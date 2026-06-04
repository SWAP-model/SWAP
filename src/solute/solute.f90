!> Solute transport (mobile concentration + decay + uptake + drainage).
module solute_mod
   use error_mod, only: fatalerr_collected
   implicit none
   private
   public :: solute_seed, solute_step, solute_cml_from_cmsy
   ! prepare_solute_dispersion is exposed for the state-grain unit test
   ! (test_solute_dispersion.pf); the other step phase helpers remain private.
   public :: prepare_solute_dispersion

contains

   subroutine solute_seed(state)
      use swap_state_mod, only: swap_state_t
      implicit none
      type(swap_state_t), intent(inout) :: state

      call determine_initial_solute_profile(state)
      call derive_solute_concentrations(state)
   end subroutine solute_seed

   !> Narrative orchestrator for one solute transport step. Each phase of the
   !> sub-stepped time loop is delegated to a named private helper; per-substep
   !> working arrays are passed explicitly.
   subroutine solute_step(state)
      use swap_array_dimensions, only: mabbc, macp
      use array_utils,           only: afgen
      use swap_state_mod,        only: swap_state_t

      implicit none

      type(swap_state_t), intent(inout) :: state

      real(8) :: cfluxt
      real(8) :: isqdra
      real(8) :: tcumsol
      real(8), dimension(macp) :: thetav, dispr1, vpore2

      associate (sol  => state%solute,        &
                 mesh => state%mesh,          &
                 soil => state%soilwater,     &
                 time => state%timecontrol)

         ! === Solute rate variables =============================================

         ! FIXME(swbr): SWBR=1 (mixed-reservoir aquifer breakthrough) is BROKEN and
         ! unsupported. The breakthrough blocks below read sol%bdenskfsatporos(i) with
         ! i == numnod+1 (the do-loop counter's terminal value, used outside the loop)
         ! — an out-of-range index: the array is only filled 1..numnod, so the term
         ! divides by zero and every concentration becomes NaN. This bug is present in
         ! SWAP 4.2.0 itself (confirmed: swap420gf emits NaN for all CONC columns under
         ! SWBR=1), so there is no valid parity reference to port against. The
         ! breakthrough math is kept below for a future fix, but is quarantined here.
         ! Fixing it (likely a single aquifer coefficient bdens*kfsat+poros, not a
         ! soil-node array element) and validating against the SWAP theory manual is
         ! required before SWBR=1 can be re-enabled.
         if (sol%swbr .eq. 1) then
            call state%diag%fatal('solute_step',                                     &
               'SWBR=1 (mixed-reservoir aquifer breakthrough) is not supported: '//  &
               'broken in SWAP 4.2.0 (out-of-range bdenskfsatporos read -> division '// &
               'by zero -> NaN concentrations). See FIXME in solute_step.')
            return
         end if

         ! Cohort resets (ADR 0033).
         if (time%flZeroIntr) call sol%reset_intermediate()
         if (time%flZeroCumu) then
            call sol%reset_cumulative()
            ! Mass-balance baseline rebase: samini is zeroed by reset_cumulative
            ! (cumulative cohort); physics requires anchoring it to the current
            ! profile mass (sampro) for the next balance period.
            sol%samini = sol%sampro
         end if

         sol%isqbot = 0.0d0
         sol%isqtop = 0.0d0
         isqdra     = 0.0d0

         ! Boundary concentration (afgen table).
         if (sol%swbotbc .eq. 2) then
            sol%cseep = afgen(sol%cseeptab, mabbc*2, time%t1900 + time%dt)
         end if

         call prepare_solute_dispersion(state, thetav, dispr1, vpore2)

         tcumsol = 0.0d0
         do while ((time%dt - tcumsol) .gt. 1.0d-8)

            sol%dtsolu = min(sol%dtsolu, (time%dt - tcumsol))
            sol%dtsolu = max(sol%dtsolu, time%dtmin)
            tcumsol    = tcumsol + sol%dtsolu

            call solute_surface_flux(state, cfluxt)
            call update_solute_compartments(state, thetav, dispr1, vpore2, cfluxt, isqdra)
            call solute_aquifer_breakthrough(state, isqdra)
            call solute_bottom_flux(state)

         end do

         ! Current solute flux at bottom of soil column.
         if (soil%q(mesh%numnod + 1) .gt. 0.0d0) then
            sol%isqbot = soil%q(mesh%numnod + 1) * sol%cseep
         else
            sol%isqbot = soil%q(mesh%numnod + 1) * sol%cml(mesh%numnod)
         end if

         call solute_balance(state)

      end associate

      return
   end subroutine solute_step

   !> Maximum solute time step + per-node dispersion working arrays.
   subroutine prepare_solute_dispersion(state, thetav, dispr1, vpore2)
      use swap_array_dimensions, only: macp
      use swap_state_mod,        only: swap_state_t
      implicit none
      type(swap_state_t), intent(inout) :: state
      ! thetav is fully filled 1..numnod; dispr1/vpore2 only 1..numnod-1 (the
      ! i<numnod guard below), so they take intent(inout) rather than the
      ! full-definition intent(out) contract.
      real(8), intent(out)   :: thetav(:)
      real(8), intent(inout) :: dispr1(:), vpore2(:)

      integer :: i
      real(8) :: dispr, vpore, dummy
      real(8), dimension(macp) :: diffus

      associate (sol  => state%solute,        &
                 soil => state%soilwater,     &
                 mesh => state%mesh,          &
                 time => state%timecontrol)

         ! Maximum solute time step.
         sol%dtsolu = time%dt
         do i = 1, mesh%numnod
            thetav(i) = mesh%inpola(i + 1)*soil%theta(i) + mesh%inpolb(i)*soil%theta(i + 1)
            diffus(i) = sol%ddiffwcs(i) * thetav(i)**2.33d0
            if (i .lt. mesh%numnod) then
               vpore     = abs(soil%q(i + 1))/thetav(i)
               dispr1(i) = diffus(i) + sol%ldis(mesh%layer(i)) * vpore
               vpore2(i) = vpore*vpore
            end if

            dispr = diffus(i) + sol%ldis(mesh%layer(i))*abs(soil%q(i))/soil%theta(i)
            if (dispr .lt. 1.0d-8) dispr = 1.0d-8
            dummy      = mesh%dz(i)*mesh%dz(i)/(2.0d0*dispr)
            sol%dtsolu = min(sol%dtsolu, dummy)
         end do

      end associate
   end subroutine prepare_solute_dispersion

   !> Solute flux at the soil surface (ponded reservoir mixing).
   subroutine solute_surface_flux(state, cfluxt)
      use swap_state_mod, only: swap_state_t
      implicit none
      type(swap_state_t), intent(inout) :: state
      real(8), intent(out) :: cfluxt

      real(8) :: ArMpSs   ! macropore-retired (always 0; ADR 0040)

      associate (sol  => state%solute,        &
                 soil => state%soilwater,     &
                 atmo => state%atmosphere)

         ! Macropore-area at the surface (ADR 0040: always 0).
         ArMpSs = 0.d0

         ! Solute flux at the soil surface.
         sol%csurf = (atmo%nird*sol%cirr + atmo%nraidt*sol%cpre)*sol%dtsolu + sol%csurf
         if (soil%qtop .lt. -1.d-6) then
            sol%cpond  = sol%csurf / (soil%pond - soil%qtop*sol%dtsolu)
            cfluxt     = soil%qtop*(1.0d0 - ArMpSs)*sol%cpond*sol%dtsolu
            sol%csurf  = sol%csurf + cfluxt
            sol%isqtop = soil%qtop*(1.0d0 - ArMpSs)*sol%cpond
         else
            sol%cpond = 0.0d0
            cfluxt    = 0.0d0
         end if

      end associate
   end subroutine solute_surface_flux

   !> Per-compartment mass balance: convection/dispersion, decomposition, root
   !> uptake, lateral drainage, conservation + isotherm recovery.
   subroutine update_solute_compartments(state, thetav, dispr1, vpore2, cfluxt, isqdra)
      use swap_state_mod, only: swap_state_t
      implicit none
      type(swap_state_t), intent(inout) :: state
      real(8), intent(in)    :: thetav(:), dispr1(:), vpore2(:)
      real(8), intent(inout) :: cfluxt
      real(8), intent(inout) :: isqdra

      integer :: i, level
      real(8) :: cmlav, cfluxb, cdrtot, ctrans, crot, dispr
      real(8) :: ftemp, ftheta, decact

      real(8), parameter :: vsmall = 1.0d-15

      associate (sol  => state%solute,        &
                 soil => state%soilwater,     &
                 mesh => state%mesh,          &
                 time => state%timecontrol,   &
                 drai => state%drainage,      &
                 heat => state%heat)

         ! Per-compartment mass balance.
         do i = 1, mesh%numnod

            ! Convective + dispersive fluxes.
            if (i .lt. mesh%numnod) then
               cmlav  = mesh%inpola(i + 1) * sol%cml(i) + mesh%inpolb(i) * sol%cml(i + 1)
               dispr  = dispr1(i) + 0.5d0 * sol%dtsolu*vpore2(i)
               cfluxb = (soil%q(i + 1)*cmlav +                                    &
                         thetav(i) * dispr * (sol%cml(i + 1) - sol%cml(i))/mesh%disnod(i + 1)) &
                        * sol%dtsolu
            else
               if (soil%q(i + 1) .gt. 0.0d0) then
                  cfluxb = soil%q(i + 1)*sol%cseep*sol%dtsolu
               else
                  cfluxb = soil%q(i + 1)*sol%cml(i)*sol%dtsolu
               end if
            end if

            ! Solute decomposition. NB (faithful to SWAP 4.2.0): with no soil-
            ! temperature model (flTemperature=.false.) the temperature factor is
            ! 0, which disables decomposition entirely (not a reference-temp rate).
            if (time%flTemperature) then
               if (heat%tsoil(i) .lt. 35.0d0) then
                  ftemp = exp(sol%gampar*(heat%tsoil(i) - 20.0d0))
               else
                  ftemp = exp(sol%gampar*15.0d0)
               end if
            else
               ftemp = 0.0d0
            end if
            ftheta = min(1.0d0, (soil%theta(i)/sol%rtheta)**sol%bexp)
            decact = sol%decpotfdepth(i) * ftemp * ftheta
            ctrans = decact*soil%theta(i)*sol%cml(i) +                            &
                     decact*sol%bdenskfcref(i)*((sol%cml(i)/sol%cref)**sol%frexp)
            sol%dectot   = sol%dectot   + ctrans*sol%dtsolu*mesh%dz(i)
            sol%imdectot = sol%imdectot + ctrans*sol%dtsolu*mesh%dz(i)

            ! Solute uptake by plant roots.
            crot         = sol%tscf*soil%qrot(i)*sol%cml(i)/mesh%dz(i)
            sol%rottot   = sol%rottot   + sol%tscf*soil%qrot(i)*sol%cml(i)*sol%dtsolu
            sol%imrottot = sol%imrottot + sol%tscf*soil%qrot(i)*sol%cml(i)*sol%dtsolu

            ! Lateral drainage.
            cdrtot = 0.0d0
            if (allocated(drai%qdra)) then
               do level = 1, drai%nrlevs
                  if (drai%qdra(level, i) .gt. 0.0d0) then
                     cdrtot = cdrtot + drai%qdra(level, i)*sol%cml(i)/mesh%dz(i)
                  else
                     cdrtot = cdrtot + drai%qdra(level, i)*sol%cdrain/mesh%dz(i)
                  end if
               end do
            end if

            ! Cumulative solute to lateral drainage.
            isqdra      = isqdra      + cdrtot*mesh%dz(i)*sol%dtsolu
            sol%sqdra   = sol%sqdra   + cdrtot*mesh%dz(i)*sol%dtsolu
            sol%imsqdra = sol%imsqdra + cdrtot*mesh%dz(i)*sol%dtsolu

            ! Conservation equation for the substance.
            sol%cmsy(i) = sol%cmsy(i) + (cfluxb - cfluxt) / mesh%dz(i) +          &
                          (-ctrans - crot - cdrtot) * sol%dtsolu

            ! Iterate to recover cml from cmsy with the Freundlich isotherm.
            if (sol%cmsy(i) .lt. vsmall) then
               sol%cmsy(i) = 0.0d0
               sol%cml(i)  = 0.0d0
            else
               sol%cml(i) = solute_cml_from_cmsy(sol%cmsy(i), soil%theta(i), &
                               sol%bdenskf(i), sol%frexp, sol%cref, sol%cml(i))
            end if

            cfluxt = cfluxb
         end do

      end associate
   end subroutine update_solute_compartments

   !> Aquifer breakthrough + flux to surface water (SWBR=1). QUARANTINED.
   subroutine solute_aquifer_breakthrough(state, isqdra)
      use swap_state_mod, only: swap_state_t
      implicit none
      type(swap_state_t), intent(inout) :: state
      real(8), intent(in) :: isqdra   ! cumulative solute-to-drainage from start-of-step through current substep (NOT a per-substep rate)

      ! i replicates the original post-loop counter value (numnod+1). The
      ! sol%bdenskfsatporos(i) reads below are the documented OOB/div-by-zero bug;
      ! this helper is unreachable at runtime — solute_step fatal-errors on
      ! SWBR=1 before the time loop. Kept verbatim for a future, validated fix.
      integer :: i

      associate (sol  => state%solute,        &
                 mesh => state%mesh,          &
                 surf => state%surfacewater)

         i = mesh%numnod + 1

         ! Aquifer breakthrough.
         ! FIXME(swbr): quarantined — unreachable; solute_step fatal-errors on
         ! SWBR=1 above. sol%bdenskfsatporos(i) here reads i==numnod+1 (OOB).
         if (sol%swbr .eq. 1) then
            if (surf%qdrtot .gt. 0.0d0) then
               sol%cdrain = sol%cdrain + sol%dtsolu/sol%bdenskfsatporos(i) *          &
                            ((isqdra - surf%qdrtot*sol%cdrain)/sol%daquif -       &
                             sol%decsat*sol%cdrain*sol%bdenskfsatporos(i))
            else
               sol%cdrain = sol%cdrain + sol%dtsolu/sol%bdenskfsatporos(i) *          &
                            (isqdra/sol%daquif - sol%decsat*sol%cdrain*sol%bdenskfsatporos(i))
            end if
            sol%cseep = sol%cdrain
         end if

         ! Flux to surface water from the aquifer.
         if (sol%swbr .eq. 1) then
            sol%sqsur = sol%sqsur + surf%qdrtot*sol%cdrain*sol%dtsolu
         end if

      end associate
   end subroutine solute_aquifer_breakthrough

   !> Flux through bottom of soil profile.
   subroutine solute_bottom_flux(state)
      use swap_state_mod, only: swap_state_t
      implicit none
      type(swap_state_t), intent(inout) :: state

      associate (sol  => state%solute,        &
                 soil => state%soilwater,     &
                 mesh => state%mesh)

         ! Flux through bottom of soil profile.
         if (soil%qbot .gt. 0.0d0) then
            sol%sqbot   = sol%sqbot   + soil%qbot*sol%cseep*sol%dtsolu
            sol%imsqbot = sol%imsqbot + soil%qbot*sol%cseep*sol%dtsolu
         else
            sol%sqbot   = sol%sqbot   + soil%qbot*sol%cml(mesh%numnod)*sol%dtsolu
            sol%imsqbot = sol%imsqbot + soil%qbot*sol%cml(mesh%numnod)*sol%dtsolu
         end if

      end associate
   end subroutine solute_bottom_flux

   !> Solute mass-balance components (profile mass + flux totals + residual).
   subroutine solute_balance(state)
      use swap_state_mod, only: swap_state_t
      implicit none
      type(swap_state_t), intent(inout) :: state

      integer :: i

      associate (sol  => state%solute,        &
                 mesh => state%mesh,          &
                 time => state%timecontrol,   &
                 atmo => state%atmosphere)

         ! Solute balance components.
         sol%sampro = 0.0d0
         do i = 1, mesh%numnod
            sol%sampro = sol%sampro + sol%cmsy(i) * mesh%dz(i)
         end do
         sol%sampro = sol%sampro + sol%csurf

         sol%sqprec    = sol%sqprec    + atmo%nraidt * sol%cpre * time%dt
         sol%imsqprec  = sol%imsqprec  + atmo%nraidt * sol%cpre * time%dt
         sol%sqirrig   = sol%sqirrig   + atmo%nird   * sol%cirr * time%dt
         sol%imsqirrig = sol%imsqirrig + atmo%nird   * sol%cirr * time%dt

         sol%solbal = sol%sampro - sol%sqprec - sol%sqirrig - sol%sqbot +            &
                      sol%sqdra + sol%dectot + sol%rottot - sol%samini

      end associate
   end subroutine solute_balance

   !> Recover mobile concentration cml from total cmsy via the Freundlich
   !> isotherm. Linear shortcut when frexp ~ 1; otherwise fixed-point iterate
   !> seeded from cml_guess. Caller handles the cmsy < vsmall zeroing.
   pure function solute_cml_from_cmsy(cmsy, theta, bdenskf, frexp, cref, cml_guess) result(cml)
      use, intrinsic :: iso_fortran_env, only: real64
      real(real64), intent(in) :: cmsy       ! total (dissolved+adsorbed) conc (M/L3 soil)
      real(real64), intent(in) :: theta      ! volumetric water content (-)
      real(real64), intent(in) :: bdenskf    ! bdens*kf (-)
      real(real64), intent(in) :: frexp      ! Freundlich exponent (-)
      real(real64), intent(in) :: cref       ! reference concentration (M/L3)
      real(real64), intent(in) :: cml_guess  ! previous cml, iteration seed (M/L3)
      real(real64)             :: cml

      real(real64), parameter :: rer    = 1.0d-3
      real(real64), parameter :: vsmall = 1.0d-15
      real(real64) :: old, dummy
      logical      :: differ

      if (abs(frexp - 1.0d0) .lt. 0.001d0) then
         cml = cmsy / (theta + bdenskf)
      else
         cml = cml_guess
         if (cml .lt. vsmall) cml = vsmall
         differ = .true.
         do while (differ)
            old   = cml
            dummy = bdenskf*(cml/cref)**(frexp - 1.0d0)
            cml   = cmsy/(theta + dummy)
            if (abs(cml - old) .lt. rer*cml) differ = .false.
         end do
      end if
   end function solute_cml_from_cmsy

   !> Initial solute profile from input concentrations (AFGEN table).
   subroutine determine_initial_solute_profile(state)
      use swap_array_dimensions, only: mabbc, macp
      use array_utils,           only: afgen
      use swap_state_mod,        only: swap_state_t
      implicit none
      type(swap_state_t), intent(inout) :: state
      integer :: i
      real(8) :: tab(mabbc*2)
      associate (sol => state%solute, soil => state%soilwater, mesh => state%mesh)
         if (soil%swinco .ne. 3 .and. allocated(sol%cml_init) .and. allocated(sol%zc_init)) then
            do i = 1, sol%nconc
               tab(i*2)     = sol%cml_init(i)
               tab(i*2 - 1) = abs(sol%zc_init(i))
            end do
            do i = 1, mesh%numnod
               sol%cml(i) = afgen(tab, macp*2, abs(mesh%z(i)))
            end do
         end if
      end associate
   end subroutine determine_initial_solute_profile

   !> Per-node derived coefficients + initial cmsy / mass-balance baseline.
   subroutine derive_solute_concentrations(state)
      use swap_state_mod, only: swap_state_t
      implicit none
      type(swap_state_t), intent(inout) :: state
      integer :: i
      associate (sol => state%solute, soil => state%soilwater, mesh => state%mesh)
         sol%samini = 0.0d0
         do i = 1, mesh%numnod
            sol%bdenskf(i)         = soil%bdens(mesh%layer(i))*sol%kf(mesh%layer(i))
            sol%bdenskfcref(i)     = sol%bdenskf(i)*sol%cref
            sol%bdenskfsatporos(i) = soil%bdens(mesh%layer(i))*sol%kfsat + sol%poros
            sol%cmsy(i)            = soil%theta(i)*sol%cml(i) +                          &
                                     sol%bdenskfcref(i)*(sol%cml(i)/sol%cref)**sol%frexp
            sol%samini             = sol%samini + sol%cmsy(i) * mesh%dz(i)
            sol%ddiffwcs(i)        = sol%ddif / (soil%thetsl(mesh%layer(i))**2)
            sol%decpotfdepth(i)    = sol%decpot(mesh%layer(i))*sol%fdepth(mesh%layer(i))
         end do
         sol%sampro = sol%samini
      end associate
   end subroutine derive_solute_concentrations

end module solute_mod
