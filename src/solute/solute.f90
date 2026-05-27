!> Solute transport (mobile concentration + decay + uptake + drainage).
module solute_mod
   use error_mod, only: fatalerr_collected
   implicit none
   private
   public :: solute_seed, solute_step

contains

   subroutine solute_seed(state)
      use swap_array_dimensions, only: mabbc, macp
      use array_utils,           only: afgen
      use swap_state_mod,        only: swap_state_t
      use, intrinsic :: iso_fortran_env, only: real64
      use solute_kernels_mod,    only: bdenskf_coeff, bdenskfsatporos_coeff, &
                                       ddiffwcs_coeff, decpotfdepth_coeff

      implicit none

      type(swap_state_t), intent(inout) :: state

      integer :: level, i
      real(8) :: cmlav, ftemp, ftheta, decact, cfluxt, cfluxb
      real(8) :: cdrtot, ctrans, crot, dispr, old, dummy, vpore
      real(8) :: isqdra, tab(mabbc*2)
      real(8) :: tcumsol
      real(8) :: ArMpSs   ! macropore-retired (always 0; ADR 0040)
      logical :: differ
      real(8), parameter :: rer    = 1.0d-3
      real(8), parameter :: vsmall = 1.0d-15

      associate (sol  => state%solute,        &
                 soil => state%soilwater,     &
                 mesh => state%mesh,          &
                 time => state%timecontrol,   &
                 drai => state%drainage,      &
                 surf => state%surfacewater,  &
                 atmo => state%atmosphere,    &
                 heat => state%heat)

         ! === Initialise solute rate/state variables ============================

         ! Determine initial solute profile from input concentrations.
         if (soil%swinco .ne. 3 .and. allocated(sol%cml_init) .and. allocated(sol%zc_init)) then
            do i = 1, sol%nconc
               tab(i*2)     = sol%cml_init(i)
               tab(i*2 - 1) = abs(sol%zc_init(i))
            end do
            do i = 1, mesh%numnod
               sol%cml(i) = afgen(tab, macp*2, abs(mesh%z(i)))
            end do
         end if

         ! Derived solute concentrations + time-invariant coefficients.
         sol%samini = 0.0d0
         do i = 1, mesh%numnod
            sol%bdenskf(i)         = bdenskf_coeff(soil%bdens(mesh%layer(i)), sol%kf(mesh%layer(i)))
            sol%bdenskfcref(i)     = sol%bdenskf(i)*sol%cref
            sol%bdenskfsatporos(i) = bdenskfsatporos_coeff(soil%bdens(mesh%layer(i)), sol%kfsat, sol%poros)
            sol%cmsy(i)            = soil%theta(i)*sol%cml(i) +                          &
                                     sol%bdenskfcref(i)*(sol%cml(i)/sol%cref)**sol%frexp
            sol%samini             = sol%samini + sol%cmsy(i) * mesh%dz(i)
            sol%ddiffwcs(i)        = ddiffwcs_coeff(sol%ddif, soil%thetsl(mesh%layer(i)))
            sol%decpotfdepth(i)    = decpotfdepth_coeff(sol%decpot(mesh%layer(i)), sol%fdepth(mesh%layer(i)))
         end do
         sol%sampro = sol%samini

      end associate

      return
   end subroutine solute_seed

   subroutine solute_step(state)
      use swap_array_dimensions, only: mabbc, macp
      use array_utils,           only: afgen
      use swap_state_mod,        only: swap_state_t
      use, intrinsic :: iso_fortran_env, only: real64
      use solute_kernels_mod,    only: solute_cml_from_cmsy, solute_ftemp, solute_ftheta, solute_decomp_ctrans

      implicit none

      type(swap_state_t), intent(inout) :: state

      integer :: level, i
      real(8) :: cmlav, ftemp, ftheta, decact, cfluxt, cfluxb
      real(8) :: cdrtot, ctrans, crot, dispr, dummy, vpore
      real(8) :: isqdra, tab(mabbc*2)
      real(8) :: tcumsol
      real(8) :: ArMpSs   ! macropore-retired (always 0; ADR 0040)
      real(8), dimension(macp) :: thetav, diffus, dispr1, vpore2

      real(8), parameter :: vsmall = 1.0d-15

      associate (sol  => state%solute,        &
                 soil => state%soilwater,     &
                 mesh => state%mesh,          &
                 time => state%timecontrol,   &
                 drai => state%drainage,      &
                 surf => state%surfacewater,  &
                 atmo => state%atmosphere,    &
                 heat => state%heat)

         ! === Solute rate variables =============================================

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

         ! Macropore-area at the surface (ADR 0040: always 0).
         ArMpSs = 0.d0

         ! Boundary concentration (afgen table).
         if (sol%swbotbc .eq. 2) then
            sol%cseep = afgen(sol%cseeptab, mabbc*2, time%t1900 + time%dt)
         end if

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

         tcumsol = 0.0d0
         do while ((time%dt - tcumsol) .gt. 1.0d-8)

            sol%dtsolu = min(sol%dtsolu, (time%dt - tcumsol))
            sol%dtsolu = max(sol%dtsolu, time%dtmin)
            tcumsol    = tcumsol + sol%dtsolu

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

               ! Solute decomposition.
               ftemp  = solute_ftemp(heat%tsoil(i), sol%gampar, time%flTemperature)
               ftheta = solute_ftheta(soil%theta(i), sol%rtheta, sol%bexp)
               decact = sol%decpotfdepth(i) * ftemp * ftheta
               ctrans = solute_decomp_ctrans(decact, soil%theta(i), sol%cml(i), &
                                             sol%bdenskfcref(i), sol%cref, sol%frexp)
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

            ! Aquifer breakthrough.
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

            ! Flux through bottom of soil profile.
            if (soil%qbot .gt. 0.0d0) then
               sol%sqbot   = sol%sqbot   + soil%qbot*sol%cseep*sol%dtsolu
               sol%imsqbot = sol%imsqbot + soil%qbot*sol%cseep*sol%dtsolu
            else
               sol%sqbot   = sol%sqbot   + soil%qbot*sol%cml(mesh%numnod)*sol%dtsolu
               sol%imsqbot = sol%imsqbot + soil%qbot*sol%cml(mesh%numnod)*sol%dtsolu
            end if

         end do

         ! Current solute flux at bottom of soil column.
         if (soil%q(mesh%numnod + 1) .gt. 0.0d0) then
            sol%isqbot = soil%q(mesh%numnod + 1) * sol%cseep
         else
            sol%isqbot = soil%q(mesh%numnod + 1) * sol%cml(mesh%numnod)
         end if

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

      return
   end subroutine solute_step

end module solute_mod
