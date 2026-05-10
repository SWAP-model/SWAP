! File VersionID:
!   $Id: solute.f90 341 2017-09-29 18:12:25Z kroes006 $
! ----------------------------------------------------------------------
module solute_mod
   use error_mod, only: fatalerr_collected
   implicit none
   private
   public :: solute
   public :: solute_init

contains

      subroutine solute (task, state)
! ----------------------------------------------------------------------
!     date               : december 2007; code update: June, 2019
!     purpose            : calculation of solute concentrations
! ----------------------------------------------------------------------
      use Variables
      use array_utils, only: afgen
      use swap_state_mod, only: swap_state_t
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none

!     SS-SWST Phase 2 Task 5: read qdra / qdrtot from state%surfacewater.
!     SS-SLST Phase 2 Task 8: solute writes ONLY to typed state (dual-write dropped).
      type(swap_state_t), intent(inout) :: state

!     local variables
      integer level,i,task
      real(8) cmlav,ftemp,ftheta,decact,cfluxt,cfluxb
      real(8) cdrtot,ctrans,crot,dispr,old,dummy,vpore
      real(8) isqdra,tab(mabbc*2)
      real(8) tcumsol
      logical differ
!     work arrays for intermediate calculations (recomputed each timestep)
      real(8), dimension(macp) :: thetav, diffus, dispr1, vpore2, ddiffwcs, bdenskf, bdenskfcref, bdenskfsatporos, decpotfdepth

      ! Small constants for numerical stability
      real(8), parameter :: rer = 1.0d-3
      real(8), parameter :: vsmall = 1.0d-15

! ----------------------------------------------------------------------


      select case (task)
      case (1)

! === initialize Solute rate/state variables ===========================

      associate( &
         cml    => state%solute%cml,                      &
         cmsy   => state%solute%cmsy,                     &
         samini => state%solute%cumulative%samini,         &
         sampro => state%solute%sampro                     &
      )

! --- determine initial solute profile from input concentrations
      if (swinco.ne.3) then
        do i = 1, nconc
          tab(i*2) = cml(i)
          tab(i*2-1) = abs(zc(i))
        end do
        do i = 1, numnod
          cml(i) = afgen(tab,macp*2,abs(z(i)))
        end do
      endif

! --- determine derived solute concentrations
      samini = 0.0d0
!      samaq = 0.0d0
      do i = 1, numnod
         bdenskf(i)         = bdens(layer(i))*kf(layer(i))
         bdenskfcref(i)     = bdenskf(i)*cref
         bdenskfsatporos(i) = bdens(layer(i))*kfsat + poros
         cmsy(i) = (theta(i)*cml(i) + bdenskfcref(i)*(cml(i)/cref)**frexp)
         samini = samini + cmsy(i) * dz(i)
         ddiffwcs(i) = ddif / (thetsl(layer(i))**2)
         decpotfdepth(i) = decpot(layer(i))*fdepth(layer(i))
      end do
      sampro = samini

      end associate  ! cml, cmsy, samini, sampro (case 1)

      case (2)

! === calculate Solute rate variables ========================

      associate( &
         cml      => state%solute%cml,                        &
         cmsy     => state%solute%cmsy,                       &
         csurf    => state%solute%cumulative%csurf,           &
         cpond    => state%solute%cpond,                      &
         cdrain   => state%solute%cdrain,                     &
         cseep    => state%solute%cseep,                      &
         dtsolu   => state%solute%dtsolu,                     &
         isqbot   => state%solute%isqbot,                     &
         isqtop   => state%solute%isqtop,                     &
         samini   => state%solute%cumulative%samini,          &
         sampro   => state%solute%sampro,                     &
         solbal   => state%solute%solbal,                     &
         dectot   => state%solute%cumulative%dectot,          &
         imdectot => state%solute%intermediate%imdectot,      &
         rottot   => state%solute%cumulative%rottot,          &
         imrottot => state%solute%intermediate%imrottot,      &
         sqprec   => state%solute%cumulative%sqprec,          &
         imsqprec => state%solute%intermediate%imsqprec,      &
         sqirrig  => state%solute%cumulative%sqirrig,         &
         imsqirrig=> state%solute%intermediate%imsqirrig,     &
         sqbot    => state%solute%cumulative%sqbot,           &
         imsqbot  => state%solute%intermediate%imsqbot,       &
         sqdra    => state%solute%cumulative%sqdra,           &
         imsqdra  => state%solute%intermediate%imsqdra,       &
         sqsur    => state%solute%cumulative%sqsur            &
      )

! --- reset cumulative solute fluxes
      ! SS-CRR Phase B Task B4: cohort-owned reset; see solute_state_mod.
      if (flzerointr) call state%solute%intermediate%reset()
      if (flzerocumu) then
         call state%solute%cumulative%reset()
         ! Rebase mass-balance baseline. samini is in the cumulative cohort
         ! and was zeroed by reset(); physics requires anchoring it to
         ! current profile mass (sampro) for the next balance period. See
         ! the solbal equation later in case(2).
         samini = sampro
      end if

      isqbot = 0.0d0
      isqtop = 0.0d0
      isqdra = 0.0d0

!     set value of macropore area at soil surface
      ArMpSs = 0.d0                                           !     set value of macropore area at soil surface
      if (FlMacropore .and. Z_Tp.gt.-1.d-8) ArMpSs = ArMpTp

! --- boundary concentrations
      if (swbotbc .eq. 2) then
        cseep = afgen (cseeptab,mabbc*2,t1900+dt)
      endif

! --- determine maximum timestep
      dtsolu = dt
      do i = 1,numnod
        thetav(i) = inpola(i+1)*theta(i)+inpolb(i)*theta(i+1)
        diffus(i) = ddiffwcs(i) * thetav(i)**2.33d0
        if (i < numnod) then
            vpore     = abs(q(i+1))/thetav(i)
            dispr1(i) = diffus(i) + ldis(layer(i)) * vpore
            vpore2(i) = vpore*vpore
         end if

        dispr = diffus(i)+ldis(layer(i))*abs(q(i))/theta(i)
        if (dispr.lt.1.0d-8) dispr = 1.0d-8
        !dummy = dz(i)*dz(i)*theta(i)/(2.0d0*dispr)
        dummy = dz(i)*dz(i)/(2.0d0*dispr)
        !dummy = 1.2d0*dz(i)*dz(i)/(2.0d0*dispr)
        dtsolu = min(dtsolu,dummy)
      enddo

      tcumsol = 0.0d0
      ! SS-DRST Task 3: qdra read from state%drainage; qdrtot remains in state%surfacewater
      associate(qdra   => state%drainage%qdra, &
                qdrtot => state%surfacewater%qdrtot)
      do while ((dt-tcumsol).gt.1.0d-8)

! ---    time step and cumulative time
         dtsolu  = min(dtsolu,(dt-tcumsol))
         dtsolu  = max(dtsolu,dtmin)
         tcumsol = tcumsol + dtsolu

! --- solute flux at soil surface
         csurf = (nird*cirr + nraidt*cpre)*dtsolu + csurf   ! gr cm-2
         ! SS-BND Phase 2 Task B-2.3: qtop read from state%soilwater (boundary home).
         if (state%soilwater%qtop.lt.-1.d-6) then
            cpond  = csurf / (pond-state%soilwater%qtop*dtsolu)             ! gr cm-3
            cfluxt = state%soilwater%qtop*(1.0d0-ArMpSs)*cpond*dtsolu       ! gr cm-2
            csurf  = csurf + cfluxt                                          ! gr cm-2
            isqtop = state%soilwater%qtop*(1.0d0-ArMpSs)*cpond              ! gr cm-2 d-1
         else
            cpond  = 0.0d0
            cfluxt = 0.0d0
         endif

! --- calculate mass balance for each compartment

         do i = 1,numnod

! --- convective and dispersive fluxes
            if (i .lt. numnod) then
               cmlav = inpola(i+1) * cml(i) + inpolb(i) * cml(i+1)
               !thetav = inpola(i+1)*theta(i)+inpolb(i)*theta(i+1)
               !vpore = abs(q(i+1))/thetav
               !diffus = ddif*(thetav**2.33d0)/(thetsl(layer(i))**2)
               !dispr = diffus + ldis(layer(i)) * vpore + 0.5d0 * dtsolu*vpore*vpore
               dispr = dispr1(i) + 0.5d0 * dtsolu*vpore2(i)
               cfluxb = (q(i+1)*cmlav + thetav(i) * dispr * (cml(i+1)-cml(i))/disnod(i+1))*dtsolu
            else
               if (q(i+1).gt.0.0d0) then
                  cfluxb = q(i+1)*cseep*dtsolu
               else
                  cfluxb = q(i+1)*cml(i)*dtsolu
               endif
            endif

! --- solute decomposition
            if (fltemperature) then
               ! SS-HEAT Phase 2 Task 6: read tsoil from state%heat
               if (state%heat%tsoil(i) .lt. 35.0d0) then
                  ftemp = exp(gampar*(state%heat%tsoil(i)-20.0d0))
               else
                  ftemp = exp(gampar*15.0d0)
               endif
            else
              ftemp = 0.0d0
            endif
            ftheta = min(1.0d0,(theta(i)/rtheta)**bexp)
            decact = decpotfdepth(i) * ftemp * ftheta
            ctrans = decact*theta(i)*cml(i) + decact*bdenskfcref(i)*((cml(i)/cref)**frexp)
            dectot = dectot + ctrans*dtsolu*dz(i)
            imdectot = imdectot + ctrans*dtsolu*dz(i)

! --- solute uptake by plant roots
            crot   = tscf*qrot(i)*cml(i)/dz(i)
            rottot = rottot + tscf*qrot(i)*cml(i)*dtsolu
            imrottot = imrottot + tscf*qrot(i)*cml(i)*dtsolu

! --- lateral drainage
            ! SS-SWST Phase 2 Task 11: qdra read from state (global dropped).
            cdrtot = 0.0d0
            if (allocated(state%drainage%qdra)) then
            do level = 1,nrlevs
               if (state%drainage%qdra(level,i) .gt. 0.0d0) then
                  cdrtot = cdrtot+state%drainage%qdra(level,i)*cml(i)/dz(i)
               else
                  cdrtot = cdrtot+state%drainage%qdra(level,i)*cdrain/dz(i)
               endif
            enddo
            end if

! --- cumulative amount of solutes to lateral drainage
            isqdra = isqdra + cdrtot*dz(i)*dtsolu
            sqdra  = sqdra + cdrtot*dz(i)*dtsolu
            imsqdra  = imsqdra + cdrtot*dz(i)*dtsolu

! --- conservation equation for the substance
            cmsy(i) = cmsy(i) + (cfluxb-cfluxt) / dz(i) + (-ctrans-crot-cdrtot) * dtsolu

! --- iteration procedure for calculation of cml
            differ = .true.
            if (cmsy(i).lt.vsmall) then
               cmsy(i) = 0.0d0
               cml(i)  = 0.0d0
            else
               if (abs(frexp-1.0d0).lt.0.001d0) then
                  cml(i) = cmsy(i) / (theta(i) + bdenskf(i))
               else
                  if (cml(i).lt.vsmall) cml(i) = vsmall
                  do while (differ)
                     old    = cml(i)
                     dummy  = bdenskf(i)*(cml(i)/cref)**(frexp-1.0d0)
                     cml(i) = cmsy(i)/(theta(i)+dummy)
                     if (abs(cml(i)-old).lt.rer*cml(i)) differ = .false.
                  enddo
               endif
            endif

!           make top flux next compartment equal to current bottom flux
            cfluxt = cfluxb

! ---    next compartment
         enddo

! ---    solute balance in aquifer for breakthrough curve
         ! SS-SWST Phase 2 Task 11: qdrtot read from state (global dropped).
         if (swbr .eq. 1) then
            if (state%surfacewater%qdrtot .gt. 0.0d0) then
               cdrain = cdrain + dtsolu/bdenskfsatporos(i) *         &
     &         ( (isqdra - state%surfacewater%qdrtot*cdrain)/daquif - decsat*cdrain*bdenskfsatporos(i) )
            else
               cdrain = cdrain + dtsolu/bdenskfsatporos(i) *         &
     &                           ( isqdra/daquif - decsat*cdrain*bdenskfsatporos(i) )
            endif
            cseep = cdrain
         endif

! --- flux to surface water from aquifer
         if (swbr .eq. 1) then
            sqsur = sqsur + state%surfacewater%qdrtot*cdrain*dtsolu
         endif

! --- flux through bottom of soil profile
         ! SS-BND Phase 2 Task B-2.3: qbot read from state%soilwater (boundary home).
         if (state%soilwater%qbot .gt. 0.0d0) then
            sqbot = sqbot + state%soilwater%qbot*cseep*dtsolu
            imsqbot = imsqbot + state%soilwater%qbot*cseep*dtsolu
         else
            sqbot = sqbot + state%soilwater%qbot*cml(numnod)*dtsolu
            imsqbot = imsqbot + state%soilwater%qbot*cml(numnod)*dtsolu
         endif

! --- continue with next solute time step
      end do
      end associate  ! qdra, qdrtot from state%surfacewater

! --- current solute flux at bottom of soil column
      if (q(numnod+1) .gt. 0.0d0) then
        isqbot = q(numnod+1) * cseep
      else
        isqbot = q(numnod+1) * cml(numnod)
      endif

! === calculate solute balance components ========================

! --- total amount in soil profile
      sampro = 0.0d0
      do i = 1,numnod
        sampro = sampro + cmsy(i) * dz(i)
      enddo
      sampro = sampro + csurf
!      if (swbr .eq. 1) samaq = cdrain*poros*daquif

! --- add time step fluxes to total cumulative values
      sqprec = sqprec + nraidt * cpre * dt
      imsqprec = imsqprec + nraidt * cpre * dt
      sqirrig = sqirrig + nird * cirr * dt
      imsqirrig = imsqirrig + nird * cirr * dt

! --- cumulative solute balance
      solbal = sampro - sqprec - sqirrig - sqbot + sqdra + dectot + rottot - samini

      end associate  ! cml, cmsy, csurf, cpond, ..., sqsur (case 2)

      case default
         call fatalerr_collected ('Solute', 'Illegal value for TASK')
      end select

      return
      end subroutine solute

!> Lifecycle init for solute typed state. Allocates per-node arrays
!! from numnod and seeds them from the config-time-populated legacy
!! globals. Called from swap_main once per simulation, after
!! config_to_variables has populated the legacy globals.
!!
!! Mirrors drainage_init's pattern (ADR 0031). The two-stage cml
!! seeding (discovery hazard #7) is preserved: config-time seed in
!! the legacy global persists; this routine copies it into the typed
!! state at run init.
subroutine solute_init(state)
   use, intrinsic :: iso_fortran_env, only: real64
   use swap_state_mod, only: swap_state_t
   use Variables, only: numnod, cml, cmsy
   implicit none
   type(swap_state_t), intent(inout) :: state

   if (.not. allocated(state%solute%cml)) then
      allocate(state%solute%cml(numnod))
   end if
   if (.not. allocated(state%solute%cmsy)) then
      allocate(state%solute%cmsy(numnod))
   end if

   state%solute%cml(:)  = cml(1:numnod)
   state%solute%cmsy(:) = cmsy(1:numnod)
end subroutine solute_init

end module solute_mod
