!> @file agetracer.f90
!! AgeTracer feature — extracted from solute.f90 during ADR 0032
!! to declutter the solute compute module.
!!
!! Status: DEAD. flAgeTracer is never set to .true. anywhere in the
!! codebase. The runtime body is preserved verbatim for source
!! completeness but is gated behind a stub-error guard.
!!
!! Reactivation checklist (DO NOT skip — these decisions defer until
!! AgeTracer is genuinely needed):
!!
!!   1. Define agetracer_state_t. The 12 AgeTracer-specific globals
!!      currently in variables.f90 (Ageirr, Agedrain, Agepre,
!!      Agepond, Agepondm1, icAgetopupw, icAgetopdwn, icAgeBot,
!!      icAgeDra, icAgeRot, icAgeSur, AgeGwl1m) are the target field
!!      set.
!!
!!   2. Resolve the cml dual-use (discovery hazard #1). AgeTracer
!!      currently overwrites state%solute%cml(i) with age storage at
!!      the end of each task=2 call. Either:
!!      (a) Add a separate state%agetracer%cml_age(:) field and
!!          stop overwriting state%solute%cml; OR
!!      (b) Document the dual-use sequencing rules formally and
!!          retain the overwrite (faster but more fragile).
!!
!!   3. Wire flAgeTracer from typed config or external trigger
!!      (currently never assigned anywhere).
!!
!!   4. Update outage's flAgeTracer guard (in swapoutput.f90) to
!!      ungate output once the runtime is back online.
!!
!!   5. Remove the stub-error guard at the top of AgeTracer below.
!!
!! See docs/superpowers/specs/2026-05-10-state-migration-solute-discovery.md
!! Section 8 hazard #4 + the 2026-05-10 design spec D5.

module agetracer_mod
   use error_mod, only: fatalerr_collected
   implicit none
   private
   public :: AgeTracer

contains

      subroutine AgeTracer (task, state)
! ----------------------------------------------------------------------
!     date               : Oct 2010
!     purpose            : Ageing according to Goode (1996):
!                        : "Direct simulation of groundwater age, WRR vol.32, p 289-296"
! ----------------------------------------------------------------------
! --- STUB-ERROR GUARD (ADR 0032) ---
!     AgeTracer is currently inert: flAgeTracer is never set to .true.
!     anywhere in the codebase. The runtime body is preserved verbatim
!     below for future reactivation, but it is unreachable because this
!     guard fires first. Before re-enabling, work through the
!     reactivation checklist in the file header above.
! ---
      use Variables
      use array_utils, only: afgen
      use swap_state_mod, only: swap_state_t
      implicit none

!     SS-SWST Phase 2 Task 5: read qdra from state%surfacewater.
      type(swap_state_t), intent(in) :: state

!     global
      integer task

!     local variables
      integer level,i
      real(8) Agemlav,thetav,Agefluxt,Agefluxb
      real(8) Agedrtot,Agerot,dispr,diffus,dummy
      real(8) vpore,isqdra,tab(mabbc*2)
      real(8) tcumsol
      ! [GR-BH Task 36] ArMpSs retired from variables.f90 — local (always 0.d0, ADR 0040 complete)
      real(8) ArMpSs
!
!      real(8) Ageevp
      real(8) Ageml(macp)    ! Array with age mass solute concentration (M/L3 water)
      real(8) Agemsy(macp)   ! Array with dissolved solute concentration (M/L3 soil volume)
      real(8) Agesurf
      real(8) AgeProd
      real(8) sum0,sum1,deltaz,zzbot,zztop

! Note: Ageirr, Agedrain, Agepre, Agepond, Agepondm1, icAgetopupw, icAgetopdwn, ArMpSs
! are now module-level variables in variables.f90 (synced via state%solute)
! This enables multi-instance execution

      if (flAgeTracer) then
         call fatalerr_collected('AgeTracer', &
            'AgeTracer feature is currently inert. The runtime body '// &
            'was preserved during ADR 0032 (solute migration) for '// &
            'future reactivation. Before re-enabling, work through the '// &
            'reactivation checklist in src/solute/agetracer.f90 (defines '// &
            'agetracer_state_t, resolves cml dual-use, wires flAgeTracer).')
         return
      end if

! [Original AgeTracer body from solute.f90:304-547 preserved verbatim
!  below this line. Body is unreachable in current build because the
!  stub-error guard above always exits when flAgeTracer is .true.,
!  and the caller in swap.f90 already gates on flAgeTracer being .true.
!  before calling — so this body executes only if flAgeTracer is .false.,
!  in which case it drops straight to the end select and returns.]

! SS-SWC Phase 2 S-2.8: theta/thetsl/q/thetm1/gwl/nodgwl/thetas/pond
!  read from state%soilwater (reader cutover; body is preserved but unreachable).
! SS-TC TC-12: dt read via state%timecontrol tc_* alias.
      ! [GR-BH Audit 31] numnod/dz/disnod/zbotcp/ztopcp aliased via state%mesh;
      !   nrlevs aliased via state%drainage
      ! [GR-BH Task 35] z and layer added to associate (globals deleted)
      associate( &
         tc_dt     => state%timecontrol%dt,    &  ! TC-12
         sw_theta  => state%soilwater%theta,   &
         sw_thetm1 => state%soilwater%thetm1,  &
         sw_thetas => state%soilwater%thetas,  &
         sw_thetsl => state%soilwater%thetsl,  &
         sw_q      => state%soilwater%q,       &
         sw_gwl    => state%soilwater%gwl,     &
         sw_nodgwl => state%soilwater%nodgwl,  &
         sw_pond   => state%soilwater%pond,    &
         numnod    => state%mesh%numnod,       &  ! [GR-BH Audit 31]
         dz        => state%mesh%dz,           &  ! [GR-BH Audit 31]
         z         => state%mesh%z,            &  ! [GR-BH Task 35]
         disnod    => state%mesh%disnod,       &  ! [GR-BH Audit 31]
         zbotcp    => state%mesh%zbotcp,       &  ! [GR-BH Audit 31]
         ztopcp    => state%mesh%ztopcp,       &  ! [GR-BH Audit 31]
         layer     => state%mesh%layer,        &  ! [GR-BH Task 35]
         nrlevs    => state%drainage%nrlevs    &  ! [GR-BH Audit 31]
      )

! ----------------------------------------------------------------------

      select case (task)
      case (1)

! === initialize Solute rate/state variables ===========================

! --- determine initial solute profile from input concentrations
      if (swinco.ne.3) then
        do i = 1, nconc
          tab(i*2) = cml(i)
          tab(i*2-1) = abs(zc(i))
        end do
        do i = 1, numnod
          Ageml(i) = afgen(tab,macp*2,abs(z(i)))
        end do
      else
        do i = 1, numnod
          Ageml(i) = cml(i)
        end do
      endif
! --- boundary conditions (ALL SET TO 0 ON INPUT)
      Agepre = cpre
      Ageirr = cirr
      Agedrain = cdrain
      Agepond = cpre
      Agepondm1 = Agepond

! --- determine derived solute concentrations
      samini = 0.0d0
      do i = 1,numnod
         Agemsy(i) = sw_theta(i)*Ageml(i)
         samini = samini + Agemsy(i) * dz(i)
      end do

      case (2)

! === calculate Solute rate variables ========================

! --- reset cumulative solute fluxes
      if (state%timecontrol%flZeroCumu) then
        icAgesur    = 0.0d0
        icAgetopdwn = 0.0d0
        icAgetopupw = 0.0d0
        icAgerot    = 0.0d0
        do level = 1,nrlevs
           icAgedra(level) = 0.0d0
        end do
        icAgebot = 0.0d0
      endif

      isqbot = 0.0d0
      isqtop = 0.0d0
      isqdra = 0.0d0

!     set value of macropore area at soil surface
      ArMpSs = 0.d0                                           !     set value of macropore area at soil surface
      if (FlMacropore .and. Z_Tp.gt.-1.d-8) ArMpSs = ArMpTp

! --- determine maximum timestep
      dtsolu = tc_dt  ! TC-12
      do i = 1,numnod
        ! SS-SWC Phase 2 S-2.8: theta/thetsl/q read from state%soilwater
        thetav = inpola(i+1)*sw_theta(i)+inpolb(i)*sw_theta(i+1)
        diffus = ddif * (thetav**2.33d0)/(sw_thetsl(layer(i))**2)
        dispr = diffus+ldis(layer(i))*abs(sw_q(i))/sw_theta(i)
        if (dispr.lt.1.0d-8) dispr = 1.0d-8
        dummy = dz(i)*dz(i)*sw_theta(i)/2.0/dispr
        dtsolu = min(dtsolu,dummy)
      enddo

      tcumsol = 0.0
      ! SS-DRST Task 3: qdra read from state%drainage
      ! [SS-BMI2 Task 4] dtmin aliased from state%timecontrol
      associate(qdra  => state%drainage%qdra, &
                dtmin => state%timecontrol%dtmin)
      do while ((tc_dt-tcumsol).gt.1.0d-8)  ! TC-12

! ---    time step and cumulative time
         dtsolu = min(dtsolu,(tc_dt-tcumsol))  ! TC-12
         dtsolu = max(dtsolu,dtmin)
         tcumsol   = tcumsol + dtsolu

! --- solute flux at soil surface
         ! SS-ATM Phase 2 Task A-2.4: nraidt read from state%atmosphere (atmosphere home).
         Agesurf = (nird*Ageirr + state%atmosphere%nraidt*Agepre)*dtsolu +               &
     &                                       state%soilwater%pondm1*Agepondm1              ! [SS-SWC S-2.12B] gr cm-2
         ! SS-BND Phase 2 Task B-2.3: qtop/runots read from state%soilwater (boundary home).
         ! SS-SWC Phase 2 S-2.8: pond read from state%soilwater
         if (state%soilwater%qtop.lt.-1.d-6) then
            Agepond  = Agesurf / (sw_pond-state%soilwater%qtop*dtsolu)                        ! gr cm-3
            Agefluxt = state%soilwater%qtop*(1.0d0-ArMpSs)*Agepond*dtsolu                  ! gr cm-2
            Agesurf  = Agesurf + Agefluxt                                                   ! gr cm-2
            isqtop   = state%soilwater%qtop*(1.0d0-ArMpSs)*Agepond                         ! gr cm-2 d-1
            icAgetopdwn = icAgetopdwn + state%soilwater%qtop*(1.0d0-ArMpSs)*0.5d0*      &  ! gr cm-2
     &                    (Agepond+Agepondm1)  * dtsolu
         else
            Agepond  = 0.0d0
            Agefluxt = 0.0d0
         endif
         icAgesur    = icAgesur + 0.5d0*(Agepond+Agepondm1) * state%soilwater%runots

! --- calculate mass balance for each compartment

         do i = 1,numnod

! --- convective and dispersive fluxes
            if (i .lt. numnod) then
               Agemlav = inpola(i+1) * Ageml(i) + inpolb(i) * Ageml(i+1)
               ! SS-SWC Phase 2 S-2.8: theta/thetsl/q read from state%soilwater
               thetav = inpola(i+1)*sw_theta(i)+inpolb(i)*sw_theta(i+1)
               vpore = abs(sw_q(i+1))/thetav
               diffus = ddif*(thetav**2.33d0)/(sw_thetsl(layer(i))**2)
               dispr = diffus + ldis(layer(i)) * vpore +                &
     &                         0.5d0 * dtsolu*vpore*vpore
               Agefluxb = (sw_q(i+1)*Agemlav +thetav * dispr               &
     &                  *(Ageml(i+1)-Ageml(i))/disnod(i+1))*dtsolu
            else
               ! SS-SWC Phase 2 S-2.8: q(numnod+1) read from state%soilwater
               if (sw_q(i+1).gt.0.0d0) then
                  Agefluxb = sw_q(i+1)*Agedrain*dtsolu
               else
                  Agefluxb = sw_q(i+1)*Ageml(i)*dtsolu
               endif
            endif

! --- Age tracer uptake by plant roots
            ! SS-CRP Phase 2 C-2.3: qrot read from state%soilwater
            Agerot = state%soilwater%qrot(i)*Ageml(i)/dz(i)
            rottot = rottot + state%soilwater%qrot(i)*Ageml(i)*dtsolu

! --- Age tracer uptake by soil evaporation ???
!            Ageevp = evso??*Ageml(i)/dz(i)

! --- lateral drainage
            ! SS-SWST Phase 2 Task 11: qdra read from state (global dropped).
            Agedrtot = 0.0d0
            if (allocated(state%drainage%qdra)) then
            do level = 1,nrlevs
               if (state%drainage%qdra(level,i) .gt. 0.0d0) then
                  Agedrtot = Agedrtot+state%drainage%qdra(level,i)*Ageml(i)/dz(i)
               else
                  Agedrtot = Agedrtot+state%drainage%qdra(level,i)*Agedrain/dz(i)
               endif
            enddo
            end if

! --- cumulative amount of solutes to lateral drainage
            isqdra = isqdra + Agedrtot*dz(i)*dtsolu
            sqdra = sqdra + Agedrtot*dz(i)*dtsolu

! --- zero order production
            ! SS-SWC Phase 2 S-2.8: theta/thetm1 read from state%soilwater
            AgeProd = 1.0d0 * 0.5d0*(sw_theta(i)+sw_thetm1(i))

! --- conservation equation for the substance
            Agemsy(i) = Agemsy(i) + (Agefluxb-Agefluxt) / dz(i) +       &
     &                  (-Agerot-Agedrtot+AgeProd) * dtsolu
            Ageml(i) = Agemsy(i) / sw_theta(i)

!           make top flux next compartment equal to current bottom flux
            Agefluxt = Agefluxb

! ---    next compartment
         enddo

!        age of effluent terms
         ! SS-SWC Phase 2 S-2.8: q(1) read from state%soilwater
         icAgetopupw = icAgetopupw + max(0.0d0,sw_q(1))*                   &
     &                              0.5d0*(Ageml(1)+cml(1))*dtsolu
         do i = 1,numnod
            ! SS-CRP Phase 2 C-2.3: qrot read from state%soilwater
            icAgerot = icAgerot + state%soilwater%qrot(i)*              &
     &                              0.5d0*(Ageml(i)+cml(i))*dtsolu
            do level = 1,nrlevs
               icAgedra(level) = icAgedra(level) + qdra(level,i)*       &
     &                              0.5d0*(Ageml(i)+cml(i))*dtsolu
            end do
         end do
         ! SS-SWC Phase 2 S-2.8: q(numnod+1) read from state%soilwater
         icAgebot = icAgebot - min(sw_q(numnod+1),0.0d0)*                  &
     &                         0.5d0*(Ageml(numnod)+cml(numnod))*dtsolu

      end do
      end associate  ! qdra from state%surfacewater

!     age of variable 1m-plane of groundwater
      ! SS-SWC Phase 2 S-2.8: nodgwl/gwl/thetas read from state%soilwater
      i = sw_nodgwl+1
      zzbot = zbotcp(i)
      zztop = ztopcp(i)
      sum0 = 0.0d0
      sum1 = 0.0d0
      do while(sw_gwl-100.0d0.lt.zzbot .and. i.le.numnod)
         if(sw_gwl.lt.zztop) zztop = sw_gwl
         if(sw_gwl-100.0d0.gt.zzbot) zzbot = sw_gwl-100.0d0
         deltaz = zztop - zzbot
         sum1   = sum1 + deltaz * Ageml(i) * sw_thetas(i)
         sum0   = sum0 + deltaz * sw_thetas(i)
         i = i + 1
         zztop = ztopcp(i)
         zzbot = zbotcp(i)
      end do
      Agegwl1m = -9999.9d0
      if(sum0.gt.1.0d-12) Agegwl1m = sum1/sum0


! --- initial values next timestep
      Agepondm1 = Agepond

! --- Save for output only
      do i = 1, numnod
        cml(i) = Ageml(i)
      end do

      case default
         call fatalerr_collected ('AgeTracer', 'Illegal value for TASK')
      end select

      end associate  ! sw_theta, sw_thetm1, sw_thetas, sw_thetsl, sw_q, sw_gwl, sw_nodgwl, sw_pond; numnod/dz/disnod/zbotcp/ztopcp/nrlevs [GR-BH Audit 31]

      return
      end subroutine AgeTracer

end module agetracer_mod
