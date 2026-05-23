module surfacewater_mod
   use error_mod, only: fatalerr_collected
!! Module for calculating surface water balance and drainage fluxes
!!
!! This module handles the surface water system dynamics including:
!!
!! * Lateral drainage fluxes from soil to surface water
!! * Surface water level management (automatic weirs, fixed weirs)
!! * Water supply and discharge calculations
!! * Distribution of drainage over soil compartments
!!
!! The module supports both input-based and simulated surface water levels
!! and includes functionality for macropore drainage and extended drainage systems.
!!
!! @note The surface water calculations are sensitive to timestep size and may
!! trigger automatic timestep reduction if oscillations occur.
!! @endnote
      use distribute_drainage, only: DIVDRA
      use drainage_mod, only: bocodre
      implicit none
      public :: SurfaceWater, surfacewater_year_reset
      contains

subroutine SurfaceWater(task, state, request_smaller_dt)
      !! Main driver subroutine for surface water calculations
      !!
      !! This subroutine is called with different task numbers to perform
      !! different stages of the surface water calculation:
      !!
      !! * Task 1: Initialization - read drainage input data and initialize variables
      !! * Task 2: Calculate lateral drainage fluxes to surface water
      !! * Task 3: Calculate surface water balance
      !!
      !! The subroutine manages the partitioning of drainage fluxes over soil
      !! compartments and handles both primary and secondary drainage systems.
      ! [SS-GR-FINAL B10] blanket use Variables narrowed; madr/mawlp → swap_array_dimensions
      use swap_array_dimensions, only: madr, mawlp
      ! [GR-DRA 2026-05-23] swdislay/swtopdislay/fTopDisLay retired — aliased from state%drainage below.
      ! [GR-DRA 2026-05-23] wlptab retired — read via state%surfacewater%wlptab.
      ! state%cfg%surface_water%swsrf/state%cfg%surface_water%swsec retired (→state%cfg%surface_water%X); pilot
      use array_utils, only: afgen
      use swap_state_mod, only: swap_state_t
      implicit none
      integer,            intent(in)    :: task
      type(swap_state_t), intent(inout) :: state
      logical,            intent(out)   :: request_smaller_dt

!     local
      integer level, node
      ! ADR 0031 Phase 2 Task 5: zTopDisLay declared local (was global ztopdislay).
      real(8) zCum,zTopDisLay(madr),difzTopDisLay(madr),ratio,ratiodz,sumqdr(madr),dh
      integer nodeTopDisLay(madr)
      character(len=300) messag

      request_smaller_dt = .false.

      associate( &
         tc_dt    => state%timecontrol%dt,    &  ! TC-8: SurfaceWater TC reader cutover
         tc_t1900 => state%timecontrol%t1900, &  ! TC-8
         ms_numnod    => state%mesh%numnod,             &  ! GR-BH Task 28
         ms_dz        => state%mesh%dz,                 &  ! GR-BH Task 28
         ms_layer     => state%mesh%layer,              &  ! GR-BH Task 28
         dr_nrlevs    => state%drainage%nrlevs,         &  ! GR-BH Task 28
         dr_swdivd    => state%drainage%swdivd,         &  ! GR-BH Task 28
         dr_swnrsrf   => state%drainage%swnrsrf,        &  ! GR-BH Task 28
         dr_swdivdinf  => state%drainage%swdivdinf,     &  ! GR-BH Task 28
         dr_swtopnrsrf => state%drainage%swtopnrsrf,    &  ! GR-BH Task 28
         dr_FacDpthInf => state%drainage%FacDpthInf,   &  ! GR-BH Task 28
         dr_zbotdr    => state%drainage%zbotdr,         &  ! GR-BH Task 28
         dr_L         => state%drainage%L,              &  ! GR-BH Task 28
         dr_owltab    => state%drainage%owltab,         &  ! GR-BH Task 28
         sw_ksatfit   => state%soilwater%ksatfit,       &  ! GR-BH Task 28
         sw_ksatexm   => state%soilwater%ksatexm,       &  ! GR-BH Task 28
         sw_cofani    => state%soilwater%cofani,        &  ! GR-BH Task 28
         swdislay     => state%drainage%swdislay,       &  ! [GR-DRA 2026-05-23]
         swtopdislay  => state%drainage%swtopdislay,    &
         ftopdislay   => state%drainage%ftopdislay      )

! ----------------------------------------------------------------------
      select case (task)
      case (1)

! === initialization — HOISTED ==========================================
! state%surfacewater%init(...) is now called directly from swap_init's S4 pass
! (spec docs/superpowers/specs/2026-05-13-state-init-pilot-surfacewater-design.md).
! Stub kept to preserve the task=1 dispatcher signature; case removal is a
! separate cleanup follow-up.

      return

      case (2)

! === lateral drainage fluxes to surface water ========================

! --- reset intermediate surface water and drainage fluxes
      ! SS-CRR Phase A Task A4: cohort-owned reset; see surfacewater_state_mod.
      if (state%timecontrol%flZeroIntr) call state%surfacewater%reset_intermediate()

! --- reset cumulative surface water and drainage fluxes
      ! SS-CRR Phase A correction Task 2: cumulative cohort partitioned by activity gate.
      ! SurfaceWater(2) is the canonical owner for both cohorts under swdra=2 (the only
      ! swdra value at which this code runs — see flSurfaceWater wiring in timecontrol).
      if (state%timecontrol%flZeroCumu) then
         call state%surfacewater%reset_cumulative_drainage()
         call state%surfacewater%reset_cumulative_reservoir()
      end if

! --- calculate lateral drainage
      ! SS-DRST Phase 2 Task 4: bocodre writes state%drainage%qdrain directly; no bridge sync.
      call bocodre (dh, state)

! --- partition drainage flux over compartments

      ! ADR 0031 Phase 2 Task 5: qdra global deleted; state%drainage%qdra is
      ! the sole working array throughout.  Sync loops removed.
      ! GR-BH Task 28: nrlevs/numnod/dz/swdivd/ksatfit/ksatexm/layer/cofani/l/zbotdr/
      !   Swdivdinf/Swnrsrf/SwTopnrsrf/FacDpthInf/owltab → state aliases.
      do level=1,dr_nrlevs
         do node = 1,ms_numnod
            state%drainage%qdra(level,node) = 0.0d0
         end do
      end do

      if (dr_swdivd.eq.1) then
!cD        do level=1,nrlevs
!cD           qdrain_old(level) = qdrain(level)
!cD        end do
         ! SS-SWC Phase 2 S-2.8: gwl read from state%soilwater
         call divdra (ms_numnod,dr_nrlevs,ms_dz,sw_ksatfit,sw_ksatexm,state%soilwater%fluseksatexm,    &  ! [SS-SWC S-2.12B]
            ms_layer,sw_cofani,state%soilwater%gwl,dr_L,state%drainage%qdrain,state%drainage%qdra,dr_swdivdinf,dr_swnrsrf, &
     &      dr_swtopnrsrf,dr_zbotdr,tc_dt,dr_FacDpthInf,dr_owltab,state%drainage%nowltab,tc_t1900)  ! [GR-DRA 2026-05-23]

!        redistribute qdrain with new top boundary for discharge layers
         if(swdislay.eq.2) then
            do level=1,dr_nrlevs
               if(swtopdislay(level).eq.1)  then
                  zTopDisLay(level) = fTopDisLay(level) * state%soilwater%gwl  +        &
     &                       (1.0d0-fTopDisLay(level)) * (state%soilwater%gwl-dh)
               end if
            end do
         end if
         if(swdislay.eq.1 .or. swdislay.eq.2) then
            do level=1,dr_nrlevs
               if(swtopdislay(level).eq.1)  then
!                 find node nr of new top of discharge layer
                  nodeTopDisLay(level) = 1
                  zCum               = - ms_dz(1)
                  do while (zTopDisLay(level) .lt. zCum)
                     nodeTopDisLay(level) = nodeTopDisLay(level) + 1
                     zCum              = zCum - ms_dz(nodeTopDisLay(level))
                  enddo
!                 saturated part (difzTopDisLay(lev)) of compartment containing waterlevel
                  difzTopDisLay(level) = zTopDisLay(level) - zCum
                  ratiodz =                                             &
     &                     difzTopDisLay(level)/ms_dz(nodeTopDisLay(level))
                  sumqdr(level) =                                       &
     &                        ratiodz * state%drainage%qdra(level,nodeTopDisLay(level))
                  do node = nodeTopDisLay(level)+1,ms_numnod
                     sumqdr(level) =  sumqdr(level) + state%drainage%qdra(level,node)
                  end do
                  if( dabs(sumqdr(level)) .lt. 1.0d-8)then
                     ratio = 1.0d0
                  else
                     ratio = state%drainage%qdrain(level)/sumqdr(level)
                  end if
!                 redistribute drainwater fluxes
                  do node = 1,nodeTopDisLay(level)-1
                     state%drainage%qdra(level,node) = 0.0d0
                  end do
                  state%drainage%qdra(level,nodeTopDisLay(level)) =                    &
     &                state%drainage%qdra(level,nodeTopDisLay(level)) * ratio * ratiodz
                  do node = nodeTopDisLay(level)+1,ms_numnod
                     state%drainage%qdra(level,node) = state%drainage%qdra(level,node)* ratio
                  end do
               endif
            end do
         endif

! --- error handling
!cD        do level=1,nrlevs
!cD           qdrain_new(level) = 0.0d0
!cD           do node = 1,numnod
!cD              qdrain_new(level) = qdrain_new(level) + state%drainage%qdra(level,node)
!cD           end do
!cD           if (abs(qdrain_new(level)-qdrain_old(level)).gt.0.001) then
!cD              write(messag,55)
!cD    &        ' SwDislay qdrain-diff,datetime',date,t1900,'level',level,
!cD    &        ' qdrain_new = ',qdrain_new(level),' qdrain_old = ',
!cD    &         qdrain_old(level),'(cm/d)'
!cD  55          format(a,a11,f10.4,a,i3,a,f10.5,a,f10.5,a)
!c!D              call fatalerr ('Surfacwater',messag)
!cD              call warn ('Surfacwater',messag,logf,swscre)
!cD           endif
!cD        end do

      else
! --- drainage flux through lowest compartment
        do level = 1,dr_nrlevs
           state%drainage%qdra(level,ms_numnod) = state%drainage%qdrain(level)
        end do
      endif

      ! SS-SWST Phase 2 Task 11: qdrtot global write dropped; only state written.
      state%surfacewater%qdrtot = 0.0d0
      do level=1,dr_nrlevs
          state%surfacewater%qdrtot = state%surfacewater%qdrtot + state%drainage%qdrain(level)
      end do

      return

      case (3)

! === surface water balance ========================

      if (state%cfg%surface_water%swsrf .eq. 3) then
        state%surfacewater%wlp = afgen (state%surfacewater%wlptab,2*mawlp,tc_t1900-1.0d0+tc_dt)  ! [TC-8]
      endif
      if (state%cfg%surface_water%swsec.eq.2) then
! ---    water level of secondary system is simulated
         call wlevbal (state, request_smaller_dt)
      elseif (state%cfg%surface_water%swsec.eq.1) then
! ---    water level of secondary system is input
         call wballev (state)
      endif

      case default
         call fatalerr_collected ('SurfaceWater', 'Illegal value for TASK')
      end select

      end associate  ! tc_dt/tc_t1900 => state%timecontrol [TC-8]

      return
      end

      SUBROUTINE WLEVBAL (state, request_smaller_dt)
      !! Calculate surface water level from water balance (simulated level)
      !!
      !! This subroutine determines the surface water level based on a complete
      !! water balance calculation. It handles:
      !!
      !! * Target level determination (automatic or fixed weir)
      !! * Storage calculation at target levels
      !! * Supply and discharge calculations
      !! * Overflow conditions
      !! * Ponding limitations
      !!
      !! The routine uses an iterative procedure to balance incoming drainage
      !! with discharge through weirs or open channels.
      !!
      !! @warning May trigger timestep reduction if water level oscillations
      !! exceed threshold or if ponding becomes excessive.
      !! @endwarning
      !! @note
      !! ----------------------------------------------------------------------
      !!     UpDate             : 20080109
      !!     Date               : 19990929
      !!     Purpose            : calculate surf. water level from water balance
      !!
      !! --- Set target level wlstar:
      !! --- SWMAN = 1: HBWEIR
      !! --- SWMAN = 2: from table (4e #2) and within adjustment period and
      !! --- taking into account the maximum drop rate.
      !! --- Calculate level for maximum supply wlstara ( = wlstar-wldip)
      !!
      !! --- Calculate storage at wlstar and at wlstara
      !! ---   If wlstara is above deepest bottom level then water supply
      !! ---   capacity is set to maximum value otherwise both swsttara and
      !! ---   wsmax are set to zero
      !!
      !! --- Calculate max. new storage: old + incoming/outgoing fluxes
      !! --- 1 System falls dry: set supply to maximum,
      !!       set wls at bottom of deepest channel
      !! --- 2 System will not become full: set supply to maximum
      !!       calculate new level from storage
      !! --- 3 System becomes full under maximum supply conditions,
      !! ---   Determine how, first without supply
      !! --- 3a wlstara cannot be reached: calculate required supply,
      !! ---    to reach wlstara and set wls to wlstara
      !! --- 3b system can be filled above wlstara, however no discharge:
      !! ---    wls is calculated from storage and is somewhere between
      !! ---    wlstar and wlstara, no supply, no discharge
      !! --- 3c wlstar can be reached :wsupp = 0;
      !! ---    Now check whether there are automatic weirs for keeping the
      !! ---    level at target value or that the discharge relationship
      !! ---    determines new level:
      !! --- 3c1 SWMAN=2: calculate discharge
      !! ---     check whether there is enough discharge capacity:
      !! ---     in case SWQHR = 1: calculate discap
      !! ---     in case SWQHR = 2: find discap from table
      !! ---     if sufficient capacity wls = wlstar otherwise set overflow
      !! --- 3c2 SWMAN=1 or overflow
      !! ---     Calculate highest possible discharge, if still unsufficient
      !! ---     to handle present drain fluxes: stop - system overflow
      !! ---     otherwise start iteration procedure to determine new level,
      !! ---     storage and discharge:
      !!
      !! ***********************************************
      !! --- Iteration Procedure:
      !! --- Establish lower and upper bounds: hbweir(imper) and +100cm
      !! --- Calculate storage and discharge for point halfway: swsti, wdisi
      !! --- Calculate new storage: if higher than swsti then adjust lower
      !! --- bound to point halfway otherwise adapt upper bound to point halfway
      !! --- Continue until upper and lower bounds converge (< 0.001 cm)
      !! --- Ready: update wls, swst and wdis
      !! ***********************************************
      !!
      !!     Subroutines called :
      !!     Functions called   : swstlev
      !!     File usage         :
      !!@endnote
      ! SS-SWST Phase 2 Task 11: wlstar global removed; use sw_wlstar (state alias) throughout.
      ! SS-BND Phase 2 Task B-2.4: runots removed from use clause; read via state%soilwater%runots.
      ! SS-SWC Phase 2 S-2.8: gwl,pond,THETA,THETAS,H removed from use-list; read from state%soilwater.
      ! [SS-TC TC-14] T retired — read via state%timecontrol%t
      ! GR-BH Task 28: zbotdr/NUMNOD/DZ off variables → state%drainage/mesh aliases.
      ! [GR-DRA 2026-05-23] NRPRI/nmper/impend/wldip/intwl/osswlm/wscap/dropr
      ! retired — aliased from state%surfacewater below.
      ! [GR-DRA 2026-05-23] swman/hbweir/wlsman/gwlcrit/nphase/VCRIT/NODHD/HCRIT/SWQHR/QQHTAB
      ! retired — aliased from state%surfacewater below.
      use variables, only: &   ! [SS-GR-FINAL B10] residuals — all DEFERRED
         ! [GR-DRA 2026-05-23] alphaw/betaw retired — aliased from state%surfacewater below.
         ! rsro/pondmx retired (→state%surfacewater%X)

         ! DEFERRED: QRapDra — rapid drainage flux runtime state; Phase C3
         QRapDra  ! [TC-8: dropped tcum,dt,t1900,fldtmin]
      use swap_state_mod, only: swap_state_t
      use surfacewater_utils, only: wlevst, swstlev, qhtab
      use swap_log, only: log_warn
      IMPLICIT NONE

      type(swap_state_t), intent(inout) :: state
      logical,            intent(inout) :: request_smaller_dt

! --- local
      INTEGER iphase,NODE,Intday,imper
      real(8) wlstx,swsttar,dvmax,swstmax,wsupp,wdis,wlstarb
      real(8) wover,discap,wlsl,wlsu,wlsi,swsti,wdisi,swstn
      real(8) wprod1,wprod2,oscil,wlstara
      real(8) swsttara,rday,wsmax
      character(len=200) messag
      character(len=19)  datetime
      logical :: fl_early_return
! removed the blanket save statement to avoid issues in parallel runs
!-----------------------------------------------------------------------

      fl_early_return = .false.

      associate( &
         sw_wls    => state%surfacewater%wls,    &
         sw_wlstar => state%surfacewater%wlstar, &
         sw_swst   => state%surfacewater%swst,   &
         sw_wlsbak => state%surfacewater%wlsbak, &
         sw_overfl => state%surfacewater%overfl, &
         sw_numadj => state%surfacewater%numadj, &
         sw_hwlman => state%surfacewater%hwlman, &
         sw_vtair  => state%surfacewater%vtair,  &
         sw_imper  => state%surfacewater%imper,  &
         sw_cqdrd  => state%surfacewater%cqdrd,  &
         sw_cwsupp => state%surfacewater%cwsupp, &
         sw_cwout  => state%surfacewater%cwout,  &
         tc_dt      => state%timecontrol%dt,     &  ! TC-8: WLEVBAL TC reader cutover
         tc_t1900   => state%timecontrol%t1900,  &  ! TC-8
         tc_tcum    => state%timecontrol%tcum,   &  ! TC-8
         tc_fldtmin  => state%timecontrol%fldtmin,  &  ! TC-8
         swscre      => state%timecontrol%swscre,   &  ! [SS-BMI2 Task 4]
         dr_zbotdr   => state%drainage%zbotdr,      &  ! GR-BH Task 28
         ms_numnod   => state%mesh%numnod,          &  ! GR-BH Task 28
         ms_dz       => state%mesh%dz,              &  ! GR-BH Task 28
         nrpri  => state%surfacewater%nrpri,        &  ! [GR-DRA 2026-05-23]
         nmper  => state%surfacewater%nmper,        &
         impend => state%surfacewater%impend,       &
         wldip  => state%surfacewater%wldip,        &
         intwl  => state%surfacewater%intwl,        &
         osswlm => state%surfacewater%osswlm,       &
         wscap  => state%surfacewater%wscap,        &
         dropr  => state%surfacewater%dropr,        &
         swqhr   => state%surfacewater%swqhr,       &  ! [GR-DRA 2026-05-23] weir cluster
         swman   => state%surfacewater%swman,       &
         hbweir  => state%surfacewater%hbweir,      &
         wlsman  => state%surfacewater%wlsman,      &
         gwlcrit => state%surfacewater%gwlcrit,     &
         nphase  => state%surfacewater%nphase,      &
         vcrit   => state%surfacewater%vcrit,       &
         nodhd   => state%surfacewater%nodhd,       &
         hcrit   => state%surfacewater%hcrit,       &
         alphaw  => state%surfacewater%alphaw,      &
         betaw   => state%surfacewater%betaw        )

! --- resetting of flag for overflowing of automatic weir
      ! overfl global write dropped: only sw_overfl (state alias) used henceforth.
      sw_overfl = .false.

! --- memorizing previous target level
      ! SS-SWST Phase 2 Task 11: wlstar global removed; sw_wlstar (state alias) is authoritative.
      wlstarb = sw_wlstar

! --- determine which management period the model is in:
      imper = 0
      sw_imper = 0
 100  imper = imper + 1
      sw_imper = imper

! --- error handling
      if (imper .gt. nmper) then
        messag = ' sw-level oscillation at '//datetime//                &
     &        '       advise: reduction of dtmax !'
        messag = 'error sw-management periods(IMPER), more than defined'
        call fatalerr_collected ('Wlevbal',messag)
      endif

      if (tc_t1900-1.d0+0.1d-10 .gt. impend(imper)) goto 100  ! [TC-8]

! --- determine the target sw-level:
      if (swman(imper) .eq. 1) then

! --- In the case of a fixed weir the 'target level' is set to the
! --- weir crest, for later use in calculations to determine whether
! --- there is any outflow at all (see below):
        sw_wlstar = hbweir(imper)
      else

! ---  For automatic weir, determine the target level of the surface
! ---  water, dependent on groundwater level:

! ---   only adjust it if new subperiod, with length intwl(imper) has
!       started  (or if it is is the first call)
        rday = (state%timecontrol%t+1.0D0)/intwl(imper)
        intday = int(rday)

        if (abs(rday-1.0*intday).lt.0.00001d0 .or. tc_tcum.lt.1.0d-10) then  ! [TC-8]

          iphase = nphase(imper)
          ! SS-SWC Phase 2 S-2.8: gwl read from state%soilwater
          do while(state%soilwater%gwl.gt.gwlcrit(imper,iphase).and.iphase.gt.1)
            iphase = iphase-1
          enddo

! --- compare total air volume with VCRIT, adapt iphase
          ! VTAIR global write dropped; sw_vtair (state alias) used as accumulator.
          ! SS-SWC Phase 2 S-2.8: THETAS/THETA read from state%soilwater
          sw_vtair = 0.0d0
          do NODE = 1,ms_numnod
            sw_vtair = sw_vtair + (state%soilwater%thetas(NODE)-state%soilwater%theta(NODE)) &
     &              *abs(ms_dz(NODE))
          enddo
          do while (sw_vtair.lt.VCRIT(imper,iphase).AND.IPHASE.gt.1)
            iphase = iphase - 1
          enddo

! --- compare H(nodhd(imper)) with HCRIT, adapt iphase
          ! SS-SWC Phase 2 S-2.8: H read from state%soilwater
          do while (state%soilwater%h(nodhd(imper)).gt.hcrit(imper,iphase) &
     &                            .and.iphase.gt.1)
            iphase = iphase - 1
          enddo
          ! hwlman global write dropped; sw_hwlman (state alias) set directly.
          sw_hwlman = state%soilwater%h(nodhd(imper))

          wlstx = wlsman(imper,iphase)
        else

! --- use old level
          wlstx = wlstarb
        endif

! ---   if the level must drop, then do not let it drop at more than
! ---   the specified rate:
        if (wlstx .lt. sw_wlstar .and. dropr(imper) .gt. 0.001d0) then
          sw_wlstar = sw_wlstar - dropr(imper)*tc_dt  ! [TC-8]
          if (sw_wlstar .lt. wlstx) sw_wlstar = wlstx
        else
          sw_wlstar = wlstx
        endif
      endif

! --- counter of adjustments
      if (abs(sw_wlstar-wlstarb) .gt. 0.00001d0) then
         ! numadj global write dropped; sw_numadj (state alias) accumulates directly.
         sw_numadj = sw_numadj + 1
      endif

! --- storage for the 'target level'
      swsttar = swstlev(state, sw_wlstar)

! --- level and storage for "max. level for supply"
      wlstara = sw_wlstar - wldip(imper)
      if (wlstara .gt. (dr_zbotdr(1+nrpri)+1.d-4)) then
         swsttara = swstlev(state, wlstara)
         wsmax = wscap(imper)
      else
         swsttara = 0.0d0
         wsmax = 0.0d0
      endif

! --- determine whether the system will become full (target level or
! --- level of weir crest):
      dvmax = (state%drainage%qdrd + QRapDra + wsmax) * tc_dt + state%soilwater%runots  ! [TC-8]
      swstmax = sw_swst + dvmax

      if (swstmax .lt. 1.0d-7) then
! ---   storage decreases to zero, then the surface water system
! ---   falls dry; set surface water supply to maximum:

        if (swstmax .lt. -0.1d0) then
          messag = 'error algorithm for sw falling dry'
          call fatalerr_collected ('Wlevbal',messag)
        endif

        wsupp = wsmax
        wdis = 0.0d0
        sw_swst = 0.0d0
        sw_wls = dr_zbotdr(nrpri+1)

      elseif (swstmax .ge. 0.0d0 .and. swstmax .lt. swsttara) then
! --- system will not become full - set supply to max. capacity
        wsupp = wsmax
        wdis = 0.0d0
        sw_swst = swstmax

! --- calculate new level from storage
        sw_wls = wlevst(state, sw_swst)
      else

! --- determine how system will become full: with or without needing
!     surface water supply; try first without any supply:
        dvmax = (state%drainage%qdrd + QRapDra) * tc_dt + state%soilwater%runots  ! [TC-8]
        swstmax = sw_swst + dvmax
        if (swstmax .le. swsttara) then

! --- apparently supply is needed for reaching target level, system
!     is made full up to level wlstara, because supply is controllable:
          wsupp = (swsttara - sw_swst - (state%drainage%qdrd+QRapDra)*tc_dt-state%soilwater%runots)/tc_dt  ! [TC-8]
          wdis = 0.0d0
          sw_swst = swsttara
          sw_wls = wlstara
        elseif (swstmax .le. swsttar) then

! --- apparently drainage is more than sufficient for filling system
!     to a level above the target level FOR SUPPLY, but not enough
!     for generating discharge; so calculate level from storage:
          wsupp = 0.0d0
          wdis = 0.0d0
          sw_swst = swstmax
          sw_wls = wlevst(state, sw_swst)
        else

! --- drainage water is more than sufficient for reaching target
!     level;  now see whether there are automatic weirs for keeping
!     the level at target value, or that the discharge relationship
!     determines new level:
          wsupp = 0.
          if (swman(imper) .eq. 2) then

! --- the outflow equals the drainage flux, plus the storage
!     excess (or deficit !) in the target situation compared to
!     the actual situation:
            wdis = (sw_swst-swsttar + (state%drainage%qdrd+QRapDra)*tc_dt + state%soilwater%runots )/tc_dt  ! [TC-8]

! --- now check whether the weir has enough discharge capacity
!     at this water level
            if (SWQHR.eq.1) then
              wover = sw_wls - hbweir(imper)
              discap = alphaw(imper) * (wover**betaw(imper))
            elseif (SWQHR.eq.2) then

! --- interpolate QH table
              ! SS-SWST Phase 2 Task 11: pass imper explicitly (no longer a global).
              ! GR-UTILS Task 12: pass state%surfacewater as first arg.
              discap = qhtab(state%surfacewater, sw_wlstar, imper)
            endif
            if (discap .gt. wdis) then
              sw_wls = sw_wlstar
              sw_swst = swsttar
              ! overfl global drops dropped; sw_overfl (state alias) is the signal.
              sw_overfl = .false.
            else
              sw_overfl = .true.
            endif
          endif
          if (swman(imper) .eq. 1 .or. sw_overfl) then

! --- determine level from q-h relationship, and also account
!     for change in storage (see below)
!     check first that the system does not overflow
            if (SWQHR.eq.1) then
              wover = state%surfacewater%sttab(1,1) - hbweir(imper)
              discap = alphaw(imper) * (wover**betaw(imper))
            elseif (SWQHR.eq.2) then
              discap = state%surfacewater%qqhtab(imper,1)
            endif

! ---       error handling
            swstn = sw_swst + (state%drainage%qdrd + QRapDra - discap)*tc_dt + state%soilwater%runots  ! [TC-8]
            if ( swstn .gt. state%surfacewater%sttab(1,2) ) then
              messag = 'surface water system has overflowed!'
              call fatalerr_collected ('Wlevbal',messag)
            endif

! --- iteration procedure for determining new level, storage,
!     and discharge
            wlsl = hbweir(imper)
            wlsu = state%surfacewater%sttab(1,1)

! --- find storage and discharge for intermediate point
 700        wlsi = (wlsl + wlsu) * 0.5
            swsti = swstlev(state, wlsi)
            if (SWQHR.eq.1) then
              wdisi = alphaw(imper)*(wlsi-hbweir(imper))**betaw(imper)
            else
              ! SS-SWST Phase 2 Task 11: pass imper explicitly (no longer a global).
              ! GR-UTILS Task 12: pass state%surfacewater as first arg.
              wdisi = qhtab(state%surfacewater, wlsi, imper)
            endif
            swstn = sw_swst + (state%drainage%qdrd + QRapDra - wdisi)*tc_dt + state%soilwater%runots  ! [TC-8]
            if (swstn .lt. swsti) then
              wlsu = wlsi
            else
              wlsl = wlsi
            endif
            if ((wlsu - wlsl) .gt. 0.001d0) then

! --- continue iteration procedure:
              goto 700
            else

! --- updating of sw-parameters:
              sw_wls = wlsi
              sw_swst = swstn
              wdis = wdisi
            endif
          endif
        endif
      endif

!        ponding in case of extended drainage may limit timestep

      ! SS-SWC Phase 2 S-2.8: pond read from state%soilwater
      if (sw_wls.gt.state%surfacewater%pondmx .or. state%soilwater%pond.gt.state%surfacewater%pondmx) then
        if(tc_dt .gt. 0.02*state%surfacewater%rsro) then  ! [TC-8]
          request_smaller_dt = .true.
        end if
        fl_early_return = .true.
      endif

      if (.not. fl_early_return) then
! --- updating registration of last four levels, for noting oscillation
      ! wlsbak global writes dropped; sw_wlsbak (state alias) is the ring buffer.
      sw_wlsbak(1) = sw_wlsbak(2)
      sw_wlsbak(2) = sw_wlsbak(3)
      sw_wlsbak(3) = sw_wlsbak(4)
      sw_wlsbak(4) = sw_wls
      wprod1 = (sw_wlsbak(2)-sw_wlsbak(1))*(sw_wlsbak(3)-sw_wlsbak(2))
      wprod2 = (sw_wlsbak(3)-sw_wlsbak(2))*(sw_wlsbak(4)-sw_wlsbak(3))
      if (wprod1.lt.0.0d0 .and. wprod2.lt.0.0d0) then
        oscil = abs(sw_wlsbak(3)-sw_wlsbak(2))
        if (oscil .gt. osswlm) then
           if (.not.tc_fldtmin ) then  ! [TC-8]
              request_smaller_dt = .true.
              fl_early_return = .true.
           else
              call dtdpst                                               &
     &        ('year-month-day,hour:minute:seconds',tc_t1900,datetime)  ! [TC-8]
              messag = ' sw-level oscillation at '//datetime//          &
     &        '       advise: reduction of dtmax !'
              call log_warn('Wlevbal', messag)
              call fatalerr_collected ('Wlevbal',messag)
           endif
        endif
      endif
      endif  ! .not. fl_early_return

      if (.not. fl_early_return) then
! --- cumulative terms (global accumulations dropped; state aliases are authoritative):
      sw_cqdrd  = sw_cqdrd  + state%drainage%qdrd*tc_dt  ! [TC-8]
      sw_cwsupp = sw_cwsupp + wsupp*tc_dt                ! [TC-8]
      sw_cwout  = sw_cwout  + wdis*tc_dt                 ! [TC-8]
      endif  ! .not. fl_early_return

      end associate

      return
      end

! ----------------------------------------------------------------------
      subroutine WBALLEV (state)
      !! Close surface water balance using given input surface water levels
      !!
      !! This subroutine handles the case where surface water levels are
      !! provided as input rather than calculated. It:
      !!
      !! * Fetches the water level from input time series
      !! * Calculates surface water storage at previous and current time
      !! * Determines supply or discharge as the residual of the water balance
      !! * Updates cumulative water balance terms
      !!
      !! The difference between storage at time t and t+dt, accounting for
      !! drainage fluxes, determines whether supply or discharge occurred.
      !!@note
      !!     Date               : 29/9/99
      !!     Purpose            : close surface water balance using
      !!                          given (input) surface water levels
      !!
      !! --- Compare storage at T + drainage fluxes with storage at T+DT
      !! --- The difference will be either discharge or supply
      !! --- Update totals (of discharge, supply and qdrd)
      !!
      !!     Subroutines called :
      !!     Functions called   : swstlev
      !!     File usage         :
      !!     Differences SWAP/SWAPS: None
      !!@endnote
      ! SS-BND Phase 2 Task B-2.4: runots removed from use clause; read via state%soilwater%runots.
      ! [GR-DRA 2026-05-23] wlstab retired — read via state%surfacewater%wlstab.
      use variables, only: &   ! [SS-GR-FINAL B10] residuals — all DEFERRED
         ! DEFERRED: QRapDra — rapid drainage flux runtime state; Phase C3
         QRapDra  ! [TC-8: dropped dt,t1900]
      use swap_state_mod, only: swap_state_t
      use array_utils, only: afgen
      use surfacewater_utils, only: swstlev
      use swap_array_dimensions, only: mawls
      IMPLICIT NONE

      type(swap_state_t), intent(inout) :: state

! --- local
      real(8) swstold,swstrest,wdis,wsupp

! ----------------------------------------------------------------------

      associate( &
         sw_wls    => state%surfacewater%wls,    &
         sw_wlsold => state%surfacewater%wlsold, &
         sw_swst   => state%surfacewater%swst,   &
         sw_cqdrd  => state%surfacewater%cqdrd,  &
         sw_cwsupp => state%surfacewater%cwsupp, &
         sw_cwout  => state%surfacewater%cwout,  &
         tc_dt    => state%timecontrol%dt,    &  ! TC-8: WBALLEV TC reader cutover
         tc_t1900 => state%timecontrol%t1900  )  ! TC-8

! --- wlsold gets w-level of previous time step
      ! wlsold global write dropped; sw_wlsold (state alias) is the signal.
      sw_wlsold = sw_wls

! --- fetch new level from input series
      sw_wls = AFGEN (state%surfacewater%wlstab,2*MAWLS,tc_t1900-1.d0+tc_dt)  ! [TC-8]

! --- determine surface water storage for level(t-dt) and level(t)
      swstold = swstlev(state, sw_wlsold)
      sw_swst = swstlev(state, sw_wls)

! --- determine from the surface water storages and the qdrain whether
! --- supply has taken place during period (t)-(t+dt) or water has
! --- been discharged
      swstrest = swstold + (state%drainage%qdrd + QRapDra)*tc_dt + state%soilwater%runots - sw_swst  ! [TC-8]

! --- if supply was needed, set discharge to zero
      if (swstrest.le.0.0d0) then
        wdis = 0.0d0
        wsupp = -swstrest/tc_dt  ! [TC-8]

! --- if discharge has taken place, set supply to zero
      else
        wdis = swstrest/tc_dt  ! [TC-8]
        wsupp = 0.0d0
      endif

! --- cumulation of water balance terms (global accumulations dropped; state aliases authoritative):
      sw_cqdrd  = sw_cqdrd  + state%drainage%qdrd*tc_dt  ! [TC-8]
      sw_cwsupp = sw_cwsupp + wsupp*tc_dt                ! [TC-8]
      sw_cwout  = sw_cwout  + wdis*tc_dt                 ! [TC-8]

      end associate

      return
      end


!> Year-boundary reset for surface-water cumulative state.
!! Captures the current swst as swstini for the new year.
subroutine surfacewater_year_reset(sw)
   use surfacewater_state_mod, only: surfacewater_state_t
   implicit none
   type(surfacewater_state_t), intent(inout) :: sw
   sw%swstini = sw%swst
end subroutine surfacewater_year_reset

end module surfacewater_mod
