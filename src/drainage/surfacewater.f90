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
      use swap_array_dimensions, only: madr, mawlp
      use array_utils,           only: afgen
      use swap_state_mod,        only: swap_state_t
      implicit none

      integer,            intent(in)    :: task
      type(swap_state_t), intent(inout) :: state
      logical,            intent(out)   :: request_smaller_dt

      integer :: level, node
      real(8) :: zCum, zTopDisLay(madr), difzTopDisLay(madr), ratio, ratiodz, sumqdr(madr), dh
      integer :: nodeTopDisLay(madr)
      character(len=300) :: messag

      request_smaller_dt = .false.

      associate (mesh => state%mesh,         &
                 drai => state%drainage,     &
                 soil => state%soilwater,    &
                 surf => state%surfacewater, &
                 time => state%timecontrol,  &
                 sw_cfg => state%cfg%surface_water)

         select case (task)
         case (1)
            ! Init — hoisted: state%surfacewater%init(...) called directly from swap_init
            ! (spec 2026-05-13-state-init-pilot-surfacewater-design.md). Stub kept to
            ! preserve dispatcher signature; case removal is a separate cleanup.
            return

         case (2)
            ! Lateral drainage fluxes to surface water.

            ! Reset intermediate surface-water and drainage fluxes (cohort-owned).
            if (time%flZeroIntr) call surf%reset_intermediate()

            ! Reset cumulative cohorts (ADR 0042). SurfaceWater(2) is the canonical
            ! owner for both cohorts under swdra=2 (the only swdra value at which
            ! this code runs — see flSurfaceWater wiring in timecontrol).
            if (time%flZeroCumu) then
               call surf%reset_cumulative_drainage()
               call surf%reset_cumulative_reservoir()
            end if

            ! Lateral drainage.
            call bocodre(dh, state)

            ! Partition drainage flux over compartments.
            do level = 1, drai%nrlevs
               do node = 1, mesh%numnod
                  drai%qdra(level, node) = 0.0d0
               end do
            end do

            if (drai%swdivd .eq. 1) then
               call divdra(mesh%numnod, drai%nrlevs, mesh%dz, soil%ksatfit, soil%ksatexm, &
                           soil%fluseksatexm, mesh%layer, soil%cofani, soil%gwl,          &
                           drai%L, drai%qdrain, drai%qdra,                                &
                           drai%swdivdinf, drai%swnrsrf, drai%swtopnrsrf, drai%zbotdr,    &
                           time%dt, drai%FacDpthInf, drai%owltab, drai%nowltab, time%t1900)

               ! Redistribute qdrain with the new top boundary for discharge layers.
               if (drai%swdislay .eq. 2) then
                  do level = 1, drai%nrlevs
                     if (drai%swtopdislay(level) .eq. 1) then
                        zTopDisLay(level) = drai%fTopDisLay(level) * soil%gwl + &
                                            (1.0d0 - drai%fTopDisLay(level)) * (soil%gwl - dh)
                     end if
                  end do
               end if
               if (drai%swdislay .eq. 1 .or. drai%swdislay .eq. 2) then
                  do level = 1, drai%nrlevs
                     if (drai%swtopdislay(level) .eq. 1) then
                        ! Find node number of the new top of the discharge layer.
                        nodeTopDisLay(level) = 1
                        zCum = -mesh%dz(1)
                        do while (zTopDisLay(level) .lt. zCum)
                           nodeTopDisLay(level) = nodeTopDisLay(level) + 1
                           zCum = zCum - mesh%dz(nodeTopDisLay(level))
                        end do
                        ! Saturated fraction of the partial compartment.
                        difzTopDisLay(level) = zTopDisLay(level) - zCum
                        ratiodz = difzTopDisLay(level)/mesh%dz(nodeTopDisLay(level))
                        sumqdr(level) = ratiodz * drai%qdra(level, nodeTopDisLay(level))
                        do node = nodeTopDisLay(level) + 1, mesh%numnod
                           sumqdr(level) = sumqdr(level) + drai%qdra(level, node)
                        end do
                        if (dabs(sumqdr(level)) .lt. 1.0d-8) then
                           ratio = 1.0d0
                        else
                           ratio = drai%qdrain(level)/sumqdr(level)
                        end if
                        ! Redistribute drain-water fluxes.
                        do node = 1, nodeTopDisLay(level) - 1
                           drai%qdra(level, node) = 0.0d0
                        end do
                        drai%qdra(level, nodeTopDisLay(level)) = &
                           drai%qdra(level, nodeTopDisLay(level)) * ratio * ratiodz
                        do node = nodeTopDisLay(level) + 1, mesh%numnod
                           drai%qdra(level, node) = drai%qdra(level, node) * ratio
                        end do
                     end if
                  end do
               end if
            else
               ! Drainage flux through lowest compartment only.
               do level = 1, drai%nrlevs
                  drai%qdra(level, mesh%numnod) = drai%qdrain(level)
               end do
            end if

            surf%qdrtot = 0.0d0
            do level = 1, drai%nrlevs
               surf%qdrtot = surf%qdrtot + drai%qdrain(level)
            end do

            return

         case (3)
            ! Surface water balance.
            if (sw_cfg%swsrf .eq. 3) then
               surf%wlp = afgen(surf%wlptab, 2*mawlp, time%t1900 - 1.0d0 + time%dt)
            end if
            if (sw_cfg%swsec .eq. 2) then
               ! Secondary system: water level is simulated.
               call wlevbal(state, request_smaller_dt)
            elseif (sw_cfg%swsec .eq. 1) then
               ! Secondary system: water level is input.
               call wballev(state)
            end if

         case default
            call fatalerr_collected('SurfaceWater', 'Illegal value for TASK')
         end select

      end associate

      return
      end subroutine SurfaceWater

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
      ! [GR-DRA 2026-05-23] alphaw/betaw retired — aliased from state%surfacewater below.
      ! [GR-DRA 2026-05-23] QRapDra retired — read via state%drainage%QRapDra.
      ! rsro/pondmx retired (→state%surfacewater%X)
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

      if (tc_t1900-1.d0+0.1d-10 .gt. impend(imper)) goto 100

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

        if (abs(rday-1.0*intday).lt.0.00001d0 .or. tc_tcum.lt.1.0d-10) then

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
          sw_wlstar = sw_wlstar - dropr(imper)*tc_dt
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
      dvmax = (state%drainage%qdrd + state%drainage%QRapDra + wsmax) * tc_dt + state%soilwater%runots
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
        dvmax = (state%drainage%qdrd + state%drainage%QRapDra) * tc_dt + state%soilwater%runots
        swstmax = sw_swst + dvmax
        if (swstmax .le. swsttara) then

! --- apparently supply is needed for reaching target level, system
!     is made full up to level wlstara, because supply is controllable:
          wsupp = (swsttara - sw_swst - (state%drainage%qdrd+state%drainage%QRapDra)*tc_dt-state%soilwater%runots)/tc_dt
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
            wdis = (sw_swst-swsttar + (state%drainage%qdrd+state%drainage%QRapDra)*tc_dt + state%soilwater%runots )/tc_dt

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
            swstn = sw_swst + (state%drainage%qdrd + state%drainage%QRapDra - discap)*tc_dt + state%soilwater%runots
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
            swstn = sw_swst + (state%drainage%qdrd + state%drainage%QRapDra - wdisi)*tc_dt + state%soilwater%runots
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
        if(tc_dt .gt. 0.02*state%surfacewater%rsro) then
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
           if (.not.tc_fldtmin ) then
              request_smaller_dt = .true.
              fl_early_return = .true.
           else
              call dtdpst                                               &
     &        ('year-month-day,hour:minute:seconds',tc_t1900,datetime)
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
      sw_cqdrd  = sw_cqdrd  + state%drainage%qdrd*tc_dt
      sw_cwsupp = sw_cwsupp + wsupp*tc_dt
      sw_cwout  = sw_cwout  + wdis*tc_dt
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
      ! [GR-DRA 2026-05-23] QRapDra retired — read via state%drainage%QRapDra.
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
      sw_wls = AFGEN (state%surfacewater%wlstab,2*MAWLS,tc_t1900-1.d0+tc_dt)

! --- determine surface water storage for level(t-dt) and level(t)
      swstold = swstlev(state, sw_wlsold)
      sw_swst = swstlev(state, sw_wls)

! --- determine from the surface water storages and the qdrain whether
! --- supply has taken place during period (t)-(t+dt) or water has
! --- been discharged
      swstrest = swstold + (state%drainage%qdrd + state%drainage%QRapDra)*tc_dt + state%soilwater%runots - sw_swst

! --- if supply was needed, set discharge to zero
      if (swstrest.le.0.0d0) then
        wdis = 0.0d0
        wsupp = -swstrest/tc_dt

! --- if discharge has taken place, set supply to zero
      else
        wdis = swstrest/tc_dt
        wsupp = 0.0d0
      endif

! --- cumulation of water balance terms (global accumulations dropped; state aliases authoritative):
      sw_cqdrd  = sw_cqdrd  + state%drainage%qdrd*tc_dt
      sw_cwsupp = sw_cwsupp + wsupp*tc_dt
      sw_cwout  = sw_cwout  + wdis*tc_dt

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
