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
      use swap_state_mod,      only: swap_state_t
      use surfacewater_utils,  only: wlevst, swstlev, qhtab
      use swap_log,            only: log_warn
      implicit none

      type(swap_state_t), intent(inout) :: state
      logical,            intent(inout) :: request_smaller_dt

      integer :: iphase, node, intday, imper
      real(8) :: wlstx, swsttar, dvmax, swstmax, wsupp, wdis, wlstarb
      real(8) :: wover, discap, wlsl, wlsu, wlsi, swsti, wdisi, swstn
      real(8) :: wprod1, wprod2, oscil, wlstara
      real(8) :: swsttara, rday, wsmax
      character(len=200) :: messag
      character(len=19)  :: datetime
      logical :: fl_early_return

      fl_early_return = .false.

      associate (mesh => state%mesh,         &
                 drai => state%drainage,     &
                 soil => state%soilwater,    &
                 surf => state%surfacewater, &
                 time => state%timecontrol)

         surf%overfl = .false.

         ! Remember the previous target level.
         wlstarb = surf%wlstar

         ! Determine the management period the model is in.
         imper = 0
         surf%imper = 0
100      imper = imper + 1
         surf%imper = imper

         if (imper .gt. surf%nmper) then
            messag = ' sw-level oscillation at '//datetime//                &
       &           '       advise: reduction of dtmax !'
            messag = 'error sw-management periods(IMPER), more than defined'
            call fatalerr_collected('Wlevbal', messag)
         end if

         if (time%t1900 - 1.d0 + 0.1d-10 .gt. surf%impend(imper)) goto 100

         ! Determine the target surface-water level.
         if (surf%swman(imper) .eq. 1) then
            ! Fixed weir: target = weir crest (used downstream to decide outflow).
            surf%wlstar = surf%hbweir(imper)
         else
            ! Automatic weir: target depends on groundwater state.
            ! Adjust only when a new subperiod (length intwl) has started
            ! (or on the very first call).
            rday   = (time%t + 1.0d0)/surf%intwl(imper)
            intday = int(rday)

            if (abs(rday - 1.0*intday) .lt. 0.00001d0 .or. time%tcum .lt. 1.0d-10) then

               iphase = surf%nphase(imper)
               do while (soil%gwl .gt. surf%gwlcrit(imper, iphase) .and. iphase .gt. 1)
                  iphase = iphase - 1
               end do

               ! Compare total air volume with VCRIT, adapt iphase.
               surf%vtair = 0.0d0
               do node = 1, mesh%numnod
                  surf%vtair = surf%vtair + (soil%thetas(node) - soil%theta(node)) &
       &                       *abs(mesh%dz(node))
               end do
               do while (surf%vtair .lt. surf%vcrit(imper, iphase) .and. iphase .gt. 1)
                  iphase = iphase - 1
               end do

               ! Compare H(nodhd(imper)) with HCRIT, adapt iphase.
               do while (soil%h(surf%nodhd(imper)) .gt. surf%hcrit(imper, iphase) &
       &                  .and. iphase .gt. 1)
                  iphase = iphase - 1
               end do
               surf%hwlman = soil%h(surf%nodhd(imper))

               wlstx = surf%wlsman(imper, iphase)
            else
               wlstx = wlstarb
            end if

            ! If the level must drop, don't let it drop faster than dropr.
            if (wlstx .lt. surf%wlstar .and. surf%dropr(imper) .gt. 0.001d0) then
               surf%wlstar = surf%wlstar - surf%dropr(imper)*time%dt
               if (surf%wlstar .lt. wlstx) surf%wlstar = wlstx
            else
               surf%wlstar = wlstx
            end if
         end if

         ! Counter of target-level adjustments.
         if (abs(surf%wlstar - wlstarb) .gt. 0.00001d0) surf%numadj = surf%numadj + 1

         ! Storage at the target level.
         swsttar = swstlev(state, surf%wlstar)

         ! Level and storage for "max. level for supply".
         wlstara = surf%wlstar - surf%wldip(imper)
         if (wlstara .gt. (drai%zbotdr(1 + surf%nrpri) + 1.d-4)) then
            swsttara = swstlev(state, wlstara)
            wsmax    = surf%wscap(imper)
         else
            swsttara = 0.0d0
            wsmax    = 0.0d0
         end if

         ! Determine whether the system becomes full (target / weir crest).
         dvmax   = (drai%qdrd + drai%QRapDra + wsmax)*time%dt + soil%runots
         swstmax = surf%swst + dvmax

         if (swstmax .lt. 1.0d-7) then
            ! Storage decreases to zero — system falls dry; set supply to max.
            if (swstmax .lt. -0.1d0) then
               messag = 'error algorithm for sw falling dry'
               call fatalerr_collected('Wlevbal', messag)
            end if
            wsupp    = wsmax
            wdis     = 0.0d0
            surf%swst = 0.0d0
            surf%wls  = drai%zbotdr(surf%nrpri + 1)

         elseif (swstmax .ge. 0.0d0 .and. swstmax .lt. swsttara) then
            ! System won't become full — supply at max capacity.
            wsupp    = wsmax
            wdis     = 0.0d0
            surf%swst = swstmax
            surf%wls  = wlevst(state, surf%swst)
         else
            ! Determine whether supply is needed; try first without any supply.
            dvmax   = (drai%qdrd + drai%QRapDra)*time%dt + soil%runots
            swstmax = surf%swst + dvmax
            if (swstmax .le. swsttara) then
               ! Supply is needed; fill up to wlstara (supply is controllable).
               wsupp = (swsttara - surf%swst - (drai%qdrd + drai%QRapDra)*time%dt - soil%runots)/time%dt
               wdis  = 0.0d0
               surf%swst = swsttara
               surf%wls  = wlstara
            elseif (swstmax .le. swsttar) then
               ! Drainage is enough to fill above the supply target but not enough
               ! for discharge — calculate level from storage.
               wsupp    = 0.0d0
               wdis     = 0.0d0
               surf%swst = swstmax
               surf%wls  = wlevst(state, surf%swst)
            else
               ! Drainage is sufficient to reach the target level; check whether
               ! an automatic weir keeps the level at target, or the discharge
               ! relation determines the new level.
               wsupp = 0.
               if (surf%swman(imper) .eq. 2) then
                  ! Outflow = drainage flux + (target − actual) storage difference.
                  wdis = (surf%swst - swsttar + (drai%qdrd + drai%QRapDra)*time%dt + soil%runots)/time%dt

                  ! Does the weir have enough discharge capacity at this water level?
                  if (surf%swqhr .eq. 1) then
                     wover  = surf%wls - surf%hbweir(imper)
                     discap = surf%alphaw(imper) * (wover**surf%betaw(imper))
                  elseif (surf%swqhr .eq. 2) then
                     discap = qhtab(surf, surf%wlstar, imper)
                  end if
                  if (discap .gt. wdis) then
                     surf%wls    = surf%wlstar
                     surf%swst   = swsttar
                     surf%overfl = .false.
                  else
                     surf%overfl = .true.
                  end if
               end if
               if (surf%swman(imper) .eq. 1 .or. surf%overfl) then
                  ! Determine the level from the q-h relationship and account
                  ! for the change in storage. First check the system doesn't overflow.
                  if (surf%swqhr .eq. 1) then
                     wover  = surf%sttab(1, 1) - surf%hbweir(imper)
                     discap = surf%alphaw(imper) * (wover**surf%betaw(imper))
                  elseif (surf%swqhr .eq. 2) then
                     discap = surf%qqhtab(imper, 1)
                  end if

                  swstn = surf%swst + (drai%qdrd + drai%QRapDra - discap)*time%dt + soil%runots
                  if (swstn .gt. surf%sttab(1, 2)) then
                     messag = 'surface water system has overflowed!'
                     call fatalerr_collected('Wlevbal', messag)
                  end if

                  ! Bisection iteration for new level, storage and discharge.
                  wlsl = surf%hbweir(imper)
                  wlsu = surf%sttab(1, 1)

700               wlsi  = (wlsl + wlsu) * 0.5
                  swsti = swstlev(state, wlsi)
                  if (surf%swqhr .eq. 1) then
                     wdisi = surf%alphaw(imper)*(wlsi - surf%hbweir(imper))**surf%betaw(imper)
                  else
                     wdisi = qhtab(surf, wlsi, imper)
                  end if
                  swstn = surf%swst + (drai%qdrd + drai%QRapDra - wdisi)*time%dt + soil%runots
                  if (swstn .lt. swsti) then
                     wlsu = wlsi
                  else
                     wlsl = wlsi
                  end if
                  if ((wlsu - wlsl) .gt. 0.001d0) then
                     goto 700
                  else
                     surf%wls  = wlsi
                     surf%swst = swstn
                     wdis      = wdisi
                  end if
               end if
            end if
         end if

         ! Ponding in extended drainage may limit the timestep.
         if (surf%wls .gt. surf%pondmx .or. soil%pond .gt. surf%pondmx) then
            if (time%dt .gt. 0.02*surf%rsro) request_smaller_dt = .true.
            fl_early_return = .true.
         end if

         if (.not. fl_early_return) then
            ! Update last-four-levels ring buffer to detect oscillation.
            surf%wlsbak(1) = surf%wlsbak(2)
            surf%wlsbak(2) = surf%wlsbak(3)
            surf%wlsbak(3) = surf%wlsbak(4)
            surf%wlsbak(4) = surf%wls
            wprod1 = (surf%wlsbak(2) - surf%wlsbak(1))*(surf%wlsbak(3) - surf%wlsbak(2))
            wprod2 = (surf%wlsbak(3) - surf%wlsbak(2))*(surf%wlsbak(4) - surf%wlsbak(3))
            if (wprod1 .lt. 0.0d0 .and. wprod2 .lt. 0.0d0) then
               oscil = abs(surf%wlsbak(3) - surf%wlsbak(2))
               if (oscil .gt. surf%osswlm) then
                  if (.not. time%fldtmin) then
                     request_smaller_dt = .true.
                     fl_early_return    = .true.
                  else
                     call dtdpst                                            &
       &                ('year-month-day,hour:minute:seconds', time%t1900, datetime)
                     messag = ' sw-level oscillation at '//datetime//       &
       &                      '       advise: reduction of dtmax !'
                     call log_warn('Wlevbal', messag)
                     call fatalerr_collected('Wlevbal', messag)
                  end if
               end if
            end if
         end if

         if (.not. fl_early_return) then
            ! Cumulative terms.
            surf%cqdrd  = surf%cqdrd  + drai%qdrd*time%dt
            surf%cwsupp = surf%cwsupp + wsupp*time%dt
            surf%cwout  = surf%cwout  + wdis*time%dt
         end if

      end associate

      return
      end subroutine WLEVBAL

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
      use swap_state_mod,        only: swap_state_t
      use array_utils,           only: afgen
      use surfacewater_utils,    only: swstlev
      use swap_array_dimensions, only: mawls
      implicit none

      type(swap_state_t), intent(inout) :: state

      real(8) :: swstold, swstrest, wdis, wsupp

      associate (drai => state%drainage,    &
                 soil => state%soilwater,   &
                 surf => state%surfacewater, &
                 time => state%timecontrol)

         ! Memorise the previous-step water level.
         surf%wlsold = surf%wls

         ! Fetch new level from the input series.
         surf%wls = afgen(surf%wlstab, 2*MAWLS, time%t1900 - 1.d0 + time%dt)

         ! Storage for level(t-dt) and level(t).
         swstold  = swstlev(state, surf%wlsold)
         surf%swst = swstlev(state, surf%wls)

         ! Decide whether supply or discharge took place over the timestep.
         swstrest = swstold + (drai%qdrd + drai%QRapDra)*time%dt + soil%runots - surf%swst

         if (swstrest .le. 0.0d0) then
            wdis  = 0.0d0
            wsupp = -swstrest/time%dt
         else
            wdis  = swstrest/time%dt
            wsupp = 0.0d0
         end if

         ! Cumulative balance terms.
         surf%cqdrd  = surf%cqdrd  + drai%qdrd*time%dt
         surf%cwsupp = surf%cwsupp + wsupp*time%dt
         surf%cwout  = surf%cwout  + wdis*time%dt

      end associate

      return
      end subroutine WBALLEV


!> Year-boundary reset for surface-water cumulative state.
!! Captures the current swst as swstini for the new year.
subroutine surfacewater_year_reset(sw)
   use surfacewater_state_mod, only: surfacewater_state_t
   implicit none
   type(surfacewater_state_t), intent(inout) :: sw
   sw%swstini = sw%swst
end subroutine surfacewater_year_reset

end module surfacewater_mod
