!> @brief Soil water balance calculations module
!>
!> @details
!> This module contains subroutines and functions for water balance and flux calculations
!> in the soil profile, including groundwater level determination, water table searches,
!> flux calculations between compartments, and mass balance checking.
!>
!> Public procedures:
!> - calcgwl: Calculate groundwater level
!> - level: Calculate water level from pressure head
!> - watertable: Search for watertable and perched watertable
!> - fluxes: Calculate fluxes between compartments
!> - integral: Calculate intermediate and cumulative fluxes
!> - checkmassbal: Check mass balance per output period
!> - watstor: Calculate water storage in soil profile
module soilwaterbalance_mod
   use error_mod, only: fatalerr_collected
   use swap_state_mod, only: swap_state_t
    implicit none
    private
    ! `watertable` extracted to src/soil/dormant/watertable.f90 (no live callers).
    ! `checkmassbal` extracted to src/soil/dormant/checkmassbal.f90 (no live callers).
    ! `level` retained — calcgwl calls it.
    public :: calcgwl, level, fluxes, integral, watstor
contains
      !> @brief Calculate groundwater level
      !>
      !> @details
      !> Searches for the watertable and perched watertable (if existing) in the soil profile.
      !> The groundwater level is determined based on pressure heads in the soil compartments.
      !> For flux calculations, the profile must extend below the groundwater level.
      !>
      !> @note
      !> Date: July 2002, updated April 2008
      !>
      !> Update: An unsaturated zone embedded in a saturated soil column should contain
      !> at least a total of 'CritAir' cm of air to be recognized as really unsaturated.
      !>
      !> SAVE statement removed - all local variables are reset at start of each call
      !> (legacy code that was unnecessary)
      !> @endnote
      subroutine calcgwl (state)
      use swap_log, only: log_debug, log_warn, to_str
      implicit none
      type(swap_state_t), intent(inout) :: state
      ! local
      integer   i, node, nodhlp, nodheq1
      integer   nodgwlflcpzo_loc
      real(8)   gwlflcpzo_loc
      logical   flsat,flunsat
      character(len=200) messag
      character(len=19) datexti

      associate (mesh => state%mesh, soil => state%soilwater, time => state%timecontrol)

      ! set initial values
      soil%gwl    = 999.0d0
      soil%pegwl  = 999.0d0
      flsat     = .false.
      soil%nodgwl = mesh%numnod+1
      nodhlp    = mesh%numnod
      nodheq1   = mesh%numnod

      ! search for groundwater table
      if (soil%h(mesh%numnod).ge.0.0d0) flsat  = .true.

      node = mesh%numnod
      nodgwlflcpzo_loc = mesh%numnod + 1
      soil%nodgwlflcpzo = mesh%numnod + 1
      gwlflcpzo_loc = soil%gwl
      soil%gwlflcpzo    = soil%gwl
      do while (flsat .and. node.gt.1)
         node = node - 1
         if(soil%swbotb_runtime.eq.1)then
            if (soil%h(node) .lt. 0.0d0) then
               soil%gwl = mesh%z(node+1) + soil%h(node+1) / (soil%h(node+1)-soil%h(node)) * mesh%disnod(node+1)
               flsat   =.false.
               soil%nodgwl = node
            endif
         else
            if (soil%h(node) .lt. 1.0d0 .and. nodheq1.eq.mesh%numnod) nodheq1 = node

            if (soil%h(node) .lt. 0.0d0) then
               ! Macropore branch removed (ADR 0040; always .false.).
               flsat  = .false.
               soil%nodgwl = node
               soil%gwl    = level (state,1,node,nodheq1)
            endif
         endif
      end do

      ! whole profile saturated, then add ponding layer to groundwater level
      if (flsat)then
         if(soil%h(1) .gt. 0.0d0)then
            if (soil%pond .lt. 1.d-8) then
               soil%gwl = min(mesh%z(1)+soil%h(1),soil%pond)
            else
               soil%gwl = soil%pond
            endif
         else
            soil%gwl = 0.0d0
         end if
         soil%nodgwl = 1
         ! if (flmacropore) block dropped (ADR 0040; always .false.).
      endif

      ! search for perched groundwater table

      ! first, search for first saturated compartment (i) above groundwater level
      i = nodhlp
      flunsat = .true.
      do while (flunsat .and. i.ge.1)
         if (soil%h(i).ge.0.0d0) flunsat = .false.
         i = i - 1
      enddo

      ! if saturated compartment above gwl exists, then find perched groundwater table
      if (i.ne.0) then
         flsat  = .true.
         soil%bpegwl = i
         node   = soil%bpegwl
         nodheq1 = soil%bpegwl

         do while (flsat .and. node.gt.1)
            node = node - 1

            if (soil%h(node) .lt. 1.0d0 .and. nodheq1.eq.soil%bpegwl) nodheq1 = node

            if (soil%h(node) .lt. 0.0d0) then
               ! Macropore branch removed (ADR 0040; always .false.).
               flsat = .false.
               soil%npegwl = node
               soil%pegwl  = level (state,1,node,nodheq1)
            endif
         end do

         ! whole profile saturated, then add ponding layer to perched groundwater level
         if (flsat)then
            if(soil%h(1) .gt. 0.0d0)then
               if (soil%pond .lt. 1.d-8) then
                  soil%pegwl = min(mesh%z(1)+soil%h(1),soil%pond)
               else
                  soil%pegwl = soil%pond
               endif
            else
               soil%pegwl = 0.0d0
            end if
            soil%npegwl = 1
         endif
      else
         soil%bpegwl = -1
         soil%npegwl = -1
      endif

      ! fatal error if gwl below profile and flux has to be calculated
      if ((soil%swbotb_runtime.eq.3.or.soil%swbotb_runtime.eq.4)  &
     &    .and.soil%gwl.gt.998.0d0) then
          messag = 'The groundwater level descends below the lower' &
     &     //' boundary. This conflicts with bottom boundary' &
     &     //' condition 3 and 4. Extend soil profile!'
         call fatalerr_collected ('calcgwl',messag)
      endif

      ! warning error if there is inconsistency between defined gwl and soil physics
      if (soil%swbotb_runtime.eq.1 .and.  &
     &    (soil%gwlinp .ge.mesh%z(1) .or. soil%gwl.gt.998.0d0)) then
         call dtdpst('year-month-day,hour:minute:seconds',time%t1900,datexti)
         write(messag,'(6a)')                                           &
     &         'No groundwater level because unsaturation at bottom ',  &
     &         'compartment ( ', datexti,  ' ). ',                      &
     &         'This is caused by inconsistency between ',              &
     &         'given gwl and soil physical parameters '
         call log_warn('Calcgwl', messag)
      endif

      end associate

      return
      end subroutine calcgwl


      !> @brief Calculate water level from pressure head
      !>
      !> @details
      !> Method 1 (swoptlev=1): GWL = elevation head where h = 0.
      !> Method 2 (swoptlev=2): GWL = average of elevation heads at h = -1 and h = +1.
      function level (state,swoptlev,node,nodheq1)
      implicit none

      type(swap_state_t), intent(in) :: state
      integer node, nodheq1, swoptlev
      integer i
      real(8) levm1, levp1
      real(8) level

      associate (mesh => state%mesh, soil => state%soilwater)

      if (swoptlev.eq.1) then
         if (soil%h(node+1).ge.0.0d0)then
            level = mesh%z(node+1) + soil%h(node+1) / (soil%h(node+1)-soil%h(node)) * mesh%disnod(node+1)
         else
            level = mesh%zbotcp(node) - soil%h(node)
            level = min(mesh%z(node),max(mesh%zbotcp(node),level))
         end if

      elseif (swoptlev.eq.2) then
         i = nodheq1
         if (nodheq1.eq.mesh%numnod) then
            levp1 = mesh%z(i) - 0.5d0 * mesh%dz(i)
         else
            levp1 = mesh%z(i) - (mesh%z(i) - mesh%z(i+1)) * (1.d0-soil%h(i)) / (soil%h(i+1)-soil%h(i))
         endif
         i = node
         do while (soil%h(i).gt.-1.d0 .and. i.gt.1)
            i = i - 1
         enddo
         if (i.eq.1 .and. soil%h(1).gt.-1.d0 .and. node.gt.2) then
            levm1 = mesh%z(node+1) + soil%h(node+1) / (soil%h(node+1)-soil%h(node)) * mesh%disnod(node+1)
            levp1 = levm1
         else
            levm1 = mesh%z(i+1) + (mesh%z(i) - mesh%z(i+1)) * (1.d0+soil%h(i+1)) / (soil%h(i+1)-soil%h(i))
         endif
         level = (levp1 + levm1) / 2.d0
      endif

      end associate

      return
    end function level

      ! `watertable` extracted to src/soil/dormant/watertable.f90
      ! (no live callers; CritUndSatVol path was never ported to the TOML pipeline).
      ! See dormant module header for reactivation checklist.

      !> @brief Calculate fluxes between compartments
      !>
      !> @details
      !> Calculates the water fluxes between soil compartments based on water content
      !> changes, root extraction, drainage, and boundary conditions. Fluxes are computed
      !> from volume changes per compartment and accumulated for mass balance tracking.
      !>
      !> @note
      !> Date: 29/9/99
      !> @endnote
      subroutine fluxes (state)
      ! qimmob/QExcMpMtx/QMaPo inlined as 0 (fingered-flow + ADR 0040 macropore terms).
      use swap_state_mod, only: swap_state_t
      implicit none

      type(swap_state_t), intent(inout) :: state
      integer i,level

      associate (mesh => state%mesh,         &
                 soil => state%soilwater,    &
                 drai => state%drainage,     &
                 surf => state%surfacewater, &
                 time => state%timecontrol)

      ! determine qbot if not specified
      if (soil%swbotb_runtime .eq. 5 .or. soil%swbotb_runtime .eq. 7 .or.                         &
     &    soil%swbotb_runtime .eq. 8 .or. soil%swbotb_runtime .eq. -2 .or.                        &
     &    (soil%swbotb_runtime .eq. 1 .and. soil%fllowgwl)) then
        ! QMaPo retired-zero inlined (ADR 0040 macropore outflow term).
        soil%qbot = soil%qtop + soil%qrosum + surf%qdrtot + (soil%volact-soil%volm1)/time%dt - soil%qssdisum
      endif

      ! calculate fluxes (cm/d) from changes in volume per compartment
      i = mesh%numnod+1
      soil%q(i)              = soil%qbot
      soil%inq(i)            = soil%inq(i) + soil%q(i)*time%dt
      do i = mesh%numnod,1,-1
        ! qimmob (fingered-flow flux) and QExcMpMtx (ADR 0040 macropore exchange) retired-zero inlined.
        soil%q(i) = - (soil%theta(i)-soil%thetm1(i))*soil%FrArMtrx(i)*mesh%dz(i)/time%dt +  &
     &                soil%q(i+1)-soil%qrot(i)+soil%qssdi(i)

        if (allocated(drai%qdra)) then
          do level=1,drai%nrlevs
             soil%q(i) = soil%q(i) - drai%qdra(level,i)
          enddo
        end if
        soil%inq(i) = soil%inq(i) + soil%q(i)*time%dt
      end do

      end associate

      return
      end

      !> @brief Calculate intermediate and cumulative fluxes
      !>
      !> @details
      !> Calculates and accumulates water fluxes over timesteps, including:
      !> - Root extraction (actual and potential transpiration)
      !> - Soil evaporation (potential and reduced)
      !> - Drainage fluxes at multiple levels
      !> - Bottom boundary fluxes
      !> - Interception, precipitation, runoff, and runon
      !> - Computes both intermediate totals and cumulative values
      !> - Tracks water balance errors for compensation
      !>
      !> @note
      !> Date: November 2004
      !> @endnote
      subroutine integral (state)
      ! qdraincomp dropped (per-comp write, no reader; orphan).
      use iso_fortran_env, only: real64
      implicit none

      type(swap_state_t), intent(inout) :: state
      integer node,level
      real(8) qrotts,qdrats,ptrats,pevats,revats,qbotts

      associate (mesh => state%mesh,         &
                 soil => state%soilwater,    &
                 drai => state%drainage,     &
                 surf => state%surfacewater, &
                 atmo => state%atmosphere,   &
                 crop => state%crop,         &
                 time => state%timecontrol)

      if (time%flZeroIntr) then
      endif

      ! potential transpiration of this timestep
      ptrats = atmo%ptra * time%dt

      ! potential soil evaporation of this timestep
      pevats = atmo%peva * time%dt

      ! reduced soil evaporation of this timestep
      revats = soil%reva * time%dt

      ! flux lower boundary of this timestep
      qbotts = soil%qbot*time%dt

      ! total root extraction of this timestep
      qrotts = soil%qrosum * time%dt

      ! total drainage flux of this timestep
      qdrats = surf%qdrtot * time%dt

      ! determine daily actual transpiration
      if (time%flDayStart) soil%tra = 0.0_real64
      soil%tra = soil%tra + qrotts

      ! add time step fluxes to intermediate totals
      soil%iqrot = soil%iqrot + qrotts
      do node = 1, crop%common%noddrz
        soil%inqrot(node) = soil%inqrot(node) + soil%qrot(node) * time%dt
        soil%qpotrot_day(node) = soil%qpotrot_day(node) + soil%qpotrot(node) * time%dt
        soil%qredtot_day(node) = soil%qredtot_day(node) + (soil%qredwet(node) + soil%qreddry(node) + soil%qredsol(node) + soil%qredfrs(node)) * time%dt
      end do
      do node = 1,mesh%numnod
        soil%inqssdi(node) = soil%inqssdi(node) + soil%qssdi(node) * time%dt
        soil%iqssdi = soil%iqssdi + soil%qssdi(node) * time%dt
      end do
      soil%iqredwet = soil%iqredwet + soil%qredwetsum*time%dt
      soil%iqreddry = soil%iqreddry + soil%qreddrysum*time%dt
      soil%iqredsol = soil%iqredsol + soil%qredsolsum*time%dt
      soil%iqredfrs = soil%iqredfrs + soil%qredfrssum*time%dt
      soil%iqredwet_day = soil%iqredwet_day + soil%qredwetsum*time%dt
      soil%iqreddry_day = soil%iqreddry_day + soil%qreddrysum*time%dt
      soil%iqredsol_day = soil%iqredsol_day + soil%qredsolsum*time%dt
      soil%iqredfrs_day = soil%iqredfrs_day + soil%qredfrssum*time%dt
      soil%iptra_day = soil%iptra_day + atmo%ptra * time%dt
      soil%ies0 = soil%ies0 + 0.1d0*crop%es0*time%dt
      soil%iet0 = soil%iet0 + 0.1d0*crop%et0*time%dt
      soil%iew0 = soil%iew0 + 0.1d0*crop%ew0*time%dt

      surf%iqdra = surf%iqdra + qdrats + drai%QRapDra*time%dt
      ! per-comp qdraincomp write dropped (no reader anywhere in src/).
      do node = 1,mesh%numnod
        do level = 1,drai%nrlevs
          if (allocated(surf%inqdra) .and. allocated(drai%qdra)) then
            surf%inqdra(level,node) = surf%inqdra(level,node) + drai%qdra(level,node)*time%dt
            if (drai%qdra(level,node) > 0.0d0) then
               surf%inqdra_out(level,node) = surf%inqdra_out(level,node) + drai%qdra(level,node)*time%dt
            else
               surf%inqdra_in(level,node)  = surf%inqdra_in(level,node) - drai%qdra(level,node)*time%dt
            end if
          end if
        end do
      end do

      soil%iintc = soil%iintc + (atmo%aintcdt+crop%gird-atmo%nird)*time%dt

      atmo%intr%iptra = atmo%intr%iptra + ptrats
      atmo%intr%ipeva = atmo%intr%ipeva + pevats
      atmo%intr%ievap = atmo%intr%ievap + revats
      soil%iruno = soil%iruno + soil%runots
      soil%irunon = soil%irunon + soil%runon*time%dt
      soil%iprec = soil%iprec + (atmo%graidt+crop%gird)*time%dt
      atmo%intr%igrai = atmo%intr%igrai + atmo%graidt*time%dt
      soil%igird = soil%igird + crop%gird*time%dt
      atmo%intr%inrai = atmo%intr%inrai + atmo%nraidt*time%dt
      soil%inird = soil%inird + atmo%nird*time%dt
      soil%iqbot = soil%iqbot + qbotts
      if (soil%q(1) < 0.0d0) then
         soil%iqtdo = soil%iqtdo - soil%q(1)*time%dt
      else
         soil%iqtup = soil%iqtup + soil%q(1)*time%dt
      end if
      do node = 1, mesh%numnod+1
         if (soil%q(node) < 0.0d0) then
            soil%iqdo(node) = soil%iqdo(node) - soil%q(node)*time%dt
         else
            soil%iqup(node) = soil%iqup(node) + soil%q(node)*time%dt
         end if
      end do

      ! add time step fluxes to total cumulative values
      soil%cqssdi = soil%cqssdi + soil%qssdisum*time%dt
      soil%cqrot  = soil%cqrot  + qrotts
      surf%cqdra = surf%cqdra + qdrats
      atmo%cumu%cptra = atmo%cumu%cptra + ptrats
      atmo%cumu%cpeva = atmo%cumu%cpeva + pevats
      atmo%cumu%cevap = atmo%cumu%cevap + revats
      if (soil%runots.lt.0.0d0) then
        soil%cinund = soil%cinund - soil%runots
      else if (soil%runots.gt.0.0d0) then
        soil%crunoff = soil%crunoff + soil%runots
      endif
      soil%irunoCN   = soil%irunoCN   + atmo%Runoff_CN*time%dt
      soil%crunoffCN = soil%crunoffCN + atmo%Runoff_CN*time%dt

      atmo%cumu%caintc = atmo%cumu%caintc + (atmo%aintcdt+crop%gird-atmo%nird)*time%dt

      atmo%cumu%cgrai = atmo%cumu%cgrai + atmo%graidt*time%dt
      atmo%cumu%cnrai = atmo%cumu%cnrai + atmo%nraidt*time%dt
!      cnrai = cgrai - caintc
      soil%cgird = soil%cgird + crop%gird*time%dt
      soil%cnird = soil%cnird + atmo%nird*time%dt

      if (qbotts.lt.0.0d0) then
        soil%cqbotdo = soil%cqbotdo - qbotts
      else if (qbotts.gt.0.0d0) then
        soil%cqbotup = soil%cqbotup + qbotts
      endif
      soil%cqbot = soil%cqbot + qbotts
      if (allocated(surf%cqdrain)) then
        do level = 1,drai%nrlevs
          ! infiltration
          if (drai%qdrain(level).lt.0.0d0) then
            surf%cqdrainin(level) = surf%cqdrainin(level) - drai%qdrain(level)*time%dt
          ! drainage
          else if (drai%qdrain(level).gt.0.0d0) then
            surf%cqdrainout(level) = surf%cqdrainout(level) + drai%qdrain(level)*time%dt
          endif
          surf%cqdrain(level) = surf%cqdrain(level) + drai%qdrain(level)*time%dt
        enddo
      end if

      ! rain on the ponding surface
      soil%cqprai = soil%cqprai + atmo%nraidt*time%dt
      soil%crunon = soil%crunon + soil%runon*time%dt
      if (soil%q(1).lt.0.0d0) then
        soil%cqtdo = soil%cqtdo - soil%q(1)*time%dt
      else if (soil%q(1).gt.0.0d0) then
        soil%cqtup = soil%cqtup + soil%q(1)*time%dt
      endif

      ! compensate water balance error of this time step during remaining day part
      ! cumulative water balance error
      if (state%cfg%meteo%snow%swsnow.eq.0) then
        soil%wbalance = atmo%cumu%cnrai + soil%cnird           &
     &        + soil%crunon - soil%crunoff                             &
     &        - soil%cqrot - atmo%cumu%cevap                              &
     &        - surf%cqdra                                          &
     &        + soil%cqbot + soil%volini                                   &
     &        - soil%volact + soil%pondini                                      &
     &        - soil%pond + soil%cqssdi
      else
         soil%wbalance = soil%cqprai + soil%cnird          &
     &        + atmo%cumu%cmelt                                                           &
     &        + soil%crunon - soil%crunoff                            &
     &        - soil%cqrot - atmo%cumu%cevap                              &
     &        - surf%cqdra                                          &
     &        + soil%cqbot + soil%volini                                   &
     &        - soil%volact + soil%pondini                                      &
     &        - soil%pond + soil%cqssdi
      endif

      ! if (FlMacropore) wbalance block dropped (ADR 0040; always .false.).

      end associate

      return
      end

      ! `checkmassbal` extracted to
      ! src/soil/dormant/checkmassbal.f90 (no live callers; ANIMO/PEARL
      ! coupling dispatch site was never ported to TOML pipeline).
      ! Globals it held hostage (NumNodNew, outfil, pathwork, DZNew,
      ! CritDevMasBal, dev_cmb) now have no other consumer in this
      ! module and remain in variables.f90 only for other call sites.

      !> @brief Calculate water storage in soil profile
      !>
      !> @details
      !> Calculates the total water storage in the soil profile by summing
      !> water content over all compartments, accounting for compartment thickness
      !> and matrix fraction. Updates both previous and current storage values.
      !>
      !> @note
      !> Date: 29/9/99
      !>
      !> Differences SWAP/SWAPS: SWAPS has extra parameters
      !> @endnote
      subroutine watstor (state)
      use swap_state_mod, only: swap_state_t
      IMPLICIT NONE

      type(swap_state_t), intent(inout) :: state
      INTEGER i

      associate (mesh => state%mesh, soil => state%soilwater)

      ! update soil profile water storage
      soil%volm1  = soil%volact
      soil%volact = 0.0d0
      do 10 i = 1,mesh%numnod
        soil%volact = soil%volact + soil%theta(i)*mesh%dz(i)*soil%FrArMtrx(i)
 10   continue

      end associate

      return
      end

end module