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
    public :: calcgwl, level, watertable, fluxes, integral, checkmassbal, watstor
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
      ! All bare-global imports retired — calcgwl reads/writes through the
      ! sub-record associate (mesh/soil/time) below.
      use swap_log, only: log_debug, log_warn, to_str
      implicit none
      ! SS-BND B-2.7: state added to read soil%gwlinp (gwlinp global retired)
      ! SS-SWC S-1.6: intent(in) -> intent(inout) to allow gwl/nodgwl/pegwl/bpegwl/npegwl/gwlflcpzo/nodgwlflcpzo dual-writes
      ! SS-SWC S-2.4: h/pond reads cut over to soil; level/watertable pass state
      type(swap_state_t), intent(inout) :: state
      ! local
      integer   i, node, nodhlp, nodheq1
      integer   nodgwlflcpzo_loc
      real(8)   gwlflcpzo_loc
      logical   flsat,flunsat
      character(len=200) messag
      character(len=19) datexti

      ! S-2.4 ASSOCIATE: alias soil arrays for h/pond reads
      ! [GR-BH C4] mesh globals aliased via mesh
      ! [GR-BH Audit 31] soil%swbotb_runtime aliased via soil%swbotb_runtime
      associate (mesh => state%mesh, soil => state%soilwater, time => state%timecontrol)

      ! set initial values — [SS-SWC S-2.12B] legacy half-writes dropped
      soil%gwl    = 999.0d0                                     ! S-1.6/S-2.12B
      soil%pegwl  = 999.0d0                        ! S-1.6/S-2.12B
      flsat     = .false.
      soil%nodgwl = mesh%numnod+1                       ! S-1.6/S-2.12B
      nodhlp    = mesh%numnod
      nodheq1   = mesh%numnod

      ! search for groundwater table
      if (soil%h(mesh%numnod).ge.0.0d0) flsat  = .true.      ! S-2.4 read cutover: h -> soil%h

      node = mesh%numnod
      nodgwlflcpzo_loc = mesh%numnod + 1
      soil%nodgwlflcpzo = mesh%numnod + 1               ! S-1.6/S-2.12B
      gwlflcpzo_loc = soil%gwl
      soil%gwlflcpzo    = soil%gwl                   ! S-1.6/S-2.12B
      do while (flsat .and. node.gt.1)
         node = node - 1
         if(soil%swbotb_runtime.eq.1)then
            if (soil%h(node) .lt. 0.0d0) then                  ! S-2.4 read cutover
               soil%gwl = mesh%z(node+1) + soil%h(node+1) / (soil%h(node+1)-soil%h(node)) * mesh%disnod(node+1)  ! S-2.4/S-2.12B
               flsat   =.false.
               soil%nodgwl = node                  ! S-1.6/S-2.12B
            endif
         else
            if (soil%h(node) .lt. 1.0d0 .and. nodheq1.eq.mesh%numnod) nodheq1 = node  ! S-2.4

            if (soil%h(node) .lt. 0.0d0) then                  ! S-2.4 read cutover
               ! [SS-GR-CROPRT A2] flmacropore branch removed (ADR 0040; always .false.)
               flsat  = .false.
               soil%nodgwl = node                ! S-1.6/S-2.12B
               soil%gwl    = level (state,1,node,nodheq1)     ! S-2.4/S-2.12B
            endif
         endif
      end do

      ! whole profile saturated, then add ponding layer to groundwater level
      if (flsat)then
         if(soil%h(1) .gt. 0.0d0)then                          ! S-2.4 read cutover
            if (soil%pond .lt. 1.d-8) then                     ! S-2.4 read cutover: pond -> soil%pond
               soil%gwl = min(mesh%z(1)+soil%h(1),soil%pond)            ! S-2.4/S-2.12B
            else
               soil%gwl = soil%pond                              ! S-2.4/S-2.12B
            endif
         else
            soil%gwl = 0.0d0
         end if
         soil%nodgwl = 1                           ! S-1.6/S-2.12B
         ! [SS-GR-CROPRT A2] if (flmacropore) block dropped (ADR 0040; always .false.)
      endif

      ! search for perched groundwater table

      ! first, search for first saturated compartment (i) above groundwater level
      i = nodhlp
      flunsat = .true.
      do while (flunsat .and. i.ge.1)
         if (soil%h(i).ge.0.0d0) flunsat = .false.             ! S-2.4 read cutover
         i = i - 1
      enddo

      ! if saturated compartment above gwl exists, then find perched groundwater table
      if (i.ne.0) then
         flsat  = .true.
         soil%bpegwl = i                           ! S-1.6/S-2.12B
         node   = soil%bpegwl
         nodheq1 = soil%bpegwl

         do while (flsat .and. node.gt.1)
            node = node - 1

            if (soil%h(node) .lt. 1.0d0 .and. nodheq1.eq.soil%bpegwl) nodheq1 = node  ! S-2.4

            if (soil%h(node) .lt. 0.0d0) then                  ! S-2.4 read cutover
               ! [SS-GR-CROPRT A2] flmacropore branch removed (ADR 0040; always .false.)
               flsat = .false.
               soil%npegwl = node                ! S-1.6/S-2.12B
               soil%pegwl  = level (state,1,node,nodheq1)  ! S-2.4/S-2.12B
            endif
         end do

         ! whole profile saturated, then add ponding layer to perched groundwater level
         if (flsat)then
            if(soil%h(1) .gt. 0.0d0)then                       ! S-2.4 read cutover
               if (soil%pond .lt. 1.d-8) then                  ! S-2.4 read cutover
                  soil%pegwl = min(mesh%z(1)+soil%h(1),soil%pond)  ! S-2.4/S-2.12B
               else
                  soil%pegwl = soil%pond            ! S-2.4/S-2.12B
               endif
            else
               soil%pegwl = 0.0d0                  ! S-2.12B
            end if
            soil%npegwl = 1                        ! S-1.6/S-2.12B
         endif
      else
         soil%bpegwl = -1                          ! S-1.6/S-2.12B
         soil%npegwl = -1                          ! S-1.6/S-2.12B
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
      !> Calculates the water level (elevation head) from pressure head using one of two methods:
      !> - Method 1 (swoptlev=1): Groundwater level equals elevation head where h = 0
      !> - Method 2 (swoptlev=2): Groundwater level equals average of elevation heads at h = -1 and h = +1
      !>
      !> @param[in] swoptlev Option for level calculation method (1 or 2)
      !> @param[in] node Node number for calculation
      !> @param[in] nodheq1 Node number where h equals 1
      !> @return Water level (elevation head) in cm
      !>
      !> @note
      !> Date: April 2008
      !> @endnote
      ! SS-SWC S-2.4: state added as first arg so h reads come from soil%h
      function level (state,swoptlev,node,nodheq1)
      ! [SS-SWC S-2.12B] h retired from use clause; read via soil%h (associate below)
      ! [GR-BH C4] mesh%numnod/mesh%disnod/mesh%dz/mesh%z/mesh%zbotcp migrated to mesh; use variables no longer needed here
      implicit none

      type(swap_state_t), intent(in) :: state
      integer node, nodheq1, swoptlev
      integer i
      real(8) levm1, levp1
      real(8) level

      ! S-2.4 ASSOCIATE: alias soil%h for reads inside level
      ! [GR-BH C4] mesh globals aliased via mesh
      associate (mesh => state%mesh, soil => state%soilwater)

      if (swoptlev.eq.1) then
         ! groundwater level equals elevation head where h = 0
         if (soil%h(node+1).ge.0.0d0)then                     ! S-2.4 read cutover
            level = mesh%z(node+1) + soil%h(node+1) / (soil%h(node+1)-soil%h(node)) * mesh%disnod(node+1)  ! S-2.4
         else
            level = mesh%zbotcp(node) - soil%h(node)               ! S-2.4
            level = min(mesh%z(node),max(mesh%zbotcp(node),level))
         end if

      elseif (swoptlev.eq.2) then
         ! groundwater level equals average of elevation heads of h = -1 and h = +1
         ! elevation head of h = +1
         i = nodheq1
         if (nodheq1.eq.mesh%numnod) then
            levp1 = mesh%z(i) - 0.5d0 * mesh%dz(i)
         else
            levp1 = mesh%z(i) - (mesh%z(i) - mesh%z(i+1)) * (1.d0-soil%h(i)) / (soil%h(i+1)-soil%h(i))  ! S-2.4
         endif
         ! elevation head of h = -1
         i = node
         do while (soil%h(i).gt.-1.d0 .and. i.gt.1)          ! S-2.4 read cutover
            i = i - 1
         enddo
         if (i.eq.1 .and. soil%h(1).gt.-1.d0 .and. node.gt.2) then   ! S-2.4
            ! no compartment with pressure head < -1 cm in top of profile:
            ! use elevation head of h = 0 as estimation for groundwater level
            levm1 = mesh%z(node+1) + soil%h(node+1) / (soil%h(node+1)-soil%h(node)) * mesh%disnod(node+1)  ! S-2.4
            levp1 = levm1
         else
            levm1 = mesh%z(i+1) + (mesh%z(i) - mesh%z(i+1)) * (1.d0+soil%h(i+1)) / (soil%h(i+1)-soil%h(i))  ! S-2.4
         endif
         ! groundwater level = average of levp1 and levm1
         level = (levp1 + levm1) / 2.d0
      endif

      end associate  ! soil%h (S-2.4); mesh%numnod/mesh%disnod/mesh%dz/mesh%z/mesh%zbotcp [GR-BH C4]

      return
    end function level

      !> @brief Search for watertable and perched watertable
      !>
      !> @details
      !> Searches for the watertable and perched watertable (if existing) using
      !> a criterion based on total unsaturated volume. An unsaturated zone embedded
      !> in saturated soil must contain at least CritUndSatVol cm of air to be recognized
      !> as truly unsaturated.
      !>
      !> @param[inout] node Starting node for search
      !> @param[out] nodlev Node containing water level
      !> @param[out] nodhlp Deepest unsaturated node
      !> @param[in] nodheq1 Node where h equals 1
      !> @param[in] CritUndSatVol Critical unsaturated volume threshold (cm)
      !> @param[inout] flsat Saturation flag
      !> @param[out] waterlevel Calculated water level (cm)
      !>
      !> @note
      !> Date: April 2008
      !>
      !> NODGWL is NOT node with GWL, but DEEPEST UNSATURATED NODE
      !> @endnote
      ! SS-SWC S-2.4: state added as first arg so h/Theta/ThetaS reads come from soil
      subroutine watertable (state,node,nodlev,nodhlp,nodheq1,CritUndSatVol,flsat,waterlevel)
      ! [SS-SWC S-2.12B] h/Theta/ThetaS retired — read via soil (associate below)
      ! [GR-BH C4] mesh%numnod/mesh%dz/mesh%z migrated to mesh; use variables no longer needed here
      implicit none

      type(swap_state_t), intent(in) :: state
      integer node,nodhlp,nodheq1,nodlev
      logical flsat
      real(8) waterlevel
      integer i
      real(8) CritUndSatVol, TotUndSatVol
      logical flsat2

      ! S-2.4 ASSOCIATE: alias soil arrays for h/Theta/ThetaS reads
      ! [GR-BH C4] mesh globals aliased via mesh
      associate (mesh => state%mesh, soil => state%soilwater)

      TotUndSatVol = 0.0d0
      flsat2 = .false.
      i = node
      do while (TotUndSatVol.lt.CritUndSatVol .and. .not.flsat2 .and. i.ge.1)
         TotUndSatVol = TotUndSatVol + (soil%thetas(i) - soil%theta(i)) * mesh%dz(i)  ! S-2.4 read cutover
         if (soil%h(i).gt.-1.d-7) flsat2 = .true.                              ! S-2.4 read cutover
         i = i - 1
      enddo
!
      if (i.eq.0 .or. TotUndSatVol.gt.CritUndSatVol-1.d-8) then
         flsat = .false.
!!!!!!!! NODGWL is NOT node with GWL, but DEEPEST UNSATURATED NODE !!!!!!!!!!!!!
         nodlev = node
         nodhlp = i
      elseif(flsat2) then
         node   = i + 1
      endif
!
      if (.not.flsat) then
! find groundwater level containing node
         if (CritUndSatVol.gt.0.d0) then
            waterlevel = level(state,1,node,nodheq1)          ! S-2.4: pass state to level
         else
            waterlevel = level(state,2,node,nodheq1)          ! S-2.4: pass state to level
         endif
         i = max(node-2,1)
         do while(mesh%z(i)-0.5d0*mesh%dz(i).gt.waterlevel .and. i.gt.2 .and. i.lt.mesh%numnod)
            i = i + 1
         enddo
         nodlev = min(max(i,1),mesh%numnod)
      endif

      end associate  ! soil%h, soil%theta, soil%thetas (S-2.4); mesh%numnod/mesh%dz/mesh%z [GR-BH C4]

      return
    end subroutine watertable

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
      ! SS-SWST Phase 2 Task 11 A2: state added to fluxes() so qdra/qdrtot read from state.
      ! SS-SWC S-2.4: theta/thetm1/FrArMtrx/volact/volm1/fllowgwl/q/inq reads cut over to soil
      subroutine fluxes (state)
      ! [SS-SWC S-2.12B] q, inq, thetm1, theta, volact, volm1, FrArMtrx, fllowgwl retired from variables
      ! [GR-BH C4] mesh%numnod/mesh%dz migrated to mesh
      ! [GR-BH Audit 31] soil%swbotb_runtime/drai%nrlevs removed from use-list; read via state%
      ! [SS-GR-FINAL B9] DEFERRED: qimmob/QExcMpMtx/QMaPo — soil immobile water / macropore fluxes; Phase C3/D
      !   qssdi/qssdisum — SSDI source term arrays; needs irrigation_state_t; Phase C3
      use variables, only: qimmob, QExcMpMtx, QMaPo, qssdi, qssdisum
      ! [SS-TC TC-6] dt read cut over to time%dt
      use swap_state_mod, only: swap_state_t
      implicit none

      type(swap_state_t), intent(inout) :: state
      integer i,level

      ! S-2.4 ASSOCIATE: alias soil arrays/scalars for read/write cutover
      ! [SS-TC TC-6] time%dt aliases time%dt
      ! [GR-BH C4] mesh globals aliased via mesh
      ! [GR-BH Audit 31] soil%swbotb_runtime/drai%nrlevs aliased via state%
      associate (mesh => state%mesh,         &
                 soil => state%soilwater,    &
                 drai => state%drainage,     &
                 surf => state%surfacewater, &
                 time => state%timecontrol)

      ! determine qbot if not specified
      ! SS-BND B-2.7: qtop and qbot read/written via soil (globals retired)
      if (soil%swbotb_runtime .eq. 5 .or. soil%swbotb_runtime .eq. 7 .or.                         &
     &    soil%swbotb_runtime .eq. 8 .or. soil%swbotb_runtime .eq. -2 .or.                        &
     &    (soil%swbotb_runtime .eq. 1 .and. soil%fllowgwl)) then        ! S-2.4 read cutover: fllowgwl -> soil%fllowgwl
        ! SS-CRP Phase 2 Task C-2.2: qrosum -> soil%qrosum
        soil%qbot = soil%qtop + soil%qrosum + surf%qdrtot - QMaPo + (soil%volact-soil%volm1)/time%dt - qssdisum  ! S-2.4, TC-6
      endif

      ! calculate fluxes (cm/d) from changes in volume per compartment
      ! [SS-SWC S-2.12B] all legacy half-writes dropped — soil is canonical
      i = mesh%numnod+1
      soil%q(i)              = soil%qbot               ! S-1.6/S-2.12B
      soil%inq(i)            = soil%inq(i) + soil%q(i)*time%dt             ! S-2.12B, TC-6
      do i = mesh%numnod,1,-1
        soil%q(i) = - (soil%theta(i)-soil%thetm1(i)+qimmob(i))*soil%FrArMtrx(i)*mesh%dz(i)/time%dt +  &
     &                soil%q(i+1)-soil%qrot(i)+QExcMpMtx(i)+qssdi(i)        ! S-1.6/S-2.12B, TC-6

        if (allocated(drai%qdra)) then
          do level=1,drai%nrlevs
             soil%q(i) = soil%q(i) - drai%qdra(level,i)
          enddo
        end if
        soil%inq(i) = soil%inq(i) + soil%q(i)*time%dt                       ! S-1.6/S-2.12B, TC-6
      end do

      end associate  ! soil%theta etc (S-2.4/S-2.12B/TC-6); mesh%numnod/mesh%dz [GR-BH C4]; soil%swbotb_runtime/drai%nrlevs [GR-BH Audit 31]

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
      Use Variables
      ! [SS-TC TC-6] dt and fldaystart reads cut over to time via tc_* aliases below
      use iso_fortran_env, only: real64        ! [SS-SWC S-2.12B] needed for 0.0_real64 literal below
      implicit none

      type(swap_state_t), intent(inout) :: state
      integer node,level
      real(8) qrotts,qdrats,ptrats,pevats,revats,qbotts

      ! TC-6 ASSOCIATE: alias TC fields for dense dt / flDayStart reads in integral body
      ! [GR-BH C4] mesh%numnod aliased via mesh
      ! [GR-BH Audit 31] drai%nrlevs aliased via drai
      associate (mesh => state%mesh,         &
                 soil => state%soilwater,    &
                 drai => state%drainage,     &
                 surf => state%surfacewater, &
                 atmo => state%atmosphere,   &
                 crop => state%crop,         &
                 time => state%timecontrol)

      if (time%flZeroIntr) then
        ! SS-ATM Phase 2 Task A-2.2 (D6): igrai/inrai removed — canonical reset
        ! is atmo%intr%reset() invoked in meteoday's ResetMetFlx (A-2.1).
        ! [SS-SWC S-2.12B] iprec/igird/inird retired — soil%reset_intermediate() handles
      endif

      ! potential transpiration of this timestep
      ! SS-ATM Phase 2 Task A-2.2: ptra read from atmo (atmosphere home).
      ptrats = atmo%ptra * time%dt                         ! TC-6

      ! potential soil evaporation of this timestep
      ! SS-ATM Phase 2 Task A-2.2: peva read from atmo (atmosphere home).
      pevats = atmo%peva * time%dt                         ! TC-6

      ! reduced soil evaporation of this timestep
      ! SS-BND Phase 2 Task B-2.2: reva read from soil (boundary home).
      revats = soil%reva * time%dt                          ! TC-6

      ! flux lower boundary of this timestep
      ! SS-BND Phase 2 Task B-2.2: qbot read from soil (boundary home).
      qbotts = soil%qbot*time%dt                            ! TC-6

      ! total root extraction of this timestep
      ! SS-CRP Phase 2 Task C-2.2: qrosum -> soil%qrosum
      qrotts = soil%qrosum * time%dt                        ! TC-6

      ! total drainage flux of this timestep
      ! SS-SWST Phase 2 Task 11 A2: qdrtot removed from globals; read from state.
      qdrats = surf%qdrtot * time%dt                     ! TC-6

      ! determine daily actual transpiration — [SS-SWC S-2.12B] legacy half-writes dropped
      if (time%flDayStart) soil%tra = 0.0_real64              ! S-2.12B, TC-6: fldaystart -> time%flDayStart
      soil%tra = soil%tra + qrotts          ! S-1.5/S-2.12B

      ! add time step fluxes to intermediate totals
      soil%iqrot = soil%iqrot + qrotts      ! S-1.5/S-2.12B
      do node = 1,noddrz
        soil%inqrot(node) = soil%inqrot(node) + soil%qrot(node) * time%dt  ! S-2.12B, TC-6
        soil%qpotrot_day(node) = soil%qpotrot_day(node) + soil%qpotrot(node) * time%dt  ! S-2.12B, TC-6
        soil%qredtot_day(node) = soil%qredtot_day(node) + (soil%qredwet(node) + soil%qreddry(node) + soil%qredsol(node) + soil%qredfrs(node)) * time%dt  ! S-2.12B, TC-6
      end do
      do node = 1,mesh%numnod
        soil%inqssdi(node) = soil%inqssdi(node) + qssdi(node) * time%dt    ! S-2.12B, TC-6
        soil%iqssdi = soil%iqssdi + qssdi(node) * time%dt                  ! S-2.12B, TC-6
      end do
      soil%iqredwet = soil%iqredwet + soil%qredwetsum*time%dt   ! S-2.12B, TC-6
      soil%iqreddry = soil%iqreddry + soil%qreddrysum*time%dt   ! S-2.12B, TC-6
      soil%iqredsol = soil%iqredsol + soil%qredsolsum*time%dt   ! S-2.12B, TC-6
      soil%iqredfrs = soil%iqredfrs + soil%qredfrssum*time%dt   ! S-2.12B, TC-6
      soil%iqredwet_day = soil%iqredwet_day + soil%qredwetsum*time%dt  ! S-2.12B, TC-6
      soil%iqreddry_day = soil%iqreddry_day + soil%qreddrysum*time%dt  ! S-2.12B, TC-6
      soil%iqredsol_day = soil%iqredsol_day + soil%qredsolsum*time%dt  ! S-2.12B, TC-6
      soil%iqredfrs_day = soil%iqredfrs_day + soil%qredfrssum*time%dt  ! S-2.12B, TC-6
      soil%iptra_day = soil%iptra_day + atmo%ptra * time%dt           ! S-2.12B, TC-6
      soil%ies0 = soil%ies0 + 0.1d0*crop%es0*time%dt                                   ! S-2.12B, TC-6
      soil%iet0 = soil%iet0 + 0.1d0*crop%et0*time%dt                                   ! S-2.12B, TC-6
      soil%iew0 = soil%iew0 + 0.1d0*crop%ew0*time%dt                                   ! S-2.12B, TC-6

      ! SS-SWST Phase 2 Task 7: iqdra/inqdra* accumulated directly into state; global dropped.
      ! ADR 0031 Phase 2 Task 5: qdra global deleted; read from drai%qdra.
      surf%iqdra = surf%iqdra + qdrats + drai%QRapDra*time%dt  ! TC-6
      do node = 1,mesh%numnod
        qdraincomp(node) = 0.d0
        do level = 1,drai%nrlevs
          if (allocated(surf%inqdra) .and. allocated(drai%qdra)) then
            surf%inqdra(level,node) = surf%inqdra(level,node) + drai%qdra(level,node)*time%dt  ! TC-6
            if (drai%qdra(level,node) > 0.0d0) then
               surf%inqdra_out(level,node) = surf%inqdra_out(level,node) + drai%qdra(level,node)*time%dt  ! TC-6
            else
               surf%inqdra_in(level,node)  = surf%inqdra_in(level,node) - drai%qdra(level,node)*time%dt  ! TC-6
            end if
          end if
          if (allocated(drai%qdra)) then
            qdraincomp(node) = drai%qdra(level,node) + qdraincomp(node)
          end if
        end do
      end do

      ! [SS-SWC S-2.12B] all legacy half-writes dropped — soil is canonical
      ! SS-ATM Phase 2 Task A-2.2: aintcdt read from atmo (atmosphere home).
      soil%iintc = soil%iintc + (atmo%aintcdt+crop%gird-atmo%nird)*time%dt  ! S-2.12B, TC-6

      atmo%intr%iptra = atmo%intr%iptra + ptrats
      atmo%intr%ipeva = atmo%intr%ipeva + pevats
      atmo%intr%ievap = atmo%intr%ievap + revats
      ! SS-BND Phase 2 Task B-2.2: runots read from soil (boundary home).
      soil%iruno = soil%iruno + soil%runots                  ! S-2.12B
      soil%irunon = soil%irunon + soil%runon*time%dt              ! S-2.12B, TC-6
      ! SS-ATM Phase 2 Task A-2.2: graidt/nraidt read from atmo (atmosphere home).
      soil%iprec = soil%iprec + (atmo%graidt+crop%gird)*time%dt       ! S-2.12B, TC-6
      atmo%intr%igrai = atmo%intr%igrai + atmo%graidt*time%dt             ! TC-6
      soil%igird = soil%igird + crop%gird*time%dt                                  ! S-2.12B, TC-6
      atmo%intr%inrai = atmo%intr%inrai + atmo%nraidt*time%dt             ! TC-6
      soil%inird = soil%inird + atmo%nird*time%dt                                  ! S-2.12B, TC-6
      soil%iqbot = soil%iqbot + qbotts                                   ! S-2.12B
      if (soil%q(1) < 0.0d0) then
         soil%iqtdo = soil%iqtdo - soil%q(1)*time%dt              ! S-2.12B, TC-6
      else
         soil%iqtup = soil%iqtup + soil%q(1)*time%dt              ! S-2.12B, TC-6
      end if
      do node = 1, mesh%numnod+1
         if (soil%q(node) < 0.0d0) then
            soil%iqdo(node) = soil%iqdo(node) - soil%q(node)*time%dt  ! S-2.12B, TC-6
         else
            soil%iqup(node) = soil%iqup(node) + soil%q(node)*time%dt  ! S-2.12B, TC-6
         end if
      end do

      ! add time step fluxes to total cumulative values
      soil%cqssdi = soil%cqssdi + qssdisum*time%dt                            ! S-2.12B, TC-6
      soil%cqrot  = soil%cqrot  + qrotts                                 ! S-2.12B
      ! SS-SWST Phase 2 Task 7: cqdra accumulated directly into state; global dropped.
      surf%cqdra = surf%cqdra + qdrats
      atmo%cumu%cptra = atmo%cumu%cptra + ptrats
      atmo%cumu%cpeva = atmo%cumu%cpeva + pevats
      atmo%cumu%cevap = atmo%cumu%cevap + revats
      if (soil%runots.lt.0.0d0) then
        soil%cinund = soil%cinund - soil%runots              ! S-2.12B
      else if (soil%runots.gt.0.0d0) then
        soil%crunoff = soil%crunoff + soil%runots            ! S-2.12B
      endif
      soil%irunoCN   = soil%irunoCN   + atmo%Runoff_CN*time%dt
      soil%crunoffCN = soil%crunoffCN + atmo%Runoff_CN*time%dt

      ! SS-ATM Phase 2 Task A-2.2: aintcdt/graidt/nraidt read from atmo (atmosphere home).
      atmo%cumu%caintc = atmo%cumu%caintc + (atmo%aintcdt+crop%gird-atmo%nird)*time%dt  ! TC-6

      atmo%cumu%cgrai = atmo%cumu%cgrai + atmo%graidt*time%dt             ! TC-6
      atmo%cumu%cnrai = atmo%cumu%cnrai + atmo%nraidt*time%dt             ! TC-6
!      cnrai = cgrai - caintc
      soil%cgird = soil%cgird + crop%gird*time%dt                                  ! S-2.12B, TC-6
      soil%cnird = soil%cnird + atmo%nird*time%dt                                  ! S-2.12B, TC-6

      if (qbotts.lt.0.0d0) then
        soil%cqbotdo = soil%cqbotdo - qbotts                            ! S-2.12B
      else if (qbotts.gt.0.0d0) then
        soil%cqbotup = soil%cqbotup + qbotts                            ! S-2.12B
      endif
      soil%cqbot = soil%cqbot + qbotts                                   ! S-2.12B
      ! SS-SWST Phase 2 Task 7: cqdrain/in/out accumulated directly into state; globals dropped.
      if (allocated(surf%cqdrain)) then
        do level = 1,drai%nrlevs
          ! infiltration
          if (drai%qdrain(level).lt.0.0d0) then
            surf%cqdrainin(level) = surf%cqdrainin(level) - drai%qdrain(level)*time%dt  ! TC-6
          ! drainage
          else if (drai%qdrain(level).gt.0.0d0) then
            surf%cqdrainout(level) = surf%cqdrainout(level) + drai%qdrain(level)*time%dt  ! TC-6
          endif
          surf%cqdrain(level) = surf%cqdrain(level) + drai%qdrain(level)*time%dt  ! TC-6
        enddo
      end if

      ! rain on the ponding surface — [SS-SWC S-2.12B] legacy half-writes dropped
      soil%cqprai = soil%cqprai + atmo%nraidt*time%dt  ! S-2.12B, TC-6
      soil%crunon = soil%crunon + soil%runon*time%dt    ! S-2.12B, TC-6
      if (soil%q(1).lt.0.0d0) then
        soil%cqtdo = soil%cqtdo - soil%q(1)*time%dt     ! S-2.12B, TC-6
      else if (soil%q(1).gt.0.0d0) then
        soil%cqtup = soil%cqtup + soil%q(1)*time%dt     ! S-2.12B, TC-6
      endif

      ! compensate water balance error of this time step during remaining day part
      ! cumulative water balance error
      ! SS-SWC S-2.4: cnird/crunon/crunoff/cqrot/cqbot/volini/volact/PondIni/pond/cqssdi/cqprai
      !               read cutover to soil flat fields
      ! [SS-SWC S-2.12B] write directly to soil%wbalance — legacy global retired
      if (swsnow.eq.0) then
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

      ! [SS-GR-CROPRT A2] if (FlMacropore) wbalance block dropped (ADR 0040; always .false.)

      end associate  ! time%dt, time%flDayStart (TC-6); mesh%numnod [GR-BH C4]; drai%nrlevs [GR-BH Audit 31]

      return
      end

      !> @brief Check mass balance per output period
      !>
      !> @details
      !> Performs comprehensive mass balance checking for different subsystems:
      !> 1. Ponding layer (surface water and snow)
      !> 2. Total soil profile
      !> 3. Individual soil compartments
      !> 4. Macropore domains (Dm1 and Dm2, if macropores are enabled)
      !>
      !> Compares all water balance terms and writes deviations exceeding CritDevMasBal
      !> threshold to output file. Used for validation and debugging of water balance calculations.
      !>
      !> @param[inout] flopenfiledev Flag indicating if deviation file is open
      !> @param[in] inqdranew Drainage fluxes per level and compartment
      !> @param[in] iqexcmtxdm1cpnew Exchange from matrix to macropore domain 1 per compartment
      !> @param[in] iqexcmtxdm2cpnew Exchange from matrix to macropore domain 2 per compartment
      !> @param[in] inqnew Fluxes between compartments
      !> @param[in] iqoutdrrapcpnew Outflow from domains to drains/rapid drainage per compartment
      !> @param[in] inqrotnew Root extraction per compartment
      !> @param[in] ithetabegnew Initial water content per compartment
      !> @param[in] thetanew Current water content per compartment
      !>
      !> @note
      !> Date: 26-jun-2003
      !>
      !> Purpose: Checking of mass balance per period OutPer for ANIMO/PEARL output
      !>
      !> File usage: outfil
      !>
      !> SAVE removed - dev_cmb now in variables.f90 module
      !> @endnote
      subroutine checkmassbal (flopenfiledev,inqdranew,iqexcmtxdm1cpnew,iqexcmtxdm2cpnew,inqnew,iqoutdrrapcpnew,inqrotnew,ithetabegnew,thetanew,state)
      ! [SS-SWC S-2.12B] igird/inird/IPondBeg/iruno/irunon/pond retired — read from soil
      ! [GR-BH Audit 31] drai%nrlevs removed from use-list; read via drai%nrlevs
      ! [SS-GR-FINAL B9] DEFERRED: all checkmassbal globals; Phase C3/D
      !   NumNodNew/DZNew — regridding temporaries; Phase C3
      !   outfil/pathwork — output file path globals; Phase C3
      !   CritDevMasBal — mass balance criterion; Phase C3
      !   dev_cmb — mass balance device file handle; Phase C3
      ! [SS-GR-CROPRT A2] IcTopMp/FlMacropore/IQInTopVrt*/IQInTop*/IWaSrDm*/WaSrDm* dropped — retired (ADR 0040)
      use variables, only: NumNodNew, outfil, pathwork, DZNew,    &
                           CritDevMasBal, dev_cmb
      ! [SS-TC TC-6] DayCum read cut over to time%daycum
      use swap_state_mod, only: swap_state_t
      use swap_array_dimensions, only: macp, madr
      use file_io_mod, only: file_open
      implicit none

      ! -   global
      real(8) IQExcMtxDm1CpNew(macp), IQExcMtxDm2CpNew(macp)
      real(8) inqdraNew(Madr,macp)
      real(8) inqNew(macp+1), IQOutDrRapCpNew(macp), inqrotNew(macp)
      real(8) IThetaBegNew(MaCp),thetaNew(macp)
      logical FlOpenFileDev
      type(swap_state_t), intent(in) :: state  ! [SS-ATM A-2.6] for retired igrai/inrai/ievap/igSnow/isnrai/isubl

      ! local
      integer Level, ic
      real(8) DevMasBalDm1,DevMasBalDm2, DevMasBalCmp(MaCp)
      real(8) DevMasBalPnd, DevMasBalPrf, IQExcMtxDm1
      real(8) IQExcMtxDm2,IQInTopLatDm,  IQInTopPreDm, IQOutDrRap
      real(8) Qdra(MaCp), QdraPrf, QrotPrf, SrDif
      real(8) WaSr(MaCp), WaSrBeg(MaCp), WaSrPrf, WaSrPrfBeg  
      character(len=300) filnam
      logical FlWriteDevCmp(MaCp), FlWriteDev, FlWriteDevDm1 
      logical FlWriteDevDm2, FlWriteDevPnd, FlWriteDevPrf

      associate (soil => state%soilwater,    &
                 atmo => state%atmosphere,   &
                 drai => state%drainage,     &
                 time => state%timecontrol)

      ! Checking of mass balances of sub systems per period OutPer
      FlWriteDev= .false.
      FlWriteDevPnd = .false.
      FlWriteDevPrf = .false.
      do ic= 1, NumNodNew
         FlWriteDevCmp(ic) = .false.
      enddo
      FlWriteDevDm1 = .false.
      FlWriteDevDm2 = .false.

      ! 1) Ponding layer
      ! SS-ATM A-2.6: Ssnow retired — read from atmo%ssnow
      ! [SS-SWC S-2.12B] IPondBeg/pond -> soil
      SrDif = soil%IPondBeg-soil%pond + atmo%ISsnowBeg-atmo%ssnow
      IQInTopPreDm= 0.d0
      IQInTopLatDm= 0.d0
      ! [SS-GR-CROPRT A2] if (FlMacropore) IQInTopPreDm/IQInTopLatDm block dropped (ADR 0040)

      ! Deviation mass balance Ponding layer in cm
      ! [SS-SWC S-2.12B] igird/inird/iruno/irunon -> soil
      DevMasBalPnd = atmo%intr%igrai + atmo%intr%igsnow + soil%igird + soil%irunon + inqNew(1) + SrDif &
     &             - (atmo%intr%igrai-atmo%intr%inrai-atmo%intr%isnrai + soil%igird-soil%inird + atmo%intr%isubl + atmo%intr%ievap + soil%iruno) &
     &             - IQInTopPreDm - IQInTopLatDm

      ! Check mass balance against criteria
      if (abs(DevMasBalPnd).gt.CritDevMasBal) then
         FlWriteDev = .true.
         FlWriteDevPnd = .true.
      endif

      ! 2) Total Soil Profile
      WaSrPrfBeg = 0.d0 
      WaSrPrf    = 0.d0 
      QrotPrf    = 0.d0 
      QdraPrf    = 0.d0 
      IQExcMtxDm1= 0.d0
      IQExcMtxDm2= 0.d0
      do ic = 1, numnodnew
         WaSrPrfBeg= WaSrPrfBeg + dzNew(ic)*IThetaBegNew(ic)
         WaSrPrf   = WaSrPrf    + dzNew(ic)*ThetaNew(ic)
         QrotPrf   = QrotPrf    + inqrotNew(ic)
         do level=1,drai%nrlevs   ! [GR-BH Audit 31]
            QdraPrf = QdraPrf + InqdraNew(level,ic)
         enddo
         ! [SS-GR-CROPRT A2] if (FlMacropore) IQExcMtxDm1/2 accumulation dropped (ADR 0040)
      enddo
      SrDif= WaSrPrfBeg - WaSrPrf

      ! Deviation mass balance total Profile in cm
      DevMasBalPrf = inqNew(NumNodNew+1) + SrDif + IQExcMtxDm1 + &
     &               IQExcMtxDm2 - (inqNew(1) + QrotPrf + QdraPrf)

      ! Check mass balance against criteria
      if (abs(DevMasBalPrf).gt.CritDevMasBal) then
         FlWriteDev = .true.
         FlWriteDevPrf = .true.
      endif

      ! 3) Individual Soil Compartments
      do 100 ic = 1, numnodnew
         SrDif = 0.d0
         Qdra(ic)= 0.d0 
         WaSrBeg(ic)= dzNew(ic) * IThetaBegNew(ic)
         WaSr(ic)   = dzNew(ic) * ThetaNew(ic)
         SrDif = WaSrBeg(ic) - WaSr(ic) 
         do level=1,drai%nrlevs   ! [GR-BH Audit 31]
            Qdra(ic) = Qdra(ic) + inqdraNew(level,ic)
         enddo

         ! Deviation mass balance Soil Compartments in cm
         DevMasBalCmp(ic) = inqNew(ic+1) + SrDif &
     &                    - (inqNew(ic) + inqrotNew(ic) + Qdra(ic))
         ! [SS-GR-CROPRT A2] if (FlMacropore) DevMasBalCmp IQExc adjustment dropped (ADR 0040)

         ! Check mass balance against criteria
         if (abs(DevMasBalCmp(ic)).gt.CritDevMasBal) then
            FlWriteDev = .true.
            FlWriteDevCmp(ic) = .true.
         endif
 100  continue

      ! [SS-GR-CROPRT A2] Section 4 Macropore domains Dm1/Dm2 dropped (ADR 0040; FlMacropore always .false.)

      ! In case of deviations of mass balance open file 'xxxxx.dwb.csv'

      if (FlWriteDev .and. .not.FlOpenFileDev) then
         filnam = trim(pathwork)//trim(outfil)//'.dwb'
         call file_open(dev_cmb, filnam, 'replace', 'write')
         write(dev_cmb,1)
         ! [SS-GR-CROPRT A2] if (FlMacropore) write(dev_cmb,2) dropped (ADR 0040)
         FlOpenFileDev = .true.
      endif

      ! Write deviations of water balance Top system
      ! SS-ATM A-2.6: igrai/inrai/igsnow/isnrai/isubl/ievap retired — read from atmo%intr
      ! [SS-SWC S-2.12B] igird/irunon/inird/iruno/Pond/IPondBeg -> soil
      if (FlWriteDevPnd) write(dev_cmb,3) time%daycum, DevMasBalPnd, &  ! TC-6: daycum -> time
     &    atmo%intr%igrai, atmo%intr%igsnow, soil%igird, soil%irunon, atmo%intr%isnrai, &
     &    atmo%intr%igrai-atmo%intr%inrai, soil%igird-soil%inird, &
     &    atmo%intr%isubl, atmo%intr%ievap, soil%iruno, inqNew(1), soil%pond, soil%IPondBeg, atmo%ssnow, &
     &    atmo%ISsnowBeg,IQInTopPreDm, IQInTopLatDm

      ! Write deviations of water balance whole Profile
      if (FlWriteDevPrf) write(dev_cmb,4) time%daycum, DevMasBalPrf, &  ! TC-6
     &    inqNew(1), inqNew(NumNodNew+1), QrotPrf, QdraPrf, WaSrPrf, &
     &    WaSrPrfBeg, IQExcMtxDm1, IQExcMtxDm2

      ! Write deviations of water balance of Individual Soil Compartments
      do ic= 1, numnodnew
         if (FlWriteDevCmp(ic)) write(dev_cmb,5) time%daycum,ic,DevMasBalCmp(ic), &  ! TC-6
     &      inqNew(ic), inqNew(ic+1), inqrotNew(ic), Qdra(ic), WaSr(ic), &
     &      WaSrBeg(ic), IQExcMtxDm1CpNew(ic), IQExcMtxDm2CpNew(ic)
      enddo

      ! [SS-GR-CROPRT A2] Macropore Domain Dm1/Dm2 write statements dropped (ADR 0040; FlMacropore always .false.)
    1 format(' DEVIATIONS WATERBALANCE for different subsystems: 1. Pon'&
     &'d.layer; 2. Whole profile; 3. Compartment; (optional: Macrop.Dom'&
     &'.: 4. Dom1; 5. Dom2)',/,                                         &
     &' Relevant terms of waterbalance per subsystem',                  &
     &' (all terms in cm):',//,                                         &
     &' DayCum, 1. PONDLAY., DevMasBalAbs, IgRai, IgSnow, IgIrd, IRunon'&
     &', SnowFall,IntcpRai, IntcpIrd, ISubl, IEvap, IRuno, InQTop,   ', &
     &'Pond, IPondBeg, Ssnow, ISsnowBeg, IQInTopPreDm, IQInTopLatDm,',/,&
     &' , 2. PROFILE, DevMasBalPrf, InQTop, InQBot, QrotPrf, QdraPrf,', &
     &' WaSrPrf, WaSrPrfBeg, InQExcMtxDm1, InQExcMtxDm2',/,             &
     &' , 3. COMPno, DevMasBalCmp, InQNew(top), InQNew(bot),',          &
     &' InQrotNew, Qdra, WaSr, WaSrBeg, InQExcMtxDm1CpNew,',            &
     &' InQExcMtxDm2CpNew')
    2 format(' , 4. MPDOM1, DevMasBalDm1, IQInTopPre/VrtDm1,',          &
     &' IQInTopLatDm1, InQExcMtxDm1, WaSrDm1, IWaSrDm1Beg, InQOutDrRap'/&
     &' , 5. MPDOM2, DevMasBalDm2, IQInTopPre/VrtDm2, IQInTopLatDm2,    &
     & InQExcMtxDm2, WaSrDm2, IWaSrDm2Beg')
    3 format(i5,',',' Pondlay. : ',18(',',f12.8))
    4 format(i5,',',' Profile : ',9(',',f12.8))
    5 format(i5,',',' Comp',i3,': ',9(',',f12.8))
    6 format(i5,',',' MpDom1 : ',7(',',f12.8))
    7 format(i5,',',' MpDom2 : ',6(',',f12.8))

      end associate

      return
      end

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
      ! SS-SWC S-2.4: theta/FrArMtrx/volact reads cut over to soil
      subroutine watstor (state)
      ! [SS-SWC S-2.12B] volm1/volact/theta/FrArMtrx retired from variables; all via soil
      ! [GR-BH C4] mesh%numnod/mesh%dz migrated to mesh; use variables no longer needed here
      use swap_state_mod, only: swap_state_t
      IMPLICIT NONE

      type(swap_state_t), intent(inout) :: state
      INTEGER i

      ! [GR-BH C4] mesh globals aliased via mesh
      associate (mesh => state%mesh, soil => state%soilwater)

      ! update soil profile water storage — [SS-SWC S-2.12B] legacy half-writes dropped
      soil%volm1  = soil%volact
      soil%volact = 0.0d0
      do 10 i = 1,mesh%numnod
        soil%volact = soil%volact + soil%theta(i)*mesh%dz(i)*soil%FrArMtrx(i)
 10   continue

      end associate  ! mesh%numnod/mesh%dz [GR-BH C4]

      return
      end

end module