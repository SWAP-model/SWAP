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
      ! [SS-SWC S-2.12B] retired globals removed from use clause; all reads/writes via state%soilwater
      ! [GR-BH C4] numnod/z/disnod migrated to state%mesh
      use variables, only: logf,swbotb,flmacropore,CritUndSatVol
      ! [SS-TC TC-6] t1900 read cut over to state%timecontrol%t1900
      use swap_log, only: log_debug, to_str
      implicit none
      ! SS-BND B-2.7: state added to read state%soilwater%gwlinp (gwlinp global retired)
      ! SS-SWC S-1.6: intent(in) -> intent(inout) to allow gwl/nodgwl/pegwl/bpegwl/npegwl/gwlflcpzo/nodgwlflcpzo dual-writes
      ! SS-SWC S-2.4: h/pond reads cut over to state%soilwater; level/watertable pass state
      type(swap_state_t), intent(inout) :: state
      ! local
      integer   i, node, nodhlp, nodheq1
      integer   nodgwlflcpzo_loc
      real(8)   gwlflcpzo_loc
      logical   flsat,flunsat
      character(len=200) messag
      character(len=19) datexti

      ! S-2.4 ASSOCIATE: alias state%soilwater arrays for h/pond reads
      ! [GR-BH C4] mesh globals aliased via state%mesh
      associate( sw_h   => state%soilwater%h,   &
                 sw_pond => state%soilwater%pond, &
                 sw_gwl  => state%soilwater%gwl,  &
                 numnod  => state%mesh%numnod,     &  ! [GR-BH C4]
                 z       => state%mesh%z,          &  ! [GR-BH C4]
                 disnod  => state%mesh%disnod       )  ! [GR-BH C4]

      ! set initial values — [SS-SWC S-2.12B] legacy half-writes dropped
      sw_gwl    = 999.0d0                                     ! S-1.6/S-2.12B
      state%soilwater%pegwl  = 999.0d0                        ! S-1.6/S-2.12B
      flsat     = .false.
      state%soilwater%nodgwl = numnod+1                       ! S-1.6/S-2.12B
      nodhlp    = numnod
      nodheq1   = numnod

      ! search for groundwater table
      if (sw_h(numnod).ge.0.0d0) flsat  = .true.      ! S-2.4 read cutover: h -> sw_h

      node = numnod
      nodgwlflcpzo_loc = numnod + 1
      state%soilwater%nodgwlflcpzo = numnod + 1               ! S-1.6/S-2.12B
      gwlflcpzo_loc = sw_gwl
      state%soilwater%gwlflcpzo    = sw_gwl                   ! S-1.6/S-2.12B
      do while (flsat .and. node.gt.1)
         node = node - 1
         if(swbotb.eq.1)then
            if (sw_h(node) .lt. 0.0d0) then                  ! S-2.4 read cutover
               sw_gwl = z(node+1) + sw_h(node+1) / (sw_h(node+1)-sw_h(node)) * disnod(node+1)  ! S-2.4/S-2.12B
               flsat   =.false.
               state%soilwater%nodgwl = node                  ! S-1.6/S-2.12B
            endif
         else
            if (sw_h(node) .lt. 1.0d0 .and. nodheq1.eq.numnod) nodheq1 = node  ! S-2.4

            if (sw_h(node) .lt. 0.0d0) then                  ! S-2.4 read cutover
               if (.not.flmacropore) then
                  flsat  = .false.
                  state%soilwater%nodgwl = node                ! S-1.6/S-2.12B
                  sw_gwl    = level (state,1,node,nodheq1)     ! S-2.4/S-2.12B
               elseif (flmacropore) then
                  if (sw_gwl.gt.990.0d0) then
                     state%soilwater%nodgwl = node             ! S-1.6/S-2.12B
                     sw_gwl    = level (state,2,node,nodheq1)  ! S-2.4/S-2.12B
                  endif
                  call watertable (state,node,nodgwlflcpzo_loc,nodhlp,nodheq1,0.0d0,flsat,gwlflcpzo_loc)  ! S-2.4
                  state%soilwater%nodgwlflcpzo = nodgwlflcpzo_loc  ! S-1.6/S-2.12B (watertable out-arg)
                  state%soilwater%gwlflcpzo    = gwlflcpzo_loc     ! S-1.6/S-2.12B (watertable out-arg)
               endif
            endif
         endif
      end do

      ! whole profile saturated, then add ponding layer to groundwater level
      if (flsat)then
         if(sw_h(1) .gt. 0.0d0)then                          ! S-2.4 read cutover
            if (sw_pond .lt. 1.d-8) then                     ! S-2.4 read cutover: pond -> sw_pond
               sw_gwl = min(z(1)+sw_h(1),sw_pond)            ! S-2.4/S-2.12B
            else
               sw_gwl = sw_pond                              ! S-2.4/S-2.12B
            endif
         else
            sw_gwl = 0.0d0
         end if
         state%soilwater%nodgwl = 1                           ! S-1.6/S-2.12B
         if (flmacropore) then
            state%soilwater%nodgwlflcpzo = 1                  ! S-1.6/S-2.12B
            state%soilwater%gwlflcpzo    = sw_gwl             ! S-1.6/S-2.12B
         endif
      endif

      ! search for perched groundwater table

      ! first, search for first saturated compartment (i) above groundwater level
      i = nodhlp
      flunsat = .true.
      do while (flunsat .and. i.ge.1)
         if (sw_h(i).ge.0.0d0) flunsat = .false.             ! S-2.4 read cutover
         i = i - 1
      enddo

      ! if saturated compartment above gwl exists, then find perched groundwater table
      if (i.ne.0) then
         flsat  = .true.
         state%soilwater%bpegwl = i                           ! S-1.6/S-2.12B
         node   = state%soilwater%bpegwl
         nodheq1 = state%soilwater%bpegwl

         do while (flsat .and. node.gt.1)
            node = node - 1

            if (sw_h(node) .lt. 1.0d0 .and. nodheq1.eq.state%soilwater%bpegwl) nodheq1 = node  ! S-2.4

            if (sw_h(node) .lt. 0.0d0) then                  ! S-2.4 read cutover
               if (.not.flmacropore) then
                  flsat = .false.
                  state%soilwater%npegwl = node                ! S-1.6/S-2.12B
                  state%soilwater%pegwl  = level (state,1,node,nodheq1)  ! S-2.4/S-2.12B
               elseif (flmacropore) then
                  ! [SS-SWC S-2.12B] use local out args for watertable interface
                  call watertable (state,node,state%soilwater%npegwl,nodhlp,nodheq1,CritUndSatVol,flsat,state%soilwater%pegwl)  ! S-2.4
               endif
            endif
         end do

         ! whole profile saturated, then add ponding layer to perched groundwater level
         if (flsat)then
            if(sw_h(1) .gt. 0.0d0)then                       ! S-2.4 read cutover
               if (sw_pond .lt. 1.d-8) then                  ! S-2.4 read cutover
                  state%soilwater%pegwl = min(z(1)+sw_h(1),sw_pond)  ! S-2.4/S-2.12B
               else
                  state%soilwater%pegwl = sw_pond            ! S-2.4/S-2.12B
               endif
            else
               state%soilwater%pegwl = 0.0d0                  ! S-2.12B
            end if
            state%soilwater%npegwl = 1                        ! S-1.6/S-2.12B
         endif
      else
         state%soilwater%bpegwl = -1                          ! S-1.6/S-2.12B
         state%soilwater%npegwl = -1                          ! S-1.6/S-2.12B
      endif

      end associate  ! sw_h, sw_pond, sw_gwl (S-2.4/S-2.12B); numnod/z/disnod [GR-BH C4]

      ! fatal error if gwl below profile and flux has to be calculated
      if ((swbotb.eq.3.or.swbotb.eq.4).and.state%soilwater%gwl.gt.998.0d0) then
          messag = 'The groundwater level descends below the lower' &
     &     //' boundary. This conflicts with bottom boundary' &
     &     //' condition 3 and 4. Extend soil profile!'
         call fatalerr_collected ('calcgwl',messag)
      endif

      ! warning error if there is inconsistency between defined gwl and soil physics
      if (swbotb.eq.1 .and. (state%soilwater%gwlinp .ge.state%mesh%z(1) .or. state%soilwater%gwl.gt.998.0d0)) then  ! [GR-BH C4]
         ! determine date and date-time
         call dtdpst('year-month-day,hour:minute:seconds',state%timecontrol%t1900,datexti)  ! TC-6: t1900 -> state%timecontrol
         write(messag,'(6a)')                                           &
     &         'No groundwater level because unsaturation at bottom ',  &
     &         'compartment ( ', datexti,  ' ). ',                      &
     &         'This is caused by inconsistency between ',              &
     &         'given gwl and soil physical parameters '
         call warn ('Calcgwl',messag,logf,state%timecontrol%swscre)  ! [SS-BMI2 Task 4]
      endif

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
      ! SS-SWC S-2.4: state added as first arg so h reads come from state%soilwater%h
      function level (state,swoptlev,node,nodheq1)
      ! [SS-SWC S-2.12B] h retired from use clause; read via state%soilwater%h (associate below)
      ! [GR-BH C4] numnod/disnod/dz/z/zbotcp migrated to state%mesh; use variables no longer needed here
      implicit none

      type(swap_state_t), intent(in) :: state
      integer node, nodheq1, swoptlev
      integer i
      real(8) levm1, levp1
      real(8) level

      ! S-2.4 ASSOCIATE: alias state%soilwater%h for reads inside level
      ! [GR-BH C4] mesh globals aliased via state%mesh
      associate( sw_h   => state%soilwater%h,   &
                 numnod => state%mesh%numnod,     &  ! [GR-BH C4]
                 disnod => state%mesh%disnod,     &  ! [GR-BH C4]
                 dz     => state%mesh%dz,         &  ! [GR-BH C4]
                 z      => state%mesh%z,          &  ! [GR-BH C4]
                 zbotcp => state%mesh%zbotcp       )  ! [GR-BH C4]

      if (swoptlev.eq.1) then
         ! groundwater level equals elevation head where h = 0
         if (sw_h(node+1).ge.0.0d0)then                     ! S-2.4 read cutover
            level = z(node+1) + sw_h(node+1) / (sw_h(node+1)-sw_h(node)) * disnod(node+1)  ! S-2.4
         else
            level = zbotcp(node) - sw_h(node)               ! S-2.4
            level = min(z(node),max(zbotcp(node),level))
         end if

      elseif (swoptlev.eq.2) then
         ! groundwater level equals average of elevation heads of h = -1 and h = +1
         ! elevation head of h = +1
         i = nodheq1
         if (nodheq1.eq.numnod) then
            levp1 = z(i) - 0.5d0 * dz(i)
         else
            levp1 = z(i) - (z(i) - z(i+1)) * (1.d0-sw_h(i)) / (sw_h(i+1)-sw_h(i))  ! S-2.4
         endif
         ! elevation head of h = -1
         i = node
         do while (sw_h(i).gt.-1.d0 .and. i.gt.1)          ! S-2.4 read cutover
            i = i - 1
         enddo
         if (i.eq.1 .and. sw_h(1).gt.-1.d0 .and. node.gt.2) then   ! S-2.4
            ! no compartment with pressure head < -1 cm in top of profile:
            ! use elevation head of h = 0 as estimation for groundwater level
            levm1 = z(node+1) + sw_h(node+1) / (sw_h(node+1)-sw_h(node)) * disnod(node+1)  ! S-2.4
            levp1 = levm1
         else
            levm1 = z(i+1) + (z(i) - z(i+1)) * (1.d0+sw_h(i+1)) / (sw_h(i+1)-sw_h(i))  ! S-2.4
         endif
         ! groundwater level = average of levp1 and levm1
         level = (levp1 + levm1) / 2.d0
      endif

      end associate  ! sw_h (S-2.4); numnod/disnod/dz/z/zbotcp [GR-BH C4]

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
      ! SS-SWC S-2.4: state added as first arg so h/Theta/ThetaS reads come from state%soilwater
      subroutine watertable (state,node,nodlev,nodhlp,nodheq1,CritUndSatVol,flsat,waterlevel)
      ! [SS-SWC S-2.12B] h/Theta/ThetaS retired — read via state%soilwater (associate below)
      ! [GR-BH C4] numnod/dz/z migrated to state%mesh; use variables no longer needed here
      implicit none

      type(swap_state_t), intent(in) :: state
      integer node,nodhlp,nodheq1,nodlev
      logical flsat
      real(8) waterlevel
      integer i
      real(8) CritUndSatVol, TotUndSatVol
      logical flsat2

      ! S-2.4 ASSOCIATE: alias state%soilwater arrays for h/Theta/ThetaS reads
      ! [GR-BH C4] mesh globals aliased via state%mesh
      associate( sw_h      => state%soilwater%h,     &
                 sw_theta  => state%soilwater%theta,  &
                 sw_thetas => state%soilwater%thetas, &
                 numnod    => state%mesh%numnod,       &  ! [GR-BH C4]
                 dz        => state%mesh%dz,           &  ! [GR-BH C4]
                 z         => state%mesh%z              )  ! [GR-BH C4]

      TotUndSatVol = 0.0d0
      flsat2 = .false.
      i = node
      do while (TotUndSatVol.lt.CritUndSatVol .and. .not.flsat2 .and. i.ge.1)
         TotUndSatVol = TotUndSatVol + (sw_thetas(i) - sw_theta(i)) * dz(i)  ! S-2.4 read cutover
         if (sw_h(i).gt.-1.d-7) flsat2 = .true.                              ! S-2.4 read cutover
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
         do while(z(i)-0.5d0*dz(i).gt.waterlevel .and. i.gt.2 .and. i.lt.numnod)
            i = i + 1
         enddo
         nodlev = min(max(i,1),numnod)
      endif

      end associate  ! sw_h, sw_theta, sw_thetas (S-2.4); numnod/dz/z [GR-BH C4]

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
      ! SS-SWC S-2.4: theta/thetm1/FrArMtrx/volact/volm1/fllowgwl/q/inq reads cut over to state%soilwater
      subroutine fluxes (state)
      ! [SS-SWC S-2.12B] q, inq, thetm1, theta, volact, volm1, FrArMtrx, fllowgwl retired from variables
      ! [GR-BH C4] numnod/dz migrated to state%mesh
      use variables, only: qimmob,swbotb,QExcMpMtx,QMaPo,nrlevs,qssdi, qssdisum
      ! [SS-TC TC-6] dt read cut over to state%timecontrol%dt
      use swap_state_mod, only: swap_state_t
      implicit none

      type(swap_state_t), intent(inout) :: state
      integer i,level

      ! S-2.4 ASSOCIATE: alias state%soilwater arrays/scalars for read/write cutover
      ! [SS-TC TC-6] tc_dt aliases state%timecontrol%dt
      ! [GR-BH C4] mesh globals aliased via state%mesh
      associate( sw_theta    => state%soilwater%theta,   &
                 sw_thetm1   => state%soilwater%thetm1,  &
                 sw_FrArMtrx => state%soilwater%FrArMtrx, &
                 sw_q        => state%soilwater%q,       &
                 sw_inq      => state%soilwater%inq,     &
                 sw_volact   => state%soilwater%volact,  &
                 sw_volm1    => state%soilwater%volm1,   &
                 sw_fllowgwl => state%soilwater%fllowgwl, &
                 tc_dt       => state%timecontrol%dt,     &  ! TC-6
                 numnod      => state%mesh%numnod,         &  ! [GR-BH C4]
                 dz          => state%mesh%dz               )  ! [GR-BH C4]

      ! determine qbot if not specified
      ! SS-BND B-2.7: qtop and qbot read/written via state%soilwater (globals retired)
      if (swbotb .eq. 5 .or. swbotb .eq. 7 .or.                         &
     &    swbotb .eq. 8 .or. swbotb .eq. -2 .or.                        &
     &    (swbotb .eq. 1 .and. sw_fllowgwl)) then        ! S-2.4 read cutover: fllowgwl -> sw_fllowgwl
        ! SS-CRP Phase 2 Task C-2.2: qrosum -> state%soilwater%qrosum
        state%soilwater%qbot = state%soilwater%qtop + state%soilwater%qrosum + state%surfacewater%qdrtot - QMaPo + (sw_volact-sw_volm1)/tc_dt - qssdisum  ! S-2.4, TC-6
      endif

      ! calculate fluxes (cm/d) from changes in volume per compartment
      ! [SS-SWC S-2.12B] all legacy half-writes dropped — state%soilwater is canonical
      i = numnod+1
      sw_q(i)              = state%soilwater%qbot               ! S-1.6/S-2.12B
      sw_inq(i)            = sw_inq(i) + sw_q(i)*tc_dt             ! S-2.12B, TC-6
      do i = numnod,1,-1
        sw_q(i) = - (sw_theta(i)-sw_thetm1(i)+qimmob(i))*sw_FrArMtrx(i)*dz(i)/tc_dt +  &
     &                sw_q(i+1)-state%soilwater%qrot(i)+QExcMpMtx(i)+qssdi(i)        ! S-1.6/S-2.12B, TC-6

        if (allocated(state%drainage%qdra)) then
          do level=1,nrlevs
             sw_q(i) = sw_q(i) - state%drainage%qdra(level,i)
          enddo
        end if
        sw_inq(i) = sw_inq(i) + sw_q(i)*tc_dt                       ! S-1.6/S-2.12B, TC-6
      end do

      end associate  ! sw_theta etc (S-2.4/S-2.12B/TC-6); numnod/dz [GR-BH C4]

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
      ! [SS-TC TC-6] dt and fldaystart reads cut over to state%timecontrol via tc_* aliases below
      use iso_fortran_env, only: real64        ! [SS-SWC S-2.12B] needed for 0.0_real64 literal below
      implicit none

      type(swap_state_t), intent(inout) :: state
      integer node,level
      real(8) qrotts,qdrats,ptrats,pevats,revats,qbotts

      ! TC-6 ASSOCIATE: alias TC fields for dense dt / flDayStart reads in integral body
      ! [GR-BH C4] numnod aliased via state%mesh
      associate( tc_dt         => state%timecontrol%dt,         &  ! TC-6
                 tc_flDayStart => state%timecontrol%flDayStart,  &  ! TC-6
                 numnod        => state%mesh%numnod               )  ! [GR-BH C4]

      if (state%timecontrol%flZeroIntr) then
        ! SS-ATM Phase 2 Task A-2.2 (D6): igrai/inrai removed — canonical reset
        ! is state%atmosphere%intr%reset() invoked in meteoday's ResetMetFlx (A-2.1).
        ! [SS-SWC S-2.12B] iprec/igird/inird retired — state%soilwater%reset_intermediate() handles
      endif

      ! potential transpiration of this timestep
      ! SS-ATM Phase 2 Task A-2.2: ptra read from state%atmosphere (atmosphere home).
      ptrats = state%atmosphere%ptra * tc_dt                         ! TC-6

      ! potential soil evaporation of this timestep
      ! SS-ATM Phase 2 Task A-2.2: peva read from state%atmosphere (atmosphere home).
      pevats = state%atmosphere%peva * tc_dt                         ! TC-6

      ! reduced soil evaporation of this timestep
      ! SS-BND Phase 2 Task B-2.2: reva read from state%soilwater (boundary home).
      revats = state%soilwater%reva * tc_dt                          ! TC-6

      ! flux lower boundary of this timestep
      ! SS-BND Phase 2 Task B-2.2: qbot read from state%soilwater (boundary home).
      qbotts = state%soilwater%qbot*tc_dt                            ! TC-6

      ! total root extraction of this timestep
      ! SS-CRP Phase 2 Task C-2.2: qrosum -> state%soilwater%qrosum
      qrotts = state%soilwater%qrosum * tc_dt                        ! TC-6

      ! total drainage flux of this timestep
      ! SS-SWST Phase 2 Task 11 A2: qdrtot removed from globals; read from state.
      qdrats = state%surfacewater%qdrtot * tc_dt                     ! TC-6

      ! determine daily actual transpiration — [SS-SWC S-2.12B] legacy half-writes dropped
      if (tc_flDayStart) state%soilwater%tra = 0.0_real64              ! S-2.12B, TC-6: fldaystart -> tc_flDayStart
      state%soilwater%tra = state%soilwater%tra + qrotts          ! S-1.5/S-2.12B

      ! add time step fluxes to intermediate totals
      state%soilwater%iqrot = state%soilwater%iqrot + qrotts      ! S-1.5/S-2.12B
      do node = 1,noddrz
        state%soilwater%inqrot(node) = state%soilwater%inqrot(node) + state%soilwater%qrot(node) * tc_dt  ! S-2.12B, TC-6
        state%soilwater%qpotrot_day(node) = state%soilwater%qpotrot_day(node) + state%soilwater%qpotrot(node) * tc_dt  ! S-2.12B, TC-6
        state%soilwater%qredtot_day(node) = state%soilwater%qredtot_day(node) + (state%soilwater%qredwet(node) + state%soilwater%qreddry(node) + state%soilwater%qredsol(node) + state%soilwater%qredfrs(node)) * tc_dt  ! S-2.12B, TC-6
      end do
      do node = 1,numnod
        state%soilwater%inqssdi(node) = state%soilwater%inqssdi(node) + qssdi(node) * tc_dt    ! S-2.12B, TC-6
        state%soilwater%iqssdi = state%soilwater%iqssdi + qssdi(node) * tc_dt                  ! S-2.12B, TC-6
      end do
      state%soilwater%iqredwet = state%soilwater%iqredwet + state%soilwater%qredwetsum*tc_dt   ! S-2.12B, TC-6
      state%soilwater%iqreddry = state%soilwater%iqreddry + state%soilwater%qreddrysum*tc_dt   ! S-2.12B, TC-6
      state%soilwater%iqredsol = state%soilwater%iqredsol + state%soilwater%qredsolsum*tc_dt   ! S-2.12B, TC-6
      state%soilwater%iqredfrs = state%soilwater%iqredfrs + state%soilwater%qredfrssum*tc_dt   ! S-2.12B, TC-6
      state%soilwater%iqredwet_day = state%soilwater%iqredwet_day + state%soilwater%qredwetsum*tc_dt  ! S-2.12B, TC-6
      state%soilwater%iqreddry_day = state%soilwater%iqreddry_day + state%soilwater%qreddrysum*tc_dt  ! S-2.12B, TC-6
      state%soilwater%iqredsol_day = state%soilwater%iqredsol_day + state%soilwater%qredsolsum*tc_dt  ! S-2.12B, TC-6
      state%soilwater%iqredfrs_day = state%soilwater%iqredfrs_day + state%soilwater%qredfrssum*tc_dt  ! S-2.12B, TC-6
      state%soilwater%iptra_day = state%soilwater%iptra_day + state%atmosphere%ptra * tc_dt           ! S-2.12B, TC-6
      state%soilwater%ies0 = state%soilwater%ies0 + 0.1d0*es0*tc_dt                                   ! S-2.12B, TC-6
      state%soilwater%iet0 = state%soilwater%iet0 + 0.1d0*et0*tc_dt                                   ! S-2.12B, TC-6
      state%soilwater%iew0 = state%soilwater%iew0 + 0.1d0*ew0*tc_dt                                   ! S-2.12B, TC-6

      ! SS-SWST Phase 2 Task 7: iqdra/inqdra* accumulated directly into state; global dropped.
      ! ADR 0031 Phase 2 Task 5: qdra global deleted; read from state%drainage%qdra.
      state%surfacewater%iqdra = state%surfacewater%iqdra + qdrats + QRapDra*tc_dt  ! TC-6
      do node = 1,numnod
        qdraincomp(node) = 0.d0
        do level = 1,nrlevs
          if (allocated(state%surfacewater%inqdra) .and. allocated(state%drainage%qdra)) then
            state%surfacewater%inqdra(level,node) = state%surfacewater%inqdra(level,node) + state%drainage%qdra(level,node)*tc_dt  ! TC-6
            if (state%drainage%qdra(level,node) > 0.0d0) then
               state%surfacewater%inqdra_out(level,node) = state%surfacewater%inqdra_out(level,node) + state%drainage%qdra(level,node)*tc_dt  ! TC-6
            else
               state%surfacewater%inqdra_in(level,node)  = state%surfacewater%inqdra_in(level,node) - state%drainage%qdra(level,node)*tc_dt  ! TC-6
            end if
          end if
          if (allocated(state%drainage%qdra)) then
            qdraincomp(node) = state%drainage%qdra(level,node) + qdraincomp(node)
          end if
        end do
      end do

      ! [SS-SWC S-2.12B] all legacy half-writes dropped — state%soilwater is canonical
      ! SS-ATM Phase 2 Task A-2.2: aintcdt read from state%atmosphere (atmosphere home).
      state%soilwater%iintc = state%soilwater%iintc + (state%atmosphere%aintcdt+gird-nird)*tc_dt  ! S-2.12B, TC-6

      state%atmosphere%intr%iptra = state%atmosphere%intr%iptra + ptrats
      state%atmosphere%intr%ipeva = state%atmosphere%intr%ipeva + pevats
      state%atmosphere%intr%ievap = state%atmosphere%intr%ievap + revats
      ! SS-BND Phase 2 Task B-2.2: runots read from state%soilwater (boundary home).
      state%soilwater%iruno = state%soilwater%iruno + state%soilwater%runots                  ! S-2.12B
      state%soilwater%irunon = state%soilwater%irunon + state%soilwater%runon*tc_dt              ! S-2.12B, TC-6
      ! SS-ATM Phase 2 Task A-2.2: graidt/nraidt read from state%atmosphere (atmosphere home).
      state%soilwater%iprec = state%soilwater%iprec + (state%atmosphere%graidt+gird)*tc_dt       ! S-2.12B, TC-6
      state%atmosphere%intr%igrai = state%atmosphere%intr%igrai + state%atmosphere%graidt*tc_dt             ! TC-6
      state%soilwater%igird = state%soilwater%igird + gird*tc_dt                                  ! S-2.12B, TC-6
      state%atmosphere%intr%inrai = state%atmosphere%intr%inrai + state%atmosphere%nraidt*tc_dt             ! TC-6
      state%soilwater%inird = state%soilwater%inird + nird*tc_dt                                  ! S-2.12B, TC-6
      state%soilwater%iqbot = state%soilwater%iqbot + qbotts                                   ! S-2.12B
      if (state%soilwater%q(1) < 0.0d0) then
         state%soilwater%iqtdo = state%soilwater%iqtdo - state%soilwater%q(1)*tc_dt              ! S-2.12B, TC-6
      else
         state%soilwater%iqtup = state%soilwater%iqtup + state%soilwater%q(1)*tc_dt              ! S-2.12B, TC-6
      end if
      do node = 1, numnod+1
         if (state%soilwater%q(node) < 0.0d0) then
            state%soilwater%iqdo(node) = state%soilwater%iqdo(node) - state%soilwater%q(node)*tc_dt  ! S-2.12B, TC-6
         else
            state%soilwater%iqup(node) = state%soilwater%iqup(node) + state%soilwater%q(node)*tc_dt  ! S-2.12B, TC-6
         end if
      end do

      ! add time step fluxes to total cumulative values
      state%soilwater%cqssdi = state%soilwater%cqssdi + qssdisum*tc_dt                            ! S-2.12B, TC-6
      state%soilwater%cqrot  = state%soilwater%cqrot  + qrotts                                 ! S-2.12B
      ! SS-SWST Phase 2 Task 7: cqdra accumulated directly into state; global dropped.
      state%surfacewater%cqdra = state%surfacewater%cqdra + qdrats
      state%atmosphere%cumu%cptra = state%atmosphere%cumu%cptra + ptrats
      state%atmosphere%cumu%cpeva = state%atmosphere%cumu%cpeva + pevats
      state%atmosphere%cumu%cevap = state%atmosphere%cumu%cevap + revats
      if (state%soilwater%runots.lt.0.0d0) then
        state%soilwater%cinund = state%soilwater%cinund - state%soilwater%runots              ! S-2.12B
      else if (state%soilwater%runots.gt.0.0d0) then
        state%soilwater%crunoff = state%soilwater%crunoff + state%soilwater%runots            ! S-2.12B
      endif
      state%soilwater%irunoCN = state%soilwater%irunoCN + Runoff_CN*tc_dt                         ! S-2.12B, TC-6
      state%soilwater%crunoffCN = state%soilwater%crunoffCN + Runoff_CN*tc_dt                     ! S-2.12B, TC-6

      ! SS-ATM Phase 2 Task A-2.2: aintcdt/graidt/nraidt read from state%atmosphere (atmosphere home).
      state%atmosphere%cumu%caintc = state%atmosphere%cumu%caintc + (state%atmosphere%aintcdt+gird-nird)*tc_dt  ! TC-6

      state%atmosphere%cumu%cgrai = state%atmosphere%cumu%cgrai + state%atmosphere%graidt*tc_dt             ! TC-6
      state%atmosphere%cumu%cnrai = state%atmosphere%cumu%cnrai + state%atmosphere%nraidt*tc_dt             ! TC-6
!      cnrai = cgrai - caintc
      state%soilwater%cgird = state%soilwater%cgird + gird*tc_dt                                  ! S-2.12B, TC-6
      state%soilwater%cnird = state%soilwater%cnird + nird*tc_dt                                  ! S-2.12B, TC-6

      if (qbotts.lt.0.0d0) then
        state%soilwater%cqbotdo = state%soilwater%cqbotdo - qbotts                            ! S-2.12B
      else if (qbotts.gt.0.0d0) then
        state%soilwater%cqbotup = state%soilwater%cqbotup + qbotts                            ! S-2.12B
      endif
      state%soilwater%cqbot = state%soilwater%cqbot + qbotts                                   ! S-2.12B
      ! SS-SWST Phase 2 Task 7: cqdrain/in/out accumulated directly into state; globals dropped.
      if (allocated(state%surfacewater%cqdrain)) then
        do level = 1,nrlevs
          ! infiltration
          if (state%drainage%qdrain(level).lt.0.0d0) then
            state%surfacewater%cqdrainin(level) = state%surfacewater%cqdrainin(level) - state%drainage%qdrain(level)*tc_dt  ! TC-6
          ! drainage
          else if (state%drainage%qdrain(level).gt.0.0d0) then
            state%surfacewater%cqdrainout(level) = state%surfacewater%cqdrainout(level) + state%drainage%qdrain(level)*tc_dt  ! TC-6
          endif
          state%surfacewater%cqdrain(level) = state%surfacewater%cqdrain(level) + state%drainage%qdrain(level)*tc_dt  ! TC-6
        enddo
      end if

      ! rain on the ponding surface — [SS-SWC S-2.12B] legacy half-writes dropped
      state%soilwater%cqprai = state%soilwater%cqprai + state%atmosphere%nraidt*tc_dt  ! S-2.12B, TC-6
      state%soilwater%crunon = state%soilwater%crunon + state%soilwater%runon*tc_dt    ! S-2.12B, TC-6
      if (state%soilwater%q(1).lt.0.0d0) then
        state%soilwater%cqtdo = state%soilwater%cqtdo - state%soilwater%q(1)*tc_dt     ! S-2.12B, TC-6
      else if (state%soilwater%q(1).gt.0.0d0) then
        state%soilwater%cqtup = state%soilwater%cqtup + state%soilwater%q(1)*tc_dt     ! S-2.12B, TC-6
      endif

      ! compensate water balance error of this time step during remaining day part
      ! cumulative water balance error
      ! SS-SWC S-2.4: cnird/crunon/crunoff/cqrot/cqbot/volini/volact/PondIni/pond/cqssdi/cqprai
      !               read cutover to state%soilwater flat fields
      ! [SS-SWC S-2.12B] write directly to state%soilwater%wbalance — legacy global retired
      if (swsnow.eq.0) then
        state%soilwater%wbalance = state%atmosphere%cumu%cnrai + state%soilwater%cnird           &
     &        + state%soilwater%crunon - state%soilwater%crunoff                             &
     &        - state%soilwater%cqrot - state%atmosphere%cumu%cevap                              &
     &        - state%surfacewater%cqdra                                          &
     &        + state%soilwater%cqbot + state%soilwater%volini                                   &
     &        - state%soilwater%volact + state%soilwater%pondini                                      &
     &        - state%soilwater%pond + state%soilwater%cqssdi
      else
         state%soilwater%wbalance = state%soilwater%cqprai + state%soilwater%cnird          &
     &        + state%atmosphere%cumu%cmelt                                                           &
     &        + state%soilwater%crunon - state%soilwater%crunoff                            &
     &        - state%soilwater%cqrot - state%atmosphere%cumu%cevap                              &
     &        - state%surfacewater%cqdra                                          &
     &        + state%soilwater%cqbot + state%soilwater%volini                                   &
     &        - state%soilwater%volact + state%soilwater%pondini                                      &
     &        - state%soilwater%pond + state%soilwater%cqssdi
      endif

      if (FlMacropore) state%soilwater%wbalance = state%soilwater%wbalance - cQMpOutDrRap -            &
     &                  (WaSrDm1 + WaSrDm2 - WaSrDm1Ini - WaSrDm2Ini)

      end associate  ! tc_dt, tc_flDayStart (TC-6); numnod [GR-BH C4]

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
      ! [SS-SWC S-2.12B] igird/inird/IPondBeg/iruno/irunon/pond retired — read from state%soilwater
      use variables, only: nrlevs,NumNodNew,IcTopMp,FlMacropore,outfil,pathwork,DZNew,    &
                           CritDevMasBal,IQInTopVrtDm1,IQInTopLatDm1,IQInTopVrtDm2, &
                           IQInTopLatDm2,ISsnowBeg,IWaSrDm1Beg,IWaSrDm2Beg,WaSrDm1,WaSrDm2, &
                           dev_cmb
      ! [SS-TC TC-6] DayCum read cut over to state%timecontrol%daycum
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
      ! SS-ATM A-2.6: Ssnow retired — read from state%atmosphere%ssnow
      ! [SS-SWC S-2.12B] IPondBeg/pond -> state%soilwater
      SrDif = state%soilwater%IPondBeg-state%soilwater%pond + ISsnowBeg-state%atmosphere%ssnow
      IQInTopPreDm= 0.d0
      IQInTopLatDm= 0.d0
      if (FlMacropore .and. IcTopMp.eq.1) then
         IQInTopPreDm= IQInTopVrtDm1 + IQInTopVrtDm2
         IQInTopLatDm= IQInTopLatDm1 + IQInTopLatDm2
      endif

      ! Deviation mass balance Ponding layer in cm
      ! [SS-SWC S-2.12B] igird/inird/iruno/irunon -> state%soilwater
      DevMasBalPnd = state%atmosphere%intr%igrai + state%atmosphere%intr%igsnow + state%soilwater%igird + state%soilwater%irunon + inqNew(1) + SrDif &
     &             - (state%atmosphere%intr%igrai-state%atmosphere%intr%inrai-state%atmosphere%intr%isnrai + state%soilwater%igird-state%soilwater%inird + state%atmosphere%intr%isubl + state%atmosphere%intr%ievap + state%soilwater%iruno) &
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
         do level=1,nrlevs
            QdraPrf = QdraPrf + InqdraNew(level,ic)
         enddo
         if (FlMacropore) then
            IQExcMtxDm1= IQExcMtxDm1 + IQExcMtxDm1CpNew(ic)
            IQExcMtxDm2= IQExcMtxDm2 + IQExcMtxDm2CpNew(ic)
         endif
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
         do level=1,nrlevs
            Qdra(ic) = Qdra(ic) + inqdraNew(level,ic)
         enddo

         ! Deviation mass balance Soil Compartments in cm
         DevMasBalCmp(ic) = inqNew(ic+1) + SrDif &
     &                    - (inqNew(ic) + inqrotNew(ic) + Qdra(ic))
         if (FlMacropore)  DevMasBalCmp(ic) = DevMasBalCmp(ic) + &
     &                     IQExcMtxDm1CpNew(ic) + IQExcMtxDm2CpNew(ic)

         ! Check mass balance against criteria
         if (abs(DevMasBalCmp(ic)).gt.CritDevMasBal) then
            FlWriteDev = .true.
            FlWriteDevCmp(ic) = .true.
         endif
 100  continue

      ! 4) Macropore domains Dm1 and Dm2
      if (FlMacropore) then
         IQOutDrRap= 0.d0
         do ic = 1, numnodnew
            IQOutDrRap= IQOutDrRap + IQOutDrRapCpNew(ic)
         enddo

         ! Deviation mass balance Macropore Domains in cm
         SrDif = IWaSrDm1Beg - WaSrDm1
         DevMasBalDm1= IQInTopLatDm1 + SrDif - (IQExcMtxDm1  + IQOutDrRap)
         if (IcTopMp.eq.1) DevMasBalDm1= DevMasBalDm1 + IQInTopVrtDm1

         SrDif = IWaSrDm2Beg - WaSrDm2
         DevMasBalDm2= IQInTopLatDm2 + SrDif - IQExcMtxDm2
         if (IcTopMp.eq.1) DevMasBalDm2= DevMasBalDm2 + IQInTopVrtDm2

         ! Check mass balance against criteria
         if (abs(DevMasBalDm1).gt.CritDevMasBal) then
            FlWriteDev = .true.
            FlWriteDevDm1 = .true.
         endif
         if (abs(DevMasBalDm2).gt.CritDevMasBal) then
            FlWriteDev = .true.
            FlWriteDevDm2 = .true.
         endif
      endif

      ! In case of deviations of mass balance open file 'xxxxx.dwb.csv'

      if (FlWriteDev .and. .not.FlOpenFileDev) then
         filnam = trim(pathwork)//trim(outfil)//'.dwb'
         call file_open(dev_cmb, filnam, 'replace', 'write')
         write(dev_cmb,1)
         if (FlMacropore) write(dev_cmb,2)
         FlOpenFileDev = .true.
      endif

      ! Write deviations of water balance Top system
      ! SS-ATM A-2.6: igrai/inrai/igsnow/isnrai/isubl/ievap retired — read from state%atmosphere%intr
      ! [SS-SWC S-2.12B] igird/irunon/inird/iruno/Pond/IPondBeg -> state%soilwater
      if (FlWriteDevPnd) write(dev_cmb,3) state%timecontrol%daycum, DevMasBalPnd, &  ! TC-6: daycum -> state%timecontrol
     &    state%atmosphere%intr%igrai, state%atmosphere%intr%igsnow, state%soilwater%igird, state%soilwater%irunon, state%atmosphere%intr%isnrai, &
     &    state%atmosphere%intr%igrai-state%atmosphere%intr%inrai, state%soilwater%igird-state%soilwater%inird, &
     &    state%atmosphere%intr%isubl, state%atmosphere%intr%ievap, state%soilwater%iruno, inqNew(1), state%soilwater%pond, state%soilwater%IPondBeg, state%atmosphere%ssnow, &
     &    ISsnowBeg,IQInTopPreDm, IQInTopLatDm

      ! Write deviations of water balance whole Profile
      if (FlWriteDevPrf) write(dev_cmb,4) state%timecontrol%daycum, DevMasBalPrf, &  ! TC-6
     &    inqNew(1), inqNew(NumNodNew+1), QrotPrf, QdraPrf, WaSrPrf, &
     &    WaSrPrfBeg, IQExcMtxDm1, IQExcMtxDm2

      ! Write deviations of water balance of Individual Soil Compartments
      do ic= 1, numnodnew
         if (FlWriteDevCmp(ic)) write(dev_cmb,5) state%timecontrol%daycum,ic,DevMasBalCmp(ic), &  ! TC-6
     &      inqNew(ic), inqNew(ic+1), inqrotNew(ic), Qdra(ic), WaSr(ic), &
     &      WaSrBeg(ic), IQExcMtxDm1CpNew(ic), IQExcMtxDm2CpNew(ic)
      enddo

      ! Write deviations of water balance Macropore Domains
      if (FlWriteDevDm1) write(dev_cmb,6) state%timecontrol%daycum, DevMasBalDm1, &  ! TC-6
     &   IQInTopVrtDm1, IQInTopLatDm1, IQExcMtxDm1, WaSrDm1, &
     &   IWaSrDm1Beg, IQOutDrRap
      if (FlWriteDevDm2) write(dev_cmb,7) state%timecontrol%daycum, DevMasBalDm2, &  ! TC-6
     &   IQInTopVrtDm2, IQInTopLatDm2, IQExcMtxDm2, WaSrDm2, &
     &   IWaSrDm2Beg
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
      ! SS-SWC S-2.4: theta/FrArMtrx/volact reads cut over to state%soilwater
      subroutine watstor (state)
      ! [SS-SWC S-2.12B] volm1/volact/theta/FrArMtrx retired from variables; all via state%soilwater
      ! [GR-BH C4] numnod/dz migrated to state%mesh; use variables no longer needed here
      use swap_state_mod, only: swap_state_t
      IMPLICIT NONE

      type(swap_state_t), intent(inout) :: state
      INTEGER i

      ! [GR-BH C4] mesh globals aliased via state%mesh
      associate( numnod => state%mesh%numnod, &  ! [GR-BH C4]
                 dz     => state%mesh%dz       )  ! [GR-BH C4]

      ! update soil profile water storage — [SS-SWC S-2.12B] legacy half-writes dropped
      state%soilwater%volm1  = state%soilwater%volact
      state%soilwater%volact = 0.0d0
      do 10 i = 1,numnod
        state%soilwater%volact = state%soilwater%volact + state%soilwater%theta(i)*dz(i)*state%soilwater%FrArMtrx(i)
 10   continue

      end associate  ! numnod/dz [GR-BH C4]

      return
      end

end module