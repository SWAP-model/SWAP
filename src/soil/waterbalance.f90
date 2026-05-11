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
      use variables, only: disnod,logf,swscre,swbotb,flmacropore,numnod,h,z,pond,t1900,  &
                           gwl,nodgwl,bpegwl,npegwl,pegwl,nodgwlflcpzo,gwlflcpzo,CritUndSatVol
      use swap_log, only: log_debug, to_str
      implicit none
      ! SS-BND B-2.7: state added to read state%soilwater%gwlinp (gwlinp global retired)
      ! SS-SWC S-1.6: intent(in) -> intent(inout) to allow gwl/nodgwl/pegwl/bpegwl/npegwl/gwlflcpzo/nodgwlflcpzo dual-writes
      ! SS-SWC S-2.4: h/pond reads cut over to state%soilwater; level/watertable pass state
      type(swap_state_t), intent(inout) :: state
      ! local
      integer   i, node, nodhlp, nodheq1
      logical   flsat,flunsat
      character(len=200) messag
      character(len=19) datexti

      ! S-2.4 ASSOCIATE: alias state%soilwater arrays for h/pond reads
      associate( sw_h => state%soilwater%h, sw_pond => state%soilwater%pond )

      ! set initial values
      gwl       = 999.0d0
      state%soilwater%gwl    = 999.0d0                        ! S-1.6 dual-write
      pegwl     = 999.0d0
      state%soilwater%pegwl  = 999.0d0                        ! S-1.6 dual-write
      flsat     = .false.
      nodgwl    = numnod+1
      state%soilwater%nodgwl = numnod+1                       ! S-1.6 dual-write
      nodhlp    = numnod
      nodheq1   = numnod

      ! search for groundwater table
      if (sw_h(numnod).ge.0.0d0) flsat  = .true.      ! S-2.4 read cutover: h -> sw_h

      node = numnod
      nodgwlflcpzo = numnod + 1
      state%soilwater%nodgwlflcpzo = numnod + 1               ! S-1.6 dual-write
      gwlflcpzo    = gwl
      state%soilwater%gwlflcpzo    = gwl                      ! S-1.6 dual-write
      do while (flsat .and. node.gt.1)
         node = node - 1
         if(swbotb.eq.1)then
            if (sw_h(node) .lt. 0.0d0) then                  ! S-2.4 read cutover
               gwl = z(node+1) + sw_h(node+1) / (sw_h(node+1)-sw_h(node)) * disnod(node+1)  ! S-2.4
               state%soilwater%gwl    = gwl                   ! S-1.6 dual-write
               flsat   =.false.
               nodgwl  = node
               state%soilwater%nodgwl = node                  ! S-1.6 dual-write
            endif
         else
            if (sw_h(node) .lt. 1.0d0 .and. nodheq1.eq.numnod) nodheq1 = node  ! S-2.4

            if (sw_h(node) .lt. 0.0d0) then                  ! S-2.4 read cutover
               if (.not.flmacropore) then
                  flsat  = .false.
                  nodgwl = node
                  state%soilwater%nodgwl = node                ! S-1.6 dual-write
                  gwl    = level (state,1,node,nodheq1)        ! S-2.4: pass state to level
                  state%soilwater%gwl    = gwl                 ! S-1.6 dual-write
               elseif (flmacropore) then
                  if (gwl.gt.990.0d0) then
                     nodgwl = node
                     state%soilwater%nodgwl = node             ! S-1.6 dual-write
                     gwl    = level (state,2,node,nodheq1)     ! S-2.4: pass state to level
                     state%soilwater%gwl    = gwl              ! S-1.6 dual-write
                  endif
                  call watertable (state,node,nodgwlflcpzo,nodhlp,nodheq1,0.0d0,flsat,gwlflcpzo)  ! S-2.4
                  state%soilwater%nodgwlflcpzo = nodgwlflcpzo  ! S-1.6 dual-write (watertable out-arg)
                  state%soilwater%gwlflcpzo    = gwlflcpzo     ! S-1.6 dual-write (watertable out-arg)
               endif
            endif
         endif
      end do

      ! whole profile saturated, then add ponding layer to groundwater level
      if (flsat)then
         if(sw_h(1) .gt. 0.0d0)then                          ! S-2.4 read cutover
            if (sw_pond .lt. 1.d-8) then                     ! S-2.4 read cutover: pond -> sw_pond
               gwl = min(z(1)+sw_h(1),sw_pond)               ! S-2.4
            else
               gwl = sw_pond                                  ! S-2.4
            endif
         else
            gwl = 0.0d0
         end if
         state%soilwater%gwl    = gwl                         ! S-1.6 dual-write
         nodgwl = 1
         state%soilwater%nodgwl = 1                           ! S-1.6 dual-write
         if (flmacropore) then
            nodgwlflcpzo = 1
            state%soilwater%nodgwlflcpzo = 1                  ! S-1.6 dual-write
            gwlflcpzo    = gwl
            state%soilwater%gwlflcpzo    = gwl                ! S-1.6 dual-write
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
         bpegwl = i
         state%soilwater%bpegwl = i                           ! S-1.6 dual-write
         node   = bpegwl
         nodheq1 = bpegwl

         do while (flsat .and. node.gt.1)
            node = node - 1

            if (sw_h(node) .lt. 1.0d0 .and. nodheq1.eq.bpegwl) nodheq1 = node  ! S-2.4

            if (sw_h(node) .lt. 0.0d0) then                  ! S-2.4 read cutover
               if (.not.flmacropore) then
                  flsat = .false.
                  npegwl = node
                  state%soilwater%npegwl = node                ! S-1.6 dual-write
                  pegwl  = level (state,1,node,nodheq1)        ! S-2.4: pass state to level
                  state%soilwater%pegwl  = pegwl               ! S-1.6 dual-write
               elseif (flmacropore) then
                  call watertable (state,node,npegwl,nodhlp,nodheq1,CritUndSatVol,flsat,pegwl)  ! S-2.4
                  state%soilwater%npegwl = npegwl              ! S-1.6 dual-write (watertable out-arg)
                  state%soilwater%pegwl  = pegwl               ! S-1.6 dual-write (watertable out-arg)
               endif
            endif
         end do

         ! whole profile saturated, then add ponding layer to perched groundwater level
         if (flsat)then
            if(sw_h(1) .gt. 0.0d0)then                       ! S-2.4 read cutover
               if (sw_pond .lt. 1.d-8) then                  ! S-2.4 read cutover
                  pegwl = min(z(1)+sw_h(1),sw_pond)          ! S-2.4
               else
                  pegwl = sw_pond                             ! S-2.4
               endif
            else
               pegwl = 0.0d0
            end if
            state%soilwater%pegwl  = pegwl                    ! S-1.6 dual-write
            npegwl = 1
            state%soilwater%npegwl = 1                        ! S-1.6 dual-write
         endif
      else
         bpegwl = -1
         state%soilwater%bpegwl = -1                          ! S-1.6 dual-write
         npegwl = -1
         state%soilwater%npegwl = -1                          ! S-1.6 dual-write
      endif

      end associate  ! sw_h, sw_pond (S-2.4)

      ! fatal error if gwl below profile and flux has to be calculated
      if ((swbotb.eq.3.or.swbotb.eq.4).and.gwl.gt.998.0d0) then
          messag = 'The groundwater level descends below the lower' &
     &     //' boundary. This conflicts with bottom boundary' &
     &     //' condition 3 and 4. Extend soil profile!'
         call fatalerr_collected ('calcgwl',messag)
      endif

      ! warning error if there is inconsistency between defined gwl and soil physics
      if (swbotb.eq.1 .and. (state%soilwater%gwlinp .ge.z(1) .or. gwl.gt.998.0d0)) then
         ! determine date and date-time
         call dtdpst('year-month-day,hour:minute:seconds',t1900,datexti)
         write(messag,'(6a)')                                           &
     &         'No groundwater level because unsaturation at bottom ',  &
     &         'compartment ( ', datexti,  ' ). ',                      &
     &         'This is caused by inconsistency between ',              &
     &         'given gwl and soil physical parameters '
         call warn ('Calcgwl',messag,logf,swscre)
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
      use variables, only: numnod, disnod, dz, h, z, zbotcp
      implicit none

      type(swap_state_t), intent(in) :: state
      integer node, nodheq1, swoptlev
      integer i
      real(8) levm1, levp1
      real(8) level

      ! S-2.4 ASSOCIATE: alias state%soilwater%h for reads inside level
      associate( sw_h => state%soilwater%h )

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

      end associate  ! sw_h (S-2.4)

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
      use variables, only: numnod,dz,h,Theta,ThetaS,z
      implicit none

      type(swap_state_t), intent(in) :: state
      integer node,nodhlp,nodheq1,nodlev
      logical flsat
      real(8) waterlevel
      integer i
      real(8) CritUndSatVol, TotUndSatVol
      logical flsat2

      ! S-2.4 ASSOCIATE: alias state%soilwater arrays for h/Theta/ThetaS reads
      associate( sw_h     => state%soilwater%h,     &
                 sw_theta => state%soilwater%theta,  &
                 sw_thetas => state%soilwater%thetas )

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

      end associate  ! sw_h, sw_theta, sw_thetas (S-2.4)

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
      ! SS-CRP Phase 2 Task C-2.2: qrot/qrosum removed from use variables; read from state%soilwater.
      use variables, only: q,dt,inq,numnod,thetm1,theta,dz,qimmob,volact,volm1,swbotb,     &
                           FrArMtrx,QExcMpMtx,QMaPo,nrlevs,fllowgwl,qssdi, qssdisum
      use swap_state_mod, only: swap_state_t
      implicit none

      type(swap_state_t), intent(inout) :: state
      integer i,level

      ! S-2.4 ASSOCIATE: alias state%soilwater arrays/scalars for read cutover
      associate( sw_theta   => state%soilwater%theta,   &
                 sw_thetm1  => state%soilwater%thetm1,  &
                 sw_FrArMtrx => state%soilwater%FrArMtrx, &
                 sw_q       => state%soilwater%q,       &
                 sw_inq     => state%soilwater%intr%inq, &
                 sw_volact  => state%soilwater%volact,  &
                 sw_volm1   => state%soilwater%volm1,   &
                 sw_fllowgwl => state%soilwater%fllowgwl )

      ! determine qbot if not specified
      ! SS-BND B-2.7: qtop and qbot read/written via state%soilwater (globals retired)
      if (swbotb .eq. 5 .or. swbotb .eq. 7 .or.                         &
     &    swbotb .eq. 8 .or. swbotb .eq. -2 .or.                        &
     &    (swbotb .eq. 1 .and. sw_fllowgwl)) then        ! S-2.4 read cutover: fllowgwl -> sw_fllowgwl
        ! SS-CRP Phase 2 Task C-2.2: qrosum -> state%soilwater%qrosum
        state%soilwater%qbot = state%soilwater%qtop + state%soilwater%qrosum + state%surfacewater%qdrtot - QMaPo + (sw_volact-sw_volm1)/dt - qssdisum  ! S-2.4
      endif

      ! calculate fluxes (cm/d) from changes in volume per compartment
      i = numnod+1
      q(i) = state%soilwater%qbot
      sw_q(i)              = q(i)                             ! S-1.6 dual-write (via ASSOCIATE)
      inq(i) = sw_inq(i) + q(i)*dt                           ! S-2.4 read cutover: inq(i) -> sw_inq(i)
      state%soilwater%intr%inq(i) = inq(i)                   ! S-1.6 dual-write (snapshot after accumulation)
      do i = numnod,1,-1
        ! SS-CRP Phase 2 Task C-2.2: qrot(i) -> state%soilwater%qrot(i)
        ! S-2.4: theta/thetm1/FrArMtrx -> sw_theta/sw_thetm1/sw_FrArMtrx; q(i+1) -> sw_q(i+1)
        q(i) = - (sw_theta(i)-sw_thetm1(i)+qimmob(i))*sw_FrArMtrx(i)*dz(i)/dt +  &
     &                sw_q(i+1)-state%soilwater%qrot(i)+QExcMpMtx(i)+qssdi(i)

        if (allocated(state%drainage%qdra)) then
          do level=1,nrlevs
             q(i) = q(i) - state%drainage%qdra(level,i)
          enddo
        end if
        sw_q(i)              = q(i)                          ! S-1.6 dual-write
        inq(i) = sw_inq(i) + q(i)*dt                        ! S-2.4 read cutover: inq(i) -> sw_inq(i)
        state%soilwater%intr%inq(i) = inq(i)                ! S-1.6 dual-write (snapshot after accumulation)
      end do

      end associate  ! sw_theta etc (S-2.4)

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
      implicit none

      type(swap_state_t), intent(inout) :: state
      integer node,level
      real(8) qrotts,qdrats,ptrats,pevats,revats,qbotts
             
      if (flzerointr) then
        ! SS-ATM Phase 2 Task A-2.2 (D6): igrai/inrai removed — canonical reset
        ! is state%atmosphere%intr%reset() invoked in meteoday's ResetMetFlx (A-2.1).
        iprec = 0.d0
        igird = 0.d0
        inird = 0.d0
      endif

      ! potential transpiration of this timestep
      ! SS-ATM Phase 2 Task A-2.2: ptra read from state%atmosphere (atmosphere home).
      ptrats = state%atmosphere%ptra * dt

      ! potential soil evaporation of this timestep
      ! SS-ATM Phase 2 Task A-2.2: peva read from state%atmosphere (atmosphere home).
      pevats = state%atmosphere%peva * dt

      ! reduced soil evaporation of this timestep
      ! SS-BND Phase 2 Task B-2.2: reva read from state%soilwater (boundary home).
      revats = state%soilwater%reva * dt

      ! flux lower boundary of this timestep
      ! SS-BND Phase 2 Task B-2.2: qbot read from state%soilwater (boundary home).
      qbotts = state%soilwater%qbot*dt

      ! total root extraction of this timestep
      ! SS-CRP Phase 2 Task C-2.2: qrosum -> state%soilwater%qrosum
      qrotts = state%soilwater%qrosum * dt

      ! total drainage flux of this timestep
      ! SS-SWST Phase 2 Task 11 A2: qdrtot removed from globals; read from state.
      qdrats = state%surfacewater%qdrtot * dt

      ! determine daily actual transpiration
      if (fldaystart) tra = 0.0d0
      tra = tra + qrotts
      state%soilwater%intr%tra = state%soilwater%intr%tra + qrotts     ! S-1.5 dual-write (per-day)

      ! add time step fluxes to intermediate totals
      iqrot = iqrot + qrotts
      state%soilwater%intr%iqrot = state%soilwater%intr%iqrot + qrotts  ! S-1.5 dual-write
      do node = 1,noddrz
        ! SS-CRP Phase 2 Task C-2.2: qrot/qpotrot/qredwet/qreddry/qredsol/qredfrs -> state%soilwater
        inqrot(node) = inqrot(node) + state%soilwater%qrot(node) * dt
        state%soilwater%intr%inqrot(node) = state%soilwater%intr%inqrot(node) + state%soilwater%qrot(node) * dt  ! S-1.5 dual-write
        qpotrot_day(node) = qpotrot_day(node) + state%soilwater%qpotrot(node) * dt
        state%soilwater%intr%qpotrot_day(node) = state%soilwater%intr%qpotrot_day(node) + state%soilwater%qpotrot(node) * dt  ! S-1.5 dual-write (per-day)
        qredtot_day(node) = qredtot_day(node) + (state%soilwater%qredwet(node) + state%soilwater%qreddry(node) + state%soilwater%qredsol(node) + state%soilwater%qredfrs(node)) * dt
        state%soilwater%intr%qredtot_day(node) = state%soilwater%intr%qredtot_day(node) + (state%soilwater%qredwet(node) + state%soilwater%qreddry(node) + state%soilwater%qredsol(node) + state%soilwater%qredfrs(node)) * dt  ! S-1.5 dual-write (per-day)
      end do
      do node = 1,numnod
        inqssdi(node) = inqssdi(node) + qssdi(node) * dt
        state%soilwater%intr%inqssdi(node) = state%soilwater%intr%inqssdi(node) + qssdi(node) * dt  ! S-1.5 dual-write
        iqssdi = iqssdi + qssdi(node) * dt
        state%soilwater%intr%iqssdi = state%soilwater%intr%iqssdi + qssdi(node) * dt               ! S-1.5 dual-write
      end do
      ! SS-CRP Phase 2 Task C-2.2: qredXXXsum -> state%soilwater%qredXXXsum
      iqredwet = iqredwet + state%soilwater%qredwetsum*dt
      state%soilwater%intr%iqredwet = state%soilwater%intr%iqredwet + state%soilwater%qredwetsum*dt  ! S-1.5 dual-write
      iqreddry = iqreddry + state%soilwater%qreddrysum*dt
      state%soilwater%intr%iqreddry = state%soilwater%intr%iqreddry + state%soilwater%qreddrysum*dt  ! S-1.5 dual-write
      iqredsol = iqredsol + state%soilwater%qredsolsum*dt
      state%soilwater%intr%iqredsol = state%soilwater%intr%iqredsol + state%soilwater%qredsolsum*dt  ! S-1.5 dual-write
      iqredfrs = iqredfrs + state%soilwater%qredfrssum*dt
      state%soilwater%intr%iqredfrs = state%soilwater%intr%iqredfrs + state%soilwater%qredfrssum*dt  ! S-1.5 dual-write
      iqredwet_day = iqredwet_day + state%soilwater%qredwetsum*dt
      state%soilwater%intr%iqredwet_day = state%soilwater%intr%iqredwet_day + state%soilwater%qredwetsum*dt  ! S-1.5 dual-write (per-day)
      iqreddry_day = iqreddry_day + state%soilwater%qreddrysum*dt
      state%soilwater%intr%iqreddry_day = state%soilwater%intr%iqreddry_day + state%soilwater%qreddrysum*dt  ! S-1.5 dual-write (per-day)
      iqredsol_day = iqredsol_day + state%soilwater%qredsolsum*dt
      state%soilwater%intr%iqredsol_day = state%soilwater%intr%iqredsol_day + state%soilwater%qredsolsum*dt  ! S-1.5 dual-write (per-day)
      iqredfrs_day = iqredfrs_day + state%soilwater%qredfrssum*dt
      state%soilwater%intr%iqredfrs_day = state%soilwater%intr%iqredfrs_day + state%soilwater%qredfrssum*dt  ! S-1.5 dual-write (per-day)
      ! SS-ATM Phase 2 Task A-2.2: ptra read from state%atmosphere (atmosphere home).
      iptra_day    = iptra_day    + state%atmosphere%ptra * dt
      state%soilwater%intr%iptra_day = state%soilwater%intr%iptra_day + state%atmosphere%ptra * dt  ! S-1.5 dual-write (per-day)
      ies0 = ies0 + 0.1d0*es0*dt
      state%soilwater%intr%ies0 = state%soilwater%intr%ies0 + 0.1d0*es0*dt  ! S-1.5 dual-write
      iet0 = iet0 + 0.1d0*et0*dt
      state%soilwater%intr%iet0 = state%soilwater%intr%iet0 + 0.1d0*et0*dt  ! S-1.5 dual-write
      iew0 = iew0 + 0.1d0*ew0*dt
      state%soilwater%intr%iew0 = state%soilwater%intr%iew0 + 0.1d0*ew0*dt  ! S-1.5 dual-write

      ! SS-SWST Phase 2 Task 7: iqdra/inqdra* accumulated directly into state; global dropped.
      ! ADR 0031 Phase 2 Task 5: qdra global deleted; read from state%drainage%qdra.
      state%surfacewater%intermediate%iqdra = state%surfacewater%intermediate%iqdra + qdrats + QRapDra*dt
      do node = 1,numnod
        qdraincomp(node) = 0.d0
        do level = 1,nrlevs
          if (allocated(state%surfacewater%intermediate%inqdra) .and. allocated(state%drainage%qdra)) then
            state%surfacewater%intermediate%inqdra(level,node) = state%surfacewater%intermediate%inqdra(level,node) + state%drainage%qdra(level,node)*dt
            if (state%drainage%qdra(level,node) > 0.0d0) then
               state%surfacewater%intermediate%inqdra_out(level,node) = state%surfacewater%intermediate%inqdra_out(level,node) + state%drainage%qdra(level,node)*dt
            else
               state%surfacewater%intermediate%inqdra_in(level,node)  = state%surfacewater%intermediate%inqdra_in(level,node) - state%drainage%qdra(level,node)*dt
            end if
          end if
          if (allocated(state%drainage%qdra)) then
            qdraincomp(node) = state%drainage%qdra(level,node) + qdraincomp(node)
          end if
        end do
      end do

      ! SS-ATM Phase 2 Task A-2.2: aintcdt read from state%atmosphere (atmosphere home).
      iintc = iintc + (state%atmosphere%aintcdt+gird-nird)*dt
      state%soilwater%intr%iintc = state%soilwater%intr%iintc + (state%atmosphere%aintcdt+gird-nird)*dt  ! S-1.5 dual-write

      state%atmosphere%intr%iptra = state%atmosphere%intr%iptra + ptrats
      state%atmosphere%intr%ipeva = state%atmosphere%intr%ipeva + pevats
      state%atmosphere%intr%ievap = state%atmosphere%intr%ievap + revats
      ! SS-BND Phase 2 Task B-2.2: runots read from state%soilwater (boundary home).
      iruno = iruno + state%soilwater%runots
      state%soilwater%intr%iruno = state%soilwater%intr%iruno + state%soilwater%runots  ! S-1.5 dual-write
      ! SS-SWC S-2.4: runon read cutover to state%soilwater%runon
      irunon = irunon + state%soilwater%runon*dt              ! S-2.4 read cutover: runon -> state%soilwater%runon
      state%soilwater%intr%irunon = state%soilwater%intr%irunon + state%soilwater%runon*dt  ! S-1.5 dual-write
      ! SS-ATM Phase 2 Task A-2.2: graidt/nraidt read from state%atmosphere (atmosphere home).
      iprec = iprec + (state%atmosphere%graidt+gird)*dt
      state%soilwater%intr%iprec = state%soilwater%intr%iprec + (state%atmosphere%graidt+gird)*dt  ! S-1.5 dual-write
      state%atmosphere%intr%igrai = state%atmosphere%intr%igrai + state%atmosphere%graidt*dt
      igird = igird + gird*dt
      state%soilwater%intr%igird = state%soilwater%intr%igird + gird*dt  ! S-1.5 dual-write
      state%atmosphere%intr%inrai = state%atmosphere%intr%inrai + state%atmosphere%nraidt*dt
      inird = inird + nird*dt
      state%soilwater%intr%inird = state%soilwater%intr%inird + nird*dt  ! S-1.5 dual-write
      iqbot = iqbot + qbotts
      state%soilwater%intr%iqbot = state%soilwater%intr%iqbot + qbotts  ! S-1.5 dual-write
      ! SS-SWC S-2.4: q(1)/q(node) read cutover to state%soilwater%q
      if (state%soilwater%q(1) < 0.0d0) then                ! S-2.4 read cutover: q(1) -> state%soilwater%q(1)
         iqtdo = iqtdo - state%soilwater%q(1)*dt             ! S-2.4
         state%soilwater%intr%iqtdo = state%soilwater%intr%iqtdo - state%soilwater%q(1)*dt  ! S-1.5 dual-write
      else
         iqtup = iqtup + state%soilwater%q(1)*dt             ! S-2.4
         state%soilwater%intr%iqtup = state%soilwater%intr%iqtup + state%soilwater%q(1)*dt  ! S-1.5 dual-write
      end if
      do node = 1, numnod+1
         if (state%soilwater%q(node) < 0.0d0) then          ! S-2.4 read cutover: q(node) -> state%soilwater%q
            iqdo(node) = iqdo(node) - state%soilwater%q(node)*dt          ! S-2.4
            state%soilwater%intr%iqdo(node) = state%soilwater%intr%iqdo(node) - state%soilwater%q(node)*dt  ! S-1.5 dual-write
         else
            iqup(node) = iqup(node) + state%soilwater%q(node)*dt          ! S-2.4
            state%soilwater%intr%iqup(node) = state%soilwater%intr%iqup(node) + state%soilwater%q(node)*dt  ! S-1.5 dual-write
         end if
      end do

      ! add time step fluxes to total cumulative values
      cqssdi = cqssdi + qssdisum*dt
      state%soilwater%cumu%cqssdi = state%soilwater%cumu%cqssdi + qssdisum*dt  ! S-1.5 dual-write
      cqrot = cqrot + qrotts
      state%soilwater%cumu%cqrot = state%soilwater%cumu%cqrot + qrotts          ! S-1.5 dual-write
      ! SS-SWST Phase 2 Task 7: cqdra accumulated directly into state; global dropped.
      state%surfacewater%drainage_cumulative%cqdra = state%surfacewater%drainage_cumulative%cqdra + qdrats
      state%atmosphere%cumu%cptra = state%atmosphere%cumu%cptra + ptrats
      state%atmosphere%cumu%cpeva = state%atmosphere%cumu%cpeva + pevats
      state%atmosphere%cumu%cevap = state%atmosphere%cumu%cevap + revats
      if (state%soilwater%runots.lt.0.0d0) then
        cinund = cinund - state%soilwater%runots
        state%soilwater%cumu%cinund = state%soilwater%cumu%cinund - state%soilwater%runots  ! S-1.5 dual-write
      else if (state%soilwater%runots.gt.0.0d0) then
        crunoff = crunoff + state%soilwater%runots
        state%soilwater%cumu%crunoff = state%soilwater%cumu%crunoff + state%soilwater%runots  ! S-1.5 dual-write
      endif
      irunoCN = irunoCN + Runoff_CN*dt
      state%soilwater%intr%irunoCN = state%soilwater%intr%irunoCN + Runoff_CN*dt  ! S-1.5 dual-write (intr)
      crunoffCN = crunoffCN + Runoff_CN*dt
      state%soilwater%cumu%crunoffCN = state%soilwater%cumu%crunoffCN + Runoff_CN*dt  ! S-1.5 dual-write

      ! SS-ATM Phase 2 Task A-2.2: aintcdt/graidt/nraidt read from state%atmosphere (atmosphere home).
      state%atmosphere%cumu%caintc = state%atmosphere%cumu%caintc + (state%atmosphere%aintcdt+gird-nird)*dt

      state%atmosphere%cumu%cgrai = state%atmosphere%cumu%cgrai + state%atmosphere%graidt*dt
      state%atmosphere%cumu%cnrai = state%atmosphere%cumu%cnrai + state%atmosphere%nraidt*dt
!      cnrai = cgrai - caintc
      cgird = cgird + gird*dt
      state%soilwater%cumu%cgird = state%soilwater%cumu%cgird + gird*dt  ! S-1.5 dual-write
      cnird = cnird + nird*dt
      state%soilwater%cumu%cnird = state%soilwater%cumu%cnird + nird*dt  ! S-1.5 dual-write

      if (qbotts.lt.0.0d0) then
        cqbotdo = cqbotdo - qbotts
        state%soilwater%cumu%cqbotdo = state%soilwater%cumu%cqbotdo - qbotts  ! S-1.5 dual-write
      else if (qbotts.gt.0.0d0) then
        cqbotup = cqbotup + qbotts
        state%soilwater%cumu%cqbotup = state%soilwater%cumu%cqbotup + qbotts  ! S-1.5 dual-write
      endif
      cqbot = cqbot + qbotts
      state%soilwater%cumu%cqbot = state%soilwater%cumu%cqbot + qbotts  ! S-1.5 dual-write
      ! SS-SWST Phase 2 Task 7: cqdrain/in/out accumulated directly into state; globals dropped.
      if (allocated(state%surfacewater%drainage_cumulative%cqdrain)) then
        do level = 1,nrlevs
          ! infiltration
          if (state%drainage%qdrain(level).lt.0.0d0) then
            state%surfacewater%drainage_cumulative%cqdrainin(level) = state%surfacewater%drainage_cumulative%cqdrainin(level) - state%drainage%qdrain(level)*dt
          ! drainage
          else if (state%drainage%qdrain(level).gt.0.0d0) then
            state%surfacewater%drainage_cumulative%cqdrainout(level) = state%surfacewater%drainage_cumulative%cqdrainout(level) + state%drainage%qdrain(level)*dt
          endif
          state%surfacewater%drainage_cumulative%cqdrain(level) = state%surfacewater%drainage_cumulative%cqdrain(level) + state%drainage%qdrain(level)*dt
        enddo
      end if

      ! rain on the ponding surface
      ! SS-ATM Phase 2 Task A-2.2: nraidt read from state%atmosphere (atmosphere home).
      cqprai = cqprai + state%atmosphere%nraidt*dt
      state%soilwater%cumu%cqprai = state%soilwater%cumu%cqprai + state%atmosphere%nraidt*dt  ! S-1.5 dual-write
      ! SS-SWC S-2.4: runon read cutover to state%soilwater%runon
      crunon = crunon + state%soilwater%runon*dt               ! S-2.4 read cutover: runon -> state%soilwater%runon
      state%soilwater%cumu%crunon = state%soilwater%cumu%crunon + state%soilwater%runon*dt  ! S-1.5 dual-write
      ! SS-SWC S-2.4: q(1) read cutover to state%soilwater%q(1)
      if (state%soilwater%q(1).lt.0.0d0) then                 ! S-2.4 read cutover
        cqtdo = cqtdo - state%soilwater%q(1)*dt               ! S-2.4
        state%soilwater%cumu%cqtdo = state%soilwater%cumu%cqtdo - state%soilwater%q(1)*dt  ! S-1.5 dual-write
      else if (state%soilwater%q(1).gt.0.0d0) then            ! S-2.4
        cqtup = cqtup + state%soilwater%q(1)*dt               ! S-2.4
        state%soilwater%cumu%cqtup = state%soilwater%cumu%cqtup + state%soilwater%q(1)*dt  ! S-1.5 dual-write
      endif

      ! compensate water balance error of this time step during remaining day part
      ! cumulative water balance error
      ! SS-SWC S-2.4: cnird/crunon/crunoff/cqrot/cqbot/volini/volact/PondIni/pond/cqssdi/cqprai
      !               read cutover to state%soilwater%cumu / state%soilwater flat fields
      if (swsnow.eq.0) then
        ! SS-ATM A-2.6: cnrai/cevap retired; read from state%atmosphere%cumu
        wbalance = state%atmosphere%cumu%cnrai + state%soilwater%cumu%cnird           &  ! S-2.4
     &        + state%soilwater%cumu%crunon - state%soilwater%cumu%crunoff             &  ! S-2.4
     &        - state%soilwater%cumu%cqrot - state%atmosphere%cumu%cevap              &  ! S-2.4
     &        - state%surfacewater%drainage_cumulative%cqdra                          &
     &        + state%soilwater%cumu%cqbot + state%soilwater%volini                   &  ! S-2.4
     &        - state%soilwater%volact + state%soilwater%pondini                      &  ! S-2.4
     &        - state%soilwater%pond + state%soilwater%cumu%cqssdi                       ! S-2.4
      else
         ! SS-ATM Phase 2 Task A-2.2: cmelt read from state%atmosphere%cumu (atmosphere home).
         ! SS-ATM A-2.6: cevap retired — read from state%atmosphere%cumu%cevap
         wbalance = state%soilwater%cumu%cqprai + state%soilwater%cumu%cnird          &  ! S-2.4
     &        + state%atmosphere%cumu%cmelt                                           &
     &        + state%soilwater%cumu%crunon - state%soilwater%cumu%crunoff            &  ! S-2.4
     &        - state%soilwater%cumu%cqrot - state%atmosphere%cumu%cevap             &  ! S-2.4
     &        - state%surfacewater%drainage_cumulative%cqdra                          &
     &        + state%soilwater%cumu%cqbot + state%soilwater%volini                   &  ! S-2.4
     &        - state%soilwater%volact + state%soilwater%pondini                      &  ! S-2.4
     &        - state%soilwater%pond + state%soilwater%cumu%cqssdi                       ! S-2.4
      endif

      if (FlMacropore) wbalance = wbalance - cQMpOutDrRap -            &
     &                  (WaSrDm1 + WaSrDm2 - WaSrDm1Ini - WaSrDm2Ini)

      state%soilwater%wbalance = wbalance                      ! S-1.6 dual-write

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
      use variables, only: DayCum,nrlevs,NumNodNew,IcTopMp,FlMacropore,outfil,pathwork,DZNew,    &
                           CritDevMasBal,igird,inird,IPondBeg,IQInTopVrtDm1,IQInTopLatDm1,IQInTopVrtDm2, &
                           IQInTopLatDm2,ISsnowBeg,iruno,irunon,pond,IWaSrDm1Beg,IWaSrDm2Beg,WaSrDm1,WaSrDm2, &
                           dev_cmb
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
      SrDif = IPondBeg-Pond + ISsnowBeg-state%atmosphere%ssnow
      IQInTopPreDm= 0.d0
      IQInTopLatDm= 0.d0
      if (FlMacropore .and. IcTopMp.eq.1) then
         IQInTopPreDm= IQInTopVrtDm1 + IQInTopVrtDm2
         IQInTopLatDm= IQInTopLatDm1 + IQInTopLatDm2
      endif

      ! Deviation mass balance Ponding layer in cm
      ! SS-ATM A-2.6: igrai/inrai/ievap/igSnow/isnrai/isubl retired — read from state%atmosphere
      DevMasBalPnd = state%atmosphere%intr%igrai + state%atmosphere%intr%igsnow + igird + irunon + inqNew(1) + SrDif &
     &             - (state%atmosphere%intr%igrai-state%atmosphere%intr%inrai-state%atmosphere%intr%isnrai + igird-inird + state%atmosphere%intr%isubl + state%atmosphere%intr%ievap + iruno) &
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
      if (FlWriteDevPnd) write(dev_cmb,3) daycum, DevMasBalPnd, &
     &    state%atmosphere%intr%igrai, state%atmosphere%intr%igsnow, igird, irunon, state%atmosphere%intr%isnrai, &
     &    state%atmosphere%intr%igrai-state%atmosphere%intr%inrai, igird-inird, &
     &    state%atmosphere%intr%isubl, state%atmosphere%intr%ievap, iruno, inqNew(1), Pond, IPondBeg, state%atmosphere%ssnow, &
     &    ISsnowBeg,IQInTopPreDm, IQInTopLatDm

      ! Write deviations of water balance whole Profile
      if (FlWriteDevPrf) write(dev_cmb,4) daycum, DevMasBalPrf, &
     &    inqNew(1), inqNew(NumNodNew+1), QrotPrf, QdraPrf, WaSrPrf, &
     &    WaSrPrfBeg, IQExcMtxDm1, IQExcMtxDm2

      ! Write deviations of water balance of Individual Soil Compartments
      do ic= 1, numnodnew
         if (FlWriteDevCmp(ic)) write(dev_cmb,5) daycum,ic,DevMasBalCmp(ic), &
     &      inqNew(ic), inqNew(ic+1), inqrotNew(ic), Qdra(ic), WaSr(ic), &
     &      WaSrBeg(ic), IQExcMtxDm1CpNew(ic), IQExcMtxDm2CpNew(ic)
      enddo

      ! Write deviations of water balance Macropore Domains
      if (FlWriteDevDm1) write(dev_cmb,6) daycum, DevMasBalDm1, &
     &   IQInTopVrtDm1, IQInTopLatDm1, IQExcMtxDm1, WaSrDm1, &
     &   IWaSrDm1Beg, IQOutDrRap
      if (FlWriteDevDm2) write(dev_cmb,7) daycum, DevMasBalDm2, &
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
      use variables, only: volm1,volact,numnod,theta,dz,FrArMtrx
      use swap_state_mod, only: swap_state_t
      IMPLICIT NONE

      type(swap_state_t), intent(inout) :: state               ! S-1.6: state added for volm1/volact dual-write
      INTEGER i

      ! update soil profile water storage
      ! S-2.4: volact read from state%soilwater%volact (read cutover)
      volm1 = state%soilwater%volact                          ! S-2.4 read cutover: volact -> state%soilwater%volact
      state%soilwater%volm1  = state%soilwater%volact         ! S-1.6 dual-write
      volact = 0.0d0
      state%soilwater%volact = 0.0d0                         ! S-2.4: reset state alongside global
      do 10 i = 1,numnod
        ! S-2.4: theta/FrArMtrx reads cut over to state%soilwater%theta/FrArMtrx
        volact = volact + state%soilwater%theta(i)*dz(i)*state%soilwater%FrArMtrx(i)  ! S-2.4
 10   continue
      state%soilwater%volact = volact                         ! S-1.6 dual-write

      return
      end

end module