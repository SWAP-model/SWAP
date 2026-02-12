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
    implicit none
    private
    public :: calcgwl, level, watertable, fluxes, integral, checkmassbal, watstor
   public :: calcgwl_state, fluxes_state, integral_state, watstor_state
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
      subroutine calcgwl ()
      use variables, only: disnod,logf,swscre,swbotb,flmacropore,numnod,gwlinp,h,z,pond,t1900,  &
                           gwl,nodgwl,bpegwl,npegwl,pegwl,nodgwlflcpzo,gwlflcpzo,CritUndSatVol
      use swap_log, only: log_debug, to_str
      implicit none
      ! local
      integer   i, node, nodhlp, nodheq1
      logical   flsat,flunsat
      character(len=200) messag
      character(len=19) datexti

      ! set initial values
      gwl       = 999.0d0
      pegwl     = 999.0d0
      flsat     = .false.
      nodgwl    = numnod+1
      nodhlp    = numnod
      nodheq1   = numnod

      ! search for groundwater table
      if (h(numnod).ge.0.0d0) flsat  = .true.

      node = numnod
      nodgwlflcpzo = numnod + 1
      gwlflcpzo    = gwl
      do while (flsat .and. node.gt.1)
         node = node - 1 
         if(swbotb.eq.1)then
            if (h(node) .lt. 0.0d0) then
               gwl = z(node+1) + h(node+1) / (h(node+1)-h(node)) * disnod(node+1)
               flsat   =.false.
               nodgwl  = node
            endif
         else
            if (h(node) .lt. 1.0d0 .and. nodheq1.eq.numnod) nodheq1 = node

            if (h(node) .lt. 0.0d0) then
               if (.not.flmacropore) then
                  flsat  = .false.
                  nodgwl = node
                  gwl    = level (1,node,nodheq1)
               elseif (flmacropore) then
                  if (gwl.gt.990.0d0) then
                     nodgwl = node
                     gwl    = level (2,node,nodheq1)
                  endif
                  call watertable (node,nodgwlflcpzo,nodhlp,nodheq1,0.0d0,flsat,gwlflcpzo)
               endif
            endif
         endif
      end do

      ! whole profile saturated, then add ponding layer to groundwater level
      if (flsat)then
         if(h(1) .gt. 0.0d0)then
            if (pond .lt. 1.d-8) then
               gwl = min(z(1)+h(1),pond)
            else
               gwl = pond
            endif
         else
            gwl = 0.0d0
         end if
         nodgwl = 1
         if (flmacropore) then
            nodgwlflcpzo = 1
            gwlflcpzo    = gwl
         endif         
      endif         
 
      ! search for perched groundwater table

      ! first, search for first saturated compartment (i) above groundwater level
      i = nodhlp
      flunsat = .true.
      do while (flunsat .and. i.ge.1)
         if (h(i).ge.0.0d0) flunsat = .false.
         i = i - 1 
      enddo

      ! if saturated compartment above gwl exists, then find perched groundwater table
      if (i.ne.0) then
         flsat  = .true.
         bpegwl = i
         node   = bpegwl
         nodheq1 = bpegwl

         do while (flsat .and. node.gt.1)
            node = node - 1

            if (h(node) .lt. 1.0d0 .and. nodheq1.eq.bpegwl) nodheq1 = node

            if (h(node) .lt. 0.0d0) then
               if (.not.flmacropore) then
                  flsat = .false.
                  npegwl = node
                  pegwl  = level (1,node,nodheq1)
               elseif (flmacropore) then
                  call watertable (node,npegwl,nodhlp,nodheq1,CritUndSatVol,flsat,pegwl)
               endif
            endif
         end do

         ! whole profile saturated, then add ponding layer to perched groundwater level
         if (flsat)then
            if(h(1) .gt. 0.0d0)then
               if (pond .lt. 1.d-8) then
                  pegwl = min(z(1)+h(1),pond)
               else
                  pegwl = pond
               endif
            else
               pegwl = 0.0d0
            end if
            npegwl = 1
         endif 
      else
         bpegwl = -1
         npegwl = -1
      endif

      ! fatal error if gwl below profile and flux has to be calculated
      if ((swbotb.eq.3.or.swbotb.eq.4).and.gwl.gt.998.0d0) then
          messag = 'The groundwater level descends below the lower' &
     &     //' boundary. This conflicts with bottom boundary' &
     &     //' condition 3 and 4. Extend soil profile!'
         call fatalerr ('calcgwl',messag)
      endif

      ! warning error if there is inconsistency between defined gwl and soil physics
      if (swbotb.eq.1 .and. (gwlinp .ge.z(1) .or. gwl.gt.998.0d0)) then
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
      function level (swoptlev,node,nodheq1)
      use variables, only: numnod, disnod, dz, h, z, zbotcp
      implicit none
      
      integer node, nodheq1, swoptlev
      integer i
      real(8) levm1, levp1
      real(8) level

      if (swoptlev.eq.1) then
         ! groundwater level equals elevation head where h = 0
         if (h(node+1).ge.0.0d0)then
            level = z(node+1) + h(node+1) / (h(node+1)-h(node)) * disnod(node+1)
         else   
            level = zbotcp(node) - h(node)
            level = min(z(node),max(zbotcp(node),level))
         end if

      elseif (swoptlev.eq.2) then
         ! groundwater level equals average of elevation heads of h = -1 and h = +1
         ! elevation head of h = +1
         i = nodheq1
         if (nodheq1.eq.numnod) then
            levp1 = z(i) - 0.5d0 * dz(i)
         else
            levp1 = z(i) - (z(i) - z(i+1)) * (1.d0-h(i)) / (h(i+1)-h(i))
         endif
         ! elevation head of h = -1
         i = node
         do while (h(i).gt.-1.d0 .and. i.gt.1)
            i = i - 1
         enddo
         if (i.eq.1 .and. h(1).gt.-1.d0 .and. node.gt.2) then
            ! no compartment with pressure head < -1 cm in top of profile:
            ! use elevation head of h = 0 as estimation for groundwater level
            levm1 = z(node+1) + h(node+1) / (h(node+1)-h(node)) * disnod(node+1)
            levp1 = levm1
         else
            levm1 = z(i+1) + (z(i) - z(i+1)) * (1.d0+h(i+1)) / (h(i+1)-h(i))
         endif
         ! groundwater level = average of levp1 and levm1
         level = (levp1 + levm1) / 2.d0
      endif

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
      subroutine watertable (node,nodlev,nodhlp,nodheq1,CritUndSatVol,flsat,waterlevel)
      use variables, only: numnod,dz,h,Theta,ThetaS,z
      implicit none

      integer node,nodhlp,nodheq1,nodlev
      logical flsat
      real(8) waterlevel
      integer i
      real(8) CritUndSatVol, TotUndSatVol
      logical flsat2
      TotUndSatVol = 0.0d0
      flsat2 = .false. 
      i = node
      do while (TotUndSatVol.lt.CritUndSatVol .and. .not.flsat2 .and. i.ge.1)
         TotUndSatVol = TotUndSatVol + (ThetaS(i) - Theta(i)) * dz(i)
         if (h(i).gt.-1.d-7) flsat2 = .true.
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
            waterlevel = level(1,node,nodheq1)
         else
            waterlevel = level(2,node,nodheq1)
         endif
         i = max(node-2,1)
         do while(z(i)-0.5d0*dz(i).gt.waterlevel .and. i.gt.2 .and. i.lt.numnod)
            i = i + 1
         enddo
         nodlev = min(max(i,1),numnod)
      endif

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
      subroutine fluxes ()
      use variables, only: q,qbot,dt,inq,numnod,thetm1,theta,dz,qrot,qdra,qimmob,qtop,qrosum,qdrtot,volact,volm1,swbotb,     &
                           FrArMtrx,QExcMpMtx,QMaPo,nrlevs,fllowgwl,qssdi, qssdisum
      implicit none

      integer i,level

      ! determine qbot if not specified
      if (swbotb .eq. 5 .or. swbotb .eq. 7 .or.                         &
     &    swbotb .eq. 8 .or. swbotb .eq. -2 .or.                        &
     &    (swbotb .eq. 1 .and. fllowgwl)) then
        qbot = qtop + qrosum + qdrtot - QMaPo + (volact-volm1)/dt - qssdisum
      endif

      ! calculate fluxes (cm/d) from changes in volume per compartment
      i = numnod+1
      q(i) = qbot
      inq(i) = inq(i) + q(i)*dt
      do i = numnod,1,-1
        q(i) = - (theta(i)-thetm1(i)+qimmob(i))*FrArMtrx(i)*dz(i)/dt +  &
     &                q(i+1)-qrot(i)+QExcMpMtx(i)+qssdi(i)
     
        do level=1,nrlevs
           q(i) = q(i) - qdra(level,i)
        enddo
        inq(i) = inq(i) + q(i)*dt
      end do

      return
      end

      !> State-aware wrapper for `calcgwl`
      !!
      !! Executes the legacy groundwater-level calculation and synchronizes
      !! selected soil water balance outputs.
      !!
      !! @param[inout] state SWAP model state container
      subroutine calcgwl_state(state)
      use swap_state_mod, only: swap_state_t
      use swap_state_sync, only: soilwaterbalance_outputs_from_variables
      implicit none

      type(swap_state_t), intent(inout) :: state

      call calcgwl()
      call soilwaterbalance_outputs_from_variables(state%soil, state%numnod)
      end subroutine calcgwl_state

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
      subroutine integral 
      Use Variables
      implicit none

      integer node,level
      real(8) qrotts,qdrats,ptrats,pevats,revats,qbotts
             
      if (flzerointr) then
        igrai = 0.d0
        inrai = 0.d0
        iprec = 0.d0
        igird = 0.d0
        inird = 0.d0
      endif

      ! potential transpiration of this timestep
      ptrats = ptra * dt

      ! potential soil evaporation of this timestep
      pevats = peva * dt

      ! reduced soil evaporation of this timestep
      revats = reva * dt

      ! flux lower boundary of this timestep
      qbotts = qbot*dt

      ! total root extraction of this timestep
      qrotts = qrosum * dt

      ! total drainage flux of this timestep
      qdrats = qdrtot * dt

      ! determine daily actual transpiration
      if (fldaystart) tra = 0.0d0
      tra = tra + qrotts

      ! add time step fluxes to intermediate totals
      iqrot = iqrot + qrotts
      do node = 1,noddrz
        inqrot(node) = inqrot(node) + qrot(node) * dt
        qpotrot_day(node) = qpotrot_day(node) + qpotrot(node) * dt
        qredtot_day(node) = qredtot_day(node) + (qredwet(node) + qreddry(node) + qredsol(node) + qredfrs(node)) * dt
      end do
      do node = 1,numnod
        inqssdi(node) = inqssdi(node) + qssdi(node) * dt
        iqssdi = iqssdi + qssdi(node) * dt
      end do
      iqredwet = iqredwet + qredwetsum*dt
      iqreddry = iqreddry + qreddrysum*dt
      iqredsol = iqredsol + qredsolsum*dt
      iqredfrs = iqredfrs + qredfrssum*dt
      iqredwet_day = iqredwet_day + qredwetsum*dt
      iqreddry_day = iqreddry_day + qreddrysum*dt
      iqredsol_day = iqredsol_day + qredsolsum*dt
      iqredfrs_day = iqredfrs_day + qredfrssum*dt
      iptra_day    = iptra_day    + ptra * dt
      ies0 = ies0 + 0.1d0*es0*dt
      iet0 = iet0 + 0.1d0*et0*dt
      iew0 = iew0 + 0.1d0*ew0*dt

      iqdra = iqdra + qdrats + QRapDra*dt
      do node = 1,numnod
        qdraincomp(node) = 0.d0
        do level = 1,nrlevs
          inqdra(level,node) = inqdra(level,node)+qdra(level,node)*dt
          if (qdra(level,node) > 0.0d0) then
             inqdra_out(level,node) = inqdra_out(level,node) + qdra(level,node)*dt
          else
             inqdra_in(level,node)  = inqdra_in(level,node) - qdra(level,node)*dt
          end if
          qdraincomp(node) = qdra(level,node) + qdraincomp(node)
        end do
      end do

      iintc = iintc + (aintcdt+gird-nird)*dt

      iptra = iptra + ptrats
      ipeva = ipeva + pevats
      ievap = ievap + revats
      iruno = iruno + runots
      irunon = irunon + runon*dt
      iprec = iprec + (graidt+gird)*dt
      igrai = igrai + graidt*dt
      igird = igird + gird*dt
      inrai = inrai + nraidt*dt
      inird = inird + nird*dt
      iqbot = iqbot + qbotts
      if (q(1) < 0.0d0) then
         iqtdo = iqtdo - q(1)*dt
      else
         iqtup = iqtup + q(1)*dt
      end if
      do node = 1, numnod+1
         if (q(node) < 0.0d0) then
            iqdo(node) = iqdo(node) - q(node)*dt
         else
            iqup(node) = iqup(node) + q(node)*dt
         end if
      end do

      ! add time step fluxes to total cumulative values
      cqssdi = cqssdi + qssdisum*dt
      cqrot = cqrot + qrotts
      cqdra = cqdra + qdrats
      cptra = cptra + ptrats
      cpeva = cpeva + pevats
      cevap = cevap + revats
      if (runots.lt.0.0d0) then
        cinund = cinund - runots
      else if (runots.gt.0.0d0) then
        crunoff = crunoff + runots
      endif
      irunoCN = irunoCN + Runoff_CN*dt
      crunoffCN = crunoffCN + Runoff_CN*dt

      caintc = caintc + (aintcdt+gird-nird)*dt

      cgrai = cgrai + graidt*dt
      cnrai = cnrai + nraidt*dt
!      cnrai = cgrai - caintc
      cgird = cgird + gird*dt
      cnird = cnird + nird*dt

      if (qbotts.lt.0.0d0) then
        cqbotdo = cqbotdo - qbotts
      else if (qbotts.gt.0.0d0) then
        cqbotup = cqbotup + qbotts
      endif
      cqbot = cqbot + qbotts
      do level = 1,nrlevs
        ! infiltration
        if (qdrain(level).lt.0.0d0) then
          cqdrainin(level) = cqdrainin(level) - qdrain(level)*dt
        ! drainage
        else if (qdrain(level).gt.0.0d0) then
          cqdrainout(level) = cqdrainout(level) + qdrain(level)*dt
        endif      
        cqdrain(level) = cqdrain(level) + qdrain(level)*dt
      enddo

      ! rain on the ponding surface
      cqprai = cqprai + nraidt*dt
      crunon = crunon + runon*dt
      if (q(1).lt.0.0d0) then
        cqtdo = cqtdo - q(1)*dt
      else if (q(1).gt.0.0d0) then
        cqtup = cqtup + q(1)*dt
      endif

      ! compensate water balance error of this time step during remaining day part
      ! cumulative water balance error
      if (swsnow.eq.0) then 
        wbalance = cnrai + cnird + crunon - crunoff - cqrot - cevap     &
     &        - cqdra + cqbot + volini - volact + PondIni - pond + cqssdi
      else
         wbalance = cqprai + cnird + cmelt + crunon - crunoff           &
     &        - cqrot - cevap - cqdra                                   &
     &        + cqbot + volini - volact + PondIni - pond + cqssdi
      endif

      if (FlMacropore) wbalance = wbalance - cQMpOutDrRap -            &
     &                  (WaSrDm1 + WaSrDm2 - WaSrDm1Ini - WaSrDm2Ini)

      return
      end

      !> State-aware wrapper for `fluxes`
      !!
      !! Executes the legacy flux integration between compartments and
      !! synchronizes selected soil water balance outputs.
      !!
      !! @param[inout] state SWAP model state container
      subroutine fluxes_state(state)
      use swap_state_mod, only: swap_state_t
      use swap_state_sync, only: soilwaterbalance_outputs_from_variables
      implicit none

      type(swap_state_t), intent(inout) :: state

      call fluxes()
      call soilwaterbalance_outputs_from_variables(state%soil, state%numnod)
      end subroutine fluxes_state

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
      subroutine checkmassbal (flopenfiledev,inqdranew,iqexcmtxdm1cpnew,iqexcmtxdm2cpnew,inqnew,iqoutdrrapcpnew,inqrotnew,ithetabegnew,thetanew)
      use variables, only: DayCum,nrlevs,NumNodNew,IcTopMp,FlMacropore,outfil,pathwork,DZNew,    &
                           CritDevMasBal,ievap,igird,igrai,igSnow,inird,inrai,IPondBeg,IQInTopVrtDm1,IQInTopLatDm1,IQInTopVrtDm2, &
                           IQInTopLatDm2,ISsnowBeg,iruno,irunon,isnrai,iSubl,pond,Ssnow,IWaSrDm1Beg,IWaSrDm2Beg,WaSrDm1,WaSrDm2, &
                           dev_cmb
      use swap_array_dimensions, only: macp, madr
      implicit none

      ! -   global
      real(8) IQExcMtxDm1CpNew(macp), IQExcMtxDm2CpNew(macp)
      real(8) inqdraNew(Madr,macp)
      real(8) inqNew(macp+1), IQOutDrRapCpNew(macp), inqrotNew(macp)
      real(8) IThetaBegNew(MaCp),thetaNew(macp)
      logical FlOpenFileDev

      ! local
      integer Level, getun, ic
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
      SrDif = IPondBeg-Pond + ISsnowBeg-Ssnow
      IQInTopPreDm= 0.d0
      IQInTopLatDm= 0.d0
      if (FlMacropore .and. IcTopMp.eq.1) then
         IQInTopPreDm= IQInTopVrtDm1 + IQInTopVrtDm2
         IQInTopLatDm= IQInTopLatDm1 + IQInTopLatDm2
      endif

      ! Deviation mass balance Ponding layer in cm
      DevMasBalPnd = igrai + igsnow + igird + irunon + inqNew(1) + SrDif &
     &             - (igrai-inrai-isnrai + igird-inird + isubl + ievap + iruno) &
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
         dev_cmb = getun (20,90)
         call fopens(dev_cmb,filnam,'new','del')
         write(dev_cmb,1)
         if (FlMacropore) write(dev_cmb,2)
         FlOpenFileDev = .true.
      endif

      ! Write deviations of water balance Top system
      if (FlWriteDevPnd) write(dev_cmb,3) daycum, DevMasBalPnd, &
     &    igrai, igsnow, igird, irunon, isnrai, igrai-inrai,igird-inird, &
     &    isubl,ievap, iruno, inqNew(1), Pond, IPondBeg, Ssnow, &
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

      !> State-aware wrapper for `integral`
      !!
      !! Executes the legacy cumulative/intermediate flux accounting and
      !! synchronizes selected soil water balance outputs.
      !!
      !! @param[inout] state SWAP model state container
      subroutine integral_state(state)
      use swap_state_mod, only: swap_state_t
      use swap_state_sync, only: soilwaterbalance_outputs_from_variables
      implicit none

      type(swap_state_t), intent(inout) :: state

      call integral()
      call soilwaterbalance_outputs_from_variables(state%soil, state%numnod)
      end subroutine integral_state

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
      subroutine watstor ()
      use variables, only: volm1,volact,numnod,theta,dz,FrArMtrx 
      IMPLICIT NONE

      INTEGER i

      ! update soil profile water storage
      volm1 = volact
      volact = 0.0d0
      do 10 i = 1,numnod
        volact = volact+theta(i)*dz(i)*FrArMtrx(i)
 10   continue

      return
      end

      !> State-aware wrapper for `watstor`
      !!
      !! Executes the legacy water-storage update and synchronizes selected
      !! soil water balance outputs.
      !!
      !! @param[inout] state SWAP model state container
      subroutine watstor_state(state)
      use swap_state_mod, only: swap_state_t
      use swap_state_sync, only: soilwaterbalance_outputs_from_variables
      implicit none

      type(swap_state_t), intent(inout) :: state

      call watstor()
      call soilwaterbalance_outputs_from_variables(state%soil, state%numnod)
      end subroutine watstor_state

end module