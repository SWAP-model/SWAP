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
      use Variables
      use array_utils, only: afgen
      use surfacewater_init_mod, only: surfacewater_init
      use swap_state_mod, only: swap_state_t
      implicit none
      integer,            intent(in)    :: task
      type(swap_state_t), intent(inout) :: state
      logical,            intent(out)   :: request_smaller_dt

!     local
      integer level, node
      real(8) zCum,difzTopDisLay(madr),ratio,ratiodz,sumqdr(madr),dh
      integer nodeTopDisLay(madr)
      character(len=300) messag

      request_smaller_dt = .false.

! ----------------------------------------------------------------------
      select case (task)
      case (1)

! === initialization ===================================================

! --- read input data
      call surfacewater_init (state)

      ! hwlman and vtair global writes dropped: only output reads them,
      ! via state%surfacewater%hwlman / state%surfacewater%vtair.
      state%surfacewater%hwlman = 0.0d0
      state%surfacewater%vtair  = 0.0d0

!   - In case of macropores: initialise drainage basis for rapid drainage through macropores
      if (state%surfacewater%flInitDraBas) then
         if (NumLevRapDra.gt.nrlevs) then
            messag = ' NUMLEVRAPDRA greater then NRLEVS'
            call fatalerr_collected('MacroRead',messag)
         endif
!
         if (swdtyp(NumLevRapDra).eq.1) then
            ZDraBas = zbotdr(NumLevRapDra)      ! drain tube
            state%surfacewater%ZDraBas = zbotdr(NumLevRapDra)
         elseif (Swsec.eq.1) then
            ZDraBas = afgen (wlstab,2*maowl,t1900) ! open drain, surf.wat. level input
            state%surfacewater%ZDraBas = ZDraBas
         elseif (Swsec.eq.2) then
            ZDraBas = WlStar                    ! open drain, srf.wat. level simulated
            state%surfacewater%ZDraBas = WlStar
         endif
!
         flInitDraBas = .false.
         state%surfacewater%flInitDraBas = .false.
!
         Return
!
      endif

      return

      case (2)

! === lateral drainage fluxes to surface water ========================

! --- reset intermediate surface water and drainage fluxes
      if (flzerointr) then
        do node = 1,numnod
          do level = 1,nrlevs
            inqdra(level,node)                        = 0.0d0
            state%surfacewater%inqdra(level,node)     = 0.0d0
            inqdra_in(level,node)                     = 0.0d0
            state%surfacewater%inqdra_in(level,node)  = 0.0d0
            inqdra_out(level,node)                    = 0.0d0
            state%surfacewater%inqdra_out(level,node) = 0.0d0
          enddo
        enddo
        iqdra = 0.0d0
        state%surfacewater%iqdra = 0.0d0
      endif

! --- reset cumulative surface water and drainage fluxes
      if (flzerocumu) then
        ! cqdrd/cwsupp/cwout global writes dropped: only output reads them,
        ! via state%surfacewater%*.  State resets remain authoritative.
        state%surfacewater%cqdrd  = 0.0d0
        state%surfacewater%cwsupp = 0.0d0
        state%surfacewater%cwout  = 0.0d0
        cqdra = 0.0d0
        state%surfacewater%cqdra = 0.0d0
        do level = 1,nrlevs
          cqdrain(level) = 0.0d0
          state%surfacewater%cqdrain(level) = 0.0d0
          cqdrainin(level) = 0.0d0
          state%surfacewater%cqdrainin(level) = 0.0d0
          cqdrainout(level) = 0.0d0
          state%surfacewater%cqdrainout(level) = 0.0d0
        enddo
      endif

! --- no drainage at all
      if (gwl.gt.998.0d0) then
        do level = 1,nrlevs
          qdrain(level) = 0.0d0
        end do
        return
      endif

! --- calculate lateral drainage
      call bocodre (dh)

! --- partition drainage flux over compartments

      do level=1,nrlevs
         do node = 1,numnod
            qdra(level,node) = 0.0d0
            state%surfacewater%qdra(level,node) = 0.0d0
         end do
      end do

      if (swdivd.eq.1) then
!cD        do level=1,nrlevs
!cD           qdrain_old(level) = qdrain(level)
!cD        end do
         call divdra (numnod,nrlevs,dz,ksatfit,ksatexm,fluseksatexm,    &
            layer,cofani,gwl,l,qdrain,qdra,Swdivdinf,Swnrsrf,           &
     &      SwTopnrsrf,Zbotdr,dt,FacDpthInf,owltab,t1900)

!        redistribute qdrain with new top boundary for discharge layers
         if(swdislay.eq.2) then
            do level=1,nrlevs
               if(swtopdislay(level).eq.1)  then
                  zTopDisLay(level) = fTopDisLay(level) * gwl  +        &
     &                       (1.0d0-fTopDisLay(level)) * (gwl-dh)
               end if
            end do
         end if
         if(swdislay.eq.1 .or. swdislay.eq.2) then
            do level=1,nrlevs
               if(swtopdislay(level).eq.1)  then
!                 find node nr of new top of discharge layer
                  nodeTopDisLay(level) = 1
                  zCum               = - dz(1)
                  do while (zTopDisLay(level) .lt. zCum)
                     nodeTopDisLay(level) = nodeTopDisLay(level) + 1
                     zCum              = zCum - dz(nodeTopDisLay(level))
                  enddo
!                 saturated part (difzTopDisLay(lev)) of compartment containing waterlevel
                  difzTopDisLay(level) = zTopDisLay(level) - zCum
                  ratiodz =                                             &
     &                     difzTopDisLay(level)/dz(nodeTopDisLay(level))
                  sumqdr(level) =                                       &
     &                        ratiodz * qdra(level,nodeTopDisLay(level))
                  do node = nodeTopDisLay(level)+1,numnod
                     sumqdr(level) =  sumqdr(level) + qdra(level,node)
                  end do
                  if( dabs(sumqdr(level)) .lt. 1.0d-8)then
                     ratio = 1.0d0
                  else
                     ratio = qdrain(level)/sumqdr(level)
                  end if
!                 redistribute drainwater fluxes
                  do node = 1,nodeTopDisLay(level)-1
                     qdra(level,node) = 0.0d0
                  end do
                  qdra(level,nodeTopDisLay(level)) =                    &
     &                qdra(level,nodeTopDisLay(level)) * ratio * ratiodz
                  do node = nodeTopDisLay(level)+1,numnod
                     qdra(level,node) = qdra(level,node)* ratio
                  end do
               endif
            end do
         endif

! --- error handling
!cD        do level=1,nrlevs
!cD           qdrain_new(level) = 0.0d0
!cD           do node = 1,numnod
!cD              qdrain_new(level) = qdrain_new(level) + qdra(level,node)
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

         ! dual-write qdra after divdra (and possible redistribution)
         do node = 1, numnod
           do level = 1, nrlevs
             state%surfacewater%qdra(level,node) = qdra(level,node)
           end do
         end do

      else
! --- drainage flux through lowest compartment
        do level = 1,nrlevs
           qdra(level,numnod) = qdrain(level)
           state%surfacewater%qdra(level,numnod) = qdrain(level)
        end do
      endif

      qdrtot = 0.0d0
      do level=1,nrlevs
          qdrtot = qdrtot + qdrain(level)
      end do
      state%surfacewater%qdrtot = qdrtot

      return

      case (3)

! === surface water balance ========================

      if (swsrf .eq. 3) then
        wlp = afgen (wlptab,2*mawlp,t1900-1.0d0+dt)
      endif
      if (swsec.eq.2) then
! ---    water level of secondary system is simulated
         call wlevbal (state, request_smaller_dt)
      elseif (swsec.eq.1) then
! ---    water level of secondary system is input
         call wballev (state)
      endif

      case default
         call fatalerr_collected ('SurfaceWater', 'Illegal value for TASK')
      end select

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
      use variables, only: tcum,NRPRI,impend,nmper,swman,wls,wlstar,hbweir,gwl,wlsman,gwlcrit,nphase,dropr,wscap,   &
                           dt,runots,QRapDra,qdrd,swst,zbotdr,alphaw,betaw,osswlm,T,NUMNOD,THETAS,THETA,DZ,VCRIT,NODHD,HCRIT, &
                           H,SWQHR,QQHTAB,wldip,intwl,t1900,logf,swscre,fldecdt,fldtmin,rsro,pond,pondmx,imper,sttab
      use swap_state_mod, only: swap_state_t
      use surfacewater_utils, only: wlevst, swstlev, qhtab
      IMPLICIT NONE

      type(swap_state_t), intent(inout) :: state
      logical,            intent(inout) :: request_smaller_dt

! --- local
      INTEGER iphase,NODE,Intday
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
         sw_sttab  => state%surfacewater%sttab,  &
         sw_cqdrd  => state%surfacewater%cqdrd,  &
         sw_cwsupp => state%surfacewater%cwsupp, &
         sw_cwout  => state%surfacewater%cwout)

! --- resetting of flag for overflowing of automatic weir
      ! overfl global write dropped: only sw_overfl (state alias) used henceforth.
      sw_overfl = .false.

! --- memorizing previous target level
      wlstarb = wlstar

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

      if (t1900-1.d0+0.1d-10 .gt. impend(imper)) goto 100

! --- determine the target sw-level:
      if (swman(imper) .eq. 1) then

! --- In the case of a fixed weir the 'target level' is set to the
! --- weir crest, for later use in calculations to determine whether
! --- there is any outflow at all (see below):
        wlstar = hbweir(imper)
        sw_wlstar = hbweir(imper)
      else

! ---  For automatic weir, determine the target level of the surface
! ---  water, dependent on groundwater level:

! ---   only adjust it if new subperiod, with length intwl(imper) has
!       started  (or if it is is the first call)
        rday = (T+1.0D0)/intwl(imper)
        intday = int(rday)

        if (abs(rday-1.0*intday).lt.0.00001d0 .or. tcum.lt.1.0d-10) then

          iphase = nphase(imper)
          do while(gwl.gt.gwlcrit(imper,iphase).and.iphase.gt.1)
            iphase = iphase-1
          enddo

! --- compare total air volume with VCRIT, adapt iphase
          ! VTAIR global write dropped; sw_vtair (state alias) used as accumulator.
          sw_vtair = 0.0d0
          do NODE = 1,NUMNOD
            sw_vtair = sw_vtair + (THETAS(NODE)-THETA(NODE))            &
     &              *abs(DZ(NODE))
          enddo
          do while (sw_vtair.lt.VCRIT(imper,iphase).AND.IPHASE.gt.1)
            iphase = iphase - 1
          enddo

! --- compare H(nodhd(imper)) with HCRIT, adapt iphase
          do while (h(nodhd(imper)).gt.hcrit(imper,iphase)              &
     &                            .and.iphase.gt.1)
            iphase = iphase - 1
          enddo
          ! hwlman global write dropped; sw_hwlman (state alias) set directly.
          sw_hwlman = h(nodhd(imper))

          wlstx = wlsman(imper,iphase)
        else

! --- use old level
          wlstx = wlstarb
        endif

! ---   if the level must drop, then do not let it drop at more than
! ---   the specified rate:
        if (wlstx .lt. wlstar .and. dropr(imper) .gt. 0.001d0) then
          wlstar = wlstar - dropr(imper)*dt
          if (wlstar .lt. wlstx) wlstar = wlstx
        else
          wlstar = wlstx
        endif
        sw_wlstar = wlstar
      endif

! --- counter of adjustments
      if (abs(wlstar-wlstarb) .gt. 0.00001d0) then
         ! numadj global write dropped; sw_numadj (state alias) accumulates directly.
         sw_numadj = sw_numadj + 1
      endif

! --- storage for the 'target level'
      swsttar = swstlev(wlstar)

! --- level and storage for "max. level for supply"
      wlstara = wlstar - wldip(imper)
      if (wlstara .gt. (zbotdr(1+nrpri)+1.d-4)) then
         swsttara = swstlev(wlstara)
         wsmax = wscap(imper)
      else
         swsttara = 0.0d0
         wsmax = 0.0d0
      endif

! --- determine whether the system will become full (target level or
! --- level of weir crest):
      dvmax = (qdrd + QRapDra + wsmax) * dt + runots
      swstmax = swst + dvmax

      if (swstmax .lt. 1.0d-7) then
! ---   storage decreases to zero, then the surface water system
! ---   falls dry; set surface water supply to maximum:

        if (swstmax .lt. -0.1d0) then
          messag = 'error algorithm for sw falling dry'
          call fatalerr_collected ('Wlevbal',messag)
        endif

        wsupp = wsmax
        wdis = 0.0d0
        swst = 0.0d0
        sw_swst = 0.0d0
        wls = zbotdr(nrpri+1)
        sw_wls = wls

      elseif (swstmax .ge. 0.0d0 .and. swstmax .lt. swsttara) then
! --- system will not become full - set supply to max. capacity
        wsupp = wsmax
        wdis = 0.0d0
        swst = swstmax
        sw_swst = swstmax

! --- calculate new level from storage
        wls = wlevst(swst)
        sw_wls = wls
      else

! --- determine how system will become full: with or without needing
!     surface water supply; try first without any supply:
        dvmax = (qdrd + QRapDra) * dt + runots
        swstmax = swst + dvmax
        if (swstmax .le. swsttara) then

! --- apparently supply is needed for reaching target level, system
!     is made full up to level wlstara, because supply is controllable:
          wsupp = (swsttara - swst - (qdrd+QRapDra)*dt-runots)/dt
          wdis = 0.0d0
          swst = swsttara
          sw_swst = swsttara
          wls = wlstara
          sw_wls = wlstara
        elseif (swstmax .le. swsttar) then

! --- apparently drainage is more than sufficient for filling system
!     to a level above the target level FOR SUPPLY, but not enough
!     for generating discharge; so calculate level from storage:
          wsupp = 0.0d0
          wdis = 0.0d0
          swst = swstmax
          sw_swst = swstmax
          wls = wlevst(swst)
          sw_wls = wls
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
            wdis = (swst-swsttar + (qdrd+QRapDra)*dt + runots )/dt

! --- now check whether the weir has enough discharge capacity
!     at this water level
            if (SWQHR.eq.1) then
              wover = wls - hbweir(imper)
              discap = alphaw(imper) * (wover**betaw(imper))
            elseif (SWQHR.eq.2) then

! --- interpolate QH table
              discap = qhtab(wlstar)
            endif
            if (discap .gt. wdis) then
              wls = wlstar
              sw_wls = wlstar
              swst = swsttar
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
              wover = sttab(1,1) - hbweir(imper)
              discap = alphaw(imper) * (wover**betaw(imper))
            elseif (SWQHR.eq.2) then
              discap = QQHTAB(imper,1)
            endif

! ---       error handling
            swstn = swst + (qdrd + QRapDra - discap)*dt + runots
            if ( swstn .gt. sttab(1,2) ) then
              messag = 'surface water system has overflowed!'
              call fatalerr_collected ('Wlevbal',messag)
            endif

! --- iteration procedure for determining new level, storage,
!     and discharge
            wlsl = hbweir(imper)
            wlsu = sttab(1,1)

! --- find storage and discharge for intermediate point
 700        wlsi = (wlsl + wlsu) * 0.5
            swsti = swstlev(wlsi)
            if (SWQHR.eq.1) then
              wdisi = alphaw(imper)*(wlsi-hbweir(imper))**betaw(imper)
            else
              wdisi = qhtab(wlsi)
            endif
            swstn = swst + (qdrd + QRapDra - wdisi)*dt + runots
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
              wls = wlsi
              sw_wls = wlsi
              swst = swstn
              sw_swst = swstn
              wdis = wdisi
            endif
          endif
        endif
      endif

!        ponding in case of extended drainage may limit timestep

      if (wls.gt.pondmx .or. pond.gt.pondmx) then
        if(dt .gt. 0.02*rsro) then
          request_smaller_dt = .true.
          fldecdt = .true.   ! transitional dual-write
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
           if (.not.fldtmin ) then
              request_smaller_dt = .true.
              fldecdt = .true.   ! transitional dual-write
              fl_early_return = .true.
           else
              call dtdpst                                               &
     &        ('year-month-day,hour:minute:seconds',t1900,datetime)
              messag = ' sw-level oscillation at '//datetime//          &
     &        '       advise: reduction of dtmax !'
              call warn ('Wlevbal',messag,logf,swscre)
              call fatalerr_collected ('Wlevbal',messag)
           endif
        endif
      endif
      endif  ! .not. fl_early_return

      if (.not. fl_early_return) then
! --- cumulative terms (global accumulations dropped; state aliases are authoritative):
      sw_cqdrd  = sw_cqdrd  + qdrd*dt
      sw_cwsupp = sw_cwsupp + wsupp*dt
      sw_cwout  = sw_cwout  + wdis*dt
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
      use variables, only: wls,wlstab,swst,dt,runots,QRapDra,qdrd,t1900
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
         sw_cwout  => state%surfacewater%cwout)

! --- wlsold gets w-level of previous time step
      ! wlsold global write dropped; sw_wlsold (state alias) is the signal;
      ! wls global kept because drainage.f90 (bocodre) reads it.
      sw_wlsold = wls

! --- fetch new level from input series
      wls = AFGEN (WLSTAB,2*MAWLS,t1900-1.d0+DT)
      sw_wls = wls

! --- determine surface water storage for level(t-dt) and level(t)
      swstold = swstlev(sw_wlsold)
      swst = swstlev(wls)
      sw_swst = swst

! --- determine from the surface water storages and the qdrain whether
! --- supply has taken place during period (t)-(t+dt) or water has
! --- been discharged
      swstrest = swstold + (qdrd + QRapDra)*dt + runots - swst

! --- if supply was needed, set discharge to zero
      if (swstrest.le.0.0d0) then
        wdis = 0.0d0
        wsupp = -swstrest/dt

! --- if discharge has taken place, set supply to zero
      else
        wdis = swstrest/dt
        wsupp = 0.0d0
      endif

! --- cumulation of water balance terms (global accumulations dropped; state aliases authoritative):
      sw_cqdrd  = sw_cqdrd  + qdrd*dt
      sw_cwsupp = sw_cwsupp + wsupp*dt
      sw_cwout  = sw_cwout  + wdis*dt

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
