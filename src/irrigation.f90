! File VersionID:
!   $Id: irrigation.f90 372 2018-03-13 10:01:20Z heine003 $
! ----------------------------------------------------------------------
      subroutine irrigation (task)
! ----------------------------------------------------------------------
!     Date               : November 2004
!     Purpose            : evaluate and schedule irrigations
!                        : Modified July 2019; 
!                             - tcs=5 obsolete and replaced by tcs=7 (theta) and tcs=8 (presh)
!                             - some calculations only once during initialization
! ----------------------------------------------------------------------
! --  global variables
      use variables
      implicit none

! --  local variables
      integer irr,node,nodsen,task,tcs,tcsfix,dcslim,dcs
      integer ifnd,i,datea(6),getun2,irgdayfix
      integer endirr(2),startirr(2)
      integer yearendcrp, yearstacrp
      real(8) frlow,phlo,phhi,phme,awlh,awmh,awah,cdef
      real(8) wclo,wcme,wchi,wcac,tps1,tps2,tps3,tps4,tps5,depl,phcrit
      real(8) watcon,dps1,dps2,afgen,Tred
      real(8) dvstage(7),trel(7),raw(7),taw(7),dwa(7),hcri(7),tcri(7)
      real(8) di(7),fid(7),irgdepmax,irgdepmin,irgthreshold
      real(4) fsec
      real(8) tstairryrx,tendirryrx, grai_red
      real(8), dimension(maho) :: wclos, wcmes, wchis
      logical flIrriTime
      logical rdinqr
      character(len=80) filnam
      character(len=200) messag

!     save local variables
      save    

!     dcs(1)  = Amount of under- or over-irrigation (L) in case of a scheduled irrigation event
!     dcs(2)  = Prescribed fixed irrigation depth (L) for each scheduled irrigation event

! ----------------------------------------------------------------------
      select case (task)
      case (1)

! === initialization calculated irrigation ============================
!     At start of new crop
         dayfix = 366

! ---    open crop file
         filnam = trim(pathcrop)//trim(cropfil(icrop))//'.crp'
         irr = getun2 (10,90,2)
         call rdinit(irr,logf,filnam)

! ---    determine whether irrigation scheduling is applied
         call rdsinr ('schedule',0,1,schedule)

         if (schedule .eq. 1) then

            flIrrigate = .true.          

! ---       no scheduled irrigation before this date
            call rdfinr('startirr',1,31,startirr,2,2)
            call rdfinr('endirr',1,31,endirr,2,2)
            datea(1) = 1900
            datea(2) = startirr(2)
            datea(3) = startirr(1)        
            fsec = 0.0
            call dtardp (datea, fsec, tstairrig)
            datea(2) = endirr(2)
            datea(3) = endirr(1)        
            call dtardp (datea, fsec, tendirrig)

            call dtdpar (cropstart(icrop),datea,fsec)
            yearstacrp = datea(1)
            call dtdpar (cropend(icrop),datea,fsec)
            yearendcrp = datea(1)
            if (yearendcrp.gt.yearstacrp) then
               datea(1) = yearstacrp
               datea(2) = startirr(2)
               datea(3) = startirr(1)        
               fsec     = 0.0
               call dtardp (datea, fsec, tstairryrx)
               datea(1) = yearendcrp
               datea(2) = endirr(2)
               datea(3) = endirr(1)        
               call dtardp (datea, fsec, tendirryrx)
            end if
            
! ---       application method
            call rdsinr ('isuas',0,1,isuas)
! ---       solute concentration irrigation water
            call rdsdor ('cirrs',0.0d0,100.0d0,cirrs)

! ---       timing criteria
            call rdsinr ('tcs',1,8,tcs)
            if (tcs == 5) then
               write (logf,'(A)') ' Option tcs=5 (together with phormc) is obsolete and cannot be used anymore.'
               write (logf,'(A)') ' Use option tcs=7 for pressure head (previously: tcs=5 + phormc=0) or '
               write (logf,'(A)') '            tcs=8 for water content (previously: tcs=5 + phormc=1).'
               write (*,'(A)') ' Option tcs=5 (together with phormc) is obsolete and cannot be used anymore.'
               write (*,'(A)') ' Use option tcs=7 for pressure head (previously: tcs=5 + phormc=0) or '
               write (*,'(A)') '            tcs=8 for water content (previously: tcs=5 + phormc=1).'
               call fatalerr ('Irrigation','Option tcs=5 (together with phormc) is obsolete and cannot be used anymore.')
            end if

! -1-       timing - allowable daily stress
!           minimum of ratio actual/potential transpiration
            if (tcs.eq.1) then
               call rdador('dvs_tc1',0.0d0,2.0d0,dvstage,7,ifnd)
               call rdfdor('trel',0.0d0,1.0d0,trel,7,ifnd)
               do i = 1, ifnd
                  treltab(i*2) = trel(i)
                  treltab(i*2-1) = dvstage(i)
               end do
            end if

! -2-       timing - depletion of readily available water
!           minimum fraction of readily available water
            if (tcs.eq.2) then
               call rdador('dvs_tc2',0.0d0,2.0d0,dvstage,7,ifnd)
               call rdfdor('raw',0.0d0,1.0d0,raw,7,ifnd)
               do i = 1, ifnd
                  rawtab(i*2) = raw(i)
                  rawtab(i*2-1) = dvstage(i)
               end do
            end if

! -3-       timing - depletion of totally available water
!           minimum fraction of totally available water
            if (tcs.eq.3) then
               call rdador('dvs_tc3',0.0d0,2.0d0,dvstage,7,ifnd)
               call rdfdor('taw',0.0d0,1.0d0,taw,7,ifnd)
               do i = 1, ifnd
                  tawtab(i*2) = taw(i)
                  tawtab(i*2-1) = dvstage(i)
               end do
            end if

! -4-       timing - allowable depletion amount            
!           maximum amount of water (L) depleted below field capacity
            if (tcs.eq.4) then
                call rdador('dvs_tc4',0.0d0,2.0d0,dvstage,7,ifnd)
                call rdfdor('dwa',0.0d0,500.0d0,dwa,7,ifnd)
                do i = 1, ifnd
                  dwatab(i*2) = dwa(i)
                  dwatab(i*2-1) = dvstage(i)
               end do
            end if

! -5-       timing - critical press. head or moist. content at sensor depth
!              OBSOLETE replaced by tcs=7 or tcs=8; see below  (july, 2019)

! -6-       Timing - fixed intervals, weekly with threshold (below it: no irrigation)
            if (tcs.eq.6) then
               call rdsdor ('irgthreshold',0.0d0,20.0d0,irgthreshold)
               dayfix = 366
            end if

! -7-       Timing - fixed intervals
            call rdsinr ('tcsfix',0,1,tcsfix)
            if (tcsfix.eq.1) then
               if (tcs.eq.6) then   
                  messag = 'Timing criteria: conflict with fixed intervals tcsfix=1 AND tcs=6 not allowed, adapt input !'
                  call fatalerr ('Irrigation',messag)
               end if
               call rdsinr('irgdayfix',1,366,irgdayfix)
            end if

! -7-       timing - critical pressure head at sensor depth (dcrit; read later)
            swcirrthres = 0
            if (tcs == 7) then
               call rdador ('dvs_tc7',0.0d0,2.0d0,dvstage,7,ifnd)
               call rdfdor ('value_tc7',-1.0d6,100.0d0,hcri,7,ifnd)
               do i = 1, ifnd
                  hcritab(i*2)   = hcri(i)
                  hcritab(i*2-1) = dvstage(i)
               end do
            end if
            
! -8-       timing - critical water content at sensor depth (dcrit; read later)
            if (tcs == 8) then
               call rdador ('dvs_tc8',0.0d0,2.0d0,dvstage,7,ifnd)
               call rdfdor ('value_tc8',0.0d0,1.0d0,tcri,7,ifnd)
               do i = 1, ifnd
                  tcritab(i*2)   = tcri(i)
                  tcritab(i*2-1) = dvstage(i)
               end do
            end if
            
            if (tcs == 7 .or. tcs == 8) then
! ---          sensor depth for options 7 or 8
               call rdsdor ('dcrit',-100.0d0,0.0d0,dcrit)
! ---          determine compartment number of sensor depth
               i = 1
               do while (zbotcp(i) .gt. (dcrit + 1.0d-5))
                  i = i+1
               end do
               nodsen = i

! -7b, 8b -    in case of solutes: allow overirrigation when conc exceeds concthreshold
               cirrthres = 0.0d0
               perirrsurp = 0.0d0
               if (rdinqr('swcirrthres')) then
                  call rdsinr ('swcirrthres',0,1,swcirrthres)
                  if (swcirrthres.eq.1) then
                     call rdsdor ('cirrthres', 0.0d0,100.0d0,cirrthres)
                     call rdsdor ('perirrsurp',0.0d0,100.0d0,perirrsurp)
                  end if
               end if

            end if

! ---       depth criteria
            call rdsinr ('dcs',1,2,dcs)

! ---       field capacity as input for timing and depth options
            if (tcs.eq.2 .or. tcs.eq.3 .or. tcs.eq.4 .or. dcs.eq.1) then
               call rdsdor('phFieldCapacity',-500.0d0,0.0d0,phlo)
            else
               phlo = -100.0d0
            end if

! -1-       depth - back to field capacity
            if (dcs.eq.1) then
               call rdador('dvs_dc1',0.0d0,2.0d0,dvstage,7,ifnd)
               call rdfdor('di',-100.0d0,+100.0d0,di,7,ifnd)
               raithreshold = 0.0d0
               if (rdinqr('raithreshold')) then
                  call rdsdor('raithreshold',0.0d0,1000.0d0,raithreshold)
               end if
               do i = 1, ifnd
                  ditab(i*2)   = di(i)
                  ditab(i*2-1) = dvstage(i)
               end do
            end if

! -2-       depth - fixed depth
            if (dcs.eq.2) then
               call rdador('dvs_dc2',0.0d0,2.0d0,dvstage,7,ifnd)
               call rdfdor('fid',0.0d0,400.0d0,fid,7,ifnd)
               do i = 1, ifnd
                  fidtab(i*2)   = fid(i)
                  fidtab(i*2-1) = dvstage(i)
               end do
            end if

! ---       depth - limited by min and max
            call rdsinr ('dcslim',0,1,dcslim)
            if (dcslim.eq.1) then
               call rdsdor('irgdepmin',0.0d0,100.0d0,irgdepmin)
               call rdsdor('irgdepmax',irgdepmin,1.0d7,irgdepmax)
            else
               irgdepmin = 0.0d0
               irgdepmax = 1.0d7
            end if

! ---       for some options some constants can be computed once (per crop) and stored
!            phlo = -100.0d0 ! field capacity became input (from swap32(6))
            ! read hlim values here, since reading crop file for crop initialization occurs one time step (or day) later (?)
            call rdsdor ('hlim3h',-10000.0d0,100.0d0,hlim3h)               
            call rdsdor ('hlim3l',-10000.0d0,100.0d0,hlim3l)               
            call rdsdor ('hlim4' ,-20000.0d0,100.0d0,hlim4)                
            phme = (hlim3l+hlim3h)*0.5d0
            phhi = hlim4
            do i = 1, numlay
               node     = botcom(i)    ! bottom most compartment per soil layer
               wclos(i) = watcon(node,phlo)
               wcmes(i) = watcon(node,phme)
               wchis(i) = watcon(node,phhi)
            end do

         end if ! if (schedule .eq. 1)

! -      close input file
         close(irr)
         return
      
      case (2)

! ===    determine irrigation rates and states  =============================================
!        daily
      
! ---    reset intermediate soil water fluxes
         if (flzerointr) then
            igird = 0.0d0
            inird = 0.0d0
         end if

! ---    reset cumulative soil water fluxes
         if (flzerocumu) then
            cgird = 0.0d0
            cnird = 0.0d0
         end if

         gird = 0.0d0
         irrigevent = 0

! ---    fixed irrigations events
         if (swirfix .eq. 1) then
            if (abs(irdate(nirri) - t1900) .lt. 1.d-3) then
               gird = irdepth(nirri)
               cirr = irconc(nirri)
               isua = irtype(nirri)
               nirri = nirri + 1
               irrigevent = 1
            end if
         end if

! ---    scheduling mode - current timing and depth criterion

!        scheduled timing within desired period ?
         if (schedule.eq.1) then
            !call dtdpar (cropstart(icrop),datea,fsec)
            !yearstacrp = datea(1)
            !call dtdpar (cropend(icrop),datea,fsec)
            !yearendcrp = datea(1)
            flIrriTime = .false.
            if (yearendcrp.gt.yearstacrp) then
               !datea(1) = yearstacrp
               !datea(2) = startirr(2)
               !datea(3) = startirr(1)        
               !fsec = 0.0
               !call dtardp (datea, fsec, tstairryrx)
               !datea(1) = yearendcrp
               !datea(2) = endirr(2)
               !datea(3) = endirr(1)        
               !call dtardp (datea, fsec, tendirryrx)
               if ( (t1900-tstairryrx).gt.1.0d-3 .and. (t1900-tendirryrx).le.1.0d-3 ) then
                  flIrriTime = .true.
               end if
            else
               if ((t-tstairrig).ge.-1.0d-3.and.(t-tendirrig).le.1.0d-3) then
                  flIrriTime = .true.
               end if
            end if
         end if

         if (schedule.eq.1 .and. irrigevent.eq.0 .and. flCropCalendar .and. .not. flCropHarvest .and. flIrriTime) then
            cirr = cirrs
            isua = isuas

! ---       determine water holding capacity, readily available water, 
! ---       actual available water and water deficit
            frlow = (ztopcp(noddrz) + rd) / dz(noddrz)
            awlh = 0.0d0; awmh = 0.0d0; awah = 0.0d0; cdef = 0.0d0
            do node = 1,noddrz
               wclo = wclos(layer(node))*dz(node);       if (node.eq.noddrz) wclo = wclo*frlow
               wcme = wcmes(layer(node))*dz(node);       if (node.eq.noddrz) wcme = wcme*frlow
               wchi = wchis(layer(node))*dz(node);       if (node.eq.noddrz) wchi = wchi*frlow
               wcac = watcon(node,h(node))*dz(node);     if (node.eq.noddrz) wcac = wcac*frlow
               awlh = awlh+(wclo-wchi)
               awmh = awmh+(wcme-wchi)
               awah = awah+(wcac-wchi)
               cdef = cdef+(wclo-wcac) 
            end do

! -1-       timing - allowable daily stress - only under dry stress circumstances
            if (tcs.eq.1) then
               tps1 = afgen(treltab,14,dvs)
! ---          transpiration fraction due to drought and salinity stress
               if (iptra_day .gt. 1.d-10) then
                  Tred = 1.0d0 - (iqreddry_day + iqredsol_day) / iptra_day
               else
                  Tred = 1.0d0
               end if
               if (Tred .lt. tps1) irrigevent = 2
            end if

! -2-       timing - depletion of readily available water (fraction)
            if (tcs.eq.2) then
! ---          compare readily available water and actual available water
               tps2 = afgen(rawtab,14,dvs)
               depl = tps2*(awlh-awmh)
               if (depl.gt.awlh) depl=awlh 
               if (awah .lt. (awlh-depl)) irrigevent = 2
            end if

! -3-       timing - depletion of totally available water (fraction)
            if (tcs.eq.3) then
! ---          compare totally available water and actual available water
               tps3 = afgen(tawtab,14,dvs)
               depl = tps3*awlh
               if (awah.lt.(awlh-depl)) irrigevent = 2
            end if

! -4-       timing - allowable amount of depletion
            if (tcs.eq.4) then
! ---          check if depletion amount has been exceeded                
               tps4 = afgen(dwatab,14,dvs)
               if ((awlh-awah).gt.(tps4*0.1d0)) irrigevent = 2
            end if

! -5-       timing - critical pressure head or moisture content exceeded
!             OBSOLETE: replaced by tcs=7 or tcs=8

! -6-       Timing - fixed irrigation time (weekly during crop growth)
            if (tcs.eq.6) then

!           (weekly) irrigation only when deficit is higher then threshold
!              cdef (cm) en IrgThreshold (mm)
               dayfix = dayfix + 1
               if (dayfix.ge.7)  then
                  dayfix = 0
                  if (10.0d0*cdef.gt.irgthreshold) then
                     irrigevent = 2
                  end if
               end if
            end if

! -7-       timing - critical pressure head at dcrit (node=nodsen) exceeded
            if (tcs.eq.7) then
! ---          calculation of critical pressure head
               tps5 = afgen(hcritab,14,dvs)
! PG/JK start  15-feb-2010
! originally not intended to simulate paddy rice fields,
! but made applicable for paddy by changing the statement:
               phcrit = tps5        ! old statement was: phcrit = -abs(tps5)
! PG/JK end    15-feb-2010
! ---          compare critical pressure head and actual pressure head
               if (h(nodsen).le.phcrit) irrigevent = 2
            end if

! -8-       timing - critical watercontent at dcrit (node=nodsen) exceeded
            if (tcs.eq.8) then
               tps5 = afgen(tcritab,14,dvs)
! ---          compare critical water content and actual water content
               if (theta(nodsen).le.tps5) irrigevent = 2
               !phcrit = prhead(nodsen,disnod(nodsen),tps5,cofgen,h)
            end if

! -9-       Timing - fixed intervals
            if (tcsfix.eq.1) then
               if (irrigevent.eq.2 .and. (dayfix .ge. irgdayfix)) then
                  irrigevent = 2
                  dayfix = 1
               else
                  irrigevent = 0
                  if (dayfix .lt. irgdayfix) dayfix = dayfix + 1
               end if
            end if

! ---       depth - back to field capacity [cm]
            if ((irrigevent.eq.2).and.(dcs.eq.1)) then
! ---       correct for over- or under irrigation
               dps1 = afgen(ditab,14,dvs)
! PG/JK start  15-feb-2010
! option to reduce irrigation on rainy (> raithreshold) day
! raithreshold =     ! threshold (cm/d) to define rainy days;  used to reduce irrigation
               grai_red = 0.0d0
               if (grai .gt. raithreshold) grai_red = grai
               gird = max (0.0d0,cdef+dps1*0.1d0-grai_red) 
! PG/JK start  15-feb-2010
            end if

! ---       depth - fixed depth [cm]
            if ((irrigevent.eq.2).and.(dcs.eq.2)) then
               dps2 = afgen(fidtab,14,dvs)
               gird = dps2*0.1d0
            end if

! ---       depth - limited depth [cm]
            if ((irrigevent.eq.2).and.(dcslim.eq.1)) then
               gird = max(gird,irgdepmin*0.1d0)
               gird = min(gird,irgdepmax*0.1d0)
            end if

! ---       in case of solutes: allow overirrigation when conc exceeds concthreshold
            if (swsolu.eq.1 .and.irrigevent.eq.2 .and.swcirrthres.eq.1) then
               if (cml(nodsen).gt.cirrthres) then
                  gird = gird + 0.01d0*perirrsurp*gird
               end if
            end if

         end if

         if (irrigevent .ne. 0) flIrrigationOutput = .true.

      case default
         call fatalerr ('Irrigation', 'Illegal value for TASK')
      end select

      return
   end

subroutine SSDI_irrigation(iTask)

use variables, only: swpfile, logf, mairg, numnod, tend, tstart, t1900, zbotcp, irrigevent, qssdi, qssdisum, dt_SSDI_event,   &
                     h, theta, iptra_day, iqreddry_day, iqredsol_day

implicit none
! global
integer, intent(in) :: iTask

! local, saved
integer,                   save :: swssdi             ! switch indicating if ssdi is active in current simulation
integer, dimension(2),     save :: nod_ssdi           ! upper and lower nodes where ssdi takes place
integer,                   save :: ssdi_schedule      ! switch indicating type of scedule (0 = fixed dates; 1 = following internal schedule)
! specific is ssdi_schedule = 1
integer,                   save :: ssdi_sched_type    ! type of internal schedule: 1 = Tact/Tpot; 2 = pressure head; 3 = water content
integer,                   save :: nod_ssdi_sensor    ! node where sensor is located (if ssdi_sched_type > 1)
real(8),                   save :: ssdi_threshold     ! Threshold value for ssdi_sched_type
real(8),                   save :: ssdi_threshold_z   ! Depth for threshold value
real(8),                   save :: ssdi_amount        ! Amount of scheduled irirgation (input as mm, converted to cm)
real(8),                   save :: ssdi_appl_rate     ! Application rate (inputs as mm/h, converted to cm/d)
integer,                   save :: sw_interval        ! Switch indicating whether (1) or not (0)there must be a minimum number of days between succesive SSDI applications
integer,                   save :: days_interval      ! Miminmum number of days between succesive SSDI applications (d)
integer,                   save :: days_counter       ! Number of days since previsou SSDI application (d)

! specific is ssdi_schedule = 0
integer,                   save :: nirri              ! ssdi counter; entry point in array of ssdi dates
real(8), dimension(mairg), save :: ssdi_date          ! Array with fixed irrigation dates
real(8), dimension(mairg), save :: ssdi_rate_f        ! Array with fixed irrigation rates (cm/d; input as mm/d)
real(8), dimension(mairg), save :: ssdi_amount_f      ! Array with fixed irrigation amounts (cm; input as mm)

! local, help
integer                         :: i, j, swp, ifnd
real(8)                         :: Tred
real(8), dimension(2)           :: ssdi_z
character(len=132)              :: ssdi_file

! functions
integer :: getun2
logical :: rdinqr

   select case  (iTask)
   case (1)
      swp = getun2(10, 90, 2)
      call rdinit(swp, logf, swpfile)
         swssdi = 0
         if (rdinqr('swssdi')) call rdsinr ('swssdi', 0, 1, swssdi)
         if (swssdi == 1)      call rdscha ('ssdi_file', ssdi_file)
      close(swp)
      
      if (swssdi == 0) then
         nod_ssdi = 0
         qssdi    = 0.0
         nirri    = 1
         dt_SSDI_event = 1.0d0
         return
      else
         call read_ssdi_input()
      end if
      
      ! check if ssdi_date is ascending, and determine initial entry point nirri
      nirri = 0
      if (ssdi_schedule == 0) then
         nirri = 1
         do i = 1, ifnd-1
            if (ssdi_date(i) >= ssdi_date(i+1)) call fatalerr ('SSDI_irrigation', 'ssdi_date not in ascending order')
            if (t1900 >= ssdi_date(i)) nirri = i
         end do
         if (t1900 >= ssdi_date(ifnd)) nirri = ifnd
      end if

      ! determine layer number of sensor (if applicable)
      nod_ssdi_sensor = 0
      if (ssdi_schedule == 1 .AND. ssdi_sched_type > 1) then
         i = 1
         do while (zbotcp(i) .gt. (ssdi_threshold_z + 1.0d-5))
            i = i + 1
         end do
         nod_ssdi_sensor = i
      end if
      
      ! determine layer number where SSDI takes place
      do j = 1, 2
         i = 1
         do while (zbotcp(i) .gt. (ssdi_z(j) + 1.0d-5))
            i = i + 1
         end do
         nod_ssdi(j) = i
      end do

      ! redefine application rate: uniformly spread over all nodes
      if (ssdi_schedule == 0) then
         ssdi_amount_f = ssdi_amount_f/dble((nod_ssdi(2) - nod_ssdi(1) + 1))
      else
         ssdi_amount = ssdi_amount/dble((nod_ssdi(2) - nod_ssdi(1) + 1))
      end if
      
      ! initialize
      qssdi         = 0.0d0
      dt_SSDI_event = 1.0d0
      
   case (2)
      ! no subsurface drip irrigation: return
      if (swssdi == 0) return
      
      irrigevent      = 0
      qssdi(1:numnod) = 0.0d0
      dt_SSDI_event   = 1.0d0
      qssdisum        = 0.0d0
      
      if (ssdi_schedule == 0) then
         ! check if today is a day with ssdi
         if (abs(ssdi_date(nirri) - t1900) .lt. 1.d-3) then
            irrigevent                     = 2
            dt_SSDI_event                  = ssdi_amount_f(nirri) / ssdi_rate_f(nirri)
            qssdi(nod_ssdi(1):nod_ssdi(2)) = ssdi_rate_f(nirri)
            nirri                          = nirri + 1
            qssdisum                       = qssdisum + sum(qssdi(nod_ssdi(1):nod_ssdi(2)))
         end if
      else
         ! scheduling based on exceedance of a certain threshold
         if (ssdi_sched_type == 1) then
            ! transpiration fraction due to drought and salinity stress
            if (iptra_day .gt. 1.d-10) then
               Tred = 1.0d0 - (iqreddry_day + iqredsol_day) / iptra_day
            else
               Tred = 1.0d0
            end if
            if (Tred .lt. ssdi_threshold) irrigevent = 2
            
         else if (ssdi_sched_type == 2) then
            if (h(nod_ssdi_sensor) <= ssdi_threshold) irrigevent = 2
            
         else if (ssdi_sched_type == 3) then
            if (theta(nod_ssdi_sensor) <= ssdi_threshold) irrigevent = 2

         end if
         
         if (sw_interval == 1) then
            if (irrigevent == 2 .and. (days_counter >= days_interval)) then
               irrigevent   = 2
               days_counter = 1
            else
               irrigevent   = 0
               if (days_counter < days_interval) days_counter = days_counter + 1
            end if
         end if

         if (irrigevent == 2) then
            dt_SSDI_event                  = ssdi_amount / ssdi_appl_rate
            qssdi(nod_ssdi(1):nod_ssdi(2)) = ssdi_appl_rate
            qssdisum                       = qssdisum + sum(qssdi(nod_ssdi(1):nod_ssdi(2)))
         end if
         
      end if

   case (9)
      ! special: reset scheduled irrigation at end of irrigation event
      irrigevent      = 0
      qssdi(1:numnod) = 0.0d0
      dt_SSDI_event   = 1.0d0
      qssdisum        = 0.0d0
      
   case default
      call fatalerr ('SSDI_irrigation', 'Illegal value for iTask')
   end select
   
!---------------------------
   contains

   subroutine read_ssdi_input()
   real(8) :: dummy
   logical :: rdinar
   
   call rdinit(swp, logf, ssdi_file)
      call rdsinr('ssdi_schedule', 0, 1, ssdi_schedule)
      if (ssdi_schedule == 0) then
         call rdatim('ssdi_date',                      ssdi_date,     mairg, ifnd)
         call rdfdor('ssdi_rate_f',    0.0d0, 100.0d0, ssdi_rate_f,   mairg, ifnd)     ! mm/h
         call rdfdor('ssdi_amount_f',  0.0d0, 100.0d0, ssdi_amount_f, mairg, ifnd)     ! mm
         if (.NOT. rdinar('ssdi_z')) then
            ! application at single depth (single compartment)
            call rdsdor('ssdi_z', -100.0d0, 0.0d0, dummy)                ! cm
            ssdi_z(1) = dummy
            ssdi_z(2) = dummy
         else
            ! application over depth interval btween two depth levels (multiple consecutive compartments)
            call rdfdor('ssdi_z', -100.0d0, 0.0d0, ssdi_z, 2, 2)          ! cm
         end if
         ! convert mm/h or mm irrigation to cm/d or cm
         ssdi_rate_f(1:ifnd)   = ssdi_rate_f(1:ifnd)*0.1d0*24.0d0
         ssdi_amount_f(1:ifnd) = ssdi_amount_f(1:ifnd)*0.1d0
      else
         call rdsinr('ssdi_sched_type', 1, 3, ssdi_sched_type)
         if (ssdi_sched_type == 1) call rdsdor('threshold_Tred',     0.0d0,   1.0d0, ssdi_threshold)
         if (ssdi_sched_type == 2) call rdsdor('threshold_presh',   -1.0d7,   0.0d0, ssdi_threshold)
         if (ssdi_sched_type == 3) call rdsdor('threshold_watc',     0.0d0,   1.0d0, ssdi_threshold)
         if (ssdi_sched_type >  1) call rdsdor('threshold_depth', -100.0d0,   0.0d0, ssdi_threshold_z)      ! cm
                                   call rdsdor('ssdi_amount',        0.0d0, 100.0d0, ssdi_amount)           ! mm
                                   call rdsdor('ssdi_appl_rate',     0.0d0, 100.0d0, ssdi_appl_rate)        ! mm/h  !!!
         if (.NOT. rdinar('ssdi_z')) then
            ! application at single depth (single compartment)
            call rdsdor('ssdi_z', -100.0d0, 0.0d0, dummy)                ! cm
            ssdi_z(1) = dummy
            ssdi_z(2) = dummy
         else
            ! application over depth interval btween two depth levels (multiple consecutive compartments)
            call rdfdor('ssdi_z', -100.0d0, 0.0d0, ssdi_z, 2, 2)          ! cm
         end if
         ! Restriction on number of days between two successive SSDI appliocations
         call rdsinr('sw_interval', 0, 1, sw_interval)
         if (sw_interval == 0) then
            days_interval = 1
         else
            call rdsinr('days_interval', 1, 366, days_interval)
         end if
         days_counter = 366

         ! convert mm irrigation to cm
         ssdi_amount = ssdi_amount*0.1d0
         
         ! convert mm/h irrigation to cm/d
         ssdi_appl_rate = ssdi_appl_rate*0.1d0*24.0d0

      end if
   close(swp)
         
!  at least one date must be within simulation period
   if (ssdi_schedule == 0) call checkdate(ifnd, ssdi_date, tend, tstart, 'irdate', 'SSDI_irrigation//swssdi=1')
   
   end subroutine read_ssdi_input
   
end subroutine SSDI_irrigation
