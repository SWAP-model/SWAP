! File VersionID:
!   $Id: readswap.f90 379 2018-05-16 07:02:44Z heine003 $
! ----------------------------------------------------------------------
      subroutine readswap
! ----------------------------------------------------------------------
!     Date               : April 2014   
!     Purpose            : read main input file .SWP
! ----------------------------------------------------------------------
      use variables
      use doln
      implicit none

      integer posarg,numchar,mxcrop,idum
      integer swp,i,datea(6),getun,getun2,runf,swrunon
      integer swirgfil,datefix(2),ifnd,swyrvar,swmonth
      integer swerror,swbbcfile
      integer lay,swdc,ini,tss
      integer isublay(macp),bbc,il, sol, j
      integer SwMacro
      integer swvapor(maho)
      integer swIrrigate

      real(4) fsec
      real(8) outdate1,rainflux(30),raintime(30)
      real(8) ores(maho),osat(maho),alfa(maho),npar(maho)
      real(8) alfa_2(maho),npar_2(maho),omega_1(maho)       ! for bi-modal MvG
      real(8) h0(maho), ha(maho), apar(maho), omega_K(maho) ! for PDI
      real(8) maxdepth,lexp(maho),alfaw(maho)
      real(8) dates(mabbc),gwlevel(mabbc),cseeparr(mabbc)
      real(8) haquif(mabbc),hbottom(mabbc),ttop(mabbc)
      real(8) hhtab(mabbc),qhtab(mabbc),qboti(mabbc),tbot(mabbc)
      real(8) hthr,help,term1
      real(8) headtab(matab),thetatab(matab),conductab(matab)
      real(8) dydx(matab),pondmxtb(mairg),datepmx(mairg),sigma(matab)
      real(8) dum, dummy
      character(len=200) filenamesophy(maho)
      character(len=200) filnam,rufil
      character(len=80)  filtext
      character(len=400) swpfilnam,logfilnam
      character(len=32)  bbcfil,irgfil,swpfil,tsoilfile
      character(len=11)  tmp
      character(len=400) messag

      logical   toscr, flsat, rdinqr,rdinar

      real(8) nihil,vsmall
      parameter (vsmall = 1.d-16)
      parameter (nihil = 1.d-24)

      character(len=256)                  :: SWAP_SOILS_DB     ! name of database file with soil physical parameters
      character(len=40), dimension(maho)  :: SOILNAMES         ! names of soil from soils database to be used
      integer                             :: idb               ! for unit number to open SWAP_SOILS_DB
      character(len=1)                    :: cdum              ! dummy

!     pressure head where interpolation is used to calculate K from VG and ksatexm
      data hthr  /-2.0d0/
! ----------------------------------------------------------------------

! --- delete existing file Swap.ok
      call delfil ('Swap.ok',.false.)

! --- write message running to screen
      write (*,'(/,a)') '  running swap ....'

! --- path and filename of executable through argument command line
      PosArg = 1
      Call Get_Command_Argument (PosArg,swpfil,NumChar)
      if (NumChar.lt.1) swpfil = 'swap'
      project = swpfil

      if(NumChar.gt.3) then
         if(swpfil(NumChar-3:NumChar).eq.'.swp' .or.                    &
     &                       swpfil(NumChar-3:NumChar).eq.'.SWP') then
            swpfilnam = trim(swpfil)
         else
            swpfilnam = trim(swpfil)//'.swp'
         end if
      else 
        swpfilnam = trim(swpfil)//'.swp'
      end if
      Numchar = len(trim(swpfilnam))
      logfilnam = swpfilnam(1:NumChar-4)//'_swap.log'
      swpfile = trim(swpfilnam)

! --- open log file
      logf = getun (20,99)
      call fopens(logf,logfilnam,'new','del')

! --  write head to logfile
      filtext = 'Main input variables'
      call writehead(logf,1,swpfilnam,filtext,project)
      write (logf,*)
      
! --- open swp file
      swp = getun2 (10,90,2)
      call rdinit(swp,logf,swpfilnam)

! -   environment
      call rdscha ('project',project)
      call rdscha ('pathwork',pathwork)
      call rdscha ('pathatm',pathatm)
      call rdscha ('pathcrop',pathcrop)
      call rdscha ('pathdrain',pathdrain)
      call rdsinr ('swscre',0,2,swscre)
      call rdsinr ('swerror',0,1,swerror)
      if (swerror .eq. 1) then
        toscr = .true.
      else
        toscr = .false.
      endif
      call messini(toscr,.true.,logf)

! --  Shared simulation
      flSwapShared = .false.
      if(rdinqr('flSwapShared')) then
        call rdslog ('flSwapShared',flSwapShared)
        if(flSwapShared) then
          messag = 'Simulation with shared data exchange (flSwapShared)'
          call warn ('Readswap',messag,logf,swscre)
        endif
      endif

! -   simulation period
      call rdstim('tstart',tstart)
      call dtdpar (tstart+0.1d0, datea, fsec)
      iyear = datea(1)
      imonth = datea(2)
      call dtdpar (tstart, datea, fsec)
      call rdstim('tend',tend)

! -   check begin and end date of simulation
      if ((tend - tstart) .lt. 0.0d0) then
        messag = 'The end date of the simulation'//                     &
     &                    ' should be larger than the begin date!'
        call fatalerr ('readswap',messag)
      endif

! -   Output dates for balances
      call rdsinr ('swyrvar',0,1,swyrvar)
      if (swyrvar .eq. 0) then 
        call rdfinr('datefix',1,31,datefix,2,2)
        datea(1) = iyear
        datea(2) = datefix(2)
        datea(3) = datefix(1)
        fsec = 0.0
        call dtardp (datea, fsec, outdate1)
        if (outdate1 .lt. tstart) then
          datea(1) = datea(1) + 1
          call dtardp (datea, fsec, outdate1)
        endif
        i = 1
        outdat(i) = outdate1
        datea(1) = datea(1) + 1
        call dtardp (datea, fsec, outdate1)
        do while (outdate1 .lt. (tend + 0.1d0))
          i = i + 1
          outdat(i) = outdate1
          datea(1) = datea(1) + 1
          call dtardp (datea, fsec, outdate1)
        end do
      else
        call rdatim ('outdat',outdat,maout,ifnd)
      endif

! -   Intermediate output dates
      call rdsinr ('nprintday',1,1440,nprintday)

! --- output each time interval dt
      flprintdt = .false.
      if(rdinqr('flprintdt')) then
        call rdslog ('flprintdt',flprintdt)
      endif

      call rdsinr ('swmonth',0,1,swmonth)
      if (swmonth .eq. 1) then 
        datea(1) = iyear
        datea(2) = imonth
        if (datea(2) .lt. 12) then
          datea(2) = datea(2) + 1          
        else
          datea(1) = datea(1) + 1          
          datea(2) = 1          
        endif
        datea(3) = 1
        fsec = 0.0
        call dtardp (datea, fsec, outdate1)
        i = 0
        do while ((outdate1 - 1.d0) .lt. (tend + 0.1d0))
          i = i + 1
          outdatint(i) = outdate1 - 1.d0
          if (datea(2) .lt. 12) then
            datea(2) = datea(2) + 1          
          else
            datea(1) = datea(1) + 1          
            datea(2) = 1          
          endif
          call dtardp (datea, fsec, outdate1)
        enddo
        period = 0
        swres  = 0
        swodat = 0
      else
        call rdsinr ('period',0,366,period)
        call rdsinr ('swres',0,1,swres)
        
        call rdsinr ('swodat',0,1,swodat)
        if (swodat.eq.1) then  
          call rdatim ('outdatint',outdatint,maout,ifnd)
        endif  
      endif

! -   output files
      call rdscha ('outfil',outfil)
      call rdsinr ('swheader',0,1,swheader)

      call rdsinr ('swafo',0,3,swafo)
      call rdsinr ('swaun',0,2,swaun)
      if (swaun.ge.1 .or. swafo.ge.1) then
        call rdsinr ('swdiscrvert',0,1,swdiscrvert)
        if (swdiscrvert.eq.1) then

! -       read values for numnodNew and dzNew
          call rdsinr ('numnodNew',1,macp,numnodNew)
          call rdfdor ('dzNew',1.0d-6,5.0d2,dzNew,macp,numnodNew)
        endif
        call rdsdor ('CritDevMasBal',                                   &
     &                1.0d-30, 1.0d0,CritDevMasBal)
      endif
      
      swcsv = 0
      if(rdinqr('swcsv')) then
        call rdsinr ('swcsv',0,1,swcsv)
        if (swcsv == 1) call RDscha ('InList_csv', InList_csv)
      endif
      swcsv_tz = 0
      if(rdinqr('swcsv_tz')) then
        call rdsinr ('swcsv_tz',0,1,swcsv_tz)
        if (swcsv_tz == 1) call RDscha ('InList_csv_tz', InList_csv_tz)
        tz_z1_z2(1:2) = 99.9d0
        if (RDinqr('tz_z1_z2')) call RDfdor('tz_z1_z2', -1.0d4, 0.0d0, tz_z1_z2, 2, 2)
      endif

      swbal = 1
      if(rdinqr('swbal')) then
        call rdsinr ('swbal',0,1,swbal)
      endif

      swwba = 1
      if(rdinqr('swwba')) then
        call rdsinr ('swwba',0,1,swwba)
      endif

      swsba = 1
      if(rdinqr('swsba')) then
        call rdsinr ('swsba',0,1,swsba)
      endif
      
      swini = 1
      if(rdinqr('swini')) then
        call rdsinr ('swini',0,1,swini)
      endif

      swend = 1
      if(rdinqr('swend')) then
        call rdsinr ('swend',0,2,swend)
      endif
      
      swinc = 1
      if(rdinqr('swinc')) then
        call rdsinr ('swinc',0,1,swinc)
      endif
      swstr = 1
      if(rdinqr('swstr')) then
        call rdsinr ('swstr',0,1,swstr)
      endif
      swcrp = 1
      if(rdinqr('swcrp')) then
        call rdsinr ('swcrp',0,1,swcrp)
      endif
      swrum = 0
      if(rdinqr('swrum')) then
        call rdsinr ('swrum',0,1,swrum)
      endif
      swirg = 1
      if(rdinqr('swirg')) then
        call rdsinr ('swirg',0,1,swirg)
      endif

      call rdsinr ('swvap',0,1,swvap)
      call rdsinr ('swate',0,1,swtem)
      call rdsinr ('swblc',0,1,swblc)
      call rdsinr ('swbma',0,1,swbma)

! --- output of recharge/storage info for ModFlow
      swoutputmodflow = 0
      if (rdinqr('swoutputmodflow')) then
        call rdsinr ('swoutputmodflow',0,1,swoutputmodflow)
        if(swoutputmodflow.eq.1) then
          messag = 'simulation with addtional output for ModFlow'
          call warn ('Readswap',messag,logf,swscre)
        endif
      endif

! -   meteo
      call rdscha ('metfil',metfil)
      call rdsdor ('lat',-90.0d0,90.0d0,lat)
      


!     type of weather data
      call rdsinr ('swetr',0,1,swetr)
      
      if (swetr.eq.0) then
        
        call rdsinr ('swmetdetail',0,1,swmetdetail)
        
        call rdsdor ('alt',-400.0d0,3000.0d0,alt)
        call rdsdor ('altw',0.0d0,99.0d0,altw)
        
        if(rdinqr('swdivide')) then
          call rdsinr ('swdivide',0,1,swdivide)
        endif
        
        angstroma = 0.25d0
        if(rdinqr('angstroma')) then
          call rdsdor ('angstroma',0.0d0,1.0d0,angstroma)
        endif
        
        angstromb = 0.50d0
        if(rdinqr('angstromb')) then
          call rdsdor ('angstromb',0.0d0,1.0d0,angstromb)
        endif
        
      else

        if (rdinqr('swmetdetail')) then
          call rdsinr ('swmetdetail',0,1,swmetdetail)
          if(swmetdetail.eq.1) then
            swmetdetail = 0
            messag = 'SWMETDETAIL is set to 0, SWMETDETAIL = 1 only available in case of SWETR = 1'
            call warn ('Readswap',messag,logf,swscre)
          endif
        endif  
        
        if (rdinqr('swdivide')) then
          call rdsinr ('swdivide',0,1,swdivide)
          if(swdivide.eq.1) then
            swdivide = 0
            messag = 'SWDIVIDE is set to 0, SWDIVIDE = 1 only available in case of SWETR = 0'
            call warn ('Readswap',messag,logf,swscre)
          endif
        endif  
        
      endif   
      

!     detailed meteo input as option
      if (swmetdetail.eq.0) then
        
        call rdsinr ('swrain',0,3,swrain)
        call rdsinr ('swetsine',0,1,swetsine)

      else  
      
        call rdsinr ('nmetdetail',1,96,nmetdetail)
        
        if (rdinqr('swrain')) then
          call rdsinr ('swrain',0,3,swrain)
          if(swrain.ge.1) then
            swrain = 0
            messag = 'SWRAIN is set to 0, SWRAIN = 1, 2 or 3 only available in case of SWMETDETAIL = 0'
            call warn ('Readswap',messag,logf,swscre)
          endif
        endif  

        if (rdinqr('swetsine')) then
          call rdsinr ('swetsine',0,1,swetsine)
          if(swetsine.eq.1) then
            swetsine = 0
            messag = 'SWETSINE is set to 0, SWETSINE = 1 only available in case of SWMETDETAIL = 0'
            call warn ('Readswap',messag,logf,swscre)
          endif
        endif  

      endif
      
!     data for actual rainfall intensity
      if (swrain .eq. 1) then
        call rdador ('time',0.d0,366.d0,raintime,30,ifnd)
        call rdfdor ('rainflux',0.d0,1000.d0,rainflux,30,ifnd)
        do i = 1, ifnd
          raintab(i*2) = rainflux(i)
          raintab(i*2-1) = raintime(i)
        end do
      endif

      if (swrain .eq. 3) then
        call rdscha ('rainfil',rainfil)
      endif

! -   special case: if METFIL is provided wit hextension .MET, then all weather dta will be erad at once
      ! only possible if: SWMETDETAIL = 0 and SWRAIN = 0 or 2
      call lowerc (metfil)
      swMetFilAll = 0
      if (index(trim(metfil),".met") > 0) then
         swMetFilAll = 1
         if (swmetdetail == 1 .OR. swrain == 1 .OR. swrain == 3) then
            messag = 'Extension .met in metfil discarded because not allowed in combination with SWETSINE = 1 or SWRAIN = 1 or 3'
            call warn ('Readswap',messag,logf,swscre)
            idum = index(metfil,".met")
            metfil = trim(metfil(1:idum-1))
            swMetFilAll = 0
         end if
      end if

! -   crop rotation scheme

! -   switch crop simulation
      swcrop = 1   ! default: with crop simulation
      if(rdinqr('swcrop')) then
        call rdsinr ('swcrop',0,1,swcrop)
      endif

      if (swcrop .eq. 1) then

!       crop calendar
        mxcrop = 0
        
        call rdatim ('cropstart',cropstart,macrop,mxcrop)
        call rdftim ('cropend',cropend,macrop,mxcrop)
        call rdfcha ('cropfil',cropfil,macrop,mxcrop)
        call rdfinr ('croptype',1,3,croptype,macrop,mxcrop)
        
        do i = 1,mxcrop
          write(tmp,'(i11)') i
          tmp = adjustl(tmp)
          if ((cropstart(i+1)-cropend(i)).lt.0.5d0                      &
     &                                           .and. i.lt.mxcrop) then
            write(tmp,'(i11)') i+1
            tmp = adjustl(tmp)
            messag = 'The begin date of crop number '//trim(tmp)//      &
     &      ' should be larger than the end date of the former crop!'
            call fatalerr ('readswap',messag)
          endif
! -       check combination detailed crop with LAT>60.
          if(croptype(i).ne.1 .and. abs(lat).gt.66.5d0) then
            messag = 'Fatal combination for crop number '//trim(tmp)//  &
     &      'detailed crop module within polar circle (LAT>66.5)!'
            call fatalerr ('readswap',messag)
          endif
        end do

! ---   Rooting depth limitation (by soil conditions)
        call rdsdor ('rds',1.d0,5000.0d0,rdmax)

! ---   Nitrogen in Crop
! -     flag for nutrients in crop and soil
        flCropNut = .false.
        if(rdinqr('flCropNut')) then
          call rdslog ('flCropNut',flCropNut)
        endif

        flCropReadFile = .true.
        flCropOpenFile = .true.      
      endif
!

! --- fixed irrigation events
      swirgfil = 0
      irgfil = ' ' 
      call rdsinr ('swirfix',0,1,swirfix)
      if (swirfix .eq. 1) then
        do i = 1,mairg
          irdate(i) = 0.d0
        end do
        call rdsinr ('swirgfil',0,1,swirgfil)
        if (swirgfil .eq. 0) then
          call rdatim ('irdate',irdate,mairg,ifnd)
          call rdfdor ('irdepth',0.d0,1000.d0,irdepth,mairg,ifnd)
          call rdfdor ('irconc',0.d0,1000.d0,irconc,mairg,ifnd)
          call rdfinr ('irtype',0,1,irtype,mairg,ifnd)
! -     at least one date must be within simulation period
          call checkdate(ifnd,irdate,tend,tstart,'irdate',              &
     &                   'readswap/swirfix=1      ')
! ---   convert mm irrigation to cm
          do i = 1, ifnd
            irdepth(i) = irdepth(i) / 10.d0
          end do
        else if (swirgfil .eq. 1) then
          call rdscha ('irgfil',irgfil)
        endif
      endif

! =================================================================
!     Section soil profile

! --- initial water conditions
      call rdsinr ('swinco',1,3,swinco)
      if (swinco.eq.1) then
        call rdador ('zi',-1.d5,0.d0,zi,macp,ifnd)
        call rdfdor ('h',-1.d10,1.d4,h,macp,ifnd)
        nhead = ifnd
      elseif (swinco.eq.2) then
        call rdsdor ('gwli',-10000.0d0,1000.0d0,gwli)
      elseif (swinco.eq.3) then
!       read data from file with results from previous simulation
        call rdscha ('inifil',inifil)
      endif

! --  drainage resistance of surface runoff
      call rdsdor ('rsro',0.001d0,1.0d0,rsro)
      call rdsdor ('rsroexp',0.01d0,10.0d0,rsroexp)
      swuseCN = 0
      if (rdinqr('swuseCN')) then
         call rdsinr('swuseCN',0,1,swuseCN)
         if (swuseCN == 1 .AND. (swmetdetail == 1 .OR. swrain > 0)) then
            messag = 'You cannot use CN method in combination with swmetdetail == 1 .OR. swrain > 0'
            call fatalerr ('readswap',messag)
         end if
         !!!call rdsdor('CNref',1.0d-2,100.0d0,CNref)
         call RDatim ('CNtimTAB', CNtimTAB, mayrs*5, iCNtab)
         call RDfdor ('CNrefTAB', 0.0d0, 100.0d0, CNrefTAB, mayrs*5, iCNtab)
         wc_cor = 0  ! default
         if (RDinqr('wc_cor')) call rdsinr('wc_cor',0,2,wc_cor)
      end if

! --- ponding
      swpondmx = 0
      if(rdinqr('swpondmx')) then
        call rdsinr ('swpondmx',0,1,swpondmx)
        if(swpondmx.eq.1) then
          messag = 'Simulation with time dependent ponding-threshold'
          call warn ('Readswap',messag,logf,swscre)
          call rdatim ('datepmx',datepmx,mairg,ifnd)
! -     at least one date must be within simulation period
          call checkdate(ifnd,datepmx,tend,tstart,'datepmx ',           &
     &                      'readswap/swpondmx=1')
          call rdfdor('pondmxtb',0.0d0,1000.d0,pondmxtb,mairg,ifnd)
          do i = 1, ifnd
            pondmxtab(i*2-1) = datepmx(i)
            pondmxtab(i*2)   = pondmxtb(i)
          enddo
          pondmx = pondmxtb(1)
        endif
      endif
      if(swpondmx.eq.0) then
        call rdsdor ('pondmx',0.0d0,1000.0d0,pondmx)
      endif

! --  drainage resistance of surface runoff
      call rdsdor ('rsro',0.001d0,1.0d0,rsro)
      call rdsdor ('rsroexp',0.01d0,10.0d0,rsroexp)
! --  check combination of time-dependent ponding and runoff-resistance
      if(swpondmx.eq.1 .and. rsro.lt.1.0d-02) then
        messag = 'Fatal error: time-dependent threshold for ponding '// &
     &       ' (SWPONDMX=1) requires resistance (RSRO) > 0.01 d-1!'
        call fatalerr ('readswap',messag)
      endif

! --- soil evaporation
      call rdsinr ('swcfbs',0,1,swcfbs)
      if (swcfbs.eq.1) then
        call rdsdor ('cfbs',0.0d0,1.5d0,cfbs)
      else
        cfbs = 1.0d0
      endif
      if (swdivide.eq.1) call rdsdor ('rsoil',0.0d0,1.0d3,rsoil)
      call rdsinr ('swredu',0,2,swredu)
      if (swredu.eq.1) then
        call rdsdor ('rsigni',0.0d0,1.0d0,rsigni)  
        if(rdinqr('cofredbl')) then
          call rdsdor ('cofredbl',0.0d0,1.0d0,cofred)
        else
          call rdsdor ('cofred',0.0d0,1.0d0,cofred)
        endif
      elseif (swredu.eq.2) then
        if(rdinqr('cofredbo')) then
          call rdsdor ('cofredbo',0.0d0,1.0d0,cofred)
        else
          call rdsdor ('cofred',0.0d0,1.0d0,cofred)
        endif
      endif
      cfevappond = 1.25d0
      if(rdinqr('cfevappond')) then
        call rdsdor ('cfevappond',0.0d0,3.0d0,cfevappond)
      endif
 
! --- vertical discretization of soil profile
      call rdainr ('isoillay',1,maho,isoillay,macp,ifnd)
      call rdfinr ('isublay',1,macp,isublay,macp,ifnd)
      call rdfdor ('hsublay',0.d0,10000.d0,hsublay,macp,ifnd)
      call rdfinr ('ncomp',0,macp,ncomp,macp,ifnd)
      nsublay = ifnd
      numlay = isoillay(ifnd)

      ! set thickness of compartments in sublayer (hcomp)
      do i = 1, nsublay
         hcomp(i) = hsublay(i) / dble(ncomp(i))
         if (hcomp(i) .gt. 1000.d0) then
            messag = 'HCOMP exceeding maximum thickness [0..1000,cm]: adjust NCOMP'
            call fatalerr('Readswap',messag)
         end if
      end do

! --- check input soil profile: isoillay and isublay
      if (isoillay(1).ne.1) then
        messag = 'Adjust vert.discr.soil profile: ISOILLAY(1) must be 1'
        call fatalerr('Readswap',messag)
      end if
      if (isublay(1).ne.1) then
        messag = 'Adjust vert.discr. soil profile: ISUBLAY(1) must be 1'
        call fatalerr('Readswap',messag)
      end if
      do i = 2, nsublay
        if ((isoillay(i) < isoillay(i-1)) .or.                          &
     &      (isoillay(i) - isoillay(i-1)) > 1) then
          messag = 'Adjust vert.discr.soil profile: ISOILLAY must be '//&
     &               'in increasing order'
          call fatalerr('Readswap',messag)
        end if
      end do

! --- check input vertical discretization of soil profile
      do i = 1, nsublay
        if (abs(hsublay(i) - hcomp(i)*ncomp(i)) .gt. 1.d-6) lay = 1
        if (i .gt. 1)then
          if( (isoillay(i) .ne. isoillay(i-1) .or.                      &
     &        isoillay(i) .ne. isoillay(i-1) + 1)) lay = 1
        end if
      end do

! --- Soil Hydraulic relation: as MVG-functions or as Tables
      swsophy = 0
      iHWCKmodel(1:numlay) = 1    ! default MvG analytical functions
      if (rdinqr('iHWCKmodel')) then
         call rdfinr ('iHWCKmodel',1,11,iHWCKmodel,maho,numlay)
         ! 1 = uni-modal MvG (default)
         ! 2 = exponential wc(h) and K(h) relationships for testing
         ! 3 = bi-modal MvG relationships
         ! 4-11 = 8 versions of PDI model
      end if
! -   Tables for each soil layer
      if(rdinqr('swsophy')) then
        call rdsinr ('swsophy',0,1,swsophy)
        if(swsophy.eq.1) then
!         not allowed: table-option combined with output with 
!         adjusted vertical discrectization (swdiscrvert=1)
          if(swdiscrvert.eq.1) then
            messag = 'combi of sophys-tables (swsophy=1) and swdiscrvert=1 not allowed'
            call fatalerr ('Readswap',messag)
          endif
          if (numlay .gt. 1) then
            call rdacha ('filenamesophy',filenamesophy,maho,numlay)
          else
            call rdscha ('filenamesophy',filenamesophy(1))
          endif
          messag = 'simulation with tables for soil hydraulic relations'
          call warn ('Readswap',messag,logf,swscre)
        endif
      endif

! --- Switch/flag to reduce capillary rise in compartment below rootzone
      swcaprise = .false.
      if(rdinqr('swcaprise')) then
        call rdslog ('swcaprise',swcaprise)
      endif
      swcapriseoutput = .false.
      if(rdinqr('swcapriseoutput')) then
        call rdslog ('swcapriseoutput',swcapriseoutput)
      endif

! ----  bulk density required for oxygen stress and solute transport
      if(rdinqr('bdens')) then
        call rdfdor ('bdens',100.0d0,10000.0d0,bdens,maho,numlay)
      endif

! --- hysteresis
      call rdsinr ('swhyst',0,2,swhyst)
      if (swhyst.gt.0) then
        call rdsdor ('tau',0.0d0,1.0d0,tau)
      endif

!     combination of soilphysical tables and hysteresis is not possible
      if (swhyst.eq.1 .and. swsophy.eq.1) then
         messag = 'Combination of hysteresis (swhyst=1) and tabulated ' &
     &//'soil physics (swsophy=1) is not possible !'
          call fatalerr('Readswap',messag)
      end if
      
! --- MVG-functions: parameters of functions of each soil layer
      alfa_2(1:maho) = 0.0d0
      npar_2(1:maho) = 1.0d0
      omega_1(1:maho) = 1.0d0
      h0(1:maho) = 0.0d0
      ha(1:maho) = 0.0d0
      apar(1:maho) = 0.0d0
      omega_K(1:maho) = 0.1d0

      if(swsophy.eq.0) then
!##: MH start
         if (RDinqr('SWAP_SOILS_DB')) then
            call RDscha ('SWAP_SOILS_DB', SWAP_SOILS_DB)          ! name of the soils database
            call RDfcha ('SoilNames', SoilNames, maho, numlay)    ! names of the soils per layer

            close (swp)   ! we must close main input file, so that soils database can be accessed by RDinit
            
            idb = getun (120, 200)
            flksatexm  =.false.
            call RDinit (idb, logf, SWAP_SOILS_DB)
!              store all data in database
               call GTSOIL (1,idb,SWAP_SOILS_DB,cdum,dum,dum,dum,dum,dum,dum,dum,dum,dum,flksatexm) ! initialize; most arguments are dummy
               do i = 1, numlay
!                 get values per soil layer
                  call GTSOIL (2,idb,SWAP_SOILS_DB,SoilNames(i),Ores(i),Osat(i),Npar(i),Alfa(i),AlfaW(i),   &
                               Lexp(i),Ksatfit(i),KsatExm(i),H_enpr(i),flksatexm)
               end do
            close (idb)

            call rdinit(swp,logf,swpfilnam)     ! open main input file again
!## MH: end

         else
            call rdfdor ('osat',0.d0,1.0d0,osat,maho,numlay)
            call rdfdor ('ores',0.d0,1.0d0,ores,maho,numlay)
            call rdfdor ('alfa',1.d-4,100.d0,alfa,maho,numlay)
            if (swhyst.gt.0) then
              call rdfdor ('alfaw',1.d-4,100.d0,alfaw,maho,numlay)
            endif
            call rdfdor ('npar',1.001d0,9.d0,npar,maho,numlay)
            call rdfdor ('lexp',-25.d0,25d0,lexp,maho,numlay)
            call rdfdor ('h_enpr',-40.d0,0.d0,h_enpr,maho,numlay)
!           to allow downward compatibility from Swap3.2.23
            if(rdinqr('ksatfit')) then
                 call rdfdor ('ksatfit',1.d-5,1.d5,ksatfit,maho,numlay)
            else
                 call rdfdor ('ksat',1.d-5,1.d5,ksatfit,maho,numlay)
            endif
            flksatexm  =.false.
            if(rdinqr('ksatexm')) then
                 call rdfdor ('ksatexm',1.d-5,1.d5,ksatexm,maho,numlay)
                 flksatexm  =.true.
            endif
!            if (iHWCKmodel == 3) then
            if (any(iHWCKmodel(1:numlay) ==  3) .or. any(iHWCKmodel(1:numlay) ==  6) .or. &
                any(iHWCKmodel(1:numlay) ==  7) .or. any(iHWCKmodel(1:numlay) == 10) .or. &
                any(iHWCKmodel(1:numlay) == 11)) then
               call rdfdor ('alfa_2',1.d-4,100.d0,alfa_2,maho,numlay)
               call rdfdor ('npar_2',1.001d0,9.d0,npar_2,maho,numlay)
               call rdfdor ('omega_1',1.0d-4,0.9999d0,omega_1,maho,numlay)
            end if
            if (any(iHWCKmodel(1:numlay) == 5) .or. any(iHWCKmodel(1:numlay) == 7)) then
               call rdfdor ('h0',-5.0d7,-1.0d5,h0,maho,numlay); h0(1:numlay) = -h0(1:numlay)
            end if
            if (any(iHWCKmodel(1:numlay) ==  8) .or. any(iHWCKmodel(1:numlay) ==  9) .or. &
                any(iHWCKmodel(1:numlay) == 10) .or. any(iHWCKmodel(1:numlay) == 11)) then
               call rdfdor ('h0',-5.0d7,-1.0d5,h0,maho,numlay); h0(1:numlay) = -h0(1:numlay)
               call rdfdor ('ha',-1.0d5,0.0d0,ha,maho,numlay);  ha(1:numlay) = -ha(1:numlay)
               call rdfdor ('apar',-5.0d0,0.0d0,apar,maho,numlay)
               call rdfdor ('omega_K',1.0d-8,0.1d0,omega_K,maho,numlay)
               call rdfinr ('swvapor',0,1,swvapor,maho,numlay)
            end if
         end if      ! if (RDinqr('SWAP_SOILS_DB')) then

         if (flksatexm) then
            messag = 'simulation with additonal Ksat value (Ksatexm)'
            call warn ('Readswap',messag,logf,swscre)
         endif

!        assign sophy-values to paramvg (input of cofgen)
         paramvg = 0.0d0
         do lay = 1,numlay
            paramvg(1,lay) = ores(lay)
            paramvg(2,lay) = osat(lay)
            paramvg(3,lay) = ksatfit(lay)
            paramvg(4,lay) = alfa(lay)
            paramvg(5,lay) = lexp(lay)
            paramvg(6,lay) = npar(lay)
            paramvg(7,lay) = 1.d0 - (1.d0 / paramvg(6,lay))
            if (swhyst.eq.0) then
              paramvg(8,lay) = alfa(lay)
            else
              paramvg(8,lay) = alfaw(lay)
            endif
            paramvg(9,lay) = h_enpr(lay)
            paramvg(10,lay) = -999.d0    ! Ksatexm
            if (flksatexm) then
               if (ksatexm(lay) > ksatfit(lay)) then                          !## MH: 20180111 - only effectively use Ksatexm in case Ksatem > Ksatfit
                  paramvg(10,lay) = ksatexm(lay)
                  help = abs(hthr * paramvg(4,lay))**paramvg(6,lay)
                  help = (1.0d0 + help) ** paramvg(7,lay)
                  relsatthr(lay) = 1.0d0/help
                  term1 = ( 1.0d0-relsatthr(lay)**(1.0/paramvg(7,lay)) )**paramvg(7,lay)
                  ksatthr(lay)=paramvg(3,lay)*(relsatthr(lay)**paramvg(5,lay)) * (1.0d0-term1)*(1.0d0-term1)
                  if(ksatthr(lay).ge.ksatexm(lay)) then
                     write(messag,'(a,i2)') 'ksatexm < ksatthr for layer ',lay
                     call fatalerr ('readswap',messag)
                  endif
               endif
            endif
            paramvg(13,lay) = alfa_2(lay)
            paramvg(14,lay) = npar_2(lay)
            paramvg(15,lay) = 1.d0 - 1.d0 / npar_2(lay)
            paramvg(16,lay) = omega_1(lay)
            paramvg(17,lay) = 1.0d0 - omega_1(lay)
            paramvg(18,lay) = h0(lay)
            paramvg(19,lay) = ha(lay)
            paramvg(20,lay) = apar(lay)
            paramvg(21,lay) = omega_K(lay)
         end do
      endif

!     determine/set BiModal and NoVap
      do i = 1, numlay
         BiModal(i) = .false.
         if (iHWCKmodel(i) ==  3 .OR. iHWCKmodel(i) ==  6 .OR. iHWCKmodel(i) ==  7 .OR.   &
             iHWCKmodel(i) == 10 .OR. iHWCKmodel(i) == 11 ) BiModal(i) = .true.
         NoVap(i) = .true.
         if (iHWCKmodel(i) ==  8 .OR. iHWCKmodel(i) ==  9 .OR. iHWCKmodel(i) == 10 .OR. iHWCKmodel(i) == 11 ) then
            if (swvapor(i) == 0) then
               NoVap(i) = .true.
            else
               NoVap(i) = .false.
            end if
         end if
      end do
      

! --- preferential flow due to macropores
!     first check whether Macropore option is choosen
      flMacroPore = .false.
      call rdsinr ('swmacro',0,1,swmacro)
      if (swmacro.eq.1) then
         flMacroPore = .true.
      endif
      flInitDraBas = .false.
!
      if (flMacroPore) then
!- a. PARAMETERS FOR GEOMETRY MACROPORES
         MaxDepth= 0.d0
         do il = 1, nsublay
            MaxDepth= MaxDepth - hsublay(il)
         enddo
         call rdsdor('Z_Ah',maxdepth,0.d0,Z_Ah)
         call rdsdor('Z_Ic',maxdepth,0.d0,Z_Ic)
         call rdsdor('Z_St',maxdepth,0.d0,Z_St)
         Z_Tp = 0.d0                               ! Adaptation for GEM 
         if(rdinqr('Z_Tp')) then                   !
           call rdsdor('Z_Tp',Z_St,0.d0,Z_Tp)      !
           if (Z_Tp.lt.0.d0 .and. Z_Tp.gt.-dz(1)) Z_Tp = -dz(1) !
         endif                                     ! Adaptation for GEM
         call rdsdor('VlMpStSs',0.0d0,0.5d0,VlMpStSs)
         call rdsdor('PpIcSs',0.0d0,0.99d0,PpIcSs)
         call rdsinr('NumSbDm',0,(MaDm-2),NumSbDm)
         if (NumSbDm.eq.0 .and. PpIcSs.gt.0.d0) then
            messag = ' NumSbDm .eq.0 .and. PpIcSs.gt.0.d0'
            call fatalerr('MacroRead',messag)
         endif
         call rdsdor('DiPoMi',0.1d0,1000.0d0,DiPoMi)
         call rdsdor('DiPoMa',0.1d0,1000.0d0,DiPoMa)
         call rdsdor('PndmxMp',0.0d0,10.0d0,PndmxMp)    
!   - optional: 
         call rdsdor('PowM',0.0d0,100.0d0,PowM)    ! default 1.0
         call rdsdor('Rzah',0.0d0,1.0d0,Rzah)      ! default 0.0
         call rdsdor('Spoint',0.0d0,1.0d0,Spoint)  ! default 1.0
         call rdsinr('SwPowM',0,1,SwPowM)          ! default 0
         call rdsdor('PndmxMp',0.0d0,10.0d0,PndmxMp)    
         Z_MB50 = 0.5d0*(Z_Ic+Z_St)
         if(rdinqr('Z_MB50')) then
           call rdsdor('Z_MB50',Z_St,Z_Ic,Z_MB50)
         endif
!- b. PARAMETERS FOR SHRINKAGE CHARACTERISTICS
         call rdfinr('SwSoilShr',0,2,SwSoilShr,MaHo,numlay)
         call rdfinr('SwShrInp',1,3,SwShrInp,MaHo,numlay)
         call rdfdor('ThetCrMP',0.0d0,1.0d0,ThetCrMp,MaHo,numlay)
         call rdfdor('GeomFac',0.0d0,10.0d0,GeomFac,MaHo,numlay)
         call rdfdor('ShrParA',0.0d0,10.0d0,ShrParA,MaHo,numlay)
         call rdfdor('ShrParB',-10.0d0,100.0d0,ShrParB,MaHo,numlay)
         call rdfdor('ShrParC',-10.0d0,100.0d0,ShrParC,MaHo,numlay)
         call rdfdor('ShrParD',0.0d0,100.0d0,ShrParD,MaHo,numlay)
         call rdfdor('ShrParE',-10.0d0,10.0d0,ShrParE,MaHo,numlay)
         call rdsdor('ZnCrAr',Z_Ah,0.0d0,ZnCrAr)
!- c. PARAMETERS FOR SORPTIVITY
         call rdfinr('SwSorp',1,2,SwSorp,MaHo,numlay)
         call rdfdor('SorpFacParl',0.0d0,1.0d2,SorpFacParl,MaHo,numlay)
         call rdfdor('SorpMax',0.0d0,100.0d0,SorpMax,MaHo,numlay)
         call rdfdor('SorpAlfa',-10.0d0,10.0d0,SorpAlfa,MaHo,numlay)
!- d. SHAPE FACTOR for saturated exchange between macropores and matrix
         call rdsdor('ShapeFacMp',0.0d0,100.0d0,ShapeFacMp)
!- e. CRITICAL value for undersaturation volume
         call rdsdor('CritUndSatVol',0.0d0,10.0d0,CritUndSatVol)
         call rdsinr('SwDarcy',0,1,SwDarcy)
!- f. PARAMETERS FOR RAPID DRAINAGE
!       only possible when at least one drainege level!!!
         call rdsinr('SwDrRap',0,1,SwDrRap)
         if (SwDrRap.eq.1) then
!          at this moment for 1 system only, may be extended to other systems 
            call rdsdor ('RapDraResRef',0.d0,1.d10,RapDraResRef(1))
            do il= 2,madr
               RapDraResRef(il) = RapDraResRef(1)
            enddo
            call rdsdor('RapDrareaExp',0.0d0,100.0d0,RapDraReaExp)
            call rdsinr('NumLevRapDra',1,5,NumLevRapDra)
!
!-    set flag for initialisation of drainage basis rapid drainage
            flInitDraBas = .true.
         endif
      endif

! -   snow and frost conditions 
      call rdsinr ('swsnow',0,1,swsnow)
      if (swsnow.eq.1) then
        call rdsdor ('snowinco',0.0d0,1000.0d0,snowinco)
        call rdsdor ('teprrain',0.0d0,10.0d0,teprrain)
        call rdsdor ('teprsnow',-10.0d0,0.0d0,teprsnow)
        call rdsdor ('snowcoef',0.0d0,10.0d0,snowcoef)
        swsublim = 0
        if(rdinqr('swsublim')) then
          call rdsinr ('swsublim',0,1,swsublim)
        endif
      endif
      call rdsinr ('swfrost',0,1,swfrost)
      if (swfrost.eq.1) then
        call rdsdor ('tfroststa',-10.0d0,5.0d0,tfroststa)
        call rdsdor ('tfrostend',-10.0d0,5.0d0,tfrostend)
      endif
!     combination of frost and macropore-flow is not possible (yet)
      if (swfrost.eq.1 .and. flMacroPore) then
         messag = 'Combination of frost (swfrost=1) and macropore-flow '&
     &//'(swmacro=1) is not operational !'
          call fatalerr('Readswap',messag)
      end if

!     combination of snow and Et variation during the day (SwEtSine)
      if (swmetdetail.eq.1 .and. swsnow.eq.1) then
        messag = 'In case of snow the sublimation is not simulated'//   &
     & 'correctly if short time meteorological records are used'//      &
     & '(swmetdetail=1 and swsnow=1), please adapt input !'
        call warn ('Readswap',messag,logf,swscre)
      end if

! --- parameters numerical scheme
      call rdsdor ('dtmin', 1.0d-7,0.1d0,dtmin)
      call rdsdor ('dtmax', dtmin, 1.0d0,dtmax)
      call rdsdor ('gwlconv',1.0d-5,1000.0d0,gwlconv)
      call rdsdor ('CritDevPondDt',1.0d-6,1.0d-01,CritDevPondDt)
! Convergence criteria (optional input)
      CritDevh1Cp = 1.0d-2
      if(rdinqr('CritDevh1Cp')) then
        call rdsdor ('CritDevh1Cp',1.0d-10,1.0d3,CritDevh1Cp)
      endif
      CritDevh2Cp = 1.0d-1
      if(rdinqr('CritDevh2Cp')) then
        call rdsdor ('CritDevh2Cp',1.0d-10,1.0d3,CritDevh2Cp)
      endif
! flag to generate additional output about convergence-warnings from subr Headcalc
      fldumpconvcrit = .false.
      if(rdinqr('fldumpconvcrit')) then
        call rdslog ('fldumpconvcrit',fldumpconvcrit)
      endif
! Maximum number of iterations [5,100]
      call rdsinr ('MaxIt',5,100,MaxIt)
! Maximum number of back track cycles within an iteration cycle [1,10]
      call rdsinr ('MaxBackTr',1,10,MaxBackTr)
! Maximum number of iterations: no input
      msteps = 100000000
! Switch for mean of hydraulic conductivity 
!  SwkMean=1,2:unweighted arithmic mean,weighted arithmic mean
!  SwkMean=3,4:unweighted geometric mean, weighted geometric mean
!  SwkMean=5,6:unweighted harmonic mean, weighted harmonic mean
      call rdsinr ('SWkmean',1,6,SWkmean)
! Switch for implicit solution with hydraulic conductivity: 0 = explicit, 1 = implicit
      call rdsinr ('SwkImpl',0,1,SwkImpl)
! Maximum cputime, introduced to be able to interrupt (near) endless iterations
      flMaxIterTime = .false.
      if(rdinqr('flMaxIterTime')) then
        call rdslog ('flMaxIterTime',flMaxIterTime)
        MaxIterTime = 2419200                ! 4 weeks = 60*60*24*7*4=2419200 secs)
        if(flMaxIterTime) then
          call rdsinr ('MaxIterTime',1,2419200,MaxIterTime)
        endif
      endif

! =================================================================
!     Lateral drainage section
!     extended or basic drainage
      call rdsinr ('swdra',0,2,swdra)
      if (swdra .ne. 0) call rdscha ('drfil',drfil)
      if (SwDrRap.eq.1 .and. SwDra.eq.0) then
          messag = ' There are no drainage levels, so rapid drainage'// &
     &'is not possible !'
          call fatalerr('MacroRead',messag)
      endif

!     runon from external source (field)
      call rdsinr ('swrunon',0,1,swrunon)
      if (swrunon .eq. 1) then
        call rdscha ('rufil',rufil)
      endif

!     output-options drainage fluxes, surface reservoir 
      if (swdra.eq.2) then
        call rdsinr ('swdrf',0,1,swdrf)
        call rdsinr ('swswb',0,1,swswb)
      endif


! =================================================================
!     Bottom boundary section

! --- Initialise
      do i = 1,2*mabbc
        gwltab(i) = 0.0D0
        hbotab(i) = 0.0D0
        qbotab(i) = 0.0D0
        haqtab(i) = 0.0D0
        cseeptab(i) = 0.0D0
      end do

! --- option for input of bottom boundary condition
      call rdsinr ('swbbcfile',0,1,swbbcfile)
      if (swbbcfile .eq. 1) then
        call rdscha ('bbcfil',bbcfil)
      endif

! =================================================================
!     Section heat flow

! --- Switch whether simulation includes heat simulation or not
      call rdsinr ('swhea',0,1,swhea)

      if (swhea .eq. 1) then
! ---   analytical or numerical method 
        call rdsinr ('swcalt',1,2,swcalt)

        if (swcalt.eq.1) then
          call rdsdor ('tampli',0.0d0,50.0d0,tampli)
          call rdsdor ('tmean',-10.0d0,30.0d0,tmean)
          call rdsdor ('timref',0.0d0,366.0d0,timref)
          call rdsdor ('ddamp',1.0d0,500.0d0,ddamp)
        else
          call rdador ('zh',-1.0d5,0.0d0,zh,macp,ifnd)
          call rdfdor ('tsoil',-50.0d0,50.0d0,tsoil,macp,ifnd)
          call rdfdor ('psand',0.0d0,1.0d0,psand,maho,numlay)
          call rdfdor ('psilt',0.0d0,1.0d0,psilt,maho,numlay)
          call rdfdor ('pclay',0.0d0,1.0d0,pclay,maho,numlay)
          call rdfdor ('orgmat',0.0d0,1.0d0,orgmat,maho,numlay)
          nheat = ifnd

!   -     top boundary temperature
          call rdsinr ('SwTopbHea',1,2,SwTopbHea)
          Tsoilfile = ' '
          if (swtopbhea .eq. 2) then
            call rdscha ('TSoilFile',tsoilfile)
          endif

!   -     bottom boundary temperature
          call rdsinr ('SwBotbHea',1,2,SwBotbHea)
          if (SwBotbHea.eq.2) then
             call rdatim ('datet',dates,mabbc,ifnd)
! -     at least one date must be within simulation period
             call checkdate(ifnd,dates,tend,tstart,'datet ',            &
     &                      'readswap/swbotbhea=2    ')
             call rdfdor('tbot',-50.0d0,50.d0,tbot,mabbc,ifnd)
!
             do i = 1, ifnd
               tembtab(i*2)   = tbot(i)
               tembtab(i*2-1) = dates(i)
             enddo
!
           endif
!
        endif
      endif

! -   fatal error when frost or snow is simulated without heat flow
      if ((swfrost.eq.1 .or. swsnow.eq.1) .and. swhea.eq.0) then
        messag = 'In case of snow or frost the soil heat flow should '//&
     &  'be simulated! Adapt .swp input file.'
        call fatalerr ('readswap',messag)
      endif

! =================================================================
!     Section solute transport

! --- Switch whether simulation includes solute transport or not
      call rdsinr ('swsolu',0,1,swsolu)

      if (swsolu .eq. 1) then

! ---   top boundary and initial condition
        call rdsdor ('cpre',0.0d0,100.0d0,cpre)
        call rdsdor ('cdrain',0.0d0,100.0d0,cdrain)

!       bottom condition with 3 options:
!       SWBOTBC = 0: lateral (cdrain as input), bottom set to lateral (cseep=cdrain)
!       SWBOTBC = 1: separate input of bottom conc as cseep
!       SWBOTBC = 2: separate input of dynamic bottom conc as (datec,cseeparr)
        swbotbc = 0
        if (rdinqr('swbotbc')) then
          call rdsinr ('swbotbc',0,2,swbotbc)
        endif
        if (swbotbc.eq.0) then
          cseep = cdrain
        elseif (swbotbc.eq.1) then
          call rdsdor ('cseep',0.0d0,100.0d0,cseep)
        elseif (swbotbc.eq.2) then
          call rdatim ('datec',dates,mabbc,ifnd)
          call rdfdor ('cseeparr',0.0d0,100.0d0,cseeparr,mabbc,ifnd)
          do i = 1, ifnd
            cseeptab(i*2) = cseeparr(i)
            cseeptab(i*2-1) = dates(i)
          end do
        endif

        if (swinco.ne.3) then
          call rdador ('zc',-1.0d5,0.0d0,zc,macp,ifnd)
          call rdfdor ('cml',0.0d0,1000.0d0,cml,macp,ifnd)
          nconc = ifnd
        endif

! ---   diffusion, dispersion and solute uptake by roots
        call rdsdor ('ddif',0.0d0,10.0d0,ddif)
        if (rdinar('ldis')) then
          call rdfdor ('ldis',0.0d0,100.d0,ldis,maho,numlay)
        else
          call rdsdor ('ldis',0.0d0,100.d0,ldis(1))
        endif
        call rdsdor ('tscf',0.0d0,10.0d0,tscf)

! ---   sorption
        call rdsinr ('swsp',0,1,swsp)
        if (swsp .eq. 1) then
          call rdsdor ('frexp',0.0d0,   10.0d0,frexp)  
          call rdsdor ('cref', 0.0d0, 1000.0d0,cref)
          if (rdinar('kf')) then
            call rdfdor ('kf',0.0d0,10000.0d0,kf,maho,numlay)
          else
            call rdsdor ('kf',0.0d0,10000.0d0,kf(1))
          endif
          if (rdinar('bdens')) then
            call rdfdor ('bdens',100.0d0,10000.0d0,bdens,maho,numlay)
          else
            call rdsdor ('bdens',100.0d0,10000.0d0,bdens(1))
          endif
        else
          cref  = 1.0d0
          frexp = 1.0d0
          do lay = 1,numlay
            kf(lay) = 0.0d0
          end do
        endif

! ---   decomposition
        call rdsinr ('swdc',0,1,swdc)
        if (swdc .eq. 1) then
          if (rdinar('decpot')) then
            call rdfdor ('decpot',0.0d0,10.0d0,decpot,maho,numlay)
          else
            call rdsdor ('decpot',0.0d0,10.0d0,decpot(1))
          endif
          call rdsdor ('gampar',0.0d0,0.5d0,gampar)
          call rdsdor ('rtheta',0.0d0,0.4d0,rtheta)
          call rdsdor ('bexp',0.0d0,2.0d0,bexp)
          if (rdinar('fdepth')) then
            call rdfdor ('fdepth',0.0d0,1.0d0,fdepth,maho,numlay)
          else
            call rdsdor ('fdepth',0.0d0,1.0d0,fdepth(1))
          endif
        else 
          gampar = 0.0d0
          rtheta = 0.5d0
          bexp = 0.0d0
          do lay = 1,numlay
            decpot(lay) = 0.0d0
            fdepth(lay) = 0.0d0
          end do
        endif

! ---   breakthrough
        call rdsinr ('swbr',0,1,swbr)
        if (swbr .eq. 0) then
          daquif = 100.0d0
          poros = 1.0d0
          kfsat = 0.0d0
          decsat = 0.0d0
        else
          call rdsdor ('daquif',0.0d0,10000.0d0,daquif)
          call rdsdor ('poros',0.0d0,0.6d0,poros)
          call rdsdor ('kfsat',0.0d0,100.0d0,kfsat)
          call rdsdor ('decsat',0.0d0,10.0d0,decsat)
          call rdsdor ('cdraini',0.0d0,100.0d0,cdrain)
          cseep = cdrain
        endif

      endif

! =================================================================
!     Section Ageing according to Goode (1996): 
!      "Direct simulation of groundwater age, WRR vol.32, p 289-296"
      flAgeTracer = .false.
      if (rdinqr('flAgeTracer')) then
        call rdslog ('flAgeTracer',flAgeTracer)
        if (flAgeTracer) then

! ---     top and bottom boundary and initial condition
          cpre = 0.0d0
          cirr = 0.0d0
          cdrain = 0.0d0
          cseep = 0.0d0
          
          if (swinco.ne.3) then
            call rdador ('zc',-1.0d5,0.0d0,zc,macp,ifnd)
            call rdfdor ('cml',0.0d0,1000.0d0,cml,macp,ifnd)
            nconc = ifnd
          endif

! ---     diffusion, dispersion and uniform solute uptake by roots
          call rdsdor ('ddif',0.0d0,10.0d0,ddif)
          if (rdinar('ldis')) then
            call rdfdor ('ldis',0.0d0,100.d0,ldis,maho,numlay)
          else
            call rdsdor ('ldis',0.0d0,100.d0,ldis(1))
          endif
          tscf = 1.0d0

!         warnings and errors
          messag = 'simulation with Age Tracer option (flAgeTracer)'
          call warn ('Readswap',messag,logf,swscre)
          if(swSolu.eq.1) then
            write(messag,'(3a)')                                        &
     &      'Combination of solute transport and Ageing is not allowed',&
     &      '(flAgeTracer=.true.  AND   swSolu=1 is not allowed !' 
            call fatalerr ('Readswap',messag)
          endif
          if(swSnow.eq.1) then
            write(messag,'(3a)')                                        &
     &       'Combination of snow and Ageing is not allowed',           &
     &       '(flAgeTracer=.true.  AND   swSnow=1 is not allowed !'
            call fatalerr ('Readswap',messag)
          endif
        endif
      endif

! =================================================================
!     Section bottom boundary condition

      if (swbbcfile .eq. 1) then
        close (swp)
        bbc = getun2 (10,90,2)
        filnam = trim(pathwork)//trim(bbcfil)//'.bbc'
        call rdinit(bbc,logf,filnam)
      endif

! --- option for bottom boundary condition
      call rdsinr ('swbotb',1,8,swbotb)

! --- given groundwaterlevel
      if (swbotb.eq.1) then
        call rdatim ('date1',dates,mabbc,ifnd)
        call rdfdor ('gwlevel',-10000.0d0,1000.d0,gwlevel,mabbc,ifnd)

! -     check: at least one date must be within simulation period
        call checkdate(ifnd,dates,tend,tstart,'date1 ',                 &
     &                      'readswap//swbotb=1      ')
! -     check if gwlevel is above soil surface, eliminate certain combinations
        flsat = .false.
        do i = 1,ifnd
          if (gwlevel(i) .gt. -0.5d0*(hcomp(1)))  flsat = .true.
        enddo
        if(flsat) then
           if (swKimpl.eq.1) then
              write(messag,'(3a)')                                      &
     &         'Groundwaterlevel above soil surface, combined with ',   &
     &         'Implicit Conductivity (swbotb=1 AND swKimpl=1 AND  ',   &
     &         'gwl>z(1))  is not allowed !'
              call fatalerr ('Readswap',messag)
           endif
           if (flMacroPore) then
              write(messag,'(3a)')                                      &
     &         'Groundwaterlevel above soil surface, combined with ',   &
     &         'MacroPore flow (swbotb=1 AND swMacro=1 AND gwl>z(1)) ', &
     &         'is not allowed !' 
              call fatalerr ('Readswap',messag)
           endif
           if (swfrost.eq.1) then
              write(messag,'(3a)')                                      &
     &         'Groundwaterlevel above soil surface, combined with ',   &
     &         'frost conditions (swbotb=1 AND swforst=1 AND gwl>z(1))',&
     &         ' not well tested, be aware of balance errors !' 
              call warn ('Readswap',messag,logf,swscre)
           endif
        endif

! -     store values in gwltab
        do i = 1,ifnd
          gwltab(i*2) = gwlevel(i) 
          gwltab(i*2-1) = dates(i)
        enddo
      endif

! --- regional bottom flux is given                             
      if (swbotb.eq.2) then
        call rdsinr ('sw2',1,2,sw2)
        if (sw2 .eq. 1) then
          call rdsdor ('sinave',-10.0d0,10.0d0,sinave)
          call rdsdor ('sinamp',-10.0d0,10.0d0,sinamp)
          call rdsdor ('sinmax',0.d0,366.d0,sinmax)
        else
! -       read tabular data
          call rdatim ('date2',dates,mabbc,ifnd)
          call rdfdor ('qbot2',-100.0d0,100.0d0,qboti,mabbc,ifnd)
! -     at least one date must be within simulation period
          call checkdate(ifnd,dates,tend,tstart,'date2 ',               &
     &                      'readswap//swbotb=2      ')
! -       fill qbotab table
          do i = 1,ifnd
            qbotab(i*2) = qboti(i) 
            qbotab(i*2-1) = dates(i)
          enddo
        endif
      endif

! --- calculated flux through the bottom of the profile
      if (swbotb.eq.3) then

!       Switch for implicit solution with lower boundary option 3 (Cauchy): 0 = explicit, 1 = implicit
        call rdsinr ('swbotb3Impl',0,1,swbotb3Impl)

        call rdsdor ('shape',0.0d0,1.0d0,shape)
        if (swbotb3Impl.eq.1 .and. abs(shape-1.0d0).gt.1.0d-7) then
          write(messag,'(a)')                                           &
     &    ' Possible lower boundary inconsistency using SwBotb=3: ',    &
     &    '  Combination of swbotb3Impl AND shape not equal 1.0',       &
     &    '  This is not recommended. Suggestion is: swbotb3Impl=0'
          call warn ('readswap',messag,logf,swscre)
        endif

        call rdsdor ('hdrain',-1.0d4,0.0d0,hdrain)
        call rdsdor ('rimlay',0.0d0,1.0d5,rimlay)

!       Switch to suppress addition of vertical resistance 
!                      between bottom of model and groundwater level
        call rdsinr ('SwBotb3ResVert ',0,1,SwBotb3ResVert)

        call rdsinr ('sw3',1,2,sw3)
        if (sw3 .eq. 1) then
          call rdsdor ('aqave',-10000.0d0,1000.0d0,aqave)
          call rdsdor ('aqamp',0.0d0,1000.0d0, aqamp)
          call rdsdor ('aqtmax',0.0d0,366.d0,aqtmax)
          call rdsdor ('aqper',0.0d0,366.0d0,aqper)
        else
! -       read tabular data
          call rdatim ('date3',dates,mabbc,ifnd)
          call rdfdor ('haquif',-10000.0d0,1000.d0,haquif,mabbc,ifnd)
! -     at least one date must be within simulation period
          call checkdate(ifnd,dates,tend,tstart,'date3 ',               &
     &                   'readswap//swbotb=3      ')
! -       fill haqtab table
          do i = 1,ifnd
            haqtab(i*2) = haquif(i) 
            haqtab(i*2-1) = dates(i)
          enddo
        endif
        call rdsinr ('sw4',0,1,sw4)
        if (sw4 .eq. 1) then
          if (swbotb3Impl.eq.1) then
            write(messag,'(3a)')                                        &
     &       ' Implicit solution of Cauchy, combined with fluxes ',     &
     &       'is active !(swbotb3Impl=1 AND sw4=1)'
            call warn ('readswap',messag,logf,swscre)
          endif
! -       read tabular data
          call rdatim ('date4',dates,mabbc,ifnd)
          call rdfdor ('qbot4',-100.0d0,100.d0,qboti,mabbc,ifnd)
! -     at least one date must be within simulation period
          call checkdate(ifnd,dates,tend,tstart,'date4 ',               &
     &                   'readswap//swbotb=3/sw4=1')
! -       fill qbotab table
          do i = 1,ifnd
            qbotab(i*2) = qboti(i) 
            qbotab(i*2-1) = dates(i)
          enddo
        endif
      endif

! --- flux-groundwater level relationship 
      if (swbotb.eq.4) then
        call rdsinr ('swqhbot',1,2,swqhbot)
        if (swqhbot.eq.1) then
          call rdsdor ('cofqha',-100.0d0,100.0d0,cofqha)
          call rdsdor ('cofqhb',-1.0d0,1.0d0,cofqhb)
          swcofqhc = 0
          cofqhc = 0.0d0
          if (rdinqr('cofqhc')) then
            swcofqhc = 1
            call rdsdor ('cofqhc',-10.0d0,10.0d0,cofqhc)
          endif
        else if (swqhbot.eq.2) then
          call rdador ('qtab',-100.0d0,100.d0,qhtab,mabbc,ifnd)
          call rdfdor ('htab', -1.0d4, 0.0d0, hhtab,mabbc,ifnd)
          do i = 1,ifnd
            qbotab(i*2) = qhtab(i) 
            qbotab(i*2-1) = abs(hhtab(i))
          enddo
        endif
      endif

! --- pressure head of lowest compartment is given 
      if (swbotb.eq.5) then
        call rdatim ('date5',dates,mabbc,ifnd)
        call rdfdor('hbot5',-1.0d10,1000.d0,hbottom,mabbc,ifnd)
! -     at least one date must be within simulation period
        call checkdate(ifnd,dates,tend,tstart,'date5 ',                 &
     &                 'readswap//swbotb=5      ')
! -     store pressure head values in hbotab
        do i = 1, ifnd
          hbotab(i*2) = hbottom(i) 
          hbotab(i*2-1) = dates(i)
        end do
      endif

! --- free outflow (lysimeter, suction candle possible) 
      if (swbotb.eq.8) then
        if (rdinqr('hplate')) then
          call rdsdor('hplate',-1000.d0,0.0d0,hplate)
        else
          hplate = 0.0d0
        endif     
      endif

      if (swbbcfile .eq. 0) then
        close (swp)
      elseif (swbbcfile .eq. 1) then
        close (bbc)
      endif

! -   Tables for each soil layer
      if(swsophy.eq.1) then
          ientrytablay = 0
          do lay = 1,numlay
            filnam = filenamesophy(lay)
            sol = getun2 (50,90,2)
            call rdinit(sol,logf,filnam)
            call rdador ('headtab',-1.0d15,1.0d15,headtab,matab,ifnd)
            call rdfdor ('thetatab',0.0d0,1.0d0,thetatab,matab,ifnd)
            call rdfdor                                                 &
     &                ('conductab',0.0d0,10000.0d0,conductab,matab,ifnd)
            close (sol)
            numtablay(lay) = ifnd

!           verify incremental sequence of values
            do i = 2,numtablay(lay)
              if((headtab(i)-headtab(i-1)) .le. nihil) then
                messag = ' No incremental values for head soil physic'  &
     &          //'in input file '//trim(filnam)
              call fatalerr('Readswap',messag)
              endif
              if((thetatab(i)-thetatab(i-1)) .le. vsmall) then
                messag = ' No incremental values for theta soil physic' &
     &          //'in input file '//trim(filnam)
              call fatalerr('Readswap',messag)
              endif
              if((conductab(i)-conductab(i-1)) .le. nihil**2) then
                messag = ' No incremental values for conduc soil physic'&
     &          //'in input file '//trim(filnam)
              call fatalerr('Readswap',messag)
              endif
            enddo

!           preprocess table
!---- sorting of arrays to an ascending sequence

!## MH start
            ! ln-transformation of h and K if global flag do_ln_trans = true
            ! works only if all headtab values are < 0
            if (do_ln_trans) then
               do i = 1, numtablay(lay)
                  headtab(i)   = -dlog(-headtab(i)+1.0d0)     ! minus, so that table remains in stricly increasing order
                  conductab(i) =  dlog(conductab(i))
               end do
            end if
!## MH end

            Do i=1,numtablay(lay)-1
               Do j=i+1,numtablay(lay)
                  If(headtab(i) .gt. headtab(j))Then
                     Dum          = headtab(j)
                     headtab(j)   = headtab(i)
                     headtab(i)   = Dum
                     Dum          = thetatab(j)
                     thetatab(j)  = thetatab(i)
                     thetatab(i)  = Dum
                     Dum          = conductab(j)
                     conductab(j) = conductab(i)
                     conductab(i) = Dum
                  End If
               End Do
            End Do
            ientrytablay(lay,0)=numtablay(lay)
            do i = 1,numtablay(lay)

!## MH start
              !if(headtab(i).gt.-1.0d-1)then
              !   j=0
              !else
              !   j = int(1000*(log10(-headtab(i))+1.d0))+1
              !end if
              if (do_ln_trans) then
                  dummy = -(dexp(-headtab(i))-1.0d0)
              else
                  dummy = headtab(i)
              end if
              if(dummy.gt.-1.0d-5)then
                 j=0
              else
                 j = int(1000*(log10(-dummy)+1.d0))+4001    ! +4001 to ensure proper functioning for near-zero entries of pressure head (say -1.0d-5 and lower)
              end if
!## MH end

              ientrytablay(lay,j) = i
              sptablay(1,lay,i) = headtab(i)
              sptablay(2,lay,i) = thetatab(i)
              sptablay(3,lay,i) = conductab(i)
            enddo
            ientrytablay(lay,1) = 0
            do j=matabentries-1,1,-1
               if(ientrytablay(lay,j) .eq. 0)                           &
     &            ientrytablay(lay,j) = ientrytablay(lay,j+1)
!## MH: additional check should be: all integers between 0 and n-1 should be included in table ientrytablay for best performance ???
            end do
            call PreProcTabulatedFunction(1,                            &
     &                       numtablay(lay),headtab,thetatab,dydx,sigma)
            do i = 1,numtablay(lay)
              sptablay(4,lay,i) = dydx(i)
              sptablay(6,lay,i) = sigma(i)   !## MH: new
            enddo
            call PreProcTabulatedFunction(2,                            &
     &                      numtablay(lay),headtab,conductab,dydx,sigma)
            do i = 1,numtablay(lay)
              sptablay(5,lay,i) = dydx(i)
              sptablay(7,lay,i) = sigma(i)   !## MH: new
            enddo
!           tables must have values for a head=0
            if (sptablay(1,lay,numtablay(lay)) .lt. -1.d-20) then
               messag = ' No values for head=0 in tabulated soil physic'&
     &          //'in input file '//trim(filnam)
              call fatalerr('Readswap',messag)
            end if
          enddo
      endif


! =================================================================
!     Read data of drainage input file
      if (swdra.eq.1) call rddrb ()

! =================================================================
! -   read data of fixed irrigation from separate file
      if (swirgfil .eq. 1 .and. swirfix .eq. 1) then
        filnam = trim(pathwork)//trim(irgfil)//'.irg'
        irg = getun2 (10,90,2)
        call rdinit(irg,logf,filnam)
        call rdatim ('irdate',irdate,mairg,ifnd)
        call rdfdor ('irdepth',0.d0,1000.d0,irdepth,mairg,ifnd)
        call rdfdor ('irconc',0.d0,1000.d0,irconc,mairg,ifnd)
        call rdfinr ('irtype',0,1,irtype,mairg,ifnd)
        close (irg)
! -     at least one date must be within simulation period
        call checkdate(ifnd,irdate,tend,tstart,'irdate',                &
     &                 'readswap//swirgfil=1     ')
! ---   convert mm irrigation to cm
        do i = 1, ifnd
          irdepth(i) = irdepth(i) / 10.d0
        end do
      endif

! --- read data from file with results from previous simulation
! --- set initial ponding conditions
      if (swinco.eq.3) then
        ini = getun2 (10,90,2)
        call rdinit(ini,logf,inifil)
!        write (messag,107)                                              &
!     &          '*  I/O of variables from file (SWINCO=3)',             &
!     &          '   should be considered carefully',                    &
!     &          '   this option is incomplete and poorly tested!'
 107    format (3a)
!        call warn ('readswap',messag,logf,swscre)

        call rdsdor ('ssnow',0.0d0,1000.0d0,ssnow)
        if (swsnow.ne.1) then
          write (messag,107)                                            &
     &         'No Simulation of snow, therefore : Initial',            &
     &         '  storage of snow set to 0.0',                          &
     &         '  neglecting ssnow-value from file '//trim(inifil)
          call warn ('readswap',messag,logf,swscre)
          ssnow = 0.0d0
        endif
        call rdsdor ('slw',0.0d0,1000.0d0,slw)
        call rdsdor ('pond',0.0d0,100.0d0,pond)
        call rdador ('z_h',-1.0d5,0.0d0,zi,macp,ifnd)
        call rdfdor ('h',-1.0d10,1.0d4,h,macp,ifnd)
        nhead = ifnd
        if (swhea.eq.1 .and. swcalt.eq.2) then
          call rdador ('z_Tsoil',-1.0d5,0.0d0,zh,macp,ifnd)
          call rdfdor ('Tsoil',-50.0d0,50.0d0,Tsoil,macp,ifnd)
        endif
        if (swsolu.eq.1 .or. flAgeTracer) then
          call rdador ('z_Cml',-1.0d5,0.0d0,zc,macp,ifnd)
          call rdfdor ('Cml',0.0d0,1000000.0d0,cml,macp,ifnd)
        endif
!       surface water level
        if(swdra.eq.2) then
          if(rdinqr('wls')) then
            call rdsdor ('wls',-1000.0d0,1000.0d0,wls)
            messag = 'Initial wls read from file (SWINCO=3)'//          &
     &      ' not implemented yet'
            call warn ('readswap',messag,logf,swscre)
          endif
        endif
!       write soil evaporation reservoirs
        if(swredu.eq.1) then
          messag = 'Initial ldwet read from file (SWINCO=3)'//          &
     &     ' if absent, then default of 0.0 is assumed'
          call warn ('reduceva',messag,logf,swscre)
          if(rdinqr('ldwet')) then
!           Time after significant rainfall (d)  (Black)
            call rdsdor ('ldwet',0.0d0,300.0d0,ldwet)
          else
            ldwet = 1.0d0
          endif
        endif
        if(swredu.eq.2) then
          messag = 'Initial spev and saev read from file (SWINCO=3)'//  &
     &     ' if absent, then default of 0.0 is assumed'
          call warn ('reduceva',messag,logf,swscre)
          if(rdinqr('spev')) then
!           Rainfall excess (cm) (Boesten/Stroosnijder)
            call rdsdor ('spev',0.0d0,1000.0d0,spev)
          else
            spev = 0.0d0
          endif
          if(rdinqr('saev')) then
!           Cumulative actual evaporation (cm) (Boesten/Stroosnijder)
            call rdsdor ('saev',0.0d0,1000.0d0,saev)
          else
            saev = 0.0d0
          endif
        endif
!       length of final timestep (d)
        if(rdinqr('dt')) then
          !call rdsdor ('dt',dtmin,dtmax,dt) ! Aanpassing H.M. Mulder 2016-10-31 probleem bij TTUTIL tijdens inlezen dt (strenge controle niet nodig vanwege ini-file geproduceerd met SWAP)
          call rdsdou ('dt',dt)
          if (dt .lt. dtmin) dt = dtmin
          if (dt .gt. dtmax) dt = dtmax
        endif

        call rdsinr ('swIrrigate',0,1,swIrrigate)
        if (swIrrigate .eq. 0) then
          flIrrigate = .false.
        else
          flIrrigate = .true.
        endif
        
!       close file
        close(ini)

!         macropore variables
        if(flMacroPore) then
          write (messag,107)                                            &
     &          '*  I/O of macropore variables from file ',             &
     &          '   (SWINCO=3) not implemented yet !'
          call warn ('readswap',messag,logf,swscre) 
        endif

      endif

! --- runon from external file
      if (swrunon .eq. 1) then
        runf = getun2 (10,90,2)
        call rdinit(runf,logf,trim(rufil))
        call rdador ('runoff',0.d0,1000.d0,runonarr,maday,ifnd)
        flrunon = .true.
        close (runf)
      else
        flrunon = .false.
      endif

! --- read soil surface temperatures
      if (swtopbhea .eq. 2) then
        filnam = trim(pathwork)//trim(tsoilfile)//'.tss'
        tss = getun2 (10,90,2)
        call rdinit(tss,logf,filnam)
        call rdatim ('datet',dates,mabbc,ifnd)
        call rdfdor('ttop',-50.0d0,50.d0,ttop,mabbc,ifnd)
        do i = 1, ifnd
           temtoptab(i*2)   = ttop(i)
           temtoptab(i*2-1) = dates(i)
        enddo
        close (tss)
      endif

!     final checks (typically for bi-modal situation
      if (any(iHWCKmodel(1:numlay) > 1)) then
         if (swhyst == 1) then
            messag = 'Combination of hysteresis (swhyst=1) and iHWCKmodel > 1 is not (yet) possible!'
            call fatalerr ('readswap', messag)
         end if
         if (swsophy == 1) then
            messag = 'Combination of tabulated input (swsophy=1) and iHWCKmodel > 1 is impossible!'
            call fatalerr ('readswap', messag)
         end if
         if (swdiscrvert == 1) then
            messag = 'Combination of output coarse discretization (swdiscrvert=1) and iHWCKmodel > 1 is not (yet) possible!'
            call fatalerr ('readswap', messag)
         end if
         !if (swkimpl == 1) then
         !   messag = 'Combination of swkimpl=1 and iHWCKmodel > 1 is not (yet) possible!'
         !   call fatalerr ('readswap', messag)
         !end if
         if (flMacroPore) then
            messag = 'Combination of macro-pore option (flMacroPore) and iHWCKmodel > 1 is impossible!'
            call fatalerr ('readswap', messag)
         end if
         if (minval(h_enpr) < 0.0d0) then
            messag = 'User-input of h_enpr < 0 not used in combination with iHWCKmodel > 1.'
            call warning ('readswap', messag)
         end if
      end if

! --- Special case: all meteo data in a single file; handled separately      
      if (swMetFilAll == 1) then
         call MeteoInOneFile (1, idum)
      end if
      
! --- copy content of key-file to log-file
      write (logf,14)  
 14   format('*',70('-'),'*',/,' Echo of input file:',/)
      call copfl2(swp,swpfilnam,logf,.true.)

      return
      end

! ----------------------------------------------------------------------
      subroutine rddrb ()
! ----------------------------------------------------------------------
!     Date               : July 2002
!     Purpose            : read input data for basic drainage    
! ----------------------------------------------------------------------
      use variables, only: swallo,swdtyp,flMacroPore,dramet,dra,ipos,nrlevs,logf,swdivd,l,zbotdr,owltab,drares,infres,qdrtab,basegw,     &
                           wetper,khtop,khbot,kvtop,kvbot,zintf,entres,geofac,drfil,pathdrain,numlay,cofani,swnrsrf,swscre,cofintfl,expintfl,     &
                           swdislay,swtopdislay,ztopdislay,ftopdislay,shape,SwTopnrsrf,swdivdinf,FacDpthInf,swliminf,nowltab
      implicit none
      include 'arrays.fi'

! ----------------------------------------------------------------------
! --- local
      integer   ifnd,i,getun2, swintfl
      real(8)   datowl(maowl),level(maowl),qdrain(25),gwl(25)
      character(len=80)  filnam
      character(len=200) messag
      logical   rdinqr,rdinar

! ----------------------------------------------------------------------

! --- open file with drainage data
      filnam = trim(pathdrain)//trim(drfil)//'.dra'
      dra = getun2 (10,90,2)
      call rdinit(dra,logf,filnam)

! --- method to establish drainage/infiltration fluxes
      call rdsinr ('dramet',1,3,dramet) 

      if (dramet.ne.3) nrlevs = 1

! --- division of drainage fluxes
      call rdsinr ('swdivd',0,1,swdivd)
      if (swdivd.eq.0) then
          write(messag,'(3a)')                                          &
     &    ' Variabel SWDIVD=0 in input file : ',trim(filnam),           &
     &    ' this is not recommended and may cause numerical instability'
          call warn ('rddrb',messag,logf,swscre)
      endif
      if (swdivd.eq.1 .and. dramet.eq.1) then
          write(messag,'(3a,i3)')                                       &
     &    ' Variabel SWDIVD=1 and DRAMET=1 in inputfile :',trim(filnam),&
     &    ' this is not allowed for drainage method (dramet)= ',dramet
          call fatalerr ('Rddrb',messag)
      endif
      if (swdivd .eq. 1) then
        if(rdinqr('swdivdinf')) then
          call rdsinr ('swdivdinf',0,1,swdivdinf)
            if (swdivdinf.eq.1 .and. dramet.ne.3) then
              write(messag,'(3a,i3)')                                   &
     &        ' Variabel SWDIVDINF=1 and DRAMET not 1 in inputfile :',  &
     &        trim(filnam),' this option is not allowed for drainage ', &
     &        ' method (dramet)= ',dramet
              call fatalerr ('Rddrb',messag)
            endif
          call rdsdor ('FacDpthInf',0.d0,1.d0,FacDpthInf)
        else
          swdivdinf = 0
        endif
        if (rdinar('cofani')) then
          call rdfdor ('cofani',1.d-4,1000.d0,cofani,maho,numlay)
        else
          call rdsdor ('cofani',1.d-4,1000.d0,cofani(1))
        endif
      endif


! --- input table of drainage flux as function of groundwater level
      if (dramet.eq.1) then  
        call rdsdor ('lm1',1.0d0,1000.d0,l(1))
        call rdador ('gwl',-10000.0d0,10.0d0,gwl,25,ifnd)
        call rdfdor ('qdrain',-100.d0,1000.0d0,qdrain,25,ifnd)

        l(1) = 100.0d0*l(1)
        do  i = 1,50
          qdrtab(i) = 0.0d0
        end do
        do i = 1,ifnd
          qdrtab(i*2-1) = abs(gwl(i))
          qdrtab(i*2) = qdrain(i)
        end do

!   - in case of rapid drainage due to macropore flow: read zbotdr
        if (flMacroPore) then
           call rdsinr ('swdtyp' ,1,2,swdtyp(1))
           call rdsdor ('zdrabas',-1000.0d0,0.0d0,zbotdr(1))
        endif

! --- input drainage formula of Hooghoudt or Ernst
      elseif (dramet.eq.2) then
   
! ---   read drain characteristics
        call rdsdor ('lm2',1.0d0,1000.0d0,l(1))
        l(1) = 100.0d0*l(1)
        call rdsdor ('shape',0.0d0,1.0d0,shape)
        call rdsdor ('wetper',0.0d0,1000.0d0,wetper(1))
        call rdsdor ('zbotdr',-1000.0d0,0.0d0,zbotdr(1))
        call rdsdor ('entres',0.0d0,1000.0d0,entres)

! ---   read profile characteristics
        call rdsinr ('ipos',1,5,ipos)
        call rdsdor ('basegw',-1.0d4,0.0d0,basegw)
        call rdsdor ('khtop',0.0d0,1000.0d0,khtop)

        if (ipos.ge.3) then
          call rdsdor ('khbot',0.0d0,1000.0d0,khbot)
          call rdsdor ('zintf',-1.0d4,0.0d0,zintf)
        endif
        if (ipos.ge.4) then
          call rdsdor ('kvtop',0.0d0,1000.0d0,kvtop)
          call rdsdor ('kvbot',0.0d0,1000.0d0,kvbot)  
        endif
        if (ipos.eq.5) then
          call rdsdor ('geofac',0.0d0,100.0d0,geofac)
        endif

! --- drainage and infiltration resistance
      elseif (dramet.eq.3) then
   
        call rdsinr ('nrlevs',1,Madr,nrlevs)
        if (nrlevs.gt.5) then
          write(messag,'(3a)')                                          &
     &    ' Number of Drainage levels >5 in inputfile : ',trim(filnam), &
     &    '   part of the output is limited to 5 levels'
          call warn ('rddrb',messag,logf,swscre)
        endif

! --- type of highest drainage level
        call rdsinr ('swintfl',0,1,swintfl)
!   - 'swintfl' is only used as input variabele, its function is taken over in the entire code by 'swintfl'    
        swnrsrf = swintfl
        if (swnrsrf.eq.1) then
          call rdsdor ('cofintflb',0.01d0,10.0d0,cofintfl)     
          call rdsdor ('expintflb',0.1d0,1.0d0,expintfl) 
        endif

! kroes 20080707 : allow disabling of option
        if (swdivd .eq. 1) then
          if (swnrsrf.gt.0)  call rdsinr ('SwTopnrsrf',0,1,SwTopnrsrf)
        endif

! ---   drainage level 1
        if (nrlevs.ge.1) then
          call rdsdor ('drares1',1.0d0,1.0d5,drares(1))
          call rdsdor ('infres1',0.0d0,1.0d5,infres(1))
          call rdsinr ('swallo1',1,3,swallo(1))
          if (swdivd .eq. 1) then
            call rdsdor ('l1',1.0d0,100000.0d0,l(1))
            l(1) = 100.0d0*l(1)
          endif
          call rdsdor ('zbotdr1',-10000.0d0,0.0d0,zbotdr(1))
          call rdsinr ('swdtyp1',1,2,swdtyp(1))
          ! always read surface water levels to allow adaptive drainage
          call rdatim ('datowl1',datowl,maowl,ifnd)
          call rdfdor ('level1',-10000.0d0,200.0d0,level,maowl,ifnd)
          ! store values in table
          do i = 1,ifnd
            owltab(1,i*2) = level(i)
            owltab(1,i*2-1) = datowl(i)
          end do
          nowltab(1) = ifnd
        endif 
! ---   drainage level 2
        if (nrlevs.ge.2) then
          call rdsdor ('drares2',1.0d0,1.0d5,drares(2))
          call rdsdor ('infres2',0.0d0,1.0d5,infres(2))
          call rdsinr ('swallo2',1,3,swallo(2))
          if (swdivd .eq. 1) then
            call rdsdor ('l2',1.0d0,100000.0d0,l(2))
            l(2) = 100.0d0*l(2)
          endif
          call rdsdor ('zbotdr2',-10000.0d0,0.0d0,zbotdr(2))
          call rdsinr ('swdtyp2',1,2,swdtyp(2))
          ! always read surface water levels to allow adaptive drainage
          call rdatim ('datowl2',datowl,maowl,ifnd)
          call rdfdor ('level2',-10000.0d0,200.0d0,level,maowl,ifnd)
          ! store values in table
          do i = 1,ifnd
            owltab(2,i*2) = level(i)
            owltab(2,i*2-1) = datowl(i)
          end do
          nowltab(2) = ifnd
        endif 
! ---   drainage level 3
        if (nrlevs.ge.3) then
          call rdsdor ('drares3',1.0d0,1.0d5,drares(3))
          call rdsdor ('infres3',0.0d0,1.0d5,infres(3))
          call rdsinr ('swallo3',1,3,swallo(3))
          if (swdivd .eq. 1) then
            call rdsdor ('l3',1.0d0,100000.0d0,l(3))
            l(3) = 100.0d0*l(3)
          endif
          call rdsdor ('zbotdr3',-10000.0d0,0.0d0,zbotdr(3))
          call rdsinr ('swdtyp3',1,2,swdtyp(3))
          ! always read surface water levels to allow adaptive drainage
          call rdatim ('datowl3',datowl,maowl,ifnd)
          call rdfdor ('level3',-10000.0d0,200.0d0,level,maowl,ifnd)
          ! store values in table
          do i = 1,ifnd
            owltab(3,i*2) = level(i)
            owltab(3,i*2-1) = datowl(i)
          end do
          nowltab(3) = ifnd
        endif 
! ---   drainage level 4
        if (nrlevs.ge.4) then
          call rdsdor ('drares4',1.0d0,1.0d5,drares(4))
          call rdsdor ('infres4',0.0d0,1.0d5,infres(4))
          call rdsinr ('swallo4',1,3,swallo(4))
          if (swdivd .eq. 1) then
            call rdsdor ('l4',1.0d0,100000.0d0,l(4))
            l(4) = 100.0d0*l(4)
          endif
          call rdsdor ('zbotdr4',-10000.0d0,0.0d0,zbotdr(4))
          call rdsinr ('swdtyp4',1,2,swdtyp(4))
          ! always read surface water levels to allow adaptive drainage
          call rdatim ('datowl4',datowl,maowl,ifnd)
          call rdfdor ('level4',-10000.0d0,200.0d0,level,maowl,ifnd)
          ! store values in table
          do i = 1,ifnd
            owltab(4,i*2) = level(i)
            owltab(4,i*2-1) = datowl(i)
          end do
          nowltab(4) = ifnd
        endif 
! ---   drainage level 5
        if (nrlevs.ge.5) then
          call rdsdor ('drares5',1.0d0,1.0d5,drares(5))
          call rdsdor ('infres5',0.0d0,1.0d5,infres(5))
          call rdsinr ('swallo5',1,3,swallo(5))
          if (swdivd .eq. 1) then
            call rdsdor ('l5',1.0d0,100000.0d0,l(5))
            l(5) = 100.0d0*l(5)
          endif
          call rdsdor ('zbotdr5',-10000.0d0,0.0d0,zbotdr(5))
          call rdsinr ('swdtyp5',1,2,swdtyp(5))
          ! always read surface water levels to allow adaptive drainage
          call rdatim ('datowl1',datowl,maowl,ifnd)
          call rdfdor ('level1',-10000.0d0,200.0d0,level,maowl,ifnd)
          ! store values in table
          do i = 1,ifnd
            owltab(5,i*2) = level(i)
            owltab(5,i*2-1) = datowl(i)
          end do
          nowltab(5) = ifnd
        endif

!       kro 20170906 : Switch for limit of infiltration head to the waterdepth in the channel
!                     (0 = nolimit, 1 = limitation)
        swliminf = 1
        if(rdinqr('swliminf')) then
          call rdsinr ('swliminf',0,1,swliminf)
        endif 

      endif 

!     top of model dicharge layer, determined by factor or direct input
      if (swdivd .eq. 1) then
        call rdsinr ('swdislay',0,2,swdislay)
        if (swdislay .eq. 1) then
          call rdfinr ('swtopdislay',0,1,swtopdislay,madr,nrlevs)
          call rdfdor('ztopdislay',-1.0d4,0.0d0,ztopdislay,madr,nrlevs)
        elseif (swdislay .eq. 2) then
          call rdfinr ('swtopdislay',0,1,swtopdislay,madr,nrlevs)
          call rdfdor('ftopdislay',0.0d0,1.0d0,ftopdislay,madr,nrlevs)
        endif
      endif


! --- close file with drainage data
      close (dra)         

      return
      end
      
! ----------------------------------------------------------------------
      subroutine readcropfixed (icrop,crpfil,lcc,swhydrlift)
! ----------------------------------------------------------------------
!     Update             : October 2012
!     Update             : July 2009
!     date               : July 2002             
!     purpose            : get crop parameters from cropfile
! ----------------------------------------------------------------------
      use variables, only: pathcrop,idev,tsumea,tsumam,tbase,kdif,kdir,gctb,swgc,cftb,cfeictb,swcf,                                                 &
                           hlim1,hlim2u,hlim2l,hlim3h,hlim3l,hlim4,rsc,adcrh,adcrl,logf,schedule,cumdens,chtb,                                      &
                           albedo,swetr,flsolute,swdrought,wiltpoint,flIrrigate,flIrg1Start,dayfix,rootradius,rootcoefa,rsw,                        &
                           q10_root,q10_microbial,specific_resp_humus,c_mroot,srl,dry_mat_cont_roots,air_filled_root_por,spec_weight_root_tissue,   &
                           var_a,f_senes,swhea,swcalt,swoxygen,swtopsub,nrstaring,oxygenslope,oxygenintercept,swoxygentype,rooteff,swsalinity,      &
                           criterhr,stephr,kroot,rxylem,taccur,kstem,wrtb,mrftb,swrootradius,root_radiusO2,bdens,inifil,                            &
                           nofd,daycrop,dvs,tsum,t1900,tstart,swinco,cropstart,saltmax,saltslope,salthead,                                          &
                           swjarvis,swcompensate,swstressor,alphacrit,dcritrtz,                                                                     &
                           swWrtNonox,aeratecrit,swinter,cofab,pfreetb,pstemtb,scanopytb,avprectb,avevaptb,fimin,siccaplai,sicact,dvsend,swharv,    &
                           swrd,rdtb,rdi,rri,rdc,swdmi2rd,swrdc,rdctb,rd,rdpot
      implicit  none
      include  'arrays.fi'
      
      integer              icrop,lcc,swhydrlift
      character(len=*)     crpfil
      
! --- locals
      integer              crp,getun2,ifnd,i,ini
      integer              swIrrigate
      real(8)              sum,depth,rootdis(202)
      real(8)              dvsinput(magrs),cfinput(magrs),chinput(magrs)
      real(8)              tinter(magrs),pfree(magrs),pstem(magrs)
      real(8)              scanopy(magrs),avprec(magrs),avevap(magrs)
      character(len=200)   message,filnam
      real(8)              afgen
      logical              rdinqr
! ----------------------------------------------------------------------

! --- initialise and start reading
      filnam = trim(pathcrop)//trim(crpfil)//'.crp'
      crp = getun2 (10,90,2)
      call rdinit(crp,logf,filnam)

! --- phenology
      call rdsinr ('idev',1,2,idev)
      if (idev.eq.1) then
        call rdsinr ('lcc',1,366,lcc)
      elseif (idev.eq.2) then
        call rdsdor ('tsumea',0.0d0,10000.0d0,tsumea)
        call rdsdor ('tsumam',0.0d0,10000.0d0,tsumam)
        call rdsdor ('tbase',-10.0d0, 30.0d0,tbase)
      endif

! --- assimilation                        
      call rdsdor ('kdif',0.0d0,2.0d0,kdif)
      call rdsdor ('kdir',0.0d0,2.0d0,kdir)
     
! --- LAI or soil cover fraction 
      call rdsinr ('swgc',1,2,swgc)
      if (swgc.eq.1) then
        call rdador ('gctb',0.0d0,12.0d0,gctb,(2*magrs),ifnd)
      elseif (swgc.eq.2) then
        call rdador ('gctb',0.0d0,2.0d0,gctb,(2*magrs),ifnd)
        do i = 2,ifnd,2
          if (gctb(i).gt.1.0d0) then
            message = 'Fatal Error in crop-(.crp)-file: SoilCover > 1 '
            call fatalerr ('ReadCropFixed',message)
          endif
        enddo
      endif

! --- Crop factor or crop height
      call rdsinr ('swcf',1,3,swcf)

! --- check use of crop factors in case of ETref
      if (swetr.eq.1 .and. swcf.eq.2) then
        message = 'If ETref is used (SWETR = 1), always define crop '// &
     &           'factors (SWCF = 1 or 3)' 
        call fatalerr ('ReadCropFixed',message)
      endif

      if (swcf.eq.1 .or. swcf.eq.3) then
! ---   crop factor is input
        call rdador ('dvs',0.0d0,2.0d0,dvsinput,(magrs),ifnd)
        call rdfdor ('cf',0.0d0,2.0d0,cfinput,(magrs),ifnd)
! ---   store values in cftb
        do i = 1,ifnd
          cftb(i*2) = cfinput(i) 
          cftb(i*2-1) = dvsinput(i)
        enddo
        if (swcf.eq.3) then
! ---     crop factor of wet crop is input
          call rdador ('dvs',0.0d0,2.0d0,dvsinput,(magrs),ifnd)
          call rdfdor ('cfw',0.0d0,2.0d0,cfinput,(magrs),ifnd)
! ---     store values in cfwtb with crop factor of wet crop is input
          do i = 1,ifnd
            cfeictb(i*2) = cfinput(i) 
            cfeictb(i*2-1) = dvsinput(i)
          enddo
        endif
        chtb = -99.99d0
      else
! ---   crop height is input
        call rdador ('dvs',0.0d0,2.0d0,dvsinput,(magrs),ifnd)
        call rdfdor ('ch',0.0d0,1.0d4,chinput,(magrs),ifnd)
! ---   store values in chtb
        do i = 1,ifnd
          chtb(i*2) = chinput(i) 
          chtb(i*2-1) = dvsinput(i)
        enddo
        cftb = -99.99d0
      endif

! --- reflection coefficient and crop resistance
      if (swcf.eq.1 .or. swcf.eq.3) then
! ---   use standard values for ETref
        albedo = 0.23d0
        rsc = 70.0d0
        rsw = 0.0d0
      else
! ---   use crop specific values
        call rdsdor ('albedo',0.0d0,1.0d0,albedo)
        call rdsdor ('rsc',0.0d0,1.0d6,rsc)
        call rdsdor ('rsw',0.0d0,1.0d6,rsw)
      endif

! --- rooting

      ! switch development root density
      swrdc = 0
      
      ! root density
      call rdador ('rdctb',0.0d0,100.0d0,rdctb,22,ifnd)
      
      ! switch development root extension
      swrd = 1
      if(rdinqr('swrd')) then
        call rdsinr ('swrd',1,3,swrd)
      endif
      
      ! root extension depends on development stage
      if (swrd.eq.1) then
        
        call rdador ('rdtb',0.0d0,1000.0d0,rdtb,(2*magrs),ifnd)
          
      ! root extension depends on maximum daily increase
      elseif (swrd.eq.2) then
        
        call rdsdor ('rdi',0.0d0,1000.0d0,rdi)
        call rdsdor ('rri',0.0d0,100.0d0,rri)
        call rdsdor ('rdc',0.0d0,1000.0d0,rdc)

        ! rooting depth influenced by dry matter increase (transpiration)
        swdmi2rd = 0
        if(rdinqr('swdmi2rd')) then
          call rdsinr ('swdmi2rd',0,1,swdmi2rd)
        endif
      
      ! root extension on available root biomass  
      elseif (swrd.eq.3) then
        
         message = 'Root extension based on available root biomass'//&
     &  ' is not possible with simple cropgrowth module.'
         call fatalerr ('readcropfixed',message)

      endif

!---  timing of harvest
      dvsend = 2.0d0
      if(rdinqr('dvsend')) then
        call rdsdor ('dvsend',0.0d0,3.0d0,dvsend)
      endif
      swharv  = 0
      if(rdinqr('swharv')) then
        call rdsinr ('swharv',0,1,swharv)
      endif
      
! --- oxygen stress
      swoxygen = 1
      if(rdinqr('swoxygen')) then
        call rdsinr ('swoxygen',0,2,swoxygen)
! -     fatal error when physical oxygen stress is simulated without numerical heat flow
        if (swoxygen.eq.2 .and. (swhea.eq.0 .or. swcalt.eq.1)) then
          message = 'In case oxygen stress is calculated according to'//&
     &   ' physical approach (SwOxygen = 2), soil heat flow should be'//&
     &   ' numerically simulated: SwCalT=2! Adapt .swp input file.'
          call fatalerr ('readcropfixed',message)
        endif
! -     fatal error when physical oxygen stress is simulated without realistic bdens value
        if (swoxygen.eq.2 .and. bdens(1).lt.100.d0) then
          message = 'In case oxygen stress is calculated according to'//&
     &   ' physical approach (SwOxygen = 2), bulk density must have'//  &
     &   ' realistic values; adjust BDENS-value(s) in .swp input file.'
          call fatalerr ('readcropfixed',message)
        endif
      endif

      if (swoxygen.eq.1) then                                          !
! ---   oxygen stress according to Feddes
        call rdsdor ('hlim1' ,-100.0d0,100.0d0,hlim1)                   !
        call rdsdor ('hlim2u',-1000.0d0,100.0d0,hlim2u)                 !
        call rdsdor ('hlim2l',-1000.0d0,100.0d0,hlim2l)                 !
      endif

      if (swoxygen.eq.2) then                                          
! ---   oxygen stress according to Bartholomeus

        swoxygentype = 1
        if(rdinqr('swoxygentype')) then
          call rdsinr ('swoxygentype',1,2,swoxygentype)
        endif

        if (swoxygentype.eq.1) then
! ---      use physical processes
           call rdsdor ('q10_root',1.0d0,4.0d0,q10_root)
           call rdsdor ('q10_microbial',1.0d0,4.0d0,q10_microbial)
           call rdsdor ('specific_resp_humus',0.0d0,1.0d0,              &
     &                   specific_resp_humus)
           call rdsdor ('c_mroot',0.0d0,1.0d0,c_mroot)
           call rdsdor ('srl',0.0d0,1.0d10,srl)
           call rdsdor ('f_senes',0.0d0,1.0d0,f_senes)         
           call rdador ('wrtb',0.0d0,10.0d0,wrtb,(2*magrs),ifnd)
           call rdador ('mrftb',0.0d0,100.0d0,mrftb,(2*magrs),ifnd)
           call rdsinr ('swrootradius',1,2,swrootradius)
           if (swrootradius .eq. 1) then
             call rdsdor ('dry_mat_cont_roots',0.0d0,1.0d0,             &
     &                   dry_mat_cont_roots)
             call rdsdor ('air_filled_root_por',0.0d0,1.0d0,            &
     &                   air_filled_root_por)
             call rdsdor ('spec_weight_root_tissue',0.0d0,1.0d5,        &
     &                   spec_weight_root_tissue)
             call rdsdor ('var_a',0.0d0,1.0d0,var_a)
           else
             call rdsdor ('root_radiusO2',1.0d-6,0.1d0,root_radiusO2)
           endif
        else
! ---      use reproduction functions
           call rdsinr ('SwTopSub',1,2,SwTopSub)
           call rdsinr ('NrStaring',1,18,NrStaring)
           call oxygen_dat (SwTopSub,NrStaring,OxygenSlope,             &
     &                      OxygenIntercept)
        endif
      endif

!
!     Check oxygen stress for root(zone) development
      swWrtNonox = 0
      if(rdinqr('swWrtNonox')) then
        call rdsinr ('swWrtNonox',0,1,swWrtNonox)
      endif
      
      aeratecrit = 0.0001d0
      if (swWrtNonox .eq. 1) then
        call rdsdor ('aeratecrit',1.d-4,1.d0,aeratecrit)
      endif
      
! --- drought stress
      swdrought = 1
      if(rdinqr('swdrought')) then
        call rdsinr ('swdrought',1,2,swdrought)                         
      endif

      if (swdrought.eq.1) then                                          
         
         ! drought stress according to Feddes
         call rdsdor ('hlim3h',-10000.0d0,100.0d0,hlim3h)               
         call rdsdor ('hlim3l',-10000.0d0,100.0d0,hlim3l)               
         call rdsdor ('hlim4' ,-20000.0d0,100.0d0,hlim4)                
         call rdsdor ('adcrh',0.0d0,5.0d0,adcrh)                        
         call rdsdor ('adcrl',0.0d0,5.0d0,adcrl)                        
      
      else
         ! drought stress according to De Jong van Lier
         call rdsdor ('wiltpoint',-1.0d8,-1.0d2,wiltpoint)
         call rdsdor ('kstem',1.0d-10,1.0d1,kstem)               
         call rdsdor ('rxylem',1.0d-4,1.d0,rxylem)               
         call rdsdor ('rootradius',1.0d-4,1.0d0,rootradius)            
         call rdsdor ('kroot',1.0d-10,1.0d10,kroot)            
         call rdsdor ('rootcoefa',0.0d0,1.0d0,rootcoefa)
         call rdsinr ('swhydrlift',0,1,swhydrlift)
         call rdsdor ('rooteff',0.0d0,1.0d0,rooteff)                 
         call rdsdor ('stephr',0.0d0,10.d0,stephr)                 
         call rdsdor ('criterhr',0.0d0,10.d0,criterhr)                 
         call rdsdor ('taccur',0.0d-5,10.d-2,taccur)                 
      endif                                                             


! --- salt stress
      if (flsolute) then

        if(rdinqr('swsalinity')) then
          call rdsinr ('swsalinity',0,2,swsalinity)
        endif
 
        if (swsalinity .eq. 1) then
          
          ! input for Maas and Hoffman salt reduction function
          call rdsdor ('saltmax',0.0d0,100.0d0,saltmax)
          call rdsdor ('saltslope',0.0d0,1.0d0,saltslope)
        elseif (swsalinity .eq. 2) then
          
          ! osmotic head salinity stress concept
          call rdsdor ('salthead',0.0d0,1.0d3,salthead)
          if (swdrought.eq.1) then
           message = 'In case salinity stress is calculated with'//     &
     &     ' osmotic head (SwSalinity = 2), the drought stress'//       &
     &     ' should be calculated according to De Jong van Lier:'//     &
     &     ' SwDrought = 2! Adapt .crp input file.'
           call fatalerr ('readcropfixed',message)
          endif
        endif
      endif

! --- compensation of root water uptake stress (-)
      swcompensate = 0
      swjarvis = 0
      if (rdinqr('swcompensate')) then
        call rdsinr ('swcompensate',0,2,swcompensate)
      else
        if (rdinqr('swjarvis')) then
          write(message,'(a)') 'SWJARVIS is deprecated, use SWCOMPENSATE instead'
          call warn ('readcropfixed',message,logf,0)
          call rdsinr ('swjarvis',0,4,swjarvis)
          if (swjarvis .gt. 0) then
            if (swjarvis .ne. 4) then
                write(message,'(a)') 'SWJARVIS is applied to all forms of stresses'
                call warn ('readcropfixed',message,logf,0)
            endif
            swcompensate = 1
          endif
        endif
      endif
      
      ! check compensation in combination with drought stress according to De Jong van Lier
      if (swcompensate .gt. 0 .and. swdrought .eq.2) then
           message = 'Compensation of root water uptake stress'//       &
     &     ' is not allowed in case drought stress is calculated'//     &
     &     ' according to De Jong van Lier. Adapt .crp input file.'
           call fatalerr ('readcropfixed',message)
      endif
      
      ! selection of stressors to compensate (default: all stressors)
      swstressor = 1
      if (swcompensate .gt. 0) then
        if (rdinqr('swstressor')) then
          call rdsinr ('swstressor',1,5,swstressor)
        endif
      endif

      ! compensated root water uptake according to Jarvis (1989)
      if (swcompensate .eq. 1) then
          
        ! criticial stress index for compensation of root water uptake (-)
        alphacrit = 1.0d0
        call rdsdor ('alphacrit',0.2d0,1.0d0,alphacrit)
      
      ! compensated root water uptake according to Walsum (2020)
      elseif (swcompensate .eq. 2) then
         
        ! threshold rootzone depth for compensation of root water uptake (cm)
        dcritrtz = 0.0d0
        call rdsdor ('dcritrtz',0.02d0,100.0d0,dcritrtz)

      endif

! --- interception
      call rdsinr ('swinter',0,3,swinter)
      if (swinter .eq. 1) then
        call rdsdor ('cofab',0.0d0,2.0d0,cofab)
      else if (swinter .eq. 2) then
        call rdador ('t',0.d0,366.d0,tinter,(magrs),ifnd)
        call rdfdor ('pfree',0.d0,1.d0,pfree,(magrs),ifnd)
        call rdfdor ('pstem',0.d0,1.d0,pstem,(magrs),ifnd)
        call rdfdor ('scanopy',0.d0,10.d0,scanopy,(magrs),ifnd)
        call rdfdor ('avprec',0.d0,100.d0,avprec,(magrs),ifnd)
        call rdfdor ('avevap',0.d0,10.d0,avevap,(magrs),ifnd)
        do i = 1, ifnd
          pfreetb(i*2) = pfree(i)
          pfreetb(i*2-1) = tinter(i)
          pstemtb(i*2) = pstem(i)
          pstemtb(i*2-1) = tinter(i)
          scanopytb(i*2) = scanopy(i)
          scanopytb(i*2-1) = tinter(i)
          avprectb(i*2) = avprec(i)
          avprectb(i*2-1) = tinter(i)
          avevaptb(i*2) = avevap(i)
          avevaptb(i*2-1) = tinter(i)
        end do
      else if (swinter .eq. 3) then
        call rdsdor ('fimin',0.0d0,1.0d0,fimin)
        call rdsdor ('siccaplai',0.0d0,1.0d-2,siccaplai)
      endif

! --- determine whether irrigation scheduling is applied
      call rdsinr ('schedule',0,1,schedule)

      if (schedule.eq.1 .and. swdrought.eq.2) then
! ---    read limiting pressure heads for irrigation scheduling
         call rdsdor ('hlim3h',-10000.0d0,100.0d0,hlim3h)               
         call rdsdor ('hlim3l',-10000.0d0,100.0d0,hlim3l)               
         call rdsdor ('hlim4' ,-20000.0d0,100.0d0,hlim4)                
      endif

! --- close file with crop data
      close (crp)

      ! irrigation
      if (schedule .eq. 1) then
          call irrigation(1)
          if (flIrg1Start) call IrrigationOutput(1)
      end if

! --- CALCULATE NORMALIZED CUMULATIVE ROOT DENSITY FUNCTION
      if (swdrought.eq.1) then 
! ---   root water extraction according to Feddes function  


! ---   specify array ROOTDIS with root density distribution
        do i = 0,100
          depth = 0.01d0 * dble(i)
          rootdis(i*2+1) = depth
          rootdis(i*2+2) = afgen(rdctb,22,depth)
        enddo

! ---   calculate cumulative root density function
        do i = 1,202,2
! ---     relative depths
          cumdens(i) = rootdis(i)
        enddo
        sum = 0.d0
        cumdens(2) = 0.d0
        do i = 4,202,2
! ---     cumulative root density
          sum = sum + (rootdis(i-2)+rootdis(i)) * 0.5d0                 &
     &               * (cumdens(i-1)-cumdens(i-3))
          cumdens(i) = sum
        enddo

! ---   normalize cumulative root density function to one
        do i = 2,202,2
          cumdens(i) = cumdens(i) / sum
        enddo
      endif

! --- read crop data of former day from *.END file
      if (t1900 - tstart .lt. 1.d-3 .and. swinco .eq. 3 .and. .not.     &
     &   abs(t1900 - cropstart(icrop)) .lt. 1.d-3) then

!       open inifile
        ini = getun2 (10,90,2)
        call rdinit(ini,logf,inifil)

        ! irrigation
        call rdsinr ('swIrrigate',0,1,swIrrigate)
        if (swIrrigate .eq. 0) then
          flIrrigate = .false.
        else
          flIrrigate = .true.
        endif
        if (flIrrigate) then
          call rdsinr ('dayfix',0,366,dayfix)
        endif
        
        ! crop development
        call rdsdor ('rd',0.0d0,1.0d4,rd)
        call rdsdor ('rdpot',0.0d0,1.0d4,rdpot)
        if (swinter.eq.3) then
          call rdsdor ('sicact',0.0d0,1.0d-2,sicact)
        endif
        
        call rdsinr ('nofd',0,7,nofd)
        call rdsinr ('daycrop',0,366,daycrop)
        call rdsdor ('dvs',0.0d0,2.0d0,dvs)
        call rdsdor ('tsum',0.0d0,1.0d4,tsum)

!       close inifile
        close(ini)
        
      endif

      return
      end

! ----------------------------------------------------------------------
      subroutine readwofost (icrop,crpfil,swhydrlift,swsoybean,mg,dvsi,dvrmax1,dvrmax2,flrfphotoveg, &
                             tmaxdvr,tmindvr,toptdvr,popt,pcrt,flphenodayl,fradeceasedlvtosoil)

! ----------------------------------------------------------------------
!     Update:  March 2017
!     Update:  December 2009
!     purpose: read parameters for wofost crop growth routine
! ----------------------------------------------------------------------
      use variables, only: pathcrop,swcf,cftb,idsl,dlo,dlc,tsumea,tsumam,dtsmtb,tdwi,laiem,rgrlai,slatb,spa,                                                      &
                           ssa,span,tbase,kdif,kdir,eff,amaxtb,tmpftb,tmnftb,cvl,cvo,cvr,cvs,q10,rml,rmo,rmr,rms,rfsetb,frtb,fltb,fstb,fotb,perdl,rdrrtb,         &
                           rdrstb,hlim1,hlim2u,hlim2l,hlim3h,hlim3l,hlim4,rsc,adcrh,adcrl,logf,schedule,cumdens,chtb,albedo,                                      &
                           swetr,flsolute,cfeictb,swdrought,wiltpoint,flIrrigate,flIrg1Start,dayfix,rootradius,rootcoefa,rsw,                                     &
                           q10_microbial,specific_resp_humus,srl,dry_mat_cont_roots,air_filled_root_por,                                                          &
                           spec_weight_root_tissue,var_a,swhea,swcalt,swoxygen,swtopsub,nrstaring,oxygenslope,                                                    &
                           oxygenintercept,swoxygentype,rooteff,swsalinity,criterhr,stephr,kroot,rxylem,taccur,kstem,                                             &
                           swrootradius,root_radiusO2,bdens,dvs,flanthesis,tsum,ilvold,ilvoldpot,wrt,wrtpot,tadw,                                                 &
                           tadwpot,wst,wstpot,wso,wsopot,wlv,wlvpot,laiexp,laiexppot,lai,laipot,laimax,dwrt,dwrtpot,dwlv,dwlvpot,dwst,dwstpot,dwlvSoil,dwlvCrop,  &
                           gasst,gasstpot,mrest,mrestpot,cwdm,cwdmpot,sla,slapot,lvage,lvagepot,lv,lvpot,inifil,daycrop,                                          &
                           dvsend,swharv,swwrtnonox,aeratecrit,nofd,atmin7,flCropCalendar,flCropEmergence,flCropHarvest,swinco,                                   &
                           t1900,tstart,cropstart,verndvs,vernsat,vernbase,vernrtb,swbulb,fbltb,pld,plwti,remoc,saltmax,saltslope,salthead,                       &
                           swjarvis,swcompensate,swstressor,alphacrit,dcritrtz,                                                                                   &
                           flco2,co2amaxtb,co2efftb,co2tratb,co2year,co2ppm,swpotrelmf,relmf,swinter,cofab,pfreetb,pstemtb,scanopytb,                             &
                           avprectb,avevaptb,fimin,siccaplai,sicact,swrd,rdtb,rdi,rri,rdc,swdmi2rd,rlwtb,wrtmax,swrdc,rdctb,rd,rdpot,glaiex,glaiexpot,tsumgerm,   &
                           sowdelay,prepdelay
       
      implicit none
      include  'arrays.fi'

      integer icrop,swhydrlift
      character(len=*)   crpfil

! --- only for soybean
      real(8) mg,dvsi,dvrmax1,dvrmax2,tmaxdvr,tmindvr,toptdvr
      integer swsoybean
      logical flrfphotoveg
      logical flphenodayl     ! Flag to allow input of POPT and PCRT or using empirical relation from Setiyono et al
      real(8) popt            ! optimal daylength for phenological development (hr)
      real(8) pcrt            ! critical daylength for phenological development (hr)

! --- CO2 impact
      integer   uco2
      character(len=200) atmofil
      integer   swco2

! --- Fraction  Deceased Leaves falling onto Soil 
      real(8)   FraDeceasedLvToSoil
      
! --- locals
      integer              crp,getun2,ifnd,i,ini,count
      integer              swCropHarvest,swCropEmergence,swAnthesis,swIrrigate
      real(8)              tinter(magrs),pfree(magrs),pstem(magrs)
      real(8)              scanopy(magrs),avprec(magrs),avevap(magrs)
      real(8)              dvsinput(magrs),cfinput(magrs),chinput(magrs)
      real(8)              laiinput(magrs),cfeicinput(magrs)
      real(8)              sum,afgen,depth,rootdis(202)
      character(len=200)   message,filnam
      logical              rdinqr

! ----------------------------------------------------------------------

! --- initialise and start reading
      filnam = trim(pathcrop)//trim(crpfil)//'.crp'
      crp = getun2 (10,90,2)
      call rdinit(crp,logf,filnam)

! --- crop factor or crop height
      call rdsinr ('swcf',1,3,swcf)

! --- check use of crop factors in case of ETref
      if (swetr.eq.1 .and. swcf.eq.2) then
        message = 'If ETref is used (SWETR = 1), always define crop '// &
     &           'factors (SWCF = 1 or 3)' 
        call fatalerr ('ReadWofost',message)
      endif

      if (swcf.eq.1) then
! ---   crop factor is input
        call rdador ('dvs',0.0d0,2.0d0,dvsinput,(magrs),ifnd)
        call rdfdor ('cf',0.0d0,2.0d0,cfinput,(magrs),ifnd)
! ---   store values in cftb
        do i = 1,ifnd
          cftb(i*2) = cfinput(i) 
          cftb(i*2-1) = dvsinput(i)
        enddo
        chtb = -99.99d0
      elseif (swcf.eq.2)then
! ---   crop height is input
        call rdador ('dvs',0.0d0,2.0d0,dvsinput,(magrs),ifnd)
        call rdfdor ('ch',0.0d0,1.0d4,chinput,(magrs),ifnd)
! ---   store values in chtb
        do i = 1,ifnd
          chtb(i*2) = chinput(i) 
          chtb(i*2-1) = dvsinput(i)
        enddo
        cftb = -99.99d0
      elseif (swcf.eq.3)then
! ---   LAI-dependent crop factors for dual crop coefficient
        call rdador ('lai',0.0d0,10.0d0,laiinput,36,ifnd)
        call rdfdor ('cf',0.0d0,2.0d0,cfinput,36,ifnd)
        call rdfdor ('cfeic',0.0d0,2.0d0,cfeicinput,36,ifnd)

        do i = 1,ifnd
          cftb(i*2)        = cfinput(i) 
          cftb(i*2-1)      = laiinput(i)
          cfeictb(i*2)     = cfeicinput(i) 
          cfeictb(i*2-1)   = laiinput(i)
        enddo

! ---   LAI-dependent crop height
        call rdfdor ('ch',0.0d0,1.0d2,chinput,36,ifnd)

        do i = 1,ifnd
          chtb(i*2)     = chinput(i)
          chtb(i*2-1)   = laiinput(i)
        enddo
      endif
      
! --- interception
      call rdsinr ('swinter',0,3,swinter)
      if (swinter .eq. 1) then
        call rdsdor ('cofab',0.0d0,2.0d0,cofab)
      else if (swinter .eq. 2) then
        call rdador ('t',0.d0,366.d0,tinter,(magrs),ifnd)
        call rdfdor ('pfree',0.d0,1.d0,pfree,(magrs),ifnd)
        call rdfdor ('pstem',0.d0,1.d0,pstem,(magrs),ifnd)
        call rdfdor ('scanopy',0.d0,10.d0,scanopy,(magrs),ifnd)
        call rdfdor ('avprec',0.d0,100.d0,avprec,(magrs),ifnd)
        call rdfdor ('avevap',0.d0,10.d0,avevap,(magrs),ifnd)
        do i = 1, ifnd
          pfreetb(i*2) = pfree(i)
          pfreetb(i*2-1) = tinter(i)
          pstemtb(i*2) = pstem(i)
          pstemtb(i*2-1) = tinter(i)
          scanopytb(i*2) = scanopy(i)
          scanopytb(i*2-1) = tinter(i)
          avprectb(i*2) = avprec(i)
          avprectb(i*2-1) = tinter(i)
          avevaptb(i*2) = avevap(i)
          avevaptb(i*2-1) = tinter(i)
        end do
      else if (swinter .eq. 3) then
        call rdsdor ('fimin',0.0d0,1.0d0,fimin)
        call rdsdor ('siccaplai',0.0d0,1.0d-2,siccaplai)
      endif

! --- reflection coefficient and crop resistance
      if (swcf.eq.1 .or. swcf.eq.3) then
! ---   use standard values for ETref
        albedo = 0.23d0
        rsc = 70.0d0
        rsw = 0.0d0
      else
! ---   use crop specific values
        call rdsdor ('albedo',0.0d0,1.0d0,albedo)
        call rdsdor ('rsc',0.0d0,1.0d6,rsc)
        call rdsdor ('rsw',0.0d0,1.0d6,rsw)
      endif

! --- phenology
! --- only for soybean:  maturity group parameters
      swsoybean = 0
      if(rdinqr('swsoybean')) then
        call rdsinr ('swsoybean',0,1,swsoybean)
      endif
      if (swsoybean.eq.1) then
        call rdsdor ('mg',0.1d0,9.0d0,mg)
        call rdsdor ('dvsi',0.0d0,2.0d0,dvsi)
        call rdsdor ('dvrmax1',0.0d0,1.0d0,dvrmax1)
        call rdsdor ('dvrmax2',0.0d0,1.0d0,dvrmax2)
        call rdsdor ('tmaxdvr',0.0d0,45.0d0,tmaxdvr)
        call rdsdor ('tmindvr',0.0d0,tmaxdvr,tmindvr)
        call rdsdor ('toptdvr',tmindvr,tmaxdvr,toptdvr)
        flrfphotoveg = .true.
        if(rdinqr('flrfphotoveg')) then
           call rdslog('flrfphotoveg',flrfphotoveg)
        endif
        flphenodayl = .false.
        if(rdinqr('flphenodayl')) then
           call rdslog('flphenodayl',flphenodayl)
           if(flphenodayl) then
              call rdsdor ('POPT',0.0d0,24.0d0,POPT)
              call rdsdor ('PCRT',0.0d0,24.0d0,PCRT)
           endif
        endif
      endif
! --- for non-soybean crops:
!     idsl =  Switch for crop development before anthesis: 
!             0 depends on temperature; 
!             1 depends on temperature and day length; 
!             2 depends on temperature, day length and vernalisation factor
      if (swsoybean.eq.0) then
        call rdsinr ('idsl',0,2,idsl)
        if (idsl.eq.1.or.idsl.eq.2) then
          call rdsdor ('dlo',0.0d0,24.0d0,dlo)
          call rdsdor ('dlc',0.0d0,24.0d0,dlc)
        endif
        call rdsdor ('tsumea',0.0d0,10000.0d0,tsumea)
        call rdsdor ('tsumam',0.0d0,10000.0d0,tsumam)
        call rdador ('dtsmtb',0.0d0,100.0d0,dtsmtb,30,ifnd)
      endif

!     vernalisation
      if (idsl.eq.2) then
        call rdsdor ('verndvs',0.0d0,0.4d0,verndvs)              ! critical development stage after which the effect of vernalisation is halted [-]
        call rdsdor ('vernsat',0.0d0,200.0d0,vernsat)            ! saturated vernalisation requirement [d]
        call rdsdor ('vernbase',0.0d0,200.0d0,vernbase)          ! base vernalisation requirement [d]
        call rdador ('vernrtb',-100.0d0,100.0d0,vernrtb,30,ifnd) ! rate of vernalisation as function of tav [days/degrees]
      endif

! --- only for bulb crops (tulips etc..)
      swbulb = 0
      if(rdinqr('swbulb')) then
        call rdsinr ('swbulb',0,1,swbulb)
      endif
      if(swbulb.eq.1) then
         call rdador ('fbltb',0.0d0,3.0d0,fbltb,30,ifnd)
!        parameters for remobilisation of carbohydrates from planted material
         call rdsdor ('plwti',0.0d0,10000.0d0,plwti)
         call rdsdor ('pld',0.0d0,10000000.0d0,pld)
         call rdsdor ('remoc',-1.0d0,0.0d0,remoc)
      endif

!---  timing of harvest
      call rdsdor ('dvsend',0.0d0,3.0d0,dvsend)
      swharv  = 0
      if(rdinqr('swharv')) then
        call rdsinr ('swharv',0,1,swharv)
      endif

! --- initial
      call rdsdor ('tdwi',0.0d0,10000.0d0,tdwi)
      call rdsdor ('laiem',0.0d0,10.0d0,laiem)
      call rdsdor ('rgrlai',0.0d0,1.0d0,rgrlai)

! --- green area
      call rdador ('slatb',0.0d0,2.0d0,slatb,30,ifnd)
      call rdsdor ('spa',0.0d0,1.0d0,spa)
      call rdsdor ('ssa',0.0d0,1.0d0,ssa)
      call rdsdor ('span',0.0d0,366.0d0,span)
      call rdsdor ('tbase',-10.0d0,30.0d0,tbase)

! --- assimilation
      call rdsdor ('kdif',0.0d0,2.0d0,kdif)
      call rdsdor ('kdir',0.0d0,2.0d0,kdir)
      call rdsdor ('eff',0.0d0,10.0d0,eff)
      call rdador ('amaxtb',0.0d0,100.0d0,amaxtb,30,ifnd)
      call rdador ('tmpftb',-10.0d0,50.0d0,tmpftb,30,ifnd)
      call rdador ('tmnftb',-10.0d0,50.0d0,tmnftb,30,ifnd)

! --- conversion of assimilates into biomass
      call rdsdor ('cvl',0.0d0,1.0d0,cvl)
      call rdsdor ('cvo',0.0d0,1.0d0,cvo)
      call rdsdor ('cvr',0.0d0,1.0d0,cvr)
      call rdsdor ('cvs',0.0d0,1.0d0,cvs)

! --- maintenance respiration
      call rdsdor ('q10',0.0d0,5.0d0,q10)
      call rdsdor ('rml',0.0d0,1.0d0,rml)
      call rdsdor ('rmo',0.0d0,1.0d0,rmo)
      call rdsdor ('rmr',0.0d0,1.0d0,rmr)
      call rdsdor ('rms',0.0d0,1.0d0,rms)
      call rdador ('rfsetb',0.0d0,3.0d0,rfsetb,30,ifnd)

! --- partitioning
      call rdador ('frtb',0.0d0,3.0d0,frtb,30,ifnd)
      call rdador ('fltb',0.0d0,3.0d0,fltb,30,ifnd)
      call rdador ('fstb',0.0d0,3.0d0,fstb,30,ifnd)
      call rdador ('fotb',0.0d0,3.0d0,fotb,30,ifnd)

! --- death rates
      call rdsdor ('perdl',0.0d0,3.0d0,perdl)
      call rdador ('rdrrtb',0.0d0,3.0d0,rdrrtb,30,ifnd)
      call rdador ('rdrstb',0.0d0,3.0d0,rdrstb,30,ifnd)

! --- oxygen stress
      swoxygen = 1
      if(rdinqr('swoxygen')) then
        call rdsinr ('swoxygen',0,2,swoxygen)
! -     fatal error when physical oxygen stress is simulated without numerical heat flow
        if (swoxygen.eq.2 .and. (swhea.eq.0 .or. swcalt.eq.1)) then
          message = 'In case oxygen stress is calculated according to'//&
     &   ' physical approach (SwOxygen = 2), soil heat flow should be'//&
     &   ' numerically simulated: SwCalT=2! Adapt .swp input file.'
          call fatalerr ('readwofost',message)
        endif
! -     fatal error when physical oxygen stress is simulated without realistic bdens value
        if (swoxygen.eq.2 .and. bdens(1).lt.100.0d0) then
          message = 'In case oxygen stress is calculated according to'//&
     &   ' physical approach (SwOxygen = 2), bulk density must have'//  &
     &   ' realistic values; adjust BDENS-value(s) in .swp input file.'
          call fatalerr ('readwofost',message)
        endif
      endif

      if (swoxygen.eq.1) then
! ---   oxygen stress according to Feddes
        call rdsdor ('hlim1' ,-100.0d0,100.0d0,hlim1)
        call rdsdor ('hlim2u',-1000.0d0,100.0d0,hlim2u)
        call rdsdor ('hlim2l',-1000.0d0,100.0d0,hlim2l)
      endif

      if (swoxygen.eq.2) then                                          
! ---   oxygen stress according to Bartholomeus

        swoxygentype = 1
        if(rdinqr('swoxygentype')) then
          call rdsinr ('swoxygentype',1,2,swoxygentype)
        endif

        if (swoxygentype.eq.1) then
! ---      use physical processes
           call rdsdor ('q10_microbial',1.0d0,4.0d0,q10_microbial)
           call rdsdor ('specific_resp_humus',0.0d0,1.0d0,              &
     &                   specific_resp_humus)                        
           call rdsdor ('srl',0.0d0,1.0d10,srl)                        
           call rdsinr ('swrootradius',1,2,swrootradius)
           if (swrootradius.eq.1) then
              call rdsdor ('dry_mat_cont_roots',0.0d0,1.0d0,            &
     &                   dry_mat_cont_roots)                        
              call rdsdor ('air_filled_root_por',0.0d0,1.0d0,           &
     &                   air_filled_root_por)                        
              call rdsdor ('spec_weight_root_tissue',0.0d0,1.0d5,       &
     &                   spec_weight_root_tissue)
              call rdsdor ('var_a',0.0d0,1.0d0,var_a)
           else 
              call rdsdor ('root_radiusO2',0.0d0,1.0d0,root_radiusO2)
           endif  
        else
! ---      use reproduction functions
           call rdsinr ('SwTopSub',1,2,SwTopSub)
           call rdsinr ('NrStaring',1,18,NrStaring)
           call oxygen_dat (SwTopSub,NrStaring,OxygenSlope,             &
     &                      OxygenIntercept)
        endif
      endif

!     Check oxygen stress for root(zone) development
      swWrtNonox = 0
      if(rdinqr('swWrtNonox')) then
        call rdsinr ('swWrtNonox',0,1,swWrtNonox)
      endif

      aeratecrit = 0.0001d0
      if (swWrtNonox .eq. 1) then
        call rdsdor ('aeratecrit',1.d-4,1.d0,aeratecrit)
      endif

! --- drought stress
      swdrought = 1
      if(rdinqr('swdrought')) then
        call rdsinr ('swdrought',1,2,swdrought)                         
      endif

      if (swdrought.eq.1) then                                          
         ! drought stress according to Feddes
         call rdsdor ('hlim3h',-10000.0d0,100.0d0,hlim3h)               
         call rdsdor ('hlim3l',-10000.0d0,100.0d0,hlim3l)               
         call rdsdor ('hlim4' ,-20000.0d0,100.0d0,hlim4)                
         call rdsdor ('adcrh',0.0d0,5.0d0,adcrh)                        
         call rdsdor ('adcrl',0.0d0,5.0d0,adcrl)                        
      
      else 

         ! drought stress according to De Jong van Lier
         call rdsdor ('wiltpoint',-1.0d8,-1.0d2,wiltpoint)
         call rdsdor ('kstem',1.0d-10,1.0d1,kstem)               
         call rdsdor ('rxylem',1.0d-4,1.d0,rxylem)               
         call rdsdor ('rootradius',1.0d-4,1.0d0,rootradius)            
         call rdsdor ('kroot',1.0d-10,1.0d10,kroot)            
         call rdsdor ('rootcoefa',0.0d0,1.0d0,rootcoefa)
         call rdsinr ('swhydrlift',0,1,swhydrlift)
         call rdsdor ('rooteff',0.0d0,1.0d0,rooteff)                 
         call rdsdor ('stephr',0.0d0,10.d0,stephr)                 
         call rdsdor ('criterhr',0.0d0,10.d0,criterhr)                 
         call rdsdor ('taccur',0.0d-5,10.d-2,taccur)                 
      endif

! --- salt stress
      if (flsolute) then

        if(rdinqr('swsalinity')) then
          call rdsinr ('swsalinity',0,2,swsalinity)
        endif
 
        if (swsalinity .eq. 1) then

          ! input for Maas and Hoffman salt reduction function
          call rdsdor ('saltmax',0.0d0,100.0d0,saltmax)
          call rdsdor ('saltslope',0.0d0,1.0d0,saltslope)
        
        elseif (swsalinity .eq. 2) then

          ! osmotic head salinity stress concept
          call rdsdor ('salthead',0.0d0,1.0d3,salthead)
          if (swdrought.eq.1) then
           message = 'In case salinity stress is calculated with'//     &
     &     ' osmotic head (SwSalinity = 2), the drought stress'//       &
     &     ' should be calculated according to De Jong van Lier:'//     &
     &     ' SwDrought = 2! Adapt .crp input file.'
           call fatalerr ('readcropfixed',message)
          endif
        endif
      endif

! --- compensation of root water uptake stress (-)
      swcompensate = 0
      swjarvis = 0
      if (rdinqr('swcompensate')) then
        call rdsinr ('swcompensate',0,2,swcompensate)
      else
        if (rdinqr('swjarvis')) then
          write(message,'(a)') 'SWJARVIS is deprecated, use SWCOMPENSATE instead'
          call warn ('readcropfixed',message,logf,0)
          call rdsinr ('swjarvis',0,4,swjarvis)
          if (swjarvis .gt. 0) then
            if (swjarvis .ne. 4) then
                write(message,'(a)') 'SWJARVIS is applied to all forms of stresses'
                call warn ('readcropfixed',message,logf,0)
            endif
            swcompensate = 1
          endif
        endif
      endif
      
      ! check compensation in combination with drought stress according to De Jong van Lier
      if (swcompensate .gt. 0 .and. swdrought .eq.2) then
           message = 'Compensation of root water uptake stress'//       &
     &     ' is not allowed in case drought stress is calculated'//     &
     &     ' according to De Jong van Lier. Adapt .crp input file.'
           call fatalerr ('readwofost',message)
      endif
      
      ! selection of stressors to compensate (default: all stressors)
      swstressor = 1
      if (swcompensate .gt. 0) then
        if (rdinqr('swstressor')) then
          call rdsinr ('swstressor',1,5,swstressor)
        endif
      endif

      ! compensated root water uptake according to Jarvis (1989)
      if (swcompensate .eq. 1) then
          
        ! criticial stress index for compensation of root water uptake (-)
        alphacrit = 1.0d0
        call rdsdor ('alphacrit',0.2d0,1.0d0,alphacrit)
      
      ! compensated root water uptake according to Walsum (2020)
      elseif (swcompensate .eq. 2) then
         
        ! threshold rootzone depth for compensation of root water uptake (cm)
        dcritrtz = 0.0d0
        call rdsdor ('dcritrtz',0.02d0,100.0d0,dcritrtz)

      endif

! --- management factor to account for other forms of stress
      relmf = 1.0d0
      if(rdinqr('relmf')) then
        call rdsdor ('relmf',0.0d0,1.0d0,relmf) 
      endif

! --- switch to account for difference between theoretical potential and attainable yield
      swpotrelmf = 1
      if(rdinqr('swpotrelmf')) then
        call rdsinr ('swpotrelmf',1,2,swpotrelmf)
      endif
      
! --- rooting
      
      ! switch development root density
      swrdc = 0
      if (rdinqr ('swrdc')) then
        call rdsinr ('swrdc',0,1,swrdc)
      endif
      
      ! root density
      call rdador ('rdctb',0.0d0,100.0d0,rdctb,22,ifnd)
      
      ! switch development root extension
      swrd = 2
      if (rdinqr ('swrd')) then
        call rdsinr ('swrd',1,3,swrd)
      endif
      
      ! root extension depends on development stage
      if (swrd.eq.1) then
        
        call rdador ('rdtb',0.0d0,1000.0d0,rdtb,(2*magrs),ifnd)
          
      ! root extension depends on maximum daily increase
      elseif (swrd.eq.2) then
        
        call rdsdor ('rdi',0.0d0,1000.0d0,rdi)
        call rdsdor ('rri',0.0d0,100.0d0,rri)
        call rdsdor ('rdc',0.0d0,1000.0d0,rdc)

        ! rooting depth influenced by dry matter increase
        swdmi2rd = 0
        if(rdinqr('swdmi2rd')) then
          call rdsinr ('swdmi2rd',0,1,swdmi2rd)
        endif
      
      ! root extension on available root biomass  
      elseif (swrd.eq.3) then
        
        call rdador ('rlwtb',0.0d0,5000.0d0,rlwtb,22,ifnd)
        call rdsdor ('wrtmax',0.0d0,100000.0d0,wrtmax)
      
      endif

! --- determine whether irrigation scheduling is applied
      call rdsinr ('schedule',0,1,schedule)

      if (schedule.eq.1 .and. swdrought.eq.2) then
! ---    read limiting pressure heads for irrigation scheduling
         call rdsdor ('hlim3h',-10000.0d0,100.0d0,hlim3h)               
         call rdsdor ('hlim3l',-10000.0d0,100.0d0,hlim3l)               
         call rdsdor ('hlim4' ,-20000.0d0,100.0d0,hlim4)                
      endif

!---  Fraction of deceased leaves entering the soil system
      FraDeceasedLvToSoil = 0.0d0
      if(rdinqr('FraDeceasedLvToSoil')) then
         call rdsdor ('FraDeceasedLvToSoil',0.d0,1.0d0,FraDeceasedLvToSoil) 
      endif
      
! --- assimilation correction due to CO2 impact      
      flco2 = .false.
      swco2 = 0
      if(rdinqr('swco2')) then
        call rdsinr ('swco2',0,1,swco2)
      endif
      if (swco2.eq.1) flco2 = .true.

!     Read CO2 data      
      if(flco2) then
         
         CALL rdadou ('CO2AMAXTB', CO2AMAXTB, 30, ifnd) 
         CALL rdadou ('CO2EFFTB', CO2EFFTB, 30, ifnd) 
         CALL rdadou ('CO2TRATB', CO2TRATB, 30, ifnd) 

         if(rdinqr('atmofil')) then
            call rdscha ('atmofil',atmofil)
            filnam = trim(pathcrop)//trim(atmofil)//'.co2'
         else
            filnam = trim(pathcrop)//'Atmospheric.co2'
         end if
         
         close (crp)
         
         uco2 = getun2 (30,90,2)
         call rdinit(uco2,logf,filnam)
         CALL rdainr ('CO2year', 1000, 3000, CO2year, mayrs, ifnd)
         CALL rdfdor ('CO2ppm',  10.0d0, 1000.0d0, CO2ppm, mayrs, ifnd)
         close(uco2)
      
      else
         
          close (crp)
      
      endif

      ! irrigation
      if (schedule .eq. 1) then
          call irrigation(1)
          if (flIrg1Start) call IrrigationOutput(1)
      end if
      
      
! --- CALCULATE NORMALIZED CUMULATIVE ROOT DENSITY FUNCTION

      if (swdrought.eq.1) then

! ---   specify array ROOTDIS with root density distribution
        do i = 0,100
          depth = 0.01d0 * dble(i)
          rootdis(i*2+1) = depth
          rootdis(i*2+2) = afgen(rdctb,22,depth)
        enddo

! ---   calculate cumulative root density function
        do i = 1,202,2
! ---     relative depths
          cumdens(i) = rootdis(i)
        enddo
        sum = 0.d0
        cumdens(2) = 0.d0
        do i = 4,202,2
! ---     cumulative root density
          sum = sum + (rootdis(i-2)+rootdis(i)) * 0.5d0                 &
     &               * (cumdens(i-1)-cumdens(i-3))
          cumdens(i) = sum
        enddo

! ---   normalize cumulative root density function to one
        do i = 2,202,2
          cumdens(i) = cumdens(i) / sum
        enddo
      endif 

! --- read crop data of former day from *.END file
      if (t1900 - tstart .lt. 1.d-3 .and. swinco .eq. 3 .and. .not.     &
     &   abs(t1900 - cropstart(icrop)) .lt. 1.d-3) then
         
        if (flCropCalendar) then
        
          ini = getun2 (10,90,2)
          call rdinit(ini,logf,inifil)
          
          if (rdinqr('swcropharvest')) then
          
            call rdsinr ('swcropharvest',0,1,swcropharvest)
            if (swcropharvest .eq. 0) then
              flCropHarvest = .false.
            else
              flCropHarvest = .true.
            endif
            
            call rdsinr ('swCropEmergence',0,1,swCropEmergence)
            if (swCropEmergence .eq. 0) then
              flCropEmergence = .false.
            else
              flCropEmergence = .true.
            endif
          
            call rdsinr ('PrepDelay',0,366,PrepDelay)
            call rdsinr ('SowDelay',0,366,SowDelay)
            call rdsdor ('tsumgerm',0.0d0,1.0d4,tsumgerm)
            
            if (.not. flCropHarvest) then
              
              ! irrigation
              call rdsinr ('swIrrigate',0,1,swIrrigate)
              if (swIrrigate .eq. 0) then
                flIrrigate = .false.
              else
                flIrrigate = .true.
              endif
              if (flIrrigate) then
                call rdsinr ('dayfix',0,366,dayfix)
              endif

              ! crop development
              call rdsdor ('rd',0.0d0,1.0d4,rd)
              call rdsdor ('rdpot',0.0d0,1.0d4,rdpot)
              if (swinter.eq.3) then
                call rdsdor ('sicact',0.0d0,1.0d-2,sicact)
              endif
              call rdsdor ('dvs',0.0d0,3.0d0,dvs)               
              call rdsinr ('daycrop',0,366,daycrop)
              call rdsinr ('swanthesis',0,1,swanthesis)               
              call rdsdor ('tsum',0.0d0,1.0d4,tsum)               
              call rdsinr ('ilvold',0,366,ilvold)               
              call rdsinr ('ilvoldpot',0,366,ilvoldpot)               
              call rdsdor ('wrt',0.0d0,1.0d8,wrt)               
              call rdsdor ('wrtpot',0.0d0,1.0d8,wrtpot)               
              call rdsdor ('tadw',0.0d0,1.0d8,tadw)               
              call rdsdor ('tadwpot',0.0d0,1.0d8,tadwpot)               
              call rdsdor ('wst',0.0d0,1.0d8,wst)               
              call rdsdor ('wstpot',0.0d0,1.0d8,wstpot)               
              call rdsdor ('wso',0.0d0,1.0d8,wso)               
              call rdsdor ('wsopot',0.0d0,1.0d8,wsopot)               
              call rdsdor ('wlv',0.0d0,1.0d8,wlv)               
              call rdsdor ('wlvpot',0.0d0,1.0d8,wlvpot)               
              call rdsdor ('laiexp',0.0d0,1.0d8,laiexp)
              call rdsdor ('laiexppot',0.0d0,1.0d8,laiexppot)
              call rdsdor ('glaiex',0.0d0,1.0d8,glaiex)
              call rdsdor ('glaiexpot',0.0d0,1.0d8,glaiexpot)
              call rdsdor ('lai',0.0d0,1.0d8,lai)
              call rdsdor ('laipot',0.0d0,1.0d8,laipot)
              call rdsdor ('laimax',0.0d0,1.0d8,laimax)
              call rdsdor ('dwrt',0.0d0,1.0d8,dwrt)
              call rdsdor ('dwrtpot',0.0d0,1.0d8,dwrtpot)
              call rdsdor ('dwlv',0.0d0,1.0d8,dwlv)
              call rdsdor ('dwlvpot',0.0d0,1.0d8,dwlvpot)
              call rdsdor ('dwst',0.0d0,1.0d8,dwst)
              call rdsdor ('dwstpot',0.0d0,1.0d8,dwstpot)
              call rdsdor ('dwlvSoil',0.0d0,1.0d8,dwlvSoil)
              call rdsdor ('dwlvCrop',0.0d0,1.0d8,dwlvCrop)
              call rdsdor ('gasst',0.0d0,1.0d10,gasst)
              call rdsdor ('gasstpot',0.0d0,1.0d10,gasstpot)
              call rdsdor ('mrest',0.0d0,1.0d10,mrest)
              call rdsdor ('mrestpot',0.0d0,1.0d10,mrestpot)
              call rdsdor ('cwdm',0.0d0,1.0d10,cwdm)
              call rdsdor ('cwdmpot',0.0d0,1.0d10,cwdmpot)
              call rdsinr ('nofd',0,366,nofd)               
              call rdfdor ('atmin7',-50.d0,60.d0,atmin7,7,7)
              
              count = max(ilvold,ilvoldpot)
              call rdfdor ('sla',0.d0,1.0d3,sla,366,count)
              call rdfdor ('slapot',0.d0,1.0d3,slapot,366,count)
              call rdfdor ('lvage',0.d0,1.0d3,lvage,366,count)
              call rdfdor ('lvagepot',0.d0,1.0d3,lvagepot,366,count)
              call rdfdor ('lv',0.d0,1.0d3,lv,366,count)
              call rdfdor ('lvpot',0.d0,1.0d3,lvpot,366,count)

              if (swanthesis .eq. 0) then
                flanthesis = .false.
              else
                flanthesis = .true.
              endif
            else

              daycrop = 0
              flCropEmergence = .false.

            endif
          else

            daycrop = 0
            flCropEmergence = .false.

          endif
        endif

!       close inifile
        close(ini)

      endif

      return
      end

    
! ----------------------------------------------------------------------
      subroutine readarablelandgerm (icrop,crpfil)

! ----------------------------------------------------------------------
!     Update:  December 2017
!     purpose: read parameters for crop germiniation
! ----------------------------------------------------------------------
      use variables, only: pathcrop,inifil,flCropPrep,zPrep,hPrep,MaxPrepDelay,PrepDelay,flCropSow,zSow,hSow,zTempSow,TempSow,MaxSowDelay,SowDelay,SwHea,&
                           flCropGerm,tsumemeopt,tbasem,teffmx,hdrygerm,hwetgerm,zgerm,agerm,cgerm,bgerm,tsumgerm,                                 &     
                           flCropEmergence,flCropHarvest,swinco,t1900,tstart,cropstart,daycrop,logf
      implicit none

      integer           icrop
      character(len=*)  crpfil
      
      integer swCropEmergence,swCropHarvest
      integer swCropPrep,swCropSow,swCropGerm
      integer swPrep,swSow,swGerm

      integer crp,ini,getun2
      character(len=200) message,filnam
      logical rdinqr

! ----------------------------------------------------------------------

! --- initialise and start reading
      filnam = trim(pathcrop)//trim(crpfil)//'.crp'
      crp = getun2 (10,90,2)
      call rdinit(crp,logf,filnam)

! --- Preparation before crop growth (default = 0)
      SwPrep = 0
      if (rdinqr('SwPrep')) then
        call rdsinr ('SwPrep',0,1,swPrep)
      endif

      flCropPrep = .true.
      if (SwPrep .eq. 1) then

        flCropPrep = .false.  
        PrepDelay  = 0

        call rdsdor ('zPrep',-1.d2,0.0d0,zPrep)
        call rdsdor ('hPrep',-2.d2,0.0d0,hPrep)
        call rdsinr ('MaxPrepDelay',1,366,MaxPrepDelay)

      endif

! --- Sowing before crop growth (default = 0)
      SwSow = 0
      if (rdinqr('SwSow')) then
        call rdsinr ('SwSow',0,1,swSow)
      endif
      
      flCropSow = .true.
      if (SwSow .eq. 1) then

        ! fatal error when sowing is simulated without heat flow
        if (SwHea.eq.0) then
          message = 'In case sowing is calculated, soil heat flow'//&
     &   ' should be simulated: SWHEA=1! Adapt .swp input file.'
          call fatalerr ('ReadArablelandGerm',message)
        endif

        flCropSow = .false.
        SowDelay  = 0

        call rdsdor ('zSow',-1.0d2,0.0d0,zSow)
        call rdsdor ('hSow',-2.0d2,0.0d0,hSow)
        call rdsdor ('zTempSow',-1.0d2,0.0d0,zTempSow)
        call rdsdor ('TempSow',0.0d0,30.0d0,TempSow)
        call rdsinr ('MaxSowDelay',1,366,MaxSowDelay)
      
      endif
      
!---  Simulation of germination (default = 0)
      SwGerm = 0
      if (rdinqr('SwGerm')) then
        call rdsinr ('SwGerm',0,2,swGerm)
      endif

      flCropGerm = .true.
      
      if (SwGerm.ge.1) then

        flCropGerm = .false.  

        call rdsdor ('tsumemeopt',0.0d0,1.0d3,tsumemeopt)
        call rdsdor ('tbasem',-20.0d0,4.0d1,tbasem)
        call rdsdor ('teffmx',0.0d0,4.0d1,teffmx)
      
        agerm = -99.d0
        
        if (SwGerm.eq.2) then

          call rdsdor ('hdrygerm',-1000.0d0,-1.0d-2,hdrygerm)
          call rdsdor ('hwetgerm',-100.0d0,-1.0d-2,hwetgerm)
          zgerm = -1.0d1
          if (rdinqr('zgerm')) then
            call rdsdor ('zgerm',-1.0d2,0.0d0,zGerm)
          endif
          call rdsdor ('agerm',1.0d0,1.0d3,agerm)
          if(rdinqr('bgerm')  .or.  rdinqr('cgerm')) then
            write(message,'(2a)')                                        &
          &    ' Variables bgerm and cgerm are determined by SWAP ',         &
          &    ' they are no longer input (see additional doc)'
            call warn ('ReadArablelandGerm',message,logf,0)
          endif
          cgerm = - (tsumemeopt - agerm*log10(-hdrygerm))
          bgerm =   (tsumemeopt + agerm*log10(-hwetgerm))
        endif
      endif
      
! --- Set crop emergence
      flCropEmergence = .false.
      if (flCropPrep .and. flCropSow .and. flCropGerm) then
        flCropEmergence = .true.
      endif
      
      close (crp)
      
! --- read crop data of former day from *.END file
      if (t1900 - tstart .lt. 1.d-3 .and. swinco .eq. 3 .and. .not.     &
     &  abs(t1900 - cropstart(icrop)) .lt. 1.d-3) then
         
        ini = getun2 (10,90,2)
        call rdinit(ini,logf,inifil)

        if (rdinqr('swcropharvest')) then
            
          call rdsinr ('swCropPrep',0,1,swCropPrep)
          if (swCropPrep .eq. 0) then
            flCropPrep = .false.
          else
            flCropPrep = .true.
          endif

          call rdsinr ('swCropSow',0,1,swCropSow)
          if (swCropSow .eq. 0) then
            flCropSow = .false.
          else
            flCropSow = .true.
          endif

          call rdsinr ('swCropGerm',0,1,swCropGerm)
          if (swCropGerm .eq. 0) then
            flCropGerm = .false.
          else
            flCropGerm = .true.
          endif
          
          call rdsinr ('swCropEmergence',0,1,swCropEmergence)
          if (swCropEmergence .eq. 0) then
            flCropEmergence = .false.
          else
            flCropEmergence = .true.
          endif

          call rdsinr ('swCropHarvest',0,1,swCropHarvest)
          if (swcropharvest .eq. 0) then
            flCropHarvest = .false.
          else
            flCropHarvest = .true.
          endif
          
          call rdsinr ('daycrop',0,366,daycrop)
          
          call rdsinr ('PrepDelay',0,366,PrepDelay)
          call rdsinr ('SowDelay',0,366,SowDelay)
          call rdsdor ('tsumgerm',0.0d0,1.0d4,tsumgerm)
          
        else

          daycrop         = 0
          tsumgerm        = 0.d0
          flCropPrep      = .false.
          PrepDelay       = 0
          flCropSow       = .false.
          SowDelay        = 0
          flCropEmergence = .false.
        
        endif
        
        close(ini)
      endif      
      
      return
      end

! ----------------------------------------------------------------------
      subroutine readgrass (icrop,crpfil,swharvest,dmharvest,daylastharvest,dmlastharvest,swdmmow,    &
                            maxdaymow,swlossmow,swlossgrz,swdmgrz, &
                            maxdaygrz,dmgrazing,lsdb,tagprest,swhydrlift)

! ----------------------------------------------------------------------
!     Update             : August 2014   
!     date               : november 2004   
!     purpose            : read parameters for grass growth routine
! ----------------------------------------------------------------------
      use variables, only: pathcrop,tdwi,laiem,rgrlai,slatb,ssa,span,tbase,kdif,kdir,eff,amaxtb,tmpftb,tmnftb,cvl,cvr,cvs,                                       &
                           q10,rml,rmr,rms,rfsetb,frtb,fltb,fstb,perdl,rdrrtb,rdrstb,hlim1,hlim2u,hlim2l,hlim3h,hlim3l,hlim4,rsc,adcrh,adcrl,                    &
                           logf,schedule,cumdens,flsolute,SeqGrazMow,dateharvest,dmmowtb,daygrowth,daygrowthpot,DelayRegrowthTab,                                &
                           lossmowtab,lossgrztab,dmgrztb,lsda,DaysGrazingtab,UptGrazingtab,LossGrazingtab,swcf,swetr,cftb,chtb,cfeictb,                          &
                           swdrought,wiltpoint,flIrrigate,flIrg1Start,dayfix,rootradius,rootcoefa,rsw,q10_microbial,specific_resp_humus,                         &
                           srl,dry_mat_cont_roots,air_filled_root_por,spec_weight_root_tissue,var_a,                                                             &
                           swhea,swcalt,swoxygen,swtopsub,nrstaring,oxygenslope,oxygenintercept,swoxygentype,rooteff,swsalinity,                                 &
                           criterhr,stephr,kroot,rxylem,taccur,kstem,swrootradius,root_radiusO2,bdens,swwrtnonox,aeratecrit,                                     &
                           dvs,flanthesis,tsum,ilvold,ilvoldpot,wrt,wrtpot,tadw,tadwpot,wst,wstpot,wso,wsopot,wlv,wlvpot,laiexp,laiexppot,laimax,lai,laipot,dwrt,&
                           dwrtpot,dwlv,dwlvpot,dwst,dwstpot,gasst,gasstpot,mrest,mrestpot,cwdm,cwdmpot,sla,slapot,lvage,lvagepot,lv,                            &
                           lvpot,inifil,daycrop,nofd,atmin7,rid,flgrazingpot,idregr,idregrpot,tagp,tagppot,mowrest,dewrest,                                      &
                           tagpt,tagptpot,iharvest,idaysgraz,iseqgm,iseqgmpot,flgrazing,idaysgrazpot,                                                            &
                           flCropCalendar,flCropEmergence,flCropHarvest,t1900,tstart,tend,swinco,cropstart,zmow,zgrz,                                            &
                           cuptgraz,cuptgrazpot,saltmax,saltslope,salthead,flco2,co2amaxtb,co2efftb,co2tratb,co2year,co2ppm,                                     &
                           swjarvis,swcompensate,swstressor,alphacrit,dcritrtz,                                                                                  &                 
                           albedo,swpotrelmf,relmf,swinter,cofab,pfreetb,pstemtb,scanopytb,avprectb,avevaptb,                                                    &
                           fimin,siccaplai,sicact,swrd,rdtb,rdi,rri,rdc,swdmi2rd,rlwtb,wrtmax,swrdc,rdctb,rd,rdpot,                                              &
                           swtsum,tsumtime,tsumtemp,tsumdepth,glaiex,glaiexpot
      implicit none
      include  'arrays.fi'

      integer           icrop,swhydrlift,swharvest,swdmmow,swlossmow,swlossgrz,swdmgrz
      integer           daylastharvest,maxdaymow,maxdaygrz
      real(8)           dmharvest,dmlastharvest,dmgrazing,lsdb(100),tagprest
      character(len=*)  crpfil
                            
! locals
      integer   crp,getun2,i,nrofSeqGM,ini,count,ifnd
      integer   swcropemergence,swcropharvest,swgrazing,swgrazingpot,swanthesis
      integer   swIrrigate
      integer   daydelay(100)
      real(8)   dnrinput(magrs),cfinput(magrs),chinput(magrs)
      real(8)   laiinput(magrs),cfeicinput(magrs),depth,rootdis(202)
      real(8)   tinter(magrs),pfree(magrs),pstem(magrs)
      real(8)   scanopy(magrs),avprec(magrs),avevap(magrs)
      real(8)   hlossmow(100),lossmow(100)
      real(8)   hlossgrz(100),lossgrz(100)
      real(8)   DaysGrazing(100),UptGrazing(100),LossGrazing(100),dmmowdelay(100)
      real(8)   afgen,sum
      logical   swMow,swGrz,swDew,RDinqr

      character(len=200) message,filnam
      
      ! CO2 impact
      integer   uco2
      character(len=200) atmofil
      integer   swco2
      
      
! ----------------------------------------------------------------------

! --- initialise and start reading
      filnam = trim(pathcrop)//trim(crpfil)//'.crp'
      crp = getun2 (10,90,2)
      call rdinit(crp,logf,filnam)

! --- ET related params  ---------

! --- crop factor or crop height
      call rdsinr ('swcf',1,3,swcf)

! --- check use of crop factors in case of ETref
      if (swetr.eq.1 .and. swcf.eq.2) then
        message = 'If ETref is used (SWETR = 1), always define crop '// &
     &           'factors (SWCF = 1 or 3)' 
        call fatalerr ('ReadGrass',message)
      endif

      if (swcf.eq.1) then
! ---   crop factor is input
        call rdador ('dnr',0.0d0,366.0d0,dnrinput,(magrs),ifnd)
        call rdfdor ('cf',0.0d0,2.0d0,cfinput,(magrs),ifnd)
! ---   store values in cftb
        do i = 1,ifnd
          cftb(i*2) = cfinput(i) 
          cftb(i*2-1) = dnrinput(i)
        enddo
        chtb = -99.99d0
      elseif (swcf.eq.2) then
! ---   crop height is input
        call rdador ('dnr',0.0d0,366.0d0,dnrinput,(magrs),ifnd)
        call rdfdor ('ch',0.0d0,1.0d4,chinput,(magrs),ifnd)
! ---   store values in chtb
        do i = 1,ifnd
          chtb(i*2) = chinput(i) 
          chtb(i*2-1) = dnrinput(i)
        enddo
        cftb = -99.99d0
      elseif (swcf.eq.3) then
! ---   LAI-dependent crop factors for dual crop coefficient
        call rdador ('lai',0.0d0,10.0d0,laiinput,36,ifnd)
        call rdfdor ('cf',0.0d0,2.0d0,cfinput,36,ifnd)
        call rdfdor ('cfeic',0.0d0,2.0d0,cfeicinput,36,ifnd)

        do i = 1,ifnd
          cftb(i*2)        = cfinput(i) 
          cftb(i*2-1)      = laiinput(i)
          cfeictb(i*2)     = cfeicinput(i) 
          cfeictb(i*2-1)   = laiinput(i)
        enddo

! ---   LAI-dependent crop height
        call rdfdor ('ch',0.0d0,1.0d2,chinput,36,ifnd)

        do i = 1,ifnd
          chtb(i*2)     = chinput(i)
          chtb(i*2-1)   = laiinput(i)
        enddo
      endif

! --- interception
      call rdsinr ('swinter',0,3,swinter)
      if (swinter .eq. 1) then
        call rdsdor ('cofab',0.0d0,2.0d0,cofab)
      else if (swinter .eq. 2) then
        call rdador ('t',0.d0,366.d0,tinter,(magrs),ifnd)
        call rdfdor ('pfree',0.d0,1.d0,pfree,(magrs),ifnd)
        call rdfdor ('pstem',0.d0,1.d0,pstem,(magrs),ifnd)
        call rdfdor ('scanopy',0.d0,10.d0,scanopy,(magrs),ifnd)
        call rdfdor ('avprec',0.d0,100.d0,avprec,(magrs),ifnd)
        call rdfdor ('avevap',0.d0,10.d0,avevap,(magrs),ifnd)
        do i = 1, ifnd
          pfreetb(i*2) = pfree(i)
          pfreetb(i*2-1) = tinter(i)
          pstemtb(i*2) = pstem(i)
          pstemtb(i*2-1) = tinter(i)
          scanopytb(i*2) = scanopy(i)
          scanopytb(i*2-1) = tinter(i)
          avprectb(i*2) = avprec(i)
          avprectb(i*2-1) = tinter(i)
          avevaptb(i*2) = avevap(i)
          avevaptb(i*2-1) = tinter(i)
        end do
      else if (swinter .eq. 3) then
        call rdsdor ('fimin',0.0d0,1.0d0,fimin)
        call rdsdor ('siccaplai',0.0d0,1.0d-2,siccaplai)
      endif
      
! --- reflection coefficient and crop resistance
      if (swcf.eq.1 .or. swcf.eq.3) then
! ---   use standard values for ETref
        albedo = 0.23d0
        rsc = 70.0d0
        rsw = 0.0d0
      else
! ---   use crop specific values
        call rdsdor ('albedo',0.0d0,1.0d0,albedo)
        call rdsdor ('rsc',0.0d0,1.0d6,rsc)
        call rdsdor ('rsw',0.0d0,1.0d6,rsw)
      endif


! --- crop growth related params ---------

! --- initial
      call rdsdor ('tdwi',0.0d0,10000.0d0,tdwi)
      call rdsdor ('laiem',0.0d0,10.0d0,laiem)
      call rdsdor ('rgrlai',0.0d0,1.0d0,rgrlai)

! --- yearly start of growth (0=no delay;1=TSUM200;2=TSOIL)
      call rdsinr ('swtsum',0,2,swtsum)
      if (swtsum.eq.2) then
        call rdsdor ('tsumtemp',0.0d0,20.0d0,tsumtemp)
        call rdsinr ('tsumtime',1,20,tsumtime)
        call rdsdor ('tsumdepth',0.0d0,100.0d0,tsumdepth)
      endif

! --- green area
      call rdador ('slatb',0.0d0,366.0d0,slatb,30,ifnd)
      call rdsdor ('ssa',0.0d0,1.0d0,ssa)
      call rdsdor ('span',0.0d0,366.0d0,span)
      call rdsdor ('tbase',-10.0d0,30.0d0,tbase)

! --- assimilation
      call rdsdor ('kdif',0.0d0,2.0d0,kdif)
      call rdsdor ('kdir',0.0d0,2.0d0,kdir)
      call rdsdor ('eff',0.0d0,10.0d0,eff)
      call rdador ('amaxtb',0.0d0,366.0d0,amaxtb,30,ifnd)
      call rdador ('tmpftb',-10.0d0,50.0d0,tmpftb,30,ifnd)
      call rdador ('tmnftb',-10.0d0,50.0d0,tmnftb,30,ifnd)

! --- conversion of assimilates into biomass
      call rdsdor ('cvl',0.0d0,1.0d0,cvl)
      call rdsdor ('cvr',0.0d0,1.0d0,cvr)
      call rdsdor ('cvs',0.0d0,1.0d0,cvs)

! --- maintenance respiration
      call rdsdor ('q10',0.0d0,5.0d0,q10)
      call rdsdor ('rml',0.0d0,1.0d0,rml)
      call rdsdor ('rmr',0.0d0,1.0d0,rmr)
      call rdsdor ('rms',0.0d0,1.0d0,rms)
      call rdador ('rfsetb',0.0d0,366.0d0,rfsetb,30,ifnd)

! --- partitioning
      call rdador ('frtb',0.0d0,366.0d0,frtb,30,ifnd)
      call rdador ('fltb',0.0d0,366.0d0,fltb,30,ifnd)
      call rdador ('fstb',0.0d0,366.0d0,fstb,30,ifnd)

! --- death rates
      call rdsdor ('perdl',0.0d0,3.0d0,perdl)
      call rdador ('rdrrtb',0.0d0,366.0d0,rdrrtb,30,ifnd)
      call rdador ('rdrstb',0.0d0,366.0d0,rdrstb,30,ifnd)

! --- oxygen stress
      swoxygen = 1
      if(rdinqr('swoxygen')) then
        call rdsinr ('swoxygen',0,2,swoxygen)
! -     fatal error when physical oxygen stress is simulated without numerical heat flow
        if (swoxygen.eq.2 .and. (swhea.eq.0 .or. swcalt.eq.1)) then
          message = 'In case oxygen stress is calculated according to'//&
     &   ' physical approach (SwOxygen = 2), soil heat flow should be'//&
     &   ' numerically simulated: SwCalT=2! Adapt .swp input file.'
          call fatalerr ('readgrass',message)
        endif
! -     fatal error when physical oxygen stress is simulated without realistic bdens value
        if (swoxygen.eq.2 .and. bdens(1).lt.100.d0) then
          message = 'In case oxygen stress is calculated according to'//&
     &   ' physical approach (SwOxygen = 2), bulk density must have'//  &
     &   ' realistic values; adjust BDENS-value(s) in .swp input file.'
          call fatalerr ('readgrass',message)
        endif
      endif

      if (swoxygen.eq.1) then                                          !
! ---   oxygen stress according to Feddes
        call rdsdor ('hlim1' ,-100.0d0,100.0d0,hlim1)                   !
        call rdsdor ('hlim2u',-1000.0d0,100.0d0,hlim2u)                 !
        call rdsdor ('hlim2l',-1000.0d0,100.0d0,hlim2l)                 !
      endif

      if (swoxygen.eq.2) then                                          
! ---   oxygen stress according to Bartholomeus

        swoxygentype = 1
        if(rdinqr('swoxygentype')) then
          call rdsinr ('swoxygentype',1,2,swoxygentype)
        endif

        if (swoxygentype.eq.1) then                                 
! ---      use physical processes
           call rdsdor ('q10_microbial',1.0d0,4.0d0,q10_microbial)  
           call rdsdor ('specific_resp_humus',0.0d0,1.0d0,              &
     &                   specific_resp_humus)                        
           call rdsdor ('srl',0.0d0,1.0d10,srl)                        
           call rdsinr ('swrootradius',1,2,swrootradius)
           if (swrootradius.eq.1) then
              call rdsdor ('dry_mat_cont_roots',0.0d0,1.0d0,            &
     &                   dry_mat_cont_roots)                        
              call rdsdor ('air_filled_root_por',0.0d0,1.0d0,           &
     &                   air_filled_root_por)                        
              call rdsdor ('spec_weight_root_tissue',0.0d0,1.0d5,       &
     &                   spec_weight_root_tissue)                   
              call rdsdor ('var_a',0.0d0,1.0d0,var_a)
           else
              call rdsdor ('root_radiusO2',0.0d0,1.0d0,root_radiusO2)
           endif  
        else
! ---      use reproduction functions
           call rdsinr ('SwTopSub',1,2,SwTopSub)
           call rdsinr ('NrStaring',1,18,NrStaring)
           call oxygen_dat (SwTopSub,NrStaring,OxygenSlope,             &
     &                      OxygenIntercept)
        endif
      endif

!     Growth of roots during oxygen stress 
      swWrtNonox = 0
      if(rdinqr('swWrtNonox')) then
        call rdsinr ('swWrtNonox',0,1,swWrtNonox)
      endif
      
      aeratecrit = 0.0001d0
      if (swWrtNonox .eq. 1) then
        call rdsdor ('aeratecrit',1.d-4,1.d0,aeratecrit)
      endif
      
! --- drought stress
      swdrought = 1
      if(rdinqr('swdrought')) then
        call rdsinr ('swdrought',1,2,swdrought)                         
      endif

      if (swdrought.eq.1) then                                          

         ! drought stress according to Feddes
         call rdsdor ('hlim3h',-10000.0d0,100.0d0,hlim3h)               
         call rdsdor ('hlim3l',-10000.0d0,100.0d0,hlim3l)               
         call rdsdor ('hlim4' ,-20000.0d0,100.0d0,hlim4)                
         call rdsdor ('adcrh',0.0d0,5.0d0,adcrh)                        
         call rdsdor ('adcrl',0.0d0,5.0d0,adcrl)                        
      
      else

         ! drought stress according to De Jong van Lier
         call rdsdor ('wiltpoint',-1.0d8,-1.0d2,wiltpoint)
         call rdsdor ('kstem',1.0d-10,1.0d1,kstem)               
         call rdsdor ('rxylem',1.0d-4,1.d0,rxylem)               
         call rdsdor ('rootradius',1.0d-4,1.0d0,rootradius)            
         call rdsdor ('kroot',1.0d-10,1.0d10,kroot)            
         call rdsdor ('rootcoefa',0.0d0,1.0d0,rootcoefa)
         call rdsinr ('swhydrlift',0,1,swhydrlift)
         call rdsdor ('rooteff',0.0d0,1.0d0,rooteff)                 
         call rdsdor ('stephr',0.0d0,10.d0,stephr)                 
         call rdsdor ('criterhr',0.0d0,10.d0,criterhr)                 
         call rdsdor ('taccur',0.0d-5,10.d-2,taccur)                 
      endif

! --- salt stress
      if (flsolute) then

        if(rdinqr('swsalinity')) then
          call rdsinr ('swsalinity',0,2,swsalinity)
        endif
 
        if (swsalinity .eq. 1) then
          ! input for Maas and Hoffman salt reduction function
          call rdsdor ('saltmax',0.0d0,100.0d0,saltmax)
          call rdsdor ('saltslope',0.0d0,1.0d0,saltslope)
        elseif (swsalinity .eq. 2) then
          ! osmotic head salinity stress concept
          call rdsdor ('salthead',0.0d0,1.0d3,salthead)
          if (swdrought.eq.1) then
           message = 'In case salinity stress is calculated with'//     &
     &     ' osmotic head (SwSalinity = 2), the drought stress'//       &
     &     ' should be calculated according to De Jong van Lier:'//     &
     &     ' SwDrought = 2! Adapt .crp input file.'
           call fatalerr ('readgrass',message)
          endif
        endif
      endif

! --- compensation of root water uptake stress (-)
      swcompensate = 0
      swjarvis = 0
      if (rdinqr('swcompensate')) then
        call rdsinr ('swcompensate',0,2,swcompensate)
      else
        if (rdinqr('swjarvis')) then
          write(message,'(a)') 'SWJARVIS is deprecated, use SWCOMPENSATE instead'
          call warn ('readcropfixed',message,logf,0)
          call rdsinr ('swjarvis',0,4,swjarvis)
          if (swjarvis .gt. 0) then
            if (swjarvis .ne. 4) then
                write(message,'(a)') 'SWJARVIS is applied to all forms of stresses'
                call warn ('readcropfixed',message,logf,0)
            endif
            swcompensate = 1
          endif
        endif
      endif
      
      ! check compensation in combination with drought stress according to De Jong van Lier
      if (swcompensate .gt. 0 .and. swdrought .eq.2) then
           message = 'Compensation of root water uptake stress'//       &
     &     ' is not allowed in case drought stress is calculated'//     &
     &     ' according to De Jong van Lier. Adapt .crp input file.'
           call fatalerr ('readgrass',message)
      endif
      
      ! selection of stressors to compensate (default: all stressors)
      swstressor = 1
      if (swcompensate .gt. 0) then
        if (rdinqr('swstressor')) then
          call rdsinr ('swstressor',1,5,swstressor)
        endif
      endif

      ! compensated root water uptake according to Jarvis (1989)
      if (swcompensate .eq. 1) then
          
        ! criticial stress index for compensation of root water uptake (-)
        alphacrit = 1.0d0
        call rdsdor ('alphacrit',0.2d0,1.0d0,alphacrit)
      
      ! compensated root water uptake according to Walsum (2020)
      elseif (swcompensate .eq. 2) then
         
        ! threshold rootzone depth for compensation of root water uptake (cm)
        dcritrtz = 0.0d0
        call rdsdor ('dcritrtz',0.02d0,100.0d0,dcritrtz)

      endif
      
! --- rooting
      
      ! switch development root density
      swrdc = 0
      if (rdinqr ('swrdc')) then
        call rdsinr ('swrdc',0,1,swrdc)
      endif
      
      ! root density
      call rdador ('rdctb',0.0d0,100.0d0,rdctb,22,ifnd)

      ! switch development root extension
      swrd = 3
      if(rdinqr('swrd')) then
        call rdsinr ('swrd',1,3,swrd)
      endif
      
      ! root extension depends on development stage
      if (swrd.eq.1) then
        
        call rdador ('rdtb',0.0d0,1000.0d0,rdtb,(2*magrs),ifnd)
          
      ! root extension depends on maximum daily increase
      elseif (swrd.eq.2) then
        
        call rdsdor ('rdi',0.0d0,1000.0d0,rdi)
        call rdsdor ('rri',0.0d0,100.0d0,rri)
        call rdsdor ('rdc',0.0d0,1000.0d0,rdc)

        ! rooting depth influenced by dry matter increase
        swdmi2rd = 0
        if(rdinqr('swdmi2rd')) then
          call rdsinr ('swdmi2rd',0,1,swdmi2rd)
        endif
      
      ! root extension on available root biomass  
      elseif (swrd.eq.3) then
        
        call rdador ('rlwtb',0.0d0,5000.0d0,rlwtb,22,ifnd)
        call rdsdor ('wrtmax',0.0d0,100000.0d0,wrtmax)
      
      endif

! --- management factor to account for other forms of stress (pests, diseases, nutrients, etc ..)
      relmf = 1.0d0
      if(rdinqr('relmf')) then
        call rdsdor ('relmf',0.0d0,1.0d0,relmf) 
      endif

! --- switch to account for difference between theoretical potential and attainable yield
      swpotrelmf = 1
      if(rdinqr('swpotrelmf')) then
        call rdsinr ('swpotrelmf',1,2,swpotrelmf)
      endif

! --- Management of mowing and grazing
! -   0. Define periods with mowing or grazing 

!     read (yearly) sequence of periods with mowing-grazing (1=grazing, 2=mowing, 3=grazing with dewooling)
      call rdainr ('SeqGrazMow',1,3,SeqGrazMow,366,nrofSeqGM)

      swMow = .false.      
      swGrz = .false.
      swDew = .false.
      do i = 1,nrofSeqGM
         if (SeqGrazMow(i) .eq. 1) then
            swGrz = .true.
         elseif (SeqGrazMow(i) .eq. 2) then
            swMow = .true. 
         elseif (SeqGrazMow(i) .eq. 3) then
            swGrz = .true.
            swDew = .true.
         endif
      end do
      
! -   1. Grazing settings

      ! trigger of harvest
      if (swGrz) then
         call rdsinr ('swharvest',1,2,swharvest)
         if(swharvest .eq. 1) then
            
           ! trigger mowing by fixed or flexible dry matter threshold
           call rdsinr ('swdmgrz',1,2,swdmgrz)
           if (swdmgrz .eq. 1) then
             call rdsdor ('dmgrazing',0.0d0,100000.0d0,dmgrazing)
           elseif (swdmgrz .eq. 2) then    
             call rdador ('dmgrztb',0.0d0,100000.0d0,dmgrztb,20,ifnd)
             call rdsinr ('maxdaygrz',1,366,maxdaygrz)
           endif
             
           ! Losses due to insufficient pressure head during grazing
           call rdsinr ('swlossgrz',0,1,swlossgrz)
           if (swlossgrz .eq. 1) then
             call rdsdor ('zgrz',-1.0d2,0.0d0,zgrz)
             call rdador ('hlossgrz',-1000.0d0,0.0d0,hlossgrz,100,ifnd)
             call rdador ('lossgrz',0.d0,1.d0,lossgrz,100,ifnd)
           
             ! store values in lossgrztab
             do i = 1,ifnd
               lossgrztab(i*2-1) = hlossgrz(i)
               lossgrztab(i*2)   = lossgrz(i)
             enddo

           endif
           
         elseif(swharvest.eq.2) then
           
           ! trigger mowing by date
           call rdatim ('dateharvest',dateharvest,999,ifnd)
           dateharvest(ifnd + 1) = tend + 1.d0
         
         endif

         if (swDew) then
            call rdsdor ('dewrest',0.0d0,100000.0d0,dewrest)
         endif
        
         call rdsdor ('tagprest',0.0d0,100000.0d0,tagprest)

!        LiveStock Density basal with days and uptake during grazing
         call rdador ('LSDb',0.0d0,1000.0d0,LSDb,100,ifnd)
         call rdfdor ('daysgrazing',0.0d0,366.0d0,daysgrazing,100,ifnd)
         call rdfdor ('uptgrazing',0.0d0,1000.0d0,uptgrazing,100,ifnd)
         call rdfdor ('lossgrazing',0.0d0,1000.0d0,lossgrazing,100,ifnd)

!        store values in daysgrazingtab, uptgrazingtab and lossgrazingtab 
         do i = 1,ifnd
           daysgrazingtab(i*2)   = daysgrazing(i)
           daysgrazingtab(i*2-1) = LSDb(i)
         enddo
         do i = 1,ifnd
           uptgrazingtab(i*2)   = uptgrazing(i)
           uptgrazingtab(i*2-1) = LSDb(i)
         enddo
         do i = 1,ifnd
           lossgrazingtab(i*2)   = lossgrazing(i)
           lossgrazingtab(i*2-1) = LSDb(i)
         enddo

!        LiveStock Density actual
         call rdador ('LSDa',0.0d0,1000.0d0,LSDa,366,ifnd)

!        verify length of array, should correspond to periods with grazing
         if (ifnd.ne.nrofSeqGM) then
          message = 'Dynamic grassland, input for LSDa'//               &
     &   ' (LiveStock Density actual); nr of values must be equal'//    &
     &   ' to the nr of values in SEQGRAZMOW! Adapt .crp input file.'
         endif
      
      endif
      
! -   2. Mowing settings
      if (swMow) then
        
        ! Remaining yield after mowing
        call rdsdor ('mowrest',0.0d0,5000.0d0,mowrest)
        
        ! trigger of harvest
        call rdsinr ('swharvest',1,2,swharvest)
        if(swharvest.eq.1) then

          ! trigger mowing by fixed or flexible dry matter threshold
          call rdsinr ('swdmmow',1,2,swdmmow)
          if (swdmmow .eq. 1) then
            call rdsdor ('dmharvest',0.0d0,100000.0d0,dmharvest)
            call rdsinr ('daylastharvest',1,366,daylastharvest)
            call rdsdor ('dmlastharvest',0.0d0,100000.0d0,dmlastharvest)
          elseif (swdmmow .eq. 2) then    
            call rdador ('dmmowtb',0.0d0,100000.0d0,dmmowtb,20,ifnd)
            call rdsinr ('maxdaymow',1,366,maxdaymow)
          endif

          ! Losses due to insufficient pressure head during mowing
          call rdsinr ('swlossmow',0,1,swlossmow)
          if (swlossmow .eq. 1) then
            call rdsdor ('zmow',-1.0d2,0.0d0,zmow)
            call rdador ('hlossmow',-1000.0d0,0.0d0,hlossmow,100,ifnd)
            call rdador ('lossmow',0.d0,1.d0,lossmow,100,ifnd)
          
            ! store values in LossMowTab
            do i = 1,ifnd
              lossmowtab(i*2-1) = hlossmow(i)
              lossmowtab(i*2)   = lossmow(i)
            enddo
          endif

        elseif(swharvest.eq.2) then

           ! trigger mowing by date
           call rdatim ('dateharvest',dateharvest,999,ifnd)
           dateharvest(ifnd + 1) = tend + 1.d0

        endif
        
        ! days delay of regrowth after mowing
        call rdainr ('daydelay',0,366,daydelay,100,ifnd)
        call rdador ('dmmowdelay',0.0d0,100000.0d0,dmmowdelay,100,ifnd)
        
        ! store values in DelayRegrowthTab
        do i = 1,ifnd
          DelayRegrowthTab(i*2)   = dble(daydelay(i))
          DelayRegrowthTab(i*2-1) = dmmowdelay(i)
        enddo
      
      endif

! --- determine whether irrigation scheduling is applied
      call rdsinr ('schedule',0,1,schedule)

      if (schedule.eq.1 .and. swdrought.eq.2) then
! ---    read limiting pressure heads for irrigation scheduling
         call rdsdor ('hlim3h',-10000.0d0,100.0d0,hlim3h)               
         call rdsdor ('hlim3l',-10000.0d0,100.0d0,hlim3l)               
         call rdsdor ('hlim4' ,-20000.0d0,100.0d0,hlim4)                
      endif

! --- assimilation correction due to CO2 impact      
      flco2 = .false.
      swco2 = 0
      if(rdinqr('swco2')) then
        call rdsinr ('swco2',0,1,swco2)
      endif
      if (swco2.eq.1) flco2 = .true.

!     Read CO2 data      
      if(flco2) then
         
         CALL rdadou ('CO2AMAXTB', CO2AMAXTB, 30, ifnd) 
         CALL rdadou ('CO2EFFTB', CO2EFFTB, 30, ifnd) 
         CALL rdadou ('CO2TRATB', CO2TRATB, 30, ifnd) 

         if(rdinqr('atmofil')) then
            call rdscha ('atmofil',atmofil)
            filnam = trim(pathcrop)//trim(atmofil)//'.co2'
         else
            filnam = trim(pathcrop)//'Atmospheric.co2'
         end if
         
         close (crp)
         
         uco2 = getun2 (30,90,2)
         call rdinit(uco2,logf,filnam)
         CALL rdainr ('CO2year', 1000, 3000, CO2year, mayrs, ifnd)
         CALL rdfdor ('CO2ppm',  10.0d0, 1000.0d0, CO2ppm, mayrs, ifnd)
         close(uco2)
      
      else
         
          close (crp)
      
      endif

      ! irrigation
      if (schedule .eq. 1) then
          call irrigation(1)
          if (flIrg1Start) call IrrigationOutput(1)
      end if
      
! --- CALCULATE NORMALIZED CUMULATIVE ROOT DENSITY FUNCTION

      if (swdrought.eq.1) then 

! ---   specify array ROOTDIS with root density distribution
        do i = 0,100
          depth = 0.01d0 * dble(i)
          rootdis(i*2+1) = depth
          rootdis(i*2+2) = afgen(rdctb,22,depth)
        enddo

! ---   calculate cumulative root density function
        do i = 1,202,2
! ---     relative depths
          cumdens(i) = rootdis(i)
        enddo
        sum = 0.d0
        cumdens(2) = 0.d0
        do i = 4,202,2
! ---     cumulative root density
          sum = sum + (rootdis(i-2)+rootdis(i)) * 0.5d0                 &
     &               * (cumdens(i-1)-cumdens(i-3))
          cumdens(i) = sum
        enddo

! ---   normalize cumulative root density function to one
        do i = 2,202,2
          cumdens(i) = cumdens(i) / sum
        enddo

      endif

! --- read crop data of former day from *.END file

      if (t1900 - tstart .lt. 1.d-3 .and. swinco .eq. 3 .and. .not.     &
     &   abs(t1900 - cropstart(icrop)) .lt. 1.d-3) then

        if (flCropCalendar) then  
          
          ini = getun2 (10,90,2)
          call rdinit(ini,logf,inifil)

          if (rdinqr('swcropharvest')) then
          
            call rdsinr ('swcropharvest',0,1,swcropharvest)
            if (swcropharvest .eq. 0) then
              flCropHarvest = .false.
            else
              flCropHarvest = .true.
            endif
            
            call rdsinr ('swCropEmergence',0,1,swCropEmergence)
            if (swCropEmergence .eq. 0) then
              flCropEmergence = .false.
            else
              flCropEmergence = .true.
            endif
          
            if (.not. flCropHarvest) then          
            
              ! irrigation
              call rdsinr ('swIrrigate',0,1,swIrrigate)
              if (swIrrigate .eq. 0) then
                flIrrigate = .false.
              else
                flIrrigate = .true.
              endif
              if (flIrrigate) then
                call rdsinr ('dayfix',0,366,dayfix)
              endif
                
              ! crop development
              call rdsdor ('rd',0.0d0,1.0d4,rd)               
              call rdsdor ('rdpot',0.0d0,1.0d4,rdpot)               
              if (swinter.eq.3) then
                call rdsdor ('sicact',0.0d0,1.0d-2,sicact)
              endif
              call rdsdor ('dvs',0.0d0,3.0d0,dvs)               
              call rdsinr ('daycrop',0,366,daycrop)               
              call rdsinr ('swanthesis',0,1,swanthesis)               
              call rdsdor ('tsum',0.0d0,1.0d4,tsum)               
              call rdsinr ('ilvold',0,366,ilvold)               
              call rdsinr ('ilvoldpot',0,366,ilvoldpot)               
              call rdsdor ('wrt',0.0d0,1.0d8,wrt)               
              call rdsdor ('wrtpot',0.0d0,1.0d8,wrtpot)               
              call rdsdor ('tadw',0.0d0,1.0d8,tadw)               
              call rdsdor ('tadwpot',0.0d0,1.0d8,tadwpot)               
              call rdsdor ('wst',0.0d0,1.0d8,wst)               
              call rdsdor ('wstpot',0.0d0,1.0d8,wstpot)               
              call rdsdor ('wso',0.0d0,1.0d8,wso)               
              call rdsdor ('wsopot',0.0d0,1.0d8,wsopot)               
              call rdsdor ('wlv',0.0d0,1.0d8,wlv)               
              call rdsdor ('wlvpot',0.0d0,1.0d8,wlvpot)               
              call rdsdor ('laiexp',0.0d0,1.0d8,laiexp)
              call rdsdor ('laiexppot',0.0d0,1.0d8,laiexppot)
              call rdsdor ('glaiex',0.0d0,1.0d8,glaiex)
              call rdsdor ('glaiexpot',0.0d0,1.0d8,glaiexpot)
              call rdsdor ('lai',0.0d0,1.0d8,lai)
              call rdsdor ('laipot',0.0d0,1.0d8,laipot)
              call rdsdor ('dwrt',0.0d0,1.0d8,dwrt)
              call rdsdor ('dwrtpot',0.0d0,1.0d8,dwrtpot)
              call rdsdor ('dwlv',0.0d0,1.0d8,dwlv)
              call rdsdor ('dwlvpot',0.0d0,1.0d8,dwlvpot)
              call rdsdor ('dwst',0.0d0,1.0d8,dwst)
              call rdsdor ('dwstpot',0.0d0,1.0d8,dwstpot)
              call rdsdor ('gasst',0.0d0,1.0d10,gasst)
              call rdsdor ('gasstpot',0.0d0,1.0d10,gasstpot)
              call rdsdor ('mrest',0.0d0,1.0d10,mrest)
              call rdsdor ('mrestpot',0.0d0,1.0d10,mrestpot)
              call rdsdor ('cwdm',0.0d0,1.0d10,cwdm)
              call rdsdor ('cwdmpot',0.0d0,1.0d10,cwdmpot)
              call rdsinr ('nofd',0,366,nofd)               
              call rdfdor ('atmin7',-50.d0,60.d0,atmin7,7,7)
              call rdsdor ('rid',0.0d0,1.0d4,rid)
              call rdsinr ('idregr',0,366,idregr)               
              call rdsinr ('idregrpot',0,366,idregrpot)               
              call rdsdor ('laimax',0.0d0,1.0d2,laimax)
              call rdsdor ('tagp',0.0d0,1.0d5,tagp)
              call rdsdor ('tagppot',0.0d0,1.0d5,tagppot)
              call rdsdor ('tagpt',0.0d0,1.0d5,tagpt)
              call rdsdor ('tagptpot',0.0d0,1.0d5,tagptpot)
              call rdsinr ('daygrowth',0,366,daygrowth)
              call rdsinr ('daygrowthpot',0,366,daygrowthpot)
              call rdsdor ('cuptgraz',0.0d0,1.0d5,cuptgraz)
              call rdsdor ('cuptgrazpot',0.0d0,1.0d5,cuptgrazpot)
              call rdsdor ('tsum',0.0d0,1.0d4,tsum)
              call rdsinr ('swcropharvest',0,1,swcropharvest)               
              call rdsinr ('iseqgm',0,366,iseqgm)               
              call rdsinr ('iseqgmpot',0,366,iseqgmpot)               
              call rdsinr ('swgrazing',0,1,swgrazing)
              call rdsinr ('swgrazingpot',0,1,swgrazingpot)
              call rdsinr ('iharvest',0,366,iharvest)
              call rdsinr ('idaysgraz',0,366,idaysgraz)
              call rdsinr ('idaysgrazpot',0,366,idaysgrazpot)
              
              count = max(ilvold,ilvoldpot)
              call rdfdor ('sla',0.d0,1.0d3,sla,366,count)
              call rdfdor ('slapot',0.d0,1.0d3,slapot,366,count)
              call rdfdor ('lvage',0.d0,1.0d3,lvage,366,count)
              call rdfdor ('lvagepot',0.d0,1.0d3,lvagepot,366,count)
              call rdfdor ('lv',0.d0,1.0d5,lv,366,count)
              call rdfdor ('lvpot',0.d0,1.0d5,lvpot,366,count)
              
              if (swanthesis .eq. 0) then
                flanthesis = .false.
              else
                flanthesis = .true.
              endif
              if (swgrazing .eq. 0) then
                flgrazing = .false.
              else
                flgrazing = .true.
              endif
              if (swgrazingpot .eq. 0) then
                flgrazingpot = .false.
              else
                flgrazingpot = .true.
              endif
        
            else

              daycrop = 0
              flCropEmergence = .false.

            endif
          else

            daycrop = 0
            flCropEmergence = .false.

          endif
        endif

!       close inifile
        close(ini)

      endif

      return
      end

! ----------------------------------------------------------------------
      subroutine rddre (wls1,wlp1)

! ----------------------------------------------------------------------
!     UpDate             : 20080109
!     Date               : 20010605                   
!     Purpose            : reading multilevel drainage characteristics 
!                          and specification of the surface water system
!                          for a period up to one year;  
!
! --- 1 Reading extended drainage input from .dra-file
! --- 2 Initializations
!
! --- Initializations:
! -1- wlp1: water level in primary system   (SWSRF = 3)
! -2- wls1: water level in secondary system (SWSRF = 2 or 3, SWSEC = 1 or 2) 
! -3- HBWEIR(IMPER) in case of table discharge relation (SWQHR = 2)
! -4- NUMADJ: number of target level adjustments
! -5- sttab(22,2) table with storage as a function of water level
! ---    sttab(i,1) contains levels:
! ---    i=1: +100 cm; i=2: 0 cm; i=22: bottom of deepest dr. medium
! ---    sttab(i,2) contains storage expressed as surface layer 
! ----------------------------------------------------------------------
      use variables, only: drfil,pathdrain,nrlevs,nrpri,nrsec,l,zbotdr,widthr,taludr,rdrain,rinfi,rentry,rexit,gwlinf,swdtyp,    &
                           wlptab,swsec,osswlm,nmper,wlstar,impend,swman,wlsman,wscap,swqhr,hbweir,alphaw,betaw,nqh,hqhtab,qqhtab,dropr,         &
                           gwlcrit,wlstab,sttab,swstini,swst,wlsbak,swsrf,nphase,hcrit,hdepth,vcrit,nodhd,numnod,dz,wldip,numadj,logf,intwl,       &
                           swnrsrf,rsurfdeep,rsurfshallow,t1900,cofintfl,expintfl,swscre,swdivd,cofani,numlay,tstart,tend,swdislay,swtopdislay,ztopdislay,&
                           ftopdislay,SwTopnrsrf,swdivdinf,FacDpthInf
      implicit none
      include 'arrays.fi'

! --- global
      real(8)   wls1,wlp1
      
! --- local
      integer   dra,level(Madr),itab,i,ilev,getun2
      integer   iph,nrman1,nrman2,node,imper,imperb,imperi,ifnd
      integer   imper_4b(mamp),imper_4c(mamp),imper_4d(mamp)
      integer   imper_4e1(mamp*mamte),imper_4e2(mamp*mamte)
      integer   imptab(mamte),impphase(mamp*mamte), nmper2
      real(8)   wlp(mawlp),wls(mawlp),zb, hhtab(mamp),qhtab(mamp)
      real(8)   wlsman_t(mamp*mamte),gwlcrit_t(mamp*mamte)
      real(8)   hcrit_t(mamp*mamte),vcrit_t(mamp*mamte)
      real(8)   dates(mabbc)
      real(8)   wdepth,wbreadth,wvolum,dep,swstlev,sofcu,altcu,afgen
      logical   flweir(mamp),exists(mamp),flzero(mamp),rdinqr,rdinar
      character(len=80)  filnam
      character(len=200) messag
      real(8)   small
      data      small     /0.0001d0/

! ----------------------------------------------------------------------

! --- open file with extended drainage data
      filnam = trim(pathdrain)//trim(drfil)//'.dra'
      dra = getun2 (10,90,2)
      call rdinit(dra,logf,filnam)

! --- division of drainage fluxes
      call rdsinr ('swdivd',0,1,swdivd)
      if (swdivd.eq.0) then
          write(messag,'(3a)')                                          &
     &    ' Variabel SWDIVD=0 in input file : ',trim(filnam),           &
     &    ' this is not recommended and may cause numerical instability'
          call warn ('rddre',messag,logf,swscre)
      endif
      if (swdivd .eq. 1) then
        if(rdinqr('swdivdinf')) then
          call rdsinr ('swdivdinf',0,1,swdivdinf)
          call rdsdor ('FacDpthInf',0.d0,1.d0,FacDpthInf)
        else
          swdivdinf = 0
        endif
        if (rdinar('cofani')) then
          call rdfdor ('cofani',1.d-4,1000.d0,cofani,maho,numlay)
        else
          call rdsdor ('cofani',1.d-4,1000.d0,cofani(1))
        endif
      endif

! --- altitude of control unit (relative to reference level)
      call rdsdor ('altcu',-3.0d5,3.0d5,altcu)

! --- part 1

! --  number of drainage levels
      call rdsinr ('nrsrf',1,Madr,nrlevs)
      if (nrlevs.gt.5) then
          write(messag,'(3a)')                                          &
     &    ' Number of Drainage levels >5 in inputfile : ',trim(filnam), &
     &    '   part of the output is limited to 5 levels'
          call warn ('rddre',messag,logf,swscre)
      endif

! --  top of model dicharge layer, determined by factor or direct input
      if (swdivd .eq. 1) then
        call rdsinr ('swdislay',0,2,swdislay)
        if (swdislay .eq. 1) then
          call rdfinr ('swtopdislay',0,1,swtopdislay,madr,nrlevs)
          call rdfdor('ztopdislay',-1.0d4,0.0d0,ztopdislay,madr,nrlevs)
        elseif (swdislay .eq. 2) then
          call rdfinr ('swtopdislay',0,1,swtopdislay,madr,nrlevs)
          call rdfdor('ftopdislay',0.0d0,1.0d0,ftopdislay,madr,nrlevs)
        endif
      endif


! --- characteristics of each drainage level 
      call rdfinr ('lev',1,Madr,level,Madr,nrlevs)
      call rdfinr ('swdtyp',0,1,swdtyp,Madr,nrlevs)
      call rdfdor ('l',1.0d0,100000.0d0,l,Madr,nrlevs)
      call rdfdor ('zbotdre',(altcu-1.0d3),(altcu-1.0d-2),              &
     &              zbotdr,Madr,nrlevs)
      call rdfdor ('gwlinf',-10000.0d0,0.0d0,gwlinf,Madr,nrlevs)
      call rdfdor ('rdrain',1.0d0,1.0d5,rdrain,Madr,nrlevs)
      call rdfdor ('rinfi',1.0d0,1.0d5,rinfi,Madr,nrlevs)
      call rdfdor ('rentry',0.0d0,10.0d0,rentry,Madr,nrlevs)
      call rdfdor ('rexit',0.0d0,10.0d0,rexit,Madr,nrlevs)
      call rdfdor ('widthr',0.0d0,10000.0d0,widthr,Madr,nrlevs)
      call rdfdor ('taludr',1.0d-2,5.0d0,taludr,Madr,nrlevs)

! -   conversions and security checks....
      do i = 1, nrlevs
!       conversions
        l(i) = l(i)*100.0d0
        zbotdr(i) = zbotdr(i) -altcu
        if (swdivd.eq.1 .and. swdislay.eq.1) then
          ztopdislay(i) = ztopdislay(i) - altcu
          if(ztopdislay(i).lt.zbotdr(i)) then
            messag = 'Ztopdislay cannot be below Zbotdre, verify input!' 
            call fatalerr ('Rddra',messag)
          endif
        endif
! -     security checks....
! -     0 - verify levels-index
! -     1 - levels must be ordered, starting with the deepest
! -     2 - zbotdr always below surface level
! -     3 - gwlinf must be below bottom of deepest drainage medium
! -     4 - widthr must be greater than zero for non-tube drain systems
        if (level(i).ne.i) then
          messag = 'Drainage level index is not consistent' 
          call fatalerr ('Rddra',messag)
        endif
        if (i .gt. 1) then
          if (zbotdr(i) .lt. zbotdr(i-1)) then
            messag = 'Levels must be ordered, starting with deepest' 
            call fatalerr ('Rddra',messag)
          endif
        endif
        if (gwlinf(i).gt.zbotdr(i)) then
            messag = 'Gwlinf should be lower than Zbotdr !' 
            call fatalerr ('Rddra',messag)
        endif
        if (swdtyp(i).eq.0 .and.widthr(i).lt.small) then
            write(messag,'(a,f10.5)')                                   &
     &              'widthr must be greater than',small
            call fatalerr ('Rddra',messag)
        endif

      enddo

! --- type of highest drainage level
      call rdsinr ('swnrsrf',0,2,swnrsrf)
      if (swnrsrf.eq.1) then
        call rdsdor ('rsurfdeep',0.001d0,1000.0d0,rsurfdeep)
        call rdsdor ('rsurfshallow',0.001d0,1000.0d0,rsurfshallow)
      else if (swnrsrf.eq.2) then
        call rdsdor ('cofintfl',0.01d0,10.0d0,cofintfl)     
        call rdsdor ('expintfl',0.1d0,1.0d0,expintfl) 
      endif


! kroes 20080707 : allow disabling of option
      if (swdivd .eq. 1) then
        if (swnrsrf.gt.0)  call rdsinr ('SwTopnrsrf',0,1,SwTopnrsrf)
      endif


! --- part 2a
      call rdsinr ('swsrf',1,3,swsrf)
      if (swsrf.eq.3) then
        nrpri = 1
        nrsec = nrlevs-nrpri
      elseif (swsrf.eq.2) then
        nrpri = 0
        nrsec = nrlevs-nrpri
      elseif (swsrf.eq.1) then
! ---   no surface water system...ready
        close (dra)
        return
      endif

      if (swdtyp(1+nrpri).ne.0) then
        messag = 'Deepest sec. level must be open.' 
        call fatalerr ('Rddra',messag)
      endif

! --- part 2b

! --- read table with water levels in the primary system (SWSRF=3)
      if (swsrf.eq.3) then
! ---   init table
        do 8 i = 1,2*mawlp
          wlptab(i) = 0.0d0
    8   continue
! -     read and store
        call rdatim ('date1',dates,mawlp,ifnd)
! -     at least one date must be within simulation period
        call checkdate(ifnd,dates,tend,tstart,'date1 ',                 &
     &                 'readswap//swsrf=3        ')
        call rdfdor                                                     &
     &        ('wlp',(altcu-1000.0d0),(altcu+200.0d0),wlp,mawlp,ifnd)
        do i = 1, ifnd
          wlptab(i*2) = wlp(i) - altcu
          wlptab(i*2-1) = dates(i)
        enddo
! ---   set initial value
        wlp1 = afgen (wlptab,2*mawlp,t1900-1.d0)

! ---   Ready in case of only a primary system ??
!       if (NRSEC.LT.1) then
!          CLOSE(DRE)
!          Return
!       endif 

      endif

! --- part 2c

      call rdsinr ('swsec',1,2,swsec)

! --- part 3

      if (swsec.eq.1) then
! ---   surface water level of secondary system is input
! ---   position file pointer

! ---   init table
        do 60 i = 1,2*mawls
          wlstab(i) = 0.0d0
   60   continue
! -     read and store
        call rdatim ('date2',dates,mawlp,ifnd)
! -     at least one date must be within simulation period
        call checkdate(ifnd,dates,tend,tstart,'date2 ',                 &
     &                 'readswap//swsec=1        ')

        call rdfdor                                                     &
     &         ('wls',(altcu-1000.0d0),(altcu+200.0d0),wls,mawlp,ifnd)
        do i = 1, ifnd
          wlstab(i*2) = wls(i) -altcu
          wlstab(i*2-1) = dates(i)
        enddo
! ---   set initial value
        wls1 = afgen (wlstab,2*mawls,t1900-1.0d0)

! ---   part 4a

! --- surface water level of secondary system is simulated
      elseif (swsec.eq.2) then
        call rdsdor ('wlact',dble(zbotdr(1+nrpri)+altcu)                &
     &       ,dble(altcu),wls1)
        wls1 = wls1 - altcu
        call rdsdor ('osswlm',0.0d0,10.0d0,osswlm)

! ---   part 4b

        call rdsinr ('nmper',1,mamp,nmper)
        wlstar = wls1

        call rdfinr ('imper_4b',1,nmper,imper_4b,mamp,nmper)
        call rdftim ('impend',impend,mamp,nmper)
        call rdfinr ('swman',1,2,swman,mamp,nmper)
        call rdfdor ('wscap',0.0d0,10.0d0,wscap,mamp,nmper)
        call rdfdor ('wldip',0.0d0,100.0d0,wldip,mamp,nmper)
        call rdfinr ('intwl',1,31,intwl,mamp,nmper)

! -     for each type of management: count number of periods 
        nrman1 = 0
        nrman2 = 0
        do imper = 1, nmper
          if (swman(imper).eq.2 .and. intwl(imper).lt.1) then
            messag = 'intwl (management interval) must be > = 1d'
            call fatalerr ('Rddra',messag)
          endif
          wldip(imper) = abs(wldip(imper))
          if (swman(imper).eq.1) then 
            nrman1 = nrman1+1
          elseif (swman(imper).eq.2) then
            nrman2 = nrman2+1
          else
            messag = 'SWMAN is out of range.'
            call fatalerr ('Rddra',messag)
          endif
        enddo
        if ((nrman1+nrman2).ne.nmper) then
          messag = 'nrman1+nrman2 does not match nmper'
          call fatalerr ('Rddra',messag)
        endif

! ---   type of discharge relationship
        call rdsinr ('swqhr',1,2,swqhr)

        IF (SWQHR.EQ.1) THEN

! ---     part 4c

          call rdsdor ('sofcu',0.1d0 ,100000.0d0,sofcu)

          call rdfinr ('imper_4c',1,nmper,imper_4c,mamp,nmper)
          zb = min(zbotdr(1),zbotdr(2))
          call rdfdor                                                   &
     &           ('hbweir',(altcu+zb),(altcu+100.0d0),hbweir,mamp,nmper)
          call rdfdor ('alphaw',0.1d0,50.0d0,alphaw,mamp,nmper)
          call rdfdor ('betaw',0.5d0,3.0d0,betaw,mamp,nmper)

! ---     convert and check values
!         initialise
          imperb = 0
          imperi = 0
          do imper = 1, nmper
            hbweir(imper) = hbweir(imper)-ALTCU
! ---       correction for units
            alphaw(imper) = alphaw(imper)* (8.64d0* 100.0d0**           &
     &                      (1.0d0-betaw(imper))/SOFCU)     
            if (hbweir(imper).lt.zbotdr(1+NRPRI)) then
              messag = 'Weir crest level is below bottom of'            &
     &                 //' deepest channel of secondary system.'
              call fatalerr ('Rddra',messag)
            endif
! --- check for target level above channel bottom when supply is 
! --- attempted (system may never become dry in this case)
            if (swman(imper).eq.1 .and. wscap(imper).gt.1.d-7 .and.     &
     &         (hbweir(imper)-wldip(imper)).lt.                         &
     &         (zbotdr(1+NRPRI)+1.d-4)) then
              messag = 'HBWEIR/WLDIP !'                                 &
     &                 //' Supply not possible =< zbotdr !'
              call fatalerr ('Rddra',messag)
            endif

            if (imper.ne.imperb) then
              imperb = imper
              imperi = imperi + 1
            else
              messag = '4c imper not unique'
              call fatalerr ('Rddra',messag)
            endif
          enddo

! ---     check number of records.. 
          if (imperi.ne.nmper) then
              messag = 'part 4c - not enough records'
              call fatalerr ('Rddra',messag)
          endif

        elseif (swqhr.eq.2) then

! ---     part 4d

! --      initialise
          imperb = 0
          imperi = 0
          do 34 imper = 1,nmper
            nqh (imper) = 0
            flweir (imper) = .false.
            exists (imper) = .false.
            flzero (imper) = .false.
34        continue

! --      read and store
          call rdainr ('imper_4d',1,nmper,imper_4d,mamp,ifnd)
          call rdfinr ('imptab',1,mamte,imptab,mamte,ifnd)
          call rdfdor                                                   &
     &         ('htab',(altcu-1000.0d0),(altcu+100.0d0),hhtab,mamp,ifnd)
          call rdfdor ('qtab',0.0d0,500.0d0,qhtab,mamp,ifnd)

! ---     convert
          do i = 1, ifnd
            hqhtab(imper_4d(i),imptab(i)) = hhtab(i) - Altcu
            qqhtab(imper_4d(i),imptab(i)) = qhtab(i)
          enddo

! ---     check values
          do i = 1, ifnd
            imper = imper_4d(i)
            itab = imptab(i)
            nqh(imper) = nqh(imper)+1
            exists(imper) = .true.
            if (qqhtab(imper,itab).lt.0.000001d0) flzero(imper)=.true.
            if (hqhtab(imper,itab).lt.zbotdr(1+nrpri)) then
              messag = 'Level in q-h table below bottom of'             &
     &                 //' deepest channel of secondary system.'
              call fatalerr ('Rddra',messag)
            endif
            if (imper.ne.imperb) then
              imperb = imper
              imperi = imperi + 1
            endif

! -         establish hbweir (at level where qtab = 0) 
            if (qqhtab(imper,itab).lt.0.0001d0.and..not.flweir(imper))  &
     &      then
              hbweir(imper) = hqhtab(imper,itab)
              flweir(imper) = .true.
            endif 

! -         consistency checks
            if (nqh(imper).ne.itab) then
              messag = 'qh-table / imper - itab mismatch'
              call fatalerr ('Rddra',messag)
            endif
            if (itab.eq.1) then
              if (abs(hqhtab(imper,itab)-100.0d0).gt.0.00001d0) then
              messag = 'First value in htab should be altcu+100.0'
              call fatalerr ('Rddra',messag)
              endif
            endif
            if (itab.gt.1) then
              if ((hqhtab(imper,itab).ge.hqhtab(imper,itab-1)).or.      &
     &           (qqhtab(imper,itab).gt.qqhtab(imper,itab-1))) then
                 messag = 'qh-table - no descending values'
                 call fatalerr ('Rddra',messag)
              endif
            endif

          enddo

! ---     check number of periods.. 
          if (imperi.ne.nmper) then
            messag = 'part 4d - number of periods incorrect.'
            call fatalerr ('Rddra',messag)
          endif

! ---     check that QQHTAB goes down to zero
          do imper = 1,nmper
            if (exists(imper).and..not.flzero(imper)) then
              messag = 'qqhtab not going down to zero.'
              call fatalerr ('Rddra',messag)
            endif
          enddo

        endif

! ---   part 4e1

        if (nrman2.gt.0) then

          imperb = 0
          imperi = 0

! ---     read table with drop rates (length of table must be equal 
!            to the nr of management periods referred to in tabel 4e2
          call rdainr('imper_4e1',1,nmper,imper_4e1,(mamp*mamte),nmper2)
          call rdfdor('dropr',0.0d0,100.0d0,dropr,(mamp*mamte),nmper2)
          call rdfdor('hdepth',-100d0,0.0d0,hdepth,(mamp*mamte),nmper2)

          do i = 1, nmper2
            imper = imper_4e1(i)
            hdepth(i) = -abs(hdepth(i))

! ---       determine compartment number related to hdepth
            dep = 0.0d0
            node = 1
            do while (hdepth(i).lt.(dep-1.0d-6).and.node.le.numnod) 
              nodhd(imper) = node
              dep = dep-dz(node)
              node = node + 1
            enddo

            if (swman(imper) .ne. 2) then
              messag = '#4e1 swman - imper mismatch.'
              call fatalerr ('Rddra',messag)
            endif
            if (imper.ne.imperb) then
              imperb = imper
              imperi = imperi + 1
            else
              messag = 'Two drop rates at the same period.'
              call fatalerr ('Rddra',messag)
            endif
          enddo
          if (imperi .ne. nrman2) then
            messag = '#4e1 number of periods for drop rate incorrect.'
            call fatalerr ('Rddra',messag)
          endif

! ---   part 4e2

          imperb = 0
          imperi = 0

          do 54 imper = 1,nmper
            nphase(imper) = 0
   54     continue

! --      read and store
          call rdainr('imper_4e2',1,nmper,imper_4e2,(mamte*mamp),ifnd)
          call rdfinr('impphase',1,mamte,impphase,(mamte*mamp),ifnd)
          call rdfdor                                                   &
     &       ('wlsman',(altcu-500.0d0),altcu,wlsman_t,(mamp*mamte),ifnd)
          call rdfdor                                                   &
     &       ('gwlcrit',-500.0d0,0.0d0,gwlcrit_t,(mamp*mamte),ifnd)
          call rdfdor('hcrit',-1000.0d0,0.0d0,hcrit_t,(mamp*mamte),ifnd)
          call rdfdor('vcrit',0.0d0,20.0d0,vcrit_t,(mamp*mamte),ifnd)

! ---     convert
          do i = 1, ifnd
            wlsman(imper_4e2(i),impphase(i)) = wlsman_t(i) - Altcu
            gwlcrit(imper_4e2(i),impphase(i)) = gwlcrit_t(i)
            hcrit(imper_4e2(i),impphase(i)) = hcrit_t(i)
            vcrit(imper_4e2(i),impphase(i)) = vcrit_t(i)
          enddo

! ---     check values
          do i = 1, ifnd
            imper = imper_4e2(i)
            iph = impphase(i)

            nphase(imper) = nphase(imper) + 1
            if (wlsman(imper,iph).lt.zbotdr(1+nrpri)) then
              messag = 'Level of automatic weir below bottom'           &
     &                 //' deepest channel of secondary system.'
              call fatalerr ('Rddra',messag)
            endif
            if (swman(imper) .ne. 2) then
              messag = '#4e2 swman - imper mismatch'
              call fatalerr ('Rddra',messag)
            endif
            if (imper.ne.imperb) then
              imperb = imper
              imperi = imperi + 1
            endif
          enddo

          if (imperi .ne. nrman2) then
            messag = '#4e2 inconsistency between the nr of periods'     &
     &            //'with automatic weir (tabel 4e (IMPER4e2)) and'     &
     &            //' swman in tabel4b (IMPER_4b)'
            call fatalerr ('Rddra',messag)
          endif

! ---     consistency checks WLSMAN, GWLCRIT, HCRIT and VCRIT
          do 446 imper = 1,nmper
            if (swman(imper).eq.2) then
              if (abs(gwlcrit(imper,1)) .gt. 0.01d0) then
                messag = '#4e2 - gwlcrit(1) must be 0.'
                call fatalerr ('rddra',messag)
              endif
              if (abs(vcrit(imper,1)) .gt. 0.01d0) then
                messag = '#4e2 - vcrit(1) must be 0.'
                call fatalerr ('rddra',messag)
              endif
              if (abs(hcrit(imper,1)) .gt. 0.01d0) then
                messag = '#4e2 - hcrit(1) must be 0.'
                call fatalerr ('rddra',messag)
              endif
              if (hbweir(imper).gt. (wlsman(imper,1)-0.99999d0)) then
                messag = '#4e2 - HBWEIR within 1 cm of wlsman(1)'
                call fatalerr ('Rddra',messag)
              endif

              do 448 iph = 2,nphase(imper)
                if (wlsman(imper,iph).lt.wlsman(imper,iph-1)) then
                  messag = '#4e2 - wlsman inconsistent'
                  call fatalerr ('rddra',messag)
                endif
                if (gwlcrit(imper,iph).gt.gwlcrit(imper,iph-1)) then
                  messag = '#4e2 - gwlcrit inconsistent'
                  call fatalerr ('rddra',messag)
                endif
                if (hcrit(imper,iph).gt.hcrit(imper,iph-1)) then
                  messag = '#4e2 - hcrit inconsistent'
                  call fatalerr ('rddra',messag)
                endif
                if (vcrit(imper,iph).lt.vcrit(imper,iph-1)) then
                  messag = '#4e2 - vcrit inconsistent'
                  call fatalerr ('rddra',messag)
                endif
448           continue
            endif 
446       continue
        endif
      endif
! ----------------------------------------------------------------------
! --- initialize counter for number of target level adjustments
      numadj = 0

! --- sttab(i,1) contains depths

! --- sw-levels, to be used in piece-wise linear functions (first level)
! --- is 100 cm above soil surface to allow for situations with ponding)
      sttab(1,1) = 100.0d0
      sttab(2,1) =   0.0d0
      do 100 i = 3,22
! ---   layer between surface level and deepest
! ---   drain/channel bottom is divided into 20 compartments
        sttab(i,1) = zbotdr(1+nrpri)*(i-2)/20.0d0
100   continue

! --- sttab(i,2) contains storage expressed in cm

! --- calculation of surface water storage (in cm, i.e. volume
! --- per unit area), as a function of sw-level, only for open channels:
      do 104 i = 1,22
        sttab(i,2) = 0.0d0
        do 108 ilev = 1+nrpri,nrlevs
          if (swdtyp(ilev).eq.0.and.sttab(i,1).gt.zbotdr(ilev)) then
! ---       for levels above soil surface the volume-increment is
! ---       computed for a rectangle, and not a trapezium
            if (sttab(i,1) .le. 0.0d0) then
              wdepth = sttab(i,1)-zbotdr(ilev)
              wvolum = wdepth*(widthr(ilev)+wdepth/taludr(ilev))
            else
              wdepth = -zbotdr(ilev)
              wvolum = wdepth*(widthr(ilev)+wdepth/taludr(ilev))
              wbreadth = widthr(ilev) + 2*wdepth/taludr(ilev)
              wdepth = sttab(i,1)
              wvolum = wvolum + wbreadth*wdepth
            endif
            sttab(i,2) = sttab(i,2)+wvolum/l(ilev)
          endif
108     continue
104   continue

! --- initial storage swstini
      swstini = swstlev (wls1)
      swst = swstini

! --- initialize memorization of wls for most recent 4 timesteps 
      do i=1,4
        wlsbak(i) = 0.0d0
      enddo

! --- close input file with lateral boundary conditions
      CLOSE (DRA)         

      RETURN
      END

! ----------------------------------------------------------------------
      subroutine checkdate(ifnd,dates,tend,tstart,namedat,topic) 
! ----------------------------------------------------------------------
!     Date               : April 2006   
!     Purpose            : check range of input dates to range of simulation period
! ----------------------------------------------------------------------
      implicit none

! --- global
      integer   ifnd
      real(8)   dates(ifnd), tend, tstart
      character(len=200) messag
      character(len=5)   namedat
      character(len=*)   topic

! ----------------------------------------------------------------------
! --- local
      integer i
      logical fldaterr

! ----------------------------------------------------------------------

!   - at least one input date must be within simulation period or 
!     simulation period should be completely within range of input dates
      fldaterr = .true.
      i = 1
      do while (fldaterr .and. i.le.ifnd)
         if(dates(i).gt.tstart-1.d-6.and.dates(i).lt.tend+1.d-6)        &
     &      fldaterr=.false.
         i = i + 1
      enddo
      if(fldaterr) then
         if(dates(1).lt.tstart+1.d-6 .and. dates(ifnd).gt.tend-1.d-6)   &
     &      fldaterr=.false.
      endif
      if(fldaterr) then
        messag = 'Fatal '//namedat//                                    &
     &           ', no input date within simulation period'
        call fatalerr(topic,messag)
      endif

      return
   end

!
! Version taken from FUSSIM, and adapted for use in SWAP
! December 2017, Marius Heinen
!
!-----------------------------------------------------------------------*
! Subroutine GTSOIL                                                     *
!                                                                       *
! Author  : Marius Heinen, Kees Rappoldt AB-DLO, Haren                  *
! Date    : July 1995                                                   *
! Purpose : Get parameters of the Van Genuchten-Mualem functions for    *
!           hydraulic properties for the specified soil name from soil  *
!           database.                                                   *
!                                                                       *
! Oct. 1999: two additional variables for extended van Genuchten-Mualem *
!            functions, i.e. TKK and TWCK. see explanation in:          *
! Heinen M., 1999. Extension to the van Genuchten-Mualem description    *
!    of the hydraulic properties in FUSSIM2. AB, internal note, 7 p.    *
!                                                                       *
!                                                                       *
!  FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)      *
!  name   type   meaning                                  units  class  *
!  ----   ----   -------                                  -----  -----  *
!  IUNIT  I      unit number used (and IUNIT+1)             -      I    *
!  FILNAM CH*    Name of input file                         -      I    *
!  SNAME  CH*    name of soil for which the parameters                  *
!                have to be returned                        -      I    *
!  KS     R8     hydraulic conductivity at saturation     cm/d     O    *
!  N      R8     curve shape parameter                      -      O    *
!  L      R8     curve shape parameter                      -      O    *
!  WCS    R8     saturated volumetric water content         -      O    *
!  WCR    R8     residual volumetric water content          -      O    *
!  ALPHA  R8     curve shape parameter                    1/cm     O    *
!  ALPHAW R8     curve shape parameter of main wetting                  *
!                curve                                    1/cm     O    *
!  KSAT   R8     KsatExm                                  cm/d     O    *
!  HE     R8     H_enpr                                   cm       O    *
!                                                                       *
! Functions and subroutines used:                                       *
! from TTUTIL : RD* routines, UPPERC, IFINDC                            *
!                                                                       *
!-----------------------------------------------------------------------*
SUBROUTINE GTSOIL (iTask,IUNIT,FILNAM,SNAME,WCR,WCS,ND,ALPHA,ALPHAW,L,KS,KSAT,HE,flksatexm)

IMPLICIT NONE
! global parameters
INTEGER                                :: IUNIT,iTask
real(8)                                :: KS,ND,L,WCS,WCR,ALPHA,ALPHAW,KSAT,HE
CHARACTER(len=*)                       :: FILNAM,SNAME
logical                                :: flksatexm

! local parameters
INTEGER                                :: IL,I,IS
INTEGER,             PARAMETER         :: IDECL = 500
real(8),             dimension(IDECL)  :: TKS,TND,TL,TWCS,TKSAT,THE,TWCR,TALP,TALPW
CHARACTER(len=40),   dimension(IDECL)  :: TBNAME
CHARACTER(len=40)                      :: LNAME
logical                                :: RDinqr

! functions
INTEGER                                :: IFINDC
SAVE

select case (iTask)

case (1)
!  read all soil names and parameters input section ; analyse input file
!   CALL STRIP (FILNAM,' ')
   CALL RDINIT (IUNIT, 0, FILNAM)
!     get values from file
!     SoilName Ores Osat Npar Alfa AlfaW Lexp Ksat H_enpr KsatExm
      CALL RDACHA ('SoilName', TBNAME, IDECL, IL)
      CALL RDFDOR ('Ores'    ,   0.0d0  ,   1.0d0 , TWCR  , IDECL, IL)
      CALL RDFDOR ('Osat'    ,   0.0d0  ,   1.0d0 , TWCS  , IDECL, IL)
      CALL RDFDOR ('Npar'    ,   1.001d0,   9.0d0 , TND   , IDECL, IL)
      CALL RDFDOR ('Alfa'    ,   1.0d-4 , 100.0d0 , TALP  , IDECL, IL)
      CALL RDFDOR ('AlfaW'   ,   1.0d-4 , 100.0d0 , TALPW , IDECL, IL)
      CALL RDFDOR ('Lexp'    , -25.d0   ,  25.0d0 , TL    , IDECL, IL)
      if (RDinqr('Ksatfit')) then      ! to allow downward compatibility
         CALL RDFDOR ('Ksatfit' ,   1.0d-5 ,   1.0d5 , TKS   , IDECL, IL)
      else
         CALL RDFDOR ('Ksat'    ,   1.0d-5 ,   1.0d5 , TKS   , IDECL, IL)
      end if
      if (RDinqr('Ksatexm')) then
         flksatexm = .true.
         CALL RDFDOR ('KsatExm' ,  -1.0d5  ,   1.0d5 , TKSAT , IDECL, IL)
      end if
      CALL RDFDOR ('H_enpr'  , -40.0d0  ,   0.0d0 , THE   , IDECL, IL)
      DO I = 1, IL
         CALL UPPERC (TBNAME(I))
      END DO
   CLOSE (IUNIT)

case (2)
!  check name length
   I = len_trim (SNAME)
   IF (I.GT.LEN(LNAME)) CALL FATALERR ('gtsoil','Soil name too long')

!  convert SNAME to local and all names to uppercase
   LNAME = SNAME
   CALL UPPERC (LNAME)

!  find position of LNAME in TBNAME and store as IS; check
   IS = IFINDC (TBNAME, IDECL, 1, IL, LNAME)
   IF (IS.EQ.0) CALL FATALERR ('gtsoil','Unknown soil')

!  set output variables
   WCR    = TWCR (IS)
   WCS    = TWCS (IS)
   ALPHA  = TALP (IS)
   ALPHAW = TALPW (IS)
   ND     = TND (IS)
   L      = TL (IS)
   KS     = TKS (IS)       ! Ksatfit
   KSAT   = TKSAT (IS)     ! Ksatexm
   HE     = THE (IS)       ! air-entry value

case default
   call fatalerr ('gtsoil','Illegal value for iTask')
end select

! ready
RETURN
END SUBROUTINE GTSOIL
