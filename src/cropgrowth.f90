! File VersionID:
!   $Id: cropgrowth.f90 380 2018-05-28 14:01:08Z heine003 $
! ----------------------------------------------------------------------
      subroutine CropGrowth(task) 
! ----------------------------------------------------------------------
!     UpDate             : May 2014   
!     Date               : Aug 2004   
!     Purpose            : Call proper crop routines for initialization,
!                          calculation of rate/state variables and output
! ----------------------------------------------------------------------

      use variables
      implicit none

      integer task
      integer i, node
      real(8) sumtmin
      real(8) afgen
      
      ! assimilation
      real(8) dayl,cosld,sinld
      real(8) effc,amax
      real(8) dtgapot,dtga
      ! only for bulb crops (tulips etc..)
      real(8) respmo,decrmo,remo,factblb
      include 'params.fi'

      select case (task)

      case (1)

! === initialization ==========================================================

      ! find active crop
      icrop = 1
      flCropCalendar = .false.
      do while (.not. flCropCalendar) 
        
        if (cropstart(icrop) .lt. 1.d0) exit
        
        if (t1900 - cropstart(icrop) .gt. -tiny                         &
     &                 .and. t1900 - cropend(icrop) .lt. tiny) then
          flCropCalendar = .true.
        else
          icrop = icrop + 1
        endif
      enddo
      
! --- bare soil condition  ----------------------------------------------------
      if (.not. flCropEmergence .or. flCropHarvest) then          
        call nocrop ()
      endif

! --- check crop emergence ----------------------------------------------------

      ! reset if new crop
      if (flCropCalendar) then
        if (dabs(t1900 - cropstart(icrop)) .lt. tiny) then
          call InitializeCrop
          flCropReadFile  = .true.
          flCropEmergence = .true.
          if (croptype(icrop) .le. 2) then
            flCropEmergence = .false.
          endif
        endif
      endif

! --- Preparation, Sowing and Germination of arable crop growth ---------------
      if (flCropCalendar .and. .not. flCropHarvest .and. croptype(icrop) .le. 2) then
        
        ! check crop preparation, sowing and germination (of previous day)
        if (.not. flCropEmergence) then
          if (flCropPrep .and. flCropSow .and. flCropGerm) then
            swinco          = -99
            flCropReadFile  = .true.
            flCropEmergence = .true.
          endif
        endif
        
        ! Initialize preparation, sowing and germination
        if (.not. flCropEmergence) then
          if (flCropReadFile) then
            call ArableLandGerm(1)
          endif
        endif

        if (.not. flCropEmergence .and. .not. flCropHarvest) then
          
          ! Preparation before crop growth
          if (.not. flCropPrep) then
            call ArableLandGerm(2)
          endif
          
          ! Sowing before crop growth
          if (flCropPrep .and. .not. flCropSow) then
            call ArableLandGerm(3)
          endif

          ! Germination of arable crop growth
          if (flCropPrep .and. flCropSow) then
            call ArableLandGerm(4)
          endif

        endif
      
      endif

! --- Initialization crop conditions ------------------------------------------
      
      if (flCropCalendar .and. .not. flCropHarvest) then            
        
        if (flCropReadFile) then
          
          ! fixed crop development
          if (croptype(icrop) .eq. 1 .and. flCropEmergence) call CropFixed(1)
          
          ! detailed crop growth
          if (croptype(icrop) .eq. 2 .and. flCropEmergence) call Wofost(1)
          
          ! detailed grass growth
          if (croptype(icrop) .eq. 3) call Grass(1)
          
          flCropReadFile = .false.

        endif

        ! update crop daynumber
        daycrop = daycrop + 1
        
        ! open crp-file
        if (swcrp.eq.1) call CropOutput(1)
        
        ! set correction of CO2 impact
        call FacCO2()
        
      endif

      ! set running average of minimum temperature (only for detailed crop growth)
      if (flCropEmergence .and. croptype(icrop).ge.2) then
        nofd = min(nofd+1, 7)
        sumtmin = 0.0d0
        do i = nofd,2,-1
          atmin7(i) = atmin7(i-1)
          sumtmin = sumtmin + atmin7(i)
        end do
        i = 1
        atmin7(i) = tmn
        sumtmin = sumtmin + atmin7(i)
        tmnr = sumtmin / nofd
      else
        nofd = 0
      endif

      ! determine lowest compartment containing roots
      node = 1
      do while (zbotcp(node) .gt. (-rd + 1.d-8))
        node = node + 1
      end do
      noddrz = node        
      
      ! calculate potential and actual assimilation
      if (flCropEmergence .and. croptype(icrop).ge.2) then
          
! check DAYNR during the day!!!!!!          
          
        ! phenological development rate 
        call astro (daynr+1,lat,rad,dayl,daylp,sinld,cosld,difpp,atmtr,dsinbe)
      
        ! only for bulb crops (tulips etc..)
        if(swbulb.eq.1) then
          ! remobilisation of carbohydrates from planted material
          if (plwt.le.(0.0002d0*pld)) then
            ! no remobilisation at minimum weight motherbulb
            respmo = 0.0d0
            remo = 0.0d0
          else
            ! decrease weight mother organ starts at emergence.
            ! decrease consists of respiration and remobilisation
            decrmo = plwt-(plwt*(2.71828d0**remoc))
            respmo = 0.025d0*(q10**((tavd-25.0d0)/10.0d0))*plwt
            if(respmo.lt.decrmo) then
              remo = decrmo - respmo
            else
              remo = 0.0d0
              respmo = decrmo
            end if
            ! weight motherbulb decreases by remobilisation and respiration
            plwt = plwt - remo - respmo
          endif
        endif

        ! daily gross assimilation
        effc = fco2eff * eff
        if (croptype(icrop) .eq. 2) amax = fco2amax * afgen (amaxtb,30,dvs) * afgen (tmpftb,30,tavd)
        if (croptype(icrop) .eq. 3) amax = fco2amax * afgen (amaxtb,30,dble(daycrop)) * afgen (tmpftb,30,tavd)        

        ! potential assimilation
        call totass (dayl,amax,effc,laipot,kdif,rad,difpp,dsinbe,sinld,cosld,dtgapot)

        ! correction for low minimum temperature
        dtgapot = dtgapot * afgen (tmnftb,30,tmnr)

        ! potential assimilation in kg ch2o per ha
        pgasspot = dtgapot * 30.0d0/44.0d0
        
        ! only for bulb crops (tulips etc..)
        if(swbulb.eq.1) then
          ! assimilation is raised with remobilisation from motherbulb
          ! using a factor of 1.11 given by De Ruijter et al.(1993)
          factblb  = 1.11d0
          pgasspot = pgasspot + remo*factblb
        endif
        
        ! reduction due to limited attainable maximum yield
        if (swpotrelmf.eq.2) pgasspot = pgasspot * relmf
        
        
        ! actual assimilation
        call totass (dayl,amax,effc,lai,kdif,rad,difpp,dsinbe,sinld,cosld,dtga)

        ! correction for low minimum temperature
        dtga = dtga * afgen (tmnftb,30,tmnr)

        ! actual assimilation in kg ch2o per ha
        pgass = dtga * 30.0d0/44.0d0
        ! only for bulb crops (tulips etc..)
        if(swbulb.eq.1) then
          ! assimilation is raised with remobilisation from motherbulb
          ! using a factor of 1.11 given by De Ruijter et al.(1993)
          factblb = 1.11d0
          pgass   = pgass + remo*factblb
        endif
      
        ! reduction due to limited attainable maximum yield
        pgass = pgass * relmf
        
        ! nitrogen stress reduction of pgass
        if (flCropNut) then
          call NUTRIE (NLUE,WLV,WST,DVS,ANLV,ANST,NMXLV,NMAXLV,NMAXST,  &
     &      NMAXRT,LRNR,LSNR,NNI,RNFLV,RNFST,FRNX,FSTR)
          pgass = pgass * FSTR
        endif

      endif  

      return

      case (2)

! === calculation of potential crop rate and state variables =================

      if (flCropHarvest) return

! --- detailed crop growth -------------------------------------------------
      if (croptype(icrop) .eq. 2) then
        if (flCropEmergence) then
          call Wofost(2)
        endif
      endif  
! --- detailed grass growth  -----------------------------------------------      
      if (croptype(icrop) .eq. 3) then
        call Grass(2)
      endif

      return

      case (3)

! === calculation of actual crop rate and state variables ==================
      
      if (flCropHarvest) return
      
! --- fixed crop development -----------------------------------------------
      if (croptype(icrop).eq.1) then
        if(flCropEmergence) then
          call CropFixed(3)
        endif
      endif  
! --- detailed crop growth -------------------------------------------------
      if (croptype(icrop).eq.2) then
        if (flCropEmergence) then
          call Wofost(3)
         endif
      endif
! --- detailed grass growth  -----------------------------------------------
      if (croptype(icrop).eq.3) then
        call Grass(3)
      endif

      return

      case (4)

! === harvest of crop ======================================================
      
      if (flCropHarvest) return
     
      if (croptype(icrop).le.2 .and. flCropEmergence)then

        ! Check flHarvestDay
        if (swharv.eq.0) then
          if (dabs(t1900 - cropend(icrop) - 1.d0) .lt. 1.0d-3) then
            flHarvestDay = .true.
          endif
        else
          if (dvs.ge.dvsend .or. dabs(t1900 - cropend(icrop) - 1.d0) .lt. 1.0d-3) then
            flHarvestDay = .true.
          endif
        endif
        
        if (flCropEmergence .or. flHarvestDay) then
          if (croptype(icrop).eq.1) then
            call CropFixed(4) 
          endif  
          if (croptype(icrop).eq.2) then
            call Wofost(4) 
          endif  
        endif

      endif

! --- check timing of harvest
      
! --- fixed crop development -----------------------------------------------
      if (croptype(icrop).eq.1)then
        if (flHarvestDay) then
          flCropEmergence = .false.
          flCropHarvest   = .true.
        endif
      endif

! --- detailed crop growth -------------------------------------------------
      if (croptype(icrop).eq.2)then
        if (flHarvestDay) then
          flCropEmergence = .false.
          flCropHarvest   = .true.
        endif
      endif

! --- detailed grass growth ------------------------------------------------
      if (croptype(icrop).eq.3)then
        if (dabs(t1900 - cropend(icrop) - 1.d0) .lt. 1.0d-3) then
          flCropEmergence = .false.
          flCropHarvest   = .true.
        endif
      endif
      
      return      
      
      case default
         call fatalerr ('CropGrowth', 'Illegal value for TASK')
      end select

      return

      end
!
! ----------------------------------------------------------------------
      subroutine cropfixed (task)
! ----------------------------------------------------------------------
!     date               : august 2004                           
!     purpose            : simple crop growth routine for swap 
! ----------------------------------------------------------------------
      use variables
      implicit none

! --- local variables
      integer   i,task,lcc,swhydrlift
      real(8)   dummy,afgen,dtsum,dvr,watcon
      
! --- rooting
      real(8)   rrpot,rr
      
      include 'params.fi'
      save
! ----------------------------------------------------------------------

      select case (task)
      case (1)

! === initialization ===================================================
      
! --- read crop data
      call readcropfixed (icrop,cropfil(icrop),lcc,swhydrlift)
     
! --- maximum rooting depth
      if (swrd.eq.1) then
        rdm = rdmax
      else
        rdm = min(rdmax,rdc)
      endif

! --- skip next initialization if crop parameters are read from *.END file
      if (t1900 - tstart .gt. tiny .or. swinco .ne. 3 .or.              &
     &  dabs(t1900 - cropstart(icrop)) .lt. tiny) then
        
        dvs = 0.0d0

! --- actual rooting depth
        if (swrd.eq.1) then
          rd = afgen (rdtb,22,dvs)
          rd = min(rd,rdm)
        else
          rd = min(rdi,rdm)
        endif
        rdpot = rd

      endif

! --- initial lai or sc
      lai = afgen (gctb,(2*magrs),dvs)
      if (swgc.eq.2) then
        gc = lai
        lai = lai*3.0d0
      endif

! --- initial crop factor or crop height
      cf = afgen (cftb,(2*magrs),dvs)
      ch = afgen (chtb,(2*magrs),dvs)
      if (swcf.eq.3) then
        cfeic = afgen (cfeictb,(2*magrs),dvs)
      endif

! --- initial storage on canopy
      if (swinter.eq.3) then
        siccapact = siccaplai*lai
      endif

! --- initial dry weight of roots at soil surface; oxygen module
      W_root_ss = afgen (wrtb,(2*magrs),dvs)

! --- initial ratio root total respiration / maintenance respiration; oxygen module
      max_resp_factor = afgen (mrftb,(2*magrs),dvs)

! --- initialize matric flux potential and hleaf
      if (swdrought .eq. 2) then                                        
        call MatricFlux(1,h(1),1,dummy)
        if (swhydrlift .eq. 1) then
          flhydrlift = .true.
        else
          flhydrlift = .false.
        endif
        do i = 1,numnod
         twilt(i) = watcon(i,wiltpoint)
         hroot(i) = h(i)
        enddo
        hleaf = -2000.d0
      endif                                                             

      return          

      case (2)
      continue

! === calculate potential rate and state variables ======================
      case (3)

! === calculate actual rate and state variables ======================

! --- increase in temperature sum
      dtsum = max (0.0d0,tav-tbase)

! --- development rate
      if (idev.eq.1) then
        dvr = 2.0/lcc
      elseif (idev.eq.2) then
        if (dvs.lt.1.0d0) then
          dvr = dtsum/tsumea
        else
          dvr = dtsum/tsumam
        endif
      endif

! --- water stress
      if(dabs(ptra).lt.nihil) then
        reltr = 1.0d0
      else
        reltr = max(min(tra/ptra,1.0d0),0.0d0)
      endif

! ----integrals of the crop --------------------------------------------

! --- phenological development stage
      dvs = min(dvs+dvr,2.d0)
      tsum = tsum + dtsum

! --- leaf area index or soil cover fraction    
      lai = afgen (gctb,(2*magrs),dvs)
      if (swgc.eq.2) then
        gc = lai
        lai = lai*3.0d0
      endif

! --- crop factor or crop height
      cf        = afgen (cftb,(2*magrs),dvs)
      ch        = afgen (chtb,(2*magrs),dvs)
      if (swcf.eq.3) then
        cfeic     = afgen (cfeictb,(2*magrs),dvs)
      endif
      
! --- update canopy storage capacity
      if (swinter.eq.3) then
        siccapact = siccaplai*lai
      endif

! --- dry weight of roots at soil surface; oxygen module
      W_root_ss = afgen (wrtb,(2*magrs),dvs)

! --- ratio root total respiration / maintenance respiration; oxygen module
      max_resp_factor = afgen (mrftb,(2*magrs),dvs)

      case (4)
          
! --- root extension
      if (swrd.eq.1) then
        rdpot = afgen (rdtb,22,dvs)
        rdpot = min(rdpot,rdm)
        rd    = rdpot
      else
        rrpot = min (rdm-rdpot,rri)
        if (ptra.lt.nihil) rrpot = 0.0d0
        rdpot = rdpot + rrpot

        rr = min (rdm-rd,rri)
        if (ptra.lt.nihil .or. flWrtNonox)     rr = 0.0d0
        if (swdmi2rd.eq.1 .and. ptra.ge.nihil) rr = rr * tra/ptra
        rd = rd + rr
      endif

      return
      
      case default
         call fatalerr ('CropFixed', 'Illegal value for TASK')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine cropoutput(task) 
! ----------------------------------------------------------------------
!     Date               : Aug 2004   
!     Purpose            : open and write crop output files 
! ----------------------------------------------------------------------

      use variables
      implicit none

! --- local variables ------------------
      integer task,getun   !,numcrop
      character(len=200) messag
      character(len=160) filnam,filtext

      select case (task)
      case (1)

! === open crop output file and write headers =====================

! --- open crop output file
      if (flCropOpenFile) then

! ---   open crop output file and write general header (*.crp)
        if (trim(outfil).eq.trim(cropfil(1))) then
          Messag = 'The name of the input crop-file (''//trim(cropfil'//&
     &   '(icrop))//'') cannot be equal to the name of'                 &
     &   //'the output crop-file '//trim(outfil)//' Adjust a filename !'
          call fatalerr ('crops',messag)
        endif
        filnam = trim(pathwork)//trim(outfil)//'.crp'
        crp = getun (20,90)
        call fopens(crp,filnam,'new','del')
        filtext = 'output data of simple or detailed crop growth model'
        call writehead (crp,1,filnam,filtext,project)

! ---   write header fixed crop growth
        if (croptype(icrop) .eq. 1) call OutCropFixed(1)

! ---   write header detailed crop growth
        if (croptype(icrop) .eq. 2) call OutWofost(1)

! ---   write header detailed grass growth
        if (croptype(icrop) .eq. 3) call OutGrass(1)

        flCropOpenFile = .false.

      else

! ---   header for second and subsequent crops

! ---   write header fixed crop growth
        if (croptype(icrop).eq.1 .and. swheader.eq.1) call OutCropFixed(1)

! ---   write header detailed crop growth 
        if (croptype(icrop).eq.2 .and. swheader.eq.1) call OutWofost(1)

! ---   write header detailed grass growth
        if (croptype(icrop).eq.3 .and. swheader.eq.1) call OutGrass(1)

      endif

      return

      case (2)

! --- write actual data ----------------------------------------------------

! --- fixed crop file
      if (croptype(icrop) .eq. 1) call OutCropFixed(2)

! --- detailed crop growth 
      if (croptype(icrop) .eq. 2) call OutWofost(2)
        
! --- detailed grass growth
      if (croptype(icrop) .eq. 3) call OutGrass(2)
      
      return

      case (3)
! --- close crop output file ------------------------------------------------

      close (crp)

      case default
         call fatalerr ('CropOutput', 'Illegal value for TASK')
      end select

      return
      end 

! ----------------------------------------------------------------------
      subroutine nocrop ()
! ----------------------------------------------------------------------

      use variables, only: rd,rdpot,lai,laipot,cf,ch,albedo,rsc,tsum,dvs,               &
                           cwdmpot,cwdm,wsopot,wso,wlvpot,wlv,wstpot,wst,wrtpot,wrt
      implicit none

      rd      = 0.0d0
      rdpot   = 0.0d0
      lai     = 0.0d0
      laipot  = 0.d0
      cf      = 0.0d0
      ch      = 0.0d0
      albedo  = 0.23d0
      rsc     = 70.d0
      tsum    = 0.d0
      dvs     = 0.d0
      cwdmpot = 0.d0
      cwdm    = 0.d0
      wsopot  = 0.d0
      wso     = 0.d0
      wlvpot  = 0.d0
      wlv     = 0.d0
      wstpot  = 0.d0
      wst     = 0.d0
      wrtpot  = 0.d0
      wrt     = 0.d0
      
      return
      end

! ----------------------------------------------------------------------
      subroutine ArableLandGerm(task)
! ----------------------------------------------------------------------
!     update             : December 2017
!     date               : December 2017
!     purpose            : Crop growth
! ----------------------------------------------------------------------
      use variables
      
      implicit none
 
      include 'params.fi'
      integer  task,node
      real(8)  drz1,hrz1,pFz1
      real(8)  tsumemesub      
      
      select case (task)
          
! === Initialization preparation, sowing and germination ======================
      
      case (1)
        
        call readarablelandgerm(icrop,cropfil(icrop))
      
        if (.not. flCropEmergence) flCropReadFile = .false.
    
        return
      
! === Preparation before crop growth ==========================================
      
      case (2)  
      
        node   = 1
        dhPrep = h(node) - hPrep
        drz1   = -1.d0 * zPrep - dz(node)
        do while (drz1 .gt. 0.d0)
          node   = node + 1
          dhPrep = max(dhPrep,h(node) - hPrep)
          drz1   = drz1 - dz(node)
        enddo
        
        flCropPrep = .true.
        if (dhPrep .gt. 0.d0) then
          if (PrepDelay .lt. MaxPrepDelay) then
            dvs        = -0.3d0
            flCropPrep = .false.
            PrepDelay  = PrepDelay + 1
          endif
        endif
        
        SowDelay = PrepDelay
        
        return        

! === Sowing before crop growth ==========================================
      
      case (3)  
      
        node   = 1
        dhSow  = h(node) - hSow
        drz1   = -1.d0 * zSow - dz(node)
        do while (drz1 .gt. 0.d0)
          node   = node + 1
          dhSow  = max(dhSow,h(node) - hSow)
          drz1   = drz1 - dz(node)
        enddo
          
        node     = 1
        drz1     = -1.d0 * zTempSow - dz(node)
        do while (drz1 .gt. 0.d0)
          node   = node + 1
          drz1 = drz1 - dz(node)
        enddo
        
        dtempSow = min(tsoil(node) - TempSow,0.d0)
        
        flCropSow = .true.
        if (dtempSow .lt. 0.d0 .or. dhSow.gt.0.d0) then
          if (SowDelay .lt. MaxSowDelay) then
            dvs       = -0.2d0
            flCropSow = .false.
            SowDelay  = SowDelay + 1
          endif
        endif
        
        return        
        
! === Simulate germination ==============================================
      
      case (4)

        ! Optimal situation in case germination only depends on temperature (swgerm = 1)
        if (agerm .lt. 0.d0) then
          
          tsumemesub = tsumemeopt
            
        ! Germination depends on temperature and hydrological conditions (swgerm = 2)    
        else
          
          ! ---   calculate average pressure head of rootzone ---
          if (dabs(zgerm-0.d0) .lt. small) then
            hrz1 = h(1)
          else
            node = 0
            hrz1 = 0.0d0
            drz1 = zgerm * (-1.d0)
            do while (drz1 .gt. 0.d0)
              node = node + 1
              if (drz1 - dz(node) .ge. 0.d0) then
                hrz1 = hrz1+h(node)*dz(node)/(zgerm*(-1.d0))
              else
                hrz1 = hrz1+h(node)*drz1/(zgerm*(-1.d0))
              endif
              drz1 = drz1 - dz(node)
            enddo
          endif
        
          ! --- simulate germination time ---
          pFz1 = DLOG10(MAX(1.0d0,-hrz1))
          if (hrz1 .lt. hdrygerm) then
            ! Dry situation
            tsumemesub = agerm * pFz1 - cgerm
          elseif (hrz1 .ge. hdrygerm .and. hrz1 .le. hwetgerm) then       
            ! Optimal situation
            tsumemesub = tsumemeopt
          else
            ! Wet situation
            tsumemesub = -agerm * pFz1 + bgerm
          endif

        endif
        
        ! Update of tsumgerm, for the time step of 1 day
        if (tav .gt. TBASEM)then
          if( tav .lt. TEFFMX) then
            if(tsumemesub.lt.0.1d0) then
              tsumgerm = tsumgerm + (tav-TBASEM)
            else
              tsumgerm = tsumgerm +(tsumemeopt/tsumemesub)*(tav-TBASEM)
            endif
          else
            if(tsumemesub.lt.0.1d0) then
              tsumgerm = tsumgerm + (TEFFMX-TBASEM)
            else
              tsumgerm = tsumgerm +(tsumemeopt/tsumemesub)*(TEFFMX-TBASEM)
            endif
          endif
        endif
            
        ! Delay growth until tsumgerm is reached
        flCropGerm = .true.
        if (tsumgerm .lt. tsumemeopt) then
          dvs = -0.1d0 * max(1.d0 - (tsumgerm / tsumemeopt), 0.d0)
          flCropGerm = .false.
        else
          dvs = 0.d0
        endif
        
        return
        
      case default
        call fatalerr ('ArableLandGerm', 'Illegal value for TASK')
      end select

      return
      end
    
! ----------------------------------------------------------------------
      subroutine FacCO2()
! ----------------------------------------------------------------------
!     update             : February 2018
!     date               : ?
!     purpose            : Assimilation correction for CO2 changes in
!                          atmosphere (Lintul4) added by Iwan Supit
! ----------------------------------------------------------------------
      use variables

      implicit none
 
      integer   ifindi,indexyr
      real(8)   CO2,afgen
      character(len=200) messag
      
      ! initialize CO2 impact
      fco2amax = 1.0d0 ! factor to correct AMAX for CO2
      fco2eff  = 1.0d0 ! factor to correct EFF for CO2
      fco2tra  = 1.0d0 ! factor to correct TRA for CO2

      ! correction of CO2 impact
      if(flco2) then
        indexyr = ifindi (CO2year, mayrs, 1, mayrs, iyear)
        if (indexyr.lt.1 .or. indexyr.gt.mayrs) then
          Messag ='Input if CO2year or CO2ppm inconsistent, correct'
          call fatalerr ('wofost',messag)
        endif
        CO2 = CO2ppm(indexyr)
        fco2amax = afgen(CO2AMAXTB,30,CO2)
        fco2eff = afgen(CO2EFFTB,30,CO2)    
        fco2tra = afgen(CO2TRATB,30,CO2)    
      endif

      return
      end
    
! ----------------------------------------------------------------------
      subroutine wofost(task)
! ----------------------------------------------------------------------
!     update             : march 2015
!     date               : october 2004
!     purpose            : detailed crop growth routine
! ----------------------------------------------------------------------
      use variables
      use wofost_soil_interface

      implicit none
 
      integer   i1,task,swhydrlift,i

      real(8)   afgen,asrc,ccheck,cvf
      real(8)   laicr,lasum,mres
      real(8)   dalv,delt,dmi
      real(8)   drrt,drst,dslv,dteff,dtsum,dvr
      real(8)   dvred,fl,fo,fr,fs,drlv
      real(8)   fysdel,gass,gla,grlv,grrt,grst,grso
      real(8)   gwso,gwst,rmres
      real(8)   slat,teff,twlv,twst
      real(8)   lasumpot,watcon
      real(8)   gasspot,rmrespot,mrespot,asrcpot,dmipot
      real(8)   admipot,grrtpot,drrtpot,gwrtpot,grlvpot
      real(8)   dslvpot,restpot,dalvpot,drlvpot,gwsopot
      real(8)   glasolpot,slatpot,glapot,grstpot,drstpot,gwstpot
      real(8)   dslvtpot,twlvpot,twstpot
      real(8)   dummy
      character(len=200) messag
      
! --- rooting
      real(8)   rrpot,rr
      
! --- only for bulb crops (tulips etc..)
      real(8)   admi,grbl,gwbl,grblpot,gwblpot
      real(8)   twbl,twblpot

! --- n-p-k use
      real(8) NLAI   !,LRNR, LSNR
      real(8) NMAXSO, NPART, NFIXF !,NLUE
      real(8) NSLA, RNFRT, TCNT !,RNFLV,RNFST
      real(8) DVSNLT, DVSNT, RDRNS, FNTRT !, FRNX
      integer getun2
      real(8) INSW
      !real(8) NMAXLV,NMAXST,NMAXRT,NNI,FSTR
      !real(8) NMXLV(30)
      integer ILNMXL
      real(8) Fstress
      integer nut

      real(8) NDEMTO
      real(8) ANRT,ANSO,ANLVI,ANSTI,ANRTI,ANSOI !,ANLV,ANST
      real(8) RNLDRT,RNLDST,RNLDLV,NDEML,NDEMS,NDEMR,NDEMSO
      real(8) NSUPSO,ATN,RNSO,NLIMIT,NUPTR
      real(8) NFIXTR,ATNLV,ATNST,ATNRT,RNTLV,RNTST,RNTRT
      real(8) RNULV,RNUST,RNURT,NUPTT,NFIXTT
      real(8) RNLV,RNST,RNRT
      real(8) NLOSSL,NLOSSR,NLOSSS  
      real(8) NBALAN

      real(8) NdemandBioFix
      real(8) ombalan,wlvt0,wstt0,wsot0,wrtt0,storagediff,drso
      real(8) FraHarLosOrm_lv,FraHarLosOrm_st,FraHarLosOrm_so
      real(8) FraDeceasedLvToSoil
      real(8) HarLosOrm_rt, HarLosOrm_lv, HarLosOrm_st, HarLosOrm_so
      real(8) HarLosNit_rt, HarLosNit_lv, HarLosNit_st, HarLosNit_so
      real(8) HarLosOrm_dwst, HarLosOrm_dwso, HarLosOrm_dwlv

      real(8) HarLosNit_dwst, HarLosNit_dwso, HarLosNit_dwlv
      real(8) idwlvCrop,idwlvSoil, NLOSSLDeceasedLvToSoil
      character(len=200) filnam

! --- only for soybean
      integer swsoybean
      logical flrfphotoveg
      real(8) mg,dvsi,dvrmax1,dvrmax2,tmaxdvr,tmindvr,toptdvr
      real(8) rfmgphotop, rfmgtemp
      logical flphenodayl
      real(8) popt
      real(8) pcrt

! --- only for bulb crops (tulips etc..)
      real(8) wstem

      parameter (delt=1.0d0)
      include 'params.fi'

! --- only for vernalisation
      logical flvernalised
      real(8) r,vern,vernfac,vernrate,interpol

      save
! ----------------------------------------------------------------------
      select case (task)
      
      case (1)

! === initialization ====================================================
      
! --- read general crop data
      call readwofost (icrop,cropfil(icrop),swhydrlift,swsoybean,mg,dvsi,dvrmax1,dvrmax2, &
                       flrfphotoveg,tmaxdvr,tmindvr,toptdvr,popt,pcrt,flphenodayl,FraDeceasedLvToSoil)

! --- if crop based on calendar is still active, but already harvested
      if (flCropHarvest) return

! --- n-p-k 
      if( flCropNut) then
!        initialise and start reading
         filnam = trim(pathcrop)//trim(cropfil(icrop))//'.crp'
         nut = getun2 (10,90,2)
         call rdinit(nut,logf,filnam)

         CALL rdsdou ('LRNR', LRNR)
         CALL rdsdou ('LSNR', LSNR)
         CALL rdsdou ('NLAI', NLAI)
         CALL rdsdou ('NLUE', NLUE)
         CALL rdsdou ('NMAXSO', NMAXSO)
         CALL rdsdou ('NPART', NPART)
         CALL rdsdou ('NFIXF', NFIXF)
         CALL rdsdou ('NSLA', NSLA)
         CALL rdsdou ('RNFLV', RNFLV)
         CALL rdsdou ('RNFRT', RNFRT)
         CALL rdsdou ('RNFST', RNFST)
         CALL rdsdou ('TCNT', TCNT) 
         CALL rdsdou ('DVSNLT', DVSNLT)
         CALL rdsdou ('DVSNT', DVSNT)
         CALL rdsdou ('RDRNS', RDRNS)
         CALL rdsdou ('FNTRT', FNTRT)
         CALL rdsdou ('FRNX', FRNX)
         CALL rdadou ('NMXLV', NMXLV, 30, ILNMXL)
!        Read harvest losses (fractions of leaves, stems, storage organs)
         call rdsdor ('FraHarLosOrm_lv',0.d0,1.0d0,FraHarLosOrm_lv) 
         call rdsdor ('FraHarLosOrm_st',0.d0,1.0d0,FraHarLosOrm_st) 
         call rdsdor ('FraHarLosOrm_so',0.d0,1.0d0,FraHarLosOrm_so) 

! -      close input file
         close(nut)

!        open output files and write header 
         if (icrop.eq.1) then
            call outbalcropOM1(1,pathwork,outfil,project,date,daycrop,  &
     &         t,dvs,tsum,gass,mres,fr,fl,fs,fo,dmi,cvf,ccheck)
            call outbalcropOM2(1,pathwork,outfil,project,date,daycrop,  &
     &         t,dvs,tsum,storagediff,wlv,wst,wso,wrt,delt,         &
     &         grlv,grst,grso,grrt,drlv,drst,drso,drrt,ombalan)
            call outbalcropN(1,pathwork,outfil,project,date,daycrop,    &
     &         t,dvs,tsum,nuptt,nfixtt,anlvi,ansti,anrti,ansoi,anlv,&
     &         anst,anrt,anso,nlossl,nlossr,nlosss,nbalan,nni)
         endif
      endif

! --- maximum rooting depth
      if (swrd.eq.1) then
        rdm = rdmax
      elseif (swrd.eq.2) then
        rdm = min(rdmax,rdc)
      elseif (swrd.eq.3) then
        rdc = afgen (rlwtb,22,wrtmax)
        rdm = min(rdmax,rdc)
      endif

! --- skip next initialization if crop parameters are read from *.END file
      if (t1900 - tstart .gt. tiny .or. swinco .ne. 3 .or.             &
     &   dabs(t1900 - cropstart(icrop)) .lt. tiny) then
      
        dvs = 0.0d0
        flAnthesis = .false.
        tsum = 0.0d0
        fr = afgen (frtb,30,dvs)
        fl = afgen (fltb,30,dvs)
        fs = afgen (fstb,30,dvs)
        fo = afgen (fotb,30,dvs)
! --- only for bulb crops (tulips etc..)
        if(swbulb.eq.1) then
           fbl = afgen (fbltb,30,dvs)
           plwt = plwti
        endif
        sla(1) = afgen (slatb,30,dvs)
        lvage(1) = 0.0d0
        ilvold = 1
        slapot(1) = afgen (slatb,30,dvs)
        lvagepot(1) = 0.0d0
        ilvoldpot = 1

! ---   initial state variables of the crop
        wrt = fr*tdwi
        wrtmin = wrt / 10000 ! minimum root weigth at relative depth is set to 1% of the initial value
        wrtpot = wrt
        tadw = (1.0d0-fr)*tdwi
        tadwpot = tadw
        wst = fs*tadw
        wstpot = wst
        wso = fo*tadw
        wsopot = wso
        wlv = fl*tadw
        wlvpot = wlv
! --- only for bulb crops (tulips etc..)
        if(swbulb.eq.1) then
!          blad bij opkomst is ondergronds: lai vanuit ingelezen laiem,
!          sla(l) aangepast aan initieel bladgewicht en laiem,
!          stengelgewicht bij opkomst niet meegenomen bij lai-berekening
           sla(1) = laiem / wlv
           wstem = wst
           wbl = fbl*tadw
           wblpot = wbl
        else
!          KRO-BOO-20160403: intro because comparison with Wofost
           laiem = wlv*sla(1)
        endif
        lv(1) = wlv
        lvpot(1) = wlv
        lasum = laiem     
        lasumpot = laiem     
        laiexp = laiem
        laiexppot = laiem
        glaiex = 0.0d0
        glaiexpot = 0.0d0
        laimax = laiem
! --- only for bulb crops (tulips etc..)
        if(swbulb.eq.1) then
            lai = lasum+ssa*(wst-wstem)+spa*wso 
            dwbl = 0.0d0
            dwblpot = 0.0d0
        else
            lai = lasum+ssa*wst+spa*wso 
        endif
        laipot = lai 
        dwrt = 0.0d0
        dwrtpot = 0.0d0
        dwlv = 0.0d0
        dwlvCrop = 0.0d0
        dwlvSoil = 0.0d0
        dwlvpot = 0.0d0
        dwso = 0.0d0
        dwst = 0.0d0
        dwstpot = 0.0d0
        if(flCropNut) then
          WLVt0 = wlv
          WSTt0 = wst
          WSOt0 = wso
          WRTt0 = wrt
        endif
        
! ---   n-p-k 
        if( flCropNut) then
!******************************************************************
!         initial maximum nutrient concentrations in plant organs 
!         per kg biomass [kg N kg-1 dry biomass] at sowing added IS
!******************************************************************        
          call nutrsow(anlv,anst,anrt,anso)
!******************************************************************
!         initial maximum nutrient concentrations in plant organs 
!         per kg biomass [kg N kg-1 dry biomass] at emergence added IS
!******************************************************************
          call nutremrg(nmxlv,lsnr,lrnr,wlv,wst,wrt,                    &
     &      anlv,anst,anrt,anso,anlvi,ansti,anrti,ansoi,dvs)
        endif

! --- actual rooting depth
        if (swrd.eq.1) then
          rd = afgen (rdtb,22,dvs)
          rd = min(rd,rdm)
        elseif (swrd.eq.2) then
          rd = min(rdi,rdm)
        elseif (swrd.eq.3) then
          rdi = afgen (rlwtb,22,wrt)
          rd = min(rdi,rdm)
        endif
        rdpot = rd
        
! --- initial summation variables of the crop
        gasst = 0.0d0
        gasstpot = 0.0d0
        mrest = 0.0d0 
        mrestpot = 0.0d0 
        cwdm = 0.0d0
        cwdmpot = 0.0d0
! --- only for vernalisation
        vern = 0.0d0             ! vernalisation state (d)
        flvernalised = .FALSE.   ! crop not vernalised (-)

! --- end skip above initialization if crop parameters are read from *.END file
      endif
      
! --- set crop height and cropfactor
      if (swcf.ne.3) then
        cf = afgen (cftb,(2*magrs),dvs)
        ch = afgen (chtb,(2*magrs),dvs)
      else
        cf        = afgen (cftb,(2*magrs),lai)
        cfeic     = afgen (cfeictb,(2*magrs),lai)
        ch        = afgen(chtb,(2*magrs),lai)
      endif
      
! --- initial storage on canopy
      if (swinter.eq.3) then
        siccapact = siccaplai*lai
      endif

! --- initialize matric flux potential and hleaf
      if (swdrought .eq. 2) then                                        
        call MatricFlux(1,h(1),1,dummy)
        if (swhydrlift .eq. 1) then
          flhydrlift = .true.
        else
          flhydrlift = .false.
        endif
        do i = 1,numnod
         twilt(i) = watcon(i,wiltpoint)
         hroot(i) = h(i)
        enddo
        hleaf = -2000.d0
      endif                                                             

! -      n-p-k 
      if( flCropNut) then
        call nutrinit    (nlossl,nlossr,nlosss,                         &
     &                    nuptt,rnlv,rnst,rnrt,rnso,                    &
     &                    rnldlv,rnldst,rnldrt,                         &
     &                    nfixtt,nni,NLOSSLDeceasedLvToSoil)
      endif


      return

      case (2)

! === calculate potential rate and state variables =====================

! --- rates of change of the crop variables ----------------------------

! --- increase in temperature sum
      dtsum = afgen (dtsmtb,30,tav)

! --- phenological development rate for potential AND actual crops
      if (swsoybean.eq.0) then
! ---   standard crops
        if (dvs.lt.1.0d0) then     
! ---     vegetative phase
          dvred = 1.0d0
          vernfac = 1.0d0
          vernrate = 0.0d0
          if (idsl.ge.1) then
             dvred = max(0.0d0,min(1.0d0,(daylp-dlc)/(dlo-dlc)))
          endif
          if (idsl.eq.2) then
!            vernalisation rate,based on routines from pyWofost (Allard de Wit, 2015)
             if(.not.flvernalised) then
                if(dvs.lt.verndvs) then
                   vernrate = afgen (vernrtb,30,tav)
                   r = (vern - vernbase) / (vernsat - vernbase)
                   vernfac = interpol(0.0d0,1.0d0,r)
                else
                   flvernalised = .true.
                endif
             endif
          endif
          dvr = vernfac * dvred*dtsum/tsumea
        else
! ---     generative phase
          dvr = dtsum/tsumam
        endif    

      else if (swsoybean.eq.1) then
! ---   soybean
        call mgtemprf(tav,toptdvr,tmindvr,tmaxdvr,rfmgtemp)
        call mgphotoprf(mg,daynr,lat,popt,pcrt,flphenodayl,rfmgphotop)
        if (dvs.lt.1.0d0) then 
! ---     vegetative phase
          if(flrfphotoveg) then
             dvr = dvrmax1 * rfmgphotop * rfmgtemp
          else
             dvr = dvrmax1 * rfmgtemp
          endif
        else
! ---     generative phase
          dvr = dvrmax2 * rfmgphotop * rfmgtemp
        endif
      endif

!     adjust development stage for realistic TSUM
      if (dvs.ge.1.d0 .and. (.not. flAnthesis)) then
        flAnthesis = .true.
        dvs = 1.0d0
      end if
      
! === daily dry matter production 

      gasspot = pgasspot

! --- respiration and partitioning of carbohydrates between growth and
! --- maintenance respiration
      if(swbulb.eq.1) then
        rmrespot = (rmr*wrtpot+rml*wlvpot+rms*wstpot+rms*wblpot+        &
     &           rmo*wsopot)* afgen(rfsetb,30,dvs)
      else
        rmrespot = (rmr*wrtpot+rml*wlvpot+rms*wstpot+rmo*wsopot)*       &
     &           afgen(rfsetb,30,dvs)
      endif
      teff = q10**((tav-25.0d0)/10.0d0)
      mrespot = dmin1(gasspot,rmrespot*teff)
      asrcpot = gasspot - mrespot

! --- partitioning factors
      fr = afgen(frtb,30,dvs)
      fl = afgen(fltb,30,dvs)
      fs = afgen(fstb,30,dvs)
      fo = afgen(fotb,30,dvs)
! --- only for bulb crops (tulips etc..)
      if(swbulb.eq.1) then
         fbl = afgen(fbltb,30,dvs)
      endif
! --- check on partitioning
      call chckprt(dvs,fr,fl,fs,fo,fbl)    

! --- conversion factor 
      if(swbulb.eq.1) then
!       only for bulb crops (tulips etc..)
        cvf = 1.0d0/((fl/cvl+fs/cvs+fbl/cvs+fo/cvo)*(1.0d0-fr)+fr/cvr)
      else
        cvf = 1.0d0/((fl/cvl+fs/cvs+fo/cvo)*(1.0d0-fr)+fr/cvr)
      endif
      dmipot = cvf*asrcpot
! --- check on carbon balance
      call chckcbl(dvs,cvf,dmipot,fr,fl,fs,fbl,fo,gasspot,mrespot,      &
     &             ccheck)

! == = growth rate by plant organ

! --- growth rate roots and aerial parts
      admipot = (1.0d0-fr)*dmipot
      grrtpot = fr*dmipot
      ! in case of SWRD = 3: after reaching maximum live weight of wrtmax, the
      ! growth of the roots is balanced by the death of root tissue
      if (swrd.eq.3 .and. wrtpot.gt.wrtmax) then
        drrtpot = grrtpot
        drrtpot = max(drrtpot,wrtpot*afgen (rdrrtb,30,dvs))
      else  
        drrtpot = wrtpot*afgen (rdrrtb,30,dvs)
      endif  
      gwrtpot = grrtpot - drrtpot

! --- weight of new leaves
      grlvpot = fl*admipot

! --- death of leaves due to water stress or high lai
      laicr = 3.2d0/kdif
      dslvpot = wlvpot*max(0.0d0,min(0.03d0,0.03d0*                     &
     &           (laipot-laicr)/laicr))

! --- death of leaves due to exceeding life span:

! --- first: leaf death due to water stress or high lai is imposed 
! ---        on array until no more leaves have to die or all leaves
! ---        are gone

      restpot = dslvpot*delt
      i1 = ilvoldpot

      do while (restpot.gt.lvpot(max(i1,1)).and.i1.ge.1)
        restpot = restpot - lvpot(i1) 
        i1 = i1-1
      enddo

! --- then: check if some of the remaining leaves are older than span,
! ---       sum their weights

      dalvpot = 0.0d0
      if (lvagepot(max(i1,1)).gt.span .and. restpot.gt.0.0d0            &
     &                   .and.i1.ge.1) then
        dalvpot = lvpot(i1) - restpot
        restpot = 0.0d0
        i1 = i1-1
      endif

      do while (i1.ge.1.and.lvagepot(max(i1,1)).gt.span)
        dalvpot = dalvpot+lvpot(i1)
        i1 = i1-1
      enddo

      dalvpot = dalvpot/delt

! --- finally: calculate total death rate leaves
      drlvpot = dslvpot + dalvpot

! --- physiologic ageing of leaves per time step
      fysdel = max (0.0d0,(tav-tbase)/(35.0d0-tbase))

! --- specific leaf area valid for current timestep
      slatpot = afgen (slatb,30,dvs)

! --- calculation of specific leaf area in case of exponential growth:
! --- leaf area not to exceed exponential growth curve
      if (laiexppot.lt.6.0d0) then
        dteff = max (0.0d0,tav-tbase)
! ---   increase in leaf area during exponential growth
        glaiexpot = laiexppot*rgrlai*dteff
! ---   source-limited increase in leaf area
        glasolpot = grlvpot*slatpot
! ---   actual increase is determined by lowest value
        glapot = min (glaiexpot,glasolpot)
! ---   slat will be modified in case gla equals glaiex
        if (grlvpot.gt.0.0d0) slatpot = glapot/grlvpot
      endif  

! --- growth rate stems
      grstpot = fs*admipot
! --- death rate stems
      drstpot = afgen (rdrstb,30,dvs)*wstpot
! --- net growth rate stems
      gwstpot = grstpot - drstpot

! --- only for bulb crops (tulips etc..)
      if(swbulb.eq.1) then
! --    growth rate flowers
        grblpot = fbl*admipot
        if(dvs.ge.1.0d0) then
          grblpot = 0.0d0
          drblpot = wblpot/delt
        endif
        gwblpot = grblpot - drblpot
      endif

! --- growth rate storage organs
      gwsopot = fo*admipot

! ----integrals of the crop --------------------------------------------

! --- leaf death (due to water stress or high lai) is imposed on array 
! --- untill no more leaves have to die or all leaves are gone

      dslvtpot = dslvpot*delt
      i1 = ilvoldpot
      do while (dslvtpot.gt.0.and.i1.ge.1)
        if (dslvtpot.ge.lvpot(i1)) then
          dslvtpot = dslvtpot-lvpot(i1)
          lvpot(i1) = 0.0d0
          i1 = i1-1
        else
          lvpot(i1) = lvpot(i1)-dslvtpot
          dslvtpot = 0.0d0
        endif
      enddo

! --- leaves older than span die
      do while (lvagepot(max(i1,1)).gt.span.and.i1.ge.1)
        lvpot(i1) = 0.0d0
        i1 = i1-1
      enddo

! --- oldest class with leaves
      ilvoldpot = i1

! --- shifting of contents, updating of physiological age
      do i1 = ilvoldpot,1,-1
        lvpot(i1+1) = lvpot(i1)
        slapot(i1+1) = slapot(i1)
        lvagepot(i1+1) = lvagepot(i1)+fysdel*delt
      enddo
      ilvoldpot = ilvoldpot + 1

! --- new leaves in class 1
      lvpot(1) = grlvpot*delt
      slapot(1) = slatpot
      lvagepot(1) = 0.0d0 

! --- calculation of new leaf area and weight
      lasumpot = 0.0d0
      wlvpot = 0.0d0
      do i1 = 1,ilvoldpot
        lasumpot = lasumpot + lvpot(i1)*slapot(i1)
        wlvpot = wlvpot + lvpot(i1)
      enddo

! --- leaf area index in case of exponential growth
      laiexppot = laiexppot+glaiexpot*delt

! --- dry weight of living plant organs
      wrtpot = wrtpot + gwrtpot*delt
      wstpot = wstpot + gwstpot*delt
      wsopot = wsopot + gwsopot*delt
! --- only for bulb crops (tulips etc..)
      if(swbulb.eq.1) then
        wblpot = wblpot + gwblpot*delt
      endif

! --- total above ground biomass
      tadwpot = wlvpot + wstpot + wsopot
! --- only for bulb crops (tulips etc..)
      if(swbulb.eq.1) then
         tadwpot = tadwpot + wblpot
      endif

! --- dry weight of dead plant organs (roots,leaves & stems)
      dwrtpot = dwrtpot + drrtpot*delt
      dwlvpot = dwlvpot + drlvpot*delt
      dwstpot = dwstpot + drstpot*delt
! --- only for bulb crops (tulips etc..)
      if(swbulb.eq.1) then
        dwblpot = dwblpot + drblpot*delt
      endif

! --- dry weight of dead and living plant organs
      twlvpot = wlvpot + dwlvpot
      twstpot = wstpot + dwstpot
      cwdmpot = twlvpot + twstpot + wsopot
! --- only for bulb crops (tulips etc..)
      if(swbulb.eq.1) then
        twblpot = wblpot + dwblpot
        cwdmpot = cwdmpot + twblpot
      endif

! --- total gross assimilation and maintenance respiration
      gasstpot = gasspot + gasstpot
      mrestpot = mrespot + mrestpot

! --- leaf area index
      laipot = lasumpot + ssa*wstpot + spa*wsopot
!     prevent immediate lai reduction at emergence
!     KRO-BOO-20160403: suppressed because deviates from Wofost
!      laipot = max(laipot, laiem)

! --- vernalisation state (d)
      if(idsl.eq.2) then
          vern = vern + vernrate
          if(.not.flvernalised .and. vern.ge.vernsat) then
              flvernalised = .true.
          else
              if(flvernalised .and. vern.lt.vernsat) then
                 write(messag,'(2a,i6,2a)') ' critical DVS,',           &
     &           ' for vernalised reached, day = ', daycum,             &
     &           ' but vernalisation requirements not yet fulfilled ',  &
     &           ' forcing vernalization now'
                 call warn ('wofost',messag,logf,swscre)
              endif
          endif
      endif

      return

      case (3)

! === calculate actual rate and state variables =====================
! === with optional calculation of (flCropNut) water AND nutrient stress

! --- rates of change of the crop variables ----------------------------
 
! --- water stress reduction of pgass to gass
      if(dabs(ptra).lt.nihil) then
        reltr = 1.0d0
      else
        reltr = max(0.0d0,min(1.0d0,tra/ptra))
      endif

! --- nitrogen stress reduction of pgass to gass
      if (flCropNut) then
        reltr = min(reltr,fstr)
        fstr  = reltr
      end if
      gass = pgass * reltr

! --- respiration and partitioning of carbohydrates between growth and
! --- maintenance respiration
! --  only for bulb crops (tulips etc..)
      if(swbulb.eq.1) then
        rmres = (rmr*wrt+rml*wlv+rms*wst+rms*wbl+rmo*wso)*              &
     &           afgen(rfsetb,30,dvs)
      else
        rmres = (rmr*wrt+rml*wlv+rms*wst+rmo*wso)*afgen(rfsetb,30,dvs)
      endif
      
      mres = dmin1(gass,rmres*teff)
      asrc = gass-mres

! --- partitioning factors
      fr = afgen(frtb,30,dvs)
      fl = afgen(fltb,30,dvs)
      fs = afgen(fstb,30,dvs)
      fo = afgen(fotb,30,dvs)
! --- only for bulb crops (tulips etc..)
      if(swbulb.eq.1) then
         fbl = afgen(fbltb,30,dvs)
      endif
! --- check on partitioning
      call chckprt(dvs,fr,fl,fs,fo,fbl)    

      if( flCropNut) then
!********************************************************************         
!         partitioning correction as influenced by water and N stress
!         Note: the partioning depends only on the Nitrogen stress,
!         not on the P and K stress. Personal communication Joost Wolf        
!         added IS
!******************************************************************** 
          CALL SUBPAR (reltr,NPART,NNI,FR,FL,FS,FO)
      endif
      
! --- conversion factor 
      if(swbulb.eq.1) then
!       only for bulb crops (tulips etc..)
        cvf = 1.0d0/((fl/cvl+fs/cvs+fbl/cvs+fo/cvo)*(1.0d0-fr)+fr/cvr)
      else
        cvf = 1.0d0/((fl/cvl+fs/cvs+fo/cvo)*(1.0d0-fr)+fr/cvr)
      endif
! --- dry matter increase
      dmi = cvf*asrc
! --- check on carbon balance
      call chckcbl(dvs,cvf,dmi,fr,fl,fs,fbl,fo,gass,mres,ccheck)

! --- growth rate by plant organ

! --- growth rate roots and aerial parts
      call relgrwt(dmi,fr,fl,fs,fo,grrt,grlv,grst,grso,admi)
      if (swrd.eq.3 .and. flWrtNonox) grrt = 0.d0
      
! --- death of leaves due to water stress or high lai or nitrogen stress
      call deaths(flcropnut,wlv,kdif,lai,NNI,perdl,rdrns,reltr,dslv)

! --- death of leaves due to exceeding life span:
      call deatha(dslv,delt,ilvold,lv,lvage,span,i1,dalv)

! --- death rate leaves as result of death due to water stress or high lai and 
!                                    death due to exceeding life span
      drlv = dslv+dalv

! --- death rate stems
      drst = wst * afgen (rdrstb,30,dvs)

! --- death rate roots
      ! in case of SWRD = 3: after reaching maximum live weight of wrtmax, the
      ! growth of the roots is balanced by the death of root tissue
      if (swrd.eq.3 .and. wrt.gt.wrtmax) then
        drrt = grrt
        drrt = max(drrt,wrt*afgen (rdrrtb,30,dvs))
      else  
        drrt = wrt*afgen (rdrrtb,30,dvs)
      endif  

! --- net growth rate stems, roots, storage organs
      gwst = grst - drst
      gwrt = grrt - drrt
      drso = 0.0d0    ! death rate of storage organs is assumed to be 0
      gwso = grso - drso

! --- only for bulb crops (tulips etc..)
      if(swbulb.eq.1) then
! --    growth rate flowers
        grbl = fbl*admi
        if(dvs.ge.1.0d0) then
          grbl = 0.0d0
          drbl = wbl/delt
        endif
        gwbl = grbl - drbl
      endif

! --- specific leaf area valid for current timestep
      if(flCropNut) then
!       nutrient and water stress
        slat = afgen (slatb,30,dvs)*EXP(-NSLA * (1.0d0-NNI))
      else
        slat = afgen (slatb,30,dvs)
      endif
!
!     Do not allow slat higher than slatpot; slatpot can be limited by exponential growth
      slat = min(slat,slatpot)   ! pvw

! --- calculation of specific leaf area in case of exponential growth:
! --- leaf area not to exceed exponential growth curve
! --- FSTR is actual stress: water and nutrient 
      Fstress = reltr
      if (flcropnut) then
         Fstress = FSTR
         if ((DVS .LT. 0.2d0).AND.(LAI .LT. 0.75d0)) then
           Fstress = reltr * EXP(-NLAI* (1.0d0 - NNI))
         endif
      endif
      call GLAI(Fstress,LAIEXP,GLAIEX,tav,TBASE,RGRLAI,GRLV,SLAT,GLA)


! ---- UPDATE STATES: integrals of the crop --------------------------------------------

! --- phenological development stage
      dvs = min(dvs+dvr*delt,dvsend)
      tsum = tsum + dtsum*delt

! --- leaf death (due to water stress or high lai) is imposed on array 
! --- untill no more leaves have to die or all leaves are gone
      call lvdth(delt,dslv,span,ilvold,lvage,lv,i1)

! --- oldest class with leaves
      ilvold = i1

! --- shifting of contents, updating of physiological age
      call lvshft(delt,fysdel,ilvold,grlv,slat,lv,lvage,sla)

      ilvold = ilvold+1

! --- calculation of new leaf area and weight
      call lvwgli(ilvold,lv,sla,lasum,wlv)

! --- leaf area index in case of exponential growth
      laiexp = laiexp+glaiex*delt
      
! --- dry weight of living plant organs
      wrt = wrt+gwrt*delt
      wst = wst+gwst*delt
      wso = wso+gwso*delt
! --- only for bulb crops (tulips etc..)
      if(swbulb.eq.1) then
        wbl = wbl + gwbl*delt
      endif

! --- total above ground biomass
      tadw = wlv+wst+wso
! --- only for bulb crops (tulips etc..)
      if(swbulb.eq.1) then
        tadw = tadw + wbl
      endif

! --- dry weight of dead plant organs (roots,leaves & stems)
      dwrt = dwrt + drrt*delt
      dwlv = dwlv + drlv*delt
      dwst = dwst + drst*delt
      dwso = dwso + drso*delt   ! dummy, because drso is assumed to be 0
! --- only for bulb crops (tulips etc..)
      if(swbulb.eq.1) then
        dwbl = dwbl + drbl*delt
      endif

!     split dwlv
      idwlvCrop = (1.0d0-FraDeceasedLvToSoil) * drlv*delt
      idwlvSoil = FraDeceasedLvToSoil * drlv*delt
      dwlvCrop = dwlvCrop + idwlvCrop
      dwlvSoil = dwlvSoil + idwlvSoil
      
! --- dry weight of dead and living plant organs
!     twrt = wrt+dwrt
      twlv = wlv+dwlv
      twst = wst+dwst
      cwdm = twlv+twst+wso
! --- only for bulb crops (tulips etc..)
      if(swbulb.eq.1) then
        twbl = wbl + dwbl
        cwdm = cwdm + twbl
      endif

! --- total gross assimilation and maintenance respiration
      gasst = gass + gasst
      mrest = mres + mrest

! --- leaf area index
      lai = lasum+ssa*wst+spa*wso
! --- determine maximum lai
      laimax = max (lai,laimax)
! --- determine minimum lai to prevent dying straight after 
!       emergence when growth is slowed down due to low temperature
!     KRO-BOO-20160403: suppressed because deviates from Wofost
!      lai = max(lai, laiem)

! --- update state variables for nutrient stress

!     Calling the subroutine for N losses of leaves, roots and stem storage
!     organs (kg N ha-1 d-1)
      if(flCropNut)then 

!        Calling the subroutines for N demand of leaves, roots and stem storage
!        organs (kg N ha-1 d-1)
         CALL NDEMND(WLV,WST,WRT,WSO,NMAXLV,NMAXST,                     &
     &                               NMAXRT,NMAXSO,ANLV,ANST,ANRT,ANSO, &
     &                               TCNT,NDEML,NDEMS,NDEMR,NDEMSO)

!        Total N demand (kg N ha-1)

         NDEMTO = MAX (0.0d0,(NDEML + NDEMS + NDEMR))

!        Nutrient uptake limiting factor (-) at low moisture conditions in the
!        rooted soil layer before anthesis. After anthesis/DVSNLT there is no
!        nutrient uptake from the soil
         NLIMIT = INSW(DVS-DVSNLT,INSW(reltr-0.01d0,0.0d0,1.0d0),0.d0)
         NdemandSoil = (1.d0-NFIXF) * NDEMTO * NLIMIT
         NdemandBioFix =  NFIXF * NDEMTO * NLIMIT

      end if
         
      return

      case (4)

      if(flCropNut) then

!        Total N uptake (kg N ha-1 d-1) from soil and by biological N fixation         
         NUPTR = (MAX(0.d0, MIN(NdemandSoil, NsupplySoil) ))/DELT
         NFIXTR = (MAX(0.d0, NdemandBioFix))/DELT

!        Calling the subroutine to estimate the translocatable nutrients in leaves, stem, roots and
!        storage organs (kg N ha-1)
         CALL NTRLOC(ANLV,ANST,ANRT,WLV,WST,WRT,RNFLV,RNFST,RNFRT,      &
     &                  FNTRT,ATNLV,ATNST,ATNRT,ATN)

!        N supply to the storage organs (kg N ha-1 d-1)      
         NSUPSO = INSW (DVS-DVSNT,0.0d0,ATN/TCNT)

!        Rate of N uptake in grains (kg N ha-1 d-1)
         RNSO =  MIN (NDEMSO,NSUPSO)

!        Calling the subroutine to calculate nutrient translocation from leaves, stem, and roots (kg N ha-1 d-1)
         CALL NTRANS(RNSO,ATNLV,ATNST,ATNRT,ATN,RNTLV,RNTST,RNTRT)

!        Calling the subroutine to compute the partitioning of the total
!        nutrient uptake rate (NUPTR) over the leaves, stem and roots (kg N ha-1 d-1)
         CALL RNUSUB(NDEML,NDEMS,NDEMR,NUPTR,                           &
     &                  NFIXTR,NDEMTO, RNULV,RNUST,RNURT)
        
!        Calling routine to calculate nutrient losses due to dying leaves, stems           
!        and roots (kg N ha-1 d-1)    
         CALL RNLD(DRLV,DRRT,DRST,RNFLV,RNFRT,RNFST,                    &
     &                RNLDLV,RNLDRT,RNLDST)


!   - ---Rate of change of N in crop organs   
         RNLV = RNULV - RNTLV - RNLDLV
         RNST = RNUST - RNTST - RNLDST
         RNRT = RNURT - RNTRT - RNLDRT

!-----   Total N  uptake by crop over time (kg N ha-1) from soil and by biological fixation
         NUPTT = NUPTT + NUPTR*DELT
         NFIXTT= NFIXTT+ NFIXTR*DELT

!-----   Actual N amount in various living organs and total living N amount(kg N ha-1)
         ANLV =  max(0.0d0, (ANLV + RNLV*DELT) )
         ANST =  max(0.0d0, (ANST + RNST*DELT) )
         ANRT =  max(0.0d0, (ANRT + RNRT*DELT) )
         ANSO =  ANSO + RNSO*DELT
!!!         NLIVT=  ANLV + ANST + ANRT + ANSO

!-----   N losses from leaves, roots and stems due to senescence and total N loss (kg N ha-1)
         NLOSSL =  NLOSSL + RNLDLV*DELT
         NLOSSR =  NLOSSR + RNLDRT*DELT
         NLOSSS =  NLOSSS + RNLDST*DELT
!!!         NLOSST =  NLOSSL + NLOSSR + NLOSSS

!----    total N  in living and dead roots
!!!         NROOT= ANRT + NLOSSR

!       increment values of dead weight of plant organs,
!       to be used in the soil nutrient submodel
         idwrt = drrt*delt
         idwlv = idwlvSoil
         idwst = 0.0d0
         idwso = 0.0d0
         iNLOSSR =  RNLDRT*DELT
         iNLOSSL =  FraDeceasedLvToSoil * RNLDLV*DELT
         NLOSSLDeceasedLvToSoil =  NLOSSLDeceasedLvToSoil + iNLOSSL
         NLOSSL =  NLOSSL - iNLOSSL
         iNLOSSS =  0.0d0
         iNLOSSO =  0.0d0
         HarLosOrm_rt = 0.0d0; HarLosOrm_lv = 0.0d0; HarLosOrm_st = 0.0d0
         HarLosOrm_dwlv = 0.0d0; HarLosOrm_dwst = 0.0d0
         HarLosOrm_so = 0.0d0; HarLosOrm_tot = 0.0d0 
         HarLosNit_rt = 0.0d0; HarLosNit_lv = 0.0d0 
         HarLosNit_st = 0.0d0; HarLosNit_so = 0.0d0 
!!         HarLosNit_dwrt = 0.0d0; HarLosOrm_dwrt = 0.0d0
         HarLosNit_dwst = 0.0d0; HarLosNit_dwso = 0.0d0; HarLosNit_dwlv = 0.0d0 
!        during the last day of the crop period: add the weight of living roots 
!        to the dead roots and reset living weight to zero
         if (flHarvestDay .or. (dvs.ge.dvsend) .or.                     &
     &                 dabs(t1900-1.0d0-cropend(icrop)).lt.1.0d-3 ) then
            HarLosOrm_rt = wrt
            HarLosOrm_dwlv =  FraHarLosOrm_lv * dwlv
            HarLosOrm_lv   = FraHarLosOrm_lv * wlv + HarLosOrm_dwlv
            HarLosOrm_dwst =  FraHarLosOrm_st * dwst
            HarLosOrm_st   = FraHarLosOrm_st * wst + HarLosOrm_dwst
            HarLosOrm_dwso =  FraHarLosOrm_so * dwso
            HarLosOrm_so   = FraHarLosOrm_so * wso + HarLosOrm_dwso
            HarLosOrm_tot = HarLosOrm_rt + FraHarLosOrm_lv * wlv +      & 
     &             FraHarLosOrm_st * wst + FraHarLosOrm_so * wso
!ckro_sup_20170714 : suppressed because it will happen after harvest
!            wrt = wrt - HarLosOrm_rt
!            wlv = wlv - FraHarLosOrm_lv * wlv
!            wst = wst - FraHarLosOrm_st * wst
!            wso = wso - FraHarLosOrm_so * wso
            idwrt = idwrt + HarLosOrm_rt
            idwlv = idwlv + HarLosOrm_lv
            idwst = idwst + HarLosOrm_st
            idwso = idwso + HarLosOrm_so
            HarLosNit_rt = ANRT
            HarLosNit_dwlv = FraHarLosOrm_lv * NLOSSL 
            HarLosNit_lv = FraHarLosOrm_lv * ANLV + HarLosNit_dwlv 
            HarLosNit_dwst = FraHarLosOrm_st * NLOSSS 
            HarLosNit_st = FraHarLosOrm_st * ANST + HarLosNit_dwst 
            HarLosNit_dwso = FraHarLosOrm_so * 0.0d0 
            HarLosNit_so = FraHarLosOrm_so * ANSO + HarLosNit_dwso 
            iNLOSSL = iNLOSSL + HarLosNit_lv
            iNLOSSS = iNLOSSS + HarLosNit_st
            iNLOSSO = iNLOSSO + HarLosNit_so
            iNLOSSR = iNLOSSR + HarLosNit_rt
            ANLV = ANLV - FraHarLosOrm_lv * ANLV
            ANST = ANST - FraHarLosOrm_st * ANST
            ANSO = ANSO - FraHarLosOrm_so * ANSO
            ANRT = ANRT - HarLosNit_rt
        end if

! ----- CHECK and WRITE MASS BALANCE: dry matter of crop

!       output of OM balance1: from air to partitioning (kg/ha DM CH2O)
        call outbalcropOM1(2,pathwork,outfil,project,date,daycrop,      &
     &       t,dvs,tsum,gass,mres,fr,fl,fs,fo,dmi,cvf,ccheck)

! -     OM balance2: storage difference(kg/ha DM CH2O)
        storagediff = (wlv+wst+wso+wrt) - (wlvt0+wstt0+wsot0+wrtt0)
        ombalan = storagediff - ( (grlv+grst+grso+grrt)*delt -          &
     &            (drlv+drst+drso+drrt)*delt ) 
!ckro_20171002 harvest losses happen after harvest and should not be 
!    in this balance, therefore this commented out of source code.
!        storagediff = HarLosOrm_rt + HarLosOrm_lv + HarLosOrm_st +      &
!     &  HarLosOrm_so -(HarLosOrm_dwlv + HarLosOrm_dwst + HarLosOrm_dwso)
        if (dabs(ombalan) .ge. 1.0d0) then
           write(messag,'(a,f8.3)')                                     &
     &      ' Warning Wofost: OM balance2 not 0, OMBAL = ',ombalan
           call warn ('wofost',messag,logf,swscre)
!     &     ' OM balance2 not 0, simulation stopped OMBAL=',ombalan
!           call fatalerr ('wofost',messag)
        endif
!       output of OM balance2
        call outbalcropom2(2,pathwork,outfil,project,date,daycrop,      &
     &         t,dvs,tsum,storagediff,wlv,wst,wso,wrt,delt,             &
     &         grlv,grst,grso,grrt,drlv,drst,drso,drrt,ombalan)

! ----- CHECK and WRITE MASS BALANCE: nitrogen of crop  

        NBALAN =  NUPTT + NFIXTT + (ANLVI+ANSTI+ANRTI+ANSOI)            &
     &      - (ANLV+ANST+ANRT+ANSO) - (NLOSSL+NLOSSR+NLOSSS)            &
     &      - NLOSSLDeceasedLvToSoil                                    &
     &      - (HarLosNit_lv+HarLosNit_st+HarLosNit_so+HarLosNit_rt)     &
     &      +  HarLosNit_dwlv + HarLosNit_dwst + HarLosNit_dwso

!       output of N balance
        call outbalcropN(2,pathwork,outfil,project,date,daycrop,        &
     &         t,dvs,tsum,NUPTT,NFIXTT,ANLVI,ANSTI,ANRTI,ANSOI,ANLV,&
     &         ANST,ANRT,ANSO,NLOSSL,NLOSSR,NLOSSS,NBALAN,nni)
        IF (dabs(NBALAN) .GE. 1.0d-03) then
           write(messag,'(1a,i6,a,f8.3)') ' Nitrogen balance not 0,'//  &
     &     ' simulation stopped, day = ', daycum,' NBAL=',NBALAN
!*           call fatalerr ('wofost',messag)
           call warn ('CropGrowth_Wofost',messag,logf,swscre)
        endif
 
        if (flHarvestDay .or. (dvs.ge.dvsend) .or.                      &
     &                 dabs(t1900-1.0d0-cropend(icrop)).lt.1.0d-3 ) then
          gwst  = 0.0d0
          gwrt  = 0.0d0
          gwso  = 0.0d0
          grlv  = 0.0d0
          NdemandSoil = 0.0d0
          HarLosOrm_tot = 0.0d0
        endif
      endif

! --- update normalized cumulative root density based on root extraction or stress (cumdens)
      if (swrdc.eq.1) call update_rootdistribution()
      
! --- root extension
      if (swrd.eq.1) then
      
        rdpot = afgen (rdtb,22,dvs)
        rdpot = min(rdpot,rdm)
        rd    = rdpot
      
      elseif (swrd.eq.2) then
        
        rrpot = min (rdm-rdpot,rri)
        if (fr.le.0.0d0 .or. pgasspot.lt.1.0d0) rrpot = 0.0d0
        rdpot = rdpot + rrpot

        rr = min (rdm-rd,rri)
        if (fr.le.0.0d0 .or. pgass.lt.1.0d0 .or. flWrtNonox) rr = 0.0d0
        if (swdmi2rd.eq.1 .and. pgass.ge.1.0d0)              rr = rr * gass/pgass
        rd = rd + rr
        
      elseif (swrd.eq.3) then
        rdpot = afgen (rlwtb,22,wrtpot)
        rdpot = min(rdpot,rdm)
        rd = afgen (rlwtb,22,wrt)
        rd = min(rd,rdm)
      endif

! --- crop factor or crop height
      if (swcf.ne.3) then
        cf = afgen (cftb,(2*magrs),dvs)
        ch = afgen (chtb,(2*magrs),dvs)
      else
        cf        = afgen (cftb,72,lai)
        cfeic     = afgen (cfeictb,72,lai)
        ch        = afgen(chtb,72,lai)
      endif

! --- update canopy storage capacity
      if (swinter.eq.3) then
        siccapact = siccaplai*lai
      endif
      
! --- update states of dry matter organs
      wlvt0 = wlv
      wstt0 = wst
      wsot0 = wso
      wrtt0 = wrt

      case default
         call fatalerr ('Wofost', 'Illegal value for TASK')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine grass(task)
! ----------------------------------------------------------------------
!     Date               : November 2004
!     Purpose            : detailed grass growth routine 
! ----------------------------------------------------------------------
      use variables
      implicit none
 
      integer   i1,task
      integer   idelaypot,idelay,i,swhydrlift

      real(8)   laicr,lasum,mres,grazlivinglv,grazlivinglvpot
      real(8)   admi,afgen,asrc,ccheck,cvf
      real(8)   dalv,delt,dmi
      real(8)   drrt,drst,dslv,dslv1,dslv2,dslvt,dteff
      real(8)   fcheck,fl,fr,fs,drlv
      real(8)   fysdel,gass,gla,glasol,grlv,grrt,grst
      real(8)   gwst,rest,rmres
      real(8)   slat,teff,twlv,twst
      real(8)   lasumpot,drst1,drst2
      real(8)   drst1pot,drst2pot
      real(8)   gasspot,rmrespot,mrespot,asrcpot,dmipot
      real(8)   admipot,grrtpot,drrtpot,gwrtpot,grlvpot,dslv1pot
      real(8)   dslv2pot,dslvpot,restpot,dalvpot,drlvpot
      real(8)   glasolpot,slatpot,glapot,grstpot,drstpot,gwstpot
      real(8)   dslvtpot,twlvpot,twstpot,tagpspot,tagps
      real(8)   dummy
      real(8)   dmharvest,dmlastharvest,dmgrazing
      real(8)   lsdb(100)
      real(8)   watcon,uptgraz,tagprest,lossgraz
      real(8)   uptgrazpot,lossgrazpot
      integer   daylastharvest,swharvest,daysgrazpot,daysgraz
      integer   swdmmow,swdmgrz
      character(len=11) tmp
      character(len=200) messag
      
! --- rooting
      real(8)   rrpot,rr
      
      integer   maxdaymow,maxdaygrz

      logical   flGrassGrowth
      logical   flDewoolingpot,flDewooling
      character(len=11) ::  dateGrassGrowth
      
!     In case of swlossmow = 1 or swlossgrz = 1: check work-ablility
      integer   nodmow,nodgrz
      integer   swlossmow,swlossgrz
      real(8)   fralossmow
      real(8)   fralossgrz
      real(8)   drz1
      logical   flearlyhrvendpot,flearlyhrvendact
      
      parameter (delt=1.0d0)
      include 'params.fi'

      save
! ----------------------------------------------------------------------
      select case (task)
      case (1)

! === initialization at start of crop =========================================

! --- read grass input data
      call readgrass (icrop,cropfil(icrop),swharvest,dmharvest,daylastharvest,dmlastharvest,swdmmow,maxdaymow, &
                      swlossmow,swlossgrz,swdmgrz,maxdaygrz,dmgrazing,LSDb,tagprest,swhydrlift)

! --- sequence of harvest by mowing, dewooling and grazing
      seqgrazmowpot = seqgrazmow

! --- development stage (not used by Grassland, instead Daynrs are used)
      dvs = -99.99d0

! --- maximum rooting depth
      if (swrd.eq.1) then
        rdm = rdmax
      elseif (swrd.eq.2) then
        rdm = min(rdmax,rdc)
      elseif (swrd.eq.3) then
        rdc = afgen (rlwtb,22,wrtmax)
        rdm = min(rdmax,rdc)
      endif

! --- skip next initialization if crop parameters are read from *.END file
      if (t1900 - tstart .gt. tiny .or. swinco .ne. 3 .or.             &
     &   dabs(t1900 - cropstart(icrop)) .lt. tiny) then

        iseqgm = 1
        iseqgmpot = iseqgm

! ---   initial values of crop parameters
        rid = dble(daycrop)
        fr = afgen (frtb,30,rid)
        fl = afgen (fltb,30,rid)
        fs = afgen (fstb,30,rid)
        sla(1) = afgen (slatb,30,rid)
        lvage(1) = 0.d0
        ilvold = 1
        idregr = 0
        slapot(1) = afgen (slatb,30,rid)
        lvagepot(1) = 0.d0
        ilvoldpot = 1
        idregrpot = 0

! ---   initial state variables of the crop
        wrt = fr*tdwi
        wrtmin = wrt / 10000 ! minimum root weigth at relative depth is set to 1% of the initial value
        wrtpot = wrt
        wst = fs*(1.0d0-fr)*tdwi
        wstpot = wst
        wlv = laiem/sla(1)
        wlvpot = wlv
        
!     KRO-BOO-20160403: intro because comparison with Wofost
        laiem = wlv*sla(1)  ! is not input !
        lv(1) = wlv
        lvpot(1) = lv(1)
        lasum = laiem
        lasumpot = lasum     
        glaiex = 0.0d0
        glaiexpot = 0.0d0
        laiexp = laiem
        laiexppot = laiem
        laimax = laiem
        lai = lasum+ssa*wst
        laipot = lai
        dwrt = 0.d0
        dwrtpot = dwrt
        dwlv = 0.d0
        dwlvpot = dwlv
        dwst = 0.d0
        dwstpot = dwst

        daygrowth    = 0
        daygrowthpot = 0

! ---   actual rooting depth
        if (swrd.eq.1) then
          rd = afgen (rdtb,22,rid)
          rd = min(rd,rdm)
        elseif (swrd.eq.2) then
          rd = min(rdi,rdm)
        elseif (swrd.eq.3) then
          rdi = afgen (rlwtb,22,wrt)
          rd = min(rdi,rdm)
        endif
        rdpot = rd
        
! ---   initial summation variables of the crop
        tagp = wlv+wst
        tagppot = tagp
        tagpt = 0.0d0
        tagptpot = 0.0d0
        cuptgraz = 0.0d0
        cuptgrazpot = 0.0d0
        tsum = 0.0d0
        
        cropstartpot     = rid
        cropstartact     = rid
        flhrvendpot      = .false.
        flearlyhrvendpot = .false.
        
        if (swtsum.eq.0) then
          flGrassGrowth = .true.
        else
          flGrassGrowth = .false.  
        endif
        if (swtsum.eq.2) then
          call sumttd('initial',flGrassGrowth,dateGrassGrowth)
        endif

! --- end skip above initialization if crop parameters are read from *.END file
      endif

      if (swcf.ne.3) then
        cf = afgen (cftb,(2*magrs),rid)
        ch = afgen (chtb,(2*magrs),rid)
      else
        cf        = afgen (cftb,(2*magrs),lai)
        cfeic     = afgen (cfeictb,(2*magrs),lai)
        ch        = afgen(chtb,(2*magrs),lai)
      endif

! --- initial storage on canopy
      if (swinter.eq.3) then
        siccapact = siccaplai*lai
      endif
      
! --- initialize matric flux potential and hleaf
      if (swdrought .eq. 2) then                                        
        call MatricFlux(1,h(1),1,dummy)
        if (swhydrlift .eq. 1) then
          flhydrlift = .true.
        else
          flhydrlift = .false.
        endif
        do i = 1,numnod
         twilt(i) = watcon(i,wiltpoint)
         hroot(i) = h(i)
        enddo
        hleaf = -2000.d0
      endif                                                             

! --- harvest
!     initialise 
      if (swharvest.eq.2) then
        iharvest = 1
        do while (t1900 .gt. dateharvest(iharvest))
          iharvest = iharvest + 1
        enddo
      endif      
      
! --- Find node for monitoring work-ability
      if (swlossmow .eq. 1) then
         
        ! Find node for monitoring work-ability during mowing
        nodmow = 1
        drz1       = -1.d0 * zmow - dz(nodmow)
        do while (drz1 .gt. 0.d0)
          nodmow   = nodmow + 1
          drz1 = drz1 - dz(nodmow)
        enddo
      
      endif   
      
      if (swlossgrz .eq. 1) then
        
        ! Find node and layer for monitoring work-ability at start of grazing
        nodgrz = 1
        drz1       = -1.d0 * zgrz - dz(nodgrz)
        do while (drz1 .gt. 0.d0)
          nodgrz   = nodgrz + 1
          drz1 = drz1 - dz(nodgrz)
        enddo
         
      endif
      
      return

      case (2)

! === calculate potential rate and state variables ======================================

! --- rates of change of the grass variables ---------------------------------------------

      rid = dble(daycrop)
      
! --- check end of harvest
      if (flhrvendpot) then
        if (flearlyhrvendpot) then
          cropstartpot  = rid - 1.d0
        else
          cropstartpot  = rid
        endif
        pmowdm        = 0.d0
        pgrzdm        = 0.d0
        plossdm       = 0.d0
      endif
      flhrvendpot      = .false.
      flearlyhrvendpot = .false.
      
! --- grass growth initiated by tsum from 1st day of calendar year
      tsum = tsum + max(0.0d0,tav)
      if (.not. flGrassGrowth) then
        
        ! grass growth initiated by tsum
        if (swtsum.eq.1) then
          if (tsum.ge.200.d0) then
            flGrassGrowth = .true.
          endif
        endif
        
        ! grass growth initiated by temperature, time and depth
        if (swtsum.eq.2) then
          if (dateGrassGrowth.eq.'undefined') call sumttd('dynamic',flGrassGrowth,dateGrassGrowth)
        endif
      
        ! check if grass growth has started
        if (flGrassGrowth) then
          cropstartpot = rid
          cropstartact = rid
        endif

      endif
      
! --- skip in case of: tsum<tsum200, or 3 criteria (tsummttd), or regrowth
      if (flGrassGrowth .and. daycrop.ge.idregrpot) then

! ===   daily dry matter production ===
        
        gasspot = pgasspot
        
! ---   respiration and partitioning of carbohydrates between growth and
! ---   maintenance respiration
        rmrespot=(rmr*wrtpot+rml*wlvpot+rms*wstpot)*afgen(rfsetb,30,rid)
        teff = q10**((tav-25.0d0)/10.0d0)
        mrespot = min (gasspot,rmrespot*teff)
        asrcpot = gasspot-mrespot

! ---   partitioning factors
        fr = afgen(frtb,30,rid)
        fl = afgen(fltb,30,rid)
        fs = afgen(fstb,30,rid)
! ---   check on partitioning
        fcheck = fr+(fl+fs)*(1.0d0-fr) - 1.0d0
        if (dabs(fcheck).gt.0.0001d0) then
          write(tmp,'(f6.3)') rid
          tmp = adjustl (tmp)
          Messag ='The sum of partitioning factors for leaves, stems'// &
     &    ' and storage organs is not equal to one at time '            &
     &    //trim(tmp)//'.'
          call fatalerr ('grass_pot',messag)
        endif

! ---   dry matter increase
        cvf = 1.0d0/((fl/cvl+fs/cvs)*(1.0d0-fr)+fr/cvr)
        dmipot = cvf*asrcpot

! ---   check on carbon balance
        ccheck = (gasspot-mrespot-(fr+(fl+fs)*(1.0d0-fr))*dmipot/cvf)   &
     &         /max(0.0001d0,gasspot)      
        if (dabs(ccheck).gt.0.0001d0) then
          Messag ='The carbon balance is not correct'
          call fatalerr ('grass_pot',messag)
        endif


! ===   growth rate by plant organ ===

! ---   growth rate roots and aerial parts

        grrtpot = fr*dmipot
        ! in case of SWRD = 3: after reaching maximum live weight of wrtmax, the
        ! growth of the roots is balanced by the death of root tissue
        if (swrd.eq.3 .and. wrtpot.gt.wrtmax) then
          drrtpot = grrtpot
          drrtpot = max(drrtpot,wrtpot*afgen (rdrrtb,30,rid))
        else  
          drrtpot = wrtpot*afgen (rdrrtb,30,rid)
        endif  
        gwrtpot = grrtpot - drrtpot

! ---   growth rate leaves

! ---   weight of new leaves
        admipot = (1.0d0-fr)*dmipot
        grlvpot = fl*admipot

! ---   death of leaves due to water stress or high lai
        dslv1pot = 0.0d0
        laicr = 3.2d0/kdif
        dslv2pot=wlvpot*max(0.0d0,                                      &
     &                  min(0.03d0,0.03d0*(laipot-laicr)/laicr))
        dslvpot = max (dslv1pot,dslv2pot) 

! ---   death of leaves due to exceeding life span;
! ---   leaf death is imposed on array until no more leaves have
! ---   to die or all leaves are gone

        restpot = dslvpot*delt
        i1 = ilvoldpot

        do while (restpot.gt.lvpot(max(i1,1)).and.i1.ge.1)
          restpot = restpot-lvpot(i1) 
          i1 = i1-1
        enddo

! ---   check if some of the remaining leaves are older than span,
! ---   sum their weights

        dalvpot = 0.0d0
        if (lvagepot(max(i1,1)).gt.span.and.restpot.gt.0.and.           &
     &                          i1.ge.1) then
          dalvpot = lvpot(i1)-restpot
          restpot = 0.0d0
          i1 = i1-1
        endif

        do while (i1.ge.1.and.lvagepot(max(i1,1)).gt.span)
          dalvpot = dalvpot+lvpot(i1)
          i1 = i1-1
        enddo

        dalvpot = dalvpot/delt

! ---   death rate leaves and growth rate living leaves
        drlvpot   = dslvpot+dalvpot

! ---   leaf area not to exceed exponential growth curve
        slatpot = afgen (slatb,30,rid)
        if (laiexppot.lt.6.0d0) then
          dteff = max (0.0d0,tav-tbase)
          glaiexpot = laiexppot*rgrlai*dteff
! ---   source-limited increase in leaf area
          glasolpot = grlvpot*slatpot
          glapot = min (glaiexpot,glasolpot)
! ---   adjustment of specific leaf area of youngest leaf class
          if (grlvpot.gt.0.0d0) slatpot = glapot/grlvpot
        endif  

! ---   growth rate stems
        grstpot = fs*admipot
! ---   death of stems due to water stress is zero in case of potential growth
        drst1pot = 0.0d0
! ---   death of stems due to ageing
        drst2pot = afgen (rdrstb,30,rid)*wstpot
        drstpot = (drst1pot+drst2pot)/delt 
        gwstpot = grstpot-drstpot

! ----  integrals of the crop --------------------------------------------

!       set growing period after previous harvest        
        daygrowthpot = daygrowthpot + 1

!       Check trigger to start mowing event
        if (seqgrazmowpot(iseqgmpot) .eq. 2) then

          flharvestpot = .false.
            
          ! use dry matter threshold
          if (swharvest .eq. 1) then 
      
            ! use of fixed threshold
            if (swdmmow .eq. 1) then   
              if (tagppot .gt. dmharvest .or. (daynr .gt. daylastharvest  &
     &          .and. tagppot .gt. dmlastharvest)) then
                flharvestpot = .true.
              endif
            
            ! use of flexible threshold
            elseif (swdmmow .eq. 2) then 
              dmharvest = afgen(dmmowtb,20,rid)
              if (tagppot .gt. dmharvest .or.                           &
     &                 (daygrowthpot .gt. maxdaymow .and. iseqgmpot .gt. 1)) then
                flharvestpot = .true.
              endif
            endif
          
          ! use fixed dates
          elseif (swharvest .eq. 2) then  
            if(t1900 .gt. dateharvest(iharvest)) then
              flharvestpot = .true.
            endif
          endif

!         In case mowing is triggered: Growth is initialized again and the weight of the sward is stored
          if (flharvestpot) then
            iseqgmpot = iseqgmpot + 1
            slapot(1) = afgen (slatb,30,rid)
            fl = afgen (fltb,30,rid)
            fs = afgen (fstb,30,rid)
            wlvpot = mowrest / (1.d0 + (fs/fl))
            wstpot = fs/fl*wlvpot
            dwlvpot = 0.0d0
            dwstpot = 0.0d0
            lvagepot(1) = 0.0d0
            ilvoldpot = 1
            lasumpot = wlvpot * slapot(1)
            laiexppot = lasumpot
            lvpot(1) = wlvpot
            
            gwstpot = 0.0d0
            gwrtpot = 0.0d0
            drlvpot = 0.0d0
            drstpot = 0.0d0
            drrtpot = 0.0d0
            
            daygrowthpot = 0
            
!           losses due to treading
            fralossmow = 0.d0
            if (swlossmow.eq.1) then
              fralossmow = afgen(lossmowtab,200,h(nodmow))
            end if

!           harvest
            tagpspot = max(0.0d0,(tagppot-(wlvpot+dwlvpot+wstpot+dwstpot)))
            tagptpot = tagptpot + tagpspot * (1.d0 - FraLossMow)

            cropendpot  = rid
            flhrvendpot = .true.
            pmowdm      = tagpspot * (1.d0 - FraLossMow)
            plossdm     = tagpspot * FraLossMow
            
!           set regrowth delay
            idelaypot = int(afgen(DelayRegrowthTab,200,tagpspot))
            idregrpot = daycrop + idelaypot

          endif          
          
!       Check trigger to start grazing event          
        else if (seqgrazmowpot(iseqgmpot) .eq. 1 .or. seqgrazmowpot(iseqgmpot) .eq. 3) then

          flharvestpot = .false.
          
          if (.not. flgrazingpot) then
              
            ! use dry matter threshold
            if (swharvest .eq. 1) then 
            
              ! use of fixed threshold
              if (swdmgrz .eq. 1) then   
                if (tagppot .gt. dmgrazing) then
                  flharvestpot = .true.
                endif
              
              ! use of flexible threshold
              elseif (swdmgrz .eq. 2) then 
                dmgrazing = afgen(dmgrztb,20,rid)
                if (tagppot .gt. dmgrazing .or.                           &
     &                 (daygrowthpot .gt. maxdaygrz .and. iseqgmpot .gt. 1)) then
                  flharvestpot = .true.
                endif
              endif
              
            ! use fixed dates
            elseif (swharvest .eq. 2) then  
              if(t1900 .gt. dateharvest(iharvest)) then
                flharvestpot = .true.
              endif
            endif
          endif

!         In case grazing is triggered (or still occurs):
          if (flharvestpot .or. flgrazingpot) then
          
!           Amount of grazing kg/ha DM based on livestock density (Handboek Melkveehouderij 2013)
            uptgrazpot = lsda(iseqgmpot) *                                 &
     &                         afgen(uptgrazingtab,200,lsda(iseqgmpot))

!           Amount of shoots lost (kg/ha DM) due to droppings and treading during grazing  
            lossgrazpot = lsda(iseqgmpot) *                                &
     &                         afgen(lossgrazingtab,200,lsda(iseqgmpot))

!           Extra losses due to treading in case pressure head is insufficient
            fralossgrz = 0.d0
            if (swlossgrz.eq.1) then
              fralossgrz = afgen(lossgrztab,200,h(nodgrz))
            end if
            lossgrazpot = lossgrazpot + tagppot * fralossgrz

!           Initialise Count nr of days with grazing
            if(.not. flgrazingpot) then
              daygrowthpot   = 0
              idaysgrazpot   = 0
              flDewoolingpot = .false.
            endif
            
!           verify if uptake is possible: tagprest should remain after grazing
            if ((tagppot - uptgrazpot - lossgrazpot) .gt. tagprest) then
              
              flgrazingpot = .true.
              cuptgrazpot  = cuptgrazpot + uptgrazpot
          
!             distribute grazing over stems and leaves (living and dead parts)
              wstpot  = wstpot  - (uptgrazpot+lossgrazpot) * wstpot  / tagppot
              dwstpot = dwstpot - (uptgrazpot+lossgrazpot) * dwstpot / tagppot
              dwlvpot = dwlvpot - (uptgrazpot+lossgrazpot) * dwlvpot / tagppot
              grazlivinglvpot =   (uptgrazpot+lossgrazpot) * wlvpot  / tagppot
          
!             reduce leave weights
              i1 = ilvoldpot
              do while (grazlivinglvpot .gt. 0 .and. i1 .ge. 1)
                if (grazlivinglvpot .ge. lvpot(i1)) then
                  grazlivinglvpot = grazlivinglvpot - lvpot(i1)
                  lvpot(i1) = 0.0d0
                  i1 = i1 - 1
                else
                  lvpot(i1) = lvpot(i1) - grazlivinglvpot
                  grazlivinglvpot = 0.d0
                endif
              enddo
          
!             harvest during total grazing event
              cropendpot = rid
              pgrzdm     = pgrzdm + uptgrazpot
              plossdm    = tagppot * fralossgrz
              
!             Check number of days with grazing
              daysgrazpot  = int(afgen(daysgrazingtab,200,lsda(iseqgmpot)))
              idaysgrazpot = idaysgrazpot + 1
              if(idaysgrazpot .eq. daysgrazpot) then
                flgrazingpot = .false.
                flhrvendpot  = .true.
                if (seqgrazmowpot(iseqgmpot) .eq. 3) then
                  flDewoolingpot  = .true.
                endif
                daygrowthpot = 0
                iseqgmpot = iseqgmpot + 1
              endif

!           Also end grazing when not enough grass remains on the field
            elseif (flgrazingpot .or. swharvest .eq. 2) then
              flgrazingpot     = .false.
              flhrvendpot      = .true.
              flearlyhrvendpot = .true.
              if (seqgrazmowpot(iseqgmpot) .eq. 3 .and. tagppot .gt. dewrest) then
                flDewoolingpot   = .true.
                flearlyhrvendpot = .false.
              endif
              daygrowthpot = 0
              iseqgmpot = iseqgmpot + 1
            endif

!           Assumption: no delay in regrowth during and after grazing (without dewooling)
            idregrpot = daycrop

!           Dewooling after grazing event            
            if (flDewoolingpot) then

              slapot(1) = afgen (slatb,30,rid)
              fl = afgen (fltb,30,rid)
              fs = afgen (fstb,30,rid)
              wlvpot = dewrest / (1.d0 + (fs/fl))
              wstpot = fs/fl*wlvpot
              dwlvpot = 0.0d0
              dwstpot = 0.0d0
              lvagepot(1) = 0.0d0
              ilvoldpot = 1
              lasumpot = wlvpot * slapot(1)
              laiexppot = lasumpot
              lvpot(1) = wlvpot
              
              gwstpot = 0.0d0
              gwrtpot = 0.0d0
              drlvpot = 0.0d0
              drstpot = 0.0d0
              drrtpot = 0.0d0
              
!             Assumption: one day delay in regrowth after grazing
              idregrpot = daycrop + 1
              
            endif
            
          endif
          
        endif
        
        if (daycrop .ge. idregrpot) then

! ---     physiologic ageing of leaves per time step
          fysdel = max (0.0d0,(tav-tbase)/(35.0d0-tbase))

! ---     leaf death is imposed on array untill no more leaves have to die or all leaves are gone

          dslvtpot = dslvpot*delt
          i1 = ilvoldpot
           do while (dslvtpot.gt.0.and.i1.ge.1)
            if (dslvtpot.ge.lvpot(i1)) then
              dslvtpot = dslvtpot-lvpot(i1)
              lvpot(i1) = 0.0d0
              i1 = i1-1
            else
              lvpot(i1) = lvpot(i1)-dslvtpot
              dslvtpot = 0.0d0
            endif
          enddo

          if(i1.gt.0) then
            do while (lvagepot(max(i1,1)) .gt. span .and. i1 .ge. 1)
              lvpot(i1) = 0.0d0
              i1 = i1-1
            enddo
          endif
          ilvoldpot = i1

! ---     shifting of contents, integration of physiological age
          do i1 = ilvoldpot,1,-1
            lvpot(i1+1) = lvpot(i1)
            slapot(i1+1) = slapot(i1)
            lvagepot(i1+1) = lvagepot(i1)+fysdel*delt
          enddo
          ilvoldpot = ilvoldpot+1

! ---     new leaves in class 1
          lvpot(1) = grlvpot*delt
          slapot(1) = slatpot
          lvagepot(1) = 0.d0

! ---     calculation of new leaf area and weight
          lasumpot = 0.d0
          wlvpot = 0.d0
          do i1 = 1,ilvoldpot
            lasumpot = lasumpot+lvpot(i1)*slapot(i1)
            wlvpot = wlvpot+lvpot(i1)
          enddo

          laiexppot = laiexppot+glaiexpot*delt

        endif

! ---   dry weight of living plant organs
        wrtpot = wrtpot+gwrtpot*delt
        wstpot = wstpot+gwstpot*delt

! ---   dry weight of dead plant organs (roots,leaves & stems)
        dwrtpot = dwrtpot+drrtpot*delt
        dwlvpot = dwlvpot+drlvpot*delt
        dwstpot = dwstpot+drstpot*delt

! ---   dry weight of dead and living plant organs
        twlvpot = wlvpot+dwlvpot
        twstpot = wstpot+dwstpot
        tagppot = twlvpot+twstpot

! ---   leaf area index
        laipot = lasumpot+ssa*wstpot
!       prevent immediate lai reduction at emergence
!       KRO-BOO-20160403: suppressed because deviates from Wofost
!       laipot = max(laipot, laiem)

        ! root extension
        if (swrd.eq.1) then
          rdpot = afgen (rdtb,22,rid)
          rdpot = min(rdpot,rdm)
        elseif (swrd.eq.2) then
          rrpot = min (rdm-rdpot,rri)
          if (fr.le.0.0d0 .or. pgasspot.lt.1.0d0) rrpot = 0.0d0
          rdpot = rdpot + rrpot
        elseif (swrd.eq.3) then
          rdpot = afgen (rlwtb,22,wrtpot)
          rdpot = min(rdpot,rdm)
        endif          

      endif

      return

      case (3)

! === calculate actual rate and state variables ======================================

! --- check end of harvest
      if (flhrvendact) then
        if (flearlyhrvendact) then
          cropstartact  = rid - 1.d0
        else
          cropstartact  = rid
        endif
        mowdm        = 0.d0
        grzdm        = 0.d0
        lossdm       = 0.d0
      endif
      flhrvendact      = .false.
      flearlyhrvendact = .false.

! --- rates of change of the crop variables ---------------------------------------------
      
! --- skip in case of: tsum<tsum200, or 3 criteria (tsummttd), or regrowth
      if (flGrassGrowth .and. daycrop .ge. idregr) then

! ===   daily dry matter production ===

! ---   water stress reduction of pgass to gass
        if(dabs(ptra).lt.nihil) then
          reltr = 1.0d0
        else
          reltr = max(0.0d0,min(1.0d0,tra/ptra))
        endif
        gass = pgass * reltr

! ---   respiration and partitioning of carbohydrates between growth and
! ---   maintenance respiration
        rmres = (rmr*wrt+rml*wlv+rms*wst)*afgen(rfsetb,30,rid)
        teff = q10**((tav-25.0d0)/10.0d0)
        mres = min (gass,rmres*teff)
        asrc = gass-mres

! ---   partitioning factors (relevant for restart)
        fr = afgen(frtb,30,rid)
        fl = afgen(fltb,30,rid)
        fs = afgen(fstb,30,rid)

! ---   dry matter increase
        cvf = 1.0d0/((fl/cvl+fs/cvs)*(1.0d0-fr)+fr/cvr)
        dmi = cvf*asrc
! ---   check on carbon balance
        ccheck = (gass-mres-(fr+(fl+fs)*(1.0d0-fr))*dmi/cvf)            &
     &         /max(0.0001d0,gass)      
        if (dabs(ccheck).gt.0.0001d0) then
          Messag ='The carbon balance is not correct'
          call fatalerr ('grass_act',messag)
        endif

! ===   growth rate by plant organ ===

! ---   growth rate roots and aerial parts
        ! in case of SWRD = 3: after reaching maximum live weight of wrtmax, the
        ! growth of the roots is balanced by the death of root tissue
        grrt = fr*dmi
        if (swrd.eq.3 .and. flWrtNonox) grrt = 0.d0
        if (swrd.eq.3 .and. wrt.gt.wrtmax) then
          drrt = grrt
          drrt = max(drrt,wrt*afgen (rdrrtb,30,rid))
        else  
          drrt = wrt*afgen (rdrrtb,30,rid)
        endif  
        gwrt = grrt-drrt

! ---   growth rate leaves

! ---   weight of new leaves
        admi = (1.0d0-fr)*dmi        
        grlv = fl*admi

! ---   death of leaves due to water stress or high lai
        dslv1 = wlv*(1.0d0-reltr)*perdl
        laicr = 3.2d0/kdif
        dslv2 = wlv*max(0.0d0,min(0.03d0,0.03d0*(lai-laicr)/laicr))
        dslv = max (dslv1,dslv2) 

! ---   death of leaves due to exceeding life span;
! ---   leaf death is imposed on array until no more leaves have
! ---   to die or all leaves are gone

        rest = dslv*delt
        i1 = ilvold

        do while (rest.gt.lv(max(i1,1)).and.i1.ge.1)
          rest = rest-lv(i1) 
          i1 = i1-1
        enddo

! ---   check if some of the remaining leaves are older than span,
! ---   sum their weights

        dalv = 0.0d0
        if (lvage(max(i1,1)).gt.span.and.rest.gt.0.and.i1.ge.1) then
          dalv = lv(i1)-rest
          rest = 0.0d0
          i1 = i1-1
        endif

        do while (i1.ge.1.and.lvage(max(i1,1)).gt.span)
          dalv = dalv+lv(i1)
          i1 = i1-1
        enddo

        dalv = dalv/delt

! ---   death rate leaves and growth rate living leaves
        drlv   = dslv+dalv

! ---   physiologic ageing of leaves per time step
        slat = afgen (slatb,30,rid)

! ---   leaf area not to exceed exponential growth curve
        if (laiexp.lt.6.0d0) then
          dteff = max (0.0d0,tav-tbase)
          glaiex = laiexp*rgrlai*dteff
! ---     source-limited increase in leaf area
          glasol = grlv*slat
          gla = min (glaiex,glasol)
! ---     adjustment of specific leaf area of youngest leaf class
          if (grlv.gt.0.0d0) slat = gla/grlv
        endif  

! ---   growth rate stems
        grst = fs*admi
! ---   death of stems due to water stress
        drst1 = wst*(1.0d0-reltr)*perdl
! ---   death of stems due to ageing
        drst2 = afgen (rdrstb,30,rid)*wst
        drst = (drst1+drst2)/delt 
        gwst = grst-drst

! ----  integrals of the crop --------------------------------------------

!       set growing period after previous harvest        
        daygrowth = daygrowth + 1

!       Check trigger to start mowing event
        if (seqgrazmow(iseqgm) .eq. 2) then
          
          flharvest = .false.   
            
          ! use dry matter threshold
          if (swharvest .eq. 1) then 
      
            ! use of fixed threshold
            if (swdmmow .eq. 1) then   
              if (tagp .gt. dmharvest .or. (daynr .gt. daylastharvest  &
     &          .and. tagp .gt. dmlastharvest)) then
                flharvest = .true.
              endif
            
            ! use of flexible threshold
            elseif (swdmmow .eq. 2) then 
              dmharvest = afgen(dmmowtb,20,rid)
              if (tagp .gt. dmharvest .or.                           &
     &                 (daygrowth .gt. maxdaymow .and. iseqgm .gt. 1)) then
                flharvest = .true.
              endif
            endif
          
          ! use fixed dates
          elseif (swharvest .eq. 2) then  
            if(t1900 .gt. dateharvest(iharvest)) then
              iharvest = iharvest + 1
              flharvest = .true.
            endif
          endif
          
!       In case mowing is triggered: Growth is initialized again and the weight of the sward is stored
        if (flharvest) then
          iseqgm = iseqgm + 1
          sla(1) = afgen (slatb,30,rid)
          fl = afgen (fltb,30,rid)
          fs = afgen (fstb,30,rid)
          wlv = mowrest / (1.d0 + (fs/fl))
          wst = fs/fl*wlv
          dwlv = 0.0d0
          dwst = 0.0d0
          lvage(1) = 0.0d0
          ilvold = 1
          lasum = wlv * sla(1)
          laiexp = lasum
          lv(1) = wlv

          gwst = 0.0d0
          gwrt = 0.0d0
          drlv = 0.0d0
          drst = 0.0d0
          drrt = 0.0d0

          daygrowth = 0

!         losses due to treading
          FraLossMow = 0.d0
          if (swlossmow.eq.1) then
            FraLossMow = afgen(lossmowtab,200,h(nodmow))
          end if
          
!         harvest
          tagps = max (0.0d0,(tagp-(wlv+dwlv+wst+dwst)))
          tagpt = tagpt + tagps * (1.d0 - fralossmow)

          cropendact  = rid
          flhrvendact = .true.
          mowdm   = tagps * (1.d0 - FraLossMow)
          lossdm  = tagps * FraLossMow
          
! ---     set regrowth delay
          idelay = int(afgen(DelayRegrowthTab,200,tagps))
          idregr = daycrop + idelay

        endif
          
!       Check trigger to start grazing event          
        else if (seqgrazmow(iseqgm) .eq. 1 .or. seqgrazmow(iseqgm) .eq. 3) then

          flharvest = .false.
            
          if (.not. flgrazing) then
            
            ! use dry matter threshold
            if (swharvest .eq. 1) then 
            
              ! use of fixed threshold
              if (swdmgrz .eq. 1) then   
                if (tagp .gt. dmgrazing) then
                  flharvest = .true.
                endif
              
              ! use of flexible threshold
              elseif (swdmgrz .eq. 2) then 
                dmgrazing = afgen(dmgrztb,20,rid)
                if (tagp .gt. dmgrazing .or.                           &
     &            (daygrowth .gt. maxdaygrz .and. iseqgm .gt. 1)) then
                  flharvest = .true.
                endif
              endif
            
            ! use fixed dates
            elseif (swharvest .eq. 2) then  
              if(t1900 .gt. dateharvest(iharvest)) then
                iharvest = iharvest + 1
                flharvest = .true.
              endif
            endif
          endif

!         In case grazing is triggered (or still occurs):
          if (flharvest .or. flgrazing) then
          
!           Amount of grazing kg/ha DM based on livestock density (Handboek Melkveehouderij 2013)
            uptgraz = lsda(iseqgm)*afgen(uptgrazingtab,200,lsda(iseqgm))            

!           Amount of shoots lost (kg/ha DM) due to droppings and treading during grazing  
            lossgraz = lsda(iseqgm) *                                &
     &                         afgen(lossgrazingtab,200,lsda(iseqgm))

!           Extra losses due to treading in case pressure head is insufficient
            fralossgrz = 0.d0
            if (swlossgrz.eq.1) then
              fralossgrz = afgen(lossgrztab,200,h(nodgrz))
            end if
            lossgraz = lossgraz + tagp * fralossgrz

!           Initialise Count nr of days with grazing
            if(.not. flgrazing) then
              daygrowth   = 0
              idaysgraz   = 0
              flDewooling = .false.
            endif
            
!           verify if uptake is possible: tagprest should remain after grazing
            if ((tagp - uptgraz - lossgraz) .gt. tagprest) then
              
              flgrazing = .true.
              cuptgraz  = cuptgraz + uptgraz

!             distribute grazing over stems and leaves (living and dead parts)
              wst  = wst  -  (uptgraz+lossgraz) * wst  / tagp
              dwst = dwst -  (uptgraz+lossgraz) * dwst / tagp
              dwlv = dwlv -  (uptgraz+lossgraz) * dwlv / tagp
              grazlivinglv = (uptgraz+lossgraz) * wlv  / tagp
          
!             reduce leave weights
              i1 = ilvold
              do while (grazlivinglv .gt. 0 .and. i1 .ge. 1)
                if (grazlivinglv .ge. lv(i1)) then
                  grazlivinglv = grazlivinglv - lv(i1)
                  lv(i1) = 0.0d0
                  i1 = i1 - 1
                else
                  lv(i1) = lv(i1) - grazlivinglv
                  grazlivinglv = 0.d0
                endif
              enddo
          
!             harvest during total grazing event
              cropendact = rid
              grzdm      = grzdm + uptgraz
              lossdm     = tagp * fralossgrz
              
!             Check number of days with grazing
              daysgraz  = int(afgen(daysgrazingtab,200,lsda(iseqgm)))
              idaysgraz = idaysgraz + 1
              if(idaysgraz .eq. daysgraz) then
                flgrazing   = .false.
                flhrvendact = .true.
                if (seqgrazmow(iseqgm) .eq. 3) then
                  flDewooling  = .true.
                endif
                daygrowth = 0
                iseqgm = iseqgm + 1
              endif

!           Also end grazing when not enough grass remains on the field
            elseif (flgrazing .or. swharvest .eq. 2) then
              flgrazing        = .false.
              flhrvendact      = .true.
              flearlyhrvendact = .true.
              if (seqgrazmow(iseqgm) .eq. 3 .and. tagp .gt. dewrest) then
                flDewooling      = .true.
                flearlyhrvendact = .false.
              endif
              daygrowth = 0
              iseqgm = iseqgm + 1
            endif

!           Assumption: no delay in regrowth during and after grazing (without dewooling)
            idregr = daycrop

!           Dewooling after grazing event            
            if (flDewooling) then

              sla(1) = afgen (slatb,30,rid)
              fl = afgen (fltb,30,rid)
              fs = afgen (fstb,30,rid)
              wlv = dewrest / (1.d0 + (fs/fl))
              wst = fs/fl*wlv
              dwlv = 0.0d0
              dwst = 0.0d0
              lvage(1) = 0.0d0
              ilvold = 1
              lasum = wlv * sla(1)
              laiexp = lasum
              lv(1) = wlv
    
              gwst = 0.0d0
              gwrt = 0.0d0
              drlv = 0.0d0
              drst = 0.0d0
              drrt = 0.0d0
    
!             Assumption: one day delay in regrowth after grazing
              idregr = daycrop + 1

            endif
            
          endif
          
        endif

        if (daycrop .ge. idregr) then

! ---     physiologic ageing of leaves per time step
          fysdel = max (0.0d0,(tav-tbase)/(35.0d0-tbase))

! ---     leaf death is imposed on array untill no more leaves have to die or all leaves are gone

          dslvt = dslv*delt
          i1 = ilvold
          do while (dslvt.gt.0.and.i1.ge.1)
            if (dslvt.ge.lv(i1)) then
              dslvt = dslvt-lv(i1)
              lv(i1) = 0.0d0
              i1 = i1-1
            else
              lv(i1) = lv(i1)-dslvt
              dslvt = 0.0d0
            endif
          enddo

          if(i1.gt.0) then
            do while (lvage(max(i1,1)).gt.span.and.i1.ge.1)
              lv(i1) = 0.0d0
              i1 = i1-1
            enddo
          endif
          ilvold = i1

! ---     shifting of contents, integration of physiological age
          do i1 = ilvold,1,-1
            lv(i1+1) = lv(i1)
            sla(i1+1) = sla(i1)
            lvage(i1+1) = lvage(i1)+fysdel*delt
          enddo
          ilvold = ilvold+1

! ---     new leaves in class 1
          lv(1) = grlv*delt
          sla(1) = slat
          lvage(1) = 0.d0 

! ---     calculation of new leaf area and weight
          lasum = 0.d0
          wlv = 0.d0
          do i1 = 1,ilvold
            lasum = lasum+lv(i1)*sla(i1)
            wlv = wlv+lv(i1)
          enddo

          laiexp = laiexp+glaiex*delt

        endif

! ---   dry weight of living plant organs
        wrt = wrt+gwrt*delt
        wst = wst+gwst*delt

! ---   dry weight of dead plant organs (roots,leaves & stems)
        dwrt = dwrt+drrt*delt
        dwlv = dwlv+drlv*delt
        dwst = dwst+drst*delt

! ---   dry weight of dead and living plant organs
        twlv = wlv+dwlv
        twst = wst+dwst
        tagp = twlv+twst

! ---   leaf area index
        lai = lasum+ssa*wst
        laimax = max (lai,laimax)

! ---   update normalized cumulative root density based on root extraction or stress (cumdens)
        if (swrdc .eq. 1) call update_rootdistribution()
        
        ! root extension
        if (swrd.eq.1) then
          rd = afgen (rdtb,22,rid)
          rd = min(rd,rdm)
        elseif (swrd.eq.2) then
          rr = min (rdm-rd,rri)
          if (fr.le.0.0d0 .or. pgass.lt.1.0d0 .or. flWrtNonox) rr = 0.0d0
          if (swdmi2rd.eq.1 .and. pgass.ge.1.0d0)              rr = rr * gass/pgass
          rd = rd + rr
        elseif (swrd.eq.3) then
          rd = afgen (rlwtb,22,wrt)
          rd = min(rd,rdm)
        endif            

! ---   set crop height and cropfactor
        if (swcf.ne.3) then
          cf = afgen (cftb,(2*magrs),rid)
          ch = afgen (chtb,(2*magrs),rid)
        else
          cf = afgen (cftb,(2*magrs),lai)
          cfeic = afgen (cfeictb,(2*magrs),lai)
          ch = afgen(chtb,(2*magrs),lai)
        endif

! ---   update canopy storage capacity
        if (swinter.eq.3) then
          siccapact = siccaplai*lai
        endif

      endif

      return

      case default
         call fatalerr ('Grass', 'Illegal value for TASK')
      end select

      return
      end

! ----------------------------------------------------------------------
      subroutine update_rootdistribution()

! ----------------------------------------------------------------------
!     Date: May 2021
!     Purpose: dynamic root distribution
!              The normalized cumulative root density is modified by
!              growth of root biomass based on relative root water 
!              extraction or transpiration reduction (uncompensated).
! ----------------------------------------------------------------------
      
      ! import global variables
      use variables, only: date, noddrz, zbotcp, ztopcp, cumdens,             &
                     qpotrot_day, qredtot_day,                                &
                     wrt, gwrt, wrtmin
      ! local      
      implicit none
 
      integer   node, i
      real(8)   top,bot
      real(8)   rd_noddrz
      real(8)   rel_qrot_day, rel_qred_day, sum
      real(8)   wrttot, wrtdis(202), qrotdis(202), qreddis(202)
      logical   found
      
! --- update normalized cumulative root density based on root extraction or stress (cumdens)
!      if (swrdc .eq. 1) then
          
        ! root extraction of each compartment since start of the day
        rel_qrot_day = 0.d0
        rel_qred_day = 0.d0
        do node = 1,noddrz
          rel_qrot_day = rel_qrot_day + 1 - qredtot_day(node) / qpotrot_day(node)
          rel_qred_day = rel_qred_day + qredtot_day(node)
        enddo
        
        if ((gwrt .gt. 0.d0 .and. rel_qrot_day .gt. 0.d0) .or. (gwrt .lt. 0.d0 .and. rel_qred_day .gt. 0.d0)) then
        
          ! distribution roots and root extraction at relative depth
          ! root extraction and root weight based on previous day
          rd_noddrz = abs(zbotcp(noddrz))
          node = 1
          do i = 4,202,2
            
            ! root distribution of previous day
            wrtdis(i) = (cumdens(i) - cumdens(i-2)) * (wrt - gwrt)
          
            ! determine optimal extraction and maximum reduction at relative depth
            found = .false.
            qrotdis(i) = 0.d0
            qreddis(i) = 0.d0
            top = - cumdens(i-3) * rd_noddrz
            bot = - cumdens(i-1) * rd_noddrz
            do while (.not. found)
              if (bot .ge. zbotcp(node)) then
                qrotdis(i) = qrotdis(i) + (1 - qredtot_day(node) / qpotrot_day(node)) / (ztopcp(node) - zbotcp(node)) * (top - bot)
                qreddis(i) = qreddis(i) + qredtot_day(node) / (ztopcp(node) - zbotcp(node)) * (top - bot)
                found = .true.
              else
                qrotdis(i) = qrotdis(i) + (1 - qredtot_day(node) / qpotrot_day(node)) / (ztopcp(node) - zbotcp(node)) * (top - zbotcp(node))
                qreddis(i) = qreddis(i) + qredtot_day(node) / (ztopcp(node) - zbotcp(node)) * (top - zbotcp(node))
                top = zbotcp(node)
                node = node + 1
              end if
            end do
          end do
          
          ! update relative root weight
          wrttot = 0.d0
          if (gwrt .gt. 0.d0) then
            do i = 4,202,2
              wrtdis(i) = max(wrtmin, wrtdis(i) + (qrotdis(i) / rel_qrot_day) * gwrt)
              wrttot = wrttot + wrtdis(i)
            end do
          elseif (gwrt .lt. 0.d0) then
            do i = 4,202,2
              wrtdis(i) = max(wrtmin, wrtdis(i) + (qreddis(i) / rel_qred_day) * gwrt)
              wrttot = wrttot + wrtdis(i)
            end do
          end if
          
          ! update normalized cumulative root density distribution
          sum = 0.d0
          do i = 4,202,2
            sum = sum + wrtdis(i)
            cumdens(i) = sum / wrttot
          end do

        end if
  
        ! TEMPORARY OUTPUT  DELETE
        do i = 2,202,2
          write(777,*) trim(date), ",", cumdens(i-1), ",", cumdens(i)
        end do

        do node = 1,noddrz
          write(888,'(a11,",",i4,3(",",f15.5))') trim(date), node, qpotrot_day(node), qredtot_day(node)
        end do
        ! TEMPORARY OUTPUT  DELETE
        
!      end if
    
      return
      end
    
! ----------------------------------------------------------------------
      subroutine totass (dayl,amax,eff,lai,kdif,avrad,difpp,            &
     &                   dsinbe,sinld,cosld,dtga)

!*  Purpose: This routine calculates the daily total gross CO2
!*           assimilation by performing a Gaussian integration over
!*           time. At three different times of the day, irradiance is
!*           computed and used to calculate the instantaneous canopy
!*           assimilation, whereafter integration takes place. More
!*           information on this routine is given by Spitters et al.
!*           (1988).

!*  FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)
!*  name   type meaning                                    units  class
!*  ----   ---- -------                                    -----  -----
!*  DAYL    R8  Astronomical daylength (base = 0 degrees)     h      O
!*  AMAX    R8  Assimilation rate at light saturation      kg CO2/   I
!*                                                        ha leaf/h   
!*  EFF     R8  Initial light use efficiency              kg CO2/J/  I
!*                                                        ha/h m2 s   
!*  LAI     R8  Leaf area index                             ha/ha    I
!*  KDIF    R8  Extinction coefficient for diffuse light             I
!*  AVRAD   R8  Daily shortwave radiation                  J m-2 d-1 I
!*  DIFPP   R8  Diffuse irradiation perpendicular to direction of
!*              light                                      J m-2 s-1 I
!*  DSINBE  R8  Daily total of effective solar height         s      I
!*  SINLD   R8  Seasonal offset of sine of solar height       -      I
!*  COSLD   R8  Amplitude of sine of solar height             -      I
!*  DTGA    R8  Daily total gross assimilation           kg CO2/ha/d O

!*  FATAL ERROR CHECKS: none
!*  SUBROUTINES and FUNCTIONS called : ASSIM
!*  FILE usage : none

!*  Authors: Daniel van Kraalingen 
!*  Date   : April 1991

!*  Modification: Implementated in Swap3.2.41: 
!*                several small adjustments (R4->R8, small caps)
!*  Author      : Joop Kroes
!*  Date        : May 2014

      implicit none

!*     formal parameters
      real(8) dayl,amax,eff,lai,kdif,avrad,difpp,dsinbe,sinld,cosld,dtga

!*     local parameters
      integer i1
      real(8) hour,pi,sinb,par,pardif,pardir,fgros
      real(8) xgauss(3),wgauss(3)

      parameter (pi=3.1415926d0)
      save

!**
!*     gauss points and weights are stored in an array
      data xgauss /0.1127017d0, 0.5000000d0, 0.8872983d0/
      data wgauss /0.2777778d0, 0.4444444d0, 0.2777778d0/

!*     calculation of assimilation is done only when it will not be zero
!*     (AMAX >0, LAI >0)
      dtga  = 0.0d0
      if (amax.gt.0.0d0.and.lai.gt.0.0d0) then
         do 10 i1=1,3
            hour   = 12.0d0+0.5d0*dayl*xgauss(i1)
            sinb   = max(0.0d0,                                         &
     &                  sinld+cosld*dcos(2.0d0*pi*(hour+12.0d0)/24.0d0))
            par    = 0.5d0*avrad*sinb*(1.0d0+0.4d0*sinb)/dsinbe
            pardif = min(par,sinb*difpp)
            pardir = par-pardif
            call assim (amax,eff,lai,kdif,sinb,pardir,pardif,fgros)
            dtga = dtga+fgros*wgauss(i1)
10       continue
         dtga = dtga*dayl
      end if

      return
      end

      subroutine assim (amax,eff,lai,kdif,sinb,pardir,pardif,fgros)

!*     Chapter 13 in documentation WOFOST Version 4.1 (1988)

!*     This routine calculates the gross CO2 assimilation rate of
!*     the whole crop, FGROS, by performing a Gaussian integration
!*     over depth in the crop canopy. At three different depths in
!*     the canopy, i.e. for different values of LAI, the
!*     assimilation rate is computed for given fluxes of photosynthe-
!*     tically active radiation, whereafter integration over depth
!*     takes place. More information on this routine is given by
!*     Spitters et al. (1988). The input variables SINB, PARDIR
!*     and PARDIF are calculated in routine TOTASS.

!*     Subroutines and functions called: none.
!*     Called by routine TOTASS.

!*     Author: D.W.G. van Kraalingen, 1986

!*  Modification: Implementated in Swap3.2.41: 
!*                several small adjustments (R4->R8, small caps)
!*  Author      : Joop Kroes
!*  Date        : May 2014

!*  FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)
!*  name   type meaning                                    units  class
!*  ----   ---- -------                                    -----  -----
!*  AMAX    R8  Maximum CO2 assimilation rate              kg/ha/hr  I
!*  EFF     R8  Light use efficiency of a leaf         kg CO2 / J adsorbed  I
!*  LAI     R8  Leaf area index                               -      I
!*  KDIF    R8  Extinction coefficient for diffuse visible light -   I
!*  SINB    R8  ...........nog invullen ..............               -      I
!*  PARDIR  R8  ...........nog invullen ..............                 -      I
!*  DIFPP   R8  Diffuse irradiation perpendicular to direction of     
!*              light                                      J m-2 s-1 I
!*  PARDIF  R8  ...........nog invullen ..............             -      I
!*  FGROS   R8  ...........nog invullen ..............           s      O
!**    
!*13.1 declarations
      implicit none

!*     formal parameters
      real(8) amax,eff,lai,kdif,sinb,pardir,pardif,fgros

!*     local parameters
      integer i
      real(8) scv,refh,refs,kdirbl,kdirt,laic,visdf,vist,visd,visshd
      real(8) fgrsh,vispp,fgrsun,fslla,fgl
      real(8) xgauss(3),wgauss(3)

      save

!*     initialize GAUSS array and scattering coefficient
      data xgauss /0.1127017d0, 0.5000000d0, 0.8872983d0/
      data wgauss /0.2777778d0, 0.4444444d0, 0.2777778d0/
      data scv /0.2d0/

!*13.2 extinction coefficients KDIF,KDIRBL,KDIRT
      refh   = (1.0d0-dsqrt(1.0d0-scv))/(1.0d0+dsqrt(1.0d0-scv))
      refs   = refh*2.0d0/(1.0d0+1.6d0*sinb)
      kdirbl = (0.5d0/sinb)*kdif/(0.8d0*dsqrt(1.0d0-scv))
      kdirt  = kdirbl*dsqrt(1.0d0-scv)

!*13.3 three-point Gaussian integration over LAI
      fgros  = 0.0d0
      do 10 i=1,3
         laic   = lai*xgauss(i)
!*        absorbed diffuse radiation (VISDF),light from direct
!*        origine (VIST) and direct light(VISD)
         visdf  = (1.0d0-refs)*pardif*kdif  *exp (-kdif  *laic)
         vist   = (1.0d0-refs)*pardir*kdirt *exp (-kdirt *laic)
         visd   = (1.0d0-scv) *pardir*kdirbl*exp (-kdirbl*laic)
!*        absorbed flux in W/m2 for shaded leaves and assimilation
         visshd = visdf+vist-visd
         fgrsh  = amax*(1.0d0-exp(-visshd*eff/max(2.0d0,amax)))
!*        direct light absorbed by leaves perpendicular on direct
!*        beam and assimilation of sunlit leaf area
         vispp  = (1.0d0-scv)*pardir/sinb
         if (vispp.le.0.0d0) then
            fgrsun = fgrsh
         else
            fgrsun = amax*(1.0d0-(amax-fgrsh)                           &
     &          *(1.0d0-exp (-vispp*eff/max(2.0d0,amax)))/ (eff*vispp))
         end if
!*        fraction of sunlit leaf area (FSLLA) and local
!*        assimilation rate (FGL)
         fslla  = exp (-kdirbl*laic)
         fgl    = fslla*fgrsun+(1.0d0-fslla)*fgrsh
!*        integration
         fgros  = fgros+fgl*wgauss(i)
10    continue

      fgros  = fgros*lai
      return
      end

      subroutine astro (iday,lat,avrad,                                 &
     &                  dayl,daylp,sinld,cosld,difpp,atmtr,dsinbe)

!*  Purpose: This subroutine calculates astronomic daylength,
!*           diurnal radiation characteristics such as the atmospheric
!*           transmission, diffuse radiation etc.. This routine has
!*           been modified so that it uses arrays to hold some input
!*           output variables for faster processing 

!*  FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)
!*  name   type meaning                                    units  class
!*  ----   ---- -------                                    -----  -----
!*  IDAY    I4  Day number (Jan 1st = 1)                      -      I
!*  LAT     R8  Latitude of the site                       degrees   I
!*  AVRAD   R8  Daily shortwave radiation                  J m-2 d-1 I
!*  DAYL    R8  Astronomical daylength (base = 0 degrees)     h      O
!*  DAYLP   R8  Astronomical daylength (base =-4 degrees)     h      O
!*  SINLD   R8  Seasonal offset of sine of solar height       -      O
!*  COSLD   R8  Amplitude of sine of solar height             -      O
!*  DIFPP   R8  Diffuse irradiation perpendicular to direction of     
!*              light                                      J m-2 s-1 O
!*  ATMTR   R8  Daily atmospheric transmission                -      0
!*  DSINBE  R8  Daily total of effective solar height         s      O

!*  FATAL ERROR CHECKS: none
!*  SUBROUTINES and FUNCTIONS called : none
!*  FILE usage : none

!*  Authors: Daniel van Kraalingen
!*  Date   : April 1991

!*  Modification: Include checks for 0<=daylength<=24 hour
!*                Remove caching of results
!*  Author      : Allard de Wit
!*  Date        : January 2011

!*  Modification: Implementated in Swap3.2.41: 
!*                several small adjustments (R4->R8, small caps)
!*  Author      : Joop Kroes
!*  Date        : May 2014

      implicit none
!*     formal parameters
      integer iday
      real(8) lat,avrad,dayl,daylp,sinld,cosld,difpp,atmtr,dsinbe

!*     local parameters
      real(8) pi,angle,rad
      real(8) dec,sc,aob,aob_corr,angot,dsinb,frdif

      parameter (pi=3.1415926d0, angle=-4.0d0, rad=0.0174533d0)

!*     Error check on latitude
      if (dabs(lat).gt.90.d0) call fatalerr                             &
     &   ('astro','lat > 90 or lat < -90')

!*     Declination and solar constant for this day
      dec = -asin(dsin(23.45d0*rad)*dcos(2.d0*pi*dble(iday+10)/365.0d0))
      sc  = 1370.d0*(1.d0+0.033d0*dcos(2.d0*pi*dble(iday)/365.d0))

!*     calculation of daylength from intermediate variables
!*     SINLD, COSLD and AOB
      sinld = dsin(rad*lat)*dsin(dec)
      cosld = dcos(rad*lat)*dcos(dec)
      aob = sinld/cosld

!*     For very high latitudes and days in summer and winter a limit is  
!*     inserted to avoid math errors when daylength reaches 24 hours in 
!*     summer or 0 hours in winter.

!*     Calculate solution for base=0 degrees
      if (dabs(aob).le.1.0d0) then
         dayl  = 12.0d0*(1.d0+2.d0*asin(aob)/pi)
!*        integrals of sine of solar height
         dsinb  = 3600.d0*                                              &
     &            (dayl*sinld+24.d0*cosld*dsqrt(1.d0-aob**2)/pi)
         dsinbe = 3600.d0*                                              &
     &            (dayl*(sinld+0.4d0*(sinld**2+cosld**2*0.5d0))+  &
     &     12.d0*cosld*(2.d0+3.d0*0.4d0*sinld)*dsqrt(1.d0-aob**2)/pi)
      else
         if (aob.gt.1.0d0)  dayl = 24.0d0
         if (aob.lt.-1.0d0) dayl =  0.0d0
!*        integrals of sine of solar height      
         dsinb  = 3600.d0*(dayl*sinld)
         dsinbe = 3600.d0*                                              &
     &            (dayl*(sinld+0.4d0*(sinld**2+cosld**2*0.5d0)))
      endif

!*     Calculate solution for base=-4 (ANGLE) degrees
      aob_corr = (-dsin(angle*rad)+sinld)/cosld
      if (dabs(aob_corr).le.1.0d0) then 
         daylp = 12.0d0*(1.d0+2.d0*asin(aob_corr)/pi)
      else
         if (aob_corr.gt.1.0d0)  daylp = 24.0d0
         if (aob_corr.lt.-1.0d0) daylp =  0.0d0
      endif

!*     extraterrestrial radiation and atmospheric transmission
      angot  = sc*dsinb
!*     Check for DAYL=0 as in that case the angot radiation is 0 as well
      if (dayl.gt.0.0d0) then
          atmtr = avrad/angot
      else
          atmtr = 0.0d0
      endif

!*     estimate fraction diffuse irradiation
      if (atmtr.gt.0.75d0) frdif = 0.23d0
      if (atmtr.le.0.75d0.and.atmtr.gt.0.35d0)                          &
     &  frdif = 1.33d0-1.46d0*atmtr
      if (atmtr.le.0.35d0.and.atmtr.gt.0.07d0)                          &
     &  frdif = 1.d0-2.3d0*(atmtr-0.07d0)**2
      if (atmtr.le.0.07d0) frdif = 1.d0

      difpp = frdif*atmtr*0.5d0*sc

      RETURN
      END

! ----------------------------------------------------------------------
      subroutine outbalcropN(task,pathwork,outfil,project,date,daycrop, &
     &             t,dvs,tsum,NUPTT,NFIXTT,ANLVI,ANSTI,ANRTI,ANSOI,ANLV,&
     &               ANST,ANRT,ANSO,NLOSSL,NLOSSR,NLOSSS,NBALAN,NNI) 
! ----------------------------------------------------------------------
!     Date               : March 2015  
!     Purpose            : open and write crop output N balance files 
! ----------------------------------------------------------------------
      implicit none

! --- global variables ------------------
      character(len=11) date
      character(len=*) outfil,pathwork,project
      integer task,daycrop
      real(8) t,dvs,tsum   !,laipot,lai,cf,rdpot,rd,ch,crt0,crt1
!      real(8) cwdmpot,cwdm,wsopot,wso,wstpot,wst,wlvpot,wlv,wrtpot,wrt
      real(8) NUPTT,NFIXTT,ANLVI,ANSTI,ANRTI,ANSOI,ANLV
      real(8) ANST,ANRT,ANSO,NLOSSL,NLOSSR,NLOSSS,NBALAN,NNI
! --- local variables ------------------
      character(len=1) comma
!      character(len=200) messag
      integer   getun,nba
      character(len=160) filnam,filtext

      save    nba

      comma = ',' 
    
      select case (task)
      case (1)

! === open output file and write headers =====================

      filnam = trim(pathwork)//trim(outfil)//'.nba'
      nba = getun (20,90)
      call fopens(nba,filnam,'new','del')
      filtext = 'output of N-balance of detailed crop growth model'
      call writehead (nba,1,filnam,filtext,project)

      write (nba,100)
 100    format ('*',/,                                                  &
     & '*             day     day      -    grC   kg/ha   kg/ha   kg/', &
     & 'ha   kg/ha   kg/ha   kg/ha   kg/ha   kg/ha   kg/ha   kg/ha   ', &
     & 'kg/ha   kg/ha   kg/ha  kg/ha',/,                                &
     & '      Date, Daynr, Daycrp,   DVS,  TSUM,  NUPTT, NFIXTT,  ANLVI'&
     & ,',  ANSTI,  ANRTI,  ANSOI,   ANLV,   ANST,   ANRT',             &
     & ',   ANSO, NLOSSL, NLOSSR, NLOSSS, NBALAN, NNI')
      return

      case (2)
! --- write dynamic data ----------------------------------------------------
      write (nba,200) date,comma,nint(t),comma,daycrop,comma,dvs,comma, &
     & tsum,comma,NUPTT,comma,NFIXTT,comma,ANLVI,comma,ANSTI,comma,     &
     & ANRTI,comma,ANSOI,comma,ANLV,comma,ANST,comma,ANRT,comma,        &
     & ANSO,comma,NLOSSL,comma,NLOSSR,comma,NLOSSS,comma,NBALAN,        &
     & comma,NNI
 200  format (a11,a1,i5,a1,i7,a1,f6.2,a1,f6.0, 15(a1,f7.2) )

      return

      case (3)
! --- close crop output file ------------------------------------------------

      close (nba)

      case default
         call fatalerr ('OutbalCropN', 'Illegal value for TASK')
      end select

      return
      end 

! ----------------------------------------------------------------------
      subroutine outbalcropOM1(task,pathwork,outfil,project,date,       &
     &      daycrop,t,dvs,tsum,gass,mres,fr,fl,fs,fo,dmi,cvf,ccheck)
! ----------------------------------------------------------------------
!     Date               : March 2015  
!     Purpose            : open and write crop output OM balance files 
! ----------------------------------------------------------------------
      implicit none
      include 'params.fi'

! --- global variables ------------------
      character(len=11) date
      character(len=*) outfil,pathwork,project
      integer task,daycrop
      real(8) t,dvs,tsum,gass,mres,fr,fl,fs,fo,dmi,cvf,ccheck
! --- local variables ------------------
      character(len=1) comma
!      character(len=200) messag
      integer   getun,om1
      character(len=160) filnam,filtext
      real(8)   OMroot,OMleaves,OMstems,OMstorage  !,Cccheck

      save    om1

      comma = ',' 

      select case (task)
      case (1)

! === open output file and write headers =====================

      filnam = trim(pathwork)//trim(outfil)//'.om1'
      om1 = getun (20,90)
      call fopens(om1,filnam,'new','del')
      filtext = 'output of OM1-balance (kg/ha DM increase per time'//   &
     & 'step) of detailed crop growth model'
      call writehead (om1,1,filnam,filtext,project)

      write (om1,100)
 100    format ('*',/,                                                  &
     & '*             day     day      -   grCd   kg/ha   kg/ha',       &
     & '   kg/ha   kg/ha   kg/ha   kg/ha   kg/ha       -   kg/ha',/,    &
     & '      Date, Daynr, Daycrp,   DVS,  TSUM,   gass,   mres,OMroots'&
     & ,',OMleaves,OMstems,OMstorage, dmi,    cvf,OMcheck')
      return

      case (2)

! --- write dynamic data ----------------------------------------------------
      if(cvf.lt.nihil) then
        OMroot    = 0.0d0 
        OMleaves  = 0.0d0
        OMstems   = 0.0d0
        OMstorage = 0.0d0
      else
        OMroot    = fr*dmi/cvf
        OMleaves  = fl*(1.0d0-fr)*dmi/cvf
        OMstems   = fs*(1.0d0-fr)*dmi/cvf
        OMstorage = fo*(1.0d0-fr)*dmi/cvf
      endif  

      write (om1,200) date,comma,nint(t),comma,daycrop,comma,dvs,comma, &
     & tsum,comma,gass,comma,mres,comma,OMroot,comma,OMleaves,comma,    &
     & OMstems,comma,OMstorage,comma,dmi,comma,cvf,comma,ccheck
 200  format (a11,a1,i5,a1,i7,a1,f6.2,a1,f6.0, 9(a1,f7.2) )

      return

      case (3)
! --- close crop output file ------------------------------------------------

      close (om1)

      case default
         call fatalerr ('OutbalCropOM1', 'Illegal value for TASK')
      end select

      return
      end 

      subroutine outbalcropom2(task,pathwork,outfil,project,date,       &
     &         daycrop,t,dvs,tsum,storagediff,wlv,wst,wso,wrt,      &
     &         delt,grlv,grst,grso,grrt,drlv,drst,drso,drrt,ombalan)
! ----------------------------------------------------------------------
!     Date               : March 2015  
!     Purpose            : open and write crop output OM balance files 
! ----------------------------------------------------------------------
      implicit none

! --- global variables ------------------
      character(len=11) date
      character(len=*) outfil,pathwork,project
      integer task,daycrop
      real(8) t,dvs,tsum,storagediff,wlv,wst,wso,wrt,delt
      real(8) grlv,grst,grso,grrt,drlv,drst,drso,drrt,ombalan

! --- local variables ------------------
      character(len=1) comma
!      character(len=200) messag
      integer   getun,om2
      character(len=160) filnam,filtext

      save    om2

      comma = ',' 

      select case (task)
      case (1)

! === open output file and write headers =====================

      filnam = trim(pathwork)//trim(outfil)//'.om2'
      om2 = getun (20,90)
      call fopens(om2,filnam,'new','del')
      filtext = 'output of OM2-balance (kg/ha DM, cumulative and '//    &
     & 'increments) of detailed crop growth model'
      call writehead (om2,1,filnam,filtext,project)

      write (om2,100)
 100    format ('*',/,                                                  &
     & '*             day     day      - degday   kg/ha   kg/ha ',      &
     & '  kg/ha   kg/ha   kg/ha   kg/ha   kg/ha   kg/ha   kg/ha ',      &
     & '  kg/ha kg/ha/d kg/ha/d kg/ha/d kg/ha/d',/,                     &
     & '      Date, Daynr, Daycrp,   DVS,  TSUM,storagediff,wlv,',      &
     & '    wst,    wso,    wrt,   grlv,   grst,   grso,   grrt,',      &
     & '   drlv,   drst,   drso,   drrt,ombalan')

      return

      case (2)

! --- write dynamic data ----------------------------------------------------
      write (om2,200) date,comma,nint(t),comma,daycrop,comma,dvs,comma, &
     & tsum,comma,storagediff,comma,wlv,comma,wst,comma,wso,comma,      &
     & wrt,comma,grlv*delt,comma,grst*delt,comma,grso*delt,comma,       &
     & grrt*delt,comma,drlv*delt,comma,drst*delt,comma,drso*delt,       &
     & comma,drrt*delt,comma,ombalan
 200  format (a11,a1,i5,a1,i7,a1,f6.2,a1,f6.0, 14(a1,f7.1) )

      return

      case (3)
! --- close crop output file ------------------------------------------------

      close (om2)

      case default
         call fatalerr ('OutbalCropOM2', 'Illegal value for TASK')
      end select

      return
      end 


      SUBROUTINE CHCKPRT(DVS,FR,FL,FS,FO,FBL)       
! ----------------------------------------------------------------------
!     Last modified      : April 2015
!       based on routines needed for LINTUL4 model (May 2011, Joost Wolf)
!       and extra routines for Wofost (August 2012, Iwan Supit)
!
!     Purpose            : Checks the partitioning factors, and interrupt in case of error
!     Interface parameters, class: I=input,O=output,I/O=input/output
!     class type parameter description (unit)
!       I    R8  FR        Fraction of total dry matter partitioned to the roots (-)
!       I    R8  FL        Fraction of total dry matter partitioned to the leaves (-)
!       I    R8  FS        Fraction of total dry matter partitioned to the stems (-)
!       I    R8  FO        Fraction of total dry matter partitioned to the storage organs (-)
!       I    R8  FBL       Fraction of total dry matter partitioned to the bulbs (-)
! ----------------------------------------------------------------------
      implicit none
! --- global
      real(8)   DVS,FR,FL,FS,FO,FBL
! --- local
      character(len=300) messag
      real(8)   FCHECK

!*     check on partitioning
      fcheck = fr+(fbl+fl+fs+fo)*(1.0d0-fr) - 1.0d0
      if (dabs (fcheck).gt.0.0001d0) then
!        write (messag,'(a,f5.2,/,3(a,g12.5),/,2(a,g12.5))')             &
        write (messag,'(a,f5.2,3(a,g12.5),2(a,g12.5))')                 &
     &      ' error in partitioning functions, dvs= ',dvs,              &
     &      ' fcheck = ',fcheck,' fr = ',fr,' fl = ',fl,                &
     &      ' fs = ',fs,' fo = ',fo
        call fatalerr ('wofost',messag)
      end if
      
      return
      end
 
      subroutine chckcbl(dvs,cvf,dmi,fr,fl,fs,fbl,fo,gass,mres,ccheck)
! ----------------------------------------------------------------------
!     Last modified      : July 2017
!       based on routines needed for LINTUL4 model (May 2011, Joost Wolf)
!       and extra routines for Wofost (August 2012, Iwan Supit)
!
!     Purpose            : Checks the carbon balance, and interrupt in case of error
!     Interface parameters, class: I=input,O=output,I/O=input/output
!     class type parameter description (unit)
!       I    R8  DVS       development stage (-)
!       I    R8  CVF       conversion factor (-)
!       I    R8  DMI       dry matter increase (-)
!       I    R8  FBL       Fraction of total dry matter partitioned to the bulbs (-)
!       I    R8  FR        Fraction of total dry matter partitioned to the roots (-)
!       I    R8  FL        Fraction of total dry matter partitioned to the leaves (-)
!       I    R8  FS        Fraction of total dry matter partitioned to the stems (-)
!       I    R8  FO        Fraction of total dry matter partitioned to the storage organs (-)
!       I    R8  GASS      Gross assimilation (-)
!       I    R8  MRES      maintenance respiration (-)
! ----------------------------------------------------------------------
      implicit none
! --- global
      real(8)   dvs,cvf,dmi,fbl,fr,fl,fs,fo,gass,mres,ccheck
! --- local
      character(len=300) messag

!     check on c-balance
      ccheck = (gass-mres-(fr+(fbl+fl+fs+fo)*(1.0d0-fr))*dmi/cvf)       &
     &       /max (0.0001d0,gass)
      if (dabs (ccheck).gt.0.0001d0) then
        write (messag,'(a,i3,/,3(a,g12.5),/,a,4g12.5,/,2(a,g12.5))')    &
     &     ' carbon flows nog balanced on day ',dvs,                    &
     &     ' ccheck = ',ccheck,' gass = ',gass,' mres = ',mres,         &
     &    ' fr,fbl,l,s,o = ',fr,fbl,fl,fs,fo,' dmi = ',dmi,' dvf = ',cvf
        call fatalerr ('wofost',messag)
      end if    
      return
      end

      subroutine relgrwt(dmi,fr,fl,fs,fo,grrt,grlv,grst,grso,admi)
! ----------------------------------------------------------------------
!     Last modified      : April 2015
!       based on routines needed for LINTUL4 model (May 2011, Joost Wolf)
!       and extra routines for Wofost (August 2012, Iwan Supit)
!
!     Purpose            : To calculate relative growth rate of roots, stems leaves
!                          and storage organs  
! Interface parameters, class: I=input,O=output,I/O=input/output
! ===== ==== =======  =============================================  ==============
! Class Type Name     Description                                    Unit
! ===== ==== =======  =============================================  ==============
!   I    R8  DSLV     Death rate leaves                              kg ha-1 d-1 DM
!   I    R8  DMI      total dry matter increase                      kg ha-1 d-1 DM
!   I    R8  FR       Fraction of total dry matter partitioned to the roots (-)
!   I    R8  FL       Fraction of total dry matter partitioned to the leaves (-)
!   I    R8  FS       Fraction of total dry matter partitioned to the stems (-)
!   I    R8  FO       Fraction of total dry matter partitioned to the storage organs (-)
!   O    R8  GRRT     Growth rate roots                              kg ha-1 d-1 DM
!   O    R8  GRLV     Growth rate leaves                             kg ha-1 d-1 DM
!   O    R8  GRST     Growth rate stems                              kg ha-1 d-1 DM
!   O    R8  GRSO     Growth rate storage organs                     kg ha-1 d-1 DM
!   O    R8  ADMI     total above ground dry matter increase         kg ha-1 d-1 DM
! ===== ==== =======  =============================================  ==============
!   O    R8  ADMI      above ground dry matter increase (-)
      implicit none
! --- global
      real(8)   DMI,FR,FL,FS,FO,GRRT,GRLV,GRST,GRSO,ADMI
! --- local
!     save
            
      ADMI = (1.0d0-FR)*DMI
      GRRT = FR*DMI
      GRLV = FL*ADMI
      GRST = FS*ADMI
      GRSO = FO*ADMI
      
      RETURN
      END

      subroutine deaths(flcropnut,wlv,kdif,lai,NNI,perdl,rdrns,         &
     &                  reltr,dslv)
! ----------------------------------------------------------------------
!     Last modified      : April 2015
!       based on routines needed for LINTUL4 model (May 2011, Joost Wolf)
!       and extra routines for Wofost (August 2012, Iwan Supit)
!
!     Purpose            : Compute the relative death rate leaves due
!                          to stress (kg DM ha-1 d-1)   
!
!     Interface parameters, class: I=input,O=output,I/O=input/output
! ===== ==== =======  =============================================  ==============
! Class Type Name     Description                                    Unit
! ===== ==== =======  =============================================  ==============
!   I    L   flCropNut Flag indicating simulation of nutrient stress -
!   I    R8  WLV       Dry weight of living leaves                   kg ha-1 DM
!   O    R8  reltr
!   I    R8  KDIF      Extinction coeff. for diffuse visible light   -
!   I    R8  LAI       Leaf Area Index
!   I    R8  NNI  
!   I    R8  PERDL     Max.rel.death rate of leaves due to water strs -
!   I    R8  RDRNS     Max.rel.death rate of leaves due to nitrogen strs -
!   O    R8  DSLV      Death rate leaves                             kg ha-1 d-1 DM
! ===== ==== =======  =============================================  ==============
!   -    R8  DSLV1     Death rate leaves due to water stress         kg ha-1 d-1 DM
!   -    R8  DSLV1     Death rate leaves due to self-shading         kg ha-1 d-1 DM
!   -    R8  LAICR     
! ===== ==== =======  =============================================  ==============
      implicit none
! --- global
      logical   flCropNut
      real(8)   DSLV,KDIF,LAI,NNI,PERDL,RDRNS,reltr,WLV
! --- local
      real(8)   DSLV1,DSLV2,LAICR
!      save     
      
!     death rate of leaves due to water stress
      DSLV1 = WLV*(1.d0-reltr)*PERDL
      
!     death rate of leaves due high LAI
      LAICR = 3.2d0/KDIF
      DSLV2 = WLV*max(0.0d0, min(0.03d0, (0.03d0*(LAI-LAICR)/LAICR)))
      DSLV  = MAX (DSLV1, DSLV2)
      
!     death rate increase due to nutrient shortage
      IF(flCropNut .AND. NNI.LT.1.0d0) THEN
         DSLV = DSLV + WLV*RDRNS * (1.0d0-NNI)
      END IF 
      
      RETURN
      END
      
      subroutine deatha(dslv,delt,ilvold,lv,lvage,span,i1,dalv)
! ----------------------------------------------------------------------
!     Last modified      : April 2015
!       based on routines needed for LINTUL4 model (May 2011, Joost Wolf)
!       and extra routines for Wofost (August 2012, Iwan Supit)
!
!     Purpose            : To compute the relative death rate leaves due              * 
!                          to ageing (kg DM ha-1 d-1)   

!     Interface parameters, class: I=input,O=output,I/O=input/output
! ===== ==== =======  =============================================  ==============
! Class Type Name     Description                                    Unit
! ===== ==== =======  =============================================  ==============
!   I    R8  DELT      
!   I    R8  REST   
!   I    I   i1vold      
!   I    R8  SPAN  
!   I    R8  LVAGE
!   I/O  I   11      
!   O    R8  DALV      Death rate leaves due to ageing (kg/ha/d DM)
! ===== ==== =======  =============================================  ==============
      implicit none
! --- global
      integer   i1,ilvold
      real(8)   DSLV,DELT,SPAN,DALV,lv(366),lvage(366)
! --- local
      real(8)   REST
            
! --- first: leaf death due to water stress or high lai is imposed on array
! ---        until no more leaves have to die or all leaves are gone

      rest = dslv*delt
      i1 = ilvold

      do while (rest.gt.lv(max(i1,1)).and.i1.ge.1)
        rest = rest-lv(i1) 
        i1 = i1-1
      enddo

! --- then: check if some of the remaining leaves are older than span,
! ---       sum their weights

      dalv = 0.0d0
      if (lvage(max(i1,1)).gt.span .and. rest.gt.0.0d0 .and.i1.ge.1)then
        dalv = lv(i1)-rest
        rest = 0.0d0
        i1 = i1-1
      endif

      do while (i1.ge.1.and.lvage(max(i1,1)).gt.span)
        dalv = dalv+lv(i1)
        i1 = i1-1
      enddo

      dalv = dalv/delt
      
      RETURN
      END      

      SUBROUTINE GLAI(Fstress,LAIEXP,GLAIEX,TEMP,TBASE,RGRLAI,GRLV,SLAT,GLA)
! ----------------------------------------------------------------------
!     Last modified      : April 2015
!       based on routines needed for LINTUL4 model (May 2011, Joost Wolf)
!       and extra routines for Wofost (August 2012, Iwan Supit)
!
!     Purpose            : exponential, sink limited leave increase
!                          and adjustment of specific leaf area of youngest 
!                          leaf class. Adjust for water and nitrogen stress
!
!     Interface parameters, class: I=input,O=output,I/O=input/output
! ===== ==== =======  =============================================  ==============
! Class Type Name     Description                                    Unit
! ===== ==== =======  =============================================  ==============
!   I    R8  Fstress
!   I/O  R8  LAIEXP
!   I/O  R8  GLAIEX
!   O    R8  GLA       
!   I/O  R8  SLAT     
! ===== ==== =======  =============================================  ==============
! ===== ==== =======  =============================================  ==============
      implicit none
! --- global
      real(8)   Fstress,LAIEXP,GLAIEX,TEMP,TBASE,RGRLAI,GRLV,SLAT,GLA
! --- local
      real(8)   DTEFF,GLASOL
!      save     

      IF (LAIEXP.LT.6.0d0) THEN
         DTEFF  = MAX (0.d0,TEMP-TBASE)
         GLAIEX = Fstress * LAIEXP*RGRLAI*DTEFF
!*        source-limited increase in leaf area
         GLASOL = GRLV*SLAT
!*        sink-limited increase in leaf area
         GLA    = MIN (GLAIEX, GLASOL)
!*        adjustment of specific leaf area of youngest leaf class
         IF (GRLV.GT.0.d0) SLAT = GLA/GRLV
      END IF
      
      RETURN
      END
      SUBROUTINE LVDTH(DELT,DSLV,SPAN,ILVOLD,LVAGE,LV,I1)
! ----------------------------------------------------------------------
!     Last modified      : April 2015
!       based on routines needed for LINTUL4 model (May 2011, Joost Wolf)
!       and extra routines for Wofost (August 2012, Iwan Supit)
!
!     Purpose            : Impose leave death on LV array 
!
!     Interface parameters, class: I=input,O=output,I/O=input/output
! ===== ==== =======  =============================================  ==============
! Class Type Name     Description                                    Unit
! ===== ==== =======  =============================================  ==============
!   I    R8  DELT
!   I    R8  DSLV       
!   I    I   SPAN
!   I    I   ILVOLD
!   I    R8  LVAGE    
!   I/O  R8  LV     
!   O    I   I1
! ===== ==== =======  =============================================  ==============
! ===== ==== =======  =============================================  ==============
      implicit none
! --- global
      integer   ilvold,i1
      real(8)   delt,dslv,lv(366),lvage(366),span
! --- local
      real(8)   dslvt

      save     

! --- remaining leaves
      dslvt = dslv*delt
      i1 = ilvold
      do while ((dslvt.gt.0.0d0) .and. (i1.ge.1))
        if (dslvt.ge.lv(i1)) then
          dslvt = dslvt-lv(i1)
          lv(i1) = 0.0d0
          i1 = i1-1
        else
          lv(i1) = lv(i1)-dslvt
          dslvt = 0.0d0
        endif
      enddo

! --- leaves older than span die
      do while (lvage(max(i1,1)).gt.span.and.i1.ge.1)
        lv(i1) = 0.0d0
        i1 = i1-1
      enddo

      RETURN
      END

      SUBROUTINE LVSHFT(DELT,FYSDEL,ILVOLD,grlv,slat,LV,LVAGE,SLA)
! ----------------------------------------------------------------------
!     Last modified      : April 2015
!       based on routines needed for LINTUL4 model (May 2011, Joost Wolf)
!       and extra routines for Wofost (August 2012, Iwan Supit)
!
!     Purpose            : Shift contents of LV, LVAGE, SLA tables with one day  
!
!     Interface parameters, class: I=input,O=output,I/O=input/output
! ===== ==== =======  =============================================  ==============
! Class Type Name     Description                                    Unit
! ===== ==== =======  =============================================  ==============
!   I    R8  DELT
!   I    R8  FYSDEL       
!   I    I   ILVOLD
!   I/O  R8  LV     
!   I/O  R8  LVAGE     
!   I/O  R8  SLA     
! ===== ==== =======  =============================================  ==============
! ===== ==== =======  =============================================  ==============
      implicit none
! --- global
      integer   ilvold
      real(8)   delt,fysdel,grlv,slat,lv(366),lvage(366),sla(366)
! --- local
      integer   i1

!      save     

!     Shift contents of LV, LVAGE, SLA tables with one day
      do i1 = ilvold,1,-1
        lv(i1+1) = lv(i1)
        sla(i1+1) = sla(i1)
        lvage(i1+1) = lvage(i1)+fysdel*delt
      enddo

!     new leaves in class 1
      lv(1) = grlv*delt
      sla(1) = slat
      lvage(1) = 0.0d0 
     
      return
      end
      SUBROUTINE LVWGLI(ILVOLD,LV,SLA,LASUM,WLV)
! ----------------------------------------------------------------------
!     Last modified      : April 2015
!       based on routines needed for LINTUL4 model (May 2011, Joost Wolf)
!       and extra routines for Wofost (August 2012, Iwan Supit)
!
!     Purpose            : Calculate new wlv and lai 
!
!     Interface parameters, class: I=input,O=output,I/O=input/output
! ===== ==== =======  =============================================  ==============
! Class Type Name     Description                                    Unit
! ===== ==== =======  =============================================  ==============
!   I    I   ILVOLD
!   I    R8  LV     
!   I    R8  SLA     
!   O    R8  LASUM     
!   O    R8  WLV     
! ===== ==== =======  =============================================  ==============
! ===== ==== =======  =============================================  ==============
      implicit none
! --- global
      integer   ilvold
      real(8)   lv(366),sla(366),lasum,wlv
! --- local
      integer    i1
!      save     

! --- calculation of new leaf area and weight
      lasum = 0.0d0
      wlv = 0.0d0
      do i1 = 1,ilvold
        lasum = lasum+lv(i1)*sla(i1)
        wlv = wlv+lv(i1)
      enddo
  
      return
      end

! --- for soybean (swsoybean=1): temperature and photoperiodicity
      subroutine mgtemprf(tav,toptdvr,tmindvr,tmaxdvr,rfmgtemp)
! ----------------------------------------------------------------------
!     Last modified      : Sept 2015
!       based on routines needed from pyWofost for soybean (Allard de Wit, 2015)
!
!     Purpose            : temperature reduction factor for soybean (short day)
!       approach and parameters based on Setiyono et al. doi 10.1016/j.fcr.2006.07.011
!       http://digitalcommons.unl.edu/agronomyfacpub/112
!
!     Interface parameters, class: I=input,O=output,I/O=input/output
! ===== ==== =======  =============================================  ==============
! Class Type Name     Description                                    Unit
! ===== ==== =======  =============================================  ==============
!   I    R8  tav        
!   I    R8  toptdvr       
!   I    R8  tmindvr
!   I    R8  tmaxdvr     
!   O    R8  rfmgtemp    
! ===== ==== =======  =============================================  ==============
! ===== ==== =======  =============================================  ==============
      implicit none
! --- global
      real(8)   tav, toptdvr, tmindvr, tmaxdvr, rfmgtemp
! --- local
      real(8)   alpha, p1, p2, p3, p4

      alpha = log(2.0d0)/(log((tmaxdvr-tmindvr)/(toptdvr-tmindvr)))
      if(tav.lt.tmindvr .or. tav.gt.tmaxdvr) then
        rfmgtemp = 0.0d0
      else
        p1 = 2.0d0 * (tav - tmindvr)**alpha
        p2 = (toptdvr - tmindvr)**alpha
        p3 = (tav - tmindvr)**(2.0d0*alpha)
        p4 = (toptdvr - tmindvr)**(2.0d0*alpha)
        rfmgtemp = (p1 * p2 - p3) / p4

      endif

      return
      end

      subroutine mgphotoprf(mg,iday,lat,popt,pcrt,flphenodayl,rfmgphotop)
! ----------------------------------------------------------------------
!     Last modified      : Sept 2015
!       based on routines needed from pyWofost for soybean (Allard de Wit, 2015)
!
!     Purpose            : Photoperiod reduction factor for soybean (short day)
!       approach and parameters based on Setiyono et al. doi 10.1016/j.fcr.2006.07.011
!       http://digitalcommons.unl.edu/agronomyfacpub/112
!
!     Interface parameters, class: I=input,O=output,I/O=input/output
! ===== ==== =======  =============================================  ==============
! Class Type Name     Description                                    Unit
! ===== ==== =======  =============================================  ==============
!   I    R8  mg       maturity group 
!   I    I   iday     daynr      
!   I    R8  lat      lattitude
!   I    L   flphenodayl Flag to allow input of POPT and PCRT or using 
!                        empirical relation from Setiyono et al
!   I    R8  popt     optimal daylength for phenological developm.   hr
!   I    R8  pcrt     critical daylength for phenological developm.  hr
!   O    R8  rfmgphotop  reduction factor for photoperiodicity       -  
! ===== ==== =======  =============================================  ==============
! ===== ==== =======  =============================================  ==============
      implicit none
! --- global
      integer   iday
      logical   flphenodayl
      real(8)   lat, popt, pcrt
      real(8)   mg, rfmgphotop
! --- local
      real(8)   alpha,m,p0,p1,p2
      real(8)   dec,daylp,pi,rad,sinld,cosld,aob
      parameter (pi=3.1415926d0, rad=0.0174533d0)

! astronomic daylength according to solar elevation angle of -0.833 day
!*     Declination and solar constant for this day
      dec = -asin(dsin(23.45d0*rad)*dcos(2.d0*pi*dble(iday+10)/365.0d0))
      sinld = dsin(rad*lat)*dsin(dec)
      cosld = dcos(rad*lat)*dcos(dec)
      aob = sinld/cosld
      if (dabs(aob).le.1.0d0) then 
         daylp = 12.0d0*(1.d0+2.d0*asin(aob)/pi)
      else
         if (aob.gt.1.0d0)  daylp = 24.0d0
         if (aob.lt.-1.0d0) daylp =  0.0d0
      endif

! First determine Popt and Pcrt based on maturity group rating
      m = 3.0d0
      if(flphenodayl) then
        continue            ! popt and pcrt  are input
      else
        popt = 12.759d0 - 0.388d0*mg - 0.058d0*mg**2
        pcrt = 27.275d0 - 0.493d0*mg - 0.066d0*mg**2
      endif
      alpha = log(2.0d0)/log(((pcrt - popt)/m) + 1.0d0)
      p0 = (pcrt - popt)/m

      if(daylp.lt.popt) then
        rfmgphotop = 1.0d0
      else
        if (daylp.gt.pcrt) then
          rfmgphotop = 0.0d0
        else
          p1 = (daylp - popt)/m + 1.0d0
          p2 = (pcrt - daylp)/(pcrt - popt)
          rfmgphotop = (p1*(p2**p0))**alpha
        endif
      endif

!      write(99,*)  iday, daylp, rfmgphotop
      
      return
      end

      subroutine sumttd(task,flGrassGrowth,dateGrassGrowth)
! ----------------------------------------------------------------------
!     Last modified      : Jan 2016
!     Author             : Joop Kroes
!
!     Purpose            : Suppress grass growth as long as 3 criteria are not met: 
!                          temperature, time and depth
!
!     Interface parameters, class: I=input,O=output,I/O=input/output
!     class type parameter description (unit)
!       I    C   condition case: 'initial' or 'dynamic' (-)
!       I    R8  tsumtemp  temperature limit to initiate grass growth  [0.0..20.0 grC, R]
!       I    I   tsumtime  time (nrs of sequential days) with temp above tsumtemp for grass growth [1..20 days, I]
!       I    R8  tsumdepth depth at which temp above tsumtemp for grass growth [0.0..100.0 cm below soil surface, R]
!       I    R8  z         depth of a node (L)
!       I    R8  tsoil     Array with soil temperatures (oC) for each compartment
!       O    L   flGrassGrowth flag indicating grass growth (suppressed=.false. when criteria are not met) [.true .or. .false. -, L]
! ----------------------------------------------------------------------
      use Variables
      implicit none

! --- global
      character(len=*), intent(in) :: task
      logical, intent(out)         :: flGrassGrowth      ! flag indicating grass growth (suppressed=.false. when criteria are not met) [.true .or. .false. -, L]

! --- local
      integer    :: tsumtimecum      ! cumulative, from 1-jan, time (nrs of sequential days) with temp above tsumtemp for grass growth [1..20 days, I]
      logical    :: fltsumtemp       ! flag to indicate if temperature criteria is met 
      logical    :: fltsimprev       ! flag to indicate if temperature criterion is met during simulated previous day
      logical    :: fltsimcount      ! flag to indicate nr of contiuous simulated days that temperature criteria is met
      integer    :: cmpcrit, node
      character(len=11), intent(out) ::  dateGrassGrowth            ! date of start of GrassGrowth
!     output
      integer    :: uo, getun   !, idum, ios
      character(len=160)  :: filnam, filtext
      character(len=1)  :: comma
      logical    :: flexist, flopened
      
      save

      comma = ',' 

      
      select case(task)

      case('initial')
          tsumtimecum = 0
          fltsumtemp = .false.
          fltsimprev = .false.
          ! find depth of compartment with critical soil temperature
          cmpcrit = 1
          do node = 1, numnod
             if (z(node) .le. (-1.0d0*tsumdepth)) then
                cmpcrit = node
                exit
             endif
          enddo
          ! no growth as long as 3 criteria are not met
          flGrassGrowth = .false.
          dateGrassGrowth = 'undefined'
          
          ! === open output file and write headers
          filnam = trim(pathwork)//trim(outfil)//'.ttd'
          flopened = .false.
          if(uo.gt.0) then
              !  Inquiry by Unit
              inquire (uo, opened=flopened, exist=flexist)
          endif
          if (.not.flexist .and. .not.flopened) then
              uo = getun (20,90)
              call fopens(uo,filnam,'new','del')
              filtext = 'output of subr sumttd'
              call writehead (uo,1,filnam,filtext,project)
              write (uo,100)
 100          format (' Date,z(cmpcrit),tsoil(cmpcrit),',               &
     &           'fltsumtemp,fltsimprev,fltsimcount')
          endif

      case('dynamic')
          ! temperature and depth criterium
          if(tsoil(cmpcrit).ge.tsumtemp) then
              fltsumtemp = .true.
          else
              fltsumtemp = .false.
          endif
          ! timing criterium: set flag for continuous days that exceed critical temperature
          if(fltsumtemp .and. fltsimprev) then
              tsumtimecum = tsumtimecum + 1
          else
              tsumtimecum = 0
              fltsimprev = .false.
          endif
          !  timing criterium: count nr of days exceeding critical temperature
          if(tsumtimecum.ge.tsumtime) then
              fltsimcount = .true.
          else
              fltsimcount = .false.
          endif
          ! no growth as long as 3 criteria are not met
          flGrassGrowth = .false.
          if(fltsumtemp .and. fltsimprev .and. fltsimcount) then
              call dtdpst ('year-month-day',t1900,dateGrassGrowth)
              flGrassGrowth = .true.
          endif

          ! timing criteria for next timestep
          if(fltsumtemp) then
              fltsimprev = .true.
          endif
          
          ! === write output 
          write (uo,200) Date,comma,z(cmpcrit),comma,tsoil(cmpcrit),    &
     &           comma,fltsumtemp,comma,fltsimprev,comma,fltsimcount
 200      format (a11,2(a1,f7.2),3(a1,i3))

      end select

      return
      end
