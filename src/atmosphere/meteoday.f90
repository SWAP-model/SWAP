! File VersionID:
!   $Id: meteoday.f90 372 2018-03-13 10:01:20Z heine003 $
!
!     This file contains the following subroutines, in order of calling:
!     1. MeteoDay          : main routine                       ; called in SWAP
!     2. ResetMetFlx       : resets meteorological fluxes       ; called in MeteoDay
!     3. ProcessMeteoDays  : processes meteo input data per day ; called in MeteoDay
!     4-6.                 : calculate interception according to:
!     4. VonHHBraden         - Von Hoyningen-Hune and Braden    ; called in ProcessMeteoDays (optional)
!     5. Gash                - Gash                             ; called in ProcessMeteoDays (optional)
!     6. Ruttervw            - Rutter (adapted)                 ; called in ProcessMeteoDays (optional)
!     7. DivIntercep       : divides interception into rain part
!                                              & sprinkling part; called in ProcessMeteoDays
!     8. Reduceva          : calc. reduction of soil evaporation; called in ProcessMeteoDays, ProcessMeteoTsteps, ETSine

module MeteoVars
!  for availability in both ReadMeteoDay and ProcessMeteoDay
   integer   count,first,i,irecord,last,ndayparts
   real(8)   arain(96),awind(96),restint,interc,Edirectpond
   real(8)   aintc,dttp,eintc,etr,gctp,hum,netrainflux,rainflux
   real(8)   sumtav,svp,wfrac,win,Edirect,Tdirect,Tdirectwet
end module MeteoVars


! SUBROUTINE 1.
! ----------------------------------------------------------------------
      subroutine ReadMeteoDay
! ----------------------------------------------------------------------
!     Last modified      : March 2014
!     Purpose            : returns meteorological fluxes of current day
!                        : or of parts of a day (detailed meteo input)
! ----------------------------------------------------------------------
      use variables
      use MeteoVars
      implicit none

! --- local
      character(len=11)  detdate
      character(len=3)   ext
      character(len=200) filnam
      character(len=300) messag

! ----------------------------------------------------------------------

      call ResetMetFlx ()

! 1: Check whether meteo data are available of today; pass on weather of today
! 1.0 Daily Meteo 0000000000000000000000000000000000000000000000000000000 Daily Meteo
!
      if (swmetdetail.eq.0) then

!   - Check availability of meteo data of today
        if (daymeteo.lt.daynrfirst .or. daymeteo.gt.daynrlast) then
          messag ='In meteo file no meteo data are'//                   &
     &    ' available for '//date//'. First adapt meteo file!'
          call fatalerr ('meteo',messag)
        end if

!   - Pass on weather values of today
        rad  = arad(daymeteo+1-daynrfirst)
        tmn  = atmn(daymeteo+1-daynrfirst)
        tmx  = atmx(daymeteo+1-daynrfirst)
        hum  = ahum(daymeteo+1-daynrfirst)
        win  = awin(daymeteo+1-daynrfirst)
        grai = arai(daymeteo+1-daynrfirst)
        etr  = aetr(daymeteo+1-daynrfirst)

!   - If hum is missing or tav cannot be calculated: set rh at -99.0
        rh = 1.0d0
        if (hum.lt.-98.0d0 .or. tmn.lt.-98.0d0 .or. tmx.lt.-98.0d0)     &
     &      rh=-99.0d0

!   - Calculate 24h average temperature
        tav = (tmx+tmn)*0.5d0
!   - Calculate average day temperature
        tavd = (tmx+tav)*0.5d0

        if (rh.ge.-98.0d0) then
!   - Calculate saturated vapour pressure [kpa]
        svp = 0.3055d0*(dexp(17.27d0*tmn/(tmn+237.3d0)) +               &
     &                  dexp(17.27d0*tmx/(tmx+237.3d0)))
!   - Calculate relative humidity [fraction]
          rh = min(hum/svp,1.0d0)
        endif

!   - CFO file for PEARL: save meteo variables of today for output
        out_rad = real(rad)
        out_tmn = real(tmn)
        out_tmx = real(tmx)
        out_hum = real(hum)
        out_win = real(win)
        out_etr = real(etr)*0.001
        if (swrain.eq.2) then
          out_wet = real(wet(daymeteo+1-daynrfirst))
        else
          out_wet = -1.0
        endif
!
! end 1 Daily Meteo 00000000000000000000000000000000000000000000000000000 Daily Meteo
!!
! 1.1 Detailed Meteo 1111111111111111111111111111111111111111111111111111 Detailed Meteo
!
      elseif (swmetdetail.eq.1) then

!   - Check availability of meteo data of today
!     + compose filename meteorological file for use in warnings
        write (ext,'(i3.3)') mod(yearmeteo,1000)
        filnam = trim(pathatm)//trim(metfil)//'.'//trim(ext)

        do i = 1, nmetdetail
          irectotal = irectotal + 1
          if (i .ne. detrecord(irectotal)) then
            messag='In meteo file '//trim(filnam)//' record number(s)'//&
     &      ' are not correct at '//date//'. First adapt meteo file!'
            call fatalerr ('meteo',messag)
          end if
          call dtdpst('year-month-day',                               &
     &                 dettime(irectotal)+0.1d0,detdate)
          call dtdpst('year-month-day',t1900+0.1d0,date)
          if (detdate .ne. date) then
            messag ='In meteo file '//trim(filnam)//' the amount of '// &
     &      'records deviate near '//date//'. First adapt meteo file!'
            call fatalerr ('meteo',messag)
          end if

!   - Pass on weather records of today
          arad(i)  = detrad(irectotal)
          ahum(i)  = dethum(irectotal)
          atav(i)  = dettav(irectotal)
          awind(i) = detwind(irectotal)
          arain(i) = detrain(irectotal) * 0.1d0 ! convert from mm to cm
        enddo
      endif
!
! end 1 Detailed Meteo 11111111111111111111111111111111111111111111111111 Detailed Meteo
! end 1.

! 2: Rain and Snow
! 2.0 Daily Meteo 0000000000000000000000000000000000000000000000000000000 Daily Meteo
!
      if (swmetdetail.eq.0) then
!   - Rain: convert daily sum of precipitation from mm to cm
        grai = grai * 0.1d0

!   - Snow: determine whether the precipitations falls as rain or snow
        if (swsnow.eq.1) then
!     >>> snow calculations
          if (tav.gt.TePrRain) then
            gsnow = 0.0d0
          elseif (tav.lt.TePrSnow) then
            gsnow = grai
          else
!     + distribution over snow and rain according to linear relation between
!       the two transition temperatures TePrRain and TePrSnow
            gsnow = grai * (TePrRain-tav) / (TePrRain-TePrSnow)
          endif

!     + rain fall rate on snowpack
          if (ssnow.gt.1.0d-6) then
            snrai = grai - gsnow
          else
            snrai = 0.0d0
          endif

!     + fraction of precipitation that is not snow nor rain on snowpack
          if (grai.gt.0.d0) then
             fprecnosnow = 1.0d0 - (gsnow + snrai) / grai
          else
             fprecnosnow = 0.0d0
          endif
        else
!     >>> no snow calculations
          fprecnosnow = 1.0d0
        endif
!
! end Daily Meteo 0000000000000000000000000000000000000000000000000000000 Daily Meteo
!!
! 2.1 Detailed Meteo 1111111111111111111111111111111111111111111111111111 Detailed Meteo
!
      elseif (swmetdetail.eq.1) then
!   - Total daily sum of precipitation (cm) today
        grai = 0.0d0
        do i = 1, nmetdetail
          grai = grai + arain(i)
        enddo
!   - Remaining amount of interception
        restint = 0.0d0

!   - Combination of detailed meteo and snow is not included:
!     omit snow and rain fall rate on snowpack
        gsnow = 0.0d0
        ssnow = 0.0d0
        snrai = 0.0d0
      endif
!
! end 2 Detailed Meteo 11111111111111111111111111111111111111111111111111 Detailed Meteo
! end 2.
      return
      end subroutine ReadMeteoDay

! SUBROUTINE 2.
! ----------------------------------------------------------------------
      subroutine ResetMetFlx ()
! ----------------------------------------------------------------------
!     Last modified      : February 2014
!     Purpose            : reset intermediate and cumulative meteorological
!                          fluxes
!     Interface:
!       I   - flzerointr,flzerocumu,caintc,cgrai,cnrai,igrai,inrai,iprec
!       O   - caintc,cgrai,cnrai,igrai,inrai,iprec
! ----------------------------------------------------------------------
      use variables, only: flzerointr,flzerocumu,caintc,cgrai,cnrai,igrai,inrai,iprec
      implicit none

! --- local

! --- reset cumulative intermediate fluxes
      if (flzerointr) then
        iprec = 0.0d0
        igrai = 0.0d0
        inrai = 0.0d0
      endif

! --- reset cumulative meteorological fluxes
      if (flzerocumu) then
        cgrai = 0.0d0
        cnrai = 0.0d0
        caintc = 0.0d0
      endif

      return
      end subroutine ResetMetFlx


! SUBROUTINE 3.
! ----------------------------------------------------------------------
      subroutine ProcessMeteoDay
! ----------------------------------------------------------------------
!     Last modified      : February 2014
!     Purpose            : processing of
! ----------------------------------------------------------------------

      use Variables
      use MeteoVars
      implicit none
      include  'params.fi'
      real(8)  rcs,afgen
      data     rcs/0.15d0/

! 3: Interception: VonHHBraden and Gash
! 3.x Both Daily and Detailed xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Both
!
!   - Calculation of Interception and Net rain & Net irrigation depth [cm]
      if ((lai .lt. 1.d-3) .or. (grai+gird .lt. 1.d-5) .or.             &
     &    (swinter.eq.0) .or. (gsnow.gt.0.0d0) .or.(ssnow.gt.0.0d0))then

!   - No vegetation, rainfall/irrigation or interception calculation
          aintc = 0.d0

      else if (swinter .eq. 1) then
!   - Calculate interception, method Von Hoyningen-Hune and Braden
         call VonHHBraden (aintc)
      else if (swinter .eq. 2) then
!   - Calculate interception, method Gash (1995)
         call Gash (aintc)
      end if

!   - Divide interception into rain part and irrigation part and
!     calculate net rain (nraida) and net sprinkling irrigation (nird)
      if (swinter.ne.3)                                                 &
     &   call DivIntercep (aintc)
!
! end 3 Both xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Both
! end 3.

!
!======================== L O O P  over dayparts =============================
!
!   - Set number of repetitions within one day
      if (swmetdetail.eq.0) then
!   + in case of daily meteo input: only one time per day
         ndayparts = 1
      elseif (swmetdetail.eq.1) then
         ndayparts = nmetdetail
      endif

      do 1000 irecord = 1, ndayparts
!
! 4: Calculate evapotranspiration: et0, ew0, es0
! 4.x Both Daily and Detailed xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Both
!
!   - Reference evapotranspiration has been specified
         if (swmetdetail.eq.0 .and. swetr.eq.1) then
           if (.not. flCropEmergence) then
!     + no crop
             et0 = 0.0d0
             ew0 = 0.0d0
             es0 = etr
             if (swcfbs.eq.1) es0 = cfbs*etr
           else
!     + crop is present
             if (swcf.eq.1 .or. swcf.eq.3) then
                et0 = cf*etr
                if (swcf .eq. 1) then
                   ew0 = cf*etr
                else
                   ew0 = cfeic*etr
                endif
             endif
             es0 = etr
             if (swcfbs.eq.1) es0 = cfbs*etr
           endif
!
!   - Reference evapotranspiration must be calculated
         elseif (swmetdetail.eq.1 .or. swetr.eq.0) then

           if (swmetdetail.eq.1) then
!     + define weather variables of current record
             rad = arad(irecord) / metperiod     ! from j/m2/period to j/m2/d
             tav = atav(irecord)
             hum = ahum(irecord)
             win = awind(irecord)
           endif

!   - Calculate evapotranspiration using Penman-Monteith: et0, ew0, es0 (mm/d)
!          in case of daily meteo(swmetdetail = 0) irecord is always 1
           call PenMon (logf,swscre,daynr,lat,alt,Altw,angstroma,       &
     &          angstromb,rcs,rad,tav,hum,win,rsc,es0,et0,ew0,swcf,ch,  &
     &          flCropEmergence,daylp,flmetdetail,irecord,              &
     &          nmetdetail,albedo,tmn,tmx,rsw,difpp,dsinbe,atmtr,       &
     &          Edirect,Tdirect,Tdirectwet,rsoil,swdivide,kdif,kdir,    &
     &          lai,Edirectpond)

           if (.not. flCropEmergence) then
!     + no crop
             if (swcfbs .eq. 1) then
               if (swcf .eq. 1) then
                 es0 = cfbs*et0
               else
                 es0 = cfbs*es0
               endif
             endif
             et0 = 0.0d0
             if (swmetdetail.eq.1 .and. (swcf.eq.1 .or. swcf.eq.3)) then
               if (swcf.eq.1) then
                 ew0 = cf*ew0
               else
                 ew0 = cfeic*ew0
               endif
             endif
           else
!     + crop is present
             if (swcfbs .eq. 1) then
               if (swcf .eq. 1) then
                 es0 = cfbs*et0
               else
                 es0 = cfbs*es0
               endif
             endif
             if (swcf.eq.1 .or. swcf.eq.3) then
               et0 = cf*et0
               if (swcf.eq.1) then
                 ew0 = cf*ew0
               else
                 ew0 = cfeic*ew0
               endif
             endif
           endif
!
        endif
!
! end 4 Both xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Both
! end 4.

! 5: Interception option NHI, adapted Rutter model
! 5.x Both Daily and Detailed xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Both
!
        if (swinter .eq. 3) then
!   - Set parameter values
           if (swmetdetail.eq.0) then      ! if swmetdetail = 0, siccapact is set in cropgrowth module
             if (croptype(icrop).eq.1 .and. swgc.eq.2) then
               gctp  = gc
             else
               gctp  = 1.0d0 - dexp(-1.0d0*kdir*kdif*lai)
               if (gctp .lt. 1.0d-5) siccapact=0.
             endif
             dttp = 1.0d0                 ! value of 1 d required for the daily meteo option
           elseif (swmetdetail.eq.1) then
             siccapact = afgen(siccaptb,(2*magrs),t)
             if (croptype(icrop).eq.1 .and. swgc.eq.2) then
               gctp  = gc
             elseif (croptype(icrop).eq.1 .and. swgc.eq.1) then
               gctp  = 1.0d0 - dexp(-1.0d0*kdir*kdif*lai)
             endif
             if (gctp .lt. 1.0d-5) siccapact=0.
             dttp = dt
           endif

!   - Calculate interception, method Rutter
           call ruttervw (gctp,aintc,eintc)

!   - Divide interception into rain part and irrigation part and
!     calculate net rain (nraida) and net sprinkling irrigation (nird)
           call DivIntercep(aintc)
         endif
!
! end 5 Both xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Both
! end 5.

! 6: Fraction of the day or period the crop is wet
! 6.x Both Daily and Detailed xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Both
!
!   - Calculate fraction of the day or period the crop is wet
         if (swmetdetail.eq.0) then
!        + fraction of the day the crop is wet
           if (ew0.lt.0.0001d0) then
             wfrac = 0.0d0
           else
             if (ew0.lt.0.0001d0) then
               wfrac = 0.0d0
             else
               if (swinter .ne. 3) then
                 if (swdivide .eq. 0) then
                   wfrac = max(min(aintc*10.0d0/ew0,1.0d0),0.0d0)
                 else
                   if(tdirectwet.gt.nihil) then
                     wfrac = max(min(aintc*10.0d0/tdirectwet,1.0d0),0.0d0)
                   else
                     wfrac = 0.0d0
                   endif
                 endif
               else
                 wfrac = max(min(eintc*10.0d0/ew0,1.0d0),0.0d0)
               endif
             endif
           endif
!     + fraction of the period the crop is wet
         elseif (swmetdetail.eq.1) then
           if (grai .lt. 1.0d-12) then
             interc = 0.0d0
             wfrac  = 0.0d0
           else
             interc = restint + aintc * arain(irecord) / grai
             if (ew0.lt.0.0001d0) then
               wfrac = 0.0d0
             else
               if (swcf.ne.3) then
                 wfrac = max(min(interc*10.0d0/ew0/metperiod,1.d0),0.d0)
               else
                 wfrac = max(min(eintc/ew0,1.0d0),0.0d0)
               endif
             endif
           endif
!   - Remaining amount of interception for swmetdetail = 1
           restint = max(interc - wfrac * metperiod * ew0 * 0.1d0, 0.d0)
         endif
!
! end 6 Both xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Both
! end 6.

! 7: Potential soil evaporation & transpiration
! 7.x Both Daily and Detailed xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Both
!
!   - Potential soil evaporation (peva) [cm/d]
         peva = max(0.0d0, (es0*dexp(-1.0d0*kdir*kdif*lai)*0.1d0))
         if (swcf.ne.3 .or. (swmetdetail.eq.0 .and. swinter.ne.3))      &
     &   peva = max(0.0d0,(1.0d0-wfrac)*peva)

!   - Alternative for peva (simple model, soil cover fraction specified)
         if (flCropCalendar .and. .not.flCropHarvest) then
           if (croptype(icrop).eq.1 .and. swgc.eq.2) then
             peva = (1.0d0-gc)*es0*0.1d0
             if (swcf.ne.3 .or. (swmetdetail.eq.0 .and. swinter.ne.3))  &
     &       peva = (1.0d0-wfrac)*peva
           endif
         endif

!   - Adapt peva in case of ponding
         if (pond .gt. 1.0d-10) then
           if (SwETr.eq.0 .and. es0.gt.1.0d-8) then
             peva = ew0/es0 * peva
           elseif (es0.gt.1.0d-8) then
             if (swcfbs .eq. 1 .and. cfbs .gt. small) then
               peva = cfevappond * peva / cfbs
             else
             peva = cfevappond * peva
           endif
         endif
         endif

!   - Potential soil evaporation [cm/d] according to PMdirect
         if (swdivide .eq. 1) then
           if (pond .gt. 1.0d-10) then
             peva = Edirectpond*0.1d0
           else
             peva = Edirect*0.1d0
           endif
         endif

!   - Potential transpiration (ptra) [cm/d]
         if (swcf .ne. 3) then
           ptra = ((1.0d0-wfrac)*et0-peva*10.0d0)*0.1d0
         else
           ptra = (1.0d0-wfrac)*et0*0.1d0
         endif
         ptra = max(ptra,(1.01d0*nihil))

!   - Potential transpiration [cm/d] according to PMdirect
         if (swdivide .eq. 1) then
           ptra = (1.0d0-wfrac) * Tdirect * 0.1d0
           ptra = max(ptra,(1.01d0*nihil))
         endif

!   - Correction of potential transpiration as a function of atmospheric CO2 concentration
         if (flCO2 .and. flCropEmergence) then
            ptra = fco2tra * ptra
         endif

! end 7 Both xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Both
! end 7.

! 8: results Detailed weather records
! 8.1 Detailed Meteo 1111111111111111111111111111111111111111111111111111 Detailed Meteo
!
!   - Result of detailed weather records (cm/d)
         if (swmetdetail.eq.1) then
           tpot(irecord) = ptra
           epot(irecord) = peva
           if (grai .lt. 1.0d-12) then
             grain(irecord) = 0.0d0
             nrain(irecord) = 0.0d0
           else
             grain(irecord) = arain(irecord) / metperiod
             nrain(irecord) = arain(irecord) / metperiod * nraida /grai
           endif
         endif
!
! end 8 Detailed Meteo 11111111111111111111111111111111111111111111111111 Detailed Meteo
! end 8.
!
1000  continue
!======================== E N D  L O O P  over dayparts ================================


! 9: Actual daily rain and snow fluxes and soil evaporation rate of today for Daily Meteo
! 9.0 Daily Meteo 0000000000000000000000000000000000000000000000000000000 Daily Meteo
!
      if (swmetdetail.eq.0) then

!   - Finterception: ratio net / gross rain flux; net rainflux = gross - interception
         if (grai.gt.1.d-5) then
!     + finterception is exclusively meant for dividing rain flux into interception part and net rain part;
!       not for sprinkler irrigation!
           finterception = nraida / grai
           if (aintc.lt.1.0d-5) finterception = 1.0d0
         else
           finterception = 1.0d0
         endif

!   - In case of daily precipitation sum: set actual gross and net rainflux, and interception on TIMESTEP basis
         if (swrain.eq.0) then
           rainflux    = fprecnosnow * grai
           netrainflux = finterception * rainflux
           graidt  = rainflux
           nraidt  = netrainflux
           if (swuseCN == 1) then
               call CNmethod(2)
               nraidt = nraidt - Runoff_CN
           end if
           aintcdt = rainflux - netrainflux  ! aintcdt involves ONLY interception of RAIN
         endif

!   - Soil evaporation rate of today
         if (.not. fletsine) then
           call reduceva (1,nraida)
         endif

!   - Save daily potential values for use in ETSine
         ptraday = ptra
         pevaday = peva

!   - Calculate atmospheric demand [cm]
         atmdem = et0*0.1d0

      endif
!
! end 9 Daily Meteo 00000000000000000000000000000000000000000000000000000 Daily Meteo
! end 9.

! 10: Set daily weather values for Detailed Meteo
! 10.1 Detailed Meteo 111111111111111111111111111111111111111111111111111 Detailed Meteo
!
      if (swmetdetail.eq.1) then
!   - Average temperature of today
         sumtav = 0.d0
         do i = 1, nmetdetail
           sumtav = sumtav + atav(i)
         enddo
         tav = sumtav * metperiod

!   - Minimum and maximum temperature of today
         tmx = -50.d0
         tmn = 99.d0
         do i = 1, nmetdetail
           tmx = max(tmx,atav(i))
           tmn = min(tmn,atav(i))
         enddo

!   - Calculate saturated vapour pressure [kpa]
        svp = 0.3055d0*(dexp(17.27d0*tmn/(tmn+237.3d0)) +               &
     &                  dexp(17.27d0*tmx/(tmx+237.3d0)))
!   - Calculate relative humidity [fraction]
         rh = min(hum/svp,1.0d0)

!   - Average temperature between 6 and 18 hour
         sumtav = 0.d0
         count = 0
         first = int(0.25/metperiod) + 1
         last = int(0.75/metperiod)
         do i = first, last
           sumtav = sumtav + atav(i)
           count = count + 1
         enddo
         tavd = sumtav / count

!   - Daily radiation (J/m2/d) and atmospheric demand (cm/d)
         rad = 0.d0
         atmdem = 0.d0
         do i = 1,nmetdetail
           rad = rad + arad(i)
           atmdem = atmdem + tpot(i)
         enddo

!   - Fluxes of current time step (start of the day)
         ptra = tpot(1)
         peva = epot(1)
         graidt = grain(1)
         nraidt = nrain(1)
         aintcdt = graidt - nraidt    ! aintcdt involves ONLY interception of RAIN

      endif
!
! end 10 Detailed Meteo 1111111111111111111111111111111111111111111111111 Detailed Meteo
! end 10.

      return

      end subroutine ProcessMeteoDay


! SUBROUTINE 4.
! ----------------------------------------------------------------------
      subroutine VonHHBraden (aintc)
! ----------------------------------------------------------------------
!     Last modified      : July 2012
!     Purpose            : calculate interception, method Von Hoyningen-Hune and Braden
!     Interface:
!       I   - callroutine,logf,grai,gird,kdif,kdir,lai
!       O   - aintc
! ----------------------------------------------------------------------
      use variables, only: grai,gird,kdif,kdir,cofab,lai,isua
      implicit none

! --- global
      real(8)   aintc              ! Amount of rainfall interception during current day (L) in cm/d

! --- local
      real(8)   rpd                ! intercepted precipitation (rain+irrig) in mm
      real(8)   cofbb              ! Interception coefficient b Von Hoyningen-Hune and Braden (L)

!     intercepted precipitation (rain+irrig) in mm
      rpd = grai*10.0d0
      if (isua.eq.0) rpd = (grai+gird)*10.0d0

!     exponential relation between soil cover and lai
      cofbb = 1.0d0 - dexp(-1.0d0*kdif*kdir*lai)
      cofbb = min(cofbb,1.0d0)

!     interception: evaporation of intercepted precipitation in cm
      if (cofab.gt.0.000001d0) then
        aintc = (cofab*lai*(1.0d0-(1/(1.0d0+rpd*cofbb/                  &
     &    (cofab*lai)))))*0.1d0
      else
        aintc = 0.0d0
      endif

! --- check data
!      if (aintc .gt. rpd) then
!         messag = ' Error-message from module '//trim(callroutine)//
!     &        ' Value for interception : ',aintc,
!     &        ' is larger than value for precipitation : ',rpd,
!     &        ' verify source code or input of meteo'
!         call fatalerr ('VonHHBraden',messag)
!      endif

      return
      end subroutine VonHHBraden


! SUBROUTINE 5.
! ----------------------------------------------------------------------
      subroutine Gash (aintc)
! ----------------------------------------------------------------------
!     Last modified      : July 2012
!     Purpose            : calculate interception for forests according to Gash (1995)
!     Interface:
!       I   - callroutine,logf,rai,gird,avevaptb,avprectb,
!             pfreetb,pstemtb,scanopytb,isua,t
!       O   - aintc
! ----------------------------------------------------------------------
      use variables, only: grai,gird,avevaptb,avprectb,pfreetb,pstemtb,scanopytb,isua,t
      implicit none
      include 'arrays.fi'

! --- global
      real(8)   aintc              ! Amount of rainfall interception during current day (L) in cm

! --- local
      real(8)   afgen
      real(8)   avevap             ! Gash interception model: average evaporation intensity during shower (-) as function of time (T)
      real(8)   avprec             ! Gash interception model: average rainfall intensity (-) as function of time (T)
      real(8)   cGash              ! slope of dPi/dPgross before saturation of canopy
      real(8)   pfree              ! Gash interception model: free throughfall coefficient (-) as function of time (T)
      real(8)   pstem              ! Gash interception model: stem flow coefficient (-) as function of time (T)
      real(8)   psatcan            ! Amount of rainfall to saturate canopy
      real(8)   rpd                ! intercepted precipitation (rain+irrig) in cm
      real(8)   scanopy            ! Gash interception model: storage capacity of canopy (-) as function of time (T)

!     intercepted precipitation (rain+irrig) in cm
      if (isua.eq.0) then
        rpd = grai+gird
      else
        rpd = grai
      endif

! --- calculate interception for forests according to Gash (1995)
      pfree = afgen(pfreetb,(2*magrs),t)
      pstem = afgen(pstemtb,(2*magrs),t)
      cGash = 1.d0-pfree-pstem
      scanopy = afgen(scanopytb,(2*magrs),t) / cGash
      avprec = afgen(avprectb,(2*magrs),t)
      avevap = afgen(avevaptb,(2*magrs),t) / cGash

! --- amount of rainfall to saturate canopy
      if ( (1.0d0 - avevap/avprec) .gt. 1.0d-4) then
        psatcan = -avprec*scanopy/avevap *                              &
     &              dlog(1.0d0 - avevap/avprec)
      else
        psatcan = avprec*scanopy/avevap
      endif

!     interception: evaporation of intercepted precipitation in cm
      if (grai .lt. psatcan) then
        aintc = cGash * rpd
      else
        aintc = cGash * ( psatcan +                                     &
     &            avevap*cGash / avprec * (rpd - psatcan) )
      endif

      return
      end subroutine Gash


! SUBROUTINE 6a.
! ----------------------------------------------------------------------
      subroutine ruttervw (gctp,aintc,eintc)
! ----------------------------------------------------------------------
!
! Program   ruttervw(logf,dt,cs,sicact,siccapact,fimin,ew0,grai,aintc,eintc)
! Author    Paul van Walsum
! Date      08/06/2012
! Purpose:  This routine simulates the interception process using the adapted
!           Rutter method of Van walsum & Supit (21012)
!
! Formal   parameters:  (I=input,O=output,C=control,IN=init,T=time)
! name     type meaning                                     units  class
! ----     ---- -------                                     -----  -----
! logf      I4  Internal number of logbook output file *.LOG   -      I
! dt        R8  time step                                      d      I
! gctp      R8  soil cover                                     -      I
! sicact    R8  storage on vegetaiton canopy                   cm    I/O
! siccapact R8  canopy storage capacity                        cm     I
! fimin     R8  startup fraction of reduction factor act/pot   -      I
! ew0       R8  potential transpiration of wet canaopy         cm/d   I
! grai      R8  gross rainfall                                 cm/d   I
! aintc     R8  intercepted rainfall                           cm/d   O
! eintc     R8  interception evaporation                       cm/d   O
!
!
      use variables, only: logf,dt,sicact,siccapact,fimin,ew0,grai
      implicit none
!
! ----------------------------------------------------------------------
!     globals
      real(8) gctp,aintc,eintc
!
!     locals
      integer(4) nuk_i4,ibd_i4(1),ib_i4
      real(4)   dc_r4,dtsw_r4,csk_r4(1),vxick_r4(1),fecmnk_r4(1)
      real(4)   ETw0_r4(1),Pgdtsw_r4(1),Sic_r4(1),Sicolddtsw_r4(1)
      real(4)   Picdtsw_r4(1),Eicdtsw_r4(1),tcap_r4(1),beta_r4(1)
      real(4)   zeta_r4(1),fricdtsw_r4(1)
!
!     The use of these arg's is to keep the msw1eic routine the same as in metaswap
      nuk_i4       = 1
      ibd_i4(1)    = 1
      ib_i4        = logf
      dc_r4        = 1.0e-4
      dtsw_r4      = REAL(dt)
      csk_r4(1)    = REAL(gctp)
      vxick_r4(1)  = REAL(siccapact)
      fecmnk_r4(1) = REAL(fimin)
      ETw0_r4(1)   = REAL(ew0*0.1d0)
      Pgdtsw_r4(1) = REAL(grai)
      Sic_r4(1)    = REAL(sicact)
!
      call msw1eic(nuk_i4,ibd_i4,dc_r4,dtsw_r4,csk_r4,vxick_r4,         &
     &  fecmnk_r4,ETw0_r4,Pgdtsw_r4,Sic_r4,Sicolddtsw_r4,Picdtsw_r4,    &
     &  Eicdtsw_r4,tcap_r4,beta_r4,zeta_r4,fricdtsw_r4,ib_i4)
!
      sicact = DBLE(Sic_r4(1))
      aintc  = DBLE(Picdtsw_r4(1))
      eintc  = DBLE(Eicdtsw_r4(1))
!
      return
      end subroutine ruttervw

! Subroutine 6b.
!**************************************************************************************
!** FILE:
!**    MSW1EIC.FOR
!**
!** COPYRIGHT: 2009
!**    Alterra
!**
!**    This PROGRAM, or parts thereof, may not be reproduced,
!**    modified or transferred to third parties without the
!**    written permission
!**
!**------------------------------------------------------------------------------------
!** PROGRAM:
!**    MSW1PENMON
!**
!** DESCRIPTION:
!**     Interception simulation with Sparse Gash concept, modified for relationship
!**     with saturation degree of canopy
!***************************************************************************************
!!
      SUBROUTINE  msw1eic(nuk,ibd,dc,dtsw,csk,vxick,fecmnk,ETw0,Pgdtsw, &
     &               Sic,Sicolddtsw,Picdtsw,Eicdtsw,tcap,beta,zeta,     &
     &               fricdtsw,ib)

!
!---  DECLARATIONS --------------------------------------------------------------------
!
      IMPLICIT     NONE
!
!     Global variables
!     formal parameters  : (i = input, o = output)
!     nuk           number of SVATs  ................................................i
!     ibd(k)        0/1 for existence of SVAT........................................i
!     dc            near-zero real ..................................................i
!     dtsw          time step [T]....................................................i
!     csk(k)        Soil cover [L2/L2]...............................................i
!     vxick(k)      Interception capacity of canopy [L3/L2]..........................i
!     fecmnk(k)     Minimum relative canopy evaporation factor [L3/L2]...............i
!     ETwo(k)       Evaporation from a wet canopy [L3/L2]............................i
!     Pgdtsw(k)     Gross rainfall + sprinkling [L]..................................i
!     Sic(k)        Interception storage of SVAT[L]........................ ........i/o
!     Sicolddtsw(k) Interception storage of SVAT, old value of time step [L]........i/o
!     Picdtsw(k)    Intercepted precipitation of SVAT [L]............................o
!     Eicdtsw(k)    Interception evaporation of SVAT [L].............................o
!     tcap(k)       Time to full interception reservoir [T] .........................o
!     beta(k)       Coefficient of differential equation [1/T].......................o
!     zeta(k)       Coefficient of differential equation [L/t].......................o
!     fricdtsw(k)   Fraction of time used by interception evaporation [-]............o
!     ib            Unit number of log file

!     Local variables
!     k             index for SVATs

      INTEGER(4)    nuk
      INTEGER(4)    ibd(1)
      real(4)       dc
      real(4)       dtsw
      real(4)       csk(1)
      real(4)       vxick(1)
      real(4)       fecmnk(1)
      real(4)       ETw0(1)
      real(4)       Pgdtsw(1)
      real(4)       Sic(1)
      real(4)       Sicolddtsw(1)
      real(4)       Picdtsw(1)
      real(4)       Eicdtsw(1)
      real(4)       tcap(1)
      real(4)       beta(1)
      real(4)       zeta(1)
      real(4)       fricdtsw(1)
      INTEGER(4)    ib
!
!-------------------------------------------------------------------------------------
!     Local variables
      INTEGER(4)    k
!-------------------------------------------------------------------------------------
      SAVE
!-------------------------------------------------------------------------------------
!
!
!$OMP PARALLEL DO
!$OMP&  DEFAULT(SHARED)
!$OMP&  PRIVATE(k)
      DO k=1,nuk
      IF (ibd(k) .GE. 1) THEN
!
!       Check that non zero interception capacity of canopy is accompanied by
!       non-zero soil cover:
        IF (vxick(k) .GT. dc) THEN
          IF (csk(k) .LT. dc) THEN
             WRITE(ib,9199) k
             WRITE(*,9199) k
             STOP
          ENDIF
        ENDIF
9199    FORMAT(' Interception capacity >0, but soil cover = 0, k =',i10)
!
!       Check first for shortcut to zero interception storage
        IF ( vxick(k) .LT. dc .OR.                                      &
     &        (Sic(k) .LT. dc .AND. Pgdtsw(k) .LT. dc) ) THEN
!
!         Interception capacity to zero, or nothing doing
          IF (Sic(k) .GT. dc) THEN
!
!           remaining interception water assumed to evaporate as vegetation dies off
            Eicdtsw(k) = Sic(k)/dtsw
            Sic(k)     = 0.
          ELSE
            Eicdtsw(k) = 0.
          ENDIF
          Picdtsw(k)   = 0.
        ELSE
!
!         Starting point of complete calculation
          Sicolddtsw(k) = Sic(k)
!
!         Check first for shortcut of full reservoir that stays full
          IF ( Sicolddtsw(k) .GE. (vxick(k) - 2*dc) .AND.               &
     &           (csk(k)*Pgdtsw(k)) .GE. ETw0(k) ) THEN
!
!           Reservoir starts full and stays full
            Sic(k)      = vxick(k)
            Eicdtsw(k)  = ETw0(k)  !    csk implicit in faeic
          ELSE
!
!           First check for beta=0 in differential equation
            IF (ETw0(k) .LT. dc .OR. fecmnk(k) .GT. 0.99999) THEN
!
!             Beta=0 (see below for expression)
              Sic(k) = Sicolddtsw(k) + (csk(k)*Pgdtsw(k) - ETw0(k))*dtsw
              Sic(k) = MIN(Sic(k),vxick(k))
              Sic(k) = MAX(Sic(k),0.)
              IF (Sic(k) .LT. vxick(k)) THEN
!
!               Reservoir does not become full, can become empty; evaporation is
!               obtained from balance because evaporation stops in the case of
!               the reservioir ending empty
                Eicdtsw(k) = (Sicolddtsw(k) - Sic(k))/dtsw +            &
     &                                           csk(k)*Pgdtsw(k)
              ELSE
!
!               Reservoir becomes full, evaporation during whole interval
!               at full rate (fecmn=1.) or zero (ETw0=0.)
                Eicdtsw(k) = ETw0(k)
              ENDIF
            ELSE
!
!             Auxiliary parameters for solving linear differential equation
              beta(k)  = (1.0 - fecmnk(k))*ETw0(k)/(vxick(k))
!
              zeta(k)  = csk(k)*Pgdtsw(k)- fecmnk(k)*ETw0(k)
!
!             First calculation of new Sic is tentative
              Sic(k)    = (Sicolddtsw(k) - zeta(k)/beta(k))*            &
     &                    EXP(-beta(k)*dtsw) + zeta(k)/beta(k)
!
!             Value can have become negative because formula did not take into
!             account that evaporation rate goes to zero when Sic goes to zero
              Sic(k)    = MAX(Sic(k),0.)
!
!             Final calculation of new Sic and Eictdtsw
              IF (Sic(k) .LT. vxick(k)) THEN
!
!               Reservoir does not become full, Ec from the balance
                Eicdtsw(k) = (Sicolddtsw(k)-Sic(k))/dtsw +              &
     &                                           csk(k)*Pgdtsw(k)
              ELSE
!
!               Reservoir becomes full; find out when this happens: tcap
                Sic(k)     = vxick(k)
                tcap(k)   = (1./beta(k))*                               &
     &                       LOG((Sicolddtsw(k)-zeta(k)/beta(k))/       &
     &                       (vxick(k) - zeta(k)/beta(k)))
!
!               From t=0 to tcap obtain Eic from a balance, rest of time Eic has potential rate
                Eicdtsw(k) = (1./dtsw)*( Sicolddtsw(k) - Sic(k) +       &
     &                                  csk(k)*Pgdtsw(k)*tcap(k) +      &
     &                                  ETw0(k)*(dtsw - tcap(k)) )
              ENDIF
            ENDIF
!
          ENDIF
!
!         Intercepted precipitation from balance
          Picdtsw(k) = (Sic(k)-Sicolddtsw(k))/dtsw + Eicdtsw(k) !
        ENDIF
!
!
!       Fraction of time that interception evaporation is active
        IF (ETw0(k) .GT. dc) THEN
          fricdtsw(k) = Eicdtsw(k)/ETw0(k)
        ELSE
          fricdtsw(k) = 0.
        ENDIF
!
      ENDIF
      ENDDO
!$OMP END PARALLEL DO

      RETURN
      END
!      end subroutine Ruttervw


! SUBROUTINE 7.
! ----------------------------------------------------------------------
      subroutine DivIntercep (aintc)
! ----------------------------------------------------------------------
!     Last modified      : februari 2014
!     Purpose            : divides interception into rain part and
!                          irrigation part and subsequently calculates
!                          net rain and net sprinkling irrigation
!     Interface:
!       I   - isua,aintc,gird,grai,gsnow,snrai
!       O   - nird,nraida
! ----------------------------------------------------------------------
      use variables, only: isua,gird,grai,gsnow,snrai,nird,nraida
      implicit none

! --- global
!   - In
      real(8)   aintc

!   - Out

! --- divide interception into rain part and irrigation part
!     and calculate net rain and net sprinkling irrigation
      if (aintc.lt.0.001d0) then
        nraida = grai - gsnow - snrai
        nird = gird
      else
         if (isua.eq.0) then
           nraida = grai-aintc*(grai/(grai+gird))
           nird = gird-aintc*(gird/(grai+gird))
         else
           nraida = grai-aintc
           nird = gird
         endif
      endif

      return
      end subroutine DivIntercep


! SUBROUTINE 8.
! ----------------------------------------------------------------------
      subroutine reduceva (task,nrai)
! ----------------------------------------------------------------------
!     date               : February 2014
!     purpose            : calculates reduction of soil evaporation
! ----------------------------------------------------------------------
      use variables, only: swredu,fldaystart,cofred,dt,empreva,    &
     &               ldwet,nird,peva,pond,rsigni,spev,saev
      implicit none

! --- global
      integer task
      real(8) nrai

! ----------------------------------------------------------------------
! --- local
      real(8) saevm1

! ----------------------------------------------------------------------
      select case (task)
      case (1)

! === reduction potential soil evaporation on daily basis ==============

      if (pond .gt. 1.0d-10) then
! ---   in case of ponding no reduction ---
        empreva = peva
        ldwet = 0.0d0
        spev = 0.0d0
        saev = 0.0d0
      elseif (swredu.eq.1) then
! ---   reduction with black model ---
        if ((nrai+nird) .gt. rsigni) ldwet = 0.0d0
        ldwet = ldwet + 1.0d0
        empreva = cofred * (dsqrt(ldwet)-dsqrt(ldwet-1.0d0))
        empreva = min(empreva,peva)
      elseif (swredu.eq.2) then
! ---   reduction with boesten and stroosnijder model ---
        if ((nrai+nird) .lt. peva) then
          spev = spev + (peva - (nrai+nird))
          saevm1 = saev
          if (spev .lt. cofred**2) then
            saev = spev
          else
            saev = cofred * dsqrt(spev)
          endif
          empreva = nrai + nird + saev - saevm1
        else
          empreva = peva
          saev = max (0.0d0,saev - (nrai+nird-peva))
          if (saev .lt. (cofred**2)) then
            spev = saev
          else
            spev = (saev/cofred)**2
          endif
        endif
      endif

      return

      case (2)

! === reduction potential soil evaporation for every time step ==============

      if (pond .gt. 1.0d-10) then
! ---   in case of ponding no reduction ---
        empreva = peva
        ldwet = 0.0d0
        spev = 0.0d0
        saev = 0.0d0
      elseif (swredu.eq.1) then
! ---   reduction with black model ---
        if (fldaystart .and. (nrai+nird) .gt. rsigni) ldwet = 0.0d0
        empreva = (cofred * (dsqrt(ldwet+dt)-dsqrt(ldwet))) / dt
        empreva = min(empreva,peva)
        ldwet = ldwet + dt
      elseif (swredu.eq.2) then
! ---   reduction with boesten and stroosnijder model ---
        if ((nrai+nird) .lt. peva) then
          spev = spev + (peva - (nrai+nird)) * dt
          saevm1 = saev
          if (spev .lt. cofred**2) then
            saev = spev
          else
            saev = cofred * dsqrt(spev)
          endif
          empreva = ((nrai + nird)*dt + saev - saevm1) / dt
        else
          empreva = peva
          saev = max (0.0d0,saev - (nrai+nird-peva)*dt)
          if (saev .lt. (cofred**2)) then
            spev = saev
          else
            spev = (saev/cofred)**2
          endif
        endif
      endif

      case default
         call fatalerr ('reduceva', 'Illegal value for TASK')
      end select

      return
      end subroutine Reduceva


subroutine CNmethod(Itask)
use variables, only: CNref, CNdry, CNwet, ThetaRef, theta, nraidt, Runoff_CN, zbotcp, dz, numnod, t1900, wc_cor, iCNtab, CNtimTAB, CNrefTAB, melt, wc10
implicit none
! global
integer, intent(in)  :: Itask
! local
integer              :: i
integer, save        :: Nod10, iCN   !, t1900_old
real(8)              :: wc1, wc2, CN, S, Ia, watcon
real(8), save        :: Z10

select case (Itask)
! initialization; calculate and store some constants
case (1)

   iCN = 0
   ! check if times in CNtimeTAB are in ascending order
   ! set initial position in CNtimTAB
   do i = 2, iCNtab
      if (CNtimTAB(i) < CNtimTAB(i-1)) call fatalerr ('CNmethod', 'CNtimTAB not in ascending order')
      if (t1900 >= CNtimTAB(i-1) .and. t1900 < CNtimTAB(i)) iCN = i-1
   end do
   ! error if start time t1900 not in CNtimTAB
   if (iCN == 0) call fatalerr ('CNmethod', 'Start time of simulation not present in CNtimTAB')
   
!  to be replaced by average for layer 0-10 cm
   do i = 1, numnod
      if (zbotcp(i) < -10.0d0) then
         Nod10 = i-1
         Z10 = -zbotcp(Nod10)
         exit
      end if
   end do
   ThetaRef = 0.0d0
   do i = 1, Nod10
      if (wc_cor == 1) then
         wc1 = watcon(i,-100.0d0)
         wc2 = watcon(i,-16000.0d0)
         ThetaRef = ThetaRef + (wc1+wc2)*0.5d0*dz(i)
      else if (wc_cor == 2) then
         wc1 = watcon(i,0.0d0)
         wc2 = watcon(i,-16000.0d0)
         ThetaRef = ThetaRef + (wc1+wc2)*0.5d0*dz(i)
      end if
   end do
   ThetaRef = ThetaRef/Z10
   
   !!!t1900_old = int(t1900) - 1 ! for testing intermediate output
   
! dynamic part: calculate runoff   
case (2)
   
   ! see if t1900 has moved ahead in CNtimTAB; iCN can nver exceed last entry
   if (iCN < iCNtab .and. t1900 >= CNtimTAB(iCN+1)) iCN = iCN + 1
   CNref = CNrefTAB(iCN)
   CN    = CNref
   CNdry =  4.2d0*CNref/(10.0d0-0.058d0*CNref)
   CNwet = 23.0d0*CNref/(10.0d0+0.13d0*CNref)

   if (wc_cor > 0) then
      wc10 = 0.0d0
      do i = 1, Nod10
         wc10 = wc10 + theta(i)*dz(i)
      end do
      wc10 = wc10/Z10
      if (wc10 < ThetaRef) then
         CN = CNdry + wc10/ThetaRef*(CNref-CNdry)
      else
         CN = CNref + (wc10-ThetaRef)/ThetaRef*(CNwet-CNref)
      end if
   end if
   S = 2540d0/CN-25.4d0    ! in cm
   Ia = 0.2d0*S
!   Ia = 0.3d0*S
   if (nraidt+melt > Ia) then
      Runoff_CN = (nraidt+melt-Ia)**2/(nraidt+melt-Ia+S)
   else
      Runoff_CN = 0.0d0
   end if
   
   ! for testing intermediate output
   !!!if (int(t1900) > t1900_old) then
   !!!   write (123, '(7F20.6)') t1900, nraidt, runoff_cn, cn, wc10, thetaref, melt
   !!!   t1900_old = int(t1900)
   !!!end if
   
case default
   call fatalerr ('CNmethod', 'Illegal Itask option')
end select

end subroutine CNmethod
