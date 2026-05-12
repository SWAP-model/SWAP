! File VersionID:
!   $Id: timecontrol.f90 374 2018-03-21 13:12:23Z heine003 $
! ----------------------------------------------------------------------
      subroutine TimeControl(task, state)
      use error_mod, only: fatalerr_collected
! ----------------------------------------------------------------------
!     Date               : Aug 2004
!     Purpose            : Handles time variables, switches and flags
! ----------------------------------------------------------------------

      use variables
      use timestep_control_mod, only: fldecdt
      use irrigation_mod, only: SSDI_irrigation
      use swap_state_mod, only: swap_state_t            ! [SS-SWC S-2.12B]
      implicit none

      ! [SS-SWC S-2.12B] state added — needed to pass to SSDI_irrigation(9)
      type(swap_state_t), intent(inout) :: state
! ----------------------------------------------------------------------
!     DAYNR  : = daynumber relative to start of calendar year
!     DAYCUM : = daynumber relative to start of simulation
!     DAYCROP: = daynumber relative to emergence of crop
!     DAYMETEO: = day number for which meteorological data should be read
!     IMONTH : = actual month number
!     IYEAR  : = actual year
!
!     T      : = time relative to start of calendar year
!     TCUM   : = time relative to start of simulation
!     TEND   : = time (relative to 1900) at which simulation ends
!     TSTART : = time (relative to 1900) at which simulation starts
!     t1900  : = time relative to 1900
!
!     flCultivate : = true if crops are grown on the soil profile
!     flDayEnd : = true if time level equals end of day
!     flRunEnd : = true if time level equals TEND
!     flRunStart  : = true if time level equals TSTART
! ----------------------------------------------------------------------
!     local
      character(len=200) messag
      character(len=11)  tmp
      character(len=80)  filtext
      integer       task, itask
      real(8)       :: dtCrit, dtfletsine
      real(8)       :: dtRestDay


!     critical time-interval 
      dtCrit = 1.d-8
      dtfletsine = 0.05d0


      ! Bind former SAVE locals to explicit, synchronizable module variables
      associate( datea => tc_datea, &
           nextyear => tc_nextyear, &
           flprevious => tc_flprevious, &
           flTnext => tc_flTnext, &
           fsec => tc_fsec, &
           tchange => tc_tchange, &
           dtEvent => tc_dtEvent, &
           tEvent => tc_tEvent, &
           tcumold => tc_tcumold, &
           dtprevious => tc_dtprevious, &
           tmptimestart => tc_tmptimestart, &
           tmptimeend => tc_tmptimeend )

      itask = task
      if (itask.eq.2 .and. (fldecdt .or. fldecdtmin)) itask = 3

      select case (task)
      case (1)

! === initialization ===================================================

! --- initialize flags ----------------------------
      fldecdt = .false.
      fldecdtmin = .false.
      state%timecontrol%fldecdtmin = fldecdtmin
      flRunEnd = .false.
      state%timecontrol%flRunEnd = flRunEnd
      flDayStart = .true.
      state%timecontrol%flDayStart = flDayStart
      fldtmin = .false.
      state%timecontrol%fldtmin = fldtmin
      flZeroIntr = .true.
      flZeroCumu = .true.
      floutput = .false.
      state%timecontrol%floutput = floutput
      flbaloutput = .false.
      state%timecontrol%flbaloutput = flbaloutput
      flheader = .false.
      state%timecontrol%flheader = flheader
      flheadirg = .false.
      state%timecontrol%flheadirg = flheadirg
      flIrg1Start = .true.
      state%timecontrol%flIrg1Start = flIrg1Start
      flUpdMetDet = .true.
      state%timecontrol%flUpdMetDet = flUpdMetDet
      flYearStart = .true.
      state%timecontrol%flYearStart = flYearStart

      if (nprintday .gt. 1 .or. flprintdt) then
        flprintshort = .true.
        period = 1
      else
        flprintshort = .false.
      endif
      state%timecontrol%flprintshort = flprintshort
      if (swmetdetail .eq. 1) then
        flmetdetail = .true.
        metperiod = 1.0d0 / dble(nmetdetail)
        state%timecontrol%metperiod = metperiod
        wrecord = 0
        state%timecontrol%wrecord = wrecord
      else
        flmetdetail = .false.
      endif
      state%timecontrol%flmetdetail = flmetdetail
      if (swrain.gt.0) then
         flrainintens = .true.
      else
         flrainintens = .false.
      endif
      state%timecontrol%flrainintens = flrainintens
      if (flmetdetail .or. flrainintens) then
         flmeteodt = .true.
      else
         flmeteodt = .false.
      endif
      state%timecontrol%flmeteodt = flmeteodt
      fletsine = .false.
      if (swetsine .eq. 1) fletsine = .true.
      state%timecontrol%fletsine = fletsine
      if (swirfix.eq.1) flIrrigate = .true.
      state%timecontrol%flIrrigate = flIrrigate
      flDrain = .false.
      if (swdra .eq. 1) flDrain = .true.
      state%timecontrol%flDrain = flDrain
      flSurfaceWater = .false.
      if (swdra .eq. 2) flSurfaceWater = .true.
      state%timecontrol%flSurfaceWater = flSurfaceWater
      flTemperature = .false.
      if (swhea .eq. 1) flTemperature = .true.
      state%timecontrol%flTemperature = flTemperature
      flSnow = .false.
      if (swsnow .eq. 1) flSnow = .true.
      state%timecontrol%flSnow = flSnow
      flSolute = .false.
      if (swsolu .eq. 1) flSolute = .true.
      state%timecontrol%flSolute = flSolute

! --- initialize counters ----------------------------
      nirri = 1
      isteps = 0
      state%timecontrol%isteps = isteps
      ioutdat = 1
      state%timecontrol%ioutdat = ioutdat
      ioutdatint = 1
      state%timecontrol%ioutdatint = ioutdatint
      cntper = 0
      state%timecontrol%cntper = cntper
      outper = 0.0d0
      state%timecontrol%outper = outper
      tcumold = 0.0d0
      state%timecontrol%tcumold = tcumold
      nprintcount = 1
      state%timecontrol%nprintcount = nprintcount

! --- set main time variable of SWAP ------------------------
      t1900 = tstart
      state%timecontrol%t1900 = t1900

! --- determine time from beginning of calendar year
      call dtdpar (t1900, datea, fsec)
      state%timecontrol%datea = datea
      state%timecontrol%fsec = fsec
      datea(1) = iyear
      datea(2) = 1
      datea(3) = 1
      state%timecontrol%datea = datea
      fsec = 0.0
      state%timecontrol%fsec = fsec
      call dtardp (datea,fsec,timjan1)
      state%timecontrol%timjan1 = timjan1
      t = tstart - timjan1
      state%timecontrol%t = t
      tEvent = 0.0d0
      state%timecontrol%tEvent = tEvent
      tcum = 0.d0
      state%timecontrol%tcum = tcum
      daynr = nint(t)
      state%timecontrol%daynr = daynr
      daycum = 0
      state%timecontrol%daycum = daycum
      daymeteo = daynr + 1
      state%timecontrol%daymeteo = daymeteo

! --- determine year,month and day number of current day
      call dtdpar (t1900+0.1d0, datea, fsec)
      state%timecontrol%datea = datea
      state%timecontrol%fsec = fsec
      iyear = datea(1)
      state%timecontrol%iyear = iyear
      imonth = datea(2)
      state%timecontrol%imonth = imonth
      yearmeteo = iyear
      state%timecontrol%yearmeteo = yearmeteo

! --- determine date of current day
      call dtdpst ('year-month-day',t1900+0.1d0,date)
      state%timecontrol%date = date

! --- output to screen
      if (swscre .eq. 2) then
        filtext = 'Screen output of daynumbers'
        call writehead (5,1,'screen',filtext,project)
        call dtdpst ('year-month-day',tstart,date)
        state%timecontrol%date = date
        write (*,'(2x,2a)') 'First day of simulation:  ',date
        call dtdpst ('year-month-day',tend,date)
        state%timecontrol%date = date
        write (*,'(2x,2a)') 'Last day of simulation:   ',date
        write (*,'(/,a,/)') '            date  daynr  daycum'
      endif

! -   set crop number
      icrop = 1
      do while (.not. flCropCalendar) 
          
        if (cropstart(icrop) .lt. 1.d0) exit
          
        if (abs(tstart - cropstart(icrop)) .lt. 1.d-3) then
          flCropCalendar = .true.
          if (flCropCalendar) then
            if (tstart - cropstart(icrop) .lt. -1.d-3 .and.             &
     &                                              swinco .ne. 3) then
              messag = 'The start of simulation (tstart) begins in '//  &
     &        'crop growing season with swinco 1 or 2'
              call fatalerr_collected ('readswap',messag)
            endif
          end if
        else
          icrop = icrop + 1
        end if
      enddo
      
! --- detailed meteo data needed for crop growth?
      swmeteo = 1
      if (flCropCalendar) then
        if (icrop .gt. 0) then
          if(croptype(icrop).ge.2) then
             swmeteo = 2
          endif
        endif
      endif
      state%timecontrol%swmeteo = swmeteo

! --- initialize
      dtEvent = 1.0d0

! --- limit initial dt in case of short time interval
      if (flprintshort .and. .not.flprintdt) then
        if (dtEvent .gt. dble(nprintcount)/dble(nprintday))then
           dtEvent = dble(nprintcount)/dble(nprintday)
         end if
      endif

! --- intial timestep for rainfall intensities
      if (swmetdetail.eq.0 .and. swrain.gt.0) then
         dtEvent = min(dtEvent, dtmin)
      endif

! --- limit initial dt in case of detailed meteorological input
       if (flmetdetail) then
        if (dtEvent .gt. metperiod)then
           dtEvent = min(dtEvent, metperiod)
         end if
      endif
      state%timecontrol%dtEvent = dtEvent
      tEvent = dtEvent
      state%timecontrol%tEvent = tEvent

! --- initial time step
      if (swinco.eq.3) then
        if (dt.lt.dtmin) then
          messag = 'Initial dt read from file (SWINCO=3)'//             &
     &    ' if absent, then default is assumed'
          call warn ('soilwater',messag,logf,swscre)
          dt = dsqrt(dtmin*dtmax)
        endif
      else
        dt = dsqrt(dtmin*dtmax)
      endif
      state%timecontrol%dt = dt
      dtprevious = dt
      state%timecontrol%dtprevious = dtprevious
      if(dt+dtCrit .gt. dtEvent)then
         flprevious = 2
         dt = dtEvent
         state%timecontrol%dt = dt
         flTnext = .true.
      else
         flprevious = 1
         flTnext = .false.
      end if
      state%timecontrol%flprevious = flprevious
      state%timecontrol%flTnext = flTnext

! --- in case of sine wave of ET, limit dt and dtmax
      if (fletsine) then
        dt = min(dt, dtfletsine)
        dtmax = min(dtmax, dtfletsine)
        state%timecontrol%dt = dt
      endif

! --- set value for dtold, used in headcalc for macropore-iteration
      dtold = dt
      state%timecontrol%dtold = dtold

      return

      case (2)

! === next time step ===================================================

! 2.1  check maximum number of time steps during this day
      isteps = isteps + 1
      state%timecontrol%isteps = isteps
      if (isteps .gt. msteps) then
        write(tmp,'(i11)') daynr
        tmp = adjustl(tmp)
        messag ='The maximum number of time steps for a day is exceeded'&
     &    //' at daynumber '//trim(tmp)//'. Check input for numerical'  &
     &    //' solution of Richards equation'
        call fatalerr_collected ('timer',messag)
      endif


! 2.2 update time variables

      t = t + dt             ! relative to yyyy0101:00:00:00
      state%timecontrol%t = t
      tcum = tcum + dt
      state%timecontrol%tcum = tcum
      t1900 = tstart + tcum
      state%timecontrol%t1900 = t1900

! 2.3  flag assignments after first time step of a day

      if (flDayStart) then

! ---   set crop conditions
        if (flCropCalendar) then
          flCropOutput = .true.
          if (flCropHarvest) then
            flCropOutput = .false.
          endif
        endif

!! --    length crop
!        if (flCropCalendar) then
!          daycrop = daycrop + 1
!        endif

! ---   set flags for reset intermediate and cumulative fluxes
        flZeroCumu = .false.
        if (flZeroIntr) then
          outper = 0.0d0
          state%timecontrol%outper = outper
          flZeroIntr = .false.
        endif

! ---   set flags for output
        if (floutput) then
         floutput = .false.
         state%timecontrol%floutput = floutput
         if (swheader .eq. 1) then
           flheader = .false.
           state%timecontrol%flheader = flheader
         endif
        endif
        if (flbaloutput) then
         flbaloutput = .false.
         state%timecontrol%flbaloutput = flbaloutput
         if (swheader .eq. 1) then
           flheader = .true.
           state%timecontrol%flheader = flheader
           flheadirg = .true.
           state%timecontrol%flheadirg = flheadirg
         endif
        endif

! 2.4  determine year,month and day number (only during first time step of a day)
        call dtdpar (t1900, datea, fsec)
        state%timecontrol%datea = datea
        state%timecontrol%fsec = fsec
        iyearm1 = iyear
        state%timecontrol%iyearm1 = iyearm1
        iyear = datea(1)
        state%timecontrol%iyear = iyear
        imonth = datea(2)
        state%timecontrol%imonth = imonth

! ---   determine date of current day
        call dtdpst ('year-month-day',t1900,date)
        state%timecontrol%date = date

! ---   update day numbers
        daynr = daynr+1
        state%timecontrol%daynr = daynr
        daycum = daycum + 1
        state%timecontrol%daycum = daycum
        cntper = cntper + 1
        state%timecontrol%cntper = cntper

! 2.5  in case of detailed meteorological input, reset weather record

        if (iyear .ne. iyearm1) then
! ---     reset daynumber and time because new calender year has started
          daynr = 1
          state%timecontrol%daynr = daynr
          t = dt
          state%timecontrol%t = t

! ---     in case SWRES = 1 reset counter for periodic output to 1
          if (swres.eq.1 .and. period.ne.0) then
            cntper = 1
            state%timecontrol%cntper = cntper
          endif
        endif


      endif


! 2.6  update  logicals for indication of start / end of day
      if ( dble(daycum) - tcum .lt. dtCrit) then
        flDayEnd = .true.
        state%timecontrol%flDayEnd = flDayEnd
        flDayStart = .true.
        state%timecontrol%flDayStart = flDayStart
        flTnext = .true.
        state%timecontrol%flTnext = flTnext
        if (flmetdetail) then
          wrecord = 1
          state%timecontrol%wrecord = wrecord
        endif
      else
        flDayEnd = .false.
        state%timecontrol%flDayEnd = flDayEnd
        flDayStart = .false.
        state%timecontrol%flDayStart = flDayStart
      endif


! 2.7  maximum size of time interval (dtEvent)

      if(flTnext)then

! 2.7.1  remaining part of a day
        dtEvent = dble(int(tcum+1.0d0+dtCrit)) - tcum

! 2.7.2  printing more than one times a day may limit timestep
         if (flprintshort  .and. .not.flprintdt) then
           if (dble(nprintcount)/dble(nprintday)-tcum .lt. dtCrit) then
              dtEvent = min(dtEvent,1.0d0/dble(nprintday))
           else
              dtEvent = min(dtEvent,                                    &
     &                   (dble(nprintcount)/dble(nprintday)) - tcum)
           endif
         endif

! 2.7.3  input of detailed meteo may limit timestep
         if (flmetdetail) then
            tchange = dble(int(t + dtCrit)) + dble(wrecord) * metperiod
            state%timecontrol%tchange = tchange
            if ((t+dtEvent) .gt. tchange) dtEvent = tchange - t
         end if

! 2.7.4  precipitation event may limit timestep
         if (swmetdetail.eq.0 .and. swrain.gt.0) then

!        next rainevent! Set new values
           if (raintimearray(rainrec).lt.tcum+dtCrit) then
!        new rain event valid
             rainrec = rainrec + 1
             state%timecontrol%rainrec = rainrec
           endif

           dtEvent = min(dtevent,dtEventRain)
           dtEvent = max(dtEvent,dtmin)
         endif
         state%timecontrol%dtEvent = dtEvent

! 2.8  set end of time interval (determined by I/O)
         tEvent = tEvent + dtEvent
         state%timecontrol%tEvent = tEvent
         flTnext = .false.
         state%timecontrol%flTnext = flTnext

      endif

! 2.9 determine next time step, based on numerical performance
      if (flDayStart) then
        if(flprevious.eq.2)then
           dt = max(dt, dsqrt(dtmin*dtmax),dtprevious)
        else
           dt = max(dt, dsqrt(dtmin*dtmax))
        end if
        dtprevious = dt
        state%timecontrol%dtprevious = dtprevious
      else
        if(flprevious .eq. 2)then
           dt = dtprevious
        else
           if (numbit.le.3)     dt = min(dt*2.0d0,DtMax)
           if (numbit.ge.MaxIt) dt = max(dt*0.5d0,DtMin)
           dtprevious = dt
           state%timecontrol%dtprevious = dtprevious
         endif
      endif
      state%timecontrol%dt = dt
      flprevious = 1
      state%timecontrol%flprevious = flprevious

      if ( tcum + dt - tEvent .gt. dtCrit) then
         dt = tEvent - tcum
         state%timecontrol%dt = dt
         flprevious = 2
         state%timecontrol%flprevious = flprevious
         flTnext = .true.
         state%timecontrol%flTnext = flTnext
      endif

!     SSDI: adapt dt as to not pass dt_SSDI_event end time of the day
      if (dt_SSDI_event < 1.0d0) then
         dt = max(dtmin, min(dt, dble(int(tcum) + dt_SSDI_event) - tcum))
         state%timecontrol%dt = dt
      end if

! 2.10  test last time step of the day: limit dt if it exceeds end of day
      dtRestDay = dble(int(tcum+1.0d0+dtCrit)) - tcum
      if (dtRestDay-dt.lt.1.d-6) then
        dt = dtRestDay
      endif

! JK20131230: when dtmin is large then the timestep-closure may not be correct and errors may occur
! resulting in water balance errors see email Paul v Walsum 210131225. Elimination of next statement
! is not the correct solution
      dt = max(dt,dtmin)
      dt = min(dt,dtmax)
      state%timecontrol%dt = dt

! 2.11 set flags and variables

! --- in case of output during a day
      if (flprintshort) then
        floutputshort = .false.
        state%timecontrol%floutputshort = floutputshort
        flzerointr = .false.
! ---   determine whether output is required
        if (flprintdt) then
           outper = tcum - tcumold
           state%timecontrol%outper = outper
           tcumold = tcum
           state%timecontrol%tcumold = tcumold
           if (abs(outdatint(ioutdatint) - t1900 + 1.d0).lt.1.d-3) then
              floutputshort = .true.
              state%timecontrol%floutputshort = floutputshort
              flzerointr = .true.
           endif
        else
           if (tcum+dtCrit .gt. dble(nprintcount)/dble(nprintday)) then
              floutputshort = .true.
              state%timecontrol%floutputshort = floutputshort
              flzerointr = .true.
              outper = tcum - tcumold
              state%timecontrol%outper = outper
              tcumold = tcum
              state%timecontrol%tcumold = tcumold
           endif
! ---      update counter nprintcount for printing
           do while (tcum+dtCrit .gt. dble(nprintcount)/dble(nprintday))
              nprintcount = nprintcount + 1
           end do
           state%timecontrol%nprintcount = nprintcount
         endif
      endif

! --- in case of detailed meteorological input
      if (flmetdetail) then
        tchange = dble(int(t + dtCrit)) + dble(wrecord) * metperiod
        state%timecontrol%tchange = tchange
        if ((tchange - t) .lt. dtCrit) then
! ---     update actual weather record and fluxes

          flUpdMetDet = .true.
          state%timecontrol%flUpdMetDet = flUpdMetDet
        endif
      endif

!     SSDI: end of subsurface irirgation event reached; reset
      if (flSSDI .and. tcum - int(tcum) + dtCrit > dt_SSDI_event) then
         call SSDI_irrigation(9, state)  ! [SS-SWC S-2.12B]
      end if

! --- update fldtmin
      if (dt .gt. (1.0d0+dtCrit)*dtmin) then
        fldtmin = .false.
        state%timecontrol%fldtmin = fldtmin
      endif


! --- procedure when day is finished ----------------------------------------------
      if (flDayEnd) then

!! --    length crop
!        if (flCropCalendar) then
!          daycrop = daycrop + 1
!        endif

! ---   length output period
        if (.not.flprintshort) then
          outper = outper + 1.0d0
          state%timecontrol%outper = outper
        endif

! ---   write daynumber to screen
        if (swscre .eq. 2) then
          write(*,'("+ ",4x,a11,i6,i8)') date,daynr,daycum
        endif

! ---   end of run?
        if ((tend - t1900 + 1.d0) .lt. 1.d-3) then
          flRunEnd = .true.
          state%timecontrol%flRunEnd = flRunEnd
          floutput = .true.
          state%timecontrol%floutput = floutput
          flbaloutput = .true.
          state%timecontrol%flbaloutput = flbaloutput
          ioutdat = ioutdat + 1
          state%timecontrol%ioutdat = ioutdat
          return
        endif

! ---   in case no end of run, determine whether today output should be written
        if (cntper .eq. period .and. .not.flprintdt) then
          cntper = 0
          state%timecontrol%cntper = cntper
          floutput = .true.
          state%timecontrol%floutput = floutput
          flzerointr = .true.
        endif

        if(flprintdt) then
          if (abs(outdatint(ioutdatint) - t1900 + dt) .lt. 1.d-3) then
            floutput = .true.
            state%timecontrol%floutput = floutput
            flzerointr = .true.
            ioutdatint = ioutdatint + 1
            state%timecontrol%ioutdatint = ioutdatint
          endif
        else
          if (abs(outdatint(ioutdatint) - t1900 + 1.d0) .lt. 1.d-3) then
            floutput = .true.
            state%timecontrol%floutput = floutput
            flzerointr = .true.
            ioutdatint = ioutdatint + 1
            state%timecontrol%ioutdatint = ioutdatint
          endif
        endif
        if (abs(outdat(ioutdat) - t1900 + 1.d0) .lt. 1.d-3) then
! ---     output of water and solute balances
          floutput = .true.
          state%timecontrol%floutput = floutput
          flbaloutput = .true.
          state%timecontrol%flbaloutput = flbaloutput
          flzerointr = .true.
          flzerocumu = .true.
          ioutdat = ioutdat + 1
          state%timecontrol%ioutdat = ioutdat
        endif


! ---   reset flags for next day
        fldtmin = .false.
        state%timecontrol%fldtmin = fldtmin

! ---   reset counters for next day
        isteps = 0
        state%timecontrol%isteps = isteps

! ---   determine daynumber and switch for reading meteorological data
        call dtdpar (t1900 + 0.1d0, datea, fsec)
        state%timecontrol%datea = datea
        state%timecontrol%fsec = fsec
!        call dtdpar (t1900 , datea, fsec)
        nextyear = datea(1)
        state%timecontrol%nextyear = nextyear
        if (nextyear .eq. iyear) then
          daymeteo = daynr + 1
          state%timecontrol%daymeteo = daymeteo
        else
          yearmeteo = nextyear
          state%timecontrol%yearmeteo = yearmeteo
! ---     set flag for new meteo year
          flYearStart = .true.
          state%timecontrol%flYearStart = flYearStart

          daymeteo = 1
          state%timecontrol%daymeteo = daymeteo
! ---     detailed meteo data needed for crop growth?
          swmeteo = 1
          if (flCropCalendar) then
            if (icrop .gt. 0) then
              if (croptype(icrop).ge.2) then
                swmeteo = 2
              endif
            endif
          endif
          state%timecontrol%swmeteo = swmeteo
        endif

      endif

      dtprevious = dt
      state%timecontrol%dtprevious = dtprevious

      return

      case (3)

! === reduce time step ===================================================


! --- decrease time step in case of no convergence in headcalc
      if (fldecdt) then
        if (dt .gt. 3.0*dtmin) then
          dt = dt / 3.0
!         force dt to equal multiple dtmin to prevent very small dt-values at end of day
!          dt = dtmin * dble(max(1,int(dt/dtmin)))
        else
          dt = dtmin
          fldtmin = .true.
          state%timecontrol%fldtmin = fldtmin
        endif
        state%timecontrol%dt = dt
        fldecdt = .false.
        flprevious = 1
        state%timecontrol%flprevious = flprevious
        dtprevious = dt
        state%timecontrol%dtprevious = dtprevious
        flTnext = .false.
        state%timecontrol%flTnext = flTnext

        return
      endif

! --- decrease time step to dtmin if required by boundtop
      if (fldecdtmin) then
        dt = dtmin
        state%timecontrol%dt = dt
        fldtmin = .true.
        state%timecontrol%fldtmin = fldtmin
        fldecdtmin = .false.
        state%timecontrol%fldecdtmin = fldecdtmin
        flprevious = 1
        state%timecontrol%flprevious = flprevious
        dtprevious = dt
        state%timecontrol%dtprevious = dtprevious
        return
      endif

! --- decrease in case of Macropores
      if (flMacroPore .and. FlDecMpRat) then
        dt = dsqrt(dtmin*dtmax)
        state%timecontrol%dt = dt
        dtprevious = dt
        state%timecontrol%dtprevious = dtprevious
        return
      endif

      case (9)
!        special: at end of day, the possible initial time step for next day may be too large: adapt if necessary
         if (dt_SSDI_event < 1.0d0) then
            dt = min(dt, dt_SSDI_event)
            state%timecontrol%dt = dt
         end if
      
      case default
         call fatalerr_collected ('TimeControl', 'Illegal value for TASK')
      end select

      end associate

      return
      end 


      subroutine IterTime(task, state)
      use error_mod, only: fatalerr_collected
! ----------------------------------------------------------------------
!     date    : 20080303
!     update  : 20170223: intro part2 - to interrupt (near) endless simulations
!     purpose : statistics of timing and numerical iterations
! ----------------------------------------------------------------------
! --- global variables
      use variables
      use swap_state_mod, only: swap_state_t            ! [SS-TC TC-2]
      implicit none

! --- local variables
      integer task, i, j, timediff
      type(swap_state_t), intent(inout) :: state
      character(len=400) messag
      real(4)       ::   tmptimeinterrupt
      ! Use module variables (tc_tmptimestart/tc_tmptimeend) for persistence
      ! to keep IterTime multi-instance safe.

      select case (task)

      case (1)
! --- part1 - initial values
      call cpu_time(tc_tmptimestart)
      state%timecontrol%tmptimestart = tc_tmptimestart
      return

      case (2)
! --- part2 - calculate intermediate time to be able to interrupt (near) endless simulations
      call cpu_time(tmptimeinterrupt)
      timediff = int(tmptimeinterrupt)-MaxIterTime
!     fatal error if cpu time exceeds input value MaxIterTime
      if(timediff.gt.0) then
        write(messag,'(a,i10,3a)')                                      &
     &     'The maximum cpu time of ',MaxIterTime,' (secs)',            &
     &     ' was exceeded.  Therefore simulation was interrupted'
        call fatalerr_collected ('IterTime',messag)
      endif
      return

      case (3)
! --- part3 - write statistics
      write(logf,'(/,a20)')      'Iteration statistics'
      write(logf,'(/,a29,i4)')   'Maximum number of iterations:',MaxIt
      write(logf,'(/,a35/,a35)') 'It Numb  No of Hits  Tot BTr cycles', &
     &                           '-------  ----------  --------------'
      do i=1,100                                       
        if(itnumb(i,1).gt.0)                                            &
     &     write(logf,'(i7,2x,i10,4x,i10)')i,(itnumb(i,j),j=1,2)
      end do

      call cpu_time(tc_tmptimeend)
      state%timecontrol%tmptimeend = tc_tmptimeend
      write(logf,'(/,a12,f12.2,a4)')                                    &
         &           ' Run-time: ',tc_tmptimeend-tc_tmptimestart,' sec'
     

      case default
         call fatalerr_collected ('IterTime', 'Illegal value for TASK')
      end select

      return
      end
