!> @file timecontrol_mod.f90
!! SS-TCM: module form of the legacy subroutine TimeControl + IterTime.
!! Seven named lifecycle procedures replace the magic-int dispatch
!! (task=1/2/3/9 + IterTime task=1/2/3). State threaded explicitly;
!! each procedure carries its own associate block. flZeroIntr and
!! flZeroCumu are owned by state%timecontrol — bare globals
!! consumed by Task 11 readers retire in Task 12.
module timecontrol_mod
   use swap_state_mod, only: swap_state_t
   implicit none
   private
   public :: timecontrol_init, timecontrol_advance, &
             timecontrol_reduce_dt, timecontrol_day_end
   public :: itertime_init, itertime_check, itertime_close

contains

   subroutine timecontrol_init(state)
      use variables, only: dtmin, dtmax, period, swscre, &
                            tend, tstart, &
                            flprintdt, nprintday, logf, flCropCalendar, &
                            swirfix, swsnow, swdra, &
                            swhea, swsolu, swetsine, swrain, swmetdetail, &
                            flzerointr, flzerocumu, &
                            nmetdetail, nirri, swinco, icrop, &
                            cropstart, croptype, project
      use timestep_control_mod, only: fldecdt
      use error_mod, only: fatalerr_collected
      implicit none
      type(swap_state_t), intent(inout) :: state

      ! local variables
      character(len=200) :: messag
      character(len=80)  :: filtext
      real(8)            :: dtCrit, dtfletsine

      dtCrit = 1.d-8
      dtfletsine = 0.05d0

      ! Associate block — copied verbatim from timecontrol.f90:54-114.
      ! Each TC field aliased; bare names below resolve to state%timecontrol%*.
      associate( datea => state%timecontrol%datea, &
           nextyear => state%timecontrol%nextyear, &
           flprevious => state%timecontrol%flprevious, &
           flTnext => state%timecontrol%flTnext, &
           fsec => state%timecontrol%fsec, &
           tchange => state%timecontrol%tchange, &
           dtEvent => state%timecontrol%dtEvent, &
           tEvent => state%timecontrol%tEvent, &
           tcumold => state%timecontrol%tcumold, &
           dtprevious => state%timecontrol%dtprevious, &
           tmptimestart => state%timecontrol%tmptimestart, &
           tmptimeend => state%timecontrol%tmptimeend, &
           iyear => state%timecontrol%iyear, &
           iyearm1 => state%timecontrol%iyearm1, &
           imonth => state%timecontrol%imonth, &
           dt => state%timecontrol%dt, &
           dtold => state%timecontrol%dtold, &
           daynr => state%timecontrol%daynr, &
           daycum => state%timecontrol%daycum, &
           daymeteo => state%timecontrol%daymeteo, &
           yearmeteo => state%timecontrol%yearmeteo, &
           t => state%timecontrol%t, &
           t1900 => state%timecontrol%t1900, &
           tcum => state%timecontrol%tcum, &
           timjan1 => state%timecontrol%timjan1, &
           outper => state%timecontrol%outper, &
           cntper => state%timecontrol%cntper, &
           isteps => state%timecontrol%isteps, &
           ioutdat => state%timecontrol%ioutdat, &
           ioutdatint => state%timecontrol%ioutdatint, &
           nprintcount => state%timecontrol%nprintcount, &
           rainrec => state%timecontrol%rainrec, &
           wrecord => state%timecontrol%wrecord, &
           swmeteo => state%timecontrol%swmeteo, &
           date => state%timecontrol%date, &
           metperiod => state%timecontrol%metperiod, &
           flDayStart => state%timecontrol%flDayStart, &
           flDayEnd => state%timecontrol%flDayEnd, &
           flRunEnd => state%timecontrol%flRunEnd, &
           flYearStart => state%timecontrol%flYearStart, &
           floutput => state%timecontrol%floutput, &
           floutputshort => state%timecontrol%floutputshort, &
           flbaloutput => state%timecontrol%flbaloutput, &
           flheader => state%timecontrol%flheader, &
           flheadirg => state%timecontrol%flheadirg, &
           flIrg1Start => state%timecontrol%flIrg1Start, &
           flUpdMetDet => state%timecontrol%flUpdMetDet, &
           fldecdtmin => state%timecontrol%fldecdtmin, &
           fldtmin => state%timecontrol%fldtmin, &
           fldtreduce => state%timecontrol%fldtreduce, &
           flprintshort => state%timecontrol%flprintshort, &
           flmetdetail => state%timecontrol%flmetdetail, &
           flmeteodt => state%timecontrol%flmeteodt, &
           flrainintens => state%timecontrol%flrainintens, &
           fletsine => state%timecontrol%fletsine, &
           flIrrigate => state%timecontrol%flIrrigate, &
           flDrain => state%timecontrol%flDrain, &
           flSurfaceWater => state%timecontrol%flSurfaceWater, &
           flTemperature => state%timecontrol%flTemperature, &
           flSnow => state%timecontrol%flSnow, &
           flSolute => state%timecontrol%flSolute )

! === initialization ===================================================

! [SS-TC TC-14] iyear/imonth/dt seeded into state%timecontrol from
!   tc_iyear_init_buf / tc_imonth_init_buf / tc_dt_init_buf by swap.f90
!   before calling TimeControl(1). No-op writes here removed.

! --- initialize flags ----------------------------
      fldecdt = .false.
      fldecdtmin = .false.
      flRunEnd = .false.
      flDayStart = .true.
      fldtmin = .false.
      flZeroIntr = .true.
      flZeroCumu = .true.
      state%timecontrol%flZeroIntr = .true.   ! [SS-TCM transition] dual-write; readers cut over Task 11
      state%timecontrol%flZeroCumu = .true.   ! [SS-TCM transition] dual-write; readers cut over Task 11
      floutput = .false.
      flbaloutput = .false.
      flheader = .false.
      flheadirg = .false.
      flIrg1Start = .true.
      flUpdMetDet = .true.
      flYearStart = .true.

      if (nprintday .gt. 1 .or. flprintdt) then
        flprintshort = .true.
        period = 1
      else
        flprintshort = .false.
      endif
      if (swmetdetail .eq. 1) then
        flmetdetail = .true.
        metperiod = 1.0d0 / dble(nmetdetail)
        wrecord = 0
      else
        flmetdetail = .false.
      endif
      if (swrain.gt.0) then
         flrainintens = .true.
      else
         flrainintens = .false.
      endif
      if (flmetdetail .or. flrainintens) then
         flmeteodt = .true.
      else
         flmeteodt = .false.
      endif
      fletsine = .false.
      if (swetsine .eq. 1) fletsine = .true.
      if (swirfix.eq.1) flIrrigate = .true.
      flDrain = .false.
      if (swdra .eq. 1) flDrain = .true.
      flSurfaceWater = .false.
      if (swdra .eq. 2) flSurfaceWater = .true.
      flTemperature = .false.
      if (swhea .eq. 1) flTemperature = .true.
      flSnow = .false.
      if (swsnow .eq. 1) flSnow = .true.
      flSolute = .false.
      if (swsolu .eq. 1) flSolute = .true.

! --- initialize counters ----------------------------
      nirri = 1
      isteps = 0
      ioutdat = 1
      ioutdatint = 1
      cntper = 0
      outper = 0.0d0
      tcumold = 0.0d0
      nprintcount = 1

! --- set main time variable of SWAP ------------------------
      t1900 = tstart

! --- determine time from beginning of calendar year
      call dtdpar (t1900, datea, fsec)
      datea(1) = iyear
      datea(2) = 1
      datea(3) = 1
      fsec = 0.0
      call dtardp (datea,fsec,timjan1)
      t = tstart - timjan1
      tEvent = 0.0d0
      tcum = 0.d0
      daynr = nint(t)
      daycum = 0
      daymeteo = daynr + 1

! --- determine year,month and day number of current day
      call dtdpar (t1900+0.1d0, datea, fsec)
      iyear = datea(1)
      imonth = datea(2)
      yearmeteo = iyear

! --- determine date of current day
      call dtdpst ('year-month-day',t1900+0.1d0,date)

! --- output to screen
      if (swscre .eq. 2) then
        filtext = 'Screen output of daynumbers'
        call writehead (5,1,'screen',filtext,project)
        call dtdpst ('year-month-day',tstart,date)
        write (*,'(2x,2a)') 'First day of simulation:  ',date
        call dtdpst ('year-month-day',tend,date)
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
      tEvent = dtEvent

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
      dtprevious = dt
      if(dt+dtCrit .gt. dtEvent)then
         flprevious = 2
         dt = dtEvent
         flTnext = .true.
      else
         flprevious = 1
         flTnext = .false.
      end if

! --- in case of sine wave of ET, limit dt and dtmax
      if (fletsine) then
        dt = min(dt, dtfletsine)
        dtmax = min(dtmax, dtfletsine)
      endif

! --- set value for dtold, used in headcalc for macropore-iteration
      dtold = dt

      end associate
   end subroutine timecontrol_init

   subroutine timecontrol_advance(state)
      use variables, only: dtmin, dtmax, period, outdat, outdatint, tend, tstart, &
                            nprintday, swheader, swres, swscre, &
                            flCropCalendar, msteps, flzerointr, flzerocumu, &
                            flprintdt, swrain, swmetdetail, &
                            flCropHarvest, flCropOutput, croptype, &
                            raintimearray, dtEventRain, dt_SSDI_event, flSSDI, &
                            numbit, MaxIt, icrop
      use timestep_control_mod, only: fldecdt
      use irrigation_mod, only: SSDI_irrigation
      use error_mod, only: fatalerr_collected
      implicit none
      type(swap_state_t), intent(inout) :: state

      ! local variables
      character(len=200) :: messag
      character(len=11)  :: tmp
      real(8)            :: dtCrit, dtRestDay

      dtCrit = 1.d-8

      ! Associate block — copied verbatim from timecontrol.f90:54-114
      associate( datea => state%timecontrol%datea, &
           nextyear => state%timecontrol%nextyear, &
           flprevious => state%timecontrol%flprevious, &
           flTnext => state%timecontrol%flTnext, &
           fsec => state%timecontrol%fsec, &
           tchange => state%timecontrol%tchange, &
           dtEvent => state%timecontrol%dtEvent, &
           tEvent => state%timecontrol%tEvent, &
           tcumold => state%timecontrol%tcumold, &
           dtprevious => state%timecontrol%dtprevious, &
           tmptimestart => state%timecontrol%tmptimestart, &
           tmptimeend => state%timecontrol%tmptimeend, &
           iyear => state%timecontrol%iyear, &
           iyearm1 => state%timecontrol%iyearm1, &
           imonth => state%timecontrol%imonth, &
           dt => state%timecontrol%dt, &
           dtold => state%timecontrol%dtold, &
           daynr => state%timecontrol%daynr, &
           daycum => state%timecontrol%daycum, &
           daymeteo => state%timecontrol%daymeteo, &
           yearmeteo => state%timecontrol%yearmeteo, &
           t => state%timecontrol%t, &
           t1900 => state%timecontrol%t1900, &
           tcum => state%timecontrol%tcum, &
           timjan1 => state%timecontrol%timjan1, &
           outper => state%timecontrol%outper, &
           cntper => state%timecontrol%cntper, &
           isteps => state%timecontrol%isteps, &
           ioutdat => state%timecontrol%ioutdat, &
           ioutdatint => state%timecontrol%ioutdatint, &
           nprintcount => state%timecontrol%nprintcount, &
           rainrec => state%timecontrol%rainrec, &
           wrecord => state%timecontrol%wrecord, &
           swmeteo => state%timecontrol%swmeteo, &
           date => state%timecontrol%date, &
           metperiod => state%timecontrol%metperiod, &
           flDayStart => state%timecontrol%flDayStart, &
           flDayEnd => state%timecontrol%flDayEnd, &
           flRunEnd => state%timecontrol%flRunEnd, &
           flYearStart => state%timecontrol%flYearStart, &
           floutput => state%timecontrol%floutput, &
           floutputshort => state%timecontrol%floutputshort, &
           flbaloutput => state%timecontrol%flbaloutput, &
           flheader => state%timecontrol%flheader, &
           flheadirg => state%timecontrol%flheadirg, &
           flIrg1Start => state%timecontrol%flIrg1Start, &
           flUpdMetDet => state%timecontrol%flUpdMetDet, &
           fldecdtmin => state%timecontrol%fldecdtmin, &
           fldtmin => state%timecontrol%fldtmin, &
           fldtreduce => state%timecontrol%fldtreduce, &
           flprintshort => state%timecontrol%flprintshort, &
           flmetdetail => state%timecontrol%flmetdetail, &
           flmeteodt => state%timecontrol%flmeteodt, &
           flrainintens => state%timecontrol%flrainintens, &
           fletsine => state%timecontrol%fletsine, &
           flIrrigate => state%timecontrol%flIrrigate, &
           flDrain => state%timecontrol%flDrain, &
           flSurfaceWater => state%timecontrol%flSurfaceWater, &
           flTemperature => state%timecontrol%flTemperature, &
           flSnow => state%timecontrol%flSnow, &
           flSolute => state%timecontrol%flSolute )

! === next time step ===================================================

! 2.1  check maximum number of time steps during this day
      isteps = isteps + 1
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
      tcum = tcum + dt
      t1900 = tstart + tcum

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
        state%timecontrol%flZeroCumu = .false.   ! [SS-TCM transition] dual-write; readers cut over Task 11
        if (flZeroIntr) then
          outper = 0.0d0
          flZeroIntr = .false.
          state%timecontrol%flZeroIntr = .false.   ! [SS-TCM transition] dual-write; readers cut over Task 11
        endif

! ---   set flags for output
        if (floutput) then
         floutput = .false.
         if (swheader .eq. 1) then
           flheader = .false.
         endif
        endif
        if (flbaloutput) then
         flbaloutput = .false.
         if (swheader .eq. 1) then
           flheader = .true.
           flheadirg = .true.
         endif
        endif

! 2.4  determine year,month and day number (only during first time step of a day)
        call dtdpar (t1900, datea, fsec)
        iyearm1 = iyear
        iyear = datea(1)
        imonth = datea(2)

! ---   determine date of current day
        call dtdpst ('year-month-day',t1900,date)

! ---   update day numbers
        daynr = daynr+1
        daycum = daycum + 1
        cntper = cntper + 1

! 2.5  in case of detailed meteorological input, reset weather record

        if (iyear .ne. iyearm1) then
! ---     reset daynumber and time because new calender year has started
          daynr = 1
          t = dt

! ---     in case SWRES = 1 reset counter for periodic output to 1
          if (swres.eq.1 .and. period.ne.0) then
            cntper = 1
          endif
        endif


      endif


! 2.6  update  logicals for indication of start / end of day
      if ( dble(daycum) - tcum .lt. dtCrit) then
        flDayEnd = .true.
        flDayStart = .true.
        flTnext = .true.
        if (flmetdetail) then
          wrecord = 1
        endif
      else
        flDayEnd = .false.
        flDayStart = .false.
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
            if ((t+dtEvent) .gt. tchange) dtEvent = tchange - t
         end if

! 2.7.4  precipitation event may limit timestep
         if (swmetdetail.eq.0 .and. swrain.gt.0) then

!        next rainevent! Set new values
           if (raintimearray(rainrec).lt.tcum+dtCrit) then
!        new rain event valid
             rainrec = rainrec + 1
           endif

           dtEvent = min(dtevent,dtEventRain)
           dtEvent = max(dtEvent,dtmin)
         endif

! 2.8  set end of time interval (determined by I/O)
         tEvent = tEvent + dtEvent
         flTnext = .false.

      endif

! 2.9 determine next time step, based on numerical performance
      if (flDayStart) then
        if(flprevious.eq.2)then
           dt = max(dt, dsqrt(dtmin*dtmax),dtprevious)
        else
           dt = max(dt, dsqrt(dtmin*dtmax))
        end if
        dtprevious = dt
      else
        if(flprevious .eq. 2)then
           dt = dtprevious
        else
           if (numbit.le.3)     dt = min(dt*2.0d0,DtMax)
           if (numbit.ge.MaxIt) dt = max(dt*0.5d0,DtMin)
           dtprevious = dt
         endif
      endif
      flprevious = 1

      if ( tcum + dt - tEvent .gt. dtCrit) then
         dt = tEvent - tcum
         flprevious = 2
         flTnext = .true.
      endif

!     SSDI: adapt dt as to not pass dt_SSDI_event end time of the day
      if (dt_SSDI_event < 1.0d0) then
         dt = max(dtmin, min(dt, dble(int(tcum) + dt_SSDI_event) - tcum))
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

! 2.11 set flags and variables

! --- in case of output during a day
      if (flprintshort) then
        floutputshort = .false.
        flzerointr = .false.
        state%timecontrol%flZeroIntr = .false.   ! [SS-TCM transition] dual-write; readers cut over Task 11
! ---   determine whether output is required
        if (flprintdt) then
           outper = tcum - tcumold
           tcumold = tcum
           if (abs(outdatint(ioutdatint) - t1900 + 1.d0).lt.1.d-3) then
              floutputshort = .true.
              flzerointr = .true.
              state%timecontrol%flZeroIntr = .true.   ! [SS-TCM transition] dual-write; readers cut over Task 11
           endif
        else
           if (tcum+dtCrit .gt. dble(nprintcount)/dble(nprintday)) then
              floutputshort = .true.
              flzerointr = .true.
              state%timecontrol%flZeroIntr = .true.   ! [SS-TCM transition] dual-write; readers cut over Task 11
              outper = tcum - tcumold
              tcumold = tcum
           endif
! ---      update counter nprintcount for printing
           do while (tcum+dtCrit .gt. dble(nprintcount)/dble(nprintday))
              nprintcount = nprintcount + 1
           end do
         endif
      endif

! --- in case of detailed meteorological input
      if (flmetdetail) then
        tchange = dble(int(t + dtCrit)) + dble(wrecord) * metperiod
        if ((tchange - t) .lt. dtCrit) then
! ---     update actual weather record and fluxes

          flUpdMetDet = .true.
        endif
      endif

!     SSDI: end of subsurface irirgation event reached; reset
      if (flSSDI .and. tcum - int(tcum) + dtCrit > dt_SSDI_event) then
         call SSDI_irrigation(9, state)  ! [SS-SWC S-2.12B]
      end if

! --- update fldtmin
      if (dt .gt. (1.0d0+dtCrit)*dtmin) then
        fldtmin = .false.
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
        endif

! ---   write daynumber to screen
        if (swscre .eq. 2) then
          write(*,'("+ ",4x,a11,i6,i8)') date,daynr,daycum
        endif

! ---   end of run?
        if ((tend - t1900 + 1.d0) .lt. 1.d-3) then
          flRunEnd = .true.
          floutput = .true.
          flbaloutput = .true.
          ioutdat = ioutdat + 1
          return
        endif

! ---   in case no end of run, determine whether today output should be written
        if (cntper .eq. period .and. .not.flprintdt) then
          cntper = 0
          floutput = .true.
          flzerointr = .true.
          state%timecontrol%flZeroIntr = .true.   ! [SS-TCM transition] dual-write; readers cut over Task 11
        endif

        if(flprintdt) then
          if (abs(outdatint(ioutdatint) - t1900 + dt) .lt. 1.d-3) then
            floutput = .true.
            flzerointr = .true.
            state%timecontrol%flZeroIntr = .true.   ! [SS-TCM transition] dual-write; readers cut over Task 11
            ioutdatint = ioutdatint + 1
          endif
        else
          if (abs(outdatint(ioutdatint) - t1900 + 1.d0) .lt. 1.d-3) then
            floutput = .true.
            flzerointr = .true.
            state%timecontrol%flZeroIntr = .true.   ! [SS-TCM transition] dual-write; readers cut over Task 11
            ioutdatint = ioutdatint + 1
          endif
        endif
        if (abs(outdat(ioutdat) - t1900 + 1.d0) .lt. 1.d-3) then
! ---     output of water and solute balances
          floutput = .true.
          flbaloutput = .true.
          flzerointr = .true.
          state%timecontrol%flZeroIntr = .true.   ! [SS-TCM transition] dual-write; readers cut over Task 11
          flzerocumu = .true.
          state%timecontrol%flZeroCumu = .true.   ! [SS-TCM transition] dual-write; readers cut over Task 11
          ioutdat = ioutdat + 1
        endif


! ---   reset flags for next day
        fldtmin = .false.

! ---   reset counters for next day
        isteps = 0

! ---   determine daynumber and switch for reading meteorological data
        call dtdpar (t1900 + 0.1d0, datea, fsec)
!        call dtdpar (t1900 , datea, fsec)
        nextyear = datea(1)
        if (nextyear .eq. iyear) then
          daymeteo = daynr + 1
        else
          yearmeteo = nextyear
! ---     set flag for new meteo year
          flYearStart = .true.

          daymeteo = 1
! ---     detailed meteo data needed for crop growth?
          swmeteo = 1
          if (flCropCalendar) then
            if (icrop .gt. 0) then
              if (croptype(icrop).ge.2) then
                swmeteo = 2
              endif
            endif
          endif
        endif

      endif

      dtprevious = dt

      end associate
   end subroutine timecontrol_advance

   subroutine timecontrol_reduce_dt(state)
      type(swap_state_t), intent(inout) :: state
      ! Body filled in Task 7 (migrated from timecontrol.f90 case (3)).
   end subroutine timecontrol_reduce_dt

   subroutine timecontrol_day_end(state)
      type(swap_state_t), intent(inout) :: state
      ! Body filled in Task 8 (migrated from timecontrol.f90 case (9)).
   end subroutine timecontrol_day_end

   subroutine itertime_init(state)
      type(swap_state_t), intent(inout) :: state
      ! Body filled in Task 9 (migrated from IterTime case (1)).
   end subroutine itertime_init

   subroutine itertime_check(state)
      type(swap_state_t), intent(inout) :: state
      ! Body filled in Task 9 (migrated from IterTime case (2)).
   end subroutine itertime_check

   subroutine itertime_close(state)
      type(swap_state_t), intent(inout) :: state
      ! Body filled in Task 9 (migrated from IterTime case (3)).
   end subroutine itertime_close

end module timecontrol_mod
