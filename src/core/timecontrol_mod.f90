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
      type(swap_state_t), intent(inout) :: state
      ! Body filled in Task 6 (migrated from timecontrol.f90 case (2)).
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
