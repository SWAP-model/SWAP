!> @file timecontrol_mod.f90
!! SS-TCM: module form of the legacy subroutine TimeControl + IterTime.
!! Seven named lifecycle procedures replace the magic-int dispatch
!! (task=1/2/3/9 + IterTime task=1/2/3). State threaded explicitly;
!! each procedure carries its own associate block. flZeroIntr and
!! flZeroCumu are owned by state%timecontrol; bare globals
!! flzerointr/flzerocumu retire in Task 12 (variables.f90/initialize.f90).
module timecontrol_mod
   use swap_state_mod, only: swap_state_t
   implicit none
   private
   public :: timecontrol_init, timecontrol_advance, &
             timecontrol_reduce_dt, timecontrol_day_end
   public :: itertime_init, itertime_check, itertime_close

contains

   subroutine timecontrol_init(state)
      ! [GR-CROP C1] flCropCalendar/icrop: dual-write to both legacy global + state%crop%common%X
      ! [SS-GR-FINAL B8] DEFERRED: all symbols — init-path config/flags; no state home yet; Phase C3
      use swap_log, only: log_warn
      ! [GR-CROP 2026-05-25] flCropCalendar/icrop/cropstart/project → state/config reads.
      ! [GR-TIME 2026-05-25] swirfix/swsnow/swhea/swsolu/swetsine/nirri now read
      ! from state%cfg/state%crop%irrigation; bare-global use-variables retired.
      use error_mod, only: fatalerr_collected
      implicit none
      type(swap_state_t), intent(inout) :: state

      ! local variables
      character(len=200) :: messag
      character(len=80)  :: filtext
      real(8)            :: dtCrit, dtfletsine

      dtCrit = 1.d-8
      dtfletsine = 0.05d0

      associate (time      => state%timecontrol,   &
                 soil      => state%soilwater,     &
                 crop      => state%crop,          &
                 meteo_cfg => state%cfg%meteo,     &
                 crop_cfg  => state%cfg%crop)

! === initialization ===================================================

! --- initialize flags ----------------------------
      time%fldecdt = .false.
      time%fldecdtmin = .false.
      time%flRunEnd = .false.
      time%flDayStart = .true.
      time%fldtmin = .false.
      time%flZeroIntr = .true.   ! reset gate: intermediate accumulators
      time%flZeroCumu = .true.   ! reset gate: cumulative accumulators
      time%floutput = .false.
      time%flbaloutput = .false.
      time%flheader = .false.
      time%flheadirg = .false.
      time%flIrg1Start = .true.
      time%flUpdMetDet = .true.
      time%flYearStart = .true.

      if (time%nprintday .gt. 1 .or. time%flprintdt) then
        time%flprintshort = .true.
        time%period = 1
      else
        time%flprintshort = .false.
      endif
      if (meteo_cfg%swmetdetail .eq. 1) then
        time%flmetdetail = .true.
        time%metperiod = 1.0d0 / dble(meteo_cfg%nmetdetail)
        time%wrecord = 0
      else
        time%flmetdetail = .false.
      endif
      if (meteo_cfg%swrain.gt.0) then
         time%flrainintens = .true.
      else
         time%flrainintens = .false.
      endif
      if (time%flmetdetail .or. time%flrainintens) then
         time%flmeteodt = .true.
      else
         time%flmeteodt = .false.
      endif
      time%fletsine = .false.
      if (meteo_cfg%swetsine .eq. 1) time%fletsine = .true.
      if (state%cfg%irrigation%swirfix .eq. 1) time%flIrrigate = .true.
      ! flDrain and flSurfaceWater are seeded by timecontrol_state_init
      ! (called from swap_init_body via state%timecontrol%init) which
      ! reads config_drain%swdra directly — no state%cfg read needed here.
      time%flTemperature = .false.
      if (state%cfg%heat%swhea .eq. 1) time%flTemperature = .true.
      time%flSnow = .false.
      if (meteo_cfg%snow%swsnow .eq. 1) time%flSnow = .true.
      time%flSolute = .false.
      if (state%cfg%solute%swsolu .eq. 1) time%flSolute = .true.

! --- initialize counters ----------------------------
      time%isteps = 0
      time%ioutdat = 1
      time%ioutdatint = 1
      time%cntper = 0
      time%outper = 0.0d0
      time%tcumold = 0.0d0
      time%nprintcount = 1

! --- set main time variable of SWAP ------------------------
      time%t1900 = time%tstart

! --- determine time from beginning of calendar year
      call dtdpar (time%t1900, time%datea, time%fsec)
      time%datea(1) = time%iyear
      time%datea(2) = 1
      time%datea(3) = 1
      time%fsec = 0.0
      call dtardp (time%datea, time%fsec, time%timjan1)
      time%t = time%tstart - time%timjan1
      time%tEvent = 0.0d0
      time%tcum = 0.d0
      time%daynr = nint(time%t)
      time%daycum = 0
      time%daymeteo = time%daynr + 1

! --- determine year,month and day number of current day
      call dtdpar (time%t1900+0.1d0, time%datea, time%fsec)
      time%iyear = time%datea(1)
      time%imonth = time%datea(2)
      time%yearmeteo = time%iyear

! --- determine date of current day
      call dtdpst ('year-month-day', time%t1900+0.1d0, time%date)

! --- output to screen
      if (time%swscre .eq. 2) then
        filtext = 'Screen output of daynumbers'
        call writehead (5,1,'screen',filtext,state%cfg%general%project)
        call dtdpst ('year-month-day', time%tstart, time%date)
        write (*,'(2x,2a)') 'First day of simulation:  ', time%date
        call dtdpst ('year-month-day', time%tend, time%date)
        write (*,'(2x,2a)') 'Last day of simulation:   ', time%date
        write (*,'(/,a,/)') '            date  daynr  daycum'
      endif

! -   set crop number
      crop%common%icrop = 1
      do while (.not. crop%common%flCropCalendar)

        if (crop_cfg%rotation_start(crop%common%icrop) .lt. 1.d0) exit

        if (abs(time%tstart - crop_cfg%rotation_start(crop%common%icrop)) .lt. 1.d-3) then
          crop%common%flCropCalendar = .true.
          if (time%tstart - crop_cfg%rotation_start(crop%common%icrop) .lt. -1.d-3 .and. &
     &                                              soil%swinco .ne. 3) then
            messag = 'The start of simulation (tstart) begins in '//  &
     &      'crop growing season with swinco 1 or 2'
            call fatalerr_collected ('readswap',messag)
          endif
        else
          crop%common%icrop = crop%common%icrop + 1
        end if
      enddo

! --- detailed meteo data needed for crop growth?
      time%swmeteo = 1
      if (crop%common%flCropCalendar) then
        if (crop%common%icrop .gt. 0) then
          if (crop%common%croptype(crop%common%icrop).ge.2) then
             time%swmeteo = 2
          endif
        endif
      endif

! --- initialize
      time%dtEvent = 1.0d0

! --- limit initial dt in case of short time interval
      if (time%flprintshort .and. .not.time%flprintdt) then
        if (time%dtEvent .gt. dble(time%nprintcount)/dble(time%nprintday))then
           time%dtEvent = dble(time%nprintcount)/dble(time%nprintday)
         end if
      endif

! --- intial timestep for rainfall intensities
      if (meteo_cfg%swmetdetail.eq.0 .and. meteo_cfg%swrain.gt.0) then
         time%dtEvent = min(time%dtEvent, time%dtmin)
      endif

! --- limit initial dt in case of detailed meteorological input
       if (time%flmetdetail) then
        if (time%dtEvent .gt. time%metperiod)then
           time%dtEvent = min(time%dtEvent, time%metperiod)
         end if
      endif
      time%tEvent = time%dtEvent

! --- initial time step
      if (soil%swinco.eq.3) then
        if (time%dt.lt.time%dtmin) then
          messag = 'Initial dt read from file (SWINCO=3)'//             &
     &    ' if absent, then default is assumed'
          call log_warn('soilwater', messag)
          time%dt = dsqrt(time%dtmin*time%dtmax)
        endif
      else
        time%dt = dsqrt(time%dtmin*time%dtmax)
      endif
      time%dtprevious = time%dt
      if (time%dt + dtCrit .gt. time%dtEvent) then
         time%flprevious = 2
         time%dt = time%dtEvent
         time%flTnext = .true.
      else
         time%flprevious = 1
         time%flTnext = .false.
      end if

! --- in case of sine wave of ET, limit dt and dtmax
      if (time%fletsine) then
        time%dt = min(time%dt, dtfletsine)
        time%dtmax = min(time%dtmax, dtfletsine)
      endif

! --- set value for dtold, used in headcalc for macropore-iteration
      time%dtold = time%dt

      end associate
   end subroutine timecontrol_init

   subroutine timecontrol_advance(state)
      ! [GR-CROP Phase B] raintimearray migrated → state%atmosphere%raintimearray via associate.
      ! [GR-CROP C1] flCropCalendar/flCropOutput/icrop: reads from state%crop%common%X; writes dual to state+legacy
      ! [SS-GR-FINAL B8] DEFERRED: all residual symbols; Phase C3
      ! [GR-CROP 2026-05-25] flCropCalendar/flCropOutput/flCropHarvest/icrop → state%crop%common.
      ! [GR-TIME 2026-05-25] outdat/outdatint migrated to state%timecontrol;
      ! flSSDI bare global replaced with `state%cfg%irrigation%swssdi == 1`.

      use irrigation_mod, only: ssdi_irrigation_reset
      use error_mod, only: fatalerr_collected
      implicit none
      type(swap_state_t), intent(inout) :: state

      ! local variables
      character(len=200) :: messag
      character(len=11)  :: tmp
      real(8)            :: dtCrit, dtRestDay

      dtCrit = 1.d-8

      associate (time      => state%timecontrol,   &
                 soil      => state%soilwater,     &
                 crop      => state%crop,          &
                 atmo      => state%atmosphere,    &
                 meteo_cfg => state%cfg%meteo)

! === next time step ===================================================

! 2.1  check maximum number of time steps during this day
      time%isteps = time%isteps + 1
      if (time%isteps .gt. time%msteps) then
        write(tmp,'(i11)') time%daynr
        tmp = adjustl(tmp)
        messag ='The maximum number of time steps for a day is exceeded'&
     &    //' at daynumber '//trim(tmp)//'. Check input for numerical'  &
     &    //' solution of Richards equation'
        call fatalerr_collected ('timer',messag)
      endif


! 2.2 update time variables

      time%t = time%t + time%dt           ! relative to yyyy0101:00:00:00
      time%tcum = time%tcum + time%dt
      time%t1900 = time%tstart + time%tcum

! 2.3  flag assignments after first time step of a day

      if (time%flDayStart) then

! ---   set crop conditions  [GR-CROP 2026-05-25] state-only
        if (crop%common%flCropCalendar) then
          crop%common%flCropOutput = .true.
          if (crop%common%flCropHarvest) then
            crop%common%flCropOutput = .false.
          endif
        endif

!! --    length crop
!        if (flCropCalendar) then
!          daycrop = daycrop + 1
!        endif

! ---   set flags for reset intermediate and cumulative fluxes
        time%flZeroCumu = .false.
        if (time%flZeroIntr) then
          time%outper = 0.0d0
          time%flZeroIntr = .false.
        endif

! ---   set flags for output
        if (time%floutput) then
         time%floutput = .false.
         if (time%swheader .eq. 1) then
           time%flheader = .false.
         endif
        endif
        if (time%flbaloutput) then
         time%flbaloutput = .false.
         if (time%swheader .eq. 1) then
           time%flheader = .true.
           time%flheadirg = .true.
         endif
        endif

! 2.4  determine year,month and day number (only during first time step of a day)
        call dtdpar (time%t1900, time%datea, time%fsec)
        time%iyearm1 = time%iyear
        time%iyear = time%datea(1)
        time%imonth = time%datea(2)

! ---   determine date of current day
        call dtdpst ('year-month-day', time%t1900, time%date)

! ---   update day numbers
        time%daynr = time%daynr + 1
        time%daycum = time%daycum + 1
        time%cntper = time%cntper + 1

! 2.5  in case of detailed meteorological input, reset weather record

        if (time%iyear .ne. time%iyearm1) then
! ---     reset daynumber and time because new calender year has started
          time%daynr = 1
          time%t = time%dt

! ---     in case SWRES = 1 reset counter for periodic output to 1
          if (time%swres.eq.1 .and. time%period.ne.0) then
            time%cntper = 1
          endif
        endif


      endif


! 2.6  update  logicals for indication of start / end of day
      if (dble(time%daycum) - time%tcum .lt. dtCrit) then
        time%flDayEnd = .true.
        time%flDayStart = .true.
        time%flTnext = .true.
        if (time%flmetdetail) then
          time%wrecord = 1
        endif
      else
        time%flDayEnd = .false.
        time%flDayStart = .false.
      endif


! 2.7  maximum size of time interval (dtEvent)

      if (time%flTnext) then

! 2.7.1  remaining part of a day
        time%dtEvent = dble(int(time%tcum+1.0d0+dtCrit)) - time%tcum

! 2.7.2  printing more than one times a day may limit timestep
         if (time%flprintshort .and. .not.time%flprintdt) then
           if (dble(time%nprintcount)/dble(time%nprintday) - time%tcum .lt. dtCrit) then
              time%dtEvent = min(time%dtEvent, 1.0d0/dble(time%nprintday))
           else
              time%dtEvent = min(time%dtEvent,                                    &
     &                   (dble(time%nprintcount)/dble(time%nprintday)) - time%tcum)
           endif
         endif

! 2.7.3  input of detailed meteo may limit timestep
         if (time%flmetdetail) then
            time%tchange = dble(int(time%t + dtCrit)) + dble(time%wrecord) * time%metperiod
            if ((time%t + time%dtEvent) .gt. time%tchange) time%dtEvent = time%tchange - time%t
         end if

! 2.7.4  precipitation event may limit timestep
         if (meteo_cfg%swmetdetail.eq.0 .and. meteo_cfg%swrain.gt.0) then

!        next rainevent! Set new values
           if (atmo%raintimearray(time%rainrec) .lt. time%tcum + dtCrit) then
!        new rain event valid
             time%rainrec = time%rainrec + 1
           endif

           time%dtEvent = min(time%dtEvent, atmo%dtEventRain)
           time%dtEvent = max(time%dtEvent, time%dtmin)
         endif

! 2.8  set end of time interval (determined by I/O)
         time%tEvent = time%tEvent + time%dtEvent
         time%flTnext = .false.

      endif

! 2.9 determine next time step, based on numerical performance
      if (time%flDayStart) then
        if (time%flprevious .eq. 2) then
           time%dt = max(time%dt, dsqrt(time%dtmin*time%dtmax), time%dtprevious)
        else
           time%dt = max(time%dt, dsqrt(time%dtmin*time%dtmax))
        end if
        time%dtprevious = time%dt
      else
        if (time%flprevious .eq. 2) then
           time%dt = time%dtprevious
        else
           if (soil%numbit.le.3)        time%dt = min(time%dt*2.0d0, time%DtMax)  ! [GR-SOIL 2026-05-24]
           if (soil%numbit.ge.time%MaxIt) time%dt = max(time%dt*0.5d0, time%DtMin)  ! [GR-SOIL 2026-05-24]
           time%dtprevious = time%dt
         endif
      endif
      time%flprevious = 1

      if (time%tcum + time%dt - time%tEvent .gt. dtCrit) then
         time%dt = time%tEvent - time%tcum
         time%flprevious = 2
         time%flTnext = .true.
      endif

!     SSDI: adapt dt as to not pass dt_SSDI_event end time of the day
      if (crop%irrigation%dt_SSDI_event < 1.0d0) then
         time%dt = max(time%dtmin, min(time%dt, dble(int(time%tcum) + crop%irrigation%dt_SSDI_event) - time%tcum))
      end if

! 2.10  test last time step of the day: limit dt if it exceeds end of day
      dtRestDay = dble(int(time%tcum+1.0d0+dtCrit)) - time%tcum
      if (dtRestDay - time%dt .lt. 1.d-6) then
        time%dt = dtRestDay
      endif

! JK20131230: when dtmin is large then the timestep-closure may not be correct and errors may occur
! resulting in water balance errors see email Paul v Walsum 210131225. Elimination of next statement
! is not the correct solution
      time%dt = max(time%dt, time%dtmin)
      time%dt = min(time%dt, time%dtmax)

! 2.11 set flags and variables

! --- in case of output during a day
      if (time%flprintshort) then
        time%floutputshort = .false.
        time%flZeroIntr = .false.
! ---   determine whether output is required
        if (time%flprintdt) then
           time%outper = time%tcum - time%tcumold
           time%tcumold = time%tcum
           if (abs(time%outdatint(time%ioutdatint) - time%t1900 + 1.d0).lt.1.d-3) then
              time%floutputshort = .true.
              time%flZeroIntr = .true.
           endif
        else
           if (time%tcum + dtCrit .gt. dble(time%nprintcount)/dble(time%nprintday)) then
              time%floutputshort = .true.
              time%flZeroIntr = .true.
              time%outper = time%tcum - time%tcumold
              time%tcumold = time%tcum
           endif
! ---      update counter nprintcount for printing
           do while (time%tcum + dtCrit .gt. dble(time%nprintcount)/dble(time%nprintday))
              time%nprintcount = time%nprintcount + 1
           end do
         endif
      endif

! --- in case of detailed meteorological input
      if (time%flmetdetail) then
        time%tchange = dble(int(time%t + dtCrit)) + dble(time%wrecord) * time%metperiod
        if ((time%tchange - time%t) .lt. dtCrit) then
! ---     update actual weather record and fluxes

          time%flUpdMetDet = .true.
        endif
      endif

!     SSDI: end of subsurface irirgation event reached; reset
      if (state%cfg%irrigation%swssdi == 1 .and. &
          time%tcum - int(time%tcum) + dtCrit > crop%irrigation%dt_SSDI_event) then
         call ssdi_irrigation_reset(state)  ! [SS-SWC S-2.12B]
      end if

! --- update fldtmin
      if (time%dt .gt. (1.0d0+dtCrit)*time%dtmin) then
        time%fldtmin = .false.
      endif


! --- procedure when day is finished ----------------------------------------------
      if (time%flDayEnd) then

!! --    length crop
!        if (flCropCalendar) then
!          daycrop = daycrop + 1
!        endif

! ---   length output period
        if (.not.time%flprintshort) then
          time%outper = time%outper + 1.0d0
        endif

! ---   write daynumber to screen
        if (time%swscre .eq. 2) then
          write(*,'("+ ",4x,a11,i6,i8)') time%date, time%daynr, time%daycum
        endif

! ---   end of run?
        if ((time%tend - time%t1900 + 1.d0) .lt. 1.d-3) then
          time%flRunEnd = .true.
          time%floutput = .true.
          time%flbaloutput = .true.
          time%ioutdat = time%ioutdat + 1
          return
        endif

! ---   in case no end of run, determine whether today output should be written
        if (time%cntper .eq. time%period .and. .not.time%flprintdt) then
          time%cntper = 0
          time%floutput = .true.
          time%flZeroIntr = .true.
        endif

        if (time%flprintdt) then
          if (abs(time%outdatint(time%ioutdatint) - time%t1900 + time%dt) .lt. 1.d-3) then
            time%floutput = .true.
            time%flZeroIntr = .true.
            time%ioutdatint = time%ioutdatint + 1
          endif
        else
          if (abs(time%outdatint(time%ioutdatint) - time%t1900 + 1.d0) .lt. 1.d-3) then
            time%floutput = .true.
            time%flZeroIntr = .true.
            time%ioutdatint = time%ioutdatint + 1
          endif
        endif
        if (abs(time%outdat(time%ioutdat) - time%t1900 + 1.d0) .lt. 1.d-3) then
! ---     output of water and solute balances
          time%floutput = .true.
          time%flbaloutput = .true.
          time%flZeroIntr = .true.
          time%flZeroCumu = .true.
          time%ioutdat = time%ioutdat + 1
        endif


! ---   reset flags for next day
        time%fldtmin = .false.

! ---   reset counters for next day
        time%isteps = 0

! ---   determine daynumber and switch for reading meteorological data
        call dtdpar (time%t1900 + 0.1d0, time%datea, time%fsec)
!        call dtdpar (t1900 , datea, fsec)
        time%nextyear = time%datea(1)
        if (time%nextyear .eq. time%iyear) then
          time%daymeteo = time%daynr + 1
        else
          time%yearmeteo = time%nextyear
! ---     set flag for new meteo year
          time%flYearStart = .true.

          time%daymeteo = 1
! ---     detailed meteo data needed for crop growth?
          time%swmeteo = 1
          if (crop%common%flCropCalendar) then
            if (crop%common%icrop .gt. 0) then
              if (crop%common%croptype(crop%common%icrop).ge.2) then
                time%swmeteo = 2
              endif
            endif
          endif
        endif

      endif

      time%dtprevious = time%dt

      end associate
   end subroutine timecontrol_advance

   subroutine timecontrol_reduce_dt(state)
      implicit none
      type(swap_state_t), intent(inout) :: state

      associate (time => state%timecontrol)

! === reduce time step ===================================================

! --- decrease time step in case of no convergence in headcalc
      if (time%fldecdt) then
        if (time%dt .gt. 3.0*time%dtmin) then
          time%dt = time%dt / 3.0
!         force dt to equal multiple dtmin to prevent very small dt-values at end of day
!          dt = dtmin * dble(max(1,int(dt/dtmin)))
        else
          time%dt = time%dtmin
          time%fldtmin = .true.
        endif
        time%fldecdt = .false.
        time%flprevious = 1
        time%dtprevious = time%dt
        time%flTnext = .false.

        return
      endif

! --- decrease time step to dtmin if required by boundtop
      if (time%fldecdtmin) then
        time%dt = time%dtmin
        time%fldtmin = .true.
        time%fldecdtmin = .false.
        time%flprevious = 1
        time%dtprevious = time%dt
        return
      endif

      end associate
   end subroutine timecontrol_reduce_dt

   subroutine timecontrol_day_end(state)
      implicit none
      type(swap_state_t), intent(inout) :: state

!     Special: at end of day, the possible initial time step for next day may be too large; adapt if necessary.
      if (state%crop%irrigation%dt_SSDI_event < 1.0d0) then
         state%timecontrol%dt = min(state%timecontrol%dt, state%crop%irrigation%dt_SSDI_event)
      end if
   end subroutine timecontrol_day_end

   subroutine itertime_init(state)
      implicit none
      type(swap_state_t), intent(inout) :: state
      call cpu_time(state%timecontrol%tmptimestart)
   end subroutine itertime_init

   subroutine itertime_check(state)
      use error_mod, only: fatalerr_collected
      implicit none
      type(swap_state_t), intent(inout) :: state
      real(4) :: tmptimeinterrupt
      integer :: timediff
      character(len=400) :: messag

      call cpu_time(tmptimeinterrupt)
      timediff = int(tmptimeinterrupt) - state%timecontrol%MaxIterTime
      if (timediff > 0) then
         write(messag,'(a,i10,3a)') &
            'The maximum cpu time of ', state%timecontrol%MaxIterTime, ' (secs)', &
            ' was exceeded.  Therefore simulation was interrupted'
         call fatalerr_collected('IterTime', messag)
      end if
   end subroutine itertime_check

   subroutine itertime_close(state)
      use swap_log,  only: log_info, to_str
      implicit none
      type(swap_state_t), intent(inout) :: state
      integer :: i

      associate (time => state%timecontrol, &
                 soil => state%soilwater)

         call log_info('itertime', 'Iteration statistics')
         call log_info('itertime', 'Maximum number of iterations: ' // to_str(time%MaxIt))
         call log_info('itertime', 'It Numb  No of Hits  Tot BTr cycles')
         if (allocated(soil%Itnumb)) then
            do i = 1, 100
               if (soil%Itnumb(i,1) > 0) &
                  call log_info('itertime', to_str(i) // '  ' // to_str(soil%Itnumb(i,1)) &
                     // '  ' // to_str(soil%Itnumb(i,2)))
            end do
         end if

         call cpu_time(time%tmptimeend)
         call log_info('itertime', 'Run-time: ' // &
                       to_str(time%tmptimeend - time%tmptimestart) // ' sec')

      end associate
   end subroutine itertime_close

end module timecontrol_mod
