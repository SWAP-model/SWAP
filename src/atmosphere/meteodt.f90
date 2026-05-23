!> Meteorological data processing module for sub-daily time steps
!!
!! This module provides routines for processing meteorological data at sub-daily
!! time steps, including rainfall event processing, detailed meteorological input
!! handling, and diurnal distribution of evapotranspiration.
!!
!! ## Module Contents
!!
!! ### Public Routines
!! - [[MeteoDT]]: Main coordinator for sub-daily meteorological processing
!!
!! ### Private Routines
!! - [[ProcessRainEvents]]: Process rainfall events for an entire calendar year
!! - [[ProcessMeteoTsteps]]: Update meteorological fluxes for current time step
!! - [[ETSine]]: Distribute potential ET according to diurnal sine wave
!!
!! ## Calling Sequence
!! 1. [[MeteoDT]] (main routine) - called in SWAP and ReadMeteo (optional)
!! 2. [[ProcessRainEvents]] - processes input data on rain events; called in MeteoDT (optional)
!! 3. [[ProcessMeteoTsteps]] - processes meteo input data per dt; called in MeteoDT (optional)
!! 4. [[ETSine]] - distributes potential transpiration & evaporation according to sine wave; called in MeteoDT (optional)
module meteodt_mod

   use, intrinsic :: iso_fortran_env, only: real64
   use swap_state_mod, only: swap_state_t
   implicit none
   private
   public :: MeteoDT

contains

   !> Main coordinator for sub-daily meteorological processing
  !!
  !! This subroutine orchestrates the processing of meteorological data at sub-daily
  !! time steps. It handles:
  !! - Rain event processing at the beginning of each year (for swrain 1-3)
  !! - Updating meteorological fluxes at each time step
  !! - Distributing potential ET according to diurnal sine wave (if enabled)
  !!
  !! ## Processing Sequence
  !! 1. **Year initialization**: Process rain events for entire year (if flYearStart=true)
  !! 2. **Time step processing**: Update rain/meteo fluxes for current time step (if flMeteoDT=true)
  !! 3. **ET distribution**: Apply sine wave distribution to potential ET (if flETSine=true)
  !!
  !! Control flags determine which processing steps are active:
  !! - `flRainIntens`: Use rainfall intensities/durations (swrain 1-3)
  !! - `flMeteoDT`: Process meteorological data at sub-daily time steps
  !! - `flETSine`: Distribute ET according to photoperiod-based sine wave
  !!
  !! @note
  !! **Original documentation:**
  !! Last modified: February 2014
  !! Purpose: Returns meteorological fluxes of current day or of parts of a day
  !! (detailed meteo input)
  !! @endnote
   subroutine MeteoDT(state)
      implicit none

      type(swap_state_t), intent(inout) :: state

      associate (time => state%timecontrol)

      ! Beginning of year: process rain events
      if (time%flYearStart .and. time%flrainintens) then
         call ProcessRainEvents(state)
         time%flYearStart = .false.
      end if

      ! Per timestep: update rain or meteo fluxes
      if (time%flmeteodt) call ProcessMeteoTsteps(state)

      ! Distribute potential transpiration and evaporation according to sine wave
      if (time%fletsine) call ETSine(state)

      end associate

   end subroutine MeteoDT

   !> Process rainfall events for an entire calendar year
  !!
  !! Converts daily or event-based rainfall input into time-stamped rain flux arrays
  !! that can be used for sub-daily time step calculations. The processing differs
  !! depending on the rainfall input mode (swrain).
  !!
  !! ## Processing by Input Mode
  !!
  !! ### swrain=1 (Daily sums + mean intensities)
  !! - Uses lookup table (raintab) to determine rainfall intensity for each day
  !! - Calculates duration from: duration = amount / intensity
  !! - Creates timestamped events with uniform intensity during each rain period
  !!
  !! ### swrain=2 (Daily sums + durations)
  !! - Uses specified wet period durations (wet array)
  !! - Distributes daily amount uniformly over specified duration
  !!
  !! ### swrain=3 (Event-based input)
  !! - Processes pre-defined rain events from input file
  !! - Calculates daily totals (arai) by accumulating events
  !! - Handles events spanning multiple days by weighted allocation
  !!
  !! ## Output Arrays
  !! Creates two parallel arrays for rain flux time series:
  !! - `raintimearray(i)`: Time of event closure (days since start)
  !! - `rainfluxarray(i)`: Rainfall intensity during period [i-1, i] (cm/d)
  !!
  !! The flux at time i applies to the interval from time i-1 to time i.
  !!
  !! @warning Must be called at the beginning of each calendar year
  !! @warning Arrays sized to accommodate entire year of events
  !!
  !! @note
  !! **Original documentation:**
  !! Last modified: February 2014
  !! Purpose: Process rain events of one calendar year
  !! Interface:
  !! - I: swrain, yearmeteo, dtmin, raintab, tcum, tend, tstart, wet,
  !!      nmrain, timjan1, rainamount, rainrec
  !! - O: arai, rainfluxarray, raintimearray
  !! @endnote
   subroutine ProcessRainEvents(state)
      use array_utils,           only: afgen
      use swap_array_dimensions, only: mrain
      implicit none

      type(swap_state_t), intent(inout) :: state

      integer :: i, iendyear, j, l, nlack, nn, rday, rdaya(367), rdayold
      real(8) :: araihlp(367), day(mrain), rainam(mrain), rainflux
      real(8) :: raintime, ratimar(mrain), tendyear, vsmall, wght, wwet(368)
      vsmall = 1.0d-8

      associate (atmo  => state%atmosphere,    &
                 time  => state%timecontrol,   &
                 meteo => state%cfg%meteo)

      ! === Process rain events on yearly basis ===

      ! For rain options 1 and 2: convert daily rain quantities and intensities or durations
      ! into rain events by creating raintime and rainflux arrays conform rain option 3
      if (meteo%swrain .eq. 1 .or. meteo%swrain .eq. 2) then
         time%rainrec = 1
         do i = 1, atmo%nmrain
            if (atmo%raintimearray(i + 1) .gt. time%tstart - vsmall) then
               time%rainrec = time%rainrec + 1
               ! Beginning (00:00) of days of current year within simulation period
               day(time%rainrec)    = atmo%raintimearray(i + 1) - time%timjan1 + 1.d0
               rainam(time%rainrec) = 0.1d0 * atmo%rainamount(i)   ! mm → cm
               wwet(time%rainrec)   = atmo%wet(i)
            end if
         end do

         ! Set first record of raintime and rainflux (= 0)
         atmo%raintimearray(1) = time%tcum + time%dtmin
         rainam(1)             = 0.d0
         atmo%rainfluxarray(1) = 0.d0

         ! Set rest of records of raintime and rainflux (only when rainam[ount] > 0)
         atmo%nmrain  = time%rainrec + 1
         time%rainrec = 0
         do i = 2, atmo%nmrain
            if (rainam(i) .gt. vsmall) then
               if (meteo%swrain .eq. 1) then
                  ! Mean rainfall intensities are specified
                  rainflux = afgen(meteo%raintab, 60, day(i))
                  raintime = min(0.99d0, rainam(i)/rainflux)
               elseif (meteo%swrain .eq. 2) then
                  ! Rainfall durations are specified
                  raintime = wwet(i)
               end if

               if (i .eq. 2) then
                  time%rainrec = time%rainrec + 1
               else
                  ! First raintime of a day: closure of last period of former day with rain = 0
                  time%rainrec = time%rainrec + 2
                  atmo%raintimearray(time%rainrec) = real(i - 2, real64) + time%tcum
                  atmo%rainfluxarray(time%rainrec) = 0.d0
               end if
               ! Second raintime of a day: closure of first period of the day, rain = rainam
               atmo%raintimearray(time%rainrec + 1) = real(i - 2, real64) + time%tcum + raintime
               atmo%rainfluxarray(time%rainrec + 1) = rainam(i)/raintime
            end if
         end do

         ! Extend array with records at end of current year
         tendyear = 365.d0
         if (mod(time%yearmeteo, 4) .eq. 0) tendyear = 366.d0
         atmo%raintimearray(time%rainrec + 2) = time%tcum + tendyear + time%dtmin
         atmo%rainfluxarray(time%rainrec + 2) = 0.d0

      elseif (meteo%swrain .eq. 3) then
         ! Rain events: 1) calculate daily values, 2) fill raintimearray and rainfluxarray

         ! Total amount of rain per meteo day arai — init array with sum of rain
         do i = 1, 366
            araihlp(i) = 0.d0
         end do

         ! Less rain days than meteo days? Fill gap with dummies
         rdayold = int(atmo%raintimearray(1) - time%timjan1) + 1
         nlack = rdayold - 1
         do j = 1, nlack
            rdaya(j) = j
            araihlp(j) = 0.d0
         end do
         rdaya(j) = rdayold
         araihlp(j) = 0.d0

         ! Fill array of daily sums of rain with real values
         do i = 1, atmo%nmrain
            rday = int(atmo%raintimearray(i) - time%timjan1) + 1
            if (rday .gt. rdayold) then
               do l = 1, rday - rdayold - 1
                  j = j + 1
                  rdaya(j) = rdaya(j - 1) + 1
                  araihlp(j) = 0.d0
               end do
               ! For an event exceeding current day, weight parts between days
               wght = (1.d0 - (atmo%raintimearray(i - 1) - real(int(atmo%raintimearray(i - 1)), real64))) / &
                      (atmo%raintimearray(i) - atmo%raintimearray(i - 1))
               araihlp(j) = araihlp(j) + atmo%rainamount(i)*wght
               rdayold = rday
               j = j + 1
               rdaya(j) = rday
               araihlp(j) = araihlp(j) + atmo%rainamount(i)*(1.d0 - wght)
            else
               if (i .gt. 1) araihlp(j) = araihlp(j) + atmo%rainamount(i)
            end if
         end do

         ! Rain days missing at the end of the year? Fill gap with dummies
         iendyear = 365
         if (mod(time%yearmeteo, 4) .eq. 0) iendyear = 366
         nlack = iendyear - j
         do i = 1, nlack
            rdaya(j + i) = j + i
            araihlp(j + i) = 0.d0
         end do

         ! Save help array araihlp into arai
         do i = 1, 366
            atmo%arai(i) = araihlp(i)
         end do

         ! Build raintimearray and rainam for calculating rainfluxarray
         ! Find time gap without rain events at the beginning of the year
         time%rainrec = 1
         ratimar(1) = atmo%raintimearray(1) - time%tstart
         i = 1
         do while (ratimar(i) .lt. vsmall)
            i = i + 1
            ratimar(i) = atmo%raintimearray(i) - time%tstart
         end do
         do while (atmo%rainamount(i) .lt. vsmall .and. atmo%rainamount(i + 1) .lt. vsmall)
            i = i + 1
            ratimar(i) = atmo%raintimearray(i) - time%tstart
         end do

         ! Fill arrays with real values of rain events
         nn = i
         do i = nn, atmo%nmrain
            ratimar(i + 1) = atmo%raintimearray(i + 1) - time%tstart
            time%rainrec = time%rainrec + 1
            atmo%raintimearray(time%rainrec) = ratimar(i)
            rainam(time%rainrec)             = 0.1d0 * atmo%rainamount(i)   ! mm → cm
         end do
         atmo%nmrain = time%rainrec

         ! Set first record of arrays
         atmo%raintimearray(1) = time%tcum + time%dtmin
         atmo%rainfluxarray(1) = 0.0d0

         ! rainfluxes (cm/d): Flx(t1) = P(t1) / (T(t1)−T(t0)) over interval T(t0)→T(t1)
         do i = 2, atmo%nmrain
            atmo%rainfluxarray(i) = rainam(i) / (atmo%raintimearray(i) - atmo%raintimearray(i - 1))
         end do

         ! Set final values of raintimearray and corresponding rainfluxarray
         atmo%raintimearray(atmo%nmrain + 1) = max(time%tend + 1.1d0 - time%tstart, &
                                                   atmo%raintimearray(atmo%nmrain) + 1.d0)
         atmo%rainfluxarray(atmo%nmrain + 1) = 0.d0
         atmo%nmrain = atmo%nmrain + 1

      end if

      ! For swrain = 1-3: determine start rain record
      time%rainrec = 1

      end associate
   end subroutine ProcessRainEvents

   !> Update meteorological fluxes for current time step
  !!
  !! Handles two types of sub-daily meteorological processing:
  !! 1. **Rainfall intensity mode** (flrainintens=true): Updates precipitation fluxes
  !!    based on rain event time series
  !! 2. **Detailed meteorology mode** (flmetdetail=true): Updates all meteorological
  !!    variables including ET and precipitation from detailed input records
  !!
  !! ## Rainfall Intensity Processing (swrain 1-3)
  !! - Retrieves rainfall flux from event arrays for current time step
  !! - Applies snow/rain partitioning (fprecnosnow)
  !! - Applies interception correction (finterception)
  !! - Calculates time until next rain event (dtEventRain)
  !!
  !! ## Detailed Meteorology Processing (swmetdetail=1)
  !! - Updates to next meteorological record when needed (flUpdMetDet=true)
  !! - Sets potential transpiration (ptra) and evaporation (peva)
  !! - Sets gross and net rainfall rates (graidt, nraidt)
  !! - Calculates actual soil evaporation considering surface wetness
  !!
  !! The routine is called at every time step when sub-daily processing is active.
  !!
  !! @note
  !! **Original documentation:**
  !! Last modified: February 2014
  !! Purpose: Calculations of meteo variables on time step basis
  !! (in case of precipitation intensities [swrain 1-3] or detailed meteo input)
  !! @endnote
   subroutine ProcessMeteoTsteps(state)
      use et_mod, only: reduceva_dt
      implicit none

      type(swap_state_t), intent(inout) :: state

      associate (atmo => state%atmosphere, time => state%timecontrol)

      if (time%flrainintens) then
         ! Per timestep: set precipitation fluxes from the rain-event array
         atmo%graidt  = atmo%fprecnosnow * atmo%rainfluxarray(time%rainrec)
         atmo%nraidt  = atmo%finterception * atmo%graidt
         atmo%aintcdt = (1.0d0 - atmo%finterception) * atmo%graidt

         ! Minimum time step length until next rain event
         ! (tcum + dt = time at end of current timestep)
         atmo%dtEventRain = atmo%raintimearray(time%rainrec) - (time%tcum + time%dt)

      else if (time%flmetdetail) then

         if (time%flUpdMetDet) then
            ! Per meteo time interval: update actual meteo record and set fluxes
            time%wrecord  = time%wrecord + 1
            atmo%ptra     = atmo%tpot(time%wrecord)
            atmo%peva     = atmo%epot(time%wrecord)
            atmo%graidt   = atmo%grain(time%wrecord)
            atmo%nraidt   = atmo%nrain(time%wrecord)
            atmo%aintcdt  = atmo%graidt - atmo%nraidt

            time%flUpdMetDet = .false.
         end if

         ! Per timestep: actual soil-evaporation rate
         call reduceva_dt(atmo%nraida, state)

      end if

      end associate
   end subroutine ProcessMeteoTsteps

   !> Distribute potential evapotranspiration according to diurnal sine wave
  !!
  !! Calculates time-varying potential transpiration and evaporation rates that follow
  !! a sine wave pattern during daylight hours. This provides a more realistic representation
  !! of diurnal ET dynamics compared to constant daily rates.
  !!
  !! ## Mathematical Approach
  !! The sine wave is centered at solar noon and extends over the photoperiodic daylength:
  !! - Before sunrise: ET = 0
  !! - During daylight: ET follows sine curve (maximum at noon)
  !! - After sunset: ET = 0
  !!
  !! The sine wave integral over the day equals the daily total (ptraday, pevaday),
  !! ensuring mass conservation.
  !!
  !! ## Implementation Details
  !! At day start (fldaystart=true):
  !! - Calculates photoperiodic daylength (daylp) using astronomical relationships
  !! - Determines sunrise time (tsunrise_atm) and sunset time (tsunset_atm)
  !!
  !! At each time step:
  !! - Determines fraction of daily ET for current time interval
  !! - Accounts for time steps spanning sunrise or sunset transitions
  !! - Calculates instantaneous rates: ptra, peva (cm/d)
  !! - Updates actual soil evaporation considering surface moisture
  !!
  !! ## State Variables
  !! - `tsunrise_atm`: Time of sunrise (fraction of day, 0-1)
  !! - `tsunset_atm`: Time of sunset (fraction of day, 0-1)
  !!
  !! Both are module-level variables in variables.f90
  !!
  !! @note
  !! **Original documentation:**
  !! Last modified: October 2008
  !! Purpose: Distributes potential transpiration and evaporation according to
  !! sine wave during photoperiodic daylight
  !! Note: tsunrise_atm and tsunset_atm are now module-level in variables.f90
  !! @endnote
   subroutine ETSine(state)
      ! DEFERRED — astro-derived scalars (rad, daylp, difpp, atmtr, dsinbe)
      ! plus tsunrise_atm/tsunset_atm remain bare globals. Migration belongs
      ! with the astro-caching arc when/if it's prioritized; per-day perf
      ! impact is negligible.
      use variables, only: rad, daylp, difpp, atmtr, dsinbe, tsunrise_atm, tsunset_atm
      use et_mod,    only: reduceva_dt
      implicit none

      type(swap_state_t), intent(inout) :: state

      real(8) :: daytime, dayl, sinld, cosld, fraction
      real(8), parameter :: pi = 3.14159265d0

      associate (atmo => state%atmosphere, time => state%timecontrol)

      if (time%flDayStart) then
         ! Photoperiodic daylength + sunrise/sunset (computed once/day)
         call astro(time%daynr, state%cfg%meteo%lat, rad, dayl, daylp, sinld, cosld, difpp, atmtr, dsinbe)
         tsunrise_atm = 0.5d0 - daylp/48.d0
         tsunset_atm  = 0.5d0 + daylp/48.d0
      end if

      ! Time as fraction of the day
      daytime = time%t1900 + time%dt - int(time%t1900)

      ! Fraction of daily ET in this timestep (sine-wave over the photoperiod)
      if (daytime .lt. tsunrise_atm) then
         fraction = 0.d0
      elseif (daytime .gt. tsunrise_atm .and. (daytime - time%dt) .lt. tsunrise_atm) then
         fraction = 0.5d0*(cos(pi/2.d0 + (tsunrise_atm - 0.5d0)/(tsunset_atm - tsunrise_atm)*pi) &
                         - cos(pi/2.d0 + (daytime      - 0.5d0)/(tsunset_atm - tsunrise_atm)*pi))
      elseif ((daytime - time%dt) .gt. tsunrise_atm .and. daytime .lt. tsunset_atm) then
         fraction = 0.5d0*(cos(pi/2.d0 + (daytime - time%dt - 0.5d0)/(tsunset_atm - tsunrise_atm)*pi) &
                         - cos(pi/2.d0 + (daytime           - 0.5d0)/(tsunset_atm - tsunrise_atm)*pi))
      elseif (daytime .gt. tsunset_atm .and. (daytime - time%dt) .lt. tsunset_atm) then
         fraction = 0.5d0*(cos(pi/2.d0 + (daytime - time%dt - 0.5d0)/(tsunset_atm - tsunrise_atm)*pi) &
                         - cos(pi/2.d0 + (tsunset_atm       - 0.5d0)/(tsunset_atm - tsunrise_atm)*pi))
      else
         fraction = 0.d0
      end if

      ! Set E and T fluxes
      atmo%peva = atmo%pevaday * fraction / time%dt
      atmo%ptra = atmo%ptraday * fraction / time%dt

      ! Actual soil-evaporation rate of current moment
      call reduceva_dt(atmo%nraida, state)

      end associate
   end subroutine ETSine

end module meteodt_mod
