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
      ! [SS-TC TC-14] yearmeteo,timjan1,rainrec read/written via state%timecontrol (ADR 0041)
      ! GR-ATM C7: wet(i) migrated → state%atmosphere%wet(i).
      ! [GR-CROP Phase B] nmrain/rainamount/rainfluxarray/raintimearray/arai migrated →
      !   state%atmosphere%X via associate aliases.  swrain/raintab remain narrow
      !   use variables until config arg is threaded (Arc 8).
      use array_utils, only: afgen
      use swap_array_dimensions, only: mrain
      implicit none

      type(swap_state_t), intent(inout) :: state
        !! Simulation state (for TC reader cutover — tcum via state%timecontrol)

      ! --- local
      integer i, iendyear, j, l, nlack, nn, rday, rdaya(367), rdayold
      real(8) araihlp(367), day(mrain), rainam(mrain), rainflux
      real(8) raintime, ratimar(mrain), tendyear, vsmall, wght, wwet(368)
      vsmall = 1.0d-8

      ! [SS-TC TC-14] alias TC fields directly so bare names below resolve to state%timecontrol
      ! [SS-BMI2 Task 4] tstart, tend, dtmin added to associate — retire globals in Task 5
      ! [GR-CROP Phase B] rain timing arrays aliased to state%atmosphere%X
      associate( tc_tcum => state%timecontrol%tcum, &
                 yearmeteo => state%timecontrol%yearmeteo, &
                 timjan1 => state%timecontrol%timjan1, &
                 rainrec => state%timecontrol%rainrec, &
                 tstart  => state%timecontrol%tstart, &
                 tend    => state%timecontrol%tend, &
                 dtmin   => state%timecontrol%dtmin, &
                 nmrain        => state%atmosphere%nmrain,        &  ! [GR-CROP Phase B]
                 rainamount    => state%atmosphere%rainamount,    &  ! [GR-CROP Phase B]
                 rainfluxarray => state%atmosphere%rainfluxarray, &  ! [GR-CROP Phase B]
                 raintimearray => state%atmosphere%raintimearray, &  ! [GR-CROP Phase B]
                 arai          => state%atmosphere%arai           )  ! [GR-CROP Phase B]

      ! === Process rain events on yearly basis ===

      ! For rain options 1 and 2: convert daily rain quantities and intensities or durations
      ! into rain events by creating raintime and rainflux arrays conform rain option 3
      if (state%cfg%meteo%swrain .eq. 1 .or. state%cfg%meteo%swrain .eq. 2) then
         rainrec = 1
         do i = 1, nmrain
            if (raintimearray(i + 1) .gt. tstart - vsmall) then
               rainrec = rainrec + 1
               ! Beginning (00:00) of days of current year within simulation period
               day(rainrec) = raintimearray(i + 1) - timjan1 + 1.d0
               rainam(rainrec) = 0.1d0*rainamount(i)  ! convert from mm to cm
               wwet(rainrec) = state%atmosphere%wet(i)   ! GR-ATM C7: wet→state%atmosphere%wet
            end if
         end do

         ! Set first record of raintime and rainflux (= 0)
         raintimearray(1) = tc_tcum + dtmin
         rainam(1) = 0.d0
         rainfluxarray(1) = 0.d0

         ! Set rest of records of raintime and rainflux (only when rainam[ount] > 0)
         nmrain = rainrec + 1
         rainrec = 0
         do i = 2, nmrain
            if (rainam(i) .gt. vsmall) then
               if (state%cfg%meteo%swrain .eq. 1) then
                  ! Mean rainfall intensities are specified
                  rainflux = afgen(state%cfg%meteo%raintab, 60, day(i))
                  raintime = min(0.99d0, rainam(i)/rainflux)

               elseif (state%cfg%meteo%swrain .eq. 2) then
                  ! Rainfall durations are specified
                  raintime = wwet(i)
               end if

               if (i .eq. 2) then
                  rainrec = rainrec + 1
               else
                  ! First raintime of a day: closure of last period of former day with rain = 0
                  rainrec = rainrec + 2
                  raintimearray(rainrec) = real(i - 2, real64) + tc_tcum
                  rainfluxarray(rainrec) = 0.d0
               end if
               ! Second raintime of a day: closure of first period of the day, rain = rainam
               raintimearray(rainrec + 1) = real(i - 2, real64) + tc_tcum + raintime
               rainfluxarray(rainrec + 1) = rainam(i)/raintime

            end if
         end do

         ! Extend array with records at end of current year
         tendyear = 365.d0
         if (mod(yearmeteo, 4) .eq. 0) tendyear = 366.d0
         raintimearray(rainrec + 2) = tc_tcum + tendyear + dtmin
         rainfluxarray(rainrec + 2) = 0.d0

         ! In case of rain events: 1) calculate daily values, 2) fill raintimearray and rainfluxarray
      elseif (state%cfg%meteo%swrain .eq. 3) then

         ! Total amount of rain per meteo day arai
         ! Initialize array with sum of rain
         do i = 1, 366
            araihlp(i) = 0.d0
         end do

         ! Less rain days than meteo days? Fill gap with dummies
         rdayold = int(raintimearray(1) - timjan1) + 1  ! first day with rain record of the year
         nlack = rdayold - 1
         do j = 1, nlack
            rdaya(j) = j
            araihlp(j) = 0.d0
         end do
         rdaya(j) = rdayold
         araihlp(j) = 0.d0

         ! Fill array of daily sums of rain with real values
         do i = 1, nmrain
            rday = int(raintimearray(i) - timjan1) + 1
            if (rday .gt. rdayold) then
               do l = 1, rday - rdayold - 1
                  j = j + 1
                  rdaya(j) = rdaya(j - 1) + 1
                  araihlp(j) = 0.d0
               end do
               ! In case of rain event exceeding current day, calculate weights for assigning parts to current and next day
               wght = (1.d0 - (raintimearray(i - 1) - real(int(raintimearray(i - 1)), real64)))/ &
                      (raintimearray(i) - raintimearray(i - 1))
               araihlp(j) = araihlp(j) + rainamount(i)*wght
               rdayold = rday
               j = j + 1
               rdaya(j) = rday
               araihlp(j) = araihlp(j) + rainamount(i)*(1.d0 - wght)
            else
               if (i .gt. 1) then
                  araihlp(j) = araihlp(j) + rainamount(i)
               end if
            end if
         end do

         ! Rain days missing at the end of the year? Fill gap with dummies
         iendyear = 365
         if (mod(yearmeteo, 4) .eq. 0) iendyear = 366
         nlack = iendyear - j
         do i = 1, nlack
            rdaya(j + i) = j + i
            araihlp(j + i) = 0.d0
         end do

         ! Save help array araihlp into arai
         do i = 1, 366
            arai(i) = araihlp(i)
         end do

         ! Assign values to raintimearray and rainam array for calculating rainfluxarray
         ! Find time gap without rain events at the beginning of the year
         rainrec = 1
         ratimar(1) = raintimearray(1) - tstart
         i = 1
         do while (ratimar(i) .lt. vsmall)
            i = i + 1
            ratimar(i) = raintimearray(i) - tstart
         end do
         do while (rainamount(i) .lt. vsmall .and. rainamount(i + 1) .lt. vsmall)
            i = i + 1
            ratimar(i) = raintimearray(i) - tstart
         end do

         ! Fill arrays with real values of rain events
         nn = i
         do i = nn, nmrain
            ratimar(i + 1) = raintimearray(i + 1) - tstart
            rainrec = rainrec + 1
            raintimearray(rainrec) = ratimar(i)
            rainam(rainrec) = 0.1d0*rainamount(i)  ! convert from mm to cm
         end do
         nmrain = rainrec

         ! Set first record of arrays
         raintimearray(1) = tc_tcum + dtmin
         rainfluxarray(1) = 0.0d0

         ! Calculate rainfluxes (cm/d) and fill rainfluxarray
         ! Flx(t1) = P(t1) / (T(t1)-T(t0))  = counts for time interval T(t0) -> T(t1)
         ! Flx = flux, P = quantity of rain, T = time
         do i = 2, nmrain
            rainfluxarray(i) = rainam(i)/(raintimearray(i) - raintimearray(i - 1))
         end do

         ! Set final values of raintimearray and corresponding rainfluxarray
         raintimearray(nmrain + 1) = dmax1(tend + 1.1d0 - tstart, raintimearray(nmrain) + 1.d0)
         rainfluxarray(nmrain + 1) = 0.d0
         nmrain = nmrain + 1

      end if

      ! For swrain = 1-3: determine start rain record
      rainrec = 1

      ! [GR-CROP Phase B] rain timing arrays written directly via state%atmosphere aliases —
      ! no bulk mirror needed (state IS the target).

      end associate  ! tc_tcum => state%timecontrol [TC-9] + state%atmosphere rain timing [GR-CROP Phase B]

      return
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
      ! SS-TC TC-9: flrainintens,tcum,dt,flmetdetail,flUpdMetDet removed from bare use variables;
      !             reads/writes via state%timecontrol.
      ! [SS-GR-ATM B28] bare use variables replaced with narrow only: list.
      !   Migrated reads: tpot/epot/grain/nrain → state%atmosphere%X (pure reads, no dual-write needed).
      ! [GR-CROP Phase B] rainfluxarray/raintimearray migrated → state%atmosphere%X.
      !   DEFERRED: finterception, dtEventRain — not yet in state schema (Arc 8).
      use variables, only: &   ! [SS-GR-FINAL B11] DEFERRED
         ! DEFERRED: finterception — fraction interception config; Phase C3
         finterception, &
         ! DEFERRED: dtEventRain — rain event timestep runtime state; Phase C3
         dtEventRain  ! [GR-CROP Phase B]
      use et_mod, only: reduceva_dt
      implicit none

      type(swap_state_t), intent(inout) :: state
        !! Simulation state (passed through to reduceva for atmosphere dual-writes)

      ! [SS-TC TC-14] alias TC fields directly so bare names below resolve to state%timecontrol
      ! [GR-CROP Phase B] rainfluxarray, raintimearray aliased to state%atmosphere%X
      associate( &
        tc_flrainintens => state%timecontrol%flrainintens,  &
        tc_tcum         => state%timecontrol%tcum,          &
        tc_dt           => state%timecontrol%dt,            &
        tc_flmetdetail  => state%timecontrol%flmetdetail,   &
        tc_flUpdMetDet  => state%timecontrol%flUpdMetDet,   &
        rainrec         => state%timecontrol%rainrec,       &
        wrecord         => state%timecontrol%wrecord,       &
        rainfluxarray   => state%atmosphere%rainfluxarray,  &  ! [GR-CROP Phase B]
        raintimearray   => state%atmosphere%raintimearray   )  ! [GR-CROP Phase B]

      ! === Precipitation intensities ===

      if (tc_flrainintens) then
         ! Per time step: set precipitation fluxes for current time step
         state%atmosphere%graidt  = state%atmosphere%fprecnosnow*rainfluxarray(rainrec)
         state%atmosphere%nraidt  = finterception*state%atmosphere%graidt
         state%atmosphere%aintcdt = (1.d0 - finterception)*state%atmosphere%graidt

         ! Calculate minimum time step length for occurrence of next rain event
         ! (tcum + dt = time at end of current timestep)
         dtEventRain = raintimearray(rainrec) - (tc_tcum + tc_dt)

         ! === Detailed meteo ===

      elseif (tc_flmetdetail) then

         if (tc_flUpdMetDet) then
            ! Per meteo time interval: update actual meteo record and set fluxes
            ! for current time of detailed meteo input
            wrecord = wrecord + 1
            state%atmosphere%ptra = state%atmosphere%tpot(wrecord)   ! [SS-GR-ATM B28] tpot → state%atmosphere%tpot
            state%atmosphere%peva = state%atmosphere%epot(wrecord)   ! [SS-GR-ATM B28] epot → state%atmosphere%epot
            state%atmosphere%graidt  = state%atmosphere%grain(wrecord)   ! [SS-GR-ATM B28] grain → state%atmosphere%grain
            state%atmosphere%nraidt  = state%atmosphere%nrain(wrecord)   ! [SS-GR-ATM B28] nrain → state%atmosphere%nrain
            state%atmosphere%aintcdt = state%atmosphere%graidt - state%atmosphere%nraidt

            tc_flUpdMetDet = .false.   ! [SS-TC TC-14] legacy flUpdMetDet write retired
         end if

         ! Per time step: calculate soil evaporation rate of current time step
         call reduceva_dt(state%atmosphere%nraida, state)

      end if

      end associate  ! tc_flrainintens,...,raintimearray,rainfluxarray => state%timecontrol/atmosphere [TC-9, GR-CROP Phase B]

      return
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
      ! SS-TC TC-9: fldaystart,daynr,t1900,dt removed from bare use variables;
      !             reads via state%timecontrol.
      ! [SS-GR-ATM B28] bare use variables replaced with narrow only: list.
      !   DEFERRED: lat → config%meteo%lat (requires config arg, Arc 8 candidate);
      !             rad, daylp, difpp, atmtr, dsinbe, tsunrise_atm, tsunset_atm —
      !             not yet in state/config schema (Arc 8+).
      use variables, only: &   ! [SS-GR-FINAL B11] DEFERRED
         ! DEFERRED: lat — site latitude; config%meteo%lat; Phase C3
         ! DEFERRED: rad/daylp/difpp/atmtr/dsinbe/tsunrise_atm/tsunset_atm — meteo derived scalars; Phase C3
         rad, daylp, difpp, atmtr, dsinbe, tsunrise_atm, tsunset_atm
      use et_mod, only: reduceva_dt
      implicit none

      type(swap_state_t), intent(inout) :: state
        !! Simulation state (passed through to reduceva for atmosphere dual-writes)

      ! --- local
      real(8) daytime, pi, dayl, sinld, cosld, fraction
      data pi/3.14159265d0/    ! number pi [-]

      ! SS-TC TC-9: tc_* aliases for fldaystart, daynr, t1900, dt.
      associate( &
        tc_fldaystart => state%timecontrol%flDayStart,  &  ! TC-9
        tc_daynr      => state%timecontrol%daynr,       &  ! TC-9
        tc_t1900      => state%timecontrol%t1900,       &  ! TC-9
        tc_dt         => state%timecontrol%dt           )  ! TC-9

      if (tc_fldaystart) then
         ! Determine duration photoperiodic daylight in hours
         call astro(tc_daynr, state%cfg%meteo%lat, rad, dayl, daylp, sinld, cosld, difpp, atmtr, dsinbe)
         ! Determine tsunrise_atm, tsunset_atm and daytime
         tsunrise_atm = 0.5d0 - daylp/48.d0
         tsunset_atm = 0.5d0 + daylp/48.d0
      end if

      ! Set time as fraction of the day
      daytime = tc_t1900 + tc_dt - int(tc_t1900)

      ! Determine fraction of fluxes according to sine wave during this time step
      if (daytime .lt. tsunrise_atm) then
         fraction = 0.d0
      elseif (daytime .gt. tsunrise_atm .and. (daytime - tc_dt) .lt. tsunrise_atm) then
         fraction = 0.5d0*(dcos(pi/2.d0 + (tsunrise_atm - 0.5d0)/ &
                                (tsunset_atm - tsunrise_atm)*pi) - dcos(pi/2.d0 + (daytime - 0.5d0)/ &
                                                                        (tsunset_atm - tsunrise_atm)*pi))
      elseif ((daytime - tc_dt) .gt. tsunrise_atm .and. (daytime) .lt. tsunset_atm) then
         fraction = 0.5d0*(dcos(pi/2.d0 + (daytime - tc_dt - 0.5d0)/ &
                                (tsunset_atm - tsunrise_atm)*pi) - dcos(pi/2.d0 + (daytime - 0.5d0)/ &
                                                                        (tsunset_atm - tsunrise_atm)*pi))
      elseif (daytime .gt. tsunset_atm .and. (daytime - tc_dt) .lt. tsunset_atm) then
         fraction = 0.5d0*(dcos(pi/2.d0 + (daytime - tc_dt - 0.5d0)/ &
                                (tsunset_atm - tsunrise_atm)*pi) - dcos(pi/2.d0 + (tsunset_atm - 0.5d0)/ &
                                                                        (tsunset_atm - tsunrise_atm)*pi))
      else
         fraction = 0.d0
      end if

      ! Set E and T fluxes
      state%atmosphere%peva = state%atmosphere%pevaday*fraction/tc_dt
      state%atmosphere%ptra = state%atmosphere%ptraday*fraction/tc_dt

      ! Actual soil evaporation rate of current moment
      ! SS-ATM A-2.6: nraida retired — read from state%atmosphere%nraida
      call reduceva_dt(state%atmosphere%nraida, state)

      end associate  ! tc_fldaystart, tc_daynr, tc_t1900, tc_dt => state%timecontrol [TC-9]

      return
   end subroutine ETSine

end module meteodt_mod
