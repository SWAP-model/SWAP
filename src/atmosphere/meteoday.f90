!> Module for meteorological data input/output operations
!!
!! This module handles reading and managing meteorological forcing data for the SWAP model.
!! It supports two modes of meteorological input:
!! - **Daily meteorology** (swmetdetail=0): Single daily values for radiation, temperature,
!!   humidity, wind speed, precipitation, and reference ET
!! - **Detailed meteorology** (swmetdetail=1): Sub-daily timestep data for high-resolution
!!   simulations
!!
!! The module also manages precipitation partitioning into rain and snow based on temperature
!! thresholds and maintains cumulative meteorological flux accounting.
!!
!! @author Original SWAP development team
!! @date Last modified March 2014, refactored February 2026
module meteo_process_mod
   use error_mod, only: fatalerr_collected
   use swap_state_mod, only: swap_state_t
   use swap_config_mod, only: swap_config_t   ! [SS-GR-ATM B23] config added for meteo switches
   implicit none
   private

   public :: ReadMeteoDay
   public :: ResetMetFlx

contains

  !> Read meteorological data for the current simulation day
  !!
  !! Retrieves meteorological forcing data from pre-loaded arrays and performs
  !! initial processing including:
  !! - Data availability checking and error handling
  !! - Unit conversions (mm to cm for precipitation, period averages for detailed meteo)
  !! - Temperature-based calculations (24h average, daytime average)
  !! - Vapor pressure and relative humidity computations
  !! - Rain/snow partitioning based on temperature thresholds
  !!
  !! ## Operating Modes
  !!
  !! ### Daily Meteo (swmetdetail=0)
  !! - Reads single daily values for: radiation, min/max temperature, humidity, wind, 
  !!   precipitation, reference ET
  !! - Calculates saturated vapor pressure using Tetens equation
  !! - Partitions precipitation into rain/snow using transition temperatures
  !!
  !! ### Detailed Meteo (swmetdetail=1)
  !! - Reads sub-daily records (typically hourly or 3-hourly)
  !! - Validates record numbers and timestamps
  !! - Converts precipitation from mm/period to cm
  !! - Snow calculations disabled for detailed mode
  !!
  !! ## Temperature-Based Snow Partitioning
  !! When snow calculations are enabled (swsnow=1):
  !! - \( T_{av} > T_{rain} \): All precipitation as rain
  !! - \( T_{av} < T_{snow} \): All precipitation as snow
  !! - \( T_{snow} < T_{av} < T_{rain} \): Linear interpolation between rain and snow
  !!
  !! @warning Requires meteorological data to be pre-loaded into module variables
  !! @warning Simulation must start within the time range of available meteo data
  !!
  !! @note
  !! **Original documentation:**
  !! Last modified: March 2014
  !! Purpose: Returns meteorological fluxes of current day or of parts of a day 
  !! (detailed meteo input)
  !! @endnote
  subroutine ReadMeteoDay(state, config)
      ! [SS-GR-ATM B23] use variables dropped; symbols → state%atmosphere/config%meteo
      ! SS-TC TC-9: date,t1900 removed from only-list; reads/writes via state%timecontrol.
      ! [SS-TC TC-14] yearmeteo, daymeteo retired — read via state%timecontrol
      ! [SS-GR-ATM B23] DEFERRED symbols (not yet in state); [SS-GR-FINAL B11] reviewed
      use variables, only: &
          rad, tmn, tmx,                                  &  ! B23 DEFERRED — daily scalars
          tav,                                            &  ! B23 DEFERRED — dual-write (consumed by snow.f90/swapoutput.f90; tavd/rh dropped B.5)
          pathatm, metfil,                                &  ! B23 DEFERRED — filename strings
          detrecord, dettime, detrad, dethum, dettav,     &  ! B23 DEFERRED — detail arrays
          detrain, detwind, irectotal                        ! B23 DEFERRED — detail arrays + counter
      ! [SS-GR-ATM B.5] tavd retired from import: no legacy consumers remain after cropgrowth migration
      ! [SS-GR-ATM B.5] rh retired from import: no legacy consumers outside init seed (swap_mod A12)
      ! [SS-GR-ATM B.5] out_rad/tmn/tmx/hum/win/etr/wet retired from import: no consumers outside swap_mod init seed
      use MeteoVars
      use precipitation_mod, only: PartitionPrecipitation
      implicit none

      type(swap_state_t),   intent(inout) :: state   !! [SS-ATM] threaded for atmosphere dual-writes
      type(swap_config_t),  intent(in)    :: config  !! [SS-GR-ATM B23] for meteo config switches

    ! --- local
    ! [SS-ATM A-2.6] grai/gsnow/snrai/ssnow/fprecnosnow are transitional locals fed to
    ! PartitionPrecipitation; state%atmosphere%X is the canonical write target.
    real(8) :: grai, gsnow, snrai, ssnow, fprecnosnow
    character(len=11)  detdate
    character(len=3)   ext
    character(len=200) filnam
    character(len=300) messag

    ! ----------------------------------------------------------------------

    ! [SS-TC TC-14] alias TC fields directly so bare names resolve to state%timecontrol
    associate( &
      tc_t1900 => state%timecontrol%t1900, &
      tc_date  => state%timecontrol%date,  &
      yearmeteo => state%timecontrol%yearmeteo, &
      daymeteo  => state%timecontrol%daymeteo )

    call ResetMetFlx (state)

    ! 1: Check whether meteo data are available of today; pass on weather of today
    ! 1.0 Daily Meteo 0000000000000000000000000000000000000000000000000000000 Daily Meteo

    if (config%meteo%swmetdetail.eq.0) then

      ! Check availability of meteo data of today
      if (daymeteo.lt.state%atmosphere%daynrfirst .or. daymeteo.gt.state%atmosphere%daynrlast) then
        messag ='In meteo file no meteo data are'// &
                ' available for '//tc_date//'. First adapt meteo file!'
        call fatalerr_collected ('meteo',messag)
      end if

      ! Pass on weather values of today
      rad  = state%atmosphere%arad(daymeteo+1-state%atmosphere%daynrfirst)
      tmn  = state%atmosphere%atmn(daymeteo+1-state%atmosphere%daynrfirst)
      tmx  = state%atmosphere%atmx(daymeteo+1-state%atmosphere%daynrfirst)
      hum  = state%atmosphere%ahum(daymeteo+1-state%atmosphere%daynrfirst)
      win  = state%atmosphere%awin(daymeteo+1-state%atmosphere%daynrfirst)
      grai = state%atmosphere%arai(daymeteo+1-state%atmosphere%daynrfirst)  ! [SS-ATM A-2.6] local; PartitionPrecipitation converts mm->cm and writes state
      etr  = state%atmosphere%aetr(daymeteo+1-state%atmosphere%daynrfirst)

      ! If hum is missing or tav cannot be calculated: set rh at -99.0
      state%atmosphere%rh = 1.0d0
      if (hum.lt.-98.0d0 .or. tmn.lt.-98.0d0 .or. tmx.lt.-98.0d0) &
        state%atmosphere%rh=-99.0d0

      ! Calculate 24h average temperature
      state%atmosphere%Tav = (tmx+tmn)*0.5d0
      tav = state%atmosphere%Tav   ! [SS-GR-ATM B23] dual-write — legacy tav consumed by snow.f90/swapoutput.f90
      ! Calculate average day temperature
      state%atmosphere%tavd = (tmx+state%atmosphere%Tav)*0.5d0   ! [SS-GR-ATM B23] direct state write
      ! [SS-GR-ATM B.5] tavd dual-write to legacy global RETIRED: no crop consumers remain

      if (state%atmosphere%rh.ge.-98.0d0) then
        ! Calculate saturated vapour pressure [kpa]
        svp = 0.3055d0*(dexp(17.27d0*tmn/(tmn+237.3d0)) + &
                        dexp(17.27d0*tmx/(tmx+237.3d0)))
        ! Calculate relative humidity [fraction]
        state%atmosphere%rh = min(hum/svp,1.0d0)   ! [SS-GR-ATM B23] direct state write
      endif
      ! [SS-GR-ATM B.5] rh dual-write to legacy global RETIRED: no crop consumers remain

      ! CFO file for PEARL: save meteo variables of today for output
      state%atmosphere%out_rad = real(rad, kind=8)          ! [SS-GR-ATM B23] direct state write
      state%atmosphere%out_tmn = real(tmn, kind=8)          ! [SS-GR-ATM B23] direct state write
      state%atmosphere%out_tmx = real(tmx, kind=8)          ! [SS-GR-ATM B23] direct state write
      state%atmosphere%out_hum = real(hum, kind=8)          ! [SS-GR-ATM B23] direct state write
      state%atmosphere%out_win = real(win, kind=8)          ! [SS-GR-ATM B23] direct state write
      state%atmosphere%out_etr = real(etr, kind=8)*0.001d0  ! [SS-GR-ATM B23] direct state write
      if (config%meteo%swrain.eq.2) then
        state%atmosphere%out_wet = real(state%atmosphere%wet(daymeteo+1-state%atmosphere%daynrfirst), kind=8)  ! [SS-GR-ATM B23]
      else
        state%atmosphere%out_wet = -1.0d0
      endif
      ! [SS-GR-ATM B.5] out_rad/tmn/tmx/hum/win/etr/wet legacy dual-writes RETIRED.
      ! state%atmosphere%out_* are now the sole write targets; legacy globals have no consumers.

    ! end 1 Daily Meteo 00000000000000000000000000000000000000000000000000000 Daily Meteo

    ! 1.1 Detailed Meteo 1111111111111111111111111111111111111111111111111111 Detailed Meteo

    elseif (config%meteo%swmetdetail.eq.1) then

      ! Check availability of meteo data of today
      ! Compose filename meteorological file for use in warnings
      write (ext,'(i3.3)') mod(yearmeteo,1000)
      filnam = trim(pathatm)//trim(metfil)//'.'//trim(ext)

      do i = 1, config%meteo%nmetdetail
        irectotal = irectotal + 1
        if (i .ne. detrecord(irectotal)) then
          messag='In meteo file '//trim(filnam)//' record number(s)'// &
                 ' are not correct at '//tc_date//'. First adapt meteo file!'
          call fatalerr_collected ('meteo',messag)
        end if
        call dtdpst('year-month-day', &
                    dettime(irectotal)+0.1d0,detdate)
        call dtdpst('year-month-day',tc_t1900+0.1d0,tc_date)
        if (detdate .ne. tc_date) then
          messag ='In meteo file '//trim(filnam)//' the amount of '// &
                  'records deviate near '//tc_date//'. First adapt meteo file!'
          call fatalerr_collected ('meteo',messag)
        end if

        ! Pass on weather records of today
        state%atmosphere%arad(i)  = detrad(irectotal)
        state%atmosphere%ahum(i)  = dethum(irectotal)
        state%atmosphere%atav(i)  = dettav(irectotal)
        awind(i) = detwind(irectotal)
        arain(i) = detrain(irectotal) * 0.1d0 ! convert from mm to cm
      enddo
    endif

    ! end 1 Detailed Meteo 11111111111111111111111111111111111111111111111111 Detailed Meteo
    ! end 1.

    ! Call precipitation partitioning module
    ! [SS-ATM A-2.6] ssnow seeded from state (intent inout for read in swmetdetail==0 path);
    ! state%atmosphere%ssnow is written directly by PartitionPrecipitation in swmetdetail==1.
    ssnow = state%atmosphere%ssnow
    call PartitionPrecipitation(config%meteo%swmetdetail, config%meteo%snow%swsnow, &
                                state%atmosphere%Tav, state%atmosphere%teprrain, &
                                state%atmosphere%teprsnow, &
                                ssnow, config%meteo%nmetdetail, arain, grai, gsnow, snrai, &
                                fprecnosnow, restint, state)

    end associate  ! tc_t1900, tc_date => state%timecontrol [TC-9]

    return
  end subroutine ReadMeteoDay


  !> Reset intermediate and cumulative meteorological flux counters
  !!
  !! Manages the zeroing of meteorological flux accumulators based on control flags.
  !! This allows for flexible reporting periods (e.g., daily, seasonal, annual totals).
  !!
  !! ## Reset Operations
  !!
  !! ### Intermediate fluxes (flzerointr = .true.)
  !! Typically reset at shorter intervals (e.g., daily):
  !! - `iprec`: Intermediate precipitation
  !! - `igrai`: Intermediate gross rainfall
  !! - `inrai`: Intermediate net rainfall
  !!
  !! ### Cumulative fluxes (flzerocumu = .true.)
  !! Typically reset at longer intervals (e.g., seasonal, annual):
  !! - `cgrai`: Cumulative gross rainfall
  !! - `cnrai`: Cumulative net rainfall  
  !! - `caintc`: Cumulative interception
  !!
  !! @note Control flags `flzerointr` and `flzerocumu` are set by the main 
  !! controller based on the reporting schedule
  !!
  !! @note
  !! **Original documentation:**
  !! Last modified: February 2014
  !! Purpose: Reset intermediate and cumulative meteorological fluxes
  !! Interface: I - flzerointr, flzerocumu, caintc, cgrai, cnrai, igrai, inrai, iprec
  !!            O - caintc, cgrai, cnrai, igrai, inrai, iprec
  !! @endnote
  subroutine ResetMetFlx (state)
      ! [SS-SWC S-2.12B] iprec retired — state%soilwater%reset_intermediate() handles it
      implicit none

      type(swap_state_t), intent(inout) :: state  !! [SS-ATM A-2.1] cohort reset() dispatch

    ! --- local

    ! Reset cumulative intermediate fluxes
    if (state%timecontrol%flZeroIntr) then
      ! [SS-ATM A-2.6] canonical reset for all 8 intr fields; legacy igrai/inrai retired
      call state%atmosphere%intr%reset()
    endif

    ! Reset cumulative meteorological fluxes
    if (state%timecontrol%flZeroCumu) then
      ! [SS-ATM A-2.6] canonical reset for all 10 cumu fields; legacy cgrai/cnrai/caintc retired
      call state%atmosphere%cumu%reset()
    endif

    return
  end subroutine ResetMetFlx

end module meteo_process_mod

!> Main meteorological processing coordinator module
!!
!! This module serves as the primary coordinator for daily meteorological processing
!! in the SWAP model. It orchestrates the calculation of evapotranspiration components,
!! interception, and partitioning between soil evaporation and transpiration.
!!
!! The module integrates functionality from several specialized modules:
!! - meteo_process_mod: Reading and resetting meteorological data
!! - interception_mod: Various interception calculation methods
!! - et_mod: Evapotranspiration calculations (Penman-Monteith)
!! - runoff_mod: Surface runoff calculations (SCS Curve Number method)
!!
!! ## Main Functionality
!! The ProcessMeteoDay subroutine handles:
!! - Daily and sub-daily meteorological data processing
!! - Interception calculations using multiple methods (Von Hoyningen-Hune & Braden, Gash, Rutter)
!! - Reference evapotranspiration via Penman-Monteith or user input
!! - Partitioning of atmospheric demand into soil evaporation and transpiration
!! - Wet fraction calculations for crop canopy
!! - Atmospheric CO2 corrections for transpiration
!!
!! @author Original SWAP development team
!! @date Last modified February 2014, refactored February 2026
module meteo_mod

  use meteo_process_mod, only: ReadMeteoDay, ResetMetFlx
  use interception_mod, only: VonHHBraden, Gash, ruttervw, DivIntercep
  use et_mod, only: PenMon, reduceva, pm_inputs_t, pm_outputs_t
  use runoff_mod, only: CNmethod
  use swap_state_mod, only: swap_state_t
  use swap_config_mod, only: swap_config_t   ! [SS-GR-ATM B24] config added for meteo switches

  implicit none

  private
  public :: ProcessMeteoDay

contains

  !> Process daily meteorological data and calculate evapotranspiration components
  !!
  !! Main orchestrator for daily meteorological processing. Performs:
  !! - Interception calculations (Von Hoyningen-Hune & Braden, Gash, or Rutter methods)
  !! - Reference evapotranspiration calculation or input processing
  !! - Potential soil evaporation and transpiration partitioning
  !! - Wet fraction calculations for crop canopy
  !! - Sub-daily timestep handling for detailed meteorology
  !!
  !! ## Processing Steps
  !! 1. Calculate interception and net rainfall/irrigation
  !! 2. Loop over day parts (1 for daily, nmetdetail for sub-daily)
  !! 3. Calculate ET0, EW0, ES0 via Penman-Monteith or use specified values
  !! 4. Apply interception for adapted Rutter model (if selected)
  !! 5. Calculate wet fraction of canopy
  !! 6. Partition into potential soil evaporation and transpiration
  !! 7. Apply atmospheric CO2 corrections if enabled
  !! 8. Aggregate results for detailed meteorology
  !!
  !! @note
  !! **Original documentation:**
  !! Last modified: February 2014
  !! Supports both daily (swmetdetail=0) and detailed (swmetdetail=1) meteorology
  !! Uses module variables from Variables and MeteoVars
  !! @endnote
  subroutine ProcessMeteoDay(state, config)
    ! [SS-GR-ATM B24] use variables partially retired; symbols → state%atmosphere/state%crop/config%meteo
    ! SS-TC TC-9: daynr,t,dt,flmetdetail,fletsine removed from only-list; reads via state%timecontrol.
    ! [SS-TC TC-14] metperiod retired — read via state%timecontrol
    ! [SS-GR-ATM B24] DEFERRED symbols (not yet in state/config); [SS-GR-FINAL B11] reviewed
    use variables, only: &
        ! B24 DEFERRED — crop factor scalar; cf/ch/cfeic retired
        rad,                                               &  ! B24 DEFERRED — daily radiation scalar
        logf,                                              &  ! B24 DEFERRED — Arc 9 (logging)
        angstroma, angstromb,              &  ! B24 DEFERRED — config ET params
        daylp, tmn, tmx, difpp,                            &  ! B24 DEFERRED — ET calculation params; albedo/rsc/rsw retired
        dsinbe, atmtr, rsoil,                              &  ! B24 DEFERRED — ET calculation params
        ! swinter retired — see state%crop%common%swinter
        croptype, gc,                                      &  ! B24 DEFERRED — crop/config fields; swgc retired
        ! [GR-CROP C3] icrop/flCropCalendar → state%crop%common%X
        flCropHarvest, cfevappond, flco2,                  &  ! B24 DEFERRED; fco2tra retired
        siccaptb,                                          &  ! B24 DEFERRED — interception table
        swusecn,                                           &  ! B24 DEFERRED — config switch (CN)
        finterception,                                     &  ! B24 DEFERRED — not yet in state
        tav                                                   ! B24 DEFERRED — dual-write (consumed by snow.f90/swapoutput.f90; tavd/rh dropped B.5)
     ! [SS-SWC S-2.12B] pond retired — read via state%soilwater%pond
    use swap_array_dimensions, only: magrs
    use MeteoVars
    use array_utils, only: afgen
    use runoff_mod, only: CNmethod
    use interception_mod, only: VonHHBraden, Gash, ruttervw, msw1eic, DivIntercep
    use swap_constants, only: nihil, small
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none

    type(swap_state_t),  intent(inout) :: state
      !! Simulation state (passed through to reduceva for atmosphere dual-writes)
    type(swap_config_t), intent(in)    :: config
      !! [SS-GR-ATM B24] for meteo config switches

    real(8)  rcs
    data     rcs/0.15d0/

    type(pm_inputs_t)  :: pmi
    type(pm_outputs_t) :: pmo

    ! SS-ATM Phase 1 Task A-1.8: ASSOCIATE aliases for atmosphere flat-scalar dual-writes
    ! SS-TC TC-9: tc_* aliases for daynr, t, dt, flmetdetail added.
    associate( &
       at_peva       => state%atmosphere%peva,        &
       at_ptra       => state%atmosphere%ptra,        &
       at_atmdem     => state%atmosphere%atmdem,      &
       at_grai       => state%atmosphere%grai,        &
       at_pevaday    => state%atmosphere%pevaday,     &
       at_ptraday    => state%atmosphere%ptraday,     &
       tc_daynr       => state%timecontrol%daynr,       &  ! TC-9
       tc_t           => state%timecontrol%t,           &  ! TC-9
       tc_dt          => state%timecontrol%dt,          &  ! TC-9
       tc_flmetdetail => state%timecontrol%flmetdetail, &  ! TC-9
       tc_fletsine    => state%timecontrol%fletsine,    &  ! TC-9
       metperiod      => state%timecontrol%metperiod,   &  ! [SS-TC TC-14]
       swscre         => state%timecontrol%swscre       )  ! [SS-BMI2 Task 4]

    ! === Section 3: Interception calculations ===

    ! Calculation of interception and net rain & net irrigation depth [cm]
    if ((state%crop%lai .lt. 1.d-3) .or. (state%atmosphere%grai+state%crop%gird .lt. 1.d-5) .or. &
        (state%crop%common%swinter.eq.0) .or. (state%atmosphere%gsnow.gt.0.0d0) .or.(state%atmosphere%ssnow.gt.0.0d0)) then

      ! No vegetation, rainfall/irrigation or interception calculation
      aintc = 0.d0

    else if (state%crop%common%swinter .eq. 1) then
      ! Calculate interception, method Von Hoyningen-Hune and Braden
      ! SS-ATM A-2.6: grai retired — pass state%atmosphere%grai explicitly
      ! SS-GR-ATM B8: state added for crop/atmosphere fields
      call VonHHBraden (aintc, state%atmosphere%grai, state)
    else if (state%crop%common%swinter .eq. 2) then
      ! Calculate interception, method Gash (1995)
      ! SS-ATM A-2.6: grai retired — pass state%atmosphere%grai explicitly
      ! SS-TC TC-11: state added for t via state%timecontrol%t
      call Gash (aintc, state%atmosphere%grai, state)
    end if

    ! Divide interception into rain part and irrigation part and
    ! calculate net rain (nraida) and net sprinkling irrigation (nird)
    if (state%crop%common%swinter.ne.3) &
      call DivIntercep (aintc, state)

    ! === LOOP over dayparts ===

    ! Set number of repetitions within one day
    if (config%meteo%swmetdetail.eq.0) then
      ! In case of daily meteo input: only one time per day
      ndayparts = 1
    elseif (config%meteo%swmetdetail.eq.1) then
      ndayparts = config%meteo%nmetdetail
    endif

    do 1000 irecord = 1, ndayparts

      ! === Section 4: Calculate evapotranspiration (et0, ew0, es0) ===
      ! Reference evapotranspiration has been specified
      if (config%meteo%swmetdetail.eq.0 .and. config%meteo%swetr.eq.1) then
        if (.not. state%crop%flCropEmergence) then
          ! no crop
          state%crop%et0 = 0.0d0
          state%crop%ew0 = 0.0d0
          state%crop%es0 = etr
          if (state%crop%swcfbs.eq.1) state%crop%es0 = state%crop%cfbs*etr
        else
          ! crop is present
          if (state%crop%swcf.eq.1 .or. state%crop%swcf.eq.3) then
            state%crop%et0 = state%crop%common%cf*etr
            if (state%crop%swcf .eq. 1) then
              state%crop%ew0 = state%crop%common%cf*etr
            else
              state%crop%ew0 = state%crop%fixed%cfeic*etr
            endif
          endif
          state%crop%es0 = etr
          if (state%crop%swcfbs.eq.1) state%crop%es0 = state%crop%cfbs*etr
        endif

      ! Reference evapotranspiration must be calculated
      elseif (config%meteo%swmetdetail.eq.1 .or. config%meteo%swetr.eq.0) then

        if (config%meteo%swmetdetail.eq.1) then
          ! Define weather variables of current record
          rad = state%atmosphere%arad(irecord) / metperiod     ! from j/m2/period to j/m2/d
          state%atmosphere%Tav = state%atmosphere%atav(irecord)
          hum = state%atmosphere%ahum(irecord)
          win = awind(irecord)
        endif

        ! Calculate evapotranspiration using Penman-Monteith: et0, ew0, es0 (mm/d)
        ! in case of daily meteo (swmetdetail = 0) irecord is always 1
        ! Pack PM inputs.
        pmi%daynr           = tc_daynr
        pmi%irecord         = irecord
        pmi%nmetdetail      = config%meteo%nmetdetail
        pmi%flmetdetail     = tc_flmetdetail
        pmi%flCropEmergence = state%crop%flCropEmergence
        pmi%swcf            = state%crop%swcf
        pmi%swdivide        = state%cfg%meteo%swdivide

        pmi%lat  = state%cfg%meteo%lat
        pmi%alt  = state%cfg%meteo%alt
        pmi%altw = state%cfg%meteo%altw
        pmi%a    = angstroma
        pmi%b    = angstromb
        pmi%rcs  = rcs

        pmi%rad    = rad
        pmi%tav    = state%atmosphere%Tav
        pmi%tmn    = tmn
        pmi%tmx    = tmx
        pmi%hum    = hum
        pmi%win    = win
        pmi%atmtr  = atmtr
        pmi%difpp  = difpp
        pmi%dsinbe = dsinbe
        pmi%daylp  = daylp

        pmi%rsc    = state%crop%common%rsc
        pmi%rsw    = state%crop%common%rsw
        pmi%ch     = state%crop%common%ch
        pmi%albedo = state%crop%common%albedo
        pmi%kdif   = state%crop%kdif
        pmi%kdir   = state%crop%kdir
        pmi%lai    = state%crop%lai

        pmi%rsoil  = rsoil

        call PenMon(pmi, pmo, logf, swscre)

        ! Unpack PM outputs to existing state / locals.
        state%crop%es0 = pmo%es0
        state%crop%et0 = pmo%et0
        state%crop%ew0 = pmo%ew0
        Edirect        = pmo%Edirect
        Tdirect        = pmo%Tdirect
        Tdirectwet     = pmo%Tdirectwet
        Edirectpond    = pmo%Edirectpond

        if (.not. state%crop%flCropEmergence) then
          ! no crop
          if (state%crop%swcfbs .eq. 1) then
            if (state%crop%swcf .eq. 1) then
              state%crop%es0 = state%crop%cfbs*state%crop%et0
            else
              state%crop%es0 = state%crop%cfbs*state%crop%es0
            endif
          endif
          state%crop%et0 = 0.0d0
          if (config%meteo%swmetdetail.eq.1 .and. (state%crop%swcf.eq.1 .or. state%crop%swcf.eq.3)) then
            if (state%crop%swcf.eq.1) then
              state%crop%ew0 = state%crop%common%cf*state%crop%ew0
            else
              state%crop%ew0 = state%crop%fixed%cfeic*state%crop%ew0
            endif
          endif
        else
          ! crop is present
          if (state%crop%swcfbs .eq. 1) then
            if (state%crop%swcf .eq. 1) then
              state%crop%es0 = state%crop%cfbs*state%crop%et0
            else
              state%crop%es0 = state%crop%cfbs*state%crop%es0
            endif
          endif
          if (state%crop%swcf.eq.1 .or. state%crop%swcf.eq.3) then
            state%crop%et0 = state%crop%common%cf*state%crop%et0
            if (state%crop%swcf.eq.1) then
              state%crop%ew0 = state%crop%common%cf*state%crop%ew0
            else
              state%crop%ew0 = state%crop%fixed%cfeic*state%crop%ew0
            endif
          endif
        endif

      endif
      ! [SS-GR-ATM B24] et0/ew0/es0 written directly to state%crop above

      ! === Section 5: Interception option NHI (adapted Rutter model) ===

      if (state%crop%common%swinter .eq. 3) then
        ! Set parameter values
        if (config%meteo%swmetdetail.eq.0) then  ! if swmetdetail = 0, siccapact is set in cropgrowth module
          if (croptype(state%crop%common%icrop).eq.1 .and. state%crop%common%swgc.eq.2) then   ! [GR-CROP C3]
            gctp  = gc
          else
            gctp  = 1.0d0 - dexp(-1.0d0*state%crop%kdir*state%crop%kdif*state%crop%lai)
            if (gctp .lt. 1.0d-5) then
              state%atmosphere%siccapact = 0.0d0   ! [SS-GR-ATM B24] direct state write
            endif
          endif
          dttp = 1.0d0  ! value of 1 d required for the daily meteo option
        elseif (config%meteo%swmetdetail.eq.1) then
          state%atmosphere%siccapact = afgen(siccaptb,(2*magrs),tc_t)   ! [SS-GR-ATM B24] direct state write
          if (croptype(state%crop%common%icrop).eq.1 .and. state%crop%common%swgc.eq.2) then
            gctp  = gc
          elseif (croptype(state%crop%common%icrop).eq.1 .and. state%crop%common%swgc.eq.1) then
            gctp  = 1.0d0 - dexp(-1.0d0*state%crop%kdir*state%crop%kdif*state%crop%lai)
          endif
          if (gctp .lt. 1.0d-5) then
            state%atmosphere%siccapact = 0.0d0   ! [SS-GR-ATM B24] direct state write
          endif
          dttp = tc_dt
        endif

        ! Calculate interception, method Rutter
        call ruttervw (gctp,aintc,eintc,state)

        ! Divide interception into rain part and irrigation part and
        ! calculate net rain (nraida) and net sprinkling irrigation (nird)
        call DivIntercep(aintc, state)
      endif

      ! === Section 6: Fraction of the day or period the crop is wet ===

      ! Calculate fraction of the day or period the crop is wet
      if (config%meteo%swmetdetail.eq.0) then
        ! Fraction of the day the crop is wet
        if (state%crop%ew0.lt.0.0001d0) then
          wfrac = 0.0d0
        else
          if (state%crop%ew0.lt.0.0001d0) then
            wfrac = 0.0d0
          else
            if (state%crop%common%swinter .ne. 3) then
              if (state%cfg%meteo%swdivide .eq. 0) then
                wfrac = max(min(aintc*10.0d0/state%crop%ew0,1.0d0),0.0d0)
              else
                if(tdirectwet.gt.nihil) then
                  wfrac = max(min(aintc*10.0d0/tdirectwet,1.0d0),0.0d0)
                else
                  wfrac = 0.0d0
                endif
              endif
            else
              wfrac = max(min(eintc*10.0d0/state%crop%ew0,1.0d0),0.0d0)
            endif
          endif
        endif
      ! Fraction of the period the crop is wet
      elseif (config%meteo%swmetdetail.eq.1) then
        if (state%atmosphere%grai .lt. 1.0d-12) then
          interc = 0.0d0
          wfrac  = 0.0d0
        else
          interc = restint + aintc * arain(irecord) / state%atmosphere%grai
          if (state%crop%ew0.lt.0.0001d0) then
            wfrac = 0.0d0
          else
            if (state%crop%swcf.ne.3) then
              wfrac = max(min(interc*10.0d0/state%crop%ew0/metperiod,1.d0),0.d0)
            else
              wfrac = max(min(eintc/state%crop%ew0,1.0d0),0.0d0)
            endif
          endif
        endif
        ! Remaining amount of interception for swmetdetail = 1
        restint = max(interc - wfrac * metperiod * state%crop%ew0 * 0.1d0, 0.d0)
      endif

      ! === Section 7: Potential soil evaporation & transpiration ===

      ! Potential soil evaporation (peva) [cm/d]
      at_peva = max(0.0d0, (state%crop%es0*dexp(-1.0d0*state%crop%kdir*state%crop%kdif*state%crop%lai)*0.1d0))
      if (state%crop%swcf.ne.3 .or. (config%meteo%swmetdetail.eq.0 .and. state%crop%common%swinter.ne.3)) then
        at_peva = max(0.0d0,(1.0d0-wfrac)*at_peva)
      end if

      ! Alternative for peva (simple model, soil cover fraction specified)
      if (state%crop%common%flCropCalendar .and. .not.flCropHarvest) then   ! [GR-CROP C3]
        if (croptype(state%crop%common%icrop).eq.1 .and. state%crop%common%swgc.eq.2) then
          at_peva = (1.0d0-gc)*state%crop%es0*0.1d0
          if (state%crop%swcf.ne.3 .or. (config%meteo%swmetdetail.eq.0 .and. state%crop%common%swinter.ne.3)) then
            at_peva = (1.0d0-wfrac)*at_peva
          end if
        endif
      endif

      ! Adapt peva in case of ponding — [SS-SWC S-2.12B] state%soilwater%pond
      if (state%soilwater%pond .gt. 1.0d-10) then
        if (config%meteo%swetr.eq.0 .and. state%crop%es0.gt.1.0d-8) then
          at_peva = state%crop%ew0/state%crop%es0 * at_peva
        elseif (state%crop%es0.gt.1.0d-8) then
          if (state%crop%swcfbs .eq. 1 .and. state%crop%cfbs .gt. small) then
            at_peva = cfevappond * at_peva / state%crop%cfbs
          else
            at_peva = cfevappond * at_peva
          endif
        endif
      endif

      ! Potential soil evaporation [cm/d] according to PMdirect
      if (state%cfg%meteo%swdivide .eq. 1) then
        if (state%soilwater%pond .gt. 1.0d-10) then  ! [SS-SWC S-2.12B]
          at_peva = Edirectpond*0.1d0
        else
          at_peva = Edirect*0.1d0
        endif
      endif

      ! Potential transpiration (ptra) [cm/d]
      if (state%crop%swcf .ne. 3) then
        at_ptra = ((1.0d0-wfrac)*state%crop%et0-at_peva*10.0d0)*0.1d0
      else
        at_ptra = (1.0d0-wfrac)*state%crop%et0*0.1d0
      endif
      at_ptra = max(at_ptra,(1.01d0*nihil))

      ! Potential transpiration [cm/d] according to PMdirect
      if (state%cfg%meteo%swdivide .eq. 1) then
        at_ptra = (1.0d0-wfrac) * Tdirect * 0.1d0
        at_ptra = max(at_ptra,(1.01d0*nihil))
      endif

      ! Correction of potential transpiration as a function of atmospheric CO2 concentration
      if (flCO2 .and. state%crop%flCropEmergence) then
        at_ptra = state%crop%wofost%fco2tra * at_ptra
      endif

      ! === Section 8: Results for detailed weather records ===

      ! Result of detailed weather records (cm/d)
      if (config%meteo%swmetdetail.eq.1) then
        state%atmosphere%tpot(irecord) = at_ptra
        state%atmosphere%epot(irecord) = at_peva
        if (state%atmosphere%grai .lt. 1.0d-12) then
          state%atmosphere%grain(irecord) = 0.0d0
          state%atmosphere%nrain(irecord) = 0.0d0
        else
          state%atmosphere%grain(irecord) = arain(irecord) / metperiod
          state%atmosphere%nrain(irecord) = arain(irecord) / metperiod * state%atmosphere%nraida / state%atmosphere%grai
        endif
      endif

1000  continue
    ! === END LOOP over dayparts ===

    ! === Section 9: Actual daily rain/snow fluxes and soil evaporation for Daily Meteo ===

    if (config%meteo%swmetdetail.eq.0) then

      ! Finterception: ratio net / gross rain flux; net rainflux = gross - interception
      if (state%atmosphere%grai.gt.1.d-5) then
        ! finterception is exclusively meant for dividing rain flux into interception part
        ! and net rain part; not for sprinkler irrigation!
        finterception = state%atmosphere%nraida / state%atmosphere%grai
        if (aintc.lt.1.0d-5) finterception = 1.0d0
      else
        finterception = 1.0d0
      endif

      ! In case of daily precipitation sum: set actual gross and net rainflux,
      ! and interception on TIMESTEP basis
      if (config%meteo%swrain.eq.0) then
        rainflux    = state%atmosphere%fprecnosnow * state%atmosphere%grai
        netrainflux = finterception * rainflux
        state%atmosphere%graidt  = rainflux
        state%atmosphere%nraidt  = netrainflux
        if (swuseCN == 1) then
          ! SS-ATM A-2.6: state added — CNmethod needs nraidt/melt from state%atmosphere
          call CNmethod(2, state)
          state%atmosphere%nraidt = state%atmosphere%nraidt - state%atmosphere%Runoff_CN   ! [SS-GR-ATM B24]
        end if
        state%atmosphere%aintcdt = rainflux - state%atmosphere%nraidt  ! aintcdt involves ONLY interception of RAIN
      endif

      ! Soil evaporation rate of today
      if (.not. tc_fletsine) then
        call reduceva (1, state%atmosphere%nraida, state)
      endif

      ! Save daily potential values for use in ETSine
      at_ptraday = at_ptra
      at_pevaday = at_peva

      ! Calculate atmospheric demand [cm]
      at_atmdem = state%crop%et0*0.1d0

    endif

    ! === Section 10: Set daily weather values for Detailed Meteo ===

    if (config%meteo%swmetdetail.eq.1) then
      ! Average temperature of today
      sumtav = 0.d0
      do i = 1, config%meteo%nmetdetail
        sumtav = sumtav + state%atmosphere%atav(i)
      enddo
      state%atmosphere%Tav = sumtav * metperiod
      tav = state%atmosphere%Tav   ! [SS-GR-ATM B24] dual-write — legacy tav consumed by snow.f90/swapoutput.f90

      ! Minimum and maximum temperature of today
      tmx = -50.d0
      tmn = 99.d0
      do i = 1, config%meteo%nmetdetail
        tmx = max(tmx,state%atmosphere%atav(i))
        tmn = min(tmn,state%atmosphere%atav(i))
      enddo

      ! Calculate saturated vapour pressure [kpa]
      svp = 0.3055d0*(dexp(17.27d0*tmn/(tmn+237.3d0)) + &
                      dexp(17.27d0*tmx/(tmx+237.3d0)))
      ! Calculate relative humidity [fraction]
      state%atmosphere%rh = min(hum/svp,1.0d0)   ! [SS-GR-ATM B24] direct state write
      ! [SS-GR-ATM B.5] rh dual-write to legacy global RETIRED: no crop consumers remain

      ! Average temperature between 6 and 18 hour
      sumtav = 0.d0
      count = 0
      first = int(0.25/metperiod) + 1
      last = int(0.75/metperiod)
      do i = first, last
        sumtav = sumtav + state%atmosphere%atav(i)
        count = count + 1
      enddo
      state%atmosphere%tavd = sumtav / count   ! [SS-GR-ATM B24] direct state write
      ! [SS-GR-ATM B.5] tavd dual-write to legacy global RETIRED: no crop consumers remain

      ! Daily radiation (J/m2/d) and atmospheric demand (cm/d)
      rad = 0.d0
      at_atmdem = 0.d0
      do i = 1,config%meteo%nmetdetail
        rad = rad + state%atmosphere%arad(i)
        at_atmdem = at_atmdem + state%atmosphere%tpot(i)
      enddo

      ! Fluxes of current time step (start of the day)
      at_ptra = state%atmosphere%tpot(1)
      at_peva = state%atmosphere%epot(1)
      state%atmosphere%graidt  = state%atmosphere%grain(1)
      state%atmosphere%nraidt  = state%atmosphere%nrain(1)
      state%atmosphere%aintcdt = state%atmosphere%graidt - state%atmosphere%nraidt  ! aintcdt involves ONLY interception of RAIN

    endif

    end associate  ! at_peva/at_ptra/at_atmdem/at_pevaday/at_ptraday + tc_daynr/tc_t/tc_dt/tc_flmetdetail/tc_fletsine [TC-9]

  end subroutine ProcessMeteoDay

end module meteo_mod










