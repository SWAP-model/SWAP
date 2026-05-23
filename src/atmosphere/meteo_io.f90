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
      use precipitation_mod, only: PartitionPrecipitation
      implicit none

      type(swap_state_t),   intent(inout) :: state   !! [SS-ATM] threaded for atmosphere dual-writes
      type(swap_config_t),  intent(in)    :: config  !! [SS-GR-ATM B23] for meteo config switches

    ! --- local
    ! [SS-ATM A-2.6] grai/gsnow/snrai/ssnow/fprecnosnow are transitional locals fed to
    ! PartitionPrecipitation; state%atmosphere%X is the canonical write target.
    real(8) :: grai, gsnow, snrai, ssnow, fprecnosnow
    ! [GR-ATM-CLEAN Phase D.2] formerly module MeteoVars members
    integer :: i               ! sub-daily record loop counter
    real(8) :: hum, win, etr   ! within-call scratch (daily branch)
    real(8) :: svp             ! within-call saturated vapor pressure
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
        svp = 0.3055d0*(exp(17.27d0*tmn/(tmn+237.3d0)) + &
                        exp(17.27d0*tmx/(tmx+237.3d0)))
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
        state%atmosphere%awind_subdaily(i) = detwind(irectotal)
        state%atmosphere%arain_subdaily(i) = detrain(irectotal) * 0.1d0 ! convert from mm to cm
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
                                ssnow, config%meteo%nmetdetail, state%atmosphere%arain_subdaily, &
                                grai, gsnow, snrai, &
                                fprecnosnow, state%atmosphere%restint, state)

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

    ! Reset intermediate fluxes + snapshot intermediate-period baseline.
    ! Mirrors the soilwater coordinator pattern (soilhydraulics.f90:1190-1196).
    if (state%timecontrol%flZeroIntr) then
      call state%atmosphere%intr%reset()
      state%atmosphere%ISsnowBeg = state%atmosphere%ssnow
    endif

    ! Reset cumulative fluxes + snapshot cumulative-period baseline.
    if (state%timecontrol%flZeroCumu) then
      call state%atmosphere%cumu%reset()
      state%atmosphere%snowinco = state%atmosphere%ssnow
    endif

    return
  end subroutine ResetMetFlx

end module meteo_process_mod
