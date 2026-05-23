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
   use error_mod,       only: fatalerr_collected
   use swap_state_mod,  only: swap_state_t
   use swap_config_mod, only: swap_config_t
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
      ! DEFERRED — rad/tmn/tmx/tav (daily meteo scalars) dual-written to bare
      ! globals because downstream consumers (meteo_orchestrator, cropgrowth)
      ! still read them. Migration to state%atmosphere is a coordinated
      ! multi-file commit; pathatm/metfil/det* arrays + irectotal stay too.
      use variables, only: rad, tmn, tmx, tav,                  &
                           pathatm, metfil,                      &
                           detrecord, dettime, detrad, dethum,   &
                           dettav, detrain, detwind, irectotal
      use precipitation_mod, only: PartitionPrecipitation
      implicit none

      type(swap_state_t),  intent(inout) :: state
      type(swap_config_t), intent(in)    :: config

      integer :: i, today_idx
      real(8) :: hum, win, etr, svp
      character(len=11)  :: detdate
      character(len=3)   :: ext
      character(len=200) :: filnam
      character(len=300) :: messag

      associate (atmo => state%atmosphere, time => state%timecontrol)

      call ResetMetFlx(state)

      ! ===== Daily meteo =====
      if (config%meteo%swmetdetail .eq. 0) then

         ! Check availability of meteo data of today
         if (time%daymeteo .lt. atmo%daynrfirst .or. time%daymeteo .gt. atmo%daynrlast) then
            messag = 'In meteo file no meteo data are available for ' &
                     // time%date // '. First adapt meteo file!'
            call fatalerr_collected('meteo', messag)
         end if

         ! Pass on weather values of today
         today_idx = time%daymeteo + 1 - atmo%daynrfirst
         rad  = atmo%arad(today_idx)
         tmn  = atmo%atmn(today_idx)
         tmx  = atmo%atmx(today_idx)
         hum  = atmo%ahum(today_idx)
         win  = atmo%awin(today_idx)
         etr  = atmo%aetr(today_idx)

         ! If hum is missing or tav cannot be calculated: set rh at -99.0
         atmo%rh = 1.0d0
         if (hum .lt. -98.0d0 .or. tmn .lt. -98.0d0 .or. tmx .lt. -98.0d0) atmo%rh = -99.0d0

         ! 24h average + day temperature
         atmo%Tav  = (tmx + tmn) * 0.5d0
         tav       = atmo%Tav   ! dual-write — bare 'tav' still consumed in meteo_orchestrator subdaily
         atmo%tavd = (tmx + atmo%Tav) * 0.5d0

         if (atmo%rh .ge. -98.0d0) then
            ! Saturated vapour pressure [kPa]
            svp = 0.3055d0 * (exp(17.27d0*tmn/(tmn+237.3d0)) + exp(17.27d0*tmx/(tmx+237.3d0)))
            atmo%rh = min(hum/svp, 1.0d0)
         endif

         ! CFO output snapshot (PEARL coupling) — sole writers to atmo%out_*
         atmo%out_rad = real(rad, kind=8)
         atmo%out_tmn = real(tmn, kind=8)
         atmo%out_tmx = real(tmx, kind=8)
         atmo%out_hum = real(hum, kind=8)
         atmo%out_win = real(win, kind=8)
         atmo%out_etr = real(etr, kind=8) * 0.001d0
         if (config%meteo%swrain .eq. 2) then
            atmo%out_wet = real(atmo%wet(today_idx), kind=8)
         else
            atmo%out_wet = -1.0d0
         endif

      ! ===== Detailed meteo =====
      elseif (config%meteo%swmetdetail .eq. 1) then

         ! Compose filename meteorological file for use in warnings
         write (ext, '(i3.3)') mod(time%yearmeteo, 1000)
         filnam = trim(pathatm) // trim(metfil) // '.' // trim(ext)

         do i = 1, config%meteo%nmetdetail
            irectotal = irectotal + 1
            if (i .ne. detrecord(irectotal)) then
               messag = 'In meteo file '// trim(filnam) // ' record number(s)' &
                        // ' are not correct at ' // time%date // '. First adapt meteo file!'
               call fatalerr_collected('meteo', messag)
            end if
            call dtdpst('year-month-day', dettime(irectotal) + 0.1d0, detdate)
            call dtdpst('year-month-day', time%t1900       + 0.1d0, time%date)
            if (detdate .ne. time%date) then
               messag = 'In meteo file ' // trim(filnam) // ' the amount of ' &
                        // 'records deviate near ' // time%date // '. First adapt meteo file!'
               call fatalerr_collected('meteo', messag)
            end if

            ! Pass on weather records of today
            atmo%arad(i)           = detrad(irectotal)
            atmo%ahum(i)           = dethum(irectotal)
            atmo%atav(i)           = dettav(irectotal)
            atmo%awind_subdaily(i) = detwind(irectotal)
            atmo%arain_subdaily(i) = detrain(irectotal) * 0.1d0   ! mm → cm
         enddo
      endif

      ! Partition today's precipitation (reads + writes via state%atmosphere).
      call PartitionPrecipitation(state, config)

      end associate
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
  subroutine ResetMetFlx(state)
      implicit none

      type(swap_state_t), intent(inout) :: state

      associate (atmo => state%atmosphere, time => state%timecontrol)

      ! Reset intermediate fluxes + snapshot intermediate-period baseline.
      ! Mirrors the soilwater coordinator pattern (soilhydraulics.f90:1190-1196).
      if (time%flZeroIntr) then
         call atmo%intr%reset()
         atmo%ISsnowBeg = atmo%ssnow
      endif

      ! Reset cumulative fluxes + snapshot cumulative-period baseline.
      if (time%flZeroCumu) then
         call atmo%cumu%reset()
         atmo%snowinco = atmo%ssnow
      endif

      end associate
   end subroutine ResetMetFlx

end module meteo_process_mod
