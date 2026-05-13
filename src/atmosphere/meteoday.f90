!> Module containing shared meteorological variables
!!
!! This module provides shared variables for communication between meteorological
!! processing routines (ReadMeteoDay and ProcessMeteoDay). Variables are used to
!! store temporary meteorological data, intermediate calculations, and results
!! during the processing of daily or sub-daily meteorological inputs.
!!
!! ## Variable Categories
!! - Loop control: count, first, i, irecord, last, ndayparts
!! - Meteorological arrays: arain(96), awind(96)
!! - Interception: restint, interc, aintc, eintc
!! - Flux components: netrainflux, rainflux, Edirect, Tdirect, Tdirectwet, Edirectpond
!! - Meteorological scalars: etr, gctp, hum, svp, wfrac, win, dttp, sumtav
!!
!! @author Original SWAP development team
!! @date Refactored February 2026
module MeteoVars
   implicit none
   
   integer   count           !! Loop counter for meteorological records
   integer   first           !! First record index in processing range
   integer   i               !! General loop counter
   integer   irecord         !! Current record index
   integer   last            !! Last record index in processing range
   integer   ndayparts       !! Number of timesteps per day (1 for daily, >1 for sub-daily)
   
   real(8)   arain(96)       !! Array of precipitation values (cm)
   real(8)   awind(96)       !! Array of wind speed values (m/s)
   real(8)   restint         !! Remaining interception from previous timestep (cm)
   real(8)   interc          !! Current interception amount (cm)
   real(8)   Edirectpond     !! Direct evaporation from ponded water (mm/d)
   real(8)   aintc           !! Daily interception (cm)
   real(8)   dttp            !! Timestep duration for interception calculations (d)
   real(8)   eintc           !! Evaporated interception (cm)
   real(8)   etr             !! Reference evapotranspiration (mm/d)
   real(8)   gctp            !! Ground cover at current timestep (fraction)
   real(8)   hum             !! Humidity or vapor pressure (kPa)
   real(8)   netrainflux     !! Net rainfall flux after interception (cm)
   real(8)   rainflux        !! Gross rainfall flux (cm)
   real(8)   sumtav          !! Sum of temperatures for averaging
   real(8)   svp             !! Saturated vapor pressure (kPa)
   real(8)   wfrac           !! Wet fraction of crop canopy (fraction)
   real(8)   win             !! Wind speed (m/s)
   real(8)   Edirect         !! Direct soil evaporation (mm/d)
   real(8)   Tdirect         !! Direct transpiration (mm/d)
   real(8)   Tdirectwet      !! Direct transpiration from wet canopy (mm/d)
end module MeteoVars

!> Module for surface runoff calculation using the SCS Curve Number method
!!
!! This module implements the USDA Soil Conservation Service (SCS) Curve Number (CN) method
!! for estimating direct surface runoff from rainfall events. The method accounts for:
!! - Time-varying CN values through lookup tables
!! - Soil moisture corrections (dry, normal, and wet conditions)
!! - Snowmelt contribution to runoff
!!
!! The CN method relates runoff to rainfall through the empirical equation:
!! \[ Q = \frac{(P - I_a)^2}{P - I_a + S} \]
!! where \(Q\) is runoff depth, \(P\) is precipitation depth, \(I_a\) is initial abstraction,
!! and \(S\) is maximum potential retention.
!!
!! @author Original SWAP development team
!! @date Refactored February 2026
module runoff_mod
   use error_mod, only: fatalerr_collected
   use soilhydraulics_utils, only: watcon
   use swap_state_mod, only: swap_state_t  ! [SS-ATM A-2.6] nraidt/melt retired to state%atmosphere
   implicit none
   private
   
   public :: CNmethod

   ! Physical constants
   real(8), parameter :: DEPTH_10CM = 10.0d0           !! Reference depth for moisture correction (cm)
   real(8), parameter :: H_FIELD_CAPACITY = -100.0d0   !! Pressure head at field capacity (cm)
   real(8), parameter :: H_WILTING_POINT = -16000.0d0  !! Pressure head at wilting point (cm)
   real(8), parameter :: IA_RATIO = 0.2d0              !! Initial abstraction ratio (dimensionless)

contains

  subroutine CNmethod(Itask, state)
    !> Calculate surface runoff using the SCS Curve Number method
    !!
    !! This subroutine operates in two modes controlled by the Itask parameter:
    !! - Itask=1: Initialization - validates time series, identifies top soil layer (0-10 cm),
    !!            and calculates reference water content for moisture corrections
    !! - Itask=2: Dynamic calculation - computes runoff for the current timestep using
    !!            moisture-adjusted CN values
    !!
    !! ## Moisture Correction Options (wc_cor)
    !! - 0: No moisture correction (CN = CNref)
    !! - 1: Field capacity based (\(h = -100\) cm and \(h = -16000\) cm)
    !! - 2: Saturation based (\(h = 0\) cm and \(h = -16000\) cm)
    !!
    !! ## References
    !! USDA-NRCS National Engineering Handbook, Part 630 Hydrology
    !!
    !! @warning CNref values > 170 will cause numerical issues in CNdry calculation
    ! [SS-ATM A-2.6] nraidt/melt retired from variables; state added to read from state%atmosphere
    ! [SS-SWC S-2.12B] theta retired — read via state%soilwater%theta
    ! SS-TC TC-9: t1900 removed from only-list; read via state%timecontrol.
    ! [SS-GR-ATM B22] CN symbols → state%atmosphere%X (nod10_cn/icn_atm/z10_cn added to atmosphere_state)
    ! GR-BH: numnod, zbotcp, dz migrated to state%mesh%X
    use soilhydraulics_utils, only: watcon
    implicit none
    ! global
    integer, intent(in)  :: Itask
    type(swap_state_t), intent(inout) :: state  ! [SS-ATM A-2.6] for retired nraidt/melt
    ! local
    integer              :: i
    real(8)              :: wc1, wc2, CN, S, Ia
    ! Note: nod10_cn/icn_atm/z10_cn migrated to state%atmosphere (B22)

    ! SS-TC TC-9: t1900 read via state%timecontrol (tc_* alias).
    associate( tc_t1900 => state%timecontrol%t1900 )  ! TC-9

    select case (Itask)
    ! initialization; calculate and store some constants
    case (1)

      state%atmosphere%icn_atm = 0
      ! check if times in CNtimeTAB are in ascending order
      ! set initial position in CNtimTAB
      do i = 2, state%atmosphere%iCNtab
          if (state%atmosphere%CNtimTAB(i) < state%atmosphere%CNtimTAB(i-1)) call fatalerr_collected ('CNmethod', 'CNtimTAB not in ascending order')
          if (tc_t1900 >= state%atmosphere%CNtimTAB(i-1) .and. tc_t1900 < state%atmosphere%CNtimTAB(i)) state%atmosphere%icn_atm = i-1
      end do
      ! error if start time t1900 not in CNtimTAB
      if (state%atmosphere%icn_atm == 0) call fatalerr_collected ('CNmethod', 'Start time of simulation not present in CNtimTAB')

    !  to be replaced by average for layer 0-10 cm
      do i = 1, state%mesh%numnod
          if (state%mesh%zbotcp(i) < -DEPTH_10CM) then
            state%atmosphere%nod10_cn = i-1
            state%atmosphere%z10_cn = -state%mesh%zbotcp(state%atmosphere%nod10_cn)
            exit
          end if
      end do
      state%atmosphere%ThetaRef = 0.0d0
      do i = 1, state%atmosphere%nod10_cn
          if (state%atmosphere%wc_cor == 1) then
            wc1 = watcon(H_FIELD_CAPACITY, &
                          state%soilwater%vg_params(i), &
                          state%soilwater%iHWCKmodel(state%soilwater%layer(i)), &
                          i, state%soilwater)                      ! [SS-GR-UTILS Task 5]
            wc2 = watcon(H_WILTING_POINT, &
                          state%soilwater%vg_params(i), &
                          state%soilwater%iHWCKmodel(state%soilwater%layer(i)), &
                          i, state%soilwater)                      ! [SS-GR-UTILS Task 5]
            state%atmosphere%ThetaRef = state%atmosphere%ThetaRef + (wc1+wc2)*0.5d0*state%mesh%dz(i)
          else if (state%atmosphere%wc_cor == 2) then
            wc1 = watcon(0.0d0, &
                          state%soilwater%vg_params(i), &
                          state%soilwater%iHWCKmodel(state%soilwater%layer(i)), &
                          i, state%soilwater)                      ! [SS-GR-UTILS Task 5]
            wc2 = watcon(H_WILTING_POINT, &
                          state%soilwater%vg_params(i), &
                          state%soilwater%iHWCKmodel(state%soilwater%layer(i)), &
                          i, state%soilwater)                      ! [SS-GR-UTILS Task 5]
            state%atmosphere%ThetaRef = state%atmosphere%ThetaRef + (wc1+wc2)*0.5d0*state%mesh%dz(i)
          end if
      end do
      state%atmosphere%ThetaRef = state%atmosphere%ThetaRef/state%atmosphere%z10_cn
      
      !!!t1900_old = int(t1900) - 1 ! for testing intermediate output
      
    ! dynamic part: calculate runoff   
    case (2)
      
      ! see if t1900 has moved ahead in CNtimTAB; icn_atm can never exceed last entry
      !  if (icn_atm < iCNtab .and. t1900 >= CNtimTAB(icn_atm+1)) icn_atm = icn_atm + 1
      ! Update position in CN time series if time has advanced (do while is more efficient if time steps are large and CN time series is long)
      do while (state%atmosphere%icn_atm < state%atmosphere%iCNtab .and. tc_t1900 >= state%atmosphere%CNtimTAB(state%atmosphere%icn_atm + 1))
        state%atmosphere%icn_atm = state%atmosphere%icn_atm + 1
      end do
      state%atmosphere%CNref = state%atmosphere%CNrefTAB(state%atmosphere%icn_atm)
      CN    = state%atmosphere%CNref
      state%atmosphere%CNdry =  4.2d0*state%atmosphere%CNref/(10.0d0-0.058d0*state%atmosphere%CNref)
      state%atmosphere%CNwet = 23.0d0*state%atmosphere%CNref/(10.0d0+0.13d0*state%atmosphere%CNref)

      if (state%atmosphere%wc_cor > 0) then
          state%atmosphere%wc10 = 0.0d0
          do i = 1, state%atmosphere%nod10_cn
            state%atmosphere%wc10 = state%atmosphere%wc10 + state%soilwater%theta(i)*state%mesh%dz(i)  ! [SS-SWC S-2.12B]
          end do
          state%atmosphere%wc10 = state%atmosphere%wc10/state%atmosphere%z10_cn
          if (state%atmosphere%wc10 < state%atmosphere%ThetaRef) then
            CN = state%atmosphere%CNdry + state%atmosphere%wc10/state%atmosphere%ThetaRef*(state%atmosphere%CNref-state%atmosphere%CNdry)
          else
            CN = state%atmosphere%CNref + (state%atmosphere%wc10-state%atmosphere%ThetaRef)/state%atmosphere%ThetaRef*(state%atmosphere%CNwet-state%atmosphere%CNref)
          end if
      end if
      S = 2540d0/CN-25.4d0    ! in cm
      Ia = 0.2d0*S
    !   Ia = 0.3d0*S
      ! SS-ATM A-2.6: nraidt/melt retired — read from state%atmosphere
      if (state%atmosphere%nraidt+state%atmosphere%melt > Ia) then
          state%atmosphere%Runoff_CN = (state%atmosphere%nraidt+state%atmosphere%melt-Ia)**2/ &
                      (state%atmosphere%nraidt+state%atmosphere%melt-Ia+S)
      else
          state%atmosphere%Runoff_CN = 0.0d0
      end if
      
      ! for testing intermediate output
      !!!if (int(t1900) > t1900_old) then
      !!!   write (123, '(7F20.6)') t1900, nraidt, runoff_cn, cn, wc10, thetaref, melt
      !!!   t1900_old = int(t1900)
      !!!end if
      
    case default
      call fatalerr_collected ('CNmethod', 'Illegal Itask option')
    end select

    end associate  ! tc_t1900 => state%timecontrol [TC-9]

  end subroutine CNmethod

end module runoff_mod

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
  subroutine ReadMeteoDay(state)
      ! use variables
      ! SS-TC TC-9: date,t1900 removed from only-list; reads/writes via state%timecontrol.
      ! [SS-TC TC-14] yearmeteo, daymeteo retired — read via state%timecontrol
      use variables, only: out_tmn, out_tmx, out_hum, out_win, out_etr, out_wet, swrain, wet, rh, tav, tavd, out_rad, arai, atmx, ahum, aetr, arad, teprrain, teprsnow, &
                        detrecord, nmetdetail, dettime, detrad, dethum, dettav, atav, swmetdetail, daynrfirst, daynrlast, rad, tmn, tmx, pathatm, awin, atmn, metfil, detrain, swsnow, irectotal, detwind
      use MeteoVars
      use precipitation_mod, only: PartitionPrecipitation
      implicit none

      type(swap_state_t), intent(inout) :: state  !! [SS-ATM] threaded for atmosphere dual-writes

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
    
    if (swmetdetail.eq.0) then

      ! Check availability of meteo data of today
      if (daymeteo.lt.daynrfirst .or. daymeteo.gt.daynrlast) then
        messag ='In meteo file no meteo data are'// &
                ' available for '//tc_date//'. First adapt meteo file!'
        call fatalerr_collected ('meteo',messag)
      end if

      ! Pass on weather values of today
      rad  = arad(daymeteo+1-daynrfirst)
      tmn  = atmn(daymeteo+1-daynrfirst)
      tmx  = atmx(daymeteo+1-daynrfirst)
      hum  = ahum(daymeteo+1-daynrfirst)
      win  = awin(daymeteo+1-daynrfirst)
      grai = arai(daymeteo+1-daynrfirst)  ! [SS-ATM A-2.6] local; PartitionPrecipitation converts mm->cm and writes state
      etr  = aetr(daymeteo+1-daynrfirst)

      ! If hum is missing or tav cannot be calculated: set rh at -99.0
      rh = 1.0d0
      if (hum.lt.-98.0d0 .or. tmn.lt.-98.0d0 .or. tmx.lt.-98.0d0) &
        rh=-99.0d0

      ! Calculate 24h average temperature
      tav = (tmx+tmn)*0.5d0
      ! Calculate average day temperature
      tavd = (tmx+tav)*0.5d0
      state%atmosphere%tavd = tavd   ! [SS-GR-ATM A5.5] runtime dual-write

      if (rh.ge.-98.0d0) then
        ! Calculate saturated vapour pressure [kpa]
        svp = 0.3055d0*(dexp(17.27d0*tmn/(tmn+237.3d0)) + &
                        dexp(17.27d0*tmx/(tmx+237.3d0)))
        ! Calculate relative humidity [fraction]
        rh = min(hum/svp,1.0d0)
      endif
      state%atmosphere%rh = rh   ! [SS-GR-ATM A5.5] runtime dual-write

      ! CFO file for PEARL: save meteo variables of today for output
      out_rad = real(rad)
      state%atmosphere%out_rad = real(out_rad, kind=8)   ! [SS-GR-ATM A5.5] runtime dual-write (real4→real64)
      out_tmn = real(tmn)
      state%atmosphere%out_tmn = real(out_tmn, kind=8)   ! [SS-GR-ATM A5.5] runtime dual-write (real4→real64)
      out_tmx = real(tmx)
      state%atmosphere%out_tmx = real(out_tmx, kind=8)   ! [SS-GR-ATM A5.5] runtime dual-write (real4→real64)
      out_hum = real(hum)
      state%atmosphere%out_hum = real(out_hum, kind=8)   ! [SS-GR-ATM A5.5] runtime dual-write (real4→real64)
      out_win = real(win)
      state%atmosphere%out_win = real(out_win, kind=8)   ! [SS-GR-ATM A5.5] runtime dual-write (real4→real64)
      out_etr = real(etr)*0.001
      state%atmosphere%out_etr = real(out_etr, kind=8)   ! [SS-GR-ATM A5.5] runtime dual-write (real4→real64)
      if (swrain.eq.2) then
        out_wet = real(wet(daymeteo+1-daynrfirst))
      else
        out_wet = -1.0
      endif
      state%atmosphere%out_wet = real(out_wet, kind=8)   ! [SS-GR-ATM A5.5] runtime dual-write (real4→real64)
    
    ! end 1 Daily Meteo 00000000000000000000000000000000000000000000000000000 Daily Meteo
    
    ! 1.1 Detailed Meteo 1111111111111111111111111111111111111111111111111111 Detailed Meteo
    
    elseif (swmetdetail.eq.1) then

      ! Check availability of meteo data of today
      ! Compose filename meteorological file for use in warnings
      write (ext,'(i3.3)') mod(yearmeteo,1000)
      filnam = trim(pathatm)//trim(metfil)//'.'//trim(ext)

      do i = 1, nmetdetail
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
        arad(i)  = detrad(irectotal)
        ahum(i)  = dethum(irectotal)
        atav(i)  = dettav(irectotal)
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
    call PartitionPrecipitation(swmetdetail, swsnow, tav, TePrRain, TePrSnow, &
                                ssnow, nmetdetail, arain, grai, gsnow, snrai, &
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
  use et_mod, only: PenMon, reduceva
  use runoff_mod, only: CNmethod
  use swap_state_mod, only: swap_state_t

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
  subroutine ProcessMeteoDay(state)
    ! use Variables
    ! SS-TC TC-9: daynr,t,dt,flmetdetail,fletsine removed from only-list; reads via state%timecontrol.
    ! [SS-TC TC-14] metperiod retired — read via state%timecontrol
    use variables, only: lai, gird, swinter, swmetdetail, nmetdetail, swetr, flCropEmergence, et0, ew0, es0, swcf, swcfbs, cfbs, &
    cf, cfeic, rad, arad, tav, atav, ahum, logf, lat, alt, altw, angstroma, angstromb, rsc, ch, daylp, albedo, tmn, tmx, rsw, difpp, &
    dsinbe, atmtr, rsoil, swdivide, kdif, kdir, croptype, swgc, gc, siccapact, siccaptb, icrop, flcropcalendar, &
     flCropHarvest, cfevappond, flco2, fco2tra, tpot, epot, grain, nrain, finterception, swrain, &
     swusecn, runoff_cn, rh, tavd
     ! [SS-SWC S-2.12B] pond retired — read via state%soilwater%pond
    use swap_array_dimensions, only: magrs
    use MeteoVars
    use array_utils, only: afgen
    use runoff_mod, only: CNmethod
    use interception_mod, only: VonHHBraden, Gash, ruttervw, msw1eic, DivIntercep
    use swap_constants, only: nihil, small
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none

    type(swap_state_t), intent(inout) :: state
      !! Simulation state (passed through to reduceva for atmosphere dual-writes)

    real(8)  rcs
    data     rcs/0.15d0/

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
    if ((lai .lt. 1.d-3) .or. (state%atmosphere%grai+gird .lt. 1.d-5) .or. &
        (swinter.eq.0) .or. (state%atmosphere%gsnow.gt.0.0d0) .or.(state%atmosphere%ssnow.gt.0.0d0)) then

      ! No vegetation, rainfall/irrigation or interception calculation
      aintc = 0.d0

    else if (swinter .eq. 1) then
      ! Calculate interception, method Von Hoyningen-Hune and Braden
      ! SS-ATM A-2.6: grai retired — pass state%atmosphere%grai explicitly
      ! SS-GR-ATM B8: state added for crop/atmosphere fields
      call VonHHBraden (aintc, state%atmosphere%grai, state)
    else if (swinter .eq. 2) then
      ! Calculate interception, method Gash (1995)
      ! SS-ATM A-2.6: grai retired — pass state%atmosphere%grai explicitly
      ! SS-TC TC-11: state added for t via state%timecontrol%t
      call Gash (aintc, state%atmosphere%grai, state)
    end if

    ! Divide interception into rain part and irrigation part and
    ! calculate net rain (nraida) and net sprinkling irrigation (nird)
    if (swinter.ne.3) &
      call DivIntercep (aintc, state)

    ! === LOOP over dayparts ===

    ! Set number of repetitions within one day
    if (swmetdetail.eq.0) then
      ! In case of daily meteo input: only one time per day
      ndayparts = 1
    elseif (swmetdetail.eq.1) then
      ndayparts = nmetdetail
    endif

    do 1000 irecord = 1, ndayparts

      ! === Section 4: Calculate evapotranspiration (et0, ew0, es0) ===
      ! Reference evapotranspiration has been specified
      if (swmetdetail.eq.0 .and. swetr.eq.1) then
        if (.not. flCropEmergence) then
          ! no crop
          et0 = 0.0d0
          ew0 = 0.0d0
          es0 = etr
          if (swcfbs.eq.1) es0 = cfbs*etr
        else
          ! crop is present
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

      ! Reference evapotranspiration must be calculated
      elseif (swmetdetail.eq.1 .or. swetr.eq.0) then

        if (swmetdetail.eq.1) then
          ! Define weather variables of current record
          rad = arad(irecord) / metperiod     ! from j/m2/period to j/m2/d
          tav = atav(irecord)
          hum = ahum(irecord)
          win = awind(irecord)
        endif

        ! Calculate evapotranspiration using Penman-Monteith: et0, ew0, es0 (mm/d)
        ! in case of daily meteo (swmetdetail = 0) irecord is always 1
        call PenMon (logf,swscre,tc_daynr,lat,alt,Altw,angstroma, &
                     angstromb,rcs,rad,tav,hum,win,rsc,es0,et0,ew0,swcf,ch, &
                     flCropEmergence,daylp,tc_flmetdetail,irecord, &
                     nmetdetail,albedo,tmn,tmx,rsw,difpp,dsinbe,atmtr, &
                     Edirect,Tdirect,Tdirectwet,rsoil,swdivide,kdif,kdir, &
                     lai,Edirectpond)

        if (.not. flCropEmergence) then
          ! no crop
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
          ! crop is present
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

      endif
      ! [SS-GR-ATM A5.5] mirror final et0/ew0/es0 values to state%crop
      state%crop%et0 = et0
      state%crop%ew0 = ew0
      state%crop%es0 = es0

      ! === Section 5: Interception option NHI (adapted Rutter model) ===

      if (swinter .eq. 3) then
        ! Set parameter values
        if (swmetdetail.eq.0) then  ! if swmetdetail = 0, siccapact is set in cropgrowth module
          if (croptype(icrop).eq.1 .and. swgc.eq.2) then
            gctp  = gc
          else
            gctp  = 1.0d0 - dexp(-1.0d0*kdir*kdif*lai)
            if (gctp .lt. 1.0d-5) then
              siccapact=0.
              state%atmosphere%siccapact = siccapact   ! [SS-GR-ATM A5.5] dual-write
            endif
          endif
          dttp = 1.0d0  ! value of 1 d required for the daily meteo option
        elseif (swmetdetail.eq.1) then
          siccapact = afgen(siccaptb,(2*magrs),tc_t)
          state%atmosphere%siccapact = siccapact   ! [SS-GR-ATM A5.5] dual-write
          if (croptype(icrop).eq.1 .and. swgc.eq.2) then
            gctp  = gc
          elseif (croptype(icrop).eq.1 .and. swgc.eq.1) then
            gctp  = 1.0d0 - dexp(-1.0d0*kdir*kdif*lai)
          endif
          if (gctp .lt. 1.0d-5) then
            siccapact=0.
            state%atmosphere%siccapact = siccapact   ! [SS-GR-ATM A5.5] dual-write
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
      if (swmetdetail.eq.0) then
        ! Fraction of the day the crop is wet
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
      ! Fraction of the period the crop is wet
      elseif (swmetdetail.eq.1) then
        if (state%atmosphere%grai .lt. 1.0d-12) then
          interc = 0.0d0
          wfrac  = 0.0d0
        else
          interc = restint + aintc * arain(irecord) / state%atmosphere%grai
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
        ! Remaining amount of interception for swmetdetail = 1
        restint = max(interc - wfrac * metperiod * ew0 * 0.1d0, 0.d0)
      endif

      ! === Section 7: Potential soil evaporation & transpiration ===

      ! Potential soil evaporation (peva) [cm/d]
      at_peva = max(0.0d0, (es0*dexp(-1.0d0*kdir*kdif*lai)*0.1d0))
      if (swcf.ne.3 .or. (swmetdetail.eq.0 .and. swinter.ne.3)) then
        at_peva = max(0.0d0,(1.0d0-wfrac)*at_peva)
      end if

      ! Alternative for peva (simple model, soil cover fraction specified)
      if (flCropCalendar .and. .not.flCropHarvest) then
        if (croptype(icrop).eq.1 .and. swgc.eq.2) then
          at_peva = (1.0d0-gc)*es0*0.1d0
          if (swcf.ne.3 .or. (swmetdetail.eq.0 .and. swinter.ne.3)) then
            at_peva = (1.0d0-wfrac)*at_peva
          end if
        endif
      endif

      ! Adapt peva in case of ponding — [SS-SWC S-2.12B] state%soilwater%pond
      if (state%soilwater%pond .gt. 1.0d-10) then
        if (SwETr.eq.0 .and. es0.gt.1.0d-8) then
          at_peva = ew0/es0 * at_peva
        elseif (es0.gt.1.0d-8) then
          if (swcfbs .eq. 1 .and. cfbs .gt. small) then
            at_peva = cfevappond * at_peva / cfbs
          else
            at_peva = cfevappond * at_peva
          endif
        endif
      endif

      ! Potential soil evaporation [cm/d] according to PMdirect
      if (swdivide .eq. 1) then
        if (state%soilwater%pond .gt. 1.0d-10) then  ! [SS-SWC S-2.12B]
          at_peva = Edirectpond*0.1d0
        else
          at_peva = Edirect*0.1d0
        endif
      endif

      ! Potential transpiration (ptra) [cm/d]
      if (swcf .ne. 3) then
        at_ptra = ((1.0d0-wfrac)*et0-at_peva*10.0d0)*0.1d0
      else
        at_ptra = (1.0d0-wfrac)*et0*0.1d0
      endif
      at_ptra = max(at_ptra,(1.01d0*nihil))

      ! Potential transpiration [cm/d] according to PMdirect
      if (swdivide .eq. 1) then
        at_ptra = (1.0d0-wfrac) * Tdirect * 0.1d0
        at_ptra = max(at_ptra,(1.01d0*nihil))
      endif

      ! Correction of potential transpiration as a function of atmospheric CO2 concentration
      if (flCO2 .and. flCropEmergence) then
        at_ptra = fco2tra * at_ptra
      endif

      ! === Section 8: Results for detailed weather records ===

      ! Result of detailed weather records (cm/d)
      if (swmetdetail.eq.1) then
        tpot(irecord) = at_ptra
        epot(irecord) = at_peva
        if (state%atmosphere%grai .lt. 1.0d-12) then
          grain(irecord) = 0.0d0
          nrain(irecord) = 0.0d0
        else
          grain(irecord) = arain(irecord) / metperiod
          nrain(irecord) = arain(irecord) / metperiod * state%atmosphere%nraida / state%atmosphere%grai
        endif
      endif

1000  continue
    ! === END LOOP over dayparts ===

    ! === Section 9: Actual daily rain/snow fluxes and soil evaporation for Daily Meteo ===

    if (swmetdetail.eq.0) then

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
      if (swrain.eq.0) then
        rainflux    = state%atmosphere%fprecnosnow * state%atmosphere%grai
        netrainflux = finterception * rainflux
        state%atmosphere%graidt  = rainflux
        state%atmosphere%nraidt  = netrainflux
        if (swuseCN == 1) then
          ! SS-ATM A-2.6: state added — CNmethod needs nraidt/melt from state%atmosphere
          call CNmethod(2, state)
          state%atmosphere%nraidt = state%atmosphere%nraidt - Runoff_CN
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
      at_atmdem = et0*0.1d0

    endif

    ! === Section 10: Set daily weather values for Detailed Meteo ===

    if (swmetdetail.eq.1) then
      ! Average temperature of today
      sumtav = 0.d0
      do i = 1, nmetdetail
        sumtav = sumtav + atav(i)
      enddo
      tav = sumtav * metperiod

      ! Minimum and maximum temperature of today
      tmx = -50.d0
      tmn = 99.d0
      do i = 1, nmetdetail
        tmx = max(tmx,atav(i))
        tmn = min(tmn,atav(i))
      enddo

      ! Calculate saturated vapour pressure [kpa]
      svp = 0.3055d0*(dexp(17.27d0*tmn/(tmn+237.3d0)) + &
                      dexp(17.27d0*tmx/(tmx+237.3d0)))
      ! Calculate relative humidity [fraction]
      rh = min(hum/svp,1.0d0)
      state%atmosphere%rh = rh   ! [SS-GR-ATM A5.5] runtime dual-write

      ! Average temperature between 6 and 18 hour
      sumtav = 0.d0
      count = 0
      first = int(0.25/metperiod) + 1
      last = int(0.75/metperiod)
      do i = first, last
        sumtav = sumtav + atav(i)
        count = count + 1
      enddo
      tavd = sumtav / count
      state%atmosphere%tavd = tavd   ! [SS-GR-ATM A5.5] runtime dual-write

      ! Daily radiation (J/m2/d) and atmospheric demand (cm/d)
      rad = 0.d0
      at_atmdem = 0.d0
      do i = 1,nmetdetail
        rad = rad + arad(i)
        at_atmdem = at_atmdem + tpot(i)
      enddo

      ! Fluxes of current time step (start of the day)
      at_ptra = tpot(1)
      at_peva = epot(1)
      state%atmosphere%graidt  = grain(1)
      state%atmosphere%nraidt  = nrain(1)
      state%atmosphere%aintcdt = state%atmosphere%graidt - state%atmosphere%nraidt  ! aintcdt involves ONLY interception of RAIN

    endif

    end associate  ! at_peva/at_ptra/at_atmdem/at_pevaday/at_ptraday + tc_daynr/tc_t/tc_dt/tc_flmetdetail/tc_fletsine [TC-9]

  end subroutine ProcessMeteoDay

end module meteo_mod










