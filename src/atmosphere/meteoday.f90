! File VersionID:
!   $Id: meteoday.f90 372 2018-03-13 10:01:20Z heine003 $
!
!     This file contains the following subroutines, in order of calling:
!     1. MeteoDay          : main routine                       ; called in SWAP
!     2. ResetMetFlx       : resets meteorological fluxes       ; called in MeteoDay
!     3. ProcessMeteoDays  : processes meteo input data per day ; called in MeteoDay
!     4-6.                 : calculate interception according to:
!     4. VonHHBraden         - Von Hoyningen-Hune and Braden    ; called in ProcessMeteoDays (optional)
!     5. Gash                - Gash                             ; called in ProcessMeteoDays (optional)
!     6. Ruttervw            - Rutter (adapted)                 ; called in ProcessMeteoDays (optional)
!     7. DivIntercep       : divides interception into rain part
!                                              & sprinkling part; called in ProcessMeteoDays
!     8. Reduceva          : calc. reduction of soil evaporation; called in ProcessMeteoDays, ProcessMeteoTsteps, ETSine
! //[ ] write a module around the CNmethod, leave it open as runoff_mod for future additions to runoff calculation (e.g. Green-Ampt, TOPMODEL, etc.)
! //[ ] add more documentation to the CNmethod subroutine and runoff_mod.
! New dependency hierchy in the meteoday.f90 file:
! Level 1 (No dependencies on other meteo modules):
!   ├─ meteo_process_mod
!   ├─ interception_mod
!   ├─ et_mod
!   └─ runoff_mod

! Level 2 (Uses Level 1):
!   └─ precipitation_mod

! Level 3 (Orchestrator):
!   └─ meteo_mod (uses all above)

module MeteoVars
!  for availability in both ReadMeteoDay and ProcessMeteoDay
   integer   count,first,i,irecord,last,ndayparts
   real(8)   arain(96),awind(96),restint,interc,Edirectpond
   real(8)   aintc,dttp,eintc,etr,gctp,hum,netrainflux,rainflux
   real(8)   sumtav,svp,wfrac,win,Edirect,Tdirect,Tdirectwet
end module MeteoVars

module runoff_mod
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
   use soilhydraulics_utils, only: watcon
   implicit none
   private
   
   public :: CNmethod

   ! Physical constants
   real(8), parameter :: DEPTH_10CM = 10.0d0           !! Reference depth for moisture correction (cm)
   real(8), parameter :: H_FIELD_CAPACITY = -100.0d0   !! Pressure head at field capacity (cm)
   real(8), parameter :: H_WILTING_POINT = -16000.0d0  !! Pressure head at wilting point (cm)
   real(8), parameter :: IA_RATIO = 0.2d0              !! Initial abstraction ratio (dimensionless)

contains

  subroutine CNmethod(Itask)
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
    use variables, only: CNref, CNdry, CNwet, ThetaRef, theta, nraidt, Runoff_CN, zbotcp, dz, numnod, t1900, wc_cor, iCNtab, CNtimTAB, CNrefTAB, melt, wc10, &
                        nod10_cn, icn_atm, z10_cn
    use soilhydraulics_utils, only: watcon
    implicit none
    ! global
    integer, intent(in)  :: Itask
    ! local
    integer              :: i
    real(8)              :: wc1, wc2, CN, S, Ia
    ! Note: Nod10, iCN, Z10 are now module-level in variables.f90 as nod10_cn, icn_atm, z10_cn

    select case (Itask)
    ! initialization; calculate and store some constants
    case (1)

      icn_atm = 0
      ! check if times in CNtimeTAB are in ascending order
      ! set initial position in CNtimTAB
      do i = 2, iCNtab
          if (CNtimTAB(i) < CNtimTAB(i-1)) call fatalerr ('CNmethod', 'CNtimTAB not in ascending order')
          if (t1900 >= CNtimTAB(i-1) .and. t1900 < CNtimTAB(i)) icn_atm = i-1
      end do
      ! error if start time t1900 not in CNtimTAB
      if (icn_atm == 0) call fatalerr ('CNmethod', 'Start time of simulation not present in CNtimTAB')
      
    !  to be replaced by average for layer 0-10 cm
      do i = 1, numnod
          if (zbotcp(i) < -DEPTH_10CM) then
            nod10_cn = i-1
            z10_cn = -zbotcp(nod10_cn)
            exit
          end if
      end do
      ThetaRef = 0.0d0
      do i = 1, nod10_cn
          if (wc_cor == 1) then
            wc1 = watcon(i,H_FIELD_CAPACITY)
            wc2 = watcon(i,H_WILTING_POINT)
            ThetaRef = ThetaRef + (wc1+wc2)*0.5d0*dz(i)
          else if (wc_cor == 2) then
            wc1 = watcon(i,0.0d0)
            wc2 = watcon(i,H_WILTING_POINT)
            ThetaRef = ThetaRef + (wc1+wc2)*0.5d0*dz(i)
          end if
      end do
      ThetaRef = ThetaRef/z10_cn
      
      !!!t1900_old = int(t1900) - 1 ! for testing intermediate output
      
    ! dynamic part: calculate runoff   
    case (2)
      
      ! see if t1900 has moved ahead in CNtimTAB; icn_atm can never exceed last entry
      !  if (icn_atm < iCNtab .and. t1900 >= CNtimTAB(icn_atm+1)) icn_atm = icn_atm + 1
      ! Update position in CN time series if time has advanced (do while is more efficient if time steps are large and CN time series is long)
      do while (icn_atm < iCNtab .and. t1900 >= CNtimTAB(icn_atm + 1))
        icn_atm = icn_atm + 1
      end do
      CNref = CNrefTAB(icn_atm)
      CN    = CNref
      CNdry =  4.2d0*CNref/(10.0d0-0.058d0*CNref)
      CNwet = 23.0d0*CNref/(10.0d0+0.13d0*CNref)

      if (wc_cor > 0) then
          wc10 = 0.0d0
          do i = 1, nod10_cn
            wc10 = wc10 + theta(i)*dz(i)
          end do
          wc10 = wc10/z10_cn
          if (wc10 < ThetaRef) then
            CN = CNdry + wc10/ThetaRef*(CNref-CNdry)
          else
            CN = CNref + (wc10-ThetaRef)/ThetaRef*(CNwet-CNref)
          end if
      end if
      S = 2540d0/CN-25.4d0    ! in cm
      Ia = 0.2d0*S
    !   Ia = 0.3d0*S
      if (nraidt+melt > Ia) then
          Runoff_CN = (nraidt+melt-Ia)**2/(nraidt+melt-Ia+S)
      else
          Runoff_CN = 0.0d0
      end if
      
      ! for testing intermediate output
      !!!if (int(t1900) > t1900_old) then
      !!!   write (123, '(7F20.6)') t1900, nraidt, runoff_cn, cn, wc10, thetaref, melt
      !!!   t1900_old = int(t1900)
      !!!end if
      
    case default
      call fatalerr ('CNmethod', 'Illegal Itask option')
    end select

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
   implicit none
   private
   
   public :: ReadMeteoDay
   public :: ResetMetFlx

contains

  subroutine ReadMeteoDay
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
      use variables
      use MeteoVars
      use precipitation_mod, only: PartitionPrecipitation
      implicit none

    ! --- local
    character(len=11)  detdate
    character(len=3)   ext
    character(len=200) filnam
    character(len=300) messag

    ! ----------------------------------------------------------------------

    call ResetMetFlx ()

    ! 1: Check whether meteo data are available of today; pass on weather of today
    ! 1.0 Daily Meteo 0000000000000000000000000000000000000000000000000000000 Daily Meteo
    !
      if (swmetdetail.eq.0) then

    !   - Check availability of meteo data of today
        if (daymeteo.lt.daynrfirst .or. daymeteo.gt.daynrlast) then
          messag ='In meteo file no meteo data are'//                   &
     &    ' available for '//date//'. First adapt meteo file!'
          call fatalerr ('meteo',messag)
        end if

    !   - Pass on weather values of today
        rad  = arad(daymeteo+1-daynrfirst)
        tmn  = atmn(daymeteo+1-daynrfirst)
        tmx  = atmx(daymeteo+1-daynrfirst)
        hum  = ahum(daymeteo+1-daynrfirst)
        win  = awin(daymeteo+1-daynrfirst)
        grai = arai(daymeteo+1-daynrfirst)
        etr  = aetr(daymeteo+1-daynrfirst)

    !   - If hum is missing or tav cannot be calculated: set rh at -99.0
        rh = 1.0d0
        if (hum.lt.-98.0d0 .or. tmn.lt.-98.0d0 .or. tmx.lt.-98.0d0)     &
     &      rh=-99.0d0

    !   - Calculate 24h average temperature
        tav = (tmx+tmn)*0.5d0
    !   - Calculate average day temperature
        tavd = (tmx+tav)*0.5d0

        if (rh.ge.-98.0d0) then
    !   - Calculate saturated vapour pressure [kpa]
        svp = 0.3055d0*(dexp(17.27d0*tmn/(tmn+237.3d0)) +               &
     &                  dexp(17.27d0*tmx/(tmx+237.3d0)))
    !   - Calculate relative humidity [fraction]
          rh = min(hum/svp,1.0d0)
        endif

    !   - CFO file for PEARL: save meteo variables of today for output
        out_rad = real(rad)
        out_tmn = real(tmn)
        out_tmx = real(tmx)
        out_hum = real(hum)
        out_win = real(win)
        out_etr = real(etr)*0.001
        if (swrain.eq.2) then
          out_wet = real(wet(daymeteo+1-daynrfirst))
        else
          out_wet = -1.0
        endif
    !
    ! end 1 Daily Meteo 00000000000000000000000000000000000000000000000000000 Daily Meteo
    !!
    ! 1.1 Detailed Meteo 1111111111111111111111111111111111111111111111111111 Detailed Meteo
    !
      elseif (swmetdetail.eq.1) then

    !   - Check availability of meteo data of today
    !     + compose filename meteorological file for use in warnings
        write (ext,'(i3.3)') mod(yearmeteo,1000)
        filnam = trim(pathatm)//trim(metfil)//'.'//trim(ext)

        do i = 1, nmetdetail
          irectotal = irectotal + 1
          if (i .ne. detrecord(irectotal)) then
            messag='In meteo file '//trim(filnam)//' record number(s)'//&
     &      ' are not correct at '//date//'. First adapt meteo file!'
            call fatalerr ('meteo',messag)
          end if
          call dtdpst('year-month-day',                               &
     &                 dettime(irectotal)+0.1d0,detdate)
          call dtdpst('year-month-day',t1900+0.1d0,date)
          if (detdate .ne. date) then
            messag ='In meteo file '//trim(filnam)//' the amount of '// &
     &      'records deviate near '//date//'. First adapt meteo file!'
            call fatalerr ('meteo',messag)
          end if

    !   - Pass on weather records of today
          arad(i)  = detrad(irectotal)
          ahum(i)  = dethum(irectotal)
          atav(i)  = dettav(irectotal)
          awind(i) = detwind(irectotal)
          arain(i) = detrain(irectotal) * 0.1d0 ! convert from mm to cm
        enddo
      endif
    !
    ! end 1 Detailed Meteo 11111111111111111111111111111111111111111111111111 Detailed Meteo
    ! end 1.

   ! Call precipitation partitioning module
   call PartitionPrecipitation(swmetdetail, swsnow, tav, TePrRain, TePrSnow, &
                               ssnow, nmetdetail, arain, grai, gsnow, snrai, &
                               fprecnosnow, restint)
      return
  end subroutine ReadMeteoDay


  subroutine ResetMetFlx ()
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
      use variables, only: flzerointr,flzerocumu,caintc,cgrai,cnrai,igrai,inrai,iprec
      implicit none

    ! --- local

    ! --- reset cumulative intermediate fluxes
      if (flzerointr) then
        iprec = 0.0d0
        igrai = 0.0d0
        inrai = 0.0d0
      endif

    ! --- reset cumulative meteorological fluxes
      if (flzerocumu) then
        cgrai = 0.0d0
        cnrai = 0.0d0
        caintc = 0.0d0
      endif

      return
  end subroutine ResetMetFlx

end module meteo_process_mod

module meteo_mod
  !> Main meteorological processing coordinator
  ! Subroutines:
  ! - ProcessMeteoDay  ! Main daily processing orchestrator
  
  use meteo_process_mod, only: ReadMeteoDay, ResetMetFlx
  use interception_mod, only: VonHHBraden, Gash, ruttervw, DivIntercep
  use et_mod, only: PenMon, reduceva
  use runoff_mod, only: CNmethod
  
  implicit none
  
  private
  public :: ProcessMeteoDay

contains

  !> Process daily meteorological data and calculate evapotranspiration components
  !>
  !> Main orchestrator for daily meteorological processing. Performs:
  !> - Interception calculations (Von Hoyningen-Hune & Braden, Gash, or Rutter methods)
  !> - Reference evapotranspiration calculation or input processing
  !> - Potential soil evaporation and transpiration partitioning
  !> - Wet fraction calculations for crop canopy
  !> - Sub-daily timestep handling for detailed meteorology
  !>
  !> ## Processing Steps
  !> 1. Calculate interception and net rainfall/irrigation
  !> 2. Loop over day parts (1 for daily, nmetdetail for sub-daily)
  !> 3. Calculate ET0, EW0, ES0 via Penman-Monteith or use specified values
  !> 4. Apply interception for adapted Rutter model (if selected)
  !> 5. Calculate wet fraction of canopy
  !> 6. Partition into potential soil evaporation and transpiration
  !> 7. Apply atmospheric CO2 corrections if enabled
  !> 8. Aggregate results for detailed meteorology
  !>
  !> @note
  !> Last modified: February 2014
  !> Supports both daily (swmetdetail=0) and detailed (swmetdetail=1) meteorology
  !> Uses module variables from Variables and MeteoVars
  !> @endnote
  subroutine ProcessMeteoDay
    use Variables
    use MeteoVars
    use array_utils, only: afgen
    use runoff_mod, only: CNmethod
    use interception_mod, only: VonHHBraden, Gash, ruttervw, msw1eic, DivIntercep
    implicit none
    include  'params.fi'
    real(8)  rcs
    data     rcs/0.15d0/

    ! === Section 3: Interception calculations ===

    ! Calculation of interception and net rain & net irrigation depth [cm]
    if ((lai .lt. 1.d-3) .or. (grai+gird .lt. 1.d-5) .or. &
        (swinter.eq.0) .or. (gsnow.gt.0.0d0) .or.(ssnow.gt.0.0d0)) then

      ! No vegetation, rainfall/irrigation or interception calculation
      aintc = 0.d0

    else if (swinter .eq. 1) then
      ! Calculate interception, method Von Hoyningen-Hune and Braden
      call VonHHBraden (aintc)
    else if (swinter .eq. 2) then
      ! Calculate interception, method Gash (1995)
      call Gash (aintc)
    end if

    ! Divide interception into rain part and irrigation part and
    ! calculate net rain (nraida) and net sprinkling irrigation (nird)
    if (swinter.ne.3) &
      call DivIntercep (aintc)

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
        call PenMon (logf,swscre,daynr,lat,alt,Altw,angstroma, &
                     angstromb,rcs,rad,tav,hum,win,rsc,es0,et0,ew0,swcf,ch, &
                     flCropEmergence,daylp,flmetdetail,irecord, &
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

      ! === Section 5: Interception option NHI (adapted Rutter model) ===

      if (swinter .eq. 3) then
        ! Set parameter values
        if (swmetdetail.eq.0) then  ! if swmetdetail = 0, siccapact is set in cropgrowth module
          if (croptype(icrop).eq.1 .and. swgc.eq.2) then
            gctp  = gc
          else
            gctp  = 1.0d0 - dexp(-1.0d0*kdir*kdif*lai)
            if (gctp .lt. 1.0d-5) siccapact=0.
          endif
          dttp = 1.0d0  ! value of 1 d required for the daily meteo option
        elseif (swmetdetail.eq.1) then
          siccapact = afgen(siccaptb,(2*magrs),t)
          if (croptype(icrop).eq.1 .and. swgc.eq.2) then
            gctp  = gc
          elseif (croptype(icrop).eq.1 .and. swgc.eq.1) then
            gctp  = 1.0d0 - dexp(-1.0d0*kdir*kdif*lai)
          endif
          if (gctp .lt. 1.0d-5) siccapact=0.
          dttp = dt
        endif

        ! Calculate interception, method Rutter
        call ruttervw (gctp,aintc,eintc)

        ! Divide interception into rain part and irrigation part and
        ! calculate net rain (nraida) and net sprinkling irrigation (nird)
        call DivIntercep(aintc)
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
        if (grai .lt. 1.0d-12) then
          interc = 0.0d0
          wfrac  = 0.0d0
        else
          interc = restint + aintc * arain(irecord) / grai
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
      peva = max(0.0d0, (es0*dexp(-1.0d0*kdir*kdif*lai)*0.1d0))
      if (swcf.ne.3 .or. (swmetdetail.eq.0 .and. swinter.ne.3)) &
        peva = max(0.0d0,(1.0d0-wfrac)*peva)

      ! Alternative for peva (simple model, soil cover fraction specified)
      if (flCropCalendar .and. .not.flCropHarvest) then
        if (croptype(icrop).eq.1 .and. swgc.eq.2) then
          peva = (1.0d0-gc)*es0*0.1d0
          if (swcf.ne.3 .or. (swmetdetail.eq.0 .and. swinter.ne.3)) &
            peva = (1.0d0-wfrac)*peva
        endif
      endif

      ! Adapt peva in case of ponding
      if (pond .gt. 1.0d-10) then
        if (SwETr.eq.0 .and. es0.gt.1.0d-8) then
          peva = ew0/es0 * peva
        elseif (es0.gt.1.0d-8) then
          if (swcfbs .eq. 1 .and. cfbs .gt. small) then
            peva = cfevappond * peva / cfbs
          else
            peva = cfevappond * peva
          endif
        endif
      endif

      ! Potential soil evaporation [cm/d] according to PMdirect
      if (swdivide .eq. 1) then
        if (pond .gt. 1.0d-10) then
          peva = Edirectpond*0.1d0
        else
          peva = Edirect*0.1d0
        endif
      endif

      ! Potential transpiration (ptra) [cm/d]
      if (swcf .ne. 3) then
        ptra = ((1.0d0-wfrac)*et0-peva*10.0d0)*0.1d0
      else
        ptra = (1.0d0-wfrac)*et0*0.1d0
      endif
      ptra = max(ptra,(1.01d0*nihil))

      ! Potential transpiration [cm/d] according to PMdirect
      if (swdivide .eq. 1) then
        ptra = (1.0d0-wfrac) * Tdirect * 0.1d0
        ptra = max(ptra,(1.01d0*nihil))
      endif

      ! Correction of potential transpiration as a function of atmospheric CO2 concentration
      if (flCO2 .and. flCropEmergence) then
        ptra = fco2tra * ptra
      endif

      ! === Section 8: Results for detailed weather records ===

      ! Result of detailed weather records (cm/d)
      if (swmetdetail.eq.1) then
        tpot(irecord) = ptra
        epot(irecord) = peva
        if (grai .lt. 1.0d-12) then
          grain(irecord) = 0.0d0
          nrain(irecord) = 0.0d0
        else
          grain(irecord) = arain(irecord) / metperiod
          nrain(irecord) = arain(irecord) / metperiod * nraida /grai
        endif
      endif

1000  continue
    ! === END LOOP over dayparts ===

    ! === Section 9: Actual daily rain/snow fluxes and soil evaporation for Daily Meteo ===

    if (swmetdetail.eq.0) then

      ! Finterception: ratio net / gross rain flux; net rainflux = gross - interception
      if (grai.gt.1.d-5) then
        ! finterception is exclusively meant for dividing rain flux into interception part
        ! and net rain part; not for sprinkler irrigation!
        finterception = nraida / grai
        if (aintc.lt.1.0d-5) finterception = 1.0d0
      else
        finterception = 1.0d0
      endif

      ! In case of daily precipitation sum: set actual gross and net rainflux,
      ! and interception on TIMESTEP basis
      if (swrain.eq.0) then
        rainflux    = fprecnosnow * grai
        netrainflux = finterception * rainflux
        graidt  = rainflux
        nraidt  = netrainflux
        if (swuseCN == 1) then
          call CNmethod(2)
          nraidt = nraidt - Runoff_CN
        end if
        aintcdt = rainflux - netrainflux  ! aintcdt involves ONLY interception of RAIN
      endif

      ! Soil evaporation rate of today
      if (.not. fletsine) then
        call reduceva (1,nraida)
      endif

      ! Save daily potential values for use in ETSine
      ptraday = ptra
      pevaday = peva

      ! Calculate atmospheric demand [cm]
      atmdem = et0*0.1d0

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

      ! Daily radiation (J/m2/d) and atmospheric demand (cm/d)
      rad = 0.d0
      atmdem = 0.d0
      do i = 1,nmetdetail
        rad = rad + arad(i)
        atmdem = atmdem + tpot(i)
      enddo

      ! Fluxes of current time step (start of the day)
      ptra = tpot(1)
      peva = epot(1)
      graidt = grain(1)
      nraidt = nrain(1)
      aintcdt = graidt - nraidt    ! aintcdt involves ONLY interception of RAIN

    endif

  end subroutine ProcessMeteoDay
end module meteo_mod










