!> Evaporation and evapotranspiration calculations
!>
!! This module provides routines for calculating potential evaporation and 
!! evapotranspiration rates using the Penman-Monteith approach, and for
!! reducing soil evaporation based on atmospheric conditions.
!>
!! @note
!! Subroutines:
!! - PenMon: Penman-Monteith ET calculation
!! - reduceva: Soil evaporation reduction (Black/Boesten-Stroosnijder)
!! @endnote
module et_mod
   use error_mod, only: fatalerr_collected
   use swap_state_mod, only: swap_state_t
    implicit none
    private

  public :: PenMon, PenMon_calc, reduceva
contains
  !> Pure Penman-Monteith calculation (no I/O, fully deterministic)
  !!
  !! Performs all ET calculations without side effects. This pure version
  !! enables compiler optimizations, parallelization, and easier testing.
  !!
  !! @note All warnings are returned as optional status codes instead of being logged
  pure subroutine PenMon_calc(daynr, lat, alt, altw, a, b, rcs, rad, &
                              tav, hum, win, rsc, swcf, ch, flCropEmergence, &
                              daylp, flmetdetail, irecord, nmetdetail, &
                              albedo, tmn, tmx, rsw, difpp, dsinbe, atmtr, &
                              rsoil, swdivide, kdif, kdir, lai, &
                              es0, et0, ew0, Edirect, Tdirect, Tdirectwet, &
                              Edirectpond, warning_code)
      use swap_constants, only: vlarge, small, KARMAN_CONSTANT, &
                                GRASS_HEIGHT_CM, MEASUREMENT_HEIGHT_CM, &
                                BARE_SOIL_HEIGHT_CM, PI, ALBEDO_PONDING
      implicit none
      
      ! Input parameters
      integer, intent(in) :: daynr
        !! Day number (January 1st = 1) [-]
      real(8), intent(in) :: lat
        !! Latitude [deg, decimal degrees, N=+, S=-]
      real(8), intent(in) :: alt
        !! Altitude above mean sea level [m]
      real(8), intent(in) :: altw
        !! Altitude of wind speed measurement [m]
      real(8), intent(in) :: a, b
        !! Angstrom coefficients [-]
      real(8), intent(in) :: rcs
        !! Reflection coefficient soil [-]
      real(8), intent(in) :: rad
        !! Incoming short wave radiation [J/m2/d]
      real(8), intent(in) :: tav
        !! Average temperature (24 hour) [C]
      real(8), intent(in) :: hum
        !! Vapour pressure [kPa]
      real(8), intent(in) :: win
        !! Wind speed at 2 m height [m/s]
      real(8), intent(in) :: rsc
        !! Minimum canopy resistance of dry crop [s/m]
      integer, intent(in) :: swcf
        !! Switch: use crop factor (=1) or crop height (=2)
      real(8), intent(in) :: ch
        !! Crop height [cm]
      logical, intent(in) :: flCropEmergence
        !! Crop emergence flag
      real(8), intent(in) :: daylp
        !! Day length parameter
      logical, intent(in) :: flmetdetail
        !! Flag for detailed meteorological data
      integer, intent(in) :: irecord
        !! Current record number
      integer, intent(in) :: nmetdetail
        !! Number of detailed meteorological records
      real(8), intent(in) :: albedo
        !! Reflection coefficient crop [-]
      real(8), intent(in) :: tmn, tmx
        !! Minimum and maximum temperature [C]
      real(8), intent(in) :: rsw
        !! Canopy resistance of intercepted water [s/m]
      real(8), intent(in) :: difpp
        !! Diffuse radiation parameter
      real(8), intent(in) :: dsinbe
        !! Solar elevation parameter
      real(8), intent(in) :: atmtr
        !! Daily atmospheric transmission
      real(8), intent(in) :: rsoil
        !! Soil resistance in PMdirect [s/m]
      integer, intent(in) :: swdivide
        !! Switch for direct partitioning method (1=PMdirect)
      real(8), intent(in) :: kdif, kdir
        !! Diffuse and direct light extinction coefficients [-]
      real(8), intent(in) :: lai
        !! Leaf area index [-]
      
      ! Output variables
      real(8), intent(out) :: es0
        !! Potential evaporation rate from a wet bare soil [mm/d]
      real(8), intent(out) :: et0
        !! Potential transpiration rate from a dry crop [mm/d]
      real(8), intent(out) :: ew0
        !! Potential transpiration rate from a wet crop [mm/d]
      real(8), intent(out) :: Edirect
        !! Direct evaporation [mm/d]
      real(8), intent(out) :: Tdirect
        !! Direct transpiration [mm/d]
      real(8), intent(out) :: Tdirectwet
        !! Direct transpiration from wet crop [mm/d]
      real(8), intent(out) :: Edirectpond
        !! Direct evaporation from ponding layer [mm/d]
      integer, intent(out), optional :: warning_code
        !! Warning code: 0=none, 1=polar circle 0hrs, 2=polar circle 24hrs
      
      ! Local variables - atmospheric
      real(8) :: lambda, delta, ea, ed, vpd, gamma, palt, rho, cp
      real(8) :: tavk, tmnk, tmxk, tkv
      
      ! Local variables - radiation
      real(8) :: rns, rnc, rnw, rnp, rnl, relssd
      real(8) :: gs, gc, gw
      real(8) :: sinld, cosld, dayl, radial, dec, aob
      real(8) :: sunrise, sunset, startrec, endrec
      
      ! Local variables - aerodynamic
      real(8) :: chplant, zmeasw
      real(8) :: ud, rac, ras, raw, rss
      real(8) :: zm, zh, d, zom, zoh
      real(8) :: dgrass, zomgrass, zact, dact, zomact, fact, fmeas
      
      ! Local variables - Penman-Monteith terms
      real(8) :: gammos, gammoc, gammow, gammop
      real(8) :: etaers, etaerc, etaerw, etaerp
      real(8) :: etrads, etradc, etradw, etradp
      
      ! Local variables - PMdirect
      real(8) :: vcover, laieff
      
      ! Initialize outputs
      es0 = 0.0d0
      et0 = 0.0d0
      ew0 = 0.0d0
      Edirect = 0.0d0
      Tdirect = 0.0d0
      Tdirectwet = 0.0d0
      Edirectpond = 0.0d0
      if (present(warning_code)) warning_code = 0
      
      ! ========================================================================
      ! 1. PREPROCESSING: Unit conversions and crop height determination
      ! ========================================================================
      
      zmeasw = 100.0d0 * altw  ! Convert wind measurement height to cm
      
      ! Determine effective crop height
      if (.not. flCropEmergence) then
          chplant = GRASS_HEIGHT_CM
      else
          if (swcf == 1 .or. swcf == 3) then
              chplant = GRASS_HEIGHT_CM
          else
              chplant = max(ch, 0.1d0)
          endif
      endif
      
      ! ========================================================================
      ! 2. ATMOSPHERIC PROPERTIES
      ! ========================================================================
      
      ! Temperature conversions [K]
      tavk = tav + 273.15d0
      tmnk = tmn + 273.15d0
      tmxk = tmx + 273.15d0
      
      ! Atmospheric pressure at elevation [kPa]
      palt = 101.3d0 * ((tavk - 0.0065d0*alt) / tavk)**5.26d0
      
      ! Latent heat of vaporization [MJ/kg]
      lambda = 2.501d0 - 0.002361d0*tav
      
      ! Saturation vapour pressure [kPa]
      if (flmetdetail) then
          ea = 0.611d0 * exp(17.27d0*tav / (tav + 237.3d0))
      else
          ea = 0.3055d0 * (exp(17.27d0*tmn / (tmn + 237.3d0)) + &
                          exp(17.27d0*tmx / (tmx + 237.3d0)))
      endif
      
      ! Measured vapour pressure (capped at saturation)
      ed = min(hum, ea)
      
      ! Vapour pressure deficit [kPa]
      vpd = ea - ed
      
      ! Slope of vapour pressure curve [kPa/C]
      delta = 4098.0d0 * ea / (tav + 237.3d0)**2
      
      ! Psychrometric constant [kPa/C]
      gamma = 0.00163d0 * palt / lambda
      
      ! Atmospheric density [kg/m3]
      tkv = tavk / (1.0d0 - 0.378d0*ed/palt)
      rho = 3.486d0 * palt / tkv
      
      ! Specific heat of moist air [kJ/kg/C]
      cp = 622.0d0 * gamma * lambda / palt
      
      ! ========================================================================
      ! 3. WIND SPEED AND AERODYNAMIC RESISTANCE
      ! ========================================================================
      
      ! Day wind speed [m/s], avoid zero
      ud = max(win, 0.0001d0)
      
      ! Adjust wind speed for height differences
      if (chplant > MEASUREMENT_HEIGHT_CM .or. zmeasw > MEASUREMENT_HEIGHT_CM) then
          dgrass = 2.0d0/3.0d0 * GRASS_HEIGHT_CM
          zomgrass = 0.123d0 * GRASS_HEIGHT_CM
          fmeas = log((1.0d4 - dgrass) / zomgrass) / &
                  log((zmeasw - dgrass) / zomgrass)
          
          zact = max(chplant, 200.0d0)
          dact = 2.0d0/3.0d0 * chplant
          zomact = 0.123d0 * chplant
          fact = log((zact - dact) / zomact) / log((1.0d4 - dact) / zomact)
          
          ud = ud * fact * fmeas
      endif
      
      ! Aerodynamic parameters for crop
      zm = max(chplant, MEASUREMENT_HEIGHT_CM)
      zh = zm
      d = 2.0d0/3.0d0 * chplant
      zom = 0.123d0 * chplant
      zoh = 0.1d0 * zom
      
      ! Aerodynamic resistance for crop (dry and wet) [s/m]
      rac = log((zm - d)/zom) * log((zh - d)/zoh) / KARMAN_CONSTANT**2 / ud
      raw = rac
      
      ! Aerodynamic resistance for bare soil [s/m]
      d = 2.0d0/3.0d0 * BARE_SOIL_HEIGHT_CM
      zom = 0.123d0 * BARE_SOIL_HEIGHT_CM
      zoh = 0.1d0 * zom
      ras = log((zm - d)/zom) * log((zh - d)/zoh) / KARMAN_CONSTANT**2 / ud
      
      ! Surface resistance of soil [s/m]
      if (swdivide == 1) then
          rss = rsoil  ! PMdirect partitioning
      else
          rss = 0.0d0
      endif
      
      ! Modified psychrometric constants [kPa/C]
      gammos = gamma * (1.0d0 + rss/ras)
      gammoc = gamma * (1.0d0 + rsc/rac)
      gammow = gamma * (1.0d0 + rsw/raw)
      
      ! ========================================================================
      ! 4. RADIATION CALCULATIONS
      ! ========================================================================
      
      ! Net shortwave radiation [MJ/m2/d]
      rns = (1.0d0 - rcs) * rad / 1.0d6
      rnc = (1.0d0 - albedo) * rad / 1.0d6
      rnw = (1.0d0 - albedo) * rad / 1.0d6
      rnp = (1.0d0 - ALBEDO_PONDING) * rad / 1.0d6
      
      ! Extraterrestrial radiation and daylength
      if (flmetdetail) then
          ! Sub-daily: calculate astronomical parameters
          radial = PI / 180.0d0
          dec = -asin(sin(23.45d0*radial) * &
                      cos(2.0d0*PI*dble(daynr+10) / 365.0d0))
          
          sinld = sin(radial*lat) * sin(dec)
          cosld = cos(radial*lat) * cos(dec)
          aob = sinld / cosld
          
          ! Daylength calculation with polar circle handling
          if (aob < -1.0d0) then
              dayl = 0.0d0
              if (present(warning_code)) warning_code = 1  ! Polar circle, 0 hours
          else if (aob > 1.0d0) then
              dayl = 24.0d0
              if (present(warning_code)) warning_code = 2  ! Polar circle, 24 hours
          else
              dayl = 12.0d0 * (1.0d0 + 2.0d0*asin(aob)/PI)
          endif
          
          sunrise = 0.5d0 - dayl / 48.0d0
          sunset = 0.5d0 + dayl / 48.0d0
          startrec = dble(real(irecord-1) / real(nmetdetail))
          endrec = dble(real(irecord) / real(nmetdetail))
      else
          ! Daily: use provided parameters
          ! Note: sinld, cosld would come from astro() call in wrapper
          sinld = dsinbe
          cosld = sqrt(max(0.0d0, 1.0d0 - sinld**2))
          startrec = 0.0d0
          endrec = 1.0d0
          sunrise = 0.0d0
          sunset = 1.0d0
      endif
      
      ! Net longwave radiation [MJ/m2/d]
      relssd = max(min((atmtr - a) / b, 1.0d0), 0.0d0)
      rnl = 4.9d-9 * 0.5d0 * (tmxk**4 + tmnk**4) * &
            (0.34d0 - 0.14d0*sqrt(ed)) * (0.1d0 + 0.9d0*relssd)
      
      ! Soil heat flux [MJ/m2/d]
      if (flmetdetail) then
          if ((startrec + endrec)/2.0d0 > sunrise .and. &
              (startrec + endrec)/2.0d0 < sunset) then
              ! Daytime
              gs = 0.1d0 * (rns - rnl)
              gc = 0.1d0 * (rnc - rnl)
              gw = 0.1d0 * (rnw - rnl)
          else
              ! Nighttime
              gs = 0.5d0 * (rns - rnl)
              gc = 0.5d0 * (rnc - rnl)
              gw = 0.5d0 * (rnw - rnl)
          endif
      else
          ! Daily: negligible net flux
          gs = 0.0d0
          gc = 0.0d0
          gw = 0.0d0
      endif
      
      ! ========================================================================
      ! 5. PENMAN-MONTEITH EQUATION - STANDARD METHOD
      ! ========================================================================
      
      ! Aerodynamic term [mm/d]
      etaers = (86.4d0/lambda) * (1.0d0/(delta + gammos)) * (rho*cp*vpd/ras)
      etaerc = (86.4d0/lambda) * (1.0d0/(delta + gammoc)) * (rho*cp*vpd/rac)
      etaerw = (86.4d0/lambda) * (1.0d0/(delta + gammow)) * (rho*cp*vpd/raw)
      
      ! Radiation term [mm/d]
      etrads = delta/(delta + gammos) * (rns - rnl - gs) / lambda
      etradc = delta/(delta + gammoc) * (rnc - rnl - gc) / lambda
      etradw = delta/(delta + gammow) * (rnw - rnl - gw) / lambda
      
      ! Total potential rates [mm/d]
      es0 = max(0.0d0, etaers + etrads)
      et0 = max(0.0d0, etaerc + etradc)
      ew0 = max(0.0d0, etaerw + etradw)
      
      ! ========================================================================
      ! 6. PENMAN-MONTEITH DIRECT PARTITIONING (PMdirect)
      ! ========================================================================
      
      if (swdivide == 1) then
          ! Vegetation cover fraction
          vcover = 1.0d0 - exp(-kdif * kdir * lai)
          
          ! Adjust aerodynamic resistances for partial cover
          if (vcover > 1.0d-6) then
              rac = rac / vcover
          else
              rac = 1.0d12
          endif
          raw = rac
          
          if ((1.0d0 - vcover) > 1.0d-6) then
              ras = ras / (1.0d0 - vcover)
          else
              ras = 1.0d12
          endif
          
          ! Effective LAI for resistance scaling
          laieff = lai / (0.3d0*lai + 1.2d0)
          
          ! Modified psychrometric constants with resistance scaling
          gammos = vlarge
          if (ras > small) gammos = gamma * (1.0d0 + rss/ras)
          
          gammoc = vlarge
          if ((rac*laieff) > small) gammoc = gamma * (1.0d0 + rsc/(rac*laieff))
          
          gammow = vlarge
          if ((raw*laieff) > small) gammow = gamma * (1.0d0 + rsw/(raw*laieff))
          
          gammop = vlarge
          if (ras > small) gammop = gamma
          
          ! Aerodynamic terms [mm/d]
          etaers = (86.4d0/lambda) * (1.0d0/(delta + gammos)) * (rho*cp*vpd/ras)
          etaerc = (86.4d0/lambda) * (1.0d0/(delta + gammoc)) * (rho*cp*vpd/rac)
          etaerw = (86.4d0/lambda) * (1.0d0/(delta + gammow)) * (rho*cp*vpd/raw)
          etaerp = (86.4d0/lambda) * (1.0d0/(delta + gammop)) * (rho*cp*vpd/ras)
          
          ! Radiation terms [mm/d] weighted by cover fraction
          etrads = delta/(delta + gammos) * (rns - rnl - gs) * (1.0d0 - vcover) / lambda
          etradc = delta/(delta + gammoc) * (rnc - rnl - gc) * vcover / lambda
          etradw = delta/(delta + gammow) * (rnw - rnl - gw) * vcover / lambda
          etradp = delta/(delta + gammop) * (rnp - rnl - gs) * (1.0d0 - vcover) / lambda
          
          ! Direct partitioned rates [mm/d]
          Edirect = max(0.0d0, etaers + etrads)
          Tdirect = max(0.0d0, etaerc + etradc)
          Tdirectwet = max(0.0d0, etaerw + etradw)
          Edirectpond = max(0.0d0, etaerp + etradp)
      endif

  end subroutine PenMon_calc

  !> Penman-Monteith evapotranspiration calculation (wrapper with I/O)
  !!
  !! This is a thin wrapper around PenMon_calc that handles validation,
  !! warnings, and calls to external routines (astro). Legacy interface
  !! is preserved for backward compatibility.
  subroutine PenMon(logf, swscre, daynr, lat, alt, altw, a, b, rcs, rad, &
                    tav, hum, win, rsc, es0, et0, ew0, swcf, ch, flCropEmergence, &
                    daylp, flmetdetail, irecord, nmetdetail, albedo, tmn, tmx, rsw, &
                    difpp, dsinbe, atmtr, Edirect, Tdirect, Tdirectwet, rsoil, &
                    swdivide, kdif, kdir, lai, Edirectpond)
      implicit none
      
      ! I/O parameters
      integer, intent(in) :: logf
        !! Internal number of logbook output file
      integer, intent(in) :: swscre
        !! Switch of screen display: 0=none, 1=summary, 2=daynumber
      
      ! All other parameters (same as PenMon_calc)
      integer, intent(in) :: daynr, swcf, irecord, nmetdetail, swdivide
      real(8), intent(in) :: lat, alt, altw, albedo, tmn, tmx, ch, difpp, dsinbe
      real(8), intent(in) :: a, b, rcs, rsc, rsw, hum, rad, tav, win, atmtr
      real(8), intent(in) :: rsoil, kdif, kdir, lai, daylp
      real(8), intent(out) :: es0, et0, ew0
      real(8), intent(out) :: Edirect, Tdirect, Tdirectwet, Edirectpond
      logical, intent(in) :: flCropEmergence, flmetdetail
      
      ! Local variables
      integer :: warning_code
      real(8) :: dayl, sinld, cosld
      character(len=200) :: messag
      
      ! Call astro() for daily radiation if needed
      if (.not. flmetdetail) then
          call astro(daynr, lat, rad, dayl, daylp, sinld, cosld, difpp, &
                    atmtr, dsinbe)
      endif
      
      ! Call pure calculation core
      call PenMon_calc(daynr, lat, alt, altw, a, b, rcs, rad, &
                      tav, hum, win, rsc, swcf, ch, flCropEmergence, &
                      daylp, flmetdetail, irecord, nmetdetail, &
                      albedo, tmn, tmx, rsw, difpp, dsinbe, atmtr, &
                      rsoil, swdivide, kdif, kdir, lai, &
                      es0, et0, ew0, Edirect, Tdirect, Tdirectwet, &
                      Edirectpond, warning_code)
      
      ! Handle warnings
      if (warning_code == 1) then
          messag = 'Warning: latitude above polar circle, daylength = 0hrs'
          call warn('Astro', messag, logf, swscre)
      else if (warning_code == 2) then
          messag = 'Warning: latitude within polar circle, daylength = 24hrs'
          call warn('Astro', messag, logf, swscre)
      endif

  end subroutine PenMon

    !> Black's evaporation reduction model
    !!
    !! Reduces soil evaporation based on square root of time since last rainfall.
    !! Formulation: E = β(√t - √(t-Δt)) / Δt
    !!
    subroutine black_reduction(nrai, nird, peva, cofred, rsigni, &
                               ldwet, empreva, dt, fldaystart, task)
        real(8), intent(in)    :: nrai
          !! Net rainfall [cm]
        real(8), intent(in)    :: nird
          !! Net irrigation [cm]
        real(8), intent(in)    :: peva
          !! Potential evaporation [cm/d]
        real(8), intent(in)    :: cofred
          !! Reduction coefficient β [cm/√d]
        real(8), intent(in)    :: rsigni
          !! Significant rainfall threshold [cm]
        real(8), intent(inout) :: ldwet
          !! Days since last significant rainfall [d]
        real(8), intent(out)   :: empreva
          !! Empirical potential evaporation [cm/d]
        real(8), intent(in)    :: dt
          !! Time step [d]
        logical, intent(in)    :: fldaystart
          !! Flag for start of day
        integer, intent(in)    :: task
          !! Task selector (1=daily, 2=timestep)
        
        ! Reset counter if significant rainfall occurred
        if (task == 1) then
            ! Daily: reset at any time significant rain occurs
            if ((nrai + nird) > rsigni) ldwet = 0.0d0
            ldwet = ldwet + 1.0d0
            
            ! Calculate reduction using daily formulation
            empreva = cofred * (sqrt(ldwet) - sqrt(ldwet - 1.0d0))
            
        else if (task == 2) then
            ! Sub-daily: only reset at start of day
            if (fldaystart .and. (nrai + nird) > rsigni) ldwet = 0.0d0
            
            ! Calculate reduction rate using sub-daily formulation
            empreva = (cofred * (sqrt(ldwet + dt) - sqrt(ldwet))) / dt
            
            ! Update state for next time step
            ldwet = ldwet + dt
        end if
        
        ! Limit to potential evaporation
        empreva = min(empreva, peva)
        
    end subroutine black_reduction


    !> Boesten-Stroosnijder evaporation reduction model
    !!
    !! Reduces soil evaporation based on cumulative deficit and actual evaporation.
    !! Uses two regimes:
    !! - Stage 1 (spev < β²): Actual evaporation = potential evaporation deficit
    !! - Stage 2 (spev ≥ β²): Actual evaporation = β√(deficit)
    !!
    subroutine boesten_stroosnijder_reduction(nrai, nird, peva, cofred, &
                                              spev, saev, empreva, dt)
        real(8), intent(in)    :: nrai
          !! Net rainfall [cm]
        real(8), intent(in)    :: nird
          !! Net irrigation [cm]
        real(8), intent(in)    :: peva
          !! Potential evaporation [cm/d]
        real(8), intent(in)    :: cofred
          !! Reduction coefficient β [cm/√d]
        real(8), intent(in)    :: dt
          !! Time step [d]
        real(8), intent(inout) :: spev
          !! Cumulative potential evaporation deficit [cm]
        real(8), intent(inout) :: saev
          !! Cumulative actual evaporation [cm]
        real(8), intent(out)   :: empreva
          !! Empirical potential evaporation [cm/d]
        
        real(8) :: water_input, deficit, saev_old
        real(8) :: cofred_squared
        
        water_input = nrai + nird
        cofred_squared = cofred**2
        
        if (water_input < peva) then
            ! Water input less than potential evaporation: deficit accumulates
            ! Increase deficit by unmet potential evaporation
            deficit = (peva - water_input) * dt
            spev = spev + deficit
            
            ! Calculate cumulative actual evaporation
            saev_old = saev
            if (spev < cofred_squared) then
                ! Stage 1: Energy-limited (actual = potential)
                saev = spev
            else
                ! Stage 2: Soil-limited (square root relationship)
                saev = cofred * sqrt(spev)
            end if
            
            ! Evaporation rate for this time step
            empreva = (water_input * dt + saev - saev_old) / dt
            
        else
            ! Water input exceeds potential evaporation: reduce deficit
            empreva = peva
            
            ! Reduce cumulative actual evaporation by excess water
            saev = max(0.0d0, saev - (water_input - peva) * dt)
            
            ! Recalculate deficit from actual evaporation
            if (saev < cofred_squared) then
                ! Stage 1 regime
                spev = saev
            else
                ! Stage 2 regime (inverse function)
                spev = (saev / cofred)**2
            end if
        end if
        
    end subroutine boesten_stroosnijder_reduction

      !> Soil evaporation reduction calculation
      !>
      !! Calculates reduction of potential soil evaporation using either the Black model
      !! or the Boesten-Stroosnijder model.
      !! 
      !! @note
      !! Date: February 2014
      !! Uses module variables from 'variables' module:
      !! - swredu: switch for reduction method (1=Black, 2=Boesten-Stroosnijder)
      !! - fldaystart: flag for start of day
      !! - cofred: coefficient for reduction calculation
      !! - dt: time step
      !! - empreva: empirical potential evaporation (output)
      !! - ldwet: number of days since last significant rainfall
      !! - nird: irrigation amount
      !! - peva: potential evaporation
      !! - pond: ponding depth
      !! - rsigni: threshold rainfall amount
      !! - spev, saev: state variables for Boesten-Stroosnijder method
      !! @endnote
      subroutine reduceva (task, nrai, state)
      ! [SS-SWC S-2.12B] pond retired — read via state%soilwater%pond
      use variables, only: swredu,fldaystart,cofred,dt,    &
     &               nird,rsigni
      implicit none

        ! Arguments
        integer, intent(in) :: task
          !! Task selector: 1 = daily basis, 2 = timestep basis
        real(8), intent(in) :: nrai
          !! Rainfall amount [mm]
        type(swap_state_t), intent(inout) :: state
          !! Simulation state (atmosphere fields dual-written here)

        ! Local variables
        real(8) :: timestep
        real(8), parameter :: POND_THRESHOLD = 1.0d-10  ! Minimum ponding depth [cm]

        ! Validate task
        if (task /= 1 .and. task /= 2) then
            call fatalerr_collected('reduceva', 'Illegal value for TASK')
        end if

        if (task == 1) then
            timestep = 1.0d0  ! Daily
        else
            timestep = dt     ! Sub-daily
        end if

        associate( &
            at_empreva => state%atmosphere%empreva, &
            at_ldwet   => state%atmosphere%ldwet,   &
            at_spev    => state%atmosphere%spev,     &
            at_saev    => state%atmosphere%saev      &
        )

        ! Check for ponding (no reduction needed)
        if (state%soilwater%pond > POND_THRESHOLD) then  ! [SS-SWC S-2.12B]
            at_empreva = state%atmosphere%peva
            at_ldwet   = 0.0d0
            at_spev    = 0.0d0
            at_saev    = 0.0d0
            return
        end if

        ! Apply selected reduction model
        select case (swredu)
        case (1)
            ! Black model
            call black_reduction(nrai, nird, state%atmosphere%peva, cofred, rsigni, &
                                at_ldwet, at_empreva, timestep, fldaystart, task)
        case (2)
            ! Boesten-Stroosnijder model
            call boesten_stroosnijder_reduction(nrai, nird, state%atmosphere%peva, cofred, &
                                              at_spev, at_saev, at_empreva, timestep)
        case default
            call fatalerr_collected('reduceva', 'Unknown reduction method SWREDU')
        end select

        end associate
      end subroutine reduceva
end module et_mod