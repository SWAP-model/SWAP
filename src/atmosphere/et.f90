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
   use, intrinsic :: iso_fortran_env, only: real64
   use error_mod, only: fatalerr_collected
   use swap_state_mod, only: swap_state_t
   use atmosphere_constants_mod, only: POND_THRESHOLD_CM
   use swap_log, only: log_warn
    implicit none
    private

  public :: PenMon, PenMon_calc, reduceva_daily, reduceva_dt
  public :: pm_inputs_t, pm_outputs_t

  type :: pm_inputs_t
      ! Time / control
      integer :: daynr           = 0
      integer :: irecord         = 1
      integer :: nmetdetail      = 1
      logical :: flmetdetail     = .false.
      logical :: flCropEmergence = .false.
      integer :: swcf            = 1
      integer :: swdivide        = 0

      ! Site
      real(real64) :: lat  = 0.0_real64
      real(real64) :: alt  = 0.0_real64
      real(real64) :: altw = 0.0_real64
      real(real64) :: a    = 0.0_real64
      real(real64) :: b    = 0.0_real64
      real(real64) :: rcs  = 0.0_real64

      ! Atmospheric forcing
      real(real64) :: rad    = 0.0_real64
      real(real64) :: tav    = 0.0_real64
      real(real64) :: tmn    = 0.0_real64
      real(real64) :: tmx    = 0.0_real64
      real(real64) :: hum    = 0.0_real64
      real(real64) :: win    = 0.0_real64
      real(real64) :: atmtr  = 0.0_real64
      real(real64) :: difpp  = 0.0_real64
      real(real64) :: dsinbe = 0.0_real64
      real(real64) :: daylp  = 0.0_real64

      ! Crop / canopy
      real(real64) :: rsc    = 0.0_real64
      real(real64) :: rsw    = 0.0_real64
      real(real64) :: ch     = 0.0_real64
      real(real64) :: albedo = 0.0_real64
      real(real64) :: kdif   = 0.0_real64
      real(real64) :: kdir   = 0.0_real64
      real(real64) :: lai    = 0.0_real64

      ! PMdirect-only
      real(real64) :: rsoil  = 0.0_real64
  end type pm_inputs_t

  type :: pm_outputs_t
      real(real64) :: es0         = 0.0_real64
      real(real64) :: et0         = 0.0_real64
      real(real64) :: ew0         = 0.0_real64
      real(real64) :: Edirect     = 0.0_real64
      real(real64) :: Tdirect     = 0.0_real64
      real(real64) :: Tdirectwet  = 0.0_real64
      real(real64) :: Edirectpond = 0.0_real64
      integer      :: warning_code = 0   ! 0=ok, 1=polar 0hrs, 2=polar 24hrs
  end type pm_outputs_t

contains
  !> Pure Penman-Monteith calculation (no I/O, fully deterministic)
  !!
  !! Performs all ET calculations without side effects. This pure version
  !! enables compiler optimizations, parallelization, and easier testing.
  !!
  !! @note Warnings are returned via outputs%warning_code; see pm_outputs_t.
  pure subroutine PenMon_calc(inputs, outputs)
      use swap_constants, only: vlarge, small, KARMAN_CONSTANT, &
                                GRASS_HEIGHT_CM, MEASUREMENT_HEIGHT_CM, &
                                BARE_SOIL_HEIGHT_CM, PI, ALBEDO_PONDING
      implicit none

      type(pm_inputs_t),  intent(in)  :: inputs
      type(pm_outputs_t), intent(out) :: outputs

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

      ! Outputs are auto-initialised by intent(out) on the derived type
      ! (defaults in pm_outputs_t).

      ! ========================================================================
      ! 1. PREPROCESSING: Unit conversions and crop height determination
      ! ========================================================================

      zmeasw = 100.0d0 * inputs%altw  ! Convert wind measurement height to cm

      ! Determine effective crop height
      if (.not. inputs%flCropEmergence) then
          chplant = GRASS_HEIGHT_CM
      else
          if (inputs%swcf == 1 .or. inputs%swcf == 3) then
              chplant = GRASS_HEIGHT_CM
          else
              chplant = max(inputs%ch, 0.1d0)
          endif
      endif

      ! ========================================================================
      ! 2. ATMOSPHERIC PROPERTIES
      ! ========================================================================

      ! Temperature conversions [K]
      tavk = inputs%tav + 273.15d0
      tmnk = inputs%tmn + 273.15d0
      tmxk = inputs%tmx + 273.15d0

      ! Atmospheric pressure at elevation [kPa]
      palt = 101.3d0 * ((tavk - 0.0065d0*inputs%alt) / tavk)**5.26d0

      ! Latent heat of vaporization [MJ/kg]
      lambda = 2.501d0 - 0.002361d0*inputs%tav

      ! Saturation vapour pressure [kPa]
      if (inputs%flmetdetail) then
          ea = 0.611d0 * exp(17.27d0*inputs%tav / (inputs%tav + 237.3d0))
      else
          ea = 0.3055d0 * (exp(17.27d0*inputs%tmn / (inputs%tmn + 237.3d0)) + &
                          exp(17.27d0*inputs%tmx / (inputs%tmx + 237.3d0)))
      endif

      ! Measured vapour pressure (capped at saturation)
      ed = min(inputs%hum, ea)

      ! Vapour pressure deficit [kPa]
      vpd = ea - ed

      ! Slope of vapour pressure curve [kPa/C]
      delta = 4098.0d0 * ea / (inputs%tav + 237.3d0)**2

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
      ud = max(inputs%win, 0.0001d0)

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
      if (inputs%swdivide == 1) then
          rss = inputs%rsoil  ! PMdirect partitioning
      else
          rss = 0.0d0
      endif

      ! Modified psychrometric constants [kPa/C]
      gammos = gamma * (1.0d0 + rss/ras)
      gammoc = gamma * (1.0d0 + inputs%rsc/rac)
      gammow = gamma * (1.0d0 + inputs%rsw/raw)

      ! ========================================================================
      ! 4. RADIATION CALCULATIONS
      ! ========================================================================

      ! Net shortwave radiation [MJ/m2/d]
      rns = (1.0d0 - inputs%rcs) * inputs%rad / 1.0d6
      rnc = (1.0d0 - inputs%albedo) * inputs%rad / 1.0d6
      rnw = (1.0d0 - inputs%albedo) * inputs%rad / 1.0d6
      rnp = (1.0d0 - ALBEDO_PONDING) * inputs%rad / 1.0d6

      ! Extraterrestrial radiation and daylength
      if (inputs%flmetdetail) then
          ! Sub-daily: calculate astronomical parameters
          radial = PI / 180.0d0
          dec = -asin(sin(23.45d0*radial) * &
                      cos(2.0d0*PI*real(inputs%daynr+10, real64) / 365.0d0))

          sinld = sin(radial*inputs%lat) * sin(dec)
          cosld = cos(radial*inputs%lat) * cos(dec)
          aob = sinld / cosld

          ! Daylength calculation with polar circle handling
          if (aob < -1.0d0) then
              dayl = 0.0d0
              outputs%warning_code = 1   ! Polar circle, 0 hours
          else if (aob > 1.0d0) then
              dayl = 24.0d0
              outputs%warning_code = 2   ! Polar circle, 24 hours
          else
              dayl = 12.0d0 * (1.0d0 + 2.0d0*asin(aob)/PI)
          endif

          sunrise = 0.5d0 - dayl / 48.0d0
          sunset = 0.5d0 + dayl / 48.0d0
          startrec = real(real(inputs%irecord-1) / real(inputs%nmetdetail), real64)
          endrec = real(real(inputs%irecord) / real(inputs%nmetdetail), real64)
      else
          ! Daily: use provided parameters
          ! Note: sinld, cosld would come from astro() call in wrapper
          sinld = inputs%dsinbe
          cosld = sqrt(max(0.0d0, 1.0d0 - sinld**2))
          startrec = 0.0d0
          endrec = 1.0d0
          sunrise = 0.0d0
          sunset = 1.0d0
      endif

      ! Net longwave radiation [MJ/m2/d]
      relssd = max(min((inputs%atmtr - inputs%a) / inputs%b, 1.0d0), 0.0d0)
      rnl = 4.9d-9 * 0.5d0 * (tmxk**4 + tmnk**4) * &
            (0.34d0 - 0.14d0*sqrt(ed)) * (0.1d0 + 0.9d0*relssd)

      ! Soil heat flux [MJ/m2/d]
      if (inputs%flmetdetail) then
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
      outputs%es0 = max(0.0d0, etaers + etrads)
      outputs%et0 = max(0.0d0, etaerc + etradc)
      outputs%ew0 = max(0.0d0, etaerw + etradw)

      ! ========================================================================
      ! 6. PENMAN-MONTEITH DIRECT PARTITIONING (PMdirect)
      ! ========================================================================

      if (inputs%swdivide == 1) then
          ! Vegetation cover fraction
          vcover = 1.0d0 - exp(-inputs%kdif * inputs%kdir * inputs%lai)

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
          laieff = inputs%lai / (0.3d0*inputs%lai + 1.2d0)

          ! Modified psychrometric constants with resistance scaling
          gammos = vlarge
          if (ras > small) gammos = gamma * (1.0d0 + rss/ras)

          gammoc = vlarge
          if ((rac*laieff) > small) gammoc = gamma * (1.0d0 + inputs%rsc/(rac*laieff))

          gammow = vlarge
          if ((raw*laieff) > small) gammow = gamma * (1.0d0 + inputs%rsw/(raw*laieff))

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
          outputs%Edirect     = max(0.0d0, etaers + etrads)
          outputs%Tdirect     = max(0.0d0, etaerc + etradc)
          outputs%Tdirectwet  = max(0.0d0, etaerw + etradw)
          outputs%Edirectpond = max(0.0d0, etaerp + etradp)
      endif

  end subroutine PenMon_calc

  !> Penman-Monteith evapotranspiration calculation (wrapper with I/O)
  !!
  !! Thin wrapper around PenMon_calc that calls astro() for the daily
  !! branch and forwards outputs%warning_code to log_warn().
  subroutine PenMon(inputs, outputs)
      implicit none

      type(pm_inputs_t),  intent(in)    :: inputs
      type(pm_outputs_t), intent(inout) :: outputs

      ! Local variables
      real(8) :: dayl, sinld, cosld
      character(len=200) :: messag
      type(pm_inputs_t)  :: inputs_local

      ! Make a mutable working copy so astro() can backfill the
      ! daily-branch astronomical fields (dayl/sinld/cosld are local;
      ! daylp/difpp/atmtr/dsinbe land in inputs_local).
      inputs_local = inputs

      ! Call astro() for daily radiation if needed
      if (.not. inputs_local%flmetdetail) then
          call astro(inputs_local%daynr, inputs_local%lat, inputs_local%rad, &
                     dayl, inputs_local%daylp, sinld, cosld, &
                     inputs_local%difpp, inputs_local%atmtr, inputs_local%dsinbe)
      endif

      ! Call pure calculation core
      call PenMon_calc(inputs_local, outputs)

      ! Handle warnings
      if (outputs%warning_code == 1) then
          messag = 'Warning: latitude above polar circle, daylength = 0hrs'
          call log_warn('Astro', messag)
      else if (outputs%warning_code == 2) then
          messag = 'Warning: latitude within polar circle, daylength = 24hrs'
          call log_warn('Astro', messag)
      endif

  end subroutine PenMon

    !> Black's evaporation reduction model
    !!
    !! Reduces soil evaporation based on square root of time since last rainfall.
    !! Formulation: E = β(√t - √(t-Δt)) / Δt
    !!
    pure subroutine black_reduction(nrai, nird, peva, cofred, rsigni, &
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
    pure subroutine boesten_stroosnijder_reduction(nrai, nird, peva, cofred, &
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
      !> Daily-basis soil evaporation reduction (formerly reduceva(task=1, ...))
      subroutine reduceva_daily(nrai, state)
        implicit none
        real(8),            intent(in)    :: nrai
        type(swap_state_t), intent(inout) :: state
        call reduceva_apply(1, nrai, state, daily=.true.)
      end subroutine reduceva_daily

      !> Sub-daily soil evaporation reduction (formerly reduceva(task=2, ...))
      subroutine reduceva_dt(nrai, state)
        implicit none
        real(8),            intent(in)    :: nrai
        type(swap_state_t), intent(inout) :: state
        call reduceva_apply(2, nrai, state, daily=.false.)
      end subroutine reduceva_dt

      !> Shared body for reduceva_daily / reduceva_dt (private helper).
      !!
      !! task_flag is forwarded to black_reduction (which still uses it to
      !! gate daily-vs-dt formulas internally). The daily logical flag
      !! selects the timestep value (1.0 d vs tc_dt).
      subroutine reduceva_apply(task_flag, nrai, state, daily)
      ! [SS-SWC S-2.12B] pond retired — read via state%soilwater%pond
      ! SS-TC TC-11: dt, fldaystart read via state%timecontrol tc_* aliases.
      ! [SS-GR-ATM B2] DEFERRED — swredu/cofred/rsigni config paths verified
      !   (config%meteo%evaporation%swredu, cofredbl/bo, rsigni) but caller chain
      !   (ProcessMeteoDay, MeteoDT, ProcessMeteoTsteps, ReadMeteoYear) lacks config arg;
      !   threading deferred to Tasks 22-29 (meteoday/meteodt migration).
      use variables, only: &   ! [SS-GR-FINAL B11] residuals — all DEFERRED
         ! DEFERRED: swredu/cofred — ET reduction switch/coefficient; config%meteo%evaporation; Phase C3
         swredu, cofred, &
         ! DEFERRED: rsigni — significant radiation threshold; config; Phase C3
         rsigni
      implicit none

        ! Arguments
        integer, intent(in) :: task_flag
          !! Forwarded to black_reduction: 1 = daily, 2 = sub-daily
        real(8), intent(in) :: nrai
          !! Rainfall amount [mm]
        type(swap_state_t), intent(inout) :: state
          !! Simulation state (atmosphere fields dual-written here)
        logical, intent(in) :: daily
          !! .true. selects daily timestep (1 d); .false. uses tc_dt

        ! Local variables
        real(8) :: timestep

        ! SS-TC TC-11: dt, flDayStart read via state%timecontrol tc_* aliases.
        associate( &
            tc_dt        => state%timecontrol%dt,         &  ! TC-11
            tc_flDayStart => state%timecontrol%flDayStart, &  ! TC-11
            at_empreva => state%atmosphere%empreva, &
            at_ldwet   => state%atmosphere%ldwet,   &
            at_spev    => state%atmosphere%spev,     &
            at_saev    => state%atmosphere%saev      &
        )

        if (daily) then
            timestep = 1.0d0  ! Daily
        else
            timestep = tc_dt  ! Sub-daily  ! TC-11
        end if

        ! Check for ponding (no reduction needed)
        if (state%soilwater%pond > POND_THRESHOLD_CM) then  ! [SS-SWC S-2.12B]
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
            call black_reduction(nrai, state%atmosphere%nird, state%atmosphere%peva, cofred, rsigni, &
                                at_ldwet, at_empreva, timestep, tc_flDayStart, task_flag)  ! TC-11
        case (2)
            ! Boesten-Stroosnijder model
            call boesten_stroosnijder_reduction(nrai, state%atmosphere%nird, state%atmosphere%peva, cofred, &
                                              at_spev, at_saev, at_empreva, timestep)
        case default
            call fatalerr_collected('reduceva_apply', 'Unknown reduction method SWREDU')
        end select

        end associate
      end subroutine reduceva_apply
end module et_mod