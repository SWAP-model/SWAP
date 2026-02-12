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
    implicit none
    private

  public :: PenMon, reduceva, reduceva_state
contains
      !> Penman-Monteith evapotranspiration calculation
      !>
      !! Calculates potential evaporation and transpiration rates from a bare soil surface,
      !! a dry crop canopy, and a wet crop canopy based on the Penman-Monteith approach.
      !>
      !! @param[in]     logf              Internal number of logbook output file
      !! @param[in]     swscre            Switch of screen display: 0 = no display; 1 = summary water balance; 2 = daynumber
      !! @param[in]     daynr             Day number (January 1st = 1) [-]
      !! @param[in]     lat               Latitude [deg, decimal degrees, N=+, S=-]
      !! @param[in]     alt               Altitude above mean sea level [m]
      !! @param[in]     altw              Altitude of wind speed measurement [m]
      !! @param[in]     a                 First Angstrom coefficient [-]
      !! @param[in]     b                 Second Angstrom coefficient [-]
      !! @param[in]     rcs               Reflection coefficient soil [-]
      !! @param[in]     rad               Incoming short wave radiation [J/m2/d]
      !! @param[in]     tav               Average temperature (24 hour) [C]
      !! @param[in]     hum               Vapour pressure [kPa]
      !! @param[in]     win               Wind speed at 2 m height [m/s]
      !! @param[in]     rsc               Minimum canopy resistance of dry crop [s/m]
      !! @param[out]    es0               Potential evaporation rate from a wet bare soil [mm/d]
      !! @param[out]    et0               Potential transpiration rate from a dry crop [mm/d]
      !! @param[out]    ew0               Potential transpiration rate from a wet crop [mm/d]
      !! @param[in]     swcf              Switch use crop factor (=1) or crop height (=2)
      !! @param[in]     ch                Crop height [cm]
      !! @param[in]     flCropEmergence   Logical flag for crop emergence
      !! @param[in]     daylp             Day length parameter
      !! @param[in]     flmetdetail       Logical flag for detailed meteorological data
      !! @param[in]     irecord           Current record number
      !! @param[in]     nmetdetail        Number of detailed meteorological records
      !! @param[in]     albedo            Reflection coefficient crop [-]
      !! @param[in]     tmn               Minimum temperature [C]
      !! @param[in]     tmx               Maximum temperature [C]
      !! @param[in]     rsw               Canopy resistance of intercepted water [s/m]
      !! @param[in]     difpp             Diffuse radiation parameter
      !! @param[in]     dsinbe            Solar elevation parameter
      !! @param[in]     atmtr             Daily atmospheric transmission
      !! @param[out]    Edirect           Direct evaporation [mm/d]
      !! @param[out]    Tdirect           Direct transpiration [mm/d]
      !! @param[out]    Tdirectwet        Direct transpiration from wet crop [mm/d]
      !! @param[in]     rsoil             Soil resistance in PMdirect [s/m]
      !! @param[in]     swdivide          Switch for direct partitioning method
      !! @param[in]     kdif              Diffuse light extinction coefficient [-]
      !! @param[in]     kdir              Direct light extinction coefficient [-]
      !! @param[in]     lai               Leaf area index [-]
      !! @param[out]    Edirectpond       Direct evaporation from ponding layer [mm/d]
      !>
      !! @note
      !! Date: 14/01/99
      !! Local variables:
      !! - ckarman: von Karman constant [-]
      !! - zm: height of wind speed [cm]
      !! - zh: height of temperature and humidity measurement [cm]
      !! - d: zero displacement of wind profile [cm]
      !! - zom: roughness parameter for momentum [cm]
      !! - zoh: roughness parameter for heat and vapour [cm]
      !! - zmeasw: altitude of wind speed measurement [m]
      !! - Vcover: vegetation cover [-]
      !! @endnote
      subroutine PenMon (logf,swscre,daynr,lat,alt,altw,a,b,rcs,rad,    &
     & tav,hum,win,rsc,es0,et0,ew0,swcf,ch,flCropEmergence,             &
     & daylp,flmetdetail,irecord,nmetdetail,albedo,tmn,tmx,rsw,difpp,   &
     & dsinbe,atmtr,Edirect,Tdirect,Tdirectwet,rsoil,swdivide,kdif,kdir,&
     & lai,Edirectpond)
use swap_constants, only: vlarge, small
      implicit none

      ! Global variables
      integer daynr,swcf,logf,swscre,irecord,nmetdetail,swdivide
      real(8) lat,alt,altw,albedo,tmn,tmx,ch,difpp,dsinbe
      real(8) a,b,rcs,rsc,rsw,es0,et0,ew0,hum,rad,tav,win,atmtr
      real(8) Edirect,Tdirect,Tdirectwet,rsoil,kdif,kdir,lai
      real(8) Edirectpond
      logical flCropEmergence,flmetdetail
      
      ! Local variables
      real(8) lambda,cosld,dayl,delta,ea,vcover
      real(8) ed,etaerc,etaers,etaerw,etradc,etrads,etradw
      real(8) gamma,gammoc,gammos,gammow,palt,rac,ras,raw,cp
      real(8) relssd,rho,rnc,rnl,rns,rnw,rss,sinld,tavk,tkv
      real(8) ud,vpd,daylp,gs,gc,gw,radial,dec,aob,tmnk,tmxk
      real(8) ckarman, chgrass, chsoil, d, dgrass, chplant
      real(8) zm, zmeasw, zmeash, zom, zomgrass, zh,zoh
      real(8) zact, dact, zomact, fact, fmeas
      real(8) sunrise,sunset,startrec,endrec,pi,laieff
      real(8) albpond,rnp,gammop,etaerp,etradp
      character(len=200) messag
      
      ! Local parameters
      data    chgrass  /12.0d0/      ! Height of reference crop grassland [cm]
      data    zmeash   /200.0d0/     ! Default height of humidity and temperature measurement [cm]
      data    ckarman  /0.41d0/      ! von Karman constant [-]
      data    chsoil   /0.1d0/       ! Nihil crop height for a wet bare soil [cm]
      data    pi       /3.141593d0/  ! Number pi
      data    albpond  /0.08d0/      ! Albedo of ponding layer


      ! Conversion to cm
      zmeasw = 100.0d0 * altw

      ! Avoid zero crop height
      if (.not. flCropEmergence) then
        chplant = chgrass
      else 
        if (swcf.eq.1 .or. swcf.eq.3) then
          chplant = chgrass
        else
          chplant = max (ch,0.1d0)
        endif
      endif

      ! Conversion of temperature from [C] to [K]
      tavk = tav+273.15d0
      tmnk = tmn+273.15d0
      tmxk = tmx+273.15d0

      ! Atmospheric pressure at elevation alt [kPa]
      palt = 101.3d0*((tavk-0.0065d0*alt)/tavk)**5.26d0

      ! Latent heat of vaporization [MJ/kg]
      lambda = 2.501d0-0.002361d0*tav

      ! Saturation vapour pressure [kPa]
      if (flmetdetail) then
        ea = 0.611d0*exp(17.27d0*tav/(tav+237.3d0))
      else
        ea = 0.3055d0*(exp(17.27d0*tmn/(tmn+237.3d0)) +                 &
     &                 exp(17.27d0*tmx/(tmx+237.3d0)))
      endif

      ! Measured vapour pressure not to exceed saturated vapour pressure
      ed = min(hum,ea)
      
      ! Vapour pressure deficit [kPa]
      vpd = ea-ed

      ! Slope vapour pressure curve [kPa/C]
      delta = 4098.0d0*ea/(tav+237.3d0)**2

      ! Psychrometric constant [kPa/C]
      gamma = 0.00163d0*palt/lambda

      ! Atmospheric density [kg/m3]
      tkv = tavk/(1.0d0-0.378d0*ed/palt)
      rho = 3.486d0*palt/tkv

      ! Specific heat moist air [kJ/kg/C]
      cp = 622.0d0*gamma*lambda/palt

      ! Aerodynamic resistance [s/m] - soil, crop & wet crop
      
      ! Day wind [m/s] for daily records, avoid zero windspeed
      ud = max (win, 0.0001d0)

      ! Adjust wind speed if crop height deviates from measurement height of wind speed
      ! assuming equal wind speed at 100 meter (1.0d4 cm) above the soil surface
      ! (all length-units in cm)
      if (chplant.gt.zmeash .or. zmeasw.gt.zmeash) then
        dgrass = 2.0d0/3.0d0 * chgrass
        zomgrass = 0.123d0 * chgrass
        fmeas = log((1.0d4-dgrass)/zomgrass) /                          &
     &         log((zmeasw-dgrass)/zomgrass) 
        zact = max(chplant,200.0d0)
        dact = 2.0d0/3.0d0 * chplant
        zomact = 0.123d0 * chplant
       fact = log((zact-dact)/zomact) / log((1.0d4-dact)/zomact)
       ud = ud * fact * fmeas
      endif

      ! Constants to determine aerodynamic resistance
      zm = max(chplant,zmeash)     ! Measured height of wind speed measurement [cm]
      zh = zm                      ! Height of humidity and temperature measurement [cm]
      d = 2.0d0/3.0d0 * chplant    ! Zero displacement of wind profile [cm]
      zom = 0.123d0 * chplant      ! Roughness parameter for momentum [cm]
      zoh = 0.1d0 * zom            ! Roughness parameter for heat and vapour [cm]
      
      ! Aerodynamic resistance for dry and wet crop
      rac = log ((zm - d)/zom) * log((zh - d)/zoh) / ckarman**2/ud 
      raw = rac

      ! Aerodynamic resistance for bare wet soil with a crop height set to 0.1 cm
      d = 2.0d0/3.0d0 * chsoil    ! Zero displacement of wind profile at low height [cm]
      zom = 0.123d0 * chsoil      ! Roughness parameter for momentum at low height [cm]
      zoh = 0.1d0 * zom           ! Roughness parameter for heat and vapour [cm]
      ras = log ((zm - d)/zom) * log((zh - d)/zoh) / ckarman**2/ud

      ! Surface resistance of wet soil [s/m]
      if (swdivide .eq. 1) then
        ! Apply specified soil resistance for PMdirect partitioning
        rss = rsoil
      else
        rss = 0.d0
      endif

      ! Modified psychrometric constant [kPa/C] - soil, crop & wet crop 
      gammos = gamma*(1.0d0+rss/ras)
      gammoc = gamma*(1.0d0+rsc/rac)
      gammow = gamma*(1.0d0+rsw/raw)

      ! Net short wave radiation [MJ/m**2/d] - soil, crop, wet crop & pond layer
      rns = (1.0d0-rcs)*rad/1000000.0d0
      rnc = (1.0d0-albedo)*rad/1000000.0d0
      rnw = (1.0d0-albedo)*rad/1000000.0d0
      rnp = (1.0d0-albpond)*rad/1000000.0d0

      ! Procedure to derive extraterrestrial radiation [MJ/m2/day]
      if (flmetdetail) then

        ! Declination of the sun as a function of daynr
        radial = pi/180.d0
        dec = -asin(dsin(23.45d0*radial)*                                &
     &                        dcos(2.d0*pi*dble(daynr+10)/365.0d0))
        
        ! Some intermediate variables
        sinld = dsin(radial*lat)*dsin(dec)
        cosld = dcos(radial*lat)*dcos(dec)
        aob = sinld/cosld
        
        ! Calculation of astronomical daylength
        if (aob.lt.-1.0d0) then
          messag='Warning: latitude above polar circle, daylength= 0hrs'
          call warn ('Astro',messag,logf,swscre)
        else if (aob.gt.1.0d0) then
          messag='Warning: latitude within polar circle,daylength=24hrs'
          call warn ('Astro',messag,logf,swscre)
        else
          dayl  = 12.0d0*(1.0d0+2.0d0*asin(aob)/pi)
        endif

        ! Extraterrestrial radiation of current period (dso) [J/m2/day]
        sunrise = 0.5d0 - dayl / 48.d0
        sunset = 0.5d0 + dayl / 48.d0
        startrec = dble(real(irecord-1)/real(nmetdetail))
        endrec = dble(real(irecord)/real(nmetdetail))

      else
        ! Just daily extraterrestrial radiation (one record per day)
        call astro(daynr,lat,rad,dayl,daylp,sinld,cosld,difpp,          &
     &             atmtr,dsinbe)
      endif

      ! Net long wave radiation (MJ/m2/d)
      relssd = max(min((atmtr-a)/b,1.0d0),0.0d0)
      rnl = 4.9d-9* 0.5d0 * (tmxk**4 + tmnk**4) *                       &
     &      (0.34d0-0.14d0*dsqrt(ed))*(0.1d0+0.9d0*relssd)
 
      ! Soil heat flux [MJ/m2/d]
      if (flmetdetail) then
        if ((startrec+endrec)/2.d0 .gt. sunrise .and.                   &
     &      (startrec+endrec)/2.d0 .lt. sunset) then
          ! Daytime period
          gs = 0.1d0 * (rns - rnl)
          gc = 0.1d0 * (rnc - rnl)
          gw = 0.1d0 * (rnw - rnl)
        else
          ! Nighttime period
          gs = 0.5d0 * (rns - rnl)
          gc = 0.5d0 * (rnc - rnl)
          gw = 0.5d0 * (rnw - rnl)
        endif
      else
        ! Daily record, net flux negligable
        gs = 0.d0
        gc = 0.d0
        gw = 0.d0
      endif

      ! Aerodynamic term of the PM equation [mm/d] - soil, crop & wet crop
      etaers = (86.4d0/lambda)*(1.0d0/(delta+gammos))*(rho*cp*vpd/ras)
      etaerc = (86.4d0/lambda)*(1.0d0/(delta+gammoc))*(rho*cp*vpd/rac)
      etaerw = (86.4d0/lambda)*(1.0d0/(delta+gammow))*(rho*cp*vpd/raw)

      ! Radiation term of the PM equation [mm/d] - soil, crop & wet crop
      etrads = delta/(delta+gammos)*(rns-rnl-gs)*1.0d0/lambda      
      etradc = delta/(delta+gammoc)*(rnc-rnl-gc)*1.0d0/lambda      
      etradw = delta/(delta+gammow)*(rnw-rnl-gw)*1.0d0/lambda      

      ! Sum of both terms [mm/d] - soil, crop & wet crop      
      es0 = max (0.0d0,etaers+etrads)
      et0 = max (0.0d0,etaerc+etradc)
      ew0 = max (0.0d0,etaerw+etradw)

      ! PMdirect: potential transpiration Tdirect and potential evaporation Edirect
      if (swdivide .eq. 1) then

        ! Determine vegetation cover Vcover
        Vcover = 1.0d0 - exp(-1.0d0*kdif*kdir*lai)

        ! Adjust aerodynamic resistances
        if (Vcover .gt. 1.d-6) then
          rac = rac / Vcover
        else
          rac = 1.d12
        endif
        raw = rac

        if ((1.d0 - Vcover) .gt. 1.d-6) then
          ras = ras / (1.d0 - Vcover)
        else
          ras = 1.d12
        endif

        ! Effective LAI
        LAIeff = lai / (0.3d0*lai + 1.2d0)

        ! Modified psychrometric constant [kPa/C] - crop, wet crop & pond layer 
        gammos = vlarge
        if (ras .gt. small) gammos = gamma*(1.0d0+rss/ras)
        gammoc = vlarge
        if ((rac*LAIeff) .gt. small) gammoc = gamma*(1.0d0+rsc/(rac*LAIeff))
        gammow = vlarge
        if ((raw*LAIeff) .gt. small) gammow = gamma*(1.0d0+rsw/(raw*LAIeff))
        gammop = vlarge
        if (ras .gt. small) gammop = gamma

        ! Aerodynamic term of the PM equation [mm/d] - crop, wet crop & pond layer
        etaers = (86.4d0/lambda)*(1.0d0/(delta+gammos))*(rho*cp*vpd/ras)
        etaerc = (86.4d0/lambda)*(1.0d0/(delta+gammoc))*(rho*cp*vpd/rac)
        etaerw = (86.4d0/lambda)*(1.0d0/(delta+gammow))*(rho*cp*vpd/raw)
        etaerp = (86.4d0/lambda)*(1.0d0/(delta+gammop))*(rho*cp*vpd/ras)

        ! Radiation term of the PM equation [mm/d] - crop, wet crop & pond layer
        etrads = delta/(delta+gammos)*(rns-rnl-gs)*(1.d0-vcover)        &
     &           *1.0d0/lambda
        etradc = delta/(delta+gammoc)*(rnc-rnl-gc)*vcover*1.0d0/lambda      
        etradw = delta/(delta+gammow)*(rnw-rnl-gw)*vcover*1.0d0/lambda      
        etradp = delta/(delta+gammop)*(rnp-rnl-gs)*(1.d0-vcover)        &
     &           *1.0d0/lambda

        ! Sum of both terms [mm/d] - crop, wet crop & pond layer
        Edirect = max (0.0d0,etaers+etrads)
        Tdirect = max (0.0d0,etaerc+etradc)
        Tdirectwet = max (0.0d0,etaerw+etradw)
        Edirectpond = max (0.0d0,etaerp+etradp)

      endif

      return
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
      subroutine reduceva (task,nrai)
      use variables, only: swredu,fldaystart,cofred,dt,empreva,    &
     &               ldwet,nird,peva,pond,rsigni,spev,saev
      implicit none

        ! Arguments
        integer, intent(in) :: task
          !! Task selector: 1 = daily basis, 2 = timestep basis
        real(8), intent(in) :: nrai
          !! Rainfall amount [mm]

        ! Local variables
        real(8) :: timestep
        real(8), parameter :: POND_THRESHOLD = 1.0d-10  ! Minimum ponding depth [cm]

        ! Validate task
        if (task /= 1 .and. task /= 2) then
            call fatalerr('reduceva', 'Illegal value for TASK')
        end if

        if (task == 1) then
            timestep = 1.0d0  ! Daily
        else
            timestep = dt     ! Sub-daily
        end if

        ! Check for ponding (no reduction needed)
        if (pond > POND_THRESHOLD) then
            empreva = peva
            ldwet = 0.0d0
            spev = 0.0d0
            saev = 0.0d0
            return
        end if

        ! Apply selected reduction model
        select case (swredu)
        case (1)
            ! Black model
            call black_reduction(nrai, nird, peva, cofred, rsigni, &
                                ldwet, empreva, timestep, fldaystart, task)
        case (2)
            ! Boesten-Stroosnijder model
            call boesten_stroosnijder_reduction(nrai, nird, peva, cofred, &
                                              spev, saev, empreva, timestep)
        case default
            call fatalerr('reduceva', 'Unknown reduction method SWREDU')
        end select
      end subroutine reduceva
      !> State-aware wrapper for `reduceva`
      !!
      !! Bridges legacy module-variable implementation with explicit model state.
      !! The wrapper restores legacy variables from `state`, executes the
      !! original routine, and snapshots the updated values back into `state`.
      !!
      !! @param[inout] state SWAP model state container
      !! @param[in]    task  Task selector: 1=daily, 2=timestep
      !! @param[in]    nrai  Rainfall amount [mm]
      subroutine reduceva_state(state, task, nrai)
      use swap_state_mod, only: swap_state_t
      use swap_state_sync, only: atmosphere_state_from_variables
      implicit none

      type(swap_state_t), intent(inout) :: state
      integer,            intent(in)     :: task
      real(8),            intent(in)     :: nrai

      call reduceva(task, nrai)
      call atmosphere_state_from_variables(state%atm)

      end subroutine reduceva_state
end module et_mod