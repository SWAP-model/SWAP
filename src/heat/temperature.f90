! File VersionID:
!   $Id: temperature.f90 362 2018-01-08 13:08:33Z kroes006 $

!> Module for soil temperature simulation
!!
!! This module handles the calculation of soil temperature profiles in SWAP using
!! either analytical or numerical methods. It includes:
!! - Analytical solution based on sinusoidal temperature waves with depth damping
!! - Numerical solution using the de Vries heat transport model
!! - Calculation of soil thermal properties (heat capacity and conductivity)
!!
!! ## Temperature Calculation Methods
!!
!! **Analytical (swcalt=1):**
!! - Uses sinusoidal function with exponential damping
!! - Suitable for homogeneous soil with constant thermal properties
!! - Fast computation, no numerical stability concerns
!!
!! **Numerical (swcalt=2):**
!! - Solves heat transport equation using finite differences
!! - Uses de Vries model for thermal properties
!! - Accounts for spatial variability in soil properties and water content
!! - Requires boundary conditions (top and/or bottom)
!!
!! ## Boundary Conditions
!!
!! **Top boundary (swtopbhea):**
!! - 1: Air temperature
!! - 2: Prescribed soil surface temperature
!!
!! **Bottom boundary (swbotbhea):**
!! - 1: No heat flow (insulated)
!! - 2: Prescribed temperature
!!
!! @author Original SWAP development team
!! @date Last modified January 2018
module temperature_mod
   use error_mod, only: fatalerr_collected
  implicit none
  private
  public :: temperature_seed, temperature_step, devries

contains

  !> Seed (initialize) soil temperature profile
  !!
  !! Phase-2 initialization: computes initial tsoil profile and (for numerical
  !! method) fquartz/fclay/forg from the soilwater sub-state.  Called once after
  !! SoilWater(1) so that soil%thetas/psand/psilt/pclay/orgmat are available.
  !!
  !! @note
  !! **Original documentation:**
  !! Date: November 2004
  !! Purpose: Calculate soil temperatures (initialization)
  !! @endnote
  subroutine temperature_seed(state, config)
      ! GR-ATM C4: Tav/atav → state%atmosphere; use variables dropped.
      use swap_state_mod,        only: swap_state_t
      use swap_config_mod,       only: swap_config_t
      use array_utils,           only: afgen
      use swap_array_dimensions, only: macp, mabbc
      implicit none

    ! Arguments
    type(swap_state_t),  intent(inout) :: state
    !! Typed simulation state — heat compute writes only state%heat%* (Task 8: dual-write dropped)
    type(swap_config_t), intent(in)    :: config
    !! Typed simulation configuration

    ! Local variables
    integer i,lay, nheat_loc
    real(8) tab(mabbc*2),dummy,gmineral

    associate (heat => state%heat,           &
               mesh => state%mesh,           &
               soil => state%soilwater,     &
               atmo => state%atmosphere,    &
               time => state%timecontrol,   &
               heat_cfg => config%heat)

         ! === Initialization ===

         ! Initial temperature profile.
         if (heat_cfg%swcalt .eq. 1) then
            ! Analytical solution.
            call heat_analytical_profile(state, config)
         else
            ! Numerical solution: use specified initial soil temperatures.
            if (config%soil%swinco .ne. 3 .and. allocated(heat_cfg%tsoil_init)) then
               nheat_loc = size(heat_cfg%tsoil_init, 1)
               do i = 1, nheat_loc
                  tab(i*2)     = heat_cfg%tsoil_init(i, 2)         ! temp — col 2
                  tab(i*2 - 1) = dabs(heat_cfg%tsoil_init(i, 1))   ! depth — col 1
               end do
               do i = 1, mesh%numnod
                  heat%tsoil(i) = afgen(tab, macp*2, dabs(mesh%z(i)))
               end do
            end if
         end if

         if (heat_cfg%swcalt .eq. 2) then
            ! Initialise dry bulk density and volume fractions (sand, clay, OM).
            do i = 1, mesh%numnod
               lay      = mesh%layer(i)
               dummy    = soil%orgmat(lay)/(1.0d0 - soil%orgmat(lay))
               gmineral = (1.0d0 - soil%thetas(i)) / (0.370d0 + 0.714d0*dummy)
               heat%fquartz(i) = (soil%psand(lay) + soil%psilt(lay))*gmineral/2.7d0
               heat%fclay(i)   = soil%pclay(lay)*gmineral/2.7d0
               heat%forg(i)    = dummy*gmineral/1.4d0
            end do
         end if

    end associate

    return
  end subroutine temperature_seed

  !> Advance soil temperature one time step
  !!
  !! Dynamic calculation: updates tsoil at each time step using either the
  !! analytical sinusoidal solution or the de Vries numerical heat transport.
  !!
  !! @note
  !! **Original documentation:**
  !! Date: November 2004
  !! Purpose: Calculate soil temperatures (time step)
  !! @endnote
  subroutine temperature_step(state, config)
      ! GR-ATM C4: Tav/atav → state%atmosphere; use variables dropped.
      use swap_state_mod,        only: swap_state_t
      use swap_config_mod,       only: swap_config_t
      use array_utils,           only: afgen
      use numericalsolvers_mod,  only: tridag
      use swap_array_dimensions, only: macp, mabbc
      implicit none

    ! Arguments
    type(swap_state_t),  intent(inout) :: state
    !! Typed simulation state — heat compute writes only state%heat%* (Task 8: dual-write dropped)
    type(swap_config_t), intent(in)    :: config
    !! Typed simulation configuration

    ! Local variables
    integer i,lay, ierror, nheat_loc
    real(8) tmpold(macp),tab(mabbc*2),ttab(mabbc*2),btab(mabbc*2),dummy,gmineral
    real(8) thoma(macp),thomb(macp),thomc(macp),thomf(macp)
    real(8) theave(macp),heacnd(macp),heacap_loc(macp)
    real(8) heaconbot,qhbot
    real(8) apar, dzsnw, heaconsnw, Rosnw
    character(len=200) messag

    associate (heat => state%heat,           &
               mesh => state%mesh,           &
               soil => state%soilwater,     &
               atmo => state%atmosphere,    &
               time => state%timecontrol,   &
               heat_cfg => config%heat)

         ! === Soil temperature rate and state variables ===

         if (heat_cfg%swcalt .eq. 2) then
            ! Numerical solution.

            ! Top boundary condition.
            if (heat_cfg%swtopbhea .eq. 2) then
               ! Use prescribed soil-surface temperatures.
               ttab = 0.0d0
               if (allocated(heat_cfg%temtoptab)) then
                  do i = 1, min(size(heat_cfg%temtoptab, 1), mabbc)
                     ttab(2*i - 1) = heat_cfg%temtoptab(i, 1)
                     ttab(2*i)     = heat_cfg%temtoptab(i, 2)
                  end do
               end if
               heat%tetop = afgen(ttab, 2*mabbc, time%t1900 + time%dt)
            elseif (dabs(atmo%ssnow) .gt. 1.0d-10) then
               ! Snow layer present — air temperature can't be used directly;
               ! calculate the temperature at the soil/snow interface.
               Rosnw     = 170.0d0
               heaconsnw = 2.86d-6 * 864.0d0 * Rosnw**2.d0
               dzsnw     = atmo%ssnow / 0.170d0
               if (heat%heacon(1) .lt. 1.d-10) heat%heacon(1) = 100.0d0
               apar = (0.5d0*heaconsnw*mesh%dz(1)) / (heat%heacon(1)*dzsnw)
               if (time%flmetdetail) then
                  heat%tetop = (heat%tsoil(1) + apar*atmo%atav(time%wrecord)) / (1.d0 + apar)
               else
                  heat%tetop = (heat%tsoil(1) + apar*atmo%Tav) / (1.d0 + apar)
               end if
            else
               if (time%flmetdetail) then
                  heat%tetop = atmo%atav(time%wrecord)
               else
                  heat%tetop = atmo%Tav
               end if
            end if

            ! Bottom boundary condition.
            if (heat_cfg%swbotbhea .eq. 1) then
               ! Zero heat flow through profile bottom.
               heat%tebot = heat%tsoil(mesh%numnod)
            elseif (heat_cfg%swbotbhea .eq. 2) then
               ! Prescribed bottom temperature.
               btab = 0.0d0
               if (allocated(heat_cfg%tembtab)) then
                  do i = 1, min(size(heat_cfg%tembtab, 1), mabbc)
                     btab(2*i - 1) = heat_cfg%tembtab(i, 1)
                     btab(2*i)     = heat_cfg%tembtab(i, 2)
                  end do
               end if
               heat%tebot = afgen(btab, 2*mabbc, time%t1900 + time%dt)
            end if

            ! Save old temperature profile.
            do i = 1, mesh%numnod
               tmpold(i) = heat%tsoil(i)
            end do

            ! Heat conductivity and capacity.
            do i = 1, mesh%numnod
               theave(i) = 0.5d0 * (soil%theta(i) + soil%thetm1(i))
            end do
            ! heacap_loc is a local workspace (macp-sized) so devries' explicit-shape args fit.
            call devries(mesh%numnod, theave, heacap_loc, heacnd, &
                         heat%fquartz, heat%fclay, heat%forg, soil%thetas)
            heat%heacon(1) = heacnd(1)
            do i = 2, mesh%numnod
               heat%heacon(i) = 0.5d0 * (heacnd(i) + heacnd(i - 1))
            end do
            heat%heacap(1:mesh%numnod) = heacap_loc(1:mesh%numnod)

            ! Build tridiagonal coefficients (node 1: surface-fixed temperature).
            i = 1
            thoma(i) = -time%dt * heat%heacon(i)     / (mesh%dz(i) * mesh%disnod(i))
            thomc(i) = -time%dt * heat%heacon(i + 1) / (mesh%dz(i) * mesh%disnod(i + 1))
            thomb(i) = heat%heacap(i) - thoma(i) - thomc(i)
            thomf(i) = heat%heacap(i) * tmpold(i) - thoma(i) * heat%tetop

            ! 2 < node < numnod.
            do i = 2, mesh%numnod - 1
               thoma(i) = -time%dt * heat%heacon(i)     / (mesh%dz(i) * mesh%disnod(i))
               thomc(i) = -time%dt * heat%heacon(i + 1) / (mesh%dz(i) * mesh%disnod(i + 1))
               thomb(i) = heat%heacap(i) - thoma(i) - thomc(i)
               thomf(i) = heat%heacap(i) * tmpold(i)
            end do

            ! node = numnod.
            i = mesh%numnod
            if (heat_cfg%swbotbhea .eq. 1) then
               ! Zero heat flux through bottom.
               qhbot    = 0.0d0
               thoma(i) = -time%dt * heat%heacon(i) / (mesh%dz(i) * mesh%disnod(i))
               thomb(i) = heat%heacap(i) - thoma(i)
               thomf(i) = heat%heacap(i) * tmpold(i) - (qhbot * time%dt)/mesh%dz(i)
            elseif (heat_cfg%swbotbhea .eq. 2) then
               ! Prescribed bottom temperature.
               heaconBot = heacnd(i)
               thoma(i)  = -time%dt * heat%heacon(i) / (mesh%dz(i) * mesh%disnod(i))
               thomc(i)  = -time%dt * heaconBot      / (mesh%dz(i) * 0.5d0 * mesh%dz(i))
               thomb(i)  = heat%heacap(i) - thoma(i) - thomc(i)
               thomf(i)  = heat%heacap(i) * tmpold(i) - thomc(i) * heat%tebot
            end if

            ! Solve for the new temperature profile.
            call tridag(mesh%numnod, thoma, thomb, thomc, thomf, heat%tsoil, ierror)
            if (ierror .ne. 0) then
               messag = 'During a call from Temperature an error occured in TriDag'
               call fatalerr_collected('Temperature', messag)
            end if
         else
            ! Analytical solution profile.
            call heat_analytical_profile(state, config)
         end if

    end associate

    return
  end subroutine temperature_step

  !> Calculate soil heat capacity and conductivity using de Vries model
  !!
  !! Implements the full de Vries (1963) model for soil thermal properties.
  !! The model accounts for soil composition (mineral, organic, water, air)
  !! and calculates thermal properties as weighted averages of component properties.
  !!
  !! ## Heat Capacity
  !!
  !! Calculated as volume-weighted average of component heat capacities:
  !! \[ C = \sum_i (f_i \rho_i c_i) \]
  !! where \(f_i\) is volume fraction, \(\rho_i\) is density, and \(c_i\) is
  !! specific heat of component i.
  !!
  !! ## Thermal Conductivity
  !!
  !! Calculated using shape factor theory with different formulations for:
  !!
  !! **Dry conditions (θ ≤ 0.02):**
  !! - Air is the main transport medium
  !! - Includes empirical correction factor (1.25)
  !! - Uses air-based weighting factors
  !!
  !! **Wet conditions (θ ≥ 0.05):**
  !! - Liquid water is the main transport medium
  !! - Uses water-based weighting factors
  !!
  !! **Intermediate (0.02 < θ < 0.05):**
  !! - Linear interpolation between dry and wet conductivities
  !!
  !! ## Shape Factors
  !!
  !! Shape factors (G) describe the geometry of soil particles and affect
  !! heat flow paths. They are moisture-dependent for air:
  !! - At θ > 0.02: G = 0.333 - (f_air/θ_sat) × 0.298
  !! - At θ < 0.02: Linear scaling from 0.333 to 0.013
  !!
  !! ## Units
  !!
  !! - Input θ: Volumetric water content (m³/m³)
  !! - Output HeaCap: Heat capacity (J/cm³/K)
  !! - Output HeaCon: Thermal conductivity (J/cm/K/d)
  !!
  !! @note
  !! **Original documentation:**
  !! Purpose: Calculate soil heat capacity and conductivity for each
  !! compartment by full de Vries model
  !!
  !! Description: de Vries model for soil heat capacity and thermal conductivity.
  !! Heat capacity is calculated as average of heat capacities for each soil
  !! component. Thermal conductivity is calculated as weighted average of
  !! conductivities for each component. If theta > 0.05 liquid water is assumed
  !! to be the main transport medium in calculating the weights. If theta < 0.02
  !! air is assumed to be the main transport medium (there is also an empirical
  !! adjustment to the conductivity). For 0.02 < theta < 0.05 conductivity is
  !! interpolated.
  !!
  !! See: Heat and water transfer at the bare soil surface, H.F.M Ten Berge
  !! (pp 48-54 and Appendix 2)
  !!
  !! Input:
  !! - NumNod: number of compartments (-)
  !! - theta: average volumetric soil moisture (m³/m³) — local arg, not VARIABLES
  !! - thetas_in: saturated vol. moisture per node — caller supplies from state%soilwater%thetas [SS-SWC S-2.10]
  !! - fquartz_in, fclay_in, forg_in: volume fractions of sand, clay and org. matter
  !!   (passed explicitly; callers supply from state%heat to avoid stale global reads)
  !!
  !! Output:
  !! - HeaCap: heat capacity (J/m³/K)
  !! - HeaCon: thermal conductivity (W/m/K)
  !! @endnote
  subroutine Devries (numnod_in, theta,HeaCap,HeaCon,fquartz_in,fclay_in,forg_in,thetas_in)
    use swap_array_dimensions, only: macp
    implicit none

    ! Arguments
    integer, intent(in) :: numnod_in
    !! Number of nodes (compartments)
    real(8) theta(macp)
    !! Average water content (m³/m³) - different from theta in VARIABLES
    real(8) HeaCap(MACP)
    !! Output: Heat capacity (J/m³/K)
    real(8) HeaCon(MACP)
    !! Output: Thermal conductivity (W/m/K)
    real(8), intent(in) :: fquartz_in(*)
    !! Volume fraction of quartz per node (caller supplies from state%heat%fquartz)
    real(8), intent(in) :: fclay_in(*)
    !! Volume fraction of clay per node
    real(8), intent(in) :: forg_in(*)
    !! Volume fraction of organic matter per node
    real(8), intent(in) :: thetas_in(*)
    !! Saturated water content per node (caller supplies from state%soilwater%thetas)  [SS-SWC S-2.10]

    ! Local variables
    integer Node
    real(8) kaw
    real(8) fAir(MACP)
    real(8) HeaConDry,HeaConWet
    real(8) GAir,GAirdry

    ! Physical constants - Specific heats (J/kg/K)
    real(8), parameter :: cQuartz =  800.0d0
    real(8), parameter :: cClay   =  900.0d0  
    real(8), parameter :: cWat    = 4180.0d0
    real(8), parameter :: cAir    = 1010.0d0
    real(8), parameter :: cOrg    = 1920.0d0

    ! Physical constants - Density (kg/m³)
    real(8), parameter :: dQuartz = 2660.0d0
    real(8), parameter :: dClay   = 2650.0d0
    real(8), parameter :: dWat    = 1000.0d0
    real(8), parameter :: dAir    =    1.2d0
    real(8), parameter :: dOrg    = 1300.0d0

    ! Physical constants - Thermal conductivities (W/m/K)
    real(8), parameter :: kQuartz = 8.8d0
    real(8), parameter :: kClay   = 2.92d0
    real(8), parameter :: kWat    = 0.57d0
    real(8), parameter :: kAir    = 0.025d0
    real(8), parameter :: kOrg    = 0.25d0

    ! Shape factors for particle geometry
    real(8), parameter :: GQuartz = 0.14d0
    real(8), parameter :: GClay   = 0.125d0
    real(8), parameter :: GWat    = 0.14d0
    real(8), parameter :: GOrg    = 0.5d0

    ! Moisture thresholds
    real(8), parameter :: thetaDry = 0.02d0
    real(8), parameter :: thetaWet = 0.05d0

    ! Weighting factors
    real(8), parameter :: kaa = 1.0d0
    real(8), parameter :: kww = 1.0d0

    ! Weights for each component in conductivity calculations
    real(8), parameter :: kqw = 0.66d0 / (1.0d0 + ((kQuartz/kWat)-1.0d0) * GQuartz) + 0.33d0 / &
                                (1.0d0 + ((kQuartz/kWat) - 1.0d0) * (1.0d0 - 2.0d0 * GQuartz))
    real(8), parameter :: kcw = 0.66d0 / (1.0d0 + ((kClay/kWat) - 1.0d0) * GClay) + 0.33d0 / &
                                (1.0d0 + ((kClay/kWat) - 1.0d0) * (1.0d0 - 2.0d0 * GClay))
    real(8), parameter :: kow = 0.66d0 / (1.0d0 + ((kOrg/kWat) - 1.0d0) * GOrg) + 0.33d0 / &
                                (1.0d0 + ((kOrg/kWat) - 1.0d0) * (1.0d0 - 2.0d0 * GOrg))
    real(8), parameter :: kwa = 0.66d0 / (1.0d0 + ((kWat/kAir) - 1.0d0) * GWat) + 0.33d0 / &
                                (1.0d0 + ((kWat/kAir) - 1.0d0) * (1.0d0 - 2.0d0 * GWat))
    real(8), parameter :: kqa = 0.66d0 / (1.0d0 + ((kQuartz/kAir)-1.0d0) * GQuartz) + 0.33d0 / &
                                (1.0d0 + ((kQuartz/kAir) - 1.0d0) * (1.0d0 - 2.0d0 * GQuartz))
    real(8), parameter :: kca = 0.66d0 / (1.0d0 + ((kClay/kAir) - 1.0d0) * GClay) + 0.33d0 / &
                                (1.0d0 + ((kClay/kAir) - 1.0d0) * (1.0d0 - 2.0d0 * GClay))
    real(8), parameter :: koa = 0.66d0 / (1.0d0 + ((kOrg/kAir) - 1.0d0) * GOrg) + 0.33d0 / &
                                (1.0d0 + ((kOrg/kAir) - 1.0d0) * (1.0d0 - 2.0d0 * GOrg))

    ! Additional derived constants
    real(8), parameter :: kAirDIVkWat = kAir/kWat
    real(8), parameter :: cdQuartz    = dQuartz*cQuartz
    real(8), parameter :: cdClay      = dClay*cClay
    real(8), parameter :: cdWat       = dWat*cWat
    real(8), parameter :: cdAir       = dAir*cAir
    real(8), parameter :: cdOrg       = dOrg*cOrg
    real(8), parameter :: kqaXkQuartz = kqa*kQuartz
    real(8), parameter :: kcaXkClay   = kca*kClay
    real(8), parameter :: kaaXkAir    = kaa*kAir
    real(8), parameter :: koaXkOrg    = koa*kOrg
    real(8), parameter :: kwaXkWat    = kwa*kWat
    real(8), parameter :: kqwXkQuartz = kqw*kQuartz
    real(8), parameter :: kcwXkClay   = kcw*kClay
    real(8), parameter :: kowXkOrg    = kow*kOrg
    real(8), parameter :: kwwXkWat    = kww*kWat

    do Node = 1,numnod_in

      ! (1) Air fraction and related parameters
      fAir(Node) = thetas_in(Node) - theta(Node)              ! [SS-SWC S-2.10]

      ! Determine shape factor of air
      if (theta(node) .gt. thetadry) then
        GAir = 0.333d0 - fair(node)/thetas_in(node)*0.298d0   ! [SS-SWC S-2.10]
      else
        GAirdry = 0.333d0 - fair(node)/thetas_in(node)*0.298d0 ! [SS-SWC S-2.10]
        GAir = 0.013d0 + theta(node)/thetaDry*(GAirdry - 0.013d0)
      endif

      ! Determine weighting factor air - water
      kaw = 0.66d0 / (1.0d0 + ((kAirDIVkWat) - 1.0d0) * GAir) + 0.33d0/ &
            (1.0d0 + ((kAirDIVkWat) - 1.0d0) * (1.0d0 - 2.0d0 * GAir))

      ! (2) Heat capacity (W/m³/K) is average of heat capacities for
      ! all components (multiplied by density for correct units)
      HeaCap(Node) = fquartz_in(Node)*cdQuartz + fclay_in(Node)*cdClay + &
                     theta(Node)*cdWat + fAir(Node)*cdAir + forg_in(Node)*cdOrg

      ! (3) Thermal conductivity (W/m/K) is weighted average of
      ! conductivities of all components

      ! (3.1) Dry conditions (include empirical correction factor 1.25)
      if (theta(Node).LE.thetaDry) Then
        HeaCon(Node) = heacon_dry(theta(Node), fquartz_in(Node), fclay_in(Node), fAir(Node), forg_in(Node))

      ! (3.2) Wet conditions
      else if (theta(Node).GE.thetaWet) Then
        HeaCon(Node) = heacon_wet(theta(Node), fquartz_in(Node), fclay_in(Node), fAir(Node), forg_in(Node), kaw)

      ! (3.3) Intermediate conditions (interpolate between dry and wet)
      else
        ! (3.3.1) Conductivity for theta = 0.02
        HeaConDry = heacon_dry(thetaDry, fquartz_in(Node), fclay_in(Node), fAir(Node), forg_in(Node))

        ! (3.3.2) Conductivity for theta = 0.05
        HeaConWet = heacon_wet(thetaWet, fquartz_in(Node), fclay_in(Node), fAir(Node), forg_in(Node), kaw)

        ! (3.3.3) Interpolate
        HeaCon(Node) = HeaConDry + (theta(Node)-thetaDry) * &
                       (HeaConWet - HeaConDry) / &
                       (thetaWet - thetaDry)
      end if

      ! Conversion of capacity from J/m³/K to J/cm³/K
      HEACAP(NODE) = HEACAP(NODE)*1.0d-6

      ! Conversion of conductivity from W/m/K to J/cm/K/d
      HEACON(NODE) = HEACON(NODE)*864.0d0

    end do

    return

  contains

    !> de Vries weighted conductivity, dry-soil weights (incl. 1.25 correction).
    real(8) function heacon_dry(theta_val, fq, fc, fa, fo)
       real(8), intent(in) :: theta_val, fq, fc, fa, fo
       heacon_dry = 1.25d0 * &
                    (fq*kqaXkQuartz + fc*kcaXkClay + fa*kaaXkAir + fo*koaXkOrg + theta_val*kwaXkWat) / &
                    (kqa*fq + kca*fc + kaa*fa + koa*fo + kwa*theta_val)
    end function heacon_dry

    !> de Vries weighted conductivity, wet-soil weights (kaw is node-specific air-water weight).
    real(8) function heacon_wet(theta_val, fq, fc, fa, fo, kaw_val)
       real(8), intent(in) :: theta_val, fq, fc, fa, fo, kaw_val
       heacon_wet = (fq*kqwXkQuartz + fc*kcwXkClay + fa*kaw_val*kAir + fo*kowXkOrg + theta_val*kwwXkWat) / &
                    (kqw*fq + kcw*fc + kaw_val*fa + kow*fo + kww*theta_val)
    end function heacon_wet

  end subroutine Devries

  !> Analytical soil-temperature profile: damped sinusoidal wave with depth.
  subroutine heat_analytical_profile(state, config)
     use swap_state_mod,  only: swap_state_t
     use swap_config_mod, only: swap_config_t
     implicit none
     type(swap_state_t),  intent(inout) :: state
     type(swap_config_t), intent(in)    :: config
     integer :: i
     associate (heat => state%heat, mesh => state%mesh, &
                time => state%timecontrol, heat_cfg => config%heat)
        do i = 1, mesh%numnod
           heat%tsoil(i) = heat_cfg%tmean + heat_cfg%tampli *                       &
                           (dsin(0.0172d0*(time%daynr - heat_cfg%timref + 91.0d0) + &
                                 mesh%z(i)/heat_cfg%ddamp))                         &
                           / dexp(-mesh%z(i)/heat_cfg%ddamp)
        end do
     end associate
  end subroutine heat_analytical_profile

end module temperature_mod

