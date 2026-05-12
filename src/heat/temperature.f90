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
  public :: temperature, devries, heat_init

contains

  !> Calculate soil temperature profile
  !!
  !! Main driver routine for soil temperature simulation. Operates in two modes:
  !! - Task 1: Initialization (set initial temperature profile)
  !! - Task 2: Dynamic calculation (update temperatures at each time step)
  !!
  !! ## Initialization (task=1)
  !!
  !! Sets initial temperature profile using:
  !! - **Analytical method** (swcalt=1): Sinusoidal wave with depth damping
  !! - **Numerical method** (swcalt=2): Interpolation from specified depths (zh, tsoil)
  !!
  !! For numerical method, also initializes soil composition:
  !! - Dry bulk density
  !! - Volume fractions of quartz, clay, and organic matter
  !!
  !! ## Dynamic Calculation (task=2)
  !!
  !! **Analytical method:**
  !! Updates temperature using sinusoidal function based on day number
  !!
  !! **Numerical method:**
  !! 1. Set boundary conditions (top and bottom temperatures)
  !! 2. Calculate thermal properties using de Vries model
  !! 3. Set up tridiagonal system of equations
  !! 4. Solve for new temperature profile using Thomas algorithm
  !!
  !! The numerical solution accounts for:
  !! - Time-varying soil water content
  !! - Snow layer insulation (if present)
  !! - Spatial variability in thermal properties
  !!
  !! @note
  !! **Original documentation:**
  !! Date: November 2004
  !! Purpose: Calculate soil temperatures
  !! @endnote
  subroutine temperature(task, state)
      use variables
      use array_utils, only: afgen
      use numericalsolvers_mod, only: tridag
      use swap_state_mod, only: swap_state_t
      implicit none

    ! Arguments
    integer,            intent(in)    :: task
    !! Task selector: 1=initialization, 2=calculation
    type(swap_state_t), intent(inout) :: state
    !! Typed simulation state — heat compute writes only state%heat%* (Task 8: dual-write dropped)

    ! Local variables
    integer i,lay, ierror
    real(8) tmpold(macp),tab(mabbc*2),dummy,gmineral
    real(8) thoma(macp),thomb(macp),thomc(macp),thomf(macp)
    real(8) theave(macp),heacnd(macp),heacap_loc(macp)
    real(8) heaconbot,qhbot
    real(8) apar, dzsnw, heaconsnw, Rosnw
    character(len=200) messag

    ! SS-TC TC-11: daynr, t1900, dt read via state%timecontrol tc_* aliases.
    associate( &
        tc_daynr         => state%timecontrol%daynr,    &  ! TC-11
        tc_t1900         => state%timecontrol%t1900,    &  ! TC-11
        tc_dt            => state%timecontrol%dt,       &  ! TC-11
        ht_tsoil         => state%heat%tsoil,           &
        ht_heacap        => state%heat%heacap,          &
        ht_heacon        => state%heat%heacon,          &
        ht_fquartz       => state%heat%fquartz,         &
        ht_fclay         => state%heat%fclay,           &
        ht_forg          => state%heat%forg,            &
        ht_tetop         => state%heat%tetop,           &
        ht_tebot         => state%heat%tebot,           &
        sw_theta         => state%soilwater%theta,      &  ! [SS-SWC S-2.10]
        sw_thetm1        => state%soilwater%thetm1,     &  ! [SS-SWC S-2.10]
        sw_thetas        => state%soilwater%thetas)        ! [SS-SWC S-2.10]

    select case (task)
    case (1)

      ! === Initialization ===

      ! Determine initial temperature profile

      if (swcalt.eq.1) then
        ! Analytical solution
        do i = 1,numnod
          ht_tsoil(i) = tmean+tampli*(dsin(0.0172d0*(tc_daynr-timref+91.0d0)+ &  ! TC-11
                        z(i)/ddamp)) / dexp(-z(i)/ddamp)
        enddo
      else
        ! Numerical solution, use specified soil temperatures
        if (swinco.ne.3) then
          do i = 1, nheat
            tab(i*2)   = tsoil(i)   ! reads config-loaded initial profile (not compute state)
            tab(i*2-1) = dabs(zh(i))
          end do
          do i = 1, numnod
            ht_tsoil(i) = afgen(tab,macp*2,dabs(z(i)))
          end do
        end if
      endif

      if (swcalt.eq.2) then
        ! Initialize dry bulk density and volume fractions sand, clay and organic matter
        do i = 1, numnod
          lay = layer(i)
          dummy = orgmat(lay)/(1.0d0 - orgmat(lay))
          gmineral = (1.0d0 - sw_thetas(i)) / (0.370d0 + 0.714d0*dummy)   ! [SS-SWC S-2.10]
          ht_fquartz(i) = (psand(lay) + psilt(lay))*gmineral/2.7d0
          ht_fclay(i)   = pclay(lay)*gmineral/2.7d0
          ht_forg(i)    = dummy*gmineral/1.4d0
        end do
      endif

      return

    case (2)

      ! === Soil temperature rate and state variables ===

      if (swcalt .eq. 2) then
        ! Numerical solution

        ! Set top boundary condition
        if (swtopbhea .eq. 2) then
          ! Use specified soil surface temperatures as top boundary condition
          ht_tetop = afgen (temtoptab,2*mabbc,tc_t1900+tc_dt)  ! TC-11
        ! SS-ATM A-2.6: ssnow retired — read from state%atmosphere%ssnow
        elseif (dabs(state%atmosphere%ssnow).gt.1.0d-10) then
          ! Air temperature cannot be used with a snow layer,
          ! calculate temperature on soil-snow interface
          Rosnw = 170.0d0
          heaconsnw = 2.86d-6 * 864.0d0 * Rosnw**2.d0
          dzsnw = state%atmosphere%ssnow / 0.170d0
          if (ht_heacon(1).lt.1.d-10) ht_heacon(1) = 100.0d0
          apar = (0.5d0*heaconsnw*dz(1)) / (ht_heacon(1)*dzsnw)
          if (state%timecontrol%flmetdetail) then
            ht_tetop = (ht_tsoil(1) + apar*atav(state%timecontrol%wrecord)) / (1.d0+apar)
          else
            ht_tetop = (ht_tsoil(1) + apar*Tav) / (1.d0+apar)
          endif
        else
          if (state%timecontrol%flmetdetail) then
            ht_tetop = atav(state%timecontrol%wrecord)
          else
            ht_tetop = Tav
          endif
        endif

        ! Set bottom boundary condition
        if (SwBotbHea.eq.1) then
          ! No heat flow through bottom of profile assumed
          ht_tebot = ht_tsoil(Numnod)
        elseif (SwBotbHea.eq.2) then
          ! Bottom temperature is prescribed
          ht_tebot = afgen (tembtab,2*mabbc,tc_t1900+tc_dt)  ! TC-11
        endif

        ! Save old temperature profile
        do i = 1,numnod
          tmpold(i) = ht_tsoil(i)
        enddo

        ! Compute heat conductivity and capacity
        do i = 1,numnod
          theave(i) = 0.5d0 * (sw_theta(i) + sw_thetm1(i))   ! [SS-SWC S-2.10]
        enddo

        ! Calculate nodal heat capacity and thermal conductivity
        ! heacap_loc is a local workspace (macp-sized) so devries explicit-shape args are satisfied
        call devries(theave,heacap_loc,heacnd,ht_fquartz,ht_fclay,ht_forg,sw_thetas)  ! [SS-SWC S-2.10]
        ht_heacon(1) = heacnd(1)
        do i = 2,numnod
          ht_heacon(i) = 0.5d0 * (heacnd(i) + heacnd(i-1))
        enddo
        ht_heacap(1:numnod) = heacap_loc(1:numnod)

        ! Calculate new temperature profile using tridiagonal solver

        ! Calculation of coefficients for node = 1 (temperature fixed at soil surface)
        i = 1
        thoma(i) = - tc_dt * ht_heacon(i) / (dz(i) * disnod(i))    ! TC-11
        thomc(i) = - tc_dt * ht_heacon(i+1) / (dz(i) * disnod(i+1)) ! TC-11
        thomb(i) = ht_heacap(i) - thoma(i) - thomc(i)
        thomf(i) = ht_heacap(i) * tmpold(i) - thoma(i) * ht_tetop

        ! Calculation of coefficients for 2 < node < numnod
        do i = 2,numnod-1
          thoma(i) = - tc_dt * ht_heacon(i) / (dz(i) * disnod(i))    ! TC-11
          thomc(i) = - tc_dt * ht_heacon(i+1) / (dz(i) * disnod(i+1)) ! TC-11
          thomb(i) = ht_heacap(i) - thoma(i) - thomc(i)
          thomf(i) = ht_heacap(i) * tmpold(i)
        enddo

        ! Calculation of coefficients for node = numnod
        i = numnod
        if (SwBotbHea.eq.1) then
          ! No heat flow through bottom of profile assumed
          qhbot = 0.0d0
          thoma(i) = - tc_dt * ht_heacon(i) / (dz(i) * disnod(i))    ! TC-11
          thomb(i) = ht_heacap(i) - thoma(i)
          thomf(i) = ht_heacap(i) * tmpold(i) - (qhbot * tc_dt)/dz(i)  ! TC-11
        elseif (SwBotbHea.eq.2) then
          ! Bottom temperature is prescribed
          heaconBot = heacnd(i)
          thoma(i)  = - tc_dt * ht_heacon(i) / (dz(i) * disnod(i))    ! TC-11
          thomc(i)  = - tc_dt * heaconBot / (dz(i) * 0.5d0 * dz(i))   ! TC-11
          thomb(i)  = ht_heacap(i) - thoma(i) - thomc(i)
          thomf(i)  = ht_heacap(i) * tmpold(i) - thomc(i) * ht_tebot
        endif

        ! Solve for temperature profile; result written directly to state
        call tridag (numnod, thoma, thomb, thomc, thomf, ht_tsoil,ierror)
        if(ierror.ne.0)then
          messag = 'During a call from Temperature an error occured in TriDag'
          call fatalerr_collected ('Temperature',messag)
        end if
      else

        ! Analytical solution temperature profile
        do i = 1,numnod
          ht_tsoil(i) = tmean+tampli*(dsin(0.0172d0*(tc_daynr-timref+91.0d0)+ &  ! TC-11
                        z(i)/ddamp)) / dexp(-z(i)/ddamp)
        enddo

      endif

    case default
      call fatalerr_collected ('Temperature', 'Illegal value for TASK')
    end select

    end associate

    return
  end subroutine Temperature

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
  subroutine Devries (theta,HeaCap,HeaCon,fquartz_in,fclay_in,forg_in,thetas_in)
    use variables, only: NumNod
    use swap_array_dimensions, only: macp
    implicit none

    ! Arguments
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

    do Node = 1,NumNod

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
        HeaCon(Node) = 1.25d0 * &
                       (fquartz_in(Node)*kqaXkQuartz + &
                        fclay_in(Node)*kcaXkClay + &
                        fAir(Node)*kaaXkAir + &
                        forg_in(Node)*koaXkOrg + &
                        theta(Node)*kwaXkWat) / &
                       (kqa * fquartz_in(Node) + kca * fclay_in(Node) + kaa * fAir(Node) + &
                        koa * forg_in(Node) + kwa * theta(Node))

      ! (3.2) Wet conditions
      else if (theta(Node).GE.thetaWet) Then
        HeaCon(Node) = (fquartz_in(Node)*kqwXkQuartz + &
                        fclay_in(Node)*kcwXkClay + &
                        fAir(Node)*kaw*kAir + &
                        forg_in(Node)*kowXkOrg + &
                        theta(Node)*kwwXkWat) / &
                       (kqw * fquartz_in(Node) + kcw * fclay_in(Node) + kaw * fAir(Node) + &
                        kow * forg_in(Node) + kww * theta(Node))

      ! (3.3) Intermediate conditions (interpolate between dry and wet)
      else
        ! (3.3.1) Conductivity for theta = 0.02
        HeaConDry = 1.25d0 * &
                       (fquartz_in(Node)*kqaXkQuartz + &
                        fclay_in(Node)*kcaXkClay + &
                        fAir(Node)*kaaXkAir + &
                        forg_in(Node)*koaXkOrg + &
                        thetaDry*kwaXkWat) / &
                       (kqa * fquartz_in(Node) + kca * fclay_in(Node) + kaa * fAir(Node) + &
                        koa * forg_in(Node) + kwa * thetaDry)

        ! (3.3.2) Conductivity for theta = 0.05
        HeaConWet = (fquartz_in(Node)*kqwXkQuartz + &
                     fclay_in(Node)*kcwXkClay + &
                     fAir(Node)*kaw*kAir + &
                     forg_in(Node)*kowXkOrg + &
                     thetaWet*kwwXkWat) / &
                    (kqw * fquartz_in(Node) + kcw * fclay_in(Node) + kaw * fAir(Node) + &
                     kow * forg_in(Node) + kww * thetaWet)

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
  end subroutine Devries

  !> Lifecycle init for heat typed state. Allocates per-node arrays
  !! from numnod and seeds them. Called from swap_main once per
  !! simulation, after config_to_variables has populated the legacy
  !! globals. Mirrors drainage_init / solute_init patterns.
  subroutine heat_init(state)
    use, intrinsic :: iso_fortran_env, only: real64
    use swap_state_mod, only: swap_state_t
    use Variables, only: numnod
    implicit none
    type(swap_state_t), intent(inout) :: state

    if (.not. allocated(state%heat%tsoil))   allocate(state%heat%tsoil(numnod))
    if (.not. allocated(state%heat%heacap))  allocate(state%heat%heacap(numnod))
    if (.not. allocated(state%heat%heacon))  allocate(state%heat%heacon(numnod))
    if (.not. allocated(state%heat%rfcp))    allocate(state%heat%rfcp(numnod))
    if (.not. allocated(state%heat%fquartz)) allocate(state%heat%fquartz(numnod))
    if (.not. allocated(state%heat%fclay))   allocate(state%heat%fclay(numnod))
    if (.not. allocated(state%heat%forg))    allocate(state%heat%forg(numnod))

    state%heat%tsoil   = 0.0_real64
    state%heat%heacap  = 0.0_real64
    state%heat%heacon  = 0.0_real64
    state%heat%rfcp    = 1.0_real64    ! NB: 1.0 not 0.0 — matches legacy initial value
    state%heat%fquartz = 0.0_real64
    state%heat%fclay   = 0.0_real64
    state%heat%forg    = 0.0_real64

    ! Scalars (tetop, tebot, zfrostbot, zfrosttop, nodfrostbot) keep
    ! their type defaults (zero) — no explicit reset needed here.
  end subroutine heat_init

end module temperature_mod

