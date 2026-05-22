!> Shared physical constants for the atmosphere subsystem.
!! Single source of truth for magic numbers that previously appeared
!! inline in snow.f90, runoff.f90, et.f90. Names are uppercase + SI-unit-suffixed.
module atmosphere_constants_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   public

   ! ----- Thermodynamics -----
   !> Specific heat capacity of liquid water [J/kg/K]
   real(real64), parameter :: SPECIFIC_HEAT_WATER = 4180.0_real64

   !> Latent heat of melting for water/ice [J/kg]
   real(real64), parameter :: LATENT_HEAT_MELTING = 333580.0_real64

   !> Reference snow temperature for melt-energy calculation [deg C]
   real(real64), parameter :: SNOW_TEMPERATURE_C = 0.0_real64

   ! ----- Snow -----
   !> Maximum liquid-water storage in snow as a fraction of total water storage
   real(real64), parameter :: SNOW_LIQUID_WATER_FRACTION = 0.07_real64

   !> Soil-surface temperature above which fresh snow cannot accumulate [deg C]
   real(real64), parameter :: SOIL_SURFACE_FREEZE_THRESHOLD_C = 0.5_real64

   ! ----- Ponding -----
   !> Minimum ponding depth at which the soil surface is treated as wet [cm]
   real(real64), parameter :: POND_THRESHOLD_CM = 1.0e-10_real64

   ! ----- Curve Number / runoff -----
   !> Reference top-soil depth for the SCS-CN moisture correction [cm]
   real(real64), parameter :: DEPTH_10CM_CM = 10.0_real64

   !> Pressure head at field capacity used in CN ThetaRef calc [cm]
   real(real64), parameter :: H_FIELD_CAPACITY_CM = -100.0_real64

   !> Pressure head at wilting point used in CN ThetaRef calc [cm]
   real(real64), parameter :: H_WILTING_POINT_CM = -16000.0_real64

   !> Initial-abstraction ratio (Ia/S) in the SCS-CN runoff equation
   real(real64), parameter :: INITIAL_ABSTRACTION_RATIO = 0.2_real64

end module atmosphere_constants_mod
