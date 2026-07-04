!> @file exchange_state.f90
!! Cross-compartment exchange records (T2-A / ADR 0053).
!!
!! One typed, directional record per coupling surface — the single channel a
!! compartment talks through, instead of reaching into a sibling's `state%X`
!! fields. Aggregated under `swap_exchange_t`, held as `state%exchange`. Copy
!! semantics: a record carries the crossing *values* (populated at the boundary),
!! so the surface is a real contract and an external component (BMI/XMI, no shared
!! memory) can be a drop-in producer/consumer.
!!
!! Grown incrementally, surface by surface, each cutover byte-identical. See
!! dev-docs/2026-07-03-arc-t2a-exchange-records-design.md for the full field set.
module exchange_state_mod
   use iso_fortran_env, only: real64
   implicit none
   private

   public :: crop_water_exchange_t, heat_soil_exchange_t, drain_soil_exchange_t
   public :: atmos_soil_exchange_t, swap_exchange_t

   !> Crop <-> atmosphere + soil-water (the transpiration/water-relations
   !! surface). The full crop<->SWAP interface an external crop model must
   !! satisfy: it feeds SWAP's ET + root sink and reads back the water-stress
   !! signal. Populated so far: the stress feedback (potential + actual
   !! transpiration); the ET/canopy inputs and root sink are added in subsequent
   !! byte-identical steps.
   integer, parameter :: ROOTDENS_LEN = 202   !! matches crop%common%cumdens

   type :: crop_water_exchange_t
      ! --- from crop: canopy state for the ET step (set at day start; stable per day) ---
      real(real64) :: lai            = 0.0_real64          !! leaf area index (-)
      real(real64) :: crop_height    = 0.0_real64          !! ch (cm)
      real(real64) :: interc_demand  = 0.0_real64          !! gird — gross irrigation depth (cm)
      real(real64) :: co2_transp_fac = 1.0_real64          !! fco2tra — CO2 transpiration factor (-)
      ! --- from crop: root distribution (set before the sink; stable per day) ---
      real(real64) :: rooting_depth = 0.0_real64          !! rd (cm)
      integer      :: root_nodes    = 0                   !! noddrz
      real(real64) :: root_density(ROOTDENS_LEN) = 0.0_real64  !! cumdens (afgen table)
      ! --- from atmosphere: potential transpiration (set by the ET step) ---
      real(real64) :: pot_transp = 0.0_real64   !! ptra (cm/d)
      ! --- from soil-water: actual transpiration (finalized at day end) ---
      real(real64) :: act_transp = 0.0_real64   !! tra  (cm/d)
   end type crop_water_exchange_t

   !> Heat <-> soil-water. Heat provides the frozen-fraction / soil-temperature
   !! for the hydraulic-conductivity adjustment read in the Richards solve;
   !! soil-water provides the water content for the de Vries thermal properties.
   type :: heat_soil_exchange_t
      ! --- from heat -> soil-water (read in the Richards solve) ---
      real(real64), allocatable :: rfcp(:)        !! reduced-frozen-fraction on K
      real(real64), allocatable :: tsoil(:)       !! soil temperature (degC)
      ! --- from soil-water -> heat (de Vries thermal props) ---
      real(real64), allocatable :: theta(:)       !! current water content
      real(real64), allocatable :: theta_prev(:)  !! thetm1 (previous step)
   end type heat_soil_exchange_t

   !> Drainage <-> soil-water. Soil-water provides the groundwater level + ponding
   !! that drive the lateral-drainage computation; drainage provides the per-level
   !! drainage flux (the sink applied in the soil-water flux balance).
   type :: drain_soil_exchange_t
      ! --- from soil-water -> drainage ---
      real(real64) :: gwl  = 0.0_real64   !! groundwater level
      real(real64) :: pond = 0.0_real64   !! surface ponding
      ! --- from drainage -> soil-water (the sink) ---
      integer                   :: nrlevs = 0   !! number of drainage levels
      real(real64), allocatable :: qdra(:,:)    !! drainage flux per level/node
   end type drain_soil_exchange_t

   !> Atmosphere <-> soil-water (top boundary). Atmosphere provides the net
   !! surface flux (rain/irrigation/snowmelt + potential evaporation) that drives
   !! the top-boundary condition; soil-water provides the ponding that feeds back
   !! into the potential-evaporation adjustment. (The atmo%cumu/intr water-balance
   !! accumulation is output accounting, NOT coupling — excluded.)
   type :: atmos_soil_exchange_t
      ! --- from atmosphere -> soil-water (top-boundary flux) ---
      real(real64) :: net_rain = 0.0_real64   !! nraidt
      real(real64) :: irrig    = 0.0_real64   !! nird
      real(real64) :: snowmelt = 0.0_real64   !! melt
      real(real64) :: pot_evap = 0.0_real64   !! peva
      real(real64) :: emp_evap = 0.0_real64   !! empreva
      integer      :: evap_reduce_method = 0  !! swredu
      ! --- from soil-water -> atmosphere ---
      real(real64) :: pond = 0.0_real64       !! ponding (peva adjustment)
   end type atmos_soil_exchange_t

   !> Aggregate of the cross-compartment exchange records. One field per surface.
   type :: swap_exchange_t
      type(crop_water_exchange_t) :: crop_water
      type(heat_soil_exchange_t)  :: heat_soil
      type(drain_soil_exchange_t) :: drain_soil
      type(atmos_soil_exchange_t) :: atmos_soil
   end type swap_exchange_t

end module exchange_state_mod
