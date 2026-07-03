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

   public :: crop_water_exchange_t, heat_soil_exchange_t, swap_exchange_t

   !> Crop <-> atmosphere + soil-water (the transpiration/water-relations
   !! surface). The full crop<->SWAP interface an external crop model must
   !! satisfy: it feeds SWAP's ET + root sink and reads back the water-stress
   !! signal. Populated so far: the stress feedback (potential + actual
   !! transpiration); the ET/canopy inputs and root sink are added in subsequent
   !! byte-identical steps.
   integer, parameter :: ROOTDENS_LEN = 202   !! matches crop%common%cumdens

   type :: crop_water_exchange_t
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

   !> Aggregate of the cross-compartment exchange records. One field per surface.
   type :: swap_exchange_t
      type(crop_water_exchange_t) :: crop_water
      type(heat_soil_exchange_t)  :: heat_soil
   end type swap_exchange_t

end module exchange_state_mod
