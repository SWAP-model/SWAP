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

   public :: crop_water_exchange_t, swap_exchange_t

   !> Crop <-> atmosphere + soil-water (the transpiration/water-relations
   !! surface). The full crop<->SWAP interface an external crop model must
   !! satisfy: it feeds SWAP's ET + root sink and reads back the water-stress
   !! signal. Populated so far: the stress feedback (potential + actual
   !! transpiration); the ET/canopy inputs and root sink are added in subsequent
   !! byte-identical steps.
   type :: crop_water_exchange_t
      ! --- from atmosphere: potential transpiration (set by the ET step) ---
      real(real64) :: pot_transp = 0.0_real64   !! ptra (cm/d)
      ! --- from soil-water: actual transpiration (finalized at day end) ---
      real(real64) :: act_transp = 0.0_real64   !! tra  (cm/d)
   end type crop_water_exchange_t

   !> Aggregate of the cross-compartment exchange records. One field per surface.
   type :: swap_exchange_t
      type(crop_water_exchange_t) :: crop_water
   end type swap_exchange_t

end module exchange_state_mod
