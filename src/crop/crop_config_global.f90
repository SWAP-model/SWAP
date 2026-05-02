!> Module-level pointer to the parsed crop config.
!!
!! Transitional element of the .crp port (ADR 0016). The runtime path
!! reads through this pointer to access the per-rotation cache stored
!! in `crop_config_t.rotation_fixed(:)` etc., without changing the
!! signatures of every legacy crop subroutine. After Phase 4 of the
!! .crp port and the config-passing follow-on spec, computation subs
!! take typed config + state as explicit arguments and this module
!! goes away entirely.
!!
!! Lifecycle: set by config_to_variables at the end of its crop block;
!! valid for the duration of the simulation (the underlying config
!! lives in a local of core/swap.f90's iTask=1 block, which encloses
!! all simulation initialization). Readers must check `associated(...)`
!! defensively.
module crop_config_global_mod
   use crop_config_mod, only: crop_config_t
   implicit none
   private

   public :: crop_config_global

   type(crop_config_t), pointer :: crop_config_global => null()
end module crop_config_global_mod
