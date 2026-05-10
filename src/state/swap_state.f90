!> @file swap_state.f90
!! Top-level aggregator of all per-subsystem typed state records.
!! Phase 1 of the state-migration arc starts with only the
!! surface-water subsystem; subsequent migration arcs (water-balance,
!! crop, drainage, atmosphere, …) add their own fields here.
!!
!! Threaded through subroutine signatures from `swap_main` down,
!! replacing the implicit shared state held in `variables.f90`.

module swap_state_mod
   use surfacewater_state_mod, only: surfacewater_state_t
   use drainage_state_mod,     only: drainage_state_t
   use solute_state_mod,       only: solute_state_t
   use heat_state_mod,         only: heat_state_t
   implicit none
   private
   public :: swap_state_t

   type :: swap_state_t
      type(surfacewater_state_t) :: surfacewater
      type(drainage_state_t)     :: drainage
      type(solute_state_t)       :: solute
      type(heat_state_t)         :: heat
   end type swap_state_t

end module swap_state_mod
