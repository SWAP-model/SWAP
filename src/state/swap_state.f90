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
   implicit none
   private
   public :: swap_state_t

   type :: swap_state_t
      type(surfacewater_state_t) :: surfacewater
      ! Subsequent migration arcs add fields here.
   end type swap_state_t

end module swap_state_mod
