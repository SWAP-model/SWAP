!> State-aware wrapper module for legacy `CropGrowth`.
!!
!! Provides an explicit-state entry point while keeping the original
!! `CropGrowth` routine unchanged.
module cropgrowth_state_mod
   implicit none
   private

   public :: CropGrowth_state

contains

!> State-aware wrapper for `CropGrowth`.
!!
!! Executes the legacy task-based crop growth routine and synchronizes crop
!! outputs back into the explicit state container.
!!
!! @param[inout] state SWAP model state container
!! @param[in]    task  Legacy task selector
   subroutine CropGrowth_state(state, task)
      use swap_state_mod, only: swap_state_t
      use swap_state_sync, only: crop_state_from_variables
      implicit none

      type(swap_state_t), intent(inout) :: state
      integer,            intent(in)    :: task

      call CropGrowth(task)
      call crop_state_from_variables(state%crop, state%ncrop)
   end subroutine CropGrowth_state

end module cropgrowth_state_mod
