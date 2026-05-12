!> @file swap_mod.f90
!! SS-DRV Phase 1: module form of the legacy `subroutine swap`.
!! Three named lifecycle procedures replace the (iCaller, iTask) dispatch.
!! Time loop lives in the caller. State and config are threaded explicitly.
module swap_mod
   use swap_state_mod,  only: swap_state_t
   use swap_config_mod, only: swap_config_t
   implicit none
   private
   public :: swap_init, swap_run_step, swap_close

contains

   subroutine swap_init(config_file, state, config)
      character(len=*),            intent(in)  :: config_file
      type(swap_state_t),          intent(out) :: state
      type(swap_config_t), target, intent(out) :: config
      ! Body filled in Task 3.
   end subroutine swap_init

   subroutine swap_run_step(state, config)
      type(swap_state_t),  intent(inout) :: state
      type(swap_config_t), intent(in)    :: config
      ! Body filled in Task 4.
   end subroutine swap_run_step

   subroutine swap_close(state, config)
      type(swap_state_t),  intent(inout) :: state
      type(swap_config_t), intent(in)    :: config
      ! Body filled in Task 5.
   end subroutine swap_close

end module swap_mod
