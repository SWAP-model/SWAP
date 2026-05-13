!> @file timecontrol_mod.f90
!! SS-TCM: module form of the legacy subroutine TimeControl + IterTime.
!! Seven named lifecycle procedures replace the magic-int dispatch
!! (task=1/2/3/9 + IterTime task=1/2/3). State threaded explicitly;
!! each procedure carries its own associate block. flZeroIntr and
!! flZeroCumu are owned by state%timecontrol — bare globals
!! consumed by Task 11 readers retire in Task 12.
module timecontrol_mod
   use swap_state_mod, only: swap_state_t
   implicit none
   private
   public :: timecontrol_init, timecontrol_advance, &
             timecontrol_reduce_dt, timecontrol_day_end
   public :: itertime_init, itertime_check, itertime_close

contains

   subroutine timecontrol_init(state)
      type(swap_state_t), intent(inout) :: state
      ! Body filled in Task 5 (migrated from timecontrol.f90 case (1)).
   end subroutine timecontrol_init

   subroutine timecontrol_advance(state)
      type(swap_state_t), intent(inout) :: state
      ! Body filled in Task 6 (migrated from timecontrol.f90 case (2)).
   end subroutine timecontrol_advance

   subroutine timecontrol_reduce_dt(state)
      type(swap_state_t), intent(inout) :: state
      ! Body filled in Task 7 (migrated from timecontrol.f90 case (3)).
   end subroutine timecontrol_reduce_dt

   subroutine timecontrol_day_end(state)
      type(swap_state_t), intent(inout) :: state
      ! Body filled in Task 8 (migrated from timecontrol.f90 case (9)).
   end subroutine timecontrol_day_end

   subroutine itertime_init(state)
      type(swap_state_t), intent(inout) :: state
      ! Body filled in Task 9 (migrated from IterTime case (1)).
   end subroutine itertime_init

   subroutine itertime_check(state)
      type(swap_state_t), intent(inout) :: state
      ! Body filled in Task 9 (migrated from IterTime case (2)).
   end subroutine itertime_check

   subroutine itertime_close(state)
      type(swap_state_t), intent(inout) :: state
      ! Body filled in Task 9 (migrated from IterTime case (3)).
   end subroutine itertime_close

end module timecontrol_mod
