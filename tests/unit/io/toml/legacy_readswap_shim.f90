!> Thin wrapper around the legacy readswap entry point, used only by
!! parity tests.  Takes a case directory and populates a fresh
!! swap_state_t via the legacy code path.
!!
!! NOT intended for production.  Lives under tests/ so it never ships.
!!
!! Implementation notes (Task 26 exploration):
!!
!!   The legacy readswap subroutine reads its .swp filename from
!!   Get_Command_Argument(1), not from a parameter.  A direct call
!!   therefore requires a workaround:
!!
!!     1. chdir into case_dir so relative file opens work.
!!     2. The command-line argument is already set by the test runner
!!        (or by a pre-call to set_command_argument if available), OR
!!     3. Set the global `swpfile` / `swpfilnam` directly before calling
!!        key sub-routines (more invasive).
!!
!!   For Task 27 (parity test) both approaches are viable; the JSON-
!!   snapshot fallback is the safer default if command-argument wiring
!!   proves fragile under pFUnit.
module legacy_readswap_shim_mod
   use swap_state_mod, only: swap_state_t, swap_state_init
   implicit none
   private
   public :: legacy_readswap_case

contains

   !> Call the legacy readswap reader for *case_dir*/*swp_filename* and
   !! populate *state*.
   !!
   !! Current status: **stub** — populates a minimal state so the shim
   !! compiles and the parity-test scaffold can be written (Task 27).
   !! Full wiring (chdir + readswap + state_from_variables) is deferred
   !! until the command-argument workaround is settled.
   subroutine legacy_readswap_case(case_dir, swp_filename, state)
      character(len=*),   intent(in)    :: case_dir
      character(len=*),   intent(in)    :: swp_filename
      type(swap_state_t), intent(inout) :: state

      ! Suppress unused-argument warnings while the stub is in place.
      ! Remove when real implementation is added.
      associate(cd => case_dir, sf => swp_filename)
      end associate

      ! Stub: initialise with placeholder dimensions so callers have a
      ! valid (if empty) state object.
      call swap_state_init(state, 1, 1)

      ! TODO Task 27: replace the stub with the real call sequence:
      !
      !   integer :: istat
      !   call chdir(trim(case_dir), istat)   ! or use iso_c_binding chdir
      !   call readswap(state)                 ! fills legacy `variables` globals
      !   call swap_state_init(state, numnod, numlay)
      !   call state_from_variables(state)    ! propagate globals -> state

   end subroutine legacy_readswap_case

end module legacy_readswap_shim_mod
