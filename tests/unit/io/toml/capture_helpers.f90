! tests/unit/io/toml/capture_helpers.f90
!
! Temporary instrumentation used during SS-A literal-value capture.
! Each parity test gets a few `call print_capture(...)` lines, runs
! once to emit literals to stderr, then both the lines and this
! module are deleted. The module exists only for the duration of SS-A.
module capture_helpers_mod
   use iso_fortran_env, only: real64
   implicit none
   private
   public :: print_capture_real, print_capture_int, print_capture_str

   character(len=*), parameter :: CAP_FILE = '/tmp/cap_literals.txt'
contains
   subroutine print_capture_real(name, val)
      character(len=*), intent(in) :: name
      real(real64),     intent(in) :: val
      integer :: u, ios
      open(newunit=u, file=CAP_FILE, status='unknown', position='append', &
           action='write', iostat=ios)
      if (ios == 0) then
         write(u, '("EXPECT ", a, " = ", es24.16)') name, val
         close(u)
      end if
   end subroutine
   subroutine print_capture_int(name, val)
      character(len=*), intent(in) :: name
      integer,          intent(in) :: val
      integer :: u, ios
      open(newunit=u, file=CAP_FILE, status='unknown', position='append', &
           action='write', iostat=ios)
      if (ios == 0) then
         write(u, '("EXPECT ", a, " = ", i0)') name, val
         close(u)
      end if
   end subroutine
   subroutine print_capture_str(name, val)
      character(len=*), intent(in) :: name
      character(len=*), intent(in) :: val
      integer :: u, ios
      open(newunit=u, file=CAP_FILE, status='unknown', position='append', &
           action='write', iostat=ios)
      if (ios == 0) then
         write(u, '("EXPECT ", a, " = ", a)') name, trim(val)
         close(u)
      end if
   end subroutine
end module
