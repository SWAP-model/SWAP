!> Test helper: chdir, getcwd, and file staging for pFUnit tests that
!! need to set up a case directory before calling legacy readers.
!!
!! Lives under tests/ — never shipped with the production binary.
module chdir_helper_mod
   use iso_c_binding, only: c_char, c_int, c_size_t, c_null_char
   implicit none
   private

   public :: chdir_to
   public :: get_cwd
   public :: stage_swp_template

   interface
      function c_chdir(path) bind(C, name="chdir") result(stat)
         import :: c_char, c_int
         character(kind=c_char), dimension(*), intent(in) :: path
         integer(c_int) :: stat
      end function c_chdir

      function c_getcwd(buf, size) bind(C, name="getcwd") result(res)
         import :: c_char, c_size_t, c_int
         character(kind=c_char), dimension(*), intent(out) :: buf
         integer(c_size_t), value, intent(in) :: size
         integer(c_int) :: res     ! actually a char*; nonzero on success
      end function c_getcwd
   end interface

contains

   !> Change working directory. Aborts (error stop) on failure to keep
   !! the test setup deterministic.
   subroutine chdir_to(path)
      character(len=*), intent(in) :: path
      integer(c_int) :: stat
      stat = c_chdir(trim(path) // c_null_char)
      if (stat /= 0) then
         write(*,'(A)') "chdir_to FAILED for: " // trim(path)
         error stop "chdir_helper: chdir failed"
      end if
   end subroutine chdir_to

   !> Read the current working directory into the provided buffer.
   !! Buffer should be large (e.g. 1024 chars).
   subroutine get_cwd(buf)
      character(len=*), intent(out) :: buf
      character(kind=c_char) :: cbuf(len(buf) + 1)
      integer(c_int) :: res
      integer :: i, n
      buf = ""
      cbuf = c_null_char
      res = c_getcwd(cbuf, int(len(buf), c_size_t))
      if (res == 0) then
         error stop "chdir_helper: getcwd failed"
      end if
      n = len(buf)
      do i = 1, n
         if (cbuf(i) == c_null_char) exit
         buf(i:i) = cbuf(i)
      end do
   end subroutine get_cwd

   !> Copy `<src_template>` to `<dst_basename>.swp` in the current
   !! working directory. The legacy SWAP convention uses
   !! `swap_linux.swp.template` -> `swap.swp` (so dst_basename = "swap").
   subroutine stage_swp_template(src_template, dst_basename)
      character(len=*), intent(in) :: src_template
      character(len=*), intent(in) :: dst_basename

      integer :: src_unit, dst_unit, ios
      character(len=4096) :: line
      character(len=:), allocatable :: dst_path

      dst_path = trim(dst_basename) // ".swp"

      open(newunit=src_unit, file=trim(src_template), status='old', &
           action='read', iostat=ios)
      if (ios /= 0) then
         write(*,'(A)') "stage_swp_template: cannot open source: " // trim(src_template)
         error stop "chdir_helper: stage source open failed"
      end if

      open(newunit=dst_unit, file=dst_path, status='replace', &
           action='write', iostat=ios)
      if (ios /= 0) then
         close(src_unit)
         write(*,'(A)') "stage_swp_template: cannot open dest: " // dst_path
         error stop "chdir_helper: stage dest open failed"
      end if

      do
         read(src_unit, '(A)', iostat=ios) line
         if (ios /= 0) exit
         write(dst_unit, '(A)') trim(line)
      end do

      close(src_unit)
      close(dst_unit)
   end subroutine stage_swp_template

end module chdir_helper_mod
