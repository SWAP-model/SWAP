!> @file swap_c_strings_mod.f90
!! Shared C <-> Fortran string helpers for the C-ABI facades (arc T1-G′ sub-arc 2).
!!
!! Previously each facade (swap_bmi_mod / swap_capi_mod / swap_xmi_mod) carried its
!! own copies. `c_to_f_string` was identical in all three. `f_to_c_string` came in
!! two contracts — BMI's writes at most `max_len-1` chars then a NUL (bounded);
!! XMI's writes `len_trim` chars then a NUL with no bound (the caller guarantees a
!! large-enough buffer). Both are kept, exposed through one generic `f_to_c_string`
!! that resolves by argument count, so every existing call site is unchanged.
module swap_c_strings_mod
   use iso_c_binding, only: c_char, c_null_char
   implicit none
   private

   public :: c_to_f_string, f_to_c_string

   interface f_to_c_string
      module procedure f_to_c_string_bounded
      module procedure f_to_c_string_unbounded
   end interface f_to_c_string

contains

   !> Copy a NUL-terminated C string into a blank-padded Fortran string.
   subroutine c_to_f_string(c_str, f_str)
      character(kind=c_char), intent(in)  :: c_str(*)
      character(len=*),       intent(out) :: f_str
      integer :: i
      f_str = ' '
      do i = 1, len(f_str)
         if (c_str(i) == c_null_char) exit
         f_str(i:i) = c_str(i)
      end do
   end subroutine c_to_f_string

   !> Bounded form (BMI): writes at most max_len-1 chars then a NUL terminator.
   subroutine f_to_c_string_bounded(f_str, c_buf, max_len)
      character(len=*),       intent(in)  :: f_str
      character(kind=c_char), intent(out) :: c_buf(*)
      integer,                intent(in)  :: max_len
      integer :: i, n
      n = min(len_trim(f_str), max_len - 1)
      do i = 1, n
         c_buf(i) = f_str(i:i)
      end do
      c_buf(n + 1) = c_null_char
   end subroutine f_to_c_string_bounded

   !> Unbounded form (XMI): writes len_trim chars then a NUL; the caller must
   !! provide a buffer large enough (matches xmipy's fixed-size buffers).
   subroutine f_to_c_string_unbounded(f_str, c_str)
      character(len=*),       intent(in)  :: f_str
      character(kind=c_char), intent(out) :: c_str(*)
      integer :: i
      do i = 1, len_trim(f_str)
         c_str(i) = f_str(i:i)
      end do
      c_str(len_trim(f_str) + 1) = c_null_char
   end subroutine f_to_c_string_unbounded

end module swap_c_strings_mod
