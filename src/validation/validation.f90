!> Shared validator primitives used by per-section config validators.
!!
!! Each subroutine checks one invariant and appends an error to the
!! passed `error_collection_t` on failure. No state. No side effects
!! other than the append (which auto-logs). Pass `context` so the
!! error identifies the offending field.
module validation_mod
   use iso_fortran_env, only: real64
   use error_mod, only: error_collection_t,          &
                        ERR_VALIDATION_OUT_OF_RANGE, &
                        ERR_VALIDATION_ENUM,         &
                        ERR_VALIDATION_CROSS_FIELD
   implicit none
   private

   public :: check_int_range
   public :: check_real_range
   public :: check_int_enum
   public :: check_not_empty
   public :: check_nonnegative_real
   public :: check_ordered_pair

contains

   !> Verify low <= value <= high.
   subroutine check_int_range(value, low, high, context, errors)
      integer,                   intent(in)    :: value, low, high
      character(len=*),          intent(in)    :: context
      type(error_collection_t),  intent(inout) :: errors
      character(len=64) :: msg
      if (value < low .or. value > high) then
         write(msg, '("value ",I0," outside [",I0,",",I0,"]")') value, low, high
         call errors%append(ERR_VALIDATION_OUT_OF_RANGE, trim(msg), context)
      end if
   end subroutine check_int_range

   !> Verify low <= value <= high (real64).
   subroutine check_real_range(value, low, high, context, errors)
      real(real64),              intent(in)    :: value, low, high
      character(len=*),          intent(in)    :: context
      type(error_collection_t),  intent(inout) :: errors
      character(len=128) :: msg
      if (value < low .or. value > high) then
         write(msg, '("value ",ES12.5," outside [",ES12.5,",",ES12.5,"]")') value, low, high
         call errors%append(ERR_VALIDATION_OUT_OF_RANGE, trim(msg), context)
      end if
   end subroutine check_real_range

   !> Verify value is one of the allowed integers.
   subroutine check_int_enum(value, allowed, context, errors)
      integer,                   intent(in)    :: value
      integer,                   intent(in)    :: allowed(:)
      character(len=*),          intent(in)    :: context
      type(error_collection_t),  intent(inout) :: errors
      character(len=128) :: msg
      integer :: i
      do i = 1, size(allowed)
         if (value == allowed(i)) return
      end do
      write(msg, '("value ",I0," not in allowed set")') value
      call errors%append(ERR_VALIDATION_ENUM, trim(msg), context)
   end subroutine check_int_enum

   !> Verify the string (trimmed) is non-empty.
   subroutine check_not_empty(value, context, errors)
      character(len=*),          intent(in)    :: value
      character(len=*),          intent(in)    :: context
      type(error_collection_t),  intent(inout) :: errors
      if (len_trim(value) == 0) then
         call errors%append(ERR_VALIDATION_OUT_OF_RANGE, "empty string", context)
      end if
   end subroutine check_not_empty

   !> Verify value >= 0.
   subroutine check_nonnegative_real(value, context, errors)
      real(real64),              intent(in)    :: value
      character(len=*),          intent(in)    :: context
      type(error_collection_t),  intent(inout) :: errors
      character(len=64) :: msg
      if (value < 0.0_real64) then
         write(msg, '("value ",ES12.5," is negative")') value
         call errors%append(ERR_VALIDATION_OUT_OF_RANGE, trim(msg), context)
      end if
   end subroutine check_nonnegative_real

   !> Verify low_val <= high_val (cross-field real64 check).
   subroutine check_ordered_pair(low_val, high_val, low_name, high_name, context, errors)
      real(real64),              intent(in)    :: low_val, high_val
      character(len=*),          intent(in)    :: low_name, high_name
      character(len=*),          intent(in)    :: context
      type(error_collection_t),  intent(inout) :: errors
      character(len=256) :: msg
      if (low_val > high_val) then
         write(msg, '(A," (",ES12.5,") exceeds ",A," (",ES12.5,")")') &
              low_name, low_val, high_name, high_val
         call errors%append(ERR_VALIDATION_CROSS_FIELD, trim(msg), context)
      end if
   end subroutine check_ordered_pair

end module validation_mod
