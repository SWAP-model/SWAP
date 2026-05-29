!> Shared array-decoding helpers for the TOML reader family.
!!
!! Centralises three primitives that were verbatim-duplicated across
!! the TOML readers:
!!
!!   read_table_2d        -- TOML array-of-arrays  -> real(real64)(:,:)
!!   read_real_array_1d   -- flat TOML real array   -> real(real64)(:)
!!   read_int_array_1d    -- flat TOML integer array -> integer(:)
!!   read_real_pair_array -- flat even-length real array (dvs,value pairs) -> real(real64)(:)
!!
!! Formerly each reader owned a private copy.  Wiring in Phase 4f
!! (polish arc W9) eliminates the duplication.
module toml_array_helpers_mod
   use iso_fortran_env, only: real64
   use tomlf, only: toml_table, toml_array, get_value, len
   use error_mod, only: error_collection_t, ERR_PARSE_TYPE_MISMATCH, ERR_PARSE_ROW_SHAPE
   implicit none
   private

   public :: read_table_2d
   public :: read_real_array_1d
   public :: read_int_array_1d
   public :: read_real_pair_array

contains

   !> Decode a TOML array-of-arrays at sec[key] into a (nrows, ncols)
   !! real(real64) allocatable. Absent key leaves table unallocated.
   !! Empty array (`key = []`) yields a 0-row allocation. Ragged or
   !! wrong-width inner arrays append a parse-type-mismatch error and
   !! leave the table unallocated.
   subroutine read_table_2d(sec, key, table, expected_cols, context, errors)
      type(toml_table), pointer, intent(in)    :: sec
      character(len=*),          intent(in)    :: key
      real(real64), allocatable, intent(out)   :: table(:,:)
      integer,                   intent(in)    :: expected_cols
      character(len=*),          intent(in)    :: context
      type(error_collection_t),  intent(inout) :: errors

      type(toml_array), pointer :: outer, inner
      integer :: nrows, i, j, stat, n_inner
      real(real64) :: val

      if (.not. associated(sec)) return

      outer => null()
      call get_value(sec, key, outer, requested=.false., stat=stat)
      if (.not. associated(outer)) return

      nrows = len(outer)
      if (nrows == 0) then
         allocate(table(0, expected_cols))
         return
      end if

      allocate(table(nrows, expected_cols))
      table = 0.0_real64

      do i = 1, nrows
         inner => null()
         call get_value(outer, i, inner, stat=stat)
         if (stat /= 0 .or. .not. associated(inner)) then
            call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                               "row not an array", context)
            if (allocated(table)) deallocate(table)
            return
         end if
         n_inner = len(inner)
         if (n_inner /= expected_cols) then
            call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                               "row width mismatch", context)
            if (allocated(table)) deallocate(table)
            return
         end if
         do j = 1, expected_cols
            call get_value(inner, j, val, stat=stat)
            if (stat /= 0) then
               call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                                  "non-real cell", context)
               if (allocated(table)) deallocate(table)
               return
            end if
            table(i, j) = val
         end do
      end do
   end subroutine read_table_2d

   !> Decode a flat TOML real array at sec[key] into a 1-D real(real64)
   !! allocatable. Absent key leaves arr unallocated. Empty array
   !! (`key = []`) yields a 0-element allocation. Non-real cells append
   !! a parse-type-mismatch error and leave arr unallocated.
   subroutine read_real_array_1d(sec, key, arr, context, errors)
      type(toml_table), pointer, intent(in)    :: sec
      character(len=*),          intent(in)    :: key
      real(real64), allocatable, intent(out)   :: arr(:)
      character(len=*),          intent(in)    :: context
      type(error_collection_t),  intent(inout) :: errors

      type(toml_array), pointer :: outer
      integer :: n, i, stat
      real(real64) :: val

      if (.not. associated(sec)) return

      outer => null()
      call get_value(sec, key, outer, requested=.false., stat=stat)
      if (.not. associated(outer)) return

      n = len(outer)
      if (n == 0) then
         allocate(arr(0))
         return
      end if

      allocate(arr(n))
      arr = 0.0_real64

      do i = 1, n
         call get_value(outer, i, val, stat=stat)
         if (stat /= 0) then
            call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                               "non-real cell", context)
            if (allocated(arr)) deallocate(arr)
            return
         end if
         arr(i) = val
      end do
   end subroutine read_real_array_1d

   !> Decode a flat TOML integer array at sec[key] into a 1-D integer
   !! allocatable. Absent key leaves arr unallocated. Empty array yields
   !! a 0-element allocation. Non-integer cells append a parse-type-mismatch
   !! error and leave arr unallocated.
   subroutine read_int_array_1d(sec, key, arr, context, errors)
      type(toml_table), pointer, intent(in)    :: sec
      character(len=*),          intent(in)    :: key
      integer, allocatable,      intent(out)   :: arr(:)
      character(len=*),          intent(in)    :: context
      type(error_collection_t),  intent(inout) :: errors

      type(toml_array), pointer :: outer
      integer :: n, i, stat, val

      if (.not. associated(sec)) return

      outer => null()
      call get_value(sec, key, outer, requested=.false., stat=stat)
      if (.not. associated(outer)) return

      n = len(outer)
      if (n == 0) then
         allocate(arr(0))
         return
      end if

      allocate(arr(n))
      arr = 0

      do i = 1, n
         call get_value(outer, i, val, stat=stat)
         if (stat /= 0) then
            call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                               "non-int cell", context)
            if (allocated(arr)) deallocate(arr)
            return
         end if
         arr(i) = val
      end do
   end subroutine read_int_array_1d

   !> Decode a flat TOML real array at sec[key] into a 1-D real(real64)
   !! allocatable. Enforces even length (dvs, value) pair constraint.
   !! Absent key leaves arr unallocated. Odd-length arrays append a
   !! parse-row-shape error and leave arr unallocated. Non-real cells
   !! append a parse-type-mismatch error and leave arr unallocated.
   subroutine read_real_pair_array(sec, key, arr, context, errors)
      type(toml_table), pointer, intent(in)    :: sec
      character(len=*),          intent(in)    :: key
      real(real64), allocatable, intent(out)   :: arr(:)
      character(len=*),          intent(in)    :: context
      type(error_collection_t),  intent(inout) :: errors

      type(toml_array), pointer :: a
      integer :: stat, i, n
      real(real64) :: v

      if (.not. associated(sec)) return
      call get_value(sec, key, a, requested=.false., stat=stat)
      if (stat /= 0 .or. .not. associated(a)) return
      n = len(a)
      if (mod(n, 2) /= 0) then
         call errors%append(ERR_PARSE_ROW_SHAPE, &
            'expected even-length (dvs,value) pair array', context)
         return
      end if
      allocate(arr(n))
      do i = 1, n
         call get_value(a, i, v, stat=stat)
         if (stat /= 0) then
            call errors%append(ERR_PARSE_TYPE_MISMATCH, &
               'non-real cell in pair array', context)
            deallocate(arr)
            return
         end if
         arr(i) = v
      end do
   end subroutine read_real_pair_array

end module toml_array_helpers_mod
