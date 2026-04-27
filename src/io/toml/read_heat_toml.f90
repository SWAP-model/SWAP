!> Reader for the [heat] section of swap.toml.
!!
!! Populates a heat_config_t from a fully-formed swap.toml document. The
!! per-layer texture fractions (`psand`, `pclay`, `porg`) are encoded as flat
!! TOML arrays of real, decoded into 1-D allocatable real(real64) arrays.
!! The initial soil-temperature table (`tsoil_init`) is encoded as a TOML
!! array of arrays of two reals (depth_cm, temp_C) and decoded into a
!! 2-D allocatable real(real64) array.
!!
!! Mechanical pattern (mirrors read_bottom_boundary_toml): scalars one per
!! line, arrays/tables one block per array. No switch-conditional skipping
!! — the validator handles cross-field consistency. If [heat] is absent
!! the reader returns silently and leaves `config` at defaults.
module read_heat_toml_mod
   use iso_fortran_env, only: real64
   use tomlf, only: toml_table, toml_array, get_value, len
   use heat_config_mod, only: heat_config_t
   use toml_field_helpers_mod, only: get_table,                         &
                                     get_optional_int_with_default,     &
                                     get_optional_real_with_default
   use error_mod, only: error_collection_t, ERR_PARSE_TYPE_MISMATCH
   implicit none
   private

   public :: read_heat_toml

contains

   subroutine read_heat_toml(doc_root, config, errors)
      type(toml_table), pointer, intent(in)    :: doc_root
      type(heat_config_t),       intent(inout) :: config
      type(error_collection_t),  intent(inout) :: errors

      type(toml_table), pointer :: sec

      call get_table(doc_root, 'heat', sec, 'heat', errors)
      if (.not. associated(sec)) return

      ! Switches.
      call get_optional_int_with_default(sec, 'swhea',     config%swhea,     0, &
                                         'heat.swhea',     errors)
      call get_optional_int_with_default(sec, 'swcalt',    config%swcalt,    0, &
                                         'heat.swcalt',    errors)
      call get_optional_int_with_default(sec, 'swtopbhea', config%swtopbhea, 0, &
                                         'heat.swtopbhea', errors)
      call get_optional_int_with_default(sec, 'swbotbhea', config%swbotbhea, 0, &
                                         'heat.swbotbhea', errors)

      ! Frost params (scalars).
      call get_optional_real_with_default(sec, 'tfroststa', config%tfroststa, 0.0_real64, &
                                          'heat.tfroststa', errors)
      call get_optional_real_with_default(sec, 'tfrostend', config%tfrostend, 0.0_real64, &
                                          'heat.tfrostend', errors)

      ! Per-layer texture fractions (1-D arrays).
      call read_array_1d(sec, 'psand', config%psand, 'heat.psand', errors)
      call read_array_1d(sec, 'pclay', config%pclay, 'heat.pclay', errors)
      call read_array_1d(sec, 'porg',  config%porg,  'heat.porg',  errors)

      ! Initial soil-temperature table (2-D, depth/temp pairs).
      call read_table_2d(sec, 'tsoil_init', config%tsoil_init, 2, &
                         'heat.tsoil_init', errors)
   end subroutine read_heat_toml

   !> Decode a flat TOML array at sec[key] into a 1-D real(real64)
   !! allocatable. Absent key leaves arr unallocated. Empty array
   !! (`key = []`) yields a 0-element allocation. Non-real cells append
   !! a parse-type-mismatch error and leave arr unallocated.
   subroutine read_array_1d(sec, key, arr, context, errors)
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
   end subroutine read_array_1d

   !> Decode a TOML array-of-arrays at sec[key] into a (nrows, ncols)
   !! real(real64) allocatable. Absent key leaves table unallocated.
   !! Empty array (`key = []`) yields a 0-row allocation. Ragged or
   !! wrong-width inner arrays append a parse-type-mismatch error and
   !! leave the table unallocated.
   !!
   !! Local copy of the helper from read_bottom_boundary_toml — that
   !! helper is private to its module, so duplicating here keeps the
   !! readers decoupled (Phase 4d Task 7).
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

end module read_heat_toml_mod
