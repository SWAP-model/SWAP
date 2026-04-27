!> Reader for the [bottom_boundary] section of swap.toml.
!!
!! Populates a bottom_boundary_config_t from a fully-formed swap.toml
!! document. Tables (`swc_table`, `qbot_table`, `cofqha_table`) are encoded
!! as TOML arrays of arrays of two reals and decoded into 2-D allocatable
!! real(real64) arrays.
!!
!! Mechanical pattern (mirrors read_cropwofost_toml): scalars read one per
!! line, tables one block per table. No switch-conditional skipping — the
!! validator handles cross-field consistency. If [bottom_boundary] is
!! absent the reader returns silently and leaves `config` at defaults.
module read_bottom_boundary_toml_mod
   use iso_fortran_env, only: real64
   use tomlf, only: toml_table, toml_array, get_value, len
   use bottom_boundary_config_mod, only: bottom_boundary_config_t
   use toml_field_helpers_mod, only: get_table,                         &
                                     get_optional_int_with_default,     &
                                     get_optional_real_with_default,    &
                                     get_optional_string_with_default
   use error_mod, only: error_collection_t, ERR_PARSE_TYPE_MISMATCH
   implicit none
   private

   public :: read_bottom_boundary_toml

contains

   subroutine read_bottom_boundary_toml(doc_root, config, errors)
      type(toml_table), pointer,       intent(in)    :: doc_root
      type(bottom_boundary_config_t),  intent(inout) :: config
      type(error_collection_t),        intent(inout) :: errors

      type(toml_table), pointer :: sec

      call get_table(doc_root, 'bottom_boundary', sec, 'bottom_boundary', errors)
      if (.not. associated(sec)) return

      ! Always-present switch.
      call get_optional_int_with_default(sec, 'swbotb', config%swbotb, 0, &
                                         'bottom_boundary.swbotb', errors)

      ! SWBOTB=1 external file reference.
      call get_optional_string_with_default(sec, 'bbcfil', config%bbcfil, '', &
                                            'bottom_boundary.bbcfil', errors)

      ! SWBOTB=3 (Cauchy) scalars.
      call get_optional_int_with_default(sec,  'shape',  config%shape,  0, &
                                         'bottom_boundary.shape',  errors)
      call get_optional_real_with_default(sec, 'hdrain', config%hdrain, 0.0_real64, &
                                          'bottom_boundary.hdrain', errors)
      call get_optional_real_with_default(sec, 'rimlay', config%rimlay, 0.0_real64, &
                                          'bottom_boundary.rimlay', errors)
      call get_optional_real_with_default(sec, 'aqave',  config%aqave,  0.0_real64, &
                                          'bottom_boundary.aqave',  errors)
      call get_optional_real_with_default(sec, 'aqamp',  config%aqamp,  0.0_real64, &
                                          'bottom_boundary.aqamp',  errors)
      call get_optional_real_with_default(sec, 'aqper',  config%aqper,  0.0_real64, &
                                          'bottom_boundary.aqper',  errors)
      call get_optional_real_with_default(sec, 'aqtmax', config%aqtmax, 0.0_real64, &
                                          'bottom_boundary.aqtmax', errors)

      ! SWBOTB=5 scalars.
      call get_optional_real_with_default(sec, 'hbot',   config%hbot,   0.0_real64, &
                                          'bottom_boundary.hbot',   errors)
      call get_optional_real_with_default(sec, 'rhobot', config%rhobot, 0.0_real64, &
                                          'bottom_boundary.rhobot', errors)

      ! Tables.
      call read_table_2d(sec, 'swc_table',    config%swc_table,    2, &
                         'bottom_boundary.swc_table',    errors)
      call read_table_2d(sec, 'qbot_table',   config%qbot_table,   2, &
                         'bottom_boundary.qbot_table',   errors)
      call read_table_2d(sec, 'cofqha_table', config%cofqha_table, 2, &
                         'bottom_boundary.cofqha_table', errors)
   end subroutine read_bottom_boundary_toml

   !> Decode a TOML array-of-arrays at sec[key] into a (nrows, ncols)
   !! real(real64) allocatable. Absent key leaves table unallocated.
   !! Empty array (`key = []`) yields a 0-row allocation. Ragged or
   !! wrong-width inner arrays append a parse-type-mismatch error and
   !! leave the table unallocated.
   !!
   !! Local copy of the helper from read_cropwofost_toml — the helper is
   !! private to that module, so duplicating here keeps the readers
   !! decoupled (Phase 4d Task 4).
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

end module read_bottom_boundary_toml_mod
