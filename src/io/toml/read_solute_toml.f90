!> Reader for the [solute] section of swap.toml.
!!
!! Populates a solute_config_t from a fully-formed swap.toml document.
!! Scalars are read one per line; the per-layer decomposition table
!! (`pertabsolu`) is encoded as a TOML array of arrays of two reals
!! (depth, factor) and decoded into a 2-D allocatable real(real64)
!! array.
!!
!! Mechanical pattern (mirrors read_heat_toml / read_bottom_boundary_toml):
!! no switch-conditional skipping — the validator handles cross-field
!! consistency. If [solute] is absent the reader returns silently and
!! leaves `config` at defaults (the validator's `swsolu=0` sentinel
!! short-circuits cleanly in that case).
module read_solute_toml_mod
   use iso_fortran_env, only: real64
   use tomlf, only: toml_table, toml_array, get_value, len
   use solute_config_mod, only: solute_config_t
   use toml_field_helpers_mod, only: get_table,                         &
                                     get_optional_int_with_default,     &
                                     get_optional_real_with_default
   use error_mod, only: error_collection_t, ERR_PARSE_TYPE_MISMATCH
   implicit none
   private

   public :: read_solute_toml

contains

   subroutine read_solute_toml(doc_root, config, errors)
      type(toml_table), pointer, intent(in)    :: doc_root
      type(solute_config_t),     intent(inout) :: config
      type(error_collection_t),  intent(inout) :: errors

      type(toml_table), pointer :: sec

      call get_table(doc_root, 'solute', sec, 'solute', errors)
      if (.not. associated(sec)) return

      ! Switches.
      call get_optional_int_with_default(sec, 'swsolu',   config%swsolu,   0, &
                                         'solute.swsolu',   errors)
      call get_optional_int_with_default(sec, 'swbotbc',  config%swbotbc,  0, &
                                         'solute.swbotbc',  errors)
      call get_optional_int_with_default(sec, 'swsoltyp', config%swsoltyp, 0, &
                                         'solute.swsoltyp', errors)
      call get_optional_int_with_default(sec, 'swdc',     config%swdc,     0, &
                                         'solute.swdc',     errors)

      ! Concentrations and dispersion.
      call get_optional_real_with_default(sec, 'cdrain', config%cdrain, 0.0_real64, &
                                          'solute.cdrain', errors)
      call get_optional_real_with_default(sec, 'cseep',  config%cseep,  0.0_real64, &
                                          'solute.cseep',  errors)
      call get_optional_real_with_default(sec, 'tscf',   config%tscf,   0.0_real64, &
                                          'solute.tscf',   errors)
      call get_optional_real_with_default(sec, 'ldis',   config%ldis,   0.0_real64, &
                                          'solute.ldis',   errors)

      ! Root-uptake.
      call get_optional_real_with_default(sec, 'rtheta', config%rtheta, 0.0_real64, &
                                          'solute.rtheta', errors)
      call get_optional_real_with_default(sec, 'bexp',   config%bexp,   0.0_real64, &
                                          'solute.bexp',   errors)

      ! Salinity scalars.
      call get_optional_real_with_default(sec, 'ecmax',  config%ecmax,  0.0_real64, &
                                          'solute.ecmax',  errors)
      call get_optional_real_with_default(sec, 'ecslop', config%ecslop, 0.0_real64, &
                                          'solute.ecslop', errors)

      ! Per-layer decomposition table.
      call read_table_2d(sec, 'pertabsolu', config%pertabsolu, 2, &
                         'solute.pertabsolu', errors)
   end subroutine read_solute_toml

   !> Decode a TOML array-of-arrays at sec[key] into a (nrows, ncols)
   !! real(real64) allocatable. Absent key leaves table unallocated.
   !! Empty array (`key = []`) yields a 0-row allocation. Ragged or
   !! wrong-width inner arrays append a parse-type-mismatch error and
   !! leave the table unallocated.
   !!
   !! Local copy of the helper from read_bottom_boundary_toml /
   !! read_heat_toml — those helpers are private to their modules, so
   !! duplicating here keeps the readers decoupled (per Phase 4d
   !! pattern; see Task 7's note on read_heat_toml).
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

end module read_solute_toml_mod
