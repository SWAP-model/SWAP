!> Reader for the [irrigation] section of swap.toml plus a per-crop
!! [irrigation_schedule] sub-section helper.
!!
!! Two public entry points (Phase 4d Task 11):
!!   * `read_irrigation_toml(doc_root, config, errors)` — top-level
!!     fixed-irrigation block (.swp side). Reads `swirfix`, `irgfil`,
!!     and the optional 4-column `fixed_events` table.
!!   * `read_irrigation_schedule_from_section(section, schedule, errors)`
!!     — fills an `irrigation_schedule_t` from a per-crop
!!     `[irrigation_schedule]` table (.crp side). The per-crop readers
!!     will call this helper once the `schedule` field is wired onto
!!     each crop type in Task 12.
!!
!! Mechanical pattern (mirrors read_bottom_boundary_toml /
!! read_heat_toml): scalars one per line, tables one block per table.
!! No switch-conditional skipping — the validator owns required-ness.
!! If the relevant section is absent, the routine returns silently and
!! leaves `config` / `schedule` at defaults (swirfix=0, schedule=0).
module read_irrigation_toml_mod
   use iso_fortran_env, only: real64
   use tomlf, only: toml_table, toml_array, get_value, len
   use irrigation_config_mod, only: irrigation_config_t, irrigation_schedule_t
   use toml_field_helpers_mod, only: get_table,                         &
                                     get_optional_int_with_default,     &
                                     get_optional_real_with_default,    &
                                     get_optional_string_with_default
   use error_mod, only: error_collection_t, ERR_PARSE_TYPE_MISMATCH
   implicit none
   private

   public :: read_irrigation_toml
   public :: read_irrigation_schedule_from_section

contains

   subroutine read_irrigation_toml(doc_root, config, errors)
      type(toml_table), pointer,   intent(in)    :: doc_root
      type(irrigation_config_t),   intent(inout) :: config
      type(error_collection_t),    intent(inout) :: errors

      type(toml_table), pointer :: sec

      call get_table(doc_root, 'irrigation', sec, 'irrigation', errors)
      if (.not. associated(sec)) return

      ! Always-present switch.
      call get_optional_int_with_default(sec, 'swirfix', config%swirfix, 0, &
                                         'irrigation.swirfix', errors)

      ! External fixed-irrigation file (alternative to inline table).
      call get_optional_string_with_default(sec, 'irgfil', config%irgfil, '', &
                                            'irrigation.irgfil', errors)

      ! Inline fixed-events: 4 columns (date / depth / conc / type).
      call read_table_2d(sec, 'fixed_events', config%fixed_events, 4, &
                         'irrigation.fixed_events', errors)
   end subroutine read_irrigation_toml

   !> Populate an `irrigation_schedule_t` from a per-crop
   !! `[irrigation_schedule]` sub-section. Caller resolves the section
   !! pointer (typically `get_table(crop_section, 'irrigation_schedule', ...)`)
   !! and may pass an unassociated pointer to signal absence — in that
   !! case the routine returns silently with `schedule` at defaults.
   subroutine read_irrigation_schedule_from_section(section, schedule, errors)
      type(toml_table), pointer,    intent(in)    :: section
      type(irrigation_schedule_t),  intent(inout) :: schedule
      type(error_collection_t),     intent(inout) :: errors

      if (.not. associated(section)) return

      ! Switches and timing-/depth-criterion enums.
      call get_optional_int_with_default(section, 'schedule',      schedule%schedule,      0, &
                                         'irrigation_schedule.schedule', errors)
      call get_optional_int_with_default(section, 'startirr_day',  schedule%startirr_day,  0, &
                                         'irrigation_schedule.startirr_day', errors)
      call get_optional_int_with_default(section, 'startirr_month', schedule%startirr_month, 0, &
                                         'irrigation_schedule.startirr_month', errors)
      call get_optional_int_with_default(section, 'endirr_day',    schedule%endirr_day,    0, &
                                         'irrigation_schedule.endirr_day', errors)
      call get_optional_int_with_default(section, 'endirr_month',  schedule%endirr_month,  0, &
                                         'irrigation_schedule.endirr_month', errors)
      call get_optional_real_with_default(section, 'cirrs',        schedule%cirrs,         0.0_real64, &
                                          'irrigation_schedule.cirrs', errors)
      call get_optional_int_with_default(section, 'isuas',         schedule%isuas,         0, &
                                         'irrigation_schedule.isuas', errors)
      call get_optional_int_with_default(section, 'tcs',           schedule%tcs,           0, &
                                         'irrigation_schedule.tcs', errors)
      call get_optional_int_with_default(section, 'dcs',           schedule%dcs,           0, &
                                         'irrigation_schedule.dcs', errors)

      ! tcs-conditional tables (validator enforces required-ness per branch).
      call read_table_2d(section, 'trel_table', schedule%trel_table, 2, &
                         'irrigation_schedule.trel_table', errors)
      call read_table_2d(section, 'raw_table',  schedule%raw_table,  2, &
                         'irrigation_schedule.raw_table',  errors)
      call read_table_2d(section, 'taw_table',  schedule%taw_table,  2, &
                         'irrigation_schedule.taw_table',  errors)
      call read_table_2d(section, 'dwa_table',  schedule%dwa_table,  2, &
                         'irrigation_schedule.dwa_table',  errors)
      call read_table_2d(section, 'hcri_table', schedule%hcri_table, 2, &
                         'irrigation_schedule.hcri_table', errors)
      call read_table_2d(section, 'tcri_table', schedule%tcri_table, 2, &
                         'irrigation_schedule.tcri_table', errors)

      ! Threshold / overrun scalars.
      call get_optional_real_with_default(section, 'irgthreshold',    schedule%irgthreshold,    0.0_real64, &
                                          'irrigation_schedule.irgthreshold', errors)
      call get_optional_real_with_default(section, 'dcrit',           schedule%dcrit,           0.0_real64, &
                                          'irrigation_schedule.dcrit', errors)
      call get_optional_int_with_default(section,  'swcirrthres',     schedule%swcirrthres,     0, &
                                         'irrigation_schedule.swcirrthres', errors)
      call get_optional_real_with_default(section, 'cirrthres',       schedule%cirrthres,       0.0_real64, &
                                          'irrigation_schedule.cirrthres', errors)
      call get_optional_real_with_default(section, 'perirrsurp',      schedule%perirrsurp,      0.0_real64, &
                                          'irrigation_schedule.perirrsurp', errors)
      call get_optional_int_with_default(section,  'tcsfix',          schedule%tcsfix,          0, &
                                         'irrigation_schedule.tcsfix', errors)
      call get_optional_int_with_default(section,  'irgdayfix',       schedule%irgdayfix,       0, &
                                         'irrigation_schedule.irgdayfix', errors)
      call get_optional_real_with_default(section, 'phfieldcapacity', schedule%phfieldcapacity, 0.0_real64, &
                                          'irrigation_schedule.phfieldcapacity', errors)

      ! dcs-conditional tables.
      call read_table_2d(section, 'di_table',  schedule%di_table,  2, &
                         'irrigation_schedule.di_table',  errors)
      call read_table_2d(section, 'fid_table', schedule%fid_table, 2, &
                         'irrigation_schedule.fid_table', errors)

      call get_optional_real_with_default(section, 'raithreshold', schedule%raithreshold, 0.0_real64, &
                                          'irrigation_schedule.raithreshold', errors)
      call get_optional_int_with_default(section,  'dcslim',       schedule%dcslim,       0, &
                                         'irrigation_schedule.dcslim', errors)
      call get_optional_real_with_default(section, 'irgdepmin',    schedule%irgdepmin,    0.0_real64, &
                                          'irrigation_schedule.irgdepmin', errors)
      call get_optional_real_with_default(section, 'irgdepmax',    schedule%irgdepmax,    0.0_real64, &
                                          'irrigation_schedule.irgdepmax', errors)
   end subroutine read_irrigation_schedule_from_section

   !> Decode a TOML array-of-arrays at sec[key] into a (nrows, ncols)
   !! real(real64) allocatable. Absent key leaves table unallocated.
   !! Empty array (`key = []`) yields a 0-row allocation. Ragged or
   !! wrong-width inner arrays append a parse-type-mismatch error and
   !! leave the table unallocated.
   !!
   !! Local copy of the helper from read_bottom_boundary_toml /
   !! read_heat_toml — those copies are private to their modules, so
   !! duplicating here keeps the readers decoupled (Phase 4d Task 11).
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

end module read_irrigation_toml_mod
