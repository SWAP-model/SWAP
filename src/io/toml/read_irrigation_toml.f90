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
   use tomlf, only: toml_table, toml_array, toml_datetime, get_value, len
   use irrigation_config_mod, only: irrigation_config_t, irrigation_schedule_t
   use read_irrigation_ssdi_toml_mod, only: read_irrigation_ssdi_toml
   use toml_field_helpers_mod, only: get_table, get_array_of_tables,    &
                                     get_optional_int_with_default,     &
                                     get_optional_real_with_default,    &
                                     get_optional_string_with_default,  &
                                     parse_date_to_days1900
   use toml_array_helpers_mod, only: read_table_2d
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

      ! Phase 4f cleanup: CSV companion file (preferred replacement for
      ! the legacy .irg path). Path is relative to swap.toml and is
      ! resolved + read by the strangler adapter, not at parse time.
      call get_optional_string_with_default(sec, 'fixed_events_file', &
                                            config%fixed_events_file, '', &
                                            'irrigation.fixed_events_file', errors)

      ! Inline fixed-events: array-of-tables with named keys
      ! (date / depth / conc / type). Stored internally as a (nrows, 4)
      ! real(real64) table where col 1 is days-since-1900 (legacy axis)
      ! so the validator and downstream adapter remain shape-stable.
      call read_fixed_events(sec, config%fixed_events, errors)

      ! Delegate [irrigation.ssdi] parsing to the sibling reader.
      call read_irrigation_ssdi_toml(sec, config%ssdi, errors)
   end subroutine read_irrigation_toml

   !> Decode `[[irrigation.fixed_events]]` array-of-tables into the legacy
   !! (nrows, 4) real(real64) layout. Each row's `date` field is a TOML
   !! date literal that gets converted to days-since-1900 via
   !! `parse_date_to_days1900`. Missing keys append a parse error and
   !! leave the table unallocated (matches the legacy malformed-table
   !! semantics from `read_table_2d`).
   subroutine read_fixed_events(sec, table, errors)
      type(toml_table), pointer, intent(in)    :: sec
      real(real64), allocatable, intent(out)   :: table(:,:)
      type(error_collection_t),  intent(inout) :: errors

      type(toml_array), pointer :: arr
      type(toml_table), pointer :: row
      type(toml_datetime)       :: dtv
      integer  :: i, n, stat
      real(real64) :: depth_val, conc_val
      integer  :: type_val

      if (.not. associated(sec)) return

      call get_array_of_tables(sec, 'fixed_events', arr, &
                               'irrigation.fixed_events', errors)
      if (.not. associated(arr)) return

      n = len(arr)
      if (n == 0) then
         allocate(table(0, 4))
         return
      end if

      allocate(table(n, 4))
      table = 0.0_real64

      do i = 1, n
         row => null()
         call get_value(arr, i, row, stat=stat)
         if (stat /= 0 .or. .not. associated(row)) then
            call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                               "row not a table", 'irrigation.fixed_events')
            if (allocated(table)) deallocate(table)
            return
         end if

         call get_value(row, 'date', dtv, stat=stat)
         if (stat /= 0) then
            call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                               "missing or non-date 'date' key", &
                               'irrigation.fixed_events')
            if (allocated(table)) deallocate(table)
            return
         end if
         table(i, 1) = parse_date_to_days1900(dtv)

         call get_value(row, 'depth', depth_val, stat=stat)
         if (stat /= 0) then
            call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                               "missing or non-real 'depth' key", &
                               'irrigation.fixed_events')
            if (allocated(table)) deallocate(table)
            return
         end if
         table(i, 2) = depth_val

         call get_value(row, 'conc', conc_val, stat=stat)
         if (stat /= 0) then
            call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                               "missing or non-real 'conc' key", &
                               'irrigation.fixed_events')
            if (allocated(table)) deallocate(table)
            return
         end if
         table(i, 3) = conc_val

         call get_value(row, 'type', type_val, stat=stat)
         if (stat /= 0) then
            call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                               "missing or non-int 'type' key", &
                               'irrigation.fixed_events')
            if (allocated(table)) deallocate(table)
            return
         end if
         table(i, 4) = real(type_val, kind=real64)
      end do
   end subroutine read_fixed_events

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

end module read_irrigation_toml_mod
