!> @file read_irrigation_ssdi_toml.f90
!! Parses the optional [irrigation.ssdi] block into
!! irrigation_config_t%ssdi. Block omitted when irrigation.swssdi = 0;
!! reader leaves defaults intact in that case. Sibling to
!! read_irrigation_toml.f90 to keep that file focused.
module read_irrigation_ssdi_toml_mod

   use, intrinsic :: iso_fortran_env, only: real64
   use tomlf, only: toml_table, toml_array, get_value, len
   use toml_field_helpers_mod, only: get_optional_int_with_default,       &
                                     get_optional_real_with_default,      &
                                     get_optional_string_with_default
   use irrigation_config_mod, only: irrigation_ssdi_t,                    &
                                    irrigation_ssdi_fixed_t,              &
                                    irrigation_ssdi_scheduled_t
   use error_mod, only: error_collection_t, ERR_PARSE_TYPE_MISMATCH
   implicit none
   private

   public :: read_irrigation_ssdi_toml

contains

   !> Populate `ssdi` from the optional [irrigation.ssdi] table on `irr_sec`.
   !! Missing block is benign — leaves defaults.
   subroutine read_irrigation_ssdi_toml(irr_sec, ssdi, errors)
      type(toml_table), pointer, intent(in)    :: irr_sec
      type(irrigation_ssdi_t),   intent(inout) :: ssdi
      type(error_collection_t),  intent(inout) :: errors

      type(toml_table), pointer :: ssdi_tbl, fixed_tbl, sched_tbl
      type(toml_array), pointer :: zarr
      integer :: stat, n, k
      real(real64) :: zval

      if (.not. associated(irr_sec)) return

      ssdi_tbl => null()
      call get_value(irr_sec, 'ssdi', ssdi_tbl, requested=.false., stat=stat)
      if (stat /= 0 .or. .not. associated(ssdi_tbl)) return

      call get_optional_int_with_default(ssdi_tbl, 'schedule', ssdi%schedule, 0, &
                                         'irrigation.ssdi.schedule', errors)

      ! ssdi_z is a 2-element array; both elements equal for single-depth.
      zarr => null()
      call get_value(ssdi_tbl, 'ssdi_z', zarr, requested=.false., stat=stat)
      if (stat == 0 .and. associated(zarr)) then
         n = len(zarr)
         if (n == 1) then
            call get_value(zarr, 1, zval, stat=stat)
            if (stat == 0) then
               ssdi%ssdi_z(1) = zval
               ssdi%ssdi_z(2) = zval
            end if
         else if (n == 2) then
            do k = 1, 2
               call get_value(zarr, k, zval, stat=stat)
               if (stat == 0) ssdi%ssdi_z(k) = zval
            end do
         else
            call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                               'ssdi_z must be a 1- or 2-element array', &
                               'irrigation.ssdi.ssdi_z')
         end if
      end if

      ! [irrigation.ssdi.fixed]
      fixed_tbl => null()
      call get_value(ssdi_tbl, 'fixed', fixed_tbl, requested=.false., stat=stat)
      if (stat == 0 .and. associated(fixed_tbl)) then
         call get_optional_string_with_default(fixed_tbl, 'events_file', &
                                               ssdi%fixed%events_file, '', &
                                               'irrigation.ssdi.fixed.events_file', errors)
      end if

      ! [irrigation.ssdi.scheduled]
      sched_tbl => null()
      call get_value(ssdi_tbl, 'scheduled', sched_tbl, requested=.false., stat=stat)
      if (stat == 0 .and. associated(sched_tbl)) then
         call read_scheduled_fields(sched_tbl, ssdi%scheduled, errors)
      end if
   end subroutine read_irrigation_ssdi_toml


   subroutine read_scheduled_fields(sec, sched, errors)
      type(toml_table), pointer,         intent(in)    :: sec
      type(irrigation_ssdi_scheduled_t), intent(inout) :: sched
      type(error_collection_t),          intent(inout) :: errors

      call get_optional_int_with_default(sec, 'sched_type', sched%sched_type, 0, &
                                         'irrigation.ssdi.scheduled.sched_type', errors)
      call get_optional_real_with_default(sec, 'threshold', sched%threshold, 0.0_real64, &
                                          'irrigation.ssdi.scheduled.threshold', errors)
      call get_optional_real_with_default(sec, 'threshold_depth', sched%threshold_depth, 0.0_real64, &
                                          'irrigation.ssdi.scheduled.threshold_depth', errors)
      call get_optional_real_with_default(sec, 'ssdi_amount', sched%ssdi_amount, 0.0_real64, &
                                          'irrigation.ssdi.scheduled.ssdi_amount', errors)
      call get_optional_real_with_default(sec, 'ssdi_appl_rate', sched%ssdi_appl_rate, 0.0_real64, &
                                          'irrigation.ssdi.scheduled.ssdi_appl_rate', errors)
      call get_optional_int_with_default(sec, 'sw_interval', sched%sw_interval, 0, &
                                         'irrigation.ssdi.scheduled.sw_interval', errors)
      call get_optional_int_with_default(sec, 'days_interval', sched%days_interval, 1, &
                                         'irrigation.ssdi.scheduled.days_interval', errors)
   end subroutine read_scheduled_fields

end module read_irrigation_ssdi_toml_mod
