!> Reader for a type-3 (grass / WOFOST grass) .crp.toml file.
!!
!! Same shape as read_cropfixed_toml plus grass-specific [mowing] and
!! [grazing] sections. Phase 4c-a covers core schema; tables and the
!! per-event mowing/grazing arrays are added during parity-test iteration.
!! Phase 4d Task 16 extended the [mowing] and [grazing] sections with
!! per-event fields (mowing_dates / mowing_heights / lsdb) and the
!! associated DM-threshold scalars. Reader stays mechanical: every new
!! field is optional and absence leaves config at defaults; the validator
!! (Task 15) owns required-when-switch-on logic.
module read_cropgrass_toml_mod
   use iso_fortran_env, only: real64
   use tomlf, only: toml_table, toml_array, get_value, len
   use cropgrass_config_mod, only: cropgrass_config_t
   use toml_field_helpers_mod, only: get_table,                    &
                                     get_optional_int_with_default, &
                                     get_optional_real_with_default
   use read_irrigation_toml_mod, only: read_irrigation_schedule_from_section
   use error_mod, only: error_collection_t, ERR_PARSE_TYPE_MISMATCH
   implicit none
   private

   public :: read_cropgrass_toml

contains

   subroutine read_cropgrass_toml(doc, config, errors)
      type(toml_table), pointer,  intent(in)    :: doc
      type(cropgrass_config_t),   intent(inout) :: config
      type(error_collection_t),   intent(inout) :: errors

      type(toml_table), pointer :: ph, light, root, ws, salt, inter, mow, graz, irr_sched

      call get_table(doc, 'phenology', ph, 'phenology', errors)
      if (associated(ph)) then
         call get_optional_int_with_default(ph,  'idev',  config%idev,  2, 'phenology.idev',  errors)
         call get_optional_int_with_default(ph,  'lcc',   config%lcc,   0, 'phenology.lcc',   errors)
         call get_optional_real_with_default(ph, 'tbase', config%tbase, 0.0_real64, 'phenology.tbase', errors)
         call get_optional_real_with_default(ph, 'tsum1', config%tsum1, 0.0_real64, 'phenology.tsum1', errors)
         call get_optional_real_with_default(ph, 'tsum2', config%tsum2, 0.0_real64, 'phenology.tsum2', errors)
      end if

      call get_table(doc, 'light', light, 'light', errors)
      if (associated(light)) then
         call get_optional_real_with_default(light, 'kdif', config%kdif, 0.0_real64, 'light.kdif', errors)
         call get_optional_real_with_default(light, 'kdir', config%kdir, 0.0_real64, 'light.kdir', errors)
         call get_optional_real_with_default(light, 'eff',  config%eff,  0.0_real64, 'light.eff',  errors)
         call get_optional_real_with_default(light, 'amax', config%amax, 0.0_real64, 'light.amax', errors)
      end if

      call get_table(doc, 'root', root, 'root', errors)
      if (associated(root)) then
         call get_optional_real_with_default(root, 'rdi', config%rdi, 0.0_real64, 'root.rdi', errors)
         call get_optional_real_with_default(root, 'rri', config%rri, 0.0_real64, 'root.rri', errors)
         call get_optional_real_with_default(root, 'rdc', config%rdc, 0.0_real64, 'root.rdc', errors)
      end if

      call get_table(doc, 'water_stress', ws, 'water_stress', errors)
      if (associated(ws)) then
         call get_optional_real_with_default(ws, 'hlim1',  config%hlim1,  0.0_real64, 'ws.hlim1',  errors)
         call get_optional_real_with_default(ws, 'hlim2u', config%hlim2u, 0.0_real64, 'ws.hlim2u', errors)
         call get_optional_real_with_default(ws, 'hlim2l', config%hlim2l, 0.0_real64, 'ws.hlim2l', errors)
         call get_optional_real_with_default(ws, 'hlim3h', config%hlim3h, 0.0_real64, 'ws.hlim3h', errors)
         call get_optional_real_with_default(ws, 'hlim3l', config%hlim3l, 0.0_real64, 'ws.hlim3l', errors)
         call get_optional_real_with_default(ws, 'hlim4',  config%hlim4,  0.0_real64, 'ws.hlim4',  errors)
         call get_optional_real_with_default(ws, 'adcrh',  config%adcrh,  0.0_real64, 'ws.adcrh',  errors)
         call get_optional_real_with_default(ws, 'adcrl',  config%adcrl,  0.0_real64, 'ws.adcrl',  errors)
         call get_optional_real_with_default(ws, 'rsc',    config%rsc,    0.0_real64, 'ws.rsc',    errors)
      end if

      call get_table(doc, 'salinity', salt, 'salinity', errors)
      if (associated(salt)) then
         call get_optional_real_with_default(salt, 'ecmax',  config%ecmax,  0.0_real64, 'salt.ecmax',  errors)
         call get_optional_real_with_default(salt, 'ecslop', config%ecslop, 0.0_real64, 'salt.ecslop', errors)
      end if

      call get_table(doc, 'interception', inter, 'interception', errors)
      if (associated(inter)) then
         call get_optional_real_with_default(inter, 'cofab', config%cofab, 0.0_real64, 'inter.cofab', errors)
      end if

      call get_table(doc, 'mowing', mow, 'mowing', errors)
      if (associated(mow)) then
         call get_optional_int_with_default(mow, 'swharv',         config%swharv,         0,          'mowing.swharv',         errors)
         call get_optional_int_with_default(mow, 'nmow',           config%nmow,           0,          'mowing.nmow',           errors)
         call get_optional_int_with_default(mow, 'swdmmow',        config%swdmmow,        0,          'mowing.swdmmow',        errors)
         call get_optional_int_with_default(mow, 'maxdaymow',      config%maxdaymow,      0,          'mowing.maxdaymow',      errors)
         call get_optional_real_with_default(mow, 'dmharvest',     config%dmharvest,      0.0_real64, 'mowing.dmharvest',      errors)
         call get_optional_real_with_default(mow, 'daylastharvest', config%daylastharvest, 0.0_real64, 'mowing.daylastharvest', errors)
         call get_optional_real_with_default(mow, 'dmlastharvest', config%dmlastharvest,  0.0_real64, 'mowing.dmlastharvest',  errors)
         call read_array_1d(mow, 'mowing_dates',   config%mowing_dates,   'mowing.mowing_dates',   errors)
         call read_array_1d(mow, 'mowing_heights', config%mowing_heights, 'mowing.mowing_heights', errors)
      end if

      call get_table(doc, 'grazing', graz, 'grazing', errors)
      if (associated(graz)) then
         call get_optional_int_with_default(graz, 'swgraz',      config%swgraz,      0,          'grazing.swgraz',      errors)
         call get_optional_int_with_default(graz, 'nstart_graz', config%nstart_graz, 0,          'grazing.nstart_graz', errors)
         call get_optional_int_with_default(graz, 'nstop_graz',  config%nstop_graz,  0,          'grazing.nstop_graz',  errors)
         call get_optional_int_with_default(graz, 'maxdaygrz',   config%maxdaygrz,   0,          'grazing.maxdaygrz',   errors)
         call get_optional_int_with_default(graz, 'swdmgrz',     config%swdmgrz,     0,          'grazing.swdmgrz',     errors)
         call get_optional_real_with_default(graz, 'dmgrazing',  config%dmgrazing,   0.0_real64, 'grazing.dmgrazing',   errors)
         call get_optional_real_with_default(graz, 'tagprest',   config%tagprest,    0.0_real64, 'grazing.tagprest',    errors)
         call read_array_1d(graz, 'lsdb', config%lsdb, 'grazing.lsdb', errors)
      end if

      call get_table(doc, 'irrigation_schedule', irr_sched, 'irrigation_schedule', errors)
      call read_irrigation_schedule_from_section(irr_sched, config%schedule, errors)
   end subroutine read_cropgrass_toml

   !> Decode a flat TOML array at sec[key] into a 1-D real(real64)
   !! allocatable. Absent key leaves arr unallocated. Empty array
   !! (`key = []`) yields a 0-element allocation. Non-real cells append
   !! a parse-type-mismatch error and leave arr unallocated.
   !!
   !! Local copy of the helper from read_heat_toml — that helper is
   !! private to its module, so duplicating here keeps the readers
   !! decoupled (Phase 4d Task 16).
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

end module read_cropgrass_toml_mod
