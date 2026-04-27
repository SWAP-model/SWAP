!> Reader for a type-3 (grass / WOFOST grass) .crp.toml file.
!!
!! Same shape as read_cropfixed_toml plus grass-specific [mowing] and
!! [grazing] sections. Phase 4c-a covers core schema; tables and the
!! per-event mowing/grazing arrays are added during parity-test iteration.
module read_cropgrass_toml_mod
   use iso_fortran_env, only: real64
   use tomlf, only: toml_table
   use cropgrass_config_mod, only: cropgrass_config_t
   use toml_field_helpers_mod, only: get_table,                    &
                                     get_optional_int_with_default, &
                                     get_optional_real_with_default
   use read_irrigation_toml_mod, only: read_irrigation_schedule_from_section
   use error_mod, only: error_collection_t
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
         call get_optional_int_with_default(mow, 'swharv', config%swharv, 0, 'mowing.swharv', errors)
         call get_optional_int_with_default(mow, 'nmow',   config%nmow,   0, 'mowing.nmow',   errors)
      end if

      call get_table(doc, 'grazing', graz, 'grazing', errors)
      if (associated(graz)) then
         call get_optional_int_with_default(graz, 'swgraz',      config%swgraz,      0,         'grazing.swgraz',      errors)
         call get_optional_int_with_default(graz, 'nstart_graz', config%nstart_graz, 0, 'grazing.nstart_graz', errors)
         call get_optional_int_with_default(graz, 'nstop_graz',  config%nstop_graz,  0, 'grazing.nstop_graz',  errors)
      end if

      call get_table(doc, 'irrigation_schedule', irr_sched, 'irrigation_schedule', errors)
      call read_irrigation_schedule_from_section(irr_sched, config%schedule, errors)
   end subroutine read_cropgrass_toml

end module read_cropgrass_toml_mod
