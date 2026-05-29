!> Reader for a type-1 (fixed crop) .crp.toml file.
!!
!! Given the ROOT table of a loaded .crp.toml document, populates a
!! cropfixed_config_t. Phase 1 of the .crp port (Phase 4f) extends this
!! to a 1:1 match with legacy readcropfixed: ~30 new scalars + 6 flat-pair
!! tables (gctb, cftb, chtb, rdtb, rdctb). cfeictb is declared on the
!! schema but not read here — it is only relevant when swcf=3, which
!! is stub-errored at validate time per ADR 0015.
module read_cropfixed_toml_mod
   use iso_fortran_env, only: real64
   use tomlf, only: toml_table
   use cropfixed_config_mod, only: cropfixed_config_t
   use toml_field_helpers_mod, only: get_table,                    &
                                     get_optional_int_with_default, &
                                     get_optional_real_with_default
   use toml_array_helpers_mod, only: read_real_pair_array
   use read_irrigation_toml_mod, only: read_irrigation_schedule_from_section
   use error_mod, only: error_collection_t
   implicit none
   private

   public :: read_cropfixed_toml

contains

   subroutine read_cropfixed_toml(doc, config, errors)
      type(toml_table), pointer,  intent(in)    :: doc
      type(cropfixed_config_t),   intent(inout) :: config
      type(error_collection_t),   intent(inout) :: errors

      type(toml_table), pointer :: prep, harv, ph, light, lai, cf, root, &
                                    ox, dr, sa, cp, inter, sched, irr_sched

      ! Preparation, sowing, germination
      call get_table(doc, 'preparation', prep, 'preparation', errors)
      if (associated(prep)) then
         call get_optional_int_with_default(prep, 'swprep', config%swprep, 0, 'preparation.swprep', errors)
         call get_optional_int_with_default(prep, 'swsow',  config%swsow,  0, 'preparation.swsow',  errors)
         call get_optional_int_with_default(prep, 'swgerm', config%swgerm, 0, 'preparation.swgerm', errors)
      end if

      ! Harvest
      call get_table(doc, 'harvest', harv, 'harvest', errors)
      if (associated(harv)) then
         call get_optional_real_with_default(harv, 'dvsend', config%dvsend, 2.0_real64, 'harvest.dvsend', errors)
         call get_optional_int_with_default (harv, 'swharv', config%swharv, 0,           'harvest.swharv', errors)
      end if

      ! Phenology
      call get_table(doc, 'phenology', ph, 'phenology', errors)
      if (associated(ph)) then
         call get_optional_int_with_default (ph, 'idev',   config%idev,   1, 'phenology.idev',   errors)
         call get_optional_int_with_default (ph, 'lcc',    config%lcc,    0, 'phenology.lcc',    errors)
         call get_optional_real_with_default(ph, 'tsumea', config%tsumea, 0.0_real64, 'phenology.tsumea', errors)
         call get_optional_real_with_default(ph, 'tsumam', config%tsumam, 0.0_real64, 'phenology.tsumam', errors)
         call get_optional_real_with_default(ph, 'tbase',  config%tbase,  0.0_real64, 'phenology.tbase',  errors)
      end if

      ! Light
      call get_table(doc, 'light', light, 'light', errors)
      if (associated(light)) then
         call get_optional_real_with_default(light, 'kdif', config%kdif, 0.0_real64, 'light.kdif', errors)
         call get_optional_real_with_default(light, 'kdir', config%kdir, 0.0_real64, 'light.kdir', errors)
      end if

      ! LAI / SCF table
      call get_table(doc, 'lai', lai, 'lai', errors)
      if (associated(lai)) then
         call get_optional_int_with_default(lai, 'swgc', config%swgc, 1, 'lai.swgc', errors)
         call read_real_pair_array(lai, 'gctb', config%gctb, 'lai.gctb', errors)
      end if

      ! Crop factor / height
      call get_table(doc, 'crop_factor', cf, 'crop_factor', errors)
      if (associated(cf)) then
         call get_optional_int_with_default(cf, 'swcf', config%swcf, 1, 'crop_factor.swcf', errors)
         call read_real_pair_array(cf, 'cftb', config%cftb, 'crop_factor.cftb', errors)
         call read_real_pair_array(cf, 'chtb', config%chtb, 'crop_factor.chtb', errors)
         call get_optional_real_with_default(cf, 'albedo', config%albedo, 0.23_real64, 'crop_factor.albedo', errors)
         call get_optional_real_with_default(cf, 'rsc',    config%rsc,    0.0_real64,  'crop_factor.rsc',    errors)
         call get_optional_real_with_default(cf, 'rsw',    config%rsw,    0.0_real64,  'crop_factor.rsw',    errors)
      end if

      ! Root
      call get_table(doc, 'root', root, 'root', errors)
      if (associated(root)) then
         call get_optional_int_with_default (root, 'swrd',     config%swrd,     1,           'root.swrd',     errors)
         call get_optional_int_with_default (root, 'swdmi2rd', config%swdmi2rd, 0,           'root.swdmi2rd', errors)
         call get_optional_int_with_default (root, 'swrdc',    config%swrdc,    0,           'root.swrdc',    errors)
         call read_real_pair_array(root, 'rdtb',  config%rdtb,  'root.rdtb',  errors)
         call read_real_pair_array(root, 'rdctb', config%rdctb, 'root.rdctb', errors)
         call get_optional_real_with_default(root, 'rdi', config%rdi, 0.0_real64, 'root.rdi', errors)
         call get_optional_real_with_default(root, 'rri', config%rri, 0.0_real64, 'root.rri', errors)
         call get_optional_real_with_default(root, 'rdc', config%rdc, 0.0_real64, 'root.rdc', errors)
      end if

      ! Oxygen stress
      call get_table(doc, 'oxygen_stress', ox, 'oxygen_stress', errors)
      if (associated(ox)) then
         call get_optional_int_with_default (ox, 'swoxygen',   config%swoxygen,   0, 'oxygen_stress.swoxygen',   errors)
         call get_optional_int_with_default (ox, 'swwrtnonox', config%swwrtnonox, 0, 'oxygen_stress.swwrtnonox', errors)
         call get_optional_real_with_default(ox, 'aeratecrit', config%aeratecrit, 1.0e-4_real64, 'oxygen_stress.aeratecrit', errors)
         call get_optional_real_with_default(ox, 'max_resp_factor', config%max_resp_factor, 1.0_real64, 'oxygen_stress.max_resp_factor', errors)
         call get_optional_real_with_default(ox, 'hlim1',      config%hlim1,  0.0_real64, 'oxygen_stress.hlim1',  errors)
         call get_optional_real_with_default(ox, 'hlim2u',     config%hlim2u, 0.0_real64, 'oxygen_stress.hlim2u', errors)
         call get_optional_real_with_default(ox, 'hlim2l',     config%hlim2l, 0.0_real64, 'oxygen_stress.hlim2l', errors)
      end if

      ! Drought stress
      call get_table(doc, 'drought_stress', dr, 'drought_stress', errors)
      if (associated(dr)) then
         call get_optional_int_with_default (dr, 'swdrought', config%swdrought, 1,           'drought_stress.swdrought', errors)
         call get_optional_real_with_default(dr, 'hlim3h',    config%hlim3h, 0.0_real64, 'drought_stress.hlim3h', errors)
         call get_optional_real_with_default(dr, 'hlim3l',    config%hlim3l, 0.0_real64, 'drought_stress.hlim3l', errors)
         call get_optional_real_with_default(dr, 'hlim4',     config%hlim4,  0.0_real64, 'drought_stress.hlim4',  errors)
         call get_optional_real_with_default(dr, 'adcrh',     config%adcrh,  0.0_real64, 'drought_stress.adcrh',  errors)
         call get_optional_real_with_default(dr, 'adcrl',     config%adcrl,  0.0_real64, 'drought_stress.adcrl',  errors)
      end if

      ! Salinity stress (subordinate fields read but stub-errored at validate)
      call get_table(doc, 'salinity_stress', sa, 'salinity_stress', errors)
      if (associated(sa)) then
         call get_optional_int_with_default (sa, 'swsalinity', config%swsalinity, 0,           'salinity_stress.swsalinity', errors)
         call get_optional_real_with_default(sa, 'saltmax',    config%saltmax,    0.0_real64,  'salinity_stress.saltmax',    errors)
         call get_optional_real_with_default(sa, 'saltslope',  config%saltslope,  0.0_real64,  'salinity_stress.saltslope',  errors)
         call get_optional_real_with_default(sa, 'salthead',   config%salthead,   0.0_real64,  'salinity_stress.salthead',   errors)
         call get_optional_real_with_default(sa, 'ecmax',      config%ecmax,      0.0_real64,  'salinity_stress.ecmax',      errors)
         call get_optional_real_with_default(sa, 'ecslop',     config%ecslop,     0.0_real64,  'salinity_stress.ecslop',     errors)
      end if

      ! Compensation (stub-errored at validate)
      call get_table(doc, 'compensation', cp, 'compensation', errors)
      if (associated(cp)) then
         call get_optional_int_with_default (cp, 'swcompensate', config%swcompensate, 0, 'compensation.swcompensate', errors)
         call get_optional_int_with_default (cp, 'swstressor',   config%swstressor,   1, 'compensation.swstressor',   errors)
         call get_optional_real_with_default(cp, 'alphacrit',    config%alphacrit,    1.0_real64,  'compensation.alphacrit', errors)
         call get_optional_real_with_default(cp, 'dcritrtz',     config%dcritrtz,     0.0_real64,  'compensation.dcritrtz',  errors)
      end if

      ! Interception
      call get_table(doc, 'interception', inter, 'interception', errors)
      if (associated(inter)) then
         call get_optional_int_with_default (inter, 'swinter', config%swinter, 1,           'interception.swinter', errors)
         call get_optional_real_with_default(inter, 'cofab',   config%cofab,   0.0_real64,  'interception.cofab',   errors)
      end if

      ! Scheduling top-level switch
      call get_table(doc, 'scheduling', sched, 'scheduling', errors)
      if (associated(sched)) then
         call get_optional_int_with_default(sched, 'schedule', config%schedule_switch, 0, 'scheduling.schedule', errors)
      end if

      ! Per-crop irrigation schedule (existing field - kept)
      call get_table(doc, 'irrigation_schedule', irr_sched, 'irrigation_schedule', errors)
      call read_irrigation_schedule_from_section(irr_sched, config%schedule, errors)
   end subroutine read_cropfixed_toml

end module read_cropfixed_toml_mod
