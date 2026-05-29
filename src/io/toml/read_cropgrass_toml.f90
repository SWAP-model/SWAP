!> Reader for a type-3 (grass / WOFOST grass) .crp.toml file.
!!
!! Same shape as read_cropfixed_toml plus grass-specific [mowing] and
!! [grazing] sections. Phase 4c-a covers core schema; tables and the
!! per-event mowing/grazing arrays are added during parity-test iteration.
!! Phase 4d Task 16 extended the [mowing] and [grazing] sections with
!! per-event fields (mowing_dates / mowing_heights / lsdb) and the
!! associated DM-threshold scalars. Phase 3 (Task 3) extends the parser
!! with all remaining sections: [crop_state], [green_area], [assimilation],
!! [crop_factor], [oxygen_stress] (with [oxygen_stress.bartholomeus]),
!! [drought_stress], [compensation], [management], and [co2]. Reader stays
!! mechanical: every new field is optional and absence leaves config at
!! defaults; the validator owns required-when-switch-on logic.
module read_cropgrass_toml_mod
   use iso_fortran_env, only: real64
   use tomlf, only: toml_table, toml_array, get_value, len
   use cropgrass_config_mod, only: cropgrass_config_t
   use toml_field_helpers_mod, only: get_table,                    &
                                     get_optional_int_with_default, &
                                     get_optional_real_with_default
   use toml_array_helpers_mod, only: read_real_array_1d
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
      type(toml_table), pointer :: cs, ga, assim, cf_sec, oxy, barto, dro, comp, mgmt, co2_sec

      ! [phenology] — pre-existing scalars + Phase 3 extensions
      call get_table(doc, 'phenology', ph, 'phenology', errors)
      if (associated(ph)) then
         call get_optional_int_with_default(ph,  'idev',  config%idev,  2, 'phenology.idev',  errors)
         call get_optional_int_with_default(ph,  'lcc',   config%lcc,   0, 'phenology.lcc',   errors)
         call get_optional_real_with_default(ph, 'tbase', config%tbase, 0.0_real64, 'phenology.tbase', errors)
         call get_optional_real_with_default(ph, 'tsum1', config%tsum1, 0.0_real64, 'phenology.tsum1', errors)
         call get_optional_real_with_default(ph, 'tsum2', config%tsum2, 0.0_real64, 'phenology.tsum2', errors)
         ! Phase 3 additions
         call get_optional_int_with_default(ph, 'swtsum',    config%swtsum,    1, 'phenology.swtsum',    errors)
         if (config%swtsum == 2) then
            call get_optional_real_with_default(ph, 'tsumtemp',  config%tsumtemp,  0.0_real64, 'phenology.tsumtemp',  errors)
            call get_optional_real_with_default(ph, 'tsumdepth', config%tsumdepth, 0.0_real64, 'phenology.tsumdepth', errors)
            call get_optional_int_with_default(ph,  'tsumtime',  config%tsumtime,  0,          'phenology.tsumtime',  errors)
         end if
      end if

      ! [light] — pre-existing
      call get_table(doc, 'light', light, 'light', errors)
      if (associated(light)) then
         call get_optional_real_with_default(light, 'kdif', config%kdif, 0.0_real64, 'light.kdif', errors)
         call get_optional_real_with_default(light, 'kdir', config%kdir, 0.0_real64, 'light.kdir', errors)
         call get_optional_real_with_default(light, 'eff',  config%eff,  0.0_real64, 'light.eff',  errors)
         call get_optional_real_with_default(light, 'amax', config%amax, 0.0_real64, 'light.amax', errors)
      end if

      ! [crop_state] — Phase 3 new section
      call get_table(doc, 'crop_state', cs, 'crop_state', errors)
      if (associated(cs)) then
         call get_optional_real_with_default(cs, 'tdwi',   config%tdwi,   1000.0_real64, 'crop_state.tdwi',   errors)
         call get_optional_real_with_default(cs, 'laiem',  config%laiem,  0.63_real64,   'crop_state.laiem',  errors)
         call get_optional_real_with_default(cs, 'rgrlai', config%rgrlai, 0.007_real64,  'crop_state.rgrlai', errors)
      end if

      ! [green_area] — Phase 3 new section
      call get_table(doc, 'green_area', ga, 'green_area', errors)
      if (associated(ga)) then
         call get_optional_real_with_default(ga, 'ssa',  config%ssa,  0.0_real64,  'green_area.ssa',  errors)
         call get_optional_real_with_default(ga, 'span', config%span, 30.0_real64, 'green_area.span', errors)
         call read_real_array_1d(ga, 'slatb', config%slatb, 'green_area.slatb', errors)
      end if

      ! [assimilation] — Phase 3 new section (light interception, biomass conversion,
      !                  respiration, partitioning, death rates)
      call get_table(doc, 'assimilation', assim, 'assimilation', errors)
      if (associated(assim)) then
         ! Light interception (grass TOML puts these here, not under [light])
         call get_optional_real_with_default(assim, 'kdif', config%kdif, 0.0_real64, 'assimilation.kdif', errors)
         call get_optional_real_with_default(assim, 'kdir', config%kdir, 0.0_real64, 'assimilation.kdir', errors)
         call get_optional_real_with_default(assim, 'eff',  config%eff,  0.0_real64, 'assimilation.eff',  errors)
         call get_optional_real_with_default(assim, 'cvl',   config%cvl,   0.685_real64,  'assimilation.cvl',   errors)
         call get_optional_real_with_default(assim, 'cvr',   config%cvr,   0.694_real64,  'assimilation.cvr',   errors)
         call get_optional_real_with_default(assim, 'cvs',   config%cvs,   0.662_real64,  'assimilation.cvs',   errors)
         call get_optional_real_with_default(assim, 'q10',   config%q10,   2.0_real64,    'assimilation.q10',   errors)
         call get_optional_real_with_default(assim, 'rml',   config%rml,   0.03_real64,   'assimilation.rml',   errors)
         call get_optional_real_with_default(assim, 'rmr',   config%rmr,   0.015_real64,  'assimilation.rmr',   errors)
         call get_optional_real_with_default(assim, 'rms',   config%rms,   0.015_real64,  'assimilation.rms',   errors)
         call get_optional_real_with_default(assim, 'perdl', config%perdl, 0.05_real64,   'assimilation.perdl', errors)
         call read_real_array_1d(assim, 'amaxtb',  config%amaxtb,  'assimilation.amaxtb',  errors)
         call read_real_array_1d(assim, 'tmpftb',  config%tmpftb,  'assimilation.tmpftb',  errors)
         call read_real_array_1d(assim, 'tmnftb',  config%tmnftb,  'assimilation.tmnftb',  errors)
         call read_real_array_1d(assim, 'rfsetb',  config%rfsetb,  'assimilation.rfsetb',  errors)
         call read_real_array_1d(assim, 'frtb',    config%frtb,    'assimilation.frtb',    errors)
         call read_real_array_1d(assim, 'fltb',    config%fltb,    'assimilation.fltb',    errors)
         call read_real_array_1d(assim, 'fstb',    config%fstb,    'assimilation.fstb',    errors)
         call read_real_array_1d(assim, 'rdrrtb',  config%rdrrtb,  'assimilation.rdrrtb',  errors)
         call read_real_array_1d(assim, 'rdrstb',  config%rdrstb,  'assimilation.rdrstb',  errors)
      end if

      ! [crop_factor] — Phase 3 new section
      call get_table(doc, 'crop_factor', cf_sec, 'crop_factor', errors)
      if (associated(cf_sec)) then
         call get_optional_int_with_default(cf_sec,  'swcf',   config%swcf,   2,           'crop_factor.swcf',   errors)
         call get_optional_real_with_default(cf_sec, 'albedo', config%albedo, 0.23_real64, 'crop_factor.albedo', errors)
         call get_optional_real_with_default(cf_sec, 'rsw',    config%rsw,    0.0_real64,  'crop_factor.rsw',    errors)
         call get_optional_int_with_default(cf_sec,  'swinter', config%swinter, 1, 'crop_factor.swinter', errors)
         call read_real_array_1d(cf_sec, 'cftb',  config%cftb,  'crop_factor.cftb',  errors)
         call read_real_array_1d(cf_sec, 'chtb',  config%chtb,  'crop_factor.chtb',  errors)
         call read_real_array_1d(cf_sec, 'rdctb', config%rdctb, 'crop_factor.rdctb', errors)
      end if

      ! [root] — pre-existing scalars + Phase 3 extensions
      call get_table(doc, 'root', root, 'root', errors)
      if (associated(root)) then
         call get_optional_real_with_default(root, 'rdi', config%rdi, 0.0_real64, 'root.rdi', errors)
         call get_optional_real_with_default(root, 'rri', config%rri, 0.0_real64, 'root.rri', errors)
         call get_optional_real_with_default(root, 'rdc', config%rdc, 0.0_real64, 'root.rdc', errors)
         ! Phase 3 additions
         call get_optional_int_with_default(root,  'swrd',     config%swrd,     2,             'root.swrd',     errors)
         call get_optional_int_with_default(root,  'swdmi2rd', config%swdmi2rd, 0,             'root.swdmi2rd', errors)
         call get_optional_int_with_default(root,  'swrdc',    config%swrdc,    0,             'root.swrdc',    errors)
         call get_optional_real_with_default(root, 'wrtmax',   config%wrtmax,   3000.0_real64, 'root.wrtmax',   errors)
         call read_real_array_1d(root, 'rdtb',  config%rdtb,  'root.rdtb',  errors)
         call read_real_array_1d(root, 'rlwtb', config%rlwtb, 'root.rlwtb', errors)
      end if

      ! [water_stress] — pre-existing
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

      ! [oxygen_stress] — Phase 3 new section with [oxygen_stress.bartholomeus] sub-table
      call get_table(doc, 'oxygen_stress', oxy, 'oxygen_stress', errors)
      if (associated(oxy)) then
         call get_optional_int_with_default(oxy,  'swoxygen',    config%swoxygen,    1,          'oxygen_stress.swoxygen',    errors)
         call get_optional_int_with_default(oxy,  'swwrtnonox',  config%swwrtnonox,  0,          'oxygen_stress.swwrtnonox',  errors)
         call get_optional_real_with_default(oxy, 'aeratecrit',  config%aeratecrit,  1.0e-4_real64, 'oxygen_stress.aeratecrit', errors)

         call get_table(oxy, 'bartholomeus', barto, 'oxygen_stress.bartholomeus', errors)
         if (associated(barto)) then
            call get_optional_int_with_default(barto,  'swoxygentype',            config%swoxygentype,           1,               'oxygen_stress.bartholomeus.swoxygentype',            errors)
            call get_optional_real_with_default(barto, 'q10_microbial',           config%q10_microbial,          2.8_real64,      'oxygen_stress.bartholomeus.q10_microbial',           errors)
            call get_optional_real_with_default(barto, 'specific_resp_humus',     config%specific_resp_humus,    1.6e-3_real64,   'oxygen_stress.bartholomeus.specific_resp_humus',     errors)
            call get_optional_real_with_default(barto, 'srl',                     config%srl,                    383571.0_real64, 'oxygen_stress.bartholomeus.srl',                     errors)
            call get_optional_int_with_default(barto,  'swrootradius',            config%swrootradius,           2,               'oxygen_stress.bartholomeus.swrootradius',            errors)
            call get_optional_real_with_default(barto, 'dry_mat_cont_roots',      config%dry_mat_cont_roots,     0.075_real64,    'oxygen_stress.bartholomeus.dry_mat_cont_roots',      errors)
            call get_optional_real_with_default(barto, 'air_filled_root_por',     config%air_filled_root_por,    0.05_real64,     'oxygen_stress.bartholomeus.air_filled_root_por',     errors)
            call get_optional_real_with_default(barto, 'spec_weight_root_tissue', config%spec_weight_root_tissue, 1.0e3_real64,   'oxygen_stress.bartholomeus.spec_weight_root_tissue', errors)
            call get_optional_real_with_default(barto, 'var_a',                   config%var_a,                  4.175e-10_real64,'oxygen_stress.bartholomeus.var_a',                   errors)
            call get_optional_real_with_default(barto, 'root_radiuso2',           config%root_radiusO2,          0.000075_real64, 'oxygen_stress.bartholomeus.root_radiuso2',           errors)
         end if
      end if

      ! [drought_stress] — Phase 3 new section
      call get_table(doc, 'drought_stress', dro, 'drought_stress', errors)
      if (associated(dro)) then
         call get_optional_int_with_default(dro, 'swdrought', config%swdrought, 1, 'drought_stress.swdrought', errors)
      end if

      ! [salinity] — pre-existing
      call get_table(doc, 'salinity', salt, 'salinity', errors)
      if (associated(salt)) then
         call get_optional_real_with_default(salt, 'ecmax',    config%ecmax,    0.0_real64, 'salt.ecmax',    errors)
         call get_optional_real_with_default(salt, 'ecslop',   config%ecslop,   0.0_real64, 'salt.ecslop',   errors)
         call get_optional_int_with_default(salt,  'swsalinity', config%swsalinity, 0,      'salt.swsalinity', errors)
      end if

      ! [compensation] — Phase 3 new section
      call get_table(doc, 'compensation', comp, 'compensation', errors)
      if (associated(comp)) then
         call get_optional_int_with_default(comp,  'swcompensate', config%swcompensate, 0,          'compensation.swcompensate', errors)
         call get_optional_int_with_default(comp,  'swstressor',   config%swstressor,   1,          'compensation.swstressor',   errors)
         call get_optional_real_with_default(comp, 'alphacrit',    config%alphacrit,    1.0_real64, 'compensation.alphacrit',    errors)
         call get_optional_real_with_default(comp, 'dcritrtz',     config%dcritrtz,     0.0_real64, 'compensation.dcritrtz',     errors)
      end if

      ! [interception] — pre-existing
      call get_table(doc, 'interception', inter, 'interception', errors)
      if (associated(inter)) then
         call get_optional_real_with_default(inter, 'cofab', config%cofab, 0.0_real64, 'inter.cofab', errors)
      end if

      ! [management] — Phase 3 new section
      call get_table(doc, 'management', mgmt, 'management', errors)
      if (associated(mgmt)) then
         call get_optional_real_with_default(mgmt, 'mowrest',    config%mowrest,    700.0_real64, 'management.mowrest',    errors)
         call get_optional_real_with_default(mgmt, 'dewrest',    config%dewrest,    850.0_real64, 'management.dewrest',    errors)
         call get_optional_int_with_default(mgmt,  'swpotrelmf', config%swpotrelmf, 1,            'management.swpotrelmf', errors)
         call get_optional_real_with_default(mgmt, 'relmf',      config%relmf,      1.0_real64,   'management.relmf',      errors)
         call read_array_1d_int(mgmt, 'seqgrazmow', config%seqgrazmow, config%nseqgrazmow, &
                                'management.seqgrazmow', errors)
         call read_real_array_1d(mgmt, 'dmmowtb',     config%dmmowtb,     'management.dmmowtb',     errors)
         call read_real_array_1d(mgmt, 'dmmowdelay',  config%dmmowdelay,  'management.dmmowdelay',  errors)
         call read_real_array_1d(mgmt, 'dmgrztb',     config%dmgrztb,     'management.dmgrztb',     errors)
         call read_real_array_1d(mgmt, 'lsda',        config%lsda,        'management.lsda',        errors)
         call read_real_array_1d(mgmt, 'daysgrazing', config%daysgrazing, 'management.daysgrazing', errors)
         call read_real_array_1d(mgmt, 'uptgrazing',  config%uptgrazing,  'management.uptgrazing',  errors)
         call read_real_array_1d(mgmt, 'lossgrazing', config%lossgrazing, 'management.lossgrazing', errors)
         call get_optional_int_with_default(mgmt, 'swlossmow', config%swlossmow, 0, 'management.swlossmow', errors)
         call get_optional_int_with_default(mgmt, 'swlossgrz', config%swlossgrz, 0, 'management.swlossgrz', errors)
      end if

      ! [co2] — Phase 3 new section
      call get_table(doc, 'co2', co2_sec, 'co2', errors)
      if (associated(co2_sec)) then
         call get_optional_int_with_default(co2_sec, 'swco2', config%swco2, 0, 'co2.swco2', errors)
      end if

      ! [mowing] — pre-existing
      call get_table(doc, 'mowing', mow, 'mowing', errors)
      if (associated(mow)) then
         call get_optional_int_with_default(mow, 'swharv',          config%swharv,         0,          'mowing.swharv',         errors)
         call get_optional_int_with_default(mow, 'nmow',            config%nmow,           0,          'mowing.nmow',           errors)
         call get_optional_int_with_default(mow, 'swdmmow',         config%swdmmow,        0,          'mowing.swdmmow',        errors)
         call get_optional_int_with_default(mow, 'maxdaymow',       config%maxdaymow,      0,          'mowing.maxdaymow',      errors)
         call get_optional_real_with_default(mow, 'dmharvest',      config%dmharvest,      0.0_real64, 'mowing.dmharvest',      errors)
         call get_optional_real_with_default(mow, 'daylastharvest', config%daylastharvest, 0.0_real64, 'mowing.daylastharvest', errors)
         call get_optional_real_with_default(mow, 'dmlastharvest',  config%dmlastharvest,  0.0_real64, 'mowing.dmlastharvest',  errors)
         call read_real_array_1d(mow, 'mowing_dates',   config%mowing_dates,   'mowing.mowing_dates',   errors)
         call read_real_array_1d(mow, 'mowing_heights', config%mowing_heights, 'mowing.mowing_heights', errors)
      end if

      ! [grazing] — pre-existing
      call get_table(doc, 'grazing', graz, 'grazing', errors)
      if (associated(graz)) then
         call get_optional_int_with_default(graz, 'swgraz',      config%swgraz,      0,          'grazing.swgraz',      errors)
         call get_optional_int_with_default(graz, 'nstart_graz', config%nstart_graz, 0,          'grazing.nstart_graz', errors)
         call get_optional_int_with_default(graz, 'nstop_graz',  config%nstop_graz,  0,          'grazing.nstop_graz',  errors)
         call get_optional_int_with_default(graz, 'maxdaygrz',   config%maxdaygrz,   0,          'grazing.maxdaygrz',   errors)
         call get_optional_int_with_default(graz, 'swdmgrz',     config%swdmgrz,     0,          'grazing.swdmgrz',     errors)
         call get_optional_real_with_default(graz, 'dmgrazing',  config%dmgrazing,   0.0_real64, 'grazing.dmgrazing',   errors)
         call get_optional_real_with_default(graz, 'tagprest',   config%tagprest,    0.0_real64, 'grazing.tagprest',    errors)
         call read_real_array_1d(graz, 'lsdb', config%lsdb, 'grazing.lsdb', errors)
      end if

      ! [irrigation_schedule] — pre-existing
      call get_table(doc, 'irrigation_schedule', irr_sched, 'irrigation_schedule', errors)
      call read_irrigation_schedule_from_section(irr_sched, config%schedule, errors)
   end subroutine read_cropgrass_toml

   !> Decode a flat TOML integer array at sec[key] into a 1-D integer
   !! allocatable. Also updates the companion count scalar. Absent key
   !! leaves arr unallocated and count unchanged. Empty array yields a
   !! 0-element allocation. Non-integer cells append a parse-type-mismatch
   !! error and leave arr unallocated.
   subroutine read_array_1d_int(sec, key, arr, count, context, errors)
      type(toml_table), pointer, intent(in)    :: sec
      character(len=*),          intent(in)    :: key
      integer, allocatable,      intent(out)   :: arr(:)
      integer,                   intent(inout) :: count
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
         count = 0
         return
      end if

      allocate(arr(n))
      arr = 0
      count = n

      do i = 1, n
         call get_value(outer, i, val, stat=stat)
         if (stat /= 0) then
            call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                               "non-int cell", context)
            if (allocated(arr)) deallocate(arr)
            count = 0
            return
         end if
         arr(i) = val
      end do
   end subroutine read_array_1d_int

end module read_cropgrass_toml_mod
