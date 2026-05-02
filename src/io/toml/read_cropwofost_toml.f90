!> Reader for a type-2 (WOFOST) .crp.toml file.
!!
!! Populates a cropwofost_config_t from a fully-formed .crp.toml document.
!! Tables are encoded as TOML arrays of arrays (e.g., dtsmtb = [[x,y],...])
!! and decoded into 2-D allocatable real(real64) arrays. The reader is
!! mechanical: scalars read one-per-line, tables one-block-per-table. No
!! switch-conditional skipping - validators flag inconsistent state.
module read_cropwofost_toml_mod
   use iso_fortran_env, only: real64
   use tomlf, only: toml_table, toml_array, get_value, len
   use cropwofost_config_mod, only: cropwofost_config_t
   use toml_field_helpers_mod, only: get_table,                         &
                                     get_optional_int_with_default,     &
                                     get_optional_real_with_default,    &
                                     get_optional_string_with_default
   use read_irrigation_toml_mod, only: read_irrigation_schedule_from_section
   use error_mod, only: error_collection_t, ERR_PARSE_TYPE_MISMATCH
   implicit none
   private

   public :: read_cropwofost_toml

contains

   subroutine read_cropwofost_toml(doc_root, config, errors)
      type(toml_table), pointer,  intent(in)    :: doc_root
      type(cropwofost_config_t),  intent(inout) :: config
      type(error_collection_t),   intent(inout) :: errors

      type(toml_table), pointer :: prep, sow, germ, harv, cf, ph, init, ga,  &
                                   asm, conv, resp, part, deth, root, oxy,    &
                                   drou, salt, comp, inter, co2t, mgmt,       &
                                   irr_sched, soy, bul, nut

      call get_table(doc_root, 'preparation', prep, 'preparation', errors)
      if (associated(prep)) then
         call get_optional_int_with_default(prep,  'swprep',       config%preparation%swprep,       0,           'preparation.swprep',       errors)
         call get_optional_real_with_default(prep, 'zprep',        config%preparation%zprep,        0.0_real64,  'preparation.zprep',        errors)
         call get_optional_real_with_default(prep, 'hprep',        config%preparation%hprep,        0.0_real64,  'preparation.hprep',        errors)
         call get_optional_int_with_default(prep,  'maxprepdelay', config%preparation%maxprepdelay, 0,           'preparation.maxprepdelay', errors)
      end if

      call get_table(doc_root, 'sowing', sow, 'sowing', errors)
      if (associated(sow)) then
         call get_optional_int_with_default(sow,  'swsow',       config%sowing%swsow,       0,          'sowing.swsow',       errors)
         call get_optional_real_with_default(sow, 'zsow',        config%sowing%zsow,        0.0_real64, 'sowing.zsow',        errors)
         call get_optional_real_with_default(sow, 'hsow',        config%sowing%hsow,        0.0_real64, 'sowing.hsow',        errors)
         call get_optional_real_with_default(sow, 'ztempsow',    config%sowing%ztempsow,    0.0_real64, 'sowing.ztempsow',    errors)
         call get_optional_real_with_default(sow, 'tempsow',     config%sowing%tempsow,     0.0_real64, 'sowing.tempsow',     errors)
         call get_optional_int_with_default(sow,  'maxsowdelay', config%sowing%maxsowdelay, 0,          'sowing.maxsowdelay', errors)
      end if

      call get_table(doc_root, 'germination', germ, 'germination', errors)
      if (associated(germ)) then
         call get_optional_int_with_default(germ,  'swgerm',     config%germination%swgerm,     0,          'germination.swgerm',     errors)
         call get_optional_real_with_default(germ, 'tsumemeopt', config%germination%tsumemeopt, 0.0_real64, 'germination.tsumemeopt', errors)
         call get_optional_real_with_default(germ, 'tbasem',     config%germination%tbasem,     0.0_real64, 'germination.tbasem',     errors)
         call get_optional_real_with_default(germ, 'teffmx',     config%germination%teffmx,     0.0_real64, 'germination.teffmx',     errors)
         call get_optional_real_with_default(germ, 'hdrygerm',   config%germination%hdrygerm,   0.0_real64, 'germination.hdrygerm',   errors)
         call get_optional_real_with_default(germ, 'hwetgerm',   config%germination%hwetgerm,   0.0_real64, 'germination.hwetgerm',   errors)
         call get_optional_real_with_default(germ, 'zgerm',      config%germination%zgerm,      0.0_real64, 'germination.zgerm',      errors)
         call get_optional_real_with_default(germ, 'agerm',      config%germination%agerm,      0.0_real64, 'germination.agerm',      errors)
      end if

      call get_table(doc_root, 'harvest', harv, 'harvest', errors)
      if (associated(harv)) then
         call get_optional_real_with_default(harv, 'dvsend', config%harvest%dvsend, 0.0_real64, 'harvest.dvsend', errors)
         call get_optional_int_with_default(harv,  'swharv', config%harvest%swharv, 0,          'harvest.swharv', errors)
      end if

      call get_table(doc_root, 'crop_factor', cf, 'crop_factor', errors)
      if (associated(cf)) then
         call get_optional_int_with_default(cf,  'swcf',   config%crop_factor%swcf,   0,          'crop_factor.swcf',   errors)
         call get_optional_real_with_default(cf, 'albedo', config%crop_factor%albedo, 0.0_real64, 'crop_factor.albedo', errors)
         call get_optional_real_with_default(cf, 'rsc',    config%crop_factor%rsc,    0.0_real64, 'crop_factor.rsc',    errors)
         call get_optional_real_with_default(cf, 'rsw',    config%crop_factor%rsw,    0.0_real64, 'crop_factor.rsw',    errors)
         call read_table_2d(cf, 'cftb', config%crop_factor%cftb, 2, 'crop_factor.cftb', errors)
         call read_table_2d(cf, 'chtb', config%crop_factor%chtb, 2, 'crop_factor.chtb', errors)
      end if

      call get_table(doc_root, 'phenology', ph, 'phenology', errors)
      if (associated(ph)) then
         call get_optional_int_with_default(ph,  'idsl',     config%phenology%idsl,     0,          'phenology.idsl',     errors)
         call get_optional_real_with_default(ph, 'tsumea',   config%phenology%tsumea,   0.0_real64, 'phenology.tsumea',   errors)
         call get_optional_real_with_default(ph, 'tsumam',   config%phenology%tsumam,   0.0_real64, 'phenology.tsumam',   errors)
         call get_optional_real_with_default(ph, 'dlo',      config%phenology%dlo,      0.0_real64, 'phenology.dlo',      errors)
         call get_optional_real_with_default(ph, 'dlc',      config%phenology%dlc,      0.0_real64, 'phenology.dlc',      errors)
         call get_optional_real_with_default(ph, 'vernsat',  config%phenology%vernsat,  0.0_real64, 'phenology.vernsat',  errors)
         call get_optional_real_with_default(ph, 'vernbase', config%phenology%vernbase, 0.0_real64, 'phenology.vernbase', errors)
         call get_optional_real_with_default(ph, 'verndvs',  config%phenology%verndvs,  0.0_real64, 'phenology.verndvs',  errors)
         call read_table_2d(ph, 'dtsmtb', config%phenology%dtsmtb, 2, 'phenology.dtsmtb', errors)
         call read_table_2d(ph, 'verntb', config%phenology%verntb, 2, 'phenology.verntb', errors)
      end if

      call get_table(doc_root, 'initial', init, 'initial', errors)
      if (associated(init)) then
         call get_optional_real_with_default(init, 'tdwi',   config%initial%tdwi,   0.0_real64, 'initial.tdwi',   errors)
         call get_optional_real_with_default(init, 'laiem',  config%initial%laiem,  0.0_real64, 'initial.laiem',  errors)
         call get_optional_real_with_default(init, 'rgrlai', config%initial%rgrlai, 0.0_real64, 'initial.rgrlai', errors)
      end if

      call get_table(doc_root, 'green_area', ga, 'green_area', errors)
      if (associated(ga)) then
         call get_optional_real_with_default(ga, 'spa',   config%green_area%spa,   0.0_real64, 'green_area.spa',   errors)
         call get_optional_real_with_default(ga, 'ssa',   config%green_area%ssa,   0.0_real64, 'green_area.ssa',   errors)
         call get_optional_real_with_default(ga, 'span',  config%green_area%span,  0.0_real64, 'green_area.span',  errors)
         call get_optional_real_with_default(ga, 'tbase', config%green_area%tbase, 0.0_real64, 'green_area.tbase', errors)
         call read_table_2d(ga, 'slatb', config%green_area%slatb, 2, 'green_area.slatb', errors)
      end if

      call get_table(doc_root, 'assimilation', asm, 'assimilation', errors)
      if (associated(asm)) then
         call get_optional_real_with_default(asm, 'kdif', config%assimilation%kdif, 0.0_real64, 'assimilation.kdif', errors)
         call get_optional_real_with_default(asm, 'kdir', config%assimilation%kdir, 0.0_real64, 'assimilation.kdir', errors)
         call get_optional_real_with_default(asm, 'eff',  config%assimilation%eff,  0.0_real64, 'assimilation.eff',  errors)
         call read_table_2d(asm, 'amaxtb', config%assimilation%amaxtb, 2, 'assimilation.amaxtb', errors)
         call read_table_2d(asm, 'tmpftb', config%assimilation%tmpftb, 2, 'assimilation.tmpftb', errors)
         call read_table_2d(asm, 'tmnftb', config%assimilation%tmnftb, 2, 'assimilation.tmnftb', errors)
      end if

      call get_table(doc_root, 'conversion', conv, 'conversion', errors)
      if (associated(conv)) then
         call get_optional_real_with_default(conv, 'cvl', config%conversion%cvl, 0.0_real64, 'conversion.cvl', errors)
         call get_optional_real_with_default(conv, 'cvo', config%conversion%cvo, 0.0_real64, 'conversion.cvo', errors)
         call get_optional_real_with_default(conv, 'cvr', config%conversion%cvr, 0.0_real64, 'conversion.cvr', errors)
         call get_optional_real_with_default(conv, 'cvs', config%conversion%cvs, 0.0_real64, 'conversion.cvs', errors)
      end if

      call get_table(doc_root, 'respiration', resp, 'respiration', errors)
      if (associated(resp)) then
         call get_optional_real_with_default(resp, 'q10', config%respiration%q10, 0.0_real64, 'respiration.q10', errors)
         call get_optional_real_with_default(resp, 'rml', config%respiration%rml, 0.0_real64, 'respiration.rml', errors)
         call get_optional_real_with_default(resp, 'rmo', config%respiration%rmo, 0.0_real64, 'respiration.rmo', errors)
         call get_optional_real_with_default(resp, 'rmr', config%respiration%rmr, 0.0_real64, 'respiration.rmr', errors)
         call get_optional_real_with_default(resp, 'rms', config%respiration%rms, 0.0_real64, 'respiration.rms', errors)
         call read_table_2d(resp, 'rfsetb', config%respiration%rfsetb, 2, 'respiration.rfsetb', errors)
      end if

      call get_table(doc_root, 'partitioning', part, 'partitioning', errors)
      if (associated(part)) then
         call read_table_2d(part, 'frtb', config%partitioning%frtb, 2, 'partitioning.frtb', errors)
         call read_table_2d(part, 'fltb', config%partitioning%fltb, 2, 'partitioning.fltb', errors)
         call read_table_2d(part, 'fstb', config%partitioning%fstb, 2, 'partitioning.fstb', errors)
         call read_table_2d(part, 'fotb', config%partitioning%fotb, 2, 'partitioning.fotb', errors)
      end if

      call get_table(doc_root, 'death', deth, 'death', errors)
      if (associated(deth)) then
         call get_optional_real_with_default(deth, 'perdl', config%death%perdl, 0.0_real64, 'death.perdl', errors)
         call read_table_2d(deth, 'rdrrtb', config%death%rdrrtb, 2, 'death.rdrrtb', errors)
         call read_table_2d(deth, 'rdrstb', config%death%rdrstb, 2, 'death.rdrstb', errors)
      end if

      call get_table(doc_root, 'root', root, 'root', errors)
      if (associated(root)) then
         call get_optional_int_with_default(root,  'swrd',     config%root%swrd,     0,          'root.swrd',     errors)
         call get_optional_real_with_default(root, 'rdi',      config%root%rdi,      0.0_real64, 'root.rdi',      errors)
         call get_optional_real_with_default(root, 'rri',      config%root%rri,      0.0_real64, 'root.rri',      errors)
         call get_optional_real_with_default(root, 'rdc',      config%root%rdc,      0.0_real64, 'root.rdc',      errors)
         call get_optional_int_with_default(root,  'swdmi2rd', config%root%swdmi2rd, 0,          'root.swdmi2rd', errors)
         call get_optional_int_with_default(root,  'swrdc',   config%root%swrdc,   0,           'root.swrdc',    errors)
         call get_optional_real_with_default(root, 'wrtmax',   config%root%wrtmax,   0.0_real64, 'root.wrtmax',   errors)
         call read_table_2d(root, 'rdtb',  config%root%rdtb,  2, 'root.rdtb',  errors)
         call read_table_2d(root, 'rlwtb', config%root%rlwtb, 2, 'root.rlwtb', errors)
         call read_table_2d(root, 'rdctb', config%root%rdctb, 2, 'root.rdctb', errors)
      end if

      call get_table(doc_root, 'oxygen_stress', oxy, 'oxygen_stress', errors)
      if (associated(oxy)) then
         call get_optional_int_with_default(oxy,  'swoxygen',                config%oxygen_stress%swoxygen,                0,          'oxygen_stress.swoxygen',                errors)
         call get_optional_int_with_default(oxy,  'swwrtnonox',              config%oxygen_stress%swwrtnonox,              0,          'oxygen_stress.swwrtnonox',              errors)
         call get_optional_real_with_default(oxy, 'aeratecrit',              config%oxygen_stress%aeratecrit,              0.0_real64, 'oxygen_stress.aeratecrit',              errors)
         call get_optional_real_with_default(oxy, 'hlim1',                   config%oxygen_stress%hlim1,                   0.0_real64, 'oxygen_stress.hlim1',                   errors)
         call get_optional_real_with_default(oxy, 'hlim2u',                  config%oxygen_stress%hlim2u,                  0.0_real64, 'oxygen_stress.hlim2u',                  errors)
         call get_optional_real_with_default(oxy, 'hlim2l',                  config%oxygen_stress%hlim2l,                  0.0_real64, 'oxygen_stress.hlim2l',                  errors)
         call get_optional_real_with_default(oxy, 'q10_microbial',           config%oxygen_stress%q10_microbial,           0.0_real64, 'oxygen_stress.q10_microbial',           errors)
         call get_optional_real_with_default(oxy, 'specific_resp_humus',     config%oxygen_stress%specific_resp_humus,     0.0_real64, 'oxygen_stress.specific_resp_humus',     errors)
         call get_optional_real_with_default(oxy, 'srl',                     config%oxygen_stress%srl,                     0.0_real64, 'oxygen_stress.srl',                     errors)
         call get_optional_int_with_default(oxy,  'swrootradius',            config%oxygen_stress%swrootradius,            0,          'oxygen_stress.swrootradius',            errors)
         call get_optional_real_with_default(oxy, 'dry_mat_cont_roots',      config%oxygen_stress%dry_mat_cont_roots,      0.0_real64, 'oxygen_stress.dry_mat_cont_roots',      errors)
         call get_optional_real_with_default(oxy, 'air_filled_root_por',     config%oxygen_stress%air_filled_root_por,     0.0_real64, 'oxygen_stress.air_filled_root_por',     errors)
         call get_optional_real_with_default(oxy, 'spec_weight_root_tissue', config%oxygen_stress%spec_weight_root_tissue, 0.0_real64, 'oxygen_stress.spec_weight_root_tissue', errors)
         call get_optional_real_with_default(oxy, 'var_a',                   config%oxygen_stress%var_a,                   0.0_real64, 'oxygen_stress.var_a',                   errors)
         call get_optional_real_with_default(oxy, 'root_radiusO2',           config%oxygen_stress%root_radiusO2,           0.0_real64, 'oxygen_stress.root_radiusO2',           errors)
      end if

      call get_table(doc_root, 'drought_stress', drou, 'drought_stress', errors)
      if (associated(drou)) then
         call get_optional_int_with_default(drou,  'swdrought', config%drought_stress%swdrought, 0,          'drought_stress.swdrought', errors)
         call get_optional_real_with_default(drou, 'hlim3h',    config%drought_stress%hlim3h,    0.0_real64, 'drought_stress.hlim3h',    errors)
         call get_optional_real_with_default(drou, 'hlim3l',    config%drought_stress%hlim3l,    0.0_real64, 'drought_stress.hlim3l',    errors)
         call get_optional_real_with_default(drou, 'hlim4',     config%drought_stress%hlim4,     0.0_real64, 'drought_stress.hlim4',     errors)
         call get_optional_real_with_default(drou, 'adcrh',     config%drought_stress%adcrh,     0.0_real64, 'drought_stress.adcrh',     errors)
         call get_optional_real_with_default(drou, 'adcrl',     config%drought_stress%adcrl,     0.0_real64, 'drought_stress.adcrl',     errors)
      end if

      call get_table(doc_root, 'salinity', salt, 'salinity', errors)
      if (associated(salt)) then
         call get_optional_int_with_default(salt,  'swsalinity', config%salinity%swsalinity, 0,          'salinity.swsalinity', errors)
         call get_optional_real_with_default(salt, 'saltmax',    config%salinity%saltmax,    0.0_real64, 'salinity.saltmax',    errors)
         call get_optional_real_with_default(salt, 'saltslope',  config%salinity%saltslope,  0.0_real64, 'salinity.saltslope',  errors)
         call get_optional_real_with_default(salt, 'salthead',   config%salinity%salthead,   0.0_real64, 'salinity.salthead',   errors)
      end if

      call get_table(doc_root, 'compensate', comp, 'compensate', errors)
      if (associated(comp)) then
         call get_optional_int_with_default(comp,  'swcompensate', config%compensate%swcompensate, 0,          'compensate.swcompensate', errors)
         call get_optional_int_with_default(comp,  'swstressor',   config%compensate%swstressor,   0,          'compensate.swstressor',   errors)
         call get_optional_real_with_default(comp, 'alphacrit',    config%compensate%alphacrit,    0.0_real64, 'compensate.alphacrit',    errors)
         call get_optional_real_with_default(comp, 'dcritrtz',     config%compensate%dcritrtz,     0.0_real64, 'compensate.dcritrtz',     errors)
      end if

      call get_table(doc_root, 'interception', inter, 'interception', errors)
      if (associated(inter)) then
         call get_optional_int_with_default(inter,  'swinter', config%interception%swinter, 0,          'interception.swinter', errors)
         call get_optional_real_with_default(inter, 'cofab',   config%interception%cofab,   0.0_real64, 'interception.cofab',   errors)
         call read_table_2d(inter, 'gashtb', config%interception%gashtb, 6, 'interception.gashtb', errors)
      end if

      call get_table(doc_root, 'co2', co2t, 'co2', errors)
      if (associated(co2t)) then
         call get_optional_int_with_default(co2t,    'swco2',   config%co2%swco2,   0,    'co2.swco2',   errors)
         call get_optional_string_with_default(co2t, 'atmofil', config%co2%atmofil, '',   'co2.atmofil', errors)
         call read_table_2d(co2t, 'co2amaxtb', config%co2%co2amaxtb, 2, 'co2.co2amaxtb', errors)
         call read_table_2d(co2t, 'co2efftb',  config%co2%co2efftb,  2, 'co2.co2efftb',  errors)
         call read_table_2d(co2t, 'co2tratb',  config%co2%co2tratb,  2, 'co2.co2tratb',  errors)
      end if

      call get_table(doc_root, 'management', mgmt, 'management', errors)
      if (associated(mgmt)) then
         call get_optional_real_with_default(mgmt, 'fraharlosorm_lv',     config%management%fraharlosorm_lv,     0.0_real64, 'management.fraharlosorm_lv',     errors)
         call get_optional_real_with_default(mgmt, 'fraharlosorm_st',     config%management%fraharlosorm_st,     0.0_real64, 'management.fraharlosorm_st',     errors)
         call get_optional_real_with_default(mgmt, 'fraharlosorm_so',     config%management%fraharlosorm_so,     0.0_real64, 'management.fraharlosorm_so',     errors)
         call get_optional_real_with_default(mgmt, 'fradeceasedlvtosoil', config%management%fradeceasedlvtosoil, 0.0_real64, 'management.fradeceasedlvtosoil', errors)
         call get_optional_int_with_default(mgmt,  'swpotrelmf',          config%management%swpotrelmf,          0,          'management.swpotrelmf',          errors)
         call get_optional_real_with_default(mgmt, 'relmf',               config%management%relmf,               0.0_real64, 'management.relmf',               errors)
      end if

      ! ----------------------------------------------------------------
      ! Soybean variant (optional section; absent in most cases)
      ! ----------------------------------------------------------------
      call get_table(doc_root, 'soybean', soy, 'soybean', errors)
      if (associated(soy)) then
         call get_optional_int_with_default(soy, 'swsoybean',    config%soybean%swsoybean,    0,          'soybean.swsoybean',    errors)
         call get_optional_real_with_default(soy, 'mg',          config%soybean%mg,           0.0_real64, 'soybean.mg',           errors)
         call get_optional_real_with_default(soy, 'dvsi',        config%soybean%dvsi,         0.0_real64, 'soybean.dvsi',         errors)
         call get_optional_real_with_default(soy, 'dvrmax1',     config%soybean%dvrmax1,      0.0_real64, 'soybean.dvrmax1',      errors)
         call get_optional_real_with_default(soy, 'dvrmax2',     config%soybean%dvrmax2,      0.0_real64, 'soybean.dvrmax2',      errors)
         call get_optional_real_with_default(soy, 'tmaxdvr',     config%soybean%tmaxdvr,      0.0_real64, 'soybean.tmaxdvr',      errors)
         call get_optional_real_with_default(soy, 'tmindvr',     config%soybean%tmindvr,      0.0_real64, 'soybean.tmindvr',      errors)
         call get_optional_real_with_default(soy, 'toptdvr',     config%soybean%toptdvr,      0.0_real64, 'soybean.toptdvr',      errors)
         call get_optional_real_with_default(soy, 'popt',        config%soybean%popt,         0.0_real64, 'soybean.popt',         errors)
         call get_optional_real_with_default(soy, 'pcrt',        config%soybean%pcrt,         0.0_real64, 'soybean.pcrt',         errors)
      end if

      ! ----------------------------------------------------------------
      ! Bulb crops (optional section; absent in most cases)
      ! ----------------------------------------------------------------
      call get_table(doc_root, 'bulb', bul, 'bulb', errors)
      if (associated(bul)) then
         call get_optional_int_with_default(bul,  'swbulb', config%bulb%swbulb, 0,          'bulb.swbulb', errors)
         call get_optional_real_with_default(bul, 'pld',    config%bulb%pld,    0.0_real64, 'bulb.pld',    errors)
         call get_optional_real_with_default(bul, 'plwti',  config%bulb%plwti,  0.0_real64, 'bulb.plwti',  errors)
         call get_optional_real_with_default(bul, 'remoc',  config%bulb%remoc,  0.0_real64, 'bulb.remoc',  errors)
         call read_table_2d(bul, 'fbltb', config%bulb%fbltb, 2, 'bulb.fbltb', errors)
      end if

      ! ----------------------------------------------------------------
      ! Nutrient model (optional section; absent in most cases)
      ! ----------------------------------------------------------------
      call get_table(doc_root, 'nutrient', nut, 'nutrient', errors)
      if (associated(nut)) then
         call get_optional_real_with_default(nut, 'lrnr',   config%nutrient%lrnr,   0.0_real64, 'nutrient.lrnr',   errors)
         call get_optional_real_with_default(nut, 'lsnr',   config%nutrient%lsnr,   0.0_real64, 'nutrient.lsnr',   errors)
         call get_optional_real_with_default(nut, 'nlai',   config%nutrient%nlai,   0.0_real64, 'nutrient.nlai',   errors)
         call get_optional_real_with_default(nut, 'nlue',   config%nutrient%nlue,   0.0_real64, 'nutrient.nlue',   errors)
         call get_optional_real_with_default(nut, 'nmaxso', config%nutrient%nmaxso, 0.0_real64, 'nutrient.nmaxso', errors)
         call get_optional_real_with_default(nut, 'npart',  config%nutrient%npart,  0.0_real64, 'nutrient.npart',  errors)
         call get_optional_real_with_default(nut, 'nfixf',  config%nutrient%nfixf,  0.0_real64, 'nutrient.nfixf',  errors)
         call get_optional_real_with_default(nut, 'nsla',   config%nutrient%nsla,   0.0_real64, 'nutrient.nsla',   errors)
         call get_optional_real_with_default(nut, 'rnflv',  config%nutrient%rnflv,  0.0_real64, 'nutrient.rnflv',  errors)
         call get_optional_real_with_default(nut, 'rnfrt',  config%nutrient%rnfrt,  0.0_real64, 'nutrient.rnfrt',  errors)
         call get_optional_real_with_default(nut, 'rnfst',  config%nutrient%rnfst,  0.0_real64, 'nutrient.rnfst',  errors)
         call get_optional_real_with_default(nut, 'tcnt',   config%nutrient%tcnt,   0.0_real64, 'nutrient.tcnt',   errors)
         call get_optional_real_with_default(nut, 'dvsnlt', config%nutrient%dvsnlt, 0.0_real64, 'nutrient.dvsnlt', errors)
         call get_optional_real_with_default(nut, 'dvsnt',  config%nutrient%dvsnt,  0.0_real64, 'nutrient.dvsnt',  errors)
         call get_optional_real_with_default(nut, 'rdrns',  config%nutrient%rdrns,  0.0_real64, 'nutrient.rdrns',  errors)
         call get_optional_real_with_default(nut, 'fntrt',  config%nutrient%fntrt,  0.0_real64, 'nutrient.fntrt',  errors)
         call get_optional_real_with_default(nut, 'frnx',   config%nutrient%frnx,   0.0_real64, 'nutrient.frnx',   errors)
         call read_table_2d(nut, 'nmxlv', config%nutrient%nmxlv, 2, 'nutrient.nmxlv', errors)
      end if

      call get_table(doc_root, 'irrigation_schedule', irr_sched, 'irrigation_schedule', errors)
      call read_irrigation_schedule_from_section(irr_sched, config%schedule, errors)
   end subroutine read_cropwofost_toml

   !> Decode a TOML array-of-arrays at sec[key] into a (nrows, ncols)
   !! real(real64) allocatable. Absent key leaves table unallocated.
   !! Ragged or wrong-width inner arrays append a parse-type-mismatch error
   !! and leave table unallocated.
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

end module read_cropwofost_toml_mod
