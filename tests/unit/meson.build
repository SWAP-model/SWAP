# ============================================================================
# Unit test build configuration
# ============================================================================
# Only sources that exist at the rescue baseline (commit e256bc0) are
# referenced here. Drift-era forward declarations for state modules and
# TOML readers that were never actually committed have been removed.
#
# At the rescue baseline, every .pf file under tests/unit/ referenced at
# least one module that was never committed (e.g. swap_config_mod,
# time_state_mod, irrigation_state_mod, toml_helpers_mod). The .pf files
# were therefore removed as part of Phase 1 Task 2. Git history preserves
# them if Phase 4 wants to revive individual suites once the missing
# state modules are reintroduced alongside each module's TDD sub-task.
#
# The pFUnit wiring itself (dependency discovery, generator, driver) is
# kept in place so that Phase 4 can re-add .pf files one at a time with
# a minimal diff here.
# ============================================================================

test_base_sources = [
    '../../src/core/arrays.f90',
    '../../src/core/constants.f90',
    '../../src/core/swap_log.f90',
    '../../src/core/dtutil.f90',
    '../../src/utils/arrayutils.f90',
]

# ============================================================================
# pFUnit test suite (gfortran only)
# ============================================================================
# The .pf files list is intentionally empty at the rescue baseline — see
# the header comment. The block below discovers pFUnit and sets up the
# preprocessor generator so Phase 4 can add .pf files with a one-line
# change to pf_files (plus the matching ADD_TEST_SUITE line in
# testSuites.inc and the required source modules in pfunit_base_sources).
# ============================================================================
if get_option('enable_pfunit')
    python3 = import('python').find_installation('python3')

    pfunit_root_result = run_command('sh', '-c', 'echo ${PFUNIT_ROOT:-}', check: false)
    pfunit_root = pfunit_root_result.stdout().strip()

    if pfunit_root == ''
        pfunit_root = join_paths(meson.current_source_dir(),
            '../../tests/pFUnit/build/install_gfortran/PFUNIT-4.15')
        if not fs.exists(pfunit_root)
            error('pFUnit not found. Set PFUNIT_ROOT to a PFUNIT-4.15 install, or build one under tests/pFUnit/build/install_gfortran/PFUNIT-4.15.')
        endif
    endif

    message('Using pFUnit from: ' + pfunit_root)

    pfunit_dep = dependency('PFUNIT',
        method: 'cmake',
        required: true,
        cmake_args: ['-DCMAKE_PREFIX_PATH=' + pfunit_root]
    )

    funitproc = join_paths(pfunit_root, 'bin', 'funitproc')
    if not fs.exists(funitproc)
        error('funitproc not found at: ' + funitproc)
    endif

    pfunit_pp = generator(python3,
        output: '@BASENAME@.F90',
        arguments: [funitproc, '@INPUT@', '@OUTPUT@']
    )

    driver_src = join_paths(pfunit_root, 'include', 'driver.F90')

    # Additional sources beyond test_base_sources required by surviving
    # .pf files. Empty at baseline; grows as Phase 4 re-adds suites.
    pfunit_extra_sources = [
        '../../src/atmosphere/precipitation.f90',
        '../../src/core/variables.f90',
        '../../src/crop/wofost_soil_declarations.f90',
        '../../src/crop/wofost_soil_interface.f90',
        '../../src/atmosphere/et.f90',
        '../../src/error/error.f90',
        '../../src/error/fatalerr.f90',
        '../../src/validation/validation.f90',
        '../../src/config/general_config.f90',
        '../../src/config/simulation_config.f90',
        '../../src/config/meteorology_config.f90',
        '../../src/config/drainage_config.f90',
        '../../src/config/soil_config.f90',
        '../../src/config/bottom_boundary_config.f90',
        '../../src/config/heat_config.f90',
        '../../src/config/irrigation_config.f90',
        '../../src/config/solute_config.f90',
        '../../src/config/surface_water_config.f90',
        '../../src/config/nutrients_config.f90',
        '../../src/config/cropfixed_config.f90',
        '../../src/config/cropgrass_config.f90',
        '../../src/config/cropwofost_config.f90',
        '../../src/config/crop_config.f90',
        '../../src/config/output_csv_config.f90',
        '../../src/crop/crop_config_global.f90',
        '../../src/crop/cropfixed_init.f90',
        '../../src/crop/cropwofost_init.f90',
        '../../src/crop/cropgrass_init.f90',
        '../../src/config/swap_config.f90',
        '../../src/io/csv_reader.f90',
        '../../src/io/toml/path_helpers.f90',
        '../../src/io/toml/toml_field_helpers.f90',
        '../../src/io/toml/read_general_toml.f90',
        '../../src/io/toml/read_simulation_toml.f90',
        '../../src/io/toml/read_meteorology_toml.f90',
        '../../src/io/toml/read_drainage_toml.f90',
        '../../src/io/toml/read_bottom_boundary_toml.f90',
        '../../src/io/toml/read_heat_toml.f90',
        '../../src/io/toml/read_irrigation_toml.f90',
        '../../src/io/toml/read_solute_toml.f90',
        '../../src/io/toml/read_soil_toml.f90',
        '../../src/io/toml/read_soil_tillage_toml.f90',
        '../../src/io/toml/read_irrigation_ssdi_toml.f90',
        '../../src/io/toml/read_surface_water_toml.f90',
        '../../src/io/toml/read_cropfixed_toml.f90',
        '../../src/io/toml/read_cropgrass_toml.f90',
        '../../src/io/toml/read_cropwofost_toml.f90',
        '../../src/io/toml/read_crop_toml.f90',
        '../../src/io/toml/read_output_csv_toml.f90',
        '../../src/io/toml/read_nutrients_toml.f90',
        '../../src/io/toml/load_swap_config.f90',
        '../../src/io/toml/load_swap_config_string.f90',
        '../../src/io/toml/write_swap_config.f90',
        '../../src/io/toml/config_to_variables.f90',
        'io/toml/chdir_helper.f90',
        # swap_mod.f90 (added Task 3) needs the full compute-module chain compiled in the test binary
        '../../src/soil/WC_K_models_04_11.f90',
        '../../src/utils/soilhydraulicsutils.f90',
        '../../src/soil/sptabulated.f90',
        '../../src/crop/oxygenstress.f90',
        '../../src/utils/surfacewaterutils.f90',
        '../../src/state/surfacewater_state.f90',
        '../../src/state/drainage_state.f90',
        '../../src/state/solute_state.f90',
        '../../src/state/heat_state.f90',
        '../../src/state/hydraulic_params_mod.f90',
        '../../src/state/soilwater_state.f90',
        '../../src/state/atmosphere_state.f90',
        '../../src/state/tillage_state.f90',
        '../../src/state/timecontrol_state.f90',
        '../../src/state/swap_state.f90',
        '../../src/core/timestep_control_mod.f90',
        '../../src/core/initialize.f90',
        '../../src/soil/soilgrid.f90',
        '../../src/soil/waterbalance.f90',
        '../../src/soil/soilhydraulics.f90',
        '../../src/io/swap_csv_output.f90',
        '../../src/io/meteo_buffer_mod.f90',
        '../../src/io/readmeteo.f90',
        '../../src/io/swapoutput.f90',
        '../../src/atmosphere/interception.f90',
        '../../src/atmosphere/meteoday.f90',
        '../../src/atmosphere/meteodt.f90',
        '../../src/atmosphere/snow.f90',
        '../../src/heat/temperature.f90',
        '../../src/heat/frozencond.f90',
        '../../src/boundary/boundtop.f90',
        '../../src/boundary/boundbottom.f90',
        '../../src/drainage/drainage.f90',
        '../../src/drainage/divdra.f90',
        '../../src/drainage/surfacewater.f90',
        '../../src/solute/solute.f90',
        '../../src/solute/agetracer.f90',
        '../../src/crop/wofost_soil_parameters.f90',
        '../../src/crop/wofost_soil_rateconstants.f90',
        '../../src/crop/wofost_soil_orgmatn.f90',
        '../../src/crop/wofost_soil_watern.f90',
        '../../src/crop/wofost_soil_amendments.f90',
        '../../src/crop/wofost_soil_cropresidues.f90',
        '../../src/crop/wofost_soil_balancecheck.f90',
        '../../src/crop/wofostnut.f90',
        '../../src/crop/rootextraction.f90',
        '../../src/crop/tillage.f90',
        '../../src/crop/management_soil.f90',
        '../../src/crop/cropgrowth.f90',
        '../../src/utils/sharedsimulation.f90',
        '../../src/utils/numericalsolvers.f90',
        '../../src/core/swap_mod.f90',
        '../../src/crop/irrigation.f90',
        '../../src/io/file_io.f90',
        '../../src/core/timecontrol_mod.f90',
    ]

    pfunit_base_sources = test_base_sources + pfunit_extra_sources

    pf_files = [
        'utils/test_arrayutils.pf',
        'atmosphere/test_precipitation.pf',
        'atmosphere/test_et.pf',
        'error/test_error.pf',
        'validation/test_validation.pf',
        'config/test_general_config.pf',
        'config/test_simulation_config.pf',
        'config/test_meteorology_config.pf',
        'config/test_drainage_config.pf',
        'config/test_soil_config.pf',
        'config/test_bottom_boundary_config.pf',
        'config/test_heat_config.pf',
        'config/test_irrigation_config.pf',
        'config/test_solute_config.pf',
        'config/test_solute_phase0_fields.pf',
        'config/test_heat_phase0_fields.pf',
        'config/test_boundary_phase0_fields.pf',
        'config/test_surface_water_config.pf',
        'config/test_crop_config.pf',
        'config/test_output_csv_config.pf',
        'config/test_swap_config.pf',
        'config/test_cropfixed_config.pf',
        'config/test_cropgrass_config.pf',
        'config/test_cropwofost_config.pf',
        'io/test_csv_reader.pf',
        'io/test_file_io.pf',
        'io/toml/test_path_helpers.pf',
        'io/toml/test_toml_field_helpers.pf',
        'io/toml/test_read_general_toml.pf',
        'io/toml/test_read_simulation_toml.pf',
        'io/toml/test_read_meteorology_toml.pf',
        'io/toml/test_read_drainage_toml.pf',
        'io/toml/test_read_bottom_boundary_toml.pf',
        'io/toml/test_read_heat_toml.pf',
        'io/toml/test_read_irrigation_toml.pf',
        'io/toml/test_read_solute_toml.pf',
        'io/toml/test_read_soil_toml.pf',
        'io/toml/test_read_surface_water_toml.pf',
        'io/toml/test_read_output_csv_toml.pf',
        'io/toml/test_read_crop_toml.pf',
        'io/toml/test_read_cropfixed_toml.pf',
        'io/toml/test_read_cropgrass_toml.pf',
        'io/toml/test_read_cropwofost_toml.pf',
        'io/toml/test_config_to_variables.pf',
        'io/toml/test_load_swap_config.pf',
        'io/toml/test_hupselbrook_loads.pf',
        'io/toml/test_hupselbrook_roundtrip.pf',
        'io/toml/test_all_cases_smoke.pf',
        'io/toml/test_chdir_helper.pf',
        'io/toml/test_hupselbrook_parity.pf',
        'io/toml/test_grassgrowth_parity.pf',
        'io/toml/test_oxygenstress_parity.pf',
        'io/toml/test_salinitystress_parity.pf',
        'io/toml/test_surfacewater_parity.pf',
        # DEFERRED — paired with the test_macroporeflow_parity_suite
        # entry in testSuites.inc; see ADR 0011 / 0018.
        # 'io/toml/test_macroporeflow_parity.pf',
        'io/toml/test_apply_soil_tillage.pf',
        'io/toml/test_soil_tillage_validate.pf',
        'io/toml/test_read_soil_tillage_toml.pf',
        'io/toml/test_apply_irrigation_ssdi.pf',
        'io/toml/test_irrigation_ssdi_validate.pf',
        'io/toml/test_read_irrigation_ssdi_toml.pf',
        'io/toml/test_apply_cropwofost_nutrient.pf',
        'io/toml/test_nutrients_config.pf',
        'state/test_surfacewater_state.pf',
        'state/test_surfacewater_cumulative.pf',
        'state/test_drainage_state.pf',
        'state/test_solute_state.pf',
        'state/test_solute_cumulative.pf',
        'state/test_heat_state.pf',
        'state/test_soilwater_state.pf',
        'state/test_swap_state.pf',
        'state/test_atmosphere_state.pf',
        'state/test_tillage_state.pf',
        'state/test_timecontrol_state.pf',
        'crop/test_cropfixed_init.pf',
        'crop/test_cropwofost_init.pf',
        'crop/test_cropgrass_init.pf',
        'core/test_swap_mod.pf',
        'core/test_timecontrol_mod.pf',
    ]

    # When pf_files is non-empty, build and register the unit-swap-tests
    # executable. At baseline this is skipped; Phase 4 will flip the
    # switch by adding at least one .pf file above.
    if pf_files.length() > 0
        pp_sources = pfunit_pp.process(pf_files)

        unit_tests = executable('unit-swap-tests',
            sources: [pp_sources, driver_src] + pfunit_base_sources,
            dependencies: [pfunit_dep, tomlf_dep],
            include_directories: src_inc,
            fortran_args: ['-D_TEST_SUITES="testSuites.inc"'],
            install: false
        )

        # Run the binary directly. The per-suite-isolation wrapper
        # (run_pfunit.sh) was retired with the legacy readers — there is
        # no shared global state to leak between suites anymore.
        test('unit-swap-tests',
            unit_tests,
            suite: 'unit-pfunit',
            workdir: meson.project_source_root(),
            timeout: 120
        )
    else
        message('No pFUnit .pf files present — skipping unit-swap-tests executable. This is expected at the Phase 1 rescue baseline.')
    endif
endif
