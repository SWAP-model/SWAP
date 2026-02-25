# ============================================================================
# Shared base sources (no TOML dependency)
# ============================================================================
test_base_sources = [
    '../../src/core/arrays.f90',
    '../../src/utils/arrayutils.f90',
    '../../src/core/swap_log.f90',
    '../../src/core/time_state.f90',
    '../../src/atmosphere/snow_state.f90',
    '../../src/atmosphere/atmosphere_state.f90',
    '../../src/boundary/boundary_state.f90',
    '../../src/drainage/drainage_state.f90',
    '../../src/drainage/surfacewater_state.f90',
    '../../src/heat/heat_state.f90',
    '../../src/heat/frozencond_state.f90',
    '../../src/solute/solute_state.f90',
    '../../src/soil/soil_state.f90',
    '../../src/macropore/macropore_state.f90',
    '../../src/crop/irrigation_state.f90',
    '../../src/crop/tillage_state.f90',
    '../../src/crop/oxygenstress_state.f90',
    '../../src/crop/management_soil_state.f90',
    '../../src/crop/rootextraction_state.f90',
    '../../src/crop/cropgrowth_state.f90',
    '../../src/crop/wofost_state.f90',
    '../../src/core/simulation_config.f90',
    '../../src/core/output_config.f90',
    '../../src/core/swap_config.f90',
    '../../src/core/swap_io.f90',
    '../../src/core/swap_state.f90'
]

# ============================================================================
# State lifecycle tests moved to pFUnit suite (see pFUnit section below)
# ============================================================================

# ============================================================================
# I/O reader tests (TOML dependency, no pFUnit)
# Each reader gets its own executable with minimal sources to keep
# link times fast and failure attribution clear.
# ============================================================================

# readdra_toml: minimal — only what drainage reading needs
readdra_toml_test_exe = executable('test_readdra_toml',
    [
        '../../src/core/swap_log.f90',
        '../../src/drainage/drainage_state.f90',
        '../../src/io/readdra_toml.f90',
        'io/test_readdra_toml.f90'
    ],
    dependencies: [tomlf_dep, ttutil_dep],
    include_directories: src_inc,
    install: false
)
test('test_readdra_toml', readdra_toml_test_exe,
    suite: 'unit-standalone',
    workdir: meson.project_source_root()
)

# readcrop_toml: needs crop states + crop_toml_state
readcrop_toml_test_exe = executable('test_readcrop_toml',
    test_base_sources + [
        '../../src/utils/tomlutils.f90',
        '../../src/io/readcrop_toml.f90',
        'io/test_readcrop_toml.f90'
    ],
    dependencies: [tomlf_dep, ttutil_dep],
    include_directories: src_inc,
    install: false
)
test('test_readcrop_toml', readcrop_toml_test_exe,
    suite: 'unit-standalone',
    workdir: meson.project_source_root()
)

# readcrop_toml section readers: standalone test for selected subroutines
readcrop_sections_test_exe = executable('test_readcrop_sectionreaders',
    test_base_sources + [
        '../../src/utils/tomlutils.f90',
        '../../src/io/readcrop_toml.f90',
        'io/test_readcrop_sectionreaders.f90'
    ],
    dependencies: [tomlf_dep, ttutil_dep],
    include_directories: src_inc,
    install: false
)
test('test_readcrop_sectionreaders', readcrop_sections_test_exe,
    suite: 'unit-standalone',
    workdir: meson.project_source_root()
)

# ============================================================================
# pFUnit test suite (gfortran only — not compatible with ifx)
# ============================================================================
if get_option('enable_pfunit')
    python3 = import('python').find_installation('python3')

    pfunit_root_result = run_command('sh', '-c', 'echo ${PFUNIT_ROOT:-}', check: false)
    pfunit_root = pfunit_root_result.stdout().strip()

    if pfunit_root == ''
        pfunit_root = join_paths(meson.current_source_dir(),
            '../../pFUnit/build/install_gfortran/PFUNIT-4.15')
        if not fs.exists(pfunit_root)
            error('pFUnit not found at: ' + pfunit_root)
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

    # pFUnit sources mirror test_base_sources exactly
    pfunit_base_sources = test_base_sources + [
        '../../src/utils/tomlutils.f90',
        '../../src/io/readdra_toml.f90',
        '../../src/io/readcrop_toml.f90',
        '../../src/io/readswap_toml.f90',
    ]

    pf_files = files(
        'io/test_readswaptoml_full.pf',
        'io/test_sectionreaders_suite.pf',
        'state_lifecycles/test_state_lifecycles_suite.pf',
    )
    pp_sources = pfunit_pp.process(pf_files)

    unit_tests = executable('unit-swap-tests',
        sources: [pp_sources, driver_src] + pfunit_base_sources,
        dependencies: [pfunit_dep, ttutil_dep, tomlf_dep],
        include_directories: src_inc,
        fortran_args: ['-D_TEST_SUITES="testSuites.inc"'],
        install: false
    )

    test('unit-swap-tests', unit_tests,
        suite: 'unit-pfunit',
        workdir: meson.project_source_root()
    )
endif
