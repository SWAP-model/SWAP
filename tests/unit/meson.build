# ============================================================================
# Standalone unit tests (no pFUnit)
# ============================================================================
test_base_sources = [
    '../../src/core/swap_log.f90',
    '../../src/core/swap_state_mod.f90'
]

# State module tests
state_tests = {
    'test_state': ['core/test_state_standalone.f90'],
    'test_soil_state': ['soil/test_soil_state.f90'],
    'test_drainage_state': ['drainage/test_drainage_state.f90'],
    'test_solute_state': ['solute/test_solute_state.f90'],
    'test_heat_state': ['heat/test_heat_state.f90'],
    'test_atmosphere_state': ['atmosphere/test_atmosphere_state.f90'],
    'test_macropore_state': ['macropore/test_macropore_state.f90'],
    'test_tillage_state': ['crop/test_tillage_state.f90']
}

foreach test_name, test_files : state_tests
    executable(test_name, test_base_sources + test_files,
        include_directories: src_inc,
        install: false
    )
endforeach

# State sync test (requires additional sources)
executable('test_state_sync',
    test_base_sources + [
        '../../src/core/variables.f90',
        '../../src/crop/wofost_soil_declarations.f90',
        '../../src/crop/wofost_soil_interface.f90',
        '../../src/core/swap_state_sync.f90',
        'core/test_state_sync.f90'
    ],
    include_directories: src_inc,
    install: false
)

# ============================================================================
# pFUnit test suite
# ============================================================================
if get_option('enable_pfunit')
    python3 = import('python').find_installation('python3')

    # Get PFUNIT_ROOT (default to local gfortran build)
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

    # Find pFUnit dependency
    pfunit_dep = dependency('PFUNIT',
        method: 'cmake',
        required: true,
        cmake_args: ['-DCMAKE_PREFIX_PATH=' + pfunit_root]
    )

    # funitproc preprocessor
    funitproc = join_paths(pfunit_root, 'bin', 'funitproc')
    if not fs.exists(funitproc)
        error('funitproc not found at: ' + funitproc)
    endif

    # Generator for .pf → .F90
    pfunit_pp = generator(python3,
        output: '@BASENAME@.F90',
        arguments: [funitproc, '@INPUT@', '@OUTPUT@']
    )

    # Preprocess test files
    pp_sources = pfunit_pp.process(
        'test_functions.pf',
        'test_ttutil.pf',
        'test_van_genuchten.pf',
        'test_datetime.pf'
    )

    # Build test executable
    driver_src = join_paths(pfunit_root, 'include', 'driver.F90')

    unit_tests = executable('unit-swap-tests',
        sources: [pp_sources, driver_src],
        dependencies: [pfunit_dep, ttutil_dep],
        include_directories: src_inc,
        fortran_args: ['-D_TEST_SUITES="testSuites.inc"'],
        install: false
    )

    # Register with Meson test runner
    test('unit-swap-tests', unit_tests)
endif
