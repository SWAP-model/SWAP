# Minimal pFUnit glue without relying on the missing Meson pFUnit module.
# Requires PFUNIT_ROOT to point at the pFUnit install prefix that contains
# bin/funitproc and cmake/PFUNITConfig.cmake

python3 = import('python').find_installation('python3')

# Get PFUNIT_ROOT environment variable, defaulting to gfortran build in tests/pFUnit
pfunit_root_result = run_command('sh', '-c', 'echo ${PFUNIT_ROOT:-}', check: false)
pfunit_root = pfunit_root_result.stdout().strip()

if pfunit_root == ''
  # Default to local gfortran build
  pfunit_root = join_paths(meson.current_source_dir(), '../../tests/pFUnit/build/install_gfortran/PFUNIT-4.15')
  if not fs.exists(pfunit_root)
    error('pFUnit not found. Build it first: cd tests/pFUnit && FC=gfortran cmake -B build -DCMAKE_INSTALL_PREFIX=$PWD/build/install_gfortran -DCMAKE_BUILD_TYPE=Release -DSKIP_MPI=YES -DSKIP_OPENMP=YES && cmake --build build -j$(nproc) && cmake --install build')
  endif
endif

message('Using pFUnit from: ' + pfunit_root)

# Find pFUnit using CMake with CMAKE_PREFIX_PATH
pfunit_dep = dependency(
  'PFUNIT',
  method: 'cmake',
  required: true,
  cmake_args: ['-DCMAKE_PREFIX_PATH=' + pfunit_root],
)

funitproc = join_paths(pfunit_root, 'bin', 'funitproc')
if not fs.exists(funitproc)
  error('funitproc not found at ' + funitproc + '; check PFUNIT_ROOT or rebuild pFUnit')
endif

# Generator to preprocess .pf files to .F90
pfunit_pp = generator(
  python3,
  output: '@BASENAME@.F90',
  arguments: [funitproc, '@INPUT@', '@OUTPUT@'],
)

# Preprocess test files
pp_sources = pfunit_pp.process(
  'test_functions.pf',
  'test_ttutil.pf',
  'test_van_genuchten.pf',
  'test_datetime.pf'
)

# Copy pFUnit driver template and customize it
driver_src = join_paths(pfunit_root, 'include', 'driver.F90')

# Add SWAP source files needed for testing (standalone pure functions)
swap_pure_functions = files('swap_pure_functions.f90')

# Build test executable with pFUnit driver
unit_tests = executable(
  'unit-swap-tests',
  sources: [
    swap_pure_functions,
    pp_sources,
    driver_src,
  ],
  dependencies: [pfunit_dep, ttutil_dep],
  include_directories: src_inc,
  fortran_args: ['-D_TEST_SUITES="testSuites.inc"'],
  install: false,
)

# Register test with Meson
test('unit-swap-tests', unit_tests)
