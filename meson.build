project('swap', 'fortran', version : '4.2.0', default_options : ['warning_level=0'])

fc = meson.get_compiler('fortran')
fs = import('fs')
src_inc = include_directories(
    'src',
    'src/core',
    'src/io',
    'src/soil',
    'src/atmosphere',
    'src/crop',
    'src/drainage',
    'src/boundary',
    'src/macropore',
    'src/solute',
    'src/heat',
    'src/utils'
)

# Intel Fortran compiler settings
if fc.get_id() == 'intel'
    add_project_arguments([
        '-O2',
        '-fpscomp', 'general',
        '-warn', 'declarations',
        '-warn', 'unused',
        '-warn', 'uncalled', 
        '-warn', 'interfaces',
        '-init=zero',
        '-fpe0',
        '-fp-model', 'source'
    ], language : 'fortran')
    
    # NOTE: Windows Intel Fortran (ifx) static linking remains problematic.
    # Despite multiple attempts (/MT, -static-intel, /NODEFAULTLIB, etc.),
    # the executable still requires Intel runtime DLLs (libmmd.dll).
    # This appears to be a known limitation of ifx on Windows.
    # Use MinGW cross-compilation for fully static Windows executables.
    
    # Commented out Windows Intel build - uncomment to re-enable testing
    # if host_machine.system() == 'windows'
    #     add_project_arguments([
    #         '/MT',                           # Static multithreaded C runtime
    #     ], language : 'fortran')
    #     add_project_link_arguments([
    #         '-static-intel'                  # Attempt static Intel libraries
    #     ], language : 'fortran')
    # else
    if host_machine.system() != 'windows'
        # Linux: Statically link Intel runtime libraries
        # System libraries (libc, libm) remain dynamically linked
        add_project_link_arguments([
            '-static-intel',    # Statically link Intel-specific libraries
        ], language : 'fortran')
    endif
    
# GCC/gfortran compiler settings for Linux
elif fc.get_id() == 'gcc' and host_machine.system() == 'linux'
    add_project_arguments([
        '-O2',
        '-ffree-line-length-none',  # Remove line length limit
        '-Wno-line-truncation',
        '-Wno-compare-reals', 
        '-Wno-missing-include-dirs',  # Ignore missing include directories
        '-Wno-unused-variable',       # Ignore unused variable warnings  
        '-Wno-unused-dummy-argument', # Ignore unused dummy argument warnings
        '-std=legacy',                # Allow legacy Fortran features
    ], language : 'fortran')
    
# MinGW Fortran compiler settings for Windows cross-compilation 
elif fc.get_id() == 'gcc' and host_machine.system() == 'windows'
    add_project_arguments([
        '-O2',
        '-ffree-line-length-none',  # Remove line length limit
        '-Wno-line-truncation',
        '-Wno-compare-reals', 
        '-Wno-missing-include-dirs',  # Ignore missing include directories
        '-Wno-unused-variable',       # Ignore unused variable warnings  
        '-Wno-unused-dummy-argument', # Ignore unused dummy argument warnings
        '-std=legacy',
        '-mconsole',
        '-static-libgfortran',
        '-static-libgcc'
    ], language : 'fortran')
    add_project_link_arguments(['-static', '-mconsole'], language : 'fortran')
endif

# Source files grouped by subfolder (organized by functional module)
sources = [
    # Core module - framework and state management
    'src/core/swap_log.f90',        # Logging (compile first, no dependencies)
    'src/core/swap_state_mod.f90',  # State types (depends on swap_log)
    'src/core/variables.f90',
    'src/core/swap_state_sync.f90', # State sync (depends on swap_state_mod, variables)
    'src/core/initialize.f90',
    'src/core/timecontrol.f90',
    
    # Soil module - physics, hydraulics, numerical solution
    'src/soil/sptabulated.f90',
    'src/soil/WC_K_models_04_11.f90',
    'src/soil/tridag.f90',
    'src/soil/functions.f90',
    'src/soil/hysteresis.f90',
    'src/soil/calcgrid.f90',
    'src/soil/calcgwl.f90',
    'src/soil/watstor.f90',
    'src/soil/fluxes.f90',
    'src/soil/headcalc.f90',
    'src/soil/soilwater.f90',
    'src/soil/integral.f90',
    'src/soil/checkmassbal.f90',
    'src/soil/convertdiscrvert.f90',
    
    # I/O module - input parsing and output generation
    'src/io/swap_csv_output.f90',
    'src/io/readswap.f90',
    'src/io/readmeteo.f90',
    'src/io/swapoutput.f90',
    
    # Atmosphere module - meteorology and ET
    'src/atmosphere/penmon.f90',
    'src/atmosphere/meteoday.f90',
    'src/atmosphere/meteodt.f90',
    'src/atmosphere/snow.f90',
    
    # Heat module - soil temperature
    'src/heat/temperature.f90',
    'src/heat/frozencond.f90',
    
    # Boundary module - top and bottom conditions
    'src/boundary/boundtop.f90',
    'src/boundary/boundbottom.f90',
    
    # Drainage module - lateral drainage and surface water
    'src/drainage/drainage.f90',
    'src/drainage/divdra.f90',
    'src/drainage/surfacewater.f90',
    
    # Macropore module - preferential flow
    'src/macropore/macropore.f90',
    'src/macropore/macrorate.f90',
    'src/macropore/macroporeoutput.f90',
    
    # Solute module - solute transport
    'src/solute/solute.f90',
    
    # Crop module - growth, nutrients, management
    'src/crop/wofost_soil_declarations.f90',
    'src/crop/wofost_soil_interface.f90',
    'src/crop/wofost_soil_parameters.f90',
    'src/crop/wofost_soil_rateconstants.f90',
    'src/crop/wofost_soil_orgmatn.f90',
    'src/crop/wofost_soil_watern.f90',
    'src/crop/wofost_soil_amendments.f90',
    'src/crop/wofost_soil_cropresidues.f90',
    'src/crop/wofost_soil_balancecheck.f90',
    'src/crop/wofostnut.f90',
    'src/crop/oxygenstress.f90',
    'src/crop/rootextraction.f90',
    'src/crop/irrigation.f90',
    'src/crop/tillage.f90',
    'src/crop/management_soil.f90',
    'src/crop/cropgrowth.f90',
    
    # Utils module - shared simulation utilities
    'src/utils/sharedexchange.f90',
    'src/utils/sharedsimulation.f90',
    
    # Core module - main model routines (must come last due to dependencies)
    'src/core/swap.f90',
    'src/core/swap_main.f90'
]

# Download ttutil library dependency from GitHub releases
ttutil_version = '4.2.7'
ttutil_lib_dir = join_paths(meson.current_build_dir(), 'lib')

# Determine platform-specific library name and path
if host_machine.system() == 'windows' and fc.get_id() != 'gcc'
    # Native Windows compilation with Intel Fortran
    ttutil_filename = 'libttutil@0@-windows.a'.format(ttutil_version)
    ttutil_lib_path = join_paths(ttutil_lib_dir, 'libttutil.a')
    ttutil_url = 'https://github.com/SWAP-model/ttutil/releases/download/v@0@/@1@'.format(ttutil_version, ttutil_filename)
    # Download ttutil library if not present
    if not fs.exists(ttutil_lib_path)
        run_command('mkdir', '-p', ttutil_lib_dir, check: true)
        run_command('curl', '-L', ttutil_url, '-o', ttutil_lib_path, check: true)
    endif
elif host_machine.system() == 'windows' and fc.get_id() == 'gcc'
    # MinGW cross-compilation 
    ttutil_filename = 'libttutil@0@-mingw.a'.format(ttutil_version)
    ttutil_lib_path = join_paths(ttutil_lib_dir, 'libttutil.a')
    ttutil_url = 'https://github.com/SWAP-model/ttutil/releases/download/v@0@/@1@'.format(ttutil_version, ttutil_filename)
    # Download ttutil library if not present
    if not fs.exists(ttutil_lib_path)
        run_command('mkdir', '-p', ttutil_lib_dir, check: true)
        run_command('curl', '-L', ttutil_url, '-o', ttutil_lib_path, check: true)
    endif
elif fc.get_id() == 'gcc' and host_machine.system() == 'linux'
    # Linux gfortran - use local ttutil from sibling directory
    ttutil_lib_path = join_paths(meson.current_source_dir(), '../ttutil/builddir_gfortran/libttutil427_gfortran.a')
    if not fs.exists(ttutil_lib_path)
        error('ttutil gfortran library not found. Build it first: cd ../ttutil && pixi run build-linux-gfortran')
    endif
else
    # Linux Intel compilation
    ttutil_filename = 'libttutil@0@-linux.a'.format(ttutil_version)
    ttutil_lib_path = join_paths(ttutil_lib_dir, 'libttutil.a')
    ttutil_url = 'https://github.com/SWAP-model/ttutil/releases/download/v@0@/@1@'.format(ttutil_version, ttutil_filename)
    # Download ttutil library if not present
    if not fs.exists(ttutil_lib_path)
        run_command('mkdir', '-p', ttutil_lib_dir, check: true)
        run_command('curl', '-L', ttutil_url, '-o', ttutil_lib_path, check: true)
    endif
endif

# Create ttutil dependency
ttutil_dep = declare_dependency(
    link_args: [ttutil_lib_path]
)

# Build fully static executable with downloaded TTUTIL library
if host_machine.system() == 'windows'
    if fc.get_id() == 'gcc'
        # Windows cross-compilation with MinGW
        link_args = ['-static']
        executable_name = 'swap'
    else
        # Windows native compilation with Intel
        # No additional link args needed - static libs already specified above
        link_args = ['/SUBSYSTEM:CONSOLE']
        executable_name = 'swap'
    endif
else
    # Linux native compilation with Intel
    link_args = ['-static']
    executable_name = 'swap'
endif

executable(executable_name, sources,
    dependencies: [ttutil_dep],
    link_args : link_args,
    include_directories: src_inc,
    install : false
)

# Build standalone state module test (no pFUnit required)
state_test_sources = [
    'src/core/swap_log.f90',
    'src/core/swap_state_mod.f90',
    'tests/unit/core/test_state_standalone.f90',
]

executable('test_state', state_test_sources,
    include_directories: src_inc,
    install: false
)

# Build soil state sync test (tests soil_state_t structure)
soil_state_test_sources = [
    'src/core/swap_log.f90',
    'src/core/swap_state_mod.f90',
    'tests/unit/soil/test_soil_state.f90',
]

executable('test_soil_state', soil_state_test_sources,
    include_directories: src_inc,
    install: false
)

# Build drainage state test (tests drainage_state_t and surfacewater_state_t)
drainage_state_test_sources = [
    'src/core/swap_log.f90',
    'src/core/swap_state_mod.f90',
    'tests/unit/drainage/test_drainage_state.f90',
]

executable('test_drainage_state', drainage_state_test_sources,
    include_directories: src_inc,
    install: false
)

# Build solute state test (tests solute_state_t)
solute_state_test_sources = [
    'src/core/swap_log.f90',
    'src/core/swap_state_mod.f90',
    'tests/unit/solute/test_solute_state.f90',
]

executable('test_solute_state', solute_state_test_sources,
    include_directories: src_inc,
    install: false
)

# Build heat state test (tests heat_state_t)
heat_state_test_sources = [
    'src/core/swap_log.f90',
    'src/core/swap_state_mod.f90',
    'tests/unit/heat/test_heat_state.f90',
]

executable('test_heat_state', heat_state_test_sources,
    include_directories: src_inc,
    install: false
)

# Build atmosphere state test (tests atmosphere_state_t)
atmosphere_state_test_sources = [
    'src/core/swap_log.f90',
    'src/core/swap_state_mod.f90',
    'tests/unit/atmosphere/test_atmosphere_state.f90',
]

executable('test_atmosphere_state', atmosphere_state_test_sources,
    include_directories: src_inc,
    install: false
)

# Build macropore state test (tests macropore_state_t)
macropore_state_test_sources = [
    'src/core/swap_log.f90',
    'src/core/swap_state_mod.f90',
    'tests/unit/macropore/test_macropore_state.f90',
]

executable('test_macropore_state', macropore_state_test_sources,
    include_directories: src_inc,
    install: false
)

# Build state sync test (tests SAVE variable migration)
state_sync_test_sources = [
    'src/core/swap_log.f90',
    'src/core/swap_state_mod.f90',
    'src/core/variables.f90',
    'src/core/swap_state_sync.f90',
    'tests/unit/core/test_state_sync.f90',
]

executable('test_state_sync', state_sync_test_sources,
    include_directories: src_inc,
    install: false
)

# Build SWAP BMI shared library for Python (requires swap_bmi.f90)
# Uncomment when BMI interface is implemented:
# bmi_sources = sources + ['src/core/swap_bmi.f90']
# shared_library('swap',
#     bmi_sources,
#     dependencies: [ttutil_dep],
#     include_directories: src_inc,
#     install: true,
#     install_dir: get_option('libdir')
# )

# Optional pFUnit unit tests (off by default)
if get_option('enable_pfunit')
    subdir('tests/unit')
endif