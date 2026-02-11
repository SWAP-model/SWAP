project('swap', 'fortran', version: '4.2.0', default_options: ['warning_level=0'])

# ============================================================================
# Compiler setup
# ============================================================================
fc = meson.get_compiler('fortran')
fs = import('fs')

src_inc = include_directories(
    'src', 'src/core', 'src/io', 'src/soil', 'src/atmosphere',
    'src/crop', 'src/drainage', 'src/boundary', 'src/macropore',
    'src/solute', 'src/heat', 'src/utils'
)

# Platform/compiler detection
is_windows = host_machine.system() == 'windows'
is_intel = fc.get_id() == 'intel'
is_gcc = fc.get_id() == 'gcc'

# ============================================================================
# Compiler flags
# ============================================================================
if is_intel
    add_project_arguments([
        '-O2', '-fpscomp', 'general', '-warn', 'declarations',
        '-warn', 'unused', '-warn', 'uncalled', '-warn', 'interfaces',
        '-init=zero', '-fpe0', '-fp-model', 'source'
    ], language: 'fortran')
    
    if not is_windows
        add_project_link_arguments(['-static-intel'], language: 'fortran')
    endif

elif is_gcc
    gcc_flags = [
        '-O2', '-ffree-line-length-none', '-Wno-line-truncation',
        '-Wno-compare-reals', '-Wno-missing-include-dirs',
        '-Wno-unused-variable', '-Wno-unused-dummy-argument', '-std=legacy'
    ]
    
    if is_windows
        gcc_flags += ['-mconsole', '-static-libgfortran', '-static-libgcc']
        add_project_link_arguments(['-static', '-mconsole'], language: 'fortran')
    endif
    
    add_project_arguments(gcc_flags, language: 'fortran')
endif

# ============================================================================
# Source files
# ============================================================================
sources = [
    # Core
    'src/core/arrays.f90',
    'src/core/constants.f90',
    'src/core/swap_log.f90',
    'src/core/swap_state_mod.f90',
    'src/core/variables.f90',
    'src/core/swap_state_sync.f90',
    'src/core/initialize.f90',
    'src/core/timecontrol.f90',
    
    # Soil
    'src/soil/sptabulated.f90',
    'src/soil/WC_K_models_04_11.f90',
    'src/soil/tridag.f90',
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
    
    # I/O
    'src/io/swap_csv_output.f90',
    'src/io/readswap.f90',
    'src/io/readmeteo.f90',
    'src/io/swapoutput.f90',
    
    # Atmosphere & Heat
    'src/atmosphere/penmon.f90',
    'src/atmosphere/interception.f90',
    'src/atmosphere/et.f90',
    'src/atmosphere/precipitation.f90',
    'src/atmosphere/meteoday.f90',
    'src/atmosphere/meteodt.f90',
    'src/atmosphere/snow.f90',
    'src/heat/temperature.f90',
    'src/heat/frozencond.f90',
    
    # Boundaries & Drainage
    'src/boundary/boundtop.f90',
    'src/boundary/boundbottom.f90',
    'src/drainage/drainage.f90',
    'src/drainage/divdra.f90',
    'src/drainage/surfacewater.f90',
    
    # Macropore & Solute
    'src/macropore/macropore.f90',
    'src/macropore/macrorate.f90',
    'src/macropore/macroporeoutput.f90',
    'src/solute/solute.f90',
    
    # Crop
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
    
    # Utils
    'src/utils/sharedexchange.f90',
    'src/utils/sharedsimulation.f90',
    'src/utils/arrayutils.f90',
    'src/utils/surfacewaterutils.f90',
    'src/utils/soilhydraulicsutils.f90',
    
    # Main (last)
    'src/core/swap.f90',
    'src/core/swap_main.f90'
]

# ============================================================================
# External dependency: ttutil
# ============================================================================
ttutil_version = '4.2.7'
ttutil_lib_dir = join_paths(meson.current_build_dir(), 'lib')

# Determine platform-specific library
if is_windows and not is_gcc
    ttutil_filename = 'libttutil@0@-windows.a'.format(ttutil_version)
elif is_windows and is_gcc
    ttutil_filename = 'libttutil@0@-mingw.a'.format(ttutil_version)
elif is_gcc
    # Local gfortran build
    ttutil_lib_path = join_paths(meson.current_source_dir(), 
        '../ttutil/builddir_gfortran/libttutil427_gfortran.a')
    if not fs.exists(ttutil_lib_path)
        error('ttutil gfortran library not found. Build it first: cd ../ttutil && pixi run build-linux-gfortran')
    endif
else
    ttutil_filename = 'libttutil@0@-linux.a'.format(ttutil_version)
endif

# Download prebuilt library (if not using local gfortran)
if not (is_gcc and not is_windows)
    ttutil_lib_path = join_paths(ttutil_lib_dir, 'libttutil.a')
    ttutil_url = 'https://github.com/SWAP-model/ttutil/releases/download/v@0@/@1@'.format(
        ttutil_version, ttutil_filename)
    
    if not fs.exists(ttutil_lib_path)
        run_command('mkdir', '-p', ttutil_lib_dir, check: true)
        run_command('curl', '-L', ttutil_url, '-o', ttutil_lib_path, check: true)
    endif
endif

ttutil_dep = declare_dependency(link_args: [ttutil_lib_path])

# ============================================================================
# Main executable
# ============================================================================
if (is_windows and is_gcc) or not is_windows
    link_args = ['-static']
else
    link_args = []
endif

executable('swap', sources,
    dependencies: [ttutil_dep],
    link_args: link_args,
    include_directories: src_inc,
    install: false
)
# ============================================================================
# Tests
# ============================================================================
subdir('tests/unit')

