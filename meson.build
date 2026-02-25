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
    'src/atmosphere/atmosphere_state.f90',
    'src/boundary/boundary_state.f90',
    'src/drainage/drainage_state.f90',
    'src/drainage/surfacewater_state.f90',
    'src/heat/heat_state.f90',
    'src/solute/solute_state.f90',
    'src/soil/soil_state.f90',
    'src/core/swap_state_mod.f90',
    'src/core/variables.f90',
    'src/core/swap_state_sync.f90',
    'src/macropore/macropore_state.f90',
    'src/core/initialize.f90',
    'src/core/timecontrol.f90',
    
    # Soil
    'src/soil/soilgrid.f90',
    'src/soil/waterbalance.f90',
    'src/soil/soilhydraulics.f90',
    'src/soil/sptabulated.f90',
    'src/soil/WC_K_models_04_11.f90',
    
    # I/O
    'src/io/swap_csv_output.f90',
    'src/io/readdrainagetoml.f90',
    'src/io/readswaptoml.f90',
    'src/io/readswap.f90',
    'src/io/readmeteo.f90',
    'src/io/swapoutput.f90',
    'src/io/macroporeoutput.f90',
    
    # Atmosphere & Heat
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
    'src/crop/cropgrowth_state.f90',
    
    # Utils
    'src/utils/sharedexchange.f90',
    'src/utils/sharedsimulation.f90',
    'src/utils/arrayutils.f90',
    'src/utils/surfacewaterutils.f90',
    'src/utils/soilhydraulicsutils.f90',
    'src/utils/numericalsolvers.f90',
    
    # Main (last)
    'src/core/swap.f90',
    'src/core/swap_main.f90'
]

# ============================================================================
# External dependency: ttutil
# ============================================================================
if is_gcc
    # Use ttutil as a Meson subproject (built from source with gfortran)
    ttutil_prj = subproject('ttutil', default_options: ['default_library=static'])
    ttutil_dep = ttutil_prj.get_variable('ttutil_dep')
else
    # Use prebuilt library for Intel / Windows
    ttutil_version = '4.2.7'
    ttutil_lib_dir = join_paths(meson.current_build_dir(), 'lib')

    if is_windows and not is_gcc
        ttutil_filename = 'libttutil@0@-windows.a'.format(ttutil_version)
    elif is_windows and is_gcc
        ttutil_filename = 'libttutil@0@-mingw.a'.format(ttutil_version)
    else
        ttutil_filename = 'libttutil@0@-linux.a'.format(ttutil_version)
    endif

    ttutil_lib_path = join_paths(ttutil_lib_dir, 'libttutil.a')
    ttutil_url = 'https://github.com/SWAP-model/ttutil/releases/download/v@0@/@1@'.format(
        ttutil_version, ttutil_filename)

    if not fs.exists(ttutil_lib_path)
        run_command('mkdir', '-p', ttutil_lib_dir, check: true)
        run_command('curl', '-L', ttutil_url, '-o', ttutil_lib_path, check: true)
    endif

    ttutil_dep = declare_dependency(link_args: [ttutil_lib_path])
endif

# ============================================================================
# External dependency: toml-f
# ============================================================================

tomlf_prj = subproject(
    'toml-f',
    version: '>=0.2',
    default_options: [
        'default_library=static',
    ],
)
tomlf_dep = tomlf_prj.get_variable('tomlf_dep')

# ============================================================================
# Main executable
# ============================================================================
if (is_windows and is_gcc) or not is_windows
    link_args = ['-static']
else
    link_args = []
endif

executable('swap', sources,
    dependencies: [ttutil_dep, tomlf_dep],
    link_args: link_args,
    include_directories: src_inc,
    install: false
)
# ============================================================================
# Tests
# ============================================================================
if get_option('enable_unit_tests')
    subdir('tests/unit')
endif

