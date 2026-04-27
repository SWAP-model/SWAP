project('swap', 'fortran', version: '4.2.0', default_options: ['warning_level=0'])

# ============================================================================
# Compiler setup
# ============================================================================
fc = meson.get_compiler('fortran')
fs = import('fs')

src_inc = include_directories(
    'src', 'src/core', 'src/io', 'src/io/toml', 'src/soil', 'src/atmosphere',
    'src/crop', 'src/drainage', 'src/boundary', 'src/macropore',
    'src/solute', 'src/heat', 'src/utils', 'src/config'
)

# Platform detection. The rescue spec locks the supported compiler
# to GCC/gfortran (see docs/adr/0001-gfortran-first.md). Reject any
# other compiler up front rather than silently producing a different
# binary.
is_windows = host_machine.system() == 'windows'

if fc.get_id() != 'gcc'
    error('Unsupported Fortran compiler: @0@. Only gfortran is supported during the rescue (see docs/adr/0001-gfortran-first.md).'.format(fc.get_id()))
endif

# ============================================================================
# Compiler flags
# ============================================================================
gfortran_flags = [
    '-O2',
    '-ffree-line-length-none',
    '-Wno-line-truncation',
    '-Wno-compare-reals',
    '-Wno-missing-include-dirs',
    '-Wno-unused-variable',
    '-Wno-unused-dummy-argument',
    '-std=legacy',
    '-finit-local-zero',
]

if get_option('enable_coverage')
    gfortran_flags += ['--coverage', '-O0', '-fprofile-arcs', '-ftest-coverage']
    add_project_link_arguments('--coverage', language: 'fortran')
endif

if is_windows
    gfortran_flags += ['-mconsole', '-static-libgfortran', '-static-libgcc']
    add_project_link_arguments(['-static', '-mconsole'], language: 'fortran')
endif

add_project_arguments(gfortran_flags, language: 'fortran')

# ============================================================================
# Source files
# ============================================================================
sources = [
    # Core
    'src/core/arrays.f90',
    'src/core/constants.f90',
    'src/core/swap_log.f90',
    'src/error/error.f90',
    'src/validation/validation.f90',
    'src/config/general_config.f90',
    'src/config/simulation_config.f90',
    'src/config/meteorology_config.f90',
    'src/config/drainage_config.f90',
    'src/config/soil_config.f90',
    'src/config/bottom_boundary_config.f90',
    'src/config/heat_config.f90',
    'src/config/cropfixed_config.f90',
    'src/config/cropgrass_config.f90',
    'src/config/cropwofost_config.f90',
    'src/config/crop_config.f90',
    'src/config/swap_config.f90',
    'src/core/variables.f90',
    'src/core/initialize.f90',
    'src/core/timecontrol.f90',
    
    # Soil
    'src/soil/soilgrid.f90',
    'src/soil/waterbalance.f90',
    'src/soil/soilhydraulics.f90',
    'src/soil/sptabulated.f90',
    'src/soil/WC_K_models_04_11.f90',
    
    # I/O
    'src/io/toml/path_helpers.f90',
    'src/io/toml/toml_field_helpers.f90',
    'src/io/toml/read_general_toml.f90',
    'src/io/toml/read_simulation_toml.f90',
    'src/io/toml/read_meteorology_toml.f90',
    'src/io/toml/read_drainage_toml.f90',
    'src/io/toml/read_bottom_boundary_toml.f90',
    'src/io/toml/read_heat_toml.f90',
    'src/io/toml/read_soil_toml.f90',
    'src/io/toml/read_cropfixed_toml.f90',
    'src/io/toml/read_cropgrass_toml.f90',
    'src/io/toml/read_cropwofost_toml.f90',
    'src/io/toml/read_crop_toml.f90',
    'src/io/toml/load_swap_config.f90',
    'src/io/toml/write_swap_config.f90',
    'src/io/swap_csv_output.f90',
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
# Vendored as a meson subproject and built from source with gfortran.
# The Intel/ifx download-prebuilt-binary path has been removed; see
# docs/adr/0001-gfortran-first.md.
ttutil_prj = subproject('ttutil', default_options: ['default_library=static'])
ttutil_dep = ttutil_prj.get_variable('ttutil_dep')

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
# Always statically link with gfortran — see docs/adr/0001-gfortran-first.md.
link_args = ['-static']

executable('swap', sources,
    dependencies: [ttutil_dep, tomlf_dep],
    link_args: link_args,
    include_directories: src_inc,
    install: false
)
# ============================================================================
# Tests
# ============================================================================
if get_option('enable_pfunit')
    subdir('tests/unit')
endif

