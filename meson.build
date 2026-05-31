project('swap', 'fortran', version: '4.2.0', default_options: ['warning_level=0'])

# ============================================================================
# Compiler setup
# ============================================================================
fc = meson.get_compiler('fortran')
fs = import('fs')

src_inc = include_directories(
    'src', 'src/core', 'src/driver', 'src/bindings', 'src/io', 'src/io/csv', 'src/io/toml',
    'src/soilwater', 'src/atmosphere',
    'src/crop', 'src/drainage',
    'src/solute', 'src/heat', 'src/utils', 'src/config', 'src/state'
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
# Source partitions
# ============================================================================
# Modern files: compiled with -std=f2018 -Wall -Wextra.
# Cross-standard linkage works fine — gfortran .mod files produced by
# legacy-mode compilation are consumed without issue by f2018-mode
# compilation. The -std= flag only affects parsing of the file being
# compiled, not what it can `use`.
# Core modern modules shared by every facade. Excludes the C-ABI facade
# modules (swap_capi_mod / swap_bmi_mod / swap_xmi_mod), which export
# OVERLAPPING bind(C) names (initialize/update/finalize/get_current_time/...)
# and therefore must NOT be linked into the same shared library together.
modern_core_sources = [
    'src/driver/swap_mod.f90',
    'src/driver/swap_ensemble_mod.f90',
    'src/bindings/bmi_constants_mod.f90',
    'src/core/timecontrol_mod.f90',
]

# Standalone-BMI facade: singleton (capi_state) backed. Drives libswap_bmi.so
# and the executable. swap_bmi_mod and swap_capi_mod share the C names.
modern_sources = modern_core_sources + [
    'src/bindings/swap_capi_mod.f90',
    'src/bindings/swap_bmi_mod.f90',
]

# XMI coupling facade: ensemble backed. Sole definer of the coupling C-ABI
# (re-implements the BMI lifecycle/time/info C names over the ensemble) so it
# CANNOT coexist with swap_capi_mod/swap_bmi_mod in one library. Drives
# libswap_xmi.so — the library imod_coupler / xmipy loads.
modern_xmi_sources = modern_core_sources + [
    'src/bindings/swap_xmi_mod.f90',
]

# Legacy files: project-wide -std=legacy (set by add_project_arguments).
# Keep every existing entry from the prior `sources` list MINUS the two
# modern files above and MINUS src/core/swap_main.f90 (which becomes the
# executable's own source). swap.f90 and sharedexchange.f90 are already
# deleted in Task 6, so they should not be present here.
sources = [
    # Core
    'src/core/arrays.f90',
    'src/core/constants.f90',
    'src/core/swap_log.f90',
    'src/core/dtutil.f90',
    'src/error/error.f90',
    'src/error/fatalerr.f90',
    'src/validation/validation.f90',
    'src/config/general_config.f90',
    'src/config/simulation_config.f90',
    'src/config/meteorology_config.f90',
    'src/config/drainage_config.f90',
    'src/config/soil_config.f90',
    'src/config/bottom_boundary_config.f90',
    'src/config/heat_config.f90',
    'src/config/irrigation_config.f90',
    'src/config/solute_config.f90',
    'src/config/surface_water_config.f90',
    'src/config/nutrients_config.f90',
    'src/config/cropfixed_config.f90',
    'src/config/cropgrass_config.f90',
    'src/config/cropwofost_config.f90',
    'src/config/crop_config.f90',
    'src/config/output_csv_config.f90',
    'src/crop/cropfixed_init.f90',
    'src/crop/cropwofost_init.f90',
    'src/crop/cropgrass_init.f90',
    'src/config/swap_config.f90',
    # [GR-IO 2026-05-25 Phase 6 Steps 5-6] variables.f90, arrays.fi, and
    #   initialize.f90 all deleted — bare-globals module + its zero-init
    #   shim are fully retired. All defaults live in derived-type component
    #   defaults on state%X / config%X. Size parameters live in
    #   src/core/arrays.f90 (swap_array_dimensions).
    # Soil
    'src/soilwater/waterbalance.f90',
    'src/soilwater/soilhydraulics.f90',
    # [GR-SOIL 2026-05-24] sptabulated moved to src/soil/dormant/sptabulated.f90 (TSPACK + EvalTabulatedFunction)
    # [GR-SOIL 2026-05-27] do_ln_trans.f90 deleted; module doln folded back into dormant sptabulated.f90
    'src/soilwater/WC_K_models_04_11.f90',

    # I/O
    'src/io/csv/csv_common.f90',
    'src/io/csv/meteo_csv.f90',
    'src/io/csv/drainage_csv.f90',
    'src/io/csv/irrigation_csv.f90',
    'src/io/csv/nutrients_csv.f90',
    'src/io/csv/boundary_csv.f90',
    'src/io/csv/soil_init_csv.f90',
    'src/io/csv_reader.f90',
    'src/io/csv_writer.f90',
    'src/io/output_registry.f90',
    'src/io/file_io.f90',
    'src/io/toml/path_helpers.f90',
    'src/io/toml/toml_field_helpers.f90',
    'src/io/toml/toml_array_helpers.f90',
    'src/io/toml/read_general_toml.f90',
    'src/io/toml/read_simulation_toml.f90',
    'src/io/toml/read_meteorology_toml.f90',
    'src/io/toml/read_drainage_toml.f90',
    'src/io/toml/read_bottom_boundary_toml.f90',
    'src/io/toml/read_heat_toml.f90',
    'src/io/toml/read_irrigation_toml.f90',
    'src/io/toml/read_irrigation_ssdi_toml.f90',
    'src/io/toml/read_solute_toml.f90',
    'src/io/toml/read_soil_toml.f90',
    'src/io/toml/read_soil_tillage_toml.f90',
    'src/io/toml/read_surface_water_toml.f90',
    'src/io/toml/read_cropfixed_toml.f90',
    'src/io/toml/read_cropgrass_toml.f90',
    'src/io/toml/read_cropwofost_toml.f90',
    'src/io/toml/read_crop_toml.f90',
    'src/io/toml/load_crop_rotation_files.f90',
    'src/io/toml/read_output_csv_toml.f90',
    'src/io/toml/read_nutrients_toml.f90',
    'src/io/toml/load_swap_config.f90',
    'src/io/toml/load_swap_config_string.f90',
    'src/io/toml/write_swap_config.f90',
    'src/io/csv_aggregates.f90',
    'src/io/csv_output.f90',
    'src/io/meteo_buffer_mod.f90',
    'src/io/readmeteo.f90',
    'src/io/file_headers.f90',

    # Atmosphere & Heat
    'src/atmosphere/atmosphere_constants.f90',
    'src/atmosphere/interception.f90',
    'src/atmosphere/et.f90',
    'src/atmosphere/precipitation.f90',
    'src/atmosphere/runoff.f90',
    'src/atmosphere/meteo_io.f90',
    'src/atmosphere/meteo_orchestrator.f90',
    'src/atmosphere/meteodt.f90',
    'src/atmosphere/snow.f90',
    'src/heat/temperature.f90',
    'src/heat/frozencond.f90',

    # Boundaries & Drainage
    'src/soilwater/boundtop.f90',
    'src/soilwater/boundbottom.f90',
    'src/drainage/drainage.f90',
    'src/drainage/divdra.f90',
    'src/drainage/surfacewater.f90',

    # Solute
    'src/solute/solute.f90',
    # [GR-SOLUTE 2026-05-24] agetracer dropped — body moved to src/solute/dormant/agetracer.f90 (not compiled)

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
    'src/crop/cropgrowth_helpers.f90',
    'src/crop/cropwofost_runtime.f90',
    'src/crop/cropgrass_runtime.f90',
    'src/crop/cropfixed_runtime.f90',
    'src/crop/cropgrowth.f90',

    # State
    'src/state/surfacewater_state.f90',
    'src/state/drainage_state.f90',
    'src/state/solute_state.f90',
    'src/state/heat_state.f90',
    'src/state/hydraulic_params_mod.f90',
    'src/state/soilwater_state.f90',
    'src/state/atmosphere_state.f90',
    'src/state/tillage_state.f90',
    'src/state/timecontrol_state.f90',
    'src/state/mesh_state.f90',
    'src/state/crop_common_state.f90',
    'src/state/crop_fixed_state.f90',
    'src/state/crop_wofost_state.f90',
    'src/state/crop_grass_state.f90',
    'src/state/crop_irrigation_state.f90',
    'src/state/crop_oxygen_state.f90',
    'src/state/crop_state.f90',
    'src/state/nutrients_state.f90',
    'src/state/swap_state.f90',

    # Utils
    # [GR-IO 2026-05-25 Phase 5] sharedsimulation.f90 moved to src/utils/dormant/
    #   — feature was hardcoded-disabled via flSwapShared (.false.) with no
    #   TOML wiring; call sites in swap_mod.f90 deleted with this commit.
    'src/utils/arrayutils.f90',
    'src/utils/surfacewaterutils.f90',
    'src/utils/soilhydraulicsutils.f90',
    'src/utils/numericalsolvers.f90',
]

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
# Static libraries (compile-time grouping)
# ============================================================================
swap_legacy = static_library('swap_legacy',
    sources: sources,
    dependencies: [tomlf_dep],
    include_directories: src_inc,
)

swap_legacy_dep = declare_dependency(
    link_with: swap_legacy,
    include_directories: src_inc,
)

swap_modern = static_library('swap_modern',
    sources: modern_sources,
    dependencies: [tomlf_dep, swap_legacy_dep],
    include_directories: src_inc,
    fortran_args: ['-std=f2018', '-Wall', '-Wextra'],
)

# ============================================================================
# Main executable
# ============================================================================
# Always statically link with gfortran — see docs/adr/0001-gfortran-first.md.
link_args = ['-static']

executable('swap',
    sources: 'src/driver/swap_main.f90',
    link_with: [swap_modern, swap_legacy],
    dependencies: [tomlf_dep],
    fortran_args: ['-std=f2018'],
    include_directories: src_inc,
    link_args: link_args,
    install: false,
)

# ============================================================================
# BMI shared library (Python / imod_coupler entry point)
# ============================================================================
# link_whole forces all object files from the static libs into the .so
# (equivalent to -Wl,--whole-archive), so every BMI symbol is exported.
# Each member library was already compiled with its appropriate -std= flag.
swap_bmi_lib = shared_library('swap_bmi',
    link_whole: [swap_modern, swap_legacy],
    dependencies: [tomlf_dep],
    include_directories: src_inc,
    install: false,
)

# ============================================================================
# XMI shared library (imod_coupler / xmipy coupling entry point)
# ============================================================================
# Ensemble-backed coupling C-ABI. Built from modern_xmi_sources, which swaps
# the singleton BMI/CAPI facades for swap_xmi_mod. swap_xmi is the SOLE definer
# of the overlapping bind(C) names in THIS library, so the link is clean. This
# is the library the coupler loads (libswap_xmi.so), kept separate from
# libswap_bmi.so (standalone singleton BMI + its tests).
swap_modern_xmi = static_library('swap_modern_xmi',
    sources: modern_xmi_sources,
    dependencies: [tomlf_dep, swap_legacy_dep],
    include_directories: src_inc,
    fortran_args: ['-std=f2018', '-Wall', '-Wextra'],
)

swap_xmi_lib = shared_library('swap_xmi',
    link_whole: [swap_modern_xmi, swap_legacy],
    dependencies: [tomlf_dep],
    include_directories: src_inc,
    install: false,
)
# ============================================================================
# Tests
# ============================================================================
if get_option('enable_pfunit')
    subdir('tests/unit')
endif

subdir('tests/bmi')
subdir('tests/cffi-demo')
subdir('tests/coupling')
