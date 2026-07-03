project('swap', 'fortran', version: '4.2.0', default_options: ['warning_level=0'])

# ============================================================================
# Compiler setup
# ============================================================================
fc = meson.get_compiler('fortran')
fs = import('fs')

src_inc = include_directories(
    'src', 'src/core', 'src/timecontrol', 'src/driver', 'src/bindings', 'src/io', 'src/io/csv', 'src/io/toml',
    'src/soilwater', 'src/atmosphere',
    'src/crop', 'src/crop/fixed', 'src/crop/grass', 'src/crop/wofost', 'src/drainage',
    'src/solute', 'src/heat', 'src/config', 'src/state'
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
# Modern modules — one list, one facade set (ADR 0050 sub-arc 2). The former
# BMI-vs-XMI bind(C) name collision is resolved: swap_bmi_mod is the sole,
# mode-aware definer of the shared lifecycle/metadata C names (single-column
# vs coupled-ensemble, both backed by swap_ensemble_mod), and swap_xmi_mod
# keeps only the XMI-specific verbs (prepare/solve/get_value_ptr/...). All
# facades therefore link into ONE shared library: libswap.so.
modern_sources = [
    'src/driver/swap_mod.f90',
    'src/driver/swap_ensemble_mod.f90',
    'src/bindings/bmi_constants_mod.f90',
    'src/bindings/swap_c_strings_mod.f90',
    'src/bindings/swap_var_registry_mod.f90',
    'src/timecontrol/timecontrol_mod.f90',
    'src/bindings/swap_capi_mod.f90',
    'src/bindings/swap_bmi_mod.f90',
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
    'src/error/diagnostics.f90',
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
    'src/config/cropfixed_config.f90',
    'src/config/cropgrass_config.f90',
    'src/config/cropwofost_config.f90',
    'src/config/crop_config.f90',
    'src/config/output_csv_config.f90',
    'src/crop/fixed/cropfixed_init.f90',
    'src/crop/wofost/cropwofost_init.f90',
    'src/crop/grass/cropgrass_init.f90',
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
    'src/io/csv/boundary_csv.f90',
    'src/io/csv/soil_init_csv.f90',
    'src/io/csv_reader.f90',
    'src/io/csv_writer.f90',
    'src/io/output_registry.f90',
    'src/io/file_io.f90',
    'src/io/toml/path_helpers.f90',
    'src/io/toml/config_source.f90',
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
    'src/crop/oxygenstress.f90',
    'src/crop/rootextraction.f90',
    'src/crop/irrigation.f90',
    'src/crop/tillage.f90',
    'src/crop/cropgrowth_helpers.f90',
    'src/crop/wofost/cropwofost_runtime.f90',
    'src/crop/grass/cropgrass_runtime.f90',
    'src/crop/fixed/cropfixed_runtime.f90',
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
    'src/state/swap_state.f90',
    'src/state/results_state.f90',

    # Utils
    # [GR-IO 2026-05-25 Phase 5] sharedsimulation.f90 moved to src/utils/dormant/
    #   — feature was hardcoded-disabled via flSwapShared (.false.) with no
    #   TOML wiring; call sites in swap_mod.f90 deleted with this commit.
    'src/core/arrayutils.f90',
    'src/drainage/surfacewaterutils.f90',
    'src/soilwater/soilhydraulicsutils.f90',
    'src/core/numericalsolvers.f90',
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
# The one shared library: libswap.so (BMI + CAPI + XMI — ADR 0050)
# ============================================================================
# link_whole forces all object files from the static libs into the .so
# (equivalent to -Wl,--whole-archive), so every C-ABI symbol is exported.
# Every consumer (cffi demos, ctypes orchestrator, xmipy/imod_coupler) loads
# this library. Target name 'libswap' + empty prefix yields libswap.so while
# staying distinct from the 'swap' executable target.
swap_lib = shared_library('libswap',
    link_whole: [swap_modern, swap_legacy],
    dependencies: [tomlf_dep],
    include_directories: src_inc,
    name_prefix: '',
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
