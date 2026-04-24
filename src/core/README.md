# src/core

## Responsibility

The framework layer. Owns the program entry points (`swap_main.f90` driver
and `swap.f90` with the `swap(iTask)` dispatcher), the aggregated
`swap_state_t` which composes every domain's state type
(`swap_state_mod.f90`), initialization (`initialize.f90`), time control
(`timecontrol.f90`), a rescue-era legacy bridge that synchronizes the
new types with the old `variables.f90` common-block-style module
(`swap_state_sync.f90`), shared constants (`constants.f90`), and
array-size parameters (`arrays.f90`, `arrays.fi`, `params.fi`,
`description.fi`). Logging is provided by `swap_log.f90`. Every physics
subdir depends on `src/core/`.

## Public interface

- `swap_state_t`, `swap_state_init`, `swap_state_finalize` — top-level
  aggregated model state.
- `io_handles_t`, `time_state_t` plus every domain's `<domain>_state_t`
  re-exported from `swap_state_mod`.
- `swap_state_sync` — bridges `swap_state_t` against `variables.f90`
  during the rescue period (see
  [../../docs/state-management.md](../../docs/state-management.md)).
- `swap_log`: `log_init`, `log_close`, `log_debug`, `log_info`,
  `log_warn`, `log_error`, `log_set_level`, `log_message`, `to_str`.
  See [../../docs/logging.md](../../docs/logging.md) for usage and
  thread-safety constraints.
- `config_to_state` — to be added in Phase 4a Task 25; converts TOML config
  to initial/boundary/parameters state types.
- `swap_constants`, `swap_array_dimensions` — shared parameters.
- `variables` — legacy common-block module preserved for the rescue
  bridge; new code must not add fields to it.

## Dependencies

Imports from every physics subdir in order to compose `swap_state_t`
(`atmosphere_state_mod`, `boundary_state_mod`, `cropgrowth_state_mod`,
`drainage_state_mod`, `heat_state_mod`, `macropore_state_mod`,
`soil_state_mod`, `solute_state_mod`, `surfacewater_state_mod`,
`snow_mod`, `meteodt_mod`, `meteo_mod`, `meteo_process_mod`,
`runoff_mod`, `irrigation_mod`, `tillage_mod` via `management_soil_mod`,
`rootextraction_mod`, `drainage_mod`, `surfacewater_mod`,
`macropore_mod`, `macroporeoutput_mod`, `frozencond_mod`,
`solute_mod`, `readswaptoml_mod`).

This is intentional for the aggregator — `swap_state_t` owns one
component per domain. Outside of the aggregator, core code only uses
`swap_log` and the shared constants.
