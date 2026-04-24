# src/config

## Responsibility

Typed configuration hierarchy (Phase 4a) consumed by the composite
TOML reader, validator stages, and the temporary config-to-state
adapter. One file per section plus an aggregate in `swap_config.f90`.

## Public interface

- `general_config_t` — project id, paths, screen/error switches.
- `simulation_config_t` — dates, output timing.
- `meteorology_config_t` — lat/alt, ET switches, Angstrom, interception.
- `drainage_config_t` — drainage switches, basic/extended tables.
- `soil_config_t` — profile, hydraulic params, initial conditions.
- `crop_config_t` — rotation, per-crop file refs.
- `swap_config_t` — aggregates all of the above.

Each type exposes type-bound procedures `validate(errors)` and
`finalize(errors)`. The aggregate delegates to each section then
runs cross-section invariants.

## Dependencies

`error_mod`, `validation_mod`. Reads no files; produces no
mutation outside its own fields.
