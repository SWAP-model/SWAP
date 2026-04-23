# src/drainage

## Responsibility

Lateral drainage and surface-water dynamics. `drainage.f90` computes
drainage fluxes under the standard options (basic Hooghoudt/Ernst via
`bocodrb`, extended table-driven via `bocodre`); `divdra.f90`
distributes a total drainage flux across drainage levels.
`surfacewater.f90` advances the surface-water reservoir coupled to the
drainage fluxes and surface runoff. State types are
`drainage_state_t` and `surfacewater_state_t` with their usual
init / finalize / reset-cumulative / reset-intermediate lifecycle.

## Public interface

- `drainage_state_t`, `drainage_state_init`, `drainage_state_finalize`,
  `drainage_state_reset_cumulative`, `drainage_state_reset_intermediate`.
- `surfacewater_state_t` with the same lifecycle quartet.
- `drainage`, `drainage_state`, `bocodrb`, `bocodre` — drainage fluxes.
- `distribute_drainage` (module providing DIVDRA) — drainage-level
  distribution.
- `SurfaceWater`, `SurfaceWater_state` — surface-water reservoir.

## Dependencies

Imports from:

- `src/core/` — `swap_constants`, `swap_array_dimensions`,
  `swap_state_mod`, `swap_state_sync`, `variables`.
- `src/utils/` — `array_utils`, `surfacewater_utils`.

Does not depend on other physics subdirs.
