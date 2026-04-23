# src/heat

## Responsibility

Soil heat transport and frozen-soil handling. `temperature.f90` solves
the soil temperature profile (tridiagonal advection-diffusion,
conductivity via De Vries mixing). `frozencond.f90` derives frozen-soil
hydraulic conductivity and boundary adjustments when part of the
profile is below 0 °C. State is `heat_state_t` with the standard
lifecycle quartet.

## Public interface

- `heat_state_t`, `heat_state_init`, `heat_state_finalize`.
- `temperature`, `temperature_state`, `devries` — soil-temperature
  solver and conductivity.
- `FrozenCond`, `FrozenBounds` (+ `_state` twins) — frozen-soil
  conductivity and boundary terms.

## Dependencies

Imports from:

- `src/core/` — `swap_array_dimensions`, `swap_state_mod`,
  `swap_state_sync`, `variables`.
- `src/utils/` — `array_utils`, `numericalsolvers_mod`.
- `src/drainage/` — `distribute_drainage` (frozen boundary interacts
  with drainage levels).

Does not depend on `src/soil/` directly; coupling is via
`swap_state_t`.
