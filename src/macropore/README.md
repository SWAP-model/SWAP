# src/macropore

## Responsibility

Dual-domain preferential flow through shrinkage cracks and biopores.
`macropore.f90` advances the macropore water balance and exchange with
the soil matrix (including the shrink sub-process). `macrorate.f90`
computes macropore flow rates. The aggregated `macropore_state_t` and
its lifecycle routines live in `macropore_state.f90`. Macropore-
specific output is written from `src/io/macroporeoutput.f90`.

## Public interface

- `macropore_state_t`, `macropore_state_init`, `macropore_state_finalize`,
  `macropore_state_reset_cumulative`, `macropore_state_reset_intermediate`.
- `macropore`, `macropore_state`, `shrink` — macropore water balance.
- `macrorate` — macropore flow-rate calculation.

## Dependencies

Imports from:

- `src/core/` — `swap_array_dimensions`, `swap_state_mod`,
  `swap_state_sync`, `Variables`.
- `src/utils/` — `soilhydraulics_utils`.
- `src/soil/` — `soilwaterbalance_mod` (mass balance integration).
- `src/drainage/` — `drainage_mod` (drainage fluxes influence macropore
  outflow).

Does not depend on `src/atmosphere/`, `src/crop/`, `src/heat/`, or
`src/solute/` directly.
