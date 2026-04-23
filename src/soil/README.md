# src/soil

## Responsibility

Soil water dynamics. Hydraulic property models — Mualem-Van Genuchten
and the PDI / extended variants — live in `WC_K_models_04_11.f90`.
Soil-property tabulation is in `sptabulated.f90` (spline helpers in the
`doln` / `doTSPACK` / `TSPACK` modules). Grid generation and vertical
re-discretization sit in `soilgrid.f90`. The Richards solver, head
calculation, hysteresis, and soil-water state variables are in
`soilhydraulics.f90`. Water-balance integration (water-table, fluxes,
integral checks, water storage) is in `waterbalance.f90`. The
aggregated `soil_state_t` is declared in `soil_state.f90`.

## Public interface

- `soil_state_t`, `soil_state_init`, `soil_state_finalize`,
  `soil_state_reset_cumulative`, `soil_state_reset_intermediate`.
- `functionvalue_04_11` — hydraulic-model evaluator.
- `calcgrid`, `convertdiscrvert` (+ `_state` twins) — grid and
  re-discretization.
- `headcalc`, `soilwater`, `soilwaterstatevar`, `hysteresis`
  (+ `_state` twins for `soilwater`, `soilwaterstatevar`) — Richards
  solver and state-variable update.
- `calcgwl`, `level`, `watertable`, `fluxes`, `integral`,
  `checkmassbal`, `watstor` (+ `_state` twins for the main ones) —
  water-balance integration.
- `TSPBI`, `TSVAL1`, `my_HVAL`, `my_HPVAL` — spline helpers from
  `TSPACK`.

## Dependencies

Imports from:

- `src/core/` — `swap_constants`, `swap_array_dimensions`, `swap_log`,
  `swap_state_mod`, `swap_state_sync`, `variables`.
- `src/utils/` — `array_utils`, `numericalsolvers_mod`,
  `soilhydraulics_utils`.
- `src/boundary/` — `boundbottom_mod`, `boundtop_mod` (boundary fluxes
  enter the Richards solve).
- `src/crop/` — `rootextraction_mod` (root sink terms).
- `src/macropore/` — `macropore_mod` (dual-domain exchange).

Does not depend on `src/atmosphere/`, `src/heat/`, `src/drainage/`,
or `src/solute/` directly.
