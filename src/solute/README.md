# src/solute

## Responsibility

Reactive solute transport: advection-dispersion through the soil
column, with simple first-order reactions and an optional age-tracer
sub-model (`solute.f90`). The aggregated `solute_state_t` with its
lifecycle routines is declared in `solute_state.f90`.

## Public interface

- `solute_state_t`, `solute_state_init`, `solute_state_finalize`,
  `solute_state_reset_cumulative`, `solute_state_reset_intermediate`.
- `solute`, `solute_state` — transport driver.
- `AgeTracer`, `AgeTracer_state` — age-tracer sub-model.

## Dependencies

Imports from:

- `src/core/` — `swap_array_dimensions`, `swap_state_mod`,
  `swap_state_sync`, `Variables`.
- `src/utils/` — `array_utils`.

Reads soil water fluxes and concentrations via `swap_state_t`; does
not import `src/soil/` directly at the rescue baseline.
