# src/soil

## Responsibility

Soil water dynamics: the Richards-equation solver, the soil hydraulic
property functions, and water-balance integration.

Three compiled modules:

- **`soilhydraulics.f90`** (`soilhydraulics_mod`) — the Richards solver
  and the soil-water lifecycle. `headcalc` solves the mixed-form Richards
  equation by Newton–Raphson (its per-iteration residual is built by the
  private `headcalc_residual`). `soilwater_seed` initialises the profile
  (decomposed into the narrative helpers `init_soil_misc`,
  `populate_hydraulic_params`, `apply_soil_initial_conditions`,
  `compute_initial_node_hydraulics`); `soilwater_step` advances one step;
  `soilwater_update` recomputes conductivities, storage, fluxes, and
  cumulative balances; `soilwater_save_state`/`soilwater_restore_state`
  snapshot/restore for adaptive time-stepping; `hysteresis` updates the
  van Genuchten scanning curves.

- **`WC_K_models_04_11.f90`** (`WC_K_models_04_11`) — the soil hydraulic
  property library. Given a pressure head `h`, `functionvalue_04_11`
  returns water content θ(h) (`iType=1`), hydraulic conductivity K(h)
  (`iType=2`), or differential moisture capacity C(h)=dθ/dh (`iType=3`),
  for the Mualem–van Genuchten family and the Peters–Durner–Iden (PDI)
  extension, uni- and bi-modal, with optional air-entry / adsorptive /
  film+vapour terms (models 4–11). Called per node, per Newton iteration
  (via `watcon`/`hconduc`/`moiscap` in `soilhydraulics_utils`). **The
  module holds no mutable state — it is reentrant**: all curve parameters
  arrive in a `vanGenuchten_params_t` argument and all scratch is local.

- **`waterbalance.f90`** (`soilwaterbalance_mod`) — water-balance
  integration: `calcgwl`/`level` (groundwater + perched water table from
  the head profile), `fluxes` (inter-compartment fluxes from volume
  change + sinks/boundaries), `integral` (intermediate + cumulative flux
  accumulation and the water-balance error), `watstor` (profile storage).

The vertical grid is no longer generated here — `calcgrid` was folded
into `state%mesh%init` (`src/state/mesh_state.f90`). The aggregated
soil-water state record `soilwater_state_t` lives in
`src/state/soilwater_state.f90`; the hydraulic-parameter type
`vanGenuchten_params_t` in `src/state/hydraulic_params_mod.f90`.

## Public interface

- `soilhydraulics_mod`: `headcalc`, `soilwater_seed`, `soilwater_step`,
  `soilwater_update`, `soilwater_save_state`, `soilwater_restore_state`,
  `hysteresis`.
- `WC_K_models_04_11`: `functionvalue_04_11` (sole public symbol).
- `soilwaterbalance_mod`: `calcgwl`, `level`, `fluxes`, `integral`,
  `watstor`.

## Dormant (not compiled)

`dormant/` holds code preserved for possible reactivation but excluded
from the build (no live callers in the TOML pipeline):

- `sptabulated.f90` — tabulated soil-physics path (TSPACK spline library
  + `EvalTabulatedFunction`/`PreProcTabulatedFunction`); also now hosts
  the `doln` module (`do_ln_trans` log-transform flag).
- `watertable.f90` — `CritUndSatVol` perched-table search.
- `checkmassbal.f90` — ANIMO/PEARL mass-balance coupling dispatch.
- `regrid.f90` — `convertdiscrvert` vertical re-discretisation.

Each carries a header with its reactivation checklist.

## Dependencies

Imports from:

- `src/state/` — `swap_state_mod` (the `soilwater`/`mesh`/etc.
  sub-records are reached through it) and `hydraulic_params_mod`
  (`vanGenuchten_params_t`).
- `src/core/` — `swap_constants`, `swap_array_dimensions`, `swap_log`,
  `error_mod`.
- `src/utils/` — `array_utils`, `numericalsolvers_mod`,
  `soilhydraulics_utils`.
- `src/boundary/` — `boundbottom_mod`, `boundtop_mod` (boundary fluxes
  enter the Richards solve).
- `src/crop/` — `rootextraction_mod` (root sink terms).

No dependency on the retired bare-globals `variables` module, on
`src/macropore/` (retired, ADR 0040), or directly on
`src/atmosphere/`, `src/heat/`, `src/drainage/`, `src/solute/`.
