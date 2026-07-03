# src/crop

## Responsibility

Crop growth and crop-driven processes. Three crop models dispatched from
`cropgrowth.f90`, each in its own sub-package (ADR 0048):

- `fixed/` — simple fixed crop (`cropfixed_*`).
- `grass/` — detailed grass (`cropgrass_*`).
- `wofost/` — WOFOST crop growth (`cropwofost_*`).

> **Nutrients detached (ADR 0052).** The WOFOST-N / ANIMO-derived soil-N
> subsystem (the `wofost_soil_*` family, `wofostnut.f90`, `management_soil.f90`,
> `flCropNut`, and the `[nutrients]`/`[cropwofost.nutrient]` config) was removed
> per the ADR 0051 core-vs-component boundary. Detailed nutrient modelling is
> external ANIMO, fed by SWAP hydrology (the parked `.afo` output, ADR 0009).
> The legacy implementation remains on `legacy/swap-4.2.0`.

At the `crop/` root sit the dispatcher (`cropgrowth.f90`,
`cropgrowth_helpers.f90`) and the cross-mode processes: root water uptake with
Feddes-style reduction (`rootextraction.f90`), oxygen-stress response
(`oxygenstress.f90`), irrigation (`irrigation.f90`, incl. SSDI sub-surface
drip), and tillage (`tillage.f90`). The aggregated state is `cropgrowth_state_t`
(includes irrigation and tillage).

## Public interface

- `CropGrowth_state` — driver for the active crop model.
- `RootExtraction`, `MatricFlux`, `RootExtraction_state` — root water
  uptake.
- `irrigation`, `SSDI_irrigation` (+ `_state` twins) — irrigation.
- `DoTillage`, `DoTillage_state`, `swtill` — tillage.
- `O2_pars`, `oxygenstress_mod` — oxygen-stress parameters and routines.

## Dependencies

Imports from `src/core/` (constants, array dims, array/numerical helpers),
`src/state/` + `src/config/` (the typed records), and `src/soilwater/` (spline
and hydraulics helpers used by the WOFOST soil sub-model). Does not depend on
`src/atmosphere/`, `src/drainage/`, or `src/solute/` directly; inter-domain
flows are communicated through `swap_state_t`.
