# src/crop

## Responsibility

Crop growth and crop-driven processes. Three crop models dispatched from
`cropgrowth.f90`, each in its own sub-package (ADR 0048):

- `fixed/` — simple fixed crop (`cropfixed_*`).
- `grass/` — detailed grass (`cropgrass_*`).
- `wofost/` — full WOFOST, **including all nutrient dynamics**: plant growth
  (`cropwofost_*`), the WOFOST soil sub-model (nitrogen cycling, amendments,
  residues, organic matter, rate constants — the `wofost_soil_*` family),
  `wofostnut.f90`, and the soil-management wrapper `management_soil.f90`.

At the `crop/` root sit the dispatcher (`cropgrowth.f90`,
`cropgrowth_helpers.f90`) and the cross-mode processes: root water uptake with
Feddes-style reduction (`rootextraction.f90`), oxygen-stress response
(`oxygenstress.f90`), irrigation (`irrigation.f90`, incl. SSDI sub-surface
drip), and tillage (`tillage.f90`). The aggregated state is `cropgrowth_state_t`
(includes irrigation, tillage, and WOFOST soil components).

## Public interface

- `CropGrowth_state` — driver for the active crop model.
- `RootExtraction`, `MatricFlux`, `RootExtraction_state` — root water
  uptake.
- `irrigation`, `SSDI_irrigation` (+ `_state` twins) — irrigation.
- `DoTillage`, `DoTillage_state`, `swtill` — tillage.
- `SoilManagement`, `SoilManagement_state` — soil-management wrapper.
- `O2_pars`, `oxygenstress_mod` — oxygen-stress parameters and routines.
- `Wofost_Soil_Declarations`, `Wofost_Soil_Interface` — WOFOST soil
  module entry points.

## Dependencies

Imports from `src/core/` (constants, array dims, array/numerical helpers),
`src/state/` + `src/config/` (the typed records), and `src/soilwater/` (spline
and hydraulics helpers used by the WOFOST soil sub-model). Does not depend on
`src/atmosphere/`, `src/drainage/`, or `src/solute/` directly; inter-domain
flows are communicated through `swap_state_t`.
