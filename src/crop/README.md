# src/crop

## Responsibility

Crop growth and crop-driven processes. Three crop models dispatched from
`cropgrowth.f90`: simple fixed crop, detailed grass, and full WOFOST.
Root water uptake with Feddes-style reduction lives in
`rootextraction.f90`; oxygen-stress response in `oxygenstress.f90`.
Management processes include irrigation (`irrigation.f90`, incl. SSDI
sub-surface drip), tillage (`tillage.f90`), and soil management wrapper
(`management_soil.f90`). The WOFOST soil sub-model — nitrogen cycling,
amendments, residues, organic matter, rate constants — is carried in
the `wofost_soil_*` family of files. The aggregated state is
`cropgrowth_state_t` (includes irrigation, tillage, and WOFOST soil
components).

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

Imports from:

- `src/core/` — `swap_constants`, `swap_array_dimensions`,
  `swap_state_mod`, `swap_state_sync`, `variables`.
- `src/utils/` — `array_utils`, `soilhydraulics_utils`.
- `src/soil/` — `doln` (spline helpers used by WOFOST soil).

Does not depend on `src/atmosphere/`, `src/drainage/`, `src/macropore/`,
or `src/solute/` directly; inter-domain flows are communicated through
`swap_state_t`.
