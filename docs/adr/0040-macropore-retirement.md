---
title: "ADR 0040 — Macropore subsystem retirement"
date: 2026-05-12
status: accepted
---

# ADR 0040: Macropore subsystem retirement

**Status:** accepted
**Date:** 2026-05-12
**Branch:** `refactor/macropore-retire`

## Context

Macropore flow has been deferred consistently since ADR 0010 (orphan
infrastructure) and ADR 0011 (case 3 excluded from `check-full`) on
2026-04-27. ADR 0038 explicitly deferred `cQMpLatSs` (D11) to "the future
macropore arc." Across every state-rescue arc since — soilwater core
(ADR 0038), atmosphere (ADR 0037), tillage (ADR 0039) and seven earlier
arcs — the same pattern has repeated: macropore globals are touched only
to keep them compiling, never to integrate with the new state types.

Today's macropore footprint on the rescue branch:

- `src/macropore/macropore.f90` — 2245 LoC physics module.
- `src/macropore/macrorate.f90` — 911 LoC rate solver.
- `src/io/macroporeoutput.f90` — 309 LoC output writer.
- `src/config/macropore_config.f90` — orphan typed config (never wired
  to `config_to_variables`).
- `src/io/toml/read_macropore_toml.f90` — never created.
- Roughly 80 macropore-specific globals in `src/core/variables.f90` plus
  zero-initialization scaffolding in `src/core/initialize.f90`.
- `tests/swap-cases/3.macroporeflow` — long-running regression case
  archived but excluded.
- `tests/unit/config/test_macropore_config.pf` — pFUnit suite for the
  orphan typed config.
- Cross-subsystem `flMacroPore`-guarded branches and `FrArMtrx` reads in
  `soilhydraulics`, `waterbalance`, `boundtop`, `boundbottom`, `swap.f90`,
  `timecontrol`, `tillage`, `solute`, `agetracer`, `surfacewater`,
  `swap_csv_output`, `swapoutput`, `soilgrid`, `drainage`.

The cost of keeping the macropore code on the rescue branch (~4500 LoC
of unused physics plus ~80 globals lugged through every arc, plus a
permanently-disabled regression case slot, plus stale typed-config
infrastructure) now exceeds the cost of deleting it and re-implementing
later as a fresh feature.

## Decision

**Retire the entire macropore subsystem from the rescue branch.** Delete
the physics modules, output writer, typed config, regression case, pFUnit
suite, and almost all macropore-specific globals. Force `flMacroPore`
permanently `.false.` so the remaining cross-subsystem `if (flMacroPore)`
branches dead-code naturally. Default the soilwater `FrArMtrx` field to
`1.0` (whole-matrix) since the only writer (`MACROPORE.MACROGEOM`) is
gone. The legacy SWAP 4.2.0 macropore implementation is preserved on
branch `legacy/swap-4.2.0`; the pre-rescue snapshot is preserved on
`archive/main-pre-rescue`. Re-implementation is deferred to a future
"feature" arc that will use the legacy branch as TDD oracle and target
a modern architecture (macropore state type, cohort-aware accumulators,
TOML-native config, no globals).

### Decisions D1–D8

- **D1. Delete macropore physics, output, config, tests, regression case.**
  Removed: `src/macropore/macropore.f90`, `src/macropore/macrorate.f90`,
  `src/macropore/README.md`, `src/io/macroporeoutput.f90`,
  `src/config/macropore_config.f90`, `tests/unit/config/test_macropore_config.pf`,
  `tests/swap-cases/3.macroporeflow/` (was untracked). `meson.build`,
  `tests/unit/meson.build`, `tests/unit/testSuites.inc` updated.

- **D2. Retire ~70 macropore globals from `variables.f90`.** Retired
  set includes all cumulative accumulators (`cQMp*`, `iQMp*`,
  `iQExc*` — except scaffolding kept for soilgrid refinement),
  parameters (`NumSbDm`, `SwDarcy`, `SwDrRap`, `SwPowM`, `SwSorp`,
  `DiPoMa`, `DiPoMi`, `FKcovlay`, `GeomFac`, `PndmxMp`, `PowM`,
  `PpIcSs`, `Rzah`, `ShapeFacMp`, `SorpAlfa`, `SorpMax`, `SorpFacParl`,
  `ShrPar[A-E]`, `Spoint`, `VlMpStSs`, `Z_Ah`, `Z_Ic`, `Z_MB50`, `Z_St`,
  `ZDiPoMa`, `ZnCrAr`, `KsMpSs`, `KsatCovLay`, `PpDmCp`, `PpIcTpMp`),
  and work arrays (`ICpBt*`, `ICpTp*`, `ICpSat*`, `NnCrAr`, `ArMpTpDm`,
  `AwlCorFac`, `FrMpWalWet`, `KDCrRlRef`, `QExcMtxDmCp`,
  `QIn*SatDmCp`, `QInTop*Dm`, `QOut*DmCp`, `SorpDmCp`, `ThtSrpRefDmCp`,
  `TimAbsCumDmCp`, `VlMpDm*`, `WaSrMp*`, `ZBtDm`, `ZWaLevDm`,
  `flBegin`, `flDraTub`, `FlEndSrpEvt`, `NumDm`, `SwBma`, `SubsidCp`,
  `VlMpDyCp`, `VlMpStCp`, `WaLevDm1`, `VlMp`, `VlMpDm1`, `VlMpDm2`,
  `swmacro`, `ICpBtDmPot`, `iQInTopLatDm`, `iQOutDrRapCp` and the
  remainder of the legacy macropore block).

- **D3. Keep ~25 placeholders as "retired-zero" globals.** Cross-subsystem
  code still references them under permanently-dead `if (flMacroPore)`
  branches, dead `FlMacropore .and. Z_Tp.gt.-1.d-8` guards, or as
  zero-valued additive terms. These compile and evaluate to zero:
  `ArMpSs`, `ArMpTp`, `Z_Tp`, `CritUndSatVol`, `cQMpLatSs`,
  `cQMpOutDrRap`, `dFdhMp`, `iQMpOutDrRap`, `iQInTopLatDm1/2`,
  `iQInTopVrtDm1/2`, `IWaSrDm1/2Beg`, `WaSrDm1/2`, `WaSrDm1/2Ini`,
  `DiPoCp`, `IAvFrMpWlWtDm1/2`, `iQExcMtxDm1/2Cp`, `iQOutDrRapCp`,
  `VlMpStDm1/2`, `IcTopMP`, `IDecMpRat`, `QExcMpMtx`, `QMaPo`, `QRapDra`,
  `SwSoilShr`, `ThetCrMp`, `FlDecMpRat`, `flmacropore`. A future
  re-implementation arc will retire these alongside the typed
  `macropore_state_t`. (Note: `NumLevRapDra`, `RapDraReaExp`,
  `RapDraResRef` are drainage-subsystem fields that were merely grouped
  next to macropore in the legacy `variables.f90` — they stay live.)

- **D4. `flMacroPore = .false.` permanent in `initialize.f90`.** No code
  path can set it true; the soil_config validator still rejects
  `soil.swmacro=1` with a stub-error explaining macropore physics is
  retired.

- **D5. `FrArMtrx` defaults to `1.0` in `soilwater_init`** (was `0.0`).
  Matrix-area fraction is unity when macropore is gone; readers like
  `soilhydraulics.headcalc` (multiplications) and `waterbalance.watstor`
  (theta-weighting) collapse cleanly. The only writer
  (`MACROPORE.MACROGEOM`) is deleted with the module.

- **D6. swap.f90 / soilhydraulics.f90 dead-coded.** Five `call MACROPORE(...)`
  sites and three `MacroPoreOutput(...)` sites commented out with
  `[MACRO-RETIRE 2026-05-12]` markers. `boundtop.f90` macropore
  overland-flow branch deleted outright (cleanest path; no live readers).

- **D7. Legacy preserved.** Branches `legacy/swap-4.2.0` (full SWAP 4.2.0
  original macropore code) and `archive/main-pre-rescue` (pre-rescue
  snapshot) remain read-only. Future re-implementation references them
  as TDD oracles.

- **D8. Regression set is now 5 cases (was 5).** `macroporeflow` was
  already excluded from `check-full` since ADR 0011, so no change in
  case count; the case directory and pFUnit parity stubs that referenced
  it are removed.

## Consequences

**Code removed** (rough estimate):

- `src/macropore/macropore.f90` — 2245 LoC
- `src/macropore/macrorate.f90` — 911 LoC
- `src/io/macroporeoutput.f90` — 309 LoC
- `src/config/macropore_config.f90` — ~190 LoC
- `tests/unit/config/test_macropore_config.pf` — orphan suite
- `src/core/variables.f90` macropore block — ~115 lines of declarations
- `src/core/initialize.f90` macropore inits — ~100 lines

Total: ~3870 LoC deleted plus orphan tests, plus archived regression
case directory. Net reduction even after counting added `[retired-zero]`
placeholders and ADR 0040 itself.

**Subsystems unblocked:**

- ADR 0010 (orphan macropore_config_t infrastructure) — closed permanent.
- ADR 0011 (case 3 deferred) — closed permanent.
- ADR 0038 D11 (`cQMpLatSs` deferral) — closed permanent.

**Re-implementation path** (when macropore physics returns as a feature):

1. Branch from `legacy/swap-4.2.0` (or cherry-pick the macropore tree)
   for the TDD oracle.
2. Design a `macropore_state_t` flat type aggregated under
   `swap_state_t` — cumulative cohorts modelled on
   `soilwater_state_t.cumu` if needed.
3. Wire `macropore_config_t` (now deleted; re-create or restore from
   git) directly into `config_to_variables.f90`. No globals.
4. TDD against legacy: run a reference case under `legacy/swap-4.2.0`,
   capture annual stats, build the modern pipeline byte-identical
   against that fixture.
5. Re-add `3.macroporeflow` to the regression matrix.
6. Retire the ~25 retired-zero placeholders from `variables.f90`.

## Verification

- `pixi run check-full`: 5/5 byte-identical (hupselbrook, grassgrowth,
  oxygenstress, salinitystress, surfacewater).
- `pixi run -e test test-pfunit`: 727/0/1 (was 736/0/1; nine macropore
  suites retired with `test_macropore_config.pf` and the macropore
  state-test FrArMtrx assertion adjusted to expect `1.0`).
- Strategy B compile-driven discovery: 7 compile passes — 1 to surface
  the retired-symbol references, 5 to address them in
  `waterbalance`, `soilhydraulics`, `boundtop`, `swap.f90`, and the
  shrinkage/drainage globals; the final link required removing three
  surviving `call MACROPORE(...)` sites and two `use macropore_mod` /
  `use macroporeoutput_mod` clauses.
