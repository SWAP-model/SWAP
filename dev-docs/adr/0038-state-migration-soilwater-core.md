---
title: "ADR 0038 — Soil-water core state-type migration"
date: 2026-05-11
status: accepted
---

# ADR 0038: Soil-water core state-type migration

**Status:** accepted (Phases 0 + 1 + 2 complete)
**Date:** 2026-05-11
**Migration #:** 8 of N — FINAL coupling-surface arc of the soil-water 4-arc decomposition
**Branch:** `refactor/soilwater-core-state`

## Context

The three preceding coupling-surface arcs established the structural foundations this arc consummates. Boundary (ADR 0035) introduced `soilwater_state_t` and wired `soilwater_init`. Crop-uptake (ADR 0036) added 22 per-node arrays and bumped the init signature to `(sw, numnod, nlay)`. Atmosphere (ADR 0037) plumbed `state` into `DoTillage`, `checkmassbal`, `CNmethod`, `Consolidate_Bdens`, and the full atmosphere subsystem. By the time this arc opens, almost every entry point needed for the Richards interior (`soilhydraulics`, `headcalc`, `fluxes`, `integral`, `calcgwl`, `boundtop`, `boundbottom`, `tillage`) already has `state` in scope. The residual plumbing additions were the smallest of any arc despite owning the largest field set.

This arc is architecturally significant for three reasons. First, it **closes the soil-water 4-arc decomposition** — boundary → crop-uptake → atmosphere → core — resolving all open deferrals from prior arcs. Second, it **extends a pre-existing state record**: `soilwater_state_t` already carried 35 fields from ADRs 0035/0036; this arc adds two cohort sub-records into that same record (not into a fresh type), which is a new pattern in the migration series. Third, it validated the largest and most important methodological advance of any arc: **compile-driven retirement (Strategy B)**, described below.

The four soil-water arcs:

- **#5 Boundary** (ADR 0035) — top/bottom coupling-surface fields: `qtop`, `qbot`, 10 others. First population of `soilwater_state_t`.
- **#6 Crop-uptake** (ADR 0036) — root-sink fields, per-node stress arrays, JvL machinery. Per-node arrays debut.
- **#7 Atmosphere** (ADR 0037) — potential-evaporation, precipitation, snow, and interception fields. Cohort pattern debuted in a new state record.
- **#8 Soil-water core** (this arc) — Richards-equation interior: `h`, `theta`, `q`, `kmean`, `pond`, `gwl`, cohort accumulators.

## Decision

**Migrate 74 soil-water-core-owned globals** (the largest single-arc retirement in the migration series) into the existing `soilwater_state_t` in `src/state/soilwater_state.f90` by extending it with 30 flat instantaneous fields, two new cohort sub-records, and a per-day extension of the intermediate cohort. The completed `state%soilwater` record is the canonical soil-water state for all future arcs.

Three phases. Phase 0 was documentation-only (zero new config fields needed). Phase 1 extended the state type and cohort definitions, grew `soilwater_init` to allocate all new per-node allocatables, and installed dual-writes across the home tree and co-writer files. Phase 2 migrated all external readers (~400+ sites across 23 files), performed cohort reset consolidation, resolved prior-arc deferrals, dropped dual-writes, and retired the 74 globals via the compile-driven Strategy B pivot.

### Decisions (D1–D15)

- **D1. Home tree: `soilhydraulics.f90`, `waterbalance.f90`, `soilgrid.f90`, `WC_K_models_04_11.f90`, `soilhydraulicsutils.f90`.** All high-traffic routines already had `state` from prior arcs; residual plumbing additions (`SoilWaterStateVar(task)`, `hysteresis`, `watstor`, `calcgwl` intent-flip) were mechanical.

- **D2. State-type extension — TWO cohorts added to PRE-EXISTING record.** `soilwater_state_t` is EXTENDED (not replaced) with `soilwater_intermediate_t` and `soilwater_cumulative_t`. The flat fields already in the record (boundary 12 + crop-uptake 23) remain. New additions: 30 flat instantaneous fields (per-node arrays `theta`, `h`, `q`, `k`, `kmean`, `dimoca`, `cofgen`, `FrArMtrx`, `fluseksatexm`, plus scalars `pond`, `gwl`, `volact`, `wbalance`, etc.) plus the two cohorts. `soilwater_intermediate_t` carries 22 fields total (14 scalars + 3 arrays + 5 per-day scalars + 2 per-day arrays) with TWO type-bound methods. `soilwater_cumulative_t` carries 14 scalar fields with one type-bound `reset()`.

- **D3. No `swap_state.f90` change.** `soilwater_state_t` is already aggregated in `swap_state_t` as `state%soilwater`; no new component or `use` statement needed.

- **D4. `soilwater_init` signature unchanged.** Stable at `soilwater_init(sw, numnod, nlay)` from crop-uptake C-1.2. Body extended to allocate 24+ new per-node/per-layer allocatables including the cohort arrays (`inq`, `inqrot`, `inqssdi`, `iqdo`, `iqup`, `IThetaBeg`, `qpotrot_day`, `qredtot_day`).

- **D5. Cohort reset consolidation.** Replaced `soilhydraulics.f90:1083–1158` (~75 scattered zero-assignments) with three type-bound procedure calls:
  ```fortran
  if (flzerointr) call state%soilwater%intr%reset()
  if (flDayStart) call state%soilwater%intr%reset_per_day()
  if (flzerocumu) call state%soilwater%cumu%reset()
  ```
  The `pondini = pond` and `volini = volact` rebases immediately after `flzerocumu` remain as inline non-zero rebases (ADR 0033 pattern). ~55 scattered reset lines → 3 call lines.

- **D6. Resolve boundary D5 — `pond` migration.** `pond` → `state%soilwater%pond`. `DoTillage` and `boundtop` were already state-plumbed; 8 pond-write sites in `boundtop.f90` retargeted. Tillage internals `Change_MvGpars`/`Adapt_WC_H` threaded state. `et.f90:reduceva` pond-read marker `[SS-ATM] reads legacy pond` converted to real cutover. Mini-sim writeback retargeted.

- **D7. Resolve boundary D6 — `gwl` migration.** `gwl` → `state%soilwater%gwl`. `calcgwl` intent-flip from `intent(in)` to `intent(inout)`; ~12 write sites retargeted (`gwl`, `nodgwl`, `pegwl`, `bpegwl`, `npegwl`, `gwlflcpzo`, `nodgwlflcpzo`). External readers (drainage, surfacewater, agetracer, et, frozencond, output) standard cutover.

- **D8. Resolve boundary D7 — `kmean` full-array migration.** `kmean(:)` → `state%soilwater%kmean(:)`. All 19+ write sites across `boundtop`, `boundbottom`, and `soilhydraulics` retargeted. All callers state-aware post-arc-7.

- **D9. Resolve heat ADR 0034 `hconduc` `tsoil_node` sentinel.** 14+ call sites that previously passed `0.0d0` or `10.d0` or omitted `tsoil_node` updated to pass `state%heat%tsoil(node)`. Sites: `soilhydraulics.f90` ~10, `boundbottom.f90`, `boundtop.f90`, `rootextraction.f90`, `tillage.f90`, `swapoutput.f90`, `macropore.f90`. All callers had `state` in scope post-arc-7.

- **D10. `swapoutput` mini-sim writeback fully retargeted.** Boundary B-2.6 had retargeted `qbot`; this arc retargets the remaining 4 fields (`gwl`, `pond`, `theta(:)`, `h(:)`) at `swapoutput.f90:3865–3948`. Also: `SoilWaterStateVar(task)` gained a `state` arg (3 call sites).

- **D11. `cQMpLatSs` init-reset DEFERRED.** `soilhydraulics.f90:889` `cQMpLatSs = 0.0d0` stays as-is. Macropore arc territory.

- **D12. `swbotb=-2` runtime mutation KEPT AS LEGACY.** Boundary D12 carry-forward; boundary deviation documented.

- **D13. Grid dimensions stay legacy globals.** `numnod`, `numlay`, `dz`, `z`, `disnod`, `layer`, `botcom`, `nod1lay`, etc. per heat ADR 0034 precedent. Future `grid_t` arc.

- **D14. Multi-owner cumulative: `cgird`/`cnird` in `soilwater_cumulative_t`.** `irrigation.f90:93–95` continues subset-reset until the irrigation arc takes ownership. Same documented multi-owner pattern as atmosphere D10 `nird`/`gird` deferral.

- **D15. Compile-driven retirement (Strategy B) as the Phase 2 global-retirement method.** See dedicated subsection below. 28 iterative compile passes resolved all 74 legacy references.

### Deferred (out of scope)

- **Macropore arc territory** — `cQMpLatSs`, `qimmob`, `QExcMpMtx`, `ArMpSs`/`ArMpTp`, `dFdhMp`, `IcTopMp`, and FrArMtrx co-writes from the macropore-driven path.
- **Future atmosphere REFACTOR arc** — structural reorganisation of `meteoday.f90`/`meteodt.f90`. State-migration already complete for atmosphere fields.
- **Irrigation arc** — `qssdi`/`qssdisum` write-side; `cgird`/`cnird` ownership transfer; `nird`/`gird` global retirement.
- **Grid-dimensions arc** — `numnod`, `numlay`, `dz`, `z`, `disnod`, `layer`, `botcom`, etc.

### Compile-driven retirement (Strategy B)

The original S-2.12 plan and a single-file S-2.12a attempt both broke check-full with a hupselbrook hang (process at 99.9% CPU, no output). Root cause: `soilhydraulicsutils.f90` (containing `hconduc`, `dhconduc`, `watcon`, `moiscap`) and `WC_K_models_04_11.f90` both had `use variables, only: cofgen, fluseksatexm` and are called from **inside the Newton-Raphson iteration** in `headcalc`. Once the legacy half-writes stopped, these routines got stale zero-init values for those arrays → Richards equation diverged.

The fix for `soilhydraulicsutils.f90` and `WC_K_models_04_11.f90` was a **module-level pointer pattern**: `soilhydraulics.f90` calls a small `bind_state_targets` helper at init time that points the utility module's module-level pointers at the `state%soilwater%cofgen` and `state%soilwater%fluseksatexm` arrays. The utilities read through the pointers; no state argument is threaded into their 81+ call sites. Correctness preserved; the pointer bind happens once per run.

The broader **Strategy B methodology** that made this tractable:

1. Comment out all 74 globals in `variables.f90` first (with `! [SS-SWC] retired 2026-05-11` markers).
2. The compiler enumerates every remaining legacy reference as a compile error.
3. Each compile error is a precise instruction: drop the legacy write line, OR migrate a missed read to state, OR drop a use-clause import.
4. Iterate: fix the batch of errors, recompile, fix the next batch.
5. 28 iterative compile passes resolved everything in one commit.

This eliminates the "which-file-first" ordering problem that caused the prior two attempts to introduce NR-coupled stale reads. The compiler becomes the exhaustive inventory tool; there is no manual gap.

**8 routines gained state arguments** during compile-driven discovery — the largest set in any arc: `update_rootdistribution`, `SSDI_irrigation`, `TimeControl`, `handle_exchange`, `outsoilphys`, `capriseoutput`, `outshrinkchar`, `MacroPoreOutput`.

**Transient buffer pattern** for two pre-init config-adapter writes: `config_to_variables` runs before `soilwater_init`, so the state allocatables are not yet allocated. Module-level transient buffers in the adapter hold the values until after init when they are drained into state. Mirrors atmosphere A-2.6's `ssnow`/`ldwet` handling; generalizes to any pre-init writer.

## Consequences

- **74 owned globals retired** from `variables.f90` with `[SS-SWC] retired 2026-05-11` provenance markers. Breakdown: 17 per-node arrays + 18 flat scalars + 22 intermediate cohort fields (14 scalars + 3 arrays + 5 per-day) + 14 cumulative cohort scalars + 8 per-day fields folded into `intr` cohort.
- **`soilwater_state_t` is now the complete soil-water state record.** Total fields across all four arcs: boundary 12 + crop-uptake 23 + atmosphere 40 (in separate `atmosphere_state_t`) + core 74 = 149 fields across 4 arcs retired from `variables.f90`.
- **Two new cohorts added to a PRE-EXISTING state record** — first time in the migration series that cohorts are retrofitted into an existing flat record (contrast: atmosphere debuted cohorts in a fresh type; surfacewater retrofitted via ADR 0033). The retrofit is clean because the flat fields were already state-owned.
- **`soilwater_intermediate_t` has TWO type-bound reset methods:** `reset()` under `flzerointr` gate (all 22 fields) and `reset_per_day()` under `flDayStart` gate (8 per-day subset). Mild ADR 0033 extension — single cohort, two gates, two methods.
- **Prior-arc deferrals all resolved:** `pond` (boundary D5), `gwl` (boundary D6), `kmean` full array (boundary D7), `hconduc tsoil_node` sentinel (heat ADR 0034 residual).
- **Mini-sim writeback fully retargeted.** Boundary B-2.6 covered `qbot`; this arc covers `gwl`/`pond`/`theta`/`h`. `swapoutput.f90:3865–3948` entirely state-side.
- **8 routines gained state args** via compile-driven discovery (largest set in any arc).
- **~400+ read sites migrated** across 23 external reader files.
- **~150+ Phase 1 dual-writes** dropped at Phase 2 completion.
- **Module-level pointer pattern** introduced for `soilhydraulicsutils.f90` and `WC_K_models_04_11.f90` — avoids threading state into 81+ callers of the Newton-Raphson utility functions.
- **pFUnit:** 730 passing, 0 failures, 1 disabled.
- **check-full: 5/5 byte-identical** at every commit. Hupselbrook canary verified: 1.18s (no NR-divergence).

### Known residual issues

- **`cgird`/`cnird` multi-owner** — `irrigation.f90:93–95` continues subset-reset pending irrigation arc. Documented with multi-owner comment in `soilwater_state.f90`.
- **`cQMpLatSs`** — remains at `soilhydraulics.f90:889` as legacy init-zero. Macropore arc territory.
- **`swdrought=2`** (JvL path), **`swbotb=2/4/8`**, **`swsublim=1`** paths have no integration-level regression coverage. All stub-error or are inactive in TOML regression cases; zero runtime risk (documented in predecessor ADRs).

## References

- Predecessors: ADR 0030 (surfacewater — cohort retrofit pilot), ADR 0031 (drainage), ADR 0032 (solute), ADR 0033 (cumulative reset cohorts), ADR 0034 (heat — flat-layout precedent), ADR 0035 (boundary — first coupling-surface arc; D5/D6/D7 deferrals resolved here), ADR 0036 (crop-uptake — per-node arrays), ADR 0037 (atmosphere — third arc; state-arg windfalls)
- Phase 0: documentation-only audit (no new config fields)
- Phase 1: `fd26432` (soilwater_state_t + 2 cohorts) → `8ded56d` (soilwater_init growth) → `ef32b31` (SoilWater(1) init) → `23239a9` (headcalc theta/h) → `63273a5` (headcalc q/k/kmean/dimoca + cohort reset) → `a6c64da` (waterbalance cumu+intr) → `cb35737` (waterbalance fluxes/calcgwl) → `4a59dc4` (boundtop pond+kmean(1)) → `5a63a54` (boundbottom kmean(NN+1)) → `583d88f` (tillage pond/cofgen/theta/h) → `48ece46` (macropore FrArMtrx)
- Phase 2: `2a873ea` (cohort reset consolidation) → `08a0377` (hconduc tsoil_node sentinel) → `24e057e` (soilhydraulics) → `18802fd` (waterbalance) → `571eaea` (boundtop/boundbottom) → `3eaa80e` (tillage) → `ce79e91` (crop subsystem) → `d9b73d2` (solute/drainage/surfacewater) → `38c1d41` (macropore) → `39ee449` (heat/frozencond) → `ed8fd3a` (output + mini-sim) → `ad02a1b` (retire 74 globals — Strategy B)
