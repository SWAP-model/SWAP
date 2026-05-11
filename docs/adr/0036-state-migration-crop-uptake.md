---
title: "ADR 0036 — Crop water uptake state-type migration"
date: 2026-05-11
status: accepted
---

# ADR 0036: Crop water uptake state-type migration

**Status:** accepted (Phases 0 + 1 + 2 complete)
**Date:** 2026-05-11
**Migration #:** 6 of N — second coupling-surface arc of the soil-water 4-arc decomposition
**Branch:** `refactor/crop-uptake-state`

## Context

Migration #5 (boundary, ADR 0035) carved the top/bottom coupling-surface scalars into `state%soilwater` as its first 12 fields. It anticipated this arc explicitly: `soilwater_init` was placed at `swap.f90:183` (after `CalcGrid()`, before `DoTillage(1)`) precisely to allow later arcs to add per-node arrays without an init-order gap. This arc is that anticipated next step.

The four soil-water arcs decompose by coupling surface rather than by implementation phase:

- **#5 Boundary** (ADR 0035) — top/bottom coupling fields: `qtop`, `qbot`, 10 others.
- **#6 Crop-uptake** (this arc) — root-sink fields: `qrot`, stress components, JvL microscopic machinery.
- **#7 Atmosphere** — potential-evaporation, precipitation, and snow fields.
- **#8 Soil-water core** — Richards-equation interior: `h`, `theta`, `q`, `kmean`, `pond`, `gwl`, cumulatives.

`RootExtraction`, `JongvanLier`, `JongvanLierLoop`, `MatricFlux`, and `OxygenStress` already carried `state intent(in)` from the heat arc (SS-HEAT Task 6) and the solute arc (SS-SLST), providing a **state-arg windfall** symmetric with boundary's entry-point inheritance from SS-HEAT Task 9. The arc is field-carve + co-writer dual-writes + reader cutover — no signature surgery.

## Decision

**Migrate 23 crop-uptake-owned globals** into a flat extension of the existing `soilwater_state_t` in `src/state/soilwater_state.f90`. This is the first population of `state%soilwater` with per-node allocatable arrays.

Three phases matching the standard playbook shape. Phase 0 was a documentation-only coverage audit (zero new typed-config fields needed). Phase 1 extended the state type and installed dual-writes. Phase 2 migrated external readers, dropped dual-writes, and retired the globals.

### Decisions (D1–D12)

- **D1. Pilot scope** — Home file `rootextraction.f90`; co-writer `cropgrowth.f90` Task-1 init paths (CropFixed/Wofost/Grass); external readers: `soilhydraulics.f90` (10 sites), `waterbalance.f90` (~12), `swapoutput.f90` (~14), `solute.f90` + `agetracer.f90` (3 each), `cropgrowth.f90` (5 `flWrtNonox` compute sites). `swap_csv_output.f90` reads only downstream cumulatives — migrated only if compile-driven.
- **D2. State-type shape** — Flat extension of `soilwater_state_t`; no new cohort sub-records. All 23 owned fields are instantaneous or init-once (`mfluxtable`). Matches ADR 0034 (heat) + ADR 0035 (boundary) flat-layout precedent.
- **D3. No `swap_state.f90` change** — `soilwater_state_t` extends in place; `swap_state_t` already holds `type(soilwater_state_t) :: soilwater` (ADR 0035). Zero structural change to `swap_state.f90`.
- **D4. `soilwater_init` signature bump to `(sw, numnod, nlay)`** — Allocates 12 per-node arrays to `numnod` and `mfluxtable(nlay,801)`. Zeroes all scalars. Call site at `swap.f90:184` updated to `call soilwater_init(state%soilwater, numnod, numlay)`. Both dimensions are available immediately after `CalcGrid()`. **First time `state%soilwater` allocates per-node arrays.** Future arcs (atmosphere, soil-water-core) will use this signature unchanged.
- **D5. `CropGrowth(1, state)` plumbing** — Signature extended from `(task, tsoil)` to `(task, tsoil, state)`. All four call sites in `swap.f90` pass state. Co-writer dual-writes for `hroot`, `hleaf`, and `mfluxtable` installed in the CropGrowth dispatcher after the three Task-1 init dispatches; the CropFixed/Wofost/Grass init paths themselves were not state-threaded — minimal scope.
- **D6. `mfluxtable` init — design deviation** — Full relocation to `soilwater_init` was blocked: `ksatfit`, `wiltpoint`, and `cofgen` (per-layer soil hydraulic params needed to build the lookup) are populated by `SoilHydraulics(1)` inside `SoilWater(1, state)` at `swap.f90:190`, which runs *after* `soilwater_init` at `swap.f90:184`. Fallback: allocate `mfluxtable(nlay,801)` in `soilwater_init`; build via `MatricFlux(task=1)` called from the `CropGrowth(1)` dispatcher — the single Task-1 caller that has state. The build itself moved up from the three CropFixed/Wofost/Grass init paths to the dispatcher; initialization sequence preserved.
- **D7. JvL scalars/arrays on state** — `Tactual`, `alpJvLier`, `rmax(:)` have zero external readers; kept on `state%soilwater` for symmetry with the rest of the JvL field group. ASSOCIATE prefix `cw_` (crop-water) used inside `rootextraction.f90` and `cropgrowth.f90` write paths.
- **D8. `flWrtNonox` in `state%soilwater`** — Sole writer is `RootExtraction`; five external readers in `cropgrowth.f90` gate `rr=0`/`grrt=0` in compute task=3 blocks. Consistent with the owned-and-touched rule from boundary. Semantically more crop than soil-water; a future crop-state arc may relocate it. Documented here for that arc.
- **D9. Macropore co-write — no action** — Confirmed: macropore reads downstream cumulative `inqrot` (out-of-scope) but does not write any of the 23 owned fields.
- **D10. Cumulatives out of scope** — `cqrot`, `iqrot`, `inqrot`, `iqredwet/dry/sol/frs`, `qpotrot_day`, `qredtot_day`, `iptra_day` owned by `waterbalance.f90` / `soilhydraulics.f90`; assigned to the soil-water-core arc (#8). After Phase 2 cutover, `waterbalance.f90:418–435` reads `state%soilwater%X` for instantaneous values while legacy cumulative writes remain.
- **D11. Phase 0 coverage audit** — Grep of five TOML regression cases confirmed: `swdrought=2` (JvL) is **not covered** in any case; `swoxygen=1` (Feddes oxygen stress) is covered; `swfrost=1` is covered. JvL path (lines 313–760 of `rootextraction.f90`) has no integration-level regression coverage. The 14 JvL fields are also gated by `swdrought=2` which stub-errors at parse time — their migration is code-review hygiene, not regression risk (see D12). No new fixtures added in this arc; gap documented.
- **D12. Compile-driven Phase 2.7** — Expected 2–4 hidden readers; **found 5 hidden readers in 1 file**: `flWrtNonox` read at 5 sites in `cropgrowth.f90` (CropFixed task=3, CropWofost task=3/3, CropGrass task=3/3) — all in large compute subroutines not visible in the pre-flight grep. Migrated via optional `state` arg threaded from the `CropGrowth` dispatcher to the variant subroutines; bare-name fallback handles edge cases. Only 2 compile iterations required (vs boundary's 8 iterations) — the pattern is maturing.

### Deferred (out of scope)

- **`cqrot`, `iqrot`, `inqrot`, `iqredwet/dry/sol/frs`** and related day/period sums — cumulative cohort; owned by `waterbalance.f90` and assigned to the soil-water-core arc (#8).
- **`pond` / `gwl` / `kmean`** — deferred per boundary D5/D6/D7; untouched here.
- **`hconduc` / `tsoil_node` sentinel residual** from SS-HEAT — soil-water-core arc territory.
- **Three duplicated JvL init blocks** in CropFixed/Wofost/Grass — refactoring to a shared helper is out of scope; noted for a future cropgrowth refactor arc.
- **`swdrought=2` regression fixture** — no new test cases added in this arc; coverage gap documented.

## Consequences

- **23 owned globals retired** from `variables.f90` with `[SS-CRP] retired 2026-05-11` provenance markers: 12 per-node arrays (`qrot`, `qpotrot`, `qredwet`, `qreddry`, `qredsol`, `qredfrs`, `mflux`, `mroot`, `hroot`, `rootrho`, `rootphi`, `rmax`), 1 two-dimensional table (`mfluxtable(nlay,801)`), 5 primary scalars (`qrosum`, `qredwetsum`, `qreddrysum`, `qredsolsum`, `qredfrssum`), 4 JvL scalars (`Tactual`, `alpJvLier`, `hleaf`, `Hxylem`), 1 flag (`flWrtNonox`).
- **`state%soilwater` now carries 35 fields** (12 boundary + 23 crop-uptake). First arc to allocate per-node arrays in the record.
- **`soilwater_init(sw, numnod, nlay)` signature** is stable for the two remaining arcs (atmosphere, soil-water-core); no further signature changes required.
- **`CropGrowth(1, state)` plumbed** — future crop-state arcs inherit this threading without additional surgery.
- **~131 dual-write legacy writes dropped** at Phase 2.7; 5 compile-surfaced hidden readers in `cropgrowth.f90` resolved via optional state arg.
- **~47 external read sites migrated** across 7 files (soilhydraulics, waterbalance, swapoutput, solute, agetracer, cropgrowth).
- **JvL path (`swdrought=2`) stub-errored at parse time** — 14 of 23 fields are unreachable at runtime through TOML. Migration is byte-identical-free by construction; risk is code-review hygiene only.
- pFUnit: 695 passing, 0 failures, 1 disabled throughout (4 new tests for per-node allocation lifecycle).
- check-full: 5/5 byte-identical at every commit.

### Known residual issues

- **`swdrought=2` JvL path** has no integration-level regression coverage. The path stub-errors at TOML parse time, so regression risk is zero; but the large JvL compute block (lines 313–760 of `rootextraction.f90`) is untested beyond unit tests.
- **`flWrtNonox` semantic placement** — currently in `state%soilwater` for pragmatic ownership reasons; a future crop-state arc should relocate it to `state%crop` or equivalent.
- **`pond` / `gwl` / `kmean`** remain as legacy globals; plumbed in soil-water-core arc (#8).

## References

- Discovery: `docs/superpowers/specs/2026-05-10-state-migration-crop-uptake-discovery.md`
- Design: `docs/superpowers/specs/2026-05-10-state-migration-crop-uptake-design.md`
- Plan: `docs/superpowers/plans/2026-05-10-crop-uptake-state-migration.md`
- Predecessors: ADR 0030 (surfacewater pilot), ADR 0031 (drainage), ADR 0032 (solute), ADR 0033 (cumulative reset cohorts), ADR 0034 (heat — flat-layout precedent), ADR 0035 (boundary — first coupling-surface arc + soilwater_state_t definition)
- Phase 0: `fd5c39d` (coverage audit) → `6c047a5` (design + plan)
- Phase 1: `0f0dc44` (soilwater_state_t extension) → `da7cf01` (soilwater_init signature bump) → `1c89ba3` (CropGrowth plumbing) → `9d5a2dd` (rootextraction dual-write)
- Phase 2: `af464b5` (soilhydraulics) → `10eb214` (waterbalance) → `73cdc79` (solute/agetracer) → `f4e1676` (output) → `6893d88` (drop dual-write + retire globals)
