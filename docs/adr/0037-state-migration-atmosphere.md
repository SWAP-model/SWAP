---
title: "ADR 0037 — Atmosphere subsystem state-type migration"
date: 2026-05-11
status: accepted
---

# ADR 0037: Atmosphere subsystem state-type migration

**Status:** accepted (Phases 0 + 1 + 2 complete)
**Date:** 2026-05-11
**Migration #:** 7 of N — third coupling-surface arc of the soil-water 4-arc decomposition
**Branch:** `refactor/atmosphere-state`

## Context

Migrations #5 (boundary, ADR 0035) and #6 (crop-uptake, ADR 0036) established the flat-layout precedent for instantaneous scalar fields. Atmosphere is the first arc in the soil-water decomposition where the owned set spans **five reset cadences** (instantaneous, per-day, per-event, intermediate `flzerointr`, cumulative `flzerocumu`), making the ADR 0033 cohort pattern mandatory rather than optional. This arc debuts the cohort pattern in a **freshly-introduced state record** — surfacewater (ADR 0030) retrofitted cohorts onto an existing flat type; here they are introduced at type-creation.

The four soil-water arcs:

- **#5 Boundary** (ADR 0035) — top/bottom coupling-surface fields: `qtop`, `qbot`, 10 others.
- **#6 Crop-uptake** (ADR 0036) — root-sink fields: `qrot`, stress arrays, JvL machinery.
- **#7 Atmosphere** (this arc) — potential-evaporation, precipitation, snow, and interception fields.
- **#8 Soil-water core** — Richards-equation interior: `h`, `theta`, `q`, `kmean`, `pond`, `gwl`, cumulatives.

The primary write sites for `peva`, `ptra`, and `atmdem` live in `meteoday.f90` and `meteodt.f90` (the meteo orchestrator files). Those files were touched additively with dual-writes in Phase 1 without being refactored — a future atmosphere REFACTOR arc handles their structural reorganisation. The arc's entry points were partially state-plumbed from prior arcs (snow.f90 from SS-HEAT; rootextraction/cropgrowth from ADR 0036), providing a partial windfall that offset the larger external read surface.

## Decision

**Migrate 40 atmosphere-owned globals** into a new `atmosphere_state_t` in `src/state/atmosphere_state.f90`, aggregated under `swap_state_t` as `state%atmosphere`. The record introduces **two cohort sub-records** at type-creation per ADR 0033: `atmosphere_intermediate_t` (8 fields) and `atmosphere_cumulative_t` (10 fields), each with a type-bound `reset()`. The remaining 22 fields sit flat on `atmosphere_state_t`. All 40 fields are scalars; no per-node arrays.

Three phases following the standard playbook shape. Phase 0 was a documentation-only audit (0 new config fields needed — `swsublim` and `spev`/`saev` seeds were already covered). Phase 1 introduced the state type, wired `atmosphere_init`, promoted `snow.f90`'s state arg, and installed dual-writes across all 7 home-tree files plus additive patches to `meteoday.f90` and `meteodt.f90`. Phase 2 migrated all external readers, consolidated reset sites, dropped dual-writes, and retired the globals.

### Decisions (D1–D13)

- **D1. Pilot scope — full atmosphere home tree (7 files).** `atmosphere_constants.f90`, `et.f90`, `interception.f90`, `meteoday.f90`, `meteodt.f90`, `precipitation.f90`, `snow.f90`. The discovery's exclusion of `meteoday.f90` and `meteodt.f90` was reversed per design: they are the canonical writers of `peva`, `ptra`, and `atmdem`. Including them as standard dual-write targets eliminates an "excluded file that must be touched" awkwardness. The future refactor arc is free to restructure without invalidating these writes.

- **D2. State-type shape — two cohorts + 22 flat scalars.** `atmosphere_state_t` carries 11 instantaneous scalars (`peva`, `ptra`, `empreva`, `melt`, `subl`, `slw`, `ssnow`, `snowinco`, `graidt`, `nraidt`, `aintcdt`), 9 per-day scalars (`grai`, `nraida`, `atmdem`, `pevaday`, `ptraday`, `gsnow`, `snrai`, `fprecnosnow`, `sicact`), 3 per-event scalars (`ldwet`, `spev`, `saev`), plus `type(atmosphere_intermediate_t) :: intr` (8 fields: `igrai`, `inrai`, `ipeva`, `iptra`, `ievap`, `igsnow`, `isubl`, `isnrai`) and `type(atmosphere_cumulative_t) :: cumu` (10 fields: `cgrai`, `cnrai`, `caintc`, `cpeva`, `cptra`, `cevap`, `cgsnow`, `csubl`, `csnrai`, `cmelt`). ASSOCIATE prefix `at_` used in large compute bodies.

- **D3. `swap_state_t` gains `type(atmosphere_state_t) :: atmosphere`.** `swap_state_mod` adds `use atmosphere_state_mod`. No change to `soilwater_state_t` or `soilwater_init`.

- **D4. `atmosphere_init(state%atmosphere)` wired at `swap.f90:190`.** Zeroes all 22 flat scalars; cohort fields self-initialise to zero via default component initialisation in the type definition. No `numnod`/`nlay` parameters — all atmosphere fields are scalars. Insertion point: after `soilwater_init(state%soilwater, numnod, numlay)`.

- **D5. Cohort reset consolidation — 3 scattered reset blocks collapsed to 2 `reset()` calls.** Legacy resets were spread across `meteoday.f90:412–423` (cgrai/cnrai/caintc + igrai/inrai), `soilhydraulics.f90:1112–1142` (ipeva/iptra/ievap + cpeva/cptra/cevap), and `snow.f90:87–98` (igsnow/isubl/isnrai + cgsnow/csubl/csnrai/cmelt). Consolidated to `call state%atmosphere%intr%reset()` and `call state%atmosphere%cumu%reset()` at `meteoday.f90:ResetMetFlx`. Scattered inline zeroing removed.

- **D6. `igrai`/`inrai` double-reset resolution.** Both `waterbalance.f90:388–389` and `meteoday.f90:414–415` previously zeroed `igrai`/`inrai` under `flzerointr`. With the cohort-owned `reset()`, the duplicate waterbalance reset block is removed. Single canonical reset via cohort `intr%reset()` at `meteoday.f90:ResetMetFlx`.

- **D7. `reduceva` partial migration — 5 atmosphere-owned reads/writes only.** Discovery hazard #6: `reduceva` reads 12 globals; 5 are atmosphere-owned (`empreva`, `ldwet`, `peva`, `spev`, `saev`). Migrated in this arc. The other 7 (`swredu`, `fldaystart`, `cofred`, `dt`, `rsigni`, `nird`, `pond`) stay as legacy reads. New signature: `subroutine reduceva(task, nrai, state)`.

- **D8. `pond` read in `et.f90:reduceva` — deferred.** `reduceva` reads legacy global `pond` at line 637. Kept as legacy global per boundary D5 deferral. Marked with `! [SS-ATM] reads legacy pond — soil-water-core arc migrates`.

- **D9. `snow.f90` signature promotion — `optional intent(in)` → `intent(inout)`.** Snow writes `ssnow`, `peva`, `empreva`, `melt`, `subl` (including zero-out paths). Promotion required for write-back. Callers `swap.f90:216` and `swap.f90:298` unchanged.

- **D10. `nird`/`gird` ownership — deferred.** `interception.f90:DivIntercep` co-writes `nird` by partitioning; canonical owner is `irrigation.f90`. Deferred to future irrigation arc.

- **D11. `precipitation.f90:127` `ssnow = 0.0d0` mutation — kept verbatim.** Becomes `state%atmosphere%ssnow = 0.0_real64`. Existing code comment preserved. No structural fix in this arc.

- **D12. Phase 0 — documentation only.** `swsublim` confirmed covered in `meteorology_config.snow`. `spev`/`saev` initial seeds confirmed in `soil_config.initial`. Zero new config fields; one audit commit.

- **D13. Compile-driven Phase 2.6 — 13 hidden readers across 11 files.** More than crop-uptake's 2, reflecting atmosphere's wider external footprint. Four routines gained new `state` arguments: `DoTillage`, `CNmethod`, `checkmassbal`, `Consolidate_Bdens`. These are natural Phase 2.6 byproducts from atmosphere's broader coupling surface.

### Deferred (out of scope)

- **Future atmosphere REFACTOR arc** — structural reorganisation of `meteoday.f90`/`meteodt.f90` orchestration. State-migration only here; dual-writes land as additive patches.
- **`nird`/`gird` ownership** — irrigation-owned; DivIntercep co-writes `nird` but the canonical owner is `irrigation.f90`. Future irrigation arc.
- **`pond` in `et.f90:reduceva`** — boundary D5 deferral persists; stays legacy global until soil-water-core arc (#8).
- **Soil-water-core cumulatives** (`cqrot`, `iqrot`, `inqrot`, `iqredwet/dry/sol/frs`) — arc #8 territory.
- **`swdrought=2` / `swsublim=1`** — unreached by check-full. Migration is byte-identical-free for those paths (analogous to heat's `swcalt=1` and crop-uptake's `swdrought=2`).

## Consequences

- **40 owned globals retired** from `variables.f90` with `[SS-ATM] retired 2026-05-11` provenance markers. Breakdown: 22 flat scalars (11 instantaneous + 9 per-day + 3 per-event) + 8 intermediate cohort fields + 10 cumulative cohort fields. Largest single-arc retirement to date (heat 11, boundary 12, crop-uptake 23, atmosphere 40).
- **`atmosphere_state_t` with two cohorts** debuts the ADR 0033 pattern at type-creation. surfacewater retrofitted cohorts; atmosphere introduces them in a freshly-designed record.
- **`atmosphere_init(state%atmosphere)` wired at `swap.f90:190`** — scalar-only init, no dimension parameters.
- **~80+ Phase 1 dual-writes** across 7 atmosphere home-tree files (et, interception, precipitation, snow, meteoday, meteodt) plus co-writer touchpoints in waterbalance and soilhydraulics.
- **~155 Phase 2 reads migrated** across 9 reader files (soilhydraulics, waterbalance, cropgrowth, rootextraction, solute, agetracer, swapoutput, swap_csv_output, plus crop subsystem).
- **Cohort reset consolidation** (D5): `meteoday.f90:ResetMetFlx` now invokes `state%atmosphere%intr%reset()` and `state%atmosphere%cumu%reset()`. Three scattered zero-blocks dropped. soilhydraulics atmosphere zero-writes dropped.
- **`igrai`/`inrai` duplicate-reset in `waterbalance.f90`** removed (D6). Single canonical reset via cohort `intr%reset()`.
- **13 compile-driven hidden readers** across 11 files resolved in Phase 2.6. Architectural additions: `DoTillage`, `CNmethod`, `checkmassbal`, `Consolidate_Bdens` gained `state` args.
- **New playbook lesson surfaced (A-2.1):** pre-flight dual-write coverage check catches non-zero write gaps before reader cutover. `nraidt`/`aintcdt` gap caught after regression broke; A-2.2 and A-2.5 applied the lesson. 6 new lessons added to playbook.
- pFUnit: 706 passing, 0 failures, 1 disabled.
- check-full: 5/5 byte-identical at every commit.

### Known residual issues

- **`swdrought=2` / `swsublim=1` paths** have no integration-level regression coverage. Both stub-error or are inactive in TOML regression cases; zero runtime risk.
- **`pond` in `et.f90:reduceva`** stays as legacy global; plumbed in soil-water-core arc (#8).
- **`nird`/`gird`** remain as legacy globals; irrigation arc territory.

## References

- Discovery: `docs/superpowers/specs/2026-05-11-state-migration-atmosphere-discovery.md`
- Design: `docs/superpowers/specs/2026-05-11-state-migration-atmosphere-design.md`
- Plan: `docs/superpowers/plans/2026-05-11-atmosphere-state-migration.md`
- Predecessors: ADR 0030 (surfacewater pilot — cohort retrofit), ADR 0031 (drainage), ADR 0032 (solute), ADR 0033 (cumulative reset cohorts — pattern debuted here at creation), ADR 0034 (heat — flat-layout precedent), ADR 0035 (boundary — first coupling-surface arc), ADR 0036 (crop-uptake — second coupling-surface arc)
- Phase 0: `4731129` (audit)
- Phase 1: `d191f4b` (atmosphere_state_t + cohorts) → `b114d77` (aggregate + init) → `4583d68` (snow intent bump) → `b8cd862` (et dual-write) → `03f6731` (interception) → `69bf0a0` (precipitation) → `bda8918` (snow) → `eb83767` (meteoday) → `99e467e` (meteodt)
- Phase 2: `56a93b0` (soilhydraulics) → `e0080fe` (waterbalance + dedupe) → `3f15b9d` (crop files) → `4b8645d` (solute/drainage/surfacewater) → `1993d45` (output) → `06732c6` (drop dual-write + retire globals)
