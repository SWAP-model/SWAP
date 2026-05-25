---
title: "ADR 0039 — Tillage state-type migration"
date: 2026-05-12
status: accepted
---

# ADR 0039: Tillage state-type migration

**Status:** accepted (all tasks complete)
**Date:** 2026-05-12
**Migration #:** 9 of N — smallest arc to date; first arc after the soil-water 4-arc decomposition
**Branch:** `refactor/tillage-state`

## Context

After the soil-water 4-arc decomposition closed (boundary ADR 0035, crop-uptake ADR 0036,
atmosphere ADR 0037, soil-water core ADR 0038), tillage is the natural next finishing arc.
The subsystem is unusually self-contained: single home file (446 LoC), zero external compute
readers of `till_*` fields, and all entry points (`DoTillage`, `Change_MvGpars`, `Adapt_WC_H`,
`Consolidate_Bdens`) already carry `state intent(inout)` as a complete windfall from prior
arcs. Phase 0 is zero — ADR 0021 completed the TOML port; typed `soil_tillage_t` already
covers every input parameter.

Discovery revealed that the 31 `till_*` globals divide into two distinct categories: 18
config-constants (Groups A+B) that are immutable after init and already represented in typed
`soil_tillage_t`, and 13 genuine runtime-state fields (Groups C+D+E) that change during the
simulation. Only the 13 runtime-state fields migrate in this arc.

## Decision

**Migrate 13 runtime-state `till_*` globals** (Groups C+D+E) into a new flat `tillage_state_t`,
aggregated as `state%tillage` under `swap_state_t`. Groups A+B (18 config-constant fields
already in typed `soil_tillage_t`) are deferred to a future config-consolidation arc. The
H-4 latent bug in `set_iTill` is fixed in-arc. Bdens and ParamVG cross-subsystem coupling
stays as legacy globals.

### Decisions D1–D10

- **D1. Home tree: `src/crop/tillage.f90` only.** Single file, 446 LoC. Co-writer:
  `config_to_variables.f90` (Groups A+B legacy writes unchanged in this arc).

- **D2. State-type shape — flat, no cohorts.** All 13 migrated fields are instantaneous:
  no `flzerointr`/`flzerocumu` cadence anywhere in tillage. Matches heat (ADR 0034) and
  boundary (ADR 0035) flat-layout precedent.

- **D3. Aggregation.** `swap_state.f90` adds `type(tillage_state_t) :: tillage` alongside
  existing soilwater, atmosphere, heat, etc.

- **D4. Scope: 13 runtime-state fields.**
  - **Group C — 7 per-layer allocatables** (sized to `numlay`): `till_Rho_tillage`,
    `till_Rho_cons`, `till_Rho_last`, `till_K_R_cons`, `till_Rho_match`, `till_N_match`,
    `till_Slope_match`.
  - **Group D — 3 per-step scalars:** `till_sumDWC`, `till_sumAvail1`, `till_sumAvail2`.
  - **Group E — 3 init-once geometry/cursor fields:** `till_MaxNumSoilHo`,
    `till_MaxNumSoilCP`, `till_iTill`.

- **D5. Groups A+B (18 config-constants) deferred.** Already represented in typed
  `soil_tillage_t`. Migrating them to `state%tillage` would duplicate typed config. A
  future config-consolidation arc retires them.

- **D6. Bdens/ParamVG NOT in scope.** Tillage writes `Bdens` and `ParamVG` as coupling
  side-effects; `soilhydraulics.f90`, `solute.f90`, and `oxygenstress.f90` read them.
  Cross-subsystem ownership clarification is its own arc. Stay as legacy globals.

- **D7. `tillage_init(state, numlay)` placement.** New subroutine in
  `src/state/tillage_state.f90`. Call inserted at `swap.f90:228`, between
  `atmosphere_init` (line 214) and `DoTillage(1)` (line 227+). Allocates 7 per-layer
  Group C arrays to `numlay`. Group E geometry scalars (`MaxNumSoilHo`, `MaxNumSoilCP`)
  computed by `det_MNSH` inside `DoTillage(1)` — assigned into `state%tillage` there.

- **D8. Pre-init pattern for Group E cursor.** `DoTillage(1)` seeds `state%tillage%iTill`,
  `MaxNumSoilHo`, `MaxNumSoilCP` after `tillage_init` has run (allocation-safe). Groups
  A+B remain as legacy globals written by `config_to_variables.f90` in the pre-init window;
  no transient buffer needed since they stay legacy throughout this arc.

- **D9. H-4 latent bug fix.** `set_iTill` (tillage.f90:405–414) contained a tautological
  condition: `t1900 < Date_tillage(i-1)` on both sides of the `.and.`, so the branch is
  always false and `iTill` never advances past event 1. Fixed in T-3: second index changed
  from `i-1` to `i`. Observationally a no-op in the regression suite (all 5 cases have
  `swtill=0`, so `flTillage=.false.` and this code never executes), but corrects multi-event
  tillage behavior for any user with `swtill=1` and more than one scheduled event.

- **D10. Compile-driven retirement (Strategy B).** Applied from ADR 0038 playbook lesson #1:
  comment out 13 `till_*` globals in `variables.f90`; compiler enumerates all remaining
  legacy references as compile errors; iterate fix-recompile. Only 2 compile passes were
  needed — smallest ever (boundary 8, atmosphere 13, soil-water 28). T-4 reader cutover was
  a no-op: zero external readers found during Strategy B, confirming subsystem
  self-containment.

### Deferred (out of scope)

- **Groups A+B config-constants** — 6 scalar parameters (`till_swtill`, `till_i_n_model`,
  `till_iRedist`, `till_Max_Z_tillage`, `till_Ntill`, `till_Ntypes`) and 12 event/type-table
  arrays. Already in typed `soil_tillage_t`; future config-consolidation arc.
- **Bdens/ParamVG cross-subsystem coupling** — tillage writes them as side-effects but
  soilhydraulics, solute, and oxygenstress read them. Cross-subsystem ownership
  clarification deferred.
- **DoTillage(4) closure** — `continue` body, never called from `swap.f90`; placeholder.
- **Hardcoded unit numbers** in DoTillage(3) and Adapt_WC_H (units 124, 222–226, 333,
  444) — pre-existing output issue, not introduced by this arc.

## Consequences

- **13 runtime-state `till_*` globals retired** from `variables.f90` with
  `[SS-TIL] retired 2026-05-12` provenance markers. Breakdown: 7 per-layer allocatable
  arrays (Group C) + 3 per-step scalars (Group D) + 3 init-once geometry/cursor integers
  (Group E).
- **NEW flat `tillage_state_t`** in `src/state/tillage_state.f90`; aggregated under
  `swap_state_t` as `state%tillage`.
- **`tillage_init(state, numlay)`** wired at `swap.f90:228` (between `atmosphere_init` and
  `DoTillage(1)`).
- **Only 2 compile-driven passes** (Strategy B) — smallest retirement in the series.
  T-4 reader cutover was a no-op: tillage subsystem fully self-contained, zero external
  readers of `till_*` runtime fields across the codebase.
- **H-4 latent bug fixed in-arc.** `set_iTill`'s tautological condition corrected; byte-
  identical in regression (all 5 cases swtill=0).
- **6 tasks total** (T-1 through T-6) — smallest task count in the migration series.
- pFUnit: 736 passing, 0 failures, 1 disabled.
- check-full: 5/5 byte-identical at every commit.

## References

- Discovery: `docs/superpowers/specs/2026-05-12-state-migration-tillage-discovery.md`
- Design: `docs/superpowers/specs/2026-05-12-state-migration-tillage-design.md`
- Plan: `docs/superpowers/plans/2026-05-12-tillage-state-migration.md`
- Predecessors: ADR 0034 (heat — flat-layout precedent), ADR 0035 (boundary — first
  coupling-surface arc), ADR 0036 (crop-uptake), ADR 0037 (atmosphere — state-arg windfalls
  that pre-plumbed DoTillage/Change_MvGpars/Adapt_WC_H/Consolidate_Bdens), ADR 0038
  (soil-water core — Strategy B origin)
- Arc commits: `2fb05de` (T-1) → `8151570` (T-2) → `41dbdbb` (T-3) → `4f9e008` (T-4)
  → `c84c349` (T-5)
