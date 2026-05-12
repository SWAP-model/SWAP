# Tillage State-Type Migration — Design

**Date:** 2026-05-12
**Status:** accepted (pending implementation)
**Migration #:** 9 of N — small finishing arc following the soil-water 4-arc decomposition
**Branch:** `refactor/tillage-state`
**Discovery:** `docs/superpowers/specs/2026-05-12-state-migration-tillage-discovery.md`
**Predecessor ADRs:** 0035 (boundary), 0036 (crop-uptake), 0037 (atmosphere), 0038 (soil-water core)
**Successor ADR:** 0039

## Goal

Carve the 13 genuine runtime-state `till_*` globals (Groups C+D+E) out of `variables.f90` into a new flat `tillage_state_t`, aggregated as `state%tillage` under `swap_state_t`. This is the smallest arc to date: single home file (446 LoC), self-contained subsystem (8 external read/write sites across 3 files), Phase 0 is zero (ADR 0021 already completed the TOML port), and all entry points already carry `state` as `intent(inout)` from prior arcs. Fix the `set_iTill` tautological-condition latent bug in-arc (H-4). Groups A+B (18 config-constant `till_*` fields already represented in typed `soil_tillage_t`) are deferred to a future config-consolidation arc.

## Out of scope

- **Groups A+B (18 config-constants):** `till_swtill`, `till_i_n_model`, `till_iRedist`, `till_Max_Z_tillage`, `till_Ntill`, `till_Ntypes` plus the 12 event/type-table arrays. Already represented in typed `soil_tillage_t` (`soil_config.f90`). Future config-consolidation arc.
- **Bdens and ParamVG:** written by tillage as cross-subsystem side-effect coupling but consumed by `soilhydraulics.f90`, `solute.f90`, and `oxygenstress.f90`. Cross-subsystem ownership clarification deferred; stay as legacy globals (H-1).
- **DoTillage(4) closure:** never called from `swap.f90`; `continue` body only. Leave as-is.
- **Hardcoded unit numbers in DoTillage(3) and Adapt_WC_H (H-2, H-5):** writes to units 222/224/226/124/333/444 without consistent `TEST` guards. Pre-existing issue; fix as non-arc housekeeping before or during T-3 if the output section is touched anyway.

## Context

After the soil-water 4-arc decomposition (boundary ADR 0035, crop-uptake ADR 0036, atmosphere ADR 0037, soil-water core ADR 0038), the main soil-water global namespace is retired. Tillage is the natural next finishing arc: it is unusually self-contained (no external compute readers of `till_*` fields outside tillage.f90 itself; `config_to_variables.f90` is the sole co-writer; no output subsystem readers). `DoTillage`, `Change_MvGpars`, `Adapt_WC_H`, and `Consolidate_Bdens` already carry `state intent(inout)` — a complete state-arg windfall from prior arcs. Phase 0 is zero: ADR 0021 ported all tillage input to TOML and the typed `soil_tillage_t` config record already covers every input parameter. Discovery Section 2 confirms all 31 `till_*` globals divide cleanly into 18 config-constants (Groups A+B, not in scope) and 13 runtime-state fields (Groups C+D+E, this arc).

## Decisions

### D1. Pilot scope: tillage home tree

Home file: `src/crop/tillage.f90`. Co-writer: `src/io/toml/config_to_variables.f90` for the pre-init seed (Groups A+B legacy globals remain; no change to the adapter in this arc).

### D2. State-type shape — flat

All 13 migrated fields are instantaneous (no `flzerointr`/`flzerocumu` cadence). `tillage_state_t` is a flat record. No cohorts. Matches the heat (ADR 0034) and boundary (ADR 0035) flat-layout precedent.

### D3. Aggregation

`swap_state.f90` adds `type(tillage_state_t) :: tillage` alongside the existing `surfacewater`, `drainage`, `solute`, `heat`, `soilwater`, `atmosphere` fields.

### D4. Scope: 13 runtime-state fields (Groups C+D+E)

- **Group C — 7 per-layer allocatables (sized to `NumLay`):**
  `till_Rho_tillage`, `till_Rho_cons`, `till_Rho_last`, `till_K_R_cons`, `till_Rho_match`, `till_N_match`, `till_Slope_match`
- **Group D — 3 per-step scalars:**
  `till_sumDWC`, `till_sumAvail1`, `till_sumAvail2`
- **Group E — 3 init-once geometry/cursor fields:**
  `till_MaxNumSoilHo`, `till_MaxNumSoilCP`, `till_iTill`

### D5. Groups A+B (18 config-constants) deferred

Already in typed `soil_tillage_t`. Migrating them to `state%tillage` would duplicate them. A config-consolidation arc (or the "retire all legacy globals" final arc) is the correct venue. Deferring keeps this arc at 6 tasks.

### D6. Bdens/ParamVG NOT in scope

Cross-subsystem side-effect writes by tillage. Stay as legacy globals. See H-1 in discovery.

### D7. `tillage_init(state%tillage, numlay)` placement

New subroutine in `src/state/tillage_state.f90` (or in `tillage.f90` as a public subroutine). Insert call at `swap.f90` between `atmosphere_init` (line 214) and `DoTillage(1)` (line 227). Allocates the 7 per-layer Group C arrays to `numlay`. Group E geometry scalars (`MaxNumSoilHo`, `MaxNumSoilCP`) are computed by `det_MNSH` during `DoTillage(1)` — not at `tillage_init` time. Signature: `tillage_init(state, numlay)` where `state` is `type(tillage_state_t), intent(out)`.

### D8. Pre-init pattern

`config_to_variables.f90` continues writing the 18 Groups A+B `till_*` config-constant legacy globals (pre-init window; `state%tillage` not yet allocated). `DoTillage(1)` seeds the Group E cursor/geometry fields (`till_iTill`, `till_MaxNumSoilHo`, `till_MaxNumSoilCP`) into `state%tillage` after `tillage_init` has run. Same pre-init pattern as atmosphere A-2.6's transient-buffer approach, but simpler — Groups A+B stay as legacy globals throughout this arc (no buffer needed for Group C/D since they are freshly allocated and computed by DoTillage(1) itself).

### D9. H-4 latent bug fix

`set_iTill` (tillage.f90 lines 405–414) contains a tautological condition: `t1900 >= Date_tillage(i-1) .and. t1900 < Date_tillage(i-1)`. The same array index on both sides means the branch is always false; `iTill` never advances past event 1. Fix in-arc (T-3): change second index from `i-1` to `i`. Safety: if none of the 5 check-full regression cases use tillage (flTillage false everywhere), the fix is byte-identical by construction. Verify before committing T-3.

### D10. Compile-driven Phase 2 (Strategy B)

Apply Strategy B from soil-water core (ADR 0038 playbook lesson #1): comment out the 13 `till_*` runtime-state globals in `variables.f90`; the compiler enumerates remaining legacy references. Modest scope — expect 5–10 hidden readers (vs soil-water's ~80+). Discovery found 8 external sites across 3 files; undercount by 20–50% is typical (expect 3–5 additional sites).

## Phasing

**Phase 0:** None. ADR 0021 completed TOML port. Zero config gaps.

**Phase 1 (T-1 through T-3):**
- T-1: Create `src/state/tillage_state.f90` with `tillage_state_t` (13 fields). Add `state%tillage` to `swap_state_t`. pFUnit lifecycle tests.
- T-2: Add `tillage_init` and wire at `swap.f90` between `atmosphere_init` and `DoTillage(1)`.
- T-3: Dual-write all 13 runtime fields inside `tillage.f90`. Fix H-4 `set_iTill` bug. ASSOCIATE with `tl_` prefix.

**Phase 2 (T-4 through T-5):**
- T-4: Reader cutover — migrate the 8 external read sites in 3 files to `state%tillage%X`.
- T-5: Strategy B compile-driven retirement — comment out 13 `till_*` globals; iterate compile-fail → fix.

**Close-out (T-6):** ADR 0039 + playbook update + merge prep.

## Testing

**pFUnit:** new `tests/unit/state/test_tillage_state.pf` — default-value tests for 10 scalars/integers, allocation lifecycle for 7 per-layer allocatables, aggregator access via `swap_state_t`. Register in `testSuites.inc` and `meson.build`. All existing tests remain green at every commit.

**check-full:** 5/5 byte-identical at every commit. H-4 fix is byte-identical if no regression case uses multi-event tillage (flTillage likely false in all 5 cases — verify before T-3 commit).

## Non-goals

- Rewriting tillage physics. ASSOCIATE preserves variable names; compute bodies change only in dual-write mechanics.
- Retiring Groups A+B `till_*` config-constant globals. Future arc.
- Moving `Bdens`/`ParamVG` into any typed state. Cross-subsystem arc.
- Adding regression coverage for multi-event tillage scenarios.

## ADR 0039 stub

Records: tillage as migration #9; flat `tillage_state_t` (13 fields: 7 per-layer allocatables + 3 per-step scalars + 3 init-once geometry/cursor); Groups A+B deferral rationale; Bdens/ParamVG stay-legacy rationale; H-4 `set_iTill` fix; Strategy B compile-driven retirement; no Phase 0. Cross-references discovery, design, plan, ADRs 0035–0038.
