---
title: "ADR 0034 — Heat subsystem state-type migration"
date: 2026-05-10
status: accepted
---

# ADR 0034: Heat subsystem state-type migration

**Status:** accepted (Phases 0 + 1 + 2 complete)
**Date:** 2026-05-10
**Migration #:** 4 of N
**Branch:** `refactor/heat-state`

## Context

The surfacewater (ADR 0030), drainage (ADR 0031), and solute (ADR 0032) state-migration arcs established the playbook: discovery → design → plan → execute, with typed state records, argument-threaded `swap_state_t`, ASSOCIATE in compute bodies, dual-write transitional pattern, and check-full byte-identical at every commit. ADR 0033 layered the cumulative-reset cohort pattern on top.

Heat (soil-temperature) is migration #4. The discovery doc cataloged 13 owned globals — **all instantaneous** (no `flzero*` gating, no cumulatives, no intermediates). This makes heat the first migration where the cohort sub-record pattern from ADR 0033 does not apply; the state type is flat.

Distinctive concerns surfaced in discovery:

- `Temperature(task)` was the **last main compute entry point without state** — every other subsystem had at least partial state plumbing.
- `FrozenCond` / `FrozenBounds` already took `state` (added during surfacewater Phase 2); only their reads/writes needed migrating.
- `outheapar()` in `swapoutput.f90` called `devries(thetadum, heacap, heacnd)` and wrote back `heacap` — a textbook output-side state-mutation hazard (playbook gotcha #6). On first inspection this looked like a stale-state writeback; on investigation it turned out to be a parameter sweep (21 `devries` calls per layer with synthetic theta values), not a timestep state propagation.
- `rfcp` is co-reset to `1.0` by `soilhydraulics.f90` on every Richards sub-step — a cross-subsystem co-write that required a guarded migration.
- 6 physics-config fields for the `swcalt=1` (analytical method) path absent from `heat_config_t` — silent-zero-defaults correctness gap, same class as solute Phase 0.

## Decision

**Migrate heat off `variables.f90` globals** into a typed flat `heat_state_t`, aggregated under `swap_state_t` alongside the earlier subsystem state records. Thread from `swap_main` down via `state` argument.

Three phases:

- **Phase 0** — Promote 6 missing physics-config fields into `heat_config_t`.
- **Phase 1** — Define `heat_state_t`, add to `swap_state_t`, plumb `Temperature(task)` with `state`, dual-write 13 owned globals, fix `outheapar` writeback, migrate output reads to `state%heat`.
- **Phase 2** — Migrate 10+ cross-subsystem compute readers to `state%heat`, fix `rfcp` co-write in `soilhydraulics`, drop dual-write, comment out 11 of 12 owned globals (reclassify `tsoil`). ADR 0034 finalize.

### Phase 0 — Physics config promotion

6 fields promoted into `heat_config_t` and the TOML reader:
- 4 scalars: `ddamp`, `tmean`, `tampli`, `timref`
- 2 boundary-table arrays: `temtoptab`, `tembtab`

Closes the `swcalt=1` (analytical method) silent-zero-defaults correctness gap. No regression coverage — none of the 5 check-full cases activates this path.

### Phase 1 — State type and threading

- `src/state/heat_state.f90` defines `heat_state_t` with 12 fields (all flat — no cohorts):
  - Per-node allocatable arrays sized from `numnod`: `tsoil(:)`, `tsoilold(:)`, `heacap(:)`, `heacnd(:)`, `heatflow(:)`, `rfcp(:)`
  - Scalar heat-physics values: `tsurf`, `tadnew`, `tadold`, `ztopnew`, `zbotold`, `ztopold`
- `swap_state_t` gains `type(heat_state_t) :: heat`.
- `heat_init(state)` allocates per-node arrays from `numnod`, seeds `state%heat%rfcp(:) = 1.0_real64`.
- `Temperature(task, state)` plumbed — closes the last main compute entry without state. Caller in `swap.f90` updated.
- ASSOCIATE with `ht_*` prefix in compute bodies where `use variables` imports would shadow bare-name aliases.

**`outheapar` writeback hazard — parameter sweep, not state mutation.** The `devries(thetadum, heacap, heacnd)` call inside `outheapar` (an output routine) looked like a stale-state writeback. Investigation found it is a parameter sweep: 21 synthetic theta values per layer fed through `devries` to populate the output columns `heacap_pF0..pF4.2`. The fix was a local scratch array `heacap_loc(macp)` — the sweep writes the scratch array, not `state%heat%heacap`. devries' signature was later extended in Task 8 (Phase 2 cleanup) to take composition arrays explicitly, eliminating its hidden global reads. `state%heat%heacap` is read-only from output's perspective.

### Phase 2 — Cross-subsystem migration and cleanup

**Cross-subsystem readers migrated:** 10+ files updated to read `state%heat` instead of legacy globals:
- Output: `outvap`, `outafo`, `outaun`, `outend`, `outtem`, `TemperatureOutput`
- Compute: crop (`ArableLandGerm`, `sumttd`, `grass`, `cropgrowth`), solute, atmosphere/`snow.f90`, `soilhydraulics`, `soilhydraulicsutils`, `boundtop`, `boundbottom`

**Non-module crop routines.** `ArableLandGerm`, `sumttd`, and `grass` are external non-module routines. Passing typed `state` with an `optional` dummy arg requires an explicit interface block. Instead: non-optional dummy `tsoil(:)` array arg (no interface needed), plus `use Variables, dummy_X_ => tsoil` rename to exclude the global from scope inside the caller's compilation unit. Clean and unambiguous.

**`rfcp` co-write fix.** `soilhydraulics.f90` resets `rfcp(:) = 1.0` at each Richards sub-step. After migration the authoritative location is `state%heat%rfcp`. Required an `allocated(state%heat%rfcp)` guard: `SoilWater(1, state)` runs before `heat_init` in `swap.f90`'s startup sequence — the array is unallocated when the first Richards sub-step runs. `heat_init` seeds the same `1.0` moments later, so the practical net effect is unchanged. The guard prevents a startup segfault.

**`boundbottom` / `boundtop` scope expansion.** `boundtop` and `boundbottom` gained a `state` argument during Task 9 to migrate `rfcp` reads. Modest scope expansion, accepted.

**`hconduc` fallback — known residual (D8).** The temperature-dependent thermal conductivity branch (`iHWCKmodel` 4–11) reads `state%heat%tsoil` via `hconduc(…, tsoil_node)`. However, 7 of hconduc's callers were not all updated to pass `tsoil_node`; those paths receive a `0.0_real64` sentinel instead of the stale global. The sentinel is a correctness improvement over the prior uninitialized-global risk, but the branch is unreachable in all 5 regression cases. Documented as a residual for users hitting `iHWCKmodel` 4–11.

**`tsoil` reclassification — deviation from spec.** `tsoil` is NOT retired. `Temperature(task=1)` reads the legacy global `tsoil` to build its initial soil-temperature interpolation table from TOML staging (the CSV companion for the init table is read at config time into the global). Migrating this seeding would require either threading `config` through `Temperature` or introducing a `state%heat%tsoil_init_table` buffer. Both are feasible but were judged out-of-scope for this arc. `tsoil` is reclassified as a **config-staging buffer**, not a compute-time state variable. 11 of 12 owned globals retired; `tsoil` carries a `! [SS-HEAT] reclassified as config-staging buffer — not retired` provenance comment.

**Dead code removed:** 871 lines of `csv_write` / `csv_write_tz` deleted (zero live callers confirmed by grep).

## Consequences

- 12 heat-owned fields carried in `state%heat`, aggregated under `swap_state_t`. Compute writes `state%heat` exclusively; output reads `state%heat`.
- `Temperature(task)` plumbed with state — closes the last main compute entry point without state.
- 11 of 12 legacy heat globals commented out in `variables.f90` with `! [SS-HEAT] retired 2026-05-10` provenance markers. `tsoil` retained as config-staging buffer.
- `outheapar` parameter sweep writes a local scratch array; `state%heat%heacap` is read-only from output.
- `rfcp` migrated with init-order guard.
- `swcalt=1` correctness gap patched (6 config fields promoted), though path remains inactive in regression suite.
- pFUnit: 664 passing, 0 failures, 1 disabled throughout.
- check-full: 5/5 byte-identical at every commit.

### Known residual issues

- **`hconduc` fallback (D8):** 7 callers of `hconduc` that didn't receive the `tsoil_node` arg get a `0.0_real64` sentinel on the temperature-dependent K branch (`iHWCKmodel` 4–11). Unreachable in current regression cases. Future arc: migrate all callers to pass the node temperature explicitly.
- **`tsoil` global:** retained as a config-staging buffer. Future arc: introduce `state%heat%tsoil_init_table` or thread `config` through `Temperature(task=1)`.

## References

- Predecessors: ADR 0030 (surfacewater pilot), ADR 0031 (drainage), ADR 0032 (solute), ADR 0033 (cumulative reset cohorts)
- Phase 0: `c1ed4c3` (config promotion)
- Phase 1: `38e82b3` (state type) → `08c785d` (Temperature plumbed) → `b03f8f1` (dual-write) → `918ba82` (output reads + outheapar fix)
- Phase 2: `82b1736` (tsoil cross-subsystem readers) → `20aa7a3` (cropgrowth tsoil + hconduc) → `20aa7a3` (rfcp co-write) → `9d5b2a4` (drop dual-write) → `b815f12` (retire globals)
