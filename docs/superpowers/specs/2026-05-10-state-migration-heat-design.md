# Heat State-Type Migration — Design

**Date:** 2026-05-10
**Status:** accepted (pending implementation)
**Migration #:** 4 of N
**Branch:** `refactor/heat-state`
**Discovery:** `docs/superpowers/specs/2026-05-10-state-migration-heat-discovery.md`
**Predecessor playbook:** ADR 0030 + state-migration-playbook.md (Phase A correction lesson on activity gates)
**Successor ADR:** 0034

## Goal

Migrate the heat (soil-temperature) subsystem off `variables.f90` globals into a typed `heat_state_t` aggregated under `swap_state_t`. Promote 6 missing config fields to close the silent-zero-defaults correctness gap on `swcalt=1` runs. Fix the `outheapar` output-side writeback hazard. Establish heat as the simplest migration so far — small, contained, no cumulatives, no cohorts needed.

## Out of scope

- Adding regression coverage for `swcalt=1` (analytical method). None of the 5 cases activates it.
- Fixing the 10 ungated `tsoil` readers in compute files — they consume all-zero temperatures when `flTemperature=false`. Pre-existing bug pattern, not migration-introduced. Document as known-issue.
- Migrating `rfcp` to `state%heat` is in scope, but the broader soilhydraulics state migration is not.

## Context

Heat is migration #4. The discovery doc cataloged 13 owned globals, **all instantaneous** (no `flzero*` gating). There are NO cumulative or intermediate fields. The cohort sub-record pattern from ADR 0033 doesn't apply — the state type is flat.

Distinctive concerns:
- `Temperature(task)` is the **only main compute entry not yet plumbed with state** — every other subsystem we've migrated had at least partial state plumbing from earlier arcs. Plumbing this is Phase 1's pivot.
- `FrozenCond` and `FrozenBounds` already take `state` (added during surfacewater Phase 2 Task 5). They need migration of their reads/writes to `state%heat`.
- `outheapar()` in `swapoutput.f90` writes `heacap` from inside an output routine via a `devries(thetadum, heacap, heacnd)` call — classic writeback hazard (playbook gotcha #6).
- `rfcp` is co-reset to `1.0` by `soilhydraulics.f90:884` (`soilwater` init) on every Richards timestep. Multi-subsystem co-write.
- 6 config fields missing for the `swcalt=1` (analytical) path — same class as solute Phase 0.

## Decisions

### D1. Pilot scope: heat home tree

The four home files: `temperature.f90`, `frozencond.f90`, `heat_config.f90`, `read_heat_toml.f90`. Plus call-site adjustments in compute consumers.

### D2. State-type shape — flat (no cohorts)

`heat_state_t` in `src/state/heat_state.f90`. All 13 fields are instantaneous, so no `intermediate` / `cumulative` sub-records. Plain flat layout:

```fortran
type :: heat_state_t
   ! Per-node soil-temperature arrays (allocatable, sized numnod)
   real(real64), allocatable :: tsoil(:)        ! current soil temperature per node
   real(real64), allocatable :: tsoilold(:)     ! previous-step temperature (if applicable)
   ! Heat-physics scalars and arrays the discovery enumerates
   real(real64), allocatable :: heacap(:)       ! per-node heat capacity (computed each step)
   real(real64), allocatable :: heacnd(:)       ! per-node heat conductivity
   ! ... (final field set comes from discovery Section 2)

   real(real64), allocatable :: rfcp(:)         ! reduction factor for cold-period (per node)
   ! ... etc

end type heat_state_t
```

Implementer derives the precise list from discovery Section 2's 13 owned globals. None go in cohorts because none are flag-gated cumulative.

### D3. State-type aggregation

`src/state/swap_state.f90` adds `type(heat_state_t) :: heat`.

### D4. Argument-threading: plumb `Temperature(task)`

`Temperature(task)` is currently `subroutine Temperature(task)` with no state. Add `type(swap_state_t), intent(inout) :: state` per the established pattern. Update its caller (`swap.f90`).

`FrozenCond(state, ...)` and `FrozenBounds(state, ...)` already take state — Phase 1 just migrates their `state%heat%*` reads/writes from globals.

ASSOCIATE in compute bodies. Use `ht_*` prefix for shadow-safe aliasing per the playbook (sw_, dr_, sl_, now ht_).

### D5. `rfcp` migration includes soilhydraulics co-write fix

`rfcp` lives at `soilhydraulics.f90:884` as a per-step reset to `1.0`. After migration, `state%heat%rfcp` is the authoritative location.

`soilhydraulics.f90` already takes `state` as an argument (from surfacewater Phase 2 Task 5 + drainage Phase 2). Update its `rfcp = 1.0d0` line to `state%heat%rfcp = 1.0_real64` (or `state%heat%rfcp(:) = 1.0_real64` if rfcp is per-node — verify shape).

This is a Phase 2 task (cross-subsystem reader migration covers it), not Phase 1.

### D6. `outheapar` writeback fix — move `devries` call to a callback

`outheapar()` in `swapoutput.f90` is an OUTPUT routine that calls `devries(thetadum, heacap, heacnd)` and writes to `heacap`. After migration, `state%heat%heacap` is in the typed state — output should be **read-only**.

Two options:
- **(A)** Introduce a `heat_param_init(state%heat, ...)` (or similar) callback in `temperature_mod` that runs the `devries` computation. Call from compute side; `outheapar` reads `state%heat%heacap`.
- **(B)** Inline the `devries` call where it physically belongs (probably `Temperature(task=1)` init), and `outheapar` becomes read-only.

Recommend **(B)** — the `devries` computation is part of heat-init physics, not output formatting. Move it to where `Temperature(task=1)` runs. Verify by inspection that this doesn't change observable output (the values should be identical).

### D7. Phase 0 — promote 6 missing config fields

Per discovery hazard #4: `swcalt=1` (analytical method) is uncovered by TOML pipeline. Missing config fields:
- 4 scalars: `ddamp`, `tmean`, `tampli`, `timref`
- 2 boundary-table arrays: `temtoptab`, `tembtab`

Same class as solute Phase 0. Promote upfront — patches the silent-defaults correctness gap. No regression coverage (the activity is stub-errored on the TOML path; future arc enables `swcalt=1`).

### D8. 10 ungated `tsoil` readers — documented as out-of-scope known issue

Compute files in `crop/`, `solute/`, `atmosphere/snow.f90`, `soilhydraulics`, `boundtop`, `boundbottom` read `tsoil` without checking `flTemperature`. When heat is disabled, they read all-zero soil temperatures. **Pre-existing bug pattern** — not introduced by migration.

Per the surfacewater migration's pattern of accepting hazards as out-of-scope (e.g., the 12 deferred `swman=2` config fields), document this in ADR 0034 as a known issue. Future arc fixes the readers either via gating or via making `flTemperature=true` the default.

## Phasing

Three phases bundled into one overarching plan (heat is small enough to plan together, like solute):

- **Phase 0** — Promote 4 scalar + 2 array config fields into `heat_config_t`. Patches `swcalt=1` correctness gap.
- **Phase 1** — Define `heat_state_t`, add to `swap_state_t`, plumb `Temperature(task)` with state, dual-write 13 owned globals, migrate output reads + the `outheapar` writeback fix. Drop dual-write within heat-home files (when feasible).
- **Phase 2** — Cross-subsystem reader migration (10 compute-reader files), `rfcp` co-write fix in `soilhydraulics`, drop remaining dual-writes, comment out 13 globals from `variables.f90`. ADR 0034.

## Testing

- **pFUnit:** new `test_heat_state` (default values + array allocation lifecycle). New tests for the 6 promoted config fields' validators. Migrate any existing heat tests that read globals.
- **check-full:** 5/5 byte-identical at every commit. The `swcalt=1` path is inactive in regression cases, so check-full proves the migration doesn't break the active path.

## Non-goals (explicit)

- Rewriting heat physics. ASSOCIATE preserves variable names; bodies of `Temperature`, `FrozenCond`, `FrozenBounds` should not change semantically.
- Adding regression coverage for `swcalt=1`.
- Fixing the 10 ungated `tsoil` readers (consumer-side bug, separate workstream).
- Migrating other subsystems' globals just because they live near heat (e.g., `soilhydraulics` only touched for the `rfcp` co-write fix; broader soilhydraulics state migration is its own arc).

## ADR 0034

Records:
- Decision: heat is migration #4; flat state type (no cohorts) — first migration without cohort sub-records.
- `Temperature(task)` plumbed with state — closes the last main compute entry without state.
- Phase 0 physics config patch (6 fields).
- `outheapar` writeback hazard fixed by moving `devries` to compute side.
- `rfcp` migration with soilhydraulics co-write fix.
- 10 ungated `tsoil` readers documented as known issue.
- Cross-references to discovery doc, ADRs 0030–0033, plan.
