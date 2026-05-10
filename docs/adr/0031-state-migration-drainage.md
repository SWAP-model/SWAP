---
title: "ADR 0031 — Drainage state-type migration"
date: 2026-05-10
status: accepted (draft — Phase 1 complete; Phase 2 pending)
---

# ADR 0031: Drainage state-type migration

**Status:** accepted (Phases 1 + 2 complete)
**Date:** 2026-05-10
**Migration #:** 2 of N (subsystem-by-subsystem state migration umbrella)
**Branch:** `refactor/surfacewater-state` — continuing the umbrella migration branch; merge to `development` after all subsystems land

## Context

The surface-water state-type migration (ADR 0030) established the playbook: discovery → design → plan → execute, with state-type definition, argument-threaded `swap_state_t`, ASSOCIATE in compute, dual-write transitional pattern, check-full byte-identical at every commit.

ADR 0030 acknowledged a piece of technical debt: `qdra` and `qdrain` were temporarily classified into `surfacewater_state_t` because surfacewater compute writes them. The drainage subsystem is the architecturally correct owner — surfacewater consumes the result. Drainage is the natural next migration; it inherits this cleanup task plus its own owned globals.

The drainage discovery doc (`docs/superpowers/specs/2026-05-09-state-migration-drainage-discovery.md`) cataloged 6 owned globals (`qdrd`, `qdrain`, `qdra`, `drainl`, `wetper`, `ztopdislay`), 52 borrowed globals, 6 external reader files, 9 cross-subsystem hazards. The owned set is much smaller than surface-water's (29 owned), but the migration is more involved because of the `qdra`/`qdrain` re-classification and `divdra`'s explicit-shape `(Madr, macp)` argument constraint that has kept the legacy `qdra` global alive.

## Decision

**Migrate drainage off `variables.f90` globals** into a typed `drainage_state_t`, aggregated under `swap_state_t` alongside `surfacewater_state_t`. Move `qdra` and `qdrain` from `surfacewater_state_t` to `drainage_state_t`. Modernize `divdra` to take assumed-shape arrays so the legacy `qdra(Madr, macp)` global can be deleted in Phase 2. Continue using ASSOCIATE blocks in compute bodies, the dual-write transitional pattern, and check-full as the integration gate.

Two phases:

- **Phase 1** (commits `406b55e..53d2776`, complete): define `drainage_state_t`, add to `swap_state_t`, move `qdra`/`qdrain` from `surfacewater_state_t` to `drainage_state_t` (mechanical rename of consumers), introduce `drainage_init(state)` allocating per-level arrays, modernize `divdra` to assumed-shape, dual-write the 4 remaining owned fields (`qdrd`, `drainl`, `wetper`, `ztopdislay`), drop the dual-write — drainage state is authoritative for those 4 fields.
- **Phase 2** (pending): owner-rule relocations (`swdislay==2` `ztopdislay` rule from surfacewater into drainage — tentatively a no-op since discovery showed no active write site), migrate cross-subsystem readers (`frozencond`'s `qdrain` writeback, etc.), drop `qdra` and `qdrain` dual-writes, fix the `geofac` adapter ordering bug at `config_to_variables.f90:431`, delete the 6 owned globals from `variables.f90`. ADR 0031 finalize.

### Phase 1 outcomes (completed 2026-05-10)

- `src/state/drainage_state.f90` defines `drainage_state_t` with 6 fields (1 scalar + 5 allocatable arrays). `swap_state_t` aggregates it.
- `qdra` (the only field that was actually in `surfacewater_state_t`) moved to `drainage_state_t`. `qdrain` was never in `surfacewater_state_t` to begin with — only the qdra rename was needed across 10 files. Consumers in drainage, surfacewater, frozencond, soilhydraulics, waterbalance, solute, output now read `state%drainage%qdra`.
- `drainage_init(state)` introduced as the lifecycle init entry, called from swap_main BEFORE `SurfaceWater(1)`. Allocates the per-level arrays from drainage-config dimensions (`nrlevs`, `numnod`). `surfacewater_init` no longer allocates `qdra`.
- `divdra` (in `module distribute_drainage`) modernized: argument arrays `FluxDr` and `FluxDrComp` (mapped to `qdrain`/`qdra` at call sites) changed from explicit-shape to assumed-shape. Body needed no `size()` replacements because all loop bounds were already scalar parameters (`NumDrHlp`, `NumComp`). Local temporaries still use `Madr`/`macp` constants — out of scope for this arc.
- Drainage compute (`bocodre`) dual-writes `qdrd` and `drainl`/`wetper` to state; `ztopdislay` was found to have no per-timestep write sites (the D7 hazard was inactive in current source — simplifies Phase 2). After Task 7 dropped the dual-write, drainage state is authoritative for these 4 fields; the legacy globals are unwritten by compute (still declared in variables.f90).
- Discovery audit gap caught by check-full: Task 6 audit declared "zero output readers" of `qdrd`, but `WLEVBAL` and `WBALLEV` in `surfacewater.f90` are *compute* readers of `qdrd`. Task 7's check-full failure (GWL drift ~0.36-0.91cm) surfaced this; both readers were migrated to `state%drainage%qdrd` to recover. Lesson: the cross-subsystem reader inventory (Section 3.5) should also include compute consumers, not just output.
- Tests: 622 pFUnit (was 617 after surfacewater migration; +4 drainage_state + +1 swap_state-has-drainage). All green.
- check-full: 5/5 byte-identical at every commit.

**Phase 2 outcomes (completed 2026-05-10, commits `a1180f9..d248cbb`):**

- 6 drainage-owned globals commented out (with provenance markers) in `variables.f90` and removed from `initialize.f90`: `qdrd`, `qdrain`, `qdra`, `drainl`, `wetper`, `ztopdislay`. `state%drainage%*` is the sole authoritative location.
- Cross-subsystem reader migrations:
  - `frozencond.f90:FrozenBounds` reads/writes `state%drainage%qdrain` and `state%drainage%qdra` (Task 1).
  - `surfacewater.f90:WLEVBAL` reads `state%drainage%qdrain` (Task 2). The `divdra` call site at line 143 passes typed state (Task 2).
  - `waterbalance.f90:integral` reads `state%drainage%qdrain` (Task 2) and `state%drainage%qdra` (caught in Task 5).
  - `drainage.f90:drainage()` divdra call passes typed state (Task 3); `bocodrb` and `bocodre` write `state%drainage%qdrain` directly via ASSOCIATE (Task 4); `drainage_init` signature gains `config` arg, reads from typed config instead of legacy globals (Task 5).
- Adapter bug fix: `config_to_variables.f90:431` now gates the `surface_runoff%geofac` write on `ipos /= 5` (Task 5). Zero observable effect on the 5 regression cases (none uses `ipos=5`); correctness-only fix preserving both intended TOML field behaviors.
- pFUnit: 622 throughout Phase 2 (no count change). check-full: 5/5 byte-identical at every commit.

### Hidden migration gaps caught during Phase 2 (Section 3.5 addendum content)

The Phase 2 plan's pre-flight inventory was incomplete despite the Phase 1 lessons-learned addendum to Section 3.5. Task 5's "zero remaining readers" verification step caught four hidden gaps:

1. `surfacewater.f90` swdislay redistribution used global `qdra` as a working buffer; also `zTopDisLay` was a global there. Migrated to `state%drainage%qdra` and a subroutine-local respectively.
2. `waterbalance.f90:integral` had a `qdra` read for `inqdra` accumulation that the inventory missed.
3. `drainage_init` read `drainl`, `wetper`, `ztopdislay`, `qdrd` from globals to seed state. After global deletion, it had to read from `config` directly — signature gained a `config` argument.
4. `bocodrb` had a `wetper(1)` read of the global; migrated to `state%drainage%wetper(1)`.

The catch was Task 5's verification step doing the grep and discovering the issues before the global declarations were deleted. **Lesson reinforced:** verification grep must run BEFORE the deletion, and any non-zero count is an unmigrated reader to handle in the same task.

## Architectural learnings (Phases 1 + 2)

- **Cross-subsystem reader inventory must be exhaustive AND categorize by intent.** Phase 1's lesson was "include compute readers, not just output." Phase 2 added: "include working-buffer reads (e.g., `swdislay` redistribution using global qdra as a temporary), reads inside _init routines (e.g., drainage_init seeding from globals), and reads passed to call sites (e.g., bocodrb's static-geometry reads)." The grep-then-classify pattern is the safety net.
- **`divdra`'s assumed-shape signature (Phase 1 Task 4) was a clean unblock.** Phase 2 just switched call sites without further refactoring of divdra. The investment paid off.
- **The geofac schema duplication (two TOML fields mapping to one legacy global) was fixed inline.** The deeper schema reconciliation (do we need both fields?) is documented as future work — none of the regression cases exercise the bug.
- **drainage_init signature change (config arg) is a precedent for future state init routines.** Init routines that previously read globals to seed state now take config directly. This breaks the "init reads globals" convention but is architecturally cleaner — config is the authoritative source for config-derived initial values.
- **The dual-write transitional pattern works, but verification gaps compound.** Phase 1's missed compute readers AND Phase 2's missed working-buffer/init reads were both caught only at the integration gate (Task 7 of Phase 1, Task 5 of Phase 2). The pattern: each phase's first integration gate is where invisible debt surfaces. Plan accordingly — the integration-gate task usually exceeds its planned scope.

## Architectural learnings (carried forward to migration #3+)

- **Cross-subsystem reader inventory must include COMPUTE readers, not just output.** Surface-water's playbook update (Section 3.5) caught this kind of issue retroactively; drainage Phase 1 still missed compute readers of `qdrd` because the discovery's Section 3.5 categorization was output-biased. Future subsystem discoveries: scan compute files (frozencond, solute, waterbalance, etc.) explicitly for reads of owned globals. The grep template:
  ```bash
  grep -rEn "use variables.*\b<owned-var>\b" src/ --include="*.f90" | grep -v "src/core/variables.f90 ; src/state/ ; src/io/"
  ```
- **`divdra`-style modernization is cheap when argument names happen to differ from the global names.** `FluxDr` and `FluxDrComp` had no name collision with the legacy globals; the only changes were the bounds declarations. If divdra's arguments had been called `qdra` directly, the body would have needed more care.
- **Allocation-point migrations (e.g., from surfacewater_init to drainage_init) are best done in the same task as the field move.** Task 3 bundled the state-type declaration move + drainage_init introduction + consumer rename. Doing them separately would have left an unallocated state%drainage%qdra during a transitional commit.
- **The dual-write pattern's effectiveness depends on the audit being thorough.** Task 6's "zero output readers" claim was right but incomplete (compute readers existed too). The pattern still works — check-full catches the gap — but the rework is wasted effort. Better Section 3.5 categorization upfront prevents this.

## References

- Discovery: `docs/superpowers/specs/2026-05-09-state-migration-drainage-discovery.md`
- Design: `docs/superpowers/specs/2026-05-10-state-migration-drainage-design.md`
- Plan (Phase 1): `docs/superpowers/plans/2026-05-10-drainage-state-phase1.md`
- Predecessor: ADR 0030 (surface-water pilot established the playbook)
- Phase 1 commit chain: `406b55e` (state type) → `c3472e2` (aggregator) → `f66509e` (qdra move) → `4b29b3b` (divdra modernized) → `c07bb7b` (dual-write) → `53d2776` (drop dual-write) → `469252d` (ADR draft)
- Phase 2 plan: `docs/superpowers/plans/2026-05-10-drainage-state-phase2.md`
- Phase 2 commit chain: `a1180f9` (frozencond) → `1d06be5` (surfacewater + waterbalance + qdrain dual-write add) → `723fb7c` (drainage divdra call) → `d3644a3` (drop qdra/qdrain dual-write) → `d248cbb` (delete globals + geofac fix)
