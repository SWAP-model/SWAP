---
title: "ADR 0030 — Surface-water state-type migration (pilot)"
date: 2026-05-09
status: accepted (draft — finalized at end of Phase 2)
---

# ADR 0030: Surface-water state-type migration (pilot)

**Status:** accepted (Phase 1 complete; Phase 2 pending — consequences section grows once globals leave variables.f90)
**Date:** 2026-05-09
**Branch:** `refactor/surfacewater-state` — do not merge to `development` until Phase 2 lands and verification is green

## Context

`variables.f90` is the legacy module-level state container for SWAP — 1314 LoC, ~1200 globals, declared `save`, accessed by every compute and output subroutine via `use variables, only: …`. This single shared mutable surface prevents multi-instance use, blocks Python embedding via a clean state object, and couples subsystems implicitly through name-level globals.

The plan is to migrate subsystem-by-subsystem onto typed state types, threaded as arguments. Each subsystem's compute, output, and config converge on a single state slice. `variables.f90` shrinks per arc until empty.

Surface-water was chosen as the pilot subsystem because it's peripheral (bottom-of-column, outside the crop/water-balance feedback loops), has dedicated home files (~1433 LoC), and was small enough to be tractable while being substantial enough to validate the playbook.

The discovery doc (`docs/superpowers/specs/2026-05-08-state-migration-surfacewater-discovery.md`) cataloged 29 owned globals, 54 borrowed globals, 6 entry points, and 8 cross-subsystem hazards. The design doc (`docs/superpowers/specs/2026-05-09-state-migration-surfacewater-design.md`) resolved the hazards and split the work into two phases.

## Decision

**Migrate surface-water off `variables.f90` globals** into a typed `surfacewater_state_t`, aggregated under a top-level `swap_state_t` and threaded through subroutine signatures from `swap_main` down. Use `ASSOCIATE` blocks inside leaf compute routines to minimize churn (variable names preserved as aliases; `sw_*` prefix where bare names would shadow `use variables` imports).

Output routines read from the typed state; compute writes the typed state. Replace the `fldecdt` global write (the "decrease timestep" signal) with an `intent(out) :: request_smaller_dt` parameter on `SurfaceWater(task=2)`. Move `outswb`'s year-end `swstini = swst` writeback (the output-layer-mutates-state hazard) into a dedicated `surfacewater_year_reset(state%surfacewater)` callback invoked from `swap_main`.

Two phases:

- **Phase 1** — define the typed state, thread it through every entry point, dual-write state and globals during compute, redirect output reads to typed state, drop dual-write from surface-water-home files. check-full byte-identical at every task. **Phase 1 complete.**
- **Phase 2** — owner-rule relocation (`qdrain` rule to drainage; `l(Madr)` m→cm to drainage config load), removal of dead `SurfaceWater(2)` call from `swapoutput.f90`, migration of cross-subsystem readers in `drainage.f90` / `surfacewaterutils` utility funcs / `solute.f90` / `frozencond.f90` / `soilgrid.f90` / `waterbalance.f90:integral`, removal of remaining surface-water owned globals from `variables.f90`, migration of `timecontrol.f90`'s `fldecdt` readers and removal of the `fldecdt` global. **Phase 2 pending.**

Phase 0 (promote 12 missing config fields for `swman=2` / `swqhr=2` / `swsec=1` branches) was originally planned and **deferred**: those branches are stub-errored, the regression cases never exercise them, and none of the 29 state fields overlap with the 12 config inputs. Track as its own future arc.

### Reusable playbook

The migration produced a four-phase template for subsequent subsystems:

1. **Discovery** — populate the discovery template (`docs/superpowers/specs/2026-05-08-state-migration-surfacewater-discovery.md` is the worked example): source files, owned globals, borrowed globals, entry points, internal call graph, output coupling, config inputs, hazards, test surface, summary statistics.
2. **Design** — resolve each hazard, define the state-type fields, define the lifecycle, decide phasing.
3. **Plan** — mechanical task decomposition with check-full as the integration gate at every commit. Dual-write transitional pattern.
4. **Execute** — subagent-driven implementation, two-stage review, byte-identical check-full proves the migration is functionally complete.

The discovery template needs one extension: a "Section 3.5: external readers of owned globals" — globals that THIS subsystem writes but OTHER compute subsystems read. Phase 1 of the surface-water pilot found these only at execution time (e.g., `drainage.f90` reads `wls`, `swst`, `iqdra`, `cqdra`, `cqdrain`, `qdra`, `ZDraBas`, `flInitDraBas`; `solute.f90` reads `qdra`; `waterbalance.f90:integral` accumulates `iqdra`/`cqdra`/`cqdrain`). Future discovery should pre-catalog these so the design phase can plan their migration explicitly.

## Consequences

**Phase 1 outcomes:**

- 29 surface-water-owned globals are now mirrored in `state%surfacewater` (the typed state record).
- 14 of those 29 globals have their compute writes fully migrated — globals are now read-only from the perspective of the surface-water home files.
- 15 globals still have dual-writes because of cross-subsystem readers in `drainage.f90`, `surfacewaterutils.f90` utility functions, `soilgrid.f90`, `frozencond.f90`, `solute.f90`, `waterbalance.f90:integral`. Phase 2 migrates those readers and drops the remaining dual-writes.
- Output routines (`csv_out`, `set_values`, `fill_values`, `outwba`, `outinc`, `outbal`, `outdrf`, `outend`, `outswb`, `outage`, `outblc`, `AgeTracerOutput`, `swapoutput`, `surfacewateroutput`, `soilwateroutput`) read from `state%surfacewater` instead of legacy globals. The output layer is read-only against the typed state.
- `outswb`'s year-end `swstini = swst` writeback (Phase 1 hazard 4) replaced by `surfacewater_year_reset(state%surfacewater)` callback in `swap_main`.
- `fldecdt` global write is preserved because `timecontrol.f90` still reads it; Phase 2 migrates that reader and drops the global. `request_smaller_dt` (intent(out)) is the new typed signal.
- `OutputModflow` (a previously-undocumented surface-water consumer found during Phase 1) holds a `SAVE`-local `state_om`; this is a transitional shim acceptable for Phase 1.
- pFUnit: 614 → 616 tests (added 2 state-type construction suites). All green.
- check-full: 5/5 byte-identical at every task boundary.

**Phase 2 consequences (to be filled after completion):** TBD — remaining global deletions, cross-subsystem reader migrations, owner-rule relocations, and the actual `variables.f90` line-count drop.

**Architectural consequences (carried into subsequent subsystem migrations):**

- `swap_state_t` is the central state surface. Each subsystem migration adds one field. Eventually it becomes the Python embedding's primary state object.
- `ASSOCIATE` is the right tool for low-churn leaf-routine refactoring, but `use variables, only:` shadowing forces an alias prefix (`sw_*`) when the same name is imported. Future subsystems will face the same constraint until their borrowed-globals are also migrated.
- `intent(out)` parameters replace shared "signal" globals (e.g., `fldecdt` → `request_smaller_dt`). This is the recommended pattern for cross-subsystem signals during migration.
- The discovery template should be extended (Section 3.5) to pre-catalog external readers of owned globals — this was the largest source of surprise scope in Phase 1.

## References

- Discovery: `docs/superpowers/specs/2026-05-08-state-migration-surfacewater-discovery.md`
- Design: `docs/superpowers/specs/2026-05-09-state-migration-surfacewater-design.md`
- Plan (Phase 1): `docs/superpowers/plans/2026-05-09-surfacewater-state-phase1.md`
- ADR 0024 — `dtutil.f90` as TTutil-API compatibility shim (the predecessor cleanup that established subsystem boundaries)
- ADR 0029 — `dtutil` de-shim from physics layer (pending; orthogonal cleanup)
- Phase 1 commit chain: `acae2af` (state type) → `62d8347` (aggregator) → `8da337e` (init dual-write) → `eb8adde` (entry refactor) → `a8ee668` (output reads round 1) → `0c16410` (output reads round 2 / waterbalance dual-write add) → `846a7d7` (drop dual-write from surface-water home files)
