---
title: "ADR 0030 — Surface-water state-type migration (pilot)"
date: 2026-05-09
status: accepted (draft — finalized at end of Phase 2)
---

# ADR 0030: Surface-water state-type migration (pilot)

**Status:** accepted (Phases 1 + 2 complete)
**Date:** 2026-05-09
**Branch:** `refactor/surfacewater-state` — ready to merge to `development`

## Context

`variables.f90` is the legacy module-level state container for SWAP — 1314 LoC, ~1200 globals, declared `save`, accessed by every compute and output subroutine via `use variables, only: …`. This single shared mutable surface prevents multi-instance use, blocks Python embedding via a clean state object, and couples subsystems implicitly through name-level globals.

The plan is to migrate subsystem-by-subsystem onto typed state types, threaded as arguments. Each subsystem's compute, output, and config converge on a single state slice. `variables.f90` shrinks per arc until empty.

Surface-water was chosen as the pilot subsystem because it's peripheral (bottom-of-column, outside the crop/water-balance feedback loops), has dedicated home files (~1433 LoC), and was small enough to be tractable while being substantial enough to validate the playbook.

The discovery doc cataloged 29 owned globals, 54 borrowed globals, 6 entry points, and 8 cross-subsystem hazards. The design doc resolved the hazards and split the work into two phases.

## Decision

**Migrate surface-water off `variables.f90` globals** into a typed `surfacewater_state_t`, aggregated under a top-level `swap_state_t` and threaded through subroutine signatures from `swap_main` down. Use `ASSOCIATE` blocks inside leaf compute routines to minimize churn (variable names preserved as aliases; `sw_*` prefix where bare names would shadow `use variables` imports).

Output routines read from the typed state; compute writes the typed state. Replace the `fldecdt` global write (the "decrease timestep" signal) with an `intent(out) :: request_smaller_dt` parameter on `SurfaceWater(task=2)`. Move `outswb`'s year-end `swstini = swst` writeback (the output-layer-mutates-state hazard) into a dedicated `surfacewater_year_reset(state%surfacewater)` callback invoked from `swap_main`.

Two phases:

- **Phase 1** — define the typed state, thread it through every entry point, dual-write state and globals during compute, redirect output reads to typed state, drop dual-write from surface-water-home files. check-full byte-identical at every task. **Phase 1 complete.**
- **Phase 2** — owner-rule relocation (`qdrain` rule to drainage; `l(Madr)` m→cm to drainage config load), removal of dead `SurfaceWater(2)` call from `swapoutput.f90`, migration of cross-subsystem readers in `drainage.f90` / `surfacewaterutils` utility funcs / `solute.f90` / `frozencond.f90` / `soilgrid.f90` / `waterbalance.f90:integral`, removal of remaining surface-water owned globals from `variables.f90`, migration of `timecontrol.f90`'s `fldecdt` readers and removal of the `fldecdt` global. **Phase 2 pending.**

Phase 0 (promote 12 missing config fields for `swman=2` / `swqhr=2` / `swsec=1` branches) was originally planned and **deferred**: those branches are stub-errored, the regression cases never exercise them, and none of the 29 state fields overlap with the 12 config inputs. Track as its own future arc.

### Reusable playbook

The migration produced a four-phase template for subsequent subsystems:

1. **Discovery** — populate the discovery template: source files, owned globals, borrowed globals, entry points, internal call graph, output coupling, config inputs, hazards, test surface, summary statistics.
2. **Design** — resolve each hazard, define the state-type fields, define the lifecycle, decide phasing.
3. **Plan** — mechanical task decomposition with check-full as the integration gate at every commit. Dual-write transitional pattern.
4. **Execute** — subagent-driven implementation, two-stage review, byte-identical check-full proves the migration is functionally complete.

The discovery template needs one extension: a "Section 3.5: external readers of owned globals" — globals that THIS subsystem writes but OTHER compute subsystems read. Phase 1 of the surface-water pilot found these only at execution time (e.g., `drainage.f90` reads `wls`, `swst`, `iqdra`, `cqdra`, `cqdrain`, `qdra`, `ZDraBas`, `flInitDraBas`; `solute.f90` reads `qdra`; `waterbalance.f90:integral` accumulates `iqdra`/`cqdra`/`cqdrain`). Future discovery should pre-catalog these so the design phase can plan their migration explicitly.

**Cohort partitioning by activity gate** (lesson added 2026-05-10 from the SS-CRR Phase A correction). When cumulatives nested in one subsystem's state are accumulated under different activity flags (e.g., surfacewater's `flSurfaceWater` vs drainage's `fldrain`), split the cohort by flag — not by name prefix or by data type. The original Task A5 baked a single 7-field `surfacewater_cumulative_t` and treated drainage's strict-subset reset as a "mid-accumulation hazard" papered over with an inline comment. The real reason was that `flSurfaceWater = .true.` only when `swdra=2`, so the 3 reservoir fields never accumulate under `swdra=1` and need no reset site in `Drainage()`. The fix split the cohort into `surfacewater_drainage_cumulative_t` (gated by `fldrain`) and `surfacewater_reservoir_cumulative_t` (gated by `flSurfaceWater`). The owner of each cohort is the subsystem whose activity flag gates that cohort's accumulation, and only the owner calls `reset()`. **Heuristic for future migrations:** during discovery, enumerate the activity flags that gate each cumulative field's accumulation paths. If two fields in the same subsystem have different gates, they belong in different cohorts. See ADR 0033 for the full pattern.

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

**Phase 2 outcomes (completed 2026-05-09):**

- 28 of 29 surface-water-owned globals removed from `variables.f90` and `initialize.f90`. The 29th (`qdra`) is retained as a globally-declared array because `divdra` consumes it as an explicit-shape `(Madr, macp)` argument; the state%surfacewater%qdra mirror is the live read source for compute and output, while the global serves only as the divdra argument-binding scaffold. This residual global is a candidate for cleanup when divdra's signature is modernized in a future arc.
- `fldecdt` (timestep-decrease signal) deleted from `variables.f90` and `initialize.f90`. Replaced by `intent(out) :: request_smaller_dt` on `SurfaceWater(2/3)`, propagated to a swap-main-local flag (`fldecdt`) hosted in a new tiny `src/core/timestep_control_mod.f90` module rather than threaded as an argument — `headcalc` writes the flag on Richards non-convergence, which would have ballooned argument lists if pure-arg threading were forced.
- 13 cross-subsystem reader files migrated to `state%surfacewater` for surface-water-owned reads: `drainage.f90`, `macropore.f90`, `macrorate.f90`, `frozencond.f90`, `solute.f90`, `soilgrid.f90`, `management_soil.f90`, `soilhydraulics.f90`, `waterbalance.f90` (`integral` and `fluxes`), `surfacewaterutils.f90` (`wlevst`, `swstlev`, `qhtab`, `runoff`), `swap.f90`, `swapoutput.f90` (including `outafo`/`outaun`/`outend`/`outwba`/`outinc`/`outbal`/`outdrf`/`outage`/`outblc`/`OutputModflow`/`AgeTracerOutput`/`SurfaceWaterOutput`/`soilwateroutput`), `swap_csv_output.f90` (`csv_out`/`set_values`/`fill_values`), `boundtop.f90`, `timecontrol.f90`.
- Owner-rule relocations completed:
  - `qdrain(:) = 0` rule (when `gwl > 998`) moved from `SurfaceWater(2)` in `surfacewater.f90` into `Drainage` (and a guarded variant in `bocodre`) per spec D7. surface-water no longer mutates `qdrain`.
  - `l(Madr)` m→cm conversion moved from `surfacewater_init` in-place mutation to TOML-read-time in `read_drainage_toml.f90` per spec D6. The typed config now holds cm internally; the legacy in-place global mutation is gone. A double-conversion path in the adapter (DRAMET=3+swdivd=1) was discovered and eliminated as a side effect.
- `OutputModflow`'s SAVE-local `state_om` confirmed as the correct architectural shape (Phase 2 Task 8) — it's a finite-difference `dV/dGWL` perturbation experiment, not a timestep-tracking loop. Documentation added; code unchanged.
- `variables.f90` LoC: 1314 → ~1340 (the 28 surface-water declarations are commented out with `! Moved to surfacewater_state_t%X` markers preserving migration provenance; net file size is similar, but the live declarations dropped accordingly). The commented-out block is itself a candidate for cleanup in a follow-up sweep once subsequent subsystems migrate.
- Tests: 614 → 617 pFUnit tests passing (added 5 surfacewater_state, 2 swap_state, 1 read_drainage_toml conversion). Zero regressions across check-full's 5 cases.

**Phase 2 architectural learnings carried forward:**

- The 16-of-29 globals had cross-subsystem readers (drainage, macropore, soilhydraulics, frozencond, solute, soilgrid, management_soil, waterbalance, swapoutput, swap_csv_output, boundtop, timecontrol). The original Phase 1 discovery only enumerated this subsystem's writes, not external reads — Section 3.5 of the discovery doc retroactively captures these for future subsystem migrations.
- ASSOCIATE name shadowing of `use Variables` (no `only:`) imports works in gfortran but is fragile; the `sw_*` prefix pattern is safer.
- `intent(inout)` for non-trivial subroutines is contagious upward — `headcalc` had to gain inout because it calls `MACROPORE` which now takes inout state. Plan ahead for the call chain.
- Some globals (e.g., `qdra`) are structurally tangled with explicit-shape array arguments (`divdra`) and cannot be cleanly removed without modernizing the consumer signatures — a separate future arc.
- `timestep_control_mod` as a tiny dedicated module is a useful pattern when a flag is written deep in the call chain — argument-threading would have required signature changes in 4-5 routines, while the module is one new file with two readers.

**Architectural consequences (carried into subsequent subsystem migrations):**

- `swap_state_t` is the central state surface. Each subsystem migration adds one field. Eventually it becomes the Python embedding's primary state object.
- `ASSOCIATE` is the right tool for low-churn leaf-routine refactoring, but `use variables, only:` shadowing forces an alias prefix (`sw_*`) when the same name is imported. Future subsystems will face the same constraint until their borrowed-globals are also migrated.
- `intent(out)` parameters replace shared "signal" globals (e.g., `fldecdt` → `request_smaller_dt`). This is the recommended pattern for cross-subsystem signals during migration.
- The discovery template should be extended (Section 3.5) to pre-catalog external readers of owned globals — this was the largest source of surprise scope in Phase 1.
- The `config` argument proposed in design D4 (`SurfaceWater(task, state, config, request_smaller_dt)`) was NOT implemented. The actual signature is `SurfaceWater(task, state, request_smaller_dt)` — surfacewater config is still accessed via the legacy `variables` globals that the input adapter populates. Adding `config` would have required threading `swap_config_t` through every subsystem entry point in parallel with `state`; deferred to the broader config-passing refactor (ADR 0016 future direction). Subsequent subsystem migrations may apply or postpone this same way.

## Post-merge follow-up items (tracked as a punch-list for the next touch of this code)

1. Narrow the bare `use Variables` in `SurfaceWater(task)` dispatcher body (`src/drainage/surfacewater.f90:36`) to a focused `use variables, only: …` clause. The branch corrected this pattern elsewhere; only the home-file dispatcher was missed. One-commit cleanup, no functional change.
2. Plan the `qdra` global removal as a sub-task of the drainage subsystem migration. When `divdra` signature is modernized to take an allocatable array (or `state%surfacewater%qdra` directly), the global can be deleted.
3. Extend the migration-playbook discovery template's Section 3.5 guidance with the `sw_*` ASSOCIATE prefix convention for subroutines that carry `use Variables, only:` clauses — to prevent latent shadow risks.

## References

- ADR 0024 — `dtutil.f90` as TTutil-API compatibility shim (the predecessor cleanup that established subsystem boundaries)
- ADR 0029 — `dtutil` de-shim from physics layer (pending; orthogonal cleanup)
- Phase 1 commit chain: `acae2af` (state type) → `62d8347` (aggregator) → `8da337e` (init dual-write) → `eb8adde` (entry refactor) → `a8ee668` (output reads round 1) → `0c16410` (output reads round 2 / waterbalance dual-write add) → `846a7d7` (drop dual-write from surface-water home files)
