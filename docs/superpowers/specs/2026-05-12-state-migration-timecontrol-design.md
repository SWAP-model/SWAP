---
title: TimeControl state-migration — Design
date: 2026-05-12
status: accepted
adr: 0041 (pending)
discovery: docs/superpowers/specs/2026-05-12-state-migration-timecontrol-discovery.md
predecessor: ADR 0039 (tillage)
---

# TimeControl State-Type Migration — Design

## Goal

Carve 61 runtime-state TimeControl globals (clock/schedule/print-flag machinery) into a new `state%timecontrol` typed record. Migrate the largest reader inventory of any arc (~415 sites across 22 files) using Strategy B compile-driven retirement. Defer 18 config-constants and 2 cross-subsystem reset gates (flZeroIntr/flZeroCumu). Preserve byte-identical check-full at every commit.

## Out of scope

- **18 config-constants** (tstart, tend, dtmin, dtmax, nprintday, period, swres, swodat, swheader, swscre, outdat, outdatint, msteps, MaxIt, MaxIterTime, nmetdetail, swmetdetail, swrain, swetsine, swdra, swhea, swsnow, swsolu, swirfix) — already in typed config; deferred to future config-consolidation arc. Same precedent as tillage Groups A+B (ADR 0039).
- **flZeroIntr / flZeroCumu** (2 reset gates) — read in 25 files / 28 sites. Cross-subsystem reset orchestration is its own arc. TimeControl OWNS the gates (writes them in TC(case=2)); other subsystems CONSUME them via `if (flZeroIntr) ...` patterns. This arc leaves both as legacy globals; a future "reset-orchestration arc" will consolidate.
- **dtEventRain** — written by meteodt, read by TC. Cross-ownership stays legacy.
- **macropore TimeControl reads** — N/A (subsystem retired per ADR 0040).
- **Output state for retired writers** — N/A (deleted per ADR 0009 Phase 5+).

## Context

TimeControl is the most cross-cutting arc to date in field reach: `dt` (255 reads), `t1900` (120), `period` (90) — these are referenced from every physics file. The arc would have been larger before the recent cleanup (ADR 0009 Phase 5+ purge + macropore retirement removed ~85 reader sites). Post-purge, `waterbalance.f90` is the heaviest single reader (69 sites), overtaking the once-dominant swapoutput.f90.

The compute side (TimeControl(task=N) routine) is already largely state-aware from prior arcs — `state` is passed via `intent(in)`. This arc promotes the intent to `inout` and plumbs `IterTime` (the only TC-adjacent routine without state). The bulk of work is reader cutover, not signature surgery.

## Decisions

- **D1. Scope.** Migrate 61 runtime-state globals to a new `timecontrol_state_t`. Defer 18 config-constants + 2 flzero gates + dtEventRain (see "Out of scope").

- **D2. Shape — flat layout.** TimeControl OWNS reset gates rather than consumes them. Per-day scatter is single-purpose and doesn't justify a `reset_per_day()` cohort method. Init-once flags fit alongside per-step flags. Match atmosphere/tillage precedent for non-cohort arcs.

- **D3. Aggregation.** `swap_state.f90` adds `type(timecontrol_state_t) :: timecontrol`.

- **D4. No `timecontrol_init`.** Init happens inline in `TimeControl(case=1)`. This subsystem's "initial state" is fundamentally tied to the runtime entry point, unlike the per-node-array allocators in soilwater/atmosphere. The aggregator type is allocated with default scalars; TC(case=1) populates from config.

- **D5. `TimeControl(state)` intent bump.** Promote from `intent(in)` (from SS-SWC S-2.12B) to `intent(inout)`. Update single call site in swap.f90.

- **D6. `IterTime(state)` plumbing.** `IterTime(task)` currently has no state arg but writes `tc_tmptimestart` / `tc_tmptimeend`. Add `state` as `intent(inout)`. Single call site (swap.f90 timestep loop).

- **D7. Pre-init pattern for `iyear`, `imonth`, `dt`.** `config_to_variables.f90` writes these three legacy globals before any state allocator runs (no state struct exists at that point). Use the transient-buffer pattern from atmosphere A-2.6: keep the writes, drain into state during TC(case=1).

- **D8. DLL re-init path** (swap.f90:592-598). External-caller re-init mutates `iyear`. Add direct `state%timecontrol%iyear` write at the same site; legacy write stays until retirement (Strategy B).

- **D9. swapoutput co-writes `flheader` / `flheadirg` / `flIrg1Start`.** Mechanical dual-write at the surviving output write sites; these are output-state flags that TC reads in formatting decisions.

- **D10. flZeroIntr / flZeroCumu — leave as legacy.** Per "Out of scope": these gates are owned by TC but consumed by 25 other files (28 sites). Migrating them in this arc would force unrelated reader plumbing into the same commit. A future "reset-orchestration arc" handles them.

- **D11. Strategy B compile-driven retirement.** Following the success of SS-SWC S-2.12B (28 passes), atmosphere A-2.6 (13 passes), soil-water core retirement, and tillage T-5 (2 passes), this arc uses Strategy B for the final retirement step. Estimated 10-18 compile passes given cross-subsystem reach.

- **D12. Apply atmosphere pre-flight dual-write coverage check.** Before each reader cutover task, verify each migrated field has non-zero state writes (atmosphere lesson: nraidt/aintcdt regression).

- **D13. Heaviest reader is waterbalance.f90 (69 sites).** It needs its own dedicated task. Post-purge, this overtook swapoutput.f90 (53).

- **D14. Monolithic arc, not split.** 15 tasks fits one arc. Decomposition by file would create artificial boundaries (e.g., what file-group ordering matters for byte-identical?). The compile-driven approach surfaces ordering as needed.

## Phasing

- **Phase 1** — state-type creation, aggregation, signature bumps, home-tree dual-write. Field-by-field state writes in TimeControl, IterTime, and DLL re-init path.
- **Phase 2** — reader cutover across 22 files, ordered by site count (heaviest first), with pre-flight dual-write checks. Strategy B retirement of 61 globals; compile-driven discovery cleans up surprises.
- **Phase 3** — ADR 0041, playbook, merge.

## Testing

- pFUnit: new `test_timecontrol_state.pf` with ~8 tests (defaults, aggregator access, IterTime state write, init lifecycle).
- check-full: 5/5 byte-identical at every commit.
- Hupselbrook canary: < 5 seconds runtime confirms no NR-divergence regression.

## Non-goals

- Do NOT touch the 18 config-constants. They have typed homes already; cleanup is a separate arc.
- Do NOT migrate flZeroIntr/flZeroCumu reset orchestration. Future reset-orchestration arc.
- Do NOT touch dtEventRain ownership.
- Do NOT delete unused file-unit globals in variables.f90 (afo/aun/vap/bal/...). Separate cleanup pass.

## ADR 0041 stub

The arc will land ADR 0041 capturing:
- 61 runtime globals retired; 18 config-constants and 2 reset-gates deferred.
- TimeControl/IterTime intent bumps.
- DLL re-init pre-init pattern (transient buffer for iyear/imonth/dt).
- Strategy B compile-driven retirement (4th application).
- New playbook lessons (if any).
