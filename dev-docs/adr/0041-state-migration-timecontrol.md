---
title: "ADR 0041 — TimeControl state-type migration"
date: 2026-05-12
status: accepted
---

# ADR 0041: TimeControl state-type migration

**Status:** accepted (all tasks complete)
**Date:** 2026-05-12
**Migration #:** 10 of N — most cross-cutting arc to date; 4th Strategy B application
**Branch:** `refactor/timecontrol-state`

## Context

TimeControl is the central time-stepping driver of SWAP: it advances the simulation
clock (`t1900`, `tcum`, `t`), maintains day/year/month/cumulative counters
(`daynr`, `daycum`, `iyear`, `imonth`), governs the variable-step Newton–Raphson
timestep (`dt`, `dtprevious`, `dtEvent`, `tEvent`), and emits the cross-subsystem
gate flags (`flDayStart`, `flDayEnd`, `flRunEnd`, `flYearStart`, `flOutput`,
`flBalOutput`, …). Every other physics file reads at least one of these fields.

This is the **10th state-migration arc** in the rescue-arc series (after surface
water, drainage, solute, heat, boundary, crop-uptake, atmosphere, soil-water core,
tillage, and macropore retirement). It is the **most cross-cutting arc to date**:
TimeControl state reaches every physics file via `dt`/`t1900`/`daynr`/`daycum`/
`period`. The field inventory pre-purge was ~500 external read sites; after the ADR
0009 Phase 5+ output-writer purge and macropore retirement (ADR 0040), the
post-purge confirmed count is **~415 sites across 22 files** — still the largest
reader inventory of any migration arc by a large margin.

Unlike soilwater-core (ADR 0038), where `dt` and soil-physical-property arrays fed
the Richards NR loop via utility-module reads that Strategy B alone could not catch,
**TimeControl reads are syntactically explicit** — all 22 reader files import via
`use variables, only: <field>` or bare `use variables`. Strategy B compile errors
surface all of them with no hidden NR-coupled hang risk.

Prior arcs had progressively plumbed `state` into every entry point, so
`TimeControl(task, state)` already took `state intent(in)` since SS-SWC. This arc
promotes the intent to `inout` and adds `state` to `IterTime`.

## Decision

**Migrate 62 runtime-state TimeControl globals** (16 clock/calendar + 19
timestep/schedule + 27 runtime-evaluated booleans) into a new flat
`timecontrol_state_t`, aggregated as `state%timecontrol` under `swap_state_t`.
Defer 18 config-constants (already in typed `simulation_config_t`), 2 reset gates
(`flZeroIntr`/`flZeroCumu` — future reset-orchestration arc), and `dtEventRain`
(meteodt cross-ownership). Apply Strategy B compile-driven retirement as the master
retirement method (4th application of the pattern).

### Decisions D1–D14

- **D1. Scope: 62 runtime-state globals.** Groups C (16 clock+calendar), D (19
  timestep+schedule), E (27 runtime-evaluated booleans). One extra boolean surfaced
  during Strategy B compile passes versus the 61 estimated in the design.

- **D2. Flat layout — no cohort sub-records.** TimeControl OWNS the `flZeroIntr`/
  `flZeroCumu` gates; it does not consume them. All owned fields are either per-step
  mutated, per-day mutated, monotonic counters, or init-once derived booleans. No
  `flzero*`-gated accumulation cluster exists. Matches heat (ADR 0034), boundary
  (ADR 0035), and tillage (ADR 0039) flat precedents.

- **D3. Aggregation.** `swap_state.f90` adds `type(timecontrol_state_t) :: timecontrol`
  alongside existing soilwater, atmosphere, heat, drainage, solute, tillage, etc.

- **D4. No separate `timecontrol_init` routine.** TimeControl already has a
  case=1 dispatcher that IS the init entry point (lines 73–258). State initialization
  happens inline in `TimeControl(case=1)` — same as the existing case-1 body. A
  separate `*_init` subroutine is only warranted when there is no case-1 dispatcher
  (atmosphere precedent). This avoids splitting init logic across two routines.

- **D5. `TimeControl(state)` intent bump.** Promoted from `intent(in)` (SS-SWC
  S-2.12B windfall) to `intent(inout)`, enabling writes to `state%timecontrol`.
  All six call sites in `swap.f90` and `swapoutput.f90` updated.

- **D6. `IterTime(task, state)` plumbing.** `IterTime(task)` previously had no
  state argument; it read/wrote `tc_tmptimestart`/`tc_tmptimeend` via the legacy
  namespace. Added `state intent(inout)` as the second argument. All three call
  sites in `swap.f90` updated.

- **D7. Pre-init transient-buffer pattern for `iyear`, `imonth`, `dt`.** These
  three fields are written by `config_to_variables.f90` before any state init runs
  (`config_to_variables` executes at `swap.f90:184`; `TimeControl(1, state)` at
  `swap.f90:191`). Module-level buffers `tc_iyear_init_buf`, `tc_imonth_init_buf`,
  `tc_dt_init_buf` in `config_to_variables.f90` receive the writes; `TimeControl
  (case=1)` drains them into state. Generalizes from soilwater-core's
  `h_init_buf`/`pondini_init_buf` (ADR 0038 lesson #4).

- **D8. DLL re-init path (`swap.f90:592–598`).** The external-caller re-init block
  re-derives `iyear` from `Tstart` and calls `TimeControl(1, state)` again. Added
  direct `state%timecontrol%iyear = datea(1)` write at the same site (legacy write
  retained until Strategy B retirement); TC case=1 seeds everything else from config.

- **D9. `swapoutput.f90` co-writes (`flheader`, `flheadirg`, `flIrg1Start`).**
  Mechanical dual-writes at the surviving output write sites. SwapOutput already
  carries `state` (SS-SWC windfall), so the additions are single-line additions.

- **D10. `flZeroIntr` / `flZeroCumu` left as legacy globals.** These gates are
  owned by TimeControl (written in TC case=2) but consumed by 25+ files at 28
  combined sites — every subsystem's `*_state_t%intr%reset()` / `%cumu%reset()`
  call chain reads them. Migrating them here would force unrelated reader plumbing
  into this arc. A future "reset-orchestration arc" will convert them into return
  values from `state%timecontrol%advance()`.

- **D11. Strategy B compile-driven retirement (4th application).** Comment out all
  62 runtime-state globals in `variables.f90`; compiler enumerates every remaining
  legacy reference as a compile error; iterate fix-recompile. 8 compile passes to
  full clean link. Compare: boundary 8 passes, atmosphere 13, soil-water core 28,
  tillage 2, TimeControl 8. Pattern at full maturity — estimate of 10–18 proved
  high; 22 files but most had already received reader cutover in dedicated tasks.

- **D12. Pre-flight dual-write coverage check (atmosphere lesson #1).** Before
  each reader cutover task, verified each migrated field had at least one non-zero
  `state%timecontrol%<field> =` write beyond zero-init. Applied at TC-6 through
  TC-13 to prevent the nraidt/aintcdt-style regression (ADR 0037 finding).

- **D13. Largest reader inventory: ~415 sites, 22 files.** Reader cutover
  distributed across 9 dedicated tasks (TC-6 through TC-13), ordered by site count:
  waterbalance.f90 (69), swapoutput.f90 (53), surfacewater.f90 (34),
  soilhydraulics.f90 (33), readmeteo.f90 (32), cropgrowth.f90 (31), meteodt.f90
  (25), meteoday.f90 (22), initialize.f90 (20), et.f90 (19), drainage.f90 (18),
  irrigation.f90 (18), boundtop.f90 (17), config_to_variables.f90 (15),
  temperature.f90 (12), swap.f90 (11), tillage.f90 (10), cropgrass_init.f90 (10),
  solute.f90 (9), management_soil.f90 (9) and smaller files.

- **D14. Monolithic arc — not split by field group.** Groups C/D/E are
  tightly-coupled (dt/t1900/daynr/period are read together in most physics files);
  splitting into sub-arcs would amplify dual-write maintenance across PRs with no
  architectural payoff. 15 tasks shipped as one arc.

## Out of scope

- **18 config-constants** (`tstart`, `tend`, `dtmin`, `dtmax`, `nprintday`, `period`,
  `swres`, `swodat`, `swheader`, `swscre`, `outdat`, `outdatint`, `msteps`, `MaxIt`,
  `MaxIterTime`, `nmetdetail`, `swmetdetail`, `swrain`, `swetsine`, `swdra`, `swhea`,
  `swsnow`, `swsolu`, `swirfix`) — already have typed homes in `simulation_config_t`,
  `meteorology_config_t`, and per-subsystem config types. Future config-consolidation
  arc. Same precedent as tillage Groups A+B (ADR 0039 D5).

- **`flZeroIntr` / `flZeroCumu` (2 reset gates)** — TC writes these but 25 files
  (28 sites) consume them as the canonical reset-gate signal for every subsystem's
  cohort `reset()` calls. A future "reset-orchestration arc" converts them to
  return values from `state%timecontrol%advance()`.

- **`dtEventRain`** — written by `meteodt.f90`, read by TimeControl at line 392.
  Cross-ownership with the atmosphere/meteodt subsystem; stays as legacy global.
  Will migrate when a future meteo refactor delivers `state%atmosphere%dtEventRain`.

- **`fldecdt`** — already in `timestep_control_mod` (since SS-SWST); not a
  `variables.f90` global and not in scope.

## Consequences

**62 runtime-state globals retired** from `variables.f90` with `[SS-TC]` provenance
markers. Breakdown: 16 clock+calendar (Group C) + 19 timestep+schedule (Group D) +
27 runtime-evaluated booleans (Group E; one extra surfaced in Strategy B passes).

**NEW flat `timecontrol_state_t`** in `src/state/timecontrol_state.f90` — 62 fields
with explicit defaults; no cohort sub-records. Aggregated under `swap_state_t` as
`state%timecontrol`.

**`TimeControl(state)` intent(in) → intent(inout).** `IterTime` gained `state
intent(inout)` as second argument. Both are the only TimeControl signature changes
in this arc.

**Pre-init transient-buffer pattern** (`tc_iyear_init_buf`, `tc_imonth_init_buf`,
`tc_dt_init_buf`) for the three fields seeded by `config_to_variables.f90` before
the state allocator runs. Generalizes the soilwater-core `h_init_buf` pattern.
Drained at the top of `TimeControl(case=1)`.

**ASSOCIATE constraint encountered in `boundtop.f90` (PONDRUNOFF block).** Fortran
ASSOCIATE blocks require single-entry/single-exit semantics. `PONDRUNOFF` has
multiple early `return` statements inside the associate scope, causing ICE. Fixed
by declaring a local variable seeded from state at routine top, using the local
throughout. Same pattern required in `management_soil.f90` select-case-with-early-
returns. Future arcs should audit for this before using ASSOCIATE in such routines.

**4th application of Strategy B compile-driven retirement.** 8 compile passes
(estimated 10–18; actual lower because 9 reader-cutover tasks pre-migrated the
heavy sites before retirement). Compare: boundary 8, atmosphere 13, soil-water
core 28, tillage 2. Pattern at full maturity.

**No NR-coupled hidden-reader hang.** TimeControl reads are syntactically explicit
`use variables, only: <field>` imports in all 22 reader files. Strategy B compile
errors surface everything directly — no utility-module pointer pattern needed
(contrast soilwater-core's `soilhydraulicsutils.f90` NR-coupled hang risk, ADR 0038
lesson #2).

**File counts and metrics:**

- Tasks: 15 (TC-1 through TC-15), comparable to atmosphere (17); lighter Phase 1
  (no cohort design) but heavier Phase 2 (22 reader files vs atmosphere's 9).
- Strategy B compile passes: 8.
- pFUnit: 733 passing, 0 failures, 0 disabled (net: new `test_timecontrol_state.pf`
  suite added).
- check-full: 5/5 byte-identical at every commit throughout the arc.
- Hupselbrook canary: 1.14 s (stable; no NR-divergence regression).

## Cross-references

- Predecessors:
  - ADR 0035 (boundary — flat-layout precedent)
  - ADR 0037 (atmosphere — pre-flight dual-write coverage check lesson)
  - ADR 0038 (soil-water core — Strategy B origin; transient-buffer pattern)
  - ADR 0039 (tillage — config-constant vs runtime-state distinction)
  - ADR 0040 (macropore retirement — reduced reader inventory from ~500 to ~415)
- Deferred work:
  - `flZeroIntr`/`flZeroCumu` — future reset-orchestration arc
  - 18 config-constants — future config-consolidation arc
  - `dtEventRain` — future meteodt/atmosphere refactor
- Arc commits: `ef90ffa` (TC-2) → `2d7182c` (TC-3) → `255108f` (TC-4) →
  `50fad7d` (TC-5) → `8612a2e` (TC-6) → `2bd402e` (TC-7) → `aa7da57` (TC-8) →
  `7217a3f` (TC-9) → `d95f092` (TC-10) → `334ae86` (TC-11) → `3d3eff3` (TC-12) →
  `858fd67` (TC-13) → `c0d503a` (TC-14)
