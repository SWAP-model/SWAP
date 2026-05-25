---
title: "TimeControl Modernization — Design Spec"
date: 2026-05-13
status: draft
arc: SS-TCM (timecontrol modernization)
builds-on: SS-DRV Phase 1 (driver-modernization arc, 2026-05-12)
relates-to: ADR 0041 (TimeControl state migration, 2026-05-12)
---

# TimeControl Modernization — Named Procedures + Cleanup

## Context

`src/core/timecontrol.f90` (913 lines) hosts two free subroutines that the new `swap_mod` driver calls via integer dispatch:

- `TimeControl(task, state)` — magic ints 1 / 2 / 3 / 9
  - case (1): initialization (~250 lines)
  - case (2): advance one timestep (~410 lines)
  - case (3): reduce dt under non-convergence (~50 lines)
  - case (9): end-of-day SSDI dt cap (~7 lines)
- `IterTime(task, state)` — magic ints 1 / 2 / 3
  - case (1) / (2) / (3): init / per-day budget check / final stats

These are the same anti-patterns the SS-DRV arc just retired in `subroutine swap`: free subroutine, integer dispatch, `use variables` (legacy globals at line 11), `-std=legacy` compilation. The hardest prerequisite work — extracting 61 runtime fields into `state%timecontrol` — is **already complete** under ADR 0041 / Task TC-14 (committed 2026-05-12). What remains is the procedure-level refactor plus two distinct cleanup opportunities the migration exposed.

### Why now

1. **swap_mod already calls TimeControl 4× and IterTime 3× via magic ints.** After SS-DRV, `swap_run_step` reads `call TimeControl(2, state)` and `call TimeControl(3, state)` — the caller has no way to express intent except through opaque integers. Named procedures restore readability at the call site.
2. **State migration just shipped.** The state-threading argument was the precondition for any procedural cleanup. It's done. The window to layer named-procedure work on top is now.
3. **Two pieces of dead code surfaced.** Both are visible only after the migration; both would require touching this file regardless of arc.
4. **No external interface to design.** Unlike SS-DRV's BMI layer, TimeControl is internal. Scope is purely Fortran-side refactor. Smaller arc than SS-DRV (estimated 60–70% the size).

### What the migration exposed

**Dead code A — pervasive dual-writes.** Inside the existing `associate` block at `timecontrol.f90:54-114`, every TC field is aliased: `flRunEnd => state%timecontrol%flRunEnd`. The body then does:

```fortran
flRunEnd = .false.
state%timecontrol%flRunEnd = flRunEnd      ! same memory — redundant
```

Both writes go to the same location. The second line is a leftover from when fields were being dual-tracked during the migration. **Approximate count: ~120 lines of pure deletion** across case (1) and case (2).

**Dead code B — implicit dispatch override at line 117.**

```fortran
itask = task
if (itask.eq.2 .and. (fldecdt .or. fldecdtmin)) itask = 3
```

A caller passing `task=2` silently becomes task=3 under non-convergence flags. This is exactly the magic the SS-DRV arc killed in `subroutine swap`. **Mitigation:** audit the callers in `swap_mod` — they already gate `TimeControl(3, ...)` calls behind explicit `if (fldecdt .or. (flMacroPore .and. FlDecMpRat))` checks. If the audit confirms the override is either (a) never firing, or (b) always firing through a path the caller has already chosen, the override is redundant and gets deleted.

## Decision

Convert `TimeControl` and `IterTime` to named procedures in `module timecontrol_mod`. Delete the dual-write redundancy. Delete the implicit dispatch override (preconditioned on a one-task audit). Update the 7 call sites in `swap_mod`. The module joins `swap_modern` so the new file compiles with `-std=f2018 -Wall -Wextra`. Regression byte-for-byte parity is non-negotiable.

`flZeroIntr` / `flZeroCumu` are included in this arc. Audit shows the surface is small: all 8 writes are inside `timecontrol.f90` (lines 138, 139, 419, 423, 613, 623, 629, 705, 712, 720, 731, 732) plus 2 init-zero writes in `initialize.f90:28-29`. Eight reader files consume them, and the recently-shipped state-migration arcs already thread `state` into all of them. Migrating the two globals into `state%timecontrol` is a mechanical extension of the same pattern this arc is already applying.

## Architecture

### What goes

- `subroutine TimeControl(task, state)` at `src/core/timecontrol.f90:4` — converted to four named procedures inside `module timecontrol_mod`. The legacy free subroutine is deleted.
- `subroutine IterTime(task, state)` at `src/core/timecontrol.f90:853` — converted to three named procedures inside the same module. The legacy free subroutine is deleted.
- The implicit override at `timecontrol.f90:117` — deleted after the caller audit confirms the semantics are preserved at the call site.
- The bare-global `flZeroIntr` and `flZeroCumu` declarations in `variables.f90:117-118` — deleted. Their init-zero writes in `initialize.f90:28-29` — deleted (the new state fields default to `.false.`, so the explicit init is redundant).
- Dual-write `state%timecontrol%X = X` lines paired with associate-alias writes — deleted throughout case (1) and case (2) bodies (~120 lines).
- `select case (task)` dispatch — deleted; replaced by direct procedure calls.

### What stays

- 61 of the 63 fields in `state%timecontrol` (`src/state/timecontrol_state.f90`) — unchanged. **Two new fields are added: `flZeroIntr` and `flZeroCumu`** (both `logical, default .false.`, joining the Group E "runtime-evaluated boolean flags" cluster).
- All runtime semantics (every `case (N)` body, with the dead code stripped) — preserved verbatim.
- The `associate` block aliasing TC fields — each named procedure carries its OWN associate block scoped to the fields it touches. This is cleaner than one shared block: smaller alias lists per procedure, no orphaned aliases from cases that don't apply, and each procedure's signature reflects only what it reads/writes. `timecontrol_advance` has the largest list (~30 aliases, same fields the current case (2) touches); `timecontrol_day_end` has just 2–3 (`dt`).
- The `use variables` import for non-migrated globals (`dtmin`, `dtmax`, `flMacroPore`, `FlDecMpRat`, `dt_SSDI_event`, `flMaxIterTime`, `msteps`, `project`, output-config switches) — scoped with `only:`. Wholesale `use variables` retirement is a separate arc.

### New module layout

```fortran
module timecontrol_mod
   use swap_state_mod, only: swap_state_t
   implicit none
   private
   public :: timecontrol_init, timecontrol_advance, &
             timecontrol_reduce_dt, timecontrol_day_end
   public :: itertime_init, itertime_check, itertime_close
contains
   subroutine timecontrol_init(state)
   subroutine timecontrol_advance(state)
   subroutine timecontrol_reduce_dt(state)
   subroutine timecontrol_day_end(state)
   subroutine itertime_init(state)
   subroutine itertime_check(state)
   subroutine itertime_close(state)
end module timecontrol_mod
```

Each procedure-body is the corresponding case body from the legacy subroutine, with dual-writes removed and the override deleted. All procedures take `state` as `intent(inout)` (every case mutates `state%timecontrol`).

### Reader-site updates for flZeroIntr / flZeroCumu

Eight files read the two globals. All eight already accept `state` in their relevant subroutine signatures (state-migration arcs threaded it through). Each migration is a localized change to the `use variables` line and the flag access path.

| File | Current pattern | New pattern |
|---|---|---|
| `src/atmosphere/meteoday.f90:440,448,454` | `use variables, only: flzerointr, flzerocumu` | drop from use list; read `state%timecontrol%flZeroIntr` / `flZeroCumu` |
| `src/drainage/surfacewater.f90:103,109` | `if (flzerointr) ...` (state in scope) | `if (state%timecontrol%flZeroIntr) ...` |
| `src/drainage/drainage.f90:429,502,510` | bare globals in arg list / body | `state%timecontrol%flZeroIntr / flZeroCumu` |
| `src/atmosphere/snow.f90:106,114` | bare globals | `state%timecontrol%flZeroIntr / flZeroCumu` |
| `src/soil/soilhydraulics.f90:1136,1148` | bare globals | `state%timecontrol%flZeroIntr / flZeroCumu` |
| `src/soil/waterbalance.f90:443` | bare global | `state%timecontrol%flZeroIntr` |
| `src/solute/agetracer.f90:159` | bare global | `state%timecontrol%flZeroCumu` |
| `src/solute/solute.f90:128,129` | bare globals (state in scope) | `state%timecontrol%flZeroIntr / flZeroCumu` |

If any reader subroutine turns out NOT to have `state` in its signature (a discovery during implementation), threading it through is a localized signature change — the SS-* state migration has established the pattern. Risk: low.

### Caller updates in `swap_mod`

Seven call sites total. Mapping:

| Current call | New call |
|---|---|
| `call TimeControl(1, state)` (swap_init) | `call timecontrol_init(state)` |
| `call TimeControl(2, state)` (swap_run_step, end-of-step) | `call timecontrol_advance(state)` |
| `call TimeControl(3, state)` (swap_run_step, dt-reduction loop) | `call timecontrol_reduce_dt(state)` |
| `call TimeControl(9, state)` (swap_run_step, day-end SSDI) | `call timecontrol_day_end(state)` |
| `call IterTime(1, state)` (swap_init) | `call itertime_init(state)` |
| `call IterTime(2, state)` (swap_run_step, flMaxIterTime gate) | `call itertime_check(state)` |
| `call IterTime(3, state)` (swap_close) | `call itertime_close(state)` |

`swap_mod`'s procedure-scope `use` lists get one new entry each:
`use timecontrol_mod, only: timecontrol_init, ...` (with the relevant subset per procedure).

### Override audit (precondition)

Before deleting the line-117 override, an audit task instruments the legacy subroutine:

```fortran
itask = task
if (itask.eq.2 .and. (fldecdt .or. fldecdtmin)) then
   write(*,'(a,i0,a,L1,a,L1)') &
      '[TC-OVERRIDE] task=2->3 fired; fldecdt=', fldecdt, ' fldecdtmin=', fldecdtmin
   itask = 3
end if
```

Run `pixi run check-full` and inspect the log for `[TC-OVERRIDE]` lines. Three possible outcomes:

1. **No firings.** The override is dead code in current scenarios → safe to delete.
2. **Firings only on a path the caller already gates.** The caller would have called `timecontrol_reduce_dt` explicitly under named-procedure dispatch → safe to delete.
3. **Firings on a path the caller does NOT gate.** Real semantic dependency → either preserve by having `timecontrol_advance` internally invoke `timecontrol_reduce_dt` under the same condition, or fix the caller.

The audit task is mechanical (instrument, run check-full, grep, revert the instrumentation, commit findings as a note). It is the first non-trivial task in the arc.

### Build wiring

The new file `src/core/timecontrol_mod.f90` joins `modern_sources` in `meson.build`:

```python
modern_sources = [
    'src/core/swap_mod.f90',
    'src/core/swap_bmi_mod.f90',
    'src/core/timecontrol_mod.f90',
]
```

The legacy `src/core/timecontrol.f90` is deleted from the `sources` list. The `swap_modern` static library carries the new file under `-std=f2018 -Wall -Wextra`. Cross-standard linkage to `swap_legacy` modules (`variables`, `irrigation_mod`, `timestep_control_mod`) works the same way it does for `swap_mod` today.

### Test coverage

Existing tests already exercise TimeControl end-to-end through the swap lifecycle:

- `test_swap_mod_suite` — 3 lifecycle tests on hupselbrook (touch all 7 TimeControl/IterTime call sites)
- `check-full` regression — 5 cases compare CSV outputs byte-for-byte against fixtures

This arc adds:

- `test_timecontrol_mod_suite` — a single smoke test calling `timecontrol_init` → `timecontrol_advance` once → check that `daynr` advanced by 1 and `flDayStart`/`flDayEnd` flipped correctly. Pattern follows `test_swap_mod.pf`.

Full BMI / variable-registry expansion is out of scope (Phase 2 of the BMI arc).

## Verification gates

Per the per-task regression discipline (memory `feedback_per_task_regression_gate.md`), every implementer subagent ends with:

1. `pixi run build-linux`
2. `pixi run test-pfunit`
3. `pixi run check-fast`

Plus, after the cutover task, `pixi run check-full` (all 5 cases). Byte-for-byte parity is non-negotiable — the cleanup is supposed to be a pure refactor.

## Out of scope (deferred follow-on work)

- **`use variables` wholesale retirement**: the legacy globals import shrinks (`only:` clause and minus `flzerointr`/`flzerocumu`) but doesn't disappear. The remaining ~9 symbols (`dtmin`, `dtmax`, `flMacroPore`, `FlDecMpRat`, `dt_SSDI_event`, `flMaxIterTime`, `msteps`, `project`, output-config switches) still come from globals.
- **TimeControl unit test expansion**: only one smoke test is added. Comprehensive per-procedure unit tests are deferred — the existing regression cases provide coverage.
- **`-std=f2018` migration of `state/timecontrol_state.f90`**: the state record is already standards-compliant but lives in `swap_legacy` for build-graph reasons. Promoting it to `swap_modern` is a separate cleanup.

## Scope and effort

Approximately 13 implementer tasks, in this order:

1. Pre-flight baseline (mechanical).
2. Override audit (instrument, run check-full, capture findings, revert instrumentation).
3. Add `flZeroIntr` / `flZeroCumu` fields to `timecontrol_state_t` (state schema change; build still uses the bare globals — readers untouched in this task).
4. Create `timecontrol_mod.f90` skeleton with 7 empty procedures + failing pFUnit smoke test. New module joins `modern_sources` in `meson.build`.
5. Migrate case (1) → `timecontrol_init`; strip its dual-writes; writes go to `state%timecontrol%flZeroIntr/flZeroCumu` AND the bare globals (interim parity).
6. Migrate case (2) → `timecontrol_advance`; same dual-write-during-transition strategy.
7. Migrate case (3) → `timecontrol_reduce_dt`.
8. Migrate case (9) → `timecontrol_day_end`.
9. Migrate `IterTime` cases (1/2/3) → three named procedures.
10. Cutover: update the 7 callers in `swap_mod`; delete legacy `TimeControl` / `IterTime` from the build.
11. Migrate the 8 reader files to `state%timecontrol%flZeroIntr/flZeroCumu`; stop writing the bare globals from `timecontrol_mod` (drop the transitional dual-writes).
12. Retire the bare globals: delete `flzerointr` / `flzerocumu` from `variables.f90:117-118` and `initialize.f90:28-29`.
13. Final verification (`check-full` + retirement grep + arc-complete marker).

Total estimated commit count: 17–21 (TDD + fix-loop cadence similar to SS-DRV, slightly larger due to the flag migration).

## References

- SS-DRV Phase 1 spec (`docs/superpowers/specs/2026-05-12-driver-modernization-design.md`) — the pattern being applied.
- ADR 0041 — TimeControl state migration (the prerequisite that just shipped).
- Memory `feedback_per_task_regression_gate.md` — verification discipline.
- Memory `feedback_verify_before_committing.md` — pFUnit alone misses global-default regressions.
- `src/core/timecontrol.f90:117` — the implicit override that needs the audit.
- `src/core/timecontrol.f90:54-114` — the existing associate block (aliases preserved per-procedure in the new module).
