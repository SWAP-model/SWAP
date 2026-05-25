---
title: "SS-TCM Override Audit — TimeControl line 117 instrumentation findings"
date: 2026-05-13
status: final
relates-to: 2026-05-13-timecontrol-modernization-design.md
---

# Findings

## Method

Instrumented `src/core/timecontrol.f90:116-117` with a `write(*,...)` before the
`itask = 3` redirect. Ran `pixi run check-full` (5 regression cases, ~50 s).
Reverted the instrumentation. Inspected the captured log at `/tmp/tcm-override-audit.log`.

## Total TC-OVERRIDE firings during check-full

**Count: 0**

The override (`if (itask.eq.2 .and. (fldecdt .or. fldecdtmin)) itask = 3`) never
fired across any of the 5 check-full regression cases.

## Per-case breakdown

N/A — count is zero; no entries in the log.

## First few entries

```
(none)
```

## Caller-side gating in swap_mod

The two relevant call sites in `swap_run_step` (`src/core/swap_mod.f90`):

1. `call TimeControl(3, state)` — line 318, inside the inner `do while (fldtreduce)`
   loop, explicitly gated by:
   ```fortran
   if (fldecdt .or. (flMacroPore .and. FlDecMpRat)) then
      call SoilWaterStateVar(2, state)
      call TimeControl(3, state)
      fldtreduce = .true.
   end if
   ```
   This call is **always** made when a dt-reduction is needed; the gate already
   handles both `fldecdt` and macropore-rate triggers.

2. `call TimeControl(2, state)` — line 337, **ungated**, at the end of every
   timestep after the `do while` loop has exited.

## Analysis

**Outcome A — Never fires.**

The caller (`swap_run_step`) already routes all dt-reduction scenarios explicitly
through `TimeControl(3, ...)` before the ungated `TimeControl(2, ...)` is reached.
By the time `task=2` is called on line 337, `fldecdt` and `fldecdtmin` are always
`false` — the inner loop would not have exited otherwise. The override at line 117
is therefore dead code in all exercised paths and confirmed to never trigger during
the full regression suite.

## Recommendation

**Delete the override at Task 10 cutover.**

Do NOT add an internal redirect in `timecontrol_advance`. The caller-side gating in
`swap_run_step` is fully sufficient:

- `timecontrol_reduce_dt` will be called explicitly when `fldecdt .or. (flMacroPore
  .and. FlDecMpRat)` is true (caller gate, swap_mod.f90:316).
- `timecontrol_advance` will be called unconditionally at the end of the timestep,
  by which point neither flag is set.

No internal guard is needed in `timecontrol_advance`.
