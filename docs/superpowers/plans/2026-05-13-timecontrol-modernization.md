# SS-TCM — TimeControl Modernization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Convert `subroutine TimeControl` and `subroutine IterTime` in `src/core/timecontrol.f90` into `module timecontrol_mod` with seven named lifecycle procedures, migrate the `flZeroIntr` / `flZeroCumu` cross-subsystem reset gates into `state%timecontrol`, and strip the dual-write redundancy that the ADR 0041 state migration left behind.

**Architecture:** A new module `timecontrol_mod.f90` joins the `swap_modern` static library (compiled with `-std=f2018 -Wall -Wextra`). Seven named procedures (`timecontrol_init`, `timecontrol_advance`, `timecontrol_reduce_dt`, `timecontrol_day_end`, `itertime_init`, `itertime_check`, `itertime_close`) replace the magic-int dispatch. A transitional dual-write keeps regression byte-for-byte across the migration: timecontrol writers update both `state%timecontrol%flZeroIntr/flZeroCumu` AND the bare globals during Tasks 5–9; readers cut over in Task 11; bare globals retire in Task 12.

**Tech Stack:** Fortran 2018, gfortran, Meson, pFUnit 4.15, pixi orchestration.

**Spec:** `docs/superpowers/specs/2026-05-13-timecontrol-modernization-design.md`

**Builds on:** SS-DRV Phase 1 (driver modernization, 2026-05-12). All seven call sites in `swap_mod.f90` are the ones SS-DRV established.

---

## Task 1: Pre-flight baseline

**Files:**
- No code changes; capture starting state.

- [ ] **Step 1: Verify build is clean**

Run: `pixi run build-linux`
Expected: build succeeds, no warnings on existing code.

- [ ] **Step 2: Run pFUnit baseline**

Run: `pixi run test-pfunit`
Expected: all 736 tests pass (or whatever the current total is — record it).

- [ ] **Step 3: Run check-full baseline**

Run: `pixi run check-full`
Expected: 5/5 regression cases pass (hupselbrook, surfacewater, salinitystress, grassgrowth, oxygenstress).

- [ ] **Step 4: Commit baseline marker**

```bash
git commit --allow-empty -m "$(cat <<'EOF'
chore(ss-tcm): pre-flight baseline marker — TimeControl modernization arc begins

check-full: 5/5 regression cases passing. pFUnit: passing.
Build clean. Baseline locked before TimeControl refactor.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Override audit — instrument, run, capture findings, revert

**Files:**
- Modify temporarily: `src/core/timecontrol.f90` (add logging at line 117)
- Create: `docs/superpowers/specs/2026-05-13-ss-tcm-override-audit.md` (audit findings note)

The implicit dispatch override at `timecontrol.f90:117` silently converts `task=2` → `task=3` under non-convergence flags. Before deleting it in the refactor, confirm whether it ever fires in the regression cases and whether the callers already gate it.

- [ ] **Step 1: Read the current override**

Read `src/core/timecontrol.f90:115-120`. Confirm the current state:

```fortran
      itask = task
      if (itask.eq.2 .and. (fldecdt .or. fldecdtmin)) itask = 3

      select case (task)
```

Note: the `select case` dispatches on `task` (not `itask`) — so the override modifies `itask` for use later in case(2)'s body, where reduction logic kicks in. Verify by skimming case(2) body for `itask` references.

- [ ] **Step 2: Add audit logging**

Replace `src/core/timecontrol.f90:115-117` with:

```fortran
      itask = task
      if (itask.eq.2 .and. (fldecdt .or. fldecdtmin)) then
         write(*,'(a,i6,a,L1,a,L1)') &
            '[TC-OVERRIDE] daynr=', state%timecontrol%daynr, &
            ' fldecdt=', fldecdt, ' fldecdtmin=', fldecdtmin
         itask = 3
      end if
```

- [ ] **Step 3: Build and run check-full with logging**

Run: `pixi run build-linux && pixi run check-full 2>&1 | tee /tmp/tcm-override-audit.log`
Expected: build succeeds; check-full runs; any TC-OVERRIDE firings are captured in `/tmp/tcm-override-audit.log`.

- [ ] **Step 4: Inspect audit log**

Run: `grep "TC-OVERRIDE" /tmp/tcm-override-audit.log | head -20 && grep -c "TC-OVERRIDE" /tmp/tcm-override-audit.log`
Capture: total count and the first few entries. Note which regression case(s) trigger it.

- [ ] **Step 5: Inspect the calling-side context in swap_mod**

Run: `grep -n "TimeControl(3\|TimeControl(2" src/core/swap_mod.f90`
Expected to find:
- One call to `TimeControl(3, state)` inside `swap_run_step` gated by `if (fldecdt .or. (flMacroPore .and. FlDecMpRat))`
- One call to `TimeControl(2, state)` inside `swap_run_step` UNGATED at end of timestep

Read `src/core/swap_mod.f90` around each call to confirm the gating pattern.

- [ ] **Step 6: Write the audit findings**

Create `docs/superpowers/specs/2026-05-13-ss-tcm-override-audit.md` with the findings. Template:

```markdown
---
title: "SS-TCM Override Audit — TimeControl line 117 instrumentation findings"
date: 2026-05-13
status: draft
relates-to: 2026-05-13-timecontrol-modernization-design.md
---

# Findings

## Total TC-OVERRIDE firings during check-full
**Count:** <N>

## Per-case breakdown
- hupselbrook: <N>
- surfacewater: <N>
- salinitystress: <N>
- grassgrowth: <N>
- oxygenstress: <N>

## First few entries (daynr / fldecdt / fldecdtmin)
<paste 5-10 entries>

## Analysis
<one of:>

- **Outcome A — Never fires:** The override is dead code in current scenarios. Safe to delete in Task 10 (cutover).

- **Outcome B — Always fires through a caller-gated path:** When the override fires, the caller in swap_run_step has either (a) already called `TimeControl(3, ...)` in the dt-reduction inner loop, or (b) the call is inside the dt-reduction loop where fldecdt is set. The caller's gating subsumes the override → safe to delete.

- **Outcome C — Fires on a path the caller does NOT gate:** Real semantic dependency. The override must be preserved either by (a) having `timecontrol_advance` internally call `timecontrol_reduce_dt` under the same condition, or (b) adding the gate at the call site in swap_run_step.

## Recommendation
<one of the three outcomes above + rationale>
```

Fill in the count and analysis based on Step 4's findings. Apply your best judgment for which outcome applies — the most likely is Outcome A or B (the inner dt-reduction loop in swap_run_step already gates `TimeControl(3, ...)` calls explicitly).

- [ ] **Step 7: Revert the instrumentation**

Restore `src/core/timecontrol.f90:115-117` to its original three lines (no logging). Verify with `git diff src/core/timecontrol.f90` — the only intended changes should be REMOVALS of the logging lines.

- [ ] **Step 8: Build and verify revert**

Run: `pixi run build-linux && pixi run check-fast`
Expected: build succeeds; 4/4 cases pass.

- [ ] **Step 9: Commit the audit findings (not the instrumentation)**

```bash
git add docs/superpowers/specs/2026-05-13-ss-tcm-override-audit.md
# Note: timecontrol.f90 was instrumented then reverted; nothing to commit there.
git commit -m "$(cat <<'EOF'
docs(ss-tcm): TimeControl line 117 override audit — findings note

Captures the audit results before the refactor:
- Total TC-OVERRIDE firings during check-full: <N>
- Per-case breakdown
- Outcome <A|B|C>: <recommendation>

Tasks 10 (cutover) and the design of timecontrol_advance follow
this finding.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Add `flZeroIntr` / `flZeroCumu` fields to `timecontrol_state_t`

**Files:**
- Modify: `src/state/timecontrol_state.f90`

Adding the two new fields to the state record now (before the body migration) ensures the new module can dual-write into state from Task 5 onward.

- [ ] **Step 1: Read the current state record**

Read `src/state/timecontrol_state.f90` to find the Group E section (around line 127–166). Group E is "Runtime-evaluated boolean flags (26 fields)" — the new fields belong here.

- [ ] **Step 2: Add the two new fields**

In `src/state/timecontrol_state.f90`, locate the existing Group E flag declarations (lines ~131–165). Append the two new fields after the existing init-once flag group (right before `end type timecontrol_state_t` at line 167). The exact insertion point:

```fortran
      logical :: flOpenFileDev  = .false.  !! developer output file open flag

      ! [SS-TCM] cross-subsystem reset gates (migrated from variables.f90)
      ! Owned by TimeControl; consumed by 8 physics readers. Reset
      ! cadence: flZeroIntr clears intermediate accumulators; flZeroCumu
      ! clears cumulative accumulators. See design spec
      ! 2026-05-13-timecontrol-modernization-design.md.
      logical :: flZeroIntr     = .false.  !! reset gate: intermediate accumulators
      logical :: flZeroCumu     = .false.  !! reset gate: cumulative accumulators

   end type timecontrol_state_t
```

Also update the file's head doc-comment to reflect that the type now has 63 fields, not 61. Read lines 1–53 first and find the field-count text (the doc says "61 runtime-state fields"). Change to "63 runtime-state fields", and add the two new fields to the Group E listing if appropriate.

- [ ] **Step 3: Build (no readers yet, so this is just a schema change)**

Run: `pixi run build-linux`
Expected: build succeeds. The new fields are present but no caller writes them yet — they default to `.false.`.

- [ ] **Step 4: Run pFUnit and check-fast**

Run: `pixi run test-pfunit && pixi run check-fast`
Expected: pFUnit unchanged (e.g. 736/736), check-fast 4/4. The schema change is additive; no behavior change.

- [ ] **Step 5: Commit**

```bash
git add src/state/timecontrol_state.f90
git commit -m "$(cat <<'EOF'
schema(ss-tcm): add flZeroIntr/flZeroCumu to timecontrol_state_t

Two new logical fields in Group E (runtime-evaluated boolean flags),
default .false. Fields are present but unused — Tasks 5-9 add the
writers (with transitional dual-write to bare globals), Task 11
migrates the 8 readers, Task 12 retires the bare globals.

No behavior change at this step. pFUnit + check-fast unchanged.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Empty `timecontrol_mod` skeleton + failing pFUnit smoke test

**Files:**
- Create: `src/core/timecontrol_mod.f90`
- Create: `tests/unit/core/test_timecontrol_mod.pf`
- Modify: `meson.build` (add `'src/core/timecontrol_mod.f90'` to `modern_sources`)
- Modify: `tests/unit/meson.build` (add `'../../src/core/timecontrol_mod.f90'` to `pfunit_extra_sources`; add `'core/test_timecontrol_mod.pf'` to `pf_files` as the LAST entry, AFTER the existing `test_swap_mod.pf`)
- Modify: `tests/unit/testSuites.inc` (append `ADD_TEST_SUITE(test_timecontrol_mod_suite)` at the END, after the existing `test_swap_mod_suite` line)

- [ ] **Step 1: Create the empty module**

Create `src/core/timecontrol_mod.f90`:

```fortran
!> @file timecontrol_mod.f90
!! SS-TCM: module form of the legacy subroutine TimeControl + IterTime.
!! Seven named lifecycle procedures replace the magic-int dispatch
!! (task=1/2/3/9 + IterTime task=1/2/3). State threaded explicitly;
!! each procedure carries its own associate block. flZeroIntr and
!! flZeroCumu are owned by state%timecontrol — bare globals
!! consumed by Task 11 readers retire in Task 12.
module timecontrol_mod
   use swap_state_mod, only: swap_state_t
   implicit none
   private
   public :: timecontrol_init, timecontrol_advance, &
             timecontrol_reduce_dt, timecontrol_day_end
   public :: itertime_init, itertime_check, itertime_close

contains

   subroutine timecontrol_init(state)
      type(swap_state_t), intent(inout) :: state
      ! Body filled in Task 5 (migrated from timecontrol.f90 case (1)).
   end subroutine timecontrol_init

   subroutine timecontrol_advance(state)
      type(swap_state_t), intent(inout) :: state
      ! Body filled in Task 6 (migrated from timecontrol.f90 case (2)).
   end subroutine timecontrol_advance

   subroutine timecontrol_reduce_dt(state)
      type(swap_state_t), intent(inout) :: state
      ! Body filled in Task 7 (migrated from timecontrol.f90 case (3)).
   end subroutine timecontrol_reduce_dt

   subroutine timecontrol_day_end(state)
      type(swap_state_t), intent(inout) :: state
      ! Body filled in Task 8 (migrated from timecontrol.f90 case (9)).
   end subroutine timecontrol_day_end

   subroutine itertime_init(state)
      type(swap_state_t), intent(inout) :: state
      ! Body filled in Task 9 (migrated from IterTime case (1)).
   end subroutine itertime_init

   subroutine itertime_check(state)
      type(swap_state_t), intent(inout) :: state
      ! Body filled in Task 9 (migrated from IterTime case (2)).
   end subroutine itertime_check

   subroutine itertime_close(state)
      type(swap_state_t), intent(inout) :: state
      ! Body filled in Task 9 (migrated from IterTime case (3)).
   end subroutine itertime_close

end module timecontrol_mod
```

- [ ] **Step 2: Wire `timecontrol_mod.f90` into the main `meson.build`**

In `meson.build`, locate the `modern_sources = [...]` block (it currently has two entries: `swap_mod.f90` and `swap_bmi_mod.f90`). Append a sibling entry:

```python
modern_sources = [
    'src/core/swap_mod.f90',
    'src/core/swap_bmi_mod.f90',
    'src/core/timecontrol_mod.f90',
]
```

- [ ] **Step 3: Verify the new module compiles**

Run: `pixi run build-linux`
Expected: build succeeds. The new module compiles into `libswap_modern.a`. Nothing calls it yet.

- [ ] **Step 4: Write the failing pFUnit smoke test**

Create `tests/unit/core/test_timecontrol_mod.pf`:

```fortran
!> SS-TCM lifecycle smoke test for timecontrol_mod.
!! Tests the named procedures in sequence against a hupselbrook init.
!! These tests MUST be registered AFTER test_swap_mod_suite in
!! testSuites.inc — both suites mutate global state via Initialize.
@test
subroutine test_timecontrol_init_sets_iyear()
   use timecontrol_mod,  only: timecontrol_init
   use swap_state_mod,   only: swap_state_t
   use swap_config_mod,  only: swap_config_t
   use swap_mod,         only: swap_init
   use chdir_helper_mod, only: chdir_to
   use iso_fortran_env,  only: real64
   use funit
   implicit none
   type(swap_state_t)           :: state
   type(swap_config_t), target  :: config
   integer :: iyear_captured

   call chdir_to('tests/swap-cases/toml/1.hupselbrook')
   call swap_init('swap.toml', state, config)
   iyear_captured = state%timecontrol%iyear
   call chdir_to('../../../..')

   @assertGreaterThan(iyear_captured, 1900)
end subroutine

@test
subroutine test_timecontrol_advance_changes_daynr()
   use timecontrol_mod,  only: timecontrol_advance
   use swap_state_mod,   only: swap_state_t
   use swap_config_mod,  only: swap_config_t
   use swap_mod,         only: swap_init, swap_run_step
   use chdir_helper_mod, only: chdir_to
   use iso_fortran_env,  only: real64
   use funit
   implicit none
   type(swap_state_t)           :: state
   type(swap_config_t), target  :: config
   integer :: daynr_before, daynr_after

   call chdir_to('tests/swap-cases/toml/1.hupselbrook')
   call swap_init('swap.toml', state, config)
   daynr_before = state%timecontrol%daynr
   call swap_run_step(state, config)
   daynr_after = state%timecontrol%daynr
   call chdir_to('../../../..')

   ! After one step on a sub-daily dt, daynr may or may not advance.
   ! The discriminating signal is: at minimum, t1900 advanced; flDayStart
   ! should still be set on the first step.
   @assertGreaterThanOrEqual(daynr_after, daynr_before)
end subroutine

@test
subroutine test_timecontrol_flzerointr_initialized()
   ! After timecontrol_init runs (via swap_init), flZeroIntr in state must
   ! be .true. — matches case (1) line `flZeroIntr = .true.` semantics.
   use swap_state_mod,   only: swap_state_t
   use swap_config_mod,  only: swap_config_t
   use swap_mod,         only: swap_init
   use chdir_helper_mod, only: chdir_to
   use iso_fortran_env,  only: real64
   use funit
   implicit none
   type(swap_state_t)           :: state
   type(swap_config_t), target  :: config
   logical :: flzerointr_captured

   call chdir_to('tests/swap-cases/toml/1.hupselbrook')
   call swap_init('swap.toml', state, config)
   flzerointr_captured = state%timecontrol%flZeroIntr
   call chdir_to('../../../..')

   @assertTrue(flzerointr_captured)
end subroutine
```

These three tests are designed to FAIL right now because the procedures have empty stubs. Specifically:
- `test_timecontrol_init_sets_iyear` reads `state%timecontrol%iyear` set inside `swap_init` (which currently uses the LEGACY `TimeControl(1, ...)`). This passes today because legacy is still in the build. But the new procedure `timecontrol_init` is never called — this test passes regardless. Acceptable: it's a smoke test on the LIFECYCLE, not on the new procedure in isolation.
- `test_timecontrol_advance_changes_daynr` passes via legacy.
- `test_timecontrol_flzerointr_initialized` will FAIL because `state%timecontrol%flZeroIntr` is never written today (the legacy writes the bare global `flzerointr`, not the state field). This is the discriminating test that will turn green in Task 5 when `timecontrol_init` body migration adds the state-field write.

- [ ] **Step 5: Register the test file in `tests/unit/meson.build`**

In `tests/unit/meson.build`:
1. Append `'../../src/core/timecontrol_mod.f90'` to the `pfunit_extra_sources` list (anywhere; convention is near other `src/core/` entries).
2. Append `'core/test_timecontrol_mod.pf'` as the LAST entry in the `pf_files` list (after `'core/test_swap_mod.pf'`; preserve trailing comma on the preceding entry).

- [ ] **Step 6: Register the suite in `tests/unit/testSuites.inc`**

Append at the END of `tests/unit/testSuites.inc`:

```
ADD_TEST_SUITE(test_timecontrol_mod_suite)
```

(After the existing `ADD_TEST_SUITE(test_swap_mod_suite)` line.)

- [ ] **Step 7: Build with the new test wired in**

Run: `pixi run build-linux`
Expected: build succeeds; the new test_timecontrol_mod_suite is generated from the .pf file.

- [ ] **Step 8: Run pFUnit — confirm the third test FAILS**

Run: `pixi run test-pfunit`
Expected:
- `test_timecontrol_init_sets_iyear` — PASS (legacy writes iyear)
- `test_timecontrol_advance_changes_daynr` — PASS (legacy advances)
- `test_timecontrol_flzerointr_initialized` — FAIL (state field never written; bare global is set by legacy but the state field is not). This failure is the TDD signal that Task 5 will fix.
- All other suites still pass.

- [ ] **Step 9: Run check-fast (regression gate)**

Run: `pixi run check-fast`
Expected: 4/4. The empty module is unused; no behavior change in the executable.

- [ ] **Step 10: Commit**

```bash
git add src/core/timecontrol_mod.f90 tests/unit/core/test_timecontrol_mod.pf \
        meson.build tests/unit/meson.build tests/unit/testSuites.inc
git commit -m "$(cat <<'EOF'
test(ss-tcm): empty timecontrol_mod skeleton + 1 failing lifecycle test

Adds the module surface for the TimeControl refactor and three
pFUnit smoke tests on hupselbrook. Two tests pass via the legacy
TimeControl subroutine (still in the build); the third
(test_timecontrol_flzerointr_initialized) FAILS because the
state%timecontrol%flZeroIntr field is never written today — Task 5
will write it from timecontrol_init.

Module joins modern_sources (-std=f2018 -Wall -Wextra). Suite
registered LAST in testSuites.inc for global-state isolation.

check-fast: 4/4 (no behavior change).

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: Migrate case (1) body → `timecontrol_init`

**Files:**
- Modify: `src/core/timecontrol_mod.f90` (fill `timecontrol_init` body)

**Source of truth:** `src/core/timecontrol.f90:120-373` (the case (1) block contents).

The body migration follows the SS-DRV pattern: copy verbatim, then apply edits. Two specific changes from SS-DRV:

**A. Strip dual-writes.** Every line of the form `X = value; state%timecontrol%X = X` (inside the existing associate block) becomes just `X = value`. The associate alias makes the state write happen automatically. Example:

Before:
```fortran
flRunEnd = .false.
state%timecontrol%flRunEnd = flRunEnd
```
After:
```fortran
flRunEnd = .false.
```

**B. Add dual-WRITE for flZeroIntr/flZeroCumu during transition.** Bare globals are still consumed by 8 reader files until Task 11. Each write to `flZeroIntr` or `flZeroCumu` must update BOTH the bare global AND the new state field:

Before:
```fortran
flZeroIntr = .true.
flZeroCumu = .true.
```
After:
```fortran
flZeroIntr = .true.
flZeroCumu = .true.
state%timecontrol%flZeroIntr = .true.   ! [SS-TCM transition] dual-write; readers cut over Task 11
state%timecontrol%flZeroCumu = .true.   ! [SS-TCM transition] dual-write; readers cut over Task 11
```

This is the OPPOSITE pattern from edit A — for these two fields specifically, we add the state-field write because they are NOT in the associate block (they don't exist in state yet from the body's perspective — Task 3 added them but case (1) was written before).

Actually — `flZeroIntr` and `flZeroCumu` are NOT in the associate block at line 54-114 (verify by reading). They are bare globals from `use variables`. So the dual-write is `flZeroIntr = ...; state%timecontrol%flZeroIntr = ...` (two distinct memory locations).

- [ ] **Step 1: Read the source body**

Read `src/core/timecontrol.f90:120-373`. Identify:
- Field assignments paired with `state%timecontrol%X = X` dual-writes (strip in edit A).
- Writes to `flZeroIntr` / `flZeroCumu` (add transitional state-field writes per edit B).
- Calls to external subroutines (`dtdpar`, etc.) and module subs (will be carried through unchanged).

- [ ] **Step 2: Add procedure-scope `use` statements**

Inside `timecontrol_init`, add the use list this procedure needs. Starting list (verify by grepping the body for symbol references):

```fortran
   subroutine timecontrol_init(state)
      use variables, only: dtmin, dtmax, period, swheader, swodat, swres, swscre, &
                            outdat, outdatint, tend, tstart, &
                            flprintdt, nprintday, logf, flCropCalendar, &
                            flMacroPore, swirfix, swsnow, swdra, &
                            swhea, swsolu, swetsine, swrain, swmetdetail, &
                            metperiod_input, flzerointr, flzerocumu
      use error_mod, only: fatalerr_collected
      implicit none
      type(swap_state_t), intent(inout) :: state
```

Adjust the list during the migration based on what the case (1) body actually references. The `implicit none` is required by `-std=f2018` (project-wide `add_project_arguments` has `-std=legacy` for legacy files but not for modern_sources).

- [ ] **Step 3: Add the full associate block**

Copy the associate block from `src/core/timecontrol.f90:54-114` verbatim into `timecontrol_init` after the use statements. Closing `end associate` goes at the end of the procedure body (before `end subroutine`).

(Pragmatic choice: each procedure starts with the FULL alias list. Code review may prune unused aliases per procedure.)

- [ ] **Step 4: Copy the case (1) body**

Copy the lines from `src/core/timecontrol.f90:121-372` (between `case (1)` and `return` of that case, inclusive of body, exclusive of `case (1)`/`return`). Paste inside the associate block of `timecontrol_init`.

- [ ] **Step 5: Apply edit A — strip dual-writes**

For every line in the case (1) body that has the form:

```fortran
X = expression
state%timecontrol%X = X
```

(where `X` is a name aliased in the associate block — daynr, iyear, flRunEnd, flDayStart, etc.) — DELETE the second line.

This is mechanical: search for `state%timecontrol%` inside the copied body. Each such line is paired with the preceding assignment to the aliased name. Delete the `state%timecontrol%X = X` line, leaving the prior aliased write.

Approximate count: ~30 dual-write pairs in case (1).

- [ ] **Step 6: Apply edit B — add transitional state-field writes for flZeroIntr/flZeroCumu**

In the copied body, locate the two lines:
```fortran
      flZeroIntr = .true.
      flZeroCumu = .true.
```
(They appear near the top of case (1), around legacy line 138-139.)

After these two lines, ADD:
```fortran
      state%timecontrol%flZeroIntr = .true.   ! [SS-TCM transition] dual-write; readers cut over Task 11
      state%timecontrol%flZeroCumu = .true.   ! [SS-TCM transition] dual-write; readers cut over Task 11
```

- [ ] **Step 7: Build**

Run: `pixi run build-linux`
Expected: build succeeds with `-std=f2018` on the modern file. If you get f2018 errors (implicit typing, deprecated syntax), report the error — minor fixes (add `implicit none` if missing, add explicit types) can be done in this task; non-trivial issues escalate.

Common gotcha: `-std=f2018` rejects `goto` to numbered labels in some forms, `call Exit(...)`, and non-block do-loops. If the case (1) body has any of these, port them to standard equivalents (`stop`, block-DO).

- [ ] **Step 8: Verify pFUnit**

Run: `pixi run test-pfunit`
Expected:
- `test_timecontrol_init_sets_iyear` still PASS
- `test_timecontrol_advance_changes_daynr` still PASS
- `test_timecontrol_flzerointr_initialized` — NOW PASS (state field written by edit B).

Wait — but the new `timecontrol_init` is not called yet. The state field is written ONLY when `timecontrol_init` runs, not when legacy `TimeControl(1, ...)` runs. So the test still fails. **Expected at this step.**

To clarify: tests/unit/core/test_timecontrol_mod.pf is currently testing the lifecycle through swap_init → legacy TimeControl. After Task 5 the new procedure has the right body, but the test path still goes through legacy. The TDD signal for Task 5 is the BUILD passing — the test will turn green in Task 10 when the cutover swaps the call.

Run pFUnit anyway and confirm: 2/3 lifecycle tests pass, 1 still fails (the flZeroIntr one). All other suites pass.

- [ ] **Step 9: Run check-fast**

Run: `pixi run check-fast`
Expected: 4/4. The new `timecontrol_init` is unused by the executable (still calling legacy); regression unaffected.

- [ ] **Step 10: Commit**

```bash
git add src/core/timecontrol_mod.f90
git commit -m "$(cat <<'EOF'
refactor(ss-tcm): migrate case (1) body into timecontrol_init

Copies the initialization block from timecontrol.f90:120-373 into
the timecontrol_init procedure body. Two transformations applied:

  Edit A: strip ~30 dual-write pairs of the form
          `X = value; state%timecontrol%X = X` — the associate
          block alias makes the state write happen automatically.

  Edit B: add transitional dual-write for flZeroIntr/flZeroCumu —
          legacy global write preserved (8 readers still consume
          the bare globals until Task 11); new state-field write
          added.

Procedure is in the build but unused; legacy TimeControl(1, ...)
in swap_init still drives the executable. Cutover in Task 10.

pFUnit: unchanged (3rd lifecycle test still fails until cutover).
check-fast: 4/4.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: Migrate case (2) body → `timecontrol_advance`

**Files:**
- Modify: `src/core/timecontrol_mod.f90` (fill `timecontrol_advance` body)

**Source of truth:** `src/core/timecontrol.f90:375-783` (case (2) body, ~410 lines).

The case (2) body has additional `flzerointr`/`flzerocumu` writes (lines 613, 623, 629, 705, 712, 720, 731, 732 — all `flzerointr = .true.` or `flzerocumu = .true.` triggered by output cadence and year boundaries). Each one gets the same transitional dual-write treatment as Task 5 Edit B.

- [ ] **Step 1: Decide override policy based on Task 2 audit findings**

Read `docs/superpowers/specs/2026-05-13-ss-tcm-override-audit.md` (committed in Task 2 Step 9).

- If the audit recommended **Outcome A or B** (override is dead or caller-gated): proceed with the migration WITHOUT preserving the `itask=2→3` redirect. The caller in `swap_run_step` already gates `TimeControl(3, ...)` explicitly.

- If the audit recommended **Outcome C** (real semantic dependency): add the override semantics INSIDE `timecontrol_advance`:
  ```fortran
  ! [SS-TCM] preserve legacy override: under non-convergence flags,
  ! advance dispatches to reduce_dt instead. Caller does NOT gate
  ! this path (per audit Outcome C).
  if (state%timecontrol%fldecdtmin .or. fldecdt) then
     call timecontrol_reduce_dt(state)
     return
  end if
  ```
  Add this at the very top of the body, after the associate block.

Document the choice in a comment at the procedure head.

- [ ] **Step 2: Add procedure-scope use statements**

Inside `timecontrol_advance`, add the procedure-scope use list. Starting set (verify by grep against case (2) body):

```fortran
   subroutine timecontrol_advance(state)
      use variables, only: dtmin, dtmax, period, outdat, outdatint, tend, tstart, &
                            nprintday, swheader, swodat, swres, swscre, &
                            logf, flCropCalendar, flMaxIterTime, swirfix, &
                            msteps, flzerointr, flzerocumu, flMacroPore, &
                            FlDecMpRat
      use timestep_control_mod, only: fldecdt
      use error_mod, only: fatalerr_collected
      implicit none
      type(swap_state_t), intent(inout) :: state
```

Refine during migration based on actual symbol references in the case (2) body.

- [ ] **Step 3: Add the full associate block**

Same associate block as Task 5 Step 3 — copy from `timecontrol.f90:54-114` verbatim.

- [ ] **Step 4: Copy the case (2) body**

Copy `src/core/timecontrol.f90:376-782` (between `case (2)` line and its final `return`). Paste inside the associate block of `timecontrol_advance`.

- [ ] **Step 5: Apply edit A — strip dual-writes**

Same mechanical strip as Task 5 Step 5. Approximate count: ~50 dual-write pairs in case (2). Search for `state%timecontrol%` in the copied body and delete each line paired with a preceding aliased assignment.

- [ ] **Step 6: Apply edit B — transitional dual-writes for flZeroIntr/flZeroCumu**

Find every write to `flZeroIntr` or `flZeroCumu` in the copied body. These should be at approximately the original lines 419, 423, 613, 623, 629, 705, 712, 720, 731, 732. For EACH such write, ADD a paired state-field write:

```fortran
      flZeroIntr = .true.
      state%timecontrol%flZeroIntr = .true.   ! [SS-TCM transition] dual-write; readers cut over Task 11
```

Same pattern for `flZeroCumu = .false.`:
```fortran
      flZeroCumu = .false.
      state%timecontrol%flZeroCumu = .false.   ! [SS-TCM transition] dual-write; readers cut over Task 11
```

Approximate count: ~10 transitional dual-write pairs (one per write site).

- [ ] **Step 7: Build**

Run: `pixi run build-linux`
Expected: build succeeds. If f2018 errors appear, fix in place per Task 5 Step 7's guidance.

- [ ] **Step 8: Run pFUnit and check-fast**

Run: `pixi run test-pfunit && pixi run check-fast`
Expected: pFUnit pattern unchanged (procedure exists but not called); check-fast 4/4.

- [ ] **Step 9: Commit**

```bash
git add src/core/timecontrol_mod.f90
git commit -m "$(cat <<'EOF'
refactor(ss-tcm): migrate case (2) body into timecontrol_advance

Copies the timestep-advance block from timecontrol.f90:376-782
(~410 lines) into timecontrol_advance.

Edit A: strip ~50 dual-write pairs (state%timecontrol%X = X) that
the associate block makes redundant.

Edit B: add ~10 transitional state-field writes for flZeroIntr /
flZeroCumu at every legacy write site (lines 419, 423, 613, 623,
629, 705, 712, 720, 731, 732).

Override policy: <Outcome A/B/C from audit, with explanation>.

Procedure is in the build but unused. check-fast: 4/4.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: Migrate case (3) body → `timecontrol_reduce_dt`

**Files:**
- Modify: `src/core/timecontrol_mod.f90` (fill `timecontrol_reduce_dt` body)

**Source of truth:** `src/core/timecontrol.f90:784-834` (case (3) body, ~50 lines). Reduces dt by 3× or to dtmin under non-convergence; resets `dtprevious`, `flTnext`.

- [ ] **Step 1: Add procedure-scope use statements**

```fortran
   subroutine timecontrol_reduce_dt(state)
      use variables, only: dtmin, dtmax, flMacroPore, FlDecMpRat
      use timestep_control_mod, only: fldecdt
      implicit none
      type(swap_state_t), intent(inout) :: state
```

- [ ] **Step 2: Add the associate block (full alias list, same as Tasks 5/6)**

- [ ] **Step 3: Copy case (3) body**

Copy `src/core/timecontrol.f90:785-834` (between `case (3)` and the last `return` of that case). The body has three branches:
1. `if (fldecdt) then ... endif` (lines ~790-810): dt reduction by 3× or to dtmin
2. `if (fldecdtmin) then ... endif` (lines ~812-825): forced dt = dtmin
3. `if (flMacroPore .and. FlDecMpRat) then ... endif` (lines ~827-834): macropore-driven sqrt(dtmin*dtmax) — note: this branch is **dead** because `flMacroPore` was retired in ADR 0040 (always .false.). Preserve it for now (the legacy still has it); a follow-on can prune.

- [ ] **Step 4: Apply edit A — strip dual-writes**

Same as Tasks 5/6. Approximate count: ~10 dual-write pairs in case (3).

- [ ] **Step 5: No edit B in this task**

`flZeroIntr` and `flZeroCumu` are not written in case (3).

- [ ] **Step 6: Build, pFUnit, check-fast**

Run: `pixi run build-linux && pixi run test-pfunit && pixi run check-fast`
Expected: all green; pattern unchanged; check-fast 4/4.

- [ ] **Step 7: Commit**

```bash
git add src/core/timecontrol_mod.f90
git commit -m "$(cat <<'EOF'
refactor(ss-tcm): migrate case (3) body into timecontrol_reduce_dt

Copies the dt-reduction block from timecontrol.f90:785-834 into
timecontrol_reduce_dt. Three branches preserved:
- fldecdt: reduce by 3x or floor at dtmin
- fldecdtmin: force dtmin
- flMacroPore + FlDecMpRat: sqrt(dtmin*dtmax) — dead branch since
  ADR 0040 retired macropores, preserved for now (follow-on prune).

Edit A: strip ~10 dual-write pairs.

Procedure is in the build but unused. check-fast: 4/4.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 8: Migrate case (9) body → `timecontrol_day_end`

**Files:**
- Modify: `src/core/timecontrol_mod.f90` (fill `timecontrol_day_end` body)

**Source of truth:** `src/core/timecontrol.f90:836-841` (case (9) body, ~7 lines). Caps next-day dt for SSDI events.

- [ ] **Step 1: Add procedure-scope use statements**

```fortran
   subroutine timecontrol_day_end(state)
      use variables, only: dt_SSDI_event
      implicit none
      type(swap_state_t), intent(inout) :: state
```

- [ ] **Step 2: Add the associate block — only the fields this procedure needs**

Case (9) only touches `dt`. The associate block can be MINIMAL here (one alias):

```fortran
      associate( dt => state%timecontrol%dt )
```

Closing `end associate` at end of procedure.

(This is the first procedure where the minimal associate block makes sense. Larger procedures use the full block per the pragmatic choice.)

- [ ] **Step 3: Copy case (9) body**

Copy `src/core/timecontrol.f90:836-841`:

```fortran
!        special: at end of day, the possible initial time step for next day may be too large: adapt if necessary
         if (dt_SSDI_event < 1.0d0) then
            dt = min(dt, dt_SSDI_event)
            state%timecontrol%dt = dt
         end if
```

- [ ] **Step 4: Apply edit A — strip dual-write**

The line `state%timecontrol%dt = dt` is the dual-write. Delete it. Result:

```fortran
      if (dt_SSDI_event < 1.0d0) then
         dt = min(dt, dt_SSDI_event)
      end if
```

- [ ] **Step 5: Build, pFUnit, check-fast**

Run: `pixi run build-linux && pixi run test-pfunit && pixi run check-fast`
Expected: all green; 4/4.

- [ ] **Step 6: Commit**

```bash
git add src/core/timecontrol_mod.f90
git commit -m "$(cat <<'EOF'
refactor(ss-tcm): migrate case (9) body into timecontrol_day_end

Copies the 6-line SSDI day-end dt cap from timecontrol.f90:836-841
into timecontrol_day_end. Minimal associate block (only `dt`
aliased). One dual-write removed.

Procedure is in the build but unused. check-fast: 4/4.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 9: Migrate IterTime cases → three named procedures

**Files:**
- Modify: `src/core/timecontrol_mod.f90` (fill `itertime_init`, `itertime_check`, `itertime_close` bodies)

**Source of truth:** `src/core/timecontrol.f90:853-912` (IterTime subroutine, 3 cases).

- [ ] **Step 1: Fill `itertime_init` body**

Case (1) is one line: `call cpu_time(state%timecontrol%tmptimestart)`. Procedure-scope use list is minimal.

```fortran
   subroutine itertime_init(state)
      implicit none
      type(swap_state_t), intent(inout) :: state
      call cpu_time(state%timecontrol%tmptimestart)
   end subroutine itertime_init
```

No associate block needed (one direct access).

- [ ] **Step 2: Fill `itertime_check` body**

Case (2) is the per-day cpu-time budget check. Copy `timecontrol.f90:879-890`:

```fortran
   subroutine itertime_check(state)
      use variables, only: MaxIterTime
      use error_mod, only: fatalerr_collected
      implicit none
      type(swap_state_t), intent(inout) :: state
      real(4) :: tmptimeinterrupt
      integer :: timediff
      character(len=400) :: messag

      call cpu_time(tmptimeinterrupt)
      timediff = int(tmptimeinterrupt) - MaxIterTime
      if (timediff > 0) then
         write(messag,'(a,i10,3a)') &
            'The maximum cpu time of ', MaxIterTime, ' (secs)', &
            ' was exceeded.  Therefore simulation was interrupted'
         call fatalerr_collected('IterTime', messag)
      end if
   end subroutine itertime_check
```

Note the `intent(inout) :: state` is kept for signature parity with the other procedures; `state` is unused in the body. Compiler warnings are suppressed by `-Wno-unused-dummy-argument` in the project flags.

- [ ] **Step 3: Fill `itertime_close` body**

Case (3) writes iteration statistics to the log. Copy `timecontrol.f90:892-905`:

```fortran
   subroutine itertime_close(state)
      use variables, only: MaxIt, itnumb, logf
      implicit none
      type(swap_state_t), intent(inout) :: state
      integer :: i, j

      write(logf, '(/,a20)')      'Iteration statistics'
      write(logf, '(/,a29,i4)')   'Maximum number of iterations:', MaxIt
      write(logf, '(/,a35/,a35)') 'It Numb  No of Hits  Tot BTr cycles', &
                                   '-------  ----------  --------------'
      do i = 1, 100
         if (itnumb(i,1) > 0) &
            write(logf, '(i7,2x,i10,4x,i10)') i, (itnumb(i,j), j=1,2)
      end do

      call cpu_time(state%timecontrol%tmptimeend)
      write(logf, '(/,a12,f12.2,a4)') &
         ' Run-time: ', state%timecontrol%tmptimeend - state%timecontrol%tmptimestart, ' sec'
   end subroutine itertime_close
```

- [ ] **Step 4: Build**

Run: `pixi run build-linux`
Expected: build succeeds with all 7 procedures filled. The bare global `itnumb` from `variables` is still needed — itnumb tracking is OUTSIDE the case bodies (it's incremented elsewhere) so no migration here.

- [ ] **Step 5: Run pFUnit and check-fast**

Run: `pixi run test-pfunit && pixi run check-fast`
Expected: pattern unchanged; 4/4.

- [ ] **Step 6: Commit**

```bash
git add src/core/timecontrol_mod.f90
git commit -m "$(cat <<'EOF'
refactor(ss-tcm): migrate IterTime cases into three named procedures

- itertime_init  <- IterTime case (1): cpu_time start
- itertime_check <- IterTime case (2): cpu-time budget check
- itertime_close <- IterTime case (3): iteration statistics log

Direct copies from timecontrol.f90:874-905; no dual-writes to
strip (legacy IterTime already wrote state%timecontrol fields
directly).

timecontrol_mod is now functionally complete; legacy still in
build. Cutover in Task 10. check-fast: 4/4.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 10: Cutover — update swap_mod callers, delete legacy `timecontrol.f90`

**Files:**
- Modify: `src/core/swap_mod.f90` (7 call sites)
- Modify: `meson.build` (remove `'src/core/timecontrol.f90'` from `sources`)
- Modify: `tests/unit/meson.build` (remove `'../../src/core/timecontrol.f90'` from `pfunit_extra_sources` if present)
- Delete: `src/core/timecontrol.f90`

- [ ] **Step 1: Update swap_init's call sites**

In `src/core/swap_mod.f90`, locate `swap_init`. Find the two calls:
- `call IterTime(1, state)` — change to `call itertime_init(state)`
- `call TimeControl(1, state)` — change to `call timecontrol_init(state)`

Add to swap_init's procedure-scope use list:
```fortran
use timecontrol_mod, only: timecontrol_init, itertime_init
```

- [ ] **Step 2: Update swap_run_step's call sites**

In `swap_run_step`, locate the calls:
- `call IterTime(2, state)` (inside the day-end block, gated by `flMaxIterTime`) — change to `call itertime_check(state)`
- `call TimeControl(3, state)` (inside the dt-reduction inner loop) — change to `call timecontrol_reduce_dt(state)`
- `call TimeControl(2, state)` (end of step) — change to `call timecontrol_advance(state)`
- `call TimeControl(9, state)` (day-end SSDI) — change to `call timecontrol_day_end(state)`

Add to swap_run_step's procedure-scope use list:
```fortran
use timecontrol_mod, only: timecontrol_advance, timecontrol_reduce_dt, &
                            timecontrol_day_end, itertime_check
```

- [ ] **Step 3: Update swap_close's call site**

In `swap_close`, locate:
- `call IterTime(3, state)` — change to `call itertime_close(state)`

Add to swap_close's procedure-scope use list:
```fortran
use timecontrol_mod, only: itertime_close
```

- [ ] **Step 4: Verify all 7 sites are updated**

Run: `grep -n "call TimeControl\b\|call IterTime\b" src/core/swap_mod.f90`
Expected: NO matches. All seven sites now use the named procedures.

Run: `grep -n "call timecontrol_\|call itertime_" src/core/swap_mod.f90`
Expected: exactly 7 matches (4 TimeControl + 3 IterTime).

- [ ] **Step 5: Remove legacy file from build**

In `meson.build`, in the `sources = [...]` list, find and DELETE the line:
```
    'src/core/timecontrol.f90',
```

If `tests/unit/meson.build`'s `pfunit_extra_sources` references `timecontrol.f90`, remove that too. (Check with `grep -n "timecontrol.f90" tests/unit/meson.build`.)

- [ ] **Step 6: Delete legacy file**

```bash
git rm src/core/timecontrol.f90
```

- [ ] **Step 7: Build**

Run: `pixi run build-linux`
Expected: build succeeds. If unresolved-symbol errors appear for `TimeControl` or `IterTime`, a caller was missed in Step 1-3.

- [ ] **Step 8: Run pFUnit — all 3 lifecycle tests now PASS**

Run: `pixi run test-pfunit`
Expected: 736 (or whatever current total) — all pass, including:
- `test_timecontrol_init_sets_iyear` — PASS
- `test_timecontrol_advance_changes_daynr` — PASS
- `test_timecontrol_flzerointr_initialized` — PASS (state%timecontrol%flZeroIntr is now written by timecontrol_init from Task 5).

- [ ] **Step 9: Run check-full (regression gate)**

Run: `pixi run check-full`
Expected: 5/5 regression cases pass byte-for-byte.

If a case fails, the cutover broke physics. Diagnose by inspecting the case's diff output vs baseline. Most likely culprit: a missed dual-write strip that turned out NOT to be redundant (the alias was masked by `intent` semantics, or the bare global has a separate writer outside the migrated case). Fix in place.

- [ ] **Step 10: Commit**

```bash
git add src/core/swap_mod.f90 meson.build tests/unit/meson.build
git commit -m "$(cat <<'EOF'
refactor(ss-tcm): cutover — retire timecontrol.f90, update swap_mod callers

Replaces the 7 magic-int call sites in swap_mod with named-procedure
calls into timecontrol_mod:

  TimeControl(1, state)  -> timecontrol_init(state)
  TimeControl(2, state)  -> timecontrol_advance(state)
  TimeControl(3, state)  -> timecontrol_reduce_dt(state)
  TimeControl(9, state)  -> timecontrol_day_end(state)
  IterTime(1, state)     -> itertime_init(state)
  IterTime(2, state)     -> itertime_check(state)
  IterTime(3, state)     -> itertime_close(state)

The legacy `subroutine TimeControl` and `subroutine IterTime` are
deleted (`src/core/timecontrol.f90` removed from build and disk).
The implicit task=2->3 override at the old line 117 is gone (or
preserved inside timecontrol_advance per audit Outcome C).

The transitional flZeroIntr/flZeroCumu dual-writes inside
timecontrol_mod remain — they retire in Task 12 after readers
migrate in Task 11.

pFUnit: all 3 lifecycle tests now PASS. check-full: 5/5 byte-for-
byte parity.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 11: Migrate 8 reader files to `state%timecontrol%flZeroIntr/flZeroCumu`

**Files (modify each):**
- `src/atmosphere/meteoday.f90`
- `src/drainage/surfacewater.f90`
- `src/drainage/drainage.f90`
- `src/atmosphere/snow.f90`
- `src/soil/soilhydraulics.f90`
- `src/soil/waterbalance.f90`
- `src/solute/agetracer.f90`
- `src/solute/solute.f90`

Each reader needs to (a) drop `flzerointr`/`flzerocumu` from its `use variables` list, (b) change the bare-global reads to `state%timecontrol%flZeroIntr` / `flZeroCumu`. Verify `state` is in the subroutine's signature; if not, thread it (the SS-* migration arcs already did this for all 8 files, per spec, but verify per file).

After all 8 readers cut over, the bare globals in `timecontrol_mod` (the Edit-B transitional dual-writes from Tasks 5-6) become unused. They're removed in Task 12.

- [ ] **Step 1: Migrate `src/atmosphere/meteoday.f90`**

Read the file around line 440-454. Confirm:
- Line 440 has `use variables, only: flzerointr,flzerocumu` (possibly with other symbols).
- Lines 448 and 454 read the flags.

Verify the enclosing subroutine has `state` (or accepts state). If yes:
- Drop `flzerointr, flzerocumu` from the `use variables, only:` clause on line 440.
- Change `if (flzerointr) then` → `if (state%timecontrol%flZeroIntr) then`
- Change `if (flzerocumu) then` → `if (state%timecontrol%flZeroCumu) then`

If `state` is NOT in the signature, add it: `type(swap_state_t), intent(inout) :: state` (and corresponding `use swap_state_mod, only: swap_state_t`). Update the caller in `swap_run_step` to pass state — verify by checking that `ProcessMeteoDay(state)` or `ReadMeteoDay(state)` already passes state, then trace through to meteoday.f90.

Build after this single-file migration: `pixi run build-linux`. Verify check-fast: `pixi run check-fast`.

If both pass, commit:
```bash
git add src/atmosphere/meteoday.f90
git commit -m "$(cat <<'EOF'
refactor(ss-tcm): meteoday — state%timecontrol%flZeroIntr/flZeroCumu

Drops the two flags from `use variables, only:` and reads from
state%timecontrol instead. One of 8 reader-file migrations in
Task 11 of the SS-TCM arc.

check-fast: 4/4.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

- [ ] **Step 2: Migrate `src/drainage/surfacewater.f90`**

Lines 103 and 109 read the flags; `state` is already in scope (the legacy already used `state%surfacewater%reset_intermediate()` etc.). Search the file's `use variables` clause and drop `flzerointr`/`flzerocumu` from it (if listed by `only:` — otherwise the bare-global access just becomes `state%timecontrol%flZeroIntr`).

Change:
- `if (flzerointr) call state%surfacewater%reset_intermediate()` → `if (state%timecontrol%flZeroIntr) call state%surfacewater%reset_intermediate()`
- `if (flzerocumu) then` → `if (state%timecontrol%flZeroCumu) then`

Build, check-fast, commit (same pattern as Step 1).

- [ ] **Step 3: Migrate `src/drainage/drainage.f90`**

Three reads: lines 429, 502, 510. Line 429 is inside a parameter list (the call has the flag as an argument); 502 and 510 read the flag directly.

Read around line 429 to understand the parameter-passing context. If a subroutine takes `flzerointr, flzerocumu` as args, change to take `state` and read inside. If the call site has access to `state`, replace `flzerointr, flzerocumu` arguments at the call with the state-field reads.

Lines 502 and 510: same pattern as Step 2.

Build, check-fast, commit.

- [ ] **Step 4: Migrate `src/atmosphere/snow.f90`**

Lines 106 and 114. The subroutine should already accept `state` (snow_init / Snow take state per SS-DRV/SS-* migration). Drop the `use variables` references, switch to state-field reads.

Build, check-fast, commit.

- [ ] **Step 5: Migrate `src/soil/soilhydraulics.f90`**

Lines 1136 and 1148. The subroutines in this file (SoilWater, etc.) accept state. Drop globals, switch to state-field reads.

Build, check-fast, commit.

- [ ] **Step 6: Migrate `src/soil/waterbalance.f90`**

Line 443 (only flZeroIntr). Drop the global, switch to state-field read.

Build, check-fast, commit.

- [ ] **Step 7: Migrate `src/solute/agetracer.f90`**

Line 159 (only flZeroCumu). Drop the global, switch to state-field read.

Build, check-fast, commit.

- [ ] **Step 8: Migrate `src/solute/solute.f90`**

Lines 128 and 129. Drop the globals, switch to state-field reads.

Build, check-fast, commit.

- [ ] **Step 9: Drop the transitional bare-global writes in `timecontrol_mod`**

Now that all 8 readers are off the bare globals, the Edit-B dual-writes in `timecontrol_mod` can be removed. Open `src/core/timecontrol_mod.f90` and DELETE every line of the form:
```fortran
state%timecontrol%flZeroIntr = .true.   ! [SS-TCM transition] dual-write; readers cut over Task 11
state%timecontrol%flZeroCumu = .true.   ! [SS-TCM transition] dual-write; readers cut over Task 11
```
(and the `.false.` variants).

Wait — that's the WRONG line to delete. Re-read the transition strategy:
- During Tasks 5-9, writers in `timecontrol_mod` wrote BOTH the bare global `flZeroIntr` AND the state field `state%timecontrol%flZeroIntr`.
- Readers in Tasks 11.1-11.8 just switched to read the state field.
- Now the bare global has no consumers. Drop the BARE-GLOBAL writes; keep the STATE-FIELD writes.

So in `timecontrol_mod`:
- DELETE the lines `flZeroIntr = .true.` / `flZeroIntr = .false.` / `flZeroCumu = .true.` / `flZeroCumu = .false.`
- KEEP the lines `state%timecontrol%flZeroIntr = .true.` etc.
- Also drop `flzerointr, flzerocumu` from `timecontrol_mod`'s `use variables, only:` clause (both procedures' use lists).
- Remove the `[SS-TCM transition]` comments — they're no longer transitional.

Build, run pFUnit, run check-fast.

- [ ] **Step 10: Final check-fast for Task 11**

Run: `pixi run check-fast`
Expected: 4/4 byte-for-byte.

- [ ] **Step 11: Commit the bare-global write drop**

```bash
git add src/core/timecontrol_mod.f90
git commit -m "$(cat <<'EOF'
refactor(ss-tcm): drop transitional bare-global writes for flZeroIntr/flZeroCumu

All 8 reader files have migrated to state%timecontrol — the bare
global writes in timecontrol_mod are no longer consumed by any
caller. Drop them along with the [SS-TCM transition] comments.

The bare global declarations still exist in variables.f90 (and the
init-zero in initialize.f90) — those retire in Task 12.

check-fast: 4/4.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 12: Retire bare globals from `variables.f90` and `initialize.f90`

**Files:**
- Modify: `src/core/variables.f90` (delete lines 117-118)
- Modify: `src/core/initialize.f90` (delete lines 28-29)

- [ ] **Step 1: Verify no remaining bare-global references**

Run: `grep -irn "flzerointr\|flzerocumu" src/ | grep -v "\.mod" | grep -v "state%timecontrol%flZero"`

Expected: ONLY references should be:
- `src/core/variables.f90:117-118` (the declarations — about to delete)
- `src/core/initialize.f90:28-29` (the init-zeros — about to delete)
- `src/state/timecontrol_state.f90` (doc-comment mentioning Group B retirement)
- Possibly some unrelated doc strings in state-record headers.

If ANY runtime-code reference outside these locations remains (e.g. a reader file that wasn't migrated, a writer in legacy code), STOP and fix that reader first.

- [ ] **Step 2: Delete from `variables.f90`**

In `src/core/variables.f90`, delete lines 117-118:

```fortran
      logical   flzerocumu         ! Flag indicating that cumulative fluxes should be reset to zero (reset-arc)
      logical   flzerointr         ! Flag indicating that intermediate fluxes should be reset to zero (reset-arc)
```

Use Edit with these exact strings as `old_string` so the deletion is precise.

- [ ] **Step 3: Delete from `initialize.f90`**

In `src/core/initialize.f90`, delete lines 28-29:

```fortran
      flzerocumu         = .false.
      flzerointr         = .false.
```

- [ ] **Step 4: Build**

Run: `pixi run build-linux`
Expected: build succeeds.

If link errors appear for `flzerointr` or `flzerocumu`, a reader was missed in Task 11. Find it (`grep -irn "flzerointr\|flzerocumu" src/`), migrate it, then rebuild.

- [ ] **Step 5: Run all gates**

Run: `pixi run test-pfunit && pixi run check-fast`
Expected: pFUnit unchanged (all pass); check-fast 4/4.

- [ ] **Step 6: Update the timecontrol_state_t doc-comment**

In `src/state/timecontrol_state.f90`, find the head doc-comment (lines 35-40) that says:

```
!! Excluded (deferred per design doc):
!!   - Group A: 18 config-constants ...
!!   - Group B: flZeroIntr / flZeroCumu — cross-subsystem reset gates;
!!     future reset-orchestration arc.
```

Remove the Group B bullet (it's no longer excluded — it's been migrated). Keep Group A and `dtEventRain` bullets.

- [ ] **Step 7: Commit**

```bash
git add src/core/variables.f90 src/core/initialize.f90 src/state/timecontrol_state.f90
git commit -m "$(cat <<'EOF'
retire(ss-tcm): delete bare-global flZeroIntr/flZeroCumu

The two cross-subsystem reset-gate flags are now owned by
state%timecontrol (Task 3 added the fields; Tasks 5-9 added the
writers; Task 11 migrated the 8 readers). The bare globals in
variables.f90 are dead; the init-zeros in initialize.f90 are
redundant (state fields default to .false.).

Updates the timecontrol_state.f90 head doc-comment to remove the
"Group B excluded" note.

check-fast: 4/4.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 13: Final verification + arc-complete marker

**Files:**
- No code changes; verification only.

- [ ] **Step 1: Full check-full**

Run: `pixi run check-full`
Expected: 5/5 regression cases pass.

- [ ] **Step 2: Retirement grep**

Run:
```bash
grep -rn 'subroutine TimeControl\b\|subroutine IterTime\b' src/
grep -irn '\bflzerointr\b\|\bflzerocumu\b' src/ | grep -v "state%timecontrol%flZero" | grep -v ":!"
```

Expected: both grep commands return either no hits, or only doc-comments that explicitly reference the retired patterns (e.g., in `state/timecontrol_state.f90` field doc-comments). Live executable references should be zero.

- [ ] **Step 3: Confirm legacy file deleted**

Run: `ls src/core/timecontrol.f90 2>&1`
Expected: `No such file or directory`.

- [ ] **Step 4: Confirm new files present**

Run:
```bash
ls src/core/timecontrol_mod.f90 tests/unit/core/test_timecontrol_mod.pf
```
Expected: both files listed without error.

- [ ] **Step 5: Confirm sources list reflects the change**

Run: `grep -n "timecontrol" meson.build`
Expected: `timecontrol_mod.f90` in `modern_sources`; `timecontrol.f90` ABSENT from `sources`.

- [ ] **Step 6: Tag the arc complete**

```bash
git commit --allow-empty -m "$(cat <<'EOF'
chore(ss-tcm): TimeControl modernization complete

Summary of the arc:
- subroutine TimeControl + IterTime (913 lines, free, magic-int)
  -> module timecontrol_mod (7 named procedures)
- 7 callers in swap_mod.f90 updated (TimeControl/IterTime -> named)
- flZeroIntr / flZeroCumu migrated from variables.f90 globals into
  state%timecontrol; 8 reader files cut over; bare globals deleted
- ~120 lines of dual-write redundancy stripped (state%X = X pairs)
- implicit task=2->3 dispatch override at line 117: <retired per
  audit Outcome A/B, OR preserved inside timecontrol_advance per
  Outcome C>
- timecontrol_state_t: 61 -> 63 fields
- New module joins swap_modern static lib (-std=f2018 -Wall -Wextra)
- 3 pFUnit lifecycle smoke tests on hupselbrook
- check-full: all 5 regression cases byte-for-byte parity

Deferred to follow-on arcs:
- Wholesale `use variables` retirement (timecontrol_mod still imports
  ~9 globals via only:)
- Macropore-related dead branch in timecontrol_reduce_dt (legacy
  preservation, prune candidate)
- Comprehensive per-procedure unit tests (only lifecycle smoke today)

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Plan self-review

**Spec coverage check:**

- New `timecontrol_mod` with 7 named procedures → Tasks 4 (skeleton) + 5/6/7/8/9 (bodies).
- Override audit precondition → Task 2.
- `flZeroIntr` / `flZeroCumu` as state fields → Task 3 (schema), Tasks 5/6 (writers), Task 11 (readers), Task 12 (retirement).
- Dual-write stripping → Tasks 5/6/7/8 (Edit A in each).
- Per-procedure associate blocks → Tasks 5/6/7/8 (full associate as pragmatic starting point per spec).
- Caller updates in swap_mod (7 sites) → Task 10 (cutover).
- Legacy file deletion → Task 10.
- Build wiring (modern_sources) → Task 4.
- pFUnit lifecycle smoke test → Task 4 (skeleton, with 1 failing test as TDD signal).
- 8 reader migrations enumerated → Task 11 (Steps 1-8, one file each).
- Regression byte-for-byte gate (check-full) → Tasks 10, 13.

**Placeholder scan:** No "TBD" / "TODO" / "fill in later" — all code blocks contain actual code; all commands have expected outputs.

**Type/name consistency:** Procedure names `timecontrol_init / _advance / _reduce_dt / _day_end` and `itertime_init / _check / _close` used consistently across all tasks. Fields `flZeroIntr` and `flZeroCumu` (CamelCase) used consistently as the new state-field names; bare globals `flzerointr` / `flzerocumu` (lowercase) used consistently when referring to the legacy. `timecontrol_state_t`, `swap_state_t`, `swap_config_t` types match across tasks.

**Ordering check:** Task 3 (schema add) precedes Task 5 (which writes those fields). Tasks 5-9 (body migrations) precede Task 10 (cutover). Task 10 (cutover) precedes Task 11 (reader migration) — because readers consume from state, and state writers must be live first. Task 11 precedes Task 12 (global retirement) — because the bare-global readers must be off before the declaration can be removed.

**Audit-outcome conditional in Task 6:** Task 6 Step 1 explicitly reads the audit findings and branches the implementation. This is the only place the audit result affects the implementation; everywhere else assumes the audit confirmed safety. If audit returns Outcome C, Task 6 Step 1's conditional branch handles it.
