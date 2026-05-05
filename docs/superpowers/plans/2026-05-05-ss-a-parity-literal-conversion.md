# SS-A — Per-parity-suite literal-value capture + helper-file drop

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Convert each parity test from `legacy_global == config%X`
assertions to `<literal> == config%X` assertions, where `<literal>`
is the value the legacy reader currently writes. After SS-A, no
parity test invokes `readswap` or any legacy helper; the parity-test
helper files are deleted.

**Sub-spec parent:**
`docs/superpowers/specs/2026-05-05-legacy-readers-physical-deletion-design.md`
(SS-A).

**Tech stack:** Fortran 2008, pFUnit, pixi, meson.

**Effort:** M.

---

## File structure

| File | Change |
|---|---|
| `tests/unit/io/toml/test_hupselbrook_parity.pf` | **Modify:** convert ~12 assertions per test, ~6 tests in suite. Drop `parity_helpers_mod` import; replace `load_both_for_hupselbrook` with `load_swap_config` + `validate` + `finalize`. |
| `tests/unit/io/toml/test_grassgrowth_parity.pf` | Same shape, case 2 inputs. |
| `tests/unit/io/toml/test_oxygenstress_parity.pf` | Same shape, case 4 inputs. |
| `tests/unit/io/toml/test_salinitystress_parity.pf` | Same shape, case 5 inputs. Slight quirk: this suite also tests SWSOLU=1 + ldis_array path. |
| `tests/unit/io/toml/test_surfacewater_parity.pf` | Same shape, case 6 inputs. Quirk: drives `rddre` via `surfacewater_init` (SS-3 era). |
| `tests/unit/io/toml/test_legacy_crop_helper.pf` | **Decide & modify** — see Task 7. This suite tests the helper itself, so the conversion is different. |
| `tests/unit/io/toml/parity_helpers.f90` | **Delete** at end of SS-A (after all parity files converted). |
| `tests/unit/io/toml/legacy_crop_helper.f90` | **Delete** with parity_helpers (no callers after Task 7). |
| `tests/unit/io/toml/readswap_stubs.f90` | **Delete** with parity_helpers (no callers). |
| `tests/unit/meson.build` | **Modify:** drop the three deleted files from `pfunit_extra_sources`. |

---

## Conventions

1. **Capture the literal once, in a one-off run.** For each parity
   test: instrument with a `write(0,*) 'EXPECT', '<field>', '=',
   <legacy_global>` line, run the test in TAP mode, grep stderr for
   the `EXPECT` lines, then delete the instrumentation. Captured
   values become `real(real64), parameter` declarations at the top of
   each test.
2. **One commit per parity suite.** After each conversion: build →
   `pixi run -e test test-pfunit -f <suite>` → check-full → commit.
3. **Spot-check captured literals.** Each captured value should be
   non-zero (most physics fields are positive) and should match the
   value documented in the `READING NOTE` comments inside the test
   (e.g. "kdif = 0.60 from grassd.crp"). If a captured literal is 0
   or contradicts the comment, investigate before committing.
4. **Helper-file deletion is the LAST commit** of SS-A — only after
   all six parity files compile cleanly without `parity_helpers_mod`.
5. **Preserve test names and docstrings.** Conversion changes the
   data source, not the test intent.

---

## Task 1: Capture instrumentation utility

**Files:**
- New: `tests/unit/io/toml/capture_helpers.f90` (a tiny module with a
  `print_capture(field_name, value)` helper that writes to stderr in
  the standard format).
- Modify: `tests/unit/meson.build` — add the new file to
  `pfunit_extra_sources`.

- [ ] **Step 1: write the capture helper**

```fortran
! tests/unit/io/toml/capture_helpers.f90
!
! Temporary instrumentation used during SS-A literal-value capture.
! Each parity test gets a few `call print_capture(...)` lines, runs
! once to emit literals to stderr, then both the lines and this
! module are deleted. The module exists only for the duration of SS-A.
module capture_helpers_mod
   use iso_fortran_env, only: real64, error_unit
   implicit none
   private
   public :: print_capture_real, print_capture_int, print_capture_str
contains
   subroutine print_capture_real(name, val)
      character(len=*), intent(in) :: name
      real(real64),     intent(in) :: val
      write(error_unit, '("EXPECT ", a, " = ", es24.16)') name, val
   end subroutine
   subroutine print_capture_int(name, val)
      character(len=*), intent(in) :: name
      integer,          intent(in) :: val
      write(error_unit, '("EXPECT ", a, " = ", i0)') name, val
   end subroutine
   subroutine print_capture_str(name, val)
      character(len=*), intent(in) :: name
      character(len=*), intent(in) :: val
      write(error_unit, '("EXPECT ", a, " = ", a)') name, trim(val)
   end subroutine
end module
```

- [ ] **Step 2: register in meson.build**

Add `'io/toml/capture_helpers.f90'` to `pfunit_extra_sources` in
`tests/unit/meson.build` (anywhere in the list).

- [ ] **Step 3: build to verify**

Run: `pixi run -e test build-linux 2>&1 | tail -5`
Expected: build succeeds.

(No commit yet — this helper module is temporary scaffold; commit it
together with the first parity conversion in Task 2.)

---

## Task 2: Convert `test_hupselbrook_parity.pf` (smallest case)

**Files:**
- Modify: `tests/unit/io/toml/test_hupselbrook_parity.pf`
- Modify: `tests/unit/meson.build` (already done in Task 1)

- [ ] **Step 1: instrument the suite**

For each `@assertEqual(legacy_global, config%X, ...)` in the file,
insert a `call print_capture_*('legacy.<name>', legacy_global)` line
BEFORE the assertion. Add `use capture_helpers_mod, only:
print_capture_real, print_capture_int` to each test's `use` clauses.

- [ ] **Step 2: capture the literals**

Run: `./builddir/tests/unit/unit-swap-tests --tap /tmp/cap.tap -f
test_hupselbrook_parity 2> /tmp/cap.stderr < /dev/null`

Then: `grep '^EXPECT' /tmp/cap.stderr | sort -u > /tmp/literals.txt`

Inspect `/tmp/literals.txt`. Each line is `EXPECT legacy.<name> =
<value>`. These are the ground-truth literals.

- [ ] **Step 3: rewrite assertions to literals**

For each `@assertEqual(legacy_global, config%X, ...)`:
1. Look up the captured literal in `/tmp/literals.txt`.
2. Replace the assertion with `@assertEqual(<literal>, config%X,
   ...)` (preserve the tolerance).
3. Drop the now-unused `legacy_global` from the `use variables, only:
   ...` clause.

Drop the `call load_both_for_hupselbrook(config, errors)` line.
Replace with:

```fortran
call load_swap_config(  &
   'tests/swap-cases/toml/1.hupselbrook/swap.toml', config, errors)
call config%validate(errors)
call config%finalize(errors)
@assertFalse(errors%has_fatals())
```

Drop the `use parity_helpers_mod, only: ...` import.

Drop the `use capture_helpers_mod, ...` import (no longer needed).

Drop the `call print_capture_*(...)` instrumentation lines (they
served their purpose).

- [ ] **Step 4: verify**

Run: `pixi run -e test build-linux 2>&1 | tail -3`
Expected: build succeeds.

Run: `./builddir/tests/unit/unit-swap-tests --tap /tmp/v.tap -f
test_hupselbrook_parity < /dev/null > /dev/null 2>&1; grep -c '^ ok'
/tmp/v.tap; grep -c '^ not ok' /tmp/v.tap`
Expected: same number of `ok` as before SS-A; zero `not ok`.

- [ ] **Step 5: commit**

```
git add tests/unit/io/toml/test_hupselbrook_parity.pf \
        tests/unit/io/toml/capture_helpers.f90 \
        tests/unit/meson.build
git commit -m "test(parity): convert hupselbrook to literal assertions (SS-A)

Replace `legacy_global == config%X` parity assertions with captured-
once literals. Test no longer invokes readswap; load_swap_config is
called directly. Adds tests/unit/io/toml/capture_helpers.f90 as a
temporary instrumentation module used during the capture pass —
deleted at the end of SS-A.
"
```

---

## Tasks 3-6: Convert grassgrowth / oxygenstress / salinitystress / surfacewater parity suites

Same shape as Task 2, applied to:

- [ ] **Task 3:** `test_grassgrowth_parity.pf`
- [ ] **Task 4:** `test_oxygenstress_parity.pf`
- [ ] **Task 5:** `test_salinitystress_parity.pf`
- [ ] **Task 6:** `test_surfacewater_parity.pf`

For each, follow Task 2's Steps 1-5. **One commit per suite.** Verify
gates after each.

Notes:

- **Salinitystress** (Task 5) has additional state assertions for
  SWSOLU=1 / ldis_array. Capture and convert those too.
- **Surfacewater** (Task 6) drives `rddre` via the helper for the
  surface-water management block. The capture must include the
  per-period management arrays (`impend(:)`, `swman(:)`, etc.). The
  test currently asserts on `min(size(impend), size(config%impend))`;
  preserve that pattern.
- **Macroporeflow** (case 3) is already disabled in
  `testSuites.inc` — skip.

---

## Task 7: Decide on `test_legacy_crop_helper.pf`

**Files:**
- Modify or delete: `tests/unit/io/toml/test_legacy_crop_helper.pf`

- [ ] **Step 1: read the suite to understand its purpose**

`tests/unit/io/toml/test_legacy_crop_helper.pf` (4 tests) tests the
helper itself (`read_legacy_wofost`, `read_legacy_cropfixed`,
`read_legacy_grass`, `flatten_table_basic`). Its purpose was to
verify that `legacy_crop_helper.f90` correctly drove the legacy
readers under per-rotation conditions.

- [ ] **Step 2: decide**

Two options:

**Option A — delete the suite entirely.** The helper itself is
deleted in SS-A's final task; tests for a deleted helper have no
purpose. Simplest.

**Option B — convert each test to a TOML-equivalent assertion.** For
example, `test_read_legacy_wofost_populates_variables` asserts that
the legacy reader populates `kdif /= 0`; the equivalent modern test
would assert that `load_swap_config` + `cropwofost_init_from_config`
populates `kdif /= 0`. Adds incremental coverage of the per-crop init
modules.

Recommended: **Option A**. The per-crop init modules are already
covered by `test_cropwofost_init.pf`, `test_cropfixed_init.pf`, and
`test_cropgrass_init.pf`. Re-testing them under "what the legacy
helper would have done" adds no signal.

- [ ] **Step 3: execute the chosen option**

If **Option A**: `git rm tests/unit/io/toml/test_legacy_crop_helper.pf`
and remove the corresponding `pf_files` entry from
`tests/unit/meson.build` and the `ADD_TEST_SUITE(test_legacy_crop_helper_suite)`
line from `tests/unit/testSuites.inc`.

- [ ] **Step 4: verify**

Build + pFUnit + check-full.

- [ ] **Step 5: commit**

```
git commit -m "test: drop test_legacy_crop_helper suite (SS-A)

The suite tested the legacy_crop_helper module which is being
deleted in SS-A's final commit. The per-crop init modules
(cropfixed_init, cropwofost_init, cropgrass_init) have their own
unit-test coverage; this suite added no independent signal."
```

---

## Task 8: Drop the helper files + capture_helpers + retire run_pfunit.sh

**Files:**
- Delete: `tests/unit/io/toml/parity_helpers.f90`
- Delete: `tests/unit/io/toml/legacy_crop_helper.f90`
- Delete: `tests/unit/io/toml/readswap_stubs.f90`
- Delete: `tests/unit/io/toml/capture_helpers.f90` (no longer
  referenced after Tasks 2-6)
- Modify: `tests/unit/meson.build` (drop the four files from
  `pfunit_extra_sources`)

**NOTE:** SS-A drops the helper files only. The `run_pfunit.sh` wrapper
and the `meson.build` `test()` invocation that uses it stay until SS-B
(after `readswap.f90` itself is gone — only then is per-suite isolation
truly unnecessary).

- [ ] **Step 1: confirm no remaining callers**

```bash
grep -rn "use parity_helpers_mod\|use legacy_crop_helper_mod\|use capture_helpers_mod\|readswap_stubs" tests/ src/ --include='*.f90' --include='*.pf' 2>/dev/null
```

Expected: zero matches. If any remain, return to the corresponding
parity-conversion task and fix.

- [ ] **Step 2: delete the files**

```bash
git rm tests/unit/io/toml/parity_helpers.f90 \
       tests/unit/io/toml/legacy_crop_helper.f90 \
       tests/unit/io/toml/readswap_stubs.f90 \
       tests/unit/io/toml/capture_helpers.f90
```

- [ ] **Step 3: update meson.build**

Edit `tests/unit/meson.build`: remove the four lines from
`pfunit_extra_sources`.

- [ ] **Step 4: verify**

```
pixi run -e test build-linux 2>&1 | tail -3
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:"
pixi run -e test check-full 2>&1 | grep -E "passed|failed"
```

Expected: all green; failure count unchanged from current baseline
(0 failures).

- [ ] **Step 5: commit**

```
git commit -m "test: drop legacy parity helper files (SS-A close)

Helpers (parity_helpers.f90, legacy_crop_helper.f90,
readswap_stubs.f90, capture_helpers.f90) had no callers after the
per-parity-suite literal-value conversion (commits TBD). Drop them
from src list and meson.build.

This closes SS-A. SS-B (code deletion of readswap.f90 + dead
case(1) blocks) is now unblocked.
"
```

---

## Definition of done (SS-A)

- All 5 parity test suites (1, 2, 4, 5, 6) assert against literals
  instead of legacy globals.
- `test_legacy_crop_helper.pf` either converted (Option B) or
  deleted (Option A).
- `parity_helpers.f90`, `legacy_crop_helper.f90`,
  `readswap_stubs.f90`, `capture_helpers.f90` deleted from
  `tests/unit/io/toml/`.
- `meson.build` no longer references these files.
- `grep -rn "call readswap\|call read_legacy_\|load_both_for_" tests/`
  returns zero hits in `*.pf` files.
- `pixi run -e test test-pfunit` exits 0 (`Ok: 1, Fail: 0`).
- `pixi run -e test check-full` exits 0 (5 passed, 0 failed).
- One commit per parity suite + one helper-file-drop commit
  (~7 commits total).
