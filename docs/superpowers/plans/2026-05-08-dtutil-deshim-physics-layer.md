# dtutil De-Shim from Physics Layer Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove all dtutil-shape calls from the physics layer (`src/crop/`, `src/soil/`, `src/drainage/`, `src/atmosphere/`, `src/macropore/`), replacing them with native helpers in non-physics modules. `dtutil.f90` itself stays for the I/O layer.

**Architecture:** Two new exports — `format_iso_date` in a new `src/utils/date_format_mod.f90`, `index_in_sorted_int` added to existing `src/utils/array_utils.f90`. Physics files swap their dtutil calls for these. Output strings must be byte-identical to existing dtdpst output so check-full baselines stay green.

**Tech Stack:** Fortran 2008, gfortran, meson+ninja+pixi, pFUnit, check-full byte-identical regression.

**Spec:** `docs/superpowers/specs/2026-05-08-dtutil-deshim-physics-layer-design.md`

---

## File Structure

**Created:**
- `src/utils/date_format_mod.f90` — `format_iso_date(t1900, with_time)` function
- `tests/unit/utils/test_date_format.pf` — pFUnit tests
- `tests/unit/utils/test_index_in_sorted_int.pf` — pFUnit tests (or merge into existing array_utils test if there is one)
- `docs/adr/0029-dtutil-deshim-physics-layer.md`

**Modified:**
- `src/utils/array_utils.f90` — add `index_in_sorted_int` function
- `src/soil/soilhydraulics.f90` — replace 4 `dtdpst` calls
- `src/soil/waterbalance.f90` — replace 1 `dtdpst` call
- `src/crop/cropgrowth.f90` — replace 1 `dtdpst` call + 1 `ifindi` call (and remove its `integer ifindi` declaration)
- `src/crop/irrigation.f90` — delete 4 commented-out `dtdpar`/`dtardp` lines
- `src/atmosphere/meteoday.f90` — replace 2 `dtdpst` calls
- `src/drainage/surfacewater.f90` — replace 1 `dtdpst` call
- `tests/unit/testSuites.inc` — register the two new test suites
- `meson.build` (or wherever utils source files are listed) — add `date_format_mod.f90` to the build
- `docs/adr/index.md` — add 0029 entry

**No changes to:**
- `src/core/dtutil.f90` (stays for I/O layer)
- I/O-layer files (`src/io/*`, `src/utils/surfacewaterutils.f90`, `src/core/timecontrol.f90`, `src/core/swap.f90`)
- Any other code

---

### Task 1: Create date_format_mod with TDD

This is the new helper that replaces all physics-layer `dtdpst` calls. Output strings must be byte-identical to dtutil's `dtdpst` for the same inputs.

**Files:**
- Create: `src/utils/date_format_mod.f90`
- Create: `tests/unit/utils/test_date_format.pf`
- Modify: `tests/unit/testSuites.inc` (add `ADD_TEST_SUITE(test_date_format_suite)`)
- Modify: meson source list (find where `src/utils/array_utils.f90` is listed; add `date_format_mod.f90` next to it)

- [ ] **Step 1: Write the failing pFUnit tests**

Create `tests/unit/utils/test_date_format.pf`:

```fortran
! Tests for src/utils/date_format_mod.f90 — native ISO date formatting
! that replaces dtdpst calls in physics-layer code. Output strings
! must be byte-identical to dtutil's dtdpst with the matching format
! string, so check-full byte-identical regression baselines stay green.

@test
subroutine test_format_iso_date_no_time()
   use funit
   use iso_fortran_env, only: real64
   use date_format_mod, only: format_iso_date
   character(len=:), allocatable :: s

   ! 2003-06-20 noon. dtdpst('year-month-day', t1900, ...) of the same
   ! input produces '2003-06-20'.
   s = format_iso_date(37791.5_real64, .false.)
   @assertEqual('2003-06-20', s)
end subroutine test_format_iso_date_no_time

@test
subroutine test_format_iso_date_with_time()
   use funit
   use iso_fortran_env, only: real64
   use date_format_mod, only: format_iso_date
   character(len=:), allocatable :: s

   ! 2003-06-20 12:00:00. dtdpst('year-month-day,hour:minute:seconds', t1900, ...)
   ! of the same input produces '2003-06-20,12:00:00'.
   s = format_iso_date(37791.5_real64, .true.)
   @assertEqual('2003-06-20,12:00:00', s)
end subroutine test_format_iso_date_with_time

@test
subroutine test_format_iso_date_leap_day()
   use funit
   use iso_fortran_env, only: real64
   use date_format_mod, only: format_iso_date
   character(len=:), allocatable :: s

   ! 2000-02-29 (Y2K leap day): t1900 = 36585.0
   s = format_iso_date(36585.0_real64, .false.)
   @assertEqual('2000-02-29', s)
end subroutine test_format_iso_date_leap_day

@test
subroutine test_format_iso_date_year_boundary()
   use funit
   use iso_fortran_env, only: real64
   use date_format_mod, only: format_iso_date
   character(len=:), allocatable :: s

   ! 1999-12-31 23:59:59 → just before Y2K rollover.
   ! t1900 = 36525 + 0.99998842592... but to get exact second use
   ! 36525 + (23*3600 + 59*60 + 59) / 86400 d
   s = format_iso_date(36525.0_real64 + 86399.0_real64/86400.0_real64, .true.)
   @assertEqual('1999-12-31,23:59:59', s)
end subroutine test_format_iso_date_year_boundary

@test
subroutine test_format_iso_date_epoch()
   use funit
   use iso_fortran_env, only: real64
   use date_format_mod, only: format_iso_date
   character(len=:), allocatable :: s

   ! t1900 = 1.0 must yield 1900-01-01 (per dtutil.f90 docstring:
   ! "DPDTTM counts days since 1900-01-01 00:00 (1900-01-01 = 1.0)").
   s = format_iso_date(1.0_real64, .false.)
   @assertEqual('1900-01-01', s)
end subroutine test_format_iso_date_epoch
```

Then add `ADD_TEST_SUITE(test_date_format_suite)` to `tests/unit/testSuites.inc` (suite name follows the pFUnit-generated `<file_basename>_suite` convention).

- [ ] **Step 2: Run tests — confirm FAIL**

Run: `pixi run -e test test-pfunit 2>&1 | tail -20`
Expected: 5 failures with "module date_format_mod not found" or similar (the module doesn't exist yet). All previously-passing tests still pass.

- [ ] **Step 3: Implement date_format_mod**

Create `src/utils/date_format_mod.f90`:

```fortran
!> @file date_format_mod.f90
!! Native ISO date formatter — replaces dtutil's dtdpst calls in the
!! physics layer. dtutil itself stays for I/O-layer use; this module
!! is the boundary helper that physics calls instead.
!!
!! Date arithmetic mirrors dtutil.f90's DPDTTM convention: t1900
!! counts days since 1900-01-01 00:00 (so 1900-01-01 = 1.0). The
!! output is byte-identical to dtutil's dtdpst for the corresponding
!! TTutil format strings — required so check-full regression
!! baselines pass byte-for-byte.
!!
!! See ADR 0024 (the de-shim future-direction note) and ADR 0029
!! (this arc's decision record).
module date_format_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   public :: format_iso_date

contains

   !> Format a t1900 timestamp as an ISO date string.
   !!
   !!   with_time = .false. → 'YYYY-MM-DD'           (10 chars)
   !!   with_time = .true.  → 'YYYY-MM-DD,HH:MM:SS'  (19 chars)
   !!
   !! Output is byte-identical to dtdpst(format, t1900, str) where
   !! format is 'year-month-day' or 'year-month-day,hour:minute:seconds'.
   function format_iso_date(t1900, with_time) result(out)
      real(real64), intent(in) :: t1900
      logical,      intent(in) :: with_time
      character(len=:), allocatable :: out
      integer :: datea(6)
      real(real64) :: fsec

      ! Reuse dtutil's date-arithmetic to guarantee byte-identical output.
      ! dtutil.f90 exposes DTDPAR(DPDTTM, DATEA, FSEC) → fills DATEA.
      call dtdpar(t1900, datea, fsec)

      if (with_time) then
         allocate(character(len=19) :: out)
         write(out, '(i4.4,"-",i2.2,"-",i2.2,",",i2.2,":",i2.2,":",i2.2)') &
              datea(1), datea(2), datea(3), datea(4), datea(5), datea(6)
      else
         allocate(character(len=10) :: out)
         write(out, '(i4.4,"-",i2.2,"-",i2.2)') &
              datea(1), datea(2), datea(3)
      end if
   end function format_iso_date

end module date_format_mod
```

Note the call to `dtdpar` — that's allowed because `date_format_mod` is in `src/utils/`, not the physics layer. The boundary rule per ADR 0024 is: physics doesn't call dtutil. Helpers in `src/utils/` may.

Add `src/utils/date_format_mod.f90` to the meson source list. Run `grep -rn "array_utils.f90" meson.build src/meson.build src/utils/meson.build 2>/dev/null | head` first to find the right list.

- [ ] **Step 4: Run tests — confirm PASS**

Run: `pixi run -e test test-pfunit 2>&1 | tail -20`
Expected: 5 new tests pass. Total = previous baseline + 5. Zero failures.

- [ ] **Step 5: Commit**

```bash
git add src/utils/date_format_mod.f90 tests/unit/utils/test_date_format.pf tests/unit/testSuites.inc
# Plus the meson.build that registered the new source
git status
git add <the meson file you edited>
git commit -m "$(cat <<'EOF'
feat(utils): SS-deshim Task 1 — native format_iso_date helper (ADR 0029)

Replaces dtdpst for physics-layer callers. Output is byte-identical
to dtutil's dtdpst with 'year-month-day' or
'year-month-day,hour:minute:seconds' formats so existing check-full
baselines stay green. Internally delegates date arithmetic to
dtutil's DTDPAR (allowed — date_format_mod lives in src/utils, not
the physics layer).

Spec: docs/superpowers/specs/2026-05-08-dtutil-deshim-physics-layer-design.md
EOF
)"
```

---

### Task 2: Add index_in_sorted_int to array_utils

Replaces the single `ifindi` call site in cropgrowth.f90.

**Files:**
- Modify: `src/utils/array_utils.f90` — add `index_in_sorted_int` function
- Create: `tests/unit/utils/test_index_in_sorted_int.pf`
- Modify: `tests/unit/testSuites.inc` (register `test_index_in_sorted_int_suite`)

- [ ] **Step 1: Inspect existing array_utils for style**

Run: `head -40 src/utils/array_utils.f90`
Use the visible style (module name, public list, function signature shape) for the new addition.

- [ ] **Step 2: Write failing pFUnit tests**

Create `tests/unit/utils/test_index_in_sorted_int.pf`:

```fortran
! Tests for index_in_sorted_int — replaces ifindi in cropgrowth.f90:993.
! ifindi(ILIS, ILDEC, IST, IEND, IINP) returns the index in ILIS
! (between IST and IEND, inclusive) where IINP appears, or 0 if not
! found. ILDEC is the array dimension. We mirror that semantics.

@test
subroutine test_index_in_sorted_int_hit()
   use funit
   use array_utils, only: index_in_sorted_int
   integer :: arr(5) = [1900, 1950, 2000, 2050, 2100]
   integer :: idx

   idx = index_in_sorted_int(arr, 2000, 1, 5)
   @assertEqual(3, idx)
end subroutine test_index_in_sorted_int_hit

@test
subroutine test_index_in_sorted_int_first_position()
   use funit
   use array_utils, only: index_in_sorted_int
   integer :: arr(5) = [1900, 1950, 2000, 2050, 2100]
   integer :: idx

   idx = index_in_sorted_int(arr, 1900, 1, 5)
   @assertEqual(1, idx)
end subroutine test_index_in_sorted_int_first_position

@test
subroutine test_index_in_sorted_int_last_position()
   use funit
   use array_utils, only: index_in_sorted_int
   integer :: arr(5) = [1900, 1950, 2000, 2050, 2100]
   integer :: idx

   idx = index_in_sorted_int(arr, 2100, 1, 5)
   @assertEqual(5, idx)
end subroutine test_index_in_sorted_int_last_position

@test
subroutine test_index_in_sorted_int_miss()
   use funit
   use array_utils, only: index_in_sorted_int
   integer :: arr(5) = [1900, 1950, 2000, 2050, 2100]
   integer :: idx

   idx = index_in_sorted_int(arr, 1975, 1, 5)
   @assertEqual(0, idx)
end subroutine test_index_in_sorted_int_miss

@test
subroutine test_index_in_sorted_int_subrange()
   use funit
   use array_utils, only: index_in_sorted_int
   integer :: arr(5) = [1900, 1950, 2000, 2050, 2100]
   integer :: idx

   ! Search only positions 2..4. Target 1900 is at position 1 → outside
   ! range → must return 0.
   idx = index_in_sorted_int(arr, 1900, 2, 4)
   @assertEqual(0, idx)
end subroutine test_index_in_sorted_int_subrange
```

Add `ADD_TEST_SUITE(test_index_in_sorted_int_suite)` to `tests/unit/testSuites.inc`.

- [ ] **Step 3: Run tests — confirm FAIL**

Run: `pixi run -e test test-pfunit 2>&1 | tail -20`
Expected: 5 new failures referencing missing `index_in_sorted_int`.

- [ ] **Step 4: Implement index_in_sorted_int in array_utils**

In `src/utils/array_utils.f90`, add the function (matching the existing module's style — likely public list at top, contains block):

```fortran
   !> Linear search for `target` in `values` between positions `lo`
   !! and `hi` (inclusive, 1-based). Returns the index where found,
   !! or 0 if not present in the range. Mirrors dtutil's IFINDI.
   pure function index_in_sorted_int(values, target, lo, hi) result(idx)
      integer, intent(in) :: values(:)
      integer, intent(in) :: target
      integer, intent(in) :: lo, hi
      integer :: idx
      integer :: i

      idx = 0
      do i = lo, hi
         if (values(i) == target) then
            idx = i
            return
         end if
      end do
   end function index_in_sorted_int
```

Add `index_in_sorted_int` to the module's `public ::` list.

- [ ] **Step 5: Run tests — confirm PASS**

Run: `pixi run -e test test-pfunit 2>&1 | tail -20`
Expected: zero failures. Test count = previous + 5.

- [ ] **Step 6: Commit**

```bash
git add src/utils/array_utils.f90 tests/unit/utils/test_index_in_sorted_int.pf tests/unit/testSuites.inc
git commit -m "$(cat <<'EOF'
feat(utils): SS-deshim Task 2 — index_in_sorted_int (replaces ifindi)

Adds a pure native-Fortran linear-search helper to array_utils with
ifindi-equivalent semantics. Replaces the single ifindi call in
cropgrowth.f90:993 (Task 5 below) without taking a dtutil dependency
inside crop/.

Spec: docs/superpowers/specs/2026-05-08-dtutil-deshim-physics-layer-design.md
EOF
)"
```

---

### Task 3: Replace dtdpst calls in src/soil/

Migrate `src/soil/soilhydraulics.f90` (4 calls) and `src/soil/waterbalance.f90` (1 call) to use `format_iso_date`.

**Files:**
- Modify: `src/soil/soilhydraulics.f90` lines 76-77, 392, 726, 772-773
- Modify: `src/soil/waterbalance.f90` line 170

- [ ] **Step 1: Migrate soilhydraulics.f90**

For each call site, replace:

```fortran
call dtdpst('year-month-day,hour:minute:seconds', t1900, datetime)
```

with:

```fortran
datetime = format_iso_date(t1900, .true.)
```

And replace `'year-month-day'` calls similarly with `format_iso_date(<arg>, .false.)`. Note the `t1900+1.001d0` and `t1900+0.1d0` forms — pass the same arithmetic expression as the first argument unchanged.

The four sites:
- Line 76-77: `call dtdpst('year-month-day,hour:minute:seconds', t1900, datetime)` → `datetime = format_iso_date(t1900, .true.)`
- Line 392: `call dtdpst('year-month-day', t1900+1.001d0, datetmp)` → `datetmp = format_iso_date(t1900+1.001d0, .false.)`
- Line 726: `call dtdpst('year-month-day', t1900+1.001d0, datetmp)` → `datetmp = format_iso_date(t1900+1.001d0, .false.)`
- Line 772-773: `call dtdpst('year-month-day,hour:minute:seconds', t1900, datetime)` → `datetime = format_iso_date(t1900, .true.)`

(Confirm exact arguments by reading the file before editing — line numbers are approximate.)

Add to the file's existing `use` block at the top:

```fortran
   use date_format_mod, only: format_iso_date
```

If the destination strings (`datetime`, `datetmp`) are declared as fixed-length characters, ensure their declared length is at least 19 to hold the longest output. The existing declarations are likely `character(len=19)` already (TTutil convention) — confirm, and increase if smaller.

- [ ] **Step 2: Migrate waterbalance.f90**

Line 170:

```fortran
call dtdpst('year-month-day,hour:minute:seconds', t1900, datexti)
```

becomes:

```fortran
datexti = format_iso_date(t1900, .true.)
```

Add `use date_format_mod, only: format_iso_date` to the file's `use` block.

- [ ] **Step 3: Build & run check-full**

Run: `pixi run check-full`
Expected: `Results: 5 passed, 0 failed`. All 5 cases byte-identical because `format_iso_date` produces the same strings as `dtdpst`.

- [ ] **Step 4: Run pFUnit suite**

Run: `pixi run -e test test-pfunit 2>&1 | tail -5`
Expected: zero failures, test count unchanged from Task 2.

- [ ] **Step 5: Commit**

```bash
git add src/soil/soilhydraulics.f90 src/soil/waterbalance.f90
git commit -m "$(cat <<'EOF'
refactor(soil): SS-deshim Task 3 — replace dtdpst with format_iso_date

5 dtdpst call sites in soilhydraulics.f90 (×4) and waterbalance.f90
(×1) migrated to date_format_mod%format_iso_date. Byte-identical
output preserves check-full baselines.

Spec: docs/superpowers/specs/2026-05-08-dtutil-deshim-physics-layer-design.md
EOF
)"
```

---

### Task 4: Replace dtdpst calls in src/atmosphere/ and src/drainage/

Migrate `src/atmosphere/meteoday.f90` (2 calls) and `src/drainage/surfacewater.f90` (1 call).

**Files:**
- Modify: `src/atmosphere/meteoday.f90` lines 347-349
- Modify: `src/drainage/surfacewater.f90` line 587

- [ ] **Step 1: Migrate meteoday.f90**

Two adjacent `dtdpst` calls (around lines 347-349):

```fortran
call dtdpst('year-month-day', dettime(irectotal)+0.1d0, detdate)
call dtdpst('year-month-day', t1900+0.1d0, date)
```

become:

```fortran
detdate = format_iso_date(dettime(irectotal)+0.1d0, .false.)
date    = format_iso_date(t1900+0.1d0, .false.)
```

Add `use date_format_mod, only: format_iso_date` at the top.

- [ ] **Step 2: Migrate surfacewater.f90**

Line 587 — read the actual file to see the exact line shape (it's a continuation: `call dtdpst &`). Convert the multi-line call into a single-line assignment:

```fortran
<destvar> = format_iso_date(<t1900-arg>, <with_time-bool>)
```

Confirm the format string used to know whether `with_time` is true. If `'year-month-day'` then `.false.`; if `'year-month-day,hour:minute:seconds'` then `.true.`.

Add `use date_format_mod, only: format_iso_date` at the top.

- [ ] **Step 3: Build & run check-full + pFUnit**

```bash
pixi run check-full && pixi run -e test test-pfunit 2>&1 | tail -5
```

Expected: 5/5 check-full pass, zero pFUnit failures.

- [ ] **Step 4: Commit**

```bash
git add src/atmosphere/meteoday.f90 src/drainage/surfacewater.f90
git commit -m "$(cat <<'EOF'
refactor(atmosphere,drainage): SS-deshim Task 4 — replace dtdpst with format_iso_date

3 dtdpst call sites migrated. Byte-identical output preserves
check-full baselines.

Spec: docs/superpowers/specs/2026-05-08-dtutil-deshim-physics-layer-design.md
EOF
)"
```

---

### Task 5: Replace dtdpst + ifindi in cropgrowth.f90 and clean up irrigation.f90

Last live physics-layer dtutil callers: cropgrowth.f90 (1 dtdpst + 1 ifindi) and 4 commented-out dead lines in irrigation.f90.

**Files:**
- Modify: `src/crop/cropgrowth.f90` lines 982, 993, 4630
- Modify: `src/crop/irrigation.f90` lines 113-127 (delete 4 commented-out lines)

- [ ] **Step 1: Migrate the dtdpst call in cropgrowth.f90:4630**

Replace:

```fortran
call dtdpst('year-month-day', t1900, dateGrassGrowth)
```

with:

```fortran
dateGrassGrowth = format_iso_date(t1900, .false.)
```

- [ ] **Step 2: Migrate the ifindi call in cropgrowth.f90:993**

Replace:

```fortran
indexyr = ifindi(CO2year, mayrs, 1, mayrs, iyear)
```

with:

```fortran
indexyr = index_in_sorted_int(CO2year, iyear, 1, mayrs)
```

**IMPORTANT:** read the file context — confirm the actual ifindi argument order. The signature is `ifindi(ILIS, ILDEC, IST, IEND, IINP)` per dtutil.f90 docstring. Map carefully:
- `CO2year` → `values` (`ILIS`)
- `iyear` → `target` (`IINP` — last argument in ifindi, third in our helper)
- `1` → `lo` (`IST`)
- `mayrs` → `hi` (`IEND`) — but ifindi has `mayrs` for both `ILDEC` AND `IEND`. The `ILDEC` argument is the array dimension, which is `size(CO2year)` or just the last valid index. In `index_in_sorted_int`, `values` is assumed-shape so we don't need a separate dimension argument; only pass `lo` and `hi`.

Final form:

```fortran
indexyr = index_in_sorted_int(CO2year, iyear, 1, mayrs)
```

Then delete the local declaration `integer ifindi, indexyr` at line 982 (drop the `ifindi` part — keep `indexyr` since it's still a local variable):

Before:
```fortran
integer ifindi, indexyr
```

After:
```fortran
integer indexyr
```

- [ ] **Step 3: Update use list in cropgrowth.f90**

Add to the `use` block at the appropriate location (probably top of the subroutine that contains line 993):

```fortran
   use array_utils, only: index_in_sorted_int
```

And at the location containing line 4630:

```fortran
   use date_format_mod, only: format_iso_date
```

If both calls are in the same subroutine they can share one `use` block. Otherwise add separately.

- [ ] **Step 4: Clean up irrigation.f90 dead lines**

Locate and delete lines 113-127 of `src/crop/irrigation.f90` — the 4 commented-out `!call dtdpar` and `!call dtardp` lines. Read the file context first to confirm the lines are still there and still dead.

Remove only the commented-out dtutil calls; preserve any surrounding live code or non-dtutil comments.

- [ ] **Step 5: Build & run check-full + pFUnit**

```bash
pixi run check-full && pixi run -e test test-pfunit 2>&1 | tail -5
```

Expected: 5/5 check-full pass, zero pFUnit failures.

- [ ] **Step 6: Verify zero physics-layer dtutil references**

Run:

```bash
grep -rEn "\b(dtdpst|dtardp|dtdpar|dtnow|dtleap|lowerc|upperc|addstr|words|decrea|ifindi)\b" \
    src/crop/ src/soil/ src/drainage/ src/atmosphere/ src/macropore/ \
    --include="*.f90" 2>/dev/null
```

Expected: ZERO matches. If any remain, identify the file and address it before committing.

- [ ] **Step 7: Commit**

```bash
git add src/crop/cropgrowth.f90 src/crop/irrigation.f90
git commit -m "$(cat <<'EOF'
refactor(crop): SS-deshim Task 5 — last physics-layer dtutil callers go (ADR 0029)

cropgrowth.f90: dtdpst → format_iso_date, ifindi → index_in_sorted_int.
irrigation.f90: delete 4 commented-out dtdpar/dtardp dead lines.

After this commit, zero dtutil-shape calls remain in the physics
layer (crop, soil, drainage, atmosphere, macropore). dtutil.f90
itself stays for I/O-layer callers.

Spec: docs/superpowers/specs/2026-05-08-dtutil-deshim-physics-layer-design.md
EOF
)"
```

---

### Task 6: ADR 0029 + final verification

**Files:**
- Create: `docs/adr/0029-dtutil-deshim-physics-layer.md`
- Modify: `docs/adr/index.md`

- [ ] **Step 1: Read ADR 0028 for house style**

Run: `cat docs/adr/0028-nutrients-N3-runtime-activation.md | head -30`
Match the YAML frontmatter, heading style, section structure.

- [ ] **Step 2: Author docs/adr/0029-dtutil-deshim-physics-layer.md**

Match the structure of ADR 0028. Required content:

```markdown
# ADR 0029: dtutil de-shim from physics layer

**Status:** accepted
**Date:** 2026-05-08
**Predecessor:** ADR 0024 (dtutil.f90 as TTutil-API compatibility shim)

## Context

ADR 0024 introduced `src/core/dtutil.f90` as a same-signature shim
for ~11 TTutil utility functions. The future-direction note in that
ADR called out hoisting these calls out of the physics layer so
physics subroutines receive parsed/validated inputs and don't
format dates or split strings. Inspection of the actual physics-
layer footprint (src/crop, src/soil, src/drainage, src/atmosphere,
src/macropore) found 9 live `dtdpst` calls, 1 live `ifindi` call,
and 4 commented-out dead `dtdpar`/`dtardp` lines — much smaller
than the headline "~75 call sites" figure, which is dominated by
I/O-layer callers.

## Decision

Replace all physics-layer dtutil calls with native helpers in
non-physics modules:

- `format_iso_date(t1900, with_time)` in a new
  `src/utils/date_format_mod.f90` — output is byte-identical to
  `dtdpst('year-month-day', …)` and
  `dtdpst('year-month-day,hour:minute:seconds', …)` so check-full
  baselines stay green.
- `index_in_sorted_int(values, target, lo, hi)` added to existing
  `src/utils/array_utils.f90` — pure ifindi-equivalent.

`dtutil.f90` itself stays. The I/O layer (`src/io/`,
`src/utils/surfacewaterutils.f90`, `src/core/timecontrol.f90`,
`src/core/swap.f90`) continues to use it. Removing dtutil entirely
is out of scope.

`date_format_mod` internally calls `dtdpar` (dtutil's date-
arithmetic). That is allowed — `src/utils/` is not the physics
layer, and reusing the proven arithmetic guarantees byte-identical
output.

## Consequences

- The physics layer now satisfies the architectural rule of ADR
  0024: zero dtutil-shape calls. A `grep` across `src/crop`,
  `src/soil`, `src/drainage`, `src/atmosphere`, `src/macropore`
  for any of the 11 dtutil-shape functions returns nothing.
- check-full byte-identical regression baselines unchanged.
- This is one of the prerequisites for the variables.f90
  dismantling arc (state encapsulation): physics modules now have
  one fewer cross-cutting dependency to thread through their
  signatures when they migrate to typed-state-bag arguments.
- `dtutil.f90` is now the *I/O-layer's* date/string helper, not
  a project-wide one. A future arc may rename it accordingly or
  fold it into a more focused module.

## References

- ADR 0024 — dtutil.f90 compatibility shim (the future-direction
  note this arc resolves)
- ADR 0023 — TTutil retirement
- Spec: docs/superpowers/specs/2026-05-08-dtutil-deshim-physics-layer-design.md
```

- [ ] **Step 3: Add the index entry**

In `docs/adr/index.md`, after the ADR 0028 line, append:

```markdown
- [ADR 0029 — dtutil de-shim from physics layer](0029-dtutil-deshim-physics-layer.html) — Native `format_iso_date` and `index_in_sorted_int` helpers replace `dtdpst`/`ifindi` in physics; `dtutil.f90` retained for I/O layer. Zero dtutil-shape calls in `src/crop`, `src/soil`, `src/drainage`, `src/atmosphere`, `src/macropore`. Resolves the future-direction note in ADR 0024.
```

- [ ] **Step 4: Final verification**

```bash
pixi run -e test test-pfunit 2>&1 | tail -10
pixi run check-full
grep -rEn "\b(dtdpst|dtardp|dtdpar|dtnow|dtleap|lowerc|upperc|addstr|words|decrea|ifindi)\b" \
    src/crop/ src/soil/ src/drainage/ src/atmosphere/ src/macropore/ \
    --include="*.f90" 2>/dev/null
```

Expected:
- pFUnit: zero failures, test count = baseline + 10 (5 from Task 1 + 5 from Task 2).
- check-full: `Results: 5 passed, 0 failed`.
- grep: zero matches.

- [ ] **Step 5: Commit**

```bash
git add docs/adr/0029-dtutil-deshim-physics-layer.md docs/adr/index.md
git commit -m "$(cat <<'EOF'
docs(adr): ADR 0029 — dtutil de-shim from physics layer

Records the decision and result of removing dtutil-shape calls from
src/crop, src/soil, src/drainage, src/atmosphere, src/macropore.
dtutil.f90 itself stays for I/O-layer use. Resolves the
future-direction note in ADR 0024.

Spec: docs/superpowers/specs/2026-05-08-dtutil-deshim-physics-layer-design.md
EOF
)"
```

---

## Self-Review Notes

- **Spec coverage:** Section 1 (goal) → all tasks. Section 2 (out-of-scope) — dtutil.f90 stays, I/O-layer untouched, error-reporting refactor deferred — none of the tasks violate. Section 3 (decision) → Task 1 (`format_iso_date`), Task 2 (`index_in_sorted_int`), Task 5 (commented-out cleanup). Section 4 (byte-identical output) → enforced by Task 1 tests + check-full at every migration task. Section 5 (testing) → pFUnit tests in Tasks 1–2, check-full integration in Tasks 3–6. Section 6 (ADR 0029) → Task 6.
- **TDD discipline:** Tasks 1 and 2 both follow red-green: failing test first, then implementation, then verify green. Tasks 3–5 are mechanical migrations gated by byte-identical check-full — the test-first burden was paid in Tasks 1–2 by proving the helpers produce the right strings.
- **Byte-identical output is the integration gate.** Every migration task ends with `pixi run check-full` and the commit only lands if `Results: 5 passed, 0 failed`. If any case fails, the helper output diverged from dtdpst — STOP, debug, do not commit.
- **No `tests/swap-cases/` edits.** Working tree should stay clean of fixture changes throughout.
- **Branch policy:** all commits go to `development`. Do not merge to main.
