# Meteo Typed-CSV-Records Pilot — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the untyped `state%atmosphere%metcsv_dat(:,:) / metcsv_det(:,:) / raincsv_dat(:,:)` caches with three typed table records (`meteo_daily_table_t`, `meteo_detail_table_t`, `rain_events_table_t`) defined in a new `src/io/csv/meteo_csv.f90` module. Establish the typed-CSV pattern that six more CSV families will follow.

**Architecture:** Each table type owns its schema (column names, units), a `load(path, errors)` method that wraps the existing generic `csv_reader_mod%read_csv_table`, an inline `validate` step (called from `load`, never separately), a `year_window(year, i1, i2)` slice helper, and an `is_loaded :: logical` flag. Date column stays `real(real64)` days-since-1900 (no `date_t` rework in this pilot). Consumers (`atmosphere_state_init`, `readmeteo.f90`'s `MeteoCSVYear` / `MeteoCSVDetYear` / `ReadRainEvents` / `read_meteo_from_external_buffer_year`) iterate typed rows by named field — no more column indices.

**Tech Stack:** Fortran 2008 derived types with type-bound procedures, gfortran 13+, Meson build system, pFUnit 4.15 unit tests, pixi test environment.

**Regression gate (explicit):** This refactor is **byte-identical**. Both `pixi run -e test check-fast` (4 fast cases) and `pixi run -e test check-full` (6 cases) must pass with zero output diff against the gfortran-4.2.0 baseline. The new pFUnit tests for the typed tables are mandatory and registered in `testSuites.inc` (see memory: `pFUnit tests need testSuites.inc registration` — a `.pf` not listed there runs dark).

**State-schema rebuild note:** Task 7 changes the field layout of `atmosphere_state_t`. Per the saved memory `state-schema-changes-require-clean-rebuild`, the implementer **must** `rm -rf builddir` before rebuilding after Task 7's edit — incremental Meson builds do not propagate `.mod` dependency changes across the `swap_modern↔swap_legacy` boundary.

---

## File Structure

**New files:**
- `src/io/csv/meteo_csv.f90` — three typed table types + their methods. ~280 lines.
- `tests/unit/io/csv/test_meteo_csv.pf` — pFUnit tests for all three types. ~250 lines.
- `tests/unit/io/fixtures/csv_meteo_daily.csv` — 5-row fixture covering normal days + winter wrap. ~10 lines.
- `tests/unit/io/fixtures/csv_meteo_detail.csv` — 3-row sub-daily fixture. ~6 lines.
- `tests/unit/io/fixtures/csv_meteo_rain_events.csv` — 4-row fixture spanning two years. ~7 lines.
- `tests/unit/io/fixtures/csv_meteo_daily_bad_temp.csv` — fixture with `tmin > tmax` for validation test. ~4 lines.

**Modified files:**
- `src/state/atmosphere_state.f90` — fields replaced (`metcsv_dat`/`nmetcsv` etc. → `meteo`/`meteo_detail`/`rain_events` typed tables); `atmosphere_state_init` shrinks (three CSV-load blocks collapse to three method calls).
- `src/io/readmeteo.f90` — `MeteoCSVYear`, `MeteoCSVDetYear`, `ReadRainEvents`, `read_meteo_from_external_buffer_year` reach typed fields by name instead of column index.
- `tests/unit/meson.build` — register new `.pf` file and `meteo_csv.f90` source.
- `tests/unit/testSuites.inc` — register new test suite.
- `tests/unit/state/test_atmosphere_state.pf` — update any assertion that read `metcsv_dat` directly (none expected per the file as of this plan, but verify).

**Touch-but-don't-rewrite:**
- `src/state/atmosphere_state.f90:241-265` (the per-year extracted arrays `arad/atmn/atmx/...`) — these are populated *from* the cache by `MeteoCSVYear` and stay. Only the cache representation changes.

---

## Task 0: Pre-flight & worktree

**Files:** none — environment only.

- [ ] **Step 1: Create an isolated worktree.**

Run:
```bash
cd /home/zawadzkim/Code/swap
git worktree add -b meteo-typed-csv ../swap-meteo-typed-csv development
cd ../swap-meteo-typed-csv
```
Expected: new worktree at `/home/zawadzkim/Code/swap-meteo-typed-csv`, on branch `meteo-typed-csv`.

- [ ] **Step 2: Verify clean baseline regression.**

Run:
```bash
pixi run -e test check-fast
```
Expected: PASS. If it fails on baseline, stop and surface — the development branch is broken independently of this work.

- [ ] **Step 3: Snapshot a baseline output for byte-diff later.**

Run:
```bash
mkdir -p .scratch
pixi run -e test pytest tests/regression/test_output_regression.py -k hupselbrook -v > .scratch/baseline_hupselbrook.txt 2>&1
```
Expected: PASS recorded. (We'll re-run the same after Task 9 and diff.)

---

## Task 1: Typed table types — declarations only

**Files:**
- Create: `src/io/csv/meteo_csv.f90`
- Test: `tests/unit/io/csv/test_meteo_csv.pf`
- Modify: `tests/unit/meson.build`
- Modify: `tests/unit/testSuites.inc`

- [ ] **Step 1: Write a failing pFUnit test that the three types exist with the expected fields.**

Create `tests/unit/io/csv/test_meteo_csv.pf` with:
```fortran
! Tests for typed CSV record tables in meteo_csv_mod.
! Pilot for the typed-CSV-records pattern; six more families to follow.

@test
subroutine test_types_exist_and_default()
   use funit
   use iso_fortran_env, only: real64
   use meteo_csv_mod, only: meteo_daily_table_t, meteo_detail_table_t, rain_events_table_t

   type(meteo_daily_table_t)  :: daily
   type(meteo_detail_table_t) :: detail
   type(rain_events_table_t)  :: rain

   ! is_loaded defaults to false; rows unallocated.
   @assertFalse(daily%is_loaded)
   @assertFalse(detail%is_loaded)
   @assertFalse(rain%is_loaded)
   @assertFalse(allocated(daily%rows))
   @assertFalse(allocated(detail%rows))
   @assertFalse(allocated(rain%rows))
end subroutine test_types_exist_and_default
```

- [ ] **Step 2: Register the new test file in the unit-tests build.**

Edit `tests/unit/meson.build`. Find the `pf_files = [` block at line 207. Add this line in the `io/` group (alphabetical near line 234):
```meson
        'io/csv/test_meteo_csv.pf',
```

Add the new source to `pfunit_extra_sources` (search the file for `pfunit_extra_sources = [`; add):
```meson
        '../../src/io/csv/meteo_csv.f90',
```

- [ ] **Step 3: Register the suite in `testSuites.inc`.**

Edit `tests/unit/testSuites.inc`. Add at the end of the io-tests group (near the `test_csv_reader_suite` line):
```
ADD_TEST_SUITE(test_meteo_csv_suite)
```

- [ ] **Step 4: Update the project meson.build include path for the new directory.**

Edit `meson.build` at the top of the file (the `src_inc = include_directories(...)` block near line 8). Add `'src/io/csv'` to the list:
```meson
src_inc = include_directories(
    'src', 'src/core', 'src/io', 'src/io/csv', 'src/io/toml', 'src/soil', 'src/atmosphere',
    'src/crop', 'src/drainage', 'src/boundary',
    'src/solute', 'src/heat', 'src/utils', 'src/config', 'src/state'
)
```

- [ ] **Step 5: Run test to verify it fails (module not found).**

Run:
```bash
rm -rf builddir
pixi run -e test meson setup builddir
pixi run -e test meson compile -C builddir 2>&1 | tail -20
```
Expected: compile failure with "Cannot open module file 'meteo_csv_mod.mod'" or similar.

- [ ] **Step 6: Implement the three types and the module.**

Create `src/io/csv/meteo_csv.f90`:
```fortran
!> Typed CSV record tables for meteorological forcing files.
!!
!! Three families:
!!   * meteo_daily_table_t   — 9-column daily meteo (ADR 0014 canonical schema)
!!   * meteo_detail_table_t  — 7-column sub-daily meteo (swmetdetail=1)
!!   * rain_events_table_t   — 2-column rain events (swrain=3)
!!
!! Each table owns its schema, a `load(path, errors)` method that wraps the
!! generic csv_reader, an inline validate step, and a year-window slice
!! helper. The `is_loaded` flag distinguishes "not requested" from "empty".
!!
!! Date column is real(real64) days-since-JD-1900 (csv_reader's convention).
!! A future arc may introduce a typed date_t; out of scope here.
module meteo_csv_mod
   use iso_fortran_env, only: real64
   use error_mod,       only: error_collection_t
   implicit none
   private

   public :: meteo_daily_table_t
   public :: meteo_detail_table_t
   public :: rain_events_table_t

   !> One daily meteo record. Units match csv_reader output (no scaling).
   type :: meteo_daily_row_t
      real(real64) :: date  = 0.0_real64  !! days since JD 1900
      real(real64) :: rad   = 0.0_real64  !! kJ m-2 d-1 (NOT yet J/m2/d)
      real(real64) :: tmin  = 0.0_real64  !! degC
      real(real64) :: tmax  = 0.0_real64  !! degC
      real(real64) :: hum   = 0.0_real64  !! kPa
      real(real64) :: wind  = 0.0_real64  !! m s-1
      real(real64) :: rain  = 0.0_real64  !! mm d-1
      real(real64) :: etref = 0.0_real64  !! mm d-1 (-99.9 = compute internally)
      real(real64) :: wet   = 0.0_real64  !! fraction (-99.9 = missing)
   end type meteo_daily_row_t

   type :: meteo_daily_table_t
      type(meteo_daily_row_t), allocatable :: rows(:)
      logical :: is_loaded = .false.
   contains
      procedure :: load        => meteo_daily_table_load
      procedure :: year_window => meteo_daily_table_year_window
   end type meteo_daily_table_t

   !> One sub-daily record. Column 1 is datetime (continuous fractional days
   !! since JD 1900), not whole-day; column 2 is the within-day record index.
   type :: meteo_detail_row_t
      real(real64) :: datetime = 0.0_real64
      integer      :: record   = 0
      real(real64) :: rad      = 0.0_real64  !! kJ m-2 (per-period total)
      real(real64) :: temp     = 0.0_real64  !! degC
      real(real64) :: hum      = 0.0_real64  !! kPa
      real(real64) :: wind     = 0.0_real64  !! m s-1
      real(real64) :: rain     = 0.0_real64  !! mm (per-period total)
   end type meteo_detail_row_t

   type :: meteo_detail_table_t
      type(meteo_detail_row_t), allocatable :: rows(:)
      logical :: is_loaded = .false.
   contains
      procedure :: load        => meteo_detail_table_load
      procedure :: year_window => meteo_detail_table_year_window
   end type meteo_detail_table_t

   !> One rain event (sub-daily rain accumulation row).
   type :: rain_event_row_t
      real(real64) :: datetime = 0.0_real64  !! days since JD 1900 (fractional)
      real(real64) :: amount   = 0.0_real64  !! mm (raw, no unit conversion)
   end type rain_event_row_t

   type :: rain_events_table_t
      type(rain_event_row_t), allocatable :: rows(:)
      logical :: is_loaded = .false.
   contains
      procedure :: load        => rain_events_table_load
      procedure :: year_window => rain_events_table_year_window
   end type rain_events_table_t

contains

   ! Stubs — real implementations land in Tasks 2-4.
   subroutine meteo_daily_table_load(self, path, errors)
      class(meteo_daily_table_t), intent(inout) :: self
      character(len=*),           intent(in)    :: path
      type(error_collection_t),   intent(inout) :: errors
      ! TASK 2
   end subroutine meteo_daily_table_load

   subroutine meteo_daily_table_year_window(self, year, i1, i2)
      class(meteo_daily_table_t), intent(in)  :: self
      integer,                    intent(in)  :: year
      integer,                    intent(out) :: i1, i2
      i1 = 0; i2 = 0  ! TASK 2
   end subroutine meteo_daily_table_year_window

   subroutine meteo_detail_table_load(self, path, errors)
      class(meteo_detail_table_t), intent(inout) :: self
      character(len=*),            intent(in)    :: path
      type(error_collection_t),    intent(inout) :: errors
      ! TASK 3
   end subroutine meteo_detail_table_load

   subroutine meteo_detail_table_year_window(self, year, i1, i2)
      class(meteo_detail_table_t), intent(in)  :: self
      integer,                     intent(in)  :: year
      integer,                     intent(out) :: i1, i2
      i1 = 0; i2 = 0  ! TASK 3
   end subroutine meteo_detail_table_year_window

   subroutine rain_events_table_load(self, path, errors)
      class(rain_events_table_t), intent(inout) :: self
      character(len=*),           intent(in)    :: path
      type(error_collection_t),   intent(inout) :: errors
      ! TASK 4
   end subroutine rain_events_table_load

   subroutine rain_events_table_year_window(self, year, i1, i2)
      class(rain_events_table_t), intent(in)  :: self
      integer,                    intent(in)  :: year
      integer,                    intent(out) :: i1, i2
      i1 = 0; i2 = 0  ! TASK 4
   end subroutine rain_events_table_year_window

end module meteo_csv_mod
```

- [ ] **Step 7: Run test to verify it passes.**

Run:
```bash
pixi run -e test meson compile -C builddir
pixi run -e test meson test -C builddir test_meteo_csv_suite
```
Expected: PASS for `test_types_exist_and_default` and the suite count incremented (e.g., `OK (1 tests)` for this suite, or larger if integrated into the unit-swap-tests binary).

- [ ] **Step 8: Commit.**

```bash
git add src/io/csv/meteo_csv.f90 tests/unit/io/csv/test_meteo_csv.pf tests/unit/meson.build tests/unit/testSuites.inc meson.build
git commit -m "feat(io): introduce typed meteo CSV table types (skeleton)

First step of the typed-CSV-records arc. Defines meteo_daily_table_t /
meteo_detail_table_t / rain_events_table_t with row record types and
empty load/year_window stubs. Pilot establishes the pattern; six more
CSV families will follow.

Date column stays real64 days-since-1900 (ADR 0014 convention); typed
date_t is a separate arc. is_loaded flag distinguishes 'not requested'
from 'requested-but-empty'.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 2: `meteo_daily_table_t%load` + `%year_window`

**Files:**
- Modify: `src/io/csv/meteo_csv.f90`
- Test: `tests/unit/io/csv/test_meteo_csv.pf`
- Create: `tests/unit/io/fixtures/csv_meteo_daily.csv`
- Create: `tests/unit/io/fixtures/csv_meteo_daily_bad_temp.csv`

- [ ] **Step 1: Create the daily fixture.**

Write `tests/unit/io/fixtures/csv_meteo_daily.csv`:
```
date,rad,tmin,tmax,hum,wind,rain,etref,wet
1980-12-30,5.0,-2.0,4.0,0.8,2.5,1.2,1.0,0.5
1980-12-31,6.0,-1.0,5.0,0.9,2.0,0.0,1.1,0.0
1981-01-01,7.0, 0.5,6.0,0.9,1.8,3.4,1.3,1.0
1981-01-02,8.0, 1.0,7.0,0.85,1.5,0.0,1.2,0.0
1981-12-31,4.5,-3.0,2.0,0.8,3.0,2.0,0.9,1.0
```

- [ ] **Step 2: Create the bad-temp fixture (for validation test).**

Write `tests/unit/io/fixtures/csv_meteo_daily_bad_temp.csv`:
```
date,rad,tmin,tmax,hum,wind,rain,etref,wet
1981-01-01,7.0,10.0,5.0,0.9,1.8,3.4,1.3,1.0
```
(tmin=10, tmax=5 — invalid.)

- [ ] **Step 3: Write failing tests for `load` and `year_window`.**

Append to `tests/unit/io/csv/test_meteo_csv.pf`:
```fortran
@test
subroutine test_daily_load_happy()
   use funit
   use iso_fortran_env, only: real64
   use meteo_csv_mod, only: meteo_daily_table_t
   use error_mod,     only: error_collection_t

   type(meteo_daily_table_t) :: daily
   type(error_collection_t)  :: errors

   call daily%load('tests/unit/io/fixtures/csv_meteo_daily.csv', errors)

   @assertFalse(errors%has_errors())
   @assertTrue(daily%is_loaded)
   @assertTrue(allocated(daily%rows))
   @assertEqual(5, size(daily%rows))

   ! Row 3 = 1981-01-01.  csv_reader returns days-since-1900 = JD - 2415020.
   ! JD(1981-01-01) = 2444606 → 29586.
   @assertEqual(29586.0_real64, daily%rows(3)%date,  1.0_real64)
   @assertEqual(7.0_real64,     daily%rows(3)%rad,   1.0e-6_real64)
   @assertEqual(0.5_real64,     daily%rows(3)%tmin,  1.0e-6_real64)
   @assertEqual(6.0_real64,     daily%rows(3)%tmax,  1.0e-6_real64)
   @assertEqual(3.4_real64,     daily%rows(3)%rain,  1.0e-6_real64)
   @assertEqual(1.0_real64,     daily%rows(3)%wet,   1.0e-6_real64)
end subroutine test_daily_load_happy

@test
subroutine test_daily_validate_tmin_gt_tmax()
   use funit
   use meteo_csv_mod, only: meteo_daily_table_t
   use error_mod,     only: error_collection_t

   type(meteo_daily_table_t) :: daily
   type(error_collection_t)  :: errors

   call daily%load('tests/unit/io/fixtures/csv_meteo_daily_bad_temp.csv', errors)

   @assertTrue(errors%has_errors())
   @assertFalse(daily%is_loaded)
end subroutine test_daily_validate_tmin_gt_tmax

@test
subroutine test_daily_year_window_1981()
   use funit
   use meteo_csv_mod, only: meteo_daily_table_t
   use error_mod,     only: error_collection_t

   type(meteo_daily_table_t) :: daily
   type(error_collection_t)  :: errors
   integer :: i1, i2

   call daily%load('tests/unit/io/fixtures/csv_meteo_daily.csv', errors)
   @assertFalse(errors%has_errors())

   call daily%year_window(1981, i1, i2)
   ! 1981 rows are #3, #4, #5 in the fixture.
   @assertEqual(3, i1)
   @assertEqual(5, i2)
end subroutine test_daily_year_window_1981

@test
subroutine test_daily_year_window_missing_year()
   use funit
   use meteo_csv_mod, only: meteo_daily_table_t
   use error_mod,     only: error_collection_t

   type(meteo_daily_table_t) :: daily
   type(error_collection_t)  :: errors
   integer :: i1, i2

   call daily%load('tests/unit/io/fixtures/csv_meteo_daily.csv', errors)
   call daily%year_window(1999, i1, i2)
   ! No 1999 rows — sentinel (0, 0).
   @assertEqual(0, i1)
   @assertEqual(0, i2)
end subroutine test_daily_year_window_missing_year

@test
subroutine test_daily_load_missing_file()
   use funit
   use meteo_csv_mod, only: meteo_daily_table_t
   use error_mod,     only: error_collection_t

   type(meteo_daily_table_t) :: daily
   type(error_collection_t)  :: errors

   call daily%load('tests/unit/io/fixtures/does_not_exist.csv', errors)
   @assertTrue(errors%has_errors())
   @assertFalse(daily%is_loaded)
end subroutine test_daily_load_missing_file
```

- [ ] **Step 4: Run tests to verify they fail.**

Run:
```bash
pixi run -e test meson compile -C builddir
pixi run -e test meson test -C builddir test_meteo_csv_suite -v 2>&1 | tail -30
```
Expected: 5 new tests, 4 failures (load/year_window stubs) + 1 may pass for the missing-file case (errors already populated by csv_reader if called).

- [ ] **Step 5: Implement `meteo_daily_table_load`.**

Replace the stub body in `src/io/csv/meteo_csv.f90`:
```fortran
   subroutine meteo_daily_table_load(self, path, errors)
      use csv_reader_mod, only: read_csv_table
      use error_mod,      only: ERR_PARSE_TYPE_MISMATCH
      class(meteo_daily_table_t), intent(inout) :: self
      character(len=*),           intent(in)    :: path
      type(error_collection_t),   intent(inout) :: errors

      real(real64), allocatable :: tbl(:,:)
      character(len=5) :: hdr(9)
      integer :: n, r
      character(len=200) :: msg

      hdr = [character(len=5) :: 'date ', 'rad  ', 'tmin ', 'tmax ', &
             'hum  ', 'wind ', 'rain ', 'etref', 'wet  ']
      call read_csv_table(trim(path), hdr, tbl, errors)
      if (errors%has_fatals()) return

      n = size(tbl, 1)
      if (allocated(self%rows)) deallocate(self%rows)
      allocate(self%rows(n))
      do r = 1, n
         self%rows(r)%date  = tbl(r, 1)
         self%rows(r)%rad   = tbl(r, 2)
         self%rows(r)%tmin  = tbl(r, 3)
         self%rows(r)%tmax  = tbl(r, 4)
         self%rows(r)%hum   = tbl(r, 5)
         self%rows(r)%wind  = tbl(r, 6)
         self%rows(r)%rain  = tbl(r, 7)
         self%rows(r)%etref = tbl(r, 8)
         self%rows(r)%wet   = tbl(r, 9)
      end do

      ! Inline validate — load owns this, never call separately.
      do r = 1, n
         if (self%rows(r)%tmin > self%rows(r)%tmax) then
            write(msg, '("meteo daily row ", I0, ": tmin > tmax")') r
            call errors%append(ERR_PARSE_TYPE_MISMATCH, trim(msg), 'meteo_csv')
            deallocate(self%rows)
            return
         end if
         if (self%rows(r)%rain < 0.0_real64) then
            write(msg, '("meteo daily row ", I0, ": rain < 0")') r
            call errors%append(ERR_PARSE_TYPE_MISMATCH, trim(msg), 'meteo_csv')
            deallocate(self%rows)
            return
         end if
      end do

      self%is_loaded = .true.
   end subroutine meteo_daily_table_load
```

- [ ] **Step 6: Implement `meteo_daily_table_year_window`.**

Replace the stub body:
```fortran
   subroutine meteo_daily_table_year_window(self, year, i1, i2)
      class(meteo_daily_table_t), intent(in)  :: self
      integer,                    intent(in)  :: year
      integer,                    intent(out) :: i1, i2

      integer, parameter :: jd1900 = 2415020
      integer  :: i, jday
      external :: jday
      real(real64) :: t_jan1, t_dec31

      i1 = 0; i2 = 0
      if (.not. self%is_loaded) return
      if (.not. allocated(self%rows)) return

      t_jan1  = real(jday(year,  1,  1) - jd1900, real64)
      t_dec31 = real(jday(year, 12, 31) - jd1900, real64)

      do i = 1, size(self%rows)
         if (self%rows(i)%date >= t_jan1 - 0.5_real64 .and. &
             self%rows(i)%date <= t_dec31 + 0.5_real64) then
            if (i1 == 0) i1 = i
            i2 = i
         end if
      end do
   end subroutine meteo_daily_table_year_window
```

- [ ] **Step 7: Add `jday` to the meson.build pfunit_extra_sources.**

`jday` lives in `src/io/readmeteo.f90` at the bottom. Search `tests/unit/meson.build` for `readmeteo` — if not present, add to `pfunit_extra_sources`:
```meson
        '../../src/io/readmeteo.f90',
```
If `readmeteo.f90` brings in other deps that aren't already in the test build, prefer adding `'../../src/core/dtutil.f90'` (already present per Task 0's check) and inlining a local `jday` helper inside meteo_csv. To keep this plan deterministic, **inline a private helper** inside `meteo_csv.f90` instead of dragging in readmeteo:

```fortran
   ! Private helper — same epoch and convention as csv_reader's date col.
   pure function days_since_1900(year, month, day) result(d)
      integer, intent(in) :: year, month, day
      integer :: d
      integer :: a, y, m, jd
      integer, parameter :: jd1900 = 2415020
      a = (14 - month) / 12
      y = year + 4800 - a
      m = month + 12*a - 3
      jd = day + (153*m + 2)/5 + 365*y + y/4 - y/100 + y/400 - 32045
      d = jd - jd1900
   end function days_since_1900
```
Add to the `contains` section of `meteo_csv_mod`. Replace the `jday` call in `meteo_daily_table_year_window` with:
```fortran
      t_jan1  = real(days_since_1900(year,  1,  1), real64)
      t_dec31 = real(days_since_1900(year, 12, 31), real64)
```
(No external dependency. Same convention as `jday` — Fliegel/Van Flandern formula minus the jd1900 offset.)

- [ ] **Step 8: Run tests to verify they pass.**

Run:
```bash
pixi run -e test meson compile -C builddir
pixi run -e test meson test -C builddir test_meteo_csv_suite -v 2>&1 | tail -30
```
Expected: 5/5 tests in the daily group pass (plus the type-defaults test from Task 1 = 6 total).

- [ ] **Step 9: Commit.**

```bash
git add src/io/csv/meteo_csv.f90 tests/unit/io/csv/test_meteo_csv.pf tests/unit/io/fixtures/csv_meteo_daily.csv tests/unit/io/fixtures/csv_meteo_daily_bad_temp.csv
git commit -m "feat(io): implement meteo_daily_table_t load + year_window

load() wraps csv_reader and validates inline (tmin<=tmax, rain>=0); on
validation failure rows are deallocated and is_loaded stays false.
year_window returns (i1,i2) by date column, (0,0) sentinel for missing
years. Private days_since_1900 helper avoids dragging readmeteo into
the test build.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 3: `meteo_detail_table_t%load` + `%year_window`

**Files:**
- Modify: `src/io/csv/meteo_csv.f90`
- Modify: `tests/unit/io/csv/test_meteo_csv.pf`
- Create: `tests/unit/io/fixtures/csv_meteo_detail.csv`

- [ ] **Step 1: Create the detail fixture.**

Write `tests/unit/io/fixtures/csv_meteo_detail.csv`:
```
datetime,record,rad,temp,hum,wind,rain
1981-01-01T00:00:00,1,0.5,2.0,0.9,1.5,0.1
1981-01-01T06:00:00,2,2.5,3.5,0.85,2.0,0.0
1981-01-01T12:00:00,3,8.0,5.5,0.7,3.0,0.0
```

- [ ] **Step 2: Write failing tests.**

Append to `tests/unit/io/csv/test_meteo_csv.pf`:
```fortran
@test
subroutine test_detail_load_happy()
   use funit
   use iso_fortran_env, only: real64
   use meteo_csv_mod, only: meteo_detail_table_t
   use error_mod,     only: error_collection_t

   type(meteo_detail_table_t) :: detail
   type(error_collection_t)   :: errors

   call detail%load('tests/unit/io/fixtures/csv_meteo_detail.csv', errors)

   @assertFalse(errors%has_errors())
   @assertTrue(detail%is_loaded)
   @assertEqual(3, size(detail%rows))
   @assertEqual(2,    detail%rows(2)%record)
   @assertEqual(2.5_real64, detail%rows(2)%rad,  1.0e-6_real64)
   @assertEqual(3.5_real64, detail%rows(2)%temp, 1.0e-6_real64)
end subroutine test_detail_load_happy

@test
subroutine test_detail_year_window_half_open()
   use funit
   use meteo_csv_mod, only: meteo_detail_table_t
   use error_mod,     only: error_collection_t

   type(meteo_detail_table_t) :: detail
   type(error_collection_t)   :: errors
   integer :: i1, i2

   call detail%load('tests/unit/io/fixtures/csv_meteo_detail.csv', errors)
   call detail%year_window(1981, i1, i2)

   @assertEqual(1, i1)
   @assertEqual(3, i2)
end subroutine test_detail_year_window_half_open
```

- [ ] **Step 3: Run tests to verify failures.**

Run:
```bash
pixi run -e test meson compile -C builddir
pixi run -e test meson test -C builddir test_meteo_csv_suite -v 2>&1 | tail -30
```
Expected: 2 new failures.

- [ ] **Step 4: Implement `meteo_detail_table_load` and `meteo_detail_table_year_window`.**

Replace the stub bodies in `src/io/csv/meteo_csv.f90`:
```fortran
   subroutine meteo_detail_table_load(self, path, errors)
      use csv_reader_mod, only: read_csv_table
      class(meteo_detail_table_t), intent(inout) :: self
      character(len=*),            intent(in)    :: path
      type(error_collection_t),    intent(inout) :: errors

      real(real64), allocatable :: tbl(:,:)
      character(len=8) :: hdr(7)
      integer :: n, r

      hdr = [character(len=8) :: 'datetime', 'record  ', 'rad     ', &
             'temp    ', 'hum     ', 'wind    ', 'rain    ']
      call read_csv_table(trim(path), hdr, tbl, errors)
      if (errors%has_fatals()) return

      n = size(tbl, 1)
      if (allocated(self%rows)) deallocate(self%rows)
      allocate(self%rows(n))
      do r = 1, n
         self%rows(r)%datetime = tbl(r, 1)
         self%rows(r)%record   = nint(tbl(r, 2))
         self%rows(r)%rad      = tbl(r, 3)
         self%rows(r)%temp     = tbl(r, 4)
         self%rows(r)%hum      = tbl(r, 5)
         self%rows(r)%wind     = tbl(r, 6)
         self%rows(r)%rain     = tbl(r, 7)
      end do

      self%is_loaded = .true.
   end subroutine meteo_detail_table_load

   subroutine meteo_detail_table_year_window(self, year, i1, i2)
      class(meteo_detail_table_t), intent(in)  :: self
      integer,                     intent(in)  :: year
      integer,                     intent(out) :: i1, i2

      integer :: i
      real(real64) :: t_jan1, t_jan1_next

      i1 = 0; i2 = 0
      if (.not. self%is_loaded) return
      if (.not. allocated(self%rows)) return

      ! Half-open interval [t_jan1, t_jan1_next) — sub-daily timestamps are
      ! fractional days, so the daily ±0.5 slack would mis-bucket Dec 31
      ! noon-to-midnight into the next year (mirrors legacy MeteoCSVDetYear).
      t_jan1      = real(days_since_1900(year,     1, 1), real64)
      t_jan1_next = real(days_since_1900(year + 1, 1, 1), real64)

      do i = 1, size(self%rows)
         if (self%rows(i)%datetime >= t_jan1 .and. &
             self%rows(i)%datetime <  t_jan1_next) then
            if (i1 == 0) i1 = i
            i2 = i
         end if
      end do
   end subroutine meteo_detail_table_year_window
```

- [ ] **Step 5: Run tests to verify they pass.**

Run:
```bash
pixi run -e test meson compile -C builddir
pixi run -e test meson test -C builddir test_meteo_csv_suite -v 2>&1 | tail -30
```
Expected: 8/8 tests pass.

- [ ] **Step 6: Commit.**

```bash
git add src/io/csv/meteo_csv.f90 tests/unit/io/csv/test_meteo_csv.pf tests/unit/io/fixtures/csv_meteo_detail.csv
git commit -m "feat(io): implement meteo_detail_table_t load + year_window

Half-open [jan1, jan1_next) interval matches legacy MeteoCSVDetYear
to prevent FP slack from mis-bucketing Dec 31 noon-to-midnight.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 4: `rain_events_table_t%load` + `%year_window`

**Files:**
- Modify: `src/io/csv/meteo_csv.f90`
- Modify: `tests/unit/io/csv/test_meteo_csv.pf`
- Create: `tests/unit/io/fixtures/csv_meteo_rain_events.csv`

- [ ] **Step 1: Create the rain-events fixture.**

Write `tests/unit/io/fixtures/csv_meteo_rain_events.csv`:
```
datetime,amount
1980-12-30T18:00:00,3.5
1981-01-01T06:30:00,1.2
1981-01-01T18:00:00,0.8
1981-06-15T15:00:00,12.0
```

Note: row 1 is **1980-12-30** (not 12-31) — the legacy `ReadRainEvents`
bucketing uses `>= t_jan1 - 0.5d0` which would otherwise leak a Dec-31-evening
1980 row into the 1981 bucket. The pilot must preserve that behavior for
byte-identical regression; the fixture is shaped to test the clean 1981 case.

- [ ] **Step 2: Write failing tests.**

Append to `tests/unit/io/csv/test_meteo_csv.pf`:
```fortran
@test
subroutine test_rain_load_happy()
   use funit
   use iso_fortran_env, only: real64
   use meteo_csv_mod, only: rain_events_table_t
   use error_mod,     only: error_collection_t

   type(rain_events_table_t) :: rain
   type(error_collection_t)  :: errors

   call rain%load('tests/unit/io/fixtures/csv_meteo_rain_events.csv', errors)

   @assertFalse(errors%has_errors())
   @assertTrue(rain%is_loaded)
   @assertEqual(4, size(rain%rows))
   @assertEqual(1.2_real64,  rain%rows(2)%amount, 1.0e-6_real64)
   @assertEqual(12.0_real64, rain%rows(4)%amount, 1.0e-6_real64)
end subroutine test_rain_load_happy

@test
subroutine test_rain_year_window_1981()
   use funit
   use meteo_csv_mod, only: rain_events_table_t
   use error_mod,     only: error_collection_t

   type(rain_events_table_t) :: rain
   type(error_collection_t)  :: errors
   integer :: i1, i2

   call rain%load('tests/unit/io/fixtures/csv_meteo_rain_events.csv', errors)
   call rain%year_window(1981, i1, i2)
   @assertEqual(2, i1)
   @assertEqual(4, i2)
end subroutine test_rain_year_window_1981
```

- [ ] **Step 3: Run tests to verify failures.**

Run:
```bash
pixi run -e test meson compile -C builddir
pixi run -e test meson test -C builddir test_meteo_csv_suite -v 2>&1 | tail -30
```
Expected: 2 new failures.

- [ ] **Step 4: Implement `rain_events_table_load` and `rain_events_table_year_window`.**

Replace stub bodies in `src/io/csv/meteo_csv.f90`:
```fortran
   subroutine rain_events_table_load(self, path, errors)
      use csv_reader_mod, only: read_csv_table
      class(rain_events_table_t), intent(inout) :: self
      character(len=*),           intent(in)    :: path
      type(error_collection_t),   intent(inout) :: errors

      real(real64), allocatable :: tbl(:,:)
      character(len=8) :: hdr(2)
      integer :: n, r

      hdr = [character(len=8) :: 'datetime', 'amount  ']
      call read_csv_table(trim(path), hdr, tbl, errors)
      if (errors%has_fatals()) return

      n = size(tbl, 1)
      if (allocated(self%rows)) deallocate(self%rows)
      allocate(self%rows(n))
      do r = 1, n
         self%rows(r)%datetime = tbl(r, 1)
         self%rows(r)%amount   = tbl(r, 2)
      end do

      self%is_loaded = .true.
   end subroutine rain_events_table_load

   subroutine rain_events_table_year_window(self, year, i1, i2)
      class(rain_events_table_t), intent(in)  :: self
      integer,                    intent(in)  :: year
      integer,                    intent(out) :: i1, i2

      integer :: i
      real(real64) :: t_jan1, t_dec31

      i1 = 0; i2 = 0
      if (.not. self%is_loaded) return
      if (.not. allocated(self%rows)) return

      ! Matches legacy ReadRainEvents bucketing: ±0.5d slack on year end.
      t_jan1  = real(days_since_1900(year,  1,  1), real64)
      t_dec31 = real(days_since_1900(year, 12, 31), real64) + 1.0_real64

      do i = 1, size(self%rows)
         if (self%rows(i)%datetime >= t_jan1 - 0.5_real64 .and. &
             self%rows(i)%datetime <  t_dec31 + 0.5_real64) then
            if (i1 == 0) i1 = i
            i2 = i
         end if
      end do
   end subroutine rain_events_table_year_window
```

- [ ] **Step 5: Run tests.**

Run:
```bash
pixi run -e test meson compile -C builddir
pixi run -e test meson test -C builddir test_meteo_csv_suite -v 2>&1 | tail -30
```
Expected: 10/10 tests pass.

- [ ] **Step 6: Commit.**

```bash
git add src/io/csv/meteo_csv.f90 tests/unit/io/csv/test_meteo_csv.pf tests/unit/io/fixtures/csv_meteo_rain_events.csv
git commit -m "feat(io): implement rain_events_table_t load + year_window

Bucketing convention matches legacy ReadRainEvents (±0.5d slack on
year end) so the post-migration extractor in Task 8 is byte-identical.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 5: Migrate `state%atmosphere` field layout

**Files:**
- Modify: `src/state/atmosphere_state.f90`

This is the schema change. Per the saved memory, `rm -rf builddir` is mandatory before recompiling.

- [ ] **Step 1: Replace the cache field declarations.**

In `src/state/atmosphere_state.f90`, find the block at lines 266-274 (the `nmetcsv` / `metcsv_dat` / `nmetcsv_det` / `metcsv_det` / `nraincsv` / `raincsv_dat` declarations). Replace with:
```fortran
      ! [METEO-TYPED-CSV 2026-05-28] Typed CSV table records replace the
      ! raw (:,:) caches. Per-year extractors in readmeteo.f90 use
      ! %year_window(year, i1, i2) + %rows(i1:i2)%field accessors.
      type(meteo_daily_table_t)  :: meteo         !! daily meteo CSV cache
      type(meteo_detail_table_t) :: meteo_detail  !! sub-daily meteo CSV cache (swmetdetail=1)
      type(rain_events_table_t)  :: rain_events   !! rain events CSV cache (swrain=3)
```

Add the `use` near the top of the module (under the existing `use` lines around line 5-15):
```fortran
   use meteo_csv_mod, only: meteo_daily_table_t, meteo_detail_table_t, rain_events_table_t
```

- [ ] **Step 2: Replace the three CSV-load blocks in `atmosphere_state_init`.**

In `atmosphere_state_init`, find lines 342-440 (the three `block ... end block` sections that load `metcsv_dat`, `metcsv_det`, `raincsv_dat`). Replace the entire range with:
```fortran
      ! [METEO-TYPED-CSV 2026-05-28] Three typed table loads. Each owns
      ! its schema; validation is inline in load(); is_loaded flag
      ! distinguishes "not requested" from "empty".
      block
         use error_mod, only: error_collection_t
         character(len=300) :: metfile_lc, csvpath
         type(error_collection_t) :: errs

         metfile_lc = ''
         if (allocated(config%meteo%metfile)) metfile_lc = config%meteo%metfile
         call lowerc(metfile_lc)

         if (index(trim(metfile_lc), '.csv') > 0) then
            csvpath = trim(config%general%pathatm) // trim(metfile_lc)
            call self%meteo%load(trim(csvpath), errs)
            call errs%abort_if_fatal()
         end if

         if (config%meteo%swmetdetail == 1) then
            if (allocated(config%meteo%detail_file) .and. &
                len_trim(config%meteo%detail_file) > 0) then
               csvpath = trim(config%general%pathatm) // trim(config%meteo%detail_file)
               call self%meteo_detail%load(trim(csvpath), errs)
               call errs%abort_if_fatal()
            end if
         end if

         if (config%meteo%swrain == 3 .and. allocated(config%meteo%rain_events_file)) then
            if (len_trim(config%meteo%rain_events_file) > 0) then
               csvpath = trim(config%general%pathatm) // trim(config%meteo%rain_events_file)
               call self%rain_events%load(trim(csvpath), errs)
               call errs%abort_if_fatal()
            end if
         end if
      end block
```

Net diff for this block: ~95 lines → ~30.

- [ ] **Step 3: Clean rebuild and run pFUnit.**

Run (the `rm -rf builddir` is mandatory per the saved memory `state-schema-changes-require-clean-rebuild`):
```bash
rm -rf builddir
pixi run -e test meson setup builddir
pixi run -e test meson compile -C builddir 2>&1 | tail -30
```
Expected: compile fails — `readmeteo.f90` still references `state%atmosphere%metcsv_dat` etc. **That's the next task.** Do NOT commit yet. Verify the failure is specifically in `readmeteo.f90` and nowhere else (e.g., grep for stragglers):
```bash
grep -rn 'metcsv_dat\|metcsv_det\|raincsv_dat\|nmetcsv\b\|nmetcsv_det\|nraincsv' src/ tests/unit/ | grep -v 'atmosphere_state.f90'
```
Expected: hits in `src/io/readmeteo.f90` only (and maybe a comment in `meteo_csv.f90` referencing the old names — those are fine).

If hits appear elsewhere, surface them — the migration is broader than this plan scoped for.

---

## Task 6: Migrate `readmeteo.f90` consumers

**Files:**
- Modify: `src/io/readmeteo.f90`

- [ ] **Step 1: Rewrite `MeteoCSVYear` (around lines 436-510).**

Replace the body with:
```fortran
subroutine MeteoCSVYear(ifnd, state)
use error_mod, only: fatalerr_collected
use swap_state_mod, only: swap_state_t
implicit none
integer, intent(out) :: ifnd
type(swap_state_t), intent(inout) :: state

integer  :: i, i1, i2, n
integer  :: datea(6)
real(4)  :: fsec
real(8)  :: t_jan1, tval

associate (atmo => state%atmosphere, &
           time => state%timecontrol)

call atmo%meteo%year_window(time%yearmeteo, i1, i2)

if (i1 == 0) then
   call fatalerr_collected('MeteoCSVYear', &
      'No meteo CSV records found for the requested year')
   ifnd = 0; return
end if

n = i2 - i1 + 1
ifnd = n

! Populate per-day arrays from typed rows.
do i = 1, n
   atmo%arad(i) = atmo%meteo%rows(i1+i-1)%rad * 1000.0d0   ! kJ/m2/d → J/m2/d
   atmo%atmn(i) = atmo%meteo%rows(i1+i-1)%tmin
   atmo%atmx(i) = atmo%meteo%rows(i1+i-1)%tmax
   atmo%ahum(i) = atmo%meteo%rows(i1+i-1)%hum
   atmo%awin(i) = atmo%meteo%rows(i1+i-1)%wind
   atmo%arai(i) = atmo%meteo%rows(i1+i-1)%rain
   atmo%aetr(i) = atmo%meteo%rows(i1+i-1)%etref
   atmo%wet(i)  = atmo%meteo%rows(i1+i-1)%wet
end do

! Backfill ad/am from the date column.
do i = 1, n
   call days1900_to_md(nint(atmo%meteo%rows(i1+i-1)%date), atmo%am(i), atmo%ad(i))
end do

! daynrfirst / daynrlast (unchanged logic).
datea = 0; fsec = 0.0
datea(1) = time%yearmeteo; datea(2) = 1; datea(3) = 1
call dtardp(datea, fsec, t_jan1)
time%timjan1 = t_jan1

datea(2) = atmo%am(1); datea(3) = atmo%ad(1)
call dtardp(datea, fsec, tval)
atmo%daynrfirst = nint(tval - time%timjan1 + 1.0d0)

datea(2) = atmo%am(n); datea(3) = atmo%ad(n)
call dtardp(datea, fsec, tval)
atmo%daynrlast = nint(tval - time%timjan1 + 1.0d0)

end associate
end subroutine MeteoCSVYear
```

- [ ] **Step 2: Rewrite `MeteoCSVDetYear` (around lines 519-594).**

Replace the body with:
```fortran
subroutine MeteoCSVDetYear(ifnd, state)
use error_mod, only: fatalerr_collected
use swap_state_mod, only: swap_state_t
use swap_array_dimensions, only: NMETFILE
implicit none
integer, intent(out) :: ifnd
type(swap_state_t), intent(inout) :: state

integer  :: i, i1, i2, n

associate (atmo => state%atmosphere, &
           time => state%timecontrol)

call atmo%meteo_detail%year_window(time%yearmeteo, i1, i2)

if (i1 == 0) then
   call fatalerr_collected('MeteoCSVDetYear', &
      'No sub-daily meteo CSV records found for the requested year')
   ifnd = 0; return
end if

n = i2 - i1 + 1
if (n > NMETFILE) then
   call fatalerr_collected('MeteoCSVDetYear', &
      'Sub-daily meteo CSV record count exceeds NMETFILE (17568)')
   ifnd = 0; return
end if
ifnd = n

! Allocate state-side detail arrays on first use (NMETFILE-sized for parity).
if (.not. allocated(atmo%dettime)) then
   allocate(atmo%dettime(NMETFILE))
   allocate(atmo%detrecord(NMETFILE))
   allocate(atmo%detrad(NMETFILE))
   allocate(atmo%dettav(NMETFILE))
   allocate(atmo%dethum(NMETFILE))
   allocate(atmo%detwind(NMETFILE))
   allocate(atmo%detrain(NMETFILE))
end if

do i = 1, n
   atmo%dettime(i)   = atmo%meteo_detail%rows(i1+i-1)%datetime
   atmo%detrecord(i) = atmo%meteo_detail%rows(i1+i-1)%record
   atmo%detrad(i)    = atmo%meteo_detail%rows(i1+i-1)%rad * 1000.0d0   ! kJ → J
   atmo%dettav(i)    = atmo%meteo_detail%rows(i1+i-1)%temp
   atmo%dethum(i)    = atmo%meteo_detail%rows(i1+i-1)%hum
   atmo%detwind(i)   = atmo%meteo_detail%rows(i1+i-1)%wind
   atmo%detrain(i)   = atmo%meteo_detail%rows(i1+i-1)%rain
end do

end associate
end subroutine MeteoCSVDetYear
```

- [ ] **Step 3: Rewrite the rain-events scan loop in `ReadRainEvents` (around lines 282-296).**

Find the loop:
```fortran
      do i = 1, state%atmosphere%nraincsv
         if (state%atmosphere%raincsv_dat(i,1) >= t_jan1 - 0.5d0 .and. &
     &       state%atmosphere%raincsv_dat(i,1) <  t_dec31 + 0.5d0) then
            ifnd = ifnd + 1
            state%atmosphere%raintimearray(ifnd) = state%atmosphere%raincsv_dat(i,1)
            state%atmosphere%rainamount(ifnd)    = state%atmosphere%raincsv_dat(i,2)
         end if
      end do
```
Replace with:
```fortran
      block
         integer :: i1, i2, k
         call state%atmosphere%rain_events%year_window(yearmeteo, i1, i2)
         if (i1 > 0) then
            do k = i1, i2
               ifnd = ifnd + 1
               state%atmosphere%raintimearray(ifnd) = state%atmosphere%rain_events%rows(k)%datetime
               state%atmosphere%rainamount(ifnd)    = state%atmosphere%rain_events%rows(k)%amount
            end do
         end if
      end block
```
(The surrounding `t_jan1`/`t_dec31`/`vsmall` etc. stay — they're still used by the zero-prepend and dedup logic later in the subroutine. The `i` variable shadowed by the block is intentional.)

- [ ] **Step 4: Update `read_meteo_from_external_buffer_year` ad/am backfill (around line 408-410).**

The buffer path doesn't read from the cache so it stays largely unchanged. **No edit needed** — verify by reading the function and confirming it uses `get_external_meteo_value(...)`, not `metcsv_dat`. If a reference to `metcsv_dat` is still there, it's a comment only (e.g., the canonical column-order docstring); leave it.

- [ ] **Step 5: Rebuild and run pFUnit + unit-swap-tests.**

Run:
```bash
pixi run -e test meson compile -C builddir 2>&1 | tail -20
pixi run -e test meson test -C builddir 2>&1 | tail -30
```
Expected: full compile, all unit tests pass including `test_meteo_csv_suite` and `test_atmosphere_state_suite`.

- [ ] **Step 6: Verify no stragglers.**

Run:
```bash
grep -rn 'metcsv_dat\|metcsv_det\|raincsv_dat\|nmetcsv\b\|nmetcsv_det\|nraincsv' src/ tests/unit/
```
Expected: only matches inside comments (notably the buffer-column-order docstring) and inside `meteo_csv.f90`'s "replaces" comment. Zero live code references.

- [ ] **Step 7: Commit (single combined commit for state + readmeteo migration).**

```bash
git add src/state/atmosphere_state.f90 src/io/readmeteo.f90
git commit -m "refactor(atmosphere): migrate to typed CSV record tables

state%atmosphere%metcsv_dat / metcsv_det / raincsv_dat (untyped (:,:)
caches) replaced by meteo / meteo_detail / rain_events typed tables.
atmosphere_state_init shrinks 95L -> 30L (three method calls).
MeteoCSVYear / MeteoCSVDetYear / ReadRainEvents reach typed rows by
named field; bare-column-index magic gone. Byte-identical regression.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 7: Regression gate — `check-fast`

**Files:** none — verification only.

- [ ] **Step 1: Run check-fast.**

Run:
```bash
pixi run -e test check-fast 2>&1 | tail -40
```
Expected: ALL PASS. pFUnit count includes the new `test_meteo_csv_suite` tests (look for `OK (N tests)` and verify N rose by 10 compared to baseline). 4 regression cases (hupselbrook, grassgrowth, surfacewater, soilhysteresis) pass byte-identical.

If anything fails:
- pFUnit fail → fix and re-run.
- Regression diff → STOP. The migration is not byte-identical. Examine the diff (look for diff in cumulative quantities or specific dates), bisect against the column mapping in Task 6 Step 1/2/3. The most likely culprits: wrong column → field assignment, off-by-one in `year_window` (e.g., `<` vs `<=`), or `kJ→J` conversion misplaced.

- [ ] **Step 2: Re-snapshot and diff against baseline.**

Run:
```bash
pixi run -e test pytest tests/regression/test_output_regression.py -k hupselbrook -v > .scratch/post_hupselbrook.txt 2>&1
diff .scratch/baseline_hupselbrook.txt .scratch/post_hupselbrook.txt
```
Expected: zero diff. (The harness compares output files byte-by-byte already; this is a defensive double-check.)

---

## Task 8: Regression gate — `check-full`

**Files:** none — verification only.

- [ ] **Step 1: Run check-full.**

Run:
```bash
pixi run -e test check-full 2>&1 | tail -40
```
Expected: ALL PASS. 6 regression cases (check-fast set plus `oxygenstress` and `macropore`/`salinitystress` per the dev docs) pass byte-identical. Target budget ~5 min per the developer docs.

Per the saved memory `feedback-verify-before-committing`, `check-fast` alone misses global-default regressions; `check-full` is the gate for claiming done.

- [ ] **Step 2: If check-full passes, mark this as the gate of record.**

No file edits — this step is a checkpoint. The implementer confirms in their reply that `check-full` passed byte-identical with the new tests registered.

---

## Task 9: Documentation + ADR pointer

**Files:**
- Create: `docs/adr/0044-typed-csv-record-tables.md`
- Modify: `src/io/csv/meteo_csv.f90` (add ADR reference to module header)

- [ ] **Step 1: Write the ADR.**

Create `docs/adr/0044-typed-csv-record-tables.md`:
```markdown
# ADR 0044 — Typed CSV record tables

**Status:** Accepted (2026-05-28)
**Pilot:** `src/io/csv/meteo_csv.f90` — meteo daily/detail/rain-events

## Context

Post-strangler, CSV companion files were cached in `(:,:)` arrays on state
(e.g., `state%atmosphere%metcsv_dat(:,:)`). Schemas lived in three places:
the loader's `expected_header` argument, the consumer's column indices,
and any per-column unit/scaling code. Three points to keep in sync.

## Decision

Each CSV family gets a typed table module under `src/io/csv/<family>.f90`,
defining:

- A row record type (`<family>_row_t`) with one named field per column.
- A table type (`<family>_table_t`) holding `rows(:)`, `is_loaded :: logical`,
  and methods `load(path, errors)` / `year_window(year, i1, i2)`.

Rules:
- `load` wraps the generic `csv_reader_mod` and runs validation inline.
  No separate `validate` method — callers must not get an unvalidated table.
- Unit conversions stay at the *consumer* read site (kJ→J etc.), not in
  `load`. Loader is a pure parser.
- `is_loaded` distinguishes "not requested" from "loaded but empty".
- Date columns stay `real(real64)` days-since-1900 (csv_reader convention).
  A typed `date_t` is a separate arc.

## Consequences

+ Schema lives in one place per family.
+ Consumers read by name (`rows(i)%rad`), not by column index.
+ Validation is inseparable from loading.
+ Each new CSV family is a self-contained module (~250-400 L incl. tests).

- Adding a typed table requires touching 3 places (module, meson.build,
  testSuites.inc). Same as adding any pFUnit suite.

## Pilot scope and follow-ups

Pilot: meteo (daily, detail, rain events). Established 2026-05-28.

Six more CSV families to migrate using this template:
- bottom-boundary tables (gwl, qbot, haquif, hbot, qhbot)
- drainage owltab
- irrigation fixed-events + SSDI
- nutrients amendment-events
- soil initial (h-profile, cml-profile, tsoil)

Each migration closes a chunk of `state%cfg%X` reads in its subsystem,
contributing to the broader `state%cfg` retirement arc.
```

- [ ] **Step 2: Add ADR reference to `meteo_csv.f90` header.**

In `src/io/csv/meteo_csv.f90`, prepend to the docstring:
```fortran
!> Typed CSV record tables for meteorological forcing files. (ADR 0044)
```

- [ ] **Step 3: Update CSV reader README if present.**

Run:
```bash
ls src/io/csv/README.md src/io/README.md 2>/dev/null
```
If either exists, add a one-liner pointing at `meteo_csv_mod` as the typed-table pattern reference. Otherwise create `src/io/csv/README.md`:
```markdown
# CSV layer

`csv_reader.f90` — generic header-validated CSV parser (untyped `(:,:)` output).

`meteo_csv.f90` — typed table records for meteo CSVs (pilot for ADR 0044).
Each family wraps `csv_reader` and exposes typed `rows(:)` with named fields.
```

- [ ] **Step 4: Final regression smoke + commit.**

Run:
```bash
pixi run -e test check-fast 2>&1 | tail -10
```
Expected: PASS.

```bash
git add docs/adr/0044-typed-csv-record-tables.md src/io/csv/meteo_csv.f90 src/io/csv/README.md
git commit -m "docs: ADR 0044 typed CSV record tables; pilot scope notes

Establishes the pattern documented and links the six follow-up
families. Pilot module gains an ADR reference in its header.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 10: Memory update

**Files:** none in the repo — local memory only.

- [ ] **Step 1: After Task 8 passes, write a project memory entry.**

Save `/home/zawadzkim/.claude/projects/-home-zawadzkim-Code-swap/memory/project_meteo_typed_csv_pilot_2026-05-28.md`:
```markdown
---
name: project-meteo-typed-csv-pilot-2026-05-28
description: meteo CSV caches replaced by typed tables (ADR 0044); pilot for 6 more CSV families
metadata:
  type: project
---

Pilot for typed CSV record tables landed 2026-05-28 on branch `meteo-typed-csv`.
src/io/csv/meteo_csv.f90 introduces meteo_daily_table_t / meteo_detail_table_t /
rain_events_table_t with load + year_window methods; load validates inline.
state%atmosphere fields metcsv_dat/_det/raincsv_dat retired in favour of
typed table fields. readmeteo's MeteoCSVYear / MeteoCSVDetYear / ReadRainEvents
migrated to typed rows. Byte-identical regression (check-full PASS).

Six families remain on the typed-CSV roadmap (ADR 0044): bottom-boundary,
drainage owltab, irrigation fixed/SSDI, nutrients events, soil initial
(h-profile, cml-profile, tsoil). Each follows the meteo_csv template;
shrinks its subsystem's state%cfg read footprint as a side-effect.

Connected to [[project-globals-retirement-session-2026-05-20-22]] (the
broader state%cfg arc).
```

- [ ] **Step 2: Add the pointer to MEMORY.md.**

Append to `/home/zawadzkim/.claude/projects/-home-zawadzkim-Code-swap/memory/MEMORY.md`:
```
- [Meteo typed CSV pilot 2026-05-28](project_meteo_typed_csv_pilot_2026-05-28.md) — ADR 0044 pilot: meteo CSV caches typed; 6 families to follow.
```

---

## Closing notes

- **Worktree merge:** `cd /home/zawadzkim/Code/swap && git merge meteo-typed-csv && git worktree remove ../swap-meteo-typed-csv`. Do not push without the user's go-ahead (saved feedback: development-only commits).
- **Out of scope for this pilot (deliberately):** `meteo_io.f90`'s `ReadMeteoDay` (it reads already-extracted per-day arrays, not the cache); `meteo_orchestrator.f90`; the broader `state%cfg` retirement; the date-as-typed-record question.
- **Pattern reference for follow-ups:** the six follow-up families should use this exact shape — one module per family, row + table types, `load` + `year_window`, validation inline. Don't drift the API.
