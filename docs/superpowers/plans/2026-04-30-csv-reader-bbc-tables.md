# Unified CSV Reader + BBC Tabular Companion Files — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the ad-hoc `read_csv_date_reals` reader with a single header-validated `read_csv_table` reader, and migrate every tabular bottom-boundary input from inline TOML 2D arrays to CSV companion files.

**Architecture:** New reader in `csv_reader_mod` with strict positional header validation; the first column is parsed as ISO-date or real depending on whether `expected_header(1) == 'date'`. Schema gains six `*_file` slots in `[bottom_boundary]` (one per SWBOTB sub-mode that drives a table). Sub-mode integer switches (`sw2`, `sw3`, `sw4`, `swqhbot`) are added to the schema so per-sub-mode validator rules can fire deterministically. The adapter opens each CSV at adapter time after validation; load-validate-adapt separation is preserved.

**Tech stack:** Fortran 2008 (gfortran), pFUnit unit tests, tomlf for TOML parsing, pixi-managed build/test runner.

**Spec:** `docs/superpowers/specs/2026-04-30-csv-reader-bbc-tables-design.md`

---

## File Structure

| File | Role |
|---|---|
| `src/utils/error_mod.f90` | Add 4 new error codes. |
| `src/io/csv_reader.f90` | Replace `read_csv_date_reals` with `read_csv_table`. |
| `src/config/bottom_boundary_config.f90` | Add 6 `*_file` slots + 4 sub-mode switches; remove `gwl_table`, `haquif_table`, `bbcfil`; rewrite validator. |
| `src/io/toml/read_bottom_boundary_toml.f90` | Read 6 `*_file` keys + 4 sub-mode switches; drop `gwl_table` / `haquif_table` / `bbcfil` reads. |
| `src/io/toml/config_to_variables.f90` | Adapter: per-sub-mode `read_csv_table` calls; migrate irrigation `fixed_events_file` to new reader. |
| `tests/unit/io/test_csv_reader.pf` | Rewritten — 13 tests covering happy + every error path. |
| `tests/unit/io/fixtures/csv_*.csv` | New fixture set. |
| `tests/unit/config/test_bottom_boundary_config.pf` | Update validator tests; remove tests for removed slots. |
| `tests/unit/io/toml/test_read_bottom_boundary_toml.pf` | Update reader tests for new slots. |
| `tests/swap-cases/toml/2.grassgrowth/swap.toml` | Replace inline `gwl_table` with `gwl_file`. |
| `tests/swap-cases/toml/2.grassgrowth/grassgrowth.gwl.csv` | New CSV (date,gwl) — 125 rows. |
| `tests/swap-cases/toml/4.oxygenstress/swap.toml` | Replace inline `haquif_table` with `haquif_file` + `sw3=2`. |
| `tests/swap-cases/toml/4.oxygenstress/oxygenstress.haquif.csv` | New CSV (date,haquif) — current inline rows. |
| `tests/swap-cases/toml/6.surfacewater/swap.toml` | Replace inline `haquif_table` with `haquif_file` + `sw3=2`. |
| `tests/swap-cases/toml/6.surfacewater/surfacewater.haquif.csv` | New CSV (date,haquif) — current inline rows. |
| `docs/csv-companion-files.md` | Update to reflect unified reader + BBC slot list. |

---

## Task 1: Add new error codes

**Files:**
- Modify: `src/error/error.f90` (add four `integer, parameter, public` lines)

- [ ] **Step 1: Add error code constants**

Edit `src/error/error.f90` — insert after the existing `ERR_PARSE_*` block (after `ERR_PARSE_MISSING_REQUIRED = 202`):

```fortran
   integer, parameter, public :: ERR_PARSE_MISSING_HEADER     = 203
   integer, parameter, public :: ERR_PARSE_HEADER_MISMATCH    = 204
   integer, parameter, public :: ERR_PARSE_ROW_SHAPE          = 205
   integer, parameter, public :: ERR_IO_OPEN_FAILED           = 102
   integer, parameter, public :: ERR_VALIDATION_REQUIRED        = 304
```

(`ERR_VALIDATION_REQUIRED = 304` slots into the validation 30x range; `ERR_IO_OPEN_FAILED = 102` slots into the I/O 10x range.)

- [ ] **Step 2: Build to confirm error_mod compiles**

Run: `pixi run -e test build`
Expected: clean build.

- [ ] **Step 3: Commit**

```bash
git add src/error/error.f90
git commit -m "feat(error): add codes for CSV header / row-shape / required-file errors

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: New CSV reader fixtures

**Files:**
- Create: `tests/unit/io/fixtures/csv_date_gwl.csv`
- Create: `tests/unit/io/fixtures/csv_htab_qtab.csv`
- Create: `tests/unit/io/fixtures/csv_irrigation.csv`
- Create: `tests/unit/io/fixtures/csv_gwl_date_swapped.csv`
- Create: `tests/unit/io/fixtures/csv_uppercase_header.csv`
- Create: `tests/unit/io/fixtures/csv_extra_column.csv`
- Create: `tests/unit/io/fixtures/csv_no_header.csv`
- Create: `tests/unit/io/fixtures/csv_bad_real.csv`
- Create: `tests/unit/io/fixtures/csv_bad_date.csv`
- Create: `tests/unit/io/fixtures/csv_short_row.csv`
- Create: `tests/unit/io/fixtures/csv_header_only.csv`
- Create: `tests/unit/io/fixtures/csv_with_comments.csv`
- Delete: `tests/unit/io/fixtures/csv_happy.csv`, `tests/unit/io/fixtures/csv_malformed.csv`

- [ ] **Step 1: Write each fixture**

Use the `Write` tool for each file with the exact contents below. (No trailing blank line on any file.)

`csv_date_gwl.csv`:
```
date,gwl
1980-04-24,-88.0
1980-05-02,-121.0
1980-05-07,-98.0
```

`csv_htab_qtab.csv`:
```
htab,qtab
-100.0,0.0
-200.0,-0.5
-300.0,-1.0
```

`csv_irrigation.csv`:
```
# Sample CSV — three rows.
date,depth,conc,type
2012-04-17,14.40,0.401,1
2012-04-18,14.40,0.401,1
2012-04-19,14.40,0.401,1
```

`csv_gwl_date_swapped.csv` (header columns flipped):
```
gwl,date
-88.0,1980-04-24
```

`csv_uppercase_header.csv`:
```
DATE,GWL
1980-04-24,-88.0
```

`csv_extra_column.csv`:
```
date,gwl,note
1980-04-24,-88.0,seed
```

`csv_no_header.csv` (data rows only — no header line at all):
```
1980-04-24,-88.0
1980-04-25,-90.0
```

`csv_bad_real.csv`:
```
date,gwl
1980-04-24,-88.0
1980-04-25,not_a_number
1980-04-26,-90.0
```

`csv_bad_date.csv`:
```
date,gwl
1980-04-24,-88.0
not-a-date,-90.0
```

`csv_short_row.csv`:
```
date,gwl
1980-04-24,-88.0
1980-04-25
```

`csv_header_only.csv`:
```
date,gwl
```

`csv_with_comments.csv`:
```
# leading comment
# another comment
date,gwl

1980-04-24,-88.0

# inline comment between rows
1980-04-25,-90.0
```

- [ ] **Step 2: Delete obsolete fixtures**

```bash
git rm tests/unit/io/fixtures/csv_happy.csv tests/unit/io/fixtures/csv_malformed.csv
```

- [ ] **Step 3: Commit**

```bash
git add tests/unit/io/fixtures/
git commit -m "test(csv-reader): add fixture set for header-validated reader

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: Failing reader tests (TDD red)

**Files:**
- Modify (full rewrite): `tests/unit/io/test_csv_reader.pf`

- [ ] **Step 1: Replace the test file**

Overwrite `tests/unit/io/test_csv_reader.pf` with:

```fortran
! Tests for csv_reader_mod (Phase 4f cleanup — unified reader).

@test
subroutine test_date_keyed_happy()
   use funit
   use iso_fortran_env, only: real64
   use csv_reader_mod, only: read_csv_table
   use error_mod, only: error_collection_t
   real(real64), allocatable :: table(:,:)
   type(error_collection_t)  :: errors

   call read_csv_table('tests/unit/io/fixtures/csv_date_gwl.csv', &
                       ['date', 'gwl '], table, errors)

   @assertFalse(errors%has_errors())
   @assertTrue(allocated(table))
   @assertEqual(3, size(table, 1))
   @assertEqual(2, size(table, 2))
   ! 1980-04-24 = days-since-1900 = 29334.
   @assertEqual(29334.0_real64, table(1, 1), 1.0_real64)
   @assertEqual(-88.0_real64,  table(1, 2), 1.0e-6_real64)
end subroutine

@test
subroutine test_real_keyed_happy()
   use funit
   use iso_fortran_env, only: real64
   use csv_reader_mod, only: read_csv_table
   use error_mod, only: error_collection_t
   real(real64), allocatable :: table(:,:)
   type(error_collection_t)  :: errors

   call read_csv_table('tests/unit/io/fixtures/csv_htab_qtab.csv', &
                       ['htab', 'qtab'], table, errors)

   @assertFalse(errors%has_errors())
   @assertEqual(3, size(table, 1))
   @assertEqual(2, size(table, 2))
   @assertEqual(-100.0_real64, table(1, 1), 1.0e-6_real64)
   @assertEqual(0.0_real64,    table(1, 2), 1.0e-6_real64)
end subroutine

@test
subroutine test_wide_table_happy()
   use funit
   use iso_fortran_env, only: real64
   use csv_reader_mod, only: read_csv_table
   use error_mod, only: error_collection_t
   real(real64), allocatable :: table(:,:)
   type(error_collection_t)  :: errors

   call read_csv_table('tests/unit/io/fixtures/csv_irrigation.csv', &
                       ['date ', 'depth', 'conc ', 'type '], table, errors)

   @assertFalse(errors%has_errors())
   @assertEqual(3, size(table, 1))
   @assertEqual(4, size(table, 2))
   @assertEqual(14.40_real64, table(1, 2), 1.0e-6_real64)
   @assertEqual(0.401_real64, table(1, 3), 1.0e-6_real64)
   @assertEqual(1.0_real64,   table(1, 4), 1.0e-6_real64)
end subroutine

@test
subroutine test_header_wrong_order()
   use funit
   use iso_fortran_env, only: real64
   use csv_reader_mod, only: read_csv_table
   use error_mod, only: error_collection_t, ERR_PARSE_HEADER_MISMATCH
   real(real64), allocatable :: table(:,:)
   type(error_collection_t)  :: errors

   call read_csv_table('tests/unit/io/fixtures/csv_gwl_date_swapped.csv', &
                       ['date', 'gwl '], table, errors)

   @assertTrue(errors%has_errors())
   @assertEqual(ERR_PARSE_HEADER_MISMATCH, errors%items(1)%code)
   @assertFalse(allocated(table))
end subroutine

@test
subroutine test_header_case_mismatch()
   use funit
   use iso_fortran_env, only: real64
   use csv_reader_mod, only: read_csv_table
   use error_mod, only: error_collection_t, ERR_PARSE_HEADER_MISMATCH
   real(real64), allocatable :: table(:,:)
   type(error_collection_t)  :: errors

   call read_csv_table('tests/unit/io/fixtures/csv_uppercase_header.csv', &
                       ['date', 'gwl '], table, errors)

   @assertTrue(errors%has_errors())
   @assertEqual(ERR_PARSE_HEADER_MISMATCH, errors%items(1)%code)
end subroutine

@test
subroutine test_header_extra_column()
   use funit
   use iso_fortran_env, only: real64
   use csv_reader_mod, only: read_csv_table
   use error_mod, only: error_collection_t, ERR_PARSE_HEADER_MISMATCH
   real(real64), allocatable :: table(:,:)
   type(error_collection_t)  :: errors

   call read_csv_table('tests/unit/io/fixtures/csv_extra_column.csv', &
                       ['date', 'gwl '], table, errors)

   @assertTrue(errors%has_errors())
   @assertEqual(ERR_PARSE_HEADER_MISMATCH, errors%items(1)%code)
end subroutine

@test
subroutine test_header_missing()
   use funit
   use iso_fortran_env, only: real64
   use csv_reader_mod, only: read_csv_table
   use error_mod, only: error_collection_t, ERR_PARSE_HEADER_MISMATCH
   real(real64), allocatable :: table(:,:)
   type(error_collection_t)  :: errors

   call read_csv_table('tests/unit/io/fixtures/csv_no_header.csv', &
                       ['date', 'gwl '], table, errors)

   @assertTrue(errors%has_errors())
   ! First line "1980-04-24,-88.0" is treated as the header attempt
   ! and fails strict-positional comparison.
   @assertEqual(ERR_PARSE_HEADER_MISMATCH, errors%items(1)%code)
end subroutine

@test
subroutine test_missing_file()
   use funit
   use iso_fortran_env, only: real64
   use csv_reader_mod, only: read_csv_table
   use error_mod, only: error_collection_t, ERR_IO_OPEN_FAILED
   real(real64), allocatable :: table(:,:)
   type(error_collection_t)  :: errors

   call read_csv_table('tests/unit/io/fixtures/does_not_exist.csv', &
                       ['date', 'gwl '], table, errors)

   @assertTrue(errors%has_errors())
   @assertEqual(ERR_IO_OPEN_FAILED, errors%items(1)%code)
   @assertFalse(allocated(table))
end subroutine

@test
subroutine test_malformed_real()
   use funit
   use iso_fortran_env, only: real64
   use csv_reader_mod, only: read_csv_table
   use error_mod, only: error_collection_t, ERR_PARSE_TYPE_MISMATCH
   real(real64), allocatable :: table(:,:)
   type(error_collection_t)  :: errors

   call read_csv_table('tests/unit/io/fixtures/csv_bad_real.csv', &
                       ['date', 'gwl '], table, errors)

   @assertTrue(errors%has_errors())
   @assertEqual(ERR_PARSE_TYPE_MISMATCH, errors%items(1)%code)
   @assertFalse(allocated(table))
end subroutine

@test
subroutine test_malformed_date()
   use funit
   use iso_fortran_env, only: real64
   use csv_reader_mod, only: read_csv_table
   use error_mod, only: error_collection_t, ERR_PARSE_TYPE_MISMATCH
   real(real64), allocatable :: table(:,:)
   type(error_collection_t)  :: errors

   call read_csv_table('tests/unit/io/fixtures/csv_bad_date.csv', &
                       ['date', 'gwl '], table, errors)

   @assertTrue(errors%has_errors())
   @assertEqual(ERR_PARSE_TYPE_MISMATCH, errors%items(1)%code)
end subroutine

@test
subroutine test_row_shape_mismatch()
   use funit
   use iso_fortran_env, only: real64
   use csv_reader_mod, only: read_csv_table
   use error_mod, only: error_collection_t, ERR_PARSE_ROW_SHAPE
   real(real64), allocatable :: table(:,:)
   type(error_collection_t)  :: errors

   call read_csv_table('tests/unit/io/fixtures/csv_short_row.csv', &
                       ['date', 'gwl '], table, errors)

   @assertTrue(errors%has_errors())
   @assertEqual(ERR_PARSE_ROW_SHAPE, errors%items(1)%code)
end subroutine

@test
subroutine test_empty_table()
   use funit
   use iso_fortran_env, only: real64
   use csv_reader_mod, only: read_csv_table
   use error_mod, only: error_collection_t
   real(real64), allocatable :: table(:,:)
   type(error_collection_t)  :: errors

   call read_csv_table('tests/unit/io/fixtures/csv_header_only.csv', &
                       ['date', 'gwl '], table, errors)

   @assertFalse(errors%has_errors())
   @assertTrue(allocated(table))
   @assertEqual(0, size(table, 1))
   @assertEqual(2, size(table, 2))
end subroutine

@test
subroutine test_comments_and_blanks()
   use funit
   use iso_fortran_env, only: real64
   use csv_reader_mod, only: read_csv_table
   use error_mod, only: error_collection_t
   real(real64), allocatable :: table(:,:)
   type(error_collection_t)  :: errors

   call read_csv_table('tests/unit/io/fixtures/csv_with_comments.csv', &
                       ['date', 'gwl '], table, errors)

   @assertFalse(errors%has_errors())
   @assertEqual(2, size(table, 1))
   @assertEqual(2, size(table, 2))
end subroutine
```

Note on `expected_header` literals: pFUnit/Fortran array literals require uniform string lengths, so we right-pad with trailing spaces (e.g. `'gwl '` to match `'date'`). The reader compares with `trim()` on both sides, so the padding is invisible.

- [ ] **Step 2: Build and run the test (must fail)**

Run: `pixi run -e test test-pfunit -- --filter test_csv_reader`
Expected: compile failure ("read_csv_table is not defined") OR all tests fail. Either is the red state.

---

## Task 4: Implement `read_csv_table` (TDD green)

**Files:**
- Modify (full rewrite): `src/io/csv_reader.f90`

- [ ] **Step 1: Replace the reader**

Overwrite `src/io/csv_reader.f90`:

```fortran
!> Unified CSV reader for SWAP tabular companion files.
!!
!! One public sub: `read_csv_table(path, expected_header, table, errors)`.
!! Validates the file header strict-positional (lowercase, trimmed) against
!! the caller-supplied expected_header. Column 1 is parsed as ISO-date when
!! `expected_header(1) == 'date'`, else as real64. Subsequent columns are
!! always real64. Comment (`#`) and blank lines are skipped. Errors flow
!! through error_collection_t and leave `table` unallocated on failure.
module csv_reader_mod
   use iso_fortran_env, only: real64, iostat_end
   use error_mod, only: error_collection_t,         &
                        ERR_IO_OPEN_FAILED,         &
                        ERR_IO_READ_FAILED,         &
                        ERR_PARSE_TYPE_MISMATCH,    &
                        ERR_PARSE_MISSING_HEADER,   &
                        ERR_PARSE_HEADER_MISMATCH,  &
                        ERR_PARSE_ROW_SHAPE
   implicit none
   private

   public :: read_csv_table

contains

   subroutine read_csv_table(path, expected_header, table, errors)
      character(len=*),          intent(in)    :: path
      character(len=*),          intent(in)    :: expected_header(:)
      real(real64), allocatable, intent(out)   :: table(:,:)
      type(error_collection_t),  intent(inout) :: errors

      integer :: unit, ios, ncols, nrows, irow
      character(len=4096) :: line
      logical :: file_exists, header_done, date_keyed
      character(len=64) :: irow_str

      ncols = size(expected_header)
      date_keyed = (ncols >= 1) .and. (to_lower(trim(expected_header(1))) == 'date')

      inquire(file=path, exist=file_exists)
      if (.not. file_exists) then
         call errors%append(ERR_IO_OPEN_FAILED, &
            "cannot open CSV: " // trim(path), 'csv_reader')
         return
      end if

      open(newunit=unit, file=path, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         call errors%append(ERR_IO_OPEN_FAILED, &
            "cannot open CSV: " // trim(path), 'csv_reader')
         return
      end if

      ! Pass 1: locate header, count data rows.
      header_done = .false.
      nrows = 0
      do
         read(unit, '(A)', iostat=ios) line
         if (ios == iostat_end) exit
         if (ios /= 0) then
            call errors%append(ERR_IO_READ_FAILED, &
               "read error: " // trim(path), 'csv_reader')
            close(unit)
            return
         end if
         if (is_skippable(line)) cycle
         if (.not. header_done) then
            call validate_header(line, expected_header, path, errors)
            if (errors%has_fatals()) then
               close(unit)
               return
            end if
            header_done = .true.
            cycle
         end if
         nrows = nrows + 1
      end do

      if (.not. header_done) then
         call errors%append(ERR_PARSE_MISSING_HEADER, &
            trim(path) // ": missing header line", 'csv_reader')
         close(unit)
         return
      end if

      allocate(table(nrows, ncols))
      table = 0.0_real64
      if (nrows == 0) then
         close(unit)
         return
      end if

      ! Pass 2: parse data rows.
      rewind(unit)
      header_done = .false.
      irow = 0
      do
         read(unit, '(A)', iostat=ios) line
         if (ios == iostat_end) exit
         if (ios /= 0) then
            call errors%append(ERR_IO_READ_FAILED, &
               "read error: " // trim(path), 'csv_reader')
            close(unit)
            if (allocated(table)) deallocate(table)
            return
         end if
         if (is_skippable(line)) cycle
         if (.not. header_done) then
            header_done = .true.
            cycle
         end if
         irow = irow + 1
         call parse_row(line, irow, ncols, date_keyed, path, table, errors)
         if (errors%has_fatals()) then
            close(unit)
            if (allocated(table)) deallocate(table)
            return
         end if
      end do

      close(unit)
   end subroutine read_csv_table

   pure function is_skippable(line) result(skip)
      character(len=*), intent(in) :: line
      logical :: skip
      character(len=:), allocatable :: trimmed
      skip = .false.
      trimmed = adjustl(line)
      if (len_trim(trimmed) == 0) then
         skip = .true.
         return
      end if
      if (trimmed(1:1) == '#') skip = .true.
   end function is_skippable

   subroutine validate_header(line, expected, path, errors)
      character(len=*),         intent(in)    :: line
      character(len=*),         intent(in)    :: expected(:)
      character(len=*),         intent(in)    :: path
      type(error_collection_t), intent(inout) :: errors
      character(len=:), allocatable :: rest, field
      integer :: i, comma, n
      n = size(expected)
      rest = adjustl(line)
      do i = 1, n
         if (i < n) then
            comma = index(rest, ',')
            if (comma == 0) then
               call errors%append(ERR_PARSE_HEADER_MISMATCH, &
                  trim(path) // ": header column count mismatch", 'csv_reader')
               return
            end if
            field = trim(adjustl(rest(1:comma-1)))
            rest  = adjustl(rest(comma+1:))
         else
            comma = index(rest, ',')
            if (comma /= 0) then
               ! Extra column after the last expected one.
               call errors%append(ERR_PARSE_HEADER_MISMATCH, &
                  trim(path) // ": header has extra columns", 'csv_reader')
               return
            end if
            field = trim(adjustl(rest))
         end if
         if (to_lower(field) /= to_lower(trim(expected(i)))) then
            call errors%append(ERR_PARSE_HEADER_MISMATCH, &
               trim(path) // ": header column '" // field // &
               "' does not match expected '" // trim(expected(i)) // "'", &
               'csv_reader')
            return
         end if
      end do
   end subroutine validate_header

   subroutine parse_row(line, irow, ncols, date_keyed, path, table, errors)
      character(len=*),         intent(in)    :: line
      integer,                  intent(in)    :: irow, ncols
      logical,                  intent(in)    :: date_keyed
      character(len=*),         intent(in)    :: path
      real(real64),             intent(inout) :: table(:,:)
      type(error_collection_t), intent(inout) :: errors
      character(len=:), allocatable :: rest, field
      character(len=64) :: irow_str
      integer :: j, comma, ios
      real(real64) :: val, days
      rest = adjustl(line)
      do j = 1, ncols
         if (j < ncols) then
            comma = index(rest, ',')
            if (comma == 0) then
               write(irow_str, '(I0)') irow
               call errors%append(ERR_PARSE_ROW_SHAPE, &
                  trim(path) // ":row " // trim(irow_str) // &
                  ": expected " // int_str(ncols) // " fields", 'csv_reader')
               return
            end if
            field = trim(adjustl(rest(1:comma-1)))
            rest  = adjustl(rest(comma+1:))
         else
            comma = index(rest, ',')
            if (comma /= 0) then
               write(irow_str, '(I0)') irow
               call errors%append(ERR_PARSE_ROW_SHAPE, &
                  trim(path) // ":row " // trim(irow_str) // &
                  ": expected " // int_str(ncols) // " fields", 'csv_reader')
               return
            end if
            field = trim(adjustl(rest))
         end if

         if (j == 1 .and. date_keyed) then
            if (.not. parse_iso_date(field, days)) then
               write(irow_str, '(I0)') irow
               call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                  trim(path) // ":row " // trim(irow_str) // " col 1: '" // &
                  field // "' not an ISO date", 'csv_reader')
               return
            end if
            table(irow, 1) = days
         else
            read(field, *, iostat=ios) val
            if (ios /= 0) then
               write(irow_str, '(I0)') irow
               call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                  trim(path) // ":row " // trim(irow_str) // " col " // &
                  int_str(j) // ": '" // field // "' not a real number", &
                  'csv_reader')
               return
            end if
            table(irow, j) = val
         end if
      end do
   end subroutine parse_row

   pure function to_lower(s) result(out)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: out
      integer :: i, c
      do i = 1, len(s)
         c = iachar(s(i:i))
         if (c >= iachar('A') .and. c <= iachar('Z')) c = c + 32
         out(i:i) = achar(c)
      end do
   end function to_lower

   pure function int_str(i) result(s)
      integer, intent(in) :: i
      character(len=:), allocatable :: s
      character(len=32) :: buf
      write(buf, '(I0)') i
      s = trim(buf)
   end function int_str

   function parse_iso_date(s, days) result(ok)
      character(len=*), intent(in)  :: s
      real(real64),     intent(out) :: days
      logical :: ok
      integer :: y, m, d, ios, jd, jd1900
      days = 0.0_real64
      ok = .false.
      if (len_trim(s) < 10) return
      if (s(5:5) /= '-' .or. s(8:8) /= '-') return
      read(s(1:4),  '(I4)', iostat=ios) y; if (ios /= 0) return
      read(s(6:7),  '(I2)', iostat=ios) m; if (ios /= 0) return
      read(s(9:10), '(I2)', iostat=ios) d; if (ios /= 0) return
      if (m < 1 .or. m > 12) return
      if (d < 1 .or. d > 31) return
      jd     = julian_day(y, m, d)
      jd1900 = 2415020
      days   = real(jd - jd1900, kind=real64)
      ok = .true.
   end function parse_iso_date

   pure function julian_day(y, m, d) result(jd)
      integer, intent(in) :: y, m, d
      integer :: jd, a, yy, mm
      a  = (14 - m) / 12
      yy = y + 4800 - a
      mm = m + 12 * a - 3
      jd = d + (153 * mm + 2) / 5 + 365 * yy + yy / 4 - yy / 100 + yy / 400 - 32045
   end function julian_day

end module csv_reader_mod
```

- [ ] **Step 2: Run the new tests (must pass)**

Run: `pixi run -e test test-pfunit -- --filter test_csv_reader`
Expected: 13 PASS.

- [ ] **Step 3: Commit**

```bash
git add src/io/csv_reader.f90 tests/unit/io/test_csv_reader.pf
git commit -m "feat(csv-reader): unified header-validated read_csv_table

Replaces read_csv_date_reals (heuristic header skip) with strict
positional header validation. Column 1 is date or real depending on
expected_header(1).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 5: Migrate irrigation `fixed_events_file` to new reader

**Files:**
- Modify: `src/io/toml/config_to_variables.f90` (around lines 728-749 — the `fixed_events_file` block)

The current adapter calls `read_csv_date_reals(.., 3, csv_table, ..)`. Migrate to `read_csv_table` with `expected_header=['date','depth','conc','type']`.

- [ ] **Step 1: Replace the call site**

Find the block starting `if (len_trim(config%irrigation%fixed_events_file) > 0) then` and replace its body. The new body:

```fortran
         if (len_trim(config%irrigation%fixed_events_file) > 0) then
            block
               use csv_reader_mod, only: read_csv_table
               use error_mod, only: error_collection_t
               real(8), allocatable :: csv_table(:,:)
               type(error_collection_t) :: csv_errs
               integer :: k_csv, nrows_csv
               character(len=5) :: irrig_header(4)
               irrig_header(1) = 'date '
               irrig_header(2) = 'depth'
               irrig_header(3) = 'conc '
               irrig_header(4) = 'type '
               call read_csv_table( &
                  trim(pathwork)//trim(config%irrigation%fixed_events_file), &
                  irrig_header, csv_table, csv_errs)
               call csv_errs%abort_if_fatal()
               if (allocated(csv_table)) then
                  nrows_csv = size(csv_table, 1)
                  do k_csv = 1, min(nrows_csv, size(irdate))
                     irdate(k_csv)  = csv_table(k_csv, 1)
                     irdepth(k_csv) = csv_table(k_csv, 2) / 10.0d0  ! mm -> cm
                     irconc(k_csv)  = csv_table(k_csv, 3)
                     irtype(k_csv)  = nint(csv_table(k_csv, 4))
                  end do
               end if
            end block
         end if
```

- [ ] **Step 2: Run the salinitystress regression**

Run: `pixi run -e test regression -- -k salinitystress`
Expected: PASS — header `date,depth,conc,type` matches `salinitystress.irg.csv`, output unchanged.

- [ ] **Step 3: Commit**

```bash
git add src/io/toml/config_to_variables.f90
git commit -m "refactor(adapter): migrate fixed_events_file to read_csv_table

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 6: BBC schema — add `*_file` slots and sub-mode switches; remove obsolete slots

**Files:**
- Modify: `src/config/bottom_boundary_config.f90`

- [ ] **Step 1: Add new fields, remove old**

Within `type :: bottom_boundary_config_t`:

Add (insert below `swbotb` declaration):
```fortran
      ! Sub-mode switches (legacy SW2/SW3/SW4/SWQHBOT) — see core/variables.f90.
      ! sw2: 1=sine, 2=table; sw3: 1=sine, 2=table; sw4: 0=no extra flux,
      ! 1=include extra flux; swqhbot: 1=exponential, 2=tabular.
      ! Validators key on (swbotb, sw_x) combinations.
      integer :: sw2     = 1
      integer :: sw3     = 1
      integer :: sw4     = 0
      integer :: swqhbot = 1

      ! Phase 4f cleanup: per-sub-mode CSV companion file paths.
      character(len=:), allocatable :: gwl_file     ! SWBOTB=1
      character(len=:), allocatable :: qbot2_file   ! SWBOTB=2 + sw2=2
      character(len=:), allocatable :: haquif_file  ! SWBOTB=3 + sw3=2
      character(len=:), allocatable :: qbot4_file   ! SWBOTB=3 + sw4=1
      character(len=:), allocatable :: qhbot_file   ! SWBOTB=4 + swqhbot=2
      character(len=:), allocatable :: hbot5_file   ! SWBOTB=5
```

Remove:
```fortran
      character(len=:), allocatable :: bbcfil
      real(real64), allocatable :: gwl_table(:,:)
      real(real64), allocatable :: haquif_table(:,:)
```

(Leave `swc_table`, `qbot_table`, `cofqha_table`, scalar `hbot`/`rhobot`/`shape`/etc untouched — they are not in scope of this spec.)

- [ ] **Step 2: Rewrite the validator**

Replace the body of `bottom_boundary_config_validate` with:

```fortran
   subroutine bottom_boundary_config_validate(self, errors)
      use error_mod, only: ERR_VALIDATION_REQUIRED
      class(bottom_boundary_config_t), intent(in)    :: self
      type(error_collection_t),         intent(inout) :: errors
      integer :: i

      if (self%swbotb == 0) return

      call check_int_enum(self%swbotb, [(i, i=1, 8)], 'bottom_boundary.swbotb', errors)
      call check_int_enum(self%sw2,     [1, 2], 'bottom_boundary.sw2',     errors)
      call check_int_enum(self%sw3,     [1, 2], 'bottom_boundary.sw3',     errors)
      call check_int_enum(self%sw4,     [0, 1], 'bottom_boundary.sw4',     errors)
      call check_int_enum(self%swqhbot, [1, 2], 'bottom_boundary.swqhbot', errors)

      select case (self%swbotb)
      case (1)
         if (.not. has_file(self%gwl_file)) then
            call errors%append(ERR_VALIDATION_REQUIRED, &
               "bottom_boundary.gwl_file required when swbotb=1", &
               'bottom_boundary')
         end if
      case (2)
         if (self%sw2 == 2 .and. .not. has_file(self%qbot2_file)) then
            call errors%append(ERR_VALIDATION_REQUIRED, &
               "bottom_boundary.qbot2_file required when swbotb=2 and sw2=2", &
               'bottom_boundary')
         end if
      case (3)
         call check_real_range(self%shape, 0.0_real64, 2.0_real64, &
                               'bottom_boundary.shape', errors)
         call check_real_range(self%hdrain, -1.0e4_real64, 0.0_real64, &
                               'bottom_boundary.hdrain', errors)
         call check_real_range(self%rimlay, 0.0_real64, 1.0e5_real64, &
                               'bottom_boundary.rimlay', errors)
         call check_real_range(self%aqave, -1.0e4_real64, 1.0e3_real64, &
                               'bottom_boundary.aqave', errors)
         call check_real_range(self%aqamp, 0.0_real64, 1.0e3_real64, &
                               'bottom_boundary.aqamp', errors)
         call check_real_range(self%aqper, 0.0_real64, 366.0_real64, &
                               'bottom_boundary.aqper', errors)
         call check_real_range(self%aqtmax, 0.0_real64, 366.0_real64, &
                               'bottom_boundary.aqtmax', errors)
         call check_int_enum(self%swbotb3impl, [0, 1], &
                             'bottom_boundary.swbotb3impl', errors)
         if (self%sw3 == 2 .and. .not. has_file(self%haquif_file)) then
            call errors%append(ERR_VALIDATION_REQUIRED, &
               "bottom_boundary.haquif_file required when swbotb=3 and sw3=2", &
               'bottom_boundary')
         end if
         if (self%sw4 == 1 .and. .not. has_file(self%qbot4_file)) then
            call errors%append(ERR_VALIDATION_REQUIRED, &
               "bottom_boundary.qbot4_file required when swbotb=3 and sw4=1", &
               'bottom_boundary')
         end if
      case (4)
         if (self%swqhbot == 2 .and. .not. has_file(self%qhbot_file)) then
            call errors%append(ERR_VALIDATION_REQUIRED, &
               "bottom_boundary.qhbot_file required when swbotb=4 and swqhbot=2", &
               'bottom_boundary')
         end if
      case (5)
         call check_real_range(self%hbot, -1.0e10_real64, 1.0e3_real64, &
                               'bottom_boundary.hbot', errors)
         call check_real_range(self%rhobot, -1.0e4_real64, 1.0e4_real64, &
                               'bottom_boundary.rhobot', errors)
         if (.not. has_file(self%hbot5_file)) then
            call errors%append(ERR_VALIDATION_REQUIRED, &
               "bottom_boundary.hbot5_file required when swbotb=5", &
               'bottom_boundary')
         end if
      case (6, 7, 8)
         ! No required tables.
      end select
   end subroutine bottom_boundary_config_validate

   pure function has_file(slot) result(yes)
      character(len=:), allocatable, intent(in) :: slot
      logical :: yes
      yes = allocated(slot)
      if (yes) yes = len_trim(slot) > 0
   end function has_file
```

Delete the local `check_table_2d` helper — no longer used.

Note `case (4)` previously had no validator action; we now key on `swqhbot`. Note also the old `swbotb=1 requires bbcfil/swc_table/gwl_table` rule is gone — replaced by the strict `gwl_file` requirement. The intermediate alternative slots (`swc_table`) stay in the schema but are no longer accepted as fulfillment for SWBOTB=1.

- [ ] **Step 3: Build to confirm config compiles**

Run: `pixi run -e test build`
Expected: clean.

- [ ] **Step 4: Update unit tests**

Edit `tests/unit/config/test_bottom_boundary_config.pf`:

For each test, where the old code did `c%bbcfil = "case.bbc"` → swap for `c%gwl_file = "case.gwl.csv"`. Where it set `c%gwl_table = reshape(...)` or `c%haquif_table = reshape(...)` → swap for the corresponding `c%*_file = "..."` and `c%sw3 = 2` (etc) where relevant.

Concrete changes (line numbers are approximate, locate by content):

- `test_bb_swbotb1_with_bbcfil_passes` → rename to `test_bb_swbotb1_with_gwl_file_passes`. Set `c%gwl_file = "case.gwl.csv"`. Drop `c%bbcfil`.
- `test_bb_swbotb1_with_swc_table_passes` → DELETE (alternative-slot logic is gone).
- `test_bb_swbotb1_missing_both_fails` → rename to `test_bb_swbotb1_missing_gwl_file_fails`. Sets only `c%swbotb = 1` and asserts `ERR_VALIDATION_REQUIRED`.
- `test_bb_swbotb1_both_fails` → DELETE (no longer applicable).
- `test_bb_swbotb2_with_qbot_table_passes` → rename to `test_bb_swbotb2_sw2_2_with_qbot2_file_passes`. Set `c%swbotb = 2; c%sw2 = 2; c%qbot2_file = "case.qbot.csv"`.
- `test_bb_swbotb2_missing_qbot_fails` → rename to `test_bb_swbotb2_sw2_2_missing_qbot2_file_fails`. Same setup minus `qbot2_file`. Assert `ERR_VALIDATION_REQUIRED`.
- `test_bb_swbotb3_happy_path` — strip any `c%haquif_table = ...` allocation; explicitly set `c%sw3 = 1; c%sw4 = 2` so neither sub-rule fires.
- New test `test_bb_swbotb3_sw3_2_requires_haquif_file` — sets `swbotb=3` plus all required scalars, `sw3=2`, no `haquif_file` → asserts `ERR_VALIDATION_REQUIRED`.
- New test `test_bb_swbotb3_sw3_2_with_haquif_file_passes` — same setup with `c%haquif_file = "case.haquif.csv"`.
- New test `test_bb_swbotb4_swqhbot_2_requires_qhbot_file`, `test_bb_swbotb4_swqhbot_2_with_qhbot_file_passes`.
- `test_bb_swbotb5_happy_path` → add `c%hbot5_file = "case.hbot.csv"`.
- New test `test_bb_swbotb5_missing_hbot5_file_fails`.

(Where a test references `ERR_VALIDATION_OUT_OF_RANGE` for the table-missing case, change to `ERR_VALIDATION_REQUIRED`.)

- [ ] **Step 5: Run unit tests**

Run: `pixi run -e test test-pfunit -- --filter test_bottom_boundary_config`
Expected: all PASS.

- [ ] **Step 6: Commit**

```bash
git add src/config/bottom_boundary_config.f90 tests/unit/config/test_bottom_boundary_config.pf
git commit -m "feat(bbc): add per-sub-mode *_file slots and validator rules

Adds gwl_file/qbot2_file/haquif_file/qbot4_file/qhbot_file/hbot5_file
plus sw2/sw3/sw4/swqhbot sub-mode switches. Removes legacy bbcfil and
inline gwl_table/haquif_table slots. Validator requires the matching
file per (swbotb, sub_switch) combination.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 7: BBC TOML reader — populate new slots

**Files:**
- Modify: `src/io/toml/read_bottom_boundary_toml.f90`

- [ ] **Step 1: Add reader calls for new slots**

In `read_bottom_boundary_toml`, append the sub-mode switch reads after the existing `swbotb3impl` line:

```fortran
      call get_optional_int_with_default(sec, 'sw2',     config%sw2,     1, &
                                         'bottom_boundary.sw2',     errors)
      call get_optional_int_with_default(sec, 'sw3',     config%sw3,     1, &
                                         'bottom_boundary.sw3',     errors)
      call get_optional_int_with_default(sec, 'sw4',     config%sw4,     0, &
                                         'bottom_boundary.sw4',     errors)
      call get_optional_int_with_default(sec, 'swqhbot', config%swqhbot, 1, &
                                         'bottom_boundary.swqhbot', errors)
```

Replace the existing `bbcfil` read and the two `read_date_real_table` calls (`haquif_table`, `gwl_table`) with six file slot reads:

```fortran
      call get_optional_string_with_default(sec, 'gwl_file',    config%gwl_file,    '', &
                                            'bottom_boundary.gwl_file',    errors)
      call get_optional_string_with_default(sec, 'qbot2_file',  config%qbot2_file,  '', &
                                            'bottom_boundary.qbot2_file',  errors)
      call get_optional_string_with_default(sec, 'haquif_file', config%haquif_file, '', &
                                            'bottom_boundary.haquif_file', errors)
      call get_optional_string_with_default(sec, 'qbot4_file',  config%qbot4_file,  '', &
                                            'bottom_boundary.qbot4_file',  errors)
      call get_optional_string_with_default(sec, 'qhbot_file',  config%qhbot_file,  '', &
                                            'bottom_boundary.qhbot_file',  errors)
      call get_optional_string_with_default(sec, 'hbot5_file',  config%hbot5_file,  '', &
                                            'bottom_boundary.hbot5_file',  errors)
```

Delete the `read_date_real_table` subroutine (no longer used). Also delete the `read_table_2d` calls for `swc_table`/`qbot_table`/`cofqha_table` ONLY IF the matching schema slots are still present — they are; keep those calls untouched.

- [ ] **Step 2: Update reader unit tests**

Edit `tests/unit/io/toml/test_read_bottom_boundary_toml.pf`. Any test that previously authored a `[bottom_boundary]` block with `gwl_table = [...]` or `haquif_table = [...]` or `bbcfil = "..."` must be rewritten to use `gwl_file = "..."` / `haquif_file = "..."` / `sw3 = 2` etc, and assert on the new slots. Search the file for those keywords and rework each affected test.

- [ ] **Step 3: Run unit tests**

Run: `pixi run -e test test-pfunit -- --filter test_read_bottom_boundary_toml`
Expected: all PASS.

- [ ] **Step 4: Commit**

```bash
git add src/io/toml/read_bottom_boundary_toml.f90 tests/unit/io/toml/test_read_bottom_boundary_toml.pf
git commit -m "feat(bbc-reader): read *_file slots and sub-mode switches

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 8: BBC adapter — read CSVs at adapter time

**Files:**
- Modify: `src/io/toml/config_to_variables.f90` (the `case (1)`, `case (3)`, `case (5)` blocks under `select case (swbotb)`; add new `case (2)`, `case (4)`)

- [ ] **Step 1: Replace the SWBOTB select-case body**

Replace the entire block starting at `select case (swbotb)` and ending at its `end select` (currently roughly lines 602-654) with:

```fortran
      select case (swbotb)
      case (1)
         block
            use csv_reader_mod, only: read_csv_table
            use error_mod, only: error_collection_t
            real(real64), allocatable :: csv_table(:,:)
            type(error_collection_t)  :: csv_errs
            integer :: k, nrows
            character(len=4) :: hdr(2)
            hdr(1) = 'date'
            hdr(2) = 'gwl '
            call read_csv_table( &
               trim(pathwork)//trim(config%bottom_boundary%gwl_file), &
               hdr, csv_table, csv_errs)
            call csv_errs%abort_if_fatal()
            if (allocated(csv_table)) then
               nrows = size(csv_table, 1)
               do k = 1, nrows
                  gwltab(k*2 - 1) = csv_table(k, 1)
                  gwltab(k*2)     = csv_table(k, 2)
               end do
            end if
         end block
      case (2)
         sw2 = config%bottom_boundary%sw2
         if (config%bottom_boundary%sw2 == 2) then
            block
               use csv_reader_mod, only: read_csv_table
               use error_mod, only: error_collection_t
               real(real64), allocatable :: csv_table(:,:)
               type(error_collection_t)  :: csv_errs
               integer :: k, nrows
               character(len=4) :: hdr(2)
               hdr(1) = 'date'
               hdr(2) = 'qbot'
               call read_csv_table( &
                  trim(pathwork)//trim(config%bottom_boundary%qbot2_file), &
                  hdr, csv_table, csv_errs)
               call csv_errs%abort_if_fatal()
               if (allocated(csv_table)) then
                  nrows = size(csv_table, 1)
                  do k = 1, nrows
                     qbotab(k*2 - 1) = csv_table(k, 1)
                     qbotab(k*2)     = csv_table(k, 2)
                  end do
               end if
            end block
         end if
      case (3)
         shape       = config%bottom_boundary%shape
         hdrain      = config%bottom_boundary%hdrain
         rimlay      = config%bottom_boundary%rimlay
         aqave       = config%bottom_boundary%aqave
         aqamp       = config%bottom_boundary%aqamp
         aqper       = config%bottom_boundary%aqper
         aqtmax      = config%bottom_boundary%aqtmax
         swbotb3impl = config%bottom_boundary%swbotb3impl
         sw3         = config%bottom_boundary%sw3
         sw4         = config%bottom_boundary%sw4
         if (config%bottom_boundary%sw3 == 2) then
            block
               use csv_reader_mod, only: read_csv_table
               use error_mod, only: error_collection_t
               real(real64), allocatable :: csv_table(:,:)
               type(error_collection_t)  :: csv_errs
               integer :: k, nrows
               character(len=6) :: hdr(2)
               hdr(1) = 'date  '
               hdr(2) = 'haquif'
               call read_csv_table( &
                  trim(pathwork)//trim(config%bottom_boundary%haquif_file), &
                  hdr, csv_table, csv_errs)
               call csv_errs%abort_if_fatal()
               if (allocated(csv_table)) then
                  nrows = size(csv_table, 1)
                  do k = 1, nrows
                     haqtab(k*2 - 1) = csv_table(k, 1)
                     haqtab(k*2)     = csv_table(k, 2)
                  end do
               end if
            end block
         end if
         if (config%bottom_boundary%sw4 == 1) then
            block
               use csv_reader_mod, only: read_csv_table
               use error_mod, only: error_collection_t
               real(real64), allocatable :: csv_table(:,:)
               type(error_collection_t)  :: csv_errs
               integer :: k, nrows
               character(len=4) :: hdr(2)
               hdr(1) = 'date'
               hdr(2) = 'qbot'
               call read_csv_table( &
                  trim(pathwork)//trim(config%bottom_boundary%qbot4_file), &
                  hdr, csv_table, csv_errs)
               call csv_errs%abort_if_fatal()
               if (allocated(csv_table)) then
                  nrows = size(csv_table, 1)
                  do k = 1, nrows
                     qbotab(k*2 - 1) = csv_table(k, 1)
                     qbotab(k*2)     = csv_table(k, 2)
                  end do
               end if
            end block
         end if
      case (4)
         swqhbot = config%bottom_boundary%swqhbot
         if (config%bottom_boundary%swqhbot == 2) then
            block
               use csv_reader_mod, only: read_csv_table
               use error_mod, only: error_collection_t
               real(real64), allocatable :: csv_table(:,:)
               type(error_collection_t)  :: csv_errs
               integer :: k, nrows
               character(len=4) :: hdr(2)
               hdr(1) = 'htab'
               hdr(2) = 'qtab'
               call read_csv_table( &
                  trim(pathwork)//trim(config%bottom_boundary%qhbot_file), &
                  hdr, csv_table, csv_errs)
               call csv_errs%abort_if_fatal()
               ! Legacy unpack pattern from readswap.f90:1418-1419 — for the
               ! q(h) curve, qbotab(odd) = abs(htab) and qbotab(even) = qtab.
               if (allocated(csv_table)) then
                  nrows = size(csv_table, 1)
                  do k = 1, nrows
                     qbotab(k*2 - 1) = abs(csv_table(k, 1))
                     qbotab(k*2)     = csv_table(k, 2)
                  end do
               end if
            end block
         end if
      case (5)
         hbot   = config%bottom_boundary%hbot
         rhobot = config%bottom_boundary%rhobot
         block
            use csv_reader_mod, only: read_csv_table
            use error_mod, only: error_collection_t
            real(real64), allocatable :: csv_table(:,:)
            type(error_collection_t)  :: csv_errs
            integer :: k, nrows
            character(len=4) :: hdr(2)
            hdr(1) = 'date'
            hdr(2) = 'hbot'
            call read_csv_table( &
               trim(pathwork)//trim(config%bottom_boundary%hbot5_file), &
               hdr, csv_table, csv_errs)
            call csv_errs%abort_if_fatal()
            if (allocated(csv_table)) then
               nrows = size(csv_table, 1)
               do k = 1, nrows
                  hbotab(k*2 - 1) = csv_table(k, 1)
                  hbotab(k*2)     = csv_table(k, 2)
               end do
            end if
         end block
      end select
```

Notes on legacy globals (verified in `src/core/variables.f90` and `src/io/readswap.f90`):
- `gwltab(mabbc*2)` — date/gwl interleaved (SWBOTB=1).
- `qbotab(mabbc*2)` — used by SWBOTB=2 (sw2=2) for date/qbot, SWBOTB=3 (sw4=1) added flux on top, AND SWBOTB=4 (swqhbot=2) for the q(h) curve packed as `(abs(h), q)`.
- `haqtab(mabbc*2)` — date/aquifer-head (SWBOTB=3 sw3=2).
- `hbotab(mabbc*2)` — date/hbot (SWBOTB=5).
- Sub-mode legacy switches (`sw2`, `sw3`, `sw4`, `swqhbot`) are SWAP-wide globals; assigning them at the adapter level here matches readswap.f90 lines 1308-1417.

- [ ] **Step 2: Build**

Run: `pixi run -e test build`
Expected: clean. If the build complains about unknown symbols (`qbot4tab`, etc.), adjust those lines per the verification grep above. The case-2/4/5 paths are dormant in the regression suite.

- [ ] **Step 3: Commit**

```bash
git add src/io/toml/config_to_variables.f90
git commit -m "refactor(adapter): bbc reads via read_csv_table per sub-mode

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 9: Migrate grassgrowth case data

**Files:**
- Create: `tests/swap-cases/toml/2.grassgrowth/grassgrowth.gwl.csv`
- Modify: `tests/swap-cases/toml/2.grassgrowth/swap.toml`

- [ ] **Step 1: Extract the inline `gwl_table` rows into the CSV**

Open `tests/swap-cases/toml/2.grassgrowth/swap.toml`, locate the `gwl_table = [...]` block. Each inline row has the form `[YYYY-MM-DD, value],`. Rewrite into CSV format:

```
date,gwl
1980-04-24,-88.0
1980-05-02,-121.0
...
```

(125 rows total. Use a one-shot script to extract — the inline format is `  [YYYY-MM-DD, NN.N],` so a sed/awk to strip brackets and the trailing comma is fine. Verify line count matches the inline count before and after.)

- [ ] **Step 2: Update swap.toml**

In the `[bottom_boundary]` block, replace:

```toml
swbotb = 1
gwl_table = [
  [1980-04-24,  -88.0],
  ...
]
```

with:

```toml
swbotb   = 1
gwl_file = "grassgrowth.gwl.csv"
```

- [ ] **Step 3: Run grassgrowth regression**

Run: `pixi run -e test regression -- -k grassgrowth`
Expected: PASS at 1e-2 cm tolerance.

- [ ] **Step 4: Commit**

```bash
git add tests/swap-cases/toml/2.grassgrowth/grassgrowth.gwl.csv tests/swap-cases/toml/2.grassgrowth/swap.toml
git commit -m "test(grassgrowth): migrate gwl_table to gwl_file CSV companion

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 10: Migrate oxygenstress case data

**Files:**
- Create: `tests/swap-cases/toml/4.oxygenstress/oxygenstress.haquif.csv`
- Modify: `tests/swap-cases/toml/4.oxygenstress/swap.toml`

- [ ] **Step 1: Extract `haquif_table` into CSV**

Header: `date,haquif`. Rows: each inline row.

- [ ] **Step 2: Update swap.toml**

In `[bottom_boundary]`, replace `haquif_table = [...]` with:

```toml
sw3         = 2
haquif_file = "oxygenstress.haquif.csv"
```

(Keep `swbotb = 3`, `shape`, `hdrain`, `rimlay`, `swbotb3impl = 1` — they remain.)

- [ ] **Step 3: Run regression**

Run: `pixi run -e test regression -- -k oxygenstress`
Expected: PASS (no collateral; same numeric content via CSV).

- [ ] **Step 4: Commit**

```bash
git add tests/swap-cases/toml/4.oxygenstress/
git commit -m "test(oxygenstress): migrate haquif_table to haquif_file CSV companion

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 11: Migrate surfacewater case data

**Files:**
- Create: `tests/swap-cases/toml/6.surfacewater/surfacewater.haquif.csv`
- Modify: `tests/swap-cases/toml/6.surfacewater/swap.toml`

Same procedure as Task 10 for the surfacewater inline `haquif_table`. Add `sw3 = 2` and `haquif_file = "surfacewater.haquif.csv"`; remove the inline table.

- [ ] **Step 1: Extract CSV.**
- [ ] **Step 2: Update swap.toml.**
- [ ] **Step 3: Run regression** — `pixi run -e test regression -- -k surfacewater` — PASS.
- [ ] **Step 4: Commit**

```bash
git add tests/swap-cases/toml/6.surfacewater/
git commit -m "test(surfacewater): migrate haquif_table to haquif_file CSV companion

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 12: Update documentation

**Files:**
- Modify: `docs/csv-companion-files.md`

- [ ] **Step 1: Rewrite the doc**

Replace existing content of `docs/csv-companion-files.md` with content reflecting:

1. Reader contract — single sub `read_csv_table(path, expected_header, table, errors)`. Header is required, validated strict-positional case-insensitive, comma-separated, trimmed.
2. Authoring rules — column 1 may be `date` (ISO) or any real-keyed name. Columns 2..N are reals. `#` and blank lines tolerated.
3. Schema — list the 6 BBC slots and the irrigation slot:

| Slot | Section | Header | Required when |
|---|---|---|---|
| `gwl_file` | `[bottom_boundary]` | `date,gwl` | `swbotb = 1` |
| `qbot2_file` | `[bottom_boundary]` | `date,qbot` | `swbotb = 2`, `sw2 = 2` |
| `haquif_file` | `[bottom_boundary]` | `date,haquif` | `swbotb = 3`, `sw3 = 2` |
| `qbot4_file` | `[bottom_boundary]` | `date,qbot` | `swbotb = 3`, `sw4 = 1` |
| `qhbot_file` | `[bottom_boundary]` | `htab,qtab` | `swbotb = 4`, `swqhbot = 2` |
| `hbot5_file` | `[bottom_boundary]` | `date,hbot` | `swbotb = 5` |
| `fixed_events_file` | `[irrigation]` | `date,depth,conc,type` | `swirfix = 1` (long-form) |

4. Errors — `ERR_IO_OPEN_FAILED`, `ERR_PARSE_MISSING_HEADER`, `ERR_PARSE_HEADER_MISMATCH`, `ERR_PARSE_ROW_SHAPE`, `ERR_PARSE_TYPE_MISMATCH`, `ERR_VALIDATION_REQUIRED`.
5. Path resolution — basename in working directory, matching `*.crp.toml`. Staging — regression / `run_case.sh` glob `*.csv` and copy alongside `swap.toml`.

- [ ] **Step 2: Commit**

```bash
git add docs/csv-companion-files.md
git commit -m "docs: refresh csv-companion-files for unified reader

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 13: Full verification

- [ ] **Step 1: Run unit suite**

Run: `pixi run -e test test-pfunit`
Expected: green.

- [ ] **Step 2: Run full regression**

Run: `pixi run -e test check-full`
Expected: hupselbrook, grassgrowth, oxygenstress, salinitystress, surfacewater all PASS at 1e-2 cm tolerance. Pre-existing oxygenstress / salinitystress drift unchanged.

- [ ] **Step 3: Search for stragglers**

Run: `git grep -n 'gwl_table\|haquif_table\|bbcfil\|read_csv_date_reals' -- src tests`
Expected: zero matches in src/, zero in tests/swap-cases/toml/. Any matches in committed artifacts (other than this plan and the spec) are bugs to fix before claiming done.

- [ ] **Step 4: Final commit if any cleanup needed**

```bash
git add -A
git status
# only proceed if changes are intentional cleanup
git commit -m "chore: post-migration sweep"
```
