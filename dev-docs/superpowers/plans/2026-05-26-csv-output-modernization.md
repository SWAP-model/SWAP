# CSV Output Modernization (IO-OUT) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace SWAP's half-retired output layer with a symmetric, registry-driven CSV path — a reusable `csv_writer` primitive, a declarative `output_registry`, one `csv_output` module that samples `state` once for both the two result CSVs and the live C-API streams — and delete the `swapoutput.f90` graveyard.

**Architecture:** `csv_writer.f90` mirrors `csv_reader.f90` (numeric matrix out). `output_registry.f90` holds the `{name, unit, kind}` catalogue and resolves the user `inlist`. `csv_output.f90` (renamed from `swap_csv_output.f90`) keeps centralized sampling + aggregation, writes through `csv_writer`, and fills the live BMI rows. `swapoutput.f90` is deleted; the four inert C-API streams (`temperature/solute/surfacewater/agetracer`) are dropped.

**Tech Stack:** Fortran 2008 (gfortran), Meson build, pFUnit unit tests, Python regression harness (`tests/regression/test_output_regression.py`). Build/test via `pixi run -e test <task>`.

**Spec:** `dev-docs/superpowers/specs/2026-05-26-csv-output-modernization-design.md`

---

## Conventions for every task

- **Build:** `pixi run -e test build-linux`
- **Fast gate:** `pixi run -e test check-fast` (build + pFUnit + 4 regression cases: hupselbrook, surfacewater, salinitystress, grassgrowth)
- **Full gate:** `pixi run -e test check-full` (build + pFUnit + all 5 cases)
- **Single pFUnit run:** `pixi run -e test test-pfunit` (whole unit-pfunit suite; there is no per-test target)
- **pFUnit registration (load-bearing):** a new `tests/unit/*.pf` file must be added in THREE places or it compiles but never runs ("dark"): (1) the source-list and (2) the `pf_files` list in `tests/unit/meson.build`, AND (3) an `ADD_TEST_SUITE(test_<name>_suite)` line in `tests/unit/testSuites.inc`. After adding a test, confirm the printed pFUnit total (`OK (N tests)`) actually increased.
- **Commits:** `development` branch only, one focused change each, message tagged `IO-OUT/<phase>`.
- **State-schema rebuild:** any task that edits `src/state/*_state.f90` MUST `rm -rf builddir` before building (Meson does not propagate `.mod` deps across the swap_modern↔swap_legacy boundary). Only Phase E does this.
- **Format invariant:** the CSV number formatter (E-notation when `|x| < 1e-4` or `> 1e4`, else fixed 5-decimal) is load-bearing for the `1e-2` regression tolerance. It must be carried verbatim, never re-derived.

---

## File Structure

| File | Responsibility | Phase |
|------|----------------|-------|
| `src/io/csv_writer.f90` (new) | Generic CSV writer primitive: open / meta / header / row / flush / close + `fmt_real`. Symmetric to `csv_reader.f90`. | A |
| `tests/unit/io/test_csv_writer.pf` (new) | Round-trip (write→read-back) + formatter golden-string tests. | A |
| `src/io/output_registry.f90` (new) | `out_var_t` catalogue + `resolve_inlist` + accessors. Single source for value-array layout and header. | B |
| `tests/unit/io/test_output_registry.pf` (new) | inlist resolution: filter/order, unknown-name error, default fallback. | B |
| `src/io/csv_output.f90` (renamed from `swap_csv_output.f90`) | Orchestrates init/step/finalize; aggregation; sampling; writes via `csv_writer`; fills live BMI rows. Merges the `_tz` module. | C, E |
| `tests/unit/io/test_csv_aggregates.pf` (new) | `dstor`/`baldev` math against a constructed state. | C |
| `scripts/output_parity_check.py` (new) | Proves `.crp`/`.snw` columns are reproduced in `result_output.csv` within tol. | D |
| `src/io/swapoutput.f90` (deleted) | — | E |
| `src/core/swap_mod.f90` (modify) | Replace `*Output` task calls with `csv_output_init/step/finalize`. | E |
| `src/core/swap_capi_mod.f90` (modify) | Drop 4 inert stream cases. | E |
| `src/crop/cropgrowth_helpers.f90` (modify) | Drop `OutCropFixed/OutWofost/OutGrass` calls. | E |
| `src/state/{heat,solute,surfacewater}_state.f90` (modify) | Remove the dropped streams' `output_row`/`*_columns`/`*_n_cols` fields. | E |
| `meson.build`, `tests/unit/meson.build` (modify) | Register new sources/tests; drop `swapoutput.f90`. | A–E |

---

## Phase A — `csv_writer` primitive

Standalone new module; nothing consumes it yet. De-risks every later phase.

### Task A1: Create `csv_writer` module and round-trip test

**Files:**
- Create: `src/io/csv_writer.f90`
- Create: `tests/unit/io/test_csv_writer.pf`
- Modify: `meson.build` (add to `swap_io_sources`, near line 117 next to `csv_reader.f90`)
- Modify: `tests/unit/meson.build` (add source near line 104; add test near line 233)

- [ ] **Step 1: Write the failing round-trip test**

Create `tests/unit/io/test_csv_writer.pf`:

```fortran
! Tests for csv_writer_mod — symmetric counterpart of csv_reader_mod.

@test
subroutine test_roundtrip_numeric()
   use funit
   use iso_fortran_env, only: real64
   use csv_writer_mod,  only: csv_writer_t
   use csv_reader_mod,  only: read_csv_table
   use error_mod,       only: error_collection_t
   type(csv_writer_t)        :: w
   real(real64), allocatable :: table(:,:)
   type(error_collection_t)  :: errors
   character(len=*), parameter :: path = 'test_csv_writer_rt.csv'

   ! write a known 3x2 matrix (no meta, no units -> reader can parse it back)
   call w%open(path, errors)
   @assertFalse(errors%has_errors())
   call w%header(['htab', 'qtab'])
   call w%row([-100.0_real64, 0.0_real64])
   call w%row([ -50.0_real64, 1.5_real64])
   call w%row([   0.0_real64, 9.0_real64])
   call w%close()

   call read_csv_table(path, ['htab', 'qtab'], table, errors)
   @assertFalse(errors%has_errors())
   @assertEqual(3, size(table, 1))
   @assertEqual(2, size(table, 2))
   @assertEqual(-100.0_real64, table(1, 1), 1.0e-4_real64)
   @assertEqual(   1.5_real64, table(2, 2), 1.0e-4_real64)
   @assertEqual(   9.0_real64, table(3, 2), 1.0e-4_real64)
end subroutine

@test
subroutine test_leading_string_and_trailing_comma()
   use funit
   use iso_fortran_env, only: real64
   use csv_writer_mod,  only: csv_writer_t
   use error_mod,       only: error_collection_t
   type(csv_writer_t)       :: w
   type(error_collection_t) :: errors
   character(len=*), parameter :: path = 'test_csv_writer_lead.csv'
   integer :: u, ios
   character(len=256) :: hdr, rec

   call w%open(path, errors)
   call w%header(['DATETIME', 'RAIN    '])
   call w%row([12.5_real64], leading='1980-04-24')
   call w%close()

   open(newunit=u, file=path, status='old', action='read')
   read(u, '(A)', iostat=ios) hdr   ! header line
   read(u, '(A)', iostat=ios) rec   ! data line
   close(u)
   ! no trailing comma, leading datetime present
   @assertEqual('1980-04-24,12.50000', trim(adjustl(rec)))
end subroutine

@test
subroutine test_fmt_real_switches_to_exponential()
   use funit
   use iso_fortran_env, only: real64
   use csv_writer_mod,  only: fmt_real_for_test
   ! |x|<1e-4 -> E-notation; mid-range -> 5-decimal F.
   @assertTrue(index(trim(fmt_real_for_test(1.0e-6_real64)), 'E') > 0)
   @assertTrue(index(trim(fmt_real_for_test(1.0e6_real64)),  'E') > 0)
   @assertEqual('1.23450', trim(adjustl(fmt_real_for_test(1.2345_real64))))
end subroutine
```

- [ ] **Step 2: Register the new source + test, then run to verify failure**

In `tests/unit/meson.build` add after line 104 (`'../../src/io/csv_reader.f90',`):
```
        '../../src/io/csv_writer.f90',
```
and after line 233 (`'io/test_csv_reader.pf',`):
```
        'io/test_csv_writer.pf',
```

Run: `pixi run -e test test-pfunit`
Expected: FAIL — `Cannot open module file 'csv_writer_mod'` / unresolved `csv_writer_t`.

- [ ] **Step 3: Implement `csv_writer.f90`**

Create `src/io/csv_writer.f90`:

```fortran
!> Generic CSV writer for SWAP result files — symmetric counterpart of
!! csv_reader_mod. A csv_writer_t owns one open unit. The number formatter
!! (fmt_real) preserves the legacy E/F selection so byte output is stable.
module csv_writer_mod
   use iso_fortran_env, only: real64
   use error_mod,       only: error_collection_t, ERR_IO_OPEN_FAILED
   use file_io_mod,     only: file_open
   implicit none
   private

   public :: csv_writer_t
   public :: fmt_real_for_test   ! thin alias exported only for unit tests

   type :: csv_writer_t
      integer :: unit  = -1
      integer :: ncols = 0
   contains
      procedure :: open   => csv_writer_open
      procedure :: meta   => csv_writer_meta
      procedure :: header => csv_writer_header
      procedure :: row    => csv_writer_row
      procedure :: flush  => csv_writer_flush
      procedure :: close  => csv_writer_close
   end type csv_writer_t

contains

   subroutine csv_writer_open(self, path, errors)
      class(csv_writer_t),      intent(inout) :: self
      character(len=*),         intent(in)    :: path
      type(error_collection_t), intent(inout) :: errors
      integer :: ios
      call file_open(self%unit, path, 'unknown', 'readwrite', iostat=ios)
      if (ios /= 0) then
         call errors%append(ERR_IO_OPEN_FAILED, &
            "cannot open CSV for write: " // trim(path), 'csv_writer')
         self%unit = -1
      end if
   end subroutine csv_writer_open

   !> Write '*'-prefixed metadata comment lines (Project, File content, ...).
   subroutine csv_writer_meta(self, lines)
      class(csv_writer_t), intent(inout) :: self
      character(len=*),    intent(in)    :: lines(:)
      integer :: i
      do i = 1, size(lines)
         write(self%unit, '(A)') '* ' // trim(lines(i))
      end do
   end subroutine csv_writer_meta

   !> Write the column-name row and (optionally) a unit row. Sets ncols.
   subroutine csv_writer_header(self, names, units)
      class(csv_writer_t), intent(inout) :: self
      character(len=*),    intent(in)    :: names(:)
      character(len=*),    intent(in), optional :: units(:)
      self%ncols = size(names)
      write(self%unit, '(A)') join(names)
      if (present(units)) write(self%unit, '(A)') join(units)
   end subroutine csv_writer_header

   !> Write one data row: optional leading string (datetime), then values.
   subroutine csv_writer_row(self, values, leading)
      class(csv_writer_t), intent(inout) :: self
      real(real64),        intent(in)    :: values(:)
      character(len=*),    intent(in), optional :: leading
      character(len=:), allocatable :: line
      integer :: j
      line = ''
      if (present(leading)) line = trim(leading) // ','
      do j = 1, size(values)
         line = line // trim(adjustl(fmt_real(values(j))))
         if (j < size(values)) line = line // ','
      end do
      write(self%unit, '(A)') line
   end subroutine csv_writer_row

   !> Drain the runtime I/O buffer to disk without closing. Called at year
   !! boundaries so a crash mid-run leaves a near-complete file. No-op if
   !! the writer is not open (e.g. headless mode never opened a unit).
   subroutine csv_writer_flush(self)
      class(csv_writer_t), intent(inout) :: self
      if (self%unit /= -1) flush(self%unit)
   end subroutine csv_writer_flush

   subroutine csv_writer_close(self)
      class(csv_writer_t), intent(inout) :: self
      if (self%unit /= -1) close(self%unit)
      self%unit = -1
   end subroutine csv_writer_close

   pure function join(fields) result(s)
      character(len=*), intent(in) :: fields(:)
      character(len=:), allocatable :: s
      integer :: i
      s = ''
      do i = 1, size(fields)
         s = s // trim(adjustl(fields(i)))
         if (i < size(fields)) s = s // ','
      end do
   end function join

   !> Legacy E/F selector lifted verbatim from swap_csv_output.f90's `what_form`
   !! (thresholds t1=1.0d-4, t2=1.0d4, num_d=5 decimals). DO NOT re-derive.
   function fmt_real(x) result(buf)
      real(real64), intent(in) :: x
      character(len=30) :: buf
      character(len=20) :: form
      real(real64), parameter :: t1 = 1.0d-4, t2 = 1.0d4
      if (x /= 0.0d0 .and. (abs(x) < t1 .or. abs(x) > t2)) then
         form = '(1pE12.5)'
      else
         form = '(F12.5)'
      end if
      write(buf, form) x
   end function fmt_real

   function fmt_real_for_test(x) result(buf)
      real(real64), intent(in) :: x
      character(len=30) :: buf
      buf = fmt_real(x)
   end function fmt_real_for_test

end module csv_writer_mod
```

> NOTE: before committing, open `swap_csv_output.f90:473` (`what_form`) and confirm `fmt_real` above reproduces its exact `form` strings and thresholds for representative values. If `what_form` uses different widths/parameters, copy them verbatim — the regression tolerance depends on it.

- [ ] **Step 4: Add the source to the production build**

In `meson.build` after line 117 (`'src/io/csv_reader.f90',`):
```
    'src/io/csv_writer.f90',
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `pixi run -e test test-pfunit`
Expected: PASS, including the three `test_csv_writer` cases. (If `test_fmt_real_switches_to_exponential` or `test_leading_string_and_trailing_comma` fail, reconcile `fmt_real` with `what_form` and the trailing-comma rule — the legacy writer strips the final comma; `csv_writer_row` achieves this by not appending after the last value.)

- [ ] **Step 6: Commit**

```bash
git add src/io/csv_writer.f90 tests/unit/io/test_csv_writer.pf meson.build tests/unit/meson.build
git commit -m "feat(io): IO-OUT/A — csv_writer primitive symmetric to csv_reader"
```

---

## Phase B — `output_registry`

Extract the hardcoded variable catalogue into a declarative module. No file-format change; `result_output.csv` byte output must stay identical.

### Task B1: Create `output_registry` with `resolve_inlist` + test

**Files:**
- Create: `src/io/output_registry.f90`
- Create: `tests/unit/io/test_output_registry.pf`
- Modify: `meson.build`, `tests/unit/meson.build`

- [ ] **Step 1: Read the current catalogue**

Open `src/io/swap_csv_output.f90:94-198` and copy the full ordered list of variable names + units + their kind (scalar / per-node template `X[` / per-subregion template `Y[`). This list is the source of truth for the registry order; preserving order keeps the value-array indices stable.

- [ ] **Step 2: Write the failing test**

Create `tests/unit/io/test_output_registry.pf`:

```fortran
! Tests for output_registry_mod.

@test
subroutine test_var_count_nonzero()
   use funit
   use output_registry_mod, only: var_count
   @assertTrue(var_count() > 0)
end subroutine

@test
subroutine test_resolve_inlist_filters_and_orders()
   use funit
   use output_registry_mod, only: resolve_inlist
   use error_mod,           only: error_collection_t
   integer, allocatable     :: sel(:)
   type(error_collection_t) :: errors
   ! RAIN and GWL both exist; selection returns their registry indices in
   ! canonical (registry) order regardless of inlist order.
   call resolve_inlist('GWL, RAIN', sel, errors)
   @assertFalse(errors%has_errors())
   @assertEqual(2, size(sel))
   @assertTrue(sel(1) < sel(2))   ! canonical order, not inlist order
end subroutine

@test
subroutine test_resolve_inlist_unknown_name_errors()
   use funit
   use output_registry_mod, only: resolve_inlist
   use error_mod,           only: error_collection_t
   integer, allocatable     :: sel(:)
   type(error_collection_t) :: errors
   call resolve_inlist('RAIN, NOTAVAR', sel, errors)
   @assertTrue(errors%has_errors())
end subroutine
```

> Adjust `RAIN`/`GWL` if those exact names are not in the catalogue from Step 1; pick any two scalar names that are, where the first listed appears later in the registry than the second.

- [ ] **Step 3: Register source + test; run to verify failure**

`tests/unit/meson.build`: add `'../../src/io/output_registry.f90',` after the `csv_writer.f90` line, and `'io/test_output_registry.pf',` after the `test_csv_writer.pf` line.

Run: `pixi run -e test test-pfunit`
Expected: FAIL — `output_registry_mod` not found.

- [ ] **Step 4: Implement `output_registry.f90`**

Create `src/io/output_registry.f90`:

```fortran
!> Declarative catalogue of SWAP output variables. The module-level registry
!! is the single source of truth for both the value-array layout (registry
!! index == value index) and the emitted CSV header/units.
module output_registry_mod
   use error_mod, only: error_collection_t, ERR_PARSE_TYPE_MISMATCH
   implicit none
   private

   integer, parameter, public :: OUT_SCALAR = 1, OUT_NODE = 2, &
                                  OUT_SUBREGION = 3, OUT_PROFILE = 4

   type, public :: out_var_t
      character(len=32) :: name = ''
      character(len=16) :: unit = ''
      integer           :: kind = OUT_SCALAR
   end type out_var_t

   public :: var_count, var_name, var_unit, var_kind, resolve_inlist

   ! Catalogue, in canonical order. Ported verbatim from the former inline
   ! list in swap_csv_output.f90:94-198. (Fill from Step 1.)
   type(out_var_t), parameter :: REGISTRY(*) = [ &
      out_var_t('RAIN',  'cm',   OUT_SCALAR), &
      out_var_t('EACT',  'cm',   OUT_SCALAR), &
      ! ... every remaining entry from the source list, same order ...
      out_var_t('GWL',   'cm',   OUT_SCALAR)  &
   ]

contains

   pure integer function var_count()
      var_count = size(REGISTRY)
   end function var_count

   pure function var_name(i) result(s)
      integer, intent(in) :: i
      character(len=32)   :: s
      s = REGISTRY(i)%name
   end function var_name

   pure function var_unit(i) result(s)
      integer, intent(in) :: i
      character(len=16)   :: s
      s = REGISTRY(i)%unit
   end function var_unit

   pure integer function var_kind(i)
      integer, intent(in) :: i
      var_kind = REGISTRY(i)%kind
   end function var_kind

   !> Parse a comma-separated inlist, match each token (its base name, before
   !! any '[') against the registry, and return matched indices in registry
   !! order. Unknown tokens append a fatal error.
   subroutine resolve_inlist(inlist, sel, errors)
      character(len=*),         intent(in)    :: inlist
      integer, allocatable,     intent(out)   :: sel(:)
      type(error_collection_t), intent(inout) :: errors
      logical :: want(size(REGISTRY))
      character(len=:), allocatable :: rest, tok, base
      integer :: comma, br, i, n
      want = .false.
      rest = adjustl(inlist)
      do while (len_trim(rest) > 0)
         comma = index(rest, ',')
         if (comma > 0) then
            tok = trim(adjustl(rest(1:comma-1))); rest = adjustl(rest(comma+1:))
         else
            tok = trim(adjustl(rest)); rest = ''
         end if
         if (len_trim(tok) == 0) cycle
         br = index(tok, '[')                  ! strip [..] node/subregion selector
         if (br > 0) then; base = tok(1:br-1); else; base = tok; end if
         n = match(base)
         if (n == 0) then
            call errors%append(ERR_PARSE_TYPE_MISMATCH, &
               "unknown output variable: '" // trim(base) // "'", 'output_registry')
         else
            want(n) = .true.
         end if
      end do
      sel = pack([(i, i=1,size(REGISTRY))], want)
   end subroutine resolve_inlist

   pure integer function match(base)
      character(len=*), intent(in) :: base
      integer :: i
      match = 0
      do i = 1, size(REGISTRY)
         if (trim(REGISTRY(i)%name) == trim(base)) then
            match = i; return
         end if
      end do
   end function match

end module output_registry_mod
```

> The `REGISTRY(*)` literal must contain **every** entry from Step 1 in the same order. For per-node/per-subregion entries keep the trailing `[` in the stored name to match how `swap_csv_output` keys them; `resolve_inlist` strips `[..]` from the user token before matching, so store names without the `[` OR adjust `match` to compare against the base — pick one and be consistent with how `csv_output` looks names up in Phase C.

- [ ] **Step 5: Add to production build**

`meson.build`: add `'src/io/output_registry.f90',` after the `csv_writer.f90` line.

- [ ] **Step 6: Run unit tests**

Run: `pixi run -e test test-pfunit`
Expected: PASS for the three registry tests.

- [ ] **Step 7: Commit**

```bash
git add src/io/output_registry.f90 tests/unit/io/test_output_registry.pf meson.build tests/unit/meson.build
git commit -m "feat(io): IO-OUT/B — output_registry catalogue + inlist resolver"
```

### Task B2: Make `swap_csv_output` consume the registry

**Scope note (refined during execution):** `make_userlist`/`merge`/`det_which_vars` do more than name-matching — they expand the `shorts` aliases (`WATBAL`→17 vars, etc.) and parse node/sub-region selectors (`H[-10.0,4,5]`). `resolve_inlist` knows neither. So B2 does NOT replace that machinery. B2 only single-sources the *catalogue* (`vars%name`/`vars%unit`) from the registry, deleting the duplicated `data` block. The selection machinery stays and continues to read `vars%name`. `resolve_inlist` remains tested registry infrastructure (alias/selector-aware validation is a separate future concern). Output must stay byte-identical.

**Files:**
- Modify: `src/io/swap_csv_output.f90` (replace the `data` catalogue block with a runtime registry fill)

- [ ] **Step 1: Replace the inline `data` catalogue with a registry fill**

In `swap_csv_output.f90`: add `use output_registry_mod, only: var_count, var_name, var_unit`. Delete the `data (vars%name(i), vars%unit(i), i = 1, M) / … /` block (lines 94-198). In `csv_out` `case (1)`, before `make_userlist`, fill from the registry:

```fortran
if (var_count() /= M) call fatalerr_collected('csv_out', 'registry size /= M')
do i = 1, M
   vars%name(i) = var_name(i)   ! char(32) -> char(12): all names <= 12 chars
   vars%unit(i) = var_unit(i)   ! char(16) -> char(12): all units <= 12 chars
end do
```

`makeheader`, `merge`, `det_which_vars`, and the `shorts` alias table are unchanged — they keep reading `vars%name`/`vars%unit`, now sourced from the single registry.

- [ ] **Step 2: Build**

Run: `pixi run -e test build-linux`
Expected: compiles clean.

- [ ] **Step 3: Regression must be unchanged**

Run: `pixi run -e test check-fast`
Expected: 4/4 cases `regression ok`. If any case drifts, the registry order or a unit string diverged from the inline list — diff `result_output.csv` header against a pre-change run and fix the registry entry.

- [ ] **Step 4: Commit**

```bash
git add src/io/swap_csv_output.f90
git commit -m "refactor(io): IO-OUT/B — swap_csv_output reads catalogue from registry"
```

---

## Phase C — route writing through `csv_writer`, rename to `csv_output`

### Task C1: Extract `dstor`/`baldev` into a tested pure helper

**Files:**
- Modify: `src/io/swap_csv_output.f90` (factor `fill_values`' balance math into a pure subroutine)
- Create: `tests/unit/io/test_csv_aggregates.pf`

- [ ] **Step 1: Write the failing test**

Create `tests/unit/io/test_csv_aggregates.pf`:

```fortran
! Tests for the water-balance aggregates extracted from fill_values.

@test
subroutine test_dstor_and_baldev()
   use funit
   use iso_fortran_env, only: real64
   use csv_aggregates_mod, only: water_balance_dev
   real(real64) :: dstor, baldev
   ! storage change = (volact+pond+ssnow) - (old sums)
   ! inputs - dstor - outputs should net to baldev.
   call water_balance_dev( &
      volact=10.0_real64, pond=1.0_real64, ssnow=0.0_real64, &
      vol_old=9.0_real64, pond_old=1.0_real64, snow_old=0.0_real64, &
      rain=2.0_real64, irrig=0.0_real64, runon=0.0_real64, qssdi=0.0_real64, &
      interc=0.2_real64, runoff=0.0_real64, runoffcn=0.0_real64, &
      tact=0.5_real64, eact=0.3_real64, sublim=0.0_real64, &
      drainage=0.0_real64, qbottom=0.0_real64, &
      dstor=dstor, baldev=baldev)
   @assertEqual(1.0_real64, dstor, 1.0e-9_real64)
   @assertEqual(0.0_real64, baldev, 1.0e-9_real64)
end subroutine
```

- [ ] **Step 2: Register + run to verify failure**

`tests/unit/meson.build`: add `'io/test_csv_aggregates.pf',` after the registry test line. (The module lives inside `swap_csv_output.f90`, already a unit source.)

Run: `pixi run -e test test-pfunit`
Expected: FAIL — `csv_aggregates_mod` not found.

- [ ] **Step 3: Implement the pure helper**

In `swap_csv_output.f90`, add a small module `csv_aggregates_mod` (or a public pure subroutine in the existing module) containing exactly the `dstor`/`baldev` formulas currently in `fill_values` (`swap_csv_output.f90:509-572`):

```fortran
pure subroutine water_balance_dev(volact, pond, ssnow, vol_old, pond_old, &
      snow_old, rain, irrig, runon, qssdi, interc, runoff, runoffcn, tact, &
      eact, sublim, drainage, qbottom, dstor, baldev)
   use iso_fortran_env, only: real64
   real(real64), intent(in)  :: volact, pond, ssnow, vol_old, pond_old, snow_old
   real(real64), intent(in)  :: rain, irrig, runon, qssdi, interc, runoff
   real(real64), intent(in)  :: runoffcn, tact, eact, sublim, drainage, qbottom
   real(real64), intent(out) :: dstor, baldev
   dstor  = (volact + pond + ssnow) - (vol_old + pond_old + snow_old)
   baldev = (rain + irrig + runon + qssdi) - dstor &
          - (interc + runoff + runoffcn + tact + eact + sublim + drainage + qbottom)
end subroutine water_balance_dev
```

Then call `water_balance_dev` from `fill_values` instead of the inline arithmetic.

> Confirm the exact term list against `fill_values` lines 509-572 before finalizing — the spec quotes this formula but the live code is authoritative.

- [ ] **Step 4: Run tests**

Run: `pixi run -e test test-pfunit` then `pixi run -e test check-fast`
Expected: aggregate test PASS; 4/4 regression unchanged.

- [ ] **Step 5: Commit**

```bash
git add src/io/swap_csv_output.f90 tests/unit/io/test_csv_aggregates.pf tests/unit/meson.build
git commit -m "refactor(io): IO-OUT/C — extract tested water_balance_dev helper"
```

**C2 split (refined during execution).** The original single C2 (rename + merge two modules + route both writers through `csv_writer` + new API + rewire) is too large a byte-sensitive blast radius for one commit. It is split into **C2a** (mechanical rename + named API + rewire; inline writes unchanged) and **C2b** (route writes through `csv_writer` + flush). The **physical module merge is dropped** as unnecessary risk: `module csv_output` owns the public API and delegates the time-depth output to a renamed `module csv_output_tz` helper in the same file. **Gate reality:** the regression harness keys columns by NAME and compares annual stats within `1e-2` — it does NOT enforce byte-identity. Aim for byte-identical output, but the acceptance criterion is check-fast 4/4 with the same column NAMES present; diff the header + first data rows against a pre-change baseline to catch cosmetic drift.

### Task C2a: Rename to `csv_output` + named `init/step/finalize`; rewire caller

**Files:**
- Rename: `src/io/swap_csv_output.f90` → `src/io/csv_output.f90`; rename module `SWAP_csv_output` → `csv_output` and `SWAP_csv_output_tz` → `csv_output_tz` (both stay in the file; NOT merged).
- Modify: `src/io/swapoutput.f90` (`soilwateroutput` calls the new procedures).
- Modify: `meson.build` (source path), `tests/unit/meson.build` (the unit suite lists this source — update the path).

- [ ] **Step 1: Capture a baseline output for later diffing**

Build, run one case (e.g. hupselbrook) to produce a `result_output.csv` and (if enabled) `result_output_tz.csv`, and copy them aside (e.g. `/tmp/baseline_output.csv`). These are the byte-comparison reference for C2a and C2b.

- [ ] **Step 2: `git mv` + rename both modules**

`git mv src/io/swap_csv_output.f90 src/io/csv_output.f90`. In the file rename `module SWAP_csv_output` → `module csv_output` (and its `end module`), and `module SWAP_csv_output_tz` → `module csv_output_tz`. Update `meson.build` (`'src/io/swap_csv_output.f90'` → `'src/io/csv_output.f90'`) and `tests/unit/meson.build` (same path in the unit source list). Do NOT merge the two module bodies. Do NOT change any procedure body.

- [ ] **Step 3: Add named entry points in `module csv_output`**

Add three public procedures that wrap the existing task bodies, reading the enable flags from `state%cfg%output_csv` (as `csv_out` already reads `state%cfg%output_csv%inlist`). `csv_output_tz` is reached via `use csv_output_tz, only: csv_out_tz`:

```fortran
public :: csv_output_init, csv_output_step, csv_output_finalize
...
subroutine csv_output_init(state)
   type(swap_state_t), intent(inout) :: state
   if (state%cfg%output_csv%enabled    == 1) call csv_out(1, state)
   if (state%cfg%output_csv%enabled_tz == 1) call csv_out_tz(1, state)
end subroutine
subroutine csv_output_step(state)
   type(swap_state_t), intent(inout) :: state
   if (state%cfg%output_csv%enabled    == 1) call csv_out(2, state)
   if (state%cfg%output_csv%enabled_tz == 1) call csv_out_tz(2, state)
end subroutine
subroutine csv_output_finalize(state)
   type(swap_state_t), intent(inout) :: state
   if (state%cfg%output_csv%enabled    == 1) call csv_out(3, state)
   if (state%cfg%output_csv%enabled_tz == 1) call csv_out_tz(3, state)
end subroutine
```

Keep `csv_out`/`csv_out_tz` public for now (E-phase removes external callers). Add `use csv_output_tz, only: csv_out_tz` to `module csv_output`.

- [ ] **Step 4: Rewire the dispatch hub**

In `swapoutput.f90` `soilwateroutput`: replace `use SWAP_csv_output` + `use SWAP_csv_output_tz` with `use csv_output, only: csv_output_init, csv_output_step, csv_output_finalize`. Replace the task-1 body's `csv_out(1)`/`csv_out_tz(1)` calls with `call csv_output_init(state)`, the task-2 body with `call csv_output_step(state)`, and the task-4 body's `csv_out(3)`/`csv_out_tz(3)` with `call csv_output_finalize(state)`. (The enable-flag `if`s now live inside the wrappers, so drop them from `soilwateroutput`.)

- [ ] **Step 5: Build + regression + byte-diff**

Run `pixi run -e test check-fast` → 4/4. Then re-run the baseline case and `diff` its `result_output.csv` against `/tmp/baseline_output.csv` → expect ZERO diff (pure rename, no write-path change). Confirm `pixi run -e test test-pfunit` still `OK (762 tests)`.

- [ ] **Step 6: Commit**

```bash
git add src/io/csv_output.f90 src/io/swapoutput.f90 meson.build tests/unit/meson.build
git commit -m "refactor(io): IO-OUT/C — rename swap_csv_output->csv_output; named init/step/finalize"
```

### Task C2b: Route writes through `csv_writer_t` + per-year flush

**Files:**
- Modify: `src/io/csv_output.f90` (both `csv_out` and `csv_out_tz` write paths)

- [ ] **Step 1: Route `csv_out`'s scalar write through `csv_writer_t`**

Add `use csv_writer_mod, only: csv_writer_t` to `module csv_output`; add a module-level `type(csv_writer_t), save :: scalar_w`. In `csv_out` `case(1)`, replace the `file_open` + `makeheader` block with `scalar_w%open(filcsv, errors)` + `scalar_w%meta(meta_lines)` + `scalar_w%header(names, units)` — where `meta_lines`, `names`, `units` reproduce EXACTLY what `makeheader` wrote (read `makeheader` lines ~530-554 and build the same `*`-prefixed lines and the same column-name + unit rows from the active `vars`). In `case(2)`, replace the hand-built `line`/`write(iuncsv,…)` with: flatten the active `vars%value(1:Nnodes,j)` (for `j` with `iyes==1`) into a 1-D `real64` array in the SAME emission order, then `call scalar_w%row(values, leading=trim(tc_date_or_datetime))`. In `case(3)` replace `close(iuncsv)` with `scalar_w%close()`. Keep the headless guard (only open/write when `.not. headless`).

- [ ] **Step 2: Route `csv_out_tz` through a second writer**

Add module-level `type(csv_writer_t), save :: profile_w` and do the equivalent replacement in `csv_out_tz` (open/header in its init task, `row` per depth line in its write task, `close` in its close task).

- [ ] **Step 3: Per-year flush**

At the end of `csv_output_step` (in `module csv_output`), add:
```fortran
if (state%timecontrol%flYearStart) then
   call scalar_w%flush()
end if
```
(and `call profile_w%flush()` if the tz writer is module-accessible from `csv_output`; if `profile_w` lives in `csv_output_tz`, add a small public `csv_out_tz_flush()` there and call it). `flush` is a no-op when the unit was never opened (headless), so it is safe.

> `flYearStart` confirmed to exist: `state%timecontrol%flYearStart` (`src/state/timecontrol_state.f90:147`).

- [ ] **Step 4: Build + regression + byte-diff**

Run `pixi run -e test check-fast` → 4/4. Re-run the baseline case and `diff result_output.csv /tmp/baseline_output.csv`. Aim for ZERO diff; if there is only cosmetic whitespace drift but column names + numeric values match within tol and check-fast is 4/4, that is acceptable — DOCUMENT any diff in the commit body. If a column NAME changed or a value moved beyond tol, it is a bug (likely the flatten order or a header-name mismatch).

- [ ] **Step 5: Commit**

```bash
git add src/io/csv_output.f90
git commit -m "refactor(io): IO-OUT/C — route csv_output writes through csv_writer; per-year flush"
```

---

## Phase D — column-parity audit (`.crp`, `.snw`)

Prove the two result CSVs already carry every `.crp`/`.snw` quantity before Phase E deletes those writers.

### Task D1: Enumerate legacy columns; add any missing registry entries

**Files:**
- Modify: `src/io/output_registry.f90` (+ sampling case in `csv_output.f90`) only if a variable is missing

- [ ] **Step 1: List the legacy columns**

From `swapoutput.f90`: read `OutCropFixed` (550-620), `OutWofost` (623-732), `OutGrass` (735-811), and `SnowOutput`/`build_snow_output_row` (1373-1513). Write down every quantity each writes.

- [ ] **Step 2: Map to registry**

For each legacy column, find the registry entry that carries the same `state` quantity. Snow maps to `SNOW`/`SSNOW`/`SUBLIM` and the snow-row fields (`snrai/gsnow/ssnow/melt/subl`); crop maps to DVS/LAI/CWDM/CWSO/HEIGHT/RD/etc. List any quantity with **no** registry entry.

- [ ] **Step 3: Add the missing entries (only if Step 2 found gaps)**

For each gap, add an `out_var_t(...)` to `REGISTRY` and a sampling line in `csv_output`'s `set_values` (registry index → `state%…`). Build + `check-fast`; the new columns must not perturb existing-column values.

- [ ] **Step 4: Commit (only if entries were added)**

```bash
git add src/io/output_registry.f90 src/io/csv_output.f90
git commit -m "feat(io): IO-OUT/D — add crop/snow registry entries for parity"
```

### Task D2: Parity script as deletion evidence

**Files:**
- Create: `scripts/output_parity_check.py`

- [ ] **Step 1: Write the script**

Create `scripts/output_parity_check.py`:

```python
"""Prove result_output.csv reproduces legacy .crp/.snw columns within tol.

Run a case that still emits the legacy file, parse both, and assert per-day
agreement on the mapped columns. Exit non-zero on any mismatch. This is the
evidence gate before IO-OUT/E deletes the .crp/.snw writers.
"""
import csv, sys
from pathlib import Path

TOL = 1e-2

# (legacy column in .snw/.crp) -> (column in result_output.csv)
SNOW_MAP = {"ssnow": "SSNOW", "subl": "SUBLIM"}   # extend per Task D1 Step 2

def read_csv_skip_star(path):
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("*"):
                continue
            rows.append(line)
    return list(csv.DictReader(rows))

def check(legacy_path, result_path, colmap):
    legacy = read_csv_skip_star(legacy_path)
    result = read_csv_skip_star(result_path)
    n = min(len(legacy), len(result))
    bad = 0
    for i in range(n):
        for lk, rk in colmap.items():
            lv, rv = legacy[i].get(lk), result[i].get(rk)
            if lv in (None, "") or rv in (None, ""):
                continue
            if abs(float(lv) - float(rv)) > TOL:
                print(f"row {i}: {lk}={lv} vs {rk}={rv}")
                bad += 1
    return bad

if __name__ == "__main__":
    workdir = Path(sys.argv[1])   # a run dir containing both files
    bad = check(workdir / "result.snw", workdir / "result_output.csv", SNOW_MAP)
    if bad:
        print(f"PARITY FAIL: {bad} mismatches")
        sys.exit(1)
    print("PARITY OK")
```

- [ ] **Step 2: Run it against a fresh run**

Build, run a snow-active case so both `result.snw` and `result_output.csv` exist in one workdir, then:
Run: `python scripts/output_parity_check.py <workdir>`
Expected: `PARITY OK`. (Tune `SNOW_MAP`/add a `CROP_MAP` per the Step-1/2 findings until parity holds.)

- [ ] **Step 3: Commit**

```bash
git add scripts/output_parity_check.py
git commit -m "test(io): IO-OUT/D — parity script proving CSV reproduces .crp/.snw"
```

---

## Phase E — retire the graveyard

Removes `swapoutput.f90`, drops the 4 inert C-API streams, re-homes the live BMI rows, fixes `water_balance_row`, rewires the main loop. **State-schema change → `rm -rf builddir`. Gate is `check-full`, not `check-fast`.**

### Task E1: Re-home live BMI rows + fix `water_balance_row` in `csv_output`

**Files:**
- Modify: `src/io/csv_output.f90`

- [ ] **Step 1: Fill `water_balance_row` from `compute_aggregates`**

In `csv_output_step`, after computing the aggregates, populate `state%water_balance_row(:)` with the same 18 fields the dead `outinc`/`build_water_balance_row` defined (date, day, dcum, rain, snow, irrig, interc, runon, runoff, tpot, tact, epot, eact, drainage, qbottom, gwl, dstorage, baldev). Allocate it in `csv_output_init`, deallocate in `csv_output_finalize`. (This fixes the `swap_balance` stream that returns zeros today.)

- [ ] **Step 2: Move `build_snow_output_row` into `csv_output_step`**

Copy the 7-field assignment from `swapoutput.f90:1495-1513` into `csv_output_step`, gated on `state%timecontrol%flSnow`. Allocate `state%atmosphere%snow_output_row` in `csv_output_init` (gated on flSnow), deallocate in `csv_output_finalize`.

- [ ] **Step 3: Build + regression**

Run: `pixi run -e test check-fast`
Expected: 4/4 unchanged (CSV output untouched; only buffers now filled here).

- [ ] **Step 4: Commit**

```bash
git add src/io/csv_output.f90
git commit -m "refactor(io): IO-OUT/E — fill water_balance + snow BMI rows in csv_output"
```

### Task E2: Rewire the main loop to call `csv_output_*` directly

**Files:**
- Modify: `src/core/swap_mod.f90` (≈ lines 421, 615-629, 662-664)

- [ ] **Step 1: Replace the output dispatch calls**

- Init block (≈421-429): replace `call SwapOutput(1,…)` + `call SoilWaterOutput(1,…)` + the `if (flSnow) call SnowOutput(1,…)` etc. with a single `call csv_output_init(state, config)`.
- Per-step block (≈615-629): replace `SwapOutput(2)`/`SoilWaterOutput(2)`/`SnowOutput(2)`/… (incl. the `flOutputShort` path) with `call csv_output_step(state)`.
- Close block (≈662-673): replace `SwapOutput(3)`/`SoilWaterOutput(4)`/… with `call csv_output_finalize(state)`.
- Add `use csv_output, only: csv_output_init, csv_output_step, csv_output_finalize`; remove the `*Output` externals.

- [ ] **Step 2: Build**

Run: `pixi run -e test build-linux`
Expected: compiles; `swapoutput.f90` still present (deleted in E4) but no longer called — link still succeeds.

- [ ] **Step 3: Regression**

Run: `pixi run -e test check-fast`
Expected: 4/4 unchanged.

- [ ] **Step 4: Commit**

```bash
git add src/core/swap_mod.f90
git commit -m "refactor(core): IO-OUT/E — main loop calls csv_output_init/step/finalize"
```

### Task E3: Drop `.crp` writer call sites

**Files:**
- Modify: `src/crop/cropgrowth_helpers.f90` (≈100-141)

- [ ] **Step 1: Remove the calls**

Delete the `call OutCropFixed/OutWofost/OutGrass(1|2, state)` lines at `cropgrowth_helpers.f90:100-141`. Crop data remains in `result_output.csv` (verified Phase D) and the `crop` C-API stream (`crop_output_row`, built elsewhere in this file — leave that).

- [ ] **Step 2: Build + regression**

Run: `pixi run -e test check-fast`
Expected: compiles (the `OutCrop*` subroutines still exist until E4); 4/4 unchanged. `.crp` files simply stop being produced.

- [ ] **Step 3: Commit**

```bash
git add src/crop/cropgrowth_helpers.f90
git commit -m "refactor(crop): IO-OUT/E — stop writing legacy .crp files"
```

### Task E4: Delete `swapoutput.f90`; drop 4 inert streams + state fields

**Files:**
- Delete: `src/io/swapoutput.f90`
- Modify: `src/core/swap_capi_mod.f90` (remove `temperature`/`solute`/`surfacewater`/`agetracer` cases, ≈240-260)
- Modify: `src/state/heat_state.f90` (remove `output_row`/`output_columns`/`output_n_cols`)
- Modify: `src/state/solute_state.f90` (remove `output_row`/`agetracer_row` + their columns/n_cols)
- Modify: `src/state/surfacewater_state.f90` (remove `output_row`/columns/n_cols)
- Modify: `meson.build` (drop line 146 `'src/io/swapoutput.f90',`)

- [ ] **Step 1: Pre-flight — confirm no consumer of the dropped streams**

Run:
```bash
grep -rniE "'temperature'|'solute'|'surfacewater'|'agetracer'|get_output_row" tests/bmi tests/cffi-demo
```
Expected: no test asserts success on these four streams (they return zeros today). If one does, stop and report — do not delete that stream.

- [ ] **Step 2: Delete the file + meson entry**

```bash
git rm src/io/swapoutput.f90
```
Remove `'src/io/swapoutput.f90',` (line 146) from `meson.build`.

- [ ] **Step 3: Drop the four stream cases**

In `swap_capi_mod.f90`, delete the `case ('temperature')`, `case ('solute')`, `case ('agetracer')`, `case ('surfacewater')` arms (≈240-274). Keep `swap_balance`, `soilwater`, `snow`, `crop`, `tillage`.

- [ ] **Step 4: Remove the now-unused state fields**

In `heat_state.f90`, `solute_state.f90`, `surfacewater_state.f90`, delete the `output_row`/`agetracer_row`/`*_columns`/`*_n_cols` declarations and any `init_*_buffer`/`cleanup_*_buffer` helpers that lived in the deleted `swapoutput.f90` and referenced them. Grep to confirm no remaining references:
```bash
grep -rniE "heat%output_row|solute%output_row|solute%agetracer_row|surfacewater%output_row" src/
```
Expected: no hits.

- [ ] **Step 5: Clean rebuild (state schema changed)**

Run:
```bash
rm -rf builddir
pixi run -e test build-linux
```
Expected: compiles clean.

- [ ] **Step 6: Full gate**

Run: `pixi run -e test check-full`
Expected: 5/5 `regression ok` + pFUnit green.

- [ ] **Step 7: C-API smoke check**

Run: `pixi run -e test test-cffi-demo` and `pixi run -e test test-bmi`
Expected: pass. Optionally verify `swap_get_output_row('swap_balance')` now returns a non-zero, correctly-sized row, and the four dropped streams return `ierr=1`.

- [ ] **Step 8: Commit**

```bash
# git rm already staged the deletion; scoped add for the rest (never `git add -A`)
git add meson.build src/core/swap_capi_mod.f90 \
        src/state/heat_state.f90 src/state/solute_state.f90 src/state/surfacewater_state.f90
git commit -m "refactor(io): IO-OUT/E — delete swapoutput.f90; drop 4 inert C-API streams"
```

---

## Self-review checklist (run before handing off to execution)

- **Spec coverage:** A→csv_writer; B→registry; C→csv_writer routing + rename + aggregates; D→parity; E→deletion + stream drop + rewire + water_balance fix + snow re-home. Out-of-scope IO-IN items are explicitly deferred. ✔
- **Placeholder scan:** the only "fill from source" instructions (registry catalogue in B1, the `what_form`/`fill_values` verbatim ports) point at exact source line ranges with the invariant to preserve — they are guided extractions of existing tested code, not TBDs. ✔
- **Type consistency:** `csv_writer_t` (open/meta/header/row/close) used identically in A and C; `out_var_t`/`OUT_*`/`resolve_inlist`/`var_*` consistent B→D; `csv_output_init/step/finalize` consistent C→E; `water_balance_row`/`snow_output_row` field names match the state decls. ✔
