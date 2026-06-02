# In-Memory Results Record Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Accumulate SWAP's scalar output (the user's `inlist` columns) into a per-instance typed record in memory, so Python can pull a full-precision results array directly — without parsing the CSV file.

**Architecture:** Mirror the input side ("config-as-data"). The output step already builds, each output interval, the selected + node-expanded `vals(1:ncount)` row right before writing it to CSV (`csv_out_write`). We add a per-instance `results_state_t` on `swap_state_t` and `add_row(time, vals)` it at that same point — so the sieve (inlist selection) is already applied and compute logic is untouched. `add_row` buffers into a chunk and `flush`es every N rows into a growing in-memory store; the chunk/flush seam is the future swap-in point for a binary disk buffer. The CSV file write is left exactly as-is this pass (byte-identical); only the *build of `vals`* moves outside the headless guard so the record also populates in headless (pyswap) mode.

**Tech Stack:** Modern Fortran (gfortran, `iso_fortran_env`), pFUnit unit tests, meson build, Python `ctypes` for the C-ABI demo.

**Scope (this pass):** Scalar output only (the `result_output.csv` family). The depth–time profile output (`csv_out_tz_*`) is memory-heavy and stays on its current streaming path. The CSV file remains written by the existing code; making CSV serialize *from* the record (and the binary disk sink) are explicit follow-up passes.

**Non-negotiable:** Byte-identical regression vs `swap420gf` (`pixi run -e test check-full`). No physics/format/ordering changes to the on-disk CSV. No new module globals — the record lives on `state` (per-instance).

**Build/test commands:**
- Build: `pixi run build-linux`
- Unit tests: `pixi run -e test test-pfunit`  (look for the `OK (N tests)` count)
- Fast regression gate: `pixi run -e test check-fast`
- Full regression: `pixi run -e test check-full`
- Clean rebuild (MANDATORY after any `src/state/*.f90` schema change): `pixi run clean` then build.

---

## File Structure

- **Create** `src/state/results_state.f90` — defines `results_state_t` and its type-bound `init` / `add_row` / `flush` / `finalize`. In-memory sink, chunked (flush-every-N), scalar columns. No I/O dependency.
- **Modify** `src/state/swap_state.f90` — add `type(results_state_t) :: results` to `swap_state_t`. (schema change → clean rebuild)
- **Modify** `src/io/csv_output.f90` — `csv_out_header` initialises the record from the selected columns; `csv_out_write` builds `vals` and `add_row`s it *always* (outside the headless guard); `csv_out_close` calls `finalize`.
- **Modify** `src/bindings/swap_capi_mod.f90` — add C-ABI accessors `swap_results_shape`, `swap_results_columns`, `swap_view_results`, `swap_view_results_times`.
- **Modify** `prototype/pyswap_inmemory_demo.py` — pull the results array (full precision) and a column list; verify against the CSV.
- **Create** `tests/unit/state/test_results_state.pf` — unit tests for the record.
- **Modify** `tests/unit/testSuites.inc`, `tests/unit/meson.build`, `meson.build` — register the new source + test.

---

## Task 1: `results_state_t` skeleton — `init` + `add_row` (grow-in-memory)

**Files:**
- Create: `src/state/results_state.f90`
- Test: `tests/unit/state/test_results_state.pf`
- Modify: `meson.build` (add source), `tests/unit/meson.build` (add test + source), `tests/unit/testSuites.inc` (register suite)

- [ ] **Step 1: Write the failing test**

Create `tests/unit/state/test_results_state.pf`:

```fortran
! Tests for results_state_t — the per-instance in-memory scalar-output record.

@test
subroutine test_results_add_row_stores_values()
   use funit
   use iso_fortran_env, only: real64
   use results_state_mod, only: results_state_t

   type(results_state_t) :: r

   call r%init([character(len=24) :: 'RAIN', 'TACT'], &
               [character(len=12) :: 'cm', 'cm'], flush_every=2)
   call r%add_row(100.0_real64, [1.0_real64, 2.0_real64])
   call r%add_row(101.0_real64, [3.0_real64, 4.0_real64])
   call r%finalize()

   @assertEqual(2, r%ncols)
   @assertEqual(2, r%nrows)
   @assertEqual(100.0_real64, r%times(1), 1.0e-12_real64)
   @assertEqual(4.0_real64,   r%values(2, 2), 1.0e-12_real64)
   @assertEqual('TACT', trim(r%col_names(2)))
end subroutine test_results_add_row_stores_values
```

- [ ] **Step 2: Add the new source + test to the build, register the suite**

In `meson.build`, after the `'src/state/swap_state.f90',` line (or near the other `src/state/*` entries), add:
```
    'src/state/results_state.f90',
```
In `tests/unit/meson.build`, in `pfunit_extra_sources` near the other `'../../src/state/...'` entries add:
```
        '../../src/state/results_state.f90',
```
and in the `pf_files` list near the other `'state/...'` entries add:
```
        'state/test_results_state.pf',
```
In `tests/unit/testSuites.inc`, near the other state suites add:
```
ADD_TEST_SUITE(test_results_state_suite)
```

- [ ] **Step 3: Create a stub module so it compiles but the test fails**

Create `src/state/results_state.f90`:

```fortran
!> @file results_state.f90
!! Per-instance in-memory record of SWAP's scalar output. Holds only the
!! user's selected (inlist) columns. Rows are appended via add_row and
!! committed to a growing in-memory store by flush every `flush_every` rows;
!! the chunk/flush seam is the swap-in point for a future binary disk sink.
!! Scalar output only (depth-time profile output stays on its streaming path).
module results_state_mod
   use iso_fortran_env, only: real64
   implicit none
   private
   public :: results_state_t

   type :: results_state_t
      integer :: ncols       = 0
      integer :: nrows       = 0          !! committed rows
      integer :: flush_every = 512        !! rows per chunk before a flush
      character(len=24), allocatable :: col_names(:)
      character(len=12), allocatable :: col_units(:)
      real(real64),      allocatable :: times(:)        !! (nrows)
      real(real64),      allocatable :: values(:,:)     !! (nrows, ncols)
      ! chunk buffer (filled by add_row, drained by flush)
      integer :: chunk_n = 0
      real(real64), allocatable :: chunk_t(:)           !! (flush_every)
      real(real64), allocatable :: chunk_v(:,:)         !! (flush_every, ncols)
   contains
      procedure :: init     => results_state_init
      procedure :: add_row  => results_state_add_row
      procedure :: flush    => results_state_flush
      procedure :: finalize => results_state_finalize
   end type results_state_t

contains

   subroutine results_state_init(self, names, units, flush_every)
      class(results_state_t), intent(inout) :: self
      character(len=*),       intent(in)    :: names(:)
      character(len=*),       intent(in)    :: units(:)
      integer, optional,      intent(in)    :: flush_every
      ! STUB (RED)
   end subroutine results_state_init

   subroutine results_state_add_row(self, time, vals)
      class(results_state_t), intent(inout) :: self
      real(real64),           intent(in)    :: time
      real(real64),           intent(in)    :: vals(:)
      ! STUB (RED)
   end subroutine results_state_add_row

   subroutine results_state_flush(self)
      class(results_state_t), intent(inout) :: self
      ! STUB (RED)
   end subroutine results_state_flush

   subroutine results_state_finalize(self)
      class(results_state_t), intent(inout) :: self
      ! STUB (RED)
   end subroutine results_state_finalize

end module results_state_mod
```

- [ ] **Step 4: Run the test to verify it fails**

Run: `pixi run -e test test-pfunit 2>&1 | grep -iE "test_results|Fail:"`
Expected: `test_results_state_suite.test_results_add_row_stores_values` listed under failures (ncols/nrows are 0).

- [ ] **Step 5: Implement `init`, `add_row`, `flush`, `finalize`**

Replace the four stub bodies in `src/state/results_state.f90`:

```fortran
   subroutine results_state_init(self, names, units, flush_every)
      class(results_state_t), intent(inout) :: self
      character(len=*),       intent(in)    :: names(:)
      character(len=*),       intent(in)    :: units(:)
      integer, optional,      intent(in)    :: flush_every

      self%ncols = size(names)
      self%nrows = 0
      self%chunk_n = 0
      if (present(flush_every)) self%flush_every = max(1, flush_every)
      self%col_names = names
      self%col_units = units
      allocate(self%times(0))
      allocate(self%values(0, self%ncols))
      allocate(self%chunk_t(self%flush_every))
      allocate(self%chunk_v(self%flush_every, self%ncols))
   end subroutine results_state_init

   subroutine results_state_add_row(self, time, vals)
      class(results_state_t), intent(inout) :: self
      real(real64),           intent(in)    :: time
      real(real64),           intent(in)    :: vals(:)
      integer :: m
      m = min(size(vals), self%ncols)
      self%chunk_n = self%chunk_n + 1
      self%chunk_t(self%chunk_n)        = time
      self%chunk_v(self%chunk_n, 1:m)   = vals(1:m)
      if (m < self%ncols) self%chunk_v(self%chunk_n, m+1:self%ncols) = 0.0_real64
      if (self%chunk_n == self%flush_every) call self%flush()
   end subroutine results_state_add_row

   !> Drain the chunk buffer into the committed store. For the in-memory sink
   !! this grows `times`/`values`; a future disk sink would write the chunk to
   !! disk here and leave the in-memory store empty.
   subroutine results_state_flush(self)
      class(results_state_t), intent(inout) :: self
      real(real64), allocatable :: t2(:), v2(:,:)
      integer :: old, add
      if (self%chunk_n == 0) return
      old = self%nrows
      add = self%chunk_n
      allocate(t2(old + add))
      allocate(v2(old + add, self%ncols))
      if (old > 0) then
         t2(1:old)      = self%times
         v2(1:old, :)   = self%values
      end if
      t2(old+1:old+add)    = self%chunk_t(1:add)
      v2(old+1:old+add, :) = self%chunk_v(1:add, :)
      call move_alloc(t2, self%times)
      call move_alloc(v2, self%values)
      self%nrows   = old + add
      self%chunk_n = 0
   end subroutine results_state_flush

   subroutine results_state_finalize(self)
      class(results_state_t), intent(inout) :: self
      call self%flush()
   end subroutine results_state_finalize
```

- [ ] **Step 6: Run the test to verify it passes**

Run: `pixi run -e test test-pfunit 2>&1 | tail -3`
Expected: `OK (N tests)` with N one higher than before; no failures.

- [ ] **Step 7: Commit**

```bash
git add src/state/results_state.f90 tests/unit/state/test_results_state.pf meson.build tests/unit/meson.build tests/unit/testSuites.inc
git commit -m "feat(state): results_state_t — in-memory scalar-output record (add_row/flush)"
```

---

## Task 2: flush-every-N actually chunks (multiple flushes, partial tail)

**Files:**
- Test: `tests/unit/state/test_results_state.pf`

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/state/test_results_state.pf`:

```fortran
@test
subroutine test_results_flush_every_n_and_partial_tail()
   use funit
   use iso_fortran_env, only: real64
   use results_state_mod, only: results_state_t

   type(results_state_t) :: r
   integer :: k

   call r%init([character(len=24) :: 'A'], [character(len=12) :: '-'], flush_every=2)
   ! 5 rows with flush_every=2 => flushes after row 2 and 4; finalize flushes row 5.
   do k = 1, 5
      call r%add_row(real(k, real64), [real(10*k, real64)])
   end do
   ! before finalize: 4 rows committed, 1 buffered
   @assertEqual(4, r%nrows)
   call r%finalize()
   @assertEqual(5, r%nrows)
   @assertEqual(50.0_real64, r%values(5, 1), 1.0e-12_real64)
   @assertEqual(5.0_real64,  r%times(5),     1.0e-12_real64)
end subroutine test_results_flush_every_n_and_partial_tail
```

- [ ] **Step 2: Run the test to verify it fails or passes**

Run: `pixi run -e test test-pfunit 2>&1 | grep -iE "flush_every|Fail:|OK \("`
Expected: PASS already (the Task 1 implementation handles chunking). If it passes, this test simply locks the behavior — no code change needed. If it fails, fix `add_row`/`flush` until green. (Either way: never edit the test to match a bug.)

- [ ] **Step 3: Commit**

```bash
git add tests/unit/state/test_results_state.pf
git commit -m "test(state): lock flush-every-N chunking + partial-tail behavior"
```

---

## Task 3: Add `results` to `swap_state_t`

**Files:**
- Modify: `src/state/swap_state.f90`

- [ ] **Step 1: Add the field**

In `src/state/swap_state.f90`, add the use and the component. Near the other `use *_state_mod` lines add:
```fortran
   use results_state_mod, only: results_state_t
```
In the `type :: swap_state_t` definition, alongside the other subsystem records add:
```fortran
      type(results_state_t) :: results
```

- [ ] **Step 2: Clean rebuild (schema change) and confirm the suite still passes**

Run:
```bash
pixi run clean && pixi run -e test test-pfunit 2>&1 | tail -3
```
Expected: `OK (N tests)`, no failures (the field is inert so far).

- [ ] **Step 3: Confirm regression still byte-identical**

Run: `pixi run -e test check-fast 2>&1 | grep -iE "Results"`
Expected: `Results: 4 passed, 0 failed`.

- [ ] **Step 4: Commit**

```bash
git add src/state/swap_state.f90
git commit -m "feat(state): add per-instance results record to swap_state_t"
```

---

## Task 4: Wire `csv_output` to populate the record (byte-identical)

**Files:**
- Modify: `src/io/csv_output.f90:657-707` (`csv_out_header`), `:709-755` (`csv_out_write`), `:757-769` (`csv_out_close`)
- Test: `tests/unit/io/toml/test_config_source.pf` is unrelated — add the integration assertion in a new test file `tests/unit/io/test_results_wiring.pf`

Context for the implementer: `csv_out_write` already builds `vals(1:ncount)` (the selected, node-expanded row) but **only inside `if (.not. state%timecontrol%headless)`**. `makeheader` (line 892) builds the matching column names from `vars%head(k,j)` over the same `iyes`/`Nnodes` loop. The record must mirror that loop for names and build `vals` outside the headless guard.

- [ ] **Step 1: Write the failing integration test**

Create `tests/unit/io/test_results_wiring.pf`:

```fortran
! Integration: after a headless run, state%results holds one row per output
! step with the selected (inlist) columns — i.e. the record is populated even
! when the CSV file is suppressed.

@test
subroutine test_results_populated_headless()
   use funit
   use iso_fortran_env,   only: real64
   use swap_state_mod,    only: swap_state_t
   use swap_config_mod,   only: swap_config_t
   use swap_mod,          only: swap_init, swap_run_step, swap_close

   type(swap_state_t)  :: state
   type(swap_config_t) :: config
   character(*), parameter :: cfg = 'tests/swap-cases/toml/1.hupselbrook/swap.toml'
   integer :: guard

   state%timecontrol%headless = .true.
   call swap_init(cfg, state, config)
   ! headless is reset by init in some builds; force it for the run.
   state%timecontrol%headless = .true.

   guard = 0
   do while (.not. state%timecontrol%flRunEnd .and. guard < 1000000)
      call swap_run_step(state, config)
      guard = guard + 1
   end do
   call swap_close(state, config)

   @assertTrue(state%results%ncols > 0)
   @assertTrue(state%results%nrows > 0)
   @assertEqual(state%results%nrows, size(state%results%times))
end subroutine test_results_populated_headless
```

Register it: in `tests/unit/meson.build` `pf_files` add `'io/test_results_wiring.pf'`; in `tests/unit/testSuites.inc` add `ADD_TEST_SUITE(test_results_wiring_suite)`.

NOTE: confirm `swap_init`/`swap_run_step`/`swap_close` are the correct public names in `src/driver/swap_mod.f90` (they are at the time of writing) and that `state%timecontrol%headless` / `%flRunEnd` exist; adjust the harness if the run-loop entry points differ.

- [ ] **Step 2: Run the test to verify it fails**

Run: `pixi run -e test test-pfunit 2>&1 | grep -iE "results_populated|Fail:"`
Expected: FAIL — `ncols`/`nrows` are 0 (record never initialised/filled).

- [ ] **Step 3: Initialise the record in `csv_out_header`**

In `src/io/csv_output.f90`, inside `csv_out_header`, AFTER `call det_which_vars(state)` (line 688) and the headless file-open block, BEFORE the "store initial values" block (line 700), insert:

```fortran
      ! Build the selected column names/units in the SAME order csv_out_write
      ! flattens vals (iyes, then per node), and initialise the per-instance
      ! in-memory record. Always — independent of headless file output.
      block
         integer :: jj, kk, nc
         character(len=24) :: rnames(Mnodes*M)
         character(len=12) :: runits(Mnodes*M)
         nc = 0
         do jj = 1, M
            if (vars%iyes(jj) == 1) then
               do kk = 1, vars%Nnodes(jj)
                  nc = nc + 1
                  rnames(nc) = vars%head(kk, jj)
                  runits(nc) = vars%unit(jj)
               end do
            end if
         end do
         call state%results%init(rnames(1:nc), runits(1:nc))
      end block
```

- [ ] **Step 4: Append rows in `csv_out_write` (outside the headless guard)**

In `src/io/csv_output.f90`, replace the body of `csv_out_write` from the comment block through `end associate` (lines 726-753) with:

```fortran
      ! IO-OUT/C2b + results record: flatten the active values (selected +
      ! node-expanded, same emission order as the CSV) ONCE, append to the
      ! per-instance record (always), then write the CSV row (file only when
      ! not headless). Building vals outside the headless guard is what lets
      ! the in-memory record populate in pyswap's headless mode.
      associate (time => state%timecontrol)
         ncount = 0
         do j = 1, M
            if (vars%iyes(j) == 1) then
               do n = 1, vars%Nnodes(j)
                  ncount = ncount + 1
                  vals(ncount) = vars%value(n,j)
               end do
            end if
         end do

         call state%results%add_row(time%t1900, vals(1:ncount))

         if (.not. time%headless) then
            if (.not. time%flprintshort) then
               call scalar_w%row(vals(1:ncount), leading=trim(time%date))
            else
               call dtdpst ('year-month-day hour:minute:seconds', time%t1900, datexti)
               call scalar_w%row(vals(1:ncount), leading=trim(datexti))
            end if
         end if
      end associate
```

- [ ] **Step 5: Finalize the record in `csv_out_close`**

In `src/io/csv_output.f90`, in `csv_out_close`, after the headless-guarded `scalar_w%close()` block (line 767), add:
```fortran
      call state%results%finalize()
```

- [ ] **Step 6: Clean rebuild and run the integration test**

Run:
```bash
pixi run clean && pixi run -e test test-pfunit 2>&1 | grep -iE "results_populated|Fail:|OK \("
```
Expected: PASS; `OK (N tests)`.

- [ ] **Step 7: Verify byte-identical (the critical gate)**

Run: `pixi run -e test check-full 2>&1 | grep -iE "Results:"`
Expected: `Results: 11 passed, 0 failed, 7 known-divergence (xfail)` (the 7 xfail are pre-existing; no NEW failures). If any non-xfail case fails, STOP — the vals-build move changed the on-disk output; diff and fix before continuing.

- [ ] **Step 8: Commit**

```bash
git add src/io/csv_output.f90 tests/unit/io/test_results_wiring.pf tests/unit/meson.build tests/unit/testSuites.inc
git commit -m "feat(io): populate per-instance results record from csv_output (byte-identical)"
```

---

## Task 5: C-ABI accessors for the results record

**Files:**
- Modify: `src/bindings/swap_capi_mod.f90`

The record lives on `capi_state%results`. `capi_state` is already `save, target`, so `c_loc` of its array components is valid (same pattern as `swap_view_array`).

- [ ] **Step 1: Add the accessors**

In `src/bindings/swap_capi_mod.f90`, after `swap_get_water_balance` (line 207), add:

```fortran
   !----------------------------------------------------------------------
   ! In-memory results record accessors
   !----------------------------------------------------------------------

   function swap_results_shape(nrows, ncols) result(ierr) bind(C, name='swap_results_shape')
      integer(c_int), intent(out) :: nrows, ncols
      integer(c_int)              :: ierr
      nrows = capi_state%results%nrows
      ncols = capi_state%results%ncols
      ierr  = 0
   end function swap_results_shape

   !> Zero-copy view of the (nrows x ncols) values array. Fortran column-major:
   !! element (i,j) is at flat index (j-1)*nrows + (i-1).
   function swap_view_results(ptr, nrows, ncols) result(ierr) bind(C, name='swap_view_results')
      type(c_ptr),    intent(out) :: ptr
      integer(c_int), intent(out) :: nrows, ncols
      integer(c_int)              :: ierr
      nrows = capi_state%results%nrows
      ncols = capi_state%results%ncols
      if (allocated(capi_state%results%values) .and. nrows > 0 .and. ncols > 0) then
         ptr  = c_loc(capi_state%results%values(1,1))
         ierr = 0
      else
         ptr  = c_null_ptr
         ierr = 1
      end if
   end function swap_view_results

   !> Zero-copy view of the time axis (nrows doubles).
   function swap_view_results_times(ptr, nrows) result(ierr) bind(C, name='swap_view_results_times')
      type(c_ptr),    intent(out) :: ptr
      integer(c_int), intent(out) :: nrows
      integer(c_int)              :: ierr
      nrows = capi_state%results%nrows
      if (allocated(capi_state%results%times) .and. nrows > 0) then
         ptr  = c_loc(capi_state%results%times(1))
         ierr = 0
      else
         ptr  = c_null_ptr
         ierr = 1
      end if
   end function swap_view_results_times

   !> NUL-delimited column names packed into buf (truncated to n bytes).
   function swap_results_columns(buf, n) result(ierr) bind(C, name='swap_results_columns')
      character(kind=c_char), intent(out) :: buf(*)
      integer(c_int),  value, intent(in)  :: n
      integer(c_int)                      :: ierr
      integer :: i, k, c
      k = 0
      do c = 1, capi_state%results%ncols
         do i = 1, len_trim(capi_state%results%col_names(c))
            k = k + 1
            if (k > n) then; buf(min(k,n)) = c_null_char; ierr = 1; return; end if
            buf(k) = capi_state%results%col_names(c)(i:i)
         end do
         k = k + 1
         if (k > n) then; buf(min(k,n)) = c_null_char; ierr = 1; return; end if
         buf(k) = c_null_char
      end do
      ierr = 0
   end function swap_results_columns
```

- [ ] **Step 2: Build the library**

Run: `pixi run build-linux 2>&1 | grep -iE "error|FAILED" | head`
Expected: no output (clean build); `builddir/libswap_bmi.so` present.

- [ ] **Step 3: Commit**

```bash
git add src/bindings/swap_capi_mod.f90
git commit -m "feat(bindings): zero-copy C-ABI accessors for the in-memory results record"
```

---

## Task 6: Python demo pulls a full-precision results array

**Files:**
- Modify: `prototype/pyswap_inmemory_demo.py`

- [ ] **Step 1: Add a results pull + a self-check in the worker**

In `prototype/pyswap_inmemory_demo.py`, extend `load_lib()` argtypes with:

```python
    lib.swap_results_shape.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
    lib.swap_view_results.argtypes = [ctypes.POINTER(ctypes.c_void_p),
                                      ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
    lib.swap_results_columns.argtypes = [ctypes.c_char_p, ctypes.c_int]
```

After `run_to_end(lib)` returns in `worker`, before printing, add a results read:

```python
    import numpy as np
    nr, nc = ctypes.c_int(0), ctypes.c_int(0)
    lib.swap_results_shape(ctypes.byref(nr), ctypes.byref(nc))
    ptr = ctypes.c_void_p()
    rc = lib.swap_view_results(ctypes.byref(ptr), ctypes.byref(nr), ctypes.byref(nc))
    cols = b"\x00" * 4096
    lib.swap_results_columns(cols, len(cols))
    names = [s.decode() for s in cols.split(b"\x00") if s][:nc.value]
    arr = None
    if rc == 0 and nr.value > 0 and nc.value > 0:
        flat = (ctypes.c_double * (nr.value * nc.value)).from_address(ptr.value)
        # Fortran column-major -> shape (ncols, nrows) then transpose
        arr = np.ctypeslib.as_array(flat).reshape((nc.value, nr.value)).T.copy()
    wb["_results_nrows"] = nr.value
    wb["_results_ncols"] = nc.value
    wb["_results_cols"] = names
    if arr is not None:
        wb["_results_first_row"] = list(arr[0])
        wb["_results_last_row"] = list(arr[-1])
```

(`wb` is the dict returned by `run_to_end`; the underscore keys piggyback on the existing JSON round-trip.)

- [ ] **Step 2: Print the results summary in `main` and verify mem == disk**

In `main`, after the water-balance comparison loop, add:

```python
    print("\nResults record (in-memory):")
    print(f"  shape = {mem['_results_nrows']} rows x {mem['_results_ncols']} cols")
    print(f"  columns = {mem['_results_cols']}")
    print(f"  first row = {mem.get('_results_first_row')}")
    # the record must be identical mem vs disk too
    if mem.get("_results_first_row") and disk.get("_results_first_row"):
        rec_ok = all(abs(a - b) <= 1e-10
                     for a, b in zip(mem["_results_first_row"], disk["_results_first_row"])) \
                 and all(abs(a - b) <= 1e-10
                         for a, b in zip(mem["_results_last_row"], disk["_results_last_row"]))
        print("  record mem==disk:", "PASS" if rec_ok else "FAIL")
        ok = ok and rec_ok
```

- [ ] **Step 3: Run the demo**

Run: `python3 prototype/pyswap_inmemory_demo.py 2>&1 | tail -20`
Expected: the water-balance table, then a `Results record (in-memory)` block showing a non-zero shape, the inlist column names, and `record mem==disk: PASS`, and final `RESULT: PASS`.

- [ ] **Step 4: Commit**

```bash
git add prototype/pyswap_inmemory_demo.py
git commit -m "feat(prototype): pull full-precision in-memory results array into numpy from Python"
```

---

## Follow-ups (NOT in this plan — explicit next passes)

1. **CSV-serializes-from-record.** Add a CSV-file sink so `flush`/`finalize` writes the file *from* the record (one source of truth), retiring the parallel `scalar_w%row` path. Gate hard on `check-full` byte-identical.
2. **Binary disk-buffer sink.** Swap the `flush` target to a binary/mmap disk buffer for hundreds of columns; Python reads via the same `swap_view_results` contract (returns a mapped region instead of an in-memory array).
3. **Profile (tz) output.** Apply the same record pattern to the depth–time output, opt-in due to memory.
4. **Other ~14 CSV input companions.** Continue the decouple-to-config move (bottom-boundary, drainage owl, irrigation, nutrients, warm-restart profiles, meteo detail/rain) using `read_csv_table_text`.

---

## Self-Review notes

- **Spec coverage:** record type + add_row (Task 1), flush-every-N seam (Tasks 1–2), per-instance on state (Task 3), inlist-only sieve reused from `vars%iyes` (Task 4 — the record gets exactly the `vals` the CSV writes), in-memory full read from Python (Tasks 5–6), byte-identical gate (Task 4 Step 7). Binary buffer is explicitly deferred (follow-up 2).
- **Byte-identical risk** is confined to Task 4 Step 4 (moving the `vals` build outside the headless guard). The computation is identical; only its *location* moves. Step 7 is the gate.
- **No new globals:** the record is a component of `swap_state_t` (per-instance). The existing `vars`/`scalar_w` globals are untouched this pass (retiring them is follow-up 1).
- **Open verification at execution time:** confirm `vars%head`/`vars%unit`/`vars%Nnodes` field names and `state%timecontrol%{headless,flRunEnd,t1900,flprintshort,date}` exist as used; confirm `swap_init/run_step/close` signatures for the Task 4 harness. These were true at plan-writing time (`csv_output.f90`, `swap_mod.f90`) but re-check if the files have moved.
