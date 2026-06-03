# Results-Record File Sink Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give the per-instance in-memory results record a configurable **file sink** so large pyswap runs can flush results to disk every N rows (bounded memory) and read them back with pandas — without touching the byte-identical legacy CSV output.

**Architecture:** The record (`results_state_t`, already per-instance on `swap_state_t`) gains a *sink mode*: `memory` (default — accumulate in RAM, Python reads zero-copy at the end) or `file` (own a `csv_writer_t`, write the buffered chunk to a tagged `*_results.csv` on every flush, never grow RAM). The mode comes from config (`[output.csv].results_sink`), so it is set before the record is initialised — no post-init timing problem. The filename carries an instance *tag* (empty for single-instance now; the future multi-instance work populates it), so each instance writes its own file. The legacy `_output.csv` path is left exactly as-is, so regression stays byte-identical.

**Tech Stack:** Modern Fortran (gfortran, `iso_fortran_env`), pFUnit, meson, Python `ctypes`/`pandas` for the demo.

**Scope (this plan):** ONLY the file sink + its config + wiring + a file accessor. NOT in scope: retiring the legacy `_output.csv` writer (a separate byte-identical-risky cleanup), Parquet (a later Python-side upgrade), and multi-instance (a separate plan — this plan only makes the sink *ready* for it via the tag).

**Iron law:** Byte-identical regression vs swap420gf (`pixi run -e test check-full`). The legacy `_output.csv` is untouched; the new `_results.csv` is only written when `results_sink="file"`, which no regression case sets — so regression is unaffected by construction. Verify it anyway.

**Build/test commands:**
- Build: `pixi run build-linux`
- Unit tests: `pixi run -e test test-pfunit` (look for `OK (N tests)`)
- Fast gate: `pixi run -e test check-fast`
- Full regression: `pixi run -e test check-full`
- **Clean rebuild MANDATORY after any `src/state/*.f90` or `src/config/*` schema change:** `pixi run clean` then build.
- `.gitignore` has a `result*` glob — any new `results*`/`*_results.csv` path needs a `!negation` line to stage, and test fixtures writing `*_results.csv` should write to a temp dir.

---

## File Structure

- **Modify** `src/state/results_state.f90` — add `sink_mode`, `tag`, `filepath`, an owned `type(csv_writer_t) :: writer`, `header_written`. `init` opens the file + writes the header in file mode; `flush` branches (file → write chunk + reset; memory → grow store); `finalize` closes the writer. Adds `use csv_writer_mod` (state→io, same one-way dependency precedent set by `config%meteo%meteo` using `meteo_csv_mod`).
- **Modify** `src/config/output_csv_config.f90` — add `results_sink` field (`'memory'`/`'file'`, default `'memory'`).
- **Modify** `src/io/toml/read_output_csv_toml.f90` — read `[output.csv].results_sink`.
- **Modify** `src/state/timecontrol_state.f90` — snapshot `results_sink` + add `results_tag` (default `''`). (schema change → clean rebuild)
- **Modify** `src/io/csv_output.f90` — `csv_out_header` passes sink mode/tag/path to `state%results%init`.
- **Modify** `src/bindings/swap_capi_mod.f90` — add `swap_results_filepath` so Python finds the file in file mode.
- **Tests:** `tests/unit/state/test_results_state.pf` (file-sink unit test), `tests/unit/io/test_results_wiring.pf` (file-mode integration).

---

## Task 1: File sink on `results_state_t`

**Files:**
- Modify: `src/state/results_state.f90`
- Test: `tests/unit/state/test_results_state.pf`

Note for implementer: `csv_writer_t` (in `src/io/csv_writer.f90`, module `csv_writer_mod`) has type-bound `open(path, errors)`, `header(names[, units])`, `row(values[, leading])`, `flush()`, `close()`. `row` prepends `leading` + a comma then joins `values` with commas. We write `t1900` as the leading column.

- [ ] **Step 1: Write the failing test** — append to `tests/unit/state/test_results_state.pf`:

```fortran
@test
subroutine test_results_file_sink_writes_csv()
   use funit
   use iso_fortran_env,   only: real64
   use results_state_mod, only: results_state_t

   type(results_state_t)         :: r
   character(len=*), parameter   :: path = 'builddir/test_results_sink.csv'
   integer                       :: u, ios, nlines
   character(len=4096)           :: line

   call r%init([character(len=24) :: 'RAIN', 'TACT'], &
               [character(len=12) :: 'cm', 'cm'], flush_every=2, &
               sink_mode=1, filepath=path)
   call r%add_row(100.0_real64, [1.0_real64, 2.0_real64])
   call r%add_row(101.0_real64, [3.0_real64, 4.0_real64])
   call r%add_row(102.0_real64, [5.0_real64, 6.0_real64])
   call r%finalize()

   ! File mode does NOT grow the in-memory store.
   @assertEqual(0, size(r%values, 1))
   ! Total rows written is tracked.
   @assertEqual(3, r%nrows)

   ! Header + 3 data rows = 4 lines; header starts with the t1900 column.
   open(newunit=u, file=path, status='old', action='read', iostat=ios)
   @assertEqual(0, ios)
   nlines = 0
   do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (nlines == 0) then
         @assertTrue(index(line, 't1900') == 1)
         @assertTrue(index(line, 'RAIN') > 0)
      end if
      nlines = nlines + 1
   end do
   close(u)
   @assertEqual(4, nlines)
end subroutine test_results_file_sink_writes_csv
```

- [ ] **Step 2: Run test, verify it FAILS to COMPILE first** (init has no `sink_mode`/`filepath` args yet): `pixi run -e test test-pfunit 2>&1 | grep -iE "results_state|Error|sink_mode"`. Expected: a compile error about the new dummy arguments. That's the RED for this step; proceed to implement so it compiles and then fails on assertions, then passes.

- [ ] **Step 3: Add the fields + extend `init`/`flush`/`finalize`** in `src/state/results_state.f90`.

Add `use csv_writer_mod, only: csv_writer_t` near the top (after `use iso_fortran_env`).

Add these components to the `type :: results_state_t` (after the existing `chunk_v` field):
```fortran
      ! Sink: 0 = memory (grow times/values), 1 = file (flush chunk to disk).
      integer :: sink_mode = 0
      character(len=:), allocatable :: tag        !! instance tag in the filename
      character(len=:), allocatable :: filepath   !! file-mode output path
      type(csv_writer_t)            :: writer      !! owned writer (file mode)
      logical :: header_written = .false.
```

Replace `results_state_init` with (adds optional `sink_mode`, `tag`, `filepath`):
```fortran
   subroutine results_state_init(self, names, units, flush_every, sink_mode, tag, filepath)
      use error_mod, only: error_collection_t
      class(results_state_t), intent(inout) :: self
      character(len=*),       intent(in)    :: names(:)
      character(len=*),       intent(in)    :: units(:)
      integer, optional,      intent(in)    :: flush_every
      integer, optional,      intent(in)    :: sink_mode
      character(len=*), optional, intent(in) :: tag
      character(len=*), optional, intent(in) :: filepath
      type(error_collection_t) :: errs
      character(len=24), allocatable :: hdr(:)
      integer :: k

      self%ncols = size(names)
      self%nrows = 0
      self%chunk_n = 0
      self%header_written = .false.
      if (present(flush_every)) self%flush_every = max(1, flush_every)
      if (present(sink_mode))   self%sink_mode   = sink_mode
      self%tag      = ''
      self%filepath = ''
      if (present(tag))      self%tag      = tag
      if (present(filepath)) self%filepath = filepath
      self%col_names = names
      self%col_units = units
      allocate(self%chunk_t(self%flush_every))
      allocate(self%chunk_v(self%flush_every, self%ncols))

      if (self%sink_mode == 1) then
         ! File sink: open + write header now; do NOT keep rows in RAM.
         allocate(self%times(0))
         allocate(self%values(0, self%ncols))
         call self%writer%open(trim(self%filepath), errs)
         allocate(hdr(self%ncols + 1))
         hdr(1) = 't1900'
         do k = 1, self%ncols
            hdr(k + 1) = self%col_names(k)
         end do
         call self%writer%header(hdr)
         self%header_written = .true.
      else
         allocate(self%times(0))
         allocate(self%values(0, self%ncols))
      end if
   end subroutine results_state_init
```

Replace `results_state_flush` with (branch on sink mode):
```fortran
   subroutine results_state_flush(self)
      class(results_state_t), intent(inout) :: self
      real(real64), allocatable :: t2(:), v2(:,:)
      character(len=32) :: lead
      integer :: old, add, i
      if (self%chunk_n == 0) return
      add = self%chunk_n

      if (self%sink_mode == 1) then
         ! File sink: write each buffered row (t1900 as the leading column),
         ! count it, and drop it from RAM.
         do i = 1, add
            write(lead, '(F0.6)') self%chunk_t(i)
            call self%writer%row(self%chunk_v(i, :), leading=trim(lead))
         end do
         call self%writer%flush()
         self%nrows   = self%nrows + add
         self%chunk_n = 0
      else
         ! Memory sink: grow the in-memory store.
         old = self%nrows
         allocate(t2(old + add))
         allocate(v2(old + add, self%ncols))
         if (old > 0) then
            t2(1:old)    = self%times
            v2(1:old, :) = self%values
         end if
         t2(old+1:old+add)    = self%chunk_t(1:add)
         v2(old+1:old+add, :) = self%chunk_v(1:add, :)
         call move_alloc(t2, self%times)
         call move_alloc(v2, self%values)
         self%nrows   = old + add
         self%chunk_n = 0
      end if
   end subroutine results_state_flush
```

Replace `results_state_finalize` with:
```fortran
   subroutine results_state_finalize(self)
      class(results_state_t), intent(inout) :: self
      call self%flush()
      if (self%sink_mode == 1) call self%writer%close()
   end subroutine results_state_finalize
```

- [ ] **Step 4: Add `csv_writer_mod` to the unit-test build** — `tests/unit/meson.build` `pfunit_extra_sources` already lists `'../../src/io/csv_writer.f90'` (it is used by other suites). Confirm it is present; if not, add it before `results_state.f90`.

- [ ] **Step 5: Run the test, verify PASS** — `pixi run -e test test-pfunit 2>&1 | grep -iE "file_sink|Fail:|OK \("`. Expected: `OK (N tests)`, no failures. Existing `results_state` tests (memory mode, default `sink_mode=0`) must still pass — they call `init` without the new args, which default to memory.

- [ ] **Step 6: Commit**
```bash
git add src/state/results_state.f90 tests/unit/state/test_results_state.pf
git commit -m "feat(state): add file sink to results_state_t (flush chunk to CSV)"
```

---

## Task 2: `results_sink` config option

**Files:**
- Modify: `src/config/output_csv_config.f90`, `src/io/toml/read_output_csv_toml.f90`
- Test: `tests/unit/config/test_read_output_csv_toml.pf` (the suite is `test_read_output_csv_toml_suite`)

- [ ] **Step 1: Write the failing test** — append to `tests/unit/io/toml/test_read_output_csv_toml.pf` (READ the file first to match its `use`/setup style; it parses a TOML snippet via `read_output_csv_toml`):

```fortran
@test
subroutine test_read_output_csv_results_sink()
   use funit
   use tomlf,                 only: toml_table, toml_load
   use output_csv_config_mod, only: output_csv_config_t
   use read_output_csv_toml_mod, only: read_output_csv_toml
   use error_mod,             only: error_collection_t

   type(toml_table), allocatable, target :: doc
   type(toml_table), pointer             :: doc_ptr
   type(output_csv_config_t)             :: cfg
   type(error_collection_t)              :: errors
   character(len=*), parameter :: nl = char(10)
   character(len=*), parameter :: text = &
      '[output.csv]' // nl // 'enabled = 1' // nl // 'results_sink = "file"' // nl

   call toml_load(doc, text)
   doc_ptr => doc
   call read_output_csv_toml(doc_ptr, cfg, errors)

   @assertFalse(errors%has_fatals())
   @assertEqual('file', trim(cfg%results_sink))
end subroutine test_read_output_csv_results_sink
```
Register it if the file uses the standard pattern (it is already in `pf_files`/`testSuites.inc`, so just adding the subroutine is enough).

- [ ] **Step 2: Run, verify FAIL** — `pixi run -e test test-pfunit 2>&1 | grep -iE "results_sink|Error|Fail:"`. Expected: compile error (`results_sink` not a member) → add the field (Step 3), then assertion fail (default not read) → add reader (Step 4).

- [ ] **Step 3: Add the field** — in `src/config/output_csv_config.f90`, inside `type :: output_csv_config_t`, after the `inlist_tz` field add:
```fortran
      !! In-memory results record sink: 'memory' (accumulate in RAM) or
      !! 'file' (flush every N rows to <outfil>[<tag>]_results.csv).
      character(len=:), allocatable :: results_sink
```
And in `output_csv_config_finalize` (or wherever defaults are applied — READ the file), ensure a default: if `.not. allocated(self%results_sink)` then `self%results_sink = 'memory'`. If there is no finalize body that runs always, set the default in the reader (Step 4) instead.

- [ ] **Step 4: Read it** — in `src/io/toml/read_output_csv_toml.f90`, after the existing `inlist`/`inlist_tz` reads, add:
```fortran
      call get_optional_string_with_default(csv_sec, 'results_sink', config%results_sink, &
                                            'memory', 'output.csv.results_sink', errors)
```
(Use the same helper the other string fields use; READ the file to confirm `get_optional_string_with_default` is in scope.)

- [ ] **Step 5: Run, verify PASS** — `pixi run -e test test-pfunit 2>&1 | grep -iE "results_sink|Fail:|OK \("`. Expected PASS.

- [ ] **Step 6: Commit**
```bash
git add src/config/output_csv_config.f90 src/io/toml/read_output_csv_toml.f90 tests/unit/io/toml/test_read_output_csv_toml.pf
git commit -m "feat(config): add [output.csv].results_sink (memory|file)"
```

---

## Task 3: Snapshot sink mode + tag onto timecontrol state

**Files:**
- Modify: `src/state/timecontrol_state.f90`

- [ ] **Step 1: Add the fields** — in `src/state/timecontrol_state.f90`, in the type (near `csv_inlist`, line ~226) add:
```fortran
      character(len=16)  :: results_sink = 'memory'  !! snapshot of config%output_csv%results_sink
      character(len=64)  :: results_tag  = ''        !! instance tag for the _results.csv filename
```

- [ ] **Step 2: Seed `results_sink` from config** — in the same file's init/seed routine (near line ~315 where `csv_inlist` is seeded), add:
```fortran
      if (allocated(config_output_csv%results_sink)) self%results_sink = config_output_csv%results_sink
```
`results_tag` stays `''` (the future multi-instance work sets it per column).

- [ ] **Step 3: Clean rebuild + confirm suite green**
```bash
pixi run clean && pixi run -e test test-pfunit 2>&1 | tail -3
```
Expected: `OK (N tests)`, no failures (fields are inert so far).

- [ ] **Step 4: Confirm byte-identical**
```bash
pixi run -e test check-fast 2>&1 | grep -iE "Results:"
```
Expected: `Results: 4 passed, 0 failed`.

- [ ] **Step 5: Commit**
```bash
git add src/state/timecontrol_state.f90
git commit -m "feat(state): snapshot results_sink + results_tag onto timecontrol"
```

---

## Task 4: Wire `csv_out_header` to configure the record's sink

**Files:**
- Modify: `src/io/csv_output.f90` (the record-init `block` in `csv_out_header`, added by the prior arc)

- [ ] **Step 1: Read** the record-init `block` in `csv_out_header` (it builds `rnames`/`runits` and calls `call state%results%init(rnames(1:nc), runits(1:nc))`).

- [ ] **Step 2: Pass sink mode/tag/path to `init`** — replace that `call state%results%init(...)` line with:
```fortran
         block
            integer :: sink_mode
            character(len=:), allocatable :: rfile
            sink_mode = 0
            if (trim(state%timecontrol%results_sink) == 'file') sink_mode = 1
            rfile = trim(state%timecontrol%pathwork) // trim(state%timecontrol%outfil) // &
                    trim(state%timecontrol%results_tag) // '_results.csv'
            call state%results%init(rnames(1:nc), runits(1:nc), &
                                    sink_mode=sink_mode, tag=trim(state%timecontrol%results_tag), &
                                    filepath=rfile)
         end block
```
(Keep the surrounding `rnames`/`runits` build exactly as-is.)

- [ ] **Step 3: Clean rebuild + run suite**
```bash
pixi run clean && pixi run -e test test-pfunit 2>&1 | tail -3
```
Expected: `OK (N tests)`, no failures.

- [ ] **Step 4: Byte-identical gate** — the legacy `_output.csv` path is untouched and no regression case sets `results_sink="file"`, so:
```bash
pixi run -e test check-full 2>&1 | grep -iE "Results:"
```
Expected: `Results: 11 passed, 0 failed, 7 known-divergence (xfail)`. If any non-xfail fails, STOP and investigate.

- [ ] **Step 5: Commit**
```bash
git add src/io/csv_output.f90
git commit -m "feat(io): configure results record sink (memory/file) from config"
```

---

## Task 5: `swap_results_filepath` C-ABI accessor

**Files:**
- Modify: `src/bindings/swap_capi_mod.f90`

- [ ] **Step 1: Add the accessor** — after `swap_results_columns` (the prior arc's accessor), add:
```fortran
   !> Path of the _results.csv written in file-sink mode; empty in memory mode.
   function swap_results_filepath(buf, n) result(ierr) bind(C, name='swap_results_filepath')
      character(kind=c_char), intent(out) :: buf(*)
      integer(c_int),  value, intent(in)  :: n
      integer(c_int)                      :: ierr
      integer :: i, m
      m = 0
      if (allocated(capi_state%results%filepath)) m = len_trim(capi_state%results%filepath)
      do i = 1, min(m, n - 1)
         buf(i) = capi_state%results%filepath(i:i)
      end do
      buf(min(m, n - 1) + 1) = c_null_char
      ierr = 0
   end function swap_results_filepath
```

- [ ] **Step 2: Build + symbol check**
```bash
pixi run build-linux 2>&1 | grep -iE "error|FAILED" | head
nm -D builddir/libswap_bmi.so | grep swap_results_filepath
```
Expected: clean build; the symbol is exported (`T`).

- [ ] **Step 3: Commit**
```bash
git add src/bindings/swap_capi_mod.f90
git commit -m "feat(bindings): swap_results_filepath accessor for file-sink mode"
```

---

## Task 6: End-to-end file-sink verification (pandas-readable)

**Files:**
- Modify: `prototype/pyswap_inmemory_demo.py`

- [ ] **Step 1: Add a file-sink check** — READ the demo, then add a new function that runs hupselbrook in-memory with `results_sink="file"` injected into the TOML and confirms the `_results.csv` has one row per output step matching the memory-mode record. Append before `main`'s final `RESULT:` print:

```python
def run_file_sink_check():
    import csv, tempfile, os as _os
    lib = load_lib()
    lib.swap_results_filepath.argtypes = [ctypes.c_char_p, ctypes.c_int]
    lib.swap_set_headless(1)
    lib.swap_clear_config_files()
    for name in COMPANIONS:
        with open(_os.path.join(CASE, name), "rb") as fh:
            content = fh.read()
        lib.swap_attach_config_file(name.encode(), len(name), content, len(content))
    with open(_os.path.join(CASE, "swap.toml"), "rb") as fh:
        toml = fh.read()
    toml += b'\n[output.csv]\nresults_sink = "file"\n'
    rc = lib.swap_initialize_from_toml_string(toml, len(toml))
    assert rc == 0, f"file-sink init rc={rc}"
    run_to_end(lib)
    nr, nc = ctypes.c_int(0), ctypes.c_int(0)
    lib.swap_results_shape(ctypes.byref(nr), ctypes.byref(nc))
    pbuf = ctypes.create_string_buffer(4096)
    lib.swap_results_filepath(pbuf, len(pbuf))
    path = pbuf.value.decode()
    rows = list(csv.reader(open(path)))
    return {"nrows_record": nr.value, "file": path,
            "file_header": rows[0], "file_datarows": len(rows) - 1}
```

And in `main`, after the record block, call it (in a subprocess for singleton isolation — mirror `run_mode`; simplest is to add a `worker` branch `filesink` that prints the dict as JSON and a `run_mode('filesink')` in main). Print:
```python
    fs = run_mode("filesink")
    print("\nFile sink:")
    print(f"  wrote {fs['file_datarows']} data rows to {os.path.basename(fs['file'])}")
    print(f"  header[0] = {fs['file_header'][0]}")
    fs_ok = fs["file_datarows"] == mem["_results_nrows"] and fs["file_header"][0] == "t1900"
    print("  file rows == record rows:", "PASS" if fs_ok else "FAIL")
    ok = ok and fs_ok
```
(Add a `filesink` branch to `worker` that returns `run_file_sink_check()` as the sentinel JSON.)

- [ ] **Step 2: Run the demo**
```bash
python3 prototype/pyswap_inmemory_demo.py 2>&1 | tail -25
```
Expected: the existing water-balance + record blocks, then a `File sink:` block showing N data rows written to `_results.csv`, `header[0] = t1900`, `file rows == record rows: PASS`, and final `RESULT: PASS`.

- [ ] **Step 3: Commit**
```bash
git add prototype/pyswap_inmemory_demo.py
git commit -m "feat(prototype): verify results_sink=file writes a pandas-readable _results.csv"
```

---

## Follow-ups (NOT in this plan)

1. **Multi-instance handle core** — run N columns/process; the ensemble sets each instance's `results_tag` so the file sinks don't collide. The per-instance record + tagged file sink built here are the readiness for it.
2. **Retire the legacy `_output.csv` writer** — make the record the single source for standalone too (byte-identical, or regenerate fixtures). Separate, byte-identical-risky.
3. **Parquet sink** — for efficiency at hundreds of columns, have Python pull chunks (`swap_view_results` + a reset call) and write Parquet with `pyarrow`. Keeps Fortran simple.

---

## Self-Review notes

- **Spec coverage:** file sink on the record (Task 1); config option memory|file (Task 2); per-instance snapshot + tag for multi-instance readiness (Task 3); wiring from config (Task 4); Python finds the file (Task 5); end-to-end pandas-readable verification (Task 6). Memory mode (the existing default) is unchanged — existing `results_state` tests pass with default args.
- **Byte-identical:** legacy `_output.csv` write path is untouched; the new `_results.csv` only appears when `results_sink="file"`, which no regression case sets. Task 4 Step 4 is the gate.
- **No new globals:** the writer is a component of the per-instance record; the tag lives on `state%timecontrol`. Multi-instance-ready.
- **Layering:** `results_state` (state) now `use`s `csv_writer_mod` (io) — a one-way state→io dependency, consistent with the `config%meteo%meteo`→`meteo_csv_mod` precedent set in the prior arc.
- **Open verification at execution time:** confirm `csv_writer_t` bindings (`open/header/row/flush/close`), `get_optional_string_with_default` availability in `read_output_csv_toml`, the `output_csv_config` default site, and the `timecontrol_state` seed routine field names — all true at plan-writing time; re-check if files moved.
