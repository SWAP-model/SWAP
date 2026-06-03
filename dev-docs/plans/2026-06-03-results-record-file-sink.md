# Output Sinks: In-Memory Record vs Streaming CSV — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

> **Supersedes** the earlier "record file sink" draft. Rationale: the existing per-step CSV writer already *streams* to disk (OS-buffered, periodically flushed) — it does not hold rows in RAM, so it is already memory-bounded and byte-identical. Re-implementing CSV writing inside the record would be reinventing a working writer and risking the byte-identical regression. Instead: one gather feeds two purpose-built sinks, chosen by config.

**Goal:** Make the in-memory results record **optional** (a config switch), and **expose the streaming CSV file's path to Python** — so the three use cases each get the efficient sink without redundant work or byte-identical risk.

**Architecture:** The per-step value gather is unchanged. It feeds (a) the existing streaming CSV writer when a file is wanted, and/or (b) the in-memory record when NumPy arrays are wanted. A new `[output.csv].results_in_memory` flag (default off) decides whether the record is built/populated at all — so standalone and huge runs skip it (bounded memory, the file is the buffer), while pyswap turns it on for zero-copy arrays. The existing CSV writer is untouched.

**Tech Stack:** Modern Fortran (gfortran, `iso_fortran_env`), pFUnit, meson, Python `ctypes`/`pandas`.

**Three use cases this delivers:**

| Use case | `results_in_memory` | CSV file (`csv_enabled` + not headless) | Results to Python |
|---|---|---|---|
| Standalone | 0 (off) | on | the CSV file (as today) |
| pyswap, modest run | 1 (on) | off (headless) | zero-copy NumPy via `swap_view_results` |
| pyswap, huge run | 0 (off) | on | read the streamed CSV via `swap_output_filepath` + pandas |

**Scope (this plan):** the `results_in_memory` switch + the CSV-path accessor + demo proof. NOT in scope: any change to the CSV format/writer (untouched), Parquet (a later Python-side upgrade off the same file/arrays), and multi-instance de-globalization of the writer (a separate plan — the record is already per-instance; the global writer is the remaining blocker, handled there).

**Iron law:** Byte-identical regression vs swap420gf. The CSV writer is untouched; `results_in_memory` defaults to 0, so regression cases (which set neither it) build no record and write the identical CSV. Verify with `pixi run -e test check-full`.

**Build/test commands:**
- Build: `pixi run build-linux`
- Unit tests: `pixi run -e test test-pfunit` (look for `OK (N tests)`)
- Fast gate: `pixi run -e test check-fast`
- Full regression: `pixi run -e test check-full`
- **Clean rebuild MANDATORY after any `src/state/*.f90` or `src/config/*` schema change:** `pixi run clean` then build.

---

## File Structure

- **Modify** `src/config/output_csv_config.f90` — add `results_in_memory` (int, default 0).
- **Modify** `src/io/toml/read_output_csv_toml.f90` — read `[output.csv].results_in_memory`.
- **Modify** `src/state/timecontrol_state.f90` — snapshot `results_in_memory`. (schema change → clean rebuild)
- **Modify** `src/io/csv_output.f90` — gate the record init / `add_row` / finalize on `results_in_memory`.
- **Modify** `src/bindings/swap_capi_mod.f90` — add `swap_output_filepath` (the streamed `_output.csv` path).
- **Modify** `prototype/pyswap_inmemory_demo.py` — verify both modes.
- **Modify** `tests/unit/io/test_results_wiring.pf` — the headless integration test must set `results_in_memory=1` (record is now opt-in).

---

## Task 1: `results_in_memory` config option

**Files:**
- Modify: `src/config/output_csv_config.f90`, `src/io/toml/read_output_csv_toml.f90`
- Test: `tests/unit/io/toml/test_read_output_csv_toml.pf` (suite `test_read_output_csv_toml_suite`)

- [ ] **Step 1: Write the failing test** — append to `tests/unit/io/toml/test_read_output_csv_toml.pf` (READ the file first to match its `use`/parse style):

```fortran
@test
subroutine test_read_output_csv_results_in_memory()
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
      '[output.csv]' // nl // 'enabled = 1' // nl // 'results_in_memory = 1' // nl

   call toml_load(doc, text)
   doc_ptr => doc
   call read_output_csv_toml(doc_ptr, cfg, errors)

   @assertFalse(errors%has_fatals())
   @assertEqual(1, cfg%results_in_memory)
end subroutine test_read_output_csv_results_in_memory
```

- [ ] **Step 2: Run, verify FAIL** — `pixi run -e test test-pfunit 2>&1 | grep -iE "results_in_memory|Error|Fail:"`. Expected: compile error (`results_in_memory` not a member) → add field (Step 3); then assertion fail → add reader (Step 4).

- [ ] **Step 3: Add the field** — in `src/config/output_csv_config.f90`, inside `type :: output_csv_config_t`, after `enabled_tz` add:
```fortran
      !! Build the in-memory results record (state%results) for Python array
      !! access. Off (0) => no record (standalone / huge runs stream to the CSV
      !! file instead, keeping memory bounded). On (1) => record is populated.
      integer :: results_in_memory = 0
```

- [ ] **Step 4: Read it** — in `src/io/toml/read_output_csv_toml.f90`, after the existing `enabled` read, add (mirror the helper used for `enabled`):
```fortran
      call get_optional_int_with_default(csv_sec, 'results_in_memory', config%results_in_memory, &
                                         0, 'output.csv.results_in_memory', errors)
```
(READ the file to confirm `get_optional_int_with_default` is in scope — it is used for `enabled`.)

- [ ] **Step 5: Run, verify PASS** — `pixi run -e test test-pfunit 2>&1 | grep -iE "results_in_memory|Fail:|OK \("`. Expected PASS.

- [ ] **Step 6: Commit**
```bash
git add src/config/output_csv_config.f90 src/io/toml/read_output_csv_toml.f90 tests/unit/io/toml/test_read_output_csv_toml.pf
git commit -m "feat(config): add [output.csv].results_in_memory toggle"
```

---

## Task 2: Snapshot `results_in_memory` onto timecontrol state

**Files:**
- Modify: `src/state/timecontrol_state.f90`

- [ ] **Step 1: Add the field** — in `src/state/timecontrol_state.f90`, near `csv_enabled` (line ~224) add:
```fortran
      integer :: results_in_memory = 0  !! snapshot of config%output_csv%results_in_memory
```

- [ ] **Step 2: Seed it** — in the same file's seed routine (near line ~313 where `csv_enabled` is seeded) add:
```fortran
      self%results_in_memory = config_output_csv%results_in_memory
```

- [ ] **Step 3: Clean rebuild + suite**
```bash
pixi run clean && pixi run -e test test-pfunit 2>&1 | tail -3
```
Expected: `OK (N tests)` (field inert so far). NOTE: the existing `test_results_populated_headless` in `tests/unit/io/test_results_wiring.pf` will FAIL after Task 3 because the record becomes opt-in — that test is updated in Task 3 Step 2. For now it still passes (Task 3 not yet applied).

- [ ] **Step 4: Commit**
```bash
git add src/state/timecontrol_state.f90
git commit -m "feat(state): snapshot results_in_memory onto timecontrol"
```

---

## Task 3: Make the record opt-in (gate on `results_in_memory`)

**Files:**
- Modify: `src/io/csv_output.f90` (`csv_out_header` record-init block, `csv_out_write` `add_row`, `csv_out_close` finalize)
- Test: `tests/unit/io/test_results_wiring.pf`

- [ ] **Step 1: Update the integration test to opt in** — in `tests/unit/io/test_results_wiring.pf`, in `test_results_populated_headless`, after `call swap_init(cfg, state, config)` and the `state%timecontrol%headless = .true.` line, add:
```fortran
   state%timecontrol%results_in_memory = 1
```
(Setting it post-init works here because `csv_output_init`'s record block reads the flag during init — so set it BEFORE the run loop, immediately after `swap_init`. If `csv_output_init` already ran inside `swap_init`, instead add an assertion-friendly path: see Step 3's gate reads `state%timecontrol%results_in_memory` at init time, so this test must inject it via the config. SIMPLER: keep the test as a headless run but assert the record is EMPTY when the flag is off, and add a second test that the flag is what's needed — see Step 1b.)

- [ ] **Step 1b: Replace the test with two clear cases** — replace `test_results_populated_headless` with:

```fortran
@test
subroutine test_record_off_by_default()
   use funit
   use swap_state_mod,  only: swap_state_t
   use swap_config_mod, only: swap_config_t
   use swap_mod,        only: swap_init, swap_run_step, swap_close
   type(swap_state_t)  :: state
   type(swap_config_t) :: config
   character(*), parameter :: cfg = 'tests/swap-cases/toml/1.hupselbrook/swap.toml'
   integer :: guard
   state%timecontrol%headless = .true.
   call swap_init(cfg, state, config)
   state%timecontrol%headless = .true.
   guard = 0
   do while (.not. state%timecontrol%flRunEnd .and. guard < 1000000)
      call swap_run_step(state, config); guard = guard + 1
   end do
   call swap_close(state, config)
   ! hupselbrook does not set results_in_memory => record stays empty.
   @assertEqual(0, state%results%nrows)
end subroutine test_record_off_by_default
```

(The "on" path is exercised end-to-end by the Python demo in Task 5, which injects `results_in_memory=1` via the config string — the cleanest way to set it before `csv_output_init` runs during `swap_init`.)

- [ ] **Step 2: Run, verify FAIL** — `pixi run -e test test-pfunit 2>&1 | grep -iE "record_off|Fail:"`. Expected: FAIL — today the record is always populated, so `nrows` is 36, not 0.

- [ ] **Step 3: Gate the record in `csv_output.f90`** — three edits:

In `csv_out_header`, wrap the record-init `block` (the one building `rnames`/`runits` and calling `state%results%init`) in:
```fortran
      if (state%timecontrol%results_in_memory == 1) then
         block
            ! ... existing rnames/runits build + call state%results%init(...) ...
         end block
      end if
```

In `csv_out_write`, gate the `add_row` call:
```fortran
         if (state%timecontrol%results_in_memory == 1) &
            call state%results%add_row(time%t1900, vals(1:ncount))
```
(The `vals`/`ncount` flatten stays where it is — it is also needed by the CSV file write. Only the `add_row` is gated.)

In `csv_out_close`, gate the finalize:
```fortran
      if (state%timecontrol%results_in_memory == 1) call state%results%finalize()
```

- [ ] **Step 4: Run, verify PASS** — `pixi run -e test test-pfunit 2>&1 | grep -iE "record_off|Fail:|OK \("`. Expected: PASS.

- [ ] **Step 5: Byte-identical gate** — `results_in_memory` defaults 0, CSV writer untouched:
```bash
pixi run -e test check-full 2>&1 | grep -iE "Results:"
```
Expected: `Results: 11 passed, 0 failed, 7 known-divergence (xfail)`.

- [ ] **Step 6: Commit**
```bash
git add src/io/csv_output.f90 tests/unit/io/test_results_wiring.pf
git commit -m "feat(io): make the in-memory results record opt-in (results_in_memory)"
```

---

## Task 4: `swap_output_filepath` C-ABI accessor

**Files:**
- Modify: `src/bindings/swap_capi_mod.f90`

- [ ] **Step 1: Add the accessor** — after `swap_results_columns`, add (builds the same path `csv_out_header` opens: `pathwork // outfil // '_output.csv'`):
```fortran
   !> Path of the streaming scalar CSV (_output.csv); written when csv output is
   !! enabled and not headless. Lets a huge pyswap run read the file instead of
   !! holding results in RAM.
   function swap_output_filepath(buf, n) result(ierr) bind(C, name='swap_output_filepath')
      character(kind=c_char), intent(out) :: buf(*)
      integer(c_int),  value, intent(in)  :: n
      integer(c_int)                      :: ierr
      character(len=512) :: path
      integer :: i, m
      path = trim(capi_state%timecontrol%pathwork) // &
             trim(capi_state%timecontrol%outfil) // '_output.csv'
      m = len_trim(path)
      do i = 1, min(m, n - 1)
         buf(i) = path(i:i)
      end do
      buf(min(m, n - 1) + 1) = c_null_char
      ierr = 0
   end function swap_output_filepath
```

- [ ] **Step 2: Build + symbol check**
```bash
pixi run build-linux 2>&1 | grep -iE "error|FAILED" | head
nm -D builddir/libswap_bmi.so | grep swap_output_filepath
```
Expected: clean build; symbol exported (`T`).

- [ ] **Step 3: Commit**
```bash
git add src/bindings/swap_capi_mod.f90
git commit -m "feat(bindings): swap_output_filepath accessor for the streamed CSV"
```

---

## Task 5: Demo — verify both modes (arrays + streamed file)

**Files:**
- Modify: `prototype/pyswap_inmemory_demo.py`

- [ ] **Step 1: Memory mode opt-in** — READ the demo. In `worker("mem")`, after reading the toml bytes and BEFORE `swap_initialize_from_toml_string`, inject the flag so the record is built:
```python
        toml += b"\n[output.csv]\nresults_in_memory = 1\n"
```
(The `mem` worker is headless, so no file is written; the record is the only output.) Confirm the existing `Results record (in-memory)` block still prints a non-zero shape and `record mem==disk: PASS`. NOTE the `disk` worker (BMI `initialize(path)`) does NOT get this injection, so its record is empty — change the `mem==disk` record comparison to compare the `mem` record against the `disk` **water balance only**, OR also inject the flag into a small toml the disk worker writes; SIMPLEST: in `worker("disk")` skip the record comparison (its record is empty by design) and keep only the water-balance comparison for disk.

- [ ] **Step 2: Add a file-mode check** — add a `worker("file")` branch + a `run_mode("file")` call in `main`. The file worker runs hupselbrook **not headless**, with `results_in_memory = 0` (default) and CSV enabled, then reads the streamed `_output.csv`:
```python
def file_worker():
    lib = load_lib()
    lib.swap_output_filepath.argtypes = [ctypes.c_char_p, ctypes.c_int]
    lib.swap_set_headless(0)
    lib.swap_clear_config_files()
    for name in COMPANIONS:
        with open(os.path.join(CASE, name), "rb") as fh:
            lib.swap_attach_config_file(name.encode(), len(name), fh.read(),
                                        os.path.getsize(os.path.join(CASE, name)))
    with open(os.path.join(CASE, "swap.toml"), "rb") as fh:
        toml = fh.read()
    rc = lib.swap_initialize_from_toml_string(toml, len(toml))
    assert rc == 0, f"file init rc={rc}"
    run_to_end(lib)
    pbuf = ctypes.create_string_buffer(4096)
    lib.swap_output_filepath(pbuf, len(pbuf))
    path = pbuf.value.decode()
    import csv
    rows = list(csv.reader(open(path)))
    return {"file": os.path.basename(path), "datarows": len(rows) - 1,
            "header0": rows[0][0] if rows else ""}
```
Wire it: add `elif sys.argv[2] == "file": print(SENTINEL + json.dumps(file_worker()))` to the worker dispatch, and in `main`:
```python
    fr = run_mode("file")
    print("\nStreaming CSV (huge-run path):")
    print(f"  wrote {fr['datarows']} data rows to {fr['file']}")
    ok = ok and fr["datarows"] == mem["_results_nrows"]
    print("  file rows == record rows:", "PASS" if fr["datarows"] == mem["_results_nrows"] else "FAIL")
```

- [ ] **Step 3: Run the demo**
```bash
python3 prototype/pyswap_inmemory_demo.py 2>&1 | tail -25
```
Expected: water-balance table; `Results record (in-memory)` (memory mode, non-zero shape, PASS); `Streaming CSV (huge-run path)` showing the same row count read from `_output.csv`; final `RESULT: PASS`.

- [ ] **Step 4: Commit**
```bash
git add prototype/pyswap_inmemory_demo.py
git commit -m "feat(prototype): verify memory-mode arrays and streamed-CSV huge-run path"
```

---

## Follow-ups (NOT in this plan)

1. **Multi-instance handle core** — run N columns/process. The remaining blocker is the *global* CSV writer (`scalar_w`) and the global `vars` gather buffer; de-globalize them onto per-instance state (the record is already per-instance) and tag per-column filenames. Separate plan.
2. **Parquet** — if huge-run CSVs get too bulky, have Python pull the in-memory chunks (or post-process the CSV) into Parquet with `pyarrow`. Python-side; no Fortran change.

---

## Self-Review notes

- **Spec coverage:** standalone (record off, CSV as today — default, byte-identical); pyswap modest (record on → zero-copy arrays — Tasks 1-3 + demo); pyswap huge (record off, read streamed CSV — Task 4 + demo). All three in the table at the top map to a task.
- **Byte-identical:** the CSV writer is not touched; `results_in_memory` defaults 0 so regression builds no record and writes the identical CSV. Gate is Task 3 Step 5.
- **Efficiency:** modest runs get zero-copy arrays (no file); huge runs stream to the existing writer (bounded memory, no RAM accumulation); standalone stops building an unused record. No redundant CSV-writing code.
- **No new globals; multi-instance-ready:** the record is per-instance; this plan adds only a per-instance flag. The one remaining global (the CSV writer) is explicitly the next plan's job.
- **Re-verify at execution:** `get_optional_int_with_default` in `read_output_csv_toml`; the `csv_enabled` seed site in `timecontrol_state`; the `csv_out_header` record-init block location; `capi_state%timecontrol%{pathwork,outfil}` field names — all true at plan-writing time.
