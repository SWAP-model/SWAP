---
title: "BMI Phase 2 — Lean BMI + Python-Driven Foundations — Design Spec"
date: 2026-05-13
status: draft
arc: SS-BMI2 (BMI Phase 2)
builds-on: SS-DRV Phase 1 (2026-05-12), SS-TCM (2026-05-13)
relates-to: ADR 0023 (reruns retired — parameter sweeps move external)
---

# BMI Phase 2 — Lean BMI + Python-Driven Foundations

## Context

SS-DRV Phase 1 (2026-05-12) shipped a minimal BMI façade: `initialize` / `update` / `finalize` / `get_current_time` / `get_time_step` plus two sentinel variables in `get_value_double`. The remaining ~25 CSDMS BMI v2.0 methods are stubbed with `! BMI-STUB Phase 2` markers, returning `0`/empty.

Phase 1 also established a clean separation: `swap_modern` static lib for new `-std=f2018` code, `swap_legacy` for everything else, `libswap_bmi.so` shared library, and a working Python `cffi` hello-world. The infrastructure is in place; the variable registry and a richer Python surface are not.

Two distinct external use cases motivate Phase 2:

1. **Model coupling (imod_coupler / MODFLOW).** A standardised, narrow BMI surface that exchanges water-flux variables between SWAP soil columns and a MODFLOW groundwater model. CSDMS BMI v2.0 spec compliance is the contract.

2. **Python-driven ensemble research.** Build 10,000 SWAP simulations from Python (different parameters, same meteo), run them in parallel via `multiprocessing.Pool`, and consume outputs as numpy arrays — without any worker process touching the filesystem. This is the "pyswap as research tool" use case.

The first needs a strict, narrow BMI. The second needs a richer C ABI surface and a foundation of "all inputs and outputs are in-memory buffers, not paths". The two are independent in their consumers but share infrastructure (the shared library, the state record, the `swap_mod` driver). One arc delivers both.

This arc does NOT rewrite or modify the existing `pyswap` package. Instead, it ships a small Python demonstration script under `tests/cffi-demo/` that exercises the full Python-driven workflow (load TOML in Python, pass as a string, attach shared-memory meteo, run headless, read output buffers). A subsequent arc adapts pyswap to use the new foundations.

## Decision

Implement the full CSDMS BMI v2.0 method set on a **lean** variable registry (~6–10 coupling variables) in `swap_bmi_mod.f90`. Introduce a separate, richer C-binding module `swap_capi_mod.f90` for Python-driven workflows. Add four pieces of "in-memory I/O" infrastructure: in-memory TOML config (`initialize_from_toml_string`), in-memory meteo buffer (`swap_attach_meteo_buffer`), headless mode flag (`swap_set_headless`), and output-row state buffers populated by every `*Output(2, ...)`. Migrate 15 timecontrol-related globals from `variables.f90` into `state%timecontrol` (BMI prerequisite + per-instance isolation).

Both binding modules ship in the same shared library. The Python demo script proves the end-to-end workflow; the pyswap rewrite is deferred.

## Scope summary

Approximately 960 lines of Fortran across:

- ~150 lines: full CSDMS v2.0 metadata methods in `swap_bmi_mod` + ~6–10 coupling variables
- ~150 lines: new `swap_capi_mod` with zero-copy accessors, scalar getters/setters, `bind(C)` struct returns
- ~30 lines: `initialize_from_toml_string` wrapping toml-f's string parser
- ~100 lines: `meteo_buffer_mod` + `readmeteo` external-buffer mode
- ~50 lines: headless flag + per-`*Output(N, ...)` guards
- ~410 lines: output sink refactor — row builders + state buffers for 9 subsystems
- ~70 lines: TimeControl globals sweep — 15 fields migrated from `variables.f90` to `state%timecontrol` plus 1 new `headless` field

Plus a ~50-line Python demo script.

## Architecture

### File structure

**New files:**

```
src/core/swap_capi_mod.f90              ← pyswap-direct C surface (zero-copy + scalars + struct returns)
src/io/toml/load_swap_config_string.f90 ← in-memory TOML entry point (wraps toml-f's toml_loads)
src/io/meteo_buffer_mod.f90             ← external meteo buffer storage + mode flag
tests/cffi-demo/swap_capi.h             ← C header for the capi layer
tests/cffi-demo/swap_bmi_extended.h     ← extended BMI header (replaces Phase 1's swap_bmi.h)
tests/cffi-demo/run_ensemble.py         ← end-to-end Python demo
tests/cffi-demo/meson.build             ← register demo as meson test in `cffi-demo` suite
```

**Significantly modified files:**

```
src/core/swap_bmi_mod.f90               ← full v2.0 methods (replace `! BMI-STUB Phase 2` stubs)
src/core/swap_main.f90                  ← thread headless awareness through log_init
src/io/readmeteo.f90                    ← branch on meteo_buffer_mod's mode flag
src/state/timecontrol_state.f90         ← +15 fields (tstart, tend, dtmin/max, period, etc.) + headless
src/state/{soilwater,atmosphere,heat,solute,surfacewater,tillage}_state.f90
                                        ← +output_row(:), output_columns(:), output_n_cols
src/io/swap_csv_output.f90              ← split row-build from file-write; honour headless
src/io/swapoutput.f90                   ← same split for the main swap balance output
src/core/timecontrol_mod.f90            ← read from state%timecontrol instead of `use variables`
src/core/initialize.f90                 ← drop init-zeros for 15 migrated globals
src/core/variables.f90                  ← drop 15 declarations
meson.build                             ← add swap_capi_mod to modern_sources; rename shared lib target
```

### Module boundaries

| Module | Responsibility | Consumers |
|---|---|---|
| `swap_bmi_mod` | CSDMS v2.0 spec compliance, ~6–10 coupling variables | imod_coupler, generic BMI tools |
| `swap_capi_mod` | Rich Python accessors, no spec constraints | pyswap (future), Python research scripts |
| `meteo_buffer_mod` | "Use external meteo buffer" mode flag + buffer pointer | `readmeteo` (reader-side), `swap_capi_mod` (writer-side) |

Both binding modules link into the same `libswap_bmi.so` (rename candidate: `libswap.so` — the "bmi" suffix is now misleading). They share `bmi_state` and `bmi_config` (module-level singletons inherited from Phase 1; multi-instance handle support is Phase 3).

## `swap_bmi_mod` — lean BMI surface

### Coupling variable registry (read-only via `get_value_double`)

| BMI name | Source | Shape | Units | Used by |
|---|---|---|---|---|
| `soil_water_content` | `state%soilwater%theta` | `[numnod]` | m³/m³ | imod_coupler |
| `pressure_head` | `state%soilwater%h` | `[numnod]` | cm | imod_coupler |
| `groundwater_level` | derived (deepest unsaturated node) | scalar | cm | MODFLOW exchange |
| `bottom_flux` | bottom-boundary flux (implementer locates exact state field via grep — likely `state%soilwater%inqdra(numnod)` or similar) | scalar | cm/d | MODFLOW exchange |
| `actual_evapotranspiration` | `state%soilwater%iqrot` + atmospheric contribution (`state%atmosphere%intr%iaintc` etc.) | scalar | cm/d | reporting |
| `recharge` | derived from bottom-node flux (same source as `bottom_flux`, sign convention may differ) | scalar | cm/d | MODFLOW exchange |
| `surface_runoff` | `state%surfacewater%runoff` if present, otherwise the cumulative runoff field — implementer verifies via grep before mapping | scalar | cm/d | reporting |
| `soil_temperature` | `state%heat%tsoil` | `[numnod]` | °C | reporting |

8 read variables. The mapping from BMI name to state path is a `select case (trim(name))` table inside `bmi_get_value_double`.

### Settable variables (via `set_value_double`)

| BMI name | Target | Used by |
|---|---|---|
| `groundwater_level_imposed` | sets bottom Dirichlet BC | MODFLOW → SWAP feedback |
| `bottom_flux_imposed` | sets bottom Neumann BC | MODFLOW → SWAP feedback |

2 settable variables. Every other BMI name returns rc=2 ("not settable") to prevent accidental state corruption.

### Full CSDMS v2.0 method set

All methods get real implementations (no more `! BMI-STUB Phase 2` markers):

**Control:** `update_until(time)` — loops `swap_run_step` until `t1900 >= time`. `get_component_name()` — returns `"SWAP"`.

**Info:** `get_input_var_names`, `get_output_var_names` — return hard-coded arrays of the names above. `get_input_item_count`, `get_output_item_count` — return `2` and `8`.

**Variable info:** `get_var_grid` (returns 0 for all SWAP variables — single grid), `get_var_type` (returns `"double"`), `get_var_units` (table lookup), `get_var_itemsize` (`8` for double), `get_var_nbytes` (size × itemsize), `get_var_location` (`"node"`).

**Time:** `get_start_time` (reads `state%timecontrol%tstart` — requires the timecontrol sweep), `get_end_time` (reads `state%timecontrol%tend`), `get_time_units` (returns `"d"`).

**Grid (one grid, id=0, the 1D vertical soil profile):** `get_grid_type` (`"uniform_rectilinear"`), `get_grid_rank` (`1`), `get_grid_size` (`numnod`), `get_grid_shape` (`[numnod]`), `get_grid_spacing` (table of `dz`), `get_grid_x`/`y` (single-zero array, irrelevant for 1D), `get_grid_z` (cumulative `dz`), `get_grid_node_count` (`numnod`).

**Getters/setters:** `get_value_int` and `get_value_float` return rc=1 (no SWAP variables are int or float at the BMI surface). `set_value_int`/`float` similarly. `get_value_double` extended per the registry above.

Total: ~25 new method bodies + 1 extended (`get_value_double`). Each is 5–20 lines.

## TimeControl globals sweep (BMI prerequisite + isolation)

15 fields migrated from `variables.f90` to `state%timecontrol`, plus 1 brand-new field (`headless`) added at the same time since it shares the same per-instance lifecycle:

| Field | Type | Default | Source at init |
|---|---|---|---|
| `tstart` | `real(real64)` | `0.0` | `config%simulation%tstart` |
| `tend` | `real(real64)` | `0.0` | `config%simulation%tend` |
| `dtmin` | `real(real64)` | `0.0` | `config%simulation%dtmin` |
| `dtmax` | `real(real64)` | `0.0` | `config%simulation%dtmax` |
| `period` | `integer` | `0` | `config%simulation%period` |
| `nprintday` | `integer` | `0` | `config%output%nprintday` |
| `flprintdt` | `logical` | `.false.` | `config%output%flprintdt` |
| `swheader` | `integer` | `0` | `config%output%swheader` |
| `swodat` | `integer` | `0` | `config%output%swodat` |
| `swres` | `integer` | `0` | `config%output%swres` |
| `swscre` | `integer` | `0` | `config%output%swscre` |
| `MaxIt` | `integer` | `0` | `config%simulation%MaxIt` |
| `MaxIterTime` | `integer` | `0` | `config%simulation%MaxIterTime` |
| `msteps` | `integer` | `0` | `config%simulation%msteps` |
| `flMaxIterTime` | `logical` | `.false.` | `config%simulation%flMaxIterTime` |
| `headless` | `logical` | `.false.` | (set by `swap_set_headless`) |

The `headless` field is included here since it's also a per-instance flag that belongs in `state%timecontrol`.

**Migration pattern** (same as SS-TCM Tasks 3, 5, 11, 12):

1. Add fields to `timecontrol_state_t` (schema additive).
2. Seed in `timecontrol_init` from config — transitional dual-write: keep the bare-global write, add the state-field write.
3. Migrate readers — across `timecontrol_mod`, `swap_mod`, `swap_csv_output`, `swapoutput`, `meteodt`, `numericalsolvers`, and `swap_capi_mod` (the new module).
4. Drop the transitional bare-global writes from `timecontrol_init`.
5. Delete the 15 bare-global declarations from `variables.f90` and the corresponding init-zeros from `initialize.f90`.

Subsystem-switch integers (`swsnow`, `swdra`, `swhea`, `swsolu`, `swetsine`, `swrain`, `swmetdetail`) are NOT migrated to state — they're read in `timecontrol_init` to derive the corresponding `flSnow` / `flDrain` / etc. flags, which are already in state. The integer switches can be read directly from `config%*` at init time, removing them from `timecontrol_mod`'s use list without state-side work.

## `swap_capi_mod` — Python-direct surface

Three patterns, applied as appropriate:

### Pattern A — zero-copy array view

```fortran
subroutine swap_view_array(name, ptr, n, ierr) bind(C, name='swap_view_array')
   character(kind=c_char), intent(in)  :: name(*)
   type(c_ptr),            intent(out) :: ptr   ! → array storage
   integer(c_int),         intent(out) :: n     ! number of elements
   integer(c_int),         intent(out) :: ierr  ! 0=ok, 1=unknown, 2=not allocated
end subroutine
```

Recognised names: `theta`, `h`, `tsoil`, `inqrot`, `q`, `dz` (per-node, `[numnod]`).

### Pattern B — scalar getter / setter

```fortran
function swap_get_scalar(name, value) bind(C, name='swap_get_scalar') result(ierr)
function swap_set_scalar(name, value) bind(C, name='swap_set_scalar') result(ierr)
```

Get registry: `tstart`, `tend`, `dt`, `t1900`, `daynr`, `iptra`, `iqrot`, `gwl`, `lai`, `crop_height`, `rooting_depth`.

Set registry (allowlisted): `lai`, `crop_height`, `rooting_depth` (external crop driver workflow). All other names return rc=2.

### Pattern C — `bind(C)` struct return for derived summaries

```fortran
type, bind(C) :: swap_water_balance_t
   real(c_double) :: rain
   real(c_double) :: evap_pot
   real(c_double) :: evap_act
   real(c_double) :: transp_pot
   real(c_double) :: transp_act
   real(c_double) :: runoff
   real(c_double) :: drain
   real(c_double) :: percolation
   real(c_double) :: storage_change
   real(c_double) :: balance_error
end type

function swap_get_water_balance(summary) bind(C) result(ierr)
```

One struct + one accessor per logical summary. Ships with `swap_water_balance_t`; others (crop summary, solute summary) deferred to a follow-on.

### Output-row accessors

```fortran
function swap_get_output_row(subsystem, ptr, names_ptr, n) bind(C) result(ierr)
   character(kind=c_char), intent(in)  :: subsystem(*)   ! "soilwater" | "swap" | "temperature" | ...
   type(c_ptr),            intent(out) :: ptr            ! → row(:)
   type(c_ptr),            intent(out) :: names_ptr      ! → column names array
   integer(c_int),         intent(out) :: n              ! current row size
   integer(c_int)                      :: ierr
end function
```

Returns pointers to `state%<subsys>%output_row(:)` and `state%<subsys>%output_columns(:)`. Python wraps them as numpy views. Updated after every `update()` call when an output step fires.

### Input-buffer attachment

```fortran
subroutine swap_attach_meteo_buffer(ptr, n_days, n_cols, ierr) bind(C)
```

Sets `meteo_buffer_mod`'s mode flag to `EXTERNAL_BUFFER`. The buffer pointer + dimensions are stored. `readmeteo` checks the mode flag and either parses CSV (mode=PATH) or consumes the buffer (mode=EXTERNAL_BUFFER).

Column ordering for the meteo buffer is canonical: `[date_offset, rain, tmin, tmax, et_ref, radiation, vapor, wind]`. Documented in `tests/cffi-demo/swap_capi.h`.

### Lifecycle additions

```fortran
function swap_initialize_from_toml_string(buf, n) bind(C) result(rc)
   character(kind=c_char), intent(in) :: buf(*)
   integer(c_int),  value, intent(in) :: n
   integer(c_int)                     :: rc
end function

function swap_set_headless(flag) bind(C) result(rc)
   integer(c_int), value, intent(in) :: flag   ! 0=normal, 1=headless
   integer(c_int)                    :: rc
end function
```

`swap_set_headless` must be called BEFORE `swap_initialize_from_toml_string` to skip the output-file `open()` calls in `*Output(1, ...)`.

## In-memory TOML config

`src/io/toml/load_swap_config_string.f90`:

```fortran
subroutine load_swap_config_from_string(toml_text, config, errors)
   character(len=*),         intent(in)    :: toml_text
   type(swap_config_t),      intent(inout) :: config
   type(error_collection_t), intent(inout) :: errors

   type(toml_table), allocatable, target :: doc
   type(toml_error), allocatable         :: terr

   call toml_loads(doc, toml_text, error=terr)    ! ← string parser; no file I/O
   if (allocated(terr)) then
      call errors%append(ERR_PARSE_MALFORMED_TOML, terr%message, "(in-memory)")
      return
   end if

   ! ... call each section reader (same as load_swap_config does)
end subroutine
```

The section readers (`read_general_toml`, `read_simulation_toml`, etc.) consume `toml_table` objects, which `toml_loads` produces — so no changes to them. The new function reuses the entire downstream pipeline.

A small refactor: extract the post-`toml_load` body of `load_swap_config` (lines 43+) into a shared helper that both `load_swap_config` (file path) and `load_swap_config_from_string` (in-memory) call. Avoids duplication.

## Meteo buffer

`src/io/meteo_buffer_mod.f90`:

```fortran
module meteo_buffer_mod
   use iso_c_binding, only: c_double, c_ptr, c_f_pointer
   implicit none
   private
   public :: meteo_mode_t, meteo_mode_path, meteo_mode_external_buffer
   public :: get_meteo_mode, set_meteo_mode, attach_meteo_buffer, get_meteo_buffer

   enum, bind(C)
      enumerator :: meteo_mode_path = 0
      enumerator :: meteo_mode_external_buffer = 1
   end enum

   integer,                save :: current_mode = meteo_mode_path
   real(c_double), pointer :: external_buffer(:,:) => null()
contains
   ! ... getters/setters/attach
end module meteo_buffer_mod
```

`readmeteo.f90` (sketch of new branch):

```fortran
if (get_meteo_mode() == meteo_mode_external_buffer) then
   call read_meteo_from_buffer(get_meteo_buffer(), state, ...)
   return
end if
! ... original CSV-reading code unchanged
```

`read_meteo_from_buffer` copies row data from the C-pointer buffer into `state%atmosphere%a*(:)` arrays in the same shape the CSV reader produces. The downstream consumers (meteoday, ProcessMeteoDay) see no difference.

Buffer ownership is Python-side; Fortran only holds a pointer. The buffer must remain alive (shared_memory block kept attached) for the duration of the simulation.

## Headless mode + output sink refactor

### State additions — one named buffer per output stream

Each output stream gets its own triplet (`row`, `columns`, `n_cols`) on the most natural state subrecord. The cffi accessor `swap_get_output_row("stream_name", ...)` dispatches on the stream name and returns a `c_ptr` to the right buffer.

The 9 streams and their natural homes:

| Stream name | Buffer location | Approx cols |
|---|---|---|
| `swap_balance` | `state%water_balance_row(:)` (new top-level fields on `swap_state_t`) | ~20 |
| `soilwater` | `state%soilwater%output_row(:)` | ~30 |
| `temperature` | `state%heat%output_row(:)` | `numnod`-dependent |
| `solute` | `state%solute%output_row(:)` | ~20 |
| `agetracer` | `state%solute%agetracer_row(:)` (sibling on solute state) | ~10 |
| `snow` | `state%atmosphere%snow_output_row(:)` (named to disambiguate from any atmosphere output) | ~10 |
| `surfacewater` | `state%surfacewater%output_row(:)` | ~25 |
| `crop` | `state%crop_output_row(:)` (new top-level fields on `swap_state_t`, since there is no `crop_state_t` today) | ~15 |
| `tillage` | `state%tillage%output_row(:)` | ~5 |

Each triplet has the same shape:

```fortran
real(c_double),    allocatable :: <name>_row(:)
character(len=32), allocatable :: <name>_columns(:)
integer                        :: <name>_n_cols = 0
```

Allocated and column names set in each subsystem's `*Output(1, ...)` (init). Row values overwritten in `*Output(2, ...)`. Released in `*Output(3, ...)` (close). The naming convention keeps multiple streams on the same subrecord (e.g. `solute` and `agetracer` both live on `state%solute`) unambiguous.

9 streams total.

### Per-subsystem refactor

Each `*Output(2, state)` is split:

```fortran
subroutine SoilWaterOutput(task, state)
   ...
   if (task == 2) then
      call build_soilwater_output_row(state)         ! Fills state%soilwater%output_row
      if (.not. state%timecontrol%headless) then
         call write_soilwater_output_row(state, swma_unit)  ! CSV-write path
      end if
   end if
end subroutine

subroutine build_soilwater_output_row(state)
   ! Computes the same numeric values that the legacy code wrote.
end subroutine

subroutine write_soilwater_output_row(state, unit)
   ! Pure formatter: state buffer → CSV row. Same format as legacy.
end subroutine
```

`*Output(1, state)`:
- Always allocates `output_row` + `output_columns` (size and names are static per subsystem).
- File `open()` is gated by `if (.not. state%timecontrol%headless)`.

`*Output(3, state)`:
- Deallocates `output_row`/`output_columns`.
- File `close()` gated by `if (.not. state%timecontrol%headless)`.

### Behaviour guarantees

- **Standalone binary**: byte-for-byte identical CSV outputs. The state buffer is filled, then formatted to file using the same format string as today. Regression byte-for-byte across all 5 cases must hold.
- **Cffi/headless mode**: state buffer is filled; no file operations occur. `swap_get_output_row` returns the buffer to Python.
- **Per-instance correctness**: each simulation has its own `state%*%output_row`. Two simultaneous instances (when multi-instance arrives in Phase 3) do not collide.

## Python demo script

`tests/cffi-demo/run_ensemble.py` (~50 lines):

```python
"""Demonstrates the BMI Phase 2 in-memory workflow end-to-end.

- Reads a SWAP TOML config file in Python.
- Reads the meteo CSV in Python.
- Passes both to SWAP via cffi (TOML as string, meteo as a shared
  memory buffer).
- Runs SWAP in headless mode (no CSV files created).
- Pulls the soil water output buffer after each `update()` and
  collects it into a numpy array.
- Prints a small summary (n steps, final theta[0], shape of collected
  output array).
"""

import pathlib, sys
import cffi
import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent


def main(case_dir: pathlib.Path, lib_path: pathlib.Path) -> int:
    ffi = cffi.FFI()
    ffi.cdef((HERE / "swap_bmi_extended.h").read_text())
    ffi.cdef((HERE / "swap_capi.h").read_text())
    lib = ffi.dlopen(str(lib_path))

    # Load TOML in Python
    toml_text = (case_dir / "swap.toml").read_text()
    toml_buf = toml_text.encode("utf-8")

    # Load meteo in Python (no SharedMemory yet — Phase 2 demos a single
    # process; the multiprocessing.SharedMemory wrapper is a pyswap-side
    # concern in a follow-on arc)
    met_csv = case_dir / "hupsel.met"
    met_arr = pd.read_csv(met_csv).to_numpy(dtype=np.float64)
    met_ptr = ffi.cast("double *", ffi.from_buffer(met_arr))

    # Wire up SWAP
    assert lib.swap_set_headless(1) == 0
    assert lib.swap_attach_meteo_buffer(met_ptr, met_arr.shape[0], met_arr.shape[1]) == 0
    assert lib.swap_initialize_from_toml_string(toml_buf, len(toml_buf)) == 0

    # Run and collect
    rows = []
    while True:
        ended = ffi.new("int *")
        # ... pull state%timecontrol%flRunEnd via swap_get_scalar
        rc = lib.update()
        assert rc == 0

        n = ffi.new("int *")
        ptr = ffi.new("double **")
        names = ffi.new("char ***")
        lib.swap_get_output_row(b"soilwater\0", ptr, names, n)
        row = np.frombuffer(ffi.buffer(ptr[0], n[0] * 8), dtype=np.float64).copy()
        rows.append(row)

        # ... loop termination via swap_get_scalar("flRunEnd")

    assert lib.finalize() == 0

    out = np.vstack(rows)
    print(f"Collected {out.shape[0]} output rows × {out.shape[1]} columns")
    print(f"Final theta[0] = {out[-1, 1]:.4f}")
    return 0


if __name__ == "__main__":
    case = pathlib.Path(sys.argv[1])      # e.g. tests/swap-cases/toml/1.hupselbrook
    lib  = pathlib.Path(sys.argv[2])      # path to libswap_bmi.so
    sys.exit(main(case, lib))
```

The script proves the workflow works end-to-end with zero filesystem touches for inputs and outputs. It's deliberately simple — no shared_memory, no parallelism. Those wrap the same primitives at the Python level (pyswap rewrite, future arc).

Registered as a meson test in the `cffi-demo` suite. Compared against a small fixture (final theta[0] within tolerance of the regression baseline).

## Verification strategy

Per memory `feedback_per_task_regression_gate.md`, every implementer subagent ends with:

1. `pixi run build-linux`
2. `pixi run test-pfunit`  (expected pattern depends on task — see Plan)
3. `pixi run check-fast`

Plus, after the output sink refactor: byte-for-byte `check-full` (5/5 cases) is the strong correctness signal.

After the full arc:

4. `pixi run -e test meson test -C builddir --suite bmi` — Phase 1 hello-world still passes.
5. `pixi run -e test meson test -C builddir --suite cffi-demo` — new end-to-end Python demo passes.
6. Manual: BMI compliance smoke test against `bmipy.Bmi` abstract base class (Python script that subclasses bmipy with our cffi and verifies all 25+ methods callable). Optional — not gating.

## Out of scope (deferred follow-on work)

- **Pyswap rewrite.** The Python package adopting the new C ABI is a separate arc. The demo script proves the foundation; pyswap's API design is independent.
- **Multi-instance support (handle-based).** Phase 2 uses the module-level singleton from Phase 1. Multiple in-process SWAP instances require refactoring `bmi_state`/`bmi_config`/`meteo_buffer_mod`'s module variables into an opaque-handle pattern. Phase 3.
- **In-memory buffers for CO2, soil hydraulic tables, crop tables.** Each follows the same pattern as the meteo buffer (mode flag + pointer + reader branch) but is a separate ~50-line addition each. Add when needed.
- **Output buffer for sub-daily writes.** Current scope covers daily output cadence. Sub-daily (`floutputshort`) follows the same row-builder pattern; implement when a consumer needs it.
- **Multi-grid BMI support.** Phase 2 declares one grid (id=0, the 1D vertical soil profile). Distinct grids for surface fluxes (scalar) or crop-state (different node count) are spec-supported but not implemented here.
- **`get_value_ptr_*` zero-copy BMI variants.** The spec includes pointer-based getters. SWAP can support them via the same `c_loc` pattern used in `swap_capi_mod`, but the BMI surface stays copy-only for now (matches what bmipy expects most often).
- **TimeControl `outdat(:)` / `outdatint(:)` output-date arrays.** Allocatable arrays in `variables.f90` — same migration pattern as the 15 scalars but with allocation handling. Defer to a small mini-arc.
- **`itnumb(:,:)` iteration counter.** Written by `headcalc` / numerical solver code that doesn't currently have `state` threaded. Migration touches deeper than this arc warrants.
- **The stop-100 IEEE summary cosmetic fix.** Captured in `docs/superpowers/specs/2026-05-13-fp-exception-summary-cleanup-note.md`.

## References

- SS-DRV Phase 1 spec — `docs/superpowers/specs/2026-05-12-driver-modernization-design.md`.
- SS-TCM spec — `docs/superpowers/specs/2026-05-13-timecontrol-modernization-design.md`.
- CSDMS BMI v2.0 — `https://bmi.readthedocs.io`.
- toml-f string parser — `subprojects/toml-f/src/tomlf/de.f90` (`interface toml_loads`).
- Phase 1 BMI stubs — `src/core/swap_bmi_mod.f90` lines 86–112 (the `! BMI-STUB Phase 2` markers).
- imod_coupler / MetaSWAP exchange variables — public Deltares documentation; the registry above mirrors that naming.
