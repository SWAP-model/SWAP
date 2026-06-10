# Python orchestration & parallelization — assessment and plan

**Date:** 2026-06-11
**Goal (from the project brief):** move from the legacy single-column SWAP
executable to a fast, modern exe/library capable of true parallelization and
coupling with other models (MODFLOW), orchestrated from Python over a gridded,
netCDF-friendly interface, and drivable for PEST / sensitivity analysis.

This document assesses what already exists, names the gaps, and lays out a
phased plan. A working **first slice** ships with it:
`prototype/swap_orchestrator.py` (netCDF-driven, parallel, BMI-backed,
parameter-injectable; run `--selftest`).

---

## 1. Where we already are (better than the brief assumes)

The strangler-fig modernization has already built most of the substrate:

| Capability | Status | Where |
|---|---|---|
| Typed `swap_state_t` threaded by argument; no bare globals / `state%cfg` | done | `src/state/`, ADR 0045/0046 |
| In-memory, diskless init (push TOML + companions as bytes) | done | `swap_capi_mod.f90`: `swap_attach_config_file`, `swap_initialize_from_toml_string`, `swap_set_headless` |
| In-memory results record (zero-copy numpy view of the output matrix) | done | `swap_capi_mod.f90`: `swap_results_shape/columns/view_results/view_results_times`; opt-in `[output.csv] results_in_memory=1` |
| BMI 2.0 single-instance C-ABI (grid + get/set value + time) | done | `swap_bmi_mod.f90` (`libswap_bmi.so`) |
| XMI ensemble C-ABI for tight coupling (zero-copy exchange arrays) | done | `swap_xmi_mod.f90` (`libswap_xmi.so`): `get_value_ptr` on `gwl`, `qbot_volume`, `storage_coef` |
| In-process multi-column ensemble (sequential) | done | `swap_ensemble_mod.f90`: `columns(:)`, per-column `column_config` map, gwl injection / qbot extraction |
| SWAP↔MODFLOW6 loose coupling (MODFLOW leads the clock) | demo works | `tests/coupling/` (vendored imod_coupler fork, `swapmod.py`) |

**Design decisions already on record** (align the plan to these):
- **ADR 0043** — every procedure takes `state` explicitly, touches nothing
  outside it; in-process parallel multi-instance is the stated near-term goal.
- **Coupling design spec (2026-05-31)** — one `libswap.so` owns an ensemble
  `columns(:)`/`configs(:)` with an `N→M` `column_config` index map; the XMI
  `solve()` loops `swap_run_step` **sequentially** for now.
- **NetCDF I/O is explicitly delegated to Python/xarray** (coupling spec, line
  246). The Fortran side stays in-memory + CSV; Python owns gridded I/O.

So the missing piece is **not** the Fortran engine — it is a Python
orchestration layer above the existing C-ABI, plus retiring the last Fortran
globals that block *true* (threaded) parallelism.

---

## 2. Gaps to the goal

1. **No Python orchestration layer over the C-ABI.** `pyswap` today is a
   file-based input writer + `subprocess` runner + CSV/pandas parser
   (`pyswap/model/model.py:run`). It has **no** BMI/ctypes binding, no gridding,
   no ensembles, no netCDF. Its `run_parallel()` is embarrassingly-parallel
   subprocesses over scratch dirs.
2. **No gridded/netCDF in/out.** Nothing reads an xarray grid of per-cell
   parameters or writes a `(cell, time, var)` cube.
3. **Ensemble is Fortran-internal + sequential.** `swap_ensemble_mod` is not a
   bind(C) surface of its own (only reachable through XMI), shares one meteo
   buffer (homogeneous columns), and runs columns one-by-one.
4. **True (threaded) parallelism is blocked by Tier-3 module globals** — per the
   post-Phase-4 summary: `crop_config_global`, `SAVE` locals in
   `cropgrowth.f90:601,678` and `cropgrass_runtime.f90:94`, the
   `Wofost_Soil_Declarations` blanket-`save` module, and output file-handle
   `SAVE`s. Until these are retired, columns cannot share one process safely.
5. **No first-class PEST / sensitivity entry point.** The hooks exist (param
   injection via TOML, in-memory results readout) but nothing wires them into a
   sampling/calibration loop.

---

## 3. Two orchestration topologies (use the right one per job)

- **Offline / loose (this slice):** one process per column, columns fully
  independent, run via the **BMI single-instance** in-memory C-ABI, parallelized
  with a process pool. Sidesteps *all* Tier-3 global state (process isolation),
  supports **heterogeneous** per-column configs immediately, and is the natural
  fit for **PEST / sensitivity** (each "column" is a parameter sample) and for
  **offline gridded runs** (each column a cell, no inter-cell feedback within a
  step). This is what `prototype/swap_orchestrator.py` implements.

- **Online / tight (MODFLOW coupling):** all columns in **one** process via the
  **XMI ensemble**, stepped together each day, exchanging `gwl`/`qbot_volume`
  with MODFLOW by zero-copy pointer. Required when columns share a groundwater
  system that must converge within the timestep. Already demoed in
  `tests/coupling/`. Parallelizing *this* one needs the Tier-3 global retirement
  (Phase C) so the per-column loop can go `!$omp parallel do`.

The shared currency between both is the typed `swap_state_t` + the in-memory
C-ABI; netCDF is the Python-side interchange format for both.

---

## 4. Plan (phased; each phase independently shippable)

### Phase A — Orchestration core + netCDF I/O  ✅ first slice shipped
`prototype/swap_orchestrator.py`:
- `ColumnSpec(name, base_dir, overrides)` → `run_columns()` (process pool over
  the BMI in-memory C-ABI) → `ColumnResult` (zero-copy results matrix + water
  balance).
- `apply_overrides()` (tomlkit round-trip) injects per-column parameters without
  disturbing the rest of the TOML — the PEST/sensitivity hook.
- `results_to_dataset()` / `specs_from_dataset()` / `write_netcdf()` —
  xarray `(column, time, var)` cube in and out.
- Self-test proves a 3-column initial-GWL sweep runs in parallel, overrides take
  effect, the balance closes, and the netCDF round-trips.

**Next within A:** promote from `prototype/` into a package (e.g.
`pyswap.orchestrate` or a standalone `swap_orchestrate`), add a CLI
(`swap-orchestrate grid.nc --base case/ -o out.nc`), and an xarray meteo path
(per-column `283.csv` companions or the `swap_attach_meteo_buffer` C-ABI).

### Phase B — PEST / sensitivity front-ends
- A thin `design → specs` layer: SALib (Morris/Sobol) and a PEST `.pst`/
  run-manager adapter that maps parameters → `overrides` and result columns →
  objectives. Pure Python over Phase A; no Fortran change.
- Add `swap_set_scalar` allowlist entries (currently all rejected,
  `swap_capi_mod.f90:250`) for the handful of parameters worth setting *after*
  init (warm-start / data assimilation); most calibration stays init-time.

### Phase C — Retire Tier-3 globals → true threaded parallelism
Unblocks running the XMI ensemble columns concurrently in one process (and lets
Phase A drop to threads where fork overhead matters). Order by the post-Phase-4
blocker list:
1. `crop_config_global` → thread a `config` reference through the crop cluster.
2. `SAVE` locals in `cropgrowth.f90:601,678`, `cropgrass_runtime.f90:94` → move
   onto `crop_state_t`.
3. `Wofost_Soil_Declarations` blanket-`save` → per-instance record.
4. Output file-handle `SAVE`s → per-instance (or headless/in-memory only).
Each step is regression-gated (byte-identical) and independently verifiable with
a 2-column "same inputs ⇒ identical outputs, no cross-talk" test.

### Phase D — Promote the ensemble to a first-class gridded C-ABI
- Give `swap_ensemble_mod` its own bind(C) surface (init with **per-column**
  config paths / TOML strings + per-column meteo, step, read results) so Python
  can drive heterogeneous gridded ensembles directly (not only through XMI's
  homogeneous, `ensemble.txt`-driven path).
- Wire `results_in_memory` per column so the gridded cube is read zero-copy from
  one process once Phase C makes it safe.

### Phase E — pyswap refactor (the brief's "substantial refactor")
- Add a **library backend** to `pyswap.Model.run()`: instead of writing files +
  `subprocess`, serialize the model to a TOML string + companions and call the
  in-memory BMI C-ABI (keep the subprocess path as a fallback). Re-uses Phase A.
- Add a `pyswap.Grid` / `pyswap.Ensemble` that holds an xarray of per-cell
  component overrides over a base `Model`, and runs via Phase A/D.
- Output: return `xarray.Dataset` cubes alongside the existing pandas frames.

---

## 5. Risks / constraints discovered while building the slice

- **Byte-identity is the law.** The orchestrator changes *nothing* in the
  engine; it only drives the existing C-ABI. The one engine change made this
  session (the drainage first-step fix, separate commit) is regression-gated.
- **Process-singleton state.** The BMI/CAPI path uses one process-global
  `capi_state`, so **one process = one column**. The orchestrator therefore uses
  `multiprocessing` (spawn) — correct and safe today. Threaded single-process
  multi-column is Phase C work.
- **netCDF delegated to Python** (per the coupling spec) — so the cube assembly
  lives in `results_to_dataset`, not in Fortran. Good: keeps hot kernels I/O-free.
- **Heterogeneous meteo** per column isn't wired through the *XMI* ensemble yet
  (shared buffer); the *BMI offline* path already supports it (each column
  attaches its own `283.csv`).

---

## 6. Try-it

```bash
pixi run build-linux                                   # builds libswap_bmi.so
pixi run -e test python prototype/swap_orchestrator.py --selftest
```

The self-test runs a 3-column initial-GWL sweep in parallel, asserts the
overrides take effect and the water balance closes, writes + reloads a netCDF
cube. Extend `ColumnSpec(...).overrides` to sweep any TOML parameter; swap the
3 hand-written specs for `specs_from_dataset(xr.open_dataset("grid.nc"), base)`
to drive a real grid or a calibration design.
