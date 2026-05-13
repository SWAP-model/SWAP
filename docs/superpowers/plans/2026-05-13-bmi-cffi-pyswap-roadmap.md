---
title: "BMI / CFFI / pyswap Integration — Forward Roadmap"
date: 2026-05-13
status: draft
context: post-SS-BMI2 (working prototype shipped, cleanup arcs ahead)
---

# Forward Roadmap — BMI / CFFI / pyswap Integration

## Where we are today (post SS-BMI2, commit `4e9781c`)

Working prototype, all gates green:

- **BMI v2.0**: ~30 methods on `swap_bmi_mod` with 8 read + 2 settable coupling variables. Strict CSDMS compliance.
- **CFFI accessor surface**: `swap_capi_mod` — zero-copy array views, scalar get/set, `bind(C)` water-balance struct, 9-stream output-row dispatch, in-memory TOML config, meteo buffer attach, headless mode.
- **End-to-end demo**: `tests/cffi-demo/run_ensemble.py` proves the in-memory workflow on hupselbrook (29 388 update steps, zero filesystem touches for config, plausibility-checked output).
- **State migrations**: 16-field TimeControl sweep complete (tstart, tend, dt limits, output cadence, iteration limits, headless). swap_state now houses 9 output-row buffers.
- **Pixi tasks**: `test-bmi`, `test-cffi-demo`, `test-all` available alongside the existing `check-fast` / `check-full` gates.

The prototype is real but rough in places. The roadmap below organises the path from "working prototype" to "production-quality coupling/ensemble platform" by natural dependency order. Items within a phase are roughly interchangeable.

---

## Phase A — Hygiene cleanups (small arcs, weeks-not-months)

These are independent, can be picked up in any order, and each delivers value on its own. Recommended to clear at least the first two before Phase B starts, because they keep biting otherwise.

### A1. Meson cross-static-library `.mod` dependency declaration
**Cost:** ~half day. **Risk:** low.

The SS-BMI2 Task 2 incident (state schema additions silently breaking regression on incremental builds) showed that Meson's ninja graph doesn't propagate `.mod`-file dependencies across the `swap_modern` / `swap_legacy` static-library boundary. Today, the workaround is `rm -rf builddir` after any state schema change (memory `feedback_state_schema_clean_rebuild.md`).

The fix: declare explicit dependencies between `modern_sources` files and the state `.f90` files they consume. Likely a `declare_dependency(sources: [...])` block in `meson.build`. Once fixed, retire the memory note.

### A2. `stop 100` IEEE summary cosmetic fix
**Cost:** half hour. **Risk:** zero.

`docs/superpowers/specs/2026-05-13-fp-exception-summary-cleanup-note.md` documents this. Replace `stop 100` in `swap_main.f90` with a `bind(C)` libc `exit(100)` wrapper. Standard-compliant; suppresses the noisy IEEE flag summary that appears on direct binary invocation.

### A3. Stale doc-comment scrub
**Cost:** half hour. **Risk:** zero.

Several files (`meteoday.f90`, `drainage.f90`, and the state file head comments referencing the retired `flzerointr` / `flzerocumu`) still describe the old bare-global access pattern. Update the `!!` doc comments to reference `state%timecontrol%flZeroIntr` etc. Pure documentation, no behavior change.

### A4. Retire `ex_tlast` dead variable
**Cost:** half hour. **Risk:** zero.

Flagged in `variables.f90:62`. Originally tracked exchange time for the `handle_exchange` DLL pattern (retired in SS-DRV). Only `initialize.f90` still zeros it; nobody reads it.

### A5. Prune dead macropore branch in `timecontrol_reduce_dt`
**Cost:** half hour. **Risk:** zero (the branch is provably dead).

Per ADR 0040, `flMacroPore` is permanent `.false.`. The branch `if (flMacroPore .and. FlDecMpRat) then dt = dsqrt(...)` in `timecontrol_reduce_dt` will never fire. Delete it (preserved during SS-TCM as a verbatim copy).

### A6. Migrate `outdat(:)` / `outdatint(:)` arrays to `state%timecontrol`
**Cost:** ~1 day. **Risk:** low.

These output-schedule arrays were deferred from the SS-BMI2 sweep because they're allocatable (not scalars like the 15 fields we moved). The migration follows the same SS-* pattern but adds allocation handling at init and deallocation at close. Cleans up the last of the output-cadence globals.

### A7. Migrate `itnumb(:,:)` iteration counter to state
**Cost:** ~1–2 days. **Risk:** medium.

Written by `headcalc` / numerical solver code. Currently a bare-global `integer :: itnumb(100, 2)` in `variables.f90`. Threading state through the numerical solver path may touch deeper layers. Worth doing for parallelism (Phase C requires globals cleanup), but its own focused arc.

### A8. Globals cleanup, subsystem-by-subsystem
**Cost:** ~1–2 days per subsystem. **Total:** 1–2 weeks.

After SS-BMI2, ~34 files still `use variables`. The remaining bare-globals belong to: crop, irrigation, snow runtime fields, water-balance scalars not in cumu, and miscellany. Each subsystem can be migrated independently following the SS-TCM playbook. Required for Phase C (process-isolation works today; thread-level parallelism needs this).

---

## Phase B — BMI surface tightening (production-quality coupling)

These items complete the "best-effort" parts of SS-BMI2 to production quality. Drive them when a real consumer (imod_coupler integration test, pyswap user) surfaces a need.

### B1. Resolve BMI `set_value_double` BC plumbing
**Cost:** 2–3 days. **Risk:** medium (touches the bottom-boundary solver path).

Today `bmi_set_value_double('groundwater_level_imposed', ...)` writes `state%soilwater%h(numnod)` and `'bottom_flux_imposed'` writes `state%soilwater%qbot`. This is best-effort — depending on `swdrasur` / `swbotb` and the active boundary mode, the imposed values may or may not actually drive the next timestep correctly.

Required for real imod_coupler / MODFLOW exchange. Audit `BoundBottom` and the surrounding switches; ensure each imposed value flows through to the actual BC machinery. Add integration test against a known coupled-run baseline (need to coordinate with Deltares or MetaSWAP reference data).

### B2. Fill placeholder fields in `swap_water_balance_t`
**Cost:** ~1 day. **Risk:** low.

`swap_get_water_balance` returns placeholder zeros for `drain`, `storage_change`, `balance_error`. Find the canonical sources in `state%soilwater` (likely `cqdrain`, derived storage-change, computed balance-error) and wire them. The struct shape stays; only the assignment lines change.

### B3. Output body revival (SoluteOutput, SurfaceWaterOutput, CropOutput)
**Cost:** ~2 days each, separate arcs. **Risk:** medium per stream.

Three output streams ship as placeholders (N_COLS=1, no real values populated) because their original bodies were deleted by ADR 0009. The state buffer infrastructure is in place; values need to be assembled when a downstream consumer (BMI variable read, pyswap accessor) needs them.

Approach per stream:
1. Identify what columns the legacy output produced (git history before ADR 0009).
2. Map each column to its current state-record source.
3. Implement `build_<stream>_output_row(state)` to fill the buffer.
4. Increase N_COLS in `init_<stream>_output_buffer`.
5. Verify byte-for-byte regression remains 5/5.

Defer until someone actually needs the data.

### B4. Sub-daily output stream support
**Cost:** ~1 day. **Risk:** low.

Phase 2 covers daily output cadence. Sub-daily (`floutputshort` path) follows the same row-builder pattern but writes more frequently. Add the build path when sub-daily ensemble runs become a use case.

### B5. `get_value_ptr_*` zero-copy BMI variants
**Cost:** ~half day. **Risk:** low.

CSDMS BMI v2.0 spec includes pointer-based getters (`get_value_ptr_double`, etc.) for zero-copy access — same intent as our `swap_capi_mod`'s `swap_view_array`. Add as thin wrappers around the existing capi pattern, exposed under BMI naming. Useful when bmipy consumers prefer the pointer API.

### B6. Multi-grid BMI (separate scalar grid)
**Cost:** ~1 day. **Risk:** low.

Today all BMI variables share grid `id=0` (the 1D vertical soil column). Scalar variables (`groundwater_level`, `bottom_flux`, ET) are technically misrepresented as `numnod`-shape with only `[0]` filled. Add a second grid (`id=1`, `type="scalar"`) and dispatch in `get_var_grid` accordingly. Cleans up BMI metadata correctness; matters for tools that introspect grids.

---

## Phase C — Multi-instance + Fortran-level parallelism (Phase 3 of the original three-phase plan)

This is the architectural leap from "one SWAP instance per process" to "N SWAP instances in one process, optionally OpenMP-parallel."

### C1. Refactor `swap_capi_mod` and `swap_bmi_mod` to opaque handles
**Cost:** ~3–5 days. **Risk:** medium (interface change to existing callers).

Today both modules hold module-level `save` singletons (`capi_state`, `capi_config`, `bmi_state`, `bmi_config`). Refactor to:

```fortran
function swap_create() result(handle) bind(C, name='swap_create')
   type(c_ptr) :: handle
   type(swap_instance_t), pointer :: inst
   allocate(inst)
   handle = c_loc(inst)
end function

! Every existing entry point gains a handle argument:
function swap_initialize_from_toml_string(handle, buf, n) result(ierr) bind(C)
   type(c_ptr), value, intent(in) :: handle
   ...
```

The Python demo script + bmipy wrapper update to pass handles. Backwards-compatible by introducing a default-handle accessor for single-instance use cases.

### C2. `swap_ensemble_mod` — Fortran-side ensemble orchestrator
**Cost:** ~1 week. **Risk:** medium (requires globals cleanup — A8).

```fortran
type :: swap_ensemble_t
   type(swap_instance_t), allocatable :: columns(:)
end type
```

Provides batch operations (`ensemble%init_all`, `ensemble%step_all`, `ensemble%collect_outputs`). Optionally OpenMP-parallelises the per-column loops:

```fortran
!$omp parallel do
do i = 1, size(ens%columns)
   call swap_run_step(ens%columns(i)%state, ens%columns(i)%config)
end do
```

**Prerequisite:** every subsystem touched by `swap_run_step` must be thread-safe — i.e. NO `use variables` reads/writes of mutable globals during a step. This is what Phase A8 (globals cleanup) is for. Without A8, OpenMP races silently.

### C3. Ensemble C ABI + Python wrapper
**Cost:** ~2 days. **Risk:** low.

Expose `swap_ensemble_*` entry points (`create`, `add_column`, `step_all`, `get_outputs_for_column`, `destroy`) as `bind(C)`. Python wrapper iterates over columns, collects outputs into numpy stacks. Process-isolation (`multiprocessing.Pool`) remains the alternative for users who don't need shared memory between columns.

---

## Phase D — pyswap rewrite (Python-side arc)

This is the Python-side counterpart to the cffi foundations. Lives in the pyswap repository, not the SWAP repository — but worth scoping here so the SWAP side knows what to anticipate.

### D1. Replace ad-hoc cffi demo with class-based pyswap API
**Cost:** ~1 week. **Risk:** low (Python-only).

```python
from pyswap import SwapSim

sim = SwapSim()
sim.config.from_dict({...})         # build config in Python
sim.config.crop.rotation.append(...) # ergonomic per-field setters
sim.run(headless=True)               # internally serializes config to TOML
                                     # string + calls cffi
df = sim.outputs.soilwater           # pandas DataFrame from output_row buffer
```

Internally:
- pyswap holds a cffi handle.
- `sim.run` serializes config dataclass → TOML string → `swap_initialize_from_toml_string`.
- Output collection walks `swap_get_output_row` after each `update()`.
- Cleanup via `finalize` on `__del__` / context manager exit.

### D2. Shared-memory meteo helper
**Cost:** ~2 days.

```python
from pyswap import MeteoSharedBuffer
buf = MeteoSharedBuffer.from_csv("hupsel.met")    # reads once, puts in shm
with multiprocessing.Pool(8) as pool:
    pool.map(run_one_sim, parameter_ensemble, buf.handle)
```

Each worker attaches `buf` via name + calls `swap_attach_meteo_buffer` before init. Solves the "10 000 columns × same meteo = 10 000 wasted CSV parses" problem identified during brainstorming.

### D3. CO2 / soil / crop table shared buffers
**Cost:** ~2 days each. **Risk:** low.

Same pattern as meteo. Follow `meteo_buffer_mod` template:
- New Fortran module per input type (`co2_buffer_mod`, `soil_table_buffer_mod`, ...)
- Reader-side branch on mode flag
- cffi `swap_attach_<type>_buffer(ptr, ...)` entry point
- pyswap helper for shared-memory packing

Implement on demand — only meteo matters for typical ensembles; CO2 is small enough that the file read is fine.

### D4. Ensemble experiment API
**Cost:** ~1 week.

Once C2 (Fortran ensemble) lands, pyswap exposes a higher-level API:

```python
exp = pyswap.Experiment(base_config)
exp.sweep('soil.dtmin', [0.0001, 0.001, 0.01])
exp.sweep('crop.rotation[0].rdmax', np.linspace(50, 150, 20))
results = exp.run(workers=8, backend='multiprocess')  # or 'openmp' (C2)
```

Pivots cleanly between process-pool (each worker is a separate Python process loading its own cffi handle) and OpenMP (one process, ensemble of in-memory handles, parallel `update` loop).

---

## Phase E — External validation and adoption

This is what the foundations are for — the consumers that prove the architecture works.

### E1. imod_coupler integration test
**Cost:** ~1–2 days for the test wiring; ~1–2 weeks total counting BC plumbing iteration.

Coordinate with Deltares: take a small MetaSWAP + MODFLOW coupled scenario and re-run with our SWAP BMI in MetaSWAP's slot. Compare exchange variables (groundwater_level, recharge) timestep by timestep.

The actual integration may surface gaps in B1 (BC plumbing). Iterate.

Deliverable: SWAP listed as a supported BMI provider in imod_coupler docs.

### E2. eWaterCycle / generic BMI tool compatibility
**Cost:** ~1 day.

Run our library through `bmipy.Bmi`'s abstract base class verifier. Pass any tool that loads BMI providers (notably eWaterCycle's container ecosystem).

Deliverable: any tool that speaks BMI v2.0 can drive SWAP.

### E3. pyswap research-workflow case studies
**Cost:** open-ended.

Once pyswap has the new API (D1–D4), publish a notebook ensemble study — e.g., parameter sensitivity for a known dataset. Establishes the workflow as a research-grade alternative to "edit TOML, run binary, parse CSV."

---

## Suggested ordering for the next 1–2 months

If picking up tomorrow:

1. **A1 (Meson .mod deps)** — half day, eliminates a recurring footgun
2. **A2 + A3 + A4 + A5** — half day total, pure hygiene
3. **B1 (BC plumbing) + B2 (water balance fields)** — ~3 days, makes the BMI surface production-ready
4. **E1 (imod_coupler integration test)** — drives B1 iteration, ~1–2 weeks
5. **D1 (pyswap class API)** — Python-side, ~1 week, can run parallel to imod_coupler work
6. **A6 + A7 + A8 (globals cleanup arcs)** — sets up C2 (OpenMP ensemble)
7. **C1 + C2 + C3 (multi-instance + ensemble)** — once globals are clean
8. **D2 + D4 (pyswap ensemble + shm)** — once C lands

The split allows two streams of work — Fortran-side (A → B → C) and Python-side (D) — to proceed in parallel after Phase B opens up E1.

## Where this plan lives

This roadmap is a high-level forward plan. Individual items become formal specs + plans in `docs/superpowers/specs/` and `docs/superpowers/plans/` when picked up. Each numbered item (e.g. A1, B1) is a candidate arc — typically 1 day to 2 weeks of focused work.

## Files / artifacts inventory (post-SS-BMI2)

Reference points for anyone picking up this roadmap:

**Working entry points** (linked into `libswap_bmi.so`):
- `src/core/swap_bmi_mod.f90` — strict CSDMS BMI v2.0 surface
- `src/core/swap_capi_mod.f90` — pragmatic Python-direct surface

**Headers** (cffi consumes):
- `tests/bmi/swap_bmi.h` — BMI C declarations
- `tests/cffi-demo/swap_capi.h` — capi C declarations
- `tests/cffi-demo/swap_bmi.h` — symlink to bmi header (so cffi-demo gets both)

**Working demos**:
- `tests/bmi/hello_swap.py` — Phase 1 hello-world (`pixi run -e test test-bmi`)
- `tests/cffi-demo/run_ensemble.py` — Phase 2 end-to-end (`pixi run -e test test-cffi-demo`)

**Infrastructure modules**:
- `src/io/toml/load_swap_config_string.f90` — in-memory TOML loader
- `src/io/meteo_buffer_mod.f90` — meteo buffer mode flag + reader branch
- `src/state/timecontrol_state.f90` — 79-field record (incl. `headless` flag)

**Documentation**:
- `docs/superpowers/specs/2026-05-12-driver-modernization-design.md` (SS-DRV Phase 1)
- `docs/superpowers/specs/2026-05-13-timecontrol-modernization-design.md` (SS-TCM)
- `docs/superpowers/specs/2026-05-13-bmi-phase-2-design.md` (SS-BMI2 — this arc)
- `docs/superpowers/specs/2026-05-13-fp-exception-summary-cleanup-note.md` (A2)

**Memory notes**:
- `feedback_state_schema_clean_rebuild.md` — clean-rebuild requirement (A1 retires this once Meson deps fixed)
- `feedback_per_task_regression_gate.md` — verification discipline
