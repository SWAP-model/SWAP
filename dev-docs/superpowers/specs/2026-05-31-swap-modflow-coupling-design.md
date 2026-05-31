# SWAP → MODFLOW 6 Coupled Demo + Folded Cleanups — Design

**Date:** 2026-05-31
**Status:** Approved (brainstorming) — pending spec review → writing-plans
**Repo:** `/home/zawadzkim/Code/swap` (branch `development`)

## 1. Purpose & scope

Produce a **working, runnable example** of SWAP (unsaturated-zone column model)
coupled to MODFLOW 6 (saturated groundwater), driven by Deltares'
**imod_coupler** `swapmod` driver (which already exists on the `swap` branch).
This is the first concrete step of the forward arc toward parallelization,
Python-driven execution, and MODFLOW coupling, and it doubles as the forcing
function that pins down SWAP's library/XMI interface.

The MODFLOW conceptual model: **lateral flow between two channels**, with
**recharge delivered per cell from SWAP columns**. One SWAP column per MODFLOW
recharge cell; all columns are **homogeneous** (same soil, crop, meteo) but hold
**independent state** and receive a **per-cell groundwater head** from MODFLOW,
so their recharge response diverges — that divergence is the physics the demo
should show.

**Definition of done (smoke + qualitative):** the coupled run completes via
imod_coupler without crashing; the exchange arrays carry sensible non-zero
values; a plot shows the water table sitting between the two channel heads with
recharge-driven mounding, and cells at different heads showing different
recharge. No strict numerical assertion.

**Non-negotiables carried from CLAUDE.md:** gfortran only; `implicit none`;
`iso_fortran_env` kinds; **standalone SWAP regression stays byte-identical**
(the XMI/ensemble layer is additive — the CLI path is unchanged). No physics
changes (equations/coefficients/units/FP-ordering).

## 2. The contract (discovered, not invented)

`Deltares/imod_coupler` **`swap` branch** has a complete `swapmod` driver
(`drivers/swapmod/swapmod.py`) + `SwapWrapper` (`kernelwrappers/swap_wrapper.py`,
extends xmipy's `XmiWrapper`). MODFLOW 6 is the **leading clock**. Per daily
step the driver does:

```
heads → SWAP            swap_head[:] = mf6_head[:]
mf6.prepare_time_step;  delt = mf6.get_time_step()
swap.prepare_time_step(delt)
swap.prepare_solve(0); swap.solve(0); swap.finalize_solve(0)
storage + recharge → MF6   mf6_storage[:] = swap_storage[:]
                           mf6_recharge[:] = swap_volume[:] / delt
mf6.prepare_solve(1); for kiter: mf6.solve(1) until converged; mf6.finalize_solve(1)
both.finalize_time_step
```

SWAP is solved **once per step** (not re-solved inside the MF6 convergence
loop) → **loose sequential coupling**. Exchange is `[:] = [:]` (1:1, no
remapping) → SWAP column ordering must match MF6 recharge-cell ordering, which
is natural for a 1-row cross-section.

**C-ABI `libswap.so` must export** (via xmipy):
- XMI lifecycle: `initialize`, `prepare_time_step(dt)`, `prepare_solve(id)`,
  `solve(id)`, `finalize_solve(id)`, `finalize_time_step`, `finalize`,
  `get_current_time`, `get_end_time`, `get_version`, `report_timing_totals`.
- XMI variable-pointer protocol: `get_var_rank`/`get_var_type`/`get_var_shape`
  + `get_value_ptr` for three named arrays.

**Decision — SWAP-native names + thin fork.** SWAP exposes the three arrays
under readable native names (proposed: `gwl`, `qbot_volume`, `storage_coef`),
and we vendor a small fork of `swap_wrapper.py` (just the three `get_value_ptr`
strings) and `swapmod.py`. We own that tiny Python glue; the Fortran keeps clean
semantics.

**Current SWAP surface** (`src/core/swap_bmi_mod.f90`, `swap_capi_mod.f90`):
plain BMI (`initialize/update/update_until/get_value_double/set_value_double/…`)
over a single `capi_state`/`capi_config` singleton. **Gap:** no XMI sub-timestep
methods, no `get_value_ptr` pointer protocol, no ensemble, aborts via
`error stop`.

## 3. Approach (selected: A — minimal sequential ensemble)

One `libswap.so` owns an **ensemble**: `columns(:) : swap_state_t` (N
independent states), `configs(:) : swap_config_t` (M unique, read-only) +
`column_config(:)` index map (N → config; demo: all `1`), and a shared meteo
buffer. XMI `solve(0)` loops `swap_run_step(columns(i), configs(column_config(i)))`
**sequentially**.

Why this is small and safe (no full singleton retirement needed now):
- SWAP compute already threads state by argument.
- The remaining shared globals are safe under *sequential* looping:
  `O2_pars%current_state` is rebound per call; the meteo buffer is *meant* to be
  shared (homogeneous columns); per-column CSV file output is disabled in
  coupled mode.
- Full thread-safe multi-instance (retiring `O2_pars`/`csv_output`/`meteo_buffer`
  globals, stack-array → allocatable sweep) is the **later** parallelization arc,
  informed by this one. Not in scope here.

**Per-column config indirection is built in now** (index map + M unique
configs), so heterogeneous soils/crops are a later data change, not a
re-architecture. Demo uses M=1.

Rejected: **B** (full multi-instance retirement first — pulls a large arc into a
smoke demo); **C** (N=1 broadcast — fails the per-cell-GWL requirement, no
divergence, throwaway).

## 4. Components

### SWAP-side Fortran
- **`swap_ensemble_mod.f90`** (new) — `configs(:)`, `column_config(:)`,
  `columns(:)`, shared meteo buffer, and the three flat **`allocatable, target`**
  exchange arrays `gwl(:)`, `qbot_volume(:)`, `storage_coef(:)`. Procedures:
  `ensemble_init` (allocate N columns from configs, validate coupling
  preconditions), `ensemble_step` (sequential loop), `ensemble_finalize`.
- **`swap_xmi_mod.f90`** (new) — `bind(C)` XMI facade + variable-pointer
  protocol resolving native names → the three arrays. Reuses `swap_run_step`.
  `prepare_solve`/`finalize_solve` are thin (loose coupling); the column loop
  lives in `solve`.
- **Bottom-boundary prescribed-head seam** — apply `gwl(i)` (m→cm) as column
  *i*'s prescribed-GWL bottom BC (`swbotb=1` path) before its solve, reading the
  injected state field rather than a mandatory file; after solve read `qbot` →
  `qbot_volume(i)` (cm/d → volume via cell area, cm→m); compute
  `storage_coef(i)`. Scoped unification of today's split `qbot`/`gwl` logic.
- **Returnable errors** — the `fatalerr`/`abort_if_fatal` choke point
  (`src/error/error.f90:191`) + `dtutil.f90` raw `error stop` return a status in
  library mode instead of aborting; XMI methods return non-zero `rc`. Standalone
  CLI keeps `error stop` (it is the program, not the library).
- **Build** — meson `libswap.so` shared target (C-ABI + version symbol) + pixi
  task; CLI binary retained.

### Python-side (vendored in the demo dir)
- **`swap_wrapper.py`** (fork) — three `get_value_ptr` strings → native names.
- **`swapmod.py`** (fork) — same loop, native names; commented-out
  mapping/in-loop-re-solve bits resolved for our 1:1 sequential case.
- **`build_modflow.py`** (flopy) — see §6.
- **`run_coupled.py`** — orchestration: build MF6 inputs; emit SWAP inputs; write
  the imod_coupler config TOML (kernel `.so` paths + coupling block); run
  imod_coupler; collect + plot.

## 5. SWAP configuration (Hupselbrook with two deltas)

Base = `tests/swap-cases/toml/1.hupselbrook` (known-good). Changes for the
coupled case, both SWAP-supported:
- **Bottom boundary `swbotb = 6 → 1`** (prescribe groundwater level). In coupled
  mode the GWL is injected from `mf6_head` each step, not read from a file.
- **Drainage `swdra = 1 → 0`** (off). Lateral flow is MODFLOW's job; the
  `swap.dra.toml` is dropped from the coupled case.
- Unchanged: 2-layer Mualem-van Genuchten soil (`ksatfit ≈ 12.5 cm/d`), heat on,
  solute on, maize→potato→grass rotation, `283.csv` meteo. Shared across all N
  columns (M=1).
- Period: **full Hupselbrook 2002-2004** (~1096 daily steps).

## 6. MODFLOW 6 model (flopy, grounded to that soil)

- **DIS:** 1 layer × 1 row × N cols; `top = 0 m`, `botm ≈ −10 m`; unconfined
  (`NPF icelltype = 1`).
- **NPF:** `K ≈ 0.125 m/d` (= SWAP's ~12.5 cm/d ksat; same-soil aquifer
  assumption).
- **IC:** initial head ≈ `−0.75 m` (matches `gwli = −75 cm`).
- **STO:** transient; specific yield `sy` ≈ documented estimate.
- **CHD:** cells 1 and N (the two channels), heads `h0`/`h1` with a small
  gradient (e.g. −0.5 / −1.0 m).
- **RCHA:** interior cells — the coupling recharge package
  (`mf6_swap_recharge_pkg`).
- **TDIS:** transient, daily stress periods over the SWAP period.

## 7. Data flow (per daily step)

`MF6 head(:)` → driver `gwl[:] = mf6_head[:]` → `prepare_solve` applies
`gwl(i)` → column *i* bottom BC (m→cm) → `solve` loops columns sequentially
(each advances one day) → `qbot_volume(i)`, `storage_coef(i)` → driver
`mf6_recharge[:] = qbot_volume/delt`, `mf6_storage[:] = storage_coef` (cm→m) →
MF6 convergence loop → next day.

## 8. Key correctness risks (plan must tackle head-on)

1. **Recharge unit/area conversion** — SWAP `qbot` is cm/d; MF6 RCHA wants a
   rate; conversion through cell area + cm→m + the driver's `/delt`. Pin the
   exact convention; assert sign (downward percolation → positive recharge).
   *Highest-priority correctness item.*
2. **`storage_coef` semantics** — MetaSWAP feeds `sc1` so MF6 storage stays
   consistent with the unsaturated zone. For smoke, start with a documented
   fixed specific-yield-like constant; refine later. Do not block the demo on
   exact semantics.
3. **Bottom-BC mode** — validate at `ensemble_init` that every coupled column's
   config uses the prescribed-GWL bottom boundary (`swbotb=1`).
4. **Vertical datum alignment** — SWAP column bottom (−200 cm) vs MODFLOW water
   table; align so injected `gwl` and the SWAP profile are consistent. Setup
   detail to verify.
5. **Driver WIP** — the swap-branch driver has commented-out mapping /
   in-loop re-solve; verify it runs in our loose 1:1 sequential config; finish in
   our fork.

## 9. Error handling

The library never aborts: returnable status throughout the coupled path;
non-zero XMI `rc` → coupler reports and stops cleanly. The standalone CLI keeps
`error stop`.

## 10. Testing / definition of done

- **Standalone regression byte-identical** — `pixi run -e test check-fast`
  (and `check-full` at the end) stay green; the XMI/ensemble layer is additive,
  CLI path unchanged. This also guards the cleanups under the byte-identical law.
- **Coupled smoke** — run completes via imod_coupler; exchange arrays
  sensible/non-zero; plot shows the water table between `h0`/`h1` with
  recharge-driven mounding, and per-cell recharge differences. Qualitative.

## 11. Folded cleanups (separable sub-arcs, each ending `check-fast` green)

Structured so none can stall the coupling demo; each verified byte-identical.
1. **Drop MetaSWAP** — delete dormant `src/utils/dormant/sharedsimulation.f90`;
   delete the dead `msw1eic` Gash interception path in
   `src/atmosphere/interception.f90` (`swinter=2/3` already rejected by TOML
   validators → unreachable); drop the stale `swap_log.f90` comment. (Also
   removes the one `real(4)` island, the lone `write(*)+stop`, and a dead `!$OMP`
   block.)
2. **Delete other dead code** — empty `src/boundary/boundary_constants.f90` +
   fix its README; deprecated `soil_config%z_init`; dead `tillage.f90` debug
   writes (hardcoded units); empty `waterbalance.f90` `if` block.
3. **real(8) → real64 sweep** — scoped to subsystems this arc touches; FP-neutral;
   verified byte-identical.
4. **intent/pure annotations** — scoped to leaf kernels in touched subsystems
   (e.g. heat `Devries`, soil `WC_K_models`, drainage `divdra`/`Lev2Comp`).

## 12. Repo layout

```
tests/coupling/                 # new (alongside tests/bmi/, tests/cffi-demo/)
  build_modflow.py              # flopy model builder
  run_coupled.py                # orchestration
  imod_coupler_fork/            # vendored swap_wrapper.py + swapmod.py
  swap_config/                  # hupselbrook-derived coupled config (swbotb=1, swdra=0)
  imod_coupler.toml             # coupler config (kernel paths + coupling block)
  README.md
src/core/swap_xmi_mod.f90       # new
src/state/swap_ensemble_mod.f90 # new (or src/core/)
```

## 13. Sequencing (for writing-plans)

1. **Runnable pipe first:** ensemble + XMI facade + bottom-boundary seam +
   returnable errors + `libswap.so` build → driver fork → flopy model →
   orchestration → first coupled run.
2. **Cleanup sub-arcs** land alongside (independent; each `check-fast` green).
3. **Finalize:** `check-full` byte-identical; coupled smoke plot; README.

## 14. Out of scope (explicit)

- Full thread-safe multi-instance / stack-array → allocatable sweep (next arc).
- NetCDF I/O (delegated to Python/xarray per prior decision).
- Heterogeneous configs/forcing (machinery built; data deferred).
- Tight (within-MF6-iteration) coupling; quantitative/validated DoD.
