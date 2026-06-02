# SWAP modernization sweep — findings & roadmap (2026-06-01)

A full-subsystem audit against the forward intent: **in-process multi-instance
execution, clean model/Python interoperability, extensibility, and clarity.**
Thirteen parallel read-only agents swept every subsystem; the highest-stakes
structural claims were verified directly against source (see
[Verification](#verification)). Findings carry `file:line` evidence and are
scored by the rubric below.

This is an **analysis document — no code was changed.** Roadmap items become
future brainstorm → plan → implement cycles.

---

## 1. Executive summary

**The single most important finding: the compute kernel is already
multi-instance clean — every global that blocks parallelization lives at the
*edges* (ambient services + bindings + output), not in the physics.**
`swap_mod`'s `swap_init/swap_run_step/swap_close` thread `state` and `config`
as arguments. The blockers are a short, well-bounded list of module-level
`SAVE` singletons in five places, plus a cluster of per-rotation `SAVE` state
in the crop dispatchers. That makes the multi-instance goal *tractable* — it is
a contained edge-cleanup, not a kernel rewrite.

Top findings, ranked (severity × objective weight: parallelization /
extensibility dominant, then clarity):

| # | Finding | Subsystem | Sev | Why it matters |
|---|---------|-----------|-----|----------------|
| 1 | **Shared ambient services** — `error_mod` (`global_errors`, `library_mode`, `fatal_was_raised`) and `swap_log` (module scalars) are process-global `SAVE`. Used by *every* subsystem. | error, core | Critical | One column's fatal/log state bleeds into all others. Blocks any N-instance run. Foundational — unblocks everything else. |
| 2 | **Binding singletons** — `capi_state`/`capi_config` (`save`) and the `swap_ensemble_mod` globals (`columns`, `configs`, `gwl`, `qbot_volume`, `storage_coef`, all `save`). No instance-handle model. | bindings, driver | Critical | Hard ceiling of one SWAP per process for BMI/CAPI; one ensemble per process for XMI. This is *the* multi-instance enabler. |
| 3 | **Crop per-rotation `SAVE` state** — `grass`/`wofost`/`cropfixed` runtime dispatchers + `wofost_soil_declarations` carry `SAVE`; `oxygenstress` uses a module-level `current_state` pointer. | crop | Critical | **Latent correctness bug, not just a future blocker:** the existing XMI ensemble only works because all columns share one config (identical crop). Distinct cropped columns — the pyswap case — silently clobber each other's phenology. |
| 4 | **No in-memory results accessor** — output is write-only CSV via module-global writers (`profile_w`/`scalar_w`, `save`); BMI exposes only 8 vars of ~112. | io | High | pyswap's "run hundreds of columns, collect results in memory" is impossible without disk round-trips + filename collisions. |
| 5 | **Three overlapping C-ABI facades** (BMI, XMI, CAPI) exporting duplicate `bind(C)` names → forced into separate `.so`s. | bindings | High | Architectural confusion; redundant maintenance. Drives the binding recommendation (§3). |
| 6 | **Blocking `read(*,*)` + file writes inside compute kernels** — `oxygenstress` numerical routines, `wofost` `outbalcrop*`, crop debug units 777/888. | crop | High | Deadlocks in headless/batch/parallel runs; serializes I/O. |
| 7 | **Oversized fixed-size state arrays** — `soilwater_state` boundary tables `2*MABBC` ≈ **4.7 MB/instance**; `crop_irrigation_state` ≈ 520 KB; `crop_oxygen_state` ≈ 280 KB. | state | High | Violates the allocatable-default rule; ×hundreds of instances is prohibitive and risks stack SIGSEGV. |
| 8 | **Pervasive hardcoded switch dispatch** — `swbotb` (×2 sites), `croptype` 1/2/3 (×6 in cropgrowth), `swinter`/`swredu`/`swcalt`/`swfrost`/stress models. | all compute | Medium | Adding a boundary/crop/process means editing N call sites. Recurring extensibility tax. |
| 9 | **Oversized routines** — `OxygenStress` (432 L), `wofost` task 2/3 (563 L), `grass` runtime (1317 L), `divdra` (337 L / 27 args), `headcalc` (573 L), `cropwofost_config` (1063 L). | crop, soilwater, drainage, config | Medium | Exceeds what fits in one mental/context window; blocks safe change. |
| 10 | **Dead / stranded code** — `agetracer` (won't compile after next refactor), `SWBR=1` unreachable body, wofost bulb/soybean/CO2 dead branches, `sptabulated` (5557 L), `jongvanlier`, `oxygenrepro`, `regrid`, `FillOxygenStress1/2`. | solute, crop, soilwater | Medium | False complexity; the half-ports are active liabilities. |

---

## 2. Scoring rubric

**Severity** (weighted toward the user's priorities):

- **Critical** — blocks correct multi-instance execution *today* or is the
  foundation other parallelization work depends on. Includes latent
  correctness risks under N distinct instances.
- **High** — blocks a named use case (pyswap many-column, imod_coupler) or is a
  hard parallelization/extensibility wall, but with a contained workaround.
- **Medium** — extensibility tax or clarity blocker; slows future work, no
  correctness/scaling wall.
- **Low** — local cleanup, naming, documentation.

**Effort** — S (≤1 commit), M (a focused arc, a few commits), L (multi-arc /
cross-file sweep).

**Objective** — which forward goal each finding most serves:
`parallelization` · `extensibility` · `clarity` · `interop`.

**Non-negotiable constraint on every item:** byte-identical regression vs
`swap420gf`. None of these touch equations, coefficients, units, or FP
operation ordering. Hidden-state removal and seam extraction are behavior-
preserving by construction (ADRs 0043–0048 precedent).

---

## 3. Binding-layer assessment (centerpiece)

### 3.1 What exists today (verified)

| Facade | Backed by | Instance model | Use case | Status |
|--------|-----------|----------------|----------|--------|
| **CAPI** (`swap_capi_mod`) | `capi_state`/`capi_config` singletons | **single, process-global** | low-level ABI; BMI built on it | TOML-string init ✓, zero-copy `swap_view_array` ✓, **setters empty** (all rejected, rc=2), water-balance struct partial (drain/storage/error stubbed 0) |
| **BMI** (`swap_bmi_mod`) | aliases CAPI singletons | **single** | CSDMS-style single-column coupling | ~25 `bind(C)` methods; 8 read vars + 2 settable, single 1-D grid; **subset-BMI, no multi-instance/multi-grid** |
| **XMI** (`swap_xmi_mod`) | `swap_ensemble_mod` (array of `columns`) | **one ensemble / process** (N columns) | imod_coupler / MODFLOW6 | **proven, tested** (`test_coupled_smoke`); 3 exchange vars (`gwl`→, `qbot_volume`←, `storage_coef`←); `get_value_ptr` zero-copy; daily step |

Key structural facts (all verified in source):
- The **kernel is stateless** — `swap_init/run_step/close` take `(state, config)`.
- **BMI and XMI export the same `bind(C)` symbol names** (`initialize`,
  `update`, …) → they *cannot* live in one `.so`; hence two libraries. This is
  the tell that two facades are doing one job.
- XMI's multi-column infrastructure is real but **single-config**:
  `allocate(configs(1))`, `column_config = 1`. Per-column config indirection
  exists but is unused. `storage_coef` is a hardcoded 0.15 placeholder.
- The standalone `swap` executable links the kernel directly; **bindings are
  not involved** — so standalone is unaffected by any binding decision.

### 3.2 What each real use case actually needs

**(A) pyswap — hundreds of columns for SA / auto-calibration.** Needs:
cheap N independent instances in one process; per-instance parameter override;
results pulled into memory (no per-column CSV files); ideally parallel.
*Today:* none of this — singleton state, empty setters, disk-only output.

**(B) imod_coupler / MODFLOW6.** Needs: step-wise ensemble drive, head→flux
exchange, error propagation. *Today:* works via XMI for the homogeneous case;
gaps are per-column config, variable dt, real specific-yield, and mapping
internal fatals to the XMI return code.

**(C) standalone.** Works; must keep working. Independent of the above.

Both (A) and (B) fundamentally need the **same** thing: state+config behind an
**opaque per-instance handle** instead of a module singleton. XMI's `columns(:)`
array is a half-step toward exactly that.

### 3.3 Recommendation

**Commit to a single handle-based instance core, exposed through *two* thin
facades — XMI (coupling) and a pyswap C-API (many-column) — and stop
maintaining BMI as a Fortran artifact.**

Concretely:

1. **One instance core.** Retire `capi_state`/`capi_config` and the
   `swap_ensemble_mod` singletons into an **instance registry**: an opaque
   handle (`c_ptr` / integer id) indexes an array of `(state, config)`. Both
   facades ride on it. This is the load-bearing change; everything else is a
   thin wrapper.

2. **Keep XMI** as the canonical coupling contract — it's proven, tested, and
   matches the xmipy/imod_coupler symbol surface. Reshape it as a thin facade
   over the instance core (it already is, modulo the singletons).

3. **Build a pyswap C-API** over the same core: handle-per-column, a populated
   parameter-setter allowlist, and an **in-memory results accessor** (§ io
   findings). Not CSDMS-BMI-shaped — pragmatic and zero-copy.

4. **Drop standalone Fortran BMI.** Its 8-var single-column surface is strictly
   dominated by XMI + the pyswap C-API. If CSDMS-BMI conformance is ever needed,
   generate it in **Python** over the C-API (the bmipy/xmipy pattern) — don't
   carry a third Fortran facade. Archive `swap_bmi_mod` as reference.

5. **Standalone executable unchanged** — keeps linking the kernel directly.

This is the "**commit to two**" option the brief invited: **XMI + pyswap-CAPI
over one shared instance core.** It removes the duplicate-symbol/two-library
problem, serves both real consumers, and preserves standalone.

**Viable alternative (secondary):** keep BMI for generic CSDMS couplers, but
only if a concrete third consumer appears — and still build it over the shared
instance core, never the singleton. Cost: a third facade to maintain.

**Hard prerequisite for all of it:** the handle core is unsafe until the shared
ambient services (`error_mod`, `swap_log`) are instance-scoped (Finding #1).
That ordering drives the roadmap.

---

### 3.4 Decision & in-memory-init prototype (updated 2026-06-02)

Dialogue refined §3.3 to a sharper, more standard-aligned conclusion, plus a
first implementation slice.

**Refined binding decision.** The BMI-vs-XMI framing was the wrong axis; the
real problem is a **split instance model** — `swap_bmi_mod` and `swap_capi_mod`
share one singleton (`capi_state`/`capi_config`), while `swap_xmi_mod` is the
only ensemble-backed path. Target end-state:

- **One `libswap`** over a handle-based instance core, exposing the **BMI method
  set + XMI extensions** (the `libmf6` pattern — XMI ⊇ BMI, so one library, not
  two colliding ones). Both consumers drive it via **`xmipy`**.
- **Keep a CAPI-flavored init/override surface** for pyswap's per-instance
  *parameterization* (config override at init) — the one thing the BMI/XMI
  runtime-`set_value` model doesn't cover.
- **Drop standalone Fortran BMI** (strictly dominated); regenerate CSDMS-BMI in
  Python over the C-API if ever needed. Standalone executable unaffected.
- **pyswap reuses `xmipy` (the wrapper), not `imod_coupler` (the app):** the
  coupler is a lockstep multi-model orchestrator; pyswap's many-*independent*-
  columns problem isn't a coupling problem. A Python glue layer
  (pyswap + flopy + imod_coupler) to make coupled SWAP–MODFLOW model-building
  smooth is a worthwhile, separate interop goal.

**Diskless-init invariant.** Standalone reads disk; Python feeds memory; both
converge on the *same* `swap_config_t` → identical seeding → identical run.
Nothing downstream of the config knows the difference (keeps physics
byte-identical).

**Reality found:** the existing in-memory `load_swap_config_from_string`
consumed only the TOML *body* from memory; companions were still read from disk:
load-time TOML subfiles (`.dra.toml`, `.crp.toml`) and seed-time CSV data
(meteo `283.csv` via `atmosphere_state%init`, and all other CSV tables via the
shared `read_csv_table`).

**Implemented (committed, TDD, byte-identical):**
- `config_source_t` (`src/io/toml/config_source.f90`) — a companion-content
  provider: disk-backed (standalone) or in-memory blobs (Python). Threaded as an
  argument, no module global. (commit `7a72fc7`)
- Threaded through the **load-time** companion readers (`read_drainage_toml`,
  `load_crop_rotation_files`) → an in-memory hupselbrook load (drainage + 3 crop
  subfiles as blobs) yields the same config as disk (`dramet=2`); standalone
  passes no source → disk fallback. (commit `58b8a0e`) — pFUnit 814, check-fast
  4/4 byte-identical.

**Remaining for a runnable Python-driven prototype:**
- **M2 — seed-time CSV from memory.** Add a from-text path at the shared
  `read_csv_table` layer (unlocks in-memory for *all* CSV companions), thread an
  optional `config_source_t` through `swap_init_from_loaded_config →
  state%atmosphere%init → meteo table load`. This is the keystone to a diskless
  *init* and touches a shared layer — do it as its own TDD slice.
- **M3 — surface + driver.** A CAPI push (`swap_attach_config_file(name,
  content, n)`) holding blobs on the singleton between calls; a Python `ctypes`
  driver that pushes TOML + companions, `initialize`, loops `update`, reads
  results (`swap_view_array`/water-balance) — verified against the disk run.

## 4. Cross-cutting anti-patterns

These recur across subsystems and are best fixed as themed arcs, not per-file.

**C1 — Shared mutable globals at the edges.** The complete list that blocks
multi-instance (kernel is clean): `error_mod::{global_errors, library_mode,
fatal_was_raised}`; `swap_log::{current_level, log_unit, log_to_file, …}`;
`swap_capi_mod::{capi_state, capi_config}`; `swap_ensemble_mod::{columns,
configs, column_config, ncol, gwl, qbot_volume, storage_coef}`;
`meteo_buffer_mod::{current_mode, ext_buffer, ext_n_days, ext_n_cols}`;
`csv_output::{profile_w, scalar_w, + ~10 module arrays}`. *Bounded and
enumerable* — that's the good news.

**C2 — Per-rotation `SAVE` in crop dispatchers.** `grass`/`wofost`/`cropfixed`
runtime `save`; `wofost_soil_declarations` `save`; `oxygenstress` module
`current_state` pointer. Unlike C1 these sit *inside compute* and carry physics
state across calls → the highest-risk parallelization items (latent ensemble
correctness bug).

**C3 — I/O inside compute kernels.** Blocking `read(*,*)` in `oxygenstress`
numerical routines (deadlocks headless), `wofost` `outbalcrop*` file writes
mid-task, crop debug writes (units 777/888). Distinct from benign in-memory
`afgen` table lookups (fine to keep).

**C4 — No in-memory results path.** Output hardwired to disk; the retired
`output_row` C-API left no programmatic accessor. Blocks (A).

**C5 — Hardcoded integer-switch dispatch.** `swbotb`, `croptype` (×6),
`swinter`, `swredu`, `swcalt`, `swfrost`, `SWBR`, stress models in
`rootextraction`. Extensibility tax; candidate for behavior-neutral
strategy/registry seams (must not perturb FP ordering).

**C6 — Oversized routines & files.** See Finding #9. Decompose into
orchestrator + narrative phase-helpers (per existing convention), not micro-fns.

**C7 — Fixed-size state arrays.** `2*MABBC`/`mairg`/`macp`-sized fields that
should be allocatable-to-actual. Memory bloat × N instances + stack risk.

**C8 — Dead / stranded code.** Two flavors: cleanly-quarantined `dormant/`
(fine, decide keep-vs-delete) vs **stranded half-ports** that will break the
build on the next refactor (`agetracer`) or are unreachable bodies
(`SWBR=1`, wofost dead branches). The latter are liabilities, not dormancy.

---

## 5. Per-subsystem findings

Compact tables; full evidence in the sweep transcript. Sev/Eff/Obj as defined
above. IDs preserved from the sweep for traceability.

### 5.1 bindings / driver (centerpiece — see §3)
| ID | Title | Cat | Sev | Eff | Obj |
|----|-------|-----|-----|-----|-----|
| bindings-01 | No instance-handle model (singletons everywhere) | no-instance-handle | High | M | parallel |
| bindings-02 | Empty parameter-setter API (`swap_set_scalar` rejects all) | incomplete | Med | S–M | extens |
| bindings-03 | CSV output file-only; no in-memory accessor | io-hardwired | Med | M | interop |
| bindings-08 | Redundant BMI/XMI/CAPI facades, duplicate symbols → 2 `.so` | redundant | Med | S–L | clarity |
| driver-01 | `swap_log` global state; columns share log unit/level | shared-global | Critical | M | parallel |
| driver-02 | `error_mod` global error collection + `library_mode` flag | shared-global | Critical | M | parallel |
| driver-03/08 | Ensemble columnar-isolated but **not** error/log-isolated | shared-global | High | M–L | parallel |
| driver-04/11 | `swap_run_step` monolithic — no clean one-substep seam for BMI | no-step-seam | High | M | extens |
| driver-05 | `external :: dtleap` / unmoduled `CropGrowth` interface | tangled | High | S | clarity |
| driver-06 | `timecontrol_mod` 693 L, implicit call-order contract | oversized | High | M | clarity |
| driver-09 | `tridag` 40 KB automatic `gamma(macp)` per call (stack-heavy) | other | Med | M | parallel |
| driver-10 | `arrayutils`/`interpol` call `fatalerr_collected` → global fatal | tangled | Med | M | parallel |

### 5.2 crop (fixed / grass / wofost / oxygen / cropgrowth)
| ID | Title | Cat | Sev | Eff | Obj |
|----|-------|-----|-----|-----|-----|
| cropfixedgrass-01 | ~30 local `SAVE` vars in `grass()` dispatcher | hidden-state | Critical | M | parallel |
| cropfixedgrass-02 | `SAVE` block in `cropfixed()` dispatcher | hidden-state | High | S | parallel |
| wofost-01 | `wofost_soil_declarations` module `save` bag (47 vars) | hidden-state | Critical | L | parallel |
| wofost-02 | `wofost_soil_interface` module vars (implicit `SAVE`) cross-module coupling | hidden-state | High | M | parallel |
| wofost-03 | `SAVE` in `wofost()` dispatcher (25+ vars across tasks) | hidden-state | High | M | parallel |
| cropshared-01/07 | `oxygenstress` module `current_state` pointer (non-reentrant) | hidden-state | Critical | M | parallel |
| wofost-04 | File I/O (`outbalcrop*`) + `SAVE` units inside task loop | io-in-kernel | High | M | parallel/interop |
| cropfixedgrass-04 | Debug `write` to units 777/888 in kernel | io-in-kernel | High | S | interop |
| cropshared-05 | Blocking `read(*,*)` in `QROMBD`/`POLINTD`/`ZBREND` | io-in-kernel | High | S | interop |
| cropfixedgrass-03 / wofost-06 / cropshared-09/11 | `croptype` & stress hardcoded dispatch (≥6 sites) | extens-seam | High | M | extens |
| cropshared-02/04 / wofost-05 / cropfixedgrass-05 | Oversized routines (OxygenStress 432 L, wofost task2/3 563 L, grass 1317 L) | oversized | Med | L | clarity |
| cropshared-08/15 / wofost-08/09/12 / cropfixedgrass-07 | Dead/stranded: jongvanlier, oxygenrepro, bulb/soybean/CO2 branches, `wofost_soil_parameters`, FillOxygenStress1/2 | dead/dormant | Med–Low | S–L | clarity |

### 5.3 soilwater (incl. boundary)
| ID | Title | Cat | Sev | Eff | Obj |
|----|-------|-----|-----|-----|-----|
| soilwater-04/10 | `swbotb`/top-BC hardcoded switch ladders (×2 modules) | extens-seam | High | M | extens |
| soilwater-01/14 | `headcalc` 573 L + inline Jacobian; hysteresis 109 L | oversized | High | L | clarity |
| soilwater-02 | `DATA` statement for convergence constants (subtle shared state) | non-reentrant | Med | S | parallel |
| soilwater-05/06/12 | Dormant `sptabulated` (5557 L), `watertable`, `checkmassbal`, `regrid` | dead/dormant | Med–Low | M | clarity |
| soilwater-11 | Per-node convergence diagnostics → log spam under N instances | clarity | Low | S | parallel |
| soilwater-09 | **Positive:** BC kernels are I/O-free (clean edge separation) | — | — | — | — |
| soilwater-13 | **Positive:** `WC_K_models` confirmed reentrant | — | — | — | — |

### 5.4 atmosphere
| ID | Title | Cat | Sev | Eff | Obj |
|----|-------|-----|-----|-----|-----|
| atmosphere-01/05 | Cumulative snow flux written directly into shared cohort state | hidden-state | Critical* | M | parallel |
| atmosphere-03/04/15 | Reference-ET / sub-daily loop write `state` mid-loop | non-reentrant | High | M–L | parallel |
| atmosphere-02/08 | `swinter`/`swredu` hardcoded dispatch | extens-seam | Med–High | M | extens |
| atmosphere-06/11 | `meteo_orchestrator` 624 L; tangled `partition_peva_ptra` | oversized | Med | L | clarity |
| atmosphere-10 | `ruttervw` real(4) MetaSWAP kernel + conversion buffers | other | Med | M | clarity |
| atmosphere-12 | Static rain-event arrays, no bounds check | other | Med | M | interop |

*\*atmosphere-01 severity assumes per-instance state is shared; since these
write to `state%atmosphere`, the risk is the same per-rotation pattern as crop —
confirm against the threaded-state model during planning.*

### 5.5 drainage / surfacewater
| ID | Title | Cat | Sev | Eff | Obj |
|----|-------|-----|-----|-----|-----|
| drainage-01 | Per-timestep `allocate` inside `drainage()` | non-reentrant | Critical | M | parallel |
| drainage-02 | `surf%flInitDraBas` one-time-init flag (hidden state) | hidden-state | Critical | M | parallel |
| drainage-04/05 | `DIVDRA` 27-arg signature + 337 L body | oversized | High | L | extens/clarity |
| drainage-03/06 | `dramet`/`ipos`/`swnrsrf` nested hardcoded dispatch | extens-seam | High | L | extens |
| drainage-08/09/11 | `WLEVBAL` outcome/target split still 74–110 L, mutate 5+ fields | tangled | High | M | clarity |
| drainage-07 | `afgen` table reads in kernel (benign) | io-in-kernel | Med | S | interop |

*Note: drainage-01/02 — verify whether allocation/flag are per-`state%drainage`
(then re-entrant across instances) or truly process-shared, during planning.*

### 5.6 solute
| ID | Title | Cat | Sev | Eff | Obj |
|----|-------|-----|-----|-----|-----|
| solute-02 | `SWBR=1` quarantined: fatal guard + unreachable body | dead | Critical | L | clarity |
| solute-05 | `dormant/agetracer` imports retired globals — **won't compile after next refactor** | dead/stranded | High | S (delete) / L (port) | extens |
| solute-03/04 | `update_solute_compartments` 97 L mixes 5 kernels | tangled | High | M | extens |
| solute-01 | `dtsolu` time-step state on record across substeps | hidden-state | Med | M | parallel |
| solute-06/10 | `SWBR`/Freundlich hardcoded, no strategy seam | extens-seam | Med–Low | L–M | extens |
| solute-08 | **Positive:** state threaded as single arg (well-architected) | — | — | — | — |

### 5.7 heat
| ID | Title | Cat | Sev | Eff | Obj |
|----|-------|-----|-----|-----|-----|
| heat-01/07 | `MACP`-fixed arrays in `devries`/`FrozenBounds`; uninit tails | non-reentrant | High | M | parallel |
| heat-06 | `FrozenBounds` mixes frozen-conductivity + drainage redistribution | tangled | High | L | clarity |
| heat-04/05 | `temperature_step` 149 L; duplicated tridiagonal stencil blocks | oversized | Med | L–M | clarity |
| heat-02/08/09/14 | `swcalt`/`swtopbhea`/`swbotbhea`/`swfrost` hardcoded dispatch | extens-seam | Med | M | extens |
| heat-03/10/11 | Magic snow constants; divide-by-zero risk; uninit `heaconBot` | other | Low–Med | S | clarity |

### 5.8 config
| ID | Title | Cat | Sev | Eff | Obj |
|----|-------|-----|-----|-----|-----|
| config-03 | `surface_water_config_finalize` mutates `alphaw` in place → not re-runnable for per-instance param override | mutation-in-build | High | M | parallel/interop |
| config-08 | No override interface for pyswap sweeps (re-parse per column) | extens-seam | High | L | interop/extens |
| config-04 | Parallel rotation arrays (8 of them) instead of one record array | extens-seam | High | L | extens |
| config-02 | `cropwofost_config` 1063 L, 26 types + 50 validators | oversized | High | L | clarity |
| config-01/05/07 | Hardcoded crop-type enum + duplicated `check_*` (101×) + cross-section rules in monolith | duplication | Med | M | clarity |
| config-09 | **Positive:** `irrigation_schedule_t` factoring is exemplary | — | — | — | — |

*config-03 nuance (verified): `finalize` runs once at parse to match legacy
normalization; the config is immutable afterwards. The real cost is that a
per-instance parameter override can't simply re-run `finalize` idempotently —
relevant to the pyswap override path, not a runtime mutation bug.*

### 5.9 state
| ID | Title | Cat | Sev | Eff | Obj |
|----|-------|-----|-----|-----|-----|
| state-06 | `soilwater_state` `2*MABBC` BC tables ≈ **4.7 MB/instance** | fixed-array | High | M | parallel |
| state-01 | `crop_irrigation_state` `mairg` arrays ≈ 520 KB | fixed-array | Critical | M | parallel |
| state-02 | `crop_oxygen_state` `macp` arrays ≈ 280 KB, no per-instance sizing | fixed-array | High | M | parallel |
| state-08 | `drainage_state::owltab` `(MADR, 2*MAOWL)` ≈ 300 KB | fixed-array | Med | M | parallel |
| state-09 | BMI accessor covers 8 of ~112 vars; no `get_value_ptr` registry | bmi-gap | Med | L | interop |
| state-10 | `swap_state_t` init order is implicit (scattered in `swap_mod`) | clarity | Low | M | clarity |
| state-04 | `crop_oxygen` `ini_stress` first-call init — verify per-instance isolation | hidden-state | High | M | parallel |

### 5.10 io
| ID | Title | Cat | Sev | Eff | Obj |
|----|-------|-----|-----|-----|-----|
| io-01/04 | Module-global `csv_writer_t` (`profile_w`, `scalar_w`, `save`) | shared-file-unit | Critical | M | parallel |
| io-03 | `meteo_buffer_mod` module-global mode/buffer | shared-state | Critical | M | parallel/interop |
| io-06/10 | No in-memory results buffer; headless flag suppresses, doesn't collect | io-hardwired | High | L | interop |
| io-05 | Output filename from single `pathwork/outfil` → collide across instances | io-hardwired | High | S | interop |
| io-08/09 | `csv_output` 1291 L; new var needs edits in 3 places | oversized/seam | Med | L–M | clarity/extens |
| io-02/11 | `SAVE` locals spanning header→write; `ReadMeteoYear` first-call alloc | hidden-state | Med–Low | S–M | parallel |

---

## 6. Dormant / dead-code inventory

| Item | LOC | State | Recommendation |
|------|-----|-------|----------------|
| `soilwater/dormant/sptabulated.f90` | 5557 | Quarantined, no live caller (swsophy=1 dark) | **Archive or delete** — embeds a 3370-L TSPACK spline lib; biggest false-complexity source. Decide explicitly. |
| `solute/dormant/agetracer.f90` | ~400 | **Stranded half-port** — `use variables` of retired globals | **Delete or complete now** — will break the build on next refactor. Highest-priority dead item. |
| `crop/dormant/jongvanlier.f90` | — | Dormant (swdrought=2); `alpJvLier` ref dangling in `rootextraction:188` | Remove dangling ref; keep file dormant or delete. |
| `crop/dormant/oxygenrepro.f90` | — | Dormant (swoxygentype=2), fatal-stubbed dispatch | Keep dormant w/ checklist, or delete. |
| `soilwater/dormant/{watertable,checkmassbal,regrid}.f90` | ~29K | Cleanly quarantined w/ READMEs | Leave; annual review. |
| `solute SWBR=1` body | ~40 | Unreachable (fatal guard) | Delete body + guard; issue-track the fix. |
| wofost bulb/soybean/CO2 branches | ~600 | Dead on TOML path (guarded) | Delete; validator prevents reactivation. |
| `wofost_soil_parameters.f90` | 179 | Never called; refs retired `BDENS` | Delete. |
| `FillOxygenStress1/2` | ~200 | Called, but bodies fully commented | Delete subs + calls. |

---

## 7. Prioritized roadmap

Dependency-ordered. Each arc is a separate brainstorm → plan → implement cycle,
ending `check-fast` green, byte-identical. Arcs 3–7 are largely independent of
each other once Arc 1 lands.

**Arc 1 — Instance-safe ambient services** *(foundational; unblocks all
parallelization).* Move `error_mod` (`global_errors`, `library_mode`,
`fatal_was_raised`) and `swap_log` state off module-global onto a per-instance
context (thread via `state`, or an injected handle). Shared by every subsystem →
must come first. *Findings #1, driver-01/02/03/08/10.* **Effort M.**

**Arc 2 — Handle-based instance core + binding consolidation** *(the
multi-instance enabler; §3).* Retire `capi_state`/`capi_config` and the
`swap_ensemble_mod` singletons into an opaque-handle instance registry. Reshape
XMI as a thin facade over it; stand up the pyswap C-API skeleton; archive
Fortran BMI. *Findings #2, #5, bindings-01/08, driver-03.* **Effort L.** Depends
on Arc 1.

**Arc 3 — Crop hidden-state → `state%crop%*`** *(closes the latent ensemble
correctness bug).* Migrate `grass`/`wofost`/`cropfixed` runtime `SAVE`,
`wofost_soil_declarations`, and the `oxygenstress` `current_state` pointer onto
threaded state. Highest-risk because it's *in* the physics — guard with
byte-identical regression at every step. *Findings #3, C2, crop hidden-state
rows.* **Effort L.** Independent of Arcs 4–7.

**Arc 4 — In-memory results + I/O de-globalization** *(unblocks pyswap (A)).*
Instance-own `csv_output` writers + `meteo_buffer_mod`; add an in-memory results
accessor (water balance + profile arrays) and per-instance output namespacing;
remove in-kernel I/O (wofost `outbalcrop*`, oxygenstress blocking `read`, crop
debug units). *Findings #4, #6, C3/C4, io-01/03/05/06.* **Effort L.** Needs Arc 1.

**Arc 5 — State array right-sizing** *(N-instance memory scaling).* `2*MABBC`,
`mairg`, `macp`, `MAOWL` fixed fields → allocatable-to-actual. *Findings #7, C7,
state-01/02/06/08; heat-01/07; drainage-01.* **Effort M.** Independent.

**Arc 6 — Decomposition + dispatch seams** *(clarity + extensibility; ongoing,
parallelizable across subsystems).* Break up oversized routines
(`OxygenStress`, wofost task2/3, `grass`, `divdra`, `headcalc`, `WLEVBAL`,
`cropwofost_config`) into orchestrator + narrative helpers; extract
behavior-neutral strategy/registry seams for `swbotb`/`croptype`/stress/etc.
Add a clean one-substep advance seam for BMI/coupling (driver-04). *Findings #8,
#9, C5/C6.* **Effort L, incremental.** Independent.

**Arc 7 — Dead/stranded-code purge** *(do `agetracer` + `SWBR=1` early).* Delete
stranded half-ports and unreachable bodies; decide keep-vs-archive for the
`dormant/` set (§6). *Finding #10, C8.* **Effort S–M.** Independent; the
`agetracer` deletion should happen before Arc 1's refactor touches globals.

**Suggested sequence:** Arc 7 (agetracer/SWBR quick wins) → **Arc 1** → **Arc 2**
→ then Arcs 3/4/5 in parallel as capacity allows, with Arc 6 running
continuously as subsystems are touched.

---

## Verification

The following load-bearing claims were checked directly against source (not
taken from agent reports):

- `swap_capi_mod:24-25` — `capi_state`/`capi_config` `save` singletons ✓
- `swap_ensemble_mod:24-31` — `columns/configs/column_config/ncol` +
  `gwl/qbot_volume/storage_coef` all `save` ✓
- `error.f90:57,62,63` — `global_errors`, `library_mode`, `fatal_was_raised`
  `save` ✓
- `swap_log.f90:38-44` — module-level mutable scalars (implicit `SAVE`) ✓
- `surface_water_config.f90:170-191` — `finalize` mutates `alphaw` in place ✓
- `cropgrass_runtime:93`, `cropwofost_runtime:163`, `wofost_soil_declarations:119`
  — explicit `save` ✓; `wofost_soil_interface:10` — explicit `save` *removed* but
  module vars remain implicitly `SAVE` ✓
- `oxygenstress:76,145` — module `current_state` pointer set per call ✓
- `meteo_buffer_mod:30-33`, `csv_output:17,427` — module-global state ✓
- `tests/coupling/hupselbrook_coupled/ensemble.txt` + single shared config —
  confirms the crop-`SAVE` correctness risk is currently masked by identical
  columns ✓

Line numbers in the per-subsystem tables (§5) are from the sweep agents and are
reliable for the structural claims spot-checked above; treat individual
unchecked line numbers as ±a few lines, to be reconfirmed at planning time for
each arc.
