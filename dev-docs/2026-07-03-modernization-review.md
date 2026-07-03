---
title: SWAP modernization — state-of-the-code review and multi-arc sprint map
date: 2026-07-03
status: review
tags: [review, refactoring, isolation, efficiency, coupling, multi-instance]
---

# SWAP modernization — state-of-the-code review

**Purpose.** The strangler-fig migration is complete (`TOML → swap_config_t →
swap_state_t`, no bare globals; see
[`post-phase-4-modernization-summary.md`](post-phase-4-modernization-summary.md)).
This document surveys what is *left behind* after that work and maps it to a
spec-ready, multi-arc sprint. The stated objective for the next round is a SWAP
that is **modern, modular, extensible, efficient (incl. parallel), ensemble-capable,
and coupling-ready (MODFLOW 6)**. Every finding below is scored against that.

**Method.** Five parallel deep-dive passes (state/data-exchange, duplication/dead
code, I/O layer, Python orchestration/coupling, build/test infra) plus a hands-on
profiling pass (gprof + strace on the `hupselbrook` regression case, `~1.1 s`
single column). Claims that anchor recommendations were spot-verified against the
source, not taken from agent summaries. File:line references are current as of
this date.

**How to read it.** Findings are grouped by theme, each with an impact/effort/
payoff line. Two explicitly separated tiers close the document:

- **Tier 1 — incremental cleanup** within today's `swap_state_t`/`swap_config_t`
  architecture (dead code, dedup, init consolidation, efficiency, the last
  multi-instance blockers, CI).
- **Tier 2 — component isolation**, a longer-horizon, MODFLOW 6-inspired
  restructure (explicit exchange objects, per-component boundaries, handle-based
  multi-instance, a data-driven variable registry).

A headline summary sits up top; the arc breakdown sits at the bottom.

---

## 0. Headline

The engine is in good shape. The strangler migration delivered: typed state
threaded by argument, I/O confined to `src/io/` (verified — the hot loop issues
**zero** file operations; a full run makes 38 `write` syscalls total, all
buffered), a clean 98-suite pFUnit harness, and a byte-identical regression
oracle (`swap420gf`, 18/20 cases passing). Single-column runtime is already fast
and allocation-clean.

What remains is **clarity and capability**, in four buckets:

1. **Component isolation is incomplete.** Compute kernels still reach across
   compartment boundaries by direct field access — soil-water compute *writes*
   five atmosphere fields; crop and solute *read* heat/atmosphere state. There is
   no exchange-boundary object; everything shares one mutable `swap_state_t`. This
   is the gap between "typed state" and "coupling-ready components."

2. **A handful of module-global blockers still prevent true (threaded)
   parallelism** — the Tier-3 crop cluster (`crop_config_global`, `SAVE` locals)
   and the `Wofost_Soil_Interface` / `Wofost_Soil_Declarations` implicitly-static
   module variables shared across all crop instances.

3. **Dead weight and boilerplate**: ~9 K lines of dormant/dead code, a hand-listed
   `meson.build`, per-field BMI/XMI accessor switch-statements (6 edit sites per
   new exposed variable), and a genuinely legacy `MACP = 5000` fixed dimension
   used for per-step scratch arrays.

4. **The Python orchestration/coupling layer is a working prototype, not a
   product** — Phase A shipped but is stuck in `prototype/`; the MODFLOW 6 demo
   works but its storage exchange is knowingly wrong; three separate Python
   drivers each re-implement library loading.

None of this is a regression risk to fix — most is behavior-preserving cleanup —
but the ordering matters, because isolation (bucket 1) and de-globaling (bucket 2)
are prerequisites for the coupling/ensemble objective.

---

## 1. Component isolation & data exchange  *(the core ask)*

Today `swap_state_t` is a flat aggregate of 13 subsystem records
(`src/state/swap_state.f90:28`). Any compute routine can read or write any other
subsystem's fields by name. That was the right end state for the *migration* (kill
the globals first), but it is not component isolation — it is globals with a
`state%` prefix. Concrete breaks:

| # | Break | Location | Direction |
|---|-------|----------|-----------|
| 1 | Soil-water compute **writes** 5 atmosphere fields (`nraidt`, `nird`, `ldwet`, `spev`, `saev`) | `src/soilwater/soilhydraulics.f90:800-805` | soilwater → atmosphere (write) |
| 2 | Soil-water compute reads heat (`rfcp`, `tsoil`), drainage (`qdra`, `nrlevs`), atmosphere (`nraidt`, `nird`, `melt`) | `soilhydraulics.f90:79-85, 102-117` | soilwater ← 3 subsystems |
| 3 | Root extraction reads `heat%tsoil` directly | `src/crop/rootextraction.f90:258` | crop ← heat |
| 4 | Solute transport associates `atmosphere` + `heat` alongside soil | `src/solute/solute.f90:39-42` | solute ← atmosphere, heat |
| 5 | WOFOST↔soil nutrient exchange via 13 **module-global** variables | `src/crop/wofost/wofost_soil_interface.f90:16-22` | hidden static coupling |

Break #1 is the sharpest: a soil-water routine mutating atmosphere state means the
atmosphere subsystem cannot be tested, replaced, or reasoned about in isolation,
and the write ordering is an implicit contract. Break #5 is the worst for the
objective: those 13 `real(8)` module variables are implicitly `SAVE` (the file
comment even notes the explicit `save` was removed because "module variables are
implicitly static") — so **every column in an ensemble shares them**. That is a
data race the moment columns run concurrently.

**MODFLOW 6 comparison (Tier 2 direction).** MF6 separates *models* (own their
state), *exchanges* (own the coupled variables that cross a boundary and the logic
to move them), and a *memory manager* (a registry that owns and names every
allocatable, so nothing is reached by raw `use`). SWAP has models (the subsystem
records) but no exchanges and no registry. The Tier-2 proposal is not "adopt MF6
wholesale" — it is to introduce **typed exchange records** for the handful of real
inter-compartment surfaces (soilwater↔atmosphere flux, soilwater↔heat, crop↔heat,
crop↔soil-nutrients) that a compute routine returns/consumes explicitly, instead
of reaching into a sibling record. This is standard practice for coupled
process models of this class and is the single highest-leverage structural change
toward the coupling objective.

- **Impact:** high (foundation for isolation, testability, and safe concurrency).
- **Effort:** medium per surface (there are ~4–5 real surfaces, not dozens).
- **Payoff:** high — each closed surface is a subsystem that becomes independently
  testable and thread-safe.

---

## 2. Multi-instance / parallelism blockers  *(objective-critical)*

The orchestration plan ([`2026-06-11-python-orchestration-plan.md`](2026-06-11-python-orchestration-plan.md))
correctly identifies these as the gate to threaded parallelism (its Phase C). They
are unaddressed. In priority order:

1. **`Wofost_Soil_Interface` module variables** (`wofost_soil_interface.f90:16-22`,
   consumed in `management_soil.f90:144,324`) — 13 shared static scalars. **Verified
   still present.**
2. **`Wofost_Soil_Declarations` blanket-static nutrient module** (~100 vars) — the
   ADR 0047 W10 arc retired the state *writes* into globals, but the module-static
   working set remains a per-process singleton.
3. **`crop_config_global` + `SAVE` locals** in `cropgrowth.f90:601,678` and
   `cropgrass_runtime.f90:94` — the ADR 0043 Tier-3 forward note; the literal
   blocker to running two columns in one process.
4. **Binding singletons** `capi_state` / `capi_config` (`swap_capi_mod.f90:28-34`),
   shared into `swap_bmi_mod`, and the ensemble module arrays
   (`swap_ensemble_mod.f90:27-35`). These are *by design* for the current
   one-process-one-column model; true multi-instance needs a **handle-based** C-ABI
   (opaque instance pointer) — a Tier-2 item, but the Fortran side is already
   `state`-threaded, so the handle is mostly a binding-layer change.

Process isolation (the current `multiprocessing`-spawn orchestrator) sidesteps 1–3
today, which is why the prototype works. But the XMI ensemble loop is **sequential**
even in one process (`swap_ensemble_mod.f90:128-136`); making it `!$omp parallel do`
— the whole point of the in-process ensemble — requires 1–3 retired first.

- **Impact:** high (blocks threaded ensemble and gridded in-process coupling).
- **Effort:** medium; each step is byte-identical-gated and verifiable with a
  2-column "same inputs ⇒ identical outputs, no cross-talk" test.
- **Payoff:** high — unlocks Phase C of the orchestration plan.

---

## 3. Efficiency  *(profiled)*

**Setup.** gprof (`-pg`) and `strace -c` on `hupselbrook` (`~1.1 s`, 1461 days,
~29 K Richards sub-steps). `perf` is unavailable in this sandbox
(`perf_event_paranoid=4`).

**What is already good.** Output is fully buffered (38 `write` syscalls for the
whole run); no per-timestep file opens; no syscall-level allocation churn (7 `mmap`,
8 `brk`). The I/O-at-the-edges design holds up under measurement. Single-column
wall time is not a problem.

**The `MACP = 5000` scratch-array tax.** `MACP` (`src/core/arrays.f90:24`) is a
legacy maximum-compartment dimension. Several per-step routines declare **local
scratch arrays sized to `MACP`**, not to the actual node count (typically 40–100):

- `src/solute/solute.f90:37` — `thetav, dispr1, vpore2` (3 × 40 KB), zeroed every
  solute step; `:128` — `diffus`.
- `src/core/numericalsolvers.f90:52` — `gamma(macp)` in `tridag`, zeroed on all
  ~88 K calls.
- `src/soilwater/soilhydraulics.f90:1320-1322` (`hysteresis`),
  `src/heat/frozencond.f90:200-201` — conditional paths.

With the project-wide `-finit-local-zero` flag, each such array is memset-zeroed on
every entry. Under the distorting `-pg` build this shows as ~46 % "memset"; on the
real `-O2` binary, sizing these to `numnod` yields only **~4 %** single-run
improvement (the optimizer hoists/elides much of it). **So this is not a
single-column speed problem.** It matters for two forward reasons:

1. **Threaded ensemble scaling.** Under `!$omp parallel do` over columns, each
   thread gets its own copy of these stack arrays — cache footprint and zeroing
   cost multiply by thread count, exactly where the ensemble objective lives.
2. **Cheap, byte-identical hygiene.** Sizing scratch arrays to the actual node
   count (an automatic array bound like `state%mesh%numnod`) is a mechanical,
   regression-safe change that removes a latent scaling hazard.

**`-finit-local-zero` is load-bearing for correctness.** Rebuilding without it
**hangs** the model (non-convergence) — meaning genuine uninitialized-variable reads
are masked by the flag. This is both a latent-correctness and an efficiency debt:
the zeroing tax cannot be dropped wholesale until the real uninitialized reads are
found (via `-Wuninitialized`/`-Wmaybe-uninitialized` triage) and fixed. Worth a
scoped arc; do **not** simply remove the flag.

**Van Genuchten `pow`/`exp` is the compute floor.** `watcon`/`moiscap`/`hconduc`
(`src/soilwater/soilhydraulicsutils.f90`) are each called ~2 M times per run and
dominate genuine FLOPs via `pow`/`exp`. The dormant
`src/soilwater/dormant/sptabulated.f90` (5,557 lines, Hermite-spline tabulation of
exactly these curves) appears to be a purpose-built tabulation path — the real
performance lever, but a **byte-identity risk** (tabulation ≠ analytic to the last
ULP). Treat as a separate, opt-in, well-guarded experiment, not a cleanup.

- **Impact:** medium now, high for threaded ensembles.
- **Effort:** low (scratch sizing); medium (uninit triage); high/risky (tabulation).
- **Payoff:** scratch sizing = cheap parallel-readiness; uninit triage = unblocks
  dropping the zeroing tax.

---

## 4. Dead code & hygiene  *(DONE — arc T1-B, ADR 0049)*

> **Updated 2026-07-03 after execution.** The original estimate here ("~9 K lines
> of dead code") was wrong on two counts, and the correction is the finding.
> Delivered as arc T1-B; details in
> [`2026-07-03-arc-t1b-deadcode-sweep.md`](2026-07-03-arc-t1b-deadcode-sweep.md)
> and ADR 0049.

- **The `dormant/` modules (~9 K lines) are an intentional museum, NOT dead
  code.** All three `dormant/` dirs have a `README.md` with per-module
  reactivation checklists; meson tombstones point at them; they are build-excluded
  (zero binary cost). Policy formalized in **ADR 0049**: deletion requires an ADR.
  Kept. (The READMEs were refreshed — they had drifted, still listing the
  ADR-0047-deleted `jongvanlier.f90` and `variables.f90`/`state%cfg` restoration
  steps.)
- **`stepnr()`** (`arrayutils.f90`) — genuinely dead (self-documented "not called
  from any SWAP code", zero callers). **Removed.**
- **State fields: 10 of the 12 candidates are LIVE** (verified by full-tree
  read/write grep) — `evp`, `runonarr`/`flrunon`, `alfaw_layer`,
  `ksatexm`/`cofani`/`flksatexm`, `bpegwl`, `psilt`/`pclay`. Kept; recorded in
  ADR 0049 so they are not re-attempted. Only **`pegwl` and `npegwl`** were
  genuinely dead (write-only in `waterbalance.f90`) — **removed**, byte-identical.
- **Macropore "retired-zero placeholders"** — post-strangler these are mostly
  tombstone comments plus a couple of load-bearing always-zero locals; `FrArMtrx`
  is live (read in output). Nothing safely removable. **Out of scope.**

- **Outcome:** small LoC delta (1 function + 2 fields), high correctness value (no
  false-dead deletions), museum policy + docs now durable. `check-fast` 4/4
  byte-identical, pFUnit 833 green.

---

## 5. Duplication → shared utilities  *(medium, do selectively)*

The two duplication passes partly disagreed; reconciled view:

- **TOML readers (18 files, ~2,710 lines).** The `toml_field_helpers` /
  `toml_array_helpers` abstraction **is** used consistently (this is not
  copy-paste); what remains is that each reader hand-lists every field as a
  sequential `get_optional_*` call. A **declarative field registry** (table of
  section/key/type/default) + one generic loop could collapse a large fraction, but
  the current state is maintainable and low-risk. **Recommendation:** treat as
  optional, not urgent; pursue only if the config schema keeps growing.
- **Crop `copy_pair_table`/`copy_table`** duplicated across `cropfixed_init.f90:172`,
  `cropgrass_init.f90:456`, and inline in `cropwofost_init.f90` — extract one
  `arrayutils` helper (~20 lines saved, removes coupling).
- **`clamp_01` / `max(0,min(1,·))`** repeated 50+ times — cosmetic; a one-liner
  helper is nice-to-have, low payoff.
- **CSV table loaders** (`src/io/csv/*.f90`) are correctly *parallel-but-distinct*
  (per-family parsing) — leave alone.

- **Impact:** medium (readers) to low (rest). **Effort:** medium. **Payoff:**
  selective.

---

## 6. Two-phase init consolidation  *(medium)*

`swap_init_body` (`src/driver/swap_mod.f90:53-202`) mixes Phase-1 type-bound
`state%X%init()` with Phase-2 free `*_seed()` procedures, and the ordering encodes
implicit cross-subsystem dependencies (crop before timecontrol; `soilwater_seed`
after `heat%init`; `temperature_seed` reads soilwater). Per ADR 0043 this split is
deliberate and correct, but the *contract is opaque* — a new subsystem author must
reverse-engineer the ordering. Init signatures are also non-uniform (`crop%init`
takes 3 args, `soilwater%init` takes 7). This is the natural companion to §1: once
inter-compartment data moves through explicit exchange records, the Phase-2 seed
inputs become explicit and the ordering becomes checkable.

- **Impact:** medium (clarity, safe extension). **Effort:** medium. **Payoff:**
  medium; compounds with §1.

---

## 7. BMI/XMI binding boilerplate  *(medium)*

`src/bindings/{swap_bmi,swap_capi,swap_xmi}_mod.f90` expose variables via
hand-maintained `select case` switches: `bmi_get_value_double` (8 cases),
`bmi_set_value_double` (2), plus separate switches for names, units, and nbytes.
Adding one exposed variable touches **~6 sites**, with no compile-time check that
the name lists agree. `c_to_f_string`/`f_to_c_string` are defined **3×**. This is
the binding-layer face of "no registry" (§1's Tier-2 memory-manager point): a
single **data-driven variable table** (name, rank, units, a `state` accessor) would
drive get/set/metadata from one source and make the exposed surface
introspectable — directly useful for coupling, where the exchange var set is the
contract.

- **Impact:** medium (extensibility, coupling ergonomics). **Effort:** medium.
  **Payoff:** medium; enabling for Tier 2.

---

## 8. Python orchestration & MODFLOW 6 coupling  *(product gap)*

Verified against `prototype/`, `tests/coupling/`, and the orchestration plan:

- **Phase A shipped and works** (`prototype/swap_orchestrator.py`): netCDF-driven,
  process-pool parallel, BMI in-memory, override-injectable, balance closes,
  netCDF round-trips. But it is **stuck in `prototype/`** — no package, no CLI, no
  pyswap integration (plan Phases A-tail, E).
- **Coupling demo works but is qualitative.** `tests/coupling/run_coupled.py` drives
  a 2-channel MF6 model + 1 SWAP ensemble via XMI (MODFLOW leads the clock). The
  **storage exchange is currently disabled with cause** (`EXCHANGE_STORAGE=False`,
  `swapmod.py`): SWAP emits a constant `0.15` placeholder and MF6 `STO/SS` is
  *specific storage* (~1e-5), so writing it there is ~1e4× too large and pins the
  water table (observed). The correct fix is a physical phreatic Sy routed to
  `STO/SY` — a physics-design task (T1-I′), not a bug-fix. Recharge (the load-
  bearing exchange) is live and validated qualitatively.
- **Three Python drivers duplicate library loading and result extraction** (BMI
  orchestrator, XMI runner, coupled `SwapMod`). No unified "run a column/ensemble"
  entry point.
- **`imod_coupler_fork`** is a thin, self-contained vendored adapter (pydantic v2,
  stdlib logging) — fine for a demo, but needs governance (version pinning, CI
  against new `libmf6`) if it becomes the production path.
- **Not started:** Phase B (PEST/SALib front-ends), Phase C (Tier-3 globals →
  threading — see §2), Phase D (first-class gridded ensemble C-ABI so Python drives
  heterogeneous grids directly, not only XMI's homogeneous `ensemble.txt` path).

- **Impact:** high (this *is* the objective). **Effort:** mixed. **Payoff:** high.

---

## 9. Build & test infrastructure  *(mostly healthy; one quick win)*

- **No CI gate on PR/push** — `.github/workflows/ci.yaml` runs on version tags only;
  `check-fast`/`test-pfunit` run **only locally**. A breaking change can merge. A
  `test.yaml` running `check-fast` + `test-pfunit` on push/PR is the **single
  cheapest high-value fix in this document** (~2 h). **Do first.**
- **`meson.build` is a hand-listed monolith** (134 legacy + 6 modern sources in the
  root file). Works, no stale entries, but doesn't scale; `subdir()` per feature
  folder is a maintainability follow-up.
- **Physics kernels lack unit tests** — drainage, solute transport, soil-hydraulics
  models, Richards convergence are covered only by regression (byte-identical
  oracle), not algorithmically. Acceptable for a rescue, limiting for isolated
  debugging. Pairs naturally with §1 (isolated components are unit-testable).
- **Open regression items:** `winter` (frost-path xfail, ~0.5–0.8 cm, root cause
  suspected heat↔frost feedback), `swdrought2` (pending restore, perf-blocked),
  `swbotb=7` (uninvestigated modern-only crash — sentinel `999` fed into
  runoff/qbottom arithmetic), `swsalinity=2` (no oracle — `swap420gf` itself
  SIGSEGVs). All documented in `tests/regression/INVESTIGATION_NOTES.md`.

- **Impact:** high (CI) to medium (rest). **Effort:** low (CI). **Payoff:** high (CI).

---

## 10. Proposed sprint map

Two tiers, explicitly separated. Tier 1 is behavior-preserving and byte-identity-
gated. Tier 2 is the structural direction Tier 1 sets up.

> **Rescoped 2026-07-03** after landing T1-A/B and T1-G′ (the binding-convergence
> arc, ADR 0050). Key movements: **T1-G expanded into T1-G′** and delivered *one*
> `libswap.so` over the ensemble backing store with a variable registry — which
> **absorbs most of Tier-2's T2-D and half of T2-B**. **T1-I is re-scoped** (the
> coupling storage exchange is not "actively wrong"; it is *disabled with cause* —
> `EXCHANGE_STORAGE=False` in `swapmod.py` — so the work is a physics-design task,
> not a bug-fix). **T2-E is promoted earlier** (gridded ensemble C-ABI is now a
> natural extension of the unified ensemble surface). Two small infra items and a
> process rule were added from execution experience.
>
> **Reframed again by ADR 0051 + 0052.** A strategic boundary was drawn (ADR 0051):
> SWAP-native = vadose-zone water/heat/solute + the native table-driven "Plant";
> detailed crop growth (WOFOST/grass) and nutrients are *coupled components*, not
> embedded models. First detachment landed (ADR 0052): the WOFOST-N / ANIMO-derived
> soil-N subsystem was **deleted** (~2.7k LoC, ~85 module globals). So the §1
> crop↔soil-nutrient isolation break (#5) and the §2 nutrient blockers are **gone**,
> and **T1-E-b is moot**. Detailed nutrients are now external ANIMO (SWAP hydrology
> via the parked `.afo` output). The next structural item is T2-A's crop↔water
> exchange record — the seam an external WOFOST would plug into.

### Method rule (added from execution)

**Verify load-bearing claims from source before writing a spec.** Both executed
arcs caught material errors in agent/document claims only via code verification
(T1-B: the "9 K dead lines" were a deliberate museum + 10/12 "dead" fields were
live; T1-G′: the consumer map wrongly claimed the coupled driver calls
`get_time_step` on SWAP and missed that `report_timing_totals` never enters the
library). Every arc's discovery ends with this check.

### Tier 1 — incremental cleanup (within current architecture)

| Arc | Scope | §  | Effort | Payoff | Status / order |
|-----|-------|----|--------|--------|----------------|
| **T1-A CI gate** | `check-fast` + pFUnit on push/PR (`test.yaml`) | 9 | XS | High | **DONE** (unpushed) |
| **T1-B Dead-code sweep** | `stepnr`, write-only `pegwl`/`npegwl`; dormant-museum policy (ADR 0049); README refresh | 4 | S | Med | **DONE** |
| **T1-G′ Binding convergence** | one `libswap.so` (BMI+CAPI+XMI) over the ensemble; variable registry; C-string dedup (ADR 0050) | 7,2 | M–H | High | **DONE** |
| **T0.5 Binding-gate infra** | a `check-bindings` pixi task (the 6-path gate); re-point meson `bmi`/`cffi-demo` suites off the un-checked-out `tests/swap-cases` onto `tests/regression/cases` | 9 | XS | High | **next** |
| **T1-E Crop reentrancy** *(program, not one arc — see [arc doc](2026-07-03-arc-t1e-crop-reentrancy.md))* | `crop_config_global` already retired. **E1 (DONE):** TOTASS/ASSIM de-save. **E2:** `wofost()` ≥10 carries-state → `crop%wofost`. **E3:** `grass()` carries-state → `crop%grass`. **E4:** `cropfixed()`. **E5:** `O2_pars` state-pointer. **~~E-b~~ (DONE — deleted, not de-globalled):** the ~85 WOFOST nutrient globals were removed with the nutrient subsystem (ADR 0052), so E-b is moot | 2 | E1 S; E2–E5 M each | High | **E1 done; nutrients detached; E2 next** |
| **T1-D Uninit-read triage** | `-Wmaybe-uninitialized` audit; fix real reads behind `-finit-local-zero` | 3 | M | Med | mid |
| **T1-C Scratch-array sizing** | `MACP` locals → `numnod` (solute `:37,:128`, tridag `:52`, …), byte-identical | 3 | S | Med (parallel) | filler |
| **T1-F Init consolidation** | uniform `state%X%init`, explicit Phase-2 seed inputs | 6 | M | Med | mid |
| **T1-H′ Orchestrator → package** | promote `prototype/` to package + CLI over the single `libswap.so`; unify the 3 drivers | 8 | M | High | independent |
| **T1-I′ Phreatic Sy design** | derive a physical specific-yield from SWAP soil hydraulics, route to MF6 `STO/SY`, re-enable `EXCHANGE_STORAGE`; numeric validation case | 8 | M | High | with T2-E |
| **T1-J Small dedup** | crop `copy_table` helper; optional TOML field registry; XMI `get_value_ptr` → `NS_XMI` registry | 5,7 | S | Low | opportunistic |

### Tier 2 — component isolation (longer horizon, MODFLOW 6-inspired)

| Arc | Scope | §  | Effort | Payoff | Status |
|-----|-------|----|--------|--------|--------|
| **T2-A Exchange records** | typed exchange objects for soilwater↔atmosphere (start with the 5-field write, `soilhydraulics.f90:800-805`), soilwater↔heat, crop↔heat, crop↔soil-N | 1 | H | High | structural keystone |
| **T2-B Handle-based multi-instance** | per-instance ensemble + error state (opaque handle). Backing-store half **already done** by T1-G′ (`capi_state` is a pointer into `columns(1)`); needs T1-E first | 2 | M | High | after T1-E |
| **T2-C Threaded ensemble** | `!$omp parallel do` over columns in `swap_ensemble_mod` (needs T1-E + T2-B) | 2,8 | M | High | after T1-E/T2-B |
| **T2-D Variable registry** | **DONE via T1-G′** (`swap_var_registry_mod`); only the `NS_XMI` rider remains (folded into T1-J) | 1,7 | — | — | **DONE** |
| **T2-E Gridded ensemble C-ABI** | `bind(C)` for heterogeneous per-column config/meteo — now a natural extension of the unified `libswap.so` ensemble surface | 8 | M | High | **promoted** (after T1-E) |
| **T2-F Kernel unit tests** | algorithmic tests for now-isolated components | 9 | M | Med | after T2-A |

**Next arcs, in order:** **T0.5** (binding-gate infra, hours) → **T1-E** (Tier-3
de-global — the keystone that unblocks T2-B, T2-C, and real gridded coupling) →
**T2-E + T1-I′** (gridded ensemble surface + phreatic Sy) on the ensemble. T1-H′ is
independent and can interleave. **T2-A** (exchange records) is specced once T1-E
lands. (Push the 8 local commits and watch the first `test.yaml` run before
continuing — tests-only on push.)

---

## Appendix — verification notes

- Isolation break #1, `stepnr` deadness, `Wofost_Soil_Interface` module vars, and
  dormant-file build exclusion were spot-verified in source, not taken from agent
  summaries (per the "verify state, not reports" rule).
- Efficiency figures are from a local `hupselbrook` profile; the `-pg` "memset 46 %"
  is an instrumentation artifact (disabled inlining) — the honest single-column
  gain from the scratch-array fix is ~4 %, and the finding's value is
  parallel-scaling + hygiene, not single-run speed. Stated as such above.
- No source was changed by this review; all experimental edits were reverted and
  `git status` is clean.
</content>
