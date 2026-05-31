---
title: Post-Phase-4 modernization — capstone summary
date: 2026-05-31
status: complete
tags: [state-migration, globals-retirement, strangler-pattern, byte-identical-regression]
---

# Post-Phase-4 modernization — capstone summary

**Window:** 2026-05-12 → 2026-05-31 (~20 days)
**Scope:** ~851 commits on `development`; net `src/` delta ≈ **−5,500 lines**
(≈61,500 added / ≈67,000 deleted).
**Companion to:** [`phase-4-modernization-summary.md`](phase-4-modernization-summary.md)
(input-pipeline modernization, 2026-04-22 → 2026-05-05).

## Abstract

Phase 4 modernized SWAP's *input* pipeline: a TTutil fixed-format reader
chain (`.swp`/`.dra`/`.crp`/`.met`/`.YYY`) was replaced by a typed TOML
pipeline (`load_swap_config` → `swap_config_t` → a `config_to_variables`
adapter → legacy module-level globals). It left intact the runtime's
single largest piece of legacy architecture: `src/core/variables.f90`
— a ~1,500-declaration `save`d module of bare globals that every compute
and output routine read and mutated via `use variables`. This implicit
shared-mutable surface blocked multi-instance execution, made Python
embedding impossible against a clean state object, and coupled every
subsystem through name-level globals.

The post-Phase-4 work completed the **strangler-fig** transformation:
**TOML → typed config → typed state, with no bare-globals layer.** It
proceeded subsystem by subsystem — ten typed-state migrations (ADRs
0030–0041), then a series of per-symbol "globals retirement" sweeps —
until `variables.f90` and its companion `arrays.fi` were **deleted
outright** (2026-05-25). The `config_to_variables` adapter was then
dissolved into per-subsystem `state%X%init()` constructors, the
transitional `state%cfg` back-pointer was retired, and the last
transitional scaffolding (`legacy_state_t`) was removed. Throughout, the
verification oracle was **byte-identical regression**: the modern binary
had to reproduce a reference output bit-for-bit at every commit, first
against a self-baseline and later against an unmodified SWAP 4.2.0
recompiled with the modern toolchain (`swap420gf`). The headline result:
every runtime symbol now lives on a typed `state%X` (mutable) or
`config%X` (read-only) record; the model is materially closer to
in-process multi-instance execution and a clean Python binding.

This document is the durable record of that work; it is written to stand
alone after the 52 working specs/plans it draws from are deleted.

## Starting point (end of Phase 4)

At the close of Phase 4 (HEAD on `development`, ~2026-05-08), the
codebase had:

- A working **TOML-only runtime**. The legacy fixed-format readers were
  physically deleted (ADR 0019); production input flowed
  `TOML → swap_config_t → config_to_variables → variables.f90 globals`.
- **`variables.f90`: ~1,500 active declarations**, `save`d, read by every
  subsystem. `arrays.fi` provided the array-dimension include.
- **`config_to_variables.f90`: ~1,800 lines** of `<global> = config%X%Y`
  copy assignments — the adapter that bridged typed config into the bare
  globals.
- A pFUnit harness (~610 tests) and a `check-full` regression suite of
  **5 byte-identical cases** (hupselbrook, grassgrowth, oxygenstress,
  salinitystress, surfacewater; macroporeflow excluded per ADR 0011).
- An ADR backbone through 0028 (infrastructure, the TOML pipeline, the
  nutrient-subsystem reactivation work).

The aggregator type `swap_state_t` existed but was nearly empty — a
shell awaiting subsystem records. The post-Phase-4 work is the story of
filling it, then collapsing everything else into it.

---

## Arc-by-arc chronological narrative

The work falls into five overlapping phases:

1. **State-type migrations** (2026-05-08 → 05-12, ADRs 0030–0042) —
   carve each subsystem's runtime globals into a typed `state%X` record.
2. **Globals retirement series** (2026-05-12 → 05-25, GR-* arcs) —
   per-symbol migration of the remaining (mostly crop/config) globals,
   ending in the **deletion of `variables.f90`**.
3. **Adapter dissolution** (2026-05-25, ADR 0045 predecessor) —
   `config_to_variables` → per-subsystem `state%X%init()`.
4. **Clarity & modernization arcs** (2026-05-26 → 05-28) — task-dispatch
   retirement (ADR 0043), typed CSV records (ADR 0044), CSV-output
   modernization, narrative-helper decompositions, regression-coverage
   expansion.
5. **Architectural close-out** (2026-05-28 → 05-29, ADRs 0045–0047) —
   orchestrator dissolution, `state%cfg` retirement, post-cfg polish,
   `legacy_state_t` deletion.

### Phase 1 — State-type migrations (ADRs 0030–0042)

The unifying method was a **four-stage playbook**: *discovery* (catalog
owned globals, borrowed globals, entry points, cross-subsystem readers,
hazards) → *design* (define the state record, resolve each hazard, decide
phasing) → *plan* (mechanical task decomposition with `check-full` as the
integration gate at every commit) → *execute* (subagent-driven, with a
**dual-write transitional pattern**: compute writes both the new typed
state and the legacy global until every reader is cut over, then the
legacy write is dropped). The playbook accumulated 32 lessons in
`state-migration-playbook.md`.

Each migration added exactly one field to `swap_state_t`. The order was
chosen to climb from peripheral subsystems toward the coupled core.

#### Surface-water pilot — ADR 0030 (2026-05-08 → 05-09)

The pilot. Surface-water was chosen as peripheral (bottom-of-column,
outside the crop/water-balance feedback loop) yet substantial (~1,433
LoC home files). 29 owned globals migrated into `surfacewater_state_t`,
threaded from `swap_main` down through `ASSOCIATE` blocks in leaf
routines (a `sw_*` alias prefix avoided shadowing `use variables`
imports). Two phases: dual-write + output-read redirect (Phase 1), then
cross-subsystem reader migration across 13 files + global deletion (Phase
2). The `fldecdt` "decrease timestep" signal global was replaced by an
`intent(out) :: request_smaller_dt` parameter hosted in a tiny new
`timestep_control_mod`. One global (`qdra`) was retained as a
divdra-argument-binding scaffold pending drainage. **Outcome:** 28/29
globals retired; pFUnit 614→617; check-full 5/5 byte-identical at every
task. The pilot produced the discovery-template's critical Section 3.5
("external readers of owned globals") — the largest source of surprise
scope.

#### Drainage — ADR 0031 (2026-05-10)

6 owned globals into `drainage_state_t`. `qdra`/`qdrain` re-homed from
surfacewater (their architecturally-correct owner). `divdra` modernized
from explicit-shape `(Madr, macp)` arrays to assumed-shape, finally
freeing the legacy `qdra` global. Key lesson reinforced: the
cross-subsystem reader inventory must include **compute** readers, not
just output (a missed `qdrd` reader in `WLEVBAL` caused a ~0.4 cm GWL
drift caught only at the integration gate). A `geofac` adapter ordering
bug was fixed as a side effect. check-full 5/5 throughout.

#### Solute — ADR 0032 (2026-05-10)

25 fields into `solute_state_t`. Three distinctive findings: (a) a
**Phase-0 physics-config gap** — 14 fields (`cref`, `kf`, `frexp`,
`decpot`, …) consumed by `solute.f90` but absent from `solute_config_t`,
meaning every TOML `swsolu=1` run silently used zero defaults (a real
correctness bug, patched). (b) **AgeTracer is dead code** (`flAgeTracer`
never set true) — its 244-line body was extracted to a new
`agetracer.f90` module behind a stub-error guard, preserving it verbatim
for future reactivation while uncluttering the solute path. (c) the
two-stage `cml` seeding (config CSV → global → typed state via a new
`solute_init`) was preserved. 18 of 27 globals retired; 9 retained
(config-seed bridges + AgeTracer dead-code deps).

#### Heat — ADR 0034 (2026-05-10)

13 owned globals, **all instantaneous** — the first **flat** state record
(no cohorts). Closed the last main compute entry point without state
(`Temperature(task)`). The `outheapar` "writeback" hazard was diagnosed
as a parameter sweep, not state mutation (fixed with a local scratch
array). `tsoil` was reclassified as a config-staging buffer (not
retired). 871 lines of dead `csv_write`/`csv_write_tz` deleted. A
Phase-0 config gap (`swcalt=1` analytical method, 6 fields) was patched.

#### The soil-water four-arc decomposition (ADRs 0035–0038)

Soil-water touched ~100 globals — too large for one arc. It was
**decomposed by coupling surface**, not by phase: **boundary →
crop-uptake → atmosphere → core**. Each arc carved the fields belonging
to its surface; the final core arc handled the Richards interior. This
minimized per-commit blast radius and produced clean ADR narratives. A
deliberate design choice: `soilwater_init` was placed early (after
`CalcGrid()`, before `DoTillage(1)`) so later arcs could add per-node
arrays without an init-order gap.

- **Boundary — ADR 0035 (2026-05-10).** 12 top/bottom coupling scalars
  (`qtop`, `qbot`, `hbot`, `gwlinp`, …) — the first population of
  `soilwater_state_t`. Phase 0 promoted 8 config fields (swbotb=2 sine /
  4 q(h)-exp / 8 lysimeter). ~165 reads across 8 files; `pond`/`gwl`/
  `kmean` deferred to core.

- **Crop-uptake — ADR 0036 (2026-05-11).** 23 root-sink globals; the
  first **per-node allocatable arrays** in the record. `soilwater_init`
  signature stabilized at `(sw, numnod, nlay)`. `CropGrowth(1, state)`
  plumbed. The `mfluxtable` lookup build had to move from the three crop
  init paths up to the dispatcher (the per-layer hydraulic params it
  needs are populated later by `SoilHydraulics(1)`). JvL machinery
  (`swdrought=2`) is stub-errored at parse time — migration was
  byte-identical-free by construction.

- **Atmosphere — ADR 0037 (2026-05-11).** 40 globals — the largest arc to
  that point. First arc spanning **five reset cadences** (instantaneous,
  per-day, per-event, intermediate, cumulative), making the cohort
  pattern (ADR 0033) mandatory and debuting it **at type creation**. A
  new playbook lesson — the **pre-flight dual-write coverage check** —
  was surfaced after an `nraidt`/`aintcdt` regression: before reader
  cutover, verify each migrated field has a non-zero state write beyond
  zero-init. 13 hidden readers across 11 files surfaced during the
  compile-driven phase; 4 routines (`DoTillage`, `CNmethod`,
  `checkmassbal`, `Consolidate_Bdens`) gained state args.

- **Soil-water core — ADR 0038 (2026-05-11/12).** **74 globals — the
  largest single-arc retirement in the series** — the Richards interior
  (`h`, `theta`, `q`, `kmean`, `pond`, `gwl`, cohort accumulators). This
  arc resolved every prior soil-water deferral and validated the most
  important methodological advance of the whole effort: **Strategy B
  compile-driven retirement** (see Methodology). ~400+ read sites across
  23 files. The **module-level pointer pattern** was introduced for
  `soilhydraulicsutils.f90` and `WC_K_models_04_11.f90` (called from
  inside the Newton-Raphson loop): rather than thread state through 81+
  callers, an init-time `bind_state_targets` points module-level pointers
  at the state arrays. After all four arcs, `soilwater_state_t` carried
  ~149 fields.

#### Cumulative reset cohorts — ADR 0033 (2026-05-10) → flatten — ADR 0042 (2026-05-12)

ADR 0033 made reset cadence first-class: each subsystem record gained
nested `<subsystem>_intermediate_t` / `<subsystem>_cumulative_t`
sub-records, each with a single type-bound `reset()`. A key correction
during the surfacewater rollout established **cohort partitioning by
activity gate**: when cumulatives accumulate under different flags
(surfacewater's `flSurfaceWater` vs drainage's `fldrain`), split the
cohort by flag, not by name or type — the owner of each cohort is the
subsystem whose flag gates its accumulation. Applied across surfacewater,
solute, soilwater.

Two days later, **ADR 0042 walked the wrapper back.** In review the
cohort sub-records bought less than expected: gates were still enforced
at call sites, the `i*`/`c*` name prefixes already encoded role, the
cohort-vs-flat boundary on the parent was inconsistent, and hot loops
traversed three levels. The asymmetric semantics the cohorts captured
were re-expressed as **named reset procedures on the flat parent type**
(`reset_intermediate`, `reset_cumulative_drainage`,
`reset_cumulative_reservoir`, `reset_intermediate_per_day`), with
section-header comments restoring the self-documentation. ~270 call sites
migrated across ~22 files; 8 commits; byte-identical throughout. ADR 0042
supersedes only the cohort-pattern portion of ADR 0033 (0033 left intact
as historical record).

#### Tillage — ADR 0039 (2026-05-12)

The smallest arc: 13 runtime-state `till_*` globals, flat. Notable for
distinguishing **config-constants from runtime-state** — 18 of the 31
`till_*` globals were immutable-after-init and already in typed
`soil_tillage_t`, so they were *deferred* (only runtime state migrates).
Only **2 Strategy-B compile passes** (smallest ever); reader cutover was a
no-op (zero external readers — fully self-contained). An in-arc latent
bug fix (`set_iTill` tautological condition that never advanced past
event 1).

#### Macropore retirement — ADR 0040 (2026-05-12)

Not a migration but a **deletion**. Macropore flow had been deferred
since ADR 0010/0011; every state arc lugged ~80 macropore globals through
just to keep them compiling. The cost of carrying ~4,500 LoC of unused
physics now exceeded the cost of deleting and re-implementing later.
Deleted: `macropore.f90` (2,245 LoC), `macrorate.f90` (911),
`macroporeoutput.f90` (309), the orphan typed config, the pFUnit suite,
~115 lines of declarations + ~100 of inits — **~3,870 LoC total**.
`flMacroPore` forced permanently `.false.`; remaining cross-subsystem
`if(flMacroPore)` branches dead-code naturally; `FrArMtrx` defaults to
`1.0` (whole-matrix). ~25 placeholder "retired-zero" globals kept where
dead branches still reference them. The legacy implementation is
preserved on `legacy/swap-4.2.0` as a future TDD oracle; ADR 0040
documents the re-implementation path.

#### TimeControl — ADR 0041 (2026-05-12)

The 10th and **most cross-cutting** migration. 62 runtime-state globals
(16 clock/calendar + 19 timestep/schedule + 27 runtime booleans) into a
flat `timecontrol_state_t`. **~415 read sites across 22 files** — the
largest reader inventory of any arc (already reduced from ~500 by the
ADR 0040 macropore retirement and the ADR 0009 Phase-5+ output purge).
Reads were syntactically explicit (`use variables, only:`), so Strategy B
surfaced everything with no NR-coupled hang risk (8 compile passes). The
`flZeroIntr`/`flZeroCumu` reset gates and 18 config-constants were
deferred. A reusable lesson: the **ASSOCIATE single-entry/single-exit
constraint** — routines with early `return`s inside an associate scope
(`PONDRUNOFF`, select-case bodies in `management_soil`) cause a gfortran
ICE; the fix is a local copy seeded from state at routine top.

**State of the aggregator at the close of Phase 1 (2026-05-12):**
`swap_state_t` carried 8 typed records — `surfacewater`, `drainage`,
`solute`, `heat`, `soilwater` (149 fields across 4 arcs), `atmosphere`
(40), `tillage` (13), `timecontrol` (62). `variables.f90` was down from
~1,500 to ~700 active declarations. pFUnit ~733; check-full 5/5
byte-identical; hupselbrook canary stable at ~1.1 s.

### Phase 2 — Globals retirement series (GR-* arcs, 2026-05-12 → 05-25)

The state migrations cleared the subsystem *runtime* globals. The
remaining ~700 declarations were dominated by **crop biomass / WOFOST**
fields and **config-side switches** — harder, because crop globals are
read across many files (the `cropgrowth.f90` dispatcher + its split
runtime files + irrigation + oxygenstress) and config switches need
consumers rerouted to read `config%X` directly rather than a global copy.

This phase was a long series of subagent-driven sweeps. ~8 arcs shipped
2026-05-12 → 05-16 (GR-UTILS, GR-BH, GR-ATM, GR-CROP, GR-FINAL,
GR-CROPRT, GR-CROPWS), then per-symbol arcs through 05-25. Headline
structural win: **`cropgrowth.f90` split from 5,431 → 918 lines** (the
dispatcher), with four focused runtime files extracted
(`cropwofost_runtime` 2,044, `cropgrass_runtime` 1,444, `cropgrowth_helpers`
768, `cropfixed_runtime` 297). New state sub-records: `state%mesh`
(vertical discretization), `state%crop` (`%common` ~95 fields, `%fixed`,
`%wofost`, `%grass`), `state%nutrients` (~100 fields from
`Wofost_Soil_Declarations`).

Three methodological refinements made the long tail tractable:

- **Per-symbol-across-all-files migration** beat per-file. Subagents
  defaulted to retaining dual-writes and deferring write-site migration;
  per-file conversion breaks regression when a *different* file still
  reads the bare global (it sees a stale value). Migrating a symbol in
  *every* file at once (reads first, then the write, then the
  declaration) is the safe grain. `scripts/retire_global.py` and
  `scripts/show_consumers.py` (a locate-only categorizer) automated the
  playbook. The **GR-DVS arc** (2026-05-19, 29 commits, ~59 globals)
  proved the inline pattern on crop biomass symbols (dvs, lai, wlv, …).
  **GR-CROPWS Sweep 2** (2026-05-20 → 05-22, 31 commits) retired **137
  per-rotation crop fields**, bringing `variables.f90` to 477 active
  declarations.

- **The `state%cfg` back-pointer** (introduced 2026-05-19):
  `swap_state_t` gained `type(swap_config_t), pointer :: cfg`, set once at
  init. Top-level compute routines could read config switches via
  `state%cfg%X%Y` with no signature churn — a fast strangler accelerant
  for the hundreds of config-side switches. (It was itself retired later;
  see ADR 0046.)

- **Sweep 3 transitional bag** (2026-05-22): a `state%legacy` record (437
  allocatable-array fields) was staged as a holding pen for the
  still-active globals, to allow `variables.f90` to be deleted while
  consumers migrated incrementally. (A linker GOTPCREL overflow forced
  the fields to be `allocatable`, not fixed-size — a recurring storage
  lesson.)

**The atmosphere cleanup arc (2026-05-22 → 05-23)** brought `src/atmosphere`
to **zero `use variables` imports**, retiring ~40 bare globals and
establishing the template patterns for subsequent subsystem cleanups
(sub-record `associate` aliasing, config-out-of-compute, dual-write
retirement). It surfaced two latent bugs: `swap_log` file output was
unreachable (a `log_unit > 0` guard rejected the valid `-1`-vs-negative
unit convention, so debug logs were always 0 bytes), and `Runoff_CN`
silently accumulated zero under `swuseCN=1` (no writer after a migration;
not caught by check-fast because no case exercises CN runoff).

**The crop `use variables` sweep (2026-05-25, 52 commits)** rendered
`src/crop` locally use-variables-free, adding `state%crop%oxygen` (23
Bartholomeus fields) and `state%crop%irrigation` (14 SSDI fields), and
extracting two dormant modules (`jongvanlier.f90`, `oxygenrepro.f90`) to
release JvL-only globals. Parallel `src/io` (22 globals), drainage (13),
and timecontrol (10) sweeps cleared their residuals.

#### Terminal state — `variables.f90` + `arrays.fi` DELETED (2026-05-25)

In a final 5-commit "Phase 6," the last consumers of the bare-globals
module were eliminated, `config_to_variables` was stripped of its ~35
remaining mirror writes, `initialize.f90` was reduced to no-op stubs
(595 → 30 lines), and **`src/core/variables.f90` (~1,330 lines) and
`src/core/arrays.fi` were deleted** (commit `1938653`). Every runtime
symbol now lives on `state%X` (typed runtime state), `config%X` (typed
config), or the `swap_array_dimensions` module (size parameters as integer
`parameter`s). **The strangler pattern's bare-globals layer was gone.**
Per-step `rm -rf builddir && check-fast` was byte-identical at all 5
steps (748 pFUnit + 4/4 regression).

(The `legacy_state_t` holding bag from Sweep 3 still existed at this point;
it was deleted later, in ADR 0047.)

### Phase 3 — Adapter dissolution / seed-state arc (2026-05-25)

With the globals gone, the `config_to_variables` adapter (now 1,631 lines
of typed-config-to-typed-state copies) was decomposed in 12 tasks / 14
commits into **per-subsystem type-bound `state%X%init(config%Y, …)`
constructors**. Five subsystems' free `*_init` subs were promoted to
type-bound; five existing inits were extended. The adapter shrank to a
**302-line orchestrator** (renamed `seed_state_from_config.f90`) holding
only ~5 residual cross-subsystem writes. `swap_init` flattened to a
2-deep chain (modulo a thin CAPI shim that loads config from a string
buffer rather than a file). A runtime SIGSEGV during this arc taught that
**init order matters and is not statically checkable**: `state%crop%init`
must run before `timecontrol_init` (which reads
`state%crop%common%croptype`), discovered only on a clean rebuild.

### Phase 4 — Clarity & modernization arcs (2026-05-26 → 05-28)

#### CSV-output modernization / IO-OUT arc — (2026-05-26)

CSV output was modernized into `csv_writer` (a generic primitive
symmetric to `csv_reader`), `output_registry` (a 112-entry
single-source-of-truth catalogue), and `csv_output` (public
init/step/finalize). **`swapoutput.f90` was deleted in full**, and the
**entire per-subsystem `output_row` C-API mechanism** (`swap_get_output_row`
+ 9 streams + every `*_output_row`/`*_columns` field across 8 state
files) was retired as unused — **net −1,814 lines**. +9 crop/snow parity
columns were added before deletion to keep `result_output.csv`
byte-identical. The design rationale ("unify-strongly"): Python use is
batch (seed → run → read), with per-step BMI access deferred, so the
per-step stream zoo was retired rather than re-homed.

#### Task-dispatch retirement — ADR 0043 (2026-05-27)

The legacy `subroutine X(task)` + `select case(task)` lifecycle-dispatch
idiom (case 1 = init, 2 = step, 3 = finalize) was retired from **10
Tier-1/Tier-2 dispatchers** (irrigation, SSDI, csv_out, temperature,
solute, tillage, soilwater, surfacewater, …). The key insight: "init" is
genuinely **two phases** that cannot merge — Phase 1 *construct*
(type-bound `state%X%init(cfg)`: allocate/zero/snapshot, order-free) and
Phase 2 *seed* (a free `x_seed(state)` computing derived initial state
from *sibling* subsystems, order-dependent: soilwater → heat → solute).
Folding the seed into the type-bound init would either drop the derived
computation or break ordering. So `x_seed` stays at its original call
site (order preserved byte-for-byte), and step/sub-step cases became named
free procedures (`x_step`, `x_update`, `x_finalize`). ADR 0043
established a binding **multi-instance invariant**: every extracted
procedure takes `state` explicitly and touches nothing global (no module
vars, no `SAVE`). The arc also identified **the real parallelism blocker**
— not the dispatch idiom but module-global state in the deferred crop
cluster (`crop_config_global` pointer + `SAVE` locals in cropgrowth/
cropgrass). A **Tier-3 crop-cluster arc** (on an isolated `tier3-crop`
worktree) subsequently retired SP1–SP2d (crop_config_global,
ArableLandGerm+CropFixed, Grass, Wofost, CropGrowth parent) SAVE-free; SP3
(SoilManagement + WSN nutrient module-SAVE) and SP4 remained.

#### Typed CSV record tables — ADR 0044 (2026-05-28)

Post-strangler, CSV companion files were cached as untyped `(:,:)` arrays
on state, with the schema split across three places (loader header,
consumer column indices, unit code). ADR 0044 gave each CSV family a typed
module under `src/io/csv/<family>.f90`: a `<family>_row_t` record (one
named field per column) and a `<family>_table_t` (`rows(:)`, `is_loaded`,
`load(path, errors)` / `year_window(...)`). Rules: `load` validates inline
(no separable `validate`); unit conversions stay at the consumer; date
columns stay `real64` days-since-1900. Piloted on meteo (daily/detail/rain),
then applied to **7 families** (drainage owltab, irrigation fixed+SSDI,
nutrients amendments, bottom-boundary 5 sub-tables, soil initial profiles)
— 14 typed table types across 6 modules. A `days_since_1900` helper was
de-duplicated (4 copies). This arc closed **W3** (config mutation in
`soilwater_state_init` — `config_soil` demoted to `intent(in)`) and **W4**
(duplicate h_file CSV read). pFUnit 769 → 812.

#### Narrative-helper decompositions; regression-coverage expansion (2026-05-27 → 05-28)

A **solute pure-kernels pilot** (2026-05-27) tested extracting physics
into pure kernels, then — after feedback — was reworked into the
**narrative-helpers** pattern: decompose a big subroutine into a thin
orchestrator + named phase-helpers operating on `state`, keep genuine
multi-line algorithms as standalone `pure` functions (the Freundlich
fixed-point solver), keep control flags in the orchestrator. It surfaced
two latent bugs: a masked uninitialized-coefficient read (the seed/step
split discarded 5 derived arrays), and a div-by-zero NaN on `swbr=1` that
**exists in SWAP 4.2.0 itself** (no valid parity reference) — quarantined
with a `fatalerr_collected` guard. The pattern was then applied to the
**soil** and **drainage** cleanup arcs (2026-05-28): `headcalc` deduped
(one `headcalc_residual`), `WC_K_models` made **reentrant** (removing a
soil-side parallelization blocker), `WLEVBAL` (~335 L) split into four
phase helpers, cross-file `qdra`-redistribution dedup, and legacy
`goto`/numbered-`do` modernization across drainage.

A **regression-coverage-expansion arc** (2026-05-27) replaced the
self-baseline with **`swap420gf`** — an unmodified SWAP 4.2.0 recompiled
with the modern build's gfortran flags (zero source edits) — making the
suite a *physics* fidelity check with the compiler held constant
(verified: 4.2.0-ifx ≡ 4.2.0-gfortran, 0.000 drift). It added two cases
(`soilhysteresis`, `winter`) and, via systematic debugging, traced their
~1.9 cm / ~0.04 cm divergence not to physics or compiler but to the
**adaptive-`dt` controller desyncing on the first timestep** (modern's
first Richards solve converges in `numbit=1` vs 4.2.0's `numbit=4` from a
bit-identical init; the faithful `dt`-doubling rule then walks different
timestep sequences). With fixed `dt` the builds are bit-identical.
Decision: accept and document as `xfail` (`known_divergence`), since
fixing would perturb the core solver to replicate a sub-1e-6 quirk.

### Phase 5 — Architectural close-out (ADRs 0045–0047)

#### Orchestrator dissolution — ADR 0045 (2026-05-28)

The strangler adapter `seed_state_from_config.f90` (8 live writes amid
~290 lines of comment-graveyard) was dissolved in 9 byte-identical steps,
each folding a cross-subsystem write into the appropriate
`state%X%init()` and then **deleting the file entirely**. `swap_init_body`
became the orchestrator — init order explicit in one place, top-to-bottom.
Every state field now has a single owner; init signatures encode their
config dependencies explicitly (e.g. `state%soilwater%init(config%soil,
config%drain, config%heat, config%bottom_boundary, …)`). Closed W1+W2.
812 → 806 pFUnit (the deleted module's 6 tests).

#### `state%cfg` pointer retirement — ADR 0046 (2026-05-28)

The `state%cfg => config` back-pointer — the equivalent of bare-global
reads at the config layer — was retired. **159 live `state%cfg` reads
across 25 files → 0**, migrated via three patterns: direct config-parameter
substitution (consumer already takes `config`); snapshot-to-state-at-init
(consumer takes only `state`); or signature expansion (init-only
consumer). ~30 new snapshot fields were added across the state hierarchy.
Shipped as 7 clusters + a final delete commit (9 commits). A gfortran
gotcha: the `timecontrol_init` config parameter needed the `target`
attribute to prevent heap corruption. This brought the model materially
closer to multi-instance: with no global pointer, N concurrent
`swap_state_t` instances no longer share a config reach-through.

#### Post-cfg polish arc — ADR 0047 (2026-05-29)

A bundle of end-of-strangler housekeeping (none individually worth an arc):

- **TOML reader dedup (W9):** `toml_array_helpers.f90` extracted; 5 copies
  of `read_table_2d` + 3 of `read_array_1d` removed (**−335 lines**).
- **pathwork unification (W8):** `config%general%pathwork` became the
  canonical relative-path anchor (default = directory of the main TOML,
  was hardcoded `'./'`); two competing path conventions collapsed to one.
- **Crop subfile-loader extraction (W5):** `.crp.toml` loading/dispatch
  moved to `load_crop_rotation_files.f90`; `read_crop_toml.f90` became
  pure parsing.
- **Nutrients WSN dual-write retirement (W10):** the 8 dual-written
  `Wofost_Soil_Declarations` globals retired to `state%nutrients` — **the
  last remaining "state init writes to globals" pattern in the codebase.**
- **`legacy_state` deletion (W12):** a recount found every
  `legacy_state_t` field had zero live readers; the entire
  `src/state/legacy_state.f90` (379 lines), the `state%legacy` composition
  field, and all build registrations were **deleted**. The Sweep 3
  transitional bag is closed.
- **Obsolete-tooling cleanup:** 4 now-defunct Python scripts
  (`audit_variables.py`, `output_parity_check.py`, `retire_global.py`,
  `sweep1_config_side.py`) and the dormant `jongvanlier.f90` museum file
  deleted.

**Outcome:** the strangler pattern is fully retired end-to-end —
`variables.f90` → `legacy_state_t` → both gone; `swap_state_t` contains
only typed subsystem records, no transitional bag. Net ~−700 lines.

### Cross-cutting — Docs restructure (2026-05-25)

In parallel with the seed-state arc, the documentation site was
reorganized. Public docs under `docs/` were split into
**Tutorial / Reference / Developer / API** sections with a custom
Atlas-palette navbar and a lunr.js typeahead search. All developer
material — **42 ADRs, the superpowers specs/plans/notes, and the frozen
Phase-4 archive** — was moved to `dev-docs/` (out of FORD's scan path).
This is why the present document and its sources live under `dev-docs/`.

---

## Aggregate metrics

| Metric | Value |
|---|---|
| Window | 2026-05-12 → 2026-05-31 (~20 days) |
| Commits on `development` | ~851 |
| Net `src/` line delta | ≈ −5,500 (≈61.5k added / ≈67k deleted) |
| ADRs produced | **0030–0047** (18 ADRs) |
| `variables.f90` trajectory | ~1,500 active decls → **0 (file deleted)** |
| State-type migrations | 10 (surfacewater → timecontrol) |
| `swap_state_t` subsystem records | surfacewater, drainage, solute, heat, soilwater, atmosphere, tillage, timecontrol, mesh, crop (+oxygen/irrigation/common/fixed/wofost/grass), nutrients |
| `config_to_variables` adapter | ~1,800 L → 302 L → **deleted** |
| `state%cfg` reads | 159 → 0 (pointer deleted) |
| `legacy_state.f90` | 437-field bag → **deleted** |
| Largest deletions | `swapoutput.f90` (−1,814), macropore subsystem (−3,870), cropgrowth split (5,431 → 918) |
| pFUnit tests | ~610 → ~806–812 (suite grew with each arc) |
| Regression posture | `check-full` **5/5 byte-identical** non-xfail at every commit; 2 documented adaptive-`dt` xfails (`soilhysteresis`, `winter`) |
| Regression oracle | self-baseline → **`swap420gf`** (unmodified SWAP 4.2.0, modern toolchain) |

The −5,500 net `src/` figure *understates* the structural change: tens of
thousands of lines were rewritten (globals → typed state), and large
deletions (macropore, swapoutput, the output C-API) were partly offset by
new typed state records, init methods, typed CSV modules, and tests.

---

## Methodology & reproducible lessons

This is the most transferable part of the work — a playbook for
modernizing a legacy scientific Fortran codebase without altering its
numerics.

### The strangler-fig pattern, completed

The end-to-end transformation was **TOML → typed config → typed state**,
with the bare-globals layer (`variables.f90`) strangled out from under a
running, regression-locked system. The sequence per subsystem:
typed config schema (Phase 4) → typed state record threaded through
signatures → dual-write transitional period → reader cutover → drop the
dual-write → delete the global. The adapter (`config_to_variables`) and
the back-pointer (`state%cfg`) were themselves strangler scaffolding, and
both were eventually dissolved (ADRs 0045/0046). The final state has no
intermediate layer: every runtime symbol is on a typed `state%X` or
`config%X` record.

### Strategy B — compile-driven retirement (ADR 0038, matured by 0041)

The pivotal technique. Manual file-by-file dropping of legacy writes
caused a Newton-Raphson divergence in the soil-water core arc (utility
modules read stale zero-init arrays once their half-writes stopped). The
fix was to **invert the process**: comment out *all* the arc's globals in
`variables.f90` first, then let the compiler enumerate every remaining
reference as a compile error. Each error is a precise instruction (drop a
legacy write, migrate a missed read, or drop a use-clause import); iterate
fix-recompile. The compiler becomes the exhaustive inventory tool — there
is no manual gap and no "which-file-first" ordering problem. Compile-pass
counts measure subsystem coupling: tillage 2, boundary 8, TimeControl 8,
atmosphere 13, soil-water core 28. By ADR 0041 the pattern was at full
maturity (4th application).

A corollary technique for **NR-coupled hidden readers** (utility functions
called inside the Richards loop with many callers): the **module-level
pointer pattern** — bind module pointers at the state arrays once at init,
rather than thread state through 81+ call sites.

### The cohort pattern and its walk-back

A cautionary tale worth recording in full: a typed abstraction (nested
reset cohorts, ADR 0033) was introduced to make reset cadence first-class,
adopted across three subsystems, then **deliberately reverted** (ADR 0042)
when review showed it bought less than its cost (gates still enforced at
call sites; name prefixes already encoded role; hot-loop access paths
deepened). The asymmetric semantics it captured were re-expressed as named
procedures on the flat type. Lesson: prefer the simplest representation
that expresses the contract; a type-shaped wrapper is not free.

### Two-phase init (ADR 0043)

Initialization is genuinely two phases: **construct** (order-free,
idempotent, per-subsystem, type-bound `state%X%init(cfg)`) and **seed**
(order-dependent, reads sibling subsystems, a free `x_seed(state)` left at
its original call site to preserve ordering byte-for-byte). Conflating
them into a single `case(1)` body — or trying to fold the seed into the
type-bound init — drops the derived computation or breaks ordering. This
distinction is the backbone of the multi-instance roadmap.

### Byte-identical regression as the verification oracle

The non-negotiable gate. Every commit had to reproduce a reference output
bit-for-bit (`check-full`, later expanded against `swap420gf`). This made
the entire refactor a sequence of provably behavior-preserving steps, and
it repeatedly caught regressions that pFUnit alone missed (global-default
drift). Key disciplines:

- **Per-task regression gating, not end-of-arc.** Every implementer
  subagent ends with `check-fast`; deferring to an end-of-arc gate makes
  regressions hard to bisect. Several silent regressions were caught only
  because each commit was individually byte-identical (GR-ATM, GR-CROPWS).
- **`check-full` over pFUnit alone.** pFUnit misses global-default
  regressions; the byte-identical suite is the dominant correctness signal.
- **The hupselbrook canary.** An NR-coupled hidden-reader hang manifests
  as a 99.9%-CPU hot-spin with no output; kill at ~5 s and investigate the
  Strategy-B compile residuals.
- **`swap420gf` as physics oracle.** Recompiling unmodified 4.2.0 with the
  modern toolchain (zero source edits) holds the compiler constant, so any
  divergence is attributable to the modern *source*, not the compiler —
  and proved the `soilhysteresis`/`winter` divergence was an adaptive-`dt`
  scheduling artifact, not a physics error.

### Narrative-helpers over micro-functions

For decomposing large subroutines: a thin orchestrator + named
phase-helpers operating on `state`, with genuine multi-line algorithms as
standalone `pure` functions and control flags kept in the orchestrator —
*not* a zoo of one-line elemental kernels (which a pilot tried and was
rejected for poor readability). Applied to solute, soil (`headcalc`,
`soilwater_seed`), drainage (`WLEVBAL`, `bocodre`). A side benefit: the
exercise surfaces latent bugs (masked uninitialized reads, NaN paths).

### What blocked progress, and how it was unblocked

- **Init-order regressions.** Legacy globals doubled as init-order sync
  points; removing them exposed state reads happening before state writes
  (GR-DVS flCropX batch; the crop/timecontrol SIGSEGV). Not statically
  checkable — surfaced only on clean rebuilds. Unblocked by sequencing
  inits explicitly in `swap_init_body` and deferring fragile batches.
- **NR-coupled hidden readers.** (See Strategy B / module-pointer pattern.)
- **FP summation-order differences.** Algebraically-equal expressions
  summed in different term orders round differently under IEEE-754,
  blocking naive dedup (the `headcalc` residual). Safe to unify *only* when
  the divergent branch is dark in the regression suite — verified per case
  (e.g. SWBOTB=1 unused).
- **SAVE / module-global multi-instance blockers.** `crop_config_global`
  (a reassigned module pointer) and `SAVE` locals in the crop cluster are
  the literal blocker to in-process parallelism. Addressed incrementally
  by the Tier-3 crop arc; `WC_K_models` was made reentrant in the soil
  cleanup arc.
- **Subagent conservatism.** Subagents default to retaining dual-writes
  and deferring write-site migration; per-file conversion breaks regression
  across files. Unblocked by the per-symbol-across-all-files grain plus
  locate-only tooling (`show_consumers.py`) feeding a manual checklist.
- **Build-system gotchas.** State-schema changes don't propagate `.mod`
  deps across the swap_modern↔swap_legacy library boundary —
  `rm -rf builddir` is mandatory after any `src/state/*_state.f90` edit.
  Embedding many fixed-size arrays in a derived type triggered a linker
  GOTPCREL overflow — array fields default to `allocatable` (also the
  W10/iamend SIGSEGV lesson: fixed-size state arrays land on the stack via
  `swap_state_t` and blow smaller stack limits).
- **Fortran-specific traps.** ASSOCIATE single-entry/single-exit (early
  `return`s cause an ICE — use a local copy); CamelCase use-clause
  mismatches (case-insensitive grep misses them — surface at compile time);
  dangling `&` continuations when retiring the last symbol in a use-list;
  pFUnit's dark-test trap (a new `.pf` runs silently unless registered in
  `testSuites.inc` — verify the `OK (N tests)` count rises).

### Process conventions

- **Development-only commits**, never pushed without explicit ask;
  umbrella specs decomposed into per-task plans; CSV over inline for big
  tables; verify *state*, not subagent reports.
- **Tombstone provenance markers** (`! [SS-X] retired <date> — moved to
  state%Y`) preserve migration history at the deletion site.
- **Subagent-driven development** with a controller verifying each
  implementer's report against byte-identical `check-full` + pFUnit; the
  two-stage spec/code review used selectively (architectural tasks) and
  skipped for highly-repetitive mechanical cutovers once the pattern was
  established.

---

## Open threads / deferred work

- **Macropore re-implementation** as a fresh feature against
  `legacy/swap-4.2.0` as TDD oracle: a `macropore_state_t`, TOML-native
  config, no globals; re-add `3.macroporeflow` to the regression matrix;
  retire the ~25 "retired-zero" placeholder globals (ADR 0040).
- **AgeTracer reactivation:** un-stub `agetracer.f90`, define
  `agetracer_state_t`, narrow its wildcard `use Variables`, retire the 7–12
  dead-code-dependency globals (ADR 0032).
- **`swbr=1` solute physics:** correct the mixed-reservoir aquifer
  breakthrough NaN (the bug exists in 4.2.0 itself, so no parity reference
  — needs validation against the SWAP theory manual; possibly an ADR).
- **Config-constant consolidation:** the deferred config-constants (18
  TimeControl, 18 tillage Groups A+B, atmosphere/boundary Phase-0 fields)
  already have typed homes; an arc could remove the residual legacy
  representations. (Largely subsumed by the `state%cfg` retirement, but
  worth a closing audit.)
- **Reset-orchestration arc:** convert `flZeroIntr`/`flZeroCumu` (owned by
  TimeControl, consumed by ~25 files) into return values from a
  `state%timecontrol%advance()` (ADR 0041 D10).
- **Tier-3 crop dispatch — remaining:** SP3 (`SoilManagement` dispatch +
  the `Wofost_Soil_Declarations` blanket-`save` ~100-var module) and SP4
  (output file-handle SAVEs). These are the **last multi-instance
  blockers** in the crop cluster (ADR 0043 forward note).
- **Multi-instance / parallelization:** with globals, the adapter, and
  `state%cfg` gone, and `WC_K_models` reentrant, the path is mostly the
  remaining crop SAVE state. A possible end state is a single top-level
  `state%init(config)` wrapping construct → seed → first-day compute so N
  instances cannot mis-sequence.
- **Grid-dimensions arc:** `numnod`, `numlay`, `dz`, `z`, `disnod`,
  `layer`, … remain as intentional legacy globals pending a `grid_t`
  (partially started via `state%mesh`).
- **Regression coverage:** many `SW*` option clusters are stub-guarded on
  the TOML path (swdrought=2, swco2=1, swsalinity=2, grazing, …) — dead
  until the dropped physics is re-implemented; buildable non-crop clusters
  (swqhbot=1, swbotbhea=2, swnrsrf, swdislay, swsec=1, swqhr=2, swrain=1,
  swinco=1, swetsine) are candidates for new fixtures.
- **CSV-output follow-ups:** wire `resolve_inlist` (retire legacy
  `make_userlist`/`det_which_vars`); a unified in-memory results accessor
  for batch-Python (replacing the deleted `output_row` C-API) extensible
  to per-step BMI.
- **CAPI path resolution (W7):** hardcoded `base_dir = './'` — deferred as
  part of a BMI/CAPI arc.
- **`development` not pushed:** the entire post-Phase-4 arc (~851 commits)
  remains local on `development`; a review/PR consolidation is a separate
  decision.

---

## ADR map (post-Phase-4)

| ADR | Subject | Date |
|---|---|---|
| 0030 | State migration — surfacewater (pilot) | 2026-05-09 |
| 0031 | State migration — drainage | 2026-05-10 |
| 0032 | State migration — solute (+ AgeTracer extraction) | 2026-05-10 |
| 0033 | Cumulative reset cohorts (later superseded in part) | 2026-05-10 |
| 0034 | State migration — heat (first flat record) | 2026-05-10 |
| 0035 | State migration — boundary (soil-water arc #1) | 2026-05-10 |
| 0036 | State migration — crop-uptake (arc #2; per-node arrays) | 2026-05-11 |
| 0037 | State migration — atmosphere (arc #3; cohorts at creation) | 2026-05-11 |
| 0038 | State migration — soil-water core (arc #4; Strategy B origin) | 2026-05-11 |
| 0039 | State migration — tillage (smallest arc) | 2026-05-12 |
| 0040 | Macropore subsystem retirement (−3,870 LoC) | 2026-05-12 |
| 0041 | State migration — TimeControl (largest reader inventory) | 2026-05-12 |
| 0042 | Flatten reset cohorts (supersedes 0033 cohort portion) | 2026-05-12 |
| 0043 | Retire integer task-dispatch (two-phase init) | 2026-05-27 |
| 0044 | Typed CSV record tables | 2026-05-28 |
| 0045 | Orchestrator dissolution (`seed_state_from_config` deleted) | 2026-05-28 |
| 0046 | `state%cfg` pointer retirement (159 → 0 reads) | 2026-05-28 |
| 0047 | Post-cfg polish arc (`legacy_state` deleted; strangler fully retired) | 2026-05-29 |
