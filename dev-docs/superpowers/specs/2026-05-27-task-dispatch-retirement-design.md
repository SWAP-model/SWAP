---
title: "Task-Dispatch Retirement — Tier 1 + Tier 2 (init/step/finalize de-multiplexing)"
date: 2026-05-27
status: approved
branch: development
scope: Tier 1 + Tier 2 dispatchers (crop cluster / Tier 3 deferred)
---

# Task-Dispatch Retirement — Tier 1 + Tier 2

## 1. Problem

Several compute subroutines still carry the legacy SWAP idiom of an integer
lifecycle selector — `subroutine X(task, …)` with an internal
`select case (task)` where the cases are *lifecycle phases*, not algorithm
options:

- `case (1)` = initialization
- `case (2)` = per-step / rate computation
- `case (3)` = terminal / output

The smell is loudest at the call sites in `src/core/swap_mod.f90`, where the
orchestrator reads as a list of magic-int calls (`call SoilWater(1, state)` …
`call SoilWater(2, state)` … `call SoilWater(3, state)`) that the reader must
decode against the callee's case labels.

This must be distinguished from **legitimate** `select case` usage where the
integer encodes a genuine algorithm/option choice (e.g.
`select case (model)` in `soil/WC_K_models_04_11.f90`,
`reduceva_apply`'s daily/sub-daily mode). Those are **not** targets.

The canonical orchestration pattern already in place is the proof of the
target shape: `swap_init` / `swap_run_step` / `swap_close` (free module
procedures, time loop in `swap_main.f90`), with each subsystem exposing a
type-bound `state%X%init(...)`. The orchestrator already calls free per-step
procedures with no integer (`call Drainage(state)`, `call RootExtraction(state)`).

## 2. Scope

This spec covers the **9 independent dispatchers** (Tier 1 + Tier 2). The
crop cluster (Tier 3 — `CropGrowth`, `CropFixed`/`Wofost`/`Grass`,
`ArableLandGerm`, `SoilManagement`, `sumttd`, `Wofost_SoilRateConstants`) is
deferred to a follow-up spec because those dispatchers call **each other**
with magic ints and carry internal phase-sequence ordering that needs its
own design.

**Out of scope (confirmed not anti-patterns):**
`functionvalue_04_11(model)` (algorithm choice), `reduceva_apply(task_flag)`
(daily/sub-daily mode; mild redundancy with its `daily` arg, low priority),
`SharedSimulation(task)` (file I/O; dormant).

## 3. Target convention (the recipe)

For each dispatcher `X(task, …)`:

- **init case → free `x_init(state)` facade.** The facade internally
  `call state%X%init(state%cfg%Y, …)` (the typed, testable core) **and** runs
  the old `case (1)` runtime-seed logic. The orchestrator calls
  `call x_init(state)`. `state%cfg` is a `type(swap_config_t), pointer`
  associated at `swap_mod.f90:102` (before all inits), so the facade reaches
  every config subtree via `state%cfg%…` while taking only `state`.

  Rationale: a uniform free-proc API per subsystem
  (`x_init` / `x_step` / `x_finalize`) makes the orchestrator self-documenting
  and hides config-argument plumbing, while the type-bound `init` remains the
  typed core. The facade is a deliberate strangler layer — if it ever becomes
  a pure pass-through, deleting it is a trivial later step.

- **step / sub-step cases → free `x_step(state)`** (plus `x_update`, etc. for
  genuinely distinct phases), integer dropped.

- **finalize case → free `x_finalize(state)`.**

- **save/restore & lazy-table (Tier 2) → named pairs**, not lifecycle verbs.

The `x_init`-wraps-type-bound facade is the standard **only** for the
subsystems that have both a real init case and a partner type-bound init:
SoilWater, Temperature(heat), Solute, DoTillage(tillage). Stub-only inits
(irrigation, SSDI) get **no** facade — their init was already migrated to
config-load time, so only the step procedure remains.

## 4. Per-dispatcher migration map

Procedure names below are proposals; final names may be adjusted during
implementation but must drop the integer and read as named phases.

### Tier 1 — pure lifecycle dispatch

| # | Dispatcher | File | New free procedures | Init fold |
|---|---|---|---|---|
| 1 | `irrigation(task)` | `crop/irrigation.f90:31` | delete dead `case(1)` (no-op `return`); `case(2)` → `irrigation_step` | none (config-load init) |
| 2 | `SSDI_irrigation(iTask)` | `crop/irrigation.f90:~340` | delete dead `case(1)`; `case(2)` → `ssdi_irrigation_step` | none |
| 3 | `csv_out` / `csv_out_tz(iTask)` | `io/csv_output.f90:644 / :24` | inline the three cases into the existing `csv_output_init` / `csv_output_step` / `csv_output_finalize` wrappers; delete both dispatchers | n/a (output) |
| 4 | `Temperature(task)` | `heat/temperature.f90:82` | `temperature_init` (wraps `state%heat%init` + seed), `temperature_step` | yes |
| 5 | `Solute(task)` | `solute/solute.f90:10` | `solute_init` (wraps `state%solute%init` + seed), `solute_step` | yes |
| 6 | `DoTillage(iTask)` | `crop/tillage.f90:53` | `tillage_init` (wraps `state%tillage%init`), `tillage_step`, `tillage_output` (`case(3)`); **resolve `case(4)` liveness** (only 1/2/3 are called from the orchestrator) | yes |
| 7 | `SoilWater(task)` | `soil/soilhydraulics.f90:820` | `soilwater_init` (wraps `state%soilwater%init` + the large `case(1)` seed), `soilwater_step` (Richards/headcalc, `case(2)`), `soilwater_update` (rate+state, `case(3)`) | yes |

### Tier 2 — two-distinct-operations dispatch (looks like lifecycle, isn't)

| # | Dispatcher | File | New free procedures | Notes |
|---|---|---|---|---|
| 8 | `SoilWaterStateVar(task)` | `soil/soilhydraulics.f90:1203` | `soilwater_save_state` (`case(1)`, snapshot @ t), `soilwater_restore_state` (`case(2)`, restore on dt-reduce) | called internally `(1)` and from `swap_run_step` `(2)` |
| 9 | `MatricFlux(task)` | `crop/rootextraction.f90:329` | `matricflux_build_table` (`case(1)`; fold into soilwater init), `matric_flux(phead, node, outcome, state)` (`case(2)` lookup) | **dead runtime path** (swdrought=2 stub-errored on TOML) — compile-only confidence |

## 5. Sequencing — one dispatcher per commit, easy → hard

1. `irrigation` (near-pure deletion of dead stub)
2. `SSDI_irrigation` (same)
3. `csv_out` + `csv_out_tz` (output-only; inline into existing wrappers)
4. `Temperature`
5. `Solute`
6. `DoTillage` (resolve `case(4)`)
7. `SoilWater` (largest, 3 phases)
8. `SoilWaterStateVar` (save/restore)
9. `MatricFlux` (last; dead path)

Each step is its own commit so any regression bisects to a single dispatcher.

## 6. Verification & risk

- **Per-commit byte-identical regression** against the gfortran-4.2.0
  reference (`swap420gf`, now the regression truth) via the `check-fast`
  gate; `check-full` before finishing the arc. Each implementer/fix subagent
  ends with `check-fast` (per standing per-task regression-gate convention) —
  not deferred to an end-of-arc gate.

- **Init-order risk (#4–7).** Wrapping `state%X%init` inside the `x_init`
  facade and calling the facade at a single point *relocates* the type-bound
  init call from its current early position in `swap_init_body` to the old
  `X(1)` call site (or vice-versa). Each fold must verify:
  (a) nothing initialized between the two points is read by the relocated
  seed, and (b) nothing after the new position depends on the relocated
  type-bound init. **Fallback:** if a fold is unsafe, keep `x_init` doing
  *only* the seed and leave the type-bound `state%X%init` call where it is —
  the magic int is still eliminated, which is the actual goal.

- **MatricFlux (#9).** Dead at runtime, so regression cannot exercise it.
  Rely on compilation + careful structural split; flag the limited
  confidence explicitly in the commit.

- **State-schema rule.** If any `src/state/*_state.f90` is touched (none
  expected in this arc, but the type-bound init bodies live in state modules),
  `rm -rf builddir` before regression — incremental Meson builds don't
  propagate `.mod` deps across the swap_modern↔swap_legacy boundary.

- Subagent-driven implementers; development-only commits on `development`.

## 7. Deliverables

- This spec.
- **ADR 0043** — *Retire integer task-dispatch lifecycle multiplexing*:
  records the convention (init → `x_init` facade over type-bound `init`;
  step/finalize/sub-phases → named free procedures; option-switches
  explicitly excluded). Written alongside the first implementation commit.
- An implementation plan (via writing-plans) with 9 tasks following §5.

## 8. Success criteria

- All 9 dispatchers in §4 no longer take an integer task/iTask lifecycle
  selector; their bodies live in named free procedures.
- `swap_init_body` and `swap_run_step` contain no `X(1, …)` / `X(2, …)` /
  `X(3, …)` magic-int calls for the in-scope subsystems.
- `check-full` passes byte-identical against `swap420gf` across the
  regression suite (with the documented MatricFlux dead-path caveat).
