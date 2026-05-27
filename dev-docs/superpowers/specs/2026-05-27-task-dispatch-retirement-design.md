---
title: "Task-Dispatch Retirement — Tier 1 + Tier 2 (init/step/finalize de-multiplexing)"
date: 2026-05-27
status: complete
branch: development
scope: Tier 1 + Tier 2 dispatchers (crop cluster / Tier 3 deferred)
---

# Task-Dispatch Retirement — Tier 1 + Tier 2

> **COMPLETE (2026-05-27).** All 10 in-scope dispatchers retired over 11 commits
> (`6122a8d` ADR → `2e3768a`). `check-full`: 5 passed, 0 failed, 2 known xfails
> (`soilhysteresis`/`winter`, pre-existing adaptive-dt divergence per `cfc447e`)
> — zero new regressions; `oxygenstress` also byte-identical.
> **Note:** `SurfaceWater` was inadvertently omitted from the §4 table during
> design; it is a textbook Tier-1 dispatcher and was folded in as plan Task 9
> (`surfacewater_lateral` / `surfacewater_balance`, init stub deleted).

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

Reading the bodies revealed that initialization is genuinely **two phases**,
already latent in the code:

| Phase | What it does | Ordering | Mechanism |
|---|---|---|---|
| 1 — construct | allocate arrays, zero-fill, snapshot config scalars/tables into state | order-free (idempotent, takes no sibling state) | type-bound `state%X%init(cfg, …)` — **already exists** |
| 2 — seed | compute derived initial state from **sibling subsystems that must already be seeded** | ordered (soilwater → heat → solute) | the old `case (1)` body |

Evidence the seed is a distinct phase, not a duplicate init:
`Solute(1)` computes `cmsy`/`samini` from `soil%theta`, `soil%bdens`,
`mesh%z`; `Temperature(1)` sets the tsoil *profile* and `fquartz/fclay/forg`
from `soil%thetas`/`soil%orgmat` and is deliberately ordered *after*
`SoilWater(1)`; `heat_state_init` cannot do this (it takes no `soil` and its
`heat_cfg` arg is explicitly "reserved for future seed migration"). Folding
the seed into the type-bound init would either drop the derived computation
or break the ordering, so it stays a separate phase.

For each dispatcher `X(task, …)`:

- **Phase 1 stays put.** The type-bound `state%X%init(...)` call in
  `swap_init_body` is left exactly where it is. We *commit to it* as the
  canonical, order-free constructor.
- **init case → free `x_seed(state)`** (Phase 2): the old `case (1)` body,
  integer dropped, called at the *same site* `X(1)` is called today (so the
  cross-subsystem ordering is preserved byte-for-byte). It takes the whole
  `state` because it reads sibling subsystems.
- **step / sub-step cases → free `x_step(state)`** (plus `x_update`, etc. for
  genuinely distinct phases), integer dropped.
- **finalize case → free `x_finalize(state)`.**
- **save/restore & lazy-table (Tier 2) → named pairs**, not lifecycle verbs.

This is uniform across all subsystems: every one ends as
`state%X%init(...)` (Phase 1) + `x_seed(state)` (Phase 2) +
`x_step`/`x_update`/`x_finalize`. Stub-only inits (irrigation, SSDI) had
their Phase-1 work migrated to config-load time and have no Phase-2 seed, so
they contribute only a step procedure.

### Multi-instance invariant

In-process parallel multi-instance execution (multiple `swap_state_t` running
concurrently) is on the near-term roadmap. The decisive enabler is that
**every procedure takes `state` explicitly and touches nothing outside it —
no module-level variables, no `SAVE` locals** — so N instances each mutate
only their own `state` and the procedures are re-entrant. The Tier 1 + Tier 2
files already satisfy this (verified: no module vars, no `SAVE` in
`soilhydraulics`, `solute`, `temperature`, `tillage`, `irrigation`); the
de-multiplexing into `proc(state)` keeps them there. This invariant is
recorded in ADR 0043 and binds all future work.

> **Out-of-scope blocker (flagged, not fixed here):** the real obstacle to
> parallel in-process runs lives in the deferred crop cluster —
> `crop_config_global` (a module-global `pointer`, `crop/crop_config_global.f90:23`,
> reassigned on every `crop_state_init`, `crop_state.f90:85`) and `SAVE` locals
> in `cropgrowth.f90:601,678` / `cropgrass_runtime.f90:94`. These are tracked
> for the Tier 3 / globals effort, not this arc.

## 4. Per-dispatcher migration map

Procedure names below are proposals; final names may be adjusted during
implementation but must drop the integer and read as named phases.

### Tier 1 — pure lifecycle dispatch

Phase-1 (`state%X%init`) is unchanged in every row; the table lists the
Phase-2+ procedures extracted from the dispatcher.

| # | Dispatcher | File | New free procedures |
|---|---|---|---|
| 1 | `irrigation(task)` | `crop/irrigation.f90:31` | delete dead `case(1)` (no-op `return`); `case(2)` → `irrigation_step` (no seed — Phase-1 done at config load) |
| 2 | `SSDI_irrigation(iTask)` | `crop/irrigation.f90:~340` | delete dead `case(1)`; `case(2)` → `ssdi_irrigation_step` (no seed) |
| 3 | `csv_out` / `csv_out_tz(iTask)` | `io/csv_output.f90:644 / :24` | inline the three cases into the existing `csv_output_init` / `csv_output_step` / `csv_output_finalize` wrappers; delete both dispatchers |
| 4 | `Temperature(task)` | `heat/temperature.f90:82` | `temperature_seed` (`case(1)`), `temperature_step` (`case(2)`) |
| 5 | `Solute(task)` | `solute/solute.f90:10` | `solute_seed` (`case(1)`), `solute_step` (`case(2)`) |
| 6 | `DoTillage(iTask)` | `crop/tillage.f90:53` | `tillage_seed` (`case(1)`), `tillage_step` (`case(2)`), `tillage_output` (`case(3)`); **resolve `case(4)` liveness** (only 1/2/3 called from the orchestrator; `case(4)` body is empty `! CLOSURE`) |
| 7 | `SoilWater(task)` | `soil/soilhydraulics.f90:820` | `soilwater_seed` (`case(1)`), `soilwater_step` (Richards/headcalc, `case(2)`), `soilwater_update` (rate+state, `case(3)`) |
| 7b | `SurfaceWater(task)` | `drainage/surfacewater.f90:24` | delete `case(1)` init stub; `surfacewater_lateral` (`case(2)`), `surfacewater_balance` (`case(3)`); keeps `request_smaller_dt` out-arg *(added post-design as plan Task 9)* |

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

- **Init-order is preserved by construction.** Because Phase-1
  (`state%X%init`) is left untouched and each `x_seed` is called at the exact
  site its `X(1)` was, no init is relocated — the extraction is a pure rename
  + signature change, byte-identical by construction. The earlier
  "wrap-into-init" relocation risk is designed out. (Reason it would have been
  unsafe is documented in §3: e.g. `Temperature(1)` must run after
  `SoilWater(1)`, and soilwater arrays allocated early at `swap_mod.f90:122`
  are read before `:219`.)

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
- **ADR 0043** — *Retire integer task-dispatch lifecycle multiplexing*. Records:
  (a) the two-phase init model — type-bound `state%X%init` (construct) +
  free `x_seed(state)` (derived seed); (b) step/finalize/sub-phases → named
  free `proc(state)`; (c) the **multi-instance invariant** (procedures take
  `state`, touch nothing global); (d) option-switches explicitly excluded;
  (e) a forward note flagging `crop_config_global` + crop `SAVE` locals as the
  parallel-execution blocker, and a possible future top-level
  `state%init(config)` wrapping all phases. Written alongside the first
  implementation commit.
- An implementation plan (via writing-plans) with 9 tasks following §5.

## 8. Success criteria

- All 9 dispatchers in §4 no longer take an integer task/iTask lifecycle
  selector; their bodies live in named free procedures.
- `swap_init_body` and `swap_run_step` contain no `X(1, …)` / `X(2, …)` /
  `X(3, …)` magic-int calls for the in-scope subsystems.
- `check-full` passes byte-identical against `swap420gf` across the
  regression suite (with the documented MatricFlux dead-path caveat).
