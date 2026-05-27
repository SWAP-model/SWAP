# Task-Dispatch Retirement (Tier 1 + Tier 2) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Eliminate the legacy integer `select case (task)` lifecycle-dispatch anti-pattern from 10 independent compute subroutines, replacing magic-int calls with named free procedures, with byte-identical output preserved at every step.

**Architecture:** Two-phase init is made explicit: the order-free type-bound `state%X%init(...)` (construct) stays exactly where it is; each old `case(1)` body becomes a free `x_seed(state)` called at its current site (preserving cross-subsystem order); compute phases become `x_step`/`x_update`/etc. taking `state`. Every extracted procedure takes `state` and touches nothing global (multi-instance invariant). See `dev-docs/superpowers/specs/2026-05-27-task-dispatch-retirement-design.md`.

**Tech Stack:** Fortran 2008 (gfortran), Meson build via `pixi`, Python regression harness against the `swap420gf` reference binary.

---

## Refactoring discipline (read first)

This is a **behavior-preserving refactor**, not feature work. There is no new
behavior to test-drive; the standing invariant *is* the test:

```
pixi run -e test check-fast
```

(builds via `build-linux`, runs `test-pfunit`, then byte-compares 4 cases —
`hupselbrook surfacewater salinitystress grassgrowth` — against `swap420gf`.)

**Every task ends with `check-fast` passing before commit.** This is the
per-task regression gate; do not defer it to the end. The final task runs the
broader `check-full` (all six cases).

> **Known pre-existing xfails — do NOT treat as task regressions.** The
> `soilhysteresis` and `winter` fixtures currently fail against `swap420gf`
> due to an adaptive-dt desync (root-caused as *not physics* in commit
> `cfc447e`), unrelated to this arc. They are **not** in `check-fast`'s four
> cases (`hupselbrook surfacewater salinitystress grassgrowth`), so **every
> per-task `check-fast` must be fully clean (4/4)**. They WILL appear as
> failures in Task 12's `check-full` — that is expected and pre-existing. A
> task is only a regression if it breaks one of the four `check-fast` cases,
> or adds a *new* `check-full` failure beyond `soilhysteresis`/`winter`.

**Moving code:** when a step says "move the `case(N)` body", move the existing
lines **verbatim** (cut, don't retype) into the new procedure. Retyping a
200-line numerical body invites transcription drift and a regression you'll
spend an hour bisecting. The plan shows the new procedure *skeleton* and the
*call-site/`use`/`public`* edits — those are the parts that change.

**Per-task mechanical recipe** (every dispatcher task follows this shape):
1. In the owning module, add the new procedure(s); move each `case(N)` body into the matching one.
2. Delete the now-empty `select case` / dispatcher (or keep a thin private dispatcher only if other live callers still pass a task — none do here unless noted).
3. Update the module's `public ::` list: drop the old name, add the new ones.
4. In `src/core/swap_mod.f90`, update every `use <module>, only:` list (the names appear in both `swap_init_body`'s and `swap_run_step`'s `use` blocks) and replace each `call X(N, …)` with `call x_phase(…)`.
5. `pixi run -e test build-linux` → expect clean build.
6. `pixi run -e test check-fast` → expect all 4 cases byte-identical (PASS).
7. Commit.

**Branch:** `development`. Commits are development-only (do not push unless asked).

---

## Task 1: ADR 0043 — record the convention

**Files:**
- Create: `dev-docs/adr/0043-retire-task-dispatch.md`

- [ ] **Step 1: Write the ADR**

Use the house format (see `dev-docs/adr/0042-flatten-reset-cohorts.md`). Frontmatter `title`/`date: 2026-05-27`/`status: accepted`/`branch: development`, then sections **Context**, **Decision**, **Consequences**. The Decision section MUST state:

- The anti-pattern: a subroutine with an integer `task`/`iTask` argument whose `select case` branches are *lifecycle phases* (init/step/finalize), called from the orchestrator with magic ints.
- The replacement: two-phase init — type-bound `state%X%init(cfg, …)` (construct: allocate + zero + config snapshot; order-free) stays; the old `case(1)` body becomes a free `x_seed(state)` (derived state from sibling subsystems; called at its current site to preserve order). Compute phases → named free `proc(state)` (`x_step`, `x_update`, `x_finalize`, save/restore pairs).
- The **multi-instance invariant**: every extracted procedure takes `state` and touches nothing global — no module variables, no `SAVE` locals — so concurrent `swap_state_t` instances are re-entrant.
- Explicitly **excluded** (not the anti-pattern): genuine option-switches like `select case (model)` in `soil/WC_K_models_04_11.f90` and `reduceva_apply`'s daily/sub-daily mode.
- Forward note: the parallel-execution blocker is the deferred crop cluster — `crop_config_global` (`src/crop/crop_config_global.f90:23`) and `SAVE` locals in `cropgrowth.f90:601,678` / `cropgrass_runtime.f90:94`; a future top-level `state%init(config)` could wrap all phases as one per-instance constructor.

- [ ] **Step 2: Commit**

```bash
git add dev-docs/adr/0043-retire-task-dispatch.md
git commit -m "docs(adr): 0043 retire integer task-dispatch lifecycle multiplexing"
```

---

## Task 2: `irrigation` — delete dead init stub, rename step

**Files:**
- Modify: `src/crop/irrigation.f90` (subroutine `irrigation(task, state)` at `:31`; `public :: irrigation` at `:21`)
- Modify: `src/core/swap_mod.f90` (call at `:322` in `swap_run_step`; `use irrigation_mod` block)

`case(1)` is a dead no-op `return` (init done at config load); `case(2)` is the only live path.

- [ ] **Step 1: Replace the dispatcher with a plain step procedure**

New skeleton — drop the `task` arg and the `select case`, keep the `case(2)` body:

```fortran
subroutine irrigation_step(state)
   ! ... existing use/declares/associate from irrigation(task,state) ...
   type(swap_state_t), intent(inout) :: state
   ! <body of the old case(2)>  (move verbatim; drop the case(1) return stub)
end subroutine irrigation_step
```

- [ ] **Step 2: Update the public list**

`src/crop/irrigation.f90:21` — change `public :: irrigation` to `public :: irrigation_step`.

- [ ] **Step 3: Rewire the call site**

`src/core/swap_mod.f90` — in `swap_run_step`'s `use irrigation_mod, only:` list replace `irrigation` with `irrigation_step`; change `:322`
`if (time%flIrrigate) call irrigation(2, state)` → `if (time%flIrrigate) call irrigation_step(state)`.

- [ ] **Step 4: Build** — `pixi run -e test build-linux` → clean build.
- [ ] **Step 5: Regression** — `pixi run -e test check-fast` → 4/4 byte-identical.
- [ ] **Step 6: Commit**

```bash
git add src/crop/irrigation.f90 src/core/swap_mod.f90
git commit -m "refactor(crop): irrigation — drop task dispatch, irrigation_step(state)"
```

---

## Task 3: `SSDI_irrigation` — delete dead init stub, rename step

**Files:**
- Modify: `src/crop/irrigation.f90` (subroutine `SSDI_irrigation(iTask, state)` ~`:340`; `public :: SSDI_irrigation` at `:22`)
- Modify: `src/core/swap_mod.f90` (init call `:212` `SSDI_irrigation(1, …)`; step call `:393` `SSDI_irrigation(2, …)`; `use irrigation_mod` blocks)

`case(1)` is a dead stub (`return`); `case(2)` is the live "decide for next day" step.

- [ ] **Step 1: Replace the dispatcher**

```fortran
subroutine ssdi_irrigation_step(state)
   type(swap_state_t), intent(inout) :: state
   ! <body of the old case(2)>  (move verbatim)
end subroutine ssdi_irrigation_step
```

- [ ] **Step 2: Update the public list** — `:22` `public :: SSDI_irrigation` → `public :: ssdi_irrigation_step`.

- [ ] **Step 3: Rewire call sites in `swap_mod.f90`**

- `use irrigation_mod, only:` lists: replace `SSDI_irrigation` with `ssdi_irrigation_step`.
- `swap_init_body` `:212`: **delete** the line `if (state%cfg%irrigation%swssdi == 1) call SSDI_irrigation(1, state)` (dead init stub).
- `swap_run_step` `:393`: `if (state%cfg%irrigation%swssdi == 1) call SSDI_irrigation(2, state)` → `… call ssdi_irrigation_step(state)`.

- [ ] **Step 4: Build** — `pixi run -e test build-linux`.
- [ ] **Step 5: Regression** — `pixi run -e test check-fast` → 4/4 byte-identical.
- [ ] **Step 6: Commit**

```bash
git add src/crop/irrigation.f90 src/core/swap_mod.f90
git commit -m "refactor(crop): SSDI_irrigation — drop task dispatch, ssdi_irrigation_step(state)"
```

---

## Task 4: `csv_out` / `csv_out_tz` — split into header/write/close

**Files:**
- Modify: `src/io/csv_output.f90` (module `csv_output`: `csv_out(iTask, state)` at `:644`; wrappers `csv_output_init/step/finalize` at `:1237–:1263`)
- Modify: the file/module defining `csv_out_tz(iTask, state)` (imported via `use csv_output_tz` at `csv_output.f90:1238`; confirm its path with `grep -rn "subroutine csv_out_tz" src/io`)

The public wrappers `csv_output_init/step/finalize` already exist and are the canonical API — the magic-int lives one layer down. Split each dispatcher into three private procedures and point the wrappers at them.

- [ ] **Step 1: Split `csv_out`**

In module `csv_output`, replace `subroutine csv_out(iTask, state)` with three private procedures, moving the case bodies verbatim:

```fortran
private :: csv_out_header, csv_out_write, csv_out_close

subroutine csv_out_header(state)   ! <body of csv_out case(1)>
subroutine csv_out_write(state)    ! <body of csv_out case(2)>
subroutine csv_out_close(state)    ! <body of csv_out case(3)>
```

- [ ] **Step 2: Split `csv_out_tz` the same way** in its module: `csv_out_tz_header` / `csv_out_tz_write` / `csv_out_tz_close`.

- [ ] **Step 3: Repoint the wrappers**

`csv_output.f90:1240–1262` — replace the `call csv_out(N, state)` / `call csv_out_tz(N, state)` lines:

```fortran
! csv_output_init
if (state%cfg%output_csv%enabled    == 1) call csv_out_header(state)
if (state%cfg%output_csv%enabled_tz == 1) call csv_out_tz_header(state)
! csv_output_step  (keep the existing flYearStart flush block unchanged)
if (state%cfg%output_csv%enabled    == 1) call csv_out_write(state)
if (state%cfg%output_csv%enabled_tz == 1) call csv_out_tz_write(state)
! csv_output_finalize
if (state%cfg%output_csv%enabled    == 1) call csv_out_close(state)
if (state%cfg%output_csv%enabled_tz == 1) call csv_out_tz_close(state)
```

Update the `use csv_output_tz, only:` imports in the three wrappers to the new `_header/_write/_close` names.

- [ ] **Step 4: Build** — `pixi run -e test build-linux`.
- [ ] **Step 5: Regression** — `pixi run -e test check-fast` → 4/4 byte-identical (this task changes only output formatting paths; CSV bytes must match).
- [ ] **Step 6: Commit**

```bash
git add src/io/csv_output.f90 src/io/<csv_output_tz file>
git commit -m "refactor(io): split csv_out/csv_out_tz iTask dispatch into header/write/close"
```

---

## Task 5: `Temperature` — seed + step

**Files:**
- Modify: `src/heat/temperature.f90` (`subroutine temperature(task, state, config)` at `:82`; `public :: temperature, devries` at `:41`; cases at `:116`/`:156`)
- Modify: `src/core/swap_mod.f90` (init `:255` `Temperature(1,…)`; step `:373` `Temperature(2,…)`; `use temperature_mod` blocks)

Phase-1 (`state%heat%init`) at `swap_mod.f90:216` stays. `case(1)` is the Phase-2 seed (sets tsoil profile + fquartz/fclay/forg from soilwater — must stay after `SoilWater(1)`). `case(2)` is the step.

- [ ] **Step 1: Extract two procedures** (note `temperature` takes `config`):

```fortran
subroutine temperature_seed(state, config)
   type(swap_state_t),  intent(inout) :: state
   type(swap_config_t), intent(in)    :: config
   ! ... existing use/declares/associate ...
   ! <body of case(1)>  (move verbatim; keep the trailing `return`'s logic but no select)
end subroutine temperature_seed

subroutine temperature_step(state, config)
   type(swap_state_t),  intent(inout) :: state
   type(swap_config_t), intent(in)    :: config
   ! <body of case(2)>  (move verbatim)
end subroutine temperature_step
```

- [ ] **Step 2: Public list** — `:41` `public :: temperature, devries` → `public :: temperature_seed, temperature_step, devries`.

- [ ] **Step 3: Rewire `swap_mod.f90`**

- `use temperature_mod, only:` lists: replace `Temperature` with `temperature_seed, temperature_step`.
- `:255`: `if (time%flTemperature) call Temperature(1, state, config)` → `… call temperature_seed(state, config)`.
- `:373`: `if (time%flTemperature) call Temperature(2, state, config)` → `… call temperature_step(state, config)`.

- [ ] **Step 4: Build** — `pixi run -e test build-linux`.
- [ ] **Step 5: Regression** — `pixi run -e test check-fast` → 4/4 byte-identical.
- [ ] **Step 6: Commit**

```bash
git add src/heat/temperature.f90 src/core/swap_mod.f90
git commit -m "refactor(heat): temperature — drop task dispatch, temperature_seed/step(state,config)"
```

---

## Task 6: `Solute` — seed + step

**Files:**
- Modify: `src/solute/solute.f90` (`subroutine solute(task, state)` at `:10`; `public :: solute` at `:6`; cases `:45`/`:73`)
- Modify: `src/core/swap_mod.f90` (init `:268` `Solute(1,…)`; step `:374` `Solute(2,…)`; `use solute_mod` blocks)

Phase-1 (`state%solute%init`) at `swap_mod.f90:245` stays. `case(1)` is the Phase-2 seed (computes `cmsy`/`samini` from `soil%theta`/`soil%bdens`/`mesh%z` — must stay after `SoilWater(1)`).

- [ ] **Step 1: Extract two procedures**

```fortran
subroutine solute_seed(state)
   type(swap_state_t), intent(inout) :: state
   ! ... existing use/declares/associate ...
   ! <body of case(1)>  (move verbatim)
end subroutine solute_seed

subroutine solute_step(state)
   type(swap_state_t), intent(inout) :: state
   ! <body of case(2)>  (move verbatim)
end subroutine solute_step
```

- [ ] **Step 2: Public list** — `:6` `public :: solute` → `public :: solute_seed, solute_step`.

- [ ] **Step 3: Rewire `swap_mod.f90`**

- `use solute_mod, only:` lists: replace `solute` with `solute_seed, solute_step`.
- `:268`: `if (time%flSolute) call Solute(1, state)` → `… call solute_seed(state)`.
- `:374`: `if (time%flSolute) call Solute(2, state)` → `… call solute_step(state)`.

- [ ] **Step 4: Build** — `pixi run -e test build-linux`.
- [ ] **Step 5: Regression** — `pixi run -e test check-fast` → 4/4 byte-identical.
- [ ] **Step 6: Commit**

```bash
git add src/solute/solute.f90 src/core/swap_mod.f90
git commit -m "refactor(solute): solute — drop task dispatch, solute_seed/step(state)"
```

---

## Task 7: `DoTillage` — seed + step + output (resolve case 4)

**Files:**
- Modify: `src/crop/tillage.f90` (`subroutine DoTillage(iTask, state)` at `:29`; `public :: DoTillage` at `:25`; cases `:54`/`:92`/`:154`/`:179`; guard `if (iTask > 1 .and. cfg_soil%swtill == 0) return` at `:45`)
- Modify: `src/core/swap_mod.f90` (init `:211` `DoTillage(1,…)`; step `:324` `DoTillage(2,…)`; output `:400` `DoTillage(3,…)`; `use tillage_mod` blocks)

Phase-1 (`state%tillage%init`) at `swap_mod.f90:208` stays. `case(4)` (`! CLOSURE`) has an **empty body and no caller** — drop it. The `case default` fatalerr also goes away with the `select`. The `:45` guard (`iTask > 1 .and. swtill == 0`) protected the step/output paths — fold the `swtill == 0` early-return into `tillage_step`/`tillage_output` (they are already gated by `if (state%cfg%soil%swtill == 1)` at the call sites, so the guard is redundant at the call sites but keep an internal guard for safety parity).

- [ ] **Step 1: Extract three procedures**

```fortran
subroutine tillage_seed(state)     ! <body of case(1)>  — validation/setup
subroutine tillage_step(state)     ! <body of case(2)>
subroutine tillage_output(state)   ! <body of case(3)>  — debug writes (headless-gated)
```

Each: `type(swap_state_t), intent(inout) :: state` plus the existing `use`/declares/`associate (… )`. Delete `case(4)` and `case default`.

- [ ] **Step 2: Public list** — `:25` `public :: DoTillage` → `public :: tillage_seed, tillage_step, tillage_output`.

- [ ] **Step 3: Rewire `swap_mod.f90`**

- `use tillage_mod, only:` lists: replace `DoTillage` with `tillage_seed, tillage_step, tillage_output`.
- `:211`: `if (state%cfg%soil%swtill == 1) call DoTillage(1, state)` → `… call tillage_seed(state)`.
- `:324`: `if (state%cfg%soil%swtill == 1) call DoTillage(2, state)` → `… call tillage_step(state)`.
- `:400`: `if (state%cfg%soil%swtill == 1) call DoTillage(3, state)` → `… call tillage_output(state)`.

- [ ] **Step 4: Build** — `pixi run -e test build-linux`.
- [ ] **Step 5: Regression** — `pixi run -e test check-fast` → 4/4 byte-identical. (Default cases run with `swtill == 0`, so tillage paths are inactive — but the build must still link and the gated paths must compile.)
- [ ] **Step 6: Commit**

```bash
git add src/crop/tillage.f90 src/core/swap_mod.f90
git commit -m "refactor(crop): DoTillage — drop iTask dispatch (seed/step/output); delete dead case 4"
```

---

## Task 8: `SoilWater` — seed + step + update

**Files:**
- Modify: `src/soil/soilhydraulics.f90` (`subroutine soilwater(task, state)` at `:820`; `public :: … soilwater …` at `:12`; cases `:868`/`:1112`/`:1151`)
- Modify: `src/core/swap_mod.f90` (init `:219` `SoilWater(1,…)`; step `:358` `SoilWater(2,…)`; update `:371` `SoilWater(3,…)`; `use soilhydraulics_mod` blocks)

Phase-1 (`state%soilwater%init`) at `swap_mod.f90:122` stays. `case(1)` is the large Phase-2 seed; `case(2)` is the Richards/headcalc step; `case(3)` is the rate+state update.

- [ ] **Step 1: Extract three procedures**

```fortran
subroutine soilwater_seed(state)     ! <body of case(1)>  (large; move verbatim)
subroutine soilwater_step(state)     ! <body of case(2)>  (Richards/headcalc)
subroutine soilwater_update(state)   ! <body of case(3)>  (rate + state)
```

Each carries the full existing `use` list, declares, and the big `associate (mesh => …, soil => …, …)` block (replicate the associate in each — it aliases sub-states the body uses). `type(swap_state_t), intent(inout) :: state`.

- [ ] **Step 2: Public list** — `:12` replace `soilwater` with `soilwater_seed, soilwater_step, soilwater_update` (keep `headcalc, soilwaterstatevar, hysteresis`).

- [ ] **Step 3: Rewire `swap_mod.f90`**

- `use soilhydraulics_mod, only:` lists: replace `soilwater` with `soilwater_seed, soilwater_step, soilwater_update` (keep `SoilWaterStateVar` for now — Task 10).
- `:219`: `call SoilWater(1, state)` → `call soilwater_seed(state)`.
- `:358`: `if (.not.time%fldecdt) call SoilWater(2, state)` → `… call soilwater_step(state)`.
- `:371`: `call SoilWater(3, state)` → `call soilwater_update(state)`.

- [ ] **Step 4: Build** — `pixi run -e test build-linux`.

  > **If you touched a `src/state/*_state.f90` (you should not have for this task), `rm -rf builddir` first** — incremental Meson does not propagate `.mod` deps across the swap_modern↔swap_legacy boundary.

- [ ] **Step 5: Regression** — `pixi run -e test check-fast` → 4/4 byte-identical. (This is the highest-traffic numerical path; a non-identical result means a body was split mid-state — re-check the `associate` aliases.)
- [ ] **Step 6: Commit**

```bash
git add src/soil/soilhydraulics.f90 src/core/swap_mod.f90
git commit -m "refactor(soil): soilwater — drop task dispatch (seed/step/update)"
```

---

## Task 9: `SurfaceWater` — lateral + balance (delete init stub)

**Files:**
- Modify: `src/drainage/surfacewater.f90` (`subroutine SurfaceWater(task, state, request_smaller_dt)` at `:24`; `public :: SurfaceWater, surfacewater_year_reset` at `:21`; cases `:60`/`:66`/`:154`)
- Modify: `src/core/swap_mod.f90` (init stub `:252` `SurfaceWater(1,…)`; lateral `:353` `SurfaceWater(2,…)`; balance `:360` `SurfaceWater(3,…)`; `use surfacewater_mod` blocks)

Phase-1 (`state%surfacewater%init`) at `swap_mod.f90:248` stays. `case(1)` is a no-op stub (init hoisted) — delete it. `case(2)` = lateral drainage fluxes, `case(3)` = surface-water balance. Both produce the `request_smaller_dt` out-arg that drives the dt-reduction loop in `swap_run_step` — preserve that argument.

- [ ] **Step 1: Extract two procedures** (keep the `request_smaller_dt` out-arg and the `request_smaller_dt = .false.` reset at entry):

```fortran
subroutine surfacewater_lateral(state, request_smaller_dt)
   type(swap_state_t), intent(inout) :: state
   logical,            intent(out)   :: request_smaller_dt
   ! ... existing use/declares/associate ...
   request_smaller_dt = .false.
   ! <body of case(2)>  (move verbatim)
end subroutine surfacewater_lateral

subroutine surfacewater_balance(state, request_smaller_dt)
   type(swap_state_t), intent(inout) :: state
   logical,            intent(out)   :: request_smaller_dt
   request_smaller_dt = .false.
   ! <body of case(3)>  (move verbatim)
end subroutine surfacewater_balance
```

- [ ] **Step 2: Public list** — `:21` `public :: SurfaceWater, surfacewater_year_reset` → `public :: surfacewater_lateral, surfacewater_balance, surfacewater_year_reset`.

- [ ] **Step 3: Rewire `swap_mod.f90`**

- `use surfacewater_mod, only:` lists: replace `SurfaceWater` with `surfacewater_lateral, surfacewater_balance` (keep `surfacewater_year_reset`).
- `:252`: **delete** `if (time%flSurfaceWater) call SurfaceWater(1, state, request_smaller_dt)` (dead stub).
- `:353`: `… call SurfaceWater(2, state, request_smaller_dt)` → `… call surfacewater_lateral(state, request_smaller_dt)`.
- `:360`: `… call SurfaceWater(3, state, request_smaller_dt)` → `… call surfacewater_balance(state, request_smaller_dt)`.

- [ ] **Step 4: Build** — `pixi run -e test build-linux`.
- [ ] **Step 5: Regression** — `pixi run -e test check-fast` → 4/4 byte-identical (the `surfacewater` case exercises this directly).
- [ ] **Step 6: Commit**

```bash
git add src/drainage/surfacewater.f90 src/core/swap_mod.f90
git commit -m "refactor(drainage): SurfaceWater — drop task dispatch (lateral/balance); delete init stub"
```

---

## Task 10: `SoilWaterStateVar` — save + restore

**Files:**
- Modify: `src/soil/soilhydraulics.f90` (`subroutine SoilWaterStateVar(task, state)` at `:1203`; `public :: … soilwaterstatevar …` at `:12`; cases `:1215`/case(2))
- Modify: internal caller `src/soil/soilhydraulics.f90:1144` `call SoilWaterStateVar(1, state)`
- Modify: `src/core/swap_mod.f90` (`:364` `SoilWaterStateVar(2, state)`; `use soilhydraulics_mod` block)

`case(1)` saves the snapshot at t (`hm1`, `thetm1`, `gwlm1`, `pondm1`); `case(2)` restores on dt-reduce. Not lifecycle — name them honestly.

- [ ] **Step 1: Extract two procedures**

```fortran
subroutine soilwater_save_state(state)     ! <body of case(1)>
subroutine soilwater_restore_state(state)  ! <body of case(2)>
```

Each: `type(swap_state_t), intent(inout) :: state` + the existing `associate (mesh => …, soil => …)`.

- [ ] **Step 2: Public list** — `:12` replace `soilwaterstatevar` with `soilwater_save_state, soilwater_restore_state`.

- [ ] **Step 3: Rewire call sites**

- `src/soil/soilhydraulics.f90:1144`: `call SoilWaterStateVar(1, state)` → `call soilwater_save_state(state)` (intra-module; ensure both new procs precede/are accessible).
- `src/core/swap_mod.f90` `use soilhydraulics_mod, only:` list: replace `SoilWaterStateVar` with `soilwater_save_state, soilwater_restore_state`.
- `:364`: `call SoilWaterStateVar(2, state)` → `call soilwater_restore_state(state)`.

- [ ] **Step 4: Build** — `pixi run -e test build-linux`.
- [ ] **Step 5: Regression** — `pixi run -e test check-fast` → 4/4 byte-identical (the dt-reduction path runs in all cases).
- [ ] **Step 6: Commit**

```bash
git add src/soil/soilhydraulics.f90 src/core/swap_mod.f90
git commit -m "refactor(soil): SoilWaterStateVar — split into soilwater_save_state/restore_state"
```

---

## Task 11: `MatricFlux` — build_table + lookup (dead path)

**Files:**
- Modify: `src/crop/rootextraction.f90` (`subroutine MatricFlux(task,phead,node,outcome,state)` at `:329`; `public :: RootExtraction, MatricFlux` at `:13`; cases `:349`/case(2))
- Modify: caller `src/crop/cropgrowth.f90:302` `call MatricFlux(1, …)`
- Modify: dormant callers `src/crop/dormant/jongvanlier.f90:183,503,504,527` `call MatricFlux(2, …)`

**Caveat:** this path is **dead at runtime** (reached only when `swdrought=2`, which is stub-errored on the TOML path). Regression cannot exercise it — confidence is **compile + structural** only. State this in the commit body.

`case(1)` builds `state%soilwater%mfluxtable`; `case(2)` is the per-call lookup writing `outcome`.

- [ ] **Step 1: Extract two procedures**

```fortran
subroutine matricflux_build_table(state)
   type(swap_state_t), intent(inout) :: state
   ! <body of case(1)>  (move verbatim)
end subroutine matricflux_build_table

subroutine matric_flux(phead, node, outcome, state)
   real(8),            intent(in)    :: phead
   integer,            intent(in)    :: node
   real(8),            intent(out)   :: outcome
   type(swap_state_t), intent(inout) :: state
   ! <body of case(2)>  (move verbatim)
end subroutine matric_flux
```

(Drop the `optional` on `state` — both callers pass it.)

- [ ] **Step 2: Public list** — `:13` `public :: RootExtraction, MatricFlux` → `public :: RootExtraction, matricflux_build_table, matric_flux`.

- [ ] **Step 3: Rewire callers**

- `src/crop/cropgrowth.f90:302`: `call MatricFlux(1, state%soilwater%h(1), 1, dummy_mf_, state)` → `call matricflux_build_table(state)`.
- `src/crop/dormant/jongvanlier.f90` lines `:183,:503,:504,:527`: `call MatricFlux(2, <phead>, <node>, <outcome>, state)` → `call matric_flux(<phead>, <node>, <outcome>, state)`. Update jongvanlier's `use rootextraction_mod, only:` import (`MatricFlux` → `matric_flux`).

- [ ] **Step 4: Build** — `pixi run -e test build-linux` → clean build (this is the primary safety check for this task).
- [ ] **Step 5: Regression** — `pixi run -e test check-fast` → 4/4 byte-identical (confirms the dead path didn't accidentally activate; it should not change anything).
- [ ] **Step 6: Commit**

```bash
git add src/crop/rootextraction.f90 src/crop/cropgrowth.f90 src/crop/dormant/jongvanlier.f90
git commit -m "refactor(crop): MatricFlux — split into matricflux_build_table/matric_flux (dead path; compile-verified)"
```

---

## Task 12: Arc close-out — full regression + spec update

**Files:**
- Modify: `dev-docs/superpowers/specs/2026-05-27-task-dispatch-retirement-design.md` (mark status complete; note SurfaceWater added)

- [ ] **Step 1: Full regression gate**

Run: `pixi run -e test check-full`
Expected: the four `check-fast` cases byte-identical against `swap420gf`;
pFUnit `OK (N tests)`. **`soilhysteresis` and `winter` are expected to fail
(known pre-existing adaptive-dt xfails, commit `cfc447e`) — confirm those are
the *only* failures and that this arc introduced no new ones.**

- [ ] **Step 2: Confirm no Tier-1/2 task dispatch remains**

```bash
grep -rnE "select case *\( *(task|iTask) *\)" src/soil/soilhydraulics.f90 src/heat/temperature.f90 \
  src/solute/solute.f90 src/crop/tillage.f90 src/crop/irrigation.f90 src/drainage/surfacewater.f90 \
  src/crop/rootextraction.f90 src/io/csv_output.f90
```
Expected: **no matches**.

```bash
grep -rnE "call (SoilWater|Temperature|Solute|DoTillage|SurfaceWater|irrigation|SSDI_irrigation|SoilWaterStateVar|MatricFlux) *\( *[0-9]" src/core/swap_mod.f90
```
Expected: **no matches**.

- [ ] **Step 3: Update spec status to complete** (front-matter `status: complete`, add a line noting SurfaceWater was folded in as task 9).

- [ ] **Step 4: Commit**

```bash
git add dev-docs/superpowers/specs/2026-05-27-task-dispatch-retirement-design.md
git commit -m "docs(spec): mark task-dispatch retirement (Tier 1 + Tier 2) complete"
```

---

## Self-review notes (for the executor)

- **Spec coverage:** §4 dispatchers map to Tasks 2–11; SurfaceWater (omitted from the spec table) is Task 9 and must be confirmed in-scope with the user before starting. ADR (§7) is Task 1. Final `check-full` (§6/§8) is Task 12.
- **Naming consistency:** `*_seed` (Phase-2), `*_step` (per-step), `*_update`/`*_output`/`*_lateral`/`*_balance` (distinct phases), `soilwater_save_state`/`soilwater_restore_state` and `matricflux_build_table`/`matric_flux` (Tier-2 pairs). These names are used identically in the module, the `public` list, and the `swap_mod.f90` call sites within each task.
- **Order safety:** no `state%X%init` call is moved; every `*_seed` is called at its old `X(1)` site. Init order is preserved by construction (spec §3, §6).
- **Line numbers** are anchors from 2026-05-27; they drift as edits land. If a line doesn't match, re-grep for the `call X(N` pattern named in the step.
