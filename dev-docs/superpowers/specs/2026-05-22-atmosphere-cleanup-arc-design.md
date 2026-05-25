# Atmosphere subsystem cleanup arc (GR-ATM-CLEAN) — design

**Date:** 2026-05-22
**Scope:** `src/atmosphere/` — five distinct refactors bundled into one sequential arc.
**Goal:** Bring the atmosphere subsystem to a modernised, single-responsibility, locals-by-default shape: kill the shared-state scratch bag, split bundled-module files, replace task-selector subroutines with named init/step pairs, split the daily/sub-daily orchestrators, hoist physical constants, and sweep F77 intrinsics. Each phase is independently bisectable.

## Background

The architecture review captured in [src/atmosphere/README.md](../../../src/atmosphere/README.md#cross-cutting-modernisation-themes) identified eight modernisation themes. This arc executes themes #1, #2, #3, #6 and #7 (themes #4 and #5 are handled separately — #4 awaits its own design discussion; #5 was shipped 2026-05-22 as the PenMon derived-type refactor).

The five items being addressed:

1. **Kill the shared-state scratch bag (`MeteoVars`)** — module-level scalars used as parameter channels between `ReadMeteoDay` and `ProcessMeteoDay`. Replace with locals or `state%atmosphere` fields.
2. **Split `meteoday.f90` and stop bundling modules** — four modules in one 1036-line file. Split per-module so each translation unit has one responsibility.
3. **Replace task-selector subroutines with named init/step pairs** — `snow(task, …)`, `CNmethod(Itask, …)`, `reduceva(task, …)` dispatch on an integer flag at every call. Split into named pairs.
6. **Separate daily vs sub-daily orchestrators** — `ProcessMeteoDay` has six scattered `if swmetdetail` branches. Promote to two top-level orchestrators sharing private helpers.
7. **Smaller efficiency wins worth bundling** — hoist physical constants to `atmosphere_constants_mod`; sweep `dexp`/`dlog`/`dmin1`/`dble` to F90 generics.

The PenMon refactor that just shipped (`c29408b`) is the template: pure-where-possible, explicit-interface, derived-type-bundled arguments. This arc applies the same discipline to the rest of the subsystem.

## Coordination with other arcs

- **GR-ATM Phase C3** (in-flight globals-retirement arc) owns the DEFERRED `use variables` imports inside `et.f90`, `meteodt.f90`, etc. This arc does **not** touch those imports.
- **sweep3 `state%legacy_t`** (introduced 2026-05-22 in `75df34c`) is a transitional bag for retiring `variables.f90` globals. `MeteoVars` scratch (this arc's #1) is a *separate* concern — its members are not `variables.f90` globals. No overlap.
- **PenMon derived-type refactor** (shipped 2026-05-22 in `c29408b`) is the immediate predecessor. This arc inherits the `pm_inputs_t` / `pm_outputs_t` pattern as the style template.

## Arc shape — six phases, sequential, dependency-driven

Each phase is a self-contained refactor. After every phase the build is green and `check-fast` (pFUnit + 4 regression cases) passes. Subagent-driven execution dispatches one subagent per phase.

| Phase | Item | What it does | Files touched | Est. commits |
|-------|------|--------------|---------------|--------------|
| **A** | #2  | Split `meteoday.f90`'s four modules into four files | `src/atmosphere/{meteo_vars.f90, runoff.f90, meteo_io.f90, meteo_orchestrator.f90}`; delete `meteoday.f90`; update `meson.build` | 4 |
| **B** | #7a | Populate `atmosphere_constants_mod`; hoist physical constants from `snow.f90`, `runoff.f90`, `et.f90` | `atmosphere_constants.f90`, `snow.f90`, `runoff.f90`, `et.f90` | 2 |
| **C** | #3  | Replace `snow`, `CNmethod`, `reduceva` task-selector subroutines with named init/step pairs; migrate callers | `snow.f90`, `runoff.f90`, `et.f90`, `meteo_orchestrator.f90`, `meteodt.f90`, `swap_mod.f90` | 3 |
| **D** | #1  | Retire `MeteoVars` — scratch scalars become locals; long-lived arrays migrate to `state%atmosphere` | `meteo_vars.f90` (deleted), `meteo_io.f90`, `meteo_orchestrator.f90`, `state/atmosphere_state.f90` | 2 |
| **E** | #6  | Split `ProcessMeteoDay` into two orchestrators sharing private helpers; dispatcher stays as the public entry point | `meteo_orchestrator.f90` (heavy surgery), no caller changes | 4 |
| **F** | #7b | Sweep `dexp`/`dlog`/`dmin1`/`dble` → F90 generics across the subsystem | All `src/atmosphere/*.f90` | 1 |

Total: ~16 commits. Verification gate at the end of every phase.

## Phase A — Split `meteoday.f90`

### Current state

`meteoday.f90` (1036 lines) bundles four modules:

- `module MeteoVars` (lines 17-47) — scratch scalars + arrays
- `module runoff_mod` (lines 64-218) — `CNmethod`
- `module meteo_process_mod` (lines 234-490) — `ReadMeteoDay`, `ResetMetFlx`
- `module meteo_mod` (lines 515-1036) — `ProcessMeteoDay`

### Target

| New file | Module | Contents |
|----------|--------|----------|
| `src/atmosphere/meteo_vars.f90` | `MeteoVars` | (Phase A: identical contents to current `module MeteoVars`. Phase D deletes this file.) |
| `src/atmosphere/runoff.f90` | `runoff_mod` | `CNmethod` |
| `src/atmosphere/meteo_io.f90` | `meteo_process_mod` | `ReadMeteoDay`, `ResetMetFlx` |
| `src/atmosphere/meteo_orchestrator.f90` | `meteo_mod` | `ProcessMeteoDay` |
| `src/atmosphere/meteoday.f90` | — | **deleted** |

### Scope

- **No code changes** to module bodies. Pure file relocation.
- **`meson.build` update:** drop `meteoday.f90` from `swap_atmosphere_sources`; add the four new files.
- **No `use` statement changes anywhere else** — module names are preserved. Existing `use meteo_mod`, `use meteo_process_mod`, `use runoff_mod`, `use MeteoVars` continue to resolve to the new files unchanged.

### Why this comes first

Every later phase needs to edit one of these four modules. Splitting first means edits land in focused files rather than churning the 1036-line bundle. Reduces merge surface for subsequent phases.

## Phase B — Constants hoist

### Current state

`atmosphere_constants.f90` is a header-only stub — comment block, no `module` declaration. Magic literals are scattered across:

- `snow.f90`: `cwat = 4180.0d0`, `lm = 333580d0`, `ts = 0.0d0`, `0.07` (liquid-water cap), `0.5` (soil-surface threshold)
- `runoff.f90` (after Phase A): `DEPTH_10CM = 10.0`, `H_FIELD_CAPACITY = -100.0`, `H_WILTING_POINT = -16000.0`, `IA_RATIO = 0.2` (currently declared as file-scope `parameter`s — promote to shared module)
- `et.f90` (`reduceva`): `POND_THRESHOLD = 1.0d-10`

### Target

Populate `atmosphere_constants.f90` with a `module atmosphere_constants_mod` exporting `real(real64), parameter` constants:

```fortran
module atmosphere_constants_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   public

   ! Thermodynamics
   real(real64), parameter :: SPECIFIC_HEAT_WATER   = 4180.0_real64    ! J/kg/K
   real(real64), parameter :: LATENT_HEAT_MELTING   = 333580.0_real64  ! J/kg
   real(real64), parameter :: SNOW_TEMPERATURE_C    = 0.0_real64       ! deg C

   ! Snow
   real(real64), parameter :: SNOW_LIQUID_WATER_FRACTION = 0.07_real64 ! fraction of total water storage
   real(real64), parameter :: SOIL_SURFACE_FREEZE_THRESHOLD_C = 0.5_real64

   ! Ponding
   real(real64), parameter :: POND_THRESHOLD_CM = 1.0e-10_real64

   ! Curve Number / runoff
   real(real64), parameter :: DEPTH_10CM_CM           = 10.0_real64
   real(real64), parameter :: H_FIELD_CAPACITY_CM     = -100.0_real64
   real(real64), parameter :: H_WILTING_POINT_CM      = -16000.0_real64
   real(real64), parameter :: INITIAL_ABSTRACTION_RATIO = 0.2_real64
end module atmosphere_constants_mod
```

### Scope

- **Do not move** constants already in `swap_constants_mod` (`KARMAN_CONSTANT`, `GRASS_HEIGHT_CM`, `MEASUREMENT_HEIGHT_CM`, `BARE_SOIL_HEIGHT_CM`, `PI`, `ALBEDO_PONDING`, `vlarge`, `small`, `nihil`). These stay where they are — `et.f90` already imports them from `swap_constants`.
- **Do not migrate** any constants from other subsystems (`src/soil/`, `src/crop/`).
- **Replace inline literals** in `snow.f90`, `runoff.f90`, `et.f90` with `use atmosphere_constants_mod, only: …` imports.
- **Verify byte-identical output** via check-fast — constants are just symbolic substitution.

## Phase C — Task selectors → named init/step pairs

### Current state

Three task-selector subroutines:

| Routine | Tasks | Callers |
|---------|-------|---------|
| `snow(task, state)` | 1=init, 2=step | `swap_mod.f90:432` (init), `swap_mod.f90:550` (step) |
| `CNmethod(Itask, state)` | 1=init, 2=step | `swap_mod.f90:384` (init), `meteo_orchestrator.f90:957` (step) |
| `reduceva(task, nrai, state)` | 1=daily, 2=dt | `meteo_orchestrator.f90:965` (daily), `meteodt.f90:429` (dt), `meteodt.f90:541` (dt) |

### Target

| Old | New |
|-----|-----|
| `snow(1, state)` | `snow_init(state)` |
| `snow(2, state)` | `snow_step(state)` |
| `CNmethod(1, state)` | `cn_init(state)` |
| `CNmethod(2, state)` | `cn_step(state)` |
| `reduceva(1, nrai, state)` | `reduceva_daily(nrai, state)` |
| `reduceva(2, nrai, state)` | `reduceva_dt(nrai, state)` |

### Scope

- For each routine: extract the `case (1)` and `case (2)` bodies into the two new public subroutines. Drop the `select case` dispatcher entirely. Drop the `call fatalerr_collected('X', 'Illegal Itask')` branches — those become unreachable.
- Update every caller to use the new name.
- Delete the old `snow`, `CNmethod`, `reduceva` public exports — no deprecation shims needed (7 call sites total across the codebase: snow×2, CNmethod×2, reduceva×3, all updated in this phase).
- Module names (`snow_mod`, `runoff_mod`, `et_mod`) stay the same.

### Why split

Each old routine bundles two distinct lifecycles (init runs once; step runs every timestep). Splitting (a) self-documents at the call site, (b) lets the compiler inline the hot path, (c) eliminates the runtime `select case` branch, (d) removes the now-dead "illegal task" fatalerr path.

## Phase D — Retire `MeteoVars`

### Current state

`module MeteoVars` (in `meteo_vars.f90` after Phase A) declares 24 module-level scalars and 2 arrays, used as parameter channels between `ReadMeteoDay` and `ProcessMeteoDay`:

```fortran
integer   count, first, i, irecord, last, ndayparts
real(8)   arain(96), awind(96), restint, interc, Edirectpond
real(8)   aintc, dttp, eintc, etr, gctp, hum, svp, wfrac, win
real(8)   netrainflux, rainflux, sumtav, Edirect, Tdirect, Tdirectwet
```

### Target

Categorise each member by lifetime and relocate:

| Category | Members | Relocation |
|----------|---------|------------|
| **Pure loop locals** | `count, first, i, irecord, last, ndayparts` | Inline into the routines that use them. Delete module declarations. (Note: `irecord` is currently a *module-level* variable that the `do 1000 irecord = 1, ndayparts` loop silently mutates — same name, same storage. Promoting it to a local in `ProcessMeteoDay` makes the per-call lifetime explicit.) |
| **Within-call scratch** | `interc, aintc, eintc, dttp, gctp, etr, hum, win, svp, wfrac, netrainflux, rainflux, sumtav, Edirect, Tdirect, Tdirectwet, Edirectpond` | Locals in `ProcessMeteoDay` (still a single routine at this point; Phase E later splits it). |
| **Cross-call persistent** | `restint, arain(96), awind(96)` | Migrate to `state%atmosphere%X` (new fields). `restint` persists across iterations of the `do 1000 irecord` loop and across days; `arain`/`awind` are populated in `ReadMeteoDay` and consumed in `ProcessMeteoDay`. |

After categorisation:

1. Add fields to `state%atmosphere` (in `src/state/atmosphere_state.f90`):
   - `real(real64) :: restint = 0.0_real64`
   - `real(real64) :: arain_subdaily(96) = 0.0_real64`
   - `real(real64) :: awind_subdaily(96) = 0.0_real64`

2. Rewrite all readers/writers of these three names to use `state%atmosphere%X`.

3. Rewrite `ReadMeteoDay` and `ProcessMeteoDay` to declare the scratch scalars as locals at the routine top.

4. Delete `meteo_vars.f90` entirely. Delete `use MeteoVars` from `meteo_io.f90` and `meteo_orchestrator.f90`.

5. Update `meson.build` to drop `meteo_vars.f90` from the sources list.

### Scope

- **No physics changes.** Pure scope migration.
- **State schema change:** three new fields in `state%atmosphere`. Per project convention, this requires a **clean rebuild** (`rm -rf builddir`) — flagged in the plan.
- **Naming choice for the migrated arrays:** `arain_subdaily` / `awind_subdaily` to disambiguate from `state%atmosphere%arai` (which is the per-day total rain in the daily-meteo path) and avoid collision.

## Phase E — Daily/sub-daily orchestrator split

### Current state

`ProcessMeteoDay` (~430 lines in `meteo_orchestrator.f90` after Phase A) has ten labelled sections, with **six** scattered `if (config%meteo%swmetdetail.eq.0) … elseif (config%meteo%swmetdetail.eq.1) …` branches:

- Section 4 (compute ET): daily branch uses `etr` directly; sub-daily branch calls `PenMon`
- Section 5 (Rutter interception): both branches but with different table-vs-runtime treatment
- Section 6 (wet fraction): different formulas
- Section 7 (peva/ptra): different handling of cover fraction
- Section 8 (per-record aggregation): sub-daily only
- Section 9 (daily finalisation): daily only — `graidt/nraidt/aintcdt`, `CNmethod`, `reduceva`, daily totals
- Section 10 (sub-daily finalisation): sub-daily only — `Tav`/`tmn`/`tmx` averaging, sub-daily fluxes

Plus an outer `do 1000 irecord = 1, ndayparts` loop where `ndayparts = 1` for daily and `nmetdetail` for sub-daily.

### Target

Three public-or-private routines + four private helpers in `meteo_orchestrator.f90`:

```fortran
! Public — preserves the legacy API for all existing callers
subroutine ProcessMeteoDay(state, config)
   if (config%meteo%swmetdetail == 0) then
      call process_meteo_day_daily(state, config)
   else
      call process_meteo_day_subdaily(state, config)
   end if
end subroutine

! Private — daily-mode path
subroutine process_meteo_day_daily(state, config)
   ! Section 3 (interception, single iteration)
   ! Section 4 (compute reference ET, daily branch — etr-direct or PenMon)
   ! Sections 6 + 7 (single iteration via helpers)
   ! Section 9 (daily finalisation — graidt/nraidt, CN, reduceva_daily, atmdem)
end subroutine

! Private — sub-daily-mode path
subroutine process_meteo_day_subdaily(state, config)
   ! Section 3 (interception, once per day)
   ! do irecord = 1, nmetdetail
   !   refresh per-record meteo
   !   Section 4 (compute reference ET, sub-daily branch — PenMon)
   !   Section 5 (Rutter)
   !   Sections 6 + 7 (per record via helpers)
   !   Section 8 (record aggregation)
   ! end do
   ! Section 10 (daily totals from sub-daily records)
end subroutine
```

Private helpers extracted (used by both orchestrators):

- `apply_interception_step(state, config, aintc)` — VonHHBraden / Gash / Rutter dispatch + DivIntercep
- `compute_reference_et(state, config, irecord, pmi, pmo)` — wraps PenMon pack/call/unpack; daily-vs-subdaily `swetr` handling lives inside
- `compute_wet_fraction(state, config, aintc, eintc, wfrac)` — Section 6 logic
- `partition_peva_ptra(state, config, wfrac)` — Section 7 logic (cover correction, ponding correction, PMdirect, CO2 correction). Writes `state%atmosphere%peva` and `%ptra` directly.

### Scope

- **No caller changes.** `ProcessMeteoDay` keeps its signature; it becomes a 5-line dispatcher.
- **Builds on Phase D.** Without `MeteoVars` retired, the orchestrator split would need to thread the scratch state through arguments — adding work to the split that doesn't belong there.
- **No physics changes.** The extracted helpers contain identical arithmetic to the inlined sections.
- **Caller chain unchanged:** `meteo_io.f90`, `meteodt.f90`, `swap_mod.f90` all continue to call `ProcessMeteoDay` exactly as today.

### Risk

This is the largest phase. Mitigation:

- Each commit extracts **one** helper at a time, verifying check-fast after each.
- Final commit removes the original inline code and switches to the helper-based version.
- Byte-identical regression outputs are the gate. Any drift → revert the offending commit, diff.

## Phase F — Intrinsic sweep

### Target

Mechanical search-and-replace across `src/atmosphere/*.f90`:

| Old | New |
|-----|-----|
| `dexp(x)` | `exp(x)` |
| `dlog(x)` | `log(x)` |
| `dmin1(a, b)` | `min(a, b)` |
| `dble(x)` | `real(x, real64)` (with `iso_fortran_env` import) |

### Scope

- Single commit covering all atmosphere files.
- `real(real64)` requires `iso_fortran_env` import — add where missing.
- check-fast must be **byte-identical** (these are aliases the compiler folds identically). Any drift = bug, revert and investigate.
- Do **not** also migrate `real(8)` → `real(real64)` in declarations. That's a separate sweep (different risk profile — `real(8)` is non-portable but byte-equivalent on every compiler we use).

## Out of scope

Deliberately not in this arc:

- **Item #4 (msw1eic / ruttervw quarantine):** native rewrite vs quarantine needs its own design discussion. Separate spec.
- **AFGEN cache in `Gash`:** perf optimisation, not a clarity refactor. Defer until profile evidence.
- **Dropping dual-write transitional blocks:** delicate per-field verification work owned by the state-rescue arc (GR-ATM Phase C3).
- **PenMon `pm_atmos_baseline_t` partial precomputation:** flagged in [PenMon spec](2026-05-22-penman-monteith-derived-types-design.md#forward-compatibility-with-partial-precomputation) as future arc.
- **DEFERRED `use variables` imports** in `et.f90`, `meteodt.f90` (e.g. `swredu`, `cofred`, `nird`, `finterception`, `dtEventRain`, `lat`): owned by GR-ATM Phase C3.
- **`real(8)` → `real(real64)` declaration sweep:** separate phase if/when it's prioritised.
- **Splitting `meteodt.f90`** (which also has multiple subroutines per module): focus stays on `meteoday.f90` for this arc.

## Verification

After **every** phase:

1. `pixi run build-linux` — clean build.
2. `pixi run test-pfunit` — all pFUnit tests pass.
3. `pixi run check-fast` — pFUnit + 4 regression cases (hupselbrook, surfacewater, salinitystress, grassgrowth) match golden outputs.

Phase D requires a clean rebuild (`rm -rf builddir`) because it adds fields to `state%atmosphere` — per the project's `feedback_state_schema_clean_rebuild` playbook entry, incremental Meson builds don't propagate `.mod` dependency changes across the `swap_modern` ↔ `swap_legacy` boundary on state-schema edits.

Phases A, B, C, E, F do **not** require clean rebuilds.

Regression case match must be byte-identical for Phases A, B, F (file moves, constant substitution, intrinsic aliasing — none affect arithmetic). Phases C, D, E can have tiny floating-point drift in principle (different temporary-variable evaluation order is possible if the compiler vectorises differently), but the project's regression harness asserts to single-precision tolerance; drift beyond that signals a real bug and the offending commit reverts.

## Definition of done

After all six phases:

- [ ] `meteoday.f90` does not exist.
- [ ] `MeteoVars` module does not exist; no atmosphere file imports it.
- [ ] `snow_mod`, `runoff_mod`, `et_mod` export `*_init` / `*_step` (or `*_daily` / `*_dt`) pairs; no `task` / `Itask` integer dispatchers remain.
- [ ] `ProcessMeteoDay` is a 5-line dispatcher; `process_meteo_day_daily` and `process_meteo_day_subdaily` are the two work-doing routines.
- [ ] `atmosphere_constants_mod` exports the constants enumerated above; no inline magic literals remain in `snow.f90`, `runoff.f90`, `et.f90` for those values.
- [ ] No `dexp` / `dlog` / `dmin1` / `dble` calls in `src/atmosphere/*.f90`.
- [ ] `pixi run check-fast` is green at every phase boundary and at arc end.
- [ ] No new DEFERRED imports introduced; no DEFERRED imports newly retired (that's a different arc).
