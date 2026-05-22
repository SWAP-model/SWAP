# Atmosphere subsystem cleanup arc (GR-ATM-CLEAN) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Dispatch one subagent per **phase** (A–F). Each phase ends with a `pixi run check-fast` verification gate; do not start the next phase until the gate is green. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Bring the `src/atmosphere/` subsystem to a modernised, single-responsibility, locals-by-default shape via six dependency-ordered phases.

**Architecture:** Six sequential phases — Phase A splits `meteoday.f90`'s four modules into four files; Phase B hoists physical constants to `atmosphere_constants_mod`; Phase C replaces three task-selector subroutines with named init/step pairs; Phase D retires the `MeteoVars` scratch bag; Phase E splits `ProcessMeteoDay` into daily and sub-daily orchestrators sharing private helpers; Phase F sweeps F77 intrinsics to F90 generics. Verification gate after each phase. Each phase is independently bisectable.

**Tech Stack:** Fortran 2003+; Meson; pFUnit; pixi task runner.

**Spec:** [docs/superpowers/specs/2026-05-22-atmosphere-cleanup-arc-design.md](../specs/2026-05-22-atmosphere-cleanup-arc-design.md)

---

## Pre-arc setup

Before dispatching the first phase, confirm:

- [ ] Working tree on a branch dedicated to this arc (or on `development` if that's the project's convention — check recent commit history for the established workflow).
- [ ] `pixi run check-fast` is currently green from HEAD. If not, fix that first — do **not** start this arc on a broken baseline.
- [ ] The PenMon refactor (commit `c29408b` or later) is in the history. This arc depends on `et_mod` exporting `pm_inputs_t`/`pm_outputs_t`.

---

# Phase A — Split `meteoday.f90`

**Goal:** Move each of the four modules in `meteoday.f90` into its own file. No code changes inside the modules. After this phase, every later phase edits a focused file rather than the 1036-line bundle.

**Files at end of phase:**
- Create: `src/atmosphere/meteo_vars.f90`
- Create: `src/atmosphere/runoff.f90`
- Create: `src/atmosphere/meteo_io.f90`
- Create: `src/atmosphere/meteo_orchestrator.f90`
- Modify: `meson.build`
- Delete: `src/atmosphere/meteoday.f90`

**No `use` statement changes elsewhere** — module names (`MeteoVars`, `runoff_mod`, `meteo_process_mod`, `meteo_mod`) are preserved, so existing imports resolve to the new files automatically.

## Task A.1: Extract `module MeteoVars` → `meteo_vars.f90`

- [ ] **Step A.1.1: Identify the module bounds in source**

Read `src/atmosphere/meteoday.f90` and locate `module MeteoVars` (currently lines 17–47). Copy the entire module block (from `module MeteoVars` through `end module MeteoVars`, inclusive) plus the doc-comment block immediately above it.

- [ ] **Step A.1.2: Create `src/atmosphere/meteo_vars.f90`**

Write the new file. Contents:
- The doc-comment block from above `module MeteoVars` (if present).
- The complete `module MeteoVars … end module MeteoVars` block, **byte-identical** to what's in `meteoday.f90`.

- [ ] **Step A.1.3: Add `meteo_vars.f90` to the atmosphere sources in `meson.build`**

Find the atmosphere sources list (`grep -n "src/atmosphere" meson.build` to locate). Add `'src/atmosphere/meteo_vars.f90'` to the list. Leave `'src/atmosphere/meteoday.f90'` in place — it will be removed in Task A.4 after all four modules are extracted.

- [ ] **Step A.1.4: Verify the build still passes**

Run: `pixi run build-linux 2>&1 | tail -20`
Expected: clean build. (We now have two definitions of `module MeteoVars` — one in the original `meteoday.f90` and one in the new file. This will fail with a duplicate-symbol error.)

**If duplicate-symbol error:** Remove the `module MeteoVars` block from `meteoday.f90` (delete lines that previously held the module) — keep the other three modules intact. Re-run build. Expected: clean build now.

- [ ] **Step A.1.5: Run check-fast partial verification (build + pFUnit)**

Run: `pixi run test-pfunit 2>&1 | tail -10`
Expected: all tests pass.

- [ ] **Step A.1.6: Commit**

```bash
git add src/atmosphere/meteo_vars.f90 src/atmosphere/meteoday.f90 meson.build
git commit -m "$(cat <<'EOF'
refactor(atmosphere): extract MeteoVars into meteo_vars.f90

Phase A.1 of GR-ATM-CLEAN: split meteoday.f90's four modules into
four files. Module body byte-identical to original.

Spec: docs/superpowers/specs/2026-05-22-atmosphere-cleanup-arc-design.md
EOF
)"
```

## Task A.2: Extract `module runoff_mod` → `runoff.f90`

- [ ] **Step A.2.1: Copy `module runoff_mod` block from `meteoday.f90`**

In `meteoday.f90`, the `runoff_mod` block now starts at the line that begins `module runoff_mod` and ends at `end module runoff_mod` (lines have shifted after Task A.1; locate via `grep`).

- [ ] **Step A.2.2: Create `src/atmosphere/runoff.f90`**

Write the file containing the byte-identical `module runoff_mod … end module runoff_mod` block plus any doc comment immediately preceding it.

- [ ] **Step A.2.3: Add to `meson.build`**

Add `'src/atmosphere/runoff.f90'` to the atmosphere sources list.

- [ ] **Step A.2.4: Remove the `module runoff_mod` block from `meteoday.f90`**

Delete the entire `module runoff_mod … end module runoff_mod` block from `meteoday.f90`.

- [ ] **Step A.2.5: Build + test**

Run: `pixi run build-linux 2>&1 | tail -10`
Expected: clean.
Run: `pixi run test-pfunit 2>&1 | tail -5`
Expected: tests pass.

- [ ] **Step A.2.6: Commit**

```bash
git add src/atmosphere/runoff.f90 src/atmosphere/meteoday.f90 meson.build
git commit -m "$(cat <<'EOF'
refactor(atmosphere): extract runoff_mod (CNmethod) into runoff.f90

Phase A.2 of GR-ATM-CLEAN. Module body byte-identical to original.

Spec: docs/superpowers/specs/2026-05-22-atmosphere-cleanup-arc-design.md
EOF
)"
```

## Task A.3: Extract `module meteo_process_mod` → `meteo_io.f90`

- [ ] **Step A.3.1: Copy + create + register + delete**

Same recipe as Task A.2:
- Locate `module meteo_process_mod … end module meteo_process_mod` in `meteoday.f90`.
- Copy verbatim into new file `src/atmosphere/meteo_io.f90`.
- Add `'src/atmosphere/meteo_io.f90'` to `meson.build`.
- Delete the `module meteo_process_mod` block from `meteoday.f90`.

- [ ] **Step A.3.2: Build + test**

Run: `pixi run build-linux 2>&1 | tail -10 && pixi run test-pfunit 2>&1 | tail -5`
Expected: clean build + tests pass.

- [ ] **Step A.3.3: Commit**

```bash
git add src/atmosphere/meteo_io.f90 src/atmosphere/meteoday.f90 meson.build
git commit -m "$(cat <<'EOF'
refactor(atmosphere): extract meteo_process_mod into meteo_io.f90

Phase A.3 of GR-ATM-CLEAN. ReadMeteoDay + ResetMetFlx, body byte-identical.

Spec: docs/superpowers/specs/2026-05-22-atmosphere-cleanup-arc-design.md
EOF
)"
```

## Task A.4: Extract `module meteo_mod` → `meteo_orchestrator.f90`, delete `meteoday.f90`

- [ ] **Step A.4.1: Move final module + delete original file**

- `meteoday.f90` should now contain only `module meteo_mod … end module meteo_mod`. Verify.
- Copy `meteoday.f90` to `src/atmosphere/meteo_orchestrator.f90`.
- Delete `src/atmosphere/meteoday.f90`.

- [ ] **Step A.4.2: Update `meson.build`**

Remove `'src/atmosphere/meteoday.f90'` from the atmosphere sources list. Add `'src/atmosphere/meteo_orchestrator.f90'`.

- [ ] **Step A.4.3: Verification gate — full check-fast**

Run: `pixi run check-fast 2>&1 | tail -30`
Expected: pFUnit green + 4/4 regression cases match (hupselbrook, surfacewater, salinitystress, grassgrowth).

Byte-identical regression match expected — Phase A is a pure file split.

- [ ] **Step A.4.4: Commit + update README**

```bash
git add src/atmosphere/meteo_orchestrator.f90 src/atmosphere/meteoday.f90 meson.build
git commit -m "$(cat <<'EOF'
refactor(atmosphere): extract meteo_mod into meteo_orchestrator.f90; delete meteoday.f90

Phase A.4 of GR-ATM-CLEAN. Final extraction; meteoday.f90 deleted.
Module body byte-identical to original. check-fast green.

Spec: docs/superpowers/specs/2026-05-22-atmosphere-cleanup-arc-design.md
EOF
)"
```

Also update [src/atmosphere/README.md](../../../src/atmosphere/README.md) — the file/module table at the top should reflect the new layout. This can be a follow-up commit at the end of the phase.

---

# Phase B — Constants hoist

**Goal:** Populate `atmosphere_constants_mod` with shared physical constants; replace inline literals in `snow.f90`, `runoff.f90`, `et.f90`.

## Task B.1: Populate `atmosphere_constants_mod`

**File:** Modify `src/atmosphere/atmosphere_constants.f90` (currently a header-comment stub).

- [ ] **Step B.1.1: Replace `atmosphere_constants.f90` with the module declaration**

Write the file with:

```fortran
!> Shared physical constants for the atmosphere subsystem.
!! Single source of truth for magic numbers that previously appeared
!! inline in snow.f90, runoff.f90, et.f90. Names are uppercase + SI-unit-suffixed.
module atmosphere_constants_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   public

   ! ----- Thermodynamics -----
   !> Specific heat capacity of liquid water [J/kg/K]
   real(real64), parameter :: SPECIFIC_HEAT_WATER = 4180.0_real64

   !> Latent heat of melting for water/ice [J/kg]
   real(real64), parameter :: LATENT_HEAT_MELTING = 333580.0_real64

   !> Reference snow temperature for melt-energy calculation [deg C]
   real(real64), parameter :: SNOW_TEMPERATURE_C = 0.0_real64

   ! ----- Snow -----
   !> Maximum liquid-water storage in snow as a fraction of total water storage
   real(real64), parameter :: SNOW_LIQUID_WATER_FRACTION = 0.07_real64

   !> Soil-surface temperature above which fresh snow cannot accumulate [deg C]
   real(real64), parameter :: SOIL_SURFACE_FREEZE_THRESHOLD_C = 0.5_real64

   ! ----- Ponding -----
   !> Minimum ponding depth at which the soil surface is treated as wet [cm]
   real(real64), parameter :: POND_THRESHOLD_CM = 1.0e-10_real64

   ! ----- Curve Number / runoff -----
   !> Reference top-soil depth for the SCS-CN moisture correction [cm]
   real(real64), parameter :: DEPTH_10CM_CM = 10.0_real64

   !> Pressure head at field capacity used in CN ThetaRef calc [cm]
   real(real64), parameter :: H_FIELD_CAPACITY_CM = -100.0_real64

   !> Pressure head at wilting point used in CN ThetaRef calc [cm]
   real(real64), parameter :: H_WILTING_POINT_CM = -16000.0_real64

   !> Initial-abstraction ratio (Ia/S) in the SCS-CN runoff equation
   real(real64), parameter :: INITIAL_ABSTRACTION_RATIO = 0.2_real64

end module atmosphere_constants_mod
```

- [ ] **Step B.1.2: Ensure `meson.build` lists the file**

`grep -n atmosphere_constants meson.build` — verify it's already in the sources list (it should be — `atmosphere_constants.f90` was always a compiled file even when empty). If not, add it.

- [ ] **Step B.1.3: Build**

Run: `pixi run build-linux 2>&1 | tail -10`
Expected: clean build. The constants exist now but aren't yet used.

- [ ] **Step B.1.4: Commit**

```bash
git add src/atmosphere/atmosphere_constants.f90
git commit -m "$(cat <<'EOF'
feat(atmosphere): populate atmosphere_constants_mod with shared constants

Phase B.1 of GR-ATM-CLEAN. Adds named parameters for snow/pond/CN
physical constants previously inline in snow.f90, runoff.f90, et.f90.
Callers will be migrated in Task B.2.

Spec: docs/superpowers/specs/2026-05-22-atmosphere-cleanup-arc-design.md
EOF
)"
```

## Task B.2: Replace inline literals across callers

**Files:** Modify `src/atmosphere/snow.f90`, `src/atmosphere/runoff.f90`, `src/atmosphere/et.f90`.

- [ ] **Step B.2.1: Migrate `snow.f90`**

In `snow.f90`:

1. Add to the module's `use` block (top of `module snow_mod`):
   ```fortran
   use atmosphere_constants_mod, only: SPECIFIC_HEAT_WATER, LATENT_HEAT_MELTING, &
                                       SNOW_TEMPERATURE_C, SNOW_LIQUID_WATER_FRACTION, &
                                       SOIL_SURFACE_FREEZE_THRESHOLD_C
   ```

2. Inside `subroutine snow`, delete the local parameter declarations:
   ```fortran
   real(8), parameter :: cwat = 4180.0d0
   real(8), parameter :: lm = 333580d0
   real(8), parameter :: ts = 0.0d0
   ```

3. Replace identifier uses:
   - `cwat` → `SPECIFIC_HEAT_WATER`
   - `lm` → `LATENT_HEAT_MELTING`
   - `ts` → `SNOW_TEMPERATURE_C`

4. Find the inline literal `0.07*(at_slw + at_ssnow)` and replace with `SNOW_LIQUID_WATER_FRACTION * (at_slw + at_ssnow)`.

5. Find the inline `if (tsoil_surf .gt. 0.5d0 ...)` and replace with `if (tsoil_surf > SOIL_SURFACE_FREEZE_THRESHOLD_C ...)`.

- [ ] **Step B.2.2: Migrate `runoff.f90`**

In `runoff.f90`:

1. Add to `module runoff_mod` `use` block:
   ```fortran
   use atmosphere_constants_mod, only: DEPTH_10CM_CM, H_FIELD_CAPACITY_CM, &
                                       H_WILTING_POINT_CM, INITIAL_ABSTRACTION_RATIO
   ```

2. Delete the existing file-level `parameter` block:
   ```fortran
   real(8), parameter :: DEPTH_10CM = 10.0d0
   real(8), parameter :: H_FIELD_CAPACITY = -100.0d0
   real(8), parameter :: H_WILTING_POINT = -16000.0d0
   real(8), parameter :: IA_RATIO = 0.2d0
   ```

3. Replace identifier uses:
   - `DEPTH_10CM` → `DEPTH_10CM_CM`
   - `H_FIELD_CAPACITY` → `H_FIELD_CAPACITY_CM`
   - `H_WILTING_POINT` → `H_WILTING_POINT_CM`
   - `IA_RATIO` → `INITIAL_ABSTRACTION_RATIO`

4. Find the inline `Ia = 0.2d0*S` and replace with `Ia = INITIAL_ABSTRACTION_RATIO * S`.

- [ ] **Step B.2.3: Migrate `et.f90` (`reduceva`)**

In `et.f90`:

1. Add to `module et_mod` `use` block:
   ```fortran
   use atmosphere_constants_mod, only: POND_THRESHOLD_CM
   ```

2. Inside `subroutine reduceva`, delete the local `POND_THRESHOLD` parameter declaration.

3. Replace identifier uses of `POND_THRESHOLD` with `POND_THRESHOLD_CM`.

- [ ] **Step B.2.4: Verification gate — full check-fast**

Run: `pixi run check-fast 2>&1 | tail -30`
Expected: pFUnit green + 4/4 regression cases match **byte-identically** (named-constant substitution does not change arithmetic).

If any regression case shows drift: the substitution introduced a typo (e.g. wrong value in `atmosphere_constants_mod`). Diff `atmosphere_constants.f90` against the original inline literals (case-insensitive) and fix.

- [ ] **Step B.2.5: Commit**

```bash
git add src/atmosphere/snow.f90 src/atmosphere/runoff.f90 src/atmosphere/et.f90
git commit -m "$(cat <<'EOF'
refactor(atmosphere): consume atmosphere_constants_mod in snow/runoff/et

Phase B.2 of GR-ATM-CLEAN. Replaces inline literals + file-local
parameters with shared module constants. check-fast byte-identical.

Spec: docs/superpowers/specs/2026-05-22-atmosphere-cleanup-arc-design.md
EOF
)"
```

---

# Phase C — Task selectors → named init/step pairs

**Goal:** Replace `snow(task, …)`, `CNmethod(Itask, …)`, `reduceva(task, …)` with named init/step pairs. Update all 7 call sites.

## Task C.1: Split `snow` → `snow_init` + `snow_step`

**Files:** Modify `src/atmosphere/snow.f90`, `src/core/swap_mod.f90`.

- [ ] **Step C.1.1: Rewrite `snow.f90`**

In `module snow_mod`:

1. Change the `public` line:
   ```fortran
   public :: snow_init, snow_step
   ```
   (delete `public :: snow`)

2. Replace `subroutine snow(task, state)` with two routines. The recipe is mechanical: take the `case (1)` body verbatim for `snow_init`, the `case (2)` body verbatim for `snow_step`. Drop the `select case (task)` wrapper and the `case default … fatalerr_collected` branch.

```fortran
   subroutine snow_init(state)
      use Variables
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none
      type(swap_state_t), intent(inout) :: state

      ! === initialization ===================================================

      if (swinco .eq. 3) then
         state%atmosphere%snowinco = state%atmosphere%ssnow
      else
         state%atmosphere%ssnow = state%atmosphere%snowinco
      end if
   end subroutine snow_init

   subroutine snow_step(state)
      use Variables
      use, intrinsic :: iso_fortran_env, only: real64
      use atmosphere_constants_mod, only: SPECIFIC_HEAT_WATER, LATENT_HEAT_MELTING, &
                                          SNOW_TEMPERATURE_C, SNOW_LIQUID_WATER_FRACTION, &
                                          SOIL_SURFACE_FREEZE_THRESHOLD_C
      implicit none
      type(swap_state_t), intent(inout) :: state

      ! ... (verbatim body of former case (2), with the associate block
      !      and all dual-write aliases preserved) ...
   end subroutine snow_step
```

Subagent: copy the entire `case (2)` block (including the `associate(...)` aliases and `end associate`) verbatim into `snow_step`. Do not refactor the body — that's out of scope for this phase.

- [ ] **Step C.1.2: Update callers in `swap_mod.f90`**

In `src/core/swap_mod.f90`:

1. Line ~110: `use snow_mod, only: snow` → `use snow_mod, only: snow_init`
2. Line ~432: `if (flSnow) call Snow(1, state)` → `if (flSnow) call snow_init(state)`
3. Line ~477: `use snow_mod, only: snow` → `use snow_mod, only: snow_step`
4. Line ~550: `if (flSnow .and. tc_flDayStart) call Snow(2, state)` → `if (flSnow .and. tc_flDayStart) call snow_step(state)`

- [ ] **Step C.1.3: Verification gate**

Run: `pixi run check-fast 2>&1 | tail -30`
Expected: pFUnit green + 4/4 regression cases match.

- [ ] **Step C.1.4: Commit**

```bash
git add src/atmosphere/snow.f90 src/core/swap_mod.f90
git commit -m "$(cat <<'EOF'
refactor(snow): split snow(task) into snow_init + snow_step

Phase C.1 of GR-ATM-CLEAN. Removes the integer-task dispatcher in
snow_mod; callers in swap_mod.f90 now use the named init/step pair.

Spec: docs/superpowers/specs/2026-05-22-atmosphere-cleanup-arc-design.md
EOF
)"
```

## Task C.2: Split `CNmethod` → `cn_init` + `cn_step`

**Files:** Modify `src/atmosphere/runoff.f90`, `src/atmosphere/meteo_orchestrator.f90`, `src/core/swap_mod.f90`.

- [ ] **Step C.2.1: Rewrite `runoff.f90`**

In `module runoff_mod`:

1. Change `public :: CNmethod` to `public :: cn_init, cn_step`.

2. Replace `subroutine CNmethod(Itask, state)` with two routines extracted from `case (1)` and `case (2)`:

```fortran
   subroutine cn_init(state)
      use error_mod, only: fatalerr_collected
      use soilhydraulics_utils, only: watcon
      use atmosphere_constants_mod, only: DEPTH_10CM_CM, H_FIELD_CAPACITY_CM, &
                                          H_WILTING_POINT_CM
      implicit none
      type(swap_state_t), intent(inout) :: state

      integer :: i
      real(8) :: wc1, wc2

      associate( tc_t1900 => state%timecontrol%t1900 )
      ! ... verbatim body of former case (1) ...
      end associate
   end subroutine cn_init

   subroutine cn_step(state)
      use atmosphere_constants_mod, only: INITIAL_ABSTRACTION_RATIO
      implicit none
      type(swap_state_t), intent(inout) :: state

      integer :: i
      real(8) :: CN, S, Ia

      associate( tc_t1900 => state%timecontrol%t1900 )
      ! ... verbatim body of former case (2) ...
      end associate
   end subroutine cn_step
```

Subagent: extract the `case (1)` body verbatim into `cn_init`, the `case (2)` body verbatim into `cn_step`. Wrap each in its own `associate( tc_t1900 => state%timecontrol%t1900 )` block matching the original. Drop the outer `select case`, the `case default` fatalerr, and the `Itask` argument.

- [ ] **Step C.2.2: Update callers**

In `src/core/swap_mod.f90` line ~384:
- Old: `if (swuseCN == 1) call CNmethod(1, state)`
- New: `if (swuseCN == 1) call cn_init(state)`
- Update the `use runoff_mod, only: CNmethod` → `use runoff_mod, only: cn_init` at the top of that routine.

In `src/atmosphere/meteo_orchestrator.f90`:
- Find `use runoff_mod, only: CNmethod` — change to `use runoff_mod, only: cn_step`. (Verify there's no other use of `CNmethod` in this file before changing — `grep -n CNmethod src/atmosphere/meteo_orchestrator.f90` should show only the one call site at the former line 957.)
- Line ~957: `call CNmethod(2, state)` → `call cn_step(state)`

- [ ] **Step C.2.3: Verification gate**

Run: `pixi run check-fast 2>&1 | tail -30`
Expected: pFUnit green + 4/4 regression cases match.

- [ ] **Step C.2.4: Commit**

```bash
git add src/atmosphere/runoff.f90 src/atmosphere/meteo_orchestrator.f90 src/core/swap_mod.f90
git commit -m "$(cat <<'EOF'
refactor(runoff): split CNmethod(Itask) into cn_init + cn_step

Phase C.2 of GR-ATM-CLEAN. Removes the integer-task dispatcher in
runoff_mod; two callers (swap_mod, meteo_orchestrator) updated.

Spec: docs/superpowers/specs/2026-05-22-atmosphere-cleanup-arc-design.md
EOF
)"
```

## Task C.3: Split `reduceva` → `reduceva_daily` + `reduceva_dt`

**Files:** Modify `src/atmosphere/et.f90`, `src/atmosphere/meteo_orchestrator.f90`, `src/atmosphere/meteodt.f90`.

- [ ] **Step C.3.1: Rewrite `reduceva` in `et.f90`**

In `module et_mod`:

1. Change `public :: PenMon, PenMon_calc, reduceva` to `public :: PenMon, PenMon_calc, reduceva_daily, reduceva_dt`.

2. Replace `subroutine reduceva(task, nrai, state)` with two routines. They share a non-trivial body (the ponding short-circuit, the Black vs Boesten-Stroosnijder dispatch). Extract the shared body into a private helper to avoid duplication:

```fortran
   ! Public — daily-mode reduction
   subroutine reduceva_daily(nrai, state)
      implicit none
      real(8),            intent(in)    :: nrai
      type(swap_state_t), intent(inout) :: state
      call reduceva_apply(1, nrai, state, daily=.true.)
   end subroutine reduceva_daily

   ! Public — sub-daily reduction
   subroutine reduceva_dt(nrai, state)
      implicit none
      real(8),            intent(in)    :: nrai
      type(swap_state_t), intent(inout) :: state
      call reduceva_apply(2, nrai, state, daily=.false.)
   end subroutine reduceva_dt

   ! Private — shared body
   subroutine reduceva_apply(task_flag, nrai, state, daily)
      ! body equivalent to the former reduceva, but with `task` replaced
      ! by `task_flag` and the `if (task == 1) timestep = 1.0d0; else
      ! timestep = tc_dt` branch replaced by `if (daily) ... else ...`.
      ! Pass `task_flag` into `black_reduction` unchanged (it still uses
      ! the integer 1/2 internally to distinguish daily-vs-dt formulas).
      ...
   end subroutine reduceva_apply
```

Subagent: the cleanest implementation just keeps `reduceva_apply` private and re-uses the existing `select case (swredu)` dispatch into `black_reduction` / `boesten_stroosnijder_reduction`. `black_reduction` already takes a `task` argument — that stays, since it gates the daily vs sub-daily formula.

- [ ] **Step C.3.2: Update callers**

In `src/atmosphere/meteo_orchestrator.f90`:
- Line ~965: `call reduceva (1, state%atmosphere%nraida, state)` → `call reduceva_daily(state%atmosphere%nraida, state)`
- Find the `use et_mod, only: ...` import that currently includes `reduceva` (likely module-level around line 519); change `reduceva` → `reduceva_daily`.

In `src/atmosphere/meteodt.f90`:
- Line ~429: `call reduceva(2, state%atmosphere%nraida, state)` → `call reduceva_dt(state%atmosphere%nraida, state)`
- Line ~541: same substitution.
- Update the two `use et_mod, only: reduceva` imports (one per subroutine that calls reduceva) to `use et_mod, only: reduceva_dt`.

- [ ] **Step C.3.3: Verification gate**

Run: `pixi run check-fast 2>&1 | tail -30`
Expected: pFUnit green + 4/4 regression cases match.

- [ ] **Step C.3.4: Commit**

```bash
git add src/atmosphere/et.f90 src/atmosphere/meteo_orchestrator.f90 src/atmosphere/meteodt.f90
git commit -m "$(cat <<'EOF'
refactor(et): split reduceva(task) into reduceva_daily + reduceva_dt

Phase C.3 of GR-ATM-CLEAN. Removes the integer-task dispatcher in
et_mod; three callers (meteo_orchestrator x1, meteodt x2) updated.
Shared body in private reduceva_apply helper.

Spec: docs/superpowers/specs/2026-05-22-atmosphere-cleanup-arc-design.md
EOF
)"
```

---

# Phase D — Retire `MeteoVars`

**Goal:** Delete `meteo_vars.f90`. Categorise each member by lifetime; promote scratch to locals, persistent state to `state%atmosphere`.

**Critical:** this phase adds three fields to `state%atmosphere`. After Task D.1 lands the state-schema change, the next build **requires a clean rebuild** — `rm -rf builddir` — per the project's `feedback_state_schema_clean_rebuild` rule. Incremental Meson builds don't propagate `.mod` deps across the `swap_modern` / `swap_legacy` boundary on state-schema edits.

## Task D.1: Add new fields to `state%atmosphere`

**Files:** Modify `src/state/atmosphere_state.f90`.

- [ ] **Step D.1.1: Locate `atmosphere_state_t`**

Read `src/state/atmosphere_state.f90` and find the `type :: atmosphere_state_t` block. Identify a sensible insertion point — likely at the end of the type before `end type`.

- [ ] **Step D.1.2: Add three fields**

Insert (preserving the file's existing formatting style):

```fortran
   ! [GR-ATM-CLEAN Phase D] migrated from module MeteoVars (meteo_vars.f90).
   !> Remaining interception storage at start of timestep (cm).
   !> Persists across iterations of the sub-daily dayparts loop and across days.
   real(real64) :: restint = 0.0_real64

   !> Sub-daily precipitation per record (cm). Populated by ReadMeteoDay,
   !> consumed by ProcessMeteoDay's sub-daily branch. Allocation size = nmetdetail (<= 96).
   real(real64) :: arain_subdaily(96) = 0.0_real64

   !> Sub-daily wind speed per record (m/s). Same lifetime as arain_subdaily.
   real(real64) :: awind_subdaily(96) = 0.0_real64
```

- [ ] **Step D.1.3: Clean rebuild + verification gate**

```bash
rm -rf builddir
pixi run check-fast 2>&1 | tail -30
```

Expected: pFUnit green + 4/4 regression cases match. The new fields exist but no one reads them yet, so behaviour is unchanged.

- [ ] **Step D.1.4: Commit**

```bash
git add src/state/atmosphere_state.f90
git commit -m "$(cat <<'EOF'
feat(state): add restint + arain/awind_subdaily to atmosphere_state_t

Phase D.1 of GR-ATM-CLEAN. Pre-stages the three cross-call MeteoVars
members. Callers will be migrated in Task D.2.

Spec: docs/superpowers/specs/2026-05-22-atmosphere-cleanup-arc-design.md
EOF
)"
```

## Task D.2: Retire `meteo_vars.f90`

**Files:** Modify `src/atmosphere/meteo_io.f90`, `src/atmosphere/meteo_orchestrator.f90`; delete `src/atmosphere/meteo_vars.f90`; modify `meson.build`.

- [ ] **Step D.2.1: Migrate `meteo_io.f90` (`ReadMeteoDay`, `ResetMetFlx`)**

In `meteo_io.f90`:

1. Remove `use MeteoVars`.
2. Identify every reference inside the module body to a MeteoVars symbol. The members are:
   - **Pure loop locals** (`count, first, i, irecord, last, ndayparts`): if used in this file, declare as routine locals at the top of the routine.
   - **Within-call scratch** (`interc, aintc, eintc, dttp, gctp, etr, hum, win, svp, wfrac, netrainflux, rainflux, sumtav, Edirect, Tdirect, Tdirectwet, Edirectpond`): declare as routine locals at the top of the routine.
   - **Cross-call persistent**: replace `restint` → `state%atmosphere%restint`, `arain(i)` → `state%atmosphere%arain_subdaily(i)`, `awind(i)` → `state%atmosphere%awind_subdaily(i)`.

3. The `ResetMetFlx` subroutine is small (~20 lines, mostly the cohort `reset` dispatches) — quickly grep for any MeteoVars use; likely none.

4. `ReadMeteoDay` uses `arain`, `awind` (lines around `arain(i) = detrain(...) * 0.1d0`, `awind(i) = detwind(...)`) — switch to `state%atmosphere%arain_subdaily(i)`, `state%atmosphere%awind_subdaily(i)`.

Subagent: be thorough. Run `grep -nE '\b(count|first|i|irecord|last|ndayparts|arain|awind|restint|interc|Edirectpond|aintc|dttp|eintc|etr|gctp|hum|svp|wfrac|win|netrainflux|rainflux|sumtav|Edirect|Tdirect|Tdirectwet)\b' src/atmosphere/meteo_io.f90` to find every use. Migrate them all. Some of those names are also used as local variables in unrelated contexts (`i` especially is a generic loop counter) — declare them as locals in the routine and confirm no shadowing surprises.

- [ ] **Step D.2.2: Migrate `meteo_orchestrator.f90` (`ProcessMeteoDay`)**

In `meteo_orchestrator.f90`:

1. Remove `use MeteoVars`.
2. Run the same grep to enumerate uses.
3. Inside `ProcessMeteoDay`, at the top of the routine (alongside `rcs`, `pmi`, `pmo` declared by the PenMon refactor), add local declarations for every MeteoVars member used in the body. Subagent: prefer one block of `real(8)` declarations for the scratch scalars, one for `integer` locals.
4. Replace `restint`/`arain(...)`/`awind(...)` with `state%atmosphere%restint`/`state%atmosphere%arain_subdaily(...)`/`state%atmosphere%awind_subdaily(...)`.
5. The `do 1000 irecord = 1, ndayparts` loop: `irecord` becomes a local integer (already declared above); `ndayparts` becomes a local integer set near the loop entry.

- [ ] **Step D.2.3: Delete `meteo_vars.f90` and update `meson.build`**

```bash
rm src/atmosphere/meteo_vars.f90
```

Edit `meson.build` to remove `'src/atmosphere/meteo_vars.f90'` from the atmosphere sources list.

- [ ] **Step D.2.4: Verification gate**

```bash
rm -rf builddir
pixi run check-fast 2>&1 | tail -30
```

Expected: pFUnit green + 4/4 regression cases match.

**If regression drift:** the most likely cause is a missed `restint` / `arain` / `awind` site that's still trying to read from the (now-deleted) MeteoVars module — but the compile would have failed in that case, so drift implies an arithmetic transcription bug. Most likely a `state%atmosphere%arain_subdaily` index off-by-one, or a write-then-read ordering issue in the migrated locals. Diff against `HEAD~1` and audit.

- [ ] **Step D.2.5: Commit**

```bash
git add src/atmosphere/meteo_io.f90 src/atmosphere/meteo_orchestrator.f90 src/atmosphere/meteo_vars.f90 meson.build
git commit -m "$(cat <<'EOF'
refactor(atmosphere): retire MeteoVars; scratch scalars to locals, arrays to state

Phase D.2 of GR-ATM-CLEAN. Deletes meteo_vars.f90 entirely. Loop
counters and within-call scratch become routine locals; restint and
the arain/awind sub-daily arrays migrate to state%atmosphere.

Spec: docs/superpowers/specs/2026-05-22-atmosphere-cleanup-arc-design.md
EOF
)"
```

---

# Phase E — Daily/sub-daily orchestrator split

**Goal:** Split `ProcessMeteoDay` into `process_meteo_day_daily` and `process_meteo_day_subdaily`, sharing private helpers. `ProcessMeteoDay` becomes a 5-line dispatcher (preserves the public API).

This is the largest phase. Extract one helper at a time, verifying check-fast after each, so any regression is bisectable to a single helper extraction.

## Task E.1: Extract `apply_interception_step` helper

**File:** Modify `src/atmosphere/meteo_orchestrator.f90`.

- [ ] **Step E.1.1: Identify the interception block**

In `ProcessMeteoDay`, the interception block (Section 3 of the labelled sections) starts around the comment `! === Section 3: Interception calculations ===` and ends before the `! === LOOP over dayparts ===` comment. It contains the `VonHHBraden` / `Gash` dispatch + `DivIntercep`. The Rutter case (Section 5) calls `ruttervw + DivIntercep` *inside* the dayparts loop — it's a separate path.

Extract Section 3 into:

```fortran
   ! Private helper: handles the Section-3 interception calculation
   ! (VonHHBraden / Gash dispatch + DivIntercep). Returns aintc and updates
   ! state%atmosphere%nraida via DivIntercep.
   subroutine apply_interception_step(state, config, aintc)
      use interception_mod, only: VonHHBraden, Gash, DivIntercep
      use swap_state_mod, only: swap_state_t
      use swap_config_mod, only: swap_config_t
      implicit none
      type(swap_state_t),  intent(inout) :: state
      type(swap_config_t), intent(in)    :: config
      real(8), intent(out) :: aintc

      ! ... copy lines from the former Section 3 of ProcessMeteoDay ...
   end subroutine apply_interception_step
```

In `ProcessMeteoDay`, replace the Section-3 inline block with `call apply_interception_step(state, config, aintc)`.

- [ ] **Step E.1.2: Verification gate**

Run: `pixi run check-fast 2>&1 | tail -30`
Expected: pFUnit green + 4/4 regression cases match.

- [ ] **Step E.1.3: Commit**

```bash
git add src/atmosphere/meteo_orchestrator.f90
git commit -m "$(cat <<'EOF'
refactor(meteo_orchestrator): extract apply_interception_step helper

Phase E.1 of GR-ATM-CLEAN. Lifts Section 3 (interception calc +
DivIntercep) out of ProcessMeteoDay into a private helper. Same body,
called identically. Prepares the orchestrator split in Tasks E.4.

Spec: docs/superpowers/specs/2026-05-22-atmosphere-cleanup-arc-design.md
EOF
)"
```

## Task E.2: Extract `compute_reference_et` helper

**File:** Modify `src/atmosphere/meteo_orchestrator.f90`.

- [ ] **Step E.2.1: Identify the ET-compute block**

In `ProcessMeteoDay`, Section 4 has two top-level branches:
- `if (swmetdetail.eq.0 .and. swetr.eq.1)`: use supplied `etr` directly (no PenMon).
- `elseif (swmetdetail.eq.1 .or. swetr.eq.0)`: compute via PenMon — includes the sub-daily per-record refresh, the PenMon pack/call/unpack, and the post-call `if .not. flCropEmergence then ... else ...` adjustment.

Extract both branches into:

```fortran
   ! Private helper: computes es0/ew0/et0 (and PMdirect outputs) for one
   ! record. Handles both the etr-direct and PenMon-based paths.
   ! For sub-daily mode the caller refreshes rad/Tav/hum/win before calling.
   subroutine compute_reference_et(state, config, irecord, &
                                   rad, hum, win, tmn, tmx, &
                                   rcs, angstroma, angstromb, daylp, difpp, dsinbe, atmtr, rsoil, &
                                   logf, swscre, &
                                   pmo)
      use et_mod, only: PenMon, pm_inputs_t, pm_outputs_t
      use swap_state_mod, only: swap_state_t
      use swap_config_mod, only: swap_config_t
      implicit none
      type(swap_state_t),  intent(inout) :: state
      type(swap_config_t), intent(in)    :: config
      integer, intent(in) :: irecord, logf, swscre
      real(8), intent(in) :: rad, hum, win, tmn, tmx
      real(8), intent(in) :: rcs, angstroma, angstromb, daylp
      real(8), intent(in) :: difpp, dsinbe, atmtr, rsoil
      type(pm_outputs_t), intent(out) :: pmo

      type(pm_inputs_t) :: pmi

      ! ... pack pmi from arguments + state + config ...
      ! ... call PenMon(pmi, pmo, logf, swscre) ...
      ! ... post-call swcf/swcfbs adjustment to es0/et0/ew0 ...
      ! ... write state%crop%es0/et0/ew0 ...
   end subroutine compute_reference_et
```

Subagent: this helper subsumes the existing PenMon pack/call/unpack added by `c29408b`. Migrate that pack/call/unpack into the helper, plus the surrounding `if (config%meteo%swmetdetail.eq.0 .and. config%meteo%swetr.eq.1) ... elseif ...` selector and the post-call crop-factor adjustments.

Argument-list discipline: only pass what the helper needs. Everything else (state, config) reads from the typed-state passed in. The list above is intentionally short — many of the existing locals (`angstroma`, `daylp`, `rsc`, etc.) are read inside via `state%crop%...` or `state%cfg%meteo%...` already. Re-examine and prune.

In `ProcessMeteoDay`, replace Section 4 inline block with `call compute_reference_et(state, config, irecord, ..., pmo)`. Then unpack `pmo` to whatever locals still need it (mostly `Edirect`, `Tdirect`, `Tdirectwet`, `Edirectpond`).

- [ ] **Step E.2.2: Verification gate**

Run: `pixi run check-fast 2>&1 | tail -30`
Expected: pFUnit green + 4/4 regression cases match.

- [ ] **Step E.2.3: Commit**

```bash
git add src/atmosphere/meteo_orchestrator.f90
git commit -m "$(cat <<'EOF'
refactor(meteo_orchestrator): extract compute_reference_et helper

Phase E.2 of GR-ATM-CLEAN. Lifts Section 4 (etr-direct or PenMon
+ crop-factor adjustments) out of ProcessMeteoDay into a private
helper. Subsumes the PenMon pack/call/unpack from c29408b.

Spec: docs/superpowers/specs/2026-05-22-atmosphere-cleanup-arc-design.md
EOF
)"
```

## Task E.3: Extract `compute_wet_fraction` + `partition_peva_ptra` helpers

**File:** Modify `src/atmosphere/meteo_orchestrator.f90`.

- [ ] **Step E.3.1: Extract `compute_wet_fraction`**

In `ProcessMeteoDay`, Section 6 computes `wfrac` based on `aintc`, `eintc`, `ew0`, `swinter`, `swdivide`. The two top-level branches (`if (swmetdetail.eq.0) ... elseif (swmetdetail.eq.1)`) live inside this section.

Extract into:

```fortran
   subroutine compute_wet_fraction(state, config, aintc, eintc, &
                                   interc, restint, arain_irecord, irecord, wfrac)
      use swap_state_mod, only: swap_state_t
      use swap_config_mod, only: swap_config_t
      use swap_constants, only: nihil
      implicit none
      type(swap_state_t),  intent(in)    :: state
      type(swap_config_t), intent(in)    :: config
      real(8), intent(in)    :: aintc, eintc
      real(8), intent(inout) :: interc, restint
      real(8), intent(in)    :: arain_irecord
      integer, intent(in)    :: irecord
      real(8), intent(out)   :: wfrac

      ! ... copy lines from the former Section 6, both branches ...
   end subroutine compute_wet_fraction
```

Replace Section 6 inline with `call compute_wet_fraction(state, config, aintc, eintc, interc, state%atmosphere%restint, state%atmosphere%arain_subdaily(irecord), irecord, wfrac)`.

- [ ] **Step E.3.2: Extract `partition_peva_ptra`**

Section 7 computes `peva` and `ptra` with cover correction, ponding correction, PMdirect override, and CO2 correction. Extract:

```fortran
   subroutine partition_peva_ptra(state, config, wfrac, &
                                  Edirect, Tdirect, Edirectpond)
      use swap_state_mod, only: swap_state_t
      use swap_config_mod, only: swap_config_t
      use swap_constants, only: nihil, small
      use variables, only: cfevappond, flco2, croptype, flCropHarvest, gc
      implicit none
      type(swap_state_t),  intent(inout) :: state
      type(swap_config_t), intent(in)    :: config
      real(8), intent(in) :: wfrac, Edirect, Tdirect, Edirectpond

      ! ... copy lines from the former Section 7 ...
      ! Writes state%atmosphere%peva and state%atmosphere%ptra.
   end subroutine partition_peva_ptra
```

Replace Section 7 inline with `call partition_peva_ptra(state, config, wfrac, Edirect, Tdirect, Edirectpond)`.

- [ ] **Step E.3.3: Verification gate**

Run: `pixi run check-fast 2>&1 | tail -30`
Expected: pFUnit green + 4/4 regression cases match.

- [ ] **Step E.3.4: Commit**

```bash
git add src/atmosphere/meteo_orchestrator.f90
git commit -m "$(cat <<'EOF'
refactor(meteo_orchestrator): extract wfrac + peva/ptra partition helpers

Phase E.3 of GR-ATM-CLEAN. Lifts Sections 6 and 7 of ProcessMeteoDay
into compute_wet_fraction and partition_peva_ptra. Same body, called
identically.

Spec: docs/superpowers/specs/2026-05-22-atmosphere-cleanup-arc-design.md
EOF
)"
```

## Task E.4: Split `ProcessMeteoDay` into daily + sub-daily orchestrators

**File:** Modify `src/atmosphere/meteo_orchestrator.f90`.

- [ ] **Step E.4.1: Create the two new orchestrator subroutines**

After the helpers (which are now extracted), add two new private subroutines:

```fortran
   subroutine process_meteo_day_daily(state, config)
      use swap_state_mod, only: swap_state_t
      use swap_config_mod, only: swap_config_t
      implicit none
      type(swap_state_t),  intent(inout) :: state
      type(swap_config_t), intent(in)    :: config

      real(8) :: aintc, eintc, interc, wfrac
      real(8) :: Edirect, Tdirect, Tdirectwet, Edirectpond
      type(pm_outputs_t) :: pmo
      ! ... whatever other locals are needed (rad, hum, win, tmn, tmx — derived
      !     from state%atmosphere) ...

      ! Section 3: interception (one call)
      call apply_interception_step(state, config, aintc)

      ! Section 4: reference ET (one record)
      call compute_reference_et(state, config, irecord=1, &
                                rad=...(today), hum=..., win=..., tmn=..., tmx=..., &
                                ..., pmo=pmo)
      Edirect     = pmo%Edirect
      Tdirect     = pmo%Tdirect
      Tdirectwet  = pmo%Tdirectwet
      Edirectpond = pmo%Edirectpond

      ! Section 6 + 7
      call compute_wet_fraction(state, config, aintc, eintc, interc, &
                                state%atmosphere%restint, 0.0d0, 1, wfrac)
      call partition_peva_ptra(state, config, wfrac, Edirect, Tdirect, Edirectpond)

      ! Section 9: daily-only finalization
      ! - finterception, rainflux, graidt, nraidt, aintcdt for swrain.eq.0
      ! - CNmethod (cn_step) if swuseCN==1
      ! - reduceva_daily for actual soil evaporation
      ! - state%atmosphere%pevaday/ptraday assignment for ETSine carryover
      ! - state%atmosphere%atmdem computation
      ! ... (body verbatim from former Section 9) ...
   end subroutine process_meteo_day_daily

   subroutine process_meteo_day_subdaily(state, config)
      use swap_state_mod, only: swap_state_t
      use swap_config_mod, only: swap_config_t
      implicit none
      type(swap_state_t),  intent(inout) :: state
      type(swap_config_t), intent(in)    :: config

      integer :: irecord
      real(8) :: aintc, eintc, interc, wfrac
      real(8) :: Edirect, Tdirect, Tdirectwet, Edirectpond
      type(pm_outputs_t) :: pmo
      ! ... per-record locals refreshed each iteration ...

      ! Section 3: interception (once per day)
      call apply_interception_step(state, config, aintc)

      ! Loop over sub-daily records
      do irecord = 1, config%meteo%nmetdetail
         ! Refresh per-record meteo (rad/Tav/hum/win)
         ! ... from state%atmosphere%arad(irecord), atav, ahum, awind_subdaily ...

         call compute_reference_et(state, config, irecord, &
                                   rad=..., hum=..., win=..., ..., pmo=pmo)
         Edirect     = pmo%Edirect
         Tdirect     = pmo%Tdirect
         Tdirectwet  = pmo%Tdirectwet
         Edirectpond = pmo%Edirectpond

         ! Section 5: Rutter (sub-daily-specific, swinter==3)
         if (state%crop%common%swinter == 3) then
            ! ... siccapact from afgen, gctp calc, ruttervw call, DivIntercep ...
         end if

         ! Section 6 + 7
         call compute_wet_fraction(state, config, aintc, eintc, interc, &
                                   state%atmosphere%restint, &
                                   state%atmosphere%arain_subdaily(irecord), &
                                   irecord, wfrac)
         call partition_peva_ptra(state, config, wfrac, Edirect, Tdirect, Edirectpond)

         ! Section 8: per-record aggregation
         ! - state%atmosphere%tpot(irecord), epot(irecord), grain(irecord), nrain(irecord)
         ! ... (body verbatim from former Section 8) ...
      end do

      ! Section 10: daily totals from sub-daily records
      ! - Tav from sumtav, tmn/tmx from min/max
      ! - svp / rh
      ! - tavd (6h-18h average)
      ! - daily rad sum, atmdem sum
      ! - flux at start of day for ETSine seed
      ! ... (body verbatim from former Section 10) ...
   end subroutine process_meteo_day_subdaily
```

**Section-to-orchestrator mapping** (use as the cut sheet):

| Former section in `ProcessMeteoDay` | Identified by comment | Belongs to |
|------|-----------------------|------------|
| Section 3 (interception) | `! === Section 3: Interception ===` | both (helper) |
| Section 4 (compute ET) | `! === Section 4: Calculate evapotranspiration ===` | both (helper) |
| Section 5 (Rutter, swinter==3) | `! === Section 5: Interception option NHI ===` | sub-daily only inline; daily branch of Section 5 disappears because the `swmetdetail.eq.0` siccapact branch is no longer needed once the orchestrators are split (siccapact is set in cropgrowth for daily mode) |
| Section 6 (wet fraction) | `! === Section 6: Fraction of the day the crop is wet ===` | both (helper) |
| Section 7 (peva/ptra) | `! === Section 7: Potential soil evaporation & transpiration ===` | both (helper) |
| Section 8 (sub-daily aggregation) | `! === Section 8: Results for detailed weather records ===` | sub-daily only |
| Section 9 (daily finalisation) | `! === Section 9: Actual daily rain/snow fluxes for Daily Meteo ===` | daily only |
| Section 10 (sub-daily finalisation) | `! === Section 10: Set daily weather values for Detailed Meteo ===` | sub-daily only |

Subagent: this is the most delicate task in the arc. Recommended approach:

1. Copy the entire current `ProcessMeteoDay` body to a scratch location (a temp file, NOT the source).
2. Identify which sections / blocks belong to the daily path vs the sub-daily path. Use the existing `if (swmetdetail.eq.0)` and `if (swmetdetail.eq.1)` branches as the guide.
3. Build `process_meteo_day_daily` by concatenating: helper calls + daily-only blocks (Section 9).
4. Build `process_meteo_day_subdaily` by concatenating: helper calls (in loop) + Section 5 (Rutter) + Section 8 (per-record agg) + Section 10 (daily totals).
5. Both routines need the meteo refresh locals (`rad`, `hum`, `win`, `tmn`, `tmx`) — in daily mode these are set once from `state%atmosphere%X`; in sub-daily mode they're refreshed per iteration.

- [ ] **Step E.4.2: Reduce `ProcessMeteoDay` to a dispatcher**

```fortran
   subroutine ProcessMeteoDay(state, config)
      use swap_state_mod, only: swap_state_t
      use swap_config_mod, only: swap_config_t
      implicit none
      type(swap_state_t),  intent(inout) :: state
      type(swap_config_t), intent(in)    :: config

      if (config%meteo%swmetdetail == 0) then
         call process_meteo_day_daily(state, config)
      else
         call process_meteo_day_subdaily(state, config)
      end if
   end subroutine ProcessMeteoDay
```

Make the two new orchestrators `private` (only `ProcessMeteoDay` stays public — caller surface unchanged).

- [ ] **Step E.4.3: Final verification gate**

Run: `pixi run check-fast 2>&1 | tail -30`
Expected: pFUnit green + 4/4 regression cases match.

**If regression drift:** the most likely culprits are:
- Wrong branch picked up by the wrong orchestrator (e.g. Section 9 logic ending up in sub-daily, Section 10 ending up in daily).
- A control-flow guard dropped during the split (e.g. `if (swrain.eq.0)` guarding the rainflux/graidt assignment in Section 9 — must be preserved inside `process_meteo_day_daily`).
- Order-of-operations difference (e.g. extracting helpers before the meteo-refresh writes, causing PenMon to see stale `rad`).

Diff against the pre-split commit and walk the diff section-by-section. Do not "fix" regression values; the split is supposed to be byte-identical.

- [ ] **Step E.4.4: Commit**

```bash
git add src/atmosphere/meteo_orchestrator.f90
git commit -m "$(cat <<'EOF'
refactor(meteo_orchestrator): split ProcessMeteoDay into daily + subdaily

Phase E.4 of GR-ATM-CLEAN. ProcessMeteoDay becomes a 5-line dispatcher.
process_meteo_day_daily and process_meteo_day_subdaily are the two
work-doing routines, sharing helpers from Tasks E.1-E.3. Public API
unchanged.

Spec: docs/superpowers/specs/2026-05-22-atmosphere-cleanup-arc-design.md
EOF
)"
```

---

# Phase F — Intrinsic sweep

**Goal:** Replace F77-era intrinsics with their F90 generic counterparts across the atmosphere subsystem. Single commit.

## Task F.1: Sweep across `src/atmosphere/*.f90`

**Files:** All `.f90` files under `src/atmosphere/`.

- [ ] **Step F.1.1: Inventory current uses**

Run, for visibility:

```bash
grep -nE '\bdexp\b|\bdlog\b|\bdmin1\b|\bdble\b' src/atmosphere/*.f90
```

Expected output: a list of lines across multiple files (mostly inside `et.f90`, `interception.f90`, `meteodt.f90`).

- [ ] **Step F.1.2: Apply mechanical substitutions**

For each file in `src/atmosphere/*.f90`:

- `dexp(` → `exp(`
- `dlog(` → `log(`
- `dmin1(` → `min(`  (n-ary, but in this codebase only called with 2 args)
- `dble(` → `real(` … but watch for the second arg: `dble(x)` becomes `real(x, real64)`. This requires `real64` in scope.

Subagent: use `sed -i 's/\bdexp(/exp(/g; s/\bdlog(/log(/g; s/\bdmin1(/min(/g' src/atmosphere/*.f90` for the first three (they're plain renames). For `dble`, do it more carefully — `sed` can't add the `real64` kind argument. Use `grep -n 'dble(' src/atmosphere/*.f90` first, then do per-occurrence edits.

For each file that gets a `dble → real(.., real64)` change, add to the module's `use` block (if not already present):

```fortran
use, intrinsic :: iso_fortran_env, only: real64
```

- [ ] **Step F.1.3: Verification gate**

Run: `pixi run check-fast 2>&1 | tail -30`
Expected: pFUnit green + 4/4 regression cases match **byte-identically** (these are aliases the compiler folds identically).

**If regression drift:** the substitution was not byte-equivalent — most likely a `dble(integer)` call where the implicit conversion was already to `real(8)` and the new `real(int_var, real64)` should also produce identical bits. If drift persists, revert and investigate per-file.

- [ ] **Step F.1.4: Commit**

```bash
git add src/atmosphere
git commit -m "$(cat <<'EOF'
refactor(atmosphere): sweep dexp/dlog/dmin1/dble to F90 generics

Phase F.1 of GR-ATM-CLEAN. Mechanical substitution; check-fast
byte-identical. Closes the GR-ATM-CLEAN arc.

Spec: docs/superpowers/specs/2026-05-22-atmosphere-cleanup-arc-design.md
EOF
)"
```

---

# Arc completion

- [ ] Run a final `pixi run check-fast` to confirm the arc's end state.
- [ ] Update [src/atmosphere/README.md](../../../src/atmosphere/README.md) to reflect the new file layout (the "Module layout" table at the top) and remove items #1, #2, #3, #6, #7 from the "Cross-cutting modernisation themes" section.
- [ ] Report to the user: 6 phases complete, ~16 commits, check-fast green at every phase boundary.

**Next-arc candidates** (not in scope here):
- Item #4: msw1eic / ruttervw quarantine — needs its own design discussion.
- AFGEN cache in `Gash` — defer until profile evidence.
- `pm_atmos_baseline_t` partial precomputation for PenMon — defer.
- Drop dual-write transitional blocks — coordinate with GR-ATM Phase C3.
- `real(8)` → `real(real64)` declaration sweep — separate phase.
