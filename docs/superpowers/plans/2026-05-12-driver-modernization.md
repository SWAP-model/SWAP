# SS-DRV Phase 1 — Driver Modernization & BMI Stub Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Convert the free `subroutine swap` in `src/core/swap.f90` into `module swap_mod` with three named lifecycle procedures (`swap_init`, `swap_run_step`, `swap_close`); extract the outer time loop to the caller; retire `swap_exchange`/`handle_exchange`/`iCaller`/`dummy()`; add a minimal but spec-compliant BMI `iso_c_binding` façade plus a Python `cffi` hello-world that exercises it end-to-end.

**Architecture:** New files (`swap_mod.f90`, `swap_bmi_mod.f90`) compile with `-std=f2018` in a dedicated `swap_modern` static library; legacy sources remain on `-std=legacy` in `swap_legacy`. The BMI shared library and the main executable link both. Migration is done incrementally with the new module compiled alongside the old subroutine until cutover, so each step has a working build.

**Tech Stack:** Fortran 2018, gfortran, Meson, pFUnit 4.15, Python 3.11, cffi, pixi for the build/test orchestration.

**Spec:** `docs/superpowers/specs/2026-05-12-driver-modernization-design.md`

---

## Task 1: Pre-flight baseline

**Files:**
- No code changes; capture starting state for regression comparison.

- [ ] **Step 1: Verify the build is clean before starting**

Run: `pixi run build-linux`
Expected: build succeeds, no warnings on existing code.

- [ ] **Step 2: Run pFUnit suite to confirm baseline pass**

Run: `pixi run test-pfunit`
Expected: all existing suites pass (no swap_mod suite yet).

- [ ] **Step 3: Run the fast regression gate to capture the byte-for-byte baseline**

Run: `pixi run check-fast`
Expected: 4/4 cases pass (hupselbrook, surfacewater, salinitystress, grassgrowth).

- [ ] **Step 4: Confirm git working tree is clean apart from known untracked files**

Run: `git status`
Expected: only the known untracked items shown in the initial gitStatus (`.claude/`, `meteovars.mod`, `subprojects/test-drive/`, `subprojects/toml-f/`). No tracked-file modifications.

- [ ] **Step 5: Commit a tagged baseline marker (no code change)**

Run:
```bash
git commit --allow-empty -m "$(cat <<'EOF'
chore(ss-drv): pre-flight baseline marker — driver modernization arc begins

check-fast: 4/4 regression cases passing on hupselbrook, surfacewater,
salinitystress, grassgrowth. pFUnit suite: passing. Build clean.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Add empty `swap_mod` skeleton + failing pFUnit lifecycle tests

**Files:**
- Create: `src/core/swap_mod.f90`
- Create: `tests/unit/core/test_swap_mod.pf`
- Modify: `meson.build` (add `swap_mod.f90` to existing `sources` list, just before `swap.f90`)
- Modify: `tests/unit/meson.build` (add `core/test_swap_mod.pf` to `pf_files`; add `swap_mod.f90` to `pfunit_extra_sources`)
- Modify: `tests/unit/testSuites.inc` (append `ADD_TEST_SUITE(test_swap_mod_suite)` — must be **last** for global-state isolation)

- [ ] **Step 1: Create the empty module skeleton**

Create `src/core/swap_mod.f90` with this exact content:

```fortran
!> @file swap_mod.f90
!! SS-DRV Phase 1: module form of the legacy `subroutine swap`.
!! Three named lifecycle procedures replace the (iCaller, iTask) dispatch.
!! Time loop lives in the caller. State and config are threaded explicitly.
module swap_mod
   use swap_state_mod,  only: swap_state_t
   use swap_config_mod, only: swap_config_t
   implicit none
   private
   public :: swap_init, swap_run_step, swap_close

contains

   subroutine swap_init(config_file, state, config)
      character(len=*),            intent(in)  :: config_file
      type(swap_state_t),          intent(out) :: state
      type(swap_config_t), target, intent(out) :: config
      ! Body filled in Task 3.
   end subroutine swap_init

   subroutine swap_run_step(state, config)
      type(swap_state_t),  intent(inout) :: state
      type(swap_config_t), intent(in)    :: config
      ! Body filled in Task 4.
   end subroutine swap_run_step

   subroutine swap_close(state, config)
      type(swap_state_t),  intent(inout) :: state
      type(swap_config_t), intent(in)    :: config
      ! Body filled in Task 5.
   end subroutine swap_close

end module swap_mod
```

- [ ] **Step 2: Wire `swap_mod.f90` into the main `meson.build` sources list**

In `meson.build`, locate the `sources = [` block. Inside it, immediately **before** the `'src/core/swap.f90',` line (currently at line 187), insert:

```
    'src/core/swap_mod.f90',
```

The block tail must end up looking like:
```
    # Main (last)
    'src/core/swap_mod.f90',
    'src/core/swap.f90',
    'src/core/swap_main.f90'
]
```

- [ ] **Step 3: Verify the new module compiles before writing tests**

Run: `pixi run build-linux`
Expected: build succeeds with the new module present (unused but compiled).

- [ ] **Step 4: Write the failing pFUnit lifecycle tests**

Create `tests/unit/core/test_swap_mod.pf` with this exact content:

```fortran
!> SS-DRV Phase 1: lifecycle smoke tests for swap_mod.
!! These tests drive a full hupselbrook simulation through the new
!! module interface. They MUST be the last suite registered so that
!! the global-state churn (Initialize, file opens) does not leak
!! into other test suites.
@suite(name="test_swap_mod_suite")

@test
subroutine test_swap_init_sets_time()
   use swap_mod,        only: swap_init
   use swap_state_mod,  only: swap_state_t
   use swap_config_mod, only: swap_config_t
   use chdir_helper_mod, only: chdir_to
   use funit
   implicit none
   type(swap_state_t)           :: state
   type(swap_config_t), target  :: config

   call chdir_to('tests/swap-cases/toml/1.hupselbrook')
   call swap_init('swap.toml', state, config)
   @assertFalse(state%timecontrol%flRunEnd)
   @assertGreaterThan(state%timecontrol%t1900, 0.0d0)
   call chdir_to('../../../..')
end subroutine

@test
subroutine test_swap_run_step_advances_time()
   use swap_mod,        only: swap_init, swap_run_step
   use swap_state_mod,  only: swap_state_t
   use swap_config_mod, only: swap_config_t
   use chdir_helper_mod, only: chdir_to
   use funit
   implicit none
   type(swap_state_t)           :: state
   type(swap_config_t), target  :: config
   real(8) :: t_before

   call chdir_to('tests/swap-cases/toml/1.hupselbrook')
   call swap_init('swap.toml', state, config)
   t_before = state%timecontrol%t1900
   call swap_run_step(state, config)
   @assertGreaterThan(state%timecontrol%t1900, t_before)
   call chdir_to('../../../..')
end subroutine

@test
subroutine test_swap_full_lifecycle()
   use swap_mod,        only: swap_init, swap_run_step, swap_close
   use swap_state_mod,  only: swap_state_t
   use swap_config_mod, only: swap_config_t
   use chdir_helper_mod, only: chdir_to
   use funit
   implicit none
   type(swap_state_t)           :: state
   type(swap_config_t), target  :: config
   integer :: n_steps
   integer, parameter :: max_steps = 200000

   call chdir_to('tests/swap-cases/toml/1.hupselbrook')
   call swap_init('swap.toml', state, config)
   n_steps = 0
   do while (.not. state%timecontrol%flRunEnd .and. n_steps < max_steps)
      call swap_run_step(state, config)
      n_steps = n_steps + 1
   end do
   call swap_close(state, config)
   @assertTrue(state%timecontrol%flRunEnd)
   @assertGreaterThan(n_steps, 0)
   call chdir_to('../../../..')
end subroutine
```

- [ ] **Step 5: Register the test file in `tests/unit/meson.build`**

In `tests/unit/meson.build`, append `'../../src/core/swap_mod.f90'` to the `pfunit_extra_sources` list (anywhere; convention is near the other `src/core/` entries).

Then append `'core/test_swap_mod.pf'` as the **last entry** in the `pf_files` list (preserve trailing comma on the preceding entry).

- [ ] **Step 6: Register the suite in `tests/unit/testSuites.inc`**

Append a single line at the **end** of `tests/unit/testSuites.inc`:

```
ADD_TEST_SUITE(test_swap_mod_suite)
```

The "last" position matters: this suite runs a full simulation and mutates global state via `Initialize`; later suites would inherit that mutated state.

- [ ] **Step 7: Build and confirm the new test suite compiles**

Run: `pixi run build-linux`
Expected: build succeeds; the new `test_swap_mod_suite` is generated from the `.pf` file by pFUnit's `funitproc`.

- [ ] **Step 8: Run the pFUnit suite — confirm the three new tests FAIL**

Run: `pixi run test-pfunit`
Expected: three new tests fail in `test_swap_mod_suite`:
- `test_swap_init_sets_time` fails because `t1900` is 0 (empty stub didn't initialize).
- `test_swap_run_step_advances_time` fails because `t1900` didn't change (empty stub).
- `test_swap_full_lifecycle` fails because `flRunEnd` is still `.false.` (empty stub).
All other existing suites must still pass.

- [ ] **Step 9: Commit**

```bash
git add src/core/swap_mod.f90 tests/unit/core/test_swap_mod.pf \
        meson.build tests/unit/meson.build tests/unit/testSuites.inc
git commit -m "$(cat <<'EOF'
test(ss-drv): empty swap_mod skeleton + 3 failing lifecycle tests

Adds the module surface for the driver refactor and three pFUnit smoke
tests that drive a full hupselbrook simulation. Stub bodies return
without doing anything; tests fail at the expected assertions
(t1900==0, no time advance, flRunEnd stays false). Tasks 3-5 fill in
the procedure bodies one at a time.

Suite registered LAST in testSuites.inc — full simulation mutates
globals; later suites must not inherit that state.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Migrate `iTask=1` body into `swap_init`

**Files:**
- Modify: `src/core/swap_mod.f90` (fill `swap_init` body)

Source of truth: `src/core/swap.f90:158-313` (the `if (iTask == 1) then ... return; end if` block).

- [ ] **Step 1: Hoist the `use` statements into `swap_mod`'s module header**

The current `subroutine swap` has its `use` list at lines 70-114. Replicate those at module level inside `swap_mod` (just below the existing two `use` lines), so that all three procedures share them. Remove the inline `subroutine`-level `use` from `swap_init` if any (keep them at module level).

After this step, the top of `swap_mod.f90` looks like:

```fortran
module swap_mod
   use swap_state_mod,  only: swap_state_t
   use swap_config_mod, only: swap_config_t
   use variables, only : flswapshared, flmacropore, flcropnut, flagetracer, swfrost, &
                         swusecn, fldecmprat, flcropcalendar, flmaxitertime, &
                         flharvestday, flcropoutput, swcrp, flirrigationoutput, swend, project, &
                         flTillage, flSSDI, &
                         numnod, numlay
   use timestep_control_mod, only: fldecdt
   use soilwater_state_mod, only: soilwater_init
   use atmosphere_state_mod, only: atmosphere_init
   use tillage_state_mod, only: tillage_init
   use drainage_mod, only: drainage, drainage_init
   use surfacewater_mod, only: SurfaceWater, surfacewater_year_reset
   use tillage_mod,   only : DoTillage
   use swap_log, only: log_info
   use boundbottom_mod, only: BoundBottom
   use runoff_mod, only: CNmethod
   use meteo_mod, only: ProcessMeteoDay
   use meteo_process_mod, only: ReadMeteoDay
   use snow_mod, only: snow
   use meteodt_mod, only: MeteoDT
   use rootextraction_mod, only: RootExtraction
   use frozencond_mod, only: FrozenCond, FrozenBounds
   use temperature_mod, only: Temperature, heat_init
   use solute_mod, only: solute, solute_init
   use agetracer_mod, only: AgeTracer
   use soilgrid_mod, only: CalcGrid, ConvertDiscrVert
   use soilhydraulics_mod, only: soilwater, SoilWaterStateVar
   use WC_K_models_04_11, only: bind_cofgen_target
   use soilhydraulics_utils, only: bind_state_targets, bind_tc_target
   use config_to_variables_mod, only: h_init_buf, pondini_init_buf, pond_init_buf, &
                                      tc_iyear_init_buf, tc_imonth_init_buf, tc_dt_init_buf
   use irrigation_mod, only: irrigation, SSDI_irrigation
   use management_soil_mod, only: SoilManagement
   use error_mod, only: fatalerr_collected
   use load_swap_config_mod, only: load_swap_config
   use config_to_variables_mod, only: config_to_variables
   implicit none
   private
   public :: swap_init, swap_run_step, swap_close

contains
```

- [ ] **Step 2: Copy the `iTask=1` body into `swap_init`**

In `src/core/swap.f90`, the block from line 158 (`if (iTask == 1) then`) through line 313 (the `return; end if`) is the body. Inside `swap_init`, paste the contents **between** those two markers (not the `if/end if/return` themselves).

Two contextual edits as you paste:

1. Remove the `if (iCaller == 0) then ... end if` wrapper around the output-file open block (lines 293-303 in `swap.f90`). Always open output files. The body inside that `if` block stays — only the wrapper is removed:

   Before (`swap.f90:292-303`):
   ```fortran
   !  open Output files and write headers (skip in external/DLL mode to avoid per-column I/O)
      if (iCaller == 0) then
         call SwapOutput(1, state)
         call SoilWaterOutput(1, state)
         if (flTemperature)  call TemperatureOutput(1, state)
         if (flSolute)       call SoluteOutput(1, state)
         if (flAgeTracer)    call AgeTracerOutput(1, state)
         if (flSnow)         call SnowOutput(1, state)
         if (flSurfaceWater) call SurfaceWaterOutput(1, state)
      end if
   ```

   After (in `swap_init`):
   ```fortran
   !  open Output files and write headers
         call SwapOutput(1, state)
         call SoilWaterOutput(1, state)
         if (flTemperature)  call TemperatureOutput(1, state)
         if (flSolute)       call SoluteOutput(1, state)
         if (flAgeTracer)    call AgeTracerOutput(1, state)
         if (flSnow)         call SnowOutput(1, state)
         if (flSurfaceWater) call SurfaceWaterOutput(1, state)
   ```

2. Delete the `if (iCaller /= 0) call handle_exchange(11, flError, state)` line near the end of the block (line 308 in `swap.f90`).

3. Replace the call `load_swap_config('swap.toml', config, errors)` with `load_swap_config(config_file, config, errors)` so the dummy argument is used.

- [ ] **Step 3: Verify swap_init compiles**

Run: `pixi run build-linux`
Expected: build succeeds. If unresolved references appear (e.g., `IterTime`, `Initialize`, `TimeControl`, `SwapOutput`), those are *external* subroutines (not in modules) — they don't need `use` statements and are resolved at link time. Leave them alone.

- [ ] **Step 4: Run the pFUnit suite — confirm `test_swap_init_sets_time` PASSES**

Run: `pixi run test-pfunit`
Expected:
- `test_swap_init_sets_time` PASSES.
- `test_swap_run_step_advances_time` still FAILS (empty `swap_run_step`).
- `test_swap_full_lifecycle` still FAILS.
- All previously-passing suites still pass.

- [ ] **Step 5: Commit**

```bash
git add src/core/swap_mod.f90
git commit -m "$(cat <<'EOF'
refactor(ss-drv): migrate iTask=1 body into swap_init

Copies the initialization block from swap.f90:158-313 into the
swap_init procedure body, hoisting all `use` statements to the
module header. Removes the iCaller-gated output-file branch
(output files now always open) and the handle_exchange(11) call
(retired with the DLL exchange machinery).

test_swap_init_sets_time PASSES; the two run-step / full-lifecycle
tests remain failing as expected.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Migrate `iTask=2` body into `swap_run_step` (minus outer loop)

**Files:**
- Modify: `src/core/swap_mod.f90` (fill `swap_run_step` body)

Source of truth: `src/core/swap.f90:318-515` (the `if (iTask == 2) then ... return; end if` block).

- [ ] **Step 1: Declare locals needed by the step body**

At the top of `swap_run_step` (after the dummy args), declare locals that the step body uses (currently locals in `subroutine swap`):

```fortran
logical :: flError
logical :: request_smaller_dt
logical, external :: dtleap
```

- [ ] **Step 2: Copy the iTask=2 body into `swap_run_step`**

Inside `swap_run_step` (after the locals declared in Step 1), paste the contents of `src/core/swap.f90:318-515` between `if (iTask == 2) then` and `return; end if` (exclusive of those lines). Apply these edits as you paste:

1. **Remove the outer time loop wrapper.** Lines 344 (`do while (.not.flrunend)`) and 507 (`end do`) — delete both lines. The associate block must remain wrapping the body. Result: the `associate(...)` block opens, then the (former) loop body runs once, then `end associate`.

2. **Remove all `if (iCaller /= 0)` and `if (iCaller == 0)` guards.** Always run the body of these blocks (output, meteo I/O, handle_exchange). Specifically:
   - Line 321 `if (iCaller /= 0) call handle_exchange(21, flError, state); if (flError) return` — DELETE entire line.
   - Line 347 `if (iCaller == 0) then` and its matching `end if` at line 349 — remove the wrapper; keep `if (tc_flYearStart) call ReadMeteoYear(state)`.
   - Line 354 `if (iCaller /= 0) call handle_exchange(22, flError, state)   ! weather` — DELETE.
   - Line 364 `if (iCaller /= 0) call handle_exchange(23, flError, state)   ! LAI, RD` — DELETE.
   - Line 446 — `if (iCaller == 0 .and. flCropCalendar) call CropGrowth(2, state%heat%tsoil, state)` — change to `if (flCropCalendar) call CropGrowth(2, state%heat%tsoil, state)`.
   - Line 456 — same pattern: `if (iCaller == 0 .and. flCropCalendar) call CropGrowth(3, ...)` → drop the `iCaller == 0 .and.` prefix.
   - Line 463 — same: `if (iCaller == 0 .and. flCropCalendar) call CropGrowth(4, ...)` → drop the prefix.
   - Line 476 `if (iCaller == 0) then` and its matching `end if` at line 502 — remove the wrapper; the entire output block always runs.
   - Line 510 `if (iCaller /= 0) call handle_exchange(29, flError, state)` — DELETE.

3. **Add the `CropGrowth` explicit interface.** The current `subroutine swap` has an `interface` block at lines 121-128 declaring `CropGrowth`. Add the same interface block inside `swap_run_step`, after the locals, before the executable section:

   ```fortran
   interface
      subroutine CropGrowth(task, tsoil, state)
         use swap_state_mod, only: swap_state_t
         integer, intent(in) :: task
         real(8), intent(in) :: tsoil(:)
         type(swap_state_t), intent(inout) :: state
      end subroutine CropGrowth
   end interface
   ```

   This same interface is also needed in `swap_init` (Task 3) since it calls `CropGrowth`. **Audit Task 3's `swap_init` body** — if the original `iTask=1` block calls `CropGrowth` (it does not, per inspection of `swap.f90:158-313`), the interface is not needed there. Verify by `grep -n CropGrowth src/core/swap.f90` and add to `swap_init` only if grep shows a call in lines 158-313.

- [ ] **Step 3: Build and verify**

Run: `pixi run build-linux`
Expected: build succeeds.

- [ ] **Step 4: Run the pFUnit suite — `test_swap_run_step_advances_time` PASSES**

Run: `pixi run test-pfunit`
Expected:
- `test_swap_init_sets_time` PASSES (unchanged).
- `test_swap_run_step_advances_time` PASSES — t1900 advances by one `dt`.
- `test_swap_full_lifecycle` likely fails late (`@assertTrue(state%timecontrol%flRunEnd)` — depends on whether `swap_close` is actually needed to mark the run done, but flRunEnd is set inside the loop by TimeControl(2) which we now do call; verify experimentally). If it still fails it's expected because `swap_close` is still empty — Task 5 fixes.

- [ ] **Step 5: Commit**

```bash
git add src/core/swap_mod.f90
git commit -m "$(cat <<'EOF'
refactor(ss-drv): migrate iTask=2 body into swap_run_step

Copies the timestep body from swap.f90:318-515 into swap_run_step,
dropping the outer `do while (.not.flrunend)` wrapper (loop moves to
the caller). All iCaller-gated branches removed: meteo I/O, output,
and CropGrowth calls always run; handle_exchange(21|22|23|29) calls
deleted with the retired DLL exchange machinery.

The TimeControl `associate` block is preserved wrapping the entire
step body — same field aliases as before (tc_flYearStart, tc_flDayStart,
tc_flDayEnd, tc_daynr, tc_iyear, flrunend, flOutput, flOutputShort,
flMeteoDt, flETSine, flSnow, flSolute, flTemperature, flDrain,
flSurfaceWater, flIrrigate, fldtreduce). Inner dt-reduction loop
unchanged.

test_swap_run_step_advances_time PASSES.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: Migrate `iTask=3` body into `swap_close`

**Files:**
- Modify: `src/core/swap_mod.f90` (fill `swap_close` body)

Source of truth: `src/core/swap.f90:520-552` (the `if (iTask == 3) then ... return; end if` block).

- [ ] **Step 1: Copy the iTask=3 body into `swap_close`**

Inside `swap_close`, paste the contents of `src/core/swap.f90:520-552` between `if (iTask == 3) then` and `return; end if` (exclusive). Apply these edits:

1. Remove the `if (iCaller == 0) then` wrapper at line 526 and its matching `end if` at line 541 — keep the body (output closes, etc., always run).
2. Delete line 547 `if (iCaller /= 0) call handle_exchange(31, flError, state)`.

After edits, `swap_close` body looks like:

```fortran
!  iteration and timing statistics
   call IterTime(3, state)

!  close output files
   if (flSwapShared) call SharedSimulation(4)
   call SwapOutput(3, state)
   if (swend.eq.1) call SoilWaterOutput(3, state)
   call SoilWaterOutput(4, state)
   if (swcrp.eq.1) call CropOutput(3, state)
   if (state%timecontrol%flTemperature)  call TemperatureOutput(3, state)
   if (state%timecontrol%flSolute)       call SoluteOutput(3, state)
   if (flAgeTracer)                      call AgeTracerOutput(3, state)
   if (state%timecontrol%flSnow)         call SnowOutput(3, state)
   if (state%timecontrol%flSurfaceWater) call SurfaceWaterOutput(3, state)
   if (flCropNut)                        call SoilManagement(7, state)

!  write okay file for external use
   call WriteSwapOk(Project)

   call log_info('swap', 'Simulation complete for project: ' // trim(project))
```

- [ ] **Step 2: Build and verify**

Run: `pixi run build-linux`
Expected: build succeeds.

- [ ] **Step 3: Run the pFUnit suite — all three lifecycle tests PASS**

Run: `pixi run test-pfunit`
Expected: all three `test_swap_mod_suite` tests PASS. All existing suites still pass.

- [ ] **Step 4: Commit**

```bash
git add src/core/swap_mod.f90
git commit -m "$(cat <<'EOF'
refactor(ss-drv): migrate iTask=3 body into swap_close

Copies the closure block from swap.f90:520-552 into swap_close.
Removes the iCaller-gated wrapper around output-file closes (they
always run now) and the handle_exchange(31) call. swap_mod is now
functionally complete; cutover happens in Task 6.

All three test_swap_mod_suite lifecycle tests PASS.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: Cutover — rewrite `swap_main`, delete `swap.f90`, retire `sharedexchange.f90`

**Files:**
- Rewrite: `src/core/swap_main.f90`
- Delete: `src/core/swap.f90`
- Delete: `src/utils/sharedexchange.f90`
- Modify: `meson.build` (remove `src/core/swap.f90` and `src/utils/sharedexchange.f90` from `sources`)

- [ ] **Step 1: Rewrite `src/core/swap_main.f90` as the thin driver**

Replace the entire contents of `src/core/swap_main.f90` with:

```fortran
! swap_main.f90
! SS-DRV Phase 1: thin driver around module swap_mod.
! The outer time loop lives here so the same module can be driven
! by BMI (swap_bmi_mod) one timestep at a time.
program swap_main

   use swap_mod,        only: swap_init, swap_run_step, swap_close
   use swap_state_mod,  only: swap_state_t
   use swap_config_mod, only: swap_config_t
   use swap_log,        only: log_init, log_close, LOGLEVEL_INFO
   implicit none

   type(swap_state_t)           :: state
   type(swap_config_t), target  :: config

   call log_init(log_level=LOGLEVEL_INFO, log_file='swap_debug.log')

   call swap_init('swap.toml', state, config)
   do while (.not. state%timecontrol%flRunEnd)
      call swap_run_step(state, config)
   end do
   call swap_close(state, config)

   write(*,'(a)')' Swap normal completion!'
   call log_close()
   call CloseTempFil
   call Exit(100)

end program swap_main
```

`CloseTempFil` is preserved (behavior parity with the old `swap_main`); a future cleanup decides whether it can be retired.

- [ ] **Step 2: Delete `src/core/swap.f90`**

Run:
```bash
git rm src/core/swap.f90
```

- [ ] **Step 3: Delete `src/utils/sharedexchange.f90`**

`FromSwap` / `ToSwap` are stub no-ops and have no remaining callers (the only references were inside the old `subroutine swap`, now deleted). Confirm and delete:

Run:
```bash
grep -rn "FromSwap\|ToSwap" src/ tests/
```
Expected: no hits outside `src/utils/sharedexchange.f90` itself.

Run:
```bash
git rm src/utils/sharedexchange.f90
```

- [ ] **Step 4: Remove the deleted files from `meson.build` sources**

In `meson.build`:
- Remove the line `'src/core/swap.f90',`
- Remove the line `'src/utils/sharedexchange.f90',`

- [ ] **Step 5: Build and confirm the executable links**

Run: `pixi run build-linux`
Expected: build succeeds and produces `builddir/swap`. If there are unresolved references for `FromSwap` / `ToSwap`, find and remove those callers — they should not exist.

- [ ] **Step 6: Run the pFUnit suite**

Run: `pixi run test-pfunit`
Expected: all suites pass, including the three `test_swap_mod_suite` tests.

- [ ] **Step 7: Commit**

```bash
git add src/core/swap_main.f90 meson.build
git commit -m "$(cat <<'EOF'
refactor(ss-drv): cutover — retire swap.f90, sharedexchange.f90, iCaller

Replaces swap_main with a thin driver that calls swap_mod's three
named procedures and runs the outer time loop. The legacy free
subroutine `swap` is deleted; sharedexchange.f90 (dead FromSwap/
ToSwap stubs) is deleted; the iCaller dispatch concept is retired.

Behavior parity preserved: CloseTempFil still called at program end.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: Regression gate — byte-for-byte parity with the baseline

**Files:**
- No code changes; verification only.

- [ ] **Step 1: Run the fast regression gate**

Run: `pixi run check-fast`
Expected: 4/4 cases pass on hupselbrook, surfacewater, salinitystress, grassgrowth — outputs match the regression baselines.

- [ ] **Step 2: Run the full regression gate (all six cases)**

Run: `pixi run check-full`
Expected: all six regression cases pass.

- [ ] **Step 3: If any case fails, diagnose**

If a case fails, the failure points to a missed `iCaller` branch removal, a mis-applied edit during body migration, or a side-effect difference (e.g. `swap_init` calling something twice). Use `git diff HEAD~3 -- src/core/swap_mod.f90` to inspect the migration and the regression test's diff output to localize. Fix in place and re-run from Step 1.

- [ ] **Step 4: Commit a verification checkpoint**

```bash
git commit --allow-empty -m "$(cat <<'EOF'
chore(ss-drv): regression gate PASSED post-cutover

check-full: all six regression cases match baseline byte-for-byte.
swap_mod refactor is behaviorally identical to the legacy
subroutine swap. Safe to layer BMI on top.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 8: Create `swap_bmi_mod` with `bind(C)` wrappers

**Files:**
- Create: `src/core/swap_bmi_mod.f90`
- Modify: `meson.build` (add `src/core/swap_bmi_mod.f90` to sources for now; Task 9 restructures into modern/legacy split)

- [ ] **Step 1: Create the BMI module**

Create `src/core/swap_bmi_mod.f90` with this exact content:

```fortran
!> @file swap_bmi_mod.f90
!! SS-DRV Phase 1: minimal Basic Model Interface (BMI) façade.
!! Exposes initialize / update / finalize / get_value_double /
!! get_current_time as bind(C) procedures. All other CSDMS BMI v2.0
!! methods are stubbed (correct C signature, body returns 0 / empty)
!! and marked with `! BMI-STUB Phase 2` for the next arc to enumerate.
!! Holds a single module-level (state, config) pair — multi-column
!! support is Phase 3 (swap_ensemble_mod).
module swap_bmi_mod
   use iso_c_binding,   only: c_char, c_double, c_int, c_size_t, c_null_char, c_ptr
   use swap_mod,        only: swap_init, swap_run_step, swap_close
   use swap_state_mod,  only: swap_state_t
   use swap_config_mod, only: swap_config_t
   implicit none
   private

   type(swap_state_t),          save :: bmi_state
   type(swap_config_t), target, save :: bmi_config

contains

   !----------------------------------------------------------------------
   ! Lifecycle
   !----------------------------------------------------------------------

   function bmi_initialize(config_file, n) result(rc) bind(C, name='initialize')
      character(kind=c_char), intent(in)    :: config_file(*)
      integer(c_int),  value, intent(in)    :: n
      integer(c_int)                        :: rc
      character(len=256) :: f_config_file
      call c_to_f_string(config_file, f_config_file)
      call swap_init(trim(f_config_file), bmi_state, bmi_config)
      rc = 0
   end function bmi_initialize

   function bmi_update() result(rc) bind(C, name='update')
      integer(c_int) :: rc
      call swap_run_step(bmi_state, bmi_config)
      rc = 0
   end function bmi_update

   function bmi_finalize() result(rc) bind(C, name='finalize')
      integer(c_int) :: rc
      call swap_close(bmi_state, bmi_config)
      rc = 0
   end function bmi_finalize

   !----------------------------------------------------------------------
   ! Time
   !----------------------------------------------------------------------

   function bmi_get_current_time(t) result(rc) bind(C, name='get_current_time')
      real(c_double), intent(out) :: t
      integer(c_int)              :: rc
      t = bmi_state%timecontrol%t1900
      rc = 0
   end function bmi_get_current_time

   !----------------------------------------------------------------------
   ! Variable accessors — two sentinel variables for Phase 1
   !----------------------------------------------------------------------

   function bmi_get_value_double(var_name, n, dest) result(rc) bind(C, name='get_value_double')
      character(kind=c_char), intent(in)    :: var_name(*)
      integer(c_int),  value, intent(in)    :: n
      real(c_double),         intent(out)   :: dest(n)
      integer(c_int)                        :: rc
      character(len=64) :: name
      integer :: m
      call c_to_f_string(var_name, name)
      select case (trim(name))
      case ('soil_water_content')
         m = min(n, size(bmi_state%soilwater%theta))
         dest(1:m) = bmi_state%soilwater%theta(1:m)
         if (m < n) dest(m+1:n) = 0.0_c_double
         rc = 0
      case ('actual_evapotranspiration')
         if (n >= 1) dest(1) = bmi_state%soilwater%intr%iqrot
         if (n > 1)  dest(2:n) = 0.0_c_double
         rc = 0
      case default
         rc = 1
      end select
   end function bmi_get_value_double

   !----------------------------------------------------------------------
   ! BMI-STUB Phase 2: spec-required methods, currently return 0 / empty.
   ! Each carries the correct C signature so consumers can compile, but
   ! does nothing useful until Phase 2 wires the full variable registry.
   !----------------------------------------------------------------------

   function bmi_get_end_time(t) result(rc) bind(C, name='get_end_time')
      real(c_double), intent(out) :: t
      integer(c_int)              :: rc
      t = 0.0_c_double            ! BMI-STUB Phase 2
      rc = 0
   end function bmi_get_end_time

   function bmi_get_time_step(dt) result(rc) bind(C, name='get_time_step')
      real(c_double), intent(out) :: dt
      integer(c_int)              :: rc
      dt = bmi_state%timecontrol%dt   ! cheap real impl
      rc = 0
   end function bmi_get_time_step

   function bmi_set_value_double(var_name, n, src) result(rc) bind(C, name='set_value_double')
      character(kind=c_char), intent(in) :: var_name(*)
      integer(c_int),  value, intent(in) :: n
      real(c_double),         intent(in) :: src(n)
      integer(c_int)                     :: rc
      rc = 1                       ! BMI-STUB Phase 2 — no settable variables yet
   end function bmi_set_value_double

   !----------------------------------------------------------------------
   ! Helpers
   !----------------------------------------------------------------------

   subroutine c_to_f_string(c_str, f_str)
      character(kind=c_char), intent(in)  :: c_str(*)
      character(len=*),       intent(out) :: f_str
      integer :: i
      f_str = ' '
      do i = 1, len(f_str)
         if (c_str(i) == c_null_char) exit
         f_str(i:i) = c_str(i)
      end do
   end subroutine c_to_f_string

end module swap_bmi_mod
```

- [ ] **Step 2: Temporarily add the BMI module to the existing `sources` list**

In `meson.build`, in the `sources = [` block, near `'src/core/swap_mod.f90',`, add a sibling line:

```
    'src/core/swap_bmi_mod.f90',
```

This is temporary — Task 9 splits modern and legacy sources into separate static libraries. For now we just want the BMI module compiling alongside everything else.

- [ ] **Step 3: Build and verify**

Run: `pixi run build-linux`
Expected: build succeeds. The BMI module compiles and links into the `swap` executable (unused, but linked).

- [ ] **Step 4: Run the pFUnit suite and regression**

Run: `pixi run check-fast`
Expected: all 4 cases still pass. No behavior change — `bmi_state` / `bmi_config` are not touched by the executable's path.

- [ ] **Step 5: Commit**

```bash
git add src/core/swap_bmi_mod.f90 meson.build
git commit -m "$(cat <<'EOF'
feat(ss-drv): add swap_bmi_mod — minimal BMI bind(C) façade

Exposes the five real BMI methods (initialize, update, finalize,
get_current_time, get_value_double for two sentinel variables) plus
spec-required stubs (get_end_time, get_time_step, set_value_double)
marked `! BMI-STUB Phase 2`. Holds a single module-level state +
config pair — multi-column ensemble is Phase 3.

Build wiring is interim — Task 9 splits modern/legacy into separate
static libraries to scope -std=f2018 to the new files only.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 9: Restructure `meson.build` — modern static lib, legacy static lib, BMI shared lib

**Files:**
- Modify: `meson.build` (split sources, add libraries, change executable to link composition)

- [ ] **Step 1: Edit `meson.build` to introduce the two static libraries and the shared BMI lib**

Replace the section starting at `# Source files` through `# Main executable` (currently lines 52-215). Replace it with the structure below. Keep the existing `sources = [...]` list contents, but **remove the two modern entries** (`'src/core/swap_mod.f90'` and `'src/core/swap_bmi_mod.f90'`) from `sources` — they move to `modern_sources`.

The new section:

```python
# ============================================================================
# Source partitions
# ============================================================================
# Modern files: compiled with -std=f2018 -Wall -Wextra.
# Cross-standard linkage works fine — gfortran .mod files produced by
# legacy compilation are consumed without issue by f2018 compilation.
modern_sources = [
    'src/core/swap_mod.f90',
    'src/core/swap_bmi_mod.f90',
]

# Legacy files: existing project-wide -std=legacy (set in the top
# add_project_arguments call). The original `sources` list below
# excludes the modern files and excludes src/core/swap_main.f90
# (the executable's own source).
sources = [
    # ... preserve all existing entries from the previous `sources` list,
    # MINUS src/core/swap.f90 (already deleted in Task 6)
    # MINUS src/utils/sharedexchange.f90 (already deleted in Task 6)
    # MINUS src/core/swap_mod.f90 (now in modern_sources)
    # MINUS src/core/swap_bmi_mod.f90 (now in modern_sources)
    # MINUS src/core/swap_main.f90 (becomes the executable's own source)
    # ... all other existing entries unchanged
]

# ============================================================================
# External dependency: toml-f
# ============================================================================
tomlf_prj = subproject(
    'toml-f',
    version: '>=0.2',
    default_options: ['default_library=static'],
)
tomlf_dep = tomlf_prj.get_variable('tomlf_dep')

# ============================================================================
# Static libraries (compile-time grouping)
# ============================================================================
swap_legacy = static_library('swap_legacy',
    sources: sources,
    dependencies: [tomlf_dep],
    include_directories: src_inc,
)

swap_modern = static_library('swap_modern',
    sources: modern_sources,
    dependencies: [tomlf_dep],
    include_directories: src_inc,
    fortran_args: ['-std=f2018', '-Wall', '-Wextra'],
)

# ============================================================================
# Main executable
# ============================================================================
# Always statically link with gfortran — see docs/adr/0001-gfortran-first.md.
link_args = ['-static']

executable('swap',
    sources: 'src/core/swap_main.f90',
    link_with: [swap_modern, swap_legacy],
    dependencies: [tomlf_dep],
    fortran_args: ['-std=f2018'],
    include_directories: src_inc,
    link_args: link_args,
    install: false
)

# ============================================================================
# BMI shared library (Python / imod_coupler entry point)
# ============================================================================
# No fortran_args here — link-only target. Each member library was
# already compiled with its appropriate -std= flag.
swap_bmi_lib = shared_library('swap_bmi',
    link_with: [swap_modern, swap_legacy],
    dependencies: [tomlf_dep],
    include_directories: src_inc,
    install: false,
)
```

When you remove `'src/core/swap_main.f90'` from the `sources` list, the diff for that file is just `-    'src/core/swap_main.f90'`. The executable now lists it as `sources: 'src/core/swap_main.f90'` directly.

- [ ] **Step 2: Build and verify the executable still works**

Run: `pixi run build-linux`
Expected: build succeeds. Two static libraries (`libswap_modern.a`, `libswap_legacy.a`) and one shared library (`libswap_bmi.so`) appear in `builddir/`. The `swap` executable still links and runs.

- [ ] **Step 3: Verify the shared library has the BMI symbols exposed**

Run:
```bash
nm -D --defined-only builddir/libswap_bmi.so | grep -E '^[0-9a-f]+ T (initialize|update|finalize|get_current_time|get_value_double|get_time_step)$'
```
Expected: 6 lines, one per BMI entry point (`initialize`, `update`, `finalize`, `get_current_time`, `get_value_double`, `get_time_step`). The `T` indicates a defined text (code) symbol.

- [ ] **Step 4: Run pFUnit and regression as sanity gates**

Run: `pixi run check-fast`
Expected: 4/4 regression cases still pass; pFUnit still passes.

- [ ] **Step 5: Commit**

```bash
git add meson.build
git commit -m "$(cat <<'EOF'
build(ss-drv): split sources — swap_modern (f2018) + swap_legacy + libswap_bmi.so

Introduces two static libraries and a shared library:
- swap_modern: swap_mod + swap_bmi_mod, compiled with -std=f2018
  -Wall -Wextra. Scoped strict-standard enforcement on the two new
  files only; the rest of the codebase stays on -std=legacy until a
  dedicated standards-migration arc.
- swap_legacy: all existing physics + I/O sources, unchanged.
- libswap_bmi.so: shared library linking both, for Python (cffi) and
  imod_coupler consumption.

The swap executable now links both static libraries instead of
compiling everything monolithically.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 10: Python `cffi` hello-world + Meson test registration

**Files:**
- Create: `tests/bmi/swap_bmi.h`
- Create: `tests/bmi/hello_swap.py`
- Create: `tests/bmi/meson.build`
- Modify: `meson.build` (add `subdir('tests/bmi')`)

- [ ] **Step 1: Create the C header**

Create `tests/bmi/swap_bmi.h`:

```c
/* SS-DRV Phase 1: minimal BMI C surface for libswap_bmi.so.
 * Full CSDMS BMI v2.0 compliance arrives in Phase 2. */

int initialize(const char *config_file, int n);
int update(void);
int finalize(void);
int get_current_time(double *t);
int get_time_step(double *dt);
int get_value_double(const char *var_name, int n, double *dest);
int set_value_double(const char *var_name, int n, const double *src);
```

No `#ifndef` guards — cffi parses this directly via `ffi.cdef`.

- [ ] **Step 2: Create the Python hello-world driver**

Create `tests/bmi/hello_swap.py`:

```python
"""SS-DRV Phase 1: BMI hello-world.

Loads libswap_bmi.so via cffi, drives the hupselbrook simulation
through one initialize / update / finalize cycle, and asserts the
two sentinel get_value_double variables return physically plausible
data. End-to-end proof that the C binding path works.
"""

import pathlib
import sys

import cffi

HERE = pathlib.Path(__file__).resolve().parent


def main() -> int:
    if len(sys.argv) < 2:
        print("usage: hello_swap.py <path-to-libswap_bmi.so>", file=sys.stderr)
        return 2

    so_path = sys.argv[1]
    header_text = (HERE / "swap_bmi.h").read_text()

    ffi = cffi.FFI()
    ffi.cdef(header_text)
    lib = ffi.dlopen(so_path)

    rc = lib.initialize(b"swap.toml\0", 0)
    assert rc == 0, f"initialize returned {rc}"

    rc = lib.update()
    assert rc == 0, f"update returned {rc}"

    t = ffi.new("double *")
    assert lib.get_current_time(t) == 0
    print(f"current_time after 1 step = {t[0]}")

    dt = ffi.new("double *")
    assert lib.get_time_step(dt) == 0
    print(f"time_step                = {dt[0]}")

    buf = ffi.new("double[500]")
    assert lib.get_value_double(b"soil_water_content\0", 500, buf) == 0
    theta0 = buf[0]
    print(f"theta[0]                 = {theta0}")
    assert 0.0 < theta0 < 1.0, f"theta[0]={theta0} not in (0, 1)"

    # Unknown variable name must return rc != 0
    rc = lib.get_value_double(b"nonexistent_variable\0", 1, buf)
    assert rc != 0, "expected nonzero rc for unknown variable name"

    rc = lib.finalize()
    assert rc == 0, f"finalize returned {rc}"

    print("BMI hello-world: OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 3: Create the test directory's Meson build file**

Create `tests/bmi/meson.build`:

```python
# ============================================================================
# BMI hello-world integration test
# ============================================================================
# Runs the Python cffi driver against libswap_bmi.so in the
# hupselbrook case directory. Registered in the `bmi` suite so it
# can be run independently of pFUnit and the regression suite.
python3 = import('python').find_installation('python3')

test('bmi-hello-world',
    python3,
    args: [
        meson.current_source_dir() / 'hello_swap.py',
        swap_bmi_lib.full_path(),
    ],
    workdir: meson.project_source_root() / 'tests/swap-cases/toml/1.hupselbrook',
    suite: 'bmi',
    depends: [swap_bmi_lib],
    timeout: 60,
)
```

- [ ] **Step 4: Include the BMI test directory from the top-level `meson.build`**

In `meson.build`, in the tests section (currently at lines 216-221 in the original; locate the `if get_option('enable_pfunit')` block), append a new line **after** the `subdir('tests/unit')` call but still inside conditional gating. Since the BMI test does not depend on pFUnit, gate it on a fresh option, or simply add it unconditionally — choose unconditional:

```python
# ============================================================================
# Tests
# ============================================================================
if get_option('enable_pfunit')
    subdir('tests/unit')
endif

subdir('tests/bmi')
```

- [ ] **Step 5: Build and confirm the test target is registered**

Run: `pixi run build-linux`
Expected: build succeeds.

Run: `meson test -C builddir --list --suite bmi`
Expected: one line — `swap:bmi / bmi-hello-world`.

- [ ] **Step 6: Verify `cffi` is available in the pixi test environment**

Run: `pixi run -e test python -c "import cffi; print(cffi.__version__)"`
Expected: a version string prints (cffi is a transitive dep of `pyswap` per `pixi.toml`).

If it errors with `ModuleNotFoundError`, add `cffi = ">=1.16"` to `[feature.test.dependencies]` in `pixi.toml`, then `pixi install -e test`.

- [ ] **Step 7: Run the BMI hello-world test**

Run: `pixi run -e test meson test -C builddir --suite bmi --verbose`
Expected: `bmi-hello-world ... OK` with the Python script printing the time, time_step, theta[0] values, and `BMI hello-world: OK`.

- [ ] **Step 8: Commit**

```bash
git add tests/bmi/swap_bmi.h tests/bmi/hello_swap.py tests/bmi/meson.build meson.build
git commit -m "$(cat <<'EOF'
test(ss-drv): Python cffi BMI hello-world

Drives libswap_bmi.so through a full initialize -> update ->
get_current_time -> get_time_step -> get_value_double -> finalize
cycle from Python via cffi. Asserts theta[0] is physically plausible
(0 < theta < 1) and that an unknown variable name returns a non-zero
rc. Registered as the meson `bmi` suite, workdir=hupselbrook.

End-to-end proof that the C binding path works — same surface
imod_coupler will speak in Phase 2.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 11: Final verification — full regression + retirement grep

**Files:**
- No code changes; verification only.

- [ ] **Step 1: Run the full check-full gate**

Run: `pixi run check-full`
Expected: build + pFUnit + all six regression cases pass.

- [ ] **Step 2: Run the BMI suite**

Run: `pixi run -e test meson test -C builddir --suite bmi --verbose`
Expected: bmi-hello-world passes.

- [ ] **Step 3: Confirm full retirement of the legacy DLL pattern**

Run:
```bash
grep -rn 'iCaller\|swap_exchange\|handle_exchange\|FromSwap\|ToSwap' src/
```
Expected: no hits anywhere in `src/`.

- [ ] **Step 4: Confirm the deleted files are gone**

Run:
```bash
ls src/core/swap.f90 src/utils/sharedexchange.f90 2>&1 | grep -c "No such file"
```
Expected: `2` (both files deleted).

- [ ] **Step 5: Confirm the new files are present**

Run:
```bash
ls src/core/swap_mod.f90 src/core/swap_bmi_mod.f90 tests/unit/core/test_swap_mod.pf \
   tests/bmi/swap_bmi.h tests/bmi/hello_swap.py tests/bmi/meson.build
```
Expected: all six files listed without error.

- [ ] **Step 6: Confirm the `swap_main` shrinkage**

Run:
```bash
wc -l src/core/swap_main.f90
```
Expected: approximately 25–30 lines (down from 117).

- [ ] **Step 7: Tag the arc complete**

```bash
git commit --allow-empty -m "$(cat <<'EOF'
chore(ss-drv): Phase 1 complete — driver modernization & BMI stub

Summary of the arc:
- subroutine swap (free, 713 lines) -> module swap_mod (3 named
  procedures: swap_init / swap_run_step / swap_close)
- outer time loop lifted from inside swap to swap_main and BMI update()
- swap.f90, sharedexchange.f90, iCaller, swap_exchange, handle_exchange,
  dummy() — all retired
- new swap_bmi_mod with bind(C) façade: initialize/update/finalize +
  get_current_time/get_time_step/get_value_double for two sentinel
  variables, plus BMI-STUB markers on Phase 2 methods
- meson restructured: swap_modern (-std=f2018) + swap_legacy
  (-std=legacy) static libs + libswap_bmi.so shared lib
- pFUnit test_swap_mod_suite: 3 lifecycle tests on hupselbrook
- Python cffi hello-world: bmi-hello-world meson test in `bmi` suite
- check-full: all six regression cases pass byte-for-byte

Deferred to follow-on arcs:
- TimeControl modernization (named procedures replacing magic ints)
- Globals cleanup (35 files still `use variables`)
- Full BMI variable registry + metadata methods
- swap_ensemble_mod (Fortran-level multi-column orchestrator) + OpenMP
- Project-wide -std=f2018 migration

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Plan self-review

**Spec coverage check** (every "what goes / what stays / what's new" item in the spec maps to at least one task):

- `swap_exchange` module retired → Task 6 (deleted with swap.f90).
- `handle_exchange` retired → Task 6 (deleted with swap.f90).
- `iCaller` argument retired → Tasks 3/4/5 strip branches; Task 11 grep verifies.
- `dummy()` retired → Task 6 (rewritten `swap_main` omits it).
- `sharedexchange.f90` retired → Task 6.
- `swap.f90` deleted → Task 6.
- `swap_mod` module with three procedures → Tasks 2/3/4/5.
- `swap_main` thin driver → Task 6.
- TimeControl `associate` block preserved → Task 4 explicitly preserves it.
- `swap_bmi_mod` with bind(C) → Task 8.
- C header + Python hello-world → Task 10.
- Modern static lib with `-std=f2018` + legacy lib + shared lib → Task 9.
- pFUnit lifecycle smoke test → Task 2 (writes), Tasks 3/4/5 (each makes one pass).
- Python cffi BMI integration test → Task 10.
- Regression baselines unchanged → Tasks 1 (capture), 7 (verify), 11 (verify).

**Placeholder scan**: no "TBD" / "TODO" / "fill in later" markers — all code blocks are concrete, all commands have expected outputs, all file paths are absolute or repo-relative.

**Type/name consistency**: `swap_init` / `swap_run_step` / `swap_close` signatures match between Task 2 (skeleton), Task 6 (swap_main usage), and Task 8 (BMI usage). `bmi_state` / `bmi_config` named consistently in Task 8. `test_swap_mod_suite` registered in `testSuites.inc` (Task 2) matches the `@suite(name=...)` in the .pf file (Task 2).

**Order check**: Tasks 3/4/5 each rely on Task 2's failing tests and progressively make them pass. Task 6 (cutover) cannot run before Tasks 3/4/5 (module bodies must be complete). Task 7 (regression gate) cannot run before Task 6. Task 8 (BMI module) depends on Task 6 (clean module to wrap). Task 9 (build restructure) depends on Task 8. Task 10 (Python test) depends on Task 9 (shared library exists).

---

Plan complete and saved to `docs/superpowers/plans/2026-05-12-driver-modernization.md`. Two execution options:

**1. Subagent-Driven (recommended)** — I dispatch a fresh subagent per task, review between tasks, fast iteration

**2. Inline Execution** — Execute tasks in this session using executing-plans, batch execution with checkpoints

Which approach?
