# SWAP → MODFLOW 6 Coupled Demo Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Produce a runnable example of SWAP coupled to MODFLOW 6 (two-channel lateral-flow cross-section, recharge per cell from N homogeneous SWAP columns) driven by imod_coupler's `swapmod` driver, and fold in the agreed quick cleanups — all with standalone SWAP regression staying byte-identical.

**Architecture:** One `libswap.so` exposes an XMI shared-library kernel wrapping an *ensemble* of N independent `swap_state_t` columns sharing one read-only `swap_config_t` (per-column config indirection built in, M=1 for the demo). The XMI `solve` loops `swap_run_step` over columns sequentially, injecting each cell's MODFLOW head as the prescribed-GWL bottom boundary (`gwlinp`, cm) and reading back the bottom flux (`qbot`, cm/d) as recharge. A flopy-built MODFLOW 6 model and a thin fork of imod_coupler's `swapmod` driver + `swap_wrapper` (SWAP-native exchange names) complete the loop. MODFLOW 6 leads the clock.

**Tech Stack:** Fortran (gfortran, `iso_fortran_env`, `iso_c_binding`, meson), pFUnit, Python (ctypes/xmipy, cffi for SWAP-only smoke), flopy, imod_coupler (swap branch fork), pixi.

**Spec:** `dev-docs/superpowers/specs/2026-05-31-swap-modflow-coupling-design.md`

**Working repo:** `/home/zawadzkim/Code/swap` (branch `development`; this work lands on a short-lived branch off `development`).

---

## Conventions used throughout

- **Byte-identical law:** after every Fortran task, run `pixi run -e test check-fast` and confirm it stays green. The XMI/ensemble layer is *additive* — the standalone CLI path (`swap_main` → `swap_run_step`) is never modified, so regression output cannot change. At milestones run `pixi run -e test check-full`.
- **Clean rebuild rule:** any change to `src/state/*.f90` schema requires `pixi run clean` before rebuild (stale `.mod` SIGSEGV). Tasks that touch state say so explicitly.
- **pFUnit registration:** every new `.pf` file MUST be added to `tests/unit/testSuites.inc`; confirm the `OK (N tests)` count rises.
- **Commit discipline:** one cohesive change per commit; conventional prefixes; never mix refactor with behavior change. End commit messages with the `Co-Authored-By` trailer.
- **Exchange-array unit conventions (locked here, referenced everywhere):**
  - `gwl(i)` — MODFLOW head in **metres**, datum = SWAP surface = MODFLOW `top` = 0. Injected as `columns(i)%soilwater%gwlinp = gwl(i) * 100.0` (m→cm; below-surface is negative).
  - `qbot_volume(i)` — recharge **depth in metres over the step** = `(-columns(i)%soilwater%qbot) * 0.01 * dt_days`. Sign flip: SWAP `qbot < 0` = downward percolation = positive recharge. The driver computes `mf6_recharge = qbot_volume / delt` → m/d rate for RCHA (no area factor — matches the swapmod driver's `[:] = [:]/delt`).
  - `storage_coef(i)` — dimensionless specific yield for MF6 `sc1`; fixed constant `0.15` for the smoke demo (documented placeholder).

---

## Branch setup

- [ ] **Step 1: Create the work branch**

Run:
```bash
cd /home/zawadzkim/Code/swap
git checkout development && git pull --ff-only 2>/dev/null; git checkout -b feat/swap-modflow-coupling
```

- [ ] **Step 2: Confirm baseline green**

Run: `pixi run -e test check-fast`
Expected: pFUnit `OK (N tests)` + 4 regression cases PASS.

---

# PHASE 0 — Cleanups (independent; do first to de-risk; each byte-identical)

### Task 0.1: Drop MetaSWAP

**Files:**
- Delete: `src/utils/dormant/sharedsimulation.f90`
- Modify: `src/atmosphere/interception.f90` (remove `msw1eic` + its wrapper)
- Modify: `src/core/swap_log.f90:135` (stale comment)
- Check: `src/atmosphere/meteo_orchestrator.f90` (must not reference the removed wrapper)

- [ ] **Step 1: Confirm `msw1eic` is unreachable**

Run:
```bash
cd /home/zawadzkim/Code/swap
grep -rn "msw1eic\|apply_interception_step" src/atmosphere/meteo_orchestrator.f90
grep -rn "swinter == 2\|swinter == 3\|swinter.eq.2\|swinter.eq.3" src/config/*.f90
```
Expected: the config validators reject `swinter=2/3` (Gash) on the TOML path, so the `msw1eic` Gash branch is dead. `meteo_orchestrator` calls `DivIntercep`/`VonHHBraden`/`Gash` — confirm whether `Gash` (the public name) routes into `msw1eic`. If `Gash` is the live Von-Hoyningen path and `msw1eic` is only the dead Sparse-Gash kernel, only `msw1eic` is removed.

- [ ] **Step 2: Remove `msw1eic` and its wrapper from `interception.f90`**

Edit `src/atmosphere/interception.f90`:
- Remove `msw1eic` from the module `public ::` list (line ~19).
- Delete the wrapper subroutine that calls `msw1eic` (the `apply_interception_step`-adjacent block, ~lines 135–229) **only if** it is the dead Gash path; keep `DivIntercep`, `VonHHBraden`, `Gash`, `ruttervw`.
- Delete the `msw1eic` subroutine itself (~lines 231–380, the `MSW1EIC.FOR` port with `real(4)`, `write+stop`, dead `!$OMP`).

If Step 1 shows the wrapper is reachable via a live `swinter` value, STOP and report — do not remove a live physics path.

- [ ] **Step 3: Delete the dormant shared-simulation file + stale comment**

Run:
```bash
git rm src/utils/dormant/sharedsimulation.f90
```
Edit `src/core/swap_log.f90:135` — remove the clause referencing `msw1eic`'s quarantined error path (reword to drop the MetaSWAP mention).

- [ ] **Step 4: Verify build + regression byte-identical**

Run: `pixi run clean && pixi run build-linux && pixi run -e test check-fast`
Expected: builds clean; pFUnit `OK (N tests)`; 4 regression cases PASS (byte-identical — `msw1eic` never ran).

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "chore(metaswap): drop dormant shared-simulation + dead msw1eic Gash kernel

sharedsimulation.f90 (obsolete file-based MetaSWAP coupling) deleted;
msw1eic Sparse-Gash interception (swinter=2/3, already rejected by TOML
validators -> unreachable) removed, taking the lone real(4) island,
write(*)+stop host-killer and dead !\$OMP block with it. We replace
MetaSWAP with direct MODFLOW coupling. Byte-identical.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 0.2: Delete other dead code

**Files:**
- Delete: `src/boundary/boundary_constants.f90` (empty, not in build)
- Modify: `src/boundary/README.md` (drop the reference to it)
- Modify: `src/config/soil_config.f90` (remove deprecated `z_init`)
- Modify: `src/crop/tillage.f90` (remove dead debug writes to hardcoded units)
- Modify: `src/soil/waterbalance.f90` (remove empty `if` block)

- [ ] **Step 1: Confirm `z_init` is unconsumed**

Run: `grep -rn "z_init" src/`
Expected: only the declaration + validation in `soil_config.f90` (commented DEPRECATED, "no longer populated or consumed"). If consumed anywhere, STOP.

- [ ] **Step 2: Remove the dead items**

- `git rm src/boundary/boundary_constants.f90`; edit `src/boundary/README.md` to delete the line pointing at it.
- In `src/config/soil_config.f90`: delete the `z_init(:)` field declaration and any `check_*`/validation referencing it.
- In `src/crop/tillage.f90`: delete the `write(222,…)`, `write(226,…)`, `write(333,…)`, `write(444,…)` debug statements (hardcoded units; the `TEST`-gated ones are already dead) and the dead `if(TEST)`/`if(TEST2)` blocks.
- In `src/soil/waterbalance.f90`: delete the empty `if (time%flZeroIntr) then / endif` block (no body).

- [ ] **Step 3: Rebuild + regression (z_init touches state-adjacent config → clean rebuild)**

Run: `pixi run clean && pixi run build-linux && pixi run -e test check-fast`
Expected: clean build; `OK (N tests)`; 4 cases PASS byte-identical.

- [ ] **Step 4: Commit**

```bash
git add -A
git commit -m "chore(cleanup): delete dead code (boundary_constants, z_init, tillage debug, empty if)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

# PHASE 1 — Returnable errors (prerequisite: a library must never abort the host)

### Task 1.1: Library-mode abort policy in `error_mod`

**Files:**
- Modify: `src/error/error.f90` (add a module-level abort-policy flag + status; gate the `error stop`)
- Test: `tests/unit/error/test_error_library_mode.pf` (create)
- Modify: `tests/unit/testSuites.inc` (register)

- [ ] **Step 1: Write the failing pFUnit test**

Create `tests/unit/error/test_error_library_mode.pf`:
```fortran
@test
subroutine test_library_mode_sets_status_not_abort()
   use error_mod
   use funit
   implicit none
   call global_errors%clear()
   call set_library_mode(.true.)
   call fatalerr_collected('TESTROUT', 'simulated fatal')
   @assertTrue(library_fatal_raised())
   call set_library_mode(.false.)   ! restore for other tests
   call global_errors%clear()
end subroutine

@test
subroutine test_default_mode_flag_clear()
   use error_mod
   use funit
   implicit none
   call set_library_mode(.false.)
   @assertFalse(library_fatal_raised())
end subroutine
```

- [ ] **Step 2: Register the suite**

Add to `tests/unit/testSuites.inc`:
```
ADD_TEST_SUITE(test_error_library_mode_suite)
```
(Match the existing macro style in that file; confirm the file's exact macro name by reading a neighbouring entry.)

- [ ] **Step 3: Run, expect FAIL**

Run: `pixi run -e test test-pfunit`
Expected: compile error / FAIL — `set_library_mode`, `library_fatal_raised` undefined.

- [ ] **Step 4: Implement the policy in `error.f90`**

In `src/error/error.f90`, add after the `global_errors` declaration (~line 57):
```fortran
   !> Library-embedding mode: when .true., a fatal does NOT `error stop`
   !! (which would kill a Python/MODFLOW host) but sets a sticky status
   !! flag the XMI layer can poll. Default .false. preserves CLI behaviour.
   logical, save, private :: library_mode    = .false.
   logical, save, private :: fatal_was_raised = .false.

   public :: set_library_mode, library_fatal_raised, clear_library_fatal
```
Add the procedures (in the module `contains`):
```fortran
   subroutine set_library_mode(on)
      logical, intent(in) :: on
      library_mode = on
      if (on) fatal_was_raised = .false.
   end subroutine set_library_mode

   logical function library_fatal_raised()
      library_fatal_raised = fatal_was_raised
   end function library_fatal_raised

   subroutine clear_library_fatal()
      fatal_was_raised = .false.
   end subroutine clear_library_fatal
```
Modify `error_collection_abort_if_fatal` (lines 187–192) to honour the flag:
```fortran
   subroutine error_collection_abort_if_fatal(self)
      class(error_collection_t), intent(in) :: self
      if (.not. self%has_fatals()) return
      write(error_unit, '(A)') self%summary()
      if (library_mode) then
         fatal_was_raised = .true.
         return                     ! caller (XMI) polls library_fatal_raised()
      end if
      error stop "fatal error(s) in swap input pipeline"
   end subroutine error_collection_abort_if_fatal
```

- [ ] **Step 5: Run, expect PASS**

Run: `pixi run -e test test-pfunit`
Expected: `OK (N tests)` with the count risen by 2.

- [ ] **Step 6: Regression unchanged (CLI never sets library_mode)**

Run: `pixi run -e test check-fast`
Expected: 4 cases PASS byte-identical.

- [ ] **Step 7: Commit**

```bash
git add src/error/error.f90 tests/unit/error/test_error_library_mode.pf tests/unit/testSuites.inc
git commit -m "feat(error): add library-mode abort policy (returnable status, no error stop)

Embedded/coupled runs must not be killed by the library. set_library_mode(.true.)
makes a fatal set a sticky poll-able flag instead of error stop; CLI default
unchanged (byte-identical).

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 1.2: Make the `dtutil` ADDSTR `error stop` returnable

**Files:**
- Modify: `src/core/dtutil.f90:365-369`

- [ ] **Step 1: Replace the raw `error stop` with the collected path**

In `src/core/dtutil.f90`, change the buffer-overflow guard so it routes through `fatalerr_collected` (which honours library mode) instead of a hard `error stop`:
```fortran
      IF (SIGLEN + L > LEN(STRING)) THEN
         call fatalerr_collected('ADDSTR', 'string buffer overflow')
         RETURN
      END IF
```
Add `use error_mod, only: fatalerr_collected` to the relevant procedure/module scope if not already imported (check the top of `dtutil.f90`).

- [ ] **Step 2: Rebuild + regression**

Run: `pixi run build-linux && pixi run -e test check-fast`
Expected: 4 cases PASS byte-identical (overflow path is not exercised by the cases).

- [ ] **Step 3: Commit**

```bash
git add src/core/dtutil.f90
git commit -m "feat(error): route ADDSTR overflow through fatalerr_collected (no raw error stop)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

# PHASE 2 — Bottom-boundary GWL injection seam

### Task 2.1: Coupled-GWL inject hook on soilwater state

**Files:**
- Modify: `src/state/soilwater_state.f90` (add a `flcoupled_gwl` flag + `gwl_injected` field)
- Modify: `src/boundary/boundbottom.f90:56-59` (use the injected value when coupled)
- Test: `tests/unit/boundary/test_coupled_gwl_inject.pf` (create) + register

**This task changes `src/state/*` schema → mandatory clean rebuild.**

- [ ] **Step 1: Write the failing pFUnit test**

Create `tests/unit/boundary/test_coupled_gwl_inject.pf`:
```fortran
@test
subroutine test_injected_gwl_overrides_table()
   use swap_state_mod, only: swap_state_t
   use funit
   implicit none
   type(swap_state_t) :: state
   ! Simulate coupled injection: flag on, value set
   state%soilwater%swbotb_runtime = 1
   state%soilwater%flcoupled_gwl  = .true.
   state%soilwater%gwl_injected   = -120.0d0   ! cm
   call apply_coupled_gwl(state)               ! helper under test
   @assertEqual(-120.0d0, state%soilwater%gwlinp, tolerance=1.0d-9)
end subroutine
```

- [ ] **Step 2: Register** in `tests/unit/testSuites.inc`:
```
ADD_TEST_SUITE(test_coupled_gwl_inject_suite)
```

- [ ] **Step 3: Run, expect FAIL**

Run: `pixi run -e test test-pfunit`
Expected: FAIL — `flcoupled_gwl`, `gwl_injected`, `apply_coupled_gwl` undefined.

- [ ] **Step 4: Add the state fields**

In `src/state/soilwater_state.f90`, after the `gwlinp` declaration (line ~80):
```fortran
      logical      :: flcoupled_gwl = .false.  !! when .true., gwlinp is injected externally (MODFLOW coupling)
      real(real64) :: gwl_injected  = 0.0_real64 !! externally injected groundwater level, swbotb=1 coupled (cm)
```

- [ ] **Step 5: Add `apply_coupled_gwl` and wire it into BoundBottom**

In `src/boundary/boundbottom.f90`, add a small public helper (same module) and call it in the swbotb=1 branch. Replace lines 56–59:
```fortran
         ! ---- swbotb=1: groundwater level (table or injected) ---------------
         if (soil%swbotb_runtime .eq. 1) then
            if (soil%flcoupled_gwl) then
               call apply_coupled_gwl(state)
            else
               soil%gwlinp = afgen(soil%gwltab, mabbc*2, time%t1900 + time%dt)
            end if
         end if
```
Add (module `contains`):
```fortran
   subroutine apply_coupled_gwl(state)
      type(swap_state_t), intent(inout) :: state
      associate (soil => state%soilwater)
         soil%gwlinp = soil%gwl_injected
      end associate
   end subroutine apply_coupled_gwl
```
Export it via the module `public` list. (`state`/`soil`/`time` aliasing follows the existing pattern in `BoundBottom`.)

- [ ] **Step 6: Clean rebuild + run, expect PASS**

Run: `pixi run clean && pixi run build-linux && pixi run -e test test-pfunit`
Expected: `OK (N tests)` count risen by 1.

- [ ] **Step 7: Regression byte-identical (flag defaults .false.)**

Run: `pixi run -e test check-fast`
Expected: 4 cases PASS byte-identical.

- [ ] **Step 8: Commit**

```bash
git add src/state/soilwater_state.f90 src/boundary/boundbottom.f90 tests/unit/boundary/test_coupled_gwl_inject.pf tests/unit/testSuites.inc
git commit -m "feat(boundary): coupled-GWL injection hook for swbotb=1

flcoupled_gwl + gwl_injected on soilwater state; BoundBottom uses the
injected value instead of the afgen table when coupled. Default off ->
byte-identical.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

# PHASE 3 — Ensemble core

### Task 3.1: `swap_ensemble_mod` — N columns over one shared config

**Files:**
- Create: `src/core/swap_ensemble_mod.f90`
- Modify: `meson.build` (add to `modern_sources`)
- Test: `tests/unit/core/test_swap_ensemble.pf` (create) + register

**Touches build + new module → clean rebuild.**

- [ ] **Step 1: Write the failing pFUnit test**

Create `tests/unit/core/test_swap_ensemble.pf`:
```fortran
@test
subroutine test_ensemble_init_and_step()
   use swap_ensemble_mod
   use funit
   implicit none
   integer :: rc
   character(len=*), parameter :: case_dir = &
      'tests/swap-cases/toml/1.hupselbrook'   ! note: must run from repo root
   ! Build a 3-column homogeneous ensemble from the coupled config.
   rc = ensemble_init(case_dir//'/swap_coupled.toml', ncol_in=3)
   @assertEqual(0, rc)
   @assertEqual(3, ensemble_ncol())
   ! All columns must be in prescribed-GWL coupled mode.
   @assertTrue(ensemble_all_swbotb1())
   ! Inject a head, step one day, read back a (finite) recharge depth.
   call ensemble_set_gwl(1, -0.80d0)   ! metres
   rc = ensemble_step_day()
   @assertEqual(0, rc)
   @assertTrue(ensemble_qbot_volume(1) == ensemble_qbot_volume(1))  ! not NaN
   rc = ensemble_finalize()
   @assertEqual(0, rc)
end subroutine
```
*(Note: this test reads a config that Task 7.2 creates. To keep Phase 3 self-contained, Step 1b below creates a minimal `swap_coupled.toml` now; Task 7.2 finalises it.)*

- [ ] **Step 1b: Create the coupled config (minimal, finalised in Task 7.2)**

Copy the hupselbrook TOML and apply the two deltas:
```bash
mkdir -p tests/swap-cases/toml/1.hupselbrook
cp tests/swap-cases/toml/1.hupselbrook/swap.toml tests/swap-cases/toml/1.hupselbrook/swap_coupled.toml
```
Edit `swap_coupled.toml`:
- `[bottom_boundary]` → `swbotb = 1`
- `[drainage]` → set `swdra = 0` (and remove the `file = "swap.dra.toml"` line; drainage off).

- [ ] **Step 2: Register** in `tests/unit/testSuites.inc`:
```
ADD_TEST_SUITE(test_swap_ensemble_suite)
```

- [ ] **Step 3: Run, expect FAIL** (module missing).

Run: `pixi run -e test test-pfunit`
Expected: FAIL — `swap_ensemble_mod` undefined.

- [ ] **Step 4: Implement `swap_ensemble_mod`**

Create `src/core/swap_ensemble_mod.f90`:
```fortran
!> SWAP ensemble: N independent column states sharing read-only config(s),
!! stepped sequentially. Backing store for the XMI coupling kernel. One
!! ensemble per process (matches imod_coupler's one-kernel-per-lib model).
!! Per-column config indirection (column_config -> configs) is built in;
!! the demo uses a single shared config (M=1).
module swap_ensemble_mod
   use iso_fortran_env, only: real64, int32
   use swap_state_mod,  only: swap_state_t
   use swap_config_mod, only: swap_config_t
   use swap_mod,        only: swap_init_from_loaded_config, swap_run_step, swap_close
   use load_swap_config_mod, only: load_swap_config
   use error_mod,       only: set_library_mode, library_fatal_raised, clear_library_fatal
   implicit none
   private

   public :: ensemble_init, ensemble_step_day, ensemble_finalize
   public :: ensemble_ncol, ensemble_all_swbotb1
   public :: ensemble_set_gwl, ensemble_qbot_volume, ensemble_storage_coef
   ! Exposed (target) exchange arrays for the XMI get_value_ptr layer:
   public :: gwl, qbot_volume, storage_coef

   integer, parameter :: STORAGE_COEF_DEFAULT_X100 = 15   ! sy = 0.15 (smoke placeholder)

   type(swap_state_t),  allocatable, save :: columns(:)
   type(swap_config_t), allocatable, save, target :: configs(:)
   integer,             allocatable, save :: column_config(:)
   integer,             save :: ncol = 0

   ! Exchange arrays (length ncol). target => stable address for c_loc.
   real(real64), allocatable, save, target :: gwl(:)          ! MODFLOW head, metres (driver writes)
   real(real64), allocatable, save, target :: qbot_volume(:)  ! recharge depth over step, metres (we write)
   real(real64), allocatable, save, target :: storage_coef(:) ! specific yield, - (we write)

contains

   integer function ensemble_init(config_file, ncol_in) result(rc)
      character(len=*), intent(in) :: config_file
      integer,          intent(in) :: ncol_in
      integer :: i
      rc = 0
      call set_library_mode(.true.)
      call clear_library_fatal()
      ncol = ncol_in
      allocate(configs(1))
      configs(1) = load_swap_config(config_file)   ! parse+validate+finalize one shared config
      if (library_fatal_raised()) then; rc = 1; return; end if
      allocate(column_config(ncol)); column_config = 1
      allocate(columns(ncol))
      do i = 1, ncol
         call swap_init_from_loaded_config(columns(i), configs(column_config(i)))
         if (library_fatal_raised()) then; rc = 2; return; end if
         ! Enable coupled GWL injection on every column.
         columns(i)%soilwater%flcoupled_gwl = .true.
      end do
      allocate(gwl(ncol),          source=0.0_real64)
      allocate(qbot_volume(ncol),  source=0.0_real64)
      allocate(storage_coef(ncol), source=real(STORAGE_COEF_DEFAULT_X100,real64)/100.0_real64)
      if (.not. ensemble_all_swbotb1()) rc = 3   ! coupling precondition
   end function ensemble_init

   integer function ensemble_step_day() result(rc)
      integer      :: i
      real(real64) :: dt_days
      rc = 0
      do i = 1, ncol
         ! Inject this cell's head (m -> cm, datum surface=0).
         columns(i)%soilwater%gwl_injected = gwl(i) * 100.0_real64
         call swap_run_step(columns(i), configs(column_config(i)))
         if (library_fatal_raised()) then; rc = 1; return; end if
         dt_days = columns(i)%timecontrol%dt
         ! Recharge depth over step (m): downward qbot (cm/d, negative) -> positive m.
         qbot_volume(i) = (-columns(i)%soilwater%qbot) * 0.01_real64 * dt_days
         ! storage_coef left at the init constant (smoke placeholder).
      end do
   end function ensemble_step_day

   integer function ensemble_finalize() result(rc)
      integer :: i
      rc = 0
      do i = 1, ncol
         call swap_close(columns(i), configs(column_config(i)))
      end do
      if (allocated(columns))       deallocate(columns)
      if (allocated(configs))       deallocate(configs)
      if (allocated(column_config)) deallocate(column_config)
      if (allocated(gwl))           deallocate(gwl)
      if (allocated(qbot_volume))   deallocate(qbot_volume)
      if (allocated(storage_coef))  deallocate(storage_coef)
      ncol = 0
   end function ensemble_finalize

   integer function ensemble_ncol(); ensemble_ncol = ncol; end function

   logical function ensemble_all_swbotb1()
      integer :: i
      ensemble_all_swbotb1 = (ncol > 0)
      do i = 1, ncol
         if (columns(i)%soilwater%swbotb_runtime /= 1) ensemble_all_swbotb1 = .false.
      end do
   end function

   subroutine ensemble_set_gwl(i, head_m)
      integer,      intent(in) :: i
      real(real64), intent(in) :: head_m
      gwl(i) = head_m
   end subroutine

   real(real64) function ensemble_qbot_volume(i); integer, intent(in) :: i
      ensemble_qbot_volume = qbot_volume(i); end function
   real(real64) function ensemble_storage_coef(i); integer, intent(in) :: i
      ensemble_storage_coef = storage_coef(i); end function

end module swap_ensemble_mod
```
*Verify the exact name of the config loader* (`load_swap_config` vs `load_swap_config` returning a `swap_config_t` vs a subroutine). Read `src/io/toml/` to confirm whether it's a function returning `swap_config_t` or a subroutine with `intent(out)`. If a subroutine, adapt the `configs(1) = ...` line to `call load_swap_config(config_file, configs(1), errors)` and poll `errors`/`library_fatal_raised()`.

- [ ] **Step 5: Add to build**

In `meson.build`, add `'src/core/swap_ensemble_mod.f90'` to the `modern_sources` list (lines ~60-65), placed after `swap_mod.f90` (dependency order: it `use`s `swap_mod`).

- [ ] **Step 6: Clean rebuild + run, expect PASS**

Run: `pixi run clean && pixi run build-linux && pixi run -e test test-pfunit`
Expected: `OK (N tests)` risen by 1.

- [ ] **Step 7: Regression unchanged**

Run: `pixi run -e test check-fast`
Expected: 4 cases PASS byte-identical.

- [ ] **Step 8: Commit**

```bash
git add src/core/swap_ensemble_mod.f90 meson.build tests/unit/core/test_swap_ensemble.pf tests/unit/testSuites.inc tests/swap-cases/toml/1.hupselbrook/swap_coupled.toml
git commit -m "feat(core): swap_ensemble_mod — N homogeneous columns, sequential step, exchange arrays

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

# PHASE 4 — XMI shared-library facade

### Task 4.1: `BMI_LEN*` exported integer constants

**Files:**
- Create: `src/core/bmi_constants_mod.f90`
- Modify: `meson.build` (`modern_sources`)

- [ ] **Step 1: Implement the constants module**

xmipy reads these as C global ints via `c_int.in_dll`. They must be `bind(C)` module variables (not `parameter`, so they get external symbols). Create `src/core/bmi_constants_mod.f90`:
```fortran
!> BMI string-length constants exported as C globals (xmipy reads them via
!! c_int.in_dll). Must be bind(C) module variables to have external linkage.
module bmi_constants_mod
   use iso_c_binding, only: c_int
   implicit none
   integer(c_int), bind(C, name='BMI_LENCOMPONENTNAME') :: BMI_LENCOMPONENTNAME = 256
   integer(c_int), bind(C, name='BMI_LENVERSION')       :: BMI_LENVERSION       = 256
   integer(c_int), bind(C, name='BMI_LENVARADDRESS')    :: BMI_LENVARADDRESS    = 256
   integer(c_int), bind(C, name='BMI_LENVARTYPE')       :: BMI_LENVARTYPE       = 256
   integer(c_int), bind(C, name='BMI_LENGRIDTYPE')      :: BMI_LENGRIDTYPE       = 256
   integer(c_int), bind(C, name='BMI_LENERRMESSAGE')    :: BMI_LENERRMESSAGE    = 1024
end module bmi_constants_mod
```

- [ ] **Step 2: Add to `modern_sources` in `meson.build`** (before `swap_xmi_mod`).

- [ ] **Step 3: Build**

Run: `pixi run build-linux`
Expected: compiles.

- [ ] **Step 4: Commit**

```bash
git add src/core/bmi_constants_mod.f90 meson.build
git commit -m "feat(xmi): export BMI_LEN* constants as C globals for xmipy

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 4.2: XMI lifecycle + the variable-pointer protocol

**Files:**
- Create: `src/core/swap_xmi_mod.f90`
- Modify: `meson.build` (`modern_sources`)

All functions are `function … result(rc) bind(C, name='…')` returning `c_int` (0 = success), following the existing `swap_bmi_mod.f90` pattern. xmipy passes most outputs by reference and `prepare_time_step`'s dt and the solve ids by reference (`int*`/`double*`); `update_until` is by value.

- [ ] **Step 1: Implement the XMI module**

Create `src/core/swap_xmi_mod.f90`:
```fortran
!> XMI (BMI + Deltares eXtended Model Interface) C-ABI facade over the SWAP
!! ensemble. Symbol names/signatures match what xmipy's XmiWrapper calls.
module swap_xmi_mod
   use iso_c_binding
   use swap_ensemble_mod
   use error_mod, only: library_fatal_raised
   implicit none
   private

   ! Last-error text buffer (filled when a function returns rc /= 0).
   character(len=1024), save :: last_error = ''

contains

   ! ---- lifecycle -------------------------------------------------------
   integer(c_int) function xmi_initialize(config_file) result(rc) bind(C, name='initialize')
      character(kind=c_char), intent(in) :: config_file(*)
      character(len=512) :: f_path
      integer :: ncol
      call c_to_f_string(config_file, f_path)
      ! f_path points at the SWAP working dir's coupled config; ncol comes
      ! from the sidecar ensemble descriptor 'ensemble.txt' (one int).
      ncol = read_ncol_sidecar(trim(f_path))
      if (ncol <= 0) then; last_error = 'ensemble.txt missing/invalid'; rc = 1; return; end if
      rc = ensemble_init(trim(f_path), ncol)
      if (rc /= 0) last_error = 'ensemble_init failed'
   end function xmi_initialize

   integer(c_int) function xmi_update() result(rc) bind(C, name='update')
      rc = ensemble_step_day()
   end function

   integer(c_int) function xmi_update_until(t) result(rc) bind(C, name='update_until')
      real(c_double), value, intent(in) :: t
      rc = ensemble_step_day()   ! MF6 leads; one day per call in coupled mode
   end function

   integer(c_int) function xmi_finalize() result(rc) bind(C, name='finalize')
      rc = ensemble_finalize()
   end function

   ! ---- time ------------------------------------------------------------
   integer(c_int) function xmi_get_current_time(t) result(rc) bind(C, name='get_current_time')
      real(c_double), intent(out) :: t
      t = ensemble_current_t1900(); rc = 0
   end function
   integer(c_int) function xmi_get_start_time(t) result(rc) bind(C, name='get_start_time')
      real(c_double), intent(out) :: t
      t = ensemble_start_t1900(); rc = 0
   end function
   integer(c_int) function xmi_get_end_time(t) result(rc) bind(C, name='get_end_time')
      real(c_double), intent(out) :: t
      t = ensemble_end_t1900(); rc = 0
   end function
   integer(c_int) function xmi_get_time_step(dt) result(rc) bind(C, name='get_time_step')
      real(c_double), intent(out) :: dt
      dt = 1.0_c_double; rc = 0     ! daily coupling
   end function

   ! ---- component info --------------------------------------------------
   integer(c_int) function xmi_get_version(buf) result(rc) bind(C, name='get_version')
      character(kind=c_char), intent(out) :: buf(*)
      call f_to_c_string('SWAP-modern', buf); rc = 0
   end function
   integer(c_int) function xmi_get_component_name(buf) result(rc) bind(C, name='get_component_name')
      character(kind=c_char), intent(out) :: buf(*)
      call f_to_c_string('SWAP', buf); rc = 0
   end function
   integer(c_int) function xmi_get_input_item_count(n) result(rc) bind(C, name='get_input_item_count')
      integer(c_int), intent(out) :: n
      n = 1; rc = 0                 ! gwl
   end function
   integer(c_int) function xmi_get_output_item_count(n) result(rc) bind(C, name='get_output_item_count')
      integer(c_int), intent(out) :: n
      n = 2; rc = 0                 ! qbot_volume, storage_coef
   end function

   ! ---- XMI time-step / solve control ----------------------------------
   integer(c_int) function xmi_prepare_time_step(dt) result(rc) bind(C, name='prepare_time_step')
      real(c_double), intent(in) :: dt        ! by reference (xmipy passes byref)
      rc = 0
   end function
   integer(c_int) function xmi_finalize_time_step() result(rc) bind(C, name='finalize_time_step')
      rc = 0
   end function
   integer(c_int) function xmi_prepare_solve(cid) result(rc) bind(C, name='prepare_solve')
      integer(c_int), intent(in) :: cid
      rc = 0
   end function
   integer(c_int) function xmi_solve(cid, has_converged) result(rc) bind(C, name='solve')
      integer(c_int), intent(in)  :: cid
      integer(c_int), intent(out) :: has_converged
      rc = ensemble_step_day()      ! the actual column loop
      has_converged = 1
   end function
   integer(c_int) function xmi_finalize_solve(cid) result(rc) bind(C, name='finalize_solve')
      integer(c_int), intent(in) :: cid
      rc = 0
   end function
   integer(c_int) function xmi_get_subcomponent_count(n) result(rc) bind(C, name='get_subcomponent_count')
      integer(c_int), intent(out) :: n
      n = 1; rc = 0
   end function

   ! ---- variable metadata ----------------------------------------------
   integer(c_int) function xmi_get_var_rank(name, r) result(rc) bind(C, name='get_var_rank')
      character(kind=c_char), intent(in)  :: name(*)
      integer(c_int),         intent(out) :: r
      r = 1; rc = 0                 ! all three exchange vars are rank-1
   end function
   integer(c_int) function xmi_get_var_type(name, tbuf) result(rc) bind(C, name='get_var_type')
      character(kind=c_char), intent(in)  :: name(*)
      character(kind=c_char), intent(out) :: tbuf(*)
      call f_to_c_string('double', tbuf); rc = 0   ! string must start 'double'
   end function
   integer(c_int) function xmi_get_var_shape(name, shp) result(rc) bind(C, name='get_var_shape')
      character(kind=c_char), intent(in)  :: name(*)
      integer(c_int),         intent(out) :: shp(*)
      shp(1) = ensemble_ncol(); rc = 0
   end function
   integer(c_int) function xmi_get_var_itemsize(name, sz) result(rc) bind(C, name='get_var_itemsize')
      character(kind=c_char), intent(in)  :: name(*)
      integer(c_int),         intent(out) :: sz
      sz = 8; rc = 0
   end function
   integer(c_int) function xmi_get_var_nbytes(name, nb) result(rc) bind(C, name='get_var_nbytes')
      character(kind=c_char), intent(in)  :: name(*)
      integer(c_int),         intent(out) :: nb
      nb = 8 * ensemble_ncol(); rc = 0
   end function

   ! ---- the zero-copy pointer protocol ---------------------------------
   integer(c_int) function xmi_get_value_ptr(name, ptr) result(rc) bind(C, name='get_value_ptr')
      character(kind=c_char), intent(in)  :: name(*)
      type(c_ptr),            intent(out) :: ptr
      character(len=64) :: nm
      call c_to_f_string(name, nm)
      rc = 0
      select case (trim(nm))
      case ('gwl');          ptr = c_loc(gwl(1))
      case ('qbot_volume');  ptr = c_loc(qbot_volume(1))
      case ('storage_coef'); ptr = c_loc(storage_coef(1))
      case default;          ptr = c_null_ptr; rc = 1; last_error = 'unknown var: '//trim(nm)
      end select
   end function

   ! ---- error reporting -------------------------------------------------
   integer(c_int) function xmi_get_last_bmi_error(buf) result(rc) bind(C, name='get_last_bmi_error')
      character(kind=c_char), intent(out) :: buf(*)
      call f_to_c_string(trim(last_error), buf); rc = 0
   end function

   ! ---- private C-string helpers (copy from swap_bmi_mod.f90) -----------
   subroutine c_to_f_string(c_str, f_str)
      character(kind=c_char), intent(in)  :: c_str(*)
      character(len=*),       intent(out) :: f_str
      integer :: i
      f_str = ''
      do i = 1, len(f_str)
         if (c_str(i) == c_null_char) exit
         f_str(i:i) = c_str(i)
      end do
   end subroutine
   subroutine f_to_c_string(f_str, c_str)
      character(len=*),       intent(in)  :: f_str
      character(kind=c_char), intent(out) :: c_str(*)
      integer :: i
      do i = 1, len_trim(f_str)
         c_str(i) = f_str(i:i)
      end do
      c_str(len_trim(f_str)+1) = c_null_char
   end subroutine

   integer function read_ncol_sidecar(dir_or_file) result(n)
      character(len=*), intent(in) :: dir_or_file
      integer :: u, ios
      character(len=512) :: path
      n = 0
      ! ensemble.txt sits next to the config; accept a dir or a file path.
      path = trim(dir_or_file)
      if (index(path, '.toml') > 0) path = path(1:scan(path,'/',back=.true.))
      open(newunit=u, file=trim(path)//'ensemble.txt', status='old', action='read', iostat=ios)
      if (ios /= 0) return
      read(u, *, iostat=ios) n
      close(u)
      if (ios /= 0) n = 0
   end function

end module swap_xmi_mod
```
This module references ensemble time accessors (`ensemble_current_t1900`, `ensemble_start_t1900`, `ensemble_end_t1900`) — add them to `swap_ensemble_mod` in Step 2.

- [ ] **Step 2: Add ensemble time accessors**

In `src/core/swap_ensemble_mod.f90`, add to the `public` list and `contains`:
```fortran
   public :: ensemble_current_t1900, ensemble_start_t1900, ensemble_end_t1900
   ...
   real(real64) function ensemble_current_t1900()
      ensemble_current_t1900 = columns(1)%timecontrol%t1900
   end function
   real(real64) function ensemble_start_t1900()
      ensemble_start_t1900 = columns(1)%timecontrol%tstart
   end function
   real(real64) function ensemble_end_t1900()
      ensemble_end_t1900 = columns(1)%timecontrol%tend
   end function
```
*Verify the exact field names* `tstart`/`tend`/`t1900` on `timecontrol_state_t` (the CAPI `swap_get_scalar` reads `tstart`/`tend`/`t1900` — reuse those exact names).

- [ ] **Step 3: Add `swap_xmi_mod.f90` to `modern_sources`** (after `swap_ensemble_mod.f90`).

- [ ] **Step 4: Build**

Run: `pixi run clean && pixi run build-linux`
Expected: compiles; `libswap_bmi.so` now also exports the XMI symbols.

- [ ] **Step 5: Verify the symbols are exported**

Run:
```bash
nm -D builddir/libswap_bmi.so | grep -E " T (initialize|update|solve|prepare_solve|prepare_time_step|finalize_solve|get_value_ptr|get_var_rank|get_var_type|get_var_shape)$"
nm -D builddir/libswap_bmi.so | grep -E "BMI_LEN"
```
Expected: all listed symbols present as `T` (text/defined); the `BMI_LEN*` as `B`/`D` (data).

- [ ] **Step 6: Regression unchanged**

Run: `pixi run -e test check-fast`
Expected: 4 cases PASS byte-identical.

- [ ] **Step 7: Commit**

```bash
git add src/core/swap_xmi_mod.f90 src/core/swap_ensemble_mod.f90 meson.build
git commit -m "feat(xmi): XMI C-ABI facade over the ensemble (lifecycle, solve control, get_value_ptr)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

# PHASE 5 — SWAP-only XMI smoke (no MODFLOW yet)

### Task 5.1: Drive `libswap.so` via xmipy and assert the contract

**Files:**
- Create: `tests/coupling/test_swap_xmi_smoke.py`
- Create: `tests/coupling/ensemble.txt` (contains `3`)
- Modify: `tests/coupling/meson.build` (create) + register `subdir('tests/coupling')` in root `meson.build`

- [ ] **Step 1: Ensure xmipy is available**

Confirm `xmipy` is in the pixi test env (`pixi list -e test | grep xmipy`). If absent, add to `pixi.toml` test feature deps: `xmipy = "*"` (conda-forge) and `pixi install`.

- [ ] **Step 2: Write the smoke test**

Create `tests/coupling/ensemble.txt` with a single line `3`.
Create `tests/coupling/test_swap_xmi_smoke.py`:
```python
"""SWAP-only XMI contract smoke test (no MODFLOW). Loads libswap via xmipy,
runs a few coupled-style steps, exercises get_value_ptr on the 3 exchange
arrays, asserts shapes and finite values."""
import sys, shutil
from pathlib import Path
import numpy as np
from xmipy import XmiWrapper

lib_path = Path(sys.argv[1])           # libswap_bmi.so
case_dir = Path(sys.argv[2])           # .../1.hupselbrook
work = Path(sys.argv[3])               # writable work dir (meson build dir)

# Stage a work dir: coupled config + meteo + ensemble.txt
for f in ["swap_coupled.toml", "283.csv", "grassd.crp.toml", "maizes.crp.toml", "potatod.crp.toml"]:
    shutil.copy(case_dir / f, work / f)
(work / "ensemble.txt").write_text("3\n")

swap = XmiWrapper(lib_path=str(lib_path), working_directory=str(work))
swap.initialize(str(work / "swap_coupled.toml"))

gwl          = swap.get_value_ptr("gwl")
qbot_volume  = swap.get_value_ptr("qbot_volume")
storage_coef = swap.get_value_ptr("storage_coef")
assert gwl.shape == (3,), gwl.shape
assert storage_coef.shape == (3,)
assert np.allclose(storage_coef, 0.15)

# Three coupled-style daily steps with distinct injected heads per column.
for day in range(3):
    gwl[:] = np.array([-0.5, -0.75, -1.0])      # metres
    swap.prepare_time_step(1.0)
    swap.prepare_solve(0)
    converged = swap.solve(0)
    swap.finalize_solve(0)
    swap.finalize_time_step()
    assert np.all(np.isfinite(qbot_volume)), qbot_volume

swap.finalize()
print("SWAP XMI smoke OK:", qbot_volume)
```

- [ ] **Step 3: Register the meson test**

Create `tests/coupling/meson.build`:
```python
python3 = import('python').find_installation('python3')

test('swap-xmi-smoke',
    python3,
    args: [
        meson.current_source_dir() / 'test_swap_xmi_smoke.py',
        swap_bmi_lib.full_path(),
        meson.project_source_root() / 'tests/swap-cases/toml/1.hupselbrook',
        meson.current_build_dir(),
    ],
    suite: 'coupling',
    depends: [swap_bmi_lib],
    timeout: 300,
)
```
Add to root `meson.build` (near the other `subdir('tests/...')` lines): `subdir('tests/coupling')`.

- [ ] **Step 4: Run, expect PASS**

Run: `pixi run build-linux && meson test -C builddir swap-xmi-smoke -v`
Expected: prints `SWAP XMI smoke OK: [...]` with three finite recharge depths. If `solve` returns non-zero, call `get_last_bmi_error` — fix the ensemble/config wiring (most likely the `load_swap_config` call shape from Task 3.1 Step 4, or `ncol` sidecar path).

- [ ] **Step 5: Commit**

```bash
git add tests/coupling/test_swap_xmi_smoke.py tests/coupling/ensemble.txt tests/coupling/meson.build meson.build
git commit -m "test(coupling): SWAP-only XMI smoke via xmipy (3-column ensemble)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

# PHASE 6 — flopy MODFLOW 6 model

### Task 6.1: Two-channel cross-section builder

**Files:**
- Create: `tests/coupling/build_modflow.py`
- Create: `tests/coupling/test_modflow_standalone.py`

- [ ] **Step 1: Ensure flopy + a MODFLOW 6 binary are available**

Confirm `flopy` in the pixi test env; ensure a `mf6`/`libmf6.so` is available (conda-forge `modflow6` or `flopy`'s `get-modflow`). Record the `libmf6.so` path for the coupler config (Task 7.2).

- [ ] **Step 2: Write the model builder**

Create `tests/coupling/build_modflow.py`:
```python
"""Build a 1-layer x 1-row x N-col unconfined transient MODFLOW 6 model:
two CHD 'channels' at the ends, RCHA on the interior cells (coupling target),
K grounded to SWAP's ~12.5 cm/d soil (=0.125 m/d). Lengths in metres, time in days."""
from pathlib import Path
import flopy

def build(ws: Path, ncol: int = 10, nper: int = 1096, h0: float = -0.5, h1: float = -1.0):
    ws.mkdir(parents=True, exist_ok=True)
    sim = flopy.mf6.MFSimulation(sim_name="swapmf", sim_ws=str(ws), exe_name="mf6")
    # Daily transient stress periods (perlen=1 d, single step each).
    tdis = flopy.mf6.ModflowTdis(sim, time_units="days",
                                 perioddata=[(1.0, 1, 1.0)] * nper, nper=nper)
    ims = flopy.mf6.ModflowIms(sim, complexity="SIMPLE", outer_maximum=50,
                               inner_maximum=100, linear_acceleration="BICGSTAB")
    gwf = flopy.mf6.ModflowGwf(sim, modelname="swapmf", newtonoptions="NEWTON",
                               save_flows=True)
    delr = 10.0   # cell width [m]
    dis = flopy.mf6.ModflowGwfdis(gwf, nlay=1, nrow=1, ncol=ncol,
                                  delr=delr, delc=10.0, top=0.0, botm=-10.0,
                                  length_units="meters")
    flopy.mf6.ModflowGwfic(gwf, strt=-0.75)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=0.125)        # unconfined, K=0.125 m/d
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=1e-5, sy=0.15, transient={0: True})
    # Two channels: CHD on first and last column.
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[
        [(0, 0, 0), h0], [(0, 0, ncol - 1), h1]])
    # Recharge on interior cells — the coupling package (named 'RCHA').
    rch_cells = {0: [[(0, 0, j), 0.0] for j in range(1, ncol - 1)]}
    rch = flopy.mf6.ModflowGwfrch(gwf, stress_period_data=rch_cells, pname="RCHA",
                                  maxbound=ncol - 2)
    oc = flopy.mf6.ModflowGwfoc(gwf, head_filerecord="swapmf.hds",
                                budget_filerecord="swapmf.cbc",
                                saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")])
    sim.write_simulation()
    return sim

if __name__ == "__main__":
    import sys
    build(Path(sys.argv[1]), ncol=int(sys.argv[2]) if len(sys.argv) > 2 else 10,
          nper=int(sys.argv[3]) if len(sys.argv) > 3 else 1096)
    print("wrote MODFLOW 6 model")
```
*Note:* the number of RCHA interior cells (`ncol-2`) is the coupling count N that SWAP's `ensemble.txt` must match. For `ncol=10`, N=8.

- [ ] **Step 3: Write the standalone-run test**

Create `tests/coupling/test_modflow_standalone.py`:
```python
import sys, tempfile
from pathlib import Path
import flopy
from build_modflow import build

ws = Path(tempfile.mkdtemp()) / "mf"
sim = build(ws, ncol=10, nper=3)
ret, _ = sim.run_simulation(silent=True)
assert ret, "MODFLOW standalone run failed"
hds = flopy.utils.HeadFile(ws / "swapmf.hds").get_data().squeeze()
assert hds[0] != hds[-1], "expected a lateral gradient between the two channels"
print("MODFLOW standalone OK, heads:", hds)
```

- [ ] **Step 4: Run, expect PASS**

Run: `cd tests/coupling && pixi run -e test python test_modflow_standalone.py`
Expected: prints heads with a gradient between the ends.

- [ ] **Step 5: Commit**

```bash
git add tests/coupling/build_modflow.py tests/coupling/test_modflow_standalone.py
git commit -m "feat(coupling): flopy two-channel MODFLOW 6 builder + standalone test

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

# PHASE 7 — Coupler fork + configs

### Task 7.1: Vendor the swap-branch driver + wrapper (native names)

**Files:**
- Create: `tests/coupling/imod_coupler_fork/swap_wrapper.py`
- Create: `tests/coupling/imod_coupler_fork/swapmod.py`
- Create: `tests/coupling/imod_coupler_fork/__init__.py`

- [ ] **Step 1: Vendor `swap_wrapper.py` with native names**

Create `tests/coupling/imod_coupler_fork/swap_wrapper.py`:
```python
"""Fork of imod_coupler swap-branch SwapWrapper using SWAP-native var names."""
from pathlib import Path
from typing import Union
import numpy as np
from numpy.typing import NDArray
from xmipy import XmiWrapper

class SwapWrapper(XmiWrapper):
    def __init__(self, lib_path, lib_dependency=None, working_directory=None, timing=False):
        super().__init__(lib_path, lib_dependency, working_directory, timing)

    def get_head_ptr(self) -> NDArray[np.float64]:
        return self.get_value_ptr("gwl")            # was 'dhgwmod'

    def get_volume_ptr(self) -> NDArray[np.float64]:
        return self.get_value_ptr("qbot_volume")    # was 'dvsim'

    def get_storage_ptr(self) -> NDArray[np.float64]:
        return self.get_value_ptr("storage_coef")   # was 'dsc1sim'
```

- [ ] **Step 2: Vendor `swapmod.py`**

Obtain the swap-branch `swapmod.py` and adapt: import the local `swap_wrapper`, and resolve the loose 1:1 sequential exchange (the commented `create_mapping`/in-loop re-solve are not needed — keep `update()` as: heads→SWAP, prepare/solve/finalize SWAP once, storage+recharge→MF6, MF6 convergence loop). Save to `tests/coupling/imod_coupler_fork/swapmod.py`. Add `__init__.py` exporting `SwapMod`.

*(The driver body is reproduced in the spec §2; copy it verbatim, changing the `SwapWrapper` import to the local fork and confirming `exchange_swap2mod`/`exchange_mod2swap` use `swap_head`/`swap_volume`/`swap_storage`.)*

- [ ] **Step 3: Commit**

```bash
git add tests/coupling/imod_coupler_fork/
git commit -m "feat(coupling): vendor swapmod driver + wrapper fork (SWAP-native exchange names)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 7.2: Coupled SWAP config + coupler config

**Files:**
- Finalise: `tests/swap-cases/toml/1.hupselbrook/swap_coupled.toml` (from Task 3.1)
- Create: `tests/coupling/swap_config/` (staged copy with `ensemble.txt`)
- Create: `tests/coupling/imod_coupler.toml`

- [ ] **Step 1: Finalise the coupled SWAP config**

Confirm `swap_coupled.toml` has `[bottom_boundary] swbotb = 1`, `[drainage] swdra = 0` (no `file=`), and that loading it standalone still validates:
```bash
cd tests/swap-cases/toml/1.hupselbrook
pixi run -e test python ../../../bmi/hello_swap.py ../../../../builddir/libswap_bmi.so   # uses swap.toml; sanity that the lib loads
```
(For a direct standalone check of `swap_coupled.toml`, temporarily point a hello-style script at it. swbotb=1 standalone needs a `gwltab`; in coupled mode the injection bypasses it — so a *standalone* run of `swap_coupled.toml` may warn about a missing GWL table. That is expected; the coupled path injects `gwl_injected`.)

- [ ] **Step 2: Write the coupler config TOML**

Create `tests/coupling/imod_coupler.toml`:
```toml
[[driver.coupling]]
mf6_model = "swapmf"
mf6_swap_recharge_pkg = "RCHA"

[driver.kernels.modflow6]
dll = "<ABS_PATH_TO>/libmf6.so"
work_dir = "./mf"

[driver.kernels.swap]
dll = "<ABS_PATH_TO>/builddir/libswap_bmi.so"
work_dir = "./swap"
```
(`run_coupled.py` rewrites the `dll` paths at runtime to absolute discovered paths; the literal placeholders document the schema. The exact TOML schema must match the forked `SwapModConfig` pydantic model — verify keys against `imod_coupler_fork`'s `config.py`.)

- [ ] **Step 3: Commit**

```bash
git add tests/swap-cases/toml/1.hupselbrook/swap_coupled.toml tests/coupling/imod_coupler.toml
git commit -m "feat(coupling): coupled SWAP config (swbotb=1, swdra=0) + imod_coupler config

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

# PHASE 8 — Orchestration + smoke

### Task 8.1: `run_coupled.py`

**Files:**
- Create: `tests/coupling/run_coupled.py`

- [ ] **Step 1: Write the orchestrator**

Create `tests/coupling/run_coupled.py`:
```python
"""Build the MODFLOW model, stage SWAP inputs, run the coupled SWAP<->MF6
simulation via the vendored swapmod driver, collect + plot results."""
import sys, shutil, os
from pathlib import Path
import numpy as np
import flopy

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE / "imod_coupler_fork"))
from build_modflow import build
from swapmod import SwapMod                      # vendored fork
from config import BaseConfig, SwapModConfig      # from the fork

def main(libswap: Path, libmf6: Path, case_dir: Path, run_dir: Path,
         ncol: int = 10, nper: int = 1096):
    run_dir.mkdir(parents=True, exist_ok=True)
    n_couple = ncol - 2                            # interior recharge cells = SWAP columns

    # 1. MODFLOW model
    build(run_dir / "mf", ncol=ncol, nper=nper)

    # 2. SWAP work dir
    swap_ws = run_dir / "swap"; swap_ws.mkdir(exist_ok=True)
    for f in ["swap_coupled.toml", "283.csv", "grassd.crp.toml",
              "maizes.crp.toml", "potatod.crp.toml"]:
        shutil.copy(case_dir / f, swap_ws / f)
    (swap_ws / "ensemble.txt").write_text(f"{n_couple}\n")

    # 3. Coupler config with discovered absolute dll paths
    cfg_text = (HERE / "imod_coupler.toml").read_text() \
        .replace("<ABS_PATH_TO>/libmf6.so", str(libmf6)) \
        .replace("<ABS_PATH_TO>/builddir/libswap_bmi.so", str(libswap))
    (run_dir / "imod_coupler.toml").write_text(cfg_text)

    # 4. Run the coupled loop
    os.chdir(run_dir)
    base = BaseConfig(log_level="INFO", timing=False)
    smcfg = SwapModConfig(config_dir=run_dir, **_toml_load(run_dir / "imod_coupler.toml")["driver"])
    driver = SwapMod(base, smcfg)
    driver.initialize()
    while driver.get_current_time() < driver.get_end_time():
        driver.update()
    driver.finalize()

    # 5. Collect + plot
    hds = flopy.utils.HeadFile(run_dir / "mf" / "swapmf.hds").get_alldata().squeeze()
    np.save(run_dir / "heads.npy", hds)
    _plot(hds, run_dir / "coupled_gwl.png")
    print("coupled run complete; heads shape:", hds.shape)

def _toml_load(p):
    import tomllib
    with open(p, "rb") as fh: return tomllib.load(fh)

def _plot(hds, out):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.plot(hds[-1], marker="o"); ax.set_xlabel("cell"); ax.set_ylabel("head [m]")
    ax.set_title("Final water table (two-channel + SWAP recharge)")
    fig.savefig(out, dpi=120); print("wrote", out)

if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve(),
         Path(sys.argv[3]).resolve(), Path(sys.argv[4]).resolve(),
         ncol=int(sys.argv[5]) if len(sys.argv) > 5 else 10,
         nper=int(sys.argv[6]) if len(sys.argv) > 6 else 1096)
```
*Adapt `BaseConfig`/`SwapModConfig` construction to the fork's exact API* (verify against `imod_coupler_fork/config.py` — the swap-branch `SwapModConfig.__init__(config_dir, **data)` chdirs to the config dir).

- [ ] **Step 2: First coupled run (short, to wire it up)**

Run:
```bash
cd /home/zawadzkim/Code/swap
pixi run build-linux
pixi run -e test python tests/coupling/run_coupled.py \
    builddir/libswap_bmi.so $(python -c "import flopy,pathlib;print(pathlib.Path(flopy.__file__))" >/dev/null; which mf6 | xargs dirname)/libmf6.so \
    tests/swap-cases/toml/1.hupselbrook /tmp/swapmf_run 10 30
```
Expected: completes; prints `coupled run complete`; writes `coupled_gwl.png`. Debug iteratively: if `solve` errors, call `get_last_bmi_error` via the wrapper; if recharge magnitudes look wrong, re-check the unit convention (Conventions section). Locate `libmf6.so` precisely for your environment (conda-forge `modflow6` installs it under the env `lib/`).

- [ ] **Step 3: Full run**

Run the same command with `nper=1096` (full 2002-2004). Expected: completes (slower); plot shows the water table between `h0`/`h1` with recharge-driven mounding in the interior, and interior heads varying.

- [ ] **Step 4: Commit**

```bash
git add tests/coupling/run_coupled.py
git commit -m "feat(coupling): run_coupled.py orchestration (flopy build + swapmod run + plot)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 8.2: Smoke verification + README

**Files:**
- Create: `tests/coupling/test_coupled_smoke.py`
- Create: `tests/coupling/README.md`
- Modify: `tests/coupling/meson.build` (register the smoke test, longer timeout)

- [ ] **Step 1: Write the smoke assertion test (short period)**

Create `tests/coupling/test_coupled_smoke.py`:
```python
"""Coupled smoke (qualitative): a short coupled run completes and the water
table shows a two-channel gradient with finite interior heads."""
import sys, tempfile
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
from run_coupled import main

def test_coupled_smoke(libswap, libmf6, case_dir):
    run = Path(tempfile.mkdtemp()) / "run"
    main(libswap, libmf6, case_dir, run, ncol=10, nper=30)
    hds = np.load(run / "heads.npy")
    final = hds[-1]
    assert np.all(np.isfinite(final))
    assert final[0] != final[-1], "no lateral gradient"
    assert (run / "coupled_gwl.png").exists()
    print("coupled smoke OK")

if __name__ == "__main__":
    test_coupled_smoke(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve(),
                       Path(sys.argv[3]).resolve())
```

- [ ] **Step 2: Register the meson test** in `tests/coupling/meson.build`:
```python
test('coupled-smoke',
    python3,
    args: [
        meson.current_source_dir() / 'test_coupled_smoke.py',
        swap_bmi_lib.full_path(),
        get_option('libmf6_path'),                      # provided at configure time
        meson.project_source_root() / 'tests/swap-cases/toml/1.hupselbrook',
    ],
    suite: 'coupling',
    depends: [swap_bmi_lib],
    timeout: 1200,
)
```
Add a `libmf6_path` option to `meson_options.txt`:
```
option('libmf6_path', type: 'string', value: '', description: 'Absolute path to libmf6.so for coupling test')
```
*(If wiring a meson option is awkward, keep `coupled-smoke` as a manually-run script and document it in the README instead of registering it — the smoke DoD is qualitative.)*

- [ ] **Step 3: Write the README**

Create `tests/coupling/README.md` documenting: prerequisites (xmipy, flopy, libmf6), the unit conventions (copy the Conventions block), how to run `run_coupled.py`, the two SWAP config deltas (swbotb=1, swdra=0), the `ensemble.txt` = `ncol-2` rule, and that the driver/wrapper are a vendored fork of imod_coupler's swap branch.

- [ ] **Step 4: Run the smoke**

Run: `pixi run -e test python tests/coupling/test_coupled_smoke.py builddir/libswap_bmi.so <libmf6.so> tests/swap-cases/toml/1.hupselbrook`
Expected: prints `coupled smoke OK`; PNG written.

- [ ] **Step 5: Commit**

```bash
git add tests/coupling/test_coupled_smoke.py tests/coupling/README.md tests/coupling/meson.build meson_options.txt
git commit -m "test(coupling): coupled smoke test + README

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

# PHASE 9 — Scoped polish + final gate

### Task 9.1: `real(8) → real64` on touched/new files

**Files:** the new modules (`swap_ensemble_mod.f90`, `swap_xmi_mod.f90`, `bmi_constants_mod.f90` already use `real64`/`c_*`) and any file edited in this arc that still uses `real(8)` (notably `boundbottom.f90` additions, `dtutil.f90` if touched). Do NOT sweep untouched files.

- [ ] **Step 1: Find `real(8)` in touched files only**

Run: `git diff --name-only development...HEAD -- 'src/**/*.f90' | xargs grep -ln "real(8)\|real\*8" 2>/dev/null`

- [ ] **Step 2: Replace declarations** (not literals) `real(8)` → `real(real64)`; ensure `use iso_fortran_env, only: real64`. Leave `1.0d0`-style literals alone (FP-neutral, avoids churn risk).

- [ ] **Step 3: Build + regression byte-identical**

Run: `pixi run build-linux && pixi run -e test check-fast`
Expected: 4 cases PASS byte-identical.

- [ ] **Step 4: Commit**

```bash
git commit -am "refactor(kinds): real(8) -> real64 on files touched by the coupling arc

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

### Task 9.2: `intent`/`pure` on touched leaf kernels

- [ ] **Step 1:** On procedures edited/added in this arc (e.g. `apply_coupled_gwl`, ensemble accessors), confirm explicit `intent` on every dummy arg; mark side-effect-free accessors `pure` where the compiler allows.
- [ ] **Step 2:** Build + `check-fast` byte-identical.
- [ ] **Step 3:** Commit `refactor(clarity): intent/pure on coupling-arc procedures`.

### Task 9.3: Final full-suite gate

- [ ] **Step 1: Clean rebuild + full regression**

Run: `pixi run clean && pixi run build-linux && pixi run -e test check-full`
Expected: pFUnit `OK (N tests)` (count includes the new suites); all regression cases PASS byte-identical.

- [ ] **Step 2: Re-run both coupling tests**

Run the SWAP-only XMI smoke (Task 5.1) and the coupled smoke (Task 8.2); confirm both pass and the plot is sensible.

- [ ] **Step 3: Final commit / branch ready**

```bash
git commit -am "chore(coupling): final gate — check-full byte-identical + coupling smokes green

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```
Branch `feat/swap-modflow-coupling` is ready to merge into `development` (use the finishing-a-development-branch skill).

---

## Self-review notes (gaps flagged for the implementer)

1. **`load_swap_config` shape** (Task 3.1 Step 4): confirm function-vs-subroutine and adapt — the single load-bearing unknown in the ensemble.
2. **`timecontrol_state_t` field names** (Task 4.2 Step 2): reuse `tstart`/`tend`/`t1900` exactly as the CAPI `swap_get_scalar` does.
3. **`Gash` vs `msw1eic`** (Task 0.1 Step 1): verify `Gash` (live) is distinct from the dead `msw1eic` before deleting; stop if entangled.
4. **Recharge units** (Conventions + Task 8.2): the #1 correctness risk — verify exchanged magnitudes (~mm/d) and sign (mounding under positive recharge) by inspection; adjust the conversion if MF6 RCHA semantics differ from the rate-in-m/d assumption.
5. **`SwapModConfig` API** (Tasks 7.2, 8.1): match the vendored fork's pydantic schema exactly.
6. **`libmf6.so` discovery** (Tasks 6.1, 8.1): pin the path for your env (conda-forge `modflow6`).
