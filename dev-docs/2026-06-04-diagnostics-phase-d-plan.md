# Diagnostics Phase D — per-instance graceful fatal termination — Plan

**Goal (the user's stated objective):** a critical condition during stepping must
**gracefully terminate the offending instance** — returning its collection of located
errors — **without killing the host process or the sibling columns**, uniformly across
entry points.

## Design decision (important — read this)

The design's ideal end-state was to migrate all ~109 compute-kernel `fatalerr_collected`
call sites to `state%diag%fatal(...) + early-return` and retire
`global_errors`/`fatalerr_collected`/`library_mode`. **That full per-site sweep is
deferred** as an explicit follow-up, because:
- Mechanically adding `return` after 109 physics-path guards risks the byte-identical law
  on error paths the 11 regression cases do **not** exercise (the oracle can't catch
  them), and risks latent crashes where a leaf `return` leaves `intent(out)` undefined.
- Many sites are state-less leaves (solver, `WC_K_models`, `QROMBD`) with no `state` in
  scope, needing argument threading — a large, delicate change.

**Phase D instead delivers the full capability with a small, additive, byte-identical
change:** an *active-error-sink* registered by the step driver. While an instance is
stepping, every existing `fatalerr_collected` call records into **that instance's**
`state%diag%errors` and sets its `fatal_raised`, then returns (no `error stop`); the step
driver checks `aborted()` at a boundary and the entry point translates. So all 109 sites
become per-instance and host-safe **without editing any of them**. `global_errors` /
`fatalerr_collected` / `library_mode` are retained (load-phase + global fallback);
retiring them via the per-site sweep is the documented follow-up.

Residual vs ideal: a kernel fatal still runs the remainder of the current sub-phase before
the next boundary check (the "garbage window" the design acknowledged) — but this is
**bounded per sub-phase** and strictly better than today's embedded `library_mode`
(which continued globally with a single shared flag). The CLI keeps prompt halting (it
`error stop`s at the post-step check).

**Base:** branch `diagnostics-bcde` (Phase C in). Conventions as before (commit `-c
diff.ignoreSubmodules=all --no-verify`; targeted add; pFUnit gotcha). Byte-identical:
nothing here triggers in the regression cases.

---

## Task D1 — active-error-sink in `error_mod` + redirect `fatalerr_collected`

**File:** `src/error/error.f90`. Test: `tests/unit/error/test_error.pf` (or a new
`test_error_active_sink.pf` registered in all 3 places).

`error_collection_t` already lives in `error_mod`, so error_mod can hold pointers to a
collection + a flag with NO new dependency (no diagnostics_mod import — avoids a cycle).

Add module-private state + public setters:
```fortran
   type(error_collection_t), pointer, private :: active_errors => null()
   logical,                  pointer, private :: active_fatal  => null()

   public :: set_active_error_sink, clear_active_error_sink
```
```fortran
   !> Register a per-instance sink (the stepping instance's state%diag%errors +
   !! fatal_raised). While set, fatalerr_collected records into it and returns
   !! instead of aborting — so an embedded host's other columns keep running and
   !! the host process is not killed. The driver sets this at each step and the
   !! pointers stay valid for the step's duration (state outlives the step).
   subroutine set_active_error_sink(errors, fatal_flag)
      type(error_collection_t), pointer, intent(in) :: errors
      logical,                  pointer, intent(in) :: fatal_flag
      active_errors => errors
      active_fatal  => fatal_flag
   end subroutine set_active_error_sink

   subroutine clear_active_error_sink()
      active_errors => null()
      active_fatal  => null()
   end subroutine clear_active_error_sink
```
Change `fatalerr_collected` to prefer the active sink:
```fortran
   subroutine fatalerr_collected(routine, message)
      character(len=*), intent(in) :: routine
      character(len=*), intent(in) :: message
      if (associated(active_errors)) then
         call active_errors%append(ERR_LEGACY_FATAL, message, routine)
         if (associated(active_fatal)) active_fatal = .true.
         return
      end if
      call global_errors%append(ERR_LEGACY_FATAL, message, routine)
      call global_errors%abort_if_fatal()
   end subroutine fatalerr_collected
```

**Tests** (TDD): in a new/existing pf:
- After `set_active_error_sink(ep, fp)` (with a local `type(error_collection_t),
  pointer :: ep` and `logical, pointer :: fp` pointing at a local target collection +
  logical), `call fatalerr_collected('K','boom')` → the local collection `count()==1`,
  the local flag `.true.`, and it did NOT abort (test reaches the next line). Then
  `clear_active_error_sink()`; assert a subsequent call with no sink would use the global
  path (don't actually call it — that aborts; just assert `clear` nullified by checking
  behavior indirectly, or skip). Keep the test to the active-sink path only.
- Build, `pixi run -e test test-pfunit` (count +1), commit. (No regression run needed —
  pure addition not on any stepping path yet.)

Register the new pf in `pf_files` + `pfunit_extra_sources` (it only needs
`src/error/error.f90`, already present) + `testSuites.inc` if you make a new file; or just
append @tests to the existing `tests/unit/error/test_error.pf` (simpler — no registration
needed).

---

## Task D2 — driver registers the sink + boundary check

**File:** `src/driver/swap_mod.f90` (`swap_run_step`).

1. Add the `target` attribute to the `state` dummy: `type(swap_state_t), intent(inout),
   target :: state`.
2. Add `use error_mod, only: error_collection_t, set_active_error_sink` to the routine's
   `use` list (it already uses other modules at the top of the subroutine).
3. Right AFTER the existing `call state%diag%set_simtime(...)` (and before the
   `associate`), register the sink:
   ```fortran
      block
         type(error_collection_t), pointer :: diag_errors
         logical,                  pointer :: diag_fatal
         diag_errors => state%diag%errors
         diag_fatal  => state%diag%fatal_raised
         call set_active_error_sink(diag_errors, diag_fatal)
      end block
   ```
4. Add ONE boundary check that stops the rest of the step once a fatal has been recorded.
   Place it immediately AFTER `call soilwater_update(state)` (the end of the core
   hydrology), before `if (time%flTemperature) ...`:
   ```fortra n
         if (state%diag%aborted()) return
   ```
   (Use `fortran`, not `fortra n` — that was a typo guard; write `if (state%diag%aborted()) return`.)

Clean rebuild NOT required (no schema change), but if a stale-.mod issue appears, clean.
Build, `pixi run -e test check-fast` 4/4 byte-identical (no abort fires), `test-xmi`,
commit.

---

## Task D3 — entry-point translation

Each entry point checks `state%diag%aborted()` after a step and translates. Byte-identical
(never fires in regression).

1. **CLI** `src/driver/swap_main.f90`: in the `do while (.not. state%timecontrol%flRunEnd)`
   loop, after `call swap_run_step(state, config)`:
   ```fortran
      if (state%diag%aborted()) then
         write(error_unit,'(A)') state%diag%errors%summary()
         error stop 'fatal error during simulation'
      end if
   ```
   Add `use iso_fortran_env, only: error_unit` to the program.
2. **BMI** `src/bindings/swap_bmi_mod.f90` `bmi_update`: after `call
   swap_run_step(bmi_state, bmi_config)`:
   ```fortran
      if (bmi_state%diag%aborted()) then; rc = 1; return; end if
      rc = 0
   ```
   (replace the existing `rc = 0`).
3. **C-API** `src/bindings/swap_capi_mod.f90`: find the function(s) that call
   `swap_run_step(capi_state, ...)` (grep `swap_run_step` in the file). After each, add
   `if (capi_state%diag%aborted()) then; ierr = <nonzero>; return; end if` before the
   success `ierr = 0`. If no C-API function calls swap_run_step directly (stepping goes
   through BMI's bmi_update), note that and skip.
4. **Ensemble** `src/driver/swap_ensemble_mod.f90` `ensemble_step_day`: inside the
   per-column `do i = 1, ncol` loop, after the inner `call swap_run_step(columns(i), ...)`
   (and its existing `library_fatal_raised` check), add a per-column abort handling:
   ```fortran
            if (columns(i)%diag%aborted()) then
               rc = 1
               exit            ! stop this column's day; other columns already
                               ! stepped/te-step independently — see note
            end if
   ```
   Keep the existing `library_fatal_raised()` check too (load/global path). The intent:
   an aborted column reports rc=1 and stops its own substep loop; sibling columns (other
   `i`) are independent iterations and keep running. Confirm the loop structure supports
   this (the `do i` loop continues to the next column; only the inner `do` substep loop is
   exited). If the structure makes per-column isolation awkward, implement the minimal
   safe version (set rc=1 and exit the inner substep loop for that column) and report.

Build, `pixi run -e test check-fast` 4/4 byte-identical, `test-xmi`, commit.

---

## Phase D done-when
- A fatal raised by any existing `fatalerr_collected` site DURING stepping records into the
  stepping instance's `state%diag%errors`, sets its `fatal_raised`, and does NOT
  `error stop` — the step driver stops that step at the boundary and the entry point
  translates (CLI: summary→stderr + exit; BMI/C-API: nonzero rc; ensemble: column rc=1,
  siblings continue).
- `check-full` 11/0/7 byte-identical (nothing fires in regression).
- Documented residual: full per-site `state%diag%fatal` sweep + retiring
  `global_errors`/`fatalerr_collected`/`library_mode` is a follow-up.
