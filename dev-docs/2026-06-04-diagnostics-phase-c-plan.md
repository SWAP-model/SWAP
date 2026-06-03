# Diagnostics Phase C — state-borne diagnostics_t — Plan

**Goal:** Put a per-instance `diagnostics_t` on `swap_state_t` (`state%diag`) carrying the
per-instance error collection + fatal flag + sim-time + instance-id, with leveled-emit and
`fatal`/`aborted` methods. This is the **spine Phase D needs**. Byte-identical;
`swap_state_t` schema change → **clean rebuild mandatory**.

**Scope decision (documented):** the log *sink* stays process-global (one file/level);
per-instance disambiguation is achieved by **instance-id + sim-date stamping** on every
record, not by per-instance log files. Full per-instance sink isolation is deferred as
optional polish — the high-value per-instance state is the error collection + fatal flag,
which this phase delivers. This satisfies the design's "threaded on state" for the
error/fatal spine and stamping.

**Base:** branch `diagnostics-bcde` (Phase B in). Conventions as Phase B (commit with
`git -c diff.ignoreSubmodules=all commit --no-verify`; targeted add; pFUnit gotcha).

---

## Task C1 — define `diagnostics_t` (type + methods + tests), no state change yet

**File:** `src/error/diagnostics.f90` (the foundation module — it already sits above
`error_mod`/`swap_log` and below `state`). Test: `tests/unit/error/test_diagnostics.pf`.

Add imports to the module: extend `use swap_log` to also bring `log_debug, log_info,
log_warn, log_error`; add `use error_mod, only: error_collection_t, ERR_LEGACY_FATAL`.

Add the type (public) + its TBPs:

```fortran
   !> Per-instance diagnostics carried on swap_state_t (state%diag). Owns the
   !! instance's error accumulation + fatal flag + sim-time context; leveled
   !! emits route through the process-global swap_log sink, stamped with
   !! instance-id and sim-date so interleaved multi-instance logs are
   !! attributable.
   type, public :: diagnostics_t
      integer                  :: instance_id  = 0
      character(len=11)        :: sim_date     = ''
      integer                  :: daynr        = 0
      integer                  :: daycum       = 0
      type(error_collection_t) :: errors
      logical                  :: fatal_raised = .false.
   contains
      procedure :: debug       => diag_debug
      procedure :: info        => diag_info
      procedure :: warn        => diag_warn
      procedure :: error       => diag_error
      procedure :: fatal       => diag_fatal
      procedure :: aborted     => diag_aborted
      procedure :: set_simtime => diag_set_simtime
      procedure, private :: stamp => diag_stamp
   end type diagnostics_t
```

Methods in `contains`:

```fortran
   !> Prepend "[#id] date " to a message (id omitted when 0, date when empty)
   !! so single-instance CLI logs stay clean and ensemble logs are attributable.
   function diag_stamp(self, message) result(s)
      class(diagnostics_t), intent(in) :: self
      character(len=*),     intent(in) :: message
      character(len=:), allocatable    :: s
      character(len=16) :: idbuf
      s = message
      if (len_trim(self%sim_date) > 0) s = trim(self%sim_date)//' '//s
      if (self%instance_id /= 0) then
         write(idbuf,'("[#",I0,"] ")') self%instance_id
         s = trim(idbuf)//s
      end if
   end function diag_stamp

   subroutine diag_debug(self, context, message)
      class(diagnostics_t), intent(in) :: self
      character(len=*),     intent(in) :: context, message
      call log_debug(context, self%stamp(message))
   end subroutine diag_debug

   subroutine diag_info(self, context, message)
      class(diagnostics_t), intent(in) :: self
      character(len=*),     intent(in) :: context, message
      call log_info(context, self%stamp(message))
   end subroutine diag_info

   subroutine diag_warn(self, context, message)
      class(diagnostics_t), intent(in) :: self
      character(len=*),     intent(in) :: context, message
      call log_warn(context, self%stamp(message))
   end subroutine diag_warn

   subroutine diag_error(self, context, message)
      class(diagnostics_t), intent(in) :: self
      character(len=*),     intent(in) :: context, message
      call log_error(context, self%stamp(message))
   end subroutine diag_error

   !> Record a fatal into the per-instance collection (auto-logs via log_error)
   !! and set the sticky fatal flag. Does NOT abort — the step driver checks
   !! aborted() at substep boundaries (Phase D).
   subroutine diag_fatal(self, context, message)
      class(diagnostics_t), intent(inout) :: self
      character(len=*),     intent(in)    :: context, message
      call self%errors%append(ERR_LEGACY_FATAL, self%stamp(message), context)
      self%fatal_raised = .true.
   end subroutine diag_fatal

   pure logical function diag_aborted(self)
      class(diagnostics_t), intent(in) :: self
      diag_aborted = self%fatal_raised
   end function diag_aborted

   subroutine diag_set_simtime(self, date, daynr, daycum)
      class(diagnostics_t), intent(inout) :: self
      character(len=*),     intent(in)    :: date
      integer,              intent(in)    :: daynr, daycum
      self%sim_date = date
      self%daynr    = daynr
      self%daycum   = daycum
   end subroutine diag_set_simtime
```

**Tests** (TDD, standalone — no state): add @tests:
- `test_diag_fatal_records_and_aborts`: fresh `diagnostics_t` → `aborted()` false; after
  `call d%fatal('k','msg')` → `aborted()` true and `d%errors%count() == 1`.
- `test_diag_warn_does_not_abort`: after `set_simtime('2002-01-01',1,1)` + `warn('c','m')`
  → `aborted()` false, `errors%count() == 0`.
(These emit through log_* to stdout/stderr — harmless test noise; that's fine.)

Build, `pixi run -e test test-pfunit` (OK count +2), no regression run needed (no state
change yet), commit.

---

## Task C2 — add `state%diag`, set instance-id, wire sim-time

**Files:** `src/state/swap_state.f90`, `src/driver/swap_ensemble_mod.f90`,
`src/driver/swap_mod.f90`. **Clean rebuild required.**

1. `src/state/swap_state.f90`: `use diagnostics_mod, only: diagnostics_t` and add the
   field to `swap_state_t`:
   ```fortran
      type(diagnostics_t)        :: diag
   ```
   (place it last in the type). `diagnostics_t` is small (scalars + an allocatable-backed
   collection) — stack-safe.

2. `src/driver/swap_ensemble_mod.f90` `ensemble_init`: after the per-column init loop
   (`call swap_init_from_loaded_config(columns(i), ...)`), set the per-column id:
   ```fortran
         columns(i)%diag%instance_id = i
   ```
   (single-instance CLI/BMI paths keep the default id 0 → unstamped, clean.)

3. `src/driver/swap_mod.f90` `swap_run_step(state, config)`: at the TOP of the routine
   body (after declarations, before the step work), stamp the current sim-time:
   ```fortran
      call state%diag%set_simtime(state%timecontrol%date, &
                                  state%timecontrol%daynr, state%timecontrol%daycum)
   ```
   (This makes every diag emit during the step carry the date. No behavior change to
   outputs.)

**Gate:** `pixi run clean && pixi run build-linux` (schema change → clean rebuild), then
`pixi run -e test check-fast` (4/4 byte-identical) + `pixi run -e test test-xmi`. Commit.

---

## Phase C done-when
- `state%diag` exists on every instance; ensemble columns get distinct `instance_id`;
  sim-time is stamped each step.
- `diagnostics_t` provides `fatal`/`aborted`/`errors` (the Phase D spine) + leveled emits
  with instance-id/sim-date stamping, unit-tested.
- `check-fast` 4/4 byte-identical after a clean rebuild; new tests counted.

## Out of scope (Phase D)
- Migrating kernel `fatalerr_collected` sites to `state%diag%fatal`; the substep-boundary
  `aborted()` checks; entry-point translation (CLI exit / library rc / per-instance XMI
  poll); retiring `global_errors`/`fatalerr_collected`/`FatalERR`/`library_mode`.
