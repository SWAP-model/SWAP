# Diagnostics Phase A — Foundation & Wiring — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make SWAP logging user-configurable (env + C-API), initialised in *every* entry point with embedding-safe defaults, and route WARN/ERROR to stderr — without changing kernels, physics, or the `swap_state_t` schema, byte-identical against `swap420gf`.

**Architecture:** A new pure `diagnostics_mod` (`src/error/diagnostics.f90`) defines a resolved-config type and a precedence merge (C-API > env > TOML > default). `swap_log` gains stderr routing and stdout/stderr setters. Each entry point (CLI, BMI, C-API, ensemble) builds a config via the resolver and calls `log_init`; embedded paths default `to_stdout=.false.`. The `[logging]` TOML block and the full state-borne `diagnostics_t` are later phases.

**Tech Stack:** Modern Fortran (gfortran, `iso_fortran_env`), Meson build, pFUnit unit tests. Spec: `dev-docs/2026-06-03-diagnostics-arc-design.md` (§4f config surface, §4e routing, §7 Phase A).

---

## Conventions for every task

- **Build:** `pixi run build-linux`
- **Fast unit+regression gate (must be green before each commit):** `pixi run -e test check-fast`
- **Unit suite only:** `pixi run -e test test-pfunit`
- **A new `.pf` must be registered in 3 places or it runs dark:** `tests/unit/meson.build` `pf_files` list, `tests/unit/meson.build` `pfunit_extra_sources` (for any new src module it uses), and `tests/unit/testSuites.inc` (`ADD_TEST_SUITE`). Suite name = file basename + `_suite` (e.g. `error/test_diagnostics.pf` → `test_diagnostics_suite`). Confirm the `OK (N tests)` count rises.
- **No `src/state/*` schema change in this phase**, so no forced clean rebuild; but if Meson acts stale, `pixi run clean` then rebuild.
- **Commit messages:** `feat(diagnostics): …` / `refactor(diagnostics): …`, body explains *why*. End with the Co-Authored-By trailer.
- **Targeted `git add`** (never `git add -A`).

---

## File Structure

| File | Responsibility | Action |
|------|----------------|--------|
| `src/error/diagnostics.f90` | Pure config types + precedence resolver + env reader + default configs | **Create** |
| `src/core/swap_log.f90` | Add stderr routing + stdout/stderr setters + `to_stderr` init arg | Modify |
| `src/bindings/swap_capi_mod.f90` | bind(C) `swap_set_log_level`/`swap_set_log_stdout`; init logging in C-API init paths | Modify |
| `src/bindings/swap_bmi_mod.f90` | Init logging (embedded defaults) in `bmi_initialize` | Modify |
| `src/driver/swap_ensemble_mod.f90` | Init logging (embedded defaults) in `ensemble_init` | Modify |
| `src/driver/swap_main.f90` | Route CLI `log_init` through the resolver | Modify |
| `meson.build` | Register `src/error/diagnostics.f90` in the build | Modify |
| `tests/unit/error/test_diagnostics.pf` | Unit tests for the resolver/parser/routing | **Create** |
| `tests/unit/meson.build` | Register the new `.pf` + new src module | Modify |
| `tests/unit/testSuites.inc` | `ADD_TEST_SUITE(test_diagnostics_suite)` | Modify |

---

## Task 1: `diagnostics_mod` — config types + level-name parsing

**Files:**
- Create: `src/error/diagnostics.f90`
- Modify: `meson.build` (register the module)
- Test: `tests/unit/error/test_diagnostics.pf`
- Modify: `tests/unit/meson.build`, `tests/unit/testSuites.inc`

- [ ] **Step 1: Create the module with the types and `log_level_from_name`**

Create `src/error/diagnostics.f90`:

```fortran
!> Diagnostics configuration layer (Phase A).
!!
!! Pure resolution of logging settings from layered sources with the
!! precedence  C-API > env > TOML > built-in default  (spec §4f). This
!! phase owns only the *config* (level/file/routing) + a precedence merge;
!! the per-instance state-borne diagnostics record is a later phase.
module diagnostics_mod
   use swap_log, only: LOGLEVEL_DEBUG, LOGLEVEL_INFO, LOGLEVEL_WARN, &
                       LOGLEVEL_ERROR, LOGLEVEL_NONE
   implicit none
   private

   public :: diagnostics_config_t, diag_overrides_t
   public :: log_level_from_name, merge_overrides, resolve_diagnostics_config
   public :: read_env_overrides, default_cli_config, default_embedded_config

   !> Fully-resolved logging settings handed to log_init.
   type :: diagnostics_config_t
      integer                       :: level      = LOGLEVEL_INFO
      logical                       :: to_stdout  = .true.
      logical                       :: to_stderr  = .true.
      logical                       :: timestamps = .false.
      character(len=:), allocatable :: log_file       ! unallocated => no file
   end type diagnostics_config_t

   !> A sparse set of overrides from one source. Only fields whose
   !! has_* flag is .true. override the base config in merge_overrides.
   type :: diag_overrides_t
      logical                       :: has_level      = .false.
      integer                       :: level          = LOGLEVEL_INFO
      logical                       :: has_stdout     = .false.
      logical                       :: to_stdout      = .true.
      logical                       :: has_timestamps = .false.
      logical                       :: timestamps     = .false.
      logical                       :: has_file       = .false.
      character(len=:), allocatable :: log_file
   end type diag_overrides_t

contains

   !> Map a level name (case-insensitive) to a LOGLEVEL_* value.
   !! Unknown/empty => LOGLEVEL_INFO (safe default).
   pure function log_level_from_name(name) result(level)
      character(len=*), intent(in) :: name
      integer :: level
      select case (upcase(trim(adjustl(name))))
      case ('DEBUG');           level = LOGLEVEL_DEBUG
      case ('INFO');            level = LOGLEVEL_INFO
      case ('WARN', 'WARNING'); level = LOGLEVEL_WARN
      case ('ERROR');           level = LOGLEVEL_ERROR
      case ('NONE', 'OFF');     level = LOGLEVEL_NONE
      case default;             level = LOGLEVEL_INFO
      end select
   end function log_level_from_name

   !> ASCII upper-case helper (pure, no locale).
   pure function upcase(s) result(u)
      character(len=*), intent(in) :: s
      character(len=len(s))        :: u
      integer :: i, c
      do i = 1, len(s)
         c = iachar(s(i:i))
         if (c >= iachar('a') .and. c <= iachar('z')) then
            u(i:i) = achar(c - 32)
         else
            u(i:i) = s(i:i)
         end if
      end do
   end function upcase

end module diagnostics_mod
```

- [ ] **Step 2: Register the module in the main build**

In `meson.build`, the `src/error/` group is at lines 98–99:

```
    'src/error/error.f90',
    'src/error/fatalerr.f90',
```

`diagnostics_mod` uses `swap_log` only, so it can sit right after `error.f90`. Change to:

```
    'src/error/error.f90',
    'src/error/diagnostics.f90',
    'src/error/fatalerr.f90',
```

(`src/core/swap_log.f90` is already listed earlier in the file, so the `swap_log` dependency is satisfied.)

- [ ] **Step 3: Write the failing test**

Create `tests/unit/error/test_diagnostics.pf`:

```fortran
! Tests for diagnostics_mod (Phase A): level-name parsing + precedence merge.

@test
subroutine test_log_level_from_name_known_values()
   use funit
   use diagnostics_mod, only: log_level_from_name
   use swap_log,        only: LOGLEVEL_DEBUG, LOGLEVEL_INFO, LOGLEVEL_WARN, &
                              LOGLEVEL_ERROR, LOGLEVEL_NONE
   implicit none
   @assertEqual(LOGLEVEL_DEBUG, log_level_from_name('debug'))
   @assertEqual(LOGLEVEL_INFO,  log_level_from_name('INFO'))
   @assertEqual(LOGLEVEL_WARN,  log_level_from_name('Warn'))
   @assertEqual(LOGLEVEL_WARN,  log_level_from_name('warning'))
   @assertEqual(LOGLEVEL_ERROR, log_level_from_name('ERROR'))
   @assertEqual(LOGLEVEL_NONE,  log_level_from_name('off'))
end subroutine

@test
subroutine test_log_level_from_name_unknown_defaults_to_info()
   use funit
   use diagnostics_mod, only: log_level_from_name
   use swap_log,        only: LOGLEVEL_INFO
   implicit none
   @assertEqual(LOGLEVEL_INFO, log_level_from_name('gibberish'))
   @assertEqual(LOGLEVEL_INFO, log_level_from_name(''))
end subroutine
```

- [ ] **Step 4: Register the test suite**

In `tests/unit/meson.build`, add to the `pf_files` list (near the existing `'error/test_error.pf',` line ~223):

```
        'error/test_diagnostics.pf',
```

In the `pfunit_extra_sources` list (the block that already contains `'../../src/error/error.f90',`), add directly after it:

```
        '../../src/error/diagnostics.f90',
```

In `tests/unit/testSuites.inc`, after the `ADD_TEST_SUITE(test_error_suite)` line, add:

```
ADD_TEST_SUITE(test_diagnostics_suite)
```

- [ ] **Step 5: Run the unit suite — verify the new tests pass and the count rose**

Run: `pixi run -e test test-pfunit`
Expected: build succeeds; suite runs; `OK (N tests)` with N two higher than before; the two `test_diagnostics_*` tests pass.

- [ ] **Step 6: Commit**

```bash
git add src/error/diagnostics.f90 meson.build \
        tests/unit/error/test_diagnostics.pf tests/unit/meson.build tests/unit/testSuites.inc
git commit -m "feat(diagnostics): add diagnostics_mod config types + level-name parsing

Phase A foundation: diagnostics_config_t / diag_overrides_t and a
case-insensitive log_level_from_name. Pure, no I/O, no state change.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 2: Precedence merge + resolver

**Files:**
- Modify: `src/error/diagnostics.f90`
- Test: `tests/unit/error/test_diagnostics.pf`

- [ ] **Step 1: Write the failing test (precedence)**

Append to `tests/unit/error/test_diagnostics.pf`:

```fortran
@test
subroutine test_merge_overrides_applies_only_set_fields()
   use funit
   use diagnostics_mod, only: diagnostics_config_t, diag_overrides_t, merge_overrides
   use swap_log,        only: LOGLEVEL_INFO, LOGLEVEL_DEBUG
   implicit none
   type(diagnostics_config_t) :: base, out
   type(diag_overrides_t)     :: ov
   base%level = LOGLEVEL_INFO
   base%to_stdout = .true.
   ! Only override the level; to_stdout must be untouched.
   ov%has_level = .true.
   ov%level     = LOGLEVEL_DEBUG
   out = merge_overrides(base, ov)
   @assertEqual(LOGLEVEL_DEBUG, out%level)
   @assertTrue(out%to_stdout)
end subroutine

@test
subroutine test_resolve_precedence_capi_over_env_over_toml_over_base()
   use funit
   use diagnostics_mod, only: diagnostics_config_t, diag_overrides_t, &
                              resolve_diagnostics_config
   use swap_log,        only: LOGLEVEL_INFO, LOGLEVEL_DEBUG, LOGLEVEL_WARN, LOGLEVEL_ERROR
   implicit none
   type(diagnostics_config_t) :: base, out
   type(diag_overrides_t)     :: toml, env, capi
   base%level = LOGLEVEL_INFO
   toml%has_level = .true.; toml%level = LOGLEVEL_WARN
   env%has_level  = .true.; env%level  = LOGLEVEL_DEBUG
   capi%has_level = .true.; capi%level = LOGLEVEL_ERROR
   out = resolve_diagnostics_config(base, toml, env, capi)
   @assertEqual(LOGLEVEL_ERROR, out%level)   ! capi wins

   ! With capi absent, env wins over toml.
   capi%has_level = .false.
   out = resolve_diagnostics_config(base, toml, env, capi)
   @assertEqual(LOGLEVEL_DEBUG, out%level)
end subroutine
```

- [ ] **Step 2: Run to verify it fails**

Run: `pixi run -e test test-pfunit`
Expected: compile error / link error — `merge_overrides` and `resolve_diagnostics_config` not yet defined.

- [ ] **Step 3: Implement merge + resolve**

In `src/error/diagnostics.f90`, add these procedures inside `contains` (after `log_level_from_name`):

```fortran
   !> Apply a sparse override set onto a base config (only set fields win).
   pure function merge_overrides(base, ov) result(cfg)
      type(diagnostics_config_t), intent(in) :: base
      type(diag_overrides_t),     intent(in) :: ov
      type(diagnostics_config_t)             :: cfg
      cfg = base
      if (ov%has_level)      cfg%level      = ov%level
      if (ov%has_stdout)     cfg%to_stdout  = ov%to_stdout
      if (ov%has_timestamps) cfg%timestamps = ov%timestamps
      if (ov%has_file)       cfg%log_file   = ov%log_file
   end function merge_overrides

   !> Resolve the final config with precedence  capi > env > toml > base.
   pure function resolve_diagnostics_config(base, toml, env, capi) result(cfg)
      type(diagnostics_config_t), intent(in) :: base
      type(diag_overrides_t),     intent(in) :: toml, env, capi
      type(diagnostics_config_t)             :: cfg
      cfg = merge_overrides(base, toml)
      cfg = merge_overrides(cfg,  env)
      cfg = merge_overrides(cfg,  capi)
   end function resolve_diagnostics_config
```

- [ ] **Step 4: Run to verify it passes**

Run: `pixi run -e test test-pfunit`
Expected: PASS; `OK (N tests)` count rose by 2.

- [ ] **Step 5: Commit**

```bash
git add src/error/diagnostics.f90 tests/unit/error/test_diagnostics.pf
git commit -m "feat(diagnostics): precedence merge + resolver (capi>env>toml>default)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 3: `swap_log` stderr routing

**Files:**
- Modify: `src/core/swap_log.f90`
- Test: `tests/unit/error/test_diagnostics.pf` (pure routing helper)

Rationale: today every level writes to `*` (stdout). Spec §4e wants WARN/ERROR on stderr and DEBUG/INFO on stdout only when enabled. We add a *pure decision helper* (unit-testable) and wire `log_message` to it. The file output path is unchanged, so existing `test_swap_log` (which reads the file) stays green.

- [ ] **Step 1: Write the failing test for the routing decision**

Append to `tests/unit/error/test_diagnostics.pf`:

```fortran
@test
subroutine test_log_destination_warn_error_go_to_stderr()
   use funit
   use swap_log, only: log_destination, LOGLEVEL_DEBUG, LOGLEVEL_INFO, &
                       LOGLEVEL_WARN, LOGLEVEL_ERROR
   implicit none
   logical :: out, err
   ! WARN with stderr enabled -> stderr only, not stdout
   call log_destination(LOGLEVEL_WARN, .true., .true., out, err)
   @assertFalse(out); @assertTrue(err)
   ! ERROR likewise
   call log_destination(LOGLEVEL_ERROR, .true., .true., out, err)
   @assertFalse(out); @assertTrue(err)
   ! INFO with stdout enabled -> stdout only
   call log_destination(LOGLEVEL_INFO, .true., .true., out, err)
   @assertTrue(out); @assertFalse(err)
   ! INFO with stdout disabled -> neither stream
   call log_destination(LOGLEVEL_INFO, .false., .true., out, err)
   @assertFalse(out); @assertFalse(err)
   ! WARN with stderr disabled -> neither stream
   call log_destination(LOGLEVEL_WARN, .true., .false., out, err)
   @assertFalse(out); @assertFalse(err)
end subroutine
```

- [ ] **Step 2: Run to verify it fails**

Run: `pixi run -e test test-pfunit`
Expected: compile/link error — `log_destination` not defined / not public.

- [ ] **Step 3: Implement the helper and wire `log_message`**

In `src/core/swap_log.f90`:

(a) add the stderr unit + module flag and the new publics. Change the `use` line (currently the module has no `iso_fortran_env` import — add one) by inserting after `module swap_log` / `implicit none`:

```fortran
module swap_log
    use iso_fortran_env, only: output_unit, error_unit
    implicit none
    private
```

(b) add a module flag next to `log_to_stdout`:

```fortran
    logical :: log_to_stdout = .true.
    logical :: log_to_stderr = .true.
```

(c) add to the public list (next to `log_set_level`):

```fortran
    public :: log_set_level
    public :: log_set_stdout
    public :: log_destination
```

(d) add the pure helper in `contains`:

```fortran
    pure subroutine log_destination(level, to_stdout, to_stderr, write_stdout, write_stderr)
        !> Decide which console streams a record at `level` goes to.
        !! WARN/ERROR -> stderr (if enabled); DEBUG/INFO -> stdout (if enabled).
        integer, intent(in)  :: level
        logical, intent(in)  :: to_stdout, to_stderr
        logical, intent(out) :: write_stdout, write_stderr
        if (level >= LOGLEVEL_WARN) then
            write_stdout = .false.
            write_stderr = to_stderr
        else
            write_stdout = to_stdout
            write_stderr = .false.
        end if
    end subroutine log_destination

    subroutine log_set_stdout(flag)
        !> Enable/disable console stdout output at runtime (entry points use
        !! this to silence stdout when embedded).
        logical, intent(in) :: flag
        log_to_stdout = flag
    end subroutine log_set_stdout
```

(e) replace the stdout block in `log_message` (currently):

```fortran
        ! Output to stdout
        if (log_to_stdout) then
            write(*,'(A)') trim(formatted_msg)
        end if
```

with:

```fortran
        ! Output to console: WARN/ERROR -> stderr, DEBUG/INFO -> stdout.
        block
            logical :: to_out, to_err
            call log_destination(level, log_to_stdout, log_to_stderr, to_out, to_err)
            if (to_out) write(output_unit,'(A)') trim(formatted_msg)
            if (to_err) write(error_unit, '(A)') trim(formatted_msg)
        end block
```

- [ ] **Step 4: Run the unit suite — new routing test passes, existing swap_log tests still green**

Run: `pixi run -e test test-pfunit`
Expected: PASS; `test_swap_log_*` unchanged (they read the file, not the console); `test_log_destination_*` passes.

- [ ] **Step 5: Full fast gate (byte-identical regression)**

Run: `pixi run -e test check-fast`
Expected: pFUnit green + 4 regression cases byte-identical. (Console routing does not touch compared output files.)

- [ ] **Step 6: Commit**

```bash
git add src/core/swap_log.f90 tests/unit/error/test_diagnostics.pf
git commit -m "feat(diagnostics): route WARN/ERROR to stderr in swap_log

Adds pure log_destination decision helper + log_set_stdout setter; wires
log_message to send WARN/ERROR to error_unit and DEBUG/INFO to stdout only
when enabled. File output path unchanged; byte-identical regression.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 4: `log_init` gains `to_stderr`; env reader + default configs

**Files:**
- Modify: `src/core/swap_log.f90` (optional `to_stderr` arg on `log_init`)
- Modify: `src/error/diagnostics.f90` (`read_env_overrides`, `default_cli_config`, `default_embedded_config`)
- Test: `tests/unit/error/test_diagnostics.pf` (default configs)

- [ ] **Step 1: Add the optional `to_stderr` arg to `log_init`**

In `src/core/swap_log.f90`, the `log_init` signature is:

```fortran
    subroutine log_init(log_level, log_file, to_stdout, timestamps)
        integer, intent(in), optional :: log_level
        character(len=*), intent(in), optional :: log_file
        logical, intent(in), optional :: to_stdout
        logical, intent(in), optional :: timestamps
```

Add a `to_stderr` optional argument:

```fortran
    subroutine log_init(log_level, log_file, to_stdout, timestamps, to_stderr)
        integer, intent(in), optional :: log_level
        character(len=*), intent(in), optional :: log_file
        logical, intent(in), optional :: to_stdout
        logical, intent(in), optional :: timestamps
        logical, intent(in), optional :: to_stderr
```

And inside the body, after the `to_stdout` handling block, add:

```fortran
        ! Set stderr option
        if (present(to_stderr)) then
            log_to_stderr = to_stderr
        else
            log_to_stderr = .true.
        end if
```

(Existing callers pass it positionally/by keyword without `to_stderr`; the default `.true.` preserves behaviour.)

- [ ] **Step 2: Write the failing test for default configs**

Append to `tests/unit/error/test_diagnostics.pf`:

```fortran
@test
subroutine test_default_configs_cli_vs_embedded()
   use funit
   use diagnostics_mod, only: diagnostics_config_t, default_cli_config, &
                              default_embedded_config
   use swap_log,        only: LOGLEVEL_INFO
   implicit none
   type(diagnostics_config_t) :: cli, emb
   cli = default_cli_config()
   emb = default_embedded_config()
   ! CLI: stdout on, has a default log file, INFO level.
   @assertTrue(cli%to_stdout)
   @assertEqual(LOGLEVEL_INFO, cli%level)
   @assertTrue(allocated(cli%log_file))
   @assertEqual('swap_swap.log', cli%log_file)
   ! Embedded: stdout OFF (must never pollute host stdout), no default file.
   @assertFalse(emb%to_stdout)
   @assertFalse(allocated(emb%log_file))
end subroutine
```

- [ ] **Step 3: Run to verify it fails**

Run: `pixi run -e test test-pfunit`
Expected: link error — `default_cli_config` / `default_embedded_config` not defined.

- [ ] **Step 4: Implement the env reader + default configs**

In `src/error/diagnostics.f90`, add inside `contains`:

```fortran
   !> Built-in default config for a standalone CLI run.
   pure function default_cli_config() result(cfg)
      type(diagnostics_config_t) :: cfg
      cfg%level      = LOGLEVEL_INFO
      cfg%to_stdout  = .true.
      cfg%to_stderr  = .true.
      cfg%timestamps = .false.
      cfg%log_file   = 'swap_swap.log'
   end function default_cli_config

   !> Built-in default config when embedded (BMI/C-API/XMI): never write the
   !! host's stdout; no log file unless the caller asks (env/C-API).
   pure function default_embedded_config() result(cfg)
      type(diagnostics_config_t) :: cfg
      cfg%level      = LOGLEVEL_INFO
      cfg%to_stdout  = .false.
      cfg%to_stderr  = .true.
      cfg%timestamps = .false.
      ! log_file deliberately left unallocated
   end function default_embedded_config

   !> Read logging overrides from the environment:
   !!   SWAP_LOG_LEVEL   = DEBUG|INFO|WARN|ERROR|NONE
   !!   SWAP_LOG_FILE    = <path>
   !!   SWAP_LOG_STDOUT  = 0|1 (or true/false)
   !! Impure (reads the environment); not in the pure section above.
   subroutine read_env_overrides(ov)
      type(diag_overrides_t), intent(out) :: ov
      character(len=256) :: buf
      integer :: ln, st
      call get_environment_variable('SWAP_LOG_LEVEL', buf, length=ln, status=st)
      if (st == 0 .and. ln > 0) then
         ov%has_level = .true.
         ov%level     = log_level_from_name(buf(:ln))
      end if
      call get_environment_variable('SWAP_LOG_FILE', buf, length=ln, status=st)
      if (st == 0 .and. ln > 0) then
         ov%has_file = .true.
         ov%log_file = buf(:ln)
      end if
      call get_environment_variable('SWAP_LOG_STDOUT', buf, length=ln, status=st)
      if (st == 0 .and. ln > 0) then
         ov%has_stdout = .true.
         ov%to_stdout  = (buf(1:1) == '1' .or. buf(1:1) == 't' .or. buf(1:1) == 'T')
      end if
   end subroutine read_env_overrides
```

- [ ] **Step 5: Run to verify it passes; full gate**

Run: `pixi run -e test test-pfunit` → PASS (count +1).
Run: `pixi run -e test check-fast` → green, byte-identical.

- [ ] **Step 6: Commit**

```bash
git add src/core/swap_log.f90 src/error/diagnostics.f90 tests/unit/error/test_diagnostics.pf
git commit -m "feat(diagnostics): env override reader + CLI/embedded default configs

log_init gains optional to_stderr; diagnostics_mod gains read_env_overrides
(SWAP_LOG_LEVEL/FILE/STDOUT) and default_cli_config/default_embedded_config
(embedded defaults stdout OFF).

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 5: Wire all entry points + C-API setters

This task has no new unit test (entry points are bind(C)/program glue); it is verified by build + `check-fast` byte-identical regression and a manual embedded-stdout check. Each edit resolves a config and calls `log_init`; embedded paths use `default_embedded_config`.

**Files:**
- Modify: `src/driver/swap_main.f90`
- Modify: `src/bindings/swap_bmi_mod.f90`
- Modify: `src/bindings/swap_capi_mod.f90`
- Modify: `src/driver/swap_ensemble_mod.f90`

- [ ] **Step 1: CLI — route `swap_main` through the resolver**

In `src/driver/swap_main.f90`, change the `use swap_log` line and the `log_init` call.

Replace:

```fortran
   use swap_log,        only: log_init, log_close, LOGLEVEL_INFO
```
with:
```fortran
   use swap_log,         only: log_init, log_close
   use diagnostics_mod,  only: diagnostics_config_t, diag_overrides_t, &
                               default_cli_config, read_env_overrides,  &
                               resolve_diagnostics_config
```

Replace:

```fortran
   call log_init(log_level=LOGLEVEL_INFO, log_file='swap_swap.log')
```
with:
```fortran
   block
      type(diagnostics_config_t) :: dcfg
      type(diag_overrides_t)     :: none_ov, env_ov
      call read_env_overrides(env_ov)
      ! precedence: env > default (no TOML/C-API on the CLI path yet)
      dcfg = resolve_diagnostics_config(default_cli_config(), none_ov, env_ov, none_ov)
      if (allocated(dcfg%log_file)) then
         call log_init(log_level=dcfg%level, log_file=dcfg%log_file, &
                       to_stdout=dcfg%to_stdout, to_stderr=dcfg%to_stderr, &
                       timestamps=dcfg%timestamps)
      else
         call log_init(log_level=dcfg%level, to_stdout=dcfg%to_stdout, &
                       to_stderr=dcfg%to_stderr, timestamps=dcfg%timestamps)
      end if
   end block
```

- [ ] **Step 2: BMI — initialise logging in `bmi_initialize` (embedded defaults)**

In `src/bindings/swap_bmi_mod.f90`, add to the module `use` list (after the existing uses):

```fortran
   use swap_log,        only: log_init
   use diagnostics_mod, only: diagnostics_config_t, diag_overrides_t, &
                              default_embedded_config, read_env_overrides, &
                              resolve_diagnostics_config
```

In `bmi_initialize`, insert before `call swap_init(...)`:

```fortran
      block
         type(diagnostics_config_t) :: dcfg
         type(diag_overrides_t)     :: none_ov, env_ov
         call read_env_overrides(env_ov)
         dcfg = resolve_diagnostics_config(default_embedded_config(), none_ov, env_ov, none_ov)
         if (allocated(dcfg%log_file)) then
            call log_init(log_level=dcfg%level, log_file=dcfg%log_file, &
                          to_stdout=dcfg%to_stdout, to_stderr=dcfg%to_stderr, &
                          timestamps=dcfg%timestamps)
         else
            call log_init(log_level=dcfg%level, to_stdout=dcfg%to_stdout, &
                          to_stderr=dcfg%to_stderr, timestamps=dcfg%timestamps)
         end if
      end block
```

- [ ] **Step 3: C-API — add setters + initialise logging in the C-API init paths**

In `src/bindings/swap_capi_mod.f90`, add to the module `use` list:

```fortran
   use swap_log,        only: log_init, log_set_level, log_set_stdout
   use diagnostics_mod, only: diagnostics_config_t, diag_overrides_t, &
                              default_embedded_config, read_env_overrides, &
                              resolve_diagnostics_config
```

Add two bind(C) setters in the `contains` section (after `c_to_f_string`, before the Lifecycle block). These are the **highest-precedence** C-API control surface (a host calls them after initialize):

```fortran
   function swap_set_log_level(level) result(ierr) bind(C, name='swap_set_log_level')
      integer(c_int), value, intent(in) :: level
      integer(c_int)                    :: ierr
      call log_set_level(int(level))
      ierr = 0
   end function swap_set_log_level

   function swap_set_log_stdout(flag) result(ierr) bind(C, name='swap_set_log_stdout')
      integer(c_int), value, intent(in) :: flag
      integer(c_int)                    :: ierr
      call log_set_stdout(flag /= 0)
      ierr = 0
   end function swap_set_log_stdout
```

Add a small private helper (in `contains`) that the C-API init functions call to wire logging with embedded defaults:

```fortran
   subroutine capi_init_logging()
      type(diagnostics_config_t) :: dcfg
      type(diag_overrides_t)     :: none_ov, env_ov
      call read_env_overrides(env_ov)
      dcfg = resolve_diagnostics_config(default_embedded_config(), none_ov, env_ov, none_ov)
      if (allocated(dcfg%log_file)) then
         call log_init(log_level=dcfg%level, log_file=dcfg%log_file, &
                       to_stdout=dcfg%to_stdout, to_stderr=dcfg%to_stderr, &
                       timestamps=dcfg%timestamps)
      else
         call log_init(log_level=dcfg%level, to_stdout=dcfg%to_stdout, &
                       to_stderr=dcfg%to_stderr, timestamps=dcfg%timestamps)
      end if
   end subroutine capi_init_logging
```

Then add `call capi_init_logging()` as the first statement of each C-API initialize entry point — at minimum `swap_initialize_from_toml_string` (seen at the top of the module). Grep for every `bind(C, name='…initialize…')` / `bind(C, name='swap_init…')` function in this file and add the call as their first executable statement:

Run: `grep -n "bind(C, name='" src/bindings/swap_capi_mod.f90 | grep -i 'init'`
For each such initialize function, insert `call capi_init_logging()` immediately after its declarations. (Idempotent: `log_init` safely re-opens; if a host also calls a BMI initialize, the last init wins — acceptable for Phase A.)

- [ ] **Step 4: Ensemble — initialise logging in `ensemble_init` (embedded defaults)**

In `src/driver/swap_ensemble_mod.f90`, add to the `use` list:

```fortran
   use swap_log,        only: log_init
   use diagnostics_mod, only: diagnostics_config_t, diag_overrides_t, &
                              default_embedded_config, read_env_overrides, &
                              resolve_diagnostics_config
```

In `ensemble_init`, insert immediately after `rc = 0` (before `call set_library_mode(.true.)`):

```fortran
      block
         type(diagnostics_config_t) :: dcfg
         type(diag_overrides_t)     :: none_ov, env_ov
         call read_env_overrides(env_ov)
         dcfg = resolve_diagnostics_config(default_embedded_config(), none_ov, env_ov, none_ov)
         if (allocated(dcfg%log_file)) then
            call log_init(log_level=dcfg%level, log_file=dcfg%log_file, &
                          to_stdout=dcfg%to_stdout, to_stderr=dcfg%to_stderr, &
                          timestamps=dcfg%timestamps)
         else
            call log_init(log_level=dcfg%level, to_stdout=dcfg%to_stdout, &
                          to_stderr=dcfg%to_stderr, timestamps=dcfg%timestamps)
         end if
      end block
```

(`set_library_mode` stays for now; the per-instance fatal migration that retires it is Phase D.)

- [ ] **Step 5: Build**

Run: `pixi run build-linux`
Expected: compiles clean (no unused-var / missing-symbol errors).

- [ ] **Step 6: Byte-identical regression gate**

Run: `pixi run -e test check-fast`
Expected: pFUnit green + 4 regression cases byte-identical. The CLI still writes `swap_swap.log` and `Swap normal completion!`; WARN/ERROR now go to stderr (does not affect compared output files).

- [ ] **Step 7: Manual embedded-stdout smoke check**

Confirm an embedded run does not write logging to stdout by default. With the coupling/XMI test available (`pixi run -e test test-xmi` per CLAUDE.md memory), run it and confirm no `INFO/WARN …` logger lines appear on stdout (they may appear on stderr). If `test-xmi` is unavailable in this environment, note it and rely on the `default_embedded_config()` unit test (Task 4) as the guarantee that embedded paths pass `to_stdout=.false.`.

Run: `pixi run -e test test-xmi`
Expected: passes; stdout free of logger lines.

- [ ] **Step 8: Commit**

```bash
git add src/driver/swap_main.f90 src/bindings/swap_bmi_mod.f90 \
        src/bindings/swap_capi_mod.f90 src/driver/swap_ensemble_mod.f90
git commit -m "feat(diagnostics): initialise logging in every entry point

CLI/BMI/C-API/ensemble all resolve a diagnostics config (env > default) and
call log_init; embedded paths default to_stdout=.false. so the Python/MODFLOW
host's stdout is never polluted. Adds bind(C) swap_set_log_level /
swap_set_log_stdout (highest-precedence host control). Byte-identical.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Phase A done-when

- `pixi run -e test check-full` green (run once at the end), byte-identical against `swap420gf`.
- `SWAP_LOG_LEVEL=DEBUG ./swap` raises CLI verbosity with no recompile; `SWAP_LOG_FILE`/`SWAP_LOG_STDOUT` honoured.
- BMI/C-API/XMI initialise logging; embedded default writes nothing to the host's stdout.
- WARN/ERROR appear on stderr; INFO/DEBUG on stdout only when enabled.
- New `test_diagnostics_suite` registered and counted in the `OK (N tests)` total.

## Out of scope (later phases)

- `[logging]` TOML block → start of **Phase B**.
- State-borne `diagnostics_t`, per-instance isolation → **Phase C**.
- Kernel `fatalerr` → `state%diag%fatal` + boundary checks; retiring
  `global_errors`/`fatalerr_collected`/`FatalERR`/`library_mode` → **Phase D**.
- Coverage uplift → **Phase E**.

---

## Self-review notes (author)

- **Spec coverage (Phase A slice of §7):** config layer + precedence (Tasks 1–2) ✓;
  env + C-API surface (Tasks 4–5) ✓; all-entry-point wiring + embedded stdout-off
  (Task 5) ✓; stderr split (Task 3) ✓. TOML block explicitly deferred to B (noted in
  spec). State-borne `diagnostics_t` correctly excluded (Phase C).
- **Placeholder scan:** every code/test step shows complete code; commands have
  expected output. No TBD/TODO.
- **Type consistency:** `diagnostics_config_t` / `diag_overrides_t` field names
  (`level`, `to_stdout`, `to_stderr`, `timestamps`, `log_file`, `has_*`) used
  identically across Tasks 1–5; `log_destination(level,to_stdout,to_stderr,write_stdout,
  write_stderr)`, `log_set_stdout(flag)`, `read_env_overrides(ov)`,
  `resolve_diagnostics_config(base,toml,env,capi)`, `default_cli_config()`,
  `default_embedded_config()` signatures match every call site.
