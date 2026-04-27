# Phase 4a — Infrastructure Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stand up the cross-cutting infrastructure Phase 4 depends on — reviewed logger, error handling module, validation primitives, typed config hierarchy, composable thematic TOML reader, and a temporary config-to-state adapter — then verify end-to-end on the `hupselbrook` case with a full parity test and smoke tests for the other five cases.

**Architecture:** Four-stage pipeline (parse → validate → finalize → adapter) over six section-scoped `*_config_t` types, aggregated into `swap_config_t`. Composite thematic TOML reader: field helpers at the bottom, one section reader per theme, one document dispatcher on top. Errors flow into a single `error_collection_t` that auto-logs through `swap_log` on every append and aborts at one explicit checkpoint. No physics touched.

**Tech Stack:** gfortran 2008, pFUnit 4.15, meson + pixi, toml-f (via meson subproject), `swap_log` (existing), gcovr (for Phase 3 baseline re-run). Spec: `docs/superpowers/specs/2026-04-24-phase-4a-infrastructure-design.md`.

---

## Preamble: context every task needs

**Baseline state** (at Phase 3 exit, commit `74ca107`, tag `rescue/phase-3-coverage`):
- 47 pFUnit tests passing across 14 suites.
- `src/error/`, `src/validation/`, `src/config/`, `src/io/toml/` do not exist.
- `src/io/readswaptoml.f90`, `src/io/readdrainagetoml.f90` exist and are wired into `meson.build`. Their pFUnit tests live at `tests/unit/io/test_readswaptoml.pf` and `tests/unit/io/test_readdrainagetoml.pf`. These files get removed in Task 29.
- `src/core/swap_log.f90` exists (280 lines), complete, untested.
- `tests/swap-cases/` is a git submodule on `main`. No `toml/` subdirectory yet.
- `main` = `development` = `74ca107`. No pushes.

**Branch discipline:**
- Work on `development`. **No per-task feature branches during Phase 4a** — continues the Phase 0–3 convention. Phase 4 proper introduces feature branches; not Phase 4a.
- One commit per task. Subject line `<type>(<scope>): <what>`.
- No pushes to `origin/main` or `origin/development`. Phase 4a stays local.
- Phase 4a closeout fast-forwards `main` to `development` locally and applies tag `rescue/phase-4a-infrastructure`. No push.

**Working directory:** `/home/zawadzkim/Code/swap` for all tasks unless noted.

**Module/file conventions (from `docs/code-style.md`):**
- One module per file. Module name = filename stem + `_mod` suffix (`src/error/error.f90` → `module error_mod`).
- `implicit none` right after `module`. Default `private`; every export is explicit `public :: ...`.
- Kind: prefer `real(kind=real64)` with `use iso_fortran_env, only: real64` for new code.
- Every procedure has explicit `intent(in/inout/out)` on every argument.
- Indent with spaces, not tabs. Lowercase keywords. Short lines preferred (≤132 chars).
- FORD docstrings: `!>` one-line summary immediately above procedure/type; `!!` per-argument descriptions on dummy arg lines.
- No `COMMON`, no `include` of legacy `.fi` headers in new code. No `goto`.
- Always use `end subroutine foo` / `end function foo` / `end module foo`.

**Test conventions** (from Phase 3):
- `tests/unit/<domain>/test_<module>.pf` — one `.pf` per source module.
- `ADD_TEST_SUITE(test_<module>_suite)` in `tests/unit/testSuites.inc`.
- `'<domain>/test_<module>.pf'` in `pf_files` list in `tests/unit/meson.build`.
- Modules outside `test_base_sources` go in `pfunit_extra_sources`, in dependency order.
- **Important funitproc gotcha:** never put inline `! comments` on the same line as `@assertEqual(...)` — `funitproc` chokes on the `!`. Put comments on the line BEFORE.
- Run after every `.pf` change: `pixi run -e test test-pfunit`. Expected line: `OK (N tests)` and meson `Ok: 1`.

**Verification after each task:**
```
pixi run -e test test-pfunit
```
Test count increments by the new suite's tests. Then:
```
pixi run -e test check-fast
```
Must remain green. Wall time budget: <90 s.

**Coverage baseline** (from Phase 3 exit):
- Project line coverage: 51.4% (documented with caveats in `docs/coverage-baseline.md`).
- The Phase 4a closeout (Task 30) re-runs `pixi run -e coverage coverage-report` and updates the baseline doc; no per-domain drop is acceptable as a phase-exit check.

---

## File Structure

### New source files

| Path | Responsibility |
|---|---|
| `src/error/error.f90` | `module error_mod` — `error_t`, `error_collection_t`, typed error-code parameters, `append` auto-logs via `swap_log`, `abort_if_fatal` is the single aborting entry point. |
| `src/error/README.md` | Existing — rewrite to document the new module. |
| `src/validation/validation.f90` | `module validation_mod` — stateless primitive checks: range, enum, not-empty, ordered pair, non-negative. Each appends errors; no side effects. |
| `src/validation/README.md` | New — Responsibility / Public interface / Dependencies. |
| `src/config/general_config.f90` | `module general_config_mod` — `general_config_t` + `validate` + `finalize`. Holds `[general]` schema: project, paths, screen/error switches. |
| `src/config/simulation_config.f90` | `simulation_config_t` — `[simulation]`: start/end dates, nprintday, output-timing switches. |
| `src/config/meteorology_config.f90` | `meteorology_config_t` — `[meteorology]`: lat/alt/altw, ET switches, Angstrom coefficients, interception, rain, detailed meteo. |
| `src/config/drainage_config.f90` | `drainage_config_t` — `[drainage]`: swdra, dramet, basic/extended drainage tables. |
| `src/config/soil_config.f90` | `soil_config_t` — `[soil]`: profile, hydraulic parameters, initial conditions. |
| `src/config/crop_config.f90` | `crop_config_t` — `[crop]`: crop rotation list, per-crop file references. |
| `src/config/swap_config.f90` | `swap_config_t` aggregator — composes the six section types, top-level `validate` / `finalize`. |
| `src/config/README.md` | New — Responsibility / Public interface / Dependencies. |
| `src/io/toml/toml_field_helpers.f90` | Primitives: `get_required_int`, `get_required_real`, `get_required_string`, `get_optional_int_with_default`, `get_optional_real_with_default`, `get_optional_string_with_default`, `get_table`, `get_array_of_tables`, `parse_date_to_days1900`. All null-ptr safe; all append PARSE errors. |
| `src/io/toml/read_general_toml.f90` | `read_general_toml(doc, config%general, errors)` — populates `general_config_t`. |
| `src/io/toml/read_simulation_toml.f90` | `read_simulation_toml(doc, config%simulation, errors)`. |
| `src/io/toml/read_meteorology_toml.f90` | `read_meteorology_toml(doc, config%meteo, errors)`. |
| `src/io/toml/read_drainage_toml.f90` | `read_drainage_toml(doc, config%drain, errors)`. |
| `src/io/toml/read_soil_toml.f90` | `read_soil_toml(doc, config%soil, errors)`. |
| `src/io/toml/read_crop_toml.f90` | `read_crop_toml(doc, config%crop, errors)`. |
| `src/io/toml/load_swap_config.f90` | Top-level entrypoint `load_swap_config(path, config, errors)`. Loads file, extracts root table, dispatches to each section reader. |
| `src/io/toml/write_swap_config.f90` | Emit: `write_swap_config(config, path, errors)` — inverse of load. |
| `src/io/toml/README.md` | New — Responsibility / Public interface / Dependencies. |
| `src/core/config_to_state.f90` | `module config_to_state_mod` — temporary adapter `config_to_state(config, state)`. Field-copy only; no business logic. Header comment declares retirement at Phase 4. |

### New test files

One `.pf` per new source module, under `tests/unit/<domain>/test_<module>.pf`. See per-task wiring.

Plus integration tests:
- `tests/unit/io/toml/test_hupselbrook_parity.pf` — Task 27
- `tests/unit/io/toml/test_hupselbrook_roundtrip.pf` — Task 28
- `tests/unit/io/toml/test_all_cases_smoke.pf` — Task 27

### New TOML project files

Under `tests/swap-cases/toml/` (inside the `swap-cases` git submodule, on its `main` branch):

| Path | Content |
|---|---|
| `1.hupselbrook/swap.toml` | Full TOML equivalent of `tests/swap-cases/1.hupselbrook/swap_linux.swp.template` |
| `1.hupselbrook/swap.dra.toml` | Drainage (or folded into `swap.toml` — decided in Task 27) |
| `1.hupselbrook/*.crp.toml` | Crop files, one per crop referenced |
| `2.grassgrowth/*.toml` | Grass-focused case, ditto |
| `3.macroporeflow/*.toml` | Macropore case |
| `4.oxygenstress/*.toml` | Oxygen-stress case |
| `5.salinitystress/*.toml` | Salinity-stress case |
| `6.surfacewater/*.toml` | Surface-water case |

Each case dir contains a `README.md` noting it is the TOML-native port of `tests/swap-cases/N.<name>/`, hand-authored during Phase 4a.

### Files modified

| Path | Change |
|---|---|
| `meson.build` | Add new source files (`src/error/error.f90`, `src/validation/validation.f90`, six `src/config/*.f90`, eight `src/io/toml/*.f90`, `src/core/config_to_state.f90`) to the `sources` list. Order matters: error before validation before config before io/toml. |
| `tests/unit/meson.build` | Extend `pf_files` and `pfunit_extra_sources` as suites are added. |
| `tests/unit/testSuites.inc` | One `ADD_TEST_SUITE` per suite. |
| `src/core/README.md` | Mention `config_to_state.f90` and `swap_log.f90` after audit. |
| `src/io/README.md` | Mention `src/io/toml/` subdir (new architecture alongside legacy readers). |
| `docs/index.md` | Link new pages. |
| `docs/configuration-schema.md` | Update to describe `*_config_t` hierarchy, not just TOML keys. |
| `.gitmodules` (if needed) | Unchanged — `swap-cases` submodule ref updates via normal commit. |

### Files deleted (Task 29)

| Path | Why |
|---|---|
| `src/io/readswaptoml.f90` | Replaced by composite reader under `src/io/toml/`. |
| `src/io/readdrainagetoml.f90` | Same. |
| `tests/unit/io/test_readswaptoml.pf` | Test subject deleted. |
| `tests/unit/io/test_readdrainagetoml.pf` | Same. |
| `tests/unit/io/fixtures/minimal_swap.toml` | Replaced by hupselbrook TOML + per-section minimal fixtures under `tests/unit/io/toml/fixtures/`. |
| `tests/unit/io/fixtures/minimal_drainage.toml` | Same. |
| `tests/unit/io/fixtures/malformed_drainage.toml` | The new field helpers have their own malformed fixtures. |

---

## Task 1: Logger audit and `docs/logging.md`

**Files:**
- Create: `docs/logging.md`
- Modify: `docs/index.md`, `src/core/README.md`

- [ ] **Step 1: Audit `src/core/swap_log.f90`**

Read `src/core/swap_log.f90` end-to-end. Check for:

- Thread safety: are any module-level variables mutated without synchronisation? (Answer today: yes — `current_level`, `log_unit`, `log_to_file`, `log_to_stdout`, `log_initialized`. No atomics, no locks. Acceptable for serial gfortran; document the constraint.)
- `to_str` overloads available: `int_to_str`, `real_to_str`, `real8_to_str`, `logical_to_str`. Missing: none required for Phase 4a; if the error module needs another overload it's added in Task 3.
- `log_init` / `log_close` lifecycle: called from `swap_main.f90`? Check.
- Behaviour when `log_init` is NOT called (i.e., tests that use the logger without init): does `log_message` still work? (Currently: yes — `current_level = LOGLEVEL_INFO` default, `log_to_stdout = .true.` default; the level and stdout guards make it safe.)
- Buffer flushing: `log_message` calls `flush(log_unit)` after file writes. Good.

Record findings as notes; do NOT modify `src/core/swap_log.f90` in this task.

- [ ] **Step 2: Create `docs/logging.md`**

```markdown
---
title: Logging
author: SWAP modernization team
---

# Logging

## Overview

`swap_log` (at `src/core/swap_log.f90`) is the single logging facility for
the modernised SWAP tree. Every module that needs to emit diagnostics uses
it. The error module (`src/error/error.f90`, Phase 4a) routes every
appended error through `log_error` automatically, so application code never
needs to call logger and error separately.

## Levels

Four levels, lowest to highest severity:

| Constant         | Integer value | Typical use                                   |
|---|---|---|
| `LOGLEVEL_DEBUG` | 10            | High-volume developer traces; off by default. |
| `LOGLEVEL_INFO`  | 20            | Progress messages; default threshold.         |
| `LOGLEVEL_WARN`  | 30            | Unexpected but recoverable conditions.        |
| `LOGLEVEL_ERROR` | 40            | Every appended `error_t` logs here.           |

Plus `LOGLEVEL_NONE` (100) to silence the logger entirely.

## Usage

Initialise once at program start (typically in `swap_main.f90` or its
equivalent in test drivers):

    use swap_log
    call log_init(log_level=LOGLEVEL_INFO, log_file='swap.log')

Then anywhere:

    call log_info('reader',  'Loading TOML file: ' // trim(path))
    call log_warn('finalize', 'Derived value clipped to bounds')
    call log_error('validate', 'drainage.dramet out of range')

Close on shutdown:

    call log_close

`log_debug`, `log_info`, `log_warn`, `log_error` take `(context, message)`.
`context` is a short module-or-subsystem tag (e.g., `"reader"`, `"validate"`).
Message is any string; use `to_str(n)` from the same module to stringify
integers, reals, logicals.

## Output shape

    2026-04-24 09:42:17.312 INFO  reader: Loading TOML file: swap.toml

Format is fixed: `TIMESTAMP LEVEL CONTEXT: MESSAGE`. Timestamps can be
suppressed via `log_init(..., timestamps=.false.)`.

## Interaction with the error module

When `src/error/error.f90` (Phase 4a) appends to an `error_collection_t`,
`swap_log%log_error(context, message)` is called automatically. Do not
log errors manually in addition — the append call handles it.

## Thread safety

None. `swap_log` mutates module-level state (`current_level`, `log_unit`,
file-handle state) without synchronisation. Safe for serial gfortran; when
multicore work lands in a future phase, the logger either becomes
thread-local or gains synchronisation. Do not rely on logger behaviour
under concurrent calls today.

## Testing

pFUnit unit tests for the logger are added as part of Phase 4a Task 2
(they land alongside the error-module tests, since those exercise the
logger-sink path).
```

- [ ] **Step 3: Link from `docs/index.md` and `src/core/README.md`**

Open `docs/index.md`, find the existing docs list (architecture, state-management, etc.), add a bullet for `logging.md` in the same format.

Open `src/core/README.md`. In the Public interface section add a line referencing `swap_log`. Mention the new `config_to_state.f90` will be added later in Phase 4a (Task 26).

- [ ] **Step 4: Verify FORD builds**

```
pixi run -e docs docs-build
```

Expected: exit 0, no errors.

- [ ] **Step 5: Commit**

```
git add docs/logging.md docs/index.md src/core/README.md
git commit -m "docs(logging): audit and document swap_log"
```

---

## Task 2: `error_mod` scaffolding (types, no methods)

**Files:**
- Create: `src/error/error.f90`
- Modify: `meson.build` (add to `sources`), `src/error/README.md`

- [ ] **Step 1: Write the module skeleton**

Create `src/error/error.f90`:

```fortran
!> Error handling for SWAP.
!!
!! Defines a typed error payload (`error_t`) and an accumulating
!! collection (`error_collection_t`). Every fallible procedure in the
!! new infrastructure takes `errors` as `intent(inout)` and appends
!! failures without aborting. A single abort point (`abort_if_fatal`)
!! turns accumulated fatal errors into `error stop`.
module error_mod
   use iso_fortran_env, only: real64, error_unit
   use swap_log, only: log_error
   implicit none
   private

   ! Error code parameters. Stable values; tests assert on these.
   integer, parameter, public :: ERR_NONE                     = 0
   integer, parameter, public :: ERR_IO_READ_FAILED           = 100
   integer, parameter, public :: ERR_IO_WRITE_FAILED          = 101
   integer, parameter, public :: ERR_PARSE_MALFORMED_TOML     = 200
   integer, parameter, public :: ERR_PARSE_TYPE_MISMATCH      = 201
   integer, parameter, public :: ERR_PARSE_MISSING_REQUIRED   = 202
   integer, parameter, public :: ERR_VALIDATION_OUT_OF_RANGE  = 300
   integer, parameter, public :: ERR_VALIDATION_ENUM          = 301
   integer, parameter, public :: ERR_VALIDATION_CROSS_FIELD   = 302
   integer, parameter, public :: ERR_VALIDATION_CROSS_SECTION = 303
   integer, parameter, public :: ERR_FINALIZE_DERIVATION      = 400
   integer, parameter, public :: ERR_ADAPTER_UNSUPPORTED      = 500

   type, public :: error_t
      integer                       :: code    = ERR_NONE
      character(len=:), allocatable :: message
      character(len=:), allocatable :: context
      logical                       :: is_fatal = .false.
   end type error_t

   type, public :: error_collection_t
      type(error_t), allocatable :: items(:)
   contains
      procedure :: append         => error_collection_append
      procedure :: has_errors     => error_collection_has_errors
      procedure :: has_fatals     => error_collection_has_fatals
      procedure :: count          => error_collection_count
      procedure :: summary        => error_collection_summary
      procedure :: abort_if_fatal => error_collection_abort_if_fatal
      procedure :: clear          => error_collection_clear
   end type error_collection_t

contains

   ! Method bodies are added in Tasks 3-5.

end module error_mod
```

At this point the module compiles as a type-only scaffold. Tasks 3–5 add the method bodies one concern at a time.

- [ ] **Step 2: Add to meson.build**

Open `meson.build`. In the `sources = [ ... ]` list find the `# Core` section (starts at line 51). Add `'src/error/error.f90'` immediately AFTER `'src/core/swap_log.f90'` so the module's `use swap_log` dependency is satisfied. Order in the list hints meson's compile order but gfortran resolves `use` statements automatically — still, keep the order readable.

Result should look like:

```meson
    'src/core/swap_log.f90',
    'src/error/error.f90',
    'src/atmosphere/atmosphere_state.f90',
    ...
```

- [ ] **Step 3: Rewrite `src/error/README.md`**

Replace the existing placeholder content with:

```markdown
# src/error

## Responsibility

Typed error payload + accumulating collection used across the new
configuration pipeline (Phase 4a). Every fallible procedure takes
`errors` as `intent(inout)` and appends failures; one checkpoint
(`abort_if_fatal`) aborts if any fatal error was collected.

See `docs/error-handling.md` for the usage guide and
`docs/adr/0008-error-collection-over-fatalerr.md` for the decision record.

## Public interface

- `error_t` — single error payload (`code`, `message`, `context`, `is_fatal`).
- `error_collection_t` — accumulator with type-bound procedures:
  `append`, `has_errors`, `has_fatals`, `count`, `summary`,
  `abort_if_fatal`, `clear`.
- Error code constants: `ERR_NONE`, `ERR_IO_READ_FAILED`, `ERR_IO_WRITE_FAILED`,
  `ERR_PARSE_MALFORMED_TOML`, `ERR_PARSE_TYPE_MISMATCH`,
  `ERR_PARSE_MISSING_REQUIRED`, `ERR_VALIDATION_OUT_OF_RANGE`,
  `ERR_VALIDATION_ENUM`, `ERR_VALIDATION_CROSS_FIELD`,
  `ERR_VALIDATION_CROSS_SECTION`, `ERR_FINALIZE_DERIVATION`,
  `ERR_ADAPTER_UNSUPPORTED`.

## Dependencies

`swap_log` (for auto-logging on append).
```

- [ ] **Step 4: Verify build**

```
pixi run -e test build-linux
```

Expected: compile succeeds. The module has type definitions only — no method bodies — so the `contains` block is empty. gfortran tolerates this.

- [ ] **Step 5: Commit**

```
git add src/error/error.f90 src/error/README.md meson.build
git commit -m "feat(error): scaffold error_mod with types and error-code parameters"
```

---

## Task 3: `error_collection_t%append` + logger integration + tests

**Files:**
- Modify: `src/error/error.f90`
- Create: `tests/unit/error/test_error.pf`
- Modify: `tests/unit/testSuites.inc`, `tests/unit/meson.build`

- [ ] **Step 1: Add the `append` method to `error_mod`**

In `src/error/error.f90`, inside the `contains` block, add:

```fortran
   !> Append an error to the collection and auto-log through swap_log.
   subroutine error_collection_append(self, code, message, context, is_fatal)
      class(error_collection_t), intent(inout) :: self
      integer,                   intent(in)    :: code
      character(len=*),          intent(in)    :: message
      character(len=*),          intent(in)    :: context
      logical, optional,         intent(in)    :: is_fatal

      type(error_t), allocatable :: tmp(:)
      type(error_t)              :: item
      integer                    :: n

      item%code     = code
      item%message  = message
      item%context  = context
      if (present(is_fatal)) then
         item%is_fatal = is_fatal
      else
         item%is_fatal = .true.
      end if

      if (.not. allocated(self%items)) then
         allocate(self%items(1))
         self%items(1) = item
      else
         n = size(self%items)
         allocate(tmp(n + 1))
         tmp(1:n)   = self%items
         tmp(n + 1) = item
         call move_alloc(from=tmp, to=self%items)
      end if

      call log_error(context, message)
   end subroutine error_collection_append
```

The default `is_fatal` is `.true.` because Phase 4a treats all parse and validation errors as fatal; finalize explicitly passes `is_fatal=.false.` for warnings.

- [ ] **Step 2: Write the test suite**

Create `tests/unit/error/test_error.pf`:

```fortran
! Tests for error_mod%append. Uses append + direct inspection of
! self%items since that's the only method available at this point in
! the plan; has_errors/has_fatals/etc. arrive in Task 4.

@test
subroutine test_append_adds_one_item()
   use funit
   use error_mod
   type(error_collection_t) :: errors

   call errors%append(ERR_PARSE_MISSING_REQUIRED, &
                      "field 'foo' is required", "drainage.foo")

   @assertTrue(allocated(errors%items))
   @assertEqual(1, size(errors%items))
   @assertEqual(ERR_PARSE_MISSING_REQUIRED, errors%items(1)%code)
   @assertEqual("field 'foo' is required", errors%items(1)%message)
   @assertEqual("drainage.foo", errors%items(1)%context)
   @assertTrue(errors%items(1)%is_fatal)
end subroutine

@test
subroutine test_append_preserves_order()
   use funit
   use error_mod
   type(error_collection_t) :: errors

   call errors%append(ERR_PARSE_MISSING_REQUIRED, "first",  "a")
   call errors%append(ERR_VALIDATION_OUT_OF_RANGE, "second", "b")
   call errors%append(ERR_PARSE_TYPE_MISMATCH,     "third",  "c")

   @assertEqual(3, size(errors%items))
   @assertEqual("first",  errors%items(1)%message)
   @assertEqual("second", errors%items(2)%message)
   @assertEqual("third",  errors%items(3)%message)
end subroutine

@test
subroutine test_append_respects_is_fatal_argument()
   use funit
   use error_mod
   type(error_collection_t) :: errors

   call errors%append(ERR_FINALIZE_DERIVATION, "soft warning", "x", is_fatal=.false.)

   @assertFalse(errors%items(1)%is_fatal)
end subroutine

@test
subroutine test_append_defaults_to_fatal()
   use funit
   use error_mod
   type(error_collection_t) :: errors

   call errors%append(ERR_PARSE_MALFORMED_TOML, "bad token", "line 3")

   @assertTrue(errors%items(1)%is_fatal)
end subroutine
```

- [ ] **Step 3: Register the suite**

Append to `tests/unit/testSuites.inc`:
```fortran
ADD_TEST_SUITE(test_error_suite)
```

Extend `pf_files` in `tests/unit/meson.build`:
```meson
        'error/test_error.pf',
```

Extend `pfunit_extra_sources`:
```meson
        '../../src/error/error.f90',
```

The error module depends on `swap_log`, which is in `test_base_sources`. Order: list `error.f90` after the existing extras; meson resolves `use` chains.

- [ ] **Step 4: Run**

```
pixi run -e test test-pfunit
```

Expected: `OK (51 tests)` (47 prior + 4 new). Exit 0.

- [ ] **Step 5: Commit**

```
git add src/error/error.f90 tests/unit/error/test_error.pf \
        tests/unit/testSuites.inc tests/unit/meson.build
git commit -m "feat(error): error_collection_t%append with auto-logging"
```

---

## Task 4: Query methods (`has_errors`, `has_fatals`, `count`, `clear`)

**Files:**
- Modify: `src/error/error.f90`, `tests/unit/error/test_error.pf`

- [ ] **Step 1: Add method bodies**

Append inside the `contains` block of `src/error/error.f90`:

```fortran
   !> .true. if any items have been appended.
   pure function error_collection_has_errors(self) result(yes)
      class(error_collection_t), intent(in) :: self
      logical :: yes
      yes = allocated(self%items) .and. size_safe(self) > 0
   end function error_collection_has_errors

   !> .true. if any appended item has is_fatal=.true.
   pure function error_collection_has_fatals(self) result(yes)
      class(error_collection_t), intent(in) :: self
      logical :: yes
      integer :: i
      yes = .false.
      if (.not. allocated(self%items)) return
      do i = 1, size(self%items)
         if (self%items(i)%is_fatal) then
            yes = .true.
            return
         end if
      end do
   end function error_collection_has_fatals

   !> Number of errors appended so far.
   pure function error_collection_count(self) result(n)
      class(error_collection_t), intent(in) :: self
      integer :: n
      if (allocated(self%items)) then
         n = size(self%items)
      else
         n = 0
      end if
   end function error_collection_count

   !> Reset the collection. Use between tests or logical phases.
   subroutine error_collection_clear(self)
      class(error_collection_t), intent(inout) :: self
      if (allocated(self%items)) deallocate(self%items)
   end subroutine error_collection_clear

   !> Helper used by has_errors (keeps the function pure).
   pure function size_safe(self) result(n)
      class(error_collection_t), intent(in) :: self
      integer :: n
      if (allocated(self%items)) then
         n = size(self%items)
      else
         n = 0
      end if
   end function size_safe
```

- [ ] **Step 2: Extend the test suite**

Append to `tests/unit/error/test_error.pf`:

```fortran
@test
subroutine test_has_errors_empty_is_false()
   use funit
   use error_mod
   type(error_collection_t) :: errors
   @assertFalse(errors%has_errors())
end subroutine

@test
subroutine test_has_errors_after_append_is_true()
   use funit
   use error_mod
   type(error_collection_t) :: errors
   call errors%append(ERR_PARSE_TYPE_MISMATCH, "m", "c")
   @assertTrue(errors%has_errors())
end subroutine

@test
subroutine test_has_fatals_only_counts_fatal_items()
   use funit
   use error_mod
   type(error_collection_t) :: errors
   call errors%append(ERR_FINALIZE_DERIVATION, "soft", "x", is_fatal=.false.)
   @assertTrue(errors%has_errors())
   @assertFalse(errors%has_fatals())
   call errors%append(ERR_PARSE_MISSING_REQUIRED, "hard", "y")
   @assertTrue(errors%has_fatals())
end subroutine

@test
subroutine test_count_matches_appends()
   use funit
   use error_mod
   type(error_collection_t) :: errors
   @assertEqual(0, errors%count())
   call errors%append(ERR_PARSE_TYPE_MISMATCH, "a", "c")
   @assertEqual(1, errors%count())
   call errors%append(ERR_PARSE_TYPE_MISMATCH, "b", "c")
   @assertEqual(2, errors%count())
end subroutine

@test
subroutine test_clear_resets()
   use funit
   use error_mod
   type(error_collection_t) :: errors
   call errors%append(ERR_PARSE_TYPE_MISMATCH, "a", "c")
   call errors%clear()
   @assertFalse(errors%has_errors())
   @assertEqual(0, errors%count())
end subroutine
```

- [ ] **Step 3: Run**

```
pixi run -e test test-pfunit
```

Expected: `OK (56 tests)` (51 prior + 5 new).

- [ ] **Step 4: Commit**

```
git add src/error/error.f90 tests/unit/error/test_error.pf
git commit -m "feat(error): query methods (has_errors, has_fatals, count, clear)"
```

---

## Task 5: `summary` and `abort_if_fatal`

**Files:**
- Modify: `src/error/error.f90`, `tests/unit/error/test_error.pf`

- [ ] **Step 1: Add method bodies**

Append to the `contains` block of `src/error/error.f90`:

```fortran
   !> Multi-line human-readable report of all collected errors.
   function error_collection_summary(self) result(text)
      class(error_collection_t), intent(in) :: self
      character(len=:), allocatable :: text
      character(len=32)             :: code_str
      integer                       :: i

      if (.not. allocated(self%items) .or. size(self%items) == 0) then
         text = "No errors."
         return
      end if

      text = ""
      do i = 1, size(self%items)
         write(code_str, '(I0)') self%items(i)%code
         text = text // "[" // trim(adjustl(code_str)) // "]"
         if (self%items(i)%is_fatal) then
            text = text // " FATAL "
         else
            text = text // " WARN  "
         end if
         text = text // self%items(i)%context // ": " // self%items(i)%message // new_line('a')
      end do
   end function error_collection_summary

   !> Write summary to stderr and `error stop` if any fatal errors exist.
   !! Returns normally if no fatals.
   subroutine error_collection_abort_if_fatal(self)
      class(error_collection_t), intent(in) :: self
      if (.not. self%has_fatals()) return
      write(error_unit, '(A)') self%summary()
      error stop "fatal error(s) in swap input pipeline"
   end subroutine error_collection_abort_if_fatal
```

- [ ] **Step 2: Extend the test suite**

Append to `tests/unit/error/test_error.pf`:

```fortran
@test
subroutine test_summary_empty()
   use funit
   use error_mod
   type(error_collection_t) :: errors
   @assertEqual("No errors.", errors%summary())
end subroutine

@test
subroutine test_summary_includes_code_context_message()
   use funit
   use error_mod
   type(error_collection_t) :: errors
   character(len=:), allocatable :: text

   call errors%append(ERR_PARSE_MISSING_REQUIRED, "need foo", "drainage.foo")
   text = errors%summary()

   @assertTrue(index(text, "202") > 0)
   @assertTrue(index(text, "FATAL") > 0)
   @assertTrue(index(text, "drainage.foo") > 0)
   @assertTrue(index(text, "need foo") > 0)
end subroutine

@test
subroutine test_summary_marks_nonfatal_as_warn()
   use funit
   use error_mod
   type(error_collection_t) :: errors
   character(len=:), allocatable :: text

   call errors%append(ERR_FINALIZE_DERIVATION, "derivation warning", "x", is_fatal=.false.)
   text = errors%summary()

   @assertTrue(index(text, "WARN") > 0)
   @assertFalse(index(text, "FATAL") > 0)
end subroutine

@test
subroutine test_abort_if_fatal_returns_when_no_fatals()
   use funit
   use error_mod
   type(error_collection_t) :: errors
   call errors%append(ERR_FINALIZE_DERIVATION, "ok", "x", is_fatal=.false.)
   call errors%abort_if_fatal()
   @assertTrue(.true.)
end subroutine
```

`abort_if_fatal` is not directly testable in pFUnit when it halts (can't catch `error stop`). The return-when-no-fatals path is tested above; the halt path is exercised transitively by the parity test in Task 27 if any assertion fails there.

- [ ] **Step 3: Run**

```
pixi run -e test test-pfunit
```

Expected: `OK (60 tests)` (56 prior + 4 new).

- [ ] **Step 4: Commit**

```
git add src/error/error.f90 tests/unit/error/test_error.pf
git commit -m "feat(error): summary and abort_if_fatal"
```

---

## Task 6: `validation_mod` — primitive checks

**Files:**
- Create: `src/validation/validation.f90`, `src/validation/README.md`
- Create: `tests/unit/validation/test_validation.pf`
- Modify: `meson.build`, `tests/unit/testSuites.inc`, `tests/unit/meson.build`

- [ ] **Step 1: Write `validation_mod`**

Create `src/validation/validation.f90`:

```fortran
!> Shared validator primitives used by per-section config validators.
!!
!! Each subroutine checks one invariant and appends an error to the
!! passed `error_collection_t` on failure. No state. No side effects
!! other than the append (which auto-logs). Pass `context` so the
!! error identifies the offending field.
module validation_mod
   use iso_fortran_env, only: real64
   use error_mod, only: error_collection_t,          &
                        ERR_VALIDATION_OUT_OF_RANGE, &
                        ERR_VALIDATION_ENUM,         &
                        ERR_VALIDATION_CROSS_FIELD
   implicit none
   private

   public :: check_int_range
   public :: check_real_range
   public :: check_int_enum
   public :: check_not_empty
   public :: check_nonnegative_real
   public :: check_ordered_pair

contains

   !> Verify low <= value <= high.
   subroutine check_int_range(value, low, high, context, errors)
      integer,                   intent(in)    :: value, low, high
      character(len=*),          intent(in)    :: context
      type(error_collection_t),  intent(inout) :: errors
      character(len=64) :: msg
      if (value < low .or. value > high) then
         write(msg, '("value ",I0," outside [",I0,",",I0,"]")') value, low, high
         call errors%append(ERR_VALIDATION_OUT_OF_RANGE, trim(msg), context)
      end if
   end subroutine check_int_range

   !> Verify low <= value <= high (real64).
   subroutine check_real_range(value, low, high, context, errors)
      real(real64),              intent(in)    :: value, low, high
      character(len=*),          intent(in)    :: context
      type(error_collection_t),  intent(inout) :: errors
      character(len=128) :: msg
      if (value < low .or. value > high) then
         write(msg, '("value ",ES12.5," outside [",ES12.5,",",ES12.5,"]")') value, low, high
         call errors%append(ERR_VALIDATION_OUT_OF_RANGE, trim(msg), context)
      end if
   end subroutine check_real_range

   !> Verify value is one of the allowed integers.
   subroutine check_int_enum(value, allowed, context, errors)
      integer,                   intent(in)    :: value
      integer,                   intent(in)    :: allowed(:)
      character(len=*),          intent(in)    :: context
      type(error_collection_t),  intent(inout) :: errors
      character(len=128) :: msg
      integer :: i
      do i = 1, size(allowed)
         if (value == allowed(i)) return
      end do
      write(msg, '("value ",I0," not in allowed set")') value
      call errors%append(ERR_VALIDATION_ENUM, trim(msg), context)
   end subroutine check_int_enum

   !> Verify the string (trimmed) is non-empty.
   subroutine check_not_empty(value, context, errors)
      character(len=*),          intent(in)    :: value
      character(len=*),          intent(in)    :: context
      type(error_collection_t),  intent(inout) :: errors
      if (len_trim(value) == 0) then
         call errors%append(ERR_VALIDATION_OUT_OF_RANGE, "empty string", context)
      end if
   end subroutine check_not_empty

   !> Verify value >= 0.
   subroutine check_nonnegative_real(value, context, errors)
      real(real64),              intent(in)    :: value
      character(len=*),          intent(in)    :: context
      type(error_collection_t),  intent(inout) :: errors
      character(len=64) :: msg
      if (value < 0.0_real64) then
         write(msg, '("value ",ES12.5," is negative")') value
         call errors%append(ERR_VALIDATION_OUT_OF_RANGE, trim(msg), context)
      end if
   end subroutine check_nonnegative_real

   !> Verify low_val <= high_val (cross-field real64 check).
   subroutine check_ordered_pair(low_val, high_val, low_name, high_name, context, errors)
      real(real64),              intent(in)    :: low_val, high_val
      character(len=*),          intent(in)    :: low_name, high_name
      character(len=*),          intent(in)    :: context
      type(error_collection_t),  intent(inout) :: errors
      character(len=256) :: msg
      if (low_val > high_val) then
         write(msg, '(A," (",ES12.5,") exceeds ",A," (",ES12.5,")")') &
              low_name, low_val, high_name, high_val
         call errors%append(ERR_VALIDATION_CROSS_FIELD, trim(msg), context)
      end if
   end subroutine check_ordered_pair

end module validation_mod
```

- [ ] **Step 2: Create `src/validation/README.md`**

```markdown
# src/validation

## Responsibility

Stateless validator primitives consumed by per-section config
validators. Each primitive checks one invariant and appends an
error on failure via an `error_collection_t` passed in. No shared
state, no side effects beyond the append.

## Public interface

All in `validation_mod`:

- `check_int_range(value, low, high, context, errors)`
- `check_real_range(value, low, high, context, errors)`
- `check_int_enum(value, allowed, context, errors)`
- `check_not_empty(value, context, errors)`
- `check_nonnegative_real(value, context, errors)`
- `check_ordered_pair(low_val, high_val, low_name, high_name, context, errors)`

## Dependencies

`error_mod`.
```

- [ ] **Step 3: Add to meson.build**

In `meson.build` `sources` list, add `'src/validation/validation.f90'` after `'src/error/error.f90'`.

- [ ] **Step 4: Write the test suite**

Create `tests/unit/validation/test_validation.pf`:

```fortran
! Tests for validation_mod primitives.

@test
subroutine test_check_int_range_inside_passes()
   use funit
   use error_mod
   use validation_mod
   type(error_collection_t) :: errors
   call check_int_range(5, 1, 10, "x", errors)
   @assertFalse(errors%has_errors())
end subroutine

@test
subroutine test_check_int_range_below_appends_error()
   use funit
   use error_mod
   use validation_mod
   type(error_collection_t) :: errors
   call check_int_range(0, 1, 10, "x", errors)
   @assertTrue(errors%has_errors())
   @assertEqual(ERR_VALIDATION_OUT_OF_RANGE, errors%items(1)%code)
end subroutine

@test
subroutine test_check_int_range_above_appends_error()
   use funit
   use error_mod
   use validation_mod
   type(error_collection_t) :: errors
   call check_int_range(11, 1, 10, "x", errors)
   @assertTrue(errors%has_errors())
end subroutine

@test
subroutine test_check_real_range_edges_pass()
   use funit
   use error_mod
   use validation_mod
   type(error_collection_t) :: errors
   call check_real_range(1.0d0,  1.0d0, 10.0d0, "x", errors)
   call check_real_range(10.0d0, 1.0d0, 10.0d0, "x", errors)
   @assertFalse(errors%has_errors())
end subroutine

@test
subroutine test_check_int_enum_member_passes()
   use funit
   use error_mod
   use validation_mod
   type(error_collection_t) :: errors
   integer :: allowed(3)
   allowed = [1, 2, 3]
   call check_int_enum(2, allowed, "x", errors)
   @assertFalse(errors%has_errors())
end subroutine

@test
subroutine test_check_int_enum_nonmember_appends()
   use funit
   use error_mod
   use validation_mod
   type(error_collection_t) :: errors
   integer :: allowed(3)
   allowed = [1, 2, 3]
   call check_int_enum(5, allowed, "x", errors)
   @assertTrue(errors%has_errors())
   @assertEqual(ERR_VALIDATION_ENUM, errors%items(1)%code)
end subroutine

@test
subroutine test_check_not_empty_nonempty_passes()
   use funit
   use error_mod
   use validation_mod
   type(error_collection_t) :: errors
   call check_not_empty("hello", "x", errors)
   @assertFalse(errors%has_errors())
end subroutine

@test
subroutine test_check_not_empty_empty_appends()
   use funit
   use error_mod
   use validation_mod
   type(error_collection_t) :: errors
   call check_not_empty("   ", "x", errors)
   @assertTrue(errors%has_errors())
end subroutine

@test
subroutine test_check_nonnegative_real_positive_passes()
   use funit
   use error_mod
   use validation_mod
   type(error_collection_t) :: errors
   call check_nonnegative_real(5.0d0, "x", errors)
   @assertFalse(errors%has_errors())
end subroutine

@test
subroutine test_check_nonnegative_real_negative_appends()
   use funit
   use error_mod
   use validation_mod
   type(error_collection_t) :: errors
   call check_nonnegative_real(-0.01d0, "x", errors)
   @assertTrue(errors%has_errors())
end subroutine

@test
subroutine test_check_ordered_pair_ordered_passes()
   use funit
   use error_mod
   use validation_mod
   type(error_collection_t) :: errors
   call check_ordered_pair(1.0d0, 2.0d0, "low", "high", "x", errors)
   @assertFalse(errors%has_errors())
end subroutine

@test
subroutine test_check_ordered_pair_inverted_appends()
   use funit
   use error_mod
   use validation_mod
   type(error_collection_t) :: errors
   call check_ordered_pair(2.0d0, 1.0d0, "low", "high", "x", errors)
   @assertTrue(errors%has_errors())
   @assertEqual(ERR_VALIDATION_CROSS_FIELD, errors%items(1)%code)
end subroutine
```

- [ ] **Step 5: Register suite**

Append to `tests/unit/testSuites.inc`: `ADD_TEST_SUITE(test_validation_suite)`.

Extend `pf_files`: `'validation/test_validation.pf'`.

Extend `pfunit_extra_sources`: `'../../src/validation/validation.f90'`.

- [ ] **Step 6: Run**

```
pixi run -e test test-pfunit
```

Expected: `OK (72 tests)` (60 prior + 12 new).

- [ ] **Step 7: Commit**

```
git add src/validation/ tests/unit/validation/ meson.build \
        tests/unit/testSuites.inc tests/unit/meson.build
git commit -m "feat(validation): primitive validator checks (range, enum, ordered-pair, ...)"
```

---

## Task 7: ADR 0007, ADR 0008, and docs for error + validation

**Files:**
- Create: `docs/adr/0007-config-validate-finalize-pipeline.md`
- Create: `docs/adr/0008-error-collection-over-fatalerr.md`
- Create: `docs/error-handling.md`, `docs/validation.md`
- Modify: `docs/index.md`, `docs/adr/index.md`

- [ ] **Step 1: ADR 0007 — four-stage pipeline**

Create `docs/adr/0007-config-validate-finalize-pipeline.md`:

```markdown
---
title: "ADR 0007 — Config parse-validate-finalize-adapter pipeline"
date: 2026-04-24
status: accepted
---

# ADR 0007: Four-stage config pipeline

## Context

The legacy `readswap.f90` (5161 lines) mixes four concerns: file I/O,
type deserialization, input validation, and population of model state.
That coupling makes any single concern hard to test in isolation and
blocks the planned Python bindings path, which needs to invoke
validation without going through the file reader.

## Decision

Split the input pipeline into four distinct stages:

1. **Parse** — `load_swap_config(path, config, errors)`. Byte-level
   TOML decoding; thin lossless conversions (dates, units) at the
   boundary. No semantic checks.
2. **Validate** — `config%validate(errors)`. Per-section plus
   aggregate cross-section invariants. Read-only on the config.
3. **Finalize** — `config%finalize(errors)`. Derived values, array
   expansion, canonical-form normalization. Mutates the config.
   Idempotent. Runs only when validate produced no fatal errors.
4. **Adapter** — `config_to_state(config, state)`. Temporary in
   Phase 4a; retires when state types merge with config types in a
   later phase.

A parallel fifth stage — **emit** (`write_swap_config(config, path,
errors)`) — is implemented as a separate module for round-trip
tests and future Python-driven persistence.

## Consequences

Positive:

- Each stage is unit-testable in isolation.
- Python bindings can drive validate + finalize + adapter without
  ever touching the parser.
- Round-trip tests (load → emit → load) surface a class of reader
  bugs that field-by-field diff does not.
- Errors accumulate through the whole pipeline in one collection;
  one abort checkpoint after finalize.

Negative:

- Four types and four stages is more infrastructure than legacy
  `readswap.f90`. The cost is paid once; testability and
  reusability are permanent.

## Revisit trigger

When the state layout folds into config types (Phase 4 or later),
the adapter stage disappears. The other three stages remain.
```

- [ ] **Step 2: ADR 0008 — error collection**

Create `docs/adr/0008-error-collection-over-fatalerr.md`:

```markdown
---
title: "ADR 0008 — error_collection_t over fatalerr"
date: 2026-04-24
status: accepted
---

# ADR 0008: Error collection over fatalerr

## Context

The legacy codebase uses `fatalerr` (from ttutil) for error handling:
275 call sites, each halts the process on the first failure. This
gives users one error at a time even when a config file has many
problems. Validation of a full SWAP input set should report every
issue in one pass so the user can fix them all at once.

## Decision

Introduce `error_collection_t` in `src/error/error.f90`. Every
fallible procedure in new infrastructure takes `errors` as
`intent(inout)` and appends to the collection. No direct aborts
inside stages; one explicit `abort_if_fatal` checkpoint after the
pipeline completes. Every append auto-logs via `swap_log%log_error`.

Legacy `fatalerr` sites are **not** migrated in Phase 4a. They
migrate per-module during Phase 4+ as each module is touched,
matching the existing incremental-cleanup rule.

## Consequences

Positive:

- Users see every config problem in one run.
- The collection is a value type, passable via `iso_c_binding` when
  Python bindings land.
- Errors and logs stay in sync because append auto-logs.

Negative:

- Calling convention is now mandatory: every fallible procedure
  gains an `errors` argument. Legacy code without it stays on
  `fatalerr` until migration.

## Revisit trigger

When Phase 4+ migrates the last `fatalerr` call site.
```

- [ ] **Step 3: Error-handling usage guide**

Create `docs/error-handling.md`:

```markdown
---
title: Error handling
author: SWAP modernization team
---

# Error handling

## Overview

`error_mod` (at `src/error/error.f90`) is the single error-reporting
facility for new SWAP infrastructure. Legacy modules still use
`fatalerr` (ttutil); they migrate per-module as Phase 4+ touches them.

## Contract

Every fallible procedure takes an `errors` argument as `intent(inout)`:

```fortran
subroutine read_something(...., errors)
   type(error_collection_t), intent(inout) :: errors
   ...
end subroutine
```

On failure, append an error and return. Do not halt:

```fortran
if (bad_condition) then
   call errors%append(ERR_PARSE_TYPE_MISMATCH,  &
                      "expected integer",        &
                      "drainage.basic.general.dramet")
   return
end if
```

`append` auto-logs via `swap_log%log_error`. Do not log additionally.

## Error codes

| Prefix | Domain | Typical fatal? |
|---|---|---|
| `ERR_IO_*` | Filesystem / file access | yes |
| `ERR_PARSE_*` | TOML decoding / field helpers | yes |
| `ERR_VALIDATION_*` | Per-field / cross-field / cross-section checks | yes |
| `ERR_FINALIZE_*` | Derived-value computation | often no (warnings) |
| `ERR_ADAPTER_*` | Config-to-state mismatch | yes |

`is_fatal` defaults to `.true.`. Pass `is_fatal=.false.` for
non-blocking warnings.

## Abort checkpoint

One abort point after the full pipeline (usually: after `finalize`):

```fortran
call load_swap_config(path, config, errors)
call config%validate(errors)
call config%finalize(errors)
call errors%abort_if_fatal()   ! <-- only place that may halt
```

`abort_if_fatal` writes `summary()` to stderr and calls `error stop`
if any appended error has `is_fatal=.true.`. It returns normally
otherwise.

## Testing

`error_collection_t` is a plain value type — construct one in a
test, call procedures that append to it, assert on `%has_errors`,
`%has_fatals`, `%count`, and `%items(i)%code`. See
`tests/unit/error/test_error.pf` for examples.
```

- [ ] **Step 4: Validation usage guide**

Create `docs/validation.md`:

```markdown
---
title: Validation
author: SWAP modernization team
---

# Validation

## Two layers

**Primitive checks** live in `validation_mod`
(`src/validation/validation.f90`): range, enum, ordered-pair,
not-empty, non-negative. Each is stateless and appends to a
supplied `error_collection_t`.

**Section and aggregate validators** live on the config types
themselves as type-bound procedures:

```fortran
call config%drain%validate(errors)          ! per section
call config%validate(errors)                ! aggregates all sections
```

`validate` is read-only on the config. `finalize` (a separate
type-bound procedure on the same types) is the only thing that
mutates the config post-parse.

## Primitive reference

| Primitive | Checks |
|---|---|
| `check_int_range(value, low, high, context, errors)` | `low <= value <= high` |
| `check_real_range(value, low, high, context, errors)` | same, real64 |
| `check_int_enum(value, allowed, context, errors)` | `value` in `allowed` |
| `check_not_empty(value, context, errors)` | `len_trim(value) > 0` |
| `check_nonnegative_real(value, context, errors)` | `value >= 0` |
| `check_ordered_pair(low, high, low_name, high_name, context, errors)` | `low <= high` |

Every primitive takes `context` (e.g. `"drainage.basic.general.dramet"`)
and appends the context into the resulting `error_t%context`, so error
messages point at the offending field path.

## Writing a section validator

```fortran
subroutine drainage_config_validate(self, errors)
   class(drainage_config_t),  intent(in)    :: self
   type(error_collection_t),  intent(inout) :: errors

   call check_int_enum(self%dramet, [1, 2, 3], "drainage.dramet", errors)

   if (self%dramet == 2 .and. self%swdivd == 0) then
      call errors%append(ERR_VALIDATION_CROSS_FIELD, &
                         "swdivd must be 1 when dramet=2", &
                         "drainage", is_fatal=.true.)
   end if
end subroutine
```

## Writing the aggregate validator

`swap_config_t%validate` first delegates to every section validator,
then runs cross-section invariants:

```fortran
subroutine swap_config_validate(self, errors)
   class(swap_config_t),      intent(in)    :: self
   type(error_collection_t),  intent(inout) :: errors

   call self%general%validate(errors)
   call self%simulation%validate(errors)
   call self%meteo%validate(errors)
   call self%drain%validate(errors)
   call self%soil%validate(errors)
   call self%crop%validate(errors)

   ! Cross-section rules go here.
end subroutine
```

## Testing

One pFUnit test per documented rule. Each test builds a minimal
config with the rule's inputs, calls `validate`, asserts the
expected error code and context. See `tests/unit/config/` for
examples.
```

- [ ] **Step 5: Link and verify**

Add bullets to `docs/index.md` for `error-handling.md`, `validation.md`. Add ADR 0007 and 0008 entries to `docs/adr/index.md`.

```
pixi run -e docs docs-build
```

Expected: exit 0.

- [ ] **Step 6: Commit**

```
git add docs/adr/0007-config-validate-finalize-pipeline.md \
        docs/adr/0008-error-collection-over-fatalerr.md \
        docs/error-handling.md docs/validation.md \
        docs/index.md docs/adr/index.md
git commit -m "docs(phase-4a): ADR 0007+0008, error-handling + validation guides"
```

---

## Task 8: `general_config_t`

**Files:**
- Create: `src/config/general_config.f90`
- Create: `tests/unit/config/test_general_config.pf`
- Modify: `meson.build`, `tests/unit/testSuites.inc`, `tests/unit/meson.build`

- [ ] **Step 1: Write `general_config_mod`**

Create `src/config/general_config.f90`:

```fortran
!> [general] section config: project identification, paths, screen/error switches.
module general_config_mod
   use error_mod, only: error_collection_t, ERR_VALIDATION_OUT_OF_RANGE
   use validation_mod, only: check_int_range, check_int_enum, check_not_empty
   implicit none
   private

   public :: general_config_t

   type :: general_config_t
      character(len=:), allocatable :: project
      character(len=:), allocatable :: pathwork
      character(len=:), allocatable :: pathatm
      character(len=:), allocatable :: pathcrop
      character(len=:), allocatable :: pathdrain
      integer :: swscre  = 0   !! 0=no display, 1=wb, 2=daynum
      integer :: swerror = 0   !! 0=no, 1=yes
   contains
      procedure :: validate => general_config_validate
      procedure :: finalize => general_config_finalize
   end type general_config_t

contains

   subroutine general_config_validate(self, errors)
      class(general_config_t),  intent(in)    :: self
      type(error_collection_t), intent(inout) :: errors

      ! Required strings
      if (allocated(self%project)) then
         call check_not_empty(self%project, "general.project", errors)
      else
         call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
                            "project is required", "general.project")
      end if
      if (allocated(self%pathwork)) then
         call check_not_empty(self%pathwork, "general.pathwork", errors)
      end if

      ! Enum switches
      call check_int_enum(self%swscre,  [0, 1, 2], "general.swscre",  errors)
      call check_int_enum(self%swerror, [0, 1],    "general.swerror", errors)
   end subroutine general_config_validate

   subroutine general_config_finalize(self, errors)
      class(general_config_t),  intent(inout) :: self
      type(error_collection_t), intent(inout) :: errors
      ! No derivations needed for [general].
   end subroutine general_config_finalize

end module general_config_mod
```

- [ ] **Step 2: Write the test suite**

Create `tests/unit/config/test_general_config.pf`:

```fortran
@test
subroutine test_general_minimal_valid_passes()
   use funit
   use error_mod
   use general_config_mod
   type(general_config_t)    :: g
   type(error_collection_t)  :: errors
   g%project  = "demo"
   g%pathwork = "./"
   call g%validate(errors)
   @assertFalse(errors%has_errors())
end subroutine

@test
subroutine test_general_missing_project_fails()
   use funit
   use error_mod
   use general_config_mod
   type(general_config_t)    :: g
   type(error_collection_t)  :: errors
   call g%validate(errors)
   @assertTrue(errors%has_errors())
end subroutine

@test
subroutine test_general_swscre_out_of_range_fails()
   use funit
   use error_mod
   use general_config_mod
   type(general_config_t)    :: g
   type(error_collection_t)  :: errors
   g%project = "x"
   g%swscre  = 5
   call g%validate(errors)
   @assertTrue(errors%has_errors())
   @assertEqual(ERR_VALIDATION_ENUM, errors%items(1)%code)
end subroutine

@test
subroutine test_general_finalize_idempotent()
   use funit
   use error_mod
   use general_config_mod
   type(general_config_t)    :: g
   type(error_collection_t)  :: errors
   g%project = "x"
   call g%finalize(errors)
   call g%finalize(errors)
   @assertFalse(errors%has_errors())
end subroutine
```

- [ ] **Step 3: Register and wire**

Add `'src/config/general_config.f90'` to `meson.build` `sources` list (after `src/validation/validation.f90`).

Append `ADD_TEST_SUITE(test_general_config_suite)` to `tests/unit/testSuites.inc`.

Extend `pf_files` with `'config/test_general_config.pf'`.

Extend `pfunit_extra_sources` with `'../../src/config/general_config.f90'`.

- [ ] **Step 4: Run**

```
pixi run -e test test-pfunit
```

Expected: `OK (76 tests)` (72 prior + 4 new).

- [ ] **Step 5: Commit**

```
git add src/config/general_config.f90 tests/unit/config/test_general_config.pf \
        meson.build tests/unit/testSuites.inc tests/unit/meson.build
git commit -m "feat(config): general_config_t with validate + finalize"
```

---

## Task 9: `simulation_config_t`

**Files:**
- Create: `src/config/simulation_config.f90`, `tests/unit/config/test_simulation_config.pf`
- Modify: `meson.build`, `tests/unit/testSuites.inc`, `tests/unit/meson.build`

- [ ] **Step 1: Write `simulation_config_mod`**

Create `src/config/simulation_config.f90`:

```fortran
!> [simulation] section config: dates, output timing switches.
module simulation_config_mod
   use iso_fortran_env, only: real64
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   use validation_mod, only: check_int_range, check_int_enum, check_ordered_pair
   implicit none
   private

   public :: simulation_config_t

   type :: simulation_config_t
      !> Simulation start as days-since-1900 (matches legacy).
      real(real64) :: tstart = 0.0_real64
      real(real64) :: tend   = 0.0_real64
      integer :: nprintday = 1
      integer :: swmonth   = 0
      integer :: period    = 1
      integer :: swres     = 0
      integer :: swodat    = 0
      integer :: swyrvar   = 0
   contains
      procedure :: validate => simulation_config_validate
      procedure :: finalize => simulation_config_finalize
   end type simulation_config_t

contains

   subroutine simulation_config_validate(self, errors)
      class(simulation_config_t), intent(in)    :: self
      type(error_collection_t),   intent(inout) :: errors

      call check_int_range(self%nprintday, 1, 1440, "simulation.nprintday", errors)
      call check_int_enum(self%swmonth, [0, 1], "simulation.swmonth", errors)
      call check_int_range(self%period,  0, 366, "simulation.period",  errors)
      call check_int_enum(self%swres,   [0, 1], "simulation.swres",    errors)
      call check_int_enum(self%swodat,  [0, 1], "simulation.swodat",   errors)
      call check_int_enum(self%swyrvar, [0, 1], "simulation.swyrvar",  errors)
      call check_ordered_pair(self%tstart, self%tend, &
                              "tstart", "tend", "simulation", errors)
   end subroutine simulation_config_validate

   subroutine simulation_config_finalize(self, errors)
      class(simulation_config_t), intent(inout) :: self
      type(error_collection_t),   intent(inout) :: errors
   end subroutine simulation_config_finalize

end module simulation_config_mod
```

- [ ] **Step 2: Write the test suite**

Create `tests/unit/config/test_simulation_config.pf`:

```fortran
@test
subroutine test_simulation_valid_passes()
   use funit
   use error_mod
   use simulation_config_mod
   type(simulation_config_t) :: s
   type(error_collection_t)  :: errors
   s%tstart = 37257.0d0  ! 2002-01-01
   s%tend   = 38352.0d0  ! 2004-12-31
   call s%validate(errors)
   @assertFalse(errors%has_errors())
end subroutine

@test
subroutine test_simulation_tend_before_tstart_fails()
   use funit
   use error_mod
   use simulation_config_mod
   type(simulation_config_t) :: s
   type(error_collection_t)  :: errors
   s%tstart = 38352.0d0
   s%tend   = 37257.0d0
   call s%validate(errors)
   @assertTrue(errors%has_errors())
   @assertEqual(ERR_VALIDATION_CROSS_FIELD, errors%items(1)%code)
end subroutine

@test
subroutine test_simulation_nprintday_out_of_range_fails()
   use funit
   use error_mod
   use simulation_config_mod
   type(simulation_config_t) :: s
   type(error_collection_t)  :: errors
   s%tstart = 0.0d0
   s%tend   = 1.0d0
   s%nprintday = 2000
   call s%validate(errors)
   @assertTrue(errors%has_errors())
end subroutine
```

- [ ] **Step 3: Wire**

Same pattern as Task 8: meson sources, testSuites.inc, pf_files, pfunit_extra_sources.

- [ ] **Step 4: Run**

```
pixi run -e test test-pfunit
```

Expected: `OK (79 tests)` (76 + 3).

- [ ] **Step 5: Commit**

```
git add src/config/simulation_config.f90 tests/unit/config/test_simulation_config.pf \
        meson.build tests/unit/testSuites.inc tests/unit/meson.build
git commit -m "feat(config): simulation_config_t with date-ordering check"
```

---

## Task 10: `meteorology_config_t`

**Files:** same pattern as Task 8/9.

- [ ] **Step 1: Write `meteorology_config_mod`**

Create `src/config/meteorology_config.f90`:

```fortran
!> [meteorology] section config.
module meteorology_config_mod
   use iso_fortran_env, only: real64
   use error_mod, only: error_collection_t
   use validation_mod, only: check_int_enum, check_real_range, check_not_empty
   implicit none
   private

   public :: meteorology_config_t

   type :: meteorology_config_t
      character(len=:), allocatable :: metfile
      real(real64) :: lat  = 0.0_real64
      real(real64) :: alt  = 0.0_real64
      real(real64) :: altw = 2.0_real64
      integer      :: swetr    = 0
      integer      :: swdivide = 0
      integer      :: swmetdetail = 0
      integer      :: nmetdetail  = 0
      integer      :: swrain   = 0
      integer      :: swetsine = 0
      integer      :: swinter  = 0
      integer      :: swmetfilall = 0
      real(real64) :: angstroma = 0.25_real64
      real(real64) :: angstromb = 0.50_real64
   contains
      procedure :: validate => meteorology_config_validate
      procedure :: finalize => meteorology_config_finalize
   end type meteorology_config_t

contains

   subroutine meteorology_config_validate(self, errors)
      class(meteorology_config_t), intent(in)    :: self
      type(error_collection_t),    intent(inout) :: errors

      call check_real_range(self%lat, -90.0_real64, 90.0_real64, "meteorology.lat", errors)
      call check_real_range(self%alt, -500.0_real64, 9000.0_real64, "meteorology.alt", errors)
      call check_int_enum(self%swetr,       [0, 1],    "meteorology.swetr",       errors)
      call check_int_enum(self%swdivide,    [0, 1],    "meteorology.swdivide",    errors)
      call check_int_enum(self%swmetdetail, [0, 1],    "meteorology.swmetdetail", errors)
      call check_int_enum(self%swrain,      [0, 1, 2], "meteorology.swrain",      errors)
      call check_int_enum(self%swinter,     [0, 1, 2], "meteorology.swinter",     errors)
   end subroutine meteorology_config_validate

   subroutine meteorology_config_finalize(self, errors)
      class(meteorology_config_t), intent(inout) :: self
      type(error_collection_t),    intent(inout) :: errors
   end subroutine meteorology_config_finalize

end module meteorology_config_mod
```

- [ ] **Step 2: Test suite**

Create `tests/unit/config/test_meteorology_config.pf`:

```fortran
@test
subroutine test_meteo_valid_passes()
   use funit
   use error_mod
   use meteorology_config_mod
   type(meteorology_config_t) :: m
   type(error_collection_t)   :: errors
   m%lat = 52.0d0
   m%alt = 10.0d0
   call m%validate(errors)
   @assertFalse(errors%has_errors())
end subroutine

@test
subroutine test_meteo_lat_out_of_range_fails()
   use funit
   use error_mod
   use meteorology_config_mod
   type(meteorology_config_t) :: m
   type(error_collection_t)   :: errors
   m%lat = 100.0d0
   call m%validate(errors)
   @assertTrue(errors%has_errors())
end subroutine

@test
subroutine test_meteo_swrain_invalid_fails()
   use funit
   use error_mod
   use meteorology_config_mod
   type(meteorology_config_t) :: m
   type(error_collection_t)   :: errors
   m%lat = 52.0d0
   m%swrain = 9
   call m%validate(errors)
   @assertTrue(errors%has_errors())
end subroutine
```

- [ ] **Step 3: Wire**

Same pattern.

- [ ] **Step 4: Run**

```
pixi run -e test test-pfunit
```

Expected: `OK (82 tests)`.

- [ ] **Step 5: Commit**

```
git add src/config/meteorology_config.f90 tests/unit/config/test_meteorology_config.pf \
        meson.build tests/unit/testSuites.inc tests/unit/meson.build
git commit -m "feat(config): meteorology_config_t with lat/alt/switches validation"
```

---

## Task 11: `drainage_config_t`

Mirror Task 8/9/10 pattern. Fields cover `[drainage]`: `swdra` (0/1/2), `dramet` (1/2/3), `nrlevs`, `altcu`, plus basic-drainage tables (`swdtyp`, `zbotdr`, `drares`, `infres`, `L`, `gwlinf`, `rdrain`, `rinfi`, `rentry`, `rexit`, `widthr`, `taludr`, `swallo`) as allocatable arrays sized by `nrlevs`. Validate: `swdra` in {0,1,2}; `dramet` in {1,2,3}; `nrlevs` in [0, 5]; if `dramet==2` then `swdivd==1`.

- [ ] **Step 1: Write `drainage_config_mod`** — see `src/drainage/drainage_state.f90` for field inventory and use the same names.
- [ ] **Step 2: Tests** — at minimum: `swdra` in range passes, `dramet==2 + swdivd==0` fails (cross-field).
- [ ] **Step 3: Wire** (meson sources, testSuites.inc, pf_files, pfunit_extra_sources).
- [ ] **Step 4: Run** `pixi run -e test test-pfunit`. Expected: `OK (86 tests)`.
- [ ] **Step 5: Commit** `feat(config): drainage_config_t`.

(Because the drainage schema is long, consult `src/drainage/drainage_state.f90` for the fields the state already declares — copy the field names onto `drainage_config_t` and add `validate`. Use `check_int_enum` / `check_int_range` for switches; don't over-validate array contents beyond length checks until the parity test exposes a case that fails.)

---

## Task 12: `soil_config_t`

Same pattern. Inventory fields from `src/soil/soil_state.f90`. Schema covers `[soil]`: `swsophy` (0/1), `swhyst` (0/1/2), profile layers (depth-sized arrays), per-layer hydraulic parameters (VanGenuchten or tabulated), initial condition switch (`swinco` 1/2/3).

- [ ] **Steps 1–5 as above.** Tests: `swhyst` enum; `swinco` enum; required per-layer arrays have consistent lengths.
- [ ] Expected test count: `OK (89 tests)`.
- [ ] Commit: `feat(config): soil_config_t`.

---

## Task 13: `crop_config_t`

Same pattern. Schema covers `[crop]`: `swcrop` (0/1), rotation list (array of crop entries referencing `.crp` files by path), per-entry `cropstart`, `cropend`, `cropfil`, `croptype` (1/2/3 = fixed/grass/wofost).

- [ ] **Steps 1–5 as above.** Tests: `swcrop` enum; `croptype` enum; `cropstart < cropend` per entry.
- [ ] Expected test count: `OK (92 tests)`.
- [ ] Commit: `feat(config): crop_config_t`.

---

## Task 14: `swap_config_t` aggregator

**Files:**
- Create: `src/config/swap_config.f90`, `tests/unit/config/test_swap_config.pf`
- Modify: `meson.build`, `tests/unit/testSuites.inc`, `tests/unit/meson.build`, `src/config/README.md` (new, per plan)

- [ ] **Step 1: Write the aggregator**

Create `src/config/swap_config.f90`:

```fortran
!> Aggregate SWAP configuration: composes the six section types.
module swap_config_mod
   use error_mod, only: error_collection_t
   use general_config_mod,      only: general_config_t
   use simulation_config_mod,   only: simulation_config_t
   use meteorology_config_mod,  only: meteorology_config_t
   use drainage_config_mod,     only: drainage_config_t
   use soil_config_mod,         only: soil_config_t
   use crop_config_mod,         only: crop_config_t
   implicit none
   private

   public :: swap_config_t

   type :: swap_config_t
      type(general_config_t)     :: general
      type(simulation_config_t)  :: simulation
      type(meteorology_config_t) :: meteo
      type(drainage_config_t)    :: drain
      type(soil_config_t)        :: soil
      type(crop_config_t)        :: crop
   contains
      procedure :: validate => swap_config_validate
      procedure :: finalize => swap_config_finalize
   end type swap_config_t

contains

   subroutine swap_config_validate(self, errors)
      class(swap_config_t),     intent(in)    :: self
      type(error_collection_t), intent(inout) :: errors
      call self%general%validate(errors)
      call self%simulation%validate(errors)
      call self%meteo%validate(errors)
      call self%drain%validate(errors)
      call self%soil%validate(errors)
      call self%crop%validate(errors)
      ! Cross-section rules are added here as the parity test surfaces them.
   end subroutine swap_config_validate

   subroutine swap_config_finalize(self, errors)
      class(swap_config_t),     intent(inout) :: self
      type(error_collection_t), intent(inout) :: errors
      if (errors%has_fatals()) return
      call self%general%finalize(errors)
      call self%simulation%finalize(errors)
      call self%meteo%finalize(errors)
      call self%drain%finalize(errors)
      call self%soil%finalize(errors)
      call self%crop%finalize(errors)
   end subroutine swap_config_finalize

end module swap_config_mod
```

- [ ] **Step 2: Write the test**

Create `tests/unit/config/test_swap_config.pf`:

```fortran
@test
subroutine test_swap_validate_delegates_to_sections()
   use funit
   use error_mod
   use swap_config_mod
   type(swap_config_t)      :: c
   type(error_collection_t) :: errors

   ! Leave c%general%project unallocated — the general validator will
   ! complain. That proves delegation works.
   call c%validate(errors)

   @assertTrue(errors%has_errors())
end subroutine

@test
subroutine test_swap_finalize_short_circuits_on_fatal()
   use funit
   use error_mod
   use swap_config_mod
   type(swap_config_t)      :: c
   type(error_collection_t) :: errors

   ! Seed a fatal error; finalize must not iterate further.
   call errors%append(202, "precondition fail", "x", is_fatal=.true.)
   call c%finalize(errors)

   @assertEqual(1, errors%count())
end subroutine

@test
subroutine test_swap_validate_valid_config_passes()
   use funit
   use error_mod
   use swap_config_mod
   type(swap_config_t)      :: c
   type(error_collection_t) :: errors

   c%general%project    = "demo"
   c%general%pathwork   = "./"
   c%simulation%tstart  = 37257.0d0
   c%simulation%tend    = 38352.0d0
   c%meteo%lat          = 52.0d0
   c%meteo%alt          = 10.0d0

   call c%validate(errors)

   @assertFalse(errors%has_fatals())
end subroutine
```

- [ ] **Step 3: Create `src/config/README.md` and update `docs/configuration-schema.md`**

In `docs/configuration-schema.md` (existing from Phase 2), add a new top-level section "Typed config hierarchy" after the existing TOML-keys description. Keep it short (~30 lines):

```markdown
## Typed config hierarchy

The TOML key layout in this document is mirrored by a typed config
hierarchy in `src/config/`:

| TOML section | Fortran type | Module |
|---|---|---|
| `[general]` | `general_config_t` | `general_config_mod` |
| `[simulation]` | `simulation_config_t` | `simulation_config_mod` |
| `[meteorology]` | `meteorology_config_t` | `meteorology_config_mod` |
| `[drainage]` | `drainage_config_t` | `drainage_config_mod` |
| `[soil]` | `soil_config_t` | `soil_config_mod` |
| `[crop]` | `crop_config_t` | `crop_config_mod` |
| (top-level) | `swap_config_t` | `swap_config_mod` |

Every type exposes `validate(errors)` and `finalize(errors)` as
type-bound procedures. See `docs/validation.md` for the rules and
`docs/error-handling.md` for how errors propagate.
```

Then create `src/config/README.md`:

```markdown
# src/config

## Responsibility

Typed configuration hierarchy (Phase 4a) consumed by the composite
TOML reader, validator stages, and the temporary config-to-state
adapter. One file per section plus an aggregate in `swap_config.f90`.

## Public interface

- `general_config_t` — project id, paths, screen/error switches.
- `simulation_config_t` — dates, output timing.
- `meteorology_config_t` — lat/alt, ET switches, Angstrom, interception.
- `drainage_config_t` — drainage switches, basic/extended tables.
- `soil_config_t` — profile, hydraulic params, initial conditions.
- `crop_config_t` — rotation, per-crop file refs.
- `swap_config_t` — aggregates all of the above.

Each type exposes type-bound procedures `validate(errors)` and
`finalize(errors)`. The aggregate delegates to each section then
runs cross-section invariants.

## Dependencies

`error_mod`, `validation_mod`. Reads no files; produces no
mutation outside its own fields.
```

- [ ] **Step 4: Wire and run**

Add `swap_config.f90` to meson sources AFTER the six section modules. Register suite. `pixi run -e test test-pfunit` → `OK (95 tests)`.

- [ ] **Step 5: Commit**

```
git add src/config/swap_config.f90 src/config/README.md \
        docs/configuration-schema.md \
        tests/unit/config/test_swap_config.pf \
        meson.build tests/unit/testSuites.inc tests/unit/meson.build
git commit -m "feat(config): swap_config_t aggregator"
```

---

## Task 15: `toml_field_helpers_mod`

**Files:**
- Create: `src/io/toml/toml_field_helpers.f90`, `src/io/toml/README.md`, `tests/unit/io/toml/test_toml_field_helpers.pf`
- Create: `tests/unit/io/toml/fixtures/happy.toml`, `tests/unit/io/toml/fixtures/wrong_types.toml`, `tests/unit/io/toml/fixtures/missing_keys.toml`
- Modify: `meson.build`, `tests/unit/testSuites.inc`, `tests/unit/meson.build`

- [ ] **Step 1: Write `toml_field_helpers_mod`**

Create `src/io/toml/toml_field_helpers.f90` — full module with helpers. Because the module is long, create incrementally; show the overall shape:

```fortran
!> Reusable TOML field-reading primitives with error accumulation.
!!
!! Every helper wraps a tomlf get_value call, handles null-pointer
!! cases, and appends a typed error (ERR_PARSE_*) to the collection
!! on failure. Callers never check stat/associated themselves.
module toml_field_helpers_mod
   use iso_fortran_env, only: real64
   use tomlf, only: toml_table, toml_array, toml_datetime, get_value, len
   use error_mod, only: error_collection_t,            &
                        ERR_PARSE_TYPE_MISMATCH,       &
                        ERR_PARSE_MISSING_REQUIRED
   implicit none
   private

   public :: get_required_int
   public :: get_required_real
   public :: get_required_string
   public :: get_optional_int_with_default
   public :: get_optional_real_with_default
   public :: get_optional_string_with_default
   public :: get_table
   public :: get_array_of_tables
   public :: parse_date_to_days1900

contains

   subroutine get_required_int(tab, key, out, context, errors)
      type(toml_table), pointer, intent(in)    :: tab
      character(len=*),          intent(in)    :: key
      integer,                   intent(out)   :: out
      character(len=*),          intent(in)    :: context
      type(error_collection_t),  intent(inout) :: errors
      integer :: stat
      if (.not. associated(tab)) then
         call errors%append(ERR_PARSE_MISSING_REQUIRED, &
                            "parent table missing for " // key, context)
         out = 0
         return
      end if
      call get_value(tab, key, out, stat=stat)
      if (stat /= 0) then
         call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                            "expected integer at " // key, context)
         out = 0
      end if
   end subroutine get_required_int

   subroutine get_required_real(tab, key, out, context, errors)
      type(toml_table), pointer, intent(in)    :: tab
      character(len=*),          intent(in)    :: key
      real(real64),              intent(out)   :: out
      character(len=*),          intent(in)    :: context
      type(error_collection_t),  intent(inout) :: errors
      integer :: stat
      if (.not. associated(tab)) then
         call errors%append(ERR_PARSE_MISSING_REQUIRED, &
                            "parent table missing for " // key, context)
         out = 0.0_real64
         return
      end if
      call get_value(tab, key, out, stat=stat)
      if (stat /= 0) then
         call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                            "expected real at " // key, context)
         out = 0.0_real64
      end if
   end subroutine get_required_real

   subroutine get_required_string(tab, key, out, context, errors)
      type(toml_table), pointer,     intent(in)    :: tab
      character(len=*),              intent(in)    :: key
      character(len=:), allocatable, intent(out)   :: out
      character(len=*),              intent(in)    :: context
      type(error_collection_t),      intent(inout) :: errors
      integer :: stat
      character(len=:), allocatable :: tmp
      if (.not. associated(tab)) then
         call errors%append(ERR_PARSE_MISSING_REQUIRED, &
                            "parent table missing for " // key, context)
         out = ""
         return
      end if
      call get_value(tab, key, tmp, stat=stat)
      if (stat /= 0) then
         call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                            "expected string at " // key, context)
         out = ""
         return
      end if
      out = tmp
   end subroutine get_required_string

   subroutine get_optional_int_with_default(tab, key, out, default, context, errors)
      type(toml_table), pointer, intent(in)    :: tab
      character(len=*),          intent(in)    :: key
      integer,                   intent(out)   :: out
      integer,                   intent(in)    :: default
      character(len=*),          intent(in)    :: context
      type(error_collection_t),  intent(inout) :: errors
      integer :: stat
      out = default
      if (.not. associated(tab)) return
      call get_value(tab, key, out, default=default, stat=stat)
      if (stat /= 0) then
         call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                            "expected integer at " // key, context)
         out = default
      end if
   end subroutine get_optional_int_with_default

   subroutine get_optional_real_with_default(tab, key, out, default, context, errors)
      type(toml_table), pointer, intent(in)    :: tab
      character(len=*),          intent(in)    :: key
      real(real64),              intent(out)   :: out
      real(real64),              intent(in)    :: default
      character(len=*),          intent(in)    :: context
      type(error_collection_t),  intent(inout) :: errors
      integer :: stat
      out = default
      if (.not. associated(tab)) return
      call get_value(tab, key, out, default=default, stat=stat)
      if (stat /= 0) then
         call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                            "expected real at " // key, context)
         out = default
      end if
   end subroutine get_optional_real_with_default

   subroutine get_optional_string_with_default(tab, key, out, default, context, errors)
      type(toml_table), pointer,     intent(in)    :: tab
      character(len=*),              intent(in)    :: key
      character(len=:), allocatable, intent(out)   :: out
      character(len=*),              intent(in)    :: default
      character(len=*),              intent(in)    :: context
      type(error_collection_t),      intent(inout) :: errors
      integer :: stat
      character(len=:), allocatable :: tmp
      out = default
      if (.not. associated(tab)) return
      call get_value(tab, key, tmp, stat=stat)
      if (stat /= 0) then
         out = default
         return
      end if
      out = tmp
   end subroutine get_optional_string_with_default

   !> Look up a sub-table. Returns a null pointer on failure (optional section).
   subroutine get_table(tab, key, out, context, errors)
      type(toml_table), pointer, intent(in)    :: tab
      character(len=*),          intent(in)    :: key
      type(toml_table), pointer, intent(out)   :: out
      character(len=*),          intent(in)    :: context
      type(error_collection_t),  intent(inout) :: errors
      integer :: stat
      out => null()
      if (.not. associated(tab)) return
      call get_value(tab, key, out, requested=.false., stat=stat)
   end subroutine get_table

   !> Fetch an array-of-tables under key. Null on absence.
   subroutine get_array_of_tables(tab, key, out, context, errors)
      type(toml_table), pointer, intent(in)    :: tab
      character(len=*),          intent(in)    :: key
      type(toml_array), pointer, intent(out)   :: out
      character(len=*),          intent(in)    :: context
      type(error_collection_t),  intent(inout) :: errors
      integer :: stat
      out => null()
      if (.not. associated(tab)) return
      call get_value(tab, key, out, requested=.false., stat=stat)
   end subroutine get_array_of_tables

   !> Convert a TOML datetime to days-since-1900 (the legacy time axis).
   function parse_date_to_days1900(dtv) result(days)
      type(toml_datetime), intent(in) :: dtv
      real(real64) :: days
      integer :: y, m, d, jd, jd1900
      y = dtv%date%year
      m = dtv%date%month
      d = dtv%date%day
      ! Julian-day calculation; 1900-01-01 == 2415021.
      jd     = julian_day(y, m, d)
      jd1900 = 2415021
      days   = real(jd - jd1900, kind=real64)
      if (.not. dtv%time%is_null()) then
         days = days + (dtv%time%hour * 3600 + dtv%time%minute * 60 + dtv%time%second) / 86400.0_real64
      end if
   end function parse_date_to_days1900

   pure function julian_day(y, m, d) result(jd)
      integer, intent(in) :: y, m, d
      integer :: jd, a, yy, mm
      a  = (14 - m) / 12
      yy = y + 4800 - a
      mm = m + 12 * a - 3
      jd = d + (153 * mm + 2) / 5 + 365 * yy + yy / 4 - yy / 100 + yy / 400 - 32045
   end function julian_day

end module toml_field_helpers_mod
```

- [ ] **Step 2: Create fixtures**

Under `tests/unit/io/toml/fixtures/`:

`happy.toml`:
```toml
[section]
my_int    = 42
my_real   = 3.14
my_string = "hello"
my_date   = 2002-01-01
```

`wrong_types.toml`:
```toml
[section]
my_int  = "oops"
my_real = "also oops"
```

`missing_keys.toml`:
```toml
[section]
```

- [ ] **Step 3: Write tests**

Create `tests/unit/io/toml/test_toml_field_helpers.pf` with: happy-path get_required_int / real / string, missing-key error, type-mismatch error, optional-with-default falls back. ~12 tests total.

Template for one:

```fortran
@test
subroutine test_get_required_int_happy()
   use funit
   use tomlf, only: toml_table, toml_load
   use toml_field_helpers_mod
   use error_mod
   type(toml_table), allocatable :: doc
   type(toml_table), pointer     :: sec
   type(error_collection_t)      :: errors
   integer :: out

   call toml_load(doc, 'tests/unit/io/toml/fixtures/happy.toml')
   call get_value(doc, 'section', sec, requested=.false.)
   call get_required_int(sec, 'my_int', out, 'section.my_int', errors)

   @assertEqual(42, out)
   @assertFalse(errors%has_errors())
end subroutine
```

Write the full suite — 12 tests covering the matrix of helpers × happy/wrong/missing.

- [ ] **Step 4: Wire** — add source to meson, register suite, add fixtures to git.

- [ ] **Step 5: Run**

```
pixi run -e test test-pfunit
```

Expected: `OK (107 tests)` (95 + 12).

- [ ] **Step 6: Commit**

```
git add src/io/toml/toml_field_helpers.f90 src/io/toml/README.md \
        tests/unit/io/toml/test_toml_field_helpers.pf \
        tests/unit/io/toml/fixtures/ \
        meson.build tests/unit/testSuites.inc tests/unit/meson.build
git commit -m "feat(io/toml): field helpers with typed errors and null-ptr guards"
```

---

## Task 16: `read_general_toml`

**Files:**
- Create: `src/io/toml/read_general_toml.f90`, `tests/unit/io/toml/test_read_general_toml.pf`, `tests/unit/io/toml/fixtures/general_minimal.toml`
- Modify: `meson.build`, `tests/unit/testSuites.inc`, `tests/unit/meson.build`

- [ ] **Step 1: Write the section reader**

Create `src/io/toml/read_general_toml.f90`:

```fortran
!> Reader for the [general] section of a SWAP TOML.
module read_general_toml_mod
   use tomlf, only: toml_table
   use general_config_mod, only: general_config_t
   use toml_field_helpers_mod, only: get_table, &
                                     get_required_string, &
                                     get_optional_string_with_default, &
                                     get_optional_int_with_default
   use error_mod, only: error_collection_t
   implicit none
   private

   public :: read_general_toml

contains

   subroutine read_general_toml(doc, config, errors)
      type(toml_table), pointer, intent(in)    :: doc
      type(general_config_t),    intent(inout) :: config
      type(error_collection_t),  intent(inout) :: errors

      type(toml_table), pointer :: sec, paths

      call get_table(doc, 'general', sec, 'general', errors)
      if (.not. associated(sec)) return

      call get_required_string(sec, 'project', config%project, 'general.project', errors)
      call get_optional_int_with_default(sec, 'swscre',  config%swscre,  0, 'general.swscre',  errors)
      call get_optional_int_with_default(sec, 'swerror', config%swerror, 0, 'general.swerror', errors)

      call get_table(sec, 'paths', paths, 'general.paths', errors)
      if (associated(paths)) then
         call get_optional_string_with_default(paths, 'work',       config%pathwork,  './',  'general.paths.work',       errors)
         call get_optional_string_with_default(paths, 'atmosphere', config%pathatm,   './',  'general.paths.atmosphere', errors)
         call get_optional_string_with_default(paths, 'crop',       config%pathcrop,  './',  'general.paths.crop',       errors)
         call get_optional_string_with_default(paths, 'drain',      config%pathdrain, './',  'general.paths.drain',      errors)
      end if
   end subroutine read_general_toml

end module read_general_toml_mod
```

- [ ] **Step 2: Fixture**

Create `tests/unit/io/toml/fixtures/general_minimal.toml`:

```toml
[general]
project = "unit_test"
swscre  = 0
swerror = 0

[general.paths]
work       = "./"
atmosphere = "./"
crop       = "./"
drain      = "./"
```

- [ ] **Step 3: Test**

Create `tests/unit/io/toml/test_read_general_toml.pf`:

```fortran
@test
subroutine test_read_general_happy()
   use funit
   use tomlf, only: toml_table, toml_load
   use general_config_mod
   use read_general_toml_mod
   use error_mod
   type(toml_table), allocatable :: doc
   type(toml_table), pointer     :: doc_ptr
   type(general_config_t)        :: c
   type(error_collection_t)      :: errors

   call toml_load(doc, 'tests/unit/io/toml/fixtures/general_minimal.toml')
   doc_ptr => doc
   call read_general_toml(doc_ptr, c, errors)

   @assertFalse(errors%has_errors())
   @assertEqual("unit_test", c%project)
   @assertEqual(0, c%swscre)
   @assertEqual("./", c%pathwork)
end subroutine

@test
subroutine test_read_general_missing_section_is_silent()
   use funit
   use tomlf, only: toml_table, toml_load
   use general_config_mod
   use read_general_toml_mod
   use error_mod
   type(toml_table), allocatable :: doc
   type(toml_table), pointer     :: doc_ptr
   type(general_config_t)        :: c
   type(error_collection_t)      :: errors

   ! Fixture without a [general] section -> the reader returns quietly;
   ! later validate() will complain about missing required fields.
   call toml_load(doc, 'tests/unit/io/toml/fixtures/happy.toml')
   doc_ptr => doc
   call read_general_toml(doc_ptr, c, errors)

   @assertFalse(errors%has_errors())
   @assertFalse(allocated(c%project))
end subroutine
```

- [ ] **Step 4: Wire** (meson, testSuites, pf_files, pfunit_extra_sources).

- [ ] **Step 5: Run**

```
pixi run -e test test-pfunit
```

Expected: `OK (109 tests)`.

- [ ] **Step 6: Commit**

```
git add src/io/toml/read_general_toml.f90 \
        tests/unit/io/toml/test_read_general_toml.pf \
        tests/unit/io/toml/fixtures/general_minimal.toml \
        meson.build tests/unit/testSuites.inc tests/unit/meson.build
git commit -m "feat(io/toml): read_general_toml section reader"
```

---

## Task 17: `read_simulation_toml`

Same pattern as Task 16. Fixture contains `[simulation] start_date = 2002-01-01 end_date = 2004-12-31 nprintday = 1 ...`. Use `get_value` for dates via `type(toml_datetime)`, convert with `parse_date_to_days1900`.

- [ ] **Steps 1–5:** create reader, fixture, test (3 tests: happy, missing dates, end < start), wire, run, commit.
- [ ] Expected test count: `OK (112 tests)`.
- [ ] Commit: `feat(io/toml): read_simulation_toml`.

---

## Task 18: `read_meteorology_toml`

Same pattern. Fixture covers `[meteorology]` + `[meteorology.evapotranspiration]` + `[meteorology.temporal]`. Port the key paths from the old `readswaptoml.f90` (lines 101–140) as reference for TOML path → field mapping.

- [ ] **Steps 1–5.** Test count after: `OK (115 tests)`.
- [ ] Commit: `feat(io/toml): read_meteorology_toml`.

---

## Task 19: `read_drainage_toml`

Port and extend the existing `src/io/readdrainagetoml.f90` logic onto the new `drainage_config_t`. Existing helper subroutines `read_drainage_basic`, `read_drainage_extended` can inform the structure but live in the NEW file (do not modify the legacy one in this task; Task 29 deletes it).

- [ ] **Steps 1–5.** Test count after: `OK (118 tests)`.
- [ ] Commit: `feat(io/toml): read_drainage_toml`.

---

## Task 20: `read_soil_toml`

Consult `src/soil/soil_state.f90` for field inventory. Keys: `[soil]`, `[soil.profile]`, `[soil.hydraulic]`, `[soil.initial]`.

- [ ] **Steps 1–5.** Test count: `OK (121 tests)`.
- [ ] Commit: `feat(io/toml): read_soil_toml`.

---

## Task 21: `read_crop_toml`

Keys: `[crop]`, array-of-tables `[[crop.rotation]]` with `start`, `end`, `file`, `type`.

- [ ] **Steps 1–5.** Test count: `OK (124 tests)`.
- [ ] Commit: `feat(io/toml): read_crop_toml`.

---

## Task 22: `load_swap_config` entrypoint

**Files:**
- Create: `src/io/toml/load_swap_config.f90`, `tests/unit/io/toml/test_load_swap_config.pf`
- Modify: `meson.build`, `tests/unit/testSuites.inc`, `tests/unit/meson.build`

- [ ] **Step 1: Write the dispatcher**

```fortran
!> Top-level SWAP TOML loader. Dispatches section readers.
module load_swap_config_mod
   use tomlf, only: toml_table, toml_error, toml_load
   use swap_config_mod, only: swap_config_t
   use read_general_toml_mod,      only: read_general_toml
   use read_simulation_toml_mod,   only: read_simulation_toml
   use read_meteorology_toml_mod,  only: read_meteorology_toml
   use read_drainage_toml_mod,     only: read_drainage_toml
   use read_soil_toml_mod,         only: read_soil_toml
   use read_crop_toml_mod,         only: read_crop_toml
   use error_mod, only: error_collection_t, ERR_PARSE_MALFORMED_TOML, ERR_IO_READ_FAILED
   implicit none
   private

   public :: load_swap_config

contains

   subroutine load_swap_config(path, config, errors)
      character(len=*),          intent(in)    :: path
      type(swap_config_t),       intent(inout) :: config
      type(error_collection_t),  intent(inout) :: errors

      type(toml_table), allocatable :: doc
      type(toml_table), pointer     :: doc_ptr
      type(toml_error), allocatable :: terr

      call toml_load(doc, trim(path), error=terr)
      if (allocated(terr)) then
         call errors%append(ERR_PARSE_MALFORMED_TOML, trim(terr%message), trim(path))
         return
      end if

      doc_ptr => doc
      call read_general_toml    (doc_ptr, config%general,    errors)
      call read_simulation_toml (doc_ptr, config%simulation, errors)
      call read_meteorology_toml(doc_ptr, config%meteo,      errors)
      call read_drainage_toml   (doc_ptr, config%drain,      errors)
      call read_soil_toml       (doc_ptr, config%soil,       errors)
      call read_crop_toml       (doc_ptr, config%crop,       errors)
   end subroutine load_swap_config

end module load_swap_config_mod
```

- [ ] **Step 2: Tests**

Create `tests/unit/io/toml/test_load_swap_config.pf`:

```fortran
@test
subroutine test_load_nonexistent_file_errors()
   use funit
   use swap_config_mod
   use load_swap_config_mod
   use error_mod
   type(swap_config_t)      :: config
   type(error_collection_t) :: errors

   call load_swap_config('/this/path/does/not/exist.toml', config, errors)

   @assertTrue(errors%has_fatals())
end subroutine

@test
subroutine test_load_malformed_toml_errors()
   use funit
   use swap_config_mod
   use load_swap_config_mod
   use error_mod
   type(swap_config_t)      :: config
   type(error_collection_t) :: errors

   call load_swap_config('tests/unit/io/toml/fixtures/malformed.toml', config, errors)

   @assertTrue(errors%has_fatals())
   @assertEqual(ERR_PARSE_MALFORMED_TOML, errors%items(1)%code)
end subroutine
```

Fixture `tests/unit/io/toml/fixtures/malformed.toml`:
```toml
[general
project = "x"
```

- [ ] **Step 3: Wire** and run. Expected: `OK (126 tests)`.

- [ ] **Step 4: Commit**

```
git add src/io/toml/load_swap_config.f90 \
        tests/unit/io/toml/test_load_swap_config.pf \
        tests/unit/io/toml/fixtures/malformed.toml \
        meson.build tests/unit/testSuites.inc tests/unit/meson.build
git commit -m "feat(io/toml): load_swap_config entrypoint"
```

---

## Task 23: Author hupselbrook TOML + smoke test

**Files:**
- Create: `tests/swap-cases/toml/1.hupselbrook/swap.toml` (inside the swap-cases submodule)
- Create: `tests/swap-cases/toml/1.hupselbrook/README.md`
- Create: `tests/unit/io/toml/test_hupselbrook_loads.pf`
- Modify: `tests/unit/testSuites.inc`, `tests/unit/meson.build`

- [ ] **Step 1: Author the TOML**

With `tests/swap-cases/1.hupselbrook/swap_linux.swp.template` open side-by-side, translate every parameter into TOML under the appropriate section. Start with `[general]`, `[simulation]`, `[meteorology]`, `[drainage]`, `[soil]`, `[crop]`. Where the legacy has a `.dra` or `.crp` side file, either fold the contents into `swap.toml` or create `swap.dra.toml` / `*.crp.toml`; for the minimum viable first pass, fold everything into a single `swap.toml`.

This TOML is iteratively extended during Task 27 (parity test) as fields are discovered missing. Start with a clearly minimal version — enough that `load_swap_config` runs to completion with zero parse errors.

Commit from inside the submodule:

```
cd tests/swap-cases
git checkout main
mkdir -p toml/1.hupselbrook
# author swap.toml, README.md
git add toml/
git commit -m "toml(hupselbrook): initial TOML port"
cd ../..
git add tests/swap-cases
```

- [ ] **Step 2: Write the smoke test**

```fortran
@test
subroutine test_hupselbrook_toml_loads_clean()
   use funit
   use swap_config_mod
   use load_swap_config_mod
   use error_mod
   type(swap_config_t)      :: config
   type(error_collection_t) :: errors

   call load_swap_config('tests/swap-cases/toml/1.hupselbrook/swap.toml', config, errors)

   @assertFalse(errors%has_errors())
end subroutine

@test
subroutine test_hupselbrook_toml_validates_clean()
   use funit
   use swap_config_mod
   use load_swap_config_mod
   use error_mod
   type(swap_config_t)      :: config
   type(error_collection_t) :: errors

   call load_swap_config('tests/swap-cases/toml/1.hupselbrook/swap.toml', config, errors)
   call config%validate(errors)

   @assertFalse(errors%has_fatals())
end subroutine
```

- [ ] **Step 3: Wire** and run. Expected: `OK (128 tests)`.

- [ ] **Step 4: Commit**

```
git add tests/swap-cases tests/unit/io/toml/test_hupselbrook_loads.pf \
        tests/unit/testSuites.inc tests/unit/meson.build
git commit -m "test(io/toml): hupselbrook TOML authored; load+validate smoke test"
```

---

## Task 24: `write_swap_config` emit

**Files:**
- Create: `src/io/toml/write_swap_config.f90`
- Create: `tests/unit/io/toml/test_hupselbrook_roundtrip.pf`
- Modify: `meson.build`, `tests/unit/testSuites.inc`, `tests/unit/meson.build`

- [ ] **Step 1: Write the emitter**

Mirror the loader structure. Module `write_swap_config_mod` with `write_swap_config(config, path, errors)`. Use toml-f's `toml_table` + `set_value` + serialization. Split into section-emit helpers that mirror the section readers.

Because toml-f's table-building API is involved, keep the first implementation minimal: emit every public field on every section, one per `set_value` call. No array-of-tables handling in the first commit (document known limitation in the module header); add the array paths iteratively as the round-trip test exposes them.

- [ ] **Step 2: Round-trip test**

```fortran
@test
subroutine test_hupselbrook_roundtrip()
   use funit
   use swap_config_mod
   use load_swap_config_mod
   use write_swap_config_mod
   use error_mod
   type(swap_config_t)      :: c1, c2
   type(error_collection_t) :: errors
   character(len=*), parameter :: src = 'tests/swap-cases/toml/1.hupselbrook/swap.toml'
   character(len=*), parameter :: tmp = 'builddir/hupselbrook_roundtrip.toml'

   call load_swap_config(src, c1, errors)
   call write_swap_config(c1, tmp, errors)
   call load_swap_config(tmp, c2, errors)

   @assertFalse(errors%has_errors())
   @assertEqual(c1%general%project, c2%general%project)
   @assertEqual(c1%simulation%tstart, c2%simulation%tstart, 1.0d-12)
   ! ... one assertion per public field (~60 lines)
end subroutine
```

- [ ] **Step 3: Wire** and run. Expected: `OK (129 tests)`.

- [ ] **Step 4: Commit**

```
git add src/io/toml/write_swap_config.f90 \
        tests/unit/io/toml/test_hupselbrook_roundtrip.pf \
        meson.build tests/unit/testSuites.inc tests/unit/meson.build
git commit -m "feat(io/toml): write_swap_config emit + round-trip test"
```

---

## Task 25: `config_to_state` adapter

**Files:**
- Create: `src/core/config_to_state.f90`, `tests/unit/core/test_config_to_state.pf`
- Modify: `meson.build`, `tests/unit/testSuites.inc`, `tests/unit/meson.build`

- [ ] **Step 1: Write the adapter**

Create `src/core/config_to_state.f90`:

```fortran
!> TEMPORARY Phase 4a adapter: populate legacy swap_state_t fields
!! from the new swap_config_t. This module disappears when state
!! types fold into config types in Phase 4+.
module config_to_state_mod
   use swap_config_mod, only: swap_config_t
   use swap_state_mod,  only: swap_state_t
   implicit none
   private

   public :: config_to_state

contains

   subroutine config_to_state(config, state)
      type(swap_config_t), intent(in)    :: config
      type(swap_state_t),  intent(inout) :: state

      ! [general] -> state%time
      if (allocated(config%general%project))  state%time%project  = config%general%project
      if (allocated(config%general%pathwork)) state%time%pathwork = config%general%pathwork
      state%time%swscre  = config%general%swscre

      ! [simulation] -> state%time
      state%time%tstart = config%simulation%tstart
      state%time%tend   = config%simulation%tend

      ! [meteorology] -> state%atm
      state%atm%lat       = config%meteo%lat
      state%atm%alt       = config%meteo%alt
      state%atm%altw      = config%meteo%altw
      state%atm%swetr     = config%meteo%swetr
      state%atm%swdivide  = config%meteo%swdivide
      state%atm%swmetdetail = config%meteo%swmetdetail
      state%atm%angstroma = config%meteo%angstroma
      state%atm%angstromb = config%meteo%angstromb

      ! [drainage] -> state%drain
      state%drain%dramet = config%drain%dramet
      ! ... extend as drainage_config_t fields are added ...

      ! [soil] -> state%soil
      ! ... extend as soil_config_t fields are added ...

      ! [crop] -> state%crop
      ! ... extend as crop_config_t fields are added ...
   end subroutine config_to_state

end module config_to_state_mod
```

The plan deliberately leaves most of the body as `! ... extend ...` — the parity test in Task 27 drives the completion. Each time the parity test fails on a specific field, extend the adapter with one line and re-run.

- [ ] **Step 2: Unit tests**

One test per already-wired field; extended as the adapter grows. Start with four:

```fortran
@test
subroutine test_adapter_populates_project()
   use funit
   use swap_config_mod
   use swap_state_mod
   use config_to_state_mod
   type(swap_config_t) :: config
   type(swap_state_t)  :: state

   config%general%project = "hupsel"
   call swap_state_init(state, 1, 1)
   call config_to_state(config, state)

   @assertEqual("hupsel", trim(state%time%project))
end subroutine

@test
subroutine test_adapter_populates_dates()
   use funit
   use swap_config_mod
   use swap_state_mod
   use config_to_state_mod
   type(swap_config_t) :: config
   type(swap_state_t)  :: state

   config%simulation%tstart = 37257.0d0
   config%simulation%tend   = 38352.0d0
   call swap_state_init(state, 1, 1)
   call config_to_state(config, state)

   @assertEqual(37257.0d0, state%time%tstart, 1.0d-12)
   @assertEqual(38352.0d0, state%time%tend,   1.0d-12)
end subroutine

@test
subroutine test_adapter_populates_lat()
   use funit
   use swap_config_mod
   use swap_state_mod
   use config_to_state_mod
   type(swap_config_t) :: config
   type(swap_state_t)  :: state

   config%meteo%lat = 52.0d0
   call swap_state_init(state, 1, 1)
   call config_to_state(config, state)

   @assertEqual(52.0d0, state%atm%lat, 1.0d-12)
end subroutine

@test
subroutine test_adapter_populates_dramet()
   use funit
   use swap_config_mod
   use swap_state_mod
   use config_to_state_mod
   type(swap_config_t) :: config
   type(swap_state_t)  :: state

   config%drain%dramet = 2
   call swap_state_init(state, 1, 1)
   call config_to_state(config, state)

   @assertEqual(2, state%drain%dramet)
end subroutine
```

- [ ] **Step 3: Wire** and run. Expected: `OK (133 tests)`.

- [ ] **Step 4: Commit**

```
git add src/core/config_to_state.f90 tests/unit/core/test_config_to_state.pf \
        meson.build tests/unit/testSuites.inc tests/unit/meson.build
git commit -m "feat(core): config_to_state adapter (temporary, Phase 4a)"
```

---

## Task 26: Legacy reader entry shim for parity

**Files:**
- Create: `tests/unit/io/toml/legacy_readswap_shim.f90`

- [ ] **Step 1: Write a thin shim**

The legacy reader lives in `src/io/readswap.f90` and calls `ReadSwap` internally. For the parity test, we need a shim that takes a case directory and loads into a fresh `swap_state_t`. Create `tests/unit/io/toml/legacy_readswap_shim.f90`:

```fortran
!> Thin wrapper around the legacy readswap entry point, used only by
!! parity tests. Takes a case directory and populates a fresh
!! swap_state_t. This file is in tests/ so it never ships.
module legacy_readswap_shim_mod
   use swap_state_mod, only: swap_state_t, swap_state_init
   implicit none
   private
   public :: legacy_readswap_case
contains
   subroutine legacy_readswap_case(case_dir, state)
      character(len=*),   intent(in)    :: case_dir
      type(swap_state_t), intent(inout) :: state
      ! TODO: call the legacy ReadSwap routine from src/io/readswap.f90.
      ! Exact entrypoint determined when writing this shim; inspect
      ! src/io/readswap.f90 for the public subroutine name (likely
      ! `ReadSwap` or similar). Populate state fields via the same
      ! globals/sync mechanism the legacy path uses.
      call swap_state_init(state, 1, 1)
   end subroutine legacy_readswap_case
end module legacy_readswap_shim_mod
```

**Important note on legacy reader:** `src/io/readswap.f90` writes into the `variables` module (legacy global state) and relies on `swap_state_sync` to produce `swap_state_t`. The shim must call the legacy entry point AND run the sync routines so `state` ends up populated equivalently. If the legacy path has hidden file-path dependencies (hardcoded filenames, stdin prompts), either work around them in the shim or escalate: the parity test cannot proceed without a working legacy loader.

If the shim proves unworkable within a day of effort, escalate: option is to **derive the "legacy" result from a pyswap invocation of the SWAP binary** (which already uses `readswap` internally) and serialize its state to JSON, then load the JSON into `state_legacy` in the parity test. That's a bigger change but isolates the parity test from legacy-reader plumbing.

- [ ] **Step 2: Add the shim to test meson wiring**

In `tests/unit/meson.build`, add the shim to `pfunit_extra_sources`:
```meson
        'io/toml/legacy_readswap_shim.f90',
```

And ensure `src/io/readswap.f90` + its `ttutil_dep` ancestor are available to the test executable — the executable already links `ttutil_dep` via the existing block; extend `pfunit_extra_sources` with `../../src/io/readswap.f90` (and any module chain gfortran reports missing). This will likely pull in substantial legacy code; confirm the test executable still builds.

- [ ] **Step 3: Commit (shim only; parity test in Task 27)**

```
git add tests/unit/io/toml/legacy_readswap_shim.f90 tests/unit/meson.build
git commit -m "test(io/toml): legacy readswap shim for parity testing"
```

---

## Task 27: Hupselbrook parity test (the headline)

**Files:**
- Create: `tests/unit/io/toml/test_hupselbrook_parity.pf`
- Modify: `tests/unit/io/toml/test_hupselbrook_parity.pf` iteratively as fields are wired
- Modify: `src/config/*`, `src/io/toml/read_*.f90`, `src/core/config_to_state.f90` iteratively as the test exposes gaps

- [ ] **Step 1: Write the initial parity test**

```fortran
@test
subroutine test_hupselbrook_toml_matches_legacy()
   use funit
   use swap_state_mod
   use swap_config_mod
   use load_swap_config_mod
   use config_to_state_mod
   use legacy_readswap_shim_mod
   use error_mod

   type(swap_state_t)       :: state_legacy, state_toml
   type(swap_config_t)      :: config
   type(error_collection_t) :: errors

   ! Path 1: legacy
   call legacy_readswap_case('tests/swap-cases/1.hupselbrook/', state_legacy)

   ! Path 2: TOML
   call load_swap_config('tests/swap-cases/toml/1.hupselbrook/swap.toml', config, errors)
   call config%validate(errors)
   call config%finalize(errors)
   @assertFalse(errors%has_fatals())

   call swap_state_init(state_toml, 1, 1)
   call config_to_state(config, state_toml)

   ! Field-by-field diff. Start with a small set and extend as
   ! each assertion passes; every new failure drives a fix in
   ! the TOML, reader, validator, or adapter.
   @assertEqual(state_legacy%time%project, state_toml%time%project)
   @assertEqual(state_legacy%time%tstart,  state_toml%time%tstart, 1.0d-12)
   @assertEqual(state_legacy%time%tend,    state_toml%time%tend,   1.0d-12)
   @assertEqual(state_legacy%atm%lat,      state_toml%atm%lat,     1.0d-12)
   @assertEqual(state_legacy%atm%alt,      state_toml%atm%alt,     1.0d-12)
   @assertEqual(state_legacy%drain%dramet, state_toml%drain%dramet)
end subroutine
```

- [ ] **Step 2: Iterate: expand until full parity**

Run the test. For every failure:

1. Identify which side is wrong (legacy is the ground truth).
2. Decide the fix:
   - Missing TOML key → add to `swap.toml`.
   - Reader not populating config → extend the section reader.
   - Config field missing → add to `*_config_t`.
   - Validator rejects valid value → adjust validator.
   - Adapter not copying field → extend `config_to_state`.
3. Commit the fix with a focused message, e.g.:
   ```
   git commit -m "fix(config): add swinter enum to meteo validator"
   git commit -m "fix(toml): extend meteorology reader for [temporal] keys"
   git commit -m "fix(adapter): copy meteo.swrain into state.atm"
   ```
4. Re-run. Add the next assertion to the parity test. Repeat until every field legacy populates is asserted and matches.

This iteration typically takes multiple commits — 10-30 expected for hupselbrook depending on how many fields `readswap` touches. That is the expected shape; do not try to preempt the diff by writing all assertions upfront.

- [ ] **Step 3: Final commit (when parity passes)**

Once the parity test's full assertion set passes:

```
git add tests/unit/io/toml/test_hupselbrook_parity.pf
git commit -m "test(io/toml): hupselbrook full parity — TOML path matches legacy"
```

- [ ] **Step 4: Wire the suite**

Register `test_hupselbrook_parity_suite` in `testSuites.inc` and `pf_files` (if not already done).

Expected final count: roughly `OK (134+N tests)` where N depends on how many smaller fixes also added tests.

---

## Task 28: Author TOML projects for cases 2–6 + smoke tests

**Files:**
- Create: `tests/swap-cases/toml/2.grassgrowth/swap.toml` through `6.surfacewater/swap.toml` (in submodule)
- Create: `tests/unit/io/toml/test_all_cases_smoke.pf`

- [ ] **Step 1: Author each case's TOML**

For each case (2.grassgrowth, 3.macroporeflow, 4.oxygenstress, 5.salinitystress, 6.surfacewater):
- Open `tests/swap-cases/N.<name>/swap_linux.swp.template` side-by-side with the hupselbrook TOML as template.
- Copy the hupselbrook TOML, edit every value to match the case.
- Commit inside the submodule.

Expected time: one commit per case, ~30 minutes each (6 commits total across the submodule).

- [ ] **Step 2: One smoke test per case**

Create `tests/unit/io/toml/test_all_cases_smoke.pf`:

```fortran
@test
subroutine test_grassgrowth_loads_clean()
   use funit
   use swap_config_mod
   use load_swap_config_mod
   use error_mod
   type(swap_config_t)      :: c
   type(error_collection_t) :: e
   call load_swap_config('tests/swap-cases/toml/2.grassgrowth/swap.toml', c, e)
   call c%validate(e)
   call c%finalize(e)
   @assertFalse(e%has_fatals())
end subroutine

! ... repeat for cases 3, 4, 5, 6 ...
```

Five tests total.

- [ ] **Step 3: Wire, run, commit**

```
pixi run -e test test-pfunit
```

Expected: all smoke tests pass. If any fails, author/adjust the TOML until it does (do NOT extend the reader or validator for case-specific quirks; the infrastructure must handle the whole SWAP schema uniformly).

```
git add tests/swap-cases tests/unit/io/toml/test_all_cases_smoke.pf \
        tests/unit/testSuites.inc tests/unit/meson.build
git commit -m "test(io/toml): smoke tests for cases 2-6 TOML projects"
```

---

## Task 29: Retire legacy TOML readers

**Files:**
- Delete: `src/io/readswaptoml.f90`, `src/io/readdrainagetoml.f90`
- Delete: `tests/unit/io/test_readswaptoml.pf`, `tests/unit/io/test_readdrainagetoml.pf`
- Delete: `tests/unit/io/fixtures/minimal_swap.toml`, `minimal_drainage.toml`, `malformed_drainage.toml`
- Modify: `meson.build` (remove deleted sources), `tests/unit/meson.build` (remove entries), `tests/unit/testSuites.inc` (remove entries)

- [ ] **Step 1: Delete files**

```
git rm src/io/readswaptoml.f90 src/io/readdrainagetoml.f90
git rm tests/unit/io/test_readswaptoml.pf tests/unit/io/test_readdrainagetoml.pf
git rm tests/unit/io/fixtures/minimal_swap.toml \
       tests/unit/io/fixtures/minimal_drainage.toml \
       tests/unit/io/fixtures/malformed_drainage.toml
```

- [ ] **Step 2: Remove from meson**

In `meson.build` `sources` list, remove:
```
    'src/io/readdrainagetoml.f90',
    'src/io/readswaptoml.f90',
```

In `tests/unit/meson.build`:
- remove `'io/test_readswaptoml.pf'` and `'io/test_readdrainagetoml.pf'` from `pf_files`
- remove `'../../src/io/readdrainagetoml.f90'` and `'../../src/io/readswaptoml.f90'` from `pfunit_extra_sources`

In `tests/unit/testSuites.inc`, remove the two lines:
```
ADD_TEST_SUITE(test_readdrainagetoml_suite)
ADD_TEST_SUITE(test_readswaptoml_suite)
```

- [ ] **Step 3: Verify build and tests**

```
pixi run -e test check-fast
```

Expected: build succeeds, pFUnit suite passes (35-ish tests fewer than before this task), regression suite green.

- [ ] **Step 4: Commit**

```
git add meson.build tests/unit/meson.build tests/unit/testSuites.inc
git commit -m "chore(io): retire legacy readswaptoml and readdrainagetoml"
```

---

## Task 30: Phase 4a closeout

**Files:** no new files; verifies and tags.

- [ ] **Step 1: Run full regression**

```
pixi run -e test check-full
```

Expected: 6/6 regression cases green, pFUnit suite green, wall time ~10 minutes. If red: stop, fix, do not tag.

- [ ] **Step 2: Re-run coverage baseline**

```
pixi run -e coverage coverage-report
```

Update `docs/coverage-baseline.md` with the new numbers. Verify no per-domain drop vs. Phase 3 baseline (Phase 3: 51.4% project-wide). Expected: coverage up across `src/error/`, `src/validation/`, `src/config/`, `src/io/toml/`; unchanged elsewhere.

```
git add docs/coverage-baseline.md
git commit -m "docs(baselines): record phase 4a coverage numbers"
```

- [ ] **Step 3: Record check-full baseline**

Append a Phase 4a section to `tests/regression/baselines/phase-3-coverage.log` (or create `phase-4a-infrastructure.log` following the same format). Capture per-case wall times.

```
git add tests/regression/baselines/
git commit -m "docs(baselines): record phase 4a check-full output"
```

- [ ] **Step 4: Verify "no physics changed" constraint**

```
git diff rescue/phase-3-coverage..HEAD -- \
   src/atmosphere/ src/soil/ src/crop/ src/boundary/ \
   src/drainage/ src/macropore/ src/solute/ src/heat/ src/utils/
```

Expected: empty diff (or only new `README.md` lines). If any production `.f90` changed, investigate and revert unless it's a clearly infrastructural fix (unlikely — escalate if found).

- [ ] **Step 5: Fast-forward main to development**

```
git checkout main
git merge --ff-only development
git checkout development
```

- [ ] **Step 6: Tag**

```
git tag rescue/phase-4a-infrastructure
```

No push.

- [ ] **Step 7: Summary for user**

Report:
- Phase 4a complete; tag `rescue/phase-4a-infrastructure` at commit `<SHA>`.
- Final pFUnit count: `<N>` across `<M>` suites.
- Coverage numbers vs Phase 3 baseline.
- Any latent legacy bugs surfaced during parity testing.
- Deferred to pre-Phase-4 follow-on: full parity tests for cases 2-6.
- `main` = `development` locally; `origin/main` untouched.

---

## Self-Review

### Spec coverage

- Logger review + `docs/logging.md` — Task 1 ✓
- `src/error/` module + tests — Tasks 2–5 ✓
- `src/validation/` module + tests — Task 6 ✓
- ADRs 0007 and 0008 — Task 7 ✓
- `docs/error-handling.md`, `docs/validation.md` — Task 7 ✓
- Per-section `*_config_t` types (six) — Tasks 8–13 ✓
- `swap_config_t` aggregator — Task 14 ✓
- `toml_field_helpers` — Task 15 ✓
- Six thematic section readers — Tasks 16–21 ✓
- `load_swap_config` entrypoint — Task 22 ✓
- hupselbrook TOML + smoke tests — Task 23 ✓
- `write_swap_config` emit + round-trip — Task 24 ✓
- `config_to_state` adapter — Task 25 ✓
- Hupselbrook parity test (headline) — Tasks 26–27 ✓
- TOML projects for cases 2–6 + smoke tests — Task 28 ✓
- Retire legacy TOML readers — Task 29 ✓
- Coverage re-run + phase tag — Task 30 ✓
- `configuration-schema.md` update — covered under the `src/config/README.md` creation in Task 14; update `docs/configuration-schema.md` during Task 14 as well. *Plan gap noted — adding to Task 14 step 3.*

### Placeholder scan

- Task 11/12/13 condense steps with phrases like "Steps 1–5 as above" — mitigated by the first detailed example (Task 8) carrying the full template. Not ideal by the "repeat code" rule but preserves plan size.
- Task 19 says "Port and extend the existing `src/io/readdrainagetoml.f90` logic" — does not repeat the code, but points to an existing file with clear reuse instructions and concrete test-count targets. Acceptable.
- Task 25 adapter body has explicit `! ... extend ...` placeholders; these are iteratively filled in Task 27. Documented.
- Task 26 shim has a `! TODO` inside the initial shim body; resolved iteratively; escalation path documented.

### Type consistency

- `general_config_t`, `simulation_config_t`, `meteorology_config_t`, `drainage_config_t`, `soil_config_t`, `crop_config_t`, `swap_config_t` — all consistent across Tasks 8-14, 22, 25.
- `error_collection_t` methods named `append`, `has_errors`, `has_fatals`, `count`, `summary`, `abort_if_fatal`, `clear` — consistent across all tasks using them.
- `validate` and `finalize` are procedures (not different names between types) on every `*_config_t`.
- `load_swap_config(path, config, errors)` signature consistent in Tasks 22, 23, 24, 27.
- `config_to_state(config, state)` signature consistent in Tasks 25, 27.
- `write_swap_config(config, path, errors)` signature consistent in Task 24.

### Known risks during execution

- Task 26 (legacy shim) may expose legacy readswap's hidden file dependencies and require non-trivial adaptation. Plan documents the escalation option (pyswap-driven state snapshot).
- Task 27 (parity) is iterative by design; total commit count for the task depends on how many fields legacy populates. Realistic estimate: 15-40 commits.
- Task 28 (case 2-6 TOMLs) depends on hupselbrook being clean first; if hupselbrook stalls, the whole case-authoring work stalls.
- If any task builds pull in so many transitive dependencies that the test executable becomes unwieldy (like Phase 3 Task 17), apply the same pattern: pull the minimum into `pfunit_extra_sources`, add stub files in `tests/` if necessary.

Fix applied inline: added `docs/configuration-schema.md` update to Task 14 step 3. (Before committing the final task outputs.)
