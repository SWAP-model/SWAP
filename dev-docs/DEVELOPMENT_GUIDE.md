# SWAP Modernization, Parallelization, and AI Collaboration Guide

## 1. Purpose and Scope

This document defines how to modernize and extend the SWAP (Soil Water – Atmosphere – Plant) hydrological model:

- Modern Fortran style and refactoring rules.
- Multicore / GPU‑oriented design.
- Strong preference for in‑memory operation with minimal I/O (Python/other frontends).
- Constraints on modifying core physics.
- Incremental, test‑driven workflow.
- How to collaborate effectively with AI assistants (e.g. VS Code “agents”).

All contributors (human and AI) must follow these guidelines.

> **Status (2026-05).** The hardest part of this plan is **done**: the
> strangler-fig migration to `TOML → typed config → typed state` is complete,
> the bare-globals module (`variables.f90`/`arrays.fi`) and the legacy output
> C-API (`swapoutput.f90`) are deleted, and ten subsystems plus the globals
> sweeps have moved onto typed `state%X` / `config%X` records. The SAVE/COMMON
> elimination below therefore describes a pattern that is now the *norm* for
> existing code; the live frontier is keeping new code SAVE-free and reaching
> in-process multi-instance + BMI/Python execution. See
> `dev-docs/post-phase-4-modernization-summary.md` for the full arc.

> **This is the deep reference.** For the day-to-day operating contract
> (non-negotiables, exact build/test commands, current architecture, commit
> conventions, where-to-look), see **`CLAUDE.md`** in the repo root. This guide
> is the depth behind that summary: detailed Fortran style, in-memory/Python
> design, parallelization patterns, and the code-review checklist.


## 2. Project Management and Tooling

### 2.1 Environment Management (Pixi)

- The canonical development environment is defined in `pixi.toml`.
- All commands (build, test, format, benchmarks, regression) must be runnable via `pixi run …`.
- Do not introduce ad‑hoc commands in docs that bypass Pixi; always show `pixi run <task>`.

**Core Pixi tasks** (see `pixi.toml` for the full list; the build/test gate
is also summarized in `CLAUDE.md`):

```bash
pixi run build-linux          # meson compile (auto-configures builddir)
pixi run -e test check-fast   # pFUnit + 4 regression cases
pixi run -e test check-full   # pFUnit + all regression cases
pixi run -e test test-pfunit  # unit suite only
pixi run lint                 # fprettify (lint-check for a dry-run diff)
pixi run clean                # rm -rf builddir
```

If you add new tooling (formatters, linters, test runners), expose them as new Pixi tasks and briefly document them.

### 2.2 Build System (Meson)

- Meson is the single source of truth for building SWAP.
- Do not add raw compiler command lines to documentation or comments; always show Meson + Pixi usage.

**Typical workflow:**

```bash
pixi run build-linux
pixi run -e test check-fast
```

**Meson conventions:**

- Root project file: `meson.build`.
- One `meson.build` per source directory under `src/`.
- Group sources by domain: `src/soil`, `src/crop`, `src/hydraulics`, `src/io`, etc.
- Use Meson options and features for:
  - Compile flags (OpenMP, GPU offload).
  - Optional components (Python bindings, tests).

---

## 3. Modern Fortran Principles

### 3.1 General Coding Rules

- Use free‑form Fortran (`.f90`, `.F90`).
- `implicit none` in every module and procedure.
- Prefer `use iso_fortran_env, only: real64, int32` for kinds.
- Prefer descriptive names except for simple indices (`i`, `j`, `k`).
- One logical module per file; filename should match module name where practical.

**Example:**

```fortran
use iso_fortran_env, only: real64, int32

real(real64)  :: rainfall_mm_day    ! [mm/day]
real(real64)  :: soil_moisture_vv   ! [cm^3/cm^3]
integer(int32):: n_layers
```

### 3.2 Eliminating SAVE and Hidden State

Legacy SWAP code used `SAVE`, `DATA`, and `COMMON` for persistent state. These
patterns are problematic for testability, parallel execution (multicore, GPU),
and interfacing from Python with independent runs and parameter changes.

This was the central work of the modernization and is **largely complete** —
the bare-globals module is gone and persistent state lives on typed
`*_state_t` records (see `src/state/`). The rules below are now the standing
norm: keep new code free of hidden state, and when you touch a surviving
legacy pocket, migrate it the same way.

**Principles:**

- Replace `SAVE` + `DATA` with explicit state passing via derived types
  (one `*_state_t` per subsystem, aggregated under `swap_state_t`).
- Use module‑level configuration only when a value is truly global and mostly
  read‑only — and prefer a field on `swap_config_t` over a module global.
- Keep computational kernels as stateless as possible (good candidates for
  `pure`/`elemental`).

**Bad example (legacy style):**

```fortran
subroutine update_balance(rain, balance)
  real, intent(in)  :: rain
  real, intent(out) :: balance
  real :: cum_rain
  save cum_rain
  data cum_rain /0.0/

  cum_rain = cum_rain + rain
  balance  = cum_rain
end subroutine
```

**Good example (explicit state, in‑memory):**

```fortran
module water_balance_mod
  use iso_fortran_env, only: real64
  implicit none
  private

  type, public :: water_balance_state_t
    real(real64) :: cumulative_rain = 0.0_real64
  end type water_balance_state_t

  public :: update_balance

contains

  subroutine update_balance(state, rain, balance)
    type(water_balance_state_t), intent(inout) :: state
    real(real64),                intent(in)    :: rain    ! [mm/day]
    real(real64),                intent(out)   :: balance ! [mm]
    state%cumulative_rain = state%cumulative_rain + rain
    balance               = state%cumulative_rain
  end subroutine update_balance

end module water_balance_mod
```

### 3.3 State Sharing and In‑Memory Operation

The Fortran core must support being called many times from Python (or other frontends) with:

- Shared, immutable inputs (e.g. precipitation time series, soil properties).
- Small changes in a subset of parameters between runs (e.g. 1 of 100 parameters).
- Minimal disk I/O.

**Design goals:**

- The computational core of SWAP has to remain the same for backwards compatibility with legacy way of working with the code
- The ambition is to have the core of SWAP modernizedin a way that it can be called in the "legac way", from the executable, but also through BMI from programming languages like Python (pyswap package).
- For that we should have the possibility to export the core computing library separately (perhaps with a sort of an interface reading and setting the variables from memory instead of files like ttutil).
- All "run configuration" (precipitation, soil, crop parameters, control flags) is represented as in‑memory derived types.
- The Fortran API accepts these types as arguments; no reading of configuration files inside core routines.
- Simulation outputs are returned via arguments (arrays / derived types) or via pre‑allocated buffers.

**Example state definitions:**

```fortran
module swap_state_mod
  use iso_fortran_env, only: real64, int32
  implicit none
  private

  type, public :: meteo_series_t
    integer(int32)          :: n_steps
    real(real64), allocatable :: precip_mm_day(:)
    real(real64), allocatable :: t_air_c(:)
    ! add other forcings as needed
  end type meteo_series_t

  type, public :: soil_profile_t
    integer(int32)                 :: n_layers
    real(real64), allocatable :: thickness_cm(:)
    real(real64), allocatable :: theta_init(:)
    ! etc.
  end type soil_profile_t

  type, public :: parameter_set_t
    ! This is where the “100 variables” live,
    ! with some being changed between runs
    real(real64) :: param1
    real(real64) :: param2
    ! ...
  end type parameter_set_t

end module swap_state_mod
```

**High‑level simulation interface:**

```fortran
module swap_driver_mod
  use iso_fortran_env, only: real64, int32
  use swap_state_mod
  implicit none
  private
  public :: run_swap

contains

  subroutine run_swap(meteo, soil, params, state_out)
    type(meteo_series_t), intent(in)  :: meteo
    type(soil_profile_t), intent(in)  :: soil
    type(parameter_set_t),intent(in)  :: params
    ! state_out could be a time series of key outputs
    real(real64),         intent(inout) :: state_out(:,:) 
    ! Implementation: call existing core computational routines.
    ! Minimize I/O here; all I/O is orchestrated at Python or outer level.
  end subroutine run_swap

end module swap_driver_mod
```

In this design, Python code:

1. Reads configuration files once.
2. Constructs meteo, soil, params objects (via FFI).
3. Calls `run_swap` many times, only modifying `params` for each run.
4. Collects results from `state_out` without disk I/O.

### 3.4 I/O Minimization Strategy

**Rules:**

- All core compute routines must be free of file I/O.
- I/O is allowed only in:
  - Dedicated `*_io_mod` modules.
  - High‑level "driver" layers, not in numerically hot kernels.
- Treat file I/O as a front‑end concern (e.g. Python or a thin Fortran wrapper).

**Refactoring pattern:**

1. Identify routines that both compute and read/write files.
2. Split into:
   - Pure computation routine: operates only on arguments.
   - Thin I/O routine: reads/writes files and calls computation.

**Example:**

```fortran
! Old (to refactor):
subroutine compute_evaporation_with_io(input_file, output_file)
  ! reads meteorological data, computes, writes results
end subroutine

! New:
subroutine compute_evaporation(meteo, params, evap)
  ! pure-ish kernel, no file I/O
end subroutine

subroutine evaporation_from_file(input_file, output_file)
  ! thin I/O wrapper
  ! reads into meteo, params; calls compute_evaporation; writes evap
end subroutine
```

For Python integration, the recommended pattern is to use `compute_evaporation` only and keep file handling outside Fortran.

### 3.5 Protecting Core Physics

Many SWAP routines encode complex physical equations that we do not want to change unless absolutely necessary.

**Rules:**

- Consider numerically core routines as frozen physics kernels.
- Allowed changes inside physics kernels:
  - Replace `SAVE`, `COMMON`, implicit state with explicit arguments/state types, if this does not alter equations or ordering.
  - Improve interfaces (kinds, intents, purity) without changing formulas.
  - Minor reordering only if it is guaranteed not to change results within numerical precision.
- Changes that should be avoided or heavily scrutinized:
  - Modifying equations, coefficients, or algorithmic logic.
  - Changing variable semantics or units.
  - Optimizations that change floating‑point operation ordering in sensitive parts without tests.

Where possible, wrap legacy kernels instead of rewriting them:

```fortran
! legacy_physics_mod.f90 (minimally touched)
module legacy_physics_mod
  implicit none
contains
  subroutine core_soil_flux( ... )
    ! original code, lightly cleaned (implicit none, kinds, etc.)
  end subroutine
end module legacy_physics_mod

! modern wrapper:
module soil_flux_wrapper_mod
  use legacy_physics_mod
  implicit none
contains
  subroutine soil_flux_api(state_in, params, flux_out)
    ! Translate from modern derived types to arguments expected by core_soil_flux
    ! Call core_soil_flux
    ! Translate outputs back
  end subroutine
end module soil_flux_wrapper_mod
```

---

## 4. Multicore and GPU‑Oriented Design

### 4.1 Multicore (CPU) Principles

- Target data‑parallel loops: across layers, grid cells, parameters, or ensembles.
- No shared mutable globals inside parallel loops.
- Use OpenMP or other parallel frameworks in outer loops; keep inner kernels simple and stateless.

**Example OpenMP pattern:**

```fortran
subroutine run_ensemble(n_runs, meteo, soil, params_all, outputs)
  use iso_fortran_env, only: real64, int32
  use swap_state_mod
  use swap_driver_mod
  implicit none
  integer(int32),       intent(in)  :: n_runs
  type(meteo_series_t), intent(in)  :: meteo
  type(soil_profile_t), intent(in)  :: soil
  type(parameter_set_t),intent(in)  :: params_all(n_runs)
  real(real64),         intent(out) :: outputs(:, :, :)  ! e.g. [n_runs, n_steps, n_vars]
  integer(int32)                    :: r

!$omp parallel do default(none) private(r) shared(n_runs, meteo, soil, params_all, outputs)
  do r = 1, n_runs
    call run_swap(meteo, soil, params_all(r), outputs(r, :, :))
  end do
!$omp end parallel do

end subroutine run_ensemble
```

This design supports:

- Running many parameter sets in parallel.
- Minimal I/O (none inside `run_swap`).
- Clean mapping from Python (e.g. pass `params_all` as a 1D array of structures).

### 4.2 GPU‑Friendly Patterns

While the exact GPU backend may vary, the following patterns make GPU offload easier:

- Keep numerical kernels as pure subroutines/functions with array arguments.
- Avoid allocating/deallocating inside kernels; operate on pre‑allocated arrays.
- Prefer structure‑of‑arrays for large grids if GPU memory access patterns benefit.

**Example kernel:**

```fortran
pure subroutine compute_flux(n, theta, ksat, flux)
  use iso_fortran_env, only: real64, int32
  integer(int32), intent(in)  :: n
  real(real64),   intent(in)  :: theta(n)
  real(real64),   intent(in)  :: ksat(n)
  real(real64),   intent(out) :: flux(n)
  integer(int32)              :: i

  do i = 1, n
    flux(i) = ksat(i) * theta(i)  ! simplified
  end do
end subroutine compute_flux
```

Such kernels can later be wrapped for GPU with directives or offload mechanisms, without changing their interface.

---

## 5. Incremental, Test‑Driven Workflow

Modernization must proceed in small, verifiable steps.

### 5.1 Rules for Changes

For each change (human or AI‑assisted):

1. Keep the scope small (e.g. one subroutine, one module, one I/O pathway).
2. Ensure there is at least one unit or regression test that exercises the changed behavior.
3. Ensure that we can compare outputs against baseline (for physics‑sensitive routines).
4. Avoid mixing refactoring with new features in the same change set.

### 5.2 Unit Tests Around Changes

Whenever you:

- Remove a `SAVE`.
- Split I/O from computation.
- Introduce a new derived type or state API.
- Add parallelization (OpenMP, etc.).
- Adjust the calling interface for Python.

…you must add or adapt tests that verify:

- The numerics behave as before (within tolerances).
- The new API behaves as expected (e.g. multiple runs with different `parameter_set_t` produce distinct outputs).
- Thread safety where relevant (e.g. run the same code in parallel on multiple copies of state).

**Example testing pattern (pseudocode):**

```fortran
module test_run_swap
  use iso_fortran_env, only: real64
  use swap_state_mod
  use swap_driver_mod
  implicit none
contains

  subroutine test_run_swap_different_params()
    type(meteo_series_t)  :: meteo
    type(soil_profile_t)  :: soil
    type(parameter_set_t) :: p1, p2
    real(real64), allocatable :: out1(:,:), out2(:,:)

    ! Arrange: initialize meteo, soil, p1, p2 (p2 slightly different)
    ! Allocate out1, out2

    ! Act
    call run_swap(meteo, soil, p1, out1)
    call run_swap(meteo, soil, p2, out2)

    ! Assert: out1 and out2 differ in expected metrics, but both are numerically stable
  end subroutine

end module test_run_swap
```

---

## 6. Working with Agentic AI in VS Code

### 6.1 Role of AI

AI assistants act as pair programmers and refactoring helpers, not authoritative sources. You are responsible for:

- Validating all AI‑generated code.
- Ensuring compliance with this guide.
- Writing tests for every substantive change.

### 6.2 Prompting AI Effectively

When asking an AI to refactor or write code:

1. **Provide context:**
   - The SWAP component (e.g. "soil hydraulic conductivity routine, no physics changes allowed").
   - The objective (e.g. "remove SAVE and separate I/O from computation").

2. **Provide constraints:**
   - "No modification of core equations."
   - "Use explicit state passing and derived types."
   - "No file I/O inside this routine."
   - "Thread‑safe and suitable for OpenMP across parameter sets."

3. **Provide tests or expected behavior:**
   - Baseline numerical outputs or acceptance criteria.
   - A short description of what must remain identical.

**Example prompt:**

> Refactor this SWAP subroutine to remove SAVE and make it safe for OpenMP across different parameter sets.
> Constraints:
> – Do not change any physical formulas or equations.
> – No file I/O inside the routine; all inputs and outputs must be arguments.
> – Use real(real64) and implicit none.
> Also propose a small unit test that checks that two runs with different parameter sets produce different outputs while staying numerically stable.

### 6.3 Reviewing AI Output

For any AI‑generated proposal:

1. **Check for:**
   - `implicit none`.
   - Correct intent attributes.
   - Absence of `SAVE`, `COMMON`, implicit external procedures.
   - No file I/O in computational kernels.
   - No changes to equations unless explicitly requested.

2. **Check that:**
   - Interfaces are usable from Python (arguments, no hidden state).
   - The changes are small and testable.

3. Add or update tests before accepting the change.

---

## 7. Code Review Checklist

Use this checklist for all MRs/PRs and major commits.

### 7.1 Correctness & Style

- [ ] `implicit none` in all updated modules and procedures.
- [ ] No `SAVE`, `COMMON`, or hidden global mutable state introduced or left in modernized code.
- [ ] No `GOTO` for normal control flow.
- [ ] `iso_fortran_env` kinds used instead of `REAL*8`/`INTEGER*4`.
- [ ] Clear, descriptive names and units documented in comments.

### 7.2 Physics Integrity

- [ ] Core equations are unmodified unless explicitly agreed and documented.
- [ ] Any change in physics is accompanied by justification and tests.
- [ ] Wrappers are preferred over rewriting legacy physics routines.

### 7.3 Parallel & In‑Memory Design

- [ ] Numerical kernels have no file I/O.
- [ ] State is passed via arguments / derived types, suitable for Python FFI.
- [ ] Parallel loops do not rely on shared mutable global state.
- [ ] No allocations inside tight parallel loops.

### 7.4 Testing & Incremental Workflow

- [ ] Each change is small and logically cohesive.
- [ ] New or changed behavior is covered by unit tests and/or regression tests.
- [ ] Baseline outputs (where available) have been compared and are within acceptable tolerance.
- [ ] `pixi run -e test check-fast` is green (pFUnit + regression); a new
      `.pf` suite is registered in `tests/unit/testSuites.inc`.

### 7.5 AI‑Related Checks

- [ ] AI‑generated code was fully reviewed and adapted as needed.
- [ ] AI suggestions did not violate in‑memory, low‑I/O, or no‑SAVE rules.
- [ ] Any new patterns introduced by AI are consistent with this guide.

---

## 8. Document Maintenance

- This guide lives at `dev-docs/DEVELOPMENT_GUIDE.md`; the operating-contract
  summary lives at `CLAUDE.md` in the repo root. Keep the two consistent — if
  a non-negotiable or command changes, update `CLAUDE.md`; if a style/design
  rule changes, update this guide.
- Update when:
  - Python/BMI interface requirements change.
  - GPU or parallelization strategy evolves.
  - New patterns for low‑I/O, in‑memory operation are adopted.
- Discuss and review changes to this guide like any other code change.
