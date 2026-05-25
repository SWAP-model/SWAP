---
title: "Driver Modernization & BMI Stub — Design Spec"
date: 2026-05-12
status: draft
arc: SS-DRV (driver modernization)
phase: 1 of 3 (Phase 2 = Python multiprocessing ensemble; Phase 3 = swap_ensemble_mod + OpenMP)
---

# Driver Modernization & BMI Stub

## Context

`src/core/swap.f90` is the runtime driver: a standalone subroutine that runs the full SWAP lifecycle in three phases (`iTask=1` init, `iTask=2` time loop, `iTask=3` close). It predates the project's move toward modular Fortran and carries several anti-patterns that block both clean unit testing and external coupling:

- **Not a module**: `subroutine swap` is a free subroutine. Callers must declare an explicit `interface` block (`swap_main.f90:23–30`) to call it. No module means no `use` statement, no `.mod` file, and no clean way to expose it via `iso_c_binding`.
- **Magic integer dispatch**: `iTask=1/2/3` for lifecycle phases, `iCaller=0/non-zero` for "is this a DLL call?". `handle_exchange` uses task codes `11/21/22/23/29/31`. Untestable, unreadable.
- **Bespoke DLL exchange**: `swap_exchange` module (lines 1–47) defines `swap_input`/`swap_output` derived types for external calling, with fixed-size arrays `dimension(500)`. `handle_exchange` (lines 557–711) populates them from module globals. The whole pattern is a hand-rolled precursor to BMI.
- **Time loop is buried inside**: the outer `do while (.not.flrunend)` loop sits inside `iTask=2`. A caller cannot drive one timestep at a time — required for BMI `update()`, coupling with imod_coupler, and tight integration with pyswap.
- **Compiler flag `-std=legacy`**: project-wide, permitting implicit typing, `!dec$` directives (line 63 of `swap.f90` carries one), and other patterns that wouldn't compile under modern standards.

This arc is the foundation of a three-phase modernization plan:

1. **Phase 1 (this arc)**: Module refactor + minimal BMI stub. Establishes the interface contract on which everything else builds.
2. **Phase 2 (future)**: Python `multiprocessing.Pool` ensemble using the BMI shared library. Per-process isolation makes globals safe without touching `variables.f90`.
3. **Phase 3 (future)**: `swap_ensemble_mod` Fortran orchestrator holding `N` `swap_state_t` columns, with OpenMP parallelism. Requires globals cleanup first.

The architectural choice for Phase 1 is the procedural form (`swap_init`/`swap_run_step`/`swap_close` with explicit state threading) rather than a `swap_model_t` with type-bound procedures. The decision rationale is:

- `swap_state_t` already exists and already aggregates per-instance state.
- The procedural form maps directly onto the BMI C interface (no thin OOP-to-C wrapper layer).
- Phase 3's `swap_ensemble_mod` wraps it cleanly: `type(swap_state_t), allocatable :: columns(:)` with an explicit loop.
- Switching to `swap_model_t` later remains possible if it turns out to be useful — the procedural form is a subset.

## Decision

Introduce `module swap_mod` with three named procedures replacing the magic-integer dispatch. Extract the outer time loop into the caller (`swap_main` and BMI `update()`). Retire `swap_exchange`, `handle_exchange`, `iCaller`, and the `dummy()` test harness in `swap_main`. Add a minimal but spec-compliant BMI C-binding layer (`swap_bmi_mod`) plus a Python hello-world that exercises it via cffi.

New files compile with `-std=f2018 -Wall -Wextra`; legacy sources remain on `-std=legacy`. Per-file scoping is achieved by building the new files in their own static library with target-specific `fortran_args` — the project-wide flag is unchanged.

## Architecture

### What goes

- `swap_exchange` module inside `swap.f90` (lines 1–47): `swap_input`/`swap_output` types — superseded by BMI `get_value`/`set_value`.
- `handle_exchange` subroutine inside `swap.f90` (lines 557–711): all six task-code branches — superseded by BMI lifecycle + variable accessors.
- `iCaller` argument and all `if (iCaller /= 0)` branches in `swap.f90`: concept replaced by BMI (caller is always external from `swap_mod`'s perspective).
- `dummy()` subroutine in `swap_main.f90` (lines 59–113): replaced by proper pFUnit lifecycle test and the Python BMI integration test.
- `src/utils/sharedexchange.f90`: dead-stub `FromSwap`/`ToSwap` subroutines — retired with the rest of the legacy DLL exchange machinery.

### What stays

- `src/core/swap.f90` is **deleted** after its contents migrate to `swap_mod.f90` — no parallel old/new during this arc.
- All physics subroutines (`SoilWater`, `Drainage`, `Temperature`, etc.) and the state migration done by the SS-* arcs — untouched.
- `swap_state_t` and `swap_config_t` — used as-is.
- TimeControl module and `timecontrol_state_t` — unchanged. `TimeControl(1)` stays in `swap_init`; `TimeControl(2)`, `TimeControl(3)`, `TimeControl(9)` stay inside `swap_run_step`. The `associate` block aliasing TC fields (currently `swap.f90:324–341`) moves wholesale into `swap_run_step`, wrapping the entire step body.
- Existing pFUnit suites and regression tests (`hupselbrook`, `grassgrowth`, `oxygenstress`, etc.) — untouched. They are the regression safety net for this refactor.

### New files

```
src/core/swap_mod.f90              ← module form of subroutine swap
src/core/swap_bmi_mod.f90          ← BMI C-binding façade
tests/unit/core/test_swap_mod.pf   ← pFUnit lifecycle smoke test
tests/bmi/swap_bmi.h               ← C header (cffi input)
tests/bmi/hello_swap.py            ← Python BMI hello-world
```

### Modified files

```
src/core/swap_main.f90    ← stripped to ~20 lines: thin driver around swap_mod
meson.build               ← add swap_modern static lib, swap_bmi shared lib, bmi test target
tests/unit/meson.build    ← add test_swap_mod.pf to pf_files
```

## Module interface

### `swap_mod`

```fortran
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
   end subroutine

   subroutine swap_run_step(state, config)
      type(swap_state_t),  intent(inout) :: state
      type(swap_config_t), intent(in)    :: config
   end subroutine

   subroutine swap_close(state, config)
      type(swap_state_t),  intent(inout) :: state
      type(swap_config_t), intent(in)    :: config
   end subroutine

end module swap_mod
```

**`swap_init`** = current `iTask=1` body, verbatim. Loads TOML, calls `Initialize`, `IterTime(1)`, `TimeControl(1)`, allocates state, runs all `*_init` calls. Drops the `if (iCaller /= 0)` branch around the output-file opens — always opens output files.

**`swap_run_step`** = current `iTask=2` body **minus the outer `do while (.not.flrunend)` loop**. Wraps the entire body in the existing `associate` block for TimeControl field aliasing. The inner `do while(fldtreduce)` dt-reduction loop stays where it is. All `if (iCaller /= 0)` branches inside the loop are removed (output and meteo I/O always happen).

**`swap_close`** = current `iTask=3` body, verbatim, with `if (iCaller /= 0)` branch removed.

### `swap_main.f90` after cleanup

```fortran
program swap_main
   use swap_mod,        only: swap_init, swap_run_step, swap_close
   use swap_state_mod,  only: swap_state_t
   use swap_config_mod, only: swap_config_t
   use swap_log,        only: log_init, log_close, LOGLEVEL_INFO
   implicit none

   type(swap_state_t)           :: state
   type(swap_config_t), target  :: config

   call log_init(LOGLEVEL_INFO, 'swap_debug.log')
   call swap_init('swap.toml', state, config)
   do while (.not. state%timecontrol%flRunEnd)
      call swap_run_step(state, config)
   end do
   call swap_close(state, config)
   call log_close()
   write(*,'(a)') ' Swap normal completion!'
   call Exit(100)
end program swap_main
```

The `interface` block disappears (no longer needed — `use swap_mod`). `dummy()` disappears. `CloseTempFil` stays in `swap_main` unchanged for this arc — behavior preservation. It is a candidate for retirement or folding into `swap_close` in a follow-on cleanup, but not this arc.

### `swap_bmi_mod`

Single-column instance. Module-level `state` and `config` with `save`. All procedures `bind(C)`. Returns BMI status codes (`0` = success, non-zero = error).

```fortran
module swap_bmi_mod
   use iso_c_binding,   only: c_char, c_double, c_int, c_null_char, c_ptr
   use swap_mod,        only: swap_init, swap_run_step, swap_close
   use swap_state_mod,  only: swap_state_t
   use swap_config_mod, only: swap_config_t
   implicit none
   private

   type(swap_state_t),          save :: bmi_state
   type(swap_config_t), target, save :: bmi_config

contains

   function initialize(config_file, n) result(rc) bind(C, name='initialize')
      character(kind=c_char), intent(in) :: config_file(*)
      integer(c_int), value,  intent(in) :: n
      integer(c_int)                     :: rc
   end function

   function update() result(rc) bind(C, name='update')
      integer(c_int) :: rc
   end function

   function finalize() result(rc) bind(C, name='finalize')
      integer(c_int) :: rc
   end function

   function get_value_double(var_name, n, dest) result(rc) bind(C, name='get_value_double')
      character(kind=c_char), intent(in)  :: var_name(*)
      integer(c_int), value,  intent(in)  :: n
      real(c_double),         intent(out) :: dest(n)
      integer(c_int)                      :: rc
   end function

   function get_current_time(t) result(rc) bind(C, name='get_current_time')
      real(c_double), intent(out) :: t
      integer(c_int)              :: rc
   end function

end module swap_bmi_mod
```

**Implemented for real:**
- `initialize` → `swap_init`
- `update` → `swap_run_step` (one call advances one `dt`, including inner reduction loop)
- `finalize` → `swap_close`
- `get_current_time` → returns `bmi_state%timecontrol%t1900`
- `get_value_double` → recognises two variable names: `"soil_water_content"` (returns `state%soilwater%theta(1:n)`) and `"actual_evapotranspiration"` (returns `state%soilwater%intr%iqrot`). All other names return `rc=1`.

**Deferred to Phase 2 — stub with correct C signature, body returns 0 or empty:**
- `get_end_time`, `get_time_step`, `get_start_time`
- `set_value_double`, `get_value_int`, `set_value_int`
- `get_var_type`, `get_var_units`, `get_var_nbytes`, `get_var_itemsize`, `get_var_grid`
- `get_grid_type`, `get_grid_rank`, `get_grid_size`, `get_grid_shape`, `get_grid_spacing`, `get_grid_origin`, `get_grid_x`, `get_grid_y`, `get_grid_z`
- `get_input_var_names`, `get_output_var_names`, `get_input_item_count`, `get_output_item_count`, `get_component_name`

Each stub method is marked `! BMI-STUB Phase 2` in the source so the next arc can enumerate them with `grep`.

**Helper:** `c_to_f_string` — a private ~5-line subroutine inside `swap_bmi_mod` that converts a null-terminated C character array to a Fortran `character(len=*)`. Not extracted to a shared utility — it's local to this module.

## BMI binding details

### C header

```c
/* tests/bmi/swap_bmi.h */
#ifndef SWAP_BMI_H
#define SWAP_BMI_H

int initialize(const char *config_file, int n);
int update(void);
int finalize(void);
int get_value_double(const char *var_name, int n, double *dest);
int get_current_time(double *t);

#endif
```

### Python hello-world

```python
# tests/bmi/hello_swap.py
import cffi, pathlib, sys

ffi = cffi.FFI()
header = pathlib.Path(__file__).parent / 'swap_bmi.h'
ffi.cdef(header.read_text().replace('#ifndef SWAP_BMI_H', '').replace(
    '#define SWAP_BMI_H', '').replace('#endif', ''))

so_path = sys.argv[1] if len(sys.argv) > 1 else './libswap_bmi.so'
lib = ffi.dlopen(so_path)

assert lib.initialize(b'swap.toml\0', 0) == 0, 'initialize failed'
assert lib.update() == 0,                       'update failed'

t = ffi.new('double *')
assert lib.get_current_time(t) == 0
print(f'current_time = {t[0]}')

buf = ffi.new('double[500]')
assert lib.get_value_double(b'soil_water_content\0', 500, buf) == 0
print(f'theta[0] = {buf[0]}')
assert 0.0 < buf[0] < 1.0, f'theta[0]={buf[0]} not in (0, 1)'

assert lib.finalize() == 0, 'finalize failed'
print('BMI hello-world: OK')
```

## Build system

```python
# meson.build (additions)

# Modern files: strict standard
swap_modern = static_library('swap_modern',
    sources: ['src/core/swap_mod.f90', 'src/core/swap_bmi_mod.f90'],
    fortran_args: ['-std=f2018', '-Wall', '-Wextra'],
    include_directories: src_inc,
    dependencies: [tomlf_dep]
)

# Existing legacy sources continue under project-wide -std=legacy.
# (No change to project add_project_arguments call.)
# Remove 'src/core/swap.f90' from the existing `sources` list — its body
# has migrated to swap_mod.f90. Do NOT add swap_main.f90 here; it is the
# executable's own source.
swap_legacy = static_library('swap_legacy',
    sources: sources,
    include_directories: src_inc,
    dependencies: [tomlf_dep]
)

# Executable
swap_exe = executable('swap',
    sources: 'src/core/swap_main.f90',
    link_with: [swap_modern, swap_legacy],
    fortran_args: ['-std=f2018'],
    install: true
)

# BMI shared library — no fortran_args, link-only target
swap_bmi_lib = shared_library('swap_bmi',
    link_with: [swap_modern, swap_legacy],
    install: false
)

# Python BMI integration test
python3 = import('python').find_installation('python3')
test('bmi-hello-world',
    python3,
    args: [
        meson.current_source_dir() / 'tests/bmi/hello_swap.py',
        swap_bmi_lib.full_path()
    ],
    workdir: meson.project_source_root() / 'tests/swap-cases/hupselbrook',
    suite: 'bmi',
    depends: [swap_bmi_lib]
)
```

**Cross-standard linkage**: `swap_modern.use`s modules from `swap_legacy` (`swap_state_mod`, `swap_config_mod`, `variables`, etc.). gfortran .mod files produced by legacy-mode compilation are consumed by f2018-mode compilation without issue. The `-std=` flag only affects parsing of the file being compiled, not what it can `use`.

## Testing strategy

### Layer 1 — pFUnit lifecycle smoke test

`tests/unit/core/test_swap_mod.pf` — three tests calling `swap_mod` directly, no C involved. Uses the `hupselbrook` test case via `chdir_helper` (existing pattern from `test_hupselbrook_loads.pf`).

```fortran
@test
subroutine test_swap_init_completes()
   ! After init: flRunEnd is false, t1900 > 0
end subroutine

@test
subroutine test_swap_run_step_advances_time()
   ! t1900 after one step > t1900 before
end subroutine

@test
subroutine test_swap_full_lifecycle()
   ! init → loop until flRunEnd → close, no abort
end subroutine
```

Registered in `tests/unit/meson.build` by adding `'core/test_swap_mod.pf'` to `pf_files` and the corresponding `ADD_TEST_SUITE` to `testSuites.inc`.

### Layer 2 — Python BMI integration test

`tests/bmi/hello_swap.py` registered as a Meson test in the `bmi` suite. Drives initialize → update → get_current_time → get_value_double → finalize. Asserts the data returned by `get_value_double('soil_water_content', ...)` is physically plausible (`0 < theta < 1`).

### Layer 3 — existing regression tests

`hupselbrook`, `grassgrowth`, `oxygenstress`, `salinitystress`, `surfacewater` (and any other case in `tests/regression/`) continue to run via `swap_main`. They exercise the full physics through the new module interface. **No expected-value changes** — outputs must match byte-for-byte (or within existing tolerance) with the pre-refactor baseline.

### Out of scope

- Testing BMI stub methods individually (they return 0; tests would only verify the stub status).
- Full BMI compliance validation against `bmipy.Bmi` abstract base — Phase 2 work.
- Multi-column / parallel correctness — Phase 2/3 work.
- Performance benchmarks — separate arc.

## Verification

Before this arc is considered complete:

1. `meson compile -C builddir` succeeds with no new warnings on `swap_modern` (the f2018-compiled target).
2. `meson test -C builddir --suite unit-pfunit` passes including the three new `test_swap_mod` tests.
3. `meson test -C builddir --suite bmi` passes (the Python hello-world).
4. `meson test -C builddir` full suite passes with no regressions in existing `hupselbrook`/`grassgrowth`/`oxygenstress`/etc.
5. Manual: `./builddir/swap` in `tests/swap-cases/hupselbrook/` runs to completion and writes output identical to the pre-refactor baseline.
6. `grep -rn "iCaller\|swap_exchange\|handle_exchange\|FromSwap\|ToSwap" src/` returns no hits (retirement complete).

## Out of scope (deferred follow-on work)

- **TimeControl modernization**: `timecontrol.f90` retains its `(1, state)` / `(2, state)` / `(3, state)` / `(9, state)` magic-integer dispatch. A future arc converts it to named procedures (`timecontrol_init`, `timecontrol_advance`, `timecontrol_reduce_dt`, `timecontrol_day_end`).
- **Globals cleanup**: 35 files still `use variables`. Phase 3 ensemble work requires this; a separate arc continues the SS-* state migration arc to eliminate the remaining globals.
- **Full BMI variable registry**: design and implement `get_value`/`set_value` for the full set of exchangeable variables (water balance, fluxes, crop state, …). Phase 2.
- **Full BMI metadata methods**: `get_var_type`, `get_grid_*`, etc. — implemented bodies replacing stubs. Phase 2.
- **`swap_ensemble_mod`**: Fortran-level orchestrator for N columns with OpenMP. Phase 3.
- **Project-wide `-std=f2018` upgrade**: incremental, file-by-file. A separate arc.

## References

- ADR 0001 — gfortran-first.
- ADR 0016 — config-passing refactor (typed config threaded through compute subs).
- ADR 0023 — reruns retired (parameter sweeps now external — this is exactly the BMI use case).
- ADR 0041 — TimeControl state migration (the work this arc builds on).
- CSDMS BMI specification v2.0 — `https://bmi.readthedocs.io`.
- Existing `src/core/swap.f90:1–47` (swap_exchange) and `src/core/swap.f90:557–711` (handle_exchange) for the retiring pattern.
