---
title: "Small Fix Note — Suppress gfortran IEEE FP summary at termination"
date: 2026-05-13
status: deferred
relates-to: SS-DRV Phase 1 cutover (commit effde1b)
priority: low (cosmetic — no numerical impact)
---

# Small Fix — gfortran IEEE FP summary on `stop 100`

## Symptom

Direct invocation of `./builddir/swap` (e.g. via `tests/swap-cases/run_case.sh -c hupselbrook`) prints at termination:

```
Note: The following floating-point exceptions are signalling:
      IEEE_INVALID_FLAG IEEE_DIVIDE_BY_ZERO IEEE_DENORMAL
STOP 100
```

Regression tests do NOT surface this because `test_output_regression.py` captures stderr but only prints it on a non-100 exit code.

## Root cause

SS-DRV Phase 1 Task 9 (commit `effde1b`) changed `call Exit(100)` → `stop 100` in
`src/core/swap_main.f90`. `call Exit` is a gfortran extension rejected by `-std=f2018`,
which was newly enforced on the executable target in that task.

- `call Exit(100)` invokes libc `_exit(2)` — terminates immediately, bypasses Fortran cleanup, no IEEE flag report.
- `stop 100` uses Fortran's standard termination path. gfortran's runtime inspects any IEEE FP exception flags raised during the run and prints a summary to stderr.

**The flags themselves are pre-existing.** They have been raised on every SWAP run for years (typical sources: `oxygenstress`, evapotranspiration when LAI=0, denormal-producing numerical kernels). The model produces byte-for-byte identical output regardless — guard logic catches the cases where the flag fires before the result is consumed. `check-full` confirmed 5/5 byte-for-byte parity through the SS-DRV cutover and SS-TCM arc.

## Fix (when picked up)

**Preferred — Fix A: libc exit via `bind(C)`**

Replace `stop 100` with a `bind(C)` wrapper around libc's `exit(3)`. Standard-compliant (no gfortran extension) and restores the pre-Task-9 termination behavior:

```fortran
! In swap_main.f90
use iso_c_binding, only: c_int
interface
   subroutine c_exit(status) bind(C, name='exit')
      import :: c_int
      integer(c_int), value :: status
   end subroutine c_exit
end interface
...
call c_exit(100_c_int)
```

Cost: ~8 lines in `swap_main.f90`. No regression risk (regression already passes; this only changes the termination path).

**Alternative — Fix B: project-wide `-ffpe-summary=none`**

Add to `add_project_arguments` in `meson.build`:

```python
gfortran_flags += ['-ffpe-summary=none']
```

Suppresses the summary for ALL builds. Simpler (1 line) but coarser — also affects debug builds.

**Proper fix — Fix C (separate arc): hunt down the FP edge cases**

Run with `-ffpe-trap=invalid,zero,denormal` to get tracebacks at the exact firing sites, then add explicit guards or use IEEE intrinsics. Self-documenting but a 1–2 day hygiene arc. Out of scope here.

## When to pick this up

After BMI Phase 2 ships. The cosmetic message is harmless and does not affect any consumer:
- Regression: passes (exit code 100, output byte-for-byte).
- pFUnit: tests do not parse stderr.
- imod_coupler / pyswap (future): consume the shared library, not the executable — gfortran's runtime summary only fires at executable termination.

Half-hour fix; can be batched with any small cleanup arc.
