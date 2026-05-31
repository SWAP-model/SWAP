---
title: "ADR 0018 — TTutil FatalERR shim"
date: 2026-05-05
status: accepted
---

# ADR 0018: TTutil `FatalERR` shim

## Context

TTutil's upstream `FatalERR` (`subprojects/ttutil/src/fatalerr.f90`) ends
each error path with a bare `STOP` statement. In Fortran, `STOP` without
an argument exits the process with status **0**. That has two cascading
effects:

1. **The pFUnit harness silently passes failing builds.** When any test
   invokes legacy SWAP code that hits a fatal condition (FOPENG can't
   open a file, RDDATA fails to parse, etc.), the entire test binary
   terminates with exit 0 *before* pFUnit's `finalize()` can run.
   `meson test` sees exit 0 and reports `Ok: 1, Fail: 0` regardless of
   how many `@test` cases actually failed.

2. **The legacy error path is a parallel universe to `error_mod`.**
   The modern code already routes errors through `error_collection_t` /
   `fatalerr_collected` (`src/error/error.f90`), which appends to a
   typed collection, logs through `swap_log`, and aborts via
   `error stop "..."` (exit code 1). Legacy `fatalerr` calls bypass all
   of that.

Symptoms surfaced during SS-5 verification: 6 pre-existing test failures
were silently green, and the test wrapper script `tests/unit/run_pfunit.sh`
had to invoke each suite separately and parse TAP output to recover the
true pass/fail picture.

## Decision

Provide a one-file drop-in replacement for TTutil's `FatalERR`, linked
into every SWAP binary that depends on TTutil. The shim has the exact
signature TTutil's callers expect:

```fortran
SUBROUTINE FatalERR(MODULE, MESSAG)
   use error_mod, only: fatalerr_collected
   implicit none
   character(len=*), intent(in) :: MODULE, MESSAG
   call fatalerr_collected(trim(MODULE), trim(MESSAG))
END SUBROUTINE FatalERR
```

Linker behaviour: object files are processed before libraries, so OUR
`fatalerr_` symbol wins over `libttutil`'s. Verified with
`nm builddir/tests/unit/unit-swap-tests | grep fatalerr_` — exactly one
`T fatalerr_` symbol; no duplicate.

Net effect on every caller (TTutil's own internals — `FOPENG`, `RDINIT`,
`RDDATA` — *and* the 100+ legacy SWAP source files in `src/`):

- Error message text is unchanged (still `' ERROR in <module>: <message>'`
  via `swap_log`).
- The error is appended to `error_mod`'s `global_errors` collection.
- The process exits with status 1 via `error stop` instead of status 0
  via bare `STOP`.
- The "Press <Enter>" stdin block is removed — no more `End of file`
  Fortran runtime error when stdin is closed.

Files added:
- `src/error/fatalerr_shim.f90` (free subroutine, ~10 LoC)

Files modified:
- `meson.build` — adds the shim to the production source list
- `tests/unit/meson.build` — adds the shim to the test source list

Files removed/restored:
- `tests/unit/run_pfunit.sh` — initially removed when the shim landed,
  then **restored** with a narrower scope when a separate failure mode
  surfaced: several parity suites (`hupselbrook_parity`,
  `grassgrowth_parity`, …) invoke the legacy `readswap()` and mutate
  `variables` module globals. State leaks between suites and a later
  suite fatal-errors via RDDATA on what should be valid input. Each
  suite passes when run in isolation; only the cross-suite sequence is
  broken. The wrapper invokes the binary once per suite, which avoids
  the leakage. Until the legacy reader is fully retired (umbrella spec
  SS-11), this isolation is load-bearing for the unit-test gate.
- `tests/unit/meson.build` — `test()` invokes the wrapper, which in
  turn invokes the binary per-suite.

## Consequences

- `pixi run -e test test-pfunit` now correctly reports `Fail: 1` on any
  test failure, including failures that reach the legacy code path.
- CI / pre-commit gates that depend on the unit-test exit code start
  catching failures they previously missed.
- Tests that intentionally trigger `FatalERR` (parity tests invoking
  `readswap()` on case data) now terminate the binary cleanly — no
  blocking stdin read, no exit-0 cover-up. The argv-collision issue
  with parity tests (umbrella spec follow-up: parity test refactor)
  is *separate* from this shim and remains tracked independently.

## Retirement gate — TO BE REMOVED WHEN

This shim exists only because TTutil is still linked in. The umbrella
"legacy reader retirement" spec
(SS-11 closeout) phases TTutil out of the runtime. When that work
completes:

1. No production `.f90` under `src/` calls TTutil routines (`rdinit`,
   `rdsdor`, `rdador`, `rdfdor`, …).
2. `subprojects/ttutil/` is removed from the meson dependency graph for
   the production binary.
3. `subprojects/ttutil/` may still be present for parity-test fixtures,
   in which case the shim stays in the test build only.
4. If parity tests are also retired (the umbrella SS-11 closeout
   discusses converting them to fixture-data comparisons), TTutil is
   removed entirely and this shim file (`src/error/fatalerr_shim.f90`)
   is deleted in the same change. Both meson source-list entries are
   reverted.

Until then, this shim is load-bearing for the test harness and the
production exit-code correctness.

## Alternatives considered

**A. Mechanical replace of `call fatalerr` → `call fatalerr_collected`
across `src/`** (~100 sites). Pros: explicit, no symbol shadowing.
Cons: doesn't reach TTutil's *internal* calls (FOPENG → fatalerr,
RDINIT → fatalerr, RDDATA → fatalerr), which run inside `libttutil`
and have the same silent-exit-0 problem. Would require shimming
anyway. Combine with this ADR if/when convenient as a no-behavior-
change cleanup; no rush.

**B. Patch TTutil's source.** Pros: fixes the bug at the source.
Cons: TTutil is a vendored subproject; patches would have to be
re-applied on every update. Forking TTutil is heavier than shadowing
one symbol.

**C. Wrap the test binary in a script** (`tests/unit/run_pfunit.sh`).
The shim solves the silent-exit-0 root cause, but the wrapper retains
load-bearing value for **per-suite isolation** — see the
"Files removed/restored" note above. Until the legacy `readswap()`
retires (umbrella spec SS-11), some parity suites cannot share a
process without state corruption. The wrapper trades a small amount
of complexity for a working unit-test gate.
