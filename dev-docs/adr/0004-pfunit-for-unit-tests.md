# ADR 0004 — pFUnit for Fortran unit tests

Status: accepted (2026-04-23, during rescue Phase 2; wiring finalized in Phase 1)

## Context

Options considered at rescue start:

- **pFUnit 4** — mature, CMake-based, preprocessor-generated Fortran test drivers. De-facto standard in Fortran scientific computing (NASA, GFDL, ECMWF use it).
- **test-drive** (Fortran-Lang) — lightweight, pure Fortran, no preprocessor. Easier to vendor but less feature-complete (no parameterized tests, no fixtures in the pFUnit sense).
- **Handroll** — a few dozen lines of `call assert_equal(...)` helpers. Maximum control, minimum ergonomics.

At the rescue baseline, `tests/unit/meson.build` was already wired for pFUnit 4.15, though all suites were orphaned against drift-era modules (Phase 1 removed them). The scaffolding (generator, driver, `testSuites.inc`) remained functional.

## Decision

Use pFUnit 4 for Fortran unit tests throughout the rescue and beyond.

## Consequences

- **Positive**: ecosystem standard; conventions transfer to and from other scientific Fortran projects.
- **Positive**: rich fixtures and parameterized tests — useful for state lifecycles and physics routines across parameter combinations.
- **Positive**: wiring already in place; adding a new test is `.pf` file + `ADD_TEST_SUITE()` line.
- **Negative**: pFUnit is CMake-based, not meson-native. The rescue lives with a gitlink checkout at `tests/pFUnit/` (see `docs/build-and-test.md`). Proper vendoring is deferred.
- **Negative**: preprocessor adds Python compile-time dependency; pixi manages this, so net friction is zero.

## Alternative if circumstances change

test-drive as fallback. Conversion cost: rewriting `.pf` suites in test-drive's Fortran-native API.
