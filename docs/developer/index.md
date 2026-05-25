---
title: Developer
---

# Developer

For contributors to the SWAP modernization.

## Orientation

- **[Architecture](architecture.html)** — three-phase control flow, state
  aggregation, module boundaries.
- **[Contributing](contributing.html)** — branch model, commit conventions,
  verification gates.
- **[Build and test](build-and-test.html)** — building locally, fast/full
  test split, regression harness.
- **[Code style](code-style.html)** — Fortran 2008 conventions matching
  `.fprettify.rc`.

## Subsystems

- **[State management](state-management.html)** — `config_t` / `initial_t`
  / `state_t` lifecycle and the ASSOCIATE pattern.
- **[Error handling](error-handling.html)** — error collection, calling
  convention, abort checkpoint.
- **[Logging](logging.html)** — `swap_log` facility, levels, format.
- **[Validation](validation.html)** — primitive checks, section validators,
  how to add new ones.

## Operations

- **[Branches](branches.html)** — branch model and naming.
- **[Dependency management](dependency-management.html)** — subprojects,
  pixi deps, pFUnit, how to bump.
