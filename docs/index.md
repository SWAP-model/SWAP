---
title: SWAP documentation
author: SWAP Team
---

# SWAP — Soil Water Atmosphere Plant

SWAP is a one-dimensional vertical simulation model for transport processes (water, heat, solutes) in the soil–plant–atmosphere continuum at field scale. This repository holds the modernization of SWAP 4.2.0 under active rescue-and-stabilize work.

## Start here

1. [Architecture overview](architecture.html) — three-phase control flow, state aggregation, module boundaries.
2. [Build and test](build-and-test.html) — build, fast/full test split, regression harness, fixture policy.
3. [Contributing](contributing.html) — branch model, commit conventions, verification gates.

## Reference

- [State management](state-management.html) — `config_t` / `initial_t` / `state_t`, ASSOCIATE, lifecycle.
- [Configuration schema](configuration-schema.html) — TOML input reference.
- [Error handling](error-handling.html) — error collection, calling convention, abort checkpoint.
- [Validation](validation.html) — primitive checks, section validators, aggregate rules, testing.
- [Coverage baseline](coverage-baseline.html) — Phase 3 line/branch coverage numbers and known gaps.
- [Dependency management](dependency-management.html) — subprojects, pixi deps, pFUnit, how to bump.
- [Code style](code-style.html) — Fortran 2008 conventions, names, intent; matched by `.fprettify.rc`.
- [Logging](logging.html) — `swap_log` facility, levels, initialization, format, thread-safety constraints.
- [Architecture decision records](adr/) — the non-obvious choices and why:
  - [ADR 0001 — gfortran-first](adr/0001-gfortran-first.html)
  - [ADR 0002 — single builddir](adr/0002-single-builddir.html)
  - [ADR 0003 — aggregator state over compartments for now](adr/0003-aggregator-state-over-compartments-for-now.html)
  - [ADR 0004 — pFUnit for unit tests](adr/0004-pfunit-for-unit-tests.html)
  - [ADR 0005 — licensing audit](adr/0005-licensing-audit.html)
  - [ADR 0006 — coverage tracked, not gated](adr/0006-coverage-tracked-not-gated.html)

## Current status

The modernization is in a rescue-and-stabilize workflow; see `docs/superpowers/specs/` and `docs/superpowers/plans/` for the internal planning. Status of the rescue is tracked in git tags named `rescue/phase-N-*`.

## API reference

FORD-generated API reference at [api/index.html](api/index.html). Run `pixi run -e docs docs-build` to regenerate.
