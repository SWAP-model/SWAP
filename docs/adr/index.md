---
title: Architecture Decision Records
---

# Architecture Decision Records

Records of significant architecture and process decisions made during the SWAP modernization.

- [ADR 0001 — gfortran-first](0001-gfortran-first.html) — Only gfortran is supported; non-GCC Fortran compilers are out of scope.
- [ADR 0002 — single builddir](0002-single-builddir.html) — Single out-of-tree build directory; no per-variant subdirs.
- [ADR 0003 — aggregator state over compartments for now](0003-aggregator-state-over-compartments-for-now.html) — Initial state design uses aggregator `swap_state_t` rather than per-compartment state.
- [ADR 0004 — pFUnit for unit tests](0004-pfunit-for-unit-tests.html) — pFUnit framework for lightweight Fortran unit tests.
- [ADR 0005 — licensing audit](0005-licensing-audit.html) — Repository licensed as GPL v2 to match upstream SWAP 4.2.0.
- [ADR 0006 — coverage tracked, not gated](0006-coverage-tracked-not-gated.html) — Coverage is a reference metric, not a CI gate.
- [ADR 0007 — config parse-validate-finalize-adapter pipeline](0007-config-validate-finalize-pipeline.html) — Four-stage input pipeline: parse, validate, finalize, adapter.
- [ADR 0008 — error collection over fatalerr](0008-error-collection-over-fatalerr.html) — Error collection for full-pass validation; one abort checkpoint after pipeline.
- [ADR 0009 — discontinue non-CSV output formats](0009-discontinue-non-csv-outputs.html) — Retire 18 legacy output switches; CSV becomes the only supported output path.
- [ADR 0010 — macropore module deferral](0010-macropore-deferral.html) — Macropore stays in legacy code; new TOML pipeline does not wire it. Schema kept as orphan infrastructure for future macropore work.
- [ADR 0011 — macropore case exclusion from regression](0011-macropore-exclusion-from-regression.html) — Case 3 dropped from check-full; ADR 0010's runtime-fallback clause superseded; Phase 4f's strangler-fig becomes a clean cut.
- [ADR 0012 — CSV companion input files](0012-csv-companion-input-files.html) — `read_csv_table` as the canonical reader for all tabular inputs in the TOML pathway; ISO date support, header validation, `error_collection_t` errors.
- [ADR 0013 — CSV meteorology input](0013-csv-meteorology-input.html) — TOML pathway reads daily meteo from a `date,rad,tmin,…,wet` CSV; legacy `.met`/`.YYY` files unchanged for the ASCII pathway.
- [ADR 0014 — readmeteo.f90 TTutil phase-out](0014-readmeteo-phaseout.html) — Three-step plan to delete all TTutil branches from `readmeteo.f90`: implement sub-daily CSV, delete dead TTutil code, clean up dead variables.
- [ADR 0015 — Strangler narrow-scope stub-errors](0015-strangler-narrow-scope-stub-errors.html) — Phase 4f legacy-reader ports default to narrow-scope coverage with validator-level stub-errors for branches no test case exercises. Legacy reader stays alive as parity-test fixture.
- [ADR 0016 — Per-rotation crop config cache](0016-per-rotation-crop-config-cache.html) — `.crp` port (Phases 1-4) loads all rotation crop content into parallel arrays on `crop_config_t` at config-load time and dispatches into the cache from `cropgrowth.f90`'s per-rotation init. Future direction: pass typed config + state to computation subs as explicit arguments, eliminating the legacy module-global mutation pattern.
