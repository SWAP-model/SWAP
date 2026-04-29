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
