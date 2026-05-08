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
- [ADR 0017 — Sibling-reader dispatch from per-rotation cache](0017-sibling-reader-dispatch-from-cache.html) — Multiple legacy readers may open the same `.crp` file (e.g. `readcropfixed` + `readarablelandgerm`). When a sibling reader is discovered (smoke-test signal: `FOPENG: File does not exist`), apply the same cache-driven dispatch at its callsite. Audit checklist updated to grep all `'.crp'` opens before each port phase.
- [ADR 0020 — Call-site gating convention](0020-call-site-gating-convention.html) — All optional subsystems are gated at the call site by an `flX` flag; subsystem implementations assume the flag is true. Brings tillage and SSDI in line with the existing pattern (flCropNut, flMacroPore, flSurfaceWater, …).
- [ADR 0021 — Tillage parameters ported to [soil.tillage] TOML block](0021-tillage-toml-port.html) — Tillage subsystem migrated from TTutil-based `swap.swp` reads to a typed `[soil.tillage]` sub-table; `Read_Tillage` retired.
- [ADR 0022 — SSDI parameters ported to [irrigation.ssdi] TOML block](0022-ssdi-toml-port.html) — SSDI subsystem migrated from TTutil-based `swap.swp` + `ssdi_file` reads to a typed `[irrigation.ssdi]` sub-table with schedule discriminator; `read_ssdi_input` and `checkdate.f90` retired.
- [ADR 0023 — TTutil retired from the build](0023-ttutil-retirement.html) — TTutil retired from the SWAP build (subproject + utility functions + reruns); `src/core/dtutil.f90` provides canonical native-Fortran date/string utilities; `fatalerr_shim.f90` renamed to canonical `fatalerr.f90`.
- [ADR 0024 — `dtutil.f90` as TTutil-API compatibility shim](0024-dtutil-compatibility-shim.html) — Documents the Phase E choice to introduce `src/core/dtutil.f90` as a same-signature shim for 11 TTutil utility functions, deferring ~75 call-site edits. Future direction: hoist utility calls out of physics layer; physics receives parsed/validated inputs only.
- [ADR 0025 — [nutrients] N1: crop-side nutrient adapter](0025-nutrients-N1-crop-side-adapter.html) — Wires wofost_nutrient_t typed config to legacy variables globals via apply_cropwofost_nutrient adapter; validator stub-error replaced with input checks. Runtime gate at tillage.f90:73 stays in place pending N2/N3.
- [ADR 0026 — [nutrients] N2a: soil-side initial state + SorpCoef](0026-nutrients-N2a-soil-side-initial-state.html) — Top-level `[nutrients]` TOML block carries initial pool concentrations (12 reals) and SorpCoef. Adapter populates legacy globals unconditionally; SorpCoef gets a deterministic default fixing a genuine uninitialised-variable bug. SoilManagement(7) `_nut.end` dump retired (broken without .snp template files).
- [ADR 0027 — [nutrients] N2b: timed amendments via CSV companion](0027-nutrients-N2b-timed-amendments.html) — Optional `events_file` on `[nutrients]` references a 4-column CSV (`date,material,amount_kgha,volat_fraction`); adapter stages, sorts, groups same-day dosages, populates legacy globals. `volat_fraction` range tightened to [0, 1] (legacy typo was [0, 500000]). Runtime path unchanged.
- [ADR 0028 — [nutrients] N3: runtime activation](0028-nutrients-N3-runtime-activation.html) — Drives legacy global `flCropNut` per-rotation from `cfg%nutrient%flcropnut` inside `cropwofost_init_from_config`; lifts validator stub at `cropwofost_init.f90:85` and runtime stub at `tillage.f90:73`. Closes the `[nutrients]` umbrella (N1+N2a+N2b+N3). Regression fixture deferred — no legacy nutrient-enabled `.crp` exists to mirror.
