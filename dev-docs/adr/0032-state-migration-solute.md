---
title: "ADR 0032 — Solute state-type migration"
date: 2026-05-10
status: accepted
---

# ADR 0032: Solute state-type migration

**Status:** accepted (Phases 0 + 1 + 2 complete)
**Date:** 2026-05-10
**Migration #:** 3 of N
**Branch:** `refactor/solute-state` — ready to merge to `development`

## Context

The surfacewater (ADR 0030) and drainage (ADR 0031) state-type migrations established the playbook: discovery → design → plan → execute, with state-type definition, argument-threaded `swap_state_t`, ASSOCIATE in compute, dual-write transitional pattern, and check-full byte-identical at every commit.

Solute is migration #3. The discovery doc surfaced four distinctive concerns that shaped the design:

- **AgeTracer is dead code.** `flAgeTracer` is never set to `.true.` anywhere in `src/`. The feature has been dormant. AgeTracer dual-uses `cml` (the solute concentration array) as age storage by overwriting `cml(i)` at the end of each `AgeTracer(task=2)` call.
- **Major physics config gap.** 14 fields consumed by `solute.f90` were absent from `solute_config_t` and the adapter. TOML-driven runs silently used zero defaults for `cref`, `kf`, `frexp`, `gampar`, `decpot`, `fdepth`, `ddif`, `kfsat`, `poros`, `daquif`, `decsat`, `cseeptab`, `swbr`, `cpre`. A real correctness bug for any TOML run with `swsolu=1`.
- **Two-stage `cml` seeding** (config CSV → legacy global → typed state at run init) needed preservation.
- **`ArMpSs` working-buffer co-write** by 3 subsystems (solute, soilhydraulics, boundtop) — not solute-owned, kept as global.

## Decision

Three-phase migration on a single branch:

- **Phase 0 — Physics config gap.** Promote 14 missing physics fields into `solute_config_t` + TOML reader + adapter. Patches the silent-zero-defaults correctness gap. No regression coverage (none of the 5 cases activates `swsolu=1`).
- **Phase 1 — State type, AgeTracer extraction, threading, dual-write.** Define `solute_state_t` (25 fields). Aggregate under `swap_state_t`. Extract `AgeTracer` from `solute.f90` into a new `src/solute/agetracer.f90` module with stub-error guard (preserves dead code for future reactivation). Gate `outage` behind `flAgeTracer`. Solute compute dual-writes state and globals. Output reads switch to typed state.
- **Phase 2 — Cross-subsystem reader migration + cleanup.** Migrate `rootextraction` and `irrigation` (compute readers of `cml`) to typed state. Introduce `solute_init(state)` to seed `state%solute%cml/cmsy` from config-time globals at run init. Drop solute compute's dual-write of legacy globals. Comment out 18 of 27 solute-owned global declarations. ADR 0032 finalize.

### Phase 0 outcomes

- 14 fields added to `solute_config_t`: 10 scalars (`cref`, `cpre`, `ddif`, `frexp`, `gampar`, `daquif`, `kfsat`, `decsat`, `poros`, `swbr`) + 3 per-soil-layer arrays (`kf`, `decpot`, `fdepth`) + 1 2D table (`cseeptab`).
- TOML reader extended with 14 new field reads (using existing private helpers).
- `apply_solute` adapter populates the 14 legacy globals from typed config. `cseeptab` flattens 2D `(rows, 2)` typed config into the legacy interleaved 1D `(2*k-1)=time, 2*k=concentration` layout (verified from `solute.f90:107`'s `afgen(cseeptab, mabbc*2, time)` call pattern).
- 14 new pFUnit tests covering defaults, ranges, and array allocation.

### Phase 1 outcomes

- `src/state/solute_state.f90` defines `solute_state_t` with 25 fields:
  - 2 allocatable per-node arrays (`cml(:)`, `cmsy(:)`) sized from `numnod`.
  - 7 step-level scalars (`csurf`, `cpond`, `cdrain`, `cseep`, `dtsolu`, `isqbot`, `isqtop`).
  - 8 cumulative balance scalars (`samini`, `sampro`, `samcra`, `solbal`, `dectot`, `imdectot`, `rottot`, `imrottot`).
  - 8 cumulative source/sink fluxes (`sqprec`, `imsqprec`, `sqirrig`, `imsqirrig`, `sqbot`, `imsqbot`, `sqdra`, `imsqdra`, `sqsur`, `sqrap`).
- 12 AgeTracer-specific globals stay declared in `variables.f90` (commented as "future agetracer_state_t targets" per design D5).
- `ArMpSs` stays as a shared global — co-written by solute + soilhydraulics + boundtop. Macropore migration sorts ownership later.
- `AgeTracer` extracted: 244 lines moved from `src/solute/solute.f90` to a new `src/solute/agetracer.f90` declaring `module agetracer_mod`. Stub-error guard at the top fires if `flAgeTracer = .true.` (currently never assigned). Body preserved verbatim. Reactivation checklist documented in the file header.
- `outage` gated with `if (.not. flAgeTracer) return` (defense-in-depth on top of existing call-site gating in `swap.f90`).
- Solute compute writes 25 fields to both `state%solute%*` and legacy globals. Output routines (`outvap`, `outend`, `outage`, `outbal`, `outsba`, `set_values`, `fill_values`, `csv_out_tz`) switched to read from typed state.
- 5 new pFUnit tests for `solute_state_t` (defaults, allocation, independent instances).

### Phase 2 outcomes

- Cross-subsystem compute readers migrated:
  - `rootextraction.f90:RootExtraction` reads `state%solute%cml` (Maas-Hoffman salinity stress + osmotic head matric flux).
  - `irrigation.f90` reads `state%solute%cml(nodsen)` for concentration threshold.
  - Signature propagation chain caught: `RootExtraction` → `JongvanLier` → `JongvanLierLoop` → `MatricFlux`. `MatricFlux` accepts `state` as `optional intent(in)` — task=1 init callers in cropgrowth/wofost/grass don't pass it; task=2 callers in rootextraction do.
- `solute_init(state)` introduced as the lifecycle init for the typed state. Called from `swap_main` immediately after `drainage_init`. Allocates `state%solute%cml/cmsy` from `numnod` and seeds them from the legacy `cml`/`cmsy` globals (which `config_to_variables` writes from CSV at config time). Preserves the two-stage `cml` seeding pattern (discovery hazard #7).
- Solute compute's dual-write of 25 owned fields dropped — `state%solute%*` is the sole authoritative location. ASSOCIATE blocks in `solute(task=1)` (4 fields) and `solute(task=2)` (25 fields) replace per-loop `<global> = ...` writes with direct state writes; per-node array slice mirrors no longer needed.
- 18 solute-owned globals commented out in `variables.f90` with provenance markers:
  `cpond`, `cseep`, `csurf`, `dectot`, `imdectot`, `imrottot`, `sqprec`, `imsqprec`, `sqirrig`, `imsqirrig`, `sqbot`, `imsqbot`, `imsqdra`, `sqsur`, `sqrap`, `samcra`, `sampro`, `solbal`.
- 9 globals retained:
  - `cml(macp)`, `cmsy(macp)` — config-time-seeded by `apply_solute`; `solute_init` reads them at run init to populate state. Removing breaks the seeding bridge.
  - `cdrain`, `samini`, `isqbot`, `isqtop`, `dtsolu`, `rottot`, `sqdra` — referenced by `agetracer.f90`'s dead-code body via wildcard `use Variables`. Tagged `[AgeTracer dead-code dep]`. Reactivation arc removes them when `agetracer_state_t` is defined.

## Architectural learnings (Phase 0 + 1 + 2)

- **Config schema gaps surface during state-migration discovery.** Phase 0's discovery of 14 missing physics fields was orthogonal to state migration but blocked design (TOML runs would have used wrong physics). Future subsystem migrations should grep for compute-reads of borrowed globals NOT in the typed config — those are similar gaps. The pattern: any field read by `<subsystem>.f90` but absent from `<subsystem>_config_t` is a candidate.
- **Stub-errored module extraction is a clean way to preserve dead code.** AgeTracer's 244-line body was kept verbatim in a new module with a documented reactivation checklist. The cost is one extra module file; the benefit is the solute compute path is uncluttered and the dead-code path is unambiguously inert.
- **`use Variables` wildcard imports are the hardest part of clean-up.** AgeTracer kept its wildcard import (option D from Task 9) because narrowing it to a precise `only:` list would have required 50+ names. Phase 2 Task 9 left 7 solute-owned globals declared in `variables.f90` solely because of this. Reactivation of AgeTracer must include the wildcard-narrowing as part of its scope.
- **`MatricFlux`-style optional state args** are useful when a routine is called from many call chains, some pre-migration and some post-migration. Avoids signature breakage on init-time callers that don't need state.
- **The two-stage `cml` seeding** (config global → typed state via init) shows that `_init` routines bridging config→state are a recurring pattern. The same shape was needed for `drainage_init` and now `solute_init`. Future subsystems will likely use the same.

## References

- Predecessors: ADR 0030 (surfacewater pilot), ADR 0031 (drainage)
- Phase 0: `86fb4f5` (config + reader + adapter + tests)
- Phase 1: `4c7b6d8` (state type + aggregator) → `0678200` (AgeTracer extract) → `8460dd4` (dual-write) → `e3b7e58` (output reads)
- Phase 2: `598e551` (rootextraction/irrigation) → `cc3c575` (solute_init) → `378ca96` (drop dual-write) → `ed5cf93` (comment out globals)
