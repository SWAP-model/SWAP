# Atmosphere State-Type Migration — Design

**Date:** 2026-05-11
**Status:** accepted (pending implementation)
**Migration #:** 7 of N — third of the four coupling-surface arcs (boundary → crop-uptake → **atmosphere** → soil-water core)
**Branch:** `refactor/atmosphere-state`
**Discovery:** `docs/superpowers/specs/2026-05-11-state-migration-atmosphere-discovery.md`
**Predecessor:** `docs/superpowers/specs/2026-05-10-state-migration-crop-uptake-design.md` / ADR 0036
**Successor ADR:** 0037

## Goal

Carve ~33 atmosphere-owned globals (`peva`, `ptra`, `atmdem`, `empreva`, `grai`, `nraida`, `ssnow`, `melt`, and 25 others) out of `variables.f90` into a new `atmosphere_state_t` record, aggregated under `swap_state_t` as `state%atmosphere`. The type debuts the ADR 0033 cohort pattern for a freshly-introduced state record: two nested sub-records (`atmosphere_intermediate_t` and `atmosphere_cumulative_t`) replace 3 scattered reset blocks and resolve the `igrai`/`inrai` double-reset hazard. The meteoday/meteodt files are included as full home-tree members (scope-correction from discovery exclusion plan); a future atmosphere REFACTOR arc handles their structural reorganisation. Every commit must be byte-identical on check-full.

## Out of scope

- The future atmosphere **refactor** arc (restructure meteoday.f90 / meteodt.f90 orchestration) — state-migration only here.
- Soil-water-core cumulatives (`cqrot`, `iqrot`, `inqrot`, `iqredXXX*`) — separate arc (#8).
- The 7 non-atmosphere globals read by `reduceva` (`swredu`, `fldaystart`, `cofred`, `dt`, `rsigni`, `nird`, `pond`) — stay as legacy reads; only the 5 atmosphere-owned fields (`empreva`, `ldwet`, `peva`, `spev`, `saev`) are migrated in this arc.
- `nird`/`gird` ownership — irrigation-owned; deferred to a future irrigation arc.
- `pond` migration — boundary D5 deferral persists; `reduceva` reads the legacy global until the soil-water-core arc migrates `pond` to `state%soilwater`.
- `mfluxtable` — already retired in crop-uptake C-2.6 (ADR 0036). No interaction.

## Context

The atmosphere subsystem is the third coupling-surface arc in the soil-water 4-arc decomposition. It is structurally different from the two preceding arcs (boundary: 12 instantaneous scalars; crop-uptake: 22 instantaneous fields including per-node arrays): atmosphere owns fields across **five reset cadences** — 11 instantaneous, 9 per-day, 8 intermediate (`flzerointr`-gated), 10 cumulative (`flzerocumu`-gated), and 3 per-event (tillage-reset). This is the first arc where the ADR 0033 cohort pattern applies to a freshly-introduced state record rather than retrofitted onto an existing flat type. The 10 cumulatives are currently reset across three different files (`meteoday.f90`, `soilhydraulics.f90`, `snow.f90`); the two intermediate cohorts add a fourth ownership split (`waterbalance.f90`). Packaging both cohorts under `atmosphere_state_t` collapses that scatter to two type-bound `reset()` calls and positions the soil-water-core arc to reuse the pattern for `cqrot`/`iqrot`.

The scope correction vs the discovery: `meteoday.f90` and `meteodt.f90` move from "excluded co-writers requiring additive dual-writes" to **full home routines**. They are the canonical writers of `peva`, `ptra`, and `atmdem` (the highest-traffic atmosphere outputs, read by rootextraction × 25 sites and cropgrowth × 9 sites). Including them doubles the home-tree LoC to ~3 000 but eliminates the "excluded files that must still be touched" awkwardness and resolves discovery hazard #1. The future refactor arc is free to restructure the meteo orchestration without invalidating the state writes landed here.

## Decisions

- **D1. Pilot scope — full atmosphere home tree (7 files).** `atmosphere_constants.f90` (empty stub), `et.f90`, `interception.f90`, `meteoday.f90`, `meteodt.f90`, `precipitation.f90`, `snow.f90`. The discovery's exclusion of meteoday and meteodt is REVERSED per user directive. Dual-writes in those files are standard intra-home dual-writes, not "additive patches to excluded files."

- **D2. State-type shape — TWO COHORTS + flat top-level scalars.** `atmosphere_state_t` debuts the ADR 0033 cohort pattern at type-creation (first arc not to retrofit an existing flat type). Layout:

  ```fortran
  type :: atmosphere_state_t
     ! Flat instantaneous (11 scalars):
     real(real64) :: peva, ptra, empreva, melt, subl, slw
     real(real64) :: ssnow, snowinco, graidt, nraidt, aintcdt
     ! Flat per-day (9 scalars):
     real(real64) :: grai, nraida, atmdem, pevaday, ptraday
     real(real64) :: gsnow, snrai, fprecnosnow, sicact
     ! Flat per-event / tillage-reset (3 scalars):
     real(real64) :: ldwet, spev, saev
     ! Cohort sub-records:
     type(atmosphere_intermediate_t) :: intr
     type(atmosphere_cumulative_t)   :: cumu
  end type atmosphere_state_t
  ```

  `atmosphere_intermediate_t` holds 8 fields: `igrai`, `inrai`, `ipeva`, `iptra`, `ievap`, `igsnow`, `isubl`, `isnrai`. `atmosphere_cumulative_t` holds 10 fields: `cgrai`, `cnrai`, `caintc`, `cpeva`, `cptra`, `cevap`, `cgsnow`, `csubl`, `csnrai`, `cmelt`. Each cohort has one type-bound `reset()` procedure. All 40 fields are scalars; no per-node arrays. ASSOCIATE prefix `at_` used in large compute bodies.

- **D3. Aggregation — `swap_state.f90` gains `type(atmosphere_state_t) :: atmosphere`.** `swap_state_mod` adds `use atmosphere_state_mod` and a new component. No change to the existing `soilwater_state_t` or `soilwater_init`.

- **D4. `atmosphere_init(state%atmosphere)` — new init routine; called at `swap.f90:188`.** Zeros all 22 top-level scalars and initialises both cohort sub-records. No `numnod`/`nlay` parameters — all atmosphere fields are scalars. Insertion point: immediately after `soilwater_init(state%soilwater, numnod, numlay)` at `swap.f90:188`, before `DoTillage(1)` at line 190.

- **D5. Cohort reset call sites — replace 3 scattered reset blocks with type-bound calls.** The three legacy reset sites (`meteoday.f90:412–423` for cgrai/cnrai/caintc + igrai/inrai, `soilhydraulics.f90:1112–1114` + `1140–1142` for ipeva/iptra/ievap + cpeva/cptra/cevap, `snow.f90:87–98` for igsnow/isubl/isnrai + cgsnow/csubl/csnrai/cmelt) are each replaced with:
  - `if (flzerointr) call state%atmosphere%intr%reset()` at the canonical intr-reset site (soilhydraulics.f90)
  - `if (flzerocumu) call state%atmosphere%cumu%reset()` at the canonical cumu-reset site (soilhydraulics.f90 / SoilWater(1) block)
  - The snow.f90 and meteoday.f90 inline reset blocks are removed once the cohort reset owns the fields.

- **D6. `igrai`/`inrai` double-reset resolution (discovery hazard #3).** Both `waterbalance.f90:388–389` and `meteoday.f90:414–415` currently zero `igrai`/`inrai` under `flzerointr`. With the cohort-owned `reset()`, atmosphere is the single owner. The duplicate reset writes in `waterbalance.f90` lines 388–389 and the `meteoday.f90:ResetMetFlx` block at lines 414–415 are removed. The single `call state%atmosphere%intr%reset()` at the canonical site (soilhydraulics.f90 or the SoilWater call) replaces all three.

- **D7. `reduceva` globals — address only the 5 atmosphere-owned reads/writes.** Discovery hazard #6: `reduceva` reads 12 globals; 5 are atmosphere-owned (`empreva`, `ldwet`, `peva`, `spev`, `saev`). These 5 are migrated in this arc. The other 7 (`swredu`, `fldaystart`, `cofred`, `dt`, `rsigni`, `nird`, `pond`) stay as legacy reads. New signature: `subroutine reduceva(task, nrai, state)` with `state` as `intent(inout)`. Callers: `meteoday.f90:808`, `meteodt.f90:358`, `meteodt.f90:448`.

- **D8. `pond` read in `et.f90:reduceva` (boundary D5 deferral).** Keep as legacy global read for now. Mark with: `! [SS-ATM] reads legacy pond — soil-water-core arc migrates`. The `state` arg added in D7 enables future cutover without signature surgery.

- **D9. `snow.f90` signature promotion — `optional intent(in)` → `intent(inout)`.** Snow writes `ssnow`, `melt`, `subl`, `empreva`, `peva` (zero-out paths). The current `optional intent(in)` state arg cannot cover write-backs. After promotion, callers `swap.f90:216` and `swap.f90:298` remain unchanged (both already pass full `state`).

- **D10. `nird`/`gird` ownership — defer.** `interception.f90:DivIntercep` co-writes `nird` by partitioning; `gird` is set by `irrigation.f90`. The irrigation subsystem is the canonical owner. This arc reads `nird`/`gird` in DivIntercep but does not claim ownership. Deferred to a future irrigation arc.

- **D11. `precipitation.f90:127` `ssnow = 0.0d0` mutation — keep verbatim.** After migration becomes `state%atmosphere%ssnow = 0.0_real64`. The legacy code comment ("Note: This modifies a state variable — consider refactoring") is preserved. No structural fix in this arc.

- **D12. Phase 0 — small; verify `swsublim` typed-config + `spev`/`saev` adapter seeding.** `swsublim` is used in `snow.f90:104` but its config-coverage status is unconfirmed. `spev` and `saev` initial values come from `soil_config.initial`; check adapter seeding. Expected 1–2 commits at most; likely 0 if both are already covered.

- **D13. Compile-driven Phase 2.7 expectation — 3–6 fixup commits.** Per playbook lesson #5 and crop-uptake's 2-iteration maturity signal, expect 3–6 hidden readers when dual-writes are dropped (atmosphere has a larger external read surface: ~105 sites across 9 files vs crop-uptake's ~47 across 7).

## Phasing

**Phase 0 (A-0.1)** verifies `swsublim` typed-config coverage and `spev`/`saev` adapter seeding. If gaps exist, close them with a config-field addition and an adapter seed line. If already covered, one documentation-only commit records the audit. This phase is expected to be 0–1 commit — lighter than boundary's 7-field Phase 0.

**Phase 1 (A-1.1 through A-1.9)** creates the `atmosphere_state_t` with two cohort sub-records, wires `atmosphere_init` at `swap.f90:188`, promotes `snow.f90`'s state arg to `intent(inout)`, and installs dual-writes in all 7 home-tree files. The large files (`meteoday.f90`, `meteodt.f90`) are now standard dual-write tasks, not "excluded file patches." Discovery line anchors (meteoday:706–865 for peva/ptra/atmdem; meteodt:348–349, 444–445 for peva/ptra) are the authoritative write-site inventory. The cohort reset blocks in `meteoday.f90:412–423` and `soilhydraulics.f90:1112–1142` are also dual-patched in this phase (legacy inline resets remain until Phase 2 drops them).

**Phase 2 (A-2.1 through A-2.7)** migrates the ~105 external read sites across 9 files (rootextraction, cropgrowth, tillage, boundtop, soilhydraulics, waterbalance, temperature, swapoutput, swap_csv_output), collapses duplicate reset blocks (D5, D6), drops dual-writes, retires ~33 globals from `variables.f90` with `[SS-ATM]` provenance markers, publishes ADR 0037, and updates the playbook. The double-reset removal for `igrai`/`inrai` (D6) happens in A-2.2 (waterbalance cutover) and A-1.8 (meteoday dual-write), coordinated so the single cohort-owned reset is live before the duplicate sites are deleted.

## Testing

**pFUnit:** `tests/unit/state/test_atmosphere_state.pf` (new file). Tests: scalar defaults all zero after `atmosphere_state_t()` literal; `intr%reset()` zeroes all 8 intermediate fields; `cumu%reset()` zeroes all 10 cumulative fields; two independent `atmosphere_state_t` instances share no state (no aliasing). All tests written before A-1.1 to establish failing-to-passing progression.

**check-full:** 5/5 byte-identical at every commit throughout Phase 1 and Phase 2. Dual-writes during Phase 1 are the safety net. Verification is mandatory before every commit; `pixi run check-full` failure blocks the task (per playbook lifecycle invariants and the `feedback_verify_before_committing.md` memory note — pFUnit alone misses global-default regressions).

## Non-goals

- Physics changes to `PenMon`, `reduceva`, `VonHHBraden`, `Gash`, `PartitionPrecipitation`, or `snow`. ASSOCIATE `at_` prefix preserves variable names; compute bodies change only in dual-write mechanics.
- Refactoring `meteoday.f90` or `meteodt.f90` orchestration structure — that is the future atmosphere REFACTOR arc.
- Adding regression coverage for `swetsine=1`, `swinter=3` (Rutter), `swsnow=1`, or `swredu=2` paths if currently uncovered.
- Migrating `nird`/`gird` — irrigation arc territory.
- Changing `soilwater_state_t` or `soilwater_init` — stable per ADR 0036 D4.

## ADR 0037

ADR 0037 will record: atmosphere as migration #7 (third coupling-surface arc); `atmosphere_state_t` with two cohort sub-records (`atmosphere_intermediate_t` 8 fields, `atmosphere_cumulative_t` 10 fields) plus 22 flat scalars; cohort pattern debuted at type-creation (not retrofitted); meteoday/meteodt included as home routines after scope correction; `atmosphere_init(state%atmosphere)` (scalars only, no dims) wired at `swap.f90:188`; `snow.f90` signature promotion; `igrai`/`inrai` double-reset collapsed to single cohort owner; ~33 globals retired; `pond` read and `nird`/`gird` ownership deferred. Cross-references discovery, design, plan, ADRs 0030–0036.
