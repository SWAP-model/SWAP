# Soil-Water Core State-Type Migration — Design

**Date:** 2026-05-11
**Status:** accepted (pending implementation)
**Migration #:** 8 of N — FINAL of the four coupling-surface arcs (boundary → crop-uptake → atmosphere → **soil-water core**)
**Branch:** `refactor/soilwater-core-state`
**Discovery:** `docs/superpowers/specs/2026-05-11-state-migration-soilwater-core-discovery.md`
**Predecessor:** `docs/superpowers/specs/2026-05-11-state-migration-atmosphere-design.md` / ADR 0037
**Successor ADR:** 0038

## Goal

Carve the ~74 residual soil-water-owned globals (Richards-equation interior: `h`, `theta`, `q`, `kmean`, `dimoca`, `cofgen`, `FrArMtrx`, `pond`, `gwl`, `volact`, `wbalance`, and the full intermediate + cumulative cohort families) out of `variables.f90` into the existing `soilwater_state_t` in `src/state/soilwater_state.f90`. This is the FINAL of the four coupling-surface arcs; after it lands, `state%soilwater` is the complete soil-water state record and all deferred residuals from prior arcs are resolved. This arc also resolves the three open boundary deferrals (D5 `pond`, D6 `gwl`, D7 `kmean` from ADR 0035), the heat ADR 0034 residual (`hconduc` `tsoil_node` 0.0 sentinel), and retargets the swapoutput mini-sim writeback fully to state. Every commit must be byte-identical on check-full.

## Out of scope

- **Macropore arc territory** — `cQMpLatSs`, `qimmob`, `QExcMpMtx`, `ArMpSs`/`ArMpTp`, `dFdhMp`, `IcTopMp`, and FrArMtrx co-writes from the macropore FrArMtrx-driven path. `soilhydraulics.f90:889` `cQMpLatSs = 0.0d0` init-zero left as-is.
- **Future atmosphere REFACTOR arc** — structural reorganisation of `meteoday.f90`/`meteodt.f90` orchestration. State-migration only here; atmosphere cohort writes already landed in ADR 0037.
- **Irrigation arc** — `nird`/`gird` global ownership, `qssdi`/`qssdisum` write-side. The soil-water `intr` cohort carries `igird`/`inird`/`iprec` and the `cumu` cohort carries `cgird`/`cnird` (accumulated and reset from this arc's side); irrigation arc continues subset-reset of `cgird`/`cnird` until it lands.
- **Grid-dimensions arc** — `numnod`, `numlay`, `dz`, `z`, `disnod`, `layer`, `botcom`, `nod1lay`, `ztopcp`, `zbotcp`, `inpola`, `inpolb`, `numtab`, `sptab` stay as legacy globals per heat ADR 0034 precedent. A future `grid_t` arc is the destination.

## Context

The three preceding coupling-surface arcs have already done most of the heavy structural lifting. Boundary (ADR 0035) established `soilwater_state_t`, wired `soilwater_init`, and plumbed `state` into `boundtop`/`boundbottom`. Crop-uptake (ADR 0036) bumped `soilwater_init` to `(sw, numnod, nlay)` and added 22 per-node arrays. Atmosphere (ADR 0037) plumbed `state` into `DoTillage`, `checkmassbal`, `CNmethod`, `Consolidate_Bdens`, and the entire atmosphere subsystem. As a result, almost every entry point this arc needs to write through (`soilwater`, `headcalc`, `fluxes`, `integral`, `ConvertDiscrVert`, `checkmassbal`, `boundtop`, `boundbottom`, `tillage`) already has `state` in scope. The residual plumbing additions are mechanical: `SoilWaterStateVar(task, state)`, `hysteresis(state)`, `watstor(state)`, `calcgwl` intent-flip, plus threading into tillage internals.

The windfall payoff is substantial: the 5-entry-point residual required by this arc (SoilWaterStateVar, hysteresis, watstor, level/watertable, calcgwl) is the smallest plumbing overhead of any arc, despite owning 74 fields — the largest field count by 30%. The cohort consolidation is also the most impactful: `soilhydraulics.f90:1083–1158` (the entire `flzerointr`+`flzerocumu` reset block, 75+ lines) collapses to three type-bound procedure calls.

This arc introduces a **mild ADR 0033 extension** for the per-day cohort: the 8 `flDayStart`-gated per-day fields fold into `soilwater_intermediate_t` with a SECOND type-bound method `reset_per_day()`, distinct from the existing `reset()` / `flzerointr` gate. This keeps `soilwater_state_t`'s shape parallel to `atmosphere_state_t` (two cohort sub-records, not three) while expressing both reset cadences in the type system.

## Decisions

- **D1. Pilot scope: soil-water core home tree.** Files: `soilhydraulics.f90`, `waterbalance.f90`, `soilgrid.f90`, `sptabulated.f90` (third-party; no SWAP writes), `WC_K_models_04_11.f90`, `soilhydraulicsutils.f90`. The `soilgrid.f90:ConvertDiscrVert` body reads ~10 soil-water globals via `use Variables`; retarget to `state%soilwater%intr%X` in Phase 2. Discovery Section 1 confirms all high-traffic home-tree routines already have `state`.

- **D2. State-type extension — TWO cohorts mirroring atmosphere ADR 0037.** Extend existing `soilwater_state_t` with 30 flat instantaneous + 8 per-day + 2 cohort sub-records:

  ```fortran
  type :: soilwater_intermediate_t
     ! per-node intra-period flux accumulators (3 arrays)
     real(real64), allocatable :: inq(:), inqrot(:), inqssdi(:)
     real(real64), allocatable :: iqdo(:), iqup(:), IThetaBeg(:)
     ! scalars (16)
     real(real64) :: iqrot, iqssdi
     real(real64) :: iqredwet, iqreddry, iqredsol, iqredfrs
     real(real64) :: ies0, iet0, iew0, iintc
     real(real64) :: iruno, irunoCN, irunon
     real(real64) :: iqbot, iqtdo, iqtup
     real(real64) :: IPondBeg
     real(real64) :: iprec, igird, inird
     ! per-day cohort (flDayStart gate — reset_per_day() only)
     real(real64) :: tra
     real(real64) :: iqredwet_day, iqreddry_day, iqredsol_day, iqredfrs_day, iptra_day
     real(real64), allocatable :: qpotrot_day(:), qredtot_day(:)
  contains
     procedure :: reset         => soilwater_intermediate_reset      ! flzerointr
     procedure :: reset_per_day => soilwater_per_day_reset           ! flDayStart
  end type soilwater_intermediate_t

  type :: soilwater_cumulative_t
     real(real64) :: cqssdi, cqrot, cqbot, cqbotdo, cqbotup
     real(real64) :: cinund, crunon, crunoff, crunoffCN
     real(real64) :: cqtdo, cqtup, cqprai, cgird, cnird
  contains
     procedure :: reset => soilwater_cumulative_reset
  end type soilwater_cumulative_t
  ```

  Per-day fields fold into `soilwater_intermediate_t` (not a separate `per_day_t`) per discovery Section 8 recommendation. This is the mild ADR 0033 extension: same cohort type, two distinct reset procedures, two distinct activity gates (`flzerointr` vs `flDayStart`). Field count: intermediate 22 total (14 scalars + 3 arrays + 5 per-day scalars + 2 per-day arrays); cumulative 14 scalars. The `reset()` method zeroes all 22 fields including per-day; `reset_per_day()` zeroes only the 8 per-day fields.

- **D3. Aggregation — no `swap_state.f90` change.** `soilwater_state_t` is the record already aggregated in `swap_state_t` as `state%soilwater`. This arc extends it in place; no new component, no new `use` statement in `swap_state.f90`.

- **D4. `soilwater_init` signature unchanged.** Signature is already `soilwater_init(sw, numnod, nlay)` from crop-uptake C-1.2 (ADR 0036). The body is extended to allocate the 24 new per-node/per-layer allocatables (`theta`, `thetm1`, `thetar`, `thetas`, `h`, `hm1`, `q`, `k`, `kmean`, `dimoca`, `cofgen(21,numnod)`, `FrArMtrx`, `fluseksatexm`, `indeks`, `evp`, `thetsl(numlay)`, plus cohort arrays `inq`, `inqrot`, `inqssdi`, `iqdo`, `iqup`, `IThetaBeg`, `qpotrot_day`, `qredtot_day`) and initialise scalar fields.

- **D5. Cohort reset consolidation — atmosphere ADR 0037 D5 pattern.** Replace `soilhydraulics.f90:1083–1158` (~75 scattered zero-assignments) with:
  ```fortran
  if (flzerointr) call state%soilwater%intr%reset()
  if (flDayStart) call state%soilwater%intr%reset_per_day()
  if (flzerocumu) call state%soilwater%cumu%reset()
  ```
  The `pondini = pond` and `volini = volact` rebases that immediately follow the `flzerocumu` block remain as inline non-zero rebases (same pattern as solute's `samini = sampro`, ADR 0033 Phase B). ~55 scattered reset lines collapse to 3 call lines.

- **D6. Resolve boundary D5 (`pond` deferral).** `pond` migrates this arc to `state%soilwater%pond`. `DoTillage(iTask, state)` is already state-plumbed (A-2.6 windfall). Internal subs `Change_MvGpars`, `Adapt_WC_H` do not take `state`; thread state through them. `boundtop.f90` already has state; retarget `pond` writes at 8 sites. `swapoutput.f90` mini-sim snapshot+restore retargeted. `et.f90:reduceva` pond-read marker `[SS-ATM] reads legacy pond` becomes a real cutover. `surfacewaterutils.f90` and other readers get standard reader cutover.

- **D7. Resolve boundary D6 (`gwl` deferral).** `gwl` migrates to `state%soilwater%gwl`. `calcgwl(state)` already has `intent(in)`; flip to `intent(inout)` and retarget ~12 write sites (`gwl`, `nodgwl`, `pegwl`, `bpegwl`, `npegwl`, `gwlflcpzo`, `nodgwlflcpzo`) to `state%soilwater%X`. External readers (drainage, surfacewater, divdra, agetracer, et, frozencond, output files) get standard reader cutover.

- **D8. Resolve boundary D7 (`kmean` full-array deferral).** `kmean(:)` migrates to `state%soilwater%kmean(:)`. `boundtop` writes `kmean(1)` at 3 sites; `boundbottom` writes `kmean(numnod+1)` at 2 sites; `soilhydraulics` writes interior elements at ~14 sites. All callers are state-aware post-arc-7.

- **D9. Resolve heat ADR 0034 `hconduc` `tsoil_node` sentinel.** All 14+ callers of `hconduc(node, h, th, rfcp, [tsoil_node])` that currently omit `tsoil_node` or pass a literal (`0.0d0`, `10.d0`, `1.0d0`) are updated to pass `state%heat%tsoil(node)`. Sites: `soilhydraulics.f90` ~10 sites, `boundbottom.f90:171`, `boundtop.f90:111`, `rootextraction.f90:868,873`, `tillage.f90:361`, `swapoutput.f90:2103`, `macropore.f90:1142`. All callers have `state` in scope post-arc-7. Single coordinated task.

- **D10. `swapoutput` mini-sim writeback retarget.** Boundary B-2.6 retargeted `qbot`; this arc retargets the remaining 4 fields: `gwl`, `pond`, `theta(:)`, `h(:)` snapshot+restore at `swapoutput.f90:3865–3948`. Identical mechanical pattern to B-2.6. Also: `SoilWaterStateVar(task)` gains a `state` arg (3 call sites: `soilhydraulics.f90`, `swapoutput.f90:3903`, `swapoutput.f90:3915`); `kmean(numnod+1) = k(numnod)` restore write retargeted to state.

- **D11. `cQMpLatSs` init reset DEFERRED.** `soilhydraulics.f90:889` `cQMpLatSs = 0.0d0` stays as-is. Macropore arc owns this field (discovery Section 2.8). No action this arc.

- **D12. `swbotb=-2` runtime mutation KEPT AS LEGACY.** Boundary D12 carry-forward. `boundbottom.f90:97–99` toggle unchanged. Documented deviation; no action this arc.

- **D13. Grid dimensions stay legacy globals.** `numnod`, `numlay`, `dz`, `z`, `disnod`, `layer`, `botcom`, `nod1lay`, `ztopcp`, `zbotcp`, `inpola`, `inpolb`, `numtab`, `sptab`. Matches heat ADR 0034 grid-dims decision and boundary D7/D8 carry-forward pattern.

- **D14. Multi-owner cumulative pattern: `cgird`/`cnird`.** Kept in `soilwater_cumulative_t`. `irrigation.f90:93–95` continues to call a subset-reset of these two fields until the irrigation arc takes ownership. Same documented multi-owner pattern as atmosphere's `nird`/`gird` deferral (ADR 0037 D10). No structural change in `irrigation.f90` this arc.

- **D15. Compile-driven Phase 2.12 expectation — 20–25 hidden readers.** Largest compile-driven phase of any arc to date (boundary had 4, crop-uptake 2, atmosphere 13). `theta` has ~14 external reader files; `h` ~12; `gwl`/`pond` ~9 each. Budget 10–15 iterative fixup commits. Apply the atmosphere A-2.1 pre-flight dual-write coverage check before every reader-file cutover.

## Phasing

**Phase 1 (S-1.1 through S-1.12)** creates the two cohort types, extends `soilwater_state_t` with 30 flat instantaneous + 8 per-day fields, extends `soilwater_init` to allocate all new arrays, and installs dual-writes in all home-tree and co-writer files. The SoilWaterStateVar/hysteresis/watstor/calcgwl plumbing residual is also resolved in Phase 1. Per the atmosphere arc precedent, cohort reset blocks in soilhydraulics (1083–1158) are dual-patched here; the cohort `reset()` replacement and legacy-inline-zero removal happen in Phase 2.

**Phase 2 (S-2.1 through S-2.13)** migrates all external readers (23 files, ~300–400 sites), performs cohort reset consolidation (S-2.1), resolves the hconduc sentinel (S-2.2), retargets the mini-sim writeback (inside S-2.11), drops dual-writes, retires ~74 globals from `variables.f90`, and publishes ADR 0038 + playbook update.

## Testing

**pFUnit:** `tests/unit/state/test_soilwater_state.pf` — extended for new cohort types. Tests: `soilwater_intermediate_t()` literal default to zero; `intr%reset()` zeroes all 22 fields; `intr%reset_per_day()` zeroes only the 8 per-day fields (seed non-zero values on full intr, call reset_per_day, assert per-day fields are zero and non-per-day fields are unchanged); `cumu%reset()` zeroes all 14 fields; allocatable array reset guards (`if (allocated)` — arrays are allocated before reset is called in tests); two independent `soilwater_state_t` instances share no state. Minimum ~10 new tests written before S-1.1 to establish failing-to-passing progression.

**check-full:** 5/5 byte-identical at every commit throughout Phase 1 and Phase 2. Dual-writes during Phase 1 are the safety net. `pixi run check-full` is mandatory before every commit; pFUnit alone misses global-default regressions (per `feedback_verify_before_committing.md`).

## Non-goals

- Physics changes to `headcalc`, `hysteresis`, `calcgwl`, `watstor`, `hconduc`. ASSOCIATE blocks preserve variable names; compute bodies change only in dual-write mechanics and plumbing.
- Refactoring the Richards solver algorithm, the Mualem-VG parameter structure, or the soil grid machinery.
- Adding regression coverage for `iHWCKmodel=4–11` temperature-dependent K paths, macropore paths, or the `swbotb=2/4/8` modes (still inactive in TOML regression cases; documented in ADR 0035 Phase 0).
- Migrating `qssdi`/`qssdisum` (SSDI/irrigation write-side) — they remain as legacy globals; only the soil-water-side accumulators `cqssdi`/`iqssdi`/`inqssdi` migrate here.
- Migrating `cgird`/`cnird` ownership to `irrigation.f90` — irrigation arc territory.
- Moving `evp(:)` (always-zero per-node field) — migrated semantics-only; zero overhead; preserved for legacy parity.

## ADR 0038 stub

ADR 0038 will record: soil-water core as migration #8 — FINAL coupling-surface arc; `soilwater_state_t` extended with 30 flat-instantaneous + 8 per-day + `soilwater_intermediate_t` (22 fields, two methods: `reset()` / `reset_per_day()`) + `soilwater_cumulative_t` (14 fields, `reset()`); mild ADR 0033 extension via second type-bound method on `intr`; boundary D5 (`pond`) / D6 (`gwl`) / D7 (`kmean`) deferrals resolved; heat ADR 0034 `hconduc` `tsoil_node` sentinel eliminated (14 sites); `swapoutput` mini-sim writeback fully retargeted (`gwl`/`pond`/`theta`/`h`); grid dimensions retained as legacy globals (heat ADR 0034 precedent); `cgird`/`cnird` kept in soil-water cumu cohort with multi-owner comment; `cQMpLatSs` deferred to macropore arc; ~74 globals retired; `soilwater_init` signature stable. Cross-references discovery, design, plan, ADRs 0030–0037.
