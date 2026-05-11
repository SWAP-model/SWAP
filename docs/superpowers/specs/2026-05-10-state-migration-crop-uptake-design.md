# Crop Water Uptake State-Type Migration — Design

**Date:** 2026-05-11
**Status:** accepted (pending implementation)
**Migration #:** 6 of N — second of the four coupling-surface arcs (boundary → **crop-uptake** → atmosphere → soil-water core)
**Branch:** `refactor/crop-uptake-state`
**Discovery:** `docs/superpowers/specs/2026-05-10-state-migration-crop-uptake-discovery.md`
**Predecessor:** `docs/superpowers/specs/2026-05-10-state-migration-boundary-design.md` / ADR 0035
**Successor ADR:** 0036

## Goal

Carve the ~22 crop-water-uptake globals (`qrot`, `qrosum`, four per-node stress arrays, four stress-sum scalars, the JvL microscopic-uptake state `Tactual`/`alpJvLier`/`hleaf`/`Hxylem`, six per-node JvL arrays, `mfluxtable`, and `flWrtNonox`) out of `variables.f90` into a flat extension of the existing `soilwater_state_t`. Resolve the `qrot` ownership ambiguity (flagged in the soil-water mega-discovery Section 2/7) by giving `rootextraction.f90` undisputed ownership of `qrot` as soil-water state. Both `RootExtraction` and `JongvanLier`/`JongvanLierLoop`/`MatricFlux`/`OxygenStress` already take `state` (SS-HEAT + SS-SLST windfall) — the arc is field-carve + co-writer dual-writes + reader cutover, not signature surgery. Every commit must be byte-identical on check-full.

## Out of scope

- **Cumulatives:** `cqrot`, `iqrot`, `inqrot`, `iqredwet/dry/sol/frs`, `iqredwet_day/dry_day/sol_day/frs_day`, `qpotrot_day`, `qredtot_day`, `iptra_day` — owned by `waterbalance.f90` / `soilhydraulics.f90` and assigned to the soil-water-core arc (#8).
- **`pond` / `gwl` / `kmean`** — still deferred per boundary D5/D6/D7; untouched here.
- **`hconduc` `tsoil_node` sentinel residual** from SS-HEAT — soil-water-core arc territory.
- **3 cropgrowth Task-1 init paths refactored** — only state-plumbed here (write to `state%soilwater`); no structural reorganization of the three duplicated JvL init blocks (out of scope per discovery Q8).
- **`swap_state.f90` structural change** — `soilwater_state_t` extends in place; no new sub-record is added to `swap_state_t`.

## Context

`RootExtraction(state)` runs once per outer Richards timestep (`swap.f90:302`), before `SoilWater(2)` / `headcalc`. Its outputs `qrot(:)` and `qrosum` are the primary per-node sink term consumed by Richards (10 read sites in `soilhydraulics.f90`) and the mass-balance anchor in `waterbalance.f90`. The JvL microscopic-uptake machinery (`JongvanLier`, `JongvanLierLoop`, `MatricFlux`) is self-contained: write-sites are all inside `rootextraction.f90` except for three `cropgrowth.f90` Task-1 seeds (`hroot(:) = h(:)`, `hleaf = -2000.d0`, and the `MatricFlux(1)` lookup build).

The boundary arc established `soilwater_state_t` with 12 scalar fields. This arc is the anticipated next step: discovery Section 8 and boundary ADR 0035 both explicitly named crop-uptake as the arc that adds per-node arrays. Boundary lesson #8 placed `soilwater_init` at `swap.f90:183` (after `CalcGrid()`, before `DoTillage(1)`) precisely to allow this arc's array allocations to land without an init-order gap.

## Decisions

### D1. Pilot scope

Home file: `src/crop/rootextraction.f90`. Co-writer paths: `cropgrowth.f90` Task-1 init (CropFixed lines 595–612, CropWofost lines 1325–1342, CropGrass lines 2394–2411). External readers: `soilhydraulics.f90` (10 qrot sites), `waterbalance.f90` (~12 sites), `swapoutput.f90` (~14 sites), `solute.f90` (3 sites), `agetracer.f90` (3 sites), `cropgrowth.f90` (5 `flWrtNonox` compute sites). `swap_csv_output.f90` is indirect-only (reads downstream cumulatives, not this arc's instantaneous fields) — migrates only if compile-driven.

### D2. State-type shape — flat extension

Extend the existing `soilwater_state_t` in `src/state/soilwater_state.f90`. No new cohort sub-records. All 22 owned fields are instantaneous or init-once (`mfluxtable`). Matches ADR 0034 (heat) + ADR 0035 (boundary) flat-layout precedent. Target shape (appended after the 12 existing boundary fields):

```fortran
! ── Crop-uptake arc (SS-CRP) — primary Feddes + stress path ───────────────────
real(real64), allocatable :: qrot(:)        ! per-node root sink (cm/d)
real(real64), allocatable :: qpotrot(:)     ! per-node potential uptake (cm/d)
real(real64), allocatable :: qredwet(:)     ! per-node wet-stress contribution
real(real64), allocatable :: qreddry(:)     ! per-node dry-stress contribution
real(real64), allocatable :: qredsol(:)     ! per-node salt-stress contribution
real(real64), allocatable :: qredfrs(:)     ! per-node frost-stress contribution
real(real64) :: qrosum      = 0.0_real64   ! column-sum root uptake (cm/d)
real(real64) :: qredwetsum  = 0.0_real64   ! column-sum wet-stress reduction
real(real64) :: qreddrysum  = 0.0_real64   ! column-sum dry-stress reduction
real(real64) :: qredsolsum  = 0.0_real64   ! column-sum salt-stress reduction
real(real64) :: qredfrssum  = 0.0_real64   ! column-sum frost-stress reduction
logical      :: flWrtNonox  = .false.      ! per-step non-ox flag (written by RootExtraction)

! ── JvL microscopic uptake (swdrought=2) ─────────────────────────────────────
real(real64), allocatable :: mflux(:)       ! matric-flux potential per node (L2/T)
real(real64), allocatable :: mroot(:)       ! matric-flux at root surface per node
real(real64), allocatable :: hroot(:)       ! pressure head at root surface per node (cm)
real(real64), allocatable :: rootrho(:)     ! root radial geometry factor per node
real(real64), allocatable :: rootphi(:)     ! root azimuthal geometry factor per node
real(real64), allocatable :: rmax(:)        ! maximum root uptake flux per node
real(real64), allocatable :: mfluxtable(:,:)! (numlay, 801) — init-once matric-flux lookup
real(real64) :: Tactual    = 0.0_real64    ! prior-step actual transpiration (cm/d)
real(real64) :: alpJvLier  = 0.0_real64    ! JvL alpha parameter
real(real64) :: hleaf      = 0.0_real64    ! leaf water potential (cm)
real(real64) :: Hxylem     = 0.0_real64    ! xylem water potential (cm)
```

**Total new fields: 22** (6 primary per-node + 5 primary scalars + 1 flag + 6 JvL per-node + 4 JvL scalars + 1 layer×801 lookup).

### D3. Aggregation — no swap_state.f90 change

`soilwater_state_t` extends in place. `swap_state_t` already holds `type(soilwater_state_t) :: soilwater` (ADR 0035). No structural change to `swap_state.f90`.

### D4. `soilwater_init` signature change

New signature: `subroutine soilwater_init(sw, numnod, nlay)`. Allocates all 12 per-node arrays to `numnod` and `mfluxtable` to `(nlay, 801)`. Zeros all scalars. Caller at `swap.f90:184` updates to `call soilwater_init(state%soilwater, numnod, numlay)`. Both `numnod` and `numlay` are available immediately after `CalcGrid()` at line 183 (confirmed: `swap.f90:71` imports both from `variables`). ASSOCIATE prefix `cw_` (crop-water) used inside `rootextraction.f90` and `cropgrowth.f90` write paths.

### D5. `CropGrowth(1, state)` plumbing

Extend the `CropGrowth` Task-1 call at `swap.f90:272` from `CropGrowth(1, state%heat%tsoil)` to `CropGrowth(1, state)`. Update all three Task-1 init paths in `cropgrowth.f90` to write `state%soilwater%hroot(i)`, `state%soilwater%hleaf`, and call `MatricFlux` in a way that populates `state%soilwater%mfluxtable`. The heat-arc `tsoil(:)` dummy-arg threading confirmed that module-routine cropgrowth paths already accept slice or full-state threading without interface complications.

### D6. `mfluxtable` initialization relocated to `soilwater_init`

The `MatricFlux(1)` lookup-build (currently triggered by three `cropgrowth.f90` Task-1 calls at lines 601, 1331, 2400) is relocated into `soilwater_init`. This converts the only Cat-4 init-seed from outside the home tree into a self-contained state-lifecycle operation, eliminating a cropgrowth co-write. `soilwater_init` calls the private helper equivalent of `MatricFlux(1)` (or delegates to `MatricFlux(1, h(1), 1, dummy, state%soilwater)` if the signature permits; resolve during implementation). The `cropgrowth.f90` Task-1 call sites drop their `MatricFlux(1)` calls after relocation.

### D7. JvL fields `Tactual`, `alpJvLier`, `rmax` — flat scalars/array on state

Symmetric with `qrot`, `qrosum`. Zero external readers outside `rootextraction.f90` for all three. Trivial migration: dual-write in `rootextraction.f90`, then drop. `rmax(:)` is a per-node array (allocatable).

### D8. `flWrtNonox` ownership — `state%soilwater` for this arc

`flWrtNonox` is written by `RootExtraction` (sole writer) and read by `cropgrowth.f90` at 5 sites. Consistent with the owned-and-touched rule from boundary. Document in ADR 0036 that a future crop-state arc may relocate it; semantically it is more crop than soil-water.

### D9. Macropore co-write disposition — no action

Discovery confirmed no `qrot` co-write by macropore. Discovery Section 4 / Hazard #9: macropore reads `inqrot` (downstream cumulative, out of scope) but does not write any of the 22 owned fields.

### D10. Cumulatives out of scope

`cqrot`, `iqrot`, `inqrot`, `iqredwet/dry/sol/frs`, `qpotrot_day`, `qredtot_day` belong to the soil-water-core arc (#8, cumulative cohort under `flzerocumu` gate). The `waterbalance.f90:418–435` block that accumulates from this arc's instantaneous scalars continues to read from `state%soilwater%X` after Phase 2 cutover, but the writes (into the legacy cumulatives) stay until the core arc.

### D11. Phase 0 regression coverage audit

Verify that `check-full`'s five TOML cases exercise `swdrought=2` (JvL), `swoxygen=2` (Bartholomeus), `swcompensate≥1` (Jarvis/Walsum), and `swfrost=1`. The JvL path spans lines 313–760 of `rootextraction.f90` — if uncovered, this is a notable risk. Document any gap as a known-issue in ADR 0036; do NOT add new fixtures in this arc. (Discovery Section 6 Phase 0 candidates: estimated 0–1 commit; all typed-config fields already covered.)

### D12. Compile-driven Phase 2.7 expected — 2–4 hidden readers

Boundary lesson #5 applies. After dropping dual-writes, expect 2–4 hidden readers from `swapoutput.f90` use clauses, `waterbalance.f90` use clauses, and possibly `agetracer.f90` or `solute.f90` extra reads of `qpotrot`/`qredXXX`. Plan for 2–4 compiler-driven fix-up commits.

## Phasing

**Phase 0 (C-0.1) — Regression coverage audit.** Grep the five TOML regression cases for `swdrought`, `swoxygen`, `swcompensate`, `swfrost` switch values. Document results; if any path is uncovered, note as a known-issue. No typed-config fields needed (full coverage exists). Zero or one documentation-only commit.

**Phase 1 (C-1.1–C-1.4) — State-type extension, init signature bump, co-writer plumbing, dual-write.** Extend `soilwater_state_t` with 22 fields; bump `soilwater_init` to `(sw, numnod, nlay)` and allocate all per-node arrays; relocate `mfluxtable` build to `soilwater_init`; plumb `CropGrowth(1, state)` and update three Task-1 init paths; install dual-writes throughout `rootextraction.f90`. check-full must remain byte-identical after every commit.

**Phase 2 (C-2.1–C-2.6) — Reader cutover and global retirement.** Migrate reads in `soilhydraulics.f90`, `waterbalance.f90`, `oxygenstress.f90`/`solute.f90`/`cropgrowth.f90`, output routines; drop dual-writes; retire 22 legacy globals from `variables.f90`; publish ADR 0036 and update the playbook.

## Testing

**pFUnit:** extend `tests/unit/state/test_soilwater_state.pf` with new `@test` stubs for per-node array allocation, zero-default scalars, and independent-instance isolation after the `soilwater_init(sw, numnod, nlay)` call. All existing tests must remain green at every commit.

**check-full:** 5/5 byte-identical at every commit. Dual-writes during Phase 1 are the safety net; verification is mandatory before every commit. If the Phase 0 audit reveals that swdrought=2 is uncovered, the JvL path has no integration-level regression coverage — document this gap explicitly in ADR 0036.

## Non-goals

- Rewriting `RootExtraction`, `JongvanLier`, or `JongvanLierLoop` physics. ASSOCIATE `cw_` prefix preserves variable names; compute bodies change only in dual-write mechanics.
- Refactoring the three duplicated JvL init blocks in `cropgrowth.f90` into a shared helper (noted for a future cropgrowth refactor arc).
- Adding regression coverage for `swdrought=2`, `swoxygen=2`, or `swcompensate` paths.
- Migrating `OxygenStress`-owned fields — `OxygenStress` owns no globals in the 22-field owned set; it is an external callee that receives stress output (not a writer of this arc's fields).
- Changing `swap_state.f90` structure.

## ADR 0036

Records: crop-uptake as migration #6; flat extension of `soilwater_state_t` with 22 fields (6 primary per-node + 5 primary scalars + 1 flag + 6 JvL per-node + 4 JvL scalars + 1 layer×801 lookup); `soilwater_init` signature bump to `(sw, numnod, nlay)`; `CropGrowth(1, state)` plumbing rationale; `mfluxtable` init relocation to `soilwater_init`; JvL flat-scalar layout; `flWrtNonox` pragmatic ownership in `state%soilwater` with note for future relocation; cumulatives out of scope; `rmax`/`Tactual`/`alpJvLier` kept on state for JvL symmetry. Cross-references discovery, design, plan, ADRs 0030–0035.
