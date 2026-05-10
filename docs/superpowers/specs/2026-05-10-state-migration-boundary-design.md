# Boundary-Conditions State-Type Migration — Design

**Date:** 2026-05-10
**Status:** accepted (pending implementation)
**Migration #:** 5 of N — first of the four coupling-surface arcs
**Branch:** `refactor/boundary-state`
**Discovery:** `docs/superpowers/specs/2026-05-10-state-migration-boundary-discovery.md`
**Predecessor ADRs:** 0030 (surfacewater), 0031 (drainage), 0032 (solute), 0033 (cumulative cohorts), 0034 (heat — flat-layout precedent)
**Successor ADR:** 0035

## Goal

Carve the 12 instantaneous-scalar boundary-condition fields out of `variables.f90` into a new `soilwater_state_t` (boundary subset) aggregated under `swap_state_t`. Close 7–8 Phase 0 config gaps in `bottom_boundary_config_t`. Both `boundtop` and `BoundBottom` already take `state` (from SS-HEAT Task 9) — the arc is field-carve + co-writer dual-writes + reader migration, not signature surgery. Every commit must be byte-identical on check-full.

## Out of scope

- `pond`, `gwl`, and `kmean(1)` / `kmean(numnod+1)` — deferred to the soil-water-core arc (Decisions D5–D7 below).
- The `cQMpLatSs` init reset at `soilhydraulics.f90:888` — macropore arc territory.
- Refactoring the duplicate `qbot` logic at the 12 `soilhydraulics.f90` co-write sites — the duplication is a coupling-surface artifact; leave in place.
- Cumulative bottom/top fluxes (`cqbot`, `cqbotdo`, `cqtdo`, etc.) — soil-water-core's cumulative cohort.
- The `hconduc` 14-caller `tsoil_node` sentinel residual from SS-HEAT — soil-water-core arc.

## Context

`boundtop` and `BoundBottom` are the coupling surfaces of the Richards-equation soil column with the outside world. Both are physics-thin (494 LoC total) — their content is small flux assignments per boundary-type switch. The arc work is concentrated in field-carving and co-writer plumbing, not algorithmic refactor.

Both entry points already take `state intent(in)` from SS-HEAT Task 9. The discovery cataloged 13 globals owned by boundary; one (`pond`) is deferred due to the tillage co-write hazard (discovery Section 4 / Hazard #4), leaving 12 fields for this arc. All 12 are instantaneous scalars (discovery Section 8). The flat `heat_state_t` precedent (ADR 0034) applies exactly. External readers span 14 files and ~80 read sites; co-writers number 5 files, of which `tillage.f90` and `calcgwl()` are deferred with `pond`/`gwl`.

## Decisions

### D1. Pilot scope: boundary home tree

Home tree: `src/boundary/boundtop.f90` + `boundbottom.f90`. Phase 0 config additions to `bottom_boundary_config_t` (`src/config/bottom_boundary_config.f90`). Co-writer dual-write targets: `soilhydraulics.f90`, `waterbalance.f90`, `frozencond.f90`. Reader-cutover targets: `solute.f90`, `surfacewater.f90`, `drainage.f90`, `swapoutput.f90`, `swap_csv_output.f90`.

### D2. State-type shape — flat

`soilwater_state_t` in `src/state/soilwater_state.f90`. All 12 owned fields are instantaneous scalars. No cohort sub-records in this arc. Layout (discovery Section 8):

```fortran
type :: soilwater_state_t
   ! Top-boundary (boundtop / PONDRUNOFF)
   real(real64) :: qtop           = 0.0_real64  !! top-surface flux (cm/d)
   real(real64) :: reva           = 0.0_real64  !! actual soil evaporation (cm/d)
   real(real64) :: hsurf          = 0.0_real64  !! pressure head at surface (cm)
   real(real64) :: runots         = 0.0_real64  !! runoff this step (cm)
   real(real64) :: QMpLatSs       = 0.0_real64  !! lateral macropore inflow (cm/d)
   logical      :: ftoph          = .false.     !! pressure-head top-boundary flag
   logical      :: FlRunoff       = .false.     !! runoff-potential flag
   ! Bottom-boundary (BoundBottom)
   real(real64) :: qbot           = 0.0_real64  !! bottom flux (cm/d)
   real(real64) :: qbot_nonfrozen = 0.0_real64  !! bottom flux pre-frost snapshot (cm/d)
   real(real64) :: hbot           = 0.0_real64  !! prescribed head at bottom (cm)
   real(real64) :: gwlinp         = 0.0_real64  !! prescribed gwl, swbotb=1 (cm)
   real(real64) :: deepgw         = 0.0_real64  !! deep-aquifer head, swbotb=3 (cm)
   ! Deferred (soil-water-core arc): pond, gwl, kmean(:), theta(:), h(:), cohort sub-records.
end type soilwater_state_t
```

### D3. Aggregation

`swap_state.f90` adds `type(soilwater_state_t) :: soilwater` alongside the existing `heat`, `drainage`, `solute`, `surfacewater` fields.

### D4. Argument-threading — no new plumbing at entry points

`boundtop(state)` and `BoundBottom(state)` already take `state intent(inout)` (SS-HEAT Task 9). No new plumbing needed for the home tree. Use ASSOCIATE inside compute bodies with `bw_` prefix (boundary water) — follows the `ht_` heat precedent and the subsystem-prefix convention.

### D5. `pond` deferred to soil-water-core arc

`tillage.f90:301,327` co-writes `pond` via bare `use variables`. Migrating `pond` here forces tillage state-plumbing and an init-order guard. Defer alongside whole-array soil-water fields (discovery Hazard #4).

### D6. `gwl` deferred to soil-water-core arc

`calcgwl()` in `waterbalance.f90` co-writes `gwl` with no `state` arg. Same defer logic as pond (discovery Hazard #5).

### D7. `kmean(1)` / `kmean(numnod+1)` deferred to soil-water-core arc

The entire `kmean(:)` array belongs to the soil-water-core arc; carving only the two surface entries as scalars would produce an awkward split (discovery Section 2b).

### D8. `soilwater_init(state)` placement

Insert immediately after `CalcGrid()` at `swap.f90:183`, before `DoTillage(1)`. Allocates the 12 scalar fields (trivial; no `numnod`-dependent arrays in this arc). Placement after `CalcGrid` ensures forward compatibility when later arcs add per-node arrays sized by `numnod`.

### D9. `qbot_nonfrozen` retained in state

One-write (`boundbottom.f90:179`) one-read (`frozencond.f90:216`) snapshot. Kept in state for symmetry with `qbot`. Cheap; revisit if `FrozenBounds` is refactored.

### D10. `reva` is boundary-owned

`reva` is the **actual** soil evaporation — the result of the top-boundary hydraulic-cap decision (`boundtop.f90:126,128`, sole writers). The atmosphere arc owns `peva` (potential) and `empreva` (empirical reduced). Boundary arc owns `reva` (discovery Section 10, Q9).

### D11. Mini-sim writeback preserved

`swapoutput.f90:3745–3820` snapshots `qbot`/`gwl`/`pond` and writes back after a perturbation mini-sim. Same pattern as the drainage mini-sim (already resolved in SS-DRST). Retarget writeback: `qbot` moves to `state%soilwater%qbot`; `gwl` and `pond` stay as legacy globals until the soil-water-core arc lands (deferred per D5–D6).

### D12. `swbotb=-2` runtime mutation kept as legacy config mutation

`boundbottom.f90:96` permanently sets `swbotb = -2` under oven-dry conditions (a runtime mutation of a config variable). Do not relocate in this arc. Document the deviation in ADR 0035; revisit when config-vs-state separation is hardened (discovery Hazard #6).

### D13. Phase 0 closes 7–8 scalar config gaps

`bottom_boundary_config_t` is missing `sinmax`, `sinamp`, `sinave` (swbotb=2 sine), `cofqha`, `cofqhb`, `cofqhc`, `swcofqhc` (swbotb=4 exp), and `hplate` (swbotb=8 lysimeter). Task B-0.4 verifies the adapter writes all 8 globals. Task B-0.5 audits TOML fixtures for swbotb=2/4/8 coverage; if uncovered, documents as known-issue without adding a new fixture this arc (discovery Section 6).

## Phasing

**Phase 0 (B-0.1–B-0.5) — Config gap closure.** Add 8 scalar parameters to `bottom_boundary_config_t`, extend the TOML reader, update the `config_to_variables.f90` adapter, add pFUnit tests. Verify regression-case coverage for the three newly-closed boundary modes; document any coverage gaps as known-issues. Ships standalone; must not change check-full output.

**Phase 1 (B-1.1–B-1.4) — State type, init, dual-write.** Create `soilwater_state_t`, add to `swap_state_t`, wire `soilwater_init(state)` at `swap.f90:183`. Dual-write all 12 owned fields in the home tree (`boundtop`, `BoundBottom`, `PONDRUNOFF`) and at co-writer sites (`soilhydraulics.f90` 12-site `qbot` + 1-site `qtop`, `frozencond.f90:FrozenBounds` 3-site `qbot`). ASSOCIATE with `bw_` prefix. check-full must remain byte-identical throughout.

**Phase 2 (B-2.1–B-2.7) — Reader cutover and global retirement.** Migrate all 14 external reader files from bare globals to `state%soilwater%X`. Drop dual-writes in home tree and co-writers. Comment out the 12 retired globals in `variables.f90` with `! [SS-BND] retired 2026-05-10 — moved to state%soilwater%X (ADR 0035)` provenance markers.

## Testing

**pFUnit:** new `tests/unit/state/test_soilwater_state.pf` — default-value, scalar lifecycle, independent-instance tests. New config tests per Phase 0 task (validators for each promoted field). All existing tests must remain green at every commit.

**check-full:** 5/5 byte-identical at every commit. The swbotb=2/4/8 paths are likely inactive in regression cases; check-full proves the migration does not break the active paths. Any newly uncovered modes documented as known-issues per D13.

## Non-goals

- Rewriting boundary physics. ASSOCIATE preserves variable names; bodies of `boundtop`, `BoundBottom`, `PONDRUNOFF` change only in dual-write mechanics.
- Adding regression coverage for swbotb=2 sine, swbotb=4 exp, or swbotb=8 in this arc.
- Fixing the duplicate `qbot` compute logic at the 12 `soilhydraulics.f90` sites — duplication is a coupling-surface artifact, not migration-introduced.
- Migrating `pond`, `gwl`, or `kmean` entries in this arc.
- Plumbing `tillage.f90` or `calcgwl()` with `state` in this arc.

## ADR 0035

Records: boundary as migration #5; flat `soilwater_state_t` (12 scalar fields, no cohorts — same class as ADR 0034 heat); `pond`/`gwl`/`kmean` deferral rationale; `qbot_nonfrozen` retention for symmetry; `reva` ownership classification; mini-sim writeback retarget; `swbotb=-2` runtime mutation kept as legacy config mutation with documented deviation. Cross-references discovery, design, plan, ADRs 0030–0034.
