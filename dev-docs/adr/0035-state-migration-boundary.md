---
title: "ADR 0035 — Boundary subsystem state-type migration"
date: 2026-05-10
status: accepted
---

# ADR 0035: Boundary subsystem state-type migration

**Status:** accepted (Phases 0 + 1 + 2 complete)
**Date:** 2026-05-10
**Migration #:** 5 of N — first coupling-surface arc of the soil-water 4-arc decomposition
**Branch:** `refactor/boundary-state`

## Context

Migrations 1–4 (surfacewater ADR 0030, drainage ADR 0031, solute ADR 0032, heat ADR 0034) established the full playbook. ADR 0033 added cumulative-reset cohorts. Soil-water is migration #5–#8: decomposed by coupling surface rather than carved as one oversized arc. The four arcs are:

- **#5 Boundary** (this arc) — top/bottom coupling-surface fields: `qtop`, `qbot`, and 10 others.
- **#6 Crop-uptake** — crop-root sink fields.
- **#7 Atmosphere** — potential-evaporation, precipitation, and snow fields.
- **#8 Soil-water core** — Richards-equation interior: `h`, `theta`, `q`, `kmean`, `pond`, `gwl`, cumulative cohorts.

The decomposition-by-coupling-surface strategy was adopted because the soil-water subsystem touches ~100 legacy globals — too large for one arc. Each arc carves the fields belonging to its surface; the final core arc handles the Richards interior. This minimizes blast radius per commit and produces clean ADR narratives.

`boundtop` and `BoundBottom` already took `state intent(inout)` from SS-HEAT Task 9. The arc was therefore field-carve + co-writer dual-writes + reader migration, not signature surgery.

## Decision

**Migrate 12 boundary-owned globals** into a new flat `soilwater_state_t`, aggregated under `swap_state_t`. The boundary subset is the first population of `state%soilwater`; subsequent arcs expand it.

Three phases:

- **Phase 0** — Promote 8 missing scalar physics-config fields into `bottom_boundary_config_t`.
- **Phase 1** — Define `soilwater_state_t`, add to `swap_state_t`, wire `soilwater_init`, dual-write 12 owned fields in home tree and co-writer sites.
- **Phase 2** — Migrate ~165 reads across 8 files, drop dual-writes, retire 12 globals.

### Phase 0 — Physics config promotion

8 fields promoted into `bottom_boundary_config_t` and the TOML reader, closing silent-zero-default gaps for three boundary modes:

- **swbotb=2 sine:** `sinmax`, `sinamp`, `sinave`
- **swbotb=4 q(h) exponential:** `cofqha`, `cofqhb`, `cofqhc`, `swcofqhc`
- **swbotb=8 lysimeter:** `hplate`

Same pattern as heat's Phase 0 (`swcalt=1` analytical-method gaps). The swbotb=2/4/8 paths are inactive in the five TOML regression cases; documented as a known coverage gap (B-0.5).

### Phase 1 — State type and threading

`src/state/soilwater_state.f90` defines flat `soilwater_state_t` (12 scalar fields, no cohorts):

```fortran
type :: soilwater_state_t
   ! Top-boundary (boundtop / PONDRUNOFF)
   real(real64) :: qtop           = 0.0_real64
   real(real64) :: reva           = 0.0_real64
   real(real64) :: hsurf          = 0.0_real64
   real(real64) :: runots         = 0.0_real64
   real(real64) :: QMpLatSs       = 0.0_real64
   logical      :: ftoph          = .false.
   logical      :: FlRunoff       = .false.
   ! Bottom-boundary (BoundBottom)
   real(real64) :: qbot           = 0.0_real64
   real(real64) :: qbot_nonfrozen = 0.0_real64
   real(real64) :: hbot           = 0.0_real64
   real(real64) :: gwlinp         = 0.0_real64
   real(real64) :: deepgw         = 0.0_real64
end type soilwater_state_t
```

`swap_state_t` gains `type(soilwater_state_t) :: soilwater`. `soilwater_init(state)` placed at `swap.f90:183` — after `CalcGrid()`, before `DoTillage(1)` — ensuring forward compatibility when later arcs add per-node arrays.

Dual-write installed in both home tree (`boundtop.f90`, `boundbottom.f90`) and co-writer sites (`soilhydraulics.f90` 12-site `qbot` + 1-site `qtop`; `frozencond.f90:FrozenBounds` 3-site `qbot`).

### Phase 2 — Cross-subsystem migration and cleanup

**Reader files migrated (~165 read sites, 8 files):**

- `soilhydraulics.f90` (~25 sites): qtop, qbot, reva, hsurf, ftoph, FlRunoff, hbot, gwlinp, deepgw.
- `waterbalance.f90` (~10 sites): reva, runots, qbot (integral + checkmassbal).
- `solute.f90` + `agetracer.f90` (~7 sites): qtop, qbot.
- `surfacewater.f90` + `drainage.f90` (~15 sites): runots, qtop indirect.
- `frozencond.f90` (3 sites): qbot, qbot_nonfrozen.
- `swapoutput.f90` + `swap_csv_output.f90` (~15 sites): all owned fields in output columns.

**Mini-sim writeback retargeted (D11).** `swapoutput.f90:3745–3820` snapshots and restores `qbot`/`gwl`/`pond` for the stochastic mini-sim. After migration: `qbot` writeback targets `state%soilwater%qbot`. `gwl` and `pond` remain as legacy globals until the soil-water-core arc (D5–D6).

**Compile-driven Phase 2.7 (4 hidden readers surfaced).** After dropping dual-writes, the compiler identified 4 read sites that the read-only inventory missed:

1. `macropore.f90:MACROINTEGRAL` — `QMpLatSs` reader in macropore cumulative path.
2. `waterbalance.f90:calcgwl` — `qbot` read in a flux back-calc fallback branch.
3. `config_to_variables.f90` — vestigial `hbot` write (became a no-op after Phase 0; removed).
4. `soilhydraulics.f90` convergence-check block — `qtop` read in a diagnostics path not visible during grep-based discovery.

**~68 dual-write legacy writes dropped** at Phase 2.7; 12 globals commented out in `variables.f90` with `! [SS-BND] retired 2026-05-10 — moved to state%soilwater%X (ADR 0035)` provenance markers.

### Deferred (out of scope)

- **`pond`** — co-written by `tillage.f90:301,327` (no state plumbing); deferred to soil-water-core arc (D5).
- **`gwl`** — co-written by `calcgwl()` in `waterbalance.f90` (no state plumbing); deferred to soil-water-core arc (D6).
- **`kmean(1)` / `kmean(numnod+1)`** — carving only two entries of the full `kmean(:)` array would produce an awkward split; deferred with whole array to soil-water-core arc (D7).
- **`cQMpLatSs` init reset** at `soilhydraulics.f90:888` — macropore arc territory (Hazard #7).
- **`swbotb=-2` runtime mutation** — `boundbottom.f90:96` permanently sets `swbotb=-2` under oven-dry conditions. Kept as legacy config mutation; documented deviation (D12).
- **Cumulative boundary fluxes** (`cqbot`, `cqbotdo`, `cqtdo`, etc.) — written exclusively by `waterbalance.f90:integral()`; belong to soil-water-core's cumulative cohort.

## Consequences

- 12 boundary-owned instantaneous scalar fields carried in `state%soilwater` (flat layout, no cohorts — same class as ADR 0034 heat). `state%soilwater` is the first population of what will become the full soil-water record.
- No new entry-point plumbing — `boundtop`/`BoundBottom` already took `state` (SS-HEAT Task 9 windfall). The arc was pure field-carve + reader cutover.
- 8 Phase-0 config gaps closed (`sinmax/sinamp/sinave`, `cofqha/b/c`, `swcofqhc`, `hplate`). swbotb=2/4/8 modes now safe on the TOML path; coverage in regression suite remains inactive (documented).
- ~165 reads migrated across 8 files; ~68 dual-write legacy writes dropped.
- 12 legacy globals retired from `variables.f90`.
- 4 compile-surfaced hidden readers resolved during Phase 2.7 (macropore, waterbalance, config_to_variables, soilhydraulics diagnostics).
- Mini-sim writeback at `swapoutput.f90:3745–3820` retargeted for `qbot`; `gwl` and `pond` arms unchanged pending soil-water-core arc.
- pFUnit: 685 passing, 0 failures, 1 disabled throughout.
- check-full: 5/5 byte-identical at every commit.

### Known residual issues

- **`pond` / `gwl` / `kmean`** remain as legacy globals; plumbed in soil-water-core arc (#8).
- **swbotb=2 sine, swbotb=4 exp, swbotb=8 lysimeter paths** uncovered in TOML regression suite. Phase 0 config fields are correct but untested at the integration level.
- **`swbotb=-2` runtime config mutation** at `boundbottom.f90:96` (oven-dry fallback) documented but not resolved; revisit in config-vs-state hardening arc.

## References

- Discovery: `docs/superpowers/specs/2026-05-10-state-migration-boundary-discovery.md`
- Design: `docs/superpowers/specs/2026-05-10-state-migration-boundary-design.md`
- Plan: `docs/superpowers/plans/2026-05-10-boundary-state-migration.md`
- Predecessors: ADR 0030 (surfacewater pilot), ADR 0031 (drainage), ADR 0032 (solute), ADR 0033 (cumulative reset cohorts), ADR 0034 (heat — flat-layout precedent)
- Phase 0: `38004c6` (sinmax/amp/ave) → `3071b70` (cofqha/b/c + swcofqhc) → `e857103` (hplate) → `5f5d2ed` (coverage-gap doc)
- Phase 1: `897171a` (soilwater_state_t) → `cc1b9a3` (soilwater_init) → `5dfa294` (home dual-write) → `e52938d` (co-writer dual-write)
- Phase 2: `e22025d` (soilhydraulics) → `bcb990c` (waterbalance) → `fbbc230` (solute/agetracer) → `fdaf2f8` (surfacewater/drainage) → `5a9084b` (frozencond) → `8db39d5` (output + mini-sim) → `6f06046` (drop dual-write + retire globals)
