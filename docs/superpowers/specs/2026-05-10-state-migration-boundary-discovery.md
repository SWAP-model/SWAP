# Subsystem Migration Discovery: Boundary Conditions

**Date:** 2026-05-10
**Status:** discovery (read-only inventory)
**Migration #:** 5 of N — first of the four coupling-surface arcs (boundary → crop-uptake → atmosphere → soil-water core)
**Branch:** `development`
**Scope note:** This arc carves the soil-water fields touching the **top and bottom boundary surfaces** out of legacy `variables.f90` into a minimal `state%soilwater`. Crop-uptake (#6), atmosphere (#7), and soil-water-core (#8) grow that record incrementally.
**Predecessor playbooks:**
- `docs/superpowers/specs/state-migration-playbook.md` (living document)
- `docs/superpowers/specs/2026-05-10-state-migration-soilwater-discovery.md` (superseded mega-discovery — Section 2 owned globals + Section 7 hazards remain authoritative)
- `docs/superpowers/specs/2026-05-10-state-migration-heat-discovery.md` (style/depth template)
- ADRs 0030 (surfacewater pilot), 0031 (drainage), 0032 (solute), 0033 (cumulative cohorts), 0034 (heat — flat-layout precedent)

> Read-only discovery: no code changes made. All file:line references are anchors for the design phase.

---

## 1. Big picture

### Subsystem role

`boundtop` and `boundbottom` are the **coupling surfaces** of the Richards-equation soil column with the outside world:

- **`boundtop(state)`** (`src/boundary/boundtop.f90`, 310 LoC): top-surface flux assignment.
  - Computes maximum Darcy-limited evaporation `Emax` from a `hAtm` air pressure head; sets reduced evaporation `reva`.
  - Composes `q0 = (nraidt+nird+melt)*(1-ArMpSs) + runon - reva` (net top flux from atmosphere + surface water).
  - Decides between **flux boundary** (`ftoph=.false.` → assigns `qtop`) and **pressure-head boundary** (`ftoph=.true.` → assigns `pond`, `hsurf`, `kmean(1)`).
  - In ponding mode, computes `h0max` (max ponding) and may set `FlRunoff` and seed lateral macropore inflow `QMpLatSs`.
- **`PONDRUNOFF(state)`** (same file): given `h0max` from `boundtop`, computes final `pond`, `runots`, and `hsurf`. Iterative bisection when runoff is non-linear (rsroexp ≠ 1).
- **`BoundBottom(state)`** (`src/boundary/boundbottom.f90`, 184 LoC): bottom-boundary flux/head assignment.
  - Switches on `swbotb ∈ {1..8, ±2}`: prescribed gwl (1), regional flux (sine/table, ±2), Cauchy with deep aquifer (3), q(h) exponential/table (4), prescribed head (5), zero flux (6), free drainage (7), lysimeter (8).
  - Writes `qbot`, `gwlinp` (in swbotb=1), `deepgw` (in swbotb=3), `hbot` + `kmean(numnod+1)` (in swbotb=5), and a snapshot `qbot_nonfrozen = qbot` (frost feedback handoff).

### Lines of code (home files)

```
   0  src/boundary/boundary_constants.f90   (empty placeholder)
 184  src/boundary/boundbottom.f90
 310  src/boundary/boundtop.f90
 494  total
```

About 8× smaller than soil-water-core in LoC. The boundary subsystem is **physics-thin**: most of its content is small flux assignments per boundary-type switch. The arc work is therefore concentrated in field-carving and state plumbing, not algorithmic refactor.

### Entry points (call graph from outside)

| # | Call site | File:line | Phase | Currently takes `state`? |
|---|-----------|-----------|-------|--------------------------|
| 1 | `BoundBottom(state)` | `swap.f90:303` | Top of each Richards timestep (before dt-iteration loop) | **YES** (SS-HEAT Task 9) |
| 2 | `boundtop(state)` | `soilhydraulics.f90:204` (`headcalc` init), `soilhydraulics.f90:485` (`headcalc` per-iteration retry) | **Inside** Richards Newton iteration in `headcalc` | **YES** (SS-HEAT Task 9) |
| 3 | `PONDRUNOFF(state)` | called from `boundtop` only (internal sequencing via `headcalc` block) | — | **YES** |
| 4 | `BoundBottom` mini-sim call | `swapoutput.f90:3765` — **commented out** today (loop instead passes `gwlinp` directly). | — | n/a |

**Signature status:** All three public entry points already take `state intent(in)` (from SS-HEAT Task 9). No new plumbing needed for the routines themselves — the migration is field-carve + caller updates, not signature surgery.

**Critical call-order observation:** `boundtop` runs **inside `headcalc`** (the Richards iterative solver in `soilhydraulics.f90`), not at the top of the timestep. It is called once on initial F/dFdh setup (`soilhydraulics.f90:204`) and again on iterative retry when ponding/flux switch direction flips (`soilhydraulics.f90:485`). This means `boundtop`'s writes to `qtop / pond / kmean(1) / reva / hsurf / ftoph / FlRunoff / QMpLatSs` happen inside the inner loop — and `soilhydraulics.f90:headcalc` re-reads them within the same loop iteration. There is no "after boundtop" handoff seam at the swap.f90 level.

By contrast, `BoundBottom` runs at swap.f90:303 — outside the dt-reduce loop, before `headcalc` is entered. Its writes are stable across the iteration.

---

## 2. Owned-and-touched globals (the inventory)

The discovery distinguishes (a) globals **owned** by boundary (boundary writes authoritatively each call) vs (b) globals **co-written** by boundary and soil-water (boundary contributes; soil-water mostly recomputes).

### 2a. Owned globals (boundary writes authoritatively)

These are the strongest candidates to live on `state%soilwater` carved out by this arc.

| Variable | Type / shape | `variables.f90` line | Owner-of-record | Reset cadence | Activity gate | Write sites |
|---|---|---|---|---|---|---|
| `qtop` | real(8) scalar | 921 | boundary (top) | Instantaneous (per Richards iter) | `swtopb=1` (flux/head switch) — always called | `boundtop.f90:165`; **co-write (init reset)**: `soilhydraulics.f90:878` (`qtop=0` in `SoilWater(1)`); **co-write (intermediate update)**: `soilhydraulics.f90:637` (`qtop = -kmean(1)*((hsurf-h(1))/disnod(1)+1)` in headcalc inner loop, only when `ftoph=.true.`) |
| `qbot` | real(8) scalar | 878 | boundary (bottom) | Instantaneous | `swbotb ∈ {1..8, ±2}` | `boundbottom.f90:104,107,112,143,146,152,153,155,171,174,177`; **co-write**: `soilhydraulics.f90:122,252,254,257,266,274,534,537,541,552,560,721` (during headcalc setup for various swbotb modes); `frozencond.f90:216,243,287` (FrozenBounds frost feedback); `swapoutput.f90:3818` (mini-sim writeback restore) |
| `qbot_nonfrozen` | real(8) scalar | 880 | boundary (bottom) | Instantaneous (snapshot) | always | `boundbottom.f90:179` — sole writer. Read by `frozencond.f90:216` only. Pure handoff variable. |
| `pond` | real(8) scalar | 872 | shared (boundary writes outside iteration; soilhydraulics reads + restores; tillage co-writes) | Instantaneous | always | `boundtop.f90:144,163,181,254,262,270,286,304`; **co-write**: `soilhydraulics.f90:758,1031,1033,1067,1244` (init + state-restore on dt-reduce); `tillage.f90:301,327` (tillage compaction redistribution); `swapoutput.f90:3820` (mini-sim restore) |
| `gwlinp` | real(8) scalar | 800 | boundary (bottom, swbotb=1 mode) | Instantaneous | `swbotb=1` | `boundbottom.f90:76`; **co-write**: `soilhydraulics.f90:149` (clamp `gwlinp = z(NN)` if too close to node); `swapoutput.f90:3763` (mini-sim perturbation arm) |
| `deepgw` | real(8) scalar | 782 | boundary (bottom, swbotb=3 mode) | Instantaneous | `swbotb=3` | `boundbottom.f90:121,123` — sole writers. Read by `boundbottom.f90:143` and `soilhydraulics.f90:252,254,534,537`. |
| `hbot` | real(8) scalar | 807 | boundary (bottom, swbotb=5 + lysimeter setup) | Instantaneous | `swbotb=5` or swbotb=8 (lysimeter) | `boundbottom.f90:161`; `soilhydraulics.f90:155,271,557` (set in soil-side bottom-boundary setup) |
| `reva` | real(8) scalar | 924 | boundary (top) | Instantaneous (per Richards iter via boundtop reset) | always | `boundtop.f90:126,128` — sole writers. Read at `soilhydraulics.f90:109,111`, `waterbalance.f90:396` (revats = reva*dt for cumulative ievap accumulation), `et.f90:` reads via legacy. |
| `hsurf` | real(8) scalar | 816 | boundary (top) | Instantaneous (per iter) | `ftoph=.true.` mostly | `boundtop.f90:142,162,255,263,273,291,306`; read at `soilhydraulics.f90:214,335,637`. |
| `ftoph` | logical scalar | 985 | boundary (top) | Instantaneous (per iter) | always | `boundtop.f90:141,160,167`; read by `soilhydraulics.f90:213,335,637` (gates pressure-head vs flux branch of F/dFdh). |
| `runots` | real(8) scalar | 944 | boundary (top, PONDRUNOFF) | Instantaneous (per iter) | always | `boundtop.f90:145,164,253,259,272,287,305`; **co-write**: read+accumulated in `waterbalance.f90:461,491,493` (intermediate `iruno`, cumulative `crunoff`/`cinund`); `surfacewater.f90:451,480,486,511,546,566,663,699` (reservoir balance compute). |
| `FlRunoff` | logical scalar | 636 | boundary (top) | Instantaneous (per iter) | always | `boundtop.f90:101,169`; read by `soilhydraulics.f90:210` and `boundtop.f90:228` (PONDRUNOFF gate). |
| `QMpLatSs` | real(8) scalar | 1191 | boundary (top) — macropore-coupled | Instantaneous (per iter) | `flMacroPore` | `boundtop.f90:102,182,183,184,186,236,242`; read by `macropore.f90:1544` (`cQMpLatSs = cQMpLatSs + QMpLatSs` accumulation). |

**Subtotal: 13 owned scalars + zero per-node arrays.** Boundary owns no allocatable arrays — all owned fields are scalars.

### 2b. Co-write fields (boundary writes, soil-water also writes)

These are **soil-water-owned** but boundary contributes a write. Two strategies for design:
- (i) Boundary writes via the state path (`state%soilwater%X = …`) once the field is migrated — the cleaner end state.
- (ii) Boundary returns a value through call args; soil-water assigns it. More plumbing churn.

| Variable | Type / shape | Soil-water owner write sites | Boundary write sites | Recommendation |
|---|---|---|---|---|
| `kmean(macp+1)` | real(8) array | `soilhydraulics.f90:113,133,186,189,241,262,456,474,475,523,548,1054,1056,1174,1177,1242` (Richards setup + state-restore) | `boundtop.f90:143,161,168` (top entry), `boundbottom.f90:164,166` (bottom-prescribed-head entry) | Soil-water-owned (full `kmean(:)` belongs to soil-water-core). Boundary writes only `kmean(1)` and `kmean(numnod+1)` — the two surface entries. After full migration: boundary writes through `state%soilwater%kmean(1)` and `state%soilwater%kmean(numnod+1)`. **This arc** can either (a) migrate the two entry-points as scalars `kmean_top`, `kmean_bot` on the soilwater state, or (b) defer `kmean` to the soil-water-core arc and have boundary continue to write the legacy `kmean(:)` global. Recommendation: **(b) defer** — `kmean` is a Richards-interior field; carving only entries 1 and N+1 produces an awkward split. |

**Subtotal: 1 co-write field (kmean top/bottom entries).** Defer to soil-water-core arc.

### 2c. Cumulative boundary fluxes — **NOT boundary-owned**

The cumulative bottom and top fluxes (`cqbot`, `cqbotdo`, `cqbotup`, `cqtdo`, `cqtup`, `cqprai`) are written **only by `waterbalance.f90:integral()`** (soil-water-core's `integral` routine), not by boundary code itself. They accumulate `qbot*dt` and `q(1)*dt` (the inter-compartment flux at the surface interface), not boundary-routine outputs directly.

| Cumulative field | Reset (init) | Accumulator write | Boundary write? |
|---|---|---|---|
| `cqbot` | `soilhydraulics.f90:1137` (in `flzerocumu` block) | `waterbalance.f90:511` (`cqbot = cqbot + qbotts`) | **No** |
| `cqbotdo` | `soilhydraulics.f90:1138` | `waterbalance.f90:507` | No |
| `cqbotup` | `soilhydraulics.f90:1139` | `waterbalance.f90:509` | No |
| `cqtdo` | `soilhydraulics.f90:1147` | `waterbalance.f90:530` | No |
| `cqtup` | `soilhydraulics.f90:1148` | `waterbalance.f90:532` | No |
| `cqprai` | `soilhydraulics.f90:1149` | `waterbalance.f90:527` | No |

**These belong to the soil-water-core arc** (cumulative cohort, single `flzerocumu` gate). The boundary arc does **not** touch them.

### 2d. Tabulated config inputs (read-only)

Already in typed config (`bottom_boundary_config_t` per `src/config/bottom_boundary_config.f90`):

| Legacy tabular global | Read at | Typed-config coverage |
|---|---|---|
| `gwltab(mabbc*2)` | `boundbottom.f90:76` | covered (`gwl_file` → swc_table per swbotb=1) |
| `qbotab(mabbc*2)` | `boundbottom.f90:107,146,155` | covered (`qbot2_file`, `qbot4_file`, `qhbot_file` → `qbot_table`) |
| `haqtab(mabbc*2)` | `boundbottom.f90:123` | covered (`haquif_file`) |
| `hbotab(mabbc*2)` | `boundbottom.f90:161` | covered (`hbot5_file`) |
| `pondmxtab(2*mairg)` | `boundtop.f90:224` | covered (in `soil_config_t.pondmx` + table path) |

No Phase 0 tabular gaps. The boundary arc inherits Phase 0 risk from the **scalar** parameters (Section 6).

---

## 3. External readers (5-category framework)

For each boundary-owned field (Section 2a), files outside `src/boundary/` and `src/state/` that **read** the field. Categories per playbook: (1) Output, (2) Compute, (3) Working buffer, (4) Init-seed, (5) Call-site arg.

### Per-field reader index

| Field | Output (Cat 1) | Compute (Cat 2) | Working buffer (Cat 3) | Init-seed (Cat 4) | Call-site arg (Cat 5) |
|---|---|---|---|---|---|
| `qtop` | `swapoutput.f90` (afo/bfo blocks via implicit), `swap_csv_output.f90` (QTOP column) | `soilhydraulics.f90:217,502,637` (headcalc F vector + flux update), `solute.f90:173,174,175,177` (top-boundary solute flux), `waterbalance.f90:324` (use clause + line 336 qbot back-calc fallback) | — | — | — |
| `qbot` | `swapoutput.f90:2473,2769,3228` (afo/aun/swb), `swap_csv_output.f90` (QBOT column) | `soilhydraulics.f90:122,252,…,560,721` (headcalc + state-restore), `waterbalance.f90:336,341,399,511` (cumu accumulation + checkmassbal), `solute.f90:291,292,293,295,296` (bottom solute flux), `frozencond.f90:216,243,287` (FrozenBounds frost feedback) | `swapoutput.f90:3745,3818` (mini-sim snapshot+restore) | — | — |
| `qbot_nonfrozen` | — | `frozencond.f90:216` only | — | — | — |
| `pond` | `swapoutput.f90:302,325,455,475,1093,2407,2473,2769,3228,3306,3380,3912`, `swap_csv_output.f90:243,384,470,476` | `et.f90:612,637` (potential E uses pond), `meteoday.f90:497,720,734` (interception path), `drainage.f90:666` (drainage flux), `surfacewateruts.f90:14` (qhtab), `tillage.f90:8,301,327`, `macropore.f90` (indirect), `soilhydraulics.f90:758,1031,1033,1067,1244` (init/restore), `solute.f90:174,178` (cpond), `waterbalance.f90:40` (use clause + line 169 calcgwl), `surfacewater.f90:308` (use clause) | `swapoutput.f90:3747,3820` (mini-sim snapshot+restore), `swap_csv_output.f90:384,476` (PondOld save) | `soil_config.f90:73` (default 0.0); `config_to_variables.f90` seeds from `soil.initial.pond` | — |
| `gwl` | `swapoutput.f90:3228,3723,3746,3819,3912`, `swap_csv_output.f90:242` | `drainage.f90:180,422,666` (heavy — drainage gradient), `surfacewater.f90:308`, `waterbalance.f90:40` (use clause), `soilhydraulics.f90:757,1009,1019,1028,1243` (calc+restore), `et.f90` (heavy), `frozencond.f90:277`, `boundbottom.f90:118,152,155` (within home tree), `meteoday.f90` (via use chain) | `swapoutput.f90:3746,3819` (mini-sim) | — | — |
| `gwlinp` | — | `soilhydraulics.f90:107,109,111,112,113,123,142,145,147,149,155,226,344` (headcalc swbotb=1 + dirichlet branch); `waterbalance.f90:40,169` (calcgwl gate) | `swapoutput.f90:3763` (mini-sim perturbation arm — writes `gwlinp = gwltmp + dgwl(i)`) | — | — |
| `deepgw` | — | `soilhydraulics.f90:252,254,534,537` (headcalc swbotb=3 setup) | — | — | — |
| `hbot` | — | `soilhydraulics.f90:155,228,271,557` (headcalc bottom-head setup) | — | — | — |
| `reva` | `swap_csv_output.f90` (rare) | `soilhydraulics.f90:109` (re-uses reva in init q0), `waterbalance.f90:396` (revats = reva*dt for ievap), `et.f90` (reads as the actual evap output) | — | — | — |
| `hsurf` | — | `soilhydraulics.f90:214,335,637` (headcalc top-pressure-head branch) | — | — | — |
| `ftoph` | `swapoutput.f90:1093` (final-state dump) | `soilhydraulics.f90:213,335,637,796` (gates pressure-head vs flux branch in F/dFdh) | — | — | — |
| `runots` | `swapoutput.f90`, `swap_csv_output.f90` (RUNOTS column) | `waterbalance.f90:461,490,491,492,493` (iruno + cinund/crunoff intermediate/cumu), `surfacewater.f90:309,451,480,486,511,546,566,663,699` (reservoir balance — heavy), `soilhydraulics.f90:111` (q1 in init swbotb=1 branch) | — | — | — |
| `FlRunoff` | — | `soilhydraulics.f90:210` (headcalc top-boundary mode gate) | — | — | — |
| `QMpLatSs` | `swapoutput.f90` (macropore output blocks) | `macropore.f90:1544` (cQMpLatSs cumulative accumulator); `soilhydraulics.f90` (q0 chain via integral) | — | — | — |

### Distinct external reader files (union over all owned fields)

1. `src/soil/soilhydraulics.f90` — **heaviest reader**, ~25 read sites across qtop/qbot/pond/gwl/gwlinp/deepgw/hbot/reva/hsurf/ftoph/FlRunoff. Almost all are compute sites in `headcalc`. **Co-writes qbot, pond, gwlinp, kmean.**
2. `src/soil/waterbalance.f90` — `integral`, `calcgwl`, `checkmassbal`: cumulative accumulators read qbot/qtop/runots/reva and write the cumulative cohort. **Co-writes gwl** (calcgwl).
3. `src/io/swapoutput.f90` — output (afo/aun/swb, .vap, .bal, etc.) + **mini-sim writeback (Cat 3 hazard)** at lines 3745–3820.
4. `src/io/swap_csv_output.f90` — CSV columns: QTOP, QBOT, GWL, POND, REVA, RUNOTS, FTOPH, plus PondOld working buffer.
5. `src/solute/solute.f90` — top/bottom solute flux uses qtop/qbot/pond directly.
6. `src/solute/agetracer.f90` — age-tracer surface boundary.
7. `src/heat/frozencond.f90` — **co-writes qbot** in `FrozenBounds` (already takes state).
8. `src/drainage/drainage.f90` — gwl is the primary drainage driver.
9. `src/drainage/surfacewater.f90` — runots accumulates into reservoir; gwl/pond for reservoir balance.
10. `src/utils/surfacewaterutils.f90` — qhtab path uses pond.
11. `src/atmosphere/et.f90` — reads pond, reva.
12. `src/atmosphere/meteoday.f90` — interception path reads pond.
13. `src/crop/tillage.f90` — **co-writes pond** in compaction redistribution; reads pond.
14. `src/macropore/macropore.f90` — reads QMpLatSs (cumulative accumulator).

**Total external reader files: 14.** Compare with heat (15), surfacewater pilot (12).

**Per-file read-site count (rough):**
- soilhydraulics.f90: ~25 sites
- swapoutput.f90: ~15 sites
- waterbalance.f90: ~10 sites
- swap_csv_output.f90: ~6 sites
- surfacewater.f90: ~10 sites
- solute.f90: ~7 sites
- drainage.f90: ~5 sites
- All others: 1–4 sites

**Cat 3 hazard surface:** `swapoutput.f90:3745–3820` snapshots qbot/gwl/pond/theta/h, runs a mini-sim with `gwlinp = gwltmp + dgwl(i)` perturbations, then **writes back** all five fields. After migration, this writeback targets `state%soilwater%qbot`, `state%soilwater%pond`, `state%soilwater%gwl`. Same pattern as the drainage mini-sim (already resolved in SS-DRST).

---

## 4. Co-writers

Files outside the boundary home tree that **write** Section 2 fields. Each is a Phase-1-dual-write or Phase-2-migration target.

| Co-writer file | Field(s) written | Context | State-in-scope today? |
|---|---|---|---|
| `src/soil/soilhydraulics.f90` | `qbot`, `qtop`, `pond`, `kmean(1)`, `kmean(numnod+1)`, `hbot`, `gwlinp` | `SoilWater(1)` init reset (lines 870–890), `headcalc` per-iter setup (lines 100–280, 480–560), `headcalc` flux update (lines 637, 721), `SoilWaterStateVar` restore (lines 1242–1244), `SoilWater(3)` integrand setup (line 1067) | **YES** — `soilwater(task,state)`, `headcalc(state)` already plumbed. Co-writes for qbot/qtop/pond will become `state%soilwater%X = …` once those fields migrate. |
| `src/soil/waterbalance.f90` | `gwl` (calcgwl) | `calcgwl` recomputes gwl from h profile each iteration | **NO** — `calcgwl()` is plain no-arg subroutine using bare `use variables`. Needs plumbing if gwl migrates. |
| `src/heat/frozencond.f90` | `qbot` | `FrozenBounds(state)` overrides qbot under frost (zero or redistribute to deepest drain) | **YES** — already takes state%drainage. Co-write of qbot is a borrower-write into a soil-water field. |
| `src/crop/tillage.f90` | `pond` | `Adapt_WC_H` redistributes water-content after tillage compaction; if excess > saturated, adds to pond | **NO** — bare `use variables`. Needs state plumbing if pond migrates. Same pattern as soil-water-core arc (tillage hazard). |
| `src/io/swapoutput.f90` | `qbot`, `gwl`, `pond`, `gwlinp` | Mini-sim snapshot+restore for stocoav/stocot1 output (lines 3745–3820) and `gwlinp = gwltmp + dgwl(i)` perturbation arm (line 3763) | **YES** — `state_om` and `state_main` both passed. Mini-sim writeback target updates trivially. |

**Total co-writers: 5 distinct files.** Compare heat (1: soilhydraulics rfcp reset), drainage (~3), surfacewater (~3).

**State-plumbing gap:** `calcgwl()` (in waterbalance.f90) and `tillage.f90` need `state` arg added if their respective targets (`gwl`, `pond`) migrate in this arc. Both are inherited from the soil-water mega-discovery hazards — handled below.

---

## 5. Init-order analysis

Reproduce the heat-arc lesson: an `allocated()` guard is needed when a state-allocator (e.g. `soilwater_init`) runs **after** a routine that reads the state field. The boundary arc inverts the heat layout — boundary routines run **after** SoilWater(1) but the soil-water state will be allocated INSIDE SoilWater(1), so the order is favorable.

### Init-order map (from swap.f90)

```
swap.f90:182   call CalcGrid()                    ! grid dims (numnod, dz, z, …)
swap.f90:184   if (flTillage) call DoTillage(1)   ! tillage init — writes pond, cofgen
swap.f90:188   call SoilWater(1, state)           ! INIT POINT for soilwater_state allocation
swap.f90:195   call drainage_init(state, config)
swap.f90:197   call heat_init(state)              ! allocates state%heat%rfcp etc
swap.f90:206   if (flTemperature) call Temperature(1, state)
swap.f90:303   call BoundBottom(state)            ! per timestep
swap.f90:317   call SoilWater(2, state)           ! contains call to headcalc(state)
                  └─ soilhydraulics.f90:204  call boundtop(state)   ! per Richards iter
                  └─ soilhydraulics.f90:485  call boundtop(state)   ! per retry
```

**Allocation site: inside `SoilWater(1, state)`** (analogous to heat). The soil-water state must be allocated before any of {BoundBottom, boundtop, SoilWater(2), tillage(2/3), heat, anyone reading qbot/pond/gwl/qtop} runs.

### Risk: `flTillage` runs **before** `SoilWater(1)`

`DoTillage(1)` at line 184 runs before `SoilWater(1)` at line 188. If `pond` is migrated to `state%soilwater%pond`, then `tillage.f90:301,327` writing pond would need `state%soilwater` already allocated. Resolution options:

1. **Defer pond writes in DoTillage(1) init path** if state not yet allocated — guarded by `if (allocated(state%soilwater%fields))`. Mirrors the heat `allocated()` guard pattern.
2. **Move the allocation up** — create a `soilwater_init(state)` call before `DoTillage(1)`, paralleling `heat_init` placement. The cleaner end state.
3. **Don't migrate pond in this arc** — leave it as legacy global; carve only the bottom-boundary fields. Smaller scope.

**Recommendation:** option (2) — introduce a `soilwater_init(state)` at the earliest natural point (immediately after `CalcGrid()`, before `DoTillage(1)`). This matches the heat pattern.

### Where should `soilwater_init(state)` go?

Insert at `swap.f90:183` (between `CalcGrid()` and `DoTillage(1)`). It needs `numnod` from CalcGrid for per-node arrays — but the boundary arc owns only scalars, so allocation is trivial (no `numnod`-dependent arrays in this arc's scope). Could even go before `CalcGrid()`.

---

## 6. Config / Phase 0 candidates

The bottom-boundary config (`bottom_boundary_config_t` in `src/config/bottom_boundary_config.f90`) is mature — covers `swbotb`, `sw2`, `sw3`, `sw4`, `swqhbot`, `shape`, `hdrain`, `rimlay`, `aqave/amp/per/tmax`, files for all swbotb modes. But several **scalar physics parameters** read by `boundbottom.f90` are NOT in typed config and NOT populated by the TOML adapter:

| Legacy global | `variables.f90` line | Used at | Typed-config field? | Adapter populates? |
|---|---|---|---|---|
| `sinmax` | 949 | `boundbottom.f90:104` (swbotb=2 sine) | **NO** | **NO** |
| `sinamp` | 947 | `boundbottom.f90:104` (swbotb=2 sine) | **NO** | **NO** |
| `sinave` | 948 | `boundbottom.f90:104` (swbotb=2 sine) | **NO** | **NO** |
| `cofqha` | 761 | `boundbottom.f90:152` (swbotb=4 swqhbot=1 exp) | **NO** | **NO** (note: drainage config has a `cofqha_table`; that is the swsrf=4 surface case, not boundbottom swbotb=4) |
| `cofqhb` | 762 | `boundbottom.f90:152` (swbotb=4 swqhbot=1 exp) | **NO** | **NO** |
| `cofqhc` | 763 | `boundbottom.f90:153` (swbotb=4 + swcofqhc=1) | **NO** | **NO** |
| `swcofqhc` | — (gate) | `boundbottom.f90:153` | **NO** | **NO** |
| `hplate` | 813 | `soilhydraulics.f90:271,557` (swbotb=8 lysimeter) | **NO** | **NO** |
| `gwlconv` | 798 | `soilhydraulics.f90` (headcalc convergence) | covered (`simulation.numerical.gwlconv` → `config_to_variables.f90:117`) | **YES** |
| `SwBotb3ResVert` | — (gate) | `boundbottom.f90:127,139` (swbotb=3 vertical-resistance switch) | check — likely in bottom_boundary_config | **probably yes** (it's a swbotb=3 sub-switch) |

**Phase 0 candidates for THIS arc:** 7–8 fields (sinmax/amp/ave, cofqha/b/c+swcofqhc, hplate). The `cofqha_table` in drainage config is a different field (surface-runoff qhtab, not bottom-boundary). These are **adapter gaps** identical in shape to the heat-arc analytical-method gaps (ddamp/tmean/tampli/timref).

Boundary arc Phase 0 work:
1. Add `sinmax`, `sinamp`, `sinave` to `bottom_boundary_config_t` (gated by `swbotb=2 .and. sw2=1`).
2. Add `cofqha`, `cofqhb`, `cofqhc`, `swcofqhc` to `bottom_boundary_config_t` (gated by `swbotb=4 .and. swqhbot=1`).
3. Add `hplate` — but `hplate` is used in `swbotb=8` (lysimeter) inside `soilhydraulics.f90`, not in `boundbottom.f90`. It could go in `bottom_boundary_config_t` (gated by `swbotb=8`) or in soil config. Recommendation: `bottom_boundary_config_t` (semantically a bottom-boundary parameter).
4. TOML adapter (`config_to_variables.f90`) writes these legacy globals from typed config.
5. Add unit tests + at least one regression case exercising the swbotb=2 sine + swbotb=4 exp paths if not already covered.

**Test coverage check:** The five TOML regression cases (1–5) use which `swbotb`?

Looking at the soil-water mega-discovery, this was not enumerated explicitly. The design phase should grep the TOML fixtures to confirm — but it is plausible that **swbotb=2 sine** and **swbotb=4** are uncovered in TOML regression (analogous to heat swcalt=1).

### Phase 0 candidates NOT in scope for boundary arc

`gwlconv` is already covered. `hplate` lives in soilhydraulics, not boundary — but it is a bottom-boundary parameter conceptually; include in this arc's Phase 0 to retire it cleanly.

---

## 7. Known coupling hazards (boundary-specific)

### Hazard #1 — Mini-sim writeback (qbot/gwl/pond)

**Scope: small.** `swapoutput.f90:3745–3820` snapshots qbot/gwl/pond + theta/h, runs a 2-arm perturbation mini-sim (`gwlinp = gwltmp ± 1 cm`), restores. Same pattern as drainage's resolved mini-sim. After migration: writeback targets `state%soilwater%qbot/gwl/pond`. Mechanical.

### Hazard #2 — Multi-site qbot writes (frozencond + soilhydraulics + boundbottom)

**Scope: medium.** `qbot` has three writer trees:
- `boundbottom.f90` (11 sites) — authoritative per swbotb mode.
- `soilhydraulics.f90` (12 sites in lines 122/252/254/257/266/274/534/537/541/552/560/721) — these are **duplicate writes** inside `headcalc` that re-derive qbot from the current iterate's h/kmean. They mirror `boundbottom.f90`'s logic but use intermediate iteration values.
- `frozencond.f90:216,243,287` — frost override (FrozenBounds, already plumbed).

Question for design phase: **why does soilhydraulics duplicate boundbottom's qbot logic inside headcalc?** Hypothesis: boundbottom runs once per dt iteration (before headcalc) using prior h/kmean; soilhydraulics then re-uses iteration-current h/kmean for the iterative Jacobian. This is functionally a coupling-surface duplication. The boundary arc should NOT try to unify these — leave the duplicate writes in place; just have both site sets target `state%soilwater%qbot`.

### Hazard #3 — `boundtop` called inside `headcalc`

**Scope: medium.** Unlike `BoundBottom` (called at swap.f90:303 once per timestep), `boundtop` is called at `soilhydraulics.f90:204` (per-iter setup) and `soilhydraulics.f90:485` (per-iter retry). This means `boundtop`'s writes to qtop/pond/kmean(1)/reva/hsurf/ftoph/FlRunoff are **per-Richards-iteration**. Performance: irrelevant (boundtop is cheap). Semantic: `pond` written by `boundtop` is a transient per-iteration value, not a "final pond at end of timestep" — `pond` is then re-overwritten by `soilhydraulics.f90:1031,1033` in `SoilWater(3)` and again by `soilhydraulics.f90:1244` in `SoilWaterStateVar` restore. The state field's last-writer-wins semantics carry through unchanged.

### Hazard #4 — `tillage.f90` co-writes `pond` and is NOT state-plumbed

**Scope: medium.** `tillage.f90:301,327` writes `pond` in `Adapt_WC_H`. `tillage.f90` uses bare `use variables` and has no `state` arg. If `pond` migrates in this arc, `tillage.f90` needs plumbing. This is the same hazard flagged in the soil-water mega-discovery (Hazard #6). Options:

- (a) Plumb state into `DoTillage(itask, state)` as part of this arc.
- (b) Don't migrate `pond` in this arc — defer to soil-water-core where tillage plumbing happens.

Recommendation: **(b) defer pond**. The boundary arc's pond writes are inside the Richards iteration; carving pond as a standalone migration without also touching tillage and soilhydraulics's restore path produces too many half-migrated edges. Leave pond as legacy global for this arc; migrate it in the soil-water-core arc alongside h/theta/cofgen.

This narrows the boundary arc to: **qbot, qbot_nonfrozen, qtop, hbot, gwlinp, deepgw, reva, hsurf, ftoph, runots, FlRunoff, QMpLatSs**. That is the minimum-blast-radius set.

### Hazard #5 — `calcgwl()` writes `gwl` and is NOT state-plumbed

**Scope: medium.** `waterbalance.f90:calcgwl()` (line 39, no args) writes gwl from h+pond profile each iteration. If `gwl` migrates, `calcgwl` needs `state` plumbing. Same defer-to-soil-water-core logic as Hazard #4. Recommendation: **leave gwl as legacy global** in this arc; migrate alongside gwlflcpzo/nodgwl/pegwl in soil-water-core (they all share the calcgwl call site).

### Hazard #6 — `boundbottom.f90:96` mutates `swbotb` (logical -2 fallback)

**Scope: small but notable.** When swbotb=2 detects oven-dry conditions, it sets `swbotb = -2` permanently (the change persists for the rest of the simulation). `swbotb` is technically config, not state, but this is a **runtime mutation of a config variable** that re-enters the swbotb branch logic. Future config-vs-state separation should preserve this fallback path — either (a) move `swbotb_effective` into `state%soilwater` or `state%boundary`, or (b) keep `swbotb` as legacy global with a noted mutation. The boundary arc must not break this fallback. Recommendation: leave it alone for this arc; document in ADR.

### Hazard #7 — `cqMpLatSs` in soilhydraulics.f90:888 init reset

**Scope: small.** `cQMpLatSs` is declared at variables.f90:1191 (after macropore-section, grouped with QMpLatSs). It is **not in this arc's owned set** (it's a macropore cumulative, written at `macropore.f90:1544`), but it is initialized by `soilhydraulics.f90:888` (`cQMpLatSs = 0.0d0` in `SoilWater(1)` init). Recommendation: not in boundary arc; flag for macropore arc.

### Hazard #8 — Phase 0 config gaps (sinmax/amp/ave, cofqha/b/c, hplate)

**Scope: medium.** Recapping Section 6 — 7+ scalar bottom-boundary parameters are missing from typed config and the TOML adapter. The boundary arc must close these gaps (analogous to heat's ddamp/tmean/tampli/timref).

---

## 8. Reset-cadence analysis (cohort vs flat decision)

All Section 2a owned fields are **instantaneous** (computed afresh each Richards iteration or each `BoundBottom` call). None are intermediate or cumulative.

| Cadence | Count | Cohort target |
|---|---|---|
| Instantaneous | 13 | Flat fields on `soilwater_state_t` |
| Intermediate | 0 | — |
| Cumulative | 0 | — |
| Per-day | 0 | — |

**Decision: FLAT.** For the boundary arc, `soilwater_state_t` adds 13 scalar flat fields. No cohort sub-records introduced in this arc. The soil-water-core arc will add the cumulative cohort (cqbot, cqtop, etc.) and intermediate cohort (inq, iqrot, etc.) later.

Compare:
- Heat: 13 owned, all instantaneous, flat → ADR 0034 precedent.
- Surfacewater: 25 owned, mixed I+M+C → cohort partitioning per ADR 0033.
- Drainage: similar mixed → cohort.
- Boundary: 13 owned, all instantaneous, flat → **matches heat layout**.

**Recommended `soilwater_state_t` initial shape (boundary arc only):**

```fortran
type :: soilwater_state_t

   ! Top-boundary fluxes / surface variables (from boundtop)
   real(real64) :: qtop      = 0.0_real64   !! top-surface flux (cm/d)
   real(real64) :: reva      = 0.0_real64   !! actual soil evaporation (cm/d)
   real(real64) :: hsurf     = 0.0_real64   !! pressure head at surface (cm)
   real(real64) :: runots    = 0.0_real64   !! runoff this step (cm)
   real(real64) :: qmplatss  = 0.0_real64   !! lateral macropore inflow (cm/d)
   logical      :: ftoph     = .false.      !! flag: pressure-head top boundary
   logical      :: flrunoff  = .false.      !! flag: runoff potential

   ! Bottom-boundary fluxes / variables (from BoundBottom)
   real(real64) :: qbot          = 0.0_real64   !! bottom flux (cm/d)
   real(real64) :: qbot_nonfrozen = 0.0_real64  !! bottom flux pre-frost snapshot
   real(real64) :: hbot          = 0.0_real64   !! prescribed head at bottom (cm)
   real(real64) :: gwlinp        = 0.0_real64   !! prescribed gwl, swbotb=1 (cm)
   real(real64) :: deepgw        = 0.0_real64   !! deep-aquifer head, swbotb=3 (cm)

   ! NOT in this arc (deferred to soil-water-core):
   !   pond, gwl, kmean(:), theta(:), h(:), q(:), inq(:), cohort sub-records.

end type soilwater_state_t
```

12 fields. `pond` deferred per Hazard #4 (avoids tillage plumbing). `gwl` deferred per Hazard #5 (avoids calcgwl plumbing).

**Open question:** drop `qbot_nonfrozen` and replace with a local variable in `FrozenBounds`? It is a one-write-one-read snapshot. Recommendation: include in state for symmetry; revisit if `FrozenBounds` is refactored.

---

## 9. Scope estimate

| Metric | Boundary | Heat (ADR 0034) | Surfacewater (pilot) |
|---|---|---|---|
| Owned globals | 12 (all instantaneous, all scalars) | 12 (all instantaneous, mix of scalar + array) | ~25 (mix I + M + 2 C cohorts) |
| Co-write fields | 1 (kmean top/bottom — deferred to soil-water-core) | 1 (rfcp reset by soilhydraulics) | 3 |
| External reader files | 14 | 15 | 12 |
| External read sites | ~80 | ~50 | ~80 |
| Phase 0 config candidates | 7–8 (sinmax/amp/ave, cofqha/b/c, swcofqhc, hplate) | 6 (analytical method + tables) | 2 |
| Subsystems needing state plumbing | 0 (boundtop and BoundBottom already plumbed) | 3 (Temperature, FrozenCond, FrozenBounds) | 4 (SurfaceWater tasks) |
| Cohort design complexity | 0 cohorts (flat) | 0 cohorts (flat) | 3 cohorts post-ADR 0033 |
| LoC home tree | 494 | 1119 (heat + frozencond + config + reader) | similar |

### Suggested task decomposition (boundary arc)

**Phase 0 — config gaps**

1. **B-0.1**: Add `sinmax`, `sinamp`, `sinave` to `bottom_boundary_config_t` (gate: `swbotb=2 .and. sw2=1`). TOML reader update. pFUnit tests.
2. **B-0.2**: Add `cofqha`, `cofqhb`, `cofqhc`, `swcofqhc` to `bottom_boundary_config_t` (gate: `swbotb=4 .and. swqhbot=1`). TOML reader update. pFUnit tests.
3. **B-0.3**: Add `hplate` to `bottom_boundary_config_t` (gate: `swbotb=8`). TOML reader update. pFUnit test.
4. **B-0.4**: `config_to_variables.f90` writes legacy globals from new typed-config fields.
5. **B-0.5**: Regression-case audit — confirm swbotb=2 sine + swbotb=4 + swbotb=8 covered, or add minimal fixtures.

**Phase 1 — state-type design + threading**

6. **B-1.1**: Add `soilwater_state_mod` with the 12-field flat type (Section 8 layout). Add `state%soilwater` to `swap_state_t` aggregator. pFUnit init/zero test.
7. **B-1.2**: Add `soilwater_init(state)` (or fold init into a `state%soilwater%init()` constructor). Call at `swap.f90:183` (before `DoTillage(1)`).
8. **B-1.3**: Dual-write inside boundary home tree — `boundtop.f90` and `boundbottom.f90` write both legacy global and state field for the 12 owned fields. Verify regression numerical bit-equality.
9. **B-1.4**: Dual-write at soilhydraulics co-write sites for qbot (lines 122, 252–560, 721) and qtop (line 637). FrozenBounds qbot dual-write.

**Phase 2 — external reader migration**

10. **B-2.1**: `soilhydraulics.f90` cut over all qbot/qtop/reva/hsurf/ftoph/FlRunoff/hbot/gwlinp/deepgw reads to `state%soilwater%X`.
11. **B-2.2**: `waterbalance.f90` cut over (reva, runots, qbot reads in integral).
12. **B-2.3**: `solute.f90`, `agetracer.f90` cut over (qtop, qbot, runots).
13. **B-2.4**: `surfacewater.f90`, `drainage.f90`, `surfacewaterutils.f90` cut over (runots; qtop indirect).
14. **B-2.5**: `frozencond.f90:FrozenBounds` qbot/qbot_nonfrozen read/write through state.
15. **B-2.6**: `swapoutput.f90` + `swap_csv_output.f90` output blocks cut over. Mini-sim writeback (lines 3745–3820) targets state. **Hazard #1 resolution.**
16. **B-2.7**: Retire legacy globals (`qtop`, `qbot`, `qbot_nonfrozen`, `hbot`, `gwlinp`, `deepgw`, `reva`, `hsurf`, `ftoph`, `runots`, `FlRunoff`, `QMpLatSs`) from variables.f90.

**Total tasks: ~16** (5 Phase 0 + 4 Phase 1 + 7 Phase 2). Comparable to heat (~14 tasks). The arc fits one focused work session per the playbook cadence.

---

## 10. Open questions for design phase

1. **Cohort vs flat shape.** **Recommendation: FLAT.** All 12 owned fields are instantaneous scalars. No cohort sub-records needed in this arc. Future arcs (soil-water-core) will add cumulative + intermediate cohorts.

2. **`pond` ownership — boundary arc or soil-water-core arc?** **Recommendation: defer to soil-water-core.** Tillage co-writes pond and is not state-plumbed; carving pond now would force tillage plumbing too. Defer.

3. **`gwl` ownership — boundary arc or soil-water-core arc?** **Recommendation: defer to soil-water-core.** `calcgwl()` co-writes gwl and is not state-plumbed. Defer.

4. **`kmean(1)` / `kmean(numnod+1)` — carve as scalars `kmean_top` / `kmean_bot`, or defer with whole `kmean(:)` to soil-water-core?** **Recommendation: defer.** Carving two array entries is awkward; whole `kmean` belongs to soil-water-core.

5. **`qbot_nonfrozen` — state field or local in FrozenBounds?** **Recommendation: include in state.** Symmetric with qbot; cheap.

6. **Mini-sim writeback (swapoutput.f90:3745–3820) — refactor or preserve?** **Recommendation: preserve.** Same pattern as drainage's resolved mini-sim. Just retarget writeback to `state%soilwater%X`.

7. **`swbotb = -2` runtime mutation (oven-dry fallback) — relocate to state or keep as config mutation?** **Recommendation: keep as config mutation for now.** Flag in ADR. Revisit when config-vs-state separation is hardened.

8. **Phase 0 regression coverage — do the five TOML cases exercise swbotb=2 sine, swbotb=4 exp, swbotb=8?** **Verify in design phase.** Grep TOML fixtures; if uncovered, add minimal regression fixture.

9. **`reva` — boundary-owned or atmosphere-owned?** It is the **actual** soil evaporation (boundary-computed, capped by Darcy `Emax`). The atmosphere arc owns `peva` (potential) and `empreva` (empirical reduced). `reva` is best classified as boundary-owned because it is the result of the top-boundary's hydraulic-cap decision. **Recommendation: boundary arc owns reva.**

10. **`hbot` write at `soilhydraulics.f90:155` and `:271,557` — boundary co-write or boundary's responsibility?** Both lines set `hbot` from `gwlinp` or `hplate` for the bottom-Dirichlet branch. The semantic owner is the bottom-boundary subsystem; the writes happen inside `headcalc`. **Recommendation: state field is boundary-owned; soilhydraulics writes through `state%soilwater%hbot`.** Same handling as qbot.

11. **`runots` reset cadence — is it instantaneous or intermediate?** `boundtop.f90` writes `runots` per Richards iteration (multiple times per timestep). `waterbalance.f90:461` accumulates `iruno = iruno + runots` then `:493` accumulates `crunoff = crunoff + runots`. So `runots` is **per-iteration instantaneous**, then sampled into intermediate/cumulative. **Recommendation: classify as instantaneous; the iruno/crunoff cumulatives belong to soil-water-core arc.**

12. **Should `soilwater_init(state)` go before or after `CalcGrid()`?** Since this arc's fields are all scalars, `numnod` is not needed. But future arcs add per-node arrays. **Recommendation: place after `CalcGrid()`** for forward compatibility.

13. **TOML fixtures — does the swbotb=8 lysimeter path exist in `.swp` regression?** Tests would tell us if `hplate` is exercised. **Verify in design phase.**

---

## Discovery summary

- **Owned set: 12 scalars.** All instantaneous. No cohort partitioning needed.
- **Co-writers: 5 files**, of which `tillage.f90` and `calcgwl()` are not state-plumbed (relevant only if pond/gwl migrate — both deferred).
- **External readers: 14 files**, ~80 read sites. Heaviest: soilhydraulics, waterbalance, swapoutput.
- **Signature status: BOTH boundary routines already take `state`** (SS-HEAT Task 9). No new plumbing required for entry points.
- **Phase 0 gaps: 7–8 scalar parameters** (sinmax/amp/ave, cofqha/b/c, swcofqhc, hplate) missing from `bottom_boundary_config_t` + adapter. Must close.
- **Layout: flat** (matches ADR 0034 heat precedent).
- **Hazards: 8 items.** Critical: Phase 0 gaps + tillage/calcgwl deferral decisions. Manageable: mini-sim writeback, multi-site qbot writes, boundtop-inside-headcalc semantics, swbotb=-2 mutation.
- **Suggested task count: ~16** across Phase 0/1/2.
- **Decomposition: monolithic arc.** Scope is comparable to heat; no need to split.

End of discovery.
