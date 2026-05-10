# Subsystem Migration Discovery: Soil-Water

**Date:** 2026-05-10
**Status:** superseded-by-split (kept as reference for the eventual soil-water-core arc)
**Migration #:** originally proposed as #5 — replaced by 4 coupling-surface arcs (#5 boundary, #6 crop water uptake, #7 atmosphere, #8 soil-water core)
**Branch:** `development` (no branch switch performed)

> **Superseded 2026-05-10.** This doc proposed a single ~100-global soil-water migration with internal phase split (SS-SW-A..E). User redirected to a **semantic split by coupling surface**: boundary first → crop water uptake → atmosphere → soil-water core. Each arc carves the soil fields touching its surface into `state%soilwater` incrementally. The Section 2 inventory and Section 7 hazards remain authoritative; the Phase decomposition is obsolete. The soil-water-core arc (the final pass) will lean on what's left after the three coupling-surface arcs land.

**Predecessor playbooks:**
- `docs/superpowers/specs/state-migration-playbook.md` (living document)
- ADRs 0030 (surfacewater pilot), 0031 (drainage), 0032 (solute), 0033 (cumulative cohorts), 0034 (heat — flat-layout precedent)

> Read-only discovery: no code changes made. All file:line references are anchors for the design phase.

---

## 1. Big picture

### Subsystem role

The soil-water subsystem is the **central physics core** of SWAP. It owns the Richards-equation solver, the per-step hydraulic state of the unsaturated zone (pressure heads, water contents, conductivities, fluxes), and every cumulative water-balance accumulator from soil-evaporation totals through bottom-boundary fluxes. Every other physics subsystem either feeds it (boundary conditions, root extraction, drainage, macropores, frost reduction of K) or reads from it (output, solute transport, oxygen-stress crop physiology, etc.).

It is the heaviest migration target yet — five times the owned-globals count of heat, with cumulative cohorts gated by multiple activity flags.

### Lines of code (home files)

```
   442  src/soil/soilgrid.f90
  1343  src/soil/soilhydraulics.f90
  5492  src/soil/sptabulated.f90              (mostly TSPACK — third-party tabulation library)
   810  src/soil/waterbalance.f90
   627  src/soil/WC_K_models_04_11.f90        (temperature-dependent K, iHWCKmodel 4-11)
   659  src/utils/soilhydraulicsutils.f90     (watcon/hconduc/dhconduc/moiscap/prhead/hcomean/dkmean)
 ─────
  9373  total
```

Excluding `sptabulated.f90` (TSPACK third-party tabulator), the soil home tree is ~3880 LoC of SWAP-specific code.

### Entry points called from `swap.f90`

In timestep call order (see `src/core/swap.f90`):

| # | Call site | swap.f90 line | Task arg | Phase | Currently takes `state`? |
|---|-----------|---------------|----------|-------|-------------------------|
| 1 | `CalcGrid()` | 182 | — | Init: grid geometry (dz, z, disnod, layer assignments, botcom, nod1lay) | NO |
| 2 | `SoilWater(1, state)` | 188 | 1 | Init: set up h, theta, k, kmean, dimoca, cofgen, volini, gwl from `swinco`/`swbotb` paths | YES (only state%heat%rfcp and state%surfacewater%vtair touched in current dual-write) |
| 3 | `BoundBottom(state)` | 303 | — | Per-step: bottom-boundary qbot setup (calls hconduc with state%heat%rfcp) | YES |
| 4 | `SoilWater(2, state)` | 317 | 2 | Per-step: reset intermediate/cumulative accumulators; call `SoilWaterStateVar(1)`; call `headcalc(state)` for Richards iteration | YES (passes through to headcalc) |
| 5 | `SoilWater(3, state)` | 333 | 3 | Per-step: refresh k/kmean, watstor, fluxes(state), macropore(4,state), integral(state), hysteresis | YES |
| 6 | `SoilWaterStateVar(2)` | 325 | 2 | Per-step (only if dt-reduce): restore h/theta/gwl/pond from hm1/thetm1/gwlm1/pondm1 | NO |
| 7 | `SoilWaterOutput(*, state)` | 226, 388, 401, 436, 437 | 1..4 | Output: defined in `swapoutput.f90:72` | YES |

Internal helpers (`headcalc`, `hysteresis`, `calcgwl`, `level`, `watertable`, `fluxes`, `integral`, `watstor`, `checkmassbal`, `ConvertDiscrVert`) are called from within the home tree or from `swapoutput.f90`. `fluxes(state)` and `integral(state)` already take `state` (SS-SWST Phase 2 Task 11).

### What `SoilWater(task=2)` does — sketch

```
SoilWater(2, state):
  if (flDayStart):  zero per-day reduction accumulators (iqredwet_day, iqreddry_day, …, qpotrot_day, qredtot_day)
  if (flzerointr):  zero intermediate set — inqrot, inq, iqrot, iqssdi, iqredwet/dry/sol/frs, ies0, iet0, iew0, iintc, iptra,
                    ipeva, ievap, iruno, irunoCN, iqbot, iqtdo, iqtup, irunon, iqdo(:), iqup(:), IPondBeg, IThetaBeg(:);
                    if (flMacroPore) macropore(5, state)
  if (flzerocumu):  zero cumulative set — cqssdi, cqrot, cqbot, cqbotdo, cqbotup, cptra, cpeva, cevap, cinund, crunon,
                    crunoff, crunoffCN, cqtdo, cqtup, cqprai; rebase volini=volact, pondini=pond;
                    if (flMacroPore) macropore(6, state)
  call SoilWaterStateVar(1)        ! save hm1=h, thetm1=theta, gwlm1=gwl, pondm1=pond
  call headcalc(state)             ! Newton-Raphson Richards solver — main physics
```

---

## 2. Owned globals (this subsystem writes)

This subsystem's owned globals fall into three reset-cadence groups (per playbook):

- **Instantaneous (I)** — overwritten/computed fresh each Richards step or each call to `SoilWater(1/3)`. Includes h, theta, q, gwl, pond, hatm, dimoca, k, kmean, FrArMtrx (init only), cofgen (init only), etc.
- **Intermediate (M)** — accumulated within an output period, reset on `flzerointr`. The `i*` prefix family.
- **Cumulative (C)** — accumulated across periods, reset on `flzerocumu`. The `c*` prefix family.

Gating: all soil-water cumulative fields are accumulated **unconditionally each timestep** in `integral(state)` — they are not flag-gated by drainage/macropore activity. (Macropore-specific cumulatives belong to the macropore subsystem, not soil-water; see Section 4.) This means soil-water needs **one** intermediate cohort and **one** cumulative cohort, not multiple partitioned cohorts. Contrast with surfacewater (ADR 0033 partitioning).

### Per-step / instantaneous fields

| Variable | Type / shape | Cadence | variables.f90 line | Description | Write sites (file:line) |
|---|---|---|---|---|---|
| `h(macp)` | real(8) | I | 803 | Pressure head per node | `soilhydraulics.f90:123,125,439,724,1036,1239`; `headcalc:382,435,753`; `hysteresis:1317,1333` |
| `hm1(macp)` | real(8) | I (save-state) | 812 | h at former time level | `soilhydraulics.f90:1227` (SoilWaterStateVar(1)) |
| `theta(macp)` | real(8) | I | 957 | Volumetric water content per node | `soilhydraulics.f90:112,237,443,519,715,1042,1240`; `hysteresis` (via cofgen update) |
| `thetm1(macp)` | real(8) | I (save-state) | 960 | theta at former time level | `soilhydraulics.f90:1228` |
| `thetar(macp)` | real(8) | I (init+hysteresis update) | 958 | Residual water content | `soilhydraulics.f90:949`; `hysteresis:1307,1311,1322` |
| `thetas(macp)` | real(8) | I (init+hysteresis update) | 959 | Saturated water content | `soilhydraulics.f90:950`; `hysteresis:1306,1323,1327,1328` |
| `thetsl(maho)` | real(8) | I (init only) | 961 | Saturated water content per layer | `soilhydraulics.f90:914,940,942` |
| `indeks(macp)` | integer | I (init+hysteresis) | 656 | Hysteresis branch index (+1 wet, −1 dry) | `soilhydraulics.f90:954,958`; `hysteresis:1279,1300` |
| `cofgen(21,macp)` | real(8) | I (init+hysteresis+tillage) | 760 | MvG params per node (mutated by hysteresis and tillage) | `soilhydraulics.f90:891,907–937,949,950,959,1313–1331`; **co-write**: `tillage.f90:243` |
| `dimoca(macp)` | real(8) | I | 783 | Differential moisture capacity | `soilhydraulics.f90:1048,312`; `hysteresis:1338` |
| `k(macp+1)` | real(8) | I | 852 | Hydraulic conductivity per node | `soilhydraulics.f90:130,170,176,181,238,453,467,471,520,1051,1171,1176`; capacities clamped to 1e-10 at flcaprise/upward-gradient nodes |
| `kmean(macp+1)` | real(8) | I | 855 | Mean K at interface | `soilhydraulics.f90:113,133,186,189,241,262,456,474,475,523,548,1054,1056,1174,1177,1242`; **co-write**: `boundtop.f90:143,161,168`; `boundbottom.f90:164,166` |
| `q(macp+1)` | real(8) | I | 877 | Inter-compartment water flux | `soilhydraulics.f90:880` (init reset); `waterbalance.f90:fluxes:341,344` |
| `qbot` | real(8) | I | 878 | Bottom-boundary flux | `soilhydraulics.f90:122,252,254,257,266,274,534,537,541,552,560,721`; `waterbalance.f90:fluxes:336`; **co-write**: `boundbottom.f90:104,107,143,152,155`; `frozencond.f90:216,243,287`; `swapoutput.f90:3818` (mini-sim writeback) |
| `qtop` | real(8) | I | 921 | Top-boundary flux | `soilhydraulics.f90:878,637`; **co-write**: `boundtop.f90:165` |
| `qrot(macp)` | real(8) | I | 894 | Root extraction sink per node | **External owner**: `rootextraction.f90:59,88,199,278,625,634,714,720,725,728,736,742,747,750` — soil-water *reads* qrot but does not write it. Listed here because traditional accounting put it in soil-water's globals. **Strongly consider transferring ownership to crop's root-state.** |
| `qrosum` | real(8) | I | 889 | Total root extraction | Same: written by rootextraction.f90; read by soilhydraulics fluxes (`waterbalance.f90:336`) |
| `qimmob(macp)` | real(8) | I | 887 | Mobile/immobile flux (fingered flow) | Written nowhere in home tree of soil-water; presumably macropore-owned. Read at `waterbalance.f90:344` |
| `qssdi(macp)` | real(8) | I | 900 | Subsurface drip irrigation flux | Owned by `irrigation.f90` (SSDI). Read at `soilhydraulics.f90:98,345`. Not owned by soil-water. |
| `qssdisum` | real(8) | I | 888 | Total SSDI flux | Read at `waterbalance.f90:336`. Not soil-water owned. |
| `hatm` | real(8) | I | 806 | Air pressure head near surface | `soilhydraulics.f90:869` (init=−2.75e5) |
| `gwl` | real(8) | I | 797 | Groundwater level | `waterbalance.f90:calcgwl:51,68,79,83,95,97,100`; `soilhydraulics.f90:757,1009,1019,1028,1243` |
| `gwlm1` | real(8) | I (save-state) | 801 | gwl at former time level | `soilhydraulics.f90:1230` |
| `nodgwl` | integer | I | 663 | Node directly above gwl | `waterbalance.f90:calcgwl:54,70,78,82,102` |
| `nodgwlflcpzo` | integer | I | (~1192 region) | Node above capillary-zone gwl | `waterbalance.f90:calcgwl:62,82,104` |
| `gwlflcpzo` | real(8) | I | 1193 | Capillary-zone gwl | `waterbalance.f90:calcgwl:63,85,105` |
| `bpegwl` | integer | I | 649 | Node at bottom of perched gwl | `waterbalance.f90:calcgwl:122,156` |
| `npegwl` | integer | I | 667 | Node directly above perched gwl | `waterbalance.f90:calcgwl:134,153,157` |
| `pegwl` | real(8) | I | 871 | Perched gwl | `waterbalance.f90:calcgwl:52,135,146,148,151` |
| `pond` | real(8) | I | 872 | Surface ponding depth | `soilhydraulics.f90:1031,1033,1067,1244,758`; `headcalc:758`; **co-write**: `boundtop.f90:144,163,181,254,262,270,286,304`; `tillage.f90:301,327`; `swapoutput.f90:3820` |
| `pondm1` | real(8) | I (save-state) | 874 | Ponding at former time level | `soilhydraulics.f90:1231,758` |
| `pondini` | real(8) | I (rebase on flzerocumu) | — | Initial ponding for the cumu period | `soilhydraulics.f90:1067,1156` |
| `volact` | real(8) | I | 963 | Current soil-profile storage | `waterbalance.f90:watstor:803,805` |
| `volm1` | real(8) | I (save-state) | 965 | Storage at former time level | `waterbalance.f90:watstor:802` |
| `volini` | real(8) | I (rebase on flzerocumu) | 964 | Storage at start of cumu period | `soilhydraulics.f90:1065,1155` |
| `wbalance` | real(8) | I (recomputed each integral call) | 968 | Cumulative water balance error | `waterbalance.f90:integral:538,541` |
| `FrArMtrx(macp)` | real(8) | I (init reset; macropore overwrites each step) | 1192 | Matrix horizontal-area fraction | `soilhydraulics.f90:1050,1061`; **co-write**: `macropore.f90:435,550` |
| `fluseksatexm(macp)` | logical | I (init only) | 634 | Per-node flag: use Ksatexm | `soilhydraulics.f90:925` |
| `numnod`, `numlay` | integer | I (init only) | 673, 672 | Grid dimensions | `soilgrid.f90:60,82` |
| `botcom(maho)` | integer | I (init only) | 650 | Bottom compartment per layer | `soilgrid.f90:78,83` |
| `nod1lay(maho)` | integer | I (init only) | 664 | First node of each layer | `soilgrid.f90:99` |
| `dz(macp)`, `z(macp)`, `disnod(macp+1)` | real(8) | I (init only) | 787,970,784 | Compartment thickness / node depth / inter-node distance | `soilgrid.f90:48,50,51,54,55,61` |
| `ztopcp(macp)`, `zbotcp(macp)` | real(8) | I (init only) | (1192 region) | Top/bottom depth per compartment | `soilgrid.f90:66,67,69,70` |
| `inpola(macp)`, `inpolb(macp)` | real(8) | I (init only) | — | Interpolation coefficients | `soilgrid.f90:86,88,89,91` |
| `numtab(macp)` | integer | I (init only) | 675 | Tabulated entries per compartment | `soilhydraulics.f90:895` |
| `ientrytab(macp,0:matabentries)` | integer | I (init only) | 724 | Tabulated indices per compartment | `soilhydraulics.f90:897` |
| `sptab(7,macp,matab)` | real(8) | I (init only) | 951 | Tabulated soil-physics per compartment | `soilhydraulics.f90:903` |
| `numnodNew, dzNew, …, ThetaNew, hNew, …` | mixed | I (init/output mini-sim only) | — | Working buffers for ConvertDiscrVert | `soilgrid.f90:ConvertDiscrVert:226,231–258,*` |
| `flwarn_hc`, `iwarn_hc`, `nstep_hc` | logical/integer | I (init+state-save) | 629–631 | Headcalc warn flags / step counter | `soilhydraulics.f90:78,79,699,703` |
| `fldecdt` | logical | I (per-step decision) | — | Set on Richards non-convergence | `soilhydraulics.f90:headcalc:761` |
| `fllowgwl` | logical | I (init+per-step) | 983 | gwl-below-profile flag | `soilhydraulics.f90:106,154` |
| `ldwet, spev, saev` | real(8) | I (init only under swinco≠3) | 864, 950, 945 | Boesten/Stroosnijder evap parameters | `soilhydraulics.f90:873–875` (init reset) |
| `runon` | real(8) | I (init reset) | 942 | Runon flux | `soilhydraulics.f90:877` |
| `evp(macp)` | real(8) | I (init reset) | — | Evaporation per node (currently always 0) | `soilhydraulics.f90:883` |
| `nraidt`, `nird` | real(8) | I (init reset) | 174, 255 | Net rain / net irrigation rate | `soilhydraulics.f90:870,871` |
| `cQMpLatSs` | real(8) | C (init+integral) | 1174 | Cumulative macropore lateral supply | `soilhydraulics.f90:888` (init reset). Could belong to macropore — verify ownership. |
| `tra` | real(8) | I (per-day reset) | 193 | Daily actual transpiration | `waterbalance.f90:integral:409,410` |

### Intermediate fields (`flzerointr` reset)

| Variable | Type | Cadence | variables.f90 line | Description | Write sites |
|---|---|---|---|---|---|
| `inq(macp+1)` | real(8) | M | 826 | Inter-compartment cumulative flux this period | `soilhydraulics.f90:1099,1101`; `waterbalance.f90:integral:342,352` |
| `inqrot(macp)` | real(8) | M | 831 | Root extraction this period | `soilhydraulics.f90:1097`; `integral:415` |
| `inqssdi(macp)` | real(8) | M | 832 | SSDI per node this period | `soilhydraulics.f90:1098`; `integral:420` |
| `iqrot` | real(8) | M | 842 | Total root extraction this period | `soilhydraulics.f90:1102`; `integral:413` |
| `iqssdi` | real(8) | M | 843 | Total SSDI this period | `soilhydraulics.f90:1103`; `integral:421` |
| `iqredwet, iqreddry, iqredsol, iqredfrs` | real(8) | M | 844–847 | Reduced extraction by wet/dry/salt/frost stress | `soilhydraulics.f90:1104–1107`; `integral:423–426` |
| `ies0, iet0, iew0` | real(8) | M | 819–821 | Reference Eppppppppppp / T / Ewet over period | `soilhydraulics.f90:1108–1110`; `integral:432–434` |
| `iintc` | real(8) | M | 822 | Interception this period | `soilhydraulics.f90:1111`; `integral:456` |
| `iptra, ipeva, ievap` | real(8) | M | 170,169,167 | Potential transp / pot evap / actual evap | `soilhydraulics.f90:1112–1114`; `integral:458–460` |
| `iruno, irunoCN, irunon` | real(8) | M | 848,940,849 | Runoff (Manning), runoff (CN), runon | `soilhydraulics.f90:1115,1116,1120`; `integral:461,495,462` |
| `iqbot` | real(8) | M | 837 | Bottom flux this period | `soilhydraulics.f90:1117`; `integral:468` |
| `iqtdo, iqtup` | real(8) | M | 838 | Top flux down / up | `soilhydraulics.f90:1118,1119`; `integral:470,472` |
| `iqdo(macp+1), iqup(macp+1)` | real(8) | M | 839 | Per-node down / up flux | `soilhydraulics.f90:1121,1122`; `integral:476,478` |
| `IPondBeg` | real(8) | M | — | Pond at start of period | `soilhydraulics.f90:1124` |
| `IThetaBeg(macp)` | real(8) | M | — | Theta at start of period | `soilhydraulics.f90:1126` |
| `igrai, inrai, iprec, igird, inird` | real(8) | M | 818,168,836,250,251 | Gross/net rain, total precip, gross/net irrigation | `integral:382–386` (set under `flzerointr`); accumulated in `integral:463–467` |
| `iqredwet_day, iqreddry_day, iqredsol_day, iqredfrs_day, iptra_day, qpotrot_day(:), qredtot_day(:)` | real(8) | per-day (flDayStart) | 13–17,833,834 | Per-day reduction accumulators | `soilhydraulics.f90:1083–1093`; `integral:416,417,427–431` |

The per-day group resets on `flDayStart`, not `flzerointr` — a **third cadence** distinct from intermediate. Worth a dedicated `per_day` cohort or a separate flag-gated reset.

### Cumulative fields (`flzerocumu` reset)

| Variable | Type | Cadence | variables.f90 line | Description | Write sites |
|---|---|---|---|---|---|
| `cqssdi` | real(8) | C | 774 | Cumulative SSDI | `soilhydraulics.f90:1135`; `integral:483` |
| `cqrot` | real(8) | C | 775 | Cumulative root extraction | `soilhydraulics.f90:1136`; `integral:484` |
| `cqbot` | real(8) | C | 765 | Cumulative bottom flux (net) | `soilhydraulics.f90:1137`; `integral:511` |
| `cqbotdo, cqbotup` | real(8) | C | 766,767 | Cumulative bottom down / up | `soilhydraulics.f90:1138,1139`; `integral:507,509` |
| `cptra, cpeva, cevap` | real(8) | C | 150,149,146 | Cumulative potential T / pot E / actual E | `soilhydraulics.f90:1140–1142`; `integral:487–489` |
| `cinund` | real(8) | C | 757 | Cumulative inundation | `soilhydraulics.f90:1143`; `integral:491` |
| `crunon, crunoff, crunoffCN` | real(8) | C | 781,780,940 | Cumulative runon / runoff (Manning) / runoff (CN) | `soilhydraulics.f90:1144–1146`; `integral:528,493,496` |
| `cqtdo, cqtup` | real(8) | C | 776,777 | Cumulative top down / up | `soilhydraulics.f90:1147,1148`; `integral:530,532` |
| `cqprai` | real(8) | C | 773 | Cumulative net rain on pond | `soilhydraulics.f90:1149`; `integral:527` |
| `caintc, cgrai, cnrai, cgird, cnird` | real(8) | C | 145,147,148,756,758 | Cumulative interception / gross & net rain / gross & net irrigation | `integral:498–504` (also reset by `meteoday.f90:420–422` and `irrigation.f90:93–95` — **multi-owner reset**) |
| `volini` | real(8) | C (rebase) | 964 | Storage at start of period (`volini = volact` after reset) | `soilhydraulics.f90:1155` |
| `pondini` | real(8) | C (rebase) | — | Pond at start of period (`pondini = pond` after reset) | `soilhydraulics.f90:1156` |
| `cQMpLatSs` | real(8) | C (init only) | 1174 | Cumulative macropore lateral supply | `soilhydraulics.f90:888` (init reset). Probably macropore-owned — flag for design phase. |

### Summary counts

- **Instantaneous (I)**: ~50 fields (per-node arrays, scalars, grid arrays, save-state copies, working buffers).
- **Per-day (D, reset on flDayStart)**: 7 fields (`iqred*_day`, `iptra_day`, `qpotrot_day`, `qredtot_day`).
- **Intermediate (M, reset on flzerointr)**: ~25 fields.
- **Cumulative (C, reset on flzerocumu, all gated by the same default-on condition)**: ~20 fields.

**Total owned globals: ~100** (rough; depends on how strictly we class grid/save-state/working-buffer fields).

Per-day fields are a small cohort distinct from intermediate cohort. The intermediate and cumulative cohorts can each be a flat sub-record (no further partitioning needed — all members share the same reset gate). ADR 0033 partitioning rules do not apply.

---

## 3. External readers — 5-category inventory

For each owned global, files outside `src/soil/`, `src/utils/soilhydraulicsutils.f90`, `src/state/`, `src/core/variables.f90`, and `src/core/initialize.f90` that consume soil-water owned data. Grep templates per playbook.

Distinct external reader files found via `use variables, only: …` import:

| # | File | Categories | Soil-water fields read |
|---|------|------------|-------------------------|
| 1 | `src/io/swapoutput.f90` | 1 (Output) | h, theta, q, qrot, qbot, gwl, pond, volini, volact, cofgen, dimoca, kmean, thetar, thetas, FrArMtrx, hm1, theta, cgrai/cnrai/cqrot/cqbot/etc — almost everything |
| 2 | `src/io/swap_csv_output.f90` | 1 (Output) | iqbot, iqrot, ievap, ipeva, irunon, iruno, iintc, igsnow, igird, gwl, pond, theta, h, etc. |
| 3 | `src/io/macroporeoutput.f90` | 1 (Output) | FrArMtrx (via macropore-side reads), inq, dz, gwl |
| 4 | `src/drainage/drainage.f90` | 2 (Compute) | gwl, pond, theta, h (for drainage-flux computation) |
| 5 | `src/drainage/surfacewater.f90` | 2 (Compute) | gwl, pond, theta (manage surface-water reservoir state) |
| 6 | `src/drainage/divdra.f90` | 2 (Compute) | gwl, h, theta (distribute drainage over depth) |
| 7 | `src/atmosphere/meteoday.f90` | 2 (Compute) | theta, ThetaRef, caintc, cgrai, cnrai, igrai, inrai, iprec (curve-number runoff path); **also a CO-WRITER** of caintc/cgrai/cnrai (intermediate/cumulative reset) — see Section 4. |
| 8 | `src/atmosphere/et.f90` | 2 (Compute) | theta, FrArMtrx, gwl (potential E/T inputs) |
| 9 | `src/boundary/boundtop.f90` | 2 + Co-write | h, pond, hatm, qtop, kmean, k(1), pondm1, q0; **WRITES** pond, kmean(1), qtop, reva |
| 10 | `src/boundary/boundbottom.f90` | 2 + Co-write | hbot, gwl, thetabot, kmean(numnod+1); **WRITES** qbot, kmean(numnod+1) |
| 11 | `src/heat/temperature.f90` | 2 (Compute) | theta, thetm1, thetas (heat capacity / conductivity inputs) |
| 12 | `src/heat/frozencond.f90` | 2 + Co-write | theta, thetas, qbot_nonfrozen, gwl; **WRITES** qbot |
| 13 | `src/crop/rootextraction.f90` | 2 + Co-write | theta, h, thetar, thetas, hm1, k, kmean, hwet; **WRITES** qrot, qrosum, h (clamp) |
| 14 | `src/crop/tillage.f90` | 2 + Co-write | h, theta, dz, layer; **WRITES** cofgen (Change_MvGpars), ParamVG, theta, h, pond (Adapt_WC_H) |
| 15 | `src/crop/cropgrowth.f90` | 2 (Compute) | theta, h (oxygen stress / phenology) |
| 16 | `src/crop/oxygenstress.f90` | 2 (Compute) | theta, thetas, h |
| 17 | `src/crop/management_soil.f90` | 2 (Compute) | theta, h, dz |
| 18 | `src/crop/irrigation.f90` | 2 + Co-write (cumu reset) | h, theta; **WRITES** cgird, cnird (cumu reset path — same multi-owner pattern as cgrai/cnrai) |
| 19 | `src/solute/solute.f90` | 2 (Compute) | theta, h, q, inq, gwl, FrArMtrx |
| 20 | `src/solute/agetracer.f90` | 2 (Compute) | theta, q, h |
| 21 | `src/macropore/macropore.f90` | 2 + Co-write | h, theta, cofgen, dz, FrArMtrx, layer, kmean; **WRITES** FrArMtrx (lines 435, 550) |
| 22 | `src/macropore/macrorate.f90` | 2 (Compute) | theta, h |
| 23 | `src/utils/surfacewaterutils.f90` | 2 (Compute) | pond, rsro, rsroexp (qhtab path) |
| 24 | `src/core/swap.f90` | 5 (Call-site arg) | passes `state%heat%tsoil` to CropGrowth; uses `theta` only in DLL-mode comment block at line 77 (commented out) |
| 25 | `src/io/toml/config_to_variables.f90` | 4 (Init-seed) | seeds h, theta (from h_file / initial profile), pond (from soil.initial.pond) |

**Distinct external reader files:** 25 (vs heat's 15).

**Per-field heat-map (selected hot fields):**

- `theta`: read by ~12 external files (3 output + 9 compute). Heaviest reader: drainage, divdra, solute, oxygenstress, et, meteoday, rootextraction, snow (via swapoutput), tillage, macropore, agetracer.
- `h`: similar to theta — ~10+ external readers. Compute hot path: rootextraction (matric-flux potential), tillage (Adapt_WC_H), oxygenstress, macropore.
- `gwl`: ~8 readers — drainage (heavily), output, FrozenBounds (already on state).
- `pond`: ~6 readers — boundtop (co-writer), surface_water, boundbottom, output, tillage.
- `cofgen`: ~3 readers + 2 writers — soilhydraulics (init/hysteresis), tillage (Change_MvGpars), macropore (Adapt_WC_H reads, via prhead), `ConvertDiscrVert` (rebuilds cofgenNew).
- `FrArMtrx`: 1 internal writer (soilhydraulics init), 2 external writers (macropore lines 435, 550), and many readers in waterbalance/divdra.
- `qbot`: 3 writers outside soil-water — frozencond (3 sites), boundbottom (5 sites), swapoutput (mini-sim writeback). Pattern matches drainage's pre-migration multi-write.

---

## 4. Co-writers — multi-subsystem writes (Section 4)

Files outside the soil home tree that **WRITE** soil-water owned globals. Each is a Phase 2 (or Phase 1 dual-write) candidate.

| Co-writer file | Field(s) written | Write context | Has `state` in scope today? |
|---|---|---|---|
| `src/heat/frozencond.f90:216,243,287` | `qbot` | `FrozenBounds(state)` — already plumbed with state%drainage; modifies `qbot` global directly when frozen | YES (`state%drainage`, `state%surfacewater`); does **not** read/write state's soil-water fields yet |
| `src/boundary/boundbottom.f90:104,107,143,152,155` | `qbot` | Per-step bottom-boundary: sets qbot per `swbotb` mode (sine, tab, exp, qhtab) | YES (passed `state` for rfcp; current signature `BoundBottom(state)`) |
| `src/boundary/boundbottom.f90:164,166` | `kmean(numnod+1)` | Sets bottom-interface K when prescribed-head bottom | YES (state) |
| `src/boundary/boundtop.f90:143,161,163,165,168,181,254,262,270,286,304` | `pond`, `kmean(1)`, `qtop`, `reva` | Per-step top-boundary: sets ponding height + top-interface K + qtop + reva based on `Z_Tp` / runoff mode | YES (state passed in SS-HEAT Task 9 already) |
| `src/macropore/macropore.f90:435,550` | `FrArMtrx(macp)` | When macropores active, replaces matrix-area fraction by `1 - VlMpStCp/Dz`; resets to 1.0 when off | YES (`state` is already plumbed to macropore tasks 1–6) |
| `src/crop/rootextraction.f90:59,88,108,199,278,282,515,615–636,666,714–756` | `qrot(:)`, `qrosum`, `h(node)` (clamp at L382) | Computes root-extraction sink terms. The full qrot/qrosum lifecycle lives here — soil-water *only reads* these. **Recommend transferring ownership to crop's state.** | YES (rootextraction.f90 already takes `state` for MatricFlux) |
| `src/crop/tillage.f90:209–248, 263–301, 332, 344` | `cofgen`, `ParamVG`, `theta`, `h`, `pond` | `Change_MvGpars` rewrites cofgen/ParamVG; `Adapt_WC_H` rewrites theta/h/pond after redistribution. Called at init (task 1), per-day (task 2), per-output (task 3). | NO — currently uses bare `use variables`. Needs state plumbing if theta/h/cofgen migrate. |
| `src/atmosphere/meteoday.f90:420–422` | `cgrai`, `cnrai`, `caintc` | Resets cumulative atmospheric accumulators on `flzerocumu` independently of `integral()`. **Multi-owner reset pattern** — same field reset from two subsystems. | YES (state arg available; not currently used for these) |
| `src/crop/irrigation.f90:93–95` | `cgird`, `cnird` | Resets cumulative irrigation cumus on `flzerocumu`. Same multi-owner pattern. | YES |
| `src/io/swapoutput.f90:3818–3823` | `qbot`, `gwl`, `pond`, `theta(nod)`, `h(nod)` | Output mini-simulation: snapshot/restore pattern (saves qbot/gwl/pond/theta/h, runs mini sim, restores). | YES (state passed) |

**Total co-writers: 10 distinct files.** Highest impact: tillage (cofgen/theta/h writeback path) and rootextraction (qrot ownership). Boundary co-writes (qbot/qtop/pond/kmean) are already cleanly inside the timestep iteration and have state plumbing.

---

## 5. Plumbing-state status

Per playbook Section "Argument threading": which subroutines already have `state` as a dummy arg, and which need plumbing?

### Already plumbed (✓)

| Subroutine | Signature | Source |
|---|---|---|
| `soilwater(task, state)` | `type(swap_state_t), intent(inout) :: state` | `soilhydraulics.f90:842,855` |
| `headcalc(state)` | `intent(inout)` | `soilhydraulics.f90:25,42` |
| `fluxes(state)` | `intent(in)` | `waterbalance.f90:323,329` |
| `integral(state)` | `intent(inout)` | `waterbalance.f90:373,377` |
| `ConvertDiscrVert(part,…,state)` | `intent(in)` | `soilgrid.f90:182,194` |

These were plumbed in earlier migration arcs (mostly for drainage/surfacewater drain/qdra access, and heat/rfcp).

### NOT yet plumbed (need state arg in Phase 1)

| Subroutine | Source | Why state is needed in soil-water migration |
|---|---|---|
| `SoilWaterStateVar(task)` | `soilhydraulics.f90:1213` | Writes hm1/thetm1/gwlm1/pondm1; if those migrate, state must be in scope |
| `hysteresis()` | `soilhydraulics.f90:1262` | Reads h/hm1/theta, writes thetar/thetas/cofgen/dimoca — all owned |
| `calcgwl()` | `waterbalance.f90:39` | Writes gwl, nodgwl, pegwl, bpegwl, npegwl, nodgwlflcpzo, gwlflcpzo; reads h, pond, theta |
| `level(swoptlev,node,nodheq1)` | `waterbalance.f90:199` | Reads h, z (function — could pass values explicitly) |
| `watertable(node,…)` | `waterbalance.f90:267` | Reads ThetaS, Theta, h |
| `watstor()` | `waterbalance.f90:795` | Writes volact, volm1; reads theta, dz, FrArMtrx |
| `checkmassbal(flopenfiledev,…)` | `waterbalance.f90:583` | Already takes many array args; reads pond/Ssnow/igrai/etc |
| `CalcGrid()` | `soilgrid.f90:21` | Writes numnod, numlay, dz, z, disnod, layer, botcom, nod1lay (grid) — these are likely retained as legacy globals or moved to a non-state config (see open question) |

### External callers needing plumbing

Many external files import soil-water globals via `use variables`. Phase 2 will need to convert those to `state%soilwater%X` reads (or pass scalars as args). The full list is the 25-file external reader inventory (Section 3).

**The biggest threading consequence:** `state` must reach inside `rootextraction.f90`, `boundbottom.f90`, `boundtop.f90`, `macropore.f90`, `tillage.f90`, all `swapoutput.f90` output blocks, `solute.f90`, `cropgrowth.f90`, `oxygenstress.f90`, `et.f90`, `meteoday.f90`. Most already have it for other reasons; the residual is `tillage.f90`, `meteoday.f90` (state passed but unused), `cropgrowth.f90` (state already passed for tsoil), `oxygenstress.f90`, `et.f90`.

---

## 6. Config / Phase 0 candidates

Soil and bottom-boundary configs already exist (`src/config/soil_config.f90`, `src/config/bottom_boundary_config.f90`) and are populated by `read_soil_toml.f90` and `read_bottom_boundary_toml.f90`. The legacy `.swp` reader populated these globals; today the TOML adapter mirrors them.

### Existing typed config coverage (soil_config_t)

- Grid: `swdiscrvert`, `numnodnew`, `dznew(:)`.
- Per-layer hydraulic params: `ores`, `osat`, `alfa`, `npar`, `ksatfit`, `lexp`, `alfaw`, `h_enpr`, `ksatexm`, `bdens`.
- Frost: `swfrost`, `tfroststa`, `tfrostend`, `swsublim`.
- Initial state: `ssnow`, `slw`, `pond`, `ldwet`, `dt`, `atmin7(7)`, `h_file`, `tsoil_file`, `cml_file`.
- Tillage event: `date`, `z`, `intensity`, `type_id`, `id`, `rho_cons`, `rho_tillage`, `k_R`, `rho_match`, `N_match`.
- Switches: `swsophy`, `swhyst`, `swinco`, `swmacro`, `swscal`, `swtill`.
- Profile/runon: `gwli`, `pondini`, `pondmx`, `ksatexm`, `rsoil`, `rsro`, `rsroexp`, `swrunon`.
- Tabulated model selectors: `i_n_model`, `iRedist`.

### Existing typed config coverage (bottom_boundary_config_t)

- `swbotb`, `sw2`, `sw3`, `sw4`, `swqhbot`.
- Files: `gwl_file`, `qbot2_file`, `haquif_file`, `qbot4_file`, `qhbot_file`, `hbot5_file`.
- Tables: `swc_table`, `qbot_table`, `cofqha_table`.
- Cauchy: `shape`, `hdrain`, `rimlay`, `swbotb3impl`, `aqave`, `aqamp`, `aqper`, `aqtmax`.
- Other: `hbot`, `rhobot`.

### Potential Phase 0 gaps

Compared with what `soilhydraulics.f90:headcalc` and `boundbottom.f90` actually read, the typed config looks reasonably complete. Possible gaps to investigate during the design phase (not exhaustive):

1. **`deepgw`** (variables.f90:782) — hydraulic head in aquifer (used by Cauchy/swbotb=3). Need to check if it is populated from `bottom_boundary%shape/aqave` or is a separately-set global. Quick grep: only one assignment site outside `config_to_variables.f90` — needs verification it's covered.
2. **`hplate`** (variables.f90:813) — lysimeter plate pressure (used by swbotb=8 path in headcalc). Not visible in bottom_boundary_config.f90 — probable Phase 0 candidate.
3. **`sinmax`, `sinamp`, `sinave`** (variables.f90:947–949) — bottom-flux sine parameters (used by swbotb=2 sine option). Need to verify they are in `swc_table` / `qbot_table` decomposition or are separate config fields.
4. **`cofqha`, `cofqhb`** (variables.f90:761–762) — qh exponential coefficients (used by swbotb=4 + swqhbot=1). Likely from `cofqha_table` but a quick grep would confirm.
5. **`gwlconv`, `MaxBackTr`, `MaxIt`, `CritDevh1Cp`, `CritDevh2Cp`, `CritDevPondDt`, `CritDevBalCp`, `CritDevBalTot`, `Critdz`, `CritDevMasBal`** — Richards-solver convergence parameters. Some are `parameter` constants in `headcalc`; the global ones (gwlconv, CritDevMasBal) need confirmation of TOML adapter coverage. These are Phase-0 candidates if not already covered.
6. **`hsurf, ftoph, dtmin, fldumpconvcrit`** — solver runtime controls. Generally set elsewhere; verify.

This Phase 0 set is much smaller than heat's (which had 4 missing analytical-method fields plus 2 tables). The soil/bottom-boundary configs are mature.

---

## 7. Known coupling hazards

Pre-migration warnings inherited from earlier arcs, plus new ones surfaced during this discovery.

### Hazard #1 — `hconduc` `tsoil_loc = 0.0_real64` sentinel (heat-arc residual)

**Scope: medium.** Inherited from SS-HEAT Task 9. The `hconduc` function in `src/utils/soilhydraulicsutils.f90:435` falls through to a `tsoil_loc = 0.0_real64` sentinel when `tsoil_node` is not passed. **No caller passes `tsoil_node` today** (14+ callsites in soilhydraulics, boundtop, boundbottom, rootextraction, macropore, swapoutput, tillage). The iHWCKmodel 4-11 path is therefore unreachable with anything other than `tsoil=0.0` (an incorrect physics value), and the design note records it as "unreachable in regression." The soil-water migration is the natural place to either (a) thread `state%heat%tsoil(node)` through hconduc call chains, or (b) document the iHWCKmodel 4-11 path as deprecated/unsupported.

**Recommendation:** small Phase 1 task — pass `state%heat%tsoil(node)` from the ~10 callers inside the soil home tree (where state%heat is already accessible). The remaining 4-5 callers in boundary/rootextraction/macropore/tillage need their own state plumbing; can defer.

### Hazard #2 — Multi-owner cumulative reset (cgrai, cnrai, caintc, cgird, cnird)

**Scope: medium.** The cumulative fields `cgrai`, `cnrai`, `caintc` are accumulated by `integral()` (soil-water) but reset by `meteoday.f90:420–422` on `flzerocumu`. Similarly `cgird`, `cnird` are reset by `irrigation.f90:93–95`. This is a different cross-subsystem pattern than ADR 0033 partitioning — here the cumu fields belong to soil-water's `integral()` accumulation path but the reset is owned by atmospheric / irrigation subsystems. Two options:

- (a) Move these specific cumus to the meteoday / irrigation owners' state and have soil-water `integral` accumulate via owner-side fields.
- (b) Keep them in soil-water cumulative cohort and have meteoday/irrigation call `state%soilwater%cumulative%reset_atmosphere_subset()`.

ADR 0033 Phase A finding makes (a) the cleaner choice — assign ownership to the accumulator's actual "home" subsystem (meteorology owns rainfall cumus, irrigation owns irrigation cumus). The soil-water cohort then contains only the fields that soil-water accumulates and resets.

**Recommendation:** in design phase, decide whether `cgrai/cnrai/caintc/cgird/cnird` move to a meteo/irrigation state, or stay in soil-water with subset-reset comments. Lean toward moving them — it's the cleaner end state.

### Hazard #3 — `rfcp` co-write at `soilhydraulics.f90:886`

**Scope: small.** Already exists post-heat. `if (allocated(state%heat%rfcp)) state%heat%rfcp = 1.0d0` is the only soil-water-side writeback into heat state. It runs in `SoilWater(task=1)`. No change needed for soil-water migration — this is a heat-state field being initialized at the right place.

### Hazard #4 — `dead-code` rfcp reset loop (mentioned in heat ADR 0034)

**Scope: small.** Per the task brief, ADR 0034 left a "harmless dead code" comment on the `rfcp(i) = 1.0d0` loop under `swfrost=0`. Looking at line 886 (`if (allocated(state%heat%rfcp)) state%heat%rfcp = 1.0d0`), this is the **init-time** seed (task=1 only), not a per-step reset. There is no per-step `rfcp = 1.0d0` loop in soilhydraulics today — what existed in the legacy code was already collapsed by SS-HEAT. **No action required for soil-water migration.**

### Hazard #5 — `macropore` co-writes `FrArMtrx`

**Scope: small.** `macropore.f90:435,550` writes `FrArMtrx(:)` directly each timestep. Soil-water `SoilWater(1)` also writes `FrArMtrx(i) = 1.d0` as a non-macropore fallback. If `FrArMtrx` migrates to `state%soilwater`, both writers need access. macropore.f90 already takes `state` for surfacewater/drainage — adding a `state%soilwater%FrArMtrx` write is mechanical.

### Hazard #6 — `tillage` co-writes `cofgen`, `theta`, `h`, `pond`

**Scope: large.** `tillage.f90:Change_MvGpars` rewrites the global `cofgen` array using `Bdens`/`Rho_last` ratios per layer. `Adapt_WC_H` then redistributes `theta`/`h` and adjusts `pond`. This runs:

- At init (`DoTillage(1)`) — called BEFORE `SoilWater(1, state)` in `swap.f90` ordering (line 184 vs 188). So tillage's initial cofgen modification happens before any soil-water state allocation.
- Per-day (`DoTillage(2)`) — called inside the daily loop.
- Per-output (`DoTillage(3)`) — called inside the output block.

`tillage.f90` uses bare `use variables` and has no state plumbing. The soil-water migration **must add `state` to `DoTillage(itask, state)`** before migrating cofgen/theta/h to state%soilwater. Same pattern as the heat-arc `Temperature(task, state)` plumbing extension.

**Recommendation:** treat tillage state plumbing as a Phase 1 task in the soil-water arc. Or split into a Phase 1.5 mini-arc.

### Hazard #7 — `rootextraction` owns qrot but soil-water reads it

**Scope: medium.** Today `qrot(macp)` and `qrosum` are declared in `variables.f90` (lines 894, 889) and visually grouped with soil-water fields, but **all writes** happen in `rootextraction.f90`. Soil-water reads them in `headcalc` (as sink terms) and `fluxes`/`integral` (as flux components). This is a classic ownership-misalignment hazard.

**Recommendation:** in design phase, decide whether to:
- (a) Migrate qrot/qrosum to a future `state%crop%root%qrot` (defer to crop arc).
- (b) Migrate qrot/qrosum to `state%soilwater%qrot` (soil-water arc; rootextraction continues to write directly).

Option (a) is cleaner but expands scope. Option (b) preserves the current data flow.

### Hazard #8 — `boundtop` and `boundbottom` co-write multiple soil-water fields

**Scope: medium.** Top boundary writes `pond`, `kmean(1)`, `qtop`, `reva`. Bottom boundary writes `qbot`, `kmean(numnod+1)`. Both already take `state`. The mechanical change is `pond = …` → `state%soilwater%pond = …`. The hazard is the number of distinct sites (~12 in boundtop, ~6 in boundbottom).

### Hazard #9 — `swapoutput.f90:3818–3823` mini-simulation writeback

**Scope: small.** Output runs an "om" mini-simulation that snapshots `qbot, gwl, pond, theta, h`, runs a parallel sim, and restores them from snapshot. Same pattern as Drainage's mini-sim writeback (existing). Once h/theta/qbot/gwl/pond migrate, the writeback must target `state%soilwater%X` instead of the legacy globals.

### Hazard #10 — `ConvertDiscrVert` working-buffer expansion

**Scope: small/medium.** `ConvertDiscrVert` (soilgrid.f90:178) takes a large argument list of "New" arrays (`thetaNew`, `hNew`, `inqNew`, etc.) and reads from globals (`theta, h, inq, IThetaBeg, cofgen, dz, FrArMtrx, numlay, botcom`). When these globals migrate to state, the routine signature stays the same (the "New" arrays are passed-in output buffers) but the global reads inside the body must rebind to `state%soilwater%X`. Already takes `state` for `state%surfacewater%intermediate%inqdra`. Mechanical change.

### Hazard #11 — Grid dimensions / arrays: legacy globals vs state

**Scope: medium.** `numnod`, `numlay`, `dz(:)`, `z(:)`, `disnod(:)`, `layer(:)`, `botcom(:)`, `nod1lay(:)`, `ztopcp/zbotcp`, `inpola/inpolb` are written once by `CalcGrid()` and read by ~everything. Per the task brief, these are typically retained as legacy globals or moved to a non-state config struct. The heat-discovery doc explicitly flagged the same issue (dz, z, disnod, layer are read everywhere).

**Recommendation:** keep grid dimensions as legacy globals **OR** introduce a small `grid_t` typed record in `swap_state_t` (parallel to subsystem states). The latter is cleaner long-term but expands scope. Defer the decision to the design phase; for the migration, keep them as legacy globals and revisit when refactoring the grid module is on the menu.

---

## 8. Reset-cadence analysis (cohort vs flat decision)

Reset-cadence groups for soil-water owned globals:

| Cadence | Count | Activity gate | Cohort target |
|---|---|---|---|
| Instantaneous (overwrite each step, no flag gate) | ~50 | none / `flMacroPore` for FrArMtrx | Flat fields on `soilwater_state_t` |
| Per-day (`flDayStart`) | 7 | none | `soilwater_per_day_t` cohort (small) |
| Intermediate (`flzerointr`) | ~25 | none | `soilwater_intermediate_t` cohort |
| Cumulative (`flzerocumu`) | ~20 | none | `soilwater_cumulative_t` cohort |

**All cumulatives share a single reset gate** (`flzerocumu`) — no partitioning by activity flag needed. Contrast with surfacewater's two cumulatives (drainage vs reservoir gates per ADR 0033). The cumulative cohort can be a single sub-record.

**Recommended shape:**

```fortran
type :: soilwater_state_t
   ! Instantaneous fields — direct on the state type
   real(real64), allocatable :: theta(:), h(:), q(:), k(:), kmean(:), dimoca(:), …
   real(real64) :: gwl, pond, hatm, volact, volm1, volini, pondini, wbalance, …
   integer :: nodgwl, bpegwl, npegwl, …

   ! Save-state copies (former time level)
   real(real64), allocatable :: hm1(:), thetm1(:)
   real(real64) :: gwlm1, pondm1

   ! Cohort sub-records
   type(soilwater_per_day_t)      :: per_day        ! reset on flDayStart
   type(soilwater_intermediate_t) :: intermediate   ! reset on flzerointr
   type(soilwater_cumulative_t)   :: cumulative     ! reset on flzerocumu
end type
```

This is the "mix of instantaneous + cumulative" shape from the task brief — closest analog is surfacewater pre-ADR-0033 (before partitioning).

**Key open question:** the per_day cohort is small (7 fields) — should it be (a) a dedicated cohort with `reset()`, (b) folded into intermediate (since per-day fluxes are a subset of intermediate-period fluxes), or (c) split out and reset inline without a cohort? Design phase should decide.

---

## 9. Scope estimate

| Metric | Soil-water | Heat (ADR 0034 — reference) | Surfacewater (pilot — reference) |
|---|---|---|---|
| Owned globals | ~100 (50 I + 7 D + 25 M + 20 C) | 12 (all instantaneous) | ~25 (mix I + M + 2 C cohorts) |
| External reader files | 25 | 15 | 12 |
| External read sites (rough) | 200+ | ~50 | ~80 |
| Co-writers (cross-subsystem writes) | 10 | 1 (`soilhydraulics` rfcp reset) | 3 |
| Phase 0 config candidates | 5 (hplate, sinmax/amp/ave, gwlconv, deepgw verification) | 6 (analytical-method + tables) | 2 |
| Subsystems entry-points needing state plumbing | 6 (SoilWaterStateVar, hysteresis, calcgwl, watstor, watertable, CalcGrid?) + tillage co-writer | 3 (Temperature, FrozenCond, plus FrozenBounds already plumbed) | 4 (SurfaceWater tasks) |
| Cohort design complexity | 3 cohorts (per_day, intermediate, cumulative) | 0 cohorts (flat) | 3 cohorts post-ADR 0033 |
| Recommended phase shape | Decompose | Single arc | Two phases |

**Recommended decomposition:** **split into multiple sub-arcs**. Given the scale (5–10× heat), the natural split boundaries are:

- **SS-SW-A (Phase 0 + Phase 1a):** State type design + threading. Add `soilwater_state_t` to `swap_state_t`. Plumb state through tillage. Phase 0 config gaps closed.
- **SS-SW-B (Phase 1b — instantaneous fields):** Migrate the ~50 instantaneous fields (h, theta, q, gwl, pond, …). Dual-write. Drop legacy globals once all readers cut over. Handle `hconduc` `tsoil_node` threading (Hazard #1).
- **SS-SW-C (Phase 1c — per-day + intermediate):** Migrate per_day and intermediate cohorts. Replace `flzerointr` block with `state%soilwater%intermediate%reset()`.
- **SS-SW-D (Phase 1d — cumulative):** Migrate cumulative cohort. Resolve multi-owner reset (Hazard #2) by relocating cgrai/cnrai/caintc/cgird/cnird to meteo/irrigation state, OR keep with subset reset.
- **SS-SW-E (Phase 2 — external reader migration):** Convert 25 external files to read from `state%soilwater`. This phase is itself decomposable along reader-file groupings (output, drainage, crop, atmosphere, boundary).

Five sub-arcs is plausible; the design phase should refine to match team velocity and check-full discipline.

---

## 10. Open questions for design phase

1. **Cohort vs flat shape.** Recommendation: **partitioned cohort** (1 flat + 3 cohort sub-records: per_day, intermediate, cumulative). Confirm whether per_day deserves its own cohort or folds into intermediate. Confirm cumulative is single-cohort (single reset gate).

2. **`hconduc` `tsoil_node` residual.** Fix in this arc (thread `state%heat%tsoil(node)` through hconduc call chain from soil-water home tree first) or leave for crop-arc? Recommendation: fix in this arc — the home-tree callers are easy.

3. **`qrot` / `qrosum` ownership.** Move to `state%soilwater` (option b) for arc-locality, or wait for crop arc and migrate to `state%crop%root` (option a)? Recommendation: option (b) — preserves data flow; revisit if crop arc develops.

4. **Multi-owner cumulative reset (cgrai/cnrai/caintc/cgird/cnird).** Move to meteo/irrigation state (option a) or keep with subset reset (option b)? Recommendation: option (a) — cleaner long-term, matches ADR 0033 ownership principle.

5. **Tillage state plumbing.** Do as part of SS-SW-A (recommended) or split into a pre-arc? Recommendation: include in SS-SW-A — tillage is small (~400 LoC) and the change is mechanical (add `state` arg + update 3 call sites).

6. **Grid dimensions (`numnod, numlay, dz, z, disnod, layer, botcom, nod1lay`).** Keep as legacy globals (recommended, matches heat-arc decision) or introduce a `grid_t` typed record. Recommendation: keep as legacy globals; flag as a future-arc topic.

7. **`FrArMtrx` ownership.** Soil-water-owned with macropore co-write (current pattern) or macropore-owned? Recommendation: soil-water-owned — `FrArMtrx` is fundamentally a soil-matrix property; macropore's write is a per-step override.

8. **`cQMpLatSs` ownership.** Listed under soil-water (variables.f90:1174) but obviously a macropore cumu. Recommendation: move to macropore state during this arc OR explicitly defer to a macropore arc.

9. **Phase 0 config gaps verification.** Confirm `hplate`, `sinmax/amp/ave`, `cofqha/b`, `gwlconv`, `deepgw` are properly TOML-adapter-covered before Phase 1 starts.

10. **`flzerointr`/`flzerocumu` flag ownership.** These are set by `TimeControl` — verify their reset gating doesn't interact weirdly with the new cohort `reset()` calls. Confirm reset-block ordering in `SoilWater(task=2)` versus `integral()`.

11. **`ConvertDiscrVert` legacy-global reads.** When `theta, h, inq, IThetaBeg, cofgen, dz, FrArMtrx, numlay, botcom` migrate, `ConvertDiscrVert` (already takes state) needs its body updated. Mechanical but touches a 250-line subroutine. Plan: address in SS-SW-B alongside the affected globals.

12. **`SoilWaterStateVar(2)` restore path.** When `dt`-reduce backs out, this routine restores `h, theta, gwl, pond, kmean(numnod+1)` from the m1-copies. Confirm the restore path matches the new state type's invariants (especially `kmean(numnod+1) = k(numnod)` at line 1242 — note this is a separate write of a state field outside the m1-copy set).

13. **`outheapar`-style hazards in soil-water.** Heat had an `outheapar` working-buffer hazard (output mutates state). For soil-water, the analog is `swapoutput.f90:3818–3823` mini-sim writeback. Confirm no other output sites mutate soil-water globals.

---

## Discovery summary

- **Subsystem size:** 5× heat (~100 owned globals, 25 reader files, 10 co-writers).
- **Reset cadence:** 4 groups (instantaneous + per_day + intermediate + cumulative). Cumulative is single-gate — no partitioning needed.
- **Recommended state-type layout:** flat fields + 3 cohort sub-records.
- **Recommended decomposition:** 5 sub-arcs (state-type+threading, instantaneous, per_day+intermediate, cumulative, external readers).
- **Highest-risk co-writers:** tillage (needs state plumbing), rootextraction (qrot ownership question), atmosphere/irrigation (multi-owner cumu reset).
- **Heat residual to clean up:** `hconduc` `tsoil_node` sentinel — small task in SS-SW-B.
- **Phase 0 gaps:** 5 candidates to verify; soil/bottom-boundary configs are mature.

End of discovery.
