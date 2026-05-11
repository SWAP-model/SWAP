## Subsystem Migration Discovery: Soil-Water Core (Migration #8 — FINAL coupling-surface arc)

**Date:** 2026-05-11
**Status:** discovery draft (read-only investigation, no code changes)
**Migration #:** 8 — the residual of the soil-water mega-arc after #5 boundary, #6 crop-uptake, #7 atmosphere
**Branch:** `development` (no branch switch)

> This is the FINAL of the four coupling-surface arcs that replaced the original
> single ~100-global mega-arc. The mega-discovery
> (`docs/superpowers/specs/2026-05-10-state-migration-soilwater-discovery.md`)
> Section 2 inventory and Section 7 hazards remain authoritative; the original
> Phase decomposition (SS-SW-A..E) is obsolete. After #5 (12 boundary fields),
> #6 (23 crop-uptake fields), and #7 (40 atmosphere fields) — 75 fields have
> carved out — this arc owns the **Richards interior** + the **cumulative /
> intermediate cohorts** + the deferred residuals (`pond`, `gwl`, `kmean`,
> hconduc `tsoil_node` sentinel).

**Predecessor docs read (and relied on):**
- `docs/superpowers/specs/2026-05-10-state-migration-soilwater-discovery.md` (mega-discovery, Sec 2 + Sec 7)
- `docs/adr/0035-state-migration-boundary.md` (12 fields)
- `docs/adr/0036-state-migration-crop-uptake.md` (23 fields)
- `docs/adr/0037-state-migration-atmosphere.md` (40 fields, two cohorts — `intr` + `cumu`)
- `docs/superpowers/specs/2026-05-10-state-migration-boundary-design.md`
- `docs/superpowers/specs/2026-05-10-state-migration-crop-uptake-design.md`
- `docs/superpowers/specs/2026-05-11-state-migration-atmosphere-design.md`
- `docs/superpowers/specs/state-migration-playbook.md` (20 lessons)
- `src/state/soilwater_state.f90` (current 35-field record — THIS ARC EXTENDS)
- `src/state/atmosphere_state.f90` (TWO-cohort precedent)

---

## 1. Big picture

### Subsystem role (still the heart of SWAP)

The soil-water core is the Richards-equation solver and the per-step hydraulic
interior of the unsaturated zone: pressure head `h`, water content `theta`,
inter-compartment flux `q`, hydraulic conductivity `k`/`kmean`, capacity
`dimoca`, soil-physics parameter table `cofgen`, plus the matrix-fraction array
`FrArMtrx`, the per-node "former time level" copies `hm1`/`thetm1`, the gwl
geometry helpers (`nodgwl`, `pegwl`, `bpegwl`, `npegwl`, `gwlflcpzo`,
`nodgwlflcpzo`), the storage scalars (`volact`, `volm1`, `volini`, `pondini`),
the air-pressure-head sentinel `hatm`, the deferred boundary fields (`pond`,
`gwl`), and the entire **intermediate + cumulative** cohort families
(`inq*`/`iq*`/`cq*`/`crun*`/`cinund`/`wbalance` etc.).

After the three coupling-surface arcs landed, the mega-discovery's ~100 globals
break down as:
- ✅ 12 migrated by boundary (ADR 0035) — top/bottom flux scalars
- ✅ 23 migrated by crop-uptake (ADR 0036) — qrot family + JvL fields
- ✅ 40 migrated by atmosphere (ADR 0037) — peva/ptra + intr/cumu cohorts
- ✅ 5 closed as Phase 0 typed-config — gwlconv, hplate, sinmax/amp/ave,
  cofqha/b (verified Sec 6 below)

The **residual is what this arc owns: ~60 fields** falling into 4 cadences
(flat instantaneous, flat per-day, intermediate cohort, cumulative cohort).

### LoC (home files in scope)

```
   442  src/soil/soilgrid.f90
  1344  src/soil/soilhydraulics.f90
   838  src/soil/waterbalance.f90
  5492  src/soil/sptabulated.f90              (TSPACK third-party tabulator)
   627  src/soil/WC_K_models_04_11.f90        (temperature-dep K, iHWCKmodel 4-11)
   659  src/utils/soilhydraulicsutils.f90     (watcon/hconduc/dhconduc/moiscap/prhead/hcomean/dkmean)
 ─────
  9402  total
```

Excluding `sptabulated.f90` (third-party), the SWAP-specific soil home tree is
~3910 LoC.

### Entry points called from `swap.f90` (post-arc-7 status)

| # | Call site | swap.f90 line | Task | Currently takes `state`? |
|---|-----------|---------------|------|--------------------------|
| 1 | `CalcGrid()` | 188 | — (init grid) | NO — pure CalcGrid, grid dims stay legacy |
| 2 | `soilwater_init(state%soilwater, numnod, numlay)` | 189 | — | YES (sub-record direct, to dodge circular dep) |
| 3 | `atmosphere_init(state%atmosphere)` | 190 | — | YES (sub-record) |
| 4 | `DoTillage(1, state)` | 203 | 1 | YES — A-2.6 windfall (nraida retirement) |
| 5 | `SoilWater(1, state)` | 207 | 1 | YES (state already plumbed through) |
| 6 | `drainage_init(state, config)` | 215 | — | YES |
| 7 | `solute_init(state)` | 216 | — | YES |
| 8 | `heat_init(state)` | 217 | — | YES |
| 9 | `BoundBottom(state)` (per-step) | — | — | YES (ADR 0035) |
| 10 | `SoilWater(2, state)` (per-step) | 337 | 2 | YES |
| 11 | `SoilWater(3, state)` (per-step) | 353 | 3 | YES |
| 12 | `DoTillage(2/3, state)` | 300, 409 | 2/3 | YES |

### Inside the soil home tree — current state plumbing

| Subroutine | Source | State arg? | Reads/writes through state today |
|---|---|---|---|
| `soilwater(task, state)` | `soilhydraulics.f90:842` | YES | `%soilwater` (qtop, qbot, qrot, …), `%atmosphere` (nraidt, melt, ldwet/spev/saev, intr/cumu), `%heat%rfcp`, `%surfacewater%vtair`, `%drainage%qdra` |
| `headcalc(state)` | `soilhydraulics.f90:25` | YES | as above |
| `SoilWaterStateVar(task)` | `soilhydraulics.f90:1214` | **NO** — bare globals | hm1/thetm1/gwlm1/pondm1/kmean(NN+1) — all from `use Variables` |
| `hysteresis()` | `soilhydraulics.f90:1263` | **NO** — bare globals | h/hm1/indeks/cofgen/theta/thetar/thetas/dimoca |
| `calcgwl(state)` | `waterbalance.f90:39` | YES (`intent(in)`) | only reads `%soilwater%gwlinp`; writes `gwl`, `nodgwl`, `pegwl`, etc. via legacy globals |
| `level(swoptlev,…)` | `waterbalance.f90:201` | NO | reads h, z, dz, zbotcp, disnod via use Variables (pure function) |
| `watertable(node,…)` | `waterbalance.f90:269` | NO | reads h, theta, thetas, z, dz |
| `fluxes(state)` | `waterbalance.f90:325` | YES | reads `%soilwater%qtop/qbot/qrot/qrosum`, `%surfacewater%qdrtot`, `%drainage%qdra` |
| `integral(state)` | `waterbalance.f90:379` | YES | reads `%atmosphere%*`, `%soilwater%*`, `%surfacewater%intermediate`, `%drainage`; writes legacy iq*/cq*/iqdo/iqup/wbalance |
| `watstor()` | `waterbalance.f90:823` | NO | reads theta, dz, FrArMtrx; writes volact, volm1 |
| `checkmassbal(…,state)` | `waterbalance.f90:605` | YES (intent(in)) | reads `%atmosphere%intr/cumu/ssnow` |
| `CalcGrid()` | `soilgrid.f90:21` | NO | writes grid dims (legacy retained) |
| `ConvertDiscrVert(part,…,state)` | `soilgrid.f90:178` | YES | reads `%surfacewater%intermediate%inqdra`; everything else via legacy |
| `hconduc(node,h,th,rfcp,[tsoil_node])` | `soilhydraulicsutils.f90:410` | NO — value-passed | falls through to `tsoil_loc=0.0` sentinel when caller omits `tsoil_node` (heat ADR 0034 Task 9 residual) |

**Verdict:** the high-traffic routines (`soilwater`, `headcalc`, `fluxes`,
`integral`, `ConvertDiscrVert`, `checkmassbal`) already have `state`. The
**residual plumbing this arc must add**: `SoilWaterStateVar(task, state)`,
`hysteresis(state)`, `watstor(state)`, `level(state, …)` (or pass scalars
explicitly), `watertable(state, …)`. `calcgwl` already accepts state but only
reads `%soilwater%gwlinp`; it needs to be flipped to `intent(inout)` so it can
write `state%soilwater%gwl`/`nodgwl`/`pegwl` etc.

---

## 2. Owned-globals residual inventory

Cross-referenced the mega-discovery's ~100 globals against ADRs 0035/0036/0037
and Phase 0 closures. The unmigrated residual is the inventory for THIS arc.

Categorization:
- **F-INST** — flat instantaneous on `soilwater_state_t` (overwrite each step)
- **F-DAY** — flat per-day scalars/arrays (reset on `flDayStart`)
- **INTR** — intermediate cohort `soilwater_intermediate_t` (reset on `flzerointr`)
- **CUMU** — cumulative cohort `soilwater_cumulative_t` (reset on `flzerocumu`)
- **LEGACY** — stays as legacy global (grid dims, save-state, working buffers, soil-physics tables)
- **DEFER** — to a different arc (macropore, irrigation, solute)

### 2.1 Per-node arrays (the Richards interior)

| Field | Shape | Cadence | Group | variables.f90 line | Notes / representative write |
|---|---|---|---|---|---|
| `theta(:)` | `real(8) macp` | I | **F-INST** | 957 | water content; `soilhydraulics.f90:1043,237,443,519,715,1042,1240`; co-write tillage:330; co-write swapoutput:3946 |
| `thetm1(:)` | `real(8) macp` | I (save-state) | F-INST | 960 | former time-level theta; `soilhydraulics.f90:1229` (only write) |
| `thetar(:)` | `real(8) macp` | I (init+hyster) | F-INST | 958 | per-node residual; `soilhydraulics.f90:950`; `hysteresis:1308,1312,1323` |
| `thetas(:)` | `real(8) macp` | I (init+hyster) | F-INST | 959 | per-node saturated; `soilhydraulics.f90:951`; `hysteresis:1307,1324,1328` |
| `h(:)` | `real(8) macp` | I | **F-INST** | 811 | pressure head; `soilhydraulics.f90:123,125,435,439,724,1037,1318`; `hysteresis:1318,1334`; co-write rootextraction (clamp); co-write swapoutput:3947; co-write tillage |
| `hm1(:)` | `real(8) macp` | I (save-state) | F-INST | 822 | former h; written only in `SoilWaterStateVar(1)` at line 1228; read in restore path and hysteresis |
| `q(:)` | `real(8) macp+1` | I | **F-INST** | 893 | inter-comp flux; `soilhydraulics.f90:881` (init zero); `waterbalance.f90:346,350` (fluxes); read in `integral` |
| `inq(:)` | `real(8) macp+1` | M | **INTR** | 839 | intra-period flux accumulator; `soilhydraulics.f90:1100,1102` (zero reset); `integral:347,358` (accumulate) |
| `k(:)` | `real(8) macp+1` | I | **F-INST** | 855 | hydraulic conductivity; `soilhydraulics.f90:130,170,176,181,238,453,467,471,520,1052,1172,1176` |
| `kmean(:)` | `real(8) macp+1` | I | **F-INST** | 868 | mean K at interface; `soilhydraulics.f90:113,133,186,189,241,262,456,474,475,523,548,1055,1178,1243`; co-write boundtop (k(1) interface); co-write boundbottom:171 |
| `dimoca(:)` | `real(8) macp` | I | **F-INST** | 790 | differential moisture capacity; `soilhydraulics.f90:312,1049`; `hysteresis:1339` |
| `cofgen(21,:)` | `real(8)` | I (init+hyster+tillage) | **F-INST** | 766 | per-node Mualem-VG params; `soilhydraulics.f90:891,907–937,949,950,959,1314–1330`; co-write `tillage.f90:Change_MvGpars` |
| `FrArMtrx(:)` | `real(8) macp` | I (init + macropore overwrite) | **F-INST** | 1234 | matrix-area fraction; `soilhydraulics.f90:1051,1062`; co-write `macropore.f90:435,550` |
| `fluseksatexm(:)` | `logical macp` | I (init only) | F-INST | 639 | per-node Ksatexm flag; only written `soilhydraulics.f90:926`; read in `soilhydraulicsutils:hconduc` path |
| `indeks(:)` | `integer macp` | I (init+hyster) | F-INST | 656 | hysteresis branch index +1/-1; `soilhydraulics.f90:955,959`; `hysteresis:1280,1301` |
| `IThetaBeg(:)` | `real(8) macp` | M | **INTR** | — | theta at start of intr period; `soilhydraulics.f90:1127` (assign on flzerointr) |
| `evp(:)` | `real(8) macp` | I (init-zero only) | LEGACY (or F-INST) | 800 | per-node evaporation; always zero today — kept as dead field |
| `qpotrot_day(:)` | `real(8) macp` | per-day | **F-DAY** | — | already migrated? Check: still in waterbalance:428 as legacy — write through state%soilwater%qpotrot |

**Per-node count:** 17 owned arrays (`evp` is effectively dead; flag for
removal but keep semantics-only migration). `qpotrot_day` was flagged by
mega-discovery but is currently per-day reset of `qpotrot` (already in state).
**Resolution:** `qpotrot_day` is a separate scratch accumulator distinct from
state%soilwater%qpotrot — it sums daily reduction; migrate to F-DAY cohort.

### 2.2 Per-layer arrays (init-once)

`thetsl(maho)` — written only in `soilhydraulics.f90:914,940,942`; read in
hysteresis. Init-once per simulation. Candidate F-INST or LEGACY. Recommend
**F-INST** (it's a soil-water-derived saturated water content per layer).

| Field | Shape | Cadence | Group | Notes |
|---|---|---|---|---|
| `thetsl(:)` | `real(8) maho` | I (init only) | **F-INST** | per-layer saturated water content; alternative: LEGACY |

### 2.3 Top-level scalars (instantaneous)

| Field | Type | Cadence | Group | Description |
|---|---|---|---|---|
| `pond` | real(8) | I | **F-INST** | surface ponding; soil:`1031,1033,1067,1244,758`; **co-writers**: `boundtop.f90:147,166,185,258,266,274,290,308` (8 sites), `tillage.f90:304,330`, `swapoutput.f90:3944` mini-sim |
| `pondm1` | real(8) | I (save-state) | F-INST | pond at former level; `soilhydraulics.f90:1232,758` (only writes) |
| `pondini` | real(8) | I (rebase on flzerocumu) | F-INST | pond at start of cumu period; `soilhydraulics.f90:1068,1157` |
| `gwl` | real(8) | I | **F-INST** | groundwater level; `waterbalance.f90:calcgwl:53,68,79,83,95,97,99,102`; `soilhydraulics.f90:757,1010,1019,1029,1244`; co-writer `swapoutput.f90:3943` mini-sim restore |
| `gwlm1` | real(8) | I (save-state) | F-INST | gwl at former level; `soilhydraulics.f90:1231` |
| `nodgwl` | integer | I | F-INST | node directly above gwl; `calcgwl:54,70,78,82,102` |
| `pegwl` | real(8) | I | F-INST | perched gwl; `calcgwl:53,135,146,148,151` |
| `bpegwl` | integer | I | F-INST | node at bottom of perched gwl; `calcgwl:122,156` |
| `npegwl` | integer | I | F-INST | node above perched gwl; `calcgwl:134,153,157` |
| `gwlflcpzo` | real(8) | I | F-INST | capillary-zone gwl; `calcgwl:63,85,107` |
| `nodgwlflcpzo` | integer | I | F-INST | node for above; `calcgwl:62,82,104` |
| `hatm` | real(8) | I (init only) | **F-INST** | air pressure head near surface; `soilhydraulics.f90:870` (init=−2.75e5); read in boundtop |
| `volact` | real(8) | I | **F-INST** | current soil-profile water storage; `watstor:830-833` |
| `volm1` | real(8) | I (save-state) | F-INST | volact at former level; `watstor:830` |
| `volini` | real(8) | I (rebase on flzerocumu) | F-INST | storage at cumu-period start; `soilhydraulics.f90:1066,1156` |
| `wbalance` | real(8) | I (recomputed in `integral`) | F-INST | cumulative water balance error; `integral:558,563` |
| `fllowgwl` | logical | I | F-INST | gwl-below-profile flag; `soilhydraulics.f90:106,154` |
| `runon` | real(8) | I (init reset) | F-INST | runon flux; `soilhydraulics.f90:878` |

**Scalar count:** 17 instantaneous scalars + 1 init-only sentinel (`hatm`) = 18.

### 2.4 Per-day cohort (flDayStart reset)

Mega-discovery flagged 7 per-day fields. Crop-uptake already migrated some of
the underlying state. Residual:

| Field | Type | Group | Notes |
|---|---|---|---|
| `tra` | real(8) | F-DAY | daily actual transp; `integral:420-421`; reset on `flDayStart` |
| `iqredwet_day` | real(8) | F-DAY | per-day stress accumulator; `soilhydraulics.f90:1085`; `integral:440` |
| `iqreddry_day` | real(8) | F-DAY | per-day; `1086`; `integral:441` |
| `iqredsol_day` | real(8) | F-DAY | per-day; `1087`; `integral:442` |
| `iqredfrs_day` | real(8) | F-DAY | per-day; `1088`; `integral:443` |
| `iptra_day` | real(8) | F-DAY | per-day; `1089`; `integral:445` |
| `qpotrot_day(:)` | `real(8) macp` | F-DAY | per-day per-node; `1091`; `integral:428` |
| `qredtot_day(:)` | `real(8) macp` | F-DAY | per-day per-node; `1092`; `integral:429` |

**Per-day count:** 8 fields. Per the playbook these are small enough to fold
into a `per_day_t` cohort OR into the existing `intr` cohort with a different
reset gate. Recommendation: **fold into INTR cohort with comment** describing
the `flDayStart` gate (small cohort doesn't justify its own type — keeps
soilwater_state_t shape parallel to atmosphere_state_t).

Alternative: **separate `per_day_t` cohort** for cleanliness, mirroring how
surfacewater split off its sub-records post-ADR 0033. Recommend final decision
during design phase.

### 2.5 Intermediate cohort (flzerointr reset)

| Field | Type | Notes |
|---|---|---|
| `inq(:)` | `real(8) macp+1` | already listed in 2.1; reset at `1100,1102`; accumulate in `integral` |
| `inqrot(:)` | `real(8) macp` | `1098`; `integral:427` |
| `inqssdi(:)` | `real(8) macp` | `1099`; `integral:432` — note: also written by SSDI (macropore-coupled); recommend keep here |
| `IThetaBeg(:)` | `real(8) macp` | `1127` |
| `IPondBeg` | real(8) | `1125` |
| `iqrot` | real(8) | `1103`; `integral:424` |
| `iqssdi` | real(8) | `1104`; `integral:433` |
| `iqredwet/dry/sol/frs` | real(8) ×4 | `1105–1108`; `integral:436–439` |
| `ies0, iet0, iew0` | real(8) ×3 | reference evap totals; `1109–1111`; `integral:446–448` |
| `iintc` | real(8) | `1112`; `integral:471` (note: caintc owner is atmosphere; iintc stays soil-water for now) |
| `iruno` | real(8) | `1116`; `integral:477` |
| `irunon` | real(8) | `1121`; `integral:478` |
| `irunoCN` | real(8) | `1117`; `integral:512` |
| `iqbot` | real(8) | `1118`; `integral:485` |
| `iqtdo` | real(8) | `1119`; `integral:487` |
| `iqtup` | real(8) | `1120`; `integral:489` |
| `iqdo(:)` | `real(8) macp+1` | `1122`; `integral:493` |
| `iqup(:)` | `real(8) macp+1` | `1123`; `integral:495` |
| `iprec` | real(8) | `integral:390` (zeroed on flzerointr, accumulated at :480) |
| `igird` | real(8) | `integral:391, 482` |
| `inird` | real(8) | `integral:392, 484` |

**INTR field count:** 22 fields (3 arrays + 19 scalars).

Of these:
- `igird`/`inird` are irrigation accumulators — co-reset in `irrigation.f90`;
  arguably belong to a future irrigation arc. Recommend **keep in soil-water
  INTR cohort** for this arc (defer ownership move to irrigation arc).
- `iprec` is precipitation accumulator — co-reset somewhere atmospheric? Check:
  only soil-water resets `iprec`; atmosphere arc moved `igrai/inrai/isubl/etc`
  out, but `iprec = igrai + igird` is a derived total still owned here.

### 2.6 Cumulative cohort (flzerocumu reset)

| Field | Type | Notes |
|---|---|---|
| `cqssdi` | real(8) | `soilhydraulics.f90:1136`; `integral:500` |
| `cqrot` | real(8) | `1137`; `integral:501`; widely read in swapoutput |
| `cqbot` | real(8) | `1138`; `integral:529` |
| `cqbotdo` | real(8) | `1139`; `integral:525` |
| `cqbotup` | real(8) | `1140`; `integral:527` |
| `cinund` | real(8) | `1144`; `integral:508` |
| `crunon` | real(8) | `1145`; `integral:547` |
| `crunoff` | real(8) | `1146`; `integral:510` |
| `crunoffCN` | real(8) | `1147`; `integral:513` |
| `cqtdo` | real(8) | `1148`; `integral:549` |
| `cqtup` | real(8) | `1149`; `integral:551` |
| `cqprai` | real(8) | `1150`; `integral:546` |
| `cgird` | real(8) | `irrigation.f90:93-95` + `integral:521` — co-owned with irrigation arc |
| `cnird` | real(8) | `irrigation.f90:93-95` + `integral:522` — co-owned with irrigation arc |
| `pondini` | real(8) | rebase at `1157` (already in 2.3 — appears in both groups; recommend keep flat) |
| `volini` | real(8) | rebase at `1156` (already in 2.3) |

**CUMU field count:** 14 (12 cumulative + 2 rebase-on-cumu scalars; the rebase
ones can stay flat with logic in the reset path or move into the cohort —
recommend keep flat for clarity, the rebase is `volini = volact` not zeroing).

### 2.7 Legacy globals (stays — not migrated this arc)

| Field | Why it stays |
|---|---|
| `numnod`, `numlay`, `dz(:)`, `z(:)`, `disnod(:)`, `layer(:)`, `botcom(:)`, `nod1lay(:)`, `ztopcp(:)`, `zbotcp(:)`, `inpola(:)`, `inpolb(:)` | grid dimensions — heat ADR 0034 precedent: keep legacy for cross-subsystem grid access. Move to `grid_t` is a separate future arc. |
| `numtab(:)`, `ientrytab(:,:)`, `sptab(:,:,:)` | soil-physics tabulated functions (init-only, read at every step); legacy aligns with sptabulated.f90 third-party body |
| `numnodNew`, `dzNew`, `inqNew`, `thetaNew`, `hNew`, etc. | ConvertDiscrVert working buffers (caller-owned passed-in arrays); not state |
| `flwarn_hc`, `iwarn_hc`, `nstep_hc` | headcalc warn/step counters; small, write only in soilhydraulics; not external. Could be F-INST but de-prioritized — leave legacy. |
| `fldecdt` | timestep-control flag — owned by timestep_control_mod (already a typed module) |
| `dtold` | `headcalc:702` macropore-rate stash; stays legacy. |
| `iqdrainout(:)` | SwSWST already owns this on `state%surfacewater`; not soil-water. |

### 2.8 Deferred to a different arc

| Field | Defer to | Why |
|---|---|---|
| `cQMpLatSs` (cumu) | **macropore arc** | obviously macropore-owned; `macropore.f90:1548` writes via state%soilwater%QMpLatSs read; `macropore.f90:2208` resets. Mega-discovery flagged as macropore. **THIS ARC: leave legacy; flag in design phase.** |
| `qimmob(:)` | macropore arc | written nowhere in soil home tree; read at `waterbalance:fluxes:350` |
| `QExcMpMtx(:)`, `QMaPo`, `QRapDra`, `ArMpSs`, `ArMpTp`, `dFdhMp(:)`, `QExcMpMtx(:)`, `IcTopMp`, etc. | macropore | macropore-owned |
| `qssdi(:)`, `qssdisum`, `cqssdi`, `iqssdi`, `inqssdi` | SSDI/irrigation arc | written by `irrigation.f90` (SSDI path); SOIL-WATER reads them in integral/fluxes. Recommend **keep cqssdi/iqssdi/inqssdi in soil-water cumu/intr cohorts** (they are the soil-water accumulators), `qssdi/qssdisum` defer to irrigation arc. |
| `swbotb=-2` runtime mutation | NO ACTION (kept as legacy documentation) | `boundbottom.f90:97-99` toggles swbotb between -2/2 as a "subroutine-internal mode flag"; boundary D12 deferral; keep behavior, document. |

### 2.9 Summary counts (this arc)

| Group | Count | Examples |
|---|---|---|
| **F-INST** flat instantaneous | ~30 | theta, h, q, k, kmean, dimoca, cofgen, FrArMtrx, pond, gwl, hatm, volact, volini, pondini, wbalance, runon, fllowgwl, hm1, thetm1, pondm1, gwlm1, thetar, thetas, thetsl, fluseksatexm, indeks, nodgwl, pegwl, bpegwl, npegwl, gwlflcpzo, nodgwlflcpzo, volm1 |
| **F-DAY** per-day | 8 | tra, iqredwet/dry/sol/frs_day, iptra_day, qpotrot_day(:), qredtot_day(:) |
| **INTR** intermediate cohort | 22 | inq(:), inqrot(:), inqssdi(:), IThetaBeg(:), IPondBeg, iqrot, iqssdi, iqredwet/dry/sol/frs, ies0/iet0/iew0, iintc, iruno, irunon, irunoCN, iqbot, iqtdo, iqtup, iqdo(:), iqup(:), iprec, igird, inird |
| **CUMU** cumulative cohort | 14 | cqssdi, cqrot, cqbot, cqbotdo/up, cinund, crunon, crunoff, crunoffCN, cqtdo, cqtup, cqprai, cgird, cnird |
| **LEGACY** (grid + tabulated) | ~12 | numnod, dz, z, disnod, layer, botcom, nod1lay, ztopcp, zbotcp, inpola, inpolb, numtab, sptab |
| **DEFER** (macropore/SSDI/etc.) | ~5 | cQMpLatSs, qimmob, qssdi family |

**Total this arc:** **~74 owned residual fields** to migrate (30 F-INST + 8
F-DAY + 22 INTR + 14 CUMU). The mega-discovery's "~100 globals" expands once
per-day fields (8) and intermediate sub-arrays (3) are counted as separate
fields — the residual count is consistent with `100 - 75 ≈ 25`-ish for the
**unique scalars** but the with-arrays count is 74.

---

## 3. External readers — 5-category inventory

The 23 distinct external files that import any soil-water owned residual
(grep against `theta|gwl|cofgen|qbot|pond|h|q|inq|kmean|FrArMtrx|hm1|thetm1|hatm|volact|volini|wbalance|dimoca|nodgwl|pegwl|bpegwl|npegwl|tra|runon|evp|fllowgwl` minus already-state fields).

| # | File | Cats | Owned soil-water residuals consumed (representative) |
|---|------|------|--------|
| 1 | `src/io/swapoutput.f90` | 1 Output + **co-writer (mini-sim restore)** | theta, h, q, qbot, gwl, pond, kmean, cofgen, FrArMtrx, dimoca, hm1, volact, volini, cqrot, cqbot, iqrot, etc. — **the heaviest reader** (~50 distinct fields); mini-sim writeback (lines 3865-3948) saves/restores qbot/gwl/pond/theta/h |
| 2 | `src/io/swap_csv_output.f90` | 1 Output | theta, h, gwl, pond, kmean, FrArMtrx, q, inq, iqbot, iqrot, irunon, iruno — bulk read for CSV emit |
| 3 | `src/io/macroporeoutput.f90` | 1 Output | FrArMtrx (via macropore), gwl, dz — minor |
| 4 | `src/drainage/drainage.f90` | 2 Compute | gwl, pond, h, theta — drainage flux |
| 5 | `src/drainage/surfacewater.f90` | 2 Compute | gwl, pond, theta |
| 6 | `src/drainage/divdra.f90` | 2 Compute | gwl, h, theta |
| 7 | `src/atmosphere/meteoday.f90` | 2 Compute + Co-write (legacy cumu reset; partially closed by ADR 0037) | theta, ThetaRef; cumu reset of cgrai/cnrai/caintc now owned by atmosphere |
| 8 | `src/atmosphere/et.f90` | 2 Compute | pond, gwl (potential E/T inputs) |
| 9 | `src/boundary/boundtop.f90` | 2 + **co-writer** | h, hatm, kmean(1), pond, pondm1; **writes pond + kmean(1)** at 11 sites |
| 10 | `src/boundary/boundbottom.f90` | 2 + **co-writer** | hbot, gwl, theta, kmean(numnod+1); **writes kmean(numnod+1)** at 2 sites |
| 11 | `src/heat/temperature.f90` | 2 Compute | theta, thetm1, thetas |
| 12 | `src/heat/frozencond.f90` | 2 + co-writer | theta, thetas, gwl, qbot |
| 13 | `src/crop/rootextraction.f90` | 2 + co-writer | theta, h, thetar, thetas, hm1, k, kmean (state-owned writes for qrot already in ADR 0036) |
| 14 | `src/crop/tillage.f90` | 2 + **co-writer** | h, theta, dz, layer, cofgen, pond; writes cofgen (Change_MvGpars), theta+h+pond (Adapt_WC_H) |
| 15 | `src/crop/cropgrowth.f90` | 2 Compute | theta, h |
| 16 | `src/crop/oxygenstress.f90` | 2 Compute | theta, thetas, h, dimoca, cofgen |
| 17 | `src/crop/management_soil.f90` | 2 Compute | theta, h, dz |
| 18 | `src/crop/irrigation.f90` | 2 + co-writer (cumu reset) | h, theta; resets cgird/cnird on flzerocumu (multi-owner pattern; STAYS as legacy or include subset reset in soil-water cohort) |
| 19 | `src/solute/solute.f90` | 2 Compute | theta, h, q, inq, gwl, FrArMtrx |
| 20 | `src/solute/agetracer.f90` | 2 Compute | theta, q, h, gwl, pond |
| 21 | `src/macropore/macropore.f90` | 2 + **co-writer** | h, theta, cofgen, kmean, FrArMtrx; **writes FrArMtrx** (435, 550); reads cQMpLatSs |
| 22 | `src/macropore/macrorate.f90` | 2 Compute | theta, h |
| 23 | `src/utils/surfacewaterutils.f90` | 2 Compute | pond, rsro, rsroexp |

Additional reader/binding files:
- `src/core/swap.f90` — orchestration (passes state)
- `src/utils/soilhydraulicsutils.f90` — INTERNAL home-tree (uses cofgen, numtab, sptab, …)
- `src/utils/sharedexchange.f90`, `src/utils/sharedsimulation.f90` — DLL/shared-state wrappers (probably read theta, gwl, etc. — defer audit)
- `src/io/toml/config_to_variables.f90` — init seed (4 Init-seed)

**Distinct external readers: 23 files** (vs heat 15, surfacewater 12, drainage
~10, atmosphere 9). The mega-discovery counted 25 — the 2-file delta is the
crop-uptake/atmosphere arcs having already moved some readers off legacy.

**Per-field heat-map (the hottest fields):**

| Field | Reader files | Cat |
|---|---|---|
| `theta` | ~14 external | Compute (heat, root, solute, macropore, oxygen, agetracer, divdra, management_soil, meteoday, frozencond, et, irrigation, tillage, swapoutput) |
| `h` | ~12 external | Compute (rootextraction is the heaviest user — JvL + Feddes path) |
| `gwl` | ~9 external | drainage, output, frozencond, agetracer, et, surfacewater, swapoutput, swap_csv_output, divdra |
| `pond` | ~12 external | output, drainage, surfacewater, surfacewaterutils, et, agetracer, solute, swap_csv_output, swapoutput, meteoday, tillage, boundtop |
| `cofgen` | ~6 external | tillage (write), oxygenstress, irrigation, macropore, boundbottom, swapoutput |
| `kmean` | 2 external | boundtop, boundbottom (both write) |
| `FrArMtrx` | 5 external | macropore (write), boundtop, boundbottom, swap_csv_output, swapoutput |
| `cqrot` | 1 external | swapoutput (Section 2.6 cumu cohort field) |
| `dimoca` | 1 external | oxygenstress |
| `hatm` | 1 external | boundtop (only reader; soil-water's only writer) |

**Phase 2 expected compile-driven hidden readers:** **20-25 sites** — this is
the largest reader inventory of any subsystem arc; reduction-only writes are
common and easy to convert, but mechanical site count is large.

---

## 4. Co-writers — files that WRITE soil-water owned residuals

| Co-writer file | Field(s) written | State-arg today? | Disposition |
|---|---|---|---|
| `src/boundary/boundtop.f90` | `pond` (8 sites), `kmean(1)` (3 sites), `qtop` (state-side, already migrated), `reva` (state-side, already migrated) | YES — already takes state | Co-writer for pond + kmean(1); needs `state%soilwater%pond = …` retarget |
| `src/boundary/boundbottom.f90` | `kmean(numnod+1)` (2 sites), `qbot` (state-side), `hbot` (state-side) | YES | Retarget kmean(numnod+1) writes to state%soilwater%kmean(numnod+1) |
| `src/crop/tillage.f90` | `cofgen` (Change_MvGpars), `theta`+`h`+`pond` (Adapt_WC_H) | YES — A-2.6 windfall (DoTillage now takes state) | Plumb state into Change_MvGpars / Adapt_WC_H internal subs; write through state%soilwater |
| `src/crop/rootextraction.f90` | `h(node)` clamp at ~line 382 (potential mutate); qrot/qrosum already in state per ADR 0036 | YES | retarget the h-clamp to state%soilwater%h |
| `src/macropore/macropore.f90` | `FrArMtrx(:)` (lines 435, 550) | YES (state passed for macropore tasks 1-6) | Retarget FrArMtrx writes to state%soilwater%FrArMtrx |
| `src/heat/frozencond.f90` | `qbot` (already in state per ADR 0035); reads theta, gwl | YES | Reads only — no new co-write |
| `src/io/swapoutput.f90` | mini-sim writeback: `gwl`, `pond`, `theta(:)`, `h(:)` at lines 3943-3948; `qbot` already retargeted in ADR 0035 | YES (state_main) | The mini-sim snapshot path (3865-3948) needs retarget for all 4 residual fields |
| `src/atmosphere/meteoday.f90` | Was a co-writer of cgrai/cnrai/caintc (mega-discovery flagged); ADR 0037 moved atmosphere cohort to state%atmosphere. **Remaining soil-water co-write: NONE** — meteoday only resets atmosphere cumus now. | YES | No action (resolved by atmosphere arc) |
| `src/crop/irrigation.f90` | `cgird`, `cnird` cumu reset (`irrigation.f90:93-95`) | YES | Same multi-owner pattern; recommend keep cgird/cnird in soil-water CUMU cohort and add subset-reset call from irrigation.f90, OR move ownership to a future irrigation state |

**Total co-writers: 7 distinct files** (down from mega-discovery's 10; atmosphere
ADR 0037 closed meteoday + collapsed atmosphere ownership). Highest impact:
**tillage** (mutates cofgen+theta+h+pond at 3 tasks) and **swapoutput mini-sim**
(snapshot/restore of 5 fields).

---

## 5. Init-order analysis (swap.f90)

Current init order (verified):

```
182:   ...
188:   call CalcGrid()                                          ! grid dims (numnod, numlay, dz, z, …)
189:   call soilwater_init(state%soilwater, numnod, numlay)     ! ADR 0035 + 0036
190:   call atmosphere_init(state%atmosphere)                   ! ADR 0037
203:   if (flTillage) call DoTillage(1, state)                  ! tillage init (after grid; before soilwater(1))
207:   call SoilWater(1, state)                                 ! Richards init: cofgen, h, theta, gwl, ...
215:   call drainage_init(state, config)
216:   if (flSolute) call solute_init(state)
217:   call heat_init(state)
```

**Verdict on `soilwater_init` extension:**
- `numnod` and `numlay` are passed today — they are written by `CalcGrid` at
  line 188, available at line 189. ✓
- **New per-node arrays to allocate this arc:** `theta(numnod)`, `thetm1(numnod)`,
  `thetar(numnod)`, `thetas(numnod)`, `thetsl(numlay)`, `h(numnod)`, `hm1(numnod)`,
  `q(numnod+1)`, `inq(numnod+1)`, `k(numnod+1)`, `kmean(numnod+1)`, `dimoca(numnod)`,
  `cofgen(21, numnod)`, `FrArMtrx(numnod)`, `fluseksatexm(numnod)`, `indeks(numnod)`,
  `IThetaBeg(numnod)`, `iqdo(numnod+1)`, `iqup(numnod+1)`, `inqrot(numnod)`,
  `inqssdi(numnod)`, `qpotrot_day(numnod)`, `qredtot_day(numnod)`, `evp(numnod)`.
- All bounds (`numnod`, `numlay`) available at line 189 — no init-order change needed.
- **Tillage timing concern:** `DoTillage(1, state)` at line 203 mutates cofgen
  and (later via Adapt_WC_H) theta/h/pond — BUT `SoilWater(1, state)` at line
  207 writes those same fields via `cofgen = 0.0d0` + initial-h profile. After
  migration, tillage-1 must consume `state%soilwater%cofgen` and the post-tillage
  rewrite by SoilWater(1) is reset / re-written. Verify the existing init
  ordering preserves the legacy semantics. Discovery flags this as a
  **design-phase verification point** — see hazard #8.

**Tillage windfall:** `DoTillage(iTask, state)` already accepts state thanks to
the A-2.6 ADR 0037 windfall (`nraida` retirement). The tillage-side cofgen/theta/h/pond
co-writes can now be retargeted with no further plumbing change. ✓

---

## 6. Config / Phase 0 candidates

Mega-discovery flagged 5 candidates. All 5 are **already typed-config covered**:

| Candidate | Status |
|---|---|
| `hplate` | ✓ `src/config/bottom_boundary_config.f90:78` |
| `sinmax, sinamp, sinave` | ✓ `bottom_boundary_config.f90:40-42` |
| `cofqha, cofqhb` | ✓ `bottom_boundary_config.f90:47-48` |
| `gwlconv` | ✓ `src/config/simulation_config.f90:26` |
| `deepgw` | ✓ state-side (`state%soilwater%deepgw`, ADR 0035) |

**Additional sweeps:** Richards-solver criteria (`CritDevh1Cp`, `CritDevh2Cp`,
`CritDevPondDt`, `CritDevBalCp`, `CritDevBalTot`, `CritDevMasBal`, `Critdz`):
some are local `data` constants in `headcalc` (lines 60-65 — local PARAMETER
init); `CritDevh1Cp`, `CritDevh2Cp`, `CritDevPondDt`, `CritDevMasBal` are
elsewhere — `simulation_config` likely covers them. **No Phase 0 work needed
for this arc.**

**Phase 0 candidates: 0** — all surfaced gaps from the mega-discovery were
closed during the three prior arcs.

---

## 7. Coupling hazards

### Hazard #1 — `hconduc` `tsoil_node = 0.0` sentinel (heat ADR 0034 Task 9 residual)

**Scope: small-medium.** `hconduc(node, head, theta, rfcp, [tsoil_node])`
falls through to `tsoil_loc = 0.0` when `tsoil_node` is not passed. The
iHWCKmodel 4-11 path therefore evaluates at tsoil=0°C — an unphysical sentinel,
documented in ADR 0034 as "unreachable in regression."

**14+ callers identified:**
- soilhydraulics.f90: 10 sites (113, 130, 170, 238, 262, 318 (dhconduc), 453, 520, 548, 1052, 1172)
- boundary/boundbottom.f90:171
- boundary/boundtop.f90:111
- crop/rootextraction.f90:868, 873 (passes `10.d0` literal!)
- crop/tillage.f90:361 (passes `1.0d0` literal in a write statement!)
- io/swapoutput.f90:2103 (passes `rfcpx`, no tsoil_node)
- macropore/macropore.f90:1142 (passes `Dum`)

**This arc resolves:** Thread `state%heat%tsoil(node)` into the home-tree callers
(soilhydraulics 10 sites + soilhydraulicsutils as needed). External callers
(rootextraction, tillage, swapoutput, macropore, boundtop, boundbottom) all
already have `state` in scope — mechanically pass `state%heat%tsoil(node)`.
Single-coordinated task replaces 14 sites.

**Recommendation:** dedicated task (e.g. SS-SWC Task X) for hconduc
`tsoil_node` threading — small but touches many files.

### Hazard #2 — Mini-simulation writeback (swapoutput.f90 lines 3865-3948)

**Scope: small.** The OutputModFlow mini-sim:
- Saves: `qbottmp = state%soilwater%qbot` (ADR 0035 done), `gwltmp = gwl`,
  `pondtmp = pond`, `thetatmp(:) = theta(:)`, `htmp(:) = h(:)` (lines 3866-3873)
- Runs a 2-arm SWAP perturbation with `state_om` to compute storage coefficient
- Restores: same 5 fields (lines 3941-3948)

**This arc retargets:** `gwl`, `pond`, `theta(:)`, `h(:)` snapshot+restore
must move to `state%soilwater%X`. Pattern is identical to ADR 0035's `qbot`
retarget. Mechanical — 10 sites, one block.

### Hazard #3 — Tillage co-write of `cofgen`, `theta`, `h`, `pond`

**Scope: medium.** `tillage.f90:DoTillage(iTask, state)` is now state-plumbed
(A-2.6 windfall), but internal subs `Change_MvGpars`, `Adapt_WC_H`,
`Change_Bdens`, `Consolidate_Bdens` don't take state — they `use Variables`
directly. Migration must thread state into all four internal subs OR keep them
internal to DoTillage and pass state%soilwater%cofgen / theta / h / pond explicitly.

**Plus an init-order subtlety (see Section 5):** `DoTillage(1, state)` runs at
swap.f90:203 BEFORE `SoilWater(1, state)` at line 207. SoilWater(1) overwrites
cofgen + theta + h + pond at init. Tillage(1)'s effect at init is "set the
target rho values; the actual cofgen rewrite happens later when SoilWater
runs." Verify: legacy semantics use `Bdens` from tillage's storage to inform
SoilWater(1)'s cofgen build (via PdmVG indirect read). After migration, the
data path must still work.

**Recommendation:** task to migrate tillage internal subs (Change_MvGpars +
Adapt_WC_H) — design-phase verify they see state%soilwater after the SoilWater(1)
seed.

### Hazard #4 — Macropore co-write of `FrArMtrx`

**Scope: small.** `macropore.f90:435,550` writes `FrArMtrx(:)` directly. State
is already passed to macropore tasks 1-6. Retarget to `state%soilwater%FrArMtrx`.
Mechanical 2-site change. SoilWater(1) also writes `FrArMtrx(i) = 1.d0` as
non-macropore fallback (line 1051, 1062) — retarget too.

### Hazard #5 — ConvertDiscrVert legacy-global reads (soilgrid.f90:178-441)

**Scope: medium.** `ConvertDiscrVert(part, swop, ..., state)` already accepts
state (used for state%surfacewater%intermediate%inqdra), but the body reads
**~10 soil-water owned globals via `use Variables`**:
- theta, h, inq, IThetaBeg, cofgen, dz, FrArMtrx (line 186-188)
- inqrot — to be retargeted

When these migrate, the use-Variables list shrinks to grid-only (numnod, numlay,
dz, layer, botcom, etc.) and the body reads state%soilwater%X. ~30 sites inside
a 263-line subroutine. Mechanical but voluminous.

### Hazard #6 — SoilWaterStateVar(task) needs state plumbing

**Scope: small.** Current signature: `SoilWaterStateVar(task)` with bare use
Variables — writes hm1/thetm1/gwlm1/pondm1 (task=1) and restores h/theta/gwl/
pond/kmean(numnod+1) (task=2). Three callers:
- `soilhydraulics.f90:soilwater:1161` (per-step task=1)
- `swapoutput.f90:3903, 3915` (mini-sim task=1 and task=2)

Add state arg: `SoilWaterStateVar(task, state)`. Three call sites.

The task=2 restore path includes a `kmean(numnod+1) = k(numnod)` line — an
outside-the-m1-copy-set write that needs to be retained as
`state%soilwater%kmean(numnod+1) = state%soilwater%k(numnod)`.

### Hazard #7 — `hysteresis()` needs state plumbing

**Scope: small.** `hysteresis()` writes thetar/thetas/cofgen/dimoca/h via bare
use Variables. One caller (soilhydraulics.f90:1194). Add state arg, retarget
~10 sites in body. Routine is ~80 lines — surgical.

### Hazard #8 — `calcgwl(state)` write side

**Scope: small.** Already has `type(swap_state_t), intent(in) :: state`, but
**writes gwl/nodgwl/pegwl/bpegwl/npegwl/gwlflcpzo/nodgwlflcpzo via bare globals**.
Flip to `intent(inout)`; retarget ~12 write sites to state%soilwater%X.

### Hazard #9 — `watstor()` write side

**Scope: small.** `watstor()` writes volm1+volact. Add state arg, retarget 3 sites.

### Hazard #10 — Multi-owner cumu reset (cgird/cnird)

**Scope: small.** `irrigation.f90:93-95` resets cgird/cnird on flzerocumu —
duplicate of soil-water's `integral` accumulation. Atmosphere ADR 0037
resolved the cgrai/cnrai/caintc path (moved to atmosphere cohort). For
this arc, cgird/cnird are residual; recommend **keep in soil-water CUMU
cohort with a comment** noting irrigation also resets (pre-existing).
Alternative: pre-emptively move to a future irrigation_state_t — defer to
irrigation arc.

### Hazard #11 — `swbotb=-2` runtime mutation (boundary D12 carry-forward)

**Scope: NO-OP, documentation only.** `boundbottom.f90:97-99` toggles swbotb
between -2/2 as a subroutine-internal mode flag — pre-existing pattern; kept
as-is per boundary D12 deferral. Do not touch.

### Hazard #12 — `cQMpLatSs` ownership (macropore arc territory)

**Scope: deferred.** mega-discovery flagged this as macropore-owned.
`soilhydraulics.f90:889` writes `cQMpLatSs = 0.0d0` (init reset).
`macropore.f90:1548,2208` is the actual write/reset path. **This arc: leave as
legacy global; design phase decides defer-to-macropore-arc.**

### Hazard #13 — `cnrai`/`cevap` reads in `integral` body (lines 558-565)

**Scope: NONE — already resolved.** ADR 0037 retargeted these reads to
`state%atmosphere%cumu%cnrai` / `state%atmosphere%cumu%cevap`. No further
action.

### Hazard #14 — Phase 2 compile-driven readers (estimate)

**Scope: medium-large.** The 23 external reader files + ~20 read-site count
per heavy field (theta has ~14 readers × ~10 sites each = ~140 mechanical
retarget sites in the hot fields alone). Conservative estimate: **300-400
mechanical retarget sites**. This is the largest Phase 2 of any arc.

**Cohort lesson from atmosphere (ADR 0033 → 0037):** the cohort `reset()` calls
consolidate scattered reset blocks. For soil-water this collapses:
- `soilhydraulics.f90:1083-1132` (flzerointr block, 30+ lines) → `call
  state%soilwater%intr%reset()` (single line) + the per-day mini-block can
  collapse similarly.
- `soilhydraulics.f90:1134-1158` (flzerocumu block, 25+ lines) → `call
  state%soilwater%cumu%reset()` + the volini/pondini rebase outside.

Applied broadly this is **the largest reset-block consolidation of any arc**.

---

## 8. Reset-cadence + cohort decision

Four cadences:

| Cadence | Count | Activity gate | Cohort? |
|---|---|---|---|
| Flat instantaneous | ~30 | none / FlMacropore for FrArMtrx | Direct fields on `soilwater_state_t` |
| Per-day | 8 | flDayStart | **Recommend: fold into INTR cohort with comment** (parallels heat ADR 0034's "flat" choice for a small group); OR new `per_day_t` for cleanliness |
| Intermediate | 22 | flzerointr | **NEW soilwater_intermediate_t cohort** |
| Cumulative | 14 | flzerocumu | **NEW soilwater_cumulative_t cohort** |

**Recommendation: TWO COHORTS** (intermediate + cumulative) mirroring
atmosphere ADR 0037's `intr` + `cumu` pattern. Per-day fields go into a
**third small `per_day_t` cohort** OR fold into `intr` with a comment.

The atmosphere ADR 0037 precedent (TWO cohorts) is the best template. Soil-water's
per-day cohort is small enough that designing a 3rd cohort is overkill — but
the per-day reset gate is *distinct* (`flDayStart` not `flzerointr`), so
folding into INTR would mean adding a `reset_per_day()` subroutine to
`soilwater_intermediate_t` that zeros only the 8 per-day fields. **Acceptable
trade-off.**

### Proposed `soilwater_state_t` layout (post-arc)

```fortran
module soilwater_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   public :: soilwater_state_t, soilwater_intermediate_t, soilwater_cumulative_t, soilwater_init

   type :: soilwater_intermediate_t
      ! per-node intra-period flux accumulators
      real(real64), allocatable :: inq(:), inqrot(:), inqssdi(:)
      real(real64), allocatable :: iqdo(:), iqup(:), IThetaBeg(:)
      ! scalars
      real(real64) :: iqrot, iqssdi
      real(real64) :: iqredwet, iqreddry, iqredsol, iqredfrs
      real(real64) :: ies0, iet0, iew0, iintc
      real(real64) :: iruno, irunoCN, irunon
      real(real64) :: iqbot, iqtdo, iqtup
      real(real64) :: IPondBeg
      real(real64) :: iprec, igird, inird
      ! per-day cohort (flDayStart gate, separate from flzerointr)
      real(real64) :: tra
      real(real64) :: iqredwet_day, iqreddry_day, iqredsol_day, iqredfrs_day, iptra_day
      real(real64), allocatable :: qpotrot_day(:), qredtot_day(:)
   contains
      procedure :: reset => soilwater_intermediate_reset       ! zero all (flzerointr)
      procedure :: reset_per_day => soilwater_per_day_reset    ! zero the 8 per-day (flDayStart)
   end type

   type :: soilwater_cumulative_t
      real(real64) :: cqssdi, cqrot, cqbot, cqbotdo, cqbotup
      real(real64) :: cinund, crunon, crunoff, crunoffCN
      real(real64) :: cqtdo, cqtup, cqprai, cgird, cnird
   contains
      procedure :: reset => soilwater_cumulative_reset
   end type

   type :: soilwater_state_t

      ! ----- existing 35 fields from boundary + crop-uptake arcs -----
      ! (qtop, qbot, qbot_nonfrozen, hbot, gwlinp, deepgw, reva, hsurf, ftoph, runots, QMpLatSs, FlRunoff)
      ! (qrot, qpotrot, qredwet, qreddry, qredsol, qredfrs, mflux, mroot, hroot, rootrho, rootphi, rmax, mfluxtable, ...)

      ! ----- NEW: Richards interior per-node arrays -----
      real(real64), allocatable :: theta(:), thetm1(:), thetar(:), thetas(:), thetsl(:)
      real(real64), allocatable :: h(:), hm1(:)
      real(real64), allocatable :: q(:), k(:), kmean(:)
      real(real64), allocatable :: dimoca(:)
      real(real64), allocatable :: cofgen(:,:)         ! (21, numnod)
      real(real64), allocatable :: FrArMtrx(:)
      real(real64), allocatable :: evp(:)              ! retained for completeness
      logical, allocatable      :: fluseksatexm(:)
      integer, allocatable      :: indeks(:)           ! hysteresis branch

      ! ----- NEW: scalars -----
      real(real64) :: pond     = 0.0_real64
      real(real64) :: pondm1   = 0.0_real64
      real(real64) :: pondini  = 0.0_real64
      real(real64) :: gwl      = 0.0_real64
      real(real64) :: gwlm1    = 0.0_real64
      real(real64) :: pegwl    = 999.0_real64
      real(real64) :: gwlflcpzo = 999.0_real64
      real(real64) :: hatm     = -2.75e5_real64
      real(real64) :: volact   = 0.0_real64
      real(real64) :: volm1    = 0.0_real64
      real(real64) :: volini   = 0.0_real64
      real(real64) :: wbalance = 0.0_real64
      real(real64) :: runon    = 0.0_real64
      integer :: nodgwl       = 0
      integer :: bpegwl       = -1
      integer :: npegwl       = -1
      integer :: nodgwlflcpzo = 0
      logical :: fllowgwl     = .false.

      ! ----- NEW: cohorts -----
      type(soilwater_intermediate_t) :: intr
      type(soilwater_cumulative_t)   :: cumu
   end type

contains
   subroutine soilwater_init(sw, numnod, nlay)
      ! existing 35 fields + allocations
      ! + new allocations for theta/h/q/k/kmean/dimoca/cofgen/FrArMtrx/...
      ! + cohort allocations (inq/iqdo/iqup/qpotrot_day/qredtot_day/inqrot/inqssdi/IThetaBeg)
      ! + cohort reset() calls
   end subroutine
end module
```

---

## 9. Scope estimate

| Metric | Soil-water core (THIS) | Atmosphere | Crop-uptake | Boundary |
|---|---|---|---|---|
| Owned globals (residual) | **~74** | 40 | 23 | 12 |
| **Flat instantaneous** | **~30** | 22 (11 + 9 + 3 per-event) | 18 | 12 |
| **Flat per-day** | **8** | 0 | 0 | 0 |
| **Intermediate cohort** | **22** | 8 | 0 | 0 |
| **Cumulative cohort** | **14** | 10 | 0 | 0 |
| **Stays legacy (grid + tabulated)** | ~12 | 0 | 0 | 0 |
| External readers (files) | **23** | 9 | 7 | 14 |
| External read sites (rough) | **300-400** | 155 | 47 | 80 |
| Co-writers | **7** | 3 (post-collapse) | 1 | 5 |
| Phase 0 candidates | **0** | 0 | 0 | 8 |
| Subsystem entry-points needing state plumbing | **5** (SoilWaterStateVar, hysteresis, watstor, level, watertable, plus tillage internals) | 2 | 1 | 3 |
| Cohort design complexity | **2 cohorts + per-day procedure** | 2 cohorts | 0 | 0 |
| hconduc `tsoil_node` threading | **14 sites** | — | — | — |
| Suggested task count | **22-25** | 17 | 11 | 16 |

**Suggested decomposition (22-25 tasks):**
- SS-SWC Phase 1 (state-type + init): 3 tasks
  1. Add soilwater_intermediate_t + soilwater_cumulative_t types
  2. Extend soilwater_state_t with ~30 F-INST fields (arrays + scalars)
  3. Extend soilwater_init: allocate + zero all new fields
- SS-SWC Phase 2 (call-site retarget within soil home): 8 tasks
  4. SoilWaterStateVar(task, state) plumbing
  5. hysteresis(state) plumbing
  6. watstor(state), level(state), watertable(state) plumbing
  7. calcgwl(state) intent(inout) + write retarget
  8. soilhydraulics headcalc body retarget (~80 sites)
  9. soilhydraulics soilwater(1) body retarget (~40 sites)
  10. soilhydraulics soilwater(2/3) reset block consolidation → intr%reset() + cumu%reset()
  11. waterbalance integral body retarget (~50 sites)
- SS-SWC Phase 3 (co-writer retarget): 5 tasks
  12. boundtop pond + kmean(1) writes (11 sites)
  13. boundbottom kmean(numnod+1) writes (2 sites)
  14. tillage cofgen + theta + h + pond co-writes
  15. macropore FrArMtrx writes (2 sites) + soil-water FrArMtrx writes
  16. rootextraction h-clamp write
- SS-SWC Phase 4 (external reader retarget): 6 tasks (grouped by subsystem)
  17. swapoutput.f90 reads + mini-sim writeback (heaviest, ~60 sites)
  18. swap_csv_output + macroporeoutput
  19. drainage + surfacewater + divdra
  20. crop subsystems (cropgrowth, oxygenstress, rootextraction reads, management_soil, irrigation reads, tillage reads)
  21. heat + frozencond + macropore + macrorate + solute + agetracer + et + meteoday + surfacewaterutils
- SS-SWC Phase 5 (hconduc + cleanup): 3 tasks
  22. hconduc `tsoil_node` threading (14 sites)
  23. ConvertDiscrVert body retarget (mechanical inside soilgrid)
  24. Phase 2.6 compile-driven hidden readers (expected 10-20 sites)
  25. Drop legacy globals + verify-state-not-reports

**Final task count estimate: 22-25** — largest arc to date by 30%.

---

## 10. Open questions for design phase

1. **Cohort shape.** Recommendation: **2 cohorts (intr + cumu) mirroring
   atmosphere ADR 0037**, with `intr%reset()` (flzerointr) + a separate
   `intr%reset_per_day()` (flDayStart). 8 per-day fields fold into intr type
   but get their own zero procedure. **Confirmed: this is the most-uniform
   shape.**

2. **`pond` migration — verify tillage windfall.** DoTillage(iTask, state) is
   plumbed (A-2.6). Adapt_WC_H (internal sub of tillage) writes `pond` —
   plumb state through OR pass `state%soilwater%pond` explicitly. Recommendation:
   plumb state through tillage internals.

3. **`gwl` migration — verify calcgwl windfall.** calcgwl(state) already
   exists (intent(in) for state%soilwater%gwlinp). Flip to intent(inout);
   write `state%soilwater%gwl/nodgwl/pegwl/etc.` ~12 sites.

4. **`kmean` migration — verify boundary windfall.** boundtop / boundbottom
   already have state; mechanical retarget of `kmean(1)` and `kmean(numnod+1)`.
   `kmean` is full array, lives in state%soilwater%kmean(:).

5. **`hconduc` `tsoil_node` threading.** Single coordinated task across
   14 sites. All callers have state in scope (post-arc-7).

6. **`cQMpLatSs` ownership.** Defer to macropore arc. Document, leave legacy.

7. **Mini-sim writeback retarget.** swapoutput.f90:3865-3948 must retarget
   gwl/pond/theta/h snapshot+restore (qbot already done in ADR 0035). 10
   sites in one block.

8. **Tillage init-order verification.** DoTillage(1) at swap.f90:203 runs
   BEFORE SoilWater(1) at line 207. SoilWater(1) resets cofgen=0 and writes
   cofgen from PdmVG (which tillage updated). Verify legacy semantics
   preserved when both write through state%soilwater%cofgen.

9. **Grid dimensions (numnod, dz, z, disnod, layer, botcom, nod1lay,
   ztopcp/zbotcp, inpola/inpolb).** **Keep as legacy globals** per heat ADR
   0034 precedent. Future arc may introduce a `grid_t`.

10. **`evp(:)`** — always-zero per-node field. Migrate semantics-only
    (preserve behavior) OR drop. Recommend migrate (zero overhead, preserves
    legacy parity).

11. **`cgird`/`cnird` multi-owner reset.** irrigation.f90:93-95 resets these.
    Keep in soil-water CUMU cohort; let irrigation continue to call subset-reset.
    Defer ownership move to irrigation arc.

12. **`iprec`/`igird`/`inird`** — irrigation accumulators in soil-water's
    INTR cohort. Keep here (defer to irrigation arc).

13. **`swbotb=-2` runtime mutation.** Per boundary D12: NO-OP, kept as legacy.

14. **`thetsl(maho)`** — per-layer; flat F-INST OR LEGACY? Recommend F-INST
    (it's a derived saturated water content per layer, conceptually
    soil-water state).

15. **ConvertDiscrVert body retarget.** The `IThetaBeg`, `inq`, `inqrot`
    reads inside the 264-line body need to point at state%soilwater%intr%X.
    Mechanical but voluminous.

16. **fluseksatexm migration.** init-once logical; flat F-INST. Read by
    hconduc indirectly via cofgen(10) check, plus 4 external sites — verify
    the external reads are still valid post-migration.

17. **Phase 2.6 (compile-driven) expectation.** Heat arc found ~7 hidden
    readers; soil-water expects 20-25 due to volume. Budget time for the
    "compile and see" iteration.

18. **Verification: `check-full` discipline.** Per the feedback skill —
    don't claim done until check-full passes (not just pFUnit). The cohort
    reset() consolidation has historically caused subtle global-default
    regressions.

---

## Discovery summary

- **Subsystem size:** the residual after three coupling-surface arcs is still
  the largest single arc — ~74 owned fields across 4 cadences. The
  Richards-equation core lives here.
- **Reset cadence:** 4 groups (F-INST + F-DAY + INTR + CUMU). **2 cohorts
  (intr + cumu) mirroring atmosphere ADR 0037** is the recommended shape, with
  a `reset_per_day()` method on the intr type for the 8 flDayStart fields.
- **Co-writers:** 7 distinct files. The heavy ones (tillage, swapoutput
  mini-sim, boundtop pond+kmean(1), boundbottom kmean(NN+1), macropore
  FrArMtrx) all have state in scope today.
- **Plumbing residual:** SoilWaterStateVar, hysteresis, watstor, level,
  watertable, calcgwl-intent-flip. All small.
- **hconduc `tsoil_node` sentinel:** resolved in a single 14-site coordinated
  task — all callers have state in scope.
- **Phase 0:** zero new candidates (all 5 mega-discovery items now covered).
- **Suggested task count: 22-25** — largest arc by 30%.
- **Deferred items:**
  - `cQMpLatSs` → macropore arc (legacy preserved)
  - `qssdi`/`qssdisum` → irrigation arc (cqssdi/iqssdi/inqssdi stay here as cumu/intr cohort fields)
  - `cgird`/`cnird` → kept in soil-water CUMU cohort; irrigation co-resets
  - `swbotb=-2` runtime mutation → NO-OP (boundary D12 carry-forward)
- **State-in-scope status:** the high-traffic routines all have state already
  (post-arc-7); the residual plumbing is mechanical.

End of discovery.
