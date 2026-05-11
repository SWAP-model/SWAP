## Subsystem Migration Discovery: Crop Water Uptake

**Date:** 2026-05-10
**Status:** discovery (read-only inventory)
**Migration #:** 6 of N — second of the four coupling-surface arcs (boundary → **crop-uptake** → atmosphere → soil-water core)
**Branch:** `development`
**Scope note:** This arc carves the **per-node root sink term `qrot(:)` and its scalar sum `qrosum`** (plus the four per-node stress-component arrays + their scalar sums + the JvL state `Tactual` / `alpJvLier`) out of legacy `variables.f90` into `state%soilwater`. Resolves the qrot ownership ambiguity flagged in the soil-water mega-discovery (Section 2/7) by giving `rootextraction.f90` undisputed ownership of qrot as soil-water-state.
**Predecessor playbooks:**
- `docs/superpowers/specs/state-migration-playbook.md` (lessons 1–8 boundary, 1–5 heat)
- `docs/superpowers/specs/2026-05-10-state-migration-boundary-discovery.md` (structure template + 4-arc framing)
- `docs/superpowers/specs/2026-05-10-state-migration-boundary-design.md`
- `docs/adr/0035-state-migration-boundary.md`
- `docs/superpowers/specs/2026-05-10-state-migration-soilwater-discovery.md` (superseded mega-discovery — Section 2 inventory + Section 7 hazards remain authoritative)

> Read-only discovery: no code changes. All file:line references are anchors for the design phase.

---

## 1. Big picture

### Subsystem role

`RootExtraction(state)` is the **per-timestep root water uptake compute**. It is the bridge between three subsystem worlds — atmosphere (`ptra` potential transpiration), soil-water (`h`, `theta`, `kmean`), and crop (`rd`, `noddrz`, `cumdens`, `rdctb`) — and produces the per-node sink term consumed by Richards.

Three subroutines + one parameter-table helper, all in module `rootextraction_mod` (`src/crop/rootextraction.f90`):

- **`RootExtraction(state)`** (lines 31–303): the per-timestep entry point. Resets `qrot(1:numnod)`, computes potential uptake by Feddes (`swdrought=1`) or microscopic JongvanLier (`swdrought=2`), then applies multiplicative reductions for oxygen / drought / salt / frost stress, partitions reduction into the four stress components, and (optionally) applies Jarvis/Walsum compensation.
- **`JongvanLier(state)`** (lines 313–640): microscopic root uptake compute — Newton-Raphson search on `hleaf` to match `qrosum = ptra`, called when `swdrought=2`. Writes `mflux`, `mroot`, `hroot`, `rootrho`, `rootphi`, `rmax`, `Hleaf`, `Hxylem`, `alpJvLier`. Calls `MatricFlux(2, ...)` and `JongvanLierLoop`.
- **`JongvanLierLoop(state)`** (lines 649–760): one Newton-loop pass; calls `MatricFlux(2, ...)` per node; writes `hroot`, `mroot`, `qrot`, `qrosum`.
- **`MatricFlux(task, phead, node, outcome, state)`** (lines 775–873): matric-flux-potential lookup. Task 1 builds `mfluxtable` per layer; Task 2 evaluates with optional osmotic correction (reads `state%solute%cml`). `state` is **optional intent(in)** (per playbook gotcha #3 — many call chains, including pre-state cropgrowth init).

### Signature status (the windfall)

| Routine | Takes `state`? | Source |
|---|---|---|
| `RootExtraction(state)` | **YES** (`intent(in)`) | SS-HEAT Task 6 (for `state%heat%tsoil` frost gate at line 193) and SS-SLST (for `state%solute%cml` salt gate at line 181) |
| `JongvanLier(state)` | **YES** (`intent(in)`) | Same arcs |
| `JongvanLierLoop(state)` | **YES** (`intent(in)`) | Same arcs |
| `MatricFlux(..., state)` | **YES — `optional intent(in)`** | SS-HEAT/SS-SLST; called from cropgrowth Task 1 init **without** state (lines 601/1331/2400) |
| `OxygenStress(node, rwu_factor, state)` | **YES — `optional intent(in)`** | SS-HEAT Task 6 |

**No signature surgery required.** Like boundary, this is a field-carve + reader-cutover arc, not a signature-plumbing arc. Boundary lesson #3 applies again.

### Entry points (call graph from outside)

| # | Call site | File:line | Phase |
|---|---|---|---|
| 1 | `RootExtraction(state)` | `swap.f90:302` | Top of each Richards timestep, **after** `FrozenCond` and **before** `BoundBottom` / `Drainage` / `SurfaceWater` / `SoilWater(2)`. Called **once per timestep** (NOT inside `headcalc`). |
| 2 | `MatricFlux(1, h(1), 1, dummy)` | `cropgrowth.f90:601, 1331, 2400` | Task 1 init from each cropgrowth variant (`CropFixed`, `CropWofost`, `CropGrass`) when `swdrought=2`. No `state` arg — uses `optional` absence to skip `state%solute%cml` reads. |
| 3 | `MatricFlux(2, h(...), node, outcome, state)` | `rootextraction.f90:383, 687, 688, 711` (internal) | Per-node JvL evaluation |

**Critical call-order observation:** Unlike boundary's `boundtop` (called inside `headcalc`), `RootExtraction` runs at the **outer timestep level** (swap.f90:302). Its outputs `qrot(:)` and `qrosum` are read by `soilhydraulics.f90:headcalc` (10 read sites) inside the Richards inner loop, but the writes themselves happen once per outer dt. Cleaner semantics than boundary — no per-iteration retargeting concern.

### Lines of code (home + adjacent)

```
  875 src/crop/rootextraction.f90      (home)
 1830 src/crop/oxygenstress.f90        (OxygenStress + OxygenReproFunction — co-stress callees, NOT carved here)
 4696 src/crop/cropgrowth.f90          (MatricFlux task-1 init from CropFixed/Wofost/Grass)
 1343 src/soil/soilhydraulics.f90      (heaviest external reader: 10 qrot read sites)
  816 src/soil/waterbalance.f90        (cumulative qrot accumulator; intermediate qred sums)
```

The home file is ~3× larger than boundary, but the JvL/MatricFlux machinery is largely self-contained per-node compute that does not cross subsystem lines.

---

## 2. Owned-and-touched globals (the inventory)

### 2a. Owned globals — written authoritatively by `rootextraction.f90`

| Variable | Type / shape | `variables.f90` line | Cadence | Activity gate | Write sites |
|---|---|---|---|---|---|
| `qrot(macp)` | real(8) array (per-node) | 901 | **Instantaneous** (reset to 0 at top of `RootExtraction`, written per node) | always (when `flCropEmergence`; early return otherwise) | `rootextraction.f90:59, 88, 198, 199, 278` (Feddes / overall reduction / Jarvis compensation); `rootextraction.f90:633, 634` (JvL); `rootextraction.f90:714, 720, 725, 728, 736, 742, 747, 750` (JvL inner loop). **NO co-writers outside `rootextraction.f90`.** |
| `qrosum` | real(8) scalar | 896 | Instantaneous | always | `rootextraction.f90:61, 108, 200, 282, 616, 618, 627, 666, 756` — sole writer. (Read by waterbalance, drainage, soilhydraulics indirectly via qrot loops.) |
| `qpotrot(macp)` | real(8) array | 902 | Instantaneous (set on each node before reduction) | always | `rootextraction.f90:198` — sole writer. |
| `qredwet(macp)` | real(8) array | 903 | Instantaneous | always | `rootextraction.f90:208, 219` — sole writer. |
| `qreddry(macp)` | real(8) array | 904 | Instantaneous | always | `rootextraction.f90:209, 220` — sole writer. |
| `qredsol(macp)` | real(8) array | 905 | Instantaneous | always | `rootextraction.f90:210, 221` — sole writer. |
| `qredfrs(macp)` | real(8) array | 906 | Instantaneous | always | `rootextraction.f90:211, 222` — sole writer. |
| `qredwetsum` | real(8) scalar | 897 | Instantaneous | always | `rootextraction.f90:63, 225, 286, 292` — sole writer. |
| `qreddrysum` | real(8) scalar | 898 | Instantaneous | always | `rootextraction.f90:62, 226, 287, 293` — sole writer. |
| `qredsolsum` | real(8) scalar | 899 | Instantaneous | always | `rootextraction.f90:64, 227, 288, 294` — sole writer. |
| `qredfrssum` | real(8) scalar | 900 | Instantaneous | always | `rootextraction.f90:65, 228, 289, 295` — sole writer. |
| `Tactual` | real(8) scalar | 965 | Instantaneous (snapshot of prior `qrosum`) | `swdrought=2` only | `rootextraction.f90:60` (`Tactual = qrosum`); also reads/writes inside `JongvanLier` (lines 526, 529, 538, 557, 576). **NO external readers.** |
| `alpJvLier` | real(8) scalar | 325 | Instantaneous | `swdrought=2` only | `rootextraction.f90:628, 631` — sole writer. Read internally at line 174 (`alpdry = alpJvLier`) — pure JvL → Feddes-loop handoff. **NO external readers.** |
| `mflux(macp)` | real(8) array | 871 | Instantaneous (zeroed below noddrz; set above) | `swdrought=2` only | `rootextraction.f90:344, 383, 389` — sole writer (Task 2 evaluation through `MatricFlux`). Reads `mflux(node)` at 689, 690, 715, 716, 730, 738, 752 (within home). |
| `mroot(macp)` | real(8) array | 872 | Instantaneous | `swdrought=2` only | `rootextraction.f90:345, 711, 731, 753` — sole writer. |
| `hroot(macp)` | real(8) array | 818 | Instantaneous | `swdrought=2` only | `rootextraction.f90:346, 353, 354, 411, 678, 709, 729, 751` — primary writer. **Co-write: `cropgrowth.f90:609, 1339, 2408`** seeds `hroot(i) = h(i)` in Task 1 init (`flhydrlift` path). |
| `rootrho(macp)` | real(8) array | (near 937) | Instantaneous | `swdrought=2` only | `rootextraction.f90:347, 362, 374` — sole writer. |
| `rootphi(macp)` | real(8) array | — | Instantaneous | `swdrought=2` only | `rootextraction.f90:348, 365, 377` — sole writer. |
| `rmax(macp)` | real(8) array | 937 | Instantaneous | `swdrought=2` only | `rootextraction.f90:362, 373` — sole writer. |
| `hleaf` | real(8) scalar | 815 | Instantaneous (last-iterate cache) | `swdrought=2` only | `rootextraction.f90:339, 401, 416, 418, 445, 479, 516, 525` (JvL search); **co-write: `cropgrowth.f90:611, 1341, 2410`** (init `hleaf = -2000.d0`). |
| `Hxylem` | real(8) scalar | (near hleaf) | Instantaneous | `swdrought=2` only | `rootextraction.f90:405, 419, 446, 480, 518, 531, 559, 578, 610` — sole writer. Read by swapoutput `rot` output. |
| `mfluxtable(maho,801)` | real(8) layer×801 lookup | 870 | Init-once (Task 1 lookup) | `swdrought=2` only | `rootextraction.f90:798, 816` (initialization in `MatricFlux(1)`); read at 836, 841, 856, 861. Pure layer-config-derived; never re-written per step. |
| `flWrtNonox` | logical scalar | 1360 | Instantaneous (per-iteration override) | `swWrtNonox=1` | `rootextraction.f90:154, 156` — sole writer. **External readers (Compute): `cropgrowth.f90:689, 1742, 2104, 3001, 3387`** (`rr = 0` / `grrt = 0` when flagged). |

**Subtotal: ~22 owned fields.** Cadence: all instantaneous. Shape breakdown:
- 6 per-node arrays in primary Feddes path: `qrot`, `qpotrot`, `qredwet`, `qreddry`, `qredsol`, `qredfrs`
- 6 scalar sums / scalars in primary path: `qrosum`, `qredwetsum`, `qreddrysum`, `qredsolsum`, `qredfrssum`, `flWrtNonox`
- 6 per-node arrays in JvL path: `mflux`, `mroot`, `hroot`, `rootrho`, `rootphi`, `rmax`
- 3 JvL scalars: `Tactual`, `alpJvLier`, `hleaf`, `Hxylem` (4)
- 1 layer×801 lookup (`mfluxtable`) — JvL initialization-only

### 2b. Co-write fields — written by rootextraction but with co-writers outside

The two co-write hazards:

| Field | Where co-written outside rootextraction | Plumbing status |
|---|---|---|
| `hroot(:)` | `cropgrowth.f90:609, 1339, 2408` — `cropgrowth(1)` Task-1 init seeds `hroot(i) = h(i)` for `flhydrlift` path | `cropgrowth.f90` **uses bare `use variables`**. After migration, these three lines need either (a) state plumbing through `cropgrowth(1)` or (b) a deferred-write pattern. See Hazard #2. |
| `hleaf` | `cropgrowth.f90:611, 1341, 2410` — init `hleaf = -2000.d0` | Same `cropgrowth.f90` plumbing issue. |
| `mfluxtable` | Initialized by `MatricFlux(1, ...)` called from `cropgrowth.f90:601, 1331, 2400` — passed `dummy` for outcome, **no `state` arg** | `MatricFlux` already accepts `optional state`. Call sites work as-is for Phase 2.7 once `mfluxtable` migrates — the Task-1 init does not need state. But `MatricFlux(1)` writes the lookup as a module-level global today; it would need writing into a passed-in state. See Hazard #4. |

### 2c. Cumulatives / intermediates — NOT owned by this arc

The cumulative and intermediate root-uptake fields are **owned by soil-water-core (waterbalance.f90 + soilhydraulics.f90)**, not by `rootextraction.f90`:

| Field | Type | Reset (init) | Accumulator | Owner |
|---|---|---|---|---|
| `inqrot(macp)` | array | `soilhydraulics.f90:1097` (in `flzerointr` block) | `waterbalance.f90:420` (`inqrot(node) = inqrot(node) + qrot(node) * dt`) | soil-water-core (intermediate cohort) |
| `iqrot` | scalar | `soilhydraulics.f90:1102` (in `flzerointr`) | `waterbalance.f90:418` (`iqrot = iqrot + qrotts`) | soil-water-core |
| `iqredwet/dry/sol/frs` | 4 scalars | `soilhydraulics.f90:1104–1107` | `waterbalance.f90:428–435` | soil-water-core |
| `iqredwet_day / iqreddry_day / iqredsol_day / iqredfrs_day` | 4 scalars | `soilhydraulics.f90:1084–1087` (`flDayStart`) | `waterbalance.f90:432–435` | soil-water-core |
| `qpotrot_day(macp), qredtot_day(macp)` | 2 arrays | `soilhydraulics.f90:1090–1091` (`flDayStart`) | `waterbalance.f90:421–422` | soil-water-core |
| `iptra_day` | scalar | `soilhydraulics.f90:1088` | meteo path | soil-water-core / atmosphere |
| `cqrot` | scalar (cumulative) | `soilhydraulics.f90:1136` (`flzerocumu`) | `waterbalance.f90:490` | soil-water-core (cumulative cohort) |

**These belong to the soil-water-core arc.** This arc does NOT touch them — but **note** for Section 7 that `waterbalance.f90:418–435` is a 10-line block accumulating cumulatives from this arc's owned scalars. The cutover for `qrot/qrosum/qredXXXsum/qpotrot/qredwet/qreddry/qredsol/qredfrs` reads at those lines moves to `state%soilwater%X` — the writes (to legacy iqrot/iqredXXX) stay until soil-water-core.

### 2d. Tabulated/config inputs — read-only at run time

| Legacy global | Read at | Typed-config coverage |
|---|---|---|
| `cumdens(202)` | `rootextraction.f90:88` | covered (cropfixed/wofost/grass init populates) |
| `rdctb(22)` | `rootextraction.f90:360, 372` | covered |
| `hlim1, hlim2u, hlim2l, hlim3h, hlim3l, hlim4, adcrh, adcrl` | rootextraction lines 92–98, 122–168 | covered in cropfixed/wofost/grass configs |
| `wiltpoint, kstem, rxylem, rootradius, kroot, rootcoefa, rooteff, stephr, criterhr, taccur, dcritrtz, rdm` | various lines | all covered in cropfixed/wofost/grass + simulation_config |
| `saltmax, saltslope, salthead` | 181, 849 | covered in cropfixed/wofost/grass |
| `swdrought, swoxygen, swoxygentype, swsalinity, swfrost, swcompensate, swstressor, swWrtNonox, aeratecrit, swrootradius, swhydrlift, flhydrlift, alphacrit` | gates | all covered in crop variant configs |

**No tabular Phase 0 gaps.** Crop-side parameters are well-covered by the existing per-variant typed configs (`cropfixed_config`, `cropwofost_config`, `cropgrass_config`). Phase 0 risk for this arc is **low** (see Section 6).

---

## 3. External readers (5-category framework)

For each Section 2a owned field, files outside `src/crop/` and `src/state/` that **read** the field. Categories per playbook: (1) Output, (2) Compute, (3) Working buffer, (4) Init-seed, (5) Call-site arg.

### 3.1 Per-field reader index (primary fields)

| Field | Output (Cat 1) | Compute (Cat 2) | Working buffer (Cat 3) | Init-seed (Cat 4) | Call-site arg (Cat 5) |
|---|---|---|---|---|---|
| `qrot(:)` | `swapoutput.f90:581, 630, 651, 690, 709` (afo/aun/vap outputs); `swapoutput.f90:774, 834, 845` (.rot rotation output); `swap_csv_output.f90` (indirect via inqrot) | `soilhydraulics.f90:120, 202, 222, 244, 248, 481, 507, 527, 531, 719` (Richards F vector — sink term, 10 sites in headcalc + SoilWater(3) integrand); `solute.f90:223, 224, 225` (passive uptake `crot`); `agetracer.f90:224, 225, 265` (age tracer); `waterbalance.f90:348, 420` (mass balance + inqrot accumulator) | — | `soilgrid.f90:236, 385` (regridding writes `inqrotNew` from `inqrot`; **not from `qrot` directly** — Cat 5 false hit, but `inqrot` is downstream) | — |
| `qrosum` | — | `waterbalance.f90:339, 407` (computes `qbot = qtop + qrosum + ...`; `qrotts = qrosum * dt`); `waterbalance.f90:415, 418, 490` (tra/iqrot/cqrot accumulators) | — | — | — |
| `qpotrot(:)` | — | `waterbalance.f90:421` (`qpotrot_day = qpotrot_day + qpotrot * dt`) — sole external reader | — | — | — |
| `qredwet(:)` etc. (4 arrays) | — | `waterbalance.f90:422` (sum into `qredtot_day`) | — | — | — |
| `qredwetsum` etc. (4 scalars) | `swap_csv_output.f90:245–249` (via downstream `iqredXXX`) | `waterbalance.f90:428–435` (10-line accumulator into iqredXXX_day / iqredXXX) | — | — | — |
| `Tactual` | — | **NONE outside rootextraction** | — | — | — |
| `alpJvLier` | — | **NONE outside rootextraction** | — | — | — |
| `hroot(:)` | `swapoutput.f90:775, 833, 844, 865, 875` (.rot output) | — | — | `cropgrowth.f90:609, 1339, 2408` (Task-1 init seed `hroot(i) = h(i)`) | — |
| `mflux(:)` | `swapoutput.f90:775, 835, 846` (.rot output) | — | — | — | — |
| `mroot(:)` | `swapoutput.f90:775, 834, 845` (.rot output) | — | — | — | — |
| `rootrho(:)` | `swapoutput.f90:775, 835, 846` (.rot output) | — | — | — | — |
| `rootphi(:)` | `swapoutput.f90:775, 835, 846` (.rot output) | — | — | — | — |
| `hleaf` | `swapoutput.f90:775, 831, 842, 865, 875` (.rot output) | — | — | `cropgrowth.f90:611, 1341, 2410` (init `hleaf = -2000.d0`) | — |
| `Hxylem` | `swapoutput.f90:774, 832, 843` (.rot output) | — | — | — | — |
| `rmax(:)` | — | — | — | — | — |
| `mfluxtable(:,:)` | — | (internal to MatricFlux only) | — | `cropgrowth.f90:601, 1331, 2400` (Task-1 init `call MatricFlux(1, ...)`) | — |
| `flWrtNonox` | — | `cropgrowth.f90:689, 1742, 2104, 3001, 3387` — 5 sites gate `rr = 0` / `grrt = 0` in WOFOST/grass dynamics | — | — | — |

### 3.2 Distinct external reader files (union over all owned fields)

1. `src/soil/soilhydraulics.f90` — **heaviest compute reader.** 10 read sites of `qrot(:)` across headcalc F-vector setup (lines 120, 202, 222, 244, 248) and headcalc per-iter retry (481, 507, 527, 531) and SoilWater(3) integrand (719). Pure reader of `qrot` (no co-write). Also reads `inqrot`/`iqrot`/`iqredXXX`/`cqrot` but those are downstream cumulatives (out-of-scope).
2. `src/soil/waterbalance.f90` — `integral`/`checkmassbal`. Reads `qrot, qrosum, qpotrot, qredXXX, qredXXXsum`. The `integral` function accumulates ALL of this arc's primary scalars and arrays into the soil-water-core cohort. **Mass-balance owner.**
3. `src/io/swapoutput.f90` — `.afo`, `.aun`, `.vap`, `.rot`, `.bal` outputs. Heaviest user of the JvL arrays (`mflux/mroot/hroot/rootrho/rootphi/hleaf/Hxylem`) via the `.rot` rotation-output table; also `qrot(node)` columns in `.vap`-style outputs.
4. `src/io/swap_csv_output.f90` — TACT / TREDDRY / TREDWET / TREDSOL / TREDFRS columns via downstream `iqrot`/`iqredXXX` cumulatives. Indirect — does NOT read this arc's instantaneous fields directly.
5. `src/solute/solute.f90` — 3 read sites (lines 223–225) for passive solute uptake `crot = tscf * qrot * cml / dz`.
6. `src/solute/agetracer.f90` — 3 read sites (lines 224, 225, 265) for age-tracer uptake.
7. `src/crop/cropgrowth.f90` — **mixed**: (a) **init-seed Cat 4** writes (`hroot = h`, `hleaf = -2000`, `MatricFlux(1)`) at three call sites (cropfixed Task 1, wofost Task 1, grass Task 1); (b) **compute Cat 2** reads of `flWrtNonox` at five sites (root respiration / growth gates). Note: cropgrowth IS the home `src/crop/` tree, but it is NOT the home FILE — it co-writes/reads but lives in the same directory. For migration purposes, treat as adjacent.
8. `src/soil/soilgrid.f90` — reads `inqrot` (downstream cumulative; out-of-scope) for vertical-regrid path. NOT a reader of `qrot` directly.
9. `src/crop/management_soil.f90` — reads `inqrot` (downstream; out-of-scope).

**Distinct in-scope external reader files: 7** (soilhydraulics, waterbalance, swapoutput, swap_csv_output, solute, agetracer, cropgrowth).

### 3.3 Per-file read-site count (in-scope only)

| File | Read sites (rough) |
|---|---|
| `soilhydraulics.f90` | 10 (all `qrot`) |
| `waterbalance.f90` | ~12 (`qrot`, `qrosum`, `qpotrot`, `qredXXX×4`, `qredXXXsum×4`) |
| `swapoutput.f90` | ~14 (across qrot, hleaf, Hxylem, hroot, mflux, mroot, rootrho, rootphi in .rot + 4 sites of qrot column in .afo/.aun/.vap) |
| `solute.f90` | 3 |
| `agetracer.f90` | 3 |
| `cropgrowth.f90` | 5 (`flWrtNonox`) compute + 3+3 Cat-4 init writes |
| `swap_csv_output.f90` | 0 direct (5 indirect via iqredXXX — out-of-scope) |

**Total external read sites for in-scope fields: ~47.**

**Cat 3 hazard surface: NONE.** No mini-sim writeback for qrot/qrosum. (Boundary's mini-sim at `swapoutput.f90:3745–3820` snapshots qbot/gwl/pond/theta/h — NOT qrot.)

---

## 4. Co-writers

Files outside the home file that **write** Section 2a fields. Each is a Phase-1-dual-write or Phase-2-migration target.

| Co-writer file | Field(s) written | Context | State-in-scope today? |
|---|---|---|---|
| `src/crop/cropgrowth.f90` (Task 1 init, three variants) | `hroot(:)`, `hleaf`, `mfluxtable` (indirect via `MatricFlux(1)`) | `cropfixed_cropgrowth_init` (line 595–612), wofost equivalent (line 1325–1342), grass equivalent (line 2394–2411). Triggered when `swdrought=2`. | **NO** — `cropgrowth.f90` uses bare `use variables` and Tasks 1/2 are not state-plumbed for the JvL init path. (`cropgrowth` heat-arc Task 6 added `state%heat%tsoil` reads but **not** state writes.) See Hazard #2. |
| `src/core/initialize.f90` | `qrot, qrosum, qredXXXsum (4), hroot, mflux, alphacrit` zero-init at lines 353, 368, 401, 420–425, 804 | One-time program-start zero of legacy globals; redundant after migration (state defaults to 0). | **NO** — bare `use variables`. After migration these lines are dead and can be deleted as part of Phase 2.7 retirement. |
| `src/soil/soilhydraulics.f90` (downstream resets) | Resets `inqrot`, `iqrot`, `iqredXXX`, `cqrot`, `qpotrot_day`, `qredtot_day` at lines 1084–1136 in `SoilWater(1)` `flDayStart`/`flzerointr`/`flzerocumu` blocks | These are **downstream cumulatives** (soil-water-core) — NOT in this arc's owned set. Resets stay where they are. | YES (already plumbed). |

**Total co-writers for this arc's owned set: 1 file (cropgrowth.f90)** + the trivial initialize.f90 zero-init.

**State-plumbing gap:** `cropgrowth.f90:CropFixed/CropWofost/CropGrass` need `state` plumbing if `hroot`, `hleaf`, or `mfluxtable` migrate. This is the **central hazard** of the arc — see Hazard #2.

---

## 5. Init-order analysis

### Init-order map (from swap.f90)

```
swap.f90:182   call CalcGrid()                                ! grid dims
swap.f90:184   call soilwater_init(state%soilwater)           ! seeds 12 boundary scalars (SS-BND B-1.2)
swap.f90:186   if (flTillage) call DoTillage(1)
swap.f90:190   call SoilWater(1, state)                       ! resets cqrot/iqrot/inqrot
swap.f90:197   call drainage_init(state, config)
swap.f90:198   if (flSolute) call solute_init(state)
swap.f90:199   call heat_init(state)
…
   per-day:
   swap.f90:272   call CropGrowth(1, state%heat%tsoil)        ! Task-1 init at day start; writes hroot, hleaf, mfluxtable (when swdrought=2)
…
   per-timestep:
   swap.f90:302   call RootExtraction(state)                   ! WRITES qrot/qrosum/qredXXX
   swap.f90:305   call BoundBottom(state)
   swap.f90:312   call Drainage(state)                         ! READS qrosum indirectly via waterbalance
   swap.f90:319   call SoilWater(2, state)
                    └─ headcalc reads qrot(:) in F vector
```

**Allocation site for `qrot(:)` and the other per-node arrays:** must be inside or before `soilwater_init`. Three options:

1. **Extend `soilwater_init` to take `numnod`** and allocate `qrot(:)`, `qpotrot(:)`, `qredwet/dry/sol/frs(:)`, `mflux(:)`, `mroot(:)`, `hroot(:)`, `rootrho(:)`, `rootphi(:)`, `rmax(:)` there. **This is the boundary lesson #8 forward-compat plan; the boundary commit already mentioned "later arcs (crop-uptake) will add per-node arrays".**
2. Lazy-allocate inside `RootExtraction` on first call (allocated() guard).
3. Allocate in `cropgrowth(1)` init.

**Recommendation:** option **(1) — extend `soilwater_init`**. This is exactly the case the boundary design anticipated. The `numnod` parameter is available immediately after `CalcGrid()` at swap.f90:182.

### Init-order risk: cropgrowth(1) at swap.f90:272 reads/writes `hroot`/`hleaf`/`mfluxtable`

`CropGrowth(1, ...)` runs **after** `soilwater_init` (program init → daily loop), so by the time cropgrowth Task-1 init seeds `hroot(i) = h(i)`, the state arrays will already be allocated. **No allocated() guard required at the cropgrowth call sites** if `soilwater_init(state%soilwater, numnod)` is called at swap.f90:184.

However: the seed `hroot(i) = h(i)` must write into `state%soilwater%hroot(i)`. Either (a) plumb `state` through `CropGrowth(1, state)` (currently `CropGrowth(1, state%heat%tsoil)` — already takes a state slice) or (b) defer the seed to inside `rootextraction.f90` itself with a one-time-flag.

**Recommendation:** plumb `state` into `CropGrowth(1, state)`. The signature already takes `state%heat%tsoil`; passing the full `state` simplifies the seed AND positions the cropgrowth heat-arc gymnastics (the `dummy_X_ => tsoil` rename trick) for eventual cleanup. See Hazard #2.

### Does anyone read `qrot` BEFORE `RootExtraction` first writes it?

Per swap.f90 call order:
- `soilwater_init` zero-allocates → safe.
- `SoilWater(1, state)` at line 190 runs the `flzerointr/flzerocumu` resets but reads no per-node `qrot(:)`.
- Daily loop: `CropGrowth(1, ...)` at line 272 only touches `hroot/hleaf/mfluxtable`; does NOT read `qrot`.
- `FrozenCond(state)` at line 298 does NOT read `qrot`.
- `RootExtraction(state)` at line 302 writes `qrot(1:numnod) = 0.0d0` first thing, then computes.
- `Drainage`, `SurfaceWater`, `SoilWater(2)`/headcalc all read `qrot` AFTER it has been written.

**No pre-write reader.** A zero-allocate in `soilwater_init` is sufficient; no allocated() guard at read sites needed.

---

## 6. Config / Phase 0 candidates

The crop-side parameters (h-limits, root architecture, JvL constants, salt stress) are **already covered** by the three per-variant typed configs (`cropfixed_config`, `cropwofost_config`, `cropgrass_config`) plus `simulation_config%numerical%taccur`. Adapter `config_to_variables.f90` writes these into legacy globals via the per-variant `*_init.f90` modules.

Quick coverage check (all confirmed in `src/config/`):

| Parameter | Source config | Covered |
|---|---|---|
| `hlim1, hlim2u, hlim2l, hlim3h, hlim3l, hlim4, adcrh, adcrl` | cropfixed/wofost/grass | YES |
| `swdrought, swoxygen, swoxygentype, swsalinity, swfrost, swcompensate, swstressor, swWrtNonox, aeratecrit` | crop variant configs | YES |
| `wiltpoint, kstem, rxylem, rootradius, kroot, rootcoefa, rooteff, stephr, criterhr, taccur` | cropfixed_config + simulation_config | YES |
| `saltmax, saltslope, salthead` | cropfixed/wofost/grass | YES |
| `alphacrit, dcritrtz, rdm, swhydrlift, swrootradius` | crop variant configs (alphacrit in cropfixed/wofost/grass) | YES |
| `cumdens(202), rdctb(22), mrftb(2*magrs)` | crop variant configs (table arrays) | YES |

**Phase 0 candidates for this arc: zero new fields.** All inputs are covered. No analog to boundary's sinmax/cofqha/hplate gaps.

**Caveat:** the design phase should still grep the five TOML regression cases to confirm at least one case exercises:
- `swdrought=2` (JvL microscopic uptake — JvL fields actually written)
- `swoxygen=2` (Bartholomeus path — exercises `OxygenStress`)
- `swcompensate=1 or 2` (Jarvis/Walsum compensation branch)
- `swfrost=1` (frost stress branch)

If any of those is uncovered, add a minimal regression fixture as a B-0 task. (Looking at the existing test suite: hupselbrook covers swdrought=1, Feddes baseline. The JvL path likely needs additional coverage.) **Estimated Phase 0 work: 0–1 fixture-coverage commit.**

---

## 7. Known coupling hazards

### Hazard #1 — `cropgrowth.f90` co-writes `hroot/hleaf/mfluxtable` from Task-1 init and is NOT state-plumbed for these writes

**Scope: medium-high.** This is the analog of boundary's `tillage.f90` pond co-write hazard, but with **three call sites** (CropFixed Task 1, CropWofost Task 1, CropGrass Task 1). Each site currently writes:

```fortran
if (swdrought .eq. 2) then
   call MatricFlux(1, h(1), 1, dummy)   ! init mfluxtable
   if (swhydrlift .eq. 1) then
      flhydrlift = .true.
   else
      flhydrlift = .false.
   endif
   do i = 1,numnod
      twilt(i) = watcon(i,wiltpoint)
      hroot(i) = h(i)                    ! seed hroot
   enddo
   hleaf = -2000.d0                       ! seed hleaf
endif
```

Resolution options:
- (a) **Plumb `state` through `CropGrowth(1, state)`** (currently takes `state%heat%tsoil` slice — extend to full state). Cleanest. Three writes become `state%soilwater%hroot(i) = h(i)`, etc. Affects 3 task-1 init paths + `CropGrowth(2..)` signatures.
- (b) **Defer the seed into `rootextraction.f90` itself** behind a one-time `if (first_call)` flag. Avoids cropgrowth surgery but adds state-bookkeeping.
- (c) **Keep `hroot/hleaf/mfluxtable` as legacy globals in THIS arc**; migrate them in a later arc (e.g. crop-uptake-Phase-2 or atmosphere arc when ptra/atmdem migrate).

**Recommendation:** **(a)** plumb `state` through `CropGrowth(1, state)`. This sets up the atmosphere arc (#7) and crop-state arc (eventual) for clean threading. The boundary arc deferred `pond` to avoid `tillage.f90` surgery — for crop-uptake, the cropgrowth surgery is unavoidable because the JvL fields are crop-uptake's own. Doing it now (one well-scoped commit) is cleaner than splitting.

### Hazard #2 — Heat-arc deviation: 3 non-module cropgrowth routines used `use Variables, dummy_X_ => tsoil` rename trick

**Scope: medium.** Per heat Task 7.5, `cropgrowth.f90` contains three non-module routines (`ArableLandGerm`, `sumttd`, `grass`) that received `tsoil(:)` as a non-optional dummy arg via the rename trick. This is the **pre-existing plumbing pain** flagged in the task prompt. For crop-uptake, the analogous question: does any non-module routine in `cropgrowth.f90` read `flWrtNonox`?

Grep result: `flWrtNonox` is read at `cropgrowth.f90:689, 1742, 2104, 3001, 3387` — all inside the `CropFixed`, `CropWofost`, `CropGrass` driver routines (module routines per the heat Task 7.5 reorganization). **No non-module routine reads `flWrtNonox`.** Routine signatures should already be state-aware via Task 7.5 → Task 6 cascade.

**Verify in design phase:** grep `cropgrowth.f90` for routines that take `tsoil` via the dummy-rename trick; if any reads `flWrtNonox` or `hroot`/`hleaf`, the same rename trick may be needed. Discovery's grep suggests this is NOT the case, but verify.

### Hazard #3 — `MatricFlux(1)` initializes `mfluxtable` from outside the home file

**Scope: small.** Three `MatricFlux(1, h(1), 1, dummy)` calls from cropgrowth.f90 (lines 601, 1331, 2400) **WRITE** the layer-keyed `mfluxtable(:,:)`. After migration `mfluxtable` lives on `state%soilwater%mfluxtable`. Since `MatricFlux` already accepts `optional state`, the call sites need to pass state OR the Task-1 lookup-build must be relocated to a `soilwater_init` extension.

**Recommendation:** Relocate the `mfluxtable` build to `soilwater_init` (extended with `numnod`+layer info) so the init is unambiguous and self-contained. Then `MatricFlux(2, ...)` (used per-node at runtime) reads from `state%soilwater%mfluxtable` only. This eliminates the only Cat-4 init-seed read from outside the home tree. Alternative: keep the calls in cropgrowth(1), plumb state.

### Hazard #4 — `flWrtNonox` is a non-state-style flag with cross-subsystem semantics

**Scope: medium.** `flWrtNonox` is set by `RootExtraction` (line 154–156) and read by `cropgrowth` (5 sites: rr=0/grrt=0 gates in WOFOST and grass dynamics). It is logically a **per-timestep instantaneous boolean** carrying information from root-uptake's oxygen-stress test back to crop growth. After migration:
- The flag becomes `state%soilwater%flWrtNonox` (semantically odd — it is more crop than soil-water) OR
- It moves to a future crop-state record (which does not yet exist).

**Recommendation for THIS arc:** include `flWrtNonox` in `state%soilwater` (consistent ownership: written by RootExtraction, which is the new sole writer of all 22 fields). Document in the design ADR that a future crop-state arc may relocate it. Pragmatic: matches the "owned-and-touched" rule from boundary.

### Hazard #5 — JvL fields have NO external readers except output

**Scope: small / favorable.** `Tactual`, `alpJvLier`, `rmax(:)` have **zero external readers**. `mflux/mroot/hroot/rootrho/rootphi/hleaf/Hxylem` are read only by `swapoutput.f90` (.rot output) — pure Cat 1. Phase 2 cutover for these fields is trivial (output-side only).

### Hazard #6 — `qrosum` enters mass-balance check

**Scope: small.** `waterbalance.f90:339` computes `state%soilwater%qbot = state%soilwater%qtop + qrosum + state%surfacewater%qdrtot - QMaPo + (volact-volm1)/dt - qssdisum`. After cutover this becomes `... + state%soilwater%qrosum + ...`. Mechanical; no semantic change. **Note:** this line is **already mid-migration** — qtop/qbot have moved to state, qrosum has not. Confirms the natural Phase 2 cutover order.

### Hazard #7 — Compile-driven Phase 2.7 expectation (boundary lesson #5)

**Scope: 2–6 hidden readers expected.** Per boundary lesson #5, expect 2–6 hidden readers to surface only when dual-write is dropped (compiler reveals stale `use variables, only: qrot, ...` imports). Plan for that step. Likely surprise locations: `swapoutput.f90` use clauses (~3 lines), `waterbalance.f90` use clauses, possibly `agetracer.f90` or `solute.f90` extra reads of `qpotrot/qredXXX` for diagnostic columns.

### Hazard #8 — `Tactual = qrosum` (line 60) is an inter-iteration self-reference

**Scope: small.** Line 60 reads the **prior step's** `qrosum` into `Tactual` before zeroing `qrosum` at line 61. This is JvL-only state: prior-iter actual transpiration used as initial guess in the `flstress=TRUE` branch (line 538). After migration: `state%soilwater%Tactual = state%soilwater%qrosum` (read prior, write new). Last-writer-wins semantics carry through cleanly with the assignment order preserved.

### Hazard #9 — `MatricFlux` is `public` and exported from the module

**Scope: small.** `rootextraction_mod` exports `public :: RootExtraction, MatricFlux`. After state migration, `MatricFlux` keeps `optional intent(in) :: state` so the 3 cropgrowth Task-1 callers continue to work without state during init. The state-needing branch (`swsalinity=2` at line 846) reads `state%solute%cml` only at Task=2 evaluation, when state IS passed. The optional pattern (playbook gotcha #3) handles this correctly already.

---

## 8. Reset-cadence + cohort decision

| Cadence | Count | Cohort target |
|---|---|---|
| Instantaneous | 22 | Flat fields on `soilwater_state_t` |
| Intermediate | 0 | — (downstream `inqrot/iqredXXX*` belong to soil-water-core) |
| Cumulative | 0 | — (downstream `cqrot` belongs to soil-water-core) |
| Per-day | 0 | — |
| Init-once | 1 (`mfluxtable`) | Flat field, populated by soilwater_init or cropgrowth(1) |

**Decision: FLAT.** All 22 owned fields are instantaneous (or init-once for `mfluxtable`). No cohort sub-records introduced in this arc. Matches boundary + heat precedent.

### Recommended `soilwater_state_t` shape after THIS arc

Extending the existing 12-field boundary type:

```fortran
type :: soilwater_state_t

   ! ── Boundary arc (SS-BND, already merged) ─────────────────────────
   real(real64) :: qtop, reva, hsurf, runots, QMpLatSs
   logical      :: ftoph, FlRunoff
   real(real64) :: qbot, qbot_nonfrozen, hbot, gwlinp, deepgw

   ! ── Crop-uptake arc (THIS arc) — primary Feddes+stress path ───────
   real(real64), allocatable :: qrot(:)      ! per-node root sink (cm/d)
   real(real64), allocatable :: qpotrot(:)   ! per-node potential uptake
   real(real64), allocatable :: qredwet(:)   ! per-node wet stress contribution
   real(real64), allocatable :: qreddry(:)
   real(real64), allocatable :: qredsol(:)
   real(real64), allocatable :: qredfrs(:)
   real(real64) :: qrosum       = 0.0_real64 ! column-sum root uptake (cm/d)
   real(real64) :: qredwetsum   = 0.0_real64
   real(real64) :: qreddrysum   = 0.0_real64
   real(real64) :: qredsolsum   = 0.0_real64
   real(real64) :: qredfrssum   = 0.0_real64
   logical      :: flWrtNonox   = .false.

   ! ── JvL microscopic uptake (swdrought=2) ─────────────────────────
   real(real64), allocatable :: mflux(:)
   real(real64), allocatable :: mroot(:)
   real(real64), allocatable :: hroot(:)
   real(real64), allocatable :: rootrho(:)
   real(real64), allocatable :: rootphi(:)
   real(real64), allocatable :: rmax(:)
   real(real64), allocatable :: mfluxtable(:,:)   ! (maho, 801) — init-once lookup
   real(real64) :: Tactual    = 0.0_real64
   real(real64) :: alpJvLier  = 0.0_real64
   real(real64) :: hleaf      = 0.0_real64
   real(real64) :: Hxylem     = 0.0_real64

end type soilwater_state_t
```

**Total new fields: ~22** (6 primary per-node + 5 primary scalars + 1 flag + 6 JvL per-node + 4 JvL scalars + 1 layer×801 lookup).

### `soilwater_init` signature change

Today: `subroutine soilwater_init(sw)` — no dims arg, scalars only.

After this arc: `subroutine soilwater_init(sw, numnod)` (and possibly `nlay` for `mfluxtable`).

```fortran
subroutine soilwater_init(sw, numnod, nlay)
   type(soilwater_state_t), intent(inout) :: sw
   integer, intent(in) :: numnod, nlay
   ! ... allocate all per-node arrays to numnod ...
   ! ... allocate mfluxtable to (nlay, 801) ...
   ! ... zero scalars ...
end subroutine
```

**This is the boundary lesson #8 forward-compat plan landing exactly as anticipated** ("later arcs (crop-uptake) will add per-node arrays … that must be allocated before tillage runs").

Caller at `swap.f90:184` updates to `call soilwater_init(state%soilwater, numnod, numlay)`. `numnod` is available immediately after `CalcGrid()`.

---

## 9. Scope estimate

| Metric | Crop-uptake (THIS arc) | Boundary (ADR 0035) | Heat (ADR 0034) |
|---|---|---|---|
| Owned globals | ~22 (12 inst arrays + 8 inst scalars + 1 flag + 1 init-once table) | 12 (all inst scalars) | 12 (mixed inst scalars + 7 arrays) |
| Co-write fields | 3 (`hroot`, `hleaf`, `mfluxtable` — all from `cropgrowth(1)` init) | 1 (`kmean` top/bottom, deferred) | 1 (`rfcp` reset by soilhydraulics) |
| External reader files | 7 (soilhydraulics, waterbalance, swapoutput, swap_csv_output indirect, solute, agetracer, cropgrowth) | 14 | 15 |
| External read sites | ~47 (10 soilhydraulics + 12 waterbalance + 14 swapoutput + 6 solute/agetracer + 5 cropgrowth flWrtNonox) | ~80 | ~50 |
| Subsystems needing state plumbing | 1 (`CropGrowth(1, …)` to pass state for hroot/hleaf/mfluxtable seed) | 0 (already plumbed) | 3 (Temperature, FrozenCond, FrozenBounds) |
| Phase 0 config candidates | 0 (full coverage already) | 7–8 (sinmax/cofqha/hplate) | 6 (analytical method + tables) |
| Cohort design complexity | 0 cohorts (flat) | 0 cohorts (flat) | 0 cohorts (flat) |
| LoC home tree | 875 (rootextraction.f90 only) | 494 (boundtop+boundbottom) | 1119 |

### Suggested task decomposition (crop-uptake arc)

**Phase 0 — config gaps**

- **C-0.1** (optional): Confirm five TOML cases cover `swdrought=2` (JvL), `swoxygen=2` (Bartholomeus), `swcompensate≥1` (Jarvis), `swfrost=1`. Add minimal regression fixture if any uncovered.

(No new typed-config fields needed.)

**Phase 1 — state-type extension + dual-write**

- **C-1.1**: Extend `soilwater_state_mod` with 22 new fields (Section 8 layout). Update `soilwater_init` signature to `(sw, numnod, nlay)`. pFUnit alloc/zero test.
- **C-1.2**: Update `swap.f90:184` call to `soilwater_init(state%soilwater, numnod, numlay)`. (Confirm `numlay` is set by `CalcGrid()`; otherwise pass from config.)
- **C-1.3**: Plumb `state` into `CropGrowth(1, state)` (extending the current `state%heat%tsoil` slice). Update three Task-1 init paths in `cropgrowth.f90` to write `state%soilwater%hroot`, `state%soilwater%hleaf`. Either: relocate `MatricFlux(1)` call to `soilwater_init` (cleaner) OR pass `state` to it. Decision in design phase.
- **C-1.4**: Dual-write inside `rootextraction.f90` — every legacy write of the 22 fields also writes `state%soilwater%X`. Verify regression bit-identity.

**Phase 2 — external reader migration**

- **C-2.1**: `soilhydraulics.f90` cut over 10 read sites of `qrot(:)` to `state%soilwater%qrot(:)`.
- **C-2.2**: `waterbalance.f90` cut over reads: `qrot, qrosum, qpotrot, qredwet/dry/sol/frs, qredXXXsum*4`. ~12 sites.
- **C-2.3**: `solute.f90`, `agetracer.f90` cut over `qrot(:)` reads (3 + 3 sites).
- **C-2.4**: `cropgrowth.f90` cut over `flWrtNonox` reads (5 sites) — needs state in scope at those sites. Confirm Tasks 2/3 in cropgrowth already take state from heat Task 6.
- **C-2.5**: `swapoutput.f90` cut over output blocks: `.afo/.aun/.vap` (qrot column), `.rot` (full JvL array suite), `.bal` (cqrot reads — out-of-scope cumulatives stay legacy).
- **C-2.6**: Retire legacy globals from `variables.f90`: drop `qrot`, `qrosum`, `qpotrot`, `qredwet`, `qreddry`, `qredsol`, `qredfrs`, `qredwetsum`, `qreddrysum`, `qredsolsum`, `qredfrssum`, `Tactual`, `alpJvLier`, `mflux`, `mroot`, `hroot`, `rootrho`, `rootphi`, `rmax`, `hleaf`, `Hxylem`, `mfluxtable`, `flWrtNonox`. Drop the zero-init lines from `initialize.f90` (369, 401, 420–425, 790, 804, etc.). Expect 2–6 compile-driven fixup commits (boundary lesson #5).

**Total tasks: ~10–11** (1 Phase 0 + 4 Phase 1 + 6 Phase 2). Comparable to heat (~10 tasks); smaller than boundary (~16) because no Phase 0 work and no mini-sim hazard.

---

## 10. Open questions for design phase

1. **Cohort vs flat shape.** **Recommendation: FLAT.** All 22 owned fields are instantaneous (one is init-once). No cohort sub-records needed. Matches boundary + heat precedent.

2. **Per-node array allocation — in `soilwater_init` (with `numnod`) or lazily inside `RootExtraction`?** **Recommendation: in `soilwater_init`.** This is exactly what boundary lesson #8 anticipated ("later arcs … will add per-node arrays … that must be allocated before tillage runs"). The signature change is `soilwater_init(sw, numnod, nlay)`. Lazy-alloc would force allocated() guards everywhere — strictly worse.

3. **`CropGrowth(1, …)` state plumbing for hroot/hleaf/mfluxtable seed — do it in this arc or defer?** **Recommendation: do it in this arc.** Cropgrowth already takes a state slice (`state%heat%tsoil`); extending to full `state` is incremental. Deferring forces a deferred-seed pattern (one-time flag inside RootExtraction) which is messier than the cropgrowth surgery. The boundary arc deferred its tillage analog because tillage took no state at all; here cropgrowth already takes state.

4. **`mfluxtable` allocation site — `soilwater_init` or `cropgrowth(1)`?** **Recommendation: `soilwater_init`.** This relocates the only init-time write of `mfluxtable` out of cropgrowth and into the state lifecycle, eliminating one cropgrowth co-write entirely. Requires `numlay` parameter at soilwater_init call site. The actual lookup-build (current `MatricFlux(1)` body) moves into a private helper called from `soilwater_init` when `swdrought == 2`.

5. **`flWrtNonox` semantic owner — soil-water-state or future crop-state?** **Recommendation: soil-water-state for THIS arc.** It is written by RootExtraction (the sole writer). A future crop-state arc may relocate it; document in ADR. Pragmatic ownership rule (boundary precedent).

6. **`Tactual`, `alpJvLier`, `rmax` — keep on state or refactor to JvL-local?** All three have zero external readers; they are purely intra-RootExtraction state with cross-call lifetime. **Recommendation: keep on state** for symmetry with the rest of the JvL fields (mflux, mroot, hroot, etc.) all of which are also state. Revisit in a future JvL refactor.

7. **Phase 0 regression coverage — does any TOML case exercise `swdrought=2` JvL?** **Verify in design phase.** Quick grep of `tests/` for swdrought=2. If uncovered, add minimal fixture (small per-day case). Estimated 0–1 commit. The JvL path is large (lines 313–760) and zero coverage would be a notable risk.

8. **The 3 cropgrowth Task-1 init paths — refactor to a shared helper?** Currently CropFixed (line 595–612), CropWofost (line 1325–1342), CropGrass (line 2394–2411) each contain a copy-pasted JvL init block. **Recommendation: out of scope for this migration arc.** Note for a future cropgrowth refactor.

9. **Pure-reader optimization (boundary lesson #4) — which external files can drop legacy `use variables, only: qrot,...` imports immediately after Phase 2 cutover?** **All except `cropgrowth.f90`** (which still co-writes hroot/hleaf via Task-1 init pre-migration). `soilhydraulics.f90`, `waterbalance.f90`, `solute.f90`, `agetracer.f90`, `swapoutput.f90` are pure readers of the in-scope fields and can drop the imports as soon as their reads migrate.

10. **`Hxylem` declaration anchor in variables.f90 — separate or joined with hleaf?** Minor housekeeping question. Grep variables.f90 near line 815. Mechanical.

11. **Compile-driven Phase 2.7 expected surprise count?** Boundary surfaced 4. Crop-uptake's reader surface is smaller (47 vs 80 sites) and more concentrated (3 files dominate). **Expect 2–4 compile-driven fixup commits.**

---

## Discovery summary

- **Owned set: ~22 fields.** 12 per-node arrays + 8 scalars + 1 flag + 1 init-once table. All instantaneous (`mfluxtable` is init-once). No cohort partitioning.
- **External reader files: 7 in-scope** (soilhydraulics, waterbalance, swapoutput, solute, agetracer, cropgrowth, csv-output indirect). ~47 read sites.
- **Co-writers: 1 file** (`cropgrowth.f90` — three Task-1 init paths for hroot/hleaf/mfluxtable). State plumbing needed for `CropGrowth(1, state)`.
- **Signature status: BOTH RootExtraction, JongvanLier, JongvanLierLoop, MatricFlux, OxygenStress already take `state`.** Heat Task 6 + Solute arcs already plumbed. The boundary windfall pattern repeats: this is a field-carve + reader-cutover arc, NOT a signature-surgery arc.
- **Phase 0 gaps: ZERO new typed-config fields.** Full coverage exists in cropfixed/wofost/grass + simulation_config. Possible Phase 0 work limited to fixture coverage for swdrought=2 JvL path.
- **Layout: FLAT** (matches ADR 0034 heat + ADR 0035 boundary precedent).
- **`soilwater_init` signature CHANGE: `subroutine soilwater_init(sw, numnod, nlay)`** — per-node array allocation lands here per boundary lesson #8 forward-compat plan.
- **Hazards: 9 items.** Critical: cropgrowth Task-1 state plumbing (Hazard #1). Manageable: relocating `mfluxtable` init, `flWrtNonox` semantic ownership, compile-driven Phase 2.7. No mini-sim hazard. No `tillage`-style un-plumbed co-writer hazard.
- **Suggested task count: ~10–11** across Phase 0/1/2 (lighter than boundary's 16 due to zero Phase 0 work).
- **Decomposition: monolithic arc.** Scope sits between heat and boundary; no need to split.

End of discovery.
