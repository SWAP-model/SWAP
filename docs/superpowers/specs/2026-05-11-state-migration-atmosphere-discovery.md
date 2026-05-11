## Subsystem Migration Discovery: Atmosphere

**Date:** 2026-05-11
**Status:** discovery (read-only inventory)
**Migration #:** 7 of N — third of the four coupling-surface arcs (boundary → crop-uptake → **atmosphere** → soil-water core)
**Branch:** `development`

**Scope note:** This arc carves the **atmosphere-owned top-boundary fluxes** (`peva`, `ptra`, `empreva`, `atmdem`, etc.), the **rainfall / snow / interception state** (`grai`, `nraida`, `gsnow`, `snrai`, `ssnow`, `melt`, `subl`, `slw`, `aintc`, `aintcdt`, `sicact`, `ldwet`, `spev`, `saev`, …), and the related **per-day / intermediate / cumulative aggregators** (`pevaday`, `ptraday`, `igsnow/csubl/cmelt/csnrai/cgsnow/isnrai/isubl`) out of legacy `variables.f90` into `state%soilwater`. Resolves the atmosphere-side feeder ambiguity flagged in the soil-water mega-discovery (Section 7 “multi-owner cumulative reset pattern”).

**Predecessor docs:**
- `docs/superpowers/specs/state-migration-playbook.md` (14 lessons; boundary + crop-uptake apply most directly)
- `docs/superpowers/specs/2026-05-10-state-migration-boundary-discovery.md` (structure template)
- `docs/superpowers/specs/2026-05-10-state-migration-crop-uptake-discovery.md` (most recent precedent, ~22 owned fields)
- `docs/adr/0035-state-migration-boundary.md` (qbot/reva/qtop already migrated)
- `docs/adr/0036-state-migration-crop-uptake.md` (qrot/JvL allocatables already migrated)
- `docs/adr/0033-cumulative-reset-cohorts.md` (asymmetric reset hazard)

**Excluded from scope:** `src/atmosphere/meteoday.f90`, `src/atmosphere/meteodt.f90` (user is doing a parallel meteo refactor). Note these two files ARE the primary writers of `peva`, `ptra`, `atmdem`, AND the reset sites for `cgrai`, `cnrai`, `caintc`, `igrai`, `inrai`, `iprec`. See Hazards #1 and #2.

> Read-only discovery: no code changes. All file:line references are anchors for the design phase.

---

## 1. Big picture

### Atmosphere subsystem role

The atmosphere subsystem couples the meteorology forcing (rainfall, temperature, radiation, humidity, wind) to the soil-water column through three top-boundary flux pathways:

1. **Evapotranspiration partition** — `PenMon` produces `es0`/`et0`/`ew0` (returned via arguments, no global side effects in et.f90 itself); then meteoday.f90 turns these into `peva` (potential soil evaporation), `ptra` (potential transpiration), `atmdem` (atmospheric demand). `reduceva` reduces `peva → empreva` via Black or Boesten-Stroosnijder models.
2. **Interception** — `VonHHBraden` / `Gash` / `ruttervw` compute `aintc` (interception flux per day); `DivIntercep` partitions into `nraida` (net rain) and `nird` (net irrigation), zeroing canopy storage `sicact` for Rutter.
3. **Precipitation / snow** — `PartitionPrecipitation` partitions `grai` (in mm) into rain vs snow (`gsnow`, `snrai`, `fprecnosnow`); `snow(2)` advances `ssnow`, `slw`, `melt`, `subl` with rain-on-snow energy.

### Home files (in scope)

| File | LoC | Subroutines | State arg today? |
|---|---|---|---|
| `src/atmosphere/atmosphere_constants.f90` | 3 (stub) | (none — empty module) | n/a |
| `src/atmosphere/et.f90` | 658 | `PenMon`, `PenMon_calc`, `reduceva`, `black_reduction`, `boesten_stroosnijder_reduction` | **NO** |
| `src/atmosphere/interception.f90` | 417 | `VonHHBraden`, `Gash`, `ruttervw`, `msw1eic`, `DivIntercep` | **NO** |
| `src/atmosphere/precipitation.f90` | 134 | `PartitionPrecipitation` | **NO** (already pure, args-only) |
| `src/atmosphere/snow.f90` | 186 | `snow(task, state)` | **YES** (`optional intent(in)`, added by SS-HEAT Task 6 for `state%heat%tsoil`) |
| **In-scope total** | **1 398** | | |

Excluded (user is refactoring):
| `src/atmosphere/meteoday.f90` | 884 | `ReadMeteoDay`, `ResetMetFlx`, `ProcessMeteoDay` | NO |
| `src/atmosphere/meteodt.f90` | 453 | `MeteoDT` + `ETSine` helpers | NO |

### Entry points from outside the home tree (call graph)

| # | Call site | File:line | Phase | Task arg |
|---|---|---|---|---|
| 1 | `Snow(1, state)` | `swap.f90:216` | init | task=1 |
| 2 | `Snow(2, state)` | `swap.f90:298` | per-day at flDayStart | task=2 |
| 3 | `PenMon(...)` | `meteoday.f90:580` | per-day (only from meteoday) | n/a — pure (args) |
| 4 | `VonHHBraden(aintc)` | `meteoday.f90:521` | per-day | n/a |
| 5 | `Gash(aintc)` | `meteoday.f90:524` | per-day | n/a |
| 6 | `ruttervw(gctp, aintc, eintc)` | `meteoday.f90:649` | per-day | n/a |
| 7 | `DivIntercep(aintc)` | `meteoday.f90:530, 653` | per-day | n/a |
| 8 | `reduceva(1, nraida)` | `meteoday.f90:808` | per-day | task=1 |
| 9 | `reduceva(2, nraida)` | `meteodt.f90:358, 448` | per-timestep | task=2 |
| 10 | `PartitionPrecipitation(...)` | `meteoday.f90:369` | per-day | n/a — pure (args) |

**Critical call-flow observation:** Every entry point from outside `atmosphere/` is from `meteoday.f90` or `meteodt.f90` — both **out of scope**. The home files are dispatched THROUGH the excluded meteo orchestrator. This is the inverse of crop-uptake (where the home file `rootextraction.f90` was the entry point). The home files behave as **callee leaves** under meteoday's orchestration.

**Implication:** because meteoday/meteodt are excluded but they DO write atmosphere-owned globals (peva/ptra/atmdem/wfrac/aintc/grai/nraida/cgrai/cnrai/caintc/igrai/inrai/aintcdt and reset most cumulatives), this arc is **fundamentally split across the meteo boundary**. The atmosphere typed state can be defined and populated, but the canonical writes of peva/ptra/atmdem live in `meteoday.f90` (lines 706–865) and `meteodt.f90` (lines 348–349, 444–445). Phase 1 dual-write touches the meteo files even though we don't refactor them. See Section 4 + Hazard #1.

### Signature status

Only `snow` already takes `state`. Other entry points are called by meteoday with bare argument lists. Several routines (PenMon_calc, PartitionPrecipitation) are already pure-args and need no state; only `reduceva` reads globals (`pond`, `peva`, `nird`, `swredu`, `fldaystart`, `cofred`, `dt`, `empreva`, `ldwet`, `rsigni`, `spev`, `saev`). `reduceva` will need state for `pond` (deferred boundary D5) and `empreva` (atmosphere-owned target).

---

## 2. Owned-and-touched globals (the inventory)

### 2a. Owned globals — written authoritatively by atmosphere subsystem (home + meteo)

Owned by **atmosphere** (any writer in `atmosphere/`). Co-written by meteoday/meteodt flagged in Section 4.

#### Top-boundary fluxes (peva/ptra/atmdem/empreva)

| Variable | Type | `variables.f90` line | Cadence | Activity gate | Primary write sites |
|---|---|---|---|---|---|
| `peva` | real(8) scalar | 176 | Instantaneous (rewritten per call) | always | `meteoday.f90:706, 708, 713, 715, 722, 725, 727, 735, 737`; `meteodt.f90:349, 444`; `snow.f90:109` (zeros `peva` when `swsublim==0` and snowpack exists) |
| `pevaday` | real(8) scalar | 177 | Per-day | `swetsine=1` only | `meteodt.f90:444` (read), assigned by meteoday `ETSine` path (excluded). Read by `snow.f90:107`. |
| `ptra` | real(8) scalar | 178 | Instantaneous | always | `meteoday.f90:743, 745, 747, 751, 752, 757, 864`; `meteodt.f90:348, 445` |
| `ptraday` | real(8) scalar | 179 | Per-day | `swetsine=1` only | `meteodt.f90:445` (read); set by meteoday ETSine path |
| `empreva` | real(8) scalar | 159 | Instantaneous (writer = `reduceva` IN HOME) | `swredu in {1,2}`, no ponding | `et.f90:485, 515, 543, 569, 573, 638` (via `reduceva` → `black_reduction` / `boesten_stroosnijder_reduction`); `snow.f90:108` (zero when snowpack exists) |
| `atmdem` | real(8) scalar | 140 | Per-day | always | `meteoday.f90:816, 857, 860` |

#### Reduceva state (Black + Boesten-Stroosnijder)

| Variable | Type | `variables.f90` line | Cadence | Activity gate | Primary writes |
|---|---|---|---|---|---|
| `ldwet` | real(8) scalar | 873 | Per-day update (Black) / per-step (sub-daily) | `swredu=1` | `et.f90:497, 498, 505, 511, 639`; also `soilhydraulics.f90:873` (reset on tillage). Config seed: `read_soil_toml.f90:70`. |
| `spev` | real(8) scalar | 981 | Per-step accumulator (Boesten) | `swredu=2` | `et.f90:556, 581, 583, 640`; `soilhydraulics.f90:874` reset |
| `saev` | real(8) scalar | 976 | Per-step accumulator (Boesten) | `swredu=2` | `et.f90:558, 562, 565, 576, 641`; `soilhydraulics.f90:875` reset |

#### Precipitation / interception (per-day fluxes)

| Variable | Type | line | Cadence | Activity gate | Writers |
|---|---|---|---|---|---|
| `grai` | real(8) scalar | 164 | Per-day (mm→cm convert) | always | `meteoday.f90:295` (load from arai), `precipitation.f90:73, 117, 119` (unit conversion + detail aggregation) |
| `graidt` | real(8) scalar | 165 | Instantaneous (per-step) | `flmetdetail` mostly | meteoday/meteodt-owned (excluded) |
| `nraida` | real(8) scalar | 173 | Per-day | always | `interception.f90:403, 407, 410` (in `DivIntercep`) — sole home writer |
| `nraidt` | real(8) scalar | 174 | Instantaneous | always | meteoday/meteodt-owned |
| `aintc` | real(8) scalar | (returned, not declared) | Per-day | always | local in meteoday; argument-passed via `VonHHBraden/Gash/ruttervw` (interception.f90 lines 57, 60, 126, 128, 190) |
| `aintcdt` | real(8) scalar | 132 | Instantaneous | always | meteoday/meteodt-owned |
| `eintc` | real(8) scalar | (local to ruttervw) | Per-day | swinter=3 | `interception.f90:191` |
| `sicact` | real(8) scalar | 415 | Per-day (Rutter state) | swinter=3 | `interception.f90:189` (returned from msw1eic). External reads: `swapoutput.f90:1068, 1233` (output). |
| `fprecnosnow` | real(8) scalar | 163 | Per-day | swsnow=1 | `precipitation.f90:101, 103, 110, 129` |

#### Snow pack state (snow.f90)

| Variable | Type | line | Cadence | Activity gate | Writers |
|---|---|---|---|---|---|
| `ssnow` | real(8) scalar | 1079 | Instantaneous (carry-state) | `flSnow`/`swsnow=1` | `snow.f90:74, 76, 122, 141, 156, 165`; `precipitation.f90:127` (zero when swmetdetail=1); seeded from `config_to_variables.f90:561`. Heavy external reader use (heat, output, mass-balance). |
| `snowinco` | real(8) scalar | 1077 | Once-per-cumu-reset (snapshot at zerocumu) | `flSnow` | `snow.f90:74, 99` |
| `slw` | real(8) scalar | (near 1078) | Instantaneous (carry) | `flSnow` | `snow.f90:141, 144, 153, 156, 166` |
| `melt` | real(8) scalar | 1074 | Instantaneous | `flSnow` | `snow.f90:123, 138, 157, 163`. External readers: `boundtop.f90:137, 180`, `soilhydraulics.f90:109, 639, 649, 651` (top-flux equation). |
| `subl` | real(8) scalar | 1080 | Instantaneous | `flSnow`, swsublim=0 | `snow.f90:103, 106, 124, 164` |
| `gsnow` | real(8) scalar | 1070 | Per-day | `swsnow=1` | `precipitation.f90:79, 83, 88, 108, 126`; read by `snow.f90:121, 123, 141, 171`. External: `interception.f90:403` (DivIntercep). |
| `snrai` | real(8) scalar | 1078 | Per-day | `swsnow=1` | `precipitation.f90:93, 95, 109, 128`; read by `snow.f90:131, 132, 144, 173, 177`. External: `interception.f90:403`. |
| `igsnow` | real(8) scalar | 1071 | Intermediate (period accumulator) | `flSnow` | `snow.f90:87, 171` — reset on `flzerointr`, accumulate per step. |
| `isubl` | real(8) scalar | 1073 | Intermediate | `flSnow` | `snow.f90:88, 172` |
| `isnrai` | real(8) scalar | 1072 | Intermediate | `flSnow` | `snow.f90:89, 173` |
| `cgsnow` | real(8) scalar | 1066 | Cumulative | `flSnow` | `snow.f90:95, 174` — reset on `flzerocumu`. |
| `csubl` | real(8) scalar | 1069 | Cumulative | `flSnow` | `snow.f90:96, 175` |
| `cmelt` | real(8) scalar | 1067 | Cumulative | `flSnow` | `snow.f90:98, 176` |
| `csnrai` | real(8) scalar | 1068 | Cumulative | `flSnow` | `snow.f90:97, 177` |

#### Atmosphere cumulatives & intermediates — multi-owner (the hazard cluster)

| Variable | Type | line | Cadence | Reset site | Accumulate site | Owner question |
|---|---|---|---|---|---|---|
| `cgrai` | scalar | 147 | Cumulative | `meteoday.f90:420` (`flzerocumu`) — excluded | `waterbalance.f90:512` | Reset stays legacy until meteoday refactor; this arc owns the cumulative on state. **Asymmetric reset hazard.** |
| `cnrai` | scalar | 148 | Cumulative | `meteoday.f90:421` | `waterbalance.f90:513` | Same |
| `caintc` | scalar | 145 | Cumulative | `meteoday.f90:422` | `waterbalance.f90:510` | Same |
| `cpeva` | scalar | 149 | Cumulative | `soilhydraulics.f90:1141` | `waterbalance.f90:500` | Reset by soil; classify atmosphere cumulative |
| `cptra` | scalar | 150 | Cumulative | `soilhydraulics.f90:1140` | `waterbalance.f90:499` | Same |
| `cevap` | scalar | 146 | Cumulative | `soilhydraulics.f90:1142` | `waterbalance.f90:501` | actual soil evap — atmosphere-owned conceptually but reset by soil-water. Boundary of arc. |
| `ipeva` | scalar | 169 | Intermediate | `soilhydraulics.f90:1113` | `waterbalance.f90:470` | atmosphere field |
| `iptra` | scalar | 170 | Intermediate | `soilhydraulics.f90:1112` | `waterbalance.f90:469` | atmosphere field |
| `ievap` | scalar | 167 | Intermediate | `soilhydraulics.f90:1114` | n/a (waterbalance) | actual evap (qtopd-like) |
| `igrai` | scalar | 827 | Intermediate | `waterbalance.f90:388` + `meteoday.f90:414` (DUAL reset — hazard) | `waterbalance.f90:476` | **double-reset** — both files clear it at flzerointr. Asymmetric pattern. |
| `inrai` | scalar | 168 | Intermediate | `waterbalance.f90:389` + `meteoday.f90:415` | `waterbalance.f90:478` | Same — double reset. |

> The double-reset of `igrai`/`inrai` (by both `waterbalance.f90` and `meteoday.f90.ResetMetFlx`) is real, currently safe by ordering but exactly the asymmetric-reset hazard from ADR 0033. Phase 2 must collapse to a single owner.

#### Wet-fraction / intermediate canopy

`wfrac` is local-only in meteoday (lines 662–695) — not declared in variables.f90, so excluded.

**Subtotal: ~33 owned atmosphere fields.**
Breakdown:
- 5 top-boundary fluxes (peva, pevaday, ptra, ptraday, atmdem) + 1 reduced (empreva) = 6
- 3 reduceva state (ldwet, spev, saev) = 3
- 6 precipitation/interception per-day (grai, graidt, nraida, nraidt, aintcdt, fprecnosnow) + 2 canopy (sicact, eintc) = 8 [aintc is not declared in variables — purely local/returned]
- 8 snow flux/state (ssnow, snowinco, slw, melt, subl, gsnow, snrai, + igsnow/isubl/isnrai/cgsnow/csubl/cmelt/csnrai = 7 aggregators) = 15

Note: `eintc` and `aintc` are not in variables.f90 — they pass between meteoday and interception.f90 as local args. They are NOT owned state.

### 2b. Co-write fields — atmosphere writes but ownership is elsewhere

| Field | Atmosphere writes | Real owner | Status |
|---|---|---|---|
| `nird` | `interception.f90:404, 408, 411` (DivIntercep) | irrigation subsystem (`irrigation.f90` resets cgird/cnird; `tillage.f90` reads `nraida`) | Atmosphere co-writes by partitioning. Likely defer or share. |
| `nird` reset | `soilhydraulics.f90:871` (after tillage) | mixed | Co-write tracker |
| `pond` (read only in et.f90:612, 637) | (deferred boundary — D5) | boundary (deferred) | Not written, but **read-through-deferred-global** in reduceva. **Hazard #5.** |

`nird` ownership lies more naturally with crop/irrigation (which already writes `cgird/cnird/igird/inird`). For this arc: **defer** `nird` ownership; atmosphere just co-writes it via DivIntercep partition.

### 2c. Cumulatives / intermediates — the cohort question

| Cadence | Fields | Reset gate | Reset writer |
|---|---|---|---|
| Instantaneous | peva, ptra, empreva, melt, subl, slw, ssnow, snowinco, graidt, nraidt, aintcdt | (overwritten each call) | atmosphere |
| Per-day | grai, nraida, atmdem, pevaday, ptraday, gsnow, snrai, fprecnosnow, sicact, aintc-local | (overwritten each meteoday call) | meteoday |
| Intermediate | igrai, inrai, ipeva, iptra, ievap, igsnow, isubl, isnrai, ldwet (per-step) | `flzerointr` | mixed (snow.f90, soilhydraulics.f90, waterbalance.f90, meteoday.f90) |
| Cumulative | cgrai, cnrai, caintc, cpeva, cptra, cevap, cgsnow, csubl, csnrai, cmelt, (cgird/cnird are irrigation-owned) | `flzerocumu` | mixed (snow.f90, soilhydraulics.f90, meteoday.f90) |
| Per-event | spev, saev (B-S Stroosnijder) | reset on tillage | et.f90 + soilhydraulics |

ADR 0033 cohort precedent suggests packaging the 10–11 cumulative scalars under a **`atmosphere_cumulative_t` sub-record** with `flzerocumu` reset semantics, and similarly the 8–9 intermediate scalars under `atmosphere_intermediate_t` with `flzerointr` reset. See Section 8 for the recommendation.

### 2d. Tabulated / config inputs — read-only at runtime

| Legacy global | Source config | Coverage |
|---|---|---|
| `swsnow, snowcoef, teprrain, teprsnow` | `meteorology_config.snow` | YES |
| `swredu, swcfbs, cfbs, cofredbl, cofredbo, rsigni, cfevappond` | `meteorology_config.evaporation` | YES (`cofred` filled by adapter from `cofredbl/bo`) |
| `swetr, swdivide, swmetdetail, swrain, swetsine, swinter` | `meteorology_config.meteo` | YES |
| `lat, alt, altw, angstroma, angstromb` | `meteorology_config.meteo` | YES |
| `cofab, kdif, kdir, swcf, ch, cftb, cfeic, cfeictb` | crop variant configs | YES |
| `fimin, siccapact (initial), siccaplai, siccaptb, avevaptb, avprectb, pfreetb, pstemtb, scanopytb` | crop variant configs | YES |
| `swsublim` | (search needed — only used in snow.f90:104) | partial (declared in variables, may be silent-default) — **Phase 0 candidate** |
| `ssnow (initial), ldwet (initial)` | `soil_config.initial` | YES |

**Phase 0 candidates (~1 field):** verify `swsublim` is config-driven and has no silent default; otherwise it's covered. Phase 0 work is likely **0–1 commit**.

---

## 3. External readers (5-category framework)

For each Section 2a owned field, files outside `src/atmosphere/` and `src/state/` that **read** the field. Categories: (1) Output, (2) Compute, (3) Working buffer, (4) Init-seed, (5) Call-site arg.

### 3.1 Per-field reader index (primary fluxes)

| Field | External reader files (sites) | Category |
|---|---|---|
| `peva` | `boundtop.f90:127`, `snow.f90:106` (in-home), [meteoday/meteodt out-of-scope co-writers] | Cat 2 (compute: reva clamp) |
| `ptra` | `cropgrowth.f90:669, 672, 715, 719, 722, 1712, 1715, 3004, 3007` (5 distinct read sites in 3 variants); `rootextraction.f90:73, 101, 105, 107, 108, 257, 258, 297, 298, 435, 440, 456, 460, 469, 481, 486, 517, 520, 548, 654, 656, 660, 664, 764, 786` (25 sites); `waterbalance.f90:396, 442` (mass balance); `swapoutput.f90:155` (output gate) — **~35 sites across 4 files** | Cat 2 + Cat 1 |
| `empreva` | `boundtop.f90:129` (alternative reva clamp when `swredu>0`) | Cat 2 |
| `atmdem` | `rootextraction.f90:105, 107, 108` (hlim3 dynamic switch) | Cat 2 |
| `ldwet` | `soilhydraulics.f90:873` (tillage reset); `read_soil_toml.f90:70`; `soil_config.f90:74, 428, 429`; `config_to_variables.f90:565` | Cat 4 (init-seed) |
| `spev` | `soilhydraulics.f90:874` (tillage reset) | Cat 4 |
| `saev` | `soilhydraulics.f90:875` (tillage reset) | Cat 4 |
| `grai` | (only meteoday/meteodt read it; no in-scope external readers) | none |
| `nraida` | `tillage.f90:177, 178, 181, 370, 371` (tillage Bdens compute); written + read by meteoday | Cat 2 |
| `aintcdt` | `waterbalance.f90:467, 510` (iintc/caintc accumulator) | Cat 2 |
| `sicact` | `swapoutput.f90:1068, 1233` (output) | Cat 1 |
| `ssnow` | **heavy** — 36 sites total. Highlights: `temperature.f90:165, 170` (heat conductivity through snow), `swapoutput.f90:214, 326, 382, 456, 477` (output), `swap_csv_output.f90:12, 256, 385, 470, 477` (CSV), `boundtop.f90:?`, `config_to_variables.f90:561, 569` (seed), `read_soil_toml.f90:64`, `soil_config.f90:71, 422, 423` (config) | Cat 1 + Cat 2 + Cat 4 |
| `melt` | `boundtop.f90:137, 180`, `soilhydraulics.f90:109, 639, 649, 651` (top-flux equation — net inflow includes melt) | Cat 2 |
| `gsnow` | `interception.f90:403` (DivIntercep, in-home) | Cat 2 in-home |
| `snrai` | `interception.f90:403` (DivIntercep, in-home) | Cat 2 in-home |
| `pevaday` | `snow.f90:107` (in-home) | none external |
| `ptraday` | (only meteodt) | none |
| `igsnow`, `isubl`, `isnrai`, `cgsnow`, `csubl`, `csnrai`, `cmelt` | `swapoutput.f90:381` (snow output block; csubl maps to inputs to ies0/iet0 cumulatives also) | Cat 1 |
| `cgrai`, `cnrai`, `caintc`, `cpeva`, `cptra`, `cevap` | `swapoutput.f90:212–344` heavy (.bal output composite block); `swap_csv_output.f90` indirect | Cat 1 |
| `ipeva`, `iptra`, `ievap`, `igrai`, `inrai` | `swap_csv_output.f90:8, 12, 237, 244` (EPOT/TPOT columns); `swapoutput.f90:381–558` (.inc, .csv) | Cat 1 |

### 3.2 Distinct external reader files

1. `src/crop/rootextraction.f90` — heavy `ptra`, `atmdem` reader (~28 sites total). Already takes `state` per ADR 0036. Pure reader.
2. `src/crop/cropgrowth.f90` — `ptra` reader (~9 sites across CropFixed/CropWofost/CropGrass). Takes `state` since SS-CRP C-1.3.
3. `src/crop/tillage.f90` — `nraida` reader (Bdens compute, 5 sites). Currently uses bare `use variables`.
4. `src/boundary/boundtop.f90` — `peva`, `empreva`, `melt` readers. Already state-plumbed (ADR 0035).
5. `src/soil/soilhydraulics.f90` — `melt` (4 sites top-flux), `ldwet`/`spev`/`saev` (3 reset sites), `nird` (zero), `cpeva`/`cptra`/`cevap` resets, `ipeva`/`iptra`/`ievap` resets. State-aware.
6. `src/soil/waterbalance.f90` — accumulates `ipeva`, `iptra`, `igrai`, `inrai`, `cpeva`, `cptra`, `cevap`, `caintc`, `cgrai`, `cnrai` (10+ accumulator sites). State-aware.
7. `src/heat/temperature.f90` — `ssnow` reader (heat conductivity through snowpack).
8. `src/io/swapoutput.f90` — heaviest output reader: all cumulatives, intermediates, `ssnow`, `snowinco`, `sicact`. ~30 sites.
9. `src/io/swap_csv_output.f90` — EPOT/TPOT/SSNOW columns and dstor computation. ~6 sites.
10. `src/io/macroporeoutput.f90` — commented-out `ssnow` reference (dormant).
11. `src/solute/solute.f90`, `src/solute/agetracer.f90` — `nird` reads (irrigation co-write; out-of-scope for THIS arc).
12. `src/io/toml/config_to_variables.f90` — `ssnow`, `ldwet` initial-seed.
13. `src/io/toml/read_soil_toml.f90` — same.
14. `src/config/soil_config.f90` — `ssnow`, `ldwet` typed-config holders.
15. `src/config/meteorology_config.f90` — many config fields; not runtime readers.

**Distinct in-scope external reader files: ~9–10** (rootextraction, cropgrowth, tillage, boundtop, soilhydraulics, waterbalance, temperature, swapoutput, swap_csv_output).

### 3.3 Per-file read-site estimate (in-scope only)

| File | Read sites |
|---|---|
| `rootextraction.f90` | ~28 (`ptra` × 25 + `atmdem` × 3) |
| `cropgrowth.f90` | ~9 (`ptra` only) |
| `waterbalance.f90` | ~14 (10 accumulator reads + 4 cumulative reads) |
| `swapoutput.f90` | ~30 (cumulatives + intermediates + ssnow + sicact) |
| `swap_csv_output.f90` | ~6 |
| `soilhydraulics.f90` | ~8 (melt × 4 + ldwet/spev/saev resets + cumulatives resets) |
| `boundtop.f90` | 3 (peva, empreva, melt) |
| `temperature.f90` | 2 (ssnow) |
| `tillage.f90` | 5 (nraida) |

**Total external read sites: ~105.** Larger than boundary (80) and crop-uptake (47). Dominated by the cumulative/intermediate output block.

### 3.4 Cat 3 (working-buffer / mini-sim) hazard

None observed. The legacy `swapoutput.f90` mini-sim writeback (lines ~3745+) does NOT snapshot atmosphere fluxes. Confirmed clean.

---

## 4. Co-writers

Files outside the home tree that **write** Section 2a owned fields.

| Co-writer file | Fields written | Plumbing status | State-in-scope today? |
|---|---|---|---|
| **`src/atmosphere/meteoday.f90`** (EXCLUDED) | `peva`, `ptra`, `atmdem`, `grai` (load), `wfrac` (local), `aintc` (local), `cgrai/cnrai/caintc` resets, `igrai/inrai` resets, `pevaday/ptraday` writes (ETSine), `aintcdt` | bare `use variables` | NO — and not in scope to refactor |
| **`src/atmosphere/meteodt.f90`** (EXCLUDED) | `peva`, `ptra`, `nraida` (via reduceva) | bare `use variables` | NO |
| `src/soil/soilhydraulics.f90` | `cpeva, cptra, cevap` (resets at 1140–42), `ipeva, iptra, ievap` (resets at 1112–14), `ldwet, spev, saev` (resets at 873–75 after tillage), `nird = 0` (line 871) | state-aware | YES |
| `src/soil/waterbalance.f90` | `igrai, inrai` (resets 388–89) — DOUBLE-RESETTER with meteoday | state-aware | YES |
| `src/crop/irrigation.f90` | `cgird, cnird` resets (94–95); `gird, nird` writes per event | bare `use variables` | partial — `irrigation(2, state)` is state-plumbed |
| `src/core/initialize.f90` | one-time zeros of cgrai, cnrai, caintc, cpeva, cptra, cevap, cgird, cnird, ipeva, iptra, ievap, inrai, igrai, igsnow, isnrai, isubl, csnrai, csubl, cgsnow, cmelt, aintcdt | bare `use variables` | NO — redundant after migration |

**Co-writer summary:**
- **2 EXCLUDED files (meteoday, meteodt) write the canonical peva/ptra/atmdem.** This is THE central plumbing problem of this arc — see Hazard #1.
- **1 unavoidable double-reset** (`igrai/inrai` reset by both waterbalance AND meteoday) — must collapse to single owner during Phase 2.
- **soilhydraulics** owns multiple atmosphere-cumulative resets (cpeva/cptra/cevap/ipeva/iptra/ievap/ldwet/spev/saev). Should the resets move to atmosphere state owner? See Section 8.
- **irrigation** co-writes `cgird/cnird` (irrigation-owned, not atmosphere — defer).

**Total atmosphere-owned co-writers: 5 files** (or 3 if we exclude initialize.f90 zeros and the irrigation case).

---

## 5. Init-order analysis

### Init-order map (from swap.f90, lines 184–306)

```
swap.f90:184  call TimeControl(1)
swap.f90:187  call CalcGrid()                            ! grid dims (numnod, numlay)
swap.f90:188  call soilwater_init(state%soilwater, numnod, numlay)   ! seeds 12+22 boundary/crop-uptake fields
swap.f90:190  if (flTillage) call DoTillage(1)
swap.f90:194  call SoilWater(1, state)                   ! resets cpeva/cptra/cevap/ipeva/iptra/ievap (atmosphere fields!)
swap.f90:201  call drainage_init(state, config)
swap.f90:203  call heat_init(state)
swap.f90:216  if (flSnow) call Snow(1, state)            ! ssnow init from snowinco

   per-day:
swap.f90:273    call ReadMeteoDay()                       ! reads arai/atav/atmn/atmx/ahum/awin → grai/tav/tmn/tmx/hum/win/rad
swap.f90:276    call CropGrowth(1, state%heat%tsoil, state)
swap.f90:282    if (flIrrigate) call irrigation(2, state) ! gird/nird events
swap.f90:285    call ProcessMeteoDay()                    ! [meteoday — out-of-scope] writes peva/ptra/atmdem/cgrai/cnrai/caintc resets
swap.f90:291    if (flMeteoDt .or. flETSine) call MeteoDT() ! [meteodt — out-of-scope]
swap.f90:298    if (flSnow .and. flDayStart) call Snow(2, state)  ! writes ssnow/melt/subl/slw etc.

   per-timestep:
swap.f90:306    call RootExtraction(state)                ! reads ptra, atmdem
swap.f90:309    call BoundBottom(state)
swap.f90:312    ...
swap.f90:316    Drainage(state)
swap.f90:323    SoilWater(2, state) (headcalc reads melt indirectly via boundtop)
```

**Key init-order observations:**

1. **`SoilWater(1, state)` at line 194 resets atmosphere cumulatives (cpeva/cptra/cevap/ipeva/iptra/ievap) BEFORE Snow(1) at line 216.** Currently safe because Snow(1) only seeds `ssnow` from `snowinco` and does not touch cumulatives. But the reset-ownership is wrong — atmosphere cumulatives should be reset by an atmosphere routine, not by `SoilWater(1)`.

2. **`Snow(1, state)` runs at line 216, AFTER soilwater_init but BEFORE any meteo data is loaded.** It reads `swinco` to decide between `snowinco = ssnow` or `ssnow = snowinco`. The state init (atmosphere init) should be added at line 188 alongside `soilwater_init`, OR a separate `atmosphere_init(state%atmosphere, numnod)` call inserted between 188 and 190. Per-node arrays are NOT used by atmosphere — all fields are scalars. **No allocatables needed.**

3. **`ReadMeteoDay` at line 273 writes `grai`, `tav`, etc.** Currently in meteoday (excluded). The atmosphere state extension can be populated by `meteoday` via dual-write at the Phase 1 stage; the state owns the fields, but meteoday continues to write the legacy globals until the meteo refactor.

4. **`Snow(2)` runs after `ProcessMeteoDay` and `MeteoDT`** (and at flDayStart only). `peva`, `gsnow`, `snrai`, `tav` are read inside `Snow(2)`. Reading order: `peva` is written by meteoday before snow(2) — OK.

5. **`RootExtraction(state)` reads `ptra`/`atmdem`** — these are written by meteoday at lines 743–747/816. Read-after-write order is correct in legacy code; migrated state must follow same order.

### Pre-write readers?

- `Snow(1)` reads `swinco` and `snowinco` (config-seeded) — OK.
- `Snow(2)` reads `tsoil(1)` (via state%heat) → OK now.
- `peva` is written in meteoday for the day, then snow(2) may zero it (sublimation path). The subsequent `RootExtraction` reads `ptra` only (already finalized). Order is clean.

**No pre-write reader gap.** Init-zero of atmosphere scalars in `atmosphere_init` is sufficient.

### Allocation needs

All atmosphere-owned fields are **scalars**. No per-node arrays. No allocation site needed beyond simple `=0.0d0` in the init routine.

**Recommendation:** introduce `atmosphere_state_mod` with `atmosphere_init(state%atmosphere)` (no `numnod`/`nlay` parameter). Call at swap.f90:188 immediately after `soilwater_init`.

---

## 6. Config / Phase 0 candidates

Quick coverage check of typed configs:

| Parameter | Config home | Coverage |
|---|---|---|
| `swsnow, snowcoef, teprrain, teprsnow` | meteorology.snow | YES |
| `swredu, swcfbs, cfbs, cofredbl, cofredbo, rsigni, cfevappond` | meteorology.evaporation | YES |
| `swetr, swdivide, swmetdetail, swrain, swetsine, swinter, swmetfilall, nmetdetail` | meteorology.meteo | YES |
| `lat, alt, altw, angstroma, angstromb, rsoil, rsc, rsw` | meteorology.meteo | YES |
| `cofab, kdif, kdir, swcf, ch, cftb, cfeic` | crop variant configs | YES |
| `fimin, siccapact (init), siccaplai, siccaptb` | crop variant configs | YES |
| `swsublim` | (not searched in detail) | **PHASE 0 CANDIDATE — verify** |
| `ssnow (initial), ldwet (initial)` | `soil_config.initial` | YES |
| `spev (initial), saev (initial)` | (search; if absent → silent default) | **PHASE 0 CANDIDATE** |
| `swhydrlift, swcfbs etc.` | covered | YES |

**Phase 0 candidates: ~2 fields** (swsublim, spev/saev initial). Both are minor.

---

## 7. Known coupling hazards

### Hazard #1 — meteoday/meteodt (EXCLUDED) are the canonical writers of `peva`/`ptra`/`atmdem`

**Scope: HIGH.** The atmosphere subsystem's primary outputs (`peva`, `ptra`, `atmdem`) are written by `meteoday.f90:706–865` and `meteodt.f90:348–445`, NOT by any home file. The home files only write subordinate fields (`empreva` via reduceva from et.f90; `ssnow/melt/subl` from snow.f90; `grai/gsnow/snrai` from precipitation.f90 + interception.f90).

User explicitly excludes meteoday/meteodt from this arc (parallel refactor). Two resolution strategies:

- **(a) Dual-write inside meteoday/meteodt anyway** — touch the excluded files minimally to add `state%atmosphere%peva = peva` etc. alongside legacy writes. This is a 5–10 line patch per file; does NOT collide with a meteo refactor since the dual-write is additive. **RECOMMENDED.**
- **(b) Defer peva/ptra/atmdem to the meteo refactor** — only carve out the home-file-owned fields (empreva, ssnow/melt/subl/slw, grai/gsnow/snrai/aintc-side state, ldwet/spev/saev). Reduces this arc to maybe 15 fields. Cleaner subsystem cut, but leaves the highest-traffic readers (ptra in rootextraction × 25 sites) unmigrated. The meteo refactor would then own the second carve.

**Recommendation: (a)** — touch meteoday/meteodt with additive dual-writes. The meteo refactor remains free to restructure the orchestrator without invalidating the state writes.

### Hazard #2 — Multi-owner cumulative reset pattern (ADR 0033 asymmetric-reset hazard)

**Scope: HIGH.** Atmosphere cumulatives are accumulated by `waterbalance.f90` but **reset by three different files**:
- `meteoday.f90:420–422` resets `cgrai`, `cnrai`, `caintc` under `flzerocumu` (excluded)
- `soilhydraulics.f90:1140–42` resets `cpeva`, `cptra`, `cevap` under `flzerocumu`
- `snow.f90:95–98` resets `cgsnow`, `csubl`, `csnrai`, `cmelt` under `flzerocumu` (in-home)
- `irrigation.f90:94–95` resets `cgird`, `cnird` (out-of-arc — irrigation domain)

Per ADR 0033, the cohort packaging convention is: **fields with the same reset cadence + same reset gate go into the same sub-record with a type-bound reset method.** Atmosphere cumulatives split across **3 reset sites** — these need to be unified into a single `atmosphere_cumulative_t%reset()` call. The meteoday reset (cgrai/cnrai/caintc) would move out of meteoday into an atmosphere call invoked from `SoilWater(1)` or a new `atmosphere_reset(state%atmosphere, flzerocumu, flzerointr)` invocation. This contracts meteoday on the way out (helpful for the meteo refactor).

**Same for intermediates (igrai/inrai/ipeva/iptra/ievap/igsnow/isubl/isnrai)** — 3 reset sites (meteoday + soilhydraulics + snow + waterbalance for the igrai/inrai double-reset).

### Hazard #3 — `igrai`/`inrai` double-reset

**Scope: MEDIUM.** `igrai` and `inrai` are reset at BOTH:
- `meteoday.f90:414, 415` (in `ResetMetFlx`, under `flzerointr`)
- `waterbalance.f90:388, 389` (in `integral`, under `flzerointr`)

Currently safe (both write 0.0 under the same gate), but per ADR 0033 this is exactly the asymmetric-reset pattern that breaks under cohort migration. Phase 2 must elect a single owner. Recommendation: **atmosphere** (since `grai` itself is atmosphere-owned and `nraida` is the partition output of DivIntercep).

### Hazard #4 — `nird`/`gird` ownership ambiguity (irrigation vs atmosphere)

**Scope: SMALL.** `interception.f90:DivIntercep` writes `nird` (net irrigation after interception). But `gird` is set by `irrigation.f90`. The irrigation subsystem owns the canonical write. For this arc: **defer** `nird`/`gird`; treat them as irrigation-owned passthroughs. Atmosphere reads them in DivIntercep — no state write needed in this arc.

### Hazard #5 — `pond` read in et.f90 (boundary-deferred field)

**Scope: SMALL.** `et.f90:reduceva` reads `pond` at line 637 (`if (pond > POND_THRESHOLD)`). `pond` is **deferred** per boundary arc (D5 — boundary discovery flagged it as out-of-scope, stays legacy until soil-water-core arc). Resolution: `reduceva` takes `state` as `optional intent(in)`; reads `pond` from state when available, falls back to legacy global. After soil-water-core migrates `pond` to `state%soilwater%pond`, this becomes `state%soilwater%pond`. For THIS arc: **leave `pond` as legacy global**; just plumb `state` into reduceva to enable future cutover.

### Hazard #6 — `reduceva` globals soup

**Scope: MEDIUM.** `reduceva` reads/writes `swredu, fldaystart, cofred, dt, empreva, ldwet, nird, peva, pond, rsigni, spev, saev` — 12 globals. Of these, this arc owns `empreva, ldwet, peva, spev, saev` (5). Others come from time control (`dt, fldaystart`), config (`swredu, cofred, rsigni`), irrigation (`nird`), or are deferred (`pond`). Adding `state` argument and migrating the 5 owned reads/writes is straightforward.

### Hazard #7 — Snow's read of `peva` and write-back zero

**Scope: SMALL.** `snow.f90:106–109` reads `peva`, writes `subl = peva`, then `peva = 0` and `empreva = 0` to suppress further soil evaporation when snow pack exists. After migration, `state%atmosphere%peva` is both READ and WRITTEN by snow.f90. Already takes `state` (optional intent(in)) → must promote to `intent(inout)` for the writes. This breaks the snow signature contract and ripples to swap.f90:216, 298. Minor but worth noting.

### Hazard #8 — `temperature.f90` reads `ssnow` for heat conductivity

**Scope: SMALL.** `temperature.f90:165, 170` reads `ssnow` to compute snowpack thickness in the heat solver. Pure reader; no co-write. State `temperature` already takes `state` (post-heat-arc). Just retarget the read to `state%atmosphere%ssnow`.

### Hazard #9 — Compile-driven Phase 2.7 expectation

**Scope: 3–6 hidden readers expected.** Per playbook lesson #5, the dual-write drop step will reveal stale `use variables, only: peva, ptra, ssnow, ...` imports. Likely surprise sites: `swapoutput.f90` use clauses, `swap_csv_output.f90` use clauses, `macroporeoutput.f90` (commented out but may have live imports), `tillage.f90` (`nraida`). Estimate 3–6 fixup commits.

### Hazard #10 — heat-arc finding: meteoday reads `pond` for interception path

**Scope: VERIFY.** ADR 0034 (heat) flagged that meteoday's interception path reads `pond` somewhere. Quick check: `grep -n "\bpond\b" src/atmosphere/meteoday.f90` shows pond reads at meteoday — these are out of scope (excluded subsystem), but worth noting: when the meteo refactor happens, those `pond` reads will need state-aware routing too. Not this arc.

### Hazard #11 — `ssnow = 0.0d0` in precipitation.f90:127 (silent state mutation)

**Scope: SMALL.** `precipitation.f90:127` mutates `ssnow = 0.0d0` when `swmetdetail = 1` (commented "Note: This modifies a state variable - consider refactoring"). Co-write hazard — precipitation.f90 is in-home but it's overwriting a snow-state variable. After migration this becomes `state%atmosphere%ssnow = 0.0d0` — semantically owned by atmosphere already, just untidy. Document but don't fix structurally.

---

## 8. Reset-cadence + cohort decision

Field classification by cadence (33 fields total):

| Cadence | Count | Fields |
|---|---|---|
| Instantaneous (per-call) | 11 | peva, ptra, empreva, melt, subl, slw, ssnow, snowinco, graidt, nraidt, aintcdt |
| Per-day | 9 | grai, nraida, atmdem, pevaday, ptraday, gsnow, snrai, fprecnosnow, sicact |
| Intermediate (flzerointr-reset) | 8 | igrai, inrai, ipeva, iptra, ievap, igsnow, isubl, isnrai |
| Cumulative (flzerocumu-reset) | 10 | cgrai, cnrai, caintc, cpeva, cptra, cevap, cgsnow, csubl, csnrai, cmelt |
| Per-event (tillage-reset) | 3 | ldwet, spev, saev |

### Cohort vs flat decision

Boundary + crop-uptake both used **flat** layout — all fields instantaneous, no cohort needed.

Atmosphere is **structurally different**: it owns 10 cumulatives + 8 intermediates + 3 per-event state fields. **Three reset cadences** vs boundary/crop-uptake's one. This is exactly the ADR 0033 cohort-partitioning use case.

**Three options:**

**Option A — Pure flat (33 scalars, no cohorts).** Match boundary/crop-uptake precedent. Reset cadences are encoded via `if (flzerocumu) state%atmosphere%cgrai = 0.0` style scattered through the code as today. Simple to implement, follows precedent. **Loses** the ADR 0033 cohort-reset clarity.

**Option B — Two cohorts (cumulative + intermediate) under flat top-level scalars.**
```fortran
type :: atmosphere_cumulative_t
   real(real64) :: cgrai, cnrai, caintc, cpeva, cptra, cevap, cgsnow, csubl, csnrai, cmelt
contains
   procedure :: reset => atmosphere_cumulative_reset
end type

type :: atmosphere_intermediate_t
   real(real64) :: igrai, inrai, ipeva, iptra, ievap, igsnow, isubl, isnrai
contains
   procedure :: reset => atmosphere_intermediate_reset
end type

type :: atmosphere_state_t
   ! Instantaneous + per-day scalars (20 fields flat)
   real(real64) :: peva, ptra, empreva, atmdem, ...
   ! ...
   type(atmosphere_cumulative_t)   :: cumu
   type(atmosphere_intermediate_t) :: intr
end type
```
Reset becomes `call state%atmosphere%cumu%reset()` invoked once per cumu cycle. Eliminates the 3-file reset scatter. **RECOMMENDED.** Aligns with ADR 0033, sets the pattern for the soil-water-core arc's cumulative cohort.

**Option C — Defer cumulatives to soil-water-core arc.** Carve only instantaneous + per-day fields in this arc (~20 fields), leave the 18 intermediate + cumulative scalars for soil-water-core to package with `iqrot/cqrot/inqrot` cumulatives. Smaller arc, but the meteoday-reset ownership question would persist for another arc.

**Recommendation: Option B (two cohorts).** Justification:
1. Atmosphere is the natural place to introduce the cumulative-cohort pattern (largest cumulative count of any subsystem encountered so far).
2. The meteoday refactor will benefit: the cumulative-reset block at meteoday:412–423 can be deleted in favor of `state%atmosphere%cumu%reset()` invoked from a unified driver.
3. ADR 0033 cohort-reset is already the documented convention; not using it here without a strong reason is debt.
4. Soil-water-core arc can re-use the cohort pattern (cqrot/iqrot/inqrot/iqredXXX_day, etc.).

If we go Option B, soilwater_init's role does NOT change (allocations are crop-uptake's). Atmosphere needs its own `atmosphere_init(state%atmosphere)` — pure scalars, no `numnod`/`nlay`.

### `soilwater_init` signature change

**No change.** All atmosphere fields are scalars. No new per-node arrays. The new `atmosphere_init` lives in a separate module (`atmosphere_state_mod`).

### Recommended `atmosphere_state_t` shape

```fortran
type :: atmosphere_state_t

   ! === Instantaneous (per-call, no reset gate) ===
   real(real64) :: peva     = 0.0_real64
   real(real64) :: ptra     = 0.0_real64
   real(real64) :: empreva  = 0.0_real64
   real(real64) :: melt     = 0.0_real64
   real(real64) :: subl     = 0.0_real64
   real(real64) :: slw      = 0.0_real64
   real(real64) :: ssnow    = 0.0_real64
   real(real64) :: snowinco = 0.0_real64
   real(real64) :: graidt   = 0.0_real64
   real(real64) :: nraidt   = 0.0_real64
   real(real64) :: aintcdt  = 0.0_real64

   ! === Per-day (overwritten each meteo day) ===
   real(real64) :: grai     = 0.0_real64
   real(real64) :: nraida   = 0.0_real64
   real(real64) :: atmdem   = 0.0_real64
   real(real64) :: pevaday  = 0.0_real64
   real(real64) :: ptraday  = 0.0_real64
   real(real64) :: gsnow    = 0.0_real64
   real(real64) :: snrai    = 0.0_real64
   real(real64) :: fprecnosnow = 0.0_real64
   real(real64) :: sicact   = 0.0_real64

   ! === Per-event (reset on tillage) ===
   real(real64) :: ldwet = 0.0_real64
   real(real64) :: spev  = 0.0_real64
   real(real64) :: saev  = 0.0_real64

   ! === Cohorts ===
   type(atmosphere_intermediate_t) :: intr
   type(atmosphere_cumulative_t)   :: cumu

end type atmosphere_state_t
```

Total: 22 top-level scalars + 18 fields under two cohorts = **40 fields organised under one record.**

---

## 9. Scope estimate

| Metric | Atmosphere (THIS) | Boundary (ADR 0035) | Crop-uptake (ADR 0036) |
|---|---|---|---|
| Owned globals | ~33 (11 inst + 9 per-day + 8 intr + 10 cumu + 3 per-event — fits 40 fields in record with 2 cohorts) | 12 | 22 |
| Owned-by cohort | 2 cohorts (intermediate + cumulative) | 0 cohorts (flat) | 0 cohorts (flat) |
| External reader files | ~9 (rootextraction, cropgrowth, tillage, boundtop, soilhydraulics, waterbalance, temperature, swapoutput, swap_csv_output) | 14 | 7 |
| External read sites | ~105 (heavy ptra in rootextraction; heavy cumulative reads in swapoutput) | 80 | 47 |
| Co-writers | 5 (meteoday, meteodt, soilhydraulics, waterbalance, irrigation) — 2 are EXCLUDED but must be touched | 5 | 1 |
| Phase 0 candidates | 1–2 (swsublim, spev/saev initial) | 7–8 | 0 |
| New init routine | YES — `atmosphere_init(state%atmosphere)` (scalars only, no dims) | n/a (in soilwater) | n/a (in soilwater, ext signature) |
| LoC home tree | 1 398 (et+interception+precip+snow) | 494 | 875 |

### Suggested task decomposition

**Phase 0 — config gaps**
- **A-0.1** (optional): Audit `swsublim` config coverage; add missing field to `meteorology_config.snow` if needed.
- **A-0.2** (optional): Audit `spev`/`saev` initial seeding; add to `soil_config.initial` if uncovered.

**Phase 1 — state-type extension + dual-write**
- **A-1.1**: Create `src/state/atmosphere_state.f90` with `atmosphere_state_t`, `atmosphere_intermediate_t`, `atmosphere_cumulative_t`, `atmosphere_init`, type-bound `reset()` methods. pFUnit alloc/zero/reset tests.
- **A-1.2**: Add `state%atmosphere` member to `swap_state_t` in `swap_state.f90`. Add `call atmosphere_init(state%atmosphere)` at `swap.f90:188` (after soilwater_init).
- **A-1.3**: Dual-write inside `snow.f90` — every legacy write of ssnow/melt/subl/slw/snowinco/igsnow/isubl/isnrai/cgsnow/csubl/csnrai/cmelt also writes `state%atmosphere%…`. Promote `state` arg from `optional intent(in)` → `intent(inout)` to support writes. Update swap.f90:216, 298 callers.
- **A-1.4**: Dual-write inside `precipitation.f90` — grai/gsnow/snrai/fprecnosnow + the ssnow=0 mutation. Add `state` arg to `PartitionPrecipitation`. Update meteoday:369 caller (one-line patch in excluded file — additive only).
- **A-1.5**: Dual-write inside `interception.f90` — nraida (from DivIntercep), sicact (from ruttervw). Add `state` arg to DivIntercep (and msw1eic if needed). Update meteoday:530, 653 callers.
- **A-1.6**: Dual-write inside `et.f90` — empreva/ldwet/spev/saev. Add `state` (optional intent(inout)) to `reduceva`. Update meteoday:808, meteodt:358, 448 callers.
- **A-1.7**: Dual-write inside `meteoday.f90` (EXCLUDED file — additive only) — peva/ptra/atmdem/pevaday/ptraday/grai (after mm→cm) + the cumu/intr resets (`call state%atmosphere%cumu%reset()` alongside `if (flzerocumu) cgrai=0`). ~25 lines additive.
- **A-1.8**: Dual-write inside `meteodt.f90` (EXCLUDED — additive only) — peva/ptra. ~6 lines additive.
- **A-1.9**: Dual-accumulate inside `waterbalance.f90` — igrai/inrai/ipeva/iptra/cpeva/cptra/cevap/cgrai/cnrai/caintc accumulate paths now also write the state cumu/intr. Reset of igrai/inrai (lines 388–389) also moves to state (single-owner cleanup — Hazard #3).
- **A-1.10**: Dual-reset inside `soilhydraulics.f90` — cpeva/cptra/cevap/ipeva/iptra/ievap/ldwet/spev/saev resets also clear state%atmosphere%cumu/intr (or migrate the resets out of soilhydraulics into the atmosphere cohort reset). Decision in design.

**Phase 2 — external reader migration**
- **A-2.1**: `rootextraction.f90` cut over `ptra` (25 sites) + `atmdem` (3 sites) reads to `state%atmosphere%ptra/atmdem`. Already takes `state` per ADR 0036.
- **A-2.2**: `cropgrowth.f90` cut over `ptra` reads (9 sites). Already takes `state` per ADR 0036 C-1.3.
- **A-2.3**: `tillage.f90` cut over `nraida` reads (5 sites). State plumbing: tillage takes `state` per ADR 0020.
- **A-2.4**: `boundtop.f90` cut over `peva`/`empreva`/`melt` reads (3 sites). Already state-aware (ADR 0035).
- **A-2.5**: `soilhydraulics.f90` cut over `melt` reads (4 sites in F-vector) + retire the duplicate reset block (A-1.10 moves them).
- **A-2.6**: `waterbalance.f90` cut over `ipeva`/`iptra`/`igrai`/`inrai`/cumulative reads (~10 sites) + retire dup reset (A-1.9).
- **A-2.7**: `temperature.f90` cut over `ssnow` reads (2 sites).
- **A-2.8**: `swapoutput.f90` cut over output blocks — `.bal` cumulative cluster (lines 212–344), `.inc` intermediate cluster (lines 381–558), `.csv` cluster, ssnow/snowinco/sicact direct reads. ~30 sites.
- **A-2.9**: `swap_csv_output.f90` cut over EPOT/TPOT/SSNOW + dstor (~6 sites).
- **A-2.10**: Retire legacy globals from `variables.f90`: peva, ptra, empreva, atmdem, pevaday, ptraday, grai, nraida, gsnow, snrai, fprecnosnow, sicact, ssnow, snowinco, slw, melt, subl, ldwet, spev, saev, aintcdt, igrai, inrai, ipeva, iptra, ievap, igsnow, isubl, isnrai, cgrai, cnrai, caintc, cpeva, cptra, cevap, cgsnow, csubl, csnrai, cmelt. Drop zero-inits from initialize.f90. Drop reset blocks now redundant. Expect 3–6 compile-driven fixup commits (playbook lesson #5).

**Total suggested tasks: ~17–19** (2 Phase 0 + 10 Phase 1 + 10 Phase 2). Larger than crop-uptake (11) and boundary (16) due to:
- Cohort introduction (cumu + intr) — first arc using cohort pattern.
- Two excluded files that must still be touched additively (meteoday, meteodt).
- Three reset-site mergers (Hazard #2, #3).
- Snow signature promotion (Hazard #7).

---

## 10. Open questions for design phase

1. **Cohort vs flat — Recommendation: TWO COHORTS (intermediate + cumulative) under flat top-level scalars (Option B).** Atmosphere is the natural debut for the cohort pattern: it has 10 cumulatives + 8 intermediates spread across 4 reset sites. Soil-water-core will reuse the pattern. The alternative (Option A pure flat) does not lose data semantics but does not solve the multi-owner reset hazard. Option C (defer cumulatives to soil-water-core) postpones the meteoday reset cleanup.

2. **meteoday/meteodt dual-writes — Touch them additively in this arc or defer to the meteo refactor?** **Recommendation: ADDITIVE dual-writes in this arc.** ~30 total lines across two excluded files. Adds state writes alongside legacy writes; the meteo refactor remains free. Defer is technically possible but leaves peva/ptra/atmdem untouched and forces every reader (rootextraction × 25 sites, cropgrowth × 9 sites) to wait. Bad cost/benefit.

3. **pond read in et.f90:reduceva — Defer or migrate?** **Recommendation: defer.** Add `state` arg to reduceva (optional intent(inout) for empreva/ldwet/spev/saev writes); read `pond` via legacy global until soil-water-core arc migrates `pond` (boundary D5 deferred). Equal to the boundary precedent for kmean.

4. **`reduceva` globals soup — Address in this arc?** **Recommendation: address only the 5 atmosphere-owned reads/writes (empreva, ldwet, spev, saev, peva).** Other reads (swredu, fldaystart, cofred, dt, rsigni, nird, pond) stay legacy. The signature becomes `subroutine reduceva(task, nrai, state)` with state optional (for backward compat with meteodt callsites that haven't been migrated yet).

5. **soilhydraulics-owned resets of cpeva/cptra/cevap/ipeva/iptra/ievap — move out or keep?** **Recommendation: MOVE.** Atmosphere cumulatives reset belongs in `state%atmosphere%cumu%reset()` invoked from a single owner (probably from a new `atmosphere_reset` block, OR keep the call from inside SoilWater(1) but invoke the type-bound reset). Removes ownership scatter, supports the cohort pattern. Note: cpeva/cptra/cevap are tied to mass-balance state in soilhydraulics — verify the reset gate timing aligns with `flzerocumu` and not a different gate. Quick grep showed both are `flzerocumu`-gated.

6. **`igrai`/`inrai` double-reset — Which file wins?** **Recommendation: atmosphere wins.** Atmosphere is the conceptual owner (grai is the upstream forcing). Delete the resets from waterbalance.f90:388–389 and the meteoday:414–415 block in favor of `state%atmosphere%intr%reset()`. This is Hazard #3 resolution.

7. **`nird` / `gird` — irrigation or atmosphere?** **Recommendation: irrigation-owned, defer in this arc.** DivIntercep co-writes `nird` but the canonical owner is irrigation.f90. Mark as deferred; will be migrated when irrigation arc happens (or fold into atmosphere if a future arc considers them part of the meteo/irrigation drive).

8. **Phase 0 candidates — `swsublim`, `spev`/`saev` initial?** **Recommendation: verify in design phase.** Both are minor (1–2 commits at most). Run `grep -n "swsublim" src/config/` and `grep -nE "spev|saev" src/io/toml/`. If absent, add defaults to `meteorology_config.snow` / `soil_config.initial`.

9. **`Snow` signature promotion (optional intent(in) → intent(inout)) — Ripple effects?** **Recommendation: do the promotion in A-1.3.** Caller sites `swap.f90:216, 298` pass full state, so the change is mechanical. The 3 cropgrowth heat-arc rename-tricks do NOT involve snow — no rename trick needed.

10. **`peva = 0` and `empreva = 0` in snow.f90:108–109 — Is this a hidden state mutation cliff?** **Recommendation: document but keep behavior.** snow.f90 zeroes peva to prevent further soil-evaporation calculation when snow exists. This is a deliberate cascade. After migration: `state%atmosphere%peva = 0.0` and `state%atmosphere%empreva = 0.0` — same behavior, same gate. No semantic change. Worth a one-line ADR note about state-mutating side effect from snow.

11. **`precipitation.f90:127` ssnow mutation under swmetdetail=1 — Is this real or vestigial?** **Recommendation: keep verbatim.** Comment in code says "consider refactoring." Out of scope for this arc; just retarget the write to `state%atmosphere%ssnow = 0.0d0`.

12. **Compile-driven Phase 2.7 expected surprise count?** **Estimate: 3–6 fixup commits.** Heaviest surface is swapoutput.f90's use clauses (large composite cumulative block at lines 212–344 + line 381–558). Expect 3–4 hidden import surprises there. swap_csv_output.f90 ditto.

13. **Should this arc consume mfluxtable allocation slot or leave it alone?** **Out of scope.** mfluxtable is crop-uptake state (ADR 0036). No interaction.

14. **Could the meteo refactor INSTEAD own the atmosphere state-carve?** Worth raising. If the user's meteo refactor is restructuring meteoday/meteodt, it could absorb the dual-writes within its own scope, eliminating the "touch excluded files" pattern of this arc. **Open question for the user.** Default plan assumes parallel-arc additive dual-writes; an alternative is to **defer atmosphere arc until after meteo refactor merges**, then carve from cleaner files.

---

## Discovery summary

- **Owned set: ~33 atmosphere fields** across 5 cadence classes (11 instantaneous + 9 per-day + 8 intermediate + 10 cumulative + 3 per-event). New record: `atmosphere_state_t` with two cohort sub-records (`atmosphere_cumulative_t`, `atmosphere_intermediate_t`).
- **Layout: TWO COHORTS** (cumulative + intermediate) under flat top-level scalars. **First arc using ADR 0033 cohort pattern.** Sets template for soil-water-core.
- **External reader files: ~9 in-scope** (rootextraction, cropgrowth, tillage, boundtop, soilhydraulics, waterbalance, temperature, swapoutput, swap_csv_output). ~105 read sites — heaviest is the `swapoutput.f90` cumulative block (~30 sites).
- **Co-writers: 5** — including 2 EXCLUDED files (`meteoday.f90`, `meteodt.f90`) that must be touched additively for dual-write. ~30 additive lines across both.
- **Signature status: only `snow` takes `state` (optional intent(in))**. Promotion to `intent(inout)` needed. `reduceva` needs a new optional `state` arg. `PartitionPrecipitation`, `DivIntercep`, `VonHHBraden`, `Gash`, `ruttervw` need state args.
- **`atmosphere_init(state%atmosphere)` signature: scalars only**, no `numnod`/`nlay` parameter. Called from swap.f90:188 after soilwater_init.
- **Hazards: 11 items.** Critical: (#1) meteoday/meteodt are excluded but write peva/ptra/atmdem; (#2) multi-owner cumulative reset pattern needs cohort consolidation; (#3) igrai/inrai double-reset; (#7) snow signature promotion ripple. Manageable: pond read (defer), reduceva globals, temperature ssnow read, init-order, ssnow mutation in precipitation.
- **Phase 0 gaps: 1–2 minor fields** (swsublim, spev/saev initial seeds).
- **Suggested task count: ~17–19** across Phase 0/1/2. Larger than crop-uptake; smaller than the eventual soil-water-core.
- **Compile-driven Phase 2.7 expectation: 3–6 fixup commits.**
- **Decomposition: monolithic arc** with optional pre-arc question (Q14) about whether to wait for meteo refactor merge.

End of discovery.
