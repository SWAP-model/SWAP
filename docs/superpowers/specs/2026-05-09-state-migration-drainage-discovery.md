# Subsystem Migration Discovery: Drainage

**Date:** 2026-05-09
**Status:** discovery (read-only inventory)
**Migration #:** 2 of N (surface-water was #1)
**Branch:** `refactor/surfacewater-state` (continuing the umbrella migration branch)
**Predecessor playbook:** `docs/superpowers/specs/2026-05-08-state-migration-surfacewater-discovery.md`

---

## 1. Source files

| File | LoC | Role |
|------|-----|------|
| `src/drainage/drainage.f90` | 822 | Home: three public subroutines `drainage()`, `bocodrb()`, `bocodre()`. Main compute orchestration, calls `bocodrb` or is itself called for the swdra=1 (no surface water management) path. `bocodre` handles the swsrf>=2 (full surface water management) path. |
| `src/drainage/divdra.f90` | 599 | Home: `DIVDRA()` — transmissivity-weighted flux distribution over soil compartments; `Lev2Comp()` — helper to locate soil compartment containing a given depth level. |
| `src/config/drainage_config.f90` | 231 | Typed config: `drainage_config_t` struct (TOML-sourced, validated). Also `drainage_surface_runoff_t` sub-section. |
| `src/io/toml/read_drainage_toml.f90` | 227 | TOML reader: `read_drainage_toml()` and `read_drainage_inner()`. Decodes `[drainage]` (with optional external file via `drainage.file`), per-level `[[drainage.levels]]`, and `[drainage.surface_runoff]` sub-section. Converts `L` from metres to centimetres at parse time (spec D6). |
| `tests/unit/config/test_drainage_config.pf` | 349 | pFUnit: validation rules for `drainage_config_t` (enum checks, cross-field constraints dramet/swdivd, swdra=2+dramet, altcu stub-error, swliminf cross-field). |
| `tests/unit/io/toml/test_read_drainage_toml.pf` | 310 | pFUnit: TOML reader round-trip tests (minimal case, levels, surface_runoff, external file, swliminf). |
| `tests/unit/drainage/test_surfacewater_init.pf` | (in drainage/ dir) | pFUnit: `surfacewater_init` zeroing test — NOT drainage-subsystem-owned; belongs to surface-water. Included here because it physically lives under `tests/unit/drainage/`. |

---

## 2. Owned globals (this subsystem writes)

A "global" is any variable declared in `src/core/variables.f90` and accessed via `use variables, only: …`. Owned = subsystem mutates it.

| Variable | Type | `variables.f90` line | Description (from inline comment) | Write sites (file:line) |
|----------|------|----------------------|-----------------------------------|--------------------------|
| `qdrain(Madr)` | `real(8)` | 881 | Total lateral drainage flux (L/T) for each drainage level | `drainage.f90:182,191,220,235,251,284,286,287,300,301,307,480,659,718,733,735,743,751,779,813` |
| `qdra(Madr,macp)` | `real(8)` | 880 | Lateral drainage flux per drainage level and compartment | `drainage.f90:528,530,533,542,544` (see also Note A below) |
| `qdrd` | `real(8)` | 1294 | Total flux to/from secondary surface water system (cm/d) | `drainage.f90:665,761,815` (only in `bocodre`, i.e. `swsrf>=2` path) |
| `drainl(Madr)` | `real(8)` | 784 | Drainage level (max of drain bottom and surface water level) | `drainage.f90:690,700,746,754` (only in `bocodre`) |
| `wetper(Madr)` | `real(8)` | 967 | Wet perimeter of drain for each drainage level | `drainage.f90:695,704,740` (only in `bocodre`) |
| `ztopdislay(Madr)` | `real(8)` | 976 | Depth of top of model discharge layer per drain level | `drainage.f90:497` (only in `drainage()` subroutine, `swdislay=2` path) |

**Note A — `qdra` dual-ownership situation:**
`qdra` is written by the drainage subsystem in `drainage()` (lines 528–544: the `swdivd=0` and `swdivd=1` redistribution branches). It is ALSO written by `surfacewater.f90` (lines 133–134, 181–186, 219–220) which is the surface-water subsystem home file. After surface-water Phase 2, `qdra` still exists as a global because `frozencond.f90` calls `DIVDRA` directly using the global (see Section 8, Hazard #1). The canonical value is dual-tracked: both the global `qdra` and `state%surfacewater%qdra` are kept current by drainage (line 556–560) and by frozencond (line 282).

**Note B — `qdrain` and `qdra` relationship to `state%surfacewater`:**
After surface-water Phase 2 Task 11, `qdrain` is still a global (needed by `divdra` callers: `surfacewater.f90:143`, `frozencond.f90:276`). `qdrtot` was migrated to `state%surfacewater%qdrtot` (variables.f90:884 — commented out). The cumulative arrays `cqdrain`, `cqdrainin`, `cqdrainout`, `cqdra`, `inqdra`, `inqdra_in`, `inqdra_out`, `iqdra` were all migrated to `state%surfacewater` and removed from globals (variables.f90:767–840 — all commented out). Drainage still writes these via `state%surfacewater` (not globals).

---

## 3. Borrowed globals (this subsystem reads, owned elsewhere)

Variables read by the drainage home files but not written by them. Grouped by borrowing subroutine.

| Variable | Best-guess owner subsystem | Read sites (file:line) |
|----------|----------------------------|--------------------------|
| `gwl` | soil/soilhydraulics | `drainage.f90:164,168,171,478,497,498,681,686,711,712,713` |
| `dramet` | config (set once at init) | `drainage.f90:167,255,306,432,435` |
| `zbotdr(Madr)` | config (set once at init) | `drainage.f90:171,182,224,225,239,293,307` |
| `basegw` | config (set once at init) | `drainage.f90:171` |
| `l(Madr)` | config (set once at init) | `drainage.f90:171,190,197,198,203,208,216,218,219,232,233,246,248,284,286,300,307,733,734` |
| `ipos` | config (set once at init) | `drainage.f90:187,194,223,238` |
| `khtop` | config (set once at init) | `drainage.f90:190,216` |
| `khbot` | config (set once at init) | `drainage.f90:218,232` |
| `kvtop` | config (set once at init) | `drainage.f90:230` |
| `kvbot` | config (set once at init) | `drainage.f90:231,233` |
| `entres` | config (set once at init) | `drainage.f90:190,216,219,235,251` |
| `zintf` | config (set once at init) | `drainage.f90:224,239` |
| `geofac` | config (set once at init) | `drainage.f90:248` |
| `swdtyp(Madr)` | config (set once at init) | `drainage.f90:275,294,435,688,694,700,732,739,768` |
| `owltab(Madr,2*maowl)` | config/timeseries | `drainage.f90:267,307,442,492` |
| `nowltab(Madr)` | config/timeseries | `drainage.f90:267,442,492`; `divdra.f90:444` |
| `t1900` | time / timecontrol | `drainage.f90:267,442,492,793` |
| `dt` | time / timecontrol | `drainage.f90:267,284,490,796,803` |
| `swallo(Madr)` | config (set once at init) | `drainage.f90:287,301` |
| `drares(Madr)` | config (set once at init) | `drainage.f90:286` |
| `infres(Madr)` | config (set once at init) | `drainage.f90:300` |
| `qdrtab(50)` | config (set once at init) | `drainage.f90:307` |
| `nrlevs` | config (set once at init) | `drainage.f90:257,426,457,478,480,489,550` |
| `numnod` | config/soil geometry | `drainage.f90:456,489,527,540` |
| `swnrsrf` | config (set once at init) | `drainage.f90:283,717` |
| `cofintfl` | config (set once at init) | `drainage.f90:284,718` |
| `expintfl` | config (set once at init) | `drainage.f90:284,718` |
| `shape` | config (set once at init) | `drainage.f90:168` |
| `FlMacropore` | macropore subsystem / config | `drainage.f90:275,767` |
| `NumLevRapDra` | config (set once at init) | `drainage.f90:275,426,435,442,767` |
| `swliminf` | config (set once at init) | `drainage.f90:295` |
| `swdivd` | config (set once at init) | `drainage.f90:489` |
| `swdislay` | config (set once at init) | `drainage.f90:494,502` |
| `swtopdislay(Madr)` | config (set once at init) | `drainage.f90:496,504` |
| `fTopDisLay(Madr)` | config (set once at init) | `drainage.f90:497,498` |
| `dz(macp)` | soil geometry | `drainage.f90:490,507,514` |
| `ksatfit(maho)` | soil hydraulics | `drainage.f90:490` |
| `ksatexm(maho)` | soil hydraulics | `drainage.f90:490` |
| `fluseksatexm(macp)` | soil hydraulics | `drainage.f90:490` |
| `layer(macp)` | soil geometry | `drainage.f90:490` |
| `cofani(maho)` | config (drainage anisotropy) | `drainage.f90:490` |
| `Swdivdinf` | config (set once at init) | `drainage.f90:490` |
| `SwTopnrsrf` | config (set once at init) | `drainage.f90:490` |
| `FacDpthInf` | config (set once at init) | `drainage.f90:490` |
| `madr` | array dimension constant | `drainage.f90:388` |
| `flzerointr` | time / timecontrol | `drainage.f90:455` |
| `flzerocumu` | time / timecontrol | `drainage.f90:468` |
| `swsec` | config / surfacewater | `drainage.f90:775` |
| `swsrf` | config / surfacewater | `drainage.f90:670,671,687,758,782` |
| `nrpri` | config / surfacewater | `drainage.f90:671,812` |
| `taludr(Madr)` | config (set once at init) | `drainage.f90:703` |
| `widthr(Madr)` | config (set once at init) | `drainage.f90:695,701` |
| `pond` | soil boundary / surfacewater | `drainage.f90:712` |
| `pondmx` | config (surface water) | `drainage.f90:681` |
| `wlp` | surfacewater (set by `surfacewater.f90:237`) | `drainage.f90:674` |
| `gwlinf(Madr)` | config (set once at init) | `drainage.f90:713,714` |
| `rdrain(Madr)` | config (set once at init) | `drainage.f90:721,722` |
| `rinfi(Madr)` | config (set once at init) | `drainage.f90:729` |
| `rentry(Madr)` | config (set once at init) | `drainage.f90:722` |
| `rexit(Madr)` | config (set once at init) | `drainage.f90:730` |
| `impend(mamp)` | config / surfacewater management | `drainage.f90:793` |
| `nmper` | config / surfacewater management | `drainage.f90:789` |
| `wscap(mamp)` | config / surfacewater management | `drainage.f90:796,803` |
| `rsurfdeep` | config (set once at init) | `drainage.f90:725` |
| `rsurfshallow` | config (set once at init) | `drainage.f90:727` |

---

## 3.5 External readers of owned globals (cross-subsystem reader inventory)

**Populated UPFRONT per playbook update from surface-water Phase 2 lessons-learned.**

For each drainage-owned global from Section 2, the files OUTSIDE the drainage home tree that read it.

| Owned global | External reader files |
|---|---|
| `qdrain(Madr)` | `src/drainage/surfacewater.f90:143,177,219,220,227` (compute — `divdra` call arg, redistribution, qdrtot summation); `src/soil/waterbalance.f90:513,514,516,517,519` (compute — cumulative drainage accounting); `src/heat/frozencond.f90:251,258,265,269,276,289,294,303` (compute — frozen boundary redistribution, also WRITES qdrain: see Hazard #2) |
| `qdra(Madr,macp)` | `src/drainage/surfacewater.f90:133,134,143,170,172,181,183,184,186,209,212,219,220` (compute — also WRITES qdra: see Hazard #3); `src/heat/frozencond.f90:249,250,276,282,291,292,293` (compute — also WRITES qdra: see Hazard #2); `src/soil/soilhydraulics.f90:89,92,94` (compute — reads from `state%surfacewater%qdra`); `src/soil/waterbalance.f90:347,349,443,444,446,449` (compute — reads both global `qdra` and `state%surfacewater%qdra`); `src/solute/solute.f90:131,197,199,200,202` (compute — reads from `state%surfacewater%qdra`) |
| `qdrd` | `src/drainage/surfacewater.f90:332,474,503,509,534,569,589,648,686,722,736` (compute — surface water balance, WBALLEV) |
| `drainl(Madr)` | (drainage-only — no external readers found) |
| `wetper(Madr)` | `src/io/toml/config_to_variables.f90:292` (config — sets `wetper(1)` for `dramet=2` at init time, not a runtime reader) |
| `ztopdislay(Madr)` | (drainage-only — used only in `drainage.f90` and `surfacewater.f90`; `surfacewater.f90` reads and writes the same global: see Hazard #4) |

---

## 4. Entry points (subroutines called from outside this subsystem)

| Subroutine (this subsystem) | Called from (file:line) | Purpose / lifecycle stage |
|------------------------------|--------------------------|--------------------------------------------------------------------|
| `drainage()` | `src/core/swap.f90:290` | per-step — main drainage flux computation for `swdra=1` (no surface water management) path; called inside the `fldtreduce` do-while loop |
| `drainage()` | `src/io/swapoutput.f90:3729` | per-step (inside perturbation loop) — called from `stocot1` finite-difference sensitivity routine; uses private `state_om` clone (see ADR 0030 / Hazard #5) |
| `drainage()` | `src/macropore/macropore.f90:121` | init (once) — called to initialize `state%surfacewater%ZDraBas` when `swdra=1` and macropore drainage basis not yet set |
| `bocodre()` | `src/drainage/surfacewater.f90:127` | per-step — called from `SurfaceWater(task=2)` for `swsrf>=2` (full surface water management) path |
| `DIVDRA()` | `src/drainage/surfacewater.f90:142` | per-step — called from `SurfaceWater(task=2)` to distribute fluxes over compartments (surface-water path) |
| `DIVDRA()` | `src/heat/frozencond.f90:275` | per-step — called from `FrozenBounds()` after frozen-zone zeroing to redistribute thaw-adjusted drainage (see Hazard #2) |

---

## 5. Internal call graph

```
drainage_mod (drainage.f90)
├── drainage(state)                       [public entry point]
│   ├── bocodrb(dh, state)                [internal: computes qdrain(:)]
│   │   └── afgen(...)                    [from array_utils]
│   └── DIVDRA(...)                       [from distribute_drainage]
│       └── Lev2Comp(...)                 [internal to divdra.f90]
├── bocodrb(dh, state)                    [public entry point]
│   └── afgen(...)                        [from array_utils]
└── bocodre(dh, state)                    [public entry point]
    └── fatalerr_collected(...)           [from error_mod]

distribute_drainage (divdra.f90)
└── DIVDRA(...)                           [public entry point]
    ├── Lev2Comp(...)                     [internal helper]
    └── afgen(...)                        [from array_utils]
```

---

## 6. Output coupling

What does the output layer read from drainage's state?

| Output writer (file:subroutine) | Reads which fields | Destination file (.csv, .drf, …) |
|-----------------------------------|----------------------|--------------------------------------|
| `src/io/swapoutput.f90:outdrf` (line 3037) | `state%surfacewater%cqdrain(1..5)`, `state%surfacewater%cqdrd` | `.drf` file (drainage flux per level + secondary system) |
| `src/io/swapoutput.f90` (waterbal section, ~line 974) | `state%surfacewater%cqdrain(:)`, `state%surfacewater%cqdra` | `.a` waterbalance file |
| `src/io/swapoutput.f90` (soil solute waterbal, ~line 2350) | `state%surfacewater%cqdrainin(:)`, `state%surfacewater%cqdrainout(:)` | `.blc` solute balance file |
| `src/io/swapoutput.f90` (ConvertDiscrVert / `.afr` file, ~line 1782) | `state%surfacewater%inqdra(:,:)` | `.afr` detailed flux file |
| `src/io/swap_csv_output.f90:set_values` (line 237) | `state%surfacewater%iqdra` | `swap.csv` DRAINAGE column |
| `src/io/swap_csv_output.f90:fill_values` (line 494,502,503,504) | `state%surfacewater%inqdra(:,:)`, `state%surfacewater%inqdra_in(:,:)`, `state%surfacewater%inqdra_out(:,:)` | `swap.csv` DRA / QDRA / QDIN / QDOU node-level columns |
| `src/soil/waterbalance.f90` (line 336,406) | `state%surfacewater%qdrtot` | in-memory flux calculation (not file, but output-visible via waterbalance) |
| `src/io/swapoutput.f90` (waterbal, ~line 307,362,457) | `state%surfacewater%cqdra`, `state%surfacewater%iqdra` | `.a` and summary files |
| `src/io/swapoutput.f90` (waterbal, line 578,624,645) | `qdraincomp(macp)` (global, aggregated from `qdra` in waterbalance.f90) | `.afr` node-level flux file |

---

## 7. Config inputs

| TOML key | Typed-config field path | Used at (file:line) |
|-------------------------------------|---------------------------|----------------------|
| `drainage.swdra` | `config%drain%swdra` → `swdra` | `config_to_variables.f90:269`; runtime gate `fldrain` in `swap.f90:68` |
| `drainage.dramet` | `config%drain%dramet` → `dramet` | `config_to_variables.f90:270`; `drainage.f90:167,255,306,432` |
| `drainage.drfil` | `config%drain%drfil` → `drfil` | `config_to_variables.f90:281` (legacy .dra file stem) |
| `drainage.swdivd` | `config%drain%swdivd` → `swdivd` | `config_to_variables.f90:271`; `drainage.f90:489` |
| `drainage.swdislay` | `config%drain%swdislay` → `swdislay` | `config_to_variables.f90:272`; `drainage.f90:494,502` |
| `drainage.nrlevs` | `config%drain%nrlevs` → `nrlevs` | `config_to_variables.f90:273`; finalized in `drainage_config_finalize` (forces 1 for dramet=1/2) |
| `drainage.swliminf` | `config%drain%swliminf` → `swliminf` | `config_to_variables.f90:418`; `drainage.f90:295` |
| `drainage.altcu` | `config%drain%altcu` | stub-error if nonzero (drainage_config.f90:145) |
| `drainage.basic.basegw` | `config%drain%basegw` → `basegw` | `config_to_variables.f90:274`; `drainage.f90:171` |
| `drainage.basic.entres` | `config%drain%entres` → `entres` | `config_to_variables.f90:275`; `drainage.f90:190,216,219,235,251` |
| `drainage.basic.shape` | `config%drain%shape` → `shape` | `config_to_variables.f90:294` (dramet=2 only); `drainage.f90:168` |
| `drainage.basic.lm` | `config%drain%lm` → `L(1)` (×100 cm) | `config_to_variables.f90:291` (dramet=2 only) |
| `drainage.basic.wetper` | `config%drain%wetper` → `wetper(1)` | `config_to_variables.f90:292` (dramet=2 only) |
| `drainage.basic.zbotdr` | `config%drain%zbotdr_basic` → `zbotdr(1)` | `config_to_variables.f90:293` (dramet=2 only) |
| `drainage.basic.ipos` | `config%drain%ipos` → `ipos` | `config_to_variables.f90:295` (dramet=2 only); `drainage.f90:187,194,223,238` |
| `drainage.basic.khtop` | `config%drain%khtop` → `khtop` | `config_to_variables.f90:296`; `drainage.f90:190,216` |
| `drainage.basic.khbot` | `config%drain%khbot` → `khbot` | `config_to_variables.f90:298` (ipos>=3); `drainage.f90:218,232` |
| `drainage.basic.zintf` | `config%drain%zintf` → `zintf` | `config_to_variables.f90:299` (ipos>=3); `drainage.f90:224,239` |
| `drainage.basic.kvtop` | `config%drain%kvtop` → `kvtop` | `config_to_variables.f90:302` (ipos>=4); `drainage.f90:230` |
| `drainage.basic.kvbot` | `config%drain%kvbot` → `kvbot` | `config_to_variables.f90:303` (ipos>=4); `drainage.f90:231` |
| `drainage.basic.geofac` | `config%drain%geofac` → `geofac` | `config_to_variables.f90:306` (ipos=5); `drainage.f90:248` |
| `drainage.cofani` | `config%drain%cofani(:)` → `cofani(:)` | `config_to_variables.f90:311-314`; `drainage.f90:490` |
| `drainage.levels[i].swdtyp` | `config%drain%swdtyp(i)` → `swdtyp(i)` | `config_to_variables.f90:317-320`; `drainage.f90:275,294,688,732,768` |
| `drainage.levels[i].zbotdr` | `config%drain%zbotdr(i)` → `zbotdr(i)` | `config_to_variables.f90:322-325`; `drainage.f90:171,269,688,690` |
| `drainage.levels[i].drares` | `config%drain%drares(i)` → `drares(i)` | `config_to_variables.f90:327-330`; `drainage.f90:286` |
| `drainage.levels[i].infres` | `config%drain%infres(i)` → `infres(i)` | `config_to_variables.f90:332-335`; `drainage.f90:300` |
| `drainage.levels[i].L` | `config%drain%L(i)` (cm, converted m→cm at parse) → `L(i)` | `config_to_variables.f90:337-342`; `drainage.f90:733` |
| `drainage.levels[i].gwlinf` | `config%drain%gwlinf(i)` → `gwlinf(i)` | `config_to_variables.f90:344-347`; `drainage.f90:713,714` |
| `drainage.levels[i].rdrain` | `config%drain%rdrain(i)` → `rdrain(i)` | `config_to_variables.f90:349-352`; `drainage.f90:721` |
| `drainage.levels[i].rinfi` | `config%drain%rinfi(i)` → `rinfi(i)` | `config_to_variables.f90:354-357`; `drainage.f90:729` |
| `drainage.levels[i].rentry` | `config%drain%rentry(i)` → `rentry(i)` | `config_to_variables.f90:359-362`; `drainage.f90:722` |
| `drainage.levels[i].rexit` | `config%drain%rexit(i)` → `rexit(i)` | `config_to_variables.f90:364-367`; `drainage.f90:730` |
| `drainage.levels[i].widthr` | `config%drain%widthr(i)` → `widthr(i)` | `config_to_variables.f90:369-372`; `drainage.f90:695,701` |
| `drainage.levels[i].taludr` | `config%drain%taludr(i)` → `taludr(i)` | `config_to_variables.f90:374-377`; `drainage.f90:703` |
| `drainage.levels[i].swallo` | `config%drain%swallo(i)` → `swallo(i)` | `config_to_variables.f90:379-382`; `drainage.f90:287,301` |
| `drainage.levels[i].owltab_file` | `config%drain%owltab_file(i)` → `owltab(i,:)` / `nowltab(i)` | `config_to_variables.f90:391-415`; `drainage.f90:267,307,442,492`; `divdra.f90:444` |
| `drainage.surface_runoff.swnrsrf` | `config%drain%surface_runoff%swnrsrf` → `swnrsrf` | `config_to_variables.f90:425`; `drainage.f90:283,717` |
| `drainage.surface_runoff.swtopnrsrf` | → `SwTopnrsrf` | `config_to_variables.f90:426`; `drainage.f90:490` |
| `drainage.surface_runoff.swdivdinf` | → `swdivdinf` | `config_to_variables.f90:427`; `drainage.f90:490` |
| `drainage.surface_runoff.facdpthinf` | → `FacDpthInf` | `config_to_variables.f90:428`; `drainage.f90:490` |
| `drainage.surface_runoff.cofintfl` | → `cofintfl` | `config_to_variables.f90:429`; `drainage.f90:284,718` |
| `drainage.surface_runoff.expintfl` | → `expintfl` | `config_to_variables.f90:430`; `drainage.f90:284,718` |
| `drainage.surface_runoff.geofac` | → `geofac` (surface_runoff sub-section value) | `config_to_variables.f90:431` |
| `drainage.surface_runoff.rsurfdeep` | → `rsurfdeep` | `config_to_variables.f90:442`; `drainage.f90:725` |
| `drainage.surface_runoff.rsurfshallow` | → `rsurfshallow` | `config_to_variables.f90:443`; `drainage.f90:727` |
| `drainage.surface_runoff.rapdrareaexp` | → `RapDraReaExp` | `config_to_variables.f90:444` |
| `drainage.surface_runoff.numlevrapdra` | → `NumLevRapDra` | `config_to_variables.f90:445`; `drainage.f90:275,426,435,442,767` |
| `drainage.surface_runoff.swtopdislay` | → `swtopdislay(:)` (broadcast) | `config_to_variables.f90:448-452`; `drainage.f90:496,504` |
| `drainage.surface_runoff.ftopdislay` | → `ftopdislay(:)` (broadcast) | `config_to_variables.f90:453-457`; `drainage.f90:497,498` |
| `drainage.surface_runoff.rapdraresref` | → `RapDraResRef(:)` (broadcast) | `config_to_variables.f90:458-462` |

---

## 8. Open questions / cross-subsystem hazards

1. **Dual-ownership of `qdra` global (drainage + surfacewater.f90 + frozencond.f90):**
   `qdra` is written by the drainage subsystem (`drainage.f90:528–544`), by the surface-water subsystem (`surfacewater.f90:133–134,181–186,209–212,219–220`), and by `frozencond.f90:249–250,279–282,291–293`. All three must be kept consistent. After surface-water Phase 2 the canonical value is `state%surfacewater%qdra` with the global as a shadow; but `frozencond` still reads and writes the global (confirmed by comment at `frozencond.f90:184–186`). When migrating drainage to typed state, the drainage-owned write at lines 528–544 must be changed, AND `frozencond`'s direct global write must be addressed in the same migration wave.

2. **`frozencond.f90` writes `qdrain` AND calls `DIVDRA` using the global:**
   `FrozenBounds()` (`frozencond.f90:178`) modifies `qdrain(level)` at lines 251,265,269,289,294 — after `Drainage()` has computed it. It also calls `DIVDRA(…, qdrain, qdra, …)` at line 276 passing the global arrays as explicit-shape arguments. When drainage migrates to typed state, `FrozenBounds` will need to be updated to read/write from `state%drainage%qdrain` (or equivalent), and the `DIVDRA` call signature will need adjustment. This is a hard dependency on the migration.

3. **`surfacewater.f90` writes to `qdrain` and `qdra` globals:**
   `surfacewater.f90:133–134` zeroes `qdra(level,node)` and `state%surfacewater%qdra`. Lines 181–186 redistribute `qdra`. Lines 219–220 set `qdra(level,numnod)`. The surface-water and drainage subsystems share write-ownership of both `qdrain` and `qdra`. A clean drainage state type must define which subsystem's state owns the final per-compartment values. The current design (`state%surfacewater`) is the path of least resistance but conflates the two subsystems.

4. **`ztopdislay` written by `drainage.f90` AND read/written by `surfacewater.f90`:**
   `surfacewater.f90:150–151` writes `zTopDisLay(level)` (same global `variables.f90:976`) when `swdislay=2`. This is the same global that `drainage()` writes at `drainage.f90:497`. Both subroutines contain near-identical redistribution logic (compare `drainage.f90:494–537` with `surfacewater.f90:146–204`). This is likely code duplication introduced during the earlier split between `swdra=1` and `swsrf>=2` paths. Migrating `ztopdislay` to typed state will require reconciling both write sites.

5. **Compute routine called from output (`swapoutput.f90:3729`):**
   `outdrf`'s parent routine `stocot1` calls `Drainage(state_om)` at `swapoutput.f90:3729` as part of a finite-difference perturbation experiment (ADR 0030). This is a legitimate compute call from the output file, not dead code. The `state_om` is a SAVE-local clone of the main `state`, so typed state migration must ensure `state%drainage` (when it exists) is also cloned correctly. The ADR 0030 note at `swapoutput.f90:3576–3583` documents the rationale.

6. **`wlp` (surfacewater-owned global) read by `bocodre`:**
   `drainage.f90:674` reads `wlp` which is set by `surfacewater.f90:237`. This is a cross-subsystem read of a surfacewater-owned runtime variable (not a config-time value). When the two subsystems eventually have separate state types, `wlp` must be accessible across the boundary — either via state parameter passing or by keeping it in a shared location.

7. **`qdraincomp(macp)` is not written by drainage but accumulates `qdra`:**
   `waterbalance.f90:439,449` writes `qdraincomp` from `qdra`. `swapoutput.f90:578,624,645,684,703` reads `qdraincomp` for the `.afr` file. This global sits outside the drainage home tree but its value is entirely derived from `qdra`. It is currently reset and recomputed in `waterbalance.f90` at each output step — a derived quantity that will need to be reconsidered if `qdra` moves to typed state.

8. **`swdra=2` (extended drainage / legacy `.dra` reader) not yet TOML-ported:**
   `drainage_config.f90:133–138` stub-errors when `swdra=2` with `dramet/=0`. The legacy `rddre` reader populates globals for the full extended drainage case. The TOML pipeline currently gates this out, but any migration of drainage globals must be careful not to break the legacy `swdra=2` code path that still exercises the globals directly.

9. **`geofac` collision between `drainage.basic.geofac` and `drainage.surface_runoff.geofac`:**
   The adapter writes `geofac` from `config%drain%surface_runoff%geofac` (line 431) unconditionally after writing it from `config%drain%geofac` (line 306, `ipos=5` only). The second write can silently overwrite the first. This is an existing adapter ordering bug independent of the state migration, but it touches the same global and must be noted.

---

### Phase 1 + Phase 2 lessons-learned addendum (added 2026-05-10 after Phase 2 close-out)

The original Section 3.5 categorization missed several classes of cross-subsystem reads. Phase 1 + Phase 2 caught the gaps at the integration-gate tasks (Phase 1 Task 7 — `qdrd` compute readers in WLEVBAL/WBALLEV; Phase 2 Task 5 — working-buffer reads, init-routine reads, call-site reads).

**For future subsystem-migration discoveries, Section 3.5 should categorize external readers by intent:**

1. **Output readers** — output routines that read for file I/O (e.g., `outdrf`, `outbal`, `set_values`). Most obvious; usually well-cataloged.
2. **Compute readers** — compute routines in OTHER subsystems that read this subsystem's state for their own physics (e.g., `frozencond.FrozenBounds` reading drainage's `qdrain`; `WLEVBAL` reading drainage's `qdrd`). Easy to miss; the `use variables, only: …` grep finds them.
3. **Working-buffer reads** — code paths that use the legacy global as a temporary scratch space (e.g., `swdislay` redistribution in `surfacewater.f90` using global `qdra` as a working buffer alongside its own `state%drainage%qdra` reads). These need migration to either use the typed state directly or to a subroutine-local. Hard to grep for; surface during integration-gate verification.
4. **Init-routine seed reads** — init routines that read legacy globals to populate state because, at init time, the globals carry config-set values (e.g., `drainage_init` reading `drainl`, `wetper`, `ztopdislay` to seed state). These need either a different init source (typed config arg) or to be migrated AFTER the typed config is wired.
5. **Call-site argument reads** — reads that happen because a routine takes the global by reference as an argument (e.g., `bocodrb`'s `wetper(1)` static read). These need the caller to pass the typed-state slice OR the typed-config field.

**The grep templates that catch each:**

```bash
# Output readers + compute readers (use variables clauses):
grep -rEn "use variables.*\b<owned-var>\b" src/ --include="*.f90" \
  | grep -v "src/core/variables.f90" \
  | grep -v "src/state/" \
  | grep -v "src/<home-tree>"

# Working-buffer + call-site reads (raw symbol references, post-state-write):
grep -rEn "\b<owned-var>\b" src/ --include="*.f90" \
  | grep -v "src/core/variables.f90" \
  | grep -v "src/core/initialize.f90" \
  | grep -v "src/state/" \
  | grep -v "state%" \
  | grep -v "config%"
```

The first grep finds all explicitly-imported readers. The second grep (run AFTER state-side migrations) finds any remaining bare-name references that aren't migrated yet.

**Verification sequencing:** before deleting a legacy global declaration, run BOTH greps. Any non-zero result is an unmigrated reader that must be handled in the same task as the deletion. The Phase 2 lesson: skip this verification at your peril; the integration gate's check-full failure tells you something is wrong, but the grep tells you exactly where.

---

## 9. Test surface

| Test file | What it covers |
|-----------|----------------|
| `tests/unit/config/test_drainage_config.pf` (349 LoC) | `drainage_config_t` validation: enum checks for `swdra`, `dramet`, `swdivd`, `swdislay`, `nrlevs`, `swliminf`; cross-field rules `dramet=2→swdivd=1`, `swdra=2+dramet/=0`, `altcu/=0 stub`, `swliminf=1+dramet/=3`; `drainage_surface_runoff_t` enum and range checks. |
| `tests/unit/io/toml/test_read_drainage_toml.pf` (310 LoC) | TOML reader: minimal drainage section, levels array, `surface_runoff` sub-section, external file via `drainage.file`, `swliminf` field parsing, `owltab_file` per level. |
| `tests/unit/io/toml/fixtures/` (7 fixture TOML files) | Supporting fixtures: `drainage_minimal.toml`, `drainage_levels_extended.toml`, `drainage_surface_runoff_off.toml`, `drainage_surface_runoff_full.toml`, `drainage_via_file.toml`, `drainage_external.dra.toml`, `drainage_drfil.toml`, `drainage_swliminf.toml`. |
| `tests/unit/drainage/test_surfacewater_init.pf` | Belongs to surface-water (init zeroing test); physically in drainage test dir but NOT drainage subsystem scope. |
| `tests/swap-cases/3.macroporeflow/` | Integration regression: case 3 exercises the macropore + drainage interaction (`FlMacropore`, `ZDraBas` init path, `NumLevRapDra`). |

No pFUnit tests exist for `drainage.f90` compute routines (`bocodrb`, `bocodre`, `DIVDRA`). The compute layer is covered only by integration regression cases.

---

## 10. Summary statistics

- **Total LoC across home files:** 1879 (822 + 599 + 231 + 227)
- **Owned globals (Section 2 row count):** 6 (`qdrain`, `qdra`, `qdrd`, `drainl`, `wetper`, `ztopdislay`)
- **Borrowed globals (Section 3 row count):** 52 (config-time and runtime reads)
- **External readers of owned globals (Section 3.5 distinct files):** 6 (`surfacewater.f90`, `waterbalance.f90`, `frozencond.f90`, `soilhydraulics.f90`, `solute.f90`, `config_to_variables.f90` [init-time only for wetper])
- **Entry points (Section 4 row count):** 6 call sites across 4 files (`swap.f90`, `swapoutput.f90`, `macropore.f90`, `surfacewater.f90` for `bocodre` + `divdra`)
- **Output coupling sites (Section 6 row count):** 10 (reads from `state%surfacewater` fields; drainage does not yet have its own state type)
- **Number of cross-subsystem hazards (Section 8 row count):** 9
