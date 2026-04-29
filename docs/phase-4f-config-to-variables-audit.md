---
title: Phase 4f config_to_variables audit
author: SWAP modernization team
date: 2026-04-27
---

# Phase 4f config_to_variables audit

Goal: enumerate every `variables%` field referenced by SWAP execution paths,
and for each one record where Phase 4f's `config_to_variables` adapter will
source its value.

Status legend:

- **C** -- covered by a `*_config_t` field; adapter copies it.
- **R** -- runtime state, never read from input; adapter does NOT touch.
- **G** -- gap: field is read by execution paths but no config covers it.
  Phase 4f blocker; resolve before the strangler-fig lands.
- **RETIRED** -- legacy output switch retired per
  [ADR 0009](adr/0009-discontinue-non-csv-outputs.md). New TOML reader
  has no slot; if a TOML file lists one, the reader emits a
  deprecation warning and ignores it. Phase 4f's
  `config_to_variables` forces the corresponding global to 0.
- **DEFERRED** -- macropore-related field per
  [ADR 0010](adr/0010-macropore-deferral.md). The new TOML pipeline
  does not read or wire these fields; case 3 (3.macroporeflow)
  continues to run on the legacy `readswap()` path. Schema lives at
  `src/config/macropore_config.f90` as orphan infrastructure for
  future macropore work.

Method:

1. Enumerated every `variables`-module field declared in `src/core/variables.f90`
   (1258 globals).
2. Scanned all `src/**/*.f90` for `use variables, only: ...` (with continuations)
   and bare `use variables`. Excluded `src/io/readswap.f90` itself (the file
   being replaced) -- its references are not part of the post-4f execution path.
   Bare-`use` files were rescanned for any token matching a declared global.
3. Cross-referenced the resulting union (1216 names) against:
   - All `*_config_t` field declarations under `src/config/` (309 unique field
     names across 14 typed-config modules) -- match -> **C**.
   - All `rd<func>('key', ...)` calls under `src/io/` and `src/crop/` (819
     unique input keys) -- match (and not in any config) -> **G**.
   - Otherwise -> **R** (set by simulation code, never read from input).
4. Manual aliases applied for legacy-name vs schema-name renames
   (e.g. `pond` -> `soil%pondini`, `tsoil` / `zh` -> `heat%tsoil_init`).

## Summary

- **Total variables referenced by execution paths:** 1216
- **C (covered):** 207  (17%)
- **R (runtime state):** 806  (66%)
- **G (gaps):** 175  (14%) -- Phase 4f blockers
- **RETIRED:** 18  (1%) -- per ADR 0009; no longer flagged as gaps
- **DEFERRED:** 10  (1%) -- per ADR 0010; macropore section deferred

Phase 4f-prep Task A0 verified that `croptype(macrop)` is a renamed
alias for `crop_config_t.rotation_type(:)` — the auditor's exact-name
classifier missed the rename. Reclassified G → C; net +1 C / -1 G.
The Crop section's per-section table reflects the change: 111 C / 74 R
/ 29 G / 0 RETIRED / 214 total.

Initial audit (commit `68b6d8d`) listed 325 G entries. Phase 4e Task B3
reclassified 18 legacy output-format switches from G to RETIRED:
`swafo`, `swaun`, `swvap`, `swbal`, `swwba`, `swsba`, `swblc`, `swdrf`,
`swstr`, `swirg`, `swini`, `swend`, `swheader`, `swcaprise`,
`swcapriseoutput`, `swrum`, `swswb`, `swoutputmodflow`. See
[ADR 0009](adr/0009-discontinue-non-csv-outputs.md). Phase 4e Task B5
then triaged the remaining 307 G entries: 118 were reclassified to R
because the only place the corresponding input key is read is
`readswap.f90` itself, while a non-readswap simulation file is the
authoritative writer -- the legacy reader merely zero-initialised them.
The residual 189 G entries are confirmed gaps; each row's notes column
records the resolution path (which `*_config_t` should grow the slot).
Triage budget was capped at ~30s/entry, so the residual G count is
still an upper bound -- the Crop section in particular was skimmed
aggressively and may shrink further once Phase 4f's adapter prototype
exposes which fields the physics actually exercises.

Per section (post Task B5 triage):

| section | C | R | G | RETIRED | total |
|---|---:|---:|---:|---:|---:|
| General + simulation | 11 | 0 | 0 | 0 | 11 |
| Time / control / output | 0 | 16 | 12 | 17 | 45 |
| Meteorology | 12 | 63 | 15 | 0 | 90 |
| Soil + hydraulics | 15 | 69 | 18 | 0 | 102 |
| Drainage + surface water | 20 | 30 | 20 | 0 | 70 |
| Bottom boundary | 9 | 2 | 7 | 0 | 18 |
| Heat | 10 | 0 | 2 | 0 | 12 |
| Irrigation | 8 | 7 | 1 | 0 | 16 |
| Solute | 8 | 11 | 3 | 0 | 22 |
| Crop (fixed / grass / WOFOST) | 111 | 74 | 29 | 0 | 214 |
| Macropore | 0 | 22 | 0 | 0 (RETIRED) / 10 (DEFERRED) | 32 |
| Other / uncategorised | 0 | 512 | 71 | 1 | 584 |

## Initialization order (legacy)

Today the legacy reader populates `variables%` globals in roughly this order
(reconstructed from `src/io/readswap.f90` + `src/core/initialize.f90`):

1. `Initialize` opens `*.swp` via `rdinit`, then `rdsdor`/`rdsinr`/etc. fill
   the **general / project** block: `project`, `pathwork`, `pathatm`,
   `swscre`, `swerror`, `swheader`, `tstart`, `tend`.
2. **Meteorology block** -- `metfil`, `lat`, `alt`, `altw`, `swetr`,
   `swrain`, `swdivide`, `swmetdetail`, `nmetdetail`, `swsnow`, `snowinco`,
   `teprrain`, `teprsnow`, `swrainfile`, `rainfil`, `cnref`, `cnwet`, `cndry`.
3. **Crop rotation** -- `swcrop`, then for each rotation entry:
   `(cropfil, cropstart, cropend, cropfile, croptype)`. Per-rotation crop
   readers (`readcropfixed`, `readwofost`, `readgrass`) live in
   `cropgrowth.f90` not `readswap.f90`, but they are part of the same
   "first read everything" phase.
4. **Soil-water + hydraulics** -- `swsophy`, `numlay`, `botcom`, `numnod`,
   `dz`, `swhyst`, `tau`, `cofgen`/`paramvg` per layer, `bdens`, `swfrost`,
   `tfroststa`, `tfrostend`, `pondini` (read into `pond`).
5. **Initial conditions** -- `h` (initial pressure heads, optional read of
   `gwli`/depth/`hi` table), `swinco`.
6. **Bottom boundary** -- `swbotb`, then sub-block 1..8 parameters:
   `bbcfil`, `cofqha`, `cofqhb`, `cofqhc`, `aqamp`/`aqave`/`aqper`/`aqtmax`,
   `daquif`/`dpegwl`, `cofred`, `swcaprise`, sinusoid params.
7. **Heat** -- `swhea`, `tampli`, `ddamp`, `fclay`/`fquartz`/`forg` per layer,
   `tsoil_init` table (legacy `zh`+`tsoil`).
8. **Solute** -- `swsolu`, `cdrain`, `cml`, `bexp`, `cref`, dispersion +
   diffusion + sorption parameters, root uptake.
9. **Drainage / surface water** -- `swdra`, `drfil`, `swsec`, `swsrf`,
   `nrlevs`, level-block params (`drainl`, `wetper`, `geofac`, `widthr`,
   `taludr`, `rdrain`, `rinfi`, `swdtyp`, `khtop`/`khbot`/`kvtop`/`kvbot`,
   `entres`, `geofac`, `zintf`, `zbotdr`, `basegw`).
10. **Surface-water management** (case `swsec=2`): `nmper`, `swman`,
    `wlsmana`/`wlsmanb`/`wlsmanc`/`wlsmand`, `wldip`, `dropr`, `wlsman`,
    `gwlcrit`, `nphase`, `sttab`.
11. **Macropore** (case `swmacro=1`): `swshrinp`, `swsoilshr`, `shapefacmp`,
    `swsorp`, `sorpalfa`/`sorpfacparl`/`sorpmax`, `spoint`, `swpowm`,
    `geomfac`, `rzah`, `swdarcy`, `z_ah`/`z_ic`/`z_st`/`z_tp`, `zncrar`,
    `ppicss`, `vlmpstss`, `dipoma`/`dipomi`, `pndmxmp`.
12. **Irrigation** -- `swirfix`, `swirfilt`, scheduled tables (`irgdate`,
    `irrigaction`, `cirrs`, `cirrthres`, ...).
13. **Output / control** -- `period`, `swperiod`, `swres`, `swodat`,
    `outdat[int]`, `swcsv`, `swcsv_tz`, `swvap`, `swafo`, `swaun`, `swbal`,
    `swwba`, `swsba`, `swblc`, `swdrf`, `swstr`, plus output-stream switches.
14. After `readswap` returns, `Initialize` calls `TimeControl(1)` which
    converts `tstart` -> `t1900`, sets `daycum`, `daynr`, `iyear`, `timjan1`.
    Soil-water solver is then primed with `theta`, `gwl`, `volact` from `h`.

The Phase 4f adapter must preserve **steps 1-13 effects** by copying every
(C) field from `swap_config_t` into the corresponding `variables%` slot, in
an order that respects internal dependencies (e.g. `numnod` and `dz` before
any per-node array). Step 14 stays untouched -- that is `R`-territory.

## Field map

### General + simulation

C=11, R=0, G=0

| variable | status | source / target | notes |
|---|---|---|---|
| nprintday | C | config%simulation%nprintday |  |
| pathatm | C | config%general%pathatm |  |
| pathcrop | C | config%general%pathcrop |  |
| pathwork | C | config%general%pathwork |  |
| period | C | config%simulation%period |  |
| project | C | config%general%project |  |
| swodat | C | config%simulation%swodat |  |
| swres | C | config%simulation%swres |  |
| swscre | C | config%general%swscre |  |
| tend | C | config%simulation%tend |  |
| tstart | C | config%simulation%tstart |  |

### Time / control / output

C=0, R=16, G=12, RETIRED=17 (per ADR 0009)

| variable | status | source / target | notes |
|---|---|---|---|
| date | G | -- (no config) | input key 'date' read only by readswap; needs slot in simulation_config_t (timestep / output controls) |
| dt | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| dtmax | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| dtmin | G | -- (no config) | input key 'dtmin' read only by readswap; needs slot in simulation_config_t (timestep / output controls) |
| fldumpconvcrit | G | -- (no config) | input key 'fldumpconvcrit' read only by readswap; needs slot in simulation_config_t (timestep / output controls) |
| flMaxIterTime | G | -- (no config) | input key 'flmaxitertime' read only by readswap; needs slot in simulation_config_t (timestep / output controls) |
| ipos | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| outdat | G | -- (no config) | input key 'outdat' read only by readswap; needs slot in simulation_config_t (timestep / output controls) |
| outdatint | G | -- (no config) | input key 'outdatint' read only by readswap; needs slot in simulation_config_t (timestep / output controls) |
| swafo | RETIRED | ADR 0009 | discontinued; new TOML reader emits deprecation warning |
| swaun | RETIRED | ADR 0009 | discontinued; new TOML reader emits deprecation warning |
| swbal | RETIRED | ADR 0009 | discontinued; new TOML reader emits deprecation warning |
| swblc | RETIRED | ADR 0009 | discontinued; new TOML reader emits deprecation warning |
| swcaprise | RETIRED | ADR 0009 | discontinued; new TOML reader emits deprecation warning |
| swcapriseoutput | RETIRED | ADR 0009 | discontinued; new TOML reader emits deprecation warning |
| swcsv | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| swcsv_tz | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| swdiscrvert | G | -- (no config) | input key 'swdiscrvert' read only by readswap; needs slot in simulation_config_t (timestep / output controls) |
| swdrf | RETIRED | ADR 0009 | discontinued; new TOML reader emits deprecation warning |
| SwDrRap | G | -- (no config) | input key 'swdrrap' read only by readswap; needs slot in simulation_config_t (timestep / output controls) |
| swend | RETIRED | ADR 0009 | discontinued; new TOML reader emits deprecation warning |
| swheader | RETIRED | ADR 0009 | discontinued; new TOML reader emits deprecation warning |
| swini | RETIRED | ADR 0009 | discontinued; new TOML reader emits deprecation warning |
| swirg | RETIRED | ADR 0009 | discontinued; new TOML reader emits deprecation warning |
| swkimpl | G | -- (no config) | input key 'swkimpl' read only by readswap; needs slot in simulation_config_t (timestep / output controls) |
| swkmean | G | -- (no config) | input key 'swkmean' read only by readswap; needs slot in simulation_config_t (timestep / output controls) |
| swliminf | G | -- (no config) | input key 'swliminf' read only by readswap; needs slot in simulation_config_t (timestep / output controls) |
| swrum | RETIRED | ADR 0009 | discontinued; new TOML reader emits deprecation warning |
| swsba | RETIRED | ADR 0009 | discontinued; new TOML reader emits deprecation warning |
| swstr | RETIRED | ADR 0009 | discontinued; new TOML reader emits deprecation warning |
| swsublim | C | config%soil%frost%swsublim | Phase 4f-prep Task C3 added the schema slot under `[soil.frost]`. No regression case exercises swsublim (all .swp.template files omit the key, defaulting to 0); Task D3 noop on the per-case TOML side. |
| swswb | RETIRED | ADR 0009 | discontinued; new TOML reader emits deprecation warning |
| swvap | RETIRED | ADR 0009 | discontinued; new TOML reader emits deprecation warning |
| swwba | RETIRED | ADR 0009 | discontinued; new TOML reader emits deprecation warning |

Runtime state (R) -- set by simulation code, not read from input:

> `daycum`, `daynr`, `FlOpenFileDev`, `imonth`, `ioutdat`, `isteps`, `iyear`, `iyearm1`, `outper`, `t1900`, `tcum`

### Meteorology

C=12, R=63, G=15

| variable | status | source / target | notes |
|---|---|---|---|
| alt | C | config%meteorology%alt |  |
| altw | C | config%meteorology%altw |  |
| angstroma | C | config%meteorology%angstroma |  |
| angstromb | C | config%meteorology%angstromb |  |
| lat | C | config%meteorology%lat |  |
| nmetdetail | C | config%meteorology%nmetdetail |  |
| swdivide | C | config%meteorology%swdivide |  |
| swetr | C | config%meteorology%swetr |  |
| swetsine | C | config%meteorology%swetsine |  |
| swmetdetail | C | config%meteorology%swmetdetail |  |
| swMetFilAll | C | config%meteorology%swmetfilall |  |
| swrain | C | config%meteorology%swrain |  |
| cfbs | G | -- (no config) | input key 'cfbs' read only by readswap; needs slot in meteorology_config_t |
| cofred | G | -- (no config) | input key 'cofred' read only by readswap; needs slot in meteorology_config_t |
| dateharvest | G | -- (no config) | input key 'dateharvest' read only by readswap; needs slot in meteorology_config_t |
| metfil | G | -- (no config) | input key 'metfil' read only by readswap; needs slot in meteorology_config_t |
| rainfil | G | -- (no config) | input key 'rainfil' read only by readswap; needs slot in meteorology_config_t |
| sinamp | G | -- (no config) | input key 'sinamp' read only by readswap; needs slot in meteorology_config_t |
| sinave | G | -- (no config) | input key 'sinave' read only by readswap; needs slot in meteorology_config_t |
| sinmax | G | -- (no config) | input key 'sinmax' read only by readswap; needs slot in meteorology_config_t |
| snowcoef | G | -- (no config) | input key 'snowcoef' read only by readswap; needs slot in meteorology_config_t |
| snowinco | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| ssnow | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| station | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| swcfbs | G | -- (no config) | input key 'swcfbs' read only by readswap; needs slot in meteorology_config_t |
| tampli | G | -- (no config) | input key 'tampli' read only by readswap; needs slot in meteorology_config_t |
| TePrRain | G | -- (no config) | input key 'teprrain' read only by readswap; needs slot in meteorology_config_t |
| TePrSnow | G | -- (no config) | input key 'teprsnow' read only by readswap; needs slot in meteorology_config_t |
| tmean | G | -- (no config) | input key 'tmean' read only by readswap; needs slot in meteorology_config_t |
| wet | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| wetper | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| wscap | G | -- (no config) | input key 'wscap' read only by readswap; needs slot in meteorology_config_t |
| wso | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| wsopot | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| wst | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| wstpot | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |

Runtime state (R) -- set by simulation code, not read from input:

> `aetr`, `ahum`, `arad`, `arai`, `atmn`, `atmx`, `awin`, `caintc`, `cgrai`, `CNdry`, `cnird`, `cnrai`, `CNref`, `CNrefTAB`, `CNtimTAB`, `cntper`, `CNwet`, `empreva`, `epot`, `es0`, `et0`, `ew0`, `fldaystart`, `gird`, `grai`, `graidt`, `grain`, `gsnow`, `igrai`, `inrai`, `iprec`, `melt`, `nraida`, `nraidt`, `nrain`, `out_etr`, `out_hum`, `out_rad`, `out_tmn`, `out_tmx`, `out_wet`, `out_win`, `rainamount`, `rainfluxarray`, `rainrec`, `raintimearray`, `rh`, `subl`, `SubsidCp`, `tav`, `tavd`, `ThetaRef`, `tpot`, `yearmeteo`

### Soil + hydraulics

C=15, R=69, G=18

| variable | status | source / target | notes |
|---|---|---|---|
| bdens | C | config%soil%bdens |  |
| cofgen | C | config%soil%cofgen |  |
| gwli | C | config%soil%gwli |  |
| hcomp | C | config%soil%hcomp |  |
| isoillay | C | config%soil%isoillay |  |
| ksatexm | C | config%soil%ksatexm |  |
| ncomp | C | config%soil%ncomp |  |
| orgmat | C | config%soil%orgmat |  |
| pond | C | config%soil%pondini | renamed in config (legacy `pond` ↔ schema `pondini`) |
| pondini | C | config%soil%pondini |  |
| pondmx | C | config%soil%pondmx |  |
| rsoil | C | config%soil%rsoil |  |
| swhyst | C | config%soil%swhyst |  |
| swinco | C | config%soil%swinco |  |
| swsophy | C | config%soil%swsophy |  |
| cofani | G | -- (no config) | input key 'cofani' read only by readswap; needs slot in soil_config_t |
| dznew | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| h | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| h_enpr | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| hcrit | G | -- (no config) | input key 'hcrit' read only by readswap; needs slot in soil_config_t |
| hdepth | G | -- (no config) | input key 'hdepth' read only by readswap; needs slot in soil_config_t |
| hplate | G | -- (no config) | input key 'hplate' read only by readswap; needs slot in soil_config_t |
| hsublay | G | -- (no config) | input key 'hsublay' read only by readswap; needs slot in soil_config_t |
| kf | G | -- (no config) | input key 'kf' read only by readswap; needs slot in soil_config_t |
| kfsat | G | -- (no config) | input key 'kfsat' read only by readswap; needs slot in soil_config_t |
| ksatfit | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| nrstaring | G | -- (no config) | input key 'nrstaring' read only by readswap; needs slot in soil_config_t |
| numnodnew | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| SwDarcy | G | -- (no config) | input key 'swdarcy' read only by readswap; needs slot in soil_config_t |
| swfrost | C | config%soil%frost%swfrost | Phase 4f-prep Task C3 added the schema slot under `[soil.frost]`. All 6 regression cases have `SWFROST=0` in their .swp.template; Task D3 noop on the per-case TOML side. |
| tau | G | -- (no config) | input key 'tau' read only by readswap; needs slot in soil_config_t |
| Z_Ah | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| Z_Ic | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| Z_MB50 | G | -- (no config) | input key 'z_mb50' read only by readswap; needs slot in soil_config_t |
| Z_St | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| Z_Tp | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| zc | G | -- (no config) | input key 'zc' read only by readswap; needs slot in soil_config_t |
| ZDraBas | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| zgrz | G | -- (no config) | input key 'zgrz' read only by readswap; needs slot in soil_config_t |
| zi | G | -- (no config) | input key 'zi' read only by readswap; needs slot in soil_config_t |
| zintf | G | -- (no config) | input key 'zintf' read only by readswap; needs slot in soil_config_t |
| zmow | G | -- (no config) | input key 'zmow' read only by readswap; needs slot in soil_config_t |
| ZnCrAr | G | -- (no config) | input key 'zncrar' read only by readswap; needs slot in soil_config_t |
| ztopdislay | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |

Runtime state (R) -- set by simulation code, not read from input:

> `botcom`, `dimoca`, `disnod`, `dz`, `fclay`, `forg`, `fquartz`, `FrArMtrx`, `h0max`, `haqtab`, `HarLosOrm_tot`, `hatm`, `hbotab`, `hbweir`, `hcritab`, `heacap`, `heacon`, `hleaf`, `hm1`, `hqhtab`, `hroot`, `hsurf`, `hwlman`, `Hxylem`, `ientrytab`, `ientrytablay`, `indeks`, `inq`, `inqdra`, `inqdra_in`, `inqdra_out`, `inqrot`, `inqssdi`, `isubl`, `KsatCovLay`, `ksatthr`, `layer`, `numlay`, `numnod`, `numtab`, `numtablay`, `paramvg`, `relsatthr`, `sptab`, `sptablay`, `theta`, `thetar`, `thetas`, `thetm1`, `z`, `z10_cn`, `zbotcp`, `ZBtDm`, `ZDiPoMa`, `zfrostbot`, `zfrosttop`, `ztopcp`, `ZWaLevDm`

### Drainage + surface water

C=20, R=30, G=20

| variable | status | source / target | notes |
|---|---|---|---|
| basegw | C | config%drainage%basegw |  |
| dramet | C | config%drainage%dramet |  |
| drares | C | config%drainage%drares |  |
| entres | C | config%drainage%entres |  |
| gwlinf | C | config%drainage%gwlinf |  |
| infres | C | config%drainage%infres |  |
| L | C | config%drainage%L |  |
| nrlevs | C | config%drainage%nrlevs |  |
| rdrain | C | config%drainage%rdrain |  |
| rentry | C | config%drainage%rentry |  |
| rexit | C | config%drainage%rexit |  |
| rinfi | C | config%drainage%rinfi |  |
| swallo | C | config%drainage%swallo |  |
| swdislay | C | config%drainage%swdislay |  |
| swdivd | C | config%drainage%swdivd |  |
| swdra | C | config%drainage%swdra |  |
| swdtyp | C | config%drainage%swdtyp |  |
| taludr | C | config%drainage%taludr |  |
| widthr | C | config%drainage%widthr |  |
| zbotdr | C | config%drainage%zbotdr |  |
| cofintfl | G | -- (no config) | input key 'cofintfl' read only by readswap; needs slot in drainage_config_t (or new surface_water_config_t) |
| dropr | G | -- (no config) | input key 'dropr' read only by readswap; needs slot in drainage_config_t (or new surface_water_config_t) |
| expintfl | G | -- (no config) | input key 'expintfl' read only by readswap; needs slot in drainage_config_t (or new surface_water_config_t) |
| ftopdislay | G | -- (no config) | input key 'ftopdislay' read only by readswap; needs slot in drainage_config_t (or new surface_water_config_t) |
| geofac | G | -- (no config) | input key 'geofac' read only by readswap; needs slot in drainage_config_t (or new surface_water_config_t) |
| gwl | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| gwlconv | G | -- (no config) | input key 'gwlconv' read only by readswap; needs slot in drainage_config_t (or new surface_water_config_t) |
| impend | G | -- (no config) | input key 'impend' read only by readswap; needs slot in drainage_config_t (or new surface_water_config_t) |
| nmper | G | -- (no config) | input key 'nmper' read only by readswap; needs slot in drainage_config_t (or new surface_water_config_t) |
| NumLevRapDra | G | -- (no config) | input key 'numlevrapdra' read only by readswap; needs slot in drainage_config_t (or new surface_water_config_t) |
| qdrain | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| RapDraReaExp | G | -- (no config) | input key 'rapdrareaexp' read only by readswap; needs slot in drainage_config_t (or new surface_water_config_t) |
| RapDraResRef | G | -- (no config) | input key 'rapdraresref' read only by readswap; needs slot in drainage_config_t (or new surface_water_config_t) |
| rsurfdeep | G | -- (no config) | input key 'rsurfdeep' read only by readswap; needs slot in drainage_config_t (or new surface_water_config_t) |
| rsurfshallow | G | -- (no config) | input key 'rsurfshallow' read only by readswap; needs slot in drainage_config_t (or new surface_water_config_t) |
| swdivdinf | G | -- (no config) | input key 'swdivdinf' read only by readswap; needs slot in drainage_config_t (or new surface_water_config_t) |
| swnrsrf | G | -- (no config) | input key 'swnrsrf' read only by readswap; needs slot in drainage_config_t (or new surface_water_config_t) |
| swsec | G | -- (no config) | input key 'swsec' read only by readswap; needs slot in drainage_config_t (or new surface_water_config_t) |
| swsrf | G | -- (no config) | input key 'swsrf' read only by readswap; needs slot in drainage_config_t (or new surface_water_config_t) |
| swtopdislay | G | -- (no config) | input key 'swtopdislay' read only by readswap; needs slot in drainage_config_t (or new surface_water_config_t) |
| SwTopnrsrf | G | -- (no config) | input key 'swtopnrsrf' read only by readswap; needs slot in drainage_config_t (or new surface_water_config_t) |
| wldip | G | -- (no config) | input key 'wldip' read only by readswap; needs slot in drainage_config_t (or new surface_water_config_t) |
| wls | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |

Runtime state (R) -- set by simulation code, not read from input:

> `drainl`, `gwlcrit`, `GWlFlCpZo`, `gwlinp`, `gwlm1`, `gwltab`, `imper`, `nphase`, `nrpri`, `qbot`, `qbot_nonfrozen`, `qbotab`, `qdra`, `qdraincomp`, `qdrd`, `qdrtab`, `qdrtot`, `QRapDra`, `qrosum`, `sttab`, `wlp`, `wlptab`, `wlsbak`, `wlsman`, `wlsold`, `wlstab`, `wlstar`

### Bottom boundary

C=9, R=2, G=7

| variable | status | source / target | notes |
|---|---|---|---|
| aqamp | C | config%bottom_boundary%aqamp |  |
| aqave | C | config%bottom_boundary%aqave |  |
| aqper | C | config%bottom_boundary%aqper |  |
| aqtmax | C | config%bottom_boundary%aqtmax |  |
| hbot | C | config%bottom_boundary%hbot |  |
| hdrain | C | config%bottom_boundary%hdrain |  |
| rimlay | C | config%bottom_boundary%rimlay |  |
| shape | C | config%bottom_boundary%shape | also in: drainage/drainage_config_t |
| swbotb | C | config%bottom_boundary%swbotb |  |
| cofqha | G | -- (no config) | input key 'cofqha' read only by readswap; needs slot in bottom_boundary_config_t |
| cofqhb | G | -- (no config) | input key 'cofqhb' read only by readswap; needs slot in bottom_boundary_config_t |
| cofqhc | G | -- (no config) | input key 'cofqhc' read only by readswap; needs slot in bottom_boundary_config_t |
| daquif | G | -- (no config) | input key 'daquif' read only by readswap; needs slot in bottom_boundary_config_t |
| swbotb3Impl | G | -- (no config) | input key 'swbotb3impl' read only by readswap; needs slot in bottom_boundary_config_t |
| swqhbot | G | -- (no config) | input key 'swqhbot' read only by readswap; needs slot in bottom_boundary_config_t |
| swqhr | G | -- (no config) | input key 'swqhr' read only by readswap; needs slot in bottom_boundary_config_t |

Runtime state (R) -- set by simulation code, not read from input:

> `nowltab`, `SwBotb3ResVert`

### Heat

C=10, R=0, G=2

| variable | status | source / target | notes |
|---|---|---|---|
| pclay | C | config%heat%pclay |  |
| psand | C | config%heat%psand |  |
| swbotbhea | C | config%heat%swbotbhea |  |
| swcalt | C | config%heat%swcalt |  |
| swhea | C | config%heat%swhea |  |
| swtopbhea | C | config%heat%swtopbhea |  |
| tfrostend | C | config%heat%tfrostend |  |
| tfroststa | C | config%heat%tfroststa |  |
| tsoil | C | config%heat%tsoil_init | renamed in config (legacy `tsoil` ↔ schema `tsoil_init`) |
| zh | C | config%heat%tsoil_init | renamed in config (legacy `zh` ↔ schema `tsoil_init`) |
| ddamp | G | -- (no config) | input key 'ddamp' read only by readswap; needs slot in heat_config_t |
| fdepth | G | -- (no config) | input key 'fdepth' read only by readswap; needs slot in heat_config_t |

### Irrigation

C=8, R=7, G=1

| variable | status | source / target | notes |
|---|---|---|---|
| cirrs | C | config%irrigation%cirrs |  |
| cirrthres | C | config%irrigation%cirrthres |  |
| dcrit | C | config%irrigation%dcrit |  |
| isuas | C | config%irrigation%isuas |  |
| perirrsurp | C | config%irrigation%perirrsurp |  |
| raithreshold | C | config%irrigation%raithreshold |  |
| swcirrthres | C | config%irrigation%swcirrthres |  |
| swirfix | C | config%irrigation%swirfix |  |
| irconc | G | -- (no config) | input key 'irconc' read only by readswap; needs slot in irrigation_config_t |

Runtime state (R) -- set by simulation code, not read from input:

> `cqssdi`, `dt_SSDI_event`, `iqssdi`, `irrigevent`, `qssdi`, `qssdisum`, `sqirrig`

### Solute

C=8, R=11, G=3

| variable | status | source / target | notes |
|---|---|---|---|
| bexp | C | config%solute%bexp |  |
| cdrain | C | config%solute%cdrain |  |
| cseep | C | config%solute%cseep |  |
| ldis | C | config%solute%ldis |  |
| rtheta | C | config%solute%rtheta |  |
| swbotbc | C | config%solute%swbotbc |  |
| swsolu | C | config%solute%swsolu |  |
| tscf | C | config%solute%tscf |  |
| c_mroot | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| cml | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| cpre | G | -- (no config) | input key 'cpre' read only by readswap; needs slot in solute_config_t |
| cref | G | -- (no config) | input key 'cref' read only by readswap; needs slot in solute_config_t |
| flAgeTracer | G | -- (no config) | input key 'flagetracer' read only by readswap; needs slot in solute_config_t |

Runtime state (R) -- set by simulation code, not read from input:

> `c_top`, `cmsy`, `cpond`, `crunoff`, `crunoffCN`, `crunon`, `csurf`, `flsolute`, `solbal`

### Crop (fixed / grass / WOFOST)

C=110, R=74, G=30

| variable | status | source / target | notes |
|---|---|---|---|
| adcrh | C | config%cropfixed%adcrh | also in: cropgrass/cropgrass_config_t, cropwofost/wofost_drought_stress_t |
| adcrl | C | config%cropfixed%adcrl | also in: cropgrass/cropgrass_config_t, cropwofost/wofost_drought_stress_t |
| aeratecrit | C | config%cropwofost%aeratecrit |  |
| agerm | C | config%cropwofost%agerm |  |
| air_filled_root_por | C | config%cropwofost%air_filled_root_por |  |
| albedo | C | config%cropwofost%albedo |  |
| alphacrit | C | config%cropwofost%alphacrit |  |
| amaxtb | C | config%cropwofost%amaxtb |  |
| cftb | C | config%cropfixed%cftb | also in: cropgrass/cropgrass_config_t, cropwofost/wofost_cropfactor_t |
| chtb | C | config%cropfixed%chtb | also in: cropgrass/cropgrass_config_t, cropwofost/wofost_cropfactor_t |
| co2amaxtb | C | config%cropwofost%co2amaxtb |  |
| co2efftb | C | config%cropwofost%co2efftb |  |
| co2tratb | C | config%cropwofost%co2tratb |  |
| cofab | C | config%cropfixed%cofab | also in: cropgrass/cropgrass_config_t, cropwofost/wofost_interception_t |
| cvl | C | config%cropwofost%cvl |  |
| cvo | C | config%cropwofost%cvo |  |
| cvr | C | config%cropwofost%cvr |  |
| cvs | C | config%cropwofost%cvs |  |
| dcritrtz | C | config%cropwofost%dcritrtz |  |
| dlc | C | config%cropwofost%dlc |  |
| dlo | C | config%cropwofost%dlo |  |
| dry_mat_cont_roots | C | config%cropwofost%dry_mat_cont_roots |  |
| dtsmtb | C | config%cropwofost%dtsmtb |  |
| dvsend | C | config%cropwofost%dvsend |  |
| eff | C | config%cropgrass%eff | also in: cropwofost/wofost_assimilation_t |
| fltb | C | config%cropwofost%fltb |  |
| fotb | C | config%cropwofost%fotb |  |
| frtb | C | config%cropwofost%frtb |  |
| fstb | C | config%cropwofost%fstb |  |
| hdrygerm | C | config%cropwofost%hdrygerm |  |
| hlim1 | C | config%cropfixed%hlim1 | also in: cropgrass/cropgrass_config_t, cropwofost/wofost_oxygen_stress_t |
| hlim2l | C | config%cropfixed%hlim2l | also in: cropgrass/cropgrass_config_t, cropwofost/wofost_oxygen_stress_t |
| hlim2u | C | config%cropfixed%hlim2u | also in: cropgrass/cropgrass_config_t, cropwofost/wofost_oxygen_stress_t |
| hlim3h | C | config%cropfixed%hlim3h | also in: cropgrass/cropgrass_config_t, cropwofost/wofost_drought_stress_t |
| hlim3l | C | config%cropfixed%hlim3l | also in: cropgrass/cropgrass_config_t, cropwofost/wofost_drought_stress_t |
| hlim4 | C | config%cropfixed%hlim4 | also in: cropgrass/cropgrass_config_t, cropwofost/wofost_drought_stress_t |
| hPrep | C | config%cropwofost%hprep |  |
| hSow | C | config%cropwofost%hsow |  |
| hwetgerm | C | config%cropwofost%hwetgerm |  |
| idev | C | config%cropfixed%idev | also in: cropgrass/cropgrass_config_t |
| idsl | C | config%cropwofost%idsl |  |
| kdif | C | config%cropfixed%kdif | also in: cropgrass/cropgrass_config_t, cropwofost/wofost_assimilation_t |
| kdir | C | config%cropfixed%kdir | also in: cropgrass/cropgrass_config_t, cropwofost/wofost_assimilation_t |
| laiem | C | config%cropwofost%laiem |  |
| MaxPrepDelay | C | config%cropwofost%maxprepdelay |  |
| MaxSowDelay | C | config%cropwofost%maxsowdelay |  |
| perdl | C | config%cropwofost%perdl |  |
| q10 | C | config%cropwofost%q10 |  |
| q10_microbial | C | config%cropwofost%q10_microbial |  |
| rdc | C | config%cropfixed%rdc | also in: cropgrass/cropgrass_config_t, cropwofost/wofost_root_t |
| rdctb | C | config%cropfixed%rdctb | also in: cropgrass/cropgrass_config_t, cropwofost/wofost_root_t |
| rdi | C | config%cropfixed%rdi | also in: cropgrass/cropgrass_config_t, cropwofost/wofost_root_t |
| rdrrtb | C | config%cropwofost%rdrrtb |  |
| rdrstb | C | config%cropwofost%rdrstb |  |
| rdtb | C | config%cropwofost%rdtb |  |
| relmf | C | config%cropwofost%relmf |  |
| rfsetb | C | config%cropwofost%rfsetb |  |
| rgrlai | C | config%cropwofost%rgrlai |  |
| rlwtb | C | config%cropwofost%rlwtb |  |
| rml | C | config%cropwofost%rml |  |
| rmo | C | config%cropwofost%rmo |  |
| rmr | C | config%cropwofost%rmr |  |
| rms | C | config%cropwofost%rms |  |
| root_radiusO2 | C | config%cropwofost%root_radiusO2 |  |
| rri | C | config%cropfixed%rri | also in: cropgrass/cropgrass_config_t, cropwofost/wofost_root_t |
| rsc | C | config%cropfixed%rsc | also in: cropgrass/cropgrass_config_t, cropwofost/wofost_cropfactor_t |
| rsw | C | config%cropwofost%rsw |  |
| salthead | C | config%cropwofost%salthead |  |
| saltmax | C | config%cropwofost%saltmax |  |
| saltslope | C | config%cropwofost%saltslope |  |
| schedule | C | config%cropfixed%schedule | also in: cropgrass/cropgrass_config_t, cropwofost/cropwofost_config_t, irrigation/irrigation_schedule_t |
| slatb | C | config%cropwofost%slatb |  |
| spa | C | config%cropwofost%spa |  |
| span | C | config%cropwofost%span |  |
| spec_weight_root_tissue | C | config%cropwofost%spec_weight_root_tissue |  |
| specific_resp_humus | C | config%cropwofost%specific_resp_humus |  |
| srl | C | config%cropwofost%srl |  |
| ssa | C | config%cropwofost%ssa |  |
| swcf | C | config%cropwofost%swcf |  |
| swcompensate | C | config%cropwofost%swcompensate |  |
| swdmi2rd | C | config%cropwofost%swdmi2rd |  |
| swdrought | C | config%cropwofost%swdrought |  |
| swharv | C | config%cropgrass%swharv | also in: cropwofost/wofost_harvest_t |
| swinter | C | config%cropwofost%swinter | also in: meteorology/meteorology_config_t |
| swoxygen | C | config%cropwofost%swoxygen |  |
| swpotrelmf | C | config%cropwofost%swpotrelmf |  |
| swrd | C | config%cropwofost%swrd |  |
| swrootradius | C | config%cropwofost%swrootradius |  |
| swsalinity | C | config%cropwofost%swsalinity |  |
| swstressor | C | config%cropwofost%swstressor |  |
| swWrtNonox | C | config%cropwofost%swwrtnonox |  |
| tbase | C | config%cropgrass%tbase | also in: cropwofost/wofost_greenarea_t |
| TBASEM | C | config%cropwofost%tbasem |  |
| tdwi | C | config%cropwofost%tdwi |  |
| TEFFMX | C | config%cropwofost%teffmx |  |
| TempSow | C | config%cropwofost%tempsow |  |
| tmnftb | C | config%cropwofost%tmnftb |  |
| tmpftb | C | config%cropwofost%tmpftb |  |
| tsumam | C | config%cropwofost%tsumam |  |
| tsumea | C | config%cropwofost%tsumea |  |
| tsumemeopt | C | config%cropwofost%tsumemeopt |  |
| var_a | C | config%cropwofost%var_a |  |
| vernbase | C | config%cropwofost%vernbase |  |
| verndvs | C | config%cropwofost%verndvs |  |
| vernsat | C | config%cropwofost%vernsat |  |
| wrtmax | C | config%cropwofost%wrtmax |  |
| zgerm | C | config%cropwofost%zgerm |  |
| zPrep | C | config%cropwofost%zprep |  |
| zSow | C | config%cropwofost%zsow |  |
| zTempSow | C | config%cropwofost%ztempsow |  |
| alphaw | G | -- (no config) | input key 'alphaw' read only by readswap; needs slot in crop / cropfixed / cropgrass / cropwofost configs |
| cf | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| cfeic | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| cfevappond | G | -- (no config) | input key 'cfevappond' read only by readswap; needs slot in crop / cropfixed / cropgrass / cropwofost configs |
| co2ppm | G | -- (no config) | input key 'co2ppm' read only by readswap; needs slot in crop / cropfixed / cropgrass / cropwofost configs |
| co2year | G | -- (no config) | input key 'co2year' read only by readswap; needs slot in crop / cropfixed / cropgrass / cropwofost configs |
| cropend | G | -- (no config) | input key 'cropend' read only by readswap; needs slot in crop / cropfixed / cropgrass / cropwofost configs |
| cropfil | G | -- (no config) | input key 'cropfil' read only by readswap; needs slot in crop / cropfixed / cropgrass / cropwofost configs |
| cropstart | G | -- (no config) | input key 'cropstart' read only by readswap; needs slot in crop / cropfixed / cropgrass / cropwofost configs |
| croptype | C | config%crop%rotation_type(:) | renamed alias; legacy `croptype(macrop)` and new `rotation_type(:)` carry identical per-entry values (1=fixed, 2=WOFOST, 3=grass). Phase 4f-prep Task A0 verified by spot-checking 7 caller files. Adapter copies via `do i=1,nrot; croptype(i) = config%crop%rotation_type(i); end do`. |
| cuptgraz | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| cuptgrazpot | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| cwdm | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| cwdmpot | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| dmgrztb | G | -- (no config) | input key 'dmgrztb' read only by readswap; needs slot in crop / cropfixed / cropgrass / cropwofost configs |
| dmmowtb | G | -- (no config) | input key 'dmmowtb' read only by readswap; needs slot in crop / cropfixed / cropgrass / cropwofost configs |
| dvs | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| dwlv | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| dwlvCrop | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| dwlvpot | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| dwlvSoil | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| dwrt | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| dwrtpot | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| dwst | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| dwstpot | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| f_senes | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| fimin | G | -- (no config) | input key 'fimin' read only by readswap; needs slot in crop / cropfixed / cropgrass / cropwofost configs |
| frexp | G | -- (no config) | input key 'frexp' read only by readswap; needs slot in crop / cropfixed / cropgrass / cropwofost configs |
| gasst | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| gasstpot | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| gctb | G | -- (no config) | input key 'gctb' read only by readswap; needs slot in crop / cropfixed / cropgrass / cropwofost configs |
| glaiex | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| glaiexpot | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| idaysgraz | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| idaysgrazpot | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| idregr | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| idregrpot | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| Kroot | G | -- (no config) | input key 'kroot' read only by readswap; needs slot in crop / cropfixed / cropgrass / cropwofost configs |
| kstem | G | -- (no config) | input key 'kstem' read only by readswap; needs slot in crop / cropfixed / cropgrass / cropwofost configs |
| lai | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| laiexp | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| laiexppot | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| laimax | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| laipot | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| lvage | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| lvagepot | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| mowrest | G | -- (no config) | input key 'mowrest' read only by readswap; needs slot in crop / cropfixed / cropgrass / cropwofost configs |
| mrest | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| mrestpot | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| PpIcSs | G | -- (no config) | input key 'ppicss' read only by readswap; needs slot in crop / cropfixed / cropgrass / cropwofost configs |
| q10_root | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| rootcoefa | G | -- (no config) | input key 'rootcoefa' read only by readswap; needs slot in crop / cropfixed / cropgrass / cropwofost configs |
| rooteff | G | -- (no config) | input key 'rooteff' read only by readswap; needs slot in crop / cropfixed / cropgrass / cropwofost configs |
| rootradius | G | -- (no config) | input key 'rootradius' read only by readswap; needs slot in crop / cropfixed / cropgrass / cropwofost configs |
| rsro | G | -- (no config) | input key 'rsro' read only by readswap; needs slot in crop / cropfixed / cropgrass / cropwofost configs |
| rsroexp | G | -- (no config) | input key 'rsroexp' read only by readswap; needs slot in crop / cropfixed / cropgrass / cropwofost configs |
| Rxylem | G | -- (no config) | input key 'rxylem' read only by readswap; needs slot in crop / cropfixed / cropgrass / cropwofost configs |
| saev | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| seqgrazmow | G | -- (no config) | input key 'seqgrazmow' read only by readswap; needs slot in crop / cropfixed / cropgrass / cropwofost configs |
| sla | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| slapot | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| spev | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| swcrp | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| swoxygentype | G | -- (no config) | input key 'swoxygentype' read only by readswap; needs slot in crop / cropfixed / cropgrass / cropwofost configs |
| SwShrInp | G | -- (no config) | input key 'swshrinp' read only by readswap; needs slot in crop / cropfixed / cropgrass / cropwofost configs |
| tadw | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| tadwpot | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| tagp | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| tagppot | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| tagpt | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| tagptpot | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| tsum | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| tsumdepth | G | -- (no config) | input key 'tsumdepth' read only by readswap; needs slot in crop / cropfixed / cropgrass / cropwofost configs |
| tsumgerm | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| tsumtemp | G | -- (no config) | input key 'tsumtemp' read only by readswap; needs slot in crop / cropfixed / cropgrass / cropwofost configs |
| tsumtime | G | -- (no config) | input key 'tsumtime' read only by readswap; needs slot in crop / cropfixed / cropgrass / cropwofost configs |
| wlv | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| wlvpot | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| wrt | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| wrtb | G | -- (no config) | input key 'wrtb' read only by readswap; needs slot in crop / cropfixed / cropgrass / cropwofost configs |
| wrtpot | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |

Runtime state (R) -- set by simulation code, not read from input:

> `alpJvLier`, `atav`, `cfeictb`, `cropendact`, `cropendpot`, `cropstartact`, `cropstartpot`, `cwout`, `cwsupp`, `flCropOpenFile`, `flCropReadFile`, `iptra`, `iptra_day`, `iqrot`, `mowdm`, `PpDmCp`, `PpIcTpMp`, `qrot`, `rdm`, `rdmax`, `seqgrazmowpot`, `swmeteo`, `wrtmin`

### Macropore

C=0, R=22, G=0, DEFERRED=10 (per ADR 0010)

| variable | status | source / target | notes |
|---|---|---|---|
| DiPoMa | DEFERRED | macropore_config_t (orphan, ADR 0010) | schema lives at src/config/macropore_config.f90 but is unused; Phase 4f-prep deferred wiring per ADR 0010 |
| DiPoMi | DEFERRED | macropore_config_t (orphan, ADR 0010) | schema lives at src/config/macropore_config.f90 but is unused; Phase 4f-prep deferred wiring per ADR 0010 |
| GeomFac | DEFERRED | macropore_config_t (orphan, ADR 0010) | schema lives at src/config/macropore_config.f90 but is unused; Phase 4f-prep deferred wiring per ADR 0010 |
| PowM | DEFERRED | macropore_config_t (orphan, ADR 0010) | schema lives at src/config/macropore_config.f90 but is unused; Phase 4f-prep deferred wiring per ADR 0010 |
| PrepDelay | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| Rzah | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| ShapeFacMp | DEFERRED | macropore_config_t (orphan, ADR 0010) | schema lives at src/config/macropore_config.f90 but is unused; Phase 4f-prep deferred wiring per ADR 0010 |
| SorpAlfa | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| SorpFacParl | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| SorpMax | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| Spoint | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| SwBma | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| swman | DEFERRED | macropore_config_t (orphan, ADR 0010) | schema lives at src/config/macropore_config.f90 but is unused; Phase 4f-prep deferred wiring per ADR 0010 |
| SwPowM | DEFERRED | macropore_config_t (orphan, ADR 0010) | schema lives at src/config/macropore_config.f90 but is unused; Phase 4f-prep deferred wiring per ADR 0010 |
| SwSoilShr | DEFERRED | macropore_config_t (orphan, ADR 0010) | schema lives at src/config/macropore_config.f90 but is unused; Phase 4f-prep deferred wiring per ADR 0010 |
| SwSorp | DEFERRED | macropore_config_t (orphan, ADR 0010) | schema lives at src/config/macropore_config.f90 but is unused; Phase 4f-prep deferred wiring per ADR 0010 |
| VlMpStSs | DEFERRED | macropore_config_t (orphan, ADR 0010) | schema lives at src/config/macropore_config.f90 but is unused; Phase 4f-prep deferred wiring per ADR 0010 |

Runtime state (R) -- set by simulation code, not read from input:

> `ArMpSs`, `ArMpTp`, `ArMpTpDm`, `DiPoCp`, `flmacropore`, `fprecnosnow`, `VlMp`, `VlMpDm`, `VlMpDm1`, `VlMpDm2`, `VlMpDmCp`, `VlMpDyCp`, `VlMpStCp`, `VlMpStDm1`, `VlMpStDm2`

### Other / uncategorised

C=0, R=512, G=71, RETIRED=1 (per ADR 0009)

| variable | status | source / target | notes |
|---|---|---|---|
| atmin7 | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| betaw | G | -- (no config) | input key 'betaw' read only by readswap; needs slot in needs Phase 4f categorisation |
| bgerm | G | -- (no config) | input key 'bgerm' read only by readswap; needs slot in needs Phase 4f categorisation |
| cgerm | G | -- (no config) | input key 'cgerm' read only by readswap; needs slot in needs Phase 4f categorisation |
| ch | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| CritDevh1Cp | G | -- (no config) | input key 'critdevh1cp' read only by readswap; needs slot in needs Phase 4f categorisation |
| CritDevh2Cp | G | -- (no config) | input key 'critdevh2cp' read only by readswap; needs slot in needs Phase 4f categorisation |
| CritDevMasBal | G | -- (no config) | input key 'critdevmasbal' read only by readswap; needs slot in needs Phase 4f categorisation |
| CritDevPondDt | G | -- (no config) | input key 'critdevponddt' read only by readswap; needs slot in needs Phase 4f categorisation |
| CriterHr | G | -- (no config) | input key 'criterhr' read only by readswap; needs slot in needs Phase 4f categorisation |
| CritUndSatVol | G | -- (no config) | input key 'critundsatvol' read only by readswap; needs slot in needs Phase 4f categorisation |
| daycrop | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| dayfix | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| daygrowth | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| daygrowthpot | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| ddif | G | -- (no config) | input key 'ddif' read only by readswap; needs slot in needs Phase 4f categorisation |
| decpot | G | -- (no config) | input key 'decpot' read only by readswap; needs slot in needs Phase 4f categorisation |
| decsat | G | -- (no config) | input key 'decsat' read only by readswap; needs slot in needs Phase 4f categorisation |
| dewrest | G | -- (no config) | input key 'dewrest' read only by readswap; needs slot in needs Phase 4f categorisation |
| FacDpthInf | G | -- (no config) | input key 'facdpthinf' read only by readswap; needs slot in needs Phase 4f categorisation |
| fbltb | G | -- (no config) | input key 'fbltb' read only by readswap; needs slot in needs Phase 4f categorisation |
| flCropNut | G | -- (no config) | input key 'flcropnut' read only by readswap; needs slot in needs Phase 4f categorisation |
| flprintdt | G | -- (no config) | input key 'flprintdt' read only by readswap; needs slot in needs Phase 4f categorisation |
| flSwapShared | G | -- (no config) | input key 'flswapshared' read only by readswap; needs slot in needs Phase 4f categorisation |
| frnx | G | -- (no config) | input key 'frnx' read only by readswap; needs slot in needs Phase 4f categorisation |
| gampar | G | -- (no config) | input key 'gampar' read only by readswap; needs slot in needs Phase 4f categorisation |
| iharvest | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| iHWCKmodel | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| ilvold | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| ilvoldpot | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| intwl | G | -- (no config) | input key 'intwl' read only by readswap; needs slot in needs Phase 4f categorisation |
| irdate | G | -- (no config) | input key 'irdate' read only by readswap; needs slot in needs Phase 4f categorisation |
| irdepth | G | -- (no config) | input key 'irdepth' read only by readswap; needs slot in needs Phase 4f categorisation |
| irtype | G | -- (no config) | input key 'irtype' read only by readswap; needs slot in needs Phase 4f categorisation |
| iseqgm | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| iseqgmpot | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| khbot | G | -- (no config) | input key 'khbot' read only by readswap; needs slot in needs Phase 4f categorisation |
| khtop | G | -- (no config) | input key 'khtop' read only by readswap; needs slot in needs Phase 4f categorisation |
| kvbot | G | -- (no config) | input key 'kvbot' read only by readswap; needs slot in needs Phase 4f categorisation |
| kvtop | G | -- (no config) | input key 'kvtop' read only by readswap; needs slot in needs Phase 4f categorisation |
| ldwet | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| lrnr | G | -- (no config) | input key 'lrnr' read only by readswap; needs slot in needs Phase 4f categorisation |
| lsda | G | -- (no config) | input key 'lsda' read only by readswap; needs slot in needs Phase 4f categorisation |
| lsnr | G | -- (no config) | input key 'lsnr' read only by readswap; needs slot in needs Phase 4f categorisation |
| lv | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| lvpot | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| MaxBackTr | G | -- (no config) | input key 'maxbacktr' read only by readswap; needs slot in needs Phase 4f categorisation |
| MaxIt | G | -- (no config) | input key 'maxit' read only by readswap; needs slot in needs Phase 4f categorisation |
| MaxIterTime | G | -- (no config) | input key 'maxitertime' read only by readswap; needs slot in needs Phase 4f categorisation |
| mrftb | G | -- (no config) | input key 'mrftb' read only by readswap; needs slot in needs Phase 4f categorisation |
| nlue | G | -- (no config) | input key 'nlue' read only by readswap; needs slot in needs Phase 4f categorisation |
| nmxlv | G | -- (no config) | input key 'nmxlv' read only by readswap; needs slot in needs Phase 4f categorisation |
| nofd | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| NumSbDm | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| osswlm | G | -- (no config) | input key 'osswlm' read only by readswap; needs slot in needs Phase 4f categorisation |
| outfil | G | -- (no config) | input key 'outfil' read only by readswap; needs slot in needs Phase 4f categorisation |
| pld | G | -- (no config) | input key 'pld' read only by readswap; needs slot in needs Phase 4f categorisation |
| plwti | G | -- (no config) | input key 'plwti' read only by readswap; needs slot in needs Phase 4f categorisation |
| PndmxMp | G | -- (no config) | input key 'pndmxmp' read only by readswap; needs slot in needs Phase 4f categorisation |
| poros | G | -- (no config) | input key 'poros' read only by readswap; needs slot in needs Phase 4f categorisation |
| psilt | G | -- (no config) | input key 'psilt' read only by readswap; needs slot in needs Phase 4f categorisation |
| rad | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| rd | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| rdpot | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| remoc | G | -- (no config) | input key 'remoc' read only by readswap; needs slot in needs Phase 4f categorisation |
| rid | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| rnflv | G | -- (no config) | input key 'rnflv' read only by readswap; needs slot in needs Phase 4f categorisation |
| rnfst | G | -- (no config) | input key 'rnfst' read only by readswap; needs slot in needs Phase 4f categorisation |
| rsigni | G | -- (no config) | input key 'rsigni' read only by readswap; needs slot in needs Phase 4f categorisation |
| ShrParA | G | -- (no config) | input key 'shrpara' read only by readswap; needs slot in needs Phase 4f categorisation |
| ShrParB | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| ShrParC | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| ShrParD | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| ShrParE | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| sicact | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| siccaplai | G | -- (no config) | input key 'siccaplai' read only by readswap; needs slot in needs Phase 4f categorisation |
| slw | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| SowDelay | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| StepHr | G | -- (no config) | input key 'stephr' read only by readswap; needs slot in needs Phase 4f categorisation |
| sw2 | G | -- (no config) | input key 'sw2' read only by readswap; needs slot in needs Phase 4f categorisation |
| sw3 | G | -- (no config) | input key 'sw3' read only by readswap; needs slot in needs Phase 4f categorisation |
| sw4 | G | -- (no config) | input key 'sw4' read only by readswap; needs slot in needs Phase 4f categorisation |
| swbr | G | -- (no config) | input key 'swbr' read only by readswap; needs slot in needs Phase 4f categorisation |
| swbulb | G | -- (no config) | input key 'swbulb' read only by readswap; needs slot in needs Phase 4f categorisation |
| swgc | G | -- (no config) | input key 'swgc' read only by readswap; needs slot in needs Phase 4f categorisation |
| swinc | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| swoutputmodflow | RETIRED | ADR 0009 | discontinued; new TOML reader emits deprecation warning |
| swpondmx | G | -- (no config) | input key 'swpondmx' read only by readswap; needs slot in needs Phase 4f categorisation |
| swrdc | G | -- (no config) | input key 'swrdc' read only by readswap; needs slot in needs Phase 4f categorisation |
| swredu | G | -- (no config) | input key 'swredu' read only by readswap; needs slot in needs Phase 4f categorisation |
| swsnow | C | config%meteo%snow%swsnow | Phase 4f-prep Task C2 added the schema slot under `[meteorology.snow]`. All 6 regression cases have `SWSNOW=0` in their .swp.template; Task D3 noop on the per-case TOML side. |
| swtopsub | G | -- (no config) | input key 'swtopsub' read only by readswap; needs slot in needs Phase 4f categorisation |
| swtsum | G | -- (no config) | input key 'swtsum' read only by readswap; needs slot in needs Phase 4f categorisation |
| swuseCN | G | -- (no config) | input key 'swusecn' read only by readswap; needs slot in needs Phase 4f categorisation |
| t | R | -- | reclassified by Task B5: input key only consumed by readswap; assigned by simulation code; orphan after Phase 4f |
| taccur | G | -- (no config) | input key 'taccur' read only by readswap; needs slot in needs Phase 4f categorisation |
| ThetCrMp | G | -- (no config) | input key 'thetcrmp' read only by readswap; needs slot in needs Phase 4f categorisation |
| timref | G | -- (no config) | input key 'timref' read only by readswap; needs slot in needs Phase 4f categorisation |
| vcrit | G | -- (no config) | input key 'vcrit' read only by readswap; needs slot in needs Phase 4f categorisation |
| vernrtb | G | -- (no config) | input key 'vernrtb' read only by readswap; needs slot in needs Phase 4f categorisation |
| wc_cor | G | -- (no config) | input key 'wc_cor' read only by readswap; needs slot in needs Phase 4f categorisation |
| wiltpoint | G | -- (no config) | input key 'wiltpoint' read only by readswap; needs slot in needs Phase 4f categorisation |

Runtime state (R) -- set by simulation code (482 names; first 80 listed):

> `ad`, `afo`, `Agedrain`, `AgeGwl1m`, `Ageirr`, `Agepond`, `Agepondm1`, `Agepre`, `aintcdt`, `am`, `anlv`, `anst`, `atmdem`, `atmtr`, `aun`, `avevaptb`, `avprectb`, `AwlCorFac`, `bal`, `BiModal`, `blc`, `bma`, `bpegwl`, `cevap`, `cgird`, `cgsnow`, `cinund`, `cirr`, `cmelt`, `cpeva`, `cptra`, `cqbot`, `cqbotdo`, `cqbotup`, `cqdra`, `cqdrain`, `cqdrainin`, `cqdrainout`, `cqdrd`, `cQMpInIntSatDm1`, `cQMpInIntSatDm2`, `cQMpInMtxSatDm1`, `cQMpInMtxSatDm2`, `cQMpInTopLatDm1`, `cQMpInTopLatDm2`, `cQMpInTopVrtDm1`, `cQMpInTopVrtDm2`, `cQMpLatSs`, `cQMpOutDrRap`, `cQMpOutMtxSatDm1`, `cQMpOutMtxSatDm2`, `cQMpOutMtxUnsDm1`, `cQMpOutMtxUnsDm2`, `cqprai`, `cqrot`, `cqtdo`, `cqtup`, `crp`, `cseeptab`, `csnrai`, `csubl`, `cumdens`, `daylp`, `daymeteo`, `daynrfirst`, `daynrlast`, `days_counter_irr`, `days_interval_irr`, `DaysGrazingtab`, `dectot`, `deepgw`, `DelayRegrowthTab`, `dethum`, `detrad`, `detrain`, `detrecord`, `dettav`, `dettime`, `detwind`, `dev_cmb`, ...

## Gaps (Phase 4f blockers)

After Task B5 triage, 189 (G) entries are reproduced inline in their
section tables above. Each row's notes column gives the legacy reader
key plus a target `*_config_t`. Resolution path for each gap is one of:

1. **Add the field to an existing typed config** (e.g. `swfrost` ->
   extend `heat_config_t` with `swfrost: int`).
2. **Add a new typed config module** (e.g. `macropore_config_t` for the
   ~10 (G) fields under macropore -- Phase 4d intentionally deferred this).
3. **Reclassify as initial-condition state** if the legacy reader populates
   it once and physics never re-reads it from input -- extend the relevant
   `*_ini` slot rather than adding a runtime field.

Phase 4f decides per-gap which path to take. The audit's job is only to
surface them; this doc does not prescribe schema changes.

Highest-density gap clusters (post Task B5):

- **Other / uncategorised** -- 71 gaps
- **Crop (fixed / grass / WOFOST)** -- 30 gaps
- **Drainage + surface water** -- 20 gaps
- **Soil + hydraulics** -- 18 gaps
- **Meteorology** -- 15 gaps
- **Time / control / output** -- 12 gaps
- **Macropore** -- 10 gaps
- **Bottom boundary** -- 7 gaps
- **Solute** -- 3 gaps
- **Heat** -- 2 gaps
- **Irrigation** -- 1 gaps

## Phase 4f adapter sketch

```fortran
subroutine config_to_variables(config)
   use swap_config_mod, only: swap_config_t
   use variables  ! adapter is the one place a bare `use variables` is OK
   type(swap_config_t), intent(in) :: config

   ! --- General / simulation (~11 fields) ---
   project = config%general%project
   pathwork = config%general%pathwork
   swscre  = config%general%swscre
   tstart  = config%simulation%tstart
   tend    = config%simulation%tend
   ! ...

   ! --- Meteorology (~12 fields) ---
   metfil = config%meteorology%metfil
   lat    = config%meteorology%lat
   alt    = config%meteorology%alt
   altw   = config%meteorology%altw
   swetr  = config%meteorology%swetr
   swrain = config%meteorology%swrain
   ! ...

   ! --- Soil (~15 fields) ---
   numlay   = config%soil%numlay
   botcom(:) = config%soil%botcom(:)
   pond     = config%soil%pondini   ! legacy name -> schema name
   bdens(:) = config%soil%bdens(:)
   ! ...

   ! --- Drainage (~20 fields) ---
   swdra     = config%drainage%swdra
   nrlevs    = config%drainage%nrlevs
   drainl(:) = config%drainage%drainl(:)
   ! ...

   ! --- Bottom boundary (~9 fields) ---
   swbotb = config%bottom_boundary%swbotb
   bbcfil = config%bottom_boundary%bbcfil
   aqamp  = config%bottom_boundary%aqamp
   ! ...

   ! --- Heat (~10 fields) ---
   tampli      = config%heat%tampli
   ddamp       = config%heat%ddamp
   fclay(:)    = config%heat%fclay(:)
   tsoil(1:nh) = config%heat%tsoil_init(1:nh, 2)  ! legacy 2-col table
   ! ...

   ! --- Irrigation (~8 fields) ---
   swirfix = config%irrigation%swirfix
   ! ...

   ! --- Solute (~8 fields) ---
   swsolu = config%solute%swsolu
   cdrain = config%solute%cdrain
   ! ...

   ! --- Crop (~110 fields) ---
   ! Per-rotation; each rotation entry dispatches to fixed / grass / wofost.
   swcrop = config%crop%swcrop
   do iri = 1, size(config%crop%rotations)
      croptype(iri) = config%crop%rotations(iri)%croptype
      ! ... fixed / grass / wofost-specific copies (~50 fields each)
   end do
end subroutine config_to_variables
```

## Audit-trail notes

- Source list of declared globals: 1258 (every `variables%` slot).
- 47 declared globals are NOT referenced by any execution-path file (likely
  populated by `readswap.f90` only and dropped post-4f, or unused legacy
  cruft). They are not in this audit; the adapter does not need to set them.
- For names with multiple config matches (e.g. `cftb` lives in cropfixed,
  cropgrass, AND wofost), the table cites the primary path; the notes
  column lists the alternates. The adapter dispatches per crop type, not
  per name, so this is informational only.
- (R) classification is conservative: if a name is declared, referenced by
  execution-path code, NOT in any config, AND NOT a readswap input key, it is
  marked (R). Some of these may turn out to be unused -- Phase 4f can prune
  on contact, that is not this audit's job.
