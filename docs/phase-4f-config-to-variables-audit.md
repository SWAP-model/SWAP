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
- **C (covered):** 203  (16%)
- **R (runtime state):** 688  (56%)
- **G (gaps):** 325  (26%) -- Phase 4f blockers

Per section:

| section | C | R | G | total |
|---|---:|---:|---:|---:|
| General + simulation | 11 | 0 | 0 | 11 |
| Time / control / output | 0 | 11 | 34 | 45 |
| Meteorology | 12 | 54 | 24 | 90 |
| Soil + hydraulics | 15 | 58 | 29 | 102 |
| Drainage + surface water | 20 | 27 | 23 | 70 |
| Bottom boundary | 9 | 2 | 7 | 18 |
| Heat | 10 | 0 | 2 | 12 |
| Irrigation | 8 | 7 | 1 | 16 |
| Solute | 8 | 9 | 5 | 22 |
| Crop (fixed / grass / WOFOST) | 110 | 23 | 81 | 214 |
| Macropore | 0 | 15 | 17 | 32 |
| Other / uncategorised | 0 | 482 | 102 | 584 |

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

C=0, R=11, G=34

| variable | status | source / target | notes |
|---|---|---|---|
| date | G | -- (no config) | input key 'date' read by readswap; no config covers |
| dt | G | -- (no config) | input key 'dt' read by readswap; no config covers |
| dtmax | G | -- (no config) | input key 'dtmax' read by readswap; no config covers |
| dtmin | G | -- (no config) | input key 'dtmin' read by readswap; no config covers |
| fldumpconvcrit | G | -- (no config) | input key 'fldumpconvcrit' read by readswap; no config covers |
| flMaxIterTime | G | -- (no config) | input key 'flmaxitertime' read by readswap; no config covers |
| ipos | G | -- (no config) | input key 'ipos' read by readswap; no config covers |
| outdat | G | -- (no config) | input key 'outdat' read by readswap; no config covers |
| outdatint | G | -- (no config) | input key 'outdatint' read by readswap; no config covers |
| swafo | G | -- (no config) | input key 'swafo' read by readswap; no config covers |
| swaun | G | -- (no config) | input key 'swaun' read by readswap; no config covers |
| swbal | G | -- (no config) | input key 'swbal' read by readswap; no config covers |
| swblc | G | -- (no config) | input key 'swblc' read by readswap; no config covers |
| swcaprise | G | -- (no config) | input key 'swcaprise' read by readswap; no config covers |
| swcapriseoutput | G | -- (no config) | input key 'swcapriseoutput' read by readswap; no config covers |
| swcsv | G | -- (no config) | input key 'swcsv' read by readswap; no config covers |
| swcsv_tz | G | -- (no config) | input key 'swcsv_tz' read by readswap; no config covers |
| swdiscrvert | G | -- (no config) | input key 'swdiscrvert' read by readswap; no config covers |
| swdrf | G | -- (no config) | input key 'swdrf' read by readswap; no config covers |
| SwDrRap | G | -- (no config) | input key 'swdrrap' read by readswap; no config covers |
| swend | G | -- (no config) | input key 'swend' read by readswap; no config covers |
| swheader | G | -- (no config) | input key 'swheader' read by readswap; no config covers |
| swini | G | -- (no config) | input key 'swini' read by readswap; no config covers |
| swirg | G | -- (no config) | input key 'swirg' read by readswap; no config covers |
| swkimpl | G | -- (no config) | input key 'swkimpl' read by readswap; no config covers |
| swkmean | G | -- (no config) | input key 'swkmean' read by readswap; no config covers |
| swliminf | G | -- (no config) | input key 'swliminf' read by readswap; no config covers |
| swrum | G | -- (no config) | input key 'swrum' read by readswap; no config covers |
| swsba | G | -- (no config) | input key 'swsba' read by readswap; no config covers |
| swstr | G | -- (no config) | input key 'swstr' read by readswap; no config covers |
| swsublim | G | -- (no config) | input key 'swsublim' read by readswap; no config covers |
| swswb | G | -- (no config) | input key 'swswb' read by readswap; no config covers |
| swvap | G | -- (no config) | input key 'swvap' read by readswap; no config covers |
| swwba | G | -- (no config) | input key 'swwba' read by readswap; no config covers |

Runtime state (R) -- set by simulation code, not read from input:

> `daycum`, `daynr`, `FlOpenFileDev`, `imonth`, `ioutdat`, `isteps`, `iyear`, `iyearm1`, `outper`, `t1900`, `tcum`

### Meteorology

C=12, R=54, G=24

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
| cfbs | G | -- (no config) | input key 'cfbs' read by readswap; no config covers |
| cofred | G | -- (no config) | input key 'cofred' read by readswap; no config covers |
| dateharvest | G | -- (no config) | input key 'dateharvest' read by readswap; no config covers |
| metfil | G | -- (no config) | input key 'metfil' read by readswap; no config covers |
| rainfil | G | -- (no config) | input key 'rainfil' read by readswap; no config covers |
| sinamp | G | -- (no config) | input key 'sinamp' read by readswap; no config covers |
| sinave | G | -- (no config) | input key 'sinave' read by readswap; no config covers |
| sinmax | G | -- (no config) | input key 'sinmax' read by readswap; no config covers |
| snowcoef | G | -- (no config) | input key 'snowcoef' read by readswap; no config covers |
| snowinco | G | -- (no config) | input key 'snowinco' read by readswap; no config covers |
| ssnow | G | -- (no config) | input key 'ssnow' read by readswap; no config covers |
| station | G | -- (no config) | input key 'station' read by readswap; no config covers |
| swcfbs | G | -- (no config) | input key 'swcfbs' read by readswap; no config covers |
| tampli | G | -- (no config) | input key 'tampli' read by readswap; no config covers |
| TePrRain | G | -- (no config) | input key 'teprrain' read by readswap; no config covers |
| TePrSnow | G | -- (no config) | input key 'teprsnow' read by readswap; no config covers |
| tmean | G | -- (no config) | input key 'tmean' read by readswap; no config covers |
| wet | G | -- (no config) | input key 'wet' read by readswap; no config covers |
| wetper | G | -- (no config) | input key 'wetper' read by readswap; no config covers |
| wscap | G | -- (no config) | input key 'wscap' read by readswap; no config covers |
| wso | G | -- (no config) | input key 'wso' read by readswap; no config covers |
| wsopot | G | -- (no config) | input key 'wsopot' read by readswap; no config covers |
| wst | G | -- (no config) | input key 'wst' read by readswap; no config covers |
| wstpot | G | -- (no config) | input key 'wstpot' read by readswap; no config covers |

Runtime state (R) -- set by simulation code, not read from input:

> `aetr`, `ahum`, `arad`, `arai`, `atmn`, `atmx`, `awin`, `caintc`, `cgrai`, `CNdry`, `cnird`, `cnrai`, `CNref`, `CNrefTAB`, `CNtimTAB`, `cntper`, `CNwet`, `empreva`, `epot`, `es0`, `et0`, `ew0`, `fldaystart`, `gird`, `grai`, `graidt`, `grain`, `gsnow`, `igrai`, `inrai`, `iprec`, `melt`, `nraida`, `nraidt`, `nrain`, `out_etr`, `out_hum`, `out_rad`, `out_tmn`, `out_tmx`, `out_wet`, `out_win`, `rainamount`, `rainfluxarray`, `rainrec`, `raintimearray`, `rh`, `subl`, `SubsidCp`, `tav`, `tavd`, `ThetaRef`, `tpot`, `yearmeteo`

### Soil + hydraulics

C=15, R=58, G=29

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
| cofani | G | -- (no config) | input key 'cofani' read by readswap; no config covers |
| dznew | G | -- (no config) | input key 'dznew' read by readswap; no config covers |
| h | G | -- (no config) | input key 'h' read by readswap; no config covers |
| h_enpr | G | -- (no config) | input key 'h_enpr' read by readswap; no config covers |
| hcrit | G | -- (no config) | input key 'hcrit' read by readswap; no config covers |
| hdepth | G | -- (no config) | input key 'hdepth' read by readswap; no config covers |
| hplate | G | -- (no config) | input key 'hplate' read by readswap; no config covers |
| hsublay | G | -- (no config) | input key 'hsublay' read by readswap; no config covers |
| kf | G | -- (no config) | input key 'kf' read by readswap; no config covers |
| kfsat | G | -- (no config) | input key 'kfsat' read by readswap; no config covers |
| ksatfit | G | -- (no config) | input key 'ksatfit' read by readswap; no config covers |
| nrstaring | G | -- (no config) | input key 'nrstaring' read by readswap; no config covers |
| numnodnew | G | -- (no config) | input key 'numnodnew' read by readswap; no config covers |
| SwDarcy | G | -- (no config) | input key 'swdarcy' read by readswap; no config covers |
| swfrost | G | -- (no config) | input key 'swfrost' read by readswap; no config covers |
| tau | G | -- (no config) | input key 'tau' read by readswap; no config covers |
| Z_Ah | G | -- (no config) | input key 'z_ah' read by readswap; no config covers |
| Z_Ic | G | -- (no config) | input key 'z_ic' read by readswap; no config covers |
| Z_MB50 | G | -- (no config) | input key 'z_mb50' read by readswap; no config covers |
| Z_St | G | -- (no config) | input key 'z_st' read by readswap; no config covers |
| Z_Tp | G | -- (no config) | input key 'z_tp' read by readswap; no config covers |
| zc | G | -- (no config) | input key 'zc' read by readswap; no config covers |
| ZDraBas | G | -- (no config) | input key 'zdrabas' read by readswap; no config covers |
| zgrz | G | -- (no config) | input key 'zgrz' read by readswap; no config covers |
| zi | G | -- (no config) | input key 'zi' read by readswap; no config covers |
| zintf | G | -- (no config) | input key 'zintf' read by readswap; no config covers |
| zmow | G | -- (no config) | input key 'zmow' read by readswap; no config covers |
| ZnCrAr | G | -- (no config) | input key 'zncrar' read by readswap; no config covers |
| ztopdislay | G | -- (no config) | input key 'ztopdislay' read by readswap; no config covers |

Runtime state (R) -- set by simulation code, not read from input:

> `botcom`, `dimoca`, `disnod`, `dz`, `fclay`, `forg`, `fquartz`, `FrArMtrx`, `h0max`, `haqtab`, `HarLosOrm_tot`, `hatm`, `hbotab`, `hbweir`, `hcritab`, `heacap`, `heacon`, `hleaf`, `hm1`, `hqhtab`, `hroot`, `hsurf`, `hwlman`, `Hxylem`, `ientrytab`, `ientrytablay`, `indeks`, `inq`, `inqdra`, `inqdra_in`, `inqdra_out`, `inqrot`, `inqssdi`, `isubl`, `KsatCovLay`, `ksatthr`, `layer`, `numlay`, `numnod`, `numtab`, `numtablay`, `paramvg`, `relsatthr`, `sptab`, `sptablay`, `theta`, `thetar`, `thetas`, `thetm1`, `z`, `z10_cn`, `zbotcp`, `ZBtDm`, `ZDiPoMa`, `zfrostbot`, `zfrosttop`, `ztopcp`, `ZWaLevDm`

### Drainage + surface water

C=20, R=27, G=23

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
| cofintfl | G | -- (no config) | input key 'cofintfl' read by readswap; no config covers |
| dropr | G | -- (no config) | input key 'dropr' read by readswap; no config covers |
| expintfl | G | -- (no config) | input key 'expintfl' read by readswap; no config covers |
| ftopdislay | G | -- (no config) | input key 'ftopdislay' read by readswap; no config covers |
| geofac | G | -- (no config) | input key 'geofac' read by readswap; no config covers |
| gwl | G | -- (no config) | input key 'gwl' read by readswap; no config covers |
| gwlconv | G | -- (no config) | input key 'gwlconv' read by readswap; no config covers |
| impend | G | -- (no config) | input key 'impend' read by readswap; no config covers |
| nmper | G | -- (no config) | input key 'nmper' read by readswap; no config covers |
| NumLevRapDra | G | -- (no config) | input key 'numlevrapdra' read by readswap; no config covers |
| qdrain | G | -- (no config) | input key 'qdrain' read by readswap; no config covers |
| RapDraReaExp | G | -- (no config) | input key 'rapdrareaexp' read by readswap; no config covers |
| RapDraResRef | G | -- (no config) | input key 'rapdraresref' read by readswap; no config covers |
| rsurfdeep | G | -- (no config) | input key 'rsurfdeep' read by readswap; no config covers |
| rsurfshallow | G | -- (no config) | input key 'rsurfshallow' read by readswap; no config covers |
| swdivdinf | G | -- (no config) | input key 'swdivdinf' read by readswap; no config covers |
| swnrsrf | G | -- (no config) | input key 'swnrsrf' read by readswap; no config covers |
| swsec | G | -- (no config) | input key 'swsec' read by readswap; no config covers |
| swsrf | G | -- (no config) | input key 'swsrf' read by readswap; no config covers |
| swtopdislay | G | -- (no config) | input key 'swtopdislay' read by readswap; no config covers |
| SwTopnrsrf | G | -- (no config) | input key 'swtopnrsrf' read by readswap; no config covers |
| wldip | G | -- (no config) | input key 'wldip' read by readswap; no config covers |
| wls | G | -- (no config) | input key 'wls' read by readswap; no config covers |

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
| cofqha | G | -- (no config) | input key 'cofqha' read by readswap; no config covers |
| cofqhb | G | -- (no config) | input key 'cofqhb' read by readswap; no config covers |
| cofqhc | G | -- (no config) | input key 'cofqhc' read by readswap; no config covers |
| daquif | G | -- (no config) | input key 'daquif' read by readswap; no config covers |
| swbotb3Impl | G | -- (no config) | input key 'swbotb3impl' read by readswap; no config covers |
| swqhbot | G | -- (no config) | input key 'swqhbot' read by readswap; no config covers |
| swqhr | G | -- (no config) | input key 'swqhr' read by readswap; no config covers |

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
| ddamp | G | -- (no config) | input key 'ddamp' read by readswap; no config covers |
| fdepth | G | -- (no config) | input key 'fdepth' read by readswap; no config covers |

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
| irconc | G | -- (no config) | input key 'irconc' read by readswap; no config covers |

Runtime state (R) -- set by simulation code, not read from input:

> `cqssdi`, `dt_SSDI_event`, `iqssdi`, `irrigevent`, `qssdi`, `qssdisum`, `sqirrig`

### Solute

C=8, R=9, G=5

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
| c_mroot | G | -- (no config) | input key 'c_mroot' read by readswap; no config covers |
| cml | G | -- (no config) | input key 'cml' read by readswap; no config covers |
| cpre | G | -- (no config) | input key 'cpre' read by readswap; no config covers |
| cref | G | -- (no config) | input key 'cref' read by readswap; no config covers |
| flAgeTracer | G | -- (no config) | input key 'flagetracer' read by readswap; no config covers |

Runtime state (R) -- set by simulation code, not read from input:

> `c_top`, `cmsy`, `cpond`, `crunoff`, `crunoffCN`, `crunon`, `csurf`, `flsolute`, `solbal`

### Crop (fixed / grass / WOFOST)

C=110, R=23, G=81

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
| alphaw | G | -- (no config) | input key 'alphaw' read by readswap; no config covers |
| cf | G | -- (no config) | input key 'cf' read by readswap; no config covers |
| cfeic | G | -- (no config) | input key 'cfeic' read by readswap; no config covers |
| cfevappond | G | -- (no config) | input key 'cfevappond' read by readswap; no config covers |
| co2ppm | G | -- (no config) | input key 'co2ppm' read by readswap; no config covers |
| co2year | G | -- (no config) | input key 'co2year' read by readswap; no config covers |
| cropend | G | -- (no config) | input key 'cropend' read by readswap; no config covers |
| cropfil | G | -- (no config) | input key 'cropfil' read by readswap; no config covers |
| cropstart | G | -- (no config) | input key 'cropstart' read by readswap; no config covers |
| croptype | G | -- (no config) | input key 'croptype' read by readswap; no config covers |
| cuptgraz | G | -- (no config) | input key 'cuptgraz' read by readswap; no config covers |
| cuptgrazpot | G | -- (no config) | input key 'cuptgrazpot' read by readswap; no config covers |
| cwdm | G | -- (no config) | input key 'cwdm' read by readswap; no config covers |
| cwdmpot | G | -- (no config) | input key 'cwdmpot' read by readswap; no config covers |
| dmgrztb | G | -- (no config) | input key 'dmgrztb' read by readswap; no config covers |
| dmmowtb | G | -- (no config) | input key 'dmmowtb' read by readswap; no config covers |
| dvs | G | -- (no config) | input key 'dvs' read by readswap; no config covers |
| dwlv | G | -- (no config) | input key 'dwlv' read by readswap; no config covers |
| dwlvCrop | G | -- (no config) | input key 'dwlvcrop' read by readswap; no config covers |
| dwlvpot | G | -- (no config) | input key 'dwlvpot' read by readswap; no config covers |
| dwlvSoil | G | -- (no config) | input key 'dwlvsoil' read by readswap; no config covers |
| dwrt | G | -- (no config) | input key 'dwrt' read by readswap; no config covers |
| dwrtpot | G | -- (no config) | input key 'dwrtpot' read by readswap; no config covers |
| dwst | G | -- (no config) | input key 'dwst' read by readswap; no config covers |
| dwstpot | G | -- (no config) | input key 'dwstpot' read by readswap; no config covers |
| f_senes | G | -- (no config) | input key 'f_senes' read by readswap; no config covers |
| fimin | G | -- (no config) | input key 'fimin' read by readswap; no config covers |
| frexp | G | -- (no config) | input key 'frexp' read by readswap; no config covers |
| gasst | G | -- (no config) | input key 'gasst' read by readswap; no config covers |
| gasstpot | G | -- (no config) | input key 'gasstpot' read by readswap; no config covers |
| gctb | G | -- (no config) | input key 'gctb' read by readswap; no config covers |
| glaiex | G | -- (no config) | input key 'glaiex' read by readswap; no config covers |
| glaiexpot | G | -- (no config) | input key 'glaiexpot' read by readswap; no config covers |
| idaysgraz | G | -- (no config) | input key 'idaysgraz' read by readswap; no config covers |
| idaysgrazpot | G | -- (no config) | input key 'idaysgrazpot' read by readswap; no config covers |
| idregr | G | -- (no config) | input key 'idregr' read by readswap; no config covers |
| idregrpot | G | -- (no config) | input key 'idregrpot' read by readswap; no config covers |
| Kroot | G | -- (no config) | input key 'kroot' read by readswap; no config covers |
| kstem | G | -- (no config) | input key 'kstem' read by readswap; no config covers |
| lai | G | -- (no config) | input key 'lai' read by readswap; no config covers |
| laiexp | G | -- (no config) | input key 'laiexp' read by readswap; no config covers |
| laiexppot | G | -- (no config) | input key 'laiexppot' read by readswap; no config covers |
| laimax | G | -- (no config) | input key 'laimax' read by readswap; no config covers |
| laipot | G | -- (no config) | input key 'laipot' read by readswap; no config covers |
| lvage | G | -- (no config) | input key 'lvage' read by readswap; no config covers |
| lvagepot | G | -- (no config) | input key 'lvagepot' read by readswap; no config covers |
| mowrest | G | -- (no config) | input key 'mowrest' read by readswap; no config covers |
| mrest | G | -- (no config) | input key 'mrest' read by readswap; no config covers |
| mrestpot | G | -- (no config) | input key 'mrestpot' read by readswap; no config covers |
| PpIcSs | G | -- (no config) | input key 'ppicss' read by readswap; no config covers |
| q10_root | G | -- (no config) | input key 'q10_root' read by readswap; no config covers |
| rootcoefa | G | -- (no config) | input key 'rootcoefa' read by readswap; no config covers |
| rooteff | G | -- (no config) | input key 'rooteff' read by readswap; no config covers |
| rootradius | G | -- (no config) | input key 'rootradius' read by readswap; no config covers |
| rsro | G | -- (no config) | input key 'rsro' read by readswap; no config covers |
| rsroexp | G | -- (no config) | input key 'rsroexp' read by readswap; no config covers |
| Rxylem | G | -- (no config) | input key 'rxylem' read by readswap; no config covers |
| saev | G | -- (no config) | input key 'saev' read by readswap; no config covers |
| seqgrazmow | G | -- (no config) | input key 'seqgrazmow' read by readswap; no config covers |
| sla | G | -- (no config) | input key 'sla' read by readswap; no config covers |
| slapot | G | -- (no config) | input key 'slapot' read by readswap; no config covers |
| spev | G | -- (no config) | input key 'spev' read by readswap; no config covers |
| swcrp | G | -- (no config) | input key 'swcrp' read by readswap; no config covers |
| swoxygentype | G | -- (no config) | input key 'swoxygentype' read by readswap; no config covers |
| SwShrInp | G | -- (no config) | input key 'swshrinp' read by readswap; no config covers |
| tadw | G | -- (no config) | input key 'tadw' read by readswap; no config covers |
| tadwpot | G | -- (no config) | input key 'tadwpot' read by readswap; no config covers |
| tagp | G | -- (no config) | input key 'tagp' read by readswap; no config covers |
| tagppot | G | -- (no config) | input key 'tagppot' read by readswap; no config covers |
| tagpt | G | -- (no config) | input key 'tagpt' read by readswap; no config covers |
| tagptpot | G | -- (no config) | input key 'tagptpot' read by readswap; no config covers |
| tsum | G | -- (no config) | input key 'tsum' read by readswap; no config covers |
| tsumdepth | G | -- (no config) | input key 'tsumdepth' read by readswap; no config covers |
| tsumgerm | G | -- (no config) | input key 'tsumgerm' read by readswap; no config covers |
| tsumtemp | G | -- (no config) | input key 'tsumtemp' read by readswap; no config covers |
| tsumtime | G | -- (no config) | input key 'tsumtime' read by readswap; no config covers |
| wlv | G | -- (no config) | input key 'wlv' read by readswap; no config covers |
| wlvpot | G | -- (no config) | input key 'wlvpot' read by readswap; no config covers |
| wrt | G | -- (no config) | input key 'wrt' read by readswap; no config covers |
| wrtb | G | -- (no config) | input key 'wrtb' read by readswap; no config covers |
| wrtpot | G | -- (no config) | input key 'wrtpot' read by readswap; no config covers |

Runtime state (R) -- set by simulation code, not read from input:

> `alpJvLier`, `atav`, `cfeictb`, `cropendact`, `cropendpot`, `cropstartact`, `cropstartpot`, `cwout`, `cwsupp`, `flCropOpenFile`, `flCropReadFile`, `iptra`, `iptra_day`, `iqrot`, `mowdm`, `PpDmCp`, `PpIcTpMp`, `qrot`, `rdm`, `rdmax`, `seqgrazmowpot`, `swmeteo`, `wrtmin`

### Macropore

C=0, R=15, G=17

| variable | status | source / target | notes |
|---|---|---|---|
| DiPoMa | G | -- (no config) | input key 'dipoma' read by readswap; no config covers |
| DiPoMi | G | -- (no config) | input key 'dipomi' read by readswap; no config covers |
| GeomFac | G | -- (no config) | input key 'geomfac' read by readswap; no config covers |
| PowM | G | -- (no config) | input key 'powm' read by readswap; no config covers |
| PrepDelay | G | -- (no config) | input key 'prepdelay' read by readswap; no config covers |
| Rzah | G | -- (no config) | input key 'rzah' read by readswap; no config covers |
| ShapeFacMp | G | -- (no config) | input key 'shapefacmp' read by readswap; no config covers |
| SorpAlfa | G | -- (no config) | input key 'sorpalfa' read by readswap; no config covers |
| SorpFacParl | G | -- (no config) | input key 'sorpfacparl' read by readswap; no config covers |
| SorpMax | G | -- (no config) | input key 'sorpmax' read by readswap; no config covers |
| Spoint | G | -- (no config) | input key 'spoint' read by readswap; no config covers |
| SwBma | G | -- (no config) | input key 'swbma' read by readswap; no config covers |
| swman | G | -- (no config) | input key 'swman' read by readswap; no config covers |
| SwPowM | G | -- (no config) | input key 'swpowm' read by readswap; no config covers |
| SwSoilShr | G | -- (no config) | input key 'swsoilshr' read by readswap; no config covers |
| SwSorp | G | -- (no config) | input key 'swsorp' read by readswap; no config covers |
| VlMpStSs | G | -- (no config) | input key 'vlmpstss' read by readswap; no config covers |

Runtime state (R) -- set by simulation code, not read from input:

> `ArMpSs`, `ArMpTp`, `ArMpTpDm`, `DiPoCp`, `flmacropore`, `fprecnosnow`, `VlMp`, `VlMpDm`, `VlMpDm1`, `VlMpDm2`, `VlMpDmCp`, `VlMpDyCp`, `VlMpStCp`, `VlMpStDm1`, `VlMpStDm2`

### Other / uncategorised

C=0, R=482, G=102

| variable | status | source / target | notes |
|---|---|---|---|
| atmin7 | G | -- (no config) | input key 'atmin7' read by readswap; no config covers |
| betaw | G | -- (no config) | input key 'betaw' read by readswap; no config covers |
| bgerm | G | -- (no config) | input key 'bgerm' read by readswap; no config covers |
| cgerm | G | -- (no config) | input key 'cgerm' read by readswap; no config covers |
| ch | G | -- (no config) | input key 'ch' read by readswap; no config covers |
| CritDevh1Cp | G | -- (no config) | input key 'critdevh1cp' read by readswap; no config covers |
| CritDevh2Cp | G | -- (no config) | input key 'critdevh2cp' read by readswap; no config covers |
| CritDevMasBal | G | -- (no config) | input key 'critdevmasbal' read by readswap; no config covers |
| CritDevPondDt | G | -- (no config) | input key 'critdevponddt' read by readswap; no config covers |
| CriterHr | G | -- (no config) | input key 'criterhr' read by readswap; no config covers |
| CritUndSatVol | G | -- (no config) | input key 'critundsatvol' read by readswap; no config covers |
| daycrop | G | -- (no config) | input key 'daycrop' read by readswap; no config covers |
| dayfix | G | -- (no config) | input key 'dayfix' read by readswap; no config covers |
| daygrowth | G | -- (no config) | input key 'daygrowth' read by readswap; no config covers |
| daygrowthpot | G | -- (no config) | input key 'daygrowthpot' read by readswap; no config covers |
| ddif | G | -- (no config) | input key 'ddif' read by readswap; no config covers |
| decpot | G | -- (no config) | input key 'decpot' read by readswap; no config covers |
| decsat | G | -- (no config) | input key 'decsat' read by readswap; no config covers |
| dewrest | G | -- (no config) | input key 'dewrest' read by readswap; no config covers |
| FacDpthInf | G | -- (no config) | input key 'facdpthinf' read by readswap; no config covers |
| fbltb | G | -- (no config) | input key 'fbltb' read by readswap; no config covers |
| flCropNut | G | -- (no config) | input key 'flcropnut' read by readswap; no config covers |
| flprintdt | G | -- (no config) | input key 'flprintdt' read by readswap; no config covers |
| flSwapShared | G | -- (no config) | input key 'flswapshared' read by readswap; no config covers |
| frnx | G | -- (no config) | input key 'frnx' read by readswap; no config covers |
| gampar | G | -- (no config) | input key 'gampar' read by readswap; no config covers |
| iharvest | G | -- (no config) | input key 'iharvest' read by readswap; no config covers |
| iHWCKmodel | G | -- (no config) | input key 'ihwckmodel' read by readswap; no config covers |
| ilvold | G | -- (no config) | input key 'ilvold' read by readswap; no config covers |
| ilvoldpot | G | -- (no config) | input key 'ilvoldpot' read by readswap; no config covers |
| intwl | G | -- (no config) | input key 'intwl' read by readswap; no config covers |
| irdate | G | -- (no config) | input key 'irdate' read by readswap; no config covers |
| irdepth | G | -- (no config) | input key 'irdepth' read by readswap; no config covers |
| irtype | G | -- (no config) | input key 'irtype' read by readswap; no config covers |
| iseqgm | G | -- (no config) | input key 'iseqgm' read by readswap; no config covers |
| iseqgmpot | G | -- (no config) | input key 'iseqgmpot' read by readswap; no config covers |
| khbot | G | -- (no config) | input key 'khbot' read by readswap; no config covers |
| khtop | G | -- (no config) | input key 'khtop' read by readswap; no config covers |
| kvbot | G | -- (no config) | input key 'kvbot' read by readswap; no config covers |
| kvtop | G | -- (no config) | input key 'kvtop' read by readswap; no config covers |
| ldwet | G | -- (no config) | input key 'ldwet' read by readswap; no config covers |
| lrnr | G | -- (no config) | input key 'lrnr' read by readswap; no config covers |
| lsda | G | -- (no config) | input key 'lsda' read by readswap; no config covers |
| lsnr | G | -- (no config) | input key 'lsnr' read by readswap; no config covers |
| lv | G | -- (no config) | input key 'lv' read by readswap; no config covers |
| lvpot | G | -- (no config) | input key 'lvpot' read by readswap; no config covers |
| MaxBackTr | G | -- (no config) | input key 'maxbacktr' read by readswap; no config covers |
| MaxIt | G | -- (no config) | input key 'maxit' read by readswap; no config covers |
| MaxIterTime | G | -- (no config) | input key 'maxitertime' read by readswap; no config covers |
| mrftb | G | -- (no config) | input key 'mrftb' read by readswap; no config covers |
| nlue | G | -- (no config) | input key 'nlue' read by readswap; no config covers |
| nmxlv | G | -- (no config) | input key 'nmxlv' read by readswap; no config covers |
| nofd | G | -- (no config) | input key 'nofd' read by readswap; no config covers |
| NumSbDm | G | -- (no config) | input key 'numsbdm' read by readswap; no config covers |
| osswlm | G | -- (no config) | input key 'osswlm' read by readswap; no config covers |
| outfil | G | -- (no config) | input key 'outfil' read by readswap; no config covers |
| pld | G | -- (no config) | input key 'pld' read by readswap; no config covers |
| plwti | G | -- (no config) | input key 'plwti' read by readswap; no config covers |
| PndmxMp | G | -- (no config) | input key 'pndmxmp' read by readswap; no config covers |
| poros | G | -- (no config) | input key 'poros' read by readswap; no config covers |
| psilt | G | -- (no config) | input key 'psilt' read by readswap; no config covers |
| rad | G | -- (no config) | input key 'rad' read by readswap; no config covers |
| rd | G | -- (no config) | input key 'rd' read by readswap; no config covers |
| rdpot | G | -- (no config) | input key 'rdpot' read by readswap; no config covers |
| remoc | G | -- (no config) | input key 'remoc' read by readswap; no config covers |
| rid | G | -- (no config) | input key 'rid' read by readswap; no config covers |
| rnflv | G | -- (no config) | input key 'rnflv' read by readswap; no config covers |
| rnfst | G | -- (no config) | input key 'rnfst' read by readswap; no config covers |
| rsigni | G | -- (no config) | input key 'rsigni' read by readswap; no config covers |
| ShrParA | G | -- (no config) | input key 'shrpara' read by readswap; no config covers |
| ShrParB | G | -- (no config) | input key 'shrparb' read by readswap; no config covers |
| ShrParC | G | -- (no config) | input key 'shrparc' read by readswap; no config covers |
| ShrParD | G | -- (no config) | input key 'shrpard' read by readswap; no config covers |
| ShrParE | G | -- (no config) | input key 'shrpare' read by readswap; no config covers |
| sicact | G | -- (no config) | input key 'sicact' read by readswap; no config covers |
| siccaplai | G | -- (no config) | input key 'siccaplai' read by readswap; no config covers |
| slw | G | -- (no config) | input key 'slw' read by readswap; no config covers |
| SowDelay | G | -- (no config) | input key 'sowdelay' read by readswap; no config covers |
| StepHr | G | -- (no config) | input key 'stephr' read by readswap; no config covers |
| sw2 | G | -- (no config) | input key 'sw2' read by readswap; no config covers |
| sw3 | G | -- (no config) | input key 'sw3' read by readswap; no config covers |
| sw4 | G | -- (no config) | input key 'sw4' read by readswap; no config covers |
| swbr | G | -- (no config) | input key 'swbr' read by readswap; no config covers |
| swbulb | G | -- (no config) | input key 'swbulb' read by readswap; no config covers |
| swgc | G | -- (no config) | input key 'swgc' read by readswap; no config covers |
| swinc | G | -- (no config) | input key 'swinc' read by readswap; no config covers |
| swoutputmodflow | G | -- (no config) | input key 'swoutputmodflow' read by readswap; no config covers |
| swpondmx | G | -- (no config) | input key 'swpondmx' read by readswap; no config covers |
| swrdc | G | -- (no config) | input key 'swrdc' read by readswap; no config covers |
| swredu | G | -- (no config) | input key 'swredu' read by readswap; no config covers |
| swsnow | G | -- (no config) | input key 'swsnow' read by readswap; no config covers |
| swtopsub | G | -- (no config) | input key 'swtopsub' read by readswap; no config covers |
| swtsum | G | -- (no config) | input key 'swtsum' read by readswap; no config covers |
| swuseCN | G | -- (no config) | input key 'swusecn' read by readswap; no config covers |
| t | G | -- (no config) | input key 't' read by readswap; no config covers |
| taccur | G | -- (no config) | input key 'taccur' read by readswap; no config covers |
| ThetCrMp | G | -- (no config) | input key 'thetcrmp' read by readswap; no config covers |
| timref | G | -- (no config) | input key 'timref' read by readswap; no config covers |
| vcrit | G | -- (no config) | input key 'vcrit' read by readswap; no config covers |
| vernrtb | G | -- (no config) | input key 'vernrtb' read by readswap; no config covers |
| wc_cor | G | -- (no config) | input key 'wc_cor' read by readswap; no config covers |
| wiltpoint | G | -- (no config) | input key 'wiltpoint' read by readswap; no config covers |

Runtime state (R) -- set by simulation code (482 names; first 80 listed):

> `ad`, `afo`, `Agedrain`, `AgeGwl1m`, `Ageirr`, `Agepond`, `Agepondm1`, `Agepre`, `aintcdt`, `am`, `anlv`, `anst`, `atmdem`, `atmtr`, `aun`, `avevaptb`, `avprectb`, `AwlCorFac`, `bal`, `BiModal`, `blc`, `bma`, `bpegwl`, `cevap`, `cgird`, `cgsnow`, `cinund`, `cirr`, `cmelt`, `cpeva`, `cptra`, `cqbot`, `cqbotdo`, `cqbotup`, `cqdra`, `cqdrain`, `cqdrainin`, `cqdrainout`, `cqdrd`, `cQMpInIntSatDm1`, `cQMpInIntSatDm2`, `cQMpInMtxSatDm1`, `cQMpInMtxSatDm2`, `cQMpInTopLatDm1`, `cQMpInTopLatDm2`, `cQMpInTopVrtDm1`, `cQMpInTopVrtDm2`, `cQMpLatSs`, `cQMpOutDrRap`, `cQMpOutMtxSatDm1`, `cQMpOutMtxSatDm2`, `cQMpOutMtxUnsDm1`, `cQMpOutMtxUnsDm2`, `cqprai`, `cqrot`, `cqtdo`, `cqtup`, `crp`, `cseeptab`, `csnrai`, `csubl`, `cumdens`, `daylp`, `daymeteo`, `daynrfirst`, `daynrlast`, `days_counter_irr`, `days_interval_irr`, `DaysGrazingtab`, `dectot`, `deepgw`, `DelayRegrowthTab`, `dethum`, `detrad`, `detrain`, `detrecord`, `dettav`, `dettime`, `detwind`, `dev_cmb`, ...

## Gaps (Phase 4f blockers)

All 325 (G) entries are reproduced inline in their section tables above.
Each row's notes column gives the legacy reader key. Resolution path for
each gap is one of:

1. **Add the field to an existing typed config** (e.g. `swfrost` ->
   extend `heat_config_t` with `swfrost: int`).
2. **Add a new typed config module** (e.g. `macropore_config_t` for the
   ~17 (G) fields under macropore -- Phase 4d intentionally deferred this).
3. **Reclassify as initial-condition state** if the legacy reader populates
   it once and physics never re-reads it from input -- extend the relevant
   `*_ini` slot rather than adding a runtime field.

Phase 4f decides per-gap which path to take. The audit's job is only to
surface them; this doc does not prescribe schema changes.

Highest-density gap clusters (by section):

- **Other / uncategorised** -- 102 gaps
- **Crop (fixed / grass / WOFOST)** -- 81 gaps
- **Time / control / output** -- 34 gaps
- **Soil + hydraulics** -- 29 gaps
- **Meteorology** -- 24 gaps
- **Drainage + surface water** -- 23 gaps
- **Macropore** -- 17 gaps
- **Bottom boundary** -- 7 gaps
- **Solute** -- 5 gaps
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
