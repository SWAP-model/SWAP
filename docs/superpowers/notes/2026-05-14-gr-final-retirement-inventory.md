# GR-FINAL Retirement Inventory

**Date:** 2026-05-14
**Source:** Audit of `src/core/variables.f90` against current consumer state.
**Method:** Per-symbol grep classification; consumers counted excluding `variables.f90`,
`initialize.f90`, `config_to_variables.f90`, and `swap_mod.f90` (the last is Phase C3
retirement territory — dual-write blocks are not real consumers).

---

## Summary

- **W** (Write-only, directly deletable): ~55 symbols
- **R** (Read-only-by-readers, need reader migration): ~65 symbols
- **B** (Both, need migration + adapter rewrite): ~110 symbols
- **C** (Config-loaded mirror, convert to direct sourcing): ~14 symbols

Total declarations audited: ~244 (as counted by grep -nE for type keywords)

---

## Per-category symbol lists

### W — Directly deletable

These have NO remaining consumers outside variables.f90/initialize.f90/config_to_variables.f90.
Phase D can delete declarations + adapter writes directly after Phase C clears dual-writes.

**Time & control:**
- `ex_tlast` — dead; handle_exchange retired; only initialize.f90 zeros it
- `flSwapShared` — 0 real consumers (swap_mod.f90 still has calls but those are Phase C3)
- `outdat(maout)`, `outdatint(maout)` — only c2v writes; swap_mod retired all use

**Irrigation (ADR 0009 Phase 5+ already zeroed output):**
- `swirg` — set 0 in c2v; IrrigationOutput deleted
- `irg` — file handle for *.IRG; never opened now
- `FlIrrigationOutput` — only initialize.f90 sets .false.
- `phormc` — only initialize.f90 zeros it

**Soilwater output handles (ADR 0009 Phase 5+):**
- `swafo`, `afo` — AFO output deleted
- `swaun`, `aun` — AUN output deleted
- `swbal`, `bal` — BAL output deleted
- `swblc`, `blc` — BLC output deleted
- `swsba`, `sba` — SBA output deleted
- `swwba`, `wba` — WBA output deleted
- `swvap`, `vap` — VAP output deleted
- `swstr`, `str` — STR output deleted (swstr=0 confirmed; str has 22 hits but only in file-open paths in swapoutput.f90 under `if (swstr==1)` dead branch)
- `bma` — BMA macropore balance never opened (flmacropore=.false.)
- `swini` — only initialize.f90 zeros it

**Crop flags with 0 consumers:**
- `swCrop` — only c2v writes it; all callers use state%crop instead
- `swjarvis` — deprecated switch; 0 total non-var hits
- `flHarvest`, `flHarvestpot` — 0 real consumers (0 hits outside)
- `flGrazing`, `flGrazingpot` — 0 real consumers

**Crop timing:**
- `DayGerm` — only initialize.f90 zeros it
- `hmow`, `hgrz` — 0 real consumers outside init

**Grass harvest tables (never read):**
- `DaysGrazingtab`, `UptGrazingtab`, `LossGrazingtab` — 0 real consumers

**Nitrogen:**
- `amFERT` — 0 real consumers

**Heat:**
- `nheat` — only written by initialize/c2v; heat_state.f90 declares its own nheat

**Snow:**
- `issnowbeg` — 0 real consumers

**Soilwater convergence state:**
- `Itnumb(100,2)` — 0 real consumers
- `nstep_hc` — 0 real consumers
- `swcapriseoutput` — 0 real consumers
- `qinfmax` — 0 real consumers (iqinfmax has 1 hit but isolated)

**Soil:**
- `SwSoilShr(MaHo)` — [retired-zero]; only initialize.f90 zeros it
- `ThetCrMp(MaHo)` — [retired-zero]; only initialize.f90 zeros it

**Macropore retired-zeros with 0 real consumers:**
- `dFdhMp(MaCp)` — zero-write only in soilhydraulics
- `iQInTopLatDm1`, `iQInTopLatDm2`, `iQInTopVrtDm1`, `iQInTopVrtDm2` — 0 consumers
- `iQExcMtxDm1Cp(MaCp)`, `iQExcMtxDm2Cp(MaCp)`, `iQOutDrRapCp(MaCp)` — 0 consumers
- `IcTopMP` — 0 consumers
- `IDecMpRat` — 0 consumers
- `RapDraReaExp`, `RapDraResRef(Madr)` — 0 real consumers (only c2v writes)

**Surface water:**
- `swswb`, `swdrf` — 0 real consumers
- `nrsec` — 0 real consumers
- `nqh(mamp)`, `drf`, `swb` — 0 real consumers
- `wls1_init` — only c2v writes; surfacewater_state_t owns it
- `wlstab(2*mawls)` — 0 real consumers

**Atmosphere SAVE (W after adapter writes die):**
- `atop(18,6)` — 0 real consumers (reprofunctions array never read)
- `ResultsOxygenStress(19,macp)` — 0 real consumers
- `o2_initialized` — 0 real consumers (o2_ini_stress has 3, but o2_initialized itself = 0)

---

### R — Need reader migration in Phase B

These have remaining consumers in source files that read the legacy global. Phase B migrates
them to read from `state%<subsystem>%<field>` instead.

**Meteo arrays (Phase B: meteoday/readmeteo reader cutover):**

| Symbol | Consumer count | Notes |
|--------|---------------|-------|
| `arad(366)` | 10 | used in meteoday, readmeteo, radiation calcs |
| `arai(366)` | 12 | rain array sliced daily |
| `atmn(366)` | 8 | min temp daily array |
| `atmx(366)` | 8 | max temp daily array |
| `awin(366)` | 7 | wind daily array |
| `aetr(366)` | 7 | ETref daily array |
| `ahum(366)` | 9 | humidity daily array |
| `wet(366)` | 21 | wet fraction array |
| `tav` | 29 | daily avg temperature scalar |
| `tmn` | 20 | daily min temperature |
| `tmx` | 19 | daily max temperature |
| `tmnr` | 4 | 7-day avg Tmin |
| `daylp` | 28 | photoperiodic daylength |
| `rad` | 0 excluded | scalar rad; has consumers in crop code |
| `atmin7(7)` | — | running Tmin window |
| `epot(96)`, `tpot(96)`, `grain(96)`, `nrain(96)`, `atav(96)` | various | detailed meteo arrays |
| `dethum/detrad/detrain/dettav/dettime/detwind(nmetfile)` | various | detailed arrays |
| `dtEventRain` | — | rain event timestep |
| `finterception` | — | interception ratio |
| `cfevappond` | — | pond evap ratio |
| `alt`, `altw`, `angstroma`, `angstromb` | low | station params |
| `lat` | — | latitude |
| `raintab(60)` | 2 | mean intensity table |
| `ad(mrain)`, `am(mrain)` | 13 each | meteo file day/month arrays |
| `daynrfirst`, `daynrlast` | 19/10 | meteo availability range |
| `detrecord(nmetfile)`, `irectotal`, `nmetdetail`, `nofd` | various | meteo control |

**Meteo control switches (Phase B or C depending on whether state%atmosphere gets them):**
- `swdivide` (16), `swetr` (9), `swetsine` (6), `swinter` (51), `swmetdetail` (37), `swrain` (21) — all read by meteo/ET code

**Meteo CSV cache (Phase B: meteoday slicing):**
- `nmetcsv`, `metcsv_dat`, `nmetcsv_det`, `metcsv_det`, `nraincsv`, `raincsv_dat`

**Meteo SAVE state:**
- `tsunrise_atm`, `tsunset_atm` (11/9 consumers) — ETSine time tracking
- `nod10_cn`, `icn_atm`, `z10_cn` (6/7/4) — CN method state

**Soilwater state (large category, Phase B reader cutover):**

These are genuine model parameters consumed broadly by soilwater / output code — they need
to move into a config-seeded struct or state%soilwater sub-record:

`paramvg(21,maho)` (35), `sptab/sptablay` (48/1), `iHWCKmodel(maho)` (60), `h_enpr(macp)` (35),
`bdens(maho)` (10), `ksatthr/relsatthr(maho)` (6/8), `BiModal/NoVap(maho)` (6/6),
`twilt(macp)` (11), `hcomp/hsublay(macp)` (3/4), `zi/zh(macp)` (1/4),
`cfbs` (13), `cofred` (11), `CritDevMasBal` (6), `CriterHr` (1), `gwlconv` (12), `gwli` (5),
`gwltab(mabbc*2)` (3), `haqtab/hbotab/qbotab` (2/2/6), `pondmx/pondmxtab` (12/2),
`drares/infres(Madr)` (6/6), `ftopdislay(Madr)` (5), `basegw` (4),
`shape` (14), `geofac` (8), `rimlay` (12), `entres` (8), `hdrain` (7), `rsigni` (9), `rsoil` (11),
`rsro/rsroexp` (9/6), `runonarr` (2), `khbot/khtop/kvbot/kvtop` (6/7/4/5),
`sinamp/sinave/sinmax` (6/6/6), `aqamp/aqave/aqper/aqtmax` (6 each),
`tau` (1), `StepHr/taccur` (3/14), `dznew(macp)` (11), `zintf` (9), `hplate` (10),
`cofqha/cofqhb/cofqhc` (6/6/6), `c_top(macp)` (6), `qimmob(macp)` (1), `qdraincomp(macp)` (2),
`qdrtab(50)` (2), `inpola/inpolb(macp)` (8/8), `CritDevh1Cp/CritDevh2Cp/CritDevPondDt` (2/2/1)

**Irrigation parameters (Phase B):**
- `cirrs`, `cirrthres`, `dcrit`, `ditab(14)`, `dwatab(14)`, `fidtab(14)`, `gird`, `hcritab(14)`,
  `irconc(mairg)`, `irdate(mairg)`, `irdepth(mairg)`, `nird`, `perirrsurp`, `raithreshold`,
  `rawtab/tawtab/tcritab/treltab/tstairrig/tendirrig` — all read by irrigation.f90

**SSDI persistent state (Phase B: migrate to state%soilwater or state%irrigation):**
- All `ssdi_*_irr`, `swssdi_irr`, `nod_ssdi_irr`, etc. — 3 consumers each in irrigation.f90

**Irrigation scheduling (Phase B):**
- `irrigevent` (29), `swirfix` (13), `nirri` (15), `swcirrthres` (6), `schedule` (77),
  `isua` (9), `isuas` (6), `dayfix` (7)

**Drainage / surface water (Phase B):**
- `dramet` (17), `ipos` (13), `ncomp/nsublay/numlay` (4/2/23), `nod1lay(maho)` (6),
  `numbit` (7), `numnodnew` (11), `numtab/numtablay` (20/3), `isoillay(maho)` (4),
  `botcom(maho)` (0), `nrstaring` (4), `nhead` (2)
- Surface water: `swsrf` (16), `swallo/swdtyp/swman` (7/13/20), `swqhr/swsec` (11/17),
  `nrpri/nmper` (8/45), `nphase/nodhd/intwl` (1/2/8), `nowltab(madr)` (4),
  `widthr/taludr/rdrain` (9/8/5), `rsurfdeep/rsurfshallow` (7/7), `rinfi/rentry/rexit` (5/6/6),
  `gwlinf/wlptab` (7/1), `impend/wldip/wscap/hbweir` (9/13/12/16), `osswlm/wlp` (7/3),
  `alphaw/betaw(mamp)` (12/12), `dropr/gwlcrit/hcrit/vcrit/hqhtab/qqhtab/wlsman` (various),
  `cofintfl/expintfl/owltab` (10/10/10)

**Output control (still read by swapoutput.f90):**
- `swrum` (3), `swinc` (3), `inc` (10), `swcsv` (3), `InList_csv` (5),
  `swcsv_tz` (3), `InList_csv_tz` (4), `tz_z1_z2(2)` (3),
  `sw2/sw3/sw4` (9/8/10), `swdra` (17), `swfrost` (13), `swhyst` (8),
  `swinco` (31), `swkmean` (41), `swkimpl` (5), `swliminf` (7),
  `swoxygen` (37), `swoxygentype` (13), `swpondmx` (2), `swqhbot/swcofqhc` (10/6),
  `swredu` (8), `swsophy` (24), `ientrytab/ientrytablay` (12/1), `swtopsub` (1),
  `swbotb3Impl/SwBotb3ResVert` (3/9), `swcfbs` (10), `swdiscrvert` (8),
  `swdislay/swtopdislay` (9/10), `swoutputmodflow` (0-excluded), `rot` (11),
  `dra` (1), `drfil/pathdrain` (2/3)

**Crop parameters (large — Phase B crop reader cutover):**

All active crop parameters (not already migrated to state%crop) with >0 consumers:
- Crop common: `croptype(macrop)` (52), `swcrp` (3), `crp` (17), `daycrop` (43), `icrop` (108),
  `swcf` (69), `swdrought` (49), `swgc` (15), `swcompensate` (36), `swstressor` (21),
  `swrootradius` (14), `swsalinity` (34), `swpotrelmf` (14), `cropfil/pathcrop/inifil` (3-4)
- Crop growth scalars: `dvs` (143), `lai` (86), `rd` (78), `rid` (69), etc. — all active
  crop dynamics fields with high consumer counts
- Crop config tables: `amaxtb`, `dtsmtb`, `fltb/fotb/frtb/fstb`, `slatb`, `rdrrtb/rdrstb`,
  `tmnftb/tmpftb`, `rfsetb`, `mrftb`, `rdctb/rdtb/rlwtb`, `gctb/cftb/cfeictb/chtb/wrtb`, etc.
- Grass harvest: `dmmowtb`, `dmgrztb`, `dateharvest`, `lsda`, `DelayRegrowthTab`,
  `lossmowtab/lossgrztab`, `cropstartpot/endpot/act/endact`, `pmowdm/mowdm/pgrzdm/grzdm/plossdm/lossdm`
- `flCropCalendar` (22), `flCropEmergence` (41), `flCropHarvest` (19),
  `flCropReadFile` (6), `flCropOpenFile` (2), `flCropOutput` (6)
- `flHarvestDay` (12), `flhrvendpot` (7), `flhrvendact` (6), `flanthesis` (1)
- Rooting: `rdi/rri/rdc/rdmax/rdm/rdpot` (all >14 consumers)
- CO2: `fco2amax/fco2eff/fco2tra` (4-5), `co2amaxtb/co2efftb/co2tratb` (4), `co2year/co2ppm` (1)
- Vernalisation: `verndvs/vernsat/vernbase/vernrtb`
- Bulb: `drbl/drblpot/dwbl/dwblpot/fbl/fbltb/pld/remoc/plwt/plwti/wbl/wblpot`
- Germination: `agerm/bgerm/cgerm/hdrygerm/hwetgerm/tsumgerm/tsumemeopt/TBASEM/TEFFMX`
- Nutrient parameters: `nmxlv/nlue/anlv/anst/nmaxlv/nmaxst/nmaxrt/lrnr/lsnr/nni/rnflv/rnfst/frnx`
  and `nlai/nmaxso/npart/nfixf/nsla/rnfrt/tcnt/dvsnlt/dvsnt/rdrns/fntrt/ilnmxl`
  and `fraharlosorm_lv/st/so/fstr`

**Tillage (all consumed only by tillage.f90, Phase B):**
- `till_swtill/i_n_model/iRedist/Ntill/Ntypes` (1-2 consumers each in tillage.f90)
- `till_Max_Z_tillage`, `till_Date/Z/I_tillage`, `till_Type/iType/iTT1/iTT2_Tillage`
- `till_TAB_Rho_tillage/cons/K_R_cons/Rho_match/N_match`

**Heat (Phase B: temperature.f90 reader cutover):**
- `swbotbhea` (8), `swtopbhea` (5), `swcalt` (10), `swhea` (9), `swtem` (3), `tem` (9)
- `ddamp/tampli/timref/tmean` (6 each), `tembtab/temtoptab(mabbc*2)` (8 each)
- `tfroststa/tfrostend` (14/21), `tsoil(macp)` (78-staging buffer note)
- `zh(macp)` (4)

**Snow (Phase B):**
- `snowcoef` (5), `TePrRain/TePrSnow` (4 each), `snw` (7), `swsnow` (12), `swsublim` (4)

**Solute (Phase B: solute reader cutover):**
- `nconc` (2), `swbr` (6), `swbotbc` (6), `swsolu` (12), `swsp` (0 real consumers — W?)
- `AgeGwl1m` (3), `bexp/decsat/gampar/rtheta` (6 each), `cdrain/dtsolu` (15/50-AgeTracer),
  `cirr/cpre/cref/ddif/daquif/poros/kfsat` (5-9), `cml/cmsy(macp)` (46/19),
  `cseeptab/kf/ldis/fdepth/decpot(maho)` (5-10), `frexp` (8), `salthead/saltmax/saltslope` (10/14/13),
  `samini/sqdra/rottot/isqbot/isqtop` (6-12-AgeTracer deps), `tscf` (8), `zc(macp)` (2)
- `icAgeBot/Dra/Rot/Sur` (2 each — AgeTracer accumulators)

**Age tracer (Phase B: agetracer.f90 reader cutover):**
- `flAgeTracer` (4), `Ageirr/Agedrain/Agepre/Agepond/Agepondm1` (2-9),
  `icAgetopupw/icAgetopdwn` (2 each)

**O2 stress SAVE state (Phase B: OxygenStress.f90 reader cutover):**
- `o2_w_root/w_root_z0/soil_temp/sat_water_cont/gas_filled_porosity/d_o2inwater/d_root/d_soil`
  `o2_perc_org_mat/soil_density/depth/shape_factor_microbialr/root_radius/r_microbial_z0`
  `o2_waterfilm_thickness/bunsencoeff/c_min_micro/c_macro/ctopnode` — all 1 consumer (OxygenStress.f90)
- `o2_d_soil_term1/term2/gfp100/capac_term/nmin1/mplus1` — 3 consumers each
- `o2_ini_stress` — 3 consumers

**Misc output:**
- `outfil` (34), `pathwork` (35), `project` (44), `logf` (57), `metfil` (3), `pathatm` (6)
- `es0/et0/ew0` (23/19/27) — potential ET scalars used by crop/output code

---

### B — Both migration + adapter rewrite

Most "R" symbols above are technically also "B" because config_to_variables.f90 writes them
AND runtime code reads them. The distinction matters for Phase ordering: Phase B migrates
readers first; Phase C retires the adapter writes; Phase D deletes the declaration.

Key B symbols with additional runtime mutators (written by non-adapter code during simulation):

- `daycrop`, `icrop`, `dvs`, `lai`, `rd`, `rdpot`, `rdm` — updated each day by cropgrowth.f90
- `tsum`, `tsumam`, `tsumea`, `laiexp/pot`, `glaiex/pot`, `wlv/pot`, `wrt/pot`, `wso/wsopot`, `wst/wstpot` — crop state updated daily
- `lv(366)/lvpot/lvage/lvagepot`, `sla/slapot(366)` — crop day-indexed arrays
- `daygrowth/pot`, `idregr/pot`, `iharvest`, `ilvold/pot`, `iseqgm/pot` — crop counters
- `idaysgraz/pot`, `seqgrazmow/pot(366)` — grazing counters
- `cuptgraz/pot`, `tagp/pot/t/tpot`, `tadw/pot`, `cwdm/pot`, `gasst/pot` — cumulative crop scalars
- `flCropCalendar/Emergence/Harvest/ReadFile/OpenFile/Output` — set by crop scheduler
- `flanthesis`, `flHarvestDay`, `flhrvendpot/act` — runtime flags toggled mid-simulation
- `flCropPrep/Sow/Germ` — set by crop pre-growth logic
- `PrepDelay/SowDelay`, `dhPrep/Sow`, `dtempSow` — runtime state for prep/sow delays
- `crp` (file handle, opened by cropinit) — runtime-mutated
- `rot` (file handle, opened by outrot) — runtime-mutated
- `inc` (file handle, opened by outinc) — runtime-mutated
- `tem` (file handle, opened by temperature output) — runtime-mutated
- `snw` (file handle, opened by snow output) — runtime-mutated
- `numbit` — updated each Richards iteration
- `dev_cmb` — file unit opened at runtime by checkmassbal
- `flwarn_hc/iwarn_hc` — warning state updated during headcalc
- `days_counter_irr`, `nirri_ssdi_irr` — SSDI state updated each day
- `iCNtab`, `CNdry/CNwet/ThetaRef/wc10/Runoff_CN` — CN state updated by meteoday
- `Agepond/Agepondm1`, `icAgetopupw/topdwn`, `icAgeBot/Dra/Rot/Sur` — updated by AgeTracer
- `icn_atm`, `tsunrise_atm/sunset_atm` — updated by meteoday/ETSine
- `daynrfirst/last` — updated each year by readmeteo
- All macropore retained-zero globals with in-loop writes: `cQMpLatSs`, `cQMpOutDrRap`,
  `iQMpOutDrRap`, `IAvFrMpWlWtDm1/2`, `IWaSrDm1/2Beg`, `WaSrDm1/2`, `WaSrDm1/2Ini`,
  `VlMpStDm1/2(MaCp)`, `QExcMpMtx`, `QMaPo`, `QRapDra`, `FlDecMpRat`, `flmacropore`

---

### C — Config-loaded mirrors

State mirrors that source directly from a config field. Phase C converts consumers to
read from `config%<path>` or `state%<subsystem>%<field>` directly; Phase D deletes legacy global.

| Symbol | Config path / notes |
|--------|---------------------|
| `swCrop` | config%crop%swcrop (already seeded via state%crop; global only in c2v/init) |
| `swend` | config%crop%common%swend → state%crop%common%swend; swap_mod reads it (Phase C3) |
| `swtsum/tsumtime/tsumtemp/tsumdepth` | config%crop%common (grass growth start params) |
| `swharv` | config%crop%common%swharv |
| `swWrtNonox/aeratecrit` | config%crop (oxygen stress root development) |
| `MaxPrepDelay/MaxSowDelay` | config%crop%common (prep/sow tolerance) |
| `swbulb` | config%crop%wofost (bulb switch, integer→logical conversion already in swap_mod) |
| `flCO2` | config%crop%wofost (CO2 correction flag) |
| `flCropNut` | config%crop (nutrient stress flag) |
| `flTillage` | config%tillage%swtill (call-site gate, ADR 0020) |
| `flSSDI` | config%irrigation%ssdi (call-site gate, ADR 0020) |
| `NumLevRapDra` | config%drainage%nrlevs (drain feature, not macropore) |
| `flAgeTracer` | config%solute%flagetracer |
| `flmacropore` | forced .false. in init; source of truth is compile-time constant |

---

## Phase ordering implications

**Phase B unlocks Phase D** in the following order:

1. **B-meteo** (migrate arad/arai/atmn/atmx/awin/aetr/ahum/wet/tav/tmn/tmx/tmnr and
   meteo switches) → unlocks retirement of ~25 meteo globals in Phase D.

2. **B-crop-common** (migrate croptype/swcf/swdrought/swgc/swcompensate/swstressor/swsalinity
   and the 100+ crop dynamics fields) — this is the largest Phase B task and the critical path.
   It unlocks retirement of the ~150 crop-family globals in Phase D.

3. **B-soilwater-config** (migrate paramvg/sptab/iHWCKmodel/h_enpr/bdens and soilwater
   control switches) → unlocks retirement of ~30 soilwater globals in Phase D.

4. **B-output** (migrate swrum/swinc/swcsv/swcsv_tz/sw2/sw3/sw4/swdra and remaining
   output switches) → unlocks Phase D deletion of ~15 output globals.

5. **B-solute/agetracer** (migrate cml/cmsy/cdrain/dtsolu and agetracer SAVE state)
   → unlocks solute globals in Phase D.

6. **B-O2** (migrate o2_* SAVE state to state%soilwater or new state%oxygenstress)
   → unlocks O2 globals in Phase D.

**Phase D safe-to-delete immediately** (no Phase B dependency):
All W-category symbols listed above — ~55 symbols deletable as soon as Phase C3 clears the
dual-write blocks from swap_mod.f90.

---

## Edge cases / TBD symbols

1. **`tsoil(macp)`**: Annotated as "config-staging buffer" — 78 consumers. Unusual: temperature.f90
   uses it as input at task=1 init, then state%heat%tsoil is the runtime value. Phase B must
   distinguish the init-read from runtime reads. Category: B (both config write and runtime read).

2. **`swsp`** (switch for solute sorption): 0 real consumers outside init/c2v. Should be W, but
   confirm it is not read by solute.f90 under `if (swsp==1)` branches. If those branches are dead
   due to existing dead-code removal, safe to classify W.

3. **`wc_cor`**: 4 consumers in CN method code. Appears to be a runtime state variable updated in
   meteoday.f90. Category B (written by meteoday, read by same). Needs migration.

4. **`dra`** (internal number for *.DRA input file): 1 real consumer (dra_open). Low-priority B.

5. **`iqinfmax` / `qinfmax`**: `iqinfmax` has 1 consumer; `qinfmax` has 0. Both are remnants of
   retired iqtdo/iqtup. Candidate W for `qinfmax`; `iqinfmax` may be W too if the 1 consumer is
   only an assign. Verify before Phase D.

6. **AgeTracer dead-code dependencies** (`cdrain`, `dtsolu`, `samini`, `sqdra`, `rottot`,
   `isqbot`, `isqtop`): These are annotated `[AgeTracer dead-code dep, keep until agetracer_state_t]`
   in variables.f90. Phase B agetracer_state_t migration will clean these up. Category: R (for now).

7. **`o2_initialized` vs `o2_ini_stress`**: Two distinct flags. `o2_initialized` has 0 real consumers
   (W); `o2_ini_stress` has 3 consumers in OxygenStress.f90 and is initialized via `data` statement.
   Category: o2_initialized → W; o2_ini_stress → B.

8. **Macropore `.BMA` writers** (`iQInTopLatDm1/2`, `iQInTopVrtDm1/2`, `IWaSrDm1/2Beg`,
   `WaSrDm1/2`, `WaSrDm1/2Ini`, `VlMpStDm1/2`, `IAvFrMpWlWtDm1/2`): The `.BMA` writer calls
   are in waterbalance.f90 under `if (flmacropore)` which is always .false. The write-side is dead
   but the read of the global exists. All are W (0 real consumers confirmed above).

9. **`swsba`**: Shows 0 real consumers outside vars/init/c2v (ADR 0009 deleted outsba). Confirmed W.

10. **`swend`** (category C above): swap_mod.f90 still reads it directly for `if (swend.eq.1)` and
    `if (swend.eq.2)` — it seeds `state%crop%common%swend` but also reads the raw global. Phase C3
    must move these two conditionals to read `state%crop%common%swend` before Phase D can delete it.
