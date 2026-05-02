# readwofost Line-by-Line Audit

**File:** `src/io/readswap.f90`, lines 2520–3244  
**Subroutine:** `readwofost(icrop, crpfil, swhydrlift, swsoybean, mg, dvsi, dvrmax1, dvrmax2, flrfphotoveg, tmaxdvr, tmindvr, toptdvr, popt, pcrt, flphenodayl, fradeceasedlvtosoil)`  
**Date of audit:** 2026-05-02  
**Purpose:** Classification of every non-blank, non-comment line in the ~725-line `readwofost` subroutine to inform the `.crp` → TOML port (Phase 2, Tasks 2–5). The subroutine reads the legacy `*.crp` ASCII file, validates its contents, applies hardcoded defaults, and builds runtime data structures (notably the cumulative root density function and the `.END`-file restart state). After the port, it is replaced by:

- Config validators (VALIDATE-bucket lines)
- Config finalizers / defaults (NORMALIZE-bucket lines)
- A new `cropwofost_init` module (RUNTIME-bucket lines)
- Stub/error checks for unimplemented branches (GUARDED-bucket lines)

Lines that fall into multiple buckets appear in multiple rows — this is intentional and useful.

---

## Classification Table

| Lines         | Bucket    | Notes                                                                                                                                          |
| ------------- | --------- | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| 2520-2521     | READ      | Subroutine signature — `readwofost(icrop, crpfil, swhydrlift, swsoybean, …, fradeceasedlvtosoil)`                                             |
| 2528-2543     | READ      | `use variables, only: …` — imports all module-level wofost crop state variables                                                               |
| 2544-2547     | READ      | `use irrigation_mod, only: irrigation`; `use oxygenstress_mod, only: oxygen_dat`; `use array_utils, only: afgen`; `use swap_array_dimensions` |
| 2548          | READ      | `implicit none`                                                                                                                                |
| 2550-2578     | READ      | Dummy argument and local variable declarations (`icrop`, `swhydrlift`, `swsoybean`, soybean params, CO2 locals, temporaries, `rdinqr`)         |
| 2583          | READ      | `filnam = trim(pathcrop)//trim(crpfil)//'.crp'` — build file path                                                                             |
| 2584          | READ      | `crp = getun2(10,90,2)` — get free file unit number                                                                                           |
| 2585          | READ      | `call rdinit(crp,logf,filnam)` — open `.crp` file and initialise parser                                                                       |
| 2588          | READ      | `call rdsinr('swcf',1,3,swcf)` — read crop factor / crop height switch                                                                        |
| 2591-2595     | VALIDATE  | `if (swetr.eq.1 .and. swcf.eq.2)` fatalerr — ETref requires swcf=1 or swcf=3, not swcf=2                                                     |
| 2597-2606     | READ      | `if (swcf.eq.1)`: read `dvs` + `cf` arrays into `cftb`; pack interleaved                                                                     |
| 2602-2605     | RUNTIME   | Loop: pack `dvsinput`/`cfinput` into interleaved `cftb` array                                                                                 |
| 2606          | NORMALIZE | `chtb = -99.99d0` — sentinel fill for unused crop-height table when swcf=1                                                                    |
| 2607-2616     | READ      | `elseif (swcf.eq.2)`: read `dvs` + `ch` arrays into `chtb`; pack interleaved                                                                 |
| 2612-2615     | RUNTIME   | Loop: pack `dvsinput`/`chinput` into interleaved `chtb` array                                                                                 |
| 2616          | NORMALIZE | `cftb = -99.99d0` — sentinel fill for unused crop-factor table when swcf=2                                                                    |
| 2617-2637     | GUARDED   | `elseif (swcf.eq.3)`: read `lai`, `cf`, `cfeic` + `ch` arrays; pack into `cftb`, `cfeictb`, `chtb`. LAI-dependent dual crop coefficient; not in port scope. |
| 2640          | READ      | `call rdsinr('swinter',0,3,swinter)` — interception model switch                                                                              |
| 2641-2642     | READ      | `if (swinter.eq.1)`: `call rdsdor('cofab',…)` — Von Hoyningen-Huene coefficient                                                              |
| 2643-2661     | GUARDED   | `else if (swinter.eq.2)`: read `t`, `pfree`, `pstem`, `scanopy`, `avprec`, `avevap`; pack into 5 interleaved lookup tables (Gash model). Not in port scope. |
| 2662-2665     | GUARDED   | `else if (swinter.eq.3)`: read `fimin`, `siccaplai` — storage-capacity interception. Not in port scope.                                       |
| 2668-2672     | NORMALIZE | `if (swcf.eq.1 .or. swcf.eq.3)`: `albedo = 0.23d0`, `rsc = 70.0d0`, `rsw = 0.0d0` — hardcoded ETref standard defaults                       |
| 2675-2677     | READ      | `else` (swcf=2): `call rdsdor('albedo',…)`, `call rdsdor('rsc',…)`, `call rdsdor('rsw',…)` — crop-specific values                            |
| 2682-2685     | NORMALIZE | `swsoybean = 0`; `if (rdinqr('swsoybean'))`: `call rdsinr('swsoybean',0,1,swsoybean)` — default=0; soybean switch read                        |
| 2686-2706     | GUARDED   | `if (swsoybean.eq.1)`: read `mg`, `dvsi`, `dvrmax1`, `dvrmax2`, `tmaxdvr`, `tmindvr`, `toptdvr`; optional `flrfphotoveg`, `flphenodayl`, `popt`, `pcrt`. Entire soybean phenology block guarded. |
| 2712-2721     | READ      | `if (swsoybean.eq.0)`: `call rdsinr('idsl',0,2,idsl)`; if idsl=1/2 read `dlo`, `dlc`; read `tsumea`, `tsumam`, `dtsmtb`                      |
| 2714-2717     | GUARDED   | `if (idsl.eq.1 .or. idsl.eq.2)`: read `dlo`, `dlc` — daylength dependency. Case 5 has idsl=0, so dlo/dlc are guarded.                        |
| 2724-2729     | GUARDED   | `if (idsl.eq.2)`: read `verndvs`, `vernsat`, `vernbase`, `vernrtb` — vernalisation block. Guarded (case 5 has idsl=0).                        |
| 2732-2734     | NORMALIZE | `swbulb = 0`; `if (rdinqr('swbulb'))`: `call rdsinr('swbulb',0,1,swbulb)` — default=0; bulb switch read                                      |
| 2736-2742     | GUARDED   | `if (swbulb.eq.1)`: read `fbltb`, `plwti`, `pld`, `remoc` — bulb crop parameters. Not in port scope.                                         |
| 2745          | READ      | `call rdsdor('dvsend',0,3,dvsend)` — development stage at harvest (always read, no default)                                                   |
| 2746-2748     | NORMALIZE | `swharv = 0`; `if (rdinqr('swharv'))`: `call rdsinr('swharv',0,1,swharv)` — default=0                                                        |
| 2752-2754     | READ      | `call rdsdor('tdwi',…)`, `call rdsdor('laiem',…)`, `call rdsdor('rgrlai',…)` — initial dry weight, LAI, relative LAI increase                 |
| 2757-2761     | READ      | `call rdador('slatb',…)`, `call rdsdor('spa',…)`, `call rdsdor('ssa',…)`, `call rdsdor('span',…)`, `call rdsdor('tbase',…)` — green area      |
| 2764-2769     | READ      | `call rdsdor('kdif',…)`, `call rdsdor('kdir',…)`, `call rdsdor('eff',…)`, `call rdador('amaxtb',…)`, `call rdador('tmpftb',…)`, `call rdador('tmnftb',…)` — assimilation |
| 2772-2775     | READ      | `call rdsdor('cvl',…)`, `call rdsdor('cvo',…)`, `call rdsdor('cvr',…)`, `call rdsdor('cvs',…)` — assimilate conversion efficiencies           |
| 2778-2783     | READ      | `call rdsdor('q10',…)`, `call rdsdor('rml',…)`, `call rdsdor('rmo',…)`, `call rdsdor('rmr',…)`, `call rdsdor('rms',…)`, `call rdador('rfsetb',…)` — maintenance respiration |
| 2786-2789     | READ      | `call rdador('frtb',…)`, `call rdador('fltb',…)`, `call rdador('fstb',…)`, `call rdador('fotb',…)` — partitioning tables                      |
| 2792-2794     | READ      | `call rdsdor('perdl',…)`, `call rdador('rdrrtb',…)`, `call rdador('rdrstb',…)` — death rates                                                  |
| 2797-2799     | NORMALIZE | `swoxygen = 1`; `if (rdinqr('swoxygen'))`: `call rdsinr('swoxygen',0,2,swoxygen)` — default=1                                                 |
| 2801-2806     | VALIDATE  | `if (swoxygen.eq.2 .and. (swhea.eq.0 .or. swcalt.eq.1))` fatalerr — physical O2 stress needs numerical heat flow                              |
| 2808-2813     | VALIDATE  | `if (swoxygen.eq.2 .and. bdens(1).lt.100)` fatalerr — physical O2 stress needs realistic bulk density                                         |
| 2816-2821     | READ      | `if (swoxygen.eq.1)`: read `hlim1`, `hlim2u`, `hlim2l` — Feddes anaerobiosis pressure heads                                                   |
| 2823-2856     | GUARDED   | `if (swoxygen.eq.2)` block — Bartholomeus physical O2 stress. Default `swoxygentype=1`; optional `swoxygentype` read; swoxygentype=1: read q10_microbial, specific_resp_humus, srl, swrootradius + conditionals; swoxygentype=2: read SwTopSub, NrStaring; `call oxygen_dat(…)`. Entire block out of port scope. |
| 2859          | NORMALIZE | `swWrtNonox = 0` — default: no non-oxidative root senescence switch                                                                           |
| 2860-2862     | READ      | `if (rdinqr('swWrtNonox'))`: `call rdsinr('swWrtNonox',0,1,swWrtNonox)`                                                                       |
| 2864          | NORMALIZE | `aeratecrit = 0.0001d0` — default aeration criterion                                                                                          |
| 2865-2867     | READ      | `if (swWrtNonox.eq.1)`: `call rdsdor('aeratecrit',…)` — read aeration criterion threshold                                                     |
| 2870          | NORMALIZE | `swdrought = 1` — default: drought stress according to Feddes                                                                                  |
| 2871-2873     | READ      | `if (rdinqr('swdrought'))`: `call rdsinr('swdrought',1,2,swdrought)`                                                                          |
| 2875-2881     | READ      | `if (swdrought.eq.1)`: read `hlim3h`, `hlim3l`, `hlim4`, `adcrh`, `adcrl` — Feddes drought parameters                                        |
| 2883-2897     | GUARDED   | `else` (swdrought=2): read De Jong van Lier params: `wiltpoint`, `kstem`, `rxylem`, `rootradius`, `kroot`, `rootcoefa`, `swhydrlift`, `rooteff`, `stephr`, `criterhr`, `taccur`. Not in port scope. |
| 2900-2924     | READ      | `if (flsolute)` salt stress block: optional read `swsalinity`; if swsalinity=1 read `saltmax`, `saltslope`                                    |
| 2912-2923     | GUARDED   | `elseif (swsalinity.eq.2)`: read `salthead`; validate swdrought must be 2 — osmotic head concept. Guarded.                                    |
| 2916-2923     | VALIDATE  | `if (swdrought.eq.1)` inside swsalinity=2 block: fatalerr — osmotic head salt stress requires De Jong van Lier                                |
| 2927-2928     | NORMALIZE | `swcompensate = 0; swjarvis = 0` — defaults for compensation switches                                                                         |
| 2929-2930     | READ      | `if (rdinqr('swcompensate'))`: `call rdsinr('swcompensate',0,2,swcompensate)`                                                                 |
| 2932-2944     | READ      | `else if (rdinqr('swjarvis'))` legacy path: warn, read `swjarvis`, derive `swcompensate=1`                                                    |
| 2933-2934     | VALIDATE  | `call warn(…)` — SWJARVIS is deprecated; migrate to SWCOMPENSATE                                                                              |
| 2937-2940     | VALIDATE  | `call warn(…)` — SWJARVIS applied to all stresses (when swjarvis≠4)                                                                           |
| 2941          | NORMALIZE | `swcompensate = 1` — coerce legacy swjarvis > 0 to new swcompensate=1                                                                         |
| 2947-2952     | VALIDATE  | `if (swcompensate.gt.0 .and. swdrought.eq.2)` fatalerr — compensation not allowed with De Jong van Lier drought                               |
| 2955-2960     | READ      | `swstressor = 1`; `if (swcompensate.gt.0 .and. rdinqr('swstressor'))`: `call rdsinr('swstressor',1,5,swstressor)` — default=1                 |
| 2955          | NORMALIZE | `swstressor = 1` — default: compensate all stressors                                                                                          |
| 2963-2968     | GUARDED   | `if (swcompensate.eq.1)`: `alphacrit = 1.0d0`; `call rdsdor('alphacrit',…)` — Jarvis compensation. Guarded (case 5 has swcompensate=0).       |
| 2963          | NORMALIZE | `alphacrit = 1.0d0` — default critical stress index (inside swcompensate=1 guard)                                                             |
| 2970-2976     | GUARDED   | `elseif (swcompensate.eq.2)`: `dcritrtz = 0.0d0`; `call rdsdor('dcritrtz',…)` — Walsum compensation. Not in port scope.                      |
| 2979-2982     | NORMALIZE | `relmf = 1.0d0`; `if (rdinqr('relmf'))`: `call rdsdor('relmf',…)` — management factor default=1.0                                            |
| 2985-2988     | NORMALIZE | `swpotrelmf = 1`; `if (rdinqr('swpotrelmf'))`: `call rdsinr('swpotrelmf',1,2,swpotrelmf)` — default=1                                         |
| 2993-2996     | NORMALIZE | `swrdc = 0`; `if (rdinqr('swrdc'))`: `call rdsinr('swrdc',0,1,swrdc)` — default=0                                                             |
| 2999          | READ      | `call rdador('rdctb',0,100,rdctb,22,ifnd)` — root density distribution table (always read)                                                    |
| 3002-3005     | NORMALIZE | `swrd = 2`; `if (rdinqr('swrd'))`: `call rdsinr('swrd',1,3,swrd)` — default=2 (note: differs from cropfixed default=1)                        |
| 3008-3010     | GUARDED   | `if (swrd.eq.1)`: `call rdador('rdtb',…)` — rooting depth vs DVS table. Guarded in case 5 (swrd=2).                                          |
| 3013-3023     | READ      | `elseif (swrd.eq.2)`: read `rdi`, `rri`, `rdc`; optional `swdmi2rd` — max-daily-increase root extension. Active in case 5.                   |
| 3020-3023     | NORMALIZE | `swdmi2rd = 0`; `if (rdinqr('swdmi2rd'))`: `call rdsinr('swdmi2rd',0,1,swdmi2rd)` — default=0                                                 |
| 3026-3031     | GUARDED   | `elseif (swrd.eq.3)`: read `rlwtb`, `wrtmax` — root biomass-based extension. Not in port scope.                                               |
| 3034          | READ      | `call rdsinr('schedule',0,1,schedule)` — irrigation scheduling switch (always read)                                                           |
| 3036-3041     | GUARDED   | `if (schedule.eq.1 .and. swdrought.eq.2)`: re-read `hlim3h`, `hlim3l`, `hlim4` for scheduling override. Guarded (case 5 has schedule=0).     |
| 3044-3047     | NORMALIZE | `FraDeceasedLvToSoil = 0.0d0`; `if (rdinqr('FraDeceasedLvToSoil'))`: `call rdsdor(…)` — default=0.0                                          |
| 3050-3054     | NORMALIZE | `flco2 = .false.; swco2 = 0`; `if (rdinqr('swco2'))`: `call rdsinr('swco2',0,1,swco2)`; `if (swco2.eq.1) flco2 = .true.` — default=0         |
| 3058-3083     | GUARDED   | `if (flco2)` block: read `CO2AMAXTB`, `CO2EFFTB`, `CO2TRATB`; optionally read `atmofil`; build filnam for `.co2` file; `close(crp)`; open `.co2` file and read `CO2year`, `CO2ppm`; `close(uco2)`. Entire CO2 block guarded (case 5 has swco2=0). |
| 3081          | READ      | `else`: `close(crp)` — close `.crp` file in the non-CO2 path                                                                                 |
| 3086-3089     | GUARDED   | `if (schedule.eq.1)`: `call irrigation(1)` + optional `IrrigationOutput(1)` — scheduling initialisation. Guarded.                            |
| 3094-3121     | RUNTIME   | `if (swdrought.eq.1)` block: build `cumdens` array — normalized cumulative root density function over 101 depth points (3 loops; depends only on `rdctb`) |
| 3097-3101     | RUNTIME   | Inner loop i=0,100: populate `rootdis` using `afgen(rdctb,22,depth)` at 0.01-spaced depths                                                   |
| 3104-3115     | RUNTIME   | Cumulative density computation: copy depths to odd indices; trapezoidal integration loop i=4,202,2; accumulate `sum`, store in `cumdens` even indices |
| 3118-3120     | RUNTIME   | Normalization loop i=2,202,2: `cumdens(i) = cumdens(i) / sum`                                                                                 |
| 3124-3241     | GUARDED   | `.END`-file restart block: `if (t1900-tstart.lt.1d-3 .and. swinco.eq.3 .and. .not. abs(t1900-cropstart(icrop)).lt.1d-3)` — entire block guarded (swinco=3 only). |
| 3127-3130     | GUARDED   | `if (flCropCalendar)`: `ini = getun2(…)`; `call rdinit(ini,logf,inifil)` — open `.END` file                                                   |
| 3132-3135     | GUARDED   | `if (rdinqr('swcropharvest'))`: read `swcropharvest`; derive `flCropHarvest` (swinco=3 only)                                                  |
| 3141-3146     | GUARDED   | Read `swCropEmergence`; derive `flCropEmergence` (swinco=3 only)                                                                              |
| 3148-3150     | GUARDED   | Read `PrepDelay`, `SowDelay`, `tsumgerm` (swinco=3 only)                                                                                      |
| 3155-3163     | GUARDED   | `if (.not. flCropHarvest)`: read `swIrrigate`; derive `flIrrigate`; conditionally read `dayfix` (swinco=3 only)                               |
| 3166-3218     | GUARDED   | Read full crop state variables from `.END` file: `rd`, `rdpot`, optional `sicact` (swinter=3 only), `dvs`, `daycrop`, `swanthesis`, `tsum`, `ilvold`, `ilvoldpot`, `wrt`, `wrtpot`, `tadw`, `tadwpot`, `wst`, `wstpot`, `wso`, `wsopot`, `wlv`, `wlvpot`, `laiexp`, `laiexppot`, `glaiex`, `glaiexpot`, `lai`, `laipot`, `laimax`, `dwrt`, `dwrtpot`, `dwlv`, `dwlvpot`, `dwst`, `dwstpot`, `dwlvSoil`, `dwlvCrop`, `gasst`, `gasstpot`, `mrest`, `mrestpot`, `cwdm`, `cwdmpot`, `nofd`, `atmin7`, arrays `sla`, `slapot`, `lvage`, `lvagepot`, `lv`, `lvpot`. Derive `flanthesis`. (swinco=3 only) |
| 3226-3233     | GUARDED   | `else` (flCropHarvest=true or rdinqr('swcropharvest')=false): set `daycrop=0`, `flCropEmergence=.false.` as defaults (swinco=3 path)          |
| 3239          | GUARDED   | `close(ini)` — close `.END` file (swinco=3 only)                                                                                              |
| 3243          | READ      | `return`                                                                                                                                       |
| 3244          | READ      | `end` — end of subroutine                                                                                                                      |

---

## Summary

### Lines per bucket

| Bucket    | Approx. line count | Description                                                                                                                    |
| --------- | ------------------ | ------------------------------------------------------------------------------------------------------------------------------ |
| READ      | ~180 lines         | `rd*` calls, file open/close, declarations, use statements; mandatory reads for all wofost biophysics parameters               |
| VALIDATE  | ~25 lines          | `fatalerr` / `warn` calls and surrounding `if` logic (swcf+swetr, swoxygen=2+heat/bdens, swcompensate+swdrought, swsalinity=2) |
| NORMALIZE | ~35 lines          | Defaults: swsoybean=0, swbulb=0, swharv=0, swoxygen=1, swWrtNonox=0, aeratecrit=1e-4, swdrought=1, swcompensate=0, swstressor=1, alphacrit=1.0, swrdc=0, swrd=2, FraDeceasedLvToSoil=0.0, relmf=1.0, swpotrelmf=1, swdmi2rd=0, flco2=false; plus chtb/cftb sentinel fills, albedo/rsc/rsw ETref defaults |
| RUNTIME   | ~30 lines          | cumdens build after close (3 loops); cftb/chtb packing loops (swcf=1 and swcf=2 paths)                                        |
| GUARDED   | ~200 lines         | swcf=3, swinter=2/3, swsoybean=1, idsl=1/2 (dlo/dlc), idsl=2 (vernalisation), swbulb=1, swoxygen=2, swdrought=2, swsalinity=2, swcompensate=1/2, swrd=1, swrd=3, schedule=1, swco2=1, swinco=3 warm-restart |

Total: ~725 lines of substantive content (2520–3244).

### Call count summary

- **READ:** ~60 `rd*` calls (rdsinr, rdsdor, rdador, rdfdor, rdfinr, rdinqr, rdinit, rdslog, rdscha) plus file operations (getun2×2, rdinit×2, close×2)
- **VALIDATE:** 5 `fatalerr` calls + 2 `warn` calls
- **NORMALIZE:** ~17 distinct default assignments: swsoybean, swbulb, swharv, swoxygen, swWrtNonox, aeratecrit, swdrought, swcompensate/swjarvis, swstressor, alphacrit, swrdc, swrd, swdmi2rd, FraDeceasedLvToSoil, relmf, swpotrelmf, flco2; plus chtb/cftb sentinel fills, albedo/rsc/rsw ETref defaults
- **RUNTIME:** 1 initialization action: cumdens build (3 loops); 2 table-packing loops (cftb, chtb)
- **GUARDED:** 14 conditional branches: swcf=3, swinter=2, swinter=3, swsoybean=1, idsl=1/2 (dlo/dlc), idsl=2 (vernalisation), swbulb=1, swoxygen=2, swdrought=2, swsalinity=2, swcompensate=1, swcompensate=2, swrd=1, swrd=3, schedule=1, swco2=1, swinco=3

---

## Scope for Phase 2

### Active (case 5 exercises)

| Switch         | Case 5 value | Read location       |
| -------------- | ------------ | ------------------- |
| swcf           | 2            | 2588                |
| swinter        | 1            | 2640                |
| swsoybean      | 0 (default)  | 2682-2685           |
| idsl           | 0            | 2713                |
| swbulb         | 0 (default)  | 2732-2734           |
| swoxygen       | 1 (default)  | 2797-2799           |
| swWrtNonox     | 1            | 2860-2862           |
| swdrought      | 1 (default)  | 2870-2873           |
| swsalinity     | 1            | 2902-2910           |
| swcompensate   | 0 (default)  | 2927-2930           |
| swrd           | 2 (default)  | 3002-3005           |
| swdmi2rd       | 1            | 3020-3023           |
| schedule       | 0            | 3034                |
| swco2          | 0 (default)  | 3050-3054           |
| swpotrelmf     | 2            | 2985-2988           |
| FraDeceasedLvToSoil | 0.3     | 3044-3047           |

### Guarded (stub-errored in port)

swcf=3, swinter=2/3, swsoybean=1, idsl=1/2 (dlo/dlc), idsl=2 (vernalisation), swbulb=1, swoxygen=2, swdrought=2, swsalinity=2, swcompensate=1/2, swrd=1, swrd=3, schedule=1, swco2=1, flCropNut=.true. (N-P-K; see Sibling Readers), swinco=3 (warm-restart condition may evaluate true for some rotations — see spec risk register).

---

## Sibling Readers

This section documents every `.crp` open reachable from the wofost runtime path, per ADR 0017 (sibling-reader dispatch pattern first discovered during Phase 1 Task 8 for `readarablelandgerm`).

### ADR 0017 grep results

Command run:
```
grep -n "trim(pathcrop)\|'.crp'" src/io/readswap.f90 src/crop/*.f90 \
  | grep -v "outfil\|outputfile" | sort -u
```

Relevant hits on wofost call path:

| File                          | Line | Open                                    | Guard                        |
| ----------------------------- | ---- | --------------------------------------- | ---------------------------- |
| `src/io/readswap.f90`         | 2583 | `readwofost` primary `.crp` open        | Always                       |
| `src/io/readswap.f90`         | 3273 | `readarablelandgerm` `.crp` open        | Called from wofost path (see below) |
| `src/crop/cropgrowth.f90`     | 1063 | N-P-K nutrient block `.crp` open        | `if (flCropNut)` — guarded   |
| `src/crop/irrigation.f90`     | 77   | `irrigation(1)` `.crp` open             | `if (schedule.eq.1)` — guarded |

### Sibling 1: `readarablelandgerm` (called via `ArableLandGerm`)

**Location:** `src/crop/cropgrowth.f90:134` (legacy fallback), `src/io/readswap.f90:3248-3434`

`ArableLandGerm` is called from `cropgrowth.f90` on every new wofost crop emergence check (lines 91-149), regardless of crop type (fixed or wofost). It opens the same `.crp` file as `readwofost` to read `SwPrep`, `SwSow`, `SwGerm` and their dependent parameters.

The ADR 0017 cache-hit dispatch currently fires for `rotation_type==1` (cropfixed) only. For `rotation_type==2` (wofost), it falls back to the legacy `readarablelandgerm` call (line 134). Phase 2 must extend the dispatch to cover the wofost cache.

Case 5 uses `SwPrep=0`, `SwSow=0`, `SwGerm=0` — all defaults. All three active paths are the zero-case (flCropPrep/flCropSow/flCropGerm all set to `.true.` without reading any dependent fields). This means the Phase 2 stub-error for swprep/swsow/swgerm≠0 (already coded for cropfixed at cropgrowth.f90:122-124) must be replicated for the wofost cache path.

### Sibling 2: N-P-K nutrient block in `cropgrowth.f90`

**Location:** `src/crop/cropgrowth.f90:1061-1091`

After `readwofost` returns, `cropgrowth.f90` opens the same `.crp` file a second time when `flCropNut` is true (line 1061). It reads: `LRNR`, `LSNR`, `NLAI`, `NLUE`, `NMAXSO`, `NPART`, `NFIXF`, `NSLA`, `RNFLV`, `RNFRT`, `RNFST`, `TCNT`, `DVSNLT`, `DVSNT`, `RDRNS`, `FNTRT`, `FRNX`, `NMXLV` (table), and `FraHarLosOrm_lv`, `FraHarLosOrm_st`, `FraHarLosOrm_so`.

This block is **not** inside `readwofost` and is **not** in `src/io/readswap.f90`. The plan's reference to "lines 994-1021" is incorrect — those lines do not correspond to readwofost; the actual nutrient block is at `src/crop/cropgrowth.f90:1061-1091`.

The `flCropNut` flag is set from the `.swp` file (readswap.f90:474-476), not from the `.crp` file. Case 5 does not set `flCropNut = .true.`, so this block is guarded and out of Phase 2 scope.

### Sibling 3: `irrigation(1)` in `src/crop/irrigation.f90`

**Location:** `src/crop/irrigation.f90:77-79`

When `schedule=1`, `readwofost` calls `irrigation(1)` at line 3086-3088. `irrigation(task=1)` opens the same `.crp` file to re-read `schedule`, `startirr`, `endirr`, and all irrigation scheduling parameters. This is entirely guarded by `schedule=1` and is out of Phase 2 scope (case 5 has schedule=0).

---

## Implementation Guidance

### 1. swrd default differs from cropfixed

`readwofost` sets `swrd = 2` as default (line 3002), whereas `readcropfixed` sets `swrd = 1`. The TOML schema for cropwofost must document `swrd` default as 2. Case 5 exercises the swrd=2 path (rdi/rri/rdc/swdmi2rd).

### 2. dvsend is always read (no rdinqr guard)

Unlike `readcropfixed` where dvsend is optional with a default of 2.0, `readwofost` reads `dvsend` unconditionally at line 2745 with no default. The TOML schema must require dvsend for cropwofost (not optional).

### 3. FraHarLosOrm is NOT in readwofost — it is in a sibling block

The three `FraHarLosOrm_*` fields (lv, st, so) appear in `cropgrowth.f90:1086-1088`, inside the `flCropNut` guard. They are not read by `readwofost` and should not be included in the cropwofost TOML schema unless/until the N-P-K path is ported.

### 4. swco2 opens a second file (`.co2`) inside readwofost

When `swco2=1`, the `.crp` file is closed early (line 3071) and a `.co2` file (defaulting to `Atmospheric.co2`) is opened to read `CO2year` and `CO2ppm`. The CO2 block is entirely guarded in Phase 2. The stub should refuse `swco2=1` at config validation time.

### 5. The .END-file restart block is significantly more complex than in cropfixed

The wofost `.END`-file block (lines 3124-3241) reads ~40 state variables (full crop state arrays including `sla`, `lvage`, `lv`, etc.), compared to ~7 variables in `readcropfixed`. It is structured around `flCropCalendar` and `rdinqr('swcropharvest')` guards. The swinco=3 path should be refused by the config validator for Phase 2.

### 6. Salt stress (swsalinity=1) is active in case 5 but gated on flsolute

The entire salt stress section (lines 2900-2924) is inside `if (flsolute)`. Case 5 has flsolute=.true. and swsalinity=1. The TOML crop validator cannot check this alone — the cross-file validator must gate swsalinity checks on the companion `.swp` flsolute setting.

### 7. swWrtNonox=1 and aeratecrit are both active in case 5

Unlike cropfixed (where swWrtNonox=0 is the default and aeratecrit is rarely set), case 5 has swWrtNonox=1 with aeratecrit=0.5. Both must be in the active (non-guarded) port path.

### 8. cumdens is computed post-close, swdrought=1 only

Identical pattern to cropfixed (lines 3094-3121). Belongs in `cropwofost_init`. Only computed when swdrought=1.

### 9. swjarvis legacy shim present

Lines 2932-2944 replicate the same swjarvis deprecation shim as in readcropfixed. The TOML schema should only expose `swcompensate`; swjarvis must not be replicated.
