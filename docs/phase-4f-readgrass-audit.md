# readgrass Line-by-Line Audit

**File:** `src/io/readswap.f90`, lines 3437–4270  
**Subroutine:** `readgrass(icrop, crpfil, swharvest, dmharvest, daylastharvest, dmlastharvest, swdmmow, maxdaymow, swlossmow, swlossgrz, swdmgrz, maxdaygrz, dmgrazing, lsdb, tagprest, swhydrlift)`  
**Date of audit:** 2026-05-02  
**Purpose:** Classification of every non-blank, non-comment line in the 834-line `readgrass` subroutine to inform the `.crp` → TOML port (Phase 3, Tasks 2–5). The subroutine reads the legacy `*.crp` ASCII file, validates its contents, applies hardcoded defaults, and builds runtime data structures (notably the cumulative root density function and the `.END`-file restart state). After the port, it is replaced by:

- Config validators (VALIDATE-bucket lines)
- Config finalizers / defaults (NORMALIZE-bucket lines)
- A new `cropgrass_init` module (RUNTIME-bucket lines)
- Stub/error checks for unimplemented branches (GUARDED-bucket lines)

Lines that fall into multiple buckets appear in multiple rows — this is intentional and useful.

**Line range confirmed:** The plan stated 3437–4270. Grep confirms `subroutine readgrass` at 3437 and the first `end` after it at 4270. The range matches exactly — no drift.

---

## Classification Table

| Lines         | Bucket    | Notes                                                                                                                                                           |
| ------------- | --------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 3437-3439     | READ      | Subroutine signature — `readgrass(icrop, crpfil, swharvest, dmharvest, …, swhydrlift)`                                                                        |
| 3441-3445     | READ      | Header comment block                                                                                                                                            |
| 3446-3463     | READ      | `use variables, only: …` — imports all module-level grass crop state variables                                                                                 |
| 3464-3466     | READ      | `use irrigation_mod, only: irrigation`; `use oxygenstress_mod, only: oxygen_dat`; `use array_utils, only: afgen`; `use swap_array_dimensions`                 |
| 3467          | READ      | `implicit none`                                                                                                                                                 |
| 3470-3473     | READ      | Dummy argument declarations: `icrop`, `swhydrlift`, `swharvest`, `swdmmow`, `swlossmow`, `swlossgrz`, `swdmgrz`, `daylastharvest`, `maxdaymow`, `maxdaygrz`, `dmharvest`, `dmlastharvest`, `dmgrazing`, `lsdb`, `tagprest`, `crpfil` |
| 3476-3491     | READ      | Local variable declarations: `crp`, `getun2`, `i`, `nrofSeqGM`, `ini`, `count`, `ifnd`, `swcropemergence`, `swcropharvest`, `swgrazing`, `swgrazingpot`, `swanthesis`, `swIrrigate`, `daydelay`, `dnrinput`, `cfinput`, `chinput`, `laiinput`, `cfeicinput`, `depth`, `rootdis`, `tinter`, `pfree`, `pstem`, `scanopy`, `avprec`, `avevap`, `hlossmow`, `lossmow`, `hlossgrz`, `lossgrz`, `DaysGrazing`, `UptGrazing`, `LossGrazing`, `dmmowdelay`, `sum`, `swMow`, `swGrz`, `swDew`, `RDinqr`, `message`, `filnam`, `uco2`, `atmofil`, `swco2` |
| 3501          | READ      | `filnam = trim(pathcrop)//trim(crpfil)//'.crp'` — build file path                                                                                             |
| 3502          | READ      | `crp = getun2(10,90,2)` — get free file unit number                                                                                                            |
| 3503          | READ      | `call rdinit(crp,logf,filnam)` — open `.crp` file and initialise parser                                                                                        |
| 3508          | READ      | `call rdsinr('swcf',1,3,swcf)` — read crop factor / crop height switch                                                                                         |
| 3511-3515     | VALIDATE  | `if (swetr.eq.1 .and. swcf.eq.2)` fatalerr — ETref requires swcf=1 or swcf=3, not swcf=2                                                                      |
| 3517-3526     | READ      | `if (swcf.eq.1)`: read `dnr` + `cf` arrays; pack into `cftb` interleaved; set `chtb = -99.99d0`                                                               |
| 3522-3525     | RUNTIME   | Loop: pack `dnrinput`/`cfinput` into interleaved `cftb` array                                                                                                  |
| 3526          | NORMALIZE | `chtb = -99.99d0` — sentinel fill for unused crop-height table when swcf=1                                                                                     |
| 3527-3536     | READ      | `elseif (swcf.eq.2)`: read `dnr` + `ch` arrays; pack into `chtb` interleaved; set `cftb = -99.99d0` (active in both case 4 and case 2)                        |
| 3532-3535     | RUNTIME   | Loop: pack `dnrinput`/`chinput` into interleaved `chtb` array                                                                                                  |
| 3536          | NORMALIZE | `cftb = -99.99d0` — sentinel fill for unused crop-factor table when swcf=2                                                                                     |
| 3537-3557     | GUARDED   | `elseif (swcf.eq.3)`: read `lai`, `cf`, `cfeic` + `ch` arrays; pack into `cftb`, `cfeictb`, `chtb`. LAI-dependent dual crop coefficient — not in port scope   |
| 3560          | READ      | `call rdsinr('swinter',0,3,swinter)` — interception model switch                                                                                               |
| 3561-3562     | READ      | `if (swinter.eq.1)`: `call rdsdor('cofab',…)` — Von Hoyningen-Huene coefficient (active in both cases)                                                        |
| 3563-3581     | GUARDED   | `else if (swinter.eq.2)`: read `t`, `pfree`, `pstem`, `scanopy`, `avprec`, `avevap`; pack into 5 interleaved lookup tables (Gash model) — not in port scope   |
| 3582-3585     | GUARDED   | `else if (swinter.eq.3)`: read `fimin`, `siccaplai` — storage-capacity interception — not in port scope                                                        |
| 3588-3598     | NORMALIZE | `if (swcf.eq.1 .or. swcf.eq.3)`: `albedo = 0.23d0`, `rsc = 70.0d0`, `rsw = 0.0d0` — hardcoded ETref standard defaults                                       |
| 3595-3597     | READ      | `else` (swcf=2): `call rdsdor('albedo',…)`, `call rdsdor('rsc',…)`, `call rdsdor('rsw',…)` — crop-specific values (active in both cases)                      |
| 3604          | READ      | `call rdsdor('tdwi',…)` — initial total crop dry weight                                                                                                        |
| 3605          | READ      | `call rdsdor('laiem',…)` — leaf area index at emergence                                                                                                        |
| 3606          | READ      | `call rdsdor('rgrlai',…)` — maximum relative increase in LAI                                                                                                   |
| 3609          | READ      | `call rdsinr('swtsum',0,2,swtsum)` — yearly start-of-growth switch                                                                                             |
| 3610-3614     | READ      | `if (swtsum.eq.2)`: read `tsumtemp`, `tsumtime`, `tsumdepth` — soil temperature threshold for grass growth start (guarded in cases 2 and 4 which have swtsum=1) |
| 3610-3614     | GUARDED   | `if (swtsum.eq.2)`: only swtsum=0 and swtsum=1 are exercised in test cases; swtsum=2 branch present but not in Phase 3 primary test path                       |
| 3617          | READ      | `call rdador('slatb',…)` — specific leaf area as function of day number                                                                                        |
| 3618          | READ      | `call rdsdor('ssa',…)` — specific stem area                                                                                                                    |
| 3619          | READ      | `call rdsdor('span',…)` — life span of leaves                                                                                                                  |
| 3620          | READ      | `call rdsdor('tbase',…)` — lower temperature threshold for leaf ageing                                                                                         |
| 3623          | READ      | `call rdsdor('kdif',…)` — extinction coefficient for diffuse visible light                                                                                     |
| 3624          | READ      | `call rdsdor('kdir',…)` — extinction coefficient for direct visible light                                                                                      |
| 3625          | READ      | `call rdsdor('eff',…)` — light use efficiency                                                                                                                  |
| 3626          | READ      | `call rdador('amaxtb',…)` — max CO2 assimilation rate as function of day number                                                                                |
| 3627          | READ      | `call rdador('tmpftb',…)` — reduction factor of AMAX vs average temperature                                                                                    |
| 3628          | READ      | `call rdador('tmnftb',…)` — reduction factor of AMAX vs minimum temperature                                                                                    |
| 3631          | READ      | `call rdsdor('cvl',…)` — efficiency of conversion into leaves                                                                                                  |
| 3632          | READ      | `call rdsdor('cvr',…)` — efficiency of conversion into roots                                                                                                   |
| 3633          | READ      | `call rdsdor('cvs',…)` — efficiency of conversion into stems                                                                                                   |
| 3636          | READ      | `call rdsdor('q10',…)` — temperature increase in respiration rate                                                                                              |
| 3637          | READ      | `call rdsdor('rml',…)` — maintenance respiration rate of leaves                                                                                                |
| 3638          | READ      | `call rdsdor('rmr',…)` — maintenance respiration rate of roots                                                                                                 |
| 3639          | READ      | `call rdsdor('rms',…)` — maintenance respiration rate of stems                                                                                                 |
| 3640          | READ      | `call rdador('rfsetb',…)` — reduction factor of senescence as function of day number                                                                           |
| 3643          | READ      | `call rdador('frtb',…)` — fraction partitioned to roots as function of day number                                                                              |
| 3644          | READ      | `call rdador('fltb',…)` — fraction partitioned to leaves as function of day number                                                                             |
| 3645          | READ      | `call rdador('fstb',…)` — fraction partitioned to stems as function of day number                                                                              |
| 3648          | READ      | `call rdsdor('perdl',…)` — maximum relative death rate of leaves                                                                                               |
| 3649          | READ      | `call rdador('rdrrtb',…)` — relative death rates of roots as function of day number                                                                            |
| 3650          | READ      | `call rdador('rdrstb',…)` — relative death rates of stems as function of day number                                                                            |
| 3653          | NORMALIZE | `swoxygen = 1` — default: oxygen stress according to Feddes                                                                                                    |
| 3654-3655     | READ      | `if (rdinqr('swoxygen'))`: `call rdsinr('swoxygen',0,2,swoxygen)` — read oxygen stress switch (case 4 has swoxygen=2; case 2 has swoxygen=1)                   |
| 3657-3662     | VALIDATE  | `if (swoxygen.eq.2 .and. (swhea.eq.0 .or. swcalt.eq.1))` fatalerr — physical O2 stress needs numerical heat flow                                              |
| 3664-3669     | VALIDATE  | `if (swoxygen.eq.2 .and. bdens(1).lt.100)` fatalerr — physical O2 stress needs realistic bulk density                                                         |
| 3672-3677     | READ      | `if (swoxygen.eq.1)`: read `hlim1`, `hlim2u`, `hlim2l` — Feddes anaerobiosis pressure heads (active in case 2)                                                |
| 3679-3712     | READ      | `if (swoxygen.eq.2)`: Bartholomeus oxygen stress block (active in case 4)                                                                                      |
| 3682          | NORMALIZE | `swoxygentype = 1` — default: physical processes                                                                                                               |
| 3683-3685     | READ      | `if (rdinqr('swoxygentype'))`: `call rdsinr('swoxygentype',1,2,swoxygentype)` — read Bartholomeus type switch (case 4 has swrootradius=2, swoxygentype=1)       |
| 3687-3704     | READ      | `if (swoxygentype.eq.1)`: read `q10_microbial`, `specific_resp_humus`, `srl`, `swrootradius`; conditional on `swrootradius`: read dry_mat_cont_roots / air_filled_root_por / spec_weight_root_tissue / var_a (swrootradius=1) or `root_radiusO2` (swrootradius=2) — all supported in Phase 3 |
| 3705-3711     | GUARDED   | `else` (swoxygentype=2): read `SwTopSub`, `NrStaring`; `call oxygen_dat(…)` — reproduction functions — not in port scope                                      |
| 3715          | NORMALIZE | `swWrtNonox = 0` — default: no non-oxidative root senescence switch                                                                                            |
| 3716-3718     | READ      | `if (rdinqr('swWrtNonox'))`: `call rdsinr('swWrtNonox',0,1,swWrtNonox)` — both cases have swwrtnonox=1                                                        |
| 3720          | NORMALIZE | `aeratecrit = 0.0001d0` — default aeration criterion                                                                                                           |
| 3721-3723     | READ      | `if (swWrtNonox.eq.1)`: `call rdsdor('aeratecrit',…)` — both cases have aeratecrit=0.5                                                                        |
| 3726          | NORMALIZE | `swdrought = 1` — default: drought stress according to Feddes                                                                                                  |
| 3727-3729     | READ      | `if (rdinqr('swdrought'))`: `call rdsinr('swdrought',1,2,swdrought)` — both cases have swdrought=1                                                             |
| 3731-3738     | READ      | `if (swdrought.eq.1)`: read `hlim3h`, `hlim3l`, `hlim4`, `adcrh`, `adcrl` — Feddes drought parameters (active in both cases)                                  |
| 3740-3754     | GUARDED   | `else` (swdrought=2): read De Jong van Lier params: `wiltpoint`, `kstem`, `rxylem`, `rootradius`, `kroot`, `rootcoefa`, `swhydrlift`, `rooteff`, `stephr`, `criterhr`, `taccur` — not in port scope |
| 3757-3778     | GUARDED   | `if (flsolute)` salt stress block: both cases have swsalinity=0, so entire block is guarded; even if flsolute=.true., swsalinity=0 means no reads fire. Schema must expose swsalinity=0 only |
| 3763-3766     | GUARDED   | `if (swsalinity.eq.1)`: read `saltmax`, `saltslope` — not in port scope                                                                                        |
| 3767-3777     | GUARDED   | `elseif (swsalinity.eq.2)`: read `salthead`; fatalerr if swdrought=1 — not in port scope                                                                       |
| 3781          | NORMALIZE | `swcompensate = 0` — default: no compensation                                                                                                                  |
| 3782          | NORMALIZE | `swjarvis = 0` — legacy default                                                                                                                                 |
| 3783-3784     | READ      | `if (rdinqr('swcompensate'))`: `call rdsinr('swcompensate',0,2,swcompensate)` — case 4 has swcompensate=1; case 2 has swcompensate=0                           |
| 3785-3798     | READ      | `else if (rdinqr('swjarvis'))` legacy shim: warn, read `swjarvis`, coerce swcompensate=1 if swjarvis>0                                                         |
| 3787-3788     | VALIDATE  | `call warn(…)` — SWJARVIS is deprecated; migrate to SWCOMPENSATE                                                                                               |
| 3792-3793     | VALIDATE  | `call warn(…)` — SWJARVIS is applied to all stresses (when swjarvis≠4)                                                                                         |
| 3795          | NORMALIZE | `swcompensate = 1` — coerce legacy swjarvis > 0 to new swcompensate=1                                                                                          |
| 3801-3806     | VALIDATE  | `if (swcompensate.gt.0 .and. swdrought.eq.2)` fatalerr — compensation not allowed with De Jong van Lier drought                                                |
| 3809          | NORMALIZE | `swstressor = 1` — default: compensate all stressors                                                                                                           |
| 3810-3813     | READ      | `if (swcompensate.gt.0 .and. rdinqr('swstressor'))`: `call rdsinr('swstressor',1,5,swstressor)` — case 4 has swstressor=1 (default)                           |
| 3817-3821     | READ      | `if (swcompensate.eq.1)`: `alphacrit = 1.0d0`; `call rdsdor('alphacrit',…)` — Jarvis compensation; case 4 has alphacrit=0.7 (active in Phase 3)               |
| 3820          | NORMALIZE | `alphacrit = 1.0d0` — default critical stress index (inside swcompensate=1 guard)                                                                              |
| 3824-3830     | GUARDED   | `elseif (swcompensate.eq.2)`: `dcritrtz = 0.0d0`; `call rdsdor('dcritrtz',…)` — Walsum compensation — not in port scope                                       |
| 3835-3838     | NORMALIZE | `swrdc = 0`; `if (rdinqr('swrdc'))`: `call rdsinr('swrdc',0,1,swrdc)` — both cases have swrdc=0 (default)                                                     |
| 3841          | READ      | `call rdador('rdctb',0,100,rdctb,22,ifnd)` — root density distribution table (always read)                                                                     |
| 3844          | NORMALIZE | `swrd = 3` — default: root extension depends on available root biomass (differs from cropfixed default=1 and cropwofost default=2)                              |
| 3845-3847     | READ      | `if (rdinqr('swrd'))`: `call rdsinr('swrd',1,3,swrd)` — both cases have swrd=2                                                                                 |
| 3850-3852     | GUARDED   | `if (swrd.eq.1)`: `call rdador('rdtb',…)` — rooting depth vs day number table — guarded in both cases (swrd=2)                                                 |
| 3855-3865     | READ      | `elseif (swrd.eq.2)`: read `rdi`, `rri`, `rdc`; optional `swdmi2rd` — max-daily-increase root extension (active in both cases)                                 |
| 3862-3865     | NORMALIZE | `swdmi2rd = 0`; `if (rdinqr('swdmi2rd'))`: `call rdsinr('swdmi2rd',0,1,swdmi2rd)` — both cases have swdmi2rd=1                                                 |
| 3868-3873     | GUARDED   | `elseif (swrd.eq.3)`: read `rlwtb`, `wrtmax` — root biomass-based extension — not in port scope (both cases have swrd=2)                                       |
| 3876-3879     | NORMALIZE | `relmf = 1.0d0`; `if (rdinqr('relmf'))`: `call rdsdor('relmf',…)` — both cases have relmf=0.90                                                                |
| 3882-3885     | NORMALIZE | `swpotrelmf = 1`; `if (rdinqr('swpotrelmf'))`: `call rdsinr('swpotrelmf',1,2,swpotrelmf)` — both cases have swpotrelmf=1 (default)                            |
| 3891          | READ      | `call rdainr('SeqGrazMow',1,3,SeqGrazMow,366,nrofSeqGM)` — read yearly sequence of mowing/grazing periods (integer array, both cases all-2)                    |
| 3893-3905     | RUNTIME   | Loop: set `swMow`, `swGrz`, `swDew` booleans from SeqGrazMow values                                                                                            |
| 3910-3982     | GUARDED   | `if (swGrz)` entire grazing settings block — guarded because both cases have all-mowing sequences (SeqGrazMow=all-2, swGrz=.false.); includes swharvest read (grazing), swdmgrz, dmgrazing, dmgrztb, maxdaygrz, swlossgrz, zgrz+lossgrztab, dewrest, tagprest, LSDb/daysgrazing/uptgrazing/lossgrazing tables, LSDa — all guarded |
| 3910-3911     | READ      | `call rdsinr('swharvest',1,2,swharvest)` — note this is the grazing copy; entirely guarded here since swGrz=.false. in both test cases                         |
| 3910-3982     | GUARDED   | `seqgrazmow(i) ∈ {1,3}` guards entire block; parser concern noted — SeqGrazMow is an integer array read with `rdainr`, not `rdsinr`; schema must handle array-of-int differently from scalar-int |
| 3985-3988     | READ      | `if (swMow)`: `call rdsdor('mowrest',…)` — remaining yield after mowing (active in both cases)                                                                 |
| 3991          | READ      | `call rdsinr('swharvest',1,2,swharvest)` — mowing copy; this is the active one for both test cases; case 4 has swharvest=1, case 2 has swharvest=2             |
| 3992-4003     | READ      | `if (swharvest.eq.1)`: read `swdmmow`; if swdmmow=1 read `dmharvest`, `daylastharvest`, `dmlastharvest`; if swdmmow=2 read `dmmowtb`, `maxdaymow` — case 4 has swdmmow=2, case 2 has swharvest=2 (bypasses this branch) |
| 3996-3999     | READ      | `if (swdmmow.eq.1)`: read `dmharvest`, `daylastharvest`, `dmlastharvest` — case 4 has swdmmow=2, so guarded there; active for swdmmow=1 scenario               |
| 3996-3999     | GUARDED   | `swdmmow=1` path guarded in case 4 (swdmmow=2 there); but swdmmow=1 IS a supported Phase 3 variant per spec (DM-threshold fixed)                               |
| 4000-4003     | READ      | `elseif (swdmmow.eq.2)`: read `dmmowtb`, `maxdaymow` — case 4 active path                                                                                     |
| 4006-4017     | READ      | `call rdsinr('swlossmow',0,1,swlossmow)` — both cases have swlossmow=0; `if (swlossmow.eq.1)` block for zmow+lossmowtab guarded                               |
| 4007-4017     | GUARDED   | `if (swlossmow.eq.1)`: read `zmow`, `hlossmow`, `lossmow`; pack into `lossmowtab` — treading losses during mowing — not in port scope                         |
| 4019-4025     | READ      | `elseif (swharvest.eq.2)`: `call rdatim('dateharvest',dateharvest,999,ifnd)` — read fixed harvest dates; `dateharvest(ifnd+1) = tend+1.d0` sentinel fill; case 2 active path |
| 4023          | RUNTIME   | `dateharvest(ifnd + 1) = tend + 1.d0` — sentinel fill past end of dates; this is a runtime side-effect, not a pure read                                        |
| 4028-4029     | READ      | `call rdainr('daydelay',…)`, `call rdador('dmmowdelay',…)` — delay of regrowth after mowing (both cases)                                                       |
| 4031-4035     | RUNTIME   | Loop: pack `daydelay`/`dmmowdelay` into interleaved `DelayRegrowthTab` array                                                                                   |
| 4040          | READ      | `call rdsinr('schedule',0,1,schedule)` — irrigation scheduling switch (always read; both cases have schedule=0)                                                 |
| 4042-4047     | GUARDED   | `if (schedule.eq.1 .and. swdrought.eq.2)`: re-read `hlim3h`, `hlim3l`, `hlim4` for scheduling override — guarded (schedule=0 in both cases)                   |
| 4050          | NORMALIZE | `flco2 = .false.` — default: no CO2 assimilation correction                                                                                                    |
| 4051          | NORMALIZE | `swco2 = 0` — default CO2 switch                                                                                                                               |
| 4052-4054     | READ      | `if (rdinqr('swco2'))`: `call rdsinr('swco2',0,1,swco2)` — both cases have swco2=0                                                                             |
| 4055          | NORMALIZE | `if (swco2.eq.1) flco2 = .true.` — derive flco2 from swco2                                                                                                     |
| 4058-4083     | GUARDED   | `if (flco2)` block: read `CO2AMAXTB`, `CO2EFFTB`, `CO2TRATB`; optionally read `atmofil`; build filnam for `.co2` file; `close(crp)`; open `.co2` file; read `CO2year`, `CO2ppm`; `close(uco2)` — entire CO2 block guarded (both cases have swco2=0) |
| 4079-4082     | READ      | `else`: `close(crp)` — close `.crp` file in the non-CO2 path                                                                                                   |
| 4086-4089     | GUARDED   | `if (schedule.eq.1)`: `call irrigation(1)` + optional `IrrigationOutput(1)` — scheduling initialisation — guarded (schedule=0)                                 |
| 4091-4121     | RUNTIME   | `if (swdrought.eq.1)` block: build `cumdens` array — normalised cumulative root density function over 101 depth points (3 loops; depends only on `rdctb`)       |
| 4096-4100     | RUNTIME   | Inner loop i=0,100: populate `rootdis` using `afgen(rdctb,22,depth)` at 0.01-spaced depths                                                                     |
| 4103-4114     | RUNTIME   | Cumulative density computation: copy depths to odd indices; trapezoidal integration loop i=4,202,2; accumulate `sum`, store in `cumdens` even indices           |
| 4117-4119     | RUNTIME   | Normalization loop i=2,202,2: `cumdens(i) = cumdens(i) / sum`                                                                                                  |
| 4125-4267     | GUARDED   | `.END`-file restart block: `if (t1900-tstart.lt.1d-3 .and. swinco.eq.3 .and. .not. abs(t1900-cropstart(icrop)).lt.1d-3)` — entire block guarded (swinco=3 only) |
| 4128-4131     | GUARDED   | `if (flCropCalendar)`: `ini = getun2(…)`; `call rdinit(ini,logf,inifil)` — open `.END` file                                                                    |
| 4133-4140     | GUARDED   | `if (rdinqr('swcropharvest'))`: read `swcropharvest`; derive `flCropHarvest`                                                                                    |
| 4142-4147     | GUARDED   | Read `swCropEmergence`; derive `flCropEmergence`                                                                                                                |
| 4149-4161     | GUARDED   | `if (.not. flCropHarvest)`: read `swIrrigate`; derive `flIrrigate`; conditionally read `dayfix`                                                                 |
| 4163-4248     | GUARDED   | Read full grass crop state variables from `.END` file: `rd`, `rdpot`, `dvs`, `daycrop`, `swanthesis`, `tsum`, `ilvold`, `ilvoldpot`, `wrt`, `wrtpot`, `tadw`, `tadwpot`, `wst`, `wstpot`, `wso`, `wsopot`, `wlv`, `wlvpot`, `laiexp`, `laiexppot`, `glaiex`, `glaiexpot`, `lai`, `laipot`, `dwrt`, `dwrtpot`, `dwlv`, `dwlvpot`, `dwst`, `dwstpot`, `gasst`, `gasstpot`, `mrest`, `mrestpot`, `cwdm`, `cwdmpot`, `nofd`, `atmin7`, `rid`, `idregr`, `idregrpot`, `laimax`, `tagp`, `tagppot`, `tagpt`, `tagptpot`, `daygrowth`, `daygrowthpot`, `cuptgraz`, `cuptgrazpot`, `tsum` (again), `swcropharvest`, `iseqgm`, `iseqgmpot`, `swgrazing`, `swgrazingpot`, `iharvest`, `idaysgraz`, `idaysgrazpot`, arrays `sla`, `slapot`, `lvage`, `lvagepot`, `lv`, `lvpot` |
| 4234-4248     | GUARDED   | Derive `flanthesis`, `flgrazing`, `flgrazingpot` from swanthesis/swgrazing/swgrazingpot integers                                                                 |
| 4250-4261     | GUARDED   | `else` branches (flCropHarvest=true or rdinqr('swcropharvest')=false): set `daycrop=0`, `flCropEmergence=.false.` as defaults                                  |
| 4265          | GUARDED   | `close(ini)` — close `.END` file                                                                                                                                |
| 4269          | READ      | `return`                                                                                                                                                        |
| 4270          | READ      | `end` — end of subroutine                                                                                                                                       |

---

## Summary

### Lines per bucket

| Bucket    | Approx. line count | Description                                                                                                                                    |
| --------- | ------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| READ      | ~270 lines         | `rd*` calls, file open/close, declarations, use statements; mandatory reads for all grass biophysics + management parameters                   |
| VALIDATE  | ~20 lines          | `fatalerr` / `warn` calls: swcf+swetr, swoxygen=2+heat/bdens, swcompensate+swdrought, swsalinity=2+swdrought, swjarvis deprecation warning     |
| NORMALIZE | ~25 lines          | Defaults: swoxygen=1, swWrtNonox=0, aeratecrit=1e-4, swdrought=1, swcompensate=0, swjarvis=0, swstressor=1, alphacrit=1.0, swrdc=0, swrd=3, relmf=1.0, swpotrelmf=1, swdmi2rd=0, flco2=false, swco2=0; plus chtb/cftb sentinel fills, albedo/rsc/rsw ETref defaults |
| RUNTIME   | ~35 lines          | cumdens build (3 loops); cftb/chtb packing loops (swcf=1 and swcf=2); swMow/swGrz/swDew Boolean derivation; DelayRegrowthTab packing; dateharvest sentinel fill |
| GUARDED   | ~410 lines         | swcf=3, swinter=2/3, swtsum=2 (not exercised by test cases), swoxygen=2+swoxygentype=2, swdrought=2, swsalinity≠0 (entire flsolute block), swcompensate=2, swrd=1, swrd=3, swGrz block (entire grazing settings), swlossmow=1, swlossgrz=1, seqgrazmow(i)∈{1,3}, schedule=1, swco2=1, swinco=3 warm-restart |

Total: 834 lines of substantive content (3437–4270), consistent with plan's 3437-4270 range — **no drift**.

### Call count summary

- **READ:** ~75 `rd*` calls (rdsinr, rdsdor, rdador, rdfdor, rdainr, rdatim, rdinqr, rdinit, rdscha) plus file operations (getun2×2, rdinit×2, close×2)
- **VALIDATE:** 6 `fatalerr` calls + 2 `warn` calls
- **NORMALIZE:** ~20 distinct default assignments
- **RUNTIME:** 4 initialization actions: swMow/swGrz/swDew derivation loop; cftb/chtb packing loops; DelayRegrowthTab packing loop; dateharvest sentinel fill; cumdens build (3 loops)
- **GUARDED:** 17 conditional branches

---

## Scope for Phase 3

### Active switch values (both cases exercise)

| Switch           | Case 4 value      | Case 2 value     | Read location |
| ---------------- | ----------------- | ---------------- | ------------- |
| swcf             | 2                 | 2                | 3508          |
| swinter          | 1                 | 1                | 3560          |
| swtsum           | 1                 | 1                | 3609          |
| swoxygen         | 2                 | 1                | 3654-3655     |
| swoxygentype     | 1 (default)       | N/A              | 3683-3685     |
| swrootradius     | 2                 | N/A              | 3693          |
| swWrtNonox       | 1                 | 1                | 3716-3718     |
| aeratecrit       | 0.5               | 0.5              | 3721-3723     |
| swdrought        | 1 (default)       | 1 (default)      | 3727-3729     |
| swsalinity       | 0 (not read)      | 0 (not read)     | 3757-3778 (guarded) |
| swcompensate     | 1                 | 0 (default)      | 3783-3784     |
| swstressor       | 1 (default)       | N/A              | 3810-3813     |
| alphacrit        | 0.7               | N/A              | 3821          |
| swrdc            | 0 (default)       | 0 (default)      | 3835-3838     |
| swrd             | 2                 | 2                | 3845-3847     |
| swdmi2rd         | 1                 | 1                | 3862-3865     |
| relmf            | 0.90              | 0.90             | 3876-3879     |
| swpotrelmf       | 1 (default)       | 1 (default)      | 3882-3885     |
| SeqGrazMow       | all-2 (mow only)  | all-2 (mow only) | 3891          |
| swharvest        | 1 (DM-threshold)  | 2 (fixed dates)  | 3991          |
| swdmmow          | 2 (flexible)      | N/A (swharvest=2)| 3995          |
| swlossmow        | 0                 | 0                | 4006          |
| schedule         | 0                 | 0                | 4040          |
| swco2            | 0 (default)       | 0 (default)      | 4052-4054     |

### Guarded (stub-errored in port)

swcf=3, swinter=2/3, swoxygentype=2, swdrought=2, swsalinity≠0, swcompensate=2, swrd=1, swrd=3 (swrd default in readgrass is 3, but both test cases explicitly set swrd=2), swlossmow=1, swlossgrz=1, seqgrazmow(i)∈{1,3} (any grazing or dewooling), schedule=1, swco2=1, swinco=3 warm-restart.

### Notable grass-specific differences from cropfixed/cropwofost

1. **swrd default is 3** (root biomass-based), not 1 (cropfixed) or 2 (cropwofost). Both test cases explicitly override to swrd=2. The TOML schema should default swrd to 3 but both active test cases provide explicit swrd=2 values.
2. **SeqGrazMow is an integer array** read with `rdainr`, not `rdsdor`/`rdsinr`. This requires special handling in the TOML schema — it maps to an array of integers, not a single switch. Both test cases have all-20 entries equal to 2 (mow only). The parser concern is that `rdainr` accepts a maximum of 366 values; the TOML schema should use a `[integer]` array field.
3. **swharvest is read twice** inside readgrass: once inside the `if (swGrz)` block (grazing) and once inside the `if (swMow)` block (mowing). Only the mowing copy fires for the test cases. The TOML schema needs only one `swharvest` field since the two are semantically the same switch — the legacy reader duplicates the read because mowing and grazing have separate sections.
4. **dateharvest is an array of dates** read with `rdatim` (not `rdador`). Phase 3 must handle ISO-date parsing for swharvest=2. The trailing sentinel `dateharvest(ifnd+1) = tend+1.d0` is a RUNTIME action.
5. **No `cvo`, `fotb`, `rmo`** fields — grass has no storage organ, so these wofost/cropfixed fields are absent. The schema must not inherit them.
6. **DelayRegrowthTab** is packed from two parallel arrays (`daydelay` integer array + `dmmowdelay` real array) into an interleaved real lookup table. This is a RUNTIME packing step analogous to cftb/chtb in the other crop types.

---

## Sibling Readers

This section documents every `.crp` open reachable from the grass runtime path, per ADR 0017 (sibling-reader dispatch pattern).

### ADR 0017 grep results

Command run:
```
grep -n "trim(pathcrop)\|'.crp'" src/io/readswap.f90 src/crop/*.f90 \
  | grep -v "outfil\|outputfile" | sort -u
```

All `.crp` opens in the codebase:

| File                          | Line | Open                                            | Guard                           |
| ----------------------------- | ---- | ----------------------------------------------- | ------------------------------- |
| `src/io/readswap.f90`         | 2075 | `readcropfixed` primary `.crp` open             | croptype=1                      |
| `src/io/readswap.f90`         | 2583 | `readwofost` primary `.crp` open                | croptype=2                      |
| `src/io/readswap.f90`         | 3273 | `readarablelandgerm` `.crp` open                | croptype=1 or 2 (see below)     |
| `src/io/readswap.f90`         | 3501 | `readgrass` primary `.crp` open                 | croptype=3                      |
| `src/io/readswap.f90`         | 4066 | `readgrass` CO2 `.co2` open (inside swco2=1)    | swco2=1 only — guarded          |
| `src/crop/cropgrowth.f90`     | 1109 | N-P-K nutrient block `.crp` open                | `if (flCropNut)` — wofost path only |
| `src/crop/irrigation.f90`     | 77   | `irrigation(1)` `.crp` open                     | `if (schedule.eq.1)` — guarded  |

### Sibling analysis for the grass (croptype=3) path

**`readarablelandgerm` — NOT on the grass path.**

The `ArableLandGerm` subroutine is dispatched from `cropgrowth.f90` at lines 90–178. The dispatch fires only when `.not. flCropEmergence` — i.e., during the pre-emergence stage. For croptype=3 (grass), `cropgrowth.f90` calls `Grass(1)` at line 195 inside the `if (croptype(icrop) .eq. 3)` branch, **after** the ArableLandGerm dispatch block. The ADR 0017 cache-hit dispatch at lines 92–157 covers only rotation_type=1 (cropfixed) and rotation_type=2 (wofost) — there is no croptype=3 case in that select statement. The legacy fallback `call ArableLandGerm(1)` at line 155 would fire for grass if `use_cache=.false.` and `flCropReadFile=.true.`, but the `if (.not. flCropEmergence)` guard at line 90 means it only fires during pre-emergence stages that grass does not use. In practice, grass skips the ArableLandGerm block entirely because `flCropEmergence` is set to `.true.` for grass before this point (grass uses a different emergence logic).

Conclusion: `readarablelandgerm` is NOT reachable from the grass runtime path. **No new ADR 0017 dispatch is needed for grass beyond what Phases 1 and 2 already wired.**

**N-P-K block (`cropgrowth.f90:1109`) — NOT on the grass path.**

The N-P-K block at lines 1107–1134 is inside the `wofost(task=1)` subroutine. The `Grass(task=1)` subroutine is a separate subroutine at line 2115. The `flCropNut` check at `cropgrowth.f90:312` (in the main crop growth task=2 path) applies to `croptype(icrop).ge.2`, but the N-P-K initialization at line 1107 is inside wofost's `case (1)` select block — not reachable from `grass`. Both test cases have `flCropNut=.false.` in any event.

**`irrigation(1)` — guarded by schedule=1.** Both cases have schedule=0; this sibling is guarded and out of Phase 3 scope.

**Summary:** The grass runtime path opens **exactly one** `.crp` file — the primary open inside `readgrass` at line 3501. No new ADR 0017 dispatch wiring is required for Phase 3.

---

## Implementation Guidance

### 1. swrd default=3 is unique to grass

`readgrass` sets `swrd = 3` as default (line 3844), whereas `readcropfixed` sets `swrd = 1` and `readwofost` sets `swrd = 2`. The TOML schema must document `swrd` default as 3 for cropgrass — but since both test cases explicitly set `swrd = 2`, the active port path is swrd=2 (rdi/rri/rdc/swdmi2rd). The swrd=3 path (rlwtb/wrtmax) is guarded.

### 2. SeqGrazMow parser concern

`SeqGrazMow` is read with `rdainr` (integer array, max 366 values). This is unlike any scalar switch in the previous two phases. The TOML schema must represent it as an array of integers. The Phase 3 validator must check: (a) length ≥ 1; (b) all values ∈ {1,2,3}; (c) stub-error if any value ≠ 2 (grazing and dewooling not supported). Both test cases have exactly 20 entries all equal to 2.

### 3. swharvest dual-read is a single schema field

The legacy reader reads `swharvest` twice (once in the grazing block, once in the mowing block). The TOML schema exposes it as a single field. Since the grazing read is entirely guarded in Phase 3, the validator only needs to validate the mowing swharvest.

### 4. dateharvest requires ISO-date array parsing

When `swharvest=2`, `dateharvest` is a variable-length array of calendar dates read with `rdatim`. Case 2 has 30 dates spanning 1980-1984. The TOML schema must represent this as an array of date strings (or date values). The trailing sentinel `dateharvest(ifnd+1) = tend+1.d0` is a RUNTIME side-effect that `cropgrass_init` must reproduce.

### 5. swoxygen=2 supported path is swoxygentype=1 only

The `swoxygentype=1` sub-path (physical Bartholomeus with `q10_microbial`, `specific_resp_humus`, `srl`, `swrootradius`) is fully supported. The `swoxygentype=2` sub-path (reproduction functions with `SwTopSub`, `NrStaring`, `oxygen_dat`) is guarded. Case 4 uses `swrootradius=2` (root_radiusO2 given directly), so the `swrootradius=2` path is also active.

### 6. swcompensate=1 is supported (Jarvis) — swcompensate=2 is guarded (Walsum)

Case 4 has `swcompensate=1` with `alphacrit=0.7`. Case 2 has `swcompensate=0`. Both are in scope. The `swcompensate=2` path (dcritrtz) is guarded.

### 7. cumdens post-close, swdrought=1 only

Identical pattern to cropfixed and cropwofost. Belongs in `cropgrass_init`. Only computed when swdrought=1 (which is the case for both test cases).

### 8. swjarvis legacy shim present (identical to other crop types)

Lines 3785-3798 replicate the same swjarvis deprecation shim as in readcropfixed and readwofost. The TOML schema should only expose `swcompensate`; swjarvis must not be replicated.

### 9. The .END-file restart block is structurally similar to wofost but with grass-specific fields

The grass `.END`-file block (lines 4125-4267) reads ~50 state variables including grass-specific ones (`rid`, `idregr`, `idregrpot`, `tagp`, `tagppot`, `tagpt`, `tagptpot`, `daygrowth`, `daygrowthpot`, `cuptgraz`, `cuptgrazpot`, `iseqgm`, `iseqgmpot`, `swgrazing`, `swgrazingpot`, `iharvest`, `idaysgraz`, `idaysgrazpot`). The swinco=3 path should be refused by the config validator for Phase 3.
