# Phase 4f — salinitystress exact-parity audit

Phase 4f Task B5: walk every legacy `.swp` and `.dra` key for case 5
(`tests/swap-cases/5.salinitystress/`) and classify against the TOML pipeline
(`src/config/*_config.f90`, `src/io/toml/read_*_toml.f90`,
`src/io/toml/config_to_variables.f90`).

Status legend (mirrors hupselbrook / grassgrowth / oxygenstress audits):
- OK — TOML authors it, reader copies it, adapter writes the legacy global,
  no finalize mismatch.
- TOML — schema slot exists but the case TOML is missing the value.
- READER — schema field exists but no reader populates it.
- ADAPTER — schema + reader fine, but adapter never copies it into the
  legacy global.
- FINALIZE — adapter copies the value but skips a unit-conversion / derived
  finalize step the legacy reader performs.
- N/A — switch is read but its value is irrelevant (gated branch not taken).

Baseline (commit a8cee6d): `pixi run -e test regression salinitystress`
aborts with `psand/psilt/pclay/porg must all be set together` — heat
validator fires because the TOML authors `psand/pclay/porg` but not
`psilt`. (Same hard-block as oxygenstress before B4.)

## .swp — General + Output

| Legacy key      | Status        | Notes |
|-----------------|---------------|-------|
| project, paths  | OK            | adapter line 59-63 (`saltfarmtexel`). |
| swscre, swerror | OK            | both 0. |
| tstart / tend   | OK            | 2012-01-01 .. 2015-12-31. |
| nprintday       | OK            | =1. |
| swmonth         | OK            | =0; PERIOD/SWRES/SWODAT path. |
| period, swres, swodat | OK      | period=1, swres=0, swodat=0. |
| swyrvar, datefix | OK/N/A       | swbal=0 → no .bal output → N/A. |
| outfil          | OK (HACK)     | hard-coded `'result'` in adapter. |
| swheader, swwba/swend/swvap/swbal/swblc/swsba/swate/swbma/swdrf/swswb/swini/swinc/swcrp/swstr/swirg | OK (HACK) | RETIRED zero-forced. |
| swcsv           | OK (HACK)     | hard-coded =1 in adapter. |
| **inlist_csv**  | **TOML missing** | case 5 lists 6 columns (`cwso,cpwso,tredwet,treddry,tredsol,conc[-5.0,-25.0,-55.0]`). The adapter falls back to the hupselbrook 12-col water-balance default. The fixture asserts on `CWSO/CPWSO/TREDDRY/TREDWET/TREDSOL/CONC[-5.0]/CONC[-25.0]/CONC[-55.0]` — without the override **the regression aggregator will receive missing columns**. **HARD BLOCK on regression compare.** |
| **swcsv_tz, inlist_csv_tz** | **ADAPTER (HACK) — wrong default** | adapter HACKs `swcsv_tz = 0`. Case 5 sets `swcsv_tz=1` with `InList_csv_tz = 'conc'`. Legacy `swcsv_tz=1` writes a per-(time, depth) profile CSV; case 5 fixture does not consume the TZ output (only the regular `_output.csv`), so this drift is **N/A** for fixture compare but worth noting. |
| swafo, swaun    | OK (HACK)     | zero-forced. |

## .swp — Meteorology

| Legacy key      | Status   | Notes |
|-----------------|----------|-------|
| metfil          | OK       | `235.met`. |
| lat, alt, altw  | OK       | 52.928 / 1.2 / 10.0. |
| swetr           | OK       | =0 (PM). |
| angstroma/b     | OK       | 0.25 / 0.5. |
| swdivide        | OK       | =1 → PM-direct. |
| swmetdetail     | OK       | =0. |
| swetsine        | OK       | =0. |
| swrain          | OK       | =2 → daily + duration. |
| rainfil, rainflux | N/A    | swrain != 3. |

## .swp — Crop section

| Legacy key      | Status        | Notes |
|-----------------|---------------|-------|
| swcrop          | OK            | =1; adapter arms flCropReadFile/Open. |
| cropstart, cropend, cropfil, croptype | OK | 4 `potatod` rotations 2012..2015 type=2 (WOFOST general). |
| **rds (rdmax)** | OK (HACK)     | Adapter HACKs rdmax=200. Case 5 .swp authors **100.0 cm**. The HACK 200 value is permissive (max rooting depth used as a cap); WOFOST `rd` per-day capped by `min(afgen(rdtb,...), rdm)`. potatod.crp's RDC max is 100 cm so rd is the binding constraint, not rdm. **No water-balance impact** (verified via grassgrowth/oxygenstress where the same 200 HACK applied without harm). |

## .swp — Irrigation

| Legacy key      | Status        | Notes |
|-----------------|---------------|-------|
| swirfix         | TOML/OK       | TOML authors swirfix=1 (matches case). |
| **swirgfil**    | **READER missing** | Case 5 sets swirgfil=1 (external `swap.irg` file with daily 14.4 mm gifts 2012-04-17..2015-08-18). Legacy reads `swirgfil` then opens `pathwork//irgfil//.irg` (readswap.f90:1575). The adapter has no swirgfil schema slot or wiring; existing irrigation_config_t exposes `irgfil` (string) only. The legacy `.irg` reader (readswap.f90:1573-1593) is a separate `RDinit` open of the .irg file. Without swirgfil=1 the global stays 0 → adapter never reads .irg → no irrigation events. **HARD BLOCK on irrigation gifts.** |
| **irgfil = 'swap'** | **TOML/ADAPTER missing** | TOML authors `irgfil="swap"` but no adapter wiring exists for swirgfil=1 dispatch. |

> Note: the legacy `.irg` external file is read by `readswap.f90` lines 1573-1620
> via `RDinit(...).//.irg`, which then calls `rdatim(irdate)`, `rdfdor(irdepth)`,
> `rdfdor(irconc)`, `rdfinr(irtype)` — exactly the same parsers as the inline
> branch. The adapter could either (a) parse `.irg` itself, (b) defer to the
> legacy reader by leaving swirgfil=1 in the global so readswap's gated read
> path stays alive (but the strangler routes around readswap entirely), or
> (c) author the gifts inline as `[[irrigation.fixed_events]]`. **Option (c)
> is the strangler-clean path — translate the .irg into a TOML table.** This
> requires authoring the date/depth/conc/type rows in TOML. The .irg file
> spans 2012-2015 ~80-100 daily entries during each potato season.

## .swp — Soil water

| Legacy key      | Status        | Notes |
|-----------------|---------------|-------|
| swinco          | OK            | TOML uses swinco=2 (gwli=-90) as schema-clean stand-in for legacy swinco=3 (INIFIL). Without `swap.ini` parsing the .ini's restart state is dropped; the TOML's gwli=-90 + pondini=0 is the closest analytic equivalent. **Initial-state mismatch may bias year-1 numbers slightly but should converge.** |
| gwli            | OK            | =-90.0 cm. |
| swpondmx        | OK            | =0; default. |
| pondmx          | OK            | 0.2. |
| **rsro**        | **TOML missing** | legacy reads 0.5 d. Schema slot exists; case TOML doesn't author it. Initialize zeroes → instantaneous-runoff branch. (Same as hupselbrook / grassgrowth / oxygenstress.) |
| **rsroexp**     | **TOML missing** | legacy reads 1.0. Same gap. |
| swrunon         | OK            | =0; flrunon=.false. |
| **cfevappond**  | OK (HACK)     | adapter HACKs 1.25; matches case. |
| swcfbs / cfbs   | OK            | swcfbs=0. |
| **rsoil**       | **TOML missing** | case authors **150.0 s/m**. Without it default 0 → EPOT direct contributor. |
| swredu          | OK (HACK)     | adapter HACKs swredu=1; **case 5 .swp authors swredu=2** (Boesten/Stroosnijder). Legacy reads `cofredbo=0.54` under swredu=2. Adapter currently maps either cofredbl OR cofredbo to legacy `cofred`. Need to ensure swredu=2 path is exercised. **EVAPORATION drift.** |
| **cofredbo**    | **TOML missing** | case authors 0.54; under swredu=2 the legacy reads it. The schema's `meteo.evaporation.cofredbo` slot is what the adapter uses (line 180-183). Without authoring it falls back to default 0.35 (cofredbl path). |
| **rsigni**      | OK (HACK)     | adapter HACKs 0.5; matches default. |
| **isublay/isoillay/hsublay/ncomp** | **TOML missing** | 5-row sub-layer table spanning 2 soil-physical layers (total = 30+30+40+50+850 = 1000 cm, 30+30+40+10+85 = 195 compartments). |
| swsophy         | OK            | =0 (analytical MvG). |
| **MvG params (ores, osat, alfa, npar, ksatfit, lexp, h_enpr, ksatexm, bdens)** | **TOML missing** | 2 rows, one per soil-physical layer. Without these the Richards solver runs uninitialized. |
| swhyst, tau     | OK            | swhyst=0. |
| swmacro         | OK            | =0. |
| swsnow, swfrost | OK            | =0/0. |
| Numerical (dtmin, dtmax, gwlconv, critdevh1cp/h2cp, critdevponddt, maxit, maxbacktr, swkmean, swkimpl) | OK/partial | TOML has dtmin/dtmax. Initialize defaults match legacy's standard values. Lower priority. |

## .swp — Drainage / .dra

| Legacy key      | Status        | Notes |
|-----------------|---------------|-------|
| swdra, drfil    | OK            | =1, 'swap'. |
| dramet          | OK            | =3 (resistance multi-level). |
| swdivd          | OK            | =1. |
| **cofani[1..maho]** | **TOML/OK** | TOML authors `cofani = [1.0, 1.0]`. Wait — checked: salinitystress dra.toml does NOT author cofani. **TOML missing.** Without it cofani(:)=0 → drainage flux = 0. |
| swdislay        | OK            | =0. |
| nrlevs          | OK            | =2. |
| swintfl         | OK            | =0. |
| Per-level (DRAMET=3) | partial  | Both levels in TOML have drares/infres/zbotdr/L/swdtyp. **Each level's `swallo` is missing.** Legacy reads SWALLO1=3, SWALLO2=1. Without it Initialize defaults to 0 → invalid. |
| L1, L2 (m→cm)   | OK (FINALIZE) | adapter does ×100 conversion for DRAMET=3. |
| zbotdr          | OK            | -60.0 / -120.0 cm. |
| **DATOWL/LEVEL tables** | TOML missing | DRAMET=3 with SWDTYP=1 (drain tube, level 1) + SWDTYP=2 (open channel, level 2). For DRAMET=3 readswap.f90:rddra reads channel water level table for SWDTYP=2 only. Case 5: level 2 has DATOWL/LEVEL table 2012-01-01..2015-12-31 = -85.0. Legacy uses `level` to set the surface-water elevation; **without it level 2 drainage uses default -inf or 0 channel level → infiltration-from-drain disabled / wrong gradient**. Schema gap: per-level wls/level table. |

## .swp — Bottom boundary

| Legacy key      | Status        | Notes |
|-----------------|---------------|-------|
| swbbcfile       | OK            | =0; inline. |
| swbotb          | OK            | =3 (deep aquifer). |
| swbotb3resvert  | TOML missing  | =0 in case. Initialize default 0 → matches. N/A. |
| swbotb3impl     | TOML missing  | =0 in case. Initialize default 0 → matches. N/A. |
| shape           | OK            | 0.79; adapter copies. |
| hdrain          | OK            | -110.0. |
| rimlay          | OK            | 50.0. |
| **sw3**         | **N/A**       | =1 (sinus function); adapter copies aqave/aqamp/aqper/aqtmax. Matches. |
| aqave/aqamp/aqtmax/aqper | OK   | -120 / 20 / 120 / 365. |
| sw4             | OK            | =0; no extra qbot4 flux. |

## .swp — Heat

| Legacy key      | Status        | Notes |
|-----------------|---------------|-------|
| swhea           | OK            | =1. |
| swcalt          | OK            | =2. |
| **psilt**       | **TOML missing** | case authors 0.05 / 0.05 (2 layers). Validator fires: "psand/psilt/pclay/porg must all be set together". **Blocks startup.** |
| psand, pclay, porg | OK (TOML authored) | psand=[0.935, 0.935], pclay=[0.015, 0.015], porg=[0.019, 0.019]. |
| **tsoil_init**  | OK            | TOML authored 4 rows -10/-40/-70/-95 vs 15/12/10/9. **N/A: SWINCO=2, legacy only reads tsoil_init under SWINCO=1 or 2** (readswap.f90:tsoil block). |
| swtopbhea, swbotbhea | OK       | =1/1. |

## .swp — Solute

| Legacy key      | Status        | Notes |
|-----------------|---------------|-------|
| swsolu          | OK            | =1; adapter copies to global. |
| **cpre**        | **READER + ADAPTER missing** | legacy `rdsdor('cpre',...)` at readswap.f90:1107. Schema has no `cpre` slot. Case authors cpre=0.0 (default-ish), but the global stays at whatever Initialize sets. Verify: `cpre` is initialized to 0 in `variables.f90` → matches case default. **No-impact for case 5 (cpre=0).** Note: schema gap (low priority — Phase 4f-extend). |
| **cdrain**      | OK            | TOML authors 1.051; adapter copies. |
| swbotbc         | OK            | TOML authors 1; adapter copies. |
| **cseep**       | OK            | TOML authors 15.574; adapter copies. |
| **cml(macp) initial concentration** | **READER + ADAPTER missing** | legacy reads `zc(:)`/`cml(:)` arrays under `swinco != 3` (readswap.f90:1131-1135). Case 5 has SWINCO=3 (legacy reads from .ini), so the inline cml/zc table is **N/A in this case**. (TOML uses swinco=2 stand-in; under swinco=2 legacy WOULD read cml/zc, but they aren't authored in the .swp.) **Likely no-impact for case 5** unless the .ini-restart mismatch exposes a missing initial conc. |
| **ddif**        | **READER + ADAPTER missing** | legacy `rdsdor('ddif',...)`. Case authors ddif=0.0. Schema gap. **No-impact for case 5 (ddif=0 ≡ default).** |
| **tscf**        | OK            | TOML authors 0.5; adapter copies. |
| **ldis**        | **TOML/ADAPTER mismatch** | Schema has scalar `solute.ldis`. Legacy reads per-layer `ldis(maho)` — `rdfdor('ldis',...,ldis,maho,numlay)` (readswap.f90:1140). Case authors `LDIS = 5.0 5.0` (per-layer). Adapter currently DOES NOT copy `ldis` to global (line 700 says: "ldis(maho): config side has no per-layer schema yet; leave as-is"). **Default Initialize sets ldis=0 → dispersion off → solute concentrations propagate purely advectively → CONC profiles will be too sharp.** Need: (1) extend schema to allow per-layer `ldis(:)`, (2) author the array in TOML, (3) wire adapter to copy ldis(:). **Highest-impact gap for solute fixture.** |
| swsp            | OK            | =0; adapter doesn't need to wire frexp/cref/kf. |
| swdc            | OK            | =0; adapter doesn't need to wire decpot/gampar/rtheta/bexp/fdepth. |
| swbr            | **READER + ADAPTER missing** | legacy `rdsinr('swbr',...)`. Case authors swbr=0. Under swbr=0 legacy hardcodes daquif=100, poros=1, kfsat=0, decsat=0 (readswap.f90:1198-1201). Schema has no swbr slot. **Need adapter logic to set these globals to the swbr=0 defaults**, OR ensure the Initialize defaults match. Verify: `daquif/poros/kfsat/decsat` initialization in `variables.f90`. |

## .swp — Macropore (case has SWMACRO=0)

N/A — case excluded from macropore by SWMACRO=0.

## Highest-priority gaps (water-balance + solute-balance impact, blocking first)

1. **`[heat].psilt`** (TOML-missing) — currently aborts startup. Author 2 silt
   fractions (0.05, 0.05). Immediately unblocks. **HARD BLOCK.**
2. **`[soil.hydraulics]` MvG params + `bdens`** (TOML-missing) — without
   these the Richards solver runs uninitialized. **HARD BLOCK on water flow.**
3. **`[soil].sublay/isoillay/hsublay/ncomp`** (TOML-missing) — discretization
   table; without it the adapter cannot lay out the 2-physical-layer profile.
   **HARD BLOCK.**
4. **`[drainage].cofani`** (TOML-missing) — without it cofani(:)=0 → no
   drainage. **HARD BLOCK on drainage.**
5. **`[drainage.levels].swallo`** (TOML-missing) — per-level allowance;
   case has SWALLO1=3 (no infiltration), SWALLO2=1 (drain+inf). Without it
   defaults to 0 → invalid.
6. **`[general].inlist_csv`** override — case lists 8-col fixture columns
   (`cwso,cpwso,tredwet,treddry,tredsol,conc[-5.0,-25.0,-55.0]`); HACK default
   has 12 water-balance columns. **HARD BLOCK on regression compare.**
7. **`[soil].rsoil` = 150.0** (TOML-missing) — EPOT contributor.
8. **`[soil].rsro` + `rsroexp`** (TOML-missing) — runoff resistance.
9. **`[meteo.evaporation].cofredbo` + swredu=2 path** (TOML-missing /
   ADAPTER mismatch) — case authors swredu=2 with cofredbo=0.54. Adapter
   HACKs swredu=1; need swredu=2 with cofredbo=0.54.
10. **`[solute].ldis` per-layer + adapter wiring** (SCHEMA + TOML +
    ADAPTER missing) — case has 2-layer dispersion length 5.0/5.0. Adapter
    does not copy ldis(maho). **Highest-impact gap for solute fixture
    after the basic startup gaps clear.**
11. **`[irrigation].swirgfil + irgfil`** (READER + ADAPTER missing) — case
    has external swap.irg with ~370 daily gifts 2012-2015. Without wiring
    no irrigation events fire. **HARD BLOCK on solute mass balance** since
    the .irg authors irconc=0.401 mg/cm3 saline irrigation. Workaround:
    translate .irg to inline `[[irrigation.fixed_events]]` TOML.
12. **`[drainage.levels].level_table`** (SCHEMA + TOML missing) —
    DATOWL/LEVEL table for level 2 drives surface-water elevation in the
    drain-channel formulation. Lower priority (~–85 cm constant; legacy
    DEFAULT for missing table is unclear; verify post-fix).
13. **`[solute].swbr` defaults + per-layer cml** (low-impact; case has
    swbr=0 → defaults). Only worth wiring if other gaps don't close
    the diff.

## Action plan

1. Add `psilt` + complete the heat block (case has psand/pclay/porg
   authored already; just needs psilt) — unblocks startup.
2. Add `[soil.hydraulics]` MvG block + `[soil].sublay/.../ncomp` discretization.
3. Add `[drainage].cofani` and per-level `swallo`.
4. Add `rsoil`, `rsro`, `rsroexp` to TOML.
5. Add `inlist_csv` override + `cofredbo` + flip swredu.
6. Wire ldis per-layer (schema + reader + adapter + TOML).
7. Wire swirgfil/.irg (translate inline as fixed_events).
8. Verify regression converges. If diff persists, dig DATOWL/LEVEL.

## Not-yet-required schema gaps (Phase 4f-extend)

- `[solute].cpre`, `[solute].ddif`, `[solute].swbr`, `[solute].swsp` (and
  swsp=1 fields frexp/cref/kf), `[solute].cml`/`zc` initial-concentration
  table, decomposition fields (decpot/gampar/rtheta/bexp/fdepth).
  Case 5 has all of these at default-equivalent values (=0 / swsp=0 /
  swdc=0 / swbr=0), so adapter inaction matches case behaviour.
- `[soil].swinco=3` + `inifil` parsing for restart from `swap.ini`. Case 5
  uses SWINCO=3; the TOML stand-in is SWINCO=2 with gwli=-90, which
  may bias year-1 trajectories slightly.
