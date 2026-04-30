# Phase 4f — oxygenstress exact-parity audit

Phase 4f Task B4: walk every legacy `.swp` and `.dra` key for case 4
(`tests/swap-cases/4.oxygenstress/`) and classify against the TOML pipeline
(`src/config/*_config.f90`, `src/io/toml/read_*_toml.f90`,
`src/io/toml/config_to_variables.f90`).

Status legend (same as hupselbrook / grassgrowth audits):
- OK — TOML authors it, reader copies it, adapter writes the legacy global, no
  finalize mismatch.
- TOML — schema slot exists but the case TOML is missing the value.
- READER — schema field exists but no reader populates it.
- ADAPTER — schema + reader fine, but adapter never copies it into the global.
- FINALIZE — adapter copies the value but skips a unit-conversion / derived
  finalize step the legacy reader performs.
- N/A — switch is read but its value is irrelevant (gated branch not taken).

Baseline (commit a8cee6d): `pixi run -e test regression oxygenstress` aborts
with `psand/psilt/pclay/porg must all be set together` (heat validator fires
because the case TOML authors `psand/pclay/porg` but not `psilt`).

## .swp — General + Output

| Legacy key      | Status   | Notes |
|-----------------|----------|-------|
| project, paths  | OK       | adapter line 59-63 |
| swscre, swerror | OK       | |
| tstart / tend   | OK       | 1993-01-01 .. 2002-12-31 |
| nprintday       | OK       | |
| swmonth         | OK       | =0; PERIOD/SWRES/SWODAT path. |
| period, swres, swodat | OK | adapter copies. |
| swyrvar, datefix | TOML/N/A | swbal=0 → no .bal output → N/A. |
| outfil, swheader | OK (HACK)| hard-coded in adapter. |
| swwba/swend/swvap/swbal/swblc/swsba/swate/swbma/swdrf/swswb/swini/swinc/swcrp/swstr/swirg | OK (HACK) | RETIRED zero-forced. |
| swcsv           | OK (HACK)| hard-coded =1 in adapter. |
| **inlist_csv**  | **TOML missing** | case 4 .swp lists 6 columns (`pgrassdm,grassdm,pmowdm,mowdm,treddry,tredwet`). The TOML hasn't authored `inlist_csv` under `[general]`; adapter falls back to the HACK default `'pgrassdm,grassdm,pmowdm,mowdm'` (4 cols). The csv_writer / regression fixture reads `treddry/tredwet` for oxygen stress; without the override two columns will be missing. |
| swcsv_tz, inlist_csv_tz | OK (HACK) | hard-coded zero. |
| swafo, swaun    | OK (HACK)| zero-forced. |

## .swp — Meteorology

| Legacy key      | Status   | Notes |
|-----------------|----------|-------|
| metfil          | OK       | `260.met`. |
| lat, alt, altw  | OK       | |
| swetr           | OK       | =0 (PM). |
| angstroma/b     | OK       | |
| swdivide        | OK       | =1 (PM-direct). |
| swmetdetail     | OK       | =0. |
| swetsine        | OK       | =0. |
| swrain          | OK       | =2 → rainfall + duration from .met (no extra .swp keys). |
| rainfil, rainflux | N/A    | swrain != 3. |

## .swp — Crop section

| Legacy key      | Status   | Notes |
|-----------------|----------|-------|
| swcrop          | OK       | =1; adapter sets flCropReadFile/Open. |
| cropstart, cropend, cropfil, croptype | OK | 10 grass-detailed rotations 1993-2002. |
| rds (rdmax)     | OK (HACK)| hard-coded 200. |

> Note: `SwOxygen=2` is read by the legacy crop sub-reader (cropgrowth.f90 via
> readswap.f90:2233) from `grassd.crp`. The TOML adapter does not handle the
> .crp file; legacy readers consume the staged `grassd.crp` directly. The
> SwOxygen=2 branch checks `bdens(1) >= 100` (cropgrowth.f90:3668 / readswap
> :2225). **`bdens` flow is the parity-blocking dependency.**

## .swp — Irrigation

| Legacy key      | Status   | Notes |
|-----------------|----------|-------|
| swirfix         | OK       | =0. |

## .swp — Soil water

| Legacy key      | Status        | Notes |
|-----------------|---------------|-------|
| swinco          | OK            | =2. |
| gwli            | OK            | -15.0 cm. |
| swpondmx        | ADAPTER       | =0; Initialize default 0; N/A for case. |
| pondmx          | OK            | 0.2. |
| **rsro**        | **TOML missing** | legacy reads 0.5 d. Without it Initialize zeroes → instantaneous-runoff branch in `surfacewaterutils.f90`. Drainage / runoff drift source. (Same as hupselbrook / grassgrowth.) |
| **rsroexp**     | **TOML missing** | legacy reads 1.0. |
| swrunon         | TOML/ADAPTER  | =0; Initialize default 0; N/A. |
| cfevappond      | OK (HACK)     | adapter HACK 1.25; matches case. |
| swcfbs / cfbs   | OK            | swcfbs=0. |
| **rsoil**       | **TOML missing** | case 4 .swp authors **600.0 s/m** (matches grassgrowth). Schema slot already exists. EPOT direct contributor. |
| swredu          | OK (HACK)     | =1. Matches case. |
| **cofredbl**    | **TOML missing** | case authors 0.35 (= grassgrowth). Mapped to legacy `cofred` via `config%soil%evaporation%cofredbl`. |
| **rsigni**      | OK (HACK)     | adapter HACK 0.5; matches case. |
| **isublay/isoillay/hsublay/ncomp** | **TOML missing** | 8-row sub-layer table spanning 4 soil-physical layers (total = 5+20+15+35+15+110+100+100 = 400 cm, 5+8+3+7+2+11+5+4 = 45 compartments). |
| swsophy         | OK            | =0 → MvG. |
| **MvG params (ores, osat, alfa, npar, ksatfit, lexp, h_enpr, ksatexm, bdens)** | **TOML missing** | 4 rows, one per soil-physical layer. **bdens is load-bearing** for SwOxygen=2 (must be ≥ 100). Case authors `bdens = [557, 222, 194, 137]` — all ≥ 100. |
| swhyst, tau     | OK            | swhyst=0. |
| swmacro         | OK            | =0. |
| swsnow, swfrost | OK            | =0/0. |
| Numerical (dtmin, dtmax, gwlconv, critdevh1cp, critdevh2cp, critdevponddt, maxit, maxbacktr, swkmean, swkimpl) | TOML/partial | TOML has only dtmin/dtmax in `[simulation.numerical]`. Other numerical keys exist in schema but unauthored. Initialize defaults match legacy's standard values for most; TOML authoring still desirable for parity. |

## .swp — Drainage / .dra

| Legacy key      | Status        | Notes |
|-----------------|---------------|-------|
| swdra, drfil    | OK            | =1, 'swap'. |
| dramet          | OK            | =3 (resistance multi-level). |
| swdivd          | OK            | =1. |
| **cofani[1..maho]** | **TOML missing** | case authors `1.0 1.0 1.0 1.0` (4 layers). Without it cofani(:)=0 → drainage flux = 0. (Same gap pattern as grassgrowth.) |
| swdislay        | OK            | =0. |
| nrlevs          | OK            | =2. |
| swintfl         | OK            | =0. |
| Per-level (DRAMET=3) | partial | Both levels in TOML have drares/infres/zbotdr/L/swdtyp. **Each level's `swallo` is missing.** Legacy reads SWALLO1=1 and SWALLO2=3. Without it Initialize defaults to 0 → invalid; adapter has `swliminf=1` HACK for DRAMET=3 but per-level swallo controls drain/infiltration allowance. |
| L1, L2 (m→cm)   | OK (FINALIZE) | adapter does ×100 conversion for DRAMET=3 (commit 1b00cdd era). |
| zbotdr          | OK            | -140 / -20. |

## .swp — Bottom boundary

| Legacy key      | Status        | Notes |
|-----------------|---------------|-------|
| swbbcfile       | OK            | =0; inline. |
| swbotb          | OK            | =3 (Cauchy / regional aquifer). |
| swbotb3resvert  | TOML missing  | =0 in case. Initialize default 0; N/A water-balance. (Schema slot does not exist; gated by SWBOTB=3.) |
| swbotb3impl     | TOML missing  | =1 in case. Initialize default 0 → flux solved explicitly instead of implicitly. **Affects bottom-flux computation.** Schema slot does not exist. |
| shape           | OK            | adapter copies. |
| hdrain          | OK            | -25.0. |
| rimlay          | OK            | 500.0. |
| sw3             | TOML missing  | =2 (table). Schema does not select between sinus / table; the adapter currently always copies aqave/aqamp/aqper/aqtmax (sinus), but this case uses the date table. |
| **HAQUIF table (date3, haquif)** | **MISSING** | case 4 inlines 80-row monthly aquifer-head table 1993-01-01 .. 2002-12-31 with values −36.66 / −57.47 / −60.11 / −46.70. **Drives the regional-aquifer head; controls bottom flux and thus GWL / DRAINAGE.** Highest-impact gap. Schema currently has only `cofqha_table` (SWBOTB=4); no `haquif_table` slot for SWBOTB=3 sw3=2. |
| sw4             | OK            | =0; no extra qbot4 flux. |

## .swp — Heat

| Legacy key      | Status   | Notes |
|-----------------|----------|-------|
| swhea           | OK       | =1. |
| swcalt          | OK       | =2. |
| **psilt**       | **TOML missing** | case authors 0.108 / 0.138 / 0.528 / 0.373 (4 layers). Validator fires: "psand/psilt/pclay/porg must all be set together". **Blocks startup.** |
| psand, pclay, porg | OK   | already authored. |
| tsoil_init      | OK       | already authored. |
| swtopbhea, swbotbhea | OK | =1/1. |

## .swp — Solute

| Legacy key      | Status   | Notes |
|-----------------|----------|-------|
| swsolu          | OK       | =0. |

## Highest-priority gaps (water-balance impact, blocking first)

1. **`[heat].psilt`** (TOML-missing) — currently aborts startup. Author 4 silt
   fractions, immediately unblocks.
2. **`[soil.hydraulics]` MvG params + `bdens`** (TOML-missing) — without these
   the soil hydraulics / oxygen-stress paths are uninitialized. `bdens` ≥ 100
   is required for SwOxygen=2 (cropgrowth.f90:3668). Affects every flux.
3. **`[soil].sublay/isoillay/hsublay/ncomp`** (TOML-missing) — discretization
   table; without it the adapter cannot lay out the 4-physical-layer profile
   correctly.
4. **`[drainage].cofani`** (TOML-missing) — without it cofani(:)=0 → no
   drainage. (Already known pattern.)
5. **`[drainage.levels].swallo`** (TOML-missing) — per-level allowance.
6. **`[bottom_boundary].haquif_table`** (SCHEMA + TOML missing) — the SWBOTB=3
   sw3=2 case uses a date-keyed aquifer-head table that drives the regional
   aquifer head. **Largest-impact gap once startup unblocks.** Requires schema
   slot + reader + adapter wiring (mirror `cofqha_table` for SWBOTB=4).
7. **`[soil].rsoil` = 600.0** (TOML-missing) — same as grassgrowth.
8. **`[soil].rsro` + `rsroexp`** (TOML-missing) — same as hupselbrook.
9. **`[soil].cofredbl`** (TOML-missing) — soil evaporation coefficient.
10. **`[general].inlist_csv`** override — case lists 6 columns including
    `treddry,tredwet` for oxygen-stress reporting; HACK default has only 4.
11. **`swbotb3impl`** (SCHEMA + TOML missing) — case has =1 (implicit). Adapter
    leaves the global at 0 (explicit). Likely modest-impact.
12. **Per-soil-numerical** (gwlconv, critdevh1cp/h2cp, critdevponddt, maxit,
    maxbacktr, swkmean, swkimpl) — already in schema; TOML omits. Initialize
    defaults are close to legacy defaults; lower priority.

## Action plan

1. Add `psilt` to oxygenstress TOML — unblocks startup.
2. Add `[soil.hydraulics]` MvG block + `[soil].sublay/.../ncomp` discretization.
3. Add `[drainage].cofani` and per-level `swallo`.
4. Add `[bottom_boundary].haquif_table` schema slot (and reader + adapter) +
   author the 80-row table from the .swp.
5. Add `rsoil`, `rsro`, `rsroexp`, `cofredbl` to TOML.
6. Add `inlist_csv` override.
7. Add `swbotb3impl` schema slot + adapter copy.
8. Verify regression converges.
