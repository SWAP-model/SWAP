# Phase 4f — surfacewater exact-parity audit

Phase 4f Task B6: walk every legacy `.swp` and `.dra` key for case 6
(`tests/swap-cases/6.surfacewater/`) and classify against the TOML pipeline
(`src/config/*_config.f90`, `src/io/toml/read_*_toml.f90`,
`src/io/toml/config_to_variables.f90`).

Status legend (mirrors hupselbrook / grassgrowth / oxygenstress / salinitystress
audits):
- OK — TOML authors it, reader copies it, adapter writes the legacy global,
  no finalize mismatch.
- TOML — schema slot exists but the case TOML is missing the value.
- READER — schema field exists but no reader populates it.
- ADAPTER — schema + reader fine, but adapter never copies it into the
  legacy global.
- FINALIZE — adapter copies the value but skips a unit-conversion / derived
  finalize step the legacy reader performs.
- N/A — switch is read but its value is irrelevant (gated branch not taken).

Baseline (commit a8cee6d): `pixi run -e test regression surfacewater` aborts
with

    Routine RDDATA, called by RDATIM, attempts to read ... './SWAP.BBC'
    ... ERROR in RDDATA: Variable name not in data file (date1)

The TOML hard-codes `bottom_boundary.swbotb=1` and `bbcfil="swap"`. The case's
on-disk `swap.bbc` actually authors `SWBOTB=3` with a `DATE3`/`HAQUIF` aquifer-head
table, so the SWBOTB=1 adapter branch (line 620 of `config_to_variables.f90`)
fails immediately on the missing `date1` key.

## .swp — General + Output

| Legacy key      | Status        | Notes |
|-----------------|---------------|-------|
| project, paths  | OK            | adapter line 59-63 (`wildenborch`). |
| swscre, swerror | OK            | both 0. |
| tstart / tend   | OK            | 1997-01-01 .. 1999-12-31. |
| nprintday       | OK            | =1. |
| swmonth         | OK            | =0; PERIOD/SWRES/SWODAT path. |
| period, swres, swodat | OK      | period=1, swres=0, swodat=0. |
| swyrvar, datefix | OK/N/A       | swbal=0 → no .bal output → N/A. |
| outfil          | OK (HACK)     | hard-coded `'result'` in adapter. |
| swheader, swwba/swend/swvap/swbal/swblc/swsba/swate/swbma/swdrf/swswb/swini/swinc/swcrp/swstr/swirg | OK (HACK) | RETIRED zero-forced. |
| swcsv           | OK (HACK)     | hard-coded =1 in adapter. |
| **inlist_csv**  | OK            | TOML authors `inlist_csv = 'gwl,pond'` matching case .swp. fixture asserts on **GWL/POND**. |
| swcsv_tz        | OK (HACK)     | zero-forced. |
| swafo, swaun    | OK (HACK)     | zero-forced. |

## .swp — Meteorology

| Legacy key      | Status   | Notes |
|-----------------|----------|-------|
| metfil          | OK       | `290.met`. |
| lat, alt, altw  | OK       | 52.274 / 34.8 / 10.0. |
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
| cropstart, cropend, cropfil, croptype | OK | 3 `grass` rotations 1997..1999 type=1 (cropfixed). |
| **rds (rdmax)** | OK (HACK)     | Adapter HACKs rdmax=200. Case 6 .swp authors **60.0 cm**. The HACK 200 is permissive; rd capped by `min(afgen(rdtb,...), rdm)`, and case 6 grass.crp's rdtb max is 50 cm. **No water-balance impact.** |

## .swp — Irrigation

| Legacy key      | Status | Notes |
|-----------------|--------|-------|
| swirfix         | OK     | =0; no fixed events. |

## .swp — Soil water

| Legacy key      | Status        | Notes |
|-----------------|---------------|-------|
| swinco          | OK            | =2 (hydrostatic equilibrium with gwli). |
| gwli            | OK            | =-68.0 cm. |
| swpondmx        | OK            | =0; default. |
| pondmx          | OK            | 0.2. |
| **rsro**        | **TOML missing** | legacy reads 0.5 d. Schema slot exists; case TOML doesn't author it. Initialize zeroes → instantaneous-runoff branch. |
| **rsroexp**     | **TOML missing** | legacy reads 1.0. Same gap. |
| swrunon         | OK            | =0; flrunon=.false. |
| **cfevappond**  | OK (HACK)     | adapter HACKs 1.25; matches case. |
| swcfbs / cfbs   | OK            | swcfbs=0. |
| **rsoil**       | **TOML missing** | case authors **600.0 s/m**. Without it default 0 → wrong EPOT. |
| swredu          | OK (HACK)     | adapter HACKs swredu=1; case 6 .swp also authors swredu=1. Match. |
| **cofredbl**    | **TOML missing** | case authors 0.35; under swredu=1 the legacy reads it. Adapter has slot at `meteo.evaporation.cofredbl`. Initialize default may differ. |
| **rsigni**      | OK (HACK)     | adapter HACKs 0.5; matches case. |
| **isublay/isoillay/hsublay/ncomp** | **TOML missing** | 8-row sub-layer table spanning 3 soil-physical layers (10+15+35+20+80+120+90+35 = 405 cm, 10+3+7+2+8+6+3+1 = 40 compartments). |
| swsophy         | OK            | =0 (analytical MvG). |
| **MvG params (ores, osat, alfa, npar, ksatfit, lexp, h_enpr, ksatexm)** | **TOML missing** | 3 rows, one per soil-physical layer. **HARD BLOCK on water flow.** |
| swhyst          | OK            | =0. |
| swmacro         | OK            | =0. |
| swsnow, swfrost | OK            | =0/0. |
| Numerical (dtmin, dtmax, gwlconv, critdevh1cp/h2cp, critdevponddt, maxit, maxbacktr, swkmean, swkimpl) | partial | TOML has dtmin/dtmax. Other knobs use Initialize defaults; .swp values match defaults except gwlconv=100 / critdevh1cp=0.01 / critdevh2cp=0.1 / critdevponddt=0.0001. **TOML missing**, but Initialize defaults likely match. |

## .swp — Drainage / .dra (SWDRA=2 — extended drainage with surface water)

The case uses SWDRA=2, so legacy reads `.dra` via `rddre()` in
`src/drainage/surfacewater.f90:56` — NOT through `readswap.f90`. The
strangler adapter only needs to set `swdra=2` so `flSurfaceWater=.true.` is
flipped in `timecontrol.f90:114`. Once that's set, `rddre()` reads `swap.dra`
directly from disk via the legacy `ttutil` reader — the on-disk `swap.dra`
file is the source of truth.

| Legacy key (.dra) | Status | Notes |
|-------------------|--------|-------|
| swdivd          | OK     | =1 — case authors. Read by rddre. |
| **cofani[1..maho]** | OK (legacy) | =[10.0, 10.0, 10.0]. Read by rddre directly. |
| swdislay        | OK     | =0. Read by rddre. |
| altcu           | OK     | =0.0. Read by rddre directly. |
| nrsrf           | OK     | =2 levels. Read by rddre. |
| Per-level table (LEV/SWDTYP/L/ZBOTDRE/GWLINF/RDRAIN/RINFI/RENTRY/REXIT/WIDTHR/TALUDR) | OK (legacy) | Read by rddre directly from .dra. |
| swnrsrf, swtopnrsrf | OK (legacy) | =0/0. Read by rddre. |
| swsrf           | OK (legacy) | =2 (surface-water without separate primary). Read by rddre. |
| swsec           | OK (legacy) | =2 (simulated). Read by rddre. |
| wlact, osswlm   | OK (legacy) | =-77.0 / 2.5. Read by rddre. |
| nmper           | OK (legacy) | =28. Read by rddre. |
| imper_4b/impend/swman/wscap/wldip/intwl | OK (legacy) | Per-period management table. Read by rddre. |
| swqhr           | OK (legacy) | =1 (exponential). Read by rddre. |
| sofcu           | OK (legacy) | =100.0 ha. Read by rddre. |
| imper_4c/hbweir/alphaw/betaw | OK (legacy) | Per-period weir table. Read by rddre. |

**NOTE:** The TOML authors all the surface-water keys redundantly in
`[surface_water]` and `[drainage]` blocks, and the strangler adapter at
lines 798-846 of `config_to_variables.f90` copies them into module globals.
However, `rddre()` re-reads them from `swap.dra` and overwrites the adapter's
values. **The TOML surface-water block is currently inert for surfacewater
case execution** — the on-disk `swap.dra` is what matters. The TOML/schema
work in Phase 4f-prep was preparation for an eventual rddre-replacement.
**For Phase 4f Task B6, the only requirement is that swdra=2 lands in the
global and the on-disk swap.dra remains intact.**

## .swp — Bottom boundary (HARD BLOCK)

| Legacy key      | Status        | Notes |
|-----------------|---------------|-------|
| swbbcfile       | OK            | =1; external file `swap.bbc`. |
| **bbcfil**      | **TOML wrong** | TOML has `bbcfil="swap"`. Adapter's swbotb=1 branch reads .bbc looking for `date1` key but the file is actually SWBOTB=3 with `DATE3`/`HAQUIF` table. **HARD BLOCK at startup.** |
| **swbotb (in .bbc)** | **TOML wrong** | Case TOML `swbotb=1`. Real value: =3 (deep aquifer). |
| swbotb3resvert  | TOML missing  | =0 in case. Initialize default 0 → matches. N/A. |
| swbotb3impl     | TOML missing  | =0 in case. Default 0 matches. |
| shape           | TOML missing  | =1.0. (TOML has nothing because swbotb=1 in current TOML.) |
| hdrain          | TOML missing  | =-400.0. |
| rimlay          | TOML missing  | =10.0. |
| **sw3 + DATE3/HAQUIF table** | TOML missing | =2 (table). 73 (date, haquif) rows from 1996-12-28..2000-01-14. **HARD BLOCK on bottom flux.** |
| sw4             | OK (default)  | =0; no extra qbot4 flux. |
| aqave/aqamp/aqtmax/aqper | N/A  | sw3=2, not sinus. |

The simplest fix: rewrite the TOML's `[bottom_boundary]` block to mirror
oxygenstress (case 4) — `swbotb=3` with inline `haquif_table`, populated
from the on-disk `swap.bbc`. Drop `bbcfil` since we're authoring inline.

## .swp — Heat

| Legacy key      | Status | Notes |
|-----------------|--------|-------|
| swhea           | OK     | =0. No further heat fields needed. |

## .swp — Solute

| Legacy key      | Status | Notes |
|-----------------|--------|-------|
| swsolu          | OK     | =0. No further solute fields needed. |

## Highest-priority gaps (water-balance impact, blocking first)

1. **`[bottom_boundary]` swbotb=3 + haquif_table** (TOML wrong) — currently
   aborts at startup. Rewrite the block to mirror oxygenstress: `swbotb=3`,
   `shape=1.0`, `hdrain=-400.0`, `rimlay=10.0`, and an inline `haquif_table`
   with the 73 (date, haquif_cm) rows from `swap.bbc`. Drop `bbcfil`.
   **HARD BLOCK.**
2. **`[soil.hydraulics]` MvG params** (TOML missing) — without these the
   Richards solver runs uninitialized. **HARD BLOCK on water flow.**
3. **`[soil].sublay/isoillay/hsublay/ncomp`** (TOML missing) — discretization
   table; without it the adapter cannot lay out the 3-physical-layer profile.
   **HARD BLOCK.**
4. **`[soil].rsoil`** (TOML missing) — case authors 600.0 s/m; without it
   default 0 → wrong EPOT.
5. **`[soil].rsro` + `rsroexp`** (TOML missing) — runoff resistance (0.5/1.0).
6. **`[meteo.evaporation].cofredbl`** (TOML missing) — case authors 0.35;
   without it default may diverge.

Note that for SWDRA=2 (surface water) the .dra read path is via `rddre()` in
`src/drainage/surfacewater.f90`, NOT readswap. So as long as `swdra=2` is in
the global and the on-disk `swap.dra` is unchanged, the legacy reader handles
the entire extended drainage / surface-water block. The TOML's redundant
`[surface_water]` block is currently inert for case 6.

## Action plan

1. Fix `[bottom_boundary]` to use swbotb=3 + haquif_table from `swap.bbc`
   (mirror oxygenstress). Drop the bbcfil-pointing-at-real-file pattern.
2. Author `[soil.hydraulics]` MvG block + `[soil]` discretization.
3. Author `rsoil`, `rsro`, `rsroexp`, `cofredbl` in TOML.
4. Verify regression converges; if it doesn't, dig the surface-water block
   (rddre re-reads from .dra so it should already be in sync with the legacy
   numbers).

## Not-yet-required schema gaps (Phase 4f-extend)

- `[bottom_boundary]` external `.bbc` file reader for arbitrary swbotb (the
  current adapter only handles swbotb=1 + bbcfil; swbotb=3 + bbcfil would
  require new code). Workaround for case 6: inline the haquif_table.
- `[surface_water]` adapter rddre-replacement — the TOML schema slots all
  exist (Phase 4f-prep Task D2) but the adapter assignments are currently
  shadowed by `rddre()` re-reading from `.dra`. Phase 4f-extend can drop the
  rddre call once the adapter wires altcu / nrpri / nrsec / wlstar / wls1 /
  swstini / sttab.
