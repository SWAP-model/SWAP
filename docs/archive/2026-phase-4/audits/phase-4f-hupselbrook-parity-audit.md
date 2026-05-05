# Phase 4f — hupselbrook exact-parity audit

Phase 4f's strangler swap (commit 744cbfd) routes hupselbrook end-to-end through
the new TOML pipeline. The case completes but the regression diff still shows
~0.3-0.5 cm/yr drift on EACT/EPOT/DRAINAGE and ~2 cm on GWL, which is far above
the 1e-2 cm tolerance.

This audit walks every legacy `rd*` call in `src/io/readswap.f90` and matches it
against the TOML schema (`src/config/*_config.f90`), the per-section reader
(`src/io/toml/read_*_toml.f90`), and the strangler adapter
(`src/io/toml/config_to_variables.f90`). The first table is restricted to keys
the hupselbrook `.swp` actually authors (not every legacy switch).

Status legend:
- OK — TOML authors it, reader copies it, adapter writes the legacy global, no
  finalize mismatch.
- TOML — value missing from `tests/swap-cases/toml/1.hupselbrook/swap.toml`.
- READER — schema field exists but no reader populates it.
- ADAPTER — schema + reader fine, but adapter never copies it into the global.
- FINALIZE — adapter copies the value but skips a unit-conversion or derived
  finalize step the legacy reader performs.
- N/A — switch is read but its value is irrelevant (gated branch not taken).

## hupselbrook .swp — General + Output

| Legacy key      | Status   | Notes |
|-----------------|----------|-------|
| project, paths  | OK       | adapter line 59-63 |
| swscre          | OK       | |
| swerror         | TOML     | TOML authors `swerror=0` but no schema slot. Initialize sets 0; matches. N/A. |
| tstart / tend   | OK       | |
| nprintday       | OK       | |
| swyrvar, datefix, outdat | TOML | output-cycle inputs; gated by SWBAL/SWBLC etc. swbal=0 in HACK so N/A. |
| swmonth, period, swres, swodat, outdatint | partial | period/swres/swodat in adapter; swmonth + outdat[int] not. SWMONTH=1 means monthly output ticks; csv driver does not check it, so N/A for water balance. |
| outfil          | OK (HACK)| hard-coded 'result' in adapter. |
| swheader        | OK (HACK)| zero-forced. |
| swwba/swend/swvap/swbal/swblc/swsba/swate/swbma/swdrf/swswb/swini/swinc/swcrp/swstr/swirg | OK (HACK) | RETIRED switches zero-forced (ADR 0009). |
| swcsv, inlist_csv | OK (HACK) | hard-coded in adapter. |
| swcsv_tz, inlist_csv_tz | OK (HACK) | hard-coded zero. |
| swafo, swaun, critdevmasbal, swdiscrvert, numnodnew, dznew | partial | swafo/swaun zero-forced; discretization in adapter; critdevmasbal not (output deviation gate, no water balance impact). |

## hupselbrook .swp — Meteorology

| Legacy key      | Status   | Notes |
|-----------------|----------|-------|
| metfil          | OK       | with .met-suffix derivation of swMetFilAll. |
| lat, alt, altw  | OK       | |
| swetr           | OK       | |
| angstroma/b     | OK       | |
| swdivide        | OK       | hupselbrook=1 → PM-direct. |
| swmetdetail, nmetdetail | OK | swmetdetail=0 in case 1. |
| swetsine        | OK       | |
| swrain, rainfil, rainflux | OK | swrain=0, rainfil unused. |

## hupselbrook .swp — Crop section

| Legacy key      | Status   | Notes |
|-----------------|----------|-------|
| swcrop          | OK       | adapter line 635 + sets flCropReadFile/Open. |
| cropstart, cropend, cropfil, croptype | OK | adapter writes per-rotation arrays. |
| rds (rdmax)     | OK (HACK)| hard-coded 200 cm in adapter. |

## hupselbrook .swp — Irrigation

| Legacy key      | Status   | Notes |
|-----------------|----------|-------|
| swirfix, swirgfil | OK     | swirfix=1, swirgfil=0 → fixed_events table consumed. |
| irdate, irdepth, irconc, irtype | OK | adapter at line 548 with `/10` mm→cm finalize. |

## hupselbrook .swp — Soil water

| Legacy key      | Status        | Notes |
|-----------------|---------------|-------|
| swinco          | OK            | =2 → hydrostatic eq. with gwli. |
| gwli            | OK            | -75 cm. |
| **rsro**        | **MISSING**   | legacy reads 0.5 d, Initialize zeroes it. With rsro=0, surfacewater.f90:205 hits the "instantaneous drainage" branch; legacy uses the dt/rsro*(pond-pondmx)^rsroexp formula. **DRAINAGE / runoff drift source.** |
| **rsroexp**     | **MISSING**   | legacy reads 1.0 (linear). |
| swrunon         | TOML/ADAPTER  | Legacy reads at line 1011. swrunon=0 in case 1; Initialize defaults to 0 → no observable drift, but adapter never sets it. Add for parity. |
| swpondmx        | ADAPTER       | swpondmx=0 → constant pondmx; Initialize defaults to 0; no impact for case 1. |
| pondmx          | OK            | 0.2 cm. |
| cfevappond      | OK (HACK)     | 1.25. |
| swcfbs / cfbs   | OK            | swcfbs=0; cfbs unused. |
| **rsoil**       | **TOML missing** | Legacy: 30.0 s/m (read only when swdivide=1). TOML authors no value → defaults 0.0 → PM-direct uses rss=0 → **EPOT direct contributor.** |
| swredu          | OK (HACK)     | =1 (Black). |
| rsigni, cofred  | OK            | rsigni HACKed; cofred mapped from cofredbl. |
| isoillay/isublay/hsublay/ncomp | OK | discretization. |
| swsophy / MvG params (ores, osat, alfa, npar, ksatfit, lexp, alfaw, h_enpr, ksatexm, bdens) | OK | paramvg layout matches readswap.f90:786. |
| swhyst, tau     | OK            | swhyst=0 → tau N/A. |
| swmacro         | OK            | =0; case 1 macropore-free. |
| swsnow          | OK            | =0. |
| swfrost         | OK            | =0. |
| Numerical (dtmin, dtmax, gwlconv, critdevh1cp, critdevh2cp, critdevponddt, maxit, maxbacktr, swkmean, swkimpl) | OK | all in [simulation.numerical]. |

## hupselbrook .swp — Drainage / .dra

| Legacy key      | Status        | Notes |
|-----------------|---------------|-------|
| swdra, drfil    | OK            | swdra=1, drfil='swap'. |
| dramet          | OK            | =2 (Hooghoudt). |
| swdivd          | OK            | =1. |
| cofani[1..maho] | OK            | adapter line 240. |
| swdislay        | OK            | =0. |
| lm2 (m→cm)      | OK            | adapter line 220 does 100×lm. |
| shape           | OK            | drainage block. |
| wetper, zbotdr (basic) | OK     | level-1 entry of the per-level array. |
| entres          | OK            | |
| ipos            | OK            | =2. |
| basegw, khtop   | OK            | |
| khbot/zintf/kvtop/kvbot/geofac | N/A | gated by ipos; ipos=2 so skipped. |
| nrlevs etc (DRAMET=3) | N/A | DRAMET=2 selected. |
| swintfl, cofintflb, expintflb | N/A | DRAMET=3 only. |

## hupselbrook .swp — Bottom boundary

| Legacy key      | Status   | Notes |
|-----------------|----------|-------|
| swbbcfile       | TOML/ADAPTER | =0 in case 1; Initialize default 0 → no impact. |
| swbotb          | OK       | =6 (zero flux). |

(All the SWBOTB=1..5 sub-blocks are inactive; no parity issue for case 1.)

## hupselbrook .swp — Heat

| Legacy key      | Status   | Notes |
|-----------------|----------|-------|
| swhea           | OK       | =1. |
| swcalt          | OK       | =2 (numerical). |
| **psilt**       | **MISSING** | legacy reads 0.15. config has psand/pclay/porg only. Affects soil thermal conductivity (heat module); indirect impact on water flow via temperature-dependent K. **Likely cause of small EACT drift.** |
| psand, pclay, orgmat (porg) | OK | adapter copies. |
| tsoil_init      | OK       | adapter copies col 2 element-by-element. |
| swtopbhea, swbotbhea | OK   | both =1. |

## hupselbrook .swp — Solute

| Legacy key      | Status   | Notes |
|-----------------|----------|-------|
| swsolu          | OK       | =1. |
| **cpre**        | **MISSING** | legacy reads 0.0; Initialize defaults 0; no drift but unauthored. |
| cdrain          | OK       | 0.1 mg/cm3. |
| swbotbc         | OK       | =0. |
| cseep           | OK       | 0.1. |
| **cml** initial | **MISSING** | legacy interpolates the (zc, cml) table to per-compartment initial concentrations. case 1 specifies cml=0 everywhere; no observed drift, but unauthored. |
| **ddif**        | **MISSING** | molecular diffusion coefficient; case 1 = 0.0 = Initialize default; no drift. |
| tscf            | OK       | =0.0. |
| ldis (per-layer) | partial | adapter has `ldis(1) = config%solute%ldis` (scalar) — case 1 lists per-layer 5.0 5.0; scalar 5.0 broadcast OK. |
| **swsp**, frexp, cref, kf | **MISSING** | swsp=0 in case 1 → adsorption disabled; no drift. |
| swdc            | OK       | =0. |
| swbr            | TOML/ADAPTER | swbr=0; Initialize default 0; N/A. |

## Highest-priority gaps for hupselbrook exact-match

Ranked by expected impact on the diff table (largest signal first):

1. **`rsoil`** — TOML `[soil]` doesn't author `rsoil`. Schema default 0.0; legacy
   reads 30.0 s/m. With swdivide=1, `rsoil` enters the PM-direct soil-resistance
   `rss = rsoil` path in `et.f90`. Lower rsoil → higher Edirect → higher peva →
   higher EPOT. **Highest expected impact on EPOT (~1 cm/yr).**

2. **`rsro` + `rsroexp`** — Surface-runoff drainage resistance and exponent.
   Legacy reads 0.5 d and 1.0. Initialize zeroes both. With rsro=0, surface
   runoff is instantaneous (`rsro < 1e-3` branch in `surfacewaterutils.f90`)
   instead of resistance-controlled. Affects DRAINAGE and DSTOR.

3. **`psilt`** — Heat-module silt fraction; legacy reads 0.15 alongside psand/
   pclay/orgmat. Currently no schema slot. Affects soil thermal conductivity
   (when swhea=1 + swcalt=2), with feedback through temperature-dependent K.
   Smaller (~0.1-0.3 cm/yr) but real for case 1.

4. **`swrunon`** — Initialize default already 0; no observable drift. Worth
   wiring for parity.

5. **`cpre` / `ddif` / `cml(:)` initial** — solute side, all =0 in case 1; no
   water balance signal but missing for completeness.

The bottom-boundary/heat/solute remainders are gated by switches that the case
sets to inert values; they'd matter once cases 2..6 enter the loop.

## Action plan (one commit per gap, smallest commit first)

1. Add `rsoil` to `tests/swap-cases/toml/1.hupselbrook/swap.toml` as
   `[soil].rsoil = 30.0` (schema slot already exists). Submodule first.
2. Add `[soil.runoff]` `rsro` + `rsroexp` to schema + reader + adapter +
   case 1 TOML.
3. Add `psilt` to `[heat]` array list (schema + reader + adapter + case 1
   TOML).
4. Add `swrunon` to schema + adapter (no TOML change for case 1).

After each change: rebuild, rerun hupselbrook regression, record delta.
