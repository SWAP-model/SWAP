# Reconstructing the lost base regression cases (2026-06-11)

The swap-cases `toml/` tree was lost from the remote (submodule pinned a commit
upstream no longer has). The base cases are being rebuilt as self-contained LOCAL
cases under `tests/regression/cases/<name>/{legacy,toml}/` by converting the
legacy ASCII (`tests/swap-cases/<N>.<case>/`, still present) to modern TOML.

## Status

| case | bottom | crop | extras | status |
|------|--------|------|--------|--------|
| hupselbrook | swbotb=6 | maize/potato/grass | solute, irrig | ✅ byte-identical |
| grassgrowth | swbotb=1 (gwl_file) | grass | dramet=3 drainage | ✅ byte-identical |
| oxygenstress | swbotb=3 Cauchy (haquif_file) | grass + Bartholomeus O2 | dramet=3 2-level | ✅ byte-identical |
| surfacewater | swbotb=3 Cauchy | type-1 grass | **swdra=2 surface water** (28 periods) | ✅ byte-identical |
| salinitystress | swbotb=3 | wofost potato + salinity | swsolu=1, swinco=3 warm-restart, swap.irg | ⚠️ reconstructed, ~3% solute residual (xfail) |
| macroporeflow | — | — | macropore | retired (ADR 0040, physics deleted) |

## Conversion recipe (proven on grassgrowth + oxygenstress)

1. **Dir:** `tests/regression/cases/<name>/{legacy,toml}/`. Copy the legacy
   `.swp.template`, `.crp`, `.dra`, `.bbc`, `.met` into `legacy/`.
2. **Meteo:** `python tests/regression/met_to_csv.py <case>/260.met toml/260.csv`.
3. **swap.toml** — map legacy `.swp` sections (use an existing case as template):
   - `[general]` project/swscre/swerror; `[simulation]` start/end/nprintday;
     `[simulation.output]` swmonth/period/...; `[simulation.numerical]` dt/dtmin/dtmax.
   - `[meteorology]` file/lat/alt/altw; `.evapotranspiration` swetr/swdivide/angstrom;
     `.rain` swrain/swetsine; `.interception` swinter=0; `.evaporation`
     swcfbs/swredu/cofredbl/rsigni/cfevappond.
   - `[soil]` swsophy/swhyst/swinco/swmacro/rsoil/rsro/rsroexp + discretization
     arrays `sublay/isoillay/hsublay/ncomp` (legacy ISUBLAY/ISOILLAY/HSUBLAY/NCOMP;
     hcomp = hsublay/ncomp is implicit).
   - `[soil.initial]` gwli/pondini/pondmx (+ `h_file` for swinco=3).
   - `[soil.hydraulics]` per-layer arrays ores/osat/alfa/npar/ksatfit/lexp/alfaw/
     h_enpr/ksatexm/bdens. **alfaw is required even when swhyst=0** — mirror alfa.
   - `[crop]` swcrop=1 + one `[[crop.rotation]]` per legacy CROPSTART row
     (yearly entries matter for byte-identity).
   - `[bottom_boundary]` swbotb + per-mode params:
     - swbotb=1: `gwl_file="gwlevel.csv"` (CSV header `date,gwl`).
     - swbotb=3: swbotb3resvert/swbotb3impl/shape/hdrain/rimlay/sw3 +
       `haquif_file="haquif.csv"` (header `date,haquif`) when sw3=2.
   - `[heat]` (swhea=1) swcalt/swtopbhea/swbotbhea + psand/psilt/pclay/porg +
     `tsoil_init` rows. Omit when swhea=0.
   - `[solute]` swsolu (+ cdrain/cseep/... when swsolu=1).
   - `[output.csv] inlist = "..."` — MUST set to the legacy INLIST_CSV columns,
     else the modern emits the default water-balance set and the asserted columns
     come back empty.
4. **swap.dra.toml** — `[drainage]` swdra/dramet/swdivd/swdislay/nrlevs/cofani.
   - dramet=2 (Hooghoudt): `[drainage.basic]` lm/shape/wetper/zbotdr/entres/ipos/...
   - dramet=3 (resistance): one `[[drainage.levels]]` per legacy level with
     swdtyp/swallo/zbotdr/drares/infres/L (+ `owltab_file="owltabN.csv"` header
     `date,level` for open channels, swdtyp=2).
   - **swdra=2 (surface water): NOT yet mapped** — needs the surface-water
     reservoir + management-period schema (NMPER periods, SWSRF/SWSEC, WLACT...).
5. **Crop `.crp.toml`** — adapt the matching existing crop toml (grassd/potatod)
   by applying the legacy `.crp`'s differing scalars + tables. Gotchas:
   - swrd=2 needs rdi/rri/rdc/swdmi2rd; swrd=3 needs rlwtb/wrtmax.
   - grass swharv=2 (fixed dates): legacy `dateharvest` -> `mowing_dates` as
     **day-of-year floats** + nmow; the runtime advances the year when DOY resets.
   - swoxygen=2: add `[oxygen_stress.bartholomeus]` AFTER the oxygen scalars
     (before `[drought_stress]`) so swwrtnonox/aeratecrit stay in `[oxygen_stress]`.
   - swhyst=1 needs `tau` (legacy TAU); modern defaults tau=0 (degenerate).
6. **Fixture + verify:** `python regen_reference.py <name>` (swap420gf on the
   legacy dir), then run the case; iterate until byte-identical.

## Remaining work (all 4 cases now reconstructed & registered)

All four lost cases are rebuilt. grassgrowth/oxygenstress/surfacewater are
byte-identical; salinitystress runs with a documented ~3% solute-concentration
residual (see below). Engine fix landed: swinco=3 + numerical heat now seeds
tsoil from `[heat].tsoil_init` (the tsoil_file validator had no loader).

### Former remaining work (now done)

- **salinitystress**: swinco=3 (convert swap.ini h-profile -> `h_file` CSV),
  swap.irg (many fixed IRDATE/IRDEPTH/IRCONC/IRTYPE -> `[[irrigation.fixed_events]]`),
  full WOFOST-potato conversion (saltfarmtexel cultivar, 318-line diff incl. tables),
  swsalinity=1 (saltmax=0.732, saltslope=0.0868). RISK: the wofost salinity path
  may diverge like the cropfixed swsalinity1 xfail (salinity->cml feedback) — verify;
  register known_divergence if so.
- **surfacewater**: map the swdra=2 surface-water-management TOML schema (NRSRF=2
  subsurface levels with rdrain/rinfi/rentry/rexit/widthr/taludr; SWSRF=2/SWSEC=2
  simulated water level; WLACT/OSSWLM; NMPER=28 management periods IMPER/IMPEND/
  SWMAN/WSCAP/WLDIP/INTWL). The modern surfacewater_state init (swsec==2) supports
  it; the reader schema needs to be located and the periods table converted.
