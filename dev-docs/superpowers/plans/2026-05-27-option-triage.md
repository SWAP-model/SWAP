# Phase 0 — Option-Support Triage

**Date:** 2026-05-27
Feeds the cluster fan-out in `2026-05-27-regression-coverage-expansion.md`.

## Environment facts

- `swap420` (Intel SWAP 4.2.0 Linux binary) **runs on this machine** (deps OK) and
  produces `result_output.csv`. Generate references via
  `tests/swap-cases/run_case.sh -c <case> --legacy-binary ../../tests/reference/swap420 -k`
  (runs in the **legacy ASCII** dir). NOTE: the `pixi run swap-ref` task uses `--exec`,
  which run_case.sh routes to TOML mode — wrong for a 4.2.0 binary. Use `--legacy-binary`.
- CSV columns are controlled by `INLIST_CSV` in the `.swp` (and the TOML equivalent):
  a comma-separated list of output-variable names. To assert a new-physics column, add
  its name there — but only if **both** binaries know the name.
- Output vocabulary cross-check: `SNOW`, `SSNOW`, `MELT` exist in both swap420 (`strings`)
  and the modern `output_registry`. No dedicated frost/cover-fraction output var found →
  those clusters assert via standard water-balance columns + a divergence-from-base check.

## Support classification (modern build)

`core` = source files matching `\bsw…\b` outside `src/io/toml/`; `core=0` ⇒ unsupported.

**Supported (build):** swsnow, swsublim, swfrost, swhyst, swsophy, swdrought, swcalt,
swsalinity, swinter, swgc, swqhbot, swqhr, swco2, swnrsrf, swdislay, swgerm, swharv,
swsow, swprep, swrdc, swtsum, swdivide, swcrop, swrootradius, swkmean, swdmgrz, swdmmow,
swlossgrz, swlossmow, swkimpl, swman, swtill, swrain, swmetdetail, swetsine, swinco,
swbotbc, swbotbhea, swsrf, swsec, swallo, swcompensate, swcirrthres.

**Marginal (verify at use):**
- `swdc` (solute decomposition) — `core=1`, matches only `src/config/solute_config.f90`
  (read + enum-validated, **no compute-path hit**). Treat as a read-only stub → **exclude
  from the salinity cluster** unless a compute path is confirmed during execution.
- `swgraz` (grazing) — `core=1` but has a real `if (self%swgraz==1)` branch reading
  grazing events in `cropgrass_config.f90`. Includable in grass-management; verify it
  produces grazing output when activated, else drop.

**Unsupported (excluded):** `swsp` (solute adsorption, `core=0`); macropore (ADR 0040);
`swsophy` (tabulated soil physics) — **discovered 2026-05-27**: `sptabulated.f90` is in
`src/soil/dormant/` (not compiled) and there is no TOML input path for the soil-physical
table filenames (`FILENAMESOPHY`). The triage `core=7` counted the dormant module. The
soil-tabulated cluster is dropped.

## 2026-05-27 RE-TRIAGE — stub-guard reality (binding constraint)

The grep-based support classification above **over-counted**: it credited dormant modules
and config fields as "supported." The real gate is the modern build's **runtime
stub-error guards** on the TOML path (`*_init.f90`, mirrored by config validators). A
switch-value is buildable only if it is enum-allowed AND not stub-guarded AND has a TOML
input path.

**Stub-guarded / unsupported on the TOML path (crop subsystem — all crop types):**
`swdrought=2`, `swco2=1`, `swharv=1`, `swcompensate/=0`, `swinter=2/3`, `swrd=1`,
`swcf=3`, `swrdc=1`, `swtsum=2`, `swsalinity=2` (wofost; grass rejects swsalinity/=0
entirely), grazing (`seqgrazmow/=2`/`schedule=1`), `swlossmow=1`, `swlossgrz=1`,
`swsoybean`, `swbulb`, `swoxygen=2`. Plus `swsophy` (dormant). ⇒ **These clusters are
DEAD:** drought-vanlier, grass-management, salinity-osmotic, soil-tabulated, and most of
crop-calendar.

**Genuinely supported (enum-allowed, NOT stub-guarded) — the realistic buildable set:**
- Non-crop (no stub-guards in drainage/surfacewater/soil/bottom/heat/meteo init):
  `swqhbot=1`, `swqhr=2`, `swbotbhea=2`, `swnrsrf=1/2`, `swdislay=1`, `swsec=1`,
  `swsrf` variants, `swrain=1`, `swinco=1`, `swetsine=1`. (Caveat: `swmetdetail=1` forces
  ssnow=0 — incompatible with snow; needs sub-daily meteo input.)
- Crop (not in any stub-guard list): `swgc=2` (soil cover fraction), `swgerm=1`.

Revised remaining clusters: **bottom-boundary** (swqhbot=1), **heat-bottom** (swbotbhea=2),
**extended-drainage** (swnrsrf, swdislay), **surface-water-variants** (swsec=1, swqhr=2,
swsrf), **meteo** (swrain=1, swinco=1, swetsine), **crop-modes** (swgc=2, swgerm=1).

## Finalized cluster set (SUPERSEDED by the re-triage above for clusters 4–8)

| # | Cluster | Base | Switches | Asserted columns |
|---|---|---|---|---|
| 2 | soil-hysteresis | hupselbrook | swhyst, swkimpl, swkmean(alt) | standard balance + GWL (divergence-from-base) |
| 3 | winter | cold-meteo case | swsnow, swsublim, swfrost | SNOW, SSNOW (+ balance) |
| 4 | soil-tabulated | hupselbrook | swsophy | standard balance + GWL |
| 5 | drought-vanlier | oxygenstress | swdrought=2, swrdc, swrootradius=1, swdivide=0 | TREDDRY/TREDWET (+ balance) |
| 6 | crop-calendar | fixed-crop | swsow, swprep, swharv, swgerm=1, swco2, swgc=2 | crop DM / DVS |
| 7 | grass-management | grassgrowth | swdmgrz=1, swdmmow=1, swlossgrz, swlossmow, swtsum(alt), swgraz(verify) | PGRASSDM/GRASSDM/PMOWDM/MOWDM |
| 8 | salinity-osmotic | salinitystress | swsalinity=2, swbotbc=2 (swdc dropped) | TREDSOL/CONC[...] |
| 9 | extended-drainage | surfacewater | swnrsrf, swdislay, swsec=1, swsrf(alt), swallo=2, swqhbot=1, swqhr=2, swbotbhea=2 | GWL/POND/drainage |
| 10 | meteo-detail | hupselbrook | swrain=1, swmetdetail=1, swetsine, swinco=1 | RAIN/RUNOFF/INTERC |

Clusters may still split during execution if activated switches prove mutually
incompatible in one input file (documented per case).
