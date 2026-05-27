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

## Finalized cluster set

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
