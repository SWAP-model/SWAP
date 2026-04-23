# src/io

## Responsibility

All persistent I/O. Input readers: the TOML configuration reader
(`readswaptoml.f90`) and its drainage sibling (`readdrainagetoml.f90`),
the legacy fixed-format reader kept for backwards compatibility
(`readswap.f90`), and the meteo reader (`readmeteo.f90`). Output
writers: the CSV output backend (`swap_csv_output.f90`), the main
SWAP output dispatcher (`swapoutput.f90`), and the macropore-specific
output (`macroporeoutput.f90`). The TOML readers use the `toml-f`
subproject.

## Public interface

- `ReadSwapToml_state` — TOML main-input reader entry point.
- `ReadDrainageToml_state` — drainage TOML reader entry point.
- `MacroPoreOutput` — macropore-output writer.
- `csv_out`, `csv_out_tz` — CSV output writer and its time-zone variant.

The legacy `readswap.f90` and `swapoutput.f90` use traditional Fortran
fixed-form style without explicit `public ::` lists; their entry points
are called from `src/core/swap.f90`.

## Dependencies

Imports from:

- `src/core/` — `swap_array_dimensions`, `swap_log`, `swap_state_mod`,
  `variables`.
- `src/utils/` — `array_utils`, `soilhydraulics_utils`,
  `surfacewater_utils`.
- `src/soil/` — `soilgrid_mod`, `soilhydraulics_mod`,
  `soilwaterbalance_mod`, `doln`.
- `src/crop/` — `irrigation_mod`, `oxygenstress_mod`.
- `src/drainage/` — `drainage_mod`, `drainage_state_mod`,
  `surfacewater_mod`.
- `src/macropore/` — `macropore_mod`.
- `src/atmosphere/` — `meteodt_mod`.
- `src/heat/` — `temperature_mod`, `frozencond_mod`.
- `toml-f` (subproject) — `tomlf`, `tomlf_type`.

I/O is the one subdir that legitimately depends on most physics
domains, because the output writers sample values from every domain's
state.
