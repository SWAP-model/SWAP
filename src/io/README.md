# src/io

## Responsibility

All persistent I/O. Input readers: the TOML configuration reader
(`readswaptoml.f90`) and its drainage sibling (`readdrainagetoml.f90`),
the legacy fixed-format reader kept for backwards compatibility
(`readswap.f90`), and the meteo reader (`readmeteo.f90`). Output
writers: the CSV output backend (`csv_output.f90`, writing via the
`csv_writer.f90` primitive, variable schema from `output_registry.f90`),
and the macropore-specific output (`macroporeoutput.f90`). The legacy
`swapoutput.f90` has been removed. The TOML readers use the `toml-f`
subproject.

## Public interface

- `ReadSwapToml_state` — TOML main-input reader entry point.
- `ReadDrainageToml_state` — drainage TOML reader entry point.
- `MacroPoreOutput` — macropore-output writer.
- `csv_output_init`, `csv_output_step`, `csv_output_finalize` — CSV
  output entry points (scalar `result_output.csv` and time-depth
  `result_output_tz.csv`).

The legacy `readswap.f90` uses traditional Fortran fixed-form style
without an explicit `public ::` list; its entry points are called from
`src/core/swap.f90`.

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
