# src/atmosphere

## Responsibility

Atmospheric boundary processing: daily and sub-daily meteorology ingestion
(`meteoday.f90`, `meteodt.f90`), precipitation partitioning
(`precipitation.f90`), canopy interception using the Von Hoff / Gash /
Rutter / simple MSW1 schemes (`interception.f90`), reference and actual
evapotranspiration via Penman-Monteith with soil-evaporation reduction
(`et.f90`), and snow accumulation / melt plus CN-method runoff
(`snow.f90`). The aggregated `atmosphere_state_t` holds the per-timestep
fluxes and cumulative totals that the rest of the model reads through
`swap_state_t`.

## Public interface

- `atmosphere_state_t` — aggregated atmosphere state type.
- `atmosphere_state_init` / `atmosphere_state_finalize` — lifecycle.
- `atmosphere_state_reset_cumulative` / `atmosphere_state_reset_intermediate` — reset hooks.
- `snow`, `snow_state` — snow accumulation and melt.
- `PartitionPrecipitation` — rain / snow splitting.
- `PenMon`, `reduceva`, `reduceva_state` — Penman-Monteith and soil-evaporation reduction.
- `VonHHBraden`, `Gash`, `ruttervw`, `msw1eic`, `DivIntercep` (+ `_state` twins) — interception schemes.
- `MeteoDT`, `MeteoDT_state` — sub-daily meteo.
- `ReadMeteoDay`, `ResetMetFlx`, `ProcessMeteoDay` (+ `_state` twins) — daily meteo pipeline.
- `CNmethod` — SCS curve-number runoff.

## Dependencies

Imports from:

- `src/core/` — `swap_constants`, `swap_array_dimensions`, `swap_state_mod`,
  `swap_state_sync`, `variables` (legacy bridge).
- `src/utils/` — `array_utils`, `soilhydraulics_utils`.

Does not depend on any other physics subdir — downstream consumers read
atmosphere outputs through `swap_state_t`.
