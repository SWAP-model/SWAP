# ADR 0044 — Typed CSV record tables

**Status:** Accepted (2026-05-28)
**Pilot:** `src/io/csv/meteo_csv.f90` — meteo daily/detail/rain-events

## Context

Post-strangler, CSV companion files were cached in `(:,:)` arrays on state
(e.g., `state%atmosphere%metcsv_dat(:,:)`). Schemas lived in three places:
the loader's `expected_header` argument, the consumer's column indices,
and any per-column unit/scaling code. Three points to keep in sync.

## Decision

Each CSV family gets a typed table module under `src/io/csv/<family>.f90`,
defining:

- A row record type (`<family>_row_t`) with one named field per column.
- A table type (`<family>_table_t`) holding `rows(:)`, `is_loaded :: logical`,
  and methods `load(path, errors)` / `year_window(year, i1, i2)`.

Rules:
- `load` wraps the generic `csv_reader_mod` and runs validation inline.
  No separate `validate` method — callers must not get an unvalidated table.
  Where a family has no validation baseline in the legacy code (e.g. the
  meteo detail and rain_events tables), the load procedure may omit
  validation, but must do so with an explicit inline comment naming the
  legacy procedure and the defer intent. The asymmetry is intentional:
  validation is added with byte-identical regression evidence, not eagerly.
- Unit conversions stay at the *consumer* read site (kJ→J etc.), not in
  `load`. Loader is a pure parser.
- `is_loaded` distinguishes "not requested" from "loaded but empty".
- Date columns stay `real(real64)` days-since-1900 (csv_reader convention).
  A typed `date_t` is a separate arc.

## Consequences

+ Schema lives in one place per family.
+ Consumers read by name (`rows(i)%rad`), not by column index.
+ Validation is inseparable from loading.
+ Each new CSV family is a self-contained module (~250-400 L incl. tests).

- Adding a typed table requires touching 3 places (module, meson.build,
  testSuites.inc). Same as adding any pFUnit suite.

## Pilot scope and follow-ups

Pilot: meteo (daily, detail, rain events). Established 2026-05-28.

Six more CSV families to migrate using this template:
- bottom-boundary tables (gwl, qbot, haquif, hbot, qhbot)
- drainage owltab
- irrigation fixed-events + SSDI
- nutrients amendment-events
- soil initial (h-profile, cml-profile, tsoil)

Each migration closes a chunk of `state%cfg%X` reads in its subsystem,
contributing to the broader `state%cfg` retirement arc.

## Status update 2026-05-28

- All 7 originally-scoped CSV families now using the pattern: meteo
  (daily + detail + rain), drainage (owl), irrigation (fixed + ssdi),
  nutrients (amendment), boundary (gwl + qbot + haquif + qhbot + hbot),
  soil initial (h_profile + cml_profile). Note: `tsoil` was descoped
  from the soil-initial family because the heat subsystem already handles
  its initial profile via `read_heat_toml` (`cfg_heat%tsoil_init`), not
  as a separate CSV file.
- 14 typed table types across 6 modules
  (`src/io/csv/{meteo,drainage,irrigation,nutrients,boundary,soil_init}_csv.f90`).
- `days_since_1900` helper extracted to `src/io/csv/csv_common.f90` and
  reused across all loaders; `days1900_to_md` likewise. Previously
  duplicated in `meteo_csv`, `csv_reader`, `readmeteo`, and
  `toml_field_helpers` (4 copies of the forward formula, 1 of the inverse).
- Header-array length normalized to `character(len=14)` across all loaders
  (longest column name across all families is 14: `volat_fraction` in
  nutrients_csv; padding shorter names is safe because
  `csv_reader_mod%validate_header` trims before comparison).
- W3 fix (config-mutation in soilwater_state_init) landed in Family 5:
  `config_soil` now `intent(in)`, `state%soilwater%h_init` typed table
  owns the swinco=3 h-profile.
- W4 fix (duplicate h_file CSV read) landed in Family 5.
- Next arc: `state%cfg` pointer retirement — typed CSV tables on state
  make this tractable.
- Several `seed_state_from_config.f90` residuals (~10 cross-subsystem
  writes) still standing — also part of the `state%cfg` retirement arc.
