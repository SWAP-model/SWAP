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
