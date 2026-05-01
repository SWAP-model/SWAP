---
title: "ADR 0013 — CSV meteorology input for the TOML pathway"
date: 2026-05-01
status: accepted
---

# ADR 0013: CSV meteorology input for the TOML pathway

## Context

SWAP's legacy meteorology input uses one of two formats:

1. **Per-year files** (`<stem>.<YYY>`) — one file per year, TTutil keyword
   format, columns: `station`, `dd`, `mm`, `yyyy`, `rad`, `tmin`, `tmax`,
   `hum`, `wind`, `rain`, `etref`, `wet`.
2. **All-years file** (`<stem>.met`) — same format, all years in one file,
   read via `MeteoInOneFile` which pre-loads them into SAVE arrays.

Both formats carry a freeform `station` text field and three separate date
columns (`dd`/`mm`/`yyyy`) used only to compute a Julian day number.

The TOML pathway needs a clean meteo input that:
- Follows ADR 0012's CSV companion convention.
- Replaces the three date columns with a single ISO 8601 `date` field.
- Drops the `station` freeform string (not used in computations).
- Pre-loads all years in the adapter (consistent with the adapter's
  "resolve everything before the simulation loop" contract).

## Decision

The TOML pathway accepts meteorology via a CSV file with the schema:

```
date,rad,tmin,tmax,hum,wind,rain,etref,wet
```

- `date`: ISO 8601 `YYYY-MM-DD`; parsed to days since JD 2415020 by
  `read_csv_table`.
- `rad`: daily global radiation, **kJ m⁻² d⁻¹** (same unit as legacy
  files; converted to J m⁻² d⁻¹ internally by `MeteoCSVYear`).
- All other columns: `real64`, same physical units as the legacy format.

**Detection:** if `[meteorology].file` ends in `.csv`, `swMetCSV` is set to
1 and the adapter pre-loads all rows into the module-level `metcsv_dat`
array.  The legacy `swMetFilAll` flag is left at 0.

**Year extraction:** `ReadMeteoYear` (called once per simulation year)
checks `swMetCSV` first.  When set, it calls `MeteoCSVYear` and jumps past
the TTutil reader block.  `MeteoCSVYear` scans `metcsv_dat(:,1)` for rows
whose date falls within the target year (Julian day boundaries computed from
the same `jday()` helper and JD 2415020 epoch used elsewhere), populates
`arad`, `atmn`, `atmx`, `ahum`, `awin`, `arai`, `aetr`, `wet`, `ad`, `am`,
then sets `daynrfirst`/`daynrlast`/`timjan1` via `DTARDP` — identical to
what `MeteoInOneFile(2)` and the TTutil reader do.

**File placement:** CSV files live in the TOML case directory
(`tests/swap-cases/toml/<case>/`); the regression runner stages them
automatically.

**Legacy files unchanged:** ASCII case directories keep their `.met` and
`.YYY` files.  The legacy `swap420` reference binary is unaffected.

## Sub-daily (detail) meteorology

The `datetime_keyed` branch in `read_csv_table` (`parse_iso_datetime`)
supports `YYYY-MM-DD HH:MM:SS` timestamps and returns fractional days on
the same epoch.  The detail meteo CSV path (`swmetdetail=1`) is a direct
extension of this ADR; its schema and integration are left for a follow-up
phase.

## Consequences

Positive:

- Five TOML regression cases no longer depend on `.met`/`.YYY` files;
  the TTutil meteo reader is not reachable from any active test.
- ISO date schema is unambiguous, portable, and toolable.
- The adapter pre-loads the full dataset once; no per-year file I/O in
  the simulation loop.

Negative:

- `MeteoInOneFile` and the per-year TTutil reader remain in
  `readmeteo.f90` as unreachable-from-TOML code.  They are still needed
  by the legacy `.swp` pathway.  Deletion is deferred to Phase 5+.
- The all-years CSV must cover every year in the simulation range exactly;
  `MeteoCSVYear` calls `fatalerr_collected` if no rows are found for a
  year.  Partial-year edge cases (e.g. leap years, gaps) are not handled.

## Files changed

| File | Change |
|---|---|
| `src/io/csv_reader.f90` | Added `datetime_keyed` flag and `parse_iso_datetime` |
| `src/core/variables.f90` | Added `swMetCSV`, `nmetcsv`, `metcsv_dat` |
| `src/io/readmeteo.f90` | Added `MeteoCSVYear` + `days1900_to_md` + `jday` helpers; `goto 100` bypass |
| `src/io/toml/config_to_variables.f90` | Three-way branch: `.csv` → pre-load; `.met` → `MeteoInOneFile(1)` |
| `tests/swap-cases/toml/*/swap.toml` | `file` changed from `*.met` to `*.csv` |
| `tests/swap-cases/toml/*/` | New `<station>.csv` added to each active case |

## Revisit trigger

When the legacy `.swp` pathway is retired or `MeteoInOneFile` / the
per-year reader is deleted, the `swMetCSV` / `swMetFilAll` flags and the
`MeteoInOneFile` iTask dispatch can be removed from `readmeteo.f90`.
