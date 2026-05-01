---
title: "ADR 0012 — CSV companion files for tabular inputs"
date: 2026-05-01
status: accepted
---

# ADR 0012: CSV companion files for tabular inputs

## Context

SWAP's legacy `.swp` pipeline reads tabular time-series inputs (boundary
conditions, drainage water levels, irrigation schedules, etc.) through the
TTutil `rdinit`/`rdfdor`/`rdatim` reader stack. Each file uses a bespoke
keyword-value format; the schema is implicit in the reader code.

When the TOML pathway was introduced the question arose: which file formats
carry over, and which get replaced?  Replacing all legacy formats at once
would be a large-scope change; keeping them all would anchor the TOML path
to TTutil forever.

The compromise pattern that emerged organically during Phase 4f was
**CSV companion files**: structured plain-text tables (header row +
comma-separated values) read through a single `read_csv_table` subroutine.
Callers declare the expected column headers; the reader validates them and
returns a `real(real64)` matrix.  Errors flow through `error_collection_t`
so a malformed file is reported at pipeline startup rather than mid-run.

## Decision

For every tabular input that is newly wired into the TOML pathway, use CSV
companion files and `read_csv_table` as the canonical format and reader.
The legacy TTutil readers are left in place for the `.swp` pathway and are
not called from the new TOML adapter code.

**Format rules:**

- First line: comma-separated column names, lowercase, trimmed.
- If column 1 is named `date`, values must be ISO 8601 dates
  (`YYYY-MM-DD`); the reader converts them to days since JD 2415020
  (1899-12-31), matching the internal SWAP epoch used by `DTARDP`.
- If column 1 is named `datetime`, values must be `YYYY-MM-DD HH:MM:SS`;
  the reader converts them to fractional days on the same epoch.
- All other columns are `real64`.
- Lines starting with `#` and blank lines are ignored.

**`read_csv_table` public interface:**

```fortran
subroutine read_csv_table(path, expected_header, table, errors)
   character(len=*),          intent(in)    :: path
   character(len=*),          intent(in)    :: expected_header(:)
   real(real64), allocatable, intent(out)   :: table(:,:)
   type(error_collection_t),  intent(inout) :: errors
```

Header matching is strict-positional: every column name must match exactly
(case-insensitive, trimmed). Column count must match exactly.  Any mismatch
appends a fatal error and leaves `table` unallocated.

**File placement (TOML test cases):**

CSV companion files belong in the same directory as the TOML input files
(`tests/swap-cases/toml/<case>/`). The regression test runner copies all
`*.csv` from that directory to the working directory before executing.

## Inputs migrated to CSV companions (Phase 4f)

| Input | Variable | File naming convention |
|---|---|---|
| Channel water level tables | `owltab_dat` | `<case>.owltab.csv` (per drainage level) |
| Groundwater level forcing | `gwltab_dat` | `<case>.gwl.csv` |
| Aquifer head forcing | `haquiftab_dat` | `<case>.haquif.csv` |
| Drainage flux forcing (qbot2) | `qbot2tab_dat` | `<case>.qbot2.csv` |
| Drainage flux forcing (qbot4) | `qbot4tab_dat` | `<case>.qbot4.csv` |
| Head–flux relationship (qhbot) | `qhbottab_dat` | `<case>.qhbot.csv` |
| Head at bottom boundary (hbot5) | `hbot5tab_dat` | `<case>.hbot5.csv` |
| Fixed irrigation schedule | `irg_events` | `<case>.irg.csv` |
| Meteorology (daily) | `metcsv_dat` | `<station>.csv` |

## Consequences

Positive:

- One reader (`read_csv_table`) replaces ~8 distinct TTutil call sites in
  the adapter; all share the same header-validation and error-collection
  path.
- File format is human-readable, toolable, and version-control friendly.
- The adapter never opens a file through TTutil, keeping the two pathways
  fully independent.

Negative:

- CSV is whitespace-sensitive for the header row; a trailing space in a
  column name causes a mismatch.  The `trim`+case-fold in `validate_header`
  absorbs most authoring mistakes, but the strict-positional check means
  column reordering is a breaking change.
- The legacy TTutil files (`.crp`, `.ini`, per-year `.YYY` meteo) remain
  as-is until explicitly migrated; callers must know which format applies
  based on the pathway.

## Revisit trigger

If a tabular input grows beyond simple row × column layout (nested groups,
optional sections, typed enums), a TOML sub-file should be considered
instead of extending the CSV convention.
