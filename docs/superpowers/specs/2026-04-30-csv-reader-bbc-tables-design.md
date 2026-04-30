# Unified CSV Reader + BBC Tabular Companion Files — Design

**Date:** 2026-04-30
**Phase:** 4f cleanup (post-strangler-swap, pre-`.met` retirement)
**Status:** approved

## Goal

Establish a single CSV reader as the universal entry point for tabular SWAP input files, and migrate every tabular bottom-boundary input onto it. Eliminate inline 2D arrays for bottom-boundary tables and the legacy `bbcfil` external-file slot. Lay the contract that subsequent specs (`.met`, CO₂, heat, solute, ponding) reuse without modification.

## Context

Phase 4f introduced a TOML-based input pipeline. Three artefacts shipped recently are direct precursors to this spec:

- `read_csv_date_reals(path, ncols, table, errors)` in `src/io/csv_reader.f90` — first CSV reader, no header validation, used only by irrigation `fixed_events_file`.
- Inline `gwl_table` / `haquif_table` 2D-array slots in `[bottom_boundary]` — added during the `.bbc` decommission to carry SWBOTB=1 and SWBOTB=3 sw3=2 tables in `swap.toml`.
- Legacy `bbcfil` schema slot — kept as a forward-compat marker, never wired into the new path.

The friction these create:
- Two reader signatures (one validating, one not) for what is the same problem.
- Inline 2D arrays inflate `swap.toml` (grassgrowth's `gwl_table` is 125 rows).
- Header-less CSVs make column-order mistakes silent runtime errors.
- Per-table schema slots (`gwl_table`, `haquif_table`, future `qbot2_table`, `qbot4_table`, ...) multiply.

The unified rule replaces all of this with one reader, one slot pattern, one validation surface.

## Architecture

**Single reader, single contract:**

```fortran
subroutine read_csv_table(path, expected_header, table, errors)
   character(len=*),          intent(in)    :: path
   character(len=*),          intent(in)    :: expected_header(:)
   real(real64), allocatable, intent(out)   :: table(:,:)
   type(error_collection_t),  intent(inout) :: errors
end subroutine
```

`read_csv_table` is the SWAP CSV reader. Every tabular input file goes through it: irrigation events, bottom-boundary tables, meteo forcings (next spec), CO₂ (later), heat / solute / ponding tables (later).

**Reader behaviour:**
- Open `path`. Missing file → `ERR_IO_OPEN_FAILED`.
- Skip `#`-prefixed lines and blank lines.
- First non-comment, non-blank line is the header. Compare against `expected_header` strict-positional, lowercase, no extras, no missing. Mismatch → `ERR_PARSE_HEADER_MISMATCH`.
- For each subsequent non-comment, non-blank line: split on `,`, expect `size(expected_header)` fields. Row-shape mismatch → `ERR_PARSE_ROW_SHAPE`.
- Column 1: if `expected_header(1) == 'date'`, parse as ISO `YYYY-MM-DD` → days-since-1900 (real64). Otherwise parse as real64.
- Columns 2..N: parse as real64.
- Cell parse failure → `ERR_PARSE_TYPE_MISMATCH` citing row and column.
- Output: `table(nrows, ncols)`. Empty file (header only, zero data rows) → `table(0, ncols)`, no error.
- Any error short-circuits parsing and leaves `table` unallocated.

The existing `read_csv_date_reals` is **removed**. The salinitystress irrigation path (which currently calls it) is migrated to `read_csv_table` as part of this work.

**Schema slots in `[bottom_boundary]`:**

Six new `character(len=:), allocatable :: *_file` slots:
- `gwl_file` — SWBOTB=1, header `date,gwl`
- `qbot2_file` — SWBOTB=2 + SW2=2, header `date,qbot`
- `haquif_file` — SWBOTB=3 + SW3=2, header `date,haquif`
- `qbot4_file` — SWBOTB=3 + SW4=1, header `date,qbot`
- `qhbot_file` — SWBOTB=4 + SWQHBOT=2, header `htab,qtab`
- `hbot5_file` — SWBOTB=5, header `date,hbot`

Removed slots: `gwl_table`, `haquif_table`, `bbcfil`.

**Validator:**

`bottom_boundary_config_t%validate` adds a per-sub-mode rule: if the chosen SWBOTB sub-mode requires a table, the matching `*_file` must be non-empty. Six rules total, each emitting `ERR_VALIDATE_REQUIRED` with a message naming the offending field. Validation fires before any I/O.

**Adapter:**

`populate_bottom_boundary_*` in `src/io/toml/config_to_variables.f90` is the only call site for `read_csv_table` for BBC. Each sub-mode opens its CSV, validates the table shape, and copies into the relevant legacy globals (`gwlst/gwlinp`, `qbottab`, `haquif`, `qbot4`, `qhtab/htab`, `hbot5tab`).

**Path resolution:** companion CSVs are referenced by basename and resolved from the current working directory at adapter time, matching the convention used by `*.crp.toml`. The runtime stages all companion files into the case directory before SWAP starts. Per-section path policy will be revisited later.

## Components

| Area | File | Change |
|---|---|---|
| Reader | `src/io/csv_reader.f90` | Replace `read_csv_date_reals` with `read_csv_table`. Header parser, strict positional match, per-column dispatch (date vs real for col 1). |
| Reader test | `tests/unit/io/test_csv_reader.pf` | Rewrite to cover the 13 cases listed under Testing. |
| Reader fixtures | `tests/unit/io/fixtures/csv_*.csv` | New / updated fixture set covering happy, malformed, edge cases. |
| BBC schema | `src/config/bottom_boundary_config.f90` | Add 6 `*_file` slots. Remove `gwl_table`, `haquif_table`, `bbcfil`. New validator rule per sub-mode. |
| BBC schema test | `tests/unit/config/test_bottom_boundary_config.pf` | Update validator tests for new rule; remove tests for removed slots. |
| BBC reader | `src/io/toml/read_bottom_boundary_toml.f90` | 6 `get_optional_string_with_default` calls. Remove `read_array_2d` calls for table slots; remove `bbcfil` read. |
| Adapter — BBC | `src/io/toml/config_to_variables.f90` | Rewrite `populate_bottom_boundary_*` blocks to call `read_csv_table` per sub-mode. |
| Adapter — irrigation | `src/io/toml/config_to_variables.f90` | Migrate `fixed_events_file` call from `read_csv_date_reals` to `read_csv_table` with `expected_header=['date','depth','conc','type']`. |
| Error codes | `src/utils/error_mod.f90` | Add `ERR_PARSE_MISSING_HEADER`, `ERR_PARSE_HEADER_MISMATCH`, `ERR_PARSE_ROW_SHAPE`, `ERR_VALIDATE_REQUIRED`. |
| Case data | `tests/swap-cases/toml/2.grassgrowth/swap.toml` | Replace inline 125-row `gwl_table` with `gwl_file = "grassgrowth.gwl.csv"`. |
| Case data | `tests/swap-cases/toml/2.grassgrowth/grassgrowth.gwl.csv` | New, header `date,gwl` + 125 rows extracted from the inline table. |
| Docs | `docs/csv-companion-files.md` | Update to reflect the unified reader contract and the BBC slot list. |

## Data flow

```
swap.toml  ──load_swap_config──▶  swap_config_t                    (slot: <opt>_file = "foo.csv")
                                       │
                                       ▼  validate()
                              required-file check fires
                                       │
                                       ▼  config_to_variables (adapter)
                              read_csv_table(<opt>_file, header, table, errors)
                                       │                  │
                                       ▼                  ▼
                              variables%<table_name>     errors collected; abort if any
```

`load_swap_config` does not open companion files — it only records paths. Open / read / validate happens in the adapter. This preserves the load-validate-adapt separation.

## Error handling

Reader-time errors:

| Condition | Error code | Message shape |
|---|---|---|
| File not found | `ERR_IO_OPEN_FAILED` | `cannot open CSV: <path>` |
| Header line absent | `ERR_PARSE_MISSING_HEADER` | `<path>: missing header line` |
| Header column count mismatch | `ERR_PARSE_HEADER_MISMATCH` | `<path>: expected header [date,gwl], got [date,gwl,note]` |
| Header column-name mismatch | `ERR_PARSE_HEADER_MISMATCH` | `<path>: expected header [date,gwl], got [Date,GWL]` |
| Numeric cell unparseable | `ERR_PARSE_TYPE_MISMATCH` | `<path>:row 47 col 2: 'foo' not a real number` |
| Date cell unparseable | `ERR_PARSE_TYPE_MISMATCH` | `<path>:row 47 col 1: 'foo' not an ISO date` |
| Row column count mismatch | `ERR_PARSE_ROW_SHAPE` | `<path>:row 47: expected 2 fields, got 3` |
| Empty table (header only) | accepted | returns `table(0,ncols)` |

Validator-time errors:

| Condition | Error code | Message |
|---|---|---|
| SWBOTB=1 with empty `gwl_file` | `ERR_VALIDATE_REQUIRED` | `bottom_boundary.gwl_file required when swbotb=1` |
| SWBOTB=2 + SW2=2 with empty `qbot2_file` | `ERR_VALIDATE_REQUIRED` | analogous |
| SWBOTB=3 + SW3=2 with empty `haquif_file` | `ERR_VALIDATE_REQUIRED` | analogous |
| SWBOTB=3 + SW4=1 with empty `qbot4_file` | `ERR_VALIDATE_REQUIRED` | analogous |
| SWBOTB=4 + SWQHBOT=2 with empty `qhbot_file` | `ERR_VALIDATE_REQUIRED` | analogous |
| SWBOTB=5 with empty `hbot5_file` | `ERR_VALIDATE_REQUIRED` | analogous |

All errors flow through `error_collection_t`; the run aborts at the first checkpoint after any error. No silent fallbacks.

## Testing

**Reader unit tests** (`tests/unit/io/test_csv_reader.pf`):

| Test | Fixture | Asserts |
|---|---|---|
| `test_date_keyed_happy` | `csv_date_gwl.csv` | shape `(3,2)`; col 1 = days-since-1900; col 2 values match. |
| `test_real_keyed_happy` | `csv_htab_qtab.csv` | shape `(3,2)`; col 1 = htab values; col 2 = qtab values. |
| `test_wide_table_happy` | `csv_irrigation.csv` | shape `(N,4)`; header `date,depth,conc,type`; rows correct. |
| `test_header_wrong_order` | `csv_gwl_date_swapped.csv` | `ERR_PARSE_HEADER_MISMATCH`; `table` not allocated. |
| `test_header_case_mismatch` | `csv_uppercase_header.csv` | `ERR_PARSE_HEADER_MISMATCH`. |
| `test_header_extra_column` | `csv_extra_column.csv` | `ERR_PARSE_HEADER_MISMATCH`. |
| `test_header_missing` | `csv_no_header.csv` | `ERR_PARSE_MISSING_HEADER`. |
| `test_missing_file` | `does_not_exist.csv` | `ERR_IO_OPEN_FAILED`. |
| `test_malformed_real` | `csv_bad_real.csv` | `ERR_PARSE_TYPE_MISMATCH`; cites row+col. |
| `test_malformed_date` | `csv_bad_date.csv` | `ERR_PARSE_TYPE_MISMATCH`; cites row+col. |
| `test_row_shape_mismatch` | `csv_short_row.csv` | `ERR_PARSE_ROW_SHAPE`. |
| `test_empty_table` | `csv_header_only.csv` | shape `(0,2)`; no errors. |
| `test_comments_and_blanks` | `csv_with_comments.csv` | `#` and blank lines tolerated. |

**Validator unit tests** (`tests/unit/config/test_bottom_boundary_config.pf`):

Per SWBOTB sub-mode: one positive case (file slot populated → passes) and one negative (slot empty → `ERR_VALIDATE_REQUIRED` on the right field). Six pairs total. Tests for removed slots (`gwl_table`, `haquif_table`, `bbcfil`) are deleted.

**Regression** (`pixi run -e test regression`):

- **grassgrowth** — SWBOTB=1, migrated to `gwl_file`. Must remain green at 1e-2 cm tolerance.
- **salinitystress** — irrigation path migrated to `read_csv_table`. Header unchanged; numerically equivalent to baseline.
- **hupselbrook, oxygenstress, surfacewater** — must stay green; no collateral damage. Pre-existing oxygenstress and salinitystress drift remains pre-existing.

**Build verification:** `pixi run -e test test-pfunit` (unit) and `pixi run -e test check-full` (regression) green at the end of each task.

## Out of scope

- `.met` retirement — same reader, separate spec.
- CO₂ — same reader, separate spec.
- Heat / solute / ponding tables — same reader, separate spec.
- Per-year `.YYY` met files — deferred until macropore work.
- Path-resolution policy refinement — deferred. Current convention: basename in working directory, matching `*.crp.toml`.
- Compatibility shims for legacy `.bbc` / `.irg` formats — none. Translation to CSV is per-case authoring.
- Fixing oxygenstress and salinitystress regression drift — pre-existing, Phase 4f-extend territory.
