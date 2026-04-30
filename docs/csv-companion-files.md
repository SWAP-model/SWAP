# CSV companion files

The SWAP TOML pipeline keeps short tables (a handful of rows) inline
in `swap.toml`, but cases with hundreds of rows author them as CSV
companion files staged alongside the TOML.

## When to use which

| Source                        | Inline TOML       | CSV companion    |
| ----------------------------- | ----------------- | ---------------- |
| 2-row Mualem-van Genuchten    | `osat = [...]`    | -                |
| 5-row vertical discretization | `hsublay = [...]` | -                |
| 11-row tsoil_init             | `tsoil_init = []` | -                |
| 125-row gwl_table (grass)     | (currently inline; to migrate) | `*.csv`          |
| 585-row fixed irrigation      | -                 | `*.irg.csv`      |
| Multi-year met series         | -                 | `*.csv` (TBD)    |

Heuristic: above ~30 rows, prefer CSV. The TOML stays readable and
diffs/reviews of the swap.toml don't drown in tabular data.

## Authoring a CSV file

* Column 1: ISO-like date `YYYY-MM-DD`. Decoded to days-since-1900
  by the reader.
* Columns 2..N: `real(real64)` values. Plain decimal; comma-separated.
* Header row (e.g. `date,depth,conc,type`) is permitted — its first
  field is non-numeric and non-date so the reader skips it.
* Lines starting with `#` are comments and are skipped.
* Blank lines are ignored.

Example (`salinitystress.irg.csv`):

```
date,depth,conc,type
2012-04-17,14.40,0.401,1
2012-04-18,14.40,0.401,1
...
```

## Wiring in the schema

Author one optional `*_file` slot on the relevant config type and
read it via `get_optional_string_with_default`. The path is relative
to `swap.toml` (resolved at adapter time, so parsers stay free of
side effects). The strangler adapter calls
`csv_reader_mod%read_csv_date_reals(path, ncols, table, errs)` and
unpacks `table(:, 1)` (days-since-1900) + `table(:, 2..)` into the
legacy globals.

Validators should reject the case where multiple sources are set
(e.g. inline `fixed_events`, `fixed_events_file`, and the legacy
`irgfil` are mutually exclusive).

## Staging

Both `tests/swap-cases/run_case.sh` and
`tests/regression/test_output_regression.py` glob `*.csv` from
`tests/swap-cases/toml/<case>/` and copy them into the case workdir
alongside `swap.toml`. The submodule's `.gitignore` whitelists
`toml/**/*.csv` so tracked CSV companions coexist with the existing
`*.csv` ignore for SWAP run outputs.

## Reader contract

Module: `src/io/csv_reader.f90`. One public sub:

```fortran
subroutine read_csv_date_reals(path, ncols_expected, table, errors)
   character(len=*),          intent(in)    :: path
   integer,                   intent(in)    :: ncols_expected
   real(real64), allocatable, intent(out)   :: table(:,:)
   type(error_collection_t),  intent(inout) :: errors
```

* `ncols_expected` does NOT count the date column. So a CSV with
  4 columns (`date,depth,conc,type`) is read with `ncols_expected=3`
  and `table` comes out shape `(nrows, 4)` — col 1 = days-since-1900,
  cols 2..4 = the three reals.
* Errors append to `errors`; the table is left unallocated when fatal.
* Missing files append `ERR_IO_READ_FAILED`.
* Malformed cells (non-numeric in a real column, malformed date)
  append `ERR_PARSE_TYPE_MISMATCH`.

## Current users

| Reader                       | Schema slot                        | Case               |
| ---------------------------- | ---------------------------------- | ------------------ |
| `read_csv_date_reals(.., 3)` | `[irrigation].fixed_events_file`   | salinitystress     |

Future: `[meteorology].file = "<met>.csv"` will use the same reader
once the legacy ttutil-based `.met` cache is retired.
