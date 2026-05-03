# CSV companion files

The SWAP TOML pipeline keeps short tables (a handful of rows) inline; long tables
(hundreds of rows) are authored as CSV companion files staged alongside `swap.toml`
and read by a single `read_csv_table` sub at adapter time.

## When to use which

Heuristic: above ~30 rows, prefer CSV. The TOML stays readable and diffs/reviews of
`swap.toml` don't drown in tabular data.

| Source                            | Inline TOML        | CSV companion   |
| --------------------------------- | ------------------ | --------------- |
| 2-row Mualem-van Genuchten        | `osat = [...]`     | -               |
| 5-row vertical discretization     | `hsublay = [...]`  | -               |
| 126-row grassgrowth gwl table     | -                  | `*.csv`         |
| 585-row fixed irrigation schedule | -                  | `*.irg.csv`     |
| Multi-year met series             | -                  | `*.csv`         |

## Authoring a CSV file

- First non-comment, non-blank line is the header. Strict positional, lowercase,
  comma-separated, trimmed.
- Column 1 is decoded as ISO date `YYYY-MM-DD` → days-since-1900 when the column
  name is `date`. Otherwise parsed as `real(real64)`.
- Columns 2..N are always `real(real64)`. Plain decimal; no quoting.
- `#`-prefixed lines and blank lines are skipped.
- File header MUST be lowercase. The reader rejects `DATE,GWL` with
  `ERR_PARSE_HEADER_MISMATCH`.
- Empty file (header only, no data rows) is valid and yields a 0-row table.

Example (`salinitystress.irg.csv`):

```
date,depth,conc,type
2012-04-17,14.40,0.401,1
2012-04-18,14.40,0.401,1
...
```

## Reader contract

Module: `src/io/csv_reader.f90`. One public sub:

```fortran
subroutine read_csv_table(path, expected_header, table, errors)
   character(len=*),          intent(in)    :: path
   character(len=*),          intent(in)    :: expected_header(:)
   real(real64), allocatable, intent(out)   :: table(:,:)
   type(error_collection_t),  intent(inout) :: errors
end subroutine
```

- `expected_header` is the schema declared by the caller. Length must match the
  file's header column count exactly.
- If `expected_header(1) == 'date'`, col 1 is parsed as ISO date.
- Output `table` has shape `(nrows, size(expected_header))` on success, or is
  unallocated on any error.

## Errors

| Condition                                       | Error code                  |
| ----------------------------------------------- | --------------------------- |
| File missing                                    | `ERR_IO_OPEN_FAILED`        |
| File contains only blank/comment lines          | `ERR_PARSE_MISSING_HEADER`  |
| Header column count mismatch                    | `ERR_PARSE_HEADER_MISMATCH` |
| Header column name mismatch (incl. uppercase)   | `ERR_PARSE_HEADER_MISMATCH` |
| Row column count mismatch                       | `ERR_PARSE_ROW_SHAPE`       |
| Cell parse failure (real or date)               | `ERR_PARSE_TYPE_MISMATCH`   |

## Schema slots and validator wiring

The caller wires a `*_file` slot per sub-mode and a validator rule that requires it
when the matching switch is set. Current users:

| Slot               | Section             | Header           | Required when                        |
| ------------------ | ------------------- | ---------------- | ------------------------------------ |
| `gwl_file`         | `[bottom_boundary]` | `date,gwl`       | `swbotb = 1`                         |
| `qbot2_file`       | `[bottom_boundary]` | `date,qbot`      | `swbotb = 2` AND `sw2 = 2`           |
| `haquif_file`      | `[bottom_boundary]` | `date,haquif`    | `swbotb = 3` AND `sw3 = 2`           |
| `qbot4_file`       | `[bottom_boundary]` | `date,qbot`      | `swbotb = 3` AND `sw4 = 1`           |
| `qhbot_file`       | `[bottom_boundary]` | `htab,qtab`      | `swbotb = 4` AND `swqhbot = 2`       |
| `hbot5_file`       | `[bottom_boundary]` | `date,hbot`      | `swbotb = 5`                         |
| `fixed_events_file`| `[irrigation]`      | `date,depth,conc,type` | `swirfix = 1` (long-form)      |
| `h_file`           | `[soil.initial]`    | `z,h`            | `swinco = 3`                         |
| `tsoil_file`       | `[soil.initial]`    | `z,tsoil`        | `swinco = 3` AND `heat.swhea = 1` AND `heat.swcalt = 2` |
| `cml_file`         | `[soil.initial]`    | `z,cml`          | `swinco = 3` AND `solute.swsolu = 1` |

Validators should reject cases where multiple sources are set (e.g. inline
`fixed_events`, `fixed_events_file`, and the legacy `irgfil` are mutually exclusive).

## Path resolution and staging

- Companion CSVs are referenced by basename in `swap.toml`. The adapter resolves
  paths relative to the working directory at adapter time.
- The case working directory is `tests/swap-cases/toml/<N>.<case>/`. It is
  self-contained — every file SWAP reads at runtime lives there: `swap.toml`,
  `swap.dra.toml`, `*.crp.toml`, all `*.csv` companions, and
  `swap_linux.swp.template` (staged to `swap.swp` per run). For type=1
  (cropfixed, Phase 1, case 6 surfacewater), type=2 (cropwofost,
  Phase 2, case 5 salinitystress), and type=3 (cropgrass, Phase 3,
  cases 4 oxygenstress + 2 grassgrowth) rotations the `.crp` ASCII file
  is no longer staged. Phase 4 (case 1 hupselbrook) integrates all three
  crop types in one rotation.
- `tests/swap-cases/run_case.sh` runs SWAP in that directory in-place;
  `tests/regression/test_output_regression.py` copies it to a temp dir for
  parallel-safe execution. Neither tool reads from the legacy `<N>.<case>/`
  directories — those are reserved for the legacy reference binary.

## Migration history

- `salinitystress.irg.csv` — first user (irrigation `fixed_events_file`).
- `grassgrowth.gwl.csv` — replaced inline 125-row `gwl_table` (Phase 4f cleanup).
- `oxygenstress.haquif.csv` — replaced inline `haquif_table`.
- `surfacewater.haquif.csv` — replaced inline `haquif_table`.
- `salinitystress.ini.{h,tsoil,cml}.csv` — replaced legacy ASCII `swap.ini` profile blocks (Phase 4f cleanup, post-CSV-meteo). The legacy `[soil].inifil` slot was removed; `[soil.initial]` is now the canonical SWINCO=3 schema.
- `surfacewater/grass.crp.toml` — replaced legacy ASCII `grass.crp` for type=1 (cropfixed) rotations, Phase 1 of the `.crp` port. Schema 1:1 with legacy `readcropfixed`; runtime narrow per ADR 0015. Cache-driven dispatch on `crop_config_global%rotation_loaded(icrop)` per ADR 0016, with sibling-reader dispatch around `ArableLandGerm(1)` per ADR 0017. Phases 2/3/4 will extend coverage to types 2/3 + integration via case 1 hupselbrook.
- `salinitystress/potatod.crp.toml` — replaced legacy ASCII `potatod.crp` for type=2 (cropwofost) rotations, Phase 2 of the `.crp` port. Schema 1:1 with legacy `readwofost` (extended in Phase 2 with `wofost_soybean_t`, `wofost_bulb_t`, `wofost_nutrient_t` sub-types and `swrdc` field for full coverage). Runtime narrow per ADR 0015 — soybean/bulb/nutrient/CO2/scheduling and `swdrought=2`/`swoxygen=2`/`swinter=2`/`swcompensate≠0`/`swharv=1`/`swsalinity=2`/`swrdc=1` are validator-rejected. Cache-driven dispatch via `crop_config_global%rotation_loaded(icrop)` (ADR 0016); ArableLandGerm sibling dispatch (ADR 0017) extended to type=2. Init signature `cropwofost_init_from_config(cfg, icrop, FraDeceasedLvToSoil)` mirrors Phase 1's `lcc` pattern for the local SAVE in `wofost()`. Phase 3 covers type=3 (cropgrass); Phase 4 integrates all three types via case 1 hupselbrook.
- `oxygenstress/grassd.crp.toml` and `grassgrowth/grassd.crp.toml` — replaced legacy ASCII `grassd.crp` for type=3 (cropgrass) rotations, Phase 3 of the `.crp` port. Two case files because case 4 and case 2 differ on `swoxygen` (2 Bartholomeus vs 1 Feddes), `swcompensate` (1 Jarvis vs 0), and `swharvest` (1 DM-threshold vs 2 fixed-date). Schema 1:1 with legacy `readgrass` (extended in Phase 3 with ~55 scalar fields and ~15 tables, plus a Bartholomeus sub-section under `[oxygen_stress.bartholomeus]`). Runtime narrow per ADR 0015 — `swoxygen=2 swoxygentype=2`/`swcompensate=2`/`swinter ∈ {2,3}`/`swdrought=2`/`swsalinity≠0`/`swco2=1`/`swlossgrz=1`/`swlossmow=1`/`seqgrazmow ∈ {1,3}`/`swrd=1`/`swcf=3`/`swrdc=1`/`schedule=1`/`swtsum=2` are validator-rejected (`swrd=3` later un-stub-errored for Phase 4 hupselbrook). Cache-driven dispatch via `rotation_loaded(icrop)` (ADR 0016). The 14 `intent(out)` args of legacy `readgrass` (`swharvest`, `dmharvest`, `daylastharvest`, `dmlastharvest`, `swdmmow`, `maxdaymow`, `swlossmow`, `swlossgrz`, `swdmgrz`, `maxdaygrz`, `dmgrazing`, `LSDb`, `tagprest`, `swhydrlift`) are local SAVEs in `grass()` and assigned in the dispatch block prior to calling `cropgrass_init_from_config(cfg, icrop)` — chosen over a 16-arg init signature for clarity. `dateharvest` (t1900-relative timestamps) is computed from `mowing_dates` (DOY floats) by walking simulation years using `yearmeteo`. No sibling readers per ADR 0017 audit.

## Future users (deferred)

- `[meteorology].file = "<met>.csv"` — same reader once `.met` retirement lands.
- `[heat]`, `[solute]`, `[surface_water_ponding]` tables — same reader, separate specs.
