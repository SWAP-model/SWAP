# CSV output modernization arc (IO-OUT) — design

**Date:** 2026-05-26
**Scope:** `src/io/` output side + the `*Output` wiring in `src/core/swap_mod.f90` and the `.crp` writers reached from `src/crop/cropgrowth_helpers.f90`. The input side (TOML + CSV readers) is already modern and is **out of scope** here (see "Out of scope / follow-on").
**Goal:** Bring SWAP output to a single, symmetric, registry-driven CSV path: a reusable `csv_writer` primitive (mirror of `csv_reader`), a declarative output-variable registry, one orchestrating `csv_output` module that samples `state` once and feeds both the two result CSVs *and* the live C-API/BMI streams, and the full deletion of the `swapoutput.f90` graveyard. Each phase is independently bisectable and ends green.

## Background

The input pipeline is finished: `load_swap_config` + ~20 section readers produce a typed `swap_config_t`; `csv_reader.f90` reads numeric tables (`read_csv_table`, `real64` matrix, col-1 ISO date/datetime); `config_to_variables.f90` seeds typed `state%X` via per-subsystem `init`. No legacy fixed-format readers survive (the `RDinit`/`rddre`/`rdinit` mentions in `tests/regression/test_output_regression.py` comments are stale; no such subroutines exist in `src/`).

The **output** side is half-retired, not unmodernized:

- **`swap_csv_output.f90`** is the modern keeper but a ~1080-line monolith mixing four concerns: a hardcoded ~103-variable catalogue, state sampling (`set_values`), aggregation (`fill_values`: `dstor`, `baldev`, subregion sums), and CSV formatting/writing. It produces `result_output.csv` (one row/day: scalars + subregion aggregates, via `csv_out`) and `result_output_tz.csv` (time×depth profile, via `csv_out_tz`). There is **no reusable writer primitive** — asymmetric with `csv_reader`.
- **`swapoutput.f90`** (~1960 lines) is mostly dead:
  - *Live wiring:* `soilwateroutput` is the dispatch hub the main loop calls — it now only calls `csv_out`/`csv_out_tz`. `swapoutput` itself only manages the `water_balance_row` buffer + logs.
  - *Dead, uncalled:* `outinc` (`.inc`), `outrot` (`.rot`), `outtem` (`.tem`), `outage`.
  - *Empty stubs:* `SoluteOutput`, `TemperatureOutput`, `SurfaceWaterOutput`, `AgeTracerOutput` — see "The C-API stream channel" below.
  - *Still writing legacy text files:* `SnowOutput → .snw` (from `swap_mod`) and `OutCropFixed/OutWofost/OutGrass → .crp` (from `cropgrowth_helpers.f90`).

The atmosphere cleanup arc (`c29408b`, the PenMon derived-type refactor and the GR-ATM-CLEAN arc) is the style template: explicit-interface, derived-type-bundled arguments, named init/step pairs over task-selector subroutines, one responsibility per file. This arc applies the same discipline to output.

### The C-API stream channel (load-bearing constraint)

The four empty stubs are not new scaffolding — they are the gutted remains of original per-subsystem *file* writers (`SoluteOutput` Nov 2004, etc.). ADR 0009 Phase 5+ deleted their file bodies (switches forced to 0); the SS-BMI2 work then repurposed the shells to allocate/free a per-subsystem in-memory `output_row` buffer in `state`, exposed through the C API. `swap_get_output_row(stream, …)` in `src/core/swap_capi_mod.f90:213` enumerates **9 streams**:

```
swap_balance  soilwater  temperature  solute  agetracer
snow          surfacewater  crop  tillage
```

This is a **second output channel parallel to the CSV files** — same data sampled from `state`, handed to C/BMI instead of disk. Four rows are actually *filled* today: `soilwater` (`swap_csv_output`), `crop` (`cropgrowth_helpers`), `tillage` (`tillage`), and `snow` (`build_snow_output_row`, called from `SnowOutput` at `swapoutput.f90:1428`). `water_balance_row` is *intended* to be filled but its only builder lived in `outinc`, which is now dead — so `swap_balance` currently returns a zeroed row (latent bug, fixed in Phase E). The `temperature/solute/surfacewater/agetracer` streams build nothing — their `case(2)` is an explicit no-op ("buffer stays zeroed").

**Decision (REVISED 2026-05-26 — supersedes the earlier "fold + drop 4" plan).** Investigation during execution found: (a) the `crop` stream's builder lives inside the swcrp-gated `CropOutput`, which is dead in TOML mode, so `crop` is inert too; (b) `swap_balance`'s builder was the now-dead `outinc`, so it returns zeros; (c) **no consumer reads `swap_get_output_row` at all** — not the BMI tests, not the cffi-demo POC (which uses `swap_view_array` / `swap_get_scalar` / `swap_get_water_balance`, the last reading cumulative state directly, independent of the buffers). Combined with the product direction — Python is **batch** (seed → run → read results), per-step BMI is **later**, and the goal is to **unify strongly** — the per-subsystem per-step `output_row` mechanism serves no current mode.

So Phase E **retires the entire `output_row` mechanism**: the `swap_get_output_row` C function + all 9 stream cases, every `*_output_row`/`*_columns`/`*_n_cols` state field, and all `build_*_output_row`/`init|cleanup_*_buffer` helpers. The single results-production path is registry → `csv_output` → `result_output.csv`. Kept, because they are orthogonal and used: `swap_view_array`, `swap_get_scalar`, `swap_get_water_balance` (reads cumulative state directly). **Future work (separate arc):** one unified results accessor on top of the registry/`csv_output` — returns the accumulated result table in memory for batch-Python (no file needed) and extends to per-step BMI later. This replaces the deleted `swap_get_output_row`.

## Approved decisions

1. **Writer architecture:** generic `csv_writer` primitive (symmetric to `csv_reader`) + a metadata **registry** (`{name, unit, kind}`) on top. Sampling stays **centralized** in `csv_output` keyed by registry order — *not* per-variable procedure pointers (the lighter middle option, not the full descriptor-with-getter redesign).
2. **File layout:** keep **two** files — `result_output.csv` (daily scalars + subregion, with crop & snow folded in as columns) and `result_output_tz.csv` (time×depth profile). The two row shapes are genuinely different, and the regression harness already targets `result_output.csv`.
3. **`csv_reader`:** unchanged — stays numeric-only (`real64` matrix). The real addition is its symmetric counterpart, `csv_writer`.
4. **BMI streams:** option 1 + 3 (fold into `csv_output`; drop the four inert streams).

## Regression target

`tests/regression/test_output_regression.py` compares **only `result_output.csv`** — annual aggregated stats (flux sums, state means, cumulative last-values) vs JSON fixtures, tolerance `1e-2`, across the 5 live cases (hupselbrook, grassgrowth, oxygenstress, salinitystress, surfacewater). `.crp`, `.snw`, `.tem`, `.inc`, `.rot` are **not tested**. Consequences:

- Refactors that preserve the values and column set of `result_output.csv` are safe by the harness.
- Crop and snow quantities already exist as `result_output.csv` columns (DVS, LAI, PGRASSDM, GRASSDM, SNOW, SSNOW, SUBLIM, …) — so `.crp`/`.snw` are largely redundant. But because those files are untested, a one-shot **parity script** (Phase D) must confirm column coverage before deletion.

## Verification & build discipline

- After **every** phase: build green + `check-fast` (pFUnit + the regression cases) per the per-task regression gate convention. No deferring to an end-of-arc gate.
- Phases that change `src/state/*_state.f90` (removing `output_row` fields, Phase E) require `rm -rf builddir` before rebuild — incremental Meson does not propagate `.mod` deps across the swap_modern↔swap_legacy boundary.
- CSV number formatting (E-notation when `|x| < 1e-4` or `> 1e4`, else fixed 5-decimal) is **load-bearing** for the `1e-2` tolerance and must be carried verbatim from `swap_csv_output` into `csv_writer`.
- Commits are `development`-only, one focused change each, conventional-commit style with an `IO-OUT` phase tag.

## Arc shape — five phases, sequential, dependency-driven

| Phase | What it does | Key files | Est. commits |
|-------|--------------|-----------|--------------|
| **A** | Add `csv_writer.f90` primitive + pFUnit round-trip tests (write → read-back via `csv_reader`). No behavior change. | `src/io/csv_writer.f90`, `tests/pFUnit/...`, `meson.build` | 2 |
| **B** | Extract the variable catalogue into `output_registry.f90`; `swap_csv_output` consumes it for header/unit/order + inlist resolution. | `src/io/output_registry.f90`, `src/io/swap_csv_output.f90`, `meson.build` | 2–3 |
| **C** | Route `swap_csv_output`'s writing through `csv_writer_t`; merge the `_tz` module; rename `swap_csv_output → csv_output` with `csv_output_init/step/finalize`. | `src/io/csv_output.f90` (from `swap_csv_output.f90`), `src/io/swapoutput.f90` (dispatch), `meson.build` | 3 |
| **D** | Column-parity audit: confirm every `.crp`/`.snw` quantity has a registry entry; add the missing ones; ship a parity script as evidence. | `src/io/output_registry.f90`, `scripts/output_parity_check.py` | 2 |
| **E** | Retire the graveyard: delete `swapoutput.f90` (dead routines, stubs, `.snw`/`.crp` writers), delete `.crp` calls in `cropgrowth_helpers`, fold live C-API buffer-building into `csv_output`, drop the 4 inert streams, fix `water_balance_row`, rewire `swap_mod`. | `src/io/swapoutput.f90` (deleted), `src/io/csv_output.f90`, `src/core/swap_mod.f90`, `src/core/swap_capi_mod.f90`, `src/crop/cropgrowth_helpers.f90`, `src/state/{heat,solute,surfacewater}_state.f90`, `meson.build` | 4–5 |

Total: ~13–15 commits. Verification gate at the end of every phase.

## Phase A — `csv_writer` primitive

### Target

A new `src/io/csv_writer.f90`, module `csv_writer_mod`, symmetric to `csv_reader_mod`. A derived type owns the open unit and column count:

```fortran
type :: csv_writer_t
   integer :: unit  = -1
   integer :: ncols = 0
   logical :: leading_datetime = .false.
 contains
   procedure :: open   => csv_writer_open    ! open file, write '*'-prefixed meta block
   procedure :: header => csv_writer_header   ! column-name row, then unit row
   procedure :: row    => csv_writer_row      ! optional datetime string + real64 values
   procedure :: close  => csv_writer_close
end type
```

- `open(path, meta_lines, errors)` opens via `file_io` (`newunit`), writes the `* key: value` metadata comment block (Project / File content / File name / Model version / Generated at) exactly as `swap_csv_output` does today.
- `header(names, units)` writes the column-name row and the unit row. Sets `ncols`.
- `row(values)` / `row(datetime, values)` formats one data line. The number formatter is lifted verbatim from `swap_csv_output` (E vs F selection) into a private `fmt_real` so the byte output is identical.
- `close()`.

### Scope

- New module only; nothing consumes it yet.
- `meson.build`: add `csv_writer.f90` to `swap_io_sources`.

### Tests

- pFUnit: build a known `real64` matrix + header, write it through `csv_writer_t`, read it back through `read_csv_table`, assert header match and element equality within `1e-12`. This is the symmetry contract.
- pFUnit: formatter edge cases (`0.0`, `1e-5`, `1e5`, negatives) produce the same strings the legacy inline formatter produced (golden strings copied from a current run).

### Why first

Every later phase writes through this primitive; landing it standalone with its own tests de-risks the rest.

## Phase B — `output_registry`

### Current state

`swap_csv_output.f90:94-198` hardcodes the ~103-variable catalogue inline: scalar names (RAIN, EACT, GWL, DVS, LAI, …), per-node templates (`H[`, `WC[`, `TEMP[`, `CONC[`, `RWU[`, …), and per-subregion templates (`WTOT[`, `QTRANS[`, `QDRA[`, …). `set_values` fills a flat value array whose index layout is implicit in the catalogue order.

### Target

`src/io/output_registry.f90`, module `output_registry_mod`:

```fortran
integer, parameter :: OUT_SCALAR=1, OUT_NODE=2, OUT_SUBREGION=3, OUT_PROFILE=4
type :: out_var_t
   character(len=32) :: name
   character(len=16) :: unit
   integer           :: kind
end type
```

- A module-level `registry(:)` initialized from the current catalogue — the **single source** for both the value-array layout (registry index = value index) and the emitted header/units.
- `resolve_inlist(inlist_csv, sel, errors)` — validate the user's `config%output_csv%inlist` names against the registry, return selected indices in canonical order, append a typed error on any unknown name (today an unknown name is silently dropped — tightening this is in scope).
- Accessors: `var_count()`, `var_name(i)`, `var_unit(i)`, `var_kind(i)`.

### Scope

- `swap_csv_output` switches to the registry for header/unit emission, iteration order, and inlist resolution. `set_values` writes into registry indices.
- No file-format change: `result_output.csv` byte output and column set stay identical. Regression must be green/within-tol.

### Tests

- pFUnit: `resolve_inlist` filters/orders correctly; unknown name yields a fatal error in the collection; empty inlist falls back to the documented default.

## Phase C — route through `csv_writer`, rename to `csv_output`

### Target

- Merge `SWAP_csv_output` and `SWAP_csv_output_tz` into one module `csv_output` (`src/io/csv_output.f90`, renamed from `swap_csv_output.f90`), holding two `csv_writer_t` instances (scalar + profile).
- Public API replaces the `csv_out(iTask, …)`/`csv_out_tz(iTask, …)` task-selectors with named procedures:
  - `csv_output_init(state, config)` — resolve inlist via registry, open both writers.
  - `csv_output_step(state)` — compute aggregates, sample, write one row to each writer.
  - `csv_output_finalize(state)` — close both.
- All file writing now goes through `csv_writer_t`; the inline `open`/`write`/format code in `swap_csv_output` is deleted.
- `swapoutput.f90:soilwateroutput` is updated to call the three new procedures (still the dispatch hub for this phase — it is deleted in Phase E).

### Scope

- Internals retained: `compute_aggregates` (= today's `fill_values`: `dstor`, `baldev`, subregion sums) and the scalar/profile sampling (= `set_values` + the tz logic).
- `result_output.csv` / `result_output_tz.csv` byte output unchanged. Regression green/within-tol.

### Tests

- pFUnit: `compute_aggregates` reproduces `dstor`/`baldev` for a constructed state (extract the formulas into a pure helper to make them unit-testable).
- Regression: all 5 cases within tol.

## Phase D — column-parity audit (`.crp`, `.snw`)

### Goal

`.crp` and `.snw` are untested by the harness, so before deleting their writers (Phase E) we must prove the two CSVs already carry their data.

### Scope

- Enumerate the columns `OutCropFixed`/`OutWofost`/`OutGrass` write to `.crp` and `SnowOutput` writes to `.snw`; map each to a registry entry. Add any genuinely missing variable to `output_registry` (+ its sampling case in `csv_output`).
- `scripts/output_parity_check.py`: run a case that exercises each legacy file, parse the legacy `.crp`/`.snw` and the corresponding `result_output.csv` columns, assert per-day equality within tol. Committed as the evidence that Phase E's deletions lose nothing.

### Tests

- The parity script passes for a grass case (`.crp` for grass) and a snow-active case (`.snw`).
- Regression green (additions are new columns; existing values unchanged).

## Phase E — retire the graveyard

### Scope

- **Delete `src/io/swapoutput.f90` entirely:** `swapoutput`, `soilwateroutput`, `outinc`, `outrot`, `outtem`, `outage`, `SoluteOutput`, `TemperatureOutput`, `SurfaceWaterOutput`, `AgeTracerOutput`, `SnowOutput`, `OutCropFixed`, `OutWofost`, `OutGrass` and the buffer helpers.
- **Re-home live C-API buffer-building into `csv_output_step`:** `water_balance_row` (fix: it is currently built only in the dead `outinc`, so `swap_balance` returns zeros today — wire it from `compute_aggregates`); `snow` row (currently built in the to-be-deleted `SnowOutput` — move `build_snow_output_row` into `csv_output_step`, gated on `flSnow`); `soilwater` row already built in `csv_output`. `crop`/`tillage` rows stay where they are (`cropgrowth_helpers`, `tillage` — untouched by this deletion).
- **Drop the 4 inert streams** `temperature/solute/surfacewater/agetracer`: remove their `case` arms in `swap_get_output_row`, and remove the now-unused `output_row`/`*_columns`/`*_n_cols` fields from `src/state/{heat_state,solute_state,surfacewater_state}.f90` (state schema change → clean rebuild).
- **Delete the `.crp` writers' call sites** in `src/crop/cropgrowth_helpers.f90` (the `OutCropFixed/OutWofost/OutGrass` calls at task 1/2). Crop data remains available via `result_output.csv` and the `crop` C-API stream.
- **Rewire `src/core/swap_mod.f90`:** replace the `SwapOutput`/`SoilWaterOutput`/`SnowOutput`/`TemperatureOutput`/`SoluteOutput`/`SurfaceWaterOutput` task calls (≈ lines 421–673) with `csv_output_init` (init), `csv_output_step` (per output day, incl. the `flOutputShort` path), `csv_output_finalize` (close).
- `meson.build`: drop `swapoutput.f90`.

### Pre-flight

- Grep `tests/bmi` and `tests/cffi-demo` for any consumer enumerating the four dropped streams; the dropped streams return zeros today, so risk is low, but confirm no test asserts their presence.
- `rm -rf builddir` before the rebuild in this phase (state schema changed).

### Tests

- Full `check-full` (5/5 regression + pFUnit) — not just `check-fast` — because this phase removes wiring and state fields. Per the verify-check-full-before-committing convention for subsystem retirement.
- C-API smoke: `swap_get_output_row('swap_balance')` now returns a non-zero, correctly-sized row (the latent fix); the four dropped streams return `ierr=1` (unknown stream).

## What is preserved

- The two result CSVs and their formats (extended with folded-in crop/snow columns, not changed in shape).
- The `inlist` user-selection feature (the `output_registry` adds a validated, bracket-aware resolver staged for a follow-on; the live selection still flows through `make_userlist`/`det_which_vars`).
- The `dstor`/`baldev`/subregion aggregation math (extracted to the tested `csv_aggregates_mod`).
- The orthogonal C-API accessors actually used by the POC: `swap_get_water_balance` (reads cumulative state directly), `swap_view_array`, `swap_get_scalar`.

## What is retired

- `swapoutput.f90` in full, and with it the `.inc`/`.rot`/`.tem`/`.snw`/`.crp` legacy text writers, the dead `outage` path, and the `.crp` `CropOutput`/`cropoutput` dispatch. (The 3 still-live free-standing procs — `writehead`, `WriteSwapOk` — were relocated verbatim to `src/io/file_headers.f90`; `CloseTempFil` inlined into `swap_main`.)
- **The entire per-subsystem `output_row` C-API mechanism** (per the REVISED decision above): the `swap_get_output_row` function + all 9 stream cases, every `*_output_row`/`*_columns`/`*_n_cols` field across 8 state files, and all `build_*_output_row`/`init|cleanup_*_buffer` helpers. Unused by every consumer; the unified in-memory results accessor is future work on the registry/`csv_output`.

## Out of scope / follow-on (separate arc: IO-IN cleanup)

These input-side items are real but independent; they get their own design when this arc lands:

- Stale `src/io/README.md` (references deleted `readswaptoml.f90`/`readswap.f90`/`macroporeoutput.f90`).
- `config_to_variables.f90` residual "Phase 4f HACK" inline extensions / final adapter collapse.
- `readmeteo.f90` is a misnamed cache *slicer* (zero `read(` statements); rename + clarify its split with `meteo_io.f90`.
- Centralize the meteo CSV column-order knowledge duplicated between `state/atmosphere_state.f90` and `meteo_buffer_mod.f90`.
- (Deferred, decided against for now) typed/string columns in `csv_reader`.

## Risks

- **Untested legacy files:** `.crp`/`.snw` are outside the harness — mitigated by the Phase D parity script gating Phase E's deletions.
- **Format drift:** any change to `fmt_real` could push values past the `1e-2` tolerance — mitigated by lifting the formatter verbatim and the golden-string formatter test in Phase A.
- **C-API consumers:** dropping streams could break a BMI client — mitigated by the Phase E pre-flight grep; dropped streams are already zero-valued.
- **State schema rebuilds:** forgetting `rm -rf builddir` in Phase E yields stale `.mod` and confusing failures — called out explicitly in the phase.
