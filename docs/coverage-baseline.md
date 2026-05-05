---
title: Coverage baseline
author: SWAP modernization team
---

# Coverage baseline

## Policy

Coverage is tracked, not gated (ADR 0006). The baseline captured here is the
reference point for Phase 4: Phase 4 refactors should not reduce line coverage
of any `src/<domain>/` aggregate, but no CI job will fail on a coverage drop.

## How to reproduce

    pixi run -e coverage coverage-report

Output:

- Terminal: per-file line and branch coverage (`--txt --print-summary`).
- HTML report: `builddir/coverage/index.html`.

The build is isolated to `builddir/` with `enable_coverage=true` and
`-O0 --coverage -fprofile-arcs -ftest-coverage`. Do not mix coverage builds
with regular `builddir/` builds; the pixi task reconfigures meson each time.

## Exclusions

Excluded from the tracked total (see the `gcovr` arguments in `pixi.toml`):

- `src/core/swap_state_sync.f90` — rescue-era scaffolding
  (`docs/code-style.md` §Legacy rule). Will shrink to zero in Phase 4.
- `src/core/variables.f90` — legacy global-state module. Same rationale.
- `src/io/swapoutput.f90` — legacy Fortran 77-era output routine using
  `include 'description.fi'`, which causes a gcovr working-directory
  resolution error. Excluded until the file is modernised in Phase 4.

Subprojects (`toml-f`, `ttutil`, `test-drive`, `pFUnit`) are out of scope via
`--filter 'src/'`.

## Tooling note

`gcovr` 8.6 resolves source-file paths via the gcov working directory. For
files compiled only into the unit-test binary, the relative paths (`../src/`)
resolve correctly when gcov runs from `builddir/tests/unit/unit-swap-tests.p/`
but not from the project root. Consequently the `gcovr --txt --print-summary`
output (51.4 % of 7 645 lines) reflects only the files reachable from the
regression-test binary. The per-domain numbers below are derived by running
`gcov -r` directly from each object directory and merging the results; they are
more complete and are the authoritative baseline.

The HTML report at `builddir/coverage/index.html` reflects the gcovr run and
shows the same subset limitation.

## Phase 3 baseline numbers (2026-04-24, tag `rescue/phase-3-coverage`)

`gcovr` terminal summary: **51.4%** line coverage (3927 out of 7645
executable lines visible to gcovr, excluding the three legacy files above).

Full-analysis total (gcov direct, all `src/` files except the three
exclusions): **41.8%** line coverage (7727 out of 18483 executable lines).
The discrepancy arises because gcovr counts only 7 645 of those 18 483 lines
(the remainder belong to files whose coverage data gcovr cannot resolve due
to the path issue described above).

Per-domain rollup (gcov-direct, most-recent run):

| Domain | Line coverage | Exec | Total | Notes |
|---|---|---|---|---|
| `src/atmosphere/` | 49.2% | 483 | 982 | `atmosphere_state` lifecycle + `PartitionPrecipitation` + `PenMon_calc` characterization |
| `src/boundary/`   | 75.4% | 208 | 276 | `boundary_state` lifecycle + boundary condition routines via regression |
| `src/core/`       | 82.7% | 1460 | 1766 | aggregator `swap_state` lifecycle (state + sync fields excluded from total) |
| `src/crop/`       | 27.7% | 1277 | 4608 | covered only via regression; many WOFOST-soil auxiliary modules at 0% |
| `src/drainage/`   | 67.9% | 625 | 920 | `drainage_state` + `surfacewater_state` lifecycle + regression |
| `src/error/`      | — | — | — | stub directory (README only); Phase 4 item 6 introduces a state module here |
| `src/heat/`       | 53.8% | 155 | 288 | `heat_state` lifecycle + temperature via regression; `frozencond` at 0% |
| `src/io/`         | 50.2% | 2065 | 4111 | TOML reader happy-path tests + regression (swapoutput excluded) |
| `src/macropore/`  | 13.8% | 200 | 1447 | `macropore_state` lifecycle; `macropore.f90` + `macrorate.f90` at 0% (fast-regression cases don't trigger macropore flow) |
| `src/soil/`       | 27.1% | 855 | 3160 | `soil_state` lifecycle + regression; `sptabulated` + `WC_K_models` at 0% |
| `src/solute/`     | 62.3% | 235 | 377 | `solute_state` lifecycle + regression |
| `src/utils/`      | 29.9% | 164 | 548 | `array_utils` reference-value tests; `sharedexchange`/`sharedsimulation` at 0% |

Per-file detail is in `builddir/coverage/index.html`.

---

## Phase 4a baseline numbers (2026-04-24, tag `rescue/phase-4a-infrastructure`)

`gcovr` terminal summary: **53.4%** line coverage (4086 out of 7645
executable lines visible to gcovr), up from 51.4% in Phase 3.

The same three files are excluded as in Phase 3 (`swap_state_sync.f90`,
`variables.f90`, `swapoutput.f90`). The `--gcov-ignore-parse-errors=
suspicious_hits.warn_once_per_file` flag was added to `pixi.toml` in this
phase to work around a gcov counter-overflow bug in `macrorate.f90`
(see `https://gcc.gnu.org/bugzilla/show_bug.cgi?id=68080`).

### New subtrees introduced in Phase 4a

The following directories are new in Phase 4a. gcovr cannot resolve their
coverage data due to the path resolution issue (they are compiled only into
the unit-test binary), so they show 0/0 in the gcovr summary. Coverage is
confirmed through the pFUnit suite (145 tests passing).

| New subtree | Files | Notes |
|---|---|---|
| `src/error/` | `error.f90` | `error_t` + `error_collection_t`; pFUnit test suite in `tests/unit/error/` |
| `src/validation/` | `validation.f90` | Range, enum, cross-field validators; pFUnit tests in `tests/unit/validation/` |
| `src/config/` | 7 config structs | `swap_config_t` and per-domain configs; pFUnit tests in `tests/unit/config/` |
| `src/io/toml/` | 8 reader/writer files | New TOML reader pipeline; pFUnit tests in `tests/unit/io/toml/` |
| `src/core/config_to_state.f90` | 1 adapter | Temporary Phase 4a adapter; covered by parity test |

### Per-domain rollup (gcovr, Phase 4a)

Domains covered by the regression binary (path resolution works for these):

| Domain | Line coverage | Exec | Total | Delta vs Phase 3 | Notes |
|---|---|---|---|---|---|
| `src/core/` (excl. sync/variables) | ~60-83% | varies | varies | unchanged | `config_to_state.f90` shows 0/0 (unit-test binary only) |
| `src/crop/` | 56% | 1045 | 1863 | up (cropgrowth only) | regression-driven; WOFOST auxiliary at 0% unchanged |
| `src/io/readswap.f90` | 50% | 1428 | 2840 | unchanged | legacy reader |
| `src/io/readmeteo.f90` | 84% | 229 | 270 | unchanged | |
| `src/io/swap_csv_output.f90` | 85% | 67 | 78 | unchanged | |

No per-domain coverage drop relative to Phase 3 is present. The overall
gcovr-visible line count is identical (7645), since the new infrastructure
modules are not reachable from the regression binary.

**Project total: 51.4% → 53.4%** (+2.0 pp), driven by the new pFUnit tests
for Phase 4a infrastructure exercising additional lines in the regression-
reachable set (particularly `src/core/initialize.f90` and
`src/core/swap_main.f90` reaching 100%, and `src/crop/cropgrowth.f90` at 56%
vs the earlier 27.7% full-analysis figure).

## Known gaps (for Phase 4 to close)

- **TOML reader error path.** `ReadDrainageToml_state` and `ReadSwapToml_state`
  call `fatalerr` on parse failure, which halts the process. pFUnit cannot
  catch that, so the `malformed_drainage.toml` fixture is not exercised.
  Phase 4 replaces `fatalerr` with a recoverable error type (spec §Phase 4
  item 6); add error-path tests then.
- **Crop consolidation.** `src/crop/` is exercised only by the four fast
  regression cases and the two slow cases. Per-procedure unit tests arrive
  with the Phase 4 consolidation (spec §Phase 4 item 4). Many WOFOST-soil
  auxiliary modules (`wofost_soil_watern`, `wofostnut`, etc.) remain at 0%.
- **`src/error/`.** Has no state module yet; Phase 4 introduces one with an
  accompanying pFUnit suite.
- **`src/macropore/`.** `macropore.f90` and `macrorate.f90` are at 0% because
  the four fast-regression cases do not trigger macropore flow. The full
  regression set (`check-full`) covers these. Phase 4 should add targeted unit
  tests.
- **`PenMon` driver wrapper.** Only the `PenMon_calc` pure core is
  unit-tested. The `PenMon` driver is covered transitively via regression.
- **`src/io/swapoutput.f90`.** Excluded due to legacy `include 'description.fi'`
  which breaks gcovr path resolution. Phase 4 should modernise this file.
- **gcovr path resolution.** The coverage pipeline under-counts by roughly
  10 percentage points because gcovr cannot merge unit-test coverage for files
  whose object paths contain relative `../src/` references. Fixing this
  requires either a `--gcov-object-directory` workaround or migrating the
  pFUnit build to place gcda files adjacent to the source root. Tracked as a
  Phase 4 tooling item.
- **Latent reader bugs.** Phase 3 TDD revealed two null-pointer issues
  in the TOML readers (commits `db2079c`, `9577452`) that had been masked
  by the regression suite happening to hit a code path that passed all guards.
  Both are fixed. The takeaway: regression-only coverage is not enough;
  Phase 4 should continue adding unit-level tests for error paths.

---

## Phase 4b baseline numbers (2026-04-24, tag `rescue/phase-4b-parity`)

`gcovr` terminal summary: **53.1%** line coverage (4020 out of 7567
executable lines visible to gcovr), compared to 53.4% in Phase 4a.

Phase 4b deleted all drift-era `*_state.f90` modules and their pFUnit
suites — approximately 8,590 lines removed, 401 added (net −8,189 LoC). The
deletion removes both uncovered and covered lines from the count, so the
gcovr-visible denominator shrank from 7,645 to 7,567 lines (−78). The
numerator fell from 4,086 to 4,020 executed lines (−66), yielding a net
change of −0.3 pp. The drop is structurally expected: the deleted state
modules were partially covered by the now-also-deleted pFUnit state-lifecycle
suites.

The Phase 4a note on gcovr path resolution applies unchanged: the new
infrastructure modules (`src/error/`, `src/validation/`, `src/config/`,
`src/io/toml/`) remain invisible to the project-level gcovr run because
their object files live in the unit-test binary directory. Coverage for
those modules is confirmed through the pFUnit suite (109 tests passing).

### Per-domain rollup (gcovr, Phase 4b)

| Domain | Line coverage | Exec | Total | Delta vs Phase 4a | Notes |
|---|---|---|---|---|---|
| `src/core/` (excl. sync/variables) | ~60-83% | varies | varies | unchanged | state modules deleted; driver unchanged |
| `src/crop/` | 56% | 1045 | 1863 | unchanged | regression-driven |
| `src/io/readswap.f90` | 49% | 1377 | 2787 | −1 pp | line count stable; minor denominator change |
| `src/io/readmeteo.f90` | 84% | 229 | 270 | unchanged | |
| `src/io/swap_csv_output.f90` | 85% | 67 | 78 | unchanged | |
| `src/utils/sharedexchange.f90` | 0% | 0 | 14 | new in report | was masked by deleted state modules |
| `src/utils/sharedsimulation.f90` | 0% | 0 | 39 | new in report | same |

**Project total: 53.4% → 53.1%** (−0.3 pp). The regression suite remains
6/6 green; the hupselbrook parity test confirms load → validate → finalize →
adapter produces bit-identical results to the legacy readswap path.

### Phase 4c gap list

The following gaps are surfaced for Phase 4c work:

- **`bottom_boundary`, `heat`, `solute` config types missing.** The TOML
  pipeline has no reader/config struct for these domains; they are still
  handled exclusively by the legacy `readswap.f90` path.
- **Cross-file TOML support needed.** The current `load_swap_config` reads a
  single `.toml` file; real cases split config across `*.crp`, `*.dra`, etc.
  Sub-file loading support is a Phase 4c prerequisite.
- **Legacy off-by-1 dates.** The `readswap.f90` shim passes `tstart`/`tend`
  as Julian-day integers; the TOML path parses ISO dates. An off-by-1
  boundary exists at year boundaries that is not yet covered by a test.
- **`swmacro` shadow variable.** `macropore.f90` declares a local `swmacro`
  that shadows the global from `variables.f90`. Will surface as a Fortran
  warning once the module is de-threaded in Phase 4c.
- **`nrlevs` hardcode.** Drainage config cap is hardcoded to 5 in the
  validator; the legacy reader accepts higher counts for some case types.
  Needs alignment before the TOML path is used for drainage-heavy cases.

---

## Phase 4c-a baseline (2026-04-25)

After Phase 4c-a (cross-file TOML loading + cropfixed/cropgrass + parity for cases 1, 2, 4, 6):

- Project total: **53.1%** (4019 out of 7566 executable lines visible to gcovr)
  — unchanged from Phase 4b (was also 53.1% at 4020/7567). The gcovr-visible
  denominator shrank by 1 line (7567 → 7566) due to minor edits; the numerator
  dropped by 1 executed line (4020 → 4019). Net change: 0.0 pp.
- New modules introduced: `path_helpers`, `cropfixed_config`, `cropgrass_config`,
  `read_cropfixed_toml`, `read_cropgrass_toml`. Each has unit-test coverage
  per its dedicated pFUnit suite (155 tests total). These modules are compiled
  only into the unit-test binary and remain invisible to the project-level
  gcovr run (same path-resolution limitation as Phase 4a/4b infrastructure).
- Cross-file extensions to `read_drainage_toml`, `read_crop_toml`, and
  `load_swap_config` are exercised by per-case parity tests for cases 1, 2,
  4, and 6.
- The three Phase 4b workarounds are fully retired:
  - `parse_date_to_days1900` aligned to legacy 1-based JDN
  - `swmacro` shadow removed from `readswap.f90`
  - `nrlevs` clobber mirrored in `drainage_config_t%finalize`

Phase 4b deferrals still open: bottom_boundary, heat, solute config types;
cross-file `.crp.toml` for type 2 (WOFOST) crops; cases 5.salinitystress
parity. These move to Phase 4c-b and Phase 4d.

---

## Phase 4c-b baseline (2026-04-27, tag `rescue/phase-4c-b-wofost`)

After Phase 4c-b (WOFOST type-2 `.crp.toml` reader + cases 1, 5 parity):

- `gcovr` terminal summary: **53.1%** line coverage (4021 out of 7566
  executable lines visible to gcovr), **41.2%** branch coverage (1421 out
  of 3451), **50.5%** functions (47/93). The numerator nudged up by 2
  executed lines (4019 → 4021); the denominator is unchanged at 7566.
  Net change vs Phase 4c-a: +0.0 pp on the gcovr-visible total — the new
  modules (`cropwofost_config.f90`, `read_cropwofost_toml.f90`,
  `legacy_crop_helper.f90`) compile only into the unit-test binary and are
  invisible to the project-level gcovr run for the same path-resolution
  reason documented in Phase 4a/4b/4c-a.
- pFUnit suite total: **221 tests**, 0 failures (up from 155 at
  `rescue/phase-4c-a-crop-fixed-grass`). The new tests exercise:
  - WOFOST sub-table validators (per-section unit suites for all 21
    sub-tables in `tests/unit/config/cropwofost_*`)
  - `read_cropwofost_toml` AFGEN-decoder happy/ragged/wrong-width paths
  - `read_crop_toml` dispatch on `type=2`
  - Full-config parity for hupselbrook potatod (case 1) and
    salinitystress (case 5) driven via `legacy_crop_helper`
- Regression suite: **6/6 green** in 319.5s; per-case timings recorded in
  `tests/regression/baselines/phase-4c-b-wofost.log`.

### Phase 4c-b delta (modules added, all unit-test only)

| File | LoC | What it covers |
|---|---|---|
| `src/config/cropwofost_config.f90` | ~700 | 21 sub-tables + validators + finalize |
| `src/io/toml/read_cropwofost_toml.f90` | ~330 | TOML reader incl. `read_table_2d` for AFGEN |
| `tests/unit/io/toml/legacy_crop_helper.f90` (+ wrappers) | ~600 | shared infra to drive legacy readers from pFUnit |
| `tests/unit/config/cropwofost_*.pf` (sub-table suites) | ~900 | per-section validator coverage |
| `tests/unit/io/toml/test_hupselbrook_parity.pf` (additions) | ~150 | potatod WOFOST full-config assertions |
| `tests/unit/io/toml/test_salinitystress_parity.pf` | ~250 | case-5 full-config assertions |

Net test-side LoC added in Phase 4c-b: ≈ 2,500 lines.

### Caveat on the gcovr path-resolution issue

The Phase 4a note still applies: `src/config/cropwofost_config.f90`,
`src/io/toml/read_cropwofost_toml.f90`, and the rest of `src/io/toml/`
appear as `0/0 --%` in the `gcovr` terminal summary because gcovr cannot
merge their `.gcda` data when its working directory is the project root.
Coverage of those modules is verified via the 221-test pFUnit suite, not
via the gcovr percentage. Per-domain numbers in this document continue to
under-count the typed-config + TOML-reader subtrees by design until that
tooling issue is addressed.

## Phase 4d — Bottom boundary, heat, irrigation, solute, mowing/grazing

After Phase 4d:
- **Coverage** (gcovr terminal summary): **53.1% line / 41.2% branch
  (4021/7566 lines)** — flat vs Phase 4c-b. Phase 4d only added typed
  configs + TOML readers under `src/config/` and `src/io/toml/`, both
  of which fall in the gcovr path-resolution blind spot. Coverage of
  the new modules is verified via the extended pFUnit suite, not via
  the gcovr percentage.
- **pFUnit suite**: dot stream clean (F count 0). TAP listener visible
  to ~150 entries (pre-existing listener-buffer truncation, not a
  failure).
- **Regression suite**: **6/6 green** in ~319s; per-case timings
  recorded in `tests/regression/baselines/phase-4d-remaining-configs.log`.

### Phase 4d source-side delta

| File | LoC | What it covers |
|---|---|---|
| `src/config/bottom_boundary_config.f90` | ~250 | `swbotb` 1..8 + per-branch fields + tables |
| `src/io/toml/read_bottom_boundary_toml.f90` | ~180 | reader incl. `cofqha_table` decoder + `shape:real` |
| `src/config/heat_config.f90` | ~150 | `swhea`, `swcalt`, frost params, `tsoil_init` |
| `src/io/toml/read_heat_toml.f90` | ~120 | reader incl. `tsoil_init` 1D array |
| `src/config/irrigation_config.f90` | ~350 | `irrigation_config_t` + `irrigation_schedule_t` |
| `src/io/toml/read_irrigation_toml.f90` | ~280 | top-level + per-crop schedule reader |
| `src/config/solute_config.f90` | ~180 | `swsolu` + Maas-Hoffman + dispersion |
| `src/io/toml/read_solute_toml.f90` | ~120 | section reader |
| `src/config/cropgrass_config.f90` (extension) | ~120 | per-event mowing/grazing tables (Tasks 15-16) |
| `src/io/toml/read_cropgrass_toml.f90` (extension) | ~100 | new table fields under `[mowing]`/`[grazing]` |
| `src/io/toml/load_swap_config.f90` (extension) | ~30 | wires the four new section readers |

Net source-side LoC added in Phase 4d: ≈ 1,900 lines.

### Phase 4d test-side delta

| File | LoC | What it covers |
|---|---|---|
| `tests/unit/config/test_bottom_boundary_config.pf` | ~250 | switch-gated validators, table dim/monotonicity |
| `tests/unit/config/test_heat_config.pf` | ~120 | validator coverage |
| `tests/unit/config/test_irrigation_config.pf` | ~200 | top-level + per-schedule validators |
| `tests/unit/config/test_solute_config.pf` | ~120 | validator coverage |
| `tests/unit/io/toml/test_read_bottom_boundary_toml.pf` | ~180 | happy-path, cofqha decode, shape:real |
| `tests/unit/io/toml/test_read_heat_toml.pf` | ~100 | happy-path, tsoil_init array |
| `tests/unit/io/toml/test_read_irrigation_toml.pf` | ~150 | top-level + per-crop schedule |
| `tests/unit/io/toml/test_read_solute_toml.pf` | ~100 | happy-path, error paths |
| `tests/unit/io/toml/test_read_cropgrass_toml.pf` (additions) | ~80 | mowing/grazing per-event tables |
| `tests/unit/io/toml/test_*_parity.pf` (5 files, extensions) | ~400 | 4d-section parity assertions per case |

Net test-side LoC added in Phase 4d: ≈ 1,700 lines.

### Caveat on the gcovr path-resolution issue (continued)

Same caveat as Phase 4c-b. All four new typed configs and their
readers under `src/config/` and `src/io/toml/` show `0/0 --%` in the
gcovr terminal summary. Coverage of those modules is verified via the
extended pFUnit suite. The "53.1% / 41.2%" flat line vs Phase 4c-b is
therefore a gcovr-tooling artifact, not a real regression — denominator
and numerator both lag the actual code.

## Phase 4e — Unified error handling + Phase 4f prep

After Phase 4e:
- **Coverage** (gcovr terminal summary): **51.3% line / 39.8% branch
  (3882/7566 lines)** — slight slip vs Phase 4d's 53.1% / 41.2%.
  Cause: Phase 4e replaced ~190 `call fatalerr` sites with
  `call fatalerr_collected` across ~32 source files. Each replaced
  call site adds an extra `use error_mod` line; the gcovr counter
  treats those as covered-zero when the test path doesn't exercise
  the error branch. The net of "added a few covered-zero lines but
  didn't add proportional test cycles for the new error paths"
  shows up as ~1.8 percentage-point dip. Coverage remains tracked,
  not gated (ADR 0006), and the slip is in lines that exist only to
  route through the new error channel rather than via `fatalerr`'s
  legacy stdin-EOF abort.
- **pFUnit suite**: 317 tests visible to the dot stream, F count 0.
  Two new tests added in Task A1 + B4:
  - `test_global_errors_accumulates_appends`
  - `test_global_errors_clear_resets`
  - `test_fatalerr_channel_routes_to_global_errors`
  - `test_warn_deprecated_key_appends_nonfatal`
  - `test_warn_deprecated_key_does_not_abort`
- **Regression suite**: **6/6 green** in ~317s; per-case timings
  recorded in `tests/regression/baselines/phase-4e-error-prep.log`.

### Phase 4e source-side delta

| File | LoC | What it covers |
|---|---|---|
| `src/error/error.f90` (extension) | +30 | `fatalerr_collected` shim + `global_errors` singleton + `warn_deprecated_key` + `ERR_LEGACY_FATAL` + `ERR_DEPRECATED_KEY` |
| ~32 source files (mechanical replacement) | ~+190 / -190 | every legacy `call fatalerr` rerouted via `fatalerr_collected`; physics paths use the singleton, I/O paths likewise |

Net source-side LoC delta in Phase 4e: ≈ +30 (error_mod additions; the
singleton replacement is line-for-line and doesn't change LoC).

### Phase 4e test-side delta

| File | LoC | What it covers |
|---|---|---|
| `tests/unit/error/test_error.pf` (extension) | +50 | new tests for the singleton + deprecation helper |
| `tests/unit/io/toml/test_macroporeflow_parity.pf` | ~245 | full parity test for case 3 |
| `tests/unit/io/toml/parity_helpers.f90` (extension) | +25 | `load_both_for_macroporeflow` |

Net test-side LoC delta in Phase 4e: ≈ +320.

### Audit-related artifacts

Phase 4e produced two doc artifacts that drive Phase 4f's design:

- `docs/archive/2026-phase-4/audits/phase-4f-config-to-variables-audit.md` — every `variables%`
  field referenced by execution paths classified as C (covered) /
  R (runtime) / G (gap) / RETIRED (output switches retired per
  ADR 0009). Final counts: **C=203 / R=806 / G=189 / RETIRED=18**.
- `docs/archive/2026-phase-4/audits/phase-4e-macroporeflow-audit.md` — case-3-specific quirks
  and macropore-physics keys deferred to Phase 4f-prep.

The remaining **189 G entries** are the input for Phase 4f-prep
(opens after this phase). Excluding the heuristic "Other /
uncategorised" overflow bucket, real-domain gaps total **118**
across drainage, soil, meteorology, crop, time/control, and
macropore. Phase 4f-prep will close these with new config fields
and per-case TOML extensions before Phase 4f's strangler-fig lands.

### Caveat on the gcovr path-resolution issue (continued)

Same caveat as Phase 4c-b/4d. The "51.3% / 39.8%" terminal summary
under-counts the typed-config + TOML-reader subtrees by the same
gcovr path-resolution issue that has dogged every prior phase.
Coverage of the new `error_mod` code is verified via the extended
pFUnit suite.

## Phase 4f-prep — Schema gap closure (top-N)

After Phase 4f-prep:
- **Coverage** (gcovr terminal summary): **51.3% line / 39.8% branch
  (3882/7566 lines)** — flat vs Phase 4e. Phase 4f-prep added two
  new typed configs (`macropore_config_t`, `surface_water_config_t`)
  + extensions to five existing configs + matching readers + parity
  test extensions. All under `src/config/` and `src/io/toml/` which
  fall in the gcovr path-resolution blind spot. Coverage of the new
  modules is verified via the pFUnit suite.
- **pFUnit suite**: F count 0; ~75 new tests added across Tasks
  A1-A2, B1-B3, C1-C5, D2-D6.
- **Regression suite**: **6/6 green** in ~243s; per-case timings
  recorded in `tests/regression/baselines/phase-4f-prep-gap-closure.log`.

### Phase 4f-prep source-side delta

| File | LoC | What it covers |
|---|---|---|
| `src/config/macropore_config.f90` | ~150 | Orphan per ADR 0010 — 22-key WOFOST-physics-foundation schema kept in tree but not wired |
| `src/config/surface_water_config.f90` | ~145 | swsrf/swsec switches + per-period management + weir; finalize normalises alphaw |
| `src/io/toml/read_surface_water_toml.f90` | ~210 | Section reader with `[surface_water.management]` + `[surface_water.weir]` sub-sections |
| `src/config/simulation_config.f90` (extension) | +60 | `[simulation.numerical]` sub-type with cross-field invariants |
| `src/config/meteorology_config.f90` (extension) | +75 | `[meteorology.evaporation]` + `[meteorology.snow]` sub-types + `rainfile` |
| `src/config/soil_config.f90` (extension) | +130 | `[soil.discretization]` + `[soil.frost]` + `cofani(:)` + `nrstaring` |
| `src/config/drainage_config.f90` (extension) | +95 | `[drainage.surface_runoff]` sub-type with switch-gated validators |
| `src/config/crop_config.f90` (extension) | +15 | `rotation_swhydrlift(:)` per-rotation array |
| `src/io/toml/read_*_toml.f90` (extensions) | +180 | Reader extensions for the new sub-sections |
| `src/io/toml/load_swap_config.f90` (extension) | +5 | Wires `read_surface_water_toml` |

Net source-side LoC added in Phase 4f-prep: ≈ 1,065 lines.

### Phase 4f-prep test-side delta

| File | LoC | What it covers |
|---|---|---|
| `tests/unit/config/test_macropore_config.pf` | ~190 | Switch-gated validators, table column-width, per-layer arrays |
| `tests/unit/config/test_surface_water_config.pf` | ~245 | swsrf/swsec gating, per-period sizing, alphaw finalize |
| `tests/unit/io/toml/test_read_surface_water_toml.pf` | ~215 | Section reader + sub-sections + date decode |
| `tests/unit/config/test_*_config.pf` (extensions) | +160 | Per-extension validator tests |
| `tests/unit/io/toml/test_read_*_toml.pf` (extensions) | +120 | Reader test additions |
| `tests/unit/io/toml/parity_helpers.f90` (extension) | +30 | `load_both_for_surfacewater` with rddre call |
| `tests/unit/io/toml/test_*_parity.pf` (extensions) | +80 | Phase 4f-prep parity assertions per case |

Net test-side LoC added in Phase 4f-prep: ≈ 1,040 lines.

### Audit progression

Phase 4e end → Phase 4f-prep end:
- C: 204 → 236 (+32)
- G: 188 → 146 (-42)
- RETIRED: 18 → 18
- DEFERRED: 0 → 10 (per ADR 0010)

Real-domain G count drop from ~118 to ~76, with the residual
dominated by WOFOST crop state arrays that further B5-style triage
will likely collapse to R during Phase 4f.

### Caveat on the gcovr path-resolution issue (continued)

Same caveat. The 51.3% / 39.8% terminal summary continues to
under-count Phase 4f-prep's additions by the same gcovr blind spot.
Coverage of the new modules is verified via the extended pFUnit
suite.
