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
