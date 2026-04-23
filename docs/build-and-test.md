---
title: Build and test
author: SWAP modernization team
---

# Build and test

## Quick start

One command builds the binary and runs the fast test set:

    pixi run -e test check-fast

That task configures meson (gfortran, single `builddir/`), compiles the `swap`
binary, runs the pFUnit suite (empty at the rescue baseline — see the pFUnit
section), and runs four fast regression cases: `hupselbrook`, `surfacewater`,
`salinitystress`, `grassgrowth`. Target budget: under 90 seconds on typical
hardware. Run it before every commit; anything that doesn't pass `check-fast`
does not land on `development`.

## Prerequisites

Pixi manages everything. The `[dependencies]` and `[feature.test.dependencies]`
blocks of `pixi.toml` pin Python 3.11, gfortran, meson, ninja, fprettify, the
pytest family, pandas, pyswap, and FORD into `.pixi/`. No system-level install
is required beyond pixi itself.

Install pixi from https://pixi.sh and run any `pixi run ...` command in this
repository — pixi will resolve and install the environment on first use. For
the test tasks you want the `test` feature, which is selected by
`pixi run -e test <task>`.

## Build

`pixi run build-linux` produces `builddir/swap`. It depends on `_configure`,
which runs `FC=gfortran meson setup builddir --reconfigure` so meson picks up
the gfortran toolchain regardless of what is on `PATH`. A single `builddir/` is
used for both the production binary and the pFUnit unit-test executables
(see ADR 0002 — single builddir).

Re-running `build-linux` after editing Fortran sources is cheap; meson only
reconfigures when build options change. If meson options genuinely need to
change (e.g., toggling `enable_pfunit`), invoke `pixi run _configure`
explicitly or delete `builddir/` first.

`pixi run clean` removes all known build trees: `builddir/`, `builddir_windows/`,
`builddir_windows_native/`, and `builddir_gfortran/` (the last is a legacy name
from before ADR 0002; cleaned defensively in case a stale copy exists).

Cross-compilation for Windows is available via
`pixi run build-windows-cross` — see the cross-compilation section below.

## Fast vs full test protocol

Two canonical gates govern every iteration and every phase tag:

- `pixi run -e test check-fast` — build + pFUnit + four fast regression cases
  (`hupselbrook`, `surfacewater`, `salinitystress`, `grassgrowth`). Target
  budget: under 90 seconds. **Run on every commit.**
- `pixi run -e test check-full` — build + pFUnit + all six regression cases
  (`check-fast` set plus `oxygenstress` and `macropore`). Target budget: about
  ten minutes. **Run before every phase tag and before any fast-forward to
  `main`.**

Ad-hoc subsets for debugging a single case:

    pixi run -e test regression hupselbrook
    pixi run -e test regression macropore oxygenstress

Fixture regeneration (only after a deliberate physics change):

    pixi run -e test regression --regenerate-fixtures

That rewrites every `*_expected_gfortran.json` from the current binary's
output. Regeneration is a load-bearing act — see the fixture policy section
and the contributing notes before using it.

## Regression harness internals

The harness lives in `tests/regression/test_output_regression.py`. For each
case it:

1. Creates an isolated temporary directory.
2. Copies the case's input files in (from `tests/swap-cases/<case_dir>/`),
   honoring the `input_files` map that renames templates like
   `swap_linux.swp.template` → `swap.swp`.
3. Runs `builddir/swap` inside that temp dir.
4. Parses `result_output.csv` and aggregates per year:
   - `flux_vars` — annual sum (e.g., `RAIN`, `DRAINAGE`, `TACT`);
   - `state_vars` — annual mean (e.g., `GWL`);
   - `cumul_vars` — last value per year (for outputs already cumulative in
     the CSV).
5. Compares the aggregates against the stored fixture with tolerance
   `TOL = 1e-2` (cm, per variable, per year).

Cases are declared at the top of the file in the `CASES` dict, each a
`CaseConfig(name, case_dir, fixture, flux_vars, state_vars, cumul_vars,
input_files)`. Adding a new case means adding a new `CASES` entry and a
fixture file.

Cases run in parallel through a `ProcessPoolExecutor`, one worker per CPU.
Because every case runs in its own temp dir, a failure in one case cannot
corrupt the inputs of another. On mismatch the harness prints a table of
offending `(year, variable, observed, expected, delta)` rows and exits
non-zero.

## Fixture policy

Two fixture sets coexist in `tests/regression/`:

- `*_expected.json` — historical reference produced by the **Intel ifx**
  build of the pre-rescue codebase. Preserved in git; **not** compared against
  at runtime. Keep them for diffing if you want to see how the current
  gfortran build differs from the legacy Intel result.
- `*_expected_gfortran.json` — the **live truth**. The harness compares
  against these. Produced from the gfortran + `-finit-local-zero` build at
  rescue Phase 1 via `--regenerate-fixtures`.

Known, accepted tolerances baked into the `_gfortran` fixtures:

- `oxygenstress`: `MOWDM` deviates by up to 85 in year 1995 — pre-existing
  under both compilers.
- `macropore`: `DRAINAGE` deviates by up to 21 at year 1998 under gfortran
  versus ifx.

See `tests/regression/INVESTIGATION_NOTES.md` for the open
compiler-versus-physics question and the reproduction procedure using the
archived ifx reference at `tests/reference/swap420`.

Do not delete the `*_expected.json` files. Do not regenerate
`*_expected_gfortran.json` casually — the commit that does must explain
the physics change that justifies it.

## `tests/swap-cases` submodule

`tests/swap-cases/` is a git submodule pinned to branch `main`. It carries the
canonical case input directories consumed by the regression harness.

After every regression run, `git status` on the outer repo shows `tests/swap-cases`
as ` m` (dirty working tree). That is because the harness copies inputs into
temp dirs by renaming certain template files (e.g.,
`swp.toml` → `swap.toml`) — the rename happens inside the submodule working
tree as part of staging but the **outer-repo pointer to the submodule commit
is unchanged**.

To reset the submodule working tree to its pinned commit:

    git submodule update tests/swap-cases

To bump the pinned commit intentionally:

    cd tests/swap-cases
    git pull origin main
    cd ../..
    git add tests/swap-cases
    git commit -m "chore(tests): bump swap-cases pointer"

## pFUnit

`tests/pFUnit/` is a local checkout and CMake-built install of the Goddard
pFUnit 4.15 release. The install tree at
`tests/pFUnit/build/install_gfortran/PFUNIT-4.15` is what meson consumes — its
path is surfaced via the `PFUNIT_ROOT` environment variable set in pixi's
`[activation]` block.

**The tracking is weird and worth knowing about.** `tests/pFUnit/` is a
**gitlink** (git tree entry mode `160000`) pinned at commit
`581cd32be937df1922ce50fb06ddbc5522e6ce66` in the outer repo, but it is
**not** declared in `.gitmodules`. Consequences:

- `git clone` of the outer repo does **not** automatically populate
  `tests/pFUnit/`.
- `git submodule update` has no entry to act on.
- `git add tests/pFUnit/<anything>` fails with "Pathspec is in submodule".
- The working-tree contents are whatever was checked out when pFUnit was
  first dropped in.

This is a rescue-era inheritance, not a design decision. Converting it to one
of: (a) a properly-declared submodule, (b) a meson subproject via `.wrap`, or
(c) plain in-tree contents tracked by the outer repo — is deferred to Phase
1.5 or Phase 2 of the rescue. The install works on developer machines and
nothing active depends on a fresh populate.

**Rebuilding the install from the existing local clone:**

    cd tests/pFUnit
    mkdir -p build && cd build
    cmake -DCMAKE_INSTALL_PREFIX=install_gfortran ..
    make -j
    make install

**From a fresh clone** (if `tests/pFUnit/` is empty on your machine): clone
https://github.com/Goddard-Fortran-Ecosystem/pFUnit at tag `v4.15` into
`tests/pFUnit/`, then rebuild as above.

The meson config at `tests/unit/meson.build` expects the install at
`tests/pFUnit/build/install_gfortran/PFUNIT-4.15` and errors with a pointer
back to this doc if the path is missing.

## Cross-compilation

`pixi run build-windows-cross` builds for Windows via mingw, using
`cross_mingw.txt`. The result is `builddir_windows/swap.exe`, runnable on
Windows or under wine on Linux. It is best-effort, not part of the
`check-fast` / `check-full` gates, and not exercised in CI.

## Cleaning

`pixi run clean` removes `builddir/`, `builddir_windows/`,
`builddir_windows_native/`, and `builddir_gfortran/`. Use it after pulling
changes that touch `meson.build` or `subprojects/`, or whenever a build
turns flaky for unexplained reasons.

## Common problems

- **`Unsupported Fortran compiler: intel`** — your shell has Intel oneAPI
  sourced, so meson picks up `ifx` and rejects it per ADR 0001. Fix: use
  `pixi run build-linux`, which sets `FC=gfortran` explicitly. Or, for
  ad-hoc meson invocations, run in a shell without oneAPI sourced.

- **`pFUnit not found`** —
  `tests/pFUnit/build/install_gfortran/PFUNIT-4.15` is missing. Either copy
  an install from another clone of the repo, or rebuild it per the pFUnit
  section above.

- **`No suitable tests defined.`** on `pixi run test-pfunit` — expected at
  the rescue baseline. After Phase 1's test reconciliation (step D2 of the
  Phase 1 plan), `pf_files` is empty; pFUnit compiles a test executable but
  runs zero suites. Phase 4 re-populates the pFUnit suites as state modules
  come back. This is not an error.

- **Tolerance exceeded on `macropore` DRAINAGE** — you may have inadvertently
  built with Intel instead of gfortran; `macropore` output is
  compiler-sensitive (see `tests/regression/INVESTIGATION_NOTES.md`). Clean
  and rerun with `pixi run build-linux`.
