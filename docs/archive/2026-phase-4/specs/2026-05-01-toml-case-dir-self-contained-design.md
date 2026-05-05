# TOML Case Directory Self-Contained Working Dir — Design

**Date:** 2026-05-01
**Phase:** 4f cleanup, Phase 0 of CSV meteo finalization
**Status:** approved

## Goal

Make `tests/swap-cases/toml/<N>.<case>/` the canonical self-contained working
directory for the modern (TOML-based) SWAP build, and stop the legacy ASCII
case directories from participating in modern-binary runs. Sequence this work
*before* the CSV meteo finalization tasks so subsequent work touches only one
directory tree per case.

## Context

The TOML pathway shipped through Phase 4f introduced a parallel directory tree:

- `tests/swap-cases/<N>.<case>/` — legacy ASCII inputs (`.swp`, `.met`, `.dra`,
  `.YYY`, `.crp`, `.irg`, etc.) plus a few cross-cutting helpers
  (`swap_linux.swp.template`).
- `tests/swap-cases/toml/<N>.<case>/` — TOML inputs (`swap.toml`,
  `swap.dra.toml`, `*.crp.toml`) plus CSV companions (`*.csv`).

Both `run_case.sh` and the regression runner stage modern-pipeline runs by
copying *from* the TOML dir *into* the legacy dir (or a temp dir copy of it),
because the legacy dir is the home of the runtime files SWAP still
needs but the TOML pathway does not yet provide:

- `*.crp` crop files — read by legacy crop sub-readers in `cropgrowth.f90` /
  `irrigation.f90` / `management_soil.f90` (called from the simulation loop).
- `swap_linux.swp.template` → `swap.swp` — read by legacy sub-readers
  (tillage, parts of cropgrowth) that still call `RDinit(swpfile)`.

The current model has two friction points:

1. **Two sources of truth per case.** A modern run depends on files from
   *both* directories; a developer must understand both layouts.
2. **Confused legacy semantics.** `run_case.sh` has a `--legacy` mode that
   tries to run the modern binary against ASCII-only inputs; this mode is
   already documented as "post-Phase-4f generally fails" because the modern
   binary requires `swap.toml`. The mode is dead weight.

The TOML pathway is now the only target for new development. Phase 4f-extend
will eventually port the crop sub-readers, at which point the `.crp` and
`swap_linux.swp.template` dependencies disappear. Until then, those files
need a home — and it should be the TOML directory, not the legacy one.

## Architecture

**Two paths, two binaries, two locations, no crossover at runtime:**

| Path | Binary | Working dir | Inputs |
|---|---|---|---|
| **TOML / modern** (default) | `builddir/swap` | `tests/swap-cases/toml/<N>.<case>/` | `swap.toml` + companions + staged `swap.swp` |
| **Legacy / reference** | user-supplied legacy `swap` | `tests/swap-cases/<N>.<case>/` | `swap.swp` + ASCII files |

After migration, each TOML directory is fully self-contained for the modern
binary. The legacy directories are frozen — read-only inputs for an external
legacy reference binary, never touched by the modern build's tooling.

`run_case.sh` runs SWAP **in-place** in the TOML directory (with selective
cleanup of generated files). The regression runner keeps its temp-dir model
(needed for parallel `ProcessPoolExecutor` isolation) but flips the
`copytree` source from the legacy dir to the TOML dir.

## Components

| Area | Change |
|---|---|
| `tests/swap-cases/toml/<N>.<case>/` (×6) | **Submodule migration:** copy each case's `.crp` files and `swap_linux.swp.template` from `tests/swap-cases/<N>.<case>/` into the matching `toml/<N>.<case>/`. One submodule commit per case. |
| `tests/swap-cases/<N>.<case>/` | **No change.** Frozen as legacy-binary input. |
| `tests/swap-cases/run_case.sh` | **Rewrite.** Drop `--toml` / `--legacy` / `TOML_MODE=auto\|on\|off`. Default: `cd toml/<N>.<case>/`, stage `swap.swp` from the now-local `swap_linux.swp.template`, run, cleanup. New flag `--legacy-binary <path>`: `cd <N>.<case>/`, run that binary, cleanup. Cleanup uses an explicit output-only whitelist. |
| `tests/regression/test_output_regression.py` | `shutil.copytree` source → `tests/swap-cases/toml/<N>.<case>/`. Delete the lines 290–307 block that staged TOML files on top of the legacy dir. Temp-dir + parallel execution model unchanged. |
| `docs/csv-companion-files.md` | Update the `run_case.sh` reference (currently at line 95) to reflect the new working-dir model. |
| `tests/swap-cases/run_case.sh` `--help` output | Document the new flow: TOML dir is the working dir; `--legacy-binary` runs an external binary against the legacy ASCII case dir. |

The cases in scope: `1.hupselbrook`, `2.grassgrowth`, `3.macroporeflow`,
`4.oxygenstress`, `5.salinitystress`, `6.surfacewater`. Macroporeflow is
excluded from regression per ADR 0011 but gets the migration anyway for
consistency (its TOML dir already exists).

## Data flow

```
Developer manual run (TOML, default):
  ./run_case.sh -c hupselbrook -e ../../builddir/swap
  → cd tests/swap-cases/toml/1.hupselbrook
  → cp swap_linux.swp.template swap.swp
  → exec swap                              (cwd = toml/1.hupselbrook)
  → cleanup: rm swap.swp result_*.csv *.log Swap.ok swap.ok ...

Developer manual run (legacy reference binary):
  ./run_case.sh -c hupselbrook --legacy-binary /path/to/legacy/swap420
  → cd tests/swap-cases/1.hupselbrook
  → cp swap_linux.swp.template swap.swp
  → exec /path/to/legacy/swap420           (cwd = 1.hupselbrook)
  → cleanup: rm swap.swp result_*.csv *.log ...

Regression runner (TOML only — legacy not in regression):
  copytree(tests/swap-cases/toml/<N>.<case>/ → tmp/case/)
  cp tmp/case/swap_linux.swp.template tmp/case/swap.swp
  → exec swap                              (cwd = tmp/case)
  → read tmp/case/result_output.csv
  → temp dir auto-cleaned
```

## Cleanup whitelist (in-place TOML run)

`run_case.sh` post-run cleanup deletes only files SWAP generates. The
whitelist (matched in the working dir, not recursive):

- `swap.swp` (staged from template each run)
- `result_*.csv`, `result.*`
- `*.log`, `Swap.ok`, `swap.ok`, `reruns.log`
- `*.tmp`, `*rd$*.tmp`

Files **never** touched (committed inputs):

- `*.toml`, `*.crp.toml` — TOML inputs
- `*.csv` (other than `result_*.csv`) — CSV companions
- `*.crp` — legacy crop files staged into the TOML dir
- `swap_linux.swp.template` — staging source

The cleanup uses explicit globs, not `rm *`. Any new SWAP output file
not on the whitelist will accumulate until the whitelist is updated; this
is preferable to risking input deletion.

## Error handling

- **Missing TOML directory** (`tests/swap-cases/toml/<N>.<case>/` absent or
  has no `swap.toml`): `run_case.sh` aborts with a clear message naming the
  expected path. The regression runner does the same — `RuntimeError` with
  the path.
- **Missing executable**: existing `[ ! -f "$SWAP_EXEC_ABS" ]` check stays.
- **`--legacy-binary` with missing legacy dir or binary**: error and abort.
- **Both `--legacy-binary` and `--exec` set**: error out with a message
  explaining the two flags are mutually exclusive. The intent of each is
  too different to silently pick a winner.

The migration step itself has no runtime error path — it's a one-shot
`git mv`/`cp` + commit operation in the submodule.

## Testing

**Verification after migration:**

```bash
pixi run -e test check-full
```

Expected: 5/5 regression cases green at the existing 1e-2 cm tolerance.
The change is purely organizational; numeric output must not move.

**Manual verification of `run_case.sh`:**

```bash
cd tests/swap-cases
./run_case.sh -c hupselbrook -e ../../builddir/swap
# Expected: case runs to completion; result_*.csv produced briefly then cleaned;
#           swap.toml/csv/crp/template files still present and unchanged.
```

**Manual verification of `--legacy-binary`** (if a legacy binary is
available):

```bash
./run_case.sh -c hupselbrook --legacy-binary /path/to/legacy/swap420
```

Expected: case runs in legacy dir against the supplied binary, no touches
to the TOML dir.

**No new unit tests.** This phase is shell + Python plumbing; behavior is
verified end-to-end via the existing regression suite.

## Out of scope

- **Porting legacy crop sub-readers to TOML.** When that lands (Phase
  4f-extend), the `.crp` files and `swap_linux.swp.template` in the TOML
  dirs become deletable.
- **Adding regression coverage for the legacy binary.** Legacy ASCII cases
  remain manually testable only.
- **Touching the legacy case directories.** They stay exactly as they are.
- **Path-resolution policy refinement** (e.g. `pathwork` vs `pathatm` vs
  `pathcrop` separation). Current convention: all companion files in the
  working directory. Same as before.
- **Modifying `swap.toml` content** for any case. The migration is purely
  a file-location change.

## Risks

| Risk | Mitigation |
|---|---|
| Cleanup whitelist accidentally matches a committed input | Whitelist uses explicit names (`swap.swp`, `result_*.csv`, `*.log`, …); never bare `*.csv` or `*.crp`. |
| Submodule + outer-repo commits get out of sync | Each case is one submodule commit + one outer-repo bump; verify with `git submodule status` after each pair. |
| Macroporeflow migration regresses an existing manual workflow | Macroporeflow is not in regression; its TOML dir already exists. Migration is additive. |
| Developer expects `--toml` / `--legacy` flags from old `run_case.sh` | `run_case.sh --help` clearly documents the new flags; old flags removed cleanly (no silent fallback). |
| `--legacy-binary` with a path containing spaces breaks the script | Quote `"$LEGACY_BINARY"` consistently in the shell script. |

## Sequencing

This work is **Phase 0** of the CSV meteo finalization plan
(`docs/superpowers/plans/2026-05-01-csv-meteo-finalize.md`). It must land
*before* Tasks 1–10 of that plan, so subsequent work touches only the TOML
directory tree per case.
