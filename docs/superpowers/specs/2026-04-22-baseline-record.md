---
title: Baseline record at rescue/phase-0-baseline
date: 2026-04-22
commit: e256bc0a43e349ad3692f6c97fa98d06e07943d8
status: accepted
---

# Phase 0 Baseline Record

This file locks in the observed state of the `e256bc0` baseline at the moment the Rescue & Stabilize workflow started. Future phases verify against these numbers and flag regressions. Several findings here contradict assumptions in the original spec (`docs/superpowers/specs/2026-04-22-rescue-and-stabilize-design.md`); the spec will be amended separately. The truth observed on the machine, not the assumption in the spec, is what governs Phase 1 onward.

## Git topology

- `main` = `development` = `e256bc0` (local only).
- `origin/main` untouched at `7587ca3` (drift-era tip).
- Archives preserved: `archive/wip-drifted` (`cfb4aec`), `archive/main-pre-rescue` (`7587ca3`), `archive/swaplib` (`c5e35ca`), `archive/swaplib-simple` (`ddc8565`).
- Orphan branch `legacy/swap-4.2.0` to be created in a subsequent dispatch (D6).
- Submodule `tests/swap-cases` pinned at `b1bc902af7090b2e04e4b73557e8e882549e0d56` ("initial commit for TOML based SWAP config"). After a regression run, its working tree shows ` m` because `run_case.sh` renames inputs (e.g. `swp.toml` → `swap.toml`); the outer-repo pointer is unchanged.

## Compile

- `pixi run build-linux` succeeds.
- **Compiler observed: `ifx` (Intel).** Meson detected it because the user's environment has Intel oneAPI sourced. `meson.build` also supports gfortran via its `is_gcc` branch; the compiler is not hard-locked at baseline. (This contradicts the spec's "gfortran is the only supported compiler" statement; spec to be amended.)
- Build produces `builddir/swap` (ELF 64-bit, statically linked).
- Incremental wall time on the observation run: ~4.7 s (`builddir` was already configured). A clean build would be longer; this number is not a clean-build measurement.

## Unit tests (pFUnit) — BROKEN at baseline

`pixi run test-pfunit` fails during meson configure, before any tests compile:

```
tests/unit/meson.build:44:24: ERROR: File ../../src/io/readdra_toml.f90 does not exist.
```

Root cause: `tests/unit/meson.build` references TOML reader sources (`src/io/readdra_toml.f90`, and also `readcrop_toml.f90` and `readswap_toml.f90` per subsequent references) that were introduced during the drift period and do not exist at `e256bc0`. Zero pFUnit suites execute.

This is a pre-existing, pre-rescue condition. **Phase 1 must reconcile `tests/unit/meson.build` with the source files that actually exist at this baseline** before Phase 3 (coverage audit) can do any meaningful work. Treat this as the single largest known-debt item revealed by the baseline verification.

## Regression (full)

Full suite wall time: `5m44.812s` (pixi wall-clock). Run log preserved at `tests/regression/baselines/phase-0-baseline.log`.

| Case           | Status | Time (s) | Notes                                                           |
|----------------|--------|---------:|-----------------------------------------------------------------|
| surfacewater   | pass   |     2.84 |                                                                 |
| hupselbrook    | pass   |    10.13 |                                                                 |
| salinitystress | pass   |    26.54 |                                                                 |
| grassgrowth    | pass   |    27.13 |                                                                 |
| oxygenstress   | fail   |   177.41 | MOWDM deviations only; max diff 85.0 in year 1995 (accepted)    |
| macropore      | pass   |   342.92 | Accepted; perf-regression follow-on spec tracks it              |

Summary line from the harness: `Results: 5 passed, 1 failed`.

The log's reported "Total execution time: 342.93 s" is the max per-worker time under 6-way parallelism, not the true serial wall time. The true pixi wall time is the authoritative number (`5m44.812s`).

## Accepted tolerances

1. **`oxygenstress` MOWDM deviations.** Fixture expected vs actual differ on the MOWDM variable only, with deviations up to 85.0 in year 1995. These exist pre-rescue and are accepted through Phase 4; do not treat as a regression.
2. **`macropore` runtime ~343 s.** Roughly double the pre-drift ~170 s baseline the author recalls. Accepted for this spec; addressed in a separate performance follow-on spec between Phase 4 and the compartment-state refactor.

## Licensing finding — FLAGGED

The repository-root `LICENSE` file is **LGPL v2.1**, not GPL v2. First line:
`SOURCE:https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html`. Length 179 lines.

SWAP 4.2.0 upstream (`swap_org/source_swp_4.2.0/LICENSE`) is GPL v2. A derivative work of GPL v2 cannot be relicensed to LGPL v2.1 without permission from all copyright holders (Wageningen UR for SWAP). The modernization repo is therefore in a questionable licensing state.

This does not block Phase 0 (the orphan-branch preservation uses the GPL v2 LICENSE inside `swap_org/`, which is independent). It is a Phase 2 concern: an ADR (`docs/adr/0005-licensing-audit.md`) will investigate and document the correction path. Do not treat the current root `LICENSE` as authoritative for the derivative work until that ADR exists.

## Fast test set

Phases 1–4 use `pixi run check-fast` (to be added in Phase 1) for everyday iteration:

- surfacewater, hupselbrook, salinitystress, grassgrowth (~66 s total at the observed timings).

Full set is re-run before every phase tag via `pixi run check-full` (to be added in Phase 1). Until those pixi tasks exist, equivalent verification is `pixi run build-linux && pixi run regression` (accepting the unit-test breakage noted above until Phase 1 repairs it).

## Phase 1 debt items revealed by this baseline verification

Captured here so they survive into the Phase 1 plan authoring:

1. `tests/unit/meson.build` references missing sources — reconcile before any pFUnit work.
2. Add `check-fast` / `check-full` pixi tasks.
3. Record the intended compiler policy (gfortran-first, or gfortran + Intel both officially supported) and enforce it in `meson.build` accordingly — current behavior is "whichever the environment provides".
4. Root `LICENSE` relicense/correction — Phase 2 ADR 0005.
5. Submodule `tests/swap-cases` working-tree dirtying after each regression run — document in `docs/build-and-test.md` during Phase 2.
