# Phase 4f — Strangler-fig: disconnect readswap, wire TOML Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace `call ReadSwap()` in `src/core/swap.f90` with the new TOML pipeline + a new `config_to_variables` adapter. Five non-macropore regression cases run end-to-end through the new path and produce output matching their fixtures. Macropore (case 3) is excluded from `check-full` per ADR 0011. `readswap()` itself stays intact for Phase 4f-extend.

**Architecture:** New module `config_to_variables_mod` with a single public subroutine. `swap.f90` Initialize block: `load_swap_config → validate → finalize → abort_if_fatal → config_to_variables`. `run_case.sh` gains a `--toml` flag that stages TOML files into the case dir. Macropore exclusion via dropping the case from `tests/regression/test_output_regression.py`.

**Tech Stack:** Unchanged — gfortran 2008, pFUnit 4.15, meson + pixi, toml-f.

Spec: `docs/superpowers/specs/2026-04-29-phase-4f-strangler-swp-design.md`

---

## Preamble: context every task needs

**Baseline:** Phase 4f-prep complete at commit `713643a`, tag `rescue/phase-4f-prep-gap-closure`. main and development aligned. 6/6 regression green via legacy `readswap()`.

**Branch discipline:**
- Work on `development`. No feature branches.
- One commit per task. Subject `<type>(<scope>): <what> (Phase 4f Task N)`.
- No pushes during 4f.
- Phase exit: fast-forward main; tag `rescue/phase-4f-strangler-swp`.

**Working directory:** `/home/zawadzkim/Code/swap`.

**Conventions** (unchanged from 4d/4e/4f-prep).

**Verification commands:**
- `pixi run -e test build-linux` — build the binary.
- `pixi run swap --case <name> -- --toml` — run one case via the new TOML path.
- `pixi run -e test regression <name>` — compare output to fixture.
- `pixi run -e test test-pfunit` — unit tests (parity, validators, readers).
- `pixi run -e test check-full` — final gate (5 cases after macropore exclusion).

**Critical concept:**

The **iteration loop** drives Phase 4f. Each case is added one at a time. The first run of any case via the new path is highly likely to segfault or produce wrong output — the adapter is the prime suspect. Failure modes to expect:

1. **Segfault.** The adapter assigned to a not-yet-allocated array, or skipped an assignment that physics later reads.
2. **Validator/finalize abort.** A case's TOML has a value the validator rejects; loosen the validator if legacy reader was permissive.
3. **Wrong-value output.** Adapter assignment is right but a derived/normalized field needs finalize-stage transform (like Phase 4f-prep's `alphaw` normalization).

For each failure: diagnose with `print *, '[D]', config%foo, variables%foo` style probes; fix the assignment in `config_to_variables` or the validator/finalize as appropriate; re-run.

---

## File structure

### New
| File | Responsibility |
|---|---|
| `src/io/toml/config_to_variables.f90` | The adapter. Single public subroutine `config_to_variables(config)`. ~250–350 LoC. |
| `tests/unit/io/toml/test_config_to_variables.pf` | Adapter unit tests. ~150 LoC. |
| `docs/adr/0011-macropore-exclusion-from-regression.md` | Supersedes ADR 0010's runtime-fallback clause. |
| `tests/regression/baselines/phase-4f-strangler-swp.log` | Final regression log. |

### Modify
| File | Change |
|---|---|
| `src/core/swap.f90` | Initialize block: replace `call ReadSwap()` with TOML pipeline. |
| `tests/swap-cases/run_case.sh` | Add `--toml` flag + cleanup. |
| `tests/regression/test_output_regression.py` | Drop `macroporeflow` case. |
| `meson.build` | Register `config_to_variables.f90` in `sources`. |
| `tests/unit/meson.build`, `tests/unit/testSuites.inc` | Register the new test suite. |
| `docs/adr/0010-macropore-deferral.md` | Note ADR 0011 supersession. |
| `docs/adr/index.md` | Index ADR 0011. |
| `docs/configuration-schema.md` | One-line note that binary reads `swap.toml`. |

### Delete
None in Phase 4f.

---

## Tasks

### Part A — Adapter + harness changes

- [ ] **Task A1 — Author `config_to_variables` adapter (skeleton).**
  Create `src/io/toml/config_to_variables.f90` with the module shell + a stub body that includes section-by-section assignments for every (C)-classified field per `docs/phase-4f-config-to-variables-audit.md`. Walk in audit order: general → simulation → meteorology → drainage → soil → bottom_boundary → heat → irrigation → solute → surface_water → crop. For each (C)-row, emit `variables%foo = config%<section>%foo`. For allocatable arrays, guard with `if (allocated(config%foo)) ...`. For per-rotation fields, loop. Plus zero-forcing of the 18 RETIRED switches.
  Author `tests/unit/io/toml/test_config_to_variables.pf` with ~6 tests asserting representative assignments per section.
  Register both in `meson.build` and `tests/unit/meson.build` + `testSuites.inc`.
  **Verify:** test-pfunit green; new tests pass.
  **Commit:** `feat(io/toml): add config_to_variables adapter (Phase 4f Task A1)`.

- [ ] **Task A2 — Extend `run_case.sh` with `--toml` flag.**
  Add the flag to the option parser. When set: skip `cp swap_linux.swp.template swap.swp`; copy `tests/swap-cases/toml/<case>/swap.toml`, `swap.dra.toml`, `*.crp.toml` into cwd. Cleanup: extend the `rm` line to remove the new files.
  Test the flag manually: `./run_case.sh --case hupselbrook --toml --exec ../../builddir/swap` should succeed at staging (binary will fail because swap.f90 isn't wired yet, but `swap.toml` should appear in `1.hupselbrook/` before the run).
  **Verify:** smoke-stage works.
  **Commit:** `feat(test/swap-cases): add --toml mode to run_case.sh (Phase 4f Task A2)`.

- [ ] **Task A3 — ADR 0011 + macropore exclusion from regression runner.**
  Author `docs/adr/0011-macropore-exclusion-from-regression.md`: title "Macropore case exclusion from regression suite"; body documents the supersession of ADR 0010's runtime-fallback clause and explains case 3 is preserved on disk but not exercised in `check-full` going forward.
  Update `docs/adr/0010-macropore-deferral.md` with a one-line note pointing to ADR 0011.
  Update `docs/adr/index.md` to list ADR 0011.
  In `tests/regression/test_output_regression.py`: drop the `macroporeflow` entry from the registered case list. Add a comment near the drop explaining the ADR 0011 reference.
  **Verify:** `pixi run -e test regression` runs only 5 cases.
  **Commit:** `docs(adr): macropore exclusion from regression (ADR 0011) + drop case 3 from runner (Phase 4f Task A3)`.

### Part B — `swap.f90` cutover (the strangler-fig)

- [ ] **Task B1 — Wire TOML pipeline into `swap.f90` Initialize block.**
  In `src/core/swap.f90`, find the Initialize block (around line 124). Replace `call ReadSwap()` with:
  ```fortran
  block
     use load_swap_config_mod, only: load_swap_config
     use swap_config_mod, only: swap_config_t
     use config_to_variables_mod, only: config_to_variables
     use error_mod, only: error_collection_t
     type(swap_config_t)      :: config
     type(error_collection_t) :: errors
     call load_swap_config('swap.toml', config, errors)
     call config%validate(errors)
     call config%finalize(errors)
     call errors%abort_if_fatal()
     call config_to_variables(config)
  end block
  ```
  Build the binary. Don't yet run any case — that's Task B2.
  **Verify:** `pixi run -e test build-linux` succeeds.
  **Commit:** `feat(core): wire TOML pipeline into swap.f90 Initialize (Phase 4f Task B1)`.

- [ ] **Task B2 — Iterate case 1 (hupselbrook) until output matches fixture.**
  Run `pixi run swap --case hupselbrook -- --toml`. Diagnose failures (segfault → adapter omission; validator abort → schema/value mismatch; wrong output → assignment bug or finalize derivation needed).
  For each diagnosed bug: fix `config_to_variables`, the relevant validator, the finalize routine, or (if a real schema gap surfaces) extend the typed config + reader + TOML. Each fix is its own focused commit.
  When the binary runs to completion, run `pixi run -e test regression hupselbrook`. Iterate on value mismatches.
  **Verify:** `pixi run -e test regression hupselbrook` returns 1/1 green.
  **Commit:** `test(regression): hupselbrook end-to-end via TOML path (Phase 4f Task B2)`.
  (The "commit" is the marker; per-iteration fixes are their own commits before this marker. The marker commit is small or empty — possibly just a baseline log update.)

- [ ] **Task B3 — Add case 2 (grassgrowth).**
  Run `pixi run -e test regression hupselbrook grassgrowth`. Iterate.
  **Commit:** `test(regression): grassgrowth end-to-end via TOML path (Phase 4f Task B3)`.

- [ ] **Task B4 — Add case 4 (oxygenstress).**
  Run `pixi run -e test regression hupselbrook grassgrowth oxygenstress`. Iterate.
  **Commit:** `test(regression): oxygenstress end-to-end via TOML path (Phase 4f Task B4)`.

- [ ] **Task B5 — Add case 5 (salinitystress).**
  Iterate.
  **Commit:** `test(regression): salinitystress end-to-end via TOML path (Phase 4f Task B5)`.

- [ ] **Task B6 — Add case 6 (surfacewater).**
  Iterate.
  **Commit:** `test(regression): surfacewater end-to-end via TOML path (Phase 4f Task B6)`.

### Part C — Closeout

- [ ] **Task C1 — Final check-full.**
  `pixi run -e test check-full`. Expect 5/5 cases green (macropore already excluded). pFUnit unit-test suite green; F count 0.
  Save log to `tests/regression/baselines/phase-4f-strangler-swp.log`.
  **Commit:** `docs(baselines): record phase 4f check-full output (Phase 4f Task C1)`.

- [ ] **Task C2 — Schema doc one-liner.**
  Add a small note to `docs/configuration-schema.md` at the top: "As of Phase 4f, the SWAP binary reads `swap.toml` from the current directory at startup. The legacy `swap.swp` path is no longer connected to `swap.f90`."
  **Commit:** `docs(schema): note swap.toml-only binary input (Phase 4f Task C2)`.

- [ ] **Task C3 — Tag + main fast-forward.**
  Tag `rescue/phase-4f-strangler-swp HEAD -m "Phase 4f closeout: TOML pipeline + config_to_variables adapter wired into swap.f90; readswap() still in tree but disconnected; case 3 excluded from check-full per ADR 0011; 5/5 regression green."`.
  Fast-forward main: `git checkout main && git merge --ff-only development && git checkout development`.
  No push.

---

## Risks & mitigations

| Risk | Mitigation |
|---|---|
| Adapter omits a `variables%` field; segfault in case 1 | Per-iteration commit pattern catches each individually. Audit doc enumerates the list — adapter is mechanical. |
| Validator/finalize fails on a real regression-case TOML | Loosen with a `READING NOTE`; reference legacy permissive behavior. |
| TOML schema is missing a key case N actually uses | Surface during iteration; extend schema as a focused commit. Phase 4f-extend addresses systematically. |
| `bbcfil` external file resolution doesn't work after cwd change | Path-resolution already proven in parity tests; the run_case.sh `cd` puts cwd in the case dir before binary launches. |
| Cross-file TOML loading (`drainage.file`, `crop.rotation.file`) breaks | Same as bbcfil — parity tests prove it works; first failure surfaces in Task B2. |
| RETIRED-switch zero-forcing breaks output paths | Adapter explicitly zeros all 18 per ADR 0009. CSV-only output (per ADR 0009) means the comparison against fixtures uses only the `swcsv`/`swcsv_tz` paths. |
| Per-rotation crop fields aren't fully populated by the adapter | Crop section of audit has 30 G entries (mostly false-R per Phase 4e Task B5 triage). Iterate; extend adapter as needed. |

---

## Verification gates

After each Task A: build succeeds; smoke-test passes; ADR/runner exclusion verified.
After each Task B: `regression <case_list>` green for the cases added so far.
After Task C1: `check-full` 5/5 green; pFUnit green; F count 0.
After Task C3: tag exists; main fast-forwarded.

Final audit-trail check: `git diff rescue/phase-4f-prep-gap-closure..HEAD -- src/` should show: one Initialize-block change in `swap.f90`, the new `config_to_variables.f90`, and any focused fixes to typed configs / readers / finalizers that surfaced during iteration. No physics-algorithm changes.

## File count summary

- Create: ~4 files (adapter, adapter test, ADR 0011, baseline log).
- Modify: ~5 files (swap.f90, run_case.sh, test_output_regression.py, meson configs, schema doc note).
- Per-iteration fixes (Tasks B2-B6): variable count, ~5-15 focused commits depending on what surfaces.
- LoC delta: ~+400 / ~-15.
