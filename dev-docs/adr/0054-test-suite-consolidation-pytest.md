# ADR 0054 — Test-suite consolidation & pytest standardization

**Status:** Accepted — landed 2026-07-05
**Arc:** Test-suite review (post-migration cleanup, before resuming feature work)
**Relates to:** ADR 0004 (pFUnit for unit tests) — extends it with the
coverage-responsibility map; ADR 0050 (one libswap.so / binding smokes). Does
not supersede an existing ADR.

## Context

The test suite grew organically across the whole strangler-fig migration and had
accumulated overlap, dead wiring, and low-value tests. An in-depth review
(2026-07-05) found:

- **Two Fortran layers with a clean split** — the pFUnit unit suite and the
  byte-identical Python regression suite — but the *Python* side was a patchwork
  of ad-hoc `assert`+`sys.exit` scripts (`hello_swap.py`, `run_ensemble.py`, the
  coupling `test_*.py`) driven three different ways: raw `python`, meson `test()`
  suites, and `scripts/check_bindings.sh`.
- **Dead wiring**: four pixi run-tasks and the meson `bmi`/`cffi-demo` suites all
  pointed at `tests/swap-cases`, a directory removed when the cases moved to the
  `swap-testcases` sibling repo. `check_bindings.sh` and the meson suites had
  been silently failing or duplicating each other. 5.5 MB of non-gfortran oracle
  binaries (`swap420`, `swap.exe`) were tracked but unused (the project is
  gfortran-only, ADR 0001).
- **Low-value unit tests**: ~26 `*_phase0_fields` assertions checked Fortran's
  own default zero-init / allocatable semantics (compiler behavior, not SWAP
  behavior); five loader "dispatches_X" tests were near-duplicates of the direct
  `read_*_toml` reader tests.
- **A coverage picture easily misread**: a naive "which module does a `.pf`
  `use`?" scan flags heat, soil-water physics, snow, runoff, oxygen stress as
  "untested" — but those are exactly the physics kernels the regression suite
  owns and CLAUDE.md forbids unit-testing against hand-derived values. The
  genuine non-physics gaps were a short list.

## Decision

**1. Coverage-responsibility map (who owns correctness for what).**

| Layer | Owns | Tested on |
|---|---|---|
| pFUnit unit suite | config parsing/validation, state seeding, I/O glue, pure helpers | small synthetic fixtures (per-section TOMLs, `runnable_model`) |
| Byte-identical regression (pytest) | physics kernels (soil-water, heat, solute, atmosphere, crop growth, drainage) | full-model runs vs the swap420gf 4.2.0 oracle |
| pytest binding smokes | the C-ABI facades (`swap_bmi_mod`, `swap_capi_mod`, `swap_xmi_mod`) | end-to-end load/step of `libswap.so` |
| (indirect) | legacy external I/O subroutines (`file_headers`), crop-variant state seeding | regression output / crop init tests |

Physics kernels are **deliberately not** unit-tested — the byte-identical
regression is a stronger oracle, and CLAUDE.md forbids touching physics. So
their absence from the pFUnit suite is by design, not a gap.

**2. Python harnesses standardize on pytest.** The regression runner is a
parametrized pytest module (one `test_regression[<case>]` per registered case,
`fast` marker for the check-fast subset, xdist for parallelism); its case
registry + run/aggregate/compare logic lives in `regression_harness.py` (shared
with the two regen scripts). `known_divergence`/`pending_restore` become
non-strict `xfail` markers, preserving the prior pass/xfail policy exactly. The
binding/coupling smokes become pytest functions over shared fixtures
(`tests/conftest.py`: `swap_lib`, `libmf6` skip-if-absent, `bindings_case`).
The former standalone CLI is retired; diagnostic modern-build snapshots move to
`regen_expected.py`.

**3. One home per test path — no framework overlap.** The meson `bmi`,
`cffi-demo`, and `coupling` Python-test suites are deleted; `check_bindings.sh`
(pytest-driven) is the single home of the C-ABI smokes, and `check-fast` /
`check-full` are the single home of the regression gate. Dead pixi tasks and the
unused oracle binaries are removed.

**4. Trim compiler-behavior tests.** The `*_default_zero` / `*_default_unallocated`
assertions are removed (kept: allocation-shape, TOML round-trip, validator
tests). The five loader dispatch tests collapse to one "wires all optional
sections" test (the direct `read_*_toml` tests never exercise the top-level
loader, so this coverage is kept, not dropped).

**5. Fill the genuine glue gaps.** New pFUnit tests for `csv_common_mod` (the
shared days-since-1900 date math every CSV loader depends on) and
`meteo_buffer_mod` (the coupled external-meteo buffer) — both previously at zero
coverage with clean, deterministic interfaces.

## Consequences

- pFUnit 746 → 725 tests (−30 trivia/dedup, +9 new glue coverage); regression
  18 pass / 2 xfail unchanged; `check-bindings` 5 pytest smokes.
- New `.pf` files must be registered in **both** `tests/unit/meson.build`
  (`pf_files`) and `testSuites.inc` — the meson list is explicit, not globbed.
- `pytest.ini` (repo root) registers the `fast` marker and keeps collection off
  the vendored `tests/pFUnit` tree. pytest/pytest-xdist were already in the test
  env.
- The framework provisioning is unchanged: `tests/pFUnit/` stays git-ignored and
  built on demand by `scripts/setup_pfunit.sh` (not vendored).
- Not done here (deferred, low marginal value): isolated unit tests for the C-ABI
  facades (covered by the pytest smokes) and crop-variant state seeding (covered
  by the crop init tests + regression).
