# ADR 0002 — Single `builddir/` with pFUnit on by default

Status: accepted (2026-04-22, during rescue Phase 1)

## Context

Phase 0 inherited two separate meson build directories and two configure task chains:

- `builddir/` — production: `-Denable_pfunit=false`, `-Denable_unit_tests=false`. Produced the `swap` binary.
- `builddir_gfortran/` — tests: `-Denable_pfunit=true`, `-Denable_unit_tests=true`. Produced the pFUnit test binary.

Separate directories meant every compilation unit was built twice, every `meson setup` ran twice, and an operator who built with `pixi run build-linux` silently didn't build tests. Two pixi configure tasks (`_configure-swap`, `_configure-pfunit`) plus a third (`_configure-tests`) encoded this split in the task graph. A `test-unit` task targeted standalone test executables that were all removed during Task 2 reconciliation.

## Decision

A single `builddir/` used for everything. `meson_options.enable_pfunit` defaults to `true`. The `enable_unit_tests` option is removed entirely (after Task 2 there are no non-pFUnit tests). `pixi.toml` has one `_configure` task (for `FC=gfortran meson setup builddir --reconfigure`), one `build-linux` (`meson compile -C builddir`), and one `test-pfunit` (`meson test -C builddir --suite unit-pfunit --verbose`).

Two new pixi verification gates live under `[feature.test.tasks]`:

- `check-fast`: build + pFUnit + four fast regression cases. Target < 90s. Everyday iteration.
- `check-full`: build + pFUnit + all six cases. Target ~10 min. Phase-end gate.

## Consequences

- **Positive**: half the compile time for any change that isn't test-only.
- **Positive**: one obvious way to do each operation. New contributors don't need to know about the split.
- **Positive**: pFUnit tests build by default, so a configure regression in test land is caught immediately instead of when someone remembers to run the test-specific task.
- **Positive**: `check-fast` gives a predictable, quick iteration loop that every commit can gate on.
- **Neutral**: the pFUnit test binary is empty at baseline (`pf_files = []` after Task 2's reconciliation), so the test executable still gets built but produces "No suitable tests defined.". This is not wasted work — Phase 4 re-populates suites as state modules return, and the wiring is already in place.

## Override for exotic cases

Pass `-Denable_pfunit=false` to `meson setup` to build the production binary only. This is intended for CI jobs that specifically want to exercise the non-test path, not for everyday use.
