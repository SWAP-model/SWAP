---
title: "ADR 0011 — Macropore case exclusion from regression suite"
date: 2026-04-29
status: accepted
supersedes-clause: "ADR 0010 §6 (runtime fallback)"
---

# ADR 0011: Macropore case exclusion from regression suite

## Context

[ADR 0010](0010-macropore-deferral.md) deferred the macropore module:
schema kept as orphan infrastructure under
`src/config/macropore_config.f90`, no wiring of the new TOML
pipeline. Clause §6 of ADR 0010 stated:

> "Phase 4f's strangler-fig keeps `readswap.f90` as a fallback
> for macropore-using cases. The runtime branch is:
> `if (toml_path_available) then call load_swap_config(...) else
> call readswap()`. Case 3 retains its `swap_linux.swp.template`
> and runs on the legacy reader."

Phase 4f's design has since changed. The user's directive on
2026-04-29: "Remember to exclude from now on the macropore case
study from the full test suite (as we agreed before)." With case
3 out of `check-full`, the runtime fallback is no longer needed —
the new TOML pipeline becomes the only entry point for production
runs.

## Decision

**Case 3 (`3.macroporeflow`) is removed from `check-full`.** The
case directory remains on disk (`tests/swap-cases/3.macroporeflow/`)
for archival and as a reference for future macropore work, but
`tests/regression/test_output_regression.py` no longer iterates
over it. `pixi run -e test check-full` runs 5 cases; full
regression coverage of macropore-affected physics is deferred to
the future macropore phase.

**`swap.f90`'s runtime fallback to `readswap()` is dropped.** Phase
4f's strangler-fig becomes a clean cut: the binary loads
`swap.toml` at startup, validates, finalizes, calls
`config_to_variables(config)`, and proceeds. No legacy fallback,
no runtime branch.

**`readswap.f90` itself stays in the tree** (per the broader Phase
4f design). It compiles, the test-side `legacy_crop_helper.f90`
continues to call its sub-readers from parity tests, and Phase
4f-extend will walk it line by line to add schema slots for any
inputs the 5 active test cases don't exercise. Once Phase 4f-extend
is done and the schema covers every legacy input, a follow-up
phase deletes the file.

This ADR supersedes Clause §6 of ADR 0010. The rest of ADR 0010
(macropore_config_t kept as orphan, no [macropore] schema reading,
schema preserved for future macropore work) remains in effect.

## Consequences

Positive:

- Phase 4f's `swap.f90` cutover is a clean cut, not a runtime
  branch. Less complexity in the production binary.
- `check-full` is faster (~5 minutes saved by skipping macropore's
  ~240s run).
- The new TOML pipeline becomes the canonical entry point;
  there's no "second way" via legacy `readswap()` once the
  strangler-fig lands.

Negative:

- Macropore physics regressions go uncaught until the future
  macropore phase reactivates case 3.
- Anyone reading the case directory and expecting it to "just
  work" via `pixi run swap --case macroporeflow` will see a clear
  failure: no `[macropore]` section in any TOML, validator either
  errors out or runs without macropore physics. Documenting the
  exclusion in the case directory's README is worth doing if the
  case grows users.

## Revisit trigger

When the future macropore phase opens. At that point: extend or
replace `macropore_config_t`, populate case 3's TOML with
`[macropore]`, extend the parity test, re-add case 3 to
`check-full`, and drop this ADR's exclusion clause.
