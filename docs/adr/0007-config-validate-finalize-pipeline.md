---
title: "ADR 0007 — Config parse-validate-finalize-adapter pipeline"
date: 2026-04-24
status: accepted
---

# ADR 0007: Four-stage config pipeline

## Context

The legacy `readswap.f90` (5161 lines) mixes four concerns: file I/O,
type deserialization, input validation, and population of model state.
That coupling makes any single concern hard to test in isolation and
blocks the planned Python bindings path, which needs to invoke
validation without going through the file reader.

## Decision

Split the input pipeline into four distinct stages:

1. **Parse** — `load_swap_config(path, config, errors)`. Byte-level
   TOML decoding; thin lossless conversions (dates, units) at the
   boundary. No semantic checks.
2. **Validate** — `config%validate(errors)`. Per-section plus
   aggregate cross-section invariants. Read-only on the config.
3. **Finalize** — `config%finalize(errors)`. Derived values, array
   expansion, canonical-form normalization. Mutates the config.
   Idempotent. Runs only when validate produced no fatal errors.
4. **Adapter** — `config_to_state(config, state)`. Temporary in
   Phase 4a; retires when state types merge with config types in a
   later phase.

A parallel fifth stage — **emit** (`write_swap_config(config, path,
errors)`) — is implemented as a separate module for round-trip
tests and future Python-driven persistence.

## Consequences

Positive:

- Each stage is unit-testable in isolation.
- Python bindings can drive validate + finalize + adapter without
  ever touching the parser.
- Round-trip tests (load → emit → load) surface a class of reader
  bugs that field-by-field diff does not.
- Errors accumulate through the whole pipeline in one collection;
  one abort checkpoint after finalize.

Negative:

- Four types and four stages is more infrastructure than legacy
  `readswap.f90`. The cost is paid once; testability and
  reusability are permanent.

## Revisit trigger

When the state layout folds into config types (Phase 4 or later),
the adapter stage disappears. The other three stages remain.
