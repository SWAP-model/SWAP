---
title: "ADR 0008 — error_collection_t over fatalerr"
date: 2026-04-24
status: accepted
---

# ADR 0008: Error collection over fatalerr

## Context

The legacy codebase uses `fatalerr` (from ttutil) for error handling:
275 call sites, each halts the process on the first failure. This
gives users one error at a time even when a config file has many
problems. Validation of a full SWAP input set should report every
issue in one pass so the user can fix them all at once.

## Decision

Introduce `error_collection_t` in `src/error/error.f90`. Every
fallible procedure in new infrastructure takes `errors` as
`intent(inout)` and appends to the collection. No direct aborts
inside stages; one explicit `abort_if_fatal` checkpoint after the
pipeline completes. Every append auto-logs via `swap_log%log_error`.

Legacy `fatalerr` sites are **not** migrated in Phase 4a. They
migrate per-module during Phase 4+ as each module is touched,
matching the existing incremental-cleanup rule.

## Consequences

Positive:

- Users see every config problem in one run.
- The collection is a value type, passable via `iso_c_binding` when
  Python bindings land.
- Errors and logs stay in sync because append auto-logs.

Negative:

- Calling convention is now mandatory: every fallible procedure
  gains an `errors` argument. Legacy code without it stays on
  `fatalerr` until migration.

## Revisit trigger

When Phase 4+ migrates the last `fatalerr` call site.
