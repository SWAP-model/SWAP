# src/error

## Responsibility

Placeholder for the unified error-handling module planned for Phase 4 of
the rescue effort. The directory is intentionally empty at the rescue
baseline.

Today, errors propagate through a mix of `swap_log`'s `log_error`,
bespoke `write(*,*) + stop` calls, and caller-checked return flags.
Phase 4 will consolidate these into a single `error_mod` with a typed
error payload and a consistent propagation policy. See
[../../docs/architecture.md](../../docs/architecture.md) and the
phased plan in the rescue docs for details.

## Public interface

None — the subdir currently contains no source files.

## Dependencies

None.
