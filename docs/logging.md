---
title: Logging
author: SWAP modernization team
---

# Logging

## Overview

`swap_log` (at `src/core/swap_log.f90`) is the single logging facility for
the modernised SWAP tree. Every module that needs to emit diagnostics uses
it. The error module (`src/error/error.f90`, Phase 4a) routes every
appended error through `log_error` automatically, so application code never
needs to call logger and error separately.

## Levels

Four levels, lowest to highest severity:

| Constant         | Integer value | Typical use                                   |
|---|---|---|
| `LOGLEVEL_DEBUG` | 10            | High-volume developer traces; off by default. |
| `LOGLEVEL_INFO`  | 20            | Progress messages; default threshold.         |
| `LOGLEVEL_WARN`  | 30            | Unexpected but recoverable conditions.        |
| `LOGLEVEL_ERROR` | 40            | Every appended `error_t` logs here.           |

Plus `LOGLEVEL_NONE` (100) to silence the logger entirely.

## Usage

Initialise once at program start (typically in `swap_main.f90` or its
equivalent in test drivers):

    use swap_log
    call log_init(log_level=LOGLEVEL_INFO, log_file='swap.log')

Then anywhere:

    call log_info('reader',  'Loading TOML file: ' // trim(path))
    call log_warn('finalize', 'Derived value clipped to bounds')
    call log_error('validate', 'drainage.dramet out of range')

Close on shutdown:

    call log_close

`log_debug`, `log_info`, `log_warn`, `log_error` take `(context, message)`.
`context` is a short module-or-subsystem tag (e.g., `"reader"`, `"validate"`).
Message is any string; use `to_str(n)` from the same module to stringify
integers, reals, logicals.

## Output shape

    2026-04-24 09:42:17.312 INFO  reader: Loading TOML file: swap.toml

Format is fixed: `TIMESTAMP LEVEL CONTEXT: MESSAGE`. Timestamps can be
suppressed via `log_init(..., timestamps=.false.)`.

## Interaction with the error module

When `src/error/error.f90` (Phase 4a) appends to an `error_collection_t`,
`swap_log%log_error(context, message)` is called automatically. Do not
log errors manually in addition — the append call handles it.

## Thread safety

None. `swap_log` mutates module-level state (`current_level`, `log_unit`,
file-handle state) without synchronisation. Safe for serial gfortran; when
multicore work lands in a future phase, the logger either becomes
thread-local or gains synchronisation. Do not rely on logger behaviour
under concurrent calls today.

## Testing

pFUnit unit tests for the logger are added as part of Phase 4a Task 2
(they land alongside the error-module tests, since those exercise the
logger-sink path).
