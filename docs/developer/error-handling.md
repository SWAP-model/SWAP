---
title: Error handling
author: SWAP modernization team
---

# Error handling

## Overview

`error_mod` (at `src/error/error.f90`) is the single error-reporting
facility for new SWAP infrastructure. Legacy modules still use
`fatalerr` (ttutil); they migrate per-module as Phase 4+ touches them.

## Contract

Every fallible procedure takes an `errors` argument as `intent(inout)`:

    subroutine read_something(...., errors)
       type(error_collection_t), intent(inout) :: errors
       ...
    end subroutine

On failure, append an error and return. Do not halt:

    if (bad_condition) then
       call errors%append(ERR_PARSE_TYPE_MISMATCH,  &
                          "expected integer",        &
                          "drainage.basic.general.dramet")
       return
    end if

`append` auto-logs via `swap_log%log_error`. Do not log additionally.

## Error codes

| Prefix | Domain | Typical fatal? |
|---|---|---|
| `ERR_IO_*` | Filesystem / file access | yes |
| `ERR_PARSE_*` | TOML decoding / field helpers | yes |
| `ERR_VALIDATION_*` | Per-field / cross-field / cross-section checks | yes |
| `ERR_FINALIZE_*` | Derived-value computation | often no (warnings) |
| `ERR_ADAPTER_*` | Config-to-state mismatch | yes |

`is_fatal` defaults to `.true.`. Pass `is_fatal=.false.` for
non-blocking warnings.

## Abort checkpoint

One abort point after the full pipeline (usually: after `finalize`):

    call load_swap_config(path, config, errors)
    call config%validate(errors)
    call config%finalize(errors)
    call errors%abort_if_fatal()

`abort_if_fatal` writes `summary()` to stderr and calls `error stop`
if any appended error has `is_fatal=.true.`. It returns normally
otherwise.

## Testing

`error_collection_t` is a plain value type — construct one in a
test, call procedures that append to it, assert on `%has_errors`,
`%has_fatals`, `%count`, and `%items(i)%code`. See
`tests/unit/error/test_error.pf` for examples.
