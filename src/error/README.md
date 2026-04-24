# src/error

## Responsibility

Typed error payload + accumulating collection used across the new
configuration pipeline (Phase 4a). Every fallible procedure takes
`errors` as `intent(inout)` and appends failures; one checkpoint
(`abort_if_fatal`) aborts if any fatal error was collected.

See `docs/error-handling.md` for the usage guide and
`docs/adr/0008-error-collection-over-fatalerr.md` for the decision record.

## Public interface

- `error_t` — single error payload (`code`, `message`, `context`, `is_fatal`).
- `error_collection_t` — accumulator with type-bound procedures:
  `append`, `has_errors`, `has_fatals`, `count`, `summary`,
  `abort_if_fatal`, `clear`.
- Error code constants: `ERR_NONE`, `ERR_IO_READ_FAILED`, `ERR_IO_WRITE_FAILED`,
  `ERR_PARSE_MALFORMED_TOML`, `ERR_PARSE_TYPE_MISMATCH`,
  `ERR_PARSE_MISSING_REQUIRED`, `ERR_VALIDATION_OUT_OF_RANGE`,
  `ERR_VALIDATION_ENUM`, `ERR_VALIDATION_CROSS_FIELD`,
  `ERR_VALIDATION_CROSS_SECTION`, `ERR_FINALIZE_DERIVATION`,
  `ERR_ADAPTER_UNSUPPORTED`.

## Dependencies

`swap_log` (for auto-logging on append).
