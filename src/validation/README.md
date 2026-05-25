# src/validation

## Responsibility

Stateless validator primitives consumed by per-section config validators.
Each primitive checks one invariant and appends an error on failure via an `error_collection_t` passed in.
No shared state, no side effects beyond the append.

## Public interface

All in `validation_mod`:

- `check_int_range(value, low, high, context, errors)`
- `check_real_range(value, low, high, context, errors)`
- `check_int_enum(value, allowed, context, errors)`
- `check_not_empty(value, context, errors)`
- `check_nonnegative_real(value, context, errors)`
- `check_ordered_pair(low_val, high_val, low_name, high_name, context, errors)`

## Dependencies

`error_mod`.
