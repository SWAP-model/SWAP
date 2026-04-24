# src/io/toml

## Responsibility

Composite thematic TOML reader for SWAP configuration. Provides
field-reading primitives, section readers (one per domain), a
document dispatcher (`load_swap_config`), and an emit writer
(`write_swap_config`). Phase 4a — replaces the old monolithic
`readswaptoml.f90`.

## Public interface

- `toml_field_helpers_mod` — primitives: `get_required_int`,
  `get_required_real`, `get_required_string`, the three
  optional-with-default variants, `get_table`, `get_array_of_tables`,
  `parse_date_to_days1900`. Every helper takes an `errors`
  collection, null-ptr safe, appends typed errors.
- `read_<section>_toml_mod` — per-section readers (Tasks 16-21).
- `load_swap_config_mod` — document dispatcher (Task 22).
- `write_swap_config_mod` — emit writer (Task 24).

## Dependencies

`error_mod`, `validation_mod`, `config/*` modules, `tomlf` (via
meson subproject).
