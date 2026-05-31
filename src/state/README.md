# src/state

## Responsibility

The mutable runtime data model. One typed record per subsystem (`*_state_t`:
atmosphere, soilwater, solute, heat, drainage, crop*, nutrients, tillage,
timecontrol, mesh, surfacewater), aggregated under `swap_state_t` in
`swap_state.f90`. State is threaded by argument from `driver/` down and mutated
via `associate` / direct field writes.

Kept as its own layer (not folded into the feature folders) because it is the
shared, cross-cutting data the whole simulation reads — `swap_state_t` is
threaded whole through `swap_run_step`, and kernels routinely read sibling
state. It is the data model the BMI binding exposes. See ADR 0048.

## Convention

Array fields default to `allocatable` (fixed-size only if < ~64 KB) — large
fixed arrays land on the stack via `swap_state_t` and SIGSEGV on smaller stack
limits. A clean rebuild (`rm -rf builddir`) is mandatory after any schema change
here — incremental Meson does not propagate `.mod` deps across the
swap_modern↔swap_legacy boundary.
