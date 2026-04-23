# src/boundary

## Responsibility

Top and bottom boundary conditions for the Richards-equation soil column.
`boundtop.f90` supplies the top boundary flux (ponding, runoff, and
`PONDRUNOFF`), while `boundbottom.f90` computes the bottom boundary
(free drainage, pressure head, flux, or Cauchy conditions depending on
the user's `swbotb` choice). Constants live in `boundary_constants.f90`;
`boundary_state.f90` declares `boundary_state_t` and its lifecycle
routines.

## Public interface

- `boundary_state_t` — aggregated boundary state type.
- `boundary_state_init`, `boundary_state_finalize` — lifecycle.
- `boundtop`, `boundtop_state` — top boundary flux.
- `PONDRUNOFF`, `PONDRUNOFF_state` — ponding / runoff partitioning.
- `BoundBottom`, `BoundBottom_state` — bottom boundary flux.

## Dependencies

Imports from:

- `src/core/` — `swap_array_dimensions`, `swap_log`, `swap_state_mod`,
  `swap_state_sync`, `variables`.
- `src/utils/` — `array_utils`, `soilhydraulics_utils`,
  `surfacewater_utils`.

Does not depend on other physics subdirs.
