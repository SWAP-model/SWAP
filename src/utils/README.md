# src/utils

## Responsibility

Low-level shared helpers used across all physics domains. Array
utilities including `afgen` linear interpolation and `insw` / `stepnr`
/ `interpol` (`arrayutils.f90`). Numerical solvers — tridiagonal
(`tridag`) and banded LU (`bandec` / `banbks`) — in
`numericalsolvers.f90`. Hydraulics helpers `watcon`, `moiscap`,
`hconduc`, `dhconduc`, `prhead`, `hcomean`, `dkmean` in
`soilhydraulicsutils.f90`. Surface-water helpers `wlevst`, `swstlev`,
`qhtab`, `runoff` in `surfacewaterutils.f90`. Shared-model coupling
scaffolding in `sharedexchange.f90` and `sharedsimulation.f90`.

## Public interface

- `afgen`, `stepnr`, `insw`, `interpol` — array / interpolation
  helpers.
- `tridag`, `bandec`, `banbks` — numerical solvers.
- `watcon`, `moiscap`, `hconduc`, `dhconduc`, `prhead`, `hcomean`,
  `dkmean` — soil-hydraulics helpers.
- `wlevst`, `swstlev`, `qhtab`, `runoff` — surface-water helpers.
- `sharedexchange`, `sharedsimulation` — coupling-with-external-
  models scaffolding.

## Dependencies

Imports from:

- `src/core/` — `swap_array_dimensions`, `variables` (legacy bridge;
  will shrink as the rescue progresses).
- `src/soil/` — `doln`, `WC_K_models_04_11` used by the soil-
  hydraulics helpers.
- Standard library — `iso_fortran_env`.

Does not depend on the other physics subdirs. Utils is deliberately
leaf-level — adding a dependency on `atmosphere`, `boundary`, `crop`,
`drainage`, `heat`, `macropore`, or `solute` here would be a layering
violation.
