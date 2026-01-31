# Heat Unit Tests

This folder contains unit tests for soil temperature and frost routines.

## Test Scope

Tests in this folder should cover:

- **temperature.f90**
  - Heat conduction solution
  - Thermal property calculations
  - Analytical solution (sine wave)
  - Boundary condition application
  
- **frozencond.f90**
  - Frost depth calculation
  - Conductivity reduction factor
  - Freeze-thaw dynamics

## Test Files

| File | Description |
|------|-------------|
| (to be created) | Tests for heat conduction |
| (to be created) | Tests for thermal properties |
| (to be created) | Tests for frost |

## Running Tests

```bash
# From the swap directory
pixi run test-unit
```

## Analytical Solutions

Reference cases for verification:
- Harmonic temperature wave
- Step change at surface
- De Vries thermal conductivity
