# Drainage Unit Tests

This folder contains unit tests for lateral drainage and surface water routines.

## Test Scope

Tests in this folder should cover:

- **drainage.f90**
  - Hooghoudt drainage equation
  - Ernst equation
  - Entry resistance calculations
  - Multi-level drainage
  
- **divdra.f90**
  - Flux distribution algorithms
  - Transmissivity calculations
  - Layer boundary detection
  
- **surfacewater.f90**
  - Surface water level dynamics
  - Weir discharge
  - Runoff calculations

## Test Files

| File | Description |
|------|-------------|
| (to be created) | Tests for drainage equations |
| (to be created) | Tests for flux distribution |
| (to be created) | Tests for surface water |

## Running Tests

```bash
# From the swap directory
pixi run test-unit
```

## Analytical Solutions

Reference cases for verification:
- Steady-state Hooghoudt (known exact solution)
- Falling water table (Glover-Dumm)
- Ernst equation with known layer properties
