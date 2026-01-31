# Soil Unit Tests

This folder contains unit tests for soil physics and water flow routines.

## Test Scope

Tests in this folder should cover:

### Hydraulic Functions
- **functions.f90**
  - `watcon` - water content calculation
  - `hconduc` - hydraulic conductivity
  - `moiscap` - moisture capacity
  - Mean conductivity calculations
  
- **WC_K_models_04_11.f90**
  - Mualem-van Genuchten equations
  - Bimodal retention curves
  - PDI model
  
- **hysteresis.f90**
  - Scanning curve parameters
  - Reversal detection

### Numerical Solution
- **headcalc.f90**
  - Convergence behavior
  - Boundary condition application
  
- **tridag.f90**
  - Tridiagonal solver accuracy
  - Error handling

### Grid and State
- **calcgrid.f90**
  - Node position calculations
  - Layer assignments
  
- **calcgwl.f90**
  - Groundwater level detection
  - Perched water tables

## Existing Test Files

| File | Description |
|------|-------------|
| `test_van_genuchten.pf` | Tests for van Genuchten equations |
| `test_functions.pf` | Tests for hydraulic functions |

## Running Tests

```bash
# From the swap directory
pixi run test-unit
```

## Test Cases

Reference solutions from analytical cases:
- Steady-state infiltration (Green-Ampt)
- Horizontal infiltration (Philip)
- Drainage to water table
