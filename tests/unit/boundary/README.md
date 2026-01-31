# Boundary Unit Tests

This folder contains unit tests for boundary condition routines.

## Test Scope

Tests in this folder should cover:

- **boundtop.f90**
  - Infiltration calculation
  - Evaporation limiting
  - Ponding/runoff logic
  - Macropore surface entry
  
- **boundbottom.f90**
  - All SWBOTB options (1-8)
  - Free drainage (unit gradient)
  - Aquifer coupling
  - Lysimeter condition

## Test Files

| File | Description |
|------|-------------|
| (to be created) | Tests for top boundary |
| (to be created) | Tests for bottom boundary options |

## Running Tests

```bash
# From the swap directory
pixi run test-unit
```

## Test Cases

### Top Boundary
- High intensity rainfall (ponding)
- Low intensity rainfall (direct infiltration)
- Evaporation demand vs supply

### Bottom Boundary
- Prescribed head matching
- Free drainage profile
- Seepage face conditions
