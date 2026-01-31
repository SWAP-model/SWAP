# Solute Unit Tests

This folder contains unit tests for solute transport routines.

## Test Scope

Tests in this folder should cover:

- **solute.f90**
  - Advection calculation
  - Dispersion calculation
  - Adsorption isotherm
  - Decay (first-order)
  - Root uptake of solutes
  - Mass balance

## Test Files

| File | Description |
|------|-------------|
| (to be created) | Tests for transport equation |
| (to be created) | Tests for adsorption |
| (to be created) | Tests for decay |

## Running Tests

```bash
# From the swap directory
pixi run test-unit
```

## Analytical Solutions

Reference cases for verification:
- Ogata-Banks solution (advection-dispersion)
- Pure advection (step input)
- Decay without transport
- Steady-state with decay
