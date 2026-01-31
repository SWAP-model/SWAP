# Atmosphere Unit Tests

This folder contains unit tests for atmospheric process routines.

## Test Scope

Tests in this folder should cover:

- **penmon.f90**
  - Reference ET calculation
  - Radiation balance components
  - Aerodynamic resistance
  - Canopy resistance
  
- **meteoday.f90**
  - Interception calculations
  - Von Hoyningen-Hüne/Braden method
  - Gash model
  - Rutter model
  
- **meteodt.f90**
  - Sub-daily ET distribution
  - Sine wave partitioning
  - Rain event processing
  
- **snow.f90**
  - Snow accumulation
  - Degree-day melting
  - Rain-on-snow events

## Test Files

| File | Description |
|------|-------------|
| (to be created) | Tests for Penman-Monteith |
| (to be created) | Tests for interception models |
| (to be created) | Tests for snow dynamics |

## Running Tests

```bash
# From the swap directory
pixi run test-unit
```

## Reference Data

- FAO-56 reference ET examples
- Published interception measurements
- Snow course observations
