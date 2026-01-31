# Macropore Unit Tests

This folder contains unit tests for macropore flow routines.

## Test Scope

Tests in this folder should cover:

- **macropore.f90**
  - Macropore geometry calculations
  - Volume at depth
  - Domain partitioning
  
- **macrorate.f90**
  - Surface entry flux
  - Lateral absorption
  - Rapid drainage flux
  - Water level in macropores
  
- **macroporeoutput.f90**
  - Output variable calculations
  - Mass balance consistency

## Test Files

| File | Description |
|------|-------------|
| (to be created) | Tests for macropore geometry |
| (to be created) | Tests for water flow rates |
| (to be created) | Tests for absorption |

## Running Tests

```bash
# From the swap directory
pixi run test-unit
```

## Test Cases

- Known macropore geometry scenarios
- Limiting cases (no macropores, fully cracked)
- Absorption into initially dry soil
