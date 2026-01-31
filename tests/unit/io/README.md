# I/O Unit Tests

This folder contains unit tests for input/output routines.

## Test Scope

Tests in this folder should cover:

- **readswap.f90**
  - Input file parsing
  - Parameter validation
  - Error handling for malformed input
  
- **readmeteo.f90**
  - Weather file parsing
  - Missing value detection
  - Date handling
  
- **swapoutput.f90**
  - Output file creation
  - Format verification
  - Balance output accuracy
  
- **swap_csv_output.f90**
  - CSV format correctness
  - Variable selection logic
  - Header generation

## Test Files

| File | Description |
|------|-------------|
| (to be created) | Tests for I/O functionality |

## Test Data

Test input files should be placed in the `tests/cases/` directory.

## Running Tests

```bash
# From the swap directory
pixi run test-unit
```

## Notes

- Test with both valid and invalid input files
- Verify error messages are helpful
- Check boundary conditions (empty files, very long lines, etc.)
