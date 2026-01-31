# Utils Unit Tests

This folder contains unit tests for utility and coupling routines.

## Test Scope

Tests in this folder should cover:

- **sharedexchange.f90**
  - Data export format
  - Data import parsing
  - Error handling
  
- **sharedsimulation.f90**
  - Process coordination
  - File locking
  - Synchronization

## Test Files

| File | Description |
|------|-------------|
| (to be created) | Tests for data exchange |
| (to be created) | Tests for synchronization |

## Running Tests

```bash
# From the swap directory
pixi run test-unit
```

## Notes

- These modules are primarily for coupling with external models
- Tests may require mock external processes
- Focus on interface contracts and error handling
