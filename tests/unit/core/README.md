# Core Unit Tests

This folder contains unit tests for the core framework modules.

## Test Scope

Tests in this folder should cover:

- **variables.f90**
  - Variable initialization
  - Module interface verification
  
- **initialize.f90**
  - Proper initialization of all state variables
  - Default value verification
  
- **timecontrol.f90**
  - Time step management
  - Day/year boundary handling
  - Output timing logic
  
- **swap.f90**
  - Task routing (init/dynamic/close)
  - DLL interface verification
  
- **swap_main.f90**
  - Command line argument parsing
  - Rerun functionality

## Test Files

| File | Description |
|------|-------------|
| (to be created) | Tests for core functionality |

## Running Tests

```bash
# From the swap directory
pixi run test-unit
```

## Notes

- Use pFUnit framework for Fortran unit tests
- Mock external file I/O where possible
- Focus on isolated function behavior
