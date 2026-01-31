# Core Module

This folder contains the core framework files that provide the foundation for the SWAP model.

## Files

### Include Files

| File | Description |
|------|-------------|
| `arrays.fi` | Parameter definitions for array dimensions (maximum compartments, drainage systems, crops, etc.) |
| `params.fi` | Common numerical parameters and constants (`nihil`, `tiny`, `small`, `vlarge`) |
| `description.fi` | Model version and metadata (version number, license, author information) |

### Core Modules

| File | Description |
|------|-------------|
| `variables.f90` | Central module containing all shared state variables organized by category (time/control, meteorology, irrigation, crop, soil water, macropore, surface water, heat, snow, solute). This is the main data container for the entire model. |
| `initialize.f90` | Initialization routines for all model variables. Sets default values for all state and rate variables at the start of a simulation. |
| `timecontrol.f90` | Time stepping and control logic. Manages simulation time, day/year transitions, output timing, and adaptive time step control. |
| `swap.f90` | Main SWAP model subroutine orchestrating the three-phase simulation (initialization, dynamics, closure). Contains the `swap_exchange` module for DLL interface and the main `swap` subroutine. |
| `swap_main.f90` | Entry point program for standalone SWAP execution. Handles reruns, command-line interface, and calls the main `swap` subroutine. |

## Architecture

The core module follows a three-task architecture:
1. **Task 1 - Initialization**: Read input files, set initial conditions
2. **Task 2 - Dynamic**: Time loop with iterative solving of water flow equations  
3. **Task 3 - Closure**: Write final outputs, close files, cleanup

The `variables.f90` module serves as the central data hub, with all other modules importing variables from here via `use variables`.

## Dependencies

- This module has no internal dependencies on other SWAP modules
- All other SWAP modules depend on `core/variables.f90`
