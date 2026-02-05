---
description: A description of your rule
---

# SWAP Project Architecture

## Current State (Phase 2: Integration)

### Three-Layer Architecture

┌─────────────────────────────────────────────────┐
│ Layer 1: Input/Legacy (ttutil + Variables) │
│ │
│ - Reads .swp/.crp/.dra input files │
│ - Populates Variables module (module globals) │
│ - Legacy compatibility maintained │
└─────────────────────────────────────────────────┘
↓
state_from_variables() ← Sync Layer
↓
┌─────────────────────────────────────────────────┐
│ Layer 2: Computation (State-based, IN PROGRESS)│
│ │
│ - Main time loop operates on swap_state_t │
│ - Process routines accept state arguments │
│ - NO module globals in computation │
└─────────────────────────────────────────────────┘
↓
state_to_variables() ← Sync Layer
↓
┌─────────────────────────────────────────────────┐
│ Layer 3: Output (Variables-based) │
│ │
│ - Output routines read Variables module │
│ - Legacy format compatibility │
└─────────────────────────────────────────────────┘

text

## End Goal Architecture (Phase 3: BMI)

┌─────────────────────────────────────────────────┐
│ Python BMI Interface │
│ │
│ bmi.initialize(config_dict) → state │
│ bmi.update() → calls Fortran with state │
│ bmi.get_value('theta') → reads from state │
│ bmi.finalize() → cleanup │
└─────────────────────────────────────────────────┘
↓
┌─────────────────────────────────────────────────┐
│ SWAP Core (Pure State-based) │
│ │
│ call soilwater(state%soil, state%boundary) │
│ call cropgrowth(state%crop, state%soil) │
│ call drainage(state%drain, state%soil) │
│ │
│ NO Variables module dependencies │
│ NO sync layer needed │
└─────────────────────────────────────────────────┘
↓
┌─────────────────────────────────────────────────┐
│ Legacy Wrapper (Optional) │
│ │
│ call swap_legacy(input_files) │
│ → Uses Variables + ttutil │
│ → Calls state-based core via wrapper │
│ → Backward compatibility preserved │
└─────────────────────────────────────────────────┘

text

## Key Components

### State Types (`src/core/swap_state_mod.f90`)

```fortran
type :: swap_state_t
  type(time_state_t)        :: time
  type(soil_state_t)        :: soil
  type(crop_state_t)        :: crop
  type(atmosphere_state_t)  :: atm
  type(boundary_state_t)    :: boundary
  type(drainage_state_t)    :: drain
  type(solute_state_t)      :: solute
  ! ... etc
end type

Status: ✅ Complete (Phase 1)
Sync Layer (src/core/swap_state_sync.f90)

fortran
! Variables → State
subroutine state_from_variables(state)

! State → Variables  
subroutine state_to_variables(state)

Current usage: Inside time loop (temporary)
Goal: Outside time loop only (input/output boundaries)
End state: Only in legacy wrapper (Phase 3)
Variables Module (src/core/variables.f90)

Current role: Primary data storage (module globals)
Phase 2 role: Gradually reduce usage in computation
Phase 3 role: Input/output and legacy wrapper only
Never delete: Preserved for backward compatibility
Process Subsystems

Located in:

    src/soil/ - Soil water flow, heat, transport

    src/crop/ - Crop growth, phenology, water use

    src/atmosphere/ - Meteorology, ET, interception

    src/drainage/ - Subsurface drainage, deep percolation

    src/macropore/ - Preferential flow

    src/solute/ - Solute transport

Integration order: See .continue/context/integration-patterns.md
Phase 2 Integration Progress
Completed (Phase 1)

    ✅ Module-level SAVE removal (Steps 1-10)

    ✅ State type definitions

    ✅ Sync layer implementation

    ✅ Unit tests for all state types

In Progress (Phase 2)

    🔄 Main time loop state wiring

    🔄 Process subroutine conversions

    🔄 Sync elimination from computation

Future (Phase 3)

    ⏳ Python BMI interface

    ⏳ Direct state initialization from Python

    ⏳ Legacy wrapper finalization

    ⏳ Variables module isolation

Multi-Instance Support

Goal: Multiple SWAP instances in same process

fortran
! Independent instances - no state leakage
type(swap_state_t) :: instance1, instance2

call initialize(instance1, config1)
call initialize(instance2, config2)

! Run in parallel or sequentially
call run_timestep(instance1)
call run_timestep(instance2)

Requirement: Zero module-level state (all in swap_state_t)
Test: tests/integration/test_multi_instance.f90
Build System

ALWAYS use pixi (see .continue/context/pixi-commands.md)
Key files:

    pixi.toml - Task definitions

    meson.build - Build configuration

    tests/meson.build - Test configuration

Never call directly: ninja, meson, gfortran, ifort
Always use: pixi run <task>