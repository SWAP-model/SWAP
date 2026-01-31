# SWAP State Refactoring Plan

## Objective

Replace scattered SAVE/DATA patterns with explicit state types to enable:
- Multi-instance execution (multiple SWAP models in one process)
- Thread-safe parallelization
- Future GPU offloading capability
- BMI-compatible state management

## Guiding Principles

1. **Backward Compatibility is Non-Negotiable** - External interface unchanged
2. **Incremental Migration** - One module at a time, tests pass at each step
3. **Separation of Concerns** - I/O handles separate from simulation state
4. **DATA Statement Strategy** - Constants → `parameter`, tunable → `*_config_t`
5. **Performance Optimization Deferred** - Focus on correctness first

---

## State Architecture

### Top-Level State Type

```fortran
type :: swap_state_t
    type(time_state_t)       :: time        ! Time stepping state
    type(soil_state_t)       :: soil        ! Soil water flow state
    type(atmosphere_state_t) :: atm         ! Meteorological state
    type(crop_state_t)       :: crop        ! Vegetation state
    type(drainage_state_t)   :: drain       ! Lateral drainage state
    type(boundary_state_t)   :: boundary    ! Boundary conditions state
    type(macropore_state_t)  :: macro       ! Macropore flow state
    type(solute_state_t)     :: solute      ! Solute transport state
    type(heat_state_t)       :: heat        ! Heat flow state
end type swap_state_t
```

### Separate I/O Handles

```fortran
type :: io_handles_t
    integer :: swp_unit = -1    ! Main input file
    integer :: met_unit = -1    ! Meteorology file
    integer :: crp_unit = -1    ! Crop file
    integer :: log_unit = -1    ! Log output
    integer :: csv_unit = -1    ! CSV output
    ! ... additional file handles
end type io_handles_t
```

### DATA Statement Handling

| Pattern | Strategy | Example |
|---------|----------|---------|
| Mathematical constants | `parameter` | `real, parameter :: pi = 3.14159...` |
| Gauss quadrature weights | `parameter` arrays | Fixed integration weights |
| Tunable algorithm params | `*_config_t` types | Tolerances, iteration limits |
| Lookup tables | Lazy initialization | Computed once, stored in state |

---

## Implementation Steps

### Step 1: Create State Module Foundation ✅ COMPLETED

**File:** `src/core/swap_state_mod.f90`

**Actions:**
- Define all `*_state_t` types (initially empty or with key fields)
- Define `swap_state_t` container
- Define `io_handles_t` for file handles
- Add initialization procedures

**Validation:** Module compiles, can be `use`d without errors

**Result:** Created comprehensive state module with:
- 12 sub-state types covering all model domains
- `swap_state_t` top-level container
- `io_handles_t` for file handle management
- Initialization procedures for all sub-states
- All arrays use allocatable for flexibility

---

### Step 2: Refactor Core Module + I/O Handles

**Files:** `src/core/*.f90`, `src/io/*.f90`

**Actions:**
- Move time-stepping variables from `variables.f90` to `time_state_t`
- Create `io_handles_t` and refactor file unit variables
- Update `swap.f90` main program to instantiate state

**Validation:** Build succeeds, hupselbrook test passes

---

### Step 3: Pilot - Soil Module (Highest SAVE Density)

**Files:** `src/soil/*.f90`

**Actions:**
- Define `soil_state_t` with hydraulic arrays, iteration counters
- Refactor `watfd.f90`, `soilwater.f90` to accept state argument
- Move SAVE variables (convergence tracking, previous timestep values)

**Validation:** Soil physics tests pass, water balance correct

---

### Step 4: Atmosphere Module

**Files:** `src/atmosphere/*.f90`

**Actions:**
- Define `atmosphere_state_t` (met data, snow state, ET accumulators)
- Refactor `meteo.f90`, `penman.f90`, `snow.f90`
- Handle met file reading state

**Validation:** ET calculations match reference

---

### Step 5: Crop Module (Largest, Most Complex)

**Files:** `src/crop/*.f90`

**Actions:**
- Define `crop_state_t` with WOFOST state, irrigation state
- Refactor 16 crop files incrementally
- Special attention to `cropd.f90`, `grass.f90` (heavy SAVE usage)

**Validation:** Crop growth trajectories match reference

---

### Step 6: Drainage Module

**Files:** `src/drainage/*.f90`

**Actions:**
- Define `drainage_state_t`
- Refactor lateral drainage and surface water routines

**Validation:** Drainage fluxes correct

---

### Step 7: Boundary Conditions Module

**Files:** `src/boundary/*.f90`

**Actions:**
- Define `boundary_state_t` (top/bottom BC values, groundwater state)
- Refactor `bbcgwl.f90`, `bctop.f90`

**Validation:** Boundary flux calculations correct

---

### Step 8: Macropore Module

**Files:** `src/macropore/*.f90`

**Actions:**
- Define `macropore_state_t`
- Refactor preferential flow routines

**Validation:** Macropore test case passes

---

### Step 9: Solute Module

**Files:** `src/solute/*.f90`

**Actions:**
- Define `solute_state_t` (includes pre-computed dispersion arrays)
- Handle lazy initialization of expensive arrays

**Validation:** Solute transport results match

---

### Step 10: Heat Module

**Files:** `src/heat/*.f90`

**Actions:**
- Define `heat_state_t`
- Refactor temperature and frost calculations

**Validation:** Soil temperature profiles correct

---

### Step 11: Integration - Full Model Wiring

**Actions:**
- Wire all sub-states through main time loop
- Ensure state flows correctly through all process calls
- Remove remaining module-level SAVE from `variables.f90`

**Validation:** All regression tests pass

---

### Step 12: Legacy Wrapper for Backward Compatibility

**File:** `src/core/swap_legacy.f90`

**Actions:**
- Create wrapper module with `save :: legacy_state`
- Provide legacy entry points that use module-level state
- External interface completely unchanged

```fortran
module swap_legacy
    use swap_state_mod
    implicit none
    private
    
    type(swap_state_t), save :: legacy_state
    type(io_handles_t), save :: legacy_io
    
    public :: swap_run  ! Legacy entry point
contains
    subroutine swap_run(swp_file)
        character(len=*), intent(in) :: swp_file
        call swap_run_with_state(legacy_state, legacy_io, swp_file)
    end subroutine
end module
```

**Validation:** Original executable behavior identical

---

### Step 13: Multi-Instance Validation

**Actions:**
- Write test running 2+ SWAP instances in same process
- Verify no state leakage between instances
- Document thread-safety guarantees

**Validation:** Parallel instances produce identical results to sequential

---

### Step 14: Documentation & Cleanup

**Actions:**
- Update developer documentation
- Remove dead code paths
- Document new state architecture
- Update BMI interface to use new state types

---

## Risk Mitigation

| Risk | Mitigation |
|------|------------|
| Subtle state bugs | Regression tests at every step |
| Performance regression | Profile before/after, optimize later |
| Large merge conflicts | Complete in focused sprint |
| Missed SAVE variables | grep audit + runtime testing |

---

## Future Optimization Phase (Deferred)

After correctness is verified:

1. **Standalone Procedures** - Extract hot loops to pure procedures for inlining
2. **GPU Directives** - Add OpenACC/OpenMP target pragmas to soil solver
3. **Memory Layout** - Optimize array ordering for vectorization
4. **Reduced Allocation** - Pre-allocate work arrays in state types

---

## Success Criteria

- [ ] All regression tests pass
- [ ] Multiple instances can run in same process
- [ ] No global/module-level SAVE except in legacy wrapper
- [ ] BMI interface uses new state types
- [ ] External interface 100% backward compatible
