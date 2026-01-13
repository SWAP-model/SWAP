# Understanding SAVE Variables vs State Objects in Fortran

## The Problem with SAVE Variables

### What are SAVE Variables?

In Fortran, `SAVE` is an attribute that tells the compiler to preserve a variable's value between function/subroutine calls. It's like a "static" variable in C/C++.

```fortran
! Example from SWAP's irrigation.f90
subroutine irrigation_calculation()
    integer, save :: nirri              ! Counter preserved between calls
    real(8), save :: ssdi_amount        ! Configuration value
    real(8), save :: days_counter       ! State that accumulates
    
    ! This variable is LOCAL to this subroutine but PERSISTS
    ! across multiple calls
    days_counter = days_counter + 1
end subroutine
```

### How SWAP Uses SAVE Variables

Looking at the code, SWAP has SAVE variables in multiple places:

#### 1. Module-Level SAVE (variables.f90)
```fortran
module variables
    implicit none
    save  ! <-- This line makes ALL variables in the module SAVE
    
    ! Time variables
    real(8) :: t, tcum, dt
    integer :: daynr, iyear
    
    ! State arrays
    real(8) :: theta(maxnod), h(maxnod)
    
    ! Configuration
    character(len=80) :: project
end module
```

**This is the biggest issue**: The entire `variables` module uses `save`, meaning ALL 1000+ variables are global state.

#### 2. Subroutine-Level SAVE
```fortran
! From headcalc.f90
subroutine headcalc()
    logical, save :: flwarn
    integer, save :: iwarn, nstep
    
    ! These variables persist across calls to headcalc()
    ! They track state between time steps
end subroutine
```

#### 3. In Other Modules
```fortran
! From irrigation.f90
integer, save :: swssdi             ! Is subsurface drip active?
integer, save :: days_counter       ! Days since last irrigation
real(8), save :: ssdi_threshold     ! Threshold value
```

---

## Why SAVE Variables Break Multiprocessing

### The Core Problem: Global State

SAVE variables are stored in **static memory**, which means:

```
Process Memory Layout:
┌─────────────────┐
│ Text (code)     │  ← Shared across all instances
├─────────────────┤
│ Static Data     │  ← SAVE variables live here! 
│  - variables    │     ONE COPY PER PROCESS
│  - nirri        │
│  - days_counter │
├─────────────────┤
│ Heap            │  ← Dynamic allocations
├─────────────────┤
│ Stack           │  ← Local variables
└─────────────────┘
```

### Scenario: Running 2 SWAP Instances in Same Process

```python
# What we WANT to do:
swap1 = SwapModel()  # Simulate field A
swap2 = SwapModel()  # Simulate field B

swap1.initialize(grid_size=100, year=2020)
swap2.initialize(grid_size=50,  year=2021)

# What ACTUALLY happens with SAVE variables:
# Both swap1 and swap2 access THE SAME memory location!
```

**Concrete Example**:
```fortran
! Inside variables module (simplified)
module variables
    save
    integer :: NumGrid     ! Grid size
    real(8) :: theta(500)  ! Water content array
end module

! When you do:
! swap1.set_grid(NumGrid=100) → variables::NumGrid = 100
! swap2.set_grid(NumGrid=50)  → variables::NumGrid = 50  (overwrites!)
!
! Now swap1 thinks it has 50 nodes instead of 100!
```

### Why This Fails Catastrophically

1. **Configuration Overwrite**:
   ```python
   swap1.set_parameters(theta_s=0.45)  # Set for field A
   swap2.set_parameters(theta_s=0.38)  # Overwrites field A's value!
   ```

2. **State Corruption**:
   ```python
   swap1.update()  # Advances time for swap1
   # But shared variables::t, tcum, daynr are changed
   swap2.update()  # Now using wrong time values!
   ```

3. **Array Dimension Mismatch**:
   ```python
   swap1.initialize(grid=100)  # Allocates theta(100)
   swap2.initialize(grid=200)  # Reallocates theta(200) - swap1's data lost!
   ```

---

## The State Object Solution

### Concept: Explicit State Passing

Instead of using global SAVE variables, pass state explicitly through function calls:

```fortran
! BEFORE (with SAVE):
subroutine calculate_flux()
    use variables  ! Accesses global state
    
    ! Uses global theta, h, K arrays
    flux = K(i) * (h(i) - h(i+1)) / dz
end subroutine

! AFTER (with state object):
subroutine calculate_flux(state)
    type(swap_state_t), intent(inout) :: state
    
    ! Uses state's local arrays
    flux = state%K(i) * (state%h(i) - state%h(i+1)) / state%dz
end subroutine
```

### Defining a State Object

```fortran
! NEW FILE: swap_state.f90
module swap_state
    implicit none
    
    ! Define the state type that holds ALL model state
    type :: swap_state_t
        ! === Grid Configuration ===
        integer :: NumGrid              ! Number of nodes
        real(8), allocatable :: z(:)    ! Node depths (NumGrid)
        real(8), allocatable :: dz(:)   ! Layer thicknesses (NumGrid-1)
        
        ! === Time State ===
        real(8) :: t                    ! Time since year start
        real(8) :: tcum                 ! Time since simulation start
        real(8) :: dt                   ! Current time step
        integer :: daynr                ! Day of year
        integer :: iyear                ! Current year
        
        ! === Soil State ===
        real(8), allocatable :: theta(:)  ! Water content (NumGrid)
        real(8), allocatable :: h(:)      ! Pressure head (NumGrid)
        real(8), allocatable :: K(:)      ! Hydraulic conductivity (NumGrid)
        
        ! === Soil Parameters ===
        real(8), allocatable :: thetaS(:)  ! Saturated water content
        real(8), allocatable :: thetaR(:)  ! Residual water content
        real(8), allocatable :: Ksat(:)    ! Saturated K
        real(8), allocatable :: alpha(:)   ! van Genuchten alpha
        real(8), allocatable :: n(:)       ! van Genuchten n
        
        ! === Crop State ===
        real(8) :: LAI                  ! Leaf area index
        real(8) :: root_depth           ! Rooting depth
        real(8) :: crop_height          ! Crop height
        
        ! === Fluxes (cumulative) ===
        real(8) :: cevap                ! Cumulative evaporation
        real(8) :: cptra                ! Cumulative transpiration
        real(8) :: cgrai                ! Cumulative precipitation
        
        ! === Irrigation State ===
        integer :: nirri                ! Irrigation counter
        integer :: days_counter         ! Days since last irrigation
        real(8) :: ssdi_amount          ! Irrigation amount
        
        ! === Configuration (not changed during run) ===
        character(len=80) :: project    ! Project name
        character(len=80) :: pathwork   ! Work directory
        
        ! ... hundreds more variables ...
    end type swap_state_t
    
contains
    
    ! Constructor
    function swap_state_create(NumGrid) result(state)
        integer, intent(in) :: NumGrid
        type(swap_state_t) :: state
        
        state%NumGrid = NumGrid
        
        ! Allocate arrays
        allocate(state%z(NumGrid))
        allocate(state%theta(NumGrid))
        allocate(state%h(NumGrid))
        allocate(state%K(NumGrid))
        allocate(state%thetaS(NumGrid))
        ! ... etc
        
        ! Initialize to zero
        state%z = 0.0d0
        state%theta = 0.0d0
        ! ... etc
    end function
    
    ! Destructor
    subroutine swap_state_destroy(state)
        type(swap_state_t), intent(inout) :: state
        
        if (allocated(state%z)) deallocate(state%z)
        if (allocated(state%theta)) deallocate(state%theta)
        ! ... etc
    end subroutine
    
end module swap_state
```

### Refactored SWAP Routines

```fortran
! OLD: Using global variables
subroutine headcalc()
    use variables  ! Global state
    
    ! Access global arrays directly
    do i = 1, NumGrid
        h(i) = calculate_new_head(theta(i), ...)
    end do
end subroutine

! NEW: Using state object
subroutine headcalc(state)
    type(swap_state_t), intent(inout) :: state
    
    ! Access state's arrays
    do i = 1, state%NumGrid
        state%h(i) = calculate_new_head(state%theta(i), ...)
    end do
end subroutine
```

### How It Enables Multiple Instances

```fortran
! NEW BMI INTERFACE with state objects
module swap_bmi_multi
    implicit none
    
    ! Store multiple SWAP states
    type(swap_state_t), allocatable :: states(:)
    integer :: num_instances = 0
    
contains

    function bmi_create_instance(config_file) result(instance_id) bind(C)
        character(kind=c_char), intent(in) :: config_file(*)
        integer(c_int) :: instance_id
        
        ! Create a NEW state object
        num_instances = num_instances + 1
        
        ! Reallocate array to hold more states
        ! (In real code, use better data structure)
        allocate(states(num_instances))
        
        ! Initialize THIS instance's state
        states(num_instances) = swap_state_create(NumGrid=100)
        
        instance_id = num_instances
    end function
    
    function bmi_update(instance_id) result(status) bind(C)
        integer(c_int), intent(in) :: instance_id
        integer(c_int) :: status
        
        ! Update ONLY this instance's state
        call swap_timestep(states(instance_id))
        
        status = BMI_SUCCESS
    end function
    
end module
```

### Python Usage with Multiple Instances

```python
# NOW POSSIBLE: Multiple independent instances
swap1 = SwapBMI()
swap2 = SwapBMI()

# Each gets unique instance_id
id1 = swap1.create_instance("field_a.swp")  # instance_id = 1
id2 = swap2.create_instance("field_b.swp")  # instance_id = 2

# Independent state!
swap1.set_grid_size(id1, 100)
swap2.set_grid_size(id2, 200)  # Doesn't affect swap1

# Independent simulation
for t in range(365):
    swap1.update(id1)  # Updates states[1]
    swap2.update(id2)  # Updates states[2]

# Get independent results
theta1 = swap1.get_theta(id1)  # From states[1]
theta2 = swap2.get_theta(id2)  # From states[2]
```

---

## What the Refactoring Entails

### Step-by-Step Process

#### Phase 1: Create State Type (2-3 weeks)
1. **Inventory all SAVE variables** (grep through codebase)
2. **Group by logical categories**:
   - Grid configuration
   - Time tracking
   - State variables (theta, h, K)
   - Parameters (soil, crop)
   - Cumulative fluxes
   - Output tracking
3. **Define `swap_state_t` type** with all fields
4. **Create constructor/destructor** functions

#### Phase 2: Module-by-Module Refactoring (2-3 months)
Start with leaf modules (no dependencies), work up:

```
Order of refactoring:
1. variables.f90        → Move to swap_state module
2. swap_exchange.f90    → Add state parameter
3. watstor.f90          → Add state parameter  
4. fluxes.f90           → Add state parameter
5. headcalc.f90         → Remove local SAVE, use state
6. irrigation.f90       → Remove local SAVE, use state
7. crop modules         → Add state parameter
... continue for ~50 modules
```

For each module:
```fortran
! Example: refactoring watstor.f90

! STEP 1: Add state parameter to all subroutines
subroutine watstor(state)
    type(swap_state_t), intent(inout) :: state
    
! STEP 2: Replace variable access
    ! OLD: use variables; theta(i) = ...
    ! NEW: state%theta(i) = ...
    
! STEP 3: Update all function calls
    ! OLD: call fluxes(...)
    ! NEW: call fluxes(state, ...)
end subroutine
```

#### Phase 3: Remove MODULE-Level SAVE (1 week)
```fortran
! OLD variables.f90
module variables
    implicit none
    save  ! <-- DELETE THIS LINE
    
    ! Delete all variable declarations
    ! (now in swap_state_t)
end module

! This module becomes empty or just contains
! constants/parameters
```

#### Phase 4: Update BMI Interface (1-2 weeks)
```fortran
! Modify swap_bmi.f90 to manage state objects
! Add instance creation/destruction
! Thread state through all BMI functions
```

#### Phase 5: Testing (2-4 weeks)
- Unit test each refactored module
- Integration testing
- Compare results with original SWAP
- Multi-instance testing
- Memory leak checking

---

## Concrete Example: Before and After

### BEFORE (Current SWAP)

```fortran
! ===== variables.f90 =====
module variables
    save
    integer :: NumGrid
    real(8) :: theta(500)
    real(8) :: h(500)
end module

! ===== headcalc.f90 =====
subroutine headcalc()
    use variables  ! Global state
    integer, save :: iteration_count = 0
    
    iteration_count = iteration_count + 1
    
    do i = 1, NumGrid
        h(i) = update_pressure(theta(i))
    end do
end subroutine

! ===== Python interface =====
! Can ONLY create ONE instance per process:
from swap import SwapBMI
swap = SwapBMI()  # Only one!
```

### AFTER (State Object Refactoring)

```fortran
! ===== swap_state.f90 =====
module swap_state
    type :: swap_state_t
        integer :: NumGrid
        real(8), allocatable :: theta(:)
        real(8), allocatable :: h(:)
        integer :: iteration_count  ! No longer SAVE!
    end type
    
contains
    function create(NumGrid) result(state)
        type(swap_state_t) :: state
        state%NumGrid = NumGrid
        allocate(state%theta(NumGrid))
        allocate(state%h(NumGrid))
        state%iteration_count = 0
    end function
end module

! ===== headcalc.f90 =====
subroutine headcalc(state)
    use swap_state
    type(swap_state_t), intent(inout) :: state
    
    state%iteration_count = state%iteration_count + 1
    
    do i = 1, state%NumGrid
        state%h(i) = update_pressure(state%theta(i))
    end do
end subroutine

! ===== swap_bmi_multi.f90 =====
module swap_bmi_multi
    use swap_state
    type :: bmi_instance_t
        integer :: id
        type(swap_state_t) :: state
    end type
    
    type(bmi_instance_t), allocatable :: instances(:)
contains
    function create_instance() bind(C)
        integer(c_int) :: id
        ! Create new state object
        ! Store in instances array
    end function
end module

! ===== Python interface =====
# Can create MANY instances in same process!
swap1 = SwapBMI()
swap2 = SwapBMI()
swap3 = SwapBMI()

# Each has independent state
id1 = swap1.create()  # state%NumGrid, %theta, etc.
id2 = swap2.create()  # DIFFERENT state%NumGrid, %theta, etc.
id3 = swap3.create()  # DIFFERENT state%NumGrid, %theta, etc.
```

---

## Key Differences: SAVE vs State Object

| Aspect | SAVE Variables | State Object |
|--------|----------------|--------------|
| **Memory Location** | Static (global) | Heap (per-instance) |
| **Scope** | Module/subroutine | Passed explicitly |
| **Instances** | ONE per process | MANY per process |
| **Thread Safety** | ❌ Not thread-safe | ✅ Thread-safe |
| **Multiprocessing** | ⚠️ Need separate processes | ✅ Can use threads too |
| **Memory Management** | Automatic | Manual (alloc/dealloc) |
| **Debugging** | Hard (hidden state) | Easy (explicit passing) |
| **Testing** | Hard (global state) | Easy (isolated state) |

---

## Estimated Effort Breakdown

### Development Time
- **State type definition**: 2-3 weeks (inventory + design)
- **Module refactoring**: 2-3 months (50+ modules, ~1-2 days each)
- **BMI interface updates**: 1-2 weeks
- **Testing infrastructure**: 2 weeks
- **Integration testing**: 2-4 weeks
- **Documentation**: 1-2 weeks

**Total**: 4-6 months of focused development

### Risk Factors
- **High**: Introducing bugs during refactoring
- **Medium**: Performance degradation from pointer passing
- **Medium**: Breaking existing SWAP workflows
- **Low**: Memory leaks (caught by testing)

### Lines of Code Changed
- **variables.f90**: 1133 lines → convert to type definition
- **~50 modules**: Add `state` parameter to ~200 subroutines
- **~1000 call sites**: Update to pass `state`
- **Estimated**: 5000-10000 lines touched

---

## Why It's Worth It

### Benefits
1. **True multi-instance**: Run 100s of SWAP instances in one process
2. **Thread safety**: Can use Python threading (not just multiprocessing)
3. **Memory efficiency**: Shared code, independent data
4. **Cleaner architecture**: Explicit dependencies
5. **Easier testing**: Isolated state per test
6. **Modern Fortran**: Aligns with Fortran 2003+ best practices

### Performance Considerations
```fortran
! Worry: Passing large state object is slow
! Reality: Pass by REFERENCE (pointer), not by value

subroutine calc(state)
    type(swap_state_t), intent(inout) :: state
    ! Only passes 8-byte pointer, not entire struct!
end subroutine
```

### Multiprocessing Comparison

```python
# CURRENT (with SAVE variables): Need separate processes
from multiprocessing import Process

def run_swap(config):
    swap = SwapBMI()  # Each process gets own memory space
    swap.initialize(config)
    swap.run()

if __name__ == '__main__':
    # Spawn 4 separate processes (high overhead)
    procs = [Process(target=run_swap, args=(cfg,)) 
             for cfg in configs]
    for p in procs: p.start()
    for p in procs: p.join()
```

```python
# AFTER (with state objects): Can use threads OR processes
from concurrent.futures import ThreadPoolExecutor

def run_swap(config):
    swap = SwapBMI()  # All threads share process memory
    instance_id = swap.create()  # But each has own state!
    swap.initialize(instance_id, config)
    return swap.run(instance_id)

# Use threads (much lower overhead than processes)
with ThreadPoolExecutor(max_workers=4) as executor:
    results = executor.map(run_swap, configs)
```

---

## Decision Time

The refactoring is **essential for true multi-instance capability**, but it's a **major undertaking**. 

Options:
1. **Dive in now**: 4-6 months to refactor everything
2. **Hybrid first**: Use Option 3 (temp files + multiprocessing) for 1-2 months while planning refactoring
3. **Incremental**: Refactor one module at a time, starting with most critical (variables, headcalc, irrigation)
4. **Minimal**: Just remove SAVE from a few critical subroutines, keep module-level SAVE for now

Which approach fits your timeline and risk tolerance?
