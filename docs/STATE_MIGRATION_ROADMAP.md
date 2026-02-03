# SWAP State Migration Roadmap

A visual guide to modernizing the SWAP model architecture for multi-instance execution and external coupling.

---

## The Goal

Transform SWAP from a **single-instance** model with global variables to a **multi-instance** model where each instance has its own encapsulated state.

```
BEFORE: One SWAP, One State          AFTER: Many SWAPs, Each with Own State
┌─────────────────────┐              ┌─────────┐ ┌─────────┐ ┌─────────┐
│       SWAP          │              │  SWAP   │ │  SWAP   │ │  SWAP   │
│  ┌───────────────┐  │              │ ┌─────┐ │ │ ┌─────┐ │ │ ┌─────┐ │
│  │ Global State  │  │      ──▶     │ │State│ │ │ │State│ │ │ │State│ │
│  │ (variables.f90)  │              │ └─────┘ │ │ └─────┘ │ │ └─────┘ │
│  └───────────────┘  │              └─────────┘ └─────────┘ └─────────┘
└─────────────────────┘                Cell 1      Cell 2      Cell 3
```

---

## Why This Matters

| Capability | Before | After |
|------------|--------|-------|
| Run multiple SWAP instances in parallel | ❌ | ✅ |
| Couple with MODFLOW 6 (thousands of cells) | ❌ | ✅ |
| Thread-safe parallelization | ❌ | ✅ |
| BMI-compatible state access | ❌ | ✅ |
| GPU offloading (future) | ❌ | ✅ |

---

## The Migration Path

We're taking a **safe, incremental approach** with 4 phases:

```
┌──────────────────────────────────────────────────────────────────────────┐
│                         MIGRATION TIMELINE                                │
├────────────┬────────────┬────────────┬────────────┬─────────────────────┤
│  Phase 0   │  Phase 1   │  Phase 2   │  Phase 3   │      Phase 4        │
│  (done)    │  (done)    │  (current) │  (next)    │      (future)       │
├────────────┼────────────┼────────────┼────────────┼─────────────────────┤
│  Define    │   Create   │  Add Sync  │   Use      │  Remove Globals     │
│  State     │   State    │  Bridge    │   State    │  & Sync Bridge      │
│  Types     │   Init     │            │  Directly  │                     │
└────────────┴────────────┴────────────┴────────────┴─────────────────────┘
```

---

## Phase 0: Define State Types ✅

**What**: Create derived types that mirror all module variables.

```fortran
! Before: Variables scattered in module
module variables
    integer :: daynr, iyear
    real(8) :: h(MACP), theta(MACP)
    ! ... hundreds more ...
end module

! After: Organized state types
type :: time_state_t
    integer :: daynr, iyear
end type

type :: soil_state_t
    real(8), allocatable :: h(:), theta(:)
end type

type :: swap_state_t
    type(time_state_t) :: time
    type(soil_state_t) :: soil
    ! ... all sub-states ...
end type
```

---

## Phase 1: Create State Initialization ✅

**What**: Add procedures to allocate and initialize state containers.

```fortran
subroutine swap_state_init(state, numnod, numlay)
    type(swap_state_t), intent(inout) :: state
    
    ! Allocate arrays based on grid size
    allocate(state%soil%h(numnod))
    allocate(state%soil%theta(numnod))
    ! ...
end subroutine
```

---

## Phase 2: Add Sync Bridge ← WE ARE HERE

**What**: Add synchronization between module variables and state types.

```
┌─────────────────────────────────────────────────────────────────────────┐
│                          SYNC BRIDGE PATTERN                             │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                          │
│   ┌──────────────────┐                      ┌──────────────────┐        │
│   │ Module Variables │                      │   State Types    │        │
│   │ (source of truth)│                      │    (mirror)      │        │
│   │                  │                      │                  │        │
│   │  daynr = 42      │  ───── exit ─────▶   │  daynr = 42      │        │
│   │  iyear = 2024    │       sync           │  iyear = 2024    │        │
│   │  h(1:40)         │                      │  h(1:40)         │        │
│   └──────────────────┘                      └──────────────────┘        │
│                                                                          │
│   Computation happens here                  External access here         │
│   (legacy code unchanged)                   (BMI, coupling)              │
│                                                                          │
└─────────────────────────────────────────────────────────────────────────┘
```

**How it works**:
1. Subroutines still read/write module variables (no code changes needed)
2. At subroutine exit, sync captures values into state
3. State is available for external access (BMI, coupling)

```fortran
subroutine TimeControl(task, tstate)
    ! ... existing code uses module variables ...
    
    ! At exit: capture into state
    if (present(tstate)) call time_state_from_variables(tstate)
end subroutine
```

---

## Phase 3: Use State Directly (Future)

**What**: Modify subroutines to read/write state directly instead of module variables.

```
┌─────────────────────────────────────────────────────────────────────────┐
│                      DIRECT STATE ACCESS                                 │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                          │
│   ┌──────────────────┐                      ┌──────────────────┐        │
│   │ Module Variables │                      │   State Types    │        │
│   │  (deprecated)    │  ◀── entry ───       │ (source of truth)│        │
│   │                  │      sync            │                  │        │
│   │  daynr           │                      │  daynr = 42      │        │
│   │  iyear           │                      │  iyear = 2024    │        │
│   └──────────────────┘                      └──────────────────┘        │
│         ▲                                          │                     │
│         │                                          │                     │
│         └──────── backward compatibility ──────────┘                     │
│                   (old code still works)                                 │
│                                                                          │
└─────────────────────────────────────────────────────────────────────────┘
```

```fortran
subroutine TimeControl(task, tstate)
    ! New code uses state directly
    tstate%daynr = tstate%daynr + 1
    
    ! Entry sync for backward compatibility with old code
    if (present(tstate)) call time_state_to_variables(tstate)
end subroutine
```

---

## Phase 4: Remove Globals (Future)

**What**: Delete module variables and sync procedures. State types are the only storage.

```
┌─────────────────────────────────────────────────────────────────────────┐
│                        FINAL ARCHITECTURE                                │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                          │
│                              ┌──────────────────┐                        │
│                              │   State Types    │                        │
│                              │ (only storage)   │                        │
│                              │                  │                        │
│                              │  daynr = 42      │                        │
│                              │  iyear = 2024    │                        │
│                              │  h(1:40)         │                        │
│                              └──────────────────┘                        │
│                                      │                                   │
│                    ┌─────────────────┼─────────────────┐                │
│                    ▼                 ▼                 ▼                │
│             ┌──────────┐      ┌──────────┐      ┌──────────┐            │
│             │TimeControl     │SoilWater │      │ Drainage │            │
│             │(tstate)  │      │(sstate)  │      │(dstate)  │            │
│             └──────────┘      └──────────┘      └──────────┘            │
│                                                                          │
│   All subroutines receive state as argument, no global variables        │
│                                                                          │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## Multi-Instance Execution

Once Phase 4 is complete, running multiple SWAP instances is straightforward:

```fortran
program swap_mf6_coupler
    type(swap_state_t) :: states(1000)  ! One state per MODFLOW cell
    
    ! Initialize all instances
    do i = 1, 1000
        call swap_state_init(states(i), numnod, numlay)
    end do
    
    ! Run all instances (can be parallelized!)
    !$omp parallel do
    do i = 1, 1000
        call swap(iTask=2, state=states(i))
    end do
    !$omp end parallel do
end program
```

---

## Summary

| Phase | Status | Key Change |
|-------|--------|------------|
| 0. Define Types | ✅ Done | `swap_state_t` and sub-types defined |
| 1. State Init | ✅ Done | `swap_state_init` allocates arrays |
| 2. Sync Bridge | 🔄 Current | Exit-only sync captures state |
| 3. Direct Access | ⏳ Future | Subroutines use state directly |
| 4. Remove Globals | ⏳ Future | Clean architecture, multi-instance ready |

**Current milestone**: TimeControl integrated with sync bridge, all tests passing.

---

## References

- [REFACTORING_PLAN.md](REFACTORING_PLAN.md) - Detailed implementation plan
- [REFACTORING_PITFALLS.md](REFACTORING_PITFALLS.md) - Common mistakes and solutions
- `swap_state_mod.f90` - State type definitions
- `swap_state_sync.f90` - Synchronization procedures
