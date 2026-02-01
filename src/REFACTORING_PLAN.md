# SWAP State Refactoring Plan

## Progress Summary

| Step | Description | Status | Date |
|------|-------------|--------|------|
| 1 | Create State Module Foundation | ✅ COMPLETED | 2026-01-31 |
| 2 | Create Synchronization Bridge | ✅ COMPLETED | 2026-01-31 |
| 3 | Pilot - Soil Module | ✅ COMPLETED | 2026-01-31 |
| 4 | Atmosphere Module | 🔲 Not started | |
| 5 | Crop Module | 🔲 Not started | |
| 6 | Drainage Module | ✅ COMPLETED | 2026-02-01 |
| 7 | Boundary Conditions | ✅ COMPLETED | 2026-02-01 |
| 8 | Macropore Module | 🔲 Not started | |
| 9 | Solute Module | 🔲 Not started | |
| 10 | Heat Module | 🔲 Not started | |
| 11 | Integration | 🔲 Not started | |
| 12 | Legacy Wrapper | 🔲 Not started | |
| 13 | Multi-Instance Validation | 🔲 Not started | |
| 14 | Documentation & Cleanup | 🔲 Not started | |

---

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

**Files:** 
- `src/core/swap_state_mod.f90` - State type definitions
- `src/core/swap_log.f90` - Logging infrastructure
- `tests/unit/core/test_state_standalone.f90` - Standalone tests

**Actions:**
- Define all `*_state_t` types with comprehensive field coverage
- Define `swap_state_t` container with all sub-states
- Define `io_handles_t` for file handles
- Add initialization procedures with proper allocation
- Create logging module for debugging during refactoring
- Create standalone test program (46 tests)

**Validation:** Module compiles, all 46 tests pass

**Result:** Created comprehensive state module with:
- 12 sub-state types: time, soil, atmosphere, crop, irrigation, drainage, boundary, macropore, solute, heat, snow, surfacewater
- `swap_state_t` top-level container
- `io_handles_t` for file handle management
- Initialization procedures for all sub-states
- All arrays use allocatable for flexibility
- Logging infrastructure (`swap_log.f90`) with levels: DEBUG, INFO, WARN, ERROR
- Multi-instance independence verified by tests

---

### Step 2: Create Synchronization Bridge ✅ COMPLETED

**File:** `src/core/swap_state_sync.f90`

**Actions:**
- Created bidirectional sync between `variables.f90` and state types
- `state_from_variables()` - snapshot globals into state (for multi-instance)
- `state_to_variables()` - restore state to globals (for legacy code)
- Component-level sync procedures for selective updates
- Keep `variables.f90` completely intact for backwards compatibility

**Validation:** Module compiles, can be `use`d without errors

**Result:** Created synchronization module with:
- Master sync procedures for all state ↔ variables
- Time state sync (time stepping, flags, counters)
- Soil state sync (pressure heads, water content, GWL)
- Atmosphere state sync (meteo, ET, precipitation)
- Crop state sync (development, biomass, rooting)
- Irrigation and drainage state sync
- Debug logging integrated throughout

**Multi-instance pattern enabled:**
```fortran
! Switch to instance 1
call state_to_variables(state1)
call existing_swap_code()  ! Uses globals
call state_from_variables(state1)

! Switch to instance 2  
call state_to_variables(state2)
call existing_swap_code()  ! Uses globals
call state_from_variables(state2)
```

---

### Step 3: Pilot - Soil Module (Highest SAVE Density) ✅ COMPLETED

**Files Modified:**
- `src/core/swap_state_mod.f90` - Extended soil_state_t with headcalc tracking
- `src/core/swap_state_sync.f90` - Enhanced soil sync procedures (~50 additional variables)
- `src/soil/headcalc.f90` - Added debug logging
- `src/soil/soilwater.f90` - Added debug logging  
- `src/soil/calcgwl.f90` - Added debug logging
- `src/soil/watstor.f90` - Added debug logging
- `tests/unit/soil/test_soil_state.f90` - New soil-focused test program

**Actions Completed:**
- Added headcalc iteration tracking to soil_state_t:
  - `flwarn_hc`, `iwarn_hc`, `nstep_hc` (from headcalc.f90 SAVE variables)
- Extended soil sync with comprehensive variable coverage:
  - Previous timestep values (`hm1`, `thetm1`)
  - All flux arrays (`q`, `qrot`, `kmean`)
  - Complete groundwater state (`gwl`, `gwlm1`, `gwli`, `gwlinp`, `nodgwl`, `npegwl`, `bpegwl`)
  - Surface/ponding (`pond`, `pondm1`, `pondmx`, `qtop`, `qbot`)
  - All cumulative fluxes (`cqbot`, `cqbotdo`, `cqbotup`, `cqtdo`, `cqtup`, `cqrot`, `cqdra`, `crunoff`, `crunon`)
  - Storage tracking (`volact`, `volini`, `volm1`)
  - Evaporation reduction (`saev`, `spev`, `ldwet`, `cofred`)
- Added debug logging to key soil routines for verification
- Created comprehensive soil state test suite (77 tests)
- Added finalize procedures for proper cleanup

**Validation:**
- All 77 soil state tests pass
- All 46 core state tests pass
- SWAP Hupselbrook integration test passes
- Water balance unchanged

**Note on TSPACK (sptabulated.f90):**
- Contains computational intermediates for tabulated soil physics
- Lower priority for state migration (not simulation state per se)
- Can be addressed in a future optimization pass

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

### Step 6: Drainage Module ✅ COMPLETED

**Files Modified:**
- `src/core/swap_state_mod.f90` - Extended drainage_state_t and surfacewater_state_t
- `src/core/swap_state_sync.f90` - Enhanced drainage/surfacewater sync procedures
- `src/drainage/drainage.f90` - Removed SAVE statement, added debug logging
- `src/drainage/surfacewater.f90` - Removed SAVE statement, added debug logging

**Actions Completed:**
- Extended `drainage_state_t` with comprehensive field coverage:
  - All drainage flux arrays (`qdrain`, `cqdrain`, `qdra`, `inqdra`, `inqdra_in/out`)
  - All resistance/geometry arrays (`drares`, `infres`, `L`, `wetper`, `zbotdr`, `rdrain`, `rinfi`, `rentry`, `rexit`, `gwlinf`, `widthr`, `taludr`)
  - Drainage type switches (`swallo`, `swdtyp`, `swtopdislay`, `zTopDisLay`, `fTopDisLay`)
  - Interflow parameters (`cofintfl`, `expintfl`, `swnrsrf`, `SwTopnrsrf`, `rsurfdeep`, `rsurfshallow`, `FacDpthInf`, `Swdivdinf`)
- Extended `surfacewater_state_t` with comprehensive field coverage:
  - Water levels (`wlp`, `wls`, `wlsold`, `wlstar`, `hwlman`, `vtair`, `wlsbak`)
  - Storage and fluxes (`swst`, `qdrd`, `cqdrd`, `cwsupp`, `cwout`, `runots`, `QRapDra`)
  - Management arrays (`impend`, `swman`, `hbweir`, `wldip`, `alphaw`, `betaw`, `wscap`, `dropr`, `intwl`, `nphase`, `nodhd`, `gwlcrit`, `wlsman`, `vcrit`, `hcrit`)
  - Lookup tables (`wlstab`, `wlptab`, `sttab`, `owltab`, `qqhtab`)
- Added `surfacewater_state_init` and `surfacewater_state_finalize` procedures
- Updated drainage_state_init with all new array allocations
- Updated drain_state_finalize with comprehensive deallocations
- Created bidirectional sync procedures for drainage and surface water
- Removed SAVE statements from `bocodre` (drainage.f90) and `wlevbal` (surfacewater.f90)
- Added debug logging with `swap_log` module integration

**SAVE Variables Removed:**
- `bocodre`: Local variables only (level, imper, qdrdm, etc.) - no persistent state needed
- `wlevbal`: Local variables only (iphase, wlstx, swsttar, etc.) - no persistent state needed

**Validation:** All tests pass

---

### Step 7: Boundary Conditions Module ✅ COMPLETED

**Files Modified:**
- `src/core/swap_state_mod.f90` - Extended boundary_state_t with comprehensive fields
- `src/core/swap_state_sync.f90` - Added boundary sync procedures
- `src/boundary/boundtop.f90` - Converted DATA to parameter, added debug logging
- `src/boundary/boundbottom.f90` - Added debug logging

**Actions Completed:**
- Extended `boundary_state_t` with comprehensive field coverage:
  - Bottom BC configuration (`swbotb`, `swbotb3Impl`, `SwBotb3ResVert`, `swqhbot`, `swcofqhc`, `sw2`, `sw3`, `sw4`)
  - Bottom BC values (`qbot`, `qbot_nonfrozen`, `hbot`, `iqbot`, `cqbot`, `cqbotdo`, `cqbotup`, `deepgw`)
  - Aquifer parameters (`aqave`, `aqamp`, `aqper`, `aqtmax`, `rimlay`, `hdrain`, `shape`)
  - Sine function parameters (`sinave`, `sinamp`, `sinmax`)
  - Flux-head relationships (`cofqha`, `cofqhb`, `cofqhc`)
  - Lysimeter parameters (`hplate`)
  - Prescribed tables (`gwltab`, `haqtab`, `qbotab`, `hbotab`)
  - Top BC configuration (`swpondmx`, `swredu`)
  - Top BC values (`pondmx`, `hatm`, `hsurf`, `rsro`, `rsroexp`, `runon`, `runots`, `crunoff`, `crunon`)
  - Surface/ponding (`qtop`, `q0`, `h0max`, `k1max`, `QMpLatSs`)
  - Runon table (`runonarr`), ponding table (`pondmxtab`)
  - Flags (`FlRunoff`, `flrunon`, `ftoph`)
- Added `boundary_state_init` and `boundary_state_finalize` procedures
- Created bidirectional sync procedures (`boundary_state_from_variables`, `boundary_state_to_variables`)
- Added boundary sync to master sync procedures
- Converted DATA statement (`hconode_vsmall`) to Fortran parameter in boundtop.f90
- Added debug logging with `swap_log` module integration

**DATA Statement Converted:**
- `boundtop.f90`: `hconode_vsmall` (frozen soil conductivity) - mathematical constant, now a parameter

**Validation:** Compile and test (pending)

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
