# SWAP State Integration Planning - Phase 2 Strategy

## Executive Summary

This plan outlines the phased integration of state-based computation for the SWAP model, transitioning from module-based global variables to explicit state passing. The goal is to eliminate sync calls from the main time loop while maintaining numerical accuracy and backward compatibility.

### Key Objectives
- Convert all computational subroutines to accept state arguments
- Eliminate module global dependencies in process routines
- Maintain water balance accuracy (<0.1% error)
- Preserve legacy input/output compatibility
- Prepare for Python BMI interface in Phase 3

### Total Estimated Effort
- **Phase 2 Integration**: 40-60 hours
- **Testing & Validation**: 15-20 hours
- **Documentation Updates**: 5-10 hours

### Key Risks & Mitigations
- **Physics Changes**: All changes must preserve numerical algorithms
- **Call Site Complexity**: Will require careful refactoring of 100+ call sites
- **Integration Testing**: Comprehensive regression testing required

## Detailed Milestone Plan

### Milestone 1: Core Time Loop and State Initialization
**Objective**: Establish state-based time loop framework and initialize state from legacy inputs

**Files to Modify**:
- `src/core/swap.f90` - Main time loop
- `src/core/swap_state_sync.f90` - Sync procedures
- `src/core/initialize.f90` - State initialization

**Subroutines to Refactor**:
- `swap()` - Main time loop driver
- `state_from_variables()` - Initialize state from Variables
- `state_to_variables()` - Final sync for output

**Call Sites to Update**:
- All process routines in main time loop

**Testing Requirements**:
```bash
pixi run test-unit-all
pixi run test-linux-hupselbrook
pixi run regression_all
```

**Success Criteria**:
- Main time loop operates on state object
- No sync calls inside time loop
- Water balance error < 0.1%

**Dependencies**: 
- State type definitions (completed)
- Sync layer implementation (completed)

**Effort Estimate**: 8 hours

### Milestone 2: Soil Water Subsystem Integration
**Objective**: Convert soil water routines to state-based interface

**Files to Modify**:
- `src/soil/soilwater.f90`
- `src/soil/watstor.f90`
- `src/soil/fluxes.f90`
- `src/soil/calcgwl.f90`
- `src/soil/headcalc.f90`
- `src/soil/checkmassbal.f90`

**Subroutines to Refactor**:
- `soilwater()` - Main soil water routine
- `watstor()` - Water storage calculations
- `fluxes()` - Flux calculations
- `calcgwl()` - Groundwater level calculations
- `headcalc()` - Head calculations
- `checkmassbal()` - Mass balance checks

**Call Sites to Update**:
- All calls from main time loop and parent routines

**Testing Requirements**:
```bash
pixi run test-linux-hupselbrook
pixi run test-linux-grassgrowth
pixi run regression_all
```

**Success Criteria**:
- All soil water routines accept `soil_state_t` and `boundary_state_t`
- No module dependencies in soil routines
- Water balance preserved

**Dependencies**: 
- Milestone 1 completed

**Effort Estimate**: 15 hours

### Milestone 3: Crop Growth Subsystem Integration
**Objective**: Convert crop growth routines to state-based interface

**Files to Modify**:
- `src/crop/cropgrowth.f90`
- `src/crop/wofostnut.f90`
- `src/crop/rootextraction.f90`
- `src/crop/irrigation.f90`
- `src/crop/tillage.f90`
- `src/crop/oxygenstress.f90`

**Subroutines to Refactor**:
- `cropgrowth()` - Main crop growth routine
- `wofostnut()` - WOFOST nutrient routines
- `rootextraction()` - Root water extraction
- `irrigation()` - Irrigation management
- `tillage()` - Tillage operations
- `oxygenstress()` - Oxygen stress calculations

**Call Sites to Update**:
- All calls from main time loop and parent routines

**Testing Requirements**:
```bash
pixi run test-linux-grassgrowth
pixi run test-linux-oxygenstress
pixi run regression_all
```

**Success Criteria**:
- All crop routines accept `crop_state_t`, `soil_state_t`, `atmosphere_state_t`
- No module dependencies in crop routines
- Crop growth behavior unchanged

**Dependencies**: 
- Milestone 1 completed

**Effort Estimate**: 15 hours

### Milestone 4: Drainage and Surface Water Integration
**Objective**: Convert drainage and surface water routines to state-based interface

**Files to Modify**:
- `src/drainage/drainage.f90`
- `src/drainage/surfacewater.f90`
- `src/drainage/divdra.f90`

**Subroutines to Refactor**:
- `drainage()` - Main drainage routine
- `surfacewater()` - Surface water calculations
- `divdra()` - Drainage division calculations

**Call Sites to Update**:
- All calls from main time loop and parent routines

**Testing Requirements**:
```bash
pixi run test-linux-surfacewater
pixi run test-linux-hupselbrook
pixi run regression_all
```

**Success Criteria**:
- All drainage routines accept `drainage_state_t` and `soil_state_t`
- No module dependencies in drainage routines
- Drainage behavior unchanged

**Dependencies**: 
- Milestone 1 completed

**Effort Estimate**: 10 hours

### Milestone 5: Atmosphere and Heat Integration
**Objective**: Convert atmosphere and heat routines to state-based interface

**Files to Modify**:
- `src/atmosphere/meteoday.f90`
- `src/atmosphere/meteodt.f90`
- `src/atmosphere/penmon.f90`
- `src/heat/temperature.f90`
- `src/heat/frozencond.f90`

**Subroutines to Refactor**:
- `meteoday()` - Daily meteorology calculations
- `meteodt()` - Meteorology time step calculations
- `penmon()` - Penman-Monteith calculations
- `temperature()` - Temperature calculations
- `frozencond()` - Frozen conditions calculations

**Call Sites to Update**:
- All calls from main time loop and parent routines

**Testing Requirements**:
```bash
pixi run test-linux-hupselbrook
pixi run regression_all
```

**Success Criteria**:
- All atmosphere/heat routines accept `atmosphere_state_t`
- No module dependencies in atmosphere routines
- Meteorology behavior unchanged

**Dependencies**: 
- Milestone 1 completed

**Effort Estimate**: 10 hours

### Milestone 6: Solute and Macropore Integration
**Objective**: Convert solute transport and macropore routines to state-based interface

**Files to Modify**:
- `src/solute/solute.f90`
- `src/macropore/macropore.f90`
- `src/macropore/macrorate.f90`

**Subroutines to Refactor**:
- `solute()` - Solute transport calculations
- `macropore()` - Macropore flow calculations
- `macrorate()` - Macropore rate calculations

**Call Sites to Update**:
- All calls from main time loop and parent routines

**Testing Requirements**:
```bash
pixi run test-linux-hupselbrook
pixi run regression_all
```

**Success Criteria**:
- All solute/macropore routines accept relevant state types
- No module dependencies in solute/macropore routines
- Transport behavior unchanged

**Dependencies**: 
- Milestone 1 completed

**Effort Estimate**: 8 hours

## Integration Patterns Mapping

### Pattern 1: Simple Subroutine Conversion
**Applicable to**: `checkmassbal()`, `penmon()`, `frozencond()`
**Before**:
```fortran
subroutine checkmassbal()
  use Variables, only: theta, h, gwl
  implicit none
  ! Calculation using module globals
end subroutine
```
**After**:
```fortran
subroutine checkmassbal(soil, dt)
  use iso_fortran_env, only: real64
  implicit none
  type(soil_state_t), intent(in) :: soil
  real(real64), intent(in) :: dt
  ! Calculation using soil%theta, soil%h, soil%gwl
end subroutine
```

### Pattern 2: Cross-Module Dependencies
**Applicable to**: `cropgrowth()`, `soilwater()`
**Before**:
```fortran
subroutine cropgrowth()
  use Variables, only: theta, gwl, rainfall, temp, lai, root_depth
  implicit none
  ! Uses multiple module variables
end subroutine
```
**After**:
```fortran
subroutine cropgrowth(crop, soil, atm, dt)
  use iso_fortran_env, only: real64
  implicit none
  type(crop_state_t), intent(inout) :: crop
  type(soil_state_t), intent(in) :: soil
  type(atmosphere_state_t), intent(in) :: atm
  real(real64), intent(in) :: dt
  ! Uses crop%lai, soil%theta, atm%temp, etc.
end subroutine
```

### Pattern 3: Nested Call Chain
**Applicable to**: `soilwater()` and related routines
**Before**:
```fortran
subroutine soilwater()
  use Variables, only: theta, h, gwl
  implicit none
  call fluxes()  ! Nested call
  call watstor() ! Nested call
end subroutine
```
**After**:
```fortran
subroutine soilwater(soil, boundary, dt)
  use iso_fortran_env, only: real64
  implicit none
  type(soil_state_t), intent(inout) :: soil
  type(boundary_state_t), intent(in) :: boundary
  real(real64), intent(in) :: dt
  call fluxes(soil, boundary, dt)  ! Pass state down
  call watstor(soil, dt)           ! Pass state down
end subroutine
```

### Pattern 4: Temporary Bridge (Phase 2 Only)
**Applicable to**: Routines not yet converted but called by converted routines
**Before**:
```fortran
subroutine parent_converted(state)
  type(swap_state_t), intent(inout) :: state
  call state_to_variables(state)   ! Temporary bridge
  call legacy_child_not_yet_converted()  ! Still uses Variables
  call state_from_variables(state)   ! Temporary bridge
end subroutine
```
**After**:
```fortran
subroutine parent_converted(state)
  type(swap_state_t), intent(inout) :: state
  call modern_child(state%soil, state%atm)  ! All converted
end subroutine
```

## File-by-File Refactoring Sequence

### High Priority Files (Bottom-Up)
1. `src/core/swap_state_mod.f90` - State definitions (already done)
2. `src/core/swap_state_sync.f90` - Sync procedures (already done)
3. `src/core/swap.f90` - Main time loop (Milestone 1)
4. `src/soil/checkmassbal.f90` - Simple subroutine (Milestone 2)
5. `src/soil/fluxes.f90` - Simple subroutine (Milestone 2)
6. `src/soil/headcalc.f90` - Simple subroutine (Milestone 2)
7. `src/crop/oxygenstress.f90` - Simple subroutine (Milestone 3)
8. `src/drainage/divdra.f90` - Simple subroutine (Milestone 4)
9. `src/atmosphere/penmon.f90` - Simple subroutine (Milestone 5)
10. `src/heat/frozencond.f90` - Simple subroutine (Milestone 5)

### Medium Priority Files
11. `src/soil/soilwater.f90` - Complex subroutine (Milestone 2)
12. `src/crop/cropgrowth.f90` - Complex subroutine (Milestone 3)
13. `src/drainage/drainage.f90` - Complex subroutine (Milestone 4)
14. `src/solute/solute.f90` - Complex subroutine (Milestone 6)
15. `src/macropore/macropore.f90` - Complex subroutine (Milestone 6)

### Low Priority Files
16. `src/soil/watstor.f90` - Medium complexity (Milestone 2)
17. `src/crop/wofostnut.f90` - Medium complexity (Milestone 3)
18. `src/atmosphere/meteoday.f90` - Medium complexity (Milestone 5)
19. `src/heat/temperature.f90` - Medium complexity (Milestone 5)
20. `src/macropore/macrorate.f90` - Medium complexity (Milestone 6)

## Testing Matrix

### Unit Tests
- **State Type Tests**: `test_*_state.f90` files (already done)
- **Individual Routine Tests**: Add tests for converted routines
- **Multi-instance Tests**: `test_multi_instance.f90` (already done)

### Integration Tests
- **Hupselbrook**: Verify water balance preservation
- **Grassgrowth**: Verify crop growth behavior
- **Surfacewater**: Verify drainage behavior
- **Oxygenstress**: Verify oxygen stress calculations

### Regression Validation
- **All Regression Tests**: Ensure identical numerical results
- **Water Balance Check**: Verify < 0.1% error
- **Performance Benchmarking**: Ensure < 5% performance degradation

### Performance Impact Assessment
- Run baseline performance test
- Compare with converted code
- Monitor memory usage and execution time

## Documentation Updates

### REFACTORING_PLAN.md
- Update to reflect current state of integration
- Document completed milestones
- Add new integration patterns used
- Update progress tracking

### Code Comments
- Add clear documentation of state-based interface
- Update function signatures with intent declarations
- Add comments about temporary bridges (if needed)

### Architecture Diagrams
- Update flow diagrams to show state-based data flow
- Update integration architecture diagram to show final state

## Risk Assessment

### High-Risk Changes
1. **Soil Water Routines** (`soilwater.f90`, `watstor.f90`)
   - **Risk**: Complex physics calculations
   - **Mitigation**: Verify against regression tests, test water balance

2. **Crop Growth Routines** (`cropgrowth.f90`, `wofostnut.f90`)
   - **Risk**: Complex interaction with multiple subsystems
   - **Mitigation**: Test with grassgrowth and oxygenstress test cases

3. **Drainage Routines** (`drainage.f90`, `surfacewater.f90`)
   - **Risk**: Interaction with groundwater calculations
   - **Mitigation**: Test with surfacewater test case

### Mitigation Strategies
1. **Incremental Approach**: Convert one subsystem at a time
2. **Comprehensive Testing**: Run all regression tests after each milestone
3. **Code Reviews**: Peer review of all state-based conversions
4. **Performance Monitoring**: Track performance impact throughout process

## Success Criteria

✅ **All computational subroutines accept state arguments**
✅ **Main time loop operates purely on state**
✅ **No sync calls inside main time loop**
✅ **Variables module only used for input/output**
✅ **Water balance error < 0.1%**
✅ **All regression tests pass**
✅ **Multi-instance independence verified**
✅ **Performance within 5% of baseline**
✅ **Ready for Phase 3 (BMI integration)**

This plan provides a structured approach to converting SWAP to a state-based architecture while maintaining numerical accuracy and backward compatibility. Each milestone builds upon the previous one, ensuring safe and verifiable progress.