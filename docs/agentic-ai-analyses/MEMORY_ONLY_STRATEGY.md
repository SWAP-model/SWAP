# Strategy: Memory-Only SWAP with ttutil Replacement

**Objective**: Replace ttutil entirely with Python-based configuration, enabling SWAP to run with variables set from Python without file I/O in Fortran routines. Enable Python to orchestrate multiple processes efficiently.

**Date**: 2026-01-12  
**Branch**: swaplib-simple (baseline)  
**Status**: Planning phase

---

## Current Situation Analysis

### What We Know
1. **Baseline Working**: Basic BMI with file-based initialization works (9 functions tested)
2. **Previous Attempt Failed**: Memory setters segfaulted due to complex initialization order
3. **Core Blocker**: SAVE variables create global state (one instance per process)
4. **ttutil Complexity**: Parses 160+ parameters, 12+ tables, lazy-loads crop files
5. **Python Parser**: Already created, successfully parses all SWAP configuration

### Key Technical Constraints
- **Initialization Dependencies**: SWAP modules have circular dependencies requiring specific init order
- **SAVE Variables**: Module-level state prevents multiple instances per process
- **Array Allocations**: Many arrays allocated dynamically during initialization
- **Fortran 77 Legacy**: Much code uses old-style COMMON blocks and SAVE statements

---

## Strategic Options

### OPTION 1: Incremental BMI Setters (Low Risk, Medium Effort)
**Approach**: Build memory setters incrementally, respecting SWAP's initialization order

**Steps**:
1. **Map Initialization Sequence**
   - Trace complete initialization flow in `swap()` routine
   - Document all module variable dependencies
   - Identify minimum required state for `initialize_memory()`

2. **Create Initialization Stages**
   ```fortran
   bmi_init_stage1()  ! Basic dimensions (grid size, time settings)
   bmi_init_stage2()  ! Allocate arrays
   bmi_init_stage3()  ! Soil parameters
   bmi_init_stage4()  ! Crop parameters
   bmi_init_stage5()  ! Initial conditions
   ```

3. **Replace ttutil Calls One-by-One**
   - Start with simplest parameters (scalars, switches)
   - Move to tables (soil profiles, meteorology)
   - End with complex structures (crop rotations)

**Advantages**:
- Minimal disruption to existing SWAP code
- Can test after each increment
- Preserves current architecture
- Lower risk of breaking existing functionality

**Disadvantages**:
- Still limited to one instance per process (SAVE variables)
- Doesn't address architectural debt
- May hit same initialization issues as before
- Requires deep understanding of initialization order

**Estimated Effort**: 3-4 weeks  
**Risk Level**: Medium  

---

### OPTION 2: State Object Refactoring (High Risk, High Reward)
**Approach**: Refactor SWAP to use explicit state objects instead of SAVE variables

**Steps**:
1. **Create State Type**
   ```fortran
   type :: swap_state_t
     ! Grid configuration
     integer :: NumGrid
     real(8), allocatable :: z(:)
     
     ! Soil properties
     real(8), allocatable :: theta(:), h(:), K(:)
     
     ! Time settings
     integer :: year_start, doy_start, year_end, doy_end
     real(8) :: current_time
     
     ! All other variables currently in SAVE...
   end type
   ```

2. **Pass State Through Call Chain**
   - Modify all subroutines to accept `type(swap_state_t)`
   - Remove SAVE statements from modules
   - Thread state object through entire codebase

3. **Multiple Instance Support**
   ```fortran
   ! Python can create multiple independent SWAP instances
   state1 = swap_create()
   state2 = swap_create()
   ```

**Advantages**:
- Enables true multi-instance capability
- Modernizes code architecture
- Thread-safe by design
- Enables Python multiprocessing with shared memory

**Disadvantages**:
- Major refactoring (touching 100+ files)
- High risk of introducing bugs
- Requires extensive testing
- May break existing SWAP users/workflows
- 3-6 months of work minimum

**Estimated Effort**: 3-6 months  
**Risk Level**: High  

---

### OPTION 3: Hybrid Python-Fortran Split (Pragmatic)
**Approach**: Keep ttutil/SWAP as-is, orchestrate at Python level

**Steps**:
1. **Python Configuration Layer**
   - Use existing `swap_parser.py` to read/modify configs
   - Generate temporary config files per SWAP instance
   - Each Python process gets isolated temp directory

2. **BMI Wrapper Enhancement**
   ```python
   class SwapInstance:
       def __init__(self, config_dict):
           self.temp_dir = tempfile.mkdtemp()
           self.write_config_files(config_dict)
           self.bmi = SwapBMI(f"{self.temp_dir}/model.swp")
   ```

3. **Parallel Processing**
   - Use Python `multiprocessing` with isolated temp dirs
   - Each process runs independent SWAP instance
   - Collect results at Python level

**Advantages**:
- Zero changes to SWAP/ttutil
- Works with current architecture
- Python provides flexibility
- Can start immediately
- File I/O overhead not critical for many workflows

**Disadvantages**:
- Doesn't achieve "pure memory" goal
- Temp file management overhead
- Not elegant for high-frequency coupling
- File I/O still present

**Estimated Effort**: 1-2 weeks  
**Risk Level**: Low  

---

### OPTION 4: Generated Fortran Interface (Smart Middle Ground)
**Approach**: Auto-generate BMI setters from YAML registry of SWAP parameters

**Steps**:
1. **Create Parameter Registry**
   ```yaml
   # swap_params.yaml
   parameters:
     - name: NumGrid
       type: integer
       module: variables
       required: true
       init_stage: 1
     
     - name: z
       type: real8_array
       module: variables
       depends_on: [NumGrid]
       init_stage: 2
   ```

2. **Code Generator**
   ```python
   # generate_bmi_setters.py
   # Read YAML -> Generate swap_bmi_generated.f90
   # Creates setters respecting dependencies
   ```

3. **Safe Initialization**
   - Generator ensures correct order based on `init_stage`
   - Validates dependencies
   - Creates Python bindings automatically

**Advantages**:
- Documents all parameters systematically
- Reduces manual coding errors
- Easy to maintain/extend
- Can generate Python bindings too
- Clear separation of concerns

**Disadvantages**:
- Upfront effort to create registry
- Still one instance per process (SAVE variables)
- Requires YAML maintenance
- Code generation adds complexity

**Estimated Effort**: 4-6 weeks  
**Risk Level**: Medium  

---

### OPTION 5: Fortran Wrapper Library (Isolation Pattern)
**Approach**: Create thin C-compatible wrapper that manages SWAP state externally

**Steps**:
1. **External State Management**
   ```fortran
   ! swap_wrapper.f90
   type swap_handle
     integer :: id
     ! Opaque pointer to SWAP state
   end type
   
   function swap_create() bind(C)
     ! Allocate new SWAP state in heap
     ! Copy initial values from modules
   end function
   ```

2. **State Serialization**
   - Save SWAP module state to heap before each call
   - Restore specific instance state before operations
   - Copy back after operations

3. **Python Interface**
   ```python
   swap1 = SwapWrapper.create()
   swap2 = SwapWrapper.create()
   # Each tracks independent state
   ```

**Advantages**:
- No changes to core SWAP code
- Multiple instances possible
- Minimal Fortran modifications

**Disadvantages**:
- Performance overhead (state copying)
- Complex state serialization logic
- Memory management complexity
- Still limited by SAVE variables underneath

**Estimated Effort**: 6-8 weeks  
**Risk Level**: Medium-High  

---

## Recommended Pathway

### PHASE 1: Quick Win (Weeks 1-2)
**Goal**: Prove Python orchestration works without SWAP changes

Implement **OPTION 3 (Hybrid)** first:
- Enhance `swap_parser.py` with config writing
- Create `SwapInstance` Python class with temp directories
- Test multiprocessing with current BMI
- **Deliverable**: Running N SWAP instances in parallel from Python

### PHASE 2: Foundation (Weeks 3-6)
**Goal**: Build systematic parameter interface

Implement **OPTION 4 (Generated Interface)**:
- Create YAML parameter registry (start with critical 20-30 params)
- Build code generator for BMI setters
- Test incremental replacement of ttutil calls
- **Deliverable**: Memory-only initialization for subset of parameters

### PHASE 3: Deep Refactoring (Months 2-4)
**Goal**: Enable true multi-instance capability

Begin **OPTION 2 (State Object)**:
- Refactor one module at a time (start with `variables.f90`)
- Create comprehensive test suite
- Gradual migration of SAVE variables
- **Deliverable**: Modernized SWAP core supporting multiple instances

---

## Implementation Priorities

### Critical Path Parameters (Must replace ttutil first)
1. **Time settings**: `year_start`, `year_end`, `doy_start`, `doy_end`
2. **Grid configuration**: `NumGrid`, `z` array
3. **Soil hydraulics**: `thetaS`, `thetaR`, `Ksat`, `alpha`, `n`, `lambda`
4. **Initial conditions**: Initial pressure heads, water contents
5. **Boundary conditions**: Meteorology, bottom boundary type

### Secondary Parameters (Can defer)
- Detailed output options
- Crop rotation schedules
- Management practices
- Advanced coupling options

---

## Technical Deep Dive: Why Previous Attempt Failed

### Root Cause Analysis
```fortran
! In swap_bmi_setters.f90 (previous attempt)
function set_time_settings(...) bind(C)
  ! This failed because:
  ! 1. Variables module not initialized (arrays not allocated)
  ! 2. swap_exchange module not connected
  ! 3. Time-dependent arrays not sized yet
  ! 4. No call to underlying SWAP initialization routines
end function
```

### Correct Initialization Sequence
```fortran
! What MUST happen (derived from SWAP source):
1. Call swap(iCaller=0, iTask=1)  ! Opens files, reads dimensions
2. Allocate all arrays based on dimensions
3. Read soil parameters
4. Read crop parameters  
5. Initialize state variables
6. Ready for time stepping

! Memory-only version needs:
1. bmi_set_dimensions()           ! Replaces file reading step 1
2. bmi_allocate_arrays()          ! Internal allocation
3. bmi_set_soil_params()          ! Replaces file reading step 3
4. bmi_set_crop_params()          ! Replaces file reading step 4
5. bmi_set_initial_conditions()   ! Replaces file reading step 5
6. bmi_finalize_initialization()  ! Connects all modules
7. Ready for time stepping
```

---

## Testing Strategy

### Unit Tests
- Each BMI setter tested independently
- Verify array allocation before data setting
- Check module variable consistency

### Integration Tests
- Full memory-only initialization sequence
- Compare results with file-based initialization
- Multi-instance tests (if state refactoring done)

### Performance Benchmarks
- Memory-only vs file-based startup time
- Multi-instance overhead measurement
- Parallel scaling tests (1, 2, 4, 8, 16 processes)

---

## Decision Matrix

| Option | Effort | Risk | Multi-Instance | Pure Memory | Timeline |
|--------|--------|------|----------------|-------------|----------|
| 1. Incremental BMI | Medium | Medium | ❌ | ✅ | 3-4 weeks |
| 2. State Refactor | Very High | High | ✅ | ✅ | 3-6 months |
| 3. Hybrid Python | Low | Low | ✅ | ❌ | 1-2 weeks |
| 4. Generated Interface | Medium-High | Medium | ❌ | ✅ | 4-6 weeks |
| 5. Wrapper Library | High | Medium-High | ⚠️ | ⚠️ | 6-8 weeks |

**Legend**: ✅ Full support | ⚠️ Partial support | ❌ Not supported

---

## Next Steps

### Immediate Actions (Today)
1. ✅ Review this strategy document
2. Choose initial pathway (recommend: Option 3 → Option 4 → Option 2)
3. Set up development branch
4. Define acceptance criteria for Phase 1

### This Week
1. Implement Phase 1 (Hybrid approach)
2. Test multiprocessing with existing BMI
3. Measure baseline performance
4. Start YAML parameter registry

### This Month
1. Complete parameter registry (top 30 parameters)
2. Build code generator
3. Test generated BMI setters
4. Compare memory vs file performance

---

## Risk Mitigation

### Technical Risks
- **Initialization failures**: Comprehensive tracing of SWAP init before coding
- **Memory corruption**: Valgrind testing at each step
- **Performance degradation**: Benchmark after each change
- **Breaking changes**: Maintain file-based interface alongside memory interface

### Project Risks
- **Scope creep**: Stick to phased approach, don't jump to Option 2 immediately
- **Testing burden**: Automate testing from day 1
- **Knowledge silos**: Document initialization dependencies thoroughly

---

## Questions to Resolve

1. **Performance Requirements**: What's acceptable initialization time? (file: ~0.2s, target: ?)
2. **Multi-instance Priority**: Is true multi-instance critical, or is multiprocessing sufficient?
3. **SWAP Code Modification**: Can we modify SWAP core, or must it remain unchanged?
4. **Python Version**: Python 3.9+? Any constraints?
5. **Deployment**: Single machine or HPC cluster? Affects multiprocessing approach.

---

## Success Criteria

### Phase 1 Success
- [ ] 10+ SWAP instances running in parallel
- [ ] Python orchestration working
- [ ] Results identical to file-based runs
- [ ] Performance acceptable for use case

### Phase 2 Success
- [ ] 50+ critical parameters settable from Python
- [ ] No ttutil calls for basic simulations
- [ ] Generated code compiles and passes tests
- [ ] Documentation complete

### Phase 3 Success
- [ ] SAVE variables eliminated from core modules
- [ ] Multiple instances in single process
- [ ] Comprehensive test coverage
- [ ] Performance at least 80% of file-based version

