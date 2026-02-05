# SWAP State Integration Planning - Phase 2 Strategy

## Current Architecture Status

**COMPLETED (Phase 1):**
- ✅ All module-level SAVE statements removed (Steps 1-10)
- ✅ State types defined: swap_state_t, soil_state_t, crop_state_t, etc.
- ✅ Sync layer implemented: state_to_variables() and state_from_variables()
- ✅ Variables module still exists and functional (legacy compatibility)

**CURRENT STATE:**
The SWAP model currently operates in two parallel worlds:
1. **Variables module** (src/core/variables.f90) - holds all simulation data via module globals
2. **State types** (src/core/swap_state_mod.f90) - modern state containers
3. **Sync layer** (src/core/swap_state_sync.f90) - bridges Variables ↔ State

Initial input reading MUST use Variables (ttutil/legacy input system), then sync to state.

## End Goal Architecture

┌─────────────────────────────────────────────────┐
│ Python BMI Interface (Future) │
│ → Initializes state directly │
│ → No Variables dependency │
└─────────────────────────────────────────────────┘
↓
┌─────────────────────────────────────────────────┐
│ SWAP Core (State-based) │
│ → All subroutines accept state arguments │
│ → No module globals in computation │
│ → Pure functional flow │
└─────────────────────────────────────────────────┘
↓
┌─────────────────────────────────────────────────┐
│ Legacy Wrapper (Optional) │
│ → Variables + ttutil preserved │
│ → Backwards compatibility for old inputs │
│ → Sync layer only used here │
└─────────────────────────────────────────────────┘

text

## Phase 2 Strategy: Progressive Sync Elimination

**Step-by-step approach:**

1. **Initial state** (TODAY):
   ```fortran
   ! Main time loop uses Variables everywhere
   call state_from_variables(state)  ! Sync once at start
   ! Time loop uses Variables (module globals)
   call state_to_variables(state)    ! Sync once at end

    Intermediate state (PHASE 2 TARGET):

    fortran
    ! Main time loop wired with state
    call state_from_variables(state)  ! Initialize from legacy input

    ! Time loop uses STATE (no more sync inside loop!)
    do time_step = 1, nsteps
      call soilwater(state%soil, state%boundary)      ! State-based
      call cropgrowth(state%crop, state%soil, state%atm)  
      call drainage(state%drain, state%soil)
      ! No sync needed - everything flows through state
    end do

    call state_to_variables(state)    ! Final sync for legacy output

    End state (PHASE 3 - BMI):

    fortran
    ! Python BMI initializes state directly (no Variables needed)
    call bmi_initialize(config_dict, state)  ! From Python dict/JSON

    ! Time loop purely state-based
    do time_step = 1, nsteps
      call soilwater(state%soil, state%boundary)
      call cropgrowth(state%crop, state%soil, state%atm)
      call drainage(state%drain, state%soil)
    end do

    call bmi_finalize(state)  ! Return state to Python
    ! No Variables module touched - completely independent

Planning Task

Analyze the SWAP codebase and create a detailed integration plan that:
1. Call Chain Analysis

Identify the complete dependency graph:

    Map all subroutines in the main time loop

    For each subroutine, identify:

        What Variables module globals it reads

        What Variables module globals it writes

        What state types it needs (soil_state_t, crop_state_t, boundary_state_t, etc.)

        All downstream calls it makes

2. Integration Priority Order

Determine the optimal sequence for wiring state through subroutines:

    Start with leaf subroutines (no downstream dependencies)

    Progress to parent routines once children are state-based

    Group by subsystem (soil water, crop, drainage, etc.)

    Consider testing/validation checkpoints between groups

3. Sync Elimination Strategy

For each integration milestone, specify:

    What sync calls can be removed (moved outside time loop)

    What sync calls must remain temporarily (why?)

    When we can eliminate sync entirely from computational kernels

4. Backwards Compatibility Plan

Ensure Variables module remains functional:

    What functions/subroutines must remain Variables-compatible?

    Where do we maintain the legacy wrapper boundary?

    How do we test both pathways (state-based + legacy)?

5. Testing & Validation

For each integration step, define:

    Regression test requirements (pixi run test-linux-hupselbrook)

    Water balance validation (must be identical to baseline)

    Performance impact assessment

    Multi-instance independence verification

Constraints

MUST PRESERVE:

    ✅ Numerical accuracy (water balance error <0.1%)

    ✅ Physics algorithms (no equation modifications)

    ✅ Legacy input file compatibility (Variables/ttutil path still works)

    ✅ Backward compatibility (existing SWAP workflows unaffected)

MUST ELIMINATE:

    ❌ sync calls inside main time loop

    ❌ module global dependencies in computational kernels

    ❌ Hidden state sharing between process subroutines

CAN DEFER TO PHASE 3:

    Python BMI interface implementation

    Direct state initialization from config files

    Complete Variables module removal

Deliverables

Provide a comprehensive plan including:

    File-by-file refactoring sequence (which files in what order?)

    Detailed call signature changes (before/after for key subroutines)

    Integration milestones (what constitutes a "done" checkpoint?)

    Test validation strategy (how to verify each step preserves physics?)

    Risk assessment (what are the highest-risk changes? How to mitigate?)

    Time/effort estimates (which subsystems are most complex?)

Success Criteria

The plan is complete when:

    All computational subroutines accept state arguments (no module globals)

    Main time loop operates purely on state (no sync inside loop)

    Variables module only used for: input reading, legacy wrapper, optional output

    All regression tests pass with identical numerical results

    Code is ready for Phase 3: BMI integration from Python

Focus on architectural clarity and incremental safety. Every step must be verifiable through regression testing before proceeding.