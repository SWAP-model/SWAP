# SWAP Coding Standards

## Fortran Standards

### File Naming
- Module files: `swap_modulename_mod.f90`
- Subroutine-only files: `swap_functionality.f90`
- Test files: `test_modulename.f90`

### Naming Conventions
- State types: `modulename_state_t` (e.g., `soil_state_t`)
- when adding the states in the symbol signatures, use shortened names.
- Module names: `swap_modulename` (e.g., `swap_state_mod`)
- Public procedures: `modulename_action` (e.g., `state_from_variables`)
- Local variables: lowercase with underscores

### Code Structure
```fortran
module swap_example_mod
  use iso_fortran_env, only: real64
  implicit none
  private
  
  ! Public API
  public :: example_state_t
  public :: example_init
  public :: example_compute
  
  ! Type definitions
  type :: example_state_t
    real(real64) :: variable1
    real(real64) :: variable2
  end type
  
contains
  
  subroutine example_init(state)
    type(example_state_t), intent(out) :: state
    ! Initialize state
  end subroutine
  
end module

Intent Declarations (MANDATORY)

    intent(in) - Read-only input

    intent(out) - Output (uninitialized on entry)

    intent(inout) - Modified input/output

Always specify intent for all arguments
Documentation

fortran
!> Brief description of subroutine
!!
!! Detailed description if needed
!!
!! @param[inout] state Soil state containing theta, h, gwl
!! @param[in]    dt    Time step in days
subroutine soilwater(state, dt)

Python Standards
Style

    Follow PEP 8

    Type hints required for all function signatures

    Docstrings in Google style

Example

python
def initialize_state(config: dict[str, Any]) -> SwapState:
    """Initialize SWAP state from configuration.
    
    Args:
        config: Configuration dictionary with model parameters
        
    Returns:
        Initialized SwapState object
        
    Raises:
        ValueError: If required config keys are missing
    """
    # Implementation

Git Commit Standards

Use Conventional Commits format:

text
<type>(<scope>): <description>

[optional body]

[optional footer]

Types

    feat: New feature

    fix: Bug fix

    refactor: Code restructuring (no behavior change)

    test: Adding/updating tests

    docs: Documentation changes

    build: Build system changes

    chore: Maintenance tasks

Examples

text
refactor(soil): convert soilwater to state-based interface

- Update soilwater() signature to accept soil_state_t
- Update 12 call sites in main time loop
- Remove Variables module dependencies
- Tests pass: water balance error 0.02%

Refs: #42

text
test(state): add unit tests for soil_state_t initialization

- Test state allocation
- Test sync from Variables
- Test multi-instance independence

File Organization
Generally do not modify except when removing obsolete sync code after the state is fully implemented.

    src/core/swap_state_mod.f90 - State type definitions ✅

    src/core/swap_state_sync.f90 - Sync procedures ✅

    tests/unit/*/test_*_state.f90 - State unit tests ✅

Safe to Modify (Integration work)

    Process subroutines: src/soil/*.f90, src/crop/*.f90, etc.

    Main driver: src/swap.f90, src/core/swap_driver_mod.f90

    Call sites in time loop

Gradual Deprecation

    src/core/variables.f90 - Reduce dependencies over time (never delete)
