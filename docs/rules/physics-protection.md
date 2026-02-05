# Physics Protection Rules

## CRITICAL: Never Modify Physics

These rules have HIGHEST PRIORITY and override all other considerations:

❌ **NEVER** modify equations or algorithmic logic
❌ **NEVER** change variable units or semantics
❌ **NEVER** alter numerical methods or time-stepping schemes
❌ **NEVER** change water balance calculations
❌ **NEVER** modify convergence criteria or tolerances

✅ **ONLY** modify call signatures and state passing
✅ **ONLY** refactor data structures (not algorithms)
✅ **ONLY** change how data flows between routines

## When Uncertain

**ALWAYS ASK** before making changes to:
- Fortran files in `src/soil/`, `src/crop/`, `src/drainage/`
- Any file containing mathematical equations
- Water balance, mass balance, or energy balance calculations
- Convergence loops or iterative solvers

## Examples

### ✅ SAFE: Interface change only
```fortran
! Before
subroutine soilwater()
  use Variables, only: theta, h
  theta = theta + delta_theta  ! Physics unchanged
end

! After  
subroutine soilwater(soil)
  type(soil_state_t) :: soil
  soil%theta = soil%theta + delta_theta  ! Same physics, different container
end

❌ UNSAFE: Physics modification

! Before
theta = theta + delta_theta

! After - DON'T DO THIS
theta = theta + 0.95 * delta_theta  ! Changed equation!