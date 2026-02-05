# Integration Patterns - Phase 2

## Pattern 1: Simple Subroutine Conversion

**Use when:** Subroutine has minimal dependencies, no nested calls

**Before:**
```fortran
subroutine compute_drainage()
  use Variables, only: gwl, drain_flux, drain_level
  implicit none
  
  if (gwl > drain_level) then
    drain_flux = calc_flux(gwl - drain_level)
  end if
end subroutine

After:

fortran
subroutine compute_drainage(drain, soil)
  use iso_fortran_env, only: real64
  implicit none
  type(drainage_state_t), intent(inout) :: drain
  type(soil_state_t), intent(in) :: soil
  
  if (soil%gwl > drain%level) then
    drain%flux = calc_flux(soil%gwl - drain%level)
  end if
end subroutine

Steps:

    Identify Variables dependencies → map to state types

    Add state arguments to signature

    Replace Variables::var with state%var

    Update call sites

    Remove Variables USE statement

Pattern 2: Cross-Module Dependencies

Use when: Subroutine reads from multiple subsystems

Before:

fortran
subroutine cropgrowth()
  use Variables, only: theta, gwl          ! From soil
  use Variables, only: rainfall, temp      ! From atmosphere
  use Variables, only: lai, root_depth     ! From crop
  
  growth_rate = calc_growth(theta, temp, lai)
  lai = lai + growth_rate * dt
end subroutine

After:

fortran
subroutine cropgrowth(crop, soil, atm, dt)
  use iso_fortran_env, only: real64
  implicit none
  type(crop_state_t), intent(inout) :: crop
  type(soil_state_t), intent(in) :: soil
  type(atmosphere_state_t), intent(in) :: atm
  real(real64), intent(in) :: dt
  
  growth_rate = calc_growth(soil%theta, atm%temp, crop%lai)
  crop%lai = crop%lai + growth_rate * dt
end subroutine

Pattern 3: Nested Call Chain

Use when: Parent calls children that need state

Top-level:

fortran
subroutine time_step_driver(state, dt)
  type(swap_state_t), intent(inout) :: state
  real(real64), intent(in) :: dt
  
  ! Pass relevant state portions down
  call update_atmosphere(state%atm, state%time, dt)
  call update_soil(state%soil, state%boundary, state%atm, dt)
  call update_crop(state%crop, state%soil, state%atm, dt)
  call update_drainage(state%drain, state%soil, dt)
end subroutine

Child level:

fortran
subroutine update_soil(soil, boundary, atm, dt)
  type(soil_state_t), intent(inout) :: soil
  type(boundary_state_t), intent(inout) :: boundary
  type(atmosphere_state_t), intent(in) :: atm
  real(real64), intent(in) :: dt
  
  ! Further nested calls
  call soilwater_flow(soil, boundary, dt)
  call evapotranspiration(soil, atm, dt)
end subroutine

Pattern 4: Temporary Bridge (Phase 2 Only)

Use when: Child routines not yet converted, but parent is

fortran
subroutine parent_converted(state)
  type(swap_state_t), intent(inout) :: state
  
  ! Sync state → Variables for legacy child
  call state_to_variables(state)
  
  ! Call legacy routine that still uses module globals
  call legacy_child_not_yet_converted()
  
  ! Sync Variables → state to capture changes
  call state_from_variables(state)
  
  ! Continue with state-based operations
  call modern_child(state%soil, state%atm)
end subroutine

Goal: Eliminate this pattern by converting legacy_child_not_yet_converted()
Pattern 5: State Decomposition

Use when: Large state object, but routine needs only small portion

Instead of passing entire state:

fortran
! ❌ Too much coupling
subroutine small_function(state)
  type(swap_state_t), intent(inout) :: state
  ! Only uses state%soil%theta
end subroutine

Pass only what's needed:

fortran
! ✅ Minimal coupling
subroutine small_function(soil)
  type(soil_state_t), intent(inout) :: soil
  ! Uses soil%theta
end subroutine

Integration Decision Tree

text
Is subroutine in main time loop?
├─ YES → High priority, convert now
└─ NO → Is it called by time loop routines?
    ├─ YES → Medium priority, convert after parents
    └─ NO → Low priority, defer to Phase 3

Does it modify state?
├─ YES → intent(inout)
└─ NO → intent(in)

Does it use Variables module?
├─ YES → Map to state types, convert
└─ NO → Already clean, verify call sites only

Does it have nested calls?
├─ YES → Convert children first (bottom-up)
└─ NO → Safe to convert immediately

Validation Checklist

After applying any pattern:

    All call sites updated with state arguments

    Variables USE removed (where replaced by state)

    Intent declarations specified for all arguments

    Build succeeds: pixi run build-linux

    Unit tests pass: pixi run test-unit-all

    Integration test passes: pixi run test-linux-hupselbrook

    Water balance unchanged (< 0.1% error)

    No new compiler warnings

    REFACTORING_PLAN.md updated
