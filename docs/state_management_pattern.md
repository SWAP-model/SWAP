---
title: State Management
author: Mateusz Zawadzki
---

```fortran
!> {Module_Name} State Module
!!
!! Brief description of what this module manages (e.g., "Manages soil water
!! movement, pressure head distribution, and water balance tracking").
!!
!! Key components:
!! - Configuration parameters
!! - Current state variables
!! - Cumulative and intermediate flux tracking (if applicable)
!!
!! @author SWAP Development Team
!! @date YYYY-MM-DD
module {module_name}_state_mod
    implicit none
    private

    ! Public exports
    public :: {module_name}_state_t
    public :: {module_name}_state_init
    public :: {module_name}_state_finalize
    public :: {module_name}_state_reset_cumulative    ! Only if module has cumulative fluxes
    public :: {module_name}_state_reset_intermediate  ! Only if module has intermediate fluxes

    ! ===========================================================================
    ! {Module_Name} State Type
    ! ===========================================================================
    type :: {module_name}_state_t
        ! Configuration (switches, parameters)
        integer :: switch_example = 0          !! Description of switch
        
        ! State variables
        real(8) :: state_var = 0.0d0          !! Description with units (e.g., cm, cm/d)
        
        ! Allocatable arrays (if needed)
        real(8), allocatable :: array_var(:)   !! Description with units
        
        ! Cumulative fluxes (if applicable)
        real(8) :: cumulative_flux = 0.0d0    !! Cumulative since simulation start
        
        ! Intermediate fluxes (if applicable)
        real(8) :: intermediate_flux = 0.0d0  !! Accumulated between outputs
        
        ! Flags
        logical :: flag_example = .false.     !! Description
    end type {module_name}_state_t

contains

    !> Initialize the {module_name} state.
    !!
    !! Allocates arrays (if needed) and sets all values to canonical defaults.
    !! Custom initialization parameters can be passed as optional arguments.
    !!
    !! @param[inout] state {Module_name} state container
    !! @param[in] size_param Optional: array size parameter
    subroutine {module_name}_state_init(state, size_param)
        type({module_name}_state_t), intent(inout) :: state
        integer, intent(in), optional :: size_param
        
        ! Allocate arrays if needed
        ! if (present(size_param)) then
        !     allocate(state%array_var(size_param))
        ! end if
        
        ! Reset to defaults
        call {module_name}_state_reset_all(state)
        
        ! Custom initialization (e.g., set non-zero defaults)
        ! state%some_param = sensible_default_value
        
    end subroutine {module_name}_state_init

    !> Finalize the {module_name} state.
    !!
    !! Deallocates arrays and resets values to prevent state carryover
    !! between model instances.
    !!
    !! @param[inout] state {Module_name} state container
    subroutine {module_name}_state_finalize(state)
        type({module_name}_state_t), intent(inout) :: state
        
        ! Deallocate arrays
        ! if (allocated(state%array_var)) deallocate(state%array_var)
        
        ! Reset scalars for safety
        call {module_name}_state_reset_all(state)
        
    end subroutine {module_name}_state_finalize

    !> Reset all {module_name} state fields to canonical defaults.
    !!
    !! Used by init and finalize to ensure clean state. Zeros values
    !! but does not allocate/deallocate memory.
    !!
    !! @param[inout] state {Module_name} state container
    subroutine {module_name}_state_reset_all(state)
        type({module_name}_state_t), intent(inout) :: state
        
        ! Configuration
        state%switch_example = 0
        
        ! State variables
        state%state_var = 0.0d0
        
        ! Arrays (only if allocated)
        if (allocated(state%array_var)) state%array_var = 0.0d0
        
        ! Reset fluxes via specialized procedures
        call {module_name}_state_reset_cumulative(state)
        call {module_name}_state_reset_intermediate(state)
        
        ! Flags
        state%flag_example = .false.
        
    end subroutine {module_name}_state_reset_all

    !> Reset cumulative flux accumulators.
    !!
    !! Called when flzerocumu=.true. Resets long-term cumulative fluxes
    !! tracked since simulation start or last reset.
    !!
    !! @param[inout] state {Module_name} state container
    subroutine {module_name}_state_reset_cumulative(state)
        type({module_name}_state_t), intent(inout) :: state
        
        ! Only reset cumulative fluxes (if any)
        state%cumulative_flux = 0.0d0
        
    end subroutine {module_name}_state_reset_cumulative

    !> Reset intermediate flux accumulators.
    !!
    !! Called when flzerointr=.true. Resets output-interval accumulators
    !! tracked between output times.
    !!
    !! @param[inout] state {Module_name} state container
    subroutine {module_name}_state_reset_intermediate(state)
        type({module_name}_state_t), intent(inout) :: state
        
        ! Only reset intermediate accumulators (if any)
        state%intermediate_flux = 0.0d0
        
    end subroutine {module_name}_state_reset_intermediate

end module {module_name}_state_mod
```