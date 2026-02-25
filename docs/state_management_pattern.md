---
title: State Management
author: Mateusz Zawadzki
---

In the recent refactoring effort, the monolythic, Fortran 77 style variables module has been replaced by more modular architectural pattern. In each distinct section of the model (e.g., atmosphere, crop) a module was created.

```fortran
!> {Module_Name} Module
!!
!! Brief description of what this module manages (e.g., "Manages soil water
!! movement, pressure head distribution, and water balance tracking").
!!
!! Key components:
!! - Configuration: Immutable parameters read from file, including initial
!!   condition seeds that populate state at t=0. Read from a TOML config file (replacement of .swp, etc)
!! - State: Mutable variables that evolve during simulation.
!!
!! Lifecycle summary:
!!   1. {module_name}_config_init(config, nlayers)  - allocate and zero config
!!   2. <populate config>                           - fill from file/defaults
!!   3. {module_name}_state_init(state, config)     - allocate state, seed from config
!!   4. <simulation loop>
!!   5. {module_name}_state_finalize(state)         - deallocate state
!!   6. {module_name}_config_finalize(config)       - deallocate config
!!
!! Reset contract:
!!   reset_all is always PUBLIC. The main controller and BMI both require the
!!   ability to zero state mid-simulation without re-running init.
!!   reset_all always means: set all variables to 0 / .false.
!!   It never restores initial conditions. To restore ICs, call state_init again.
!!
!! Invariants:
!!   POST config_init:    all allocatables allocated; all values zeroed/defaulted
!!   POST state_init:     all allocatables allocated; accumulators = 0; state = IC values
!!   POST reset_all:      all values 0/.false.; allocations unchanged
!!   POST state_finalize: no allocatables remain allocated; all scalars defaulted
!!
!! @author SWAP Development Team
!! @date YYYY-MM-DD
module {module_name}_mod
    implicit none
    private


    ! Public exports
    public :: {module_name}_config_t
    public :: {module_name}_state_t
    public :: {module_name}_config_init
    public :: {module_name}_config_finalize
    public :: {module_name}_state_init
    public :: {module_name}_state_finalize
    public :: {module_name}_state_reset_all          ! For BMI/warm-restart use
    public :: {module_name}_state_reset_cumulative   ! Only if module has cumulative fluxes
    public :: {module_name}_state_reset_intermediate ! Only if module has intermediate fluxes


    ! ===========================================================================
    ! Configuration Type (Immutable - defines HOW the model operates)
    ! ===========================================================================
    !> Configuration parameters for {module_name}.
    !!
    !! These parameters define the model structure, methods, physical constants,
    !! and initial condition seeds. They should not change during a simulation
    !! run. Use intent(in) when passing to procedures to enforce immutability.
    !!
    !! Initial condition seeds are fields whose only role is to populate state
    !! at t=0. They are marked with [IC seed] in their inline comment. Their
    !! names match the original SWAP variable names exactly. They are copied
    !! into state once by {module_name}_state_init and never read again during
    !! the simulation.
    !!
    !! Lifecycle: call {module_name}_config_init to allocate,
    !!            call {module_name}_config_finalize to deallocate.
    type :: {module_name}_config_t
        ! Method selection switches
        integer :: swmethod = 0               !! Method choice (e.g., 0=analytical, 1=numerical)
        integer :: swoption = 0               !! Optional feature switch

        ! Structural parameters (drive array allocation)
        integer :: nlayers = 0                !! Number of layers/compartments

        ! Physical parameters (scalars)
        real(8) :: physical_param = 0.0d0    !! Description with units (e.g., cm, m/s)

        ! Per-layer arrays (allocate in config_init, deallocate in config_finalize)
        ! Annotate each with [size: <dimension>] so allocation intent is explicit.
        real(8), allocatable :: per_level_array(:)  !! Description (units) [size: nlayers]

        ! File paths (if needed for this module)
        character(len=200) :: input_file = ''  !! Input data file path

        ! Initial condition seeds (copied into state once at state_init, then unused)
        real(8) :: var_ic = 0.0d0                    !! Description (units) [IC seed]
        real(8), allocatable :: array_ic(:)          !! Description (units) [IC seed] [size: nlayers]

    end type {module_name}_config_t


    ! ===========================================================================
    ! State Type (Mutable - evolves during simulation)
    ! ===========================================================================
    !> Runtime state for {module_name}.
    !!
    !! Contains variables that change during the simulation. Organized by
    !! physical compartments or functional groups where appropriate.
    !!
    !! Reset granularity:
    !!   reset_cumulative   - long-term mass balance accumulators (since sim start)
    !!   reset_intermediate - output-interval accumulators (between output times)
    !!   reset_all          - current state variables + both of the above
    !!
    !! Lifecycle: call {module_name}_state_init to allocate,
    !!            call {module_name}_state_finalize to deallocate.
    type :: {module_name}_state_t
        ! Current state variables
        real(8) :: state_var = 0.0d0          !! Current value with units (e.g., cm, cm/d)
        real(8), allocatable :: array_var(:)  !! Current distribution/profile [size: nlayers]

        ! Cumulative fluxes (for mass balance tracking since simulation start)
        real(8) :: cumulative_flux = 0.0d0    !! Cumulative since simulation start (cm)

        ! Intermediate fluxes (for output periods between output times)
        real(8) :: intermediate_flux = 0.0d0  !! Accumulated between outputs (cm)

        ! Flags for state-dependent behaviour
        logical :: flag_example = .false.     !! Description

    end type {module_name}_state_t


contains


    ! ===========================================================================
    ! Configuration Lifecycle
    ! ===========================================================================
    !> Allocate and initialize the {module_name} configuration.
    !!
    !! Allocates all per-layer arrays based on nlayers and zero-initializes them.
    !! Scalar members use their default initialization values from the type definition.
    !! Caller is responsible for populating members after this call.
    !!
    !! @param[out] config   Configuration container (allocated and zeroed)
    !! @param[in]  nlayers  Number of layers/compartments
    subroutine {module_name}_config_init(config, nlayers)
        type({module_name}_config_t), intent(out) :: config
        integer, intent(in) :: nlayers

        config%nlayers = nlayers

        if (nlayers > 0) then
            allocate(config%per_level_array(nlayers))
            config%per_level_array = 0.0d0
            allocate(config%array_ic(nlayers))
            config%array_ic = 0.0d0
        end if

    end subroutine {module_name}_config_init


    !> Finalize the {module_name} configuration.
    !!
    !! Deallocates all allocatable members and resets scalars to defaults.
    !! Called when the model instance is destroyed.
    !!
    !! @param[inout] config  Configuration container
    subroutine {module_name}_config_finalize(config)
        type({module_name}_config_t), intent(inout) :: config

        if (allocated(config%per_level_array)) deallocate(config%per_level_array)
        if (allocated(config%array_ic))        deallocate(config%array_ic)

        config%swmethod       = 0
        config%swoption       = 0
        config%nlayers        = 0
        config%physical_param = 0.0d0
        config%var_ic         = 0.0d0
        config%input_file     = ''

    end subroutine {module_name}_config_finalize


    ! ===========================================================================
    ! State Lifecycle
    ! ===========================================================================
    !> Initialize the {module_name} state from configuration.
    !!
    !! Allocates arrays based on config dimensions, seeds current state from
    !! config IC seed fields, and zeros all accumulators.
    !!
    !! NOTE: Add extra integer parameters when state arrays depend on dimensions
    !! not captured in config (e.g., numnod for spatial discretization):
    !!   subroutine {module_name}_state_init(state, config, numnod)
    !!
    !! @param[out] state   {Module_name} state container (initialized)
    !! @param[in]  config  Configuration parameters (immutable)
    subroutine {module_name}_state_init(state, config)
        type({module_name}_state_t),  intent(out) :: state
        type({module_name}_config_t), intent(in)  :: config

        ! Allocate arrays based on configuration
        if (config%nlayers > 0) then
            allocate(state%array_var(config%nlayers))
        end if

        ! Seed current state from config IC fields
        state%state_var = config%var_ic
        if (allocated(config%array_ic) .and. allocated(state%array_var)) then
            state%array_var = config%array_ic
        else if (allocated(state%array_var)) then
            state%array_var = 0.0d0
        end if

        ! Zero all accumulators
        call {module_name}_state_reset_cumulative(state)
        call {module_name}_state_reset_intermediate(state)

        ! Initialize flags
        state%flag_example = .false.

    end subroutine {module_name}_state_init


    !> Finalize the {module_name} state.
    !!
    !! Deallocates arrays and resets scalars to prevent state carryover
    !! between model instances. Important for BMI compliance and testing.
    !!
    !! Finalize rule (normative):
    !!   state_finalize MUST deallocate all allocatables and reset scalars inline.
    !!   It MUST NOT call reset_all. After deallocation, all allocated() guards in
    !!   reset_all evaluate to .false., so array assignments become silent no-ops.
    !!   Scalar resets must be written explicitly so finalize is self-contained
    !!   and auditable.
    !!
    !! @param[inout] state  {Module_name} state container
    subroutine {module_name}_state_finalize(state)
        type({module_name}_state_t), intent(inout) :: state

        if (allocated(state%array_var)) deallocate(state%array_var)

        state%state_var         = 0.0d0
        state%cumulative_flux   = 0.0d0
        state%intermediate_flux = 0.0d0
        state%flag_example      = .false.

    end subroutine {module_name}_state_finalize


    ! ===========================================================================
    ! Reset Procedures
    ! ===========================================================================
    !> Reset all {module_name} state fields to zero/false.
    !!
    !! Resets current state variables, cumulative, and intermediate accumulators.
    !! Does not allocate/deallocate memory.
    !! Does not restore initial conditions — call state_init for that.
    !! Use for warm-restart and BMI reset operations.
    !!
    !! @param[inout] state  {Module_name} state container
    subroutine {module_name}_state_reset_all(state)
        type({module_name}_state_t), intent(inout) :: state

        state%state_var = 0.0d0
        if (allocated(state%array_var)) state%array_var = 0.0d0

        call {module_name}_state_reset_cumulative(state)
        call {module_name}_state_reset_intermediate(state)

        state%flag_example = .false.

    end subroutine {module_name}_state_reset_all


    !> Reset cumulative flux accumulators.
    !!
    !! Called when flzerocumu=.true. or at initialization.
    !! Resets long-term cumulative fluxes tracked since simulation start.
    !!
    !! @param[inout] state  {Module_name} state container
    subroutine {module_name}_state_reset_cumulative(state)
        type({module_name}_state_t), intent(inout) :: state

        state%cumulative_flux = 0.0d0

    end subroutine {module_name}_state_reset_cumulative


    !> Reset intermediate flux accumulators.
    !!
    !! Called when flzerointr=.true. or at initialization.
    !! Resets output-interval accumulators tracked between output times.
    !!
    !! @param[inout] state  {Module_name} state container
    subroutine {module_name}_state_reset_intermediate(state)
        type({module_name}_state_t), intent(inout) :: state

        state%intermediate_flux = 0.0d0

    end subroutine {module_name}_state_reset_intermediate


end module {module_name}_mod 
```
