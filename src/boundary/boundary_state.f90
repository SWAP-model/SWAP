!> Boundary Conditions State Module
!!
!! This module defines the top and bottom boundary condition state for SWAP.
!! Manages prescribed fluxes, pressure heads, groundwater levels, aquifer
!! interactions, surface ponding, runoff, and runon processes.
!!
!! Key components:
!! - Bottom boundary configuration (8 BC types with various parameterizations)
!! - Aquifer interaction parameters and prescribed tables
!! - Top boundary configuration (ponding, runoff, runon)
!! - Cumulative and intermediate flux tracking
!!
!! @author SWAP Development Team
!! @date 2026-02-13
module boundary_state_mod
    use swap_array_dimensions, only: MADAY, MABBC, MAIRG
   implicit none
   private

   public :: boundary_state_t
   public :: boundary_state_init, boundary_state_finalize

   type :: boundary_state_t
        ! Bottom boundary - configuration
        integer :: swbotb = 0                  !! Bottom BC type (1-8)
        integer :: swbotb3Impl = 0             !! Implicit solution switch for BC type 3
        integer :: SwBotb3ResVert = 0          !! Suppress vertical resistance for BC type 3
        integer :: swqhbot = 0                 !! Flux-GWL relationship type
        integer :: swcofqhc = 0                !! Additional flux switch
        integer :: sw2 = 0                     !! Sub-switch for BC type 2
        integer :: sw3 = 0                     !! Sub-switch for BC type 3
        integer :: sw4 = 0                     !! Sub-switch for BC type 4 (extra flux)
        
        ! Bottom boundary - values
        real(8) :: qbot = 0.0d0                !! Bottom flux (cm/d)
        real(8) :: qbot_nonfrozen = 0.0d0      !! Bottom flux for non-frozen soil (cm/d)
        real(8) :: hbot = 0.0d0                !! Bottom pressure head (cm)
        real(8) :: iqbot = 0.0d0               !! Intermediate bottom flux
        real(8) :: cqbot = 0.0d0               !! Cumulative bottom flux
        real(8) :: cqbotdo = 0.0d0             !! Cumulative downward bottom flux
        real(8) :: cqbotup = 0.0d0             !! Cumulative upward bottom flux
        real(8) :: deepgw = 0.0d0              !! Hydraulic head in aquifer (cm)
        
        ! Aquifer parameters
        real(8) :: aqave = 0.0d0               !! Average aquifer head (cm)
        real(8) :: aqamp = 0.0d0               !! Aquifer head amplitude (cm)
        real(8) :: aqper = 0.0d0               !! Period of prescribed sine wave (d)
        real(8) :: aqtmax = 0.0d0              !! Time with maximum head (d)
        real(8) :: rimlay = 0.0d0              !! Aquitard resistance (d)
        real(8) :: hdrain = 0.0d0              !! Mean drainage level (cm)
        real(8) :: shape = 0.0d0               !! Shape factor (-)
        
        ! Sine function parameters for bottom flux
        real(8) :: sinave = 0.0d0              !! Average bottom flux (cm/d)
        real(8) :: sinamp = 0.0d0              !! Amplitude of bottom flux (cm/d)
        real(8) :: sinmax = 0.0d0              !! Time of maximum flux (d)
        
        ! Flux-head relationships
        real(8) :: cofqha = 0.0d0              !! Coefficient A (exp function)
        real(8) :: cofqhb = 0.0d0              !! Coefficient B
        real(8) :: cofqhc = 0.0d0              !! Coefficient C
        
        ! Lysimeter parameters
        real(8) :: hplate = 0.0d0              !! Pressure head of ceramic plate (cm)
        
        ! Prescribed tables (allocated in init)
        real(8), allocatable :: gwltab(:)      !! GWL vs time [date, gwl]
        real(8), allocatable :: haqtab(:)      !! Aquifer head vs time [date, head]
        real(8), allocatable :: qbotab(:)      !! Bottom flux vs time [date, flux]
        real(8), allocatable :: hbotab(:)      !! Bottom h vs time [date, head]
        
        ! Top boundary - configuration
        integer :: swpondmx = 0                !! Time-dependent pondmx switch
        integer :: swredu = 0                  !! Soil evaporation reduction switch
        
        ! Top boundary - values
        real(8) :: pondmx = 0.0d0              !! Max ponding depth (cm)
        real(8) :: hatm = 0.0d0                !! Atmospheric pressure head (cm)
        real(8) :: hsurf = 0.0d0               !! Surface pressure head (cm)
        real(8) :: rsro = 0.0d0                !! Runoff resistance (d)
        real(8) :: rsroexp = 0.0d0             !! Runoff exponent (-)
        real(8) :: runon = 0.0d0               !! Runon flux (cm/d)
        real(8) :: runots = 0.0d0              !! Runoff during timestep (cm)
        real(8) :: crunoff = 0.0d0             !! Cumulative runoff
        real(8) :: crunon = 0.0d0              !! Cumulative runon
        real(8) :: iruno = 0.0d0               !! Intermediate runoff
        real(8) :: irunon = 0.0d0              !! Intermediate runon
        
        ! Top boundary - surface/ponding
        real(8) :: qtop = 0.0d0                !! Surface flux (cm/d)
        real(8) :: q0 = 0.0d0                  !! Net surface flux (precip-evap) (cm/d)
        real(8) :: h0max = 0.0d0               !! Max ponding without runoff (cm)
        real(8) :: k1max = 0.0d0               !! Max conductivity at surface (cm/d)
        real(8) :: QMpLatSs = 0.0d0            !! Lateral inflow to macropores at surface (cm/d)
        
        ! Runon/ponding tables (allocated in init)
        real(8), allocatable :: runonarr(:)    !! Daily runon values (cm/d)
        real(8), allocatable :: pondmxtab(:)   !! Time-dependent pondmx table [date, pondmx]
        
        ! Inundation
        real(8) :: cinund = 0.0d0              !! Cumulative inundation (cm)
        
        ! Flags
        logical :: FlRunoff = .false.          !! Runoff potential possible
        logical :: flrunon = .false.           !! Runon exists
        logical :: ftoph = .false.             !! Pressure head prescribed at surface
   end type boundary_state_t

contains

    !> Initialize the boundary state with optional array sizing.
    !!
    !! Allocates prescribed boundary condition tables and sets all values
    !! to canonical defaults. Array sizes can be customized via optional
    !! parameters, otherwise sensible defaults are used.
    !!
    !! @param[inout] state Boundary state container
    !! @param[in] nday Optional: number of days for runon array (default: MADAY)
    !! @param[in] nbbc Optional: size for bottom BC tables (default: 2*MABBC)
    !! @param[in] npondmx Optional: size for pondmx table (default: 2*MAIRG)
    subroutine boundary_state_init(state, nday, nbbc, npondmx)
    use swap_array_dimensions, only: MADAY, MABBC, MAIRG
        
        type(boundary_state_t), intent(inout) :: state
        integer, intent(in), optional :: nday
        integer, intent(in), optional :: nbbc
        integer, intent(in), optional :: npondmx
        
        integer :: nd, nb, np
        
        ! Determine array sizes
        nd = MADAY
        if (present(nday)) nd = nday
        
        nb = 2 * MABBC
        if (present(nbbc)) nb = nbbc
        
        np = 2 * MAIRG
        if (present(npondmx)) np = npondmx
        
        ! Allocate tables with specified sizes
        allocate(state%gwltab(nb))
        allocate(state%haqtab(nb))
        allocate(state%qbotab(nb))
        allocate(state%hbotab(nb))
        allocate(state%runonarr(nd))
        allocate(state%pondmxtab(np))
        
        ! Reset everything to canonical defaults
        call boundary_state_reset_all(state)
        
    end subroutine boundary_state_init

    !> Finalize the boundary state, deallocating arrays.
    subroutine boundary_state_finalize(state)
        type(boundary_state_t), intent(inout) :: state
        
        ! Deallocate bottom BC tables
        if (allocated(state%gwltab)) deallocate(state%gwltab)
        if (allocated(state%haqtab)) deallocate(state%haqtab)
        if (allocated(state%qbotab)) deallocate(state%qbotab)
        if (allocated(state%hbotab)) deallocate(state%hbotab)
        
        ! Deallocate top BC tables
        if (allocated(state%runonarr)) deallocate(state%runonarr)
        if (allocated(state%pondmxtab)) deallocate(state%pondmxtab)
        
        ! Reset all scalar values for safety
        call boundary_state_reset_all(state)
        
    end subroutine boundary_state_finalize

    !> Reset all boundary state fields to canonical defaults.
    !!
    !! This routine is used by both initialization and finalization to avoid
    !! hidden carry-over state between model runs and instances.
    subroutine boundary_state_reset_all(state)
        type(boundary_state_t), intent(inout) :: state
        
        ! Bottom boundary - configuration
        state%swbotb = 0
        state%swbotb3Impl = 0
        state%SwBotb3ResVert = 0
        state%swqhbot = 0
        state%swcofqhc = 0
        state%sw2 = 0
        state%sw3 = 0
        state%sw4 = 0
        
        ! Bottom boundary - values
        state%qbot = 0.0d0
        state%qbot_nonfrozen = 0.0d0
        state%hbot = 0.0d0
        state%deepgw = 0.0d0
        
        ! Aquifer parameters
        state%aqave = 0.0d0
        state%aqamp = 0.0d0
        state%aqper = 0.0d0
        state%aqtmax = 0.0d0
        state%rimlay = 0.0d0
        state%hdrain = 0.0d0
        state%shape = 0.0d0
        
        ! Sine function parameters
        state%sinave = 0.0d0
        state%sinamp = 0.0d0
        state%sinmax = 0.0d0
        
        ! Flux-head relationships
        state%cofqha = 0.0d0
        state%cofqhb = 0.0d0
        state%cofqhc = 0.0d0
        
        ! Lysimeter
        state%hplate = 0.0d0
        
        ! Reset prescribed tables (if allocated)
        if (allocated(state%gwltab)) state%gwltab = 0.0d0
        if (allocated(state%haqtab)) state%haqtab = 0.0d0
        if (allocated(state%qbotab)) state%qbotab = 0.0d0
        if (allocated(state%hbotab)) state%hbotab = 0.0d0
        
        ! Top boundary - configuration
        state%swpondmx = 0
        state%swredu = 0
        
        ! Top boundary - values
        state%pondmx = 0.0d0
        state%hatm = 0.0d0
        state%hsurf = 0.0d0
        state%rsro = 0.0d0
        state%rsroexp = 0.0d0
        state%runon = 0.0d0
        state%runots = 0.0d0
        
        ! Top boundary - surface/ponding
        state%qtop = 0.0d0
        state%q0 = 0.0d0
        state%h0max = 0.0d0
        state%k1max = 0.0d0
        state%QMpLatSs = 0.0d0
        
        ! Reset tables
        if (allocated(state%runonarr)) state%runonarr = 0.0d0
        if (allocated(state%pondmxtab)) state%pondmxtab = 0.0d0
        
        ! Inundation
        state%cinund = 0.0d0
        
        ! Reset cumulative and intermediate fluxes
        call boundary_state_reset_cumulative(state)
        call boundary_state_reset_intermediate(state)
        
        ! Flags
        state%FlRunoff = .false.
        state%flrunon = .false.
        state%ftoph = .false.
        
    end subroutine boundary_state_reset_all

    !> Reset cumulative flux accumulators (called when flzerocumu=.true.)
    !!
    !! Resets long-term cumulative fluxes tracked since simulation start
    !! or last reset event.
    !!
    !! @param[inout] state Boundary state container
    subroutine boundary_state_reset_cumulative(state)
        type(boundary_state_t), intent(inout) :: state
        
        ! Bottom boundary cumulative
        state%cqbot = 0.0d0
        state%cqbotdo = 0.0d0
        state%cqbotup = 0.0d0
        
        ! Top boundary cumulative
        state%crunoff = 0.0d0
        state%crunon = 0.0d0
        state%cinund = 0.0d0
        
    end subroutine boundary_state_reset_cumulative

    !> Reset intermediate flux accumulators (called when flzerointr=.true.)
    !!
    !! Resets output-interval accumulators tracked between output times.
    !!
    !! @param[inout] state Boundary state container
    subroutine boundary_state_reset_intermediate(state)
        type(boundary_state_t), intent(inout) :: state
        
        ! Bottom boundary intermediate
        state%iqbot = 0.0d0
        
        ! Top boundary intermediate
        state%iruno = 0.0d0
        state%irunon = 0.0d0
        
    end subroutine boundary_state_reset_intermediate

end module boundary_state_mod