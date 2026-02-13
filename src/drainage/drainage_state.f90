!> Drainage State Module
!!
!! This module defines the drainage state for SWAP.
!! Manages drainage/infiltration fluxes, level-wise drainage configuration,
!! discharge-layer controls, and interflow-related parameters.
!!
!! Key components:
!! - Drainage level configuration and geometry
!! - Drainage and infiltration flux arrays
!! - Cumulative and intermediate drainage tracking
!! - Interflow configuration
!!
!! @author SWAP Development Team
!! @date 2026-02-13
module drainage_state_mod
    implicit none
    private

    public :: drainage_state_t
    public :: drainage_state_init
    public :: drainage_state_finalize
    public :: drainage_state_reset_cumulative
    public :: drainage_state_reset_intermediate

    type :: drainage_state_t
        ! Drainage fluxes per level
        real(8), allocatable :: qdrain(:)       !! Drainage flux per level (cm/d)
        real(8), allocatable :: cqdrain(:)      !! Cumulative drainage per level
        real(8), allocatable :: cqdrainin(:)    !! Cumulative infiltration per level
        real(8), allocatable :: cqdrainout(:)   !! Cumulative drainage out per level
        real(8), allocatable :: drainl(:)       !! Drainage level depth

        ! Spatial drainage arrays
        real(8), allocatable :: qdra(:,:)       !! Drainage flux (level,node)
        real(8), allocatable :: inqdra(:,:)     !! Intermediate drainage (level,node)
        real(8), allocatable :: inqdra_in(:,:)  !! Intermediate infiltration flux per level/node
        real(8), allocatable :: inqdra_out(:,:) !! Intermediate drainage out per level/node
        real(8), allocatable :: qdraincomp(:)   !! Drainage per compartment

        ! Totals
        real(8) :: qdrtot = 0.0d0               !! Total drainage flux
        real(8) :: iqdra = 0.0d0                !! Intermediate drainage
        real(8) :: cqdra = 0.0d0                !! Cumulative total lateral drainage

        ! Parameters
        integer :: nrlevs = 0                   !! Number of drainage levels
        integer :: nrpri = 0                    !! Number of primary drainage levels
        integer :: dramet = 0                   !! Drainage method switch
        integer :: swdivd = 0                   !! Distribution of drainage over profile
        integer :: swdislay = 0                 !! Discharge layer option
        real(8) :: basegw = 0.0d0               !! Impervious layer depth
        real(8) :: entres = 0.0d0               !! Drain entry resistance
        real(8) :: shape = 0.0d0                !! Shape factor

        ! Resistances and geometry per level
        real(8), allocatable :: drares(:)       !! Drainage resistance
        real(8), allocatable :: infres(:)       !! Infiltration resistance
        real(8), allocatable :: L(:)            !! Drain spacing
        real(8), allocatable :: wetper(:)       !! Wet perimeter
        real(8), allocatable :: zbotdr(:)       !! Drain bottom depth
        real(8), allocatable :: rdrain(:)       !! Drainage resistance per level
        real(8), allocatable :: rinfi(:)        !! Infiltration resistance per level
        real(8), allocatable :: rentry(:)       !! Entry resistance per level
        real(8), allocatable :: rexit(:)        !! Exit resistance per level
        real(8), allocatable :: gwlinf(:)       !! GWL below which no infiltration
        real(8), allocatable :: widthr(:)       !! Width of drain/channel
        real(8), allocatable :: taludr(:)       !! Talus slope of drain
        integer, allocatable :: swallo(:)       !! Allow drainage/infiltration per level
        integer, allocatable :: swdtyp(:)       !! Drainage type per level (0=channel, 1=tube, 2=interflow)
        integer, allocatable :: swtopdislay(:)  !! Top of discharge layer option
        real(8), allocatable :: zTopDisLay(:)   !! Top of discharge layer
        real(8), allocatable :: fTopDisLay(:)   !! Fraction for discharge layer top

        ! Interflow parameters
        real(8) :: cofintfl = 0.0d0             !! Interflow coefficient
        real(8) :: expintfl = 0.0d0             !! Interflow exponent
        integer :: swnrsrf = 0                  !! Interflow switch
        integer :: SwTopnrsrf = 0               !! Interflow top switch
        real(8) :: rsurfdeep = 0.0d0            !! Deep interflow resistance
        real(8) :: rsurfshallow = 0.0d0         !! Shallow interflow resistance
        real(8) :: FacDpthInf = 0.0d0           !! Factor for infiltration depth
        integer :: Swdivdinf = 0                !! Infiltration distribution switch

        ! Flags
        logical :: fldrain = .false.

        ! Paths
        character(len=16) :: drfil = ''
        character(len=80) :: pathdrain = ''
    end type drainage_state_t

contains

    !> Initialize the drainage state.
    !!
    !! Allocates level-wise and spatial arrays and sets all values to
    !! canonical defaults.
    !!
    !! @param[inout] state Drainage state container
    !! @param[in] nrlevs Number of drainage levels
    !! @param[in] numnod Number of soil nodes
    subroutine drainage_state_init(state, nrlevs, numnod)
        type(drainage_state_t), intent(inout) :: state
        integer, intent(in) :: nrlevs
        integer, intent(in) :: numnod

        ! Drainage fluxes per level
        allocate(state%qdrain(nrlevs))
        allocate(state%cqdrain(nrlevs))
        allocate(state%cqdrainin(nrlevs))
        allocate(state%cqdrainout(nrlevs))
        allocate(state%drainl(nrlevs))

        ! Spatial arrays
        allocate(state%qdra(nrlevs, numnod))
        allocate(state%inqdra(nrlevs, numnod))
        allocate(state%inqdra_in(nrlevs, numnod))
        allocate(state%inqdra_out(nrlevs, numnod))
        allocate(state%qdraincomp(numnod))

        ! Resistances and geometry per level
        allocate(state%drares(nrlevs))
        allocate(state%infres(nrlevs))
        allocate(state%L(nrlevs))
        allocate(state%wetper(nrlevs))
        allocate(state%zbotdr(nrlevs))
        allocate(state%rdrain(nrlevs))
        allocate(state%rinfi(nrlevs))
        allocate(state%rentry(nrlevs))
        allocate(state%rexit(nrlevs))
        allocate(state%gwlinf(nrlevs))
        allocate(state%widthr(nrlevs))
        allocate(state%taludr(nrlevs))
        allocate(state%swallo(nrlevs))
        allocate(state%swdtyp(nrlevs))
        allocate(state%swtopdislay(nrlevs))
        allocate(state%zTopDisLay(nrlevs))
        allocate(state%fTopDisLay(nrlevs))

        call drainage_state_reset_all(state)
        state%nrlevs = nrlevs

    end subroutine drainage_state_init

    !> Finalize the drainage state.
    !!
    !! Deallocates all dynamic arrays and resets scalar values to canonical
    !! defaults to prevent carry-over state between runs.
    !!
    !! @param[inout] state Drainage state container
    subroutine drainage_state_finalize(state)
        type(drainage_state_t), intent(inout) :: state

        ! Flux arrays
        if (allocated(state%qdrain)) deallocate(state%qdrain)
        if (allocated(state%cqdrain)) deallocate(state%cqdrain)
        if (allocated(state%cqdrainin)) deallocate(state%cqdrainin)
        if (allocated(state%cqdrainout)) deallocate(state%cqdrainout)
        if (allocated(state%drainl)) deallocate(state%drainl)

        ! Spatial arrays
        if (allocated(state%qdra)) deallocate(state%qdra)
        if (allocated(state%inqdra)) deallocate(state%inqdra)
        if (allocated(state%inqdra_in)) deallocate(state%inqdra_in)
        if (allocated(state%inqdra_out)) deallocate(state%inqdra_out)
        if (allocated(state%qdraincomp)) deallocate(state%qdraincomp)

        ! Resistance/geometry arrays
        if (allocated(state%drares)) deallocate(state%drares)
        if (allocated(state%infres)) deallocate(state%infres)
        if (allocated(state%L)) deallocate(state%L)
        if (allocated(state%wetper)) deallocate(state%wetper)
        if (allocated(state%zbotdr)) deallocate(state%zbotdr)
        if (allocated(state%rdrain)) deallocate(state%rdrain)
        if (allocated(state%rinfi)) deallocate(state%rinfi)
        if (allocated(state%rentry)) deallocate(state%rentry)
        if (allocated(state%rexit)) deallocate(state%rexit)
        if (allocated(state%gwlinf)) deallocate(state%gwlinf)
        if (allocated(state%widthr)) deallocate(state%widthr)
        if (allocated(state%taludr)) deallocate(state%taludr)
        if (allocated(state%swallo)) deallocate(state%swallo)
        if (allocated(state%swdtyp)) deallocate(state%swdtyp)
        if (allocated(state%swtopdislay)) deallocate(state%swtopdislay)
        if (allocated(state%zTopDisLay)) deallocate(state%zTopDisLay)
        if (allocated(state%fTopDisLay)) deallocate(state%fTopDisLay)

        call drainage_state_reset_all(state)

    end subroutine drainage_state_finalize

    !> Reset all drainage state fields to canonical defaults.
    !!
    !! Used by both initialization and finalization. Scalars are reset,
    !! and allocatable arrays are zeroed when allocated.
    !!
    !! @param[inout] state Drainage state container
    subroutine drainage_state_reset_all(state)
        type(drainage_state_t), intent(inout) :: state

        ! Array fields
        if (allocated(state%qdrain)) state%qdrain = 0.0d0
        if (allocated(state%cqdrain)) state%cqdrain = 0.0d0
        if (allocated(state%cqdrainin)) state%cqdrainin = 0.0d0
        if (allocated(state%cqdrainout)) state%cqdrainout = 0.0d0
        if (allocated(state%drainl)) state%drainl = 0.0d0
        if (allocated(state%qdra)) state%qdra = 0.0d0
        if (allocated(state%inqdra)) state%inqdra = 0.0d0
        if (allocated(state%inqdra_in)) state%inqdra_in = 0.0d0
        if (allocated(state%inqdra_out)) state%inqdra_out = 0.0d0
        if (allocated(state%qdraincomp)) state%qdraincomp = 0.0d0
        if (allocated(state%drares)) state%drares = 0.0d0
        if (allocated(state%infres)) state%infres = 0.0d0
        if (allocated(state%L)) state%L = 0.0d0
        if (allocated(state%wetper)) state%wetper = 0.0d0
        if (allocated(state%zbotdr)) state%zbotdr = 0.0d0
        if (allocated(state%rdrain)) state%rdrain = 0.0d0
        if (allocated(state%rinfi)) state%rinfi = 0.0d0
        if (allocated(state%rentry)) state%rentry = 0.0d0
        if (allocated(state%rexit)) state%rexit = 0.0d0
        if (allocated(state%gwlinf)) state%gwlinf = 0.0d0
        if (allocated(state%widthr)) state%widthr = 0.0d0
        if (allocated(state%taludr)) state%taludr = 0.0d0
        if (allocated(state%swallo)) state%swallo = 0
        if (allocated(state%swdtyp)) state%swdtyp = 0
        if (allocated(state%swtopdislay)) state%swtopdislay = 0
        if (allocated(state%zTopDisLay)) state%zTopDisLay = 0.0d0
        if (allocated(state%fTopDisLay)) state%fTopDisLay = 0.0d0

        ! Totals and parameters
        state%qdrtot = 0.0d0
        state%nrlevs = 0
        state%nrpri = 0
        state%dramet = 0
        state%swdivd = 0
        state%swdislay = 0
        state%basegw = 0.0d0
        state%entres = 0.0d0
        state%shape = 0.0d0

        ! Interflow parameters
        state%cofintfl = 0.0d0
        state%expintfl = 0.0d0
        state%swnrsrf = 0
        state%SwTopnrsrf = 0
        state%rsurfdeep = 0.0d0
        state%rsurfshallow = 0.0d0
        state%FacDpthInf = 0.0d0
        state%Swdivdinf = 0

        ! Cumulative and intermediate fields
        call drainage_state_reset_cumulative(state)
        call drainage_state_reset_intermediate(state)

        ! Flags and paths
        state%fldrain = .false.
        state%drfil = ''
        state%pathdrain = ''

    end subroutine drainage_state_reset_all

    !> Reset cumulative drainage accumulators.
    !!
    !! Called when long-term cumulative outputs are reset.
    !!
    !! @param[inout] state Drainage state container
    subroutine drainage_state_reset_cumulative(state)
        type(drainage_state_t), intent(inout) :: state

        state%cqdra = 0.0d0
        if (allocated(state%cqdrain)) state%cqdrain = 0.0d0
        if (allocated(state%cqdrainin)) state%cqdrainin = 0.0d0
        if (allocated(state%cqdrainout)) state%cqdrainout = 0.0d0

    end subroutine drainage_state_reset_cumulative

    !> Reset intermediate drainage accumulators.
    !!
    !! Called when output-interval accumulators are reset.
    !!
    !! @param[inout] state Drainage state container
    subroutine drainage_state_reset_intermediate(state)
        type(drainage_state_t), intent(inout) :: state

        state%iqdra = 0.0d0
        if (allocated(state%inqdra)) state%inqdra = 0.0d0
        if (allocated(state%inqdra_in)) state%inqdra_in = 0.0d0
        if (allocated(state%inqdra_out)) state%inqdra_out = 0.0d0

    end subroutine drainage_state_reset_intermediate

end module drainage_state_mod
  