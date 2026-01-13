! swap_state.f90
! State object definition for SWAP model refactoring
! 
! This module provides an explicit state object to replace SAVE variables,
! enabling multiple independent SWAP instances in a single process.
!
! Author: Zawadzki M.
! Date: 2026-01-12
! Phase: Proof-of-concept - minimal critical variables

module swap_state
    use iso_c_binding, only: c_int, c_double, c_ptr
    implicit none
    
    !-------------------------------------------------------------------
    ! SWAP State Type
    !-------------------------------------------------------------------
    ! This type will eventually contain ALL variables currently in the
    ! variables module. For proof-of-concept, we start with critical ones.
    
    type :: swap_state_t
        ! === BMI tracking ===
        logical :: initialized = .false.
        logical :: memory_mode = .false.
        integer :: iCaller = 0
        
        ! === Time variables (critical for BMI) ===
        real(8) :: t                  ! Time since start of calendar year (T)
        real(8) :: t1900              ! Time since 1900 (T)
        real(8) :: tcum               ! Time since start of simulation (T)
        real(8) :: dt                 ! Time step (T)
        real(8) :: tend               ! End date of simulation run
        real(8) :: tstart             ! Start date of simulation run
        
        ! === Control flags ===
        logical :: flrunend = .false.   ! Flag indicating end of run
        logical :: fldayend = .false.   ! Flag indicating end of day
        logical :: fldaystart = .false. ! Flag indicating start of day
        
        ! === Grid configuration ===
        integer :: numnod = 0          ! Number of soil nodes
        
        ! === State arrays (allocatable for dynamic sizing) ===
        ! These will be allocated during initialization
        real(8), allocatable :: theta(:)     ! Water content
        real(8), allocatable :: h(:)         ! Pressure head
        
        ! === Project info ===
        character(len=80) :: project = ''
        character(len=80) :: swpfile = ''
        
        ! === CNmethod state (meteoday.f90) ===
        integer :: cn_Nod10 = 0
        integer :: cn_iCN = 0
        real(8) :: cn_Z10 = 0.0d0
        
    end type swap_state_t
    
contains

    !-------------------------------------------------------------------
    ! Create and initialize a new SWAP state object
    !-------------------------------------------------------------------
    function swap_state_create() result(state)
        type(swap_state_t) :: state
        
        ! Initialize with default values
        state%initialized = .false.
        state%memory_mode = .false.
        state%iCaller = 0
        
        state%t = 0.0d0
        state%t1900 = 0.0d0
        state%tcum = 0.0d0
        state%dt = 0.0d0
        state%tend = 0.0d0
        state%tstart = 0.0d0
        
        state%flrunend = .false.
        state%fldayend = .false.
        state%fldaystart = .false.
        
        state%numnod = 0
        
        state%project = ''
        state%swpfile = ''
        
        ! CNmethod state
        state%cn_Nod10 = 0
        state%cn_iCN = 0
        state%cn_Z10 = 0.0d0
        
        ! Arrays will be allocated when numnod is known
        
    end function swap_state_create
    
    !-------------------------------------------------------------------
    ! Allocate state arrays after grid size is known
    !-------------------------------------------------------------------
    subroutine swap_state_allocate_arrays(state, numnod)
        type(swap_state_t), intent(inout) :: state
        integer, intent(in) :: numnod
        
        state%numnod = numnod
        
        ! Allocate arrays
        if (.not. allocated(state%theta)) then
            allocate(state%theta(numnod))
            state%theta = 0.0d0
        end if
        
        if (.not. allocated(state%h)) then
            allocate(state%h(numnod))
            state%h = 0.0d0
        end if
        
    end subroutine swap_state_allocate_arrays
    
    !-------------------------------------------------------------------
    ! Destroy state object and deallocate memory
    !-------------------------------------------------------------------
    subroutine swap_state_destroy(state)
        type(swap_state_t), intent(inout) :: state
        
        ! Deallocate arrays
        if (allocated(state%theta)) deallocate(state%theta)
        if (allocated(state%h)) deallocate(state%h)
        
        ! Reset state
        state%initialized = .false.
        state%numnod = 0
        
    end subroutine swap_state_destroy
    
    !-------------------------------------------------------------------
    ! Copy state values FROM module variables TO state object
    ! Used during transition period - eventually will be removed
    !-------------------------------------------------------------------
    subroutine swap_state_sync_from_globals(state)
        use variables
        type(swap_state_t), intent(inout) :: state
        
        ! Copy time variables
        state%t = t
        state%t1900 = t1900
        state%tcum = tcum
        state%dt = dt
        state%tend = tend
        state%tstart = tstart
        
        ! Copy flags
        state%flrunend = flrunend
        state%fldayend = fldayend
        state%fldaystart = fldaystart
        
        ! Copy project info
        state%project = project
        state%swpfile = swpfile
        
    end subroutine swap_state_sync_from_globals
    
    !-------------------------------------------------------------------
    ! Copy state values FROM state object TO module variables
    ! Used during transition period - eventually will be removed
    !-------------------------------------------------------------------
    subroutine swap_state_sync_to_globals(state)
        use variables
        type(swap_state_t), intent(in) :: state
        
        ! Copy time variables
        t = state%t
        t1900 = state%t1900
        tcum = state%tcum
        dt = state%dt
        tend = state%tend
        tstart = state%tstart
        
        ! Copy flags
        flrunend = state%flrunend
        fldayend = state%fldayend
        fldaystart = state%fldaystart
        
        ! Copy project info
        project = state%project
        swpfile = state%swpfile
        
    end subroutine swap_state_sync_to_globals
    
end module swap_state

!-------------------------------------------------------------------
! Module-level state object for SWAP core
! (Makes state accessible to all SWAP subroutines without passing)
!-------------------------------------------------------------------
module swap_core_state
    use swap_state, only: swap_state_t
    implicit none
    type(swap_state_t), save :: core_state
end module swap_core_state
