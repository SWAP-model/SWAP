! ==============================================================================
! SWAP State Module Test Program
! ==============================================================================
! Standalone test program to verify the swap_state_mod functionality.
! Run with: ./test_state
!
! This test:
!   1. Initializes the logging system
!   2. Creates a state object with realistic dimensions
!   3. Verifies arrays are allocated correctly
!   4. Tests setting and getting values
!   5. Tests I/O handles initialization
!
! Author: SWAP Development Team
! Date: 2026-01-31
! ==============================================================================

program test_state
    use swap_log
    use swap_state_mod
    implicit none
    
    ! Test state objects
    type(swap_state_t) :: state
    type(io_handles_t) :: io
    
    ! Test parameters (similar to hupselbrook case)
    integer, parameter :: NUMNOD = 40    ! 40 compartments
    integer, parameter :: NUMLAY = 5     ! 5 soil layers
    integer, parameter :: NRLEVS = 2     ! 2 drainage levels
    integer, parameter :: NCROP = 1      ! 1 crop
    
    ! Test counters
    integer :: n_tests = 0
    integer :: n_passed = 0
    integer :: n_failed = 0
    
    ! Initialize logging with debug level
    call log_init(log_level=LOGLEVEL_DEBUG, log_file='test_state.log', timestamps=.true.)
    
    call log_info('test', '========================================')
    call log_info('test', 'SWAP State Module Tests')
    call log_info('test', '========================================')
    
    ! Run test suites
    call test_state_initialization()
    call test_soil_state()
    call test_time_state()
    call test_atmosphere_state()
    call test_drainage_state()
    call test_io_handles()
    call test_multiple_instances()
    
    ! Summary
    call log_info('test', '========================================')
    call log_info('test', 'Test Summary')
    call log_info('test', '========================================')
    call log_info('test', 'Total tests: ' // trim(to_str(n_tests)))
    call log_info('test', 'Passed: ' // trim(to_str(n_passed)))
    call log_info('test', 'Failed: ' // trim(to_str(n_failed)))
    
    if (n_failed > 0) then
        call log_error('test', 'SOME TESTS FAILED!')
        call log_close()
        stop 1
    else
        call log_info('test', 'ALL TESTS PASSED!')
    end if
    
    call log_close()
    
contains

    subroutine assert_true(condition, test_name)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: test_name
        
        n_tests = n_tests + 1
        if (condition) then
            n_passed = n_passed + 1
            call log_debug('assert', 'PASS: ' // trim(test_name))
        else
            n_failed = n_failed + 1
            call log_error('assert', 'FAIL: ' // trim(test_name))
        end if
    end subroutine
    
    subroutine assert_equal_int(expected, actual, test_name)
        integer, intent(in) :: expected, actual
        character(len=*), intent(in) :: test_name
        
        n_tests = n_tests + 1
        if (expected == actual) then
            n_passed = n_passed + 1
            call log_debug('assert', 'PASS: ' // trim(test_name))
        else
            n_failed = n_failed + 1
            call log_error('assert', 'FAIL: ' // trim(test_name) // &
                          ' (expected=' // trim(to_str(expected)) // &
                          ', actual=' // trim(to_str(actual)) // ')')
        end if
    end subroutine
    
    subroutine assert_equal_real(expected, actual, tolerance, test_name)
        real(8), intent(in) :: expected, actual, tolerance
        character(len=*), intent(in) :: test_name
        
        n_tests = n_tests + 1
        if (abs(expected - actual) <= tolerance) then
            n_passed = n_passed + 1
            call log_debug('assert', 'PASS: ' // trim(test_name))
        else
            n_failed = n_failed + 1
            call log_error('assert', 'FAIL: ' // trim(test_name) // &
                          ' (expected=' // trim(to_str(expected)) // &
                          ', actual=' // trim(to_str(actual)) // ')')
        end if
    end subroutine
    
    ! ==========================================================================
    ! Test: State Initialization
    ! ==========================================================================
    subroutine test_state_initialization()
        call log_info('test', '--- Test: State Initialization ---')
        
        ! Initialize state
        call swap_state_init(state, NUMNOD, NUMLAY, NRLEVS, NCROP)
        
        ! Verify initialization flag
        call assert_true(state%initialized, 'state%initialized is true')
        
        ! Verify dimensions stored correctly
        call assert_equal_int(NUMNOD, state%numnod, 'state%numnod')
        call assert_equal_int(NUMLAY, state%numlay, 'state%numlay')
        call assert_equal_int(NRLEVS, state%nrlevs, 'state%nrlevs')
        call assert_equal_int(NCROP, state%ncrop, 'state%ncrop')
        
    end subroutine
    
    ! ==========================================================================
    ! Test: Soil State
    ! ==========================================================================
    subroutine test_soil_state()
        call log_info('test', '--- Test: Soil State ---')
        
        ! Verify arrays are allocated
        call assert_true(allocated(state%soil%h), 'soil%h allocated')
        call assert_true(allocated(state%soil%theta), 'soil%theta allocated')
        call assert_true(allocated(state%soil%k), 'soil%k allocated')
        call assert_true(allocated(state%soil%z), 'soil%z allocated')
        call assert_true(allocated(state%soil%dz), 'soil%dz allocated')
        call assert_true(allocated(state%soil%bdens), 'soil%bdens allocated')
        call assert_true(allocated(state%soil%paramvg), 'soil%paramvg allocated')
        
        ! Verify array sizes
        call assert_equal_int(NUMNOD, size(state%soil%h), 'size(soil%h)')
        call assert_equal_int(NUMNOD, size(state%soil%theta), 'size(soil%theta)')
        call assert_equal_int(NUMNOD+1, size(state%soil%k), 'size(soil%k)')
        call assert_equal_int(NUMLAY, size(state%soil%bdens), 'size(soil%bdens)')
        
        ! Verify initial values are zero
        call assert_equal_real(0.0d0, state%soil%h(1), 1.0d-12, 'soil%h(1) = 0')
        call assert_equal_real(0.0d0, state%soil%gwl, 1.0d-12, 'soil%gwl = 0')
        
        ! Test setting values
        state%soil%h(1) = -100.0d0
        state%soil%theta(1) = 0.35d0
        state%soil%gwl = -150.0d0
        
        call assert_equal_real(-100.0d0, state%soil%h(1), 1.0d-12, 'soil%h(1) after set')
        call assert_equal_real(0.35d0, state%soil%theta(1), 1.0d-12, 'soil%theta(1) after set')
        call assert_equal_real(-150.0d0, state%soil%gwl, 1.0d-12, 'soil%gwl after set')
        
    end subroutine
    
    ! ==========================================================================
    ! Test: Time State
    ! ==========================================================================
    subroutine test_time_state()
        call log_info('test', '--- Test: Time State ---')
        
        ! Test default values
        call assert_equal_real(0.0d0, state%time%dt, 1.0d-12, 'time%dt initial')
        call assert_equal_real(0.0d0, state%time%t, 1.0d-12, 'time%t initial')
        call assert_true(.not. state%time%flrunend, 'time%flrunend false initially')
        
        ! Test setting values
        state%time%dt = 0.1d0
        state%time%dtmax = 1.0d0
        state%time%dtmin = 0.001d0
        state%time%t = 100.0d0
        state%time%tstart = 0.0d0
        state%time%tend = 365.0d0
        state%time%daynr = 100
        state%time%iyear = 2026
        
        call assert_equal_real(0.1d0, state%time%dt, 1.0d-12, 'time%dt after set')
        call assert_equal_int(100, state%time%daynr, 'time%daynr after set')
        call assert_equal_int(2026, state%time%iyear, 'time%iyear after set')
        
    end subroutine
    
    ! ==========================================================================
    ! Test: Atmosphere State
    ! ==========================================================================
    subroutine test_atmosphere_state()
        call log_info('test', '--- Test: Atmosphere State ---')
        
        ! Test setting meteorological values
        state%atm%tav = 15.0d0
        state%atm%tmn = 10.0d0
        state%atm%tmx = 20.0d0
        state%atm%rad = 15000000.0d0  ! J/m2/d
        state%atm%rh = 0.75d0
        state%atm%peva = 0.3d0
        state%atm%ptra = 0.2d0
        
        call assert_equal_real(15.0d0, state%atm%tav, 1.0d-12, 'atm%tav')
        call assert_equal_real(0.75d0, state%atm%rh, 1.0d-12, 'atm%rh')
        call assert_equal_real(0.3d0, state%atm%peva, 1.0d-12, 'atm%peva')
        
    end subroutine
    
    ! ==========================================================================
    ! Test: Drainage State
    ! ==========================================================================
    subroutine test_drainage_state()
        call log_info('test', '--- Test: Drainage State ---')
        
        ! Verify arrays allocated
        call assert_true(allocated(state%drain%qdrain), 'drain%qdrain allocated')
        call assert_true(allocated(state%drain%qdra), 'drain%qdra allocated')
        
        ! Verify sizes
        call assert_equal_int(NRLEVS, size(state%drain%qdrain), 'size(drain%qdrain)')
        call assert_equal_int(NRLEVS, state%drain%nrlevs, 'drain%nrlevs')
        
        ! Test setting values
        state%drain%qdrain(1) = 0.05d0
        state%drain%qdrtot = 0.1d0
        
        call assert_equal_real(0.05d0, state%drain%qdrain(1), 1.0d-12, 'drain%qdrain(1)')
        call assert_equal_real(0.1d0, state%drain%qdrtot, 1.0d-12, 'drain%qdrtot')
        
    end subroutine
    
    ! ==========================================================================
    ! Test: I/O Handles
    ! ==========================================================================
    subroutine test_io_handles()
        call log_info('test', '--- Test: I/O Handles ---')
        
        ! Initialize I/O handles
        call io_handles_init(io)
        
        ! Verify all handles are invalid (-1) initially
        call assert_equal_int(-1, io%swp_unit, 'io%swp_unit = -1')
        call assert_equal_int(-1, io%log_unit, 'io%log_unit = -1')
        call assert_equal_int(-1, io%bal_unit, 'io%bal_unit = -1')
        call assert_true(.not. io%files_open, 'io%files_open = false')
        
    end subroutine
    
    ! ==========================================================================
    ! Test: Multiple Instances
    ! ==========================================================================
    subroutine test_multiple_instances()
        type(swap_state_t) :: state1, state2
        
        call log_info('test', '--- Test: Multiple Instances ---')
        
        ! Initialize two independent states
        call swap_state_init(state1, 20, 3, 1, 1)
        call swap_state_init(state2, 50, 7, 3, 2)
        
        ! Verify they have different dimensions
        call assert_equal_int(20, state1%numnod, 'state1%numnod')
        call assert_equal_int(50, state2%numnod, 'state2%numnod')
        call assert_equal_int(3, state1%numlay, 'state1%numlay')
        call assert_equal_int(7, state2%numlay, 'state2%numlay')
        
        ! Verify arrays are independent
        state1%soil%h(1) = -100.0d0
        state2%soil%h(1) = -200.0d0
        
        call assert_equal_real(-100.0d0, state1%soil%h(1), 1.0d-12, 'state1 h independent')
        call assert_equal_real(-200.0d0, state2%soil%h(1), 1.0d-12, 'state2 h independent')
        
        call log_info('test', 'Multiple instances are independent - key for multi-core!')
        
    end subroutine

end program test_state
