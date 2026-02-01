! ==============================================================================
! SWAP State Sync Test Program
! ==============================================================================
! Tests the synchronization of variables between the global Variables module
! and the swap_state_t type. This verifies that SAVE variables that were
! migrated to Variables module are properly captured and restored.
!
! Key tests:
!   1. Headcalc SAVE variables (flwarn_hc, iwarn_hc, nstep_hc)
!   2. Mass balance file unit (dev_cmb)
!   3. Bidirectional sync (variables <-> state)
!
! Author: SWAP Development Team
! Date: 2026-01-31
! ==============================================================================

program test_state_sync
    use swap_log
    use swap_state_mod
    use swap_state_sync
    use variables
    implicit none
    
    ! Test state objects
    type(swap_state_t) :: state1, state2
    
    ! Test parameters
    integer, parameter :: NUMNOD_TEST = 40
    integer, parameter :: NUMLAY_TEST = 5
    integer, parameter :: NRLEVS_TEST = 2
    integer, parameter :: NCROP_TEST = 1
    
    ! Test counters
    integer :: n_tests = 0
    integer :: n_passed = 0
    integer :: n_failed = 0
    
    ! Initialize logging
    call log_init(log_level=LOGLEVEL_INFO, log_file='test_sync.log')
    
    call log_info('test', '========================================')
    call log_info('test', 'SWAP State Sync Tests')
    call log_info('test', '========================================')
    
    ! Initialize global variables module dimensions
    numnod = NUMNOD_TEST
    numlay = NUMLAY_TEST
    nrlevs = NRLEVS_TEST
    
    ! Run test suites
    call test_headcalc_variables_to_state()
    call test_headcalc_state_to_variables()
    call test_headcalc_roundtrip()
    call test_multiple_instance_isolation()
    
    ! Summary
    call log_info('test', '========================================')
    call log_info('test', 'Test Summary')
    call log_info('test', '========================================')
    write(*, '(A,I3)') 'Total tests: ', n_tests
    write(*, '(A,I3)') 'Passed: ', n_passed
    write(*, '(A,I3)') 'Failed: ', n_failed
    
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
            call log_warn('assert', 'FAIL: ' // trim(test_name))
        end if
    end subroutine

    subroutine assert_equal_int(expected, actual, test_name)
        integer, intent(in) :: expected, actual
        character(len=*), intent(in) :: test_name
        
        call assert_true(expected == actual, test_name)
    end subroutine

    subroutine assert_equal_real(expected, actual, tol, test_name)
        real(8), intent(in) :: expected, actual, tol
        character(len=*), intent(in) :: test_name
        
        call assert_true(abs(expected - actual) < tol, test_name)
    end subroutine

    ! =========================================================================
    ! Test: Headcalc variables -> state sync
    ! =========================================================================
    subroutine test_headcalc_variables_to_state()
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Headcalc Variables to State ---')
        
        ! Initialize state
        call swap_state_init(state, numnod=NUMNOD_TEST, numlay=NUMLAY_TEST, &
                            nrlevs=NRLEVS_TEST, ncrop=NCROP_TEST)
        
        ! Set global variables (formerly SAVE variables in headcalc)
        flwarn_hc = .false.
        iwarn_hc = 3
        nstep_hc = 50
        
        ! Sync from variables to state
        call soil_state_from_variables(state%soil, NUMNOD_TEST, NUMLAY_TEST)
        
        ! Verify state captured the values
        call assert_true(.not. state%soil%flwarn_hc, 'flwarn_hc synced to state')
        call assert_equal_int(3, state%soil%iwarn_hc, 'iwarn_hc synced to state')
        call assert_equal_int(50, state%soil%nstep_hc, 'nstep_hc synced to state')
        
    end subroutine

    ! =========================================================================
    ! Test: State -> headcalc variables sync
    ! =========================================================================
    subroutine test_headcalc_state_to_variables()
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: State to Headcalc Variables ---')
        
        ! Initialize state
        call swap_state_init(state, numnod=NUMNOD_TEST, numlay=NUMLAY_TEST, &
                            nrlevs=NRLEVS_TEST, ncrop=NCROP_TEST)
        
        ! Set state values directly
        state%soil%flwarn_hc = .false.
        state%soil%iwarn_hc = 7
        state%soil%nstep_hc = 200
        
        ! Reset global variables to different values
        flwarn_hc = .true.
        iwarn_hc = 0
        nstep_hc = 0
        
        ! Sync from state to variables
        call soil_state_to_variables(state%soil, NUMNOD_TEST, NUMLAY_TEST)
        
        ! Verify global variables were updated
        call assert_true(.not. flwarn_hc, 'flwarn_hc restored from state')
        call assert_equal_int(7, iwarn_hc, 'iwarn_hc restored from state')
        call assert_equal_int(200, nstep_hc, 'nstep_hc restored from state')
        
    end subroutine

    ! =========================================================================
    ! Test: Roundtrip sync (variables -> state -> variables)
    ! =========================================================================
    subroutine test_headcalc_roundtrip()
        type(swap_state_t) :: state
        logical :: orig_flwarn
        integer :: orig_iwarn, orig_nstep
        
        call log_info('test', '--- Test: Headcalc Roundtrip Sync ---')
        
        ! Initialize state
        call swap_state_init(state, numnod=NUMNOD_TEST, numlay=NUMLAY_TEST, &
                            nrlevs=NRLEVS_TEST, ncrop=NCROP_TEST)
        
        ! Set original values in global variables
        flwarn_hc = .false.
        iwarn_hc = 42
        nstep_hc = 999
        
        ! Save originals
        orig_flwarn = flwarn_hc
        orig_iwarn = iwarn_hc
        orig_nstep = nstep_hc
        
        ! Sync to state
        call soil_state_from_variables(state%soil, NUMNOD_TEST, NUMLAY_TEST)
        
        ! Corrupt global variables
        flwarn_hc = .true.
        iwarn_hc = 0
        nstep_hc = 0
        
        ! Sync back from state
        call soil_state_to_variables(state%soil, NUMNOD_TEST, NUMLAY_TEST)
        
        ! Verify original values restored
        call assert_true(flwarn_hc .eqv. orig_flwarn, 'flwarn_hc roundtrip preserved')
        call assert_equal_int(orig_iwarn, iwarn_hc, 'iwarn_hc roundtrip preserved')
        call assert_equal_int(orig_nstep, nstep_hc, 'nstep_hc roundtrip preserved')
        
    end subroutine

    ! =========================================================================
    ! Test: Multiple instance isolation
    ! =========================================================================
    subroutine test_multiple_instance_isolation()
        type(swap_state_t) :: state_a, state_b
        
        call log_info('test', '--- Test: Multiple Instance Isolation ---')
        
        ! Initialize two states
        call swap_state_init(state_a, numnod=NUMNOD_TEST, numlay=NUMLAY_TEST, &
                            nrlevs=NRLEVS_TEST, ncrop=NCROP_TEST)
        call swap_state_init(state_b, numnod=NUMNOD_TEST, numlay=NUMLAY_TEST, &
                            nrlevs=NRLEVS_TEST, ncrop=NCROP_TEST)
        
        ! === Simulate Instance A running ===
        flwarn_hc = .false.
        iwarn_hc = 10
        nstep_hc = 100
        
        ! Capture state A
        call soil_state_from_variables(state_a%soil, NUMNOD_TEST, NUMLAY_TEST)
        
        ! === Simulate Instance B running ===
        flwarn_hc = .true.
        iwarn_hc = 20
        nstep_hc = 200
        
        ! Capture state B
        call soil_state_from_variables(state_b%soil, NUMNOD_TEST, NUMLAY_TEST)
        
        ! Verify states are different
        call assert_true(state_a%soil%flwarn_hc .neqv. state_b%soil%flwarn_hc, &
                        'States have different flwarn_hc')
        call assert_true(state_a%soil%iwarn_hc /= state_b%soil%iwarn_hc, &
                        'States have different iwarn_hc')
        call assert_true(state_a%soil%nstep_hc /= state_b%soil%nstep_hc, &
                        'States have different nstep_hc')
        
        ! === Switch back to Instance A ===
        call soil_state_to_variables(state_a%soil, NUMNOD_TEST, NUMLAY_TEST)
        
        ! Verify Instance A values restored
        call assert_true(.not. flwarn_hc, 'Instance A flwarn_hc restored')
        call assert_equal_int(10, iwarn_hc, 'Instance A iwarn_hc restored')
        call assert_equal_int(100, nstep_hc, 'Instance A nstep_hc restored')
        
        ! === Switch to Instance B ===
        call soil_state_to_variables(state_b%soil, NUMNOD_TEST, NUMLAY_TEST)
        
        ! Verify Instance B values restored
        call assert_true(flwarn_hc, 'Instance B flwarn_hc restored')
        call assert_equal_int(20, iwarn_hc, 'Instance B iwarn_hc restored')
        call assert_equal_int(200, nstep_hc, 'Instance B nstep_hc restored')
        
    end subroutine

end program test_state_sync
