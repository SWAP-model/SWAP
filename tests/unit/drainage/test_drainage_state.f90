!> @file test_drainage_state.f90
!> @brief Unit tests for drainage state structures
!> 
!> Tests the drainage_state_t and surfacewater_state_t types from swap_state_mod.
!> This test is self-contained and does not require Variables module arrays.
!> Key functionality tested:
!> - Drainage state initialization
!> - All drainage arrays (qdrain, qdra, resistances, etc.)
!> - Surface water state
!> - Multiple drainage instances independence

program test_drainage_state
    use swap_state_mod
    use swap_log
    implicit none
    
    integer :: tests_run, tests_passed
    
    tests_run = 0
    tests_passed = 0
    
    ! Initialize logging
    call log_init(LOGLEVEL_DEBUG)
    call log_info('test', '========================================')
    call log_info('test', 'SWAP Drainage State Tests')
    call log_info('test', '========================================')
    
    ! Run test suites
    call test_drainage_initialization()
    call test_drainage_array_values()
    call test_drainage_resistances()
    call test_drainage_interflow()
    call test_surfacewater_initialization()
    call test_surfacewater_values()
    call test_multiple_drainage_instances()
    
    ! Print summary
    call log_info('test', '========================================')
    call log_info('test', 'Test Summary')
    call log_info('test', '========================================')
    call log_info('test', 'Total tests: ' // to_str(tests_run))
    call log_info('test', 'Passed: ' // to_str(tests_passed))
    call log_info('test', 'Failed: ' // to_str(tests_run - tests_passed))
    
    if (tests_passed == tests_run) then
        call log_info('test', 'ALL TESTS PASSED!')
        stop 0
    else
        call log_info('test', 'SOME TESTS FAILED!')
        stop 1
    end if

contains

    subroutine assert_true(condition, test_name)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: test_name
        
        tests_run = tests_run + 1
        if (condition) then
            tests_passed = tests_passed + 1
            call log_debug('assert', 'PASS: ' // trim(test_name))
        else
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

    subroutine test_drainage_initialization()
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Drainage State Initialization ---')
        
        ! Initialize state with 3 drainage levels and 40 nodes
        call swap_state_init(state, numnod=40, numlay=5, nrlevs=3, ncrop=1)
        
        ! Verify drainage arrays are allocated
        call assert_true(allocated(state%drain%qdrain), 'drain%qdrain allocated')
        call assert_true(allocated(state%drain%cqdrain), 'drain%cqdrain allocated')
        call assert_true(allocated(state%drain%cqdrainin), 'drain%cqdrainin allocated')
        call assert_true(allocated(state%drain%cqdrainout), 'drain%cqdrainout allocated')
        call assert_true(allocated(state%drain%drainl), 'drain%drainl allocated')
        call assert_true(allocated(state%drain%qdra), 'drain%qdra allocated')
        call assert_true(allocated(state%drain%inqdra), 'drain%inqdra allocated')
        call assert_true(allocated(state%drain%inqdra_in), 'drain%inqdra_in allocated')
        call assert_true(allocated(state%drain%inqdra_out), 'drain%inqdra_out allocated')
        call assert_true(allocated(state%drain%qdraincomp), 'drain%qdraincomp allocated')
        
        ! Verify sizes
        call assert_equal_int(3, size(state%drain%qdrain), 'size(drain%qdrain) = 3')
        call assert_equal_int(3, size(state%drain%cqdrain), 'size(drain%cqdrain) = 3')
        call assert_equal_int(3, size(state%drain%drainl), 'size(drain%drainl) = 3')
        call assert_equal_int(3, state%drain%nrlevs, 'drain%nrlevs = 3')
        call assert_equal_int(40, size(state%drain%qdraincomp), 'size(drain%qdraincomp) = 40')
        
        ! Verify 2D array dimensions
        call assert_equal_int(3, size(state%drain%qdra, 1), 'qdra dim 1 = 3')
        call assert_equal_int(40, size(state%drain%qdra, 2), 'qdra dim 2 = 40')
        
    end subroutine test_drainage_initialization
    
    subroutine test_drainage_array_values()
        type(swap_state_t) :: state
        integer :: i
        
        call log_info('test', '--- Test: Drainage Array Values ---')
        
        call swap_state_init(state, numnod=40, numlay=5, nrlevs=3, ncrop=1)
        
        ! Verify initial values are zero
        call assert_equal_real(0.0d0, state%drain%qdrain(1), 1.0d-12, 'qdrain(1) = 0 initially')
        call assert_equal_real(0.0d0, state%drain%qdrtot, 1.0d-12, 'qdrtot = 0 initially')
        call assert_equal_real(0.0d0, state%drain%cqdra, 1.0d-12, 'cqdra = 0 initially')
        
        ! Set values for each drainage level
        do i = 1, 3
            state%drain%qdrain(i) = 0.1d0 * i
            state%drain%cqdrain(i) = 10.0d0 * i
            state%drain%drainl(i) = -50.0d0 - 25.0d0 * i
        end do
        state%drain%qdrtot = 0.6d0  ! sum of 0.1 + 0.2 + 0.3
        
        ! Verify values
        call assert_equal_real(0.1d0, state%drain%qdrain(1), 1.0d-12, 'qdrain(1) = 0.1')
        call assert_equal_real(0.2d0, state%drain%qdrain(2), 1.0d-12, 'qdrain(2) = 0.2')
        call assert_equal_real(0.3d0, state%drain%qdrain(3), 1.0d-12, 'qdrain(3) = 0.3')
        call assert_equal_real(0.6d0, state%drain%qdrtot, 1.0d-12, 'qdrtot = 0.6')
        call assert_equal_real(-75.0d0, state%drain%drainl(1), 1.0d-12, 'drainl(1) = -75')
        call assert_equal_real(-100.0d0, state%drain%drainl(2), 1.0d-12, 'drainl(2) = -100')
        
        ! Test 2D array
        state%drain%qdra(2, 30) = 0.05d0
        call assert_equal_real(0.05d0, state%drain%qdra(2, 30), 1.0d-12, 'qdra(2,30) = 0.05')
        
    end subroutine test_drainage_array_values
    
    subroutine test_drainage_resistances()
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Drainage Resistances ---')
        
        call swap_state_init(state, numnod=40, numlay=5, nrlevs=3, ncrop=1)
        
        ! Verify resistance arrays allocated
        call assert_true(allocated(state%drain%drares), 'drares allocated')
        call assert_true(allocated(state%drain%infres), 'infres allocated')
        call assert_true(allocated(state%drain%rdrain), 'rdrain allocated')
        call assert_true(allocated(state%drain%rinfi), 'rinfi allocated')
        call assert_true(allocated(state%drain%rentry), 'rentry allocated')
        call assert_true(allocated(state%drain%rexit), 'rexit allocated')
        call assert_true(allocated(state%drain%L), 'L allocated')
        call assert_true(allocated(state%drain%wetper), 'wetper allocated')
        call assert_true(allocated(state%drain%zbotdr), 'zbotdr allocated')
        call assert_true(allocated(state%drain%gwlinf), 'gwlinf allocated')
        call assert_true(allocated(state%drain%widthr), 'widthr allocated')
        call assert_true(allocated(state%drain%taludr), 'taludr allocated')
        
        ! Set and verify resistance values
        state%drain%drares(1) = 100.0d0
        state%drain%infres(1) = 200.0d0
        state%drain%L(1) = 2000.0d0
        state%drain%zbotdr(1) = -90.0d0
        
        call assert_equal_real(100.0d0, state%drain%drares(1), 1.0d-12, 'drares(1) = 100')
        call assert_equal_real(200.0d0, state%drain%infres(1), 1.0d-12, 'infres(1) = 200')
        call assert_equal_real(2000.0d0, state%drain%L(1), 1.0d-12, 'L(1) = 2000')
        call assert_equal_real(-90.0d0, state%drain%zbotdr(1), 1.0d-12, 'zbotdr(1) = -90')
        
        ! Test integer arrays
        call assert_true(allocated(state%drain%swallo), 'swallo allocated')
        call assert_true(allocated(state%drain%swdtyp), 'swdtyp allocated')
        
        state%drain%swallo(1) = 1
        state%drain%swdtyp(1) = 2
        call assert_equal_int(1, state%drain%swallo(1), 'swallo(1) = 1')
        call assert_equal_int(2, state%drain%swdtyp(1), 'swdtyp(1) = 2')
        
    end subroutine test_drainage_resistances
    
    subroutine test_drainage_interflow()
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Drainage Interflow Parameters ---')
        
        call swap_state_init(state, numnod=40, numlay=5, nrlevs=3, ncrop=1)
        
        ! Set interflow parameters
        state%drain%swnrsrf = 1
        state%drain%SwTopnrsrf = 1
        state%drain%cofintfl = 0.5d0
        state%drain%expintfl = 1.5d0
        state%drain%rsurfdeep = 500.0d0
        state%drain%rsurfshallow = 50.0d0
        state%drain%FacDpthInf = 0.25d0
        state%drain%Swdivdinf = 1
        
        call assert_equal_int(1, state%drain%swnrsrf, 'swnrsrf = 1')
        call assert_equal_int(1, state%drain%SwTopnrsrf, 'SwTopnrsrf = 1')
        call assert_equal_real(0.5d0, state%drain%cofintfl, 1.0d-12, 'cofintfl = 0.5')
        call assert_equal_real(1.5d0, state%drain%expintfl, 1.0d-12, 'expintfl = 1.5')
        call assert_equal_real(500.0d0, state%drain%rsurfdeep, 1.0d-12, 'rsurfdeep = 500')
        call assert_equal_real(50.0d0, state%drain%rsurfshallow, 1.0d-12, 'rsurfshallow = 50')
        call assert_equal_real(0.25d0, state%drain%FacDpthInf, 1.0d-12, 'FacDpthInf = 0.25')
        call assert_equal_int(1, state%drain%Swdivdinf, 'Swdivdinf = 1')
        
        ! Test discharge layer parameters
        call assert_true(allocated(state%drain%swtopdislay), 'swtopdislay allocated')
        call assert_true(allocated(state%drain%zTopDisLay), 'zTopDisLay allocated')
        call assert_true(allocated(state%drain%fTopDisLay), 'fTopDisLay allocated')
        
        state%drain%swdislay = 2
        state%drain%swtopdislay(1) = 1
        state%drain%zTopDisLay(1) = -30.0d0
        state%drain%fTopDisLay(1) = 0.5d0
        
        call assert_equal_int(2, state%drain%swdislay, 'swdislay = 2')
        call assert_equal_int(1, state%drain%swtopdislay(1), 'swtopdislay(1) = 1')
        call assert_equal_real(-30.0d0, state%drain%zTopDisLay(1), 1.0d-12, 'zTopDisLay(1) = -30')
        call assert_equal_real(0.5d0, state%drain%fTopDisLay(1), 1.0d-12, 'fTopDisLay(1) = 0.5')
        
    end subroutine test_drainage_interflow
    
    subroutine test_surfacewater_initialization()
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Surface Water State Initialization ---')
        
        call swap_state_init(state, numnod=40, numlay=5, nrlevs=3, ncrop=1)
        
        ! Verify surfacewater arrays allocated
        call assert_true(allocated(state%surfwater%wlsbak), 'wlsbak allocated')
        call assert_equal_int(4, size(state%surfwater%wlsbak), 'size(wlsbak) = 4')
        
        ! Verify initial scalar values
        call assert_equal_real(0.0d0, state%surfwater%wlp, 1.0d-12, 'wlp = 0 initially')
        call assert_equal_real(0.0d0, state%surfwater%wls, 1.0d-12, 'wls = 0 initially')
        call assert_equal_real(0.0d0, state%surfwater%swst, 1.0d-12, 'swst = 0 initially')
        call assert_equal_real(0.0d0, state%surfwater%qdrd, 1.0d-12, 'qdrd = 0 initially')
        
    end subroutine test_surfacewater_initialization
    
    subroutine test_surfacewater_values()
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Surface Water Values ---')
        
        call swap_state_init(state, numnod=40, numlay=5, nrlevs=3, ncrop=1)
        
        ! Set water levels
        state%surfwater%wlp = -50.0d0
        state%surfwater%wls = -45.0d0
        state%surfwater%wlsold = -48.0d0
        state%surfwater%wlstar = -40.0d0
        
        call assert_equal_real(-50.0d0, state%surfwater%wlp, 1.0d-12, 'wlp = -50')
        call assert_equal_real(-45.0d0, state%surfwater%wls, 1.0d-12, 'wls = -45')
        call assert_equal_real(-48.0d0, state%surfwater%wlsold, 1.0d-12, 'wlsold = -48')
        call assert_equal_real(-40.0d0, state%surfwater%wlstar, 1.0d-12, 'wlstar = -40')
        
        ! Set storage and fluxes
        state%surfwater%swst = 100.0d0
        state%surfwater%qdrd = 0.5d0
        state%surfwater%cqdrd = 50.0d0
        state%surfwater%cwsupp = 25.0d0
        state%surfwater%cwout = 30.0d0
        
        call assert_equal_real(100.0d0, state%surfwater%swst, 1.0d-12, 'swst = 100')
        call assert_equal_real(0.5d0, state%surfwater%qdrd, 1.0d-12, 'qdrd = 0.5')
        call assert_equal_real(50.0d0, state%surfwater%cqdrd, 1.0d-12, 'cqdrd = 50')
        call assert_equal_real(25.0d0, state%surfwater%cwsupp, 1.0d-12, 'cwsupp = 25')
        call assert_equal_real(30.0d0, state%surfwater%cwout, 1.0d-12, 'cwout = 30')
        
        ! Test flags
        state%surfwater%overfl = .true.
        state%surfwater%fldecdt = .true.
        call assert_true(state%surfwater%overfl, 'overfl = true')
        call assert_true(state%surfwater%fldecdt, 'fldecdt = true')
        
        ! Test water level history for oscillation check
        state%surfwater%wlsbak(1) = -50.0d0
        state%surfwater%wlsbak(2) = -48.0d0
        state%surfwater%wlsbak(3) = -46.0d0
        state%surfwater%wlsbak(4) = -45.0d0
        
        call assert_equal_real(-50.0d0, state%surfwater%wlsbak(1), 1.0d-12, 'wlsbak(1) = -50')
        call assert_equal_real(-45.0d0, state%surfwater%wlsbak(4), 1.0d-12, 'wlsbak(4) = -45')
        
    end subroutine test_surfacewater_values
    
    subroutine test_multiple_drainage_instances()
        type(swap_state_t) :: state1, state2, state3
        
        call log_info('test', '--- Test: Multiple Drainage Instances ---')
        
        ! Initialize three states with different drainage levels
        call swap_state_init(state1, numnod=20, numlay=3, nrlevs=1, ncrop=1)
        call swap_state_init(state2, numnod=40, numlay=5, nrlevs=3, ncrop=1)
        call swap_state_init(state3, numnod=60, numlay=7, nrlevs=5, ncrop=2)
        
        ! Set different values in each instance
        state1%drain%qdrain(1) = 0.1d0
        state1%drain%qdrtot = 0.1d0
        
        state2%drain%qdrain(1) = 0.2d0
        state2%drain%qdrain(2) = 0.3d0
        state2%drain%qdrain(3) = 0.4d0
        state2%drain%qdrtot = 0.9d0
        
        state3%drain%qdrain(1) = 1.0d0
        state3%drain%qdrtot = 5.0d0
        
        ! Verify independence
        call assert_equal_real(0.1d0, state1%drain%qdrain(1), 1.0d-12, 'state1 qdrain(1) independent')
        call assert_equal_real(0.2d0, state2%drain%qdrain(1), 1.0d-12, 'state2 qdrain(1) independent')
        call assert_equal_real(1.0d0, state3%drain%qdrain(1), 1.0d-12, 'state3 qdrain(1) independent')
        
        call assert_equal_real(0.1d0, state1%drain%qdrtot, 1.0d-12, 'state1 qdrtot independent')
        call assert_equal_real(0.9d0, state2%drain%qdrtot, 1.0d-12, 'state2 qdrtot independent')
        call assert_equal_real(5.0d0, state3%drain%qdrtot, 1.0d-12, 'state3 qdrtot independent')
        
        ! Verify different array sizes
        call assert_equal_int(1, size(state1%drain%qdrain), 'state1 qdrain size 1')
        call assert_equal_int(3, size(state2%drain%qdrain), 'state2 qdrain size 3')
        call assert_equal_int(5, size(state3%drain%qdrain), 'state3 qdrain size 5')
        
        call assert_equal_int(1, state1%drain%nrlevs, 'state1 nrlevs = 1')
        call assert_equal_int(3, state2%drain%nrlevs, 'state2 nrlevs = 3')
        call assert_equal_int(5, state3%drain%nrlevs, 'state3 nrlevs = 5')
        
        call log_info('test', 'Multiple drainage instances are fully independent!')
        
    end subroutine test_multiple_drainage_instances

end program test_drainage_state
