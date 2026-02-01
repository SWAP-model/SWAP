!> @file test_solute_state.f90
!> @brief Unit tests for solute state structures
!> 
!> Tests the solute_state_t type from swap_state_mod.
!> This test is self-contained and does not require Variables module arrays.
!> Key functionality tested:
!> - Solute state initialization
!> - All solute arrays (cml, cmsy, ldis, kf, decpot, etc.)
!> - Age tracer variables
!> - Multiple solute instances independence

program test_solute_state
    use swap_state_mod
    use swap_log
    implicit none
    
    integer :: tests_run, tests_passed
    
    tests_run = 0
    tests_passed = 0
    
    ! Initialize logging
    call log_init(LOGLEVEL_DEBUG)
    call log_info('test', '========================================')
    call log_info('test', 'SWAP Solute State Tests')
    call log_info('test', '========================================')
    
    ! Run test suites
    call test_solute_initialization()
    call test_solute_array_values()
    call test_solute_concentrations()
    call test_solute_transport_params()
    call test_solute_age_tracer()
    call test_multiple_solute_instances()
    
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

    subroutine test_solute_initialization()
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Solute State Initialization ---')
        
        ! Initialize state with 40 nodes and 5 layers
        call swap_state_init(state, numnod=40, numlay=5, nrlevs=2, ncrop=1)
        
        ! Verify solute arrays are allocated
        call assert_true(allocated(state%solute%cml), 'solute%cml allocated')
        call assert_true(allocated(state%solute%cmsy), 'solute%cmsy allocated')
        call assert_true(allocated(state%solute%ldis), 'solute%ldis allocated')
        call assert_true(allocated(state%solute%kf), 'solute%kf allocated')
        call assert_true(allocated(state%solute%decpot), 'solute%decpot allocated')
        call assert_true(allocated(state%solute%fdepth), 'solute%fdepth allocated')
        call assert_true(allocated(state%solute%icAgeDra), 'solute%icAgeDra allocated')
        call assert_true(allocated(state%solute%cseeptab), 'solute%cseeptab allocated')
        call assert_true(allocated(state%solute%zc), 'solute%zc allocated')
        
        ! Verify array sizes
        call assert_equal_int(40, size(state%solute%cml), 'size(solute%cml) = 40')
        call assert_equal_int(40, size(state%solute%cmsy), 'size(solute%cmsy) = 40')
        call assert_equal_int(5, size(state%solute%ldis), 'size(solute%ldis) = 5')
        call assert_equal_int(5, size(state%solute%kf), 'size(solute%kf) = 5')
        call assert_equal_int(5, size(state%solute%decpot), 'size(solute%decpot) = 5')
        call assert_equal_int(5, size(state%solute%fdepth), 'size(solute%fdepth) = 5')
        call assert_equal_int(2, size(state%solute%icAgeDra), 'size(solute%icAgeDra) = nrlevs')
        
        call swap_state_finalize(state)
    end subroutine

    subroutine test_solute_array_values()
        type(swap_state_t) :: state
        integer :: i
        
        call log_info('test', '--- Test: Solute Array Values ---')
        
        call swap_state_init(state, numnod=40, numlay=5, nrlevs=2, ncrop=1)
        
        ! Verify initial values are zero
        call assert_equal_real(0.0d0, state%solute%cml(1), 1.0d-10, 'cml(1) = 0 initially')
        call assert_equal_real(0.0d0, state%solute%sampro, 1.0d-10, 'sampro = 0 initially')
        call assert_equal_real(0.0d0, state%solute%samini, 1.0d-10, 'samini = 0 initially')
        call assert_equal_real(0.0d0, state%solute%dectot, 1.0d-10, 'dectot = 0 initially')
        
        ! Set some concentration values
        do i = 1, 40
            state%solute%cml(i) = 0.1d0 * i
            state%solute%cmsy(i) = 0.05d0 * i
        end do
        
        ! Verify values were set correctly
        call assert_equal_real(0.1d0, state%solute%cml(1), 1.0d-10, 'cml(1) = 0.1')
        call assert_equal_real(4.0d0, state%solute%cml(40), 1.0d-10, 'cml(40) = 4.0')
        call assert_equal_real(0.05d0, state%solute%cmsy(1), 1.0d-10, 'cmsy(1) = 0.05')
        call assert_equal_real(2.0d0, state%solute%cmsy(40), 1.0d-10, 'cmsy(40) = 2.0')
        
        ! Set layer-based values
        do i = 1, 5
            state%solute%ldis(i) = 5.0d0 * i
            state%solute%kf(i) = 0.1d0 * i
            state%solute%decpot(i) = 0.01d0 * i
        end do
        
        call assert_equal_real(5.0d0, state%solute%ldis(1), 1.0d-10, 'ldis(1) = 5.0')
        call assert_equal_real(25.0d0, state%solute%ldis(5), 1.0d-10, 'ldis(5) = 25.0')
        call assert_equal_real(0.1d0, state%solute%kf(1), 1.0d-10, 'kf(1) = 0.1')
        call assert_equal_real(0.05d0, state%solute%decpot(5), 1.0d-10, 'decpot(5) = 0.05')
        
        call swap_state_finalize(state)
    end subroutine

    subroutine test_solute_concentrations()
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Solute Concentrations ---')
        
        call swap_state_init(state, numnod=40, numlay=5, nrlevs=2, ncrop=1)
        
        ! Set boundary concentrations
        state%solute%cpond = 0.5d0
        state%solute%csurf = 1.0d0
        state%solute%cdrain = 0.3d0
        state%solute%cseep = 0.2d0
        state%solute%cpre = 0.01d0
        state%solute%cirr = 0.05d0
        
        call assert_equal_real(0.5d0, state%solute%cpond, 1.0d-10, 'cpond = 0.5')
        call assert_equal_real(1.0d0, state%solute%csurf, 1.0d-10, 'csurf = 1.0')
        call assert_equal_real(0.3d0, state%solute%cdrain, 1.0d-10, 'cdrain = 0.3')
        call assert_equal_real(0.2d0, state%solute%cseep, 1.0d-10, 'cseep = 0.2')
        call assert_equal_real(0.01d0, state%solute%cpre, 1.0d-10, 'cpre = 0.01')
        call assert_equal_real(0.05d0, state%solute%cirr, 1.0d-10, 'cirr = 0.05')
        
        ! Set cumulative amounts
        state%solute%sampro = 100.0d0
        state%solute%sqbot = 10.0d0
        state%solute%sqdra = 5.0d0
        state%solute%dectot = 2.0d0
        state%solute%rottot = 3.0d0
        
        call assert_equal_real(100.0d0, state%solute%sampro, 1.0d-10, 'sampro = 100')
        call assert_equal_real(10.0d0, state%solute%sqbot, 1.0d-10, 'sqbot = 10')
        call assert_equal_real(5.0d0, state%solute%sqdra, 1.0d-10, 'sqdra = 5')
        call assert_equal_real(2.0d0, state%solute%dectot, 1.0d-10, 'dectot = 2')
        call assert_equal_real(3.0d0, state%solute%rottot, 1.0d-10, 'rottot = 3')
        
        call swap_state_finalize(state)
    end subroutine

    subroutine test_solute_transport_params()
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Solute Transport Parameters ---')
        
        call swap_state_init(state, numnod=40, numlay=5, nrlevs=2, ncrop=1)
        
        ! Verify default values
        call assert_equal_real(1.0d0, state%solute%frexp, 1.0d-10, 'frexp default = 1.0')
        call assert_equal_real(1.0d0, state%solute%tscf, 1.0d-10, 'tscf default = 1.0')
        call assert_equal_real(0.7d0, state%solute%bexp, 1.0d-10, 'bexp default = 0.7')
        call assert_equal_real(0.01d0, state%solute%rtheta, 1.0d-10, 'rtheta default = 0.01')
        call assert_equal_real(1.0d0, state%solute%cref, 1.0d-10, 'cref default = 1.0')
        
        ! Set transport parameters
        state%solute%ddif = 1.2d-4
        state%solute%frexp = 0.8d0
        state%solute%tscf = 0.5d0
        state%solute%gampar = 0.08d0
        state%solute%bexp = 0.6d0
        
        call assert_equal_real(1.2d-4, state%solute%ddif, 1.0d-10, 'ddif = 1.2e-4')
        call assert_equal_real(0.8d0, state%solute%frexp, 1.0d-10, 'frexp = 0.8')
        call assert_equal_real(0.5d0, state%solute%tscf, 1.0d-10, 'tscf = 0.5')
        call assert_equal_real(0.08d0, state%solute%gampar, 1.0d-10, 'gampar = 0.08')
        call assert_equal_real(0.6d0, state%solute%bexp, 1.0d-10, 'bexp = 0.6')
        
        ! Set salt stress parameters
        state%solute%salthead = -400.0d0
        state%solute%saltmax = 3.0d0
        state%solute%saltslope = 0.1d0
        
        call assert_equal_real(-400.0d0, state%solute%salthead, 1.0d-10, 'salthead = -400')
        call assert_equal_real(3.0d0, state%solute%saltmax, 1.0d-10, 'saltmax = 3')
        call assert_equal_real(0.1d0, state%solute%saltslope, 1.0d-10, 'saltslope = 0.1')
        
        call swap_state_finalize(state)
    end subroutine

    subroutine test_solute_age_tracer()
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Solute Age Tracer ---')
        
        call swap_state_init(state, numnod=40, numlay=5, nrlevs=3, ncrop=1)
        
        ! Set age tracer values
        state%solute%AgeGwl1m = 365.0d0
        state%solute%icAgeBot = 10.0d0
        state%solute%icAgeRot = 5.0d0
        state%solute%icAgeSur = 2.0d0
        
        call assert_equal_real(365.0d0, state%solute%AgeGwl1m, 1.0d-10, 'AgeGwl1m = 365')
        call assert_equal_real(10.0d0, state%solute%icAgeBot, 1.0d-10, 'icAgeBot = 10')
        call assert_equal_real(5.0d0, state%solute%icAgeRot, 1.0d-10, 'icAgeRot = 5')
        call assert_equal_real(2.0d0, state%solute%icAgeSur, 1.0d-10, 'icAgeSur = 2')
        
        ! Set age tracer per drainage level
        state%solute%icAgeDra(1) = 100.0d0
        state%solute%icAgeDra(2) = 200.0d0
        state%solute%icAgeDra(3) = 300.0d0
        
        call assert_equal_real(100.0d0, state%solute%icAgeDra(1), 1.0d-10, 'icAgeDra(1) = 100')
        call assert_equal_real(200.0d0, state%solute%icAgeDra(2), 1.0d-10, 'icAgeDra(2) = 200')
        call assert_equal_real(300.0d0, state%solute%icAgeDra(3), 1.0d-10, 'icAgeDra(3) = 300')
        
        ! Set age tracer boundary state (now in solute_state_t)
        state%solute%Ageirr = 0.0d0
        state%solute%Agedrain = 50.0d0
        state%solute%Agepre = 0.0d0
        state%solute%Agepond = 1.0d0
        state%solute%Agepondm1 = 0.5d0
        
        call assert_equal_real(0.0d0, state%solute%Ageirr, 1.0d-10, 'Ageirr = 0')
        call assert_equal_real(50.0d0, state%solute%Agedrain, 1.0d-10, 'Agedrain = 50')
        call assert_equal_real(0.0d0, state%solute%Agepre, 1.0d-10, 'Agepre = 0')
        call assert_equal_real(1.0d0, state%solute%Agepond, 1.0d-10, 'Agepond = 1')
        call assert_equal_real(0.5d0, state%solute%Agepondm1, 1.0d-10, 'Agepondm1 = 0.5')
        
        ! Set flux-related age tracer state
        state%solute%icAgetopupw = 15.0d0
        state%solute%icAgetopdwn = 20.0d0
        state%solute%ArMpSs = 0.05d0
        
        call assert_equal_real(15.0d0, state%solute%icAgetopupw, 1.0d-10, 'icAgetopupw = 15')
        call assert_equal_real(20.0d0, state%solute%icAgetopdwn, 1.0d-10, 'icAgetopdwn = 20')
        call assert_equal_real(0.05d0, state%solute%ArMpSs, 1.0d-10, 'ArMpSs = 0.05')
        
        call swap_state_finalize(state)
    end subroutine

    subroutine test_multiple_solute_instances()
        type(swap_state_t) :: state1, state2, state3
        
        call log_info('test', '--- Test: Multiple Solute Instances ---')
        
        ! Create three independent states with different sizes
        call swap_state_init(state1, numnod=20, numlay=3, nrlevs=1, ncrop=1)
        call swap_state_init(state2, numnod=40, numlay=5, nrlevs=2, ncrop=1)
        call swap_state_init(state3, numnod=60, numlay=7, nrlevs=3, ncrop=2)
        
        ! Set different values in each
        state1%solute%cml(1) = 1.0d0
        state1%solute%sampro = 100.0d0
        state1%solute%AgeGwl1m = 10.0d0
        
        state2%solute%cml(1) = 2.0d0
        state2%solute%sampro = 200.0d0
        state2%solute%AgeGwl1m = 20.0d0
        
        state3%solute%cml(1) = 3.0d0
        state3%solute%sampro = 300.0d0
        state3%solute%AgeGwl1m = 30.0d0
        
        ! Verify independence
        call assert_equal_real(1.0d0, state1%solute%cml(1), 1.0d-10, 'state1 cml(1) independent')
        call assert_equal_real(2.0d0, state2%solute%cml(1), 1.0d-10, 'state2 cml(1) independent')
        call assert_equal_real(3.0d0, state3%solute%cml(1), 1.0d-10, 'state3 cml(1) independent')
        
        call assert_equal_real(100.0d0, state1%solute%sampro, 1.0d-10, 'state1 sampro independent')
        call assert_equal_real(200.0d0, state2%solute%sampro, 1.0d-10, 'state2 sampro independent')
        call assert_equal_real(300.0d0, state3%solute%sampro, 1.0d-10, 'state3 sampro independent')
        
        call assert_equal_real(10.0d0, state1%solute%AgeGwl1m, 1.0d-10, 'state1 AgeGwl1m independent')
        call assert_equal_real(20.0d0, state2%solute%AgeGwl1m, 1.0d-10, 'state2 AgeGwl1m independent')
        call assert_equal_real(30.0d0, state3%solute%AgeGwl1m, 1.0d-10, 'state3 AgeGwl1m independent')
        
        ! Verify array sizes are different
        call assert_equal_int(20, size(state1%solute%cml), 'state1 cml size 20')
        call assert_equal_int(40, size(state2%solute%cml), 'state2 cml size 40')
        call assert_equal_int(60, size(state3%solute%cml), 'state3 cml size 60')
        
        call assert_equal_int(3, size(state1%solute%ldis), 'state1 ldis size 3')
        call assert_equal_int(5, size(state2%solute%ldis), 'state2 ldis size 5')
        call assert_equal_int(7, size(state3%solute%ldis), 'state3 ldis size 7')
        
        call log_info('test', 'Multiple solute instances are fully independent!')
        
        call swap_state_finalize(state1)
        call swap_state_finalize(state2)
        call swap_state_finalize(state3)
    end subroutine

end program test_solute_state
