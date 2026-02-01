!> @file test_soil_state.f90
!> @brief Unit tests for soil state structures
!> 
!> Tests the soil_state_t type from swap_state_mod.
!> This test is self-contained and does not require Variables module arrays.
!> Key functionality tested:
!> - Soil state initialization
!> - All soil arrays (h, theta, k, gwl, fluxes, etc.)
!> - Multiple soil instances independence
!> - Headcalc iteration tracking

program test_soil_state
    use swap_state_mod
    use swap_log
    implicit none
    
    integer :: tests_run, tests_passed
    
    tests_run = 0
    tests_passed = 0
    
    ! Initialize logging
    call log_init(LOGLEVEL_DEBUG)
    call log_info('test', '========================================')
    call log_info('test', 'SWAP Soil State Tests')
    call log_info('test', '========================================')
    
    ! Run test suites
    call test_soil_initialization()
    call test_soil_array_values()
    call test_soil_previous_timestep()
    call test_soil_groundwater_state()
    call test_soil_flux_arrays()
    call test_soil_cumulative_fluxes()
    call test_soil_storage()
    call test_soil_evaporation_reduction()
    call test_soil_headcalc_tracking()
    call test_multiple_soil_instances()
    
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

    subroutine test_soil_initialization()
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Soil State Initialization ---')
        
        ! Initialize state
        call swap_state_init(state, numnod=40, numlay=5, nrlevs=2, ncrop=1)
        
        ! Verify soil arrays are allocated
        call assert_true(allocated(state%soil%h), 'soil%h allocated')
        call assert_true(allocated(state%soil%theta), 'soil%theta allocated')
        call assert_true(allocated(state%soil%k), 'soil%k allocated')
        call assert_true(allocated(state%soil%z), 'soil%z allocated')
        call assert_true(allocated(state%soil%dz), 'soil%dz allocated')
        call assert_true(allocated(state%soil%hm1), 'soil%hm1 allocated')
        call assert_true(allocated(state%soil%thetm1), 'soil%thetm1 allocated')
        call assert_true(allocated(state%soil%q), 'soil%q allocated')
        call assert_true(allocated(state%soil%qrot), 'soil%qrot allocated')
        call assert_true(allocated(state%soil%kmean), 'soil%kmean allocated')
        
        ! Verify sizes
        call assert_equal_int(40, size(state%soil%h), 'size(soil%h) = 40')
        call assert_equal_int(40, size(state%soil%theta), 'size(soil%theta) = 40')
        call assert_equal_int(40, size(state%soil%hm1), 'size(soil%hm1) = 40')
        call assert_equal_int(41, size(state%soil%q), 'size(soil%q) = 41 (numnod+1)')
        call assert_equal_int(40, size(state%soil%qrot), 'size(soil%qrot) = 40')
        call assert_equal_int(41, size(state%soil%kmean), 'size(soil%kmean) = 41')
        
        ! Automatic cleanup at end of scope
    end subroutine

    subroutine test_soil_array_values()
        type(swap_state_t) :: state
        integer :: i
        
        call log_info('test', '--- Test: Soil Array Values ---')
        
        call swap_state_init(state, numnod=40, numlay=5, nrlevs=2, ncrop=1)
        
        ! Verify initial values are zero
        call assert_equal_real(0.0d0, state%soil%h(1), 1.0d-10, 'soil%h(1) = 0 initially')
        call assert_equal_real(0.0d0, state%soil%theta(1), 1.0d-10, 'soil%theta(1) = 0 initially')
        call assert_equal_real(0.0d0, state%soil%gwl, 1.0d-10, 'soil%gwl = 0 initially')
        call assert_equal_real(0.0d0, state%soil%pond, 1.0d-10, 'soil%pond = 0 initially')
        
        ! Set values and verify
        do i = 1, 40
            state%soil%h(i) = -real(i, 8) * 10.0d0
            state%soil%theta(i) = 0.30d0 + real(i, 8) * 0.001d0
        end do
        
        call assert_equal_real(-10.0d0, state%soil%h(1), 1.0d-10, 'h(1) = -10')
        call assert_equal_real(-400.0d0, state%soil%h(40), 1.0d-10, 'h(40) = -400')
        call assert_equal_real(0.301d0, state%soil%theta(1), 1.0d-10, 'theta(1) = 0.301')
        call assert_equal_real(0.340d0, state%soil%theta(40), 1.0d-10, 'theta(40) = 0.340')
        
        ! Automatic cleanup at end of scope
    end subroutine

    subroutine test_soil_previous_timestep()
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Soil Previous Timestep Values ---')
        
        call swap_state_init(state, numnod=40, numlay=5, nrlevs=2, ncrop=1)
        
        ! Set current and previous timestep values
        state%soil%h(1) = -100.0d0
        state%soil%hm1(1) = -90.0d0
        state%soil%theta(1) = 0.35d0
        state%soil%thetm1(1) = 0.36d0
        state%soil%gwl = -150.0d0
        state%soil%gwlm1 = -145.0d0
        state%soil%pond = 2.0d0
        state%soil%pondm1 = 1.8d0
        
        ! Verify values are independent
        call assert_true(state%soil%h(1) /= state%soil%hm1(1), 'h and hm1 independent')
        call assert_true(state%soil%theta(1) /= state%soil%thetm1(1), 'theta and thetm1 independent')
        call assert_true(state%soil%gwl /= state%soil%gwlm1, 'gwl and gwlm1 independent')
        call assert_true(state%soil%pond /= state%soil%pondm1, 'pond and pondm1 independent')
        
        ! Verify exact values
        call assert_equal_real(-100.0d0, state%soil%h(1), 1.0d-10, 'h(1) = -100')
        call assert_equal_real(-90.0d0, state%soil%hm1(1), 1.0d-10, 'hm1(1) = -90')
        
        ! Automatic cleanup at end of scope
    end subroutine

    subroutine test_soil_groundwater_state()
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Soil Groundwater State ---')
        
        call swap_state_init(state, numnod=40, numlay=5, nrlevs=2, ncrop=1)
        
        ! Set groundwater-related values
        state%soil%gwl = -200.0d0
        state%soil%gwlm1 = -195.0d0
        state%soil%gwli = -190.0d0
        state%soil%gwlinp = -185.0d0
        state%soil%nodgwl = 25
        state%soil%npegwl = 15
        state%soil%bpegwl = 10
        
        ! Verify all groundwater fields
        call assert_equal_real(-200.0d0, state%soil%gwl, 1.0d-10, 'gwl set correctly')
        call assert_equal_real(-195.0d0, state%soil%gwlm1, 1.0d-10, 'gwlm1 set correctly')
        call assert_equal_real(-190.0d0, state%soil%gwli, 1.0d-10, 'gwli set correctly')
        call assert_equal_real(-185.0d0, state%soil%gwlinp, 1.0d-10, 'gwlinp set correctly')
        call assert_equal_int(25, state%soil%nodgwl, 'nodgwl set correctly')
        call assert_equal_int(15, state%soil%npegwl, 'npegwl set correctly')
        call assert_equal_int(10, state%soil%bpegwl, 'bpegwl set correctly')
        
        ! Automatic cleanup at end of scope
    end subroutine

    subroutine test_soil_flux_arrays()
        type(swap_state_t) :: state
        integer :: i
        
        call log_info('test', '--- Test: Soil Flux Arrays ---')
        
        call swap_state_init(state, numnod=40, numlay=5, nrlevs=2, ncrop=1)
        
        ! Set flux values
        do i = 1, 40
            state%soil%q(i) = real(i, 8) * 0.01d0
            state%soil%qrot(i) = real(i, 8) * 0.005d0
        end do
        state%soil%q(41) = 0.41d0  ! Bottom boundary flux
        state%soil%qtop = -0.3d0
        state%soil%qbot = 0.2d0
        
        ! Verify flux array values
        call assert_equal_real(0.01d0, state%soil%q(1), 1.0d-10, 'q(1) = 0.01')
        call assert_equal_real(0.40d0, state%soil%q(40), 1.0d-10, 'q(40) = 0.40')
        call assert_equal_real(0.41d0, state%soil%q(41), 1.0d-10, 'q(41) = 0.41')
        call assert_equal_real(0.20d0, state%soil%qrot(40), 1.0d-10, 'qrot(40) = 0.20')
        call assert_equal_real(-0.3d0, state%soil%qtop, 1.0d-10, 'qtop = -0.3')
        call assert_equal_real(0.2d0, state%soil%qbot, 1.0d-10, 'qbot = 0.2')
        
        ! Automatic cleanup at end of scope
    end subroutine

    subroutine test_soil_cumulative_fluxes()
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Soil Cumulative Fluxes ---')
        
        call swap_state_init(state, numnod=40, numlay=5, nrlevs=2, ncrop=1)
        
        ! Set cumulative flux values
        state%soil%cqbot = 100.0d0
        state%soil%cqbotdo = 50.0d0
        state%soil%cqbotup = 50.0d0
        state%soil%cqtdo = 80.0d0
        state%soil%cqtup = 20.0d0
        state%soil%cqrot = 30.0d0
        state%soil%cqdra = 45.0d0
        state%soil%crunoff = 15.0d0
        state%soil%crunon = 5.0d0
        
        ! Verify cumulative values
        call assert_equal_real(100.0d0, state%soil%cqbot, 1.0d-10, 'cqbot = 100')
        call assert_equal_real(50.0d0, state%soil%cqbotdo, 1.0d-10, 'cqbotdo = 50')
        call assert_equal_real(50.0d0, state%soil%cqbotup, 1.0d-10, 'cqbotup = 50')
        call assert_equal_real(80.0d0, state%soil%cqtdo, 1.0d-10, 'cqtdo = 80')
        call assert_equal_real(20.0d0, state%soil%cqtup, 1.0d-10, 'cqtup = 20')
        call assert_equal_real(30.0d0, state%soil%cqrot, 1.0d-10, 'cqrot = 30')
        call assert_equal_real(45.0d0, state%soil%cqdra, 1.0d-10, 'cqdra = 45')
        call assert_equal_real(15.0d0, state%soil%crunoff, 1.0d-10, 'crunoff = 15')
        call assert_equal_real(5.0d0, state%soil%crunon, 1.0d-10, 'crunon = 5')
        
        ! Automatic cleanup at end of scope
    end subroutine

    subroutine test_soil_storage()
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Soil Storage Tracking ---')
        
        call swap_state_init(state, numnod=40, numlay=5, nrlevs=2, ncrop=1)
        
        ! Set storage values
        state%soil%volact = 150.0d0
        state%soil%volini = 145.0d0
        state%soil%volm1 = 148.0d0
        
        ! Verify storage values
        call assert_equal_real(150.0d0, state%soil%volact, 1.0d-10, 'volact = 150')
        call assert_equal_real(145.0d0, state%soil%volini, 1.0d-10, 'volini = 145')
        call assert_equal_real(148.0d0, state%soil%volm1, 1.0d-10, 'volm1 = 148')
        
        ! Automatic cleanup at end of scope
    end subroutine

    subroutine test_soil_evaporation_reduction()
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Soil Evaporation Reduction ---')
        
        call swap_state_init(state, numnod=40, numlay=5, nrlevs=2, ncrop=1)
        
        ! Set evaporation reduction values
        state%soil%saev = 5.0d0
        state%soil%spev = 3.0d0
        state%soil%ldwet = 2.5d0
        state%soil%cofred = 0.8d0
        
        ! Verify values
        call assert_equal_real(5.0d0, state%soil%saev, 1.0d-10, 'saev = 5')
        call assert_equal_real(3.0d0, state%soil%spev, 1.0d-10, 'spev = 3')
        call assert_equal_real(2.5d0, state%soil%ldwet, 1.0d-10, 'ldwet = 2.5')
        call assert_equal_real(0.8d0, state%soil%cofred, 1.0d-10, 'cofred = 0.8')
        
        ! Automatic cleanup at end of scope
    end subroutine

    subroutine test_soil_headcalc_tracking()
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Soil Headcalc Tracking ---')
        
        call swap_state_init(state, numnod=40, numlay=5, nrlevs=2, ncrop=1)
        
        ! Verify default values for headcalc tracking
        call assert_true(state%soil%flwarn_hc, 'flwarn_hc default true')
        call assert_equal_int(0, state%soil%iwarn_hc, 'iwarn_hc default 0')
        call assert_equal_int(0, state%soil%nstep_hc, 'nstep_hc default 0')
        
        ! Set headcalc tracking values
        state%soil%flwarn_hc = .false.
        state%soil%iwarn_hc = 5
        state%soil%nstep_hc = 100
        
        ! Verify values set correctly
        call assert_true(.not. state%soil%flwarn_hc, 'flwarn_hc set to false')
        call assert_equal_int(5, state%soil%iwarn_hc, 'iwarn_hc set to 5')
        call assert_equal_int(100, state%soil%nstep_hc, 'nstep_hc set to 100')
        
        ! Automatic cleanup at end of scope
    end subroutine

    subroutine test_multiple_soil_instances()
        type(swap_state_t) :: state1, state2, state3
        
        call log_info('test', '--- Test: Multiple Soil Instances ---')
        
        ! Initialize three instances with different sizes
        call swap_state_init(state1, numnod=20, numlay=3, nrlevs=1, ncrop=1)
        call swap_state_init(state2, numnod=40, numlay=5, nrlevs=2, ncrop=1)
        call swap_state_init(state3, numnod=80, numlay=10, nrlevs=3, ncrop=2)
        
        ! Set different values in each
        state1%soil%h(1) = -100.0d0
        state1%soil%gwl = -50.0d0
        state1%soil%volact = 50.0d0
        
        state2%soil%h(1) = -200.0d0
        state2%soil%gwl = -100.0d0
        state2%soil%volact = 100.0d0
        
        state3%soil%h(1) = -300.0d0
        state3%soil%gwl = -150.0d0
        state3%soil%volact = 200.0d0
        
        ! Verify they are independent
        call assert_equal_real(-100.0d0, state1%soil%h(1), 1.0d-10, 'state1 h independent')
        call assert_equal_real(-200.0d0, state2%soil%h(1), 1.0d-10, 'state2 h independent')
        call assert_equal_real(-300.0d0, state3%soil%h(1), 1.0d-10, 'state3 h independent')
        
        call assert_equal_real(-50.0d0, state1%soil%gwl, 1.0d-10, 'state1 gwl independent')
        call assert_equal_real(-100.0d0, state2%soil%gwl, 1.0d-10, 'state2 gwl independent')
        call assert_equal_real(-150.0d0, state3%soil%gwl, 1.0d-10, 'state3 gwl independent')
        
        call assert_equal_real(50.0d0, state1%soil%volact, 1.0d-10, 'state1 volact independent')
        call assert_equal_real(100.0d0, state2%soil%volact, 1.0d-10, 'state2 volact independent')
        call assert_equal_real(200.0d0, state3%soil%volact, 1.0d-10, 'state3 volact independent')
        
        ! Verify different array sizes
        call assert_equal_int(20, size(state1%soil%h), 'state1 h size 20')
        call assert_equal_int(40, size(state2%soil%h), 'state2 h size 40')
        call assert_equal_int(80, size(state3%soil%h), 'state3 h size 80')
        
        call log_info('test', 'Multiple soil instances are fully independent!')
        
        ! Automatic cleanup at end of scope
        ! Automatic cleanup at end of scope
        ! Automatic cleanup at end of scope
    end subroutine

end program test_soil_state
