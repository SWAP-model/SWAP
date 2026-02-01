!> @file test_heat_state.f90
!> @brief Unit tests for heat state structures
!> 
!> Tests the heat_state_t type from swap_state_mod.
!> This test is self-contained and does not require Variables module arrays.
!> Key functionality tested:
!> - Heat state initialization
!> - All heat arrays (tsoil, heacap, heacon, rfcp, etc.)
!> - Thermal properties and frost parameters
!> - Multiple heat instances independence

program test_heat_state
    use swap_state_mod
    use swap_log
    implicit none
    
    integer :: tests_run, tests_passed
    
    tests_run = 0
    tests_passed = 0
    
    ! Initialize logging
    call log_init(LOGLEVEL_DEBUG)
    call log_info('test', '========================================')
    call log_info('test', 'SWAP Heat State Tests')
    call log_info('test', '========================================')
    
    ! Run test suites
    call test_heat_initialization()
    call test_heat_array_values()
    call test_heat_temperatures()
    call test_heat_thermal_properties()
    call test_heat_soil_composition()
    call test_heat_boundary_conditions()
    call test_heat_frost_state()
    call test_multiple_heat_instances()
    
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

    subroutine test_heat_initialization()
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Heat State Initialization ---')
        
        ! Initialize state with 40 nodes and 5 layers
        call swap_state_init(state, numnod=40, numlay=5, nrlevs=2, ncrop=1)
        
        ! Verify heat arrays are allocated
        call assert_true(allocated(state%heat%tsoil), 'heat%tsoil allocated')
        call assert_true(allocated(state%heat%heacap), 'heat%heacap allocated')
        call assert_true(allocated(state%heat%heacon), 'heat%heacon allocated')
        call assert_true(allocated(state%heat%rfcp), 'heat%rfcp allocated')
        call assert_true(allocated(state%heat%fclay), 'heat%fclay allocated')
        call assert_true(allocated(state%heat%forg), 'heat%forg allocated')
        call assert_true(allocated(state%heat%fquartz), 'heat%fquartz allocated')
        call assert_true(allocated(state%heat%pclay), 'heat%pclay allocated')
        call assert_true(allocated(state%heat%psand), 'heat%psand allocated')
        call assert_true(allocated(state%heat%psilt), 'heat%psilt allocated')
        call assert_true(allocated(state%heat%orgmat), 'heat%orgmat allocated')
        call assert_true(allocated(state%heat%tembtab), 'heat%tembtab allocated')
        call assert_true(allocated(state%heat%temtoptab), 'heat%temtoptab allocated')
        call assert_true(allocated(state%heat%zh), 'heat%zh allocated')
        
        ! Verify array sizes
        call assert_equal_int(40, size(state%heat%tsoil), 'size(heat%tsoil) = 40')
        call assert_equal_int(40, size(state%heat%heacap), 'size(heat%heacap) = 40')
        call assert_equal_int(40, size(state%heat%heacon), 'size(heat%heacon) = 40')
        call assert_equal_int(40, size(state%heat%rfcp), 'size(heat%rfcp) = 40')
        call assert_equal_int(40, size(state%heat%fclay), 'size(heat%fclay) = 40')
        call assert_equal_int(40, size(state%heat%forg), 'size(heat%forg) = 40')
        call assert_equal_int(40, size(state%heat%fquartz), 'size(heat%fquartz) = 40')
        call assert_equal_int(5, size(state%heat%pclay), 'size(heat%pclay) = 5')
        call assert_equal_int(5, size(state%heat%psand), 'size(heat%psand) = 5')
        call assert_equal_int(5, size(state%heat%psilt), 'size(heat%psilt) = 5')
        call assert_equal_int(5, size(state%heat%orgmat), 'size(heat%orgmat) = 5')
        
        call swap_state_finalize(state)
    end subroutine

    subroutine test_heat_array_values()
        type(swap_state_t) :: state
        integer :: i
        
        call log_info('test', '--- Test: Heat Array Values ---')
        
        call swap_state_init(state, numnod=40, numlay=5, nrlevs=2, ncrop=1)
        
        ! Verify initial values (tsoil defaults to 10.0, rfcp to 1.0)
        call assert_equal_real(10.0d0, state%heat%tsoil(1), 1.0d-10, 'tsoil(1) = 10 initially (default)')
        call assert_equal_real(0.0d0, state%heat%heacap(1), 1.0d-10, 'heacap(1) = 0 initially')
        call assert_equal_real(0.0d0, state%heat%heacon(1), 1.0d-10, 'heacon(1) = 0 initially')
        
        ! Set temperature values
        do i = 1, 40
            state%heat%tsoil(i) = 10.0d0 + 0.1d0 * i
        end do
        
        ! Verify values were set correctly
        call assert_equal_real(10.1d0, state%heat%tsoil(1), 1.0d-10, 'tsoil(1) = 10.1')
        call assert_equal_real(14.0d0, state%heat%tsoil(40), 1.0d-10, 'tsoil(40) = 14.0')
        
        call swap_state_finalize(state)
    end subroutine

    subroutine test_heat_temperatures()
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Heat Temperatures ---')
        
        call swap_state_init(state, numnod=40, numlay=5, nrlevs=2, ncrop=1)
        
        ! Set boundary temperatures
        state%heat%tetop = 15.0d0
        state%heat%tebot = 10.0d0
        
        call assert_equal_real(15.0d0, state%heat%tetop, 1.0d-10, 'tetop = 15')
        call assert_equal_real(10.0d0, state%heat%tebot, 1.0d-10, 'tebot = 10')
        
        ! Set soil temperatures
        state%heat%tsoil(1) = 14.0d0
        state%heat%tsoil(20) = 12.0d0
        state%heat%tsoil(40) = 10.5d0
        
        call assert_equal_real(14.0d0, state%heat%tsoil(1), 1.0d-10, 'tsoil(1) = 14')
        call assert_equal_real(12.0d0, state%heat%tsoil(20), 1.0d-10, 'tsoil(20) = 12')
        call assert_equal_real(10.5d0, state%heat%tsoil(40), 1.0d-10, 'tsoil(40) = 10.5')
        
        call swap_state_finalize(state)
    end subroutine

    subroutine test_heat_thermal_properties()
        type(swap_state_t) :: state
        integer :: i
        
        call log_info('test', '--- Test: Heat Thermal Properties ---')
        
        call swap_state_init(state, numnod=40, numlay=5, nrlevs=2, ncrop=1)
        
        ! Set thermal properties per compartment
        do i = 1, 40
            state%heat%heacap(i) = 2.5d6  ! J/m3/K
            state%heat%heacon(i) = 1.0d0  ! W/m/K
        end do
        
        call assert_equal_real(2.5d6, state%heat%heacap(1), 1.0d-10, 'heacap(1) = 2.5e6')
        call assert_equal_real(2.5d6, state%heat%heacap(40), 1.0d-10, 'heacap(40) = 2.5e6')
        call assert_equal_real(1.0d0, state%heat%heacon(1), 1.0d-10, 'heacon(1) = 1.0')
        call assert_equal_real(1.0d0, state%heat%heacon(40), 1.0d-10, 'heacon(40) = 1.0')
        
        call swap_state_finalize(state)
    end subroutine

    subroutine test_heat_soil_composition()
        type(swap_state_t) :: state
        integer :: i
        
        call log_info('test', '--- Test: Heat Soil Composition ---')
        
        call swap_state_init(state, numnod=40, numlay=5, nrlevs=2, ncrop=1)
        
        ! Set volume fractions per compartment
        do i = 1, 40
            state%heat%fclay(i) = 0.15d0
            state%heat%forg(i) = 0.05d0
            state%heat%fquartz(i) = 0.30d0
        end do
        
        call assert_equal_real(0.15d0, state%heat%fclay(1), 1.0d-10, 'fclay(1) = 0.15')
        call assert_equal_real(0.05d0, state%heat%forg(1), 1.0d-10, 'forg(1) = 0.05')
        call assert_equal_real(0.30d0, state%heat%fquartz(1), 1.0d-10, 'fquartz(1) = 0.30')
        
        ! Set mass fractions per layer
        do i = 1, 5
            state%heat%pclay(i) = 0.20d0
            state%heat%psand(i) = 0.50d0
            state%heat%psilt(i) = 0.30d0
            state%heat%orgmat(i) = 0.03d0
        end do
        
        call assert_equal_real(0.20d0, state%heat%pclay(1), 1.0d-10, 'pclay(1) = 0.20')
        call assert_equal_real(0.50d0, state%heat%psand(1), 1.0d-10, 'psand(1) = 0.50')
        call assert_equal_real(0.30d0, state%heat%psilt(1), 1.0d-10, 'psilt(1) = 0.30')
        call assert_equal_real(0.03d0, state%heat%orgmat(1), 1.0d-10, 'orgmat(1) = 0.03')
        
        ! Verify layer-based arrays have correct size
        call assert_equal_int(5, size(state%heat%pclay), 'size(pclay) = 5')
        call assert_equal_int(5, size(state%heat%psand), 'size(psand) = 5')
        call assert_equal_int(5, size(state%heat%psilt), 'size(psilt) = 5')
        call assert_equal_int(5, size(state%heat%orgmat), 'size(orgmat) = 5')
        
        call swap_state_finalize(state)
    end subroutine

    subroutine test_heat_boundary_conditions()
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Heat Boundary Conditions ---')
        
        call swap_state_init(state, numnod=40, numlay=5, nrlevs=2, ncrop=1)
        
        ! Set analytical solution parameters
        state%heat%swcalt = 1
        state%heat%tmean = 10.0d0
        state%heat%tampli = 8.0d0
        state%heat%timref = 91.0d0
        state%heat%ddamp = 100.0d0
        
        call assert_equal_int(1, state%heat%swcalt, 'swcalt = 1')
        call assert_equal_real(10.0d0, state%heat%tmean, 1.0d-10, 'tmean = 10')
        call assert_equal_real(8.0d0, state%heat%tampli, 1.0d-10, 'tampli = 8')
        call assert_equal_real(91.0d0, state%heat%timref, 1.0d-10, 'timref = 91')
        call assert_equal_real(100.0d0, state%heat%ddamp, 1.0d-10, 'ddamp = 100')
        
        ! Set numerical solution switches
        state%heat%swcalt = 2
        state%heat%swtopbhea = 2
        state%heat%swbotbhea = 1
        
        call assert_equal_int(2, state%heat%swcalt, 'swcalt = 2')
        call assert_equal_int(2, state%heat%swtopbhea, 'swtopbhea = 2')
        call assert_equal_int(1, state%heat%swbotbhea, 'swbotbhea = 1')
        
        call swap_state_finalize(state)
    end subroutine

    subroutine test_heat_frost_state()
        type(swap_state_t) :: state
        integer :: i
        
        call log_info('test', '--- Test: Heat Frost State ---')
        
        call swap_state_init(state, numnod=40, numlay=5, nrlevs=2, ncrop=1)
        
        ! Set frost switch
        state%heat%swfrost = 1
        call assert_equal_int(1, state%heat%swfrost, 'swfrost = 1')
        
        ! Set frost reduction factors per compartment
        do i = 1, 40
            state%heat%rfcp(i) = 1.0d0  ! No reduction initially
        end do
        
        ! Simulate frost in upper compartments
        state%heat%rfcp(1) = 0.1d0
        state%heat%rfcp(2) = 0.3d0
        state%heat%rfcp(3) = 0.5d0
        state%heat%rfcp(4) = 0.8d0
        state%heat%rfcp(5) = 1.0d0
        
        call assert_equal_real(0.1d0, state%heat%rfcp(1), 1.0d-10, 'rfcp(1) = 0.1 (frozen)')
        call assert_equal_real(0.3d0, state%heat%rfcp(2), 1.0d-10, 'rfcp(2) = 0.3')
        call assert_equal_real(0.5d0, state%heat%rfcp(3), 1.0d-10, 'rfcp(3) = 0.5')
        call assert_equal_real(1.0d0, state%heat%rfcp(5), 1.0d-10, 'rfcp(5) = 1.0 (unfrozen)')
        
        ! Set frost depth tracking
        state%heat%zfrosttop = 0.0d0
        state%heat%zfrostbot = -15.0d0
        state%heat%tfroststa = 335.0d0  ! Frost start day
        state%heat%tfrostend = 60.0d0   ! Frost end day
        state%heat%nodfrostbot = 5
        
        call assert_equal_real(0.0d0, state%heat%zfrosttop, 1.0d-10, 'zfrosttop = 0')
        call assert_equal_real(-15.0d0, state%heat%zfrostbot, 1.0d-10, 'zfrostbot = -15')
        call assert_equal_real(335.0d0, state%heat%tfroststa, 1.0d-10, 'tfroststa = 335')
        call assert_equal_real(60.0d0, state%heat%tfrostend, 1.0d-10, 'tfrostend = 60')
        call assert_equal_int(5, state%heat%nodfrostbot, 'nodfrostbot = 5')
        
        call swap_state_finalize(state)
    end subroutine

    subroutine test_multiple_heat_instances()
        type(swap_state_t) :: state1, state2, state3
        
        call log_info('test', '--- Test: Multiple Heat Instances ---')
        
        ! Create three independent states with different sizes
        call swap_state_init(state1, numnod=20, numlay=3, nrlevs=1, ncrop=1)
        call swap_state_init(state2, numnod=40, numlay=5, nrlevs=2, ncrop=1)
        call swap_state_init(state3, numnod=60, numlay=7, nrlevs=3, ncrop=2)
        
        ! Set different temperature values in each
        state1%heat%tsoil(1) = 5.0d0
        state1%heat%tetop = 8.0d0
        state1%heat%tmean = 7.0d0
        
        state2%heat%tsoil(1) = 10.0d0
        state2%heat%tetop = 15.0d0
        state2%heat%tmean = 12.0d0
        
        state3%heat%tsoil(1) = 20.0d0
        state3%heat%tetop = 25.0d0
        state3%heat%tmean = 18.0d0
        
        ! Verify independence
        call assert_equal_real(5.0d0, state1%heat%tsoil(1), 1.0d-10, 'state1 tsoil(1) independent')
        call assert_equal_real(10.0d0, state2%heat%tsoil(1), 1.0d-10, 'state2 tsoil(1) independent')
        call assert_equal_real(20.0d0, state3%heat%tsoil(1), 1.0d-10, 'state3 tsoil(1) independent')
        
        call assert_equal_real(8.0d0, state1%heat%tetop, 1.0d-10, 'state1 tetop independent')
        call assert_equal_real(15.0d0, state2%heat%tetop, 1.0d-10, 'state2 tetop independent')
        call assert_equal_real(25.0d0, state3%heat%tetop, 1.0d-10, 'state3 tetop independent')
        
        call assert_equal_real(7.0d0, state1%heat%tmean, 1.0d-10, 'state1 tmean independent')
        call assert_equal_real(12.0d0, state2%heat%tmean, 1.0d-10, 'state2 tmean independent')
        call assert_equal_real(18.0d0, state3%heat%tmean, 1.0d-10, 'state3 tmean independent')
        
        ! Verify array sizes are different
        call assert_equal_int(20, size(state1%heat%tsoil), 'state1 tsoil size 20')
        call assert_equal_int(40, size(state2%heat%tsoil), 'state2 tsoil size 40')
        call assert_equal_int(60, size(state3%heat%tsoil), 'state3 tsoil size 60')
        
        call assert_equal_int(20, size(state1%heat%rfcp), 'state1 rfcp size 20')
        call assert_equal_int(40, size(state2%heat%rfcp), 'state2 rfcp size 40')
        call assert_equal_int(60, size(state3%heat%rfcp), 'state3 rfcp size 60')
        
        call assert_equal_int(3, size(state1%heat%pclay), 'state1 pclay size 3')
        call assert_equal_int(5, size(state2%heat%pclay), 'state2 pclay size 5')
        call assert_equal_int(7, size(state3%heat%pclay), 'state3 pclay size 7')
        
        call log_info('test', 'Multiple heat instances are fully independent!')
        
        call swap_state_finalize(state1)
        call swap_state_finalize(state2)
        call swap_state_finalize(state3)
    end subroutine

end program test_heat_state
