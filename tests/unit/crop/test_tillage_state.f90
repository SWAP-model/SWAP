!> @file test_tillage_state.f90
!> @brief Unit tests for tillage state structures
!> 
!> Tests the tillage_state_t type from swap_state_mod.
!> This test is self-contained and does not require Variables module arrays.
!> Key functionality tested:
!> - Tillage state initialization
!> - All tillage arrays (density tracking, event tables, etc.)
!> - Multiple tillage state instances independence

program test_tillage_state
    use swap_state_mod
    use swap_log
    implicit none
    
    integer :: tests_run, tests_passed
    
    tests_run = 0
    tests_passed = 0
    
    ! Initialize logging
    call log_init(LOGLEVEL_DEBUG)
    call log_info('test', '========================================')
    call log_info('test', 'SWAP Tillage State Tests')
    call log_info('test', '========================================')
    
    ! Run test suites
    call test_tillage_initialization()
    call test_tillage_event_arrays()
    call test_tillage_type_arrays()
    call test_tillage_layer_arrays()
    call test_tillage_scalars()
    call test_multiple_tillage_instances()
    
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

    subroutine test_tillage_initialization()
        type(tillage_state_t) :: tstate
        integer :: numlay, maxtill, maxtypes
        
        call log_info('test', '--- Test: Tillage State Initialization ---')
        
        numlay = 10
        maxtill = 5
        maxtypes = 3
        
        ! Initialize tillage state
        call tillage_state_init(tstate, numlay, maxtill, maxtypes)
        
        ! Verify per-event arrays are allocated
        call assert_true(allocated(tstate%Date_tillage), 'Date_tillage allocated')
        call assert_true(allocated(tstate%Z_tillage), 'Z_tillage allocated')
        call assert_true(allocated(tstate%I_tillage), 'I_tillage allocated')
        call assert_true(allocated(tstate%Type_Tillage), 'Type_Tillage allocated')
        
        ! Verify sizes (Date_tillage has maxtill+1)
        call assert_equal_int(maxtill + 1, size(tstate%Date_tillage), 'Date_tillage size')
        call assert_equal_int(maxtill, size(tstate%Z_tillage), 'Z_tillage size')
        call assert_equal_int(maxtill, size(tstate%I_tillage), 'I_tillage size')
        call assert_equal_int(maxtill, size(tstate%Type_Tillage), 'Type_Tillage size')
        
        ! Verify per-type arrays are allocated
        call assert_true(allocated(tstate%iType_Tillage), 'iType_Tillage allocated')
        call assert_true(allocated(tstate%iTT1), 'iTT1 allocated')
        call assert_true(allocated(tstate%iTT2), 'iTT2 allocated')
        call assert_true(allocated(tstate%TAB_Rho_tillage), 'TAB_Rho_tillage allocated')
        call assert_true(allocated(tstate%TAB_Rho_cons), 'TAB_Rho_cons allocated')
        call assert_true(allocated(tstate%TAB_K_R_cons), 'TAB_K_R_cons allocated')
        call assert_true(allocated(tstate%TAB_Rho_match), 'TAB_Rho_match allocated')
        call assert_true(allocated(tstate%TAB_N_match), 'TAB_N_match allocated')
        
        ! Verify per-layer arrays are allocated
        call assert_true(allocated(tstate%Rho_tillage), 'Rho_tillage allocated')
        call assert_true(allocated(tstate%Rho_cons), 'Rho_cons allocated')
        call assert_true(allocated(tstate%Rho_last), 'Rho_last allocated')
        call assert_true(allocated(tstate%K_R_cons), 'K_R_cons allocated')
        call assert_true(allocated(tstate%Rho_match), 'Rho_match allocated')
        call assert_true(allocated(tstate%N_match), 'N_match allocated')
        call assert_true(allocated(tstate%Slope_match), 'Slope_match allocated')
        
        ! Verify per-layer array sizes
        call assert_equal_int(numlay, size(tstate%Rho_tillage), 'Rho_tillage size')
        call assert_equal_int(numlay, size(tstate%Rho_cons), 'Rho_cons size')
        call assert_equal_int(numlay, size(tstate%Rho_last), 'Rho_last size')
        call assert_equal_int(numlay, size(tstate%K_R_cons), 'K_R_cons size')
        call assert_equal_int(numlay, size(tstate%Rho_match), 'Rho_match size')
        call assert_equal_int(numlay, size(tstate%N_match), 'N_match size')
        call assert_equal_int(numlay, size(tstate%Slope_match), 'Slope_match size')
        
        ! Clean up
        call tillage_state_finalize(tstate)
        
    end subroutine test_tillage_initialization
    
    subroutine test_tillage_event_arrays()
        type(tillage_state_t) :: tstate
        integer :: i
        
        call log_info('test', '--- Test: Tillage Event Arrays ---')
        
        call tillage_state_init(tstate, 5, 3, 2)
        
        ! Set event values
        tstate%Date_tillage(1) = 45000.0d0  ! Some date in t1900 format
        tstate%Date_tillage(2) = 45100.0d0
        tstate%Date_tillage(3) = 45200.0d0
        tstate%Z_tillage(1) = 25.0d0        ! 25 cm tillage depth
        tstate%Z_tillage(2) = 30.0d0
        tstate%Z_tillage(3) = 25.0d0
        tstate%I_tillage(1) = 0.8d0         ! 80% intensity
        tstate%I_tillage(2) = 1.0d0         ! 100% intensity
        tstate%I_tillage(3) = 0.5d0         ! 50% intensity
        tstate%Type_Tillage(1) = 1
        tstate%Type_Tillage(2) = 2
        tstate%Type_Tillage(3) = 1
        
        ! Verify values
        call assert_equal_real(45000.0d0, tstate%Date_tillage(1), 1.0d-6, 'Date_tillage(1)')
        call assert_equal_real(45100.0d0, tstate%Date_tillage(2), 1.0d-6, 'Date_tillage(2)')
        call assert_equal_real(25.0d0, tstate%Z_tillage(1), 1.0d-6, 'Z_tillage(1)')
        call assert_equal_real(30.0d0, tstate%Z_tillage(2), 1.0d-6, 'Z_tillage(2)')
        call assert_equal_real(0.8d0, tstate%I_tillage(1), 1.0d-6, 'I_tillage(1)')
        call assert_equal_real(1.0d0, tstate%I_tillage(2), 1.0d-6, 'I_tillage(2)')
        call assert_equal_int(1, tstate%Type_Tillage(1), 'Type_Tillage(1)')
        call assert_equal_int(2, tstate%Type_Tillage(2), 'Type_Tillage(2)')
        
        call tillage_state_finalize(tstate)
        
    end subroutine test_tillage_event_arrays
    
    subroutine test_tillage_type_arrays()
        type(tillage_state_t) :: tstate
        
        call log_info('test', '--- Test: Tillage Type Arrays ---')
        
        call tillage_state_init(tstate, 5, 3, 4)
        
        ! Set type table values
        tstate%TAB_Rho_tillage(1) = 1200.0d0  ! kg/m3
        tstate%TAB_Rho_tillage(2) = 1150.0d0
        tstate%TAB_Rho_cons(1) = 1400.0d0
        tstate%TAB_Rho_cons(2) = 1350.0d0
        tstate%TAB_K_R_cons(1) = 0.01d0
        tstate%TAB_K_R_cons(2) = 0.02d0
        tstate%TAB_Rho_match(1) = 1300.0d0
        tstate%TAB_N_match(1) = 1.2d0
        
        tstate%iType_Tillage(1) = 1
        tstate%iType_Tillage(2) = 1
        tstate%iType_Tillage(3) = 2
        tstate%iType_Tillage(4) = 2
        tstate%iTT1(1) = 1
        tstate%iTT2(1) = 2
        tstate%iTT1(2) = 3
        tstate%iTT2(2) = 4
        
        ! Verify values
        call assert_equal_real(1200.0d0, tstate%TAB_Rho_tillage(1), 1.0d-6, 'TAB_Rho_tillage(1)')
        call assert_equal_real(1400.0d0, tstate%TAB_Rho_cons(1), 1.0d-6, 'TAB_Rho_cons(1)')
        call assert_equal_real(0.01d0, tstate%TAB_K_R_cons(1), 1.0d-9, 'TAB_K_R_cons(1)')
        call assert_equal_real(1300.0d0, tstate%TAB_Rho_match(1), 1.0d-6, 'TAB_Rho_match(1)')
        call assert_equal_real(1.2d0, tstate%TAB_N_match(1), 1.0d-6, 'TAB_N_match(1)')
        call assert_equal_int(1, tstate%iTT1(1), 'iTT1(1)')
        call assert_equal_int(2, tstate%iTT2(1), 'iTT2(1)')
        
        call tillage_state_finalize(tstate)
        
    end subroutine test_tillage_type_arrays
    
    subroutine test_tillage_layer_arrays()
        type(tillage_state_t) :: tstate
        integer :: i
        
        call log_info('test', '--- Test: Tillage Layer Arrays ---')
        
        call tillage_state_init(tstate, 8, 3, 2)
        
        ! Set per-layer density values
        do i = 1, 8
            tstate%Rho_tillage(i) = 1200.0d0 + i * 25.0d0
            tstate%Rho_cons(i) = 1400.0d0 + i * 20.0d0
            tstate%Rho_last(i) = 1300.0d0 + i * 22.0d0
            tstate%K_R_cons(i) = 0.01d0 + i * 0.002d0
            tstate%Rho_match(i) = 1350.0d0
            tstate%N_match(i) = 1.15d0 + i * 0.05d0
            tstate%Slope_match(i) = 0.001d0 * i
        end do
        
        ! Verify some values
        call assert_equal_real(1225.0d0, tstate%Rho_tillage(1), 1.0d-6, 'Rho_tillage(1)')
        call assert_equal_real(1300.0d0, tstate%Rho_tillage(4), 1.0d-6, 'Rho_tillage(4)')
        call assert_equal_real(1420.0d0, tstate%Rho_cons(1), 1.0d-6, 'Rho_cons(1)')
        call assert_equal_real(1322.0d0, tstate%Rho_last(1), 1.0d-6, 'Rho_last(1)')
        call assert_equal_real(0.012d0, tstate%K_R_cons(1), 1.0d-9, 'K_R_cons(1)')
        call assert_equal_real(1.20d0, tstate%N_match(1), 1.0d-6, 'N_match(1)')
        call assert_equal_real(0.001d0, tstate%Slope_match(1), 1.0d-9, 'Slope_match(1)')
        
        call tillage_state_finalize(tstate)
        
    end subroutine test_tillage_layer_arrays
    
    subroutine test_tillage_scalars()
        type(tillage_state_t) :: tstate
        
        call log_info('test', '--- Test: Tillage Scalar Values ---')
        
        call tillage_state_init(tstate, 5, 3, 2)
        
        ! Test default values
        call assert_equal_int(0, tstate%swtill, 'swtill default')
        call assert_equal_int(2, tstate%i_n_model, 'i_n_model default')
        call assert_equal_int(2, tstate%iRedist, 'iRedist default')
        call assert_equal_int(0, tstate%Ntill, 'Ntill default')
        call assert_equal_int(1, tstate%iTill, 'iTill default')
        call assert_equal_int(0, tstate%Ntypes, 'Ntypes default')
        call assert_equal_real(0.0d0, tstate%Max_Z_tillage, 1.0d-10, 'Max_Z_tillage default')
        call assert_equal_real(0.0d0, tstate%sumDWC, 1.0d-10, 'sumDWC default')
        call assert_equal_real(0.0d0, tstate%sumAvail1, 1.0d-10, 'sumAvail1 default')
        call assert_equal_real(0.0d0, tstate%sumAvail2, 1.0d-10, 'sumAvail2 default')
        
        ! Set values and verify
        tstate%swtill = 1
        tstate%Ntill = 5
        tstate%iTill = 2
        tstate%Ntypes = 3
        tstate%MaxNumSoilHo = 4
        tstate%MaxNumSoilCP = 20
        tstate%Max_Z_tillage = 30.0d0
        tstate%i_n_model = 3
        tstate%iRedist = 1
        tstate%sumDWC = 0.5d0
        tstate%sumAvail1 = 10.0d0
        tstate%sumAvail2 = 8.0d0
        
        call assert_equal_int(1, tstate%swtill, 'swtill set')
        call assert_equal_int(5, tstate%Ntill, 'Ntill set')
        call assert_equal_int(2, tstate%iTill, 'iTill set')
        call assert_equal_int(3, tstate%Ntypes, 'Ntypes set')
        call assert_equal_int(4, tstate%MaxNumSoilHo, 'MaxNumSoilHo set')
        call assert_equal_int(20, tstate%MaxNumSoilCP, 'MaxNumSoilCP set')
        call assert_equal_real(30.0d0, tstate%Max_Z_tillage, 1.0d-10, 'Max_Z_tillage set')
        call assert_equal_int(3, tstate%i_n_model, 'i_n_model set')
        call assert_equal_int(1, tstate%iRedist, 'iRedist set')
        call assert_equal_real(0.5d0, tstate%sumDWC, 1.0d-10, 'sumDWC set')
        call assert_equal_real(10.0d0, tstate%sumAvail1, 1.0d-10, 'sumAvail1 set')
        call assert_equal_real(8.0d0, tstate%sumAvail2, 1.0d-10, 'sumAvail2 set')
        
        call tillage_state_finalize(tstate)
        
    end subroutine test_tillage_scalars
    
    subroutine test_multiple_tillage_instances()
        type(tillage_state_t) :: tstate1, tstate2
        
        call log_info('test', '--- Test: Multiple Tillage State Instances ---')
        
        ! Initialize two separate instances
        call tillage_state_init(tstate1, 5, 4, 2)
        call tillage_state_init(tstate2, 8, 6, 3)
        
        ! Set different values
        tstate1%swtill = 1
        tstate1%Ntill = 4
        tstate1%Max_Z_tillage = 25.0d0
        tstate1%Rho_tillage(1) = 1200.0d0
        tstate1%Rho_cons(1) = 1400.0d0
        
        tstate2%swtill = 1
        tstate2%Ntill = 6
        tstate2%Max_Z_tillage = 35.0d0
        tstate2%Rho_tillage(1) = 1100.0d0
        tstate2%Rho_cons(1) = 1350.0d0
        
        ! Verify instances are independent
        call assert_equal_int(4, tstate1%Ntill, 'instance1 Ntill')
        call assert_equal_int(6, tstate2%Ntill, 'instance2 Ntill')
        call assert_true(tstate1%Ntill /= tstate2%Ntill, 'Ntill independent')
        
        call assert_equal_real(25.0d0, tstate1%Max_Z_tillage, 1.0d-6, 'instance1 Max_Z_tillage')
        call assert_equal_real(35.0d0, tstate2%Max_Z_tillage, 1.0d-6, 'instance2 Max_Z_tillage')
        call assert_true(abs(tstate1%Max_Z_tillage - tstate2%Max_Z_tillage) > 1.0d0, 'Max_Z_tillage independent')
        
        call assert_equal_real(1200.0d0, tstate1%Rho_tillage(1), 1.0d-6, 'instance1 Rho_tillage(1)')
        call assert_equal_real(1100.0d0, tstate2%Rho_tillage(1), 1.0d-6, 'instance2 Rho_tillage(1)')
        call assert_true(abs(tstate1%Rho_tillage(1) - tstate2%Rho_tillage(1)) > 50.0d0, 'Rho_tillage independent')
        
        ! Verify different array sizes
        call assert_equal_int(5, size(tstate1%Rho_tillage), 'instance1 Rho_tillage size')
        call assert_equal_int(8, size(tstate2%Rho_tillage), 'instance2 Rho_tillage size')
        call assert_equal_int(4, size(tstate1%Z_tillage), 'instance1 Z_tillage size')
        call assert_equal_int(6, size(tstate2%Z_tillage), 'instance2 Z_tillage size')
        
        ! Clean up
        call tillage_state_finalize(tstate1)
        call tillage_state_finalize(tstate2)
        
        ! Verify arrays deallocated
        call assert_true(.not. allocated(tstate1%Rho_tillage), 'instance1 Rho_tillage deallocated')
        call assert_true(.not. allocated(tstate2%Rho_tillage), 'instance2 Rho_tillage deallocated')
        call assert_true(.not. allocated(tstate1%Date_tillage), 'instance1 Date_tillage deallocated')
        call assert_true(.not. allocated(tstate2%Date_tillage), 'instance2 Date_tillage deallocated')
        
    end subroutine test_multiple_tillage_instances

end program test_tillage_state
