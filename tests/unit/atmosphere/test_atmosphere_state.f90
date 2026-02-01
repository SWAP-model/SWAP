! ==============================================================================
! SWAP Atmosphere State Unit Tests
! ==============================================================================
! Tests for atmosphere_state_t type to ensure:
!   - Correct initialization of all fields
!   - Proper default values
!   - Independence of multiple instances
!   - Synchronization with module variables
!
! Author: SWAP Development Team
! Date: 2026-02-01
! ==============================================================================

program test_atmosphere_state
    use swap_state_mod
    use swap_log
    implicit none
    
    integer :: passed, failed, total
    
    passed = 0
    failed = 0
    total = 0
    
    call log_init(LOGLEVEL_DEBUG)
    
    call log_info('test', '========================================')
    call log_info('test', 'SWAP Atmosphere State Tests')
    call log_info('test', '========================================')
    
    ! Test groups
    call test_atmosphere_initialization(passed, failed, total)
    call test_atmosphere_meteo_values(passed, failed, total)
    call test_atmosphere_precipitation(passed, failed, total)
    call test_atmosphere_evapotranspiration(passed, failed, total)
    call test_atmosphere_cumulative(passed, failed, total)
    call test_atmosphere_etsine_state(passed, failed, total)
    call test_atmosphere_cn_method_state(passed, failed, total)
    call test_atmosphere_flags(passed, failed, total)
    call test_multiple_instances(passed, failed, total)
    
    ! Summary
    call log_info('test', '========================================')
    call log_info('test', 'Test Summary')
    call log_info('test', '========================================')
    write(*, '(A, I4)') 'Total tests: ', total
    write(*, '(A, I4)') 'Passed: ', passed
    write(*, '(A, I4)') 'Failed: ', failed
    
    if (failed == 0) then
        call log_info('test', 'ALL TESTS PASSED!')
        stop 0
    else
        call log_error('test', 'SOME TESTS FAILED!')
        stop 1
    end if
    
contains

    subroutine assert_true(condition, msg, passed, failed, total)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: msg
        integer, intent(inout) :: passed, failed, total
        
        total = total + 1
        if (condition) then
            passed = passed + 1
            call log_debug('assert', 'PASS: ' // trim(msg))
        else
            failed = failed + 1
            call log_error('assert', 'FAIL: ' // trim(msg))
        end if
    end subroutine assert_true
    
    subroutine assert_eq_real(actual, expected, msg, passed, failed, total, tol)
        real(8), intent(in) :: actual, expected
        character(len=*), intent(in) :: msg
        integer, intent(inout) :: passed, failed, total
        real(8), intent(in), optional :: tol
        real(8) :: tolerance
        
        tolerance = 1.0d-10
        if (present(tol)) tolerance = tol
        
        total = total + 1
        if (abs(actual - expected) < tolerance) then
            passed = passed + 1
            call log_debug('assert', 'PASS: ' // trim(msg))
        else
            failed = failed + 1
            call log_error('assert', 'FAIL: ' // trim(msg))
        end if
    end subroutine assert_eq_real
    
    subroutine assert_eq_int(actual, expected, msg, passed, failed, total)
        integer, intent(in) :: actual, expected
        character(len=*), intent(in) :: msg
        integer, intent(inout) :: passed, failed, total
        
        total = total + 1
        if (actual == expected) then
            passed = passed + 1
            call log_debug('assert', 'PASS: ' // trim(msg))
        else
            failed = failed + 1
            call log_error('assert', 'FAIL: ' // trim(msg))
        end if
    end subroutine assert_eq_int
    
    subroutine test_atmosphere_initialization(passed, failed, total)
        integer, intent(inout) :: passed, failed, total
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Atmosphere State Initialization ---')
        
        call swap_state_init(state, 40, 5, 2, 1)
        
        ! Test default values for meteorological state
        call assert_eq_real(state%atm%tav, 0.0d0, 'tav = 0 initially', passed, failed, total)
        call assert_eq_real(state%atm%tavd, 0.0d0, 'tavd = 0 initially', passed, failed, total)
        call assert_eq_real(state%atm%tmn, 0.0d0, 'tmn = 0 initially', passed, failed, total)
        call assert_eq_real(state%atm%tmx, 0.0d0, 'tmx = 0 initially', passed, failed, total)
        call assert_eq_real(state%atm%tmnr, 0.0d0, 'tmnr = 0 initially', passed, failed, total)
        call assert_eq_real(state%atm%rad, 0.0d0, 'rad = 0 initially', passed, failed, total)
        call assert_eq_real(state%atm%rh, 0.0d0, 'rh = 0 initially', passed, failed, total)
        call assert_eq_real(state%atm%lat, 0.0d0, 'lat = 0 initially', passed, failed, total)
        call assert_eq_real(state%atm%alt, 0.0d0, 'alt = 0 initially', passed, failed, total)
        call assert_eq_real(state%atm%daylp, 0.0d0, 'daylp = 0 initially', passed, failed, total)
        
        call swap_state_finalize(state)
    end subroutine test_atmosphere_initialization
    
    subroutine test_atmosphere_meteo_values(passed, failed, total)
        integer, intent(inout) :: passed, failed, total
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Atmosphere Meteorological Values ---')
        
        call swap_state_init(state, 40, 5, 2, 1)
        
        ! Set meteorological values
        state%atm%tav = 15.5d0
        state%atm%tavd = 18.2d0
        state%atm%tmn = 10.0d0
        state%atm%tmx = 22.0d0
        state%atm%tmnr = 9.5d0
        state%atm%rad = 20000000.0d0  ! J/m2/d
        state%atm%rh = 0.75d0
        state%atm%lat = 52.0d0
        state%atm%alt = 10.0d0
        state%atm%daylp = 14.5d0
        
        call assert_eq_real(state%atm%tav, 15.5d0, 'tav = 15.5', passed, failed, total)
        call assert_eq_real(state%atm%tavd, 18.2d0, 'tavd = 18.2', passed, failed, total)
        call assert_eq_real(state%atm%tmn, 10.0d0, 'tmn = 10.0', passed, failed, total)
        call assert_eq_real(state%atm%tmx, 22.0d0, 'tmx = 22.0', passed, failed, total)
        call assert_eq_real(state%atm%tmnr, 9.5d0, 'tmnr = 9.5', passed, failed, total)
        call assert_eq_real(state%atm%rad, 20000000.0d0, 'rad = 2e7', passed, failed, total)
        call assert_eq_real(state%atm%rh, 0.75d0, 'rh = 0.75', passed, failed, total)
        call assert_eq_real(state%atm%lat, 52.0d0, 'lat = 52.0', passed, failed, total)
        call assert_eq_real(state%atm%alt, 10.0d0, 'alt = 10.0', passed, failed, total)
        call assert_eq_real(state%atm%daylp, 14.5d0, 'daylp = 14.5', passed, failed, total)
        
        call swap_state_finalize(state)
    end subroutine test_atmosphere_meteo_values
    
    subroutine test_atmosphere_precipitation(passed, failed, total)
        integer, intent(inout) :: passed, failed, total
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Atmosphere Precipitation ---')
        
        call swap_state_init(state, 40, 5, 2, 1)
        
        ! Set precipitation values
        state%atm%grai = 1.5d0       ! 1.5 cm/d gross rain
        state%atm%graidt = 0.02d0    ! 0.02 cm in timestep
        state%atm%nraida = 1.2d0     ! 1.2 cm/d net rain
        state%atm%nraidt = 0.016d0   ! 0.016 cm in timestep
        state%atm%finterception = 0.8d0  ! 80% reaches soil
        
        call assert_eq_real(state%atm%grai, 1.5d0, 'grai = 1.5', passed, failed, total)
        call assert_eq_real(state%atm%graidt, 0.02d0, 'graidt = 0.02', passed, failed, total)
        call assert_eq_real(state%atm%nraida, 1.2d0, 'nraida = 1.2', passed, failed, total)
        call assert_eq_real(state%atm%nraidt, 0.016d0, 'nraidt = 0.016', passed, failed, total)
        call assert_eq_real(state%atm%finterception, 0.8d0, 'finterception = 0.8', passed, failed, total)
        
        ! Additional evaporation parameters
        state%atm%empreva = 0.1d0
        state%atm%fprecnosnow = 0.95d0
        call assert_eq_real(state%atm%empreva, 0.1d0, 'empreva = 0.1', passed, failed, total)
        call assert_eq_real(state%atm%fprecnosnow, 0.95d0, 'fprecnosnow = 0.95', passed, failed, total)
        
        call swap_state_finalize(state)
    end subroutine test_atmosphere_precipitation
    
    subroutine test_atmosphere_evapotranspiration(passed, failed, total)
        integer, intent(inout) :: passed, failed, total
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Atmosphere Evapotranspiration ---')
        
        call swap_state_init(state, 40, 5, 2, 1)
        
        ! Set evapotranspiration values
        state%atm%peva = 0.3d0       ! Potential E
        state%atm%pevaday = 0.3d0    ! Daily pot E
        state%atm%ptra = 0.25d0      ! Potential T
        state%atm%ptraday = 0.25d0   ! Daily pot T
        state%atm%tra = 0.22d0       ! Actual T
        state%atm%reva = 0.28d0      ! Actual E
        state%atm%atmdem = 0.5d0     ! Atmospheric demand
        state%atm%es0 = 0.35d0       ! Pot E wet bare soil
        state%atm%et0 = 0.4d0        ! Pot T dry crop
        state%atm%ew0 = 0.45d0       ! Pot T wet crop
        
        call assert_eq_real(state%atm%peva, 0.3d0, 'peva = 0.3', passed, failed, total)
        call assert_eq_real(state%atm%pevaday, 0.3d0, 'pevaday = 0.3', passed, failed, total)
        call assert_eq_real(state%atm%ptra, 0.25d0, 'ptra = 0.25', passed, failed, total)
        call assert_eq_real(state%atm%ptraday, 0.25d0, 'ptraday = 0.25', passed, failed, total)
        call assert_eq_real(state%atm%tra, 0.22d0, 'tra = 0.22', passed, failed, total)
        call assert_eq_real(state%atm%reva, 0.28d0, 'reva = 0.28', passed, failed, total)
        call assert_eq_real(state%atm%atmdem, 0.5d0, 'atmdem = 0.5', passed, failed, total)
        call assert_eq_real(state%atm%es0, 0.35d0, 'es0 = 0.35', passed, failed, total)
        call assert_eq_real(state%atm%et0, 0.4d0, 'et0 = 0.4', passed, failed, total)
        call assert_eq_real(state%atm%ew0, 0.45d0, 'ew0 = 0.45', passed, failed, total)
        
        call swap_state_finalize(state)
    end subroutine test_atmosphere_evapotranspiration
    
    subroutine test_atmosphere_cumulative(passed, failed, total)
        integer, intent(inout) :: passed, failed, total
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Atmosphere Cumulative Values ---')
        
        call swap_state_init(state, 40, 5, 2, 1)
        
        ! Set cumulative values
        state%atm%cgrai = 150.0d0    ! Cumulative gross precip
        state%atm%cnrai = 120.0d0    ! Cumulative net precip
        state%atm%caintc = 30.0d0    ! Cumulative interception
        state%atm%cevap = 50.0d0     ! Cumulative actual E
        state%atm%cpeva = 60.0d0     ! Cumulative potential E
        state%atm%cptra = 80.0d0     ! Cumulative potential T
        
        call assert_eq_real(state%atm%cgrai, 150.0d0, 'cgrai = 150', passed, failed, total)
        call assert_eq_real(state%atm%cnrai, 120.0d0, 'cnrai = 120', passed, failed, total)
        call assert_eq_real(state%atm%caintc, 30.0d0, 'caintc = 30', passed, failed, total)
        call assert_eq_real(state%atm%cevap, 50.0d0, 'cevap = 50', passed, failed, total)
        call assert_eq_real(state%atm%cpeva, 60.0d0, 'cpeva = 60', passed, failed, total)
        call assert_eq_real(state%atm%cptra, 80.0d0, 'cptra = 80', passed, failed, total)
        
        ! Intermediate values
        state%atm%inrai = 10.0d0
        state%atm%igrai = 12.0d0
        state%atm%ievap = 5.0d0
        state%atm%ipeva = 6.0d0
        state%atm%iptra = 8.0d0
        state%atm%ies0 = 3.5d0
        state%atm%iet0 = 4.0d0
        state%atm%iew0 = 4.5d0
        
        call assert_eq_real(state%atm%inrai, 10.0d0, 'inrai = 10', passed, failed, total)
        call assert_eq_real(state%atm%igrai, 12.0d0, 'igrai = 12', passed, failed, total)
        call assert_eq_real(state%atm%ievap, 5.0d0, 'ievap = 5', passed, failed, total)
        call assert_eq_real(state%atm%ipeva, 6.0d0, 'ipeva = 6', passed, failed, total)
        call assert_eq_real(state%atm%iptra, 8.0d0, 'iptra = 8', passed, failed, total)
        call assert_eq_real(state%atm%ies0, 3.5d0, 'ies0 = 3.5', passed, failed, total)
        call assert_eq_real(state%atm%iet0, 4.0d0, 'iet0 = 4.0', passed, failed, total)
        call assert_eq_real(state%atm%iew0, 4.5d0, 'iew0 = 4.5', passed, failed, total)
        
        call swap_state_finalize(state)
    end subroutine test_atmosphere_cumulative
    
    subroutine test_atmosphere_etsine_state(passed, failed, total)
        integer, intent(inout) :: passed, failed, total
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Atmosphere ETSine State (from SAVE) ---')
        
        call swap_state_init(state, 40, 5, 2, 1)
        
        ! Test ETSine state variables (previously local SAVE in meteodt.f90)
        call assert_eq_real(state%atm%tsunrise, 0.0d0, 'tsunrise = 0 initially', passed, failed, total)
        call assert_eq_real(state%atm%tsunset, 0.0d0, 'tsunset = 0 initially', passed, failed, total)
        
        ! Set values
        state%atm%tsunrise = 0.25d0   ! 6:00 AM
        state%atm%tsunset = 0.75d0    ! 6:00 PM
        
        call assert_eq_real(state%atm%tsunrise, 0.25d0, 'tsunrise = 0.25 (6AM)', passed, failed, total)
        call assert_eq_real(state%atm%tsunset, 0.75d0, 'tsunset = 0.75 (6PM)', passed, failed, total)
        
        ! Verify typical day length relationship
        call assert_true(state%atm%tsunset > state%atm%tsunrise, &
                        'tsunset > tsunrise (day length)', passed, failed, total)
        
        call swap_state_finalize(state)
    end subroutine test_atmosphere_etsine_state
    
    subroutine test_atmosphere_cn_method_state(passed, failed, total)
        integer, intent(inout) :: passed, failed, total
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Atmosphere CN Method State (from SAVE) ---')
        
        call swap_state_init(state, 40, 5, 2, 1)
        
        ! Test CN runoff method state variables (previously local SAVE in meteoday.f90)
        call assert_eq_int(state%atm%nod10_cn, 0, 'nod10_cn = 0 initially', passed, failed, total)
        call assert_eq_int(state%atm%icn_atm, 0, 'icn_atm = 0 initially', passed, failed, total)
        call assert_eq_real(state%atm%z10_cn, 0.0d0, 'z10_cn = 0 initially', passed, failed, total)
        
        ! Set values (simulate CN method initialization)
        state%atm%nod10_cn = 5      ! Node at -10cm
        state%atm%icn_atm = 1       ! Current position in CN table
        state%atm%z10_cn = 10.0d0   ! Depth = 10 cm
        
        call assert_eq_int(state%atm%nod10_cn, 5, 'nod10_cn = 5', passed, failed, total)
        call assert_eq_int(state%atm%icn_atm, 1, 'icn_atm = 1', passed, failed, total)
        call assert_eq_real(state%atm%z10_cn, 10.0d0, 'z10_cn = 10.0', passed, failed, total)
        
        call swap_state_finalize(state)
    end subroutine test_atmosphere_cn_method_state
    
    subroutine test_atmosphere_flags(passed, failed, total)
        integer, intent(inout) :: passed, failed, total
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Atmosphere Flags ---')
        
        call swap_state_init(state, 40, 5, 2, 1)
        
        ! Test default flag values
        call assert_true(.not. state%atm%fletsine, 'fletsine = false initially', passed, failed, total)
        call assert_true(.not. state%atm%flmeteodt, 'flmeteodt = false initially', passed, failed, total)
        call assert_true(.not. state%atm%flmetdetail, 'flmetdetail = false initially', passed, failed, total)
        call assert_true(.not. state%atm%flrainintens, 'flrainintens = false initially', passed, failed, total)
        call assert_true(.not. state%atm%flupdmetdet, 'flupdmetdet = false initially', passed, failed, total)
        
        ! Set flags
        state%atm%fletsine = .true.
        state%atm%flmeteodt = .true.
        state%atm%flmetdetail = .true.
        state%atm%flrainintens = .true.
        state%atm%flupdmetdet = .true.
        
        call assert_true(state%atm%fletsine, 'fletsine = true after set', passed, failed, total)
        call assert_true(state%atm%flmeteodt, 'flmeteodt = true after set', passed, failed, total)
        call assert_true(state%atm%flmetdetail, 'flmetdetail = true after set', passed, failed, total)
        call assert_true(state%atm%flrainintens, 'flrainintens = true after set', passed, failed, total)
        call assert_true(state%atm%flupdmetdet, 'flupdmetdet = true after set', passed, failed, total)
        
        call swap_state_finalize(state)
    end subroutine test_atmosphere_flags
    
    subroutine test_multiple_instances(passed, failed, total)
        integer, intent(inout) :: passed, failed, total
        type(swap_state_t) :: state1, state2, state3
        
        call log_info('test', '--- Test: Multiple Atmosphere Instances ---')
        
        call swap_state_init(state1, 20, 3, 1, 1)
        call swap_state_init(state2, 40, 5, 2, 1)
        call swap_state_init(state3, 60, 7, 3, 2)
        
        ! Set different values in each instance
        state1%atm%tav = 10.0d0
        state2%atm%tav = 20.0d0
        state3%atm%tav = 30.0d0
        
        state1%atm%grai = 0.5d0
        state2%atm%grai = 1.0d0
        state3%atm%grai = 1.5d0
        
        state1%atm%tsunrise = 0.2d0
        state2%atm%tsunrise = 0.25d0
        state3%atm%tsunrise = 0.3d0
        
        state1%atm%nod10_cn = 3
        state2%atm%nod10_cn = 5
        state3%atm%nod10_cn = 7
        
        ! Verify independence - meteorological values
        call assert_eq_real(state1%atm%tav, 10.0d0, 'state1 tav independent', passed, failed, total)
        call assert_eq_real(state2%atm%tav, 20.0d0, 'state2 tav independent', passed, failed, total)
        call assert_eq_real(state3%atm%tav, 30.0d0, 'state3 tav independent', passed, failed, total)
        
        ! Verify independence - precipitation
        call assert_eq_real(state1%atm%grai, 0.5d0, 'state1 grai independent', passed, failed, total)
        call assert_eq_real(state2%atm%grai, 1.0d0, 'state2 grai independent', passed, failed, total)
        call assert_eq_real(state3%atm%grai, 1.5d0, 'state3 grai independent', passed, failed, total)
        
        ! Verify independence - ETSine state
        call assert_eq_real(state1%atm%tsunrise, 0.2d0, 'state1 tsunrise independent', passed, failed, total)
        call assert_eq_real(state2%atm%tsunrise, 0.25d0, 'state2 tsunrise independent', passed, failed, total)
        call assert_eq_real(state3%atm%tsunrise, 0.3d0, 'state3 tsunrise independent', passed, failed, total)
        
        ! Verify independence - CN method state
        call assert_eq_int(state1%atm%nod10_cn, 3, 'state1 nod10_cn independent', passed, failed, total)
        call assert_eq_int(state2%atm%nod10_cn, 5, 'state2 nod10_cn independent', passed, failed, total)
        call assert_eq_int(state3%atm%nod10_cn, 7, 'state3 nod10_cn independent', passed, failed, total)
        
        ! Modify state1 and verify no effect on others
        state1%atm%tav = 99.0d0
        state1%atm%tsunrise = 0.99d0
        state1%atm%nod10_cn = 99
        
        call assert_eq_real(state2%atm%tav, 20.0d0, 'state2 tav still independent', passed, failed, total)
        call assert_eq_real(state3%atm%tav, 30.0d0, 'state3 tav still independent', passed, failed, total)
        call assert_eq_real(state2%atm%tsunrise, 0.25d0, 'state2 tsunrise still independent', passed, failed, total)
        call assert_eq_real(state3%atm%tsunrise, 0.3d0, 'state3 tsunrise still independent', passed, failed, total)
        call assert_eq_int(state2%atm%nod10_cn, 5, 'state2 nod10_cn still independent', passed, failed, total)
        call assert_eq_int(state3%atm%nod10_cn, 7, 'state3 nod10_cn still independent', passed, failed, total)
        
        call log_info('test', 'Multiple atmosphere instances are fully independent!')
        
        call swap_state_finalize(state1)
        call swap_state_finalize(state2)
        call swap_state_finalize(state3)
    end subroutine test_multiple_instances

end program test_atmosphere_state
