! ==============================================================================
! SWAP Macropore State Unit Tests
! ==============================================================================
! Tests for macropore_state_t type to ensure:
!   - Correct initialization of all fields
!   - Proper allocation of arrays
!   - Proper default values
!   - Independence of multiple instances
!   - Work arrays from SAVE are properly handled
!
! Author: SWAP Development Team
! Date: 2026-02-01
! ==============================================================================

program test_macropore_state
    use swap_state_mod
    use swap_log
    implicit none
    
    integer :: passed, failed, total
    
    passed = 0
    failed = 0
    total = 0
    
    call log_init(LOGLEVEL_DEBUG)
    
    call log_info('test', '========================================')
    call log_info('test', 'SWAP Macropore State Tests')
    call log_info('test', '========================================')
    
    ! Test groups
    call test_macropore_initialization(passed, failed, total)
    call test_macropore_domain_storage(passed, failed, total)
    call test_macropore_fluxes(passed, failed, total)
    call test_macropore_cumulative(passed, failed, total)
    call test_macropore_incremental(passed, failed, total)
    call test_macropore_arrays(passed, failed, total)
    call test_macropore_work_arrays(passed, failed, total)
    call test_macropore_state_tracking(passed, failed, total)
    call test_macropore_domain_config(passed, failed, total)
    call test_macropore_flags(passed, failed, total)
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
    
    subroutine test_macropore_initialization(passed, failed, total)
        integer, intent(inout) :: passed, failed, total
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Macropore State Initialization ---')
        
        call swap_state_init(state, 40, 5, 2, 1)
        
        ! Test default values for domain storage
        call assert_eq_real(state%macro%VlMp, 0.0d0, 'VlMp = 0 initially', passed, failed, total)
        call assert_eq_real(state%macro%VlMpDm1, 0.0d0, 'VlMpDm1 = 0 initially', passed, failed, total)
        call assert_eq_real(state%macro%VlMpDm2, 0.0d0, 'VlMpDm2 = 0 initially', passed, failed, total)
        call assert_eq_real(state%macro%WaSrDm1, 0.0d0, 'WaSrDm1 = 0 initially', passed, failed, total)
        call assert_eq_real(state%macro%WaSrDm2, 0.0d0, 'WaSrDm2 = 0 initially', passed, failed, total)
        call assert_eq_real(state%macro%WaSrMp, 0.0d0, 'WaSrMp = 0 initially', passed, failed, total)
        
        ! Test flag defaults
        call assert_true(.not. state%macro%flmacropore, 'flmacropore = false initially', passed, failed, total)
        call assert_true(state%macro%flBegin, 'flBegin = true initially', passed, failed, total)
        
        call swap_state_finalize(state)
    end subroutine test_macropore_initialization
    
    subroutine test_macropore_domain_storage(passed, failed, total)
        integer, intent(inout) :: passed, failed, total
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Macropore Domain Storage ---')
        
        call swap_state_init(state, 40, 5, 2, 1)
        
        ! Set domain storage values
        state%macro%VlMp = 0.005d0
        state%macro%VlMpDm1 = 0.003d0
        state%macro%VlMpDm2 = 0.002d0
        state%macro%WaSrDm1 = 0.001d0
        state%macro%WaSrDm2 = 0.0005d0
        state%macro%WaSrDm1Ini = 0.0d0
        state%macro%WaSrDm2Ini = 0.0d0
        state%macro%WaLevDm1 = -50.0d0
        
        ! Verify values are stored
        call assert_eq_real(state%macro%VlMp, 0.005d0, 'VlMp = 0.005', passed, failed, total)
        call assert_eq_real(state%macro%VlMpDm1, 0.003d0, 'VlMpDm1 = 0.003', passed, failed, total)
        call assert_eq_real(state%macro%VlMpDm2, 0.002d0, 'VlMpDm2 = 0.002', passed, failed, total)
        call assert_eq_real(state%macro%WaSrDm1, 0.001d0, 'WaSrDm1 = 0.001', passed, failed, total)
        call assert_eq_real(state%macro%WaSrDm2, 0.0005d0, 'WaSrDm2 = 0.0005', passed, failed, total)
        call assert_eq_real(state%macro%WaLevDm1, -50.0d0, 'WaLevDm1 = -50', passed, failed, total)
        
        call swap_state_finalize(state)
    end subroutine test_macropore_domain_storage
    
    subroutine test_macropore_fluxes(passed, failed, total)
        integer, intent(inout) :: passed, failed, total
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Macropore Fluxes ---')
        
        call swap_state_init(state, 40, 5, 2, 1)
        
        ! Set flux values
        state%macro%QMaPo = 0.01d0
        state%macro%QRapDra = 0.005d0
        state%macro%QMpLatSs = 0.002d0
        state%macro%QInTopLatDm1 = 0.001d0
        state%macro%QInTopLatDm2 = 0.0008d0
        state%macro%QInTopVrtDm1 = 0.003d0
        state%macro%QInTopVrtDm2 = 0.002d0
        
        ! Verify values
        call assert_eq_real(state%macro%QMaPo, 0.01d0, 'QMaPo = 0.01', passed, failed, total)
        call assert_eq_real(state%macro%QRapDra, 0.005d0, 'QRapDra = 0.005', passed, failed, total)
        call assert_eq_real(state%macro%QMpLatSs, 0.002d0, 'QMpLatSs = 0.002', passed, failed, total)
        call assert_eq_real(state%macro%QInTopLatDm1, 0.001d0, 'QInTopLatDm1 = 0.001', passed, failed, total)
        call assert_eq_real(state%macro%QInTopLatDm2, 0.0008d0, 'QInTopLatDm2 = 0.0008', passed, failed, total)
        call assert_eq_real(state%macro%QInTopVrtDm1, 0.003d0, 'QInTopVrtDm1 = 0.003', passed, failed, total)
        call assert_eq_real(state%macro%QInTopVrtDm2, 0.002d0, 'QInTopVrtDm2 = 0.002', passed, failed, total)
        
        call swap_state_finalize(state)
    end subroutine test_macropore_fluxes
    
    subroutine test_macropore_cumulative(passed, failed, total)
        integer, intent(inout) :: passed, failed, total
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Macropore Cumulative Fluxes ---')
        
        call swap_state_init(state, 40, 5, 2, 1)
        
        ! Set cumulative flux values
        state%macro%cQMpLatSs = 10.0d0
        state%macro%cQMpOutDrRap = 5.0d0
        state%macro%cQMpInMtxSatDm1 = 3.0d0
        state%macro%cQMpInMtxSatDm2 = 2.0d0
        state%macro%cQMpOutMtxUnsDm1 = 1.5d0
        state%macro%cQMpOutMtxUnsDm2 = 1.0d0
        state%macro%cQMpInIntSatDm1 = 0.5d0
        state%macro%cQMpInIntSatDm2 = 0.3d0
        
        ! Verify values
        call assert_eq_real(state%macro%cQMpLatSs, 10.0d0, 'cQMpLatSs = 10', passed, failed, total)
        call assert_eq_real(state%macro%cQMpOutDrRap, 5.0d0, 'cQMpOutDrRap = 5', passed, failed, total)
        call assert_eq_real(state%macro%cQMpInMtxSatDm1, 3.0d0, 'cQMpInMtxSatDm1 = 3', passed, failed, total)
        call assert_eq_real(state%macro%cQMpInMtxSatDm2, 2.0d0, 'cQMpInMtxSatDm2 = 2', passed, failed, total)
        call assert_eq_real(state%macro%cQMpOutMtxUnsDm1, 1.5d0, 'cQMpOutMtxUnsDm1 = 1.5', passed, failed, total)
        call assert_eq_real(state%macro%cQMpOutMtxUnsDm2, 1.0d0, 'cQMpOutMtxUnsDm2 = 1.0', passed, failed, total)
        call assert_eq_real(state%macro%cQMpInIntSatDm1, 0.5d0, 'cQMpInIntSatDm1 = 0.5', passed, failed, total)
        call assert_eq_real(state%macro%cQMpInIntSatDm2, 0.3d0, 'cQMpInIntSatDm2 = 0.3', passed, failed, total)
        
        call swap_state_finalize(state)
    end subroutine test_macropore_cumulative
    
    subroutine test_macropore_incremental(passed, failed, total)
        integer, intent(inout) :: passed, failed, total
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Macropore Incremental Fluxes ---')
        
        call swap_state_init(state, 40, 5, 2, 1)
        
        ! Set incremental flux values
        state%macro%iQMpOutDrRap = 0.01d0
        state%macro%iQInTopLatDm1 = 0.002d0
        state%macro%iQInTopLatDm2 = 0.001d0
        state%macro%iQInTopVrtDm1 = 0.003d0
        state%macro%iQInTopVrtDm2 = 0.002d0
        state%macro%IWaSrDm1Beg = 0.0005d0
        state%macro%IWaSrDm2Beg = 0.0003d0
        
        ! Verify values
        call assert_eq_real(state%macro%iQMpOutDrRap, 0.01d0, 'iQMpOutDrRap = 0.01', passed, failed, total)
        call assert_eq_real(state%macro%iQInTopLatDm1, 0.002d0, 'iQInTopLatDm1 = 0.002', passed, failed, total)
        call assert_eq_real(state%macro%iQInTopLatDm2, 0.001d0, 'iQInTopLatDm2 = 0.001', passed, failed, total)
        call assert_eq_real(state%macro%iQInTopVrtDm1, 0.003d0, 'iQInTopVrtDm1 = 0.003', passed, failed, total)
        call assert_eq_real(state%macro%iQInTopVrtDm2, 0.002d0, 'iQInTopVrtDm2 = 0.002', passed, failed, total)
        call assert_eq_real(state%macro%IWaSrDm1Beg, 0.0005d0, 'IWaSrDm1Beg = 0.0005', passed, failed, total)
        call assert_eq_real(state%macro%IWaSrDm2Beg, 0.0003d0, 'IWaSrDm2Beg = 0.0003', passed, failed, total)
        
        call swap_state_finalize(state)
    end subroutine test_macropore_incremental
    
    subroutine test_macropore_arrays(passed, failed, total)
        integer, intent(inout) :: passed, failed, total
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Macropore Per-Compartment Arrays ---')
        
        call swap_state_init(state, 40, 5, 2, 1)
        
        ! Test that per-compartment arrays are allocated with correct size
        call assert_true(allocated(state%macro%DiPoCp), 'DiPoCp allocated', passed, failed, total)
        call assert_true(allocated(state%macro%FrArMtrx), 'FrArMtrx allocated', passed, failed, total)
        call assert_true(allocated(state%macro%VlMpDyCp), 'VlMpDyCp allocated', passed, failed, total)
        call assert_true(allocated(state%macro%VlMpStCp), 'VlMpStCp allocated', passed, failed, total)
        call assert_true(allocated(state%macro%VlMpStDm1), 'VlMpStDm1 allocated', passed, failed, total)
        call assert_true(allocated(state%macro%VlMpStDm2), 'VlMpStDm2 allocated', passed, failed, total)
        call assert_true(allocated(state%macro%QExcMpMtx), 'QExcMpMtx allocated', passed, failed, total)
        call assert_true(allocated(state%macro%SubsidCp), 'SubsidCp allocated', passed, failed, total)
        
        ! Check sizes (should be numnod = 40)
        call assert_eq_int(size(state%macro%DiPoCp), 40, 'size(DiPoCp) = 40', passed, failed, total)
        call assert_eq_int(size(state%macro%FrArMtrx), 40, 'size(FrArMtrx) = 40', passed, failed, total)
        call assert_eq_int(size(state%macro%VlMpDyCp), 40, 'size(VlMpDyCp) = 40', passed, failed, total)
        call assert_eq_int(size(state%macro%VlMpStCp), 40, 'size(VlMpStCp) = 40', passed, failed, total)
        
        ! Set some values
        state%macro%DiPoCp(1) = 0.1d0
        state%macro%DiPoCp(40) = 0.05d0
        state%macro%FrArMtrx(1) = 0.95d0
        
        call assert_eq_real(state%macro%DiPoCp(1), 0.1d0, 'DiPoCp(1) = 0.1', passed, failed, total)
        call assert_eq_real(state%macro%DiPoCp(40), 0.05d0, 'DiPoCp(40) = 0.05', passed, failed, total)
        call assert_eq_real(state%macro%FrArMtrx(1), 0.95d0, 'FrArMtrx(1) = 0.95', passed, failed, total)
        
        call swap_state_finalize(state)
    end subroutine test_macropore_arrays
    
    subroutine test_macropore_work_arrays(passed, failed, total)
        integer, intent(inout) :: passed, failed, total
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Macropore Work Arrays (from SAVE) ---')
        
        call swap_state_init(state, 40, 5, 2, 1)
        
        ! Test work arrays that were previously local SAVE variables
        call assert_true(allocated(state%macro%ICpBtDm), 'ICpBtDm allocated', passed, failed, total)
        call assert_true(allocated(state%macro%ICpTpWaSrDm), 'ICpTpWaSrDm allocated', passed, failed, total)
        call assert_true(allocated(state%macro%ArMpTpDm), 'ArMpTpDm allocated', passed, failed, total)
        call assert_true(allocated(state%macro%AwlCorFac), 'AwlCorFac allocated', passed, failed, total)
        call assert_true(allocated(state%macro%FrMpWalWet), 'FrMpWalWet allocated', passed, failed, total)
        call assert_true(allocated(state%macro%KDCrRlRef), 'KDCrRlRef allocated', passed, failed, total)
        call assert_true(allocated(state%macro%QExcMtxDmCp), 'QExcMtxDmCp allocated', passed, failed, total)
        call assert_true(allocated(state%macro%QInTopLatDm), 'QInTopLatDm allocated', passed, failed, total)
        call assert_true(allocated(state%macro%QInTopVrtDm), 'QInTopVrtDm allocated', passed, failed, total)
        call assert_true(allocated(state%macro%QOutDrRapCp), 'QOutDrRapCp allocated', passed, failed, total)
        call assert_true(allocated(state%macro%VlMpDm), 'VlMpDm allocated', passed, failed, total)
        call assert_true(allocated(state%macro%VlMpDmCp), 'VlMpDmCp allocated', passed, failed, total)
        call assert_true(allocated(state%macro%WaSrMpDm), 'WaSrMpDm allocated', passed, failed, total)
        call assert_true(allocated(state%macro%WaSrMpDmCp), 'WaSrMpDmCp allocated', passed, failed, total)
        call assert_true(allocated(state%macro%ZBtDm), 'ZBtDm allocated', passed, failed, total)
        call assert_true(allocated(state%macro%ZWaLevDm), 'ZWaLevDm allocated', passed, failed, total)
        call assert_true(allocated(state%macro%flDraTub), 'flDraTub allocated', passed, failed, total)
        call assert_true(allocated(state%macro%FlEndSrpEvt), 'FlEndSrpEvt allocated', passed, failed, total)
        
        ! Test array sizes - based on MaDm (domains) and MaCp (compartments)
        ! From arrays.fi: MaDm = 20, MaDr = 5
        call assert_eq_int(size(state%macro%ICpBtDm), 20, 'size(ICpBtDm) = 20 (MaDm)', passed, failed, total)
        call assert_eq_int(size(state%macro%ArMpTpDm), 20, 'size(ArMpTpDm) = 20 (MaDm)', passed, failed, total)
        call assert_eq_int(size(state%macro%VlMpDm), 20, 'size(VlMpDm) = 20 (MaDm)', passed, failed, total)
        call assert_eq_int(size(state%macro%KDCrRlRef), 5, 'size(KDCrRlRef) = 5 (MaDr)', passed, failed, total)
        call assert_eq_int(size(state%macro%flDraTub), 5, 'size(flDraTub) = 5 (MaDr)', passed, failed, total)
        
        ! Set some work array values
        state%macro%ICpBtDm(1) = 10
        state%macro%ICpBtDm(2) = 20
        state%macro%ArMpTpDm(1) = 0.03d0
        state%macro%ArMpTpDm(2) = 0.02d0
        state%macro%flDraTub(1) = .true.
        state%macro%flBegin = .false.
        
        call assert_eq_int(state%macro%ICpBtDm(1), 10, 'ICpBtDm(1) = 10', passed, failed, total)
        call assert_eq_int(state%macro%ICpBtDm(2), 20, 'ICpBtDm(2) = 20', passed, failed, total)
        call assert_eq_real(state%macro%ArMpTpDm(1), 0.03d0, 'ArMpTpDm(1) = 0.03', passed, failed, total)
        call assert_eq_real(state%macro%ArMpTpDm(2), 0.02d0, 'ArMpTpDm(2) = 0.02', passed, failed, total)
        call assert_true(state%macro%flDraTub(1), 'flDraTub(1) = true', passed, failed, total)
        call assert_true(.not. state%macro%flBegin, 'flBegin = false after set', passed, failed, total)
        
        call swap_state_finalize(state)
    end subroutine test_macropore_work_arrays
    
    subroutine test_macropore_state_tracking(passed, failed, total)
        integer, intent(inout) :: passed, failed, total
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Macropore State Tracking ---')
        
        call swap_state_init(state, 40, 5, 2, 1)
        
        ! Set state tracking values
        state%macro%ICpBtPerZon = 15
        state%macro%ICpSatGWl = 20
        state%macro%ICpSatPeGWl = 18
        state%macro%ICpTpPerZon = 5
        state%macro%ICpTpSatZon = 22
        state%macro%NnCrAr = 3
        
        ! Verify values
        call assert_eq_int(state%macro%ICpBtPerZon, 15, 'ICpBtPerZon = 15', passed, failed, total)
        call assert_eq_int(state%macro%ICpSatGWl, 20, 'ICpSatGWl = 20', passed, failed, total)
        call assert_eq_int(state%macro%ICpSatPeGWl, 18, 'ICpSatPeGWl = 18', passed, failed, total)
        call assert_eq_int(state%macro%ICpTpPerZon, 5, 'ICpTpPerZon = 5', passed, failed, total)
        call assert_eq_int(state%macro%ICpTpSatZon, 22, 'ICpTpSatZon = 22', passed, failed, total)
        call assert_eq_int(state%macro%NnCrAr, 3, 'NnCrAr = 3', passed, failed, total)
        
        call swap_state_finalize(state)
    end subroutine test_macropore_state_tracking
    
    subroutine test_macropore_domain_config(passed, failed, total)
        integer, intent(inout) :: passed, failed, total
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Macropore Domain Configuration ---')
        
        call swap_state_init(state, 40, 5, 2, 1)
        
        ! Set domain configuration values
        state%macro%NumDm = 2
        state%macro%NumSbDm = 4
        state%macro%IcTopMP = 3
        state%macro%NumLevRapDra = 2
        state%macro%Z_Tp = 0.0d0
        state%macro%Z_St = -50.0d0
        state%macro%Z_Ic = -100.0d0
        state%macro%Z_Ah = -25.0d0
        state%macro%ArMpTp = 0.05d0
        state%macro%ArMpSs = 0.03d0
        state%macro%KsatCovLay = 10.0d0
        state%macro%KsMpSs = 100.0d0
        state%macro%PpIcTpMp = 0.4d0
        state%macro%dtold = 0.001d0
        
        ! Verify values
        call assert_eq_int(state%macro%NumDm, 2, 'NumDm = 2', passed, failed, total)
        call assert_eq_int(state%macro%NumSbDm, 4, 'NumSbDm = 4', passed, failed, total)
        call assert_eq_int(state%macro%IcTopMP, 3, 'IcTopMP = 3', passed, failed, total)
        call assert_eq_int(state%macro%NumLevRapDra, 2, 'NumLevRapDra = 2', passed, failed, total)
        call assert_eq_real(state%macro%Z_Tp, 0.0d0, 'Z_Tp = 0', passed, failed, total)
        call assert_eq_real(state%macro%Z_St, -50.0d0, 'Z_St = -50', passed, failed, total)
        call assert_eq_real(state%macro%Z_Ic, -100.0d0, 'Z_Ic = -100', passed, failed, total)
        call assert_eq_real(state%macro%Z_Ah, -25.0d0, 'Z_Ah = -25', passed, failed, total)
        call assert_eq_real(state%macro%ArMpTp, 0.05d0, 'ArMpTp = 0.05', passed, failed, total)
        call assert_eq_real(state%macro%ArMpSs, 0.03d0, 'ArMpSs = 0.03', passed, failed, total)
        call assert_eq_real(state%macro%KsatCovLay, 10.0d0, 'KsatCovLay = 10', passed, failed, total)
        call assert_eq_real(state%macro%KsMpSs, 100.0d0, 'KsMpSs = 100', passed, failed, total)
        call assert_eq_real(state%macro%PpIcTpMp, 0.4d0, 'PpIcTpMp = 0.4', passed, failed, total)
        call assert_eq_real(state%macro%dtold, 0.001d0, 'dtold = 0.001', passed, failed, total)
        
        call swap_state_finalize(state)
    end subroutine test_macropore_domain_config
    
    subroutine test_macropore_flags(passed, failed, total)
        integer, intent(inout) :: passed, failed, total
        type(swap_state_t) :: state
        
        call log_info('test', '--- Test: Macropore Flags ---')
        
        call swap_state_init(state, 40, 5, 2, 1)
        
        ! Check initial flag values
        call assert_true(.not. state%macro%flmacropore, 'flmacropore = false initially', passed, failed, total)
        call assert_true(state%macro%flBegin, 'flBegin = true initially', passed, failed, total)
        call assert_true(.not. state%macro%FlDecMpRat, 'FlDecMpRat = false initially', passed, failed, total)
        call assert_true(.not. state%macro%flInitDraBas, 'flInitDraBas = false initially', passed, failed, total)
        
        ! Set flags
        state%macro%flmacropore = .true.
        state%macro%flBegin = .false.
        state%macro%FlDecMpRat = .true.
        state%macro%flInitDraBas = .true.
        state%macro%IDecMpRat = 5
        
        ! Verify changed values
        call assert_true(state%macro%flmacropore, 'flmacropore = true after set', passed, failed, total)
        call assert_true(.not. state%macro%flBegin, 'flBegin = false after set', passed, failed, total)
        call assert_true(state%macro%FlDecMpRat, 'FlDecMpRat = true after set', passed, failed, total)
        call assert_true(state%macro%flInitDraBas, 'flInitDraBas = true after set', passed, failed, total)
        call assert_eq_int(state%macro%IDecMpRat, 5, 'IDecMpRat = 5', passed, failed, total)
        
        call swap_state_finalize(state)
    end subroutine test_macropore_flags
    
    subroutine test_multiple_instances(passed, failed, total)
        integer, intent(inout) :: passed, failed, total
        type(swap_state_t) :: state1, state2, state3
        
        call log_info('test', '--- Test: Multiple Macropore Instances ---')
        
        ! Initialize three separate state instances with different sizes
        call swap_state_init(state1, 20, 3, 1, 1)
        call swap_state_init(state2, 40, 5, 2, 1)
        call swap_state_init(state3, 60, 7, 3, 2)
        
        ! Set different values in each instance
        state1%macro%VlMp = 0.001d0
        state1%macro%WaSrDm1 = 0.0005d0
        state1%macro%flmacropore = .true.
        state1%macro%NumDm = 1
        
        state2%macro%VlMp = 0.002d0
        state2%macro%WaSrDm1 = 0.001d0
        state2%macro%flmacropore = .true.
        state2%macro%NumDm = 2
        
        state3%macro%VlMp = 0.003d0
        state3%macro%WaSrDm1 = 0.0015d0
        state3%macro%flmacropore = .false.
        state3%macro%NumDm = 2
        
        ! Verify instances are independent - scalars
        call assert_eq_real(state1%macro%VlMp, 0.001d0, 'state1 VlMp independent', passed, failed, total)
        call assert_eq_real(state2%macro%VlMp, 0.002d0, 'state2 VlMp independent', passed, failed, total)
        call assert_eq_real(state3%macro%VlMp, 0.003d0, 'state3 VlMp independent', passed, failed, total)
        
        call assert_eq_real(state1%macro%WaSrDm1, 0.0005d0, 'state1 WaSrDm1 independent', passed, failed, total)
        call assert_eq_real(state2%macro%WaSrDm1, 0.001d0, 'state2 WaSrDm1 independent', passed, failed, total)
        call assert_eq_real(state3%macro%WaSrDm1, 0.0015d0, 'state3 WaSrDm1 independent', passed, failed, total)
        
        call assert_true(state1%macro%flmacropore, 'state1 flmacropore independent', passed, failed, total)
        call assert_true(state2%macro%flmacropore, 'state2 flmacropore independent', passed, failed, total)
        call assert_true(.not. state3%macro%flmacropore, 'state3 flmacropore independent', passed, failed, total)
        
        call assert_eq_int(state1%macro%NumDm, 1, 'state1 NumDm independent', passed, failed, total)
        call assert_eq_int(state2%macro%NumDm, 2, 'state2 NumDm independent', passed, failed, total)
        call assert_eq_int(state3%macro%NumDm, 2, 'state3 NumDm independent', passed, failed, total)
        
        ! Verify arrays are independent - different sizes based on numnod
        call assert_eq_int(size(state1%macro%DiPoCp), 20, 'state1 DiPoCp size 20', passed, failed, total)
        call assert_eq_int(size(state2%macro%DiPoCp), 40, 'state2 DiPoCp size 40', passed, failed, total)
        call assert_eq_int(size(state3%macro%DiPoCp), 60, 'state3 DiPoCp size 60', passed, failed, total)
        
        ! Modify state1 array and ensure state2/state3 are unaffected
        if (allocated(state1%macro%DiPoCp)) state1%macro%DiPoCp(1) = 0.5d0
        if (allocated(state2%macro%DiPoCp)) state2%macro%DiPoCp(1) = 0.6d0
        if (allocated(state3%macro%DiPoCp)) state3%macro%DiPoCp(1) = 0.7d0
        
        call assert_eq_real(state1%macro%DiPoCp(1), 0.5d0, 'state1 DiPoCp(1) still 0.5', passed, failed, total)
        call assert_eq_real(state2%macro%DiPoCp(1), 0.6d0, 'state2 DiPoCp(1) still 0.6', passed, failed, total)
        call assert_eq_real(state3%macro%DiPoCp(1), 0.7d0, 'state3 DiPoCp(1) still 0.7', passed, failed, total)
        
        call log_info('test', 'Multiple macropore instances are fully independent!')
        
        call swap_state_finalize(state1)
        call swap_state_finalize(state2)
        call swap_state_finalize(state3)
    end subroutine test_multiple_instances

end program test_macropore_state
