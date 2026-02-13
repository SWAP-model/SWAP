!> @file test_surfacewater_state.f90
!> @brief Unit tests for surface-water state module.
program test_surfacewater_state
    use surfacewater_state_mod, only: surfacewater_state_t, surfacewater_state_init, surfacewater_state_finalize, &
                              surfacewater_state_reset_cumulative, surfacewater_state_reset_intermediate
    implicit none

    type(surfacewater_state_t) :: state

    call test_surfacewater_init_and_dimensions(state)
    call test_surfacewater_flux_resets(state)
    call test_surfacewater_finalize(state)

    print *, 'test_surfacewater_state: PASS'

contains

    !> Validate allocation and canonical initialization.
    subroutine test_surfacewater_init_and_dimensions(sw)
        type(surfacewater_state_t), intent(inout) :: sw

        call surfacewater_state_init(sw, nmper=3, mamte=5, maowl=7, nrlevs=2)

        call assert_true(allocated(sw%wlsbak), 'wlsbak allocated')
        call assert_equal_int(4, size(sw%wlsbak), 'size(wlsbak)')

        call assert_true(allocated(sw%impend), 'impend allocated')
        call assert_equal_int(3, size(sw%impend), 'size(impend)')

        call assert_true(allocated(sw%gwlcrit), 'gwlcrit allocated')
        call assert_equal_int(3, size(sw%gwlcrit, 1), 'size(gwlcrit,1)')
        call assert_equal_int(5, size(sw%gwlcrit, 2), 'size(gwlcrit,2)')

        call assert_true(allocated(sw%owltab), 'owltab allocated')
        call assert_equal_int(2, size(sw%owltab, 1), 'size(owltab,1)')
        call assert_equal_int(14, size(sw%owltab, 2), 'size(owltab,2)')

        call assert_equal_int(3, sw%nmper, 'nmper initialized')
        call assert_equal_real(0.0d0, sw%swst, 1.0d-12, 'swst default')
        call assert_equal_real(0.0d0, sw%cqdrd, 1.0d-12, 'cqdrd default')
        call assert_true(.not. sw%overfl, 'overfl default false')

    end subroutine test_surfacewater_init_and_dimensions

    !> Validate selective reset procedures.
    subroutine test_surfacewater_flux_resets(sw)
        type(surfacewater_state_t), intent(inout) :: sw

        sw%cqdrd = 11.0d0
        sw%cwsupp = 12.0d0
        sw%cwout = 13.0d0
        sw%runots = 2.0d0
        sw%qdrd = 3.0d0

        call surfacewater_state_reset_cumulative(sw)
        call assert_equal_real(0.0d0, sw%cqdrd, 1.0d-12, 'cqdrd cumulative reset')
        call assert_equal_real(0.0d0, sw%cwsupp, 1.0d-12, 'cwsupp cumulative reset')
        call assert_equal_real(0.0d0, sw%cwout, 1.0d-12, 'cwout cumulative reset')
        call assert_equal_real(2.0d0, sw%runots, 1.0d-12, 'runots unchanged by cumulative reset')

        call surfacewater_state_reset_intermediate(sw)
        call assert_equal_real(0.0d0, sw%runots, 1.0d-12, 'runots intermediate reset')
        call assert_equal_real(3.0d0, sw%qdrd, 1.0d-12, 'qdrd unchanged by intermediate reset')

    end subroutine test_surfacewater_flux_resets

    !> Validate deallocation and scalar reset at finalize.
    subroutine test_surfacewater_finalize(sw)
        type(surfacewater_state_t), intent(inout) :: sw

        sw%wlp = -25.0d0
        sw%overfl = .true.

        call surfacewater_state_finalize(sw)

        call assert_true(.not. allocated(sw%wlsbak), 'wlsbak deallocated')
        call assert_true(.not. allocated(sw%impend), 'impend deallocated')
        call assert_true(.not. allocated(sw%owltab), 'owltab deallocated')
        call assert_equal_real(0.0d0, sw%wlp, 1.0d-12, 'wlp reset after finalize')
        call assert_true(.not. sw%overfl, 'overfl reset after finalize')

    end subroutine test_surfacewater_finalize

    !> Assert integer equality.
    subroutine assert_equal_int(expected, actual, label)
        integer, intent(in) :: expected, actual
        character(len=*), intent(in) :: label

        if (actual /= expected) then
            write(*,'(A,1X,A,1X,I0,1X,A,1X,I0)') 'Assertion failed:', trim(label), actual, '/=', expected
            error stop 1
        end if
    end subroutine assert_equal_int

    !> Assert real equality with absolute tolerance.
    subroutine assert_equal_real(expected, actual, tol, label)
        real(8), intent(in) :: expected, actual, tol
        character(len=*), intent(in) :: label

        if (abs(actual - expected) > tol) then
            write(*,'(A,1X,A,1X,ES16.8,1X,A,1X,ES16.8)') 'Assertion failed:', trim(label), actual, '/=', expected
            error stop 1
        end if
    end subroutine assert_equal_real

    !> Assert logical truth.
    subroutine assert_true(condition, label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label

        if (.not. condition) then
            write(*,'(A,1X,A)') 'Assertion failed:', trim(label)
            error stop 1
        end if
    end subroutine assert_true

end program test_surfacewater_state
