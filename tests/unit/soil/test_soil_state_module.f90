!> @file test_soil_state_module.f90
!> @brief Unit tests for standalone soil_state_mod lifecycle routines
program test_soil_state_module
    use soil_state_mod, only: soil_state_t, soil_state_init, soil_state_finalize, &
                              soil_state_reset_cumulative, soil_state_reset_intermediate
    implicit none

    type(soil_state_t) :: soil

    call test_init_allocates_and_sets_sizes()
    call test_reset_accumulators()
    call test_finalize_deallocates_and_resets()

    print *, 'test_soil_state_module: PASS'

contains

    !> Verify allocation and expected dimensions after initialization.
    subroutine test_init_allocates_and_sets_sizes()
        call soil_state_init(soil, numnod=5, numlay=2)

        call assert_true(allocated(soil%h), 'soil%h allocated')
        call assert_true(allocated(soil%q), 'soil%q allocated')
        call assert_true(allocated(soil%paramvg), 'soil%paramvg allocated')
        call assert_true(allocated(soil%cofgen), 'soil%cofgen allocated')
        call assert_true(allocated(soil%layer), 'soil%layer allocated')
        call assert_true(allocated(soil%botcom), 'soil%botcom allocated')

        call assert_equal_int(5, size(soil%h), 'size(soil%h)')
        call assert_equal_int(6, size(soil%q), 'size(soil%q)')
        call assert_equal_int(21, size(soil%paramvg, 1), 'size(soil%paramvg,1)')
        call assert_equal_int(2, size(soil%paramvg, 2), 'size(soil%paramvg,2)')
        call assert_equal_int(21, size(soil%cofgen, 1), 'size(soil%cofgen,1)')
        call assert_equal_int(5, size(soil%cofgen, 2), 'size(soil%cofgen,2)')
        call assert_equal_int(5, soil%numnod, 'soil%numnod')
        call assert_equal_int(2, soil%numlay, 'soil%numlay')
    end subroutine test_init_allocates_and_sets_sizes

    !> Verify cumulative and intermediate reset routines only clear accumulators.
    subroutine test_reset_accumulators()
        soil%cqbot = 1.0d0
        soil%cqbotdo = 2.0d0
        soil%cqrot = 3.0d0
        soil%crunoff = 4.0d0
        soil%iqbot = 5.0d0
        soil%iqdra = 6.0d0
        soil%irunon = 7.0d0

        call soil_state_reset_cumulative(soil)
        call assert_equal_real(0.0d0, soil%cqbot, 1.0d-12, 'soil%cqbot reset')
        call assert_equal_real(0.0d0, soil%cqbotdo, 1.0d-12, 'soil%cqbotdo reset')
        call assert_equal_real(0.0d0, soil%cqrot, 1.0d-12, 'soil%cqrot reset')
        call assert_equal_real(0.0d0, soil%crunoff, 1.0d-12, 'soil%crunoff reset')

        call soil_state_reset_intermediate(soil)
        call assert_equal_real(0.0d0, soil%iqbot, 1.0d-12, 'soil%iqbot reset')
        call assert_equal_real(0.0d0, soil%iqdra, 1.0d-12, 'soil%iqdra reset')
        call assert_equal_real(0.0d0, soil%irunon, 1.0d-12, 'soil%irunon reset')
    end subroutine test_reset_accumulators

    !> Verify finalize deallocates arrays and restores canonical defaults.
    subroutine test_finalize_deallocates_and_resets()
        soil%gwl = -123.0d0
        soil%pondmx = 1.5d0
        soil%swhyst = 2
        soil%flwarn_hc = .false.

        call soil_state_finalize(soil)

        call assert_true(.not. allocated(soil%h), 'soil%h deallocated')
        call assert_true(.not. allocated(soil%q), 'soil%q deallocated')
        call assert_true(.not. allocated(soil%paramvg), 'soil%paramvg deallocated')
        call assert_true(.not. allocated(soil%cofgen), 'soil%cofgen deallocated')
        call assert_true(.not. allocated(soil%layer), 'soil%layer deallocated')
        call assert_true(.not. allocated(soil%botcom), 'soil%botcom deallocated')

        call assert_equal_real(0.0d0, soil%gwl, 1.0d-12, 'soil%gwl reset')
        call assert_equal_real(0.0d0, soil%pondmx, 1.0d-12, 'soil%pondmx reset')
        call assert_equal_int(0, soil%swhyst, 'soil%swhyst reset')
        call assert_true(soil%flwarn_hc, 'soil%flwarn_hc reset true')
        call assert_equal_int(0, soil%numnod, 'soil%numnod reset')
        call assert_equal_int(0, soil%numlay, 'soil%numlay reset')
    end subroutine test_finalize_deallocates_and_resets

    !> Assert logical truth.
    subroutine assert_true(condition, label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label

        if (.not. condition) then
            write(*,'(A,1X,A)') 'Assertion failed:', trim(label)
            error stop 1
        end if
    end subroutine assert_true

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

end program test_soil_state_module
