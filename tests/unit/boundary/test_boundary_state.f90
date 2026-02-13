!> Standalone unit test for boundary state module lifecycle.
program test_boundary_state
    use boundary_state_mod, only: boundary_state_t, boundary_state_init, boundary_state_finalize
    implicit none

    type(boundary_state_t) :: b

    call boundary_state_init(b)

    if (.not. allocated(b%gwltab)) error stop 1
    if (.not. allocated(b%haqtab)) error stop 1
    if (.not. allocated(b%qbotab)) error stop 1
    if (.not. allocated(b%hbotab)) error stop 1
    if (.not. allocated(b%runonarr)) error stop 1
    if (.not. allocated(b%pondmxtab)) error stop 1

    b%swbotb = 6
    b%qbot = -0.01d0
    b%pondmx = 2.0d0

    call boundary_state_finalize(b)

    if (allocated(b%gwltab)) error stop 1
    if (allocated(b%haqtab)) error stop 1
    if (allocated(b%qbotab)) error stop 1
    if (allocated(b%hbotab)) error stop 1
    if (allocated(b%runonarr)) error stop 1
    if (allocated(b%pondmxtab)) error stop 1

    if (b%swbotb /= 0) error stop 1
    if (abs(b%qbot) > 1.0d-12) error stop 1
    if (abs(b%pondmx) > 1.0d-12) error stop 1

    print *, 'test_boundary_state: PASS'
end program test_boundary_state
