program test_readcrop_toml
    use swap_state_mod, only: swap_state_t, swap_state_init, swap_state_finalize
    use swap_config_mod, only: swap_config_t, swap_config_init, swap_config_finalize
    use readcroptoml_mod, only: ReadCropToml_state
    implicit none

    type(swap_state_t) :: state
    type(swap_config_t) :: cfg

    call swap_config_init(cfg, 5, 2, 1, 1)
    call swap_state_init(state, cfg)

    call ReadCropToml_state(state, 'tests/swap-cases/1.1.hupselbrook-toml/maizes.toml')
    if (trim(state%crop%model) /= 'fixed') error stop 1
    if (abs(state%crop%dvsend - 3.0d0) > 1.0d-12) error stop 2
    if (state%crop%swcf /= 2) error stop 3
    if (state%crop%swoxygen /= 1) error stop 4
    if (.not. allocated(state%crop%rdtb_table)) error stop 5
    if (state%crop%scheduled_irrigation%schedule /= 0) error stop 14
    if (state%crop%scheduled_irrigation%isuas /= 1) error stop 15
    if (state%crop%scheduled_irrigation%startirr(1) /= 30) error stop 16
    if (state%crop%scheduled_irrigation%startirr(2) /= 3) error stop 17
    if (state%crop%scheduled_irrigation%endirr(1) /= 31) error stop 18
    if (state%crop%scheduled_irrigation%endirr(2) /= 12) error stop 19
    if (abs(state%crop%scheduled_irrigation%cirrs - 0.0d0) > 1.0d-12) error stop 20

    call ReadCropToml_state(state, 'tests/swap-cases/1.1.hupselbrook-toml/potatod.toml')
    if (trim(state%crop%model) /= 'wofost') error stop 6
    if (abs(state%crop%rdi - 10.0d0) > 1.0d-12) error stop 7
    if (state%crop%swdrought /= 1) error stop 8
    if (.not. allocated(state%crop%rlwtb_table)) error stop 9

    call ReadCropToml_state(state, 'tests/swap-cases/1.1.hupselbrook-toml/grassd.toml')
    if (trim(state%crop%model) /= 'wofost_grass') error stop 10
    if (state%crop%swrd /= 3) error stop 11
    if (abs(state%rootextract%hlim3h + 200.0d0) > 1.0d-12) error stop 12
    if (.not. allocated(state%crop%cf_table)) error stop 13

    call swap_state_finalize(state)
    call swap_config_finalize(cfg)

end program test_readcrop_toml
