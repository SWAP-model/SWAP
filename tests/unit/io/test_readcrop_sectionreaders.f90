program test_readcrop_sectionreaders
    use tomlf, only: toml_table
    use toml_helpers_mod, only: get_opt_table
    use swap_state_mod, only: swap_state_t, swap_state_init, swap_state_finalize
    use swap_config_mod, only: swap_config_t, swap_config_init, swap_config_finalize
    use readcroptoml_mod, only: ReadCropToml_open, read_crop_irrigation
    implicit none

    type(swap_state_t) :: state
    type(swap_config_t) :: cfg
    type(toml_table), allocatable :: doc
    type(toml_table), pointer :: crop_tab

    call swap_config_init(cfg, 5, 2, 1, 1)
    call swap_state_init(state, cfg)

    call ReadCropToml_open(doc, 'tests/swap-cases/1.1.hupselbrook-toml/maizes.toml')

    call get_opt_table(doc, 'crop', crop_tab, err_label='test_readcrop_sectionreaders/crop')
    if (.not. associated(crop_tab)) error stop 1

    call read_crop_irrigation(crop_tab, state)

    if (state%crop%scheduled_irrigation%schedule /= 0) error stop 2
    if (state%crop%scheduled_irrigation%isuas /= 1) error stop 3
    if (state%crop%scheduled_irrigation%startirr(1) /= 30) error stop 4
    if (state%crop%scheduled_irrigation%startirr(2) /= 3) error stop 5
    if (state%crop%scheduled_irrigation%endirr(1) /= 31) error stop 6
    if (state%crop%scheduled_irrigation%endirr(2) /= 12) error stop 7

    call swap_state_finalize(state)
    call swap_config_finalize(cfg)

end program test_readcrop_sectionreaders
