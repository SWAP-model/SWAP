program test_readdra_toml
    use drainage_state_mod, only: drainage_config_t
    use readdrainagetoml_mod, only: ReadDrainageToml_config
    implicit none

    type(drainage_config_t) :: config

    call ReadDrainageToml_config(config, 'tests/swap-cases/1.1.hupselbrook-toml/swap.dra.toml')

    if (config%dramet <= 0) error stop 1

end program test_readdra_toml
