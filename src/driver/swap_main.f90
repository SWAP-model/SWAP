! swap_main.f90
! SS-DRV Phase 1: thin driver around module swap_mod.
! The outer time loop lives here so the same module can be driven
! by BMI (swap_bmi_mod) one timestep at a time.
program swap_main

   use swap_mod,        only: swap_init, swap_run_step, swap_close
   use swap_state_mod,  only: swap_state_t
   use swap_config_mod, only: swap_config_t
   use swap_log,         only: log_close
   use diagnostics_mod,  only: diag_overrides_t, default_cli_config, &
                               read_logging_overrides_from_file, init_logging
   implicit none

   type(swap_state_t)           :: state
   type(swap_config_t), target  :: config  ! target retained for callers that take a pointer into config (crop_config_global retired)
   logical                      :: fileopen
   type(diag_overrides_t)       :: toml_ov

   call read_logging_overrides_from_file('swap.toml', toml_ov)
   call init_logging(default_cli_config(), toml_ov)

   call swap_init('swap.toml', state, config)
   do while (.not. state%timecontrol%flRunEnd)
      call swap_run_step(state, config)
   end do
   call swap_close(state, config)

   write(*,'(a)')' Swap normal completion!'
   call log_close()

   ! IO-OUT/E: CloseTempFil inlined here (swapoutput.f90 deleted). The
   ! unit-20 close-with-DELETE is the project's own scratch-file cleanup
   ! (TTutil scratch retired with ADR 0023); no other code opens unit 20.
   inquire(unit=20, opened=fileopen)
   if (fileopen) close(20, status='DELETE')

   stop 100

end program swap_main
