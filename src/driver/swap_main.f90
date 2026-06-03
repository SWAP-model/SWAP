! swap_main.f90
! SS-DRV Phase 1: thin driver around module swap_mod.
! The outer time loop lives here so the same module can be driven
! by BMI (swap_bmi_mod) one timestep at a time.
program swap_main

   use swap_mod,        only: swap_init, swap_run_step, swap_close
   use swap_state_mod,  only: swap_state_t
   use swap_config_mod, only: swap_config_t
   use swap_log,         only: log_init, log_close
   use diagnostics_mod,  only: diagnostics_config_t, diag_overrides_t, &
                               default_cli_config, read_env_overrides,  &
                               resolve_diagnostics_config
   implicit none

   type(swap_state_t)           :: state
   type(swap_config_t), target  :: config  ! target retained for callers that take a pointer into config (crop_config_global retired)
   logical                      :: fileopen

   block
      type(diagnostics_config_t) :: dcfg
      type(diag_overrides_t)     :: none_ov, env_ov
      call read_env_overrides(env_ov)
      dcfg = resolve_diagnostics_config(default_cli_config(), none_ov, env_ov, none_ov)
      if (allocated(dcfg%log_file)) then
         call log_init(log_level=dcfg%level, log_file=dcfg%log_file, &
                       to_stdout=dcfg%to_stdout, to_stderr=dcfg%to_stderr, &
                       timestamps=dcfg%timestamps)
      else
         call log_init(log_level=dcfg%level, to_stdout=dcfg%to_stdout, &
                       to_stderr=dcfg%to_stderr, timestamps=dcfg%timestamps)
      end if
   end block

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
