! swap_main.f90
! SS-DRV Phase 1: thin driver around module swap_mod.
! The outer time loop lives here so the same module can be driven
! by BMI (swap_bmi_mod) one timestep at a time.
program swap_main

   use swap_mod,        only: swap_init, swap_run_step, swap_close
   use swap_state_mod,  only: swap_state_t
   use swap_config_mod, only: swap_config_t
   use swap_log,        only: log_init, log_close, LOGLEVEL_INFO
   implicit none

   type(swap_state_t)           :: state
   type(swap_config_t), target  :: config  ! target: crop_config_global pointer set inside swap_init
   logical                      :: fileopen

   call log_init(log_level=LOGLEVEL_INFO, log_file='swap_swap.log')

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
