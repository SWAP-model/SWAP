!> @file timestep_control_mod.f90
!! Timestep-control flags migrated out of variables.f90.
!!
!! `fldecdt` — decrease-timestep signal set by:
!!   - SurfaceWater/WLEVBAL when oscillation or ponding limit is hit
!!     (now also exposed as intent(out) :: request_smaller_dt on SurfaceWater(2/3))
!!   - headcalc when Richards-equation iteration does not converge
!!
!! Initialized in TimeControl(task=1); reset in TimeControl(task=3) after
!! reducing dt. Readers: swap_main (gates subsystem calls), timecontrol,
!! OutputModflow diagnostic loop.
!!
!! See SS-SWST Phase 2 Task 1 and ADR 0030.

module timestep_control_mod
   implicit none
   public

   logical :: fldecdt = .false.   ! Flag indicating decrease of time step

end module timestep_control_mod
