!> @file crop_irrigation_state.f90
!! SS-GR-CROP: typed crop runtime state — irrigation runtime/persistent state.
!! Populated by config_to_variables%apply_irrigation_ssdi (mode 0/1) from
!! state%cfg%irrigation%ssdi at startup; mutated at runtime by
!! src/crop/irrigation.f90 SSDI_irrigation(2)/SSDI_irrigation(9).
!!
!! [GR-CROP 2026-05-25] hosts the 16 SSDI persistent state fields that
!! were SAVE-state on irrigation.f90 module variables (later legacy
!! globals in variables.f90 module). The corresponding `*_irr` legacy
!! globals retire in the same commit.
!!
!! dt_SSDI_event lives on the legacy global side for now — cross-file
!! consumer is src/core/timecontrol_mod.f90 which is outside the crop
!! sub-arc's allow-list. A future timecontrol-side migration will rehome
!! it here (state%crop%irrigation%dt_SSDI_event).
module crop_irrigation_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   use swap_array_dimensions, only: mairg
   implicit none
   private
   public :: crop_irrigation_state_t

   type :: crop_irrigation_state_t

      ! SSDI configuration (init-once from state%cfg%irrigation%ssdi via
      ! apply_irrigation_ssdi in config_to_variables.f90)
      integer      :: swssdi            = 0          !! SSDI active (0=no, 1=yes)
      integer      :: nod_ssdi(2)       = 0          !! Upper and lower nodes for SSDI
      integer      :: ssdi_schedule     = 0          !! Schedule type (0=fixed dates, 1=internal)
      integer      :: ssdi_sched_type   = 0          !! Internal schedule type (1=Tact/Tpot, 2=h, 3=theta)
      integer      :: nod_ssdi_sensor   = 0          !! Sensor node (if ssdi_sched_type > 1)
      real(real64) :: ssdi_threshold    = 0.0_real64 !! Threshold value for scheduling
      real(real64) :: ssdi_threshold_z  = 0.0_real64 !! Depth for threshold value (cm)
      real(real64) :: ssdi_amount       = 0.0_real64 !! Amount of scheduled irrigation (cm)
      real(real64) :: ssdi_appl_rate    = 0.0_real64 !! Application rate (cm/d)
      integer      :: sw_interval       = 0          !! Switch for minimum interval
      integer      :: days_interval     = 1          !! Minimum days between applications

      ! SSDI runtime cursor/counter state
      integer      :: days_counter      = 366        !! Days since previous application
      integer      :: nirri             = 1          !! SSDI counter / entry point into ssdi_date

      ! Fixed-date schedule tables (mode 0)
      real(real64) :: ssdi_date(mairg)     = 0.0_real64 !! Fixed irrigation dates (days-since-1900)
      real(real64) :: ssdi_rate_f(mairg)   = 0.0_real64 !! Fixed irrigation rates (cm/d)
      real(real64) :: ssdi_amount_f(mairg) = 0.0_real64 !! Fixed irrigation amounts (cm)

      ! Surface fixed-irrigation events (populated by config_to_variables%apply_irrigation
      ! from state%cfg%irrigation%fixed_events / fixed_events_file when swirfix == 1).
      integer      :: nirri_fixed              = 1          !! Cursor into fixed-irrigation event arrays
      real(real64) :: irdate(mairg)            = 0.0_real64 !! Fixed irrigation dates (days-since-1900)
      real(real64) :: irdepth(mairg)           = 0.0_real64 !! Fixed irrigation depths (cm)
      real(real64) :: irconc(mairg)            = 0.0_real64 !! Fixed irrigation concentrations (M/L3)
      integer      :: irtype(mairg)            = 0          !! Fixed irrigation types (0=sprinkler, 1=surface)

   contains
      procedure :: init => crop_irrigation_state_init
   end type crop_irrigation_state_t

contains

   subroutine crop_irrigation_state_init(self)
      class(crop_irrigation_state_t), intent(inout) :: self
      ! Defaults set on type declaration; apply_irrigation_ssdi populates
      ! at startup when state%cfg%irrigation%swssdi == 1.
   end subroutine

end module crop_irrigation_state_mod
