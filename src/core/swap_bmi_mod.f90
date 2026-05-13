!> @file swap_bmi_mod.f90
!! SS-DRV Phase 1: minimal Basic Model Interface (BMI) facade.
!! Exposes initialize / update / finalize / get_value_double /
!! get_current_time as bind(C) procedures. All other CSDMS BMI v2.0
!! methods are stubbed (correct C signature, body returns 0 / empty)
!! and marked with `! BMI-STUB Phase 2` for the next arc to enumerate.
!! Holds a single module-level (state, config) pair — multi-column
!! support is Phase 3 (swap_ensemble_mod).
module swap_bmi_mod
   use iso_c_binding,   only: c_char, c_double, c_int, c_size_t, c_null_char, c_ptr
   use swap_mod,        only: swap_init, swap_run_step, swap_close
   use swap_state_mod,  only: swap_state_t
   use swap_config_mod, only: swap_config_t
   implicit none
   private

   type(swap_state_t),          save :: bmi_state
   type(swap_config_t), target, save :: bmi_config

contains

   !----------------------------------------------------------------------
   ! Lifecycle
   !----------------------------------------------------------------------

   function bmi_initialize(config_file, n) result(rc) bind(C, name='initialize')
      character(kind=c_char), intent(in)    :: config_file(*)
      integer(c_int),  value, intent(in)    :: n
      integer(c_int)                        :: rc
      character(len=256) :: f_config_file
      call c_to_f_string(config_file, f_config_file)
      call swap_init(trim(f_config_file), bmi_state, bmi_config)
      rc = 0
   end function bmi_initialize

   function bmi_update() result(rc) bind(C, name='update')
      integer(c_int) :: rc
      call swap_run_step(bmi_state, bmi_config)
      rc = 0
   end function bmi_update

   function bmi_finalize() result(rc) bind(C, name='finalize')
      integer(c_int) :: rc
      call swap_close(bmi_state, bmi_config)
      rc = 0
   end function bmi_finalize

   !----------------------------------------------------------------------
   ! Time
   !----------------------------------------------------------------------

   function bmi_get_current_time(t) result(rc) bind(C, name='get_current_time')
      real(c_double), intent(out) :: t
      integer(c_int)              :: rc
      t = bmi_state%timecontrol%t1900
      rc = 0
   end function bmi_get_current_time

   !----------------------------------------------------------------------
   ! Variable accessors — two sentinel variables for Phase 1
   !----------------------------------------------------------------------

   function bmi_get_value_double(var_name, n, dest) result(rc) bind(C, name='get_value_double')
      character(kind=c_char), intent(in)    :: var_name(*)
      integer(c_int),  value, intent(in)    :: n
      real(c_double),         intent(out)   :: dest(n)
      integer(c_int)                        :: rc
      character(len=64) :: name
      integer :: m
      call c_to_f_string(var_name, name)
      select case (trim(name))
      case ('soil_water_content')
         m = min(n, size(bmi_state%soilwater%theta))
         dest(1:m) = bmi_state%soilwater%theta(1:m)
         if (m < n) dest(m+1:n) = 0.0_c_double
         rc = 0
      case ('pressure_head')
         m = min(n, size(bmi_state%soilwater%h))
         dest(1:m) = bmi_state%soilwater%h(1:m)
         if (m < n) dest(m+1:n) = 0.0_c_double
         rc = 0
      case ('soil_temperature')
         m = min(n, size(bmi_state%heat%tsoil))
         dest(1:m) = bmi_state%heat%tsoil(1:m)
         if (m < n) dest(m+1:n) = 0.0_c_double
         rc = 0
      case ('groundwater_level')
         ! state%soilwater%gwl verified: soilwater_state.f90 line 152
         if (n >= 1) dest(1) = bmi_state%soilwater%gwl
         if (n > 1)  dest(2:n) = 0.0_c_double
         rc = 0
      case ('bottom_flux')
         ! state%soilwater%qbot verified: soilwater_state.f90 line 79
         if (n >= 1) dest(1) = bmi_state%soilwater%qbot
         if (n > 1)  dest(2:n) = 0.0_c_double
         rc = 0
      case ('actual_evapotranspiration')
         ! state%soilwater%iqrot verified: soilwater_state.f90 line 182
         if (n >= 1) dest(1) = bmi_state%soilwater%iqrot
         if (n > 1)  dest(2:n) = 0.0_c_double
         rc = 0
      case ('recharge')
         ! sign: positive = downward into aquifer (recharge), hence -qbot
         if (n >= 1) dest(1) = -bmi_state%soilwater%qbot
         if (n > 1)  dest(2:n) = 0.0_c_double
         rc = 0
      case ('surface_runoff')
         ! state%soilwater%runots verified: soilwater_state.f90 line 72
         if (n >= 1) dest(1) = bmi_state%soilwater%runots
         if (n > 1)  dest(2:n) = 0.0_c_double
         rc = 0
      case default
         rc = 1
      end select
   end function bmi_get_value_double

   !----------------------------------------------------------------------
   ! BMI-STUB Phase 2: spec-required methods, currently return 0 / empty.
   ! Each carries the correct C signature so consumers can compile, but
   ! does nothing useful until Phase 2 wires the full variable registry.
   !----------------------------------------------------------------------

   function bmi_get_end_time(t) result(rc) bind(C, name='get_end_time')
      real(c_double), intent(out) :: t
      integer(c_int)              :: rc
      t = 0.0_c_double            ! BMI-STUB Phase 2
      rc = 0
   end function bmi_get_end_time

   function bmi_get_time_step(dt) result(rc) bind(C, name='get_time_step')
      real(c_double), intent(out) :: dt
      integer(c_int)              :: rc
      dt = bmi_state%timecontrol%dt   ! cheap real impl
      rc = 0
   end function bmi_get_time_step

   function bmi_set_value_double(var_name, n, src) result(rc) bind(C, name='set_value_double')
      character(kind=c_char), intent(in) :: var_name(*)
      integer(c_int),  value, intent(in) :: n
      real(c_double),         intent(in) :: src(n)
      integer(c_int)                     :: rc
      character(len=64) :: name
      call c_to_f_string(var_name, name)
      rc = 0
      select case (trim(name))
      case ('groundwater_level_imposed')
         ! Sets bottom Dirichlet BC: overwrite last node pressure head.
         ! Full MODFLOW exchange plumbing (swdrasur/BoundBottom) is a follow-on.
         if (n >= 1) bmi_state%soilwater%h(size(bmi_state%soilwater%h)) = src(1)
      case ('bottom_flux_imposed')
         ! Sets bottom Neumann BC via state%soilwater%qbot.
         if (n >= 1) bmi_state%soilwater%qbot = src(1)
      case default
         rc = 2   ! not settable
      end select
   end function bmi_set_value_double

   !----------------------------------------------------------------------
   ! Helpers
   !----------------------------------------------------------------------

   subroutine c_to_f_string(c_str, f_str)
      character(kind=c_char), intent(in)  :: c_str(*)
      character(len=*),       intent(out) :: f_str
      integer :: i
      f_str = ' '
      do i = 1, len(f_str)
         if (c_str(i) == c_null_char) exit
         f_str(i:i) = c_str(i)
      end do
   end subroutine c_to_f_string

end module swap_bmi_mod
