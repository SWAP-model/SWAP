!> XMI (Deltares eXtended Model Interface) verbs for the SWAP ensemble —
!! prepare/solve/finalize-solve timestep control plus the zero-copy
!! get_value_ptr exchange protocol, matching what xmipy's XmiWrapper /
!! imod_coupler's SwapWrapper call.
!!
!! T1-G′ sub-arc 2 (ADR 0050): this module no longer re-implements the BMI
!! lifecycle / time / component-info / variable-metadata C names. Those 14
!! formerly-colliding symbols (initialize/update/finalize/get_current_time/…)
!! now have a single mode-aware definition in swap_bmi_mod, backed by the same
!! ensemble (coupled mode when an ensemble.txt sidecar is present). With the
!! collision gone, swap_xmi_mod links into the SAME shared library as
!! swap_bmi_mod/swap_capi_mod: libswap.so — the one library every consumer
!! (cffi, ctypes, xmipy, imod_coupler) loads.
!!
!! Scope note: the module-global error state (error_mod library_mode) and the
!! module-global ensemble (swap_ensemble_mod) are single-kernel-per-process.
!! That matches imod_coupler's one-kernel-per-library model; true
!! multi-instance support would require per-instance error + ensemble state
!! (handle-based, a future arc).
module swap_xmi_mod
   use iso_c_binding
   use swap_ensemble_mod
   use error_mod, only: library_fatal_raised
   use swap_c_strings_mod, only: c_to_f_string, f_to_c_string
   implicit none
   private

contains

   ! ---- XMI time-step / solve control ----------------------------------
   integer(c_int) function xmi_prepare_time_step(dt) result(rc) bind(C, name='prepare_time_step')
      real(c_double), intent(in) :: dt        ! by reference (xmipy passes byref)
      associate (unused => dt); end associate ! silence unused-arg warning
      rc = 0
   end function xmi_prepare_time_step

   integer(c_int) function xmi_finalize_time_step() result(rc) bind(C, name='finalize_time_step')
      rc = 0
   end function xmi_finalize_time_step

   integer(c_int) function xmi_prepare_solve(cid) result(rc) bind(C, name='prepare_solve')
      integer(c_int), intent(in) :: cid
      associate (unused => cid); end associate
      rc = 0
   end function xmi_prepare_solve

   integer(c_int) function xmi_solve(cid, has_converged) result(rc) bind(C, name='solve')
      integer(c_int), intent(in) :: cid
      integer(c_int), intent(out) :: has_converged
      associate (unused => cid); end associate
      rc = ensemble_step_day()      ! the actual column loop (one day per call)
      if (rc /= 0) ensemble_last_error = 'ensemble_step_day failed'
      has_converged = 1
   end function xmi_solve

   integer(c_int) function xmi_finalize_solve(cid) result(rc) bind(C, name='finalize_solve')
      integer(c_int), intent(in) :: cid
      associate (unused => cid); end associate
      rc = 0
   end function xmi_finalize_solve

   integer(c_int) function xmi_get_subcomponent_count(n) result(rc) bind(C, name='get_subcomponent_count')
      integer(c_int), intent(out) :: n
      n = 1
      rc = 0
   end function xmi_get_subcomponent_count

   ! ---- component / variable metadata (XMI-only names) -----------------
   integer(c_int) function xmi_get_version(buf) result(rc) bind(C, name='get_version')
      character(kind=c_char), intent(out) :: buf(*)
      call f_to_c_string('SWAP-modern', buf)
      rc = 0
   end function xmi_get_version

   integer(c_int) function xmi_get_var_rank(name, r) result(rc) bind(C, name='get_var_rank')
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int), intent(out) :: r
      associate (unused => name(1)); end associate
      r = 1                 ! all three exchange vars are rank-1
      rc = 0
   end function xmi_get_var_rank

   integer(c_int) function xmi_get_var_shape(name, shp) result(rc) bind(C, name='get_var_shape')
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int), intent(out) :: shp(*)
      associate (unused => name(1)); end associate
      shp(1) = ensemble_ncol()
      rc = 0
   end function xmi_get_var_shape

   ! ---- the zero-copy pointer protocol ---------------------------------
   integer(c_int) function xmi_get_value_ptr(name, ptr) result(rc) bind(C, name='get_value_ptr')
      character(kind=c_char), intent(in) :: name(*)
      type(c_ptr), intent(out) :: ptr
      character(len=64) :: nm
      call c_to_f_string(name, nm)
      rc = 0
      select case (trim(nm))
      case ('gwl');          ptr = c_loc(gwl(1))
      case ('qbot_volume');  ptr = c_loc(qbot_volume(1))
      case ('storage_coef'); ptr = c_loc(storage_coef(1))
      case default
         ptr = c_null_ptr
         rc = 1
         ensemble_last_error = 'unknown var: '//trim(nm)
      end select
   end function xmi_get_value_ptr

   ! ---- error reporting -------------------------------------------------
   integer(c_int) function xmi_get_last_bmi_error(buf) result(rc) bind(C, name='get_last_bmi_error')
      character(kind=c_char), intent(out) :: buf(*)
      call f_to_c_string(trim(ensemble_last_error), buf)
      rc = 0
   end function xmi_get_last_bmi_error

end module swap_xmi_mod
