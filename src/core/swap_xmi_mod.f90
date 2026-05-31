!> XMI (BMI + Deltares eXtended Model Interface) C-ABI facade over the SWAP
!! ensemble. Symbol names/signatures match what xmipy's XmiWrapper /
!! imod_coupler's SwapWrapper call.
!!
!! This module is the SOLE provider of the coupling C-ABI for the
!! libswap_xmi.so shared library. It deliberately re-implements the BMI
!! lifecycle / time / component-info / variable-metadata C names
!! (initialize/update/finalize/get_current_time/...) so they are backed by
!! the ENSEMBLE (swap_ensemble_mod), not the single capi_state singleton that
!! swap_bmi_mod / swap_capi_mod drive. Because both modules would otherwise
!! export the same bind(C) names, swap_xmi_mod and swap_bmi_mod/swap_capi_mod
!! are kept in SEPARATE shared libraries (libswap_xmi.so vs libswap_bmi.so);
!! see meson.build. The coupler loads libswap_xmi.so.
!!
!! Scope note: the module-global error state (error_mod library_mode) and the
!! module-global ensemble (swap_ensemble_mod) are single-kernel-per-process.
!! That matches imod_coupler's one-kernel-per-library model and is the
!! documented scope for this demo; true multi-instance support would require
!! per-instance error + ensemble state (handle-based, a future arc).
module swap_xmi_mod
   use iso_c_binding
   use swap_ensemble_mod
   use error_mod, only: library_fatal_raised
   implicit none
   private

   ! Last-error text buffer (filled when a function returns rc /= 0).
   character(len=1024), save :: last_error = ''

contains

   ! ---- lifecycle -------------------------------------------------------
   integer(c_int) function xmi_initialize(config_file) result(rc) bind(C, name='initialize')
      character(kind=c_char), intent(in) :: config_file(*)
      character(len=512) :: f_path
      integer :: ncol
      call c_to_f_string(config_file, f_path)
      ! f_path points at the SWAP working dir's coupled config; ncol comes
      ! from the sidecar ensemble descriptor 'ensemble.txt' (one int) next to it.
      ncol = read_ncol_sidecar(trim(f_path))
      if (ncol <= 0) then
         last_error = 'ensemble.txt missing/invalid'
         rc = 1
         return
      end if
      rc = ensemble_init(trim(f_path), ncol)
      if (rc /= 0) last_error = 'ensemble_init failed'
   end function xmi_initialize

   integer(c_int) function xmi_update() result(rc) bind(C, name='update')
      rc = ensemble_step_day()
      if (rc /= 0) last_error = 'ensemble_step_day failed'
   end function xmi_update

   integer(c_int) function xmi_update_until(t) result(rc) bind(C, name='update_until')
      real(c_double), value, intent(in) :: t
      rc = ensemble_step_day()   ! MF6 leads; one day per call in coupled mode
      if (rc /= 0) last_error = 'ensemble_step_day failed'
   end function xmi_update_until

   integer(c_int) function xmi_finalize() result(rc) bind(C, name='finalize')
      rc = ensemble_finalize()
   end function xmi_finalize

   ! ---- time ------------------------------------------------------------
   integer(c_int) function xmi_get_current_time(t) result(rc) bind(C, name='get_current_time')
      real(c_double), intent(out) :: t
      t = ensemble_current_t1900()
      rc = 0
   end function xmi_get_current_time

   integer(c_int) function xmi_get_start_time(t) result(rc) bind(C, name='get_start_time')
      real(c_double), intent(out) :: t
      t = ensemble_start_t1900()
      rc = 0
   end function xmi_get_start_time

   integer(c_int) function xmi_get_end_time(t) result(rc) bind(C, name='get_end_time')
      real(c_double), intent(out) :: t
      t = ensemble_end_t1900()
      rc = 0
   end function xmi_get_end_time

   integer(c_int) function xmi_get_time_step(dt) result(rc) bind(C, name='get_time_step')
      real(c_double), intent(out) :: dt
      dt = 1.0_c_double   ! daily coupling
      rc = 0
   end function xmi_get_time_step

   ! ---- component info --------------------------------------------------
   integer(c_int) function xmi_get_version(buf) result(rc) bind(C, name='get_version')
      character(kind=c_char), intent(out) :: buf(*)
      call f_to_c_string('SWAP-modern', buf)
      rc = 0
   end function xmi_get_version

   integer(c_int) function xmi_get_component_name(buf) result(rc) bind(C, name='get_component_name')
      character(kind=c_char), intent(out) :: buf(*)
      call f_to_c_string('SWAP', buf)
      rc = 0
   end function xmi_get_component_name

   integer(c_int) function xmi_get_input_item_count(n) result(rc) bind(C, name='get_input_item_count')
      integer(c_int), intent(out) :: n
      n = 1   ! gwl
      rc = 0
   end function xmi_get_input_item_count

   integer(c_int) function xmi_get_output_item_count(n) result(rc) bind(C, name='get_output_item_count')
      integer(c_int), intent(out) :: n
      n = 2   ! qbot_volume, storage_coef
      rc = 0
   end function xmi_get_output_item_count

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
      if (rc /= 0) last_error = 'ensemble_step_day failed'
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

   ! ---- variable metadata ----------------------------------------------
   integer(c_int) function xmi_get_var_rank(name, r) result(rc) bind(C, name='get_var_rank')
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int), intent(out) :: r
      associate (unused => name(1)); end associate
      r = 1                 ! all three exchange vars are rank-1
      rc = 0
   end function xmi_get_var_rank

   integer(c_int) function xmi_get_var_type(name, tbuf) result(rc) bind(C, name='get_var_type')
      character(kind=c_char), intent(in) :: name(*)
      character(kind=c_char), intent(out) :: tbuf(*)
      associate (unused => name(1)); end associate
      call f_to_c_string('double', tbuf)   ! string must start 'double'
      rc = 0
   end function xmi_get_var_type

   integer(c_int) function xmi_get_var_shape(name, shp) result(rc) bind(C, name='get_var_shape')
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int), intent(out) :: shp(*)
      associate (unused => name(1)); end associate
      shp(1) = ensemble_ncol()
      rc = 0
   end function xmi_get_var_shape

   integer(c_int) function xmi_get_var_itemsize(name, sz) result(rc) bind(C, name='get_var_itemsize')
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int), intent(out) :: sz
      associate (unused => name(1)); end associate
      sz = 8
      rc = 0
   end function xmi_get_var_itemsize

   integer(c_int) function xmi_get_var_nbytes(name, nb) result(rc) bind(C, name='get_var_nbytes')
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int), intent(out) :: nb
      associate (unused => name(1)); end associate
      nb = 8*ensemble_ncol()
      rc = 0
   end function xmi_get_var_nbytes

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
         last_error = 'unknown var: '//trim(nm)
      end select
   end function xmi_get_value_ptr

   ! ---- error reporting -------------------------------------------------
   integer(c_int) function xmi_get_last_bmi_error(buf) result(rc) bind(C, name='get_last_bmi_error')
      character(kind=c_char), intent(out) :: buf(*)
      call f_to_c_string(trim(last_error), buf)
      rc = 0
   end function xmi_get_last_bmi_error

   ! ---- private C-string helpers (mirror swap_bmi_mod.f90) --------------
   subroutine c_to_f_string(c_str, f_str)
      character(kind=c_char), intent(in) :: c_str(*)
      character(len=*), intent(out) :: f_str
      integer :: i
      f_str = ''
      do i = 1, len(f_str)
         if (c_str(i) == c_null_char) exit
         f_str(i:i) = c_str(i)
      end do
   end subroutine c_to_f_string

   subroutine f_to_c_string(f_str, c_str)
      character(len=*), intent(in) :: f_str
      character(kind=c_char), intent(out) :: c_str(*)
      integer :: i
      do i = 1, len_trim(f_str)
         c_str(i) = f_str(i:i)
      end do
      c_str(len_trim(f_str) + 1) = c_null_char
   end subroutine f_to_c_string

   integer function read_ncol_sidecar(dir_or_file) result(n)
      character(len=*), intent(in) :: dir_or_file
      integer :: u, ios
      character(len=512) :: path
      n = 0
      ! ensemble.txt sits next to the config; accept a dir or a file path.
      path = trim(dir_or_file)
      if (index(path, '.toml') > 0) path = path(1:scan(path, '/', back=.true.))
      open (newunit=u, file=trim(path)//'ensemble.txt', status='old', action='read', iostat=ios)
      if (ios /= 0) return
      read (u, *, iostat=ios) n
      close (u)
      if (ios /= 0) n = 0
   end function read_ncol_sidecar

end module swap_xmi_mod
