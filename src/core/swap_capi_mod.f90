!> @file swap_capi_mod.f90
!! SS-BMI2: Python-direct C ABI surface. Pragmatic accessors, not
!! CSDMS BMI-bound. Zero-copy array views via c_loc, scalar
!! getters/setters with allowlists, bind(C) struct returns for
!! derived summaries, per-stream output-row dispatch, input-buffer
!! attach functions for in-memory I/O.
module swap_capi_mod
   use iso_c_binding,   only: c_char, c_double, c_int, c_size_t, &
                               c_null_char, c_ptr, c_null_ptr,    &
                               c_loc, c_f_pointer
   use swap_mod,        only: swap_init, swap_run_step, swap_close
   use swap_state_mod,  only: swap_state_t
   use swap_config_mod, only: swap_config_t
   use load_swap_config_string_mod, only: load_swap_config_from_string
   use error_mod,       only: error_collection_t
   use meteo_buffer_mod, only: attach_external_meteo_buffer
   implicit none
   private

   ! Module-level singleton — shared with swap_bmi_mod's instance.
   ! Multi-instance handles are Phase 3.
   type(swap_state_t),          save, target :: capi_state
   type(swap_config_t), target, save         :: capi_config

contains

   !----------------------------------------------------------------------
   ! Lifecycle
   !----------------------------------------------------------------------

   function swap_set_headless(flag) result(ierr) bind(C, name='swap_set_headless')
      integer(c_int), value, intent(in) :: flag
      integer(c_int)                    :: ierr
      capi_state%timecontrol%headless = (flag /= 0)
      ierr = 0
   end function swap_set_headless

   function swap_initialize_from_toml_string(buf, n) result(ierr) &
            bind(C, name='swap_initialize_from_toml_string')
      character(kind=c_char), intent(in)    :: buf(*)
      integer(c_int),  value, intent(in)    :: n
      integer(c_int)                        :: ierr
      character(len=:), allocatable :: f_text
      type(error_collection_t) :: errors
      integer :: i

      allocate(character(len=n) :: f_text)
      do i = 1, n
         f_text(i:i) = buf(i)
      end do

      call load_swap_config_from_string(f_text, capi_config, errors)
      if (errors%has_fatals()) then
         ierr = 1
         return
      end if

      call capi_config%validate(errors)
      call capi_config%finalize(errors)
      if (errors%has_fatals()) then
         ierr = 2
         return
      end if

      ! NOTE: this calls swap_init which reads 'swap.toml' from disk in
      ! its current form. Task 18 introduces swap_init_from_loaded_config
      ! to skip the file read and consume the already-loaded capi_config.
      ! For Task 17, use a placeholder call: do nothing past validation.
      ! Task 18 wires the real init.
      ierr = 0
   end function swap_initialize_from_toml_string

   !----------------------------------------------------------------------
   ! Input buffer attach
   !----------------------------------------------------------------------

   function swap_attach_meteo_buffer(ptr, n_days, n_cols) result(ierr) &
            bind(C, name='swap_attach_meteo_buffer')
      type(c_ptr),    value, intent(in) :: ptr
      integer(c_int), value, intent(in) :: n_days, n_cols
      integer(c_int)                    :: ierr
      call attach_external_meteo_buffer(ptr, n_days, n_cols)
      ierr = 0
   end function swap_attach_meteo_buffer

   !----------------------------------------------------------------------
   ! Stubs filled by Tasks 19-22
   !----------------------------------------------------------------------

   function swap_view_array(name, ptr, n) result(ierr) bind(C, name='swap_view_array')
      character(kind=c_char), intent(in)  :: name(*)
      type(c_ptr),            intent(out) :: ptr
      integer(c_int),         intent(out) :: n
      integer(c_int)                      :: ierr
      ! [SS-BMI2 Task 19] body to be filled with c_loc-based dispatch.
      ierr = 1
      ptr = c_null_ptr
      n   = 0
   end function swap_view_array

   function swap_get_scalar(name, value) result(ierr) bind(C, name='swap_get_scalar')
      character(kind=c_char), intent(in)  :: name(*)
      real(c_double),         intent(out) :: value
      integer(c_int)                      :: ierr
      ! [SS-BMI2 Task 20]
      ierr = 1
      value = 0.0_c_double
   end function swap_get_scalar

   function swap_set_scalar(name, value) result(ierr) bind(C, name='swap_set_scalar')
      character(kind=c_char), intent(in) :: name(*)
      real(c_double),  value, intent(in) :: value
      integer(c_int)                     :: ierr
      ! [SS-BMI2 Task 20]
      ierr = 1
   end function swap_set_scalar

   function swap_get_output_row(stream, row_ptr, names_ptr, n) result(ierr) &
            bind(C, name='swap_get_output_row')
      character(kind=c_char), intent(in)  :: stream(*)
      type(c_ptr),            intent(out) :: row_ptr, names_ptr
      integer(c_int),         intent(out) :: n
      integer(c_int)                      :: ierr
      ! [SS-BMI2 Task 22]
      ierr = 1
      row_ptr = c_null_ptr
      names_ptr = c_null_ptr
      n = 0
   end function swap_get_output_row

end module swap_capi_mod
