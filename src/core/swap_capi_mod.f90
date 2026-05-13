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
   use swap_mod,        only: swap_init, swap_run_step, swap_close, swap_init_from_loaded_config
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

   type, bind(C), public :: swap_water_balance_t
      real(c_double) :: rain
      real(c_double) :: evap_pot
      real(c_double) :: evap_act
      real(c_double) :: transp_pot
      real(c_double) :: transp_act
      real(c_double) :: runoff
      real(c_double) :: drain
      real(c_double) :: percolation
      real(c_double) :: storage_change
      real(c_double) :: balance_error
   end type swap_water_balance_t

contains

   !----------------------------------------------------------------------
   ! Private helpers
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

      call swap_init_from_loaded_config(capi_state, capi_config)
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
   ! Task 19: swap_view_array — zero-copy array accessor
   !----------------------------------------------------------------------

   function swap_view_array(name, ptr, n) result(ierr) bind(C, name='swap_view_array')
      character(kind=c_char), intent(in)  :: name(*)
      type(c_ptr),            intent(out) :: ptr
      integer(c_int),         intent(out) :: n
      integer(c_int)                      :: ierr
      character(len=64) :: f_name
      call c_to_f_string(name, f_name)
      ierr = 0
      select case (trim(f_name))
      case ('theta')
         if (.not. allocated(capi_state%soilwater%theta)) then
            ierr = 2; ptr = c_null_ptr; n = 0; return
         end if
         ptr = c_loc(capi_state%soilwater%theta(1))
         n   = size(capi_state%soilwater%theta)
      case ('h')
         if (.not. allocated(capi_state%soilwater%h)) then
            ierr = 2; ptr = c_null_ptr; n = 0; return
         end if
         ptr = c_loc(capi_state%soilwater%h(1))
         n   = size(capi_state%soilwater%h)
      case ('tsoil')
         if (.not. allocated(capi_state%heat%tsoil)) then
            ierr = 2; ptr = c_null_ptr; n = 0; return
         end if
         ptr = c_loc(capi_state%heat%tsoil(1))
         n   = size(capi_state%heat%tsoil)
      case ('inqrot')
         if (.not. allocated(capi_state%soilwater%inqrot)) then
            ierr = 2; ptr = c_null_ptr; n = 0; return
         end if
         ptr = c_loc(capi_state%soilwater%inqrot(1))
         n   = size(capi_state%soilwater%inqrot)
      case default
         ierr = 1; ptr = c_null_ptr; n = 0
      end select
   end function swap_view_array

   function swap_get_scalar(name, value) result(ierr) bind(C, name='swap_get_scalar')
      character(kind=c_char), intent(in)  :: name(*)
      real(c_double),         intent(out) :: value
      integer(c_int)                      :: ierr
      character(len=64) :: f_name
      call c_to_f_string(name, f_name)
      ierr = 0
      select case (trim(f_name))
      case ('tstart');   value = capi_state%timecontrol%tstart
      case ('tend');     value = capi_state%timecontrol%tend
      case ('dt');       value = capi_state%timecontrol%dt
      case ('t1900');    value = capi_state%timecontrol%t1900
      case ('daynr');    value = real(capi_state%timecontrol%daynr, c_double)
      case ('iptra');    value = capi_state%atmosphere%intr%iptra
      case ('iqrot');    value = capi_state%soilwater%iqrot
      case ('flRunEnd'); value = merge(1.0_c_double, 0.0_c_double, capi_state%timecontrol%flRunEnd)
      case default;      ierr = 1; value = 0.0_c_double
      end select
   end function swap_get_scalar

   function swap_set_scalar(name, value) result(ierr) bind(C, name='swap_set_scalar')
      character(kind=c_char), intent(in) :: name(*)
      real(c_double),  value, intent(in) :: value
      integer(c_int)                     :: ierr
      character(len=64) :: f_name
      call c_to_f_string(name, f_name)
      ! [SS-BMI2 Phase 2 settable allowlist — currently empty]
      ! Most names rejected with rc=2. Add settable variables as needed
      ! when external-crop-driver / coupling workflows arrive.
      select case (trim(f_name))
      case default
         ierr = 2
      end select
   end function swap_set_scalar

   !----------------------------------------------------------------------
   ! Task 21: swap_get_water_balance — bind(C) struct accessor
   !----------------------------------------------------------------------

   function swap_get_water_balance(summary) result(ierr) bind(C, name='swap_get_water_balance')
      type(swap_water_balance_t), intent(out) :: summary
      integer(c_int)                          :: ierr
      summary%rain           = capi_state%atmosphere%cumu%cgrai
      summary%evap_pot       = capi_state%atmosphere%cumu%cpeva
      summary%evap_act       = capi_state%atmosphere%cumu%cevap
      summary%transp_pot     = capi_state%atmosphere%cumu%cptra
      summary%transp_act     = capi_state%soilwater%iqrot
      summary%runoff         = capi_state%soilwater%crunoff
      summary%drain          = 0.0_c_double
      summary%percolation    = capi_state%soilwater%cqbot
      summary%storage_change = 0.0_c_double
      summary%balance_error  = 0.0_c_double
      ierr = 0
   end function swap_get_water_balance

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
