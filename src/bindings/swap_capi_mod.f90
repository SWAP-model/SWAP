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
   use config_source_mod, only: config_source_t, config_source_memory
   implicit none
   private

   ! Shared singleton — swap_bmi_mod imports these via USE-rename so that
   ! the BMI lifecycle methods (update/finalize) and the CAPI accessors
   ! all operate on the same (state, config) pair.
   ! Multi-instance handles are Phase 3.
   type(swap_state_t),          save, public, target :: capi_state
   type(swap_config_t), target, save, public         :: capi_config

   ! In-memory companion blobs (TOML subfiles + CSV tables) pushed from
   ! Python via swap_attach_config_file, consumed by the in-memory init so
   ! the run reads zero files. Empty (in_memory=.false.) => disk fallback.
   type(config_source_t), save :: capi_companions

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

      if (capi_companions%in_memory) then
         call load_swap_config_from_string(f_text, capi_config, errors, source=capi_companions)
      else
         call load_swap_config_from_string(f_text, capi_config, errors)
      end if
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
   ! In-memory companion files (diskless init)
   !----------------------------------------------------------------------

   !> Drop any previously-attached companion blobs and switch the singleton
   !! into in-memory mode. Call before a batch of swap_attach_config_file.
   function swap_clear_config_files() result(ierr) bind(C, name='swap_clear_config_files')
      integer(c_int) :: ierr
      capi_companions = config_source_memory()
      ierr = 0
   end function swap_clear_config_files

   !> Register one companion file (TOML subfile or CSV table) by name with
   !! its contents, so the next swap_initialize_from_toml_string resolves it
   !! from memory instead of disk. `name` is the filename as referenced in
   !! the config (e.g. "swap.dra.toml", "283.csv").
   function swap_attach_config_file(name, n_name, content, n_content) result(ierr) &
            bind(C, name='swap_attach_config_file')
      character(kind=c_char), intent(in)    :: name(*)
      integer(c_int),  value, intent(in)    :: n_name
      character(kind=c_char), intent(in)    :: content(*)
      integer(c_int),  value, intent(in)    :: n_content
      integer(c_int)                        :: ierr
      character(len=:), allocatable :: f_name, f_content
      integer :: i

      if (.not. capi_companions%in_memory) capi_companions = config_source_memory()

      allocate(character(len=n_name) :: f_name)
      do i = 1, n_name
         f_name(i:i) = name(i)
      end do
      allocate(character(len=n_content) :: f_content)
      do i = 1, n_content
         f_content(i:i) = content(i)
      end do

      call capi_companions%add_blob(f_name, f_content)
      ierr = 0
   end function swap_attach_config_file

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

   !----------------------------------------------------------------------
   ! In-memory results record accessors
   !----------------------------------------------------------------------

   function swap_results_shape(nrows, ncols) result(ierr) bind(C, name='swap_results_shape')
      integer(c_int), intent(out) :: nrows, ncols
      integer(c_int)              :: ierr
      nrows = capi_state%results%nrows
      ncols = capi_state%results%ncols
      ierr  = 0
   end function swap_results_shape

   !> Zero-copy view of the (nrows x ncols) values array. Fortran column-major:
   !! element (i,j) is at flat index (j-1)*nrows + (i-1).
   function swap_view_results(ptr, nrows, ncols) result(ierr) bind(C, name='swap_view_results')
      type(c_ptr),    intent(out) :: ptr
      integer(c_int), intent(out) :: nrows, ncols
      integer(c_int)              :: ierr
      nrows = capi_state%results%nrows
      ncols = capi_state%results%ncols
      if (allocated(capi_state%results%values) .and. nrows > 0 .and. ncols > 0) then
         ptr  = c_loc(capi_state%results%values(1,1))
         ierr = 0
      else
         ptr  = c_null_ptr
         ierr = 1
      end if
   end function swap_view_results

   !> Zero-copy view of the time axis (nrows doubles).
   function swap_view_results_times(ptr, nrows) result(ierr) bind(C, name='swap_view_results_times')
      type(c_ptr),    intent(out) :: ptr
      integer(c_int), intent(out) :: nrows
      integer(c_int)              :: ierr
      nrows = capi_state%results%nrows
      if (allocated(capi_state%results%times) .and. nrows > 0) then
         ptr  = c_loc(capi_state%results%times(1))
         ierr = 0
      else
         ptr  = c_null_ptr
         ierr = 1
      end if
   end function swap_view_results_times

   !> NUL-delimited column names packed into buf (truncated to n bytes).
   function swap_results_columns(buf, n) result(ierr) bind(C, name='swap_results_columns')
      character(kind=c_char), intent(out) :: buf(*)
      integer(c_int),  value, intent(in)  :: n
      integer(c_int)                      :: ierr
      integer :: i, k, c
      k = 0
      do c = 1, capi_state%results%ncols
         do i = 1, len_trim(capi_state%results%col_names(c))
            k = k + 1
            if (k > n) then; buf(min(k,n)) = c_null_char; ierr = 1; return; end if
            buf(k) = capi_state%results%col_names(c)(i:i)
         end do
         k = k + 1
         if (k > n) then; buf(min(k,n)) = c_null_char; ierr = 1; return; end if
         buf(k) = c_null_char
      end do
      ierr = 0
   end function swap_results_columns

   !> Path of the streaming scalar CSV (_output.csv); written when csv output is
   !! enabled and not headless. Lets a huge pyswap run read the file instead of
   !! holding results in RAM.
   function swap_output_filepath(buf, n) result(ierr) bind(C, name='swap_output_filepath')
      character(kind=c_char), intent(out) :: buf(*)
      integer(c_int),  value, intent(in)  :: n
      integer(c_int)                      :: ierr
      character(len=512) :: path
      integer :: i, m
      path = trim(capi_state%timecontrol%pathwork) // &
             trim(capi_state%timecontrol%outfil) // '_output.csv'
      m = len_trim(path)
      do i = 1, min(m, n - 1)
         buf(i) = path(i:i)
      end do
      buf(min(m, n - 1) + 1) = c_null_char
      ierr = 0
   end function swap_output_filepath

end module swap_capi_mod
