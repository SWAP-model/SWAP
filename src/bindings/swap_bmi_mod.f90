!> @file swap_bmi_mod.f90
!! SS-BMI2: Full CSDMS BMI v2.0 facade for SWAP.
!! Exposes all ~25 CSDMS BMI v2.0 methods as bind(C) procedures.
!! Variable registry: 8 read-only coupling variables + 2 settable BCs.
!! Grid: single 1D vertical soil column (id=0, numnod nodes).
!! Holds a single module-level (state, config) pair — multi-column
!! support is Phase 3 (swap_ensemble_mod).
module swap_bmi_mod
   use iso_c_binding,   only: c_char, c_double, c_int, c_size_t, c_null_char, c_ptr
   use swap_mod,        only: swap_init, swap_run_step, swap_close
   use swap_state_mod,  only: swap_state_t
   use swap_config_mod, only: swap_config_t
   ! Import the shared singleton from swap_capi_mod so that BMI lifecycle
   ! methods (update/finalize) and the CAPI accessors all operate on the
   ! same (state, config) pair.
   use swap_capi_mod,   only: bmi_state => capi_state, bmi_config => capi_config, &
                              bmi_registry => capi_registry
   ! [GR-BH Task 24] numnod/dz removed — now read via bmi_state%mesh%numnod / bmi_state%mesh%dz
   use swap_var_registry_mod, only: build_variable_registry, NS_BMI
   use diagnostics_mod, only: diag_overrides_t, default_embedded_config, &
                              read_logging_overrides_from_file, init_logging
   implicit none
   private

contains

   !----------------------------------------------------------------------
   ! Lifecycle
   !----------------------------------------------------------------------

   function bmi_initialize(config_file, n) result(rc) bind(C, name='initialize')
      character(kind=c_char), intent(in)    :: config_file(*)
      integer(c_int),  value, intent(in)    :: n
      integer(c_int)                        :: rc
      character(len=256) :: f_config_file
      type(diag_overrides_t) :: toml_ov
      call c_to_f_string(config_file, f_config_file)
      call read_logging_overrides_from_file(trim(f_config_file), toml_ov)
      call init_logging(default_embedded_config(), toml_ov)
      call swap_init(trim(f_config_file), bmi_state, bmi_config)
      call build_variable_registry(bmi_registry, bmi_state)
      rc = 0
   end function bmi_initialize

   function bmi_update() result(rc) bind(C, name='update')
      integer(c_int) :: rc
      call swap_run_step(bmi_state, bmi_config)
      if (bmi_state%diag%aborted()) then
         rc = 1
         return
      end if
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

   function bmi_get_start_time(t) result(rc) bind(C, name='get_start_time')
      real(c_double), intent(out) :: t
      integer(c_int)              :: rc
      t = bmi_state%timecontrol%tstart
      rc = 0
   end function bmi_get_start_time

   function bmi_get_end_time(t) result(rc) bind(C, name='get_end_time')
      real(c_double), intent(out) :: t
      integer(c_int)              :: rc
      t = bmi_state%timecontrol%tend
      rc = 0
   end function bmi_get_end_time

   function bmi_get_time_step(dt) result(rc) bind(C, name='get_time_step')
      real(c_double), intent(out) :: dt
      integer(c_int)              :: rc
      dt = bmi_state%timecontrol%dt
      rc = 0
   end function bmi_get_time_step

   function bmi_get_time_units(units_buf, n) result(rc) bind(C, name='get_time_units')
      character(kind=c_char), intent(out) :: units_buf(*)
      integer(c_int),  value, intent(in)  :: n
      integer(c_int)                      :: rc
      call f_to_c_string("d", units_buf, n)
      rc = 0
   end function bmi_get_time_units

   function bmi_update_until(target_time) result(rc) bind(C, name='update_until')
      real(c_double), value, intent(in) :: target_time
      integer(c_int)                    :: rc
      do while (bmi_state%timecontrol%t1900 < target_time .and. &
                .not. bmi_state%timecontrol%flRunEnd)
         call swap_run_step(bmi_state, bmi_config)
      end do
      rc = 0
   end function bmi_update_until

   !----------------------------------------------------------------------
   ! Component info
   !----------------------------------------------------------------------

   function bmi_get_component_name(name_buf, n) result(rc) bind(C, name='get_component_name')
      character(kind=c_char), intent(out) :: name_buf(*)
      integer(c_int),  value, intent(in)  :: n
      integer(c_int)                      :: rc
      call f_to_c_string("SWAP", name_buf, n)
      rc = 0
   end function bmi_get_component_name

   function bmi_get_input_item_count(count) result(rc) bind(C, name='get_input_item_count')
      integer(c_int), intent(out) :: count
      integer(c_int)              :: rc
      count = bmi_registry%count_ns(NS_BMI, want_settable=.true.)
      rc = 0
   end function bmi_get_input_item_count

   function bmi_get_output_item_count(count) result(rc) bind(C, name='get_output_item_count')
      integer(c_int), intent(out) :: count
      integer(c_int)              :: rc
      count = bmi_registry%count_ns(NS_BMI, want_settable=.false.)
      rc = 0
   end function bmi_get_output_item_count

   function bmi_get_input_var_names(buf, n) result(rc) bind(C, name='get_input_var_names')
      ! NUL-delimited list of the BMI settable variable names (from the registry).
      character(kind=c_char), intent(out) :: buf(*)
      integer(c_int),  value, intent(in)  :: n
      integer(c_int)                      :: rc
      character(len=512) :: names
      integer :: i, nb
      call bmi_registry%pack_names(NS_BMI, .true., names, nb)
      do i = 1, min(int(n), nb)
         buf(i) = names(i:i)
      end do
      rc = 0
   end function bmi_get_input_var_names

   function bmi_get_output_var_names(buf, n) result(rc) bind(C, name='get_output_var_names')
      ! NUL-delimited list of the BMI readable variable names (from the registry).
      character(kind=c_char), intent(out) :: buf(*)
      integer(c_int),  value, intent(in)  :: n
      integer(c_int)                      :: rc
      character(len=512) :: names
      integer :: i, nb
      call bmi_registry%pack_names(NS_BMI, .false., names, nb)
      do i = 1, min(int(n), nb)
         buf(i) = names(i:i)
      end do
      rc = 0
   end function bmi_get_output_var_names

   !----------------------------------------------------------------------
   ! Variable accessors
   !----------------------------------------------------------------------

   function bmi_get_value_double(var_name, n, dest) result(rc) bind(C, name='get_value_double')
      character(kind=c_char), intent(in)    :: var_name(*)
      integer(c_int),  value, intent(in)    :: n
      real(c_double),         intent(out)   :: dest(n)
      integer(c_int)                        :: rc
      character(len=64) :: name
      integer           :: idx
      call c_to_f_string(var_name, name)
      idx = bmi_registry%find(NS_BMI, trim(name))
      if (idx == 0 .or. .not. bmi_registry%is_readable(idx)) then
         rc = 1; return
      end if
      call bmi_registry%read_into(idx, dest, int(n))
      rc = 0
   end function bmi_get_value_double

   function bmi_set_value_double(var_name, n, src) result(rc) bind(C, name='set_value_double')
      character(kind=c_char), intent(in) :: var_name(*)
      integer(c_int),  value, intent(in) :: n
      real(c_double),         intent(in) :: src(n)
      integer(c_int)                     :: rc
      character(len=64) :: name
      integer           :: idx
      call c_to_f_string(var_name, name)
      ! Settable BCs (groundwater_level_imposed -> bottom-node head;
      ! bottom_flux_imposed -> qbot) are registry entries flagged settable.
      idx = bmi_registry%find(NS_BMI, trim(name))
      if (idx == 0 .or. .not. bmi_registry%is_settable(idx)) then
         rc = 2; return   ! unknown / not settable
      end if
      call bmi_registry%write_from(idx, src, int(n))
      rc = 0
   end function bmi_set_value_double

   !----------------------------------------------------------------------
   ! Variable metadata
   !----------------------------------------------------------------------

   function bmi_get_var_type(var_name, type_buf, n) result(rc) bind(C, name='get_var_type')
      character(kind=c_char), intent(in)  :: var_name(*)
      character(kind=c_char), intent(out) :: type_buf(*)
      integer(c_int),  value, intent(in)  :: n
      integer(c_int)                      :: rc
      ! All SWAP BMI variables are double precision
      call f_to_c_string("double", type_buf, n)
      rc = 0
   end function bmi_get_var_type

   function bmi_get_var_units(var_name, units_buf, n) result(rc) bind(C, name='get_var_units')
      character(kind=c_char), intent(in)  :: var_name(*)
      character(kind=c_char), intent(out) :: units_buf(*)
      integer(c_int),  value, intent(in)  :: n
      integer(c_int)                      :: rc
      character(len=64) :: name
      integer           :: idx
      call c_to_f_string(var_name, name)
      idx = bmi_registry%find(NS_BMI, trim(name))
      if (idx == 0) then
         call f_to_c_string("", units_buf, n)
         rc = 1
         return
      end if
      call f_to_c_string(trim(bmi_registry%units(idx)), units_buf, n)
      rc = 0
   end function bmi_get_var_units

   function bmi_get_var_grid(var_name, grid_id) result(rc) bind(C, name='get_var_grid')
      character(kind=c_char), intent(in)  :: var_name(*)
      integer(c_int),         intent(out) :: grid_id
      integer(c_int)                      :: rc
      ! All SWAP variables live on the single 1D soil column grid (id=0)
      grid_id = 0
      rc = 0
   end function bmi_get_var_grid

   function bmi_get_var_itemsize(var_name, sz) result(rc) bind(C, name='get_var_itemsize')
      character(kind=c_char), intent(in)  :: var_name(*)
      integer(c_int),         intent(out) :: sz
      integer(c_int)                      :: rc
      ! All variables are 8-byte double precision
      sz = 8
      rc = 0
   end function bmi_get_var_itemsize

   function bmi_get_var_nbytes(var_name, nb) result(rc) bind(C, name='get_var_nbytes')
      character(kind=c_char), intent(in)  :: var_name(*)
      integer(c_int),         intent(out) :: nb
      integer(c_int)                      :: rc
      character(len=64) :: name
      integer           :: idx
      call c_to_f_string(var_name, name)
      idx = bmi_registry%find(NS_BMI, trim(name))
      if (idx == 0) then
         nb = 0
         rc = 1
         return
      end if
      ! Profile arrays: numnod × 8 bytes; scalars: 8 bytes (from the registry).
      nb = bmi_registry%nbytes(idx)
      rc = 0
   end function bmi_get_var_nbytes

   function bmi_get_var_location(var_name, loc_buf, n) result(rc) bind(C, name='get_var_location')
      character(kind=c_char), intent(in)  :: var_name(*)
      character(kind=c_char), intent(out) :: loc_buf(*)
      integer(c_int),  value, intent(in)  :: n
      integer(c_int)                      :: rc
      ! All SWAP variables are defined at grid nodes
      call f_to_c_string("node", loc_buf, n)
      rc = 0
   end function bmi_get_var_location

   !----------------------------------------------------------------------
   ! Grid metadata — single 1D vertical soil column (grid id = 0)
   !----------------------------------------------------------------------

   function bmi_get_grid_type(grid_id, type_buf, n) result(rc) bind(C, name='get_grid_type')
      integer(c_int),  value, intent(in)  :: grid_id
      character(kind=c_char), intent(out) :: type_buf(*)
      integer(c_int),  value, intent(in)  :: n
      integer(c_int)                      :: rc
      call f_to_c_string("uniform_rectilinear", type_buf, n)
      rc = 0
   end function bmi_get_grid_type

   function bmi_get_grid_rank(grid_id, rank) result(rc) bind(C, name='get_grid_rank')
      integer(c_int), value,  intent(in)  :: grid_id
      integer(c_int),         intent(out) :: rank
      integer(c_int)                      :: rc
      rank = 1   ! 1D vertical column
      rc = 0
   end function bmi_get_grid_rank

   function bmi_get_grid_size(grid_id, sz) result(rc) bind(C, name='get_grid_size')
      integer(c_int), value,  intent(in)  :: grid_id
      integer(c_int),         intent(out) :: sz
      integer(c_int)                      :: rc
      sz = bmi_state%mesh%numnod
      rc = 0
   end function bmi_get_grid_size

   function bmi_get_grid_shape(grid_id, shape_arr, max_n) result(rc) bind(C, name='get_grid_shape')
      integer(c_int), value,  intent(in)  :: grid_id
      integer(c_int),         intent(out) :: shape_arr(*)
      integer(c_int), value,  intent(in)  :: max_n
      integer(c_int)                      :: rc
      if (max_n >= 1) shape_arr(1) = bmi_state%mesh%numnod
      rc = 0
   end function bmi_get_grid_shape

   function bmi_get_grid_node_count(grid_id, count) result(rc) bind(C, name='get_grid_node_count')
      integer(c_int), value,  intent(in)  :: grid_id
      integer(c_int),         intent(out) :: count
      integer(c_int)                      :: rc
      count = bmi_state%mesh%numnod
      rc = 0
   end function bmi_get_grid_node_count

   function bmi_get_grid_z(grid_id, z_arr, max_n) result(rc) bind(C, name='get_grid_z')
      ! Returns cumulative depth (negative-downward convention, in cm).
      integer(c_int), value,  intent(in)  :: grid_id
      real(c_double),         intent(out) :: z_arr(*)
      integer(c_int), value,  intent(in)  :: max_n
      integer(c_int)                      :: rc
      integer :: i, m
      real(c_double) :: cumz
      m = min(max_n, bmi_state%mesh%numnod)
      cumz = 0.0_c_double
      do i = 1, m
         cumz = cumz - bmi_state%mesh%dz(i)   ! negative-downward
         z_arr(i) = cumz
      end do
      rc = 0
   end function bmi_get_grid_z

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

   subroutine f_to_c_string(f_str, c_buf, max_len)
      ! Copies a Fortran character string into a C NUL-terminated char array.
      ! Writes at most max_len-1 characters, then appends a NUL terminator.
      character(len=*),              intent(in)  :: f_str
      character(kind=c_char),        intent(out) :: c_buf(*)
      integer,                       intent(in)  :: max_len
      integer :: i, n
      n = min(len_trim(f_str), max_len - 1)
      do i = 1, n
         c_buf(i) = f_str(i:i)
      end do
      c_buf(n + 1) = c_null_char
   end subroutine f_to_c_string

end module swap_bmi_mod
