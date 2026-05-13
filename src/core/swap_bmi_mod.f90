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
   use variables,       only: numnod, dz
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
      count = 2
      rc = 0
   end function bmi_get_input_item_count

   function bmi_get_output_item_count(count) result(rc) bind(C, name='get_output_item_count')
      integer(c_int), intent(out) :: count
      integer(c_int)              :: rc
      count = 8
      rc = 0
   end function bmi_get_output_item_count

   function bmi_get_input_var_names(buf, n) result(rc) bind(C, name='get_input_var_names')
      ! NUL-delimited list of 2 settable variable names packed into buf.
      ! Format: "groundwater_level_imposed\0bottom_flux_imposed\0"
      character(kind=c_char), intent(out) :: buf(*)
      integer(c_int),  value, intent(in)  :: n
      integer(c_int)                      :: rc
      integer :: i
      character(len=46), parameter :: names = &
         "groundwater_level_imposed" // c_null_char // "bottom_flux_imposed" // c_null_char
      do i = 1, min(n, len(names))
         buf(i) = names(i:i)
      end do
      rc = 0
   end function bmi_get_input_var_names

   function bmi_get_output_var_names(buf, n) result(rc) bind(C, name='get_output_var_names')
      ! NUL-delimited list of 8 readable variable names packed into buf.
      character(kind=c_char), intent(out) :: buf(*)
      integer(c_int),  value, intent(in)  :: n
      integer(c_int)                      :: rc
      integer :: i
      character(len=145), parameter :: names = &
         "soil_water_content"         // c_null_char // &
         "pressure_head"              // c_null_char // &
         "soil_temperature"           // c_null_char // &
         "groundwater_level"          // c_null_char // &
         "bottom_flux"                // c_null_char // &
         "actual_evapotranspiration"  // c_null_char // &
         "recharge"                   // c_null_char // &
         "surface_runoff"             // c_null_char
      do i = 1, min(n, len(names))
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
      call c_to_f_string(var_name, name)
      rc = 0
      select case (trim(name))
      case ('soil_water_content')
         call f_to_c_string("m3 m-3", units_buf, n)
      case ('pressure_head')
         call f_to_c_string("cm", units_buf, n)
      case ('soil_temperature')
         call f_to_c_string("degC", units_buf, n)
      case ('groundwater_level', 'groundwater_level_imposed')
         call f_to_c_string("cm", units_buf, n)
      case ('bottom_flux', 'bottom_flux_imposed', 'recharge', &
            'actual_evapotranspiration', 'surface_runoff')
         call f_to_c_string("cm d-1", units_buf, n)
      case default
         call f_to_c_string("", units_buf, n)
         rc = 1
      end select
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
      call c_to_f_string(var_name, name)
      rc = 0
      select case (trim(name))
      case ('soil_water_content', 'pressure_head', 'soil_temperature')
         ! Profile arrays: numnod elements × 8 bytes
         nb = numnod * 8
      case ('groundwater_level', 'bottom_flux', 'actual_evapotranspiration', &
            'recharge', 'surface_runoff', &
            'groundwater_level_imposed', 'bottom_flux_imposed')
         ! Scalar variables: 1 element × 8 bytes
         nb = 8
      case default
         nb = 0
         rc = 1
      end select
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
      sz = numnod
      rc = 0
   end function bmi_get_grid_size

   function bmi_get_grid_shape(grid_id, shape_arr, max_n) result(rc) bind(C, name='get_grid_shape')
      integer(c_int), value,  intent(in)  :: grid_id
      integer(c_int),         intent(out) :: shape_arr(*)
      integer(c_int), value,  intent(in)  :: max_n
      integer(c_int)                      :: rc
      if (max_n >= 1) shape_arr(1) = numnod
      rc = 0
   end function bmi_get_grid_shape

   function bmi_get_grid_node_count(grid_id, count) result(rc) bind(C, name='get_grid_node_count')
      integer(c_int), value,  intent(in)  :: grid_id
      integer(c_int),         intent(out) :: count
      integer(c_int)                      :: rc
      count = numnod
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
      m = min(max_n, numnod)
      cumz = 0.0_c_double
      do i = 1, m
         cumz = cumz - dz(i)   ! negative-downward
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
