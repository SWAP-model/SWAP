!> Reader for the [surface_water] section of swap.toml.
!!
!! Populates a surface_water_config_t from a fully-formed swap.toml
!! document. Top-level scalars cover the surface-water switches and
!! the period count `nmper`. Per-period management arrays live under
!! [surface_water.management]; per-period weir arrays live under
!! [surface_water.weir].
!!
!! Mechanical pattern (mirrors read_bottom_boundary_toml): scalars
!! one per line; per-period arrays sized by `nmper` (read first, then
!! allocate, then loop). No switch-conditional skipping — the validator
!! handles cross-field consistency. If [surface_water] is absent the
!! reader returns silently and leaves `config` at defaults.
!!
!! Date strings in `impend` decode through parse_date_to_days1900,
!! matching the simulation reader's convention. Phase 4f-prep Task B2.
module read_surface_water_toml_mod
   use iso_fortran_env, only: real64
   use tomlf, only: toml_table, toml_array, toml_datetime, get_value, len
   use surface_water_config_mod, only: surface_water_config_t
   use toml_field_helpers_mod, only: get_table,                         &
                                     get_optional_int_with_default,     &
                                     get_optional_real_with_default,    &
                                     parse_date_to_days1900
   use error_mod, only: error_collection_t, ERR_PARSE_TYPE_MISMATCH
   implicit none
   private

   public :: read_surface_water_toml

contains

   subroutine read_surface_water_toml(doc_root, config, errors)
      type(toml_table), pointer,     intent(in)    :: doc_root
      type(surface_water_config_t),  intent(inout) :: config
      type(error_collection_t),      intent(inout) :: errors

      type(toml_table), pointer :: sec, mgmt, weir

      call get_table(doc_root, 'surface_water', sec, 'surface_water', errors)
      if (.not. associated(sec)) return

      ! Top-level switches and scalars.
      call get_optional_int_with_default (sec, 'swsrf',  config%swsrf,  1, &
                                          'surface_water.swsrf',  errors)
      call get_optional_int_with_default (sec, 'swsec',  config%swsec,  1, &
                                          'surface_water.swsec',  errors)
      call get_optional_real_with_default(sec, 'wlact',  config%wlact,  0.0_real64, &
                                          'surface_water.wlact',  errors)
      call get_optional_real_with_default(sec, 'osswlm', config%osswlm, 0.0_real64, &
                                          'surface_water.osswlm', errors)
      call get_optional_int_with_default (sec, 'nmper',  config%nmper,  0, &
                                          'surface_water.nmper',  errors)
      call get_optional_int_with_default (sec, 'swqhr',  config%swqhr,  1, &
                                          'surface_water.swqhr',  errors)
      call get_optional_real_with_default(sec, 'sofcu',  config%sofcu,  0.0_real64, &
                                          'surface_water.sofcu',  errors)

      ! Per-period management arrays.
      call get_table(sec, 'management', mgmt, 'surface_water.management', errors)
      if (associated(mgmt)) then
         call read_period_array_real(mgmt, 'impend', config%nmper, &
                                     config%impend, &
                                     'surface_water.management.impend', &
                                     errors, decode_date=.true.)
         call read_period_array_int (mgmt, 'swman',  config%nmper, &
                                     config%swman,  &
                                     'surface_water.management.swman',  errors)
         call read_period_array_real(mgmt, 'wscap',  config%nmper, &
                                     config%wscap,  &
                                     'surface_water.management.wscap',  errors)
         call read_period_array_real(mgmt, 'wldip',  config%nmper, &
                                     config%wldip,  &
                                     'surface_water.management.wldip',  errors)
         call read_period_array_int (mgmt, 'intwl',  config%nmper, &
                                     config%intwl,  &
                                     'surface_water.management.intwl',  errors)
      end if

      ! Per-period weir arrays (SWQHR=1 exponential discharge relation).
      call get_table(sec, 'weir', weir, 'surface_water.weir', errors)
      if (associated(weir)) then
         call read_period_array_real(weir, 'hbweir', config%nmper, &
                                     config%hbweir, &
                                     'surface_water.weir.hbweir', errors)
         call read_period_array_real(weir, 'alphaw', config%nmper, &
                                     config%alphaw, &
                                     'surface_water.weir.alphaw', errors)
         call read_period_array_real(weir, 'betaw',  config%nmper, &
                                     config%betaw,  &
                                     'surface_water.weir.betaw',  errors)
      end if
   end subroutine read_surface_water_toml

   !> Decode a flat 1-D real TOML array sized to nmper. Optionally each
   !! element is a TOML date literal which we decode through
   !! parse_date_to_days1900 (used for impend). If the array length
   !! differs from nmper, append a parse-type-mismatch error and skip
   !! allocation (validator will see unallocated and complain on its own).
   !! Absent key leaves arr unallocated.
   subroutine read_period_array_real(sec, key, nmper, arr, context, errors, decode_date)
      type(toml_table), pointer, intent(in)    :: sec
      character(len=*),          intent(in)    :: key
      integer,                   intent(in)    :: nmper
      real(real64), allocatable, intent(out)   :: arr(:)
      character(len=*),          intent(in)    :: context
      type(error_collection_t),  intent(inout) :: errors
      logical, optional,         intent(in)    :: decode_date

      type(toml_array), pointer :: outer
      integer :: n, i, stat
      real(real64) :: val
      type(toml_datetime) :: dtv
      logical :: as_date
      character(len=128) :: msg

      if (.not. associated(sec)) return

      as_date = .false.
      if (present(decode_date)) as_date = decode_date

      outer => null()
      call get_value(sec, key, outer, requested=.false., stat=stat)
      if (.not. associated(outer)) return

      n = len(outer)
      if (n /= nmper) then
         write(msg, '("array length=",I0," /= nmper=",I0)') n, nmper
         call errors%append(ERR_PARSE_TYPE_MISMATCH, trim(msg), context)
         return
      end if

      allocate(arr(n))
      arr = 0.0_real64

      do i = 1, n
         if (as_date) then
            call get_value(outer, i, dtv, stat=stat)
            if (stat /= 0) then
               call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                                  "expected date cell", context)
               if (allocated(arr)) deallocate(arr)
               return
            end if
            arr(i) = parse_date_to_days1900(dtv)
         else
            call get_value(outer, i, val, stat=stat)
            if (stat /= 0) then
               call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                                  "non-real cell", context)
               if (allocated(arr)) deallocate(arr)
               return
            end if
            arr(i) = val
         end if
      end do
   end subroutine read_period_array_real

   !> Integer counterpart to read_period_array_real. Length-mismatch
   !! against nmper appends ERR_PARSE_TYPE_MISMATCH and skips allocation.
   subroutine read_period_array_int(sec, key, nmper, arr, context, errors)
      type(toml_table), pointer, intent(in)    :: sec
      character(len=*),          intent(in)    :: key
      integer,                   intent(in)    :: nmper
      integer, allocatable,      intent(out)   :: arr(:)
      character(len=*),          intent(in)    :: context
      type(error_collection_t),  intent(inout) :: errors

      type(toml_array), pointer :: outer
      integer :: n, i, stat, val
      character(len=128) :: msg

      if (.not. associated(sec)) return

      outer => null()
      call get_value(sec, key, outer, requested=.false., stat=stat)
      if (.not. associated(outer)) return

      n = len(outer)
      if (n /= nmper) then
         write(msg, '("array length=",I0," /= nmper=",I0)') n, nmper
         call errors%append(ERR_PARSE_TYPE_MISMATCH, trim(msg), context)
         return
      end if

      allocate(arr(n))
      arr = 0

      do i = 1, n
         call get_value(outer, i, val, stat=stat)
         if (stat /= 0) then
            call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                               "non-integer cell", context)
            if (allocated(arr)) deallocate(arr)
            return
         end if
         arr(i) = val
      end do
   end subroutine read_period_array_int

end module read_surface_water_toml_mod
