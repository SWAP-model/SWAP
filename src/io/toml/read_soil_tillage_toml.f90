!> @file read_soil_tillage_toml.f90
!! Parses the optional [soil.tillage] block into soil_config_t%tillage.
!! Block is omitted when soil.swtill = 0; reader leaves defaults intact
!! in that case. Sibling to read_soil_toml.f90 to keep that file focused.
module read_soil_tillage_toml_mod

   use, intrinsic :: iso_fortran_env, only: real64
   use tomlf, only: toml_table, toml_array, toml_datetime, get_value, len
   use toml_field_helpers_mod, only: get_array_of_tables,                 &
                                     get_optional_int_with_default,       &
                                     get_optional_real_with_default
   use soil_config_mod, only: soil_tillage_t, soil_tillage_event_t,       &
                              soil_tillage_type_t
   use error_mod, only: error_collection_t, ERR_PARSE_TYPE_MISMATCH
   implicit none
   private

   public :: read_soil_tillage_toml

contains

   !> Populate `tillage` from the optional [soil.tillage] table on `soil_sec`.
   !! Missing block is benign — leaves defaults.
   subroutine read_soil_tillage_toml(soil_sec, tillage, errors)
      type(toml_table), pointer, intent(in)    :: soil_sec
      type(soil_tillage_t),      intent(inout) :: tillage
      type(error_collection_t),  intent(inout) :: errors

      type(toml_table), pointer :: till_tbl, item
      type(toml_array), pointer :: events_arr, types_arr
      integer :: i, n, stat

      if (.not. associated(soil_sec)) return

      ! [soil.tillage] is optional.
      till_tbl => null()
      call get_value(soil_sec, 'tillage', till_tbl, requested=.false., stat=stat)
      if (stat /= 0 .or. .not. associated(till_tbl)) return

      call get_optional_int_with_default(till_tbl, 'i_n_model', tillage%i_n_model, 2, &
                                         'soil.tillage.i_n_model', errors)
      call get_optional_int_with_default(till_tbl, 'iRedist',   tillage%iRedist,   2, &
                                         'soil.tillage.iRedist',   errors)

      ! [[soil.tillage.events]]
      events_arr => null()
      call get_array_of_tables(till_tbl, 'events', events_arr, &
                               'soil.tillage.events', errors)
      if (associated(events_arr)) then
         n = len(events_arr)
         if (n > 0) then
            allocate(tillage%events(n))
            do i = 1, n
               item => null()
               call get_value(events_arr, i, item, stat=stat)
               if (stat /= 0 .or. .not. associated(item)) cycle
               call read_event_row(item, tillage%events(i), errors)
            end do
         else
            allocate(tillage%events(0))
         end if
      end if

      ! [[soil.tillage.types]]
      types_arr => null()
      call get_array_of_tables(till_tbl, 'types', types_arr, &
                               'soil.tillage.types', errors)
      if (associated(types_arr)) then
         n = len(types_arr)
         if (n > 0) then
            allocate(tillage%types(n))
            do i = 1, n
               item => null()
               call get_value(types_arr, i, item, stat=stat)
               if (stat /= 0 .or. .not. associated(item)) cycle
               call read_type_row(item, tillage%types(i), errors)
            end do
         else
            allocate(tillage%types(0))
         end if
      end if
   end subroutine read_soil_tillage_toml


   subroutine read_event_row(row, ev, errors)
      type(toml_table), pointer, intent(in)    :: row
      type(soil_tillage_event_t), intent(out)  :: ev
      type(error_collection_t),  intent(inout) :: errors
      type(toml_datetime) :: dtv
      integer             :: stat

      ! Date as TOML local-date literal — reader stores as ISO string for
      ! round-trip fidelity; adapter parses to days-since-1900.
      call get_value(row, 'date', dtv, stat=stat)
      if (stat /= 0) then
         call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                            "missing or non-date 'date' key", &
                            'soil.tillage.events')
      else
         write(ev%date, '(i4.4,"-",i2.2,"-",i2.2)') &
            dtv%date%year, dtv%date%month, dtv%date%day
      end if

      call get_optional_real_with_default(row, 'z',         ev%z,         0.0_real64, &
                                          'soil.tillage.events.z',         errors)
      call get_optional_real_with_default(row, 'intensity', ev%intensity, 0.0_real64, &
                                          'soil.tillage.events.intensity', errors)
      call get_optional_int_with_default(row,  'type_id',   ev%type_id,   0, &
                                         'soil.tillage.events.type_id',   errors)
   end subroutine read_event_row


   subroutine read_type_row(row, ty, errors)
      type(toml_table), pointer, intent(in)    :: row
      type(soil_tillage_type_t), intent(out)   :: ty
      type(error_collection_t),  intent(inout) :: errors
      call get_optional_int_with_default(row,  'id',          ty%id,          0, &
                                         'soil.tillage.types.id',          errors)
      call get_optional_real_with_default(row, 'rho_cons',    ty%rho_cons,    0.0_real64, &
                                          'soil.tillage.types.rho_cons',    errors)
      call get_optional_real_with_default(row, 'rho_tillage', ty%rho_tillage, 0.0_real64, &
                                          'soil.tillage.types.rho_tillage', errors)
      call get_optional_real_with_default(row, 'k_R',         ty%k_R,         0.0_real64, &
                                          'soil.tillage.types.k_R',         errors)
      call get_optional_real_with_default(row, 'rho_match',   ty%rho_match,   -99.0_real64, &
                                          'soil.tillage.types.rho_match',   errors)
      call get_optional_real_with_default(row, 'N_match',     ty%N_match,     -99.0_real64, &
                                          'soil.tillage.types.N_match',     errors)
   end subroutine read_type_row

end module read_soil_tillage_toml_mod
