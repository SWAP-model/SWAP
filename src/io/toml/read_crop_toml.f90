!> Reader for the [crop] section of a SWAP TOML.
module read_crop_toml_mod
   use iso_fortran_env, only: real64
   use tomlf, only: toml_table, toml_array, toml_datetime, get_value, len
   use crop_config_mod, only: crop_config_t
   use toml_field_helpers_mod, only: get_table, get_array_of_tables,     &
                                     get_optional_int_with_default,      &
                                     get_optional_string_with_default,   &
                                     parse_date_to_days1900
   use error_mod, only: error_collection_t, ERR_PARSE_TYPE_MISMATCH
   implicit none
   private

   public :: read_crop_toml

contains

   subroutine read_crop_toml(doc, config, errors)
      type(toml_table), pointer, intent(in)    :: doc
      type(crop_config_t),       intent(inout) :: config
      type(error_collection_t),  intent(inout) :: errors

      type(toml_table), pointer :: sec, item
      type(toml_array), pointer :: rotation
      type(toml_datetime)       :: dtv
      integer :: i, n, stat
      character(len=:), allocatable :: fname

      call get_table(doc, 'crop', sec, 'crop', errors)
      if (.not. associated(sec)) return

      call get_optional_int_with_default(sec, 'swcrop', config%swcrop, 0, 'crop.swcrop', errors)

      call get_array_of_tables(sec, 'rotation', rotation, 'crop.rotation', errors)
      if (.not. associated(rotation)) return

      n = len(rotation)
      if (n == 0) return

      allocate(config%rotation_start(n), config%rotation_end(n), &
               config%rotation_file(n),  config%rotation_type(n))
      config%rotation_start = 0.0_real64
      config%rotation_end   = 0.0_real64
      config%rotation_file  = ""
      config%rotation_type  = 0

      do i = 1, n
         call get_value(rotation, i, item, stat=stat)
         if (stat /= 0 .or. .not. associated(item)) cycle

         call get_value(item, 'start', dtv, stat=stat)
         if (stat == 0) then
            config%rotation_start(i) = parse_date_to_days1900(dtv)
         else
            call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                               "expected date at crop.rotation.start", &
                               "crop.rotation.start")
         end if

         call get_value(item, 'end', dtv, stat=stat)
         if (stat == 0) then
            config%rotation_end(i) = parse_date_to_days1900(dtv)
         else
            call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                               "expected date at crop.rotation.end", &
                               "crop.rotation.end")
         end if

         call get_optional_string_with_default(item, 'file', fname, '', 'crop.rotation.file', errors)
         config%rotation_file(i) = fname

         call get_optional_int_with_default(item, 'type', config%rotation_type(i), 0, 'crop.rotation.type', errors)
      end do
   end subroutine read_crop_toml

end module read_crop_toml_mod
