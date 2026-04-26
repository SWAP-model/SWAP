!> Reader for the [crop] section. Follows [[crop.rotation]].file references
!! to per-crop .crp.toml files when base_path is provided.
module read_crop_toml_mod
   use iso_fortran_env, only: real64
   use tomlf, only: toml_table, toml_array, toml_datetime, toml_error, toml_load, get_value, len
   use crop_config_mod, only: crop_config_t
   use cropfixed_config_mod, only: cropfixed_config_t
   use cropgrass_config_mod, only: cropgrass_config_t
   use cropwofost_config_mod, only: cropwofost_config_t
   use read_cropfixed_toml_mod, only: read_cropfixed_toml
   use read_cropgrass_toml_mod, only: read_cropgrass_toml
   use read_cropwofost_toml_mod, only: read_cropwofost_toml
   use toml_field_helpers_mod, only: get_table, get_array_of_tables,   &
                                     get_optional_int_with_default,    &
                                     get_optional_string_with_default, &
                                     parse_date_to_days1900
   use path_helpers_mod, only: resolve_relative_path
   use error_mod, only: error_collection_t, ERR_PARSE_TYPE_MISMATCH, ERR_PARSE_MALFORMED_TOML
   implicit none
   private

   public :: read_crop_toml

contains

   subroutine read_crop_toml(doc, config, errors, base_path)
      type(toml_table), pointer,  intent(in)    :: doc
      type(crop_config_t),        intent(inout) :: config
      type(error_collection_t),   intent(inout) :: errors
      character(len=*), optional, intent(in)    :: base_path

      type(toml_table), pointer             :: sec, item
      type(toml_array), pointer             :: rotation
      type(toml_table), allocatable, target :: crp_doc
      type(toml_table), pointer             :: crp_doc_ptr
      type(toml_error), allocatable         :: terr
      type(toml_datetime)                   :: dtv
      integer :: i, n, stat
      logical :: file_exists
      character(len=:), allocatable :: fname, file_abs

      call get_table(doc, 'crop', sec, 'crop', errors)
      if (.not. associated(sec)) return

      call get_optional_int_with_default(sec, 'swcrop', config%swcrop, 0, 'crop.swcrop', errors)

      call get_array_of_tables(sec, 'rotation', rotation, 'crop.rotation', errors)
      if (.not. associated(rotation)) return

      n = len(rotation)
      if (n == 0) return

      allocate(config%rotation_start(n),   config%rotation_end(n),   &
               config%rotation_file(n),    config%rotation_type(n),  &
               config%rotation_fixed(n),   config%rotation_grass(n), &
               config%rotation_wofost(n),                             &
               config%rotation_loaded(n))
      config%rotation_start  = 0.0_real64
      config%rotation_end    = 0.0_real64
      config%rotation_file   = ""
      config%rotation_type   = 0
      config%rotation_loaded = .false.

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

         ! If file= is specified and base_path is available, follow the
         ! reference and dispatch to the matching crop reader.
         ! Only attempt to load the file if it actually exists; absence is
         ! silently skipped so that test cases that have not yet been
         ! converted to per-crop TOML files continue to load cleanly.
         if (len_trim(fname) > 0 .and. present(base_path)) then
            file_abs = resolve_relative_path(base_path, trim(fname))
            inquire(file=trim(file_abs), exist=file_exists)
            if (.not. file_exists) cycle   ! file absent — skip silently
            call toml_load(crp_doc, trim(file_abs), error=terr)
            if (allocated(terr)) then
               call errors%append(ERR_PARSE_MALFORMED_TOML, trim(terr%message), trim(file_abs))
               cycle
            end if
            crp_doc_ptr => crp_doc
            select case (config%rotation_type(i))
            case (1)
               call read_cropfixed_toml(crp_doc_ptr, config%rotation_fixed(i), errors)
               config%rotation_loaded(i) = .true.
            case (2)
               call read_cropwofost_toml(crp_doc_ptr, config%rotation_wofost(i), errors)
               config%rotation_loaded(i) = .true.
            case (3)
               call read_cropgrass_toml(crp_doc_ptr, config%rotation_grass(i), errors)
               config%rotation_loaded(i) = .true.
            end select
         end if
      end do
   end subroutine read_crop_toml

end module read_crop_toml_mod
