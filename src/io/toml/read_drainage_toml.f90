!> TOML reader for the [drainage] section.
!! Supports inline drainage data OR external file via `[drainage].file = "..."`.
module read_drainage_toml_mod
   use iso_fortran_env, only: real64
   use tomlf, only: toml_table, toml_array, toml_error, toml_load, get_value, len
   use drainage_config_mod, only: drainage_config_t
   use toml_field_helpers_mod, only: get_table, get_array_of_tables,   &
                                     get_optional_int_with_default,    &
                                     get_optional_real_with_default,   &
                                     get_optional_string_with_default
   use path_helpers_mod, only: resolve_relative_path
   use error_mod, only: error_collection_t, ERR_PARSE_MALFORMED_TOML
   implicit none
   private

   public :: read_drainage_toml

contains

   subroutine read_drainage_toml(doc, config, errors, base_path)
      type(toml_table), pointer,  intent(in)    :: doc
      type(drainage_config_t),    intent(inout) :: config
      type(error_collection_t),   intent(inout) :: errors
      character(len=*), optional, intent(in)    :: base_path

      type(toml_table), pointer             :: drain_tab, ext_root, ext_sec
      type(toml_table), allocatable, target :: ext_doc
      type(toml_error), allocatable         :: terr
      character(len=:), allocatable         :: file_rel, file_abs

      call get_table(doc, 'drainage', drain_tab, 'drainage', errors)
      if (.not. associated(drain_tab)) return

      call get_optional_string_with_default(drain_tab, 'file', file_rel, '', 'drainage.file', errors)
      if (len_trim(file_rel) > 0 .and. present(base_path)) then
         file_abs = resolve_relative_path(base_path, file_rel)
         call toml_load(ext_doc, trim(file_abs), error=terr)
         if (allocated(terr)) then
            call errors%append(ERR_PARSE_MALFORMED_TOML, trim(terr%message), trim(file_abs))
            return
         end if
         ext_root => ext_doc
         call get_table(ext_root, 'drainage', ext_sec, 'drainage', errors)
         if (.not. associated(ext_sec)) return
         call read_drainage_inner(ext_sec, config, errors)
      else
         call read_drainage_inner(drain_tab, config, errors)
      end if
   end subroutine read_drainage_toml

   subroutine read_drainage_inner(sec, config, errors)
      type(toml_table), pointer,  intent(in)    :: sec
      type(drainage_config_t),    intent(inout) :: config
      type(error_collection_t),   intent(inout) :: errors

      type(toml_table), pointer :: basic, item
      type(toml_array), pointer :: levels
      integer :: i, n, stat

      call get_optional_int_with_default(sec, 'swdra',    config%swdra,    0, 'drainage.swdra',    errors)
      call get_optional_int_with_default(sec, 'dramet',   config%dramet,   0, 'drainage.dramet',   errors)
      call get_optional_int_with_default(sec, 'swdivd',   config%swdivd,   0, 'drainage.swdivd',   errors)
      call get_optional_int_with_default(sec, 'swdislay', config%swdislay, 0, 'drainage.swdislay', errors)
      call get_optional_int_with_default(sec, 'nrlevs',   config%nrlevs,   0, 'drainage.nrlevs',   errors)
      call get_optional_real_with_default(sec, 'altcu',   config%altcu,    0.0_real64, 'drainage.altcu', errors)

      call get_table(sec, 'basic', basic, 'drainage.basic', errors)
      if (associated(basic)) then
         call get_optional_real_with_default(basic, 'basegw', config%basegw, 0.0_real64, 'drainage.basic.basegw', errors)
         call get_optional_real_with_default(basic, 'entres', config%entres, 0.0_real64, 'drainage.basic.entres', errors)
         call get_optional_real_with_default(basic, 'shape',  config%shape,  0.0_real64, 'drainage.basic.shape',  errors)
      end if

      call get_array_of_tables(sec, 'levels', levels, 'drainage.levels', errors)
      if (associated(levels)) then
         n = len(levels)
         if (n > 0) then
            allocate(config%swdtyp(n), config%zbotdr(n), config%drares(n), &
                     config%infres(n), config%L(n),      config%gwlinf(n), &
                     config%rdrain(n), config%rinfi(n),  config%rentry(n), &
                     config%rexit(n),  config%widthr(n), config%taludr(n), &
                     config%swallo(n))
            config%swdtyp  = 0
            config%zbotdr  = 0.0_real64
            config%drares  = 0.0_real64
            config%infres  = 0.0_real64
            config%L       = 0.0_real64
            config%gwlinf  = 0.0_real64
            config%rdrain  = 0.0_real64
            config%rinfi   = 0.0_real64
            config%rentry  = 0.0_real64
            config%rexit   = 0.0_real64
            config%widthr  = 0.0_real64
            config%taludr  = 0.0_real64
            config%swallo  = 0
            do i = 1, n
               call get_value(levels, i, item, stat=stat)
               if (stat /= 0 .or. .not. associated(item)) cycle
               call get_optional_int_with_default(item,  'swdtyp', config%swdtyp(i), 0,           'drainage.levels.swdtyp', errors)
               call get_optional_real_with_default(item, 'zbotdr', config%zbotdr(i), 0.0_real64,  'drainage.levels.zbotdr', errors)
               call get_optional_real_with_default(item, 'drares', config%drares(i), 0.0_real64,  'drainage.levels.drares', errors)
               call get_optional_real_with_default(item, 'infres', config%infres(i), 0.0_real64,  'drainage.levels.infres', errors)
               call get_optional_real_with_default(item, 'L',      config%L(i),      0.0_real64,  'drainage.levels.L',      errors)
            end do
         end if
      end if
   end subroutine read_drainage_inner

end module read_drainage_toml_mod
