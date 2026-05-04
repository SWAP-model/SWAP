!> Reader for the [general] section of a SWAP TOML.
module read_general_toml_mod
   use tomlf, only: toml_table
   use general_config_mod, only: general_config_t
   use toml_field_helpers_mod, only: get_table, &
                                     get_required_string, &
                                     get_optional_string_with_default, &
                                     get_optional_int_with_default
   use error_mod, only: error_collection_t
   implicit none
   private

   public :: read_general_toml

contains

   subroutine read_general_toml(doc, config, errors)
      type(toml_table), pointer, intent(in)    :: doc
      type(general_config_t),    intent(inout) :: config
      type(error_collection_t),  intent(inout) :: errors

      type(toml_table), pointer :: sec, paths

      call get_table(doc, 'general', sec, 'general', errors)
      if (.not. associated(sec)) return

      call get_required_string(sec, 'project', config%project, 'general.project', errors)
      call get_optional_int_with_default(sec, 'swscre',  config%swscre,  0, 'general.swscre',  errors)
      call get_optional_int_with_default(sec, 'swerror', config%swerror, 0, 'general.swerror', errors)

      call get_optional_string_with_default(sec, 'outfil', config%outfil, 'result', &
                                            'general.outfil', errors)

      call get_table(sec, 'paths', paths, 'general.paths', errors)
      if (associated(paths)) then
         call get_optional_string_with_default(paths, 'work',       config%pathwork,  './',  'general.paths.work',       errors)
         call get_optional_string_with_default(paths, 'atmosphere', config%pathatm,   './',  'general.paths.atmosphere', errors)
         call get_optional_string_with_default(paths, 'crop',       config%pathcrop,  './',  'general.paths.crop',       errors)
         call get_optional_string_with_default(paths, 'drain',      config%pathdrain, './',  'general.paths.drain',      errors)
      end if
   end subroutine read_general_toml

end module read_general_toml_mod
