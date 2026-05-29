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

   !> Read the [general] section.
   !!
   !! `base_dir` is the directory containing the main TOML file (with trailing
   !! slash). It is used as the default for `general.paths.work` when the user
   !! has not set that key, making pathwork the canonical anchor for all
   !! relative paths (pathatm/pathcrop/pathdrain then default to pathwork).
   subroutine read_general_toml(doc, config, errors, base_dir)
      type(toml_table), pointer,  intent(in)    :: doc
      type(general_config_t),     intent(inout) :: config
      type(error_collection_t),   intent(inout) :: errors
      character(len=*), optional, intent(in)    :: base_dir

      type(toml_table), pointer :: sec, paths
      character(len=:), allocatable :: work_default

      call get_table(doc, 'general', sec, 'general', errors)
      if (.not. associated(sec)) return

      call get_required_string(sec, 'project', config%project, 'general.project', errors)
      call get_optional_int_with_default(sec, 'swscre',  config%swscre,  0, 'general.swscre',  errors)
      call get_optional_int_with_default(sec, 'swerror', config%swerror, 0, 'general.swerror', errors)

      call get_optional_string_with_default(sec, 'outfil', config%outfil, 'result', &
                                            'general.outfil', errors)

      ! Default chain: pathwork → base_dir (directory of main TOML file)
      !                pathatm/pathcrop/pathdrain → pathwork
      ! This makes pathwork the single canonical anchor for relative paths.
      if (present(base_dir)) then
         work_default = base_dir
      else
         work_default = './'
      end if

      call get_table(sec, 'paths', paths, 'general.paths', errors)
      if (associated(paths)) then
         call get_optional_string_with_default(paths, 'work', config%pathwork, work_default, &
                                               'general.paths.work', errors)
         call get_optional_string_with_default(paths, 'atmosphere', config%pathatm,   config%pathwork, &
                                               'general.paths.atmosphere', errors)
         call get_optional_string_with_default(paths, 'crop',       config%pathcrop,  config%pathwork, &
                                               'general.paths.crop',       errors)
         call get_optional_string_with_default(paths, 'drain',      config%pathdrain, config%pathwork, &
                                               'general.paths.drain',      errors)
      else
         ! No [general.paths] table — populate all four from work_default.
         config%pathwork  = work_default
         config%pathatm   = work_default
         config%pathcrop  = work_default
         config%pathdrain = work_default
      end if
   end subroutine read_general_toml

end module read_general_toml_mod
