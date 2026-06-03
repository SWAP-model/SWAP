!> Reader for the [output.csv] section of a SWAP TOML.
module read_output_csv_toml_mod
   use iso_fortran_env, only: real64
   use tomlf, only: toml_table, toml_array, len
   use output_csv_config_mod, only: output_csv_config_t
   use toml_field_helpers_mod, only: get_table, &
                                     get_optional_int_with_default, &
                                     get_optional_string_with_default
   use error_mod, only: error_collection_t
   implicit none
   private

   public :: read_output_csv_toml

contains

   subroutine read_output_csv_toml(doc, config, errors)
      type(toml_table), pointer,    intent(in)    :: doc
      type(output_csv_config_t),    intent(inout) :: config
      type(error_collection_t),     intent(inout) :: errors

      type(toml_table), pointer :: out_sec, csv_sec

      call get_table(doc, 'output', out_sec, 'output', errors)
      if (.not. associated(out_sec)) return

      call get_table(out_sec, 'csv', csv_sec, 'output.csv', errors)
      if (.not. associated(csv_sec)) return

      call get_optional_int_with_default(csv_sec, 'enabled', config%enabled, 1, &
                                         'output.csv.enabled', errors)
      call get_optional_int_with_default(csv_sec, 'enabled_tz', config%enabled_tz, 0, &
                                         'output.csv.enabled_tz', errors)
      call get_optional_int_with_default(csv_sec, 'results_in_memory', config%results_in_memory, &
                                         0, 'output.csv.results_in_memory', errors)
      call get_optional_string_with_default(csv_sec, 'inlist', config%inlist, &
         'rain,irrig,interc,runoff,drainage,dstor,epot,eact,tpot,tact,qbottom,gwl', &
         'output.csv.inlist', errors)
      call get_optional_string_with_default(csv_sec, 'inlist_tz', config%inlist_tz, &
         'wc,h,conc', 'output.csv.inlist_tz', errors)
      ! tz_z1_z2 = [z_top, z_bottom] cm; default [0.0, 0.0] preserves legacy
      ! csv_out_tz behaviour (top compartment only).
      call read_tz_z1_z2(csv_sec, config%tz_z1_z2)
   end subroutine read_output_csv_toml

   subroutine read_tz_z1_z2(sec, tz_z1_z2)
      use tomlf, only: get_value
      type(toml_table), pointer, intent(in)    :: sec
      real(real64),              intent(inout) :: tz_z1_z2(2)
      type(toml_array), pointer :: arr
      integer :: stat

      call get_value(sec, 'tz_z1_z2', arr, requested=.false.)
      if (.not. associated(arr)) return
      if (len(arr) >= 1) call get_value(arr, 1, tz_z1_z2(1), stat=stat)
      if (len(arr) >= 2) call get_value(arr, 2, tz_z1_z2(2), stat=stat)
   end subroutine read_tz_z1_z2

end module read_output_csv_toml_mod
