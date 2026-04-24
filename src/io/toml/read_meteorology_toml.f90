!> Reader for the [meteorology] section of a SWAP TOML.
module read_meteorology_toml_mod
   use tomlf, only: toml_table
   use meteorology_config_mod, only: meteorology_config_t
   use toml_field_helpers_mod, only: get_table, &
                                     get_optional_string_with_default, &
                                     get_optional_int_with_default, &
                                     get_optional_real_with_default
   use error_mod, only: error_collection_t
   implicit none
   private

   public :: read_meteorology_toml

contains

   subroutine read_meteorology_toml(doc, config, errors)
      type(toml_table), pointer,   intent(in)    :: doc
      type(meteorology_config_t),  intent(inout) :: config
      type(error_collection_t),    intent(inout) :: errors

      type(toml_table), pointer :: sec, et, temporal, rain, interception

      call get_table(doc, 'meteorology', sec, 'meteorology', errors)
      if (.not. associated(sec)) return

      call get_optional_string_with_default(sec, 'file', config%metfile, '', 'meteorology.file', errors)
      call get_optional_real_with_default(sec,   'lat',  config%lat,     0.0d0, 'meteorology.lat',  errors)
      call get_optional_real_with_default(sec,   'alt',  config%alt,     0.0d0, 'meteorology.alt',  errors)
      call get_optional_real_with_default(sec,   'altw', config%altw,    2.0d0, 'meteorology.altw', errors)

      call get_table(sec, 'evapotranspiration', et, 'meteorology.evapotranspiration', errors)
      if (associated(et)) then
         call get_optional_int_with_default(et,  'swetr',      config%swetr,    0, 'meteorology.evapotranspiration.swetr',    errors)
         call get_optional_int_with_default(et,  'swdivide',   config%swdivide, 0, 'meteorology.evapotranspiration.swdivide', errors)
         call get_optional_real_with_default(et, 'angstrom_a', config%angstroma, 0.25d0, 'meteorology.evapotranspiration.angstrom_a', errors)
         call get_optional_real_with_default(et, 'angstrom_b', config%angstromb, 0.50d0, 'meteorology.evapotranspiration.angstrom_b', errors)
      end if

      call get_table(sec, 'temporal', temporal, 'meteorology.temporal', errors)
      if (associated(temporal)) then
         call get_optional_int_with_default(temporal, 'swmetdetail',  config%swmetdetail,  0, 'meteorology.temporal.swmetdetail',  errors)
         call get_optional_int_with_default(temporal, 'nmetdetail',   config%nmetdetail,   0, 'meteorology.temporal.nmetdetail',   errors)
         call get_optional_int_with_default(temporal, 'swmetfilall',  config%swmetfilall,  0, 'meteorology.temporal.swmetfilall',  errors)
      end if

      call get_table(sec, 'rain', rain, 'meteorology.rain', errors)
      if (associated(rain)) then
         call get_optional_int_with_default(rain, 'swrain',   config%swrain,   0, 'meteorology.rain.swrain',   errors)
         call get_optional_int_with_default(rain, 'swetsine', config%swetsine, 0, 'meteorology.rain.swetsine', errors)
      end if

      call get_table(sec, 'interception', interception, 'meteorology.interception', errors)
      if (associated(interception)) then
         call get_optional_int_with_default(interception, 'swinter', config%swinter, 0, 'meteorology.interception.swinter', errors)
      end if
   end subroutine read_meteorology_toml

end module read_meteorology_toml_mod
