!> Reader for a type-1 (fixed crop) .crp.toml file.
!!
!! Given the ROOT table of a loaded .crp.toml document, populates a
!! cropfixed_config_t. Phase 4c-a covers a core schema; tables (cftb,
!! chtb, rdctb) are added during parity-test iteration when needed.
module read_cropfixed_toml_mod
   use iso_fortran_env, only: real64
   use tomlf, only: toml_table
   use cropfixed_config_mod, only: cropfixed_config_t
   use toml_field_helpers_mod, only: get_table,                    &
                                     get_optional_int_with_default, &
                                     get_optional_real_with_default
   use error_mod, only: error_collection_t
   implicit none
   private

   public :: read_cropfixed_toml

contains

   subroutine read_cropfixed_toml(doc, config, errors)
      type(toml_table), pointer,  intent(in)    :: doc
      type(cropfixed_config_t),   intent(inout) :: config
      type(error_collection_t),   intent(inout) :: errors

      type(toml_table), pointer :: ph, light, root, ws, salt, inter

      call get_table(doc, 'phenology', ph, 'phenology', errors)
      if (associated(ph)) then
         call get_optional_int_with_default(ph,  'idev', config%idev, 1, 'phenology.idev', errors)
         call get_optional_int_with_default(ph,  'lcc',  config%lcc,  0, 'phenology.lcc',  errors)
      end if

      call get_table(doc, 'light', light, 'light', errors)
      if (associated(light)) then
         call get_optional_real_with_default(light, 'kdif', config%kdif, 0.0_real64, 'light.kdif', errors)
         call get_optional_real_with_default(light, 'kdir', config%kdir, 0.0_real64, 'light.kdir', errors)
      end if

      call get_table(doc, 'root', root, 'root', errors)
      if (associated(root)) then
         call get_optional_real_with_default(root, 'rdi', config%rdi, 0.0_real64, 'root.rdi', errors)
         call get_optional_real_with_default(root, 'rri', config%rri, 0.0_real64, 'root.rri', errors)
         call get_optional_real_with_default(root, 'rdc', config%rdc, 0.0_real64, 'root.rdc', errors)
      end if

      call get_table(doc, 'water_stress', ws, 'water_stress', errors)
      if (associated(ws)) then
         call get_optional_real_with_default(ws, 'hlim1',  config%hlim1,  0.0_real64, 'ws.hlim1',  errors)
         call get_optional_real_with_default(ws, 'hlim2u', config%hlim2u, 0.0_real64, 'ws.hlim2u', errors)
         call get_optional_real_with_default(ws, 'hlim2l', config%hlim2l, 0.0_real64, 'ws.hlim2l', errors)
         call get_optional_real_with_default(ws, 'hlim3h', config%hlim3h, 0.0_real64, 'ws.hlim3h', errors)
         call get_optional_real_with_default(ws, 'hlim3l', config%hlim3l, 0.0_real64, 'ws.hlim3l', errors)
         call get_optional_real_with_default(ws, 'hlim4',  config%hlim4,  0.0_real64, 'ws.hlim4',  errors)
         call get_optional_real_with_default(ws, 'adcrh',  config%adcrh,  0.0_real64, 'ws.adcrh',  errors)
         call get_optional_real_with_default(ws, 'adcrl',  config%adcrl,  0.0_real64, 'ws.adcrl',  errors)
         call get_optional_real_with_default(ws, 'rsc',    config%rsc,    0.0_real64, 'ws.rsc',    errors)
      end if

      call get_table(doc, 'salinity', salt, 'salinity', errors)
      if (associated(salt)) then
         call get_optional_real_with_default(salt, 'ecmax',  config%ecmax,  0.0_real64, 'salt.ecmax',  errors)
         call get_optional_real_with_default(salt, 'ecslop', config%ecslop, 0.0_real64, 'salt.ecslop', errors)
      end if

      call get_table(doc, 'interception', inter, 'interception', errors)
      if (associated(inter)) then
         call get_optional_real_with_default(inter, 'cofab', config%cofab, 0.0_real64, 'inter.cofab', errors)
      end if
   end subroutine read_cropfixed_toml

end module read_cropfixed_toml_mod
