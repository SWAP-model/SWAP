!> Reader for the [soil] section of a SWAP TOML.
!!
!! Phase 4a scope: top-level scalars + initial conditions. Per-layer
!! arrays (sublay, hcomp, orgmat, etc.) and the 2D cofgen array are
!! loader-skeleton only; the hupselbrook parity test drives completion
!! of per-layer paths.
module read_soil_toml_mod
   use iso_fortran_env, only: real64
   use tomlf, only: toml_table
   use soil_config_mod, only: soil_config_t
   use toml_field_helpers_mod, only: get_table, &
                                     get_optional_int_with_default, &
                                     get_optional_real_with_default
   use error_mod, only: error_collection_t
   implicit none
   private

   public :: read_soil_toml

contains

   subroutine read_soil_toml(doc, config, errors)
      type(toml_table), pointer, intent(in)    :: doc
      type(soil_config_t),       intent(inout) :: config
      type(error_collection_t),  intent(inout) :: errors

      type(toml_table), pointer :: sec, initial

      call get_table(doc, 'soil', sec, 'soil', errors)
      if (.not. associated(sec)) return

      call get_optional_int_with_default(sec, 'swsophy', config%swsophy, 0, 'soil.swsophy', errors)
      call get_optional_int_with_default(sec, 'swhyst',  config%swhyst,  0, 'soil.swhyst',  errors)
      call get_optional_int_with_default(sec, 'swinco',  config%swinco,  1, 'soil.swinco',  errors)
      call get_optional_int_with_default(sec, 'swmacro', config%swmacro, 0, 'soil.swmacro', errors)
      call get_optional_int_with_default(sec, 'swscal',  config%swscal,  0, 'soil.swscal',  errors)

      call get_optional_real_with_default(sec, 'ksatexm', config%ksatexm, 0.0_real64, 'soil.ksatexm', errors)
      call get_optional_real_with_default(sec, 'rsoil',   config%rsoil,   0.0_real64, 'soil.rsoil',   errors)
      call get_optional_int_with_default(sec, 'reva_top', config%reva_top, 0, 'soil.reva_top', errors)

      call get_table(sec, 'initial', initial, 'soil.initial', errors)
      if (associated(initial)) then
         call get_optional_real_with_default(initial, 'gwli',    config%gwli,    0.0_real64, 'soil.initial.gwli',    errors)
         call get_optional_real_with_default(initial, 'pondini', config%pondini, 0.0_real64, 'soil.initial.pondini', errors)
         call get_optional_real_with_default(initial, 'pondmx',  config%pondmx,  0.0_real64, 'soil.initial.pondmx',  errors)
      end if
   end subroutine read_soil_toml

end module read_soil_toml_mod
