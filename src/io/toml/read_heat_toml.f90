!> Reader for the [heat] section of swap.toml.
!!
!! Populates a heat_config_t from a fully-formed swap.toml document. The
!! per-layer texture fractions (`psand`, `pclay`, `porg`) are encoded as flat
!! TOML arrays of real, decoded into 1-D allocatable real(real64) arrays.
!! The initial soil-temperature table (`tsoil_init`) is encoded as a TOML
!! array of arrays of two reals (depth_cm, temp_C) and decoded into a
!! 2-D allocatable real(real64) array.
!!
!! Mechanical pattern (mirrors read_bottom_boundary_toml): scalars one per
!! line, arrays/tables one block per array. No switch-conditional skipping
!! — the validator handles cross-field consistency. If [heat] is absent
!! the reader returns silently and leaves `config` at defaults.
module read_heat_toml_mod
   use iso_fortran_env, only: real64
   use tomlf, only: toml_table
   use heat_config_mod, only: heat_config_t
   use toml_field_helpers_mod, only: get_table,                         &
                                     get_optional_int_with_default,     &
                                     get_optional_real_with_default
   use toml_array_helpers_mod, only: read_real_array_1d, read_table_2d
   use error_mod, only: error_collection_t
   implicit none
   private

   public :: read_heat_toml

contains

   subroutine read_heat_toml(doc_root, config, errors)
      type(toml_table), pointer, intent(in)    :: doc_root
      type(heat_config_t),       intent(inout) :: config
      type(error_collection_t),  intent(inout) :: errors

      type(toml_table), pointer :: sec

      call get_table(doc_root, 'heat', sec, 'heat', errors)
      if (.not. associated(sec)) return

      ! Switches.
      call get_optional_int_with_default(sec, 'swhea',     config%swhea,     0, &
                                         'heat.swhea',     errors)
      call get_optional_int_with_default(sec, 'swcalt',    config%swcalt,    0, &
                                         'heat.swcalt',    errors)
      call get_optional_int_with_default(sec, 'swtopbhea', config%swtopbhea, 0, &
                                         'heat.swtopbhea', errors)
      call get_optional_int_with_default(sec, 'swbotbhea', config%swbotbhea, 0, &
                                         'heat.swbotbhea', errors)

      ! Frost params (scalars).
      call get_optional_real_with_default(sec, 'tfroststa', config%tfroststa, 0.0_real64, &
                                          'heat.tfroststa', errors)
      call get_optional_real_with_default(sec, 'tfrostend', config%tfrostend, 0.0_real64, &
                                          'heat.tfrostend', errors)

      ! Per-layer texture fractions (1-D arrays).
      call read_real_array_1d(sec, 'psand', config%psand, 'heat.psand', errors)
      call read_real_array_1d(sec, 'psilt', config%psilt, 'heat.psilt', errors)
      call read_real_array_1d(sec, 'pclay', config%pclay, 'heat.pclay', errors)
      call read_real_array_1d(sec, 'porg',  config%porg,  'heat.porg',  errors)

      ! Initial soil-temperature table (2-D, depth/temp pairs).
      call read_table_2d(sec, 'tsoil_init', config%tsoil_init, 2, &
                         'heat.tsoil_init', errors)

      ! Phase 0 (SS-HEAT) — swcalt=1 analytical method scalars.
      call get_optional_real_with_default(sec, 'ddamp',  config%ddamp,  0.0_real64, &
                                          'heat.ddamp',  errors)
      call get_optional_real_with_default(sec, 'tmean',  config%tmean,  0.0_real64, &
                                          'heat.tmean',  errors)
      call get_optional_real_with_default(sec, 'tampli', config%tampli, 0.0_real64, &
                                          'heat.tampli', errors)
      call get_optional_real_with_default(sec, 'timref', config%timref, 0.0_real64, &
                                          'heat.timref', errors)

      ! Boundary-condition temperature tables (2-D, time/temperature pairs).
      ! Absent key leaves table unallocated; adapter flattens to interleaved
      ! 1D afgen layout in apply_heat (2*k-1=time, 2*k=value).
      call read_table_2d(sec, 'temtoptab', config%temtoptab, 2, &
                         'heat.temtoptab', errors)
      call read_table_2d(sec, 'tembtab',   config%tembtab,   2, &
                         'heat.tembtab',   errors)
   end subroutine read_heat_toml

end module read_heat_toml_mod
