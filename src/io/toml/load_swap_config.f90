!> Top-level SWAP TOML loader. Dispatches to each section reader.
module load_swap_config_mod
   use tomlf, only: toml_table, toml_error, toml_load
   use swap_config_mod, only: swap_config_t
   use read_general_toml_mod,      only: read_general_toml
   use read_simulation_toml_mod,   only: read_simulation_toml
   use read_meteorology_toml_mod,  only: read_meteorology_toml
   use read_drainage_toml_mod,     only: read_drainage_toml
   use read_soil_toml_mod,         only: read_soil_toml
   use read_bottom_boundary_toml_mod, only: read_bottom_boundary_toml
   use read_heat_toml_mod,         only: read_heat_toml
   use read_irrigation_toml_mod,   only: read_irrigation_toml
   use read_crop_toml_mod,         only: read_crop_toml
   use error_mod, only: error_collection_t, ERR_PARSE_MALFORMED_TOML
   implicit none
   private

   public :: load_swap_config

contains

   subroutine load_swap_config(path, config, errors)
      use path_helpers_mod, only: directory_of
      character(len=*),          intent(in)    :: path
      type(swap_config_t),       intent(inout) :: config
      type(error_collection_t),  intent(inout) :: errors

      type(toml_table), allocatable, target :: doc
      type(toml_table), pointer             :: doc_ptr
      type(toml_error), allocatable         :: terr
      character(len=:), allocatable         :: base_dir

      call toml_load(doc, trim(path), error=terr)
      if (allocated(terr)) then
         call errors%append(ERR_PARSE_MALFORMED_TOML, trim(terr%message), trim(path))
         return
      end if

      doc_ptr => doc
      base_dir = directory_of(trim(path))

      call read_general_toml    (doc_ptr, config%general,    errors)
      call read_simulation_toml (doc_ptr, config%simulation, errors)
      call read_meteorology_toml(doc_ptr, config%meteo,      errors)
      call read_drainage_toml   (doc_ptr, config%drain,      errors, base_path=base_dir)
      call read_soil_toml       (doc_ptr, config%soil,       errors)
      call read_bottom_boundary_toml(doc_ptr, config%bottom_boundary, errors)
      call read_heat_toml       (doc_ptr, config%heat,       errors)
      call read_irrigation_toml (doc_ptr, config%irrigation, errors)
      call read_crop_toml       (doc_ptr, config%crop,       errors, base_path=base_dir)
   end subroutine load_swap_config

end module load_swap_config_mod
