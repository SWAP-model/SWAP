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
   use read_solute_toml_mod,       only: read_solute_toml
   use read_surface_water_toml_mod, only: read_surface_water_toml
   use read_crop_toml_mod,         only: read_crop_toml
   use read_output_csv_toml_mod,   only: read_output_csv_toml
   use read_nutrients_toml_mod,    only: read_nutrients_toml
   use load_crop_rotation_files_mod, only: load_crop_rotation_files
   use config_source_mod, only: config_source_t
   use error_mod, only: error_collection_t, ERR_PARSE_MALFORMED_TOML
   implicit none
   private

   public :: load_swap_config
   public :: apply_section_readers

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

      call apply_section_readers(doc_ptr, base_dir, config, errors)

      ! Load per-rotation .crp.toml subfiles. Separate phase because:
      ! (1) general%pathwork must be populated first (by read_general_toml).
      ! (2) The per-crop readers are I/O — separated from pure config parsing
      !     per interpretation A of W5.
      call load_crop_rotation_files(config%general, config%crop, errors)
   end subroutine load_swap_config

   !> Shared post-parse pipeline: takes a parsed toml_table and calls
   !! every section reader. Used by both file-based and string-based
   !! load entry points.
   !!
   !! Execution order is load-bearing: read_general_toml must run first so
   !! that config%general%pathwork is populated before read_drainage_toml
   !! uses it as its subfile base path and before load_crop_rotation_files
   !! (called from load_swap_config) follows .crp.toml file references.
   subroutine apply_section_readers(doc_ptr, base_dir, config, errors, source)
      use path_helpers_mod, only: directory_of
      type(toml_table), pointer,           intent(in)    :: doc_ptr
      character(len=*),                    intent(in)    :: base_dir
      type(swap_config_t),                 intent(inout) :: config
      type(error_collection_t),            intent(inout) :: errors
      type(config_source_t), optional,     intent(in)    :: source

      character(len=:), allocatable :: pathwork_eff

      call read_general_toml    (doc_ptr, config%general,    errors, base_dir=base_dir)
      call read_simulation_toml (doc_ptr, config%simulation, errors)
      call read_meteorology_toml(doc_ptr, config%meteo,      errors)

      ! Guard: if [general] was absent or returned early, pathwork is not
      ! allocated; fall back to base_dir so the subfile readers still work.
      if (allocated(config%general%pathwork)) then
         pathwork_eff = config%general%pathwork
      else
         pathwork_eff = base_dir
      end if

      call read_drainage_toml   (doc_ptr, config%drain,      errors, base_path=pathwork_eff, source=source)
      call read_soil_toml       (doc_ptr, config%soil,       errors)
      call read_bottom_boundary_toml(doc_ptr, config%bottom_boundary, errors)
      call read_heat_toml       (doc_ptr, config%heat,       errors)
      call read_irrigation_toml (doc_ptr, config%irrigation, errors)
      call read_solute_toml     (doc_ptr, config%solute,     errors)
      call read_surface_water_toml(doc_ptr, config%surface_water, errors)
      call read_crop_toml       (doc_ptr, config%crop,       errors)
      call read_output_csv_toml (doc_ptr, config%output_csv, errors)
      call read_nutrients_toml  (doc_ptr, config%nutrients,  errors)
   end subroutine apply_section_readers

end module load_swap_config_mod
