!> Loads per-rotation .crp.toml subfiles after the main TOML config
!! has been parsed. Iterates config%crop%rotation_file(:) and dispatches
!! to the appropriate per-crop reader (cropfixed/cropgrass/cropwofost)
!! based on config%crop%rotation_type(i).
!!
!! Called from load_swap_config after apply_section_readers completes,
!! so config%general%pathwork is populated by then.
module load_crop_rotation_files_mod
   use tomlf, only: toml_table, toml_load, toml_error
   use crop_config_mod, only: crop_config_t
   use general_config_mod, only: general_config_t
   use read_cropfixed_toml_mod, only: read_cropfixed_toml
   use read_cropgrass_toml_mod, only: read_cropgrass_toml
   use read_cropwofost_toml_mod, only: read_cropwofost_toml
   use path_helpers_mod, only: resolve_relative_path
   use error_mod, only: error_collection_t, ERR_PARSE_MALFORMED_TOML
   implicit none
   private

   public :: load_crop_rotation_files

contains

   subroutine load_crop_rotation_files(general, crop, errors)
      type(general_config_t),    intent(in)    :: general
      type(crop_config_t),       intent(inout) :: crop
      type(error_collection_t),  intent(inout) :: errors

      type(toml_table), allocatable, target :: crp_doc
      type(toml_table), pointer             :: crp_doc_ptr
      type(toml_error), allocatable         :: terr
      integer :: i, n
      logical :: file_exists
      character(len=:), allocatable :: file_abs, pathwork_eff

      if (.not. allocated(crop%rotation_file)) return
      n = size(crop%rotation_file)
      if (n == 0) return

      if (allocated(general%pathwork)) then
         pathwork_eff = general%pathwork
      else
         pathwork_eff = "./"
      end if

      do i = 1, n
         if (len_trim(crop%rotation_file(i)) == 0) cycle

         file_abs = resolve_relative_path(pathwork_eff, trim(crop%rotation_file(i)))
         inquire(file=trim(file_abs), exist=file_exists)
         if (.not. file_exists) cycle   ! file absent — skip silently

         call toml_load(crp_doc, trim(file_abs), error=terr)
         if (allocated(terr)) then
            call errors%append(ERR_PARSE_MALFORMED_TOML, trim(terr%message), trim(file_abs))
            cycle
         end if
         crp_doc_ptr => crp_doc
         select case (crop%rotation_type(i))
         case (1)
            call read_cropfixed_toml(crp_doc_ptr, crop%rotation_fixed(i), errors)
            crop%rotation_loaded(i) = .true.
         case (2)
            call read_cropwofost_toml(crp_doc_ptr, crop%rotation_wofost(i), errors)
            crop%rotation_loaded(i) = .true.
         case (3)
            call read_cropgrass_toml(crp_doc_ptr, crop%rotation_grass(i), errors)
            crop%rotation_loaded(i) = .true.
         end select
      end do
   end subroutine load_crop_rotation_files

end module load_crop_rotation_files_mod
