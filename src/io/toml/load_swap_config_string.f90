!> @file load_swap_config_string.f90
!! SS-BMI2: in-memory TOML config loader. Mirrors load_swap_config
!! but consumes a character buffer instead of a file path. The cffi
!! entry point swap_initialize_from_toml_string calls this.
module load_swap_config_string_mod
   use tomlf, only: toml_table, toml_error, toml_loads
   use swap_config_mod, only: swap_config_t
   use load_swap_config_mod, only: apply_section_readers
   use error_mod, only: error_collection_t, ERR_PARSE_MALFORMED_TOML
   implicit none
   private
   public :: load_swap_config_from_string

contains

   subroutine load_swap_config_from_string(toml_text, config, errors)
      character(len=*),                    intent(in)    :: toml_text
      type(swap_config_t),                 intent(inout) :: config
      type(error_collection_t),            intent(inout) :: errors

      type(toml_table), allocatable, target :: doc
      type(toml_table), pointer             :: doc_ptr
      type(toml_error), allocatable         :: terr

      call toml_loads(doc, toml_text, error=terr)
      if (allocated(terr)) then
         call errors%append(ERR_PARSE_MALFORMED_TOML, trim(terr%message), "(in-memory)")
         return
      end if

      doc_ptr => doc
      call apply_section_readers(doc_ptr, ".", config, errors)
   end subroutine load_swap_config_from_string

end module load_swap_config_string_mod
