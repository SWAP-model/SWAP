!> Phase 4b Task 12: shared setup helper for hupselbrook parity tests.
!! Lives under tests/; never shipped with the production binary.
module parity_helpers_mod
   use variables
   use swap_config_mod
   use load_swap_config_mod
   use chdir_helper_mod
   use error_mod
   implicit none
   private

   public :: load_both_for_hupselbrook

   character(len=*), parameter :: CASE_DIR = &
      'tests/swap-cases/1.hupselbrook'
   character(len=*), parameter :: TOML_FILE = &
      'tests/swap-cases/toml/1.hupselbrook/swap.toml'
   character(len=*), parameter :: TEMPLATE = 'swap_linux.swp.template'

contains

   !> Run readswap() for the hupselbrook case (sets `variables` module globals)
   !! then load+validate+finalize the matching TOML config.
   subroutine load_both_for_hupselbrook(config, errors)
      type(swap_config_t),      intent(out) :: config
      type(error_collection_t), intent(out) :: errors
      character(len=1024) :: orig_cwd

      call get_cwd(orig_cwd)
      call chdir_to(CASE_DIR)
      call stage_swp_template(TEMPLATE, 'swap')
      call readswap()
      close(logf)
      call chdir_to(trim(orig_cwd))

      call load_swap_config(TOML_FILE, config, errors)
      call config%validate(errors)
      call config%finalize(errors)
   end subroutine load_both_for_hupselbrook

end module parity_helpers_mod
