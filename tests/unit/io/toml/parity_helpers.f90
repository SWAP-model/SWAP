!> Phase 4b Task 12: shared setup helper for hupselbrook parity tests.
!! Lives under tests/; never shipped with the production binary.
module parity_helpers_mod
   use iso_fortran_env, only: real64
   use variables
   use swap_config_mod
   use load_swap_config_mod
   use chdir_helper_mod
   use error_mod
   implicit none
   private

   public :: load_both_for_hupselbrook
   public :: load_both_for_salinitystress
   public :: load_both_for_macroporeflow
   public :: load_both_for_surfacewater
   public :: reset_for_next_readswap

   character(len=*), parameter :: CASE_DIR = &
      'tests/swap-cases/1.hupselbrook'
   character(len=*), parameter :: TOML_FILE = &
      'tests/swap-cases/toml/1.hupselbrook/swap.toml'
   character(len=*), parameter :: TEMPLATE = 'swap_linux.swp.template'

   character(len=*), parameter :: SAL_CASE_DIR = &
      'tests/swap-cases/5.salinitystress'
   character(len=*), parameter :: SAL_TOML_FILE = &
      'tests/swap-cases/toml/5.salinitystress/swap.toml'

   character(len=*), parameter :: MAC_CASE_DIR = &
      'tests/swap-cases/3.macroporeflow'
   character(len=*), parameter :: MAC_TOML_FILE = &
      'tests/swap-cases/toml/3.macroporeflow/swap.toml'

   character(len=*), parameter :: SW_CASE_DIR = &
      'tests/swap-cases/6.surfacewater'
   character(len=*), parameter :: SW_TOML_FILE = &
      'tests/swap-cases/toml/6.surfacewater/swap.toml'

contains

   !> Phase 4c-b: pFUnit runs every @test in a single process, but the legacy
   !! readswap() opens its log file with status='new'. Even though it asks for
   !! 'del' privilege, the file lingers when a previous test was aborted or
   !! closed mid-flight (and on success the close-with-delete only fires after
   !! readswap finishes — we may be re-entering before that). This helper
   !! force-deletes the stale `swap_swap.log` (and the legacy `Swap.ok` marker)
   !! and closes any leaked Fortran units in the getun-managed range so the
   !! next readswap() can claim a fresh logf.
   subroutine reset_for_next_readswap()
      integer :: u, ios
      logical :: is_open

      ! Defensive: close logf if a prior test left it open.
      inquire(unit=logf, opened=is_open)
      if (is_open) close(logf, iostat=ios)

      ! Forcibly delete the legacy log file via Fortran open+close-with-delete.
      ! Open as 'old' and unlink on close. iostat ignored — absence is fine.
      open(newunit=u, file='swap_swap.log', status='old', &
           action='read', iostat=ios)
      if (ios == 0) close(u, status='delete', iostat=ios)

      open(newunit=u, file='Swap.ok', status='old', &
           action='read', iostat=ios)
      if (ios == 0) close(u, status='delete', iostat=ios)

      ! Belt-and-suspenders for any other process that may have leaked a
      ! handle: try a shell-level rm so a still-open inode does not block
      ! the next status='new' open. Also nuke the rdinit scratch files
      ! (`<basename>rd$NNNNN.tmp`) and any stray `<project>_swap.log`
      ! variants that earlier direct binary invocations with -f/-t/-v
      ! flags may have produced. Test-only path.
      call execute_command_line( &
         'rm -f swap_swap.log Swap.ok *rd\$*.tmp *_swap.log 2>/dev/null', &
         wait=.true.)

      ! Close any leaked units in the getun-managed range so getun(20,99)
      ! finds a free slot. Iostat ignored — closing an already-closed unit
      ! is a no-op in gfortran.
      do u = 20, 99
         close(u, iostat=ios)
      end do

      ! Reset selected variables-module globals that readswap() does NOT
      ! always overwrite. Without this, a previous test's value leaks into
      ! the next case (e.g. case 6.surfacewater has swdra=2, where rddre()
      ! never assigns dramet — so it would inherit 3 from case 4 and the
      ! parity assertion would silently pass for the wrong reason).
      dramet = 0
   end subroutine reset_for_next_readswap

   !> Run readswap() for the hupselbrook case (sets `variables` module globals)
   !! then load+validate+finalize the matching TOML config.
   subroutine load_both_for_hupselbrook(config, errors)
      type(swap_config_t),      intent(out) :: config
      type(error_collection_t), intent(out) :: errors
      character(len=1024) :: orig_cwd

      call get_cwd(orig_cwd)
      call chdir_to(CASE_DIR)
      call stage_swp_template(TEMPLATE, 'swap')
      call reset_for_next_readswap()
      call readswap('swap')
      close(logf)
      call chdir_to(trim(orig_cwd))

      call load_swap_config(TOML_FILE, config, errors)
      call config%validate(errors)
      call config%finalize(errors)
   end subroutine load_both_for_hupselbrook

   !> Phase 4c-b Task 12: same pattern as `load_both_for_hupselbrook` but for
   !! the salinitystress (case 5) regression case.
   subroutine load_both_for_salinitystress(config, errors)
      type(swap_config_t),      intent(out) :: config
      type(error_collection_t), intent(out) :: errors
      character(len=1024) :: orig_cwd

      call get_cwd(orig_cwd)
      call chdir_to(SAL_CASE_DIR)
      call stage_swp_template(TEMPLATE, 'swap')
      call reset_for_next_readswap()
      call readswap('swap')
      close(logf)
      call chdir_to(trim(orig_cwd))

      call load_swap_config(SAL_TOML_FILE, config, errors)
      call config%validate(errors)
      call config%finalize(errors)
   end subroutine load_both_for_salinitystress

   !> Phase 4e Task C3: same pattern as `load_both_for_hupselbrook` but for
   !! the macroporeflow (case 3) regression case. Macropore-specific physics
   !! (SWMACRO=1 sub-block) is exercised by the legacy reader but has no
   !! schema slot; the parity test asserts only the schema-covered subset.
   subroutine load_both_for_macroporeflow(config, errors)
      type(swap_config_t),      intent(out) :: config
      type(error_collection_t), intent(out) :: errors
      character(len=1024) :: orig_cwd

      call get_cwd(orig_cwd)
      call chdir_to(MAC_CASE_DIR)
      call stage_swp_template(TEMPLATE, 'swap')
      call reset_for_next_readswap()
      call readswap('swap')
      close(logf)
      call chdir_to(trim(orig_cwd))

      call load_swap_config(MAC_TOML_FILE, config, errors)
      call config%validate(errors)
      call config%finalize(errors)
   end subroutine load_both_for_macroporeflow

   !> Phase 4f-prep Task D2: same pattern as `load_both_for_hupselbrook` but
   !! for the surfacewater (case 6) regression case. SWSRF=2/SWSEC=2 with
   !! 28 management periods exercises the full surface_water_config_t.
   !!
   !! READING NOTE: the surface-water management block in swap.dra is read
   !! by `rddre()` (defined in src/io/readswap.f90:4273), which is called
   !! from `SurfaceWater(1)` at simulation init — NOT from `readswap()`.
   !! To populate `variables%swsrf`, `swsec`, `nmper`, `impend(:)`, etc.,
   !! the helper invokes `rddre` directly after `readswap()` while still
   !! in the case dir. Mirrors the Phase 4c-b legacy_crop_helper pattern.
   subroutine load_both_for_surfacewater(config, errors)
      type(swap_config_t),      intent(out) :: config
      type(error_collection_t), intent(out) :: errors
      character(len=1024) :: orig_cwd
      real(real64) :: wls1, wlp1
      interface
         subroutine rddre(wls1, wlp1)
            use iso_fortran_env, only: real64
            real(real64), intent(out) :: wls1, wlp1
         end subroutine rddre
      end interface

      call get_cwd(orig_cwd)
      call chdir_to(SW_CASE_DIR)
      call stage_swp_template(TEMPLATE, 'swap')
      call reset_for_next_readswap()
      call readswap('swap')
      ! `rddre` populates surface-water management globals from swap.dra.
      ! Must run from the case dir before chdir-back.
      call rddre(wls1, wlp1)
      close(logf)
      call chdir_to(trim(orig_cwd))

      call load_swap_config(SW_TOML_FILE, config, errors)
      call config%validate(errors)
      call config%finalize(errors)
   end subroutine load_both_for_surfacewater

end module parity_helpers_mod
