!> Bottom-boundary config — populated from `[bottom_boundary]` TOML section.
!! Phase 4d Task 2: type skeleton only. Validators land in Task 3.
module bottom_boundary_config_mod
   use iso_fortran_env, only: real64
   use error_mod, only: error_collection_t
   implicit none
   private

   public :: bottom_boundary_config_t

   type :: bottom_boundary_config_t
      integer :: swbotb = 0

      ! Optional path to external file for SWBOTB=1 reference (legacy bbcfil).
      character(len=:), allocatable :: bbcfil

      ! SWBOTB=1 inline alternative
      real(real64), allocatable :: swc_table(:,:)

      ! SWBOTB=2 inline
      real(real64), allocatable :: qbot_table(:,:)

      ! SWBOTB=3 (Cauchy / regional aquifer)
      integer      :: shape   = 0
      real(real64) :: hdrain  = 0.0_real64
      real(real64) :: rimlay  = 0.0_real64
      real(real64) :: aqave   = 0.0_real64
      real(real64) :: aqamp   = 0.0_real64
      real(real64) :: aqomeg  = 0.0_real64
      real(real64), allocatable :: cofqha_table(:,:)

      ! SWBOTB=5
      real(real64) :: hbot    = 0.0_real64
      real(real64) :: rhobot  = 0.0_real64
   contains
      procedure :: validate => bottom_boundary_config_validate
      procedure :: finalize => bottom_boundary_config_finalize
   end type bottom_boundary_config_t

contains

   subroutine bottom_boundary_config_validate(self, errors)
      class(bottom_boundary_config_t), intent(in)    :: self
      type(error_collection_t),         intent(inout) :: errors
      ! TODO: Phase 4d Task 3
   end subroutine bottom_boundary_config_validate

   subroutine bottom_boundary_config_finalize(self, errors)
      class(bottom_boundary_config_t), intent(inout) :: self
      type(error_collection_t),         intent(inout) :: errors
      ! No-op for now.
   end subroutine bottom_boundary_config_finalize

end module bottom_boundary_config_mod
