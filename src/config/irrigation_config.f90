!> Irrigation config — Phase 4d Task 9 skeleton.
!! Two public types in one module: `irrigation_config_t` holds top-level
!! `.swp` fixed-irrigation parameters; `irrigation_schedule_t` is a nested
!! sub-type embedded in each per-crop config (cropfixed, cropwofost,
!! cropgrass) for `.crp` scheduling. Validators land in Task 10; readers
!! in Task 11; swap_config wiring in Task 12. Field names mirror the
!! legacy `.swp`/`.crp` keys verbatim.
module irrigation_config_mod
   use iso_fortran_env, only: real64
   use error_mod, only: error_collection_t
   implicit none
   private

   public :: irrigation_config_t
   public :: irrigation_schedule_t

   ! Top-level fixed-irrigation (.swp side).
   type :: irrigation_config_t
      integer                       :: swirfix = 0
      character(len=:), allocatable :: irgfil
      real(real64),     allocatable :: fixed_events(:,:)
   contains
      procedure :: validate => irrigation_config_validate
      procedure :: finalize => irrigation_config_finalize
   end type irrigation_config_t

   ! Nested per-crop scheduling (.crp side).
   type :: irrigation_schedule_t
      integer      :: schedule = 0
      integer      :: startirr_day = 0
      integer      :: startirr_month = 0
      integer      :: endirr_day = 0
      integer      :: endirr_month = 0
      real(real64) :: cirrs = 0.0_real64
      integer      :: isuas = 0
      integer      :: tcs = 0
      integer      :: dcs = 0
      real(real64), allocatable :: trel_table(:,:)
      real(real64), allocatable :: raw_table(:,:)
      real(real64), allocatable :: taw_table(:,:)
      real(real64), allocatable :: dwa_table(:,:)
      real(real64) :: irgthreshold = 0.0_real64
      real(real64), allocatable :: hcri_table(:,:)
      real(real64), allocatable :: tcri_table(:,:)
      real(real64) :: dcrit = 0.0_real64
      integer      :: swcirrthres = 0
      real(real64) :: cirrthres = 0.0_real64
      real(real64) :: perirrsurp = 0.0_real64
      integer      :: tcsfix = 0
      integer      :: irgdayfix = 0
      real(real64) :: phfieldcapacity = 0.0_real64
      real(real64), allocatable :: di_table(:,:)
      real(real64) :: raithreshold = 0.0_real64
      real(real64), allocatable :: fid_table(:,:)
      integer      :: dcslim = 0
      real(real64) :: irgdepmin = 0.0_real64
      real(real64) :: irgdepmax = 0.0_real64
   contains
      procedure :: validate => irrigation_schedule_validate
      procedure :: finalize => irrigation_schedule_finalize
   end type irrigation_schedule_t

contains

   subroutine irrigation_config_validate(self, errors)
      class(irrigation_config_t), intent(in)    :: self
      type(error_collection_t),   intent(inout) :: errors
      ! TODO: Phase 4d Task 10
   end subroutine irrigation_config_validate

   subroutine irrigation_config_finalize(self, errors)
      class(irrigation_config_t), intent(inout) :: self
      type(error_collection_t),   intent(inout) :: errors
      ! No-op.
   end subroutine irrigation_config_finalize

   subroutine irrigation_schedule_validate(self, errors)
      class(irrigation_schedule_t), intent(in)    :: self
      type(error_collection_t),     intent(inout) :: errors
      ! TODO: Phase 4d Task 10
   end subroutine irrigation_schedule_validate

   subroutine irrigation_schedule_finalize(self, errors)
      class(irrigation_schedule_t), intent(inout) :: self
      type(error_collection_t),     intent(inout) :: errors
      ! No-op.
   end subroutine irrigation_schedule_finalize

end module irrigation_config_mod
