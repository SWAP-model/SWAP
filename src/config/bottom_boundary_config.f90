!> Bottom-boundary config — populated from `[bottom_boundary]` TOML section.
!! Phase 4d Task 3: validators wired per `swbotb` branch (1..8).
module bottom_boundary_config_mod
   use iso_fortran_env, only: real64
   use error_mod, only: error_collection_t, ERR_VALIDATION_OUT_OF_RANGE
   use validation_mod, only: check_int_enum, check_int_range, check_real_range
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
      ! NOTE: legacy `shape` is a real (e.g. case 5 SHAPE=0.79). Phase 4d
      ! Task 20-prep: store as real(real64) and use a continuous range check.
      real(real64) :: shape   = 0.0_real64
      real(real64) :: hdrain  = 0.0_real64
      real(real64) :: rimlay  = 0.0_real64
      real(real64) :: aqave   = 0.0_real64
      real(real64) :: aqamp   = 0.0_real64
      real(real64) :: aqper   = 0.0_real64
      real(real64) :: aqtmax  = 0.0_real64
      real(real64), allocatable :: cofqha_table(:,:)

      ! SWBOTB=5
      real(real64) :: hbot    = 0.0_real64
      real(real64) :: rhobot  = 0.0_real64
   contains
      procedure :: validate => bottom_boundary_config_validate
      procedure :: finalize => bottom_boundary_config_finalize
   end type bottom_boundary_config_t

contains

   !> Local helper: verify allocatable 2D table has expected ncols and >=1 row.
   !! Skips silently when unallocated (caller decides whether absence is an error).
   subroutine check_table_2d(table, expected_cols, label, errors)
      real(real64), allocatable, intent(in)    :: table(:,:)
      integer,                   intent(in)    :: expected_cols
      character(len=*),          intent(in)    :: label
      type(error_collection_t),  intent(inout) :: errors
      character(len=128) :: msg
      integer :: nrows, ncols
      if (.not. allocated(table)) return
      nrows = size(table, 1)
      ncols = size(table, 2)
      if (ncols /= expected_cols) then
         write(msg, '("ncols=",I0," expected ",I0)') ncols, expected_cols
         call errors%append(ERR_VALIDATION_OUT_OF_RANGE, trim(msg), label)
      end if
      if (nrows < 1) then
         call errors%append(ERR_VALIDATION_OUT_OF_RANGE, "nrows<1", label)
      end if
   end subroutine check_table_2d

   subroutine bottom_boundary_config_validate(self, errors)
      class(bottom_boundary_config_t), intent(in)    :: self
      type(error_collection_t),         intent(inout) :: errors
      logical :: have_file, have_table
      integer :: i

      ! Sentinel: swbotb=0 means the [bottom_boundary] section was absent in the
      ! TOML (reader leaves config at defaults). Skip validation so existing case
      ! TOMLs without the section still load clean — Phase 4d Tasks 17-19 will
      ! add [bottom_boundary] sections to each per-case TOML, after which
      ! swbotb will always be in [1..8] and full validation runs.
      if (self%swbotb == 0) return

      call check_int_enum(self%swbotb, [(i, i=1, 8)], 'bottom_boundary.swbotb', errors)

      select case (self%swbotb)
      case (1)
         have_file  = allocated(self%bbcfil)
         if (have_file) have_file = len_trim(self%bbcfil) > 0
         have_table = allocated(self%swc_table)
         if (.not. have_file .and. .not. have_table) then
            call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
               "swbotb=1 requires bbcfil or swc_table", 'bottom_boundary')
         end if
         if (have_file .and. have_table) then
            call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
               "swbotb=1 cannot set both bbcfil and swc_table", 'bottom_boundary')
         end if
         call check_table_2d(self%swc_table, 2, 'bottom_boundary.swc_table', errors)
      case (2)
         if (.not. allocated(self%qbot_table)) then
            call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
               "swbotb=2 requires qbot_table", 'bottom_boundary')
         end if
         call check_table_2d(self%qbot_table, 2, 'bottom_boundary.qbot_table', errors)
      case (3)
         call check_real_range(self%shape, 0.0_real64, 2.0_real64, &
                               'bottom_boundary.shape', errors)
         call check_real_range(self%hdrain, -1.0e4_real64,  0.0_real64, &
                               'bottom_boundary.hdrain', errors)
         call check_real_range(self%rimlay,  0.0_real64,    1.0e5_real64, &
                               'bottom_boundary.rimlay', errors)
         call check_real_range(self%aqave,  -1.0e4_real64,  1.0e3_real64, &
                               'bottom_boundary.aqave', errors)
         call check_real_range(self%aqamp,   0.0_real64,    1.0e3_real64, &
                               'bottom_boundary.aqamp', errors)
         call check_real_range(self%aqper,   0.0_real64,    366.0_real64, &
                               'bottom_boundary.aqper', errors)
         call check_real_range(self%aqtmax,  0.0_real64,    366.0_real64, &
                               'bottom_boundary.aqtmax', errors)
         call check_table_2d(self%cofqha_table, 2, 'bottom_boundary.cofqha_table', errors)
      case (5)
         call check_real_range(self%hbot,   -1.0e10_real64, 1.0e3_real64, &
                               'bottom_boundary.hbot', errors)
         call check_real_range(self%rhobot, -1.0e4_real64,  1.0e4_real64, &
                               'bottom_boundary.rhobot', errors)
      case (4, 6, 7, 8)
         ! No additional scalar params required.
      end select
   end subroutine bottom_boundary_config_validate

   subroutine bottom_boundary_config_finalize(self, errors)
      class(bottom_boundary_config_t), intent(inout) :: self
      type(error_collection_t),         intent(inout) :: errors
      ! No-op for now.
   end subroutine bottom_boundary_config_finalize

end module bottom_boundary_config_mod
