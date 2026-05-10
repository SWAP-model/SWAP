!> Bottom-boundary config — populated from `[bottom_boundary]` TOML section.
!! Phase 4d Task 3: validators wired per `swbotb` branch (1..8).
!! Phase 4f cleanup (Tasks 6-8): replaces inline table slots with *_file paths
!! and sub-mode switches.
module bottom_boundary_config_mod
   use iso_fortran_env, only: real64
   use error_mod, only: error_collection_t, ERR_VALIDATION_OUT_OF_RANGE
   use validation_mod, only: check_int_enum, check_int_range, check_real_range
   implicit none
   private

   public :: bottom_boundary_config_t

   type :: bottom_boundary_config_t
      integer :: swbotb = 0

      ! Sub-mode switches (legacy SW2/SW3/SW4/SWQHBOT) — see core/variables.f90.
      ! sw2: 1=sine, 2=table; sw3: 1=sine, 2=table; sw4: 0=no extra flux,
      ! 1=include extra flux; swqhbot: 1=exponential, 2=tabular.
      ! Validators key on (swbotb, sw_x) combinations.
      integer :: sw2     = 1
      integer :: sw3     = 1
      integer :: sw4     = 0
      integer :: swqhbot = 1

      ! Phase 4f cleanup: per-sub-mode CSV companion file paths.
      character(len=:), allocatable :: gwl_file     ! SWBOTB=1
      character(len=:), allocatable :: qbot2_file   ! SWBOTB=2 + sw2=2
      character(len=:), allocatable :: haquif_file  ! SWBOTB=3 + sw3=2
      character(len=:), allocatable :: qbot4_file   ! SWBOTB=3 + sw4=1
      character(len=:), allocatable :: qhbot_file   ! SWBOTB=4 + swqhbot=2
      character(len=:), allocatable :: hbot5_file   ! SWBOTB=5

      ! SWBOTB=1 inline alternative (kept for backward compat with non-BBC paths)
      real(real64), allocatable :: swc_table(:,:)

      ! SWBOTB=2 sine-wave bottom flux scalars (sw2=1).
      ! Phase 0 B-0.1: promoted from legacy variables.f90 lines 947-949.
      ! Used in boundbottom.f90:104 when swbotb=2 .and. sw2=1.
      real(real64) :: sinmax = 0.0_real64  ! Day of year with maximum bottom flux
      real(real64) :: sinamp = 0.0_real64  ! Amplitude of bottom flux (L/T)
      real(real64) :: sinave = 0.0_real64  ! Average value of bottom flux (L/T)

      ! SWBOTB=4 exponential q(h) scalars (swqhbot=1).
      ! Phase 0 B-0.2: promoted from legacy variables.f90 lines 761-763, 720.
      ! Used in boundbottom.f90:152-153 when swbotb=4 .and. swqhbot=1.
      real(real64) :: cofqha   = 0.0_real64  ! Coefficient A: q = A * exp(B * |gwl|)
      real(real64) :: cofqhb   = 0.0_real64  ! Coefficient B: exponent multiplier (1/L)
      real(real64) :: cofqhc   = 0.0_real64  ! Coefficient C: additional flux term (L/T)
      integer      :: swcofqhc = 0           ! Switch: 1 = include c-term, 0 = omit

      ! SWBOTB=2 inline
      real(real64), allocatable :: qbot_table(:,:)

      ! SWBOTB=3 (Cauchy / regional aquifer)
      ! NOTE: legacy `shape` is a real (e.g. case 5 SHAPE=0.79). Phase 4d
      ! Task 20-prep: store as real(real64) and use a continuous range check.
      real(real64) :: shape   = 0.0_real64
      real(real64) :: hdrain  = 0.0_real64
      real(real64) :: rimlay  = 0.0_real64
      ! Phase 4f Task B4: numerical-solution switch for SWBOTB=3 bottom flux.
      ! Legacy default is 0 (explicit). Case 4 (oxygenstress) authors 1.
      integer      :: swbotb3impl = 0
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

   subroutine bottom_boundary_config_validate(self, errors)
      use error_mod, only: ERR_VALIDATION_REQUIRED
      class(bottom_boundary_config_t), intent(in)    :: self
      type(error_collection_t),         intent(inout) :: errors
      integer :: i

      if (self%swbotb == 0) return

      call check_int_enum(self%swbotb, [(i, i=1, 8)], 'bottom_boundary.swbotb', errors)
      call check_int_enum(self%sw2,     [1, 2], 'bottom_boundary.sw2',     errors)
      call check_int_enum(self%sw3,     [1, 2], 'bottom_boundary.sw3',     errors)
      call check_int_enum(self%sw4,     [0, 1], 'bottom_boundary.sw4',     errors)
      call check_int_enum(self%swqhbot, [1, 2], 'bottom_boundary.swqhbot', errors)

      select case (self%swbotb)
      case (1)
         if (.not. has_file(self%gwl_file)) then
            call errors%append(ERR_VALIDATION_REQUIRED, &
               "bottom_boundary.gwl_file required when swbotb=1", &
               'bottom_boundary')
         end if
      case (2)
         if (self%sw2 == 1) then
            ! Sine-wave path: validate sine parameters (Phase 0 B-0.1).
            call check_real_range(self%sinmax, 0.0_real64, 366.0_real64, &
                                  'bottom_boundary.sinmax', errors)
            call check_real_range(self%sinamp, -1.0e3_real64, 1.0e3_real64, &
                                  'bottom_boundary.sinamp', errors)
            call check_real_range(self%sinave, -1.0e3_real64, 1.0e3_real64, &
                                  'bottom_boundary.sinave', errors)
         end if
         if (self%sw2 == 2 .and. .not. has_file(self%qbot2_file)) then
            call errors%append(ERR_VALIDATION_REQUIRED, &
               "bottom_boundary.qbot2_file required when swbotb=2 and sw2=2", &
               'bottom_boundary')
         end if
      case (3)
         call check_real_range(self%shape, 0.0_real64, 2.0_real64, &
                               'bottom_boundary.shape', errors)
         call check_real_range(self%hdrain, -1.0e4_real64, 0.0_real64, &
                               'bottom_boundary.hdrain', errors)
         call check_real_range(self%rimlay, 0.0_real64, 1.0e5_real64, &
                               'bottom_boundary.rimlay', errors)
         call check_real_range(self%aqave, -1.0e4_real64, 1.0e3_real64, &
                               'bottom_boundary.aqave', errors)
         call check_real_range(self%aqamp, 0.0_real64, 1.0e3_real64, &
                               'bottom_boundary.aqamp', errors)
         call check_real_range(self%aqper, 0.0_real64, 366.0_real64, &
                               'bottom_boundary.aqper', errors)
         call check_real_range(self%aqtmax, 0.0_real64, 366.0_real64, &
                               'bottom_boundary.aqtmax', errors)
         call check_int_enum(self%swbotb3impl, [0, 1], &
                             'bottom_boundary.swbotb3impl', errors)
         if (self%sw3 == 2 .and. .not. has_file(self%haquif_file)) then
            call errors%append(ERR_VALIDATION_REQUIRED, &
               "bottom_boundary.haquif_file required when swbotb=3 and sw3=2", &
               'bottom_boundary')
         end if
         if (self%sw4 == 1 .and. .not. has_file(self%qbot4_file)) then
            call errors%append(ERR_VALIDATION_REQUIRED, &
               "bottom_boundary.qbot4_file required when swbotb=3 and sw4=1", &
               'bottom_boundary')
         end if
         if (allocated(self%cofqha_table)) then
            if (size(self%cofqha_table, 2) /= 2) then
               call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
                  "bottom_boundary.cofqha_table: expected 2 columns", &
                  'bottom_boundary.cofqha_table')
            end if
            if (size(self%cofqha_table, 1) < 1) then
               call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
                  "bottom_boundary.cofqha_table: expected at least 1 row", &
                  'bottom_boundary.cofqha_table')
            end if
         end if
      case (4)
         if (self%swqhbot == 1) then
            ! Exponential q(h) path: validate scalar coefficients (Phase 0 B-0.2).
            call check_real_range(self%cofqha, -1.0e10_real64, 1.0e10_real64, &
                                  'bottom_boundary.cofqha', errors)
            call check_real_range(self%cofqhb, -1.0e10_real64, 1.0e10_real64, &
                                  'bottom_boundary.cofqhb', errors)
            call check_real_range(self%cofqhc, -1.0e10_real64, 1.0e10_real64, &
                                  'bottom_boundary.cofqhc', errors)
            call check_int_enum(self%swcofqhc, [0, 1], &
                                'bottom_boundary.swcofqhc', errors)
         end if
         if (self%swqhbot == 2 .and. .not. has_file(self%qhbot_file)) then
            call errors%append(ERR_VALIDATION_REQUIRED, &
               "bottom_boundary.qhbot_file required when swbotb=4 and swqhbot=2", &
               'bottom_boundary')
         end if
      case (5)
         call check_real_range(self%hbot, -1.0e10_real64, 1.0e3_real64, &
                               'bottom_boundary.hbot', errors)
         call check_real_range(self%rhobot, -1.0e4_real64, 1.0e4_real64, &
                               'bottom_boundary.rhobot', errors)
         if (.not. has_file(self%hbot5_file)) then
            call errors%append(ERR_VALIDATION_REQUIRED, &
               "bottom_boundary.hbot5_file required when swbotb=5", &
               'bottom_boundary')
         end if
      case (6, 7, 8)
         ! No required tables.
      end select
   end subroutine bottom_boundary_config_validate

   pure function has_file(slot) result(yes)
      character(len=:), allocatable, intent(in) :: slot
      logical :: yes
      yes = allocated(slot)
      if (yes) yes = len_trim(slot) > 0
   end function has_file

   subroutine bottom_boundary_config_finalize(self, errors)
      class(bottom_boundary_config_t), intent(inout) :: self
      type(error_collection_t),         intent(inout) :: errors
      ! No-op for now.
   end subroutine bottom_boundary_config_finalize

end module bottom_boundary_config_mod
