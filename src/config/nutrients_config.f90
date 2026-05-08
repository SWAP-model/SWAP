!> @file nutrients_config.f90
!! Top-level [nutrients] config block: soil-side initial pool
!! concentrations + SorpCoef sorption coefficient.
!!
!! N2a of the [nutrients] umbrella (ADR 0026). Block is optional;
!! absent → present=.false., all defaults zero.
module nutrients_config_mod
   use, intrinsic :: iso_fortran_env, only: real64
   use error_mod, only: error_collection_t, ERR_VALIDATION_OUT_OF_RANGE
   use validation_mod, only: check_real_range
   implicit none
   private

   public :: nutrients_initial_t
   public :: nutrients_config_t

   !> Initial concentrations of soil organic-matter and N pools.
   !! Units: kg/m^2 for FOM/Bio/Hum (depth-integrated); kg/m^3 for
   !! cNH4/cNO3. Legacy ranges (rdsdor [0, 1000]) preserved.
   type :: nutrients_initial_t
      real(real64) :: fom(8) = 0.0_real64
      real(real64) :: bio    = 0.0_real64
      real(real64) :: hum    = 0.0_real64
      real(real64) :: cnh4   = 0.0_real64
      real(real64) :: cno3   = 0.0_real64
   end type nutrients_initial_t

   !> [nutrients] top-level block.
   type :: nutrients_config_t
      logical      :: present   = .false.
      real(real64) :: sorp_coef = 0.0_real64
      character(len=:), allocatable :: events_file       ! relative to pathwork; CSV companion (N2b, ADR 0027)
      type(nutrients_initial_t) :: initial
   contains
      procedure :: validate => nutrients_config_validate
      procedure :: finalize => nutrients_config_finalize
   end type nutrients_config_t

contains

   subroutine nutrients_config_validate(self, errors)
      class(nutrients_config_t), intent(in)    :: self
      type(error_collection_t),  intent(inout) :: errors
      integer :: i

      if (.not. self%present) return

      ! sorp_coef must be non-negative; no upper bound (literature
      ! values vary widely; high values may be intentional).
      if (self%sorp_coef < 0.0_real64) then
         call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
            'nutrients.sorp_coef: must be >= 0', &
            'nutrients.sorp_coef')
      end if

      ! Initial pool concentrations: legacy rdsdor [0, 1000] range.
      do i = 1, 8
         call check_real_range(self%initial%fom(i), 0.0_real64, 1000.0_real64, &
                               'nutrients.initial.fom', errors)
      end do
      call check_real_range(self%initial%bio,  0.0_real64, 1000.0_real64, &
                            'nutrients.initial.bio',  errors)
      call check_real_range(self%initial%hum,  0.0_real64, 1000.0_real64, &
                            'nutrients.initial.hum',  errors)
      call check_real_range(self%initial%cnh4, 0.0_real64, 1000.0_real64, &
                            'nutrients.initial.cnh4', errors)
      call check_real_range(self%initial%cno3, 0.0_real64, 1000.0_real64, &
                            'nutrients.initial.cno3', errors)
   end subroutine nutrients_config_validate

   subroutine nutrients_config_finalize(self, errors)
      class(nutrients_config_t), intent(inout) :: self
      type(error_collection_t),  intent(inout) :: errors
      ! No finalization needed; included for interface uniformity.
      return
   end subroutine nutrients_config_finalize

end module nutrients_config_mod
