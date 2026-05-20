!> @file crop_fixed_state.f90
!! SS-GR-CROP: typed crop runtime state — fixed-crop runtime.
!! Fields populated in Task A9 via audit of cropfixed_init.f90 and cropgrowth.f90.
module crop_fixed_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   use swap_array_dimensions, only: MAGRS
   implicit none
   private
   public :: crop_fixed_state_t

   type :: crop_fixed_state_t

      ! Lookup tables (populated from config at init; indexed by DVS or LAI)
      real(real64) :: cftb(2*MAGRS)    = 0.0_real64  !! crop factor or height vs. DVS (2*366)
      real(real64) :: chtb(2*MAGRS)    = 0.0_real64  !! crop height (cm) vs. DVS (2*366)
      real(real64) :: cfeictb(2*MAGRS) = 0.0_real64  !! crop factor wet (-) vs. DVS (2*366)
      real(real64) :: gctb(2*MAGRS)    = 0.0_real64  !! LAI or soil-cover fraction vs. DVS (2*366)

      ! Runtime scalar derived from table lookups each timestep
      real(real64) :: cfeic            = 0.0_real64  !! crop factor wet (-) current timestep

   contains
      procedure :: init => crop_fixed_state_init
   end type crop_fixed_state_t

contains

   subroutine crop_fixed_state_init(self)
      class(crop_fixed_state_t), intent(inout) :: self
   end subroutine

end module crop_fixed_state_mod
