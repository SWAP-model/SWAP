!> @file crop_fixed_state.f90
!! SS-GR-CROP: typed crop runtime state — fixed-crop runtime.
!! Fields populated in Task A9 via audit.
module crop_fixed_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   public :: crop_fixed_state_t

   type :: crop_fixed_state_t
      integer :: placeholder_ = 0
   contains
      procedure :: init => crop_fixed_state_init
   end type crop_fixed_state_t

contains

   subroutine crop_fixed_state_init(self)
      class(crop_fixed_state_t), intent(inout) :: self
   end subroutine

end module crop_fixed_state_mod
