!> @file crop_common_state.f90
!! SS-GR-CROP: typed crop runtime state — shared across all crop types.
!! Fields populated in Task A6 via audit.
module crop_common_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   public :: crop_common_state_t

   type :: crop_common_state_t
      integer :: placeholder_ = 0
   contains
      procedure :: init => crop_common_state_init
   end type crop_common_state_t

contains

   subroutine crop_common_state_init(self)
      class(crop_common_state_t), intent(inout) :: self
   end subroutine

end module crop_common_state_mod
