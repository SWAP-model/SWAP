!> @file crop_wofost_state.f90
!! SS-GR-CROP: typed crop runtime state — WOFOST biomass pools + flows.
!! Fields populated in Task A7 via audit.
module crop_wofost_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   public :: crop_wofost_state_t

   type :: crop_wofost_state_t
      integer :: placeholder_ = 0
   contains
      procedure :: init => crop_wofost_state_init
   end type crop_wofost_state_t

contains

   subroutine crop_wofost_state_init(self)
      class(crop_wofost_state_t), intent(inout) :: self
   end subroutine

end module crop_wofost_state_mod
