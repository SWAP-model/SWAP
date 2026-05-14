!> @file nutrients_state.f90
!! SS-GR-CROP: typed nutrients runtime state. Hosts management_soil +
!! wofost_soil_* nutrient pool and flow data. Fields populated in
!! Task A11 via audit.
module nutrients_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   public :: nutrients_state_t

   type :: nutrients_state_t
      integer :: placeholder_ = 0
   contains
      procedure :: init => nutrients_state_init
   end type nutrients_state_t

contains

   subroutine nutrients_state_init(self, nlay)
      class(nutrients_state_t), intent(inout) :: self
      integer, intent(in) :: nlay
      ! Allocation/seeding deferred to Task A11.
   end subroutine

end module nutrients_state_mod
