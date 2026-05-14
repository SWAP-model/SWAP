!> @file crop_grass_state.f90
!! SS-GR-CROP: typed crop runtime state — grass-specific runtime.
!! Fields populated in Task A8 via audit of cropgrass_init.f90 and cropgrowth.f90.
module crop_grass_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   public :: crop_grass_state_t

   type :: crop_grass_state_t

      ! Mowing / grazing event schedule
      integer      :: seqgrazmow(366)    = 0          !! sequence grazing/mowing (actual)
      integer      :: seqgrazmowpot(366) = 0          !! sequence grazing/mowing (potential)
      real(real64) :: dateharvest(999)   = 0.0_real64 !! dates of mowing/grazing events
      real(real64) :: mowrest            = 0.0_real64 !! dry weight above-ground after mowing (kg/ha)

      ! Sub-period start/end tracking (reset after each mowing/grazing event)
      real(real64) :: cropstartpot       = 0.0_real64 !! start of potential grass growth (t1900)
      real(real64) :: cropendpot         = 0.0_real64 !! end of potential grass growth (t1900)
      real(real64) :: cropstartact       = 0.0_real64 !! start of actual grass growth (t1900)
      real(real64) :: cropendact         = 0.0_real64 !! end of actual grass growth (t1900)

      ! Management factors
      integer      :: swpotrelmf         = 0          !! calculation of potential yield
      real(real64) :: relmf              = 0.0_real64 !! management factor (attainable yield)

   contains
      procedure :: init => crop_grass_state_init
   end type crop_grass_state_t

contains

   subroutine crop_grass_state_init(self)
      class(crop_grass_state_t), intent(inout) :: self
   end subroutine

end module crop_grass_state_mod
