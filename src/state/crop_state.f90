!> @file crop_state.f90
!! SS-GR-ATM: typed crop runtime state. Foundation for Arc 8 (crop
!! cluster). Arc 4 (atmosphere) introduces this module to host the
!! ~12 crop-runtime symbols read by atmosphere code; Arc 8 later
!! migrates all remaining crop readers.
module crop_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   use crop_common_state_mod, only: crop_common_state_t
   use crop_fixed_state_mod,   only: crop_fixed_state_t
   use crop_wofost_state_mod,  only: crop_wofost_state_t
   implicit none
   private
   public :: crop_state_t

   type :: crop_state_t
      real(real64) :: lai             = 0.0_real64    !! leaf area index (-)
      real(real64) :: kdif            = 0.0_real64    !! extinction for diffuse light (-)
      real(real64) :: kdir            = 0.0_real64    !! extinction for direct light (-)
      real(real64) :: cofab           = 0.0_real64    !! interception coefficient (cm)
      real(real64) :: cfbs            = 1.0_real64    !! bare-soil ET factor (-)
      integer      :: swcf            = 0             !! crop factor switch
      integer      :: swcfbs          = 0             !! bare-soil factor switch
      real(real64) :: gird            = 0.0_real64    !! gross irrigation depth (cm)
      logical      :: flCropEmergence = .false.       !! crop emerged flag
      real(real64) :: et0             = 0.0_real64    !! potential ET (cm/d)
      real(real64) :: ew0             = 0.0_real64    !! potential evap wet crop (cm/d)
      real(real64) :: es0             = 0.0_real64    !! potential evap bare soil (cm/d)
      ! [SS-GR-CROP A2] sub-record for shared crop runtime fields
      type(crop_common_state_t) :: common
      ! [SS-GR-CROP A3] sub-record for fixed-crop runtime fields
      type(crop_fixed_state_t) :: fixed
      ! [SS-GR-CROP A4] sub-record for WOFOST biomass pools + flows
      type(crop_wofost_state_t) :: wofost
   contains
      procedure :: init => crop_state_init
   end type crop_state_t

contains

   subroutine crop_state_init(self, crop_cfg)
      use crop_config_mod, only: crop_config_t
      class(crop_state_t),  intent(inout) :: self
      type(crop_config_t),  intent(in)    :: crop_cfg
      ! Defaults retained from type initializers; runtime seeders update
      ! each timestep. crop_cfg arg reserved for future seed migration.
   end subroutine crop_state_init

end module crop_state_mod
