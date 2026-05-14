!> @file crop_wofost_state.f90
!! SS-GR-CROP: typed crop runtime state — WOFOST biomass pools + flows.
!! Fields populated in Task A7 via audit of cropgrowth.f90 and cropwofost_init.f90.
module crop_wofost_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   public :: crop_wofost_state_t

   type :: crop_wofost_state_t

      ! Biomass pools: actual + potential (kg/ha DM)
      real(real64) :: wlv     = 0.0_real64   !! dry weight plant leaves (actual)
      real(real64) :: wlvpot  = 0.0_real64   !! dry weight plant leaves (potential)
      real(real64) :: wst     = 0.0_real64   !! dry weight plant stem (actual)
      real(real64) :: wstpot  = 0.0_real64   !! dry weight plant stem (potential)
      real(real64) :: wrt     = 0.0_real64   !! dry weight plant root (actual)
      real(real64) :: wrtpot  = 0.0_real64   !! dry weight plant root (potential)
      real(real64) :: wso     = 0.0_real64   !! dry weight storage organ (actual)
      real(real64) :: wsopot  = 0.0_real64   !! dry weight storage organ (potential)

      ! Aggregate biomass (kg/ha DM)
      real(real64) :: tagp    = 0.0_real64   !! dry weight dead+living grass organs (actual)
      real(real64) :: tagppot = 0.0_real64   !! dry weight dead+living grass organs (potential)
      real(real64) :: tagpt   = 0.0_real64   !! dry weight harvested grass (actual)
      real(real64) :: tagptpot = 0.0_real64  !! dry weight harvested grass (potential)
      real(real64) :: cwdm    = 0.0_real64   !! dry weight dead+living plant organs (actual)
      real(real64) :: cwdmpot = 0.0_real64   !! dry weight dead+living plant organs (potential)

      ! Assimilation rates
      real(real64) :: pgass   = 0.0_real64   !! assimilation rate (actual, kg/ha/d)
      real(real64) :: pgasspot = 0.0_real64  !! assimilation rate (potential, kg/ha/d)

      ! Death / decay flows (kg/ha/d)
      real(real64) :: dwlv    = 0.0_real64   !! death rate leaves (actual)
      real(real64) :: dwlvpot = 0.0_real64   !! death rate leaves (potential)
      real(real64) :: dwst    = 0.0_real64   !! death rate stem (actual)
      real(real64) :: dwstpot = 0.0_real64   !! death rate stem (potential)
      real(real64) :: dwrt    = 0.0_real64   !! death rate root (actual)
      real(real64) :: dwrtpot = 0.0_real64   !! death rate root (potential)
      real(real64) :: dwso    = 0.0_real64   !! death rate storage organ (actual)
      real(real64) :: dwlvCrop = 0.0_real64  !! deceased leaves remaining on plant (actual)
      real(real64) :: dwlvSoil = 0.0_real64  !! deceased leaves allocated to soil (actual)

      ! Harvest losses (kg/ha)
      real(real64) :: plossdm = 0.0_real64   !! total loss potential harvest (insufficient h)
      real(real64) :: lossdm  = 0.0_real64   !! total loss actual harvest (insufficient h)

      ! Bulb-crop fields (swbulb=1: tulips etc.)
      logical      :: swbulb  = .false.      !! enable bulb crop simulation
      real(real64) :: plwt    = 0.0_real64   !! dry weight mother bulb (kg/ha)
      real(real64) :: plwti   = 0.0_real64   !! initial dry weight mother bulb (kg/ha)
      real(real64) :: wbl     = 0.0_real64   !! dry weight living flowers (actual, kg/ha)
      real(real64) :: wblpot  = 0.0_real64   !! dry weight living flowers (potential, kg/ha)
      real(real64) :: dwbl    = 0.0_real64   !! dry weight dead flowers (actual, kg/ha)
      real(real64) :: dwblpot = 0.0_real64   !! dry weight dead flowers (potential, kg/ha)

      ! [SS-GR-CROPRT A4] FCO2 derived correction factors (default 1.0 = no CO2 correction)
      real(real64) :: fco2amax = 1.0_real64  !! CO2 correction factor for AMAX (-)
      real(real64) :: fco2eff  = 1.0_real64  !! CO2 correction factor for EFF (-)
      real(real64) :: fco2tra  = 1.0_real64  !! CO2 correction factor for TRA (-)

   contains
      procedure :: init => crop_wofost_state_init
   end type crop_wofost_state_t

contains

   subroutine crop_wofost_state_init(self)
      class(crop_wofost_state_t), intent(inout) :: self
   end subroutine

end module crop_wofost_state_mod
