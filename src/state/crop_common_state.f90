!> @file crop_common_state.f90
!! SS-GR-CROP: typed crop runtime state — shared across all crop types.
!! Fields populated in Task A6 via audit of cropgrowth.f90 and init files.
module crop_common_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   public :: crop_common_state_t

   type :: crop_common_state_t

      ! Development & calendar
      integer      :: daycrop        = 0              !! days since crop start
      real(real64) :: dvs            = 0.0_real64     !! development stage (-)
      real(real64) :: tsum           = 0.0_real64     !! temperature sum (degC)
      integer      :: swcrp          = 0              !! crop type: 1=fixed, 2=wofost, 3=grass
      integer      :: icrop          = 0              !! current crop number
      logical      :: flCropCalendar = .false.        !! crop season is active
      logical      :: flCropOutput   = .false.        !! write *.CRP output
      logical      :: flCropNut      = .false.        !! simulate crop nutrient stress
      logical      :: flHarvestDay   = .false.        !! current day is harvest day

      ! End-of-simulation control (integer switch in legacy: 0/1/2)
      integer      :: swend          = 0

      ! Root depth
      real(real64) :: rd             = 0.0_real64     !! actual rooting depth (L)
      real(real64) :: rdpot          = 0.0_real64     !! rooting depth potential run (L)
      real(real64) :: rdm            = 0.0_real64     !! max rooting depth (min soil/crop) (L)
      real(real64) :: rri            = 0.0_real64     !! max daily root depth increase (L/T)
      real(real64) :: rdi            = 0.0_real64     !! initial rooting depth (L)
      real(real64) :: rdc            = 0.0_real64     !! max crop rooting depth (L)

      ! Crop physiology summary
      real(real64) :: ch             = 0.0_real64     !! crop height (cm)
      real(real64) :: cf             = 0.0_real64     !! crop factor (-)
      real(real64) :: laipot         = 0.0_real64     !! leaf area index potential run (-)

      ! Grazing / harvest totals
      real(real64) :: cuptgraz       = 0.0_real64     !! cumul. dry weight grazed actual (kg/ha)
      real(real64) :: cuptgrazpot    = 0.0_real64     !! cumul. dry weight grazed potential (kg/ha)
      real(real64) :: HarLosOrm_tot  = 0.0_real64     !! harvest losses to soil at harvest (kg/ha DM)

   contains
      procedure :: init => crop_common_state_init
   end type crop_common_state_t

contains

   subroutine crop_common_state_init(self)
      class(crop_common_state_t), intent(inout) :: self
   end subroutine

end module crop_common_state_mod
