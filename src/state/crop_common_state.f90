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

      ! [SS-GR-CROPRT D1] swend state field retired — ADR 0009 zeroed permanently; swap_mod branches dropped

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

      ! [SS-GR-CROPRT A4] lifecycle flags
      logical      :: flCropReadFile = .false.        !! reading of input .crp file occurred
      logical      :: flCropPrep     = .false.        !! ploughing opportunity realized
      logical      :: flCropSow      = .false.        !! sowing opportunity realized
      logical      :: flCropGerm     = .false.        !! germination realized
      logical      :: flCropHarvest  = .false.        !! harvest event day flag

      ! [SS-GR-CROPRT A4] crop window dates (current crop's start/end, days since 1900)
      real(real64) :: cropstart      = 0.0_real64     !! current crop season start (d)
      real(real64) :: cropend        = 0.0_real64     !! current crop season end (d)

      ! [SS-GR-CROPRT A4] germination delay scratch
      integer      :: PrepDelay      = 0              !! delay of preparation (d)
      integer      :: SowDelay       = 0              !! delay of sowing (d)

      ! [SS-GR-CROPRT A4] runtime root zone
      integer      :: noddrz         = 0              !! compartment number at bottom root zone (-)
      real(real64) :: cumdens(202)   = 0.0_real64     !! cumul. root density as fn of rel. soil depth (-)

      ! [SS-GR-CROPRT A4] surface params
      real(real64) :: albedo         = 0.0_real64     !! crop reflection coefficient (-)
      real(real64) :: rsc            = 0.0_real64     !! minimum canopy resistance dry crop (T/L)

      ! [Sweep 2] per-rotation switches mirrored from active rotation's cfg
      integer      :: swrd           = 0              !! root depth method (1=table,2=daily,3=biomass)
      integer      :: swdmi2rd       = 0              !! transpiration limit on root depth (0/1)
      integer      :: swrdc          = 0              !! root density input (0/1)
      real(real64) :: tbase          = 0.0_real64     !! base temperature for ageing of leaves (°C)
      integer      :: idev           = 0              !! length-of-growth-period switch: 1=fixed, 2=tsum
      real(real64) :: tsumea         = 0.0_real64     !! temperature sum emergence→anthesis (°C·d)
      real(real64) :: tsumam         = 0.0_real64     !! temperature sum anthesis→maturity (°C·d)
      real(real64) :: dvsend         = 0.0_real64     !! crop development stage at harvest (-)
      integer      :: swharv         = 0              !! harvest timing switch (0=cropend, 1=maturity)
      integer      :: swgc           = 0              !! green-canopy switch (1=LAI input, 2=soil cover fraction)
      integer      :: swsalinity     = 0              !! salinity stress switch (0=none, 1=Maas-Hoffman, 2=osmotic)
      integer      :: swcompensate   = 0              !! root water uptake compensation method switch
      integer      :: swoxygen       = 0              !! oxygen stress switch (1=Feddes, 2=Bartholomeus)
      integer      :: swdrought      = 0              !! drought stress switch (1=Feddes, 2=De Jong van Lier)

   contains
      procedure :: init => crop_common_state_init
   end type crop_common_state_t

contains

   subroutine crop_common_state_init(self)
      class(crop_common_state_t), intent(inout) :: self
      self%cumdens = 0.0_real64
   end subroutine

end module crop_common_state_mod
