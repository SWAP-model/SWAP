!> @file crop_common_state.f90
!! SS-GR-CROP: typed crop runtime state — shared across all crop types.
!! Fields populated in Task A6 via audit of cropgrowth.f90 and init files.
module crop_common_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   use swap_array_dimensions, only: MAGRS
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
      integer      :: swinter        = 0              !! interception switch (0=none, 1=ag crops, 2=trees, 3=Gash)
      real(real64) :: saltmax        = 0.0_real64     !! threshold conc above which yield reduces (Maas-Hoffman)
      real(real64) :: saltslope      = 0.0_real64     !! slope of yield reduction vs concentration (Maas-Hoffman)
      real(real64) :: salthead       = 0.0_real64     !! osmotic head per concentration unit
      real(real64) :: hlim3h         = 0.0_real64     !! pressure head at low atmospheric demand (cm, Feddes)
      real(real64) :: hlim3l         = 0.0_real64     !! pressure head at high atmospheric demand (cm, Feddes)
      real(real64) :: hlim4          = 0.0_real64     !! wilting point pressure head (cm, Feddes)
      real(real64) :: adcrh          = 0.0_real64     !! high atmospheric demand threshold (cm/d)
      real(real64) :: adcrl          = 0.0_real64     !! low atmospheric demand threshold (cm/d)
      real(real64) :: hlim1          = 0.0_real64     !! pressure head above which uptake stops (cm, anaerobic)
      real(real64) :: hlim2u         = 0.0_real64     !! pressure head — optimum uptake starts, top layer (cm)
      real(real64) :: hlim2l         = 0.0_real64     !! pressure head — optimum uptake starts, sub layer (cm)
      integer      :: swWrtNonox     = 0              !! switch for crop survival under non-aerated conditions
      real(real64) :: aeratecrit     = 0.0_real64     !! critical aeration level for crop survival
      integer      :: swstressor     = 0              !! compensation stressor switch (1=all, 2=drought, 3=oxy, 4=salt, 5=frost)
      real(real64) :: alphacrit      = 0.0_real64     !! critical stress index for compensation (mutated at runtime when swcompensate==2)
      real(real64) :: dcritrtz       = 0.0_real64     !! root-zone depth threshold for Walsum compensation (cm)
      integer      :: schedule       = 0              !! per-crop irrigation scheduling switch (0=fixed, 1=scheduled)
      real(real64) :: rsw            = 0.0_real64     !! canopy resistance to intercepted water (s/m)

      ! Root tables (per-rotation config copied at init)
      real(real64) :: rdtb(2*MAGRS)  = 0.0_real64     !! root depth vs DVS table (pair, 2*MAGRS)
      real(real64) :: rdctb(22)      = 0.0_real64     !! relative root density vs relative depth (22-elem)
      real(real64) :: rlwtb(22)      = 0.0_real64     !! root depth vs root biomass (22-elem pair)
      real(real64) :: wrtmax         = 0.0_real64     !! maximum root weight (kg/ha)

      ! Physiology lookup tables (30-elem) + light use efficiency scalars
      real(real64) :: slatb(30)      = 0.0_real64     !! specific leaf area vs DVS (ha/kg)
      real(real64) :: amaxtb(30)     = 0.0_real64     !! max CO2 assimilation rate vs DVS (kg/ha/hr)
      real(real64) :: tmpftb(30)     = 0.0_real64     !! reduction factor vs daily-mean temp
      real(real64) :: tmnftb(30)     = 0.0_real64     !! reduction factor vs min daily temp
      real(real64) :: eff            = 0.0_real64     !! light use efficiency of a leaf (kg CO2 / J adsorbed)
      real(real64) :: rgrlai         = 0.0_real64     !! max relative increase in LAI (1/d)

      ! Biomass conversion and maintenance respiration (scalars + table)
      real(real64) :: cvl            = 0.0_real64     !! efficiency assimilates→leaves (kg/kg)
      real(real64) :: cvr            = 0.0_real64     !! efficiency assimilates→roots (kg/kg)
      real(real64) :: cvs            = 0.0_real64     !! efficiency assimilates→stems (kg/kg)
      real(real64) :: cvo            = 0.0_real64     !! efficiency assimilates→storage organs (kg/kg)
      real(real64) :: q10            = 0.0_real64     !! Q10 temperature factor for respiration
      real(real64) :: rmr            = 0.0_real64     !! rel. maint. respiration rate, roots (kg/kg/d)
      real(real64) :: rml            = 0.0_real64     !! rel. maint. respiration rate, leaves
      real(real64) :: rms            = 0.0_real64     !! rel. maint. respiration rate, stems
      real(real64) :: rmo            = 0.0_real64     !! rel. maint. respiration rate, storage organs
      real(real64) :: rfsetb(30)     = 0.0_real64     !! senescence-effect factor vs DVS (30-elem)

      ! Partitioning + development + leaf-death tables (30-elem each)
      real(real64) :: frtb(30)       = 0.0_real64     !! root fraction vs DVS
      real(real64) :: fltb(30)       = 0.0_real64     !! leaf fraction vs DVS
      real(real64) :: fstb(30)       = 0.0_real64     !! stem fraction vs DVS
      real(real64) :: fotb(30)       = 0.0_real64     !! storage-organ fraction vs DVS
      real(real64) :: fbltb(30)      = 0.0_real64     !! bulb fraction vs DVS (swbulb=1 only)
      real(real64) :: dtsmtb(30)     = 0.0_real64     !! daily temp-sum increment vs Tavg
      real(real64) :: rdrrtb(30)     = 0.0_real64     !! relative death rate of roots vs DVS
      real(real64) :: rdrstb(30)     = 0.0_real64     !! relative death rate of stems vs DVS

      ! Leaf/stem/pod area + senescence scalars
      real(real64) :: ssa            = 0.0_real64     !! specific stem area (ha/kg)
      real(real64) :: spa            = 0.0_real64     !! specific pod area (ha/kg)
      real(real64) :: span           = 0.0_real64     !! life span of leaves at optimum (d)
      real(real64) :: perdl          = 0.0_real64     !! max relative leaf death rate due to water stress (1/d)

      ! Bartholomeus oxygen-stress physical-model config (per rotation)
      integer      :: swoxygentype          = 0          !! 1=physical, 2=reproduction functions
      integer      :: swrootradius          = 0          !! 1=calculated, 2=given
      real(real64) :: srl                   = 0.0_real64 !! specific root length (cm/g)
      real(real64) :: dry_mat_cont_roots    = 0.0_real64 !! dry matter content of roots
      real(real64) :: air_filled_root_por   = 0.0_real64 !! air-filled root porosity
      real(real64) :: spec_weight_root_tissue = 0.0_real64 !! specific weight root tissue
      real(real64) :: var_a                 = 0.0_real64 !! Carsel-Parrish variability constant
      real(real64) :: root_radiusO2         = 0.0_real64 !! given root radius (m)

      ! Leaf-area initial + runtime scalars (per-rotation; runtime-mutated)
      real(real64) :: tdwi          = 0.0_real64       !! initial total dry weight (kg/ha)
      real(real64) :: laiem         = 0.0_real64       !! LAI at emergence (m2/m2)
      real(real64) :: laiexp        = 0.0_real64       !! current exponential-phase LAI (actual)
      real(real64) :: laiexppot     = 0.0_real64       !! current exponential-phase LAI (potential)
      real(real64) :: laimax        = 0.0_real64       !! max LAI achieved (m2/m2)
      real(real64) :: glaiex        = 0.0_real64       !! daily LAI increase, exponential (actual)
      real(real64) :: glaiexpot     = 0.0_real64       !! daily LAI increase, exponential (potential)

      ! Leaf cohort arrays (per crop day; 366 entries each, runtime-mutated)
      real(real64) :: lv(366)       = 0.0_real64       !! leaf weight by cohort day (actual, kg/ha)
      real(real64) :: lvpot(366)    = 0.0_real64       !! leaf weight by cohort day (potential, kg/ha)
      real(real64) :: lvage(366)    = 0.0_real64       !! leaf age by cohort day (actual, d)
      real(real64) :: lvagepot(366) = 0.0_real64       !! leaf age by cohort day (potential, d)
      real(real64) :: sla(366)      = 0.0_real64       !! specific leaf area by cohort day (actual)
      real(real64) :: slapot(366)   = 0.0_real64       !! specific leaf area by cohort day (potential)
      integer      :: ilvold        = 0                !! oldest-leaf day index (actual)
      integer      :: ilvoldpot     = 0                !! oldest-leaf day index (potential)

   contains
      procedure :: init => crop_common_state_init
   end type crop_common_state_t

contains

   subroutine crop_common_state_init(self)
      class(crop_common_state_t), intent(inout) :: self
      self%cumdens = 0.0_real64
   end subroutine

end module crop_common_state_mod
