!> WOFOST (type-2) crop config — sub-typed per legacy .crp Part 0–15 sectioning.
!! Phase 4c-b Task 1 skeleton: type structure + default initializers + finalize stubs.
!! Validators are added in Task 2.
module cropwofost_config_mod
   use iso_fortran_env, only: real64
   use error_mod, only: error_collection_t
   implicit none
   private

   public :: wofost_preparation_t
   public :: wofost_sowing_t
   public :: wofost_germination_t
   public :: wofost_harvest_t
   public :: wofost_cropfactor_t
   public :: wofost_phenology_t
   public :: wofost_initial_t
   public :: wofost_greenarea_t
   public :: wofost_assimilation_t
   public :: wofost_conversion_t
   public :: wofost_respiration_t
   public :: wofost_partitioning_t
   public :: wofost_death_t
   public :: wofost_root_t
   public :: wofost_oxygen_stress_t
   public :: wofost_drought_stress_t
   public :: wofost_salinity_t
   public :: wofost_compensate_t
   public :: wofost_interception_t
   public :: wofost_co2_t
   public :: wofost_management_t
   public :: cropwofost_config_t

   ! ------------------------------------------------------------------
   ! Preparation
   ! ------------------------------------------------------------------
   type :: wofost_preparation_t
      integer      :: swprep       = 0
      real(real64) :: zprep        = 0.0_real64
      real(real64) :: hprep        = 0.0_real64
      integer      :: maxprepdelay = 0
   contains
      procedure :: finalize => wofost_preparation_finalize
   end type wofost_preparation_t

   ! ------------------------------------------------------------------
   ! Sowing
   ! ------------------------------------------------------------------
   type :: wofost_sowing_t
      integer      :: swsow       = 0
      real(real64) :: zsow        = 0.0_real64
      real(real64) :: hsow        = 0.0_real64
      real(real64) :: ztempsow    = 0.0_real64
      real(real64) :: tempsow     = 0.0_real64
      integer      :: maxsowdelay = 0
   contains
      procedure :: finalize => wofost_sowing_finalize
   end type wofost_sowing_t

   ! ------------------------------------------------------------------
   ! Germination
   ! ------------------------------------------------------------------
   type :: wofost_germination_t
      integer      :: swgerm     = 0
      real(real64) :: tsumemeopt = 0.0_real64
      real(real64) :: tbasem     = 0.0_real64
      real(real64) :: teffmx     = 0.0_real64
      real(real64) :: hdrygerm   = 0.0_real64
      real(real64) :: hwetgerm   = 0.0_real64
      real(real64) :: zgerm      = 0.0_real64
      real(real64) :: agerm      = 0.0_real64
   contains
      procedure :: finalize => wofost_germination_finalize
   end type wofost_germination_t

   ! ------------------------------------------------------------------
   ! Harvest
   ! ------------------------------------------------------------------
   type :: wofost_harvest_t
      real(real64) :: dvsend = 0.0_real64
      integer      :: swharv = 0
   contains
      procedure :: finalize => wofost_harvest_finalize
   end type wofost_harvest_t

   ! ------------------------------------------------------------------
   ! Crop factor
   ! ------------------------------------------------------------------
   type :: wofost_cropfactor_t
      integer      :: swcf   = 0
      real(real64) :: albedo = 0.0_real64
      real(real64) :: rsc    = 0.0_real64
      real(real64) :: rsw    = 0.0_real64
      real(real64), allocatable :: cftb(:,:)
      real(real64), allocatable :: chtb(:,:)
   contains
      procedure :: finalize => wofost_cropfactor_finalize
   end type wofost_cropfactor_t

   ! ------------------------------------------------------------------
   ! Phenology
   ! ------------------------------------------------------------------
   type :: wofost_phenology_t
      integer      :: idsl     = 0
      real(real64) :: tsumea   = 0.0_real64
      real(real64) :: tsumam   = 0.0_real64
      real(real64) :: dlo      = 0.0_real64
      real(real64) :: dlc      = 0.0_real64
      real(real64) :: vernsat  = 0.0_real64
      real(real64) :: vernbase = 0.0_real64
      real(real64) :: verndvs  = 0.0_real64
      real(real64), allocatable :: dtsmtb(:,:)
      real(real64), allocatable :: verntb(:,:)
   contains
      procedure :: finalize => wofost_phenology_finalize
   end type wofost_phenology_t

   ! ------------------------------------------------------------------
   ! Initial
   ! ------------------------------------------------------------------
   type :: wofost_initial_t
      real(real64) :: tdwi   = 0.0_real64
      real(real64) :: laiem  = 0.0_real64
      real(real64) :: rgrlai = 0.0_real64
   contains
      procedure :: finalize => wofost_initial_finalize
   end type wofost_initial_t

   ! ------------------------------------------------------------------
   ! Green area
   ! ------------------------------------------------------------------
   type :: wofost_greenarea_t
      real(real64) :: spa   = 0.0_real64
      real(real64) :: ssa   = 0.0_real64
      real(real64) :: span  = 0.0_real64
      real(real64) :: tbase = 0.0_real64
      real(real64), allocatable :: slatb(:,:)
   contains
      procedure :: finalize => wofost_greenarea_finalize
   end type wofost_greenarea_t

   ! ------------------------------------------------------------------
   ! Assimilation
   ! ------------------------------------------------------------------
   type :: wofost_assimilation_t
      real(real64) :: kdif = 0.0_real64
      real(real64) :: kdir = 0.0_real64
      real(real64) :: eff  = 0.0_real64
      real(real64), allocatable :: amaxtb(:,:)
      real(real64), allocatable :: tmpftb(:,:)
      real(real64), allocatable :: tmnftb(:,:)
   contains
      procedure :: finalize => wofost_assimilation_finalize
   end type wofost_assimilation_t

   ! ------------------------------------------------------------------
   ! Conversion
   ! ------------------------------------------------------------------
   type :: wofost_conversion_t
      real(real64) :: cvl = 0.0_real64
      real(real64) :: cvo = 0.0_real64
      real(real64) :: cvr = 0.0_real64
      real(real64) :: cvs = 0.0_real64
   contains
      procedure :: finalize => wofost_conversion_finalize
   end type wofost_conversion_t

   ! ------------------------------------------------------------------
   ! Respiration
   ! ------------------------------------------------------------------
   type :: wofost_respiration_t
      real(real64) :: q10 = 0.0_real64
      real(real64) :: rml = 0.0_real64
      real(real64) :: rmo = 0.0_real64
      real(real64) :: rmr = 0.0_real64
      real(real64) :: rms = 0.0_real64
      real(real64), allocatable :: rfsetb(:,:)
   contains
      procedure :: finalize => wofost_respiration_finalize
   end type wofost_respiration_t

   ! ------------------------------------------------------------------
   ! Partitioning
   ! ------------------------------------------------------------------
   type :: wofost_partitioning_t
      real(real64), allocatable :: frtb(:,:)
      real(real64), allocatable :: fltb(:,:)
      real(real64), allocatable :: fstb(:,:)
      real(real64), allocatable :: fotb(:,:)
   contains
      procedure :: finalize => wofost_partitioning_finalize
   end type wofost_partitioning_t

   ! ------------------------------------------------------------------
   ! Death
   ! ------------------------------------------------------------------
   type :: wofost_death_t
      real(real64) :: perdl = 0.0_real64
      real(real64), allocatable :: rdrrtb(:,:)
      real(real64), allocatable :: rdrstb(:,:)
   contains
      procedure :: finalize => wofost_death_finalize
   end type wofost_death_t

   ! ------------------------------------------------------------------
   ! Root
   ! ------------------------------------------------------------------
   type :: wofost_root_t
      integer      :: swrd     = 0
      real(real64) :: rdi      = 0.0_real64
      real(real64) :: rri      = 0.0_real64
      real(real64) :: rdc      = 0.0_real64
      integer      :: swdmi2rd = 0
      real(real64) :: wrtmax   = 0.0_real64
      real(real64), allocatable :: rdtb(:,:)
      real(real64), allocatable :: rlwtb(:,:)
      real(real64), allocatable :: rdctb(:,:)
   contains
      procedure :: finalize => wofost_root_finalize
   end type wofost_root_t

   ! ------------------------------------------------------------------
   ! Oxygen stress
   ! ------------------------------------------------------------------
   type :: wofost_oxygen_stress_t
      integer      :: swoxygen               = 0
      integer      :: swwrtnonox             = 0
      real(real64) :: aeratecrit             = 0.0_real64
      real(real64) :: hlim1                  = 0.0_real64
      real(real64) :: hlim2u                 = 0.0_real64
      real(real64) :: hlim2l                 = 0.0_real64
      real(real64) :: q10_microbial          = 0.0_real64
      real(real64) :: specific_resp_humus    = 0.0_real64
      real(real64) :: srl                    = 0.0_real64
      integer      :: swrootradius           = 0
      real(real64) :: dry_mat_cont_roots     = 0.0_real64
      real(real64) :: air_filled_root_por    = 0.0_real64
      real(real64) :: spec_weight_root_tissue = 0.0_real64
      real(real64) :: var_a                  = 0.0_real64
      real(real64) :: root_radiusO2          = 0.0_real64
   contains
      procedure :: finalize => wofost_oxygen_stress_finalize
   end type wofost_oxygen_stress_t

   ! ------------------------------------------------------------------
   ! Drought stress
   ! ------------------------------------------------------------------
   type :: wofost_drought_stress_t
      integer      :: swdrought = 0
      real(real64) :: hlim3h    = 0.0_real64
      real(real64) :: hlim3l    = 0.0_real64
      real(real64) :: hlim4     = 0.0_real64
      real(real64) :: adcrh     = 0.0_real64
      real(real64) :: adcrl     = 0.0_real64
   contains
      procedure :: finalize => wofost_drought_stress_finalize
   end type wofost_drought_stress_t

   ! ------------------------------------------------------------------
   ! Salinity
   ! ------------------------------------------------------------------
   type :: wofost_salinity_t
      integer      :: swsalinity = 0
      real(real64) :: saltmax    = 0.0_real64
      real(real64) :: saltslope  = 0.0_real64
      real(real64) :: salthead   = 0.0_real64
   contains
      procedure :: finalize => wofost_salinity_finalize
   end type wofost_salinity_t

   ! ------------------------------------------------------------------
   ! Compensate
   ! ------------------------------------------------------------------
   type :: wofost_compensate_t
      integer      :: swcompensate = 0
      integer      :: swstressor   = 0
      real(real64) :: alphacrit    = 0.0_real64
      real(real64) :: dcritrtz     = 0.0_real64
   contains
      procedure :: finalize => wofost_compensate_finalize
   end type wofost_compensate_t

   ! ------------------------------------------------------------------
   ! Interception
   ! ------------------------------------------------------------------
   type :: wofost_interception_t
      integer      :: swinter = 0
      real(real64) :: cofab   = 0.0_real64
      real(real64), allocatable :: gashtb(:,:)
   contains
      procedure :: finalize => wofost_interception_finalize
   end type wofost_interception_t

   ! ------------------------------------------------------------------
   ! CO2
   ! ------------------------------------------------------------------
   type :: wofost_co2_t
      integer :: swco2 = 0
      character(len=:), allocatable :: atmofil
      real(real64), allocatable :: co2amaxtb(:,:)
      real(real64), allocatable :: co2efftb(:,:)
      real(real64), allocatable :: co2tratb(:,:)
   contains
      procedure :: finalize => wofost_co2_finalize
   end type wofost_co2_t

   ! ------------------------------------------------------------------
   ! Management
   ! ------------------------------------------------------------------
   type :: wofost_management_t
      real(real64) :: fraharlosorm_lv     = 0.0_real64
      real(real64) :: fraharlosorm_st     = 0.0_real64
      real(real64) :: fraharlosorm_so     = 0.0_real64
      real(real64) :: fradeceasedlvtosoil = 0.0_real64
      logical      :: swpotrelmf          = .false.
      real(real64) :: relmf               = 0.0_real64
   contains
      procedure :: finalize => wofost_management_finalize
   end type wofost_management_t

   ! ------------------------------------------------------------------
   ! Top-level
   ! ------------------------------------------------------------------
   type :: cropwofost_config_t
      type(wofost_preparation_t)     :: preparation
      type(wofost_sowing_t)          :: sowing
      type(wofost_germination_t)     :: germination
      type(wofost_harvest_t)         :: harvest
      type(wofost_cropfactor_t)      :: crop_factor
      type(wofost_phenology_t)       :: phenology
      type(wofost_initial_t)         :: initial
      type(wofost_greenarea_t)       :: green_area
      type(wofost_assimilation_t)    :: assimilation
      type(wofost_conversion_t)      :: conversion
      type(wofost_respiration_t)     :: respiration
      type(wofost_partitioning_t)    :: partitioning
      type(wofost_death_t)           :: death
      type(wofost_root_t)            :: root
      type(wofost_oxygen_stress_t)   :: oxygen_stress
      type(wofost_drought_stress_t)  :: drought_stress
      type(wofost_salinity_t)        :: salinity
      type(wofost_compensate_t)      :: compensate
      type(wofost_interception_t)    :: interception
      type(wofost_co2_t)             :: co2
      type(wofost_management_t)      :: management
   contains
      procedure :: finalize => cropwofost_config_finalize
   end type cropwofost_config_t

contains

   subroutine wofost_preparation_finalize(self, errors)
      class(wofost_preparation_t), intent(inout) :: self
      type(error_collection_t),    intent(inout) :: errors
      return
   end subroutine wofost_preparation_finalize

   subroutine wofost_sowing_finalize(self, errors)
      class(wofost_sowing_t),   intent(inout) :: self
      type(error_collection_t), intent(inout) :: errors
      return
   end subroutine wofost_sowing_finalize

   subroutine wofost_germination_finalize(self, errors)
      class(wofost_germination_t), intent(inout) :: self
      type(error_collection_t),    intent(inout) :: errors
      return
   end subroutine wofost_germination_finalize

   subroutine wofost_harvest_finalize(self, errors)
      class(wofost_harvest_t),  intent(inout) :: self
      type(error_collection_t), intent(inout) :: errors
      return
   end subroutine wofost_harvest_finalize

   subroutine wofost_cropfactor_finalize(self, errors)
      class(wofost_cropfactor_t), intent(inout) :: self
      type(error_collection_t),   intent(inout) :: errors
      return
   end subroutine wofost_cropfactor_finalize

   subroutine wofost_phenology_finalize(self, errors)
      class(wofost_phenology_t), intent(inout) :: self
      type(error_collection_t),  intent(inout) :: errors
      return
   end subroutine wofost_phenology_finalize

   subroutine wofost_initial_finalize(self, errors)
      class(wofost_initial_t),  intent(inout) :: self
      type(error_collection_t), intent(inout) :: errors
      return
   end subroutine wofost_initial_finalize

   subroutine wofost_greenarea_finalize(self, errors)
      class(wofost_greenarea_t), intent(inout) :: self
      type(error_collection_t),  intent(inout) :: errors
      return
   end subroutine wofost_greenarea_finalize

   subroutine wofost_assimilation_finalize(self, errors)
      class(wofost_assimilation_t), intent(inout) :: self
      type(error_collection_t),     intent(inout) :: errors
      return
   end subroutine wofost_assimilation_finalize

   subroutine wofost_conversion_finalize(self, errors)
      class(wofost_conversion_t), intent(inout) :: self
      type(error_collection_t),   intent(inout) :: errors
      return
   end subroutine wofost_conversion_finalize

   subroutine wofost_respiration_finalize(self, errors)
      class(wofost_respiration_t), intent(inout) :: self
      type(error_collection_t),    intent(inout) :: errors
      return
   end subroutine wofost_respiration_finalize

   subroutine wofost_partitioning_finalize(self, errors)
      class(wofost_partitioning_t), intent(inout) :: self
      type(error_collection_t),     intent(inout) :: errors
      return
   end subroutine wofost_partitioning_finalize

   subroutine wofost_death_finalize(self, errors)
      class(wofost_death_t),    intent(inout) :: self
      type(error_collection_t), intent(inout) :: errors
      return
   end subroutine wofost_death_finalize

   subroutine wofost_root_finalize(self, errors)
      class(wofost_root_t),     intent(inout) :: self
      type(error_collection_t), intent(inout) :: errors
      return
   end subroutine wofost_root_finalize

   subroutine wofost_oxygen_stress_finalize(self, errors)
      class(wofost_oxygen_stress_t), intent(inout) :: self
      type(error_collection_t),      intent(inout) :: errors
      return
   end subroutine wofost_oxygen_stress_finalize

   subroutine wofost_drought_stress_finalize(self, errors)
      class(wofost_drought_stress_t), intent(inout) :: self
      type(error_collection_t),       intent(inout) :: errors
      return
   end subroutine wofost_drought_stress_finalize

   subroutine wofost_salinity_finalize(self, errors)
      class(wofost_salinity_t), intent(inout) :: self
      type(error_collection_t), intent(inout) :: errors
      return
   end subroutine wofost_salinity_finalize

   subroutine wofost_compensate_finalize(self, errors)
      class(wofost_compensate_t), intent(inout) :: self
      type(error_collection_t),   intent(inout) :: errors
      return
   end subroutine wofost_compensate_finalize

   subroutine wofost_interception_finalize(self, errors)
      class(wofost_interception_t), intent(inout) :: self
      type(error_collection_t),     intent(inout) :: errors
      return
   end subroutine wofost_interception_finalize

   subroutine wofost_co2_finalize(self, errors)
      class(wofost_co2_t),      intent(inout) :: self
      type(error_collection_t), intent(inout) :: errors
      return
   end subroutine wofost_co2_finalize

   subroutine wofost_management_finalize(self, errors)
      class(wofost_management_t), intent(inout) :: self
      type(error_collection_t),   intent(inout) :: errors
      return
   end subroutine wofost_management_finalize

   subroutine cropwofost_config_finalize(self, errors)
      class(cropwofost_config_t), intent(inout) :: self
      type(error_collection_t),   intent(inout) :: errors

      call self%preparation%finalize(errors)
      call self%sowing%finalize(errors)
      call self%germination%finalize(errors)
      call self%harvest%finalize(errors)
      call self%crop_factor%finalize(errors)
      call self%phenology%finalize(errors)
      call self%initial%finalize(errors)
      call self%green_area%finalize(errors)
      call self%assimilation%finalize(errors)
      call self%conversion%finalize(errors)
      call self%respiration%finalize(errors)
      call self%partitioning%finalize(errors)
      call self%death%finalize(errors)
      call self%root%finalize(errors)
      call self%oxygen_stress%finalize(errors)
      call self%drought_stress%finalize(errors)
      call self%salinity%finalize(errors)
      call self%compensate%finalize(errors)
      call self%interception%finalize(errors)
      call self%co2%finalize(errors)
      call self%management%finalize(errors)
   end subroutine cropwofost_config_finalize

end module cropwofost_config_mod
