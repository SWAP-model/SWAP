!> WOFOST (type-2) crop config — sub-typed per legacy .crp Part 0–15 sectioning.
!! Phase 4c-b Task 1 skeleton: type structure + default initializers + finalize stubs.
!! Phase 4c-b Task 2: per-sub-type validators wired into top-level validate.
module cropwofost_config_mod
   use iso_fortran_env, only: real64
   use error_mod, only: error_collection_t, ERR_VALIDATION_OUT_OF_RANGE
   use validation_mod, only: check_int_enum, check_int_range, check_real_range, &
                             check_nonnegative_real, check_ordered_pair
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
      procedure :: validate => wofost_preparation_validate
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
      procedure :: validate => wofost_sowing_validate
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
      procedure :: validate => wofost_germination_validate
      procedure :: finalize => wofost_germination_finalize
   end type wofost_germination_t

   ! ------------------------------------------------------------------
   ! Harvest
   ! ------------------------------------------------------------------
   type :: wofost_harvest_t
      real(real64) :: dvsend = 0.0_real64
      integer      :: swharv = 0
   contains
      procedure :: validate => wofost_harvest_validate
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
      procedure :: validate => wofost_cropfactor_validate
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
      procedure :: validate => wofost_phenology_validate
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
      procedure :: validate => wofost_initial_validate
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
      procedure :: validate => wofost_greenarea_validate
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
      procedure :: validate => wofost_assimilation_validate
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
      procedure :: validate => wofost_conversion_validate
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
      procedure :: validate => wofost_respiration_validate
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
      procedure :: validate => wofost_partitioning_validate
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
      procedure :: validate => wofost_death_validate
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
      procedure :: validate => wofost_root_validate
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
      procedure :: validate => wofost_oxygen_stress_validate
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
      procedure :: validate => wofost_drought_stress_validate
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
      procedure :: validate => wofost_salinity_validate
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
      procedure :: validate => wofost_compensate_validate
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
      procedure :: validate => wofost_interception_validate
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
      procedure :: validate => wofost_co2_validate
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
      integer      :: swpotrelmf          = 0
      real(real64) :: relmf               = 0.0_real64
   contains
      procedure :: validate => wofost_management_validate
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
      procedure :: validate => cropwofost_config_validate
      procedure :: finalize => cropwofost_config_finalize
   end type cropwofost_config_t

contains

   !> Verify table has expected ncols, 2..15 rows, and (optionally)
   !! a strictly increasing first column. Skips entirely when unallocated.
   subroutine check_table(table, expected_cols, monotone, label, errors)
      real(real64), allocatable, intent(in)    :: table(:,:)
      integer,                   intent(in)    :: expected_cols
      logical,                   intent(in)    :: monotone
      character(len=*),          intent(in)    :: label
      type(error_collection_t),  intent(inout) :: errors
      character(len=128) :: msg
      integer :: nrows, ncols, i

      if (.not. allocated(table)) return

      nrows = size(table, 1)
      ncols = size(table, 2)

      if (ncols /= expected_cols) then
         write(msg, '("ncols=",I0," expected ",I0)') ncols, expected_cols
         call errors%append(ERR_VALIDATION_OUT_OF_RANGE, trim(msg), label)
      end if
      if (nrows < 2) then
         call errors%append(ERR_VALIDATION_OUT_OF_RANGE, "nrows<2", label)
      end if
      if (nrows > 15) then
         call errors%append(ERR_VALIDATION_OUT_OF_RANGE, "nrows>15", label)
      end if
      if (monotone .and. nrows >= 2 .and. ncols >= 1) then
         do i = 2, nrows
            if (table(i, 1) <= table(i-1, 1)) then
               call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
                                  "first column not strictly increasing", label)
               exit
            end if
         end do
      end if
   end subroutine check_table

   subroutine wofost_preparation_validate(self, errors)
      class(wofost_preparation_t), intent(in)    :: self
      type(error_collection_t),    intent(inout) :: errors
      call check_int_enum(self%swprep, [0, 1], 'wofost.preparation.swprep', errors)
      if (self%swprep == 1) then
         call check_real_range(self%zprep, -100.0_real64,    0.0_real64, 'wofost.preparation.zprep', errors)
         call check_real_range(self%hprep, -200.0_real64,    0.0_real64, 'wofost.preparation.hprep', errors)
         call check_int_range(self%maxprepdelay, 1, 366, 'wofost.preparation.maxprepdelay', errors)
      end if
   end subroutine wofost_preparation_validate

   subroutine wofost_sowing_validate(self, errors)
      class(wofost_sowing_t),   intent(in)    :: self
      type(error_collection_t), intent(inout) :: errors
      call check_int_enum(self%swsow, [0, 1], 'wofost.sowing.swsow', errors)
      if (self%swsow == 1) then
         call check_real_range(self%zsow,     -100.0_real64,   0.0_real64, 'wofost.sowing.zsow',     errors)
         call check_real_range(self%hsow,     -200.0_real64,   0.0_real64, 'wofost.sowing.hsow',     errors)
         call check_real_range(self%ztempsow, -100.0_real64,   0.0_real64, 'wofost.sowing.ztempsow', errors)
         call check_real_range(self%tempsow,    0.0_real64,   30.0_real64, 'wofost.sowing.tempsow',  errors)
         call check_int_range(self%maxsowdelay, 1, 366, 'wofost.sowing.maxsowdelay', errors)
      end if
   end subroutine wofost_sowing_validate

   subroutine wofost_germination_validate(self, errors)
      class(wofost_germination_t), intent(in)    :: self
      type(error_collection_t),    intent(inout) :: errors
      call check_int_enum(self%swgerm, [0, 1, 2], 'wofost.germination.swgerm', errors)
      if (self%swgerm == 1 .or. self%swgerm == 2) then
         call check_real_range(self%tsumemeopt, 0.0_real64, 1000.0_real64, 'wofost.germination.tsumemeopt', errors)
         call check_real_range(self%tbasem,     0.0_real64,   40.0_real64, 'wofost.germination.tbasem',     errors)
         call check_real_range(self%teffmx,     0.0_real64,   40.0_real64, 'wofost.germination.teffmx',     errors)
      end if
      if (self%swgerm == 2) then
         call check_real_range(self%hdrygerm, -1000.0_real64, -0.01_real64, 'wofost.germination.hdrygerm', errors)
         call check_real_range(self%hwetgerm,  -100.0_real64, -0.01_real64, 'wofost.germination.hwetgerm', errors)
         call check_real_range(self%zgerm,     -100.0_real64,   0.0_real64, 'wofost.germination.zgerm',    errors)
         call check_real_range(self%agerm,        1.0_real64, 1000.0_real64, 'wofost.germination.agerm',   errors)
      end if
   end subroutine wofost_germination_validate

   subroutine wofost_harvest_validate(self, errors)
      class(wofost_harvest_t),  intent(in)    :: self
      type(error_collection_t), intent(inout) :: errors
      call check_real_range(self%dvsend, 0.0_real64, 3.0_real64, 'wofost.harvest.dvsend', errors)
      call check_int_enum(self%swharv,   [0, 1], 'wofost.harvest.swharv', errors)
   end subroutine wofost_harvest_validate

   subroutine wofost_cropfactor_validate(self, errors)
      class(wofost_cropfactor_t), intent(in)    :: self
      type(error_collection_t),   intent(inout) :: errors
      call check_int_enum(self%swcf, [1, 2], 'wofost.crop_factor.swcf', errors)
      if (self%swcf == 2) then
         call check_real_range(self%albedo, 0.0_real64,        1.0_real64, 'wofost.crop_factor.albedo', errors)
         call check_real_range(self%rsc,    0.0_real64,    1.0e6_real64, 'wofost.crop_factor.rsc',    errors)
         call check_real_range(self%rsw,    0.0_real64,    1.0e6_real64, 'wofost.crop_factor.rsw',    errors)
      end if
      if (self%swcf == 1) then
         call check_table(self%cftb, 2, .true., 'wofost.crop_factor.cftb', errors)
      else
         call check_table(self%chtb, 2, .true., 'wofost.crop_factor.chtb', errors)
      end if
   end subroutine wofost_cropfactor_validate

   subroutine wofost_phenology_validate(self, errors)
      class(wofost_phenology_t), intent(in)    :: self
      type(error_collection_t),  intent(inout) :: errors
      call check_int_enum(self%idsl, [0, 1, 2], 'wofost.phenology.idsl', errors)
      call check_real_range(self%tsumea, 0.0_real64, 10000.0_real64, 'wofost.phenology.tsumea', errors)
      call check_real_range(self%tsumam, 0.0_real64, 10000.0_real64, 'wofost.phenology.tsumam', errors)
      call check_table(self%dtsmtb, 2, .true., 'wofost.phenology.dtsmtb', errors)
      if (self%idsl == 1 .or. self%idsl == 2) then
         call check_real_range(self%dlo, 0.0_real64, 24.0_real64, 'wofost.phenology.dlo', errors)
         call check_real_range(self%dlc, 0.0_real64, 24.0_real64, 'wofost.phenology.dlc', errors)
      end if
      if (self%idsl == 2) then
         call check_real_range(self%vernsat,  0.0_real64, 100.0_real64, 'wofost.phenology.vernsat',  errors)
         call check_real_range(self%vernbase, 0.0_real64, 100.0_real64, 'wofost.phenology.vernbase', errors)
         call check_real_range(self%verndvs,  0.0_real64,   0.3_real64, 'wofost.phenology.verndvs',  errors)
         call check_table(self%verntb, 2, .true., 'wofost.phenology.verntb', errors)
      end if
   end subroutine wofost_phenology_validate

   subroutine wofost_initial_validate(self, errors)
      class(wofost_initial_t),  intent(in)    :: self
      type(error_collection_t), intent(inout) :: errors
      call check_real_range(self%tdwi,   0.0_real64, 10000.0_real64, 'wofost.initial.tdwi',   errors)
      call check_real_range(self%laiem,  0.0_real64,    10.0_real64, 'wofost.initial.laiem',  errors)
      call check_real_range(self%rgrlai, 0.0_real64,     1.0_real64, 'wofost.initial.rgrlai', errors)
   end subroutine wofost_initial_validate

   subroutine wofost_greenarea_validate(self, errors)
      class(wofost_greenarea_t), intent(in)    :: self
      type(error_collection_t),  intent(inout) :: errors
      call check_real_range(self%spa,   0.0_real64,   1.0_real64, 'wofost.green_area.spa',   errors)
      call check_real_range(self%ssa,   0.0_real64,   1.0_real64, 'wofost.green_area.ssa',   errors)
      call check_real_range(self%span,  0.0_real64, 366.0_real64, 'wofost.green_area.span',  errors)
      call check_real_range(self%tbase,-10.0_real64,  30.0_real64, 'wofost.green_area.tbase', errors)
      call check_table(self%slatb, 2, .true., 'wofost.green_area.slatb', errors)
   end subroutine wofost_greenarea_validate

   subroutine wofost_assimilation_validate(self, errors)
      class(wofost_assimilation_t), intent(in)    :: self
      type(error_collection_t),     intent(inout) :: errors
      call check_real_range(self%kdif, 0.0_real64,  2.0_real64, 'wofost.assimilation.kdif', errors)
      call check_real_range(self%kdir, 0.0_real64,  2.0_real64, 'wofost.assimilation.kdir', errors)
      call check_real_range(self%eff,  0.0_real64, 10.0_real64, 'wofost.assimilation.eff',  errors)
      call check_table(self%amaxtb, 2, .true., 'wofost.assimilation.amaxtb', errors)
      call check_table(self%tmpftb, 2, .true., 'wofost.assimilation.tmpftb', errors)
      call check_table(self%tmnftb, 2, .true., 'wofost.assimilation.tmnftb', errors)
   end subroutine wofost_assimilation_validate

   subroutine wofost_conversion_validate(self, errors)
      class(wofost_conversion_t), intent(in)    :: self
      type(error_collection_t),   intent(inout) :: errors
      call check_real_range(self%cvl, 0.0_real64, 1.0_real64, 'wofost.conversion.cvl', errors)
      call check_real_range(self%cvo, 0.0_real64, 1.0_real64, 'wofost.conversion.cvo', errors)
      call check_real_range(self%cvr, 0.0_real64, 1.0_real64, 'wofost.conversion.cvr', errors)
      call check_real_range(self%cvs, 0.0_real64, 1.0_real64, 'wofost.conversion.cvs', errors)
   end subroutine wofost_conversion_validate

   subroutine wofost_respiration_validate(self, errors)
      class(wofost_respiration_t), intent(in)    :: self
      type(error_collection_t),    intent(inout) :: errors
      call check_real_range(self%q10, 0.0_real64, 5.0_real64, 'wofost.respiration.q10', errors)
      call check_real_range(self%rml, 0.0_real64, 1.0_real64, 'wofost.respiration.rml', errors)
      call check_real_range(self%rmo, 0.0_real64, 1.0_real64, 'wofost.respiration.rmo', errors)
      call check_real_range(self%rmr, 0.0_real64, 1.0_real64, 'wofost.respiration.rmr', errors)
      call check_real_range(self%rms, 0.0_real64, 1.0_real64, 'wofost.respiration.rms', errors)
      call check_table(self%rfsetb, 2, .true., 'wofost.respiration.rfsetb', errors)
   end subroutine wofost_respiration_validate

   subroutine wofost_partitioning_validate(self, errors)
      class(wofost_partitioning_t), intent(in)    :: self
      type(error_collection_t),     intent(inout) :: errors
      call check_table(self%frtb, 2, .true., 'wofost.partitioning.frtb', errors)
      call check_table(self%fltb, 2, .true., 'wofost.partitioning.fltb', errors)
      call check_table(self%fstb, 2, .true., 'wofost.partitioning.fstb', errors)
      call check_table(self%fotb, 2, .true., 'wofost.partitioning.fotb', errors)
   end subroutine wofost_partitioning_validate

   subroutine wofost_death_validate(self, errors)
      class(wofost_death_t),    intent(in)    :: self
      type(error_collection_t), intent(inout) :: errors
      call check_real_range(self%perdl, 0.0_real64, 3.0_real64, 'wofost.death.perdl', errors)
      call check_table(self%rdrrtb, 2, .true., 'wofost.death.rdrrtb', errors)
      call check_table(self%rdrstb, 2, .true., 'wofost.death.rdrstb', errors)
   end subroutine wofost_death_validate

   subroutine wofost_root_validate(self, errors)
      class(wofost_root_t),     intent(in)    :: self
      type(error_collection_t), intent(inout) :: errors
      call check_int_enum(self%swrd,     [1, 2, 3], 'wofost.root.swrd',     errors)
      call check_int_enum(self%swdmi2rd, [0, 1],    'wofost.root.swdmi2rd', errors)
      if (self%swrd == 1) then
         call check_table(self%rdtb, 2, .true., 'wofost.root.rdtb', errors)
      end if
      if (self%swrd == 2) then
         call check_real_range(self%rdi, 0.0_real64, 1000.0_real64, 'wofost.root.rdi', errors)
         call check_real_range(self%rri, 0.0_real64,  100.0_real64, 'wofost.root.rri', errors)
         call check_real_range(self%rdc, 0.0_real64, 1000.0_real64, 'wofost.root.rdc', errors)
      end if
      if (self%swrd == 3) then
         call check_table(self%rlwtb, 2, .true., 'wofost.root.rlwtb', errors)
         call check_real_range(self%wrtmax, 0.0_real64, 1.0e5_real64, 'wofost.root.wrtmax', errors)
      end if
      call check_table(self%rdctb, 2, .true., 'wofost.root.rdctb', errors)
   end subroutine wofost_root_validate

   subroutine wofost_oxygen_stress_validate(self, errors)
      class(wofost_oxygen_stress_t), intent(in)    :: self
      type(error_collection_t),      intent(inout) :: errors
      call check_int_enum(self%swoxygen,   [0, 1, 2], 'wofost.oxygen_stress.swoxygen',   errors)
      call check_int_enum(self%swwrtnonox, [0, 1],    'wofost.oxygen_stress.swwrtnonox', errors)
      call check_real_range(self%aeratecrit, 0.0001_real64, 1.0_real64, 'wofost.oxygen_stress.aeratecrit', errors)
      if (self%swoxygen == 1) then
         call check_real_range(self%hlim1,  -100.0_real64,  100.0_real64, 'wofost.oxygen_stress.hlim1',  errors)
         call check_real_range(self%hlim2u, -1000.0_real64, 100.0_real64, 'wofost.oxygen_stress.hlim2u', errors)
         call check_real_range(self%hlim2l, -1000.0_real64, 100.0_real64, 'wofost.oxygen_stress.hlim2l', errors)
      end if
      if (self%swoxygen == 2) then
         call check_real_range(self%q10_microbial,       1.0_real64,    4.0_real64, 'wofost.oxygen_stress.q10_microbial',       errors)
         call check_real_range(self%specific_resp_humus, 0.0_real64,    1.0_real64, 'wofost.oxygen_stress.specific_resp_humus', errors)
         call check_real_range(self%srl,                 0.0_real64, 1.0e10_real64, 'wofost.oxygen_stress.srl',                 errors)
         call check_int_enum(self%swrootradius, [1, 2], 'wofost.oxygen_stress.swrootradius', errors)
         if (self%swrootradius == 1) then
            call check_real_range(self%dry_mat_cont_roots,      0.0_real64,    1.0_real64, 'wofost.oxygen_stress.dry_mat_cont_roots',      errors)
            call check_real_range(self%air_filled_root_por,     0.0_real64,    1.0_real64, 'wofost.oxygen_stress.air_filled_root_por',     errors)
            call check_real_range(self%spec_weight_root_tissue, 0.0_real64, 1.0e5_real64,  'wofost.oxygen_stress.spec_weight_root_tissue', errors)
            call check_real_range(self%var_a,                   0.0_real64,    1.0_real64, 'wofost.oxygen_stress.var_a',                   errors)
         end if
         if (self%swrootradius == 2) then
            call check_real_range(self%root_radiusO2, 1.0e-6_real64, 0.1_real64, 'wofost.oxygen_stress.root_radiusO2', errors)
         end if
      end if
   end subroutine wofost_oxygen_stress_validate

   subroutine wofost_drought_stress_validate(self, errors)
      class(wofost_drought_stress_t), intent(in)    :: self
      type(error_collection_t),       intent(inout) :: errors
      call check_int_enum(self%swdrought, [1, 2], 'wofost.drought_stress.swdrought', errors)
      call check_real_range(self%hlim3h, -1.0e4_real64,    100.0_real64, 'wofost.drought_stress.hlim3h', errors)
      call check_real_range(self%hlim3l, -1.0e4_real64,    100.0_real64, 'wofost.drought_stress.hlim3l', errors)
      call check_real_range(self%hlim4,  -1.6e4_real64,    100.0_real64, 'wofost.drought_stress.hlim4',  errors)
      call check_real_range(self%adcrh,    0.0_real64,       5.0_real64, 'wofost.drought_stress.adcrh',  errors)
      call check_real_range(self%adcrl,    0.0_real64,       5.0_real64, 'wofost.drought_stress.adcrl',  errors)
   end subroutine wofost_drought_stress_validate

   subroutine wofost_salinity_validate(self, errors)
      class(wofost_salinity_t), intent(in)    :: self
      type(error_collection_t), intent(inout) :: errors
      call check_int_enum(self%swsalinity, [0, 1, 2], 'wofost.salinity.swsalinity', errors)
      if (self%swsalinity == 1) then
         call check_real_range(self%saltmax,   0.0_real64, 100.0_real64, 'wofost.salinity.saltmax',   errors)
         call check_real_range(self%saltslope, 0.0_real64,   1.0_real64, 'wofost.salinity.saltslope', errors)
      end if
      if (self%swsalinity == 2) then
         call check_real_range(self%salthead, 0.0_real64, 1000.0_real64, 'wofost.salinity.salthead', errors)
      end if
   end subroutine wofost_salinity_validate

   subroutine wofost_compensate_validate(self, errors)
      class(wofost_compensate_t), intent(in)    :: self
      type(error_collection_t),   intent(inout) :: errors
      call check_int_enum(self%swcompensate, [0, 1, 2], 'wofost.compensate.swcompensate', errors)
      if (self%swcompensate == 1 .or. self%swcompensate == 2) then
         call check_int_enum(self%swstressor, [1, 2, 3, 4, 5], 'wofost.compensate.swstressor', errors)
      end if
      if (self%swcompensate == 1) then
         call check_real_range(self%alphacrit, 0.2_real64, 1.0_real64, 'wofost.compensate.alphacrit', errors)
      end if
      if (self%swcompensate == 2) then
         call check_real_range(self%dcritrtz, 0.02_real64, 100.0_real64, 'wofost.compensate.dcritrtz', errors)
      end if
   end subroutine wofost_compensate_validate

   subroutine wofost_interception_validate(self, errors)
      class(wofost_interception_t), intent(in)    :: self
      type(error_collection_t),     intent(inout) :: errors
      call check_int_enum(self%swinter, [0, 1, 2], 'wofost.interception.swinter', errors)
      if (self%swinter == 1) then
         call check_real_range(self%cofab, 0.0_real64, 1.0_real64, 'wofost.interception.cofab', errors)
      end if
      if (self%swinter == 2) then
         call check_table(self%gashtb, 6, .true., 'wofost.interception.gashtb', errors)
      end if
   end subroutine wofost_interception_validate

   subroutine wofost_co2_validate(self, errors)
      class(wofost_co2_t),      intent(in)    :: self
      type(error_collection_t), intent(inout) :: errors
      call check_int_enum(self%swco2, [0, 1], 'wofost.co2.swco2', errors)
      if (self%swco2 == 1) then
         call check_table(self%co2amaxtb, 2, .true., 'wofost.co2.co2amaxtb', errors)
         call check_table(self%co2efftb,  2, .true., 'wofost.co2.co2efftb',  errors)
         call check_table(self%co2tratb,  2, .true., 'wofost.co2.co2tratb',  errors)
      end if
   end subroutine wofost_co2_validate

   subroutine wofost_management_validate(self, errors)
      class(wofost_management_t), intent(in)    :: self
      type(error_collection_t),   intent(inout) :: errors
      call check_real_range(self%fraharlosorm_lv,     0.0_real64, 1.0_real64, 'wofost.management.fraharlosorm_lv',     errors)
      call check_real_range(self%fraharlosorm_st,     0.0_real64, 1.0_real64, 'wofost.management.fraharlosorm_st',     errors)
      call check_real_range(self%fraharlosorm_so,     0.0_real64, 1.0_real64, 'wofost.management.fraharlosorm_so',     errors)
      call check_real_range(self%fradeceasedlvtosoil, 0.0_real64, 1.0_real64, 'wofost.management.fradeceasedlvtosoil', errors)
      call check_int_enum(self%swpotrelmf, [1, 2], 'wofost.management.swpotrelmf', errors)
      call check_real_range(self%relmf, 0.0_real64, 1.0_real64, 'wofost.management.relmf', errors)
   end subroutine wofost_management_validate

   subroutine cropwofost_config_validate(self, errors)
      class(cropwofost_config_t), intent(in)    :: self
      type(error_collection_t),   intent(inout) :: errors

      call self%preparation%validate(errors)
      call self%sowing%validate(errors)
      call self%germination%validate(errors)
      call self%harvest%validate(errors)
      call self%crop_factor%validate(errors)
      call self%phenology%validate(errors)
      call self%initial%validate(errors)
      call self%green_area%validate(errors)
      call self%assimilation%validate(errors)
      call self%conversion%validate(errors)
      call self%respiration%validate(errors)
      call self%partitioning%validate(errors)
      call self%death%validate(errors)
      call self%root%validate(errors)
      call self%oxygen_stress%validate(errors)
      call self%drought_stress%validate(errors)
      call self%salinity%validate(errors)
      call self%compensate%validate(errors)
      call self%interception%validate(errors)
      call self%co2%validate(errors)
      call self%management%validate(errors)
   end subroutine cropwofost_config_validate

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
