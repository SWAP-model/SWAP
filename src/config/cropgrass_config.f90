!> Type 3 (grass / WOFOST grass) crop config — populated from a .crp.toml file.
!! Field set is similar to cropfixed_config_t plus grass-specific management
!! (mowing, grazing, fertilizer). Per design Q2A, no shared base type.
module cropgrass_config_mod
   use iso_fortran_env, only: real64
   use error_mod, only: error_collection_t, ERR_VALIDATION_OUT_OF_RANGE, &
                        ERR_VALIDATION_CROSS_FIELD
   use validation_mod, only: check_int_enum, check_int_range, check_real_range, &
                             check_nonnegative_real, check_ordered_pair
   use irrigation_config_mod, only: irrigation_schedule_t
   implicit none
   private

   public :: cropgrass_config_t

   type :: cropgrass_config_t
      ! Phenology / development
      integer      :: idev = 2                  !! For grass typically 2 (temp-sum-based)
      integer      :: lcc  = 0
      real(real64) :: tbase = 0.0_real64
      real(real64) :: tsum1 = 0.0_real64
      real(real64) :: tsum2 = 0.0_real64

      ! Light & growth
      real(real64) :: kdif = 0.0_real64
      real(real64) :: kdir = 0.0_real64
      real(real64) :: eff  = 0.0_real64
      real(real64) :: amax = 0.0_real64

      ! Tables (pre-existing)
      real(real64), allocatable :: cftb(:)
      real(real64), allocatable :: chtb(:)
      real(real64), allocatable :: rdctb(:)

      ! Root growth
      real(real64) :: rdi = 0.0_real64
      real(real64) :: rri = 0.0_real64
      real(real64) :: rdc = 0.0_real64

      ! Water stress (Feddes) — same as fixed
      real(real64) :: hlim1 = 0.0_real64
      real(real64) :: hlim2u = 0.0_real64
      real(real64) :: hlim2l = 0.0_real64
      real(real64) :: hlim3h = 0.0_real64
      real(real64) :: hlim3l = 0.0_real64
      real(real64) :: hlim4 = 0.0_real64
      real(real64) :: adcrh = 0.0_real64
      real(real64) :: adcrl = 0.0_real64
      real(real64) :: rsc   = 0.0_real64

      ! Salinity & interception (pre-existing scalars)
      real(real64) :: ecmax  = 0.0_real64
      real(real64) :: ecslop = 0.0_real64
      real(real64) :: cofab  = 0.0_real64

      ! ====================================================================
      ! Phase 3 (.crp port) additions — new scalar and table fields
      ! ====================================================================

      ! Crop factor / height switch (Part 1)
      integer      :: swcf   = 2              !! 1=crop factor, 2=crop height, 3=wet (stub)
      real(real64) :: albedo = 0.23d0         !! Reflection coeff (when swcf=2)
      real(real64) :: rsw    = 0.0d0          !! Canopy resistance intercepted water

      ! Interception switch (Part 13)
      integer :: swinter = 1                  !! 0=none, 1=Von Hoyningen (supported), 2/3=stub

      ! Crop state initialisation (Part 2a)
      real(real64) :: tdwi   = 1000.0d0       !! Initial total crop dry weight [kg/ha]
      real(real64) :: laiem  = 0.63d0         !! Leaf area index at emergence
      real(real64) :: rgrlai = 0.007d0        !! Max relative increase in LAI

      ! Start-of-growth trigger (Part 2b)
      integer      :: swtsum    = 1           !! 0=none, 1=air-temp sum, 2=soil temp (2=stub)
      real(real64) :: tsumtemp  = 0.0d0       !! Soil temp threshold (when swtsum=2)
      real(real64) :: tsumdepth = 0.0d0       !! Soil depth for temp (when swtsum=2)
      integer      :: tsumtime  = 0           !! Consecutive days above threshold

      ! Green area (Part 3)
      real(real64) :: ssa  = 0.0d0            !! Specific stem area [ha/kg]
      real(real64) :: span = 30.0d0           !! Leaf life span [d]
      ! (tbase already declared above)
      real(real64), allocatable :: slatb(:)   !! Specific leaf area vs DNR (flat pairs)

      ! Assimilation tables (Part 4)
      real(real64), allocatable :: amaxtb(:)  !! Max CO2 assim rate vs DNR
      real(real64), allocatable :: tmpftb(:)  !! AMAX reduction vs avg temp
      real(real64), allocatable :: tmnftb(:)  !! AMAX reduction vs min temp

      ! Biomass conversion (Part 5)
      real(real64) :: cvl = 0.685d0           !! Leaves
      real(real64) :: cvr = 0.694d0           !! Roots
      real(real64) :: cvs = 0.662d0           !! Stems

      ! Maintenance respiration (Part 6)
      real(real64) :: q10 = 2.0d0             !! Q10 factor
      real(real64) :: rml = 0.03d0            !! Leaves
      real(real64) :: rmr = 0.015d0           !! Roots
      real(real64) :: rms = 0.015d0           !! Stems
      real(real64), allocatable :: rfsetb(:)  !! Senescence reduction vs DNR

      ! Partitioning (Part 7)
      real(real64), allocatable :: frtb(:)    !! Roots fraction vs DNR
      real(real64), allocatable :: fltb(:)    !! Leaves fraction vs DNR
      real(real64), allocatable :: fstb(:)    !! Stems fraction vs DNR

      ! Death rates (Part 8)
      real(real64) :: perdl = 0.05d0          !! Max death rate leaves under water stress
      real(real64), allocatable :: rdrrtb(:)  !! Root death rate vs DNR
      real(real64), allocatable :: rdrstb(:)  !! Stem death rate vs DNR

      ! Root depth and density switches (Part 9)
      integer      :: swrd     = 2            !! 1=DVS table (stub), 2=daily increase, 3=biomass (stub)
      integer      :: swdmi2rd = 0            !! 0=assimilates, 1=DM increase
      integer      :: swrdc    = 0            !! 0=unmodified, 1=modified by water extract (stub)
      real(real64) :: wrtmax   = 3000.0d0     !! Max root weight (when swrd=3)
      real(real64), allocatable :: rdtb(:)    !! Rooting depth vs DNR (when swrd=1)
      real(real64), allocatable :: rlwtb(:)   !! Rooting depth vs root weight (swrd=3)

      ! Oxygen stress additions (Part 10)
      integer      :: swoxygen   = 1          !! 0=none, 1=Feddes, 2=Bartholomeus
      integer      :: swwrtnonox = 0          !! Check aerobic conditions for root growth
      real(real64) :: aeratecrit = 1.0d-4     !! Aerobic threshold for root extension

      ! Bartholomeus physical sub-model (swoxygentype=1 supported; =2 stub)
      integer      :: swoxygentype          = 1
      real(real64) :: q10_microbial         = 2.8d0
      real(real64) :: specific_resp_humus   = 1.6d-3
      real(real64) :: srl                   = 383571.0d0
      integer      :: swrootradius          = 2  !! 1=calculated, 2=given
      real(real64) :: dry_mat_cont_roots    = 0.075d0
      real(real64) :: air_filled_root_por   = 0.05d0
      real(real64) :: spec_weight_root_tissue = 1.0d3
      real(real64) :: var_a                 = 4.175d-10
      real(real64) :: root_radiusO2         = 0.000075d0

      ! Drought stress switch (Part 11)
      integer :: swdrought = 1               !! 1=Feddes (supported), 2=De Jong (stub)

      ! Salinity stress switch (Part 12)
      integer :: swsalinity = 0             !! 0=none (others stub)

      ! Root water uptake compensation
      integer      :: swcompensate = 0      !! 0=none, 1=Jarvis (supported), 2=Walsum (stub)
      integer      :: swstressor   = 1      !! Stressors to compensate
      real(real64) :: alphacrit    = 1.0d0  !! Jarvis critical index
      real(real64) :: dcritrtz     = 0.0d0  !! Walsum threshold (stub)

      ! CO2 switch (Part 14)
      integer :: swco2 = 0                  !! 0=none (swco2=1 stub)

      ! ====================================================================
      ! End of Phase 3 scalar additions
      ! ====================================================================

      ! Mowing schedule (grass-specific) — pre-existing
      ! swharv: 0 = no scheduled mow, 1 = DM-threshold-driven, 2 = fixed-date table
      integer :: swharv = 0
      integer :: nmow   = 0                     !! Number of mowing events
      real(real64), allocatable :: dates_mowing(:)
      real(real64), allocatable :: lai_after_mow(:)

      ! Phase 4d additions for per-event mowing/grazing tables.
      ! Mowing block (when swharv=1 or 2):
      integer      :: swdmmow = 0                           !! 0=use heights, 1=DM threshold, 2=DM threshold (legacy SWDMMOW=2 in cases)
      real(real64), allocatable :: mowing_dates(:)          !! day-of-year per event (when swharv=2)
      real(real64), allocatable :: mowing_heights(:)        !! optional (when swdmmow=0)
      real(real64) :: dmharvest      = 0.0d0                !! DM threshold (when swdmmow=1)
      real(real64) :: daylastharvest = 0.0d0
      real(real64) :: dmlastharvest  = 0.0d0
      integer      :: maxdaymow = 0

      ! Grazing (grass-specific) — pre-existing
      integer :: swgraz = 0                     !! 0=no grazing, 1=scheduled
      integer :: nstart_graz = 0                !! Start day-of-year
      integer :: nstop_graz  = 0                !! Stop day-of-year

      ! Phase 4d additions for grazing per-event details (when swgraz=1):
      integer      :: maxdaygrz = 0
      real(real64) :: dmgrazing = 0.0d0
      integer      :: swdmgrz   = 0
      real(real64), allocatable :: lsdb(:)                  !! per-day stocking density
      real(real64) :: tagprest = 0.0d0                      !! threshold above-ground residue

      ! Management general (MANAGEMENT SECTION)
      integer :: nseqgrazmow = 20              !! Number of periods in SEQGRAZMOW
      integer, allocatable :: seqgrazmow(:)    !! 1=graze, 2=mow, 3=dewool (only 2 supported)
      real(real64) :: mowrest    = 700.0d0     !! Remaining DM after mowing [kg/ha]
      real(real64) :: dewrest    = 850.0d0     !! Remaining DM after dewooling (stub)
      integer      :: swpotrelmf = 1           !! 1=theoretical, 2=attainable yield
      real(real64) :: relmf      = 1.0d0       !! Relative management factor

      ! Mowing DM-threshold flexible table (when swdmmow=2)
      real(real64), allocatable :: dmmowtb(:)   !! DM threshold vs DNR (flat pairs)

      ! Mowing regrowth delay table
      real(real64), allocatable :: dmmowdelay(:) !! DM harvest vs delay days (flat pairs)

      ! Grazing flexible DM threshold table (when swdmgrz=2)
      real(real64), allocatable :: dmgrztb(:)   !! DM threshold vs DNR (flat pairs)

      ! Grazing livestock density tables (schema 1:1)
      real(real64), allocatable :: lsda(:)          !! Actual livestock density per period
      real(real64), allocatable :: daysgrazing(:)   !! Max days grazing per LSDb entry
      real(real64), allocatable :: uptgrazing(:)    !! DM uptake per LSDb entry
      real(real64), allocatable :: lossgrazing(:)   !! DM loss per LSDb entry

      ! Treading loss switches (stub-errored at validator)
      integer :: swlossmow = 0               !! 0=no losses (others stub)
      integer :: swlossgrz = 0              !! 0=no losses (others stub)

      ! Per-crop irrigation schedule (Phase 4d Task 12)
      type(irrigation_schedule_t) :: schedule
   contains
      procedure :: validate => cropgrass_config_validate
      procedure :: finalize => cropgrass_config_finalize
   end type cropgrass_config_t

contains

   subroutine cropgrass_config_validate(self, errors)
      class(cropgrass_config_t), intent(in)    :: self
      type(error_collection_t),  intent(inout) :: errors

      ! ----- Phase 3 stub-errors for unsupported runtime branches -----
      ! ADR 0015: schema accepts these values 1:1 with legacy; runtime
      ! plumbing for these branches has not been ported. Cases that
      ! need them must run via the legacy executable.
      if (self%swdrought == 2) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropgrass.swdrought=2 (De Jong van Lier) not yet supported in ' // &
            'the TOML pipeline; use the legacy executable.', 'cropgrass')
      end if
      if (self%swcompensate == 2) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropgrass.swcompensate=2 (Walsum) not yet supported in ' // &
            'the TOML pipeline; use the legacy executable.', 'cropgrass')
      end if
      if (self%swsalinity /= 0) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropgrass.swsalinity /= 0 not yet supported in the TOML ' // &
            'pipeline; use the legacy executable.', 'cropgrass')
      end if
      if (self%swinter == 2 .or. self%swinter == 3) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropgrass.swinter=2 or 3 (Gash/storage-cap) not yet supported ' // &
            'in the TOML pipeline; use the legacy executable.', 'cropgrass')
      end if
      if (self%swco2 == 1) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropgrass.swco2=1 (CO2 correction) not yet supported in the ' // &
            'TOML pipeline; use the legacy executable.', 'cropgrass')
      end if
      if (self%swlossmow == 1) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropgrass.swlossmow=1 (treading losses during mowing) not yet ' // &
            'supported in the TOML pipeline; use the legacy executable.', &
            'cropgrass')
      end if
      if (self%swlossgrz == 1) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropgrass.swlossgrz=1 (treading losses during grazing) not yet ' // &
            'supported in the TOML pipeline; use the legacy executable.', &
            'cropgrass')
      end if
      if (self%swoxygen == 2 .and. self%swoxygentype == 2) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropgrass.swoxygentype=2 (reproduction functions) not yet ' // &
            'supported in the TOML pipeline; use the legacy executable.', &
            'cropgrass')
      end if
      if (self%swrd == 1) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropgrass.swrd=1 (DVS-table root depth) not yet supported in ' // &
            'the TOML pipeline; use the legacy executable.', 'cropgrass')
      end if
      ! swrd=3 (biomass-based root extension via rlwtb/wrtmax) is now supported;
      ! cropgrass_init copies rlwtb and wrtmax to legacy globals.
      if (self%swcf == 3) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropgrass.swcf=3 (LAI-dependent dual-coeff) not yet supported ' // &
            'in the TOML pipeline; use the legacy executable.', 'cropgrass')
      end if
      if (self%swrdc == 1) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropgrass.swrdc=1 (root density modified by water extraction) ' // &
            'not yet supported in the TOML pipeline; use the legacy ' // &
            'executable.', 'cropgrass')
      end if
      if (self%swtsum == 2) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropgrass.swtsum=2 (soil-temperature start trigger) not yet ' // &
            'supported in the TOML pipeline; use the legacy executable.', &
            'cropgrass')
      end if
      if (self%schedule%schedule == 1) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropgrass.schedule=1 (per-crop irrigation scheduling) not yet ' // &
            'supported in the TOML pipeline; use the legacy executable.', &
            'cropgrass')
      end if
      if (allocated(self%seqgrazmow)) then
         block
            integer :: jseq
            do jseq = 1, size(self%seqgrazmow)
               if (self%seqgrazmow(jseq) /= 2) then
                  call errors%append(ERR_VALIDATION_CROSS_FIELD, &
                     'cropgrass.seqgrazmow: grazing (=1) and dewooling (=3) ' // &
                     'periods are not yet supported in the TOML pipeline; ' // &
                     'use the legacy executable.', 'cropgrass')
                  exit
               end if
            end do
         end block
      end if

      ! ----- Pre-existing enum / range checks -----
      call check_int_enum(self%idev,   [1, 2],       'cropgrass.idev',   errors)
      call check_int_enum(self%swharv, [0, 1, 2],    'cropgrass.swharv', errors)
      call check_int_enum(self%swgraz, [0, 1],       'cropgrass.swgraz', errors)

      call check_real_range(self%rdi, 0.0_real64, 1000.0_real64, 'cropgrass.rdi', errors)
      call check_real_range(self%rdc, 0.0_real64, 1000.0_real64, 'cropgrass.rdc', errors)
      call check_nonnegative_real(self%rri, 'cropgrass.rri', errors)

      call check_ordered_pair(self%hlim2l, self%hlim2u, 'hlim2l', 'hlim2u', 'cropgrass', errors)
      call check_ordered_pair(self%hlim3l, self%hlim3h, 'hlim3l', 'hlim3h', 'cropgrass', errors)

      call check_nonnegative_real(self%ecmax,  'cropgrass.ecmax',  errors)
      call check_nonnegative_real(self%ecslop, 'cropgrass.ecslop', errors)

      ! ----- Mowing branches (swharv = 0/1/2) -----
      if (self%swharv == 1) then
         ! DM-threshold mowing: nmow drives event budget; mowing_dates not required
         call check_int_range(self%nmow, 1, 50, 'cropgrass.nmow', errors)
         call check_int_enum(self%swdmmow, [0, 1, 2], 'cropgrass.swdmmow', errors)
         if (self%swdmmow == 1) then
            if (self%dmharvest <= 0.0_real64) then
               call check_real_range(self%dmharvest, tiny(1.0_real64), huge(1.0_real64), &
                                     'cropgrass.dmharvest', errors)
            end if
         end if
         call check_nonnegative_real(self%daylastharvest, 'cropgrass.daylastharvest', errors)
         call check_nonnegative_real(self%dmlastharvest,  'cropgrass.dmlastharvest',  errors)
         call check_int_range(self%maxdaymow, 1, 366, 'cropgrass.maxdaymow', errors)
      else if (self%swharv == 2) then
         ! Fixed-date mowing: require mowing_dates allocated and nmow == size
         call check_int_range(self%nmow, 1, 50, 'cropgrass.nmow', errors)
         call check_int_enum(self%swdmmow, [0, 1, 2], 'cropgrass.swdmmow', errors)
         if (.not. allocated(self%mowing_dates)) then
            call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
                               'required when swharv=2 (allocate event table)', &
                               'cropgrass.mowing_dates')
         else if (size(self%mowing_dates) /= self%nmow) then
            call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
                               'size must equal nmow', &
                               'cropgrass.mowing_dates')
         end if
         if (self%swdmmow == 0) then
            if (.not. allocated(self%mowing_heights)) then
               call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
                                  'required when swharv=2 and swdmmow=0', &
                                  'cropgrass.mowing_heights')
            else if (size(self%mowing_heights) /= self%nmow) then
               call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
                                  'size must equal nmow', &
                                  'cropgrass.mowing_heights')
            end if
         end if
      end if

      ! ----- Grazing branch (swgraz = 1) -----
      if (self%swgraz == 1) then
         call check_int_range(self%nstart_graz, 1, 366, 'cropgrass.nstart_graz', errors)
         call check_int_range(self%nstop_graz,  1, 366, 'cropgrass.nstop_graz',  errors)
         if (self%nstop_graz < self%nstart_graz) then
            call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
                               'must be >= nstart_graz', &
                               'cropgrass.nstop_graz')
         end if
         call check_int_range(self%maxdaygrz, 1, 366, 'cropgrass.maxdaygrz', errors)
         if (self%dmgrazing <= 0.0_real64) then
            call check_real_range(self%dmgrazing, tiny(1.0_real64), huge(1.0_real64), &
                                  'cropgrass.dmgrazing', errors)
         end if
         call check_int_enum(self%swdmgrz, [0, 1, 2], 'cropgrass.swdmgrz', errors)
         call check_nonnegative_real(self%tagprest, 'cropgrass.tagprest', errors)
      end if

      ! Optional lsdb: if allocated, sanity-check non-negative entries
      if (allocated(self%lsdb)) then
         block
            integer :: i
            do i = 1, size(self%lsdb)
               if (self%lsdb(i) < 0.0_real64) then
                  call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
                                     'stocking density must be non-negative', &
                                     'cropgrass.lsdb')
                  exit
               end if
            end do
         end block
      end if

      call self%schedule%validate(errors)
   end subroutine cropgrass_config_validate

   subroutine cropgrass_config_finalize(self, errors)
      class(cropgrass_config_t), intent(inout) :: self
      type(error_collection_t),  intent(inout) :: errors
      call self%schedule%finalize(errors)
   end subroutine cropgrass_config_finalize

end module cropgrass_config_mod
