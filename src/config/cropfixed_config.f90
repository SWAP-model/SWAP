!> Type 1 (fixed/simple) crop config — populated from a .crp.toml file.
!! Phase 1 of the .crp port (Phase 4f) extended this to a 1:1 match with
!! legacy `readcropfixed` in src/io/readswap.f90: every legacy field is
!! present in the schema, regardless of whether the parent switch is
!! supported in Phase 1's runtime. Stub-errored switches are documented
!! per ADR 0015.
module cropfixed_config_mod
   use iso_fortran_env, only: real64
   use error_mod, only: error_collection_t, ERR_VALIDATION_OUT_OF_RANGE, &
                        ERR_VALIDATION_CROSS_FIELD
   use validation_mod, only: check_int_enum, check_int_range, check_real_range, &
                             check_nonnegative_real, check_ordered_pair
   use irrigation_config_mod, only: irrigation_schedule_t
   implicit none
   private

   public :: cropfixed_config_t

   type :: cropfixed_config_t
      ! ====================================================================
      ! Existing fields (Phase 4c-a) — KEPT
      ! ====================================================================
      integer      :: idev = 1
      integer      :: lcc  = 0
      real(real64) :: kdif = 0.0_real64
      real(real64) :: kdir = 0.0_real64
      real(real64), allocatable :: cftb(:)
      real(real64), allocatable :: chtb(:)
      real(real64) :: rdi = 0.0_real64
      real(real64) :: rri = 0.0_real64
      real(real64) :: rdc = 0.0_real64
      real(real64), allocatable :: rdctb(:)
      real(real64) :: hlim1 = 0.0_real64
      real(real64) :: hlim2u = 0.0_real64
      real(real64) :: hlim2l = 0.0_real64
      real(real64) :: hlim3h = 0.0_real64
      real(real64) :: hlim3l = 0.0_real64
      real(real64) :: hlim4 = 0.0_real64
      real(real64) :: adcrh = 0.0_real64
      real(real64) :: adcrl = 0.0_real64
      real(real64) :: rsc   = 0.0_real64
      real(real64) :: ecmax  = 0.0_real64
      real(real64) :: ecslop = 0.0_real64
      real(real64) :: cofab = 0.0_real64
      type(irrigation_schedule_t) :: schedule

      ! ====================================================================
      ! Phase 1 (Phase 4f .crp port) additions
      ! ====================================================================

      ! Part 0a/b/c — preparation, sowing, germination
      integer :: swprep = 0    !! 0=no preparation, 1=preparation
      integer :: swsow  = 0    !! 0=no sowing, 1=sowing before crop growth
      integer :: swgerm = 0    !! 0=no germination, 1/2=temperature-based

      ! Part 0d — harvest
      real(real64) :: dvsend = 2.0_real64   !! Development stage at harvest
      integer      :: swharv = 0            !! 0=CROPEND-based, 1=DVS-based (stub-errored)

      ! Part 1 (idev=2 path; safe defaults for idev=1)
      real(real64) :: tsumea = 0.0_real64
      real(real64) :: tsumam = 0.0_real64
      real(real64) :: tbase  = 0.0_real64

      ! Part 3 — LAI vs SCF
      integer                   :: swgc = 1   !! 1=LAI, 2=SCF
      real(real64), allocatable :: gctb(:)    !! (dvs, lai|scf) flat pairs

      ! Part 4 — crop factor / height switch
      integer                   :: swcf = 1   !! 1=crop factor (cftb), 2=crop height (chtb), 3=wet-crop (stub-errored)
      real(real64), allocatable :: cfeictb(:) !! Wet-crop factor table (only if swcf=3 - unused in Phase 1)
      real(real64) :: albedo = 0.23_real64
      real(real64) :: rsw    = 0.0_real64

      ! Part 10 — root depth & density
      integer :: swrd     = 1   !! 1=DVS table (rdtb), 2=daily increase (stub-errored), 3=biomass (fatal-errored in legacy)
      integer :: swdmi2rd = 0   !! Only used when swrd=2 (stub-errored)
      integer :: swrdc    = 0   !! Switch development root density (legacy hard-codes to 0)
      real(real64), allocatable :: rdtb(:)    !! (dvs, rd) flat pairs (when swrd=1)

      ! Part 11 — oxygen stress
      integer :: swoxygen   = 0   !! 0=none, 1=Feddes, 2=Bartholomeus (stub-errored).
                                  !! NOTE: TOML default 0; legacy readcropfixed default was 1.
      integer :: swwrtnonox = 0   !! 0=no aerobic check, 1=check
      real(real64) :: aeratecrit = 1.0e-4_real64    !! Required when swwrtnonox=1
      real(real64) :: max_resp_factor = 1.0_real64  !! Oxygen stress max respiration factor (ratio total/maintenance)

      ! Part 12 — drought stress
      integer :: swdrought = 1   !! 1=Feddes, 2=De Jong van Lier (stub-errored)

      ! Part 13 — salinity stress
      integer :: swsalinity = 0   !! 0=none, 1=Maas-Hoffman (stub-errored), 2=osmotic head (stub-errored)
      real(real64) :: saltmax   = 0.0_real64
      real(real64) :: saltslope = 0.0_real64
      real(real64) :: salthead  = 0.0_real64

      ! Part xx — root water uptake compensation
      integer :: swcompensate = 0   !! 0=none, 1=Jarvis (stub-errored), 2=Walsum (stub-errored)
      integer :: swstressor   = 1
      real(real64) :: alphacrit = 1.0_real64
      real(real64) :: dcritrtz  = 0.0_real64

      ! Part 14 — interception
      integer :: swinter = 1   !! 0=none, 1=Von Hoyningen-Hune (supported), 2=Gash (stub-errored), 3=storage-cap (stub-errored)

      ! Part 15 — irrigation scheduling top-level switch
      integer :: schedule_switch = 0   !! 0=no scheduling, 1=apply (stub-errored)

      ! Subordinate fields under stub-errored switches (schema 1:1; runtime never reads them)
      ! swdrought=2 fields (stub-errored at validate time):
      real(real64) :: wiltpoint  = 0.0_real64
      real(real64) :: kstem      = 0.0_real64
      real(real64) :: rxylem     = 0.0_real64
      real(real64) :: rootradius = 0.0_real64
      real(real64) :: kroot      = 0.0_real64
      real(real64) :: rootcoefa  = 0.0_real64
      real(real64) :: rooteff    = 0.0_real64
      real(real64) :: stephr     = 0.0_real64
      real(real64) :: criterhr   = 0.0_real64
      real(real64) :: taccur     = 0.0_real64
      ! swoxygen=2 fields (stub-errored at validate time):
      integer      :: swoxygentype          = 1
      integer      :: swrootradius          = 1
      integer      :: swtopsub              = 1
      integer      :: nrstaring             = 1
      real(real64) :: q10_root              = 0.0_real64
      real(real64) :: q10_microbial         = 0.0_real64
      real(real64) :: specific_resp_humus   = 0.0_real64
      real(real64) :: c_mroot               = 0.0_real64
      real(real64) :: srl                   = 0.0_real64
      real(real64) :: f_senes               = 0.0_real64
      real(real64) :: dry_mat_cont_roots    = 0.0_real64
      real(real64) :: air_filled_root_por   = 0.0_real64
      real(real64) :: spec_weight_root_tissue = 0.0_real64
      real(real64) :: var_a                 = 0.0_real64
      real(real64) :: root_radiusO2         = 1.0e-3_real64
   contains
      procedure :: validate => cropfixed_config_validate
      procedure :: finalize => cropfixed_config_finalize
   end type cropfixed_config_t

contains

   subroutine cropfixed_config_validate(self, errors)
      class(cropfixed_config_t), intent(in)    :: self
      type(error_collection_t),  intent(inout) :: errors

      ! ----- Phase 1 stub-errors for unsupported runtime branches -----
      ! ADR 0015: schema accepts the value 1:1 with legacy; runtime
      ! plumbing for these branches has not been ported. Cases that
      ! need them must run via the legacy executable.
      if (self%swdrought == 2) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropfixed.swdrought=2 (De Jong van Lier) not yet supported in ' // &
            'the TOML pipeline; use the legacy executable.', 'cropfixed')
      end if
      if (self%swoxygen == 2) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropfixed.swoxygen=2 (Bartholomeus) not yet supported in ' // &
            'the TOML pipeline; use the legacy executable.', 'cropfixed')
      end if
      if (self%swcf == 3) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropfixed.swcf=3 (wet-crop factor) not yet supported in the ' // &
            'TOML pipeline; use the legacy executable.', 'cropfixed')
      end if
      if (self%swinter == 2 .or. self%swinter == 3) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropfixed.swinter=2 or 3 (Gash / storage-cap interception) ' // &
            'not yet supported in the TOML pipeline; use the legacy ' // &
            'executable.', 'cropfixed')
      end if
      ! swrd=1 (DVS table) and swrd=2 (max daily increase) are both supported.
      ! swrd=3 (root extension from available root biomass) is not possible
      ! with the simple crop module — legacy readcropfixed fatal-errors it.
      if (self%swrd == 3) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropfixed.swrd=3 (root extension based on available root ' // &
            'biomass) is not possible with the simple crop module.', &
            'cropfixed')
      end if
      ! swsalinity=1 (Maas-Hoffman) and swsalinity=2 (osmotic head) both stay
      ! rejected: the Maas-Hoffman kernel is intact, but enabling it for a
      ! simple crop fails byte-identical regression against swap420gf — the
      ! salinity→uptake→solute-concentration feedback loop diverges (the only
      ! stress whose magnitude reads sol%cml). See dev-docs investigation note.
      if (self%swsalinity == 1 .or. self%swsalinity == 2) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropfixed.swsalinity != 0 (Maas-Hoffman / osmotic head) not ' // &
            'yet supported in the TOML pipeline; use the legacy executable.', &
            'cropfixed')
      end if
      if (self%schedule_switch == 1) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropfixed.schedule=1 (per-crop irrigation scheduling) not yet ' // &
            'supported in the TOML pipeline; use the legacy executable.', &
            'cropfixed')
      end if

      ! ----- Top-level enum / range checks -----
      call check_int_enum(self%idev, [1, 2], 'cropfixed.idev', errors)
      if (self%idev == 1) then
         call check_int_range(self%lcc, 1, 366, 'cropfixed.lcc', errors)
      end if

      call check_int_enum(self%swgc,    [1, 2], 'cropfixed.swgc',    errors)
      call check_int_enum(self%swcf,    [1, 2, 3], 'cropfixed.swcf', errors)
      call check_int_enum(self%swrd,    [1, 2, 3], 'cropfixed.swrd', errors)
      call check_int_enum(self%swoxygen,    [0, 1, 2], 'cropfixed.swoxygen',    errors)
      call check_int_enum(self%swwrtnonox,  [0, 1],    'cropfixed.swwrtnonox',  errors)
      call check_int_enum(self%swdrought,   [1, 2],    'cropfixed.swdrought',   errors)
      call check_int_enum(self%swsalinity,  [0, 1, 2], 'cropfixed.swsalinity',  errors)
      call check_int_enum(self%swcompensate,[0, 1, 2], 'cropfixed.swcompensate',errors)
      call check_int_enum(self%swinter,     [0, 1, 2, 3], 'cropfixed.swinter',  errors)
      call check_int_enum(self%swharv,      [0, 1],    'cropfixed.swharv',      errors)
      call check_int_enum(self%swprep,      [0, 1],    'cropfixed.swprep',      errors)
      call check_int_enum(self%swsow,       [0, 1],    'cropfixed.swsow',       errors)
      call check_int_enum(self%swgerm,      [0, 1, 2], 'cropfixed.swgerm',      errors)
      call check_int_enum(self%schedule_switch, [0, 1], 'cropfixed.schedule', errors)

      call check_real_range(self%dvsend, 0.0_real64, 3.0_real64,    'cropfixed.dvsend', errors)
      call check_real_range(self%kdif, 0.0_real64,  2.0_real64,     'cropfixed.kdif',   errors)
      call check_real_range(self%kdir, 0.0_real64,  2.0_real64,     'cropfixed.kdir',   errors)
      call check_real_range(self%rdi,  0.0_real64, 1000.0_real64,   'cropfixed.rdi',    errors)
      call check_real_range(self%rdc,  0.0_real64, 1000.0_real64,   'cropfixed.rdc',    errors)
      call check_nonnegative_real(self%rri, 'cropfixed.rri', errors)

      ! Feddes ordering — applied only when the corresponding branch is active.
      ! hlim1 (saturation, near 0) >= hlim2u >= hlim2l >= hlim3h >= hlim3l >= hlim4 (wilting).
      if (self%swoxygen == 1) then
         call check_ordered_pair(self%hlim2l, self%hlim2u, 'hlim2l', 'hlim2u', 'cropfixed', errors)
      end if
      if (self%swdrought == 1) then
         call check_ordered_pair(self%hlim3l, self%hlim3h, 'hlim3l', 'hlim3h', 'cropfixed', errors)
      end if

      call check_nonnegative_real(self%ecmax,  'cropfixed.ecmax',  errors)
      call check_nonnegative_real(self%ecslop, 'cropfixed.ecslop', errors)

      call self%schedule%validate(errors)
   end subroutine cropfixed_config_validate

   subroutine cropfixed_config_finalize(self, errors)
      class(cropfixed_config_t), intent(inout) :: self
      type(error_collection_t),  intent(inout) :: errors
      call self%schedule%finalize(errors)
   end subroutine cropfixed_config_finalize

end module cropfixed_config_mod
