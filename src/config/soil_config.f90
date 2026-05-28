!> [soil] section config.
module soil_config_mod
   use iso_fortran_env, only: real64
   use error_mod, only: error_collection_t,          &
                        ERR_VALIDATION_CROSS_FIELD,  &
                        ERR_VALIDATION_OUT_OF_RANGE
   use validation_mod, only: check_int_enum, check_int_range, &
                             check_real_range, check_nonnegative_real
   implicit none
   private

   public :: soil_config_t
   public :: soil_discretization_t
   public :: soil_frost_t
   public :: soil_hydraulics_t
   public :: soil_initial_t
   public :: soil_tillage_event_t
   public :: soil_tillage_type_t
   public :: soil_tillage_t

   !> Optional re-discretization of the vertical grid for output reporting.
   !! When swdiscrvert == 1, dznew(:) (sized to numnodnew) carries the
   !! per-new-node thickness in cm.
   type :: soil_discretization_t
      integer :: swdiscrvert = 0       !! 0=use existing, 1=re-discretize
      integer :: numnodnew   = 0       !! count when swdiscrvert=1
      real(real64), allocatable :: dznew(:)  !! per-new-node thickness (cm)
   contains
      procedure :: validate => soil_discretization_validate
   end type soil_discretization_t

   !> Per-soil-physical-layer Mualem-van Genuchten hydraulic parameters.
   !! All arrays are sized to `numlay` (number of soil-physical layers,
   !! the deepest `isoillay` value in the discretization). Required when
   !! `swsophy = 0` (analytical MvG); the validator does not enforce
   !! presence here because the runtime fatalerrs if values are missing
   !! at the time soilhydraulics builds `paramvg`. Field names mirror the
   !! legacy `.swp` keys verbatim (Phase 4f Task B2).
   type :: soil_hydraulics_t
      real(real64), allocatable :: ores(:)     !! residual water content [0..1]
      real(real64), allocatable :: osat(:)     !! saturated water content [0..1]
      real(real64), allocatable :: alfa(:)     !! MvG alpha drying [1e-4..100 /cm]
      real(real64), allocatable :: npar(:)     !! MvG n shape [1.001..9 -]
      real(real64), allocatable :: ksatfit(:)  !! fitted Ksat [1e-5..1e5 cm/d]
      real(real64), allocatable :: lexp(:)     !! K(h) exponent [-25..25 -]
      real(real64), allocatable :: alfaw(:)    !! MvG alpha wetting (hysteresis)
      real(real64), allocatable :: h_enpr(:)   !! air-entry pressure head [-40..0 cm]
      real(real64), allocatable :: ksatexm(:)  !! measured Ksat [1e-5..1e5 cm/d]
      real(real64), allocatable :: bdens(:)    !! dry bulk density [100..1e4 mg/cm3]
   contains
      procedure :: validate => soil_hydraulics_validate
   end type soil_hydraulics_t

   !> Frost-induced flow reduction parameters.
   type :: soil_frost_t
      integer      :: swfrost   = 0
      real(real64) :: tfroststa = 0.0_real64
      real(real64) :: tfrostend = 0.0_real64
      integer      :: swsublim  = 0
   contains
      procedure :: validate => soil_frost_validate
   end type soil_frost_t

   !> Initial-state inputs consumed when soil.swinco == 3.
   !! Replaces the legacy ASCII swap.ini file. Scalars copy directly to
   !! globals (ssnow, slw, pond, pondini, ldwet, dt, atmin7); the three
   !! z-indexed profiles are read via read_csv_table as separate companion
   !! CSVs. swirrigate is metadata only — no TOML-side global consumer.
   type :: soil_initial_t
      integer       :: swirrigate = 0
      real(real64)  :: ssnow  = 0.0_real64
      real(real64)  :: slw    = 0.0_real64
      real(real64)  :: pond   = 0.0_real64
      real(real64)  :: ldwet  = 0.0_real64
      real(real64)  :: dt     = 0.0_real64
      real(real64)  :: atmin7(7) = 0.0_real64
      character(len=:), allocatable :: h_file      !! header z,h
      character(len=:), allocatable :: tsoil_file  !! header z,tsoil
      character(len=:), allocatable :: cml_file    !! header z,cml
      ! [GR-SOIL 2026-05-24] swinco=1/3 initial-head table z-axis (depths, cm).
      ! TODO(W3 2026-05-28): DEPRECATED — no longer populated or consumed.
      ! soilwater_state_init now loads h_file into state%soilwater%h_init (h_profile_table_t)
      ! and soilhydraulics.f90 reads rows(i)%z from there directly. This field is
      ! allocatable and will simply be empty; retire it in a future cleanup arc.
      real(real64), allocatable :: z_init(:)
   contains
      procedure :: validate => soil_initial_validate
   end type soil_initial_t

   !> Single tillage event row read from `[[soil.tillage.events]]`.
   type :: soil_tillage_event_t
      character(len=10) :: date      = ''       !! ISO YYYY-MM-DD (parsed to days-since-1900 by adapter)
      real(real64)      :: z         = 0.0_real64
      real(real64)      :: intensity = 0.0_real64
      integer           :: type_id   = 0
   end type soil_tillage_event_t

   !> Single tillage type row read from `[[soil.tillage.types]]`.
   type :: soil_tillage_type_t
      integer      :: id          = 0
      real(real64) :: rho_cons    = 0.0_real64
      real(real64) :: rho_tillage = 0.0_real64
      real(real64) :: k_R         = 0.0_real64
      real(real64) :: rho_match   = -99.0_real64   !! used only when i_n_model = 3
      real(real64) :: N_match     = -99.0_real64   !! used only when i_n_model = 3
   end type soil_tillage_type_t

   !> `[soil.tillage]` block container.
   type :: soil_tillage_t
      integer :: i_n_model = 2
      integer :: iRedist   = 2
      type(soil_tillage_event_t), allocatable :: events(:)
      type(soil_tillage_type_t),  allocatable :: types(:)
   end type soil_tillage_t

   type :: soil_config_t
      integer :: swsophy = 0
      integer :: swhyst  = 0
      ! Minimum pressure-head difference (cm) to flip wetting↔drying branch
      ! in the hysteresis routine. Only consumed when swhyst /= 0. Default 0
      ! retains legacy parity (legacy globals were also zero-initialised).
      real(real64) :: tau = 0.0_real64
      integer :: swinco  = 1
      integer :: swmacro = 0
      integer :: swscal  = 0
      integer :: swtill  = 0       !! Tillage-event simulation. 0=off, 1=on (activates flTillage call-site gate; ADR 0020).

      real(real64) :: gwli    = 0.0_real64
      real(real64) :: pondini = 0.0_real64
      real(real64) :: pondmx  = 0.0_real64
      real(real64) :: ksatexm = 0.0_real64
      real(real64) :: rsoil   = 0.0_real64

      ! Surface-runoff drainage resistance (legacy RSRO) and exponent
      ! (legacy RSROEXP). Read by readswap.f90:527-528 right after the
      ! initial-water section. Defaults below match the legacy schema
      ! ranges from rdsdor (0.001..1.0 d, 0.01..10.0 -). Initialize sets
      ! both to 0; the runoff equation `dt/rsro * (pond - pondmx)^rsroexp`
      ! falls back to the instantaneous branch when rsro < 1e-3, so for
      ! parity each case must author finite values matching its .swp.
      real(real64) :: rsro    = 0.0_real64
      real(real64) :: rsroexp = 0.0_real64

      ! Switch for runon (legacy SWRUNON, readswap.f90:1011). When 1, a
      ! companion .inc file feeds runon time-series into the model. Cases
      ! that don't enable runon leave it at 0; the legacy default after
      ! Initialize is also 0 so the schema default agrees.
      integer      :: swrunon = 0

      integer :: reva_top = 0
      integer :: nrstaring = 0  !! 0=user-supplied, 1..6=Staring series

      integer,      allocatable :: sublay(:)
      real(real64), allocatable :: hsublay(:)  !! per-sub-layer height (cm)
      real(real64), allocatable :: hcomp(:)
      integer,      allocatable :: ncomp(:)
      integer,      allocatable :: isoillay(:)
      real(real64), allocatable :: orgmat(:)
      real(real64), allocatable :: bdens(:)
      real(real64), allocatable :: wcontent(:)
      real(real64), allocatable :: cofgen(:,:)
      real(real64), allocatable :: cofani(:)  !! per-soil-layer anisotropy ratio

      type(soil_discretization_t) :: discretization
      type(soil_frost_t)          :: frost
      type(soil_hydraulics_t)     :: hydraulics
      type(soil_initial_t)        :: initial
      type(soil_tillage_t)        :: tillage
   contains
      procedure :: validate => soil_config_validate
      procedure :: finalize => soil_config_finalize
   end type soil_config_t

contains

   subroutine soil_config_validate(self, errors)
      class(soil_config_t),     intent(in)    :: self
      type(error_collection_t), intent(inout) :: errors

      integer :: i

      ! ----- Phase 4f-extend SS-4 stub-error: macropore deferred -----
      ! ADR 0010 keeps macropore_config_t as orphan infrastructure and
      ! the macropore reader code in readswap.f90 unwired. ADR 0011
      ! excludes case 3 (3.macroporeflow) from regression. Authoring
      ! swmacro=1 in a TOML config would let the modern binary copy the
      ! switch into the legacy global without populating any macropore
      ! state, producing silent runtime corruption. Reject it here with
      ! a clear message instead.
      if (self%swmacro == 1) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'soil.swmacro=1 (macropore physics) not yet supported in ' // &
            'the TOML pipeline; case 3 is excluded from regression per ' // &
            'ADR 0011 and the macropore module remains deferred per ' // &
            'ADR 0010.', 'soil')
      end if

      ! ADR 0020 SS-B: swtill=1 is a legitimate value (tillage
      ! activated). Validation accepts it; the call-site gate
      ! (flTillage) controls whether the tillage subsystem actually
      ! runs at runtime. ADR 0021 (candidate, future) will TOML-port
      ! the rest of the Read_Tillage rdinqr/rdsdor block currently
      ! read via TTutil from the staged swap.swp.

      call check_int_enum(self%swsophy, [0, 1],       "soil.swsophy", errors)
      call check_int_enum(self%swhyst,  [0, 1, 2],    "soil.swhyst",  errors)
      call check_int_enum(self%swinco,  [1, 2, 3],    "soil.swinco",  errors)
      call check_int_enum(self%swmacro, [0, 1],       "soil.swmacro", errors)
      call check_int_enum(self%swscal,  [0, 1],       "soil.swscal",  errors)
      call check_int_enum(self%swtill,  [0, 1],       "soil.swtill",  errors)

      call check_nonnegative_real(self%pondmx,  "soil.pondmx",  errors)
      call check_nonnegative_real(self%ksatexm, "soil.ksatexm", errors)
      call check_nonnegative_real(self%rsro,    "soil.rsro",    errors)
      call check_nonnegative_real(self%rsroexp, "soil.rsroexp", errors)
      call check_int_enum(self%swrunon, [0, 1], "soil.swrunon", errors)

      call check_int_range(self%nrstaring, 0, 6, "soil.nrstaring", errors)

      if (allocated(self%cofani)) then
         do i = 1, size(self%cofani)
            call check_real_range(self%cofani(i), 0.01_real64, 100.0_real64, &
                                  "soil.cofani", errors)
         end do
      end if

      call self%discretization%validate(errors)
      call self%frost%validate(errors)
      call self%hydraulics%validate(errors)
      ! soil.initial holds the TOML companion-CSV pathway for SWINCO=3.
      if (self%swinco == 3) then
         call self%initial%validate(errors)
      end if

      ! [soil.tillage] validation — only fires when swtill = 1.
      ! Depth check (events.z vs |zbotcp(NumNod)|) is deferred to the
      ! adapter (apply_soil_tillage) since NumNod isn't known here.
      ! Cross-section date-window check is deferred to a higher-level
      ! validator (Task 4, adapter, which has access to swap_config_t).
      if (self%swtill == 1) then
         call check_int_range(self%tillage%i_n_model, 1, 3, &
                              'soil.tillage.i_n_model', errors)
         call check_int_range(self%tillage%iRedist,   0, 2, &
                              'soil.tillage.iRedist',   errors)

         if (.not. allocated(self%tillage%events) .or. &
             size(self%tillage%events) == 0) then
            call errors%append(ERR_VALIDATION_CROSS_FIELD, &
                               'must be non-empty when soil.swtill = 1', &
                               'soil.tillage.events')
         end if
         if (.not. allocated(self%tillage%types) .or. &
             size(self%tillage%types) == 0) then
            call errors%append(ERR_VALIDATION_CROSS_FIELD, &
                               'must be non-empty when soil.swtill = 1', &
                               'soil.tillage.types')
         end if

         if (allocated(self%tillage%events)) then
            do i = 1, size(self%tillage%events)
               associate (ev => self%tillage%events(i))
                  call check_real_range(ev%intensity, 0.0_real64, 1.0_real64, &
                                        'soil.tillage.events.intensity', errors)
                  if (ev%type_id < 1) then
                     call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
                                        'must be >= 1', &
                                        'soil.tillage.events.type_id')
                  end if
                  if (allocated(self%tillage%types)) then
                     if (.not. any(self%tillage%types%id == ev%type_id)) then
                        call errors%append(ERR_VALIDATION_CROSS_FIELD, &
                                           'type_id not present in soil.tillage.types', &
                                           'soil.tillage.events.type_id')
                     end if
                  end if
                  if (i > 1) then
                     if (ev%date <= self%tillage%events(i - 1)%date) then
                        call errors%append(ERR_VALIDATION_CROSS_FIELD, &
                                           'dates must be strictly ascending', &
                                           'soil.tillage.events.date')
                     end if
                  end if
               end associate
            end do
         end if

         if (allocated(self%tillage%types)) then
            do i = 1, size(self%tillage%types)
               associate (ty => self%tillage%types(i))
                  call check_real_range(ty%rho_cons,    100.0_real64, 3000.0_real64, &
                                        'soil.tillage.types.rho_cons',    errors)
                  call check_real_range(ty%rho_tillage, 100.0_real64, 3000.0_real64, &
                                        'soil.tillage.types.rho_tillage', errors)
                  call check_real_range(ty%k_R,         1.0e-4_real64, 10.0_real64, &
                                        'soil.tillage.types.k_R',         errors)
                  if (self%tillage%i_n_model == 3) then
                     call check_real_range(ty%rho_match, 100.0_real64, 3000.0_real64, &
                                           'soil.tillage.types.rho_match', errors)
                     call check_real_range(ty%N_match,   1.001_real64, 10.0_real64, &
                                           'soil.tillage.types.N_match',   errors)
                  end if
               end associate
            end do
         end if
      end if
   end subroutine soil_config_validate

   !> Per-soil-physical-layer hydraulics validator. All arrays are
   !! optional at the schema level — presence is required only when
   !! swsophy=0, but that cross-section check belongs upstream of the
   !! per-array range checks (the runtime catches missing arrays via
   !! its own fatalerr). Here we just enforce that, when present, every
   !! array has the same length and each cell is within the legacy
   !! `rdfdor` range.
   subroutine soil_hydraulics_validate(self, errors)
      class(soil_hydraulics_t), intent(in)    :: self
      type(error_collection_t), intent(inout) :: errors
      integer :: i, n
      n = -1
      call check_size(self%ores,    n, "soil.hydraulics.ores",    errors)
      call check_size(self%osat,    n, "soil.hydraulics.osat",    errors)
      call check_size(self%alfa,    n, "soil.hydraulics.alfa",    errors)
      call check_size(self%npar,    n, "soil.hydraulics.npar",    errors)
      call check_size(self%ksatfit, n, "soil.hydraulics.ksatfit", errors)
      call check_size(self%lexp,    n, "soil.hydraulics.lexp",    errors)
      call check_size(self%alfaw,   n, "soil.hydraulics.alfaw",   errors)
      call check_size(self%h_enpr,  n, "soil.hydraulics.h_enpr",  errors)
      call check_size(self%ksatexm, n, "soil.hydraulics.ksatexm", errors)
      call check_size(self%bdens,   n, "soil.hydraulics.bdens",   errors)

      if (allocated(self%ores)) then
         do i = 1, size(self%ores)
            call check_real_range(self%ores(i),    0.0_real64,    1.0_real64, &
                                  "soil.hydraulics.ores", errors)
            call check_real_range(self%osat(i),    0.0_real64,    1.0_real64, &
                                  "soil.hydraulics.osat", errors)
            call check_real_range(self%alfa(i),    1.0e-4_real64, 100.0_real64, &
                                  "soil.hydraulics.alfa", errors)
            call check_real_range(self%npar(i),    1.001_real64,  9.0_real64, &
                                  "soil.hydraulics.npar", errors)
            call check_real_range(self%ksatfit(i), 1.0e-5_real64, 1.0e5_real64, &
                                  "soil.hydraulics.ksatfit", errors)
            call check_real_range(self%lexp(i),   -25.0_real64,   25.0_real64, &
                                  "soil.hydraulics.lexp", errors)
            call check_real_range(self%alfaw(i),   1.0e-4_real64, 100.0_real64, &
                                  "soil.hydraulics.alfaw", errors)
            call check_real_range(self%h_enpr(i), -40.0_real64,   0.0_real64, &
                                  "soil.hydraulics.h_enpr", errors)
            call check_real_range(self%ksatexm(i), 1.0e-5_real64, 1.0e5_real64, &
                                  "soil.hydraulics.ksatexm", errors)
            call check_real_range(self%bdens(i),   100.0_real64,  1.0e4_real64, &
                                  "soil.hydraulics.bdens", errors)
         end do
      end if
   end subroutine soil_hydraulics_validate

   !> Helper: on first call (n<0) record the array's size; on subsequent
   !! calls verify every other array matches. Append an error when
   !! sizes diverge so the typical 10-array layer table is sanity-checked.
   subroutine check_size(arr, n, label, errors)
      real(real64), allocatable, intent(in)    :: arr(:)
      integer,                   intent(inout) :: n
      character(len=*),          intent(in)    :: label
      type(error_collection_t),  intent(inout) :: errors
      integer :: this
      if (.not. allocated(arr)) return
      this = size(arr)
      if (n < 0) then
         n = this
      else if (this /= n) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
                            "size mismatch across soil.hydraulics arrays", label)
      end if
   end subroutine check_size

   subroutine soil_discretization_validate(self, errors)
      class(soil_discretization_t), intent(in)    :: self
      type(error_collection_t),     intent(inout) :: errors

      integer :: i

      call check_int_enum(self%swdiscrvert, [0, 1], &
                          "soil.discretization.swdiscrvert", errors)

      if (self%swdiscrvert == 1) then
         call check_int_range(self%numnodnew, 1, 1000, &
                              "soil.discretization.numnodnew", errors)

         if (.not. allocated(self%dznew)) then
            call errors%append(ERR_VALIDATION_CROSS_FIELD, &
                               "dznew required when swdiscrvert=1", &
                               "soil.discretization")
         else
            if (size(self%dznew) /= self%numnodnew) then
               call errors%append(ERR_VALIDATION_CROSS_FIELD, &
                                  "dznew size != numnodnew", &
                                  "soil.discretization")
            end if
            do i = 1, size(self%dznew)
               call check_real_range(self%dznew(i), 1.0e-6_real64, 1000.0_real64, &
                                     "soil.discretization.dznew", errors)
            end do
         end if
      end if
   end subroutine soil_discretization_validate

   subroutine soil_frost_validate(self, errors)
      class(soil_frost_t),      intent(in)    :: self
      type(error_collection_t), intent(inout) :: errors

      call check_int_enum(self%swfrost,  [0, 1], "soil.frost.swfrost",  errors)
      call check_int_enum(self%swsublim, [0, 1], "soil.frost.swsublim", errors)

      if (self%swfrost == 1) then
         call check_real_range(self%tfroststa, -10.0_real64, 0.0_real64, &
                               "soil.frost.tfroststa", errors)
         call check_real_range(self%tfrostend, -10.0_real64, 0.0_real64, &
                               "soil.frost.tfrostend", errors)
         if (self%tfrostend >= self%tfroststa) then
            call errors%append(ERR_VALIDATION_CROSS_FIELD, &
                               "tfrostend >= tfroststa", "soil.frost")
         end if
      end if
   end subroutine soil_frost_validate

   subroutine soil_initial_validate(self, errors)
      class(soil_initial_t),    intent(in)    :: self
      type(error_collection_t), intent(inout) :: errors
      integer :: i

      call check_int_enum(self%swirrigate, [0, 1], "soil.initial.swirrigate", errors)
      call check_real_range(self%ssnow, 0.0_real64, 1000.0_real64, &
                            "soil.initial.ssnow", errors)
      call check_real_range(self%slw,   0.0_real64, 1000.0_real64, &
                            "soil.initial.slw",   errors)
      call check_real_range(self%pond,  0.0_real64,  100.0_real64, &
                            "soil.initial.pond",  errors)
      call check_real_range(self%ldwet, 0.0_real64,  366.0_real64, &
                            "soil.initial.ldwet", errors)
      call check_real_range(self%dt, 1.0e-12_real64, 1.0_real64, &
                            "soil.initial.dt", errors)
      do i = 1, 7
         call check_real_range(self%atmin7(i), -50.0_real64, 50.0_real64, &
                               "soil.initial.atmin7", errors)
      end do
   end subroutine soil_initial_validate

   subroutine soil_config_finalize(self, errors)
      class(soil_config_t),     intent(inout) :: self
      type(error_collection_t), intent(inout) :: errors
   end subroutine soil_config_finalize

end module soil_config_mod
