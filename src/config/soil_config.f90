!> [soil] section config.
module soil_config_mod
   use iso_fortran_env, only: real64
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   use validation_mod, only: check_int_enum, check_int_range, &
                             check_real_range, check_nonnegative_real
   implicit none
   private

   public :: soil_config_t
   public :: soil_discretization_t
   public :: soil_frost_t
   public :: soil_hydraulics_t

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

   type :: soil_config_t
      integer :: swsophy = 0
      integer :: swhyst  = 0
      integer :: swinco  = 1
      integer :: swmacro = 0
      integer :: swscal  = 0

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

      ! Phase 4f Task B5: legacy SWINCO=3 reads initial state (h, cml,
      ! ssnow, slw, pond, Tsoil) from a previous-run .end-style file
      ! named here. Adapter reads the file directly when allocated.
      character(len=:), allocatable :: inifil

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
   contains
      procedure :: validate => soil_config_validate
      procedure :: finalize => soil_config_finalize
   end type soil_config_t

contains

   subroutine soil_config_validate(self, errors)
      class(soil_config_t),     intent(in)    :: self
      type(error_collection_t), intent(inout) :: errors

      integer :: i

      call check_int_enum(self%swsophy, [0, 1],       "soil.swsophy", errors)
      call check_int_enum(self%swhyst,  [0, 1, 2],    "soil.swhyst",  errors)
      call check_int_enum(self%swinco,  [1, 2, 3],    "soil.swinco",  errors)
      call check_int_enum(self%swmacro, [0, 1],       "soil.swmacro", errors)
      call check_int_enum(self%swscal,  [0, 1],       "soil.swscal",  errors)

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

   subroutine soil_config_finalize(self, errors)
      class(soil_config_t),     intent(inout) :: self
      type(error_collection_t), intent(inout) :: errors
   end subroutine soil_config_finalize

end module soil_config_mod
