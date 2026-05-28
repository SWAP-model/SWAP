!> @file drainage_state.f90
!! Typed state record for the drainage subsystem. Owns the 6
!! drainage-flux variables that legacy SWAP held in the
!! `variables.f90` globals module.
!!
!! `qdra` and `qdrain` were temporarily classified into
!! surfacewater_state_t during the surface-water migration pilot
!! (ADR 0030) because of write-site overlap. The drainage migration
!! arc moves them here, where they architecturally belong.
!!
!! See ADR 0031 (state-migration drainage subsystem) and the
!! 2026-05-10 design spec.

module drainage_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   use swap_array_dimensions, only: MADR
   implicit none
   private
   public :: drainage_state_t

   type :: drainage_state_t
      ! Scalar fluxes
      real(real64) :: qdrd     = 0.0_real64    ! drain-direction sub-flux (cm/d)

      ! Per-level (Madr-sized) arrays — allocated by drainage_init
      real(real64), allocatable :: qdrain(:)       ! lateral drainage flux per level (cm/d)
      real(real64), allocatable :: drainl(:)       ! drain length per level (cm)
      real(real64), allocatable :: wetper(:)       ! wetted perimeter per level (cm)
      real(real64), allocatable :: ztopdislay(:)   ! top of discharge layer per level (cm)

      ! Per-level / per-compartment array (Madr × macp)
      real(real64), allocatable :: qdra(:,:)       ! lateral drainage flux per level/compartment (cm/d)

      ! [SS-GR-BH A8] drainage geometry + switches formerly in variables.f90
      integer :: nrlevs      = 0
      integer :: swdivd      = 0
      integer :: swnrsrf     = 0
      integer :: swtopnrsrf  = 0
      integer :: swdivdinf   = 0
      real(real64) :: FacDpthInf = 0.0_real64
      real(real64), allocatable :: L(:)       !! drainage spacing per level [cm]
      real(real64), allocatable :: zbotdr(:)  !! drainage depth per level [cm]
      real(real64), allocatable :: owltab(:,:)  !! open-water level table (level × 2*maowl time-value pairs)

      ! Hooghoudt/Ernst geometry scalars (dramet=2), snapshotted from
      ! config%drain at config_to_variables. ipos selects the Ernst case
      ! and gates which of {khbot, zintf, kvtop, kvbot, geofac} are
      ! meaningful — see drainage.f90:260-323.
      real(real64) :: basegw = 0.0_real64  !! [cm] base of saturated zone (zimp floor)
      real(real64) :: entres = 0.0_real64  !! [d] drain entry resistance
      real(real64) :: geofac = 0.0_real64  !! [-] Ernst geometry factor (ipos=5)
      integer      :: ipos   = 0           !! 1..5 Ernst position-of-drain selector
      real(real64) :: khtop  = 0.0_real64  !! [cm/d] horizontal K top layer
      real(real64) :: khbot  = 0.0_real64  !! [cm/d] horizontal K bottom layer (ipos>=3)
      real(real64) :: kvtop  = 0.0_real64  !! [cm/d] vertical K top layer (ipos>=4)
      real(real64) :: kvbot  = 0.0_real64  !! [cm/d] vertical K bottom layer (ipos>=4)
      real(real64) :: zintf  = 0.0_real64  !! [cm] fine/coarse interface depth (ipos>=3)

      ! Drainage surface-runoff config (snapshotted from
      ! config%drain%surface_runoff at config_to_variables).
      real(real64) :: cofintfl     = 0.0_real64  !! interflow coefficient
      real(real64) :: expintfl     = 0.0_real64  !! interflow exponent
      integer      :: NumLevRapDra = 0           !! rapid drainage drain-level index

      ! Drainage method/level config switches (snapshotted at config_to_variables).
      integer :: dramet    = 0   !! drainage method (1=table, 2=Hooghoudt/Ernst, 3=ext.drainage)
      integer :: swliminf  = 0   !! limit infiltration head to channel water depth (0/1)
      integer :: swdislay  = 0   !! distributed-layer drainage flag (0/1/2)

      ! Per-drain-level config arrays (sized MADR, broadcast or copied
      ! from config). Guarded-alloc — config_to_variables runs first.
      integer,      allocatable :: swdtyp(:)       !! drain type per level (0=open, 1=closed, 2=q-h)
      integer,      allocatable :: swallo(:)       !! allow-drainage flag per level (1/2/3)
      integer,      allocatable :: swtopdislay(:)  !! top-discharge-layer switch per level (0/1)
      real(real64), allocatable :: ftopdislay(:)   !! top-discharge-layer factor per level [0..1]

      ! Drainage-flux-vs-groundwater-level table (dramet=1 branch).
      ! Dormant — no TOML writer; reader at drainage.f90:389 (afgen, 50 pairs).
      real(real64) :: qdrtab(50) = 0.0_real64

      ! Per-drain-level count of open-water-level table entries (owltab(:,1:2*nowltab(lev))).
      ! Written by config_to_variables when an owltab CSV is supplied; consumed by
      ! divdra + drainage for afgen extents.
      integer :: nowltab(MADR) = 0

      ! Rapid drainage flux (cm/d). Dormant — kept at zero in current pipeline
      ! (legacy [retired-zero] tag). Consumed by surfacewater task=2/3 and waterbalance.
      real(real64) :: QRapDra = 0.0_real64

      ! Per-level drainage resistance arrays (cluster 4 retirement).
      ! Snapshotted from config%drain at init so that compute_drainage_level_flux
      ! and bocodrb no longer need state%cfg%drain access.
      real(real64), allocatable :: drares(:)   !! drainage resistance per level (d)
      real(real64), allocatable :: infres(:)   !! infiltration resistance per level (d)
      real(real64), allocatable :: widthr(:)   !! channel bottom width per level (cm)
      real(real64), allocatable :: taludr(:)   !! channel side slope per level (-)
      real(real64), allocatable :: gwlinf(:)   !! gwl limit for infiltration per level (cm)
      real(real64), allocatable :: rdrain(:)   !! drainage resistance per level (d) — ext.drain
      real(real64), allocatable :: rinfi(:)    !! infiltration resistance per level (d) — ext.drain
      real(real64), allocatable :: rentry(:)   !! entry resistance per level (d) — ext.drain
      real(real64), allocatable :: rexit(:)    !! exit resistance per level (d) — ext.drain

      ! Scalar drainage config snapshots (cluster 4 retirement).
      real(real64) :: shape        = 0.0_real64  !! trapezoidal shape factor
      real(real64) :: rsurfdeep    = 0.0_real64  !! surface resistance deep (d)
      real(real64) :: rsurfshallow = 0.0_real64  !! surface resistance shallow (d)

   contains
      procedure :: init => drainage_state_init
   end type drainage_state_t

contains

   !> One-time runtime initialization for drainage state.
   !! [GR-SEED 2026-05-25 Task 7] Promoted from free subroutine drainage_init
   !! and absorbs: drainl/wetper/qdrain/qdra/ztopdislay allocation + zeroing
   !! (legacy drainage_init body); scalar seeding (dramet/swdislay/basegw/entres);
   !! DRAMET=2 ipos chain; per-level swdtyp/swallo arrays; owltab CSV pre-load;
   !! surface_runoff sub-section (cofintfl/expintfl/geofac/NumLevRapDra/
   !! swtopdislay/ftopdislay).
   subroutine drainage_state_init(self, config_drain, numnod)
      use drainage_config_mod,   only: drainage_config_t
      use swap_array_dimensions, only: MADR, MAOWL
      use, intrinsic :: iso_fortran_env, only: real64
      class(drainage_state_t),  intent(inout) :: self
      type(drainage_config_t),  intent(in)    :: config_drain
      integer,                  intent(in)    :: numnod
      integer :: i
      integer :: nr

      ! ---- [Piece A] Legacy drainage_init body — allocation + zeroing ----
      nr = config_drain%nrlevs
      if (.not. allocated(self%qdrain))     allocate(self%qdrain(nr))
      if (.not. allocated(self%drainl))     allocate(self%drainl(nr))
      if (.not. allocated(self%wetper))     allocate(self%wetper(nr))
      if (.not. allocated(self%ztopdislay)) allocate(self%ztopdislay(nr))
      if (.not. allocated(self%qdra))       allocate(self%qdra(nr, numnod))
      if (.not. allocated(self%L))          allocate(self%L(nr))
      if (.not. allocated(self%zbotdr))     allocate(self%zbotdr(nr))
      ! owltab: guard here; Piece D may have already allocated via CSV pre-load.
      if (.not. allocated(self%owltab)) then
         allocate(self%owltab(nr, 2*MAOWL))
         self%owltab = 0.0_real64
      end if

      self%drainl     = 0.0_real64
      self%wetper     = 0.0_real64
      self%ztopdislay = 0.0_real64
      self%qdrd       = 0.0_real64
      if (config_drain%dramet == 2) then
         self%wetper(1) = config_drain%wetper
      end if
      self%qdrain = 0.0_real64
      self%qdra   = 0.0_real64

      ! ---- [Piece B] Scalar seeding + DRAMET=2 ipos chain ----
      self%nrlevs   = config_drain%nrlevs
      self%swdivd   = config_drain%swdivd
      self%swnrsrf  = config_drain%surface_runoff%swnrsrf
      self%dramet   = config_drain%dramet
      self%swdislay = config_drain%swdislay
      self%basegw   = config_drain%basegw
      self%entres   = config_drain%entres

      ! DRAMET=2 (Hooghoudt/Ernst) geometry scalars.
      ! ADR 0031 Phase 2 Task 5: wetper(1) already seeded in Piece A above.
      if (config_drain%dramet == 2) then
         self%ipos  = config_drain%ipos
         self%khtop = config_drain%khtop
         if (config_drain%ipos >= 3) then
            self%khbot = config_drain%khbot
            self%zintf = config_drain%zintf
         end if
         if (config_drain%ipos >= 4) then
            self%kvtop = config_drain%kvtop
            self%kvbot = config_drain%kvbot
         end if
         if (config_drain%ipos == 5) then
            self%geofac = config_drain%geofac
         end if
      end if

      ! ---- [Piece C] Per-level swdtyp + swallo arrays ----
      if (allocated(config_drain%swdtyp)) then
         if (.not. allocated(self%swdtyp)) then
            allocate(self%swdtyp(size(config_drain%swdtyp)))
            self%swdtyp = 0
         end if
         do i = 1, size(config_drain%swdtyp)
            self%swdtyp(i) = config_drain%swdtyp(i)
         end do
      end if
      if (allocated(config_drain%swallo)) then
         if (.not. allocated(self%swallo)) then
            allocate(self%swallo(size(config_drain%swallo)))
            self%swallo = 0
         end if
         do i = 1, size(config_drain%swallo)
            self%swallo(i) = config_drain%swallo(i)
         end do
      end if
      self%swliminf = config_drain%swliminf

      ! ---- [Piece D] owltab CSV pre-load (typed loader — ADR 0044 Family 1) ----
      ! owl_events_table_t validates strictly-ascending dates (required by afgen).
      ! One fresh error_collection_t per level so errors don't accumulate across
      ! levels (per meteo-pilot lesson: don't share an accumulator across loads).
      ! Legacy interleaved layout preserved: owltab(lev, 2k-1) = date,
      ! owltab(lev, 2k) = level; nowltab(lev) = nrows for afgen extents.
      if (allocated(config_drain%owltab_file)) then
         block
            use drainage_csv_mod, only: owl_events_table_t
            use error_mod,        only: error_collection_t
            type(owl_events_table_t) :: owl_tbl
            type(error_collection_t) :: lev_errs
            integer :: lev, nrows, k
            do lev = 1, size(config_drain%owltab_file)
               if (len_trim(config_drain%owltab_file(lev)) == 0) cycle
               lev_errs = error_collection_t()
               call owl_tbl%load(trim(config_drain%owltab_file(lev)), lev_errs)
               call lev_errs%abort_if_fatal()
               if (owl_tbl%is_loaded) then
                  nrows = size(owl_tbl%rows)
                  ! Cap at MAOWL pairs (afgen table bound).
                  if (nrows > MAOWL) nrows = MAOWL
                  self%nowltab(lev) = nrows
                  do k = 1, nrows
                     self%owltab(lev, 2*k-1) = owl_tbl%rows(k)%date
                     self%owltab(lev, 2*k)   = owl_tbl%rows(k)%level
                  end do
               end if
            end do
         end block
      end if

      ! ---- [Piece E] surface_runoff sub-section ----
      self%cofintfl    = config_drain%surface_runoff%cofintfl
      self%expintfl    = config_drain%surface_runoff%expintfl
      self%swtopnrsrf  = config_drain%surface_runoff%swtopnrsrf
      self%swdivdinf   = config_drain%surface_runoff%swdivdinf
      self%FacDpthInf  = config_drain%surface_runoff%facdpthinf
      ! ADR 0031: gate the surface_runoff geofac write to avoid overwriting
      ! the ipos==5 (Ernst geometry factor) write from Piece B above.
      if (config_drain%ipos /= 5) then
         self%geofac = config_drain%surface_runoff%geofac
      end if
      self%NumLevRapDra = config_drain%surface_runoff%numlevrapdra

      ! swtopdislay(MADR) + ftopdislay(MADR): broadcast scalar config field to all levels.
      if (.not. allocated(self%swtopdislay)) then
         allocate(self%swtopdislay(MADR))
         self%swtopdislay = 0
      end if
      if (.not. allocated(self%ftopdislay)) then
         allocate(self%ftopdislay(MADR))
         self%ftopdislay = 0.0_real64
      end if
      do i = 1, size(self%swtopdislay)
         self%swtopdislay(i) = config_drain%surface_runoff%swtopdislay
      end do
      do i = 1, size(self%ftopdislay)
         self%ftopdislay(i) = config_drain%surface_runoff%ftopdislay
      end do

      ! ---- [Piece G] Per-level drainage resistance arrays (cluster 4 retirement) ----
      ! Snapshotted from config_drain so that compute_drainage_level_flux and bocodrb
      ! no longer need state%cfg%drain.
      if (allocated(config_drain%drares)) then
         if (.not. allocated(self%drares)) allocate(self%drares(nr))
         self%drares(1:size(config_drain%drares)) = config_drain%drares
      end if
      if (allocated(config_drain%infres)) then
         if (.not. allocated(self%infres)) allocate(self%infres(nr))
         self%infres(1:size(config_drain%infres)) = config_drain%infres
      end if
      if (allocated(config_drain%widthr)) then
         if (.not. allocated(self%widthr)) allocate(self%widthr(nr))
         self%widthr(1:size(config_drain%widthr)) = config_drain%widthr
      end if
      if (allocated(config_drain%taludr)) then
         if (.not. allocated(self%taludr)) allocate(self%taludr(nr))
         self%taludr(1:size(config_drain%taludr)) = config_drain%taludr
      end if
      if (allocated(config_drain%gwlinf)) then
         if (.not. allocated(self%gwlinf)) allocate(self%gwlinf(nr))
         self%gwlinf(1:size(config_drain%gwlinf)) = config_drain%gwlinf
      end if
      if (allocated(config_drain%rdrain)) then
         if (.not. allocated(self%rdrain)) allocate(self%rdrain(nr))
         self%rdrain(1:size(config_drain%rdrain)) = config_drain%rdrain
      end if
      if (allocated(config_drain%rinfi)) then
         if (.not. allocated(self%rinfi)) allocate(self%rinfi(nr))
         self%rinfi(1:size(config_drain%rinfi)) = config_drain%rinfi
      end if
      if (allocated(config_drain%rentry)) then
         if (.not. allocated(self%rentry)) allocate(self%rentry(nr))
         self%rentry(1:size(config_drain%rentry)) = config_drain%rentry
      end if
      if (allocated(config_drain%rexit)) then
         if (.not. allocated(self%rexit)) allocate(self%rexit(nr))
         self%rexit(1:size(config_drain%rexit)) = config_drain%rexit
      end if
      self%shape        = config_drain%shape
      self%rsurfdeep    = config_drain%surface_runoff%rsurfdeep
      self%rsurfshallow = config_drain%surface_runoff%rsurfshallow

      ! ---- [Piece F] L / zbotdr seeding (orchestrator-dissolution arc step 4) ----
      ! dramet=2: scalar lm (m) converted to cm for L(1); zbotdr_basic for zbotdr(1).
      ! Otherwise: copy from per-level config arrays when allocated.
      if (config_drain%dramet == 2) then
         self%L(1)      = 100.0d0 * config_drain%lm
         self%zbotdr(1) = config_drain%zbotdr_basic
      else
         if (allocated(config_drain%L)) &
            self%L(1:size(config_drain%L)) = config_drain%L
         if (allocated(config_drain%zbotdr)) &
            self%zbotdr(1:size(config_drain%zbotdr)) = config_drain%zbotdr
      end if

   end subroutine drainage_state_init

end module drainage_state_mod
