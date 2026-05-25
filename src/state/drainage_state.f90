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

      ! ---- [Piece D] owltab CSV pre-load ----
      ! owltab may already be allocated (the adapter ran config_to_variables
      ! before this init fires). Guard with if (.not. allocated) to avoid
      ! double-allocation; if already allocated the data is already in place.
      if (allocated(config_drain%owltab_file)) then
         block
            use csv_reader_mod, only: read_csv_table
            use error_mod, only: error_collection_t
            real(real64), allocatable :: csv_table(:,:)
            type(error_collection_t)  :: csv_errs
            integer :: lev, nrows, k
            character(len=6) :: hdr(2)
            hdr(1) = 'date  '
            hdr(2) = 'level '
            ! Allocate only if not already done (Piece A guard covers this).
            if (.not. allocated(self%owltab)) then
               allocate(self%owltab(config_drain%nrlevs, 2*MAOWL))
               self%owltab = 0.0_real64
            end if
            do lev = 1, size(config_drain%owltab_file)
               if (len_trim(config_drain%owltab_file(lev)) == 0) cycle
               call read_csv_table(trim(config_drain%owltab_file(lev)), hdr, csv_table, csv_errs)
               call csv_errs%abort_if_fatal()
               if (csv_errs%count() == 0) then
                  nrows = size(csv_table, 1)
                  self%nowltab(lev) = nrows
                  do k = 1, nrows
                     self%owltab(lev, 2*k-1) = csv_table(k, 1)
                     self%owltab(lev, 2*k)   = csv_table(k, 2)
                  end do
               end if
            end do
         end block
      end if

      ! ---- [Piece E] surface_runoff sub-section ----
      self%cofintfl = config_drain%surface_runoff%cofintfl
      self%expintfl = config_drain%surface_runoff%expintfl
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

   end subroutine drainage_state_init

end module drainage_state_mod
