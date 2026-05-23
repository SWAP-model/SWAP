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
   end type drainage_state_t

end module drainage_state_mod
