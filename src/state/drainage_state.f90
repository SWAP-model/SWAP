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
   end type drainage_state_t

end module drainage_state_mod
