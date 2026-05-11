!> @file soilwater_state.f90
!! Typed state record for the soil-water boundary subsystem (ADR 0035)
!! and crop water uptake subsystem (ADR 0036).
!!
!! Boundary subset (12 instantaneous scalars — ADR 0035):
!!   Top-boundary fields (from boundtop):
!!     qtop, reva, hsurf, runots, QMpLatSs, ftoph, FlRunoff
!!   Bottom-boundary fields (from BoundBottom):
!!     qbot, qbot_nonfrozen, hbot, gwlinp, deepgw
!!
!! Crop water uptake subset (22 fields — ADR 0036):
!!   Per-node arrays (allocated by soilwater_init in C-1.2):
!!     Primary Feddes+stress path: qrot, qpotrot, qredwet, qreddry, qredsol, qredfrs
!!     JvL microscopic path: mflux, mroot, hroot, rootrho, rootphi, rmax
!!     Init-once lookup table: mfluxtable
!!   Scalars: qrosum, qredwetsum, qreddrysum, qredsolsum, qredfrssum
!!   Flag: flWrtNonox
!!   JvL scalars: Tactual, alpJvLier, hleaf, Hxylem
!!
!! Excluded (config, not runtime state):
!!   - swbotb, swqhbot, swtopb — boundary-condition switches; config flags.
!!   - pond, rsro, rsroexp — surface runoff parameters; config-driven.
!!   - cumulative/intermediate fields (inqrot, iqrot, iqredXXX, cqrot) —
!!     soil-water-core arc territory.
!!
!! No type-bound reset() procedure — all fields are instantaneous; the
!! implementer simply overwrites them each step.
!!
!! See docs/superpowers/specs/2026-05-10-state-migration-boundary-design.md (ADR 0035)
!!     docs/superpowers/specs/2026-05-10-state-migration-crop-uptake-design.md (ADR 0036)

module soilwater_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   public :: soilwater_state_t, soilwater_init

   type :: soilwater_state_t

      ! Top-boundary fluxes / surface variables (from boundtop)

      real(real64) :: qtop     = 0.0_real64  !! top-surface flux (cm/d)
      real(real64) :: reva     = 0.0_real64  !! actual soil evaporation (cm/d)
      real(real64) :: hsurf    = 0.0_real64  !! pressure head at surface (cm)
      real(real64) :: runots   = 0.0_real64  !! runoff this step (cm)
      real(real64) :: QMpLatSs = 0.0_real64  !! lateral macropore inflow (cm/d)
      logical      :: ftoph    = .false.     !! flag: pressure-head top boundary
      logical      :: FlRunoff = .false.     !! flag: runoff potential

      ! Bottom-boundary fluxes / variables (from BoundBottom)

      real(real64) :: qbot           = 0.0_real64  !! bottom flux (cm/d)
      real(real64) :: qbot_nonfrozen = 0.0_real64  !! bottom flux pre-frost snapshot
      real(real64) :: hbot           = 0.0_real64  !! prescribed head at bottom (cm)
      real(real64) :: gwlinp         = 0.0_real64  !! prescribed gwl, swbotb=1 (cm)
      real(real64) :: deepgw         = 0.0_real64  !! deep-aquifer head, swbotb=3 (cm)

      ! ===========================================================================
      ! CROP WATER UPTAKE (22 fields — ADR 0036, 2026-05-11)
      ! ===========================================================================
      ! Per-node arrays (allocated by soilwater_init in C-1.2; unallocated here):

      ! Primary Feddes+stress path (6 per-node arrays):
      real(real64), allocatable :: qrot(:)      !! per-node root sink term (cm/d)
      real(real64), allocatable :: qpotrot(:)   !! per-node potential uptake before stress reduction (cm/d)
      real(real64), allocatable :: qredwet(:)   !! per-node wet-stress reduction (cm/d)
      real(real64), allocatable :: qreddry(:)   !! per-node dry-stress reduction (cm/d)
      real(real64), allocatable :: qredsol(:)   !! per-node salt-stress reduction (cm/d)
      real(real64), allocatable :: qredfrs(:)   !! per-node frost-stress reduction (cm/d)

      ! JvL microscopic uptake path (6 per-node arrays; active when swdrought=2):
      real(real64), allocatable :: mflux(:)     !! matric-flux potential per node (cm²/d)
      real(real64), allocatable :: mroot(:)     !! matric-flux potential at root surface (cm²/d)
      real(real64), allocatable :: hroot(:)     !! pressure head at root surface (cm)
      real(real64), allocatable :: rootrho(:)   !! root density per node (cm/cm³)
      real(real64), allocatable :: rootphi(:)   !! root geometry factor per node (-)
      real(real64), allocatable :: rmax(:)      !! maximum radial uptake per node (cm/d)

      ! Init-once lookup table (1 per-layer × 801 table; active when swdrought=2):
      real(real64), allocatable :: mfluxtable(:,:)  !! matric-flux lookup (nlay × 801)

      ! Primary scalars (column sums over the Feddes+stress arrays):
      real(real64) :: qrosum       = 0.0_real64  !! column-sum root uptake (cm/d)
      real(real64) :: qredwetsum   = 0.0_real64  !! column-sum wet-stress reduction (cm/d)
      real(real64) :: qreddrysum   = 0.0_real64  !! column-sum dry-stress reduction (cm/d)
      real(real64) :: qredsolsum   = 0.0_real64  !! column-sum salt-stress reduction (cm/d)
      real(real64) :: qredfrssum   = 0.0_real64  !! column-sum frost-stress reduction (cm/d)

      ! Flag:
      logical :: flWrtNonox = .false.  !! non-oxygen stress override flag (swWrtNonox=1)

      ! JvL scalars (active when swdrought=2):
      real(real64) :: Tactual   = 0.0_real64  !! prior-step actual transpiration for JvL initial guess (cm/d)
      real(real64) :: alpJvLier = 0.0_real64  !! alpha factor from JvL solve (-)
      real(real64) :: hleaf     = 0.0_real64  !! leaf water potential (cm); last-iterate cache
      real(real64) :: Hxylem    = 0.0_real64  !! xylem water potential (cm)

   end type soilwater_state_t

contains

   !> Lifecycle init for soilwater typed state.
   !! Zeros/resets all scalar fields in both the boundary subset (ADR 0035)
   !! and the crop-uptake scalar subset (ADR 0036).  Per-node allocatable
   !! arrays (qrot, qpotrot, qredwet, qreddry, qredsol, qredfrs, mflux,
   !! mroot, hroot, rootrho, rootphi, rmax, mfluxtable) are NOT allocated
   !! here — that is deferred to C-1.2 when the signature gains numnod/nlay.
   !!
   !! Takes soilwater_state_t directly (not swap_state_t) to avoid a circular
   !! dependency: soilwater_state_mod is used by swap_state_mod.
   !! Mirrors heat_init pattern but at the sub-record level.
   !! Called from swap.f90 immediately after CalcGrid(), before DoTillage(1).
   !!
   !! Design: docs/superpowers/specs/2026-05-10-state-migration-boundary-design.md D8
   !!         docs/superpowers/specs/2026-05-10-state-migration-crop-uptake-design.md
   subroutine soilwater_init(sw)
      type(soilwater_state_t), intent(inout) :: sw

      ! Top-boundary fields
      sw%qtop      = 0.0_real64
      sw%reva      = 0.0_real64
      sw%hsurf     = 0.0_real64
      sw%runots    = 0.0_real64
      sw%QMpLatSs  = 0.0_real64
      sw%ftoph     = .false.
      sw%FlRunoff  = .false.

      ! Bottom-boundary fields
      sw%qbot           = 0.0_real64
      sw%qbot_nonfrozen = 0.0_real64
      sw%hbot           = 0.0_real64
      sw%gwlinp         = 0.0_real64
      sw%deepgw         = 0.0_real64

      ! Crop-uptake scalars (ADR 0036 — C-1.1)
      ! Per-node arrays are allocated in C-1.2; scalars zeroed here.
      sw%qrosum      = 0.0_real64
      sw%qredwetsum  = 0.0_real64
      sw%qreddrysum  = 0.0_real64
      sw%qredsolsum  = 0.0_real64
      sw%qredfrssum  = 0.0_real64
      sw%flWrtNonox  = .false.
      sw%Tactual     = 0.0_real64
      sw%alpJvLier   = 0.0_real64
      sw%hleaf       = 0.0_real64
      sw%Hxylem      = 0.0_real64
   end subroutine soilwater_init

end module soilwater_state_mod
