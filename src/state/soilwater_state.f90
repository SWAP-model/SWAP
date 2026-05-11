!> @file soilwater_state.f90
!! Typed state record for the soil-water boundary subsystem (ADR 0035),
!! crop water uptake subsystem (ADR 0036), and soil-water core (ADR 0038).
!!
!! Boundary subset (12 instantaneous scalars — ADR 0035):
!!   Top-boundary fields (from boundtop):
!!     qtop, reva, hsurf, runots, QMpLatSs, ftoph, FlRunoff
!!   Bottom-boundary fields (from BoundBottom):
!!     qbot, qbot_nonfrozen, hbot, gwlinp, deepgw
!!
!! Crop water uptake subset (22 fields — ADR 0036):
!!   Per-node arrays (allocated by soilwater_init, C-1.2):
!!     Primary Feddes+stress path: qrot, qpotrot, qredwet, qreddry, qredsol, qredfrs
!!     JvL microscopic path: mflux, mroot, hroot, rootrho, rootphi, rmax
!!     Init-once lookup table: mfluxtable
!!   Scalars: qrosum, qredwetsum, qreddrysum, qredsolsum, qredfrssum
!!   Flag: flWrtNonox
!!   JvL scalars: Tactual, alpJvLier, hleaf, Hxylem
!!
!! Soil-water core subset (74 new fields — ADR 0038, S-1.1):
!!   Flat instantaneous — 17 per-node/per-layer allocatable arrays:
!!     theta, thetm1, thetar, thetas, h, hm1, q, k, kmean, dimoca,
!!     cofgen(21,:), FrArMtrx, fluseksatexm (logical), indeks (integer),
!!     evp, thetsl (per-layer)
!!   Flat instantaneous — 13 scalars/logicals:
!!     pond, pondm1, pondini, gwl, gwlm1, nodgwl, pegwl, bpegwl, npegwl,
!!     gwlflcpzo, nodgwlflcpzo, hatm, volact, volm1, volini, wbalance,
!!     runon, fllowgwl
!!   Intermediate cohort (soilwater_intermediate_t :: intr):
!!     22 non-per-day fields (6 arrays + 16 scalars) + 8 per-day fields
!!     reset() = flzerointr; reset_per_day() = flDayStart
!!   Cumulative cohort (soilwater_cumulative_t :: cumu):
!!     14 scalar fields; reset() = flzerocumu
!!
!! Excluded (config, not runtime state):
!!   - swbotb, swqhbot, swtopb — boundary-condition switches; config flags.
!!   - pond, rsro, rsroexp — surface runoff parameters; config-driven.
!!   - Grid dimensions (numnod, numlay, dz, z, disnod, …) stay legacy globals
!!     (heat ADR 0034 precedent; future grid_t arc territory).
!!   - cQMpLatSs — macropore arc territory (ADR 0038 D11).
!!
!! See docs/superpowers/specs/2026-05-10-state-migration-boundary-design.md (ADR 0035)
!!     docs/superpowers/specs/2026-05-10-state-migration-crop-uptake-design.md (ADR 0036)
!!     docs/superpowers/specs/2026-05-11-state-migration-soilwater-core-design.md (ADR 0038)

module soilwater_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   public :: soilwater_state_t, soilwater_init
   public :: soilwater_intermediate_t
   public :: soilwater_cumulative_t

   ! ---------------------------------------------------------------------------
   !> Intermediate accumulators — reset when flzerointr fires (reset()),
   !! or when flDayStart fires (reset_per_day() for the 8 per-day fields).
   !!
   !! 22 non-per-day fields:
   !!   6 allocatable arrays (per-node, allocated in soilwater_init S-1.2):
   !!     inq, inqrot, inqssdi, iqdo, iqup, IThetaBeg
   !!   16 scalars:
   !!     iqrot, iqssdi, iqredwet, iqreddry, iqredsol, iqredfrs,
   !!     ies0, iet0, iew0, iintc, iruno, irunoCN, irunon,
   !!     iqbot, iqtdo, iqtup, IPondBeg, iprec, igird, inird
   !!
   !! 8 per-day fields (reset_per_day() only; reset() also touches these):
   !!   6 per-day scalars: tra, iqredwet_day, iqreddry_day,
   !!                      iqredsol_day, iqredfrs_day, iptra_day
   !!   2 per-day arrays:  qpotrot_day(:), qredtot_day(:)
   !!
   !! Mild ADR 0033 extension: same cohort type, two distinct reset procedures,
   !! two distinct activity gates (flzerointr vs flDayStart).
   !! Mirrors atmosphere_intermediate_t (ADR 0037) pattern.
   ! ---------------------------------------------------------------------------
   type :: soilwater_intermediate_t

      ! Per-node flux accumulators (6 allocatable arrays; size numnod or numnod+1)
      ! Allocated by soilwater_init (S-1.2). Unallocated until then.
      real(real64), allocatable :: inq(:)       !< intra-period inter-comp flux (cm)
      real(real64), allocatable :: inqrot(:)    !< intra-period root uptake (cm)
      real(real64), allocatable :: inqssdi(:)   !< intra-period SSDI flux (cm)
      real(real64), allocatable :: iqdo(:)      !< intra-period downward flux (cm)
      real(real64), allocatable :: iqup(:)      !< intra-period upward flux (cm)
      real(real64), allocatable :: IThetaBeg(:) !< theta at start of intr period (-)

      ! Non-per-day scalars (20; zeroed by both reset() and reset_per_day()-indirectly)
      real(real64) :: iqrot     = 0.0_real64   !< period root uptake sum (cm)
      real(real64) :: iqssdi    = 0.0_real64   !< period SSDI flux sum (cm)
      real(real64) :: iqredwet  = 0.0_real64   !< period wet-stress reduction (cm)
      real(real64) :: iqreddry  = 0.0_real64   !< period dry-stress reduction (cm)
      real(real64) :: iqredsol  = 0.0_real64   !< period salt-stress reduction (cm)
      real(real64) :: iqredfrs  = 0.0_real64   !< period frost-stress reduction (cm)
      real(real64) :: ies0      = 0.0_real64   !< period reference soil evap (cm)
      real(real64) :: iet0      = 0.0_real64   !< period reference transpiration (cm)
      real(real64) :: iew0      = 0.0_real64   !< period reference evaporation (cm)
      real(real64) :: iintc     = 0.0_real64   !< period interception (cm)
      real(real64) :: iruno     = 0.0_real64   !< period runoff (cm)
      real(real64) :: irunoCN   = 0.0_real64   !< period CN runoff (cm)
      real(real64) :: irunon    = 0.0_real64   !< period runon (cm)
      real(real64) :: iqbot     = 0.0_real64   !< period bottom flux (cm)
      real(real64) :: iqtdo     = 0.0_real64   !< period total downward flux (cm)
      real(real64) :: iqtup     = 0.0_real64   !< period total upward flux (cm)
      real(real64) :: IPondBeg  = 0.0_real64   !< ponding at start of intr period (cm)
      real(real64) :: iprec     = 0.0_real64   !< period precipitation (cm)
      real(real64) :: igird     = 0.0_real64   !< period irrigation gross (cm)
      real(real64) :: inird     = 0.0_real64   !< period irrigation net (cm)

      ! Per-day cohort (flDayStart gate — zeroed by reset_per_day() and reset())
      real(real64) :: tra           = 0.0_real64  !< daily actual transpiration (cm)
      real(real64) :: iqredwet_day  = 0.0_real64  !< per-day wet-stress reduction (cm)
      real(real64) :: iqreddry_day  = 0.0_real64  !< per-day dry-stress reduction (cm)
      real(real64) :: iqredsol_day  = 0.0_real64  !< per-day salt-stress reduction (cm)
      real(real64) :: iqredfrs_day  = 0.0_real64  !< per-day frost-stress reduction (cm)
      real(real64) :: iptra_day     = 0.0_real64  !< per-day potential transpiration (cm)
      real(real64), allocatable :: qpotrot_day(:) !< per-day potential uptake per node (cm)
      real(real64), allocatable :: qredtot_day(:) !< per-day total reduction per node (cm)

   contains
      procedure :: reset         => soilwater_intermediate_reset      !< zeroes all 22+8 fields (flzerointr)
      procedure :: reset_per_day => soilwater_per_day_reset           !< zeroes only 8 per-day fields (flDayStart)
   end type soilwater_intermediate_t

   ! ---------------------------------------------------------------------------
   !> Cumulative accumulators — reset when flzerocumu fires (reset()).
   !!
   !! 14 scalar fields:
   !!   cqssdi, cqrot, cqbot, cqbotdo, cqbotup, cinund, crunon, crunoff,
   !!   crunoffCN, cqtdo, cqtup, cqprai, cgird, cnird
   !!
   !! No allocatable arrays. All scalars; no allocated() guards needed in reset().
   !! cgird/cnird: multi-owner with irrigation arc (ADR 0038 D14); irrigation.f90
   !! continues subset-reset until irrigation arc takes ownership.
   ! ---------------------------------------------------------------------------
   type :: soilwater_cumulative_t
      real(real64) :: cqssdi    = 0.0_real64   !< cumulative SSDI flux (cm)
      real(real64) :: cqrot     = 0.0_real64   !< cumulative root uptake (cm)
      real(real64) :: cqbot     = 0.0_real64   !< cumulative bottom flux (cm)
      real(real64) :: cqbotdo   = 0.0_real64   !< cumulative downward bottom flux (cm)
      real(real64) :: cqbotup   = 0.0_real64   !< cumulative upward bottom flux (cm)
      real(real64) :: cinund    = 0.0_real64   !< cumulative inundation (cm)
      real(real64) :: crunon    = 0.0_real64   !< cumulative runon (cm)
      real(real64) :: crunoff   = 0.0_real64   !< cumulative runoff (cm)
      real(real64) :: crunoffCN = 0.0_real64   !< cumulative CN runoff (cm)
      real(real64) :: cqtdo     = 0.0_real64   !< cumulative total downward flux (cm)
      real(real64) :: cqtup     = 0.0_real64   !< cumulative total upward flux (cm)
      real(real64) :: cqprai    = 0.0_real64   !< cumulative precipitation (cm)
      real(real64) :: cgird     = 0.0_real64   !< cumulative gross irrigation (cm)
      real(real64) :: cnird     = 0.0_real64   !< cumulative net irrigation (cm)
   contains
      procedure :: reset => soilwater_cumulative_reset                !< zeroes all 14 fields (flzerocumu)
   end type soilwater_cumulative_t

   ! ---------------------------------------------------------------------------
   !> Top-level soil-water state record. 109 fields total (35 existing + 74 new):
   !!   12 boundary (ADR 0035)
   !!   22 crop water uptake (ADR 0036)
   !!   30 flat instantaneous — 17 per-node/layer arrays + 13 scalars/flags (ADR 0038)
   !!   soilwater_intermediate_t :: intr  — 22+8 = 30 fields in cohort
   !!   soilwater_cumulative_t   :: cumu  — 14 fields in cohort
   !!
   !! Flat instantaneous arrays are unallocated until soilwater_init (S-1.2).
   !! Cohort arrays are unallocated until soilwater_init (S-1.2).
   !! ASSOCIATE prefix sw_ recommended in heavy compute bodies.
   ! ---------------------------------------------------------------------------
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
      ! Per-node arrays (allocated by soilwater_init, C-1.2):

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

      ! ===========================================================================
      ! SOIL-WATER CORE — FLAT INSTANTANEOUS (ADR 0038, S-1.1)
      ! ===========================================================================
      ! Per-node arrays (17; allocated by soilwater_init S-1.2):

      real(real64), allocatable :: theta(:)        !< volumetric water content per node (-)
      real(real64), allocatable :: thetm1(:)       !< theta at previous time level (-)
      real(real64), allocatable :: thetar(:)       !< residual water content per node (-)
      real(real64), allocatable :: thetas(:)       !< saturated water content per node (-)
      real(real64), allocatable :: h(:)            !< pressure head per node (cm)
      real(real64), allocatable :: hm1(:)          !< h at previous time level (cm)
      real(real64), allocatable :: q(:)            !< inter-compartment flux per node (cm/d)
      real(real64), allocatable :: k(:)            !< hydraulic conductivity per node (cm/d)
      real(real64), allocatable :: kmean(:)        !< mean K at node interface (cm/d)
      real(real64), allocatable :: dimoca(:)       !< differential moisture capacity per node (1/cm)
      real(real64), allocatable :: cofgen(:,:)     !< Mualem-VG parameters (21 × numnod)
      real(real64), allocatable :: FrArMtrx(:)    !< matrix-area fraction per node (-)
      logical,      allocatable :: fluseksatexm(:) !< per-node Ksatexm flag (-)
      integer,      allocatable :: indeks(:)       !< hysteresis branch index per node (+1/-1)
      real(real64), allocatable :: evp(:)          !< per-node evaporation (cm/d) — always zero; kept for legacy parity
      real(real64), allocatable :: thetsl(:)       !< saturated water content per layer (-); size numlay

      ! Flat instantaneous scalars + integers + logicals (13):

      real(real64) :: pond         = 0.0_real64   !< surface ponding depth (cm)
      real(real64) :: pondm1       = 0.0_real64   !< ponding at previous time level (cm)
      real(real64) :: pondini      = 0.0_real64   !< ponding at cumu-period start (cm)
      real(real64) :: gwl          = 0.0_real64   !< groundwater level below surface (cm)
      real(real64) :: gwlm1        = 0.0_real64   !< gwl at previous time level (cm)
      integer      :: nodgwl       = 0            !< node directly above gwl (-)
      real(real64) :: pegwl        = 0.0_real64   !< perched groundwater level (cm)
      integer      :: bpegwl       = 0            !< node at bottom of perched gwl (-)
      integer      :: npegwl       = 0            !< node above perched gwl (-)
      real(real64) :: gwlflcpzo    = 0.0_real64   !< capillary-zone gwl (cm)
      integer      :: nodgwlflcpzo = 0            !< node index for capillary-zone gwl (-)
      real(real64) :: hatm         = 0.0_real64   !< air pressure head near surface (cm); init to −2.75e5 in soilwater_init
      real(real64) :: volact       = 0.0_real64   !< current soil-profile water storage (cm)
      real(real64) :: volm1        = 0.0_real64   !< volact at previous time level (cm)
      real(real64) :: volini       = 0.0_real64   !< storage at cumu-period start (cm)
      real(real64) :: wbalance     = 0.0_real64   !< cumulative water balance error (cm)
      real(real64) :: runon        = 0.0_real64   !< runon flux this step (cm/d)
      logical      :: fllowgwl     = .false.      !< flag: gwl is below the soil profile

      ! ===========================================================================
      ! SOIL-WATER CORE — COHORT SUB-RECORDS (ADR 0038, S-1.1)
      ! ===========================================================================

      type(soilwater_intermediate_t) :: intr  !< intermediate accumulators (flzerointr + flDayStart gates)
      type(soilwater_cumulative_t)   :: cumu  !< cumulative accumulators (flzerocumu gate)

   end type soilwater_state_t

contains

   !> Lifecycle init for soilwater typed state.
   !! Zeros/resets all scalar fields in both the boundary subset (ADR 0035)
   !! and the crop-uptake scalar subset (ADR 0036).  Also allocates all 13
   !! crop-uptake per-node / per-layer allocatable arrays and zeros them:
   !!   12 per-node arrays sized numnod:
   !!     qrot, qpotrot, qredwet, qreddry, qredsol, qredfrs  (Feddes path)
   !!     mflux, mroot, hroot, rootrho, rootphi, rmax         (JvL path)
   !!   1 lookup table mfluxtable(nlay, 801) — allocated here, zeroed here.
   !!     NOTE: the table is built (filled) by MatricFlux(1) inside
   !!     CropGrowth(1) — not here.  MatricFlux(1) depends on ksatfit,
   !!     wiltpoint, and cofgen which are not yet populated when
   !!     soilwater_init runs (they are set by SoilHydraulics(1) inside
   !!     SoilWater(1) at swap.f90:190).  Relocation of the build step is
   !!     therefore BLOCKED; see ADR 0036 §mfluxtable disposition.
   !!
   !! Also allocates and zeroes all soil-water core arrays (ADR 0038, S-1.2):
   !!   14 flat per-node arrays sized numnod: theta, thetm1, thetar, thetas,
   !!     h, hm1, dimoca, FrArMtrx, evp, indeks, fluseksatexm (logical)
   !!   3 flat per-node+1 arrays sized numnod+1: q, k, kmean (flux boundaries)
   !!   1 flat 2D parameter array: cofgen(21, numnod) — Mualem-VG params
   !!   1 flat per-layer array sized nlay: thetsl
   !!   intr cohort arrays: inqrot, inqssdi, IThetaBeg (numnod),
   !!                       inq, iqdo, iqup (numnod+1),
   !!                       qpotrot_day, qredtot_day (numnod)
   !!   Non-zero default: hatm = -2.75e5_real64 (mirrors soilhydraulics.f90:870)
   !!
   !! Takes soilwater_state_t directly (not swap_state_t) to avoid a circular
   !! dependency: soilwater_state_mod is used by swap_state_mod.
   !! Mirrors heat_init pattern but at the sub-record level.
   !! Called from swap.f90 immediately after CalcGrid(), before DoTillage(1).
   !!
   !! Design: docs/superpowers/specs/2026-05-10-state-migration-boundary-design.md D8
   !!         docs/superpowers/specs/2026-05-10-state-migration-crop-uptake-design.md D4/D6
   !!         docs/superpowers/specs/2026-05-11-state-migration-soilwater-core-design.md D4
   subroutine soilwater_init(sw, numnod, nlay)
      type(soilwater_state_t), intent(inout) :: sw
      integer, intent(in) :: numnod  !! number of soil nodes (from CalcGrid)
      integer, intent(in) :: nlay    !! number of soil layers (from CalcGrid)

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

      ! Per-node arrays — Feddes + stress path (ADR 0036, C-1.2)
      allocate(sw%qrot(numnod));    sw%qrot    = 0.0_real64
      allocate(sw%qpotrot(numnod)); sw%qpotrot = 0.0_real64
      allocate(sw%qredwet(numnod)); sw%qredwet = 0.0_real64
      allocate(sw%qreddry(numnod)); sw%qreddry = 0.0_real64
      allocate(sw%qredsol(numnod)); sw%qredsol = 0.0_real64
      allocate(sw%qredfrs(numnod)); sw%qredfrs = 0.0_real64

      ! Per-node arrays — JvL microscopic path (ADR 0036, C-1.2)
      allocate(sw%mflux(numnod));   sw%mflux   = 0.0_real64
      allocate(sw%mroot(numnod));   sw%mroot   = 0.0_real64
      allocate(sw%hroot(numnod));   sw%hroot   = 0.0_real64
      allocate(sw%rootrho(numnod)); sw%rootrho = 0.0_real64
      allocate(sw%rootphi(numnod)); sw%rootphi = 0.0_real64
      allocate(sw%rmax(numnod));    sw%rmax    = 0.0_real64

      ! Lookup table — allocated here; built by MatricFlux(1) in CropGrowth(1)
      ! (build is BLOCKED from relocation: ksatfit/wiltpoint/cofgen not yet
      ! populated when this routine runs — see ADR 0036 §mfluxtable disposition)
      allocate(sw%mfluxtable(nlay, 801)); sw%mfluxtable = 0.0_real64

      ! ===========================================================================
      ! SOIL-WATER CORE — flat per-node / per-layer arrays (ADR 0038, S-1.2)
      ! ===========================================================================

      ! Per-node arrays sized numnod (10 arrays):
      allocate(sw%theta(numnod));        sw%theta        = 0.0_real64
      allocate(sw%thetm1(numnod));       sw%thetm1       = 0.0_real64
      allocate(sw%thetar(numnod));       sw%thetar       = 0.0_real64
      allocate(sw%thetas(numnod));       sw%thetas       = 0.0_real64
      allocate(sw%h(numnod));            sw%h            = 0.0_real64
      allocate(sw%hm1(numnod));          sw%hm1          = 0.0_real64
      allocate(sw%dimoca(numnod));       sw%dimoca       = 0.0_real64
      allocate(sw%FrArMtrx(numnod));     sw%FrArMtrx     = 0.0_real64
      allocate(sw%evp(numnod));          sw%evp          = 0.0_real64
      allocate(sw%indeks(numnod));       sw%indeks       = 0

      ! Per-node logical array sized numnod:
      allocate(sw%fluseksatexm(numnod)); sw%fluseksatexm = .false.

      ! Per-node arrays sized numnod+1 (flux arrays — one value per node boundary):
      ! Legacy: q(macp+1), k(macp+1), kmean(macp+1)
      allocate(sw%q(numnod+1));          sw%q            = 0.0_real64
      allocate(sw%k(numnod+1));          sw%k            = 0.0_real64
      allocate(sw%kmean(numnod+1));      sw%kmean        = 0.0_real64

      ! 2D parameter array: cofgen(21, numnod) — legacy cofgen(21, macp)
      allocate(sw%cofgen(21, numnod));   sw%cofgen       = 0.0_real64

      ! Per-layer array sized nlay (thetsl per soil layer, not per node):
      ! Legacy: thetsl(maho) where maho = max number of soil layers
      allocate(sw%thetsl(nlay));         sw%thetsl       = 0.0_real64

      ! ===========================================================================
      ! SOIL-WATER CORE — intermediate cohort arrays (ADR 0038, S-1.2)
      ! ===========================================================================

      ! Non-per-day per-node arrays sized numnod:
      allocate(sw%intr%inqrot(numnod));    sw%intr%inqrot    = 0.0_real64
      allocate(sw%intr%inqssdi(numnod));   sw%intr%inqssdi   = 0.0_real64
      allocate(sw%intr%IThetaBeg(numnod)); sw%intr%IThetaBeg = 0.0_real64

      ! Non-per-day per-node flux arrays sized numnod+1:
      ! Legacy: inq(macp+1), iqdo(macp+1), iqup(macp+1)
      allocate(sw%intr%inq(numnod+1));     sw%intr%inq       = 0.0_real64
      allocate(sw%intr%iqdo(numnod+1));    sw%intr%iqdo      = 0.0_real64
      allocate(sw%intr%iqup(numnod+1));    sw%intr%iqup      = 0.0_real64

      ! Per-day per-node arrays sized numnod:
      ! Legacy: qpotrot_day(macp), qredtot_day(macp)
      allocate(sw%intr%qpotrot_day(numnod)); sw%intr%qpotrot_day = 0.0_real64
      allocate(sw%intr%qredtot_day(numnod)); sw%intr%qredtot_day = 0.0_real64

      ! ===========================================================================
      ! NON-ZERO DEFAULTS (ADR 0038, S-1.2)
      ! ===========================================================================

      ! hatm: air pressure head near soil surface; legacy soilhydraulics.f90:870
      ! sets hatm = -2.75d+05 at SoilWater(1) task 1 init.  We mirror that here
      ! so sw%hatm is consistent from the moment soilwater_init returns.
      sw%hatm = -2.75e5_real64

   end subroutine soilwater_init

   ! ---------------------------------------------------------------------------
   !> Zero ALL fields in the intermediate cohort — called under flzerointr gate.
   !! Zeroes all 6 non-per-day allocatable arrays (if allocated), all 20
   !! non-per-day scalars, and all 8 per-day fields (6 scalars + 2 arrays).
   !! Pattern: atmosphere_intermediate_reset (ADR 0037) extended for arrays.
   ! ---------------------------------------------------------------------------
   subroutine soilwater_intermediate_reset(self)
      class(soilwater_intermediate_t), intent(inout) :: self

      ! Non-per-day allocatable arrays (allocated() guard — not yet allocated pre-S-1.2)
      if (allocated(self%inq))       self%inq       = 0.0_real64
      if (allocated(self%inqrot))    self%inqrot    = 0.0_real64
      if (allocated(self%inqssdi))   self%inqssdi   = 0.0_real64
      if (allocated(self%iqdo))      self%iqdo      = 0.0_real64
      if (allocated(self%iqup))      self%iqup      = 0.0_real64
      if (allocated(self%IThetaBeg)) self%IThetaBeg = 0.0_real64

      ! Non-per-day scalars
      self%iqrot    = 0.0_real64
      self%iqssdi   = 0.0_real64
      self%iqredwet = 0.0_real64
      self%iqreddry = 0.0_real64
      self%iqredsol = 0.0_real64
      self%iqredfrs = 0.0_real64
      self%ies0     = 0.0_real64
      self%iet0     = 0.0_real64
      self%iew0     = 0.0_real64
      self%iintc    = 0.0_real64
      self%iruno    = 0.0_real64
      self%irunoCN  = 0.0_real64
      self%irunon   = 0.0_real64
      self%iqbot    = 0.0_real64
      self%iqtdo    = 0.0_real64
      self%iqtup    = 0.0_real64
      self%IPondBeg = 0.0_real64
      self%iprec    = 0.0_real64
      self%igird    = 0.0_real64
      self%inird    = 0.0_real64

      ! Per-day fields — also zeroed by full reset (flzerointr subsumes flDayStart)
      call soilwater_per_day_reset(self)

   end subroutine soilwater_intermediate_reset

   ! ---------------------------------------------------------------------------
   !> Zero only the 8 per-day fields — called under flDayStart gate.
   !! Zeroes 6 per-day scalars + 2 per-day arrays (if allocated).
   !! Non-per-day intermediate fields are UNTOUCHED (orthogonality invariant).
   ! ---------------------------------------------------------------------------
   subroutine soilwater_per_day_reset(self)
      class(soilwater_intermediate_t), intent(inout) :: self

      ! Per-day scalars
      self%tra          = 0.0_real64
      self%iqredwet_day = 0.0_real64
      self%iqreddry_day = 0.0_real64
      self%iqredsol_day = 0.0_real64
      self%iqredfrs_day = 0.0_real64
      self%iptra_day    = 0.0_real64

      ! Per-day arrays (allocated() guard)
      if (allocated(self%qpotrot_day)) self%qpotrot_day = 0.0_real64
      if (allocated(self%qredtot_day)) self%qredtot_day = 0.0_real64

   end subroutine soilwater_per_day_reset

   ! ---------------------------------------------------------------------------
   !> Zero all 14 cumulative cohort fields — called under flzerocumu gate.
   !! No allocatable arrays in this type; no allocated() guards needed.
   !! Pattern: atmosphere_cumulative_reset (ADR 0037).
   ! ---------------------------------------------------------------------------
   subroutine soilwater_cumulative_reset(self)
      class(soilwater_cumulative_t), intent(inout) :: self
      self%cqssdi    = 0.0_real64
      self%cqrot     = 0.0_real64
      self%cqbot     = 0.0_real64
      self%cqbotdo   = 0.0_real64
      self%cqbotup   = 0.0_real64
      self%cinund    = 0.0_real64
      self%crunon    = 0.0_real64
      self%crunoff   = 0.0_real64
      self%crunoffCN = 0.0_real64
      self%cqtdo     = 0.0_real64
      self%cqtup     = 0.0_real64
      self%cqprai    = 0.0_real64
      self%cgird     = 0.0_real64
      self%cnird     = 0.0_real64
   end subroutine soilwater_cumulative_reset

end module soilwater_state_mod
