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
!!   Intermediate accumulators (flattened onto parent type):
!!     22 non-per-day fields (6 arrays + 16 scalars) + 8 per-day fields
!!     reset_intermediate() = flzerointr; reset_intermediate_per_day() = flDayStart
!!   Cumulative accumulators (flattened onto parent type):
!!     14 scalar fields; reset_cumulative() = flzerocumu
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
   use iso_c_binding, only: c_double
   use hydraulic_params_mod, only: vanGenuchten_params_t
   use swap_array_dimensions, only: MADAY, MABBC
   implicit none
   private
   public :: soilwater_state_t, soilwater_init

   ! ---------------------------------------------------------------------------
   !> Top-level soil-water state record. 109 fields total (35 existing + 74 new):
   !!   12 boundary (ADR 0035)
   !!   22 crop water uptake (ADR 0036)
   !!   30 flat instantaneous — 17 per-node/layer arrays + 13 scalars/flags (ADR 0038)
   !!   intermediate (flat) — 22+8 = 30 fields, reset_intermediate / reset_intermediate_per_day
   !!   cumulative (flat)   — 14 fields, reset_cumulative
   !!
   !! Flat instantaneous arrays are unallocated until soilwater_init (S-1.2).
   !! Intermediate per-node arrays are unallocated until soilwater_init (S-1.2).
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
      real(real64), allocatable :: qssdi(:)        !< per-node SSDI source flux (cm/d) — written by irrigation, read by waterbalance/headcalc
      real(real64), allocatable :: k(:)            !< hydraulic conductivity per node (cm/d)
      real(real64), allocatable :: kmean(:)        !< mean K at node interface (cm/d)
      real(real64), allocatable :: dimoca(:)       !< differential moisture capacity per node (1/cm)
      type(vanGenuchten_params_t), allocatable :: vg_params(:)   !< [SS-GR-UTILS] typed VG parameters, one per node

      ! [SS-GR-UTILS] Soil hydraulic property metadata (migrated from variables.f90).
      ! Populated by SoilHydraulics(1) / config_to_variables as transitional
      ! dual-writes; bare globals retire in Arc 9.
      integer                       :: swsophy    = 0   !< soil-hydraulic-property switch (0=MvG, 1=table)
      integer,         allocatable  :: numtab(:)        !< per-node table entry count
      real(real64),    allocatable  :: sptab(:,:,:)     !< soil property table (7 × numnod × matab)
      integer,         allocatable  :: ientrytab(:,:)   !< entry indices (numnod × 0:matabentries)
      integer,         allocatable  :: iHWCKmodel(:)    !< per-layer hydraulic-K model selector
      integer,         allocatable  :: layer(:)         !< per-node soil-layer index
      integer                       :: swfrost    = 0   !< frost-reduction switch (0=no, 1=yes)
      logical,         allocatable  :: BiModal(:)       !< per-layer bi-modal flag
      logical,         allocatable  :: NoVap(:)         !< per-layer no-vapor flag
      real(real64), allocatable :: FrArMtrx(:)    !< matrix-area fraction per node (-)
      logical,      allocatable :: fluseksatexm(:) !< per-node Ksatexm flag (-)
      integer,      allocatable :: indeks(:)       !< hysteresis branch index per node (+1/-1)
      real(real64), allocatable :: evp(:)          !< per-node evaporation (cm/d) — always zero; kept for legacy parity
      real(real64), allocatable :: thetsl(:)       !< saturated water content per layer (-); size numlay

      ! [SS-GR-BH A4] layer-flat fields — maho-sized, one per soil layer
      real(real64), allocatable :: ksatexm(:)    !! layer Ksat (examined extension)
      real(real64), allocatable :: ksatfit(:)    !! layer fitted Ksat
      real(real64), allocatable :: cofani(:)     !! layer anisotropy coefficient
      logical                   :: flksatexm = .false.  !! global flag: Ksatexm present in input
      real(real64), allocatable :: orgmat(:)     !! layer gravimetric organic matter
      real(real64), allocatable :: bdens(:)      !! layer dry bulk density (g/cm3)
      real(real64), allocatable :: psand(:)      !! layer sand fraction
      real(real64), allocatable :: psilt(:)      !! layer silt fraction
      real(real64), allocatable :: pclay(:)      !! layer clay fraction

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
      real(real64) :: qssdisum     = 0.0_real64   !< column-sum SSDI source flux (cm/d); set by irrigation alongside qssdi
      logical      :: fllowgwl     = .false.      !< flag: gwl is below the soil profile

      ! Runon feature (config-rooted switch + day-indexed time series).
      ! flrunon is set from config%soil%swrunon by config_to_variables;
      ! runonarr is currently dormant (no TOML wiring) — boundtop reads
      ! runonarr(daycum+1) when flrunon is true. Allocated to MADAY at init.
      logical                   :: flrunon  = .false.
      real(real64), allocatable :: runonarr(:)

      ! Initial soil-water condition switch (snapshotted from config%soil%swinco).
      ! 1 = pressure heads; 2 = hydrostatic equilibrium; 3 = warm-restart from CSV.
      integer :: swinco = 1

      ! Cauchy bottom-boundary (swbotb=3) vertical-resistance switch.
      ! Dormant — no TOML writer; always 0 in the TOML pipeline.
      ! 0 = add the (modelled-profile) vertical resistance to rimlay; 1 = use rimlay alone.
      integer :: swbotb3resvert = 0

      ! Bottom-boundary CSV-driven tables (date/value interleaved pairs,
      ! sized 2*MABBC to match legacy fixed-size globals). Populated by
      ! config_to_variables from the bottom_boundary CSV inputs:
      !   gwltab — prescribed groundwater level, swbotb=1
      !   qbotab — bottom flux table; multi-use: swbotb=2 (sw2=2),
      !            swbotb=3 (sw4=1 extra flux), swbotb=4 (swqhbot=2 q(h))
      !   haqtab — deep-aquifer head, swbotb=3 with sw3=2
      !   hbotab — bottom pressure head, swbotb=5
      real(real64), allocatable :: gwltab(:)
      real(real64), allocatable :: qbotab(:)
      real(real64), allocatable :: haqtab(:)
      real(real64), allocatable :: hbotab(:)

      ! [SS-GR-BH A5] runtime scalars formerly bare globals (boundtop/PONDRUNOFF/boundbottom)
      real(real64) :: q0           = 0.0_real64   !! surface flux (precip + runon - reva) [cm/d]
      real(real64) :: k1max        = 0.0_real64   !! max conductivity at z=0 [cm/d]
      real(real64) :: H0max        = 0.0_real64   !! max ponding pre-runoff [cm]
      integer      :: swbotb_runtime = 0          !! runtime-overridable bottom-boundary switch

      ! [GR-SOIL 2026-05-24] Richards non-convergence warning state (formerly bare globals
      ! that were "previously SAVE variables" — pinned at module scope to persist across
      ! headcalc calls). Reset to defaults each flDayStart at the top of headcalc.
      logical :: flwarn_hc = .true.   !! emit Richards non-convergence warning this day
      integer :: iwarn_hc  = 0        !! count of warnings emitted today (cap at 4)

      ! [GR-SOIL 2026-05-24] Richards iteration counter + statistics (formerly bare globals).
      ! `numbit` is the live iteration index inside headcalc's Newton-Raphson loop, read
      ! by timecontrol_advance for dt adjustment. `Itnumb` is the (100×2) histogram of
      ! iteration counts / back-tracking cycles, accumulated in headcalc and dumped at
      ! end-of-run by itertime_close.
      integer              :: numbit = 0
      integer, allocatable :: Itnumb(:,:)

      ! ===========================================================================
      ! INTERMEDIATE accumulators (reset_intermediate / gate: flzerointr)
      !   subsumes the per-day subset, which has its own reset under flDayStart.
      ! ===========================================================================

      ! Non-per-day per-node arrays (allocated by soilwater_init)
      real(real64), allocatable :: inq(:)        !< intra-period inter-comp flux (cm)
      real(real64), allocatable :: inqrot(:)     !< intra-period root uptake (cm)
      real(real64), allocatable :: inqssdi(:)    !< intra-period SSDI flux (cm)
      real(real64), allocatable :: iqdo(:)       !< intra-period downward flux (cm)
      real(real64), allocatable :: iqup(:)       !< intra-period upward flux (cm)
      real(real64), allocatable :: IThetaBeg(:)  !< theta at start of intermediate period (-)

      ! Non-per-day scalars
      real(real64) :: iqrot     = 0.0_real64
      real(real64) :: iqssdi    = 0.0_real64
      real(real64) :: iqredwet  = 0.0_real64
      real(real64) :: iqreddry  = 0.0_real64
      real(real64) :: iqredsol  = 0.0_real64
      real(real64) :: iqredfrs  = 0.0_real64
      real(real64) :: ies0      = 0.0_real64
      real(real64) :: iet0      = 0.0_real64
      real(real64) :: iew0      = 0.0_real64
      real(real64) :: iintc     = 0.0_real64
      real(real64) :: iruno     = 0.0_real64
      real(real64) :: irunoCN   = 0.0_real64
      real(real64) :: irunon    = 0.0_real64
      real(real64) :: iqbot     = 0.0_real64
      real(real64) :: iqtdo     = 0.0_real64
      real(real64) :: iqtup     = 0.0_real64
      real(real64) :: IPondBeg  = 0.0_real64
      real(real64) :: iprec     = 0.0_real64
      real(real64) :: igird     = 0.0_real64
      real(real64) :: inird     = 0.0_real64

      ! ===========================================================================
      ! PER-DAY subset (reset_intermediate_per_day / gate: flDayStart)
      !   Also zeroed as part of reset_intermediate (flzerointr subsumes flDayStart).
      ! ===========================================================================
      real(real64) :: tra           = 0.0_real64
      real(real64) :: iqredwet_day  = 0.0_real64
      real(real64) :: iqreddry_day  = 0.0_real64
      real(real64) :: iqredsol_day  = 0.0_real64
      real(real64) :: iqredfrs_day  = 0.0_real64
      real(real64) :: iptra_day     = 0.0_real64
      real(real64), allocatable :: qpotrot_day(:)
      real(real64), allocatable :: qredtot_day(:)

      ! ===========================================================================
      ! CUMULATIVE accumulators (reset_cumulative / gate: flzerocumu)
      ! ===========================================================================
      real(real64) :: cqssdi    = 0.0_real64
      real(real64) :: cqrot     = 0.0_real64
      real(real64) :: cqbot     = 0.0_real64
      real(real64) :: cqbotdo   = 0.0_real64
      real(real64) :: cqbotup   = 0.0_real64
      real(real64) :: cinund    = 0.0_real64
      real(real64) :: crunon    = 0.0_real64
      real(real64) :: crunoff   = 0.0_real64
      real(real64) :: crunoffCN = 0.0_real64
      real(real64) :: cqtdo     = 0.0_real64
      real(real64) :: cqtup     = 0.0_real64
      real(real64) :: cqprai    = 0.0_real64
      real(real64) :: cgird     = 0.0_real64
      real(real64) :: cnird     = 0.0_real64

      ! ===========================================================================
      ! [SS-BMI2] SOILWATER OUTPUT STREAM (Task 9)
      ! ===========================================================================
      !> Row buffer filled by build_soilwater_output_row; sized to the number
      !! of active columns in the user-defined csv output (swcsv=1 path).
      real(c_double),    allocatable :: output_row(:)
      character(len=32), allocatable :: output_columns(:)
      integer                        :: output_n_cols = 0

   contains
      procedure :: reset_intermediate         => soilwater_reset_intermediate
      procedure :: reset_intermediate_per_day => soilwater_reset_intermediate_per_day
      procedure :: reset_cumulative           => soilwater_reset_cumulative

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
   !!   intermediate per-node arrays: inqrot, inqssdi, IThetaBeg (numnod),
   !!                                 inq, iqdo, iqup (numnod+1),
   !!                                 qpotrot_day, qredtot_day (numnod)
   !!   Non-zero default: hatm = -2.75e5_real64 (mirrors soilhydraulics.f90:870)
   !!
   !! Takes soilwater_state_t directly (not swap_state_t) to avoid a circular
   !! dependency: soilwater_state_mod is used by swap_state_mod.
   !! Mirrors the heat_init pattern.
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
      ! [MACRO-RETIRE 2026-05-12] FrArMtrx defaults to 1.0 (whole-matrix);
      ! macropore retirement removed the only writer (MACROGEOM). ADR 0040.
      allocate(sw%FrArMtrx(numnod));     sw%FrArMtrx     = 1.0_real64
      allocate(sw%evp(numnod));          sw%evp          = 0.0_real64
      allocate(sw%indeks(numnod));       sw%indeks       = 0

      ! Per-node logical array sized numnod:
      allocate(sw%fluseksatexm(numnod)); sw%fluseksatexm = .false.

      ! Per-node arrays sized numnod+1 (flux arrays — one value per node boundary):
      ! Legacy: q(macp+1), k(macp+1), kmean(macp+1)
      allocate(sw%q(numnod+1));          sw%q            = 0.0_real64
      allocate(sw%qssdi(numnod));        sw%qssdi        = 0.0_real64  ! [GR-SOIL 2026-05-24] migrated from variables.f90
      if (.not. allocated(sw%Itnumb)) then
         allocate(sw%Itnumb(100, 2));    sw%Itnumb       = 0            ! [GR-SOIL 2026-05-24] Richards iter stats
      end if
      allocate(sw%k(numnod+1));          sw%k            = 0.0_real64
      allocate(sw%kmean(numnod+1));      sw%kmean        = 0.0_real64

      allocate(sw%vg_params(numnod))    ! [SS-GR-UTILS] components default-init from type

      ! [SS-GR-UTILS] Soil hydraulic property metadata (Task 4)
      ! Shapes mirror the legacy fixed-size globals in variables.f90:
      !   numtab(macp), sptab(7,macp,matab), ientrytab(macp,0:matabentries)
      !   iHWCKmodel(maho), layer(macp), BiModal(maho), NoVap(maho)
      ! sptab/numtab/ientrytab only allocated when swsophy=1 (tabulated path).
      ! matab=1000, matabentries=50005 — static upper bounds, runtime numnod/nlay.
      ! iHWCKmodel/layer/BiModal/NoVap allocated unconditionally (small, per-layer/node).
      allocate(sw%iHWCKmodel(nlay));             sw%iHWCKmodel = 0
      allocate(sw%layer(numnod));                sw%layer      = 0
      allocate(sw%BiModal(nlay));                sw%BiModal    = .false.
      allocate(sw%NoVap(nlay));                  sw%NoVap      = .false.
      if (sw%swsophy == 1) then
         allocate(sw%numtab(numnod));            sw%numtab     = 0
         allocate(sw%sptab(7, numnod, 1000));    sw%sptab      = 0.0_real64
         allocate(sw%ientrytab(numnod, 0:50005));sw%ientrytab  = 0
      end if

      ! Per-layer array sized nlay (thetsl per soil layer, not per node):
      ! Legacy: thetsl(maho) where maho = max number of soil layers
      allocate(sw%thetsl(nlay));         sw%thetsl       = 0.0_real64

      ! [SS-GR-BH A4] Layer-flat fields (maho-sized, one per soil layer):
      allocate(sw%ksatexm(nlay));        sw%ksatexm      = 0.0_real64
      allocate(sw%ksatfit(nlay));        sw%ksatfit      = 0.0_real64
      allocate(sw%cofani(nlay));         sw%cofani       = 0.0_real64
      sw%flksatexm = .false.
      allocate(sw%orgmat(nlay));         sw%orgmat       = 0.0_real64
      ! bdens: Pattern 9 guarded alloc (config_to_variables may have done it first).
      if (.not. allocated(sw%bdens)) then
         allocate(sw%bdens(nlay));       sw%bdens        = 0.0_real64
      end if
      allocate(sw%psand(nlay));          sw%psand        = 0.0_real64
      allocate(sw%psilt(nlay));          sw%psilt        = 0.0_real64
      allocate(sw%pclay(nlay));          sw%pclay        = 0.0_real64

      ! [SS-GR-BH A5] Runtime scalars (seeded from config in Task 7)
      sw%q0              = 0.0_real64
      sw%k1max           = 0.0_real64
      sw%H0max           = 0.0_real64
      sw%swbotb_runtime  = 0

      ! ===========================================================================
      ! SOIL-WATER CORE — intermediate-period arrays (ADR 0038, S-1.2)
      ! ===========================================================================

      ! Non-per-day per-node arrays sized numnod:
      allocate(sw%inqrot(numnod));    sw%inqrot    = 0.0_real64
      allocate(sw%inqssdi(numnod));   sw%inqssdi   = 0.0_real64
      allocate(sw%IThetaBeg(numnod)); sw%IThetaBeg = 0.0_real64

      ! Non-per-day per-node flux arrays sized numnod+1:
      ! Legacy: inq(macp+1), iqdo(macp+1), iqup(macp+1)
      allocate(sw%inq(numnod+1));     sw%inq       = 0.0_real64
      allocate(sw%iqdo(numnod+1));    sw%iqdo      = 0.0_real64
      allocate(sw%iqup(numnod+1));    sw%iqup      = 0.0_real64

      ! Per-day per-node arrays sized numnod:
      ! Legacy: qpotrot_day(macp), qredtot_day(macp)
      allocate(sw%qpotrot_day(numnod)); sw%qpotrot_day = 0.0_real64
      allocate(sw%qredtot_day(numnod)); sw%qredtot_day = 0.0_real64

      ! ===========================================================================
      ! NON-ZERO DEFAULTS (ADR 0038, S-1.2)
      ! ===========================================================================

      ! hatm: air pressure head near soil surface; legacy soilhydraulics.f90:870
      ! sets hatm = -2.75d+05 at SoilWater(1) task 1 init.  We mirror that here
      ! so sw%hatm is consistent from the moment soilwater_init returns.
      sw%hatm = -2.75e5_real64

      ! Runon time-series (dormant — no TOML writer yet). Sized to MADAY to
      ! match the legacy fixed-size global runonarr(maday).
      allocate(sw%runonarr(MADAY)); sw%runonarr = 0.0_real64
      sw%flrunon = .false.

      ! Bottom-boundary CSV tables: sized 2*MABBC to match legacy globals.
      ! Allocation is guarded because config_to_variables (which runs
      ! BEFORE soilwater_init) may have already allocated and populated
      ! these arrays from the swbotb CSV inputs.
      if (.not. allocated(sw%gwltab)) then
         allocate(sw%gwltab(2*MABBC)); sw%gwltab = 0.0_real64
      end if
      if (.not. allocated(sw%qbotab)) then
         allocate(sw%qbotab(2*MABBC)); sw%qbotab = 0.0_real64
      end if
      if (.not. allocated(sw%haqtab)) then
         allocate(sw%haqtab(2*MABBC)); sw%haqtab = 0.0_real64
      end if
      if (.not. allocated(sw%hbotab)) then
         allocate(sw%hbotab(2*MABBC)); sw%hbotab = 0.0_real64
      end if

   end subroutine soilwater_init

   !> Zero ALL intermediate-period fields — flzerointr gate.
   !! Includes the per-day subset (flzerointr subsumes flDayStart).
   subroutine soilwater_reset_intermediate(self)
      class(soilwater_state_t), intent(inout) :: self

      ! Non-per-day allocatable arrays
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
      call soilwater_reset_intermediate_per_day(self)

   end subroutine soilwater_reset_intermediate

   !> Zero only the per-day subset — flDayStart gate.
   !! Non-per-day intermediate fields are UNTOUCHED.
   subroutine soilwater_reset_intermediate_per_day(self)
      class(soilwater_state_t), intent(inout) :: self
      self%tra          = 0.0_real64
      self%iqredwet_day = 0.0_real64
      self%iqreddry_day = 0.0_real64
      self%iqredsol_day = 0.0_real64
      self%iqredfrs_day = 0.0_real64
      self%iptra_day    = 0.0_real64
      if (allocated(self%qpotrot_day)) self%qpotrot_day = 0.0_real64
      if (allocated(self%qredtot_day)) self%qredtot_day = 0.0_real64
   end subroutine soilwater_reset_intermediate_per_day

   !> Zero the 14 cumulative fields — flzerocumu gate.
   subroutine soilwater_reset_cumulative(self)
      class(soilwater_state_t), intent(inout) :: self
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
   end subroutine soilwater_reset_cumulative

end module soilwater_state_mod
