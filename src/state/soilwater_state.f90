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
!! See ADR 0035 (boundary), ADR 0036 (crop-uptake), ADR 0038 (soilwater-core).

module soilwater_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   use hydraulic_params_mod, only: vanGenuchten_params_t
   use swap_array_dimensions, only: MADAY, MABBC
   use soil_init_csv_mod,     only: h_profile_table_t
   implicit none
   private
   public :: soilwater_state_t

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
      logical      :: flcoupled_gwl = .false.  !! when .true., gwlinp is injected externally (MODFLOW coupling)
      real(real64) :: gwl_injected  = 0.0_real64 !! externally injected groundwater level, swbotb=1 coupled (cm)
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
      type(vanGenuchten_params_t), allocatable :: vg_params(:)        !< [SS-GR-UTILS] typed VG parameters, one per node
      !> [GR-CROP 2026-05-25] mutable per-layer VG parameter store.
      !! Populated by SoilHydraulics(1) from state%cfg%soil%hydraulics;
      !! mutated by tillage events (tillage.f90 Change_MvGpars); per-node
      !! vg_params(:) is then rebuilt from this layer-keyed store after
      !! each event. Replaces the legacy paramvg(21, maho) global.
      type(vanGenuchten_params_t), allocatable :: vg_params_layer(:)  !< per-layer VG parameter store (one per soil layer)

      ! [SS-GR-UTILS] Soil hydraulic property metadata (migrated from variables.f90).
      ! Populated by SoilHydraulics(1) / config_to_variables as transitional
      ! dual-writes; bare globals retire in Arc 9.
      integer                       :: swsophy    = 0   !< soil-hydraulic-property switch (0=MvG, 1=table — 1 is dormant)
      ! [GR-SOIL 2026-05-24] numtab/sptab/ientrytab retired — swsophy=1 (tabulated) path
      !   dormant; see src/soil/dormant/sptabulated.f90 for reactivation prerequisites.
      integer,         allocatable  :: iHWCKmodel(:)    !< per-layer hydraulic-K model selector
      integer,         allocatable  :: layer(:)         !< per-node soil-layer index
      integer                       :: swfrost    = 0   !< frost-reduction switch (0=no, 1=yes)
      ! [state%cfg-retirement cluster 6] config-snapshot switches for tillage compatibility checks
      integer                       :: swtill     = 0   !< tillage-event simulation switch (0=off, 1=on); snapshotted from config%soil%swtill
      integer                       :: swhyst     = 0   !< hysteresis switch (0=off, 1=on, 2=...); snapshotted from config%soil%swhyst
      integer                       :: swdiscrvert = 0  !< vertical re-discretization switch (0=no, 1=yes); snapshotted from config%soil%discretization%swdiscrvert
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
      integer      :: bpegwl       = 0            !< node at bottom of perched gwl (-)
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

      ! [W3 fix 2026-05-28] Typed initial pressure-head profile table.
      ! Populated by soilwater_state_init from the h_file CSV ONLY when swinco == 3;
      ! is_loaded == .false. for all other swinco values.  Replaces the former
      ! mutation of config%soil%initial%z_init.  Consumed by:
      !   swap_mod%swap_init_body  — copies h values into state%soilwater%h(:) (W4 fix).
      !   soilhydraulics.f90      — uses rows(i)%z for the afgen depth axis (swinco=3
      !                             consistency gate); the swinco=1 consumer of this
      !                             table is currently unreachable (dead code, see
      !                             W3/W4 comment there).
      type(h_profile_table_t) :: h_init

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

      ! Hydraulic-conductivity averaging method, snapshotted from
      ! config%simulation%numerical%swkmean at init. Default 1 mirrors the
      ! simulation_numerical_t default. Consumed by hcomean() in boundtop
      ! (cluster 4) and soilhydraulics (cluster 7).
      integer :: swkmean = 1

      ! [state%cfg-retirement cluster 7] Numerical solver control snapshots.
      ! Snapshotted from config%simulation%numerical at soilwater_state_init.
      integer      :: swkimpl                   = 0             !! K-averaging scheme: 0=explicit, 1=implicit
      integer      :: MaxBackTr                 = 3             !! max back-tracking steps per Newton iteration
      real(real64) :: gwlconv                   = 100.0_real64  !! GWL change convergence criterion (cm)
      real(real64) :: critdevh1cp               = 0.01_real64   !! relative head convergence criterion (-)
      real(real64) :: critdevh2cp               = 0.1_real64    !! absolute head convergence criterion (cm)
      real(real64) :: critdevponddt             = 1.0e-4_real64 !! pond water balance convergence criterion (cm)
      logical      :: swcaprise                 = .false.       !! cap capillary rise into root zone
      logical      :: dump_convergence_diagnostics = .false.    !! emit Richards convergence diagnostics

      ! [state%cfg-retirement cluster 7] Soil config scalar snapshots.
      real(real64) :: gwli  = 0.0_real64  !! initial groundwater level (cm); snapshotted from config%soil%gwli
      real(real64) :: tau   = 0.0_real64  !! hysteresis scanning-curve coefficient; snapshotted from config%soil%tau

      ! [state%cfg-retirement cluster 7] Bottom-boundary config snapshots.
      ! Snapshotted from config%bottom_boundary at soilwater_state_init.
      real(real64) :: hplate      = 0.0_real64  !! lysimeter plate pressure head (cm)
      real(real64) :: rimlay      = 0.0_real64  !! Cauchy bottom-boundary vertical resistance (d)
      integer      :: swbotb3impl = 0           !! swbotb=3 implementation flag (0=explicit, 1=implicit)
      integer      :: sw4         = 0           !! swbotb=3 extra groundwater flux switch (0=no, 1=yes)

      ! [state%cfg-retirement cluster 7] Per-layer wetting-curve alpha for hysteresis.
      ! Snapshotted from config%soil%hydraulics%alfaw(lay) at soilwater_state_init.
      ! Consumed by hysteresis() to determine the wetting-branch alpha per layer.
      real(real64), allocatable :: alfaw_layer(:)  !! per-layer wetting alpha for hysteresis (1/cm)

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

   contains
      procedure :: init                       => soilwater_state_init
      procedure :: reset_intermediate         => soilwater_reset_intermediate
      procedure :: reset_intermediate_per_day => soilwater_reset_intermediate_per_day
      procedure :: reset_cumulative           => soilwater_reset_cumulative

   end type soilwater_state_t

contains

   !> [GR-SEED 2026-05-25 Task 8] Type-bound lifecycle init for soilwater state.
   !! Promoted from free subroutine soilwater_init and absorbs:
   !!   [Piece A] Legacy soilwater_init body — allocation + zero-fill of all
   !!     runtime arrays (crop-uptake, soil-water core, intermediate accumulators).
   !!   [Piece B] Soil scalar seeding from config_soil:
   !!     swsophy, swinco, flrunon (swrunon→bool), bdens (from soil.bdens or
   !!     soil.hydraulics.bdens, guarded alloc), swfrost;
   !!     swinco=3 warm-restart h_file CSV → self%h_init (W3 fix 2026-05-28:
   !!     config_soil is now intent(in); config_soil%initial%z_init no longer written).
   !!   [Piece C] Bottom-boundary case dispatch (swbotb=1..5 CSV pre-loads)
   !!     via private seed_bottom_boundary(self, config_bb, pathwork) helper.
   !!
   !! Cross-subsystem writes retained in adapter (config_to_variables):
   !!   state%atmosphere%atmin7  (swinco=3 block)
   !!   state%solute%nconc/cml_init/zc_init  (optional Cml CSV, swsolu=1)
   !!   state%mesh%numlay  (isoillay derivation)
   !!   state%surfacewater%pondmx/rsro/rsroexp  (soil section)
   !!
   !! Signature: (self, config_soil, config_drain, config_heat, config_bb, config_simulation, numnod, numlay, pathwork)
   !!   config_soil        — intent(in):  read-only after W3 fix (swinco=3 h_file now → self%h_init).
   !!   config_drain       — intent(in):  drainage config (cofani source 1).
   !!   config_heat        — intent(in):  heat config (porg/psand/psilt/pclay sources).
   !!   config_bb          — intent(in):  bottom_boundary switches + file names.
   !!   config_simulation  — intent(in):  numerical solver config (swkmean snapshot).
   !!   numnod       — number of soil nodes (from CalcGrid).
   !!   numlay       — number of soil layers (from CalcGrid).
   !!   pathwork     — intent(in):  working directory prefix for CSV paths.
   !!
   !! Called from swap_mod immediately after CalcGrid(), before DoTillage(1).
   subroutine soilwater_state_init(self, config_soil, config_drain, config_heat, config_bb, config_simulation, numnod, numlay, pathwork)
      use soil_config_mod,            only: soil_config_t
      use drainage_config_mod,        only: drainage_config_t
      use heat_config_mod,            only: heat_config_t
      use bottom_boundary_config_mod, only: bottom_boundary_config_t
      use simulation_config_mod,      only: simulation_config_t
      use swap_array_dimensions,      only: maho
      class(soilwater_state_t),       intent(inout) :: self
      type(soil_config_t),            intent(in)    :: config_soil  ! intent(in) after W3 fix: no longer mutated
      type(drainage_config_t),        intent(in)    :: config_drain
      type(heat_config_t),            intent(in)    :: config_heat
      type(bottom_boundary_config_t), intent(in)    :: config_bb
      type(simulation_config_t),      intent(in)    :: config_simulation
      integer,                        intent(in)    :: numnod
      integer,                        intent(in)    :: numlay
      character(len=*),               intent(in)    :: pathwork
      integer :: i

      ! ---- [Piece A] Legacy soilwater_init body — allocation + zero-fill ----

      ! Top-boundary fields
      self%qtop      = 0.0_real64
      self%reva      = 0.0_real64
      self%hsurf     = 0.0_real64
      self%runots    = 0.0_real64
      self%QMpLatSs  = 0.0_real64
      self%ftoph     = .false.
      self%FlRunoff  = .false.

      ! Bottom-boundary fields
      self%qbot           = 0.0_real64
      self%qbot_nonfrozen = 0.0_real64
      self%hbot           = 0.0_real64
      self%gwlinp         = 0.0_real64
      self%deepgw         = 0.0_real64

      ! Crop-uptake scalars (ADR 0036 — C-1.1)
      self%qrosum      = 0.0_real64
      self%qredwetsum  = 0.0_real64
      self%qreddrysum  = 0.0_real64
      self%qredsolsum  = 0.0_real64
      self%qredfrssum  = 0.0_real64
      self%flWrtNonox  = .false.
      self%Tactual     = 0.0_real64
      self%alpJvLier   = 0.0_real64
      self%hleaf       = 0.0_real64
      self%Hxylem      = 0.0_real64

      ! Per-node arrays — Feddes + stress path (ADR 0036, C-1.2)
      allocate(self%qrot(numnod));    self%qrot    = 0.0_real64
      allocate(self%qpotrot(numnod)); self%qpotrot = 0.0_real64
      allocate(self%qredwet(numnod)); self%qredwet = 0.0_real64
      allocate(self%qreddry(numnod)); self%qreddry = 0.0_real64
      allocate(self%qredsol(numnod)); self%qredsol = 0.0_real64
      allocate(self%qredfrs(numnod)); self%qredfrs = 0.0_real64

      ! Per-node arrays — JvL microscopic path (ADR 0036, C-1.2)
      allocate(self%mflux(numnod));   self%mflux   = 0.0_real64
      allocate(self%mroot(numnod));   self%mroot   = 0.0_real64
      allocate(self%hroot(numnod));   self%hroot   = 0.0_real64
      allocate(self%rootrho(numnod)); self%rootrho = 0.0_real64
      allocate(self%rootphi(numnod)); self%rootphi = 0.0_real64
      allocate(self%rmax(numnod));    self%rmax    = 0.0_real64

      ! Lookup table — allocated here; built by MatricFlux(1) in CropGrowth(1)
      ! (build is BLOCKED from relocation: ksatfit/wiltpoint/cofgen not yet
      ! populated when this routine runs — see ADR 0036 §mfluxtable disposition)
      allocate(self%mfluxtable(numlay, 801)); self%mfluxtable = 0.0_real64

      ! ===========================================================================
      ! SOIL-WATER CORE — flat per-node / per-layer arrays (ADR 0038, S-1.2)
      ! ===========================================================================

      ! Per-node arrays sized numnod (10 arrays):
      allocate(self%theta(numnod));        self%theta        = 0.0_real64
      allocate(self%thetm1(numnod));       self%thetm1       = 0.0_real64
      allocate(self%thetar(numnod));       self%thetar       = 0.0_real64
      allocate(self%thetas(numnod));       self%thetas       = 0.0_real64
      allocate(self%h(numnod));            self%h            = 0.0_real64
      allocate(self%hm1(numnod));          self%hm1          = 0.0_real64
      allocate(self%dimoca(numnod));       self%dimoca       = 0.0_real64
      ! [MACRO-RETIRE 2026-05-12] FrArMtrx defaults to 1.0 (whole-matrix);
      ! macropore retirement removed the only writer (MACROGEOM). ADR 0040.
      allocate(self%FrArMtrx(numnod));     self%FrArMtrx     = 1.0_real64
      allocate(self%evp(numnod));          self%evp          = 0.0_real64
      allocate(self%indeks(numnod));       self%indeks       = 0

      ! Per-node logical array sized numnod:
      allocate(self%fluseksatexm(numnod)); self%fluseksatexm = .false.

      ! Per-node arrays sized numnod+1 (flux arrays — one value per node boundary):
      ! Legacy: q(macp+1), k(macp+1), kmean(macp+1)
      allocate(self%q(numnod+1));          self%q            = 0.0_real64
      allocate(self%qssdi(numnod));        self%qssdi        = 0.0_real64  ! [GR-SOIL 2026-05-24] migrated from variables.f90
      if (.not. allocated(self%Itnumb)) then
         allocate(self%Itnumb(100, 2));    self%Itnumb       = 0            ! [GR-SOIL 2026-05-24] Richards iter stats
      end if
      allocate(self%k(numnod+1));          self%k            = 0.0_real64
      allocate(self%kmean(numnod+1));      self%kmean        = 0.0_real64

      allocate(self%vg_params(numnod))    ! [SS-GR-UTILS] components default-init from type
      ! [GR-CROP 2026-05-25] per-layer VG store for tillage mutator (replaces paramvg(21, maho))
      allocate(self%vg_params_layer(numlay))

      ! [SS-GR-UTILS] Soil hydraulic property metadata (Task 4)
      ! iHWCKmodel/layer/BiModal/NoVap allocated unconditionally (small, per-layer/node).
      ! [GR-SOIL 2026-05-24] numtab/sptab/ientrytab allocations retired — swsophy=1
      !   tabulated path is dormant; see src/soil/dormant/sptabulated.f90.
      ! [GR-SOIL 2026-05-24] iHWCKmodel default 1 (uni-modal MvG) — matches
      !   HACK Phase 4f-extend constraint (config_to_variables:612 forced =1).
      !   When the TOML schema exposes per-layer iHWCKmodel, set it here from config.
      allocate(self%iHWCKmodel(numlay));             self%iHWCKmodel = 1
      allocate(self%layer(numnod));                  self%layer      = 0
      allocate(self%BiModal(numlay));                self%BiModal    = .false.
      allocate(self%NoVap(numlay));                  self%NoVap      = .false.

      ! Per-layer array sized numlay (thetsl per soil layer, not per node):
      ! Legacy: thetsl(maho) where maho = max number of soil layers
      allocate(self%thetsl(numlay));         self%thetsl       = 0.0_real64

      ! [SS-GR-BH A4] Layer-flat fields (maho-sized, one per soil layer):
      allocate(self%ksatexm(numlay));        self%ksatexm      = 0.0_real64
      allocate(self%ksatfit(numlay));        self%ksatfit      = 0.0_real64
      allocate(self%cofani(numlay));         self%cofani       = 0.0_real64
      self%flksatexm = .false.
      allocate(self%orgmat(numlay));         self%orgmat       = 0.0_real64
      ! bdens: Pattern 9 guarded alloc (Piece B soil seeding may run before/after).
      if (.not. allocated(self%bdens)) then
         allocate(self%bdens(numlay));       self%bdens        = 0.0_real64
      end if
      allocate(self%psand(numlay));          self%psand        = 0.0_real64
      allocate(self%psilt(numlay));          self%psilt        = 0.0_real64
      allocate(self%pclay(numlay));          self%pclay        = 0.0_real64

      ! [SS-GR-BH A5] Runtime scalars (seeded from config in Task 7)
      self%q0              = 0.0_real64
      self%k1max           = 0.0_real64
      self%H0max           = 0.0_real64
      self%swbotb_runtime  = 0

      ! ===========================================================================
      ! SOIL-WATER CORE — intermediate-period arrays (ADR 0038, S-1.2)
      ! ===========================================================================

      ! Non-per-day per-node arrays sized numnod:
      allocate(self%inqrot(numnod));    self%inqrot    = 0.0_real64
      allocate(self%inqssdi(numnod));   self%inqssdi   = 0.0_real64
      allocate(self%IThetaBeg(numnod)); self%IThetaBeg = 0.0_real64

      ! Non-per-day per-node flux arrays sized numnod+1:
      ! Legacy: inq(macp+1), iqdo(macp+1), iqup(macp+1)
      allocate(self%inq(numnod+1));     self%inq       = 0.0_real64
      allocate(self%iqdo(numnod+1));    self%iqdo      = 0.0_real64
      allocate(self%iqup(numnod+1));    self%iqup      = 0.0_real64

      ! Per-day per-node arrays sized numnod:
      ! Legacy: qpotrot_day(macp), qredtot_day(macp)
      allocate(self%qpotrot_day(numnod)); self%qpotrot_day = 0.0_real64
      allocate(self%qredtot_day(numnod)); self%qredtot_day = 0.0_real64

      ! ===========================================================================
      ! NON-ZERO DEFAULTS (ADR 0038, S-1.2)
      ! ===========================================================================

      ! hatm: air pressure head near soil surface; legacy soilhydraulics.f90:870
      ! sets hatm = -2.75d+05 at SoilWater(1) task 1 init.  We mirror that here
      ! so self%hatm is consistent from the moment init returns.
      self%hatm = -2.75e5_real64

      ! Runon time-series (dormant — no TOML writer yet). Sized to MADAY to
      ! match the legacy fixed-size global runonarr(maday).
      allocate(self%runonarr(MADAY)); self%runonarr = 0.0_real64

      ! ---- [Piece B] Soil scalar seeding from config_soil ----

      self%swsophy     = config_soil%swsophy
      self%swinco      = config_soil%swinco
      self%swfrost     = config_soil%frost%swfrost
      self%swtill      = config_soil%swtill
      self%swhyst      = config_soil%swhyst
      self%swdiscrvert = config_soil%discretization%swdiscrvert
      ! Legacy parses .swp `SWRUNON` into a local int; we mirror that mapping
      ! into self%flrunon (runonarr remains dormant — no TOML writer).
      self%flrunon = (config_soil%swrunon == 1)

      ! bdens from soil.bdens (if present): Pattern 9 guarded alloc.
      if (allocated(config_soil%bdens)) then
         if (.not. allocated(self%bdens)) then
            allocate(self%bdens(maho)); self%bdens = 0.0d0
         end if
         do i = 1, size(config_soil%bdens)
            self%bdens(i) = config_soil%bdens(i)
         end do
      end if

      ! bdens from soil.hydraulics.bdens (if present): overrides/supplements soil.bdens.
      if (allocated(config_soil%hydraulics%ores)) then
         if (.not. allocated(self%bdens)) then
            allocate(self%bdens(maho)); self%bdens = 0.0d0
         end if
         do i = 1, size(config_soil%hydraulics%ores)
            self%bdens(i) = config_soil%hydraulics%bdens(i)
         end do
      end if

      ! [W3 fix 2026-05-28] swinco=3 warm-restart: h_file CSV → self%h_init.
      ! No longer mutates config_soil%initial%z_init (config is now read-only here).
      ! Consumers: swap_mod copies h values; soilhydraulics.f90 reads rows(i)%z.
      ! Cross-subsystem writes retained in adapter (state%atmosphere%atmin7,
      ! state%solute%X). See [GR-SEED 2026-05-25 Task 8] adapter comment.
      if (config_soil%swinco == 3) then
         if (allocated(config_soil%initial%h_file) .and. &
             len_trim(config_soil%initial%h_file) > 0) then
            block
               use error_mod, only: error_collection_t
               type(error_collection_t) :: errs
               call self%h_init%load(trim(config_soil%initial%h_file), errs)
               call errs%abort_if_fatal()
            end block
         end if
      end if

      ! ---- [Piece C] Bottom-boundary CSV pre-loads ----
      call seed_bottom_boundary(self, config_bb, pathwork)

      ! ---- [Piece D] Layer-flats from multiple config sources ----
      ! Folded from swap_init_body STRANGLER block (W1 closure, orchestrator-
      ! dissolution arc step 5). Multi-source resolution rules:
      !   * cofani: config_drain%cofani first, config_soil%cofani overrides (soil wins).
      !   * orgmat: config_soil%orgmat first; config_heat%porg backfills when soil absent.
      if (allocated(config_soil%hydraulics%ksatexm)) &
         self%ksatexm(:) = config_soil%hydraulics%ksatexm(1:size(self%ksatexm))
      if (allocated(config_soil%hydraulics%ksatfit)) &
         self%ksatfit(:) = config_soil%hydraulics%ksatfit(1:size(self%ksatfit))
      if (allocated(config_drain%cofani)) &
         self%cofani(1:size(config_drain%cofani)) = config_drain%cofani
      if (allocated(config_soil%cofani)) &
         self%cofani(1:size(config_soil%cofani))  = config_soil%cofani   ! soil wins
      self%flksatexm = .false.
      if (allocated(config_soil%orgmat)) then
         self%orgmat(1:size(config_soil%orgmat)) = config_soil%orgmat
      else if (allocated(config_heat%porg)) then
         self%orgmat(1:min(size(config_heat%porg), size(self%orgmat))) = &
            config_heat%porg(1:min(size(config_heat%porg), size(self%orgmat)))
      end if
      if (allocated(config_heat%psand)) &
         self%psand(:) = config_heat%psand(1:size(self%psand))
      if (allocated(config_heat%psilt)) &
         self%psilt(:) = config_heat%psilt(1:size(self%psilt))
      if (allocated(config_heat%pclay)) &
         self%pclay(:) = config_heat%pclay(1:size(self%pclay))
      self%swbotb_runtime = config_bb%swbotb
      self%q0    = 0.0d0
      self%k1max = 0.0d0
      self%H0max = 0.0d0

      ! Numerical solver control — snapshotted for boundtop (cluster 4) and
      ! soilhydraulics (cluster 7). Default 1 matches simulation_numerical_t.
      self%swkmean = config_simulation%numerical%swkmean

      ! [state%cfg-retirement cluster 7] Numerical solver snapshots (soilhydraulics).
      self%swkimpl                   = config_simulation%numerical%swkimpl
      self%MaxBackTr                 = config_simulation%numerical%MaxBackTr
      self%gwlconv                   = config_simulation%numerical%gwlconv
      self%critdevh1cp               = config_simulation%numerical%critdevh1cp
      self%critdevh2cp               = config_simulation%numerical%critdevh2cp
      self%critdevponddt             = config_simulation%numerical%critdevponddt
      self%swcaprise                 = config_simulation%numerical%swcaprise
      self%dump_convergence_diagnostics = config_simulation%numerical%dump_convergence_diagnostics

      ! [state%cfg-retirement cluster 7] Soil config scalar snapshots.
      self%gwli = config_soil%gwli
      self%tau  = config_soil%tau

      ! [state%cfg-retirement cluster 7] Bottom-boundary scalar snapshots.
      self%hplate      = config_bb%hplate
      self%rimlay      = config_bb%rimlay
      self%swbotb3impl = config_bb%swbotb3impl
      self%sw4         = config_bb%sw4

      ! [state%cfg-retirement cluster 7] Per-layer wetting alpha for hysteresis.
      if (allocated(config_soil%hydraulics%alfaw)) then
         if (.not. allocated(self%alfaw_layer)) then
            allocate(self%alfaw_layer(size(config_soil%hydraulics%alfaw)))
         end if
         self%alfaw_layer(:) = config_soil%hydraulics%alfaw(:)
      end if

   end subroutine soilwater_state_init

   !> Private helper: pre-load bottom-boundary CSV tables into self based on swbotb.
   !! Called from soilwater_state_init (Piece C).
   !! [GR-SEED 2026-05-25 Task 8] Absorbed from config_to_variables bottom_boundary block.
   subroutine seed_bottom_boundary(self, config_bb, pathwork)
      use bottom_boundary_config_mod, only: bottom_boundary_config_t
      type(soilwater_state_t),        intent(inout) :: self
      type(bottom_boundary_config_t), intent(in)    :: config_bb
      character(len=*),               intent(in)    :: pathwork

      select case (config_bb%swbotb)
      case (1)
         block
            use boundary_csv_mod, only: gwl_table_t
            use error_mod, only: error_collection_t
            type(gwl_table_t)        :: loader
            type(error_collection_t) :: csv_errs
            integer :: k
            call loader%load(trim(pathwork)//trim(config_bb%gwl_file), csv_errs)
            call csv_errs%abort_if_fatal()
            if (loader%is_loaded) then
               if (.not. allocated(self%gwltab)) then
                  allocate(self%gwltab(2*MABBC))
                  self%gwltab = 0.0d0
               end if
               do k = 1, min(size(loader%rows), MABBC)
                  self%gwltab(k*2 - 1) = loader%rows(k)%date
                  self%gwltab(k*2)     = loader%rows(k)%gwl
               end do
            end if
         end block
      case (2)
         ! [GR-IO 2026-05-25 Phase 6 Step 3] sw2 legacy mirror dropped
         ! Phase 0 B-0.1: populate sine-wave scalars regardless of sw2;
         ! the gate in boundbottom.f90:104 protects the non-sine path.
         ! [GR-IO 2026-05-25 Phase 6 Step 3] sinmax/sinamp/sinave legacy mirrors dropped
         if (config_bb%sw2 == 2) then
            block
               use boundary_csv_mod, only: qbot_table_t
               use error_mod, only: error_collection_t
               type(qbot_table_t)       :: loader
               type(error_collection_t) :: csv_errs
               integer :: k
               call loader%load(trim(pathwork)//trim(config_bb%qbot2_file), csv_errs)
               call csv_errs%abort_if_fatal()
               if (loader%is_loaded) then
                  if (.not. allocated(self%qbotab)) then
                     allocate(self%qbotab(2*MABBC))
                     self%qbotab = 0.0d0
                  end if
                  do k = 1, min(size(loader%rows), MABBC)
                     self%qbotab(k*2 - 1) = loader%rows(k)%date
                     self%qbotab(k*2)     = loader%rows(k)%qbot
                  end do
               end if
            end block
         end if
      case (3)
         ! [GR-DRAIN 2026-05-25] shape legacy mirror dropped — boundbottom.f90
         ! reads bb%shape (= config%bottom_boundary%shape) directly.
         ! [GR-IO 2026-05-25 Phase 6 Step 3] hdrain/aqave/aqamp/aqper/aqtmax/sw3
         ! legacy mirrors dropped (boundbottom.f90 reads bb%X = config%bottom_boundary%X).
         ! [GR-SOIL 2026-05-24] rimlay/swbotb3impl already dropped.
         ! [GR-SOIL 2026-05-24] sw4 legacy mirror dropped — direct config read.
         if (config_bb%sw3 == 2) then
            block
               use boundary_csv_mod, only: haquif_table_t
               use error_mod, only: error_collection_t
               type(haquif_table_t)     :: loader
               type(error_collection_t) :: csv_errs
               integer :: k
               call loader%load(trim(pathwork)//trim(config_bb%haquif_file), csv_errs)
               call csv_errs%abort_if_fatal()
               if (loader%is_loaded) then
                  if (.not. allocated(self%haqtab)) then
                     allocate(self%haqtab(2*MABBC))
                     self%haqtab = 0.0d0
                  end if
                  do k = 1, min(size(loader%rows), MABBC)
                     self%haqtab(k*2 - 1) = loader%rows(k)%date
                     self%haqtab(k*2)     = loader%rows(k)%haquif
                  end do
               end if
            end block
         end if
         if (config_bb%sw4 == 1) then
            block
               use boundary_csv_mod, only: qbot_table_t
               use error_mod, only: error_collection_t
               type(qbot_table_t)       :: loader
               type(error_collection_t) :: csv_errs
               integer :: k
               call loader%load(trim(pathwork)//trim(config_bb%qbot4_file), csv_errs)
               call csv_errs%abort_if_fatal()
               if (loader%is_loaded) then
                  if (.not. allocated(self%qbotab)) then
                     allocate(self%qbotab(2*MABBC))
                     self%qbotab = 0.0d0
                  end if
                  do k = 1, min(size(loader%rows), MABBC)
                     self%qbotab(k*2 - 1) = loader%rows(k)%date
                     self%qbotab(k*2)     = loader%rows(k)%qbot
                  end do
               end if
            end block
         end if
      case (4)
         ! [GR-IO 2026-05-25 Phase 6 Step 3] swqhbot/cofqha/cofqhb/cofqhc/swcofqhc
         ! legacy mirrors dropped — boundbottom reads via bb%X.
         if (config_bb%swqhbot == 2) then
            block
               use boundary_csv_mod, only: qhbot_table_t
               use error_mod, only: error_collection_t
               type(qhbot_table_t)      :: loader
               type(error_collection_t) :: csv_errs
               integer :: k
               call loader%load(trim(pathwork)//trim(config_bb%qhbot_file), csv_errs)
               call csv_errs%abort_if_fatal()
               ! Legacy unpack pattern from readswap.f90:1418-1419 — for the
               ! q(h) curve, qbotab(odd) = abs(htab) and qbotab(even) = qtab.
               if (loader%is_loaded) then
                  if (.not. allocated(self%qbotab)) then
                     allocate(self%qbotab(2*MABBC))
                     self%qbotab = 0.0d0
                  end if
                  do k = 1, min(size(loader%rows), MABBC)
                     self%qbotab(k*2 - 1) = abs(loader%rows(k)%htab)
                     self%qbotab(k*2)     = loader%rows(k)%qtab
                  end do
               end if
            end block
         end if
      case (5)
         ! [SS-BND B-2.7] hbot global retired; state%soilwater%hbot set by boundbottom each step.
         ! hbot = config%bottom_boundary%hbot
         ! NOTE: rhobot has no legacy SWAP-wide global; the plan's spec
         ! line `rhobot = config%bottom_boundary%rhobot` was a defect.
         ! The schema slot is read for future-proofing; consumers TBD.
         block
            use boundary_csv_mod, only: hbot_table_t
            use error_mod, only: error_collection_t
            type(hbot_table_t)       :: loader
            type(error_collection_t) :: csv_errs
            integer :: k
            call loader%load(trim(pathwork)//trim(config_bb%hbot5_file), csv_errs)
            call csv_errs%abort_if_fatal()
            if (loader%is_loaded) then
               if (.not. allocated(self%hbotab)) then
                  allocate(self%hbotab(2*MABBC))
                  self%hbotab = 0.0d0
               end if
               do k = 1, min(size(loader%rows), MABBC)
                  self%hbotab(k*2 - 1) = loader%rows(k)%date
                  self%hbotab(k*2)     = loader%rows(k)%hbot
               end do
            end if
         end block
      case (6, 7)
         ! No parameters to populate for modes 6 and 7.
      case (8)
         ! [GR-SOIL 2026-05-24] hplate legacy mirror dropped — direct config read.
      end select

      ! Bottom-boundary CSV tables: size 2*MABBC to match legacy globals.
      ! Guard allocations for tables not populated by the case dispatch above.
      if (.not. allocated(self%gwltab)) then
         allocate(self%gwltab(2*MABBC)); self%gwltab = 0.0_real64
      end if
      if (.not. allocated(self%qbotab)) then
         allocate(self%qbotab(2*MABBC)); self%qbotab = 0.0_real64
      end if
      if (.not. allocated(self%haqtab)) then
         allocate(self%haqtab(2*MABBC)); self%haqtab = 0.0_real64
      end if
      if (.not. allocated(self%hbotab)) then
         allocate(self%hbotab(2*MABBC)); self%hbotab = 0.0_real64
      end if

   end subroutine seed_bottom_boundary

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
