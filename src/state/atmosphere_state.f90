!> @file atmosphere_state.f90
!! Typed state record for the atmosphere subsystem (ADR 0037).
!!
!! Debuts the ADR 0033 cohort pattern at type-creation — two nested
!! sub-records replace 3 scattered reset blocks that currently live in
!! meteoday.f90, soilhydraulics.f90, and snow.f90:
!!
!!   atmosphere_intermediate_t (8 fields, flzerointr-reset)
!!     igrai, inrai, ipeva, iptra, ievap, igsnow, isubl, isnrai
!!
!!   atmosphere_cumulative_t (10 fields, flzerocumu-reset)
!!     cgrai, cnrai, caintc, cpeva, cptra, cevap, cgsnow, csubl, csnrai, cmelt
!!
!! atmosphere_state_t holds:
!!   22 flat top-level scalars (11 instantaneous + 9 per-day + 3 per-event)
!!   + type(atmosphere_intermediate_t) :: intr
!!   + type(atmosphere_cumulative_t)   :: cumu
!!   = 40 total fields, all scalars (no per-node arrays).
!!
!! state%atmosphere%init(config) zeroes flat scalars + cohorts and copies config-derived params. Cohort
!! sub-records are already zero-defaulted at declaration; init re-zeroes
!! them explicitly for forward-compat with future arcs that might add
!! allocatable arrays.
!!
!! Excluded:
!!   - nird/gird — irrigation-owned; deferred to irrigation arc.
!!   - pond      — soil-water-core arc territory (boundary D5 deferral).
!!   - Orchestration structure of meteoday.f90/meteodt.f90 — future
!!     atmosphere REFACTOR arc.
!!
!! See ADR 0037.

module atmosphere_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   use swap_array_dimensions, only: magrs, mayrs, mrain
   use meteo_csv_mod, only: meteo_daily_table_t, meteo_detail_table_t, rain_events_table_t
   implicit none
   private
   public :: atmosphere_state_t
   public :: atmosphere_intermediate_t
   public :: atmosphere_cumulative_t

   !> Intermediate accumulators — reset when flzerointr fires.
   !! These 8 fields accumulate within the intermediate output interval;
   !! their reset block currently lives in three places (meteoday.f90,
   !! soilhydraulics.f90, waterbalance.f90). After migration, the single
   !! owner is intr%reset() called under the flzerointr gate.
   type :: atmosphere_intermediate_t
      real(real64) :: igrai  = 0.0_real64   !< intermediate gross rainfall (cm)
      real(real64) :: inrai  = 0.0_real64   !< intermediate net rainfall (cm)
      real(real64) :: ipeva  = 0.0_real64   !< intermediate potential soil evaporation (cm)
      real(real64) :: iptra  = 0.0_real64   !< intermediate potential transpiration (cm)
      real(real64) :: ievap  = 0.0_real64   !< intermediate actual evaporation (cm)
      real(real64) :: igsnow = 0.0_real64   !< intermediate gross snowfall (cm)
      real(real64) :: isubl  = 0.0_real64   !< intermediate sublimation (cm)
      real(real64) :: isnrai = 0.0_real64   !< intermediate snow-rain partition (cm)
   contains
      procedure :: reset => atmosphere_intermediate_reset
   end type atmosphere_intermediate_t

   !> Cumulative accumulators — reset when flzerocumu fires.
   !! These 10 fields accumulate over the full simulation output interval;
   !! their reset block is currently scattered across meteoday.f90,
   !! soilhydraulics.f90, and snow.f90. After migration, the single owner
   !! is cumu%reset() called under the flzerocumu gate.
   type :: atmosphere_cumulative_t
      real(real64) :: cgrai  = 0.0_real64   !< cumulative gross rainfall (cm)
      real(real64) :: cnrai  = 0.0_real64   !< cumulative net rainfall (cm)
      real(real64) :: caintc = 0.0_real64   !< cumulative actual interception (cm)
      real(real64) :: cpeva  = 0.0_real64   !< cumulative potential soil evaporation (cm)
      real(real64) :: cptra  = 0.0_real64   !< cumulative potential transpiration (cm)
      real(real64) :: cevap  = 0.0_real64   !< cumulative actual evaporation (cm)
      real(real64) :: cgsnow = 0.0_real64   !< cumulative gross snowfall (cm)
      real(real64) :: csubl  = 0.0_real64   !< cumulative sublimation (cm)
      real(real64) :: csnrai = 0.0_real64   !< cumulative snow-rain partition (cm)
      real(real64) :: cmelt  = 0.0_real64   !< cumulative snowmelt (cm)
   contains
      procedure :: reset => atmosphere_cumulative_reset
   end type atmosphere_cumulative_t

   !> Top-level atmosphere state record. 40 fields total:
   !! 22 flat scalars (11 instantaneous + 9 per-day + 3 per-event)
   !! + 8 intermediate cohort fields (intr)
   !! + 10 cumulative cohort fields (cumu).
   !!
   !! All fields are scalars — no per-node arrays.
   !! ASSOCIATE prefix at_  is used in large compute bodies.
   type :: atmosphere_state_t

      ! -----------------------------------------------------------------------
      ! Instantaneous scalars (11) — updated every timestep
      ! -----------------------------------------------------------------------
      real(real64) :: peva     = 0.0_real64  !< potential soil evaporation (cm/d)
      real(real64) :: ptra     = 0.0_real64  !< potential transpiration (cm/d)
      ! Reference-ET intermediates (mm/d): computed by the ET step from reference
      ! ET x crop factor and consumed within the same step to form peva/ptra.
      ! Atmosphere-owned (relocated off crop state, T2-A crop_water ET sub-part).
      real(real64) :: et0      = 0.0_real64  !< potential ET (mm/d)
      real(real64) :: ew0      = 0.0_real64  !< potential evaporation of wet crop (mm/d)
      real(real64) :: es0      = 0.0_real64  !< potential evaporation of bare soil (mm/d)
      real(real64) :: empreva  = 0.0_real64  !< actual soil evaporation after reduction (cm/d)
      real(real64) :: melt     = 0.0_real64  !< snowmelt this timestep (cm)
      real(real64) :: subl     = 0.0_real64  !< sublimation this timestep (cm)
      real(real64) :: slw      = 0.0_real64  !< liquid water in snowpack (cm)
      real(real64) :: ssnow    = 0.0_real64  !< snow storage (cm water equivalent)
      real(real64) :: snowinco = 0.0_real64  !< snow-in-canopy snapshot (cm w.e.)
      real(real64) :: ISsnowBeg = 0.0_real64 !< snow w.e. at start of current intermediate period (cm); peer to state%soilwater%IPondBeg

      ! Config-derived parameters/switches (copied from config at init; compute reads from state only)
      real(real64) :: snowcoef = 0.0_real64  !< snow-melt temperature coefficient (cm/d/degC)
      integer      :: swsublim = 0           !< suppress sublimation of snow (1) or compute it (0)
      ! [state%cfg-retirement cluster 7] snow physics switch; snapshotted from config%meteo%snow%swsnow
      integer      :: swsnow   = 0           !< snow physics switch: 0=off, 1=on
      integer      :: swetsine = 0           !< Tp/Ep distribution: 0=uniform, 1=sine-wave during day
      integer      :: swredu     = 1           !< ET reduction method: 1=Black, 2=Boesten-Stroosnijder
      real(real64) :: cofred     = 0.35_real64 !< ET reduction coefficient β (Black: cofredbl, B-S: cofredbo)
      real(real64) :: rsigni     = 0.5_real64  !< Significant rainfall threshold resetting Black dry counter (cm/d)
      real(real64) :: angstroma  = 0.25_real64 !< Ångström a coefficient (-) for atmospheric transmission
      real(real64) :: angstromb  = 0.50_real64 !< Ångström b coefficient (-)
      real(real64) :: cfevappond = 1.25_real64 !< pond-evap / ETref ratio (-)
      real(real64) :: rsoil      = 0.0_real64  !< soil resistance of wet soil for PMdirect (s/m)
      real(real64) :: lat        = 0.0_real64  !< geographic latitude (degrees N); snapshotted from config%meteo%lat

      ! Today's astro outputs — cached once per day in ReadMeteoDay (daily mode)
      real(real64) :: daylp        = 12.0_real64  !< photoperiodic daylength (h)
      real(real64) :: difpp        = 0.0_real64   !< diffuse irradiation perpendicular to light direction (J/m2/s)
      real(real64) :: atmtr        = 0.0_real64   !< daily atmospheric transmission (-)
      real(real64) :: dsinbe       = 0.0_real64   !< daily total of effective solar height (s)
      real(real64) :: tsunrise_atm = 0.0_real64   !< time of sunrise (fraction of day)
      real(real64) :: tsunset_atm  = 1.0_real64   !< time of sunset (fraction of day)

      ! Dormant feature flags (no writers in current TOML pipeline; gated branches are dead).
      ! Kept in state to preserve compute-routine signatures; retire fully when the
      ! relevant features (CN runoff, CO2 effects) get wired into config.
      integer :: swusecn = 0  !< Curve-Number runoff method: 0=off, 1=on (no TOML wiring yet)
      logical :: flco2   = .false.  !< CO2 effects toggle (no TOML wiring yet)
      real(real64) :: graidt   = 0.0_real64  !< gross rainfall this timestep (cm)
      real(real64) :: nraidt   = 0.0_real64  !< net rainfall this timestep (cm)
      real(real64) :: aintcdt  = 0.0_real64  !< actual interception this timestep (cm)
      real(real64) :: nird     = 0.0_real64  !< net irrigation depth after interception (cm) — peer of nraida; may re-home to state%irrigation%nird in future arc

      ! Cross-routine rain/interception runtime state (atmosphere writes, downstream reads)
      real(real64) :: finterception = 1.0_real64  !< net/gross rain ratio after interception (-); written by meteo_orchestrator per day, read in meteodt per timestep
      real(real64) :: dtEventRain   = 0.0_real64  !< time-step length until next precipitation event (d); written/read in meteodt

      ! Today's daily-meteo scalars (set by ReadMeteoDay or sub-daily aggregator)
      real(real64) :: rad = 0.0_real64  !< today's daily-total radiation (J/m2/d)
      real(real64) :: tmn = 0.0_real64  !< today's minimum air temperature (°C)
      real(real64) :: tmx = 0.0_real64  !< today's maximum air temperature (°C)

      ! -----------------------------------------------------------------------
      ! Per-day scalars (9) — updated once per day in ProcessMeteoDay
      ! -----------------------------------------------------------------------
      real(real64) :: grai        = 0.0_real64  !< gross daily rainfall (cm)
      real(real64) :: nraida      = 0.0_real64  !< net daily rainfall after interception (cm)
      real(real64) :: atmdem      = 0.0_real64  !< atmospheric demand (cm/d)
      real(real64) :: pevaday     = 0.0_real64  !< potential evaporation for the day (cm/d)
      real(real64) :: ptraday     = 0.0_real64  !< potential transpiration for the day (cm/d)
      real(real64) :: gsnow       = 0.0_real64  !< gross daily snowfall (cm)
      real(real64) :: snrai       = 0.0_real64  !< daily snow-rain partition (cm)
      real(real64) :: fprecnosnow = 0.0_real64  !< fraction of precipitation that is rain (-)
      real(real64) :: sicact      = 0.0_real64  !< actual interception storage capacity (cm)

      ! -----------------------------------------------------------------------
      ! Per-event / tillage-reset scalars (3)
      ! -----------------------------------------------------------------------
      real(real64) :: ldwet = 0.0_real64  !< leaf area index for wet canopy (-)
      real(real64) :: spev  = 0.0_real64  !< soil evaporation reduction factor (-)
      real(real64) :: saev  = 0.0_real64  !< actual soil evaporation after reduction (cm/d)

      ! -----------------------------------------------------------------------
      ! Cohort sub-records
      ! -----------------------------------------------------------------------
      type(atmosphere_intermediate_t) :: intr  !< intermediate accumulators (flzerointr-reset)
      type(atmosphere_cumulative_t)   :: cumu  !< cumulative accumulators (flzerocumu-reset)

      ! [SS-GR-ATM A3] Block 1: daily meteo input arrays (366-sized fixed)
      real(real64) :: arad(366) = 0.0_real64    !! daily radiation input
      real(real64) :: atmn(366) = 0.0_real64    !! daily min temperature
      real(real64) :: atmx(366) = 0.0_real64    !! daily max temperature
      real(real64) :: ahum(366) = 0.0_real64    !! daily humidity
      real(real64) :: awin(366) = 0.0_real64    !! daily wind speed
      real(real64) :: arai(366) = 0.0_real64    !! daily rainfall
      real(real64) :: aetr(366) = 0.0_real64    !! daily ETref
      real(real64) :: wet(366)  = 0.0_real64    !! daily wet-fraction

      ! [SS-GR-ATM A4] Block 2: sub-daily detailed arrays (96-sized; swmetdetail==1)
      real(real64) :: atav(96)  = 0.0_real64    !! sub-daily air temp
      real(real64) :: epot(96)  = 0.0_real64    !! sub-daily potential evap
      real(real64) :: tpot(96)  = 0.0_real64    !! sub-daily potential transp
      real(real64) :: grain(96) = 0.0_real64    !! sub-daily gross rain
      real(real64) :: nrain(96) = 0.0_real64    !! sub-daily net rain

      ! Per-year detail-meteo input arrays (allocated when swmetdetail==1).
      ! Populated by MeteoCSVDetYear (readmeteo.f90), consumed by ReadMeteoDay.
      integer,      allocatable :: detrecord(:)  !! record number per sub-daily entry
      real(real64), allocatable :: dethum(:)     !! humidity per sub-daily entry
      real(real64), allocatable :: detrad(:)     !! radiation per sub-daily entry (J/m2)
      real(real64), allocatable :: detrain(:)    !! precipitation per sub-daily entry (mm)
      real(real64), allocatable :: dettav(:)     !! temperature per sub-daily entry (°C)
      real(real64), allocatable :: dettime(:)    !! timestamp per sub-daily entry
      real(real64), allocatable :: detwind(:)    !! wind speed per sub-daily entry
      integer :: irectotal = 0                   !! cumulative sub-daily record counter

      ! [GR-ATM-CLEAN Phase D] migrated from module MeteoVars (meteo_vars.f90).
      !> Remaining interception storage at start of timestep (cm).
      !> Persists across iterations of the sub-daily dayparts loop and across days.
      real(real64) :: restint = 0.0_real64

      !> Sub-daily precipitation per record (cm). Populated by ReadMeteoDay,
      !> consumed by ProcessMeteoDay's sub-daily branch. Allocation size = nmetdetail (<= 96).
      real(real64) :: arain_subdaily(96) = 0.0_real64

      !> Sub-daily wind speed per record (m/s). Same lifetime as arain_subdaily.
      real(real64) :: awind_subdaily(96) = 0.0_real64

      ! [SS-GR-ATM A5] Block 3: derived meteo scalars
      real(real64) :: Tav        = 0.0_real64    !! daily mean air temp
      real(real64) :: tavd       = 0.0_real64    !! daytime mean air temp
      real(real64) :: rh         = 0.0_real64    !! relative humidity
      integer      :: daynrfirst = 0             !! first day in meteo year
      integer      :: daynrlast  = 0             !! last day in meteo year
      real(real64) :: atmin7(7)  = 0.0_real64    !! 7-day min-temp buffer
      integer      :: nofd       = 0             !! current day-of-running-avg
      real(real64) :: teprrain   = 0.0_real64    !! threshold rain temp
      real(real64) :: teprsnow   = 0.0_real64    !! threshold snow temp

      ! [SS-GR-ATM A6] Block 4: interception state/params
      real(real64) :: siccapact          = 0.0_real64
      real(real64) :: fimin              = 0.0_real64
      integer      :: isua               = 0
      real(real64) :: avevaptb(2*magrs)  = 0.0_real64  !! actual evap table
      real(real64) :: avprectb(2*magrs)  = 0.0_real64  !! actual precip table
      real(real64) :: pfreetb(2*magrs)   = 0.0_real64  !! free throughfall table
      real(real64) :: pstemtb(2*magrs)   = 0.0_real64  !! stemflow table
      real(real64) :: scanopytb(2*magrs) = 0.0_real64  !! canopy storage table

      ! [SS-GR-ATM A7] Block 5: Runoff-CN method state + tables
      ! CNref/CNdry/CNwet retired: pure within-call intermediates of cn_step (locals).
      real(real64) :: ThetaRef          = 0.0_real64
      real(real64) :: Runoff_CN         = 0.0_real64
      integer      :: wc_cor            = 0
      real(real64) :: wc10              = 0.0_real64
      integer      :: iCNtab            = 0
      real(real64) :: CNtimTAB(mayrs*5) = 0.0_real64   !! CN time table (legacy dim: mayrs*5=1000)
      real(real64) :: CNrefTAB(mayrs*5) = 0.0_real64   !! CN ref  table (legacy dim: mayrs*5=1000)
      ! [SS-GR-ATM B22] CN method internal state (formerly variables.f90 module-level)
      integer      :: nod10_cn          = 0             !! node index at ~10 cm depth for CN method
      integer      :: icn_atm           = 0             !! current position in CN time series
      real(real64) :: z10_cn            = 0.0_real64    !! depth to nod10_cn node (cm)

      ! [SS-GR-ATM A8] Block 6: daily output scalars (real(4) in variables.f90;
      ! real(real64) here — promoted for state consistency; meteoday.f90 writes
      ! these per day before CSV output; dual-write mirrors legacy global values).
      real(real64) :: out_tmn = 0.0_real64   !! min air temperature of current day (oC)
      real(real64) :: out_tmx = 0.0_real64   !! max air temperature of current day (oC)
      real(real64) :: out_hum = 0.0_real64   !! air humidity of current day (kPa)
      real(real64) :: out_win = 0.0_real64   !! average wind speed of current day (m/s)
      real(real64) :: out_etr = 0.0_real64   !! reference ET of current day (m/d)
      real(real64) :: out_wet = 0.0_real64   !! rainfall duration of current day (d)
      real(real64) :: out_rad = 0.0_real64   !! global solar radiation (kJ/m2)

      ! [SS-GR-CROP A12] rain timing — per-year reload runtime state
      integer      :: nmrain                    = 0          !! Number of rain event records (-)
      real(real64) :: rainamount(mrain)         = 0.0_real64 !! Array with short duration rainfall amounts (L)
      real(real64) :: rainfluxarray(mrain)      = 0.0_real64 !! Array with short duration rainfall intensities (L/T)
      real(real64) :: raintimearray(mrain)      = 0.0_real64 !! Array with times (T) at which rainfall intensity changes

      ! [GR-IO 2026-05-25 Phase 3] Per-year meteo-file date columns.
      ! Written by MeteoCSVYear / read_meteo_from_external_buffer_year,
      ! read by ReadMeteoYear for validation + raintimearray init.
      integer :: ad(mrain) = 0  !! day-of-month per meteo-file day
      integer :: am(mrain) = 0  !! month per meteo-file day

      ! [METEO-TYPED-CSV 2026-05-28] Typed CSV table records replace the
      ! raw (:,:) caches. Per-year extractors in readmeteo.f90 use
      ! %year_window(year, i1, i2) + %rows(i1:i2)%field accessors.
      type(meteo_daily_table_t)  :: meteo         !! daily meteo CSV cache
      type(meteo_detail_table_t) :: meteo_detail  !! sub-daily meteo CSV cache (swmetdetail=1)
      type(rain_events_table_t)  :: rain_events   !! rain events CSV cache (swrain=3)

      ! [state%cfg-retirement mop-up] rain intensity lookup table (swrain==1).
      ! Fixed-size 60-element table; snapshotted from config%meteo%raintab at init.
      real(real64) :: raintab(60) = 0.0_real64   !! rain intensity (cm/d) vs time (T) table

   contains
      procedure :: init => atmosphere_state_init
   end type atmosphere_state_t

contains

   !> Zero all flat scalars + cohort sub-records, then copy config-derived
   !! parameters into the state. Compute routines read these via state%X,
   !! never via state%cfg — the config→state boundary lives here.
   !! Type-bound init — call as state%atmosphere%init(config).
   subroutine atmosphere_state_init(self, config)
      use swap_config_mod, only: swap_config_t
      class(atmosphere_state_t), intent(inout) :: self
      type(swap_config_t),       intent(in)    :: config

      ! Instantaneous (11)
      self%peva     = 0.0_real64
      self%ptra     = 0.0_real64
      self%et0      = 0.0_real64
      self%ew0      = 0.0_real64
      self%es0      = 0.0_real64
      self%empreva  = 0.0_real64
      self%melt     = 0.0_real64
      self%subl     = 0.0_real64
      self%slw      = 0.0_real64
      self%ssnow    = 0.0_real64
      self%snowinco = 0.0_real64
      self%graidt   = 0.0_real64
      self%nraidt   = 0.0_real64
      self%aintcdt  = 0.0_real64

      ! Per-day (9)
      self%grai        = 0.0_real64
      self%nraida      = 0.0_real64
      self%atmdem      = 0.0_real64
      self%pevaday     = 0.0_real64
      self%ptraday     = 0.0_real64
      self%gsnow       = 0.0_real64
      self%snrai       = 0.0_real64
      self%fprecnosnow = 0.0_real64
      self%sicact      = 0.0_real64

      ! Per-event (3)
      self%ldwet = 0.0_real64
      self%spev  = 0.0_real64
      self%saev  = 0.0_real64

      ! Cohort sub-records
      call self%intr%reset()
      call self%cumu%reset()

      ! Config-derived snapshot — compute routines read these via state, not state%cfg.
      self%snowcoef = config%meteo%snow%snowcoef
      self%swsublim = config%soil%frost%swsublim
      self%swsnow   = config%meteo%snow%swsnow
      self%swetsine = config%meteo%swetsine
      self%swredu   = config%meteo%evaporation%swredu
      self%rsigni   = config%meteo%evaporation%rsigni
      ! cofredbl (Black) vs cofredbo (Boesten-Stroosnijder): pick by swredu.
      ! Config validation enforces that swredu=2 has a non-default cofredbo.
      if (config%meteo%evaporation%swredu == 2) then
         self%cofred = config%meteo%evaporation%cofredbo
      else
         self%cofred = config%meteo%evaporation%cofredbl
      end if
      self%angstroma  = config%meteo%angstroma
      self%angstromb  = config%meteo%angstromb
      self%cfevappond = config%meteo%evaporation%cfevappond
      self%rsoil      = config%soil%rsoil
      self%lat        = config%meteo%lat

      ! [METEO-TYPED-CSV 2026-05-28] Three typed table loads. Each owns
      ! its schema; validation is inline in load(); is_loaded flag
      ! distinguishes "not requested" from "empty".
      block
         use error_mod, only: error_collection_t
         character(len=300) :: metfile_lc, csvpath
         type(error_collection_t) :: errs_daily, errs_detail, errs_rain

         metfile_lc = ''
         if (allocated(config%meteo%metfile)) metfile_lc = config%meteo%metfile
         call lowerc(metfile_lc)

         ! Daily meteo table is now read at config-load time (config%meteo%meteo,
         ! decoupled from init). Copy it onto state — no file I/O here.
         self%meteo = config%meteo%meteo

         if (config%meteo%swmetdetail == 1) then
            if (allocated(config%meteo%detail_file) .and. &
                len_trim(config%meteo%detail_file) > 0) then
               csvpath = trim(config%general%pathatm) // trim(config%meteo%detail_file)
               call self%meteo_detail%load(trim(csvpath), errs_detail)
               call errs_detail%abort_if_fatal()
            end if
         end if

         if (config%meteo%swrain == 3 .and. allocated(config%meteo%rain_events_file)) then
            if (len_trim(config%meteo%rain_events_file) > 0) then
               csvpath = trim(config%general%pathatm) // trim(config%meteo%rain_events_file)
               call self%rain_events%load(trim(csvpath), errs_rain)
               call errs_rain%abort_if_fatal()
            end if
         end if
      end block

      ! [state%cfg-retirement mop-up] rain intensity table — ProcessRainEvents reads
      ! this via atmo%raintab (swrain==1 path).
      self%raintab(:) = config%meteo%raintab(:)

      ! Snow scalars.
      self%TePrRain = config%meteo%snow%teprrain
      self%TePrSnow = config%meteo%snow%teprsnow

      ! ---------------------------------------------------------------
      ! Warm-restart (swinco==3) fields. Folded from swap_init_body's
      ! orchestrator-dissolution arc (Steps 6+7).
      ! atmin7 is a fixed-size array in config; no allocated() guard needed.
      ! ---------------------------------------------------------------
      if (config%soil%swinco == 3) then
         self%ssnow = config%soil%initial%ssnow
         self%ldwet = config%soil%initial%ldwet
         self%slw   = config%soil%initial%slw
         self%atmin7(:) = config%soil%initial%atmin7(:)
      end if

      ! Snow override: ssnow=0 when snow physics disabled.
      if (config%meteo%snow%swsnow /= 1) self%ssnow = 0.0_real64

      ! Snow handshake: when snow physics on, snowinco mirrors ssnow at init.
      ! Under swinco==3, snowinco := ssnow (loaded value). Else ssnow := snowinco
      ! (which is 0 from the type-default since no other writer exists at init time).
      if (config%meteo%snow%swsnow == 1) then
         if (config%soil%swinco == 3) then
            self%snowinco = self%ssnow
         else
            self%ssnow = self%snowinco
         end if
      end if

   end subroutine atmosphere_state_init

   !> Zero all 8 intermediate cohort fields.
   !! Called under the flzerointr gate (replaces 3 scattered reset blocks).
   subroutine atmosphere_intermediate_reset(self)
      class(atmosphere_intermediate_t), intent(inout) :: self
      self%igrai  = 0.0_real64
      self%inrai  = 0.0_real64
      self%ipeva  = 0.0_real64
      self%iptra  = 0.0_real64
      self%ievap  = 0.0_real64
      self%igsnow = 0.0_real64
      self%isubl  = 0.0_real64
      self%isnrai = 0.0_real64
   end subroutine atmosphere_intermediate_reset

   !> Zero all 10 cumulative cohort fields.
   !! Called under the flzerocumu gate (replaces 3 scattered reset blocks).
   subroutine atmosphere_cumulative_reset(self)
      class(atmosphere_cumulative_t), intent(inout) :: self
      self%cgrai  = 0.0_real64
      self%cnrai  = 0.0_real64
      self%caintc = 0.0_real64
      self%cpeva  = 0.0_real64
      self%cptra  = 0.0_real64
      self%cevap  = 0.0_real64
      self%cgsnow = 0.0_real64
      self%csubl  = 0.0_real64
      self%csnrai = 0.0_real64
      self%cmelt  = 0.0_real64
   end subroutine atmosphere_cumulative_reset

end module atmosphere_state_mod
