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
!! state%atmosphere%init(config%meteo) zeroes all 22 flat scalars explicitly. Cohort
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
!! See ADR 0037, docs/superpowers/specs/2026-05-11-state-migration-atmosphere-design.md
!!     docs/superpowers/plans/2026-05-11-atmosphere-state-migration.md

module atmosphere_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   use iso_c_binding, only: c_double
   use swap_array_dimensions, only: magrs, mayrs, mrain
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
      real(real64) :: empreva  = 0.0_real64  !< actual soil evaporation after reduction (cm/d)
      real(real64) :: melt     = 0.0_real64  !< snowmelt this timestep (cm)
      real(real64) :: subl     = 0.0_real64  !< sublimation this timestep (cm)
      real(real64) :: slw      = 0.0_real64  !< liquid water in snowpack (cm)
      real(real64) :: ssnow    = 0.0_real64  !< snow storage (cm water equivalent)
      real(real64) :: snowinco = 0.0_real64  !< snow-in-canopy snapshot (cm w.e.)
      real(real64) :: graidt   = 0.0_real64  !< gross rainfall this timestep (cm)
      real(real64) :: nraidt   = 0.0_real64  !< net rainfall this timestep (cm)
      real(real64) :: aintcdt  = 0.0_real64  !< actual interception this timestep (cm)

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

      ! -----------------------------------------------------------------------
      ! [SS-BMI2] Snow output row buffer (SnowOutput stream).
      ! Prefix snow_output_ used to namespace for future atmosphere output streams.
      ! N = 7: t1900, daycum, snrai, gsnow, ssnow, melt, subl
      ! -----------------------------------------------------------------------
      real(c_double),    allocatable :: snow_output_row(:)
      character(len=32), allocatable :: snow_output_columns(:)
      integer                        :: snow_output_n_cols = 0

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

   contains
      procedure :: init => atmosphere_state_init
   end type atmosphere_state_t

contains

   !> Zero all 22 flat scalars and both cohort sub-records.
   !! Type-bound init — call as state%atmosphere%init(config%meteo).
   !! meteo_cfg arg reserved for future config-driven seed migration;
   !! all seeding currently done by the dual-write block in swap_mod.f90.
   !! Mirrors heat_state%init (GR-BH Task 11).
   subroutine atmosphere_state_init(self, meteo_cfg)
      use meteorology_config_mod, only: meteorology_config_t
      class(atmosphere_state_t),  intent(inout) :: self
      type(meteorology_config_t), intent(in)    :: meteo_cfg

      ! Instantaneous (11)
      self%peva     = 0.0_real64
      self%ptra     = 0.0_real64
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
