!> @file timecontrol_state.f90
!! Typed state record for the TimeControl subsystem (ADR 0041, Task TC-1).
!!
!! Holds 79 runtime-state fields owned by timecontrol.f90 and IterTime:
!!
!!   Group C — 16 runtime clock + calendar scalars:
!!     Real(8): t, t1900, tcum, timjan1
!!     Real(4): fsec
!!     Integer: daynr, daycum, iyear, iyearm1, imonth, daymeteo,
!!              yearmeteo, swmeteo, nextyear
!!     Integer array: datea(6)
!!     Character: date
!!
!!   Group D — 19 timestep / output-schedule scalars:
!!     Real(8): dt, dtold, dtEvent, tEvent, dtprevious, tchange,
!!              tcumold, metperiod
!!     Real(4): tmptimestart, tmptimeend
!!     Integer: flprevious, isteps, nprintcount, cntper, ioutdat,
!!              ioutdatint, rainrec, wrecord
!!     Logical: flTnext
!!
!!   Group E — 28 runtime-evaluated boolean flags:
!!     per-step / per-day: flDayStart, flDayEnd, flRunEnd, flYearStart,
!!       floutputshort, floutput, flbaloutput, flheader, flheadirg,
!!       flIrg1Start, flUpdMetDet, fldecdtmin, fldtmin, fldtreduce
!!     init-once (derived from config at TC(case=1)):
!!       flprintshort, flmetdetail, flmeteodt, flrainintens, fletsine,
!!       flSnow, flDrain, flSurfaceWater, flTemperature, flSolute,
!!       flIrrigate, flOpenFileDev
!!     cross-subsystem reset gates (migrated from variables.f90):
!!       flZeroIntr, flZeroCumu
!!
!! No timecontrol_init routine — init happens inline in TimeControl(case=1)
!! (D4 from design doc). Type aggregated under swap_state_t as
!! state%timecontrol.
!!
!! Excluded (deferred per design doc):
!!   - Group A: 18 config-constants (tstart, tend, dtmin, dtmax,
!!     nprintday, period, etc.) — already in simulation_config_t.
!!   - dtEventRain — cross-owned with meteodt (atmosphere writer).
!!
!! Field name mappings from legacy globals (tc_* prefix dropped):
!!   tc_datea(6)    -> datea      tc_dtEvent     -> dtEvent
!!   tc_fsec        -> fsec       tc_tEvent      -> tEvent
!!   tc_nextyear    -> nextyear   tc_dtprevious  -> dtprevious
!!   tc_flprevious  -> flprevious tc_tchange     -> tchange
!!   tc_flTnext     -> flTnext    tc_tcumold     -> tcumold
!!   tc_tmptimestart-> tmptimestart
!!   tc_tmptimeend  -> tmptimeend
!!
!! See ADR 0041 (pending),
!!     docs/superpowers/specs/2026-05-12-state-migration-timecontrol-design.md
!!     docs/superpowers/plans/2026-05-12-timecontrol-state-migration.md

module timecontrol_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   public :: timecontrol_state_t

   type :: timecontrol_state_t

      ! -----------------------------------------------------------------------
      ! Group C — Runtime clock + calendar state (16 fields)
      ! -----------------------------------------------------------------------

      ! Real(8) time scalars — mutated every timestep/every day
      real(real64) :: t       = 0.0_real64  !! time since start of calendar year (d)
      real(real64) :: t1900   = 0.0_real64  !! time since 1 Jan 1900 (d)
      real(real64) :: tcum    = 0.0_real64  !! time since start of simulation (d)
      real(real64) :: timjan1 = 0.0_real64  !! time of 1 Jan of current year (d)

      ! Real(4) calendar helper — legacy precision preserved
      real(4) :: fsec = 0.0  !! seconds fraction in date conversions (former tc_fsec)

      ! Integer calendar scalars — per-day or per-year
      integer :: daynr    = 0  !! day number of calendar year
      integer :: daycum   = 0  !! day number from start of simulation
      integer :: iyear    = 0  !! year number of calendar year
      integer :: iyearm1  = 0  !! year number of previous calendar year
      integer :: imonth   = 0  !! month number of calendar year
      integer :: daymeteo = 0  !! calendar day for which meteo data should be read
      integer :: yearmeteo = 0 !! year for which meteo data should be read
      integer :: swmeteo  = 0  !! switch: 1=daily meteo, 2=detailed meteo (crop)
      integer :: nextyear = 0  !! year counter — next calendar year (former tc_nextyear)

      ! Integer array — date decomposition work buffer (former tc_datea)
      integer :: datea(6) = 0  !! date array for calendar conversions (y/m/d/h/min/s)

      ! Character string — formatted current date
      character(len=11) :: date = ''  !! current date string 'yyyy-mm-dd '

      ! -----------------------------------------------------------------------
      ! Group D — Timestep / output-schedule machinery (19 fields)
      ! -----------------------------------------------------------------------

      ! Real(8) timestep scalars
      real(real64) :: dt         = 0.0_real64  !! current time step (d)
      real(real64) :: dtold      = 0.0_real64  !! previous time step (for macropore-iteration)
      real(real64) :: dtEvent    = 0.0_real64  !! max time step to next scheduled event (former tc_dtEvent)
      real(real64) :: tEvent     = 0.0_real64  !! time to next scheduled event within day (former tc_tEvent)
      real(real64) :: dtprevious = 0.0_real64  !! length of previous timestep (former tc_dtprevious)
      real(real64) :: tchange    = 0.0_real64  !! time of next output/control change (former tc_tchange)
      real(real64) :: tcumold    = 0.0_real64  !! tcum at previous output event (former tc_tcumold)
      real(real64) :: metperiod  = 0.0_real64  !! length of detailed meteo sub-period (d)

      ! Real(4) CPU-time watchdog (IterTime) — legacy precision preserved
      real(4) :: tmptimestart = 0.0  !! cpu_time start for IterTime (former tc_tmptimestart)
      real(4) :: tmptimeend   = 0.0  !! cpu_time end for IterTime (former tc_tmptimeend)

      ! Integer schedule counters
      integer :: flprevious  = 0  !! previous timestep interval status: 1=dt-bounded, 2=event-bounded (former tc_flprevious)
      integer :: isteps      = 0  !! number of timesteps since start of day
      integer :: nprintcount = 0  !! counter for sub-daily output index
      integer :: cntper      = 0  !! day number within current output period
      integer :: ioutdat     = 0  !! counter for balance output dates
      integer :: ioutdatint  = 0  !! counter for intermediate output dates
      ! [GR-TIME 2026-05-25] Output-date schedules — populated by the config
      ! adapter (populate_outdatint_monthly for swmonth=1; otherwise stay
      ! zero-allocated). Read by timecontrol_advance to trigger output dumps.
      real(real64), allocatable :: outdat(:)     !! Output dates for water/solute balances
      real(real64), allocatable :: outdatint(:)  !! Intermediate output dates
      ! [GR-IO 2026-05-25] Output file units — assigned by file_open in
      ! swapoutput.f90 and read by subsequent write() calls. Follow the
      ! state%crop%common%file_unit_crp precedent.
      integer :: file_unit_inc = 0  !! *.inc water-balance incremental output
      integer :: file_unit_rot = 0  !! *.rot root-water-extraction output
      integer :: file_unit_tem = 0  !! *.tem soil-temperature output
      integer :: file_unit_snw = 0  !! *.snw snowpack output
      integer :: rainrec     = 0  !! rain-event record index (for sub-daily rain)
      integer :: wrecord     = 0  !! detailed meteo record index within day

      ! Logical schedule flag
      logical :: flTnext = .false.  !! flag: time control needs to update tEvent (former tc_flTnext)

      ! Real(8) output-period accumulator
      real(real64) :: outper = 0.0_real64  !! length of actual output interval (d)

      ! -----------------------------------------------------------------------
      ! Group E — Runtime-evaluated boolean flags (28 fields)
      ! -----------------------------------------------------------------------

      ! Per-step flags (flipped by TimeControl each step/day)
      logical :: flDayStart     = .false.  !! first timestep of the day
      logical :: flDayEnd       = .false.  !! last timestep of the day
      logical :: flRunEnd       = .false.  !! simulation end reached
      logical :: flYearStart    = .false.  !! first day of a new calendar year

      ! Output schedule flags
      logical :: floutputshort  = .false.  !! sub-daily output time reached
      logical :: floutput       = .false.  !! daily output time reached
      logical :: flbaloutput    = .false.  !! balance output time reached
      logical :: flheader       = .false.  !! header should be printed in output file
      logical :: flheadirg      = .false.  !! header should be printed in irrigation output file
      logical :: flIrg1Start    = .false.  !! flag for 1st-crop irrigation output start

      ! Meteo / detailed-input flags
      logical :: flUpdMetDet    = .false.  !! detailed meteo record needs refresh

      ! Timestep-control flags
      logical :: fldecdt        = .false.  !! decrease-timestep signal (SurfaceWater oscillation, Richards non-convergence)
      logical :: fldecdtmin     = .false.  !! timestep should be reset to dtmin
      logical :: fldtmin        = .false.  !! current dt equals dtmin
      logical :: fldtreduce     = .false.  !! timestep reduction flagged by swap.f90

      ! Init-once flags (derived from config at TC(case=1), constant afterwards)
      logical :: flprintshort   = .false.  !! sub-daily output mode active
      logical :: flmetdetail    = .false.  !! detailed meteo input active
      logical :: flmeteodt      = .false.  !! meteodt sub-stepping active (flmetdetail or flrainintens)
      logical :: flrainintens   = .false.  !! rain-intensity input active (swrain>0)
      logical :: fletsine       = .false.  !! sine-wave ET correction active
      logical :: flSnow         = .false.  !! snow module active
      logical :: flDrain        = .false.  !! drainage module active (swdra==1)
      logical :: flSurfaceWater = .false.  !! surface water module active (swdra==2)
      logical :: flTemperature  = .false.  !! soil temperature module active
      logical :: flSolute       = .false.  !! solute module active
      logical :: flIrrigate     = .false.  !! fixed irrigation active (swirfix==1)
      logical :: flOpenFileDev  = .false.  !! developer output file open flag

      ! [SS-BMI2] time bounds + dt limits (migrated from variables.f90)
      real(real64) :: tstart      = 0.0_real64  !! simulation start (days since 1900)
      real(real64) :: tend        = 0.0_real64  !! simulation end (days since 1900)
      real(real64) :: dtmin       = 0.0_real64  !! minimum timestep (d)
      real(real64) :: dtmax       = 0.0_real64  !! maximum timestep (d)

      ! [SS-BMI2] output cadence + format switches (migrated from variables.f90)
      integer      :: period      = 0           !! output period length (d)
      integer      :: nprintday   = 0           !! number of output times per day
      logical      :: flprintdt   = .false.     !! sub-daily print active
      integer      :: swheader    = 0           !! header-print switch
      integer      :: swodat      = 0           !! output-date switch
      integer      :: swres       = 0           !! result-file switch
      integer      :: swscre      = 0           !! screen-write switch

      ! [SS-BMI2] iteration control (migrated from variables.f90)
      integer      :: MaxIt        = 0          !! max Richards iterations
      integer      :: MaxIterTime  = 0          !! max CPU seconds before abort
      integer      :: msteps       = 0          !! max timesteps per day
      logical      :: flMaxIterTime = .false.   !! iteration-time guard active

      ! [SS-BMI2] cffi runtime mode — set by swap_set_headless before initialize.
      ! When .true., output formatters fill state buffers but skip file open/write.
      logical      :: headless    = .false.     !! cffi headless mode (no CSV files)

      ! [SS-TCM] cross-subsystem reset gates (migrated from variables.f90)
      ! Owned by TimeControl; consumed by 8 physics readers. Reset
      ! cadence: flZeroIntr clears intermediate accumulators; flZeroCumu
      ! clears cumulative accumulators. See design spec
      ! 2026-05-13-timecontrol-modernization-design.md.
      logical :: flZeroIntr     = .false.  !! reset gate: intermediate accumulators
      logical :: flZeroCumu     = .false.  !! reset gate: cumulative accumulators

      ! [state%cfg retirement 2026-05-28] Config snapshots for runtime reads
      ! in timecontrol_advance. Snapshotted once at init; not changed during run.
      integer :: swmetdetail = 0  !! snapshot of config%meteo%swmetdetail
      integer :: swrain      = 0  !! snapshot of config%meteo%swrain
      integer :: swssdi      = 0  !! snapshot of config%irrigation%swssdi

   contains
      procedure :: init => timecontrol_state_init
   end type timecontrol_state_t

contains

   !> Seed runtime time-control state from typed config.
   !!
   !! Computes iyear/imonth from tstart via dtdpar() (mirrors readswap.f90:126-128).
   !! Allocates outdat/outdatint to MAOUT and zero-fills (consumed by
   !! timecontrol_advance to gate output dumps). When config_simulation%swmonth==1,
   !! populates outdatint with end-of-month dates via populate_outdatint_monthly
   !! and forces period/swres/swodat = 0 (legacy behaviour).
   subroutine timecontrol_state_init(self, config_simulation, config_general, config_drain)
      use simulation_config_mod, only: simulation_config_t
      use general_config_mod,    only: general_config_t
      use drainage_config_mod,   only: drainage_config_t
      use swap_array_dimensions, only: maout
      class(timecontrol_state_t),  intent(inout) :: self
      type(simulation_config_t),   intent(in)    :: config_simulation
      type(general_config_t),      intent(in)    :: config_general
      type(drainage_config_t),     intent(in)    :: config_drain
      integer :: datea_init(6)
      real    :: fsec_init

      ! General — screen-/result-file switches
      self%swscre = config_general%swscre

      ! Simulation — clock window + output cadence + numerical solver
      self%tstart    = config_simulation%tstart
      self%tend      = config_simulation%tend
      self%nprintday = config_simulation%nprintday
      self%period    = config_simulation%period
      self%swres     = config_simulation%swres
      self%swodat    = config_simulation%swodat

      ! Numerical (sub-record)
      self%dt    = config_simulation%numerical%dt
      self%dtmin = config_simulation%numerical%dtmin
      self%dtmax = config_simulation%numerical%dtmax
      self%MaxIt = config_simulation%numerical%MaxIt
      self%msteps = config_simulation%numerical%msteps
      ! Not in schema — defaults match legacy.
      self%MaxIterTime   = 0
      self%flMaxIterTime = .false.

      ! Derive iyear/imonth from tstart (mirrors readswap.f90:126-128).
      call dtdpar(self%tstart + 0.1d0, datea_init, fsec_init)
      self%iyear  = datea_init(1)
      self%imonth = datea_init(2)

      ! Output-date schedules — allocate to legacy cap, zero-fill.
      if (.not. allocated(self%outdat))    allocate(self%outdat(maout))
      if (.not. allocated(self%outdatint)) allocate(self%outdatint(maout))
      self%outdat    = 0.0_real64
      self%outdatint = 0.0_real64

      ! Monthly output: populate end-of-month dates + clobber daily-period switches.
      if (config_simulation%swmonth == 1) then
         call populate_outdatint_monthly(self%tend, self%iyear, self%imonth, &
                                         self%outdatint)
         self%period = 0
         self%swres  = 0
         self%swodat = 0
      end if

      ! Output switches: legacy forced to 0 by ADR 0009.
      self%swheader = 0

      ! Drainage/surface-water flags derived from config_drain%swdra.
      ! Set here so timecontrol_init (timecontrol_mod.f90) reads the typed
      ! value directly rather than going through state%cfg%drain%swdra.
      self%flDrain        = (config_drain%swdra == 1)
      self%flSurfaceWater = (config_drain%swdra == 2)
   end subroutine timecontrol_state_init

   !> Populate `outdatint(:)` with the end-of-month dates between
   !! `tstart` and `tend`. Mirrors `readswap.f90:181-204` (the
   !! `swmonth == 1` branch). Relocated from config_to_variables_mod
   !! as a module-private helper [GR-SEED 2026-05-25].
   subroutine populate_outdatint_monthly(tend, iyear, imonth, outdatint)
      real(real64), intent(in)    :: tend
      integer,      intent(in)    :: iyear   !! start year
      integer,      intent(in)    :: imonth  !! start month
      real(real64), intent(inout) :: outdatint(:)
      integer  :: datea_om(6), i_om
      real     :: fsec_om
      real(real64) :: outdate_om

      datea_om = 0
      datea_om(1) = iyear
      datea_om(2) = imonth
      if (datea_om(2) < 12) then
         datea_om(2) = datea_om(2) + 1
      else
         datea_om(1) = datea_om(1) + 1
         datea_om(2) = 1
      end if
      datea_om(3) = 1
      fsec_om = 0.0
      call dtardp(datea_om, fsec_om, outdate_om)
      i_om = 0
      do while ((outdate_om - 1.0d0) < (tend + 0.1d0))
         i_om = i_om + 1
         outdatint(i_om) = outdate_om - 1.0d0
         if (datea_om(2) < 12) then
            datea_om(2) = datea_om(2) + 1
         else
            datea_om(1) = datea_om(1) + 1
            datea_om(2) = 1
         end if
         call dtardp(datea_om, fsec_om, outdate_om)
      end do
   end subroutine populate_outdatint_monthly

end module timecontrol_state_mod
