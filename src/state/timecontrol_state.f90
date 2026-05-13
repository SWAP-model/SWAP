!> @file timecontrol_state.f90
!! Typed state record for the TimeControl subsystem (ADR 0041, Task TC-1).
!!
!! Holds 63 runtime-state fields owned by timecontrol.f90 and IterTime:
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

      ! [SS-TCM] cross-subsystem reset gates (migrated from variables.f90)
      ! Owned by TimeControl; consumed by 8 physics readers. Reset
      ! cadence: flZeroIntr clears intermediate accumulators; flZeroCumu
      ! clears cumulative accumulators. See design spec
      ! 2026-05-13-timecontrol-modernization-design.md.
      logical :: flZeroIntr     = .false.  !! reset gate: intermediate accumulators
      logical :: flZeroCumu     = .false.  !! reset gate: cumulative accumulators

   end type timecontrol_state_t

end module timecontrol_state_mod
