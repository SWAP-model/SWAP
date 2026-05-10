module swap_exchange

   type :: swap_input

      !sequence

      ! start and finish time; typically a single day is to be considered, which means Tend = Tstart!!!
      ! units: days since 1900-01-01
      real(8) :: tstart, tend

      ! weather data
      ! units:   oC    oC    kPa  m/s   mm/d  -    mm/d   kJ/m2/d
      real(8) :: tmin, tmax, hum, wind, rain, wet, etref, rad

      ! crop variables from external crop
      ! cropheight, rooting depth and lai of external crop; overwrites value of dummy crop
      ! units    cm  cm     m2/m2
      real(8) :: ch, zroot, lai

      ! yes (1) or no (0) crop present
      integer(8) :: icrop

   end type

   type :: swap_output

      !sequence

      real(8) :: tstart, tend

      ! actual number of nodes of whole oil profile; also actual length of return arrays
      integer(8)                          :: numnodes

      ! Potential and actual daily transpiration
      real(8)                             :: tpot, tact

      ! for error handling
      integer(8)                          :: ierrorcode

      ! return arrays: thickness of all soil layers, and per layer: volumetric water content, root water uptake
      real(8), dimension(500)             :: dz
      real(8), dimension(500)             :: wc
      real(8), dimension(500)             :: rwu

   end type

end module swap_exchange

! -----------------------------------------------------------------------------------------------------------------------
subroutine swap(iCaller, iTask, toswap, fromswap)

! The model swap can perform three major tasks (iTask):
!    1 - initialization
!    2 - dynamic (time loop)
!    3 - closure
! The input variable iCaller determines who is calling the swap model:
!    0 - the swap model is called from swap_main
!    i - (i /= 0) the swap model as DLL is called from elsewhere, and additional actions are performed regarding exhange of data

! -----------------------------------------------------------------------------------------------------------------------

! This is needed in case a swap.dll is to be constructed
!dec$ attributes dllexport :: SWAP

!     swap modules for data communication
use variables, only : flyearstart, fldaystart, flswapshared, flsurfacewater, flmacropore, fltemperature, flsnow,        &
                      flsolute, flcropnut, flirrigate, flagetracer, flrunend, flmeteodt, fletsine, swfrost, fldtreduce, &
                      swusecn, fldrain, fldecmprat, fldayend, flcropcalendar, flmaxitertime, floutput,         &
                      floutputshort, flharvestday, flcropoutput, swcrp, flirrigationoutput, swend, project, &
                      flTillage, flSSDI, &
                      daynr, iyear, numnod, numlay
use timestep_control_mod, only: fldecdt
use swap_state_mod, only: swap_state_t
use drainage_mod, only: drainage, drainage_init
use surfacewater_mod, only: SurfaceWater, surfacewater_year_reset
                      ! for debugging
!use variables, only : iqrot, iptra, cnrai, t1900, Tstart, Tend, numnod, dz, theta, dt, h, arai, rainamount, lai

use tillage_mod,   only : DoTillage
use swap_exchange
use swap_log, only: log_info
use boundbottom_mod, only: BoundBottom
use runoff_mod, only: CNmethod
use meteo_mod, only: ProcessMeteoDay
use meteo_process_mod, only: ReadMeteoDay
use snow_mod, only: snow
use meteodt_mod, only: MeteoDT
use rootextraction_mod, only: RootExtraction
use frozencond_mod, only: FrozenCond, FrozenBounds
use temperature_mod, only: Temperature, heat_init
use macropore_mod, only: MACROPORE
use macroporeoutput_mod, only: MacroPoreOutput
use solute_mod, only: solute, solute_init
use agetracer_mod, only: AgeTracer
use soilgrid_mod, only: CalcGrid, ConvertDiscrVert
use soilhydraulics_mod, only: soilwater, SoilWaterStateVar
use irrigation_mod, only: irrigation, SSDI_irrigation
use management_soil_mod, only: SoilManagement
use error_mod, only: fatalerr_collected
use swap_config_mod, only: swap_config_t
use load_swap_config_mod, only: load_swap_config
use config_to_variables_mod, only: config_to_variables
implicit none

! SS-HEAT pre-Task-8: explicit interface for CropGrowth (external non-module sub
! that now takes tsoil(:) assumed-shape arg; interface required for caller).
interface
   subroutine CropGrowth(task, tsoil)
      integer, intent(in) :: task
      real(8), intent(in) :: tsoil(:)
   end subroutine CropGrowth
end interface

! global
integer,           intent(in)              :: iCaller, iTask
type(swap_input),  intent(in),    optional :: toswap
type(swap_output), intent(out),   optional :: fromswap

! local
logical :: flError
logical :: request_smaller_dt
logical, external :: dtleap
logical, parameter :: flDailyStateSnapshot = .false.
! Phase 1 (.crp port): saved config so crop_config_global pointer remains
! valid across the iTask=1 / iTask=2 / iTask=3 call boundary.
! See crop_config_global.f90 and ADR 0016/0017.
type(swap_config_t), target, save :: config
! SS-SWST Phase 1: typed surface-water state, threaded to SurfaceWater().
! SAVE ensures the state persists across the iTask=1 / iTask=2 / iTask=3
! call boundary (same pattern as `config` above).
type(swap_state_t), save :: state

if (iCaller /= 0 .and. iTask < 3) then
   if (.not.(present(toswap)))   call fatalerr_collected ('swap', 'Argument toswap missing in DLL call.')
   if (.not.(present(fromswap))) call fatalerr_collected ('swap', 'Argument fromswap missing in DLL call.')
end if
flError = .false.

!****************************************************************************************************************************
!*****   I N I T I A L I Z A T I O N   *****
!****************************************************************************************************************************
if (iTask == 1) then

!  Initialization of all variables in Module Variables
   call Initialize

!  iteration and timing statistics
   call IterTime(1)

!  Phase 4f strangler-fig: read time-independent input via the TOML
!  pipeline + config_to_variables adapter. The legacy readswap() entry
!  point is no longer called from the runtime path (see ADR 0007); it
!  survives in src/io/readswap.f90 only as a parity-test fixture. The
!  binary expects swap.toml in the current directory; abort_if_fatal
!  terminates with a clear summary if the file is absent or fails
!  validate/finalize.
!
!  Remaining strangler-fig debt: ~10 individual HACK Phase 4f-extend
!  slots in config_to_variables.f90 for legacy globals not yet covered
!  by typed schema slots (SWREDU, RSIGNI, CFEVAPPOND, iHWCKmodel, RDS,
!  ksatexm path, etc.). Each slot is a small typed-config extension +
!  adapter wiring. The deeper follow-on (per ADR 0016) is the config-
!  passing refactor — eliminate variables-module mutation by passing
!  typed config + state to compute subs explicitly.
   block
      use error_mod, only: error_collection_t
      type(error_collection_t) :: errors
      call load_swap_config('swap.toml', config, errors)
      call config%validate(errors)
      call config%finalize(errors)
      call errors%abort_if_fatal()
      call config_to_variables(config)
   end block

!  shared simulation
   if (flSwapShared) call SharedSimulation(1)

!  initialize time variables and switches/flags
   call TimeControl(1)

!  calculate grid parameters
   call CalcGrid()

   if (flTillage) call DoTillage(1)
   if (flSSDI)    call SSDI_irrigation(1)

!  initialize SoilWater rate/state variables
   call SoilWater(1, state)
   if (swuseCN == 1) call CNmethod(1)

!  Allocate and initialise drainage state arrays.  Config is passed so
!  drainage_init can seed state%drainage%wetper(1) from config%drain%wetper
!  (dramet==2) without reading the now-deleted legacy global wetper.
!  ADR 0031 Phase 2 Task 5: drainl/wetper/ztopdislay/qdrd globals deleted.
   call drainage_init(state, config)
   if (flSolute) call solute_init(state)   ! SS-SLST Phase 2 Task 7: seed state%solute from config-populated globals
   call heat_init(state)                  ! SS-HEAT Phase 1 Task 3: allocate state%heat per-node arrays

!  initialize SurfaceWater management variables
   if (flSurfaceWater) call SurfaceWater(1, state, request_smaller_dt)

!  initialize MacroPore rate/state variables
   if (flMacroPore) call MACROPORE(1, state)

!  initialize SoilTemperature rate/state variables
   if (flTemperature) call Temperature(1, state)

!  initialize Snow rate/state variables
   ! SS-HEAT Phase 2 Task 6: pass state so Snow reads tsoil from state%heat
   if (flSnow) call Snow(1, state)

!  initialize Solute rate/state variables
   if (flSolute) call Solute(1, state)

!  initialize Ageing rate/state variables
   if (flAgeTracer) call AgeTracer(1, state)

!  Soil Management init: SoilManagement(1) was the legacy reader entry
!  point and is now a no-op (SS-C step 3). flCropNut is now driven by
!  the per-rotation typed config (ADR 0028); the SoilManagement(2..7)
!  call sites below run when a rotation has flcropnut=true.

!  open Output files and write headers (skip in external/DLL mode to avoid per-column I/O)
   if (iCaller == 0) then
      call SwapOutput(1, state)
      call SoilWaterOutput(1, state)
      if (flIrrigate)     call IrrigationOutput(1)
      if (flTemperature)  call TemperatureOutput(1, state)
      if (flSolute)       call SoluteOutput(1, state)
      if (flAgeTracer)    call AgeTracerOutput(1, state)
      if (flSnow)         call SnowOutput(1)
      if (flMacroPore)    call MacroPoreOutput(1)
      if (flSurfaceWater) call SurfaceWaterOutput(1, state)
   end if

!  Specific for exchange when called as DLL
   if (iCaller /= 0) call handle_exchange(11, flError)

   call log_info('swap', 'Initialization complete for project: ' // trim(project))

   return
end if

!****************************************************************************************************************************
!*****   D Y N A M I C   *****
!****************************************************************************************************************************
if (iTask == 2) then

!  Specific for exchange when called as DLL
   if (iCaller /= 0) call handle_exchange(21, flError); if (flError) return

!  loop with soil water time step during entire simulation period
   do while (.not.flrunend)

!     get Meteo data (skip meteo-file I/O in external/DLL mode)
   if (iCaller == 0) then
      if (flYearStart) call ReadMeteoYear
   end if

      if (flDayStart) then

!        Specific for exchange when called as DLL
         if (iCaller /= 0) call handle_exchange(22, flError)   ! weather
         !if (iCaller /= 0) call handle_exchange(23, flError)   ! LAI, RD

!        read meteo data for current day
         call ReadMeteoDay()

!        check growing season
         call CropGrowth(1, state%heat%tsoil)

!        Specific for exchange when called as DLL
         if (iCaller /= 0) call handle_exchange(23, flError)   ! LAI, RD

!        calculate Irrigation rate/state variables
         if (flIrrigate) call irrigation(2, state)

!        process Meteo data
         call ProcessMeteoDay()
         if (flTillage) call DoTillage(2)

      end if

!     process Meteo data
      if (flMeteoDt .or. flETSine) call MeteoDT()

!     shared simulation
      if (flSwapShared .and. flDayStart) call SharedSimulation(2)

!     calculate Snow: MH+MM - probably to be moved within IF-block above, prior to call ProcessMeteoDay ...
      ! SS-HEAT Phase 2 Task 6: pass state so Snow reads tsoil from state%heat
      if (flSnow .and. flDayStart) call Snow(2, state)

!     calculate reduction for conductivities for frozen conditions
      if (SwFrost.eq.1) then
         call FrozenCond(state)
      end if

!     calculate potential and actual root water extraction profile
      call RootExtraction(state)

!     determine SoilWater bottom boundary conditions
      call BoundBottom()

      fldtreduce = .true.
      do while(fldtreduce)
         fldtreduce = .false.

!        calculate drainage fluxes
         if (fldrain)                           call Drainage(state)
         ! SS-SWST Phase 2: SurfaceWater sets request_smaller_dt; propagate to fldecdt here.
         if (.not.fldecdt .and. flSurfaceWater) call SurfaceWater(2, state, request_smaller_dt)
         if (request_smaller_dt) fldecdt = .true.
         if (SwFrost.eq.1)                      call FrozenBounds(state)

!        calculate SoilWater, incl macropores (headcalc inside may also set fldecdt on non-convergence)
         if (.not.fldecdt) call SoilWater(2, state)

!        calculate surface water balance
         if (.not.fldecdt .and. flSurfaceWater) call SurfaceWater(3, state, request_smaller_dt)
         if (request_smaller_dt) fldecdt = .true.

!        update time variables and switches/flags
         if (fldecdt .or. (flMacroPore .and. FlDecMpRat))then
            call SoilWaterStateVar(2)
            call TimeControl(3)
            fldtreduce = .true.
         end if

      end do

!     calculate SoilWater rate/state variables
      call SoilWater(3, state)

!     calculate SoilTemperature rate/state variables
   if (flTemperature) call Temperature(2, state)

!     calculate Solute rate/state variables
      if (flSolute) call Solute(2, state)

!     calculate Ageing rate/state variables
      if (flAgeTracer) call AgeTracer(2, state)

!     update time variables and switches/flags
      call TimeControl(2)

!     at the end of a day,
      if (flDayEnd) then

!        update Soil nutrient status variables
         if (flCropNut) call SoilManagement(2, state)

!        calculate potential crop growth
!        this is skipped in case called externally
         if (iCaller == 0 .and. flCropCalendar) call CropGrowth(2, state%heat%tsoil)

!        amendent of crop residues from previous day
         if (flCropNut) call SoilManagement(5, state)

!        amendent of fertilizers of current day
         if (flCropNut) call SoilManagement(3, state)

!        calculate actual crop growth (calculation of actual crop rate and state variables)
!        this is skipped in case called externally, so that LAI and CF remain their input values (for printing)
         if (iCaller == 0  .and. flCropCalendar) call CropGrowth(3, state%heat%tsoil)

!        Simulate Soil Nutrient processes
         if (flCropNut) call SoilManagement(4, state)

!        harvest of crop
!        this is skipped in case called externally, so that LAI and CF remain their input values (for printing)
         if (iCaller == 0 .and. flCropCalendar) call CropGrowth(4, state%heat%tsoil)

!        timing statistics : prevent (near) endless simulations
         if (flMaxIterTime) call IterTime(2)

!        Better here: check if subsurface irrigation is required for next day,
!                     and determine if time step needs to be changed due to dt_SSDI_event
         if (flSSDI) call SSDI_irrigation(2)
         call TimeControl(9)

      end if

!     output section (skip in external/DLL mode to avoid per-column I/O)
      if (iCaller == 0) then
         if (flOutput) then
            call SwapOutput(2, state)
            call SoilWaterOutput(2, state)
            if (flTillage) call DoTillage(3)
            if (flTemperature)   call TemperatureOutput(2, state)
            if (flSolute)        call SoluteOutput(2, state)
            if (flAgeTracer)     call AgeTracerOutput(2, state)
            if (flSnow)          call SnowOutput(2)
            if (flMacroPore)     call MacroPoreOutput(2)
            if (flSurfaceWater) then
               if (daynr == merge(366, 365, dtleap(iyear))) &
                  call surfacewater_year_reset(state%surfacewater)
               call SurfaceWaterOutput(2, state)
            end if
         else
            if (flOutputShort)   call SoilWaterOutput(2, state)
         end if
         if (flDayEnd .and. (flOutput .or. flHarvestDay)) then
            if (flCropCalendar .and. flCropOutput) then
               if (swcrp.eq.1) call CropOutput(2)
            end if
         end if
         if (flIrrigationOutput)          call IrrigationOutput(2)
         if (flDayEnd .and. flCropNut)    call SoilManagement(6, state)
         if (swend.eq.2 .and. flDayEnd)   call soilwateroutput(3, state)
      end if

!    shared simulation
     if (flSwapShared .and. flDayEnd) call SharedSimulation(3)

   end do

!  Specific for exchange when called as DLL
   if (iCaller /= 0) call handle_exchange(29, flError)

   return
end if

!****************************************************************************************************************************
!*****   C L O S U R E   *****
!****************************************************************************************************************************
if (iTask == 3) then

!  iteration and timing statistics
   call IterTime(3)

!  close output files (skip in external/DLL mode)
   if (iCaller == 0) then
      if (flSwapShared) call SharedSimulation(4)
      call SwapOutput(3, state)
      if (swend.eq.1) call SoilWaterOutput(3, state)
      call SoilWaterOutput(4, state)
      if (swcrp.eq.1) call CropOutput(3)
      if (flTemperature)        call TemperatureOutput(3, state)
      if (flSolute)             call SoluteOutput(3, state)
      if (flAgeTracer)          call AgeTracerOutput(3, state)
      if (flIrrigate)           call IrrigationOutput(3)
      if (flSnow)               call SnowOutput(3)
      if (flMacroPore)          call MacroPoreOutput(3)
      if (flSurfaceWater)       call SurfaceWaterOutput(3, state)
      if (flCropNut)            call SoilManagement(7, state)
   end if

!  write okay file for external use
   call WriteSwapOk(Project)

!  Specific for exchange when called as DLL
   if (iCaller /= 0) call handle_exchange(31, flError)

   call log_info('swap', 'Simulation complete for project: ' // trim(project))

   return
end if

contains

!  routine to handle exchange with calling program
   subroutine handle_exchange(task, flError)
   use variables, only : swetr, swdivide, swmetdetail, swrain, logf
   use variables, only : t1900, iyear, Tstart, Tend, numnod, dz, theta
   use variables, only : lai, ch, rd, iptra, iqrot, inqrot, flCropCalendar, flCropEmergence, flCropHarvest
   use variables, only : arad, atmn, atmx, awin, ahum, wet, arai, aetr, rainfluxarray, raintimearray   !, rainamount
   use variables, only : ex_tlast, daynrfirst, daynrlast
   implicit none
   integer, intent(in)   :: task
   logical, intent(out)  :: flError
   ! local
   integer               :: i
   integer, dimension(6) :: datea
   real                  :: fsec

! NOTE: the optional arguments in argument list of swap cannot be saved automatically with the attribute SAVE.
!       Therefore, each time allocation is needed and basic information must be set again

   fromswap%ierrorcode = 0

!  use tasks 11-19 to handle initial aspects
   if (task == 11) then
      ! some error checking
      if (swetr /= 0) then
         fromswap%ierrorcode = -1
         write (logf, '(A)') "swetr /= 0"
      end if
      if (swdivide /= 1) then
         fromswap%ierrorcode = -2
         write (logf, '(A)') "swdivide /= 1"
      end if
      if (swmetdetail /= 0) then
         fromswap%ierrorcode = -3
         write (logf, '(A)') "swmetdetail /= 0"
      end if
      if (swrain /= 0) then
         fromswap%ierrorcode = -4
         write (logf, '(A)') "swrain /= 0"
      end if
!      if (swrain /= 0 .and. swrain /= 2) then
!         fromswap%ierrorcode = -4
!         write (logf, '(A)') "swrain /= 0 .and. swrain /= 0"
!      end if

      fromswap%tstart     = Tstart
      fromswap%tend       = Tend
      fromswap%tpot       = iptra
      fromswap%tact       = iqrot
      fromswap%numnodes   = numnod
      !allocate(fromswap%dz(numnod));  fromswap%dz(1:numnod)  = dz(1:numnod)
      !allocate(fromswap%wc(numnod));  fromswap%wc(1:numnod)  = theta(1:numnod)
      !allocate(fromswap%rwu(numnod)); fromswap%rwu(1:numnod) = inqrot(1:numnod)

      fromswap%dz(1:numnod)  = dz(1:numnod)
      fromswap%wc(1:numnod)  = theta(1:numnod)
      fromswap%rwu(1:numnod) = inqrot(1:numnod)
      ex_tlast = 0.0d0
   end if

!  use tasks 21-29 to handle dynamic aspects
   if (task == 21) then
      Tstart = toswap%tstart
      Tend   = toswap%tend

      ! check
      if (ex_tlast > 0.0d0 .and. dabs(Tstart - ex_tlast) > 1.0d-8) then
         fromswap%ierrorcode = 1
         write (logf, '(A)') 'Unexpected timing error: tstart /= tlast'
      end if
      if (dabs(Tend - Tstart) > 1.0d-8) then
         fromswap%ierrorcode = 2
         write (logf, '(A)') 'Only single day allowed: tend must equal tstart'
      end if

      ! need to re-initialize
      flrunend   = .false.
      flDayStart = .true.

      ! first set iyear for proper use in TimeControl; this allows for start any time, irrespective of tstart in swap.swp
      call dtdpar (Tstart, datea, fsec)
      iyear = datea(1)
      call TimeControl(1)

      ! External forcing mode: provide full-year availability without reading meteo files
      daynrfirst = 1
      daynrlast  = 366

   end if

   if (task == 22) then
      arad(1:366) = toswap%rad*1000.0d0               ! Convert radiation from kJ/m2/d to J/m2/d
      atmn(1:366) = toswap%tmin                       ! deg. C
      atmx(1:366) = toswap%tmax                       ! deg. C
      ahum(1:366) = toswap%hum                        ! kPa
      awin(1:366) = toswap%wind                       ! m/s
      arai(1:366) = toswap%rain                       ! mm/d
      aetr(1:366) = toswap%etref                      ! mm/d
      ! in case swrain = 2                                                       !!!!! THIS IS NOT YET WORKING PROPERLY   !!!!!
      if (swrain == 2) then
         do i = 1, 366, 2
            raintimearray(2*i-1) = dble(i-1)
            raintimearray(2*i)   = dble(i-1) + wet(i)
            rainfluxarray(2*1-1) = 0.0d0
            if (wet(i) > 0.0d0) then
               rainfluxarray(2*1)   = toswap%rain*0.1d0/wet(i)     ! from mm/d to cm/d
            else
               rainfluxarray(2*1)   = 0.0d0     ! no check on consistency that both rain and wet should be either both > 0 or both = 0
            end if
         end do
         !rainfluxarray(1:366) = toswap%rain*0.1d0     ! from mm/d to cm/d
         !rainamount(1:366)    = toswap%rain          ! mm/d
         wet(1:366)           = toswap%wet            ! [0...1]
      end if
   end if
   if (task == 23) then
      lai = toswap%lai                        ! m2/m2
      ch  = toswap%ch                         ! cm
      rd  = toswap%zroot                      ! cm

      ! set crop status
      flCropCalendar  = toswap%icrop /= 0
      flCropEmergence = toswap%icrop /= 0
      flCropHarvest   = toswap%icrop == 0

   end if

   if (task == 29) then
      fromswap%numnodes      = numnod
      fromswap%tpot          = iptra
      fromswap%tact          = iqrot
      !if(.not.allocated(fromswap%dz))  allocate(fromswap%dz(numnod));  fromswap%dz(1:numnod)  = dz(1:numnod)
      !if(.not.allocated(fromswap%wc))  allocate(fromswap%wc(numnod));  fromswap%wc(1:numnod)  = theta(1:numnod)
      !if(.not.allocated(fromswap%rwu)) allocate(fromswap%rwu(numnod)); fromswap%rwu(1:numnod) = 0.0d0
      fromswap%dz(1:numnod)  = dz(1:numnod)
      fromswap%wc(1:numnod)  = theta(1:numnod)
      fromswap%rwu(1:numnod) = 0.0d0
      ex_tlast = t1900
   end if

!  use tasks 31-39 to handle closure aspects
   if (task == 31) then
   end if

!  set return flerror
   flError = fromswap%iErrorCode /= 0

   end subroutine handle_exchange

end subroutine swap
