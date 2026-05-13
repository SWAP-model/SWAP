!> @file swap_mod.f90
!! SS-DRV Phase 1: module form of the legacy `subroutine swap`.
!! Three named lifecycle procedures replace the (iCaller, iTask) dispatch.
!! Time loop lives in the caller. State and config are threaded explicitly.
!!
!! Note: heavy compute `use` statements are scoped to each procedure rather
!! than the module header. The unit-test build only provides state/config
!! modules, and hoisting all uses would require the full compute module chain
!! at compile time even for the lightweight smoke tests.
module swap_mod
   use swap_state_mod,  only: swap_state_t
   use swap_config_mod, only: swap_config_t
   implicit none
   private
   public :: swap_init, swap_run_step, swap_close, swap_init_from_loaded_config

contains

   subroutine swap_init(config_file, state, config)
      use load_swap_config_mod, only: load_swap_config
      character(len=*),            intent(in)  :: config_file
      type(swap_state_t),          intent(out) :: state
      type(swap_config_t), target, intent(out) :: config  ! target: crop_config_global => config%crop (set inside body, Task 3)

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
      call load_swap_config(config_file, config, errors)
      call config%validate(errors)
      call config%finalize(errors)
      call errors%abort_if_fatal()
   end block

   call swap_init_from_loaded_config(state, config)

   end subroutine swap_init

   subroutine swap_init_from_loaded_config(state, config)
      use variables, only : flswapshared, flcropnut, flagetracer, swfrost, &
                            swusecn, flcropcalendar, &
                            flharvestday, flcropoutput, swcrp, flirrigationoutput, swend, project, &
                            flTillage, flSSDI, &
                            numnod, numlay, &
                            dz, z, disnod, ztopcp, zbotcp, layer
      use soilwater_state_mod, only: soilwater_init
      use atmosphere_state_mod, only: atmosphere_init
      use tillage_state_mod, only: tillage_init
      use drainage_mod, only: drainage_init
      use surfacewater_mod, only: SurfaceWater
      use tillage_mod,   only : DoTillage
      use swap_log, only: log_info
      use runoff_mod, only: CNmethod
      use snow_mod, only: snow
      use temperature_mod, only: Temperature, heat_init
      use solute_mod, only: solute, solute_init
      use agetracer_mod, only: AgeTracer
      use soilgrid_mod, only: CalcGrid
      use soilhydraulics_mod, only: soilwater
      use config_to_variables_mod, only: h_init_buf, pondini_init_buf, pond_init_buf, &
                                         tc_iyear_init_buf, tc_imonth_init_buf, tc_dt_init_buf, &
                                         config_to_variables
      use irrigation_mod, only: SSDI_irrigation
      use timecontrol_mod, only: timecontrol_init, itertime_init
      type(swap_state_t),          intent(out)   :: state
      type(swap_config_t), target, intent(inout) :: config  ! target: crop_config_global => config%crop (Task 3); inout: already loaded by caller
      logical :: request_smaller_dt   ! intent(out) dummy for SurfaceWater(1)

!  Initialization of all variables in Module Variables
   call Initialize

!  iteration and timing statistics
   call itertime_init(state)

!  config_to_variables seeds state%timecontrol from config (Task 5).
   call config_to_variables(config, state)

   ! [SS-TC TC-14] seed state%timecontrol from transient buffers before
   ! TimeControl(1) consumes them (iyear/imonth derived from tstart, dt from config).
   state%timecontrol%iyear  = tc_iyear_init_buf
   state%timecontrol%imonth = tc_imonth_init_buf
   state%timecontrol%dt     = tc_dt_init_buf

!  shared simulation
   if (flSwapShared) call SharedSimulation(1)

!  initialize time variables and switches/flags
   call timecontrol_init(state)

!  calculate grid parameters
   call CalcGrid()
   ! [SS-GR-BH A3] dual-write: populate state%mesh alongside legacy mesh globals
   call state%mesh%init(numnod, dz, z, disnod, ztopcp, zbotcp, layer)
   call soilwater_init(state%soilwater, numnod, numlay)   ! SS-CRP Phase 1 C-1.2: allocate per-node arrays + mfluxtable
   ! [SS-SWC S-2.12B] seed state%soilwater from config buffers populated by config_to_variables.
   ! Retired globals: pondini, pond, h(1..nhead). state%soilwater%h is sized numnod and
   ! receives the swinco=3 initial-profile h values; SoilHydraulics(1) consumes the rest.
   state%soilwater%pondini = pondini_init_buf
   state%soilwater%pond    = pond_init_buf
   if (allocated(h_init_buf)) then
      block
         integer :: ki
         do ki = 1, min(size(h_init_buf), size(state%soilwater%h))
            state%soilwater%h(ki) = h_init_buf(ki)
         end do
      end block
      deallocate(h_init_buf)
   end if
   call atmosphere_init(state%atmosphere)                 ! SS-ATM Phase 1 A-1.2: zero all 22 flat scalars + cohort sub-records
   ! [SS-ATM A-2.6] swinco=3 warm-restart: seed state%atmosphere directly from config (legacy globals retired)
   if (config%soil%swinco == 3) then
      if (allocated(config%soil%initial%h_file) .and. &
          len_trim(config%soil%initial%h_file) > 0) then
         state%atmosphere%ssnow = config%soil%initial%ssnow
         state%atmosphere%ldwet = config%soil%initial%ldwet
         state%atmosphere%slw   = config%soil%initial%slw
         ! spev/saev not in config; remain zero from atmosphere_init (evaporation counters reset on rain)
         if (config%meteo%snow%swsnow /= 1) state%atmosphere%ssnow = 0.0d0
      end if
   end if

   ! [SS-TC TC-14] alias TC fields used in init block
   block
   associate( &
      flSnow         => state%timecontrol%flSnow,         &
      flSolute       => state%timecontrol%flSolute,       &
      flSurfaceWater => state%timecontrol%flSurfaceWater, &
      flTemperature  => state%timecontrol%flTemperature )

   call tillage_init(state%tillage, numlay)         ! SS-TIL T-2: allocate/zero tillage state unconditionally
   if (flTillage) call DoTillage(1, state)
   if (flSSDI)    call SSDI_irrigation(1, state)  ! [SS-SWC S-2.12B]

!  Allocate and initialise heat state arrays before SoilWater(1) so that
!  hconduc can read state%heat%tsoil(node) during hydraulic-conductivity init.
!  SS-SWC S-2.2: heat_init moved earlier to satisfy mandatory tsoil_node arg.
   call heat_init(state)                  ! SS-HEAT Phase 1 Task 3: allocate state%heat per-node arrays

!  initialize SoilWater rate/state variables
   call SoilWater(1, state)
   ! SS-ATM A-2.6: state added — CNmethod signature updated for retired nraidt/melt
   if (swuseCN == 1) call CNmethod(1, state)

!  Allocate and initialise drainage state arrays.  Config is passed so
!  drainage_init can seed state%drainage%wetper(1) from config%drain%wetper
!  (dramet==2) without reading the now-deleted legacy global wetper.
!  ADR 0031 Phase 2 Task 5: drainl/wetper/ztopdislay/qdrd globals deleted.
   call drainage_init(state, config)
   if (flSolute) call solute_init(state)   ! SS-SLST Phase 2 Task 7: seed state%solute from config-populated globals

!  initialize SurfaceWater management variables
   ! S4 state init for surfacewater (spec docs/superpowers/specs/2026-05-13-state-init-pilot-surfacewater-design.md).
   if (flSurfaceWater) call state%surfacewater%init(config%surface_water, config%drain, numnod)
   ! SurfaceWater(task=1)'s case(1) is now a no-op stub; init was hoisted to the line above.
   ! Dispatcher-case removal is a separate cleanup follow-up.
   if (flSurfaceWater) call SurfaceWater(1, state, request_smaller_dt)

!  [MACRO-RETIRE 2026-05-12] MACROPORE call retired (ADR 0040). flMacroPore is permanent .false.
!  if (flMacroPore) call MACROPORE(1, state)

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

!  open Output files and write headers (always run; iCaller branch retired)
   call SwapOutput(1, state)
   call SoilWaterOutput(1, state)
!  ADR 0009 Phase 5+: IrrigationOutput deleted (swirg=0).
   if (flTemperature)  call TemperatureOutput(1, state)
   if (flSolute)       call SoluteOutput(1, state)
   if (flAgeTracer)    call AgeTracerOutput(1, state)
   if (flSnow)         call SnowOutput(1, state)
   ! [MACRO-RETIRE 2026-05-12] MacroPoreOutput retired (ADR 0040).
   if (flSurfaceWater) call SurfaceWaterOutput(1, state)
   end associate
   end block

   call log_info('swap', 'Initialization complete for project: ' // trim(project))

   end subroutine swap_init_from_loaded_config

   subroutine swap_run_step(state, config)
      use variables, only : flswapshared, flcropnut, flagetracer, swfrost, &
                            flcropcalendar, &
                            flharvestday, flcropoutput, swcrp, swend, &
                            flTillage, flSSDI
      use timestep_control_mod, only: fldecdt
      use timecontrol_mod, only: timecontrol_advance, timecontrol_reduce_dt, &
                                  timecontrol_day_end, itertime_check
      use surfacewater_mod, only: SurfaceWater, surfacewater_year_reset
      use tillage_mod, only: DoTillage
      use boundbottom_mod, only: BoundBottom
      use meteo_mod, only: ProcessMeteoDay
      use meteo_process_mod, only: ReadMeteoDay
      use snow_mod, only: snow
      use meteodt_mod, only: MeteoDT
      use rootextraction_mod, only: RootExtraction
      use frozencond_mod, only: FrozenCond, FrozenBounds
      use temperature_mod, only: Temperature
      use solute_mod, only: solute
      use agetracer_mod, only: AgeTracer
      use soilhydraulics_mod, only: soilwater, SoilWaterStateVar
      use irrigation_mod, only: irrigation, SSDI_irrigation
      use management_soil_mod, only: SoilManagement
      use drainage_mod, only: drainage
      type(swap_state_t),  intent(inout) :: state
      type(swap_config_t), intent(in)    :: config
      logical :: request_smaller_dt
      logical, external :: dtleap

      interface
         subroutine CropGrowth(task, tsoil, state)
            use swap_state_mod, only: swap_state_t
            integer, intent(in) :: task
            real(8), intent(in) :: tsoil(:)
            type(swap_state_t), intent(inout) :: state
         end subroutine CropGrowth
      end interface

!  [SS-TC TC-14] bind TC aliases for all timestep-loop fields used below
   associate( &
      tc_flYearStart => state%timecontrol%flYearStart, &
      tc_flDayStart  => state%timecontrol%flDayStart,  &
      tc_flDayEnd    => state%timecontrol%flDayEnd,    &
      tc_daynr       => state%timecontrol%daynr,       &
      tc_iyear       => state%timecontrol%iyear,       &
      flrunend       => state%timecontrol%flRunEnd,    &
      flOutput       => state%timecontrol%floutput,    &
      flOutputShort  => state%timecontrol%floutputshort,&
      flMeteoDt      => state%timecontrol%flmeteodt,   &
      flETSine       => state%timecontrol%fletsine,    &
      flSnow         => state%timecontrol%flSnow,      &
      flSolute       => state%timecontrol%flSolute,    &
      flTemperature  => state%timecontrol%flTemperature,&
      flDrain        => state%timecontrol%flDrain,     &
      flSurfaceWater => state%timecontrol%flSurfaceWater,&
      flIrrigate     => state%timecontrol%flIrrigate,  &
      fldtreduce     => state%timecontrol%fldtreduce )

!     get Meteo data
   if (tc_flYearStart) call ReadMeteoYear(state)  ! SS-TC TC-13

      if (tc_flDayStart) then  ! SS-TC TC-13

!        read meteo data for current day
         call ReadMeteoDay(state)  ! SS-ATM A-1.6: state threaded for atmosphere dual-writes

!        check growing season
         call CropGrowth(1, state%heat%tsoil, state)  ! SS-CRP C-1.3: state added for dual-write

!        calculate Irrigation rate/state variables
         if (flIrrigate) call irrigation(2, state)

!        process Meteo data
         call ProcessMeteoDay(state)
         if (flTillage) call DoTillage(2, state)

      end if

!     process Meteo data
      if (flMeteoDt .or. flETSine) call MeteoDT(state)

!     shared simulation
      if (flSwapShared .and. tc_flDayStart) call SharedSimulation(2)  ! SS-TC TC-13

!     calculate Snow: MH+MM - probably to be moved within IF-block above, prior to call ProcessMeteoDay ...
      ! SS-HEAT Phase 2 Task 6: pass state so Snow reads tsoil from state%heat
      if (flSnow .and. tc_flDayStart) call Snow(2, state)  ! SS-TC TC-13

!     calculate reduction for conductivities for frozen conditions
      if (SwFrost.eq.1) then
         call FrozenCond(state)
      end if

!     calculate potential and actual root water extraction profile
      call RootExtraction(state)

!     determine SoilWater bottom boundary conditions
      call BoundBottom(state)  ! [SS-HEAT] Task 9: state passed for rfcp access

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
         if (fldecdt) then
            call SoilWaterStateVar(2, state)
            call timecontrol_reduce_dt(state)
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
      call timecontrol_advance(state)

!     at the end of a day,
      if (tc_flDayEnd) then  ! SS-TC TC-13

!        update Soil nutrient status variables
         if (flCropNut) call SoilManagement(2, state)

!        calculate potential crop growth
         if (flCropCalendar) call CropGrowth(2, state%heat%tsoil, state)

!        amendent of crop residues from previous day
         if (flCropNut) call SoilManagement(5, state)

!        amendent of fertilizers of current day
         if (flCropNut) call SoilManagement(3, state)

!        calculate actual crop growth (calculation of actual crop rate and state variables)
         if (flCropCalendar) call CropGrowth(3, state%heat%tsoil, state)

!        Simulate Soil Nutrient processes
         if (flCropNut) call SoilManagement(4, state)

!        harvest of crop
         if (flCropCalendar) call CropGrowth(4, state%heat%tsoil, state)

!        timing statistics : prevent (near) endless simulations
         if (state%timecontrol%flMaxIterTime) call itertime_check(state)  ! [SS-BMI2 Task 4]

!        Better here: check if subsurface irrigation is required for next day,
!                     and determine if time step needs to be changed due to dt_SSDI_event
         if (flSSDI) call SSDI_irrigation(2, state)  ! [SS-SWC S-2.12B]
         call timecontrol_day_end(state)

      end if

!     output section
         if (flOutput) then
            call SwapOutput(2, state)
            call SoilWaterOutput(2, state)
            if (flTillage) call DoTillage(3, state)
            if (flTemperature)   call TemperatureOutput(2, state)
            if (flSolute)        call SoluteOutput(2, state)
            if (flAgeTracer)     call AgeTracerOutput(2, state)
            if (flSnow)          call SnowOutput(2, state)
            ! [MACRO-RETIRE 2026-05-12] MacroPoreOutput retired (ADR 0040).
            if (flSurfaceWater) then
               if (tc_daynr == merge(366, 365, dtleap(tc_iyear))) &  ! SS-TC TC-13
                  call surfacewater_year_reset(state%surfacewater)
               call SurfaceWaterOutput(2, state)
            end if
         else
            if (flOutputShort)   call SoilWaterOutput(2, state)
         end if
         if (tc_flDayEnd .and. (flOutput .or. flHarvestDay)) then  ! SS-TC TC-13
            if (flCropCalendar .and. flCropOutput) then
               if (swcrp.eq.1) call CropOutput(2, state)
            end if
         end if
!        ADR 0009 Phase 5+: IrrigationOutput deleted (swirg=0).
         if (tc_flDayEnd .and. flCropNut)    call SoilManagement(6, state)   ! SS-TC TC-13
         if (swend.eq.2 .and. tc_flDayEnd)   call soilwateroutput(3, state)  ! SS-TC TC-13

!    shared simulation
     if (flSwapShared .and. tc_flDayEnd) call SharedSimulation(3)  ! SS-TC TC-13

   end associate  ! SS-TC TC-13: tc_flYearStart, tc_flDayStart, tc_flDayEnd, tc_daynr, tc_iyear

   end subroutine swap_run_step

   subroutine swap_close(state, config)
      use variables, only : flswapshared, flcropnut, flagetracer, project, swcrp, swend
      use swap_log,  only: log_info
      use management_soil_mod, only: SoilManagement
      use timecontrol_mod, only: itertime_close
      type(swap_state_t),  intent(inout) :: state
      type(swap_config_t), intent(in)    :: config  ! unused: kept for parallel signature with swap_init/swap_run_step

!  iteration and timing statistics
   call itertime_close(state)

!  close output files (always run; iCaller branch retired)
   if (flSwapShared) call SharedSimulation(4)
   call SwapOutput(3, state)
   if (swend.eq.1) call SoilWaterOutput(3, state)
   call SoilWaterOutput(4, state)
   if (swcrp.eq.1) call CropOutput(3, state)
   ! [SS-TC TC-14] flag reads via state%timecontrol
   if (state%timecontrol%flTemperature)  call TemperatureOutput(3, state)
   if (state%timecontrol%flSolute)       call SoluteOutput(3, state)
   if (flAgeTracer)                      call AgeTracerOutput(3, state)
!  ADR 0009 Phase 5+: IrrigationOutput deleted (swirg=0).
   if (state%timecontrol%flSnow)         call SnowOutput(3, state)
   ! [MACRO-RETIRE 2026-05-12] MacroPoreOutput retired (ADR 0040).
   if (state%timecontrol%flSurfaceWater) call SurfaceWaterOutput(3, state)
   if (flCropNut)                        call SoilManagement(7, state)

!  write okay file for external use
   call WriteSwapOk(Project)

   call log_info('swap', 'Simulation complete for project: ' // trim(project))

   end subroutine swap_close

end module swap_mod
