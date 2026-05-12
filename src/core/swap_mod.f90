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
   public :: swap_init, swap_run_step, swap_close

contains

   subroutine swap_init(config_file, state, config)
      use variables, only : flswapshared, flmacropore, flcropnut, flagetracer, swfrost, &
                            swusecn, fldecmprat, flcropcalendar, flmaxitertime, &
                            flharvestday, flcropoutput, swcrp, flirrigationoutput, swend, project, &
                            flTillage, flSSDI, &
                            numnod, numlay
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
      use WC_K_models_04_11, only: bind_cofgen_target
      use soilhydraulics_utils, only: bind_state_targets, bind_tc_target
      use config_to_variables_mod, only: h_init_buf, pondini_init_buf, pond_init_buf, &
                                         tc_iyear_init_buf, tc_imonth_init_buf, tc_dt_init_buf, &
                                         config_to_variables
      use load_swap_config_mod, only: load_swap_config
      use irrigation_mod, only: SSDI_irrigation
      character(len=*),            intent(in)  :: config_file
      type(swap_state_t),          intent(out) :: state
      type(swap_config_t), target, intent(out) :: config  ! target: crop_config_global => config%crop (set inside body, Task 3)
      logical :: request_smaller_dt   ! intent(out) dummy for SurfaceWater(1)

!  Initialization of all variables in Module Variables
   call Initialize

!  iteration and timing statistics
   call IterTime(1, state)

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
      call config_to_variables(config)
   end block

   ! [SS-TC TC-14] seed state%timecontrol from transient buffers before
   ! TimeControl(1) consumes them (iyear/imonth derived from tstart, dt from config).
   state%timecontrol%iyear  = tc_iyear_init_buf
   state%timecontrol%imonth = tc_imonth_init_buf
   state%timecontrol%dt     = tc_dt_init_buf

!  shared simulation
   if (flSwapShared) call SharedSimulation(1)

!  initialize time variables and switches/flags
   call TimeControl(1, state)

!  calculate grid parameters
   call CalcGrid()
   call soilwater_init(state%soilwater, numnod, numlay)   ! SS-CRP Phase 1 C-1.2: allocate per-node arrays + mfluxtable
   ! [SS-SWC S-2.12B] bind module-level pointers in utility modules to state%soilwater
   ! so legacy `cofgen` / `fluseksatexm` reads in WC_K_models_04_11 + soilhydraulics_utils
   ! resolve to the canonical state%soilwater storage (ADR 0038).
   call bind_cofgen_target(state%soilwater%cofgen)
   call bind_state_targets(state%soilwater%cofgen, state%soilwater%fluseksatexm)
   call bind_tc_target(state%timecontrol%dt)  ! SS-TC TC-12: wire tc_dt_ptr in moiscap()
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

   end subroutine swap_init

   subroutine swap_run_step(state, config)
      type(swap_state_t),  intent(inout) :: state
      type(swap_config_t), intent(in)    :: config
      ! Body filled in Task 4.
   end subroutine swap_run_step

   subroutine swap_close(state, config)
      type(swap_state_t),  intent(inout) :: state
      type(swap_config_t), intent(in)    :: config
      ! Body filled in Task 5.
   end subroutine swap_close

end module swap_mod
