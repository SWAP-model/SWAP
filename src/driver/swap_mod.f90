!> @file swap_mod.f90
!! Module form of the legacy `subroutine swap`. Three named lifecycle
!! procedures replace the (iCaller, iTask) dispatch. The time loop lives
!! in the caller. State and config are threaded explicitly.
!!
!! Note: heavy compute `use` statements are scoped to each procedure
!! rather than the module header. The unit-test build only provides
!! state/config modules, and hoisting all uses would require the full
!! compute module chain at compile time even for the lightweight smoke
!! tests.
module swap_mod
   use swap_state_mod,  only: swap_state_t
   use swap_config_mod, only: swap_config_t
   implicit none
   private
   public :: swap_init, swap_run_step, swap_close, swap_init_from_loaded_config
   ! swap_init_from_loaded_config: kept public for swap_capi_mod (loads config from
   ! a string buffer rather than a file path, so cannot route through swap_init).
   ! Body lives in the private swap_init_body helper; both wrappers delegate to it.

contains

   subroutine swap_init(config_file, state, config)
      use load_swap_config_mod, only: load_swap_config
      character(len=*),            intent(in)  :: config_file
      type(swap_state_t),          intent(out) :: state
      type(swap_config_t), target, intent(out) :: config  ! target retained for callers that take a pointer into config (crop_config_global retired)

      ! Load + validate + finalize TOML
      block
         use error_mod, only: error_collection_t
         type(error_collection_t) :: errors
         call load_swap_config(config_file, config, errors)
         call config%validate(errors)
         call config%finalize(errors)
         call errors%abort_if_fatal()
      end block

      ! Seeding + first-day compute init
      call swap_init_body(state, config)

   end subroutine swap_init

   !> Public shim retained for swap_capi_mod: the CAPI path loads config from
   !> a TOML string buffer (not a file path) and calls this directly after
   !> validate/finalize. Delegates to the shared private swap_init_body.
   subroutine swap_init_from_loaded_config(state, config)
      type(swap_state_t),          intent(out)   :: state
      type(swap_config_t), target, intent(inout) :: config
      call swap_init_body(state, config)
   end subroutine swap_init_from_loaded_config

   subroutine swap_init_body(state, config)
      !
      ! Init pipeline: TOML load (in swap_init) → Phase-1 type-bound state%X%init
      ! calls (construct) → Phase-2 free x_seed(state) procedures (derived initial
      ! state from sibling subsystems; see ADR 0043).
      !
      ! seed_state_from_config.f90 dissolved (OD Step 9, 2026-05-28): the strangler-
      ! pattern adapter held zero live writes after OD Steps 1-8 folded every
      ! cross-subsystem seeding into the respective state%X%init. The one remaining
      ! line (state%timecontrol%init call) is now inlined directly below.
      !
      ! Tracked debt (still works; not bugs):
      !   - Phase-2 seed procedures (tillage_seed, soilwater_seed,
      !     temperature_seed, solute_seed) — per ADR 0043 these are named
      !     free procs, deliberately distinct from the type-bound Phase-1
      !     init because they read sibling subsystems in a fixed order.
      !     Folding them into type-bound init is gated on untangling that
      !     init-order coupling (future work).
      !   - The legacy readswap() entry point survives in src/io/readswap.f90
      !     only as a parity-test fixture (ADR 0007); not called at runtime.
      !
      use tillage_mod,                only: tillage_seed
      use swap_log,                   only: to_str
      use runoff_mod,                 only: cn_init
      use temperature_mod,            only: temperature_seed
      use solute_mod,                 only: solute_seed
      use soilhydraulics_mod,         only: soilwater_seed
      use timecontrol_mod,            only: timecontrol_init, itertime_init
      use csv_output,                 only: csv_output_init

      type(swap_state_t),          intent(out)   :: state
      type(swap_config_t), target, intent(inout) :: config  ! target retained for callers that take a pointer into config (crop_config_global retired)

      ! Iteration / timing statistics.
      call itertime_init(state)

      ! Timecontrol seeding (formerly the last live line in seed_state_from_config).
      ! Must run BEFORE crop%init and timecontrol_init; crop%init provides croptype
      ! which timecontrol_init reads.
      call state%timecontrol%init(config%simulation, config%general, config%drain, config%output_csv)

      ! Crop seeding must run BEFORE timecontrol_init because
      ! timecontrol_init reads state%crop%common%croptype(icrop).
      call state%crop%init(config%crop, config%meteo, config%general%pathwork)

      ! Time variables and switches/flags.
      call timecontrol_init(state, config)

      ! Grid parameters (writes directly to state%mesh).
      call state%mesh%init(config%soil)

      ! Modern type-bound state init.
      call state%soilwater%init(config%soil, config%drain, config%heat, config%bottom_boundary, &
                                config%simulation, &
                                state%mesh%numnod, state%mesh%numlay, &
                                config%general%pathwork)
      call state%crop%irrigation%init(config%irrigation, state%timecontrol%tstart, &
                                      state%timecontrol%tend, state%mesh, &
                                      config%general%pathwork)

      ! Modern atmosphere init: zero flat scalars + cohorts + snapshot
      ! config-derived params + swinco==3 warm-restart (ssnow/ldwet/slw/atmin7)
      ! + snow handshake (snowinco ↔ ssnow). All folded here (OD Steps 6+7).
      call state%atmosphere%init(config)

      ! swinco==3 warm-restart: pond/pondini/dt/h-profile from soil.initial.
      ! swinco<3: pondini/pond from soil.pondini.
      ! atmosphere warm-restart fields are handled inside state%atmosphere%init above.
      ! [W4 fix 2026-05-28] h-profile no longer re-reads the CSV here;
      ! state%soilwater%h_init is already populated by soilwater_state_init.
      if (config%soil%swinco == 3 .and. &
          allocated(config%soil%initial%h_file) .and. &
          len_trim(config%soil%initial%h_file) > 0) then
         state%soilwater%pondini = config%soil%initial%pond
         state%soilwater%pond    = config%soil%initial%pond
         ! soil.initial.dt supersedes simulation%numerical%dt for swinco=3, but
         ! must be clamped to [dtmin, dtmax]. Legacy clamps a sub-dtmin restart dt
         ! (e.g. 1e-7 < dtmin 1e-6) to dtmin AND carries that clamped value as
         ! dtprevious. timecontrol_init instead took the `dt<dtmin -> sqrt(dtmin*
         ! dtmax)` branch, leaving dtprevious = 2e-4 (= sqrt) while legacy has
         ! dtprevious = dtmin. On step 0 the event-limit branch sets dt =
         ! dtprevious, so the 2e-4 spuriously triggered an event-limit in the
         ! modern build, desyncing the adaptive-dt sequence from 4.2.0 (GWL drift
         ! ~3 cm/4yr, solute ~3%). Clamp BOTH dt and dtprevious to match legacy.
         state%timecontrol%dt = max(min(config%soil%initial%dt, &
                                        state%timecontrol%dtmax), &
                                    state%timecontrol%dtmin)
         state%timecontrol%dtprevious = state%timecontrol%dt
         ! [W4 fix 2026-05-28] No duplicate CSV read: h_profile already loaded into
         ! state%soilwater%h_init by soilwater_state_init (Piece B). Copy h values here.
         block
            integer :: ki, nrows
            if (state%soilwater%h_init%is_loaded) then
               nrows = size(state%soilwater%h_init%rows)
               do ki = 1, min(nrows, size(state%soilwater%h))
                  state%soilwater%h(ki) = state%soilwater%h_init%rows(ki)%h
               end do
            end if
         end block
      else
         state%soilwater%pondini = config%soil%pondini
         state%soilwater%pond    = config%soil%pondini   ! legacy alias for swinco<3
      end if

      associate (time => state%timecontrol)

         ! Modern type-bound init. The swtill gate is inside
         ! state%tillage%init (skips Group AB if events not allocated).
         call state%tillage%init(config%soil%tillage, time%tend, state%mesh%numlay)

         ! Phase-2 tillage seed (validation + setup), gated on swtill.
         if (config%soil%swtill == 1)       call tillage_seed(state)

         ! Heat must init BEFORE soilwater_seed so hconduc can read
         ! state%heat%tsoil(node) during hydraulic-conductivity init.
         call state%heat%init(config%heat, state%mesh%numnod)

         ! Phase-2 runtime seed: hatm, initial h/theta, initial fluxes.
         call soilwater_seed(state, config%soil%hydraulics)
         if (state%atmosphere%swusecn == 1) call cn_init(state)

         ! Modern drainage init (absorbs scalar seeding, L/zbotdr, surface_runoff fields).
         call state%drainage%init(config%drain, state%mesh%numnod)

         ! Modern solute init (gated). Reads config%solute%swsolu directly
         ! instead of via state%cfg (OD Step 8: state%cfg access retirement).
         if (config%solute%swsolu == 1) call state%solute%init(config%solute, config%soil, state%mesh%numnod)

         ! Modern surfacewater init (unconditional). Always seeds lightweight scalars
         ! (pondmx/rsro/rsroexp/swdra); heavy work (allocations, sttab, management
         ! periods) is gated on swdra==2 internally.
         call state%surfacewater%init(config%surface_water, config%drain, config%soil, state%mesh%numnod)

         if (time%flTemperature) call temperature_seed(state, config)

         if (time%flSolute) call solute_seed(state)

         ! Open output files and write headers.
         call csv_output_init(state)

      end associate

      call state%diag%info('swap', 'initialised project ' // trim(config%general%project))
      call state%diag%info('swap', 'grid: ' // trim(to_str(state%mesh%numnod)) // ' compartments, ' // &
                                              trim(to_str(state%mesh%numlay)) // ' soil layers')
      call state%diag%info('swap', 'bottom boundary: swbotb=' // trim(to_str(config%bottom_boundary%swbotb)) // &
                                   ', solute: swsolu=' // trim(to_str(config%solute%swsolu)))

   end subroutine swap_init_body

   subroutine swap_run_step(state, config)
      use csv_output,         only: csv_output_step
      use timecontrol_mod,    only: timecontrol_advance, timecontrol_reduce_dt, &
                                    timecontrol_day_end, itertime_check
      use surfacewater_mod,   only: surfacewater_lateral, surfacewater_balance, surfacewater_year_reset
      use tillage_mod,        only: tillage_step, tillage_output
      use boundbottom_mod,    only: BoundBottom
      use meteo_mod,          only: ProcessMeteoDay
      use meteo_process_mod,  only: ReadMeteoDay
      use snow_mod,           only: snow_step
      use meteodt_mod,        only: MeteoDT
      use rootextraction_mod, only: RootExtraction
      use frozencond_mod,     only: FrozenCond, FrozenBounds
      use temperature_mod,    only: temperature_step
      use solute_mod,         only: solute_step
      use soilhydraulics_mod, only: soilwater_step, soilwater_update, soilwater_restore_state
      use irrigation_mod,     only: irrigation_step, ssdi_irrigation_step
      use drainage_mod,       only: drainage
      use error_mod,          only: error_collection_t, set_active_error_sink

      type(swap_state_t),  intent(inout), target :: state
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

      call state%diag%set_simtime(state%timecontrol%date, &
                                  state%timecontrol%daynr, state%timecontrol%daycum)

      block
         type(error_collection_t), pointer :: diag_errors
         logical,                  pointer :: diag_fatal
         diag_errors => state%diag%errors
         diag_fatal  => state%diag%fatal_raised
         call set_active_error_sink(diag_errors, diag_fatal)
      end block

      associate (time => state%timecontrol, &
                 crop => state%crop)

         ! Meteo input for the year / day.
         if (time%flYearStart) call ReadMeteoYear(state, config)

         if (time%flDayStart) then
            call ReadMeteoDay(state, config)
            call CropGrowth(1, state%heat%tsoil, state)
            if (time%flIrrigate) call irrigation_step(state)
            call ProcessMeteoDay(state, config)
            if (config%soil%swtill == 1) call tillage_step(state)
         end if

         if (time%flmeteodt .or. time%fletsine) call MeteoDT(state)

         ! Snow: MH+MM — possibly belongs inside the flDayStart block above,
         ! before ProcessMeteoDay.
         if (time%flSnow .and. time%flDayStart) call snow_step(state)

         ! Conductivity reduction for frozen conditions.
         if (config%soil%frost%swfrost == 1) call FrozenCond(state, config)

         ! Potential and actual root water extraction profile.
         call RootExtraction(state)

         ! SoilWater bottom boundary conditions.
         call BoundBottom(state, config)

         ! Inner timestep-reduction loop. Two writers can request a smaller
         ! dt: SurfaceWater(2/3) sets request_smaller_dt on oscillation /
         ! ponding limit; headcalc inside SoilWater(2) sets time%fldecdt on
         ! Richards non-convergence. Either signal aborts the substep,
         ! calls timecontrol_reduce_dt, and re-runs the loop.
         time%fldtreduce = .true.
         do while (time%fldtreduce)
            time%fldtreduce = .false.

            if (time%flDrain) call Drainage(state)

            if (.not.time%fldecdt .and. time%flSurfaceWater) call surfacewater_lateral(state, request_smaller_dt)
            if (request_smaller_dt) time%fldecdt = .true.

            if (config%soil%frost%swfrost == 1) call FrozenBounds(state, config)

            if (.not.time%fldecdt) call soilwater_step(state)

            if (.not.time%fldecdt .and. time%flSurfaceWater) call surfacewater_balance(state, request_smaller_dt)
            if (request_smaller_dt) time%fldecdt = .true.

            if (time%fldecdt) then
               call soilwater_restore_state(state)
               call timecontrol_reduce_dt(state)
               time%fldtreduce = .true.
            end if
         end do

         ! SoilWater rate/state variables.
         call soilwater_update(state)

         if (state%diag%aborted()) return

         if (time%flTemperature) call temperature_step(state, config)
         if (time%flSolute)      call solute_step(state)

         ! Update time variables and switches/flags.
         call timecontrol_advance(state)

         ! End-of-day section.
         if (time%flDayEnd) then
            ! [T2-A / ADR 0053] crop<->water exchange: hand the day's transpiration
            ! to the crop through the exchange record (potential from ET, actual
            ! from the water solve), so the crop reads it there — not from
            ! state%atmosphere / state%soilwater directly.
            state%exchange%crop_water%pot_transp = state%atmosphere%ptra
            state%exchange%crop_water%act_transp = state%soilwater%tra
            ! [ADR 0052] SoilManagement (WOFOST-N soil nutrient) calls removed.
            if (crop%common%flCropCalendar) call CropGrowth(2, state%heat%tsoil, state)  ! potential growth
            if (crop%common%flCropCalendar) call CropGrowth(3, state%heat%tsoil, state)  ! actual growth
            if (crop%common%flCropCalendar) call CropGrowth(4, state%heat%tsoil, state)  ! harvest

            ! Watchdog against (near-)endless simulations.
            if (time%flMaxIterTime) call itertime_check(state)

            ! Subsurface irrigation: decide for next day, adjust dt for dt_SSDI_event.
            if (config%irrigation%swssdi == 1) call ssdi_irrigation_step(state)
            call timecontrol_day_end(state)
         end if

         ! Output.
         if (time%floutput) then
            call csv_output_step(state)
            if (config%soil%swtill == 1) call tillage_output(state)
            if (time%flSurfaceWater) then
               if (time%daynr == merge(366, 365, dtleap(time%iyear))) &
                  call surfacewater_year_reset(state%surfacewater)
            end if
         else
            if (time%floutputshort) call csv_output_step(state)
         end if

         ! [ADR 0052] SoilManagement(6) removed (WOFOST-N detached).

      end associate

   end subroutine swap_run_step

   subroutine swap_close(state, config)
      use timecontrol_mod,      only: itertime_close
      use csv_output,           only: csv_output_finalize
      use error_mod,            only: clear_active_error_sink
      type(swap_state_t),  intent(inout) :: state
      type(swap_config_t), intent(in)    :: config  ! unused: kept for parallel signature with swap_init/swap_run_step

      ! Drop the per-instance fatal sink so the module pointer can't dangle to
      ! this instance after it is closed/deallocated (it is re-registered at the
      ! top of each swap_run_step, so this is purely defensive).
      call clear_active_error_sink()

      call itertime_close(state)
      call csv_output_finalize(state)
      ! [ADR 0052] SoilManagement(7) removed (WOFOST-N detached).

      ! Okay-file for external runners.
      call WriteSwapOk(config%general%project)

      call state%diag%info('swap', 'simulation complete for project ' // trim(config%general%project))

   end subroutine swap_close

end module swap_mod
