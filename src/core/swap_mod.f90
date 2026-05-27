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
      type(swap_config_t), target, intent(out) :: config  ! target: crop_config_global => config%crop (set inside body)

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
      ! Init pipeline: TOML load (in swap_init) → seed_state_from_config
      ! adapter → modern type-bound state%X%init calls → a handful of
      ! legacy magic-int dispatchers (X(1, state)) that have not yet been
      ! converted.
      !
      ! Strangler-fig follow-ups (the code still works; these are tracked
      ! debt, not bugs):
      !   - seed_state_from_config has ~10 HACK Phase 4f-extend slots for
      !     legacy globals not yet covered by typed schema (SWREDU, RSIGNI,
      !     CFEVAPPOND, iHWCKmodel, RDS, ksatexm path, …). Each slot is a
      !     small typed-config extension + adapter wiring. Per ADR 0016
      !     the deeper follow-on is the config-passing refactor.
      !   - Several per-subsystem post-init seeding blocks live here because
      !     state%X%init allocates the arrays the seeding writes into:
      !       * soilwater layer flats (ksatexm/ksatfit/cofani/orgmat/…)
      !       * drainage scalars + L/zbotdr per-level geometry
      !       * swinco==3 warm-restart h-profile CSV re-read
      !       * atmosphere snow handshake (snowinco ↔ ssnow)
      !     Each should fold into the relevant state%X%init.
      !   - Legacy magic-int dispatchers (DoTillage(1),
      !     SoilWater(1), Temperature(1), Solute(1)) should
      !     migrate to type-bound state%X%init matching the surrounding
      !     modern calls. SurfaceWater task dispatch retired (Task 9).
      !   - The legacy readswap() entry point survives in src/io/readswap.f90
      !     only as a parity-test fixture (ADR 0007); not called at runtime.
      !
      use tillage_mod,                only: tillage_seed
      use swap_log,                   only: log_info
      use runoff_mod,                 only: cn_init
      use temperature_mod,            only: temperature_seed
      use solute_mod,                 only: solute_seed
      use soilgrid_mod,               only: CalcGrid
      use soilhydraulics_mod,         only: soilwater_seed
      use seed_state_from_config_mod, only: seed_state_from_config
      use timecontrol_mod,            only: timecontrol_init, itertime_init
      use csv_output,                 only: csv_output_init

      type(swap_state_t),          intent(out)   :: state
      type(swap_config_t), target, intent(inout) :: config  ! target: crop_config_global => config%crop

      ! Non-owning config pointer; lifetime matches state's. Top-level
      ! compute routines read switches via state%cfg%X%Y without needing
      ! a separate config arg.
      state%cfg => config

      ! Iteration / timing statistics.
      call itertime_init(state)

      ! STRANGLER — bulk seeding adapter (state%timecontrol + cross-subsystem
      ! exceptions). Shrinks as typed schema slots cover more fields.
      call seed_state_from_config(config, state)

      ! Crop seeding must run BEFORE timecontrol_init because
      ! timecontrol_init reads state%crop%common%croptype(icrop).
      call state%crop%init(config%crop, config%meteo, config%general%pathwork)

      ! Time variables and switches/flags.
      call timecontrol_init(state)

      ! Grid parameters (writes directly to state%mesh).
      call CalcGrid(state, config)

      ! Modern type-bound state init.
      call state%soilwater%init(config%soil, config%bottom_boundary, &
                                config%general%pathwork, &
                                state%mesh%numnod, state%mesh%numlay)
      call state%nutrients%init(state%mesh%numlay, config%nutrients, config%general%pathwork)
      call state%crop%irrigation%init(config%irrigation, state%timecontrol%tstart, &
                                      state%timecontrol%tend, state%mesh, &
                                      config%general%pathwork)

      ! STRANGLER — soilwater layer flats. Seeded here because
      ! state%soilwater%init allocates the nlay-sized arrays. Fold into
      ! state%soilwater%init (multi-source rules below would move with it).
      ! cofani: drain.cofani first, soil.cofani overrides (soil wins).
      ! orgmat: soil.orgmat first; heat.porg backfills when absent.
      if (allocated(config%soil%hydraulics%ksatexm)) &
         state%soilwater%ksatexm(:) = config%soil%hydraulics%ksatexm(1:size(state%soilwater%ksatexm))
      if (allocated(config%soil%hydraulics%ksatfit)) &
         state%soilwater%ksatfit(:) = config%soil%hydraulics%ksatfit(1:size(state%soilwater%ksatfit))
      if (allocated(config%drain%cofani)) &
         state%soilwater%cofani(1:size(config%drain%cofani)) = config%drain%cofani
      if (allocated(config%soil%cofani)) &
         state%soilwater%cofani(1:size(config%soil%cofani))  = config%soil%cofani
      state%soilwater%flksatexm = .false.   ! never set in adapter
      if (allocated(config%soil%orgmat)) &
         state%soilwater%orgmat(1:size(config%soil%orgmat))  = config%soil%orgmat
      if (.not. allocated(config%soil%orgmat) .and. allocated(config%heat%porg)) &
         state%soilwater%orgmat(1:min(size(config%heat%porg), size(state%soilwater%orgmat))) = &
            config%heat%porg(1:min(size(config%heat%porg), size(state%soilwater%orgmat)))
      if (allocated(config%heat%psand)) &
         state%soilwater%psand(:) = config%heat%psand(1:size(state%soilwater%psand))
      if (allocated(config%heat%psilt)) &
         state%soilwater%psilt(:) = config%heat%psilt(1:size(state%soilwater%psilt))
      if (allocated(config%heat%pclay)) &
         state%soilwater%pclay(:) = config%heat%pclay(1:size(state%soilwater%pclay))
      state%soilwater%swbotb_runtime = config%bottom_boundary%swbotb
      state%soilwater%q0    = 0.0d0
      state%soilwater%k1max = 0.0d0
      state%soilwater%H0max = 0.0d0

      ! Modern atmosphere init (zero flat scalars + cohorts + snapshot
      ! config-derived params: snowcoef/swsublim/swetsine).
      call state%atmosphere%init(config)

      ! swinco==3 warm-restart: pond/pondini/dt/h-profile/atmosphere read
      ! from soil.initial. swinco<3: pondini/pond from soil.pondini;
      ! atmosphere already zero-init from %init above.
      ! STRANGLER — the h-profile CSV re-read inside this branch depends
      ! on state%soilwater%init having allocated state%soilwater%h above.
      ! Fold into state%soilwater%init (pass h_file path through) or
      ! split out a dedicated warm-restart helper.
      if (config%soil%swinco == 3 .and. &
          allocated(config%soil%initial%h_file) .and. &
          len_trim(config%soil%initial%h_file) > 0) then
         state%soilwater%pondini = config%soil%initial%pond
         state%soilwater%pond    = config%soil%initial%pond
         ! soil.initial.dt supersedes simulation%numerical%dt for swinco=3
         state%timecontrol%dt = config%soil%initial%dt
         block
            use csv_reader_mod, only: read_csv_table
            use error_mod,      only: error_collection_t
            real(8), allocatable     :: tbl(:,:)
            type(error_collection_t) :: errs
            character(len=2)         :: hdr(2)
            integer :: nrows, ki
            hdr(1) = 'z '
            hdr(2) = 'h '
            call read_csv_table(trim(config%soil%initial%h_file), hdr, tbl, errs)
            call errs%abort_if_fatal()
            nrows = size(tbl, 1)
            do ki = 1, min(nrows, size(state%soilwater%h))
               state%soilwater%h(ki) = tbl(ki, 2)
            end do
         end block
         ! atmosphere warm-restart fields (spev/saev not in config; stay zero).
         state%atmosphere%ssnow = config%soil%initial%ssnow
         state%atmosphere%ldwet = config%soil%initial%ldwet
         state%atmosphere%slw   = config%soil%initial%slw
         if (config%meteo%snow%swsnow /= 1) state%atmosphere%ssnow = 0.0d0
      else
         state%soilwater%pondini = config%soil%pondini
         state%soilwater%pond    = config%soil%pondini   ! legacy alias for swinco<3
      end if

      associate (time => state%timecontrol)

         ! Modern type-bound init. The swtill gate is inside
         ! state%tillage%init (skips Group AB if events not allocated).
         call state%tillage%init(config%soil%tillage, time%tend, state%mesh%numlay)

         ! LEGACY-INIT — magic-int dispatchers. Migrate to type-bound init.
         if (state%cfg%soil%swtill == 1)       call tillage_seed(state)

         ! Heat must init BEFORE SoilWater(1) so hconduc can read
         ! state%heat%tsoil(node) during hydraulic-conductivity init.
         call state%heat%init(config%heat, state%mesh%numnod)

         ! Phase-2 runtime seed: hatm, initial h/theta, initial fluxes.
         call soilwater_seed(state)
         if (state%atmosphere%swusecn == 1) call cn_init(state)

         ! Modern drainage init.
         call state%drainage%init(config%drain, state%mesh%numnod)

         ! STRANGLER — drainage scalar/L/zbotdr seeding. Fold into
         ! state%drainage%init. Multi-source rules: dramet==2 uses scalars
         ! (lm in m → cm, zbotdr_basic); per-level uses arrays (already cm).
         state%drainage%nrlevs     = config%drain%nrlevs
         state%drainage%swdivd     = config%drain%swdivd
         state%drainage%swnrsrf    = config%drain%surface_runoff%swnrsrf
         state%drainage%swtopnrsrf = config%drain%surface_runoff%swtopnrsrf
         state%drainage%swdivdinf  = config%drain%surface_runoff%swdivdinf
         state%drainage%FacDpthInf = config%drain%surface_runoff%facdpthinf
         if (config%drain%dramet == 2) then
            state%drainage%L(1)      = 100.0d0 * config%drain%lm
            state%drainage%zbotdr(1) = config%drain%zbotdr_basic
         else
            if (allocated(config%drain%L)) &
               state%drainage%L(1:size(config%drain%L)) = config%drain%L
            if (allocated(config%drain%zbotdr)) &
               state%drainage%zbotdr(1:size(config%drain%zbotdr)) = config%drain%zbotdr
         end if

         ! Modern solute init (gated).
         if (state%cfg%solute%swsolu == 1) call state%solute%init(config%solute, state%mesh%numnod)

         ! Modern surfacewater init (gated).
         if (time%flSurfaceWater) call state%surfacewater%init(config%surface_water, config%drain, state%mesh%numnod)

         if (time%flTemperature) call temperature_seed(state, config)

         ! STRANGLER — snow handshake (snowinco ↔ ssnow). Fold into
         ! state%atmosphere%init or a dedicated snow_init that reads swinco.
         if (time%flSnow) then
            if (config%soil%swinco == 3) then
               state%atmosphere%snowinco = state%atmosphere%ssnow
            else
               state%atmosphere%ssnow = state%atmosphere%snowinco
            end if
         end if

         if (time%flSolute) call solute_seed(state)

         ! Open output files and write headers.
         call csv_output_init(state)

      end associate

      call log_info('swap', 'Initialization complete for project: ' // trim(state%cfg%general%project))

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
      use management_soil_mod, only: SoilManagement
      use drainage_mod,       only: drainage

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

      associate (time => state%timecontrol, &
                 crop => state%crop)

         ! Meteo input for the year / day.
         if (time%flYearStart) call ReadMeteoYear(state, config)

         if (time%flDayStart) then
            call ReadMeteoDay(state, config)
            call CropGrowth(1, state%heat%tsoil, state)
            if (time%flIrrigate) call irrigation_step(state)
            call ProcessMeteoDay(state, config)
            if (state%cfg%soil%swtill == 1) call tillage_step(state)
         end if

         if (time%flmeteodt .or. time%fletsine) call MeteoDT(state)

         ! Snow: MH+MM — possibly belongs inside the flDayStart block above,
         ! before ProcessMeteoDay.
         if (time%flSnow .and. time%flDayStart) call snow_step(state)

         ! Conductivity reduction for frozen conditions.
         if (state%cfg%soil%frost%swfrost == 1) call FrozenCond(state, config)

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

            if (state%cfg%soil%frost%swfrost == 1) call FrozenBounds(state, config)

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

         if (time%flTemperature) call temperature_step(state, config)
         if (time%flSolute)      call solute_step(state)

         ! Update time variables and switches/flags.
         call timecontrol_advance(state)

         ! End-of-day section.
         if (time%flDayEnd) then
            if (crop%common%flCropNut)      call SoilManagement(2, state)  ! nutrient status
            if (crop%common%flCropCalendar) call CropGrowth(2, state%heat%tsoil, state)  ! potential growth
            if (crop%common%flCropNut)      call SoilManagement(5, state)  ! crop-residue amendment
            if (crop%common%flCropNut)      call SoilManagement(3, state)  ! fertilizer amendment
            if (crop%common%flCropCalendar) call CropGrowth(3, state%heat%tsoil, state)  ! actual growth
            if (crop%common%flCropNut)      call SoilManagement(4, state)  ! nutrient processes
            if (crop%common%flCropCalendar) call CropGrowth(4, state%heat%tsoil, state)  ! harvest

            ! Watchdog against (near-)endless simulations.
            if (time%flMaxIterTime) call itertime_check(state)

            ! Subsurface irrigation: decide for next day, adjust dt for dt_SSDI_event.
            if (state%cfg%irrigation%swssdi == 1) call ssdi_irrigation_step(state)
            call timecontrol_day_end(state)
         end if

         ! Output.
         if (time%floutput) then
            call csv_output_step(state)
            if (state%cfg%soil%swtill == 1) call tillage_output(state)
            if (time%flSurfaceWater) then
               if (time%daynr == merge(366, 365, dtleap(time%iyear))) &
                  call surfacewater_year_reset(state%surfacewater)
            end if
         else
            if (time%floutputshort) call csv_output_step(state)
         end if

         if (time%flDayEnd .and. crop%common%flCropNut) call SoilManagement(6, state)

      end associate

   end subroutine swap_run_step

   subroutine swap_close(state, config)
      use swap_log,             only: log_info
      use management_soil_mod,  only: SoilManagement
      use timecontrol_mod,      only: itertime_close
      use csv_output,           only: csv_output_finalize
      type(swap_state_t),  intent(inout) :: state
      type(swap_config_t), intent(in)    :: config  ! unused: kept for parallel signature with swap_init/swap_run_step

      call itertime_close(state)
      call csv_output_finalize(state)
      if (state%crop%common%flCropNut) call SoilManagement(7, state)

      ! Okay-file for external runners.
      call WriteSwapOk(state%cfg%general%project)

      call log_info('swap', 'Simulation complete for project: ' // trim(state%cfg%general%project))

   end subroutine swap_close

end module swap_mod
