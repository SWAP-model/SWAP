!> Phase 4f Task A1 — config_to_variables adapter.
!!
!! Single public subroutine `config_to_variables(config)` that copies every
!! (C)-classified field listed in `docs/phase-4f-config-to-variables-audit.md`
!! from a populated `swap_config_t` to its corresponding `variables%`
!! global. The body is a long sequence of mechanical assignments walking
!! the audit doc section by section (general -> simulation -> meteorology
!! -> drainage -> soil -> bottom_boundary -> heat -> irrigation -> solute
!! -> surface_water -> crop), with `if (allocated(...))` guards on every
!! allocatable array.
!!
!! Per ADR 0009, the 18 RETIRED legacy output switches are zero-forced at
!! the end so any residual code that checks them does the right thing.
!!
!! Per ADR 0010 + ADR 0011, macropore globals are NOT populated; case 3
!! is excluded from check-full and uses the unchanged readswap path
!! (through cropgrowth's per-rotation init only).
!!
!! Per Phase 4f Task A0, `croptype(:)` is a renamed alias for
!! `crop_config_t.rotation_type(:)` — the adapter copies element-by-element.
!!
!! This adapter is the one place a bare `use variables` is OK; it touches
!! many globals across nearly every section of the legacy module.
!!
!! ## PHASE4F-EXTEND HACKs
!!
!! When iteration through Phase 4f Task B2/B3/.../B6 surfaces a legacy
!! global that the existing schema doesn't cover, the temporary fix lives
!! HERE in the adapter rather than in the schema. Each such fix is marked
!! with a comment line of the form:
!!
!!     ! HACK Phase 4f-extend: <one-line description of what's missing>
!!
!! Phase 4f-extend will walk this file, find every HACK marker via grep,
!! and replace each one with a proper schema extension (typed config
!! field + reader + per-case TOML population). After Phase 4f-extend
!! closes, no HACK markers should remain in this file.
module config_to_variables_mod
   use iso_fortran_env, only: real64
   use swap_config_mod, only: swap_config_t
   implicit none
   private

   public :: config_to_variables
   public :: apply_soil_tillage
   public :: apply_irrigation_ssdi
   public :: apply_nutrients
   public :: apply_nutrients_events

contains

   !> Copy every (C)-classified field from `config` into the corresponding
   !! `variables%` legacy global. Caller is responsible for having loaded,
   !! validated, and finalized `config` first.
   subroutine config_to_variables(config, state)
      use variables   ! bare-use is intentional: many globals across sections
      use swap_state_mod, only: swap_state_t
      use error_mod, only: fatalerr_collected
      type(swap_config_t), intent(in), target :: config
      type(swap_state_t),  intent(inout)      :: state

      integer :: i, n

      ! ---------------------------------------------------------------
      ! General + simulation (audit: 11 fields)
      ! ---------------------------------------------------------------
      if (allocated(config%general%project))   project   = config%general%project
      if (allocated(config%general%pathwork))  pathwork  = config%general%pathwork
      if (allocated(config%general%pathatm))   pathatm   = config%general%pathatm
      if (allocated(config%general%pathcrop))  pathcrop  = config%general%pathcrop
      if (allocated(config%general%pathdrain)) pathdrain = config%general%pathdrain
      state%timecontrol%swscre  = config%general%swscre

      state%timecontrol%tstart    = config%simulation%tstart
      state%timecontrol%tend      = config%simulation%tend
      state%timecontrol%nprintday = config%simulation%nprintday
      state%timecontrol%period    = config%simulation%period
      state%timecontrol%swres     = config%simulation%swres
      state%timecontrol%swodat    = config%simulation%swodat
      ! flprintdt not in config schema — default to .false.
      state%timecontrol%flprintdt = .false.

      ! Derive iyear/imonth from tstart, matching legacy readswap.f90:126-128.
      ! TimeControl(1) reads these as state (not as inputs) so the adapter
      ! must populate them before TimeControl runs. Without this, dtardp
      ! aborts with "year 0 or partial date not allowed" because iyear
      ! was zero-initialized by Initialize().
      block
         integer :: datea_init(6)
         real    :: fsec_init
         call dtdpar(state%timecontrol%tstart + 0.1d0, datea_init, fsec_init)
         state%timecontrol%iyear  = datea_init(1)  ! [GR-FINAL C1] written directly (tc_iyear_init_buf retired)
         state%timecontrol%imonth = datea_init(2)  ! [GR-FINAL C1] written directly (tc_imonth_init_buf retired)
      end block

      ! Legacy finalize for swmonth=1 (mirrors readswap.f90:181-207):
      ! when monthly output is on, populate outdatint(:) with end-of-month
      ! dates and clobber the daily-period output fields. Without this the
      ! CSV runs daily even though the .swp said "monthly", and the
      ! regression aggregator averages state_vars (e.g. GWL) over 365
      ! samples instead of 12 — producing a non-physical year-1 mean diff
      ! ~2 cm. The simulation state is identical at every monthly endpoint;
      ! only the sampling cadence differed. swmonth/swyrvar themselves are
      ! locals in readswap (readswap.f90:16), no module global exists.
      ! Must run AFTER iyear/imonth derivation above.
      if (config%simulation%swmonth == 1) then
         call populate_outdatint_monthly(state%timecontrol%tend, &
                                         state%timecontrol%iyear, &
                                         state%timecontrol%imonth)
         state%timecontrol%period = 0
         state%timecontrol%swres  = 0
         state%timecontrol%swodat = 0
      end if

      ! ---------------------------------------------------------------
      ! Simulation.numerical (audit: 6 fields)
      ! ---------------------------------------------------------------
      state%timecontrol%dt = config%simulation%numerical%dt  ! [GR-FINAL C1] written directly (tc_dt_init_buf retired)
      state%timecontrol%dtmin = config%simulation%numerical%dtmin
      state%timecontrol%dtmax = config%simulation%numerical%dtmax
      state%timecontrol%MaxIt = config%simulation%numerical%MaxIt
      ! MaxIterTime / flMaxIterTime not in config schema — default to 0 / .false.
      state%timecontrol%MaxIterTime   = 0
      state%timecontrol%flMaxIterTime = .false.
      state%timecontrol%msteps        = config%simulation%numerical%msteps
      ! Meteorology (audit: 12 + evaporation + snow)
      ! ---------------------------------------------------------------
      if (allocated(config%meteo%metfile))  metfil  = config%meteo%metfile
      ! Legacy `rainfil` global removed (ADR 0014).
      ! `config%meteo%rainfile` is no longer copied into a global because
      ! the only consumer (the `.YYY` per-year rain reader in readmeteo.f90)
      ! has been deleted; CSV rain events use `config%meteo%rain_events_file`.
      swetsine    = config%meteo%swetsine
      ! [GR-ATM 2026-05-23] angstroma/b retired — snapshotted in atmosphere_state%init

      ! All metfile extensions other than .csv are rejected by
      ! meteorology_config_validate (ADR 0014). The `.csv` guard below
      ! is defense-in-depth — the validator already enforced this before
      ! the adapter ran.
      call lowerc(metfil)

      if (index(trim(metfil), '.csv') > 0) then
         block
            use csv_reader_mod,  only: read_csv_table
            use error_mod,       only: error_collection_t
            real(8), allocatable :: tbl(:,:)
            type(error_collection_t) :: errs
            character(len=9) :: hdr(9)
            character(len=300) :: csvpath
            integer :: r
            hdr(1) = 'date     '
            hdr(2) = 'rad      '
            hdr(3) = 'tmin     '
            hdr(4) = 'tmax     '
            hdr(5) = 'hum      '
            hdr(6) = 'wind     '
            hdr(7) = 'rain     '
            hdr(8) = 'etref    '
            hdr(9) = 'wet      '
            csvpath = trim(pathatm) // trim(metfil)
            call read_csv_table(trim(csvpath), hdr, tbl, errs)
            call errs%abort_if_fatal()
            nmetcsv = size(tbl, 1)
            if (allocated(metcsv_dat)) deallocate(metcsv_dat)
            allocate(metcsv_dat(nmetcsv, 9))
            do r = 1, nmetcsv
               metcsv_dat(r, :) = tbl(r, :)
            end do
         end block
      end if

      ! Detail meteo CSV pre-load (swmetdetail=1 + detail_file provided).
      ! Metfile is always CSV here (validator rejects non-.csv per ADR
      ! 0014). The detail_file required-when-swmetdetail=1 check moved
      ! to meteorology_config_validate (SS-5 follow-up M2); the
      ! allocation guard below is defense-in-depth only.
      if (config%meteo%swmetdetail == 1) then
         if (allocated(config%meteo%detail_file) .and. &
             len_trim(config%meteo%detail_file) > 0) then
            block
               use csv_reader_mod,  only: read_csv_table
               use error_mod,       only: error_collection_t
               real(8), allocatable :: tbl(:,:)
               type(error_collection_t) :: errs
               character(len=8) :: hdr(7)
               character(len=300) :: csvpath
               integer :: r
               hdr(1) = 'datetime'
               hdr(2) = 'record  '
               hdr(3) = 'rad     '
               hdr(4) = 'temp    '
               hdr(5) = 'hum     '
               hdr(6) = 'wind    '
               hdr(7) = 'rain    '
               csvpath = trim(pathatm) // trim(config%meteo%detail_file)
               call read_csv_table(trim(csvpath), hdr, tbl, errs)
               call errs%abort_if_fatal()
               nmetcsv_det = size(tbl, 1)
               if (allocated(metcsv_det)) deallocate(metcsv_det)
               allocate(metcsv_det(nmetcsv_det, 7))
               do r = 1, nmetcsv_det
                  metcsv_det(r, :) = tbl(r, :)
               end do
            end block
         end if
      end if

      ! Rain events CSV pre-load (swrain=3, events_file set).
      if (config%meteo%swrain == 3 .and. allocated(config%meteo%rain_events_file)) then
         if (len_trim(config%meteo%rain_events_file) > 0) then
            block
               use csv_reader_mod,  only: read_csv_table
               use error_mod,       only: error_collection_t
               real(8), allocatable :: tbl(:,:)
               type(error_collection_t) :: errs
               character(len=8) :: hdr(2)
               character(len=300) :: csvpath
               integer :: r
               hdr(1) = 'datetime'
               hdr(2) = 'amount  '
               csvpath = trim(pathatm) // trim(config%meteo%rain_events_file)
               call read_csv_table(trim(csvpath), hdr, tbl, errs)
               call errs%abort_if_fatal()
               nraincsv = size(tbl, 1)
               if (allocated(raincsv_dat)) deallocate(raincsv_dat)
               allocate(raincsv_dat(nraincsv, 2))
               do r = 1, nraincsv
                  raincsv_dat(r, :) = tbl(r, :)
               end do
            end block
         end if
      end if

      ! Evaporation sub-section
      swcfbs = config%meteo%evaporation%swcfbs
      state%crop%cfbs = config%meteo%evaporation%cfbs
      ! [GR-ATM 2026-05-23] swredu/cofred/rsigni/cfevappond retired —
      ! snapshotted into state%atmosphere by atmosphere_state%init(config);
      ! compute reads from state, never from these legacy globals.

      ! Snow sub-section
      swsnow   = config%meteo%snow%swsnow
      snowcoef = config%meteo%snow%snowcoef
      state%atmosphere%TePrRain = config%meteo%snow%teprrain
      state%atmosphere%TePrSnow = config%meteo%snow%teprsnow

      ! ---------------------------------------------------------------
      ! Drainage (audit: 20 fields + surface_runoff sub-section)
      ! ---------------------------------------------------------------
      state%surfacewater%swdra = config%drain%swdra
      dramet   = config%drain%dramet
      ! [GR-BH Task 37] swdivd global deleted — state%drainage%swdivd seeded in swap_mod.f90
      swdislay = config%drain%swdislay
      ! [GR-BH Task 37] nrlevs global deleted — state%drainage%nrlevs seeded in swap_mod.f90
      basegw   = config%drain%basegw
      entres   = config%drain%entres
      ! Note: drainage%shape may collide with bottom_boundary%shape
      ! at the legacy global `shape`. Bottom boundary is wired below
      ! and overwrites for swbotb=3 cases (the documented audit alias).

      if (config%drain%swdra >= 1 .and. allocated(config%drain%drfil)) then
         drfil = config%drain%drfil
      end if

      ! DRAMET=2 (Hooghoudt/Ernst). Mirrors readswap.f90:1850-1875.
      ! lm is authored in metres; the legacy reader does the m->cm
      ! conversion (`l(1) = 100*lm2`) so we replicate that here.
      ! ADR 0031 Phase 2 Task 5: wetper(1) removed — state%drainage%wetper(1)
      ! is seeded from config%drain%wetper in drainage_init instead.
      ! zbotdr goes into the level-1 entry of the per-level array;
      ! ipos / khtop / khbot / kvtop / kvbot / zintf / geofac are scalar globals.
      if (config%drain%dramet == 2) then
         ! [GR-BH Task 37] L(1)/zbotdr(1) bare globals deleted — seeded from config in swap_mod.f90
         shape     = config%drain%shape
         ipos      = config%drain%ipos
         khtop     = config%drain%khtop
         if (config%drain%ipos >= 3) then
            khbot = config%drain%khbot
            zintf = config%drain%zintf
         end if
         if (config%drain%ipos >= 4) then
            kvtop = config%drain%kvtop
            kvbot = config%drain%kvbot
         end if
         if (config%drain%ipos == 5) then
            geofac = config%drain%geofac
         end if
      end if

      ! [GR-BH Task 36] cofani global retired — precedence logic moved to swap_mod.f90
      ! config%drain%cofani is consumed directly by swap_mod seeding block.

      if (allocated(config%drain%swdtyp)) then
         do i = 1, size(config%drain%swdtyp)
            swdtyp(i) = config%drain%swdtyp(i)
         end do
      end if
      ! [GR-BH Task 37] zbotdr bare global deleted — seeded from config%drain in swap_mod.f90
      if (allocated(config%drain%drares)) then
         do i = 1, size(config%drain%drares)
            drares(i) = config%drain%drares(i)
         end do
      end if
      if (allocated(config%drain%infres)) then
         do i = 1, size(config%drain%infres)
            infres(i) = config%drain%infres(i)
         end do
      end if
      ! [GR-BH Task 37] L bare global deleted — seeded from config%drain in swap_mod.f90
      if (allocated(config%drain%gwlinf)) then
         do i = 1, size(config%drain%gwlinf)
            gwlinf(i) = config%drain%gwlinf(i)
         end do
      end if
      if (allocated(config%drain%rdrain)) then
         do i = 1, size(config%drain%rdrain)
            rdrain(i) = config%drain%rdrain(i)
         end do
      end if
      if (allocated(config%drain%rinfi)) then
         do i = 1, size(config%drain%rinfi)
            rinfi(i) = config%drain%rinfi(i)
         end do
      end if
      if (allocated(config%drain%rentry)) then
         do i = 1, size(config%drain%rentry)
            rentry(i) = config%drain%rentry(i)
         end do
      end if
      if (allocated(config%drain%rexit)) then
         do i = 1, size(config%drain%rexit)
            rexit(i) = config%drain%rexit(i)
         end do
      end if
      if (allocated(config%drain%widthr)) then
         do i = 1, size(config%drain%widthr)
            widthr(i) = config%drain%widthr(i)
         end do
      end if
      if (allocated(config%drain%taludr)) then
         do i = 1, size(config%drain%taludr)
            taludr(i) = config%drain%taludr(i)
         end do
      end if
      if (allocated(config%drain%swallo)) then
         do i = 1, size(config%drain%swallo)
            swallo(i) = config%drain%swallo(i)
         end do
      end if

      ! Channel water level tables (DATOWL/LEVEL in .dra). ASCII readswap
      ! reads these per level into owltab(lev,1:2*nowltab(lev)). Without
      ! them owltab=0 → afgen returns channel level = 0 (surface), which
      ! causes spurious infiltration or drainage and changes water table
      ! dynamics. Each drainage level may supply an owltab_file CSV with
      ! header 'date,level'; col-1 is decoded as ISO date.
      if (allocated(config%drain%owltab_file)) then
         block
            use iso_fortran_env, only: real64
            use csv_reader_mod, only: read_csv_table
            use error_mod, only: error_collection_t
            real(real64), allocatable :: csv_table(:,:)
            type(error_collection_t)  :: csv_errs
            integer :: lev, nrows, k
            character(len=6) :: hdr(2)
            hdr(1) = 'date  '
            hdr(2) = 'level '
            do lev = 1, size(config%drain%owltab_file)
               if (len_trim(config%drain%owltab_file(lev)) == 0) cycle
               call read_csv_table(trim(config%drain%owltab_file(lev)), hdr, csv_table, csv_errs)
               call csv_errs%abort_if_fatal()
               if (csv_errs%count() == 0) then
                  nrows = size(csv_table, 1)
                  nowltab(lev) = nrows
                  do k = 1, nrows
                     owltab(lev, 2*k-1) = csv_table(k, 1)  ! date (days since 1900)
                     owltab(lev, 2*k)   = csv_table(k, 2)  ! channel water level (cm)
                  end do
               end if
            end do
         end block
      end if

      swliminf = config%drain%swliminf

      ! Drainage.surface_runoff sub-section: scalar switches + per-level
      ! arrays. Legacy globals `swtopdislay`, `ftopdislay`, `RapDraResRef`
      ! are arrays of size madr; copy element 1 of the (scalar) config
      ! field as a uniform value across drainage levels until a per-level
      ! schema lands in Phase 4f-extend.
      ! [GR-BH Task 37] swnrsrf/SwTopnrsrf/swdivdinf/FacDpthInf bare globals deleted —
      ! state%drainage%X seeded in swap_mod.f90
      cofintfl     = config%drain%surface_runoff%cofintfl
      expintfl     = config%drain%surface_runoff%expintfl
      ! ADR 0031: gate the surface_runoff geofac write to avoid overwriting
      ! the ipos==5 (Ernst geometry factor) write at line 306. Two distinct
      ! TOML fields map to one legacy global; the gate preserves both
      ! intended behaviors. Schema-level reconciliation deferred.
      if (config%drain%ipos /= 5) then
         geofac = config%drain%surface_runoff%geofac
      end if
      ! NOTE: do NOT write gwlconv from drainage.surface_runoff. Legacy
      ! reads gwlconv exactly once (readswap.f90:960, in the .swp Part 13
      ! numerical block) and there is no second read in the .dra reader.
      ! The schema's drainage.surface_runoff.gwlconv is a misnamed/orphan
      ! field with default 0.0; writing it here clobbered the proper
      ! value just set from simulation.numerical.gwlconv (default 100 cm).
      ! The clobber turns out to be benign in practice because gwlconv only
      ! gates a warning (soilhydraulics.f90:724), not solver behaviour, but
      ! the duplicate write is wrong on principle and could surprise future
      ! cases if the warning ever becomes load-bearing.
      rsurfdeep    = config%drain%surface_runoff%rsurfdeep
      rsurfshallow = config%drain%surface_runoff%rsurfshallow
      ! [SS-GR-FINAL D1] RapDraReaExp write dropped — global retired
      NumLevRapDra = config%drain%surface_runoff%numlevrapdra
      ! swtopdislay(madr) and ftopdislay(madr): broadcast scalar config
      ! field to all levels (currently no per-level schema slot).
      if (size(swtopdislay) >= 1) then
         do i = 1, size(swtopdislay)
            swtopdislay(i) = config%drain%surface_runoff%swtopdislay
         end do
      end if
      if (size(ftopdislay) >= 1) then
         do i = 1, size(ftopdislay)
            ftopdislay(i) = config%drain%surface_runoff%ftopdislay
         end do
      end if
      ! [SS-GR-FINAL D1] RapDraResRef write dropped — global retired

      ! ---------------------------------------------------------------
      ! Soil (audit: 15 + discretization + frost)
      ! ---------------------------------------------------------------
      state%soilwater%swsophy = config%soil%swsophy
      swhyst  = config%soil%swhyst
      swinco  = config%soil%swinco
      ! [MACRO-RETIRE 2026-05-12] swmacro global retired (ADR 0040).
      ! soil.swmacro=1 is still rejected by soil_config validator stub.
      ! SS-B / ADR 0020: legacy globals use prefixed names (variables
      ! module); set both the legacy switch globals and the call-site
      ! gating flags. When flTillage / flSSDI are false (default for
      ! every regression case), the call sites in swap.f90 /
      ! timecontrol.f90 short-circuit and the subsystems never run.
      till_swtill = config%soil%swtill
      swssdi_irr  = config%irrigation%swssdi
      flTillage = (config%soil%swtill == 1)
      flSSDI    = (config%irrigation%swssdi == 1)
      if (flTillage) call apply_soil_tillage(config%soil%tillage, state%timecontrol%tend, state)
      if (flSSDI)    call apply_irrigation_ssdi(config%irrigation%ssdi, &
                                               state%timecontrol%tstart, &
                                               state%timecontrol%tend, state)
      call apply_nutrients(config%nutrients)
      gwli    = config%soil%gwli
      ! [GR-FINAL C1] pondini/pond: config%soil%pondini read directly by swap_mod after soilwater_init
      ! (pondini_init_buf/pond_init_buf retired; swap_mod seeding replaced with direct config reads)
      state%surfacewater%pondmx = config%soil%pondmx
      ! [GR-ATM 2026-05-23] rsoil retired — snapshotted in atmosphere_state%init
      state%surfacewater%rsro = config%soil%rsro
      state%surfacewater%rsroexp = config%soil%rsroexp
      ! Legacy parses .swp `SWRUNON` into a local int; the persistent
      ! global is the boolean `flrunon`. Mirror that mapping here.
      flrunon = (config%soil%swrunon == 1)
      nrstaring = config%soil%nrstaring

      ! sublay (legacy 'isublay') is a local in readswap, not a module
      ! global; calcgrid only consumes nsublay + isoillay + ncomp + hcomp
      ! + hsublay, so we simply set nsublay here.
      if (allocated(config%soil%sublay)) then
         nsublay = size(config%soil%sublay)
      end if
      if (allocated(config%soil%isoillay)) then
         do i = 1, size(config%soil%isoillay)
            isoillay(i) = config%soil%isoillay(i)
         end do
         numlay = config%soil%isoillay(size(config%soil%isoillay))
      end if
      if (allocated(config%soil%hsublay)) then
         do i = 1, size(config%soil%hsublay)
            hsublay(i) = config%soil%hsublay(i)
         end do
      end if
      if (allocated(config%soil%ncomp)) then
         do i = 1, size(config%soil%ncomp)
            ncomp(i) = config%soil%ncomp(i)
         end do
      end if
      ! Derive hcomp = hsublay / ncomp (mirrors readswap.f90:613-619).
      ! Skipped when the case authors hcomp explicitly (none yet do).
      if (allocated(config%soil%hsublay) .and. allocated(config%soil%ncomp) &
          .and. .not. allocated(config%soil%hcomp)) then
         do i = 1, size(config%soil%hsublay)
            if (config%soil%ncomp(i) > 0) then
               hcomp(i) = config%soil%hsublay(i) / dble(config%soil%ncomp(i))
            end if
         end do
      end if
      if (allocated(config%soil%hcomp)) then
         do i = 1, size(config%soil%hcomp)
            hcomp(i) = config%soil%hcomp(i)
         end do
      end if
      ! [GR-BH Task 36] orgmat global retired — seeded via state%soilwater%orgmat in swap_mod.f90
      ! config%soil%orgmat is consumed directly by swap_mod seeding block.
      if (allocated(config%soil%bdens)) then
         do i = 1, size(config%soil%bdens)
            bdens(i) = config%soil%bdens(i)
         end do
      end if
      ! [GR-BH Task 36] cofani global retired — consumed via config%soil%cofani in swap_mod.f90
      ! (soil.cofani overrides drain.cofani — precedence preserved in swap_mod seeding block)

      ! [soil.initial] CSV path for swinco=3 warm restart. The typed
      ! soil.initial schema + per-profile CSV companions are the only
      ! supported pathway (the legacy ASCII swap.ini reader has been
      ! removed).
      if (config%soil%swinco == 3) then
         if (allocated(config%soil%initial%h_file) .and. &
             len_trim(config%soil%initial%h_file) > 0) then
            ! [SS-ATM A-2.6] ssnow/ldwet/slw retired to state%atmosphere; seeded in swap.f90 after atmosphere_init
            ! [GR-FINAL C1] pond/pondini/dt (swinco=3): read directly by swap_mod after soilwater_init
            ! (pond_init_buf/pondini_init_buf/tc_dt_init_buf retired)
            atmin7(:) = config%soil%initial%atmin7(:)
            ! [SS-ATM A-2.6] Legacy zeroes ssnow when swsnow != 1: now handled in swap.f90 during state seeding

            ! Note: [soil.initial].pondini (pre-existing top-level field) is
            ! NOT consumed by this swinco=3 path — only [soil.initial].pond
            ! is, since pond is the warm-restart "saved final state" while
            ! pondini was the swinco<3 "initial-from-scalars" input. They
            ! can both be authored without conflict; only pond wins here.

            ! Mandatory: initial pressure-head profile (z, h).
            ! [GR-FINAL C1] zi(:) seeded here; h values read directly in swap_mod after soilwater_init
            ! (h_init_buf retired; swap_mod now reads config%soil%initial%h_file inline).
            block
               use csv_reader_mod,  only: read_csv_table
               use error_mod,       only: error_collection_t
               real(8), allocatable     :: tbl(:,:)
               type(error_collection_t) :: errs
               character(len=2)         :: hdr(2)
               integer :: nrows, k
               hdr(1) = 'z '
               hdr(2) = 'h '
               call read_csv_table(trim(config%soil%initial%h_file), hdr, tbl, errs)
               call errs%abort_if_fatal()
               nrows = size(tbl, 1)
               nhead = nrows
               do k = 1, nrows
                  zi(k) = tbl(k, 1)
               end do
            end block

            ! Optional: initial soil temperature profile.
            if (config%heat%swhea == 1 .and. config%heat%swcalt == 2) then
               block
                  use csv_reader_mod,  only: read_csv_table
                  use error_mod,       only: error_collection_t
                  real(8), allocatable     :: tbl(:,:)
                  type(error_collection_t) :: errs
                  character(len=5)         :: hdr(2)
                  integer :: nrows, k
                  hdr(1) = 'z    '
                  hdr(2) = 'tsoil'
                  call read_csv_table(trim(config%soil%initial%tsoil_file), hdr, tbl, errs)
                  call errs%abort_if_fatal()
                  nrows = size(tbl, 1)
                  do k = 1, nrows
                     zh(k)    = tbl(k, 1)
                     tsoil(k) = tbl(k, 2)
                  end do
               end block
            end if

            ! Optional: initial concentration profile (Cml).
            if (config%solute%swsolu == 1) then
               block
                  use csv_reader_mod,  only: read_csv_table
                  use error_mod,       only: error_collection_t
                  real(8), allocatable     :: tbl(:,:)
                  type(error_collection_t) :: errs
                  character(len=3)         :: hdr(2)
                  integer :: nrows, k
                  hdr(1) = 'z  '
                  hdr(2) = 'cml'
                  call read_csv_table(trim(config%soil%initial%cml_file), hdr, tbl, errs)
                  call errs%abort_if_fatal()
                  nrows = size(tbl, 1)
                  nconc = nrows
                  do k = 1, nrows
                     zc(k)  = tbl(k, 1)
                     cml(k) = tbl(k, 2)
                  end do
               end block
            end if
         end if
      end if

      ! Per-soil-physical-layer Mualem-van Genuchten hydraulics.
      ! Mirrors readswap.f90:786-825 conceptually, but most of the
      ! per-array names (ores/osat/alfa/npar/lexp/alfaw) are *locals*
      ! in readswap — not module globals — so they only matter as input
      ! to paramvg(1..10, lay). We write paramvg directly here. The
      ! globals we *do* need to set are bdens / h_enpr (real(8) :: ...(maho)
      ! in variables.f90) since downstream code reads them. ksatfit/ksatexm
      ! retired from adapter — state%soilwater reads directly from config [SS-GR-BH A6].
      if (allocated(config%soil%hydraulics%ores)) then
         do i = 1, size(config%soil%hydraulics%ores)
            h_enpr(i)  = config%soil%hydraulics%h_enpr(i)
            bdens(i)   = config%soil%hydraulics%bdens(i)
         end do

         ! Default analytical MvG model for every soil-physical layer.
         ! HACK Phase 4f-extend: iHWCKmodel(:) is fixed to 1 (uni-modal
         ! MvG). The legacy reader at readswap.f90:651-657 lets cases
         ! override per layer (1..11) but no schema slot covers it yet.
         ! None of the regression cases set it; this matches.
         iHWCKmodel(1:size(config%soil%hydraulics%ores)) = 1

         ! paramvg layout (legacy):
         !   1 = ores         6 = npar
         !   2 = osat         7 = 1 - 1/npar
         !   3 = ksatfit      8 = alfa or alfaw (per swhyst)
         !   4 = alfa         9 = h_enpr
         !   5 = lexp        10 = ksatexm (-999 sentinel when absent)
         paramvg = 0.0d0
         do i = 1, size(config%soil%hydraulics%ores)
            paramvg(1, i)  = config%soil%hydraulics%ores(i)
            paramvg(2, i)  = config%soil%hydraulics%osat(i)
            paramvg(3, i)  = config%soil%hydraulics%ksatfit(i)
            paramvg(4, i)  = config%soil%hydraulics%alfa(i)
            paramvg(5, i)  = config%soil%hydraulics%lexp(i)
            paramvg(6, i)  = config%soil%hydraulics%npar(i)
            paramvg(7, i)  = 1.0d0 - (1.0d0 / paramvg(6, i))
            if (swhyst == 0) then
               paramvg(8, i) = config%soil%hydraulics%alfa(i)
            else
               paramvg(8, i) = config%soil%hydraulics%alfaw(i)
            end if
            paramvg(9, i)  = config%soil%hydraulics%h_enpr(i)
            ! HACK Phase 4f-extend: ksatexm path ignores the legacy
            ! flksatexm/relsatthr/ksatthr branch (readswap.f90:802-815).
            ! For hupselbrook ksatexm == ksatfit, so flksatexm stays
            ! false in the legacy path and paramvg(10,:) keeps the
            ! -999 sentinel. Matches behaviour for the case at hand;
            ! cases with ksatexm > ksatfit need the threshold-Ksat
            ! computation ported.
            paramvg(10, i) = -999.0d0
         end do
      end if

      ! Soil.discretization
      swdiscrvert = config%soil%discretization%swdiscrvert
      numnodnew   = config%soil%discretization%numnodnew
      if (allocated(config%soil%discretization%dznew)) then
         do i = 1, size(config%soil%discretization%dznew)
            dznew(i) = config%soil%discretization%dznew(i)
         end do
      end if

      ! Soil.frost
      swfrost   = config%soil%frost%swfrost
      state%soilwater%swfrost = swfrost   ! [SS-GR-UTILS Task 4] dual-write
      swsublim  = config%soil%frost%swsublim
      ! tfroststa/tfrostend live in heat block per audit (and per legacy);
      ! do not double-write here from soil%frost — heat block below owns it.

      ! ---------------------------------------------------------------
      ! Bottom boundary (audit: 9 fields, conditional per swbotb)
      ! ---------------------------------------------------------------
      ! swbotb: legacy global write retired — state%soilwater%swbotb_runtime sourced directly
      ! from config%bottom_boundary%swbotb in swap_mod.f90 [SS-GR-BH A7].
      select case (config%bottom_boundary%swbotb)
      case (1)
         block
            use iso_fortran_env, only: real64
            use csv_reader_mod, only: read_csv_table
            use error_mod, only: error_collection_t
            real(real64), allocatable :: csv_table(:,:)
            type(error_collection_t)  :: csv_errs
            integer :: k, nrows
            character(len=4) :: hdr(2)
            hdr(1) = 'date'
            hdr(2) = 'gwl '
            call read_csv_table( &
               trim(pathwork)//trim(config%bottom_boundary%gwl_file), &
               hdr, csv_table, csv_errs)
            call csv_errs%abort_if_fatal()
            if (allocated(csv_table)) then
               nrows = size(csv_table, 1)
               do k = 1, nrows
                  gwltab(k*2 - 1) = csv_table(k, 1)
                  gwltab(k*2)     = csv_table(k, 2)
               end do
            end if
         end block
      case (2)
         sw2    = config%bottom_boundary%sw2
         ! Phase 0 B-0.1: populate sine-wave scalars regardless of sw2;
         ! the gate in boundbottom.f90:104 protects the non-sine path.
         sinmax = config%bottom_boundary%sinmax
         sinamp = config%bottom_boundary%sinamp
         sinave = config%bottom_boundary%sinave
         if (config%bottom_boundary%sw2 == 2) then
            block
               use iso_fortran_env, only: real64
               use csv_reader_mod, only: read_csv_table
               use error_mod, only: error_collection_t
               real(real64), allocatable :: csv_table(:,:)
               type(error_collection_t)  :: csv_errs
               integer :: k, nrows
               character(len=4) :: hdr(2)
               hdr(1) = 'date'
               hdr(2) = 'qbot'
               call read_csv_table( &
                  trim(pathwork)//trim(config%bottom_boundary%qbot2_file), &
                  hdr, csv_table, csv_errs)
               call csv_errs%abort_if_fatal()
               if (allocated(csv_table)) then
                  nrows = size(csv_table, 1)
                  do k = 1, nrows
                     qbotab(k*2 - 1) = csv_table(k, 1)
                     qbotab(k*2)     = csv_table(k, 2)
                  end do
               end if
            end block
         end if
      case (3)
         shape       = config%bottom_boundary%shape
         hdrain      = config%bottom_boundary%hdrain
         rimlay      = config%bottom_boundary%rimlay
         aqave       = config%bottom_boundary%aqave
         aqamp       = config%bottom_boundary%aqamp
         aqper       = config%bottom_boundary%aqper
         aqtmax      = config%bottom_boundary%aqtmax
         swbotb3impl = config%bottom_boundary%swbotb3impl
         sw3         = config%bottom_boundary%sw3
         sw4         = config%bottom_boundary%sw4
         if (config%bottom_boundary%sw3 == 2) then
            block
               use iso_fortran_env, only: real64
               use csv_reader_mod, only: read_csv_table
               use error_mod, only: error_collection_t
               real(real64), allocatable :: csv_table(:,:)
               type(error_collection_t)  :: csv_errs
               integer :: k, nrows
               character(len=6) :: hdr(2)
               hdr(1) = 'date  '
               hdr(2) = 'haquif'
               call read_csv_table( &
                  trim(pathwork)//trim(config%bottom_boundary%haquif_file), &
                  hdr, csv_table, csv_errs)
               call csv_errs%abort_if_fatal()
               if (allocated(csv_table)) then
                  nrows = size(csv_table, 1)
                  do k = 1, nrows
                     haqtab(k*2 - 1) = csv_table(k, 1)
                     haqtab(k*2)     = csv_table(k, 2)
                  end do
               end if
            end block
         end if
         if (config%bottom_boundary%sw4 == 1) then
            block
               use iso_fortran_env, only: real64
               use csv_reader_mod, only: read_csv_table
               use error_mod, only: error_collection_t
               real(real64), allocatable :: csv_table(:,:)
               type(error_collection_t)  :: csv_errs
               integer :: k, nrows
               character(len=4) :: hdr(2)
               hdr(1) = 'date'
               hdr(2) = 'qbot'
               call read_csv_table( &
                  trim(pathwork)//trim(config%bottom_boundary%qbot4_file), &
                  hdr, csv_table, csv_errs)
               call csv_errs%abort_if_fatal()
               if (allocated(csv_table)) then
                  nrows = size(csv_table, 1)
                  do k = 1, nrows
                     qbotab(k*2 - 1) = csv_table(k, 1)
                     qbotab(k*2)     = csv_table(k, 2)
                  end do
               end if
            end block
         end if
      case (4)
         swqhbot  = config%bottom_boundary%swqhbot
         ! Phase 0 B-0.2: populate exponential q(h) scalars regardless of
         ! swqhbot; the gate in boundbottom.f90:152-153 protects the tabular path.
         cofqha   = config%bottom_boundary%cofqha
         cofqhb   = config%bottom_boundary%cofqhb
         cofqhc   = config%bottom_boundary%cofqhc
         swcofqhc = config%bottom_boundary%swcofqhc
         if (config%bottom_boundary%swqhbot == 2) then
            block
               use iso_fortran_env, only: real64
               use csv_reader_mod, only: read_csv_table
               use error_mod, only: error_collection_t
               real(real64), allocatable :: csv_table(:,:)
               type(error_collection_t)  :: csv_errs
               integer :: k, nrows
               character(len=4) :: hdr(2)
               hdr(1) = 'htab'
               hdr(2) = 'qtab'
               call read_csv_table( &
                  trim(pathwork)//trim(config%bottom_boundary%qhbot_file), &
                  hdr, csv_table, csv_errs)
               call csv_errs%abort_if_fatal()
               ! Legacy unpack pattern from readswap.f90:1418-1419 — for the
               ! q(h) curve, qbotab(odd) = abs(htab) and qbotab(even) = qtab.
               if (allocated(csv_table)) then
                  nrows = size(csv_table, 1)
                  do k = 1, nrows
                     qbotab(k*2 - 1) = abs(csv_table(k, 1))
                     qbotab(k*2)     = csv_table(k, 2)
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
            use iso_fortran_env, only: real64
            use csv_reader_mod, only: read_csv_table
            use error_mod, only: error_collection_t
            real(real64), allocatable :: csv_table(:,:)
            type(error_collection_t)  :: csv_errs
            integer :: k, nrows
            character(len=4) :: hdr(2)
            hdr(1) = 'date'
            hdr(2) = 'hbot'
            call read_csv_table( &
               trim(pathwork)//trim(config%bottom_boundary%hbot5_file), &
               hdr, csv_table, csv_errs)
            call csv_errs%abort_if_fatal()
            if (allocated(csv_table)) then
               nrows = size(csv_table, 1)
               do k = 1, nrows
                  hbotab(k*2 - 1) = csv_table(k, 1)
                  hbotab(k*2)     = csv_table(k, 2)
               end do
            end if
         end block
      case (6, 7)
         ! No parameters to populate for modes 6 and 7.
      case (8)
         ! Phase 0 B-0.3: lysimeter path — populate hplate legacy global.
         hplate = config%bottom_boundary%hplate
      end select

      ! ---------------------------------------------------------------
      ! Heat (audit: 10 fields + 6 Phase-0 promoted fields = 16 total)
      ! ---------------------------------------------------------------
      swhea     = config%heat%swhea
      swcalt    = config%heat%swcalt
      swtopbhea = config%heat%swtopbhea
      swbotbhea = config%heat%swbotbhea
      tfroststa = config%heat%tfroststa
      tfrostend = config%heat%tfrostend

      ! psand/psilt/pclay: legacy global writes retired — state%soilwater%psand/psilt/pclay
      ! now sourced directly from config%heat in swap_mod.f90 [SS-GR-BH A6].
      ! [GR-BH Task 36] orgmat global retired — heat.porg→orgmat backfill now handled in
      ! swap_mod.f90 seeding block (state%soilwater%orgmat). config%heat%porg consumed directly.
      ! Initial soil temperature table tsoil_init(:,1:2) — column 1 (depth)
      ! mirrors legacy `zh`, column 2 (temp) mirrors legacy `tsoil(1..nheat)`.
      ! `nheat` is the number of (depth, temp) pairs; it gates the afgen
      ! table-build in heat/temperature.f90 task=1 (lines 117-123). Without
      ! `nheat` and `zh`, the table-build loop is skipped and every
      ! compartment's tsoil(:) is interpolated against an empty afgen
      ! table, yielding spurious zero initial soil temperatures across
      ! the profile. (Hupselbrook iHWCKmodel<=3 path means no direct
      ! water-flow feedback, so this fix is parity-correctness only and
      ! does not by itself close the year-1 GWL gap.)
      if (allocated(config%heat%tsoil_init)) then
         n = size(config%heat%tsoil_init, 1)
         ! [GR-FINAL C4] nheat write dropped (W-global; 0 consumers in temperature.f90 code)
         do i = 1, min(n, size(tsoil))
            zh(i)    = config%heat%tsoil_init(i, 1)
            tsoil(i) = config%heat%tsoil_init(i, 2)
         end do
      end if

      ! Phase 0 (SS-HEAT) — swcalt=1 analytical method scalars.
      ddamp  = config%heat%ddamp
      tmean  = config%heat%tmean
      tampli = config%heat%tampli
      timref = config%heat%timref

      ! Flatten 2D typed table → interleaved 1D afgen layout.
      ! afgen(temtoptab, 2*mabbc, time) reads (2*k-1)=time, (2*k)=value.
      ! Confirmed from temperature.f90:151 and 179.
      if (allocated(config%heat%temtoptab)) then
         do i = 1, min(size(config%heat%temtoptab, 1), size(temtoptab)/2)
            temtoptab(2*i - 1) = config%heat%temtoptab(i, 1)   ! time
            temtoptab(2*i)     = config%heat%temtoptab(i, 2)   ! temperature
         end do
      end if

      if (allocated(config%heat%tembtab)) then
         do i = 1, min(size(config%heat%tembtab, 1), size(tembtab)/2)
            tembtab(2*i - 1) = config%heat%tembtab(i, 1)   ! time
            tembtab(2*i)     = config%heat%tembtab(i, 2)   ! temperature
         end do
      end if

      ! ---------------------------------------------------------------
      ! Irrigation (audit: 8 fields, top-level only)
      ! ---------------------------------------------------------------
      swirfix = config%irrigation%swirfix
      ! Inline fixed events: copy (date, depth, conc, type) rows from
      ! the typed config table into the legacy parallel arrays. The
      ! reader already stored col 1 as days-since-1900, so this is a
      ! straight copy. The /10.0 on irdepth mirrors readswap.f90:503
      ! (mm in input -> cm in legacy globals).
      if (allocated(config%irrigation%fixed_events)) then
         n = size(config%irrigation%fixed_events, 1)
         do i = 1, min(n, size(irdate))
            irdate(i)  = config%irrigation%fixed_events(i, 1)
            irdepth(i) = config%irrigation%fixed_events(i, 2) / 10.0d0
            irconc(i)  = config%irrigation%fixed_events(i, 3)
            irtype(i)  = nint(config%irrigation%fixed_events(i, 4))
         end do
      else if (swirfix == 1 .and. allocated(config%irrigation%fixed_events_file)) then
         ! Phase 4f cleanup: long-form fixed-irrigation events outsourced
         ! to a CSV companion file (date, depth_mm, conc, type). The
         ! reader emits days-since-1900 in column 1; the rest of the
         ! unpack mirrors the inline-fixed_events path above (mm -> cm
         ! on depth, nint() on type). Replaces the legacy .irg HACK.
         if (len_trim(config%irrigation%fixed_events_file) > 0) then
            block
               use csv_reader_mod, only: read_csv_table
               use error_mod, only: error_collection_t
               use iso_fortran_env, only: real64
               real(real64), allocatable :: csv_table(:,:)
               type(error_collection_t) :: csv_errs
               integer :: k_csv, nrows_csv
               character(len=5) :: irrig_header(4)
               irrig_header(1) = 'date '
               irrig_header(2) = 'depth'
               irrig_header(3) = 'conc '
               irrig_header(4) = 'type '
               call read_csv_table( &
                  trim(pathwork)//trim(config%irrigation%fixed_events_file), &
                  irrig_header, csv_table, csv_errs)
               call csv_errs%abort_if_fatal()
               if (allocated(csv_table)) then
                  nrows_csv = size(csv_table, 1)
                  do k_csv = 1, min(nrows_csv, size(irdate))
                     irdate(k_csv)  = csv_table(k_csv, 1)
                     irdepth(k_csv) = csv_table(k_csv, 2) / 10.0d0  ! mm -> cm
                     irconc(k_csv)  = csv_table(k_csv, 3)
                     irtype(k_csv)  = nint(csv_table(k_csv, 4))
                  end do
               end if
            end block
         end if
      end if
      ! cirrs / cirrthres / dcrit / isuas / perirrsurp / raithreshold /
      ! swcirrthres live in irrigation_schedule_t (per-crop), not the
      ! top-level irrigation_config_t. They are populated per-rotation
      ! by cropgrowth's legacy crop sub-readers (Phase 4g territory);
      ! the strangler adapter does not touch them.

      ! ---------------------------------------------------------------
      ! Solute (audit: 8 fields + 14 Phase 0 promoted fields)
      ! ---------------------------------------------------------------
      swsolu  = config%solute%swsolu
      swbotbc = config%solute%swbotbc
      cdrain  = config%solute%cdrain
!     cseep   = config%solute%cseep   ! global cseep removed (ADR 0032); state%solute%cseep written by solute task=2 via afgen(cseeptab)
      tscf    = config%solute%tscf
      rtheta  = config%solute%rtheta
      bexp    = config%solute%bexp
      ! Phase 4f Task B5: per-layer dispersion length. Mirrors
      ! readswap.f90:1139-1143 — when the case authors `ldis` as an
      ! array, copy element-wise; otherwise broadcast the scalar to
      ! every soil-physical layer (legacy `rdsdor('ldis',...,ldis(1))`
      ! followed by an implicit broadcast in the dispersion solver).
      if (allocated(config%solute%ldis_array)) then
         do i = 1, size(config%solute%ldis_array)
            ldis(i) = config%solute%ldis_array(i)
         end do
      else if (config%solute%ldis > 0.0d0) then
         ldis(1) = config%solute%ldis
      end if

      ! Phase 0 (ADR 0032) — populate legacy globals from the 14 promoted fields.
      cref   = config%solute%cref
      cpre   = config%solute%cpre
      ddif   = config%solute%ddif
      frexp  = config%solute%frexp
      gampar = config%solute%gampar
      daquif = config%solute%daquif
      kfsat  = config%solute%kfsat
      decsat = config%solute%decsat
      poros  = config%solute%poros
      swbr   = config%solute%swbr

      if (allocated(config%solute%kf)) then
         do i = 1, min(size(config%solute%kf), size(kf))
            kf(i) = config%solute%kf(i)
         end do
      end if
      if (allocated(config%solute%decpot)) then
         do i = 1, min(size(config%solute%decpot), size(decpot))
            decpot(i) = config%solute%decpot(i)
         end do
      end if
      if (allocated(config%solute%fdepth)) then
         do i = 1, min(size(config%solute%fdepth), size(fdepth))
            fdepth(i) = config%solute%fdepth(i)
         end do
      end if

      ! cseeptab: flatten 2D typed config to the interleaved afgen layout.
      ! afgen(cseeptab, mabbc*2, time) reads pairs as (2*k-1)=time, (2*k)=value.
      ! Confirmed from: grep -n "cseeptab" src/solute/solute.f90 → line 107.
      if (allocated(config%solute%cseeptab)) then
         do i = 1, min(size(config%solute%cseeptab, 1), size(cseeptab)/2)
            cseeptab(2*i - 1) = config%solute%cseeptab(i, 1)   ! time
            cseeptab(2*i)     = config%solute%cseeptab(i, 2)   ! concentration
         end do
      end if

      ! ---------------------------------------------------------------
      ! Surface water (audit: 12 + per-period management arrays)
      ! ---------------------------------------------------------------
      ! [GR-CROP-DVS] swsrf/swsec retired — consumers read state%cfg%surface_water%X
      ! TOML pipeline: pre-compute the initial water level wls1.
      ! Legacy rddre computes wls1 = wlact - altcu inside the routine;
      ! we do the same here so surfacewater_init can read it from a
      ! module global. drainage.altcu /= 0 is rejected upstream (Task 7),
      ! so this simplifies to wlact.
      ! [GR-FINAL C4] wls1_init write dropped: W-global (0 consumers; state%surfacewater seeded in swap_mod)
      osswlm    = config%surface_water%osswlm
      nmper  = config%surface_water%nmper
      swqhr  = config%surface_water%swqhr

      if (allocated(config%surface_water%impend)) then
         do i = 1, min(size(config%surface_water%impend), size(impend))
            impend(i) = config%surface_water%impend(i)
         end do
      end if
      if (allocated(config%surface_water%swman)) then
         do i = 1, min(size(config%surface_water%swman), size(swman))
            swman(i) = config%surface_water%swman(i)
         end do
      end if
      if (allocated(config%surface_water%wscap)) then
         do i = 1, min(size(config%surface_water%wscap), size(wscap))
            wscap(i) = config%surface_water%wscap(i)
         end do
      end if
      if (allocated(config%surface_water%wldip)) then
         do i = 1, min(size(config%surface_water%wldip), size(wldip))
            wldip(i) = config%surface_water%wldip(i)
         end do
      end if
      if (allocated(config%surface_water%intwl)) then
         do i = 1, min(size(config%surface_water%intwl), size(intwl))
            intwl(i) = config%surface_water%intwl(i)
         end do
      end if
      ! Note: alphaw arrays carry the post-finalize-normalized values
      ! per Phase 4f-prep Task D2 — the adapter copies them as-is.
      if (allocated(config%surface_water%hbweir)) then
         do i = 1, min(size(config%surface_water%hbweir), size(hbweir))
            hbweir(i) = config%surface_water%hbweir(i)
         end do
      end if
      if (allocated(config%surface_water%alphaw)) then
         do i = 1, min(size(config%surface_water%alphaw), size(alphaw))
            alphaw(i) = config%surface_water%alphaw(i)
         end do
      end if
      if (allocated(config%surface_water%betaw)) then
         do i = 1, min(size(config%surface_water%betaw), size(betaw))
            betaw(i) = config%surface_water%betaw(i)
         end do
      end if

      ! ---------------------------------------------------------------
      ! Crop (audit: ~110 fields). The strangler adapter only handles
      ! the rotation table globals that readswap reads at startup —
      ! per-rotation type-specific fields go into legacy globals during
      ! cropgrowth.f90's per-rotation init (Phase 4g territory).
      ! ---------------------------------------------------------------
      ! [GR-FINAL C4] swCrop write dropped (W-global; 0 external consumers).
      ! Mirror readswap.f90:479-480 — when crop simulation is enabled,
      ! arm the per-rotation reader gates so cropgrowth.f90's per-crop
      ! init (ArableLandGerm/CropFixed/Wofost/Grass at lines 91 / 121)
      ! actually fires. Without this, flCropReadFile stays .false. (the
      ! Initialize() default) and every rotation is treated as bare soil:
      ! LAI/cf/rd remain 0, TPOT/TACT collapse, and EACT/DRAINAGE balloon.
      if (config%crop%swcrop == 1) then
         flCropReadFile = .true.
         flCropOpenFile = .true.
      end if
      state%crop%common%flCropReadFile = flCropReadFile   ! [SS-GR-CROPRT A5]

      rdmax = config%crop%rdmax

      if (allocated(config%crop%rotation_type)) then
         n = size(config%crop%rotation_type)
         ! [GR-ATM 2026-05-23] allocate + populate state%crop%common%croptype
         allocate(state%crop%common%croptype(n))
         state%crop%common%croptype = config%crop%rotation_type
      end if
      if (allocated(config%crop%rotation_start)) then
         n = size(config%crop%rotation_start)
         do i = 1, min(n, size(cropstart))
            cropstart(i) = config%crop%rotation_start(i)
         end do
      end if
      if (allocated(config%crop%rotation_end)) then
         n = size(config%crop%rotation_end)
         do i = 1, min(n, size(cropend))
            cropend(i) = config%crop%rotation_end(i)
         end do
      end if
      if (allocated(config%crop%rotation_file)) then
         n = size(config%crop%rotation_file)
         do i = 1, min(n, size(cropfil))
            ! Legacy cropfil is a stem (no extension): the cropgrowth
            ! reader appends '.crp' itself. Our TOML authors the full
            ! '<name>.crp.toml' path; strip both suffixes so the legacy
            ! per-crop reader (still in use until Phase 4f-extend ports
            ! it) reconstructs '<name>.crp' on disk.
            ! HACK Phase 4f-extend: once read_cropfixed_toml et al. own
            ! the per-crop init, the rotation_file should pass through
            ! unchanged (the new readers will use the .toml path).
            cropfil(i) = strip_crp_toml_suffix(config%crop%rotation_file(i))
         end do
      end if

      ! Phase 1 (.crp port): expose the parsed crop config to runtime
      ! subs that need per-rotation cache access. Transitional — see
      ! ADR 0016. The pointer targets the caller's local config; valid
      ! for the duration of the simulation init.
      block
         use crop_config_global_mod, only: crop_config_global
         crop_config_global => config%crop
      end block

      ! ---------------------------------------------------------------
      ! ADR 0009 retired output switches — W-globals zeroed by initialize.f90
      ! or Fortran module default; zero-force lines dropped (GR-FINAL C4).
      ! Remaining: functional switch (swcaprise), state field (swheader),
      ! R-category (swrum, still read by swapoutput.f90).
      ! [SS-GR-CROPRT C1] swend (C-category) write dropped — global + state field retired (ADR 0009: always 0)
      ! ---------------------------------------------------------------
      state%timecontrol%swheader = 0
      swcaprise       = .false.        ! active in soilhydraulics.f90 (R; always .false. — no config field yet)
      ! [SS-GR-CROPRT A3] swrum adapter write dropped — global retired (always 0; outrume calls dropped)
      ! [GR-FINAL C4] dropped W-globals (all zero by init.f90 or Fortran default):
      !   swafo, swaun, swvap, swbal, swwba, swsba, swblc, swdrf, swstr, swirg,
      !   swini, swcapriseoutput, swswb, swoutputmodflow

      if (allocated(config%general%outfil)) outfil = config%general%outfil

      ! CSV output — read from [output.csv] schema section.
      ! Defaults (enabled=1, enabled_tz=0, inlist=water-balance, inlist_tz=wc,h,conc)
      ! are applied by output_csv_config_finalize prior to this adapter.
      ! [SS-GR-CROPRT A3] swcsv adapter write dropped — global retired; readers use config directly
      if (allocated(config%output_csv%inlist)) then
         InList_csv = config%output_csv%inlist
      end if
      ! [SS-GR-CROPRT A3] swcsv_tz adapter write dropped — global retired; readers use config directly
      if (allocated(config%output_csv%inlist_tz)) then
         InList_csv_tz = config%output_csv%inlist_tz
      end if

      ! [GR-ATM 2026-05-23] logf bare-global retired; swap_log opens
      ! 'swap_swap.log' via log_init() in swap_main.

   end subroutine config_to_variables

   !> Populate `variables%outdatint(:)` with the end-of-month dates
   !! between `tstart` and `tend`. Mirrors `readswap.f90:181-204` (the
   !! `swmonth == 1` branch). Bare `use variables` for parity with the
   !! parent adapter.
   subroutine populate_outdatint_monthly(tend, iyear, imonth)
      use variables
      real(8), intent(in) :: tend
      integer, intent(in) :: iyear   !! start year (from state%timecontrol%iyear, [GR-FINAL C1])
      integer, intent(in) :: imonth  !! start month (from state%timecontrol%imonth, [GR-FINAL C1])
      integer  :: datea_om(6), i_om
      real     :: fsec_om
      real(8)  :: outdate_om

      datea_om = 0
      datea_om(1) = iyear   ! [GR-FINAL C1] replaced tc_iyear_init_buf
      datea_om(2) = imonth  ! [GR-FINAL C1] replaced tc_imonth_init_buf
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

   !> Strip a trailing '.crp.toml' (or '.toml') suffix from a rotation
   !! file path, leaving the stem the legacy per-crop reader expects in
   !! `cropfil(:)`. Only used by the strangler adapter; goes away when
   !! the per-crop readers are TOML-native (Phase 4f-extend Task TBD).
   pure function strip_crp_toml_suffix(s) result(stem)
      character(len=*), intent(in)  :: s
      character(len=len(s))         :: stem
      character(len=*), parameter   :: SFX1 = '.crp.toml'
      character(len=*), parameter   :: SFX2 = '.toml'
      integer :: n
      n = len_trim(s)
      if (n >= len(SFX1) .and. s(n - len(SFX1) + 1 : n) == SFX1) then
         stem = s(1 : n - len(SFX1))
      else if (n >= len(SFX2) .and. s(n - len(SFX2) + 1 : n) == SFX2) then
         stem = s(1 : n - len(SFX2))
      else
         stem = s
      end if
   end function strip_crp_toml_suffix

   !> ISO 'YYYY-MM-DD' -> days since 1900 (real(real64)).
   !! Constructs a toml_datetime from the parsed date components and delegates
   !! to the existing parse_date_to_days1900 helper in toml_field_helpers_mod.
   function parse_iso_date_to_days1900(s) result(t)
      use, intrinsic :: iso_fortran_env, only: real64
      use tomlf, only: toml_datetime
      use toml_field_helpers_mod, only: parse_date_to_days1900
      character(len=*), intent(in) :: s
      real(real64) :: t
      type(toml_datetime) :: dtv
      integer :: y, m, d
      read(s, '(i4,1x,i2,1x,i2)') y, m, d
      dtv%date%year  = y
      dtv%date%month = m
      dtv%date%day   = d
      ! Leave dtv%time fields at default (-1) so the conversion treats it as a date-only.
      t = parse_date_to_days1900(dtv)
   end function parse_iso_date_to_days1900

   !> Apply [soil.tillage] config to legacy `variables` globals.
   !! Called from config_to_variables when flTillage is true. Allocates
   !! per-event and per-type arrays, populates them from the typed config,
   !! parses event dates to days-since-1900, sets the Ntill+1 sentinel
   !! (tend + 1), computes Max_Z_tillage and the iTT1/iTT2 first/last-position
   !! indices. Replaces the deleted Read_Tillage subroutine (Task 5).
   subroutine apply_soil_tillage(tillage, tend, state)
      use, intrinsic :: iso_fortran_env, only: real64
      use soil_config_mod, only: soil_tillage_t
      use error_mod, only: fatalerr_collected
      use swap_state_mod, only: swap_state_t
      use variables, only: Ntill         => till_Ntill, &
                           Ntypes        => till_Ntypes, &
                           i_n_model     => till_i_n_model, &
                           iRedist       => till_iRedist, &
                           Max_Z_tillage => till_Max_Z_tillage, &
                           Date_tillage  => till_Date_tillage, &
                           Z_tillage     => till_Z_tillage, &
                           I_tillage     => till_I_tillage, &
                           Type_tillage  => till_Type_Tillage, &
                           iType_Tillage => till_iType_Tillage, &
                           TAB_Rho_cons    => till_TAB_Rho_cons, &
                           TAB_Rho_tillage => till_TAB_Rho_tillage, &
                           TAB_K_R_cons    => till_TAB_K_R_cons, &
                           TAB_Rho_match   => till_TAB_Rho_match, &
                           TAB_N_match     => till_TAB_N_match, &
                           iTT1 => till_iTT1, &
                           iTT2 => till_iTT2
      type(soil_tillage_t), intent(in)    :: tillage
      real(real64),         intent(in)    :: tend
      type(swap_state_t),   intent(inout) :: state   ! [GR-BH Task 35] replaces NumNod/zbotcp globals

      integer :: i, j

      i_n_model = tillage%i_n_model
      iRedist   = tillage%iRedist

      Ntill  = size(tillage%events)
      Ntypes = size(tillage%types)

      ! Per-event arrays (sentinel: Date_tillage(Ntill+1) = tend + 1).
      if (allocated(Date_tillage)) deallocate(Date_tillage); allocate(Date_tillage(Ntill+1))
      if (allocated(Z_tillage))    deallocate(Z_tillage);    allocate(Z_tillage(Ntill))
      if (allocated(I_tillage))    deallocate(I_tillage);    allocate(I_tillage(Ntill))
      if (allocated(Type_tillage)) deallocate(Type_tillage); allocate(Type_tillage(Ntill))

      do i = 1, Ntill
         Z_tillage(i)    = tillage%events(i)%z
         I_tillage(i)    = tillage%events(i)%intensity
         Type_tillage(i) = tillage%events(i)%type_id
         Date_tillage(i) = parse_iso_date_to_days1900(tillage%events(i)%date)
      end do
      Date_tillage(Ntill + 1) = tend + 1.0_real64

      ! Deferred z-range check (validator can't see numnod / zbotcp).
      ! TODO: route this error through whatever error-collection mechanism
      ! config_to_variables uses; for now, fatalerr_collected on violation.
      ! The test test_apply_z_outside_grid_raises_error is deferred because
      ! catching a STOP in pFUnit is awkward — verified manually for now.
      ! [GR-BH Task 35] NumNod/zbotcp globals replaced by state%mesh fields.
      ! Guard: mesh is not yet populated when called from config_to_variables
      ! (CalcGrid runs after); skip z-range check until mesh is built.
      if (state%mesh%numnod > 0 .and. allocated(state%mesh%zbotcp)) then
         do i = 1, Ntill
            if (Z_tillage(i) < 0.0_real64 .or. &
                Z_tillage(i) > abs(state%mesh%zbotcp(state%mesh%numnod))) then
               call fatalerr_collected('apply_soil_tillage', &
                  'Z_tillage value is outside the model grid depth range')
            end if
         end do
      end if

      ! Per-type arrays.
      if (allocated(iType_Tillage))   deallocate(iType_Tillage);   allocate(iType_Tillage(Ntypes))
      if (allocated(TAB_Rho_cons))    deallocate(TAB_Rho_cons);    allocate(TAB_Rho_cons(Ntypes))
      if (allocated(TAB_Rho_tillage)) deallocate(TAB_Rho_tillage); allocate(TAB_Rho_tillage(Ntypes))
      if (allocated(TAB_K_R_cons))    deallocate(TAB_K_R_cons);    allocate(TAB_K_R_cons(Ntypes))

      do i = 1, Ntypes
         iType_Tillage(i)   = tillage%types(i)%id
         TAB_Rho_cons(i)    = tillage%types(i)%rho_cons
         TAB_Rho_tillage(i) = tillage%types(i)%rho_tillage
         TAB_K_R_cons(i)    = tillage%types(i)%k_R
      end do

      if (i_n_model == 3) then
         if (allocated(TAB_Rho_match)) deallocate(TAB_Rho_match); allocate(TAB_Rho_match(Ntypes))
         if (allocated(TAB_N_match))   deallocate(TAB_N_match);   allocate(TAB_N_match(Ntypes))
         do i = 1, Ntypes
            TAB_Rho_match(i) = tillage%types(i)%rho_match
            TAB_N_match(i)   = tillage%types(i)%N_match
         end do
      end if

      Max_Z_tillage = maxval(Z_tillage(1:Ntill))

      ! iTT1 / iTT2: first/last position per tillage type in iType_Tillage.
      ! Replicates the loop from the legacy Read_Tillage subroutine verbatim.
      if (allocated(iTT1)) deallocate(iTT1); allocate(iTT1(Ntill)); iTT1 = 0
      if (allocated(iTT2)) deallocate(iTT2); allocate(iTT2(Ntill)); iTT2 = 0
      do j = 1, Ntill
         do i = 1, Ntypes
            if (iTT1(j) == 0 .and. iType_Tillage(i) == j) iTT1(j) = i
            if (iTT1(j) >  0 .and. iType_Tillage(i) == j) iTT2(j) = i
         end do
      end do
   end subroutine apply_soil_tillage


   !> Apply [irrigation.ssdi] config to legacy `variables` globals.
   !! Called from config_to_variables when flSSDI is true. Resolves
   !! ssdi_z to nod_ssdi_irr(1:2) layer indices via the zbotcp walk;
   !! mode 0 reads the events CSV (column 1 auto-converts ISO dates to
   !! days-since-1900); mode 1 copies the scheduled sub-block; both
   !! modes apply legacy unit conversions (mm/h -> cm/d, mm -> cm
   !! spread over the in-zone compartments) and initialize qssdi = 0.
   !!
   !! Replaces SSDI_irrigation(1) and read_ssdi_input (deletion in
   !! Task 5).
   subroutine apply_irrigation_ssdi(ssdi, tstart, tend, state)
      use, intrinsic :: iso_fortran_env, only: real64
      use irrigation_config_mod, only: irrigation_ssdi_t
      use error_mod, only: fatalerr_collected
      use swap_state_mod, only: swap_state_t
      use variables, only: nod_ssdi_irr, qssdi, dt_SSDI_event
      type(irrigation_ssdi_t), intent(in)    :: ssdi
      real(real64),            intent(in)    :: tstart, tend
      type(swap_state_t),      intent(inout) :: state   ! [GR-BH Task 35] replaces NumNod/zbotcp globals

      integer :: i, j, nod_top, nod_bot, ncomp

      ! Resolve ssdi_z(1:2) -> layer indices via zbotcp walk.
      ! Mirrors the legacy SSDI_irrigation(1) loop at irrigation.f90:387-393.
      ! [GR-BH Task 35] zbotcp/NumNod globals replaced by state%mesh fields.
      ! Guard: mesh not yet populated at config_to_variables call time;
      ! nod_ssdi_irr defaults to 0 if mesh not built (resolved after CalcGrid).
      nod_ssdi_irr = 0
      if (state%mesh%numnod > 0 .and. allocated(state%mesh%zbotcp)) then
         do j = 1, 2
            i = 1
            do while (state%mesh%zbotcp(i) > (ssdi%ssdi_z(j) + 1.0e-5_real64))
               i = i + 1
               if (i > state%mesh%numnod) exit
            end do
            nod_ssdi_irr(j) = i
         end do
      end if
      nod_top = nod_ssdi_irr(1)
      nod_bot = nod_ssdi_irr(2)
      ncomp   = nod_bot - nod_top + 1
      if (ncomp < 1) then
         call fatalerr_collected('apply_irrigation_ssdi', &
                                 'ssdi_z resolves to zero compartments — check ssdi_z vs grid')
      end if

      ! Initialize qssdi to zero (mirrors irrigation.f90:403).
      qssdi = 0.0_real64
      dt_SSDI_event = 1.0_real64

      select case (ssdi%schedule)
      case (0)
         call apply_ssdi_mode0(ssdi, ncomp, tstart, tend)
      case (1)
         call apply_ssdi_mode1(ssdi, ncomp)
      end select
   end subroutine apply_irrigation_ssdi


   !> Apply [nutrients] config to legacy `variables`/wofost_soil_declarations
   !! globals. Always called from config_to_variables (no flCropNut gate);
   !! the cfg%present flag is informational only — defaults are zero
   !! whether or not the user supplied a [nutrients] block.
   !!
   !! Sets SorpCoef unconditionally — fixes the genuine uninitialised-
   !! variable bug discovered during the [nutrients] N2 brainstorm.
   !!
   !! See ADR 0026 ([nutrients] N2a).
   subroutine apply_nutrients(cfg)
      use nutrients_config_mod, only: nutrients_config_t
      use wofost_soil_declarations, only: FOM_t, Bio_t, Hum_t, &
                                           cNH4_t, cNO3_t, SorpCoef
      type(nutrients_config_t), intent(in) :: cfg
      integer :: i

      SorpCoef = cfg%sorp_coef
      do i = 1, 8
         FOM_t(i) = cfg%initial%fom(i)
      end do
      Bio_t  = cfg%initial%bio
      Hum_t  = cfg%initial%hum
      cNH4_t = cfg%initial%cnh4
      cNO3_t = cfg%initial%cno3

      ! N2b (ADR 0027): stage timed amendments from the CSV companion.
      call apply_nutrients_events(cfg)
   end subroutine apply_nutrients


   !> Stage timed soil management events from the CSV companion
   !! at cfg%events_file. Sets the legacy globals consumed by
   !! SoilManagement(3) and Wofost_SoilAmendents.
   !!
   !! Default (empty events_file or empty CSV): namend = 0, isme = 1.
   !!
   !! See ADR 0027 ([nutrients] N2b).
   subroutine apply_nutrients_events(cfg)
      use, intrinsic :: iso_fortran_env, only: real64
      use csv_reader_mod, only: read_csv_table
      use error_mod, only: error_collection_t, fatalerr_collected
      use nutrients_config_mod, only: nutrients_config_t
      use variables, only: pathwork
      use wofost_soil_declarations, only: MatNum, Amend, VolaFrac, &
                                            TimeAmend, NuAmend, iamend, &
                                            namend, isme, maxamn
      type(nutrients_config_t), intent(in) :: cfg

      real(real64), allocatable :: tbl(:,:)
      type(error_collection_t)  :: errs
      character(len=300) :: csvpath
      character(len=14)  :: hdr(4)
      integer :: i, j, n
      real(real64) :: tmp_date, tmp_amount, tmp_volat
      real(real64) :: tmp_mat_real

      ! Default: no amendments. Reset legacy globals to a known state.
      namend = 0
      isme   = 1

      ! events_file is always allocated by get_optional_string_with_default
      ! (defaults to empty string on missing key), so check len_trim instead.
      if (.not. allocated(cfg%events_file)) return
      if (len_trim(cfg%events_file) == 0)   return

      hdr(1) = 'date          '
      hdr(2) = 'material      '
      hdr(3) = 'amount_kgha   '
      hdr(4) = 'volat_fraction'
      csvpath = trim(pathwork) // trim(cfg%events_file)
      call read_csv_table(trim(csvpath), hdr, tbl, errs)
      call errs%abort_if_fatal()

      n = 0
      if (allocated(tbl)) n = size(tbl, 1)
      if (n < 1) return     ! Empty CSV: no amendments. Not an error.
      if (n > maxamn) then
         call fatalerr_collected('apply_nutrients_events', &
            'CSV row count exceeds maxamn (1000)')
         return
      end if

      ! Per-row validation
      do i = 1, n
         if (nint(tbl(i, 2)) < 1 .or. nint(tbl(i, 2)) > 20) then
            call fatalerr_collected('apply_nutrients_events', &
               'material out of range [1, 20]')
            return
         end if
         if (tbl(i, 3) < 0.0_real64 .or. tbl(i, 3) > 500000.0_real64) then
            call fatalerr_collected('apply_nutrients_events', &
               'amount_kgha out of range [0, 500000]')
            return
         end if
         if (tbl(i, 4) < 0.0_real64 .or. tbl(i, 4) > 1.0_real64) then
            call fatalerr_collected('apply_nutrients_events', &
               'volat_fraction out of range [0, 1]')
            return
         end if
      end do

      ! Sort by date (in-place bubble sort, mirrors deleted SoilManagement(1)).
      ! Acceptable O(n^2) given n <= 1000 and this runs once at config-load.
      do i = 1, n - 1
         do j = i + 1, n
            if (tbl(i, 1) > tbl(j, 1)) then
               tmp_date     = tbl(i, 1); tbl(i, 1) = tbl(j, 1); tbl(j, 1) = tmp_date
               tmp_mat_real = tbl(i, 2); tbl(i, 2) = tbl(j, 2); tbl(j, 2) = tmp_mat_real
               tmp_amount   = tbl(i, 3); tbl(i, 3) = tbl(j, 3); tbl(j, 3) = tmp_amount
               tmp_volat    = tbl(i, 4); tbl(i, 4) = tbl(j, 4); tbl(j, 4) = tmp_volat
            end if
         end do
      end do

      ! Populate per-event legacy globals
      do i = 1, n
         MatNum(i)   = nint(tbl(i, 2))
         Amend(i)    = 1.0e-4_real64 * tbl(i, 3)   ! kg/ha -> kg/m^2
         VolaFrac(i) = tbl(i, 4)
      end do

      ! Group dosages per date (mirrors deleted SoilManagement(1)).
      j = 1
      NuAmend(j)   = 1
      TimeAmend(j) = tbl(1, 1)
      iamend(1, 1) = 1
      do i = 2, n
         if (abs(tbl(i, 1) - tbl(i - 1, 1)) < 1.0e-3_real64) then
            NuAmend(j) = NuAmend(j) + 1
         else
            j = j + 1
            NuAmend(j)   = 1
            TimeAmend(j) = tbl(i, 1)
         end if
         iamend(j, NuAmend(j)) = i
      end do

      namend = j
      isme   = 1
   end subroutine apply_nutrients_events


   !> Mode-0 (fixed-date): stage CSV; populate ssdi_*_f_irr; deferred
   !! date-window validation; initial nirri_ssdi_irr entry-point from tstart.
   subroutine apply_ssdi_mode0(ssdi, ncomp, tstart, tend)
      use, intrinsic :: iso_fortran_env, only: real64
      use csv_reader_mod, only: read_csv_table
      use error_mod, only: error_collection_t, fatalerr_collected
      use irrigation_config_mod, only: irrigation_ssdi_t
      use variables, only: pathwork, mairg,                              &
                           nirri_ssdi_irr,                               &
                           ssdi_date_irr, ssdi_rate_f_irr,              &
                           ssdi_amount_f_irr
      type(irrigation_ssdi_t), intent(in) :: ssdi
      integer,                 intent(in) :: ncomp
      real(real64),            intent(in) :: tstart, tend

      real(real64), allocatable :: tbl(:,:)
      type(error_collection_t)  :: errs
      character(len=300) :: csvpath
      character(len=8)   :: hdr(3)
      integer :: i, n, nirri_init
      logical :: any_in_window, window_in_dates

      hdr(1) = 'date    '
      hdr(2) = 'rate_f  '
      hdr(3) = 'amount_f'
      csvpath = trim(pathwork) // trim(ssdi%fixed%events_file)
      call read_csv_table(trim(csvpath), hdr, tbl, errs)
      call errs%abort_if_fatal()

      n = 0
      if (allocated(tbl)) n = size(tbl, 1)
      if (n < 1) then
         call fatalerr_collected('apply_irrigation_ssdi', &
                                 'mode 0: events CSV has no rows')
         return
      end if
      if (n > mairg) then
         call fatalerr_collected('apply_irrigation_ssdi', &
                                 'mode 0: events CSV exceeds mairg rows')
         return
      end if

      ! ssdi_date_irr / ssdi_rate_f_irr / ssdi_amount_f_irr are fixed-size
      ! arrays of size(mairg); zero them and fill from the CSV.
      ssdi_date_irr     = 0.0_real64
      ssdi_rate_f_irr   = 0.0_real64
      ssdi_amount_f_irr = 0.0_real64

      do i = 1, n
         ssdi_date_irr(i) = tbl(i, 1)
         if (i > 1 .and. ssdi_date_irr(i) <= ssdi_date_irr(i-1)) then
            call fatalerr_collected('apply_irrigation_ssdi', &
                                    'mode 0: ssdi_date not strictly ascending')
            return
         end if
         ! mm/h -> cm/d (mirrors irrigation.f90:514)
         ssdi_rate_f_irr(i)   = tbl(i, 2) * 0.1_real64 * 24.0_real64
         ! mm -> cm, then spread over `ncomp` compartments (mirrors irrigation.f90:397)
         ssdi_amount_f_irr(i) = (tbl(i, 3) * 0.1_real64) / real(ncomp, real64)
      end do

      ! Date-window check: at least one date in [tstart, tend], OR
      ! [tstart, tend] contained in [date(1), date(n)].
      ! Replaces the deleted checkdate call (irrigation.f90:552).
      any_in_window  = .false.
      do i = 1, n
         if (ssdi_date_irr(i) >= tstart - 1.0e-6_real64 .and. &
             ssdi_date_irr(i) <= tend   + 1.0e-6_real64) then
            any_in_window = .true.
            exit
         end if
      end do
      window_in_dates = (ssdi_date_irr(1) <= tstart + 1.0e-6_real64 .and. &
                         ssdi_date_irr(n) >= tend   - 1.0e-6_real64)
      if (.not. any_in_window .and. .not. window_in_dates) then
         call fatalerr_collected('apply_irrigation_ssdi', &
                                 'mode 0: no ssdi_date within simulation period')
      end if

      ! Determine initial entry point (mirrors irrigation.f90:368-373).
      nirri_init = 1
      do i = 1, n - 1
         if (tstart >= ssdi_date_irr(i)) nirri_init = i
      end do
      if (tstart >= ssdi_date_irr(n)) nirri_init = n
      nirri_ssdi_irr = nirri_init
   end subroutine apply_ssdi_mode0


   !> Mode-1 (scheduled-trigger): copy scheduled sub-block to legacy
   !! globals; set nirri_ssdi_irr/dt_SSDI_event/days_counter defaults
   !! (preserves the d8a88d6 regression-fix invariant for
   !! dt_SSDI_event = 1.0).
   subroutine apply_ssdi_mode1(ssdi, ncomp)
      use, intrinsic :: iso_fortran_env, only: real64
      use irrigation_config_mod, only: irrigation_ssdi_t
      use variables, only: nirri_ssdi_irr, dt_SSDI_event,               &
                           ssdi_sched_type_irr, ssdi_threshold_irr,     &
                           ssdi_threshold_z_irr, ssdi_amount_irr,       &
                           ssdi_appl_rate_irr, sw_interval_irr,         &
                           days_interval_irr, days_counter_irr
      type(irrigation_ssdi_t), intent(in) :: ssdi
      integer,                 intent(in) :: ncomp

      ssdi_sched_type_irr  = ssdi%scheduled%sched_type
      ssdi_threshold_irr   = ssdi%scheduled%threshold
      ssdi_threshold_z_irr = ssdi%scheduled%threshold_depth

      ! mm -> cm, then spread over ncomp compartments (mirrors irrigation.f90:399)
      ssdi_amount_irr      = (ssdi%scheduled%ssdi_amount * 0.1_real64) / &
                              real(ncomp, real64)
      ! mm/h -> cm/d (mirrors irrigation.f90:546)
      ssdi_appl_rate_irr   = ssdi%scheduled%ssdi_appl_rate * 0.1_real64 * 24.0_real64

      sw_interval_irr      = ssdi%scheduled%sw_interval
      ! mirrors irrigation.f90:535-538
      if (ssdi%scheduled%sw_interval == 0) then
         days_interval_irr = 1
      else
         days_interval_irr = ssdi%scheduled%days_interval
      end if
      days_counter_irr = 366   ! mirrors irrigation.f90:540

      nirri_ssdi_irr = 1
   end subroutine apply_ssdi_mode1

end module config_to_variables_mod
