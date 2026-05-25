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
   ! [GR-SEED 2026-05-25 Task 4] apply_nutrients + apply_nutrients_events removed:
   ! bodies moved to nutrients_state_mod as seed_nutrients_from_config +
   ! load_nutrients_events; called from state%nutrients%init in swap_mod.
   ! [GR-SEED 2026-05-25 Task 5] apply_irrigation_ssdi + apply_ssdi_mode0 +
   ! apply_ssdi_mode1 removed: bodies moved to crop_irrigation_state_mod as
   ! apply_ssdi_seed (public) + apply_ssdi_mode0/mode1 (private); called from
   ! state%crop%irrigation%init in swap_mod.

contains

   !> Copy every (C)-classified field from `config` into the corresponding
   !! `variables%` legacy global. Caller is responsible for having loaded,
   !! validated, and finalized `config` first.
   subroutine config_to_variables(config, state)
      ! [GR-IO 2026-05-25 Phase 6 Step 3] bare `use variables` retired —
      ! every config→bare-global mirror write below was dead after Steps
      ! 1-2 (no remaining readers). Size parameters now come from the
      ! canonical swap_array_dimensions module; everything else writes
      ! directly to state%X / config%X.
      use swap_array_dimensions, only: macp, madr, maho, mamp, mabbc
      use swap_state_mod, only: swap_state_t
      use error_mod, only: fatalerr_collected
      type(swap_config_t), intent(inout), target :: config  ! [GR-SOIL 2026-05-24] inout: populates config%soil%initial%z_init from h_file CSV
      type(swap_state_t),  intent(inout)         :: state

      integer :: i, n

      ! ---------------------------------------------------------------
      ! General + simulation + numerical (audit: 18 fields)
      ! [GR-SEED 2026-05-25] Task 1 — absorbed by state%timecontrol%init.
      ! ---------------------------------------------------------------
      call state%timecontrol%init(config%simulation, config%general)

      ! Meteorology (audit: 12 + evaporation + snow)
      ! ---------------------------------------------------------------
      ! [GR-SEED 2026-05-25 Task 2] Meteorology CSV pre-loads + snow scalars
      ! moved to state%atmosphere%init(config) (called from swap_mod).

      ! Evaporation sub-section
      ! [GR-IO 2026-05-25 Phase 6 Step 3] swcfbs legacy mirror dropped
      ! [SS-T10] Deferred — will move into state%crop%init
      state%crop%cfbs = config%meteo%evaporation%cfbs
      ! [GR-ATM 2026-05-23] swredu/cofred/rsigni/cfevappond retired —
      ! snapshotted into state%atmosphere by atmosphere_state%init(config);
      ! compute reads from state, never from these legacy globals.

      ! ---------------------------------------------------------------
      ! Drainage (audit: 20 fields + surface_runoff sub-section)
      ! ---------------------------------------------------------------
      state%surfacewater%swdra = config%drain%swdra
      state%drainage%dramet    = config%drain%dramet
      ! [GR-BH Task 37] swdivd global deleted — state%drainage%swdivd seeded in swap_mod.f90
      state%drainage%swdislay  = config%drain%swdislay
      ! [GR-BH Task 37] nrlevs global deleted — state%drainage%nrlevs seeded in swap_mod.f90
      state%drainage%basegw = config%drain%basegw
      state%drainage%entres = config%drain%entres

      ! DRAMET=2 (Hooghoudt/Ernst). Mirrors readswap.f90:1850-1875.
      ! lm is authored in metres; the legacy reader does the m->cm
      ! conversion (`l(1) = 100*lm2`) so we replicate that here.
      ! ADR 0031 Phase 2 Task 5: wetper(1) removed — state%drainage%wetper(1)
      ! is seeded from config%drain%wetper in drainage_init instead.
      ! zbotdr goes into the level-1 entry of the per-level array;
      ! ipos / khtop / khbot / kvtop / kvbot / zintf / geofac are scalar globals.
      if (config%drain%dramet == 2) then
         ! [GR-BH Task 37] L(1)/zbotdr(1) bare globals deleted — seeded from config in swap_mod.f90
         state%drainage%ipos    = config%drain%ipos
         state%drainage%khtop   = config%drain%khtop
         if (config%drain%ipos >= 3) then
            state%drainage%khbot = config%drain%khbot
            state%drainage%zintf = config%drain%zintf
         end if
         if (config%drain%ipos >= 4) then
            state%drainage%kvtop = config%drain%kvtop
            state%drainage%kvbot = config%drain%kvbot
         end if
         if (config%drain%ipos == 5) then
            state%drainage%geofac = config%drain%geofac
         end if
      end if

      ! [GR-BH Task 36] cofani global retired — precedence logic moved to swap_mod.f90
      ! config%drain%cofani is consumed directly by swap_mod seeding block.

      if (allocated(config%drain%swdtyp)) then
         if (.not. allocated(state%drainage%swdtyp)) then
            allocate(state%drainage%swdtyp(size(config%drain%swdtyp)))
            state%drainage%swdtyp = 0
         end if
         do i = 1, size(config%drain%swdtyp)
            state%drainage%swdtyp(i) = config%drain%swdtyp(i)
         end do
      end if
      ! [GR-BH Task 37] zbotdr bare global deleted — seeded from config%drain in swap_mod.f90
      ! [GR-DRAIN 2026-05-25] Per-level drainage arrays (drares, infres, gwlinf,
      ! rdrain, rinfi, rentry, rexit, widthr, taludr) are read directly from
      ! config%drain in src/drainage/drainage.f90; no mirror writes needed.
      if (allocated(config%drain%swallo)) then
         if (.not. allocated(state%drainage%swallo)) then
            allocate(state%drainage%swallo(size(config%drain%swallo)))
            state%drainage%swallo = 0
         end if
         do i = 1, size(config%drain%swallo)
            state%drainage%swallo(i) = config%drain%swallo(i)
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
            use swap_array_dimensions, only: MAOWL
            real(real64), allocatable :: csv_table(:,:)
            type(error_collection_t)  :: csv_errs
            integer :: lev, nrows, k
            character(len=6) :: hdr(2)
            hdr(1) = 'date  '
            hdr(2) = 'level '
            ! [GR-IO 2026-05-25 Phase 5] write directly to state%drainage%owltab —
            ! bare-global `owltab` staging buffer retired. drainage_init allocates
            ! state%drainage%owltab later in startup, but we need it now; pre-allocate
            ! here on the first per-level CSV load.
            if (.not. allocated(state%drainage%owltab)) then
               allocate(state%drainage%owltab(config%drain%nrlevs, 2*MAOWL))
               state%drainage%owltab = 0.0_real64
            end if
            do lev = 1, size(config%drain%owltab_file)
               if (len_trim(config%drain%owltab_file(lev)) == 0) cycle
               call read_csv_table(trim(config%drain%owltab_file(lev)), hdr, csv_table, csv_errs)
               call csv_errs%abort_if_fatal()
               if (csv_errs%count() == 0) then
                  nrows = size(csv_table, 1)
                  state%drainage%nowltab(lev) = nrows
                  do k = 1, nrows
                     state%drainage%owltab(lev, 2*k-1) = csv_table(k, 1)  ! date (days since 1900)
                     state%drainage%owltab(lev, 2*k)   = csv_table(k, 2)  ! channel water level (cm)
                  end do
               end if
            end do
         end block
      end if

      state%drainage%swliminf = config%drain%swliminf

      ! Drainage.surface_runoff sub-section: scalar switches + per-level
      ! arrays. Legacy globals `swtopdislay`, `ftopdislay`, `RapDraResRef`
      ! are arrays of size madr; copy element 1 of the (scalar) config
      ! field as a uniform value across drainage levels until a per-level
      ! schema lands in Phase 4f-extend.
      ! [GR-BH Task 37] swnrsrf/SwTopnrsrf/swdivdinf/FacDpthInf bare globals deleted —
      ! state%drainage%X seeded in swap_mod.f90
      state%drainage%cofintfl = config%drain%surface_runoff%cofintfl
      state%drainage%expintfl = config%drain%surface_runoff%expintfl
      ! ADR 0031: gate the surface_runoff geofac write to avoid overwriting
      ! the ipos==5 (Ernst geometry factor) write at line 306. Two distinct
      ! TOML fields map to one legacy global; the gate preserves both
      ! intended behaviors. Schema-level reconciliation deferred.
      if (config%drain%ipos /= 5) then
         state%drainage%geofac = config%drain%surface_runoff%geofac
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
      ! [GR-DRAIN 2026-05-25] rsurfdeep/rsurfshallow legacy mirror writes
      ! dropped — read directly from config%drain%surface_runoff in
      ! src/drainage/drainage.f90.
      ! [SS-GR-FINAL D1] RapDraReaExp write dropped — global retired
      state%drainage%NumLevRapDra = config%drain%surface_runoff%numlevrapdra
      ! swtopdislay(MADR) and ftopdislay(MADR): broadcast scalar config
      ! field to all drain levels (currently no per-level schema slot).
      ! Allocate to MADR to match legacy fixed-size globals.
      if (.not. allocated(state%drainage%swtopdislay)) then
         allocate(state%drainage%swtopdislay(madr))
         state%drainage%swtopdislay = 0
      end if
      if (.not. allocated(state%drainage%ftopdislay)) then
         allocate(state%drainage%ftopdislay(madr))
         state%drainage%ftopdislay = 0.0d0
      end if
      do i = 1, size(state%drainage%swtopdislay)
         state%drainage%swtopdislay(i) = config%drain%surface_runoff%swtopdislay
      end do
      do i = 1, size(state%drainage%ftopdislay)
         state%drainage%ftopdislay(i) = config%drain%surface_runoff%ftopdislay
      end do
      ! [SS-GR-FINAL D1] RapDraResRef write dropped — global retired

      ! ---------------------------------------------------------------
      ! Soil (audit: 15 + discretization + frost)
      ! ---------------------------------------------------------------
      state%soilwater%swsophy = config%soil%swsophy
      ! [GR-SOIL 2026-05-24] swhyst legacy mirror dropped — direct config read.
      state%soilwater%swinco = config%soil%swinco
      ! [MACRO-RETIRE 2026-05-12] swmacro global retired (ADR 0040).
      ! soil.swmacro=1 is still rejected by soil_config validator stub.
      ! SS-B / ADR 0020: set the call-site gating flags. When flTillage /
      ! flSSDI are false (default for every regression case), the call
      ! sites in swap.f90 / timecontrol.f90 short-circuit and the
      ! subsystems never run.
      ! [GR-CROP 2026-05-25] till_swtill mirror retired — tillage reads
      ! state%cfg%soil%swtill directly.
      ! [GR-IO 2026-05-25 Phase 5] flTillage bare global retired — swap_mod and
      ! this adapter now read (config%soil%swtill == 1) directly.
      ! [GR-TIME 2026-05-25] flSSDI bare global retired — swap_mod and
      ! timecontrol_mod now read (config%irrigation%swssdi == 1) directly.
      if (config%soil%swtill == 1) call apply_soil_tillage(config%soil%tillage, state%timecontrol%tend, state)
      ! [GR-SEED 2026-05-25 Task 5] Irrigation seeding (swssdi mirror + fixed events + SSDI)
      ! moved to state%crop%irrigation%init (called from swap_mod after CalcGrid).
      ! [GR-SEED 2026-05-25 Task 4] apply_nutrients moved to state%nutrients%init
      ! (called from swap_mod). Body now in nutrients_state_mod as public
      ! seed_nutrients_from_config + load_nutrients_events helpers.
      ! [GR-SOIL 2026-05-24] gwli legacy mirror dropped — direct config read.
      ! [GR-FINAL C1] pondini/pond: config%soil%pondini read directly by swap_mod after soilwater_init
      ! (pondini_init_buf/pond_init_buf retired; swap_mod seeding replaced with direct config reads)
      state%surfacewater%pondmx = config%soil%pondmx
      ! [GR-ATM 2026-05-23] rsoil retired — snapshotted in atmosphere_state%init
      state%surfacewater%rsro = config%soil%rsro
      state%surfacewater%rsroexp = config%soil%rsroexp
      ! Legacy parses .swp `SWRUNON` into a local int; we mirror that mapping
      ! into state%soilwater%flrunon (runonarr remains dormant — no TOML writer).
      state%soilwater%flrunon = (config%soil%swrunon == 1)
      ! [GR-IO 2026-05-25 Phase 6 Step 3] nrstaring legacy mirror dropped

      ! sublay (legacy 'isublay') is a local in readswap, not a module
      ! global; calcgrid only consumes nsublay + isoillay + ncomp + hcomp
      ! [GR-SOIL 2026-05-24] sublay/isoillay/hsublay/ncomp/hcomp bare-global writes
      ! retired — CalcGrid reads them inline from config%soil%X.
      if (allocated(config%soil%isoillay)) then
         state%mesh%numlay = config%soil%isoillay(size(config%soil%isoillay))
      end if
      ! [GR-SOIL 2026-05-24] config%soil%hcomp ingest retired — CalcGrid reads inline.
      ! [GR-BH Task 36] orgmat global retired — seeded via state%soilwater%orgmat in swap_mod.f90
      ! config%soil%orgmat is consumed directly by swap_mod seeding block.
      if (allocated(config%soil%bdens)) then
         ! state%soilwater%bdens — Pattern 9 guarded alloc (soilwater_init runs later).
         if (.not. allocated(state%soilwater%bdens)) then
            allocate(state%soilwater%bdens(maho)); state%soilwater%bdens = 0.0d0
         end if
         do i = 1, size(config%soil%bdens)
            state%soilwater%bdens(i) = config%soil%bdens(i)
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
            ! [GR-IO 2026-05-25 Phase 6 Step 3] atmin7 → state%atmosphere directly
            state%atmosphere%atmin7(:) = config%soil%initial%atmin7(:)
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
               ! [GR-SOIL 2026-05-24] z_init in config; legacy zi/nhead mirror retained for compat.
               if (allocated(config%soil%initial%z_init)) deallocate(config%soil%initial%z_init)
               allocate(config%soil%initial%z_init(nrows))
               do k = 1, nrows
                  config%soil%initial%z_init(k) = tbl(k, 1)
               end do
            end block

            ! [GR-IO 2026-05-25 Phase 6 Step 3] Legacy tsoil_file CSV → zh/tsoil
            ! bare-global block dropped. temperature.f90:case(1) reads
            ! cfg_heat%tsoil_init directly (populated by read_heat_toml).

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
                  state%solute%nconc = nrows
                  if (.not. allocated(state%solute%cml_init)) then
                     allocate(state%solute%cml_init(macp)); state%solute%cml_init = 0.0d0
                  end if
                  if (.not. allocated(state%solute%zc_init)) then
                     allocate(state%solute%zc_init(macp));  state%solute%zc_init  = 0.0d0
                  end if
                  do k = 1, nrows
                     state%solute%zc_init(k)  = tbl(k, 1)
                     state%solute%cml_init(k) = tbl(k, 2)
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
         if (.not. allocated(state%soilwater%bdens)) then
            allocate(state%soilwater%bdens(maho)); state%soilwater%bdens = 0.0d0
         end if
         ! [GR-SOIL 2026-05-24] h_enpr legacy mirror dropped — vg_params carries the typed value.
         do i = 1, size(config%soil%hydraulics%ores)
            state%soilwater%bdens(i) = config%soil%hydraulics%bdens(i)
         end do

         ! [GR-SOIL 2026-05-24] iHWCKmodel legacy write retired — soilwater_init
         !   seeds `sw%iHWCKmodel = 1` directly (HACK Phase 4f-extend constraint).
         !   The legacy reader at readswap.f90:651-657 supported per-layer override
         !   (1..11) but no TOML schema slot covers it yet.

         ! [GR-CROP 2026-05-25] paramvg legacy mirror retired — tillage.f90 now
         !   mutates the typed per-layer store state%soilwater%vg_params_layer(:),
         !   populated by SoilHydraulics(1) from state%cfg%soil%hydraulics directly.
      end if

      ! Soil.discretization
      ! [GR-IO 2026-05-25 Phase 6 Step 3] swdiscrvert legacy mirror dropped
      ! [GR-IO 2026-05-25] numnodnew/dznew legacy mirror dropped — swapoutput.f90:checkDiscrVert
      ! reads config%soil%discretization%{numnodnew,dznew} directly via state%cfg.

      ! Soil.frost
      ! [GR-IO 2026-05-25 Phase 6 Step 3] swfrost legacy mirror dropped — dual-write
      ! collapsed to single state-side write.
      state%soilwater%swfrost = config%soil%frost%swfrost
      ! [GR-IO 2026-05-25 Phase 6 Step 3] swsublim legacy mirror dropped
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
               trim(config%general%pathwork)//trim(config%bottom_boundary%gwl_file), &
               hdr, csv_table, csv_errs)
            call csv_errs%abort_if_fatal()
            if (allocated(csv_table)) then
               if (.not. allocated(state%soilwater%gwltab)) then
                  allocate(state%soilwater%gwltab(2*mabbc))
                  state%soilwater%gwltab = 0.0d0
               end if
               nrows = size(csv_table, 1)
               do k = 1, nrows
                  state%soilwater%gwltab(k*2 - 1) = csv_table(k, 1)
                  state%soilwater%gwltab(k*2)     = csv_table(k, 2)
               end do
            end if
         end block
      case (2)
         ! [GR-IO 2026-05-25 Phase 6 Step 3] sw2 legacy mirror dropped
         ! Phase 0 B-0.1: populate sine-wave scalars regardless of sw2;
         ! the gate in boundbottom.f90:104 protects the non-sine path.
         ! [GR-IO 2026-05-25 Phase 6 Step 3] sinmax/sinamp/sinave legacy mirrors dropped
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
                  trim(config%general%pathwork)//trim(config%bottom_boundary%qbot2_file), &
                  hdr, csv_table, csv_errs)
               call csv_errs%abort_if_fatal()
               if (allocated(csv_table)) then
                  if (.not. allocated(state%soilwater%qbotab)) then
                     allocate(state%soilwater%qbotab(2*mabbc))
                     state%soilwater%qbotab = 0.0d0
                  end if
                  nrows = size(csv_table, 1)
                  do k = 1, nrows
                     state%soilwater%qbotab(k*2 - 1) = csv_table(k, 1)
                     state%soilwater%qbotab(k*2)     = csv_table(k, 2)
                  end do
               end if
            end block
         end if
      case (3)
         ! [GR-DRAIN 2026-05-25] shape legacy mirror dropped — boundbottom.f90
         ! reads bb%shape (= config%bottom_boundary%shape) directly.
         ! [GR-IO 2026-05-25 Phase 6 Step 3] hdrain/aqave/aqamp/aqper/aqtmax/sw3
         ! legacy mirrors dropped (boundbottom.f90 reads bb%X = config%bottom_boundary%X).
         ! [GR-SOIL 2026-05-24] rimlay/swbotb3impl already dropped.
         ! [GR-SOIL 2026-05-24] sw4 legacy mirror dropped — direct config read.
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
                  trim(config%general%pathwork)//trim(config%bottom_boundary%haquif_file), &
                  hdr, csv_table, csv_errs)
               call csv_errs%abort_if_fatal()
               if (allocated(csv_table)) then
                  if (.not. allocated(state%soilwater%haqtab)) then
                     allocate(state%soilwater%haqtab(2*mabbc))
                     state%soilwater%haqtab = 0.0d0
                  end if
                  nrows = size(csv_table, 1)
                  do k = 1, nrows
                     state%soilwater%haqtab(k*2 - 1) = csv_table(k, 1)
                     state%soilwater%haqtab(k*2)     = csv_table(k, 2)
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
                  trim(config%general%pathwork)//trim(config%bottom_boundary%qbot4_file), &
                  hdr, csv_table, csv_errs)
               call csv_errs%abort_if_fatal()
               if (allocated(csv_table)) then
                  if (.not. allocated(state%soilwater%qbotab)) then
                     allocate(state%soilwater%qbotab(2*mabbc))
                     state%soilwater%qbotab = 0.0d0
                  end if
                  nrows = size(csv_table, 1)
                  do k = 1, nrows
                     state%soilwater%qbotab(k*2 - 1) = csv_table(k, 1)
                     state%soilwater%qbotab(k*2)     = csv_table(k, 2)
                  end do
               end if
            end block
         end if
      case (4)
         ! [GR-IO 2026-05-25 Phase 6 Step 3] swqhbot/cofqha/cofqhb/cofqhc/swcofqhc
         ! legacy mirrors dropped — boundbottom reads via bb%X.
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
                  trim(config%general%pathwork)//trim(config%bottom_boundary%qhbot_file), &
                  hdr, csv_table, csv_errs)
               call csv_errs%abort_if_fatal()
               ! Legacy unpack pattern from readswap.f90:1418-1419 — for the
               ! q(h) curve, qbotab(odd) = abs(htab) and qbotab(even) = qtab.
               if (allocated(csv_table)) then
                  if (.not. allocated(state%soilwater%qbotab)) then
                     allocate(state%soilwater%qbotab(2*mabbc))
                     state%soilwater%qbotab = 0.0d0
                  end if
                  nrows = size(csv_table, 1)
                  do k = 1, nrows
                     state%soilwater%qbotab(k*2 - 1) = abs(csv_table(k, 1))
                     state%soilwater%qbotab(k*2)     = csv_table(k, 2)
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
               trim(config%general%pathwork)//trim(config%bottom_boundary%hbot5_file), &
               hdr, csv_table, csv_errs)
            call csv_errs%abort_if_fatal()
            if (allocated(csv_table)) then
               if (.not. allocated(state%soilwater%hbotab)) then
                  allocate(state%soilwater%hbotab(2*mabbc))
                  state%soilwater%hbotab = 0.0d0
               end if
               nrows = size(csv_table, 1)
               do k = 1, nrows
                  state%soilwater%hbotab(k*2 - 1) = csv_table(k, 1)
                  state%soilwater%hbotab(k*2)     = csv_table(k, 2)
               end do
            end if
         end block
      case (6, 7)
         ! No parameters to populate for modes 6 and 7.
      case (8)
         ! [GR-SOIL 2026-05-24] hplate legacy mirror dropped — direct config read.
      end select

      ! ---------------------------------------------------------------
      ! Heat (audit: 10 fields + 6 Phase-0 promoted fields = 16 total)
      ! ---------------------------------------------------------------
      ! [GR-TIME 2026-05-25] swhea legacy mirror dropped — timecontrol_mod reads
      ! config%heat%swhea directly via state%cfg.
      ! [GR-IO 2026-05-25 Phase 6 Step 3] heat scalars legacy mirrors dropped —
      ! temperature.f90 reads cfg_heat%X (= config%heat%X) directly.

      ! psand/psilt/pclay: legacy global writes retired — state%soilwater%psand/psilt/pclay
      ! now sourced directly from config%heat in swap_mod.f90 [SS-GR-BH A6].
      ! [GR-BH Task 36] orgmat global retired — heat.porg→orgmat backfill now handled in
      ! swap_mod.f90 seeding block (state%soilwater%orgmat). config%heat%porg consumed directly.
      ! [GR-IO 2026-05-25 Phase 6 Step 3] tsoil_init → zh/tsoil mirror retired —
      ! temperature.f90:case(1) builds the afgen depth-temp table directly
      ! from cfg_heat%tsoil_init (= config%heat%tsoil_init). Same for
      ! temtoptab/tembtab — temperature.f90 reads cfg_heat%X straight.

      ! Phase 0 (SS-HEAT) — swcalt=1 analytical method scalars.
      ! [GR-IO 2026-05-25 Phase 6 Step 3] ddamp/tmean/tampli/timref legacy
      ! mirrors dropped — temperature.f90 reads cfg_heat%X directly.

      ! ---------------------------------------------------------------
      ! [GR-SEED 2026-05-25 Task 5] Irrigation seeding (fixed events + SSDI)
      ! moved to state%crop%irrigation%init (called from swap_mod after CalcGrid).
      ! ---------------------------------------------------------------
      ! cirrs / cirrthres / dcrit / isuas / perirrsurp / raithreshold /
      ! swcirrthres live in irrigation_schedule_t (per-crop), not the
      ! top-level irrigation_config_t. They are populated per-rotation
      ! by cropgrowth's legacy crop sub-readers (Phase 4g territory);
      ! the strangler adapter does not touch them.

      ! ---------------------------------------------------------------
      ! [GR-SEED 2026-05-25 Task 3] Solute seeding moved to state%solute%init
      ! (called from swap_mod). Scalars, per-layer arrays (ldis/kf/decpot/fdepth)
      ! and cseeptab flatten are all performed there.
      ! ---------------------------------------------------------------

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
      state%surfacewater%osswlm = config%surface_water%osswlm
      state%surfacewater%nmper  = config%surface_water%nmper
      state%surfacewater%swqhr  = config%surface_water%swqhr

      if (allocated(config%surface_water%impend)) then
         if (.not. allocated(state%surfacewater%impend)) then
            allocate(state%surfacewater%impend(mamp))
            state%surfacewater%impend = 0.0d0
         end if
         do i = 1, min(size(config%surface_water%impend), size(state%surfacewater%impend))
            state%surfacewater%impend(i) = config%surface_water%impend(i)
         end do
      end if
      if (allocated(config%surface_water%swman)) then
         if (.not. allocated(state%surfacewater%swman)) then
            allocate(state%surfacewater%swman(mamp))
            state%surfacewater%swman = 0
         end if
         do i = 1, min(size(config%surface_water%swman), size(state%surfacewater%swman))
            state%surfacewater%swman(i) = config%surface_water%swman(i)
         end do
      end if
      if (allocated(config%surface_water%wscap)) then
         if (.not. allocated(state%surfacewater%wscap)) then
            allocate(state%surfacewater%wscap(mamp))
            state%surfacewater%wscap = 0.0d0
         end if
         do i = 1, min(size(config%surface_water%wscap), size(state%surfacewater%wscap))
            state%surfacewater%wscap(i) = config%surface_water%wscap(i)
         end do
      end if
      if (allocated(config%surface_water%wldip)) then
         if (.not. allocated(state%surfacewater%wldip)) then
            allocate(state%surfacewater%wldip(mamp))
            state%surfacewater%wldip = 0.0d0
         end if
         do i = 1, min(size(config%surface_water%wldip), size(state%surfacewater%wldip))
            state%surfacewater%wldip(i) = config%surface_water%wldip(i)
         end do
      end if
      if (allocated(config%surface_water%intwl)) then
         if (.not. allocated(state%surfacewater%intwl)) then
            allocate(state%surfacewater%intwl(mamp))
            state%surfacewater%intwl = 0
         end if
         do i = 1, min(size(config%surface_water%intwl), size(state%surfacewater%intwl))
            state%surfacewater%intwl(i) = config%surface_water%intwl(i)
         end do
      end if
      ! Note: alphaw arrays carry the post-finalize-normalized values
      ! per Phase 4f-prep Task D2 — the adapter copies them as-is.
      if (allocated(config%surface_water%hbweir)) then
         if (.not. allocated(state%surfacewater%hbweir)) then
            allocate(state%surfacewater%hbweir(mamp))
            state%surfacewater%hbweir = 0.0d0
         end if
         do i = 1, min(size(config%surface_water%hbweir), size(state%surfacewater%hbweir))
            state%surfacewater%hbweir(i) = config%surface_water%hbweir(i)
         end do
      end if
      if (allocated(config%surface_water%alphaw)) then
         if (.not. allocated(state%surfacewater%alphaw)) then
            allocate(state%surfacewater%alphaw(mamp))
            state%surfacewater%alphaw = 0.0d0
         end if
         do i = 1, min(size(config%surface_water%alphaw), size(state%surfacewater%alphaw))
            state%surfacewater%alphaw(i) = config%surface_water%alphaw(i)
         end do
      end if
      if (allocated(config%surface_water%betaw)) then
         if (.not. allocated(state%surfacewater%betaw)) then
            allocate(state%surfacewater%betaw(mamp))
            state%surfacewater%betaw = 0.0d0
         end if
         do i = 1, min(size(config%surface_water%betaw), size(state%surfacewater%betaw))
            state%surfacewater%betaw(i) = config%surface_water%betaw(i)
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
         state%crop%common%flCropReadFile = .true.   ! [GR-CROP 2026-05-25] flCropReadFile retired
         state%crop%common%flCropOpenFile = .true.   ! [GR-CROP 2026-05-25] flCropOpenFile → state%crop%common
      end if

      ! [GR-IO 2026-05-25 Phase 6 Step 3] rdmax legacy mirror dropped

      if (allocated(config%crop%rotation_type)) then
         n = size(config%crop%rotation_type)
         ! [GR-ATM 2026-05-23] allocate + populate state%crop%common%croptype
         allocate(state%crop%common%croptype(n))
         state%crop%common%croptype = config%crop%rotation_type
      end if
      ! [GR-CROP 2026-05-25] cropstart/cropend legacy global arrays retired —
      ! readers now consume config%crop%rotation_start/rotation_end directly via
      ! state%cfg%crop. No mirror copy needed.
      ! [GR-CROP 2026-05-25] cropfil legacy array retired — readers consume
      ! config%crop%rotation_file (state%cfg%crop%rotation_file) directly.

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
      ! [GR-SEED 2026-05-25] swheader = 0 moved into state%timecontrol%init.
      ! [GR-SOIL 2026-05-24] swcaprise legacy mirror dropped — now state%cfg%simulation%numerical%swcaprise.
      ! [SS-GR-CROPRT A3] swrum adapter write dropped — global retired (always 0; outrume calls dropped)
      ! [GR-FINAL C4] dropped W-globals (all zero by init.f90 or Fortran default):
      !   swafo, swaun, swvap, swbal, swwba, swsba, swblc, swdrf, swstr, swirg,
      !   swini, swcapriseoutput, swswb, swoutputmodflow

      ! [GR-IO 2026-05-25 Phase 6 Step 3] outfil legacy mirror dropped — readers use config%general%outfil

      ! [GR-IO 2026-05-25 Phase 2] CSV output — read directly from config%output_csv
      ! by swap_csv_output.f90; InList_csv/InList_csv_tz legacy mirror writes
      ! dropped along with their bare-global declarations.
      ! Defaults (enabled=1, enabled_tz=0, inlist=water-balance, inlist_tz=wc,h,conc,
      ! tz_z1_z2=[0,0]) are applied by output_csv_config_finalize prior to this adapter.

      ! [GR-ATM 2026-05-23] logf bare-global retired; swap_log opens
      ! 'swap_swap.log' via log_init() in swap_main.

   end subroutine config_to_variables

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
      ! [GR-CROP 2026-05-25] writes redirected from variables.f90 till_* legacy
      ! globals to state%tillage Group AB fields. The till_* declarations are
      ! retired in this commit; this subroutine is the only writer.
      use, intrinsic :: iso_fortran_env, only: real64
      use soil_config_mod, only: soil_tillage_t
      use error_mod, only: fatalerr_collected
      use swap_state_mod, only: swap_state_t
      type(soil_tillage_t), intent(in)    :: tillage
      real(real64),         intent(in)    :: tend
      type(swap_state_t),   intent(inout) :: state   ! [GR-BH Task 35] replaces NumNod/zbotcp globals

      integer :: i, j

      associate(tl => state%tillage)

      tl%i_n_model = tillage%i_n_model
      tl%iRedist   = tillage%iRedist

      tl%Ntill  = size(tillage%events)
      tl%Ntypes = size(tillage%types)

      ! Per-event arrays (sentinel: Date_tillage(Ntill+1) = tend + 1).
      if (allocated(tl%Date_tillage)) deallocate(tl%Date_tillage); allocate(tl%Date_tillage(tl%Ntill+1))
      if (allocated(tl%Z_tillage))    deallocate(tl%Z_tillage);    allocate(tl%Z_tillage(tl%Ntill))
      if (allocated(tl%I_tillage))    deallocate(tl%I_tillage);    allocate(tl%I_tillage(tl%Ntill))
      if (allocated(tl%Type_Tillage)) deallocate(tl%Type_Tillage); allocate(tl%Type_Tillage(tl%Ntill))

      do i = 1, tl%Ntill
         tl%Z_tillage(i)    = tillage%events(i)%z
         tl%I_tillage(i)    = tillage%events(i)%intensity
         tl%Type_Tillage(i) = tillage%events(i)%type_id
         tl%Date_tillage(i) = parse_iso_date_to_days1900(tillage%events(i)%date)
      end do
      tl%Date_tillage(tl%Ntill + 1) = tend + 1.0_real64

      ! Deferred z-range check (validator can't see numnod / zbotcp).
      ! [GR-BH Task 35] NumNod/zbotcp globals replaced by state%mesh fields.
      ! Guard: mesh is not yet populated when called from config_to_variables
      ! (CalcGrid runs after); skip z-range check until mesh is built.
      if (state%mesh%numnod > 0 .and. allocated(state%mesh%zbotcp)) then
         do i = 1, tl%Ntill
            if (tl%Z_tillage(i) < 0.0_real64 .or. &
                tl%Z_tillage(i) > abs(state%mesh%zbotcp(state%mesh%numnod))) then
               call fatalerr_collected('apply_soil_tillage', &
                  'Z_tillage value is outside the model grid depth range')
            end if
         end do
      end if

      ! Per-type arrays.
      if (allocated(tl%iType_Tillage))   deallocate(tl%iType_Tillage);   allocate(tl%iType_Tillage(tl%Ntypes))
      if (allocated(tl%TAB_Rho_cons))    deallocate(tl%TAB_Rho_cons);    allocate(tl%TAB_Rho_cons(tl%Ntypes))
      if (allocated(tl%TAB_Rho_tillage)) deallocate(tl%TAB_Rho_tillage); allocate(tl%TAB_Rho_tillage(tl%Ntypes))
      if (allocated(tl%TAB_K_R_cons))    deallocate(tl%TAB_K_R_cons);    allocate(tl%TAB_K_R_cons(tl%Ntypes))

      do i = 1, tl%Ntypes
         tl%iType_Tillage(i)   = tillage%types(i)%id
         tl%TAB_Rho_cons(i)    = tillage%types(i)%rho_cons
         tl%TAB_Rho_tillage(i) = tillage%types(i)%rho_tillage
         tl%TAB_K_R_cons(i)    = tillage%types(i)%k_R
      end do

      if (tl%i_n_model == 3) then
         if (allocated(tl%TAB_Rho_match)) deallocate(tl%TAB_Rho_match); allocate(tl%TAB_Rho_match(tl%Ntypes))
         if (allocated(tl%TAB_N_match))   deallocate(tl%TAB_N_match);   allocate(tl%TAB_N_match(tl%Ntypes))
         do i = 1, tl%Ntypes
            tl%TAB_Rho_match(i) = tillage%types(i)%rho_match
            tl%TAB_N_match(i)   = tillage%types(i)%N_match
         end do
      end if

      tl%Max_Z_tillage = maxval(tl%Z_tillage(1:tl%Ntill))

      ! iTT1 / iTT2: first/last position per tillage type in iType_Tillage.
      ! Replicates the loop from the legacy Read_Tillage subroutine verbatim.
      if (allocated(tl%iTT1)) deallocate(tl%iTT1); allocate(tl%iTT1(tl%Ntill)); tl%iTT1 = 0
      if (allocated(tl%iTT2)) deallocate(tl%iTT2); allocate(tl%iTT2(tl%Ntill)); tl%iTT2 = 0
      do j = 1, tl%Ntill
         do i = 1, tl%Ntypes
            if (tl%iTT1(j) == 0 .and. tl%iType_Tillage(i) == j) tl%iTT1(j) = i
            if (tl%iTT1(j) >  0 .and. tl%iType_Tillage(i) == j) tl%iTT2(j) = i
         end do
      end do

      end associate
   end subroutine apply_soil_tillage


   ! [GR-SEED 2026-05-25 Task 4] apply_nutrients + apply_nutrients_events bodies
   ! relocated to nutrients_state_mod as seed_nutrients_from_config +
   ! load_nutrients_events. Called from state%nutrients%init in swap_mod.

   ! [GR-SEED 2026-05-25 Task 5] apply_irrigation_ssdi + apply_ssdi_mode0 +
   ! apply_ssdi_mode1 bodies relocated to crop_irrigation_state_mod as
   ! apply_ssdi_seed (public) + apply_ssdi_mode0/mode1 (private). Called from
   ! state%crop%irrigation%init in swap_mod.

end module config_to_variables_mod
