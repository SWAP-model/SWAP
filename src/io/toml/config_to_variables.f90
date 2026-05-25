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
   ! [GR-SEED 2026-05-25 Task 9] apply_soil_tillage retired: body moved to
   ! tillage_state_mod as tillage_state_init (type-bound). parse_iso_date_to_days1900
   ! also moved there as a private helper (single caller).
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
      use swap_array_dimensions, only: macp
      use swap_state_mod, only: swap_state_t
      use error_mod, only: fatalerr_collected
      type(swap_config_t), intent(inout), target :: config  ! target: state%cfg => config pointer; inout retained for other callee mutations
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
      ! Drainage
      ! [GR-SEED 2026-05-25 Task 7] Drainage seeding moved to state%drainage%init
      ! (called from swap_mod after CalcGrid). owltab CSV pre-load now inside
      ! drainage_state_init (Piece D); adapter block deleted.
      !
      ! EXCEPTION: state%surfacewater%swdra must remain here because
      ! timecontrol_init (line ~170) reads it to set flDrain/flSurfaceWater —
      ! which runs before state%drainage%init and state%surfacewater%init.
      ! ---------------------------------------------------------------
      state%surfacewater%swdra = config%drain%swdra

      ! ---------------------------------------------------------------
      ! Soil (audit: 15 + discretization + frost)
      ! ---------------------------------------------------------------
      ! [GR-SEED 2026-05-25 Task 8] soilwater scalar seeding (swsophy/swinco/flrunon/
      ! bdens/swfrost) + swinco=3 h_file CSV (→ config%soil%initial%z_init) +
      ! bottom_boundary CSV pre-loads (swbotb=1..5) moved to state%soilwater%init
      ! (called from swap_mod after CalcGrid). Cross-subsystem writes retained below.
      ! [MACRO-RETIRE 2026-05-12] swmacro global retired (ADR 0040).
      ! [GR-CROP 2026-05-25] till_swtill mirror retired — tillage reads state%cfg%soil%swtill directly.
      ! [GR-IO 2026-05-25 Phase 5] flTillage bare global retired — swap_mod/adapter read
      ! (config%soil%swtill == 1) directly.
      ! [GR-TIME 2026-05-25] flSSDI bare global retired — swap_mod/timecontrol_mod read
      ! (config%irrigation%swssdi == 1) directly.
      ! [GR-SEED 2026-05-25 Task 9] Tillage seeding moved to state%tillage%init
      ! (called from swap_mod after CalcGrid). The swtill==1 gate now lives inside the init.
      ! [GR-SEED 2026-05-25 Task 5] Irrigation seeding moved to state%crop%irrigation%init.
      ! [GR-SEED 2026-05-25 Task 4] apply_nutrients moved to state%nutrients%init.
      ! [GR-SOIL 2026-05-24] gwli legacy mirror dropped — direct config read.
      ! [GR-FINAL C1] pondini/pond: read directly by swap_mod after soilwater_init.
      state%surfacewater%pondmx = config%soil%pondmx
      ! [GR-ATM 2026-05-23] rsoil retired — snapshotted in atmosphere_state%init
      state%surfacewater%rsro = config%soil%rsro
      state%surfacewater%rsroexp = config%soil%rsroexp
      ! [GR-IO 2026-05-25 Phase 6 Step 3] nrstaring legacy mirror dropped

      ! sublay (legacy 'isublay') is a local in readswap, not a module global;
      ! calcgrid only consumes nsublay + isoillay + ncomp + hcomp.
      ! [GR-SOIL 2026-05-24] sublay/isoillay/hsublay/ncomp/hcomp bare-global writes
      ! retired — CalcGrid reads them inline from config%soil%X.
      if (allocated(config%soil%isoillay)) then
         state%mesh%numlay = config%soil%isoillay(size(config%soil%isoillay))
      end if
      ! [GR-SOIL 2026-05-24] config%soil%hcomp ingest retired — CalcGrid reads inline.
      ! [GR-BH Task 36] orgmat/cofani: consumed via config in swap_mod seeding block.

      ! [soil.initial] swinco=3 cross-subsystem writes retained here:
      !   state%atmosphere%atmin7 (atmosphere write — cannot go to soilwater%init).
      !   state%solute%X (solute write — cannot go to soilwater%init).
      ! soilwater%init absorbs the z_init CSV and soilwater scalar seeding.
      if (config%soil%swinco == 3) then
         if (allocated(config%soil%initial%h_file) .and. &
             len_trim(config%soil%initial%h_file) > 0) then
            ! [GR-IO 2026-05-25 Phase 6 Step 3] atmin7 → state%atmosphere directly
            state%atmosphere%atmin7(:) = config%soil%initial%atmin7(:)
            ! [SS-ATM A-2.6] Legacy zeroes ssnow when swsnow != 1: now handled in swap.f90 during state seeding

            ! [GR-IO 2026-05-25 Phase 6 Step 3] Legacy tsoil_file CSV → zh/tsoil
            ! bare-global block dropped. temperature.f90:case(1) reads
            ! cfg_heat%tsoil_init directly (populated by read_heat_toml).

            ! Optional: initial concentration profile (Cml) — solute write, stays here.
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
      ! [GR-SOIL 2026-05-24] h_enpr legacy mirror dropped — vg_params carries the typed value.
      ! [GR-SOIL 2026-05-24] iHWCKmodel legacy write retired — soilwater_state_init seeds it directly.
      ! [GR-CROP 2026-05-25] paramvg legacy mirror retired — tillage.f90 now
      !   mutates state%soilwater%vg_params_layer(:) directly.
      ! [GR-SEED 2026-05-25 Task 8] bdens writes moved to soilwater_state_init (Piece B).

      ! Soil.discretization
      ! [GR-IO 2026-05-25 Phase 6 Step 3] swdiscrvert legacy mirror dropped
      ! [GR-IO 2026-05-25] numnodnew/dznew legacy mirror dropped — swapoutput.f90:checkDiscrVert
      ! reads config%soil%discretization%{numnodnew,dznew} directly via state%cfg.

      ! Soil.frost: swfrost moved to soilwater_state_init (Piece B) [GR-SEED 2026-05-25 Task 8].
      ! [GR-IO 2026-05-25 Phase 6 Step 3] swsublim legacy mirror dropped.
      ! tfroststa/tfrostend live in heat block per audit (and per legacy).

      ! ---------------------------------------------------------------
      ! Bottom boundary (audit: 9 fields, conditional per swbotb)
      ! ---------------------------------------------------------------
      ! [GR-SEED 2026-05-25 Task 8] Bottom-boundary CSV pre-loads (swbotb=1..5
      ! → gwltab/qbotab/haqtab/hbotab) moved to state%soilwater%init via
      ! private seed_bottom_boundary helper (called from swap_mod after CalcGrid).
      ! swbotb: legacy global write retired — state%soilwater%swbotb_runtime sourced
      ! directly from config%bottom_boundary%swbotb in swap_mod.f90 [SS-GR-BH A7].

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
      ! [GR-SEED 2026-05-25 Task 6] Surface-water management-period seeding
      ! moved to state%surfacewater%init (called from swap_mod).

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

   ! [GR-SEED 2026-05-25 Task 9] parse_iso_date_to_days1900 helper relocated to
   ! tillage_state_mod as a private helper (single caller was apply_soil_tillage).
   ! apply_soil_tillage body relocated to tillage_state_mod as tillage_state_init
   ! (type-bound). Both are retired from this adapter.

   ! [GR-SEED 2026-05-25 Task 4] apply_nutrients + apply_nutrients_events bodies
   ! relocated to nutrients_state_mod as seed_nutrients_from_config +
   ! load_nutrients_events. Called from state%nutrients%init in swap_mod.

   ! [GR-SEED 2026-05-25 Task 5] apply_irrigation_ssdi + apply_ssdi_mode0 +
   ! apply_ssdi_mode1 bodies relocated to crop_irrigation_state_mod as
   ! apply_ssdi_seed (public) + apply_ssdi_mode0/mode1 (private). Called from
   ! state%crop%irrigation%init in swap_mod.

end module config_to_variables_mod
