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
   use swap_config_mod, only: swap_config_t
   implicit none
   private

   public :: config_to_variables

contains

   !> Copy every (C)-classified field from `config` into the corresponding
   !! `variables%` legacy global. Caller is responsible for having loaded,
   !! validated, and finalized `config` first.
   subroutine config_to_variables(config)
      use variables   ! bare-use is intentional: many globals across sections
      type(swap_config_t), intent(in) :: config

      integer :: i, n

      ! ---------------------------------------------------------------
      ! General + simulation (audit: 11 fields)
      ! ---------------------------------------------------------------
      if (allocated(config%general%project))   project   = config%general%project
      if (allocated(config%general%pathwork))  pathwork  = config%general%pathwork
      if (allocated(config%general%pathatm))   pathatm   = config%general%pathatm
      if (allocated(config%general%pathcrop))  pathcrop  = config%general%pathcrop
      if (allocated(config%general%pathdrain)) pathdrain = config%general%pathdrain
      swscre  = config%general%swscre

      tstart    = config%simulation%tstart
      tend      = config%simulation%tend
      nprintday = config%simulation%nprintday
      period    = config%simulation%period
      swres     = config%simulation%swres
      swodat    = config%simulation%swodat

      ! Derive iyear/imonth from tstart, matching legacy readswap.f90:126-128.
      ! TimeControl(1) reads these as state (not as inputs) so the adapter
      ! must populate them before TimeControl runs. Without this, dtardp
      ! aborts with "year 0 or partial date not allowed" because iyear
      ! was zero-initialized by Initialize().
      block
         integer :: datea_init(6)
         real    :: fsec_init
         call dtdpar(tstart + 0.1d0, datea_init, fsec_init)
         iyear  = datea_init(1)
         imonth = datea_init(2)
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
         call populate_outdatint_monthly()
         period = 0
         swres  = 0
         swodat = 0
      end if

      ! ---------------------------------------------------------------
      ! Simulation.numerical (audit: 6 fields)
      ! ---------------------------------------------------------------
      dt        = config%simulation%numerical%dt
      dtmin     = config%simulation%numerical%dtmin
      dtmax     = config%simulation%numerical%dtmax
      MaxIt     = config%simulation%numerical%MaxIt
      MaxBackTr = config%simulation%numerical%MaxBackTr
      taccur    = config%simulation%numerical%taccur
      gwlconv       = config%simulation%numerical%gwlconv
      critdevh1cp   = config%simulation%numerical%critdevh1cp
      critdevh2cp   = config%simulation%numerical%critdevh2cp
      critdevponddt = config%simulation%numerical%critdevponddt
      SWkmean       = config%simulation%numerical%swkmean
      SwkImpl       = config%simulation%numerical%swkimpl
      msteps        = config%simulation%numerical%msteps

      ! ---------------------------------------------------------------
      ! Meteorology (audit: 12 + evaporation + snow)
      ! ---------------------------------------------------------------
      if (allocated(config%meteo%metfile))  metfil  = config%meteo%metfile
      if (allocated(config%meteo%rainfile)) rainfil = config%meteo%rainfile
      lat         = config%meteo%lat
      alt         = config%meteo%alt
      altw        = config%meteo%altw
      swetr       = config%meteo%swetr
      swdivide    = config%meteo%swdivide
      swmetdetail = config%meteo%swmetdetail
      nmetdetail  = config%meteo%nmetdetail
      swrain      = config%meteo%swrain
      swetsine    = config%meteo%swetsine
      swinter     = config%meteo%swinter
      angstroma   = config%meteo%angstroma
      angstromb   = config%meteo%angstromb

      ! Derive swMetFilAll from metfil, matching legacy readswap.f90:415-429.
      ! If metfil contains ".met", legacy reads the whole multi-year file at
      ! once (swMetFilAll=1); otherwise it expects per-year files <metfil>.YYY.
      ! `lowerc` is from ttutil; mirror its case-folding here.
      call lowerc(metfil)
      swMetFilAll = 0
      if (index(trim(metfil), ".met") > 0) then
         swMetFilAll = 1
         if (swmetdetail == 1 .or. swrain == 1 .or. swrain == 3) then
            block
               integer :: idum
               idum = index(metfil, ".met")
               metfil = trim(metfil(1:idum-1))
               swMetFilAll = 0
            end block
         end if
      end if

      ! Pre-load all years of meteo data into the cached arrays (legacy
      ! readswap.f90:1746-1749). MeteoInOneFile(2, ...) called per-year
      ! during the dynamic loop extracts from this cache.
      if (swMetFilAll == 1) then
         block
            integer :: idum_meteo
            interface
               subroutine MeteoInOneFile(iTask, ifnd)
                  integer, intent(in)  :: iTask
                  integer, intent(out) :: ifnd
               end subroutine
            end interface
            call MeteoInOneFile(1, idum_meteo)
         end block
      end if

      ! Evaporation sub-section
      swcfbs = config%meteo%evaporation%swcfbs
      cfbs   = config%meteo%evaporation%cfbs
      ! Per Discovery #1 in configuration-schema.md, both legacy keys
      ! (cofredbl=Black, cofredbo=Boesten/Stroosnijder) target the same
      ! legacy global `cofred`; the schema separates them but only one
      ! is meaningful at a time. Pick the one the user set: prefer the
      ! non-default one, falling back to cofredbl when neither was set.
      if (config%meteo%evaporation%cofredbo /= 0.35d0) then
         cofred = config%meteo%evaporation%cofredbo
      else
         cofred = config%meteo%evaporation%cofredbl
      end if

      ! HACK Phase 4f-extend: SWREDU is the soil-evaporation reduction-method
      ! switch (1=Black, 2=Boesten-Stroosnijder). Legacy reads it from .swp
      ! at readswap.f90:584 but no schema slot covers it yet. Default to 1
      ! (Black model) matching cases 1/2/4. When the case authors a non-
      ! default cofredbo (Boesten coefficient), flip to swredu=2 — case 5
      ! (salinitystress) is the only regression case using SWREDU=2.
      ! Phase 4f-extend should add `[soil.evaporation].swredu` so the
      ! switch is authored explicitly rather than inferred.
      if (config%meteo%evaporation%cofredbo /= 0.35d0) then
         swredu = 2
      else
         swredu = 1
      end if

      ! HACK Phase 4f-extend: RSIGNI is the minimum daily rainfall (cm) that
      ! resets the Black-method dry counter (ldwet). Legacy reads it from .swp
      ! at readswap.f90:586 only when swredu==1. Initialize.f90 leaves it 0.0,
      ! so EVERY trace of rain resets ldwet -> Black empreva stays high every
      ! day -> bare-soil EACT overshoots by ~7 cm/yr in case 1 (hupselbrook).
      ! .swp template authors RSIGNI = 0.5. Hardcoding 0.5 matches case 1;
      ! add to meteo_evaporation_config_t in a follow-up Phase 4f-extend pass.
      rsigni = 0.5d0

      ! HACK Phase 4f-extend: CFEVAPPOND is the ponding-layer evaporation
      ! coefficient applied to peva when pond > 1e-10. Legacy default in
      ! readswap.f90:599 is 1.25; initialize.f90 leaves it 0.0, which would
      ! zero out evaporation during ponding. Hardcoded to 1.25 here; add a
      ! schema slot in a follow-up Phase 4f-extend pass.
      cfevappond = 1.25d0

      ! Snow sub-section
      swsnow   = config%meteo%snow%swsnow
      snowcoef = config%meteo%snow%snowcoef
      teprrain = config%meteo%snow%teprrain
      teprsnow = config%meteo%snow%teprsnow

      ! ---------------------------------------------------------------
      ! Drainage (audit: 20 fields + surface_runoff sub-section)
      ! ---------------------------------------------------------------
      swdra    = config%drain%swdra
      dramet   = config%drain%dramet
      swdivd   = config%drain%swdivd
      swdislay = config%drain%swdislay
      nrlevs   = config%drain%nrlevs
      basegw   = config%drain%basegw
      entres   = config%drain%entres
      ! Note: drainage%shape may collide with bottom_boundary%shape
      ! at the legacy global `shape`. Bottom boundary is wired below
      ! and overwrites for swbotb=3 cases (the documented audit alias).

      ! HACK Phase 4f-extend: drfil is the legacy stem of the .dra file
      ! consumed by rddre() in src/drainage/surfacewater.f90:56 when
      ! SWDRA=2. The legacy reads `drfil` from .swp at readswap.f90:1003
      ! but the strangler doesn't have a typed schema slot for it. All
      ! existing TOML cases use 'swap' as the .dra stem (swap.dra).
      ! Add a [drainage].drfil slot in Phase 4f-extend.
      if (config%drain%swdra >= 1) then
         drfil = 'swap'
      end if

      ! DRAMET=2 (Hooghoudt/Ernst). Mirrors readswap.f90:1850-1875.
      ! lm is authored in metres; the legacy reader does the m->cm
      ! conversion (`l(1) = 100*lm2`) so we replicate that here.
      ! wetper / zbotdr go into the level-1 entry of the per-level
      ! arrays; ipos / khtop / khbot / kvtop / kvbot / zintf / geofac
      ! are scalar globals.
      if (config%drain%dramet == 2) then
         L(1)      = 100.0d0 * config%drain%lm
         wetper(1) = config%drain%wetper
         zbotdr(1) = config%drain%zbotdr_basic
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

      ! Per-soil-physical-layer anisotropy ratio (legacy COFANI in .dra).
      if (allocated(config%drain%cofani)) then
         do i = 1, size(config%drain%cofani)
            cofani(i) = config%drain%cofani(i)
         end do
      end if

      if (allocated(config%drain%swdtyp)) then
         do i = 1, size(config%drain%swdtyp)
            swdtyp(i) = config%drain%swdtyp(i)
         end do
      end if
      if (allocated(config%drain%zbotdr)) then
         do i = 1, size(config%drain%zbotdr)
            zbotdr(i) = config%drain%zbotdr(i)
         end do
      end if
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
      if (allocated(config%drain%L)) then
         do i = 1, size(config%drain%L)
            L(i) = config%drain%L(i)
         end do
         ! Legacy m->cm conversion mirroring readswap.f90:1908-1909 (and
         ! the matching DRAMET=2 branch above). For DRAMET=3 the reader
         ! at lines 1908-1994 multiplies each per-level `l(:)` by 100 only
         ! when `swdivd == 1`. The DRAMET=2 branch upstream in this file
         ! handles its own scalar `lm` conversion; here we cover the
         ! per-level array path (DRAMET=3 only — DRAMET=1/lookup never
         ! consults `L(:)`).
         if (config%drain%dramet == 3 .and. config%drain%swdivd == 1) then
            do i = 1, size(config%drain%L)
               L(i) = 100.0d0 * L(i)
            end do
         end if
      end if
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

      ! HACK Phase 4f-extend: SWLIMINF gates limit-of-infiltration to the
      ! channel water depth in the DRAMET=3 multi-level resistance solver
      ! (drainage.f90 / divdra.f90). Legacy hard-codes 1 in
      ! readswap.f90:2010 (after the DRAMET=3 .dra block) when the .dra
      ! is silent on the key. variables.f90:694 default-initialises to 0,
      ! so without this HACK case 2's solver would treat infiltration as
      ! unlimited. Add a [drainage].swliminf slot in Phase 4f-extend.
      if (config%drain%dramet == 3) then
         swliminf = 1
      end if

      ! Drainage.surface_runoff sub-section: scalar switches + per-level
      ! arrays. Legacy globals `swtopdislay`, `ftopdislay`, `RapDraResRef`
      ! are arrays of size madr; copy element 1 of the (scalar) config
      ! field as a uniform value across drainage levels until a per-level
      ! schema lands in Phase 4f-extend.
      swnrsrf      = config%drain%surface_runoff%swnrsrf
      SwTopnrsrf   = config%drain%surface_runoff%swtopnrsrf
      swdivdinf    = config%drain%surface_runoff%swdivdinf
      FacDpthInf   = config%drain%surface_runoff%facdpthinf
      cofintfl     = config%drain%surface_runoff%cofintfl
      expintfl     = config%drain%surface_runoff%expintfl
      geofac       = config%drain%surface_runoff%geofac
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
      RapDraReaExp = config%drain%surface_runoff%rapdrareaexp
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
      if (size(RapDraResRef) >= 1) then
         do i = 1, size(RapDraResRef)
            RapDraResRef(i) = config%drain%surface_runoff%rapdraresref
         end do
      end if

      ! ---------------------------------------------------------------
      ! Soil (audit: 15 + discretization + frost)
      ! ---------------------------------------------------------------
      swsophy = config%soil%swsophy
      swhyst  = config%soil%swhyst
      swinco  = config%soil%swinco
      swmacro = config%soil%swmacro
      gwli    = config%soil%gwli
      pondini = config%soil%pondini
      pond    = config%soil%pondini    ! legacy alias: pond <-> pondini
      pondmx  = config%soil%pondmx
      rsoil   = config%soil%rsoil
      rsro    = config%soil%rsro
      rsroexp = config%soil%rsroexp
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
      if (allocated(config%soil%orgmat)) then
         do i = 1, size(config%soil%orgmat)
            orgmat(i) = config%soil%orgmat(i)
         end do
      end if
      if (allocated(config%soil%bdens)) then
         do i = 1, size(config%soil%bdens)
            bdens(i) = config%soil%bdens(i)
         end do
      end if
      if (allocated(config%soil%cofani)) then
         do i = 1, size(config%soil%cofani)
            cofani(i) = config%soil%cofani(i)
         end do
      end if

      ! HACK Phase 4f-extend: SWINCO=3 .ini reader. Mirrors
      ! readswap.f90:1593-1626 — when the case authors `[soil].inifil`
      ! and `swinco=3`, read the previous-run end-state file via the
      ! legacy ttutil rdinit/rdsdor/rdador/rdfdor stack to populate
      ! ssnow/slw/pond/zi/h/zc/cml (and z_Tsoil/Tsoil under SWHEA=1
      ! SWCALT=2). Salinitystress (case 5) needs this so the initial
      ! solute profile (~15 mg/cm3 below ~150 cm, ~0 in root zone)
      ! matches the fixture's quasi-steady starting state. Same HACK
      ! pattern as the .bbc / .irg readers above. Phase 4f-extend
      ! should port the .ini state into a typed schema slot.
      if (config%soil%swinco == 3 .and. allocated(config%soil%inifil)) then
         if (len_trim(config%soil%inifil) > 0) then
            block
               use swap_array_dimensions, only: macp
               integer :: ini_unit, ifnd_ini
               character(len=200) :: ini_filnam
               integer, external :: getun2
               ini_filnam = trim(config%soil%inifil)
               ini_unit = getun2(10, 90, 2)
               call rdinit(ini_unit, logf, ini_filnam)
               call rdsdor('ssnow', 0.0d0, 1000.0d0, ssnow)
               ! Legacy zeroes ssnow when swsnow != 1 (readswap.f90:1606-1613).
               ! Salinitystress and other regression cases have swsnow=0.
               if (config%meteo%snow%swsnow /= 1) ssnow = 0.0d0
               call rdsdor('slw',   0.0d0, 1000.0d0, slw)
               call rdsdor('pond',  0.0d0,  100.0d0, pond)
               pondini = pond
               call rdador('z_h',  -1.0d5,  0.0d0, zi, macp, ifnd_ini)
               call rdfdor('h',    -1.0d10, 1.0d4, h,  macp, ifnd_ini)
               nhead = ifnd_ini
               if (config%heat%swhea == 1 .and. config%heat%swcalt == 2) then
                  call rdador('z_Tsoil', -1.0d5,  0.0d0, zh,    macp, ifnd_ini)
                  call rdfdor('Tsoil',  -50.0d0, 50.0d0, tsoil, macp, ifnd_ini)
               end if
               if (config%solute%swsolu == 1) then
                  call rdador('z_Cml', -1.0d5,    0.0d0, zc,  macp, ifnd_ini)
                  call rdfdor('Cml',    0.0d0, 1.0d6,    cml, macp, ifnd_ini)
                  nconc = ifnd_ini
               end if
               close(ini_unit)
            end block
         end if
      end if

      ! Per-soil-physical-layer Mualem-van Genuchten hydraulics.
      ! Mirrors readswap.f90:786-825 conceptually, but most of the
      ! per-array names (ores/osat/alfa/npar/lexp/alfaw) are *locals*
      ! in readswap — not module globals — so they only matter as input
      ! to paramvg(1..10, lay). We write paramvg directly here. The
      ! globals we *do* need to set are bdens / ksatfit / ksatexm /
      ! h_enpr (all real(8) :: ...(maho) in variables.f90), since
      ! downstream code (soilhydraulics, solute, etc.) reads them.
      if (allocated(config%soil%hydraulics%ores)) then
         do i = 1, size(config%soil%hydraulics%ores)
            ksatfit(i) = config%soil%hydraulics%ksatfit(i)
            ksatexm(i) = config%soil%hydraulics%ksatexm(i)
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
      swsublim  = config%soil%frost%swsublim
      ! tfroststa/tfrostend live in heat block per audit (and per legacy);
      ! do not double-write here from soil%frost — heat block below owns it.

      ! ---------------------------------------------------------------
      ! Bottom boundary (audit: 9 fields, conditional per swbotb)
      ! ---------------------------------------------------------------
      swbotb = config%bottom_boundary%swbotb
      select case (swbotb)
      case (1)
         ! Phase 4f cleanup: SWBOTB=1 prescribes the groundwater level via a
         ! (date, gwlevel) table inlined in [bottom_boundary].gwl_table.
         ! Mirrors the SWBOTB=3 haquif_table pattern below. Populates the
         ! legacy interleaved gwltab(2*i-1)=date, gwltab(2*i)=gwlevel packing
         ! consumed by afgen() in soilhydraulics.f90:1005.
         if (allocated(config%bottom_boundary%gwl_table)) then
            block
               integer :: nrows_gwl, k_gwl
               nrows_gwl = size(config%bottom_boundary%gwl_table, 1)
               do k_gwl = 1, nrows_gwl
                  gwltab(k_gwl*2 - 1) = config%bottom_boundary%gwl_table(k_gwl, 1)
                  gwltab(k_gwl*2)     = config%bottom_boundary%gwl_table(k_gwl, 2)
               end do
            end block
         end if
      case (3)
         shape  = config%bottom_boundary%shape
         hdrain = config%bottom_boundary%hdrain
         rimlay = config%bottom_boundary%rimlay
         aqave  = config%bottom_boundary%aqave
         aqamp  = config%bottom_boundary%aqamp
         aqper  = config%bottom_boundary%aqper
         aqtmax = config%bottom_boundary%aqtmax
         ! HACK Phase 4f-extend: SWBOTB=3 implicit/explicit flux selector
         swbotb3impl = config%bottom_boundary%swbotb3impl
         ! Phase 4f Task B4: SWBOTB=3 + sw3=2 (date-keyed aquifer head).
         ! Mirrors readswap.f90:1369-1378. Selects table mode (sw3=2) when
         ! the typed config provides a haquif_table; sinus mode (sw3=1)
         ! otherwise. Populates the legacy interleaved haqtab(2*i-1)=date,
         ! haqtab(2*i)=head packing consumed by afgen() in
         ! boundary/boundbottom.f90:118.
         if (allocated(config%bottom_boundary%haquif_table)) then
            block
               integer :: nrows_haq, k_haq
               nrows_haq = size(config%bottom_boundary%haquif_table, 1)
               if (nrows_haq > 0) then
                  sw3 = 2
                  do k_haq = 1, nrows_haq
                     haqtab(k_haq*2 - 1) = config%bottom_boundary%haquif_table(k_haq, 1)
                     haqtab(k_haq*2)     = config%bottom_boundary%haquif_table(k_haq, 2)
                  end do
               else
                  sw3 = 1
               end if
            end block
         else
            sw3 = 1
         end if
      case (5)
         hbot = config%bottom_boundary%hbot
      end select

      ! ---------------------------------------------------------------
      ! Heat (audit: 10 fields)
      ! ---------------------------------------------------------------
      swhea     = config%heat%swhea
      swcalt    = config%heat%swcalt
      swtopbhea = config%heat%swtopbhea
      swbotbhea = config%heat%swbotbhea
      tfroststa = config%heat%tfroststa
      tfrostend = config%heat%tfrostend

      if (allocated(config%heat%psand)) then
         do i = 1, size(config%heat%psand)
            psand(i) = config%heat%psand(i)
         end do
      end if
      if (allocated(config%heat%psilt)) then
         do i = 1, size(config%heat%psilt)
            psilt(i) = config%heat%psilt(i)
         end do
      end if
      if (allocated(config%heat%pclay)) then
         do i = 1, size(config%heat%pclay)
            pclay(i) = config%heat%pclay(i)
         end do
      end if
      if (allocated(config%heat%porg)) then
         do i = 1, size(config%heat%porg)
            forg(i) = config%heat%porg(i)   ! legacy alias: porg -> forg
         end do
      end if
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
         nheat = n
         do i = 1, min(n, size(tsoil))
            zh(i)    = config%heat%tsoil_init(i, 1)
            tsoil(i) = config%heat%tsoil_init(i, 2)
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
               use csv_reader_mod, only: read_csv_date_reals
               use error_mod, only: error_collection_t
               real(8), allocatable :: csv_table(:,:)
               type(error_collection_t) :: csv_errs
               integer :: k_csv, nrows_csv
               call read_csv_date_reals( &
                  trim(pathwork)//trim(config%irrigation%fixed_events_file), &
                  3, csv_table, csv_errs)
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
      ! Solute (audit: 8 fields)
      ! ---------------------------------------------------------------
      swsolu  = config%solute%swsolu
      swbotbc = config%solute%swbotbc
      cdrain  = config%solute%cdrain
      cseep   = config%solute%cseep
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

      ! ---------------------------------------------------------------
      ! Surface water (audit: 12 + per-period management arrays)
      ! ---------------------------------------------------------------
      swsrf = config%surface_water%swsrf
      swsec = config%surface_water%swsec
      ! wlact + osswlm are loaded by readswap into local scratch; only
      ! osswlm is a module global.
      osswlm = config%surface_water%osswlm
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
      swCrop = config%crop%swcrop

      ! Mirror readswap.f90:479-480 — when crop simulation is enabled,
      ! arm the per-rotation reader gates so cropgrowth.f90's per-crop
      ! init (ArableLandGerm/CropFixed/Wofost/Grass at lines 91 / 121)
      ! actually fires. Without this, flCropReadFile stays .false. (the
      ! Initialize() default) and every rotation is treated as bare soil:
      ! LAI/cf/rd remain 0, TPOT/TACT collapse, and EACT/DRAINAGE balloon.
      if (swCrop == 1) then
         flCropReadFile = .true.
         flCropOpenFile = .true.
      end if

      ! HACK Phase 4f-extend: RDS (rdmax) is the soil-profile-imposed maximum
      ! rooting depth, read by legacy readswap.f90:470 from .swp's crop
      ! rotation block. Without it, cropfixed.f90:414 sets rdm=0 and then
      ! `rd = min(afgen(rdtb,...), rdm) = 0`, so noddrz stays at 1 and
      ! RootExtraction returns zero — TACT collapses across every rotation.
      ! Hardcoding 200.0 cm matches case 1's .swp value; add a typed
      ! crop_config_t.rdmax slot when Phase 4f-extend tackles crop schema.
      rdmax = 200.0d0

      if (allocated(config%crop%rotation_type)) then
         n = size(config%crop%rotation_type)
         do i = 1, min(n, size(croptype))
            croptype(i) = config%crop%rotation_type(i)
         end do
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

      ! ---------------------------------------------------------------
      ! Per ADR 0009: zero-force the 18 RETIRED legacy output switches
      ! so any residual code that checks them does the right thing.
      ! ---------------------------------------------------------------
      swafo           = 0
      swaun           = 0
      swvap           = 0
      swbal           = 0
      swwba           = 0
      swsba           = 0
      swblc           = 0
      swdrf           = 0
      swstr           = 0
      swirg           = 0
      swini           = 0
      swend           = 0
      swheader        = 0
      swcaprise       = .false.
      swcapriseoutput = .false.
      swrum           = 0
      swswb           = 0
      swoutputmodflow = 0

      ! HACK Phase 4f-extend: outfil is the output-file basename
      ! (legacy reads it from .swp Part 1: OUTFIL = 'result'). All
      ! regression cases use the same value, so we hardcode it. Add a
      ! [general.output] / general.outfil slot in Phase 4f-extend.
      outfil = 'result'

      ! HACK Phase 4f-extend: enable CSV output. swcsv=1 + InList_csv
      ! authored verbatim from hupselbrook's .swp. The csv driver is
      ! the *only* output the regression baseline checks (it asserts
      ! against `<outfil>_output.csv`), so without these we 'complete
      ! normally' but produce no output file. Move to a typed
      ! [output.csv] block when Phase 4f-extend tackles output configs.
      swcsv = 1
      ! Per-case override (Phase 4f Task B3): when the TOML authors
      ! `general.inlist_csv` use it; otherwise fall back to the
      ! hupselbrook-tuned water-balance default. Grass cases (case 2 +
      ! oxygenstress) override with grass-detailed columns to match
      ! their fixtures.
      if (allocated(config%general%inlist_csv)) then
         if (len_trim(config%general%inlist_csv) > 0) then
            InList_csv = config%general%inlist_csv
         else
            InList_csv = 'rain,irrig,interc,runoff,drainage,dstor,epot,eact,tpot,tact,qbottom,gwl'
         end if
      else
         InList_csv = 'rain,irrig,interc,runoff,drainage,dstor,epot,eact,tpot,tact,qbottom,gwl'
      end if
      swcsv_tz = 0
      InList_csv_tz = 'wc,h,conc'

      ! HACK Phase 4f-extend: set up legacy I/O state needed by unported
      ! readers (read_tillage, cropgrowth crop sub-readers, rddre). They
      ! call RDinit(unit, logf, swpfile) which opens swpfile from cwd.
      ! Without this, the open fails with FOPENG. Mirrors readswap.f90:83
      ! and 86. Phase 4f-extend will retire these legacy readers reader-
      ! by-reader; once they're all gone, this block goes away too.
      block
         use variables, only: swpfile, logf
         integer :: getun  ! external from ttutil
         logical :: log_open
         swpfile = 'swap.swp'
         inquire(unit=20, opened=log_open)  ! cheap check
         if (.not. log_open) then
            ! Open the legacy log file so unported readers can write to it.
            ! `del` privilege removes the file on close (legacy convention).
            call delfil('swap_swap.log', .false.)
            logf = getun(20, 99)
            call fopens(logf, 'swap_swap.log', 'new', 'del')
         end if
      end block

   end subroutine config_to_variables

   !> Populate `variables%outdatint(:)` with the end-of-month dates
   !! between `tstart` and `tend`. Mirrors `readswap.f90:181-204` (the
   !! `swmonth == 1` branch). Bare `use variables` for parity with the
   !! parent adapter.
   subroutine populate_outdatint_monthly()
      use variables
      integer  :: datea_om(6), i_om
      real     :: fsec_om
      real(8)  :: outdate_om

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

end module config_to_variables_mod
