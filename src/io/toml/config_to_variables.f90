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

      ! ---------------------------------------------------------------
      ! Simulation.numerical (audit: 6 fields)
      ! ---------------------------------------------------------------
      dt        = config%simulation%numerical%dt
      dtmin     = config%simulation%numerical%dtmin
      dtmax     = config%simulation%numerical%dtmax
      MaxIt     = config%simulation%numerical%MaxIt
      MaxBackTr = config%simulation%numerical%MaxBackTr
      taccur    = config%simulation%numerical%taccur

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
      swMetFilAll = config%meteo%swmetfilall
      angstroma   = config%meteo%angstroma
      angstromb   = config%meteo%angstromb

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
      gwlconv      = config%drain%surface_runoff%gwlconv
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
      nrstaring = config%soil%nrstaring

      if (allocated(config%soil%isoillay)) then
         do i = 1, size(config%soil%isoillay)
            isoillay(i) = config%soil%isoillay(i)
         end do
      end if
      if (allocated(config%soil%ncomp)) then
         do i = 1, size(config%soil%ncomp)
            ncomp(i) = config%soil%ncomp(i)
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
      case (3)
         shape  = config%bottom_boundary%shape
         hdrain = config%bottom_boundary%hdrain
         rimlay = config%bottom_boundary%rimlay
         aqave  = config%bottom_boundary%aqave
         aqamp  = config%bottom_boundary%aqamp
         aqper  = config%bottom_boundary%aqper
         aqtmax = config%bottom_boundary%aqtmax
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
      ! Initial soil temperature table tsoil_init(:,2) — column 2 (temp)
      ! flattens into legacy `tsoil(numnod)` row-by-row. Legacy convention:
      ! the (depth, temp) pairs feed an interpolation against `z` to fill
      ! tsoil(macp); for parity at the assignment level we copy the table
      ! directly into the lower-indexed entries of tsoil. Real interpolation
      ! happens in physics post-init.
      if (allocated(config%heat%tsoil_init)) then
         n = size(config%heat%tsoil_init, 1)
         do i = 1, min(n, size(tsoil))
            tsoil(i) = config%heat%tsoil_init(i, 2)
         end do
      end if

      ! ---------------------------------------------------------------
      ! Irrigation (audit: 8 fields, top-level only)
      ! ---------------------------------------------------------------
      swirfix = config%irrigation%swirfix
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
      ! ldis(maho): config side has no per-layer schema yet; leave as-is.

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
            cropfil(i) = config%crop%rotation_file(i)
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

   end subroutine config_to_variables

end module config_to_variables_mod
