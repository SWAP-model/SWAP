!> Runtime initializer for WOFOST (type-2) crop rotations from typed config.
!!
!! Replaces the side-effects of readwofost(task=1) on the TOML pipeline path.
!! Reads all fields from a populated cropwofost_config_t and writes them to
!! the legacy variables-module globals that the wofost computation loop reads.
!! Also computes cumdens (normalized root density) identical to readwofost's
!! tail math (readswap.f90:3094-3121).
!!
!! Scope (Phase 2, ADR 0015): only the supported subset of switches.
!! Defense-in-depth runtime guards mirror the validator stub-errors.
!!
!! IMPORTANT — locals vs globals:
!!   The following readwofost arguments are LOCAL save variables in wofost()
!!   (cropgrowth.f90), NOT in the variables module. They are returned via
!!   intent(out) arguments so the dispatch block can assign them to its own
!!   local SAVEs:
!!     FraDeceasedLvToSoil  — 3rd arg (intent out), set from cfg%management
!!   The following are NOT set here:
!!     FraHarLosOrm_lv/st/so  (harvest/death losses)
!!       — set by apply_cropwofost_nutrient when flcropnut=true on the
!!         current rotation (ADR 0025 N1). For flcropnut=false rotations
!!         these stay at their wofost() local defaults.
!!     swco2, flco2                                  (CO2 correction — stub=0)
!!     swsoybean, mg, dvsi, dvrmax1, dvrmax2,        (soybean — stub=0)
!!       tmaxdvr, tmindvr, toptdvr, popt, pcrt,
!!       flrfphotoveg, flphenodayl
!!
!! Transitional: removed when the config-passing direction (ADR 0016 Part C)
!! lands and wofost() takes explicit (cfg, state) arguments.
module cropwofost_init_mod
   use iso_fortran_env, only: real64
   use cropwofost_config_mod, only: cropwofost_config_t
   implicit none
   private

   public :: cropwofost_init_from_config
   public :: apply_cropwofost_nutrient

contains

   subroutine cropwofost_init_from_config(cfg, icrop, FraDeceasedLvToSoil, state)
      ! [SS-GR-FINAL B7] DEFERRED: all symbols are config→globals copy targets.
      !   Retirement requires Phase C3 adapter rewrite (config_to_variables.f90 dual-writes).
      use variables, only: &
         ! ET / crop factor — DEFERRED Phase C3
         swcf, cftb, chtb, albedo, rsc, rsw,                                &
         ! Development — DEFERRED Phase C3
         idsl, tsumea, tsumam, dlo, dlc, dtsmtb,                            &
         ! Vernalisation — DEFERRED Phase C3
         verndvs, vernsat, vernbase, vernrtb,                                &
         ! Initial crop state — DEFERRED Phase C3
         tdwi, laiem, rgrlai,                                                &
         ! Green area / assimilation — DEFERRED Phase C3
         slatb, spa, ssa, span, tbase,                                       &
         kdif, kdir, eff, amaxtb, tmpftb, tmnftb,                           &
         ! Biomass conversion — DEFERRED Phase C3
         cvl, cvo, cvr, cvs,                                                 &
         ! Maintenance respiration — DEFERRED Phase C3
         q10, rml, rmo, rmr, rms, rfsetb,                                   &
         ! Partitioning / death rates — DEFERRED Phase C3
         frtb, fltb, fstb, fotb,                                             &
         perdl, rdrrtb, rdrstb,                                              &
         ! Oxygen stress — DEFERRED Phase C3
         swoxygen, swWrtNonox, aeratecrit, hlim1, hlim2u, hlim2l,            &
         ! Drought stress — DEFERRED Phase C3
         swdrought, hlim3h, hlim3l, hlim4, adcrh, adcrl,                    &
         ! Salinity — DEFERRED Phase C3
         swsalinity, saltmax, saltslope, salthead,                           &
         ! Compensation — DEFERRED Phase C3
         swcompensate, swstressor,                                            &
         ! Interception — DEFERRED Phase C3
         swinter, cofab,                                                      &
         ! Root depth — DEFERRED Phase C3
         swrd, rdi, rri, rdc, swdmi2rd, rdctb, rdtb, rlwtb, wrtmax,         &
         swrdc, cumdens,                                                      &
         ! Harvest — DEFERRED Phase C3
         dvsend, swharv,                                                      &
         relmf, swpotrelmf,                                                   &
         ! Irrigation schedule — DEFERRED Phase C3
         schedule,                                                             &
         ! Active crop dynamics (written during init) — DEFERRED Phase C3; dvs retired
         tsum, daycrop, nofd, flCropNut
      use array_utils, only: afgen
      use error_mod,   only: fatalerr_collected
      use swap_state_mod, only: swap_state_t
      implicit none

      type(cropwofost_config_t), intent(in)    :: cfg
      integer,                   intent(in)    :: icrop  ! rotation slot (reserved)
      real(real64),              intent(out)   :: FraDeceasedLvToSoil
      type(swap_state_t),        intent(inout) :: state  ! [SS-GR-ATM A5.1] runtime dual-write target

      integer      :: i
      real(real64) :: depth, sum_dens
      real(real64) :: rootdis(202)

      ! ----------------------------------------------------------------
      ! Defense-in-depth guards — the validator should have caught these,
      ! but we guard at runtime too per ADR 0015.
      ! ----------------------------------------------------------------
      if (cfg%soybean%swsoybean == 1) &
         call fatalerr_collected('cropwofost_init', &
            'swsoybean=1 not supported on TOML path; validator should have rejected.')
      if (cfg%bulb%swbulb == 1) &
         call fatalerr_collected('cropwofost_init', &
            'swbulb=1 not supported on TOML path; validator should have rejected.')
      if (cfg%co2%swco2 == 1) &
         call fatalerr_collected('cropwofost_init', &
            'swco2=1 not supported on TOML path; validator should have rejected.')
      if (cfg%schedule%schedule == 1) &
         call fatalerr_collected('cropwofost_init', &
            'schedule=1 not supported on TOML path; validator should have rejected.')
      if (cfg%drought_stress%swdrought == 2) &
         call fatalerr_collected('cropwofost_init', &
            'swdrought=2 not supported on TOML path; validator should have rejected.')
      if (cfg%oxygen_stress%swoxygen == 2) &
         call fatalerr_collected('cropwofost_init', &
            'swoxygen=2 not supported on TOML path; validator should have rejected.')
      if (cfg%interception%swinter == 2) &
         call fatalerr_collected('cropwofost_init', &
            'swinter=2 not supported on TOML path; validator should have rejected.')
      if (cfg%compensate%swcompensate /= 0) &
         call fatalerr_collected('cropwofost_init', &
            'swcompensate/=0 not supported on TOML path; validator should have rejected.')
      if (cfg%harvest%swharv == 1) &
         call fatalerr_collected('cropwofost_init', &
            'swharv=1 not supported on TOML path; validator should have rejected.')
      if (cfg%salinity%swsalinity == 2) &
         call fatalerr_collected('cropwofost_init', &
            'swsalinity=2 not supported on TOML path; validator should have rejected.')
      if (cfg%root%swrdc == 1) &
         call fatalerr_collected('cropwofost_init', &
            'swrdc=1 not supported on TOML path; validator should have rejected.')

      ! ----------------------------------------------------------------
      ! Config → globals copy (mirroring readwofost's rd* sequence)
      ! ----------------------------------------------------------------

      ! Part 1: crop factor / crop height (readwofost lines 2588-2637)
      swcf = cfg%crop_factor%swcf
      state%crop%swcf = swcf   ! [SS-GR-ATM A5.1] runtime dual-write
      if (swcf == 1) then
         if (allocated(cfg%crop_factor%cftb)) then
            block
               integer :: nr, j
               nr = size(cfg%crop_factor%cftb, 1)
               do j = 1, nr
                  cftb(j*2-1) = cfg%crop_factor%cftb(j,1)
                  cftb(j*2)   = cfg%crop_factor%cftb(j,2)
               end do
               chtb = -99.99d0
               state%crop%fixed%cftb = cftb   ! [SS-GR-CROP A5.2]
               state%crop%fixed%chtb = chtb   ! [SS-GR-CROP A5.2]
            end block
         end if
         ! ETref standard defaults for albedo/rsc/rsw
         albedo = 0.23_real64
         rsc    = 70.0_real64
         rsw    = 0.0_real64
      else if (swcf == 2) then
         albedo = cfg%crop_factor%albedo
         rsc    = cfg%crop_factor%rsc
         rsw    = cfg%crop_factor%rsw
         if (allocated(cfg%crop_factor%chtb)) then
            block
               integer :: nr, j
               nr = size(cfg%crop_factor%chtb, 1)
               do j = 1, nr
                  chtb(j*2-1) = cfg%crop_factor%chtb(j,1)
                  chtb(j*2)   = cfg%crop_factor%chtb(j,2)
               end do
               cftb = -99.99d0
               state%crop%fixed%chtb = chtb   ! [SS-GR-CROP A5.2]
               state%crop%fixed%cftb = cftb   ! [SS-GR-CROP A5.2]
            end block
         end if
      end if
      state%crop%common%albedo = albedo   ! [SS-GR-CROPRT A5]
      state%crop%common%rsc    = rsc      ! [SS-GR-CROPRT A5]

      ! Part 14: interception (readwofost line 2640-2642)
      swinter = cfg%interception%swinter
      if (swinter == 1) then
         cofab = cfg%interception%cofab
         state%crop%cofab = cofab   ! [SS-GR-ATM A5.1] runtime dual-write
      end if

      ! Part 2: phenology (soybean=0 path; readwofost lines 2712-2721)
      idsl   = cfg%phenology%idsl
      tsumea = cfg%phenology%tsumea
      tsumam = cfg%phenology%tsumam
      if (idsl == 1 .or. idsl == 2) then
         dlo = cfg%phenology%dlo
         dlc = cfg%phenology%dlc
      end if
      if (allocated(cfg%phenology%dtsmtb)) then
         block
            integer :: nr, j
            nr = size(cfg%phenology%dtsmtb, 1)
            do j = 1, nr
               dtsmtb(j*2-1) = cfg%phenology%dtsmtb(j,1)
               dtsmtb(j*2)   = cfg%phenology%dtsmtb(j,2)
            end do
         end block
      end if
      ! vernalization (idsl=2) stub-guarded above — safe defaults already in module

      ! Part 3: initial crop state (readwofost lines 2752-2754)
      tdwi   = cfg%initial%tdwi
      laiem  = cfg%initial%laiem
      rgrlai = cfg%initial%rgrlai

      ! Part 4: green area (readwofost lines 2757-2761)
      if (allocated(cfg%green_area%slatb)) then
         block
            integer :: nr, j
            nr = size(cfg%green_area%slatb, 1)
            do j = 1, nr
               slatb(j*2-1) = cfg%green_area%slatb(j,1)
               slatb(j*2)   = cfg%green_area%slatb(j,2)
            end do
         end block
      end if
      spa   = cfg%green_area%spa
      ssa   = cfg%green_area%ssa
      span  = cfg%green_area%span
      tbase = cfg%green_area%tbase

      ! Part 5: assimilation (readwofost lines 2764-2769)
      kdif = cfg%assimilation%kdif
      kdir = cfg%assimilation%kdir
      state%crop%kdif = kdif   ! [SS-GR-ATM A5.1] runtime dual-write
      state%crop%kdir = kdir   ! [SS-GR-ATM A5.1] runtime dual-write
      eff  = cfg%assimilation%eff
      if (allocated(cfg%assimilation%amaxtb)) then
         block
            integer :: nr, j
            nr = size(cfg%assimilation%amaxtb, 1)
            do j = 1, nr
               amaxtb(j*2-1) = cfg%assimilation%amaxtb(j,1)
               amaxtb(j*2)   = cfg%assimilation%amaxtb(j,2)
            end do
         end block
      end if
      if (allocated(cfg%assimilation%tmpftb)) then
         block
            integer :: nr, j
            nr = size(cfg%assimilation%tmpftb, 1)
            do j = 1, nr
               tmpftb(j*2-1) = cfg%assimilation%tmpftb(j,1)
               tmpftb(j*2)   = cfg%assimilation%tmpftb(j,2)
            end do
         end block
      end if
      if (allocated(cfg%assimilation%tmnftb)) then
         block
            integer :: nr, j
            nr = size(cfg%assimilation%tmnftb, 1)
            do j = 1, nr
               tmnftb(j*2-1) = cfg%assimilation%tmnftb(j,1)
               tmnftb(j*2)   = cfg%assimilation%tmnftb(j,2)
            end do
         end block
      end if

      ! Part 6: conversion (readwofost lines 2772-2775)
      cvl = cfg%conversion%cvl
      cvo = cfg%conversion%cvo
      cvr = cfg%conversion%cvr
      cvs = cfg%conversion%cvs

      ! Part 7: respiration (readwofost lines 2778-2783)
      q10 = cfg%respiration%q10
      rml = cfg%respiration%rml
      rmo = cfg%respiration%rmo
      rmr = cfg%respiration%rmr
      rms = cfg%respiration%rms
      if (allocated(cfg%respiration%rfsetb)) then
         block
            integer :: nr, j
            nr = size(cfg%respiration%rfsetb, 1)
            do j = 1, nr
               rfsetb(j*2-1) = cfg%respiration%rfsetb(j,1)
               rfsetb(j*2)   = cfg%respiration%rfsetb(j,2)
            end do
         end block
      end if

      ! Part 8: partitioning (readwofost lines 2786-2789)
      if (allocated(cfg%partitioning%frtb)) then
         block
            integer :: nr, j
            nr = size(cfg%partitioning%frtb, 1)
            do j = 1, nr
               frtb(j*2-1) = cfg%partitioning%frtb(j,1)
               frtb(j*2)   = cfg%partitioning%frtb(j,2)
            end do
         end block
      end if
      if (allocated(cfg%partitioning%fltb)) then
         block
            integer :: nr, j
            nr = size(cfg%partitioning%fltb, 1)
            do j = 1, nr
               fltb(j*2-1) = cfg%partitioning%fltb(j,1)
               fltb(j*2)   = cfg%partitioning%fltb(j,2)
            end do
         end block
      end if
      if (allocated(cfg%partitioning%fstb)) then
         block
            integer :: nr, j
            nr = size(cfg%partitioning%fstb, 1)
            do j = 1, nr
               fstb(j*2-1) = cfg%partitioning%fstb(j,1)
               fstb(j*2)   = cfg%partitioning%fstb(j,2)
            end do
         end block
      end if
      if (allocated(cfg%partitioning%fotb)) then
         block
            integer :: nr, j
            nr = size(cfg%partitioning%fotb, 1)
            do j = 1, nr
               fotb(j*2-1) = cfg%partitioning%fotb(j,1)
               fotb(j*2)   = cfg%partitioning%fotb(j,2)
            end do
         end block
      end if

      ! Part 9: death rates (readwofost lines 2792-2794)
      perdl = cfg%death%perdl
      if (allocated(cfg%death%rdrrtb)) then
         block
            integer :: nr, j
            nr = size(cfg%death%rdrrtb, 1)
            do j = 1, nr
               rdrrtb(j*2-1) = cfg%death%rdrrtb(j,1)
               rdrrtb(j*2)   = cfg%death%rdrrtb(j,2)
            end do
         end block
      end if
      if (allocated(cfg%death%rdrstb)) then
         block
            integer :: nr, j
            nr = size(cfg%death%rdrstb, 1)
            do j = 1, nr
               rdrstb(j*2-1) = cfg%death%rdrstb(j,1)
               rdrstb(j*2)   = cfg%death%rdrstb(j,2)
            end do
         end block
      end if

      ! Part 11: oxygen stress (readwofost lines 2797-2867)
      swoxygen   = cfg%oxygen_stress%swoxygen
      swWrtNonox = cfg%oxygen_stress%swwrtnonox
      aeratecrit = cfg%oxygen_stress%aeratecrit
      if (swoxygen == 1) then
         hlim1  = cfg%oxygen_stress%hlim1
         hlim2u = cfg%oxygen_stress%hlim2u
         hlim2l = cfg%oxygen_stress%hlim2l
      end if

      ! Part 12: drought stress (readwofost lines 2870-2881)
      swdrought = cfg%drought_stress%swdrought
      if (swdrought == 1) then
         hlim3h = cfg%drought_stress%hlim3h
         hlim3l = cfg%drought_stress%hlim3l
         hlim4  = cfg%drought_stress%hlim4
         adcrh  = cfg%drought_stress%adcrh
         adcrl  = cfg%drought_stress%adcrl
      end if

      ! Part 13: salinity stress (readwofost lines 2900-2924; gated on flsolute
      ! at runtime — here we always copy what was validated)
      swsalinity = cfg%salinity%swsalinity
      if (swsalinity == 1) then
         saltmax   = cfg%salinity%saltmax
         saltslope = cfg%salinity%saltslope
      end if

      ! Part xx: compensation (readwofost lines 2927-2976)
      swcompensate = cfg%compensate%swcompensate
      ! swstressor is only meaningful when swcompensate > 0; default=1 from variables
      ! module init. Only set when enabled to avoid clobbering with 0 default.
      if (swcompensate > 0) swstressor = cfg%compensate%swstressor

      ! Part 10: root depth and density (readwofost lines 2993-3031)
      swrdc = cfg%root%swrdc
      if (allocated(cfg%root%rdctb)) then
         block
            integer :: nr, j
            nr = size(cfg%root%rdctb, 1)
            do j = 1, nr
               rdctb(j*2-1) = cfg%root%rdctb(j,1)
               rdctb(j*2)   = cfg%root%rdctb(j,2)
            end do
         end block
      end if
      swrd = cfg%root%swrd
      select case (swrd)
      case (1)
         if (allocated(cfg%root%rdtb)) then
            block
               integer :: nr, j
               nr = size(cfg%root%rdtb, 1)
               do j = 1, nr
                  rdtb(j*2-1) = cfg%root%rdtb(j,1)
                  rdtb(j*2)   = cfg%root%rdtb(j,2)
               end do
            end block
         end if
      case (2)
         rdi      = cfg%root%rdi
         rri      = cfg%root%rri
         rdc      = cfg%root%rdc
         swdmi2rd = cfg%root%swdmi2rd
         state%crop%common%rdi = rdi   ! [SS-GR-CROP A5.2]
         state%crop%common%rri = rri   ! [SS-GR-CROP A5.2]
         state%crop%common%rdc = rdc   ! [SS-GR-CROP A5.2]
      case (3)
         if (allocated(cfg%root%rlwtb)) then
            block
               integer :: nr, j
               nr = size(cfg%root%rlwtb, 1)
               do j = 1, nr
                  rlwtb(j*2-1) = cfg%root%rlwtb(j,1)
                  rlwtb(j*2)   = cfg%root%rlwtb(j,2)
               end do
            end block
         end if
         wrtmax = cfg%root%wrtmax
      end select

      ! Harvest (readwofost lines 2745-2748)
      dvsend = cfg%harvest%dvsend
      swharv = cfg%harvest%swharv

      ! Schedule (schedule=0 only on TOML path; stub-guarded above)
      schedule = cfg%schedule%schedule

      ! Management (readwofost lines 2979-2988)
      relmf      = cfg%management%relmf
      swpotrelmf = cfg%management%swpotrelmf
      state%crop%grass%relmf      = relmf      ! [SS-GR-CROP A5.2]
      state%crop%grass%swpotrelmf = swpotrelmf ! [SS-GR-CROP A5.2]

      ! FraDeceasedLvToSoil — local SAVE in wofost(), returned via intent(out)
      ! so the dispatch block can assign it.  (FraHarLosOrm_* are set by
      ! apply_cropwofost_nutrient below when flcropnut=true on this rotation.)
      FraDeceasedLvToSoil = cfg%management%fradeceasedlvtosoil

      ! ----------------------------------------------------------------
      ! Runtime init math — cumdens computation (readwofost lines 3094-3121)
      ! Only when swdrought=1 (Feddes). Identical algorithm to readwofost.
      ! ----------------------------------------------------------------
      if (swdrought == 1) then
         ! Build rootdis array: 101 points from 0.0 to 1.0
         do i = 0, 100
            depth = 0.01d0 * dble(i)
            rootdis(i*2+1) = depth
            rootdis(i*2+2) = afgen(rdctb, 22, depth)
         end do

         ! Copy depths to odd cumdens indices
         do i = 1, 202, 2
            cumdens(i) = rootdis(i)
         end do

         ! Trapezoidal cumulative integration into even indices
         sum_dens   = 0.0d0
         cumdens(2) = 0.0d0
         do i = 4, 202, 2
            sum_dens = sum_dens + (rootdis(i-2) + rootdis(i)) * 0.5d0 &
                                * (cumdens(i-1) - cumdens(i-3))
            cumdens(i) = sum_dens
         end do

         ! Normalize to 1
         if (sum_dens > 0.0d0) then
            do i = 2, 202, 2
               cumdens(i) = cumdens(i) / sum_dens
            end do
         end if
         state%crop%common%cumdens = cumdens   ! [SS-GR-CROPRT A5]
      end if

      ! ----------------------------------------------------------------
      ! Crop state scalars — zeroed at runtime init (readwofost is called
      ! after InitializeCrop which already zeros tsum/daycrop; we mirror
      ! here for completeness and in case the call sequence changes).
      ! .END-file restart (swinco=3) is NOT reproduced; stub-errored by
      ! the validator for Phase 2.
      ! ----------------------------------------------------------------
      state%crop%common%dvs     = 0.0d0
      tsum    = 0.0d0
      daycrop = 0
      nofd    = 0
      state%crop%common%tsum    = tsum    ! [SS-GR-CROP A5.2]
      state%crop%common%daycrop = daycrop ! [SS-GR-CROP A5.2]

      ! [nutrients] N3: drive the legacy global flCropNut from the
      ! per-rotation typed config. cropwofost_init_from_config runs at
      ! every rotation start, so a sequence of rotations with mixed
      ! flcropnut values toggles the gate correctly.
      flCropNut = cfg%nutrient%flcropnut
      state%crop%common%flCropNut = flCropNut ! [SS-GR-CROP A5.2]
      if (flCropNut) call apply_cropwofost_nutrient(cfg%nutrient)

   end subroutine cropwofost_init_from_config


   !> Apply the per-rotation [wofost.nutrient] config to the legacy
   !! `variables` globals. Called from cropwofost_init_from_config (or
   !! directly from an alternative entry point) when
   !! cfg%nutrient%flcropnut = .true..
   !!
   !! Replaces the deleted rdinit/rdsdou block in cropgrowth.f90's wofost
   !! subroutine (legacy readers physical deletion arc, SS-C step 2).
   !!
   !! See ADR 0025 ([nutrients] N1).
   subroutine apply_cropwofost_nutrient(cfg)
      use cropwofost_config_mod, only: wofost_nutrient_t
      ! [SS-GR-FINAL B7] DEFERRED: nutrient config globals — config→globals copy; Phase C3
      use variables, only: &
                           ! DEFERRED: nutrient parameters; Phase C3 (nutrients_state_t migration)
                           lrnr, lsnr, nlue, rnflv, rnfst, frnx, nmxlv,            &
                           nlai, nmaxso, npart, nfixf, nsla, rnfrt, tcnt,           &
                           dvsnlt, dvsnt, rdrns, fntrt, ilnmxl,                     &
                           fraharlosorm_lv, fraharlosorm_st, fraharlosorm_so
      type(wofost_nutrient_t), intent(in) :: cfg

      integer :: n

      ! Module-level scalars (already exist in module variables)
      lrnr   = cfg%lrnr
      lsnr   = cfg%lsnr
      nlue   = cfg%nlue
      rnflv  = cfg%rnflv
      rnfst  = cfg%rnfst
      frnx   = cfg%frnx

      ! Newly-promoted module variables (Task 1)
      nlai   = cfg%nlai
      nmaxso = cfg%nmaxso
      npart  = cfg%npart
      nfixf  = cfg%nfixf
      nsla   = cfg%nsla
      rnfrt  = cfg%rnfrt
      tcnt   = cfg%tcnt
      dvsnlt = cfg%dvsnlt
      dvsnt  = cfg%dvsnt
      rdrns  = cfg%rdrns
      fntrt  = cfg%fntrt

      ! NMXLV array — copy entries; ILNMXL records the active length
      n = 0
      if (allocated(cfg%nmxlv)) n = size(cfg%nmxlv)
      ilnmxl = n
      nmxlv  = 0.0_real64
      if (n > 0) nmxlv(1:n) = cfg%nmxlv(1:n)

      ! Harvest fractions
      fraharlosorm_lv = cfg%frahar_los_orm_lv
      fraharlosorm_st = cfg%frahar_los_orm_st
      fraharlosorm_so = cfg%frahar_los_orm_so
   end subroutine apply_cropwofost_nutrient

end module cropwofost_init_mod
