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

   ! [GR-CROP 2026-05-25] cropwofost shared per-rotation config snapshots.
   ! read by cropwofost_runtime_mod%wofost (which `use`s this module).
   ! [ADR 0052] The nutrient + harvest-loss cw_* clusters were removed with the
   ! WOFOST-N detachment; only the vernalisation cluster (non-nutrient phenology)
   ! remains.
   ! Vernalisation cluster:
   real(real64), public, save :: cw_vernbase = 0.0_real64
   real(real64), public, save :: cw_verndvs  = 0.0_real64
   real(real64), public, save :: cw_vernsat  = 0.0_real64
   real(real64), public, save :: cw_vernrtb(30) = 0.0_real64  ! 2-col table flat slice (size 30 to match legacy)

contains

   subroutine cropwofost_init_from_config(cfg, icrop, FraDeceasedLvToSoil, state)
      ! [GR-CROP 2026-05-25] cropwofost_init writes phenology config (idsl/dlo/dlc)
      ! directly to state%crop%wofost; daycrop/flCropNut written to state%crop%common.
      ! File is now `use variables`-free.
      use array_utils, only: afgen
      use error_mod,   only: fatalerr_collected
      use swap_state_mod, only: swap_state_t
      implicit none

      type(cropwofost_config_t), intent(in)    :: cfg
      integer,                   intent(in)    :: icrop  ! rotation slot (reserved)
      real(real64),              intent(out)   :: FraDeceasedLvToSoil
      type(swap_state_t),        intent(inout) :: state

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
      state%crop%swcf = cfg%crop_factor%swcf
      if (state%crop%swcf == 1) then
         if (allocated(cfg%crop_factor%cftb)) then
            block
               integer :: nr, j
               nr = size(cfg%crop_factor%cftb, 1)
               do j = 1, nr
                  state%crop%fixed%cftb(j*2-1) = cfg%crop_factor%cftb(j,1)
                  state%crop%fixed%cftb(j*2)   = cfg%crop_factor%cftb(j,2)
               end do
               state%crop%fixed%chtb = -99.99d0
            end block
         end if
         ! ETref standard defaults for albedo/rsc/rsw
         state%crop%common%albedo = 0.23_real64
         state%crop%common%rsc    = 70.0_real64
         state%crop%common%rsw    = 0.0_real64
      else if (state%crop%swcf == 2) then
         state%crop%common%albedo = cfg%crop_factor%albedo
         state%crop%common%rsc    = cfg%crop_factor%rsc
         state%crop%common%rsw    = cfg%crop_factor%rsw
         if (allocated(cfg%crop_factor%chtb)) then
            block
               integer :: nr, j
               nr = size(cfg%crop_factor%chtb, 1)
               do j = 1, nr
                  state%crop%fixed%chtb(j*2-1) = cfg%crop_factor%chtb(j,1)
                  state%crop%fixed%chtb(j*2)   = cfg%crop_factor%chtb(j,2)
               end do
               state%crop%fixed%cftb = -99.99d0
            end block
         end if
      else if (state%crop%swcf == 3) then
         ! LAI-dependent dual crop coefficient: cf, cfeic (wet) and ch are all
         ! LAI-indexed; reflection coeffs take the ETref defaults (legacy:
         ! swcf=1 or 3). The runtime looks up cfeic = afgen(cfeictb, lai).
         state%crop%common%albedo = 0.23_real64
         state%crop%common%rsc    = 70.0_real64
         state%crop%common%rsw    = 0.0_real64
         block
            integer :: nr, j
            if (allocated(cfg%crop_factor%cftb)) then
               nr = size(cfg%crop_factor%cftb, 1)
               do j = 1, nr
                  state%crop%fixed%cftb(j*2-1) = cfg%crop_factor%cftb(j,1)
                  state%crop%fixed%cftb(j*2)   = cfg%crop_factor%cftb(j,2)
               end do
            end if
            if (allocated(cfg%crop_factor%cfeictb)) then
               nr = size(cfg%crop_factor%cfeictb, 1)
               do j = 1, nr
                  state%crop%fixed%cfeictb(j*2-1) = cfg%crop_factor%cfeictb(j,1)
                  state%crop%fixed%cfeictb(j*2)   = cfg%crop_factor%cfeictb(j,2)
               end do
            end if
            if (allocated(cfg%crop_factor%chtb)) then
               nr = size(cfg%crop_factor%chtb, 1)
               do j = 1, nr
                  state%crop%fixed%chtb(j*2-1) = cfg%crop_factor%chtb(j,1)
                  state%crop%fixed%chtb(j*2)   = cfg%crop_factor%chtb(j,2)
               end do
            end if
         end block
      end if

      ! Part 14: interception (readwofost line 2640-2642)
      state%crop%common%swinter = cfg%interception%swinter
      if (state%crop%common%swinter == 1) then
         state%crop%cofab = cfg%interception%cofab
      else if (state%crop%common%swinter == 2 .and. allocated(cfg%interception%gashtb)) then
         ! Gash: split the 6-column (t, pfree, pstem, scanopy, avprec, avevap)
         ! table into the five (t, value) atmosphere arrays that legacy
         ! readwofost populated as globals.
         block
            integer :: r, nr
            nr = size(cfg%interception%gashtb, 1)
            do r = 1, nr
               if (2*r > size(state%atmosphere%pfreetb)) exit
               state%atmosphere%pfreetb(2*r-1)   = cfg%interception%gashtb(r, 1)
               state%atmosphere%pfreetb(2*r)     = cfg%interception%gashtb(r, 2)
               state%atmosphere%pstemtb(2*r-1)   = cfg%interception%gashtb(r, 1)
               state%atmosphere%pstemtb(2*r)     = cfg%interception%gashtb(r, 3)
               state%atmosphere%scanopytb(2*r-1) = cfg%interception%gashtb(r, 1)
               state%atmosphere%scanopytb(2*r)   = cfg%interception%gashtb(r, 4)
               state%atmosphere%avprectb(2*r-1)  = cfg%interception%gashtb(r, 1)
               state%atmosphere%avprectb(2*r)    = cfg%interception%gashtb(r, 5)
               state%atmosphere%avevaptb(2*r-1)  = cfg%interception%gashtb(r, 1)
               state%atmosphere%avevaptb(2*r)    = cfg%interception%gashtb(r, 6)
            end do
         end block
      end if

      ! Part 2: phenology (soybean=0 path; readwofost lines 2712-2721)
      state%crop%wofost%idsl = cfg%phenology%idsl
      state%crop%common%tsumea = cfg%phenology%tsumea
      state%crop%common%tsumam = cfg%phenology%tsumam
      if (state%crop%wofost%idsl == 1 .or. state%crop%wofost%idsl == 2) then
         state%crop%wofost%dlo = cfg%phenology%dlo
         state%crop%wofost%dlc = cfg%phenology%dlc
      end if
      if (allocated(cfg%phenology%dtsmtb)) then
         block
            integer :: nr, j
            nr = size(cfg%phenology%dtsmtb, 1)
            do j = 1, nr
               state%crop%common%dtsmtb(j*2-1) = cfg%phenology%dtsmtb(j,1)
               state%crop%common%dtsmtb(j*2)   = cfg%phenology%dtsmtb(j,2)
            end do
         end block
      end if
      ! vernalization (idsl=2) stub-guarded above — safe defaults already in module

      ! Part 3: initial crop state (readwofost lines 2752-2754)
      state%crop%common%tdwi  = cfg%initial%tdwi
      state%crop%common%laiem = cfg%initial%laiem
      state%crop%common%rgrlai = cfg%initial%rgrlai

      ! Part 4: green area (readwofost lines 2757-2761)
      if (allocated(cfg%green_area%slatb)) then
         block
            integer :: nr, j
            nr = size(cfg%green_area%slatb, 1)
            do j = 1, nr
               state%crop%common%slatb(j*2-1) = cfg%green_area%slatb(j,1)
               state%crop%common%slatb(j*2)   = cfg%green_area%slatb(j,2)
            end do
         end block
      end if
      state%crop%common%spa  = cfg%green_area%spa
      state%crop%common%ssa  = cfg%green_area%ssa
      state%crop%common%span = cfg%green_area%span
      state%crop%common%tbase = cfg%green_area%tbase

      ! Part 5: assimilation (readwofost lines 2764-2769)
      state%crop%kdif = cfg%assimilation%kdif
      state%crop%kdir = cfg%assimilation%kdir
      state%crop%common%eff = cfg%assimilation%eff
      if (allocated(cfg%assimilation%amaxtb)) then
         block
            integer :: nr, j
            nr = size(cfg%assimilation%amaxtb, 1)
            do j = 1, nr
               state%crop%common%amaxtb(j*2-1) = cfg%assimilation%amaxtb(j,1)
               state%crop%common%amaxtb(j*2)   = cfg%assimilation%amaxtb(j,2)
            end do
         end block
      end if
      if (allocated(cfg%assimilation%tmpftb)) then
         block
            integer :: nr, j
            nr = size(cfg%assimilation%tmpftb, 1)
            do j = 1, nr
               state%crop%common%tmpftb(j*2-1) = cfg%assimilation%tmpftb(j,1)
               state%crop%common%tmpftb(j*2)   = cfg%assimilation%tmpftb(j,2)
            end do
         end block
      end if
      if (allocated(cfg%assimilation%tmnftb)) then
         block
            integer :: nr, j
            nr = size(cfg%assimilation%tmnftb, 1)
            do j = 1, nr
               state%crop%common%tmnftb(j*2-1) = cfg%assimilation%tmnftb(j,1)
               state%crop%common%tmnftb(j*2)   = cfg%assimilation%tmnftb(j,2)
            end do
         end block
      end if

      ! Part 6: conversion (readwofost lines 2772-2775)
      state%crop%common%cvl = cfg%conversion%cvl
      state%crop%common%cvo = cfg%conversion%cvo
      state%crop%common%cvr = cfg%conversion%cvr
      state%crop%common%cvs = cfg%conversion%cvs

      ! Part 7: respiration (readwofost lines 2778-2783)
      state%crop%common%q10 = cfg%respiration%q10
      state%crop%common%rml = cfg%respiration%rml
      state%crop%common%rmo = cfg%respiration%rmo
      state%crop%common%rmr = cfg%respiration%rmr
      state%crop%common%rms = cfg%respiration%rms
      if (allocated(cfg%respiration%rfsetb)) then
         block
            integer :: nr, j
            nr = size(cfg%respiration%rfsetb, 1)
            do j = 1, nr
               state%crop%common%rfsetb(j*2-1) = cfg%respiration%rfsetb(j,1)
               state%crop%common%rfsetb(j*2)   = cfg%respiration%rfsetb(j,2)
            end do
         end block
      end if

      ! Part 8: partitioning (readwofost lines 2786-2789)
      if (allocated(cfg%partitioning%frtb)) then
         block
            integer :: nr, j
            nr = size(cfg%partitioning%frtb, 1)
            do j = 1, nr
               state%crop%common%frtb(j*2-1) = cfg%partitioning%frtb(j,1)
               state%crop%common%frtb(j*2)   = cfg%partitioning%frtb(j,2)
            end do
         end block
      end if
      if (allocated(cfg%partitioning%fltb)) then
         block
            integer :: nr, j
            nr = size(cfg%partitioning%fltb, 1)
            do j = 1, nr
               state%crop%common%fltb(j*2-1) = cfg%partitioning%fltb(j,1)
               state%crop%common%fltb(j*2)   = cfg%partitioning%fltb(j,2)
            end do
         end block
      end if
      if (allocated(cfg%partitioning%fstb)) then
         block
            integer :: nr, j
            nr = size(cfg%partitioning%fstb, 1)
            do j = 1, nr
               state%crop%common%fstb(j*2-1) = cfg%partitioning%fstb(j,1)
               state%crop%common%fstb(j*2)   = cfg%partitioning%fstb(j,2)
            end do
         end block
      end if
      if (allocated(cfg%partitioning%fotb)) then
         block
            integer :: nr, j
            nr = size(cfg%partitioning%fotb, 1)
            do j = 1, nr
               state%crop%common%fotb(j*2-1) = cfg%partitioning%fotb(j,1)
               state%crop%common%fotb(j*2)   = cfg%partitioning%fotb(j,2)
            end do
         end block
      end if

      ! Part 9: death rates (readwofost lines 2792-2794)
      state%crop%common%perdl = cfg%death%perdl
      if (allocated(cfg%death%rdrrtb)) then
         block
            integer :: nr, j
            nr = size(cfg%death%rdrrtb, 1)
            do j = 1, nr
               state%crop%common%rdrrtb(j*2-1) = cfg%death%rdrrtb(j,1)
               state%crop%common%rdrrtb(j*2)   = cfg%death%rdrrtb(j,2)
            end do
         end block
      end if
      if (allocated(cfg%death%rdrstb)) then
         block
            integer :: nr, j
            nr = size(cfg%death%rdrstb, 1)
            do j = 1, nr
               state%crop%common%rdrstb(j*2-1) = cfg%death%rdrstb(j,1)
               state%crop%common%rdrstb(j*2)   = cfg%death%rdrstb(j,2)
            end do
         end block
      end if

      ! Part 11: oxygen stress (readwofost lines 2797-2867)
      state%crop%common%swoxygen = cfg%oxygen_stress%swoxygen
      state%crop%common%swWrtNonox = cfg%oxygen_stress%swwrtnonox
      state%crop%common%aeratecrit = cfg%oxygen_stress%aeratecrit
      if (state%crop%common%swoxygen == 1) then
         state%crop%common%hlim1  = cfg%oxygen_stress%hlim1
         state%crop%common%hlim2u = cfg%oxygen_stress%hlim2u
         state%crop%common%hlim2l = cfg%oxygen_stress%hlim2l
      else if (state%crop%common%swoxygen == 2) then
         ! Bartholomeus physical sub-model (mirrors cropgrass_init; the shared
         ! oxygenstress.f90 kernel reads these state fields). WOFOST has no
         ! swoxygentype field in its crp — legacy readwofost defaults it to 1
         ! (physical OxygenStress); swoxygentype=2 (reproduction functions) is the
         ! grass-only alternative and is dormant in rootextraction.f90 anyway.
         state%crop%common%swoxygentype        = 1
         state%crop%oxygen%q10_microbial       = cfg%oxygen_stress%q10_microbial
         state%crop%oxygen%specific_resp_humus = cfg%oxygen_stress%specific_resp_humus
         state%crop%common%srl                 = cfg%oxygen_stress%srl
         state%crop%common%swrootradius        = cfg%oxygen_stress%swrootradius
         if (cfg%oxygen_stress%swrootradius == 1) then
            state%crop%common%dry_mat_cont_roots      = cfg%oxygen_stress%dry_mat_cont_roots
            state%crop%common%air_filled_root_por     = cfg%oxygen_stress%air_filled_root_por
            state%crop%common%spec_weight_root_tissue = cfg%oxygen_stress%spec_weight_root_tissue
            state%crop%common%var_a                   = cfg%oxygen_stress%var_a
         else
            state%crop%common%root_radiusO2 = cfg%oxygen_stress%root_radiusO2
         end if
      end if

      ! Part 12: drought stress (readwofost lines 2870-2881)
      state%crop%common%swdrought = cfg%drought_stress%swdrought
      if (state%crop%common%swdrought == 1) then
         state%crop%common%hlim3h = cfg%drought_stress%hlim3h
         state%crop%common%hlim3l = cfg%drought_stress%hlim3l
         state%crop%common%hlim4  = cfg%drought_stress%hlim4
         state%crop%common%adcrh  = cfg%drought_stress%adcrh
         state%crop%common%adcrl  = cfg%drought_stress%adcrl
      end if

      ! Part 13: salinity stress (readwofost lines 2900-2924; gated on flsolute
      ! at runtime — here we always copy what was validated)
      state%crop%common%swsalinity = cfg%salinity%swsalinity
      if (state%crop%common%swsalinity == 1) then
         state%crop%common%saltmax   = cfg%salinity%saltmax
         state%crop%common%saltslope = cfg%salinity%saltslope
      end if

      ! Part xx: compensation (readwofost lines 2927-2976)
      state%crop%common%swcompensate = cfg%compensate%swcompensate
      ! swstressor is only meaningful when swcompensate > 0; default=1 from variables
      ! module init. Only set when enabled to avoid clobbering with 0 default.
      if (state%crop%common%swcompensate > 0) state%crop%common%swstressor = cfg%compensate%swstressor
      ! Legacy readwofost reads alphacrit for Jarvis (swcompensate=1) and dcritrtz
      ! for Walsum (swcompensate=2; alphacrit is then derived in RootExtraction).
      if (cfg%compensate%swcompensate == 1) state%crop%common%alphacrit = cfg%compensate%alphacrit
      if (cfg%compensate%swcompensate == 2) state%crop%common%dcritrtz  = cfg%compensate%dcritrtz

      ! Part 10: root depth and density (readwofost lines 2993-3031)
      state%crop%common%swrdc = cfg%root%swrdc
      if (allocated(cfg%root%rdctb)) then
         block
            integer :: nr, j
            nr = size(cfg%root%rdctb, 1)
            do j = 1, nr
               state%crop%common%rdctb(j*2-1) = cfg%root%rdctb(j,1)
               state%crop%common%rdctb(j*2)   = cfg%root%rdctb(j,2)
            end do
         end block
      end if
      state%crop%common%swrd = cfg%root%swrd
      select case (state%crop%common%swrd)
      case (1)
         if (allocated(cfg%root%rdtb)) then
            block
               integer :: nr, j
               nr = size(cfg%root%rdtb, 1)
               do j = 1, nr
                  state%crop%common%rdtb(j*2-1) = cfg%root%rdtb(j,1)
                  state%crop%common%rdtb(j*2)   = cfg%root%rdtb(j,2)
               end do
            end block
         end if
      case (2)
         state%crop%common%rdi = cfg%root%rdi
         state%crop%common%rri = cfg%root%rri
         state%crop%common%rdc = cfg%root%rdc
         state%crop%common%swdmi2rd = cfg%root%swdmi2rd
      case (3)
         if (allocated(cfg%root%rlwtb)) then
            block
               integer :: nr, j
               nr = size(cfg%root%rlwtb, 1)
               do j = 1, nr
                  state%crop%common%rlwtb(j*2-1) = cfg%root%rlwtb(j,1)
                  state%crop%common%rlwtb(j*2)   = cfg%root%rlwtb(j,2)
               end do
            end block
         end if
         state%crop%common%wrtmax = cfg%root%wrtmax
      end select

      ! Harvest (readwofost lines 2745-2748)
      state%crop%common%dvsend = cfg%harvest%dvsend
      state%crop%common%swharv = cfg%harvest%swharv

      ! Schedule (schedule=0 only on TOML path; stub-guarded above)
      state%crop%common%schedule = cfg%schedule%schedule

      ! Management (readwofost lines 2979-2988)
      state%crop%grass%relmf      = cfg%management%relmf
      state%crop%grass%swpotrelmf = cfg%management%swpotrelmf

      ! FraDeceasedLvToSoil — local SAVE in wofost(), returned via intent(out)
      ! so the dispatch block can assign it.  (FraHarLosOrm_* are set by
      ! apply_cropwofost_nutrient below when flcropnut=true on this rotation.)
      FraDeceasedLvToSoil = cfg%management%fradeceasedlvtosoil

      ! ----------------------------------------------------------------
      ! Runtime init math — cumdens computation (readwofost lines 3094-3121)
      ! Only when swdrought=1 (Feddes). Identical algorithm to readwofost.
      ! ----------------------------------------------------------------
      if (state%crop%common%swdrought == 1) then
         ! Build rootdis array: 101 points from 0.0 to 1.0
         do i = 0, 100
            depth = 0.01d0 * dble(i)
            rootdis(i*2+1) = depth
            rootdis(i*2+2) = afgen(state%crop%common%rdctb, 22, depth)
         end do

         ! Copy depths to odd cumdens indices
         do i = 1, 202, 2
            state%crop%common%cumdens(i) = rootdis(i)
         end do

         ! Trapezoidal cumulative integration into even indices
         sum_dens   = 0.0d0
         state%crop%common%cumdens(2) = 0.0d0
         do i = 4, 202, 2
            sum_dens = sum_dens + (rootdis(i-2) + rootdis(i)) * 0.5d0 &
                                * (state%crop%common%cumdens(i-1) - state%crop%common%cumdens(i-3))
            state%crop%common%cumdens(i) = sum_dens
         end do

         ! Normalize to 1
         if (sum_dens > 0.0d0) then
            do i = 2, 202, 2
               state%crop%common%cumdens(i) = state%crop%common%cumdens(i) / sum_dens
            end do
         end if
      end if

      ! ----------------------------------------------------------------
      ! Crop state scalars — zeroed at runtime init (readwofost is called
      ! after InitializeCrop which already zeros tsum/daycrop; we mirror
      ! here for completeness and in case the call sequence changes).
      ! .END-file restart (swinco=3) is NOT reproduced; stub-errored by
      ! the validator for Phase 2.
      ! ----------------------------------------------------------------
      state%crop%common%dvs     = 0.0d0
      state%crop%common%tsum    = 0.0d0
      state%crop%common%daycrop = 0
      state%atmosphere%nofd     = 0  ! [GR-CROP 2026-05-25] nofd retired → state%atmosphere

      ! [ADR 0052] flCropNut / apply_cropwofost_nutrient removed (WOFOST-N detached).

   end subroutine cropwofost_init_from_config


end module cropwofost_init_mod
