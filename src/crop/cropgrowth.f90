! File VersionID:
!   $Id: cropgrowth.f90 380 2018-05-28 14:01:08Z heine003 $
!
! Phase 4e Task A2 audit: all `call fatalerr` sites in this file are in
! surviving physics dispatchers (CropGrowth, cropfixed, cropoutput,
! ArableLandGerm, FacCO2, wofost, grass, astro, outbalcrop*, chckcbl).
! The crop sub-readers (readwofost, readcropfixed, readgrass) live in
! src/io/readswap.f90 — none of them are here. Phase 4e Task A4/A5
! replaces all of this file's calls via fatalerr_collected (singleton).
! [SS-GR-CROPWS A1]: cropgrowth.f90 audit — state already intent(inout) non-optional;
!   all tracked write sites carry dual-writes from prior arcs; no optional/present guards.
!   No file changes required for Phase A.
! [GR-CROPWS B3]: reads migrated — croptype(state%crop%common%icrop), icrop (block section),
!   cropstart (post-mirror), cropend (task 4), rd (root loop), dvs (task 4, post-task-3 mirror),
!   swbulb→state%crop%wofost%swbulb, plwt→state%crop%wofost%plwt, kdif→state%crop%kdif.
!   daycrop NOT migrated: InitializeCrop zeroes legacy global but does not mirror to state,
!     so state%crop%common%daycrop would be stale on first crop day.
!   lai/laipot NOT migrated (ordering not confirmed safe for all crop types).
! ----------------------------------------------------------------------
      subroutine CropGrowth(task, tsoil, state)
! ----------------------------------------------------------------------
!     UpDate             : May 2014
!     Date               : Aug 2004
!     Purpose            : Call proper crop routines for initialization,
!                          calculation of rate/state variables and output
! SS-HEAT pre-Task-8: tsoil passed as non-optional arg (from state%heat%tsoil
!   at caller); threaded down to ArableLandGerm / grass / sumttd.
! SS-CRP Phase 1 C-1.3: state added (intent inout) for dual-write of
!   hroot/hleaf/mfluxtable into state%soilwater on task=1.
!   tsoil retained: still threaded to ArableLandGerm / grass / sumttd.
! ----------------------------------------------------------------------

      ! [SS-GR-CROPRT B1] DEFERRED — all remaining variables globals for CropGrowth:
      !   icrop, flCropCalendar, cropstart, cropend, flCropEmergence, flCropHarvest,
      !   flCropReadFile, flCropPrep, flCropSow, flCropGerm: crop calendar flags, no state home
      !   swinco, croptype, swcrp, swdrought, swharv: config switches, no state home [swbulb→state B3]
      !   daycrop, rdpot, lai, laipot, cf, ch, tsum, dvs: runtime crop state (dual-write
      !     to state%crop%common%; global still needed pending Phase C global retirement)
      !   [daycrop: NOT migrated — InitializeCrop zeroes global but not state mirror]
      !   cwdmpot, cwdm, wsopot, wso, wlvpot, wlv, wstpot, wst, wrtpot, wrt: WOFOST pools
      !   tmn, lat, rad: meteo scalars, no state%atmosphere scalar home
      !   eff, amaxtb, tmpftb, tmnftb: physiology params/tables, no state home [kdif→state B3]
      !   remoc, pld, q10, pgasspot, pgass: physiology params [plwt→state B3; swbulb→state B3]
      !   flCropNut, nlue, anlv, anst, nmxlv, nmaxlv, nmaxst, nmaxrt, lrnr, lsnr, nni,
      !     rnflv, rnfst, frnx, fstr: nutrient state/params, no nutrient_state home
      !   flHarvestDay: harvest flag (dual-write to state%crop%common%; global still needed)
      !   noddrz, pathcrop, cropfil: no state home
      !   bgerm, cgerm, agerm, hprep, dhPrep, zPrep, hSow, dhSow, zSow, zTempSow,
      !     dtempSow, TempSow, MaxPrepDelay, MaxSowDelay, PrepDelay, SowDelay: germ params
      !   tsumemeopt, tsumgerm, hdrygerm, hwetgerm, zgerm, TBASEM, TEFFMX: germ thresholds
      !   atmtr, daylp, difpp, dsinbe: astro outputs used by both CropGrowth and wofost
      !   dvsend: harvest DVS threshold, no state home
      !   tsoil: config-staging buffer, renamed to avoid clash with dummy arg
      ! MIGRATED B1: relmf → state%crop%grass%relmf (read-only in CropGrowth)
      ! MIGRATED B1: swpotrelmf → state%crop%grass%swpotrelmf (read-only in CropGrowth)
      ! MIGRATED B6: fco2amax/fco2eff/fco2tra → state%crop%wofost%X (FacCO2 now writes to state;
      !   body reads at lines 459-461 use state; redundant dual-writes at old-395-397 removed)
      use variables, only: &                                             ! [SS-GR-CROPRT B1/B6] [GR-CROPWS B3]
        icrop, flCropCalendar, cropstart, cropend, flCropEmergence,         &
        flCropHarvest, flCropReadFile, flCropPrep, flCropSow, flCropGerm,   &
        swinco, croptype, daycrop, cf, ch,         &  ! tsum/rd/lai/rdpot retired
        cwdmpot, cwdm, wsopot, wlvpot, wstpot,                       &  ! wso/wst/wlv retired
        wrtpot, tmn, lat, rad,                                         &  ! wrt retired
        albedo, rsc, cumdens,                                               &
        eff, amaxtb, tmpftb, tmnftb, swdrought, swcrp, dvsend,             &  ! [GR-CROPWS B3] kdif removed (→state%crop%kdif)
        swharv, plwt, remoc, pld, q10, pgasspot, pgass,                    &  ! [GR-CROPWS B3] swbulb removed (→state%crop%wofost%swbulb)
        flCropNut, nlue, anlv, anst, nmxlv, nmaxlv, nmaxst,               &
        nmaxrt, lrnr, lsnr, nni, rnflv, rnfst, frnx, fstr, flHarvestDay,  &
        noddrz, pathcrop, cropfil, bgerm, cgerm,                           &
        agerm, hprep, dhPrep, zPrep, hSow, dhSow, zSow, zTempSow,         &
        dtempSow, TempSow, MaxPrepDelay, MaxSowDelay, PrepDelay, SowDelay, &
        tsumemeopt, tsumgerm, hdrygerm, hwetgerm, zgerm, TBASEM, TEFFMX,  &
        atmtr, daylp, difpp, dsinbe,                                       &
        dummy_tsoil_cg_ => tsoil
      !! Rename config-staging tsoil to avoid clash with dummy arg tsoil.
      !! [SS-HEAT] Task 9: tsoil retained as config-staging buffer; global is not compute state.
      use array_utils, only: afgen
      use rootextraction_mod, only: MatricFlux
      use swap_constants, only: tiny
      use error_mod, only: fatalerr_collected
      use swap_state_mod, only: swap_state_t
      ! GR-CROPWS Phase 0: helpers extracted to cropgrowth_helpers_mod
      use cropgrowth_helpers_mod, only: nocrop, ArableLandGerm, FacCO2, &
                                         CropOutput, update_rootdistribution
      ! GR-CROPWS Phase 0.2: wofost extracted to cropwofost_runtime_mod
      use cropwofost_runtime_mod, only: wofost
      ! GR-CROPWS Phase 0.3: grass extracted to cropgrass_runtime_mod
      use cropgrass_runtime_mod, only: grass
      ! GR-CROPWS Phase 0.4: cropfixed extracted to cropfixed_runtime_mod
      use cropfixed_runtime_mod, only: cropfixed
      implicit none

      ! Explicit interfaces for non-module subs that now take tsoil(:)
      ! ArableLandGerm removed — now a module procedure in cropgrowth_helpers_mod.
      ! wofost removed — now a module procedure in cropwofost_runtime_mod.
      ! grass removed — now a module procedure in cropgrass_runtime_mod.
      ! cropfixed removed — now a module procedure in cropfixed_runtime_mod.

      integer task
      real(8), intent(in) :: tsoil(:)
      !! Soil temperature array from state%heat%tsoil, passed by caller.
      type(swap_state_t), intent(inout) :: state
      !! Full typed state record; state%soilwater%hleaf/hroot/mfluxtable written on task=1.
      integer i, node
      real(8) sumtmin
      real(8) tmnr          ! running 7-day average of min temperature; local (not global) [SS-GR-FINAL B5]
      real(8) dummy_mf_   ! dummy outcome arg for MatricFlux(1) init call

      ! assimilation
      real(8) dayl,cosld,sinld
      real(8) effc,amax
      real(8) dtgapot,dtga
      ! only for bulb crops (tulips etc..)
      real(8) respmo,decrmo,remo,factblb

      ! SS-TC TC-10: t1900, daynr read via state%timecontrol tc_* aliases.
      associate( &
        tc_t1900 => state%timecontrol%t1900,  &  ! TC-10
        tc_daynr => state%timecontrol%daynr,  &  ! TC-10
        at_tavd  => state%atmosphere%tavd     &  ! [SS-GR-ATM B.5] daytime-mean temp read from state
      )

      select case (task)

      case (1)

! === initialization ==========================================================

      ! find active crop
      icrop = 1
      flCropCalendar = .false.
      do while (.not. flCropCalendar)

        if (cropstart(icrop) .lt. 1.d0) exit

        if (tc_t1900 - cropstart(icrop) .gt. -tiny                      &
     &                 .and. tc_t1900 - cropend(icrop) .lt. tiny) then
          flCropCalendar = .true.
        else
          icrop = icrop + 1
        endif
      enddo
      state%crop%common%icrop        = icrop           ! [SS-GR-CROP A5.1]
      state%crop%common%flCropCalendar = flCropCalendar  ! [SS-GR-CROP A5.1]
      ! [SS-GR-CROPRT A5] mirror current-crop window scalars
      if (flCropCalendar) then
        state%crop%common%cropstart = cropstart(icrop)
        state%crop%common%cropend   = cropend(icrop)
      end if

! --- bare soil condition  ----------------------------------------------------
      if (.not. flCropEmergence .or. flCropHarvest) then
        call nocrop (state)
        ! [SS-GR-CROP A5.1] nocrop writes state%crop%common%dvs/rd directly; mirror remaining legacy zeros
        state%crop%common%cf         = cf
        state%crop%common%ch         = ch
        state%crop%common%albedo     = albedo   ! [SS-GR-CROPRT A5]
        state%crop%common%rsc        = rsc      ! [SS-GR-CROPRT A5]
        state%crop%wofost%cwdmpot    = cwdmpot
        state%crop%wofost%cwdm       = cwdm
        state%crop%wofost%wsopot     = wsopot
        state%crop%wofost%wlvpot     = wlvpot
        state%crop%wofost%wstpot     = wstpot
        state%crop%wofost%wrtpot     = wrtpot
      endif

! --- check crop emergence ----------------------------------------------------

      ! reset if new crop
      if (flCropCalendar) then
        if (dabs(tc_t1900 - state%crop%common%cropstart) .lt. tiny) then  ! [GR-CROPWS B3]
          call InitializeCrop
          ! [SS-GR-CROPRT A5] mirror fields zeroed by InitializeCrop
          state%crop%common%flCropPrep    = flCropPrep
          state%crop%common%flCropSow     = flCropSow
          state%crop%common%flCropGerm    = flCropGerm
          state%crop%common%flCropHarvest = flCropHarvest
          state%crop%common%PrepDelay     = PrepDelay
          state%crop%common%SowDelay      = SowDelay
          state%crop%common%noddrz        = noddrz
          state%crop%common%albedo        = albedo
          state%crop%common%rsc           = rsc
          state%crop%common%cumdens       = cumdens
          flCropReadFile  = .true.
          state%crop%common%flCropReadFile = flCropReadFile   ! [SS-GR-CROPRT A5]
          flCropEmergence = .true.
          state%crop%flCropEmergence = flCropEmergence   ! [SS-GR-ATM A5.2] dual-write
          if (croptype(state%crop%common%icrop) .le. 2) then
            flCropEmergence = .false.
            state%crop%flCropEmergence = flCropEmergence   ! [SS-GR-ATM A5.2] dual-write
          endif
        endif
      endif

! --- Preparation, Sowing and Germination of arable crop growth ---------------
      if (flCropCalendar .and. .not. flCropHarvest .and. croptype(state%crop%common%icrop) .le. 2) then

        ! check crop preparation, sowing and germination (of previous day)
        if (.not. flCropEmergence) then
          if (flCropPrep .and. flCropSow .and. flCropGerm) then
            swinco          = -99
            flCropReadFile  = .true.
            state%crop%common%flCropReadFile = flCropReadFile   ! [SS-GR-CROPRT A5]
            flCropEmergence = .true.
            state%crop%flCropEmergence = flCropEmergence   ! [SS-GR-ATM A5.2] dual-write
          endif
        endif
        
        ! Initialize preparation, sowing and germination
        if (.not. flCropEmergence) then
          if (flCropReadFile) then
            ! ADR 0017: sibling-reader dispatch around legacy ArableLandGerm
            ! (which opens pathcrop//cropfil(icrop)//'.crp' via
            ! readarablelandgerm). Cache-hit path sets flCrop* flags from
            ! the typed config. Cache-miss path falls back to legacy.
            ! Phase 2 extends to type=2 (Wofost) in addition to type=1.
            ! Phase 4 extends type=2 to handle swgerm=1/2 from TOML config,
            ! mirroring legacy readarablelandgerm:3322-3358.
            ! Teardown: end of Phase 4 removes the else-branch.
            block
               use crop_config_global_mod, only: crop_config_global
               use cropwofost_config_mod, only: wofost_germination_t
               use error_mod, only: fatalerr_collected
               logical :: use_cache
               integer :: swprep_cache, swsow_cache, swgerm_cache, rot_type
               type(wofost_germination_t), pointer :: gp
               use_cache = .false.
               swprep_cache = 0
               swsow_cache  = 0
               swgerm_cache = 0
               rot_type     = 0
               nullify(gp)
               if (associated(crop_config_global)) then
                  if (allocated(crop_config_global%rotation_loaded) .and. &
                      allocated(crop_config_global%rotation_type)) then
                     if (state%crop%common%icrop >= 1 .and. state%crop%common%icrop <= size(crop_config_global%rotation_loaded)) then  ! [GR-CROPWS B3]
                        if (crop_config_global%rotation_loaded(state%crop%common%icrop)) then  ! [GR-CROPWS B3]
                           rot_type = crop_config_global%rotation_type(state%crop%common%icrop)  ! [GR-CROPWS B3]
                           select case (rot_type)
                           case (1)
                              ! type=1 cropfixed
                              if (allocated(crop_config_global%rotation_fixed)) then
                                 use_cache    = .true.
                                 swprep_cache = crop_config_global%rotation_fixed(state%crop%common%icrop)%swprep   ! [GR-CROPWS B3]
                                 swsow_cache  = crop_config_global%rotation_fixed(state%crop%common%icrop)%swsow    ! [GR-CROPWS B3]
                                 swgerm_cache = crop_config_global%rotation_fixed(state%crop%common%icrop)%swgerm   ! [GR-CROPWS B3]
                              end if
                           case (2)
                              ! type=2 wofost: cache-hit when swprep=0 AND swsow=0.
                              ! swgerm=0/1/2 are all handled in the cache body below.
                              if (allocated(crop_config_global%rotation_wofost)) then
                                 swprep_cache = crop_config_global%rotation_wofost(state%crop%common%icrop)%preparation%swprep  ! [GR-CROPWS B3]
                                 swsow_cache  = crop_config_global%rotation_wofost(state%crop%common%icrop)%sowing%swsow         ! [GR-CROPWS B3]
                                 swgerm_cache = crop_config_global%rotation_wofost(state%crop%common%icrop)%germination%swgerm   ! [GR-CROPWS B3]
                                 if (swprep_cache == 0 .and. swsow_cache == 0) then
                                    use_cache = .true.
                                    gp => crop_config_global%rotation_wofost(state%crop%common%icrop)%germination   ! [GR-CROPWS B3]
                                 end if
                              end if
                           end select
                        end if
                     end if
                  end if
               end if
               if (use_cache) then
                  ! swprep/swsow non-zero: stub-error across all types.
                  if (swprep_cache /= 0 .or. swsow_cache /= 0) then
                     call fatalerr_collected('cropgrowth/ArableLandGerm', &
                        'swprep /= 0 or swsow /= 0 not yet supported in TOML pipeline.')
                  else
                     ! Prep and sow done (both switches are 0).
                     flCropPrep = .true.
                     flCropSow  = .true.
                     PrepDelay  = 0
                     SowDelay   = 0
                     state%crop%common%flCropPrep = flCropPrep   ! [SS-GR-CROPRT A5]
                     state%crop%common%flCropSow  = flCropSow    ! [SS-GR-CROPRT A5]
                     state%crop%common%PrepDelay  = PrepDelay    ! [SS-GR-CROPRT A5]
                     state%crop%common%SowDelay   = SowDelay     ! [SS-GR-CROPRT A5]
                     if (swgerm_cache == 0) then
                        ! swgerm=0: germination and emergence are immediate.
                        flCropGerm      = .true.
                        state%crop%common%flCropGerm = flCropGerm   ! [SS-GR-CROPRT A5]
                        flCropEmergence = .true.
                        state%crop%flCropEmergence = flCropEmergence   ! [SS-GR-ATM A5.2] dual-write
                     else if (rot_type == 2) then
                        ! type=2 with swgerm=1 or 2: copy germ params from cfg,
                        ! mirror legacy readarablelandgerm:3322-3358.
                        flCropGerm      = .false.
                        state%crop%common%flCropGerm = flCropGerm   ! [SS-GR-CROPRT A5]
                        flCropEmergence = .false.
                        state%crop%flCropEmergence = flCropEmergence   ! [SS-GR-ATM A5.2] dual-write
                        tsumemeopt = gp%tsumemeopt
                        tbasem     = gp%tbasem
                        teffmx     = gp%teffmx
                        agerm      = -99.0d0
                        if (swgerm_cache == 2) then
                           hdrygerm = gp%hdrygerm
                           hwetgerm = gp%hwetgerm
                           if (gp%zgerm /= 0.0d0) then
                              zgerm = gp%zgerm
                           else
                              zgerm = -10.0d0   ! legacy default
                           end if
                           agerm = gp%agerm
                           cgerm = - (tsumemeopt - agerm * log10(-hdrygerm))
                           bgerm =   (tsumemeopt + agerm * log10(-hwetgerm))
                        end if
                     else
                        ! type=1 cropfixed with swgerm > 0 — not yet supported.
                        call fatalerr_collected('cropgrowth/ArableLandGerm', &
                           'cropfixed swgerm > 0 not yet supported in TOML pipeline.')
                     end if
                  end if
               else
                  ! ADR 0017 sibling-reader dispatch: cache-miss is now
                  ! treated as a fatal error rather than a silent legacy
                  ! fallback. Every rotation in tests/swap-cases/toml/
                  ! has a .crp.toml; missing one is a user error.
                  call fatalerr_collected('cropgrowth/ArableLandGerm', &
                     'rotation has no loaded .crp.toml — author the file or use the legacy executable.')
               end if
            end block
          endif
        endif

        if (.not. flCropEmergence .and. .not. flCropHarvest) then
          
          ! Preparation before crop growth
          if (.not. flCropPrep) then
            call ArableLandGerm(2, tsoil, state)  ! [SS-SWC S-2.7]
          endif

          ! Sowing before crop growth
          if (flCropPrep .and. .not. flCropSow) then
            call ArableLandGerm(3, tsoil, state)  ! [SS-SWC S-2.7]
          endif

          ! Germination of arable crop growth
          if (flCropPrep .and. flCropSow) then
            call ArableLandGerm(4, tsoil, state)  ! [SS-SWC S-2.7]
          endif

        endif
      
      endif

! --- Initialization crop conditions ------------------------------------------
      
      if (flCropCalendar .and. .not. flCropHarvest) then            
        
        if (flCropReadFile) then
          
          ! fixed crop development
          if (croptype(state%crop%common%icrop) .eq. 1 .and. flCropEmergence) call CropFixed(1, state)

          ! detailed crop growth
          if (croptype(state%crop%common%icrop) .eq. 2 .and. flCropEmergence) call Wofost(1, state)

          ! detailed grass growth
          if (croptype(state%crop%common%icrop) .eq. 3) call Grass(1, tsoil, state)

          ! SS-CRP Phase 2 C-2.5: init JvL state directly (legacy globals retired).
          ! CropFixed/Wofost/Grass(1) set flhydrlift and twilt; hroot/hleaf/mfluxtable
          ! are written to state here since those subs lack state access.
          if (swdrought .eq. 2) then
            state%soilwater%hleaf = -2000.d0
            state%soilwater%hroot(1:state%mesh%numnod) = state%soilwater%h(1:state%mesh%numnod)  ! [SS-SWC S-2.7] [GR-BH Task 35]
            call MatricFlux(1, state%soilwater%h(1), 1, dummy_mf_, state)  ! [SS-SWC S-2.7]
          endif

          flCropReadFile = .false.
          state%crop%common%flCropReadFile = flCropReadFile   ! [SS-GR-CROPRT A5]

        endif

        ! update crop daynumber
        daycrop = daycrop + 1
        state%crop%common%daycrop = daycrop   ! [SS-GR-CROP A5.1]

        ! open crp-file
        if (swcrp.eq.1) call CropOutput(1, state)

        ! set correction of CO2 impact — FacCO2 now writes directly to state%crop%wofost%fco2*
        call FacCO2(state)  ! [SS-GR-CROPRT B6] FacCO2 handles state write; redundant dual-write removed
        
      endif

      ! set running average of minimum temperature (only for detailed crop growth)
      ! [SS-GR-FINAL B5] nofd/atmin7 read/written via state%atmosphere directly
      if (flCropEmergence .and. croptype(state%crop%common%icrop).ge.2) then
        state%atmosphere%nofd = min(state%atmosphere%nofd+1, 7)
        sumtmin = 0.0d0
        do i = state%atmosphere%nofd,2,-1
          state%atmosphere%atmin7(i) = state%atmosphere%atmin7(i-1)
          sumtmin = sumtmin + state%atmosphere%atmin7(i)
        end do
        i = 1
        state%atmosphere%atmin7(i) = tmn
        sumtmin = sumtmin + state%atmosphere%atmin7(i)
        tmnr = sumtmin / state%atmosphere%nofd
      else
        state%atmosphere%nofd = 0
      endif

      ! determine lowest compartment containing roots
      node = 1
      do while (state%mesh%zbotcp(node) .gt. (-state%crop%common%rd + 1.d-8))  ! [GR-BH C7] [GR-CROPWS B3]
        node = node + 1
      end do
      noddrz = node
      state%crop%common%noddrz = noddrz   ! [SS-GR-CROPRT A5]

      ! calculate potential and actual assimilation
      if (flCropEmergence .and. croptype(state%crop%common%icrop).ge.2) then
          
! check DAYNR during the day!!!!!!          
          
        ! phenological development rate 
        call astro (tc_daynr+1,lat,rad,dayl,daylp,sinld,cosld,difpp,atmtr,dsinbe)

        ! only for bulb crops (tulips etc..)
        if(state%crop%wofost%swbulb) then                                     ! [GR-CROPWS B3] swbulb → state%crop%wofost%swbulb
          ! remobilisation of carbohydrates from planted material
          if (state%crop%wofost%plwt.le.(0.0002d0*pld)) then                  ! [GR-CROPWS B3] plwt → state%crop%wofost%plwt
            ! no remobilisation at minimum weight motherbulb
            respmo = 0.0d0
            remo = 0.0d0
          else
            ! decrease weight mother organ starts at emergence.
            ! decrease consists of respiration and remobilisation
            decrmo = state%crop%wofost%plwt-(state%crop%wofost%plwt*(2.71828d0**remoc))  ! [GR-CROPWS B3]
            respmo = 0.025d0*(q10**((at_tavd-25.0d0)/10.0d0))*state%crop%wofost%plwt    ! [SS-GR-ATM B.5] [GR-CROPWS B3]
            if(respmo.lt.decrmo) then
              remo = decrmo - respmo
            else
              remo = 0.0d0
              respmo = decrmo
            end if
            ! weight motherbulb decreases by remobilisation and respiration
            plwt = state%crop%wofost%plwt - remo - respmo                     ! [GR-CROPWS B3] RHS plwt → state%crop%wofost%plwt
            state%crop%wofost%plwt = plwt   ! [SS-GR-CROP A5.1]
          endif
        endif

        ! daily gross assimilation
        effc = state%crop%wofost%fco2eff * eff  ! [SS-GR-CROPRT B6] fco2eff via state
        if (croptype(state%crop%common%icrop) .eq. 2) amax = state%crop%wofost%fco2amax * afgen (amaxtb,30,state%crop%common%dvs) * afgen (tmpftb,30,at_tavd)  ! [SS-GR-ATM B.5] [SS-GR-CROPRT B6]
        if (croptype(state%crop%common%icrop) .eq. 3) amax = state%crop%wofost%fco2amax * afgen (amaxtb,30,dble(daycrop)) * afgen (tmpftb,30,at_tavd)  ! [SS-GR-ATM B.5] [SS-GR-CROPRT B6]


        ! potential assimilation
        call totass (dayl,amax,effc,state%crop%common%laipot,state%crop%kdif,rad,difpp,dsinbe,sinld,cosld,dtgapot)  ! [GR-CROPWS B3] kdif → state%crop%kdif

        ! correction for low minimum temperature
        dtgapot = dtgapot * afgen (tmnftb,30,tmnr)

        ! potential assimilation in kg ch2o per ha
        pgasspot = dtgapot * 30.0d0/44.0d0

        ! only for bulb crops (tulips etc..)
        if(state%crop%wofost%swbulb) then                                     ! [GR-CROPWS B3]
          ! assimilation is raised with remobilisation from motherbulb
          ! using a factor of 1.11 given by De Ruijter et al.(1993)
          factblb  = 1.11d0
          pgasspot = pgasspot + remo*factblb   ! RHS pgasspot is local (just computed)
        endif

        ! reduction due to limited attainable maximum yield
        if (state%crop%grass%swpotrelmf.eq.2) pgasspot = pgasspot * state%crop%grass%relmf  ! [SS-GR-CROPRT B1]
        state%crop%wofost%pgasspot = pgasspot   ! [SS-GR-CROP A5.1]


        ! actual assimilation
        call totass (dayl,amax,effc,state%crop%lai,state%crop%kdif,rad,difpp,dsinbe,sinld,cosld,dtga)  ! [GR-CROPWS B3] kdif → state%crop%kdif; state%crop%lai=legacy

        ! correction for low minimum temperature
        dtga = dtga * afgen (tmnftb,30,tmnr)

        ! actual assimilation in kg ch2o per ha
        pgass = dtga * 30.0d0/44.0d0
        ! only for bulb crops (tulips etc..)
        if(state%crop%wofost%swbulb) then                                     ! [GR-CROPWS B3]
          ! assimilation is raised with remobilisation from motherbulb
          ! using a factor of 1.11 given by De Ruijter et al.(1993)
          factblb = 1.11d0
          pgass   = pgass + remo*factblb   ! RHS pgass is local (just computed)
        endif

        ! reduction due to limited attainable maximum yield
        pgass = pgass * state%crop%grass%relmf  ! [SS-GR-CROPRT B1]

        ! nitrogen stress reduction of pgass
        if (flCropNut) then
          call NUTRIE (NLUE,state%crop%wofost%wlv,state%crop%wofost%wst,state%crop%common%dvs,ANLV,ANST,NMXLV,NMAXLV,NMAXST,  &
     &      NMAXRT,LRNR,LSNR,NNI,RNFLV,RNFST,FRNX,FSTR)
          pgass = pgass * FSTR
        endif
        state%crop%wofost%pgass = pgass   ! [SS-GR-CROP A5.1]

      endif  

      return

      case (2)

! === calculation of potential crop rate and state variables =================

      if (flCropHarvest) return

! --- detailed crop growth -------------------------------------------------
      if (croptype(state%crop%common%icrop) .eq. 2) then
        if (flCropEmergence) then
          call Wofost(2, state)
        endif
      endif
! --- detailed grass growth  -----------------------------------------------
      if (croptype(state%crop%common%icrop) .eq. 3) then
        call Grass(2, tsoil, state)
      endif

      return

      case (3)

! === calculation of actual crop rate and state variables ==================

      if (flCropHarvest) return

! --- fixed crop development -----------------------------------------------
      if (croptype(state%crop%common%icrop).eq.1) then
        if(flCropEmergence) then
          call CropFixed(3, state)
        endif
      endif
! --- detailed crop growth -------------------------------------------------
      if (croptype(state%crop%common%icrop).eq.2) then
        if (flCropEmergence) then
          call Wofost(3, state)
         endif
      endif
! --- detailed grass growth  -----------------------------------------------
      if (croptype(state%crop%common%icrop).eq.3) then
        call Grass(3, tsoil, state)
      endif

      return

      case (4)

! === harvest of crop ======================================================
      
      if (flCropHarvest) return
     
      if (croptype(state%crop%common%icrop).le.2 .and. flCropEmergence)then

        ! Check flHarvestDay
        if (swharv.eq.0) then
          if (dabs(tc_t1900 - cropend(state%crop%common%icrop) - 1.d0) .lt. 1.0d-3) then  ! [GR-CROPWS B3] icrop → state%crop%common%icrop
            flHarvestDay = .true.
            state%crop%common%flHarvestDay = flHarvestDay   ! [SS-GR-CROP A5.1]
          endif
        else
          if (state%crop%common%dvs.ge.dvsend .or. dabs(tc_t1900 - cropend(state%crop%common%icrop) - 1.d0) .lt. 1.0d-3) then  ! [GR-CROPWS B3] dvs→state (task 4: synced after task 3); icrop→state
            flHarvestDay = .true.
            state%crop%common%flHarvestDay = flHarvestDay   ! [SS-GR-CROP A5.1]
          endif
        endif
        
        if (flCropEmergence .or. flHarvestDay) then
          if (croptype(state%crop%common%icrop).eq.1) then
            call CropFixed(4, state)
          endif
          if (croptype(state%crop%common%icrop).eq.2) then
            call Wofost(4, state)
          endif
        endif

      endif

! --- check timing of harvest
      
! --- fixed crop development -----------------------------------------------
      if (croptype(state%crop%common%icrop).eq.1)then
        if (flHarvestDay) then
          flCropEmergence = .false.
          state%crop%flCropEmergence = flCropEmergence   ! [SS-GR-ATM A5.2] dual-write
          flCropHarvest   = .true.
          state%crop%common%flCropHarvest = flCropHarvest   ! [SS-GR-CROPRT A5]
        endif
      endif

! --- detailed crop growth -------------------------------------------------
      if (croptype(state%crop%common%icrop).eq.2)then
        if (flHarvestDay) then
          flCropEmergence = .false.
          state%crop%flCropEmergence = flCropEmergence   ! [SS-GR-ATM A5.2] dual-write
          flCropHarvest   = .true.
          state%crop%common%flCropHarvest = flCropHarvest   ! [SS-GR-CROPRT A5]
        endif
      endif

! --- detailed grass growth ------------------------------------------------
      if (croptype(state%crop%common%icrop).eq.3)then
        if (dabs(tc_t1900 - cropend(state%crop%common%icrop) - 1.d0) .lt. 1.0d-3) then  ! [GR-CROPWS B3] icrop → state%crop%common%icrop
          flCropEmergence = .false.
          state%crop%flCropEmergence = flCropEmergence   ! [SS-GR-ATM A5.2] dual-write
          flCropHarvest   = .true.
          state%crop%common%flCropHarvest = flCropHarvest   ! [SS-GR-CROPRT A5]
        endif
      endif

      return
      
      case default
         call fatalerr_collected ('CropGrowth', 'Illegal value for TASK')
      end select

      end associate  ! tc_t1900, tc_daynr => state%timecontrol [TC-10]
      return

      end
! ----------------------------------------------------------------------
      subroutine totass (dayl,amax,eff,lai,kdif,avrad,difpp,            &
     &                   dsinbe,sinld,cosld,dtga)

!*  Purpose: This routine calculates the daily total gross CO2
!*           assimilation by performing a Gaussian integration over
!*           time. At three different times of the day, irradiance is
!*           computed and used to calculate the instantaneous canopy
!*           assimilation, whereafter integration takes place. More
!*           information on this routine is given by Spitters et al.
!*           (1988).

!*  FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)
!*  name   type meaning                                    units  class
!*  ----   ---- -------                                    -----  -----
!*  DAYL    R8  Astronomical daylength (base = 0 degrees)     h      O
!*  AMAX    R8  Assimilation rate at light saturation      kg CO2/   I
!*                                                        ha leaf/h   
!*  EFF     R8  Initial light use efficiency              kg CO2/J/  I
!*                                                        ha/h m2 s   
!*  LAI     R8  Leaf area index                             ha/ha    I
!*  KDIF    R8  Extinction coefficient for diffuse light             I
!*  AVRAD   R8  Daily shortwave radiation                  J m-2 d-1 I
!*  DIFPP   R8  Diffuse irradiation perpendicular to direction of
!*              light                                      J m-2 s-1 I
!*  DSINBE  R8  Daily total of effective solar height         s      I
!*  SINLD   R8  Seasonal offset of sine of solar height       -      I
!*  COSLD   R8  Amplitude of sine of solar height             -      I
!*  DTGA    R8  Daily total gross assimilation           kg CO2/ha/d O

!*  FATAL ERROR CHECKS: none
!*  SUBROUTINES and FUNCTIONS called : ASSIM
!*  FILE usage : none

!*  Authors: Daniel van Kraalingen 
!*  Date   : April 1991

!*  Modification: Implementated in Swap3.2.41: 
!*                several small adjustments (R4->R8, small caps)
!*  Author      : Joop Kroes
!*  Date        : May 2014

      implicit none

!*     formal parameters
      real(8) dayl,amax,eff,lai,kdif,avrad,difpp,dsinbe,sinld,cosld,dtga

!*     local parameters
      integer i1
      real(8) hour,pi,sinb,par,pardif,pardir,fgros
      real(8) xgauss(3),wgauss(3)

      parameter (pi=3.1415926d0)
      save

!**
!*     gauss points and weights are stored in an array
      data xgauss /0.1127017d0, 0.5000000d0, 0.8872983d0/
      data wgauss /0.2777778d0, 0.4444444d0, 0.2777778d0/

!*     calculation of assimilation is done only when it will not be zero
!*     (AMAX >0, LAI >0)
      dtga  = 0.0d0
      if (amax.gt.0.0d0.and.lai.gt.0.0d0) then
         do 10 i1=1,3
            hour   = 12.0d0+0.5d0*dayl*xgauss(i1)
            sinb   = max(0.0d0,                                         &
     &                  sinld+cosld*dcos(2.0d0*pi*(hour+12.0d0)/24.0d0))
            par    = 0.5d0*avrad*sinb*(1.0d0+0.4d0*sinb)/dsinbe
            pardif = min(par,sinb*difpp)
            pardir = par-pardif
            call assim (amax,eff,lai,kdif,sinb,pardir,pardif,fgros)
            dtga = dtga+fgros*wgauss(i1)
10       continue
         dtga = dtga*dayl
      end if

      return
      end

      subroutine assim (amax,eff,lai,kdif,sinb,pardir,pardif,fgros)

!*     Chapter 13 in documentation WOFOST Version 4.1 (1988)

!*     This routine calculates the gross CO2 assimilation rate of
!*     the whole crop, FGROS, by performing a Gaussian integration
!*     over depth in the crop canopy. At three different depths in
!*     the canopy, i.e. for different values of LAI, the
!*     assimilation rate is computed for given fluxes of photosynthe-
!*     tically active radiation, whereafter integration over depth
!*     takes place. More information on this routine is given by
!*     Spitters et al. (1988). The input variables SINB, PARDIR
!*     and PARDIF are calculated in routine TOTASS.

!*     Subroutines and functions called: none.
!*     Called by routine TOTASS.

!*     Author: D.W.G. van Kraalingen, 1986

!*  Modification: Implementated in Swap3.2.41: 
!*                several small adjustments (R4->R8, small caps)
!*  Author      : Joop Kroes
!*  Date        : May 2014

!*  FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)
!*  name   type meaning                                    units  class
!*  ----   ---- -------                                    -----  -----
!*  AMAX    R8  Maximum CO2 assimilation rate              kg/ha/hr  I
!*  EFF     R8  Light use efficiency of a leaf         kg CO2 / J adsorbed  I
!*  LAI     R8  Leaf area index                               -      I
!*  KDIF    R8  Extinction coefficient for diffuse visible light -   I
!*  SINB    R8  ...........nog invullen ..............               -      I
!*  PARDIR  R8  ...........nog invullen ..............                 -      I
!*  DIFPP   R8  Diffuse irradiation perpendicular to direction of     
!*              light                                      J m-2 s-1 I
!*  PARDIF  R8  ...........nog invullen ..............             -      I
!*  FGROS   R8  ...........nog invullen ..............           s      O
!**    
!*13.1 declarations
      implicit none

!*     formal parameters
      real(8) amax,eff,lai,kdif,sinb,pardir,pardif,fgros

!*     local parameters
      integer i
      real(8) scv,refh,refs,kdirbl,kdirt,laic,visdf,vist,visd,visshd
      real(8) fgrsh,vispp,fgrsun,fslla,fgl
      real(8) xgauss(3),wgauss(3)

      save

!*     initialize GAUSS array and scattering coefficient
      data xgauss /0.1127017d0, 0.5000000d0, 0.8872983d0/
      data wgauss /0.2777778d0, 0.4444444d0, 0.2777778d0/
      data scv /0.2d0/

!*13.2 extinction coefficients KDIF,KDIRBL,KDIRT
      refh   = (1.0d0-dsqrt(1.0d0-scv))/(1.0d0+dsqrt(1.0d0-scv))
      refs   = refh*2.0d0/(1.0d0+1.6d0*sinb)
      kdirbl = (0.5d0/sinb)*kdif/(0.8d0*dsqrt(1.0d0-scv))
      kdirt  = kdirbl*dsqrt(1.0d0-scv)

!*13.3 three-point Gaussian integration over LAI
      fgros  = 0.0d0
      do 10 i=1,3
         laic   = lai*xgauss(i)
!*        absorbed diffuse radiation (VISDF),light from direct
!*        origine (VIST) and direct light(VISD)
         visdf  = (1.0d0-refs)*pardif*kdif  *exp (-kdif  *laic)
         vist   = (1.0d0-refs)*pardir*kdirt *exp (-kdirt *laic)
         visd   = (1.0d0-scv) *pardir*kdirbl*exp (-kdirbl*laic)
!*        absorbed flux in W/m2 for shaded leaves and assimilation
         visshd = visdf+vist-visd
         fgrsh  = amax*(1.0d0-exp(-visshd*eff/max(2.0d0,amax)))
!*        direct light absorbed by leaves perpendicular on direct
!*        beam and assimilation of sunlit leaf area
         vispp  = (1.0d0-scv)*pardir/sinb
         if (vispp.le.0.0d0) then
            fgrsun = fgrsh
         else
            fgrsun = amax*(1.0d0-(amax-fgrsh)                           &
     &          *(1.0d0-exp (-vispp*eff/max(2.0d0,amax)))/ (eff*vispp))
         end if
!*        fraction of sunlit leaf area (FSLLA) and local
!*        assimilation rate (FGL)
         fslla  = exp (-kdirbl*laic)
         fgl    = fslla*fgrsun+(1.0d0-fslla)*fgrsh
!*        integration
         fgros  = fgros+fgl*wgauss(i)
10    continue

      fgros  = fgros*lai
      return
      end

      subroutine astro (iday,lat,avrad,                                 &
     &                  dayl,daylp,sinld,cosld,difpp,atmtr,dsinbe)

!*  Purpose: This subroutine calculates astronomic daylength,
!*           diurnal radiation characteristics such as the atmospheric
!*           transmission, diffuse radiation etc.. This routine has
!*           been modified so that it uses arrays to hold some input
!*           output variables for faster processing 

!*  FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)
!*  name   type meaning                                    units  class
!*  ----   ---- -------                                    -----  -----
!*  IDAY    I4  Day number (Jan 1st = 1)                      -      I
!*  LAT     R8  Latitude of the site                       degrees   I
!*  AVRAD   R8  Daily shortwave radiation                  J m-2 d-1 I
!*  DAYL    R8  Astronomical daylength (base = 0 degrees)     h      O
!*  DAYLP   R8  Astronomical daylength (base =-4 degrees)     h      O
!*  SINLD   R8  Seasonal offset of sine of solar height       -      O
!*  COSLD   R8  Amplitude of sine of solar height             -      O
!*  DIFPP   R8  Diffuse irradiation perpendicular to direction of     
!*              light                                      J m-2 s-1 O
!*  ATMTR   R8  Daily atmospheric transmission                -      0
!*  DSINBE  R8  Daily total of effective solar height         s      O

!*  FATAL ERROR CHECKS: none
!*  SUBROUTINES and FUNCTIONS called : none
!*  FILE usage : none

!*  Authors: Daniel van Kraalingen
!*  Date   : April 1991

!*  Modification: Include checks for 0<=daylength<=24 hour
!*                Remove caching of results
!*  Author      : Allard de Wit
!*  Date        : January 2011

!*  Modification: Implementated in Swap3.2.41: 
!*                several small adjustments (R4->R8, small caps)
!*  Author      : Joop Kroes
!*  Date        : May 2014

      use error_mod, only: fatalerr_collected
      implicit none
!*     formal parameters
      integer iday
      real(8) lat,avrad,dayl,daylp,sinld,cosld,difpp,atmtr,dsinbe

!*     local parameters
      real(8) pi,angle,rad
      real(8) dec,sc,aob,aob_corr,angot,dsinb,frdif

      parameter (pi=3.1415926d0, angle=-4.0d0, rad=0.0174533d0)

!*     Error check on latitude
      if (dabs(lat).gt.90.d0) call fatalerr_collected                   &
     &   ('astro','lat > 90 or lat < -90')

!*     Declination and solar constant for this day
      dec = -asin(dsin(23.45d0*rad)*dcos(2.d0*pi*dble(iday+10)/365.0d0))
      sc  = 1370.d0*(1.d0+0.033d0*dcos(2.d0*pi*dble(iday)/365.d0))

!*     calculation of daylength from intermediate variables
!*     SINLD, COSLD and AOB
      sinld = dsin(rad*lat)*dsin(dec)
      cosld = dcos(rad*lat)*dcos(dec)
      aob = sinld/cosld

!*     For very high latitudes and days in summer and winter a limit is  
!*     inserted to avoid math errors when daylength reaches 24 hours in 
!*     summer or 0 hours in winter.

!*     Calculate solution for base=0 degrees
      if (dabs(aob).le.1.0d0) then
         dayl  = 12.0d0*(1.d0+2.d0*asin(aob)/pi)
!*        integrals of sine of solar height
         dsinb  = 3600.d0*                                              &
     &            (dayl*sinld+24.d0*cosld*dsqrt(1.d0-aob**2)/pi)
         dsinbe = 3600.d0*                                              &
     &            (dayl*(sinld+0.4d0*(sinld**2+cosld**2*0.5d0))+  &
     &     12.d0*cosld*(2.d0+3.d0*0.4d0*sinld)*dsqrt(1.d0-aob**2)/pi)
      else
         if (aob.gt.1.0d0)  dayl = 24.0d0
         if (aob.lt.-1.0d0) dayl =  0.0d0
!*        integrals of sine of solar height      
         dsinb  = 3600.d0*(dayl*sinld)
         dsinbe = 3600.d0*                                              &
     &            (dayl*(sinld+0.4d0*(sinld**2+cosld**2*0.5d0)))
      endif

!*     Calculate solution for base=-4 (ANGLE) degrees
      aob_corr = (-dsin(angle*rad)+sinld)/cosld
      if (dabs(aob_corr).le.1.0d0) then 
         daylp = 12.0d0*(1.d0+2.d0*asin(aob_corr)/pi)
      else
         if (aob_corr.gt.1.0d0)  daylp = 24.0d0
         if (aob_corr.lt.-1.0d0) daylp =  0.0d0
      endif

!*     extraterrestrial radiation and atmospheric transmission
      angot  = sc*dsinb
!*     Check for DAYL=0 as in that case the angot radiation is 0 as well
      if (dayl.gt.0.0d0) then
          atmtr = avrad/angot
      else
          atmtr = 0.0d0
      endif

!*     estimate fraction diffuse irradiation
      if (atmtr.gt.0.75d0) frdif = 0.23d0
      if (atmtr.le.0.75d0.and.atmtr.gt.0.35d0)                          &
     &  frdif = 1.33d0-1.46d0*atmtr
      if (atmtr.le.0.35d0.and.atmtr.gt.0.07d0)                          &
     &  frdif = 1.d0-2.3d0*(atmtr-0.07d0)**2
      if (atmtr.le.0.07d0) frdif = 1.d0

      difpp = frdif*atmtr*0.5d0*sc

      RETURN
      END

