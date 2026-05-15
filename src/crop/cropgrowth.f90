! File VersionID:
!   $Id: cropgrowth.f90 380 2018-05-28 14:01:08Z heine003 $
!
! Phase 4e Task A2 audit: all `call fatalerr` sites in this file are in
! surviving physics dispatchers (CropGrowth, cropfixed, cropoutput,
! ArableLandGerm, FacCO2, wofost, grass, astro, outbalcrop*, chckcbl).
! The crop sub-readers (readwofost, readcropfixed, readgrass) live in
! src/io/readswap.f90 — none of them are here. Phase 4e Task A4/A5
! replaces all of this file's calls via fatalerr_collected (singleton).
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
      !   swinco, croptype, swcrp, swdrought, swharv, swbulb: config switches, no state home
      !   daycrop, rd, rdpot, lai, laipot, cf, ch, tsum, dvs: runtime crop state (dual-write
      !     to state%crop%common%; global still needed pending Phase C global retirement)
      !   cwdmpot, cwdm, wsopot, wso, wlvpot, wlv, wstpot, wst, wrtpot, wrt: WOFOST pools
      !   tmn, lat, rad: meteo scalars, no state%atmosphere scalar home
      !   eff, amaxtb, tmpftb, tmnftb, kdif: physiology params/tables, no state home
      !   plwt, remoc, pld, q10, pgasspot, pgass: physiology params
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
      use variables, only: &                                                 ! [SS-GR-CROPRT B1/B6] DEFERRED
        icrop, flCropCalendar, cropstart, cropend, flCropEmergence,         &
        flCropHarvest, flCropReadFile, flCropPrep, flCropSow, flCropGerm,   &
        swinco, croptype, daycrop, rd, rdpot, lai, laipot, cf, ch, tsum,   &
        dvs, cwdmpot, cwdm, wsopot, wso, wlvpot, wlv, wstpot, wst,        &
        wrtpot, wrt, tmn, lat, rad,                                         &
        albedo, rsc, cumdens,                                               &
        eff, amaxtb, tmpftb, tmnftb, kdif, swdrought, swcrp, dvsend,       &
        swharv, swbulb, plwt, remoc, pld, q10, pgasspot, pgass,            &
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
      implicit none

      ! Explicit interfaces for non-module subs that now take tsoil(:)
      ! ArableLandGerm removed — now a module procedure in cropgrowth_helpers_mod.
      ! wofost removed — now a module procedure in cropwofost_runtime_mod.
      interface
         subroutine grass(task, tsoil, state)
            use swap_state_mod, only: swap_state_t
            integer :: task
            real(8), intent(in) :: tsoil(:)
            type(swap_state_t), intent(inout), optional :: state
         end subroutine grass
         subroutine cropfixed(task, state)
            use swap_state_mod, only: swap_state_t
            integer :: task
            type(swap_state_t), intent(inout), optional :: state
         end subroutine cropfixed
      end interface

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
        call nocrop ()
        ! [SS-GR-CROP A5.1] nocrop() has no state arg; mirror all zeroed fields here
        state%crop%lai               = lai       ! GR-ATM fix: ProcessMeteoDay reads this
        state%crop%common%rd         = rd
        state%crop%common%rdpot      = rdpot
        state%crop%common%laipot     = laipot
        state%crop%common%cf         = cf
        state%crop%common%ch         = ch
        state%crop%common%tsum       = tsum
        state%crop%common%dvs        = dvs
        state%crop%common%albedo     = albedo   ! [SS-GR-CROPRT A5]
        state%crop%common%rsc        = rsc      ! [SS-GR-CROPRT A5]
        state%crop%wofost%cwdmpot    = cwdmpot
        state%crop%wofost%cwdm       = cwdm
        state%crop%wofost%wsopot     = wsopot
        state%crop%wofost%wso        = wso
        state%crop%wofost%wlvpot     = wlvpot
        state%crop%wofost%wlv        = wlv
        state%crop%wofost%wstpot     = wstpot
        state%crop%wofost%wst        = wst
        state%crop%wofost%wrtpot     = wrtpot
        state%crop%wofost%wrt        = wrt
      endif

! --- check crop emergence ----------------------------------------------------

      ! reset if new crop
      if (flCropCalendar) then
        if (dabs(tc_t1900 - cropstart(icrop)) .lt. tiny) then
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
          if (croptype(icrop) .le. 2) then
            flCropEmergence = .false.
            state%crop%flCropEmergence = flCropEmergence   ! [SS-GR-ATM A5.2] dual-write
          endif
        endif
      endif

! --- Preparation, Sowing and Germination of arable crop growth ---------------
      if (flCropCalendar .and. .not. flCropHarvest .and. croptype(icrop) .le. 2) then

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
                     if (icrop >= 1 .and. icrop <= size(crop_config_global%rotation_loaded)) then
                        if (crop_config_global%rotation_loaded(icrop)) then
                           rot_type = crop_config_global%rotation_type(icrop)
                           select case (rot_type)
                           case (1)
                              ! type=1 cropfixed
                              if (allocated(crop_config_global%rotation_fixed)) then
                                 use_cache    = .true.
                                 swprep_cache = crop_config_global%rotation_fixed(icrop)%swprep
                                 swsow_cache  = crop_config_global%rotation_fixed(icrop)%swsow
                                 swgerm_cache = crop_config_global%rotation_fixed(icrop)%swgerm
                              end if
                           case (2)
                              ! type=2 wofost: cache-hit when swprep=0 AND swsow=0.
                              ! swgerm=0/1/2 are all handled in the cache body below.
                              if (allocated(crop_config_global%rotation_wofost)) then
                                 swprep_cache = crop_config_global%rotation_wofost(icrop)%preparation%swprep
                                 swsow_cache  = crop_config_global%rotation_wofost(icrop)%sowing%swsow
                                 swgerm_cache = crop_config_global%rotation_wofost(icrop)%germination%swgerm
                                 if (swprep_cache == 0 .and. swsow_cache == 0) then
                                    use_cache = .true.
                                    gp => crop_config_global%rotation_wofost(icrop)%germination
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
          if (croptype(icrop) .eq. 1 .and. flCropEmergence) call CropFixed(1, state)

          ! detailed crop growth
          if (croptype(icrop) .eq. 2 .and. flCropEmergence) call Wofost(1, state)

          ! detailed grass growth
          if (croptype(icrop) .eq. 3) call Grass(1, tsoil, state)

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
      if (flCropEmergence .and. croptype(icrop).ge.2) then
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
      do while (state%mesh%zbotcp(node) .gt. (-rd + 1.d-8))  ! [GR-BH C7]
        node = node + 1
      end do
      noddrz = node
      state%crop%common%noddrz = noddrz   ! [SS-GR-CROPRT A5]

      ! calculate potential and actual assimilation
      if (flCropEmergence .and. croptype(icrop).ge.2) then
          
! check DAYNR during the day!!!!!!          
          
        ! phenological development rate 
        call astro (tc_daynr+1,lat,rad,dayl,daylp,sinld,cosld,difpp,atmtr,dsinbe)

        ! only for bulb crops (tulips etc..)
        if(swbulb.eq.1) then
          ! remobilisation of carbohydrates from planted material
          if (plwt.le.(0.0002d0*pld)) then
            ! no remobilisation at minimum weight motherbulb
            respmo = 0.0d0
            remo = 0.0d0
          else
            ! decrease weight mother organ starts at emergence.
            ! decrease consists of respiration and remobilisation
            decrmo = plwt-(plwt*(2.71828d0**remoc))
            respmo = 0.025d0*(q10**((at_tavd-25.0d0)/10.0d0))*plwt  ! [SS-GR-ATM B.5] tavd→state%atmosphere%tavd
            if(respmo.lt.decrmo) then
              remo = decrmo - respmo
            else
              remo = 0.0d0
              respmo = decrmo
            end if
            ! weight motherbulb decreases by remobilisation and respiration
            plwt = plwt - remo - respmo
            state%crop%wofost%plwt = plwt   ! [SS-GR-CROP A5.1]
          endif
        endif

        ! daily gross assimilation
        effc = state%crop%wofost%fco2eff * eff  ! [SS-GR-CROPRT B6] fco2eff via state
        if (croptype(icrop) .eq. 2) amax = state%crop%wofost%fco2amax * afgen (amaxtb,30,dvs) * afgen (tmpftb,30,at_tavd)  ! [SS-GR-ATM B.5] [SS-GR-CROPRT B6]
        if (croptype(icrop) .eq. 3) amax = state%crop%wofost%fco2amax * afgen (amaxtb,30,dble(daycrop)) * afgen (tmpftb,30,at_tavd)  ! [SS-GR-ATM B.5] [SS-GR-CROPRT B6]


        ! potential assimilation
        call totass (dayl,amax,effc,laipot,kdif,rad,difpp,dsinbe,sinld,cosld,dtgapot)

        ! correction for low minimum temperature
        dtgapot = dtgapot * afgen (tmnftb,30,tmnr)

        ! potential assimilation in kg ch2o per ha
        pgasspot = dtgapot * 30.0d0/44.0d0

        ! only for bulb crops (tulips etc..)
        if(swbulb.eq.1) then
          ! assimilation is raised with remobilisation from motherbulb
          ! using a factor of 1.11 given by De Ruijter et al.(1993)
          factblb  = 1.11d0
          pgasspot = pgasspot + remo*factblb
        endif

        ! reduction due to limited attainable maximum yield
        if (state%crop%grass%swpotrelmf.eq.2) pgasspot = pgasspot * state%crop%grass%relmf  ! [SS-GR-CROPRT B1]
        state%crop%wofost%pgasspot = pgasspot   ! [SS-GR-CROP A5.1]


        ! actual assimilation
        call totass (dayl,amax,effc,lai,kdif,rad,difpp,dsinbe,sinld,cosld,dtga)

        ! correction for low minimum temperature
        dtga = dtga * afgen (tmnftb,30,tmnr)

        ! actual assimilation in kg ch2o per ha
        pgass = dtga * 30.0d0/44.0d0
        ! only for bulb crops (tulips etc..)
        if(swbulb.eq.1) then
          ! assimilation is raised with remobilisation from motherbulb
          ! using a factor of 1.11 given by De Ruijter et al.(1993)
          factblb = 1.11d0
          pgass   = pgass + remo*factblb
        endif

        ! reduction due to limited attainable maximum yield
        pgass = pgass * state%crop%grass%relmf  ! [SS-GR-CROPRT B1]

        ! nitrogen stress reduction of pgass
        if (flCropNut) then
          call NUTRIE (NLUE,WLV,WST,DVS,ANLV,ANST,NMXLV,NMAXLV,NMAXST,  &
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
      if (croptype(icrop) .eq. 2) then
        if (flCropEmergence) then
          call Wofost(2, state)
        endif
      endif
! --- detailed grass growth  -----------------------------------------------
      if (croptype(icrop) .eq. 3) then
        call Grass(2, tsoil, state)
      endif

      return

      case (3)

! === calculation of actual crop rate and state variables ==================

      if (flCropHarvest) return

! --- fixed crop development -----------------------------------------------
      if (croptype(icrop).eq.1) then
        if(flCropEmergence) then
          call CropFixed(3, state)
        endif
      endif
! --- detailed crop growth -------------------------------------------------
      if (croptype(icrop).eq.2) then
        if (flCropEmergence) then
          call Wofost(3, state)
         endif
      endif
! --- detailed grass growth  -----------------------------------------------
      if (croptype(icrop).eq.3) then
        call Grass(3, tsoil, state)
      endif

      return

      case (4)

! === harvest of crop ======================================================
      
      if (flCropHarvest) return
     
      if (croptype(icrop).le.2 .and. flCropEmergence)then

        ! Check flHarvestDay
        if (swharv.eq.0) then
          if (dabs(tc_t1900 - cropend(icrop) - 1.d0) .lt. 1.0d-3) then
            flHarvestDay = .true.
            state%crop%common%flHarvestDay = flHarvestDay   ! [SS-GR-CROP A5.1]
          endif
        else
          if (dvs.ge.dvsend .or. dabs(tc_t1900 - cropend(icrop) - 1.d0) .lt. 1.0d-3) then
            flHarvestDay = .true.
            state%crop%common%flHarvestDay = flHarvestDay   ! [SS-GR-CROP A5.1]
          endif
        endif
        
        if (flCropEmergence .or. flHarvestDay) then
          if (croptype(icrop).eq.1) then
            call CropFixed(4, state)
          endif
          if (croptype(icrop).eq.2) then
            call Wofost(4, state)
          endif
        endif

      endif

! --- check timing of harvest
      
! --- fixed crop development -----------------------------------------------
      if (croptype(icrop).eq.1)then
        if (flHarvestDay) then
          flCropEmergence = .false.
          state%crop%flCropEmergence = flCropEmergence   ! [SS-GR-ATM A5.2] dual-write
          flCropHarvest   = .true.
          state%crop%common%flCropHarvest = flCropHarvest   ! [SS-GR-CROPRT A5]
        endif
      endif

! --- detailed crop growth -------------------------------------------------
      if (croptype(icrop).eq.2)then
        if (flHarvestDay) then
          flCropEmergence = .false.
          state%crop%flCropEmergence = flCropEmergence   ! [SS-GR-ATM A5.2] dual-write
          flCropHarvest   = .true.
          state%crop%common%flCropHarvest = flCropHarvest   ! [SS-GR-CROPRT A5]
        endif
      endif

! --- detailed grass growth ------------------------------------------------
      if (croptype(icrop).eq.3)then
        if (dabs(tc_t1900 - cropend(icrop) - 1.d0) .lt. 1.0d-3) then
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
!
! ----------------------------------------------------------------------
      subroutine cropfixed (task, state)
! ----------------------------------------------------------------------
!     date               : august 2004
!     purpose            : simple crop growth routine for swap
! SS-CRP C-2.5: state added (optional, intent in) to read flWrtNonox.
! SS-TC TC-10: t1900 read via state%timecontrol tc_t1900 alias.
! SS-GR-ATM A5.1: intent changed inout to allow dual-write in cropfixed_init_from_config.
! [GR-CROP Phase B/6] narrow use variables
! [SS-GR-CROPRT B2] DEFERRED — cropfixed: all remaining variables globals:
!   magrs: array dim (could → swap_array_dimensions, deferred with rest)
!   icrop, dvs, idev, lai, tsum, cf, ch, rd, rdpot: runtime state (dual-write to state%crop%;
!     global still needed pending Phase C global retirement)
!   max_resp_factor: config field (config%crop%fixed%), needs config threading
!   swrd, rdi, rri, rdc, swgc, swcf, swinter, swdrought, swdmi2rd: switches, no state home
!   cropstart, tbase, tsumea, tsumam, rdmax, rdm: config params, no state home
!   siccapact, siccaplai: state%atmosphere%siccapact migrated; siccaplai no home
!   w_root_ss, wiltpoint, twilt, flhydrlift: JvL params, no state home
!   cftb/chtb/cfeictb: state%crop%fixed homes exist (A5.2 dual-write) but migration to
!     state reads deferred — cropfixed takes optional state, substitution needs
!     non-optional refactor or present() guards; deferred to Phase C cleanup
!   gc, cfeic, gctb, rdtb, mrftb, wrtb: fixed-crop tables/scalars, no state home
!   swinco, reltr: switches, no state home
! ----------------------------------------------------------------------
      use variables, only: magrs, icrop, dvs, idev, lai, tsum, cf, ch, &  ! [SS-GR-CROPRT B2] DEFERRED
                           rd, rdpot, max_resp_factor, swrd, rdi, rri,  &
                           rdc, swgc, swcf, swinter, swdrought, swdmi2rd, &
                           cropstart, tbase, tsumea, tsumam, rdmax, rdm, &
                           siccapact, siccaplai, w_root_ss, wiltpoint,   &
                           twilt, flhydrlift, gc, cfeic,                 &
                           gctb, cftb, chtb, cfeictb, rdtb, mrftb, wrtb, &
                           swinco, reltr
      use soilhydraulics_utils, only: watcon
      use array_utils, only: afgen
      use rootextraction_mod, only: MatricFlux
      use swap_constants, only: tiny, nihil
      use error_mod, only: fatalerr_collected
      use swap_state_mod, only: swap_state_t
      implicit none

      type(swap_state_t), intent(inout), optional :: state

! --- local variables
      integer   i,task,lcc,swhydrlift
      real(8)   dummy,dtsum,dvr

! --- rooting
      real(8)   rrpot,rr

      save
! ----------------------------------------------------------------------
      ! TC-10: t1900 read via state%timecontrol tc_t1900 alias.
      ! [SS-BMI2 Task 4] tstart added to associate
      ! [SS-GR-ATM B.5] at_tav alias for tav read migration
      associate( tc_t1900 => state%timecontrol%t1900, &  ! TC-10
                 tstart   => state%timecontrol%tstart, &  ! [SS-BMI2 Task 4]
                 at_tav   => state%atmosphere%Tav      )  ! [SS-GR-ATM B.5]

      select case (task)
      case (1)

! === initialization ===================================================
      
! --- read crop data: dispatch on per-rotation typed-config cache
!     (ADR 0016). Falls back to legacy reader for rotations whose
!     .crp.toml is not yet authored or for rotation types not yet
!     ported (Phase 1: only type=1 cropfixed; Phases 2/3 add types
!     2 and 3). Teardown: end of Phase 4 removes the else-branch.
      block
         use crop_config_global_mod, only: crop_config_global
         use cropfixed_init_mod, only: cropfixed_init_from_config
         logical :: use_cache
         use_cache = .false.
         if (associated(crop_config_global)) then
            if (allocated(crop_config_global%rotation_loaded)) then
               if (icrop >= 1 .and. icrop <= size(crop_config_global%rotation_loaded)) then
                  if (crop_config_global%rotation_loaded(icrop)) use_cache = .true.
               end if
            end if
         end if
         if (use_cache) then
            call cropfixed_init_from_config(crop_config_global%rotation_fixed(icrop), icrop, lcc, state)
            ! [SS-GR-ATM A5.1] state passed for dual-write of kdif/kdir/swcf/cofab
            ! swhydrlift is read by legacy readcropfixed only inside the
            ! swdrought=2 branch (stub-errored in Phase 1). Set to 0 here
            ! to mirror the default; Phase 2 (cropwofost) reuses this
            ! same field on its own dispatch path.
            swhydrlift = 0
         else
            ! ADR 0016 cache-miss: typed config required for type=1 rotations.
            ! No silent legacy fallback — the user must author cropfixed.crp.toml.
            call fatalerr_collected('cropgrowth/CropFixed', &
               'cropfixed rotation has no loaded .crp.toml — author the file or use the legacy executable.')
         end if
      end block

! --- maximum rooting depth
      if (swrd.eq.1) then
        rdm = rdmax
      else
        rdm = min(rdmax,rdc)
      endif
      if (present(state)) state%crop%common%rdm = rdm   ! [SS-GR-CROP A5.1]

! --- skip next initialization if crop parameters are read from *.END file
      if (tc_t1900 - tstart .gt. tiny .or. swinco .ne. 3 .or.           &
     &  dabs(tc_t1900 - cropstart(icrop)) .lt. tiny) then

        dvs = 0.0d0

! --- actual rooting depth
        if (swrd.eq.1) then
          rd = afgen (rdtb,22,dvs)
          rd = min(rd,rdm)
        else
          rd = min(rdi,rdm)
        endif
        rdpot = rd
        if (present(state)) then
          state%crop%common%dvs   = dvs    ! [SS-GR-CROP A5.1]
          state%crop%common%rd    = rd     ! [SS-GR-CROP A5.1]
          state%crop%common%rdpot = rdpot  ! [SS-GR-CROP A5.1]
        endif

      endif

! --- initial lai or sc
      lai = afgen (gctb,(2*magrs),dvs)
      if (swgc.eq.2) then
        gc = lai
        lai = lai*3.0d0
      endif
      if (present(state)) state%crop%lai = lai   ! [SS-GR-ATM A5.2] dual-write

! --- initial crop factor or crop height
      cf = afgen (cftb,(2*magrs),dvs)
      ch = afgen (chtb,(2*magrs),dvs)
      if (swcf.eq.3) then
        cfeic = afgen (cfeictb,(2*magrs),dvs)
      endif
      if (present(state)) then
        state%crop%common%cf = cf   ! [SS-GR-CROP A5.1]
        state%crop%common%ch = ch   ! [SS-GR-CROP A5.1]
        if (swcf.eq.3) state%crop%fixed%cfeic = cfeic   ! [SS-GR-CROP A5.1]
      endif

! --- initial storage on canopy
      if (swinter.eq.3) then
        siccapact = siccaplai*lai
        if (present(state)) state%atmosphere%siccapact = siccapact   ! [SS-GR-ATM A5.2] dual-write
      endif

! --- initial dry weight of roots at soil surface; oxygen module
      W_root_ss = afgen (wrtb,(2*magrs),dvs)

! --- initial ratio root total respiration / maintenance respiration; oxygen module
      max_resp_factor = afgen (mrftb,(2*magrs),dvs)

! --- initialize matric flux potential (SS-CRP C-2.5: hroot/hleaf/mfluxtable
!     init moved to CropGrowth dispatcher which has access to state).
      if (swdrought .eq. 2) then
        if (swhydrlift .eq. 1) then
          flhydrlift = .true.
        else
          flhydrlift = .false.
        endif
        do i = 1,state%mesh%numnod  ! [GR-BH C7]
         twilt(i) = watcon(wiltpoint, &
                            state%soilwater%vg_params(i), &
                            state%soilwater%iHWCKmodel(state%soilwater%layer(i)), &
                            i, state%soilwater)                    ! [SS-GR-UTILS Task 5]
        enddo
      endif

      return

      case (2)
      continue

! === calculate potential rate and state variables ======================
      case (3)

! === calculate actual rate and state variables ======================

! --- increase in temperature sum
      dtsum = max (0.0d0,at_tav-tbase)  ! [SS-GR-ATM B.5] tav→state%atmosphere%Tav

! --- development rate
      if (idev.eq.1) then
        dvr = 2.0/lcc
      elseif (idev.eq.2) then
        if (dvs.lt.1.0d0) then
          dvr = dtsum/tsumea
        else
          dvr = dtsum/tsumam
        endif
      endif

! --- water stress
      ! SS-ATM Phase 2 Task A-2.3: ptra read from state%atmosphere (atmosphere home).
      if(dabs(state%atmosphere%ptra).lt.nihil) then
        reltr = 1.0d0
      else
        reltr = max(min(state%soilwater%tra/state%atmosphere%ptra,1.0d0),0.0d0)  ! [SS-SWC S-2.7]
      endif

! ----integrals of the crop --------------------------------------------

! --- phenological development stage
      dvs = min(dvs+dvr,2.d0)
      tsum = tsum + dtsum
      if (present(state)) then
        state%crop%common%dvs  = dvs    ! [SS-GR-CROP A5.1]
        state%crop%common%tsum = tsum   ! [SS-GR-CROP A5.1]
      endif

! --- leaf area index or soil cover fraction
      lai = afgen (gctb,(2*magrs),dvs)
      if (swgc.eq.2) then
        gc = lai
        lai = lai*3.0d0
      endif
      if (present(state)) state%crop%lai = lai   ! [SS-GR-ATM A5.2] dual-write

! --- crop factor or crop height
      cf        = afgen (cftb,(2*magrs),dvs)
      ch        = afgen (chtb,(2*magrs),dvs)
      if (swcf.eq.3) then
        cfeic     = afgen (cfeictb,(2*magrs),dvs)
      endif
      if (present(state)) then
        state%crop%common%cf = cf   ! [SS-GR-CROP A5.1]
        state%crop%common%ch = ch   ! [SS-GR-CROP A5.1]
        if (swcf.eq.3) state%crop%fixed%cfeic = cfeic   ! [SS-GR-CROP A5.1]
      endif

! --- update canopy storage capacity
      if (swinter.eq.3) then
        siccapact = siccaplai*lai
        if (present(state)) state%atmosphere%siccapact = siccapact   ! [SS-GR-ATM A5.2] dual-write
      endif

! --- dry weight of roots at soil surface; oxygen module
      W_root_ss = afgen (wrtb,(2*magrs),dvs)

! --- ratio root total respiration / maintenance respiration; oxygen module
      max_resp_factor = afgen (mrftb,(2*magrs),dvs)

      case (4)
          
! --- root extension
      if (swrd.eq.1) then
        rdpot = afgen (rdtb,22,dvs)
        rdpot = min(rdpot,rdm)
        rd    = rdpot
      else
        rrpot = min (rdm-rdpot,rri)
        ! SS-ATM Phase 2 Task A-2.3: ptra read from state%atmosphere (atmosphere home).
        if (state%atmosphere%ptra.lt.nihil) rrpot = 0.0d0
        rdpot = rdpot + rrpot

        rr = min (rdm-rd,rri)
        if (state%atmosphere%ptra.lt.nihil .or.             &
     &      (present(state) .and.                           &
     &       state%soilwater%flWrtNonox)) rr = 0.0d0
        if (swdmi2rd.eq.1 .and. state%atmosphere%ptra.ge.nihil) rr = rr * state%soilwater%tra/state%atmosphere%ptra  ! [SS-SWC S-2.7]
        rd = rd + rr
      endif
      if (present(state)) then
        state%crop%common%rdpot = rdpot   ! [SS-GR-CROP A5.1]
        state%crop%common%rd    = rd      ! [SS-GR-CROP A5.1]
      endif

      return

      case default
         call fatalerr_collected ('CropFixed', 'Illegal value for TASK')
      end select

      end associate  ! tc_t1900 => state%timecontrol [TC-10]
      return
      end
! ----------------------------------------------------------------------
      subroutine grass(task, tsoil, state)
! ----------------------------------------------------------------------
!     Date               : November 2004
!     Purpose            : detailed grass growth routine
! SS-HEAT pre-Task-8: tsoil(:) non-optional dummy arg; callers pass
!   state%heat%tsoil. Threaded through to sumttd calls.
! SS-CRP C-2.5: state added (optional, intent in) to read flWrtNonox.
! SS-TC TC-10: t1900,daynr read via state%timecontrol tc_* aliases;
!   state threaded to sumttd for its own TC reads.
! SS-GR-ATM A5.1: intent changed inout to allow dual-write in cropgrass_init_from_config.
! [GR-CROP Phase B/7] narrow use variables
! ----------------------------------------------------------------------
      ! [SS-GR-CROPRT B8] DEFERRED — grass: all remaining variables globals:
      !   magrs, macp: array dims (could → swap_array_dimensions, deferred with rest)
      !   icrop, dvs, tsum, daycrop: computed in grass loop; dual-write to state%crop%common%
      !     but global canonical pending Phase C
      !   rid, tbase, tdwi, swinco: config params, no state home
      !   wlv/wst/wrt/wso pools + dwlv/dwst/dwrt: computed in grass loop (WOFOST biomass)
      !   wrtmax, wrtmin: root biomass bounds, no state home
      !   cf, ch, cfeic, lai, laipot, laiem, laiexp/pot, laimax: computed in grass loop
      !   cftb/chtb/cfeictb: state%crop%fixed homes (A5.2 dual-write) but grass has optional
      !     state — same constraint as wofost B7 / cropfixed B2; deferred to Phase C
      !   rdtb, slatb, rgrlai etc.: no state home; rd/rdpot etc.: computed in loop
      !   config switches (swrd etc.), physiology params (reltr, cvl etc.): no state home
      !   leaf arrays (lv/lvpot etc.), JvL params (twilt etc.): no state home
      !   cropstartact/endact/pot: state%crop%grass homes (A5) but written here — Phase C
      !   cuptgraz/pot, tagp/pot, tagpt/pot, seqgrazmow/pot, mowrest, dateharvest:
      !     state%crop%grass/common homes (A5) but written in grass loop — Phase C
      !   pgass/pgasspot: state%crop%wofost homes (A4 dual-write) but written here — Phase C
      !   perdl, dateharvest, lsda: output + harvest tracking, no state home
      !   tsoil: config-staging buffer, renamed to avoid clash with dummy arg
      use variables, only: &                                                ! [SS-GR-CROPRT B8] DEFERRED
        magrs, macp, icrop, dvs, rid, tsum, tbase, daycrop, tdwi, swinco, &
        wlv, wlvpot, wst, wstpot, wrt, wrtpot, wrtmax, wrtmin,           &
        dwlv, dwlvpot, dwrt, dwrtpot, dwst, dwstpot,                     &
        cf, ch, cfeic, lai, laipot, laiem, laiexp, laiexppot, laimax,    &
        cftb, chtb, cfeictb, rdtb, slatb, rgrlai, rlwtb, rfsetb,        &
        frtb, fltb, fstb, rdrrtb, rdrstb, kdif,                         &
        rd, rdpot, rdm, rdmax, rdi, rri, rdc, swrd, swrdc, swdmi2rd,    &
        swdrought, swcf, swgc, swinter, reltr,                           &
        cvl, cvr, cvs, q10, rmr, rml, rms, span, ssa, glaiex, glaiexpot, &
        lv, lvpot, lvage, lvagepot, sla, slapot, ilvold, ilvoldpot,     &
        twilt, wiltpoint, gwrt, siccapact, siccaplai,                   &
        cropstartact, cropendact, cropstartpot, cropendpot,             &
        cropstart, idaysgraz, idaysgrazpot, idregr, idregrpot,          &
        flgrazing, flgrazingpot, flharvest, flharvestpot,               &
        flhrvendact, flhrvendpot, flhydrlift,                           &
        daygrowth, daygrowthpot, grzdm, dewrest,                        &
        cuptgraz, cuptgrazpot, tagp, tagppot, tagpt, tagptpot,          &
        seqgrazmow, seqgrazmowpot, swtsum, iseqgm, iseqgmpot,           &
        iharvest, dmgrztb, dmmowtb, daysgrazingtab, uptgrazingtab,      &
        lossgrazingtab, lossgrztab, lossmowtab,                         &
        delayregrowthtab, zgrz, zmow,                                   &
        mowdm, mowrest, lossdm, plossdm, pmowdm, pgrzdm, pgass, pgasspot, &
        perdl, dateharvest, lsda,                                        &
        dummy_tsoil_gr_ => tsoil
      !! Rename config-staging tsoil to avoid clash with dummy arg tsoil.
      !! [SS-HEAT] Task 9: tsoil retained as config-staging buffer; global is not compute state.
      use array_utils, only: afgen
      use soilhydraulics_utils, only: watcon
      use rootextraction_mod, only: MatricFlux
      use swap_constants, only: tiny, nihil
      use error_mod, only: fatalerr_collected
      use swap_state_mod, only: swap_state_t
      ! GR-CROPWS Phase 0: sumttd, update_rootdistribution extracted to cropgrowth_helpers_mod
      use cropgrowth_helpers_mod, only: sumttd, update_rootdistribution

      implicit none

      type(swap_state_t), intent(inout), optional :: state

      integer   i1,task
      real(8), intent(in) :: tsoil(:)
      !! Soil temperature array from state%heat%tsoil, passed from CropGrowth.
      integer   idelaypot,idelay,i,swhydrlift

      real(8)   laicr,lasum,mres,grazlivinglv,grazlivinglvpot
      real(8)   admi,asrc,ccheck,cvf
      real(8)   dalv,delt,dmi
      real(8)   drrt,drst,dslv,dslv1,dslv2,dslvt,dteff
      real(8)   fcheck,fl,fr,fs,drlv
      real(8)   fysdel,gass,gla,glasol,grlv,grrt,grst
      real(8)   gwst,rest,rmres
      real(8)   slat,teff,twlv,twst
      real(8)   lasumpot,drst1,drst2
      real(8)   drst1pot,drst2pot
      real(8)   gasspot,rmrespot,mrespot,asrcpot,dmipot
      real(8)   admipot,grrtpot,drrtpot,gwrtpot,grlvpot,dslv1pot
      real(8)   dslv2pot,dslvpot,restpot,dalvpot,drlvpot
      real(8)   glasolpot,slatpot,glapot,grstpot,drstpot,gwstpot
      real(8)   dslvtpot,twlvpot,twstpot,tagpspot,tagps
      real(8)   dummy
      real(8)   dmharvest,dmlastharvest,dmgrazing
      real(8)   lsdb(100)
      real(8)   uptgraz,tagprest,lossgraz
      real(8)   uptgrazpot,lossgrazpot
      integer   daylastharvest,swharvest,daysgrazpot,daysgraz
      integer   swdmmow,swdmgrz
      character(len=11) tmp
      character(len=200) messag
      
! --- rooting
      real(8)   rrpot,rr
      
      integer   maxdaymow,maxdaygrz

      logical   flGrassGrowth
      logical   flDewoolingpot,flDewooling
      character(len=11) ::  dateGrassGrowth
      
!     In case of swlossmow = 1 or swlossgrz = 1: check work-ablility
      integer   nodmow,nodgrz
      integer   swlossmow,swlossgrz
      real(8)   fralossmow
      real(8)   fralossgrz
      real(8)   drz1
      logical   flearlyhrvendpot,flearlyhrvendact
      
      parameter (delt=1.0d0)

      save
! ----------------------------------------------------------------------
      ! SS-TC TC-10: t1900,daynr read via state%timecontrol tc_* aliases.
      ! [SS-BMI2 Task 4] tstart added to associate
      ! [SS-GR-ATM B.5] at_tav alias for tav read migration
      associate( &
        tc_t1900 => state%timecontrol%t1900,  &  ! TC-10
        tc_daynr => state%timecontrol%daynr,  &  ! TC-10
        tstart   => state%timecontrol%tstart, &  ! [SS-BMI2 Task 4]
        at_tav   => state%atmosphere%Tav      &  ! [SS-GR-ATM B.5]
      )

      select case (task)
      case (1)

! === initialization at start of crop =========================================

! --- read grass input data: dispatch on per-rotation typed-config cache
!     (ADR 0016). Falls back to legacy reader for rotations whose
!     .crp.toml is not yet authored. Teardown: end of Phase 4 removes
!     the else-branch.
      block
         use crop_config_global_mod, only: crop_config_global
         use cropgrass_init_mod,     only: cropgrass_init_from_config
         logical :: use_cache
         use_cache = .false.
         if (associated(crop_config_global)) then
            if (allocated(crop_config_global%rotation_loaded)) then
               if (icrop >= 1 .and. icrop <= size(crop_config_global%rotation_loaded)) then
                  if (crop_config_global%rotation_loaded(icrop)) then
                     ! Defense-in-depth: only dispatch to cache when the schema
                     ! is fully authored (case 4 + case 2 have amaxtb; the
                     ! hupselbrook skeleton does not — Phase 4 will fill it).
                     if (allocated(crop_config_global%rotation_grass(icrop)%amaxtb)) &
                        use_cache = .true.
                  end if
               end if
            end if
         end if
         if (use_cache) then
            associate(cfg => crop_config_global%rotation_grass(icrop))
               swharvest      = cfg%swharv
               dmharvest      = cfg%dmharvest
               daylastharvest = int(cfg%daylastharvest)
               dmlastharvest  = cfg%dmlastharvest
               swdmmow        = cfg%swdmmow
               maxdaymow      = cfg%maxdaymow
               swlossmow      = cfg%swlossmow
               swlossgrz      = cfg%swlossgrz
               swdmgrz        = cfg%swdmgrz
               maxdaygrz      = cfg%maxdaygrz
               dmgrazing      = cfg%dmgrazing
               LSDb           = 0.0d0   ! grazing stub-guarded; populated via daysgrazingtab/uptgrazingtab/lossgrazingtab by init
               tagprest       = cfg%tagprest
               swhydrlift     = 0       ! swdrought=2 stub-errored; mirror cropfixed/cropwofost default
               call cropgrass_init_from_config(cfg, icrop, &
                  state%timecontrol%tend, state%timecontrol%tstart, state)  ! [SS-BMI2 Task 4] [SS-GR-ATM A5.1]
            end associate
         else
            ! ADR 0016 cache-miss: typed config required for type=3 rotations.
            ! No silent legacy fallback — the user must author cropgrass.crp.toml.
            call fatalerr_collected('cropgrowth/Grass', &
               'cropgrass rotation has no loaded .crp.toml — author the file or use the legacy executable.')
         end if
      end block

! --- sequence of harvest by mowing, dewooling and grazing
      seqgrazmowpot = seqgrazmow
      if (present(state)) state%crop%grass%seqgrazmowpot = seqgrazmowpot   ! [SS-GR-CROP A5.1]

! --- development stage (not used by Grassland, instead Daynrs are used)
      dvs = -99.99d0
      if (present(state)) state%crop%common%dvs = dvs   ! [SS-GR-CROP A5.1]

! --- maximum rooting depth
      if (swrd.eq.1) then
        rdm = rdmax
      elseif (swrd.eq.2) then
        rdm = min(rdmax,rdc)
      elseif (swrd.eq.3) then
        rdc = afgen (rlwtb,22,wrtmax)
        rdm = min(rdmax,rdc)
      endif
      if (present(state)) state%crop%common%rdm = rdm   ! [SS-GR-CROP A5.1]

! --- skip next initialization if crop parameters are read from *.END file
      if (tc_t1900 - tstart .gt. tiny .or. swinco .ne. 3 .or.          &
     &   dabs(tc_t1900 - cropstart(icrop)) .lt. tiny) then

        iseqgm = 1
        iseqgmpot = iseqgm

! ---   initial values of crop parameters
        rid = dble(daycrop)
        fr = afgen (frtb,30,rid)
        fl = afgen (fltb,30,rid)
        fs = afgen (fstb,30,rid)
        sla(1) = afgen (slatb,30,rid)
        lvage(1) = 0.d0
        ilvold = 1
        idregr = 0
        slapot(1) = afgen (slatb,30,rid)
        lvagepot(1) = 0.d0
        ilvoldpot = 1
        idregrpot = 0

! ---   initial state variables of the crop
        wrt = fr*tdwi
        wrtmin = wrt / 10000 ! minimum root weigth at relative depth is set to 1% of the initial value
        wrtpot = wrt
        wst = fs*(1.0d0-fr)*tdwi
        wstpot = wst
        wlv = laiem/sla(1)
        wlvpot = wlv
        
!     KRO-BOO-20160403: intro because comparison with Wofost
        laiem = wlv*sla(1)  ! is not input !
        lv(1) = wlv
        lvpot(1) = lv(1)
        lasum = laiem
        lasumpot = lasum     
        glaiex = 0.0d0
        glaiexpot = 0.0d0
        laiexp = laiem
        laiexppot = laiem
        laimax = laiem
        lai = lasum+ssa*wst
        if (present(state)) state%crop%lai = lai   ! [SS-GR-ATM A5.2] dual-write
        laipot = lai
        dwrt = 0.d0
        dwrtpot = dwrt
        dwlv = 0.d0
        dwlvpot = dwlv
        dwst = 0.d0
        dwstpot = dwst

        daygrowth    = 0
        daygrowthpot = 0

! ---   actual rooting depth
        if (swrd.eq.1) then
          rd = afgen (rdtb,22,rid)
          rd = min(rd,rdm)
        elseif (swrd.eq.2) then
          rd = min(rdi,rdm)
        elseif (swrd.eq.3) then
          rdi = afgen (rlwtb,22,wrt)
          rd = min(rdi,rdm)
        endif
        rdpot = rd
        
! ---   initial summation variables of the crop
        tagp = wlv+wst
        tagppot = tagp
        tagpt = 0.0d0
        tagptpot = 0.0d0
        cuptgraz = 0.0d0
        cuptgrazpot = 0.0d0
        tsum = 0.0d0
        
        cropstartpot     = rid
        cropstartact     = rid
        flhrvendpot      = .false.
        flearlyhrvendpot = .false.
        ! [SS-GR-CROP A5.1] mirror grass init-time state
        if (present(state)) then
          state%crop%common%tsum        = tsum
          state%crop%common%rd          = rd
          state%crop%common%rdpot       = rdpot
          state%crop%common%laipot      = laipot
          state%crop%wofost%wrt         = wrt
          state%crop%wofost%wrtpot      = wrtpot
          state%crop%wofost%wst         = wst
          state%crop%wofost%wstpot      = wstpot
          state%crop%wofost%wlv         = wlv
          state%crop%wofost%wlvpot      = wlvpot
          state%crop%wofost%dwrt        = dwrt
          state%crop%wofost%dwrtpot     = dwrtpot
          state%crop%wofost%dwlv        = dwlv
          state%crop%wofost%dwlvpot     = dwlvpot
          state%crop%wofost%dwst        = dwst
          state%crop%wofost%dwstpot     = dwstpot
          state%crop%wofost%tagp        = tagp
          state%crop%wofost%tagppot     = tagppot
          state%crop%wofost%tagpt       = tagpt
          state%crop%wofost%tagptpot    = tagptpot
          state%crop%common%cuptgraz    = cuptgraz
          state%crop%common%cuptgrazpot = cuptgrazpot
          state%crop%grass%cropstartpot = cropstartpot
          state%crop%grass%cropstartact = cropstartact
        endif
        
        if (swtsum.eq.0) then
          flGrassGrowth = .true.
        else
          flGrassGrowth = .false.  
        endif
        if (swtsum.eq.2) then
          ! SS-HEAT pre-Task-8: pass tsoil from state%heat%tsoil (via dummy arg) to sumttd.
          ! SS-TC TC-10: pass state so sumttd can read t1900,date via state%timecontrol.
          call sumttd('initial',flGrassGrowth,dateGrassGrowth,tsoil,state)
        endif

! --- end skip above initialization if crop parameters are read from *.END file
      endif

      if (swcf.ne.3) then
        cf = afgen (cftb,(2*magrs),rid)
        ch = afgen (chtb,(2*magrs),rid)
      else
        cf        = afgen (cftb,(2*magrs),lai)
        cfeic     = afgen (cfeictb,(2*magrs),lai)
        ch        = afgen(chtb,(2*magrs),lai)
      endif
      if (present(state)) then
        state%crop%common%cf = cf   ! [SS-GR-CROP A5.1]
        state%crop%common%ch = ch   ! [SS-GR-CROP A5.1]
        if (swcf.eq.3) state%crop%fixed%cfeic = cfeic   ! [SS-GR-CROP A5.1]
      endif

! --- initial storage on canopy
      if (swinter.eq.3) then
        siccapact = siccaplai*lai
        if (present(state)) state%atmosphere%siccapact = siccapact   ! [SS-GR-ATM A5.2] dual-write
      endif

! --- initialize matric flux potential (SS-CRP C-2.5: hroot/hleaf/mfluxtable
!     init moved to CropGrowth dispatcher which has access to state).
      if (swdrought .eq. 2) then
        if (swhydrlift .eq. 1) then
          flhydrlift = .true.
        else
          flhydrlift = .false.
        endif
        do i = 1,state%mesh%numnod  ! [GR-BH C7]
         twilt(i) = watcon(wiltpoint, &
                            state%soilwater%vg_params(i), &
                            state%soilwater%iHWCKmodel(state%soilwater%layer(i)), &
                            i, state%soilwater)                    ! [SS-GR-UTILS Task 5]
        enddo
      endif

! --- harvest
!     initialise 
      if (swharvest.eq.2) then
        iharvest = 1
        do while (tc_t1900 .gt. dateharvest(iharvest))
          iharvest = iharvest + 1
        enddo
      endif      
      
! --- Find node for monitoring work-ability
      if (swlossmow .eq. 1) then
         
        ! Find node for monitoring work-ability during mowing
        nodmow = 1
        drz1       = -1.d0 * zmow - state%mesh%dz(nodmow)  ! [GR-BH C7]
        do while (drz1 .gt. 0.d0)
          nodmow   = nodmow + 1
          drz1 = drz1 - state%mesh%dz(nodmow)  ! [GR-BH C7]
        enddo
      
      endif   
      
      if (swlossgrz .eq. 1) then
        
        ! Find node and layer for monitoring work-ability at start of grazing
        nodgrz = 1
        drz1       = -1.d0 * zgrz - state%mesh%dz(nodgrz)  ! [GR-BH C7]
        do while (drz1 .gt. 0.d0)
          nodgrz   = nodgrz + 1
          drz1 = drz1 - state%mesh%dz(nodgrz)  ! [GR-BH C7]
        enddo
         
      endif
      
      return

      case (2)

! === calculate potential rate and state variables ======================================

! --- rates of change of the grass variables ---------------------------------------------

      rid = dble(daycrop)
      
! --- check end of harvest
      if (flhrvendpot) then
        if (flearlyhrvendpot) then
          cropstartpot  = rid - 1.d0
        else
          cropstartpot  = rid
        endif
        pmowdm        = 0.d0
        pgrzdm        = 0.d0
        plossdm       = 0.d0
        if (present(state)) then
          state%crop%grass%cropstartpot = cropstartpot   ! [SS-GR-CROP A5.1]
          state%crop%wofost%plossdm     = plossdm        ! [SS-GR-CROP A5.1]
        endif
      endif
      flhrvendpot      = .false.
      flearlyhrvendpot = .false.

! --- grass growth initiated by tsum from 1st day of calendar year
      tsum = tsum + max(0.0d0,at_tav)  ! [SS-GR-ATM B.5]
      if (present(state)) state%crop%common%tsum = tsum   ! [SS-GR-CROP A5.1]
      if (.not. flGrassGrowth) then
        
        ! grass growth initiated by tsum
        if (swtsum.eq.1) then
          if (tsum.ge.200.d0) then
            flGrassGrowth = .true.
          endif
        endif
        
        ! grass growth initiated by temperature, time and depth
        if (swtsum.eq.2) then
          ! SS-HEAT pre-Task-8: pass tsoil from state%heat%tsoil (via dummy arg) to sumttd.
          ! SS-TC TC-10: pass state so sumttd can read t1900,date via state%timecontrol.
          if (dateGrassGrowth.eq.'undefined') call sumttd('dynamic',flGrassGrowth,dateGrassGrowth,tsoil,state)
        endif
      
        ! check if grass growth has started
        if (flGrassGrowth) then
          cropstartpot = rid
          cropstartact = rid
        endif

      endif
      
! --- skip in case of: tsum<tsum200, or 3 criteria (tsummttd), or regrowth
      if (flGrassGrowth .and. daycrop.ge.idregrpot) then

! ===   daily dry matter production ===

        gasspot = pgasspot

! ---   respiration and partitioning of carbohydrates between growth and
! ---   maintenance respiration
        rmrespot=(rmr*wrtpot+rml*wlvpot+rms*wstpot)*afgen(rfsetb,30,rid)
        teff = q10**((at_tav-25.0d0)/10.0d0)  ! [SS-GR-ATM B.5]
        mrespot = min (gasspot,rmrespot*teff)
        asrcpot = gasspot-mrespot

! ---   partitioning factors
        fr = afgen(frtb,30,rid)
        fl = afgen(fltb,30,rid)
        fs = afgen(fstb,30,rid)
! ---   check on partitioning
        fcheck = fr+(fl+fs)*(1.0d0-fr) - 1.0d0
        if (dabs(fcheck).gt.0.0001d0) then
          write(tmp,'(f6.3)') rid
          tmp = adjustl (tmp)
          Messag ='The sum of partitioning factors for leaves, stems'// &
     &    ' and storage organs is not equal to one at time '            &
     &    //trim(tmp)//'.'
          call fatalerr_collected ('grass_pot',messag)
        endif

! ---   dry matter increase
        cvf = 1.0d0/((fl/cvl+fs/cvs)*(1.0d0-fr)+fr/cvr)
        dmipot = cvf*asrcpot

! ---   check on carbon balance
        ccheck = (gasspot-mrespot-(fr+(fl+fs)*(1.0d0-fr))*dmipot/cvf)   &
     &         /max(0.0001d0,gasspot)
        if (dabs(ccheck).gt.0.0001d0) then
          Messag ='The carbon balance is not correct'
          call fatalerr_collected ('grass_pot',messag)
        endif


! ===   growth rate by plant organ ===

! ---   growth rate roots and aerial parts

        grrtpot = fr*dmipot
        ! in case of SWRD = 3: after reaching maximum live weight of wrtmax, the
        ! growth of the roots is balanced by the death of root tissue
        if (swrd.eq.3 .and. wrtpot.gt.wrtmax) then
          drrtpot = grrtpot
          drrtpot = max(drrtpot,wrtpot*afgen (rdrrtb,30,rid))
        else  
          drrtpot = wrtpot*afgen (rdrrtb,30,rid)
        endif  
        gwrtpot = grrtpot - drrtpot

! ---   growth rate leaves

! ---   weight of new leaves
        admipot = (1.0d0-fr)*dmipot
        grlvpot = fl*admipot

! ---   death of leaves due to water stress or high lai
        dslv1pot = 0.0d0
        laicr = 3.2d0/kdif
        dslv2pot=wlvpot*max(0.0d0,                                      &
     &                  min(0.03d0,0.03d0*(laipot-laicr)/laicr))
        dslvpot = max (dslv1pot,dslv2pot) 

! ---   death of leaves due to exceeding life span;
! ---   leaf death is imposed on array until no more leaves have
! ---   to die or all leaves are gone

        restpot = dslvpot*delt
        i1 = ilvoldpot

        do while (restpot.gt.lvpot(max(i1,1)).and.i1.ge.1)
          restpot = restpot-lvpot(i1) 
          i1 = i1-1
        enddo

! ---   check if some of the remaining leaves are older than span,
! ---   sum their weights

        dalvpot = 0.0d0
        if (lvagepot(max(i1,1)).gt.span.and.restpot.gt.0.and.           &
     &                          i1.ge.1) then
          dalvpot = lvpot(i1)-restpot
          restpot = 0.0d0
          i1 = i1-1
        endif

        do while (i1.ge.1.and.lvagepot(max(i1,1)).gt.span)
          dalvpot = dalvpot+lvpot(i1)
          i1 = i1-1
        enddo

        dalvpot = dalvpot/delt

! ---   death rate leaves and growth rate living leaves
        drlvpot   = dslvpot+dalvpot

! ---   leaf area not to exceed exponential growth curve
        slatpot = afgen (slatb,30,rid)
        if (laiexppot.lt.6.0d0) then
          dteff = max (0.0d0,at_tav-tbase)  ! [SS-GR-ATM B.5]
          glaiexpot = laiexppot*rgrlai*dteff
! ---   source-limited increase in leaf area
          glasolpot = grlvpot*slatpot
          glapot = min (glaiexpot,glasolpot)
! ---   adjustment of specific leaf area of youngest leaf class
          if (grlvpot.gt.0.0d0) slatpot = glapot/grlvpot
        endif  

! ---   growth rate stems
        grstpot = fs*admipot
! ---   death of stems due to water stress is zero in case of potential growth
        drst1pot = 0.0d0
! ---   death of stems due to ageing
        drst2pot = afgen (rdrstb,30,rid)*wstpot
        drstpot = (drst1pot+drst2pot)/delt 
        gwstpot = grstpot-drstpot

! ----  integrals of the crop --------------------------------------------

!       set growing period after previous harvest        
        daygrowthpot = daygrowthpot + 1

!       Check trigger to start mowing event
        if (seqgrazmowpot(iseqgmpot) .eq. 2) then

          flharvestpot = .false.
            
          ! use dry matter threshold
          if (swharvest .eq. 1) then 
      
            ! use of fixed threshold
            if (swdmmow .eq. 1) then
              if (tagppot .gt. dmharvest .or. (tc_daynr .gt. daylastharvest  &
     &          .and. tagppot .gt. dmlastharvest)) then
                flharvestpot = .true.
              endif

            ! use of flexible threshold
            elseif (swdmmow .eq. 2) then
              dmharvest = afgen(dmmowtb,20,rid)
              if (tagppot .gt. dmharvest .or.                           &
     &                 (daygrowthpot .gt. maxdaymow .and. iseqgmpot .gt. 1)) then
                flharvestpot = .true.
              endif
            endif

          ! use fixed dates
          elseif (swharvest .eq. 2) then
            if(tc_t1900 .gt. dateharvest(iharvest)) then
              flharvestpot = .true.
            endif
          endif

!         In case mowing is triggered: Growth is initialized again and the weight of the sward is stored
          if (flharvestpot) then
            iseqgmpot = iseqgmpot + 1
            slapot(1) = afgen (slatb,30,rid)
            fl = afgen (fltb,30,rid)
            fs = afgen (fstb,30,rid)
            wlvpot = mowrest / (1.d0 + (fs/fl))
            wstpot = fs/fl*wlvpot
            dwlvpot = 0.0d0
            dwstpot = 0.0d0
            lvagepot(1) = 0.0d0
            ilvoldpot = 1
            lasumpot = wlvpot * slapot(1)
            laiexppot = lasumpot
            lvpot(1) = wlvpot
            
            gwstpot = 0.0d0
            gwrtpot = 0.0d0
            drlvpot = 0.0d0
            drstpot = 0.0d0
            drrtpot = 0.0d0
            
            daygrowthpot = 0
            
!           losses due to treading
            fralossmow = 0.d0
            if (swlossmow.eq.1) then
              fralossmow = afgen(lossmowtab,200,state%soilwater%h(nodmow))  ! [SS-SWC S-2.7]
            end if

!           harvest
            tagpspot = max(0.0d0,(tagppot-(wlvpot+dwlvpot+wstpot+dwstpot)))
            tagptpot = tagptpot + tagpspot * (1.d0 - FraLossMow)

            cropendpot  = rid
            flhrvendpot = .true.
            pmowdm      = tagpspot * (1.d0 - FraLossMow)
            plossdm     = tagpspot * FraLossMow
            
!           set regrowth delay
            idelaypot = int(afgen(DelayRegrowthTab,200,tagpspot))
            idregrpot = daycrop + idelaypot

          endif          
          
!       Check trigger to start grazing event          
        else if (seqgrazmowpot(iseqgmpot) .eq. 1 .or. seqgrazmowpot(iseqgmpot) .eq. 3) then

          flharvestpot = .false.
          
          if (.not. flgrazingpot) then
              
            ! use dry matter threshold
            if (swharvest .eq. 1) then 
            
              ! use of fixed threshold
              if (swdmgrz .eq. 1) then   
                if (tagppot .gt. dmgrazing) then
                  flharvestpot = .true.
                endif
              
              ! use of flexible threshold
              elseif (swdmgrz .eq. 2) then 
                dmgrazing = afgen(dmgrztb,20,rid)
                if (tagppot .gt. dmgrazing .or.                           &
     &                 (daygrowthpot .gt. maxdaygrz .and. iseqgmpot .gt. 1)) then
                  flharvestpot = .true.
                endif
              endif
              
            ! use fixed dates
            elseif (swharvest .eq. 2) then
              if(tc_t1900 .gt. dateharvest(iharvest)) then
                flharvestpot = .true.
              endif
            endif
          endif

!         In case grazing is triggered (or still occurs):
          if (flharvestpot .or. flgrazingpot) then
          
!           Amount of grazing kg/ha DM based on livestock density (Handboek Melkveehouderij 2013)
            uptgrazpot = lsda(iseqgmpot) *                                 &
     &                         afgen(uptgrazingtab,200,lsda(iseqgmpot))

!           Amount of shoots lost (kg/ha DM) due to droppings and treading during grazing  
            lossgrazpot = lsda(iseqgmpot) *                                &
     &                         afgen(lossgrazingtab,200,lsda(iseqgmpot))

!           Extra losses due to treading in case pressure head is insufficient
            fralossgrz = 0.d0
            if (swlossgrz.eq.1) then
              fralossgrz = afgen(lossgrztab,200,state%soilwater%h(nodgrz))  ! [SS-SWC S-2.7]
            end if
            lossgrazpot = lossgrazpot + tagppot * fralossgrz

!           Initialise Count nr of days with grazing
            if(.not. flgrazingpot) then
              daygrowthpot   = 0
              idaysgrazpot   = 0
              flDewoolingpot = .false.
            endif
            
!           verify if uptake is possible: tagprest should remain after grazing
            if ((tagppot - uptgrazpot - lossgrazpot) .gt. tagprest) then
              
              flgrazingpot = .true.
              cuptgrazpot  = cuptgrazpot + uptgrazpot
          
!             distribute grazing over stems and leaves (living and dead parts)
              wstpot  = wstpot  - (uptgrazpot+lossgrazpot) * wstpot  / tagppot
              dwstpot = dwstpot - (uptgrazpot+lossgrazpot) * dwstpot / tagppot
              dwlvpot = dwlvpot - (uptgrazpot+lossgrazpot) * dwlvpot / tagppot
              grazlivinglvpot =   (uptgrazpot+lossgrazpot) * wlvpot  / tagppot
          
!             reduce leave weights
              i1 = ilvoldpot
              do while (grazlivinglvpot .gt. 0 .and. i1 .ge. 1)
                if (grazlivinglvpot .ge. lvpot(i1)) then
                  grazlivinglvpot = grazlivinglvpot - lvpot(i1)
                  lvpot(i1) = 0.0d0
                  i1 = i1 - 1
                else
                  lvpot(i1) = lvpot(i1) - grazlivinglvpot
                  grazlivinglvpot = 0.d0
                endif
              enddo
          
!             harvest during total grazing event
              cropendpot = rid
              pgrzdm     = pgrzdm + uptgrazpot
              plossdm    = tagppot * fralossgrz
              
!             Check number of days with grazing
              daysgrazpot  = int(afgen(daysgrazingtab,200,lsda(iseqgmpot)))
              idaysgrazpot = idaysgrazpot + 1
              if(idaysgrazpot .eq. daysgrazpot) then
                flgrazingpot = .false.
                flhrvendpot  = .true.
                if (seqgrazmowpot(iseqgmpot) .eq. 3) then
                  flDewoolingpot  = .true.
                endif
                daygrowthpot = 0
                iseqgmpot = iseqgmpot + 1
              endif

!           Also end grazing when not enough grass remains on the field
            elseif (flgrazingpot .or. swharvest .eq. 2) then
              flgrazingpot     = .false.
              flhrvendpot      = .true.
              flearlyhrvendpot = .true.
              if (seqgrazmowpot(iseqgmpot) .eq. 3 .and. tagppot .gt. dewrest) then
                flDewoolingpot   = .true.
                flearlyhrvendpot = .false.
              endif
              daygrowthpot = 0
              iseqgmpot = iseqgmpot + 1
            endif

!           Assumption: no delay in regrowth during and after grazing (without dewooling)
            idregrpot = daycrop

!           Dewooling after grazing event            
            if (flDewoolingpot) then

              slapot(1) = afgen (slatb,30,rid)
              fl = afgen (fltb,30,rid)
              fs = afgen (fstb,30,rid)
              wlvpot = dewrest / (1.d0 + (fs/fl))
              wstpot = fs/fl*wlvpot
              dwlvpot = 0.0d0
              dwstpot = 0.0d0
              lvagepot(1) = 0.0d0
              ilvoldpot = 1
              lasumpot = wlvpot * slapot(1)
              laiexppot = lasumpot
              lvpot(1) = wlvpot
              
              gwstpot = 0.0d0
              gwrtpot = 0.0d0
              drlvpot = 0.0d0
              drstpot = 0.0d0
              drrtpot = 0.0d0
              
!             Assumption: one day delay in regrowth after grazing
              idregrpot = daycrop + 1
              
            endif
            
          endif
          
        endif
        
        if (daycrop .ge. idregrpot) then

! ---     physiologic ageing of leaves per time step
          fysdel = max (0.0d0,(at_tav-tbase)/(35.0d0-tbase))  ! [SS-GR-ATM B.5]

! ---     leaf death is imposed on array untill no more leaves have to die or all leaves are gone

          dslvtpot = dslvpot*delt
          i1 = ilvoldpot
           do while (dslvtpot.gt.0.and.i1.ge.1)
            if (dslvtpot.ge.lvpot(i1)) then
              dslvtpot = dslvtpot-lvpot(i1)
              lvpot(i1) = 0.0d0
              i1 = i1-1
            else
              lvpot(i1) = lvpot(i1)-dslvtpot
              dslvtpot = 0.0d0
            endif
          enddo

          if(i1.gt.0) then
            do while (lvagepot(max(i1,1)) .gt. span .and. i1 .ge. 1)
              lvpot(i1) = 0.0d0
              i1 = i1-1
            enddo
          endif
          ilvoldpot = i1

! ---     shifting of contents, integration of physiological age
          do i1 = ilvoldpot,1,-1
            lvpot(i1+1) = lvpot(i1)
            slapot(i1+1) = slapot(i1)
            lvagepot(i1+1) = lvagepot(i1)+fysdel*delt
          enddo
          ilvoldpot = ilvoldpot+1

! ---     new leaves in class 1
          lvpot(1) = grlvpot*delt
          slapot(1) = slatpot
          lvagepot(1) = 0.d0

! ---     calculation of new leaf area and weight
          lasumpot = 0.d0
          wlvpot = 0.d0
          do i1 = 1,ilvoldpot
            lasumpot = lasumpot+lvpot(i1)*slapot(i1)
            wlvpot = wlvpot+lvpot(i1)
          enddo

          laiexppot = laiexppot+glaiexpot*delt

        endif

! ---   dry weight of living plant organs
        wrtpot = wrtpot+gwrtpot*delt
        wstpot = wstpot+gwstpot*delt

! ---   dry weight of dead plant organs (roots,leaves & stems)
        dwrtpot = dwrtpot+drrtpot*delt
        dwlvpot = dwlvpot+drlvpot*delt
        dwstpot = dwstpot+drstpot*delt

! ---   dry weight of dead and living plant organs
        twlvpot = wlvpot+dwlvpot
        twstpot = wstpot+dwstpot
        tagppot = twlvpot+twstpot

! ---   leaf area index
        laipot = lasumpot+ssa*wstpot
!       prevent immediate lai reduction at emergence
!       KRO-BOO-20160403: suppressed because deviates from Wofost
!       laipot = max(laipot, laiem)

        ! root extension
        if (swrd.eq.1) then
          rdpot = afgen (rdtb,22,rid)
          rdpot = min(rdpot,rdm)
        elseif (swrd.eq.2) then
          rrpot = min (rdm-rdpot,rri)
          if (fr.le.0.0d0 .or. pgasspot.lt.1.0d0) rrpot = 0.0d0
          rdpot = rdpot + rrpot
        elseif (swrd.eq.3) then
          rdpot = afgen (rlwtb,22,wrtpot)
          rdpot = min(rdpot,rdm)
        endif

      endif

      ! [SS-GR-CROP A5.1] mirror grass case(2) potential state
      if (present(state)) then
        state%crop%wofost%wlvpot        = wlvpot
        state%crop%wofost%wrtpot        = wrtpot
        state%crop%wofost%wstpot        = wstpot
        state%crop%wofost%dwrtpot       = dwrtpot
        state%crop%wofost%dwlvpot       = dwlvpot
        state%crop%wofost%dwstpot       = dwstpot
        state%crop%wofost%tagppot       = tagppot
        state%crop%common%laipot        = laipot
        state%crop%common%rdpot         = rdpot
        state%crop%wofost%pgasspot      = pgasspot
        state%crop%wofost%plossdm       = plossdm
        state%crop%common%cuptgrazpot   = cuptgrazpot
        state%crop%wofost%tagptpot      = tagptpot
        state%crop%grass%cropstartpot   = cropstartpot
        state%crop%grass%cropendpot     = cropendpot
      endif

      return

      case (3)

! === calculate actual rate and state variables ======================================

! --- check end of harvest
      if (flhrvendact) then
        if (flearlyhrvendact) then
          cropstartact  = rid - 1.d0
        else
          cropstartact  = rid
        endif
        mowdm        = 0.d0
        grzdm        = 0.d0
        lossdm       = 0.d0
        if (present(state)) then
          state%crop%grass%cropstartact = cropstartact   ! [SS-GR-CROP A5.1]
          state%crop%wofost%lossdm      = lossdm         ! [SS-GR-CROP A5.1]
        endif
      endif
      flhrvendact      = .false.
      flearlyhrvendact = .false.

! --- rates of change of the crop variables ---------------------------------------------
      
! --- skip in case of: tsum<tsum200, or 3 criteria (tsummttd), or regrowth
      if (flGrassGrowth .and. daycrop .ge. idregr) then

! ===   daily dry matter production ===

! ---   water stress reduction of pgass to gass
        ! SS-ATM Phase 2 Task A-2.3: ptra read from state%atmosphere (atmosphere home).
        if(dabs(state%atmosphere%ptra).lt.nihil) then
          reltr = 1.0d0
        else
          reltr = max(0.0d0,min(1.0d0,state%soilwater%tra/state%atmosphere%ptra))  ! [SS-SWC S-2.7]
        endif
        gass = pgass * reltr

! ---   respiration and partitioning of carbohydrates between growth and
! ---   maintenance respiration
        rmres = (rmr*wrt+rml*wlv+rms*wst)*afgen(rfsetb,30,rid)
        teff = q10**((at_tav-25.0d0)/10.0d0)  ! [SS-GR-ATM B.5]
        mres = min (gass,rmres*teff)
        asrc = gass-mres

! ---   partitioning factors (relevant for restart)
        fr = afgen(frtb,30,rid)
        fl = afgen(fltb,30,rid)
        fs = afgen(fstb,30,rid)

! ---   dry matter increase
        cvf = 1.0d0/((fl/cvl+fs/cvs)*(1.0d0-fr)+fr/cvr)
        dmi = cvf*asrc
! ---   check on carbon balance
        ccheck = (gass-mres-(fr+(fl+fs)*(1.0d0-fr))*dmi/cvf)            &
     &         /max(0.0001d0,gass)      
        if (dabs(ccheck).gt.0.0001d0) then
          Messag ='The carbon balance is not correct'
          call fatalerr_collected ('grass_act',messag)
        endif

! ===   growth rate by plant organ ===

! ---   growth rate roots and aerial parts
        ! in case of SWRD = 3: after reaching maximum live weight of wrtmax, the
        ! growth of the roots is balanced by the death of root tissue
        grrt = fr*dmi
        if (swrd.eq.3 .and. present(state) .and.                       &
     &      state%soilwater%flWrtNonox) grrt = 0.d0
        if (swrd.eq.3 .and. wrt.gt.wrtmax) then
          drrt = grrt
          drrt = max(drrt,wrt*afgen (rdrrtb,30,rid))
        else  
          drrt = wrt*afgen (rdrrtb,30,rid)
        endif  
        gwrt = grrt-drrt

! ---   growth rate leaves

! ---   weight of new leaves
        admi = (1.0d0-fr)*dmi        
        grlv = fl*admi

! ---   death of leaves due to water stress or high lai
        dslv1 = wlv*(1.0d0-reltr)*perdl
        laicr = 3.2d0/kdif
        dslv2 = wlv*max(0.0d0,min(0.03d0,0.03d0*(lai-laicr)/laicr))
        dslv = max (dslv1,dslv2) 

! ---   death of leaves due to exceeding life span;
! ---   leaf death is imposed on array until no more leaves have
! ---   to die or all leaves are gone

        rest = dslv*delt
        i1 = ilvold

        do while (rest.gt.lv(max(i1,1)).and.i1.ge.1)
          rest = rest-lv(i1) 
          i1 = i1-1
        enddo

! ---   check if some of the remaining leaves are older than span,
! ---   sum their weights

        dalv = 0.0d0
        if (lvage(max(i1,1)).gt.span.and.rest.gt.0.and.i1.ge.1) then
          dalv = lv(i1)-rest
          rest = 0.0d0
          i1 = i1-1
        endif

        do while (i1.ge.1.and.lvage(max(i1,1)).gt.span)
          dalv = dalv+lv(i1)
          i1 = i1-1
        enddo

        dalv = dalv/delt

! ---   death rate leaves and growth rate living leaves
        drlv   = dslv+dalv

! ---   physiologic ageing of leaves per time step
        slat = afgen (slatb,30,rid)

! ---   leaf area not to exceed exponential growth curve
        if (laiexp.lt.6.0d0) then
          dteff = max (0.0d0,at_tav-tbase)  ! [SS-GR-ATM B.5]
          glaiex = laiexp*rgrlai*dteff
! ---     source-limited increase in leaf area
          glasol = grlv*slat
          gla = min (glaiex,glasol)
! ---     adjustment of specific leaf area of youngest leaf class
          if (grlv.gt.0.0d0) slat = gla/grlv
        endif  

! ---   growth rate stems
        grst = fs*admi
! ---   death of stems due to water stress
        drst1 = wst*(1.0d0-reltr)*perdl
! ---   death of stems due to ageing
        drst2 = afgen (rdrstb,30,rid)*wst
        drst = (drst1+drst2)/delt 
        gwst = grst-drst

! ----  integrals of the crop --------------------------------------------

!       set growing period after previous harvest        
        daygrowth = daygrowth + 1

!       Check trigger to start mowing event
        if (seqgrazmow(iseqgm) .eq. 2) then
          
          flharvest = .false.   
            
          ! use dry matter threshold
          if (swharvest .eq. 1) then 
      
            ! use of fixed threshold
            if (swdmmow .eq. 1) then
              if (tagp .gt. dmharvest .or. (tc_daynr .gt. daylastharvest  &
     &          .and. tagp .gt. dmlastharvest)) then
                flharvest = .true.
              endif

            ! use of flexible threshold
            elseif (swdmmow .eq. 2) then
              dmharvest = afgen(dmmowtb,20,rid)
              if (tagp .gt. dmharvest .or.                           &
     &                 (daygrowth .gt. maxdaymow .and. iseqgm .gt. 1)) then
                flharvest = .true.
              endif
            endif
          
          ! use fixed dates
          elseif (swharvest .eq. 2) then
            if(tc_t1900 .gt. dateharvest(iharvest)) then
              iharvest = iharvest + 1
              flharvest = .true.
            endif
          endif
          
!       In case mowing is triggered: Growth is initialized again and the weight of the sward is stored
        if (flharvest) then
          iseqgm = iseqgm + 1
          sla(1) = afgen (slatb,30,rid)
          fl = afgen (fltb,30,rid)
          fs = afgen (fstb,30,rid)
          wlv = mowrest / (1.d0 + (fs/fl))
          wst = fs/fl*wlv
          dwlv = 0.0d0
          dwst = 0.0d0
          lvage(1) = 0.0d0
          ilvold = 1
          lasum = wlv * sla(1)
          laiexp = lasum
          lv(1) = wlv

          gwst = 0.0d0
          gwrt = 0.0d0
          drlv = 0.0d0
          drst = 0.0d0
          drrt = 0.0d0

          daygrowth = 0

!         losses due to treading
          FraLossMow = 0.d0
          if (swlossmow.eq.1) then
            FraLossMow = afgen(lossmowtab,200,state%soilwater%h(nodmow))  ! [SS-SWC S-2.7]
          end if
          
!         harvest
          tagps = max (0.0d0,(tagp-(wlv+dwlv+wst+dwst)))
          tagpt = tagpt + tagps * (1.d0 - fralossmow)

          cropendact  = rid
          flhrvendact = .true.
          mowdm   = tagps * (1.d0 - FraLossMow)
          lossdm  = tagps * FraLossMow
          
! ---     set regrowth delay
          idelay = int(afgen(DelayRegrowthTab,200,tagps))
          idregr = daycrop + idelay

        endif
          
!       Check trigger to start grazing event          
        else if (seqgrazmow(iseqgm) .eq. 1 .or. seqgrazmow(iseqgm) .eq. 3) then

          flharvest = .false.
            
          if (.not. flgrazing) then
            
            ! use dry matter threshold
            if (swharvest .eq. 1) then 
            
              ! use of fixed threshold
              if (swdmgrz .eq. 1) then   
                if (tagp .gt. dmgrazing) then
                  flharvest = .true.
                endif
              
              ! use of flexible threshold
              elseif (swdmgrz .eq. 2) then 
                dmgrazing = afgen(dmgrztb,20,rid)
                if (tagp .gt. dmgrazing .or.                           &
     &            (daygrowth .gt. maxdaygrz .and. iseqgm .gt. 1)) then
                  flharvest = .true.
                endif
              endif
            
            ! use fixed dates
            elseif (swharvest .eq. 2) then
              if(tc_t1900 .gt. dateharvest(iharvest)) then
                iharvest = iharvest + 1
                flharvest = .true.
              endif
            endif
          endif

!         In case grazing is triggered (or still occurs):
          if (flharvest .or. flgrazing) then
          
!           Amount of grazing kg/ha DM based on livestock density (Handboek Melkveehouderij 2013)
            uptgraz = lsda(iseqgm)*afgen(uptgrazingtab,200,lsda(iseqgm))            

!           Amount of shoots lost (kg/ha DM) due to droppings and treading during grazing  
            lossgraz = lsda(iseqgm) *                                &
     &                         afgen(lossgrazingtab,200,lsda(iseqgm))

!           Extra losses due to treading in case pressure head is insufficient
            fralossgrz = 0.d0
            if (swlossgrz.eq.1) then
              fralossgrz = afgen(lossgrztab,200,state%soilwater%h(nodgrz))  ! [SS-SWC S-2.7]
            end if
            lossgraz = lossgraz + tagp * fralossgrz

!           Initialise Count nr of days with grazing
            if(.not. flgrazing) then
              daygrowth   = 0
              idaysgraz   = 0
              flDewooling = .false.
            endif
            
!           verify if uptake is possible: tagprest should remain after grazing
            if ((tagp - uptgraz - lossgraz) .gt. tagprest) then
              
              flgrazing = .true.
              cuptgraz  = cuptgraz + uptgraz

!             distribute grazing over stems and leaves (living and dead parts)
              wst  = wst  -  (uptgraz+lossgraz) * wst  / tagp
              dwst = dwst -  (uptgraz+lossgraz) * dwst / tagp
              dwlv = dwlv -  (uptgraz+lossgraz) * dwlv / tagp
              grazlivinglv = (uptgraz+lossgraz) * wlv  / tagp
          
!             reduce leave weights
              i1 = ilvold
              do while (grazlivinglv .gt. 0 .and. i1 .ge. 1)
                if (grazlivinglv .ge. lv(i1)) then
                  grazlivinglv = grazlivinglv - lv(i1)
                  lv(i1) = 0.0d0
                  i1 = i1 - 1
                else
                  lv(i1) = lv(i1) - grazlivinglv
                  grazlivinglv = 0.d0
                endif
              enddo
          
!             harvest during total grazing event
              cropendact = rid
              grzdm      = grzdm + uptgraz
              lossdm     = tagp * fralossgrz
              
!             Check number of days with grazing
              daysgraz  = int(afgen(daysgrazingtab,200,lsda(iseqgm)))
              idaysgraz = idaysgraz + 1
              if(idaysgraz .eq. daysgraz) then
                flgrazing   = .false.
                flhrvendact = .true.
                if (seqgrazmow(iseqgm) .eq. 3) then
                  flDewooling  = .true.
                endif
                daygrowth = 0
                iseqgm = iseqgm + 1
              endif

!           Also end grazing when not enough grass remains on the field
            elseif (flgrazing .or. swharvest .eq. 2) then
              flgrazing        = .false.
              flhrvendact      = .true.
              flearlyhrvendact = .true.
              if (seqgrazmow(iseqgm) .eq. 3 .and. tagp .gt. dewrest) then
                flDewooling      = .true.
                flearlyhrvendact = .false.
              endif
              daygrowth = 0
              iseqgm = iseqgm + 1
            endif

!           Assumption: no delay in regrowth during and after grazing (without dewooling)
            idregr = daycrop

!           Dewooling after grazing event            
            if (flDewooling) then

              sla(1) = afgen (slatb,30,rid)
              fl = afgen (fltb,30,rid)
              fs = afgen (fstb,30,rid)
              wlv = dewrest / (1.d0 + (fs/fl))
              wst = fs/fl*wlv
              dwlv = 0.0d0
              dwst = 0.0d0
              lvage(1) = 0.0d0
              ilvold = 1
              lasum = wlv * sla(1)
              laiexp = lasum
              lv(1) = wlv
    
              gwst = 0.0d0
              gwrt = 0.0d0
              drlv = 0.0d0
              drst = 0.0d0
              drrt = 0.0d0
    
!             Assumption: one day delay in regrowth after grazing
              idregr = daycrop + 1

            endif
            
          endif
          
        endif

        if (daycrop .ge. idregr) then

! ---     physiologic ageing of leaves per time step
          fysdel = max (0.0d0,(at_tav-tbase)/(35.0d0-tbase))  ! [SS-GR-ATM B.5]

! ---     leaf death is imposed on array untill no more leaves have to die or all leaves are gone

          dslvt = dslv*delt
          i1 = ilvold
          do while (dslvt.gt.0.and.i1.ge.1)
            if (dslvt.ge.lv(i1)) then
              dslvt = dslvt-lv(i1)
              lv(i1) = 0.0d0
              i1 = i1-1
            else
              lv(i1) = lv(i1)-dslvt
              dslvt = 0.0d0
            endif
          enddo

          if(i1.gt.0) then
            do while (lvage(max(i1,1)).gt.span.and.i1.ge.1)
              lv(i1) = 0.0d0
              i1 = i1-1
            enddo
          endif
          ilvold = i1

! ---     shifting of contents, integration of physiological age
          do i1 = ilvold,1,-1
            lv(i1+1) = lv(i1)
            sla(i1+1) = sla(i1)
            lvage(i1+1) = lvage(i1)+fysdel*delt
          enddo
          ilvold = ilvold+1

! ---     new leaves in class 1
          lv(1) = grlv*delt
          sla(1) = slat
          lvage(1) = 0.d0 

! ---     calculation of new leaf area and weight
          lasum = 0.d0
          wlv = 0.d0
          do i1 = 1,ilvold
            lasum = lasum+lv(i1)*sla(i1)
            wlv = wlv+lv(i1)
          enddo

          laiexp = laiexp+glaiex*delt

        endif

! ---   dry weight of living plant organs
        wrt = wrt+gwrt*delt
        wst = wst+gwst*delt

! ---   dry weight of dead plant organs (roots,leaves & stems)
        dwrt = dwrt+drrt*delt
        dwlv = dwlv+drlv*delt
        dwst = dwst+drst*delt

! ---   dry weight of dead and living plant organs
        twlv = wlv+dwlv
        twst = wst+dwst
        tagp = twlv+twst

! ---   leaf area index
        lai = lasum+ssa*wst
        if (present(state)) state%crop%lai = lai   ! [SS-GR-ATM A5.2] dual-write
        laimax = max (lai,laimax)

! ---   update normalized cumulative root density based on root extraction or stress (cumdens)
        if (swrdc .eq. 1) call update_rootdistribution(state)
        
        ! root extension
        if (swrd.eq.1) then
          rd = afgen (rdtb,22,rid)
          rd = min(rd,rdm)
        elseif (swrd.eq.2) then
          rr = min (rdm-rd,rri)
          if (fr.le.0.0d0 .or. pgass.lt.1.0d0 .or.                    &
     &        (present(state) .and.                                      &
     &         state%soilwater%flWrtNonox)) rr = 0.0d0
          if (swdmi2rd.eq.1 .and. pgass.ge.1.0d0)              rr = rr * gass/pgass
          rd = rd + rr
        elseif (swrd.eq.3) then
          rd = afgen (rlwtb,22,wrt)
          rd = min(rd,rdm)
        endif
        if (present(state)) state%crop%common%rd = rd   ! [SS-GR-CROP A5.1]

! ---   set crop height and cropfactor
        if (swcf.ne.3) then
          cf = afgen (cftb,(2*magrs),rid)
          ch = afgen (chtb,(2*magrs),rid)
        else
          cf = afgen (cftb,(2*magrs),lai)
          cfeic = afgen (cfeictb,(2*magrs),lai)
          ch = afgen(chtb,(2*magrs),lai)
        endif
        if (present(state)) then
          state%crop%common%cf = cf   ! [SS-GR-CROP A5.1]
          state%crop%common%ch = ch   ! [SS-GR-CROP A5.1]
          if (swcf.eq.3) state%crop%fixed%cfeic = cfeic   ! [SS-GR-CROP A5.1]
        endif

! ---   update canopy storage capacity
        if (swinter.eq.3) then
          siccapact = siccaplai*lai
          if (present(state)) state%atmosphere%siccapact = siccapact   ! [SS-GR-ATM A5.2] dual-write
        endif

      endif

      ! [SS-GR-CROP A5.1] mirror grass case(3) actual state
      if (present(state)) then
        state%crop%wofost%wlv         = wlv
        state%crop%wofost%wrt         = wrt
        state%crop%wofost%wst         = wst
        state%crop%wofost%dwrt        = dwrt
        state%crop%wofost%dwlv        = dwlv
        state%crop%wofost%dwst        = dwst
        state%crop%wofost%tagp        = tagp
        state%crop%wofost%pgass       = pgass
        state%crop%wofost%lossdm      = lossdm
        state%crop%common%cuptgraz    = cuptgraz
        state%crop%wofost%tagpt       = tagpt
        state%crop%grass%cropstartact = cropstartact
        state%crop%grass%cropendact   = cropendact
        state%crop%common%tsum        = tsum
      endif

      return

      case default
         call fatalerr_collected ('Grass', 'Illegal value for TASK')
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

