! cropgrowth_helpers.f90
! GR-CROPWS Phase 0 Commit 0.1: leaf-level helpers extracted from cropgrowth.f90.
! Subroutines: nocrop, ArableLandGerm, FacCO2, cropoutput, update_rootdistribution,
!              sumttd, init_crop_output_buffer, build_crop_output_row,
!              cleanup_crop_output_buffer
! Pure relocation — no behavior change.
! [SS-GR-CROPWS A5]: Phase A audit — ArableLandGerm already intent(inout) with inline
!   writes; FacCO2 already writes state directly (B6); nocrop() has no state arg
!   (DEFERRED). No optional/present guards. No changes required.
! [GR-CROPWS B2]: reads migrated — PrepDelay + SowDelay compound reads in ArableLandGerm.
!   All other read sites in this file were already migrated (B3/B6/B9) or are mirror-line
!   RHS reads (stays). icrop already via state%crop%common%icrop.
! ----------------------------------------------------------------------
      module cropgrowth_helpers_mod
      implicit none
      private

      public :: nocrop
      public :: ArableLandGerm
      public :: FacCO2
      public :: cropoutput
      public :: update_rootdistribution
      public :: sumttd

      contains

! ----------------------------------------------------------------------
      subroutine cropoutput(task, state)
! ----------------------------------------------------------------------
!     Date               : Aug 2004
!     Purpose            : open and write crop output files
! SS-TC TC-10: state added (intent in); threaded to OutCropFixed /
!   OutWofost / OutGrass so they can read date,t via state%timecontrol.
! [SS-BMI2] inout: init/cleanup of crop_output_row buffer at top-level state.
! [GR-CROP Phase B/5] narrow use variables
! ----------------------------------------------------------------------

      ! [GR-CROP 2026-05-25] flCropOpenFile → state%crop%common%flCropOpenFile.
      ! outfil/pathwork/project read via state%cfg%general (this file's cfg snapshot).
      ! cropfil/crp remain on the bare-global side: writers/readers exist outside
      ! the crop subsystem (swapoutput.f90 reads crp, config_to_variables.f90
      ! populates cropfil). Imported here from variables until the file-IO arc.
      use variables, only: cropfil, crp
      use error_mod, only: fatalerr_collected
      use file_io_mod, only: file_open
      use swap_state_mod, only: swap_state_t
      implicit none

! --- local variables ------------------
      integer task   !,numcrop
      ! [SS-BMI2] inout: init/cleanup of crop_output_row buffer at top-level state.
      type(swap_state_t), intent(inout) :: state
      character(len=200) messag
      character(len=160) filnam,filtext

      associate( &
        crop => state%crop%common,   &  ! crop runtime (croptype, icrop, flCropOpenFile)
        time => state%timecontrol,   &  ! time control (headless, swheader)
        cfg_gen => state%cfg%general &  ! general config (outfil, pathwork, project)
      )

      select case (task)
      case (1)

! === open crop output file and write headers =====================

      ! [SS-BMI2] allocate crop output buffer (builder always runs; dynamic by croptype)
      call init_crop_output_buffer(state)

! --- open crop output file
      if (crop%flCropOpenFile) then

         if (.not. time%headless) then
! ---   open crop output file and write general header (*.crp)
            if (trim(cfg_gen%outfil).eq.trim(cropfil(1))) then
               Messag = 'The name of the input crop-file (''//trim(cropfil'//&
     &      '(icrop))//'') cannot be equal to the name of'                   &
     &      //'the output crop-file '//trim(cfg_gen%outfil)//' Adjust a filename !'
               call fatalerr_collected ('crops',messag)
            endif
            filnam = trim(cfg_gen%pathwork)//trim(cfg_gen%outfil)//'.crp'
            call file_open(crp, filnam, 'replace', 'write')
            filtext = 'output data of simple or detailed crop growth model'
            call writehead (crp,1,filnam,filtext,cfg_gen%project)

! ---   write header fixed crop growth
            if (crop%croptype(crop%icrop) .eq. 1) call OutCropFixed(1, state)

! ---   write header detailed crop growth
            if (crop%croptype(crop%icrop) .eq. 2) call OutWofost(1, state)

! ---   write header detailed grass growth
            if (crop%croptype(crop%icrop) .eq. 3) call OutGrass(1, state)
         end if

         crop%flCropOpenFile = .false.

      else

         if (.not. time%headless) then
! ---   header for second and subsequent crops
! ---   write header fixed crop growth
            if (crop%croptype(crop%icrop).eq.1 .and. time%swheader.eq.1) call OutCropFixed(1, state)

! ---   write header detailed crop growth
            if (crop%croptype(crop%icrop).eq.2 .and. time%swheader.eq.1) call OutWofost(1, state)

! ---   write header detailed grass growth
            if (crop%croptype(crop%icrop).eq.3 .and. time%swheader.eq.1) call OutGrass(1, state)
         end if

      endif

      case (2)

! --- write actual data ----------------------------------------------------
! [SS-BMI2] build crop output row buffer (placeholder; full build deferred to crop output refactor arc)
      call build_crop_output_row(state)

      if (.not. time%headless) then
! --- fixed crop file
         if (crop%croptype(crop%icrop) .eq. 1) call OutCropFixed(2, state)

! --- detailed crop growth
         if (crop%croptype(crop%icrop) .eq. 2) call OutWofost(2, state)

! --- detailed grass growth
         if (crop%croptype(crop%icrop) .eq. 3) call OutGrass(2, state)
      end if

      case (3)
! --- close crop output file ------------------------------------------------

      ! [SS-BMI2] headless guard: .crp file was only opened when not headless
      if (.not. time%headless) close (crp)

      ! [SS-BMI2] deallocate crop output buffer
      call cleanup_crop_output_buffer(state)

      case default
         call fatalerr_collected ('CropOutput', 'Illegal value for TASK')
      end select

      end associate
      return
      end subroutine cropoutput

! ----------------------------------------------------------------------
      subroutine nocrop (state)
! ----------------------------------------------------------------------

      ! [SS-GR-CROPRT B4] DEFERRED — nocrop: pure write site (sets globals to zero defaults).
      !   dvs now writes directly to state%crop%common%dvs (legacy global retired in dvs pilot).
      !   Other symbols still legacy; CropGrowth mirrors them to state immediately after the call.
      ! [SS-GR-CROPRT B4] nocrop now writes only to state (no variables imports needed)
      use swap_state_mod, only: swap_state_t
      implicit none
      type(swap_state_t), intent(inout) :: state

      state%crop%common%rd      = 0.0d0
      state%crop%common%rdpot   = 0.0d0
      state%crop%lai     = 0.0d0
      state%crop%common%laipot  = 0.d0
      state%crop%common%cf      = 0.0d0
      state%crop%common%ch      = 0.0d0
      state%crop%common%albedo  = 0.23d0
      state%crop%common%rsc     = 70.d0
      state%crop%common%tsum = 0.d0
      state%crop%common%dvs = 0.d0
      state%crop%wofost%cwdmpot = 0.d0
      state%crop%wofost%cwdm    = 0.d0
      state%crop%wofost%wsopot  = 0.d0
      state%crop%wofost%wso = 0.d0
      state%crop%wofost%wlvpot  = 0.d0
      state%crop%wofost%wlv = 0.d0
      state%crop%wofost%wstpot  = 0.d0
      state%crop%wofost%wst = 0.d0
      state%crop%wofost%wrtpot  = 0.d0
      state%crop%wofost%wrt = 0.d0

      return
      end subroutine nocrop

! ----------------------------------------------------------------------
      subroutine ArableLandGerm(task, tsoil, state)
! ----------------------------------------------------------------------
!     update             : December 2017
!     date               : December 2017
!     purpose            : Crop growth
! SS-HEAT pre-Task-8: tsoil(:) non-optional dummy arg; callers pass
!   state%heat%tsoil. Global tsoil excluded via rename.
! SS-SWC S-2.7: state added (intent in) for soil-water-core h reader cutover.
! [GR-CROP Phase B/5] narrow use variables
! ----------------------------------------------------------------------
      ! [GR-CROP 2026-05-25] ArableLandGerm fully migrated off germ-param globals.
      !   prep/sow params (hPrep/zPrep/MaxPrepDelay/hSow/zSow/zTempSow/TempSow/MaxSowDelay):
      !     read via crop_config_global%rotation_wofost(icrop)%preparation/sowing.
      !     dhPrep/dhSow/dtempSow become local scratch (Class E — case(2)/(3) are
      !     unreachable; cropgrowth.f90:255 stub-errors swprep!=0 OR swsow!=0).
      !   germ thresholds (tsumemeopt/tbasem/teffmx/hdrygerm/hwetgerm/zgerm/agerm):
      !     read via crop_config_global%rotation_wofost(icrop)%germination.
      !   bgerm/cgerm: derived locally from agerm/hdrygerm/hwetgerm/tsumemeopt.
      !   tsumgerm: state%crop%common%tsumgerm (accumulated runtime).
      !   flCropPrep/flCropSow/flCropGerm: legacy bare globals retained for cropgrowth.f90
      !     dispatcher reads (still gated on legacy bare names); dual-write keeps them in
      !     sync with state%crop%common.
      use variables, only: flCropPrep, flCropSow, flCropGerm
      use crop_config_global_mod, only: crop_config_global
      use swap_constants, only: small
      use error_mod, only: fatalerr_collected
      use swap_state_mod, only: swap_state_t
      implicit none

      integer  task,node
      real(8), intent(in) :: tsoil(:)
      !! Soil temperature array from state%heat%tsoil.
      type(swap_state_t), intent(inout) :: state
      !! State record; state%soilwater%h used for pressure-head reads. [SS-SWC S-2.7]
      real(8)  drz1,hrz1,pFz1
      real(8)  tsumemesub
      real(8)  dhPrep, dhSow, dtempSow      ! local scratch (former bare globals)
      real(8)  l_agerm, l_bgerm, l_cgerm    ! local copies / derived
      real(8)  l_tsumemeopt, l_hdrygerm, l_hwetgerm, l_zgerm
      real(8)  l_TBASEM, l_TEFFMX

      associate( &
        crop => state%crop%common,    &
        soil => state%soilwater,      &
        mesh => state%mesh,           &
        atmo => state%atmosphere      &
      )

      select case (task)

! === Preparation before crop growth ==========================================
! ADR 0017: case(1) deleted — legacy readarablelandgerm() callsite. The
! dispatch block in InitCropGrowth (lines 90-197) is the modern entry
! point; cache-miss is fatalerr_collected, no legacy fallback.

      case (2)

        ! [GR-CROP 2026-05-25] Dead branch: cropgrowth.f90:255 stub-errors swprep!=0,
        ! so this case never fires in the TOML pipeline. Reads route directly to the
        ! wofost preparation config sub-record for completeness.
        associate(prep => crop_config_global%rotation_wofost(crop%icrop)%preparation)
        node   = 1
        dhPrep = soil%h(node) - prep%hprep
        drz1   = -1.d0 * prep%zprep - mesh%dz(node)
        do while (drz1 .gt. 0.d0)
          node   = node + 1
          dhPrep = max(dhPrep,soil%h(node) - prep%hprep)
          drz1   = drz1 - mesh%dz(node)
        enddo

        crop%flCropPrep = .true.
        if (dhPrep .gt. 0.d0) then
          if (crop%PrepDelay .lt. prep%maxprepdelay) then
            crop%dvs = -0.3d0
            crop%flCropPrep = .false.
            crop%PrepDelay  = crop%PrepDelay + 1
          endif
        endif
        flCropPrep    = crop%flCropPrep    ! dual-write — dispatcher reads bare flCropPrep

        crop%SowDelay = crop%PrepDelay
        end associate
        return

! === Sowing before crop growth ==========================================

      case (3)

        ! [GR-CROP 2026-05-25] Dead branch: cropgrowth.f90:255 stub-errors swsow!=0,
        ! so this case never fires in the TOML pipeline. Reads route directly to the
        ! wofost sowing config sub-record for completeness.
        associate(sow => crop_config_global%rotation_wofost(crop%icrop)%sowing)
        node   = 1
        dhSow  = soil%h(node) - sow%hsow
        drz1   = -1.d0 * sow%zsow - mesh%dz(node)
        do while (drz1 .gt. 0.d0)
          node   = node + 1
          dhSow  = max(dhSow,soil%h(node) - sow%hsow)
          drz1   = drz1 - mesh%dz(node)
        enddo

        node     = 1
        drz1     = -1.d0 * sow%ztempsow - mesh%dz(node)
        do while (drz1 .gt. 0.d0)
          node   = node + 1
          drz1 = drz1 - mesh%dz(node)
        enddo

        ! SS-HEAT pre-Task-8: tsoil read from dummy arg (state%heat%tsoil via caller).
        dtempSow = min(tsoil(node) - sow%tempsow,0.d0)

        crop%flCropSow = .true.
        if (dtempSow .lt. 0.d0 .or. dhSow.gt.0.d0) then
          if (crop%SowDelay .lt. sow%maxsowdelay) then
            crop%dvs = -0.2d0
            crop%flCropSow = .false.
            crop%SowDelay  = crop%SowDelay + 1
          endif
        endif
        flCropSow = crop%flCropSow    ! dual-write — dispatcher reads bare flCropSow
        end associate
        return

! === Simulate germination ==============================================

      case (4)

        associate(germ => crop_config_global%rotation_wofost(crop%icrop)%germination)
        l_agerm      = germ%agerm
        l_tsumemeopt = germ%tsumemeopt
        l_hdrygerm   = germ%hdrygerm
        l_hwetgerm   = germ%hwetgerm
        l_zgerm      = germ%zgerm
        if (l_zgerm == 0.0d0) l_zgerm = -10.0d0   ! legacy default; mirrors cropgrowth.f90:289-291
        l_TBASEM     = germ%tbasem
        l_TEFFMX     = germ%teffmx

        ! Optimal situation in case germination only depends on temperature (swgerm = 1)
        if (germ%swgerm == 1) then
          tsumemesub = l_tsumemeopt
          l_bgerm    = 0.0d0
          l_cgerm    = 0.0d0
        else
          ! swgerm == 2: temperature + hydrological conditions
          l_cgerm = - (l_tsumemeopt - l_agerm * log10(-l_hdrygerm))
          l_bgerm =   (l_tsumemeopt + l_agerm * log10(-l_hwetgerm))

          ! ---   calculate average pressure head of rootzone ---
          if (dabs(l_zgerm-0.d0) .lt. small) then
            hrz1 = soil%h(1)
          else
            node = 0
            hrz1 = 0.0d0
            drz1 = l_zgerm * (-1.d0)
            do while (drz1 .gt. 0.d0)
              node = node + 1
              if (drz1 - mesh%dz(node) .ge. 0.d0) then
                hrz1 = hrz1+soil%h(node)*mesh%dz(node)/(l_zgerm*(-1.d0))
              else
                hrz1 = hrz1+soil%h(node)*drz1/(l_zgerm*(-1.d0))
              endif
              drz1 = drz1 - mesh%dz(node)
            enddo
          endif

          ! --- simulate germination time ---
          pFz1 = DLOG10(MAX(1.0d0,-hrz1))
          if (hrz1 .lt. l_hdrygerm) then
            ! Dry situation
            tsumemesub = l_agerm * pFz1 - l_cgerm
          elseif (hrz1 .ge. l_hdrygerm .and. hrz1 .le. l_hwetgerm) then
            ! Optimal situation
            tsumemesub = l_tsumemeopt
          else
            ! Wet situation
            tsumemesub = -l_agerm * pFz1 + l_bgerm
          endif

        endif

        ! Update of tsumgerm, for the time step of 1 day
        if (atmo%Tav .gt. l_TBASEM)then
          if (atmo%Tav .lt. l_TEFFMX) then
            if(tsumemesub.lt.0.1d0) then
              crop%tsumgerm = crop%tsumgerm + (atmo%Tav-l_TBASEM)
            else
              crop%tsumgerm = crop%tsumgerm +(l_tsumemeopt/tsumemesub)*(atmo%Tav-l_TBASEM)
            endif
          else
            if(tsumemesub.lt.0.1d0) then
              crop%tsumgerm = crop%tsumgerm + (l_TEFFMX-l_TBASEM)
            else
              crop%tsumgerm = crop%tsumgerm +(l_tsumemeopt/tsumemesub)*(l_TEFFMX-l_TBASEM)
            endif
          endif
        endif

        ! Delay growth until tsumgerm is reached
        crop%flCropGerm = .true.
        if (crop%tsumgerm .lt. l_tsumemeopt) then
          crop%dvs = -0.1d0 * max(1.d0 - (crop%tsumgerm / l_tsumemeopt), 0.d0)
          crop%flCropGerm = .false.
        else
          crop%dvs = 0.d0
        endif
        flCropGerm = crop%flCropGerm   ! dual-write — dispatcher reads bare flCropGerm
        end associate

        return

      case default
        call fatalerr_collected ('ArableLandGerm', 'Illegal value for TASK')
      end select

      end associate  ! crop, soil, mesh, atmo
      return
      end subroutine ArableLandGerm

! ----------------------------------------------------------------------
      subroutine FacCO2(state)
! ----------------------------------------------------------------------
!     update             : February 2018
!     date               : ?
!     purpose            : Assimilation correction for CO2 changes in
!                          atmosphere (Lintul4) added by Iwan Supit
! SS-TC TC-10: state added (intent in); iyear read via state%timecontrol
!   tc_iyear alias.
! [GR-CROP Phase B/5] narrow use variables
! ----------------------------------------------------------------------
      ! [SS-GR-CROPRT B6] DEFERRED — FacCO2:
      !   flco2: logical flag (maps to cropwofost_config co2%swco2); no state home
      !   co2year, co2ppm: legacy CO2 table arrays, not yet in config; need migration
      !   co2amaxtb, co2efftb, co2tratb: CO2 correction tables in cropwofost_config co2%;
      !     migration deferred — not yet threaded via state or config arg
      !   mayrs: array dimension (could → swap_array_dimensions, deferred with rest)
      ! MIGRATED B6: fco2amax/fco2eff/fco2tra → written directly to state%crop%wofost%X
      !   (intent changed in→inout). Global write retained for backward compat; CropGrowth
      !   redundant dual-write at lines 395-397 removed (B6 handles it). Reads in CropGrowth
      !   body now use state%crop%wofost%X (removed from CropGrowth use variables).
      use variables, only:                                              & ! [SS-GR-CROPRT B6] flco2 retired (→ state%atmosphere%flco2)
                           co2year, mayrs, co2ppm,                       &
                           co2amaxtb, co2efftb, co2tratb
      use array_utils, only: afgen
      use error_mod, only: fatalerr_collected
      use swap_state_mod, only: swap_state_t
      implicit none

      type(swap_state_t), intent(inout) :: state  ! [SS-GR-CROPRT B6] changed in→inout for fco2 state write

      integer   ifindi,indexyr
      real(8)   CO2
      character(len=200) messag

      ! initialize CO2 impact
      state%crop%wofost%fco2amax = 1.0d0  ! factor to correct AMAX for CO2
      state%crop%wofost%fco2eff  = 1.0d0  ! factor to correct EFF for CO2
      state%crop%wofost%fco2tra  = 1.0d0  ! factor to correct TRA for CO2

      ! correction of CO2 impact
      ! TC-10: iyear read via state%timecontrol%iyear directly (single site, no ASSOCIATE needed).
      if (state%atmosphere%flco2) then
        indexyr = ifindi (CO2year, mayrs, 1, mayrs, state%timecontrol%iyear)  ! TC-10
        if (indexyr.lt.1 .or. indexyr.gt.mayrs) then
          Messag ='Input if CO2year or CO2ppm inconsistent, correct'
          call fatalerr_collected ('wofost',messag)
        endif
        CO2 = CO2ppm(indexyr)
        state%crop%wofost%fco2amax = afgen(CO2AMAXTB,30,CO2)
        state%crop%wofost%fco2eff  = afgen(CO2EFFTB,30,CO2)
        state%crop%wofost%fco2tra  = afgen(CO2TRATB,30,CO2)
      endif

      return
      end subroutine FacCO2

! ----------------------------------------------------------------------
      subroutine update_rootdistribution(state)

! ----------------------------------------------------------------------
!     Date: May 2021
!     Purpose: dynamic root distribution
!              The normalized cumulative root density is modified by
!              growth of root biomass based on relative root water
!              extraction or transpiration reduction (uncompensated).
! ----------------------------------------------------------------------

      ! [SS-SWC S-2.12B] qpotrot_day/qredtot_day retired — read via state%soilwater
      ! SS-TC TC-10: date read via state%timecontrol tc_date alias (removed from variables use).
      ! [SS-GR-CROPRT B9] DEFERRED — update_rootdistribution:
      !   gwrt: root growth rate (computed in wofost task=3), no state home
      !   wrtmin: minimum root weight at relative depth, no state home
      ! MIGRATED B9: noddrz → state%crop%common%noddrz (read-only loop bound)
      ! MIGRATED B9: cumdens → state%crop%common%cumdens (full array in state after init;
      !   odd/even indices all set at CropGrowth task=1 via full array copy line 195)
      ! MIGRATED B9: wrt → state%crop%wofost%wrt (WOFOST root biomass; read prev-day value;
      !   dual-write in wofost keeps state%crop%wofost%wrt current)
      use variables, only: gwrt, wrtmin  ! [SS-GR-CROPRT B9] DEFERRED — no state home
      use swap_state_mod, only: swap_state_t  ! [SS-SWC S-2.12B]
      ! local
      implicit none

      type(swap_state_t), intent(inout) :: state  ! [SS-SWC S-2.12B] inout for cumdens dual-write [SS-GR-CROPRT A5]

      integer   node, i
      real(8)   top,bot
      real(8)   rd_noddrz
      real(8)   rel_qrot_day, rel_qred_day, sum
      real(8)   wrttot, wrtdis(202), qrotdis(202), qreddis(202)
      logical   found

      ! [GR-CROP 2026-05-25] sub-record associate style
      associate( &
        crop => state%crop%common,    &  ! crop runtime (cumdens, noddrz, ...)
        soil => state%soilwater,      &  ! soil-water runtime (qredtot_day, qpotrot_day)
        mesh => state%mesh,           &  ! mesh discretization (zbotcp, ztopcp)
        time => state%timecontrol     &  ! time control (date)
      )

! --- update normalized cumulative root density based on root extraction or stress (cumdens)
!      if (swrdc .eq. 1) then

        ! root extraction of each compartment since start of the day
        rel_qrot_day = 0.d0
        rel_qred_day = 0.d0
        do node = 1,crop%noddrz
          rel_qrot_day = rel_qrot_day + 1 - soil%qredtot_day(node) / soil%qpotrot_day(node)
          rel_qred_day = rel_qred_day + soil%qredtot_day(node)
        enddo

        if ((gwrt .gt. 0.d0 .and. rel_qrot_day .gt. 0.d0) .or. (gwrt .lt. 0.d0 .and. rel_qred_day .gt. 0.d0)) then

          ! distribution roots and root extraction at relative depth
          ! root extraction and root weight based on previous day
          rd_noddrz = abs(mesh%zbotcp(crop%noddrz))
          node = 1
          do i = 4,202,2

            ! root distribution of previous day
            wrtdis(i) = (crop%cumdens(i) - crop%cumdens(i-2)) * (state%crop%wofost%wrt - gwrt)

            ! determine optimal extraction and maximum reduction at relative depth
            found = .false.
            qrotdis(i) = 0.d0
            qreddis(i) = 0.d0
            top = - crop%cumdens(i-3) * rd_noddrz
            bot = - crop%cumdens(i-1) * rd_noddrz
            do while (.not. found)
              if (bot .ge. mesh%zbotcp(node)) then
                qrotdis(i) = qrotdis(i) + (1 - soil%qredtot_day(node) / soil%qpotrot_day(node)) / (mesh%ztopcp(node) - mesh%zbotcp(node)) * (top - bot)
                qreddis(i) = qreddis(i) + soil%qredtot_day(node) / (mesh%ztopcp(node) - mesh%zbotcp(node)) * (top - bot)
                found = .true.
              else
                qrotdis(i) = qrotdis(i) + (1 - soil%qredtot_day(node) / soil%qpotrot_day(node)) / (mesh%ztopcp(node) - mesh%zbotcp(node)) * (top - mesh%zbotcp(node))
                qreddis(i) = qreddis(i) + soil%qredtot_day(node) / (mesh%ztopcp(node) - mesh%zbotcp(node)) * (top - mesh%zbotcp(node))
                top = mesh%zbotcp(node)
                node = node + 1
              end if
            end do
          end do

          ! update relative root weight
          wrttot = 0.d0
          if (gwrt .gt. 0.d0) then
            do i = 4,202,2
              wrtdis(i) = max(wrtmin, wrtdis(i) + (qrotdis(i) / rel_qrot_day) * gwrt)
              wrttot = wrttot + wrtdis(i)
            end do
          elseif (gwrt .lt. 0.d0) then
            do i = 4,202,2
              wrtdis(i) = max(wrtmin, wrtdis(i) + (qreddis(i) / rel_qred_day) * gwrt)
              wrttot = wrttot + wrtdis(i)
            end do
          end if

          ! update normalized cumulative root density distribution
          sum = 0.d0
          do i = 4,202,2
            sum = sum + wrtdis(i)
            crop%cumdens(i) = sum / wrttot
          end do

        end if

        ! TEMPORARY OUTPUT  DELETE
        do i = 2,202,2
          write(777,*) trim(time%date), ",", crop%cumdens(i-1), ",", crop%cumdens(i)
        end do

        do node = 1,crop%noddrz
          write(888,'(a11,",",i4,3(",",f15.5))') trim(time%date), node, soil%qpotrot_day(node), soil%qredtot_day(node)
        end do
        ! TEMPORARY OUTPUT  DELETE

!      end if

      end associate  ! crop, soil, mesh, time
      return
      end subroutine update_rootdistribution

! ----------------------------------------------------------------------
      subroutine sumttd(task,flGrassGrowth,dateGrassGrowth,tsoil,state)
! ----------------------------------------------------------------------
!     Last modified      : Jan 2016
!     Author             : Joop Kroes
!
!     Purpose            : Suppress grass growth as long as 3 criteria are not met:
!                          temperature, time and depth
!
!     Interface parameters, class: I=input,O=output,I/O=input/output
!     class type parameter description (unit)
!       I    C   condition case: 'initial' or 'dynamic' (-)
!       I    R8  tsumtemp  temperature limit to initiate grass growth  [0.0..20.0 grC, R]
!       I    I   tsumtime  time (nrs of sequential days) with temp above tsumtemp for grass growth [1..20 days, I]
!       I    R8  tsumdepth depth at which temp above tsumtemp for grass growth [0.0..100.0 cm below soil surface, R]
!       I    R8  z         depth of a node (L)
!       I    R8  tsoil     Array with soil temperatures (oC) for each compartment
!       O    L   flGrassGrowth flag indicating grass growth (suppressed=.false. when criteria are not met) [.true .or. .false. -, L]
! SS-HEAT pre-Task-8: tsoil now non-optional dummy arg; callers pass state%heat%tsoil
!   via grass's tsoil dummy arg. Global tsoil excluded via rename.
! SS-TC TC-10: state added (intent in); t1900,date read via state%timecontrol tc_* aliases.
! ----------------------------------------------------------------------
      ! [SS-GR-CROPRT B10] DEFERRED — sumttd (stub-path: swtsum=2 is stub-errored):
      !   tsumdepth, tsumtemp, tsumtime: grass growth start thresholds; config fields in
      !     cropgrass_config%tsumdepth/tsumtemp/tsumtime; migration deferred because swtsum=2
      !     path is stub-errored in cropgrass_config.f90 (line 290) — code is unreachable
      !     until swtsum=2 is implemented; no state home needed until then
      !   pathwork, outfil, project: file-path globals; [SS-GR-CROPRT B] DEFERRED file-path arc
      !   tsoil: config-staging buffer, renamed to avoid clash with dummy arg.
      use variables, only: tsumdepth, tsumtemp, tsumtime, &     ! [SS-GR-CROPRT B10] DEFERRED — stub-path (swtsum=2)
                           pathwork, outfil, project,     &     ! [SS-GR-CROPRT B10] DEFERRED — file-path globals
                           dummy_tsoil_sumttd_ => tsoil         ! config-staging buffer
      use file_io_mod, only: file_open
      use swap_state_mod, only: swap_state_t
      implicit none

! --- arguments
      character(len=*), intent(in) :: task
      logical, intent(out)         :: flGrassGrowth      ! flag indicating grass growth (suppressed=.false. when criteria are not met) [.true .or. .false. -, L]
      character(len=11), intent(out) ::  dateGrassGrowth            ! date of start of GrassGrowth
      real(8), intent(in) :: tsoil(:)
      !! Soil temperature array from state%heat%tsoil.
      type(swap_state_t), intent(in) :: state

! --- local
      integer    :: tsumtimecum      ! cumulative, from 1-jan, time (nrs of sequential days) with temp above tsumtemp for grass growth [1..20 days, I]
      logical    :: fltsumtemp       ! flag to indicate if temperature criteria is met
      logical    :: fltsimprev       ! flag to indicate if temperature criterion is met during simulated previous day
      logical    :: fltsimcount      ! flag to indicate nr of contiuous simulated days that temperature criteria is met
      integer    :: cmpcrit, node
!     output
      integer    :: uo   !, idum, ios
      character(len=160)  :: filnam, filtext
      character(len=1)  :: comma
      logical    :: flexist, flopened

      save

      ! [GR-CROP 2026-05-25] sub-record associate style
      associate( &
        mesh => state%mesh,           &  ! mesh discretization (z, numnod)
        time => state%timecontrol     &  ! time control (t1900, date)
      )

      comma = ','

      select case(task)

      case('initial')
          tsumtimecum = 0
          fltsumtemp = .false.
          fltsimprev = .false.
          ! find depth of compartment with critical soil temperature
          cmpcrit = 1
          do node = 1, mesh%numnod
             if (mesh%z(node) .le. (-1.0d0*tsumdepth)) then
                cmpcrit = node
                exit
             endif
          enddo
          ! no growth as long as 3 criteria are not met
          flGrassGrowth = .false.
          dateGrassGrowth = 'undefined'

          ! === open output file and write headers
          filnam = trim(pathwork)//trim(outfil)//'.ttd'
          flopened = .false.
          if(uo.gt.0) then
              !  Inquiry by Unit
              inquire (uo, opened=flopened, exist=flexist)
          endif
          if (.not.flexist .and. .not.flopened) then
              call file_open(uo, filnam, 'replace', 'write')
              filtext = 'output of subr sumttd'
              call writehead (uo,1,filnam,filtext,project)
              write (uo,100)
 100          format (' Date,z(cmpcrit),tsoil(cmpcrit),',               &
     &           'fltsumtemp,fltsimprev,fltsimcount')
          endif

      case('dynamic')
          ! temperature and depth criterium
          ! SS-HEAT pre-Task-8: tsoil read from dummy arg (state%heat%tsoil via caller chain).
          if(tsoil(cmpcrit).ge.tsumtemp) then
              fltsumtemp = .true.
          else
              fltsumtemp = .false.
          endif
          ! timing criterium: set flag for continuous days that exceed critical temperature
          if(fltsumtemp .and. fltsimprev) then
              tsumtimecum = tsumtimecum + 1
          else
              tsumtimecum = 0
              fltsimprev = .false.
          endif
          !  timing criterium: count nr of days exceeding critical temperature
          if(tsumtimecum.ge.tsumtime) then
              fltsimcount = .true.
          else
              fltsimcount = .false.
          endif
          ! no growth as long as 3 criteria are not met
          flGrassGrowth = .false.
          if(fltsumtemp .and. fltsimprev .and. fltsimcount) then
              call dtdpst ('year-month-day',time%t1900,dateGrassGrowth)
              flGrassGrowth = .true.
          endif

          ! timing criteria for next timestep
          if(fltsumtemp) then
              fltsimprev = .true.
          endif

          ! === write output
          write (uo,200) time%date,comma,mesh%z(cmpcrit),comma,tsoil(cmpcrit), &
     &           comma,fltsumtemp,comma,fltsimprev,comma,fltsimcount
 200      format (a11,2(a1,f7.2),3(a1,i3))

      end select

      end associate  ! mesh, time
      return
      end subroutine sumttd

! ----------------------------------------------------------------------
! [SS-BMI2] Crop output buffer helpers (canonical output-sink pattern)
! Buffer lives at top-level swap_state_t (no crop_state_t today).
! N is dynamic — varies by croptype; placeholder N=1 until full crop
! output migration arc populates the row from OutCropFixed/OutWofost/OutGrass.
! ----------------------------------------------------------------------

      subroutine init_crop_output_buffer(state)
! ----------------------------------------------------------------------
!     Allocate state%crop_output_row and set column names.
!     Called from cropoutput(1) — always runs, headless-independent.
!     Currently: N=1 placeholder. Full migration deferred to crop arc.
! ----------------------------------------------------------------------
      use swap_state_mod, only: swap_state_t
      use iso_c_binding,  only: c_double
      implicit none
      type(swap_state_t), intent(inout) :: state
      integer, parameter :: N = 1

      state%crop_output_n_cols = N
      if (.not. allocated(state%crop_output_row))     allocate(state%crop_output_row(N))
      if (.not. allocated(state%crop_output_columns)) allocate(state%crop_output_columns(N))
      state%crop_output_row     = 0.0_c_double
      state%crop_output_columns(1) = 'placeholder'
      end subroutine init_crop_output_buffer


      subroutine build_crop_output_row(state)
! ----------------------------------------------------------------------
!     Fill state%crop_output_row(:) — currently a no-op placeholder.
!     Called from cropoutput(2) — always runs, headless-independent.
!     Full build (OutCropFixed/OutWofost/OutGrass values) deferred to
!     the crop output migration arc.
! ----------------------------------------------------------------------
      use swap_state_mod, only: swap_state_t
      implicit none
      type(swap_state_t), intent(inout) :: state
      ! No-op: buffer populated when crop output migration arc runs.
      end subroutine build_crop_output_row


      subroutine cleanup_crop_output_buffer(state)
! ----------------------------------------------------------------------
!     Deallocate state%crop_output_row and reset counter.
!     Called from cropoutput(3) — always runs, headless-independent.
! ----------------------------------------------------------------------
      use swap_state_mod, only: swap_state_t
      implicit none
      type(swap_state_t), intent(inout) :: state

      if (allocated(state%crop_output_row))     deallocate(state%crop_output_row)
      if (allocated(state%crop_output_columns)) deallocate(state%crop_output_columns)
      state%crop_output_n_cols = 0
      end subroutine cleanup_crop_output_buffer

      end module cropgrowth_helpers_mod
