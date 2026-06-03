! cropgrowth_helpers.f90
! Leaf-level crop helpers extracted from cropgrowth.f90.
! Subroutines: nocrop, ArableLandGerm, FacCO2, update_rootdistribution,
!              sumttd.
!
! [GR-CROP 2026-05-25] use-variables sweep:
!   - ArableLandGerm      → crop_config_global%rotation_wofost(icrop)%X (germ
!                            params; case(2)/(3) prep/sow are dead-branch); state
!                            tsumgerm + state-only writes of flCropPrep/Sow/Germ.
!   - FacCO2              → dead-branch (state%atmosphere%flco2 is dormant);
!                            fco2*amax/eff/tra remain 1.0; fatalerr stub if flco2
!                            is ever re-enabled without restored config wiring.
!   - update_rootdistribution → state%crop%common (cumdens, noddrz);
!                            gwrt/wrtmin still bare globals (writers in
!                            cropgrass_runtime + cropwofost_runtime out of scope).
!   - sumttd              → dead-branch (swtsum=2 stub-guarded in
!                            cropgrass_config.f90:290); body retired to
!                            fatalerr stub.
! ----------------------------------------------------------------------
      module cropgrowth_helpers_mod
      implicit none
      private

      public :: nocrop
      public :: ArableLandGerm
      public :: FacCO2
      public :: update_rootdistribution
      public :: sumttd

      contains

! ----------------------------------------------------------------------
      subroutine nocrop (state)
! ----------------------------------------------------------------------
!     Purpose: zero-out crop runtime fields when no crop is active.
! ----------------------------------------------------------------------
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
!     Last update          : December 2017
!     Purpose              : Drive crop preparation / sowing / germination
!                            up to the emergence event. Reads soil pressure
!                            heads from state%soilwater, soil temperatures
!                            from the tsoil dummy arg, and germination
!                            config from crop_config_global.
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
      !   [GR-CROP 2026-05-25] flCropPrep/flCropSow/flCropGerm fully retired —
      !     state-only via state%crop%common; dispatcher reads from state.
      use swap_constants, only: small
      use error_mod, only: fatalerr_collected
      use swap_state_mod, only: swap_state_t
      implicit none

      integer  task,node
      real(8), intent(in) :: tsoil(:)
      !! Soil temperature array from state%heat%tsoil.
      type(swap_state_t), intent(inout) :: state
      !! State record; state%soilwater%h used for pressure-head reads.
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
        associate(prep => state%crop%rotation_wofost(crop%icrop)%preparation)
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

        crop%SowDelay = crop%PrepDelay
        end associate
        return

! === Sowing before crop growth ==========================================

      case (3)

        ! [GR-CROP 2026-05-25] Dead branch: cropgrowth.f90:255 stub-errors swsow!=0,
        ! so this case never fires in the TOML pipeline. Reads route directly to the
        ! wofost sowing config sub-record for completeness.
        associate(sow => state%crop%rotation_wofost(crop%icrop)%sowing)
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
        end associate
        return

! === Simulate germination ==============================================

      case (4)

        associate(germ => state%crop%rotation_wofost(crop%icrop)%germination)
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
!     Last update          : February 2018
!     Purpose              : Assimilation correction for CO2 changes in
!                            atmosphere (Lintul4, Iwan Supit).
! ----------------------------------------------------------------------
      ! [GR-CROP 2026-05-25] FacCO2 dead-branch cluster (Class E).
      !   state%atmosphere%flco2 is hard-coded .false. (atmosphere_state.f90:128 default,
      !   no TOML wiring). The `if (flco2)` block is unreachable, so the entire CO2
      !   correction lookup (co2year/co2ppm/mayrs/co2amaxtb/co2efftb/co2tratb) is dead.
      !   fco2amax/fco2eff/fco2tra remain 1.0 (no-op multipliers). When flco2 wiring is
      !   restored, point the lookup at
      !   crop_config_global%rotation_wofost(icrop)%co2 (co2amaxtb/co2efftb/co2tratb
      !   already exist on the wofost_co2_t config sub-record).
      use swap_state_mod, only: swap_state_t
      use error_mod, only: fatalerr_collected
      implicit none

      type(swap_state_t), intent(inout) :: state

      ! initialize CO2 impact (no-op while flco2=.false.)
      state%crop%wofost%fco2amax = 1.0d0  ! factor to correct AMAX for CO2
      state%crop%wofost%fco2eff  = 1.0d0  ! factor to correct EFF for CO2
      state%crop%wofost%fco2tra  = 1.0d0  ! factor to correct TRA for CO2

      ! [GR-CROP 2026-05-25] flco2 is dormant — stub-error if it ever gets enabled
      ! without rewiring the lookup against state%cfg%crop%wofost%co2 tables.
      if (state%atmosphere%flco2) then
        call fatalerr_collected('FacCO2', &
          'flco2=.true. encountered but CO2-correction lookup retired in 2026-05-25 ' // &
          'arc (co2year/co2ppm/co2*tb bare globals removed). Re-wire against ' // &
          'crop_config_global%rotation_wofost(icrop)%co2 before re-enabling.')
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

      ! [GR-CROP 2026-05-25] update_rootdistribution reads only.
      !   gwrt/wrtmin: read from state%crop%wofost%gwrt and state%crop%wofost%wrtmin.
      !     The active rotation's cropX_runtime task=3 populates them; the dispatcher
      !     selects a single rotation per step so no cross-rotation contention.
      !   noddrz/cumdens/wrt: already on state%crop%common / state%crop%wofost.
      use swap_state_mod, only: swap_state_t
      implicit none

      type(swap_state_t), intent(inout) :: state  ! inout: cumdens written via state alias

      integer   node, i
      real(8)   top,bot
      real(8)   rd_noddrz
      real(8)   rel_qrot_day, rel_qred_day, sum
      real(8)   wrttot, wrtdis(202), qrotdis(202), qreddis(202)
      logical   found

      ! [GR-CROP 2026-05-25] sub-record associate style
      associate( &
        crop => state%crop%common,    &  ! crop runtime (cumdens, noddrz, ...)
        wof  => state%crop%wofost,    &  ! crop wofost runtime (wrt, gwrt, wrtmin)
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

        if ((wof%gwrt .gt. 0.d0 .and. rel_qrot_day .gt. 0.d0) .or. (wof%gwrt .lt. 0.d0 .and. rel_qred_day .gt. 0.d0)) then

          ! distribution roots and root extraction at relative depth
          ! root extraction and root weight based on previous day
          rd_noddrz = abs(mesh%zbotcp(crop%noddrz))
          node = 1
          do i = 4,202,2

            ! root distribution of previous day
            wrtdis(i) = (crop%cumdens(i) - crop%cumdens(i-2)) * (wof%wrt - wof%gwrt)

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
          if (wof%gwrt .gt. 0.d0) then
            do i = 4,202,2
              wrtdis(i) = max(wof%wrtmin, wrtdis(i) + (qrotdis(i) / rel_qrot_day) * wof%gwrt)
              wrttot = wrttot + wrtdis(i)
            end do
          elseif (wof%gwrt .lt. 0.d0) then
            do i = 4,202,2
              wrtdis(i) = max(wof%wrtmin, wrtdis(i) + (qreddis(i) / rel_qred_day) * wof%gwrt)
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
! ----------------------------------------------------------------------
      ! [GR-CROP 2026-05-25] sumttd is dead-branch (Class E).
      !   Called only from cropgrass_runtime.f90:299 when state%crop%grass%swtsum==2.
      !   The cropgrass config validator (cropgrass_config.f90:290) stub-errors
      !   swtsum=2, so this routine is unreachable in the TOML pipeline. The body
      !   (which read tsumdepth/tsumtemp/tsumtime + pathwork/outfil/project legacy
      !   globals plus a tsoil staging-buffer rename) has been replaced with a
      !   fatalerr_collected stub. When swtsum=2 is restored, port the body's
      !   threshold reads to cropgrass_config%tsumtemp/tsumtime/tsumdepth and the
      !   I/O reads to state%cfg%general%X.
      use error_mod, only: fatalerr_collected
      use swap_state_mod, only: swap_state_t
      implicit none

! --- arguments
      character(len=*), intent(in) :: task
      logical, intent(out)         :: flGrassGrowth      ! flag indicating grass growth (suppressed=.false. when criteria are not met) [.true .or. .false. -, L]
      character(len=11), intent(out) ::  dateGrassGrowth            ! date of start of GrassGrowth
      real(8), intent(in) :: tsoil(:)
      !! Soil temperature array from state%heat%tsoil.
      type(swap_state_t), intent(in) :: state

      ! Suppress "unused dummy argument" warnings.
      if (.false.) then
         flGrassGrowth   = .false.
         dateGrassGrowth = ''
         if (size(tsoil) > 0 .or. state%crop%common%swcrp >= 0) continue
      end if

      call fatalerr_collected('sumttd', &
        'sumttd is dead-branch (swtsum=2 stub-guarded in cropgrass_config.f90:290). ' // &
        'Body retired in 2026-05-25 arc; restore via cropgrass_config thresholds + ' // &
        'config%general I/O when swtsum=2 wiring is implemented.')

      ! Defensive defaults (never reached due to fatalerr above).
      flGrassGrowth = .false.
      dateGrassGrowth = 'undefined'
      if (task == 'initial' .or. task == 'dynamic') return
      return
      end subroutine sumttd

      end module cropgrowth_helpers_mod
