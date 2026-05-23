! File VersionID:
!   $Id: rootextraction.f90 374 2018-03-21 13:12:23Z heine003 $
! ----------------------------------------------------------------------
!> Root water extraction routines for SWAP crop module.
!!
!! Groups the legacy root extraction procedures in a module so callers use
!! explicit interfaces.
module rootextraction_mod
   use error_mod, only: fatalerr_collected
   use swap_state_mod, only: swap_state_t
  implicit none
  private
  public :: RootExtraction, MatricFlux

  contains

!> Calculate root water extraction profile and stress partitioning.
!!
!! Computes root water extraction for each rooted node, including reductions due
!! to oxygen, drought, salinity, and frost stress, and optional compensation.
!!
!! @note
!! update: August 2016: microscopic uptake according to JongvanLier(2013)
!! update: August 2012: O2-stress according to Bartholomeus(2008)
!! update: February 2011: macrosopic uptake extended with compensation
!!         according to Jarvis (1989)
!! date: August 2004
!! purpose: Calculate the root water extraction rate as function of soil water
!!          pressure head and salinity concentration for each node
!! @endnote
  subroutine RootExtraction(state)
! ----------------------------------------------------------------------
!     update    : August 2016: microscopic uptake according to JongvanLier(2013)
!     update    : August 2012: O2-stress according to Bartholomeus(2008)
!     update    : February 2011: macrosopic uptake extended with
!                                compensation according to Jarvis (1989)
!     date      : August 2004
!     purpose   : Calculate the root water extraction rate as function of soil
!                 water pressure head and salinity concentration for each node
! ----------------------------------------------------------------------
      ! [SS-GR-FINAL B6] macp → swap_array_dimensions; remainder DEFERRED
      use swap_array_dimensions, only: macp
      use swap_log, only: log_warn
      use variables, only: &                                               ! [SS-GR-FINAL B6] residuals — all DEFERRED
                           ! DEFERRED: adcrh/adcrl/aeratecrit/alphacrit — O2/drought config; Phase C3
                           ! adcrh/adcrl/aeratecrit/alphacrit retired
                           ! DEFERRED: botcom/criterhr/cumdens/dcritrtz/flhydrlift — soil/crop config; Phase C3
                           botcom, criterhr, flhydrlift,                   &  ! dcritrtz/cumdens retired
                           ! DEFERRED: hlim1/hlim2l/hlim2u/hlim3h/hlim3l/hlim4 — drought stress limits; Phase C3
                           ! hlim1/hlim2l/hlim2u/hlim3h/hlim3l/hlim4 retired
                           ! DEFERRED: kroot/kstem/noddrz/oxygenintercept — crop/log globals; Phase C3
                           kroot, kstem, noddrz, oxygenintercept,   &
                           ! DEFERRED: oxygenslope/rdctb/rootcoefa/rooteff — active crop state; Phase C3; rd/rdm retired
                           oxygenslope, rootcoefa, rooteff, &  ! rdctb retired
                           ! DEFERRED: rootradius/rxylem/saltmax/saltslope/stephr — crop/solute config; Phase C3
                           rootradius, rxylem, stephr,                    &  ! saltmax/saltslope retired
                           ! DEFERRED: swcompensate/swdrought/swfrost/swoxygen — crop stress switches; Phase C3
                           swfrost,                                       &  ! swcompensate/swoxygen/swdrought retired
                           ! DEFERRED: swoxygentype/swsalinity/swstressor/swwrtnonox — crop stress switches; Phase C3
                           ! swsalinity/swwrtnonox/swstressor/swoxygentype retired
                           ! DEFERRED: taccur/twilt/wiltpoint — soil convergence/stress params; Phase C3
                           twilt, wiltpoint
      use array_utils, only: afgen
      use oxygenstress_mod, only: OxygenStress, OxygenReproFunction
      implicit none

      type(swap_state_t), intent(inout) :: state

! --- local variables
      integer node
      real(8) top,bot,hlim3,hlim2,qred
      real(8) alpdry,alpwet,alpsol,alpfrs,alptot,vsmall
      real(8) alpdrycom,alpwetcom,alpsolcom,alpfrscom,alptotcom
      real(8) rd_noddrz, redtot

      parameter (vsmall = 1.0d-14)

! --- SS-CRP Phase 2 C-2.5: reset writes — legacy globals dropped; state-only.
      state%soilwater%qrot(1:state%mesh%numnod) = 0.0d0   ! [GR-BH Task 35] numnod global deleted
      state%soilwater%Tactual = state%soilwater%qrosum
      state%soilwater%qrosum = 0.0d0
      state%soilwater%qreddrysum = 0.0d0
      state%soilwater%qredwetsum = 0.0d0
      state%soilwater%qredsolsum = 0.0d0
      state%soilwater%qredfrssum = 0.0d0

! --- skip routine if the crop has not emerged
!      if (.not.flcropEmergence) return

! --- skip routine if there are no roots
      if (state%crop%common%rd .lt. vsmall) return

! --- skip routine if transpiration rate is zero
      ! SS-ATM Phase 2 Task A-2.3: ptra read from state%atmosphere (atmosphere home).
      if (state%atmosphere%ptra .lt. 1.d-10) return

! --- SS-CRP Phase 2 C-2.5: ASSOCIATE for the main compute body (state-only).
! --- SS-ATM Phase 2 Task A-2.3: at_ptra, at_atmdem added for atmosphere reader cutover.
! --- SS-SWC S-2.7: sw_h, sw_theta added for soil-water-core reader cutover.
      associate( &
         at_ptra       => state%atmosphere%ptra,       &
         at_atmdem     => state%atmosphere%atmdem,      &
         sw_h          => state%soilwater%h,           &   ! [SS-SWC S-2.7]
         sw_theta      => state%soilwater%theta,       &   ! [SS-SWC S-2.7]
         cw_qrot       => state%soilwater%qrot,        &
         cw_qpotrot    => state%soilwater%qpotrot,     &
         cw_qredwet    => state%soilwater%qredwet,     &
         cw_qreddry    => state%soilwater%qreddry,     &
         cw_qredsol    => state%soilwater%qredsol,     &
         cw_qredfrs    => state%soilwater%qredfrs,     &
         cw_qrosum     => state%soilwater%qrosum,      &
         cw_qredwetsum => state%soilwater%qredwetsum,  &
         cw_qreddrysum => state%soilwater%qreddrysum,  &
         cw_qredsolsum => state%soilwater%qredsolsum,  &
         cw_qredfrssum => state%soilwater%qredfrssum,  &
         cw_flWrtNonox => state%soilwater%flWrtNonox   &
      )

! --- DROUGHT REDUCTION ACCORDING TO FEDDES ET AL. (1978)
      if (state%crop%common%swdrought .eq. 1) then

! --- calculate potential root extraction of the compartments
        ! 22-10-2018: bug repair signalled by Paul van Walsum: division not by rd but by depth bottom of last compartment where roots are present
        !             rd replaced by (newly calculated) rd_noddrz
        rd_noddrz = abs(state%mesh%zbotcp(noddrz))  ! [GR-BH C7]
        do node = 1,noddrz
          top = abs(state%mesh%ztopcp(node) / rd_noddrz)  ! [GR-BH C7]
          bot = abs(state%mesh%zbotcp(node) / rd_noddrz)  ! [GR-BH C7]
          cw_qrot(node) = (afgen(state%crop%common%cumdens,202,bot)-afgen(state%crop%common%cumdens,202,top))* at_ptra
        enddo

! --- calculating critical point hlim3 according to feddes
        if (at_atmdem .lt. state%crop%common%adcrl) then
          hlim3 = state%crop%common%hlim3l
        elseif (at_atmdem .le. state%crop%common%adcrh) then
          hlim3 = state%crop%common%hlim3h + ((state%crop%common%adcrh - at_atmdem) / (state%crop%common%adcrh - state%crop%common%adcrl)) * (state%crop%common%hlim3l - state%crop%common%hlim3h)
        else
          hlim3 = state%crop%common%hlim3h
        endif
      endif

! --- DROUGHT REDUCTION ACCORDING TO DE JONG VAN LIER ET AL. (2012)
      if (state%crop%common%swdrought .eq. 2) then
        call JongvanLier(state)
      endif

! === COMBINATION OF OXYGEN, DROUGHT, SALT AND FROST STRESS ====

      cw_qrosum = 0.0d0

      do 200 node = 1,noddrz
        alpdry = 1.0d0
        alpwet = 1.0d0
        alpsol = 1.0d0
        alpfrs = 1.0d0

! ---   reduction due to oxygen stress
        if (state%crop%common%swoxygen .ne. 0) then

          ! Feddes linear reduction based on pressure head
          if (state%crop%common%swoxygen .eq. 1) then

            if (node.gt.botcom(1)) then
              hlim2 = state%crop%common%hlim2l
            else
              hlim2 = state%crop%common%hlim2u
            endif
            if (sw_h(node).le.state%crop%common%hlim1.and.sw_h(node).gt.hlim2) then   ! [SS-SWC S-2.7]
              alpwet = (state%crop%common%hlim1-sw_h(node))/(state%crop%common%hlim1-hlim2)   ! [SS-SWC S-2.7]
            endif
            if (sw_h(node).gt.state%crop%common%hlim1) then           ! [SS-SWC S-2.7]
              alpwet = 0.0d0
            endif

          ! Bartholomeus non-linear reduction based on gas filled porosity
          elseif (state%crop%common%swoxygen .eq. 2) then

            ! use physical processes
            if (state%crop%common%swoxygentype .eq. 1) then
!##MH         call OxygenStress(node,alpwet,ResultsOxygenStress)
              ! SS-HEAT Phase 2 Task 6: pass state so OxygenStress reads tsoil from state%heat
              call OxygenStress(node,alpwet,state)

            ! use reproduction functions
            else
              ! SS-HEAT Phase 2 Task 6: pass tsoil from state%heat
              call OxygenReproFunction (OxygenSlope,OxygenIntercept,sw_theta,state%soilwater%thetas,state%heat%tsoil,node,state%mesh%z,state%mesh%dz,alpwet,state)  ! [SS-SWC S-2.7] [GR-BH C7] [GR-BH Task 35]
            endif

          endif

! ---     Stop root development in case of oxgenstress at noddrz
!         WOFOST: root zone remain the aim, but biomass is increasing
!         GRASS : stop root development
          cw_flWrtNonox = .false.
          if (state%crop%common%swWrtNonox .eq. 1 .and. node .eq. noddrz) then
            if (alpwet .lt. state%crop%common%aeratecrit) then
              cw_flWrtNonox = .true.
            end if
          endif

        endif

! ---   reduction due to drought stress

        ! Feddes linear reduction based on pressure head
        if (state%crop%common%swdrought .eq. 1) then
          if (sw_h(node) .lt. state%crop%common%hlim4) then         ! [SS-SWC S-2.7]
            alpdry = 0.0d0
          elseif (sw_h(node).le.hlim3) then                         ! [SS-SWC S-2.7]
            alpdry = (state%crop%common%hlim4-sw_h(node))/(state%crop%common%hlim4-hlim3)  ! [SS-SWC S-2.7]
          endif
        endif

        ! JongvanLier microscopic concept for drought
        if (state%crop%common%swdrought .eq. 2) then
          alpdry = state%soilwater%alpJvLier
        endif

! ---   reduction due to salt stress

        ! reduction according to Maas and Hoffman linear reduction function
        if (state%crop%common%swsalinity .eq. 1) then
          if (state%solute%cml(node) .gt. state%crop%common%saltmax) then
            alpsol = 1.0d0 - (state%solute%cml(node) - state%crop%common%saltmax) * state%crop%common%saltslope
            alpsol = max(0.0d0,alpsol)
          endif
        endif
! ---   mind: in case of salt stress with osmotic head, microscopic root water extraction
! ---         according to JongvanLier (2013) should be used (swsalinity = 2); in that case
! ---         salinity stress is included in drought stress and not separately specified
! ---         in output file *.STR

! ----  reduction due to frost conditions
        ! SS-HEAT Phase 2 Task 6: read tsoil from state%heat
        if (swfrost .eq.1 .and. state%heat%tsoil(node) .lt. 0.0d0) then
          alpfrs = 0.0d0
        endif

! ----  overall reduction
        cw_qpotrot(node) = cw_qrot(node)
        cw_qrot(node) = cw_qrot(node) * alpwet * alpdry * alpsol * alpfrs
        cw_qrosum = cw_qrot(node) + cw_qrosum

! ----  apportionment to different types stresses (cm)

        qred = cw_qpotrot(node) - cw_qrot(node)
        if (qred .lt. vsmall)then

          ! no stress
          cw_qredwet(node) = 0.d0
          cw_qreddry(node) = 0.d0
          cw_qredsol(node) = 0.d0
          cw_qredfrs(node) = 0.d0

        else

          ! multiplication of stressors
          alptot = (1 - alpwet) + (1 - alpdry) + (1 - alpsol) + (1 - alpfrs)

          ! contribution of each stressor (lineair approach)
          cw_qredwet(node) = (1 - alpwet) / alptot * qred
          cw_qreddry(node) = (1 - alpdry) / alptot * qred
          cw_qredsol(node) = (1 - alpsol) / alptot * qred
          cw_qredfrs(node) = (1 - alpfrs) / alptot * qred

          ! sum of each stressor (rootzone)
          cw_qredwetsum = cw_qredwetsum + cw_qredwet(node)
          cw_qreddrysum = cw_qreddrysum + cw_qreddry(node)
          cw_qredsolsum = cw_qredsolsum + cw_qredsol(node)
          cw_qredfrssum = cw_qredfrssum + cw_qredfrs(node)

        end if

200   continue

! --- compensated root water uptake according to Jarvis (1989) or Walsum (2020)
      if (state%crop%common%swcompensate .gt. 0) then

        ! compensated root water uptake according to Walsum
        if (state%crop%common%swcompensate .eq. 2) then
            state%crop%common%alphacrit = min((state%crop%common%dcritrtz + state%crop%common%rdm - rd_noddrz) / state%crop%common%rdm, 1.0d0)
        end if

        alptot = cw_qrosum / at_ptra
        qred = at_ptra - cw_qrosum
        if (abs(state%crop%common%alphacrit - 1.0d0) .ge. vsmall .and. qred .gt. vsmall .and. alptot .ge. 0.05d0) then
          ! Only compensation when rootextraction and transpiration reduction is greater than vsmall
          ! and when alptot > 0.05, i.e. when there is less than 95% stress reduction. This minimum is
          ! also important for the approximation of alp... in the next 4 lines.
          alpdry = alptot**(cw_qreddrysum/qred)
          alpwet = alptot**(cw_qredwetsum/qred)
          alpsol = alptot**(cw_qredsolsum/qred)
          alpfrs = alptot**(cw_qredfrssum/qred)

          if (state%crop%common%swstressor .eq. 1) then
            alptotcom = min(alptot / state%crop%common%alphacrit, 1.d0)
            alpdrycom = alpdry
            alpwetcom = alpwet
            alpsolcom = alpsol
            alpfrscom = alpfrs
          else
            alpdrycom = alpdry
            alpwetcom = alpwet
            alpsolcom = alpsol
            alpfrscom = alpfrs
            if (state%crop%common%swstressor .eq. 2) then
              alpdrycom = min(alpdry / state%crop%common%alphacrit, 1.d0)
            elseif (state%crop%common%swstressor .eq. 3) then
              alpwetcom = min(alpwet / state%crop%common%alphacrit, 1.d0)
            elseif (state%crop%common%swstressor .eq. 4) then
              alpsolcom = min(alpsol / state%crop%common%alphacrit, 1.d0)
            elseif (state%crop%common%swstressor .eq. 5) then
              alpfrscom = min(alpfrs / state%crop%common%alphacrit, 1.d0)
            endif
            alptotcom = alpwetcom * alpdrycom * alpsolcom * alpfrscom
          endif

          ! Change the abstraction of the roots
          do node = 1,noddrz
            cw_qrot(node) = cw_qrot(node) * alptotcom / alptot
          enddo

          ! Change the sum-parameters
          cw_qrosum = at_ptra * alptotcom
          qred = at_ptra - cw_qrosum
          if (qred .lt. vsmall) then
            ! There is no stress.
            cw_qredwetsum = 0.0d0
            cw_qreddrysum = 0.0d0
            cw_qredsolsum = 0.0d0
            cw_qredfrssum = 0.0d0
          else
            redtot = (1 - alpwetcom) + (1 - alpdrycom) + (1 - alpsolcom) + (1 - alpfrscom)
            cw_qredwetsum = (1 - alpwetcom) / redtot * qred
            cw_qreddrysum = (1 - alpdrycom) / redtot * qred
            cw_qredsolsum = (1 - alpsolcom) / redtot * qred
            cw_qredfrssum = (1 - alpfrscom) / redtot * qred
          endif
        endif

      endif

      end associate
      return

        end subroutine RootExtraction

! ----------------------------------------------------------------------
      !> Calculate microscopic root extraction after De Jong van Lier et al. (2013).
      !!
      !! @note
      !! date: August 2016
      !! purpose: Calculate the root water extraction rate according to
      !!          De Jong van Lier et al. (2013)
      !! @endnote
        subroutine JongvanLier(state)
! ----------------------------------------------------------------------
!     date      : August 2016
!     purpose   : Calculate the root water extraction rate according to
!                 De Jong van Lier et al. (2013)
! [GR-CROP Phase B/9] narrow use variables
! ----------------------------------------------------------------------
      ! [SS-GR-FINAL B6] macp → swap_array_dimensions; remainder DEFERRED (same as RootExtraction)
      use swap_array_dimensions, only: macp
      use swap_log, only: log_warn
      use variables, only: &                                               ! [SS-GR-FINAL B6] residuals — all DEFERRED
                           ! adcrh/adcrl/aeratecrit/alphacrit retired  ! DEFERRED: crop stress config; Phase C3
                           botcom, criterhr, flhydrlift,                   &  ! dcritrtz/cumdens retired  ! DEFERRED: soil/crop config
                           ! hlim1/hlim2l/hlim2u/hlim3h/hlim3l/hlim4 retired  ! DEFERRED: drought limits
                           kroot, kstem, noddrz, oxygenintercept,   &  ! DEFERRED: crop/log globals
                           oxygenslope, rootcoefa, rooteff, &  ! rdctb retired  ! rd retired  ! DEFERRED: active crop state
                           rootradius, rxylem, stephr,                    &  ! saltmax/saltslope retired  ! DEFERRED: crop/solute config
                           swfrost,                                       &  ! swcompensate/swoxygen/swdrought retired  ! DEFERRED: stress switches
                           ! swsalinity/swwrtnonox/swstressor/swoxygentype retired  ! DEFERRED: stress switches
                           twilt, wiltpoint                          ! DEFERRED: convergence/stress params
      use array_utils, only: afgen
      implicit none

      type(swap_state_t), intent(inout) :: state

! --- local variables
      integer node,counter
      real(8) Fy3
      real(8) y1,y2,y3,Fy1,Fy2,hwet,ratio
      real(8) hleafm1,hrootm1(macp)
      real(8) reldepth,rdensity,phi,meandepth
      real(8) maxstep,dHPLant,dHPlantMax
      logical flconverg,flstress
      character(len=200) messag

! --- SS-CRP Phase 2 C-2.5: ASSOCIATE for JvL per-node arrays and scalars (state-only).
! --- SS-ATM Phase 2 Task A-2.3: at_ptra added for atmosphere reader cutover.
! --- SS-SWC S-2.7: sw_h, sw_theta added for soil-water-core reader cutover.
      associate( &
         at_ptra     => state%atmosphere%ptra,    &
         sw_h        => state%soilwater%h,        &   ! [SS-SWC S-2.7]
         sw_theta    => state%soilwater%theta,    &   ! [SS-SWC S-2.7]
         cw_mflux    => state%soilwater%mflux,    &
         cw_mroot    => state%soilwater%mroot,    &
         cw_hroot    => state%soilwater%hroot,    &
         cw_rootrho  => state%soilwater%rootrho,  &
         cw_rootphi  => state%soilwater%rootphi,  &
         cw_rmax     => state%soilwater%rmax,     &
         cw_hleaf    => state%soilwater%hleaf,    &
         cw_Hxylem   => state%soilwater%Hxylem,   &
         cw_qrosum   => state%soilwater%qrosum,   &
         cw_qrot     => state%soilwater%qrot,     &
         cw_alpJvLier => state%soilwater%alpJvLier, &
         swscre       => state%timecontrol%swscre  &  ! [SS-BMI2 Task 4]
      )

! --- initialization
      phi = 3.1415926d0
      counter = 0
      flstress = .false.
      hleafm1 = cw_hleaf
      flconverg = .false.

! --- reset values below root zone to zero
      do node = noddrz+1,state%mesh%numnod  ! [GR-BH C7]
        cw_mflux(node) = 0.d0
        cw_mroot(node) = 0.d0
        cw_hroot(node) = 0.d0
        cw_rootrho(node) = 0.d0
        cw_rootphi(node) = 0.d0
      enddo

! --- give hroot a value when root zone becomes larger and store previous values
      do node = 1, noddrz
        if (abs(cw_hroot(node)) .lt. 1.0d-12) then
          cw_hroot(node) = sw_h(node)            ! [SS-SWC S-2.7]
        end if
        hrootm1(node) = cw_hroot(node)
      enddo

! --- initialization of rootrho and rootphi
      do node = 1,noddrz-1
        reldepth = -state%mesh%z(node)/state%crop%common%rd  ! [GR-BH C7]
        rdensity = afgen(state%crop%common%rdctb,22,reldepth)
        cw_rmax(node) = 1.d0/dsqrt(phi*rdensity)
        cw_rootrho(node) = 4.d0/(rootradius*rootradius-rootcoefa*cw_rmax(node)&
     &      *rootcoefa*cw_rmax(node) + 2.d0 * (cw_rmax(node)*cw_rmax(node) +     &
     &      rootradius*rootradius)*log(rootcoefa*cw_rmax(node)/rootradius))
        cw_rootphi(node) = cw_rootrho(node) * cw_rmax(node)*cw_rmax(node) *         &
     &     log(rootradius/rxylem) * 0.5d0 / kroot
      enddo
! --- last node, partly filled with roots
      node = noddrz
      meandepth = (state%mesh%ztopcp(node)-state%crop%common%rd)*0.5d0  ! [GR-BH C7]
      reldepth = meandepth/(-state%crop%common%rd)
      rdensity = afgen(state%crop%common%rdctb,22,reldepth)
      cw_rmax(node) = 1.d0/dsqrt(phi*rdensity)
      cw_rootrho(node) = 4.d0/(rootradius*rootradius-rootcoefa*cw_rmax(node)  &
     &      *rootcoefa*cw_rmax(node) + 2.d0 * (cw_rmax(node)*cw_rmax(node) +     &
     &      rootradius*rootradius)*log(rootcoefa*cw_rmax(node)/rootradius))
      cw_rootphi(node) = cw_rootrho(node) * cw_rmax(node)*cw_rmax(node) *           &
     &     log(rootradius/rxylem) * 0.5d0 / kroot

! --- calculate current matric flux potential in soil water
      do node = 1,noddrz
        sw_h(node) = min(1.0d3,max(sw_h(node),1.0d-8))  ! [SS-SWC S-2.7]
        call MatricFlux(2,sw_h(node),node,cw_mflux(node),state)  ! [SS-SWC S-2.7]
      enddo

! --- interpolate matric flux potential in lowest compartment that is partly filled with roots
      node = noddrz
      if (node.gt.1) then
         cw_mflux(node) = cw_mflux(node) + (cw_mflux(node-1)-cw_mflux(node))/       &
     &                 state%mesh%disnod(node) * (meandepth-state%mesh%z(node))  ! [GR-BH C7]
      endif

! --- determine highest pressure head in root zone
      hwet = sw_h(1)                                     ! [SS-SWC S-2.7]
      do node = 2,noddrz
        if (sw_h(node) .gt. hwet) hwet = sw_h(node)     ! [SS-SWC S-2.7]
      enddo

! --- determine whether qrosum < ptra (flstress = true)
! --- calculate root water extraction qrosum at Hleaf = wiltpoint and Tactual = Ptra
      cw_hleaf = wiltpoint
      dHPlant = at_ptra/Kstem
      if (dHPlant .gt. (hwet - wiltpoint)) flstress = .true.

      cw_Hxylem = cw_hleaf + dHPLant
      call JongvanLierLoop(state)
      if (cw_qrosum .lt. at_ptra) flstress = .true.

! --- set hroot to previous values
      do node = 1,noddrz
        cw_hroot(node) = hrootm1(node)
      enddo

! --- FLSTRESS = FALSE: DETERMINE HLEAF WITH TACTUAL = PTRA
      if (.not. flstress) then

! --- calculate root extraction at previous value of hleaf
      cw_hleaf = Hleafm1
      cw_Hxylem = cw_hleaf + dHPLant
      call JongvanLierLoop(state)

! --- check convergence
      if (abs(cw_qrosum-at_ptra) .lt. state%cfg%simulation%numerical%taccur) flconverg = .true.

! --- set initial values y3 and Fy3
      y3 = cw_hleaf
      Fy3 = at_ptra - cw_qrosum

! --- start search on correct value of pressure head in leaves
      do while (.not. flconverg)

! --- reset y1 and Fy1
      y1 = y3
      Fy1 = Fy3

      if (cw_qrosum .lt. at_ptra) then
! ---   root water extraction too small, decrease Hleaf
          y2 = y1 - StepHr*log10(max(abs(y1),1.d0))
      else
! ---   root water extraction too large, increase Hleaf
          y2 = y1 + StepHr*log10(max(abs(y1),1.d0))
      endif

! --- calculate root extraction at modified pressure head y2 in leaf
      cw_hleaf = y2
      cw_Hxylem = cw_hleaf + dHPLant
      call JongvanLierLoop(state)
      Fy2 = at_ptra - cw_qrosum

! --- determine new value of leaf pressure head with Newton Raphson algorithm
      if (abs(Fy1-Fy2).lt. 1.d-10) then
! ---   take maximum step
        if (cw_qrosum .gt. at_ptra) then
          maxstep = -0.05d0 * y1
          maxstep = max(100.d0,maxstep)
          y3 = y1 + maxstep
        else
          maxstep = 0.05d0 * (y1 - wiltpoint)
          maxstep = max(100.d0,maxstep)
          y3 = y1 - maxstep
        endif
      else
! ---   take Newton Raphson step
        y3 = y1 - (y2-y1)*Fy1 / (Fy2-Fy1)
        if ((y3-y1) .gt. 0.d0) then
          maxstep = -0.05d0 * y1
          maxstep = max(100.d0,maxstep)
          y3 = min((y1 + maxstep),y3)
        else
          maxstep = 0.05d0 * (y1 - wiltpoint)
          maxstep = max(100.d0,maxstep)
          y3 = max((y1 - maxstep),y3)
        endif
      endif
      y3 = min(y3,hwet)
      y3 = max(y3,wiltpoint)

! --- new estimated value of leaf pressure head is equal to y3!
 300  cw_hleaf = y3
      cw_Hxylem = cw_hleaf + dHPlant

! --- calculate root extraction at y3
      call JongvanLierLoop(state)
      Fy3 = at_ptra - cw_qrosum

! --- check convergence
      if (abs(cw_qrosum-at_ptra) .lt. state%cfg%simulation%numerical%taccur .or.                             &
     &       abs(y3-y1) .lt. 1.0d0) flconverg = .true.

! --- fatal error if too many iterations
      counter = counter + 1
      if (counter .gt. 1000) then
         messag = '4 Too many iterations for microscopic root'          &
     &             //' water uptake. Please adapt input!'
         call log_warn('rootextraction', messag)
         call fatalerr_collected ('rootextraction',messag)
      endif

! --- apply linear interpolation when Fy1 and Fy3 have opposite sign
      if (.not. flconverg) then
        if (((Fy1 .gt. 0.d0 .and. Fy3 .lt. 0.d0) .or.                   &
     &                     (Fy1 .lt. 0.d0 .and. Fy3 .gt. 0.d0)) .and.   &
     &                      abs(y1-y3) .gt. 1.d0) then

          y3 = y1 + (y3 - y1) * Fy1 / (Fy1 - Fy3)
          y3 = min(y3,hwet)
          y3 = max(y3,wiltpoint)
          goto 300
        endif
      endif

! --- WHILE-DO LOOP TO DETERMINE HLEAF WITHOUT DROUGHT STRESS
      enddo

      cw_qrosum = at_ptra
      cw_hleaf = y3
      dHplant = cw_qrosum/kstem
      cw_Hxylem = cw_hleaf + dHPlant

      else
! --- FLSTRESS = TRUE: DETERMINE QROSUM WITH HLEAF = WILTPOINT
      counter = 0

! --- calculate root extraction at previous value of qrosum
      cw_hleaf = wiltpoint
      dHplant = state%soilwater%Tactual/kstem
      if (dHplant .gt. (hwet - wiltpoint)) then
        dHplant = 0.98 * (hwet - wiltpoint)
        state%soilwater%Tactual = dHplant * kstem
      endif
      cw_Hxylem = cw_hleaf + dHPLant
      call JongvanLierLoop(state)

! --- check convergence
      if (abs(state%soilwater%Tactual - cw_qrosum) .lt. state%cfg%simulation%numerical%taccur) flconverg = .true.

! --- set initial values y3 and Fy3
      y3 = state%soilwater%Tactual
      Fy3 = state%soilwater%Tactual - cw_qrosum

! --- start search on correct value of qrosum
      do while (.not. flconverg)

! --- reset y1 and Fy1
      y1 = y3
      Fy1 = Fy3

      if (state%soilwater%Tactual .gt. cw_qrosum) then
! ---   root water extraction less than adopted, decrease tactual
          y2 = y1 - 0.001d0
      else
! ---   root water extraction larger than adopted, increase tactual
          y2 = y1 + 0.001d0
      endif

! --- calculate root extraction at modified pressure head y2 in leaf
      state%soilwater%Tactual = y2
      dHplant = state%soilwater%Tactual/kstem
      cw_Hxylem = cw_hleaf + dHPLant
      call JongvanLierLoop(state)
      Fy2 = state%soilwater%Tactual - cw_qrosum

! --- determine new value of tactual with Newton Raphson algorithm
        y3 = y1 - (y2-y1)*Fy1 / (Fy2-Fy1)
        maxstep = 0.05d0
        if ((y3-y1) .gt. 0.d0) then
          y3 = min((y1 + maxstep),y3)
        else
          y3 = max((y1 - maxstep),y3)
        endif
      dHPlantMax = hwet - wiltpoint
      y3 = min(y3,(dHPlantMax * Kstem))
      y3 = max(y3,0.d0)

! --- new estimated value of tactual is equal to y3!
 400  state%soilwater%Tactual = y3
      dHPlant = state%soilwater%Tactual/kstem
      cw_Hxylem = cw_hleaf + dHPlant

! --- calculate root extraction at y3
      call JongvanLierLoop(state)
      Fy3 = state%soilwater%Tactual - cw_qrosum

! --- check convergence
      if (abs(state%soilwater%Tactual - cw_qrosum) .lt. state%cfg%simulation%numerical%taccur) flconverg = .true.

! --- fatal error if too many iterations
      counter = counter + 1
      if (counter .gt. 1000) then
         messag = '5 Too many iterations for microscopic root'          &
     &             //' water uptake. Please adapt input!'
         call log_warn('rootextraction', messag)
         call fatalerr_collected ('rootextraction',messag)
      endif

! --- apply linear interpolation when Fy1 and Fy3 have opposite sign
      if (.not. flconverg) then
        if ((Fy1 .gt. 0.d0 .and. Fy3 .lt. 0.d0) .or.                    &
     &                     (Fy1 .lt. 0.d0 .and. Fy3 .gt. 0.d0)) then

          y3 = y1 + (y3 - y1) * Fy1 / (Fy1 - Fy3)
          goto 400
        endif
      endif

! --- WHILE-DO LOOP TO DETERMINE HLEAF WITH DROUGHT STRESS
      enddo

      dHplant = cw_qrosum/kstem
      cw_Hxylem = cw_hleaf + dHPlant

      endif

! --- qrot(node) has been calculated based on Jong van Lier (2013)

      cw_qrosum = 0.0d0
      do node = 1,noddrz
        cw_qrosum = cw_qrot(node) + cw_qrosum
      enddo

      if ( (at_ptra - cw_qrosum) .lt. state%cfg%simulation%numerical%taccur) then
! ---   compensate convergence error
        ratio = at_ptra / cw_qrosum
        do node = 1,noddrz
          cw_qrot(node) = cw_qrot(node) * ratio
        enddo
        cw_qrosum = at_ptra
        cw_alpJvLier = 1.d0
      else
! ---   drought reduction factor
        cw_alpJvLier = cw_qrosum / at_ptra
! ---   calculate potential qrot
        do node = 1,noddrz
          cw_qrot(node) = cw_qrot(node)/cw_alpJvLier
        enddo
      endif

      end associate
      return

  end subroutine JongvanLier

! ----------------------------------------------------------------------
!> Evaluate one microscopic uptake loop for a given xylem/leaf pressure state.
!!
!! @note
!! date: August 2016
!! purpose: Calculate microscopic root water uptake using hleaf
!! @endnote
  subroutine JongvanLierLoop(state)
! ----------------------------------------------------------------------
!     date      : August 2016
!     purpose   : Calculate microscopic root water uptake using hleaf
! [GR-CROP Phase B/9] narrow use variables
! ----------------------------------------------------------------------
      ! [SS-GR-FINAL B6] macp → swap_array_dimensions; remainder DEFERRED (same as RootExtraction)
      use swap_array_dimensions, only: macp
      use swap_log, only: log_warn
      use variables, only: &                                               ! [SS-GR-FINAL B6] residuals — all DEFERRED
                           ! adcrh/adcrl/aeratecrit/alphacrit retired  ! DEFERRED: crop stress config; Phase C3
                           botcom, criterhr, flhydrlift,                   &  ! dcritrtz/cumdens retired  ! DEFERRED: soil/crop config
                           ! hlim1/hlim2l/hlim2u/hlim3h/hlim3l/hlim4 retired  ! DEFERRED: drought limits
                           kroot, kstem, noddrz, oxygenintercept,   &  ! DEFERRED: crop/log globals
                           oxygenslope, rootcoefa, rooteff, &  ! rdctb retired  ! rd retired  ! DEFERRED: active crop state
                           rootradius, rxylem, stephr,                    &  ! saltmax/saltslope retired  ! DEFERRED: crop/solute config
                           swfrost,                                       &  ! swcompensate/swoxygen/swdrought retired  ! DEFERRED: stress switches
                           ! swsalinity/swwrtnonox/swstressor/swoxygentype retired  ! DEFERRED: stress switches
                           twilt, wiltpoint                          ! DEFERRED: convergence/stress params
      implicit none

      type(swap_state_t), intent(inout) :: state

! --- local variables
      integer counter,node,lay
      real(8) x1,x2,x3,Fx1,Fx2,qmax,conducsoil,conducroot
      real(8) step,mflux1,mflux2,depth
      character(len=200) messag

! --- SS-CRP Phase 2 C-2.5: ASSOCIATE for JvL loop per-node and scalar fields (state-only).
! --- SS-ATM Phase 2 Task A-2.3: at_ptra added for atmosphere reader cutover.
! --- SS-SWC S-2.7: sw_h, sw_theta added for soil-water-core reader cutover.
! --- SS-TC TC-12: dt added for timecontrol reader cutover.
      associate( &
         tc_dt      => state%timecontrol%dt,    &   ! TC-12
         at_ptra    => state%atmosphere%ptra,   &
         sw_h       => state%soilwater%h,       &   ! [SS-SWC S-2.7]
         sw_theta   => state%soilwater%theta,   &   ! [SS-SWC S-2.7]
         cw_qrot   => state%soilwater%qrot,    &
         cw_qrosum => state%soilwater%qrosum,  &
         cw_hroot  => state%soilwater%hroot,   &
         cw_mroot  => state%soilwater%mroot,   &
         cw_mflux  => state%soilwater%mflux,   &
         cw_rootphi => state%soilwater%rootphi, &
         cw_rootrho => state%soilwater%rootrho, &
         swscre     => state%timecontrol%swscre &  ! [SS-BMI2 Task 4]
      )

! --  initialisatie
      cw_qrosum = 0.d0
      x1 = 999.d0
      counter = 0
      ConducRoot = KRoot / rootradius / log(rootradius/rxylem)

      do node = 1,noddrz
! ---   determine h and matricflux potential at root-soil interface
        if (sw_h(node) .gt. -1.d0) then                                                ! [SS-SWC S-2.7]
! ---      very wet conditions
           lay = state%mesh%layer(node)  ! [GR-BH C7]
           ConducSoil = state%soilwater%ksatfit(lay) / (rootcoefa*state%soilwater%rmax(node)) /  &  ! [GR-BH C7]
     &                  log(rootcoefa*state%soilwater%rmax(node)/rootradius)
           cw_hroot(node) = (ConducSoil*sw_h(node) + ConducRoot*state%soilwater%Hxylem) /  &  ! [SS-SWC S-2.7]
     &                   (ConducSoil + ConducRoot)
        else
! ---      common conditions
          x3 = cw_hroot(node)
          do while (abs(x3-x1) .gt. CriterHr*log10(max(abs(x3),1.d0)))
            Step = min(abs(x3-x1),StepHr*log10(max(abs(x3),1.d0)))
            x1 = x3
            x2 = x1 - Step
            call MatricFlux(2,x1,node,mflux1,state)
            call MatricFlux(2,x2,node,mflux2,state)
            Fx1 = state%soilwater%Hxylem -x1 + cw_rootphi(node)*(cw_mflux(node)-mflux1)
            Fx2 = state%soilwater%Hxylem -x2 + cw_rootphi(node)*(cw_mflux(node)-mflux2)
            if (abs(Fx2-Fx1).lt.1.d-10) then
              x3 = x2
            else
              x3 = x1 - (x2-x1)*Fx1 / (Fx2-Fx1)
            endif
            x3 = max(x3,min(state%soilwater%Hxylem,sw_h(node)))   ! [SS-SWC S-2.7]
            x3 = min(x3,max(state%soilwater%Hxylem,sw_h(node)))   ! [SS-SWC S-2.7]
! ---       fatal error if too many iterations
            counter = counter + 1
            if (counter .gt. 50) then
              messag = '1 Too many iterations for microscopic root'    &
     &               //' water uptake. Please adapt input!'
              call log_warn('rootextraction', messag)
              call fatalerr_collected ('rootextraction',messag)
            endif
          enddo
          x1 = 999.d0
          counter = 0
          cw_hroot(node) = x3
        endif
        call MatricFlux(2,cw_hroot(node),node,cw_mroot(node),state)
! ---   calculate root water extraction flux
        if (node .lt. noddrz) then
          cw_qrot(node) = rooteff * cw_rootrho(node) *                        &
     &                 (cw_mflux(node)-cw_mroot(node)) * state%mesh%dz(node)  ! [GR-BH C7]
          if (cw_mflux(node) .gt. cw_mroot(node) ) then
! ---       water extraction, set maximum flux to 10% of available soil water
            qmax = (sw_theta(node) - twilt(node)) * state%mesh%dz(node) * 0.1d0 / tc_dt  ! [SS-SWC S-2.7] TC-12 [GR-BH C7]
            qmax = min(qmax,at_ptra)
            cw_qrot(node) = min(cw_qrot(node),qmax)
          else
! ---       possible hydraulic lift, set maximum flux to 0.1% change water content
            if (flhydrlift) then
              qmax = -0.001d0 * state%mesh%dz(node) / tc_dt  ! TC-12 [GR-BH Task 35]
              cw_qrot(node) = max(cw_qrot(node),qmax)
            else
! ---         no hydraulic lift allowed
              cw_qrot(node) = 0.d0
              cw_hroot(node) = sw_h(node)            ! [SS-SWC S-2.7]
              cw_mroot(node) = cw_mflux(node)
            endif
          endif
        else
! ---     last node, partly filled with roots
          depth = state%mesh%ztopcp(node) + state%crop%common%rd  ! [GR-BH C7]
          cw_qrot(node) = rooteff * cw_rootrho(node) *                        &
     &                 (cw_mflux(node)-cw_mroot(node)) * depth
          if (cw_mflux(node) .gt. cw_mroot(node) ) then
! ---       water extraction, set maximum flux to 10% of available soil water
            qmax = (sw_theta(node) - twilt(node)) * depth * 0.1d0 / tc_dt  ! [SS-SWC S-2.7] TC-12
            qmax = min(qmax,at_ptra)
            cw_qrot(node) = min(cw_qrot(node),qmax)
          else
! ---       possible hydraulic lift, set maximum flux to 0.1% change water content
            if (flhydrlift) then
              qmax = -0.001d0 * depth / tc_dt  ! TC-12
              cw_qrot(node) = max(cw_qrot(node),qmax)
            else
! ---         no hydraulic lift allowed
              cw_qrot(node) = 0.d0
              cw_hroot(node) = sw_h(node)            ! [SS-SWC S-2.7]
              cw_mroot(node) = cw_mflux(node)
            endif
          endif
        endif
        cw_qrosum = cw_qrot(node) + cw_qrosum
      enddo

      end associate
      return
  end subroutine JongvanLierLoop

! ----------------------------------------------------------------------
!> Initialize or evaluate the matric flux potential table.
!!
!! @param[in] task Task selector: 1 initializes lookup tables, 2 evaluates
!!                 flux potential at node pressure head.
!! @param[in] phead Pressure head [cm].
!! @param[in] node Node index.
!! @param[out] outcome Matric flux potential.
!!
!! @note
!! Date: February 2010
!! Purpose: Initialize and calculate matric flux potential
!! @endnote
  subroutine MatricFlux(task,phead,node,outcome,state)
! ----------------------------------------------------------------------
!     Date               : February 2010
!     Purpose            : Initialize and calculate matric flux potential
! SS-CRP Phase 2 C-2.5: mfluxtable retired from variables.f90.
!   task=1 writes state%soilwater%mfluxtable (state required).
!   task=2 reads state%soilwater%mfluxtable (state required).
! ----------------------------------------------------------------------

      use Variables, only: numlay, nod1lay, wiltpoint  ! [GR-CROP Phase B/9b] narrow; swsalinity/salthead retired
      use soilhydraulics_utils, only: watcon, hconduc
      implicit none

! --- local variables
      integer task,lay,count,start,node,i
      real(8) phead1,phead2,wcontent,conduc1,conduc2
      real(8) logphead,hosm,hsalt,mfluxsalt,phead,outcome
      type(swap_state_t), intent(inout), optional :: state

      select case (task)
      case (1)

! === initialization =========================================================
! --- SS-CRP Phase 2 C-2.5: write state%soilwater%mfluxtable directly.

      do lay = 1,numlay
        do count = 1,801
          state%soilwater%mfluxtable(lay,count) = 0.0d0
        enddo
      enddo

      start = int(100.d0*log10(-wiltpoint))
      do lay = 1,numlay
        phead1 = -10.d0**(dble(start)/100.d0)

!       find first Node of the Layer
        i = nod1lay(lay)

        wcontent = watcon(phead1, &
                           state%soilwater%vg_params(i), &
                           state%soilwater%iHWCKmodel(state%soilwater%layer(i)), &
                           i, state%soilwater)                     ! [SS-GR-UTILS Task 5]
        conduc1 = hconduc(phead1,wcontent,10.d0,state%heat%tsoil(i), &
                          state%soilwater%vg_params(i), &
                          state%soilwater%iHWCKmodel(state%soilwater%layer(i)), &
                          state%soilwater%fluseksatexm(i), &
                          i, state%soilwater)                      ! [SS-GR-UTILS Task 6]

        do count = start-1,1,-1
          phead2 = -10.d0**(dble(count)/100.d0)
          wcontent = watcon(phead2, &
                             state%soilwater%vg_params(i), &
                             state%soilwater%iHWCKmodel(state%soilwater%layer(i)), &
                             i, state%soilwater)                   ! [SS-GR-UTILS Task 5]
          conduc2 = hconduc(phead2,wcontent,10.d0,state%heat%tsoil(i), &
                             state%soilwater%vg_params(i), &
                             state%soilwater%iHWCKmodel(state%soilwater%layer(i)), &
                             state%soilwater%fluseksatexm(i), &
                             i, state%soilwater)                   ! [SS-GR-UTILS Task 6]
          state%soilwater%mfluxtable(lay,count) =                       &
     &         state%soilwater%mfluxtable(lay,count+1) +                &
     &         0.5d0 * (conduc1 + conduc2) * (phead2 - phead1)
          phead1 = phead2
          conduc1 = conduc2
        enddo
      enddo

      return

      case (2)

! === calculation of matric flux potential ===================================

! --- matric flux potential based on soil water pressure head
        lay = state%mesh%layer(node)  ! [GR-BH C7]
        if (phead .lt. wiltpoint) then
! ---     very dry range
          outcome = 0.0d0
        elseif (phead .gt. -1.023293d0) then
! ---     very wet range (> -10^0.01)
          outcome = state%soilwater%mfluxtable(lay,1) +                 &
     &              (phead+1.023293d0)*state%soilwater%ksatfit(lay)  ! [GR-BH C7]
        else
! ---     direct access table, with linear interpolation
          logphead = 100.d0*log10(-phead)
          count = int(logphead)
          outcome = (logphead-dble(count))*                             &
     &              state%soilwater%mfluxtable(lay,count+1) +           &
     &              (dble(count+1)-logphead)*                           &
     &              state%soilwater%mfluxtable(lay,count)
        endif

! --- correction matric flux potential for osmotic head due to salinity
      if (state%crop%common%swsalinity .eq. 2) then
          lay = state%mesh%layer(node)  ! [GR-BH C7]
!         osmotic head in cm
          hosm = state%crop%common%salthead * state%solute%cml(node)
          hsalt = wiltpoint + hosm
          if (hosm .lt. 1.d-3) then
! ---       very dry range
            mfluxsalt = 0.0d0
          elseif (hsalt .gt. -1.023293d0) then
! ---       very wet range (> -10^0.01)
            mfluxsalt = state%soilwater%mfluxtable(lay,1)
          else
! ---       direct access table, with linear interpolation
            logphead = 100.d0*log10(-hsalt)
            count = int(logphead)
            mfluxsalt = (logphead-dble(count))*                         &
     &                  state%soilwater%mfluxtable(lay,count+1) +       &
     &                  (dble(count+1)-logphead)*                       &
     &                  state%soilwater%mfluxtable(lay,count)
          endif
          outcome = outcome - mfluxsalt
          outcome = max(0.d0,outcome)
      endif

      case default
         call fatalerr_collected ('MatricFlux', 'Illegal value for TASK')
      end select

      return
  end subroutine MatricFlux

end module rootextraction_mod
