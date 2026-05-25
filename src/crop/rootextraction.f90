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
      use swap_array_dimensions, only: macp
      use swap_log, only: log_warn
      ! [GR-CROP 2026-05-25] retained — dead-branch (swoxygen=2 stub-errored in TOML);
      ! passed by reference to OxygenReproFunction.
      use variables, only: oxygenintercept, oxygenslope
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

! --- [GR-CROP 2026-05-25] sub-record aliases (crop, soil, atmo, mesh, heat, sol, cfg_soil).
      associate( &
         crop     => state%crop,                 &
         soil     => state%soilwater,            &
         atmo     => state%atmosphere,           &
         mesh     => state%mesh,                 &
         heat     => state%heat,                 &
         sol      => state%solute,               &
         cfg_soil => state%cfg%soil              &
      )

! --- reset writes — state-only.
      soil%qrot(1:mesh%numnod) = 0.0d0
      soil%Tactual    = soil%qrosum
      soil%qrosum     = 0.0d0
      soil%qreddrysum = 0.0d0
      soil%qredwetsum = 0.0d0
      soil%qredsolsum = 0.0d0
      soil%qredfrssum = 0.0d0

! --- skip routine if there are no roots
      if (crop%common%rd .lt. vsmall) return

! --- skip routine if transpiration rate is zero
      if (atmo%ptra .lt. 1.d-10) return

! --- DROUGHT REDUCTION ACCORDING TO FEDDES ET AL. (1978)
      if (crop%common%swdrought .eq. 1) then

! --- calculate potential root extraction of the compartments
        ! 22-10-2018: bug repair signalled by Paul van Walsum: division not by rd but by depth bottom of last compartment where roots are present
        !             rd replaced by (newly calculated) rd_noddrz
        rd_noddrz = abs(mesh%zbotcp(crop%common%noddrz))
        do node = 1,crop%common%noddrz
          top = abs(mesh%ztopcp(node) / rd_noddrz)
          bot = abs(mesh%zbotcp(node) / rd_noddrz)
          soil%qrot(node) = (afgen(crop%common%cumdens,202,bot)-afgen(crop%common%cumdens,202,top))* atmo%ptra
        enddo

! --- calculating critical point hlim3 according to feddes
        if (atmo%atmdem .lt. crop%common%adcrl) then
          hlim3 = crop%common%hlim3l
        elseif (atmo%atmdem .le. crop%common%adcrh) then
          hlim3 = crop%common%hlim3h + ((crop%common%adcrh - atmo%atmdem) / (crop%common%adcrh - crop%common%adcrl)) * (crop%common%hlim3l - crop%common%hlim3h)
        else
          hlim3 = crop%common%hlim3h
        endif
      endif

! --- DROUGHT REDUCTION ACCORDING TO DE JONG VAN LIER ET AL. (2012)
      if (crop%common%swdrought .eq. 2) then
        call JongvanLier(state)
      endif

! === COMBINATION OF OXYGEN, DROUGHT, SALT AND FROST STRESS ====

      soil%qrosum = 0.0d0

      do 200 node = 1,crop%common%noddrz
        alpdry = 1.0d0
        alpwet = 1.0d0
        alpsol = 1.0d0
        alpfrs = 1.0d0

! ---   reduction due to oxygen stress
        if (crop%common%swoxygen .ne. 0) then

          ! Feddes linear reduction based on pressure head
          if (crop%common%swoxygen .eq. 1) then

            if (node.gt.mesh%botcom(1)) then
              hlim2 = crop%common%hlim2l
            else
              hlim2 = crop%common%hlim2u
            endif
            if (soil%h(node).le.crop%common%hlim1.and.soil%h(node).gt.hlim2) then
              alpwet = (crop%common%hlim1-soil%h(node))/(crop%common%hlim1-hlim2)
            endif
            if (soil%h(node).gt.crop%common%hlim1) then
              alpwet = 0.0d0
            endif

          ! Bartholomeus non-linear reduction based on gas filled porosity
          elseif (crop%common%swoxygen .eq. 2) then

            ! use physical processes
            if (crop%common%swoxygentype .eq. 1) then
!##MH         call OxygenStress(node,alpwet,ResultsOxygenStress)
              call OxygenStress(node,alpwet,state)

            ! use reproduction functions
            else
              call OxygenReproFunction (OxygenSlope,OxygenIntercept,soil%theta,soil%thetas,heat%tsoil,node,mesh%z,mesh%dz,alpwet,state)
            endif

          endif

! ---     Stop root development in case of oxgenstress at noddrz
!         WOFOST: root zone remain the aim, but biomass is increasing
!         GRASS : stop root development
          soil%flWrtNonox = .false.
          if (crop%common%swWrtNonox .eq. 1 .and. node .eq. crop%common%noddrz) then
            if (alpwet .lt. crop%common%aeratecrit) then
              soil%flWrtNonox = .true.
            end if
          endif

        endif

! ---   reduction due to drought stress

        ! Feddes linear reduction based on pressure head
        if (crop%common%swdrought .eq. 1) then
          if (soil%h(node) .lt. crop%common%hlim4) then
            alpdry = 0.0d0
          elseif (soil%h(node).le.hlim3) then
            alpdry = (crop%common%hlim4-soil%h(node))/(crop%common%hlim4-hlim3)
          endif
        endif

        ! JongvanLier microscopic concept for drought
        if (crop%common%swdrought .eq. 2) then
          alpdry = soil%alpJvLier
        endif

! ---   reduction due to salt stress

        ! reduction according to Maas and Hoffman linear reduction function
        if (crop%common%swsalinity .eq. 1) then
          if (sol%cml(node) .gt. crop%common%saltmax) then
            alpsol = 1.0d0 - (sol%cml(node) - crop%common%saltmax) * crop%common%saltslope
            alpsol = max(0.0d0,alpsol)
          endif
        endif
! ---   mind: in case of salt stress with osmotic head, microscopic root water extraction
! ---         according to JongvanLier (2013) should be used (swsalinity = 2); in that case
! ---         salinity stress is included in drought stress and not separately specified
! ---         in output file *.STR

! ----  reduction due to frost conditions
        if (cfg_soil%frost%swfrost .eq.1 .and. heat%tsoil(node) .lt. 0.0d0) then
          alpfrs = 0.0d0
        endif

! ----  overall reduction
        soil%qpotrot(node) = soil%qrot(node)
        soil%qrot(node)    = soil%qrot(node) * alpwet * alpdry * alpsol * alpfrs
        soil%qrosum        = soil%qrot(node) + soil%qrosum

! ----  apportionment to different types stresses (cm)

        qred = soil%qpotrot(node) - soil%qrot(node)
        if (qred .lt. vsmall)then

          ! no stress
          soil%qredwet(node) = 0.d0
          soil%qreddry(node) = 0.d0
          soil%qredsol(node) = 0.d0
          soil%qredfrs(node) = 0.d0

        else

          ! multiplication of stressors
          alptot = (1 - alpwet) + (1 - alpdry) + (1 - alpsol) + (1 - alpfrs)

          ! contribution of each stressor (lineair approach)
          soil%qredwet(node) = (1 - alpwet) / alptot * qred
          soil%qreddry(node) = (1 - alpdry) / alptot * qred
          soil%qredsol(node) = (1 - alpsol) / alptot * qred
          soil%qredfrs(node) = (1 - alpfrs) / alptot * qred

          ! sum of each stressor (rootzone)
          soil%qredwetsum = soil%qredwetsum + soil%qredwet(node)
          soil%qreddrysum = soil%qreddrysum + soil%qreddry(node)
          soil%qredsolsum = soil%qredsolsum + soil%qredsol(node)
          soil%qredfrssum = soil%qredfrssum + soil%qredfrs(node)

        end if

200   continue

! --- compensated root water uptake according to Jarvis (1989) or Walsum (2020)
      if (crop%common%swcompensate .gt. 0) then

        ! compensated root water uptake according to Walsum
        if (crop%common%swcompensate .eq. 2) then
            crop%common%alphacrit = min((crop%common%dcritrtz + crop%common%rdm - rd_noddrz) / crop%common%rdm, 1.0d0)
        end if

        alptot = soil%qrosum / atmo%ptra
        qred = atmo%ptra - soil%qrosum
        if (abs(crop%common%alphacrit - 1.0d0) .ge. vsmall .and. qred .gt. vsmall .and. alptot .ge. 0.05d0) then
          ! Only compensation when rootextraction and transpiration reduction is greater than vsmall
          ! and when alptot > 0.05, i.e. when there is less than 95% stress reduction. This minimum is
          ! also important for the approximation of alp... in the next 4 lines.
          alpdry = alptot**(soil%qreddrysum/qred)
          alpwet = alptot**(soil%qredwetsum/qred)
          alpsol = alptot**(soil%qredsolsum/qred)
          alpfrs = alptot**(soil%qredfrssum/qred)

          if (crop%common%swstressor .eq. 1) then
            alptotcom = min(alptot / crop%common%alphacrit, 1.d0)
            alpdrycom = alpdry
            alpwetcom = alpwet
            alpsolcom = alpsol
            alpfrscom = alpfrs
          else
            alpdrycom = alpdry
            alpwetcom = alpwet
            alpsolcom = alpsol
            alpfrscom = alpfrs
            if (crop%common%swstressor .eq. 2) then
              alpdrycom = min(alpdry / crop%common%alphacrit, 1.d0)
            elseif (crop%common%swstressor .eq. 3) then
              alpwetcom = min(alpwet / crop%common%alphacrit, 1.d0)
            elseif (crop%common%swstressor .eq. 4) then
              alpsolcom = min(alpsol / crop%common%alphacrit, 1.d0)
            elseif (crop%common%swstressor .eq. 5) then
              alpfrscom = min(alpfrs / crop%common%alphacrit, 1.d0)
            endif
            alptotcom = alpwetcom * alpdrycom * alpsolcom * alpfrscom
          endif

          ! Change the abstraction of the roots
          do node = 1,crop%common%noddrz
            soil%qrot(node) = soil%qrot(node) * alptotcom / alptot
          enddo

          ! Change the sum-parameters
          soil%qrosum = atmo%ptra * alptotcom
          qred = atmo%ptra - soil%qrosum
          if (qred .lt. vsmall) then
            ! There is no stress.
            soil%qredwetsum = 0.0d0
            soil%qreddrysum = 0.0d0
            soil%qredsolsum = 0.0d0
            soil%qredfrssum = 0.0d0
          else
            redtot = (1 - alpwetcom) + (1 - alpdrycom) + (1 - alpsolcom) + (1 - alpfrscom)
            soil%qredwetsum = (1 - alpwetcom) / redtot * qred
            soil%qreddrysum = (1 - alpdrycom) / redtot * qred
            soil%qredsolsum = (1 - alpsolcom) / redtot * qred
            soil%qredfrssum = (1 - alpfrscom) / redtot * qred
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
! ----------------------------------------------------------------------
      use swap_array_dimensions, only: macp
      use swap_log, only: log_warn
      ! [GR-CROP 2026-05-25] retained — dead-branch (swdrought=2 stub-errored in TOML)
      use variables, only: kroot, kstem, rootcoefa,        &
                           rootradius, rxylem, stephr,     &
                           wiltpoint
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

! --- [GR-CROP 2026-05-25] sub-record aliases (crop, soil, atmo, mesh, cfg_sim).
      associate( &
         crop    => state%crop,                       &
         soil    => state%soilwater,                  &
         atmo    => state%atmosphere,                 &
         mesh    => state%mesh,                       &
         cfg_sim => state%cfg%simulation,             &
         noddrz  => state%crop%common%noddrz          &  ! [GR-CROP 2026-05-25] alias for state read
      )

! --- initialization
      phi = 3.1415926d0
      counter = 0
      flstress = .false.
      hleafm1 = soil%hleaf
      flconverg = .false.

! --- reset values below root zone to zero
      do node = noddrz+1,mesh%numnod
        soil%mflux(node)   = 0.d0
        soil%mroot(node)   = 0.d0
        soil%hroot(node)   = 0.d0
        soil%rootrho(node) = 0.d0
        soil%rootphi(node) = 0.d0
      enddo

! --- give hroot a value when root zone becomes larger and store previous values
      do node = 1, noddrz
        if (abs(soil%hroot(node)) .lt. 1.0d-12) then
          soil%hroot(node) = soil%h(node)
        end if
        hrootm1(node) = soil%hroot(node)
      enddo

! --- initialization of rootrho and rootphi
      do node = 1,noddrz-1
        reldepth = -mesh%z(node)/crop%common%rd
        rdensity = afgen(crop%common%rdctb,22,reldepth)
        soil%rmax(node) = 1.d0/dsqrt(phi*rdensity)
        soil%rootrho(node) = 4.d0/(rootradius*rootradius-rootcoefa*soil%rmax(node)&
     &      *rootcoefa*soil%rmax(node) + 2.d0 * (soil%rmax(node)*soil%rmax(node) +     &
     &      rootradius*rootradius)*log(rootcoefa*soil%rmax(node)/rootradius))
        soil%rootphi(node) = soil%rootrho(node) * soil%rmax(node)*soil%rmax(node) *         &
     &     log(rootradius/rxylem) * 0.5d0 / kroot
      enddo
! --- last node, partly filled with roots
      node = noddrz
      meandepth = (mesh%ztopcp(node)-crop%common%rd)*0.5d0
      reldepth = meandepth/(-crop%common%rd)
      rdensity = afgen(crop%common%rdctb,22,reldepth)
      soil%rmax(node) = 1.d0/dsqrt(phi*rdensity)
      soil%rootrho(node) = 4.d0/(rootradius*rootradius-rootcoefa*soil%rmax(node)  &
     &      *rootcoefa*soil%rmax(node) + 2.d0 * (soil%rmax(node)*soil%rmax(node) +     &
     &      rootradius*rootradius)*log(rootcoefa*soil%rmax(node)/rootradius))
      soil%rootphi(node) = soil%rootrho(node) * soil%rmax(node)*soil%rmax(node) *           &
     &     log(rootradius/rxylem) * 0.5d0 / kroot

! --- calculate current matric flux potential in soil water
      do node = 1,noddrz
        soil%h(node) = min(1.0d3,max(soil%h(node),1.0d-8))
        call MatricFlux(2,soil%h(node),node,soil%mflux(node),state)
      enddo

! --- interpolate matric flux potential in lowest compartment that is partly filled with roots
      node = noddrz
      if (node.gt.1) then
         soil%mflux(node) = soil%mflux(node) + (soil%mflux(node-1)-soil%mflux(node))/       &
     &                 mesh%disnod(node) * (meandepth-mesh%z(node))
      endif

! --- determine highest pressure head in root zone
      hwet = soil%h(1)
      do node = 2,noddrz
        if (soil%h(node) .gt. hwet) hwet = soil%h(node)
      enddo

! --- determine whether qrosum < ptra (flstress = true)
! --- calculate root water extraction qrosum at Hleaf = wiltpoint and Tactual = Ptra
      soil%hleaf = wiltpoint
      dHPlant = atmo%ptra/Kstem
      if (dHPlant .gt. (hwet - wiltpoint)) flstress = .true.

      soil%Hxylem = soil%hleaf + dHPLant
      call JongvanLierLoop(state)
      if (soil%qrosum .lt. atmo%ptra) flstress = .true.

! --- set hroot to previous values
      do node = 1,noddrz
        soil%hroot(node) = hrootm1(node)
      enddo

! --- FLSTRESS = FALSE: DETERMINE HLEAF WITH TACTUAL = PTRA
      if (.not. flstress) then

! --- calculate root extraction at previous value of hleaf
      soil%hleaf = Hleafm1
      soil%Hxylem = soil%hleaf + dHPLant
      call JongvanLierLoop(state)

! --- check convergence
      if (abs(soil%qrosum-atmo%ptra) .lt. cfg_sim%numerical%taccur) flconverg = .true.

! --- set initial values y3 and Fy3
      y3 = soil%hleaf
      Fy3 = atmo%ptra - soil%qrosum

! --- start search on correct value of pressure head in leaves
      do while (.not. flconverg)

! --- reset y1 and Fy1
      y1 = y3
      Fy1 = Fy3

      if (soil%qrosum .lt. atmo%ptra) then
! ---   root water extraction too small, decrease Hleaf
          y2 = y1 - StepHr*log10(max(abs(y1),1.d0))
      else
! ---   root water extraction too large, increase Hleaf
          y2 = y1 + StepHr*log10(max(abs(y1),1.d0))
      endif

! --- calculate root extraction at modified pressure head y2 in leaf
      soil%hleaf = y2
      soil%Hxylem = soil%hleaf + dHPLant
      call JongvanLierLoop(state)
      Fy2 = atmo%ptra - soil%qrosum

! --- determine new value of leaf pressure head with Newton Raphson algorithm
      if (abs(Fy1-Fy2).lt. 1.d-10) then
! ---   take maximum step
        if (soil%qrosum .gt. atmo%ptra) then
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
 300  soil%hleaf = y3
      soil%Hxylem = soil%hleaf + dHPlant

! --- calculate root extraction at y3
      call JongvanLierLoop(state)
      Fy3 = atmo%ptra - soil%qrosum

! --- check convergence
      if (abs(soil%qrosum-atmo%ptra) .lt. cfg_sim%numerical%taccur .or.                             &
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

      soil%qrosum = atmo%ptra
      soil%hleaf = y3
      dHplant = soil%qrosum/kstem
      soil%Hxylem = soil%hleaf + dHPlant

      else
! --- FLSTRESS = TRUE: DETERMINE QROSUM WITH HLEAF = WILTPOINT
      counter = 0

! --- calculate root extraction at previous value of qrosum
      soil%hleaf = wiltpoint
      dHplant = soil%Tactual/kstem
      if (dHplant .gt. (hwet - wiltpoint)) then
        dHplant = 0.98 * (hwet - wiltpoint)
        soil%Tactual = dHplant * kstem
      endif
      soil%Hxylem = soil%hleaf + dHPLant
      call JongvanLierLoop(state)

! --- check convergence
      if (abs(soil%Tactual - soil%qrosum) .lt. cfg_sim%numerical%taccur) flconverg = .true.

! --- set initial values y3 and Fy3
      y3 = soil%Tactual
      Fy3 = soil%Tactual - soil%qrosum

! --- start search on correct value of qrosum
      do while (.not. flconverg)

! --- reset y1 and Fy1
      y1 = y3
      Fy1 = Fy3

      if (soil%Tactual .gt. soil%qrosum) then
! ---   root water extraction less than adopted, decrease tactual
          y2 = y1 - 0.001d0
      else
! ---   root water extraction larger than adopted, increase tactual
          y2 = y1 + 0.001d0
      endif

! --- calculate root extraction at modified pressure head y2 in leaf
      soil%Tactual = y2
      dHplant = soil%Tactual/kstem
      soil%Hxylem = soil%hleaf + dHPLant
      call JongvanLierLoop(state)
      Fy2 = soil%Tactual - soil%qrosum

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
 400  soil%Tactual = y3
      dHPlant = soil%Tactual/kstem
      soil%Hxylem = soil%hleaf + dHPlant

! --- calculate root extraction at y3
      call JongvanLierLoop(state)
      Fy3 = soil%Tactual - soil%qrosum

! --- check convergence
      if (abs(soil%Tactual - soil%qrosum) .lt. cfg_sim%numerical%taccur) flconverg = .true.

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

      dHplant = soil%qrosum/kstem
      soil%Hxylem = soil%hleaf + dHPlant

      endif

! --- qrot(node) has been calculated based on Jong van Lier (2013)

      soil%qrosum = 0.0d0
      do node = 1,noddrz
        soil%qrosum = soil%qrot(node) + soil%qrosum
      enddo

      if ( (atmo%ptra - soil%qrosum) .lt. cfg_sim%numerical%taccur) then
! ---   compensate convergence error
        ratio = atmo%ptra / soil%qrosum
        do node = 1,noddrz
          soil%qrot(node) = soil%qrot(node) * ratio
        enddo
        soil%qrosum = atmo%ptra
        soil%alpJvLier = 1.d0
      else
! ---   drought reduction factor
        soil%alpJvLier = soil%qrosum / atmo%ptra
! ---   calculate potential qrot
        do node = 1,noddrz
          soil%qrot(node) = soil%qrot(node)/soil%alpJvLier
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
! ----------------------------------------------------------------------
      use swap_array_dimensions, only: macp
      use swap_log, only: log_warn
      ! [GR-CROP 2026-05-25] retained — dead-branch (swdrought=2 stub-errored in TOML)
      use variables, only: criterhr, flhydrlift, kroot,    &
                           rootcoefa, rooteff, rootradius, &
                           rxylem, stephr, twilt
      implicit none

      type(swap_state_t), intent(inout) :: state

! --- local variables
      integer counter,node,lay
      real(8) x1,x2,x3,Fx1,Fx2,qmax,conducsoil,conducroot
      real(8) step,mflux1,mflux2,depth
      character(len=200) messag

! --- [GR-CROP 2026-05-25] sub-record aliases (crop, soil, atmo, mesh, time).
      associate( &
         crop   => state%crop,                        &
         soil   => state%soilwater,                   &
         atmo   => state%atmosphere,                  &
         mesh   => state%mesh,                        &
         time   => state%timecontrol,                 &
         noddrz => state%crop%common%noddrz           &  ! [GR-CROP 2026-05-25] alias for state read
      )

! --  initialisatie
      soil%qrosum = 0.d0
      x1 = 999.d0
      counter = 0
      ConducRoot = KRoot / rootradius / log(rootradius/rxylem)

      do node = 1,noddrz
! ---   determine h and matricflux potential at root-soil interface
        if (soil%h(node) .gt. -1.d0) then
! ---      very wet conditions
           lay = mesh%layer(node)
           ConducSoil = soil%ksatfit(lay) / (rootcoefa*soil%rmax(node)) /  &
     &                  log(rootcoefa*soil%rmax(node)/rootradius)
           soil%hroot(node) = (ConducSoil*soil%h(node) + ConducRoot*soil%Hxylem) /  &
     &                   (ConducSoil + ConducRoot)
        else
! ---      common conditions
          x3 = soil%hroot(node)
          do while (abs(x3-x1) .gt. CriterHr*log10(max(abs(x3),1.d0)))
            Step = min(abs(x3-x1),StepHr*log10(max(abs(x3),1.d0)))
            x1 = x3
            x2 = x1 - Step
            call MatricFlux(2,x1,node,mflux1,state)
            call MatricFlux(2,x2,node,mflux2,state)
            Fx1 = soil%Hxylem -x1 + soil%rootphi(node)*(soil%mflux(node)-mflux1)
            Fx2 = soil%Hxylem -x2 + soil%rootphi(node)*(soil%mflux(node)-mflux2)
            if (abs(Fx2-Fx1).lt.1.d-10) then
              x3 = x2
            else
              x3 = x1 - (x2-x1)*Fx1 / (Fx2-Fx1)
            endif
            x3 = max(x3,min(soil%Hxylem,soil%h(node)))
            x3 = min(x3,max(soil%Hxylem,soil%h(node)))
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
          soil%hroot(node) = x3
        endif
        call MatricFlux(2,soil%hroot(node),node,soil%mroot(node),state)
! ---   calculate root water extraction flux
        if (node .lt. noddrz) then
          soil%qrot(node) = rooteff * soil%rootrho(node) *                        &
     &                 (soil%mflux(node)-soil%mroot(node)) * mesh%dz(node)
          if (soil%mflux(node) .gt. soil%mroot(node) ) then
! ---       water extraction, set maximum flux to 10% of available soil water
            qmax = (soil%theta(node) - twilt(node)) * mesh%dz(node) * 0.1d0 / time%dt
            qmax = min(qmax,atmo%ptra)
            soil%qrot(node) = min(soil%qrot(node),qmax)
          else
! ---       possible hydraulic lift, set maximum flux to 0.1% change water content
            if (flhydrlift) then
              qmax = -0.001d0 * mesh%dz(node) / time%dt
              soil%qrot(node) = max(soil%qrot(node),qmax)
            else
! ---         no hydraulic lift allowed
              soil%qrot(node) = 0.d0
              soil%hroot(node) = soil%h(node)
              soil%mroot(node) = soil%mflux(node)
            endif
          endif
        else
! ---     last node, partly filled with roots
          depth = mesh%ztopcp(node) + crop%common%rd
          soil%qrot(node) = rooteff * soil%rootrho(node) *                        &
     &                 (soil%mflux(node)-soil%mroot(node)) * depth
          if (soil%mflux(node) .gt. soil%mroot(node) ) then
! ---       water extraction, set maximum flux to 10% of available soil water
            qmax = (soil%theta(node) - twilt(node)) * depth * 0.1d0 / time%dt
            qmax = min(qmax,atmo%ptra)
            soil%qrot(node) = min(soil%qrot(node),qmax)
          else
! ---       possible hydraulic lift, set maximum flux to 0.1% change water content
            if (flhydrlift) then
              qmax = -0.001d0 * depth / time%dt
              soil%qrot(node) = max(soil%qrot(node),qmax)
            else
! ---         no hydraulic lift allowed
              soil%qrot(node) = 0.d0
              soil%hroot(node) = soil%h(node)
              soil%mroot(node) = soil%mflux(node)
            endif
          endif
        endif
        soil%qrosum = soil%qrot(node) + soil%qrosum
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
! Note: task=1 writes state%soilwater%mfluxtable; task=2 reads it.
! ----------------------------------------------------------------------

      ! [GR-CROP 2026-05-25] retained — dead-branch (swdrought=2 stub-errored in TOML)
      use variables, only: wiltpoint
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

      do lay = 1,state%mesh%numlay
        do count = 1,801
          state%soilwater%mfluxtable(lay,count) = 0.0d0
        enddo
      enddo

      start = int(100.d0*log10(-wiltpoint))
      do lay = 1,state%mesh%numlay
        phead1 = -10.d0**(dble(start)/100.d0)

!       find first Node of the Layer
        i = state%mesh%nod1lay(lay)

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
