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
  public :: RootExtraction, matricflux_build_table, matric_flux

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
      use array_utils, only: afgen
      use oxygenstress_mod, only: OxygenStress
      implicit none

      type(swap_state_t), intent(inout) :: state

! --- local variables
      integer node
      real(8) top,bot,hlim3,hlim2,qred
      real(8) alpdry,alpwet,alpsol,alpfrs,alptot,vsmall
      real(8) alpdrycom,alpwetcom,alpsolcom,alpfrscom,alptotcom
      real(8) rd_noddrz, redtot

      parameter (vsmall = 1.0d-14)

! --- [GR-CROP 2026-05-25] sub-record aliases (crop, soil, atmo, mesh, heat, sol).
! --- [state%cfg-retirement cluster 6] soil_cfg dropped — only field accessed was
!     soil_cfg%frost%swfrost, now snapshotted as state%soilwater%swfrost.
      associate( &
         crop     => state%crop,                 &
         soil     => state%soilwater,            &
         atmo     => state%atmosphere,           &
         mesh     => state%mesh,                 &
         heat     => state%heat,                 &
         sol      => state%solute                &
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
      ! [GR-CROP 2026-05-25] swdrought=2 (de Jong van Lier microscopic uptake)
      ! dormant — body extracted to src/crop/dormant/jongvanlier.f90.
      if (crop%common%swdrought .eq. 2) then
        call fatalerr_collected('RootExtraction', 'swdrought=2 (de Jong van Lier) path is dormant — see src/crop/dormant/jongvanlier.f90')
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
              ! [GR-CROP 2026-05-25] swoxygen=2/swoxygentype=2 (OxygenReproFunction)
              ! is dormant — Task 4 (oxygenstress.f90 sub-arc) will fully migrate
              ! the OxygenReproFunction body + oxygenintercept/oxygenslope globals.
              call fatalerr_collected('RootExtraction', 'swoxygen=2/swoxygentype=2 (OxygenReproFunction) is dormant — see src/crop/oxygenstress.f90 (Task 4 will fully migrate)')
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
        if (soil%swfrost .eq.1 .and. heat%tsoil(node) .lt. 0.0d0) then
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
  subroutine matricflux_build_table(state)
! ----------------------------------------------------------------------
!     Date               : February 2010
!     Purpose            : Initialize matric flux potential lookup table
! Note: reached only when swdrought=2 (stub-errored on TOML path — dead path at runtime).
! ----------------------------------------------------------------------

      ! [GR-CROP 2026-05-25] MatricFlux is reached only when swdrought=2 (stub-errored
      ! on TOML path). The wiltpoint legacy global is always 0.0 on TOML — use
      ! state%crop%common%hlim4 (Feddes wilting-point pressure head) as the
      ! equivalent value. Code path stays dead at runtime.
      use soilhydraulics_utils, only: watcon, hconduc
      implicit none

      type(swap_state_t), intent(inout) :: state

! --- local variables
      integer lay,count,start,i
      real(8) phead1,phead2,wcontent,conduc1,conduc2

! === initialization =========================================================

      do lay = 1,state%mesh%numlay
        do count = 1,801
          state%soilwater%mfluxtable(lay,count) = 0.0d0
        enddo
      enddo

      start = int(100.d0*log10(-state%crop%common%hlim4))   ! [GR-CROP 2026-05-25] wiltpoint → hlim4 (dead path)
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
  end subroutine matricflux_build_table

  subroutine matric_flux(phead, node, outcome, state)
! ----------------------------------------------------------------------
!     Date               : February 2010
!     Purpose            : Calculate matric flux potential
! Note: reached only when swdrought=2 (stub-errored on TOML path — dead path at runtime).
! ----------------------------------------------------------------------
      implicit none

      real(8),            intent(in)    :: phead
      integer,            intent(in)    :: node
      real(8),            intent(out)   :: outcome
      type(swap_state_t), intent(inout) :: state

! --- local variables
      integer lay,count
      real(8) logphead,hosm,hsalt,mfluxsalt

! === calculation of matric flux potential ===================================

! --- matric flux potential based on soil water pressure head
        lay = state%mesh%layer(node)  ! [GR-BH C7]
        if (phead .lt. state%crop%common%hlim4) then   ! [GR-CROP 2026-05-25] wiltpoint → hlim4
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
          hsalt = state%crop%common%hlim4 + hosm   ! [GR-CROP 2026-05-25] wiltpoint → hlim4
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

      return
  end subroutine matric_flux

end module rootextraction_mod
