module distribute_drainage
   use error_mod, only: fatalerr_collected
!> Distribute Drainage - Lateral waterflux simulation in saturated zone
!!
!! This module provides routines for simulating lateral water fluxes in the saturated zone
!! of soil profiles, with support for multiple drainage systems and infiltration scenarios.
!!
!! The module implements the distribution of drainage and infiltration fluxes over
!! model discharge layers based on transmissivity-weighted allocation. It handles:
!!
!! * Multiple drainage systems with varying spacing and depths
!! * Anisotropic hydraulic conductivity (horizontal vs vertical)
!! * Interflow (surface runoff) as highest-order drainage
!! * Separate distribution of infiltration fluxes from surface water
!! * Dynamic groundwater level effects on discharge layer geometry
!!
!!### Physical Background
!!
!! The approach divides the saturated zone into model discharge layers, with each
!! drainage system associated with a specific layer. Fluxes are distributed within
!! each layer proportionally to the local horizontal transmissivity:
!!
!! \[ T = \sum_{i} K_{h,i} \cdot \Delta z_i \]
!!
!! where \(K_{h,i}\) is the horizontal saturated hydraulic conductivity and
!! \(\Delta z_i\) is the thickness of compartment \(i\).
!!
!! The maximum depth of each discharge layer is constrained by:
!!
!! \[ D_{max} = 0.25 \cdot L \cdot \sqrt{\frac{K_v}{K_h}} + z_{wt} \]
!!
!! where \(L\) is drain spacing, \(K_v/K_h\) is the anisotropy factor, and
!! \(z_{wt}\) is the groundwater level.
!!
!!### Module History
!!
!! @history
!! - **1999-09-29**: Initial implementation (J. Kroes)
!! - **2003-07-03**: Added anisotropy factor COFANI = Khor/Kvert (J. Kroes)
!! - **2004-09-30**: Added interflow discharge layer option (R. Hendriks)
!! - **2010-10-29**: Added Lev2Comp helper subroutine
!! - **2018-01-28**: Fixed transmissivity bug for infiltration (R. Hendriks)
!! - **2026**: Converted to module structure for modernization
!! @endhistory
!!
!!### References
!!
!! The drainage flux distribution follows the approach described in the SWAP model
!! documentation (van Dam et al., 2008). The infiltration flux distribution extension
!! is documented in Hendriks (2010).
!!
!!@note
!! This module requires the `variables` module for groundwater table data (`nowltab`)
!! and the `arrays.fi` include file for array dimension parameters (MACP, MADR, MAHO, MAOWL).
!!@endnote
!!
!!
   implicit none

contains

   SUBROUTINE divdra(NumComp, NumDrain, ThickComp, ksatfit, ksatexm,   &
        &                   fluseksatexm, layer,                            &
        &                   cofani, gwlev, DistDrain, FluxDr, FluxDrComp,   &
        &                   Swdivdinf, Swnrsrf, SwTopnrsrf, Zbotdr,           &  !  Divdra, infiltration
        &                   dt, FacDpthInf, owltab, t1900)                       !  Divdra, infiltration
    !! Distribute drainage and infiltration fluxes over soil compartments
    !!
    !! This is the main workhorse subroutine that distributes lateral drainage and
    !! infiltration fluxes across soil compartments based on transmissivity weighting.
    !!
    !! The algorithm proceeds through 10 major steps:
    !!
    !! 1. Calculate horizontal and vertical saturated conductivities for each compartment
    !! 2. Find compartment containing the groundwater level
    !! 3. Handle interflow (surface runoff) as highest-order drainage if enabled
    !! 4. Calculate profile-averaged anisotropy factor and cumulative transmissivity
    !! 5. Determine maximum depth for each discharge layer
    !! 6. Sort drainage systems by priority (spacing × depth factor)
    !! 7. Calculate bottom of first-order discharge layer
    !! 8. Calculate bottoms of higher-order discharge layers
    !! 9. Apply depth constraints based on drain spacing and anisotropy
    !! 10. Distribute fluxes proportionally to transmissivity within each layer
    !!
    !! If separate infiltration distribution is enabled (`Swdivdinf=1`), infiltration
    !! fluxes are distributed with a parabolic weighting that accounts for both
    !! saturated and unsaturated zone transmissivity.
    !!
    !!@note
    !! The transmissivity-weighted distribution assumes that lateral flow is proportional
    !! to the product of hydraulic conductivity and thickness at each depth.
    !!@endnote
    !!
    !!@warning
    !! This subroutine modifies `FluxDrComp` in place. Ensure it is properly initialized
    !! before calling.
    !!@endwarning

      use variables, only: nowltab
      use array_utils, only: afgen
      use swap_array_dimensions, only: macp, madr, maho, maowl

      INTEGER DrainSequence(Madr), Icomp, idr, iidr
      INTEGER idum, jdr, NumComWatLev, NumComBotDislay(Madr), NumActDrain
      INTEGER NumDrain, NumComp, LAYER(MACP)
      INTEGER Swdivdinf, Swnrsrf, SwTopnrsrf, NumDrHlp, NumComBotIntflw !  Divdra, infiltration

      real(8) FacAniso, CondSatHorAv, CondSatVerAV, BotDisLay(Madr)
      real(8) ThickCum
      real(8) CumThickCondHor, CumThickCondVer, DistDrain(Madr), Dum1, Dum2
      real(8) FluxDr(:)          ! assumed-shape (was FluxDr(Madr))
      real(8) FluxDrComp(:,:)    ! assumed-shape (was FluxDrComp(Madr,MACP))
      real(8) FlowDrDisch(Madr)
      real(8) HelpFl(Madr)
      real(8) HelpTh, HelpTr, MaxDepthDislay(Madr), CondSatHor(MACP)
      real(8) CondSatVer(MACP), Small, ThickComp(MACP), ThickCompSatWatLev
      real(8) ThickCompBotDislay(Madr), Transmissivity(Madr)
      real(8) WatLevAv, GWLEV
      real(8) COFANI(MAHO), ThickDislay

      real(8) ksatfit(MAHO), ksatexm(MAHO), Zbotdr(Madr)
      real(8) ThickCompBotIntflw, ThickCompUndIntflw, TransmisIntflw
      logical fluseksatexm(macp)

      LOGICAL NonExistFluxDr

!   - Declarations Divdra Infiltration
      real(8) dt, FacDpthInf, owltab(Madr, 2*maowl), t1900

      integer NumComDpthInf, NumComSrfLev
      real(8) CumKD, DpthInflay, FDisInf(Madr), FDisInfmin
      real(8) FluxComWatLev, KD(Macp), RQmax, SurfLev   !, temptab(2*maowl)
      real(8) ThickCompAbvDpthInf, ThickCompBlwSrfLev
      real(8) ThickCompUnsWatLev, TransmisSat, TransmisTot
      real(8) TransmisUns

      LOGICAL flDivdInf

      DATA Small/1.0d-10/

! ----------------------------------------------------------------------
! --- paragraph 1
! --- initial calculations
!
! --- groundwater level converted to cm below surface level
      WatLevAv = -1.0d0*MIN(GWLEV, 0.0d0)

! --- calculation of CondSatHor and CondSatVer
      do 10 icomp = 1, NumComp
         if (fluseksatexm(icomp)) then
            CondSatHor(icomp) = ksatexm(LAYER(icomp))*COFANI(LAYER(icomp))
            CondSatVer(icomp) = ksatexm(LAYER(icomp))
         else
            CondSatHor(icomp) = ksatfit(LAYER(icomp))*COFANI(LAYER(icomp))
            CondSatVer(icomp) = ksatfit(LAYER(icomp))
         end if
10       continue

! --- Initialisation and test whether any drainflux exists
         NonExistFluxDr = .true.
         flDivdInf = .false.                                          !  Divdra, infiltration
         do idr = 1, NumDrain
            if (abs(FluxDr(idr)) .gt. Small) NonExistFluxDr = .false.
            do icomp = 1, NumComp
               FluxDrComp(idr, icomp) = 0.0d0
            end do
            DrainSequence(idr) = 0
!   - infiltration flux exists and Swdivdinf = 1?  Yes, then distribution of
!     infiltration fluxes is separate from distribution discharge fluxes
            if (Swdivdinf .eq. 1 .and. FluxDr(idr) .lt. -Small) flDivdInf = .true. !  Divdra, infiltration
         end do
         if (NonExistFluxDr) return

! --- paragraph 2
! --- Search for compartment with groundwaterlevel: NumComWatLev
! --- saturated part of compartment with waterlevel: ThickCompSatWatLev
         Call Lev2Comp(NumComp, WatLevAv, ThickComp, NumComWatLev,    &
        &              ThickCum, ThickCompUnsWatLev, ThickCompSatWatLev)

! --- paragraph 3
! --- RobHen 30-09-2004; in case of interflow: the bottom of the highest order
! --- drainage system (-Zbotdr(NumDrain)) represents max depth of the interflow
! --- discharge layer and top of all other model_discharge_layers
         NumDrHlp = NumDrain
!     allow disabling of option
         if ((Swnrsrf .eq. 1 .or. Swnrsrf .eq. 2) .and. Swtopnrsrf .eq. 1) then
            NumDrHlp = NumDrain - 1
            if (FluxDr(NumDrain) .gt. Small) then
! --- 3a determine bottom and total transmissivity of the interflow
! --- discharge layer
               NumComBotIntflw = NumComWatLev
               TransmisIntflw = CondSatHor(NumComBotIntflw)*              &
         &                       ThickCompSatWatLev
               do while (-ZBotDr(NumDrain) .gt. ThickCum)
                  NumComBotIntflw = NumComBotIntflw + 1
                  ThickCum = ThickCum + ThickComp(NumComBotIntflw)
                  TransmisIntflw = TransmisIntflw +                         &
          &                        CondSatHor(NumComBotIntflw)*             &
          &                        ThickComp(NumComBotIntflw)
               end do
               ThickCompUndIntflw = ThickCum - (-ZBotDr(NumDrain))
               ThickCompBotIntflw = ThickComp(NumComBotIntflw) -            &
         &                        ThickCompUndIntflw
               TransmisIntflw = TransmisIntflw -                        &
         &                        CondSatHor(NumComBotIntflw)*             &
         &                        ThickCompUndIntflw
               if (NumComBotIntflw .eq. NumComWatLev) then
                  ThickCompSatWatLev = ThickCompSatWatLev - ThickCompUndIntflw
                  ThickCompBotIntflw = ThickCompSatWatLev
               end if

! --- 3b distribute interflow fluxes as lateral fluxes over compartments
! --- of the interflow discharge layer
               idr = NumDrain
               FluxDrComp(idr, NumComWatLev) = FluxDr(idr)*                 &
         &      ThickCompSatWatLev*CondSatHor(NumComWatLev)/TransmisIntflw
               do icomp = NumComWatLev + 1, NumComBotIntflw - 1
                  FluxDrComp(idr, icomp) = FluxDr(idr)*ThickComp(icomp)*    &
           &                              CondSatHor(icomp)/TransmisIntflw
               end do
               FluxDrComp(idr, NumComBotIntflw) = FluxDr(idr)*              &
         &                      ThickCompBotIntflw*                        &
         &                      CondSatHor(NumComBotIntflw)/TransmisIntflw

! --- 3c reset variables for determinig top of all other model_discharge_layers
               WatLevAv = -ZBotDr(NumDrain)
               NumComWatLev = NumComBotIntflw
               ThickCompSatWatLev = ThickCompUndIntflw
            end if
! --- no other drainage levels besides interflow -> leave DIVDRA
            if (NumDrHlp .eq. 0) return
         end if

! --- paragraph 4
! --- Calculate overall anisotropic factor model profile and
! --- cumulative transmissivity as a function of depth
         CumThickCondHor = ThickCompSatWatLev*CondSatHor(NumComWatLev)
         CumThickCondVer = ThickCompSatWatLev/CondSatVer(NumComWatLev)
         ThickCum = ThickCompSatWatLev

         do icomp = NumComWatLev + 1, NumComp
            CumThickCondHor = CumThickCondHor + ThickComp(icomp)*         &
        &                       CondSatHor(icomp)
            CumThickCondVer = CumThickCondVer + ThickComp(icomp)/         &
        &                       CondSatVer(icomp)
            ThickCum = ThickCum + ThickComp(icomp)
         end do

         CondSatHorAv = CumThickCondHor/ThickCum
         CondSatVerAv = ThickCum/CumThickCondVer
         FacAniso = dsqrt(CondSatVerAv/CondSatHorAv)

! --- paragraph 5
! --- Calculate maximum depth per drainage system and discharge
! --- flow rate per drainage system
         do idr = 1, NumDrHlp
!                                                                       !   Divdra, infiltration
! --- In case of infiltration and switch for seperate infiltration flux distribution
!     on: set factor for adapting depth of discharge layer to depth of infiltration layer
            if (flDivdInf .and. FluxDr(idr) .lt. -small) then
!   - Constrain infiltration depth to minimum = bottom compartment with groundwaterlevel
               FDisInfmin = ThickCompSatWatLev/                         &
        &                     (0.25d0*DistDrain(idr)*FacAniso)
               FDisInf(idr) = Max(FacDpthInf, FDisInfmin)
            else
               FDisInf(idr) = 1.d0
            end if                                                          !  Divdra, infiltration
!
            MaxDepthDislay(idr) = FDisInf(idr)*                           &  !  Divdra, infiltration
        &                         0.25d0*DistDrain(idr)*FacAniso + WatLevAv
!
! --- Constrain to SWAP profile thickness
            MaxDepthDislay(idr) = MIN(MaxDepthDislay(idr),                &
        &                              Thickcum + WatLevAv)
         end do

! --- paragraph 6
! --- determine sequence of order of drainage systems

         NumActDrain = 0
         do idr = 1, NumDrHlp
            if (abs(FluxDr(idr)) .gt. small) then
               NumActDrain = NumActDrain + 1
               DrainSequence(NumActDrain) = idr
               HelpFl(NumActDrain) = FDisInf(idr)*DistDrain(idr)    !  Divdra, infiltration
            end if
         end do

         do idr = 1, NumActDrain - 1
            do iidr = idr + 1, NumActDrain
               if (HelpFl(idr) .lt. HelpFl(iidr)) then
                  dum1 = HelpFl(iidr)
                  HelpFl(iidr) = HelpFl(idr)
                  HelpFl(idr) = dum1
                  idum = DrainSequence(iidr)
                  DrainSequence(iidr) = DrainSequence(idr)
                  DrainSequence(idr) = idum
               end if
            end do
         end do

         idr = DrainSequence(NumActDrain)
         FlowDrDisch(idr) = FDisInf(idr)*abs(FluxDr(idr)*DistDrain(idr))   !  Divdra, infiltration
         do iidr = NumActDrain - 1, 1, -1
            idr = DrainSequence(iidr)
            jdr = DrainSequence(iidr + 1)
            FlowDrDisch(idr) = FlowDrDisch(jdr) +                          &
        &                      FDisInf(idr)*abs(FluxDr(idr)*DistDrain(idr)) !  Divdra, infiltration
         end do

! --- paragraph 7
! --- Bottom of 1st order model_discharge_layer

         idr = DrainSequence(1)
         Transmissivity(idr) = CumThickCondHor
         BotDisLay(idr) = ThickCum + WatLevAv
         NumComBotDislay(idr) = NumComp
         ThickCompBotDislay(idr) = ThickComp(NumComp)

! --- correction of BotDisLay(1) if D1 < 0.25 L \/(kv/kh)
         if (abs(FluxDr(idr)) .gt. Small) then
            if (BotDisLay(idr) .gt. MaxDepthDislay(idr)) then
               BotDisLay(idr) = MaxDepthDislay(idr)
! ---       determine adjusted transmissivity and compartment
! ---       number which contains bottom of discharge layer
               NumComBotDislay(idr) = NumComWatLev
               HelpTh = ThickCompSatWatLev
               Transmissivity(idr) = ThickCompSatWatLev*                 &
        &                             CondSatHor(NumComWatLev)
               ThickDislay = BotDisLay(idr) - WatlevAv
               do while (ThickDislay .gt. HelpTh)
                  NumComBotDislay(idr) = NumComBotDislay(idr) + 1
                  HelpTh = HelpTh +                          &
        &                                ThickComp(NumComBotDislay(idr))
                  Transmissivity(idr) = Transmissivity(idr) +             &
        &                                ThickComp(NumComBotDislay(idr))* &
        &                                CondSatHor(NumComBotDislay(idr))
               end do

               Transmissivity(idr) = Transmissivity(idr) -                 &
        &                          (HelpTh - ThickDislay)*                &
        &                          CondSatHor(NumComBotDislay(idr))
! ---       thickness of part of bottom compartment
               ThickCompBotDislay(idr) = ThickComp(NumComBotDislay(idr)) - &
        &                              (HelpTh - ThickDislay)
            end if
         else
            BotDisLay(idr) = WatLevAv
         end if

! --- paragraph 8
! --- Bottom of 2nd and higher order model_discharge_layers
! --- for drainage system of orders 2 to NumDrain (number of drains)
         do iidr = 2, NumActDrain
            idr = DrainSequence(iidr)
            jdr = DrainSequence(iidr - 1)
            Transmissivity(idr) = Transmissivity(jdr)*                    &
        &                         FlowDrDisch(idr)/FlowDrDisch(jdr)

! --- bottom of discharge layer of order i
!     and thickness of bottom compartment, and
!     number of bottom compartment
            icomp = NumComWatLev
            HelpTr = ThickCompSatWatLev*CondSatHor(icomp)
            HelpTh = ThickCompSatWatLev
            do while (HelpTr .lt. transmissivity(idr))
               icomp = icomp + 1
               HelpTr = HelpTr + ThickComp(icomp)*CondSatHor(icomp)
               HelpTh = HelpTh + ThickComp(icomp)
            end do
            NumComBotDislay(idr) = icomp
            ThickCompBotDislay(idr) = ThickComp(icomp) -                   &
        &             (HelpTr - Transmissivity(idr))/CondSatHor(icomp)
            BotDisLay(idr) = WatLevAv + HelpTh +                   &
     &             (Transmissivity(idr) - HelpTr)/CondSatHor(icomp)

! --- paragraph 9
!     correction of BotDisLay(idr) if D(idr) < 0.25 L(idr) \/(kv/kh)
            if (BotDisLay(idr) .gt. MaxDepthDislay(idr)) then
               BotDisLay(idr) = MaxDepthDislay(idr)

! --- transmissivity of order i, thickness
!     of bottom compartment, and number
!     of bottom compartment
               NumComBotDislay(idr) = NumComWatLev
               HelpTh = ThickCompSatWatLev
               Transmissivity(idr) = ThickCompSatWatLev*                 &
        &                             CondSatHor(NumComWatLev)
               ThickDislay = BotDisLay(idr) - WatlevAv
               do while (ThickDislay .gt. HelpTh)
                  NumComBotDislay(idr) = NumComBotDislay(idr) + 1
                  HelpTh = HelpTh + ThickComp(NumComBotDislay(idr))
                  Transmissivity(idr) = Transmissivity(idr) +             &
        &                               ThickComp(NumComBotDislay(idr))*  &
        &                               CondSatHor(NumComBotDislay(idr))
               end do
               Transmissivity(idr) = Transmissivity(idr) -                 &
        &                            (HelpTh - ThickDislay)*              &
        &                            CondSatHor(NumComBotDislay(idr))
               ThickCompBotDislay(idr) = ThickComp(NumComBotDislay(idr)) - &
        &                            (HelpTh - ThickDislay)
            end if
         end do

! --- paragraph 10
!     Distribute drainage fluxes as lateral fluxes over
!     ith order model discharge layers
         do iidr = 1, NumActDrain
            idr = DrainSequence(iidr)

            if (.not. flDivdInf .or. FluxDr(idr) .gt. 0.d0) then              !  Divdra, infiltration
               FluxDrComp(idr, NumComWatLev) = FluxDr(idr)*ThickCompSatWatLev* &
         &                               CondSatHor(NumComWatLev)/         &
         &                               Transmissivity(idr)
               do icomp = NumComWatLev + 1, NumComBotDislay(idr) - 1
                  FluxDrComp(idr, icomp) = FluxDr(idr)*ThickComp(icomp)*     &
         &                                CondSatHor(icomp)/               &
         &                                Transmissivity(idr)
               end do
               FluxDrComp(idr, NumComBotDislay(idr)) = FluxDr(idr)*          &
         &                                ThickCompBotDislay(idr)*          &
         &                                CondSatHor(NumComBotDislay(idr))/ &
         &                                Transmissivity(idr)
               if (NumComWatLev .eq. NumComBotDislay(idr))                  &
         &         FluxDrComp(idr, NumComBotDislay(idr)) = FluxDr(idr)
            end if
         end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  Divdra, Infiltration !!!!!!!!!!!!!!!!!!!!!!!!
!
!     If flDivdInf = true, then distribution of infiltration fluxes is separate
!     from distribution discharge fluxes.
         if (flDivdInf) then

!        do idr = 1, NumActDrain
            do iidr = 1, NumActDrain
               idr = DrainSequence(iidr)
               if (FluxDr(idr) .lt. -small) then

! --- Find surface water level
!     - first copy to 1-dimensional table temptab
                  !do i = 1,2*maowl
                  !   temptab(i) = owltab(idr,i)
                  !end do
                  !SurfLev = afgen(temptab,2*maowl,t1900+dt-1.d0)
                  SurfLev = afgen(owltab(idr, 1:2*nowltab(idr)), 2*nowltab(idr), t1900 + dt - 1.d0)
! --- Search for compartment with surface water level: NumComSrfLev
!   - part of NumComSrfLev below SurfLev
                  call Lev2Comp(NumComp, -1.0d0*Min(SurfLev, 0.0d0), ThickComp, &
         &                      NumComSrfLev, Dum1, Dum2, ThickCompBlwSrfLev)

! --- Depth of Infiltration layer: DepthInflay
                  DpthInflay = MaxDepthDislay(idr)
! --- Search for compartment with depth Infiltration layer: NumComDepthInf
!   - part of NumComDepthInf Above DepthInflay
                  call Lev2Comp(NumComp, DpthInflay, ThickComp,                &
         &                      NumComDpthInf, Dum1, ThickCompAbvDpthInf, Dum2)

! --- Special case: DpthInflay = bottom compartment with groundwaterlevel
                  if (NumComDpthInf .eq. NumComWatLev) then
                     ThickCompAbvDpthInf = ThickCompSatWatLev
                  end if

! --- Determine transmissivity of total infiltration zone, of its
!     unsaturated part and saturated part, and of the compartments within
                  KD(NumComDpthInf) = ThickCompAbvDpthInf*                 &
         &                            CondSatHor(NumComDpthInf)
                  TransmisTot = KD(NumComDpthInf)
                  TransmisSat = KD(NumComDpthInf)
                  do icomp = NumComDpthInf - 1, NumComSrfLev + 1, -1
                     KD(icomp) = ThickComp(icomp)*CondSatHor(icomp)
                     TransmisTot = TransmisTot + KD(icomp)
                     if (icomp .gt. NumComWatLev) then
                        TransmisSat = TransmisSat + KD(icomp)
                     end if
                  end do
                  KD(NumComSrfLev) = ThickCompBlwSrfLev*                   &
         &                           CondSatHor(NumComSrfLev)
                  if (NumComDpthInf .eq. NumComWatLev) then
!   - special case: DpthInflay = bottom compartment with groundwaterlevel
                     if (NumComSrfLev .eq. NumComWatLev) then
                        KD(NumComSrfLev) = KD(NumComSrfLev) - TransmisTot
                     else
                        TransmisTot = TransmisTot + ThickCompUnsWatLev*    &
         &                            CondSatHor(NumComWatLev)
                     end if
                  end if
!
                  TransmisTot = TransmisTot + KD(NumComSrfLev)
                  KD(NumComWatLev) = ThickCompSatWatLev*                   &
         &                           CondSatHor(NumComWatLev)
!   Bug Fixed, RobH 28-1-2018
                  if (NumComDpthInf .gt. NumComWatLev) then
                     TransmisSat = TransmisSat + KD(NumComWatLev)
                     TransmisTot = TransmisTot + KD(NumComWatLev)
                  end if
                  TransmisUns = TransmisTot - TransmisSat

! --- RQmax = max flux per unit of transmissivity
                  RQmax = 2.d0*FluxDr(idr)/TransmisTot

! --- Distribute infiltration fluxes as lateral fluxes over unsaturated compartments
                  if (NumComSrfLev .eq. NumComWatLev) then
                     FluxComWatLev = RQmax*                                &
         &              0.5d0*(ThickCompBlwSrfLev - ThickCompSatWatLev)* &
         &              CondSatHor(NumComWatLev)
                  else
                     if (TransmisUns .gt. 1.d-8) then
                        CumKD = 0.d0
                        do icomp = NumComSrfLev, NumComWatLev - 1
                           FluxDrComp(idr, icomp) = RQmax/TransmisUns*0.5d0 &
                                                   &                            *((CumKD + KD(icomp))**2 - CumKD**2)
                           CumKD = CumKD + KD(icomp)
                        end do
                        FluxComWatLev = RQmax/TransmisUns*0.5d0*((CumKD  &
         &                + ThickCompUnsWatLev*CondSatHor(NumComWatLev))**2    &
         &                                                   - CumKD**2)
                     else
                        FluxComWatLev = 0.0d0
                     end if
                  end if

! --- Distribute infiltration fluxes as lateral fluxes over saturated compartments
                  CumKD = 0.d0
                  do icomp = NumComDpthInf, NumComWatLev, -1
                     FluxDrComp(idr, icomp) = RQmax/TransmisSat*          &
         &                 0.5d0*((CumKD + KD(icomp))**2 - CumKD**2)
                     CumKD = CumKD + KD(icomp)
                  end do
! --- Add unsaturated flux to infiltration flux of compartment containing groundwater level
                  FluxDrComp(idr, NumComWatLev) =                            &
         &                      FluxDrComp(idr, NumComWatLev) + FluxComWatLev

               end if
            end do

         end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  Divdra, Infiltration !!!!!!!!!!!!!!!!!!!!!!!!

         return
         end subroutine divdra

! File VersionID:
!   $Id: divdra.f90 370 2018-02-09 13:29:10Z heine003 $
! ----------------------------------------------------------------------
         SUBROUTINE Lev2Comp(NumComp, Level, ThickComp, NumCom2Lev,  &
        &                    ThickCum, ThickCompAbvLev, ThickCompBlwLev)
    !! Find compartment containing a specified depth level
    !!
    !! This helper subroutine locates which soil compartment contains a given
    !! depth level (measured from the soil surface) and calculates how the
    !! compartment is split by that level.
    !!
    !! The routine accumulates compartment thicknesses from top to bottom until
    !! the specified level is reached, then computes the portions of that
    !! compartment lying above and below the level.
    !!
    !!@note
    !! Depths are measured positive downward from the soil surface (cm below surface).
    !!@endnote
    !!
    !!@warning
    !! If the specified level exceeds the total soil profile depth, the routine
    !! calls `fatalerr` to terminate execution.
    !!@endwarning
    !!@note
    !! ----------------------------------------------------------------------
    !!     Date               : 29/10/10
    !!     Purpose            : Find Compartment with Level and its part
    !!                          Above and Under this Level
    !!     Subroutines called : -
    !!     Functions called   : -
    !!     File usage         : -
    !! ----------------------------------------------------------------------
    !!@endnote
        use swap_array_dimensions, only: macp

            integer NumCom2Lev, NumComp

            real(8) Level, ThickComp(Macp), ThickCompAbvLev, ThickCompBlwLev
            real(8) ThickCum

! --- Find number Compartment containing Level
            NumCom2Lev = 1
            ThickCum = ThickComp(1)
            do while (Level .gt. ThickCum + 1.d-10)
               NumCom2Lev = NumCom2Lev + 1
               ThickCum = ThickCum + ThickComp(NumCom2Lev)
            end do
! RobH 28-1-2018: Extra (most likely unnecessary) protection against possible level below soil profile
            if (NumCom2Lev .gt. NumComp) call fatalerr_collected('Lev2Comp', 'NumCom2Lev.gt.Numcomp')

! --- part of compartment Below level
            ThickCompBlwLev = ThickCum - Level

! --- part of compartment Above level
            ThickCompAbvLev = ThickComp(NumCom2Lev) - ThickCompBlwLev

            return
         end subroutine Lev2Comp
         end module distribute_drainage
