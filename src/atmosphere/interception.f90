!> @brief Rainfall interception calculations using various methods
!>
!> This module provides different interception calculation methods for simulating
!> the interception of rainfall and irrigation water by vegetation canopy. It includes
!> methods ranging from simple empirical relationships to complex analytical solutions
!> based on canopy storage and evaporation dynamics.
!>
!> Available methods:
!> - Von Hoyningen-Hune and Braden: Exponential relationship between LAI and interception
!> - Gash (1995): Analytical model for forest interception
!> - Adapted Rutter method: Sparse Gash concept for detailed canopy storage simulation
!>
!> @author Original SWAP team
!> @date Last modified: July 2012
module interception_mod

  use swap_state_mod, only: swap_state_t

  public :: VonHHBraden, Gash, ruttervw, msw1eic, DivIntercep

contains

  !> @brief Calculate interception using Von Hoyningen-Hune and Braden method
  !>
  !> This method uses an exponential relationship between leaf area index (LAI)
  !> and soil cover to calculate rainfall interception. The intercepted water
  !> is assumed to evaporate directly from the canopy.
  !>
  !> @param[out] aintc Amount of rainfall interception during current day [cm/d]
  !>
  !> @note
  !> Last modified: July 2012
  !>
  !> Uses exponential relation between soil cover and LAI.
  !>
  !> Pure function — explicit args, no state dependency.
  !> @endnote
  pure function VonHHBraden(grai, gird, isua, kdif, kdir, lai, cofab) result(aintc)
    implicit none

    real(8), intent(in) :: grai   ! Gross daily rain flux [cm/d]
    real(8), intent(in) :: gird   ! Gross daily irrigation flux [cm/d]
    integer, intent(in) :: isua   ! Sprinkler switch: 0 = above canopy, /=0 = below
    real(8), intent(in) :: kdif   ! Diffuse-light extinction coefficient [-]
    real(8), intent(in) :: kdir   ! Direct-light extinction coefficient [-]
    real(8), intent(in) :: lai    ! Leaf area index [m2/m2]
    real(8), intent(in) :: cofab  ! Von Hoyningen-Hune & Braden interception coefficient [cm/d]
    real(8)             :: aintc  ! Rainfall interception this day [cm/d]

    real(8) :: rpd    ! Intercepted precipitation (rain+irrig) [mm]
    real(8) :: cofbb  ! Interception coefficient b [-]

    ! Intercepted precipitation (rain+irrig) in mm
    rpd = grai*10.0d0
    if (isua .eq. 0) rpd = (grai + gird)*10.0d0

    ! Exponential relation between soil cover and lai
    cofbb = 1.0d0 - exp(-1.0d0*kdif*kdir*lai)
    cofbb = min(cofbb, 1.0d0)

    ! Interception: evaporation of intercepted precipitation in cm
    if (cofab .gt. 0.000001d0) then
      aintc = (cofab*lai*(1.0d0 - (1/(1.0d0 + rpd*cofbb / (cofab*lai)))))*0.1d0
    else
      aintc = 0.0d0
    endif

  end function VonHHBraden

  !> @brief Calculate interception for forests according to Gash (1995)
  !>
  !> Implements the analytical Gash model for forest rainfall interception.
  !> The model distinguishes between different phases of interception (wetting,
  !> saturation, and drying) and accounts for free throughfall and stem flow.
  !>
  !> @param[out] aintc Amount of rainfall interception during current day [cm/d]
  !>
  !> @note
  !> Last modified: July 2012
  !>
  !> Reference: Gash, J.H.C. (1995). An analytical framework for estimating
  !> evaporation using rainfall and forest data.
  !>
  !> Pure function — explicit args, no state dependency. The 5 AFGEN
  !> table lookups (pfree/pstem/scanopy/avprec/avevap raw values) are
  !> performed by the caller and passed in as scalars.
  !> @endnote
  pure function Gash(grai, gird, isua, pfree, pstem, scanopy_raw, avprec_raw, avevap_raw) result(aintc)
    implicit none

    real(8), intent(in) :: grai         ! Gross daily rain flux [cm/d]
    real(8), intent(in) :: gird         ! Gross daily irrigation flux [cm/d]
    integer, intent(in) :: isua         ! Sprinkler switch: 0 = above canopy, /=0 = below
    real(8), intent(in) :: pfree        ! Free throughfall coefficient [-] — AFGEN(pfreetb, t)
    real(8), intent(in) :: pstem        ! Stem flow coefficient [-]        — AFGEN(pstemtb, t)
    real(8), intent(in) :: scanopy_raw  ! Canopy storage capacity [cm]     — AFGEN(scanopytb, t), pre-divide
    real(8), intent(in) :: avprec_raw   ! Average rainfall intensity [cm/d] — AFGEN(avprectb, t)
    real(8), intent(in) :: avevap_raw   ! Average evaporation intensity [cm/d] — AFGEN(avevaptb, t), pre-divide
    real(8)             :: aintc        ! Rainfall interception this day [cm/d]

    real(8) :: avevap   ! Average evaporation intensity, /cGash [-]
    real(8) :: cGash    ! Slope of dPi/dPgross before saturation of canopy [-]
    real(8) :: psatcan  ! Amount of rainfall to saturate canopy [cm]
    real(8) :: rpd      ! Intercepted precipitation (rain+irrig) [cm]
    real(8) :: scanopy  ! Storage capacity of canopy, /cGash [cm]

    ! Intercepted precipitation (rain+irrig) in cm
    if (isua .eq. 0) then
      rpd = grai + gird
    else
      rpd = grai
    endif

    ! Sparse-Gash sparse-canopy scaling: covered-fraction inverse.
    cGash   = 1.d0 - pfree - pstem
    scanopy = scanopy_raw / cGash
    avevap  = avevap_raw  / cGash

    ! Amount of rainfall to saturate canopy
    if ((1.0d0 - avevap/avprec_raw) .gt. 1.0d-4) then
      psatcan = -avprec_raw*scanopy/avevap * log(1.0d0 - avevap/avprec_raw)
    else
      psatcan = avprec_raw*scanopy/avevap
    endif

    ! Interception: evaporation of intercepted precipitation in cm
    if (grai .lt. psatcan) then
      aintc = cGash * rpd
    else
      aintc = cGash * (psatcan + avevap*cGash / avprec_raw * (rpd - psatcan))
    endif

  end function Gash

  !> @brief Simulate interception using adapted Rutter method (wrapper for msw1eic)
  !>
  !> This subroutine serves as a wrapper for the msw1eic routine, which implements
  !> the Sparse Gash analytical model. It converts between SWAP and MetaSWAP data
  !> structures and handles the interception storage dynamics on the vegetation canopy.
  !>
  !> @param[in]  gctp  Soil cover [-]
  !> @param[out] aintc Intercepted rainfall [cm/d]
  !> @param[out] eintc Interception evaporation [cm/d]
  !>
  !> @note
  !> Author: Paul van Walsum
  !>
  !> Date: 08/06/2012
  !>
  !> Adapted Rutter method of Van Walsum & Supit (2012).
  !>
  !> Explicit-args subroutine — no state dependency. Not pure: msw1eic
  !> (item #4 quarantine) does `write(unit) + stop` on its single error
  !> path. Wrapping pack-from-state happens at the call site.
  !> @endnote
  subroutine ruttervw(gctp, dt, siccapact, fimin, ew0, grai, sicact, aintc, eintc)
    use swap_log, only: log_unit_handle
    implicit none

    real(8), intent(in)    :: gctp       ! Soil cover [-]
    real(8), intent(in)    :: dt         ! Time step [d]
    real(8), intent(in)    :: siccapact  ! Active canopy storage capacity [cm]
    real(8), intent(in)    :: fimin      ! Minimum relative canopy evaporation factor [-]
    real(8), intent(in)    :: ew0        ! Reference wet-canopy evaporation rate [mm/d]
    real(8), intent(in)    :: grai       ! Gross daily rain flux [cm/d]
    real(8), intent(inout) :: sicact     ! Canopy storage [cm] — read on entry, updated on exit
    real(8), intent(out)   :: aintc      ! Intercepted rainfall [cm/d]
    real(8), intent(out)   :: eintc      ! Interception evaporation [cm/d]

    ! Local conversion buffers for the (real(4), metaswap-compatible) msw1eic
    integer(4) :: nuk_i4, ibd_i4(1), ib_i4
    real(4)    :: dc_r4, dtsw_r4, csk_r4(1), vxick_r4(1), fecmnk_r4(1)
    real(4)    :: ETw0_r4(1), Pgdtsw_r4(1), Sic_r4(1), Sicolddtsw_r4(1)
    real(4)    :: Picdtsw_r4(1), Eicdtsw_r4(1), tcap_r4(1), beta_r4(1)
    real(4)    :: zeta_r4(1), fricdtsw_r4(1)

    nuk_i4       = 1
    ibd_i4(1)    = 1
    ib_i4        = log_unit_handle()  ! quarantined msw1eic still expects a unit number
    dc_r4        = 1.0e-4
    dtsw_r4      = real(dt)
    csk_r4(1)    = real(gctp)
    vxick_r4(1)  = real(siccapact)
    fecmnk_r4(1) = real(fimin)
    ETw0_r4(1)   = real(ew0*0.1d0)
    Pgdtsw_r4(1) = real(grai)
    Sic_r4(1)    = real(sicact)

    call msw1eic(nuk_i4, ibd_i4, dc_r4, dtsw_r4, csk_r4, vxick_r4, &
                 fecmnk_r4, ETw0_r4, Pgdtsw_r4, Sic_r4, Sicolddtsw_r4, Picdtsw_r4, &
                 Eicdtsw_r4, tcap_r4, beta_r4, zeta_r4, fricdtsw_r4, ib_i4)

    sicact = dble(Sic_r4(1))
    aintc  = dble(Picdtsw_r4(1))
    eintc  = dble(Eicdtsw_r4(1))

  end subroutine ruttervw

  !> @brief Interception simulation with Sparse Gash concept
  !>
  !> This subroutine implements a detailed analytical solution for canopy interception
  !> based on the Sparse Gash concept. It accounts for the saturation degree of the canopy
  !> and solves a linear differential equation to track interception storage and evaporation
  !> dynamics during rainfall events.
  !>
  !> The model distinguishes between periods when the canopy is filling, saturated, or drying,
  !> and calculates the appropriate evaporation rate for each phase.
  !>
  !> @param[in]    nuk        Number of SVATs (Soil Vegetation Atmosphere Transfer units) [-]
  !> @param[in]    ibd        Existence flag for SVAT (0=inactive, 1=active) [-]
  !> @param[in]    dc         Near-zero real value for numerical comparisons [-]
  !> @param[in]    dtsw       Time step [d]
  !> @param[in]    csk        Soil cover fraction [m²/m²]
  !> @param[in]    vxick      Interception capacity of canopy [cm]
  !> @param[in]    fecmnk     Minimum relative canopy evaporation factor [-]
  !> @param[in]    ETw0       Evaporation rate from wet canopy [cm/d]
  !> @param[in]    Pgdtsw     Gross rainfall plus sprinkling irrigation [cm]
  !> @param[inout] Sic        Interception storage of SVAT [cm]
  !> @param[inout] Sicolddtsw Interception storage at start of time step [cm]
  !> @param[out]   Picdtsw    Intercepted precipitation during time step [cm]
  !> @param[out]   Eicdtsw    Interception evaporation during time step [cm]
  !> @param[out]   tcap       Time to fill interception reservoir [d]
  !> @param[out]   beta       Coefficient of differential equation [1/d]
  !> @param[out]   zeta       Coefficient of differential equation [cm/d]
  !> @param[out]   fricdtsw   Fraction of time used by interception evaporation [-]
  !> @param[in]    ib         Unit number of log file for error messages [-]
  !>
  !> @note
  !> Copyright: 2009 Alterra
  !>
  !> File: MSW1EIC.FOR
  !>
  !> This program, or parts thereof, may not be reproduced, modified or transferred
  !> to third parties without written permission.
  !>
  !> Modified for relationship with saturation degree of canopy.
  !> @endnote
  subroutine msw1eic(nuk,ibd,dc,dtsw,csk,vxick,fecmnk,ETw0,Pgdtsw, &
                     Sic,Sicolddtsw,Picdtsw,Eicdtsw,tcap,beta,zeta, &
                     fricdtsw,ib)
    implicit none

    ! Arguments
    integer(4), intent(in)    :: nuk
    integer(4), intent(in)    :: ibd(1)
    real(4),    intent(in)    :: dc
    real(4),    intent(in)    :: dtsw
    real(4),    intent(in)    :: csk(1)
    real(4),    intent(in)    :: vxick(1)
    real(4),    intent(in)    :: fecmnk(1)
    real(4),    intent(in)    :: ETw0(1)
    real(4),    intent(in)    :: Pgdtsw(1)
    real(4),    intent(inout) :: Sic(1)
    real(4),    intent(inout) :: Sicolddtsw(1)
    real(4),    intent(out)   :: Picdtsw(1)
    real(4),    intent(out)   :: Eicdtsw(1)
    real(4),    intent(out)   :: tcap(1)
    real(4),    intent(out)   :: beta(1)
    real(4),    intent(out)   :: zeta(1)
    real(4),    intent(out)   :: fricdtsw(1)
    integer(4), intent(in)    :: ib

    ! Local variables
    integer(4) :: k                ! Index for SVATs

!$OMP PARALLEL DO
!$OMP&  DEFAULT(SHARED)
!$OMP&  PRIVATE(k)
    do k=1,nuk
      if (ibd(k) .ge. 1) then

        ! Check that non-zero interception capacity has non-zero soil cover
        if (vxick(k) .gt. dc) then
          if (csk(k) .lt. dc) then
            write(ib,9199) k
            write(*,9199) k
            stop
          endif
        endif
9199    format(' Interception capacity >0, but soil cover = 0, k =',i10)

        ! Check first for shortcut to zero interception storage
        if ( vxick(k) .lt. dc .or. &
             (Sic(k) .lt. dc .and. Pgdtsw(k) .lt. dc) ) then

          ! Interception capacity zero, or nothing doing
          if (Sic(k) .gt. dc) then
            ! Remaining interception water assumed to evaporate as vegetation dies off
            Eicdtsw(k) = Sic(k)/dtsw
            Sic(k)     = 0.
          else
            Eicdtsw(k) = 0.
          endif
          Picdtsw(k)   = 0.
        else

          ! Starting point of complete calculation
          Sicolddtsw(k) = Sic(k)

          ! Check first for shortcut of full reservoir that stays full
          if ( Sicolddtsw(k) .ge. (vxick(k) - 2*dc) .and. &
               (csk(k)*Pgdtsw(k)) .ge. ETw0(k) ) then

            ! Reservoir starts full and stays full
            Sic(k)      = vxick(k)
            Eicdtsw(k)  = ETw0(k)  ! csk implicit in faeic
          else

            ! First check for beta=0 in differential equation
            if (ETw0(k) .lt. dc .or. fecmnk(k) .gt. 0.99999) then

              ! Beta=0 (see below for expression)
              Sic(k) = Sicolddtsw(k) + (csk(k)*Pgdtsw(k) - ETw0(k))*dtsw
              Sic(k) = min(Sic(k),vxick(k))
              Sic(k) = max(Sic(k),0.)
              if (Sic(k) .lt. vxick(k)) then

                ! Reservoir does not become full, can become empty
                ! Evaporation from balance (stops when reservoir empties)
                Eicdtsw(k) = (Sicolddtsw(k) - Sic(k))/dtsw + &
                             csk(k)*Pgdtsw(k)
              else

                ! Reservoir becomes full, evaporation at full rate
                ! (fecmn=1.) or zero (ETw0=0.)
                Eicdtsw(k) = ETw0(k)
              endif
            else

              ! Auxiliary parameters for solving linear differential equation
              beta(k) = (1.0 - fecmnk(k))*ETw0(k)/(vxick(k))
              zeta(k) = csk(k)*Pgdtsw(k) - fecmnk(k)*ETw0(k)

              ! First calculation of new Sic is tentative
              Sic(k) = (Sicolddtsw(k) - zeta(k)/beta(k))* &
                       exp(-beta(k)*dtsw) + zeta(k)/beta(k)

              ! Value can be negative if formula didn't account for
              ! evaporation rate going to zero when Sic goes to zero
              Sic(k) = max(Sic(k),0.)

              ! Final calculation of new Sic and Eictdtsw
              if (Sic(k) .lt. vxick(k)) then

                ! Reservoir does not become full, Ec from balance
                Eicdtsw(k) = (Sicolddtsw(k)-Sic(k))/dtsw + &
                             csk(k)*Pgdtsw(k)
              else

                ! Reservoir becomes full; find out when: tcap
                Sic(k)  = vxick(k)
                tcap(k) = (1./beta(k))* &
                          log((Sicolddtsw(k)-zeta(k)/beta(k))/ &
                              (vxick(k) - zeta(k)/beta(k)))

                ! From t=0 to tcap: Eic from balance; rest: potential rate
                Eicdtsw(k) = (1./dtsw)*( Sicolddtsw(k) - Sic(k) + &
                                        csk(k)*Pgdtsw(k)*tcap(k) + &
                                        ETw0(k)*(dtsw - tcap(k)) )
              endif
            endif

          endif

          ! Intercepted precipitation from balance
          Picdtsw(k) = (Sic(k)-Sicolddtsw(k))/dtsw + Eicdtsw(k)
        endif

        ! Fraction of time that interception evaporation is active
        if (ETw0(k) .gt. dc) then
          fricdtsw(k) = Eicdtsw(k)/ETw0(k)
        else
          fricdtsw(k) = 0.
        endif

      endif
    enddo
!$OMP END PARALLEL DO

  end subroutine msw1eic

  !> @brief Divide interception into rain and irrigation parts
  !>
  !> This subroutine partitions the total interception amount between rainfall
  !> and irrigation water, and calculates the net amounts that reach the soil surface.
  !> It accounts for snow and handles cases with and without separate irrigation tracking.
  !>
  !> @param[in] aintc Total interception amount [cm/d]
  !>
  !> @note
  !> Last modified: February 2014
  !>
  !> Input via state: state%atmosphere%isua/grai/gsnow/snrai, state%crop%gird
  !> Writes: state%atmosphere%nraida, state%atmosphere%nird
  !> @endnote
  pure subroutine DivIntercep (aintc, state)
    implicit none

    real(8),            intent(in)    :: aintc   ! Total interception [cm/d]
    type(swap_state_t), intent(inout) :: state

    associate (atmo => state%atmosphere, crop => state%crop)

      ! Divide interception into rain and irrigation parts;
      ! compute net rain (nraida) and net sprinkling irrigation (nird).
      if (aintc .lt. 0.001d0) then
        atmo%nraida = atmo%grai - atmo%gsnow - atmo%snrai
        atmo%nird   = crop%gird
      else
        if (atmo%isua .eq. 0) then
          atmo%nraida = atmo%grai - aintc*(atmo%grai/(atmo%grai+crop%gird))
          atmo%nird   = crop%gird - aintc*(crop%gird/(atmo%grai+crop%gird))
        else
          atmo%nraida = atmo%grai - aintc
          atmo%nird   = crop%gird
        endif
      endif

    end associate

  end subroutine DivIntercep

end module interception_mod
