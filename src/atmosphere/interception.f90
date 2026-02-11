!> Rainfall interception calculations using various methods
! Subroutines:
! - VonHHBraden      ! Von Hoyningen-Hune and Braden method
! - Gash             ! Gash (1995) forest interception method
! - ruttervw         ! Adapted Rutter method wrapper
! - msw1eic          ! Sparse Gash concept implementation
! - DivIntercep      ! Divide interception into rain/irrigation parts
module interception_mod

    public :: VonHHBraden, Gash, ruttervw, msw1eic, DivIntercep

contains

  !> Calculate interception using Von Hoyningen-Hune and Braden method
  !>
  !> @param[out] aintc Amount of rainfall interception during current day [cm/d]
  !>
  !> @note
  !> Last modified: July 2012
  !> Uses exponential relation between soil cover and LAI
  !> Input from variables module: grai, gird, kdif, kdir, cofab, lai, isua
  !> @endnote
  subroutine VonHHBraden (aintc)
    use variables, only: grai,gird,kdif,kdir,cofab,lai,isua
    implicit none

    ! Arguments
    real(8), intent(out) :: aintc  ! Amount of rainfall interception during current day [cm/d]

    ! Local variables
    real(8) :: rpd                 ! Intercepted precipitation (rain+irrig) [mm]
    real(8) :: cofbb               ! Interception coefficient b Von Hoyningen-Hune and Braden [-]

    ! Intercepted precipitation (rain+irrig) in mm
    rpd = grai*10.0d0
    if (isua.eq.0) rpd = (grai+gird)*10.0d0

    ! Exponential relation between soil cover and lai
    cofbb = 1.0d0 - dexp(-1.0d0*kdif*kdir*lai)
    cofbb = min(cofbb,1.0d0)

    ! Interception: evaporation of intercepted precipitation in cm
    if (cofab.gt.0.000001d0) then
      aintc = (cofab*lai*(1.0d0-(1/(1.0d0+rpd*cofbb/ &
                                  (cofab*lai)))))*0.1d0
    else
      aintc = 0.0d0
    endif

  end subroutine VonHHBraden

  !> Calculate interception for forests according to Gash (1995)
  !>
  !> @param[out] aintc Amount of rainfall interception during current day [cm/d]
  !>
  !> @note
  !> Last modified: July 2012
  !> Input from variables module: grai, gird, avevaptb, avprectb, pfreetb, pstemtb, scanopytb, isua, t
  !> @endnote
  subroutine Gash (aintc)
    use variables, only: grai,gird,avevaptb,avprectb,pfreetb,pstemtb,scanopytb,isua,t
    use array_utils, only: afgen
    use swap_array_dimensions, only: magrs
    implicit none
    ! include 'arrays.fi'

    ! Arguments
    real(8), intent(out) :: aintc  ! Amount of rainfall interception during current day [cm/d]

    ! Local variables
    real(8) :: avevap              ! Average evaporation intensity during shower [-]
    real(8) :: avprec              ! Average rainfall intensity [-]
    real(8) :: cGash               ! Slope of dPi/dPgross before saturation of canopy [-]
    real(8) :: pfree               ! Free throughfall coefficient [-]
    real(8) :: pstem               ! Stem flow coefficient [-]
    real(8) :: psatcan             ! Amount of rainfall to saturate canopy [cm]
    real(8) :: rpd                 ! Intercepted precipitation (rain+irrig) [cm]
    real(8) :: scanopy             ! Storage capacity of canopy [cm]

    ! Intercepted precipitation (rain+irrig) in cm
    if (isua.eq.0) then
      rpd = grai+gird
    else
      rpd = grai
    endif

    ! Calculate interception for forests according to Gash (1995)
    pfree = afgen(pfreetb,(2*magrs),t)
    pstem = afgen(pstemtb,(2*magrs),t)
    cGash = 1.d0-pfree-pstem
    scanopy = afgen(scanopytb,(2*magrs),t) / cGash
    avprec = afgen(avprectb,(2*magrs),t)
    avevap = afgen(avevaptb,(2*magrs),t) / cGash

    ! Amount of rainfall to saturate canopy
    if ( (1.0d0 - avevap/avprec) .gt. 1.0d-4) then
      psatcan = -avprec*scanopy/avevap * &
                dlog(1.0d0 - avevap/avprec)
    else
      psatcan = avprec*scanopy/avevap
    endif

    ! Interception: evaporation of intercepted precipitation in cm
    if (grai .lt. psatcan) then
      aintc = cGash * rpd
    else
      aintc = cGash * ( psatcan + &
                avevap*cGash / avprec * (rpd - psatcan) )
    endif

  end subroutine Gash

  !> Simulate interception using adapted Rutter method (wrapper for msw1eic)
  !>
  !> @param[in]  gctp  Soil cover [-]
  !> @param[out] aintc Intercepted rainfall [cm/d]
  !> @param[out] eintc Interception evaporation [cm/d]
  !>
  !> @note
  !> Author: Paul van Walsum
  !> Date: 08/06/2012
  !> Adapted Rutter method of Van Walsum & Supit (2012)
  !> Input from variables module: logf, dt, sicact, siccapact, fimin, ew0, grai
  !> sicact (storage on vegetation canopy [cm]) is both input and output via variables module
  !> @endnote
  subroutine ruttervw (gctp,aintc,eintc)
    use variables, only: logf,dt,sicact,siccapact,fimin,ew0,grai
    implicit none

    ! Arguments
    real(8), intent(in)  :: gctp   ! Soil cover [-]
    real(8), intent(out) :: aintc  ! Intercepted rainfall [cm/d]
    real(8), intent(out) :: eintc  ! Interception evaporation [cm/d]

    ! Local variables
    integer(4) :: nuk_i4, ibd_i4(1), ib_i4
    real(4)    :: dc_r4, dtsw_r4, csk_r4(1), vxick_r4(1), fecmnk_r4(1)
    real(4)    :: ETw0_r4(1), Pgdtsw_r4(1), Sic_r4(1), Sicolddtsw_r4(1)
    real(4)    :: Picdtsw_r4(1), Eicdtsw_r4(1), tcap_r4(1), beta_r4(1)
    real(4)    :: zeta_r4(1), fricdtsw_r4(1)

    ! Convert arguments to keep msw1eic routine compatible with metaswap
    nuk_i4       = 1
    ibd_i4(1)    = 1
    ib_i4        = logf
    dc_r4        = 1.0e-4
    dtsw_r4      = REAL(dt)
    csk_r4(1)    = REAL(gctp)
    vxick_r4(1)  = REAL(siccapact)
    fecmnk_r4(1) = REAL(fimin)
    ETw0_r4(1)   = REAL(ew0*0.1d0)
    Pgdtsw_r4(1) = REAL(grai)
    Sic_r4(1)    = REAL(sicact)

    call msw1eic(nuk_i4,ibd_i4,dc_r4,dtsw_r4,csk_r4,vxick_r4, &
                 fecmnk_r4,ETw0_r4,Pgdtsw_r4,Sic_r4,Sicolddtsw_r4,Picdtsw_r4, &
                 Eicdtsw_r4,tcap_r4,beta_r4,zeta_r4,fricdtsw_r4,ib_i4)

    sicact = DBLE(Sic_r4(1))
    aintc  = DBLE(Picdtsw_r4(1))
    eintc  = DBLE(Eicdtsw_r4(1))

  end subroutine ruttervw

  !> Interception simulation with Sparse Gash concept
  !>
  !> Modified for relationship with saturation degree of canopy
  !>
  !> @param[in]    nuk        Number of SVATs [-]
  !> @param[in]    ibd        Existence flag for SVAT (0/1) [-]
  !> @param[in]    dc         Near-zero real [-]
  !> @param[in]    dtsw       Time step [d]
  !> @param[in]    csk        Soil cover [m²/m²]
  !> @param[in]    vxick      Interception capacity of canopy [cm]
  !> @param[in]    fecmnk     Minimum relative canopy evaporation factor [-]
  !> @param[in]    ETw0       Evaporation from wet canopy [cm/d]
  !> @param[in]    Pgdtsw     Gross rainfall + sprinkling [cm]
  !> @param[inout] Sic        Interception storage of SVAT [cm]
  !> @param[inout] Sicolddtsw Interception storage at start of time step [cm]
  !> @param[out]   Picdtsw    Intercepted precipitation [cm]
  !> @param[out]   Eicdtsw    Interception evaporation [cm]
  !> @param[out]   tcap       Time to full interception reservoir [d]
  !> @param[out]   beta       Coefficient of differential equation [1/d]
  !> @param[out]   zeta       Coefficient of differential equation [cm/d]
  !> @param[out]   fricdtsw   Fraction of time used by interception evaporation [-]
  !> @param[in]    ib         Unit number of log file [-]
  !>
  !> @note
  !> Copyright: 2009 Alterra
  !> File: MSW1EIC.FOR
  !> This program, or parts thereof, may not be reproduced, modified or transferred
  !> to third parties without written permission.
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
    integer(4) :: k  ! Index for SVATs

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

  !> Divide interception into rain and irrigation parts
  !>
  !> Calculates net rain and net sprinkling irrigation after interception
  !>
  !> @param[in] aintc Total interception [cm/d]
  !>
  !> @note
  !> Last modified: February 2014
  !> Input from variables module: isua, gird, grai, gsnow, snrai
  !> Output to variables module: nird, nraida
  !> @endnote
  subroutine DivIntercep (aintc)
    use variables, only: isua,gird,grai,gsnow,snrai,nird,nraida
    implicit none

    ! Arguments
    real(8), intent(in) :: aintc  ! Total interception [cm/d]

    ! Divide interception into rain and irrigation parts
    ! and calculate net rain and net sprinkling irrigation
    if (aintc.lt.0.001d0) then
      nraida = grai - gsnow - snrai
      nird = gird
    else
      if (isua.eq.0) then
        nraida = grai-aintc*(grai/(grai+gird))
        nird = gird-aintc*(gird/(grai+gird))
      else
        nraida = grai-aintc
        nird = gird
      endif
    endif

  end subroutine DivIntercep

end module interception_mod
