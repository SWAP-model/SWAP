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

  public :: VonHHBraden, Gash, DivIntercep

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
