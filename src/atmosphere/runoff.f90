!> Module for surface runoff calculation using the SCS Curve Number method
!!
!! This module implements the USDA Soil Conservation Service (SCS) Curve Number (CN) method
!! for estimating direct surface runoff from rainfall events. The method accounts for:
!! - Time-varying CN values through lookup tables
!! - Soil moisture corrections (dry, normal, and wet conditions)
!! - Snowmelt contribution to runoff
!!
!! The CN method relates runoff to rainfall through the empirical equation:
!! \[ Q = \frac{(P - I_a)^2}{P - I_a + S} \]
!! where \(Q\) is runoff depth, \(P\) is precipitation depth, \(I_a\) is initial abstraction,
!! and \(S\) is maximum potential retention.
!!
!! @author Original SWAP development team
!! @date Refactored February 2026
module runoff_mod
   use error_mod, only: fatalerr_collected
   use soilhydraulics_utils, only: watcon
   use swap_state_mod, only: swap_state_t  ! [SS-ATM A-2.6] nraidt/melt retired to state%atmosphere
   use atmosphere_constants_mod, only: DEPTH_10CM_CM, H_FIELD_CAPACITY_CM, &
                                       H_WILTING_POINT_CM, INITIAL_ABSTRACTION_RATIO
   implicit none
   private

   public :: cn_init, cn_step

contains

  !> Initialize the SCS Curve Number runoff method (formerly CNmethod(Itask=1, ...))
  !!
  !! Validates the CNtimTAB time series ordering, identifies the top soil layer
  !! (0-10 cm), and computes a reference water content for moisture corrections.
  subroutine cn_init(state)
    ! [SS-ATM A-2.6] nraidt/melt retired from variables; state added to read from state%atmosphere
    ! [SS-SWC S-2.12B] theta retired — read via state%soilwater%theta
    ! SS-TC TC-9: t1900 removed from only-list; read via state%timecontrol.
    ! [SS-GR-ATM B22] CN symbols → state%atmosphere%X (nod10_cn/icn_atm/z10_cn added to atmosphere_state)
    ! GR-BH: numnod, zbotcp, dz migrated to state%mesh%X
    implicit none
    ! global
    type(swap_state_t), intent(inout) :: state  ! [SS-ATM A-2.6] for retired nraidt/melt
    ! local
    integer              :: i
    real(8)              :: wc1, wc2

    associate (atmo => state%atmosphere, mesh => state%mesh, &
               soil => state%soilwater, time => state%timecontrol)

      atmo%icn_atm = 0
      ! check if times in CNtimeTAB are in ascending order
      ! set initial position in CNtimTAB
      do i = 2, atmo%iCNtab
          if (atmo%CNtimTAB(i) < atmo%CNtimTAB(i-1)) call fatalerr_collected('cn_init', 'CNtimTAB not in ascending order')
          if (time%t1900 >= atmo%CNtimTAB(i-1) .and. time%t1900 < atmo%CNtimTAB(i)) atmo%icn_atm = i-1
      end do
      ! error if start time t1900 not in CNtimTAB
      if (atmo%icn_atm == 0) call fatalerr_collected('cn_init', 'Start time of simulation not present in CNtimTAB')

    !  to be replaced by average for layer 0-10 cm
      do i = 1, mesh%numnod
          if (mesh%zbotcp(i) < -DEPTH_10CM_CM) then
            atmo%nod10_cn = i-1
            atmo%z10_cn = -mesh%zbotcp(atmo%nod10_cn)
            exit
          end if
      end do
      atmo%ThetaRef = 0.0d0
      do i = 1, atmo%nod10_cn
          if (atmo%wc_cor == 1) then
            wc1 = watcon(H_FIELD_CAPACITY_CM, soil%vg_params(i), &
                         soil%iHWCKmodel(soil%layer(i)), i, soil)
            wc2 = watcon(H_WILTING_POINT_CM, soil%vg_params(i), &
                         soil%iHWCKmodel(soil%layer(i)), i, soil)
            atmo%ThetaRef = atmo%ThetaRef + (wc1+wc2)*0.5d0*mesh%dz(i)
          else if (atmo%wc_cor == 2) then
            wc1 = watcon(0.0d0, soil%vg_params(i), &
                         soil%iHWCKmodel(soil%layer(i)), i, soil)
            wc2 = watcon(H_WILTING_POINT_CM, soil%vg_params(i), &
                         soil%iHWCKmodel(soil%layer(i)), i, soil)
            atmo%ThetaRef = atmo%ThetaRef + (wc1+wc2)*0.5d0*mesh%dz(i)
          end if
      end do
      atmo%ThetaRef = atmo%ThetaRef/atmo%z10_cn

    end associate

  end subroutine cn_init


  !> Compute runoff using the SCS Curve Number method (formerly CNmethod(Itask=2, ...))
  !!
  !! Updates the CN time-series index, computes moisture-adjusted CN values
  !! (per wc_cor mode), and writes Runoff_CN into state%atmosphere.
  subroutine cn_step(state)
    implicit none
    ! global
    type(swap_state_t), intent(inout) :: state
    ! local
    integer              :: i
    real(8)              :: CN, S, Ia

    associate (atmo => state%atmosphere, mesh => state%mesh, &
               soil => state%soilwater, time => state%timecontrol)

      ! Advance icn_atm if t1900 has moved past the next entry in CNtimTAB.
      do while (atmo%icn_atm < atmo%iCNtab .and. &
                time%t1900 >= atmo%CNtimTAB(atmo%icn_atm + 1))
        atmo%icn_atm = atmo%icn_atm + 1
      end do
      atmo%CNref = atmo%CNrefTAB(atmo%icn_atm)
      CN    = atmo%CNref
      atmo%CNdry =  4.2d0*atmo%CNref/(10.0d0-0.058d0*atmo%CNref)
      atmo%CNwet = 23.0d0*atmo%CNref/(10.0d0+0.13d0*atmo%CNref)

      if (atmo%wc_cor > 0) then
          atmo%wc10 = 0.0d0
          do i = 1, atmo%nod10_cn
            atmo%wc10 = atmo%wc10 + soil%theta(i)*mesh%dz(i)
          end do
          atmo%wc10 = atmo%wc10/atmo%z10_cn
          if (atmo%wc10 < atmo%ThetaRef) then
            CN = atmo%CNdry + atmo%wc10/atmo%ThetaRef*(atmo%CNref-atmo%CNdry)
          else
            CN = atmo%CNref + (atmo%wc10-atmo%ThetaRef)/atmo%ThetaRef*(atmo%CNwet-atmo%CNref)
          end if
      end if
      S  = 2540d0/CN - 25.4d0    ! in cm
      Ia = INITIAL_ABSTRACTION_RATIO*S
      if (atmo%nraidt+atmo%melt > Ia) then
          atmo%Runoff_CN = (atmo%nraidt+atmo%melt-Ia)**2 / &
                           (atmo%nraidt+atmo%melt-Ia+S)
      else
          atmo%Runoff_CN = 0.0d0
      end if

    end associate

  end subroutine cn_step

end module runoff_mod
