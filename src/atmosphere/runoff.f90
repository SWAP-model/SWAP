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

   public :: CNmethod

contains

  subroutine CNmethod(Itask, state)
    !> Calculate surface runoff using the SCS Curve Number method
    !!
    !! This subroutine operates in two modes controlled by the Itask parameter:
    !! - Itask=1: Initialization - validates time series, identifies top soil layer (0-10 cm),
    !!            and calculates reference water content for moisture corrections
    !! - Itask=2: Dynamic calculation - computes runoff for the current timestep using
    !!            moisture-adjusted CN values
    !!
    !! ## Moisture Correction Options (wc_cor)
    !! - 0: No moisture correction (CN = CNref)
    !! - 1: Field capacity based (\(h = -100\) cm and \(h = -16000\) cm)
    !! - 2: Saturation based (\(h = 0\) cm and \(h = -16000\) cm)
    !!
    !! ## References
    !! USDA-NRCS National Engineering Handbook, Part 630 Hydrology
    !!
    !! @warning CNref values > 170 will cause numerical issues in CNdry calculation
    ! [SS-ATM A-2.6] nraidt/melt retired from variables; state added to read from state%atmosphere
    ! [SS-SWC S-2.12B] theta retired — read via state%soilwater%theta
    ! SS-TC TC-9: t1900 removed from only-list; read via state%timecontrol.
    ! [SS-GR-ATM B22] CN symbols → state%atmosphere%X (nod10_cn/icn_atm/z10_cn added to atmosphere_state)
    ! GR-BH: numnod, zbotcp, dz migrated to state%mesh%X
    use soilhydraulics_utils, only: watcon
    implicit none
    ! global
    integer, intent(in)  :: Itask
    type(swap_state_t), intent(inout) :: state  ! [SS-ATM A-2.6] for retired nraidt/melt
    ! local
    integer              :: i
    real(8)              :: wc1, wc2, CN, S, Ia
    ! Note: nod10_cn/icn_atm/z10_cn migrated to state%atmosphere (B22)

    ! SS-TC TC-9: t1900 read via state%timecontrol (tc_* alias).
    associate( tc_t1900 => state%timecontrol%t1900 )  ! TC-9

    select case (Itask)
    ! initialization; calculate and store some constants
    case (1)

      state%atmosphere%icn_atm = 0
      ! check if times in CNtimeTAB are in ascending order
      ! set initial position in CNtimTAB
      do i = 2, state%atmosphere%iCNtab
          if (state%atmosphere%CNtimTAB(i) < state%atmosphere%CNtimTAB(i-1)) call fatalerr_collected ('CNmethod', 'CNtimTAB not in ascending order')
          if (tc_t1900 >= state%atmosphere%CNtimTAB(i-1) .and. tc_t1900 < state%atmosphere%CNtimTAB(i)) state%atmosphere%icn_atm = i-1
      end do
      ! error if start time t1900 not in CNtimTAB
      if (state%atmosphere%icn_atm == 0) call fatalerr_collected ('CNmethod', 'Start time of simulation not present in CNtimTAB')

    !  to be replaced by average for layer 0-10 cm
      do i = 1, state%mesh%numnod
          if (state%mesh%zbotcp(i) < -DEPTH_10CM_CM) then
            state%atmosphere%nod10_cn = i-1
            state%atmosphere%z10_cn = -state%mesh%zbotcp(state%atmosphere%nod10_cn)
            exit
          end if
      end do
      state%atmosphere%ThetaRef = 0.0d0
      do i = 1, state%atmosphere%nod10_cn
          if (state%atmosphere%wc_cor == 1) then
            wc1 = watcon(H_FIELD_CAPACITY_CM, &
                          state%soilwater%vg_params(i), &
                          state%soilwater%iHWCKmodel(state%soilwater%layer(i)), &
                          i, state%soilwater)                      ! [SS-GR-UTILS Task 5]
            wc2 = watcon(H_WILTING_POINT_CM, &
                          state%soilwater%vg_params(i), &
                          state%soilwater%iHWCKmodel(state%soilwater%layer(i)), &
                          i, state%soilwater)                      ! [SS-GR-UTILS Task 5]
            state%atmosphere%ThetaRef = state%atmosphere%ThetaRef + (wc1+wc2)*0.5d0*state%mesh%dz(i)
          else if (state%atmosphere%wc_cor == 2) then
            wc1 = watcon(0.0d0, &
                          state%soilwater%vg_params(i), &
                          state%soilwater%iHWCKmodel(state%soilwater%layer(i)), &
                          i, state%soilwater)                      ! [SS-GR-UTILS Task 5]
            wc2 = watcon(H_WILTING_POINT_CM, &
                          state%soilwater%vg_params(i), &
                          state%soilwater%iHWCKmodel(state%soilwater%layer(i)), &
                          i, state%soilwater)                      ! [SS-GR-UTILS Task 5]
            state%atmosphere%ThetaRef = state%atmosphere%ThetaRef + (wc1+wc2)*0.5d0*state%mesh%dz(i)
          end if
      end do
      state%atmosphere%ThetaRef = state%atmosphere%ThetaRef/state%atmosphere%z10_cn
      
      !!!t1900_old = int(t1900) - 1 ! for testing intermediate output
      
    ! dynamic part: calculate runoff   
    case (2)
      
      ! see if t1900 has moved ahead in CNtimTAB; icn_atm can never exceed last entry
      !  if (icn_atm < iCNtab .and. t1900 >= CNtimTAB(icn_atm+1)) icn_atm = icn_atm + 1
      ! Update position in CN time series if time has advanced (do while is more efficient if time steps are large and CN time series is long)
      do while (state%atmosphere%icn_atm < state%atmosphere%iCNtab .and. tc_t1900 >= state%atmosphere%CNtimTAB(state%atmosphere%icn_atm + 1))
        state%atmosphere%icn_atm = state%atmosphere%icn_atm + 1
      end do
      state%atmosphere%CNref = state%atmosphere%CNrefTAB(state%atmosphere%icn_atm)
      CN    = state%atmosphere%CNref
      state%atmosphere%CNdry =  4.2d0*state%atmosphere%CNref/(10.0d0-0.058d0*state%atmosphere%CNref)
      state%atmosphere%CNwet = 23.0d0*state%atmosphere%CNref/(10.0d0+0.13d0*state%atmosphere%CNref)

      if (state%atmosphere%wc_cor > 0) then
          state%atmosphere%wc10 = 0.0d0
          do i = 1, state%atmosphere%nod10_cn
            state%atmosphere%wc10 = state%atmosphere%wc10 + state%soilwater%theta(i)*state%mesh%dz(i)  ! [SS-SWC S-2.12B]
          end do
          state%atmosphere%wc10 = state%atmosphere%wc10/state%atmosphere%z10_cn
          if (state%atmosphere%wc10 < state%atmosphere%ThetaRef) then
            CN = state%atmosphere%CNdry + state%atmosphere%wc10/state%atmosphere%ThetaRef*(state%atmosphere%CNref-state%atmosphere%CNdry)
          else
            CN = state%atmosphere%CNref + (state%atmosphere%wc10-state%atmosphere%ThetaRef)/state%atmosphere%ThetaRef*(state%atmosphere%CNwet-state%atmosphere%CNref)
          end if
      end if
      S = 2540d0/CN-25.4d0    ! in cm
      Ia = INITIAL_ABSTRACTION_RATIO*S
    !   Ia = 0.3d0*S
      ! SS-ATM A-2.6: nraidt/melt retired — read from state%atmosphere
      if (state%atmosphere%nraidt+state%atmosphere%melt > Ia) then
          state%atmosphere%Runoff_CN = (state%atmosphere%nraidt+state%atmosphere%melt-Ia)**2/ &
                      (state%atmosphere%nraidt+state%atmosphere%melt-Ia+S)
      else
          state%atmosphere%Runoff_CN = 0.0d0
      end if
      
      ! for testing intermediate output
      !!!if (int(t1900) > t1900_old) then
      !!!   write (123, '(7F20.6)') t1900, nraidt, runoff_cn, cn, wc10, thetaref, melt
      !!!   t1900_old = int(t1900)
      !!!end if
      
    case default
      call fatalerr_collected ('CNmethod', 'Illegal Itask option')
    end select

    end associate  ! tc_t1900 => state%timecontrol [TC-9]

  end subroutine CNmethod

end module runoff_mod
