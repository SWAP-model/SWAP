!> Module for soil hydraulics calculations
!!
!! This module contains core hydraulic calculations for soil water flow,
!! including pressure head calculations, water content updates, and
!! hydraulic conductivity computations.
!!
module soilhydraulics_mod
   use error_mod, only: fatalerr_collected
   ! Core hydraulic calculations
   implicit none
   private
   public :: headcalc, soilwater, soilwaterstatevar, hysteresis
contains

   !> Calculate pressure heads, water contents, and conductivities for next time step
   !!
   !! Solves Richards' equation using a Newton-Raphson iteration scheme to update
   !! soil water state variables for the next time step. Handles various boundary
   !! conditions and includes optional macropore flow.
   !!
   !! @note
   !! Date: April 2005 / Sept 2005
   !! @endnote
   !!
   subroutine headcalc(state)
      use variables
      use timestep_control_mod, only: fldecdt
      use boundbottom_mod, only: BoundBottom
      use boundtop_mod, only: boundtop, PONDRUNOFF
      use rootextraction_mod, only: RootExtraction
      use array_utils, only: afgen
      use soilhydraulics_utils, only: watcon, hconduc, moiscap, hcomean, dhconduc
      use swap_constants, only: nihil
      ! [MACRO-RETIRE 2026-05-12] macropore_mod retired (ADR 0040).
      use swap_state_mod, only: swap_state_t
      use soilhydraulics_utils, only: dkmean
      use soilwaterbalance_mod, only: calcgwl, fluxes
      use numericalsolvers_mod, only: tridag, bandec, banbks
      implicit none

      ! Arguments
      type(swap_state_t), intent(inout) :: state

      ! Local variables
      integer   i,j, itry,  MaxIt1, NN, iBackTr
      real(8)   dFdhL(macp), dFdhM(macp), dFdhU(macp), difh(macp)
      real(8)   F(macp), factor, Fmax, QMpLatSsSav,qv(macp+1)
      real(8)   factmax, sink(macp), source(macp), sum, sum1, sumold, deviat
      real(8)   hold(macp)
      character(len=200) messag
      character(len=19) datetime
      character(len=11) datetmp
      logical   flnonconv,flnonconv1(macp), flnonconv2(macp), flnonconv3
      logical   flunsatok(3)       ! Flag indicating the performance of the iteration process
      logical   flboth
      integer   nodncr
      logical   flcaprise

      ! Convergence criteria
      real(8) CritDevBalCp, CritDevBalTot, Critdz

      data    CritDevBalCp   / 1.0d-6 / 
      data    CritDevBalTot  / 1.0d-5 / 
      data    Critdz         / 1.0d-5 / 
      ! ndr removed: loop now uses nrlevs (actual drain-level count, always <= Madr=5)

      real(8) hgrad(macp+1), dkdh(macp)

      ! Function and solver variables
      integer indx(macp), ierror
      real(8) a(macp,3), a1(macp,1), b(macp), d, q1
      logical flok

      ! Note: flwarn_hc, iwarn_hc, nstep_hc moved to variables.f90 module
      ! (previously local SAVE variables - now global for multi-instance support)

      ! [SS-SWC S-1.4a/S-1.4b/S-2.3] ASSOCIATE block: dual-writes + reader cutover in headcalc
      ! [TC-8] TC fields appended for dt/t1900/fldtmin/dtold/flDayStart reader cutover
      associate( &
         sw_theta       => state%soilwater%theta,    &
         sw_h           => state%soilwater%h,        &
         sw_k           => state%soilwater%k,        &
         sw_kmean       => state%soilwater%kmean,    &
         sw_dimoca      => state%soilwater%dimoca,   &
         sw_thetm1      => state%soilwater%thetm1,   &  ! [SS-SWC S-2.3] reader cutover
         sw_hm1         => state%soilwater%hm1,      &  ! [SS-SWC S-2.3] reader cutover
         sw_FrArMtrx    => state%soilwater%FrArMtrx, &  ! [SS-SWC S-2.3] reader cutover
         sw_cofgen      => state%soilwater%cofgen,   &  ! [SS-SWC S-2.3] reader cutover
         sw_nodgwl      => state%soilwater%nodgwl,   &  ! [SS-SWC S-2.3] reader cutover
         sw_evp         => state%soilwater%evp,      &  ! [SS-SWC S-2.3] reader cutover
         sw_gwlm1       => state%soilwater%gwlm1,    &  ! [SS-SWC S-2.3] reader cutover
         sw_pondm1      => state%soilwater%pondm1,   &  ! [SS-SWC S-2.3] reader cutover
         tc_dt          => state%timecontrol%dt,        &  ! [TC-8]
         tc_t1900       => state%timecontrol%t1900,     &  ! [TC-8]
         tc_fldtmin     => state%timecontrol%fldtmin,   &  ! [TC-8]
         tc_dtold       => state%timecontrol%dtold,     &  ! [TC-8]
         tc_flDayStart  => state%timecontrol%flDayStart, &  ! [TC-8]
         MaxIt          => state%timecontrol%MaxIt,      &  ! [SS-BMI2 Task 4]
         dtmin          => state%timecontrol%dtmin,      &  ! [SS-BMI2 Task 4]
         swscre         => state%timecontrol%swscre      )  ! [SS-BMI2 Task 4]

      if (tc_flDayStart) then  ! [TC-8]
         flwarn_hc = .true.
         iwarn_hc = 0
      endif
      call dtdpst                                                       &
     &        ('year-month-day,hour:minute:seconds',tc_t1900,datetime)  ! [TC-8]

      ! Summation of sink terms (constant for the current time step)
      iBackTr   = 0
      flunsatok(1) = .false.
      flunsatok(2) = .false.
      flunsatok(3) = .false.
      ! SS-DRST Task 3: read qdra from state%drainage (loop bound = nrlevs, not ndr=5)
      do i=1,numnod
         sink(i) = sw_evp(i)                                    ! [SS-SWC S-2.3] read from state
         if (allocated(state%drainage%qdra)) then
            do j=1,nrlevs
               sink(i) = sink(i) + state%drainage%qdra(j,i)
            end do
         end if
      end do
      source(1:numnod) = qssdi(1:numnod)

      ArMpSs = 0.d0                                            ! macropore retired (ADR 0040)

      ! Groundwater level specified — [SS-SWC S-2.12B] all legacy half-writes dropped
      if(swbotb.eq.1)then
         state%soilwater%fllowgwl = .false.               ! [SS-SWC S-1.6/S-2.12B]
         if(state%soilwater%gwlinp.ge.z(1)-1.0d-4)then

            q0 = (state%atmosphere%nraidt+nird+state%atmosphere%melt)*(1.0d0-ArMpSs) + state%soilwater%runon - state%soilwater%reva  ! [SS-ATM/SS-SWC S-2.12B]
            call pondrunoff (state)
            q1 = - q0 + (state%soilwater%pond - sw_pondm1)/tc_dt + state%soilwater%runots / tc_dt  ! [SS-SWC S-2.12B] [TC-8]
            sw_theta(1) = watcon(state%soilwater%gwlinp, &
                                  state%soilwater%vg_params(1), &
                                  state%soilwater%iHWCKmodel(state%soilwater%layer(1)), &
                                  1, state%soilwater)                         ! [SS-SWC S-1.4a/S-2.12B] [SS-GR-UTILS Task 5]
            sw_kmean(1) = hconduc(state%soilwater%gwlinp,sw_theta(1),state%heat%rfcp(1),state%heat%tsoil(1), &
                                  state%soilwater%vg_params(1), &
                                  state%soilwater%iHWCKmodel(state%soilwater%layer(1)), &
                                  state%soilwater%fluseksatexm(1), &
                                  1, state%soilwater)                          ! [SS-SWC S-1.4b/S-2.12B] [SS-GR-UTILS Task 6]

            qv(1) = q1
            do i=1,numnod
               qv(i+1) = qv(i) +dz(i)*sw_FrArMtrx(i)*(sw_theta(i)-sw_thetm1(i))  &  ! [SS-SWC S-2.12B]
     &                          / tc_dt+ sink(i) - source(i) + state%soilwater%qrot(i)  ! [TC-8]
            end do
            state%soilwater%qbot = qv(numnod+1)
            sw_h(1) = state%soilwater%gwlinp + disnod(1)*(qv(1)/sw_kmean(1)+1.0d0)  ! [SS-SWC S-1.4a/S-2.12B]
            do i=2,numnod
               sw_h(i) = sw_h(i-1) + disnod(i)*(qv(i)/sw_kmean(i)+1.0d0)    ! [SS-SWC S-1.4a/S-2.12B]
            end do

            if(SwKimpl.eq.1)then
               do i=1,numnod
                  sw_k(i) = hconduc(sw_h(i),sw_theta(i),state%heat%rfcp(i),state%heat%tsoil(i), &
                                    state%soilwater%vg_params(i), &
                                    state%soilwater%iHWCKmodel(state%soilwater%layer(i)), &
                                    state%soilwater%fluseksatexm(i), &
                                    i, state%soilwater)                        ! [SS-SWC S-1.4b/S-2.12B] [SS-GR-UTILS Task 6]
                  if(i.gt.1)then
                     sw_kmean(i) = hcomean(swkmean,sw_k(i-1),sw_k(i),dz(i-1),dz(i))  ! [SS-SWC S-1.4b/S-2.12B]
                  end if
               end do
               sw_kmean(numnod+1) = sw_k(numnod)                              ! [SS-SWC S-1.4b/S-2.12B]
            end if
            call calcgwl (state)
            return
         else
            NN = 0
            do while (z(NN+1).gt.state%soilwater%gwlinp .and. NN.lt.numnod)
               NN = NN + 1
            end do
            if (z(NN+1).lt.(state%soilwater%gwlinp+nihil)) then
               ! Groundwater within soil profile
               if ((z(NN)-state%soilwater%gwlinp) .lt. 1.0d-4 .and. (NN.gt.0)) then
                  ! Difference gwlinp with node too small to calculate gradient properly
                  state%soilwater%gwlinp = z(NN)
                  NN = NN-1
               endif
            else
               ! Groundwater below soil profile
               state%soilwater%fllowgwl = .true.          ! [SS-SWC S-1.6/S-2.12B]
               state%soilwater%hbot = state%soilwater%gwlinp - z(numnod) + 0.5*dz(numnod)
            endif
         end if
      else
         NN = numnod
      end if

      ! Reset conductivities to time level t

      ! Node nr of compartment with minimized flux of capillary rise
      if (swcaprise) then
         nodncr    = max(5,noddrz)
         flcaprise = .false.
      endif
      do i = 1,numnod
         sw_k(i) = hconduc(sw_h(i),sw_theta(i),state%heat%rfcp(i),state%heat%tsoil(i), &
                           state%soilwater%vg_params(i), &
                           state%soilwater%iHWCKmodel(state%soilwater%layer(i)), &
                           state%soilwater%fluseksatexm(i), &
                           i, state%soilwater)                                 ! [SS-GR-UTILS Task 6]

         if (swcaprise) then
            ! Prevent capillary rise into the root zone !! special for experts only
            if (i .eq. nodncr) then
               if ((sw_h(i) + z(i)) .lt. (sw_h(i+1) + z(i+1))) then  ! negative potential gradient upwards
                  sw_k(i)      = 1.0D-10
                  flcaprise = .true.
               endif
            endif
            if (flcaprise .and. i .eq. nodncr+1) then
                sw_k(i) = 1.0D-10
            endif
         endif
         if(i.gt.1)then
            sw_kmean(i) = hcomean(swkmean,sw_k(i-1),sw_k(i),dz(i-1),dz(i))
         end if
      enddo
      sw_kmean(numnod+1) = sw_k(numnod)

      if(SwKimpl.eq.0)then
         do i=2,numnod
            dFdhU(i)   = - sw_kmean(i)  /disnod(i)
            dFdhL(i-1) = dFdhU(i)
         end do
      end if

      do i=2,NN
         hgrad(i) = (sw_h(i-1)-sw_h(i))/disnod(i) + 1.0d0
      end do

      F(1) = (sw_theta(1)-sw_thetm1(1))*sw_FrArMtrx(1)*dz(1)/tc_dt + sink(1) - source(1) + state%soilwater%qrot(1) + sw_kmean(2) * hgrad(2)  ! [SS-SWC S-2.3] [TC-8]

      call boundtop(state)  ! [SS-HEAT] Task 9: state passed for rfcp access

      ! [MACRO-RETIRE 2026-05-12] MACROPORE(2,...) retired (ADR 0040).

      if (state%soilwater%FlRunoff) call pondrunoff (state)

      if(state%soilwater%ftoph)then
         hgrad(1) = (state%soilwater%hsurf-sw_h(1))/disnod(1) + 1.d0
         F(1)     = F(1) - sw_kmean(1) * hgrad(1)
      else
         F(1) = F(1) + state%soilwater%qtop
      end if

      do i=2,NN-1
         F(i) = (sw_theta(i)-sw_thetm1(i))*sw_FrArMtrx(i)*dz(i)/tc_dt + sink(i) - source(i) +     &  ! [SS-SWC S-2.3] [TC-8]
     &          state%soilwater%qrot(i) - sw_kmean(i) * hgrad(i) + sw_kmean(i+1) * hgrad(i+1)
      end do

      if(swbotb.eq.1 .and. (.not.state%soilwater%fllowgwl))then
         hgrad(NN+1) = sw_h(NN)/(z(nn)-state%soilwater%gwlinp) + 1.0d0
      else if(swbotb.eq.5 .or. (swbotb.eq.1 .and. state%soilwater%fllowgwl))then
         hgrad(NN+1) = (sw_h(NN) - state%soilwater%hbot) / disnod(NN+1)  + 1.0d0
      else if(swbotb.eq.8 .and. sw_h(NN).gt. Critdz - disnod(NN+1) + hplate) then
         hgrad(NN+1) = (sw_h(NN) - hplate) / disnod(NN+1)  + 1.0d0
         flboth = .true.
      else
         flboth = .false.
      end if

      if(swbotb.eq.1 .and. (.not.state%soilwater%fllowgwl))then
         sw_theta(NN)= watcon(sw_h(NN), &
                              state%soilwater%vg_params(NN), &
                              state%soilwater%iHWCKmodel(state%soilwater%layer(NN)), &
                              NN, state%soilwater)                            ! [SS-GR-UTILS Task 5]
         sw_k(NN)    = hconduc(sw_h(NN),sw_theta(NN),state%heat%rfcp(NN),state%heat%tsoil(NN), &
                              state%soilwater%vg_params(NN), &
                              state%soilwater%iHWCKmodel(state%soilwater%layer(NN)), &
                              state%soilwater%fluseksatexm(NN), &
                              NN, state%soilwater)                             ! [SS-GR-UTILS Task 6]
         sw_kmean(NN+1) = hcomean(swkmean,sw_k(NN),sw_cofgen(3,(NN+1)),       &  ! [SS-SWC S-2.3]
     &                        dz(NN),dz(NN+1))
         F(NN) = (sw_theta(NN) - sw_thetm1(NN))*sw_FrArMtrx(NN)*dz(NN)/tc_dt +      &  ! [SS-SWC S-2.3] [TC-8]
     &           sink(NN)-source(NN)+state%soilwater%qrot(NN)-sw_kmean(NN)*hgrad(NN) +sw_kmean(NN+1)*hgrad(NN+1)
      else

         F(NN) = (sw_theta(NN) - sw_thetm1(NN))*sw_FrArMtrx(NN)*dz(NN)/tc_dt        &  ! [SS-SWC S-2.3] [TC-8]
     &         - sw_kmean(NN) * hgrad(NN) + sink(NN) - source(NN) + state%soilwater%qrot(NN)

         if(swbotb.eq.3.and.swbotb3Impl.eq.1)then ! Cauchy-relation, implemented as head boundary
            if (SwBotb3ResVert.eq.0) then
               state%soilwater%qbot = - (sw_h(NN)+z(NN)-state%soilwater%deepgw) / (disnod(NN+1)/sw_kmean(NN+1)+rimlay)
            elseif (SwBotb3ResVert.eq.1) then
               state%soilwater%qbot = - (sw_h(NN)+z(NN)-state%soilwater%deepgw) / rimlay
            endif
! ---       extra groundwater flux might be added
            if (sw4 .eq. 1) state%soilwater%qbot = state%soilwater%qbot + afgen(qbotab,mabbc*2,tc_t1900+tc_dt)  ! [TC-8]
            F(NN) = F(NN) - state%soilwater%qbot
         else if(swbotb.eq.5 .or. (swbotb.eq.1 .and. state%soilwater%fllowgwl))then ! pressure head at lower boundary specified
            F(NN) = F(NN) + sw_kmean(NN+1) * hgrad(NN+1)
         else if(swbotb.eq.7 .or. swbotb .eq. -2)then ! free drainage option
            sw_kmean(numnod+1) = hconduc(sw_h(numnod),sw_theta(numnod),state%heat%rfcp(numnod),state%heat%tsoil(numnod), &
                                         state%soilwater%vg_params(numnod), &
                                         state%soilwater%iHWCKmodel(state%soilwater%layer(numnod)), &
                                         state%soilwater%fluseksatexm(numnod), &
                                         numnod, state%soilwater)              ! [SS-GR-UTILS Task 6]
            state%soilwater%qbot = -1.0d0 * sw_kmean(numnod+1)
            F(NN) = F(NN) - state%soilwater%qbot
         ! Lysimeter option
         else if(swbotb.eq.8)then
            if (flboth) then
               state%soilwater%hbot = hplate
               F(NN) = F(NN) + sw_kmean(NN+1) * hgrad(NN+1)
            else
               state%soilwater%qbot = 0.0d0
            end if
         ! Flux bottom boundary
         else
            F(NN) = F(NN) - state%soilwater%qbot
         end if

      end if

      ! Initial estimate of F inner product
      sumold    = 0.0d0
      do i=1,NN
         sumold = sumold + F(i)*F(i)
      end do
      sumold = 0.5d0 * sumold

      ! Start iteration loop, MaxIt specified in the input
      if(tc_fldtmin)then  ! [TC-8]
         MaxIt1 = 2*MaxIt
      else
         MaxIt1 = MaxIt
      end if

      sum= 0.d0   ! For Forcheck
      Do numbit = 1,MaxIt1

         do i = 1, NN
            ! Save values of h
            hold(i)   = sw_h(i)

            ! Derivative of theta to h (differential moisture capacity),
            ! as part of main diagonal
            sw_dimoca(i) = moiscap(sw_h(i), &
                                    state%soilwater%vg_params(i), &
                                    state%soilwater%iHWCKmodel(state%soilwater%layer(i)), &
                                    state%timecontrol%dt, &
                                    i, state%soilwater)              ! [SS-GR-UTILS Task 7]
            sw_dimoca(i) = sw_dimoca(i)                              ! [SS-SWC S-1.4b]

         enddo

         if(SwKimpl.eq.1)then
            do i = 1, NN
               dkdh(i) = dhconduc(sw_h(i),sw_theta(i),sw_dimoca(i),state%heat%rfcp(i), &
                                   state%soilwater%vg_params(i), &
                                   state%soilwater%iHWCKmodel(state%soilwater%layer(i)), &
                                   i, state%soilwater)                         ! [SS-GR-UTILS Task 6]
            enddo
            do i=2,NN
               dFdhU(i)   = - sw_kmean(i)  /disnod(i)
               dFdhL(i-1) = dFdhU(i)
            end do
         end if

         if (swcaprise .and. flcaprise) then
           dkdh(nodncr+1) = 1.0D-30
         endif

         ! Jacobian matrix elements
         dFdhM(1) = sw_dimoca(1)*sw_FrArMtrx(1)*dz(1)/tc_dt - dFdhL(1)  ! [SS-SWC S-2.3] [TC-8]
         ! If the head boundary condition applies: add the k1/(0.5*dz1) term
         ! to the first element of the main diagonal
         if(state%soilwater%ftoph) dFdhM(1) = dFdhM(1) + sw_kmean(1)/disnod(1)

         do i=2,NN-1
            dFdhM(i) = sw_dimoca(i)*sw_FrArMtrx(i)*dz(i)/tc_dt - dFdhU(i)        &  ! [SS-SWC S-2.3] [TC-8]
     &                                                - dFdhL(i)
         end do

         dFdhM(NN) = sw_dimoca(NN)*sw_FrArMtrx(NN)*dz(NN)/tc_dt - dFdhU(NN)  ! [SS-SWC S-2.3] [TC-8]
         if(swbotb.eq.1 .and. (.not.state%soilwater%fllowgwl))then
            dFdhM(NN) = dFdhM(NN) + sw_kmean(NN+1)/(z(NN)-state%soilwater%gwlinp)
         else if(swbotb.eq.3.and.swbotb3Impl.eq.1)then ! Cauchy
            if (SwBotb3ResVert.eq.0) then
               dFdhM(NN) = dFdhM(NN) + 1.0d0 /                          &
     &                              (disnod(NN+1)/sw_kmean(NN+1)+rimlay)   
            elseif (SwBotb3ResVert.eq.1) then
               dFdhM(NN) = dFdhM(NN) + 1.0d0 / rimlay
            endif
         else if(swbotb.eq.5 .or. (swbotb.eq.1 .and. state%soilwater%fllowgwl))then
            dFdhM(NN) = dFdhM(NN) + sw_kmean(NN+1)/disnod(NN+1)         
         else if(swbotb.eq.7 .or. swbotb .eq. -2)then ! implicitly: sw_kmean(NN+1)
            dFdhM(NN) = dFdhM(NN) + dkdh(NN) * 0.5d0
         else if(swbotb.eq.8 .and. flboth)then
            dFdhM(NN) = dFdhM(NN) + sw_kmean(NN+1)/disnod(NN+1)
         end if

         if(SwKimpl.eq.1)then
            dFdhM(1) = dFdhM(1) + dkdh(1) * hgrad(2) *                  &
     &                 dkmean(swkmean,sw_k(1),sw_k(2),dz(1),dz(2))
            if(state%soilwater%ftoph) dFdhM(1) = dFdhM(1) - dkdh(1) * hgrad(1) * 0.5d0
            dFdhL(1) = dFdhL(1) + dkdh(2) * hgrad(2) *                  &
     &                 dkmean(swkmean,sw_k(2),sw_k(1),dz(2),dz(1)) 
            do i=2,NN-1
               dFdhU(i) = dFdhU(i) - dkdh(i-1) * hgrad(i) *             &
     &                    dkmean(swkmean,sw_k(i-1),sw_k(i),dz(i-1),dz(i)) 
               dFdhM(i) = dFdhM(i) - dkdh(i) * hgrad(i) *               &
     &                    dkmean(swkmean,sw_k(i),sw_k(i-1),dz(i),dz(i-1))     &
     &                             + dkdh(i) * hgrad(i+1) *             &
     &                    dkmean(swkmean,sw_k(i),sw_k(i+1),dz(i),dz(i+1))
               dFdhL(i) = dFdhL(i) + dkdh(i+1) * hgrad(i+1) *           &
     &                    dkmean(swkmean,sw_k(i+1),sw_k(i),dz(i+1),dz(i)) 
            end do
            dFdhU(NN) = dFdhU(NN) - dkdh(NN-1) * hgrad(NN) *            &
     &                  dkmean(swkmean,sw_k(NN-1),sw_k(NN),dz(NN-1),dz(NN)) 
            dFdhM(NN) = dFdhM(NN) - dkdh(NN) * hgrad(NN) *              &
     &                  dkmean(swkmean,sw_k(NN),sw_k(NN-1),dz(NN),dz(NN-1))

            if(swbotb.eq.1 .or. swbotb.eq.5 .or. swbotb.eq.8            &
     &         .and. flboth)then
               dFdhM(NN) = dFdhM(NN) + 0.5d0 * dkdh(NN) * hgrad(NN+1)
            end if
         end if

         ! [MACRO-RETIRE 2026-05-12] MACROPORE(3,...) retired (ADR 0040).

         ! Solve the tridiagonal matrix
         call tridag(NN, dFdhU, dFdhM, dFdhL, F, difh, ierror)

         if(ierror.ne.0)then
            call dtdpst ('year-month-day',tc_t1900+1.001d0,datetmp)  ! [TC-8]
            messag = ' Tri-band matrix in HeadCalc appeared to be'//    &
     &               ' singular at '//datetmp//                            &
     &               '   Alternative SOLVER chosen'
            call warn ('Headcalc',messag,logf,swscre)
            do i=1,NN
               a(i,1) = dFdhU(i)
               a(i,2) = dFdhM(i)
               a(i,3) = dFdhL(i)
            end do
            call bandec(a,nn,1,1,macp,3,a1,1,indx,d)
            do i=1,NN
               b(i) = F(i)
            end do
            call banbks(a,nn,1,1,macp,3,a1,1,indx,b)
            do i=1,NN
               difh(i) = b(i)
            end do
         end if

         factor  = 1.0d0
         do itry = 1,MaxBackTr
            iBackTr = iBackTr + 1
            ! Factor reduces the change of h (difh) calculated as a full
            ! Newton Raphson step
            if(tc_fldtmin .and. numbit.gt.MaxIt)then              ! [TC-8]
               factmax = 0.0d0
               do i = 1,NN
                  if(dabs( hold(i) ) .lt. 1.0d0 )then
                     factmax = max( factmax, dabs( difh(i) ) ) 
                  else
                     factmax = max( factmax, dabs( difh(i) / hold(i) ) )
                  end if
               end do
               do i = 1,NN
                  sw_h(i) = hold(i) - difh(i) * min(1.0d0, 1.0d0 / factmax)
                  sw_h(i) = sw_h(i)                                ! [SS-SWC S-1.4a]
               end do
            else
               do i = 1,NN
                  sw_h(i) = hold(i) - factor * difh(i)
                  sw_h(i) = sw_h(i)                                ! [SS-SWC S-1.4a]
               end do
            end if
            do i = 1,NN
              sw_theta(i) = watcon(sw_h(i), &
                                   state%soilwater%vg_params(i), &
                                   state%soilwater%iHWCKmodel(state%soilwater%layer(i)), &
                                   i, state%soilwater)            ! [SS-GR-UTILS Task 5]
              sw_theta(i) = sw_theta(i)                            ! [SS-SWC S-1.4a]
            enddo
            do i=2,NN
               hgrad(i) = (sw_h(i-1)-sw_h(i))/disnod(i) + 1.0d0
            end do


            if(SwKimpl.eq.1)then
               call Rootextraction(state)
               do i = 1,NN
                  sw_k(i) = hconduc(sw_h(i),sw_theta(i),state%heat%rfcp(i),state%heat%tsoil(i), &
                                    state%soilwater%vg_params(i), &
                                    state%soilwater%iHWCKmodel(state%soilwater%layer(i)), &
                                    state%soilwater%fluseksatexm(i), &
                                    i, state%soilwater)                        ! [SS-GR-UTILS Task 6]
                  sw_k(i) = sw_k(i)                                  ! [SS-SWC S-1.4b]
                  if(i.gt.1)then
                     sw_kmean(i)=hcomean(swkmean,sw_k(i-1),sw_k(i),dz(i-1),dz(i))
                     sw_kmean(i) = sw_kmean(i)                        ! [SS-SWC S-1.4b]
                  end if
               end do
               sw_kmean(NN+1) = sw_k(NN)
               sw_kmean(NN+1) = sw_k(NN)                              ! [SS-SWC S-1.4b]
            end if

            ! Prevent capillary rise into the root zone !! special for experts
            if (swcaprise) then
               flcaprise = .false.
               i = nodncr
               if ((sw_h(i) + z(i)) .lt. (sw_h(i+1) + z(i+1))) then  ! negative potential gradient upwards
                  sw_k(i)      = 1.0D-10
                  flcaprise = .true.
               endif
               if (flcaprise) then
                 sw_k(i+1) = 1.0D-10
               endif
               if(i.gt.1)then
                 sw_kmean(i)=hcomean(swkmean,sw_k(i-1),sw_k(i),dz(i-1),dz(i))
                 sw_kmean(i) = sw_kmean(i)                            ! [SS-SWC S-1.4b]
                 sw_kmean(i+1)=hcomean(swkmean,sw_k(i),sw_k(i+1),dz(i),dz(i+1))
                 sw_kmean(i+1) = sw_kmean(i+1)                        ! [SS-SWC S-1.4b]
               end if
               sw_k(i) = sw_k(i)                                      ! [SS-SWC S-1.4b] sw_k(nodncr)
               if (flcaprise) sw_k(i+1) = sw_k(i+1)                   ! [SS-SWC S-1.4b] sw_k(nodncr+1)
            endif

            ! Calculate F-function
            F(1) = (sw_theta(1) - sw_thetm1(1))*sw_FrArMtrx(1)*dz(1)/tc_dt + sink(1) - source(1) &  ! [SS-SWC S-2.3] [TC-8]
     &           + state%soilwater%qrot(1) + sw_kmean(2) * hgrad(2)

            ! [MACRO-RETIRE 2026-05-12] FlMacropore QMpLatSsSav save retired (ADR 0040).

            call boundtop(state)  ! [SS-HEAT] Task 9: state passed for rfcp access

            ! [MACRO-RETIRE 2026-05-12] MACROPORE(2,...) retired (ADR 0040).

            if (state%soilwater%FlRunoff) call pondrunoff (state)

            if(state%soilwater%ftoph)then
               hgrad(1) = (state%soilwater%hsurf-sw_h(1))/disnod(1) + 1.d0
               F(1) = F(1) - sw_kmean(1) * hgrad(1)
            else
               F(1) = F(1) + state%soilwater%qtop
            end if

            do i=2,NN-1
               F(i) = (sw_theta(i)-sw_thetm1(i))*sw_FrArMtrx(i)*dz(i)/tc_dt + sink(i) - source(i) &  ! [SS-SWC S-2.3] [TC-8]
     &              + state%soilwater%qrot(i) - sw_kmean(i)*hgrad(i)+sw_kmean(i+1)*hgrad(i+1)
            end do

            if(swbotb.eq.1 .and. (.not.state%soilwater%fllowgwl))then
               hgrad(NN+1) = sw_h(NN)/(z(nn)-state%soilwater%gwlinp) + 1.0d0
           else if(swbotb.eq.5 .or. (swbotb.eq.1 .and. state%soilwater%fllowgwl))then
               hgrad(NN+1) = (sw_h(NN) - state%soilwater%hbot) / disnod(NN+1)  + 1.0d0
            else if(swbotb.eq.8 .and. flboth)then
               hgrad(NN+1) = (sw_h(NN) - hplate) / disnod(NN+1)  + 1.0d0
            end if

            if(swbotb.eq.1 .and. (.not.state%soilwater%fllowgwl))then
               sw_theta(NN) = watcon(sw_h(NN), &
                                     state%soilwater%vg_params(NN), &
                                     state%soilwater%iHWCKmodel(state%soilwater%layer(NN)), &
                                     NN, state%soilwater)          ! [SS-GR-UTILS Task 5]
               sw_theta(NN) = sw_theta(NN)                         ! [SS-SWC S-1.4a]
               sw_k(NN)     = hconduc(sw_h(NN),sw_theta(NN),state%heat%rfcp(NN),state%heat%tsoil(NN), &
                                      state%soilwater%vg_params(NN), &
                                      state%soilwater%iHWCKmodel(state%soilwater%layer(NN)), &
                                      state%soilwater%fluseksatexm(NN), &
                                      NN, state%soilwater)                     ! [SS-GR-UTILS Task 6]
               sw_k(NN) = sw_k(NN)                                 ! [SS-SWC S-1.4b]
               sw_kmean(NN+1) = hcomean(swkmean,sw_k(NN),sw_cofgen(3,(NN+1)) &  ! [SS-SWC S-2.3]
     &                       ,dz(NN),dz(NN+1))
               sw_kmean(NN+1) = sw_kmean(NN+1)                     ! [SS-SWC S-1.4b]
               F(NN) = (sw_theta(NN) - sw_thetm1(NN))*sw_FrArMtrx(NN)*dz(NN)/tc_dt  &  ! [SS-SWC S-2.3] [TC-8]
     &            - sw_kmean(NN) * hgrad(NN) + sw_kmean(NN+1) * hgrad(NN+1)   &
     &            + sink(NN) - source(NN) + state%soilwater%qrot(NN)
            else
               F(NN) = (sw_theta(NN) - sw_thetm1(NN))*sw_FrArMtrx(NN)*dz(NN)/tc_dt  &  ! [SS-SWC S-2.3] [TC-8]
     &               - sw_kmean(NN) * hgrad(NN)                            &
     &               + sink(NN) - source(NN) + state%soilwater%qrot(NN)
               if(swbotb.eq.3.and.swbotb3Impl.eq.1)then ! Cauchy
                  if (SwBotb3ResVert.eq.0) then
                     state%soilwater%qbot = - (sw_h(NN)+z(NN)-state%soilwater%deepgw) /       &
     &                                 (disnod(NN+1)/sw_kmean(NN+1)+rimlay)
                  elseif (SwBotb3ResVert.eq.1) then
                     state%soilwater%qbot = - (sw_h(NN)+z(NN)-state%soilwater%deepgw) / rimlay
                  endif
                  ! Extra groundwater flux might be added
                  if (sw4 .eq. 1) then
                     state%soilwater%qbot = state%soilwater%qbot + afgen(qbotab,mabbc*2,tc_t1900+tc_dt)  ! [TC-8]
                  end if
                  F(NN) = F(NN) - state%soilwater%qbot
               else if(swbotb.eq.5 .or.(swbotb.eq.1 .and. state%soilwater%fllowgwl))then
                  ! Pressure head at lower boundary specified
                  F(NN) = F(NN) + sw_kmean(NN+1) * hgrad(NN+1)
               else if(swbotb.eq.7.or. swbotb .eq. -2)then ! free drainage option
                  sw_kmean(numnod+1) = hconduc(sw_h(numnod),sw_theta(numnod),state%heat%rfcp(numnod),state%heat%tsoil(numnod), &
                                               state%soilwater%vg_params(numnod), &
                                               state%soilwater%iHWCKmodel(state%soilwater%layer(numnod)), &
                                               state%soilwater%fluseksatexm(numnod), &
                                               numnod, state%soilwater)        ! [SS-GR-UTILS Task 6]
                  sw_kmean(numnod+1) = sw_kmean(numnod+1)          ! [SS-SWC S-1.4b]
                  state%soilwater%qbot = -1.0d0 * sw_kmean(numnod+1)
                  F(NN) = F(NN) - state%soilwater%qbot
               ! Lysimeter option
               else if(swbotb.eq.8)then
                  if (flboth) then
                     state%soilwater%hbot = hplate
                     F(NN) = F(NN) + sw_kmean(NN+1) * hgrad(NN+1)
                  else
                     state%soilwater%qbot = 0.0d0
                  end if
               ! Flux bottom boundary
               else
                  F(NN) = F(NN) - state%soilwater%qbot
               end if

            end if

            ! Calculate maximum deviation per compartment and new inner product
            Fmax   = 0.0d0
            sum    = 0.0d0
            sum1   = 0.0d0
!            numnodFmax = 0                            ! disabled for Forcheck 
            do i=1,NN
               if(dabs(F(i)).gt.Fmax)then
                  Fmax = dabs(F(i))
!                  numnodFmax = i                      ! disabled for Forcheck 
               end if
               sum = sum + F(i)*F(i)
               sum1=sum1 + F(i)  ! in cm/d
            end do
            sum = 0.5d0 * sum

            ! Test for iteration progress, if Newton-step is too large:
            ! reduce dh by multiplication factor
            if(sum.lt.sumold .or. Fmax.lt.CritDevBalCp) goto 1000

            factor = factor / 3.0d0


         end do
 1000    continue

         ! Check on convergence of solution

         ! Initialize flags
         ! Main flag for testing the convergence
         flnonconv = .false.

         ! Flags introduced for debugging purposes
         do i = 1,numnod
            flnonconv1(i) = .false.
            flnonconv2(i) = .false.
         end do
         flnonconv3 = .false.
         flnonconv1(1) = flnonconv1(1) ! for Forcheck

         ! Apply performance criteria per compartment
         do i = 1,NN
            ! Test for water balance deviation of soil compartments
            if( dabs(F(i)).gt. CritDevBalCp)then
                  flnonconv1(i) = .true. ; flnonconv   = .true.
            end if
            ! Test for change of pressure head
            if( dabs(hold(i)) .lt. 1.0d0)then
               if(abs( sw_h(i)-hold(i) ) .gt. CritDevh2Cp)then
                  flnonconv2(i) = .true. ; flnonconv   = .true.
               endif
            else
               if(abs( sw_h(i)-hold(i) )/abs(hold(i)) .gt.                 &
     &                CritDevh1Cp)then
                  flnonconv2(i) = .true. ; flnonconv   = .true.
               endif
            end if
            flnonconv2(1) = flnonconv2(1) ! for Forcheck

         enddo

         ! Test for waterbalance of ponding layer
         if (state%soilwater%ftoph) then
            state%soilwater%qtop = -sw_kmean(1)*((state%soilwater%hsurf - sw_h(1))/disnod(1)+1.0d0)
            if(.not.flnonconv) then
               deviat = state%soilwater%pond - sw_pondm1 + state%soilwater%reva*tc_dt - (state%atmosphere%nraidt+nird+state%atmosphere%melt)*tc_dt &  ! [SS-SWC S-2.12B] [TC-8]
     &                - state%soilwater%runon*tc_dt  +  state%soilwater%runots  - state%soilwater%qtop * tc_dt  ! [TC-8]
               if( abs(deviat) .gt. CritDevPondDt) then
                  flnonconv3 = .true. ; flnonconv   = .true.
                  flnonconv3 = flnonconv3 ! for Forcheck
               end if
            end if
         end if

         if(dabs(sum1).gt.CritDevBalTot) flnonconv   = .true.

         ! Save sum for next iteration
         sumold = sum

         if(.not.flnonconv )then !  convergence has been reached
            if(swbotb.eq.1 .and. (.not.state%soilwater%fllowgwl))then
               ! Derive vertical flux profile in order to find qbot as a
               ! lower boundary condition for the saturated part of the soil
               ! system
               qv(1) = state%soilwater%qtop
               do i=NN+1,numnod
                  sw_theta(i) = sw_cofgen(2,i)                     ! [SS-SWC S-2.3]
                  sw_theta(i) = sw_theta(i)                        ! [SS-SWC S-1.4a]
               end do
               do i=1,numnod
                 qv(i+1) = qv(i) +dz(i)*sw_FrArMtrx(i)*(sw_theta(i)-sw_thetm1(i))&  ! [SS-SWC S-2.3]
     &                          / tc_dt + sink(i) - source(i) + state%soilwater%qrot(i)  ! [TC-8]
               end do
               state%soilwater%qbot = qv(numnod+1)

               do i=NN+1,numnod
                  sw_h(i) = sw_h(i-1) + disnod(i)*(qv(i)/sw_kmean(i)+1.0d0)
                  sw_h(i) = sw_h(i)                                ! [SS-SWC S-1.4a]
               end do

            end if

            ! Calculate new groundwater level
            call calcgwl (state)

            if(swbotb.ne.1.and.abs(state%soilwater%gwl-sw_gwlm1).ge.gwlconv .AND.          &  ! [SS-SWC S-2.12B]
     &         abs(state%soilwater%gwl-999d0).gt.1.d0.and.abs(sw_gwlm1-999d0).gt.1.d0) then  ! [SS-SWC S-2.12B]
               call dtdpst ('year-month-day',tc_t1900+1.001d0,datetmp)  ! [TC-8]
                 messag = ' Change of groundwater level exceeds'//      &
     &           ' criterion at '//datetmp//'. Consider reduction of dtMin'
               call warn ('Headcalc',messag,logf,swscre)
            endif

            ! Recording of number of iteration steps needed
            itnumb(min(100,numbit),1)=itnumb(min(100,numbit),1)+1 
            itnumb(min(100,numbit),2)=itnumb(min(100,numbit),2)+iBackTr 

            return
            
         end if
      End Do

      ! Convergence could not been reached
      if (.not.tc_fldtmin ) then  ! [TC-8]
         ! Reset soil state variables
         do j = 1,numnod
            sw_h(j) = sw_hm1(j)                                    ! [SS-SWC S-2.3]
            sw_h(j) = sw_h(j)                                      ! [SS-SWC S-1.4a]
            sw_theta(j) = sw_thetm1(j)                             ! [SS-SWC S-2.3]
            sw_theta(j) = sw_theta(j)                              ! [SS-SWC S-1.4a]
         enddo
         sw_kmean(numnod+1) = sw_k(numnod)
         sw_kmean(numnod+1) = sw_k(numnod)                          ! [SS-SWC S-1.4b]
         state%soilwater%gwl  = sw_gwlm1                         ! [SS-SWC S-2.12B]
         state%soilwater%pond = sw_pondm1                        ! [SS-SWC S-2.12B]

         ! Reset and continue iteration with smaller timestep!
         fldecdt = .true.

         return

      else
         ! Write warning to screen and log file
         if (flwarn_hc .and. iwarn_hc.lt.5) then
            iwarn_hc = iwarn_hc + 1
            call dtdpst                                                 &
     &        ('year-month-day,hour:minute:seconds',tc_t1900,datetime)  ! [TC-8]
            messag = ' No convergence was reached of Richards'//        &
     &        ' equation at '//datetime//                               &
     &        ' no more than 4 warnings per date - SWAP did continue !'
            call warn ('Headcalc',messag,logf,swscre)
            if (iwarn_hc.gt.4) then
              flwarn_hc = .false.  
            endif
         endif

         if(fldumpconvcrit) then
           write(logf,'(a,a19)')   'Datetime = ',datetime
           write(logf,'(a,f14.6)') 't1900    = ', tc_t1900  ! [TC-8]
           write(logf,'(a,f10.6)') 'dtmin = ', dtmin  ! [SS-BMI2 Task 4] via state%timecontrol alias
           write(logf,'(a,f10.6)') 'dt    = ', tc_dt  ! [TC-8]
           write(logf,'(a,i3)')    'ftoph  = ', state%soilwater%ftoph
           write(logf,'(a,f10.6)') 'CritDevBalCp  = ', CritDevBalCp
           write(logf,'(a,f10.6)') 'CritDevBalTot = ', CritDevBalTot
           write(logf,'(a,f10.6)') 'CritDz        = ', CritDz
           write(logf,'(a,f10.6)') 'CritDevh1Cp   = ', CritDevh1Cp
           write(logf,'(a,f10.6)') 'CritDevh2Cp   = ', CritDevh2Cp
           write(logf,'(a,i3)') 'flnonconv  = ', flnonconv
           write(logf,'(a,i3)') 'flnonconv3 = ', flnonconv3
           write(logf,'(a,i3)') 'flunsatok(1) = ', flunsatok(1)
           write(logf,'(a,i3)') 'flunsatok(2) = ', flunsatok(2)
           write(logf,'(a,i3)') 'flunsatok(3) = ', flunsatok(3)

           write(logf,'(a,f10.6)') 'pondm1 = ', sw_pondm1             ! [SS-SWC S-2.3]
           write(logf,'(a,f10.6)') 'pond   = ', state%soilwater%pond    ! [SS-SWC S-2.12B]
           write(logf,'(a,f10.6)') 'gwlm1 = ', sw_gwlm1              ! [SS-SWC S-2.3]
           write(logf,'(a,f10.6)') 'gwl   = ', state%soilwater%gwl     ! [SS-SWC S-2.12B]
           write(logf,'(a)')                                            &
     &      'node,flnonconv1_F, F, flnonconv2_h, hm1,     h,'//         &
     &      '    thetm1,    theta'
           do j = 1,numnod
              write(logf,'(2(i4,a),f10.6,a,i3,5(a,f10.6))')             &
     &            j,',',flnonconv1(j),',',F(j),',',flnonconv2(j),',',   &
     &            sw_h(j),',',sw_hm1(j),',', sw_theta(j),',',sw_thetm1(j)  ! [SS-SWC S-2.3]
           enddo

         endif

         ! Continue without convergence !!!
         return

      endif

      end associate  ! sw_theta/.../sw_gwlm1 => state%soilwater [SS-SWC S-1.4a/b/S-2.3]

   end subroutine headcalc

   !> Calculate soil water state variables
   !!
   !! Initializes or updates soil water state variables including water content,
   !! pressure head, hydraulic conductivity, and related parameters. Handles
   !! different initialization options and boundary conditions.
   !!
   !! @param[in] task Task selector: 1=initialize, 2=calculate rates, 3=update states
   !!
   !! @note
   !! Date: Aug 2004
   !! @endnote
   !!
   subroutine soilwater(task, state)
     use doln
      use Variables
      use swap_log, only: log_info, to_str
      use array_utils, only: afgen
      use soilhydraulics_utils, only: watcon, hconduc, moiscap, hcomean
      ! [MACRO-RETIRE 2026-05-12] macropore_mod retired (ADR 0040).
      use soilwaterbalance_mod, only: calcgwl, watstor, integral, fluxes
      use swap_state_mod, only: swap_state_t
      use, intrinsic :: iso_fortran_env, only: real64  ! [SS-ATM] for atmosphere state dual-writes
      implicit none

      ! Arguments
      integer task
      type(swap_state_t), intent(inout) :: state

      ! Local variables
      integer lay,node,i,j

      real(8) tab(mabbc*2)
      character(len=200) messag

      select case (task)
      case (1)

         ! Initialize Soilwater rate/state variables

         ! [SS-SWC S-1.3] ASSOCIATE block: dual-writes mirror every legacy global
         ! write into state%soilwater for all ~30 init-path fields.
         associate(sw => state%soilwater)

         ! Initialize miscellaneous
         ! [SS-SWC S-2.12B] legacy half-writes dropped — state%soilwater is canonical
         sw%hatm = -2.75e5_real64                           ! [SS-SWC S-1.3/S-2.12B]
      state%atmosphere%nraidt = 0.0_real64
      nird = 0.0d0
      if (swinco.ne.3) then
        state%atmosphere%ldwet = 0.0_real64
        state%atmosphere%spev  = 0.0_real64
        state%atmosphere%saev  = 0.0_real64
      endif
      sw%runon = 0.0_real64                                 ! [SS-SWC S-1.3/S-2.12B]
      state%soilwater%qtop = 0.d0
      do i = 1,numnod+1
        sw%q(i) = 0.0_real64                               ! [SS-SWC S-1.3/S-2.12B]
      enddo
      sw%evp = 0.0_real64                                  ! [SS-SWC S-1.3/S-2.12B] state sized numnod; blanket zero
      ! [SS-HEAT] Task 9: legacy rfcp global retired; state%heat%rfcp is authoritative
      if (allocated(state%heat%rfcp)) state%heat%rfcp = 1.0d0
      state%surfacewater%vtair = 0.0d0
      cQMpLatSs = 0.0d0

      ! Soil physics: tabulated or MualemVanGenuchten functions
      sw%cofgen = 0.0_real64                               ! [SS-SWC S-1.3/S-2.12B]

      ! [SS-GR-UTILS Task 4] Mirror layer + iHWCKmodel into state (both branches).
      ! layer(:) set by CalcGrid; iHWCKmodel(:) set by config_to_variables.
      ! State fields allocated by soilwater_init which ran before SoilHydraulics(1).
      do node = 1, numnod
         sw%layer(node) = layer(node)
      end do
      do lay = 1, numlay
         sw%iHWCKmodel(lay) = iHWCKmodel(lay)
      end do
      ! BiModal/NoVap: only set via legacy readswap (not TOML path); stay .false.

      if(swsophy.eq.1) then
         ! Tabulated functions (h,theta,k,dthetadh,dkdtheta) tabulated
         do node = 1,numnod
          numtab(node) = numtablay(layer(node))
          do i=0,matabentries
             ientrytab(node,i) = ientrytablay(layer(node),i)
          end do
        end do
        do node = 1,numnod
          do i = 1,7
            do j = 1, numtab(node)
              sptab(i,node,j) = sptablay(i,layer(node),j)
            end do
          end do
          ! Assign values to cofgen — [SS-SWC S-2.12B] legacy half-writes dropped
          sw%cofgen(1,node) = 0.0_real64                  ! thetar
          sw%cofgen(2,node) = sptab(2,node,numtab(node))  ! thetas
          sw%cofgen(3,node) = sptab(3,node,numtab(node))  ! ksat
          if (do_ln_trans) sw%cofgen(3,node) = dexp(sw%cofgen(3,node))
          ! [SS-GR-UTILS] Mirror tabulated cofgen writes into typed vg_params
          sw%vg_params(node)%thetar = sw%cofgen(1,node)
          sw%vg_params(node)%thetas = sw%cofgen(2,node)
          sw%vg_params(node)%ksat   = sw%cofgen(3,node)
          ! [SS-GR-UTILS Task 4] Mirror numtab/ientrytab/sptab into state (swsophy=1 only).
          sw%numtab(node) = numtab(node)
          do i = 0, matabentries
             sw%ientrytab(node,i) = ientrytab(node,i)
          end do
          do i = 1, 7
             do j = 1, numtab(node)
                sw%sptab(i,node,j) = sptab(i,node,j)
             end do
          end do
        end do
        do lay = 1,numlay
          ksatfit(lay) = sw%vg_params(nod1lay(lay))%ksat   ! [SS-SWC S-2.3]
          sw%thetsl(lay) = sw%vg_params(nod1lay(lay))%thetas  ! [SS-SWC S-1.3/S-2.12B]
        end do
      else
         ! MvanG functions
         do node = 1,numnod
          lay = layer(node)
          do i = 1, 10
            sw%cofgen(i,node) = paramvg(i,lay)            ! [SS-SWC S-1.3/S-2.12B]
          end do
          ! Assign dummy value to alphaw
          sw%cofgen(8,node) = -9999.9d0                   ! [SS-SWC S-1.3/S-2.12B]
          if (sw%vg_params(node)%ksatexm > 0.0d0) sw%fluseksatexm(node) = .true.  ! [SS-SWC S-2.12B]
          sw%cofgen(11,node) = relsatthr(lay)             ! [SS-SWC S-1.3/S-2.12B]
          sw%cofgen(12,node) = ksatthr(lay)               ! [SS-SWC S-1.3/S-2.12B]
          if (iHWCKmodel(lay) ==  3 .OR. iHWCKmodel(lay) ==  6 .OR. iHWCKmodel(lay) ==  7 .OR. &
              iHWCKmodel(lay) == 10 .OR. iHWCKmodel(lay) == 11) then
             sw%cofgen(13:17,node) = paramvg(13:17,lay)   ! [SS-SWC S-1.3/S-2.12B]
          end if
          if (iHWCKmodel(lay) ==  5 .OR. iHWCKmodel(lay) ==  7) then
             sw%cofgen(18,node) = paramvg(18,lay)         ! [SS-SWC S-1.3/S-2.12B]
          end if
          if (iHWCKmodel(lay) ==  8 .OR. iHWCKmodel(lay) ==  9 .OR. &
              iHWCKmodel(lay) == 10 .OR. iHWCKmodel(lay) == 11) then
             sw%cofgen(18:21,node) = paramvg(18:21,lay)   ! [SS-SWC S-1.3/S-2.12B]
          end if
          ! [SS-GR-UTILS] Mirror analytical cofgen writes into typed vg_params
          sw%vg_params(node)%thetar          = sw%cofgen(1,node)
          sw%vg_params(node)%thetas          = sw%cofgen(2,node)
          sw%vg_params(node)%ksat            = sw%cofgen(3,node)
          sw%vg_params(node)%alpha           = sw%cofgen(4,node)
          sw%vg_params(node)%lpar            = sw%cofgen(5,node)
          sw%vg_params(node)%npar            = sw%cofgen(6,node)
          sw%vg_params(node)%mpar            = sw%cofgen(7,node)
          sw%vg_params(node)%alphaw_sentinel = sw%cofgen(8,node)
          sw%vg_params(node)%h_enpr          = sw%cofgen(9,node)
          sw%vg_params(node)%ksatexm         = sw%cofgen(10,node)
          sw%vg_params(node)%relsatthr       = sw%cofgen(11,node)
          sw%vg_params(node)%ksatthr         = sw%cofgen(12,node)
          sw%vg_params(node)%alpha_2         = sw%cofgen(13,node)
          sw%vg_params(node)%npar_2          = sw%cofgen(14,node)
          sw%vg_params(node)%mpar_2          = sw%cofgen(15,node)
          sw%vg_params(node)%omega_1         = sw%cofgen(16,node)
          sw%vg_params(node)%omega_2         = sw%cofgen(17,node)
          sw%vg_params(node)%h0              = sw%cofgen(18,node)
          sw%vg_params(node)%ha              = sw%cofgen(19,node)
          sw%vg_params(node)%apar            = sw%cofgen(20,node)
          sw%vg_params(node)%omega_k         = sw%cofgen(21,node)
        end do
        sw%thetsl = 0.0_real64                            ! [SS-SWC S-1.3/S-2.12B]
        do lay = 1, numlay
          sw%thetsl(lay) = paramvg(2,lay)                 ! [SS-SWC S-1.3/S-2.12B]
        end do
      endif

! --- saturated and residual watercontent of each node; hysteresis parameters
      do node = 1,numnod
        lay = layer(node)
        sw%thetar(node) = sw%vg_params(node)%thetar        ! [SS-SWC S-1.3/S-2.12B]
        sw%thetas(node) = sw%vg_params(node)%thetas        ! [SS-SWC S-1.3/S-2.12B]
        !!! Kroes: disable combi of swsophy=1 and swhyst=1
        if (swhyst.eq.1) then
           ! Wetting curve
           sw%indeks(node) = 1                            ! [SS-SWC S-1.3/S-2.12B]
           sw%cofgen(4,node) = paramvg(8,lay)             ! [SS-SWC S-1.3/S-2.12B]
        elseif (swhyst.eq.0.or.swhyst.eq.2) then
           ! Drying branch or simulation without hysteresis
           sw%indeks(node) = -1                           ! [SS-SWC S-1.3/S-2.12B]
           sw%cofgen(4,node) = paramvg(4,lay)             ! [SS-SWC S-1.3/S-2.12B]
        endif
      end do

      if (swinco.eq.1) then
         ! Pressure head profile is input
         ! [SS-SWC S-2.12B] legacy h(:) reads/writes retargeted to state%soilwater%h
         do i = 1, nhead
          tab(i*2) = sw%h(i)
          tab(i*2-1) = abs(zi(i))
        end do
        do i = 1, numnod
          sw%h(i) = afgen(tab,macp*2,abs(z(i)))            ! [SS-SWC S-1.3/S-2.12B]
        end do
      endif
      if (swinco.eq.2 .and. swbotb.ne.8) then
        if (abs(gwli-(z(numnod)-0.5d0*dz(numnod))) .lt.1.0d-4) then
          messag = 'Initial groundwaterlevel (SWINCO=2) is '//          &
     &    'too close to bottom of soil profile'//                       &
     &    ' must be corrected!'
          call fatalerr_collected ('soilwater',messag)
        endif
      endif
      if (swinco.eq.3) then
        if (nhead.ne.numnod) then
          messag = 'Initial data are read from file (SWINCO=3) and '//  &
     &    'number of nodes/compartments is not consistent with NUMNOD'//&
     &    'must be corrected!'
          call fatalerr_collected ('soilwater',messag)
        endif
      endif
      if (swinco.eq.1.or.swinco.eq.3) then
         ! Determine groundwater level — [SS-SWC S-2.12B] all legacy half-writes dropped
         if (sw%h(numnod) .gt. -1.d-5) then
          i = numnod
          do while ((sw%h(i) .gt. -1.d-5) .and. (i .gt. 1))
              i = i - 1
          end do
          if (sw%h(i) .lt. -1.d-5) then
            sw%gwl = z(i+1) + sw%h(i+1) / (sw%h(i+1) - sw%h(i)) * (z(i) - z(i+1))  ! [SS-SWC S-1.3/S-2.12B]
            ! Assume hydrostatic equilibrium in saturated part
            do j = i+1, numnod
              sw%h(j) = sw%gwl - z(j)                     ! [SS-SWC S-1.3/S-2.12B]
            end do
          endif
        endif
      else
         ! Pressure head profile is calculated from groundwater level
         if (swbotb.eq.1) then
          sw%gwl = afgen (gwltab,mabbc*2,state%timecontrol%t1900+state%timecontrol%dt-1.d0)   ! [SS-SWC S-1.3/S-2.12B] [TC-8]

          if(abs(sw%gwl-(z(numnod)-0.5d0*dz(numnod))) .lt.1.0d-4) then
          messag = 'Groundwaterlevel as bottom boundary (SWBOTB=1) is'//&
     &    'below or to close to bottom of soil profile'//               &
     &    ' must be corrected!'
            call fatalerr_collected ('soilwater',messag)
          endif
        else
          sw%gwl = gwli                                   ! [SS-SWC S-1.3/S-2.12B]
        endif
        if (sw%gwl.gt.0.0d0) then
          sw%pond = sw%gwl                                ! [SS-SWC S-1.3/S-2.12B]
        else
          sw%pond = 0.0_real64                            ! [SS-SWC S-1.3/S-2.12B]
        endif
        do i = 1,numnod
          sw%h(i) = sw%gwl - z(i)                         ! [SS-SWC S-1.3/S-2.12B]
        end do
      endif

      ! In case of preferential flow, adjust Van Genuchten parameters
      do i = 1, numnod
        sw%theta(i) = watcon(sw%h(i), &
                              sw%vg_params(i), &
                              sw%iHWCKmodel(sw%layer(i)), &
                              i, state%soilwater)                  ! [SS-SWC S-1.3/S-2.12B] [SS-GR-UTILS Task 5]
      end do

      ! Hydraulic conductivities, differential moisture capacities
      ! and mean hydraulic conductivities for each node
      do node = 1,numnod
        sw%dimoca(node) = moiscap(sw%h(node), &
                                  sw%vg_params(node), &
                                  sw%iHWCKmodel(sw%layer(node)), &
                                  state%timecontrol%dt, &
                                  node, state%soilwater)             ! [SS-SWC S-1.3/S-2.12B] [SS-GR-UTILS Task 7]

        sw%FrArMtrx(node) = 1.0_real64                    ! [SS-SWC S-1.3/S-2.12B]
        sw%k(node) = hconduc(sw%h(node),sw%theta(node),state%heat%rfcp(node),state%heat%tsoil(node), &
                             sw%vg_params(node), &
                             sw%iHWCKmodel(sw%layer(node)), &
                             sw%fluseksatexm(node), &
                             node, state%soilwater)                            ! [SS-SWC S-1.3/S-2.12B] [SS-GR-UTILS Task 6]

        if(node.gt.1) sw%kmean(node) = hcomean(swkmean,sw%k(node-1),sw%k(node),dz(node-1),dz(node))  ! [SS-SWC S-1.3/S-2.12B]
      end do
      sw%kmean(numnod+1) = sw%k(numnod)                   ! [SS-SWC S-1.3/S-2.12B]

      ! Initial soil water storage
      do i = 1, NumNod
         sw%FrArMtrx(i) = 1.0_real64                      ! [SS-SWC S-1.3/S-2.12B]
      enddo
      sw%volact = 0.0_real64                              ! [SS-SWC S-1.3/S-2.12B]
      call watstor (state)                                ! [SS-SWC S-1.6] state arg added
      sw%volini = sw%volact                               ! [SS-SWC S-1.3/S-2.12B]
      sw%pondini = sw%pond                                ! [SS-SWC S-1.3/S-2.12B]

      ! Initial groundwater level
      call calcgwl (state)

      call log_info('soilwater', 'Soil state initialized: gwl=' // to_str(real(sw%gwl,4)) // &  ! [SS-SWC S-2.3]
                    ' cm, numnod=' // to_str(numnod) // ', numlay=' // to_str(numlay) // &
                    ', volini=' // to_str(real(sw%volini,4)) // ' cm')              ! [SS-SWC S-2.3]

         end associate  ! sw => state%soilwater [SS-SWC S-1.3]

      return

      case (2)

         ! Calculate Soilwater rate/state variables

         ! Reset intermediate soil water fluxes — [SS-SWC S-2.12B] reset() handles all
         if (state%timecontrol%flDayStart) then  ! [TC-8]
          call state%soilwater%reset_intermediate_per_day()              ! [SS-SWC S-2.1]
      end if

      if (state%timecontrol%flZeroIntr) then
        call state%soilwater%reset_intermediate()                        ! [SS-SWC S-2.1]

        state%soilwater%IPondBeg = state%soilwater%pond     ! [SS-SWC S-1.4b/S-2.3/S-2.12B]
        do node = 1, numnod
          state%soilwater%IThetaBeg(node) = state%soilwater%theta(node)  ! [SS-SWC S-1.4b/S-2.3/S-2.12B]
        enddo

        ! [MACRO-RETIRE 2026-05-12] macropore(5,...) retired (ADR 0040).
      endif

      ! Reset cumulative soil water fluxes — [SS-SWC S-2.12B] reset() handles all
      if (state%timecontrol%flZeroCumu) then
        call state%soilwater%reset_cumulative()                         ! [SS-SWC S-2.1]

        ! [MACRO-RETIRE 2026-05-12] macropore(6,...) retired (ADR 0040).

        ! Reset initial water storage and ponding
        state%soilwater%volini = state%soilwater%volact          ! [SS-SWC S-2.12B]
        state%soilwater%pondini = state%soilwater%pond           ! [SS-SWC S-2.12B]
      endif

      ! Save state variables of time = t
      call SoilWaterStateVar(1, state)

      ! Calculate new soil water state variables
      call headcalc(state)

      return

      case (3)

         ! Update hydraulic conductivities to time level t+1
         ! [SS-SWC S-2.12B] all legacy half-writes dropped
         do i = 1,numnod
         state%soilwater%k(i) = hconduc(state%soilwater%h(i),state%soilwater%theta(i), &
                                        state%heat%rfcp(i),state%heat%tsoil(i), &
                                        state%soilwater%vg_params(i), &
                                        state%soilwater%iHWCKmodel(state%soilwater%layer(i)), &
                                        state%soilwater%fluseksatexm(i), &
                                        i, state%soilwater)                    ! [SS-SWC S-1.4b/S-2.12B] [SS-GR-UTILS Task 6]
         if(i.gt.1)then
            state%soilwater%kmean(i) = hcomean(swkmean,state%soilwater%k(i-1),state%soilwater%k(i),dz(i-1),dz(i))  ! [SS-SWC S-1.4b/S-2.12B]
         end if
      enddo
      state%soilwater%kmean(numnod+1) = state%soilwater%k(numnod)  ! [SS-SWC S-1.4b/S-2.12B]

      ! Calculate actual water content of profile
      call watstor (state)                                ! [SS-SWC S-1.6] state arg added

      ! Calculate water fluxes between soil compartments
      ! SS-SWST Phase 2 Task 11: state passed so fluxes() reads qdra/qdrtot from state.
      call fluxes (state)

      ! [MACRO-RETIRE 2026-05-12] macropore(4,...) retired (ADR 0040).

      ! Calculate cumulative fluxes
      call integral (state)

      ! Update parameters for soil water hystereses
      if (swhyst.ne.0) call hysteresis (state)

      case default
         call fatalerr_collected ('SoilWater', 'Illegal value for TASK')
      end select

      return
      end subroutine soilwater

   !> Save and reset soil water state variables
   !!
   !! Manages state variable storage for time stepping. Can save current state
   !! or reset to previous state (useful for adaptive time stepping).
   !!
   !! @param[in] task Task selector: 1=save state, 2=reset state
   !!
   !! @note
   !! Date: January 2007
   !! @endnote
   !!
   subroutine SoilWaterStateVar(task, state)

! --- global variables
      use Variables
      use swap_state_mod, only: swap_state_t
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none

      ! Arguments
      integer task
      type(swap_state_t), intent(inout) :: state

      ! Local variables
      integer i

      ! [SS-SWC S-2.12B] legacy h/hm1/theta/thetm1/gwl/gwlm1/pond/pondm1/k/kmean half-writes dropped
      select case (task)
      case (1)

         ! Save state variables of time = t
         do i = 1,numnod
        state%soilwater%hm1(i)    = state%soilwater%h(i)
        state%soilwater%thetm1(i) = state%soilwater%theta(i)
      enddo
      state%soilwater%gwlm1  = state%soilwater%gwl
      state%soilwater%pondm1 = state%soilwater%pond

      return

      case (2)

         ! Reset soil state variables
         do i = 1,numnod
        state%soilwater%h(i)     = state%soilwater%hm1(i)
        state%soilwater%theta(i) = state%soilwater%thetm1(i)
      enddo
      state%soilwater%kmean(numnod+1) = state%soilwater%k(numnod)
      state%soilwater%gwl  = state%soilwater%gwlm1
      state%soilwater%pond = state%soilwater%pondm1

      case default
         call fatalerr_collected ('SoilWaterStateVar', 'Illegal value for TASK')
      end select

      return
      end subroutine SoilWaterStateVar

   !> Check for hysteretic reversal and update model parameters
   !!
   !! Checks whether the wetting/drying direction has changed and updates
   !! Van Genuchten parameters accordingly for hysteresis modeling.
   !!
   !! @note
   !! Date: 23/10/2000
   !! @endnote
   !!
   subroutine hysteresis (state)
      ! [SS-SWC S-2.12B] h/hm1/indeks/cofgen/dimoca/theta retired — read via state%soilwater
      use variables, only: numnod,layer,tau,paramvg,disnod
      use soilhydraulics_utils, only: moiscap, prhead
      use swap_array_dimensions, only: macp
      use swap_state_mod, only: swap_state_t
      use hydraulic_params_mod, only: vanGenuchten_params_t

      implicit none

      ! Arguments
      type(swap_state_t), intent(inout) :: state

      ! Local variables
      integer node,lay,indtem(macp)
      real(8) delp,sew,sed,fvalue
      real(8) thetar(macp),thetas(macp),alfamg(macp)
      type(vanGenuchten_params_t) :: vg_hys

      ! [SS-SWC S-2.3] reader cutover: read h/hm1/theta/indeks/cofgen/dimoca from state
      associate( &
         sw_h      => state%soilwater%h,       &
         sw_hm1    => state%soilwater%hm1,     &
         sw_theta  => state%soilwater%theta,   &
         sw_indeks => state%soilwater%indeks,  &
         sw_cofgen => state%soilwater%cofgen,  &
         sw_dimoca => state%soilwater%dimoca   &
      )

      ! Check for reversal
      do node = 1,numnod
        delp = sw_hm1(node)-sw_h(node)                          ! [SS-SWC S-2.3]
        if (delp/float(sw_indeks(node)).gt.tau.and.              &  ! [SS-SWC S-2.3]
     &     sw_h(node).lt.-10.0d0 .and. sw_h(node).gt.-1.0d3) then  ! [SS-SWC S-2.3]
          indtem(node) = -sw_indeks(node)                       ! [SS-SWC S-2.3]
        else
          indtem(node) = sw_indeks(node)                        ! [SS-SWC S-2.3]
        endif
      end do

      ! Change parameters scanning curves
      do 100 node = 1,numnod
        lay = layer(node)

         ! No change
         if (indtem(node).eq.sw_indeks(node) .or.               &  ! [SS-SWC S-2.3]
     &     abs(paramvg(4,lay)-paramvg(8,lay)) .lt. 1.d-4) goto 100

         ! Relative saturation
         sew = (1.0d0+(paramvg(8,lay)*(-sw_h(node)))**paramvg(6,lay))   &  ! [SS-SWC S-2.3]
     &        **(-paramvg(7,lay))
        sed = (1.0d0+(paramvg(4,lay)*(-sw_h(node)))**paramvg(6,lay))   &  ! [SS-SWC S-2.3]
     &        **(-paramvg(7,lay))

         ! Change index
         ! [SS-SWC S-2.12B] legacy half-write dropped
         state%soilwater%indeks(node) = -1*sw_indeks(node)      ! [SS-SWC S-1.4a/S-2.12B]

         ! Update alfa, thetar and thetas — [SS-SWC S-2.12B] sw_indeks now points to fresh state value
         if (state%soilwater%indeks(node).eq.1) then
            ! Wetting branch
            alfamg(node) = paramvg(8,lay)
          thetas(node) = paramvg(2,lay)
          thetar(node) = (sw_theta(node)-thetas(node)*sew)/(1.0d0-sew)  ! [SS-SWC S-2.3]

            ! Check on thetar(node) value, if needed correction of h
            fvalue = thetar(node)
          if(thetar(node).lt.paramvg(1,lay)) thetar(node)=paramvg(1,lay)
          if(thetar(node).gt.paramvg(2,lay)) thetar(node)=paramvg(2,lay)
          ! [SS-SWC S-2.12B] legacy cofgen/h half-writes dropped
          state%soilwater%cofgen(4,node) = alfamg(node)         ! [SS-SWC S-1.4a/S-2.12B]
          state%soilwater%cofgen(1,node) = thetar(node)         ! [SS-SWC S-1.4a/S-2.12B]
          state%soilwater%thetar(node)   = thetar(node)         ! [SS-SWC S-1.4a]
          state%soilwater%cofgen(2,node) = thetas(node)         ! [SS-SWC S-1.4a/S-2.12B]
          state%soilwater%thetas(node)   = thetas(node)         ! [SS-SWC S-1.4a]
          if (abs(fvalue-thetar(node)) .gt. 1.d-10) then
             ! Build vg with locally modified thetar/thetas/alpha for wetting branch  [SS-GR-UTILS Task 8]
             vg_hys         = state%soilwater%vg_params(node)
             vg_hys%thetar  = thetar(node)
             vg_hys%thetas  = thetas(node)
             vg_hys%alpha   = alfamg(node)
             state%soilwater%h(node) = prhead(disnod(node), sw_theta(node), sw_h, &
                                              state%soilwater%iHWCKmodel(state%soilwater%layer(node)), &
                                              node, state%soilwater, vg_in=vg_hys)  ! [SS-SWC S-2.12B] [SS-GR-UTILS Task 8]
          endif
        else
           ! Drying branch
           alfamg(node) = paramvg(4,lay)
          thetar(node) = paramvg(1,lay)
          thetas(node) = thetar(node)+(sw_theta(node)-thetar(node))/sed  ! [SS-SWC S-2.3]

            ! Check on thetas(node) value, if needed correction of h
            fvalue = thetas(node)
          if(thetas(node).lt.paramvg(1,lay)) thetas(node)=paramvg(1,lay)
          if(thetas(node).gt.paramvg(2,lay)) thetas(node)=paramvg(2,lay)
          ! [SS-SWC S-2.12B] legacy cofgen/h half-writes dropped
          state%soilwater%cofgen(4,node) = alfamg(node)         ! [SS-SWC S-1.4a/S-2.12B]
          state%soilwater%cofgen(1,node) = thetar(node)         ! [SS-SWC S-1.4a/S-2.12B]
          state%soilwater%thetar(node)   = thetar(node)         ! [SS-SWC S-1.4a]
          state%soilwater%cofgen(2,node) = thetas(node)         ! [SS-SWC S-1.4a/S-2.12B]
          state%soilwater%thetas(node)   = thetas(node)         ! [SS-SWC S-1.4a]
          if (abs(fvalue-thetas(node)) .gt. 1.d-10) then
             ! Build vg with locally modified thetar/thetas/alpha for drying branch  [SS-GR-UTILS Task 8]
             vg_hys         = state%soilwater%vg_params(node)
             vg_hys%thetar  = thetar(node)
             vg_hys%thetas  = thetas(node)
             vg_hys%alpha   = alfamg(node)
             state%soilwater%h(node) = prhead(disnod(node), sw_theta(node), sw_h, &
                                              state%soilwater%iHWCKmodel(state%soilwater%layer(node)), &
                                              node, state%soilwater, vg_in=vg_hys)  ! [SS-SWC S-2.12B] [SS-GR-UTILS Task 8]
          endif
        endif

         ! Update capacity
         ! [SS-SWC S-2.12B] legacy dimoca half-write dropped
         state%soilwater%dimoca(node) = moiscap(sw_h(node), &
                                              state%soilwater%vg_params(node), &
                                              state%soilwater%iHWCKmodel(state%soilwater%layer(node)), &
                                              state%timecontrol%dt, &
                                              node, state%soilwater)  ! [SS-SWC S-1.4b/S-2.12B] [SS-GR-UTILS Task 7]
 100  continue

      end associate  ! sw_h/.../sw_dimoca => state%soilwater [SS-SWC S-2.3]

      return
      end subroutine hysteresis

end module