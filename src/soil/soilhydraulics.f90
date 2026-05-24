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
      ! [SS-GR-FINAL B9] blanket use variables → explicit only-list; all symbols DEFERRED
      ! macp/mabbc → swap_array_dimensions (dimension constants); noddrz: crop root depth node
      use swap_array_dimensions, only: macp, mabbc
      use variables, only: &
         ! [GR-SOIL 2026-05-24] qssdi migrated to state%soilwater — read via soil%qssdi
         ! DEFERRED: swkmean/SwKimpl — hydraulic conductivity averaging switches; Phase C3
         ! DEFERRED: swcaprise — capillary rise prevention switch; Phase C3
         swcaprise, &
         ! DEFERRED: MaxBackTr — max Newton-Raphson backtrack iterations; Phase C3
         ! DEFERRED: fldumpconvcrit — debug convergence dump flag; Phase C3
         fldumpconvcrit, &
         ! DEFERRED: numbit/itnumb — Richards iteration counter/stats; Phase C3/D
         numbit, itnumb, &
         ! DEFERRED: rimlay — drainage resistance (Cauchy BC); Phase C3
         rimlay, &
         ! DEFERRED: sw4 — extra-flux switch for swbotb=3; Phase C3
         sw4, &
         ! DEFERRED: gwlconv — groundwater level convergence criterion; Phase C3
         ! DEFERRED: hplate — lysimeter tensiometer plate head; Phase C3
         hplate, &
         ! DEFERRED: swbotb3Impl — Cauchy BC option; Phase C3
         swbotb3Impl, &
         ! DEFERRED: CritDevh1Cp/CritDevh2Cp/CritDevPondDt — convergence criteria; Phase C3
         ! DEFERRED: flwarn_hc/iwarn_hc — non-convergence warning state; Phase C3
         flwarn_hc, iwarn_hc, &
         ! DEFERRED: noddrz — node index at root zone bottom; Phase C3
         noddrz
      use timestep_control_mod, only: fldecdt
      use swap_log, only: log_warn, log_debug, to_str
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
      ! ndr removed: loop now uses drai%nrlevs (actual drain-level count, always <= Madr=5)

      real(8) hgrad(macp+1), dkdh(macp)

      ! Function and solver variables
      integer indx(macp), ierror
      real(8) a(macp,3), a1(macp,1), b(macp), d, q1
      ! [GR-BH Task 36] q0 and ArMpSs retired from variables.f90 — made local (ADR 0040 complete)
      real(8) q0       ! surface flux for boundary (local; used only in swbotb=1 branch)
      real(8) ArMpSs   ! macropore area fraction at soil surface (local; always 0.d0, ADR 0040)
      logical flok

      ! Note: flwarn_hc, iwarn_hc, nstep_hc moved to variables.f90 module
      ! (previously local SAVE variables - now global for multi-instance support)

      associate (mesh => state%mesh,         &
                 soil => state%soilwater,    &
                 drai => state%drainage,     &
                 heat => state%heat,         &
                 atmo => state%atmosphere,   &
                 time => state%timecontrol,  &
                 swbotb => state%soilwater%swbotb_runtime)

      if (time%flDayStart) then  ! [TC-8]
         flwarn_hc = .true.
         iwarn_hc = 0
      endif
      call dtdpst                                                       &
     &        ('year-month-day,hour:minute:seconds',time%t1900,datetime)  ! [TC-8]

      ! Summation of sink terms (constant for the current time step)
      iBackTr   = 0
      flunsatok(1) = .false.
      flunsatok(2) = .false.
      flunsatok(3) = .false.
      ! SS-DRST Task 3: read qdra from state%drainage (loop bound = drai%nrlevs, not ndr=5)
      do i=1,mesh%numnod
         sink(i) = soil%evp(i)                                    ! [SS-SWC S-2.3] read from state
         if (allocated(drai%qdra)) then
            do j=1,drai%nrlevs
               sink(i) = sink(i) + drai%qdra(j,i)
            end do
         end if
      end do
      source(1:mesh%numnod) = soil%qssdi(1:mesh%numnod)  ! [GR-SOIL 2026-05-24]

      ArMpSs = 0.d0                                            ! macropore retired (ADR 0040)

      ! Groundwater level specified — [SS-SWC S-2.12B] all legacy half-writes dropped
      if(swbotb.eq.1)then
         soil%fllowgwl = .false.               ! [SS-SWC S-1.6/S-2.12B]
         if(soil%gwlinp.ge.mesh%z(1)-1.0d-4)then

            q0 = (atmo%nraidt+atmo%nird+atmo%melt)*(1.0d0-ArMpSs) + soil%runon - soil%reva  ! [SS-ATM/SS-SWC S-2.12B]
            call pondrunoff (state)
            q1 = - q0 + (soil%pond - soil%pondm1)/time%dt + soil%runots / time%dt  ! [SS-SWC S-2.12B] [TC-8]
            soil%theta(1) = watcon(soil%gwlinp, &
                                  soil%vg_params(1), &
                                  soil%iHWCKmodel(soil%layer(1)), &
                                  1, soil)                         ! [SS-SWC S-1.4a/S-2.12B] [SS-GR-UTILS Task 5]
            soil%kmean(1) = hconduc(soil%gwlinp,soil%theta(1),heat%rfcp(1),heat%tsoil(1), &
                                  soil%vg_params(1), &
                                  soil%iHWCKmodel(soil%layer(1)), &
                                  soil%fluseksatexm(1), &
                                  1, soil)                          ! [SS-SWC S-1.4b/S-2.12B] [SS-GR-UTILS Task 6]

            qv(1) = q1
            do i=1,mesh%numnod
               qv(i+1) = qv(i) +mesh%dz(i)*soil%FrArMtrx(i)*(soil%theta(i)-soil%thetm1(i))  &  ! [SS-SWC S-2.12B]
     &                          / time%dt+ sink(i) - source(i) + soil%qrot(i)  ! [TC-8]
            end do
            soil%qbot = qv(mesh%numnod+1)
            soil%h(1) = soil%gwlinp + mesh%disnod(1)*(qv(1)/soil%kmean(1)+1.0d0)  ! [SS-SWC S-1.4a/S-2.12B]
            do i=2,mesh%numnod
               soil%h(i) = soil%h(i-1) + mesh%disnod(i)*(qv(i)/soil%kmean(i)+1.0d0)    ! [SS-SWC S-1.4a/S-2.12B]
            end do

            if(state%cfg%simulation%numerical%swkimpl.eq.1)then
               do i=1,mesh%numnod
                  soil%k(i) = hconduc(soil%h(i),soil%theta(i),heat%rfcp(i),heat%tsoil(i), &
                                    soil%vg_params(i), &
                                    soil%iHWCKmodel(soil%layer(i)), &
                                    soil%fluseksatexm(i), &
                                    i, soil)                        ! [SS-SWC S-1.4b/S-2.12B] [SS-GR-UTILS Task 6]
                  if(i.gt.1)then
                     soil%kmean(i) = hcomean(state%cfg%simulation%numerical%swkmean,soil%k(i-1),soil%k(i),mesh%dz(i-1),mesh%dz(i))  ! [SS-SWC S-1.4b/S-2.12B]
                  end if
               end do
               soil%kmean(mesh%numnod+1) = soil%k(mesh%numnod)                              ! [SS-SWC S-1.4b/S-2.12B]
            end if
            call calcgwl (state)
            return
         else
            NN = 0
            do while (mesh%z(NN+1).gt.soil%gwlinp .and. NN.lt.mesh%numnod)
               NN = NN + 1
            end do
            if (mesh%z(NN+1).lt.(soil%gwlinp+nihil)) then
               ! Groundwater within soil profile
               if ((mesh%z(NN)-soil%gwlinp) .lt. 1.0d-4 .and. (NN.gt.0)) then
                  ! Difference gwlinp with node too small to calculate gradient properly
                  soil%gwlinp = mesh%z(NN)
                  NN = NN-1
               endif
            else
               ! Groundwater below soil profile
               soil%fllowgwl = .true.          ! [SS-SWC S-1.6/S-2.12B]
               soil%hbot = soil%gwlinp - mesh%z(mesh%numnod) + 0.5*mesh%dz(mesh%numnod)
            endif
         end if
      else
         NN = mesh%numnod
      end if

      ! Reset conductivities to time level t

      ! Node nr of compartment with minimized flux of capillary rise
      if (swcaprise) then
         nodncr    = max(5,noddrz)
         flcaprise = .false.
      endif
      do i = 1,mesh%numnod
         soil%k(i) = hconduc(soil%h(i),soil%theta(i),heat%rfcp(i),heat%tsoil(i), &
                           soil%vg_params(i), &
                           soil%iHWCKmodel(soil%layer(i)), &
                           soil%fluseksatexm(i), &
                           i, soil)                                 ! [SS-GR-UTILS Task 6]

         if (swcaprise) then
            ! Prevent capillary rise into the root zone !! special for experts only
            if (i .eq. nodncr) then
               if ((soil%h(i) + mesh%z(i)) .lt. (soil%h(i+1) + mesh%z(i+1))) then  ! negative potential gradient upwards
                  soil%k(i)      = 1.0D-10
                  flcaprise = .true.
               endif
            endif
            if (flcaprise .and. i .eq. nodncr+1) then
                soil%k(i) = 1.0D-10
            endif
         endif
         if(i.gt.1)then
            soil%kmean(i) = hcomean(state%cfg%simulation%numerical%swkmean,soil%k(i-1),soil%k(i),mesh%dz(i-1),mesh%dz(i))
         end if
      enddo
      soil%kmean(mesh%numnod+1) = soil%k(mesh%numnod)

      if(state%cfg%simulation%numerical%swkimpl.eq.0)then
         do i=2,mesh%numnod
            dFdhU(i)   = - soil%kmean(i)  /mesh%disnod(i)
            dFdhL(i-1) = dFdhU(i)
         end do
      end if

      do i=2,NN
         hgrad(i) = (soil%h(i-1)-soil%h(i))/mesh%disnod(i) + 1.0d0
      end do

      F(1) = (soil%theta(1)-soil%thetm1(1))*soil%FrArMtrx(1)*mesh%dz(1)/time%dt + sink(1) - source(1) + soil%qrot(1) + soil%kmean(2) * hgrad(2)  ! [SS-SWC S-2.3] [TC-8]

      call boundtop(state)  ! [SS-HEAT] Task 9: state passed for rfcp access

      ! [MACRO-RETIRE 2026-05-12] MACROPORE(2,...) retired (ADR 0040).

      if (soil%FlRunoff) call pondrunoff (state)

      if(soil%ftoph)then
         hgrad(1) = (soil%hsurf-soil%h(1))/mesh%disnod(1) + 1.d0
         F(1)     = F(1) - soil%kmean(1) * hgrad(1)
      else
         F(1) = F(1) + soil%qtop
      end if

      do i=2,NN-1
         F(i) = (soil%theta(i)-soil%thetm1(i))*soil%FrArMtrx(i)*mesh%dz(i)/time%dt + sink(i) - source(i) +     &  ! [SS-SWC S-2.3] [TC-8]
     &          soil%qrot(i) - soil%kmean(i) * hgrad(i) + soil%kmean(i+1) * hgrad(i+1)
      end do

      if(swbotb.eq.1 .and. (.not.soil%fllowgwl))then
         hgrad(NN+1) = soil%h(NN)/(mesh%z(nn)-soil%gwlinp) + 1.0d0
      else if(swbotb.eq.5 .or. (swbotb.eq.1 .and. soil%fllowgwl))then
         hgrad(NN+1) = (soil%h(NN) - soil%hbot) / mesh%disnod(NN+1)  + 1.0d0
      else if(swbotb.eq.8 .and. soil%h(NN).gt. Critdz - mesh%disnod(NN+1) + hplate) then
         hgrad(NN+1) = (soil%h(NN) - hplate) / mesh%disnod(NN+1)  + 1.0d0
         flboth = .true.
      else
         flboth = .false.
      end if

      if(swbotb.eq.1 .and. (.not.soil%fllowgwl))then
         soil%theta(NN)= watcon(soil%h(NN), &
                              soil%vg_params(NN), &
                              soil%iHWCKmodel(soil%layer(NN)), &
                              NN, soil)                            ! [SS-GR-UTILS Task 5]
         soil%k(NN)    = hconduc(soil%h(NN),soil%theta(NN),heat%rfcp(NN),heat%tsoil(NN), &
                              soil%vg_params(NN), &
                              soil%iHWCKmodel(soil%layer(NN)), &
                              soil%fluseksatexm(NN), &
                              NN, soil)                             ! [SS-GR-UTILS Task 6]
         soil%kmean(NN+1) = hcomean(state%cfg%simulation%numerical%swkmean,soil%k(NN),soil%vg_params(NN+1)%ksat, &  ! [SS-GR-UTILS Task 15]
     &                        mesh%dz(NN),mesh%dz(NN+1))
         F(NN) = (soil%theta(NN) - soil%thetm1(NN))*soil%FrArMtrx(NN)*mesh%dz(NN)/time%dt +      &  ! [SS-SWC S-2.3] [TC-8]
     &           sink(NN)-source(NN)+soil%qrot(NN)-soil%kmean(NN)*hgrad(NN) +soil%kmean(NN+1)*hgrad(NN+1)
      else

         F(NN) = (soil%theta(NN) - soil%thetm1(NN))*soil%FrArMtrx(NN)*mesh%dz(NN)/time%dt        &  ! [SS-SWC S-2.3] [TC-8]
     &         - soil%kmean(NN) * hgrad(NN) + sink(NN) - source(NN) + soil%qrot(NN)

         if(swbotb.eq.3.and.swbotb3Impl.eq.1)then ! Cauchy-relation, implemented as head boundary
            if (soil%swbotb3resvert.eq.0) then
               soil%qbot = - (soil%h(NN)+mesh%z(NN)-soil%deepgw) / (mesh%disnod(NN+1)/soil%kmean(NN+1)+rimlay)
            elseif (soil%swbotb3resvert.eq.1) then
               soil%qbot = - (soil%h(NN)+mesh%z(NN)-soil%deepgw) / rimlay
            endif
! ---       extra groundwater flux might be added
            if (sw4 .eq. 1) soil%qbot = soil%qbot + afgen(soil%qbotab,mabbc*2,time%t1900+time%dt)  ! [TC-8]
            F(NN) = F(NN) - soil%qbot
         else if(swbotb.eq.5 .or. (swbotb.eq.1 .and. soil%fllowgwl))then ! pressure head at lower boundary specified
            F(NN) = F(NN) + soil%kmean(NN+1) * hgrad(NN+1)
         else if(swbotb.eq.7 .or. swbotb .eq. -2)then ! free drainage option
            soil%kmean(mesh%numnod+1) = hconduc(soil%h(mesh%numnod),soil%theta(mesh%numnod),heat%rfcp(mesh%numnod),heat%tsoil(mesh%numnod), &
                                         soil%vg_params(mesh%numnod), &
                                         soil%iHWCKmodel(soil%layer(mesh%numnod)), &
                                         soil%fluseksatexm(mesh%numnod), &
                                         mesh%numnod, soil)              ! [SS-GR-UTILS Task 6]
            soil%qbot = -1.0d0 * soil%kmean(mesh%numnod+1)
            F(NN) = F(NN) - soil%qbot
         ! Lysimeter option
         else if(swbotb.eq.8)then
            if (flboth) then
               soil%hbot = hplate
               F(NN) = F(NN) + soil%kmean(NN+1) * hgrad(NN+1)
            else
               soil%qbot = 0.0d0
            end if
         ! Flux bottom boundary
         else
            F(NN) = F(NN) - soil%qbot
         end if

      end if

      ! Initial estimate of F inner product
      sumold    = 0.0d0
      do i=1,NN
         sumold = sumold + F(i)*F(i)
      end do
      sumold = 0.5d0 * sumold

      ! Start iteration loop, time%MaxIt specified in the input
      if(time%fldtmin)then  ! [TC-8]
         MaxIt1 = 2*time%MaxIt
      else
         MaxIt1 = time%MaxIt
      end if

      sum= 0.d0   ! For Forcheck
      Do numbit = 1,MaxIt1

         do i = 1, NN
            ! Save values of h
            hold(i)   = soil%h(i)

            ! Derivative of theta to h (differential moisture capacity),
            ! as part of main diagonal
            soil%dimoca(i) = moiscap(soil%h(i), &
                                    soil%vg_params(i), &
                                    soil%iHWCKmodel(soil%layer(i)), &
                                    time%dt, &
                                    i, soil)              ! [SS-GR-UTILS Task 7]
            soil%dimoca(i) = soil%dimoca(i)                              ! [SS-SWC S-1.4b]

         enddo

         if(state%cfg%simulation%numerical%swkimpl.eq.1)then
            do i = 1, NN
               dkdh(i) = dhconduc(soil%h(i),soil%theta(i),soil%dimoca(i),heat%rfcp(i), &
                                   soil%vg_params(i), &
                                   soil%iHWCKmodel(soil%layer(i)), &
                                   i, soil)                         ! [SS-GR-UTILS Task 6]
            enddo
            do i=2,NN
               dFdhU(i)   = - soil%kmean(i)  /mesh%disnod(i)
               dFdhL(i-1) = dFdhU(i)
            end do
         end if

         if (swcaprise .and. flcaprise) then
           dkdh(nodncr+1) = 1.0D-30
         endif

         ! Jacobian matrix elements
         dFdhM(1) = soil%dimoca(1)*soil%FrArMtrx(1)*mesh%dz(1)/time%dt - dFdhL(1)  ! [SS-SWC S-2.3] [TC-8]
         ! If the head boundary condition applies: add the k1/(0.5*dz1) term
         ! to the first element of the main diagonal
         if(soil%ftoph) dFdhM(1) = dFdhM(1) + soil%kmean(1)/mesh%disnod(1)

         do i=2,NN-1
            dFdhM(i) = soil%dimoca(i)*soil%FrArMtrx(i)*mesh%dz(i)/time%dt - dFdhU(i)        &  ! [SS-SWC S-2.3] [TC-8]
     &                                                - dFdhL(i)
         end do

         dFdhM(NN) = soil%dimoca(NN)*soil%FrArMtrx(NN)*mesh%dz(NN)/time%dt - dFdhU(NN)  ! [SS-SWC S-2.3] [TC-8]
         if(swbotb.eq.1 .and. (.not.soil%fllowgwl))then
            dFdhM(NN) = dFdhM(NN) + soil%kmean(NN+1)/(mesh%z(NN)-soil%gwlinp)
         else if(swbotb.eq.3.and.swbotb3Impl.eq.1)then ! Cauchy
            if (soil%swbotb3resvert.eq.0) then
               dFdhM(NN) = dFdhM(NN) + 1.0d0 /                          &
     &                              (mesh%disnod(NN+1)/soil%kmean(NN+1)+rimlay)   
            elseif (soil%swbotb3resvert.eq.1) then
               dFdhM(NN) = dFdhM(NN) + 1.0d0 / rimlay
            endif
         else if(swbotb.eq.5 .or. (swbotb.eq.1 .and. soil%fllowgwl))then
            dFdhM(NN) = dFdhM(NN) + soil%kmean(NN+1)/mesh%disnod(NN+1)         
         else if(swbotb.eq.7 .or. swbotb .eq. -2)then ! implicitly: soil%kmean(NN+1)
            dFdhM(NN) = dFdhM(NN) + dkdh(NN) * 0.5d0
         else if(swbotb.eq.8 .and. flboth)then
            dFdhM(NN) = dFdhM(NN) + soil%kmean(NN+1)/mesh%disnod(NN+1)
         end if

         if(state%cfg%simulation%numerical%swkimpl.eq.1)then
            dFdhM(1) = dFdhM(1) + dkdh(1) * hgrad(2) *                  &
     &                 dkmean(state%cfg%simulation%numerical%swkmean,soil%k(1),soil%k(2),mesh%dz(1),mesh%dz(2))
            if(soil%ftoph) dFdhM(1) = dFdhM(1) - dkdh(1) * hgrad(1) * 0.5d0
            dFdhL(1) = dFdhL(1) + dkdh(2) * hgrad(2) *                  &
     &                 dkmean(state%cfg%simulation%numerical%swkmean,soil%k(2),soil%k(1),mesh%dz(2),mesh%dz(1)) 
            do i=2,NN-1
               dFdhU(i) = dFdhU(i) - dkdh(i-1) * hgrad(i) *             &
     &                    dkmean(state%cfg%simulation%numerical%swkmean,soil%k(i-1),soil%k(i),mesh%dz(i-1),mesh%dz(i)) 
               dFdhM(i) = dFdhM(i) - dkdh(i) * hgrad(i) *               &
     &                    dkmean(state%cfg%simulation%numerical%swkmean,soil%k(i),soil%k(i-1),mesh%dz(i),mesh%dz(i-1))     &
     &                             + dkdh(i) * hgrad(i+1) *             &
     &                    dkmean(state%cfg%simulation%numerical%swkmean,soil%k(i),soil%k(i+1),mesh%dz(i),mesh%dz(i+1))
               dFdhL(i) = dFdhL(i) + dkdh(i+1) * hgrad(i+1) *           &
     &                    dkmean(state%cfg%simulation%numerical%swkmean,soil%k(i+1),soil%k(i),mesh%dz(i+1),mesh%dz(i)) 
            end do
            dFdhU(NN) = dFdhU(NN) - dkdh(NN-1) * hgrad(NN) *            &
     &                  dkmean(state%cfg%simulation%numerical%swkmean,soil%k(NN-1),soil%k(NN),mesh%dz(NN-1),mesh%dz(NN)) 
            dFdhM(NN) = dFdhM(NN) - dkdh(NN) * hgrad(NN) *              &
     &                  dkmean(state%cfg%simulation%numerical%swkmean,soil%k(NN),soil%k(NN-1),mesh%dz(NN),mesh%dz(NN-1))

            if(swbotb.eq.1 .or. swbotb.eq.5 .or. swbotb.eq.8            &
     &         .and. flboth)then
               dFdhM(NN) = dFdhM(NN) + 0.5d0 * dkdh(NN) * hgrad(NN+1)
            end if
         end if

         ! [MACRO-RETIRE 2026-05-12] MACROPORE(3,...) retired (ADR 0040).

         ! Solve the tridiagonal matrix
         call tridag(NN, dFdhU, dFdhM, dFdhL, F, difh, ierror)

         if(ierror.ne.0)then
            call dtdpst ('year-month-day',time%t1900+1.001d0,datetmp)  ! [TC-8]
            messag = ' Tri-band matrix in HeadCalc appeared to be'//    &
     &               ' singular at '//datetmp//                            &
     &               '   Alternative SOLVER chosen'
            call log_warn('Headcalc', messag)
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
         do itry = 1,state%cfg%simulation%numerical%MaxBackTr
            iBackTr = iBackTr + 1
            ! Factor reduces the change of h (difh) calculated as a full
            ! Newton Raphson step
            if(time%fldtmin .and. numbit.gt.time%MaxIt)then              ! [TC-8]
               factmax = 0.0d0
               do i = 1,NN
                  if(dabs( hold(i) ) .lt. 1.0d0 )then
                     factmax = max( factmax, dabs( difh(i) ) ) 
                  else
                     factmax = max( factmax, dabs( difh(i) / hold(i) ) )
                  end if
               end do
               do i = 1,NN
                  soil%h(i) = hold(i) - difh(i) * min(1.0d0, 1.0d0 / factmax)
                  soil%h(i) = soil%h(i)                                ! [SS-SWC S-1.4a]
               end do
            else
               do i = 1,NN
                  soil%h(i) = hold(i) - factor * difh(i)
                  soil%h(i) = soil%h(i)                                ! [SS-SWC S-1.4a]
               end do
            end if
            do i = 1,NN
              soil%theta(i) = watcon(soil%h(i), &
                                   soil%vg_params(i), &
                                   soil%iHWCKmodel(soil%layer(i)), &
                                   i, soil)            ! [SS-GR-UTILS Task 5]
              soil%theta(i) = soil%theta(i)                            ! [SS-SWC S-1.4a]
            enddo
            do i=2,NN
               hgrad(i) = (soil%h(i-1)-soil%h(i))/mesh%disnod(i) + 1.0d0
            end do


            if(state%cfg%simulation%numerical%swkimpl.eq.1)then
               call Rootextraction(state)
               do i = 1,NN
                  soil%k(i) = hconduc(soil%h(i),soil%theta(i),heat%rfcp(i),heat%tsoil(i), &
                                    soil%vg_params(i), &
                                    soil%iHWCKmodel(soil%layer(i)), &
                                    soil%fluseksatexm(i), &
                                    i, soil)                        ! [SS-GR-UTILS Task 6]
                  soil%k(i) = soil%k(i)                                  ! [SS-SWC S-1.4b]
                  if(i.gt.1)then
                     soil%kmean(i)=hcomean(state%cfg%simulation%numerical%swkmean,soil%k(i-1),soil%k(i),mesh%dz(i-1),mesh%dz(i))
                     soil%kmean(i) = soil%kmean(i)                        ! [SS-SWC S-1.4b]
                  end if
               end do
               soil%kmean(NN+1) = soil%k(NN)
               soil%kmean(NN+1) = soil%k(NN)                              ! [SS-SWC S-1.4b]
            end if

            ! Prevent capillary rise into the root zone !! special for experts
            if (swcaprise) then
               flcaprise = .false.
               i = nodncr
               if ((soil%h(i) + mesh%z(i)) .lt. (soil%h(i+1) + mesh%z(i+1))) then  ! negative potential gradient upwards
                  soil%k(i)      = 1.0D-10
                  flcaprise = .true.
               endif
               if (flcaprise) then
                 soil%k(i+1) = 1.0D-10
               endif
               if(i.gt.1)then
                 soil%kmean(i)=hcomean(state%cfg%simulation%numerical%swkmean,soil%k(i-1),soil%k(i),mesh%dz(i-1),mesh%dz(i))
                 soil%kmean(i) = soil%kmean(i)                            ! [SS-SWC S-1.4b]
                 soil%kmean(i+1)=hcomean(state%cfg%simulation%numerical%swkmean,soil%k(i),soil%k(i+1),mesh%dz(i),mesh%dz(i+1))
                 soil%kmean(i+1) = soil%kmean(i+1)                        ! [SS-SWC S-1.4b]
               end if
               soil%k(i) = soil%k(i)                                      ! [SS-SWC S-1.4b] soil%k(nodncr)
               if (flcaprise) soil%k(i+1) = soil%k(i+1)                   ! [SS-SWC S-1.4b] soil%k(nodncr+1)
            endif

            ! Calculate F-function
            F(1) = (soil%theta(1) - soil%thetm1(1))*soil%FrArMtrx(1)*mesh%dz(1)/time%dt + sink(1) - source(1) &  ! [SS-SWC S-2.3] [TC-8]
     &           + soil%qrot(1) + soil%kmean(2) * hgrad(2)

            ! [MACRO-RETIRE 2026-05-12] FlMacropore QMpLatSsSav save retired (ADR 0040).

            call boundtop(state)  ! [SS-HEAT] Task 9: state passed for rfcp access

            ! [MACRO-RETIRE 2026-05-12] MACROPORE(2,...) retired (ADR 0040).

            if (soil%FlRunoff) call pondrunoff (state)

            if(soil%ftoph)then
               hgrad(1) = (soil%hsurf-soil%h(1))/mesh%disnod(1) + 1.d0
               F(1) = F(1) - soil%kmean(1) * hgrad(1)
            else
               F(1) = F(1) + soil%qtop
            end if

            do i=2,NN-1
               F(i) = (soil%theta(i)-soil%thetm1(i))*soil%FrArMtrx(i)*mesh%dz(i)/time%dt + sink(i) - source(i) &  ! [SS-SWC S-2.3] [TC-8]
     &              + soil%qrot(i) - soil%kmean(i)*hgrad(i)+soil%kmean(i+1)*hgrad(i+1)
            end do

            if(swbotb.eq.1 .and. (.not.soil%fllowgwl))then
               hgrad(NN+1) = soil%h(NN)/(mesh%z(nn)-soil%gwlinp) + 1.0d0
           else if(swbotb.eq.5 .or. (swbotb.eq.1 .and. soil%fllowgwl))then
               hgrad(NN+1) = (soil%h(NN) - soil%hbot) / mesh%disnod(NN+1)  + 1.0d0
            else if(swbotb.eq.8 .and. flboth)then
               hgrad(NN+1) = (soil%h(NN) - hplate) / mesh%disnod(NN+1)  + 1.0d0
            end if

            if(swbotb.eq.1 .and. (.not.soil%fllowgwl))then
               soil%theta(NN) = watcon(soil%h(NN), &
                                     soil%vg_params(NN), &
                                     soil%iHWCKmodel(soil%layer(NN)), &
                                     NN, soil)          ! [SS-GR-UTILS Task 5]
               soil%theta(NN) = soil%theta(NN)                         ! [SS-SWC S-1.4a]
               soil%k(NN)     = hconduc(soil%h(NN),soil%theta(NN),heat%rfcp(NN),heat%tsoil(NN), &
                                      soil%vg_params(NN), &
                                      soil%iHWCKmodel(soil%layer(NN)), &
                                      soil%fluseksatexm(NN), &
                                      NN, soil)                     ! [SS-GR-UTILS Task 6]
               soil%k(NN) = soil%k(NN)                                 ! [SS-SWC S-1.4b]
               soil%kmean(NN+1) = hcomean(state%cfg%simulation%numerical%swkmean,soil%k(NN),soil%vg_params(NN+1)%ksat &  ! [SS-GR-UTILS Task 15]
     &                       ,mesh%dz(NN),mesh%dz(NN+1))
               soil%kmean(NN+1) = soil%kmean(NN+1)                     ! [SS-SWC S-1.4b]
               F(NN) = (soil%theta(NN) - soil%thetm1(NN))*soil%FrArMtrx(NN)*mesh%dz(NN)/time%dt  &  ! [SS-SWC S-2.3] [TC-8]
     &            - soil%kmean(NN) * hgrad(NN) + soil%kmean(NN+1) * hgrad(NN+1)   &
     &            + sink(NN) - source(NN) + soil%qrot(NN)
            else
               F(NN) = (soil%theta(NN) - soil%thetm1(NN))*soil%FrArMtrx(NN)*mesh%dz(NN)/time%dt  &  ! [SS-SWC S-2.3] [TC-8]
     &               - soil%kmean(NN) * hgrad(NN)                            &
     &               + sink(NN) - source(NN) + soil%qrot(NN)
               if(swbotb.eq.3.and.swbotb3Impl.eq.1)then ! Cauchy
                  if (soil%swbotb3resvert.eq.0) then
                     soil%qbot = - (soil%h(NN)+mesh%z(NN)-soil%deepgw) /       &
     &                                 (mesh%disnod(NN+1)/soil%kmean(NN+1)+rimlay)
                  elseif (soil%swbotb3resvert.eq.1) then
                     soil%qbot = - (soil%h(NN)+mesh%z(NN)-soil%deepgw) / rimlay
                  endif
                  ! Extra groundwater flux might be added
                  if (sw4 .eq. 1) then
                     soil%qbot = soil%qbot + afgen(soil%qbotab,mabbc*2,time%t1900+time%dt)  ! [TC-8]
                  end if
                  F(NN) = F(NN) - soil%qbot
               else if(swbotb.eq.5 .or.(swbotb.eq.1 .and. soil%fllowgwl))then
                  ! Pressure head at lower boundary specified
                  F(NN) = F(NN) + soil%kmean(NN+1) * hgrad(NN+1)
               else if(swbotb.eq.7.or. swbotb .eq. -2)then ! free drainage option
                  soil%kmean(mesh%numnod+1) = hconduc(soil%h(mesh%numnod),soil%theta(mesh%numnod),heat%rfcp(mesh%numnod),heat%tsoil(mesh%numnod), &
                                               soil%vg_params(mesh%numnod), &
                                               soil%iHWCKmodel(soil%layer(mesh%numnod)), &
                                               soil%fluseksatexm(mesh%numnod), &
                                               mesh%numnod, soil)        ! [SS-GR-UTILS Task 6]
                  soil%kmean(mesh%numnod+1) = soil%kmean(mesh%numnod+1)          ! [SS-SWC S-1.4b]
                  soil%qbot = -1.0d0 * soil%kmean(mesh%numnod+1)
                  F(NN) = F(NN) - soil%qbot
               ! Lysimeter option
               else if(swbotb.eq.8)then
                  if (flboth) then
                     soil%hbot = hplate
                     F(NN) = F(NN) + soil%kmean(NN+1) * hgrad(NN+1)
                  else
                     soil%qbot = 0.0d0
                  end if
               ! Flux bottom boundary
               else
                  F(NN) = F(NN) - soil%qbot
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
         do i = 1,mesh%numnod
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
               if(abs( soil%h(i)-hold(i) ) .gt. state%cfg%simulation%numerical%critdevh2cp)then
                  flnonconv2(i) = .true. ; flnonconv   = .true.
               endif
            else
               if(abs( soil%h(i)-hold(i) )/abs(hold(i)) .gt.                 &
     &                state%cfg%simulation%numerical%critdevh1cp)then
                  flnonconv2(i) = .true. ; flnonconv   = .true.
               endif
            end if
            flnonconv2(1) = flnonconv2(1) ! for Forcheck

         enddo

         ! Test for waterbalance of ponding layer
         if (soil%ftoph) then
            soil%qtop = -soil%kmean(1)*((soil%hsurf - soil%h(1))/mesh%disnod(1)+1.0d0)
            if(.not.flnonconv) then
               deviat = soil%pond - soil%pondm1 + soil%reva*time%dt - (atmo%nraidt+atmo%nird+atmo%melt)*time%dt &  ! [SS-SWC S-2.12B] [TC-8]
     &                - soil%runon*time%dt  +  soil%runots  - soil%qtop * time%dt  ! [TC-8]
               if( abs(deviat) .gt. state%cfg%simulation%numerical%critdevponddt) then
                  flnonconv3 = .true. ; flnonconv   = .true.
                  flnonconv3 = flnonconv3 ! for Forcheck
               end if
            end if
         end if

         if(dabs(sum1).gt.CritDevBalTot) flnonconv   = .true.

         ! Save sum for next iteration
         sumold = sum

         if(.not.flnonconv )then !  convergence has been reached
            if(swbotb.eq.1 .and. (.not.soil%fllowgwl))then
               ! Derive vertical flux profile in order to find qbot as a
               ! lower boundary condition for the saturated part of the soil
               ! system
               qv(1) = soil%qtop
               do i=NN+1,mesh%numnod
                  soil%theta(i) = soil%vg_params(i)%thetas ! [SS-GR-UTILS Task 15]
                  soil%theta(i) = soil%theta(i)                        ! [SS-SWC S-1.4a]
               end do
               do i=1,mesh%numnod
                 qv(i+1) = qv(i) +mesh%dz(i)*soil%FrArMtrx(i)*(soil%theta(i)-soil%thetm1(i))&  ! [SS-SWC S-2.3]
     &                          / time%dt + sink(i) - source(i) + soil%qrot(i)  ! [TC-8]
               end do
               soil%qbot = qv(mesh%numnod+1)

               do i=NN+1,mesh%numnod
                  soil%h(i) = soil%h(i-1) + mesh%disnod(i)*(qv(i)/soil%kmean(i)+1.0d0)
                  soil%h(i) = soil%h(i)                                ! [SS-SWC S-1.4a]
               end do

            end if

            ! Calculate new groundwater level
            call calcgwl (state)

            if(swbotb.ne.1.and.abs(soil%gwl-soil%gwlm1).ge.state%cfg%simulation%numerical%gwlconv .AND.          &  ! [SS-SWC S-2.12B]
     &         abs(soil%gwl-999d0).gt.1.d0.and.abs(soil%gwlm1-999d0).gt.1.d0) then  ! [SS-SWC S-2.12B]
               call dtdpst ('year-month-day',time%t1900+1.001d0,datetmp)  ! [TC-8]
                 messag = ' Change of groundwater level exceeds'//      &
     &           ' criterion at '//datetmp//'. Consider reduction of dtMin'
               call log_warn('Headcalc', messag)
            endif

            ! Recording of number of iteration steps needed
            itnumb(min(100,numbit),1)=itnumb(min(100,numbit),1)+1 
            itnumb(min(100,numbit),2)=itnumb(min(100,numbit),2)+iBackTr 

            return
            
         end if
      End Do

      ! Convergence could not been reached
      if (.not.time%fldtmin ) then  ! [TC-8]
         ! Reset soil state variables
         do j = 1,mesh%numnod
            soil%h(j) = soil%hm1(j)                                    ! [SS-SWC S-2.3]
            soil%h(j) = soil%h(j)                                      ! [SS-SWC S-1.4a]
            soil%theta(j) = soil%thetm1(j)                             ! [SS-SWC S-2.3]
            soil%theta(j) = soil%theta(j)                              ! [SS-SWC S-1.4a]
         enddo
         soil%kmean(mesh%numnod+1) = soil%k(mesh%numnod)
         soil%kmean(mesh%numnod+1) = soil%k(mesh%numnod)                          ! [SS-SWC S-1.4b]
         soil%gwl  = soil%gwlm1                         ! [SS-SWC S-2.12B]
         soil%pond = soil%pondm1                        ! [SS-SWC S-2.12B]

         ! Reset and continue iteration with smaller timestep!
         fldecdt = .true.

         return

      else
         ! Write warning to screen and log file
         if (flwarn_hc .and. iwarn_hc.lt.5) then
            iwarn_hc = iwarn_hc + 1
            call dtdpst                                                 &
     &        ('year-month-day,hour:minute:seconds',time%t1900,datetime)  ! [TC-8]
            messag = ' No convergence was reached of Richards'//        &
     &        ' equation at '//datetime//                               &
     &        ' no more than 4 warnings per date - SWAP did continue !'
            call log_warn('Headcalc', messag)
            if (iwarn_hc.gt.4) then
              flwarn_hc = .false.  
            endif
         endif

         if(fldumpconvcrit) then
           call log_debug('Headcalc', 'Datetime = ' // datetime)
           call log_debug('Headcalc', 't1900    = ' // to_str(time%t1900))
           call log_debug('Headcalc', 'time%dtmin = ' // to_str(time%dtmin))
           call log_debug('Headcalc', 'dt    = ' // to_str(time%dt))
           call log_debug('Headcalc', 'ftoph  = ' // to_str(soil%ftoph))
           call log_debug('Headcalc', 'CritDevBalCp  = ' // to_str(CritDevBalCp))
           call log_debug('Headcalc', 'CritDevBalTot = ' // to_str(CritDevBalTot))
           call log_debug('Headcalc', 'CritDz        = ' // to_str(CritDz))
           call log_debug('Headcalc', 'critdevh1cp   = ' // to_str(state%cfg%simulation%numerical%critdevh1cp))
           call log_debug('Headcalc', 'critdevh2cp   = ' // to_str(state%cfg%simulation%numerical%critdevh2cp))
           call log_debug('Headcalc', 'flnonconv  = ' // to_str(flnonconv))
           call log_debug('Headcalc', 'flnonconv3 = ' // to_str(flnonconv3))
           call log_debug('Headcalc', 'flunsatok(1) = ' // to_str(flunsatok(1)))
           call log_debug('Headcalc', 'flunsatok(2) = ' // to_str(flunsatok(2)))
           call log_debug('Headcalc', 'flunsatok(3) = ' // to_str(flunsatok(3)))
           call log_debug('Headcalc', 'pondm1 = ' // to_str(soil%pondm1))
           call log_debug('Headcalc', 'pond   = ' // to_str(soil%pond))
           call log_debug('Headcalc', 'gwlm1 = ' // to_str(soil%gwlm1))
           call log_debug('Headcalc', 'gwl   = ' // to_str(soil%gwl))
           call log_debug('Headcalc', 'node,flnonconv1_F, F, flnonconv2_h, hm1, h, thetm1, theta')
           do j = 1, mesh%numnod
              call log_debug('Headcalc', &
                  to_str(j) // ',' // to_str(flnonconv1(j)) // ',' // to_str(F(j)) // &
                  ',' // to_str(flnonconv2(j)) // ',' // to_str(soil%h(j)) // &
                  ',' // to_str(soil%hm1(j)) // ',' // to_str(soil%theta(j)) // &
                  ',' // to_str(soil%thetm1(j)))
           enddo
         endif

         ! Continue without convergence !!!
         return

      endif

      end associate  ! soil%theta/.../soil%gwlm1 => soil [SS-SWC S-1.4a/b/S-2.3]; mesh%numnod/mesh%dz/mesh%z/mesh%disnod [GR-BH C4]; drai%nrlevs/swbotb [GR-BH Audit 31]

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
     use doln  ! provides do_ln_trans (parameter)
      ! [SS-GR-FINAL B9] blanket use Variables → explicit only-list; all symbols DEFERRED
      ! macp/mabbc/matabentries → swap_array_dimensions (dimension constants)
      use swap_array_dimensions, only: macp, mabbc, matabentries
      use variables, only: &
         ! swsophy retired (→soil%swsophy)

         ! DEFERRED: swhyst — hysteresis switch; Phase C3
         swhyst, &
         ! [GR-SOL 2026-05-24] swinco retired — read via soil%swinco
         ! DEFERRED: swkmean — mean K averaging method; Phase C3
         ! DEFERRED: paramvg(21,maho) — VanGenuchten parameters table; Phase C3
         paramvg, &
         ! DEFERRED: relsatthr/ksatthr(maho) — threshold saturations; Phase C3
         relsatthr, ksatthr, &
         ! DEFERRED: iHWCKmodel(maho) — hydraulic conductivity model switch per mesh%layer; Phase C3
         iHWCKmodel, &
         ! DEFERRED: numtab/numtablay — number of table entries per node/mesh%layer; Phase C3
         numtab, numtablay, &
         ! DEFERRED: ientrytab/ientrytablay — table entry indices; Phase C3
         ientrytab, ientrytablay, &
         ! DEFERRED: sptab/sptablay — soil property tables; Phase C3
         sptab, sptablay, &
         ! [GR-SOIL 2026-05-24] nod1lay/numlay retired — read via state%mesh%{nod1lay,numlay}
         ! DEFERRED: gwli — initial groundwater level; Phase C3
         gwli, &
         ! DEFERRED: zi/nhead — initial head table entries; Phase C3
         zi, nhead, &
         ! DEFERRED: h_enpr — air entry pressure; Phase C3
         h_enpr, &
         ! DEFERRED: cQMpLatSs — macropore lateral flux (retired-zero); Phase D
         cQMpLatSs
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

      ! [GR-BH C4] mesh globals aliased via mesh for all cases
      ! [GR-BH Audit 31] swbotb aliased via soil%swbotb_runtime
      associate (mesh => state%mesh,         &
                 soil => state%soilwater,    &
                 drai => state%drainage,     &
                 heat => state%heat,         &
                 atmo => state%atmosphere,   &
                 time => state%timecontrol,  &
                 swbotb => state%soilwater%swbotb_runtime)

      select case (task)
      case (1)

         ! Initialize Soilwater rate/state variables


         ! Initialize miscellaneous
         ! [SS-SWC S-2.12B] legacy half-writes dropped — soil is canonical
         soil%hatm = -2.75e5_real64                           ! [SS-SWC S-1.3/S-2.12B]
      atmo%nraidt = 0.0_real64
      atmo%nird   = 0.0_real64
      if (soil%swinco.ne.3) then
        atmo%ldwet = 0.0_real64
        atmo%spev  = 0.0_real64
        atmo%saev  = 0.0_real64
      endif
      soil%runon = 0.0_real64                                 ! [SS-SWC S-1.3/S-2.12B]
      soil%qtop = 0.d0
      do i = 1,mesh%numnod+1
        soil%q(i) = 0.0_real64                               ! [SS-SWC S-1.3/S-2.12B]
      enddo
      soil%evp = 0.0_real64                                  ! [SS-SWC S-1.3/S-2.12B] state sized mesh%numnod; blanket zero
      ! [SS-HEAT] Task 9: legacy rfcp global retired; heat%rfcp is authoritative
      if (allocated(heat%rfcp)) heat%rfcp = 1.0d0
      state%surfacewater%vtair = 0.0d0
      cQMpLatSs = 0.0d0

      ! Soil physics: tabulated or MualemVanGenuchten functions
      ! [SS-GR-UTILS Task 4] Mirror mesh%layer + iHWCKmodel into state (both branches).
      ! mesh%layer(:) set by CalcGrid; iHWCKmodel(:) set by config_to_variables.
      ! State fields allocated by soilwater_init which ran before SoilHydraulics(1).
      do node = 1, mesh%numnod
         soil%layer(node) = mesh%layer(node)
      end do
      do lay = 1, mesh%numlay
         soil%iHWCKmodel(lay) = iHWCKmodel(lay)
      end do
      ! BiModal/NoVap: only set via legacy readswap (not TOML path); stay .false.

      if(soil%swsophy.eq.1) then
         ! Tabulated functions (h,theta,k,dthetadh,dkdtheta) tabulated
         do node = 1,mesh%numnod
          numtab(node) = numtablay(mesh%layer(node))
          do i=0,matabentries
             ientrytab(node,i) = ientrytablay(mesh%layer(node),i)
          end do
        end do
        do node = 1,mesh%numnod
          do i = 1,7
            do j = 1, numtab(node)
              sptab(i,node,j) = sptablay(i,mesh%layer(node),j)
            end do
          end do
          ! Populate vg_params directly — [SS-GR-UTILS Task 15] cofgen removed
          soil%vg_params(node)%thetar = 0.0_real64
          soil%vg_params(node)%thetas = sptab(2, node, numtab(node))
          soil%vg_params(node)%ksat   = sptab(3, node, numtab(node))
          if (do_ln_trans) soil%vg_params(node)%ksat = dexp(soil%vg_params(node)%ksat)
          ! [SS-GR-UTILS Task 4] Mirror numtab/ientrytab/sptab into state (swsophy=1 only).
          soil%numtab(node) = numtab(node)
          do i = 0, matabentries
             soil%ientrytab(node,i) = ientrytab(node,i)
          end do
          do i = 1, 7
             do j = 1, numtab(node)
                soil%sptab(i,node,j) = sptab(i,node,j)
             end do
          end do
        end do
        do lay = 1,mesh%numlay
          soil%ksatfit(lay) = soil%vg_params(mesh%nod1lay(lay))%ksat   ! [SS-SWC S-2.3] [GR-BH Task 36] ksatfit global retired
          soil%thetsl(lay) = soil%vg_params(mesh%nod1lay(lay))%thetas  ! [SS-SWC S-1.3/S-2.12B]
        end do
      else
         ! MvanG functions
         do node = 1,mesh%numnod
          lay = mesh%layer(node)
          ! Populate vg_params directly — [SS-GR-UTILS Task 15] cofgen removed
          soil%vg_params(node)%thetar           = paramvg(1, lay)
          soil%vg_params(node)%thetas           = paramvg(2, lay)
          soil%vg_params(node)%ksat             = paramvg(3, lay)
          soil%vg_params(node)%alpha            = paramvg(4, lay)
          soil%vg_params(node)%lpar             = paramvg(5, lay)
          soil%vg_params(node)%npar             = paramvg(6, lay)
          soil%vg_params(node)%mpar             = paramvg(7, lay)
          soil%vg_params(node)%alphaw_sentinel  = -9999.9_real64
          soil%vg_params(node)%h_enpr           = paramvg(9, lay)
          soil%vg_params(node)%ksatexm          = paramvg(10, lay)
          if (soil%vg_params(node)%ksatexm > 0.0_real64) soil%fluseksatexm(node) = .true.
          soil%vg_params(node)%relsatthr        = relsatthr(lay)
          soil%vg_params(node)%ksatthr          = ksatthr(lay)
          if (iHWCKmodel(lay) ==  3 .OR. iHWCKmodel(lay) ==  6 .OR. iHWCKmodel(lay) ==  7 .OR. &
              iHWCKmodel(lay) == 10 .OR. iHWCKmodel(lay) == 11) then
             soil%vg_params(node)%alpha_2  = paramvg(13, lay)
             soil%vg_params(node)%npar_2   = paramvg(14, lay)
             soil%vg_params(node)%mpar_2   = paramvg(15, lay)
             soil%vg_params(node)%omega_1  = paramvg(16, lay)
             soil%vg_params(node)%omega_2  = paramvg(17, lay)
          end if
          if (iHWCKmodel(lay) ==  5 .OR. iHWCKmodel(lay) ==  7) then
             soil%vg_params(node)%h0 = paramvg(18, lay)
          end if
          if (iHWCKmodel(lay) ==  8 .OR. iHWCKmodel(lay) ==  9 .OR. &
              iHWCKmodel(lay) == 10 .OR. iHWCKmodel(lay) == 11) then
             soil%vg_params(node)%h0      = paramvg(18, lay)
             soil%vg_params(node)%ha      = paramvg(19, lay)
             soil%vg_params(node)%apar    = paramvg(20, lay)
             soil%vg_params(node)%omega_k = paramvg(21, lay)
          end if
        end do
        soil%thetsl = 0.0_real64                            ! [SS-SWC S-1.3/S-2.12B]
        do lay = 1, mesh%numlay
          soil%thetsl(lay) = paramvg(2,lay)                 ! [SS-SWC S-1.3/S-2.12B]
        end do
      endif

! --- saturated and residual watercontent of each node; hysteresis parameters
      do node = 1,mesh%numnod
        lay = mesh%layer(node)
        soil%thetar(node) = soil%vg_params(node)%thetar        ! [SS-SWC S-1.3/S-2.12B]
        soil%thetas(node) = soil%vg_params(node)%thetas        ! [SS-SWC S-1.3/S-2.12B]
        !!! Kroes: disable combi of swsophy=1 and swhyst=1
        if (swhyst.eq.1) then
           ! Wetting curve
           soil%indeks(node) = 1                            ! [SS-SWC S-1.3/S-2.12B]
           soil%vg_params(node)%alpha = paramvg(8,lay)      ! [SS-GR-UTILS Task 15] hysteresis: wetting alpha
        elseif (swhyst.eq.0.or.swhyst.eq.2) then
           ! Drying branch or simulation without hysteresis
           soil%indeks(node) = -1                           ! [SS-SWC S-1.3/S-2.12B]
           soil%vg_params(node)%alpha = paramvg(4,lay)      ! [SS-GR-UTILS Task 15] hysteresis: drying alpha
        endif
      end do

      if (soil%swinco.eq.1) then
         ! Pressure head profile is input
         ! [SS-SWC S-2.12B] legacy h(:) reads/writes retargeted to soil%h
         do i = 1, nhead
          tab(i*2) = soil%h(i)
          tab(i*2-1) = abs(zi(i))
        end do
        do i = 1, mesh%numnod
          soil%h(i) = afgen(tab,macp*2,abs(mesh%z(i)))            ! [SS-SWC S-1.3/S-2.12B]
        end do
      endif
      if (soil%swinco.eq.2 .and. swbotb.ne.8) then
        if (abs(gwli-(mesh%z(mesh%numnod)-0.5d0*mesh%dz(mesh%numnod))) .lt.1.0d-4) then
          messag = 'Initial groundwaterlevel (SWINCO=2) is '//          &
     &    'too close to bottom of soil profile'//                       &
     &    ' must be corrected!'
          call fatalerr_collected ('soilwater',messag)
        endif
      endif
      if (soil%swinco.eq.3) then
        if (nhead.ne.mesh%numnod) then
          messag = 'Initial data are read from file (SWINCO=3) and '//  &
     &    'number of nodes/compartments is not consistent with NUMNOD'//&
     &    'must be corrected!'
          call fatalerr_collected ('soilwater',messag)
        endif
      endif
      if (soil%swinco.eq.1.or.soil%swinco.eq.3) then
         ! Determine groundwater level — [SS-SWC S-2.12B] all legacy half-writes dropped
         if (soil%h(mesh%numnod) .gt. -1.d-5) then
          i = mesh%numnod
          do while ((soil%h(i) .gt. -1.d-5) .and. (i .gt. 1))
              i = i - 1
          end do
          if (soil%h(i) .lt. -1.d-5) then
            soil%gwl = mesh%z(i+1) + soil%h(i+1) / (soil%h(i+1) - soil%h(i)) * (mesh%z(i) - mesh%z(i+1))  ! [SS-SWC S-1.3/S-2.12B]
            ! Assume hydrostatic equilibrium in saturated part
            do j = i+1, mesh%numnod
              soil%h(j) = soil%gwl - mesh%z(j)                     ! [SS-SWC S-1.3/S-2.12B]
            end do
          endif
        endif
      else
         ! Pressure head profile is calculated from groundwater level
         if (swbotb.eq.1) then
          soil%gwl = afgen (soil%gwltab,mabbc*2,time%t1900+time%dt-1.d0)   ! [SS-SWC S-1.3/S-2.12B] [TC-8]

          if(abs(soil%gwl-(mesh%z(mesh%numnod)-0.5d0*mesh%dz(mesh%numnod))) .lt.1.0d-4) then
          messag = 'Groundwaterlevel as bottom boundary (SWBOTB=1) is'//&
     &    'below or to close to bottom of soil profile'//               &
     &    ' must be corrected!'
            call fatalerr_collected ('soilwater',messag)
          endif
        else
          soil%gwl = gwli                                   ! [SS-SWC S-1.3/S-2.12B]
        endif
        if (soil%gwl.gt.0.0d0) then
          soil%pond = soil%gwl                                ! [SS-SWC S-1.3/S-2.12B]
        else
          soil%pond = 0.0_real64                            ! [SS-SWC S-1.3/S-2.12B]
        endif
        do i = 1,mesh%numnod
          soil%h(i) = soil%gwl - mesh%z(i)                         ! [SS-SWC S-1.3/S-2.12B]
        end do
      endif

      ! In case of preferential flow, adjust Van Genuchten parameters
      do i = 1, mesh%numnod
        soil%theta(i) = watcon(soil%h(i), &
                              soil%vg_params(i), &
                              soil%iHWCKmodel(soil%layer(i)), &
                              i, soil)                  ! [SS-SWC S-1.3/S-2.12B] [SS-GR-UTILS Task 5]
      end do

      ! Hydraulic conductivities, differential moisture capacities
      ! and mean hydraulic conductivities for each node
      do node = 1,mesh%numnod
        soil%dimoca(node) = moiscap(soil%h(node), &
                                  soil%vg_params(node), &
                                  soil%iHWCKmodel(soil%layer(node)), &
                                  time%dt, &
                                  node, soil)             ! [SS-SWC S-1.3/S-2.12B] [SS-GR-UTILS Task 7]

        soil%FrArMtrx(node) = 1.0_real64                    ! [SS-SWC S-1.3/S-2.12B]
        soil%k(node) = hconduc(soil%h(node),soil%theta(node),heat%rfcp(node),heat%tsoil(node), &
                             soil%vg_params(node), &
                             soil%iHWCKmodel(soil%layer(node)), &
                             soil%fluseksatexm(node), &
                             node, soil)                            ! [SS-SWC S-1.3/S-2.12B] [SS-GR-UTILS Task 6]

        if(node.gt.1) soil%kmean(node) = hcomean(state%cfg%simulation%numerical%swkmean,soil%k(node-1),soil%k(node),mesh%dz(node-1),mesh%dz(node))  ! [SS-SWC S-1.3/S-2.12B]
      end do
      soil%kmean(mesh%numnod+1) = soil%k(mesh%numnod)                   ! [SS-SWC S-1.3/S-2.12B]

      ! Initial soil water storage
      do i = 1, mesh%numnod
         soil%FrArMtrx(i) = 1.0_real64                      ! [SS-SWC S-1.3/S-2.12B]
      enddo
      soil%volact = 0.0_real64                              ! [SS-SWC S-1.3/S-2.12B]
      call watstor (state)                                ! [SS-SWC S-1.6] state arg added
      soil%volini = soil%volact                               ! [SS-SWC S-1.3/S-2.12B]
      soil%pondini = soil%pond                                ! [SS-SWC S-1.3/S-2.12B]

      ! Initial groundwater level
      call calcgwl (state)

      call log_info('soilwater', 'Soil state initialized: gwl=' // to_str(real(soil%gwl,4)) // &  ! [SS-SWC S-2.3]
                    ' cm, numnod=' // to_str(mesh%numnod) // ', numlay=' // to_str(mesh%numlay) // &
                    ', volini=' // to_str(real(soil%volini,4)) // ' cm')              ! [SS-SWC S-2.3]


      return

      case (2)

         ! Calculate Soilwater rate/state variables

         ! Reset intermediate soil water fluxes — [SS-SWC S-2.12B] reset() handles all
         if (time%flDayStart) then  ! [TC-8]
          call soil%reset_intermediate_per_day()              ! [SS-SWC S-2.1]
      end if

      if (time%flZeroIntr) then
        call soil%reset_intermediate()                        ! [SS-SWC S-2.1]

        soil%IPondBeg = soil%pond     ! [SS-SWC S-1.4b/S-2.3/S-2.12B]
        do node = 1, mesh%numnod
          soil%IThetaBeg(node) = soil%theta(node)  ! [SS-SWC S-1.4b/S-2.3/S-2.12B]
        enddo

        ! [MACRO-RETIRE 2026-05-12] macropore(5,...) retired (ADR 0040).
      endif

      ! Reset cumulative soil water fluxes — [SS-SWC S-2.12B] reset() handles all
      if (time%flZeroCumu) then
        call soil%reset_cumulative()                         ! [SS-SWC S-2.1]

        ! [MACRO-RETIRE 2026-05-12] macropore(6,...) retired (ADR 0040).

        ! Reset initial water storage and ponding
        soil%volini = soil%volact          ! [SS-SWC S-2.12B]
        soil%pondini = soil%pond           ! [SS-SWC S-2.12B]
      endif

      ! Save state variables of time = t
      call SoilWaterStateVar(1, state)

      ! Calculate new soil water state variables
      call headcalc(state)

      return

      case (3)

         ! Update hydraulic conductivities to time level t+1
         ! [SS-SWC S-2.12B] all legacy half-writes dropped
         do i = 1,mesh%numnod
         soil%k(i) = hconduc(soil%h(i),soil%theta(i), &
                                        heat%rfcp(i),heat%tsoil(i), &
                                        soil%vg_params(i), &
                                        soil%iHWCKmodel(soil%layer(i)), &
                                        soil%fluseksatexm(i), &
                                        i, soil)                    ! [SS-SWC S-1.4b/S-2.12B] [SS-GR-UTILS Task 6]
         if(i.gt.1)then
            soil%kmean(i) = hcomean(state%cfg%simulation%numerical%swkmean,soil%k(i-1),soil%k(i),mesh%dz(i-1),mesh%dz(i))  ! [SS-SWC S-1.4b/S-2.12B]
         end if
      enddo
      soil%kmean(mesh%numnod+1) = soil%k(mesh%numnod)  ! [SS-SWC S-1.4b/S-2.12B]

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

      end associate  ! mesh%numnod/mesh%dz/mesh%z/mesh%layer [GR-BH C4]; swbotb [GR-BH Audit 31]

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
      use swap_state_mod, only: swap_state_t
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none

      integer,            intent(in)    :: task
      type(swap_state_t), intent(inout) :: state

      integer :: i

      associate (mesh => state%mesh, soil => state%soilwater)

         select case (task)
         case (1)
            ! Save state variables at time = t.
            do i = 1, mesh%numnod
               soil%hm1(i)    = soil%h(i)
               soil%thetm1(i) = soil%theta(i)
            end do
            soil%gwlm1  = soil%gwl
            soil%pondm1 = soil%pond

         case (2)
            ! Reset soil state variables.
            do i = 1, mesh%numnod
               soil%h(i)     = soil%hm1(i)
               soil%theta(i) = soil%thetm1(i)
            end do
            soil%kmean(mesh%numnod + 1) = soil%k(mesh%numnod)
            soil%gwl  = soil%gwlm1
            soil%pond = soil%pondm1

         case default
            call fatalerr_collected('SoilWaterStateVar', 'Illegal value for TASK')
         end select

      end associate

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
   subroutine hysteresis(state)
      ! Stragglers still bare-global: tau (scalar) + paramvg (per-layer VG table).
      use variables,             only: tau, paramvg
      use soilhydraulics_utils,  only: moiscap, prhead
      use swap_array_dimensions, only: macp
      use swap_state_mod,        only: swap_state_t
      use hydraulic_params_mod,  only: vanGenuchten_params_t
      implicit none

      type(swap_state_t), intent(inout) :: state

      integer :: node, lay, indtem(macp)
      real(8) :: delp, sew, sed, fvalue
      real(8) :: thetar(macp), thetas(macp), alfamg(macp)
      type(vanGenuchten_params_t) :: vg_hys

      associate (mesh => state%mesh, soil => state%soilwater, time => state%timecontrol)

         ! Check for hysteretic reversal.
         do node = 1, mesh%numnod
            delp = soil%hm1(node) - soil%h(node)
            if (delp/float(soil%indeks(node)) .gt. tau .and.          &
         &      soil%h(node) .lt. -10.0d0 .and. soil%h(node) .gt. -1.0d3) then
               indtem(node) = -soil%indeks(node)
            else
               indtem(node) =  soil%indeks(node)
            end if
         end do

         ! Adapt parameters along scanning curves.
         do 100 node = 1, mesh%numnod
            lay = mesh%layer(node)

            ! No change.
            if (indtem(node) .eq. soil%indeks(node) .or.              &
         &       abs(paramvg(4, lay) - paramvg(8, lay)) .lt. 1.d-4) goto 100

            ! Relative saturation.
            sew = (1.0d0 + (paramvg(8, lay)*(-soil%h(node)))**paramvg(6, lay))**(-paramvg(7, lay))
            sed = (1.0d0 + (paramvg(4, lay)*(-soil%h(node)))**paramvg(6, lay))**(-paramvg(7, lay))

            ! Flip the scanning index.
            soil%indeks(node) = -1*soil%indeks(node)

            if (soil%indeks(node) .eq. 1) then
               ! Wetting branch.
               alfamg(node) = paramvg(8, lay)
               thetas(node) = paramvg(2, lay)
               thetar(node) = (soil%theta(node) - thetas(node)*sew)/(1.0d0 - sew)

               fvalue = thetar(node)
               if (thetar(node) .lt. paramvg(1, lay)) thetar(node) = paramvg(1, lay)
               if (thetar(node) .gt. paramvg(2, lay)) thetar(node) = paramvg(2, lay)
               soil%vg_params(node)%alpha  = alfamg(node)
               soil%vg_params(node)%thetar = thetar(node)
               soil%thetar(node)           = thetar(node)
               soil%vg_params(node)%thetas = thetas(node)
               soil%thetas(node)           = thetas(node)
               if (abs(fvalue - thetar(node)) .gt. 1.d-10) then
                  vg_hys         = soil%vg_params(node)
                  vg_hys%thetar  = thetar(node)
                  vg_hys%thetas  = thetas(node)
                  vg_hys%alpha   = alfamg(node)
                  soil%h(node) = prhead(mesh%disnod(node), soil%theta(node), soil%h, &
                                        soil%iHWCKmodel(soil%layer(node)),           &
                                        node, soil, vg_in=vg_hys)
               end if
            else
               ! Drying branch.
               alfamg(node) = paramvg(4, lay)
               thetar(node) = paramvg(1, lay)
               thetas(node) = thetar(node) + (soil%theta(node) - thetar(node))/sed

               fvalue = thetas(node)
               if (thetas(node) .lt. paramvg(1, lay)) thetas(node) = paramvg(1, lay)
               if (thetas(node) .gt. paramvg(2, lay)) thetas(node) = paramvg(2, lay)
               soil%vg_params(node)%alpha  = alfamg(node)
               soil%vg_params(node)%thetar = thetar(node)
               soil%thetar(node)           = thetar(node)
               soil%vg_params(node)%thetas = thetas(node)
               soil%thetas(node)           = thetas(node)
               if (abs(fvalue - thetas(node)) .gt. 1.d-10) then
                  vg_hys         = soil%vg_params(node)
                  vg_hys%thetar  = thetar(node)
                  vg_hys%thetas  = thetas(node)
                  vg_hys%alpha   = alfamg(node)
                  soil%h(node) = prhead(mesh%disnod(node), soil%theta(node), soil%h, &
                                        soil%iHWCKmodel(soil%layer(node)),           &
                                        node, soil, vg_in=vg_hys)
               end if
            end if

            ! Update moisture capacity.
            soil%dimoca(node) = moiscap(soil%h(node), soil%vg_params(node),     &
                                        soil%iHWCKmodel(soil%layer(node)),      &
                                        time%dt, node, soil)
 100     continue

      end associate

      return
      end subroutine hysteresis

end module