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
      use boundbottom_mod, only: BoundBottom
      use boundtop_mod, only: boundtop, PONDRUNOFF
      use rootextraction_mod, only: RootExtraction
      use array_utils, only: afgen
      use soilhydraulics_utils, only: watcon, hconduc, moiscap, hcomean, dhconduc
      use swap_constants, only: nihil
      use macropore_mod, only: macropore
      use swap_state_mod, only: swap_state_t
      use soilhydraulics_utils, only: dkmean
      use soilwaterbalance_mod, only: calcgwl, fluxes
      use numericalsolvers_mod, only: tridag, bandec, banbks
      implicit none

      ! Arguments
      type(swap_state_t), intent(in) :: state

      ! Local variables
      integer   i,j, itry,  MaxIt1, ndr, NN, iBackTr
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
      data    ndr            / 5 / 

      real(8) hgrad(macp+1), dkdh(macp)

      ! Function and solver variables
      integer indx(macp), ierror
      real(8) a(macp,3), a1(macp,1), b(macp), d, q1
      logical flok

      ! Note: flwarn_hc, iwarn_hc, nstep_hc moved to variables.f90 module
      ! (previously local SAVE variables - now global for multi-instance support)

      if (fldaystart) then
         flwarn_hc = .true.
         iwarn_hc = 0
      endif
      call dtdpst                                                       &
     &        ('year-month-day,hour:minute:seconds',t1900,datetime)

      ! Summation of sink terms (constant for the current time step)
      iBackTr   = 0
      flunsatok(1) = .false.
      flunsatok(2) = .false.
      flunsatok(3) = .false.
      do i=1,numnod
         sink(i) = evp(i)
         do j=1,ndr
            sink(i) = sink(i) + qdra(j,i)
         end do
      end do
      source(1:numnod) = qssdi(1:numnod)

      ! Set value of macropore area at soil surface
      ArMpSs = 0.d0
      if (FlMacropore .and. Z_Tp.gt.-1.d-8) ArMpSs = ArMpTp

      ! Groundwater level specified
      if(swbotb.eq.1)then
         fllowgwl = .false.
         if(gwlinp.ge.z(1)-1.0d-4)then

            q0 = (nraidt+nird+melt)*(1.0d0-ArMpSs) + runon - reva 
            call pondrunoff (state)
            q1 = - q0 + (pond - pondm1)/dt + runots / dt
            theta(1) = watcon(1,gwlinp)
            kmean(1) = hconduc(1,gwlinp,theta(1),rfcp(1))
            ! In case of static macropores FrArMtrx < 1
            if(FlMacropore) kmean(1) = FrArMtrx(1) * kmean(1)

            qv(1) = q1
            do i=1,numnod
               qv(i+1) = qv(i) +dz(i)*FrArMtrx(i)*(theta(i)-thetm1(i))  &
     &                          / dt+ sink(i) - source(i) + qrot(i) 
            end do
            qbot = qv(numnod+1)
            h(1) = gwlinp + disnod(1)*(qv(1)/kmean(1)+1.0d0)
            do i=2,numnod
               h(i) = h(i-1) + disnod(i)*(qv(i)/kmean(i)+1.0d0)
            end do

            if(SwKimpl.eq.1)then
               do i=1,numnod
                  k(i) = hconduc(i,h(i),theta(i),rfcp(i))
                  if(FlMacropore)  k(i) = FrArMtrx(i) * k(i)
                  if(i.gt.1)then
                     kmean(i)=hcomean(swkmean,k(i-1),k(i),dz(i-1),dz(i))
                  end if
               end do
               kmean(numnod+1) = k(numnod)                  
            end if
            call calcgwl ()
            return
         else
            NN = 0
            do while (z(NN+1).gt.gwlinp .and. NN.lt.numnod)
               NN = NN + 1
            end do
            if (z(NN+1).lt.(gwlinp+nihil)) then
               ! Groundwater within soil profile
               if ((z(NN)-gwlinp) .lt. 1.0d-4 .and. (NN.gt.0)) then
                  ! Difference gwlinp with node too small to calculate gradient properly
                  gwlinp = z(NN)
                  NN = NN-1
               endif
            else
               ! Groundwater below soil profile
               fllowgwl = .true.
               hbot = gwlinp - z(numnod) + 0.5*dz(numnod)
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
         k(i) = hconduc(i,h(i),theta(i),rfcp(i))

         if (swcaprise) then
            ! Prevent capillary rise into the root zone !! special for experts only
            if (i .eq. nodncr) then
               if ((h(i) + z(i)) .lt. (h(i+1) + z(i+1))) then  ! negative potential gradient upwards
                  k(i)      = 1.0D-10
                  flcaprise = .true.
               endif
            endif
            if (flcaprise .and. i .eq. nodncr+1) then
                k(i) = 1.0D-10
            endif
         endif
         if(FlMacropore)  k(i) = FrArMtrx(i) * k(i)
         if(i.gt.1)then
            kmean(i) = hcomean(swkmean,k(i-1),k(i),dz(i-1),dz(i))
         end if
      enddo
      kmean(numnod+1) = k(numnod)

      if(SwKimpl.eq.0)then
         do i=2,numnod
            dFdhU(i)   = - kmean(i)  /disnod(i)
            dFdhL(i-1) = dFdhU(i)
         end do
      end if

      do i=2,NN
         hgrad(i) = (h(i-1)-h(i))/disnod(i) + 1.0d0
      end do

      F(1) = (theta(1)-thetm1(1))*FrArMtrx(1)*dz(1)/dt + sink(1) - source(1) + qrot(1) + kmean(2) * hgrad(2)

      call boundtop

      if (FlMacropore) then
         call MACROPORE(2)
      end if

      if (FlRunoff .or. (FlMacropore .and. Z_Tp.gt.-1.d-8))             &   ! Adaptation for GEM
     &    call pondrunoff (state)

      if(ftoph)then
         hgrad(1) = (hsurf-h(1))/disnod(1) + 1.d0
         F(1)     = F(1) - kmean(1) * hgrad(1)
      else
         F(1) = F(1) + qtop
      end if

      do i=2,NN-1
         F(i) = (theta(i)-thetm1(i))*FrArMtrx(i)*dz(i)/dt + sink(i) - source(i) +     &
     &          qrot(i) - kmean(i) * hgrad(i) + kmean(i+1) * hgrad(i+1)
      end do

      if(swbotb.eq.1 .and. (.not.fllowgwl))then
         hgrad(NN+1) = h(NN)/(z(nn)-gwlinp) + 1.0d0
      else if(swbotb.eq.5 .or. (swbotb.eq.1 .and. fllowgwl))then
         hgrad(NN+1) = (h(NN) - hbot) / disnod(NN+1)  + 1.0d0
      else if(swbotb.eq.8 .and. h(NN).gt. Critdz - disnod(NN+1) + hplate) then
         hgrad(NN+1) = (h(NN) - hplate) / disnod(NN+1)  + 1.0d0
         flboth = .true.
      else
         flboth = .false.
      end if

      if(swbotb.eq.1 .and. (.not.fllowgwl))then
         theta(NN)= watcon(NN,h(NN))
         k(NN)    = hconduc(NN,h(NN),theta(NN),rfcp(NN))
         ! In case of static macropores FrArMtrx < 1
         if(FlMacropore) k(NN) = FrArMtrx(NN) * k(NN)
         kmean(NN+1) = hcomean(swkmean,k(NN),cofgen(3,(NN+1)),          &
     &                        dz(NN),dz(NN+1))
         F(NN) = (theta(NN) - thetm1(NN))*FrArMtrx(NN)*dz(NN)/dt +      &
     &           sink(NN)-source(NN)+qrot(NN)-kmean(NN)*hgrad(NN) +kmean(NN+1)*hgrad(NN+1)
      else

         F(NN) = (theta(NN) - thetm1(NN))*FrArMtrx(NN)*dz(NN)/dt        &
     &         - kmean(NN) * hgrad(NN) + sink(NN) - source(NN) + qrot(NN) 

         if(swbotb.eq.3.and.swbotb3Impl.eq.1)then ! Cauchy-relation, implemented as head boundary
            if (SwBotb3ResVert.eq.0) then
               qbot = - (h(NN)+z(NN)-deepgw) / (disnod(NN+1)/kmean(NN+1)+rimlay)
            elseif (SwBotb3ResVert.eq.1) then
               qbot = - (h(NN)+z(NN)-deepgw) / rimlay
            endif
! ---       extra groundwater flux might be added
            if (sw4 .eq. 1) qbot = qbot + afgen(qbotab,mabbc*2,t1900+dt)
            F(NN) = F(NN) - qbot     
         else if(swbotb.eq.5 .or. (swbotb.eq.1 .and. fllowgwl))then ! pressure head at lower boundary specified
            F(NN) = F(NN) + kmean(NN+1) * hgrad(NN+1)
         else if(swbotb.eq.7 .or. swbotb .eq. -2)then ! free drainage option
            kmean(numnod+1) = hconduc(numnod,h(numnod),theta(numnod),rfcp(numnod))
            if(FlMacropore) then 
                kmean(numnod+1) = FrArMtrx(numnod) * kmean(numnod+1)
            endif
            qbot = -1.0d0 * kmean(numnod+1)
            F(NN) = F(NN) - qbot
         ! Lysimeter option
         else if(swbotb.eq.8)then
            if (flboth) then
               hbot = hplate
               F(NN) = F(NN) + kmean(NN+1) * hgrad(NN+1)
            else
               qbot = 0.0d0
            end if
         ! Flux bottom boundary
         else
            F(NN) = F(NN) - qbot
         end if

      end if

      if (FlMacropore) then
         do i= 1, NN
            F(i) = F(i) - QExcMpMtx(i)
         enddo
      endif

      ! Initial estimate of F inner product
      sumold    = 0.0d0
      do i=1,NN
         sumold = sumold + F(i)*F(i)
      end do
      sumold = 0.5d0 * sumold

      ! Start iteration loop, MaxIt specified in the input
      if(fldtmin)then
         MaxIt1 = 2*MaxIt
      else
         MaxIt1 = MaxIt
      end if

      sum= 0.d0   ! For Forcheck
      Do numbit = 1,MaxIt1

         do i = 1, NN
            ! Save values of h
            hold(i)   = h(i)

            ! Derivative of theta to h (differential moisture capacity),
            ! as part of main diagonal
            dimoca(i) = moiscap(i,h(i))

         enddo

         if(SwKimpl.eq.1)then
            do i = 1, NN
               dkdh(i)= dhconduc(i,h(i),theta(i),dimoca(i),rfcp(i))
               if(FlMacropore) dkdh(i) = FrArMtrx(i) * dkdh(i)
            enddo
            do i=2,NN
               dFdhU(i)   = - kmean(i)  /disnod(i)
               dFdhL(i-1) = dFdhU(i)
            end do
         end if

         if (swcaprise .and. flcaprise) then
           dkdh(nodncr+1) = 1.0D-30
         endif

         ! Jacobian matrix elements
         dFdhM(1) = dimoca(1)*FrArMtrx(1)*dz(1)/dt - dFdhL(1)
         ! If the head boundary condition applies: add the k1/(0.5*dz1) term
         ! to the first element of the main diagonal
         if(ftoph) dFdhM(1) = dFdhM(1) + kmean(1)/disnod(1)  

         do i=2,NN-1
            dFdhM(i) = dimoca(i)*FrArMtrx(i)*dz(i)/dt - dFdhU(i)        &
     &                                                - dFdhL(i) 
         end do

         dFdhM(NN) = dimoca(NN)*FrArMtrx(NN)*dz(NN)/dt - dFdhU(NN) 
         if(swbotb.eq.1 .and. (.not.fllowgwl))then
            dFdhM(NN) = dFdhM(NN) + kmean(NN+1)/(z(NN)-gwlinp) 
         else if(swbotb.eq.3.and.swbotb3Impl.eq.1)then ! Cauchy
            if (SwBotb3ResVert.eq.0) then
               dFdhM(NN) = dFdhM(NN) + 1.0d0 /                          &
     &                              (disnod(NN+1)/kmean(NN+1)+rimlay)   
            elseif (SwBotb3ResVert.eq.1) then
               dFdhM(NN) = dFdhM(NN) + 1.0d0 / rimlay
            endif
         else if(swbotb.eq.5 .or. (swbotb.eq.1 .and. fllowgwl))then
            dFdhM(NN) = dFdhM(NN) + kmean(NN+1)/disnod(NN+1)         
         else if(swbotb.eq.7 .or. swbotb .eq. -2)then ! implicitly: kmean(NN+1)
            dFdhM(NN) = dFdhM(NN) + dkdh(NN) * 0.5d0
         else if(swbotb.eq.8 .and. flboth)then
            dFdhM(NN) = dFdhM(NN) + kmean(NN+1)/disnod(NN+1)
         end if

         if(SwKimpl.eq.1)then
            dFdhM(1) = dFdhM(1) + dkdh(1) * hgrad(2) *                  &
     &                 dkmean(swkmean,k(1),k(2),dz(1),dz(2))
            if(ftoph) dFdhM(1) = dFdhM(1) - dkdh(1) * hgrad(1) * 0.5d0
            dFdhL(1) = dFdhL(1) + dkdh(2) * hgrad(2) *                  &
     &                 dkmean(swkmean,k(2),k(1),dz(2),dz(1)) 
            do i=2,NN-1
               dFdhU(i) = dFdhU(i) - dkdh(i-1) * hgrad(i) *             &
     &                    dkmean(swkmean,k(i-1),k(i),dz(i-1),dz(i)) 
               dFdhM(i) = dFdhM(i) - dkdh(i) * hgrad(i) *               &
     &                    dkmean(swkmean,k(i),k(i-1),dz(i),dz(i-1))     &
     &                             + dkdh(i) * hgrad(i+1) *             &
     &                    dkmean(swkmean,k(i),k(i+1),dz(i),dz(i+1))
               dFdhL(i) = dFdhL(i) + dkdh(i+1) * hgrad(i+1) *           &
     &                    dkmean(swkmean,k(i+1),k(i),dz(i+1),dz(i)) 
            end do
            dFdhU(NN) = dFdhU(NN) - dkdh(NN-1) * hgrad(NN) *            &
     &                  dkmean(swkmean,k(NN-1),k(NN),dz(NN-1),dz(NN)) 
            dFdhM(NN) = dFdhM(NN) - dkdh(NN) * hgrad(NN) *              &
     &                  dkmean(swkmean,k(NN),k(NN-1),dz(NN),dz(NN-1))

            if(swbotb.eq.1 .or. swbotb.eq.5 .or. swbotb.eq.8            &
     &         .and. flboth)then
               dFdhM(NN) = dFdhM(NN) + 0.5d0 * dkdh(NN) * hgrad(NN+1)
            end if
         end if

         if (FlMacropore .and. .not.flunsatok(3)) then
            call MACROPORE(3)

            do i= 1, nn
               dFdhM(i) = dFdhM(i) - dFdhMp(i)
            enddo

         endif

         ! Solve the tridiagonal matrix
         call tridag(NN, dFdhU, dFdhM, dFdhL, F, difh, ierror)

         if(ierror.ne.0)then
            call dtdpst ('year-month-day',t1900+1.001d0,datetmp)
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
            if(fldtmin .and. numbit.gt.MaxIt)then              
               factmax = 0.0d0
               do i = 1,NN
                  if(dabs( hold(i) ) .lt. 1.0d0 )then
                     factmax = max( factmax, dabs( difh(i) ) ) 
                  else
                     factmax = max( factmax, dabs( difh(i) / hold(i) ) )
                  end if
               end do
               do i = 1,NN
                  h(i) = hold(i) - difh(i) * min(1.0d0, 1.0d0 / factmax)
               end do
            else
               do i = 1,NN
                  h(i) = hold(i) - factor * difh(i)
               end do
            end if
            do i = 1,NN
              theta(i) = watcon(i,h(i))
            enddo
            do i=2,NN
               hgrad(i) = (h(i-1)-h(i))/disnod(i) + 1.0d0
            end do


            if(SwKimpl.eq.1)then
               call Rootextraction
               do i = 1,NN
                  k(i) = hconduc(i,h(i),theta(i),rfcp(i))
                  if(FlMacropore)  k(i) = FrArMtrx(i) * k(i)
                  if(i.gt.1)then
                     kmean(i)=hcomean(swkmean,k(i-1),k(i),dz(i-1),dz(i))
                  end if
               end do
               kmean(NN+1) = k(NN)
            end if

            ! Prevent capillary rise into the root zone !! special for experts
            if (swcaprise) then
               flcaprise = .false.
               i = nodncr
               if ((h(i) + z(i)) .lt. (h(i+1) + z(i+1))) then  ! negative potential gradient upwards 
                  k(i)      = 1.0D-10
                  flcaprise = .true.
               endif
               if (flcaprise) then
                 k(i+1) = 1.0D-10
               endif
               if(i.gt.1)then
                 kmean(i)=hcomean(swkmean,k(i-1),k(i),dz(i-1),dz(i))
                 kmean(i+1)=hcomean(swkmean,k(i),k(i+1),dz(i),dz(i+1))
               end if
            endif

            ! Calculate F-function
            F(1) = (theta(1) - thetm1(1))*FrArMtrx(1)*dz(1)/dt + sink(1) - source(1) &
     &           + qrot(1) + kmean(2) * hgrad(2)

            if (FlMacropore) QMpLatSsSav = QMpLatSs

            call boundtop

            if (FlMacropore) then
               if (.not.flunsatok(3)) then
                  call MACROPORE(2)
               else
                  QMpLatSs = QMpLatSsSav
               endif
            endif

            if (FlRunoff .or. (FlMacropore .and. Z_Tp.gt.-1.d-8))       &   ! Adaptation for GEM
     &         call pondrunoff (state)

            if(ftoph)then
               hgrad(1) = (hsurf-h(1))/disnod(1) + 1.d0
               F(1) = F(1) - kmean(1) * hgrad(1)
            else
               F(1) = F(1) + qtop
            end if

            do i=2,NN-1
               F(i) = (theta(i)-thetm1(i))*FrArMtrx(i)*dz(i)/dt + sink(i) - source(i) &
     &              + qrot(i) - kmean(i)*hgrad(i)+kmean(i+1)*hgrad(i+1)
            end do

            if(swbotb.eq.1 .and. (.not.fllowgwl))then
               hgrad(NN+1) = h(NN)/(z(nn)-gwlinp) + 1.0d0
           else if(swbotb.eq.5 .or. (swbotb.eq.1 .and. fllowgwl))then
               hgrad(NN+1) = (h(NN) - hbot) / disnod(NN+1)  + 1.0d0
            else if(swbotb.eq.8 .and. flboth)then
               hgrad(NN+1) = (h(NN) - hplate) / disnod(NN+1)  + 1.0d0
            end if

            if(swbotb.eq.1 .and. (.not.fllowgwl))then
               theta(NN) = watcon(NN,h(NN))
               k(NN)     = hconduc(NN,h(NN),theta(NN),rfcp(NN))
               ! In case of static macropores FrArMtrx < 1
               if(FlMacropore)  k(NN) = FrArMtrx(NN) * k(NN)
               kmean(NN+1) = hcomean(swkmean,k(NN),cofgen(3,(NN+1))     &
     &                       ,dz(NN),dz(NN+1))
               F(NN) = (theta(NN) - thetm1(NN))*FrArMtrx(NN)*dz(NN)/dt  &
     &            - kmean(NN) * hgrad(NN) + kmean(NN+1) * hgrad(NN+1)   &
     &            + sink(NN) - source(NN) + qrot(NN)
            else
               F(NN) = (theta(NN) - thetm1(NN))*FrArMtrx(NN)*dz(NN)/dt  &
     &               - kmean(NN) * hgrad(NN)                            &
     &               + sink(NN) - source(NN) + qrot(NN)
               if(swbotb.eq.3.and.swbotb3Impl.eq.1)then ! Cauchy
                  if (SwBotb3ResVert.eq.0) then
                     qbot = - (h(NN)+z(NN)-deepgw) /                    &
     &                                 (disnod(NN+1)/kmean(NN+1)+rimlay)
                  elseif (SwBotb3ResVert.eq.1) then
                     qbot = - (h(NN)+z(NN)-deepgw) / rimlay
                  endif
                  ! Extra groundwater flux might be added
                  if (sw4 .eq. 1) then
                     qbot = qbot + afgen(qbotab,mabbc*2,t1900+dt)
                  end if
                  F(NN) = F(NN) - qbot
               else if(swbotb.eq.5 .or.(swbotb.eq.1 .and. fllowgwl))then
                  ! Pressure head at lower boundary specified
                  F(NN) = F(NN) + kmean(NN+1) * hgrad(NN+1)
               else if(swbotb.eq.7.or. swbotb .eq. -2)then ! free drainage option
                  kmean(numnod+1) = hconduc(numnod,h(numnod),theta(numnod),rfcp(numnod))
                  if(FlMacropore) then
                     kmean(numnod+1) = FrArMtrx(numnod)*kmean(numnod+1)
                  endif
                  qbot = -1.0d0 * kmean(numnod+1)
                  F(NN) = F(NN) - qbot
               ! Lysimeter option
               else if(swbotb.eq.8)then
                  if (flboth) then
                     hbot = hplate
                       F(NN) = F(NN) + kmean(NN+1) * hgrad(NN+1)
                  else
                     qbot = 0.0d0
                  end if
               ! Flux bottom boundary
               else
                  F(NN) = F(NN) - qbot
               end if

            end if

            if (FlMacropore) then
               do i= 1, nn
                  F(i) = F(i) - QExcMpMtx(i)
               enddo
            endif

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
               if(abs( h(i)-hold(i) ) .gt. CritDevh2Cp)then
                  flnonconv2(i) = .true. ; flnonconv   = .true.
               endif
            else
               if(abs( h(i)-hold(i) )/abs(hold(i)) .gt.                 &
     &                CritDevh1Cp)then
                  flnonconv2(i) = .true. ; flnonconv   = .true.
               endif
            end if
            flnonconv2(1) = flnonconv2(1) ! for Forcheck

         enddo

         ! Test for waterbalance of ponding layer
         if (ftoph) then
            qtop = -kmean(1)*((hsurf - h(1))/disnod(1)+1.0d0)
            if(.not.flnonconv .and. (.not.FlMacropore .or. IcTopMp.gt.1)) then
               deviat = pond - pondm1 + reva*dt - (nraidt+nird+Melt)*dt &
     &                - runon*dt  +  runots  - qtop * dt
               if( abs(deviat) .gt. CritDevPondDt) then
                  flnonconv3 = .true. ; flnonconv   = .true.
                  flnonconv3 = flnonconv3 ! for Forcheck
               end if
            end if
         end if

         if(FlMacropore)then
            deviat = pond - pondm1 + reva*dt - (nraidt+nird+Melt)*dt    &
     &             - runon*dt  +  runots  - qtop * dt                   &
     &             + ArMpSs * (nraidt+nird+Melt)*dt + QMpLatSs
            if( abs(deviat) .gt. CritDevPondDt) then
               flnonconv3 = .true. ; flnonconv   = .true.
               flnonconv3 = flnonconv3 ! for Forcheck
            end if
         end if

         ! Implemented to improve iteration performance in case of macropores
         if(dt.lt. 0.01d0 .and. FlMacropore .and. .not.flnonconv3)then 
            flok = .true.
            if (dt.gt.10.d0*dtmin) then
               do i=1,nodgwl
                  if(h(i).gt.0.0d0.and.flnonconv1(i).and.flnonconv2(i)) &
     &            then
                     flok = .false.
                  end if
               end do
            else
               continue
            endif
            if(flok)then
               if(.not.flunsatok(1))then
                  flunsatok(1) = .true.
               else if(.not. flunsatok(2))then
                  flunsatok(2) = .true.
               else
                  flunsatok(3) = .true.                 
               end if
            else
               flunsatok(1) = .false.
               flunsatok(2) = .false.
               flunsatok(3) = .false.
            end if
         end if

         if(dabs(sum1).gt.CritDevBalTot) flnonconv   = .true.

         ! Save sum for next iteration
         sumold = sum

         if(.not.flnonconv )then !  convergence has been reached
        !! Replaced the Nstep with the variable from the main variables.f90 file, the
        !! build still works as before.
            if (FlMacropore) then
               
               FlDecMpRat = .false.
               if (IDecMpRat.gt.0) then
                  if (nstep_hc.lt.10) then
                     nstep_hc = nstep_hc + 1
                  endif
                  if (dt.gt.dtold .or. nstep_hc.gt.10) then
                     dtold = dt
                     nstep_hc = 0
                     IDecMpRat = IDecMpRat - 1
                  endif
               endif
            endif

            if(swbotb.eq.1 .and. (.not.fllowgwl))then
               ! Derive vertical flux profile in order to find qbot as a
               ! lower boundary condition for the saturated part of the soil
               ! system
               qv(1) = qtop
               do i=NN+1,numnod
                  theta(i) = cofgen(2,i)
               end do
               do i=1,numnod
                 qv(i+1) = qv(i) +dz(i)*FrArMtrx(i)*(theta(i)-thetm1(i))&
     &                          / dt + sink(i) - source(i) + qrot(i)
               end do
               qbot = qv(numnod+1)

               do i=NN+1,numnod
                  h(i) = h(i-1) + disnod(i)*(qv(i)/kmean(i)+1.0d0)
               end do

            end if

            ! Calculate new groundwater level
            call calcgwl ()

            if(swbotb.ne.1.and.abs(gwl-gwlm1).ge.gwlconv .AND.          &
     &         abs(gwl-999d0).gt.1.d0.and.abs(gwlm1-999d0).gt.1.d0) then
               call dtdpst ('year-month-day',t1900+1.001d0,datetmp)
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
      if (.not.fldtmin ) then
         ! Reset soil state variables
         do j = 1,numnod
            h(j) = hm1(j)
            theta(j) = thetm1(j)
         enddo
         kmean(numnod+1) = k(numnod)
         gwl    = gwlm1
         pond   = pondm1

         ! Reset and continue iteration with smaller timestep!
         fldecdt = .true.

         return

         ! In case of macropores, retry with reduction of exchange fluxes with matrix
      elseif (FlMacropore .and. IDecMpRat .lt. 3) then
         ! Reset soil state variables
         IDecMpRat  = IDecMpRat + 1
         FlDecMpRat = .true.
         dtold = dt

!         write(104,'(f10.6,i5)') t, IDecMpRat

         return

      else
         ! Write warning to screen and log file
         if (flwarn_hc .and. iwarn_hc.lt.5) then
            iwarn_hc = iwarn_hc + 1
            call dtdpst                                                 &
     &        ('year-month-day,hour:minute:seconds',t1900,datetime)
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
           write(logf,'(a,f14.6)') 't1900    = ', t1900
           write(logf,'(a,f10.6)') 'dtmin = ', dtmin
           write(logf,'(a,f10.6)') 'dt    = ', dt
           write(logf,'(a,i3)')    'ftoph  = ', ftoph
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

           write(logf,'(a,f10.6)') 'pondm1 = ', pondm1
           write(logf,'(a,f10.6)') 'pond   = ', pond
           write(logf,'(a,f10.6)') 'gwlm1 = ', gwlm1
           write(logf,'(a,f10.6)') 'gwl   = ', gwl
           write(logf,'(a)')                                            &
     &      'node,flnonconv1_F, F, flnonconv2_h, hm1,     h,'//         &
     &      '    thetm1,    theta'
           do j = 1,numnod
              write(logf,'(2(i4,a),f10.6,a,i3,5(a,f10.6))')             &
     &            j,',',flnonconv1(j),',',F(j),',',flnonconv2(j),',',   &
     &            h(j),',',hm1(j),',', theta(j),',',thetm1(j)
           enddo

         endif

         ! Continue without convergence !!!
         return

      endif

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
      use macropore_mod, only: macropore
      use soilwaterbalance_mod, only: calcgwl, watstor, integral, fluxes
      use swap_state_mod, only: swap_state_t
      implicit none

      ! Arguments
      integer task
      type(swap_state_t), intent(in) :: state

      ! Local variables
      integer lay,node,i,j

      real(8) tab(mabbc*2)
      character(len=200) messag

      select case (task)
      case (1)

         ! Initialize Soilwater rate/state variables

         ! Initialize miscellaneous
         hatm = -2.75d+05
      nraidt = 0.0d0
      nird = 0.0d0
      if (swinco.ne.3) then
        ldwet = 0.d0
        spev = 0.d0
        saev = 0.d0
      endif
      runon = 0.0d0
      qtop = 0.d0
      do i = 1,numnod+1
        q(i) = 0.0d0
      enddo
      do i = 1,macp
        evp(i) = 0.0d0
        rfcp(i) = 1.0d0
      enddo
      vtair = 0.0d0
      cQMpLatSs = 0.0d0

      ! Soil physics: tabulated or MualemVanGenuchten functions
      cofgen = 0.0d0
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
          ! Assign values to cofgen
          cofgen(1,node) = 0.0d0                          ! thetar
          cofgen(2,node) = sptab(2,node,numtab(node))     ! thetas
          cofgen(3,node) = sptab(3,node,numtab(node))     ! ksat
          if (do_ln_trans) cofgen(3,node) = dexp(cofgen(3,node))
        end do
        do lay = 1,numlay
          ksatfit(lay) = cofgen(3,nod1lay(lay))
          thetsl(lay) = cofgen(2,nod1lay(lay))
        end do
      else
         ! MvanG functions
         do node = 1,numnod
          lay = layer(node)
          do i = 1, 10
            cofgen(i,node) = paramvg(i,lay)
          end do
          ! Assign dummy value to alphaw
          cofgen(8,node) = -9999.9d0
          if (cofgen(10,node) > 0.0d0) fluseksatexm(node) = .true.
          cofgen(11,node) = relsatthr(lay)
          cofgen(12,node) = ksatthr(lay)
          if (iHWCKmodel(lay) ==  3 .OR. iHWCKmodel(lay) ==  6 .OR. iHWCKmodel(lay) ==  7 .OR. &
              iHWCKmodel(lay) == 10 .OR. iHWCKmodel(lay) == 11) then
             cofgen(13:17,node) = paramvg(13:17,lay)
          end if
          if (iHWCKmodel(lay) ==  5 .OR. iHWCKmodel(lay) ==  7) then
             cofgen(18,node) = paramvg(18,lay)
          end if
          if (iHWCKmodel(lay) ==  8 .OR. iHWCKmodel(lay) ==  9 .OR. &
              iHWCKmodel(lay) == 10 .OR. iHWCKmodel(lay) == 11) then
             cofgen(18:21,node) = paramvg(18:21,lay)
          end if
        end do
        thetsl = 0.0d0
        do lay = 1, numlay
          thetsl(lay) = paramvg(2,lay)
        end do
      endif

! --- saturated and residual watercontent of each node; hysteresis parameters
      do node = 1,numnod
        lay = layer(node)
        thetar(node) = cofgen(1,node)
        thetas(node) = cofgen(2,node)
        !!! Kroes: disable combi of swsophy=1 and swhyst=1
        if (swhyst.eq.1) then
           ! Wetting curve
           indeks(node) = 1
           cofgen(4,node) = paramvg(8,lay)
        elseif (swhyst.eq.0.or.swhyst.eq.2) then
           ! Drying branch or simulation without hysteresis
           indeks(node) = -1
           cofgen(4,node) = paramvg(4,lay)
        endif
      end do

      ! Additional Input checks
      ! ThetCrMP should be less then CofGen(2,Lay): requires additional inputcheck
      if (FlMacropore) then
         do lay = 1, NumLay
            if (SwSoilShr(lay).gt.0 .and.                               &
     &          ThetCrMp(lay).gt.CofGen(2,nod1lay(lay))) then
               messag = ' ThetCrMP.gt.ThetSat'
               call fatalerr_collected('MacroRead',messag)
            endif
         enddo
      endif

      if (swinco.eq.1) then
         ! Pressure head profile is input
         do i = 1, nhead
          tab(i*2) = h(i)
          tab(i*2-1) = abs(zi(i))
        end do
        do i = 1, numnod
          h(i) = afgen(tab,macp*2,abs(z(i)))
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
         ! Determine groundwater level
         if (h(numnod) .gt. -1.d-5) then
          i = numnod
          do while ((h(i) .gt. -1.d-5) .and. (i .gt. 1))
              i = i - 1
          end do
          if (h(i) .lt. -1.d-5) then
            gwl = z(i+1) + h(i+1) / (h(i+1) - h(i)) * (z(i) - z(i+1))
            ! Assume hydrostatic equilibrium in saturated part
            do j = i+1, numnod
              h(j) = gwl - z(j)
            end do
          endif
        endif
      else
         ! Pressure head profile is calculated from groundwater level
         if (swbotb.eq.1) then  
          gwl = afgen (gwltab,mabbc*2,t1900+dt-1.d0)

          if(abs(gwl-(z(numnod)-0.5d0*dz(numnod))) .lt.1.0d-4) then
          messag = 'Groundwaterlevel as bottom boundary (SWBOTB=1) is'//&
     &    'below or to close to bottom of soil profile'//               &
     &    ' must be corrected!'
            call fatalerr_collected ('soilwater',messag)
          endif
        else
          gwl = gwli
        endif
        if (gwl.gt.0.0d0) then 
          pond = gwl
        else
          pond = 0.0d0
        endif
        do i = 1,numnod
          h(i) = gwl - z(i)
        end do
      endif

      ! In case of preferential flow, adjust Van Genuchten parameters
      do i = 1, numnod
        theta(i) = watcon(i,h(i))
      end do

      ! Hydraulic conductivities, differential moisture capacities
      ! and mean hydraulic conductivities for each node
      do node = 1,numnod
        dimoca(node) = moiscap(node,h(node))

        FrArMtrx(node) = 1.d0
        k(node) = hconduc (node,h(node),theta(node),rfcp(node))
        if(FlMacropore)  k(node) = FrArMtrx(node) * k(node)

        if(node.gt.1) kmean(node) =  hcomean(swkmean,k(node-1),k(node),dz(node-1),dz(node))
      end do
      kmean(numnod+1) = k(numnod)

      ! Initial soil water storage
      if (.not.flMacroPore) then
         do i = 1, NumNod
            FrArMtrx(i) = 1.d0
         enddo
         volact = 0.0d0
         call watstor ()
         volini = volact
      endif
      pondini = pond

      ! Initial groundwater level
      call calcgwl ()

      call log_info('soilwater', 'Soil state initialized: gwl=' // to_str(real(gwl,4)) // &
                    ' cm, numnod=' // to_str(numnod) // ', numlay=' // to_str(numlay) // &
                    ', volini=' // to_str(real(volini,4)) // ' cm')

      return

      case (2)

         ! Calculate Soilwater rate/state variables

         ! Reset intermediate soil water fluxes
         if (flDayStart) then
          iqredwet_day = 0.0d0
          iqreddry_day = 0.0d0
          iqredsol_day = 0.0d0
          iqredfrs_day = 0.0d0
          iptra_day    = 0.0d0
          do node = 1,numnod 
            qpotrot_day(node) = 0.d0
            qredtot_day(node) = 0.d0
          enddo
      end if
      
      if (flzerointr) then
        do node = 1,numnod
          inqrot(node) = 0.0d0
          inqssdi(node) = 0.0d0
          inq(node) = 0.0d0
        enddo
        inq(numnod+1) = 0.0d0
        iqrot = 0.0d0
        iqssdi = 0.0d0
        iqredwet = 0.0d0
        iqreddry = 0.0d0
        iqredsol = 0.0d0
        iqredfrs = 0.0d0
        ies0 = 0.0d0
        iet0 = 0.0d0
        iew0 = 0.0d0
        iintc = 0.0d0
        iptra = 0.0d0
        ipeva = 0.0d0
        ievap = 0.0d0
        iruno = 0.0d0
        irunoCN = 0.0d0
        iqbot = 0.0d0
        iqtdo = 0.0d0
        iqtup = 0.0d0
        irunon = 0.0d0
        iqdo(1:numnod+1) = 0.0d0
        iqup(1:numnod+1) = 0.0d0

        IPondBeg = Pond
        do node = 1, numnod
          IThetaBeg(node) = Theta(node)
        enddo

        ! Macropore variables
        if (flMacroPore) call macropore(5)
      endif

      ! Reset cumulative soil water fluxes
      if (flzerocumu) then
        cqssdi = 0.0d0
        cqrot = 0.0d0
        cqbot = 0.0d0
        cqbotdo = 0.0d0
        cqbotup = 0.0d0
        cptra = 0.0d0
        cpeva = 0.0d0
        cevap = 0.0d0
        cinund = 0.0d0
        crunon = 0.0d0
        crunoff = 0.0d0
        crunoffCN = 0.0d0
        cqtdo = 0.0d0
        cqtup = 0.0d0
        cqprai = 0.0d0

        ! Macropore variables
        if (flMacroPore) call macropore(6)

        ! Reset initial water storage and ponding
        volini = volact
        pondini = pond
      endif

      ! Save state variables of time = t
      call SoilWaterStateVar(1)

      ! Calculate new soil water state variables
      call headcalc(state)

      return

      case (3)

         ! Update hydraulic conductivities to time level t+1
         do i = 1,numnod
         k(i) = hconduc(i,h(i),theta(i),rfcp(i))
         if(FlMacropore)  k(i) = FrArMtrx(i) * k(i)
         if(i.gt.1)then
            kmean(i) = hcomean(swkmean,k(i-1),k(i),dz(i-1),dz(i))
         end if
      enddo
      kmean(numnod+1) = k(numnod)

      ! Calculate actual water content of profile
      call watstor ()

      ! Calculate water fluxes between soil compartments
      call fluxes ()

      ! Calculation of states macropores and intermediate & cumulative values
      if (flMacroPore) call macropore(4)

      ! Calculate cumulative fluxes
      call integral

      ! Update parameters for soil water hystereses
      if (swhyst.ne.0) call hysteresis ()

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
   subroutine SoilWaterStateVar(task)

! --- global variables
      use Variables
      implicit none

      ! Local variables
      integer task, i

      select case (task)
      case (1)

         ! Save state variables of time = t
         do i = 1,numnod
        hm1(i) = h(i)
        thetm1(i) = theta(i)
      enddo
      gwlm1 = gwl
      pondm1 = pond

      return

      case (2)

         ! Reset soil state variables
         do i = 1,numnod
        h(i) = hm1(i)
        theta(i) = thetm1(i)
      enddo
      kmean(numnod+1) = k(numnod)
      gwl    = gwlm1
      pond   = pondm1

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
   subroutine hysteresis ()
      use variables, only: numnod,layer,h,hm1,indeks,tau,paramvg,cofgen,dimoca,theta,disnod
      use soilhydraulics_utils, only: moiscap, prhead
      use swap_array_dimensions, only: macp

      implicit none

      ! Local variables
      integer node,lay,indtem(macp)
      real(8) delp,sew,sed,fvalue
      real(8) thetar(macp),thetas(macp),alfamg(macp)

      ! Check for reversal
      do node = 1,numnod
        delp = hm1(node)-h(node)
        if (delp/float(indeks(node)).gt.tau.and.                        &
     &     h(node).lt.-10.0d0 .and. h(node).gt.-1.0d3) then
          indtem(node) = -indeks(node)
        else
          indtem(node) = indeks(node)
        endif
      end do

      ! Change parameters scanning curves
      do 100 node = 1,numnod
        lay = layer(node)

         ! No change
         if (indtem(node).eq.indeks(node) .or.                           &
     &     abs(paramvg(4,lay)-paramvg(8,lay)) .lt. 1.d-4) goto 100

         ! Relative saturation
         sew = (1.0d0+(paramvg(8,lay)*(-h(node)))**paramvg(6,lay))       &
     &        **(-paramvg(7,lay)) 
        sed = (1.0d0+(paramvg(4,lay)*(-h(node)))**paramvg(6,lay))       &
     &        **(-paramvg(7,lay))

         ! Change index
         indeks(node) = -1*indeks(node)

         ! Update alfa, thetar and thetas
         if (indeks(node).eq.1) then
            ! Wetting branch
            alfamg(node) = paramvg(8,lay)
          thetas(node) = paramvg(2,lay)
          thetar(node) = (theta(node)-thetas(node)*sew)/(1.0d0-sew)

            ! Check on thetar(node) value, if needed correction of h
            fvalue = thetar(node)
          if(thetar(node).lt.paramvg(1,lay)) thetar(node)=paramvg(1,lay)
          if(thetar(node).gt.paramvg(2,lay)) thetar(node)=paramvg(2,lay)
          cofgen(4,node) = alfamg(node)
          cofgen(1,node) = thetar(node)
          cofgen(2,node) = thetas(node)
          if (abs(fvalue-thetar(node)) .gt. 1.d-10) then
             h(node) = prhead(node,disnod(node),theta(node),cofgen,h)
          endif
        else
           ! Drying branch
           alfamg(node) = paramvg(4,lay)
          thetar(node) = paramvg(1,lay)
          thetas(node) = thetar(node)+(theta(node)-thetar(node))/sed

            ! Check on thetas(node) value, if needed correction of h
            fvalue = thetas(node)
          if(thetas(node).lt.paramvg(1,lay)) thetas(node)=paramvg(1,lay)
          if(thetas(node).gt.paramvg(2,lay)) thetas(node)=paramvg(2,lay)
          cofgen(4,node) = alfamg(node)
          cofgen(1,node) = thetar(node)
          cofgen(2,node) = thetas(node)
          if (abs(fvalue-thetas(node)) .gt. 1.d-10) then
             h(node) = prhead(node,disnod(node),theta(node),cofgen,h)
          endif
        endif

         ! Update capacity
         dimoca(node) = moiscap(node,h(node))
 100  continue

      return
      end subroutine hysteresis

end module