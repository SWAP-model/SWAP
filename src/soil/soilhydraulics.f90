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
   public :: headcalc, soilwater_seed, soilwater_step, soilwater_update, soilwater_save_state, soilwater_restore_state, hysteresis
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

      use swap_array_dimensions, only: macp, mabbc
      use swap_log, only: log_warn, log_debug, to_str
      use boundbottom_mod, only: BoundBottom
      use boundtop_mod, only: PONDRUNOFF
      use rootextraction_mod, only: RootExtraction
      use soilhydraulics_utils, only: watcon, hconduc, moiscap, hcomean, dhconduc
      use swap_constants, only: nihil
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
      ! local loop variable; copied to soil%numbit at top of loop body
      ! (Fortran disallows derived-type components as loop variables).
      integer   numbit_local

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
      real(8) q0       ! surface flux for boundary (local; used only in swbotb=1 branch)
      real(8) ArMpSs   ! macropore area fraction at soil surface (local; always 0.d0, ADR 0040)
      logical flok

      associate (mesh => state%mesh,         &
                 soil => state%soilwater,    &
                 drai => state%drainage,     &
                 heat => state%heat,         &
                 atmo => state%atmosphere,   &
                 time => state%timecontrol,  &
                 swbotb => state%soilwater%swbotb_runtime)

      if (time%flDayStart) then
         soil%flwarn_hc = .true.
         soil%iwarn_hc = 0
      endif
      call dtdpst                                                       &
     &        ('year-month-day,hour:minute:seconds',time%t1900,datetime)

      ! Summation of sink terms (constant for the current time step)
      iBackTr   = 0
      flunsatok(1) = .false.
      flunsatok(2) = .false.
      flunsatok(3) = .false.

      do i=1,mesh%numnod
         sink(i) = soil%evp(i)
         if (allocated(drai%qdra)) then
            do j=1,drai%nrlevs
               sink(i) = sink(i) + drai%qdra(j,i)
            end do
         end if
      end do
      source(1:mesh%numnod) = soil%qssdi(1:mesh%numnod)

      ArMpSs = 0.d0                                            ! macropore retired (ADR 0040)

      ! Groundwater level specified
      if(swbotb.eq.1)then
         soil%fllowgwl = .false.
         if(soil%gwlinp.ge.mesh%z(1)-1.0d-4)then

            q0 = (atmo%nraidt+atmo%nird+atmo%melt)*(1.0d0-ArMpSs) + soil%runon - soil%reva
            call pondrunoff (state)
            q1 = - q0 + (soil%pond - soil%pondm1)/time%dt + soil%runots / time%dt
            soil%theta(1) = watcon(soil%gwlinp, &
                                  soil%vg_params(1), &
                                  soil%iHWCKmodel(soil%layer(1)), &
                                  1, soil)
            soil%kmean(1) = hconduc(soil%gwlinp,soil%theta(1),heat%rfcp(1),heat%tsoil(1), &
                                  soil%vg_params(1), &
                                  soil%iHWCKmodel(soil%layer(1)), &
                                  soil%fluseksatexm(1), &
                                  1, soil)

            qv(1) = q1
            do i=1,mesh%numnod
               qv(i+1) = qv(i) +mesh%dz(i)*soil%FrArMtrx(i)*(soil%theta(i)-soil%thetm1(i))  &
     &                          / time%dt+ sink(i) - source(i) + soil%qrot(i)
            end do
            soil%qbot = qv(mesh%numnod+1)
            soil%h(1) = soil%gwlinp + mesh%disnod(1)*(qv(1)/soil%kmean(1)+1.0d0)
            do i=2,mesh%numnod
               soil%h(i) = soil%h(i-1) + mesh%disnod(i)*(qv(i)/soil%kmean(i)+1.0d0)  
            end do

            if(soil%swkimpl.eq.1)then
               do i=1,mesh%numnod
                  soil%k(i) = hconduc(soil%h(i),soil%theta(i),heat%rfcp(i),heat%tsoil(i), &
                                    soil%vg_params(i), &
                                    soil%iHWCKmodel(soil%layer(i)), &
                                    soil%fluseksatexm(i), &
                                    i, soil)
                  if(i.gt.1)then
                     soil%kmean(i) = hcomean(soil%swkmean,soil%k(i-1),soil%k(i),mesh%dz(i-1),mesh%dz(i))
                  end if
               end do
               soil%kmean(mesh%numnod+1) = soil%k(mesh%numnod)
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
               soil%fllowgwl = .true.
               soil%hbot = soil%gwlinp - mesh%z(mesh%numnod) + 0.5*mesh%dz(mesh%numnod)
            endif
         end if
      else
         NN = mesh%numnod
      end if

      ! Reset conductivities to time level t

      ! Node nr of compartment with minimized flux of capillary rise
      if (soil%swcaprise) then
         nodncr    = max(5, state%crop%common%noddrz)
         flcaprise = .false.
      endif
      do i = 1,mesh%numnod
         soil%k(i) = hconduc(soil%h(i),soil%theta(i),heat%rfcp(i),heat%tsoil(i), &
                           soil%vg_params(i), &
                           soil%iHWCKmodel(soil%layer(i)), &
                           soil%fluseksatexm(i), &
                           i, soil)

         if (soil%swcaprise) then
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
            soil%kmean(i) = hcomean(soil%swkmean,soil%k(i-1),soil%k(i),mesh%dz(i-1),mesh%dz(i))
         end if
      enddo
      soil%kmean(mesh%numnod+1) = soil%k(mesh%numnod)

      if(soil%swkimpl.eq.0)then
         do i=2,mesh%numnod
            dFdhU(i)   = - soil%kmean(i)  /mesh%disnod(i)
            dFdhL(i-1) = dFdhU(i)
         end do
      end if

      ! Lysimeter plate-contact flag: based on pre-iteration h(NN); computed here (not in the helper) so boundtop cannot change it mid-step.
      flboth = (swbotb .eq. 8 .and. soil%h(NN) .gt. Critdz - mesh%disnod(NN+1) + soil%hplate)

      call headcalc_residual(state, NN, sink, source, flboth, hgrad, F)

      ! Initial estimate of F inner product
      sumold    = 0.0d0
      do i=1,NN
         sumold = sumold + F(i)*F(i)
      end do
      sumold = 0.5d0 * sumold

      ! Start iteration loop, time%MaxIt specified in the input
      if(time%fldtmin)then
         MaxIt1 = 2*time%MaxIt
      else
         MaxIt1 = time%MaxIt
      end if

      sum= 0.d0   ! For Forcheck
      Do numbit_local = 1, MaxIt1
         soil%numbit = numbit_local

         do i = 1, NN
            ! Save values of h
            hold(i)   = soil%h(i)

            ! Derivative of theta to h (differential moisture capacity),
            ! as part of main diagonal
            soil%dimoca(i) = moiscap(soil%h(i), &
                                    soil%vg_params(i), &
                                    soil%iHWCKmodel(soil%layer(i)), &
                                    time%dt, &
                                    i, soil)

         enddo

         if(soil%swkimpl.eq.1)then
            do i = 1, NN
               dkdh(i) = dhconduc(soil%h(i),soil%theta(i),soil%dimoca(i),heat%rfcp(i), &
                                   soil%vg_params(i), &
                                   soil%iHWCKmodel(soil%layer(i)), &
                                   i, soil)
            enddo
            do i=2,NN
               dFdhU(i)   = - soil%kmean(i)  /mesh%disnod(i)
               dFdhL(i-1) = dFdhU(i)
            end do
         end if

         if (soil%swcaprise .and. flcaprise) then
           dkdh(nodncr+1) = 1.0D-30
         endif

         ! Jacobian matrix elements
         dFdhM(1) = soil%dimoca(1)*soil%FrArMtrx(1)*mesh%dz(1)/time%dt - dFdhL(1)
         ! If the head boundary condition applies: add the k1/(0.5*dz1) term
         ! to the first element of the main diagonal
         if(soil%ftoph) dFdhM(1) = dFdhM(1) + soil%kmean(1)/mesh%disnod(1)

         do i=2,NN-1
            dFdhM(i) = soil%dimoca(i)*soil%FrArMtrx(i)*mesh%dz(i)/time%dt - dFdhU(i)        &
     &                                                - dFdhL(i)
         end do

         dFdhM(NN) = soil%dimoca(NN)*soil%FrArMtrx(NN)*mesh%dz(NN)/time%dt - dFdhU(NN)
         if(swbotb.eq.1 .and. (.not.soil%fllowgwl))then
            dFdhM(NN) = dFdhM(NN) + soil%kmean(NN+1)/(mesh%z(NN)-soil%gwlinp)
         else if(swbotb.eq.3.and.soil%swbotb3impl.eq.1)then ! Cauchy
            if (soil%swbotb3resvert.eq.0) then
               dFdhM(NN) = dFdhM(NN) + 1.0d0 /                          &
     &                              (mesh%disnod(NN+1)/soil%kmean(NN+1)+soil%rimlay)
            elseif (soil%swbotb3resvert.eq.1) then
               dFdhM(NN) = dFdhM(NN) + 1.0d0 / soil%rimlay
            endif
         else if(swbotb.eq.5 .or. (swbotb.eq.1 .and. soil%fllowgwl))then
            dFdhM(NN) = dFdhM(NN) + soil%kmean(NN+1)/mesh%disnod(NN+1)         
         else if(swbotb.eq.7 .or. swbotb .eq. -2)then ! implicitly: soil%kmean(NN+1)
            dFdhM(NN) = dFdhM(NN) + dkdh(NN) * 0.5d0
         else if(swbotb.eq.8 .and. flboth)then
            dFdhM(NN) = dFdhM(NN) + soil%kmean(NN+1)/mesh%disnod(NN+1)
         end if

         if(soil%swkimpl.eq.1)then
            dFdhM(1) = dFdhM(1) + dkdh(1) * hgrad(2) *                  &
     &                 dkmean(soil%swkmean,soil%k(1),soil%k(2),mesh%dz(1),mesh%dz(2))
            if(soil%ftoph) dFdhM(1) = dFdhM(1) - dkdh(1) * hgrad(1) * 0.5d0
            dFdhL(1) = dFdhL(1) + dkdh(2) * hgrad(2) *                  &
     &                 dkmean(soil%swkmean,soil%k(2),soil%k(1),mesh%dz(2),mesh%dz(1))
            do i=2,NN-1
               dFdhU(i) = dFdhU(i) - dkdh(i-1) * hgrad(i) *             &
     &                    dkmean(soil%swkmean,soil%k(i-1),soil%k(i),mesh%dz(i-1),mesh%dz(i))
               dFdhM(i) = dFdhM(i) - dkdh(i) * hgrad(i) *               &
     &                    dkmean(soil%swkmean,soil%k(i),soil%k(i-1),mesh%dz(i),mesh%dz(i-1))     &
     &                             + dkdh(i) * hgrad(i+1) *             &
     &                    dkmean(soil%swkmean,soil%k(i),soil%k(i+1),mesh%dz(i),mesh%dz(i+1))
               dFdhL(i) = dFdhL(i) + dkdh(i+1) * hgrad(i+1) *           &
     &                    dkmean(soil%swkmean,soil%k(i+1),soil%k(i),mesh%dz(i+1),mesh%dz(i))
            end do
            dFdhU(NN) = dFdhU(NN) - dkdh(NN-1) * hgrad(NN) *            &
     &                  dkmean(soil%swkmean,soil%k(NN-1),soil%k(NN),mesh%dz(NN-1),mesh%dz(NN))
            dFdhM(NN) = dFdhM(NN) - dkdh(NN) * hgrad(NN) *              &
     &                  dkmean(soil%swkmean,soil%k(NN),soil%k(NN-1),mesh%dz(NN),mesh%dz(NN-1))

            if(swbotb.eq.1 .or. swbotb.eq.5 .or. swbotb.eq.8            &
     &         .and. flboth)then
               dFdhM(NN) = dFdhM(NN) + 0.5d0 * dkdh(NN) * hgrad(NN+1)
            end if
         end if

         ! Solve the tridiagonal matrix
         call tridag(NN, dFdhU, dFdhM, dFdhL, F, difh, ierror)

         if(ierror.ne.0)then
            call dtdpst ('year-month-day',time%t1900+1.001d0,datetmp)
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
         do itry = 1,soil%MaxBackTr
            iBackTr = iBackTr + 1
            ! Factor reduces the change of h (difh) calculated as a full
            ! Newton Raphson step
            if(time%fldtmin .and. soil%numbit.gt.time%MaxIt)then
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
               end do
            else
               do i = 1,NN
                  soil%h(i) = hold(i) - factor * difh(i)
               end do
            end if
            do i = 1,NN
              soil%theta(i) = watcon(soil%h(i), &
                                   soil%vg_params(i), &
                                   soil%iHWCKmodel(soil%layer(i)), &
                                   i, soil)
            enddo

            if(soil%swkimpl.eq.1)then
               call Rootextraction(state)
               do i = 1,NN
                  soil%k(i) = hconduc(soil%h(i),soil%theta(i),heat%rfcp(i),heat%tsoil(i), &
                                    soil%vg_params(i), &
                                    soil%iHWCKmodel(soil%layer(i)), &
                                    soil%fluseksatexm(i), &
                                    i, soil)
                  if(i.gt.1)then
                     soil%kmean(i)=hcomean(soil%swkmean,soil%k(i-1),soil%k(i),mesh%dz(i-1),mesh%dz(i))
                  end if
               end do
               soil%kmean(NN+1) = soil%k(NN)
            end if

            ! Prevent capillary rise into the root zone !! special for experts
            if (soil%swcaprise) then
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
                 soil%kmean(i)=hcomean(soil%swkmean,soil%k(i-1),soil%k(i),mesh%dz(i-1),mesh%dz(i))
                 soil%kmean(i+1)=hcomean(soil%swkmean,soil%k(i),soil%k(i+1),mesh%dz(i),mesh%dz(i+1))
               end if
            endif

            ! Calculate F-function (residual + hgrad + bottom-boundary dispatch)
            call headcalc_residual(state, NN, sink, source, flboth, hgrad, F)

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

         ! Apply performance criteria per compartment
         do i = 1,NN
            ! Test for water balance deviation of soil compartments
            if( dabs(F(i)).gt. CritDevBalCp)then
                  flnonconv1(i) = .true. ; flnonconv   = .true.
            end if
            ! Test for change of pressure head
            if( dabs(hold(i)) .lt. 1.0d0)then
               if(abs( soil%h(i)-hold(i) ) .gt. soil%critdevh2cp)then
                  flnonconv2(i) = .true. ; flnonconv   = .true.
               endif
            else
               if(abs( soil%h(i)-hold(i) )/abs(hold(i)) .gt. soil%critdevh1cp)then
                  flnonconv2(i) = .true. ; flnonconv   = .true.
               endif
            end if

         enddo

         ! Test for waterbalance of ponding layer
         if (soil%ftoph) then
            soil%qtop = -soil%kmean(1)*((soil%hsurf - soil%h(1))/mesh%disnod(1)+1.0d0)
            if(.not.flnonconv) then
               deviat = soil%pond - soil%pondm1 + soil%reva*time%dt - (atmo%nraidt+atmo%nird+atmo%melt)*time%dt &
     &                - soil%runon*time%dt  +  soil%runots  - soil%qtop * time%dt
               if( abs(deviat) .gt. soil%critdevponddt) then
                  flnonconv3 = .true. ; flnonconv   = .true.
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
                  soil%theta(i) = soil%vg_params(i)%thetas
               end do
               do i=1,mesh%numnod
                 qv(i+1) = qv(i) +mesh%dz(i)*soil%FrArMtrx(i)*(soil%theta(i)-soil%thetm1(i))&
     &                          / time%dt + sink(i) - source(i) + soil%qrot(i)
               end do
               soil%qbot = qv(mesh%numnod+1)

               do i=NN+1,mesh%numnod
                  soil%h(i) = soil%h(i-1) + mesh%disnod(i)*(qv(i)/soil%kmean(i)+1.0d0)
               end do

            end if

            ! Calculate new groundwater level
            call calcgwl (state)

            if(swbotb.ne.1.and.abs(soil%gwl-soil%gwlm1).ge.soil%gwlconv .AND.          &
     &         abs(soil%gwl-999d0).gt.1.d0.and.abs(soil%gwlm1-999d0).gt.1.d0) then
               call dtdpst ('year-month-day',time%t1900+1.001d0,datetmp)
                 messag = ' Change of groundwater level exceeds'//      &
     &           ' criterion at '//datetmp//'. Consider reduction of dtMin'
               call log_warn('Headcalc', messag)
            endif

            ! Recording of number of iteration steps needed
            soil%Itnumb(min(100,soil%numbit),1)=soil%Itnumb(min(100,soil%numbit),1)+1 
            soil%Itnumb(min(100,soil%numbit),2)=soil%Itnumb(min(100,soil%numbit),2)+iBackTr 

            return
            
         end if
      End Do

      ! Convergence could not been reached
      if (.not.time%fldtmin ) then
         ! Reset soil state variables
         do j = 1,mesh%numnod
            soil%h(j) = soil%hm1(j)
            soil%theta(j) = soil%thetm1(j)
         enddo
         soil%kmean(mesh%numnod+1) = soil%k(mesh%numnod)
         soil%gwl  = soil%gwlm1
         soil%pond = soil%pondm1

         ! Reset and continue iteration with smaller timestep!
         time%fldecdt = .true.

         return

      else
         ! Write warning to screen and log file
         if (soil%flwarn_hc .and. soil%iwarn_hc.lt.5) then
            soil%iwarn_hc = soil%iwarn_hc + 1
            call dtdpst                                                 &
     &        ('year-month-day,hour:minute:seconds',time%t1900,datetime)
            messag = ' No convergence was reached of Richards'//        &
     &        ' equation at '//datetime//                               &
     &        ' no more than 4 warnings per date - SWAP did continue !'
            call log_warn('Headcalc', messag)
            if (soil%iwarn_hc.gt.4) then
              soil%flwarn_hc = .false.  
            endif
         endif

         if (soil%dump_convergence_diagnostics) then
           call log_debug('Headcalc', 'Datetime = ' // datetime)
           call log_debug('Headcalc', 't1900    = ' // to_str(time%t1900))
           call log_debug('Headcalc', 'time%dtmin = ' // to_str(time%dtmin))
           call log_debug('Headcalc', 'dt    = ' // to_str(time%dt))
           call log_debug('Headcalc', 'ftoph  = ' // to_str(soil%ftoph))
           call log_debug('Headcalc', 'CritDevBalCp  = ' // to_str(CritDevBalCp))
           call log_debug('Headcalc', 'CritDevBalTot = ' // to_str(CritDevBalTot))
           call log_debug('Headcalc', 'CritDz        = ' // to_str(CritDz))
           call log_debug('Headcalc', 'critdevh1cp   = ' // to_str(soil%critdevh1cp))
           call log_debug('Headcalc', 'critdevh2cp   = ' // to_str(soil%critdevh2cp))
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

      end associate

   end subroutine headcalc

   !> Compute the Richards residual vector F (+ hgrad and bottom-boundary dispatch)
   !!
   !! Shared by headcalc's initial estimate and its back-tracking iteration loop.
   !! Computes hgrad(2..NN), the surface residual F(1) (via boundtop / pondrunoff
   !! / ftoph), the interior residual F(2..NN-1), the lower-boundary gradient
   !! hgrad(NN+1) and the bottom-boundary residual F(NN) for the selected swbotb.
   !! The lysimeter plate-contact flag flboth is decided by the caller and passed
   !! in (it depends on the pre-iteration h(NN), which boundtop does not change).
   !!
   subroutine headcalc_residual(state, NN, sink, source, flboth, hgrad, F)

      use boundtop_mod, only: boundtop, PONDRUNOFF
      use array_utils, only: afgen
      use soilhydraulics_utils, only: watcon, hconduc, hcomean
      use swap_array_dimensions, only: mabbc
      use swap_state_mod, only: swap_state_t
      implicit none

      type(swap_state_t), intent(inout) :: state
      integer,            intent(in)    :: NN
      real(8),            intent(in)    :: sink(:), source(:)
      logical,            intent(in)    :: flboth
      real(8),            intent(inout) :: hgrad(:)
      real(8),            intent(out)   :: F(:)

      integer :: i

      associate (mesh => state%mesh,         &
                 soil => state%soilwater,    &
                 heat => state%heat,         &
                 time => state%timecontrol,  &
                 swbotb => state%soilwater%swbotb_runtime)

      do i=2,NN
         hgrad(i) = (soil%h(i-1)-soil%h(i))/mesh%disnod(i) + 1.0d0
      end do

      F(1) = (soil%theta(1)-soil%thetm1(1))*soil%FrArMtrx(1)*mesh%dz(1)/time%dt + sink(1) - source(1) + soil%qrot(1) + soil%kmean(2) * hgrad(2)

      call boundtop(state)

      if (soil%FlRunoff) call pondrunoff (state)

      if(soil%ftoph)then
         hgrad(1) = (soil%hsurf-soil%h(1))/mesh%disnod(1) + 1.d0
         F(1)     = F(1) - soil%kmean(1) * hgrad(1)
      else
         F(1) = F(1) + soil%qtop
      end if

      do i=2,NN-1
         F(i) = (soil%theta(i)-soil%thetm1(i))*soil%FrArMtrx(i)*mesh%dz(i)/time%dt + sink(i) - source(i) +     &
     &          soil%qrot(i) - soil%kmean(i) * hgrad(i) + soil%kmean(i+1) * hgrad(i+1)
      end do

      if(swbotb.eq.1 .and. (.not.soil%fllowgwl))then
         hgrad(NN+1) = soil%h(NN)/(mesh%z(nn)-soil%gwlinp) + 1.0d0
      else if(swbotb.eq.5 .or. (swbotb.eq.1 .and. soil%fllowgwl))then
         hgrad(NN+1) = (soil%h(NN) - soil%hbot) / mesh%disnod(NN+1)  + 1.0d0
      else if(swbotb.eq.8 .and. flboth)then
         hgrad(NN+1) = (soil%h(NN) - soil%hplate) / mesh%disnod(NN+1)  + 1.0d0
      end if

      ! NB: term order unified from the two original residual blocks; SWBOTB=1 is dark in the regression suite, so this is byte-identical there (algebraically identical regardless).
      if(swbotb.eq.1 .and. (.not.soil%fllowgwl))then
         soil%theta(NN)= watcon(soil%h(NN), &
                              soil%vg_params(NN), &
                              soil%iHWCKmodel(soil%layer(NN)), &
                              NN, soil)
         soil%k(NN)    = hconduc(soil%h(NN),soil%theta(NN),heat%rfcp(NN),heat%tsoil(NN), &
                              soil%vg_params(NN), &
                              soil%iHWCKmodel(soil%layer(NN)), &
                              soil%fluseksatexm(NN), &
                              NN, soil)
         soil%kmean(NN+1) = hcomean(soil%swkmean,soil%k(NN),soil%vg_params(NN+1)%ksat, &
     &                        mesh%dz(NN),mesh%dz(NN+1))
         F(NN) = (soil%theta(NN) - soil%thetm1(NN))*soil%FrArMtrx(NN)*mesh%dz(NN)/time%dt +   &
                 sink(NN) - source(NN) + soil%qrot(NN) - soil%kmean(NN)*hgrad(NN) +           &
                 soil%kmean(NN+1)*hgrad(NN+1)
      else

         F(NN) = (soil%theta(NN) - soil%thetm1(NN))*soil%FrArMtrx(NN)*mesh%dz(NN)/time%dt        &
     &         - soil%kmean(NN) * hgrad(NN) + sink(NN) - source(NN) + soil%qrot(NN)

         if(swbotb.eq.3.and.soil%swbotb3impl.eq.1)then ! Cauchy-relation, implemented as head boundary
            if (soil%swbotb3resvert.eq.0) then
               soil%qbot = - (soil%h(NN)+mesh%z(NN)-soil%deepgw) / (mesh%disnod(NN+1)/soil%kmean(NN+1)+soil%rimlay)
            elseif (soil%swbotb3resvert.eq.1) then
               soil%qbot = - (soil%h(NN)+mesh%z(NN)-soil%deepgw) / soil%rimlay
            endif
! ---       extra groundwater flux might be added
            if (soil%sw4 .eq. 1) soil%qbot = soil%qbot + afgen(soil%qbotab,mabbc*2,time%t1900+time%dt)
            F(NN) = F(NN) - soil%qbot
         else if(swbotb.eq.5 .or. (swbotb.eq.1 .and. soil%fllowgwl))then ! pressure head at lower boundary specified
            F(NN) = F(NN) + soil%kmean(NN+1) * hgrad(NN+1)
         else if(swbotb.eq.7 .or. swbotb .eq. -2)then ! free drainage option
            soil%kmean(mesh%numnod+1) = hconduc(soil%h(mesh%numnod),soil%theta(mesh%numnod),heat%rfcp(mesh%numnod),heat%tsoil(mesh%numnod), &
                                         soil%vg_params(mesh%numnod), &
                                         soil%iHWCKmodel(soil%layer(mesh%numnod)), &
                                         soil%fluseksatexm(mesh%numnod), &
                                         mesh%numnod, soil)
            soil%qbot = -1.0d0 * soil%kmean(mesh%numnod+1)
            F(NN) = F(NN) - soil%qbot
         ! Lysimeter option
         else if(swbotb.eq.8)then
            if (flboth) then
               soil%hbot = soil%hplate
               F(NN) = F(NN) + soil%kmean(NN+1) * hgrad(NN+1)
            else
               soil%qbot = 0.0d0
            end if
         ! Flux bottom boundary
         else
            F(NN) = F(NN) - soil%qbot
         end if

      end if

      end associate

   end subroutine headcalc_residual

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
   !> Narrative orchestrator for soil water initialization. Each phase is
   !> delegated to a named private helper; the storage/groundwater finalize
   !> is kept inline.
   subroutine soilwater_seed(state, hyd)
      use swap_log, only: log_info, to_str
      use soilwaterbalance_mod, only: calcgwl, watstor
      use swap_state_mod, only: swap_state_t
      use soil_config_mod, only: soil_hydraulics_t
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none

      ! Arguments
      type(swap_state_t),      intent(inout) :: state
      type(soil_hydraulics_t), intent(in)    :: hyd

      ! Local variables
      integer i

      call init_soil_misc(state)
      call populate_hydraulic_params(state, hyd)
      call apply_soil_initial_conditions(state)
      call compute_initial_node_hydraulics(state)

      associate (mesh => state%mesh,         &
                 soil => state%soilwater)

      ! Initial soil water storage
      do i = 1, mesh%numnod
         soil%FrArMtrx(i) = 1.0_real64
      enddo
      soil%volact = 0.0_real64
      call watstor (state)
      soil%volini = soil%volact
      soil%pondini = soil%pond

      ! Initial groundwater level
      call calcgwl (state)

      call log_info('soilwater', 'Soil state initialized: gwl=' // to_str(real(soil%gwl,4)) // &
                    ' cm, numnod=' // to_str(mesh%numnod) // ', numlay=' // to_str(mesh%numlay) // &
                    ', volini=' // to_str(real(soil%volini,4)) // ' cm')

      end associate

      return
      end subroutine soilwater_seed

   !> Initialize miscellaneous soil water rate/state variables.
   subroutine init_soil_misc(state)
      use swap_state_mod, only: swap_state_t
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none
      type(swap_state_t), intent(inout) :: state

      integer i

      associate (mesh => state%mesh,         &
                 soil => state%soilwater,    &
                 heat => state%heat,         &
                 atmo => state%atmosphere)

         ! Initialize Soilwater rate/state variables


         ! Initialize miscellaneous
         soil%hatm = -2.75e5_real64
      atmo%nraidt = 0.0_real64
      atmo%nird   = 0.0_real64
      if (soil%swinco.ne.3) then
        atmo%ldwet = 0.0_real64
        atmo%spev  = 0.0_real64
        atmo%saev  = 0.0_real64
      endif
      soil%runon = 0.0_real64
      soil%qtop = 0.d0
      do i = 1,mesh%numnod+1
        soil%q(i) = 0.0_real64
      enddo
      soil%evp = 0.0_real64                                  ! state sized mesh%numnod; blanket zero
      if (allocated(heat%rfcp)) heat%rfcp = 1.0d0
      state%surfacewater%vtair = 0.0d0
      ! cQMpLatSs retired-zero write dropped (ADR 0040 macropore).

      end associate
      end subroutine init_soil_misc

   !> Populate per-layer and per-node Van Genuchten hydraulic parameters,
   !> saturated/residual water contents, and hysteresis branch selection.
   subroutine populate_hydraulic_params(state, hyd)
      use swap_state_mod, only: swap_state_t
      use soil_config_mod, only: soil_hydraulics_t
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none
      type(swap_state_t),      intent(inout) :: state
      type(soil_hydraulics_t), intent(in)    :: hyd

      integer lay, node

      associate (mesh => state%mesh,         &
                 soil => state%soilwater)

      ! Soil physics: tabulated or MualemVanGenuchten functions
      ! Mirror mesh%layer into state.
      ! iHWCKmodel mirror loop retired — soilwater_init seeds
      !   `sw%iHWCKmodel = 1` (HACK Phase 4f-extend: TOML pipeline forces uni-modal MvG).
      do node = 1, mesh%numnod
         soil%layer(node) = mesh%layer(node)
      end do
      ! BiModal/NoVap: only set via legacy readswap (not TOML path); stay .false.

      ! swsophy=1 init block retired (TSPACK + tabulated dispatch dormant).
      ! See src/soil/dormant/sptabulated.f90 for the original layer→node + sptab/numtab/
      ! ientrytab population logic. Reactivation requires: layer-keyed `numtablay`/
      ! `ientrytablay`/`sptablay` config wiring, restoring the state fields
      ! `state%soilwater%{numtab,ientrytab,sptab}` + their allocations, and removing
      ! the gate below.
      if (soil%swsophy.eq.1) then
         call fatalerr_collected('SoilHydraulics', &
            'swsophy=1 (tabulated soil hydraulics) is dormant — see src/soil/dormant/sptabulated.f90')
      end if
      ! MvanG functions — populate per-node `vg_params` directly from `hyd` (config) reads.
      ! paramvg index layout (legacy reference):
      !   paramvg(1,lay)=ores, (2)=osat, (3)=ksatfit, (4)=alfa, (5)=lexp, (6)=npar,
      !   (7)=1-1/npar, (8)=alfaw (read by hysteresis only), (9)=h_enpr,
      !   (10)=ksatexm sentinel, (11)=relsatthr (always 0 in TOML), (12)=ksatthr (always 0),
      !   (13..21)=bi-modal/extended params (iHWCKmodel != 1 in TOML never fires).
      ! Populate the per-layer VG store (soil%vg_params_layer)
      !   first. This is the canonical layer-keyed source used both by the per-node
      !   init below and by the tillage mutator (tillage.f90 Change_MvGpars).
      do lay = 1, mesh%numlay
         soil%vg_params_layer(lay)%thetar          = hyd%ores(lay)
         soil%vg_params_layer(lay)%thetas          = hyd%osat(lay)
         soil%vg_params_layer(lay)%ksat            = hyd%ksatfit(lay)
         soil%vg_params_layer(lay)%alpha           = hyd%alfa(lay)
         soil%vg_params_layer(lay)%lpar            = hyd%lexp(lay)
         soil%vg_params_layer(lay)%npar            = hyd%npar(lay)
         soil%vg_params_layer(lay)%mpar            = 1.0_real64 - 1.0_real64/hyd%npar(lay)
         soil%vg_params_layer(lay)%alphaw_sentinel = -9999.9_real64
         soil%vg_params_layer(lay)%h_enpr          = hyd%h_enpr(lay)
         soil%vg_params_layer(lay)%ksatexm         = -999.0_real64
      end do
      do node = 1, mesh%numnod
         lay = mesh%layer(node)
         soil%vg_params(node)%thetar           = hyd%ores(lay)
         soil%vg_params(node)%thetas           = hyd%osat(lay)
         soil%vg_params(node)%ksat             = hyd%ksatfit(lay)
         soil%vg_params(node)%alpha            = hyd%alfa(lay)
         soil%vg_params(node)%lpar             = hyd%lexp(lay)
         soil%vg_params(node)%npar             = hyd%npar(lay)
         soil%vg_params(node)%mpar             = 1.0_real64 - 1.0_real64/hyd%npar(lay)
         soil%vg_params(node)%alphaw_sentinel  = -9999.9_real64
         soil%vg_params(node)%h_enpr           = hyd%h_enpr(lay)
         ! ksatexm: HACK Phase 4f-extend left the legacy ksatexm path unported;
         ! paramvg(10,:) was forced to the -999 sentinel. Match that here.
         soil%vg_params(node)%ksatexm          = -999.0_real64
         ! `relsatthr`/`ksatthr` (legacy paramvg(11,12)) are always 0 in the TOML
         ! pipeline — the threshold-Ksat computation from readswap.f90:802-815 was
         ! not ported. The state fields default to 0, so the mirror lines are dropped.
         if (soil%iHWCKmodel(lay) == 3 .OR. soil%iHWCKmodel(lay) == 6 .OR. soil%iHWCKmodel(lay) == 7 .OR. &
             soil%iHWCKmodel(lay) == 10 .OR. soil%iHWCKmodel(lay) == 11) then
            ! bi-modal MvG (paramvg(13..17)) unreachable in TOML:
            ! config_to_variables forces iHWCKmodel=1. Reactivation needs a config
            ! sub-record for alfa_2/npar_2/mpar_2/omega_1/omega_2.
            call fatalerr_collected('SoilHydraulics', &
               'iHWCKmodel in {3,6,7,10,11} (bi-modal MvG) needs config wiring — see src/state/hydraulic_params_mod.f90')
         end if
         if (soil%iHWCKmodel(lay) == 5 .OR. soil%iHWCKmodel(lay) == 7) then
            ! iHWCKmodel=5/7 (h0 air-entry) needs paramvg(18) config.
            call fatalerr_collected('SoilHydraulics', &
               'iHWCKmodel in {5,7} needs config wiring (h0 air-entry)')
         end if
         if (soil%iHWCKmodel(lay) == 8 .OR. soil%iHWCKmodel(lay) == 9 .OR. &
             soil%iHWCKmodel(lay) == 10 .OR. soil%iHWCKmodel(lay) == 11) then
            ! iHWCKmodel=8/9/10/11 needs paramvg(18..21) config.
            call fatalerr_collected('SoilHydraulics', &
               'iHWCKmodel in {8,9,10,11} needs config wiring (h0/ha/apar/omega_k)')
         end if
      end do
      soil%thetsl = 0.0_real64
      do lay = 1, mesh%numlay
         soil%thetsl(lay) = hyd%osat(lay)
      end do

! --- saturated and residual watercontent of each node; hysteresis parameters
      do node = 1,mesh%numnod
        lay = mesh%layer(node)
        soil%thetar(node) = soil%vg_params(node)%thetar
        soil%thetas(node) = soil%vg_params(node)%thetas
        !!! Kroes: disable combi of swsophy=1 and swhyst=1
        if (soil%swhyst.eq.1) then
           ! Wetting curve
           soil%indeks(node) = 1
           soil%vg_params(node)%alpha = hyd%alfaw(lay)      ! hysteresis: wetting alpha
        elseif (soil%swhyst.eq.0.or.soil%swhyst.eq.2) then
           ! Drying branch or simulation without hysteresis
           soil%indeks(node) = -1
           soil%vg_params(node)%alpha = hyd%alfa(lay)       ! hysteresis: drying alpha
        endif
      end do

      end associate
      end subroutine populate_hydraulic_params

   !> Apply initial conditions per SWINCO: input pressure-head profile,
   !> initial groundwater level, or profile derived from groundwater level.
   subroutine apply_soil_initial_conditions(state)
      use swap_array_dimensions, only: macp, mabbc
      use array_utils, only: afgen
      use swap_state_mod, only: swap_state_t
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none
      type(swap_state_t), intent(inout) :: state

      integer i, j
      real(8) tab(mabbc*2)
      character(len=200) messag

      associate (mesh => state%mesh,         &
                 soil => state%soilwater,    &
                 time => state%timecontrol,  &
                 swbotb => state%soilwater%swbotb_runtime)

      if (soil%swinco.eq.1) then
         ! Pressure head profile is input.
         ! [W3/W4 fix 2026-05-28] Migrated from state%cfg%soil%initial%z_init to
         ! state%soilwater%h_init (typed h_profile_table_t). This block is inside
         ! a `swinco.eq.1` branch but h_init is only populated when `swinco == 3`,
         ! so it is unreachable today — preserved as-is for byte-identical
         ! regression. Investigate / retire in a follow-up arc once the swinco
         ! switch ladder is being refactored.
         if (soil%h_init%is_loaded) then
            do i = 1, size(soil%h_init%rows)
              tab(i*2)   = soil%h(i)
              tab(i*2-1) = abs(soil%h_init%rows(i)%z)
            end do
            do i = 1, mesh%numnod
              soil%h(i) = afgen(tab,macp*2,abs(mesh%z(i)))
            end do
         end if
      endif
      if (soil%swinco.eq.2 .and. swbotb.ne.8) then
        if (abs(soil%gwli-(mesh%z(mesh%numnod)-0.5d0*mesh%dz(mesh%numnod))) .lt.1.0d-4) then
          messag = 'Initial groundwaterlevel (SWINCO=2) is '//          &
     &    'too close to bottom of soil profile'//                       &
     &    ' must be corrected!'
          call fatalerr_collected ('soilwater',messag)
        endif
      endif
      if (soil%swinco.eq.3) then
        ! [W3/W4 fix 2026-05-28] consistency gate now uses state%soilwater%h_init%rows.
        if (soil%h_init%is_loaded) then
           if (size(soil%h_init%rows).ne.mesh%numnod) then
             messag = 'Initial data are read from file (SWINCO=3) and '//  &
       &      'number of nodes/compartments is not consistent with NUMNOD'//&
       &      'must be corrected!'
             call fatalerr_collected ('soilwater',messag)
           endif
        endif
      endif
      if (soil%swinco.eq.1.or.soil%swinco.eq.3) then
         ! Determine groundwater level
         if (soil%h(mesh%numnod) .gt. -1.d-5) then
          i = mesh%numnod
          do while ((soil%h(i) .gt. -1.d-5) .and. (i .gt. 1))
              i = i - 1
          end do
          if (soil%h(i) .lt. -1.d-5) then
            soil%gwl = mesh%z(i+1) + soil%h(i+1) / (soil%h(i+1) - soil%h(i)) * (mesh%z(i) - mesh%z(i+1))
            ! Assume hydrostatic equilibrium in saturated part
            do j = i+1, mesh%numnod
              soil%h(j) = soil%gwl - mesh%z(j)
            end do
          endif
        endif
      else
         ! Pressure head profile is calculated from groundwater level
         if (swbotb.eq.1) then
          soil%gwl = afgen (soil%gwltab,mabbc*2,time%t1900+time%dt-1.d0)

          if(abs(soil%gwl-(mesh%z(mesh%numnod)-0.5d0*mesh%dz(mesh%numnod))) .lt.1.0d-4) then
          messag = 'Groundwaterlevel as bottom boundary (SWBOTB=1) is'//&
     &    'below or to close to bottom of soil profile'//               &
     &    ' must be corrected!'
            call fatalerr_collected ('soilwater',messag)
          endif
        else
          soil%gwl = soil%gwli
        endif
        if (soil%gwl.gt.0.0d0) then
          soil%pond = soil%gwl
        else
          soil%pond = 0.0_real64
        endif
        do i = 1,mesh%numnod
          soil%h(i) = soil%gwl - mesh%z(i)
        end do
      endif

      end associate
      end subroutine apply_soil_initial_conditions

   !> Compute initial water contents, differential moisture capacities, and
   !> hydraulic conductivities (incl. internodal means) for each node.
   subroutine compute_initial_node_hydraulics(state)
      use soilhydraulics_utils, only: watcon, hconduc, moiscap, hcomean
      use swap_state_mod, only: swap_state_t
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none
      type(swap_state_t), intent(inout) :: state

      integer i, node

      associate (mesh => state%mesh,         &
                 soil => state%soilwater,    &
                 heat => state%heat,         &
                 time => state%timecontrol)

      ! In case of preferential flow, adjust Van Genuchten parameters
      do i = 1, mesh%numnod
        soil%theta(i) = watcon(soil%h(i), &
                              soil%vg_params(i), &
                              soil%iHWCKmodel(soil%layer(i)), &
                              i, soil)
      end do

      ! Hydraulic conductivities, differential moisture capacities
      ! and mean hydraulic conductivities for each node
      do node = 1,mesh%numnod
        soil%dimoca(node) = moiscap(soil%h(node), &
                                  soil%vg_params(node), &
                                  soil%iHWCKmodel(soil%layer(node)), &
                                  time%dt, &
                                  node, soil)

        soil%FrArMtrx(node) = 1.0_real64
        soil%k(node) = hconduc(soil%h(node),soil%theta(node),heat%rfcp(node),heat%tsoil(node), &
                             soil%vg_params(node), &
                             soil%iHWCKmodel(soil%layer(node)), &
                             soil%fluseksatexm(node), &
                             node, soil)

        if(node.gt.1) soil%kmean(node) = hcomean(soil%swkmean,soil%k(node-1),soil%k(node),mesh%dz(node-1),mesh%dz(node))
      end do
      soil%kmean(mesh%numnod+1) = soil%k(mesh%numnod)

      end associate
      end subroutine compute_initial_node_hydraulics

   subroutine soilwater_step(state)
      use swap_array_dimensions, only: macp, mabbc, matabentries
      use swap_log, only: log_info, to_str
      use array_utils, only: afgen
      use soilhydraulics_utils, only: watcon, hconduc, moiscap, hcomean
      use soilwaterbalance_mod, only: calcgwl, watstor, integral, fluxes
      use swap_state_mod, only: swap_state_t
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none

      ! Arguments
      type(swap_state_t), intent(inout) :: state

      ! Local variables
      integer lay,node,i,j

      real(8) tab(mabbc*2)
      character(len=200) messag

      associate (mesh => state%mesh,         &
                 soil => state%soilwater,    &
                 drai => state%drainage,     &
                 heat => state%heat,         &
                 atmo => state%atmosphere,   &
                 time => state%timecontrol,  &
                 swbotb => state%soilwater%swbotb_runtime)

         ! Calculate Soilwater rate/state variables

         ! Reset intermediate soil water fluxes
         if (time%flDayStart) then
          call soil%reset_intermediate_per_day()
      end if

      if (time%flZeroIntr) then
        call soil%reset_intermediate()

        soil%IPondBeg = soil%pond
        do node = 1, mesh%numnod
          soil%IThetaBeg(node) = soil%theta(node)
        enddo

        ! macropore(5,...) retired (ADR 0040).
      endif

      ! Reset cumulative soil water fluxes
      if (time%flZeroCumu) then
        call soil%reset_cumulative()

        ! macropore(6,...) retired (ADR 0040).

        ! Reset initial water storage and ponding
        soil%volini = soil%volact
        soil%pondini = soil%pond
      endif

      ! Save state variables of time = t
      call soilwater_save_state(state)

      ! Calculate new soil water state variables
      call headcalc(state)

      end associate

      return
      end subroutine soilwater_step

   subroutine soilwater_update(state)
      use swap_array_dimensions, only: macp, mabbc, matabentries
      use swap_log, only: log_info, to_str
      use array_utils, only: afgen
      use soilhydraulics_utils, only: watcon, hconduc, moiscap, hcomean
      use soilwaterbalance_mod, only: calcgwl, watstor, integral, fluxes
      use swap_state_mod, only: swap_state_t
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none

      ! Arguments
      type(swap_state_t), intent(inout) :: state

      ! Local variables
      integer lay,node,i,j

      real(8) tab(mabbc*2)
      character(len=200) messag

      associate (mesh => state%mesh,         &
                 soil => state%soilwater,    &
                 drai => state%drainage,     &
                 heat => state%heat,         &
                 atmo => state%atmosphere,   &
                 time => state%timecontrol,  &
                 swbotb => state%soilwater%swbotb_runtime)

         ! Update hydraulic conductivities to time level t+1
         do i = 1,mesh%numnod
         soil%k(i) = hconduc(soil%h(i),soil%theta(i), &
                                        heat%rfcp(i),heat%tsoil(i), &
                                        soil%vg_params(i), &
                                        soil%iHWCKmodel(soil%layer(i)), &
                                        soil%fluseksatexm(i), &
                                        i, soil)
         if(i.gt.1)then
            soil%kmean(i) = hcomean(soil%swkmean,soil%k(i-1),soil%k(i),mesh%dz(i-1),mesh%dz(i))
         end if
      enddo
      soil%kmean(mesh%numnod+1) = soil%k(mesh%numnod)

      ! Calculate actual water content of profile
      call watstor (state)

      ! Calculate water fluxes between soil compartments
      call fluxes (state)

      ! macropore(4,...) retired (ADR 0040).

      ! Calculate cumulative fluxes
      call integral (state)

      ! Update parameters for soil water hystereses
      if (soil%swhyst.ne.0) call hysteresis (state)

      end associate

      return
      end subroutine soilwater_update

   !> Save and reset soil water state variables
   !!
   !! Manages state variable storage for time stepping. Can save current state
   !! or reset to previous state (useful for adaptive time stepping).
   !!
   !! @note
   !! Date: January 2007
   !! @endnote
   !!
   subroutine soilwater_save_state(state)
      use swap_state_mod, only: swap_state_t
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none

      type(swap_state_t), intent(inout) :: state

      integer :: i

      associate (mesh => state%mesh, soil => state%soilwater)

         ! Save state variables at time = t.
         do i = 1, mesh%numnod
            soil%hm1(i)    = soil%h(i)
            soil%thetm1(i) = soil%theta(i)
         end do
         soil%gwlm1  = soil%gwl
         soil%pondm1 = soil%pond

      end associate

      return
      end subroutine soilwater_save_state

   subroutine soilwater_restore_state(state)
      use swap_state_mod, only: swap_state_t
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none

      type(swap_state_t), intent(inout) :: state

      integer :: i

      associate (mesh => state%mesh, soil => state%soilwater)

         ! Reset soil state variables.
         do i = 1, mesh%numnod
            soil%h(i)     = soil%hm1(i)
            soil%theta(i) = soil%thetm1(i)
         end do
         soil%kmean(mesh%numnod + 1) = soil%k(mesh%numnod)
         soil%gwl  = soil%gwlm1
         soil%pond = soil%pondm1

      end associate

      return
      end subroutine soilwater_restore_state

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
      use soilhydraulics_utils,  only: moiscap, prhead
      use swap_array_dimensions, only: macp
      use swap_state_mod,        only: swap_state_t
      use hydraulic_params_mod,  only: vanGenuchten_params_t
      implicit none

      type(swap_state_t), intent(inout) :: state

      integer :: node, lay, indtem(macp)
      real(8) :: delp, sew, sed, fvalue
      real(8) :: thetar(macp), thetas(macp), alfamg(macp)
      real(8) :: pvg_thetar, pvg_thetas, pvg_alfa, pvg_npar, pvg_npar_m, pvg_alfa_wet
      type(vanGenuchten_params_t) :: vg_hys

      associate (mesh => state%mesh, soil => state%soilwater, time => state%timecontrol)

         ! Check for hysteretic reversal.
         do node = 1, mesh%numnod
            delp = soil%hm1(node) - soil%h(node)
            if (delp/float(soil%indeks(node)) .gt. soil%tau .and.          &
         &      soil%h(node) .lt. -10.0d0 .and. soil%h(node) .gt. -1.0d3) then
               indtem(node) = -soil%indeks(node)
            else
               indtem(node) =  soil%indeks(node)
            end if
         end do

         ! Adapt parameters along scanning curves.
         do 100 node = 1, mesh%numnod
            lay = mesh%layer(node)

            ! Layer-keyed VG params read from state%soilwater (snapshotted from config at init).
            pvg_thetar   = soil%vg_params_layer(lay)%thetar  ! paramvg(1, lay)
            pvg_thetas   = soil%vg_params_layer(lay)%thetas  ! paramvg(2, lay)
            pvg_alfa     = soil%vg_params_layer(lay)%alpha   ! paramvg(4, lay)  drying alpha
            pvg_npar     = soil%vg_params_layer(lay)%npar    ! paramvg(6, lay)
            pvg_npar_m   = 1.0d0 - 1.0d0/pvg_npar            ! paramvg(7, lay)
            pvg_alfa_wet = soil%alfaw_layer(lay)              ! paramvg(8, lay) when swhyst /= 0

            ! No change.
            if (indtem(node) .eq. soil%indeks(node) .or.              &
         &       abs(pvg_alfa - pvg_alfa_wet) .lt. 1.d-4) goto 100

            ! Relative saturation.
            sew = (1.0d0 + (pvg_alfa_wet*(-soil%h(node)))**pvg_npar)**(-pvg_npar_m)
            sed = (1.0d0 + (pvg_alfa*(-soil%h(node)))**pvg_npar)**(-pvg_npar_m)

            ! Flip the scanning index.
            soil%indeks(node) = -1*soil%indeks(node)

            if (soil%indeks(node) .eq. 1) then
               ! Wetting branch.
               alfamg(node) = pvg_alfa_wet
               thetas(node) = pvg_thetas
               thetar(node) = (soil%theta(node) - thetas(node)*sew)/(1.0d0 - sew)

               fvalue = thetar(node)
               if (thetar(node) .lt. pvg_thetar) thetar(node) = pvg_thetar
               if (thetar(node) .gt. pvg_thetas) thetar(node) = pvg_thetas
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
               alfamg(node) = pvg_alfa
               thetar(node) = pvg_thetar
               thetas(node) = thetar(node) + (soil%theta(node) - thetar(node))/sed

               fvalue = thetas(node)
               if (thetas(node) .lt. pvg_thetar) thetas(node) = pvg_thetar
               if (thetas(node) .gt. pvg_thetas) thetas(node) = pvg_thetas
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