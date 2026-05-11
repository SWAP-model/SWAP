!> Module for soil profile top boundary conditions
!! 
!! This module determines the top (surface) boundary condition for the SWAP
!! hydrological model. It handles atmospheric demand, soil evaporation,
!! precipitation infiltration, surface ponding, runoff generation, and
!! interactions with macropore flow at the soil surface.
!!
!! The module implements a switching boundary condition that can alternate
!! between pressure head (ponding) and flux (atmospheric demand) boundary
!! types depending on soil saturation and infiltration capacity.
!!
!! @author Original SWAP team
!! @date August 2004 - June 2012
!! @date Modified February 2026 (modularization)
!! @note
!! File VersionID:
!!   $Id: boundtop.f90 368 2018-01-11 15:44:15Z heine003 $
!! @endnote
module boundtop_mod
      use variables
      use swap_log, only: log_debug, to_str
      use surfacewater_utils, only: runoff
      use swap_state_mod, only: swap_state_t
      implicit none

      private
      public :: boundtop, PONDRUNOFF

contains

   ! ----------------------------------------------------------------------
   !> Determine soil profile top boundary condition
   !!
   !! This subroutine calculates the surface boundary condition by evaluating
   !! atmospheric demand (precipitation, evaporation) against soil hydraulic
   !! properties to determine whether flux or pressure head conditions apply.
   !!
   !! ## Algorithm Overview
   !!
   !! 1. Calculate soil evaporation based on hydraulic conductivity and
   !!    atmospheric demand (limited by Darcy's law)
   !! 2. Compute net surface flux (precipitation + runon - evaporation)
   !! 3. Check if atmospheric demand condition applies (flux boundary)
   !! 4. If soil cannot accept the flux, switch to ponding (pressure boundary)
   !! 5. Calculate potential macropore infiltration if applicable
   !!
   !! ## Boundary Condition Types
   !!
   !! The subroutine sets either:
   !! - **Flux boundary** (`ftoph = .false.`): When soil can accept atmospheric flux
   !! - **Pressure head boundary** (`ftoph = .true.`): When ponding occurs
   !!
   !! ## Global Variables Modified
   !!
   !! - `ftoph`: Boundary type flag (flux=.false., pressure=.true.)
   !! - `hsurf`: Pressure head at surface [cm]
   !! - `qtop`: Surface flux [cm/d]
   !! - `pond`: Ponding height [cm]
   !! - `runots`: Surface runoff [cm/d]
   !! - `reva`: Actual soil evaporation [cm/d]
   !! - `kmean(1)`: Mean hydraulic conductivity at top boundary [cm/d]
   !! - `QMpLatSs`: Lateral overland flow into macropores at surface [cm/d]
   !!
   !! @note This subroutine operates on global state from the `variables`
   !!       module including atmospheric inputs, soil properties, and
   !!       hydraulic state variables.
   !!
   !! @warning The subroutine may return early if atmospheric demand
   !!          condition applies, leaving `qtop` undefined for ponding cases.
   !! @note
   !! ----------------------------------------------------------------------
   !! ----------------------------------------------------------------------
   !!     date               : August 2004 - June 2012
   !!     purpose            : determine soil profile top boundary condition      
   !! ----------------------------------------------------------------------
   !! @endnote
   subroutine boundtop(state)
   ! [SS-HEAT] Task 9: state added to access state%heat%rfcp (rfcp global retired)
   ! [SS-BND B-2.7] writes only state%soilwater (legacy global writes dropped)
   use soilhydraulics_utils, only: watcon, hconduc, hcomean
   implicit none

   type(swap_state_t), intent(inout) :: state

! --- local variables
      real(8) emax,ks,theatm,ksurf

! ----------------------------------------------------------------------
! --- local variables
      real(8) h0,k1Atm,p1,p2,p2Mp,q1,RsRoMp
      
      ! Hydraulic conductivity for complete frozen soils (constant)
      real(8), parameter :: hconode_vsmall = 1.0d-10

! ----------------------------------------------------------------------
! --- Initialisation

! --- runon of present day
      if (flDayStart .and. flrunon) runon = runonarr(daycum+1)
      state%soilwater%runon = runon                              ! [SS-SWC S-2.4] dual-write

      state%soilwater%FlRunoff = .false.
      state%soilwater%QMpLatSs = 0.0d0


!     S O I L   E V A P O R A T I O N

! --- Calculate hydraulic conductivity corresponding with hAtm
      if (hAtm.lt.0.0d0) Then
         TheAtm = watcon(1,dble(hatm))
         ksurf  = hconduc (1,dble(hatm),TheAtm,state%heat%rfcp(1),state%heat%tsoil(1))
         if(FlMacropore) then
            ksurf = FrArMtrx(1) * ksurf
         endif
      else

! --- This only occurs if RH is 100% in SWAPS, never used for SWAP
         kSurf = k(1)
      endif
      k1Atm = hcomean(swkmean,kSurf,k(1),dz(1),dz(1))

! --- maximum evaporation rate according to Darcy
      Emax = -k1Atm * ((hatm-h(1))/disnod(1)+1.0d0)
      
! --- determine reduced soil evaporation rate
      ! SS-ATM A-2.6: peva/empreva retired — read from state%atmosphere
      if (swredu .eq. 0) then
        state%soilwater%reva = min(state%atmosphere%peva,max(0.0d0,Emax))
      else
        state%soilwater%reva = min(state%atmosphere%empreva,max(0.0d0,Emax))
      endif

!     H I G H   A T M O S P H E R I C   D E M A N D
!     flux through ground surface based on precipitation - evaporation 
!     and remaining ponding of previous timestep
      ArMpSs = 0.d0                                           !     set value of macropore area at soil surface
      if (FlMacropore .and. Z_Tp.gt.-1.d-8) ArMpSs = ArMpTp   
      ! SS-ATM A-2.6: nraidt/melt retired — read from state%atmosphere
      q0 = (state%atmosphere%nraidt+nird+state%atmosphere%melt)*(1.0d0-ArMpSs) + runon - state%soilwater%reva
      q1 = - q0 - pondm1/dt

!     check whether the atmospheric demand condition applies
      if (q1 .ge. 0.0d0 .and. q1.gt.Emax) then
         state%soilwater%ftoph  = .true.
         state%soilwater%hsurf  = hAtm
         kmean(1) = k1Atm
         state%soilwater%kmean(1) = k1Atm
         pond     = 0.0d0
         state%soilwater%pond = 0.0d0
         state%soilwater%runots = 0.0d0
         return
      endif

!     maximum conductivity assuming saturation at ground surface (z=0)
      if(fluseksatexm(1))then
         ks = state%heat%rfcp(1)*ksatexm(1) + (1.0d0-state%heat%rfcp(1))*hconode_vsmall
      else
         ks = state%heat%rfcp(1)*ksatfit(1) + (1.0d0-state%heat%rfcp(1))*hconode_vsmall
      endif
      k1max = hcomean(swkmean,ks,k(1),dz(1),dz(1))
!     check whether application of flux=q1 will yield a pressure head >0 
!     at ground surface. If not: flux boundary condition is valid
      h0    = h(1) - disnod(1)*(q1/k1max+1.0d0)
      if (h0.le.1.0d-6) then
         state%soilwater%ftoph  = .false.
         kmean(1) = 0.0d0
         state%soilwater%kmean(1) = 0.0d0
         state%soilwater%hsurf  = 0.0d0
         pond     = 0.0d0
         state%soilwater%pond   = 0.0d0
         state%soilwater%runots = 0.0d0
         state%soilwater%qtop   = q1
      else                 ! ponding occurs
         state%soilwater%ftoph  = .true.
         kmean(1) = k1max
         state%soilwater%kmean(1) = k1max
         state%soilwater%FlRunoff = .true. ! runoff potential possible

! --- calculate max value of pond without runoff
         p1     = k1max/disnod(1) * dt
         p2     = 1.0d0/(p1+1.0d0)
         h0max  = p2 * ( pondm1 + q0*dt - k1max*dt + p1*h(1) ) 

! --- in case of macropores, calc. potential overland flow into macrop.: QMpLatSs
         if (FlMacropore .and. Z_Tp.gt.-1.d-8) then                     ! Adaptation for GEM 
            if (h0max.gt.PndmxMp) then
               ! SS-ATM A-2.6: nraidt/melt retired — read from state%atmosphere
               RsRoMp  = (h0max + (state%atmosphere%nraidt+nird+state%atmosphere%melt)*ArMpSs*dt) / KsMpSs
               p2Mp    = 1.0d0 / (p1 + 1.0d0 + dt/RsRoMp)
               pond    = (h0max - PndmxMp) * p2Mp/p2
               state%soilwater%pond = pond
               state%soilwater%QMpLatSs = pond * dt/RsRoMp
               state%soilwater%QMpLatSs = dmin1(state%soilwater%QMpLatSs,h0max)
               if (state%soilwater%QMpLatSs.lt.1.0d-7) state%soilwater%QMpLatSs = 0.0d0
            else
               state%soilwater%QMpLatSs = 0.0d0
            endif
         endif
      endif
!  
      return
      end subroutine boundtop


! ----------------------------------------------------------------------
      SUBROUTINE PONDRUNOFF (state)
! ----------------------------------------------------------------------
!     Date               : 4/5/2005
!     Purpose            : determines ponding height and calculates runoff
!     Formal parameters  : state — typed surface-water state (read-only)
!     Subroutines called : -
!     Functions called   : runoff
!     File usage         : -
! ----------------------------------------------------------------------
      use variables, only: swdra,FlMacropore,disnod,dt,h,H0max,k1max,pondm1,pondmx,q0,rsro,rsroexp, &
                           pond,swpondmx,pondmxtab,t1900
      use array_utils, only: afgen
      use surfacewater_utils, only: runoff
      use swap_state_mod, only: swap_state_t
      implicit none

! --- arguments
      type(swap_state_t), intent(inout) :: state

! --- local variables
      INTEGER i
      real(8) h0,h0min,p1,p2
      real(8) q0hlp

! ----------------------------------------------------------------------

! --  in case of time dependent ponding: determine pondmx
      if (swpondmx.eq.1) then
         pondmx = afgen (pondmxtab,2*mairg,t1900+dt)
      endif

! --- in case of Macropores: 
      if (FlMacropore) then
         if (state%soilwater%FlRunoff) then
!   - h0max is reduced with overland flow into Macropores
            q0hlp  = q0 - state%soilwater%QMpLatSs/dt
            p1     = k1max/disnod(1) * dt
            p2     = 1.0d0/(p1+1.0d0)
            h0max  = p2 * ( pondm1 + q0hlp*dt - k1max*dt + p1*h(1) )
            if (h0max.lt.-1.d-9) then
               state%soilwater%QMpLatSs = state%soilwater%QMpLatSs + h0max
               h0max = 0.d0
            endif
         else
!   - inflow excess by direct precipitation into macropores is added to ponding
            if (state%soilwater%QMpLatSs.lt.0.d0) then
               state%soilwater%QMpLatSs = 0.d0
            endif
            return
         endif
      endif
!
! --- check whether h0max, the max value of pond, yields a runoff

!      if(swdra.ne.2 .and. h0max.le.pondmx)then

      if(h0max.le.pondmx)then
         state%soilwater%runots = 0.0d0
         pond     = h0max
         state%soilwater%pond  = pond
         state%soilwater%hsurf = pond
         return
      end if

      state%soilwater%runots = runoff(state)
      if(dabs(state%soilwater%runots).lt.1.0d-6)then
!        if no runoff occurs: first estimation of pond is OK
         pond     = h0max
         state%soilwater%pond  = pond
         state%soilwater%hsurf = pond
         return
      else if(dabs(state%soilwater%runots).ge.1.0d-6 .and. swdra.ne.2 .and.             &
     &                                dabs(rsroexp-1.0d0).lt.1.0d-6)then
         p1 = k1max/disnod(1) * dt
         p2 = 1.0d0 / (p1 + 1.0d0 + dt/rsro)

         pond     = p2 * ( pondm1 + q0*dt - k1max*dt + p1*h(1) +        &
     &                     dt/rsro * pondmx )
         state%soilwater%pond   = pond
         state%soilwater%runots = runoff(state)
         state%soilwater%hsurf  = pond
         return
      else

!        if runoff occurs: find values for pond and runots iteratively

         p1 = k1max/disnod(1) * dt
         p2 = 1.0d0/(p1+1.0d0)

!        estimation of maximum ponding: ignore runoff
         h0max = p2 * ( pondm1 + q0*dt - k1max*dt + p1*h(1) )
         h0min = 0.0d0
         do i=1,30
            pond   = 0.5d0 * (h0max + h0min)
            state%soilwater%pond   = pond
            state%soilwater%runots = runoff(state)
            h0     = p2 * ( pondm1 +q0*dt -k1max*dt +p1*h(1) -state%soilwater%runots)

            if(dabs(pond-h0).lt.1.0d-6)then
               state%soilwater%hsurf = pond
               return
            else
               if(h0.gt.pond)then
                  h0min = pond
               else
                  h0max = pond
               end if
            end if
         end do
      end if

!     if convergence has not been reached: proceed with final value
      pond   = 0.5d0 * (h0max + h0min)
      state%soilwater%pond   = pond
      state%soilwater%runots = runoff(state)
      state%soilwater%hsurf  = pond

      return
      end subroutine pondrunoff

end module boundtop_mod