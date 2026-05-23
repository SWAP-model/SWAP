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
      use swap_state_mod,        only: swap_state_t
      use swap_log,              only: log_debug, to_str
      use surfacewater_utils,    only: runoff
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
      real(8) :: dt   ! [SS-TC TC-14] local copy of state%timecontrol%dt

      ! Hydraulic conductivity for complete frozen soils (constant)
      real(8), parameter :: hconode_vsmall = 1.0d-10

! ----------------------------------------------------------------------
! --- Initialisation
      dt = state%timecontrol%dt   ! [SS-TC TC-14]

! --- runon of present day
      ! [SS-SWC S-2.12B] legacy runon retired — write directly to state%soilwater%runon
      ! [SS-TC TC-14] flDayStart/daycum read via state%timecontrol
      if (state%timecontrol%flDayStart .and. state%soilwater%flrunon) &
         state%soilwater%runon = state%soilwater%runonarr(state%timecontrol%daycum+1)

      state%soilwater%FlRunoff = .false.
      state%soilwater%QMpLatSs = 0.0d0


!     S O I L   E V A P O R A T I O N

! --- Calculate hydraulic conductivity corresponding with hAtm — [SS-SWC S-2.12B] read from state
      if (state%soilwater%hatm.lt.0.0d0) Then
         TheAtm = watcon(dble(state%soilwater%hatm), &
                         state%soilwater%vg_params(1), &
                         state%soilwater%iHWCKmodel(state%soilwater%layer(1)), &
                         1, state%soilwater)                                  ! [SS-SWC S-2.5] [SS-GR-UTILS Task 5]
         ksurf  = hconduc(dble(state%soilwater%hatm),TheAtm,state%heat%rfcp(1),state%heat%tsoil(1), &
                          state%soilwater%vg_params(1), &
                          state%soilwater%iHWCKmodel(state%soilwater%layer(1)), &
                          state%soilwater%fluseksatexm(1), &
                          1, state%soilwater)                                  ! [SS-SWC S-2.5] [SS-GR-UTILS Task 6]
      else

! --- This only occurs if RH is 100% in SWAPS, never used for SWAP
         kSurf = state%soilwater%k(1)                                         ! [SS-SWC S-2.5]
      endif
      k1Atm = hcomean(state%cfg%simulation%numerical%swkmean,kSurf,state%soilwater%k(1),state%mesh%dz(1),state%mesh%dz(1))        ! [SS-SWC S-2.5] [SS-GR-BH B10]

! --- maximum evaporation rate according to Darcy
      Emax = -k1Atm * ((state%soilwater%hatm-state%soilwater%h(1))/state%mesh%disnod(1)+1.0d0)  ! [SS-SWC S-2.5] [SS-GR-BH B10]
      
! --- determine reduced soil evaporation rate
      ! SS-ATM A-2.6: peva/empreva retired — read from state%atmosphere
      ! Config validation enforces swredu ∈ {1, 2}; the legacy "no reduction"
      ! branch (swredu==0) is unreachable in TOML-configured runs.
      if (state%atmosphere%swredu .eq. 0) then
        state%soilwater%reva = min(state%atmosphere%peva, max(0.0d0, Emax))
      else
        state%soilwater%reva = min(state%atmosphere%empreva, max(0.0d0, Emax))
      endif

!     H I G H   A T M O S P H E R I C   D E M A N D
!     flux through ground surface based on precipitation - evaporation 
!     and remaining ponding of previous timestep
      ! [SS-GR-BH B10] ArMpSs assignment deleted (macropore retired ADR 0040 — factor *(1-0)=1 is identity, bit-equivalent)
      ! SS-ATM A-2.6: nraidt/melt retired — read from state%atmosphere
      state%soilwater%q0 = (state%atmosphere%nraidt+state%atmosphere%nird+state%atmosphere%melt) + state%soilwater%runon - state%soilwater%reva  ! [SS-GR-BH B10]
      q1 = - state%soilwater%q0 - state%soilwater%pondm1/dt

!     check whether the atmospheric demand condition applies
      if (q1 .ge. 0.0d0 .and. q1.gt.Emax) then
         state%soilwater%ftoph  = .true.
         state%soilwater%hsurf  = state%soilwater%hatm                       ! [SS-SWC S-2.12B]
         state%soilwater%kmean(1) = k1Atm                            ! [SS-SWC S-2.12B]
         state%soilwater%pond = 0.0d0                                ! [SS-SWC S-2.12B]
         state%soilwater%runots = 0.0d0
         return
      endif

!     maximum conductivity assuming saturation at ground surface (z=0)
      if(state%soilwater%fluseksatexm(1))then                                 ! [SS-SWC S-2.5]
         ks = state%heat%rfcp(1)*state%soilwater%ksatexm(state%mesh%layer(1)) + (1.0d0-state%heat%rfcp(1))*hconode_vsmall  ! [SS-GR-BH B10]
      else
         ks = state%heat%rfcp(1)*state%soilwater%ksatfit(state%mesh%layer(1)) + (1.0d0-state%heat%rfcp(1))*hconode_vsmall  ! [SS-GR-BH B10]
      endif
      state%soilwater%k1max = hcomean(state%cfg%simulation%numerical%swkmean,ks,state%soilwater%k(1),state%mesh%dz(1),state%mesh%dz(1))  ! [SS-GR-BH B10]
!     check whether application of flux=q1 will yield a pressure head >0
!     at ground surface. If not: flux boundary condition is valid
      h0    = state%soilwater%h(1) - state%mesh%disnod(1)*(q1/state%soilwater%k1max+1.0d0)  ! [SS-SWC S-2.5] [SS-GR-BH B10]
      if (h0.le.1.0d-6) then
         state%soilwater%ftoph  = .false.
         state%soilwater%kmean(1) = 0.0d0                            ! [SS-SWC S-2.12B]
         state%soilwater%hsurf  = 0.0d0
         state%soilwater%pond   = 0.0d0                              ! [SS-SWC S-2.12B]
         state%soilwater%runots = 0.0d0
         state%soilwater%qtop   = q1
      else                 ! ponding occurs
         state%soilwater%ftoph  = .true.
         state%soilwater%kmean(1) = state%soilwater%k1max                  ! [SS-SWC S-2.12B] [SS-GR-BH B10]
         state%soilwater%FlRunoff = .true. ! runoff potential possible

! --- calculate max value of pond without runoff
         p1     = state%soilwater%k1max/state%mesh%disnod(1) * dt           ! [SS-GR-BH B10]
         p2     = 1.0d0/(p1+1.0d0)
         state%soilwater%H0max  = p2 * ( state%soilwater%pondm1 + state%soilwater%q0*dt - state%soilwater%k1max*dt + p1*state%soilwater%h(1) )  ! [SS-GR-BH B10]

! [MACRO-RETIRE 2026-05-12] macropore overland-flow branch deleted (ADR 0040).
! Legacy block ran only when FlMacropore=.true. — see legacy/swap-4.2.0.
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
      ! [SS-SWC S-2.12B] pond retired; read/written via state%soilwater%pond
      ! [SS-TC TC-14] dt, t1900 read via state%timecontrol (ADR 0041)
      ! [SS-GR-FINAL B11] mairg → swap_array_dimensions (dimension constant)
      use swap_array_dimensions, only: mairg
      use array_utils,           only: afgen
      use surfacewater_utils,    only: runoff
      use swap_state_mod,        only: swap_state_t
      implicit none

! --- arguments
      type(swap_state_t), intent(inout) :: state

! --- local variables
      INTEGER i
      real(8) h0,h0min,p1,p2
      real(8) :: dt, t1900   ! [SS-TC TC-14] local copies of state%timecontrol fields

! ----------------------------------------------------------------------
      dt    = state%timecontrol%dt
      t1900 = state%timecontrol%t1900

! --  in case of time dependent ponding: determine pondmx
      if (state%surfacewater%swpondmx.eq.1) then
         state%surfacewater%pondmx = afgen (state%surfacewater%pondmxtab, 2*mairg, t1900+dt)
      endif

! --- check whether h0max, the max value of pond, yields a runoff

!      if(swdra.ne.2 .and. h0max.le.pondmx)then

      if(state%soilwater%H0max.le.state%surfacewater%pondmx)then             ! [SS-GR-BH B11]
         state%soilwater%runots = 0.0d0
         state%soilwater%pond  = state%soilwater%H0max                       ! [SS-SWC S-2.12B] [SS-GR-BH B11]
         state%soilwater%hsurf = state%soilwater%pond
         return
      end if

      state%soilwater%runots = runoff(state)
      if(dabs(state%soilwater%runots).lt.1.0d-6)then
!        if no runoff occurs: first estimation of pond is OK
         state%soilwater%pond  = state%soilwater%H0max                       ! [SS-SWC S-2.12B] [SS-GR-BH B11]
         state%soilwater%hsurf = state%soilwater%pond
         return
      else if(dabs(state%soilwater%runots).ge.1.0d-6 .and. state%surfacewater%swdra.ne.2 .and.             &
     &                                dabs(state%surfacewater%rsroexp-1.0d0).lt.1.0d-6)then  ! [SS-GR-BH B11]
         p1 = state%soilwater%k1max/state%mesh%disnod(1) * dt                ! [SS-GR-BH B11]
         p2 = 1.0d0 / (p1 + 1.0d0 + dt/state%surfacewater%rsro)             ! [SS-GR-BH B11]

         state%soilwater%pond     = p2 * ( state%soilwater%pondm1 + state%soilwater%q0*dt - state%soilwater%k1max*dt + p1*state%soilwater%h(1) +  &  ! [SS-SWC S-2.12B] [SS-GR-BH B11]
     &                     dt/state%surfacewater%rsro * state%surfacewater%pondmx )  ! [SS-GR-BH B11]
         state%soilwater%runots = runoff(state)
         state%soilwater%hsurf  = state%soilwater%pond
         return
      else

!        if runoff occurs: find values for pond and runots iteratively

         p1 = state%soilwater%k1max/state%mesh%disnod(1) * dt                ! [SS-GR-BH B11]
         p2 = 1.0d0/(p1+1.0d0)

!        estimation of maximum ponding: ignore runoff
         state%soilwater%H0max = p2 * ( state%soilwater%pondm1 + state%soilwater%q0*dt - state%soilwater%k1max*dt + p1*state%soilwater%h(1) )  ! [SS-SWC S-2.5] [SS-GR-BH B11]
         h0min = 0.0d0
         do i=1,30
            state%soilwater%pond   = 0.5d0 * (state%soilwater%H0max + h0min)  ! [SS-SWC S-2.12B] [SS-GR-BH B11]
            state%soilwater%runots = runoff(state)
            h0     = p2 * ( state%soilwater%pondm1 +state%soilwater%q0*dt -state%soilwater%k1max*dt +p1*state%soilwater%h(1) -state%soilwater%runots)  ! [SS-GR-BH B11]

            if(dabs(state%soilwater%pond-h0).lt.1.0d-6)then
               state%soilwater%hsurf = state%soilwater%pond
               return
            else
               if(h0.gt.state%soilwater%pond)then
                  h0min = state%soilwater%pond
               else
                  state%soilwater%H0max = state%soilwater%pond               ! [SS-GR-BH B11]
               end if
            end if
         end do
      end if

!     if convergence has not been reached: proceed with final value
      state%soilwater%pond   = 0.5d0 * (state%soilwater%H0max + h0min)       ! [SS-SWC S-2.12B] [SS-GR-BH B11]
      state%soilwater%runots = runoff(state)
      state%soilwater%hsurf  = state%soilwater%pond

      return
      end subroutine pondrunoff

end module boundtop_mod