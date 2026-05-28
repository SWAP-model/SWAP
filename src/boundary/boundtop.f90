!> Module for soil profile top boundary conditions.
!!
!! Determines the top (surface) boundary condition for the SWAP
!! hydrological model: atmospheric demand, soil evaporation,
!! precipitation infiltration, surface ponding, and runoff generation.
!!
!! The boundary alternates between pressure-head (ponding) and flux
!! (atmospheric demand) types depending on soil saturation and
!! infiltration capacity.
module boundtop_mod
   use swap_state_mod,        only: swap_state_t
   use swap_log,              only: log_debug, to_str
   use surfacewater_utils,    only: runoff
   implicit none

   private
   public :: boundtop, PONDRUNOFF

contains

   !> Determine soil-profile top boundary condition.
   !!
   !! Algorithm:
   !!   1. Compute soil evaporation from hydraulic conductivity and
   !!      atmospheric demand (limited by Darcy's law).
   !!   2. Form net surface flux (precipitation + runon - evaporation).
   !!   3. If soil can accept the flux, set a flux boundary (ftoph=.false.).
   !!   4. Otherwise switch to a pressure-head boundary (ponding).
   !!
   !! Outputs (all on state%soilwater): ftoph, hsurf, qtop, pond, runots,
   !! reva, kmean(1), QMpLatSs, q0, k1max, H0max, FlRunoff.
   subroutine boundtop(state)
      use soilhydraulics_utils, only: watcon, hconduc, hcomean
      implicit none

      type(swap_state_t), intent(inout) :: state

      real(8) :: emax, ks, theatm, ksurf
      real(8) :: h0, k1Atm, p1, p2, q1
      real(8) :: dt

      ! Hydraulic conductivity for completely frozen soils (constant)
      real(8), parameter :: hconode_vsmall = 1.0d-10

      associate (soil => state%soilwater,    &
                 atmo => state%atmosphere,   &
                 mesh => state%mesh,         &
                 heat => state%heat,         &
                 time => state%timecontrol)

         dt = time%dt

         ! Runon of present day (runonarr is currently dormant — always 0).
         if (time%flDayStart .and. soil%flrunon) &
            soil%runon = soil%runonarr(time%daycum + 1)

         soil%FlRunoff = .false.
         soil%QMpLatSs = 0.0d0

         ! ---- Soil evaporation -----------------------------------------------

         ! Hydraulic conductivity corresponding with hAtm
         if (soil%hatm .lt. 0.0d0) then
            TheAtm = watcon(dble(soil%hatm),                &
                            soil%vg_params(1),              &
                            soil%iHWCKmodel(soil%layer(1)), &
                            1, soil)
            ksurf  = hconduc(dble(soil%hatm), TheAtm,       &
                             heat%rfcp(1), heat%tsoil(1),   &
                             soil%vg_params(1),             &
                             soil%iHWCKmodel(soil%layer(1)),&
                             soil%fluseksatexm(1),          &
                             1, soil)
         else
            ! Only reached if RH is 100% in SWAPS; never used for SWAP.
            kSurf = soil%k(1)
         endif
         k1Atm = hcomean(soil%swkmean, &
                         kSurf, soil%k(1), mesh%dz(1), mesh%dz(1))

         ! Maximum evaporation rate according to Darcy
         Emax = -k1Atm * ((soil%hatm - soil%h(1)) / mesh%disnod(1) + 1.0d0)

         ! Reduced soil evaporation rate. Config validation enforces
         ! swredu ∈ {1, 2}; the legacy swredu==0 "no reduction" branch is
         ! unreachable in TOML-configured runs.
         if (atmo%swredu .eq. 0) then
            soil%reva = min(atmo%peva,    max(0.0d0, Emax))
         else
            soil%reva = min(atmo%empreva, max(0.0d0, Emax))
         endif

         ! ---- High atmospheric demand ---------------------------------------
         ! Net surface flux: precipitation + runon - evaporation, plus
         ! the residual ponding from the previous timestep.
         soil%q0 = (atmo%nraidt + atmo%nird + atmo%melt) + soil%runon - soil%reva
         q1 = -soil%q0 - soil%pondm1 / dt

         ! Atmospheric-demand condition (flux boundary at the atmosphere)
         if (q1 .ge. 0.0d0 .and. q1 .gt. Emax) then
            soil%ftoph    = .true.
            soil%hsurf    = soil%hatm
            soil%kmean(1) = k1Atm
            soil%pond     = 0.0d0
            soil%runots   = 0.0d0
            return
         endif

         ! Maximum conductivity assuming saturation at ground surface (z=0)
         if (soil%fluseksatexm(1)) then
            ks = heat%rfcp(1) * soil%ksatexm(mesh%layer(1)) &
                 + (1.0d0 - heat%rfcp(1)) * hconode_vsmall
         else
            ks = heat%rfcp(1) * soil%ksatfit(mesh%layer(1)) &
                 + (1.0d0 - heat%rfcp(1)) * hconode_vsmall
         endif
         soil%k1max = hcomean(soil%swkmean, &
                              ks, soil%k(1), mesh%dz(1), mesh%dz(1))

         ! Will applying flux=q1 yield a pressure head >0 at the surface?
         ! If not, the flux boundary is valid.
         h0 = soil%h(1) - mesh%disnod(1) * (q1 / soil%k1max + 1.0d0)
         if (h0 .le. 1.0d-6) then
            soil%ftoph    = .false.
            soil%kmean(1) = 0.0d0
            soil%hsurf    = 0.0d0
            soil%pond     = 0.0d0
            soil%runots   = 0.0d0
            soil%qtop     = q1
         else
            ! Ponding occurs
            soil%ftoph    = .true.
            soil%kmean(1) = soil%k1max
            soil%FlRunoff = .true.

            ! Maximum value of pond without runoff
            p1 = soil%k1max / mesh%disnod(1) * dt
            p2 = 1.0d0 / (p1 + 1.0d0)
            soil%H0max = p2 * (soil%pondm1 + soil%q0 * dt        &
                               - soil%k1max * dt + p1 * soil%h(1))
            ! Macropore overland-flow branch deleted per ADR 0040 — see
            ! legacy/swap-4.2.0 for the FlMacropore=.true. block.
         endif

      end associate

      return
   end subroutine boundtop


   !> Determine ponding height and calculate runoff.
   !!
   !! Inputs:  state%soilwater%H0max (preliminary pond), q0, k1max, pondm1, h(1);
   !!          state%surfacewater%pondmx, rsro, rsroexp, swdra; state%mesh%disnod(1);
   !!          state%timecontrol%dt, t1900.
   !! Outputs: state%soilwater%pond, runots, hsurf, H0max.
   subroutine PONDRUNOFF(state)
      use swap_array_dimensions, only: MAIRG
      use array_utils,           only: afgen
      implicit none

      type(swap_state_t), intent(inout) :: state

      integer :: i
      real(8) :: h0, h0min, p1, p2
      real(8) :: dt, t1900

      associate (soil => state%soilwater,    &
                 surf => state%surfacewater, &
                 mesh => state%mesh,         &
                 time => state%timecontrol)

         dt    = time%dt
         t1900 = time%t1900

         ! Time-dependent ponding (dormant — swpondmx is always 0 in TOML).
         if (surf%swpondmx .eq. 1) then
            surf%pondmx = afgen(surf%pondmxtab, 2 * MAIRG, t1900 + dt)
         endif

         ! H0max ≤ pondmx ⇒ no runoff
         if (soil%H0max .le. surf%pondmx) then
            soil%runots = 0.0d0
            soil%pond   = soil%H0max
            soil%hsurf  = soil%pond
            return
         end if

         soil%runots = runoff(state)
         if (dabs(soil%runots) .lt. 1.0d-6) then
            ! No runoff: the first estimation of pond is OK
            soil%pond  = soil%H0max
            soil%hsurf = soil%pond
            return
         else if (dabs(soil%runots) .ge. 1.0d-6 .and. surf%swdra .ne. 2 .and. &
                  dabs(surf%rsroexp - 1.0d0) .lt. 1.0d-6) then
            p1 = soil%k1max / mesh%disnod(1) * dt
            p2 = 1.0d0 / (p1 + 1.0d0 + dt / surf%rsro)

            soil%pond   = p2 * (soil%pondm1 + soil%q0 * dt - soil%k1max * dt   &
                                + p1 * soil%h(1) + dt / surf%rsro * surf%pondmx)
            soil%runots = runoff(state)
            soil%hsurf  = soil%pond
            return
         else
            ! Runoff occurs: iterate to find consistent pond and runots
            p1 = soil%k1max / mesh%disnod(1) * dt
            p2 = 1.0d0 / (p1 + 1.0d0)

            ! Initial estimate of maximum ponding (ignores runoff)
            soil%H0max = p2 * (soil%pondm1 + soil%q0 * dt - soil%k1max * dt &
                               + p1 * soil%h(1))
            h0min = 0.0d0
            do i = 1, 30
               soil%pond   = 0.5d0 * (soil%H0max + h0min)
               soil%runots = runoff(state)
               h0          = p2 * (soil%pondm1 + soil%q0 * dt - soil%k1max * dt &
                                   + p1 * soil%h(1) - soil%runots)

               if (dabs(soil%pond - h0) .lt. 1.0d-6) then
                  soil%hsurf = soil%pond
                  return
               else
                  if (h0 .gt. soil%pond) then
                     h0min = soil%pond
                  else
                     soil%H0max = soil%pond
                  end if
               end if
            end do
         end if

         ! Convergence not reached: proceed with the bisection midpoint
         soil%pond   = 0.5d0 * (soil%H0max + h0min)
         soil%runots = runoff(state)
         soil%hsurf  = soil%pond

      end associate

      return
   end subroutine PONDRUNOFF

end module boundtop_mod
