!> Module for soil-profile bottom boundary conditions.
!!
!! Determines the bottom boundary condition for the SWAP hydrological
!! model: prescribed groundwater level, regional flux, Cauchy seepage
!! to/from deep aquifer, q(h) curves, prescribed pressure head, zero
!! flux, free drainage, or lysimeter.
module boundbottom_mod
   use swap_state_mod,        only: swap_state_t
   use swap_config_mod,       only: swap_config_t
   use swap_log,              only: log_debug, log_warn, to_str
   use swap_array_dimensions, only: mabbc
   implicit none

   private
   public :: BoundBottom

contains

   !> Determine soil-profile bottom boundary condition.
   !!
   !! Branches on `state%soilwater%swbotb_runtime`:
   !!   1    Given groundwater level (interpolated from gwltab)
   !!  ±2    Regional bottom flux (sine or table) — auto-switches to
   !!        free drainage (-2) under oven-dry h(numnod) < -1.0e7
   !!   3    Cauchy seepage/infiltration from deep groundwater
   !!   4    Flux as a function of pressure head (q(h))
   !!   5    Given pressure head at bottom (interpolated from hbotab)
   !!   6    Zero flux
   !!   7    Free drainage
   !!   8    Lysimeter with free drainage
   !!
   !! Outputs (all on state%soilwater): qbot, qbot_nonfrozen, gwlinp,
   !! deepgw, hbot, kmean(numnod+1) (swbotb=5 only).
   subroutine BoundBottom(state, config)
      use array_utils,           only: afgen
      use soilhydraulics_utils,  only: watcon, hconduc
      implicit none

      type(swap_state_t),  intent(inout) :: state
      type(swap_config_t), intent(in)    :: config

      integer :: node, nodnumgwl
      real(8) :: cvalprof, gwlmean, thetabot, twopi, freq
      real(8) :: satnodgwl
      character(len=300) :: messag

      associate (soil => state%soilwater,    &
                 mesh => state%mesh,         &
                 heat => state%heat,         &
                 time => state%timecontrol,  &
                 bb   => config%bottom_boundary)

         twopi = 8.0d0 * datan(1.0d0)
         freq  = twopi / 365.0d0

         ! ---- swbotb=1: interpolated daily groundwater level ----------------
         if (soil%swbotb_runtime .eq. 1) then
            soil%gwlinp = afgen(soil%gwltab, mabbc*2, time%t1900 + time%dt)
         end if

         ! ---- |swbotb|=2: regional bottom flux ------------------------------
         if (abs(soil%swbotb_runtime) .eq. 2) then

            ! PietG (8-1-08): if the moisture content in the profile is
            ! depleted by inconsistent boundary conditions and the bottom
            ! pressure head tends to very low values, swbotb=2 is not
            ! appropriate — fall back to free drainage.
            if (soil%h(mesh%numnod) .lt. -1.0d+7) then
               if (soil%swbotb_runtime .eq. 2) then
                  write (messag, '(a)') 'Oven dry conditions in lowest'
                  write (messag, '(a)') 'compartment therefore switched to'
                  write (messag, '(a)') 'free drainage at date '
                  write (messag, *)
                  write (messag, '(a11)') time%date
                  call log_warn('BoundBottom', messag)
               end if
               soil%swbotb_runtime = -2
            else
               soil%swbotb_runtime = 2
            end if

            if (soil%swbotb_runtime .eq. 2) then
               if (bb%sw2 .eq. 1) then
                  ! Sine wave
                  soil%qbot = bb%sinave + bb%sinamp * dcos(freq * (time%t - bb%sinmax))
               else
                  ! Table
                  soil%qbot = afgen(soil%qbotab, mabbc*2, time%t1900 + time%dt)
               end if
            end if

            ! Free drainage triggered by h(numnod) < -1.0e7
            if (soil%swbotb_runtime .eq. -2) then
               soil%qbot = -1.0d0 * soil%kmean(mesh%numnod + 1)
            end if

         end if

         ! ---- swbotb=3: seepage / infiltration from/to deep groundwater -----
         if (soil%swbotb_runtime .eq. 3) then
            gwlmean = bb%hdrain + bb%shape * (soil%gwl - bb%hdrain)

            ! Hydraulic head of the deep aquifer
            if (bb%sw3 .eq. 1) then
               soil%deepgw = bb%aqave + bb%aqamp * &
                             dcos(twopi / bb%aqper * (time%t - bb%aqtmax))
            else
               soil%deepgw = afgen(soil%haqtab, mabbc*2, time%t1900 + time%dt)
            end if

            ! C-value (vertical resistance) in saturated part of the profile
            if (soil%swbotb3resvert .eq. 0) then
               ! Find the node containing the groundwater level
               node = mesh%numnod
               do while (gwlmean .gt. mesh%ztopcp(node) .and. node .gt. 1)
                  node = node - 1
               end do
               nodnumgwl = node
               satnodgwl = gwlmean - mesh%zbotcp(nodnumgwl)
               cvalprof  = satnodgwl / soil%vg_params(nodnumgwl)%ksat
               do node = nodnumgwl + 1, mesh%numnod
                  cvalprof = cvalprof + mesh%dz(node) / soil%vg_params(node)%ksat
               end do
            elseif (soil%swbotb3resvert .eq. 1) then
               cvalprof = 0.0d0
            end if

            soil%qbot = (soil%deepgw - gwlmean) / (bb%rimlay + cvalprof)

            ! Optional extra groundwater flux from a CSV
            if (bb%sw4 .eq. 1) then
               soil%qbot = soil%qbot + afgen(soil%qbotab, mabbc*2, time%t1900 + time%dt)
            end if
         end if

         ! ---- swbotb=4: flux as a function of pressure head -----------------
         if (soil%swbotb_runtime .eq. 4) then
            if (bb%swqhbot .eq. 1) then
               soil%qbot = bb%cofqha * dexp(bb%cofqhb * dabs(soil%gwl))
               if (bb%swcofqhc .eq. 1) then
                  soil%qbot = soil%qbot + bb%cofqhc
               end if
            else if (bb%swqhbot .eq. 2) then
               soil%qbot = afgen(soil%qbotab, mabbc*2, dabs(soil%gwl))
            end if
         end if

         ! ---- swbotb=5: interpolated daily bottom pressure head -------------
         if (soil%swbotb_runtime .eq. 5) then
            soil%hbot = afgen(soil%hbotab, mabbc*2, time%t1900 + time%dt)
            thetabot  = watcon(soil%hbot,                          &
                               soil%vg_params(mesh%numnod),        &
                               soil%iHWCKmodel(soil%layer(mesh%numnod)), &
                               mesh%numnod, soil)

            soil%kmean(mesh%numnod + 1) =                          &
               hconduc(soil%hbot, thetabot,                        &
                       heat%rfcp(mesh%numnod), heat%tsoil(mesh%numnod), &
                       soil%vg_params(mesh%numnod),                &
                       soil%iHWCKmodel(soil%layer(mesh%numnod)),   &
                       soil%fluseksatexm(mesh%numnod),             &
                       mesh%numnod, soil)
         end if

         ! ---- swbotb=6: zero flux -------------------------------------------
         if (soil%swbotb_runtime .eq. 6) then
            soil%qbot = 0.0d0
         end if

         ! ---- swbotb=7: free drainage ---------------------------------------
         if (soil%swbotb_runtime .eq. 7) then
            soil%qbot = -1.0d0 * soil%kmean(mesh%numnod + 1)
         end if

         ! ---- swbotb=8: lysimeter with free drainage ------------------------
         if (soil%swbotb_runtime .eq. 8) then
            soil%qbot = 0.0d0
         end if

         soil%qbot_nonfrozen = soil%qbot

      end associate

      return
   end subroutine BoundBottom

end module boundbottom_mod
