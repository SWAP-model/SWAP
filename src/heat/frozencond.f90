! File VersionID:
!   $Id: frozencond.f90 368 2018-01-11 15:44:15Z heine003 $

!> Module for frozen soil conditions and their effects on hydraulic properties
!!
!! This module handles the simulation of frozen soil conditions in SWAP, including:
!! - Calculation of reduction factors for hydraulic conductivity under frost
!! - Determination of frozen depth (top and bottom boundaries)
!! - Adjustment of drainage and bottom boundary fluxes during frost periods
!!
!! ## Frost Effects on Soil Hydraulics
!!
!! When soil temperatures drop below freezing, ice formation reduces the effective
!! pore space available for water flow. This module implements:
!! - Temperature-dependent reduction factors (rfcp) that scale from 1.0 (unfrozen)
!!   to 0.0 (fully frozen)
!! - Linear interpolation between frost start temperature (tfroststa) and
!!   frost end temperature (tfrostend)
!!
!! ## Frozen Layer Tracking
!!
!! The module identifies and tracks:
!! - **zfrostbot**: Depth of bottom of frozen layer (most negative)
!! - **zfrosttop**: Depth of top of frozen layer (least negative)
!! - **nodfrostbot**: Node number at bottom of frozen zone
!!
!! These are used to modify drainage and boundary fluxes appropriately.
!!
!! @author Original SWAP development team
!! @date Last modified January 2018
module frozencond_mod
  implicit none
  private

  public :: FrozenCond, FrozenBounds

contains

  !> Calculate reduction factors and frozen depth under frozen soil conditions
  !!
  !! Determines the spatial extent and intensity of soil freezing based on simulated
  !! soil temperatures. Calculates:
  !! 1. Reduction factors (rfcp) for each soil node based on temperature
  !! 2. Frozen zone boundaries (top and bottom depths)
  !!
  !! ## Reduction Factor Calculation
  !!
  !! For each soil node, the reduction factor rfcp is calculated as:
  !! - If T ≥ tfroststa: rfcp = 1.0 (no freezing)
  !! - If T ≤ tfrostend: rfcp = 0.0 (fully frozen)
  !! - If tfrostend < T < tfroststa: rfcp = (T - tfrostend) / (tfroststa - tfrostend)
  !!
  !! The reduction factor directly scales hydraulic conductivity, simulating the
  !! blockage of soil pores by ice.
  !!
  !! ## Frozen Depth Determination
  !!
  !! The algorithm searches from bottom to top to find:
  !! 1. **Bottom of frozen layer** (zfrostbot): Deepest point where T ≤ tfrostend
  !! 2. **Top of frozen layer** (zfrosttop): Shallowest point where T ≤ tfrostend
  !!
  !! Linear interpolation between nodes provides sub-node resolution of freeze depths.
  !!
  !! @note
  !! **Original documentation:**
  !! Date: September 2005
  !! Purpose: If soil temperatures are simulated, determine the reduction factors
  !! and frozen depth for frozen conditions
  !! @endnote
  subroutine FrozenCond(state, config)
    use swap_state_mod, only: swap_state_t
    use swap_config_mod, only: swap_config_t
    implicit none

    type(swap_state_t),  intent(inout) :: state
    type(swap_config_t), intent(in)    :: config
    !! Typed simulation state — reads tsoil/tetop from state; writes rfcp/frozen-zone to state only

    ! Local variables
    integer node
    logical flthaw

    associate (heat => state%heat,         &
               mesh => state%mesh,         &
               heat_cfg => config%heat)

      ! Calculate reduction factor for each node.
      do node = 1, mesh%numnod
         heat%rfcp(node) = 1.0d0
         if (config%soil%frost%swfrost .eq. 1) then
            if (heat%tsoil(node) .ge. heat_cfg%tfroststa) then
               heat%rfcp(node) = 1.0d0
            else if (heat%tsoil(node) .le. heat_cfg%tfrostend) then
               heat%rfcp(node) = 0.0d0
            else if (heat%tsoil(node) .lt. heat_cfg%tfroststa .and. &
                     heat%tsoil(node) .gt. heat_cfg%tfrostend) then
               heat%rfcp(node) = (heat%tsoil(node) - heat_cfg%tfrostend) / &
                                 (heat_cfg%tfroststa - heat_cfg%tfrostend)
            end if
         end if
      end do

      ! Determine frozen depth (z) and frozen node number.
      flthaw            = .true.
      heat%nodfrostbot  = -1
      heat%zfrostbot    = 0.0d0
      heat%zfrosttop    = 0.0d0

      ! Search bottom-up for the frozen zone.
      node = mesh%numnod
      do while (flthaw .and. node .gt. 1)
         node = node - 1
         if (heat%tsoil(node) .le. heat_cfg%tfrostend + 1.0d-6) then
            heat%zfrostbot = mesh%z(node + 1) + mesh%disnod(node + 1) * &
                             (heat_cfg%tfrostend - heat%tsoil(node + 1)) / &
                             (heat%tsoil(node) - heat%tsoil(node + 1))
            flthaw            = .false.
            heat%nodfrostbot  = node
         end if
      end do

      ! If a frozen zone is found, search top-down for the upper boundary.
      if (.not. flthaw) then
         flthaw = .true.
         node   = 0
         do while (flthaw .and. node .lt. heat%nodfrostbot)
            node = node + 1
            if (heat%tsoil(node) .le. heat_cfg%tfrostend + 1.0d-6) then
               if (node .eq. 1) then
                  if (heat%tetop .le. heat_cfg%tfrostend) then
                     heat%zfrosttop = 0.0d0
                  else
                     heat%zfrosttop = mesh%z(node) -                        &
                                      (mesh%z(node) - 0.0d0) *              &
                                      (heat%tsoil(node) - heat_cfg%tfrostend) / &
                                      (heat%tsoil(node) - heat%tetop)
                  end if
               else
                  heat%zfrosttop = mesh%z(node) + mesh%disnod(node) *       &
                                   (heat%tsoil(node) - heat_cfg%tfrostend) / &
                                   (heat%tsoil(node) - heat%tsoil(node - 1))
               end if
               heat%zfrosttop = min(0.0d0, heat%zfrosttop)
               flthaw         = .false.
            end if
         end do
      end if

    end associate

    return
  end subroutine FrozenCond

  !> Reduce or stop boundary fluxes under frozen soil conditions
  !!
  !! Modifies drainage and bottom boundary fluxes when soil freezing occurs.
  !! The routine accounts for:
  !! - Available air volume in frozen soil for water storage
  !! - Reduction of hydraulic conductivity in frozen zones
  !! - Redistribution of drainage fluxes among multiple drain levels
  !!
  !! ## Processing Logic
  !!
  !! ### Air Volume Check
  !! Calculates available air-filled pore space in frozen zone:
  !! - volair = Σ(θ_sat - θ) × dz for frozen compartments
  !! - If volair < 0.01 cm, fluxes are severely restricted
  !!
  !! ### Flux Modifications
  !!
  !! **Without drainage systems (swdra=0):**
  !! - If frozen zone exists and volair < 0.01: qbot = 0
  !!
  !! **With drainage systems (swdra>0):**
  !! - Drains below frozen zone (zbotdr < zfrostbot): qdrain = 0
  !! - Hydraulic conductivity scaled by reduction factors (rfcp)
  !! - For very low volair: may set qbot = 0 or redistribute to deepest active drain
  !! - Drainage distribution recalculated using modified K values
  !!
  !! ### Hydraulic Conductivity in Frozen Soil
  !! K_frozen = K_unfrozen × rfcp + (1 - rfcp) × K_min
  !! where K_min = 1.0×10⁻¹⁰ cm/d (near-zero permeability)
  !!
  !! @note
  !! **Original documentation:**
  !! Date: 20070206
  !! Purpose: Reduce or stop boundary (drainage and bottom) fluxes under
  !! frost conditions
  !! @endnote
  subroutine FrozenBounds(state, config)
    use swap_state_mod,        only: swap_state_t
    use swap_config_mod,       only: swap_config_t
    use distribute_drainage,   only: DIVDRA
    use swap_array_dimensions, only: macp
    implicit none

    type(swap_state_t),  intent(inout) :: state
    type(swap_config_t), intent(in)    :: config

    integer :: node, level, layercp(macp), leveldeepest
    real(8) :: volair, ksatcp(macp), cofanicp(macp), qdratot
    real(8) :: zdeepest, ztop
    logical :: frozencomp

    ! Hydraulic conductivity for completely frozen soils.
    real(8), parameter :: hconode_vsmall = 1.0d-10

    ! Initialise qbot from the non-frozen value.
    state%soilwater%qbot = state%soilwater%qbot_nonfrozen

    associate (heat => state%heat,         &
               mesh => state%mesh,         &
               soil => state%soilwater,    &
               surf => state%surfacewater, &
               drai => state%drainage,     &
               time => state%timecontrol)

      ! Available air volume in the frozen zone.
      node       = mesh%numnod
      volair     = 0.0d0
      frozencomp = .true.
      do while (frozencomp)
         volair = volair + (soil%thetas(node) - soil%theta(node))*mesh%dz(node)
         node   = node - 1
         if (node .eq. 0) then
            frozencomp = .false.
         else
            if (heat%rfcp(node) .le. 0.01d0) frozencomp = .false.
         end if
      end do

      ! Apply reduction when air volume is very low. Drainage reduction only
      ! when drainage systems are present.
      if (surf%swdra .eq. 0) then
         if (heat%nodfrostbot .gt. 1 .and. volair .lt. 0.01d0) then
            soil%qbot = 0.0d0
         end if
      else

         if (heat%nodfrostbot .gt. 1 .and. volair .lt. 0.01d0) then

            leveldeepest = 0
            zdeepest     = 0.0d0
            do level = 1, drai%nrlevs
               if (drai%zbotdr(level) .lt. zdeepest) then
                  leveldeepest = level
                  zdeepest     = drai%zbotdr(level)
               end if
            end do

            do node = 1, mesh%numnod
               if (soil%fluseksatexm(node)) then
                  ksatcp(node) = soil%ksatexm(mesh%layer(node))*heat%rfcp(node) + &
                                 (1.0d0 - heat%rfcp(node))*hconode_vsmall
               else
                  ksatcp(node) = soil%ksatfit(mesh%layer(node))*heat%rfcp(node) + &
                                 (1.0d0 - heat%rfcp(node))*hconode_vsmall
               end if

               cofanicp(node) = soil%cofani(mesh%layer(node))
               layercp(node)  = node
               do level = 1, drai%nrlevs
                  if (heat%zfrostbot .lt. drai%zbotdr(level)) then
                     drai%qdra(level, node) = 0.0d0
                     drai%qdrain(level)     = 0.0d0
                  end if
               end do
            end do

            qdratot = 0.0d0
            do level = 1, drai%nrlevs
               qdratot = qdratot + drai%qdrain(level)
            end do

            if (abs(qdratot) .lt. 1.0d-6) then
               if (heat%zfrostbot .lt. drai%zbotdr(leveldeepest)) then
                  soil%qbot = 0.0d0
               else
                  drai%qdrain(leveldeepest) = soil%qbot
               end if
            else
               do level = 1, drai%nrlevs
                  drai%qdrain(level) = drai%qdrain(level) * (1.0d0 + soil%qbot/qdratot)
               end do
            end if

            if (drai%swdivd .eq. 1) then
               ztop = min(soil%gwl, heat%zfrostbot)
               call divdra(mesh%numnod, drai%nrlevs, mesh%dz, ksatcp, ksatcp, soil%fluseksatexm, &
                           layercp, cofanicp, ztop, drai%L, drai%qdrain, drai%qdra,              &
                           drai%swdivdinf, drai%swnrsrf, drai%swtopnrsrf, drai%zbotdr,           &
                           time%dt, drai%facdpthinf, drai%owltab, drai%nowltab, time%t1900)
            end if
         else
            do level = 1, drai%nrlevs
               drai%qdrain(level) = 0.0d0
               do node = 1, mesh%numnod
                  drai%qdra(level, node) = drai%qdra(level, node)*heat%rfcp(node)
                  drai%qdrain(level)     = drai%qdrain(level) + drai%qdra(level, node)
               end do
            end do
         end if

         surf%qdrtot = 0.0d0
         do level = 1, drai%nrlevs
            surf%qdrtot = surf%qdrtot + drai%qdrain(level)
         end do

      end if

    end associate

    return
  end subroutine FrozenBounds

end module frozencond_mod
