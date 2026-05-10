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
  subroutine FrozenCond(state)
    use variables
    use swap_state_mod, only: swap_state_t
    implicit none

    type(swap_state_t), intent(inout) :: state
    !! Typed simulation state — dual-write: writes both legacy globals and state%heat%*

    ! Local variables
    integer node
    logical flthaw

    associate( &
        ht_rfcp        => state%heat%rfcp,        &
        ht_nodfrostbot => state%heat%nodfrostbot,  &
        ht_zfrostbot   => state%heat%zfrostbot,    &
        ht_zfrosttop   => state%heat%zfrosttop)

    ! Calculate reduction factor for each node
    do node=1,numnod
      rfcp(node) = 1.0d0
      if (swfrost.eq.1)then
        if(tsoil(node).ge.tfroststa)then
          rfcp(node) = 1.0d0
        else if(tsoil(node).le.tfrostend) then
          rfcp(node) = 0.0d0
        else if(tsoil(node).lt.tfroststa .and. &
                tsoil(node).gt.tfrostend) then
          rfcp(node) = (tsoil(node)-tfrostend)/ &
                       (tfroststa-tfrostend)
        endif
      endif
    end do
    ht_rfcp(:) = rfcp(1:numnod)                  ! dual-write: mirror to state (rfcp legacy is macp)

    ! Determine frozen depth (z) and frozen node number
    flthaw              = .true.
    nodfrostbot         = -1
    zfrostbot           = 0.0d0
    zfrosttop           = 0.0d0

    ! Search from bottom upward for frozen zone
    node = numnod
    do while (flthaw .and. node.gt.1)
      node = node - 1
      if(tsoil(node) .le. tfrostend+1.0d-6)then
        zfrostbot = z(node+1) + disnod(node+1) * &
                    (tfrostend-tsoil(node+1)) / &
                    (tsoil(node)-tsoil(node+1))
        flthaw             =.false.
        nodfrostbot        = node
      endif
    end do

    ! If frozen zone found, search from top downward for upper boundary
    if(.not.flthaw)then
      flthaw  = .true.
      node = 0
      do while (flthaw .and. node.lt.nodfrostbot)
        node = node + 1
        if(tsoil(node) .le. tfrostend+1.0d-6)then
          if(node.eq.1) then
            if(tetop.le.tfrostend) then
              zfrosttop = 0.0d0
            else
              zfrosttop = z(node) - &
                          (z(node) - 0.0d0) * &
                          (tsoil(node)-tfrostend) / &
                          (tsoil(node)-tetop)
            endif
          else
            zfrosttop = z(node) + disnod(node) * &
                        (tsoil(node)-tfrostend) / &
                        (tsoil(node)-tsoil(node-1))
          endif
          zfrosttop = min(0.0d0,zfrosttop)
          flthaw      =.false.
        endif
      end do
    end if

    ! dual-write: mirror frozen-zone scalars to state
    ht_nodfrostbot = nodfrostbot
    ht_zfrostbot   = zfrostbot
    ht_zfrosttop   = zfrosttop

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
  subroutine FrozenBounds(state)
    use variables
    use distribute_drainage, only: DIVDRA
    use swap_state_mod, only: swap_state_t
    implicit none

    ! SS-DRST Task 3: FrozenBounds modifies qdra/qdrtot after
    ! SurfaceWater has computed them.  Receive state intent(inout) so we can
    ! both read and write state%drainage%qdra / state%surfacewater%qdrtot
    ! (dual-write pattern preserved until drainage Task 7 drops legacy globals).
    type(swap_state_t), intent(inout) :: state

    ! Local variables
    integer node,level,layercp(macp),leveldeepest
    real(8) volair,ksatcp(macp),cofanicp(macp),qdratot
    real(8) zdeepest,ztop
    logical frozencomp

    ! Hydraulic conductivity for completely frozen soils (constant)
    real(8), parameter :: hconode_vsmall = 1.0d-10

    ! Initialize qbot
    qbot = qbot_nonfrozen

    ! Verify available air volume in frozen zone
    node = numnod
    volair = 0.0d0
    frozencomp = .true.
    do while (frozencomp)
      volair = volair + (thetas(node)-theta(node))*dz(node)
      node = node - 1
      if(node.eq.0)then
        frozencomp = .false.
      else
        if(rfcp(node) .le. 0.01d0)then
          frozencomp = .false.
        end if
      end if
    end do

    ! Consider reduction when volair is very low
    ! Reduction of drainage only when systems are present
    if(swdra.eq.0) then
      if(nodfrostbot.gt.1 .and. volair.lt.0.01d0)then
        qbot = 0.0d0
      endif
    else
      ! SS-DRST Phase 2 Task 1: alias state drainage fields; all reads/writes
      ! go directly to state — legacy globals are no longer touched here.
      associate(qdrain => state%drainage%qdrain, qdra => state%drainage%qdra)

      if(nodfrostbot.gt.1 .and. volair.lt.0.01d0)then

        leveldeepest = 0
        zdeepest     = 0.0d0
        do level=1,nrlevs
          if(zbotdr(level).lt.zdeepest) then
            leveldeepest = level
            zdeepest     = zbotdr(level)
          endif
        enddo

        do node=1,numnod
          if(fluseksatexm(node))then
            ksatcp(node)  = ksatexm(layer(node))*rfcp(node) + &
                            (1.0d0-rfcp(node))*hconode_vsmall
          else
            ksatcp(node)  = ksatfit(layer(node))*rfcp(node) + &
                            (1.0d0-rfcp(node))*hconode_vsmall
          endif

          cofanicp(node) = cofani(layer(node))
          layercp(node) = node
          do level=1,nrlevs
            if(zfrostbot.lt.zbotdr(level)) then
              qdra(level,node) = 0.0d0
              qdrain(level) = 0.0d0
            endif
          enddo
        enddo

        qdratot = 0.0d0
        do level = 1,nrlevs
          qdratot = qdratot + qdrain(level)
        end do

        if(abs(qdratot).lt.1.0d-6) then
          if(zfrostbot.lt.zbotdr(leveldeepest)) then
            qbot = 0.0d0
          else
            qdrain(leveldeepest) = qbot
          endif
        else
          do level = 1,nrlevs
            qdrain(level) = qdrain(level) * (1.0d0 + qbot/qdratot)
          end do
        end if

        if (swdivd.eq.1) then
          ztop = min(gwl,zfrostbot)
          call divdra (numnod,nrlevs,dz,ksatcp,ksatcp,fluseksatexm, &
                       layercp,cofanicp,ztop,L,qdrain,qdra, &
                       Swdivdinf,Swnrsrf,SwTopnrsrf,Zbotdr, &
                       dt,FacDpthInf,owltab,t1900)
        endif
      else

        do level = 1,nrlevs
          qdrain(level) = 0.0d0
          do node = 1,numnod
            qdra(level,node) = qdra(level,node)*rfcp(node)
            qdrain(level) = qdrain(level) + qdra(level,node)
          end do
        end do

      endif

      ! SS-SWST Phase 2 Task 11: qdrtot global write dropped; only state written.
      state%surfacewater%qdrtot = 0.0d0
      do level=1,nrlevs
        state%surfacewater%qdrtot = state%surfacewater%qdrtot + qdrain(level)
      end do

      end associate
    endif

    return
  end subroutine FrozenBounds

end module frozencond_mod
