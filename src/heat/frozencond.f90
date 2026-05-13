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

    associate( &
        ht_tsoil       => state%heat%tsoil,       &
        ht_tetop       => state%heat%tetop,        &
        ht_rfcp        => state%heat%rfcp,         &
        ht_nodfrostbot => state%heat%nodfrostbot,  &
        ht_zfrostbot   => state%heat%zfrostbot,    &
        ht_zfrosttop   => state%heat%zfrosttop)

    ! Calculate reduction factor for each node
    do node=1,state%mesh%numnod
      ht_rfcp(node) = 1.0d0
      if (config%soil%frost%swfrost.eq.1)then
        if(ht_tsoil(node).ge.config%heat%tfroststa)then
          ht_rfcp(node) = 1.0d0
        else if(ht_tsoil(node).le.config%heat%tfrostend) then
          ht_rfcp(node) = 0.0d0
        else if(ht_tsoil(node).lt.config%heat%tfroststa .and. &
                ht_tsoil(node).gt.config%heat%tfrostend) then
          ht_rfcp(node) = (ht_tsoil(node)-config%heat%tfrostend)/ &
                          (config%heat%tfroststa-config%heat%tfrostend)
        endif
      endif
    end do

    ! Determine frozen depth (z) and frozen node number
    flthaw              = .true.
    ht_nodfrostbot      = -1
    ht_zfrostbot        = 0.0d0
    ht_zfrosttop        = 0.0d0

    ! Search from bottom upward for frozen zone
    node = state%mesh%numnod
    do while (flthaw .and. node.gt.1)
      node = node - 1
      if(ht_tsoil(node) .le. config%heat%tfrostend+1.0d-6)then
        ht_zfrostbot = state%mesh%z(node+1) + state%mesh%disnod(node+1) * &
                       (config%heat%tfrostend-ht_tsoil(node+1)) / &
                       (ht_tsoil(node)-ht_tsoil(node+1))
        flthaw           = .false.
        ht_nodfrostbot   = node
      endif
    end do

    ! If frozen zone found, search from top downward for upper boundary
    if(.not.flthaw)then
      flthaw  = .true.
      node = 0
      do while (flthaw .and. node.lt.ht_nodfrostbot)
        node = node + 1
        if(ht_tsoil(node) .le. config%heat%tfrostend+1.0d-6)then
          if(node.eq.1) then
            if(ht_tetop.le.config%heat%tfrostend) then
              ht_zfrosttop = 0.0d0
            else
              ht_zfrosttop = state%mesh%z(node) - &
                             (state%mesh%z(node) - 0.0d0) * &
                             (ht_tsoil(node)-config%heat%tfrostend) / &
                             (ht_tsoil(node)-ht_tetop)
            endif
          else
            ht_zfrosttop = state%mesh%z(node) + state%mesh%disnod(node) * &
                           (ht_tsoil(node)-config%heat%tfrostend) / &
                           (ht_tsoil(node)-ht_tsoil(node-1))
          endif
          ht_zfrosttop = min(0.0d0,ht_zfrosttop)
          flthaw       = .false.
        endif
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

    ! SS-DRST Task 3: FrozenBounds modifies qdra/qdrtot after
    ! SurfaceWater has computed them.  Receive state intent(inout) so we can
    ! both read and write state%drainage%qdra / state%surfacewater%qdrtot
    ! (dual-write pattern preserved until drainage Task 7 drops legacy globals).
    ! SS-HEAT Task 8: rfcp/nodfrostbot/zfrostbot read from state%heat (legacy globals dropped).
    ! GR-BH Task 13: use variables dropped; bare globals replaced with state%* refs.
    type(swap_state_t),  intent(inout) :: state
    type(swap_config_t), intent(in)    :: config

    ! Local variables
    integer node,level,layercp(macp),leveldeepest
    real(8) volair,ksatcp(macp),cofanicp(macp),qdratot
    real(8) zdeepest,ztop
    logical frozencomp

    ! Hydraulic conductivity for completely frozen soils (constant)
    real(8), parameter :: hconode_vsmall = 1.0d-10

    ! Initialize qbot from non-frozen value (SS-BND B-2.7: writes only state%soilwater)
    state%soilwater%qbot = state%soilwater%qbot_nonfrozen

    ! SS-TC TC-11: dt, t1900 read via state%timecontrol tc_* aliases.
    associate( &
        tc_dt            => state%timecontrol%dt,       &  ! TC-11
        tc_t1900         => state%timecontrol%t1900,    &  ! TC-11
        ht_rfcp          => state%heat%rfcp,            &
        ht_nodfrostbot   => state%heat%nodfrostbot,     &
        ht_zfrostbot     => state%heat%zfrostbot,       &
        sw_theta         => state%soilwater%theta,      &  ! [SS-SWC S-2.10]
        sw_thetas        => state%soilwater%thetas,     &  ! [SS-SWC S-2.10]
        sw_fluseksatexm  => state%soilwater%fluseksatexm, & ! [SS-SWC S-2.10]
        sw_gwl           => state%soilwater%gwl,          & ! [SS-SWC S-2.10]
        ms_numnod        => state%mesh%numnod,             & ! GR-BH Task 13
        ms_dz            => state%mesh%dz,                 & ! GR-BH Task 13
        ms_layer         => state%mesh%layer,              & ! GR-BH Task 13
        sw_swdra         => state%surfacewater%swdra,      & ! GR-BH Task 13
        sw_ksatexm       => state%soilwater%ksatexm,       & ! GR-BH Task 13
        sw_ksatfit       => state%soilwater%ksatfit,       & ! GR-BH Task 13
        sw_cofani        => state%soilwater%cofani,        & ! GR-BH Task 13
        dr_nrlevs        => state%drainage%nrlevs,         & ! GR-BH Task 13
        dr_zbotdr        => state%drainage%zbotdr,         & ! GR-BH Task 13
        dr_swdivd        => state%drainage%swdivd,         & ! GR-BH Task 13
        dr_swdivdinf     => state%drainage%swdivdinf,      & ! GR-BH Task 13
        dr_swnrsrf       => state%drainage%swnrsrf,        & ! GR-BH Task 13
        dr_swtopnrsrf    => state%drainage%swtopnrsrf,     & ! GR-BH Task 13
        dr_facdpthinf    => state%drainage%facdpthinf,     & ! GR-BH Task 13
        dr_owltab        => state%drainage%owltab,         & ! GR-BH Task 13
        dr_L             => state%drainage%L)                ! GR-BH Task 13

    ! Verify available air volume in frozen zone
    node = ms_numnod
    volair = 0.0d0
    frozencomp = .true.
    do while (frozencomp)
      volair = volair + (sw_thetas(node)-sw_theta(node))*ms_dz(node)  ! [SS-SWC S-2.10]
      node = node - 1
      if(node.eq.0)then
        frozencomp = .false.
      else
        if(ht_rfcp(node) .le. 0.01d0)then
          frozencomp = .false.
        end if
      end if
    end do

    ! Consider reduction when volair is very low
    ! Reduction of drainage only when systems are present
    if(sw_swdra.eq.0) then
      if(ht_nodfrostbot.gt.1 .and. volair.lt.0.01d0)then
        state%soilwater%qbot = 0.0d0
      endif
    else
      ! SS-DRST Phase 2 Task 1: alias state drainage fields; all reads/writes
      ! go directly to state — legacy globals are no longer touched here.
      associate(qdrain => state%drainage%qdrain, qdra => state%drainage%qdra)

      if(ht_nodfrostbot.gt.1 .and. volair.lt.0.01d0)then

        leveldeepest = 0
        zdeepest     = 0.0d0
        do level=1,dr_nrlevs
          if(dr_zbotdr(level).lt.zdeepest) then
            leveldeepest = level
            zdeepest     = dr_zbotdr(level)
          endif
        enddo

        do node=1,ms_numnod
          if(sw_fluseksatexm(node))then                             ! [SS-SWC S-2.10]
            ksatcp(node)  = sw_ksatexm(ms_layer(node))*ht_rfcp(node) + &
                            (1.0d0-ht_rfcp(node))*hconode_vsmall
          else
            ksatcp(node)  = sw_ksatfit(ms_layer(node))*ht_rfcp(node) + &
                            (1.0d0-ht_rfcp(node))*hconode_vsmall
          endif

          cofanicp(node) = sw_cofani(ms_layer(node))
          layercp(node) = node
          do level=1,dr_nrlevs
            if(ht_zfrostbot.lt.dr_zbotdr(level)) then
              qdra(level,node) = 0.0d0
              qdrain(level) = 0.0d0
            endif
          enddo
        enddo

        qdratot = 0.0d0
        do level = 1,dr_nrlevs
          qdratot = qdratot + qdrain(level)
        end do

        if(abs(qdratot).lt.1.0d-6) then
          if(ht_zfrostbot.lt.dr_zbotdr(leveldeepest)) then
            state%soilwater%qbot = 0.0d0
          else
            ! SS-BND Phase 2 Task B-2.5: read qbot from state (dual-write keeps state current).
            qdrain(leveldeepest) = state%soilwater%qbot
          endif
        else
          do level = 1,dr_nrlevs
            ! SS-BND Phase 2 Task B-2.5: read qbot from state.
            qdrain(level) = qdrain(level) * (1.0d0 + state%soilwater%qbot/qdratot)
          end do
        end if

        if (dr_swdivd.eq.1) then
          ztop = min(sw_gwl,ht_zfrostbot)                          ! [SS-SWC S-2.10]
          call divdra (ms_numnod, dr_nrlevs, ms_dz, ksatcp, ksatcp, sw_fluseksatexm, & ! GR-BH Task 13
                       layercp, cofanicp, ztop, dr_L, qdrain, qdra, &
                       dr_swdivdinf, dr_swnrsrf, dr_swtopnrsrf, dr_zbotdr, &
                       tc_dt, dr_facdpthinf, dr_owltab, tc_t1900)  ! TC-11
        endif
      else

        do level = 1,dr_nrlevs
          qdrain(level) = 0.0d0
          do node = 1,ms_numnod
            qdra(level,node) = qdra(level,node)*ht_rfcp(node)
            qdrain(level) = qdrain(level) + qdra(level,node)
          end do
        end do

      endif

      ! SS-SWST Phase 2 Task 11: qdrtot global write dropped; only state written.
      state%surfacewater%qdrtot = 0.0d0
      do level=1,dr_nrlevs
        state%surfacewater%qdrtot = state%surfacewater%qdrtot + qdrain(level)
      end do

      end associate
    endif

    end associate

    return
  end subroutine FrozenBounds

end module frozencond_mod
