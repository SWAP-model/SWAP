module boundbottom_mod
!> Module for soil profile bottom boundary conditions
!!
!! This module determines the bottom boundary condition for the SWAP
!! hydrological model based on various options including groundwater levels,
!! regional fluxes, seepage, and drainage conditions.
!!
!! @author Original SWAP team
!! @date August 2004 / September 2005
!! @date Modified February 2026 (modularization)
! ----------------------------------------------------------------------
    use variables
    use swap_log, only: log_debug, to_str
    use swap_state_mod, only: swap_state_t
    implicit none

    private
    public :: BoundBottom

contains

    !> Determine soil profile bottom boundary conditions
    !!
    !! This subroutine calculates the bottom boundary condition for the soil
    !! profile based on the selected boundary type (`swbotb`). It handles:
    !!
    !! - Given groundwater levels (swbotb=1)
    !! - Regional bottom flux, prescribed or free drainage (swbotb=±2)
    !! - Seepage/infiltration from deep groundwater (swbotb=3)
    !! - Flux as function of pressure head (swbotb=4)
    !! - Given pressure head at bottom (swbotb=5)
    !! - Zero flux (swbotb=6)
    !! - Free drainage (swbotb=7)
    !! - Lysimeter with free drainage (swbotb=8)
    !!
    !! The subroutine updates global variables `qbot` (bottom flux) and
    !! `qbot_nonfrozen`, and may modify `kmean(numnod+1)` for pressure
    !! head boundaries.
    !!
    !! @note This subroutine operates on global state from the `variables`
    !!       module including `swbotb`, `h`, `numnod`, `qbot`, etc.
    !!
    !! @warning If oven-dry conditions (h < -1.0E7) occur at the bottom node
    !!          with swbotb=2, the boundary automatically switches to free
    !!          drainage (swbotb=-2) and a warning is issued.
    !! @note
    !! File VersionID:
    !!   $Id: boundbottom.f90 362 2018-01-08 13:08:33Z kroes006 $
    !! ----------------------------------------------------------------------
    !! ----------------------------------------------------------------------
    !!     date               : August 2004 / Sept 2005
    !!     purpose            : determine soil profile bottom boundary conditions
    !! ----------------------------------------------------------------------
    !! @endnote
    subroutine BoundBottom(state)
        ! [SS-HEAT] Task 9: state added to access state%heat%rfcp (rfcp global retired)
        use array_utils, only: afgen
        use soilhydraulics_utils, only: watcon, hconduc
        implicit none

        ! [SS-BND B-1.3] intent bumped to inout for dual-write to state%soilwater
        type(swap_state_t), intent(inout) :: state

        ! --- local variables
    integer node, nodnumgwl

        real(8) cvalprof, gwlmean, thetabot, twopi, freq
        real(8) satnodgwl
        character(len=300) messag

        ! ----------------------------------------------------------------------
        twopi = 8.0d0*datan(1.0d0)
        freq = twopi/365.0d0
        ! ----------------------------------------------------------------------
        ! --- interpolation between daily values of given groundwaterlevel
        if (swbotb .eq. 1) then
            gwlinp = afgen(gwltab, mabbc*2, t1900 + dt)
            state%soilwater%gwlinp = gwlinp             ! [SS-BND B-1.3] dual-write
        end if

        ! --- regional bottom flux is given
        if (abs(swbotb) .eq. 2) then

            ! Comment PietG (8-1-08):
            ! ---   if the moisture content in the soil profile is depleted by
            !       a combination of inconsistent boundary conditions, and
            !       the pressure head at the bottom tends to very low values,
            !       then the choice for swbotb=2 is not appropriate.
            if (h(numnod) .lt. -1.0d+7) then ! oven dry conditions at bottom
                if (swbotb .eq. 2) then
                    write (messag, '(a)') 'Oven dry conditions in lowest'
                    write (messag, '(a)') 'compartment therefore switched to'
                    write (messag, '(a)') 'free drainage at date '
                    write (messag, *)
                    write (messag, '(a11)') date
                    call warn('BoundBottom', messag, logf, 0)
                end if
                swbotb = -2
            else
                swbotb = 2
            end if

            if (swbotb .eq. 2) then
                if (sw2 .eq. 1) then
                    ! ---     sine function is used
                    qbot = sinave + sinamp*dcos(freq*(t - sinmax))
                    state%soilwater%qbot = qbot         ! [SS-BND B-1.3] dual-write
                else
                    ! ---     table is used
                    qbot = afgen(qbotab, mabbc*2, t1900 + dt)
                    state%soilwater%qbot = qbot         ! [SS-BND B-1.3] dual-write
                end if
            end if

            ! ---   free drainage assumed in case of h(numnod) < -1.0E7
            if (swbotb .eq. -2) then
                qbot = -1.0d0*kmean(numnod + 1)
                state%soilwater%qbot = qbot             ! [SS-BND B-1.3] dual-write
            end if

        end if

        ! --- seepage or infiltration from/to deep groundwater
        if (swbotb .eq. 3) then
            gwlmean = hdrain + shape*(gwl - hdrain)
            ! ---   determine hydraulic head of deep aquifer
            if (sw3 .eq. 1) then
                deepgw = aqave + aqamp*dcos(twopi/aqper*(t - aqtmax))
                state%soilwater%deepgw = deepgw         ! [SS-BND B-1.3] dual-write
            else
                deepgw = afgen(haqtab, mabbc*2, t1900 + dt)
                state%soilwater%deepgw = deepgw         ! [SS-BND B-1.3] dual-write
            end if

            ! ---   determine C-value (vertical resistance) in saturated part of modelled profile
            if (SwBotb3ResVert .eq. 0) then
                !     -   find number node with groundwater level
                node = numnod
                do while (gwlmean .gt. ztopcp(node) .and. node .gt. 1)
                    node = node - 1
                end do
                nodnumgwl = node
                satnodgwl = gwlmean - zbotcp(nodnumgwl)
                cvalprof = satnodgwl/cofgen(3, nodnumgwl)
                do node = nodnumgwl + 1, numnod
                    cvalprof = cvalprof + dz(node)/cofgen(3, node)
                end do
            elseif (SwBotb3ResVert .eq. 1) then
                cvalprof = 0.0d0
            end if
!
            qbot = (deepgw - gwlmean)/(rimlay + cvalprof)
            state%soilwater%qbot = qbot                 ! [SS-BND B-1.3] dual-write

! ---   extra groundwater flux might be added
            if (sw4 .eq. 1) then
                qbot = qbot + afgen(qbotab, mabbc*2, t1900 + dt)
                state%soilwater%qbot = qbot             ! [SS-BND B-1.3] dual-write
            end if
        end if

! --- flux calculated as function of h
        if (swbotb .eq. 4) then
            if (swqhbot .eq. 1) then
                qbot = cofqha*dexp(cofqhb*dabs(gwl))
                state%soilwater%qbot = qbot             ! [SS-BND B-1.3] dual-write
                if (swcofqhc .eq. 1) then
                    qbot = qbot + cofqhc
                    state%soilwater%qbot = qbot         ! [SS-BND B-1.3] dual-write
                end if
            else if (swqhbot .eq. 2) then
                qbot = afgen(qbotab, mabbc*2, dabs(gwl))
                state%soilwater%qbot = qbot             ! [SS-BND B-1.3] dual-write
            end if
        end if

! --- interpolation between daily values of given pressurehead
        if (swbotb .eq. 5) then
            hbot = afgen(hbotab, mabbc*2, t1900 + dt)
            state%soilwater%hbot = hbot                 ! [SS-BND B-1.3] dual-write
            thetabot = watcon(numnod, hbot)

            kmean(numnod + 1) = hconduc(numnod, hbot, thetabot, state%heat%rfcp(numnod))
            if (flMacroPore) then
                kmean(numnod + 1) = FrArMtrx(numnod)*kmean(numnod + 1)
            end if
        end if

! --- zero flux at the bottom
        if (swbotb .eq. 6) then
            qbot = 0.0d0
            state%soilwater%qbot = qbot                 ! [SS-BND B-1.3] dual-write
        end if

! --- free drainage
        if (swbotb .eq. 7) then
            qbot = -1.0d0*kmean(numnod + 1)
            state%soilwater%qbot = qbot                 ! [SS-BND B-1.3] dual-write
        end if

! --- lysimeter with free drainage
        if (swbotb .eq. 8) then
            qbot = 0.0d0
            state%soilwater%qbot = qbot                 ! [SS-BND B-1.3] dual-write
        end if

        qbot_nonfrozen = qbot
        state%soilwater%qbot_nonfrozen = qbot_nonfrozen ! [SS-BND B-1.3] dual-write

        return
    end subroutine BoundBottom

end module boundbottom_mod
