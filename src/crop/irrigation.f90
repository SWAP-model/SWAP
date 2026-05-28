! File VersionID:
!   $Id: irrigation.f90 372 2018-03-13 10:01:20Z heine003 $
! ----------------------------------------------------------------------
!> Irrigation routines for scheduled and subsurface drip irrigation.
!!
!! This module groups legacy irrigation entry points into a modern
!! module namespace while preserving behavior and routine names.
!!
!! @note Legacy header retained:
!!   Date: November 2004
!!   Purpose: evaluate and schedule irrigations
!!   Modified July 2019:
!!   - tcs=5 obsolete and replaced by tcs=7 (theta) and tcs=8 (presh)
!!   - some calculations only once during initialization
   module irrigation_mod
   use error_mod, only: fatalerr_collected
   use swap_state_mod, only: swap_state_t
   implicit none
   private

   public :: irrigation_step
   public :: ssdi_irrigation_step, ssdi_irrigation_reset

   contains

!> Evaluate and schedule surface irrigation (daily step).
!!
!! Named replacement for the former irrigation(task=2, state) dispatch path.
   subroutine irrigation_step(state)
! ----------------------------------------------------------------------
!     Date               : November 2004
!     Purpose            : evaluate and schedule irrigations
!                        : Modified July 2019;
!                             - tcs=5 obsolete and replaced by tcs=7 (theta) and tcs=8 (presh)
!                             - some calculations only once during initialization
! ----------------------------------------------------------------------
! --  global variables
      use swap_array_dimensions, only: maho
      ! [GR-CROP 2026-05-25] irrigation.f90 is `use variables`-free for the main subroutine.
      ! All previously-imported globals were either:
      !   - migrated to state (state%crop%irrigation, state%atmosphere%isua,
      !     state%crop%common%X), or
      !   - retired-zero (schedule==1 dead branch — see local declarations below).
      ! [state%cfg-retirement cluster 6] irrig_cfg (state%cfg%irrigation) retired;
      !   swirfix → state%crop%irrigation%swirfix; swsolu → state%solute%swsolu.
      use array_utils, only: afgen
      use soilhydraulics_utils, only: watcon
      implicit none

      type(swap_state_t), intent(inout) :: state

! --  local variables
      integer irr,node,nodsen,tcs,tcsfix,dcslim,dcs
      integer ifnd,i,datea(6),irgdayfix
      integer endirr(2),startirr(2)
      integer yearendcrp, yearstacrp
      integer :: irrigevent       ! [GR-CROP 2026-05-25] localized — consumed only within this call
      real(8) frlow,phlo,phhi,phme,awlh,awmh,awah,cdef
      real(8) wclo,wcme,wchi,wcac,tps1,tps2,tps3,tps4,tps5,depl,phcrit
      real(8) dps1,dps2,Tred
      real(8) dvstage(7),trel(7),raw(7),taw(7),dwa(7),hcri(7),tcri(7)
      real(8) di(7),fid(7),irgdepmax,irgdepmin,irgthreshold
      real(4) fsec
      real(8) tstairryrx,tendirryrx, grai_red
      real(8), dimension(maho) :: wclos, wcmes, wchis
      logical flIrriTime
      character(len=80) filnam
      character(len=200) messag

      ! [GR-CROP 2026-05-25] retired-zero — schedule==1 dead branch.
      !   The crop%common%schedule==1 path is gated off in every TOML run
      !   because cropfixed/cropwofost/cropgrass init reject schedule=1 via
      !   fatalerr. These legacy globals were never written by the TOML
      !   pipeline; reading the bare global returned zero. They are
      !   therefore replaced by local zero-initialized values here so the
      !   declarations in variables.f90/legacy_state.f90/initialize.f90 can
      !   retire without changing runtime behaviour.
      integer,   parameter :: swcirrthres = 0
      integer,   parameter :: isuas       = 0
      integer              :: dayfix      = 0
      real(8),   parameter :: cirrs        = 0.0d0
      real(8),   parameter :: cirrthres    = 0.0d0
      real(8),   parameter :: perirrsurp   = 0.0d0
      real(8),   parameter :: raithreshold = 0.0d0
      real(8),   parameter :: tstairrig    = 0.0d0
      real(8),   parameter :: tendirrig    = 0.0d0
      real(8),   parameter :: treltab(14)  = 0.0d0
      real(8),   parameter :: rawtab(14)   = 0.0d0
      real(8),   parameter :: tawtab(14)   = 0.0d0
      real(8),   parameter :: dwatab(14)   = 0.0d0
      real(8),   parameter :: hcritab(14)  = 0.0d0
      real(8),   parameter :: tcritab(14)  = 0.0d0
      real(8),   parameter :: ditab(14)    = 0.0d0
      real(8),   parameter :: fidtab(14)   = 0.0d0

!     dcs(1)  = Amount of under- or over-irrigation (L) in case of a scheduled irrigation event
!     dcs(2)  = Prescribed fixed irrigation depth (L) for each scheduled irrigation event

! ----------------------------------------------------------------------

! ===    determine irrigation rates and states  =============================================
!        daily

      ! [GR-CROP 2026-05-25] sub-record associate aliases.
      associate( &
         crop    => state%crop,             &
         soil    => state%soilwater,        &
         time    => state%timecontrol,      &
         atmo    => state%atmosphere,       &
         solu    => state%solute,           &
         mesh    => state%mesh,             &
         irrig   => state%crop%irrigation   )  ! [state%cfg-retirement cluster 6] irrig_cfg → irrig (state%crop%irrigation)

! ---    reset intermediate soil water fluxes — handled by state%soilwater%reset_intermediate()
         ! igird/inird/cgird/cnird zeroed via state%soilwater%reset_intermediate() in SoilWater(2)
         ! and state%soilwater%reset_cumulative() in same path.

         crop%gird = 0.0d0
         irrigevent = 0

! ---    fixed irrigations events
         if (irrig%swirfix .eq. 1) then
            associate (irr => state%crop%irrigation)
            if (abs(irr%irdate(irr%nirri_fixed) - time%t1900) .lt. 1.d-3) then
               crop%gird = irr%irdepth(irr%nirri_fixed)
               solu%cirr = irr%irconc(irr%nirri_fixed)
               atmo%isua = irr%irtype(irr%nirri_fixed)
               irr%nirri_fixed = irr%nirri_fixed + 1
               irrigevent = 1
            end if
            end associate
         end if

! ---    scheduling mode - current timing and depth criterion

!        scheduled timing within desired period ?
         if (crop%common%schedule.eq.1) then
            !call dtdpar (cropstart(icrop),datea,fsec)
            !yearstacrp = datea(1)
            !call dtdpar (cropend(icrop),datea,fsec)
            !yearendcrp = datea(1)
            flIrriTime = .false.
            if (yearendcrp.gt.yearstacrp) then
               !datea(1) = yearstacrp
               !datea(2) = startirr(2)
               !datea(3) = startirr(1)
               !fsec = 0.0
               !call dtardp (datea, fsec, tstairryrx)
               !datea(1) = yearendcrp
               !datea(2) = endirr(2)
               !datea(3) = endirr(1)
               !call dtardp (datea, fsec, tendirryrx)
               if ( (time%t1900-tstairryrx).gt.1.0d-3 .and. (time%t1900-tendirryrx).le.1.0d-3 ) then
                  flIrriTime = .true.
               end if
            else
               if ((time%t-tstairrig).ge.-1.0d-3.and.(time%t-tendirrig).le.1.0d-3) then
                  flIrriTime = .true.
               end if
            end if
         end if

         if (crop%common%schedule.eq.1 .and. irrigevent.eq.0 .and. crop%common%flCropCalendar &
                 .and. .not. crop%common%flCropHarvest .and. flIrriTime) then
            solu%cirr = cirrs
            atmo%isua = isuas

! ---       determine water holding capacity, readily available water,
! ---       actual available water and water deficit
            frlow = (mesh%ztopcp(crop%common%noddrz) + crop%common%rd) / mesh%dz(crop%common%noddrz)
            awlh = 0.0d0; awmh = 0.0d0; awah = 0.0d0; cdef = 0.0d0
            do node = 1,crop%common%noddrz
               wclo = wclos(mesh%layer(node))*mesh%dz(node);       if (node.eq.crop%common%noddrz) wclo = wclo*frlow
               wcme = wcmes(mesh%layer(node))*mesh%dz(node);       if (node.eq.crop%common%noddrz) wcme = wcme*frlow
               wchi = wchis(mesh%layer(node))*mesh%dz(node);       if (node.eq.crop%common%noddrz) wchi = wchi*frlow
               wcac = watcon(soil%h(node), &
                              soil%vg_params(node), &
                              soil%iHWCKmodel(soil%layer(node)), &
                              node, soil) * mesh%dz(node)
               if (node.eq.crop%common%noddrz) wcac = wcac*frlow
               awlh = awlh+(wclo-wchi)
               awmh = awmh+(wcme-wchi)
               awah = awah+(wcac-wchi)
               cdef = cdef+(wclo-wcac)
            end do

! -1-       timing - allowable daily stress - only under dry stress circumstances
            if (tcs.eq.1) then
               tps1 = afgen(treltab,14,crop%common%dvs)
! ---          transpiration fraction due to drought and salinity stress
               if (soil%iptra_day .gt. 1.d-10) then
                  Tred = 1.0d0 - (soil%iqreddry_day + soil%iqredsol_day) / soil%iptra_day
               else
                  Tred = 1.0d0
               end if
               if (Tred .lt. tps1) irrigevent = 2
            end if

! -2-       timing - depletion of readily available water (fraction)
            if (tcs.eq.2) then
! ---          compare readily available water and actual available water
               tps2 = afgen(rawtab,14,crop%common%dvs)
               depl = tps2*(awlh-awmh)
               if (depl.gt.awlh) depl=awlh
               if (awah .lt. (awlh-depl)) irrigevent = 2
            end if

! -3-       timing - depletion of totally available water (fraction)
            if (tcs.eq.3) then
! ---          compare totally available water and actual available water
               tps3 = afgen(tawtab,14,crop%common%dvs)
               depl = tps3*awlh
               if (awah.lt.(awlh-depl)) irrigevent = 2
            end if

! -4-       timing - allowable amount of depletion
            if (tcs.eq.4) then
! ---          check if depletion amount has been exceeded
               tps4 = afgen(dwatab,14,crop%common%dvs)
               if ((awlh-awah).gt.(tps4*0.1d0)) irrigevent = 2
            end if

! -5-       timing - critical pressure head or moisture content exceeded
!             OBSOLETE: replaced by tcs=7 or tcs=8

! -6-       Timing - fixed irrigation time (weekly during crop growth)
            if (tcs.eq.6) then

!           (weekly) irrigation only when deficit is higher then threshold
!              cdef (cm) en IrgThreshold (mm)
               dayfix = dayfix + 1
               if (dayfix.ge.7)  then
                  dayfix = 0
                  if (10.0d0*cdef.gt.irgthreshold) then
                     irrigevent = 2
                  end if
               end if
            end if

! -7-       timing - critical pressure head at dcrit (node=nodsen) exceeded
            if (tcs.eq.7) then
! ---          calculation of critical pressure head
               tps5 = afgen(hcritab,14,crop%common%dvs)
! PG/JK start  15-feb-2010
! originally not intended to simulate paddy rice fields,
! but made applicable for paddy by changing the statement:
               phcrit = tps5        ! old statement was: phcrit = -abs(tps5)
! PG/JK end    15-feb-2010
! ---          compare critical pressure head and actual pressure head
               if (soil%h(nodsen).le.phcrit) irrigevent = 2
            end if

! -8-       timing - critical watercontent at dcrit (node=nodsen) exceeded
            if (tcs.eq.8) then
               tps5 = afgen(tcritab,14,crop%common%dvs)
! ---          compare critical water content and actual water content
               if (soil%theta(nodsen).le.tps5) irrigevent = 2
               !phcrit = prhead(nodsen,disnod(nodsen),tps5,cofgen,h)
            end if

! -9-       Timing - fixed intervals
            if (tcsfix.eq.1) then
               if (irrigevent.eq.2 .and. (dayfix .ge. irgdayfix)) then
                  irrigevent = 2
                  dayfix = 1
               else
                  irrigevent = 0
                  if (dayfix .lt. irgdayfix) dayfix = dayfix + 1
               end if
            end if

! ---       depth - back to field capacity [cm]
            if ((irrigevent.eq.2).and.(dcs.eq.1)) then
! ---       correct for over- or under irrigation
               dps1 = afgen(ditab,14,crop%common%dvs)
! PG/JK start  15-feb-2010
! option to reduce irrigation on rainy (> raithreshold) day
! raithreshold =     ! threshold (cm/d) to define rainy days;  used to reduce irrigation
               grai_red = 0.0d0
               if (atmo%grai .gt. raithreshold) grai_red = atmo%grai
               crop%gird = max (0.0d0,cdef+dps1*0.1d0-grai_red)
! PG/JK start  15-feb-2010
            end if

! ---       depth - fixed depth [cm]
            if ((irrigevent.eq.2).and.(dcs.eq.2)) then
               dps2 = afgen(fidtab,14,crop%common%dvs)
               crop%gird = dps2*0.1d0
            end if

! ---       depth - limited depth [cm]
            if ((irrigevent.eq.2).and.(dcslim.eq.1)) then
               crop%gird = max(crop%gird,irgdepmin*0.1d0)
               crop%gird = min(crop%gird,irgdepmax*0.1d0)
            end if

! ---       in case of solutes: allow overirrigation when conc exceeds concthreshold
            if (solu%swsolu.eq.1 .and.irrigevent.eq.2 .and.swcirrthres.eq.1) then
               if (solu%cml(nodsen).gt.cirrthres) then
                  crop%gird = crop%gird + 0.01d0*perirrsurp*crop%gird
               end if
            end if

         end if


      end associate

      return
   end subroutine irrigation_step

!> Daily SSDI scheduling and rate assignment (named replacement for case(2)).
subroutine ssdi_irrigation_step(state)

! [GR-CROP 2026-05-25] dt_SSDI_event now lives on state%crop%irrigation —
! cross-file consumer src/core/timecontrol_mod.f90 was migrated in the
! irrigation.f90 sub-arc follow-up. SSDI_irrigation is now use-variables-free.
use swap_state_mod, only: swap_state_t

implicit none
type(swap_state_t), intent(inout) :: state

! local, help
integer                         :: irrigevent   ! [GR-CROP 2026-05-25] localized — consumed only within this call
real(8)                         :: Tred

      ! [GR-CROP 2026-05-25] sub-record associate aliases; SSDI persistent
      ! state now lives on state%crop%irrigation.
      associate( &
         soil => state%soilwater,        &
         time => state%timecontrol,      &
         mesh => state%mesh,             &
         irr  => state%crop%irrigation   )
      irrigevent      = 0
      soil%qssdi(1:mesh%numnod) = 0.0d0
      irr%dt_SSDI_event = 1.0d0
      soil%qssdisum = 0.0d0

      if (irr%ssdi_schedule == 0) then
         ! check if today is a day with ssdi
         if (abs(irr%ssdi_date(irr%nirri) - time%t1900) .lt. 1.d-3) then
            irrigevent                     = 2
            irr%dt_SSDI_event              = irr%ssdi_amount_f(irr%nirri) / irr%ssdi_rate_f(irr%nirri)
            soil%qssdi(irr%nod_ssdi(1):irr%nod_ssdi(2)) = irr%ssdi_rate_f(irr%nirri)
            irr%nirri                      = irr%nirri + 1
            soil%qssdisum = soil%qssdisum + sum(soil%qssdi(irr%nod_ssdi(1):irr%nod_ssdi(2)))
         end if
      else
         ! scheduling based on exceedance of a certain threshold
         if (irr%ssdi_sched_type == 1) then
            ! transpiration fraction due to drought and salinity stress
            if (soil%iptra_day .gt. 1.d-10) then
               Tred = 1.0d0 - (soil%iqreddry_day + soil%iqredsol_day) / soil%iptra_day
            else
               Tred = 1.0d0
            end if
            if (Tred .lt. irr%ssdi_threshold) irrigevent = 2

         else if (irr%ssdi_sched_type == 2) then
            if (soil%h(irr%nod_ssdi_sensor) <= irr%ssdi_threshold) irrigevent = 2

         else if (irr%ssdi_sched_type == 3) then
            if (soil%theta(irr%nod_ssdi_sensor) <= irr%ssdi_threshold) irrigevent = 2

         end if

         if (irr%sw_interval == 1) then
            if (irrigevent == 2 .and. (irr%days_counter >= irr%days_interval)) then
               irrigevent       = 2
               irr%days_counter = 1
            else
               irrigevent       = 0
               if (irr%days_counter < irr%days_interval) irr%days_counter = irr%days_counter + 1
            end if
         end if

         if (irrigevent == 2) then
            irr%dt_SSDI_event              = irr%ssdi_amount / irr%ssdi_appl_rate
            soil%qssdi(irr%nod_ssdi(1):irr%nod_ssdi(2)) = irr%ssdi_appl_rate
            soil%qssdisum = soil%qssdisum + sum(soil%qssdi(irr%nod_ssdi(1):irr%nod_ssdi(2)))
         end if

      end if

      end associate

end subroutine ssdi_irrigation_step

!> Reset SSDI event state (named replacement for case(9)).
subroutine ssdi_irrigation_reset(state)

use swap_state_mod, only: swap_state_t

implicit none
type(swap_state_t), intent(inout) :: state

! local, help
integer                         :: irrigevent   ! consumed only within this call

      ! special: reset scheduled irrigation at end of irrigation event
      irrigevent      = 0
      state%soilwater%qssdi(1:state%mesh%numnod) = 0.0d0
      state%crop%irrigation%dt_SSDI_event = 1.0d0
      state%soilwater%qssdisum = 0.0d0

end subroutine ssdi_irrigation_reset

end module irrigation_mod
