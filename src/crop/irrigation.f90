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

   public :: irrigation
   public :: SSDI_irrigation

   contains

!> Evaluate and schedule surface irrigation.
!!
!! @param[in] task Task selector:
!!   - 1: initialization for current crop
!!   - 2: daily irrigation decision and depth
   subroutine irrigation(task, state)
! ----------------------------------------------------------------------
!     Date               : November 2004
!     Purpose            : evaluate and schedule irrigations
!                        : Modified July 2019;
!                             - tcs=5 obsolete and replaced by tcs=7 (theta) and tcs=8 (presh)
!                             - some calculations only once during initialization
! ----------------------------------------------------------------------
! --  global variables
      use variables
      use array_utils, only: afgen
      use soilhydraulics_utils, only: watcon
      implicit none

      type(swap_state_t), intent(in) :: state

! --  local variables
      integer irr,node,nodsen,task,tcs,tcsfix,dcslim,dcs
      integer ifnd,i,datea(6),irgdayfix
      integer endirr(2),startirr(2)
      integer yearendcrp, yearstacrp
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

!     SAVE removed - persistent state now in variables.f90 module (dayfix)
!     save    

!     dcs(1)  = Amount of under- or over-irrigation (L) in case of a scheduled irrigation event
!     dcs(2)  = Prescribed fixed irrigation depth (L) for each scheduled irrigation event

! ----------------------------------------------------------------------
      select case (task)
      case (1)
      ! Legacy fixed-/calculated-irrigation init body deleted as
      ! part of legacy readers physical deletion. The only callers
      ! of irrigation(1) live in src/io/readswap.f90, which itself
      ! has zero callers in src/ or tests/ (parity tests use
      ! literal-value assertions per ADR 0019/SS-A). Modern flow
      ! reaches irrigation(2) only, gated by flIrrigate.
      return

      case (2)

! ===    determine irrigation rates and states  =============================================
!        daily

      ! SS-TC TC-12: t1900, t read via state%timecontrol tc_* aliases.
      associate( &
         tc_t1900 => state%timecontrol%t1900,  &  ! TC-12
         tc_t     => state%timecontrol%t        &  ! TC-12
      )

! ---    reset intermediate soil water fluxes — [SS-SWC S-2.12B] handled by state%soilwater%reset_intermediate()
         ! igird/inird/cgird/cnird zeroed via state%soilwater%reset_intermediate() in SoilWater(2)
         ! and state%soilwater%reset_cumulative() in same path.

         gird = 0.0d0
         irrigevent = 0

! ---    fixed irrigations events
         if (swirfix .eq. 1) then
            if (abs(irdate(nirri) - tc_t1900) .lt. 1.d-3) then  ! TC-12
               gird = irdepth(nirri)
               cirr = irconc(nirri)
               isua = irtype(nirri)
               nirri = nirri + 1
               irrigevent = 1
            end if
         end if

! ---    scheduling mode - current timing and depth criterion

!        scheduled timing within desired period ?
         if (schedule.eq.1) then
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
               if ( (tc_t1900-tstairryrx).gt.1.0d-3 .and. (tc_t1900-tendirryrx).le.1.0d-3 ) then  ! TC-12
                  flIrriTime = .true.
               end if
            else
               if ((tc_t-tstairrig).ge.-1.0d-3.and.(tc_t-tendirrig).le.1.0d-3) then  ! TC-12
                  flIrriTime = .true.
               end if
            end if
         end if

         if (schedule.eq.1 .and. irrigevent.eq.0 .and. flCropCalendar .and. .not. flCropHarvest .and. flIrriTime) then
            cirr = cirrs
            isua = isuas

! ---       determine water holding capacity, readily available water, 
! ---       actual available water and water deficit
            frlow = (state%mesh%ztopcp(noddrz) + rd) / state%mesh%dz(noddrz)  ! [GR-BH C7]
            awlh = 0.0d0; awmh = 0.0d0; awah = 0.0d0; cdef = 0.0d0
            do node = 1,noddrz
               wclo = wclos(state%mesh%layer(node))*state%mesh%dz(node);       if (node.eq.noddrz) wclo = wclo*frlow  ! [GR-BH C7]
               wcme = wcmes(state%mesh%layer(node))*state%mesh%dz(node);       if (node.eq.noddrz) wcme = wcme*frlow  ! [GR-BH C7]
               wchi = wchis(state%mesh%layer(node))*state%mesh%dz(node);       if (node.eq.noddrz) wchi = wchi*frlow  ! [GR-BH C7]
               wcac = watcon(state%soilwater%h(node), &
                              state%soilwater%vg_params(node), &
                              state%soilwater%iHWCKmodel(state%soilwater%layer(node)), &
                              node, state%soilwater) * state%mesh%dz(node)    ! [SS-SWC S-2.12B] [SS-GR-UTILS Task 5] [GR-BH Task 35]
               if (node.eq.noddrz) wcac = wcac*frlow
               awlh = awlh+(wclo-wchi)
               awmh = awmh+(wcme-wchi)
               awah = awah+(wcac-wchi)
               cdef = cdef+(wclo-wcac) 
            end do

! -1-       timing - allowable daily stress - only under dry stress circumstances
            if (tcs.eq.1) then
               tps1 = afgen(treltab,14,dvs)
! ---          transpiration fraction due to drought and salinity stress
               ! [SS-SWC S-2.12B] iptra_day/iqreddry_day/iqredsol_day -> state%soilwater
               if (state%soilwater%iptra_day .gt. 1.d-10) then
                  Tred = 1.0d0 - (state%soilwater%iqreddry_day + state%soilwater%iqredsol_day) / state%soilwater%iptra_day
               else
                  Tred = 1.0d0
               end if
               if (Tred .lt. tps1) irrigevent = 2
            end if

! -2-       timing - depletion of readily available water (fraction)
            if (tcs.eq.2) then
! ---          compare readily available water and actual available water
               tps2 = afgen(rawtab,14,dvs)
               depl = tps2*(awlh-awmh)
               if (depl.gt.awlh) depl=awlh 
               if (awah .lt. (awlh-depl)) irrigevent = 2
            end if

! -3-       timing - depletion of totally available water (fraction)
            if (tcs.eq.3) then
! ---          compare totally available water and actual available water
               tps3 = afgen(tawtab,14,dvs)
               depl = tps3*awlh
               if (awah.lt.(awlh-depl)) irrigevent = 2
            end if

! -4-       timing - allowable amount of depletion
            if (tcs.eq.4) then
! ---          check if depletion amount has been exceeded                
               tps4 = afgen(dwatab,14,dvs)
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
               tps5 = afgen(hcritab,14,dvs)
! PG/JK start  15-feb-2010
! originally not intended to simulate paddy rice fields,
! but made applicable for paddy by changing the statement:
               phcrit = tps5        ! old statement was: phcrit = -abs(tps5)
! PG/JK end    15-feb-2010
! ---          compare critical pressure head and actual pressure head
               if (state%soilwater%h(nodsen).le.phcrit) irrigevent = 2  ! [SS-SWC S-2.12B]
            end if

! -8-       timing - critical watercontent at dcrit (node=nodsen) exceeded
            if (tcs.eq.8) then
               tps5 = afgen(tcritab,14,dvs)
! ---          compare critical water content and actual water content
               if (state%soilwater%theta(nodsen).le.tps5) irrigevent = 2  ! [SS-SWC S-2.12B]
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
               dps1 = afgen(ditab,14,dvs)
! PG/JK start  15-feb-2010
! option to reduce irrigation on rainy (> raithreshold) day
! raithreshold =     ! threshold (cm/d) to define rainy days;  used to reduce irrigation
               grai_red = 0.0d0
               ! SS-ATM A-2.6: grai retired — read from state%atmosphere%grai
               if (state%atmosphere%grai .gt. raithreshold) grai_red = state%atmosphere%grai
               gird = max (0.0d0,cdef+dps1*0.1d0-grai_red) 
! PG/JK start  15-feb-2010
            end if

! ---       depth - fixed depth [cm]
            if ((irrigevent.eq.2).and.(dcs.eq.2)) then
               dps2 = afgen(fidtab,14,dvs)
               gird = dps2*0.1d0
            end if

! ---       depth - limited depth [cm]
            if ((irrigevent.eq.2).and.(dcslim.eq.1)) then
               gird = max(gird,irgdepmin*0.1d0)
               gird = min(gird,irgdepmax*0.1d0)
            end if

! ---       in case of solutes: allow overirrigation when conc exceeds concthreshold
            if (swsolu.eq.1 .and.irrigevent.eq.2 .and.swcirrthres.eq.1) then
               if (state%solute%cml(nodsen).gt.cirrthres) then
                  gird = gird + 0.01d0*perirrsurp*gird
               end if
            end if

         end if

         if (irrigevent .ne. 0) flIrrigationOutput = .true.

      end associate  ! tc_t1900, tc_t (SS-TC TC-12)

      case default
         call fatalerr_collected ('Irrigation', 'Illegal value for TASK')
      end select

      return
   end subroutine irrigation

!> Compute subsurface drip irrigation (SSDI) scheduling and rates.
!!
!! @param[in] iTask Task selector:
!!   - 1: initialization/read SSDI settings
!!   - 2: daily SSDI scheduling and rate assignment
!!   - 9: reset SSDI event state
subroutine SSDI_irrigation(iTask, state)

! [SS-SWC S-2.12B] h/theta/iptra_day/iqreddry_day/iqredsol_day retired — read via state%soilwater
! SS-TC TC-12: t1900 retired from only-list; read via state%timecontrol.
use variables, only: mairg, irrigevent, qssdi, qssdisum, dt_SSDI_event,   &  ! [GR-BH C7] numnod->state%mesh%numnod; zbotcp unused dropped
                     swssdi_irr, nod_ssdi_irr, ssdi_schedule_irr, ssdi_sched_type_irr, &
                     nod_ssdi_sensor_irr, ssdi_threshold_irr, ssdi_threshold_z_irr, &
                     ssdi_amount_irr, ssdi_appl_rate_irr, sw_interval_irr, days_interval_irr, &
                     days_counter_irr, nirri_ssdi_irr, ssdi_date_irr, ssdi_rate_f_irr, ssdi_amount_f_irr
use swap_state_mod, only: swap_state_t

implicit none
! global
integer, intent(in) :: iTask
type(swap_state_t), intent(in) :: state  ! [SS-SWC S-2.12B]

! local aliases for module variables (for minimal code changes)
integer                         :: swssdi
integer, dimension(2)           :: nod_ssdi
integer                         :: ssdi_schedule
integer                         :: ssdi_sched_type
integer                         :: nod_ssdi_sensor
real(8)                         :: ssdi_threshold
real(8)                         :: ssdi_threshold_z
real(8)                         :: ssdi_amount
real(8)                         :: ssdi_appl_rate
integer                         :: sw_interval
integer                         :: days_interval
integer                         :: days_counter
integer                         :: nirri
real(8), dimension(mairg)       :: ssdi_date
real(8), dimension(mairg)       :: ssdi_rate_f
real(8), dimension(mairg)       :: ssdi_amount_f

! local, help
integer                         :: i, j
real(8)                         :: Tred

   ! Load state from module variables at entry
   swssdi = swssdi_irr
   nod_ssdi = nod_ssdi_irr
   ssdi_schedule = ssdi_schedule_irr
   ssdi_sched_type = ssdi_sched_type_irr
   nod_ssdi_sensor = nod_ssdi_sensor_irr
   ssdi_threshold = ssdi_threshold_irr
   ssdi_threshold_z = ssdi_threshold_z_irr
   ssdi_amount = ssdi_amount_irr
   ssdi_appl_rate = ssdi_appl_rate_irr
   sw_interval = sw_interval_irr
   days_interval = days_interval_irr
   days_counter = days_counter_irr
   nirri = nirri_ssdi_irr
   ssdi_date = ssdi_date_irr
   ssdi_rate_f = ssdi_rate_f_irr
   ssdi_amount_f = ssdi_amount_f_irr

   select case  (iTask)
   case (1)
      ! [irrigation.ssdi] init was performed at config-load time by
      ! apply_irrigation_ssdi (config_to_variables.f90). Per ADR 0022,
      ! this case is now a no-op; the runtime reads the staged
      ! _irr snapshots into per-day locals at the top of
      ! SSDI_irrigation (above the select case).
      return
      
   case (2)
      ! SS-TC TC-12: t1900 read via state%timecontrol tc_* alias.
      associate(tc_t1900 => state%timecontrol%t1900)  ! TC-12
      irrigevent      = 0
      qssdi(1:state%mesh%numnod) = 0.0d0  ! [GR-BH C7]
      dt_SSDI_event   = 1.0d0
      qssdisum        = 0.0d0

      if (ssdi_schedule == 0) then
         ! check if today is a day with ssdi
         if (abs(ssdi_date(nirri) - tc_t1900) .lt. 1.d-3) then  ! TC-12
            irrigevent                     = 2
            dt_SSDI_event                  = ssdi_amount_f(nirri) / ssdi_rate_f(nirri)
            qssdi(nod_ssdi(1):nod_ssdi(2)) = ssdi_rate_f(nirri)
            nirri                          = nirri + 1
            qssdisum                       = qssdisum + sum(qssdi(nod_ssdi(1):nod_ssdi(2)))
         end if
      else
         ! scheduling based on exceedance of a certain threshold
         if (ssdi_sched_type == 1) then
            ! transpiration fraction due to drought and salinity stress
            ! [SS-SWC S-2.12B] iptra_day/iqreddry_day/iqredsol_day -> state%soilwater
            if (state%soilwater%iptra_day .gt. 1.d-10) then
               Tred = 1.0d0 - (state%soilwater%iqreddry_day + state%soilwater%iqredsol_day) / state%soilwater%iptra_day
            else
               Tred = 1.0d0
            end if
            if (Tred .lt. ssdi_threshold) irrigevent = 2
            
         else if (ssdi_sched_type == 2) then
            if (state%soilwater%h(nod_ssdi_sensor) <= ssdi_threshold) irrigevent = 2  ! [SS-SWC S-2.12B]
            
         else if (ssdi_sched_type == 3) then
            if (state%soilwater%theta(nod_ssdi_sensor) <= ssdi_threshold) irrigevent = 2  ! [SS-SWC S-2.12B]

         end if
         
         if (sw_interval == 1) then
            if (irrigevent == 2 .and. (days_counter >= days_interval)) then
               irrigevent   = 2
               days_counter = 1
            else
               irrigevent   = 0
               if (days_counter < days_interval) days_counter = days_counter + 1
            end if
         end if

         if (irrigevent == 2) then
            dt_SSDI_event                  = ssdi_amount / ssdi_appl_rate
            qssdi(nod_ssdi(1):nod_ssdi(2)) = ssdi_appl_rate
            qssdisum                       = qssdisum + sum(qssdi(nod_ssdi(1):nod_ssdi(2)))
         end if
         
      end if

      end associate  ! tc_t1900 (SS-TC TC-12)

   case (9)
      ! special: reset scheduled irrigation at end of irrigation event
      irrigevent      = 0
      qssdi(1:state%mesh%numnod) = 0.0d0  ! [GR-BH C7]
      dt_SSDI_event   = 1.0d0
      qssdisum        = 0.0d0
      
   case default
      call fatalerr_collected ('SSDI_irrigation', 'Illegal value for iTask')
   end select
   
   ! Save state back to module variables at exit
   swssdi_irr = swssdi
   nod_ssdi_irr = nod_ssdi
   ssdi_schedule_irr = ssdi_schedule
   ssdi_sched_type_irr = ssdi_sched_type
   nod_ssdi_sensor_irr = nod_ssdi_sensor
   ssdi_threshold_irr = ssdi_threshold
   ssdi_threshold_z_irr = ssdi_threshold_z
   ssdi_amount_irr = ssdi_amount
   ssdi_appl_rate_irr = ssdi_appl_rate
   sw_interval_irr = sw_interval
   days_interval_irr = days_interval
   days_counter_irr = days_counter
   nirri_ssdi_irr = nirri
   ssdi_date_irr = ssdi_date
   ssdi_rate_f_irr = ssdi_rate_f
   ssdi_amount_f_irr = ssdi_amount_f
   
end subroutine SSDI_irrigation

end module irrigation_mod
