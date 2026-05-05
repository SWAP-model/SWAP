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
   subroutine irrigation(task)
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

! --  local variables
      integer irr,node,nodsen,task,tcs,tcsfix,dcslim,dcs
      integer ifnd,i,datea(6),getun2,irgdayfix
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
      logical rdinqr
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
      
! ---    reset intermediate soil water fluxes
         if (flzerointr) then
            igird = 0.0d0
            inird = 0.0d0
         end if

! ---    reset cumulative soil water fluxes
         if (flzerocumu) then
            cgird = 0.0d0
            cnird = 0.0d0
         end if

         gird = 0.0d0
         irrigevent = 0

! ---    fixed irrigations events
         if (swirfix .eq. 1) then
            if (abs(irdate(nirri) - t1900) .lt. 1.d-3) then
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
               if ( (t1900-tstairryrx).gt.1.0d-3 .and. (t1900-tendirryrx).le.1.0d-3 ) then
                  flIrriTime = .true.
               end if
            else
               if ((t-tstairrig).ge.-1.0d-3.and.(t-tendirrig).le.1.0d-3) then
                  flIrriTime = .true.
               end if
            end if
         end if

         if (schedule.eq.1 .and. irrigevent.eq.0 .and. flCropCalendar .and. .not. flCropHarvest .and. flIrriTime) then
            cirr = cirrs
            isua = isuas

! ---       determine water holding capacity, readily available water, 
! ---       actual available water and water deficit
            frlow = (ztopcp(noddrz) + rd) / dz(noddrz)
            awlh = 0.0d0; awmh = 0.0d0; awah = 0.0d0; cdef = 0.0d0
            do node = 1,noddrz
               wclo = wclos(layer(node))*dz(node);       if (node.eq.noddrz) wclo = wclo*frlow
               wcme = wcmes(layer(node))*dz(node);       if (node.eq.noddrz) wcme = wcme*frlow
               wchi = wchis(layer(node))*dz(node);       if (node.eq.noddrz) wchi = wchi*frlow
               wcac = watcon(node,h(node))*dz(node);     if (node.eq.noddrz) wcac = wcac*frlow
               awlh = awlh+(wclo-wchi)
               awmh = awmh+(wcme-wchi)
               awah = awah+(wcac-wchi)
               cdef = cdef+(wclo-wcac) 
            end do

! -1-       timing - allowable daily stress - only under dry stress circumstances
            if (tcs.eq.1) then
               tps1 = afgen(treltab,14,dvs)
! ---          transpiration fraction due to drought and salinity stress
               if (iptra_day .gt. 1.d-10) then
                  Tred = 1.0d0 - (iqreddry_day + iqredsol_day) / iptra_day
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
               if (h(nodsen).le.phcrit) irrigevent = 2
            end if

! -8-       timing - critical watercontent at dcrit (node=nodsen) exceeded
            if (tcs.eq.8) then
               tps5 = afgen(tcritab,14,dvs)
! ---          compare critical water content and actual water content
               if (theta(nodsen).le.tps5) irrigevent = 2
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
               if (grai .gt. raithreshold) grai_red = grai
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
               if (cml(nodsen).gt.cirrthres) then
                  gird = gird + 0.01d0*perirrsurp*gird
               end if
            end if

         end if

         if (irrigevent .ne. 0) flIrrigationOutput = .true.

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
subroutine SSDI_irrigation(iTask)

use variables, only: swpfile, logf, mairg, numnod, tend, tstart, t1900, zbotcp, irrigevent, qssdi, qssdisum, dt_SSDI_event,   &
                     h, theta, iptra_day, iqreddry_day, iqredsol_day, &
                     swssdi_irr, nod_ssdi_irr, ssdi_schedule_irr, ssdi_sched_type_irr, &
                     nod_ssdi_sensor_irr, ssdi_threshold_irr, ssdi_threshold_z_irr, &
                     ssdi_amount_irr, ssdi_appl_rate_irr, sw_interval_irr, days_interval_irr, &
                     days_counter_irr, nirri_ssdi_irr, ssdi_date_irr, ssdi_rate_f_irr, ssdi_amount_f_irr

implicit none
! global
integer, intent(in) :: iTask

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
integer                         :: i, j, swp, ifnd
real(8)                         :: Tred
real(8), dimension(2)           :: ssdi_z
character(len=132)              :: ssdi_file

! functions
integer :: getun2
logical :: rdinqr

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
      ! Phase 4f-extend SS-B (ADR 0020): entry to this case implies the
      ! call-site gate flSSDI was true → swssdi=1 was set in the TOML
      ! config and copied to the global. Reads SSDI parameters from
      ! staged swap.swp via TTutil; future ADR 0021 will replace this
      ! with a TOML schema port of the SSDI block.
      swp = getun2(10, 90, 2)
      call rdinit(swp, logf, swpfile)
         if (rdinqr('swssdi')) call rdsinr ('swssdi', 0, 1, swssdi)
         if (swssdi == 1)      call rdscha ('ssdi_file', ssdi_file)
      close(swp)
      call read_ssdi_input()
      
      ! check if ssdi_date is ascending, and determine initial entry point nirri
      nirri = 0
      if (ssdi_schedule == 0) then
         nirri = 1
         do i = 1, ifnd-1
            if (ssdi_date(i) >= ssdi_date(i+1)) call fatalerr_collected ('SSDI_irrigation', 'ssdi_date not in ascending order')
            if (t1900 >= ssdi_date(i)) nirri = i
         end do
         if (t1900 >= ssdi_date(ifnd)) nirri = ifnd
      end if

      ! determine layer number of sensor (if applicable)
      nod_ssdi_sensor = 0
      if (ssdi_schedule == 1 .AND. ssdi_sched_type > 1) then
         i = 1
         do while (zbotcp(i) .gt. (ssdi_threshold_z + 1.0d-5))
            i = i + 1
         end do
         nod_ssdi_sensor = i
      end if
      
      ! determine layer number where SSDI takes place
      do j = 1, 2
         i = 1
         do while (zbotcp(i) .gt. (ssdi_z(j) + 1.0d-5))
            i = i + 1
         end do
         nod_ssdi(j) = i
      end do

      ! redefine application rate: uniformly spread over all nodes
      if (ssdi_schedule == 0) then
         ssdi_amount_f = ssdi_amount_f/dble((nod_ssdi(2) - nod_ssdi(1) + 1))
      else
         ssdi_amount = ssdi_amount/dble((nod_ssdi(2) - nod_ssdi(1) + 1))
      end if
      
      ! initialize
      qssdi         = 0.0d0
      dt_SSDI_event = 1.0d0
      
   case (2)
      ! no subsurface drip irrigation: return
      if (swssdi == 0) return
      
      irrigevent      = 0
      qssdi(1:numnod) = 0.0d0
      dt_SSDI_event   = 1.0d0
      qssdisum        = 0.0d0
      
      if (ssdi_schedule == 0) then
         ! check if today is a day with ssdi
         if (abs(ssdi_date(nirri) - t1900) .lt. 1.d-3) then
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
            if (iptra_day .gt. 1.d-10) then
               Tred = 1.0d0 - (iqreddry_day + iqredsol_day) / iptra_day
            else
               Tred = 1.0d0
            end if
            if (Tred .lt. ssdi_threshold) irrigevent = 2
            
         else if (ssdi_sched_type == 2) then
            if (h(nod_ssdi_sensor) <= ssdi_threshold) irrigevent = 2
            
         else if (ssdi_sched_type == 3) then
            if (theta(nod_ssdi_sensor) <= ssdi_threshold) irrigevent = 2

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

   case (9)
      ! special: reset scheduled irrigation at end of irrigation event
      irrigevent      = 0
      qssdi(1:numnod) = 0.0d0
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
   
!---------------------------
   contains

   !> Read SSDI configuration from SSDI input file.
   subroutine read_ssdi_input()
   real(8) :: dummy
   logical :: rdinar
   
   call rdinit(swp, logf, ssdi_file)
      call rdsinr('ssdi_schedule', 0, 1, ssdi_schedule)
      if (ssdi_schedule == 0) then
         call rdatim('ssdi_date',                      ssdi_date,     mairg, ifnd)
         call rdfdor('ssdi_rate_f',    0.0d0, 100.0d0, ssdi_rate_f,   mairg, ifnd)     ! mm/h
         call rdfdor('ssdi_amount_f',  0.0d0, 100.0d0, ssdi_amount_f, mairg, ifnd)     ! mm
         if (.NOT. rdinar('ssdi_z')) then
            ! application at single depth (single compartment)
            call rdsdor('ssdi_z', -100.0d0, 0.0d0, dummy)                ! cm
            ssdi_z(1) = dummy
            ssdi_z(2) = dummy
         else
            ! application over depth interval btween two depth levels (multiple consecutive compartments)
            call rdfdor('ssdi_z', -100.0d0, 0.0d0, ssdi_z, 2, 2)          ! cm
         end if
         ! convert mm/h or mm irrigation to cm/d or cm
         ssdi_rate_f(1:ifnd)   = ssdi_rate_f(1:ifnd)*0.1d0*24.0d0
         ssdi_amount_f(1:ifnd) = ssdi_amount_f(1:ifnd)*0.1d0
      else
         call rdsinr('ssdi_sched_type', 1, 3, ssdi_sched_type)
         if (ssdi_sched_type == 1) call rdsdor('threshold_Tred',     0.0d0,   1.0d0, ssdi_threshold)
         if (ssdi_sched_type == 2) call rdsdor('threshold_presh',   -1.0d7,   0.0d0, ssdi_threshold)
         if (ssdi_sched_type == 3) call rdsdor('threshold_watc',     0.0d0,   1.0d0, ssdi_threshold)
         if (ssdi_sched_type >  1) call rdsdor('threshold_depth', -100.0d0,   0.0d0, ssdi_threshold_z)      ! cm
                                   call rdsdor('ssdi_amount',        0.0d0, 100.0d0, ssdi_amount)           ! mm
                                   call rdsdor('ssdi_appl_rate',     0.0d0, 100.0d0, ssdi_appl_rate)        ! mm/h  !!!
         if (.NOT. rdinar('ssdi_z')) then
            ! application at single depth (single compartment)
            call rdsdor('ssdi_z', -100.0d0, 0.0d0, dummy)                ! cm
            ssdi_z(1) = dummy
            ssdi_z(2) = dummy
         else
            ! application over depth interval btween two depth levels (multiple consecutive compartments)
            call rdfdor('ssdi_z', -100.0d0, 0.0d0, ssdi_z, 2, 2)          ! cm
         end if
         ! Restriction on number of days between two successive SSDI appliocations
         call rdsinr('sw_interval', 0, 1, sw_interval)
         if (sw_interval == 0) then
            days_interval = 1
         else
            call rdsinr('days_interval', 1, 366, days_interval)
         end if
         days_counter = 366

         ! convert mm irrigation to cm
         ssdi_amount = ssdi_amount*0.1d0
         
         ! convert mm/h irrigation to cm/d
         ssdi_appl_rate = ssdi_appl_rate*0.1d0*24.0d0

      end if
   close(swp)
         
!  at least one date must be within simulation period
   if (ssdi_schedule == 0) call checkdate(ifnd, ssdi_date, tend, tstart, 'irdate', 'SSDI_irrigation//swssdi=1')
   
   end subroutine read_ssdi_input
   
end subroutine SSDI_irrigation

end module irrigation_mod
