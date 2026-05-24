! to do: default MvG parameters that are read: do they refer to BDENS or Rho_cons; or should BDENS and RhoCons be equal?
! to do: check if Rho_match differs from BDENS or Rho_cons (if not: division by zero possible)

module tillage_mod
   use error_mod, only: fatalerr_collected
   use swap_state_mod, only: swap_state_t

   ! [GR-CROP 2026-05-25] swsolu/SwDiscrvert reads cut over to direct
   ! config reads via cfg_soil/cfg_solute aliases. Declarations remain
   ! alive in variables.f90 because other consumers exist:
   !   swsolu — still consumed by src/crop/irrigation.f90,
   !             src/core/timecontrol_mod.f90
   !   SwDiscrvert — still consumed by src/soil/dormant/regrid.f90
   ! ParamVG is the last bare global; retired in the next commit.
   use variables, only: ParamVG

   implicit none

!  by default: all in this module is private (local)
   private
!  except for these public routines/functions
   public :: DoTillage

   contains

   subroutine DoTillage (iTask, state)
   ! global
   integer, intent(in)                          :: iTask                               ! Task
   type(swap_state_t), intent(inout)            :: state
   ! local (not to be saved)
   integer                                   :: i
   character(len=20)                         :: STRNG
   logical                                   :: fine
   logical, parameter                        :: TEST = .false.
   logical, parameter                        :: TEST2 = .false.

   ! Sub-record aliases (canonical associate pattern).
   associate( &
      mesh       => state%mesh,             &
      soil       => state%soilwater,        &
      time       => state%timecontrol,      &
      atmo       => state%atmosphere,       &
      tl         => state%tillage,          &
      cfg_soil   => state%cfg%soil,         &
      cfg_solute => state%cfg%solute        )

   if (iTask > 1 .and. cfg_soil%swtill == 0) return      ! no tillage to be considered: return immediately

   ! handle iTask
   select case (iTask)
   case (1)
      ! INITIALIZE

      ! [SS-BMI2] allocate tillage output buffer (builder always runs when swtill=1)
      call init_tillage_output_buffer(state)

      ! some checks: some combinations not (yet) allowed
      if (cfg_soil%swtill == 1) then
         if (cfg_soil%swhyst == 1)              call fatalerr_collected ('DoTillage', 'swhyst = 1 not allowed')
         if (cfg_solute%swsolu == 1)            call fatalerr_collected ('DoTillage', 'swsolu = 1 not (yet) allowed')
         if (state%crop%common%swoxygen == 2)   call fatalerr_collected ('DoTillage', 'swoxygen = 2 not (yet) allowed')
         if (soil%flksatexm)                    call fatalerr_collected ('DoTillage', 'flksatexm not (yet) allowed')
         if (cfg_soil%discretization%swdiscrvert == 1) &
                                                call fatalerr_collected ('DoTillage', 'SwDiscrvert = 1 not (yet) allowed')
      end if

      ! currently: require all Z_tillage = Max_Z_tillage
      do i = 1, tl%Ntill
!         if (dabs(tl%Max_Z_tillage - tl%Z_tillage(i)) > 1.0d-2) call fatalerr ('DoTillage', 'For time being: all Z_tillage must equal Max_Z_tillage')
      end do

      ! determine entry point in tabulated tillage events based on start value t1900; check if input dates are sorted
      call set_iTill (state)

      ! determine number of horizon at depth Max_Z_tillage (MaxNumSoilHo)
      call det_MNSH (state)

      ! for special case i_n_model = 3: calculate slope per soil layer (remains constant over time)
      tl%Slope_match = 0.0d0

      ! TO ADD: CHECK THAT DEPTH OF EACH TILLAGE EVENT CORRESPONDS TO BOTTOM OF SOIL HORIZON; USER MAY NEED TO DEFINE MULTIPLE SUBS-HORIZONS WITHIN A SINGLE REAL SOIL HORIZON
      ! currently: require changes in horizon number (iSoilLayer) at depth Max_Z_tillage
      fine = .false.
      do i = 1, mesh%numlay
         if (mesh%botcom(i) == tl%MaxNumSoilCP) then
            fine = .true.
            exit
         end if
      end do
      if (.not. fine) call fatalerr_collected ('DoTillage', 'Bottom of soil horizon does not coincide with tillage depth(s)')

   case (2)
      ! RATE/STATE EVENT
      tl%Rho_last(1:tl%MaxNumSoilHo) = soil%bdens(1:tl%MaxNumSoilHo)

      if (TEST) then
         ! for technical test
         call DTDPST ("YEAR-MONTHST-DAY", time%t1900, STRNG)
         if (trim(STRNG) == "2016-Apr-05") then
            soil%bdens(1) = 1000.0d0
            call Change_MvGpars(state)
            Call Adapt_WC_H (TEST, state)
         else if (trim(STRNG) == "2016-Apr-11") then
            soil%bdens(1) = 1250.0d0
            call Change_MvGpars(state)
            Call Adapt_WC_H (TEST, state)
         else if (trim(STRNG) == "2016-Apr-18") then
            soil%bdens(1) = 1325.0d0        ! no ponding occurs
            !!!soil%bdens(1) = 1406.322d0   ! in this specific test this change causes ponding
            call Change_MvGpars(state)
            Call Adapt_WC_H (TEST, state)
         end if

      else
         if (Test2) then
            ! for technical test
            call DTDPST ("YEAR-MONTHST-DAY", time%t1900, STRNG)
            if (trim(STRNG) == "2005-Jun-05") then
               tl%Rho_cons(1) = 1350.0d0
               tl%K_R_cons(1) =  10.0d0
               call Change_MvGpars(state)
               Call Adapt_WC_H (TEST, state)
            end if
            if (trim(STRNG) == "2005-Oct-30") then
               tl%Rho_cons(1) = 1900.0d0
               tl%K_R_cons(1) =    0.1d0
               call Change_MvGpars(state)
               Call Adapt_WC_H (TEST, state)
            end if
         end if
         ! normal usage
         if (tl%iTill <= tl%Ntill .and. nint(time%t1900) == nint(tl%Date_tillage(tl%iTill))) then
            call DTDPST ("YEAR-MONTHST-DAY", time%t1900, STRNG)
            call Change_Tillage_Info (tl%iTill, state)
            call Change_Bdens(state)
            tl%iTill = tl%iTill + 1
         else
            call Consolidate_Bdens(state)
         end if

         call Change_MvGpars(state)

         call DTDPST ("YEAR-MONTHST-DAY", time%t1900, STRNG)
         if (trim(STRNG) == "2005-Apr-05") then
            !ParamVG(5,1) = -2.0d0
            !ParamVG(5,1) = -1.5d0
            !ParamVG(5,1) = -1.0d0
            !ParamVG(5,1) = 0.0d0
            !ParamVG(5,1) = 1.0d0
            !ParamVG(5,1) = 2.5d0
            !ParamVG(5,1) = 5.0d0
            !ParamVG(5,1) = 10.0d0
         end if

         Call Adapt_WC_H (TEST, state)

      end if


   case (3)
      ! OUTPUT

      ! [SS-BMI2] build tillage output row buffer (always runs, headless-independent)
      call build_tillage_output_row(state)

      ! [SS-BMI2] headless guard: debug writes to units 222/224/226 gated
      if (.not. time%headless) then
         if (TEST) then
            call DTDPST ("YEAR-MONTHST-DAY", time%t1900, STRNG)
            write (222,'(A,F15.5,10(I3,F15.5))') trim(time%date), atmo%nraida, (i, soil%bdens(i), i = 1, tl%MaxNumSoilHo)
            write (224,'(A,10F15.5)') trim(time%date), soil%theta(5), soil%theta(10), soil%theta(20), soil%theta(27), soil%theta(35), atmo%nraida, tl%sumDWC, tl%sumAvail1, tl%sumAvail2
            write (226,'(A,10F15.5)') trim(time%date), &
     &         soil%vg_params(1)%thetar, soil%vg_params(1)%thetas, &
     &         soil%vg_params(1)%ksat,   soil%vg_params(1)%alpha,  &
     &         soil%vg_params(1)%lpar,   soil%vg_params(1)%npar,   &
     &         soil%vg_params(1)%mpar,   soil%vg_params(1)%alphaw_sentinel, &
     &         soil%vg_params(1)%h_enpr, soil%vg_params(1)%ksatexm
         end if
         write (222,'(A,F15.5,10(I3,F15.5))') trim(time%date), atmo%nraida, (i, soil%bdens(i), i = 1, tl%MaxNumSoilHo)
         write (226,'(A,10F15.5)') trim(time%date), &
     &      soil%vg_params(1)%thetar, soil%vg_params(1)%thetas, &
     &      soil%vg_params(1)%ksat,   soil%vg_params(1)%alpha,  &
     &      soil%vg_params(1)%lpar,   soil%vg_params(1)%npar,   &
     &      soil%vg_params(1)%mpar,   soil%vg_params(1)%alphaw_sentinel, &
     &      soil%vg_params(1)%h_enpr, soil%vg_params(1)%ksatexm
      end if

   case (4)
      ! CLOSURE
      ! [SS-BMI2] deallocate tillage output buffer
      call cleanup_tillage_output_buffer(state)

   case default
      call fatalerr_collected ('DoTillage','Illegal value for iTask')
   end select

   end associate
   end subroutine DoTillage

! **************************************************** Change_MvGpars *********************************************************
   subroutine Change_MvGpars (state)
   implicit none
   type(swap_state_t), intent(inout) :: state
   integer              :: i, node, lay
   integer, parameter   :: Delta = 4
   integer, parameter   :: DeltaMin7 = Delta - 7
   real(8), parameter   :: Omega = -3.97d0
   real(8), parameter   :: Rho_s = 2650d0      ! later as input?
   real(8)              :: wcs_last, Epsilon
   associate( &
      mesh => state%mesh,         &
      soil => state%soilwater,    &
      tl   => state%tillage)
   do i = 1 , tl%MaxNumSoilHo
      wcs_last = ParamVG(2,i)       ! help

      ! First, fill PARAMVG
      ! wcr
      ParamVG(1,i) = ParamVG(1,i) * soil%bdens(i)/tl%Rho_last(i)
      ! wcs
      ParamVG(2,i) = ParamVG(2,i) * (Rho_s - soil%bdens(i))/(Rho_s - tl%Rho_last(i))
      ! ks,fit
      ParamVG(3,i) = ParamVG(3,i) * (ParamVG(2,i)/wcs_last)**3 * (soil%bdens(i)/tl%Rho_last(i))**DeltaMin7
      ! alpha
      ParamVG(4,i) = ParamVG(4,i) * (soil%bdens(i)/tl%Rho_last(i))**Omega
      ! lambda: do not change
      !ParamVG(5,i) = ParamVG(5,i)
      ! n
      select case (tl%i_n_model)
      case (1)
         !ParamVG(6,i) = ParamVG(6,i)
         continue
      case(2)
         Epsilon = -0.97d0 + 1.28d0 * soil%psilt(i) / soil%pclay(i)
         ParamVG(6,i) = 1.0d0 + (ParamVG(6,i) - 1.0d0) * (soil%bdens(i)/tl%Rho_last(i))**Epsilon
      case(3)
         ParamVG(6,i) = dmax1(1.001d0, ParamVG(6,i) + (soil%bdens(i) - tl%Rho_last(i)) * tl%Slope_match(i))
      end select

      ! m = 1-1/n
      ParamVG(7,i) = 1.0d0 - 1.0d0/ParamVG(6,i)
      ! alpha_w: not used; do not change
      !ParamVG(8,i) = ParamVG(8,i)
      ! h_enpr: do not change
      !ParamVG(9,i) = ParamVG(9,i)
      ! ksatexm: not used; do not change
      !ParamVG(10,i) = ParamVG(10,i)
   end do

   ! Second: update vg_params (first 10 fields, indices 1-10; index 8 is sentinel, not a real paramvg slot)
   do node = 1, tl%MaxNumSoilCP
      lay = mesh%layer(node)
      soil%vg_params(node)%thetar          = ParamVG(1, lay)
      soil%vg_params(node)%thetas          = ParamVG(2, lay)
      soil%vg_params(node)%ksat            = ParamVG(3, lay)
      soil%vg_params(node)%alpha           = ParamVG(4, lay)
      soil%vg_params(node)%lpar            = ParamVG(5, lay)
      soil%vg_params(node)%npar            = ParamVG(6, lay)
      soil%vg_params(node)%mpar            = ParamVG(7, lay)
      ! ParamVG(8) is alphaw_sentinel — keep existing sentinel value, not re-read from paramvg
      soil%vg_params(node)%h_enpr          = ParamVG(9, lay)
      soil%vg_params(node)%ksatexm         = ParamVG(10, lay)
      ! relsatthr (index 11) and ksatthr (index 12) are not changed by tillage
   end do
   end associate

   end subroutine Change_MvGpars

! **************************************************** Adapt_WC_H *********************************************************
   subroutine Adapt_WC_H (TEST, state)
   use soilhydraulics_utils, only: watcon, hconduc, prhead
   implicit none

   type(swap_state_t), intent(inout) :: state
   integer                          :: i
   real(8)                          :: sumWCtmin1, sumWCt, dwc, wcr, wcs, summ, dif
   real(8), dimension(state%tillage%MaxNumSoilCP) :: wc, hold, wcold
   logical                          :: TEST

   associate( &
      mesh => state%mesh,         &
      soil => state%soilwater,    &
      heat => state%heat,         &
      time => state%timecontrol,  &
      tl   => state%tillage)

   if (TEST) then
      hold(1:tl%MaxNumSoilCP)  = soil%h(1:tl%MaxNumSoilCP)
      wcold(1:tl%MaxNumSoilCP) = soil%theta(1:tl%MaxNumSoilCP)
   end if

   select case (tl%iRedist)
   case (0)
      if (.not.TEST) call fatalerr_collected ('Adapt_WC_H', 'Option iRedist = 0 only allowed in combination with TEST option')
      continue

   case (1)
      ! keep current wc values (wc_new = wc_old) and only change h; exception: when wc_old > wcs_new: alternative redistribution required
      summ = 0.0d0
      do i = 1, tl%MaxNumSoilCP
         wcs = ParamVG(2,mesh%layer(i))
         if (soil%theta(i) < wcs) then
            soil%h(i) = prhead(mesh%disnod(i), soil%theta(i), soil%h, &
                                          soil%iHWCKmodel(soil%layer(i)), &
                                          i, soil)
         else
            summ = summ + (wcs - soil%theta(i))*mesh%dz(i)
            soil%theta(i) = wcs
            soil%h(i) = 0.0d0
         end if
      end do
      if (summ > 0.0d0) then
         do i = tl%MaxNumSoilCP, 1, -1
            wcs = ParamVG(2,mesh%layer(i))
            dif = wcs - soil%theta(i)
            if (dif > 0.0d0) then
               if (dif < summ) then
                  soil%theta(i) = wcs
                  summ = summ - dif
               else
                  soil%theta(i) = soil%theta(i) + dif
                  summ = 0.0d0
                  exit
               end if
            end if
         end do
      end if
      soil%pond = summ

   case (2)
      sumWCtmin1 = sum(soil%theta(1:tl%MaxNumSoilCP))
      sumWCt = 0.0d0
      tl%sumDWC = 0.0d0
      do i = 1, tl%MaxNumSoilCP
         wc(i) = watcon(soil%h(i), &
                         soil%vg_params(i), &
                         soil%iHWCKmodel(soil%layer(i)), &
                         i, soil)
         sumWCt = sumWCt + wc(i)
         dwc = soil%theta(i) - wc(i)
         tl%sumDWC = tl%sumDWC + dwc * mesh%dz(i)
      end do

      tl%sumAvail1 = 0.0d0
      tl%sumAvail2 = 0.0d0
      if (sumWCt < sumWCtmin1) then
         ! water to be added; same as sumDWC > 0.0
         do i = 1, tl%MaxNumSoilCP
            wcs = ParamVG(2,mesh%layer(i))
            tl%sumAvail1 = tl%sumAvail1 + (wcs - wc(i))*mesh%dz(i)
         end do
         do i = 1, tl%MaxNumSoilCP
            wcs = ParamVG(2,mesh%layer(i))
            if (tl%sumAvail1 > 0.0d0) then
               wc(i) = wc(i) + (wcs - wc(i)) * tl%sumDWC / tl%sumAvail1
               if (wc(i) > wcs) then
                  soil%pond = soil%pond + (wc(i) - wcs) * mesh%dz(i)
                  wc(i) = wcs
                  write(333,'(A,I5,F12.4)') time%date, i, soil%pond
               end if
            end if
            soil%h(i)     = prhead(mesh%disnod(i), wc(i), soil%h, &
                                               soil%iHWCKmodel(soil%layer(i)), &
                                               i, soil)
            soil%theta(i) = wc(i)
         end do
      else if (sumWCt > sumWCtmin1) then
         ! water to be removed; same as sumDWC < 0.0
         do i = 1, tl%MaxNumSoilCP
            wcr = ParamVG(1,mesh%layer(i))
            tl%sumAvail2 = tl%sumAvail2 + (wc(i) - wcr) * mesh%dz(i)
         end do
         do i = 1, tl%MaxNumSoilCP
            wcr = ParamVG(1,mesh%layer(i))
            wc(i) = wc(i) + (wc(i) - wcr) * tl%sumDWC / tl%sumAvail2
            soil%h(i)     = prhead(mesh%disnod(i), wc(i), soil%h, &
                                               soil%iHWCKmodel(soil%layer(i)), &
                                               i, soil)
            soil%theta(i) = wc(i)
         end do

      endif

   end select

   if (TEST) then
      do i = 1, tl%MaxNumSoilCP
         wcs = ParamVG(2,mesh%layer(i))
         write (444,'(I5,8(A1,F12.6))') i, ',', hold(i), ',', wcold(i), ',', soil%h(i), ',', soil%theta(i), &
                                        ',', tl%sumDWC, ',', tl%sumAvail1, ',', tl%sumAvail2, &
                                        ',', soil%theta(i)/wcs
      end do
   end if
write(124,'(A,1P,12E12.5)') time%date, soil%bdens(1), ParamVG(2,mesh%layer(1)), soil%theta(1), soil%h(1), &
   hconduc(soil%h(1),soil%theta(1),1.0d0,heat%tsoil(1), &
           soil%vg_params(1), &
           soil%iHWCKmodel(soil%layer(1)), &
           soil%fluseksatexm(1), 1, soil), ParamVG(3,mesh%layer(1)),          &
   soil%bdens(2), ParamVG(2,mesh%layer(2)), soil%theta(2), soil%h(2),                           &
   hconduc(soil%h(2),soil%theta(2),1.0d0,heat%tsoil(2), &
           soil%vg_params(2), &
           soil%iHWCKmodel(soil%layer(2)), &
           soil%fluseksatexm(2), 2, soil), ParamVG(3,mesh%layer(2))

   end associate
   end subroutine Adapt_WC_H


! **************************************************** Consolidate_Bdens *********************************************************
   subroutine Consolidate_Bdens (state)
   implicit none
   type(swap_state_t), intent(inout) :: state
   integer :: i

   associate( &
      soil => state%soilwater,    &
      atmo => state%atmosphere,   &
      time => state%timecontrol,  &
      tl   => state%tillage)
   if (tl%iTill == 1) return        ! in the beginning before first tillage event: do nothing

   forall (i=1:tl%MaxNumSoilHo) soil%bdens(i) = tl%Rho_cons(i) - (tl%Rho_cons(i) - tl%Rho_last(i)) * dexp(-tl%K_R_cons(i)*atmo%nraida*10.0d0)    ! 10: to transform nraida from cm to mm
   write(123,'(A,1P,10E12.5)') time%date, atmo%nraida, soil%bdens(1:tl%MaxNumSoilHo)
   end associate
   end subroutine Consolidate_Bdens

! **************************************************** Change_Bdens *********************************************************
   subroutine Change_Bdens (state)
   implicit none
   type(swap_state_t), intent(inout) :: state
   integer :: i, NumSoilHo
   ! to check: why this loop to determine NumSoilHo?
   NumSoilHo = 1
   associate( &
      mesh => state%mesh,         &
      soil => state%soilwater,    &
      tl   => state%tillage)
   do i = 2, mesh%numnod
      if (tl%Z_tillage(tl%iTill) > -mesh%zbotcp(i-1) .and. tl%Z_tillage(tl%iTill) <= -mesh%zbotcp(i)) then
         NumSoilHo = mesh%layer(i)
         exit
      end if
   end do
   forall (i=1:NumSoilHo) soil%bdens(i) = tl%Rho_last(i) - tl%I_tillage(tl%iTill) * (tl%Rho_last(i) - tl%Rho_tillage(i))
   end associate

   end subroutine Change_Bdens

! **************************************************** set_iTill *********************************************************
   subroutine set_iTill (state)
   implicit none
   type(swap_state_t), intent(inout) :: state
   integer :: i
   associate( &
      time => state%timecontrol,  &
      tl   => state%tillage)
   tl%iTill = 1
   if (time%t1900 <= tl%Date_tillage(1)) tl%iTill = 1
   do i = 2, tl%Ntill
      if (tl%Date_tillage(i) < tl%Date_tillage(i-1)) call fatalerr_collected ('set_iTill', 'Dates in tabulated tillage events must be sorted')
      ! H-4 bug fix: second comparison was Date_tillage(i-1) (tautological); corrected to Date_tillage(i)
      if (time%t1900 >= tl%Date_tillage(i-1) .and. time%t1900 < tl%Date_tillage(i)) tl%iTill = i-1
   end do
   end associate
   end subroutine set_iTill

! **************************************************** det_MNSH *********************************************************
   subroutine det_MNSH (state)
   implicit none
   type(swap_state_t), intent(inout) :: state
   integer :: i
   associate( &
      mesh => state%mesh,         &
      tl   => state%tillage)
   do i = 2, mesh%numnod
      if (tl%Max_Z_tillage > -mesh%zbotcp(i-1) .and. tl%Max_Z_tillage <= -mesh%zbotcp(i)) then
         tl%MaxNumSoilHo = mesh%layer(i)
         tl%MaxNumSoilCP = i
         exit
      end if
   end do
   end associate
   end subroutine det_MNSH

   subroutine Change_Tillage_Info (iTill_in, state)
   implicit none
   ! global
   integer, intent(in) :: iTill_in
   type(swap_state_t), intent(inout) :: state
   ! local
   integer :: itype, nlay

   associate(tl => state%tillage)
   itype                  = tl%Type_Tillage(iTill_in)
   nlay                   = tl%iTT2(itype) - tl%iTT1(itype) + 1
   tl%Rho_tillage(1:nlay) = tl%TAB_Rho_tillage(tl%iTT1(itype):tl%iTT2(itype))
   tl%Rho_cons(1:nlay)    = tl%TAB_Rho_cons(tl%iTT1(itype):tl%iTT2(itype))
   tl%K_R_cons(1:nlay)    = tl%TAB_K_R_cons(tl%iTT1(itype):tl%iTT2(itype))
   tl%Rho_match(1:nlay)   = tl%TAB_Rho_match(tl%iTT1(itype):tl%iTT2(itype))
   tl%N_match(1:nlay)     = tl%TAB_N_match(tl%iTT1(itype):tl%iTT2(itype))
   tl%Slope_match(1:nlay) = (ParamVG(6,1:nlay) - tl%N_match(1:nlay)) / (tl%Rho_cons(1:nlay) - tl%Rho_match(1:nlay))
   end associate
   end subroutine Change_Tillage_Info


! ----------------------------------------------------------------------
! [SS-BMI2] Tillage output buffer helpers (canonical output-sink pattern)
! N = 5: t1900, nraida, sumDWC, sumAvail1, sumAvail2
! Debug writes to units 222/224/226 are gated by headless in DoTillage(3).
! Cleanup in DoTillage(4) (closure; not currently called in production).
! ----------------------------------------------------------------------

   subroutine init_tillage_output_buffer(state)
! ----------------------------------------------------------------------
!     Allocate state%tillage%output_row and set column names.
!     Called from DoTillage(1) — always runs, headless-independent.
!     N = 5: t1900, nraida, sumDWC, sumAvail1, sumAvail2
! ----------------------------------------------------------------------
   use iso_c_binding, only: c_double
   implicit none
   type(swap_state_t), intent(inout) :: state
   integer, parameter :: N = 5

   state%tillage%output_n_cols = N
   if (.not. allocated(state%tillage%output_row))     allocate(state%tillage%output_row(N))
   if (.not. allocated(state%tillage%output_columns)) allocate(state%tillage%output_columns(N))
   state%tillage%output_row     = 0.0_c_double
   state%tillage%output_columns(1) = 'date'
   state%tillage%output_columns(2) = 'nraida'
   state%tillage%output_columns(3) = 'sumDWC'
   state%tillage%output_columns(4) = 'sumAvail1'
   state%tillage%output_columns(5) = 'sumAvail2'
   end subroutine init_tillage_output_buffer


   subroutine build_tillage_output_row(state)
! ----------------------------------------------------------------------
!     Fill state%tillage%output_row(:) with the current tillage values.
!     Called from DoTillage(3) — always runs, headless-independent.
!     Column order: t1900, nraida, sumDWC, sumAvail1, sumAvail2.
! ----------------------------------------------------------------------
   use iso_c_binding, only: c_double
   implicit none
   type(swap_state_t), intent(inout) :: state

   state%tillage%output_row(1) = real(state%timecontrol%t1900,       c_double)
   state%tillage%output_row(2) = real(state%atmosphere%nraida,       c_double)
   state%tillage%output_row(3) = real(state%tillage%sumDWC,          c_double)
   state%tillage%output_row(4) = real(state%tillage%sumAvail1,       c_double)
   state%tillage%output_row(5) = real(state%tillage%sumAvail2,       c_double)
   end subroutine build_tillage_output_row


   subroutine cleanup_tillage_output_buffer(state)
! ----------------------------------------------------------------------
!     Deallocate state%tillage%output_row and reset counter.
!     Called from DoTillage(4) — closure; not currently called in production.
! ----------------------------------------------------------------------
   implicit none
   type(swap_state_t), intent(inout) :: state

   if (allocated(state%tillage%output_row))     deallocate(state%tillage%output_row)
   if (allocated(state%tillage%output_columns)) deallocate(state%tillage%output_columns)
   state%tillage%output_n_cols = 0
   end subroutine cleanup_tillage_output_buffer

end module tillage_mod
