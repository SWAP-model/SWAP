! to do: default MvG parameters that are read: do they refer to BDENS or Rho_cons; or should BDENS and RhoCons be equal?
! to do: check if Rho_match differs from BDENS or Rho_cons (if not: division by zero possible)

module tillage_mod
   use error_mod, only: fatalerr_collected
   use swap_state_mod, only: swap_state_t  ! [SS-ATM A-2.6] nraida retired from variables to state%atmosphere

   ! SS-TC TC-12: t1900 retired from only-list; read via state%timecontrol%t1900 at each call site.
   use variables, only: date, swhyst, swsolu, swoxygen, flMacroPore, flksatexm, zbotcp, NumNod, Bdens, layer, ParamVG, &
                        NumLay, dz, disnod, botcom, psilt, pclay, SwDiscrvert, tend, &  ! [SS-SWC S-2.6] CofGen/pond/theta/h retired to state%soilwater
                        ! Tillage bridge variables with renaming (SAVE statements removed)
                        ! [SS-TIL T-5] Groups C/D/E retired from variables — reads via state%tillage
                        swtill => till_swtill, Ntill => till_Ntill, &
                        Ntypes => till_Ntypes, i_n_model => till_i_n_model, iRedist => till_iRedist, &
                        Max_Z_tillage => till_Max_Z_tillage, &
                        Date_tillage => till_Date_tillage, Z_tillage => till_Z_tillage, &
                        I_tillage => till_I_tillage, Type_Tillage => till_Type_Tillage, &
                        iType_Tillage => till_iType_Tillage, iTT1 => till_iTT1, iTT2 => till_iTT2, &
                        TAB_Rho_tillage => till_TAB_Rho_tillage, TAB_Rho_cons => till_TAB_Rho_cons, &
                        TAB_K_R_cons => till_TAB_K_R_cons, &
                        TAB_Rho_match => till_TAB_Rho_match, TAB_N_match => till_TAB_N_match
   
   implicit none

   ! SAVE removed - persistent state now in variables module with till_ prefix
   ! Variables are imported with their original names via renaming in the USE statement above

!  by default: all in this module is private (local)
   private
!  except for these public routines/functions
   public :: DoTillage
!  and except for these public variables
   public :: swtill  ! Alias to till_swtill from variables module
   
   contains

   subroutine DoTillage (iTask, state)
   ! global
   integer, intent(in)                          :: iTask                               ! Task
   type(swap_state_t), intent(inout)            :: state                               ! [SS-ATM A-2.6] for retired nraida; [SS-SWC S-1.9] inout for soilwater dual-writes
   ! local (not to be saved)
   integer                                   :: i
   character(len=20)                         :: STRNG
   logical                                   :: fine
!   logical, parameter                        :: TEST = .true.
   logical, parameter                        :: TEST = .false.
!   logical, parameter                        :: TEST2 = .true.
   logical, parameter                        :: TEST2 = .false.

   ! functions

   if (iTask > 1 .and. swtill == 0) return      ! no tillage to be considered: return immediately

   ! handle iTask
   select case (iTask)
   case (1)
      ! INITIALIZE

      ! [SS-TIL T-5] legacy allocates dropped; state%tillage arrays already allocated in tillage_init (T-1/T-2)

      ! some checks: some combinations not (yet) allowed
      if (swtill == 1) then
         if (swhyst == 1)      call fatalerr_collected ('DoTillage', 'swhyst = 1 not allowed')
         if (swsolu == 1)      call fatalerr_collected ('DoTillage', 'swsolu = 1 not (yet) allowed')
         if (swoxygen == 2)    call fatalerr_collected ('DoTillage', 'swoxygen = 2 not (yet) allowed')
         if (flMacroPore)      call fatalerr_collected ('DoTillage', 'swmacro = 1 not (yet) allowed')
         if (flksatexm)        call fatalerr_collected ('DoTillage', 'flksatexm not (yet) allowed')
         if (SwDiscrvert == 1) call fatalerr_collected ('DoTillage', 'SwDiscrvert = 1 not (yet) allowed')
      end if

      ! currently: require all Z_tillage = Max_Z_tillage
      do i = 1, Ntill
!         if (dabs(Max_Z_tillage - Z_tillage(i)) > 1.0d-2) call fatalerr ('DoTillage', 'For time being: all Z_tillage must equal Max_Z_tillage')
      end do

      ! determine entry point in tabulated tillage events based on start value t1900; check if input dates are sorted
      call set_iTill (state)                         ! [SS-TIL T-3] state for iTill dual-write

      ! determine number of horizon at depth Max_Z_tillage (MaxNumSoilHo)
      call det_MNSH (state)                          ! [SS-TIL T-3] state for MaxNumSoilHo/MaxNumSoilCP dual-write

      ! for special case i_n_model = 3: calculate slope per soil layer (remains constant over time)
      state%tillage%Slope_match = 0.0d0              ! [SS-TIL T-5] legacy Slope_match dropped; canonical via state%tillage

      ! TO ADD: CHECK THAT DEPTH OF EACH TILLAGE EVENT CORRESPONDS TO BOTTOM OF SOIL HORIZON; USER MAY NEED TO DEFINE MULTIPLE SUBS-HORIZONS WITHIN A SINGLE REAL SOIL HORIZON
      ! currently: require changes in horizon number (iSoilLayer) at depth Max_Z_tillage
      fine = .false.
      do i = 1, NumLay
         if (botcom(i) == state%tillage%MaxNumSoilCP) then
            fine = .true.
            exit
         end if
      end do
      if (.not. fine) call fatalerr_collected ('DoTillage', 'Bottom of soil horizon does not coincide with tillage depth(s)')
      
   case (2)
      ! RATE/STATE EVENT
      ! [SS-TIL T-5] legacy Rho_last dropped; canonical via state%tillage
      state%tillage%Rho_last(1:state%tillage%MaxNumSoilHo) = Bdens(1:state%tillage%MaxNumSoilHo)

      if (TEST) then
         ! for technical test
         call DTDPST ("YEAR-MONTHST-DAY", state%timecontrol%t1900, STRNG)  ! TC-12
         if (trim(STRNG) == "2016-Apr-05") then
            BDENS(1) = 1000.0d0
            call Change_MvGpars(state)              ! [SS-SWC S-1.9]
            Call Adapt_WC_H (TEST, state)           ! [SS-SWC S-1.9]
         else if (trim(STRNG) == "2016-Apr-11") then
            BDENS(1) = 1250.0d0
            call Change_MvGpars(state)              ! [SS-SWC S-1.9]
            Call Adapt_WC_H (TEST, state)           ! [SS-SWC S-1.9]
         else if (trim(STRNG) == "2016-Apr-18") then
            BDENS(1) = 1325.0d0        ! no ponding occurs
            !!!BDENS(1) = 1406.322d0   ! in this specific test this change causes ponding
            call Change_MvGpars(state)              ! [SS-SWC S-1.9]
            Call Adapt_WC_H (TEST, state)           ! [SS-SWC S-1.9]
         end if

      else
         if (Test2) then
            ! for technical test
            call DTDPST ("YEAR-MONTHST-DAY", state%timecontrol%t1900, STRNG)  ! TC-12
            if (trim(STRNG) == "2005-Jun-05") then
               ! [SS-TIL T-5] legacy Rho_cons/K_R_cons dropped; canonical via state%tillage
               state%tillage%Rho_cons(1) = 1350.0d0
               state%tillage%K_R_cons(1) =  10.0d0
               call Change_MvGpars(state)           ! [SS-SWC S-1.9]
               Call Adapt_WC_H (TEST, state)        ! [SS-SWC S-1.9]
            end if
            if (trim(STRNG) == "2005-Oct-30") then
               ! [SS-TIL T-5] legacy Rho_cons/K_R_cons dropped; canonical via state%tillage
               state%tillage%Rho_cons(1) = 1900.0d0
               state%tillage%K_R_cons(1) =    0.1d0
               call Change_MvGpars(state)           ! [SS-SWC S-1.9]
               Call Adapt_WC_H (TEST, state)        ! [SS-SWC S-1.9]
            end if
         end if
         ! normal usage  [SS-TIL T-5] iTill reads via state%tillage
         if (state%tillage%iTill <= Ntill .and. nint(state%timecontrol%t1900) == nint(Date_tillage(state%tillage%iTill))) then  ! TC-12
            call DTDPST ("YEAR-MONTHST-DAY", state%timecontrol%t1900, STRNG)  ! TC-12
            call Change_Tillage_Info (state%tillage%iTill, state)    ! [SS-TIL T-3] state for Group C dual-writes
            call Change_Bdens(state)
            state%tillage%iTill = state%tillage%iTill + 1   ! [SS-TIL T-5] legacy iTill dropped
         else
            ! SS-ATM A-2.6: pass state for retired nraida
            call Consolidate_Bdens(state)
         end if

         call Change_MvGpars(state)                 ! [SS-SWC S-1.9]

         call DTDPST ("YEAR-MONTHST-DAY", state%timecontrol%t1900, STRNG)  ! TC-12
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

         
         
         Call Adapt_WC_H (TEST, state)              ! [SS-SWC S-1.9]

      end if
      

   case (3)
      ! OUTPUT  [SS-TIL T-5] MaxNumSoilHo/sumDWC/sumAvail1/sumAvail2 reads via state%tillage
      if (TEST) then
         call DTDPST ("YEAR-MONTHST-DAY", state%timecontrol%t1900, STRNG)  ! TC-12
         write (222,'(A,F15.5,10(I3,F15.5))') trim(DATE), state%atmosphere%nraida, (i, Bdens(i), i = 1, state%tillage%MaxNumSoilHo)
         write (224,'(A,10F15.5)') trim(DATE), state%soilwater%theta(5), state%soilwater%theta(10), state%soilwater%theta(20), state%soilwater%theta(27), state%soilwater%theta(35), state%atmosphere%nraida, state%tillage%sumDWC, state%tillage%sumAvail1, state%tillage%sumAvail2  ! [SS-SWC S-2.6]
         write (226,'(A,10F15.5)') trim(DATE), (state%soilwater%cofgen(i,1), i = 1, 10)  ! [SS-SWC S-2.6]
      end if
         write (222,'(A,F15.5,10(I3,F15.5))') trim(DATE), state%atmosphere%nraida, (i, Bdens(i), i = 1, state%tillage%MaxNumSoilHo)
         write (226,'(A,10F15.5)') trim(DATE), (state%soilwater%cofgen(i,1), i = 1, 10)  ! [SS-SWC S-2.6]
      continue

   case (4)
      ! CLOSURE
      continue
      
   case default
      call fatalerr_collected ('DoTillage','Illegal value for iTask')
   end select

   end subroutine DoTillage

! **************************************************** Change_MvGpars *********************************************************
   subroutine Change_MvGpars (state)                ! [SS-SWC S-2.6] cofgen reads cut over to state%soilwater%cofgen
   implicit none
   type(swap_state_t), intent(inout) :: state
   integer              :: i, node, lay
   integer, parameter   :: Delta = 4
   integer, parameter   :: DeltaMin7 = Delta - 7
   real(8), parameter   :: Omega = -3.97d0
   real(8), parameter   :: Rho_s = 2650d0      ! later as input?
   real(8)              :: wcs_last, Epsilon
   ! [SS-TIL T-5] MaxNumSoilHo/MaxNumSoilCP/Rho_last/Slope_match read via state%tillage
   associate(tl => state%tillage)
   do i = 1 , tl%MaxNumSoilHo
      wcs_last = ParamVG(2,i)       ! help

      ! First, fill PARAMVG
      ! wcr
      ParamVG(1,i) = ParamVG(1,i) * Bdens(i)/tl%Rho_last(i)
      ! wcs
      ParamVG(2,i) = ParamVG(2,i) * (Rho_s - Bdens(i))/(Rho_s - tl%Rho_last(i))
      ! ks,fit
      ParamVG(3,i) = ParamVG(3,i) * (ParamVG(2,i)/wcs_last)**3 * (Bdens(i)/tl%Rho_last(i))**DeltaMin7
      ! alpha
      ParamVG(4,i) = ParamVG(4,i) * (Bdens(i)/tl%Rho_last(i))**Omega
      ! lambda: do not change
      !ParamVG(5,i) = ParamVG(5,i)
      ! n
      select case (i_n_model)
      case (1)
         !ParamVG(6,i) = ParamVG(6,i)
         continue
      case(2)
         Epsilon = -0.97d0 + 1.28d0 * psilt(i) / pclay(i)
         ParamVG(6,i) = 1.0d0 + (ParamVG(6,i) - 1.0d0) * (Bdens(i)/tl%Rho_last(i))**Epsilon
      case(3)
         ParamVG(6,i) = dmax1(1.001d0, ParamVG(6,i) + (Bdens(i) - tl%Rho_last(i)) * tl%Slope_match(i))
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

   ! Second: fill COFGEN
   do node = 1, tl%MaxNumSoilCP
      lay = layer(node)
      state%soilwater%cofgen(1:10,node) = ParamVG(1:10,lay)  ! [SS-SWC S-2.6] legacy CofGen write dropped
      ! CofGen(11) and CofGen(12) are not used and not need to be changed
      !CofGen(11,node) = relsatthr(lay)
      !CofGen(12,node) = ksatthr(lay)
   end do
   end associate
   !!!thetsl(1:numlay) = ParamVG(2,1:numlay)

   end subroutine Change_MvGpars
   
! **************************************************** Adapt_WC_H *********************************************************
   subroutine Adapt_WC_H (TEST, state)                ! [SS-SWC S-2.6] theta/h/pond reads cut over to state%soilwater
   use soilhydraulics_utils, only: watcon, hconduc, prhead
   implicit none

   type(swap_state_t), intent(inout) :: state
   integer                          :: i
   real(8)                          :: sumWCtmin1, sumWCt, dwc, wcr, wcs, summ, dif
   ! [SS-TIL T-5] MaxNumSoilCP/sumDWC/sumAvail1/sumAvail2 via state%tillage
   real(8), dimension(state%tillage%MaxNumSoilCP) :: wc, hold, wcold
   logical                          :: TEST

   if (TEST) then
      hold(1:state%tillage%MaxNumSoilCP)  = state%soilwater%h(1:state%tillage%MaxNumSoilCP)      ! [SS-SWC S-2.6]
      wcold(1:state%tillage%MaxNumSoilCP) = state%soilwater%theta(1:state%tillage%MaxNumSoilCP)  ! [SS-SWC S-2.6]
   end if

   select case (iRedist)
   case (0)
      if (.not.TEST) call fatalerr_collected ('Adapt_WC_H', 'Option iRedist = 0 only allowed in combination with TEST option')
      continue

   case (1)
      ! keep current wc values (wc_new = wc_old) and only change h; exception: when wc_old > wcs_new: alternative redistribution required
      summ = 0.0d0
      do i = 1, state%tillage%MaxNumSoilCP
         wcs = ParamVG(2,layer(i))
         if (state%soilwater%theta(i) < wcs) then                                                       ! [SS-SWC S-2.6]
            state%soilwater%h(i) = prhead(i, disnod(i), state%soilwater%theta(i), &                    ! [SS-SWC S-2.6]
                                          state%soilwater%cofgen, state%soilwater%h)                    ! [SS-SWC S-2.6]
         else
            summ = summ + (wcs - state%soilwater%theta(i))*dz(i)                                       ! [SS-SWC S-2.6]
            state%soilwater%theta(i) = wcs                                                              ! [SS-SWC S-2.6]
            state%soilwater%h(i) = 0.0d0                                                               ! [SS-SWC S-2.6]
         end if
      end do
      if (summ > 0.0d0) then
         do i = state%tillage%MaxNumSoilCP, 1, -1
            wcs = ParamVG(2,layer(i))
            dif = wcs - state%soilwater%theta(i)                                                        ! [SS-SWC S-2.6]
            if (dif > 0.0d0) then
               if (dif < summ) then
                  state%soilwater%theta(i) = wcs                                                        ! [SS-SWC S-2.6]
                  summ = summ - dif
               else
                  state%soilwater%theta(i) = state%soilwater%theta(i) + dif                            ! [SS-SWC S-2.6]
                  summ = 0.0d0
                  exit
               end if
            end if
         end do
      end if
      state%soilwater%pond = summ                                                                        ! [SS-SWC S-2.6]

   case (2)
      ! [SS-TIL T-5] sumDWC/sumAvail1/sumAvail2 written directly to state%tillage (legacy scalars retired)
      sumWCtmin1 = sum(state%soilwater%theta(1:state%tillage%MaxNumSoilCP))                             ! [SS-SWC S-2.6]
      sumWCt = 0.0d0
      state%tillage%sumDWC = 0.0d0
      do i = 1, state%tillage%MaxNumSoilCP
         wc(i) = watcon(i,state%soilwater%h(i))                                                         ! [SS-SWC S-2.6]
         sumWCt = sumWCt + wc(i)
         dwc = state%soilwater%theta(i) - wc(i)                                                         ! [SS-SWC S-2.6]
         state%tillage%sumDWC = state%tillage%sumDWC + dwc * dz(i)
      end do

      state%tillage%sumAvail1 = 0.0d0
      state%tillage%sumAvail2 = 0.0d0
      if (sumWCt < sumWCtmin1) then
         ! water to be added; same as sumDWC > 0.0
         do i = 1, state%tillage%MaxNumSoilCP
            wcs = ParamVG(2,layer(i))
            state%tillage%sumAvail1 = state%tillage%sumAvail1 + (wcs - wc(i))*dz(i)
         end do
         do i = 1, state%tillage%MaxNumSoilCP
            wcs = ParamVG(2,layer(i))
            if (state%tillage%sumAvail1 > 0.0d0) then
               wc(i) = wc(i) + (wcs - wc(i)) * state%tillage%sumDWC / state%tillage%sumAvail1
               if (wc(i) > wcs) then
                  state%soilwater%pond = state%soilwater%pond + (wc(i) - wcs) * dz(i)                  ! [SS-SWC S-2.6]
                  wc(i) = wcs
                  write(333,'(A,I5,F12.4)') Date, i, state%soilwater%pond                              ! [SS-SWC S-2.6]
               end if
            end if
            state%soilwater%h(i)     = prhead(i, disnod(i), wc(i), &                                   ! [SS-SWC S-2.6]
                                               state%soilwater%cofgen, state%soilwater%h)               ! [SS-SWC S-2.6]
            state%soilwater%theta(i) = wc(i)                                                            ! [SS-SWC S-2.6]
         end do
      else if (sumWCt > sumWCtmin1) then
         ! water to be removed; same as sumDWC < 0.0
         do i = 1, state%tillage%MaxNumSoilCP
            wcr = ParamVG(1,layer(i))
            state%tillage%sumAvail2 = state%tillage%sumAvail2 + (wc(i) - wcr) * dz(i)
         end do
         do i = 1, state%tillage%MaxNumSoilCP
            wcr = ParamVG(1,layer(i))
            wc(i) = wc(i) + (wc(i) - wcr) * state%tillage%sumDWC / state%tillage%sumAvail2
            state%soilwater%h(i)     = prhead(i, disnod(i), wc(i), &                                   ! [SS-SWC S-2.6]
                                               state%soilwater%cofgen, state%soilwater%h)               ! [SS-SWC S-2.6]
            state%soilwater%theta(i) = wc(i)                                                            ! [SS-SWC S-2.6]
         end do

      endif

   end select

   if (TEST) then
      do i = 1, state%tillage%MaxNumSoilCP
         wcs = ParamVG(2,layer(i))
         write (444,'(I5,8(A1,F12.6))') i, ',', hold(i), ',', wcold(i), ',', state%soilwater%h(i), ',', state%soilwater%theta(i), &  ! [SS-SWC S-2.6]
                                        ',', state%tillage%sumDWC, ',', state%tillage%sumAvail1, ',', state%tillage%sumAvail2, &
                                        ',', state%soilwater%theta(i)/wcs               ! [SS-SWC S-2.6]
      end do
   end if
write(124,'(A,1P,12E12.5)') Date, Bdens(1), ParamVG(2,layer(1)), state%soilwater%theta(1), state%soilwater%h(1), &  ! [SS-SWC S-2.6]
   hconduc(1,state%soilwater%h(1),state%soilwater%theta(1),1.0d0,state%heat%tsoil(1)), ParamVG(3,layer(1)),          & ! [SS-SWC S-2.6]
   Bdens(2), ParamVG(2,layer(2)), state%soilwater%theta(2), state%soilwater%h(2),                                     & ! [SS-SWC S-2.6]
   hconduc(2,state%soilwater%h(2),state%soilwater%theta(2),1.0d0,state%heat%tsoil(2)), ParamVG(3,layer(2))              ! [SS-SWC S-2.6]

   end subroutine Adapt_WC_H

   
! **************************************************** Consolidate_Bdens *********************************************************
   subroutine Consolidate_Bdens (state)
   ! [SS-ATM A-2.6] state added for retired nraida
   implicit none
   type(swap_state_t), intent(in) :: state
   integer :: i

   ! [SS-TIL T-5] iTill/MaxNumSoilHo/Rho_cons/Rho_last/K_R_cons read via state%tillage
   associate(tl => state%tillage)
   if (tl%iTill == 1) return        ! in the beginning before first tillage event: do nothing

   forall (i=1:tl%MaxNumSoilHo) Bdens(i) = tl%Rho_cons(i) - (tl%Rho_cons(i) - tl%Rho_last(i)) * dexp(-tl%K_R_cons(i)*state%atmosphere%nraida*10.0d0)    ! 10: to transform nraida from cm to mm
   write(123,'(A,1P,10E12.5)') Date, state%atmosphere%nraida, Bdens(1:tl%MaxNumSoilHo)
   end associate
   end subroutine Consolidate_Bdens
   
! **************************************************** Change_Bdens *********************************************************
   subroutine Change_Bdens (state)
   ! [SS-TIL T-5] state added; iTill/Rho_last/Rho_tillage read via state%tillage
   implicit none
   type(swap_state_t), intent(in) :: state
   integer :: i, NumSoilHo
   ! to check: why this loop to determine NumSoilHo?
   NumSoilHo = 1
   associate(tl => state%tillage)
   do i = 2, NumNod
      if (Z_tillage(tl%iTill) > -zbotcp(i-1) .and. Z_tillage(tl%iTill) <= -zbotcp(i)) then
         NumSoilHo = layer(i)
         exit
      end if
   end do
   forall (i=1:NumSoilHo) Bdens(i) = tl%Rho_last(i) - I_tillage(tl%iTill) * (tl%Rho_last(i) - tl%Rho_tillage(i))
   end associate

   end subroutine Change_Bdens
   
! **************************************************** set_iTill *********************************************************
   subroutine set_iTill (state)
   ! [SS-TIL T-5] legacy iTill dropped; canonical write directly to state%tillage%iTill
   implicit none
   type(swap_state_t), intent(inout) :: state
   integer :: i
   state%tillage%iTill = 1
   if (state%timecontrol%t1900 <= Date_tillage(1)) state%tillage%iTill = 1  ! TC-12
   do i = 2, Ntill
      if (Date_tillage(i) < Date_tillage(i-1)) call fatalerr_collected ('set_iTill', 'Dates in tabulated tillage events must be sorted')
      ! H-4 bug fix: second comparison was Date_tillage(i-1) (tautological); corrected to Date_tillage(i)
      if (state%timecontrol%t1900 >= Date_tillage(i-1) .and. state%timecontrol%t1900 < Date_tillage(i)) state%tillage%iTill = i-1  ! TC-12
   end do
   end subroutine set_iTill

! **************************************************** det_MNSH *********************************************************
   subroutine det_MNSH (state)
   ! [SS-TIL T-5] legacy MaxNumSoilHo/MaxNumSoilCP dropped; canonical write directly to state%tillage
   implicit none
   type(swap_state_t), intent(inout) :: state
   integer :: i
   do i = 2, NumNod
      if (Max_Z_tillage > -zbotcp(i-1) .and. Max_Z_tillage <= -zbotcp(i)) then
         state%tillage%MaxNumSoilHo = layer(i)
         state%tillage%MaxNumSoilCP = i
         exit
      end if
   end do
   end subroutine det_MNSH

   subroutine Change_Tillage_Info (iTill_in, state)
   ! [SS-TIL T-5] legacy Group C arrays dropped; write directly to state%tillage
   implicit none
   ! global
   integer, intent(in) :: iTill_in
   type(swap_state_t), intent(inout) :: state
   ! local
   integer :: itype, nlay

   associate(tl => state%tillage)
   itype               = Type_Tillage(iTill_in)
   nlay                = iTT2(itype) - iTT1(itype) + 1
   tl%Rho_tillage(1:nlay) = TAB_Rho_tillage(iTT1(itype):iTT2(itype))
   tl%Rho_cons(1:nlay)    = TAB_Rho_cons(iTT1(itype):iTT2(itype))
   tl%K_R_cons(1:nlay)    = TAB_K_R_cons(iTT1(itype):iTT2(itype))
   tl%Rho_match(1:nlay)   = TAB_Rho_match(iTT1(itype):iTT2(itype))
   tl%N_match(1:nlay)     = TAB_N_match(iTT1(itype):iTT2(itype))
   tl%Slope_match(1:nlay) = (ParamVG(6,1:nlay) - tl%N_match(1:nlay)) / (tl%Rho_cons(1:nlay) - tl%Rho_match(1:nlay))
   end associate
   end subroutine Change_Tillage_Info

end module tillage_mod
   