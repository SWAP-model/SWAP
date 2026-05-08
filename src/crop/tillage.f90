! to do: default MvG parameters that are read: do they refer to BDENS or Rho_cons; or should BDENS and RhoCons be equal?
! to do: check if Rho_match differs from BDENS or Rho_cons (if not: division by zero possible)

module tillage_mod
   use error_mod, only: fatalerr_collected

   use variables, only: t1900, date, swhyst, swsolu, swoxygen, flMacroPore, flksatexm, zbotcp, NumNod, Bdens, layer, nraida, ParamVG, CofGen, &
                        NumLay, pond, theta, h, dz, disnod, botcom, psilt, pclay, SwDiscrvert, tend, &
                        ! Tillage bridge variables with renaming (SAVE statements removed)
                        swtill => till_swtill, Ntill => till_Ntill, iTill => till_iTill, &
                        Ntypes => till_Ntypes, i_n_model => till_i_n_model, iRedist => till_iRedist, &
                        MaxNumSoilHo => till_MaxNumSoilHo, MaxNumSoilCP => till_MaxNumSoilCP, &
                        Max_Z_tillage => till_Max_Z_tillage, &
                        Date_tillage => till_Date_tillage, Z_tillage => till_Z_tillage, &
                        I_tillage => till_I_tillage, Type_Tillage => till_Type_Tillage, &
                        iType_Tillage => till_iType_Tillage, iTT1 => till_iTT1, iTT2 => till_iTT2, &
                        TAB_Rho_tillage => till_TAB_Rho_tillage, TAB_Rho_cons => till_TAB_Rho_cons, &
                        TAB_K_R_cons => till_TAB_K_R_cons, &
                        TAB_Rho_match => till_TAB_Rho_match, TAB_N_match => till_TAB_N_match, &
                        Rho_tillage => till_Rho_tillage, Rho_cons => till_Rho_cons, &
                        Rho_last => till_Rho_last, K_R_cons => till_K_R_cons, &
                        Rho_match => till_Rho_match, N_match => till_N_match, Slope_match => till_Slope_match, &
                        sumDWC => till_sumDWC, sumAvail1 => till_sumAvail1, sumAvail2 => till_sumAvail2
   
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

   subroutine DoTillage (iTask)
   ! global
   integer, intent(in)                       :: iTask                               ! Task
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
      
      if (allocated(Rho_tillage)) deallocate(Rho_tillage); allocate(Rho_tillage(NumLay))
      if (allocated(Rho_cons))    deallocate(Rho_cons);    allocate(Rho_cons(NumLay))
      if (allocated(Rho_last))    deallocate(Rho_last);    allocate(Rho_last(NumLay))
      if (allocated(K_R_cons))    deallocate(K_R_cons);    allocate(K_R_cons(NumLay))
      if (allocated(Rho_match))   deallocate(Rho_match);   allocate(Rho_match(NumLay))
      if (allocated(N_match))     deallocate(N_match);     allocate(N_match(NumLay))
      if (allocated(Slope_match)) deallocate(Slope_match); allocate(Slope_match(NumLay))

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
      call set_iTill
      
      ! determine number of horizon at depth Max_Z_tillage (MaxNumSoilHo)
      call det_MNSH
    
      ! for special case i_n_model = 3: calculate slope per soil layer (remains constant over time)
      Slope_match = 0.0d0
      
      ! TO ADD: CHECK THAT DEPTH OF EACH TILLAGE EVENT CORRESPONDS TO BOTTOM OF SOIL HORIZON; USER MAY NEED TO DEFINE MULTIPLE SUBS-HORIZONS WITHIN A SINGLE REAL SOIL HORIZON
      ! currently: require changes in horizon number (iSoilLayer) at depth Max_Z_tillage
      fine = .false.
      do i = 1, NumLay
         if (botcom(i) == MaxNumSoilCP) then
            fine = .true.
            exit
         end if
      end do
      if (.not. fine) call fatalerr_collected ('DoTillage', 'Bottom of soil horizon does not coincide with tillage depth(s)')
      
   case (2)
      ! RATE/STATE EVENT
      Rho_last(1:MaxNumSoilHo) = Bdens(1:MaxNumSoilHo)
      
      if (TEST) then
         ! for technical test
         call DTDPST ("YEAR-MONTHST-DAY", t1900, STRNG)
         if (trim(STRNG) == "2016-Apr-05") then
            BDENS(1) = 1000.0d0
            call Change_MvGpars
            Call Adapt_WC_H (TEST)
         else if (trim(STRNG) == "2016-Apr-11") then
            BDENS(1) = 1250.0d0
            call Change_MvGpars
            Call Adapt_WC_H (TEST)
         else if (trim(STRNG) == "2016-Apr-18") then
            BDENS(1) = 1325.0d0        ! no ponding occurs
            !!!BDENS(1) = 1406.322d0   ! in this specific test this change causes ponding
            call Change_MvGpars
            Call Adapt_WC_H (TEST)
         end if
         
      else
         if (Test2) then
            ! for technical test
            call DTDPST ("YEAR-MONTHST-DAY", t1900, STRNG)
            if (trim(STRNG) == "2005-Jun-05") then
               Rho_cons(1) = 1350.0d0
               K_R_cons(1) =  10.0d0
               call Change_MvGpars
               Call Adapt_WC_H (TEST)
            end if
            if (trim(STRNG) == "2005-Oct-30") then
               Rho_cons(1) = 1900.0d0
               K_R_cons(1) =    0.1d0
               call Change_MvGpars
               Call Adapt_WC_H (TEST)
            end if
         end if
         ! normal usage
         if (iTill <= Ntill .and. nint(t1900) == nint(Date_tillage(iTill))) then
            call DTDPST ("YEAR-MONTHST-DAY", t1900, STRNG)
            call Change_Tillage_Info (iTill)
            call Change_Bdens
            iTill = iTill + 1       ! set counter for next tillage event
         else
            call Consolidate_Bdens
         end if
      
         call Change_MvGpars
         
         call DTDPST ("YEAR-MONTHST-DAY", t1900, STRNG)
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

         
         
         Call Adapt_WC_H (TEST)
         
      end if
      

   case (3)
      ! OUTPUT
      if (TEST) then
         call DTDPST ("YEAR-MONTHST-DAY", t1900, STRNG)
         write (222,'(A,F15.5,10(I3,F15.5))') trim(DATE), nraida, (i, Bdens(i), i = 1, MaxNumSoilHo)
         write (224,'(A,10F15.5)') trim(DATE), theta(5), theta(10), theta(20), theta(27), theta(35), nraida, sumDWC, sumAvail1, sumAvail2
         write (226,'(A,10F15.5)') trim(DATE), (CofGen(i,1), i = 1, 10)
      end if
         write (222,'(A,F15.5,10(I3,F15.5))') trim(DATE), nraida, (i, Bdens(i), i = 1, MaxNumSoilHo)
         write (226,'(A,10F15.5)') trim(DATE), (CofGen(i,1), i = 1, 10)
      continue

   case (4)
      ! CLOSURE
      continue
      
   case default
      call fatalerr_collected ('DoTillage','Illegal value for iTask')
   end select

   end subroutine DoTillage

! **************************************************** Change_MvGpars *********************************************************
   subroutine Change_MvGpars
   implicit none
   integer              :: i, node, lay
   integer, parameter   :: Delta = 4
   integer, parameter   :: DeltaMin7 = Delta - 7
   real(8), parameter   :: Omega = -3.97d0
   real(8), parameter   :: Rho_s = 2650d0      ! later as input?
   real(8)              :: wcs_last, Epsilon
   do i = 1 , MaxNumSoilHo
      wcs_last = ParamVG(2,i)       ! help
      
      ! First, fill PARAMVG
      ! wcr
      ParamVG(1,i) = ParamVG(1,i) * Bdens(i)/Rho_last(i)
      ! wcs
      ParamVG(2,i) = ParamVG(2,i) * (Rho_s - Bdens(i))/(Rho_s - Rho_last(i))
      ! ks,fit
      ParamVG(3,i) = ParamVG(3,i) * (ParamVG(2,i)/wcs_last)**3 * (Bdens(i)/Rho_last(i))**DeltaMin7
      ! alpha
      ParamVG(4,i) = ParamVG(4,i) * (Bdens(i)/Rho_last(i))**Omega
      ! lambda: do not change
      !ParamVG(5,i) = ParamVG(5,i)
      ! n
      select case (i_n_model)
      case (1)
         !ParamVG(6,i) = ParamVG(6,i)
         continue
      case(2)
         Epsilon = -0.97d0 + 1.28d0 * psilt(i) / pclay(i)
         ParamVG(6,i) = 1.0d0 + (ParamVG(6,i) - 1.0d0) * (Bdens(i)/Rho_last(i))**Epsilon
      case(3)
         ParamVG(6,i) = dmax1(1.001d0, ParamVG(6,i) + (Bdens(i) - Rho_last(i)) * Slope_match(i))
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
   do node = 1, MaxNumSoilCP
      lay = layer(node)
      CofGen(1:10,node) = ParamVG(1:10,lay)
      ! CofGen(11) and CofGen(12) are not used and not need to be changed
      !CofGen(11,node) = relsatthr(lay)
      !CofGen(12,node) = ksatthr(lay)
   end do
   !!!thetsl(1:numlay) = ParamVG(2,1:numlay)
   
   end subroutine Change_MvGpars
   
! **************************************************** Adapt_WC_H *********************************************************
   subroutine Adapt_WC_H (TEST)
   use soilhydraulics_utils, only: watcon, hconduc, prhead
   implicit none
   
   integer                          :: i
   real(8)                          :: sumWCtmin1, sumWCt, dwc, wcr, wcs, summ, dif
   real(8), dimension(MaxNumSoilCP) :: wc, hold, wcold
   logical                          :: TEST
   
   if (TEST) then
      hold(1:MaxNumSoilCP) = h(1:MaxNumSoilCP)
      wcold(1:MaxNumSoilCP) = theta(1:MaxNumSoilCP)
   end if
   
   select case (iRedist)
   case (0)
      if (.not.TEST) call fatalerr_collected ('Adapt_WC_H', 'Option iRedist = 0 only allowed in combination with TEST option')
      continue
      
   case (1)
      ! keep current wc values (wc_new = wc_old) and only change h; exception: when wc_old > wcs_new: alternative redistribution required
      summ = 0.0d0
      do i = 1, MaxNumSoilCP
         wcs = ParamVG(2,layer(i))
         if (theta(i) < wcs) then
            h(i) = prhead(i, disnod(i), theta(i), CofGen, h)
         else
            summ = summ + (wcs - theta(i))*dz(i)
            theta(i) = wcs
            h(i) = 0.0d0
         end if
      end do
      if (summ > 0.0d0) then
         do i = MaxNumSoilCP, 1, -1
            wcs = ParamVG(2,layer(i))
            dif = wcs - theta(i)
            if (dif > 0.0d0) then
               if (dif < summ) then
                  theta(i) = wcs
                  summ = summ - dif
               else
                  theta(i) = theta(i) + dif
                  summ = 0.0d0
                  exit
               end if
            end if
         end do
      end if
      pond = summ

   case (2)
      sumWCtmin1 = sum(theta(1:MaxNumSoilCP))
      sumWCt = 0.0d0
      sumDWC = 0.0d0
      do i = 1, MaxNumSoilCP
         wc(i) = watcon(i,h(i))
         sumWCt = sumWCt + wc(i)
         dwc = theta(i) - wc(i)
         sumDWC = sumDWC + dwc * dz(i)
      end do
   
      sumAvail1 = 0.0d0
      sumAvail2 = 0.0d0
      if (sumWCt < sumWCtmin1) then
         ! water to be added; same as sumDWC > 0.0
         do i = 1, MaxNumSoilCP
            wcs = ParamVG(2,layer(i))
            sumAvail1 = sumAvail1 + (wcs - wc(i))*dz(i)
         end do
         do i = 1, MaxNumSoilCP
            wcs = ParamVG(2,layer(i))
            if (sumAvail1 > 0.0d0) then
               wc(i) = wc(i) + (wcs - wc(i)) * sumDWC / sumAvail1
               if (wc(i) > wcs) then
                  pond = pond + (wc(i) - wcs) * dz(i)
                  wc(i) = wcs
                  write(333,'(A,I5,F12.4)') Date, i, pond
               end if
            end if
            h(i) = prhead(i, disnod(i), wc(i), CofGen, h)
            theta(i) = wc(i)
         end do
      else if (sumWCt > sumWCtmin1) then
         ! water to be removed; same as sumDWC < 0.0
         do i = 1, MaxNumSoilCP
            wcr = ParamVG(1,layer(i))
            sumAvail2 = SumAvail2 + (wc(i) - wcr) * dz(i)
         end do
         do i = 1, MaxNumSoilCP
            wcr = ParamVG(1,layer(i))
            wc(i) = wc(i) + (wc(i) - wcr) * sumDWC / sumAvail2
            h(i) = prhead(i, disnod(i), wc(i), CofGen, h)
            theta(i) = wc(i)
         end do
         
      endif
   
   end select
   
   if (TEST) then
      do i = 1, MaxNumSoilCP
         wcs = ParamVG(2,layer(i))
         write (444,'(I5,8(A1,F12.6))') i, ',', hold(i), ',', wcold(i), ',', h(i), ',', theta(i), ',', sumDWC, ',', sumAvail1, ',', sumAvail2, ',', theta(i)/wcs
      end do
   end if
write(124,'(A,1P,12E12.5)') Date, Bdens(1), ParamVG(2,layer(1)), theta(1), h(1), hconduc(1,h(1),theta(1),1.0d0), ParamVG(3,layer(1)), Bdens(2), ParamVG(2,layer(2)),theta(2), h(2), hconduc(2,h(2),theta(2),1.0d0), ParamVG(3,layer(2))
   
   end subroutine Adapt_WC_H

   
! **************************************************** Consolidate_Bdens *********************************************************
   subroutine Consolidate_Bdens
   implicit none
   integer :: i

   if (iTill == 1) return        ! in the beginning before first tillage event: do nothing

   forall (i=1:MaxNumSoilHo) Bdens(i) = Rho_cons(i) - (Rho_cons(i) - Rho_last(i)) * dexp(-K_R_cons(i)*nraida*10.0d0)    ! 10: to transform nraida from cm to mm
   write(123,'(A,1P,10E12.5)') Date, nraida, Bdens(1:MaxNumSoilHo)
   end subroutine Consolidate_Bdens
   
! **************************************************** Change_Bdens *********************************************************
   subroutine Change_Bdens
   implicit none
   integer :: i, NumSoilHo
   ! to check: why this loop to determine NumSoilHo?
   NumSoilHo = 1
   do i = 2, NumNod
      if (Z_tillage(iTill) > -zbotcp(i-1) .and. Z_tillage(iTill) <= -zbotcp(i)) then
         NumSoilHo = layer(i)
         exit
      end if
   end do
   forall (i=1:NumSoilHo) Bdens(i) = Rho_last(i) - I_tillage(iTill) * (Rho_last(i) - Rho_tillage(i))
   
   end subroutine Change_Bdens
   
! **************************************************** set_iTill *********************************************************
   subroutine set_iTill
   implicit none
   integer :: i
   iTill = 1
   if (t1900 <= Date_tillage(1)) iTill = 1
   do i = 2, Ntill
      if (Date_tillage(i) < Date_tillage(i-1)) call fatalerr_collected ('set_iTill', 'Dates in tabulated tillage events must be sorted')
      if (t1900 >= Date_tillage(i-1) .and. t1900 < Date_tillage(i-1)) iTill = i-1
   end do
   end subroutine set_iTill

! **************************************************** det_MNSH *********************************************************
   subroutine det_MNSH
   implicit none
   integer :: i
   do i = 2, NumNod
      if (Max_Z_tillage > -zbotcp(i-1) .and. Max_Z_tillage <= -zbotcp(i)) then
         MaxNumSoilHo = layer(i)
         MaxNumSoilCP = i
         exit
      end if
   end do
   end subroutine det_MNSH   

   subroutine Change_Tillage_Info (iTill)
   implicit none
   ! global
   integer, intent(in) :: iTill
   ! local
   integer :: itype, nlay
   
   itype               = Type_Tillage(iTill)
   nlay                = iTT2(itype) - iTT1(itype) + 1
   Rho_tillage(1:nlay) = TAB_Rho_tillage(iTT1(itype):iTT2(itype))
   Rho_cons(1:nlay)    = TAB_Rho_cons(iTT1(itype):iTT2(itype))
   K_R_cons(1:nlay)    = TAB_K_R_cons(iTT1(itype):iTT2(itype))
   Rho_match(1:nlay)   = TAB_Rho_match(iTT1(itype):iTT2(itype))
   N_match(1:nlay)     = TAB_N_match(iTT1(itype):iTT2(itype))
   Slope_match(1:nlay) = (ParamVG(6,1:nlay) - N_match(1:nlay)) / (Rho_cons(1:nlay) - Rho_match(1:nlay))
   end subroutine Change_Tillage_Info

end module tillage_mod
   