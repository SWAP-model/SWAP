! to do: default MvG parameters that are read: do they refer to BDENS or Rho_cons; or should BDENS and RhoCons be equal?
! to do: check if Rho_match differs from BDENS or Rho_cons (if not: division by zero possible)

module tillage

   use variables, only: t1900, date, swpfile, swhyst, swsolu, swoxygen, flCropNut, flMacroPore, flksatexm, zbotcp, NumNod, Bdens, layer, nraida, ParamVG, CofGen, &
                        NumLay, pond, theta, h, dz, disnod, botcom, psilt, pclay, SwDiscrvert, tend
                        !, thetsl
   
   implicit none

   ! local (to be saved)
   integer,                            save  :: swtill                              ! switch: 0 = no tillage; 1 = tillage
   integer,                            save  :: Ntill                               ! # of tabulated tillage events
   integer,                            save  :: iTill
   integer,                            save  :: Ntypes                              ! # of tabulated tillage types
   integer,                            save  :: i_n_model                           ! switch for how to treat change in n-parameter: 
                                                                                    ! 1: do not change n; 2: change n based on silt/clay ratio (default); 3: change n based on matching point
   integer,                            save  :: iRedist                             ! type of redistribution after each change in MvG: 0 = not (do not use!); 1 = simpel; 2 = complex (default)
   integer,                            save  :: MaxNumSoilHo, MaxNumSoilCP
   real(8),                            save  :: Max_Z_tillage                       ! max. possible depth of tillage (cm)
   real(8), dimension(:), allocatable, save  :: Date_tillage, Z_tillage, I_tillage  ! tabulated informatino per tillage event
   integer, dimension(:), allocatable, save  :: Type_Tillage

   integer, dimension(:), allocatable, save  :: iType_Tillage, iTT1, iTT2
   real(8), dimension(:), allocatable, save  :: TAB_Rho_tillage, TAB_Rho_cons
   real(8), dimension(:), allocatable, save  :: TAB_K_R_cons
   real(8), dimension(:), allocatable, save  :: TAB_Rho_match, TAB_N_match
   
   real(8), dimension(:), allocatable, save  :: Rho_tillage, Rho_cons, Rho_last
   real(8), dimension(:), allocatable, save  :: K_R_cons
   real(8), dimension(:), allocatable, save  :: Rho_match, N_match, Slope_match

   real(8),                            save  :: sumDWC, sumAvail1, sumAvail2

!  by default: all in this module is private (local)
   private
!  except for these public routines/functions
   public :: DoTillage
!  and except for these public variables
   public :: swtill
   
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
      
      ! read input data
      call Read_Tillage
      
      ! no tillage required; leave DoTillage immediately
      if (swtill == 0) return
      
      ! some checks: some combinations not (yet) allowed
      if (swtill == 1) then
         if (swhyst == 1)      call fatalerr ('DoTillage', 'swhyst = 1 not allowed')
         if (swsolu == 1)      call fatalerr ('DoTillage', 'swsolu = 1 not (yet) allowed')
         if (swoxygen == 2)    call fatalerr ('DoTillage', 'swoxygen = 2 not (yet) allowed')
         if (flCropNut)        call fatalerr ('DoTillage', 'flCropNut = 1 not (yet) allowed')
         if (flMacroPore)      call fatalerr ('DoTillage', 'swmacro = 1 not (yet) allowed')
         if (flksatexm)        call fatalerr ('DoTillage', 'flksatexm not (yet) allowed')
         if (SwDiscrvert == 1) call fatalerr ('DoTillage', 'SwDiscrvert = 1 not (yet) allowed')
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
      if (.not. fine) call fatalerr ('DoTillage', 'Bottom of soil horizon does not coincide with tillage depth(s)')
      
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
      call fatalerr ('DoTillage','Illegal value for iTask')
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
   implicit none
   
   integer                          :: i
   real(8)                          :: sumWCtmin1, sumWCt, dwc, wcr, wcs, summ, dif
   real(8), dimension(MaxNumSoilCP) :: wc, hold, wcold
   real(8)                          :: watcon, prhead, hconduc   ! functions called
   logical                          :: TEST
   
   if (TEST) then
      hold(1:MaxNumSoilCP) = h(1:MaxNumSoilCP)
      wcold(1:MaxNumSoilCP) = theta(1:MaxNumSoilCP)
   end if
   
   select case (iRedist)
   case (0)
      if (.not.TEST) call fatalerr ('Adapt_WC_H', 'Option iRedist = 0 only allowed in combination with TEST option')
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
      if (Date_tillage(i) < Date_tillage(i-1)) call fatalerr ('set_iTill', 'Dates in tabulated tillage events must be sorted')
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

! **************************************************** Read_Tillage *********************************************************
   subroutine Read_Tillage
   implicit none
   ! local (not to be saved)
   integer     :: IunIn, tmin, tmax, i, j
   ! functions
   integer     :: getun
   logical     :: RDinqr

   Max_Z_tillage = 0.0d0
   IunIn = getun (200, 900)
   call RDinit (IunIn, 0, swpfile)
      swtill = 0                                                     ! default: do not consider tillage
      if (RDinqr('swtill')) call RDsinr ('swtill', 0, 1, swtill)
      
      if (swtill == 1) then
         ! switch for how to treat change in n-parameter:
         i_n_model = 2
         if (RDinqr('i_n_model')) call RDsinr ('i_n_model', 1, 3, i_n_model)
         
         ! type of redistribution
         iRedist = 2 ! default
         if (RDinqr('iRedist')) call RDsinr ('iRedist', 0, 2, iRedist)
         
         ! get length of tabulated input data (number of tillage times, Ntill)
         call RDinne ('Date_tillage', Ntill)
         if (Ntill < 1) call fatalerr ('Read_Tillage', 'You must supply tabulated data for tillage events')
         
         ! allocate arrays for tabulated input
         if (allocated(Date_tillage)) deallocate(Date_tillage); allocate (Date_tillage(Ntill+1))
         if (allocated(Z_tillage))    deallocate(Z_tillage);    allocate (Z_tillage(Ntill))
         if (allocated(I_tillage))    deallocate(I_tillage);    allocate (I_tillage(Ntill))
         if (allocated(Type_tillage)) deallocate(Type_tillage); allocate (Type_tillage(Ntill))
         
         call RDftim ('Date_tillage', Date_tillage, Ntill, Ntill); Date_Tillage(Ntill+1) = tend + 1.0d0
         call RDfdor ('Z_tillage',    0.0d0, -zbotcp(NumNod), Z_tillage, Ntill, Ntill)   ! still to check: depth must be equal to bottom of a soil horizon
         call RDfdor ('I_tillage',    0.0d0, 1.0d0, I_tillage, Ntill, Ntill)
         call RDfint ('Type_tillage', Type_tillage, Ntill, Ntill)
         tmin = minval(Type_tillage); if (tmin < 1) call fatalerr ('Read_Tillage','Type_Tillage should be > 0')
         tmax = maxval(Type_tillage); if (tmax < 1) call fatalerr ('Read_Tillage','Type_Tillage should be > 0')

         ! length of second tabel entries: Ntypes
         call RDinne ('iType_Tillage', Ntypes)
         if (Ntypes < 1) call fatalerr ('Read_Tillage', 'You must supply tabulated data for tillage events (Ntypes)')

         ! allocate arrays for tabulated input
         if (allocated(iType_Tillage))    deallocate(iType_Tillage);    allocate(iType_Tillage(Ntypes))
         if (allocated(TAB_Rho_cons))     deallocate(TAB_Rho_cons);     allocate(TAB_Rho_cons(Ntypes))
         if (allocated(TAB_Rho_tillage))  deallocate(TAB_Rho_tillage);  allocate(TAB_Rho_tillage(Ntypes))
         if (allocated(TAB_K_R_cons))     deallocate(TAB_K_R_cons);     allocate(TAB_K_R_cons(Ntypes))
         if (allocated(TAB_Rho_match))    deallocate(TAB_Rho_match);    allocate(TAB_Rho_match(Ntypes))
         if (allocated(TAB_N_match))      deallocate(TAB_N_match);      allocate(TAB_N_match(Ntypes))
 
         call RDfinr ('iType_Tillage', tmin,     tmax,     iType_Tillage,     Ntypes, Ntypes)
         call RDfdor ('Rho_cons',      100.0d0,  3000.0d0, TAB_Rho_cons,      Ntypes, Ntypes)
         call RDfdor ('Rho_tillage',   100.0d0,  3000.0d0, TAB_Rho_tillage,   Ntypes, Ntypes)
         call RDfdor ('k_R',             1.0d-4,   10.0d0, TAB_K_R_cons,      Ntypes, Ntypes)
         
         if (i_n_model == 3) then
            call RDfdor('Rho_match', 100.0d0, 3000.0d0, TAB_Rho_match, Ntypes, Ntypes)
            call RDfdor('N_match',   1.001d0,   10.0d0, TAB_N_match,   Ntypes, Ntypes)
         end if
         
         ! store first and last position per tillage type in vector iType_Tillage
         ! assumes consitent supply of input: e.g., tmin = 1, tmax > 1 and all intermediate values are present
         if (allocated(iTT1)) deallocate(iTT1); allocate(iTT1(Ntill)); iTT1 = 0
         if (allocated(iTT2)) deallocate(iTT2); allocate(iTT2(Ntill)); iTT2 = 0
         do j = 1, Ntill
            do i = 1, Ntypes
               if (iTT1(j) == 0 .and. iType_Tillage(i) == j) iTT1(j) = i
               if (iTT1(j) >  0 .and. iType_Tillage(i) == j) iTT2(j) = i
            end do
! check if NumLay is exceeded; and if iTT2-iTT1+1 corresponds with number of layer within Z_tillage
         end do
         
         ! special case: consolidation according to time
         !use_K_T = RDinqr('k_T')
         !if (use_K_T) call RDfdor ('k_T', 1.0d-4, 10.0d0, K_T_cons, maho, NumLay)
         
         Max_Z_tillage = maxval(Z_tillage(1:Ntill))
         
      end if
   close (IunIn)
   end subroutine Read_Tillage
   
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
   
end module tillage
   