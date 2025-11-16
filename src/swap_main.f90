! swap_main.f90
! Marius Heinen, March 2021, WENR
!
! As of version 4.1.75, the original main (swap.f90) is split up in a main caller (this file), and
! a subroutine swap (swap.f90), which is the actual model.
! The reason is that all code, except the swap_main, can then also be used for
! preparing a DLL-version of swap, which then can be used by others.
! The model has been split into three tasks:
!    1 - initialization
!    2 - dynamic (time loop)
!    3 - closure
! Therefore, the model will be called three times in a row.

! ----------------------------------------------------------------------
program swap_main
! ----------------------------------------------------------------------

use variables, only: logf
implicit none

! because subroutine swap has optional arguments, we must define the interface
interface
   subroutine swap(iCaller, iTask, toswap, fromswap)
      use swap_exchange
      integer,           intent(in)              :: iCaller, iTask
      type(swap_input),  intent(in),    optional :: toswap
      type(swap_output), intent(out),   optional :: fromswap
   end subroutine swap
end interface

! local
integer              :: iTask
integer, parameter   :: iCaller = 0                   ! Who is calling swap? 0 = swap-main; > 0 external model is calling swap as DLL (iCaller > 0 for testing only)
integer, save        :: iset, insets, iun1, iun2

! functions
integer              :: getun

! open logfile and read rerun file
iun1 = getun (400,900)
iun2 = getun (iun1+1,900)
call fopens (iun1,'reruns.log','new','del')
call rdsets (iun2,iun1,'reruns.dat',insets)
if (insets == 0) write (iun1,'(a)') 'No reruns defined.'

! reruns (if supplied; else this loop is performed only once)
do iset = 0, insets

! select rerun set
   call rdfrom (iset,.true.)

!  Initialize swap
   iTask = 1
   if (iCaller == 0) call swap(iCaller, iTask)
   if (iCaller /= 0) call dummy(iTask)

!  Dynamic call to swap
   iTask = 2
   if (iCaller == 0) call swap(iCaller, iTask)
   if (iCaller /= 0) call dummy(iTask)

!  Close swap
   iTask = 3
   if (iCaller == 0) call swap(iCaller, iTask)
   if (iCaller /= 0) call dummy(iTask)

end do
close (iun2)

! write message on screen
write(*,'(a)')' Swap normal completion!'
close(logf)
call CloseTempFil
Call Exit(100)

   contains
   subroutine dummy(iTask)
   use swap_exchange
   implicit none
   integer, intent(in)  :: iTask
   integer              :: i, itel
   integer, save        :: itstart, itend
   type(swap_input)     :: testin
   type(swap_output)    :: testout

   select case (iTask)
   case(1)
      testin%tstart = 0.0d0; testin%tend = 0.0d0
      testin%rain = 0.0d0; testin%tmin = 0.0d0; testin%tmax = 0.0d0; testin%rad = 0.0d0; testin%hum = 0.0d0; testin%wind = 0.0d0; testin%etref = 0.0d0; testin%wet = 0.0d0
      testin%icrop = 0; testin%lai = 0.0d0; testin%ch = 0.0d0; testin%zroot = 0.0d0
      call swap(iCaller, iTask, toswap = testin, fromswap = testout)
      if (testout%ierrorcode /= 0) call fatalerr('swap_main', 'error')
      itstart = nint(testout%tstart + 365.0d0)
      itend   = nint(testout%tend   + 365.0d0)

   case(2)
      itel = 0
      do i = itstart, itend, 1
         itel = itel + 1
         testin%tstart = dble(i)
         testin%tend   = dble(i)    ! very important: if for example one day is to be considered then Tend must be equal to Tstart, since in TimeControl timing is based on Tend+1
         testin%rain   =   1.0d0    ! mm/d
         testin%tmin   =  15.0d0    ! deg. C
         testin%tmax   =  25.0d0    ! deg. C
         testin%rad    =   1.0d4    ! kJ/m2/d
         testin%hum    =   0.75d0   ! kPa
         testin%wind   =   5.0d0    ! m/s
         testin%etref  = -99.9d0    ! mm/d
         testin%wet    =   1.0d-1   ! [0...1]
         if (itel > 100 .and. itel < 275) then
            testin%icrop  =   1
            testin%lai    =   2.0d0
            testin%ch     =  40.0d0
            testin%zroot  =  25.0d0
         else
            testin%icrop  =   0
            testin%lai    =   0.0d0
            testin%ch     =   0.0d0
            testin%zroot  =   0.0d0
         end if
         call swap(iCaller, iTask, toswap = testin, fromswap = testout)
         if (testout%ierrorcode /= 0) call fatalerr('swap_main', 'error')
         !write (*,'(i10, f8.3,i10)') i, testout%tact/testout%tpot, testout%ierrorcode
      end do

   case(3)
      call swap(iCaller, iTask)
   end select

   end subroutine dummy

end program swap_main

