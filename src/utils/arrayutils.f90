module array_utils
   use error_mod, only: fatalerr_collected
    !> Module containing utility functions for interpolation and table lookups
    !!
    !! @author Original SWAP team
    !! @date February 2026 (modularization)
    use iso_fortran_env, only: real64
    implicit none
    private
    public :: afgen, stepnr, insw, interpol

contains
    ! ----------------------------------------------------------------------
   function afgen(table, iltab, x)
      !> Linear interpolation in table with length iltab for a given value of the independent variable x.
      implicit none
      
      ! Arguments
      integer, intent(in) :: iltab
      real(real64), intent(in) :: table(iltab)
      real(real64), intent(in) :: x
      real(real64) :: afgen
      
      ! Local variables
      integer :: i
      real(real64) :: slope
      
      if (x <= table(1)) then
         ! Argument less or equal to first x in table
         afgen = table(2)
         return
      end if
      
      do i = 3, iltab-1, 2
         if (table(i) >= x) then
            ! Argument between first and last x in table, interpolation
            slope = (table(i+1) - table(i-1)) / (table(i) - table(i-2))
            afgen = table(i-1) + (x - table(i-2)) * slope
            return
         end if
         
         if (table(i) < table(i-2)) then
            ! Table partly filled, argument larger than last x in table
            afgen = table(i-1)
            return
         end if
      end do
      
      ! Table fully filled, argument larger than last x in table
      afgen = table(iltab)
      
   end function afgen

    ! ----------------------------------------------------------------------

   function stepnr(array, length, x)
      !> Find step number in array
      !! @warning
      !! This function is not called from any SWAP code.
      !!
      implicit none
      
      ! Arguments
      integer, intent(in) :: length
      real(real64), intent(in) :: array(length)
      real(real64), intent(in) :: x
      integer :: stepnr
      
      ! Local variables
      integer :: i
      
      if (x <= array(1)) then
         ! Argument less or equal to first x in array
         stepnr = 1
         return
      end if
      
      do i = 2, length
         if (array(i) > x .or. array(i) < array(i-1)) then
            ! Array partly filled, argument larger than last x in array or
            ! argument between first and last x in table
            stepnr = i - 1
            return
         end if
      end do
      
      ! Array fully filled, argument larger than last x in array
      stepnr = length
      
   end function stepnr

    ! ----------------------------------------------------------------------
   function insw(x1, x2, x3)
      !> Switch routine taken from TTUTIL, returns x2 if x1 < 0.0, otherwise returns x3
      !! Called from: cropgrowth.f90
      implicit none
      
      ! Arguments
      real(real64), intent(in) :: x1
      real(real64), intent(in) :: x2
      real(real64), intent(in) :: x3
      real(real64) :: insw
      
      if (x1 < 0.0_real64) then
         insw = x2
      else
         insw = x3
      end if
      
   end function insw


    ! ----------------------------------------------------------------------
   function interpol(mn, mx, x)
      !> Interpolation routine taken from TTUTIL, returns x if x is between mn and mx, otherwise returns the nearest bound
      !!
      !! Called from: cropgrowth.f90
      implicit none
      
      ! Arguments
      real(real64), intent(in) :: mn
      real(real64), intent(in) :: mx
      real(real64), intent(in) :: x
      real(real64) :: interpol
      
      ! Local variables
      character(len=52) :: messag
      
      if (mx < mn) then
         ! Minimum is larger than maximum, should generate an error
         write(messag, '(2(a,g12.5))') 'argument error, min = ', mn, ', max = ', mx
         call fatalerr_collected('limit', messag)
      end if
      
      if (x < mn) then
         ! x below allowed range; return lower bound
         interpol = mn
      else if (x <= mx) then
         ! x in range; return x
         interpol = x
      else
         ! x above allowed range; return upper bound
         interpol = mx
      end if
      
   end function interpol

end module array_utils
! ----------------------------------------------------------------------