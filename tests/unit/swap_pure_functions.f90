! Extracted pure functions from SWAP for unit testing
! These functions have no module dependencies

! ----------------------------------------------------------------------
real(8) function afgen(table,iltab,x)
! ----------------------------------------------------------------------
!     source             : kees rappoldt, 1/86
!     purpose            : linear interpolation in table with
!                          length iltab for a given value of the
!                          independent variable x.
! ----------------------------------------------------------------------
    implicit none

    integer i,iltab
    real(8) table(iltab),slope,x
! ----------------------------------------------------------------------
    if (table(1).ge.x)  goto 40
    do 10 i = 3,iltab-1,2
      if (table(i).ge.x) goto 30
      if (table(i).lt.table(i-2)) goto 20
10  continue
! --- table fully filled, argument larger then last x in table
    afgen = table(iltab)
    return
! --- table partly filled, argument larger then last x in table
20  afgen = table(i-1)
    return
! --- argument between first and last x in table, interpolation
30  slope = (table(i+1)-table(i-1))/(table(i)-table(i-2))
    afgen = table(i-1) + (x-table(i-2))*slope
    return
! --- argument less or equal to first x in table
40  afgen = table(2)
    return
end
