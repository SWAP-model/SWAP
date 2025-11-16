! File VersionID:
!   $Id: watstor.f90 298 2016-07-25 20:07:31Z kroes006 $
! ----------------------------------------------------------------------
      subroutine watstor ()
! ----------------------------------------------------------------------
!     Date               : 29/9/99
!     Purpose            : calc. water storage in the soil profile        
!     Subroutines called : -                                           
!     Functions called   : -                                           
!     File usage         : -                                           
!     Differences SWAP/SWAPS: 
!     - SWAPS: extra parameters
! ----------------------------------------------------------------------
      use variables, only: volm1,volact,numnod,theta,dz,FrArMtrx 
      IMPLICIT NONE

! --- global

! --- local
      INTEGER i
! ----------------------------------------------------------------------
! --- update soil profile water storage
      volm1 = volact
      volact = 0.0d0
      do 10 i = 1,numnod
        volact = volact+theta(i)*dz(i)*FrArMtrx(i)
 10   continue

      return
      end

