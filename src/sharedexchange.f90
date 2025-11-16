! File VersionID:
!   $Id: sharedexchange.f90 334 2017-07-06 10:39:10Z kroes006 $
! ----------------------------------------------------------------------
      subroutine FromSwap(task) 
! ----------------------------------------------------------------------
!     UpDate             : July 2017
!     Date               : July 2009
!     Purpose            : write data to shared files
! ----------------------------------------------------------------------
!      use Variables
      implicit none

! --- global variables ------------------
      integer task
! --- local variables ------------------
!      save  unss     
! ----------------------------------------------------------------------

      select case (task)
      case (1)
      continue
!     section meant to initialize, open or access exchange file

      return

      case (2)
      continue
!     section meant to process intermediate SWAP results (variables/parameters) and writing to the exchange file

      return

      case (3)
      continue
!     section meant to terminate the exchange process and close the exchange file

      case default
         call fatalerr ('FromSWAP', 'Illegal value for TASK')
      end select

      return
      end
! ----------------------------------------------------------------------
      subroutine ToSwap(task) 
! ----------------------------------------------------------------------
!     UpDate             : July 2017
!     Date               : July 2009
!     Purpose            : write data to shared files
! ----------------------------------------------------------------------
!      use Variables
      implicit none

! --- global variables ------------------
      integer task
! --- local variables ------------------
!      save  unss     
! ----------------------------------------------------------------------

      select case (task)
      case (1)
      continue
!     section meant to initialize, open or access exchange file

      return

      case (2)
      continue
!     section meant to read from the exchange file and process the info to SWAP variables/parameters

      return

      case (3)
      continue
!     section meant to terminate the exchnge process and close the exchange file

      case default
         call fatalerr ('ToSWAP', 'Illegal value for TASK')
      end select

      return
      end
