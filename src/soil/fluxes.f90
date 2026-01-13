! File VersionID:
!   $Id: fluxes.f90 298 2016-07-25 20:07:31Z kroes006 $
! ----------------------------------------------------------------------
      subroutine fluxes ()
! ----------------------------------------------------------------------
!     date               : 29/9/99
!     purpose            : calculates the fluxes between compartments
! ----------------------------------------------------------------------
      use variables, only: q,qbot,dt,inq,numnod,thetm1,theta,dz,qrot,qdra,qimmob,qtop,qrosum,qdrtot,volact,volm1,swbotb,     &
                           FrArMtrx,QExcMpMtx,QMaPo,nrlevs,fllowgwl,qssdi, qssdisum
      implicit none

! --- global

! ----------------------------------------------------------------------
! --- local
      integer i,level

! ----------------------------------------------------------------------

! --- determine qbot if not specified
      if (swbotb .eq. 5 .or. swbotb .eq. 7 .or.                         &
     &    swbotb .eq. 8 .or. swbotb .eq. -2 .or.                        &
     &    (swbotb .eq. 1 .and. fllowgwl)) then
        qbot = qtop + qrosum + qdrtot - QMaPo + (volact-volm1)/dt - qssdisum
      endif

! --- calculate fluxes (cm/d) from changes in volume per compartment
      i = numnod+1
      q(i) = qbot
      inq(i) = inq(i) + q(i)*dt
      do i = numnod,1,-1
        q(i) = - (theta(i)-thetm1(i)+qimmob(i))*FrArMtrx(i)*dz(i)/dt +  &
     &                q(i+1)-qrot(i)+QExcMpMtx(i)+qssdi(i)
     
        do level=1,nrlevs
           q(i) = q(i) - qdra(level,i)
        enddo
        inq(i) = inq(i) + q(i)*dt
      end do

      return
      end

