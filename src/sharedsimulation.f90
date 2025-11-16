! File VersionID:
!   $Id: sharedsimulation.f90 362 2018-01-08 13:08:33Z kroes006 $
! ----------------------------------------------------------------------
      subroutine SharedSimulation(task)
! ----------------------------------------------------------------------
!     UpDate             : July 2017
!     Date               : July 2009
!     Purpose            : open and write data to shared files
! ----------------------------------------------------------------------

!      use Variables
      implicit none

! --- global variables ------------------
      integer task
! --- local variables ------------------
      integer PosArg
      integer(4) delay
      logical    flhold
      integer    unss, ID_Shared, getun, count, IDread
      character(len=3)   strIDss
      character(len=80)  strFINA
      character(len=200) messag

      save  unss, ID_Shared
! ----------------------------------------------------------------------
!      data  delay/1/                       ! delay in secs
      data  delay/100/                       ! delay in msecs


      select case (task)

      case (1)
      count = Command_Argument_Count()
      if(count.ne.3) then
        messag = 'Argument of executable-call is not correct for '//    &
     &                    ' Shared simulation !'
        call fatalerr ('readswap',messag)
      endif
      PosArg         = 2
      Call Get_Command_Argument (number=PosArg,value=strIDss)
      read(strIDss,'(i3.3)')ID_Shared
      PosArg         = 3
      Call Get_Command_Argument (number=PosArg,value=strFINA)
      unss  = getun (20,99)
      open(unit=unss,file=strFINA,status='unknown',                     &
     &     action='READWRITE') !,share='DENYNONE')
!     open shared data file
      call FromSwap(task)
      call ToSwap(task)
      return

      case (2)
      flhold = .true.
      do while (flhold)
!        call sleepqq(delay)             ! delay in milisecs
        call sleep(delay)             ! delay in secs
        rewind(unss)
        read(unss,'(i4)')IDread
        if(IDread.eq.ID_Shared) flhold=.false.
      end do
!     read New data
      call ToSwap(2)
      return

      case (3)
      rewind(unss)
      write(unss,'(i4)')-1*ID_Shared
!     write New data
      call FromSwap(2)
      return

      case (4)
! === close Shared Directive file ===========================
      call ToSwap(3)
      call FromSwap(3)
      close (unss)

      case default
         call fatalerr ('SharedSimulation', 'Illegal value for TASK')
      end select

      return
      end
