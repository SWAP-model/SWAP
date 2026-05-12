! File VersionID:
!   $Id: sharedsimulation.f90 362 2018-01-08 13:08:33Z kroes006 $
! ----------------------------------------------------------------------
      subroutine SharedSimulation(task)
      use error_mod, only: fatalerr_collected
      use file_io_mod, only: file_open
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
      integer    unss, ID_Shared, count, IDread
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
        call fatalerr_collected ('readswap',messag)
      endif
      PosArg         = 2
      Call Get_Command_Argument (number=PosArg,value=strIDss)
      read(strIDss,'(i3.3)')ID_Shared
      PosArg         = 3
      Call Get_Command_Argument (number=PosArg,value=strFINA)
      call file_open(unss, strFINA, 'unknown', 'readwrite')
!     open shared data file
!     FromSwap(1) / ToSwap(1): retired stubs — inlined as no-ops (SS-DRV Task 6)
      continue
      continue
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
!     ToSwap(2): retired stub — inlined as no-op (SS-DRV Task 6)
      continue
      return

      case (3)
      rewind(unss)
      write(unss,'(i4)')-1*ID_Shared
!     write New data
!     FromSwap(2): retired stub — inlined as no-op (SS-DRV Task 6)
      continue
      return

      case (4)
! === close Shared Directive file ===========================
!     ToSwap(3) / FromSwap(3): retired stubs — inlined as no-ops (SS-DRV Task 6)
      continue
      continue
      close (unss)

      case default
         call fatalerr_collected ('SharedSimulation', 'Illegal value for TASK')
      end select

      return
      end
