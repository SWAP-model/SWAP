! file_headers.f90
! IO-OUT/E: standalone output-file header writers, relocated here when
! swapoutput.f90 was deleted (the output_row C-API mechanism was retired).
! These are free-standing (external) subroutines — callers reference them by
! implicit external interface (no `use`), matching their prior home.
!
! writehead    : writes the 5-line SWAP header to formatted/unformatted output
!                files (used by csv_output, cropwofost_runtime, timecontrol).
! WriteSwapOk  : creates the Swap.ok end-of-simulation marker (used by swap_mod).

! ----------------------------------------------------------------------
      subroutine writehead(outf,ftype,filnam,filtext,project)
!-----------------------------------------------------------------------
!      Date    : January 2006
!      Purpose : Writes a header to output files of SWAP

!-----------------------------------------------------------------------
      implicit none

!     global
      integer       outf,ftype
      character(len=*) project,filnam,filtext

!     local
      integer       date_time(6)
      real(8)       dpactualtime
      real(4)       dum
      character(len=80)  model_id,dtstring,String
      character(len=132) Version
!-----------------------------------------------------------------------

! --- Version nr of model, to appear in all output files
!     together with revision nr
      include 'description.fi'

      model_id = 'Swap '//trim(Version)

! --- get actual time
      call dtnow (date_time)
      dum = 0.0
      call dtardp (date_time, dum, dpactualtime)
      call dtdpst ('year-month-day hour:min:sec', dpactualtime,       &
     &             dtstring)

      if (outf.eq.5) then
! ---   write to screen
        write (*,16)  trim(project)
        write (*,17)  trim(filtext)
        write (*,19)  trim(model_id)
        write (*,20)  trim(dtstring)
      else if (ftype.eq.1) then
!       write formatted output file
        write (outf,16)  trim(project)
        write (outf,17)  trim(filtext)
        write (outf,18)  trim(filnam)
        write (outf,19)  trim(model_id)
        write (outf,20)  trim(dtstring)
      else
!       write unformatted output file (header has fixed length of 80 characters)
        write (String,'(a80)')  '* Project:       '//trim(project)
        write (outf)  adjustl(String)
        write (String,'(a80)')  '* File content:  '//trim(filtext)
        write (outf)  adjustl(String)
        write (String,'(a80)')  '* File name:     '//trim(filnam)
        write (outf)  adjustl(String)
        write (String,'(a80)')  '* Model version: '//trim(model_id)
        write (outf)  adjustl(String)
        write (String,'(a80)')  '* Generated at:  '//trim(dtstring)
        write (outf)  adjustl(String)
      endif

 16   format('* Project:       ',a)
 17   format('* File content:  ',a)
 18   format('* File name:     ',a)
 19   format('* Model version: ',a)
 20   format('* Generated at:  ',a)

      return
      end


! ----------------------------------------------------------------------
      subroutine WriteSwapOk(Project)
!-----------------------------------------------------------------------
!      Date    : December 2004
!      Purpose : Create file Swap.ok at end of simulation
!-----------------------------------------------------------------------
      use file_io_mod, only: file_open
      implicit none

      integer   cexf
      character(len=80) project

!     create file Swap.ok to let environment programs verify termination
      call file_open(cexf,'Swap.ok','replace','write')
      call writehead (cexf,1,'Swap.ok',                                 &
     &  'this header only: simulation succesfully terminated',project)
      close (cexf)

      return
      end
