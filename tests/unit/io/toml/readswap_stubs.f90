!> No-op stubs for readswap.f90 dependencies that are not needed for the
!! parity test but are referenced by the linker because readswap calls them.
!!
!! Only included in the test executable; never shipped with the production
!! binary.
!!
!! Stubs provided:
!!   writehead        - writes SWAP header to output files (not needed in test)
!!   MeteoInOneFile   - reads meteorology from a single combined file (not
!!                      exercised by the hupselbrook .met yearly-file format)
!!   IrrigationOutput - writes irrigation output files (not opened in test)

subroutine writehead(outf, ftype, filnam, filtext, project)
   implicit none
   integer,          intent(in) :: outf, ftype
   character(len=*), intent(in) :: filnam, filtext, project
end subroutine writehead

subroutine MeteoInOneFile(iTask, ifnd)
   implicit none
   integer, intent(in)  :: iTask
   integer, intent(out) :: ifnd
   ifnd = 0
end subroutine MeteoInOneFile

subroutine IrrigationOutput(task)
   implicit none
   integer, intent(in) :: task
end subroutine IrrigationOutput
