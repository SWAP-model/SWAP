!> @file checkdate.f90
!! Date-range validator extracted from the deleted src/io/readswap.f90.
!!
!! Verifies that at least one date in the input list falls within the
!! simulation period, or that the simulation period is fully contained
!! within the input-date range. Used by read_ssdi_input (swssdi=1) when
!! the SSDI subsystem is enabled. Other readswap call sites went away
!! with the legacy reader.
!!
!! Original author: April 2006 SWAP team. Behaviour preserved verbatim
!! to keep the swssdi=1 path bit-identical pending ADR 0021 (TOML port
!! of the SSDI block).
      subroutine checkdate(ifnd, dates, tend, tstart, namedat, topic)
      implicit none

      integer   ifnd
      real(8)   dates(ifnd), tend, tstart
      character(len=200) messag
      character(len=5)   namedat
      character(len=*)   topic

      integer i
      logical fldaterr

      fldaterr = .true.
      i = 1
      do while (fldaterr .and. i .le. ifnd)
         if (dates(i) .gt. tstart - 1.d-6 .and. dates(i) .lt. tend + 1.d-6) &
            fldaterr = .false.
         i = i + 1
      enddo
      if (fldaterr) then
         if (dates(1) .lt. tstart + 1.d-6 .and. dates(ifnd) .gt. tend - 1.d-6) &
            fldaterr = .false.
      endif
      if (fldaterr) then
         messag = 'Fatal '//namedat//', no input date within simulation period'
         call fatalerr(topic, messag)
      endif

      return
      end subroutine checkdate
