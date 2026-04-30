!> Unified CSV reader for SWAP tabular companion files.
!!
!! One public sub: `read_csv_table(path, expected_header, table, errors)`.
!! Validates the file header strict-positional (lowercase, trimmed) against
!! the caller-supplied expected_header. Column 1 is parsed as ISO-date when
!! `expected_header(1) == 'date'`, else as real64. Subsequent columns are
!! always real64. Comment (`#`) and blank lines are skipped. Errors flow
!! through error_collection_t and leave `table` unallocated on failure.
module csv_reader_mod
   use iso_fortran_env, only: real64, iostat_end
   use error_mod, only: error_collection_t,         &
                        ERR_IO_OPEN_FAILED,         &
                        ERR_IO_READ_FAILED,         &
                        ERR_PARSE_TYPE_MISMATCH,    &
                        ERR_PARSE_MISSING_HEADER,   &
                        ERR_PARSE_HEADER_MISMATCH,  &
                        ERR_PARSE_ROW_SHAPE
   implicit none
   private

   public :: read_csv_table

contains

   subroutine read_csv_table(path, expected_header, table, errors)
      character(len=*),          intent(in)    :: path
      character(len=*),          intent(in)    :: expected_header(:)
      real(real64), allocatable, intent(out)   :: table(:,:)
      type(error_collection_t),  intent(inout) :: errors

      integer :: unit, ios, ncols, nrows, irow
      character(len=4096) :: line
      logical :: file_exists, header_done, date_keyed

      ncols = size(expected_header)
      date_keyed = (ncols >= 1) .and. (to_lower(trim(expected_header(1))) == 'date')

      inquire(file=path, exist=file_exists)
      if (.not. file_exists) then
         call errors%append(ERR_IO_OPEN_FAILED, &
            "cannot open CSV: " // trim(path), 'csv_reader')
         return
      end if

      open(newunit=unit, file=path, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         call errors%append(ERR_IO_OPEN_FAILED, &
            "cannot open CSV: " // trim(path), 'csv_reader')
         return
      end if

      ! Pass 1: locate header, count data rows.
      header_done = .false.
      nrows = 0
      do
         read(unit, '(A)', iostat=ios) line
         if (ios == iostat_end) exit
         if (ios /= 0) then
            call errors%append(ERR_IO_READ_FAILED, &
               "read error: " // trim(path), 'csv_reader')
            close(unit)
            return
         end if
         if (is_skippable(line)) cycle
         if (.not. header_done) then
            call validate_header(line, expected_header, path, errors)
            if (errors%has_fatals()) then
               close(unit)
               return
            end if
            header_done = .true.
            cycle
         end if
         nrows = nrows + 1
      end do

      if (.not. header_done) then
         call errors%append(ERR_PARSE_MISSING_HEADER, &
            trim(path) // ": missing header line", 'csv_reader')
         close(unit)
         return
      end if

      allocate(table(nrows, ncols))
      table = 0.0_real64
      if (nrows == 0) then
         close(unit)
         return
      end if

      ! Pass 2: parse data rows.
      rewind(unit)
      header_done = .false.
      irow = 0
      do
         read(unit, '(A)', iostat=ios) line
         if (ios == iostat_end) exit
         if (ios /= 0) then
            call errors%append(ERR_IO_READ_FAILED, &
               "read error: " // trim(path), 'csv_reader')
            close(unit)
            if (allocated(table)) deallocate(table)
            return
         end if
         if (is_skippable(line)) cycle
         if (.not. header_done) then
            header_done = .true.
            cycle
         end if
         irow = irow + 1
         call parse_row(line, irow, ncols, date_keyed, path, table, errors)
         if (errors%has_fatals()) then
            close(unit)
            if (allocated(table)) deallocate(table)
            return
         end if
      end do

      close(unit)
   end subroutine read_csv_table

   pure function is_skippable(line) result(skip)
      character(len=*), intent(in) :: line
      logical :: skip
      character(len=:), allocatable :: trimmed
      skip = .false.
      trimmed = adjustl(line)
      if (len_trim(trimmed) == 0) then
         skip = .true.
         return
      end if
      if (trimmed(1:1) == '#') skip = .true.
   end function is_skippable

   subroutine validate_header(line, expected, path, errors)
      character(len=*),         intent(in)    :: line
      character(len=*),         intent(in)    :: expected(:)
      character(len=*),         intent(in)    :: path
      type(error_collection_t), intent(inout) :: errors
      character(len=:), allocatable :: rest, field
      integer :: i, comma, n
      n = size(expected)
      rest = adjustl(line)
      do i = 1, n
         if (i < n) then
            comma = index(rest, ',')
            if (comma == 0) then
               call errors%append(ERR_PARSE_HEADER_MISMATCH, &
                  trim(path) // ": header column count mismatch", 'csv_reader')
               return
            end if
            field = trim(adjustl(rest(1:comma-1)))
            rest  = adjustl(rest(comma+1:))
         else
            comma = index(rest, ',')
            if (comma /= 0) then
               ! Extra column after the last expected one.
               call errors%append(ERR_PARSE_HEADER_MISMATCH, &
                  trim(path) // ": header has extra columns", 'csv_reader')
               return
            end if
            field = trim(adjustl(rest))
         end if
         if (field /= trim(expected(i))) then
            call errors%append(ERR_PARSE_HEADER_MISMATCH, &
               trim(path) // ": header column '" // field // &
               "' does not match expected '" // trim(expected(i)) // "'", &
               'csv_reader')
            return
         end if
      end do
   end subroutine validate_header

   subroutine parse_row(line, irow, ncols, date_keyed, path, table, errors)
      character(len=*),         intent(in)    :: line
      integer,                  intent(in)    :: irow, ncols
      logical,                  intent(in)    :: date_keyed
      character(len=*),         intent(in)    :: path
      real(real64),             intent(inout) :: table(:,:)
      type(error_collection_t), intent(inout) :: errors
      character(len=:), allocatable :: rest, field
      character(len=64) :: irow_str
      integer :: j, comma, ios
      real(real64) :: val, days
      rest = adjustl(line)
      do j = 1, ncols
         if (j < ncols) then
            comma = index(rest, ',')
            if (comma == 0) then
               write(irow_str, '(I0)') irow
               call errors%append(ERR_PARSE_ROW_SHAPE, &
                  trim(path) // ":row " // trim(irow_str) // &
                  ": expected " // int_str(ncols) // " fields", 'csv_reader')
               return
            end if
            field = trim(adjustl(rest(1:comma-1)))
            rest  = adjustl(rest(comma+1:))
         else
            comma = index(rest, ',')
            if (comma /= 0) then
               write(irow_str, '(I0)') irow
               call errors%append(ERR_PARSE_ROW_SHAPE, &
                  trim(path) // ":row " // trim(irow_str) // &
                  ": expected " // int_str(ncols) // " fields", 'csv_reader')
               return
            end if
            field = trim(adjustl(rest))
         end if

         if (j == 1 .and. date_keyed) then
            if (.not. parse_iso_date(field, days)) then
               write(irow_str, '(I0)') irow
               call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                  trim(path) // ":row " // trim(irow_str) // " col 1: '" // &
                  field // "' not an ISO date", 'csv_reader')
               return
            end if
            table(irow, 1) = days
         else
            read(field, *, iostat=ios) val
            if (ios /= 0) then
               write(irow_str, '(I0)') irow
               call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                  trim(path) // ":row " // trim(irow_str) // " col " // &
                  int_str(j) // ": '" // field // "' not a real number", &
                  'csv_reader')
               return
            end if
            table(irow, j) = val
         end if
      end do
   end subroutine parse_row

   pure function to_lower(s) result(out)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: out
      integer :: i, c
      do i = 1, len(s)
         c = iachar(s(i:i))
         if (c >= iachar('A') .and. c <= iachar('Z')) c = c + 32
         out(i:i) = achar(c)
      end do
   end function to_lower

   pure function int_str(i) result(s)
      integer, intent(in) :: i
      character(len=:), allocatable :: s
      character(len=32) :: buf
      write(buf, '(I0)') i
      s = trim(buf)
   end function int_str

   function parse_iso_date(s, days) result(ok)
      character(len=*), intent(in)  :: s
      real(real64),     intent(out) :: days
      logical :: ok
      integer :: y, m, d, ios, jd, jd1900
      days = 0.0_real64
      ok = .false.
      if (len_trim(s) < 10) return
      if (s(5:5) /= '-' .or. s(8:8) /= '-') return
      read(s(1:4),  '(I4)', iostat=ios) y; if (ios /= 0) return
      read(s(6:7),  '(I2)', iostat=ios) m; if (ios /= 0) return
      read(s(9:10), '(I2)', iostat=ios) d; if (ios /= 0) return
      if (m < 1 .or. m > 12) return
      if (d < 1 .or. d > 31) return
      jd     = julian_day(y, m, d)
      jd1900 = 2415020
      days   = real(jd - jd1900, kind=real64)
      ok = .true.
   end function parse_iso_date

   pure function julian_day(y, m, d) result(jd)
      integer, intent(in) :: y, m, d
      integer :: jd, a, yy, mm
      a  = (14 - m) / 12
      yy = y + 4800 - a
      mm = m + 12 * a - 3
      jd = d + (153 * mm + 2) / 5 + 365 * yy + yy / 4 - yy / 100 + yy / 400 - 32045
   end function julian_day

end module csv_reader_mod
