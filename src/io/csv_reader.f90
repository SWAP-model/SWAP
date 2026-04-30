!> Minimal CSV reader for SWAP "long-table" companion files.
!!
!! The TOML pipeline keeps short tables (handful of rows) inline in
!! `swap.toml`, but cases with hundreds of rows (fixed irrigation,
!! prescribed groundwater levels, met series) author them as CSV
!! companion files alongside `swap.toml`. This module is the single
!! reader for those CSVs.
!!
!! Conventions (kept minimal on purpose):
!!   * One public subroutine: `read_csv_date_reals`.
!!   * Field separator: comma. No TSV, no quoted strings, no escapes.
!!   * Column 1 is an ISO-like date `YYYY-MM-DD`; converted to
!!     days-since-1900 (the legacy time axis) using the same
!!     algorithm as `parse_date_to_days1900`.
!!   * Columns 2..(1+ncols_expected) are real(real64).
!!   * Lines starting with `#` are CSV comments and are skipped.
!!   * A header row whose first field is non-numeric, non-date is
!!     also skipped (one allowed; subsequent header-like rows
!!     append a parse error).
!!   * Errors flow through `error_collection_t` — never `write` /
!!     `stop` from this module.
!!
!! Output `table(:,:)` shape is `(nrows, 1 + ncols_expected)` where
!! col 1 is days-since-1900 and col 2.. are the real values, mirroring
!! the layout of `irrigation_config_t%fixed_events` and
!! `bottom_boundary_config_t%gwl_table` (so the adapter unpacks the
!! same way for inline-TOML and CSV-companion paths).
module csv_reader_mod
   use iso_fortran_env, only: real64, iostat_end
   use error_mod, only: error_collection_t,        &
                        ERR_IO_READ_FAILED,        &
                        ERR_PARSE_TYPE_MISMATCH
   implicit none
   private

   public :: read_csv_date_reals

contains

   !> Read a CSV file with one date column followed by `ncols_expected`
   !! real columns. On success `table` is allocated to
   !! `(nrows, 1 + ncols_expected)` with col 1 = days-since-1900.
   !!
   !! Errors append to `errors` and leave `table` unallocated.
   subroutine read_csv_date_reals(path, ncols_expected, table, errors)
      character(len=*),          intent(in)    :: path
      integer,                   intent(in)    :: ncols_expected
      real(real64), allocatable, intent(out)   :: table(:,:)
      type(error_collection_t),  intent(inout) :: errors

      integer :: unit, ios, nrows, irow, header_seen
      character(len=4096) :: line
      logical :: file_exists

      ! Guard: file must exist.
      inquire(file=path, exist=file_exists)
      if (.not. file_exists) then
         call errors%append(ERR_IO_READ_FAILED, &
            "CSV file not found: " // trim(path), 'csv_reader')
         return
      end if

      open(newunit=unit, file=path, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         call errors%append(ERR_IO_READ_FAILED, &
            "open failed for CSV file: " // trim(path), 'csv_reader')
         return
      end if

      ! First pass: count data rows (skip blanks / comments / header).
      nrows = 0
      header_seen = 0
      do
         read(unit, '(A)', iostat=ios) line
         if (ios == iostat_end) exit
         if (ios /= 0) then
            call errors%append(ERR_IO_READ_FAILED, &
               "read error scanning CSV: " // trim(path), 'csv_reader')
            close(unit)
            return
         end if
         if (is_skippable(line, header_seen)) cycle
         nrows = nrows + 1
      end do

      if (nrows == 0) then
         allocate(table(0, 1 + ncols_expected))
         close(unit)
         return
      end if

      allocate(table(nrows, 1 + ncols_expected))
      table = 0.0_real64

      ! Second pass: parse data rows.
      rewind(unit)
      irow = 0
      header_seen = 0
      do
         read(unit, '(A)', iostat=ios) line
         if (ios == iostat_end) exit
         if (ios /= 0) then
            call errors%append(ERR_IO_READ_FAILED, &
               "read error parsing CSV: " // trim(path), 'csv_reader')
            close(unit)
            if (allocated(table)) deallocate(table)
            return
         end if
         if (is_skippable(line, header_seen)) cycle
         irow = irow + 1
         call parse_data_row(line, irow, ncols_expected, table, path, errors)
         if (errors%has_fatals()) then
            close(unit)
            if (allocated(table)) deallocate(table)
            return
         end if
      end do

      close(unit)
   end subroutine read_csv_date_reals

   !> Decide whether a raw line should be skipped (blank / `#` comment /
   !! one-shot header row). `header_seen` is incremented when a header
   !! row is consumed so a second non-numeric row triggers a parse error
   !! at the caller (caught by parse_data_row).
   function is_skippable(line, header_seen) result(skip)
      character(len=*), intent(in)    :: line
      integer,          intent(inout) :: header_seen
      logical :: skip
      character(len=:), allocatable   :: trimmed, first_field

      skip = .false.
      trimmed = adjustl(line)
      if (len_trim(trimmed) == 0) then
         skip = .true.
         return
      end if
      if (trimmed(1:1) == '#') then
         skip = .true.
         return
      end if
      ! One-shot header heuristic: first field neither parses as date nor
      ! as real. Only allowed once (the first such row).
      first_field = leading_field(trimmed)
      if (.not. looks_like_date(first_field) .and. &
          .not. looks_like_real(first_field)) then
         if (header_seen == 0) then
            header_seen = 1
            skip = .true.
            return
         end if
      end if
   end function is_skippable

   !> Extract the substring up to (but not including) the first comma.
   !! Trims leading/trailing whitespace.
   function leading_field(s) result(field)
      character(len=*), intent(in) :: s
      character(len=:), allocatable :: field
      integer :: comma
      comma = index(s, ',')
      if (comma == 0) then
         field = trim(adjustl(s))
      else
         field = trim(adjustl(s(1:comma-1)))
      end if
   end function leading_field

   !> Heuristic: looks like `YYYY-MM-DD` (10 chars, dashes at positions
   !! 5 and 8, digits elsewhere). Not a full validator — the real parse
   !! happens in parse_iso_date.
   function looks_like_date(s) result(yes)
      character(len=*), intent(in) :: s
      logical :: yes
      integer :: i
      yes = .false.
      if (len_trim(s) < 10) return
      if (s(5:5) /= '-' .or. s(8:8) /= '-') return
      do i = 1, 10
         if (i == 5 .or. i == 8) cycle
         if (s(i:i) < '0' .or. s(i:i) > '9') return
      end do
      yes = .true.
   end function looks_like_date

   !> Heuristic: starts with a digit, sign, or decimal point. Good
   !! enough for screening header strings; the actual parse uses
   !! `read(... , *)` which is the source of truth.
   function looks_like_real(s) result(yes)
      character(len=*), intent(in) :: s
      logical :: yes
      character :: c
      yes = .false.
      if (len_trim(s) == 0) return
      c = s(1:1)
      if ((c >= '0' .and. c <= '9') .or. c == '-' .or. c == '+' .or. c == '.') yes = .true.
   end function looks_like_real

   !> Parse a single comma-separated data row into row `irow` of `table`.
   !! Appends a fatal error on any malformed cell.
   subroutine parse_data_row(line, irow, ncols_expected, table, path, errors)
      character(len=*),         intent(in)    :: line
      integer,                  intent(in)    :: irow
      integer,                  intent(in)    :: ncols_expected
      real(real64),             intent(inout) :: table(:,:)
      character(len=*),         intent(in)    :: path
      type(error_collection_t), intent(inout) :: errors

      character(len=:), allocatable :: rest, field
      character(len=64) :: irow_str
      integer :: comma, j, ios
      real(real64) :: val
      real(real64) :: days

      rest = adjustl(line)
      ! Field 1: date.
      comma = index(rest, ',')
      if (comma == 0) then
         write(irow_str, '(I0)') irow
         call errors%append(ERR_PARSE_TYPE_MISMATCH, &
            "row " // trim(irow_str) // ": expected at least one comma in " // &
            trim(path), 'csv_reader')
         return
      end if
      field = trim(adjustl(rest(1:comma-1)))
      if (.not. parse_iso_date(field, days)) then
         write(irow_str, '(I0)') irow
         call errors%append(ERR_PARSE_TYPE_MISMATCH, &
            "row " // trim(irow_str) // ": malformed date '" // field // "' in " // &
            trim(path), 'csv_reader')
         return
      end if
      table(irow, 1) = days
      rest = adjustl(rest(comma+1:))

      ! Fields 2..(1+ncols_expected): reals.
      do j = 1, ncols_expected
         if (j < ncols_expected) then
            comma = index(rest, ',')
            if (comma == 0) then
               write(irow_str, '(I0)') irow
               call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                  "row " // trim(irow_str) // ": missing column in " // trim(path), &
                  'csv_reader')
               return
            end if
            field = trim(adjustl(rest(1:comma-1)))
            rest  = adjustl(rest(comma+1:))
         else
            ! Last column: take the rest of the line (may have trailing
            ! comma to drop — strip it).
            comma = index(rest, ',')
            if (comma == 0) then
               field = trim(adjustl(rest))
            else
               field = trim(adjustl(rest(1:comma-1)))
            end if
         end if
         read(field, *, iostat=ios) val
         if (ios /= 0) then
            write(irow_str, '(I0)') irow
            call errors%append(ERR_PARSE_TYPE_MISMATCH, &
               "row " // trim(irow_str) // ": non-numeric cell '" // field // &
               "' in " // trim(path), 'csv_reader')
            return
         end if
         table(irow, 1 + j) = val
      end do
   end subroutine parse_data_row

   !> Parse `YYYY-MM-DD` -> days-since-1900 (1900-01-01 == 0).
   !! Returns .false. on any parse failure.
   function parse_iso_date(s, days) result(ok)
      character(len=*), intent(in)  :: s
      real(real64),     intent(out) :: days
      logical :: ok
      integer :: y, m, d, ios
      integer :: jd, jd1900

      days = 0.0_real64
      ok = .false.
      if (len_trim(s) < 10) return
      if (s(5:5) /= '-' .or. s(8:8) /= '-') return
      read(s(1:4), '(I4)', iostat=ios)  y; if (ios /= 0) return
      read(s(6:7), '(I2)', iostat=ios)  m; if (ios /= 0) return
      read(s(9:10), '(I2)', iostat=ios) d; if (ios /= 0) return
      if (m < 1 .or. m > 12) return
      if (d < 1 .or. d > 31) return
      jd     = julian_day(y, m, d)
      jd1900 = 2415020
      days   = real(jd - jd1900, kind=real64)
      ok = .true.
   end function parse_iso_date

   !> Same Julian-day algorithm used by toml_field_helpers_mod's private
   !! helper. Inlined here to keep this module independent of tomlf.
   pure function julian_day(y, m, d) result(jd)
      integer, intent(in) :: y, m, d
      integer :: jd, a, yy, mm
      a  = (14 - m) / 12
      yy = y + 4800 - a
      mm = m + 12 * a - 3
      jd = d + (153 * mm + 2) / 5 + 365 * yy + yy / 4 - yy / 100 + yy / 400 - 32045
   end function julian_day

end module csv_reader_mod
