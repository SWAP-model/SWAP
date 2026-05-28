!> Shared helpers for typed CSV record table loaders. (ADR 0044)
!!
!! Centralises the days-since-1900 date arithmetic that was previously
!! duplicated across meteo_csv, csv_reader, readmeteo, and
!! toml_field_helpers.  All implementations used the same
!! Fliegel/Van Flandern Gregorian-to-Julian-Day formula with the
!! epoch constant jd1900 = 2415020 (= Julian Day of 1899-12-31).
!!
!! Public interface
!! ----------------
!!   days_since_1900(year, month, day)  -- forward: Gregorian -> t1900
!!   days1900_to_md(d1900, mm, dd)      -- inverse: t1900 -> (month, day)
!!
!! The inverse subroutine was originally in src/io/readmeteo.f90;
!! the forward function was duplicated in meteo_csv_mod (as
!! days_since_1900), csv_reader_mod (as julian_day + inline offset),
!! and toml_field_helpers_mod (as julian_day + parse_date_to_days1900).
module csv_common_mod
   implicit none
   private

   public :: days_since_1900
   public :: days1900_to_md

   integer, parameter :: jd1900 = 2415020   !! JD of 1899-12-31

contains

   !> Convert a Gregorian calendar date to days since the SWAP epoch
   !! (JD 2415020 = 1899-12-31, so 1900-01-01 -> 1).
   !!
   !! Formula: Fliegel & Van Flandern (1968) proleptic Gregorian calendar.
   pure function days_since_1900(year, month, day) result(d)
      integer, intent(in) :: year, month, day
      integer :: d
      integer :: a, y, m, jd
      a  = (14 - month) / 12
      y  = year + 4800 - a
      m  = month + 12 * a - 3
      jd = day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045
      d  = jd - jd1900
   end function days_since_1900

   !> Convert days-since-jd1900 back to (month, day).
   !! Does not return the year; callers that need the year should
   !! use year_from_t1900 in cropgrass_init_mod or reconstruct via
   !! the forward function.
   !!
   !! Algorithm: inverse Julian Day by the Gregorian calendar
   !! proleptic algorithm (same epoch as days_since_1900).
   subroutine days1900_to_md(d1900, mm, dd)
      integer, intent(in)  :: d1900
      integer, intent(out) :: mm, dd
      integer :: jd, a, b, c, d, e, m
      jd = d1900 + jd1900
      a  = jd + 32044
      b  = (4 * a + 3) / 146097
      c  = a - (146097 * b) / 4
      d  = (4 * c + 3) / 1461
      e  = c - (1461 * d) / 4
      m  = (5 * e + 2) / 153
      dd = e - (153 * m + 2) / 5 + 1
      mm = m + 3 - 12 * (m / 10)
   end subroutine days1900_to_md

end module csv_common_mod
