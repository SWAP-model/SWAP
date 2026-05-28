!> Typed CSV record tables for meteorological forcing files.
!!
!! Three families:
!!   * meteo_daily_table_t   — 9-column daily meteo (ADR 0014 canonical schema)
!!   * meteo_detail_table_t  — 7-column sub-daily meteo (swmetdetail=1)
!!   * rain_events_table_t   — 2-column rain events (swrain=3)
!!
!! Each table owns its schema, a `load(path, errors)` method that wraps the
!! generic csv_reader, an inline validate step, and a year-window slice
!! helper. The `is_loaded` flag distinguishes "not requested" from "empty".
!!
!! Date column is real(real64) days-since-JD-1900 (csv_reader's convention).
!! A future arc may introduce a typed date_t; out of scope here.
module meteo_csv_mod
   use iso_fortran_env, only: real64
   use error_mod,       only: error_collection_t
   implicit none
   private

   public :: meteo_daily_table_t
   public :: meteo_detail_table_t
   public :: rain_events_table_t
   public :: meteo_daily_row_t
   public :: meteo_detail_row_t
   public :: rain_event_row_t

   !> One daily meteo record. Units match csv_reader output (no scaling).
   type :: meteo_daily_row_t
      real(real64) :: date  = 0.0_real64  !! days since JD 1900
      real(real64) :: rad   = 0.0_real64  !! kJ m-2 d-1 (NOT yet J/m2/d)
      real(real64) :: tmin  = 0.0_real64  !! degC
      real(real64) :: tmax  = 0.0_real64  !! degC
      real(real64) :: hum   = 0.0_real64  !! kPa
      real(real64) :: wind  = 0.0_real64  !! m s-1
      real(real64) :: rain  = 0.0_real64  !! mm d-1
      real(real64) :: etref = 0.0_real64  !! mm d-1 (-99.9 = compute internally)
      real(real64) :: wet   = 0.0_real64  !! fraction (-99.9 = missing)
   end type meteo_daily_row_t

   type :: meteo_daily_table_t
      type(meteo_daily_row_t), allocatable :: rows(:)
      logical :: is_loaded = .false.
   contains
      procedure :: load        => meteo_daily_table_load
      procedure :: year_window => meteo_daily_table_year_window
   end type meteo_daily_table_t

   !> One sub-daily record. Column 1 is datetime (continuous fractional days
   !! since JD 1900), not whole-day; column 2 is the within-day record index.
   type :: meteo_detail_row_t
      real(real64) :: datetime = 0.0_real64
      integer      :: record   = 0
      real(real64) :: rad      = 0.0_real64  !! kJ m-2 (per-period total)
      real(real64) :: temp     = 0.0_real64  !! degC
      real(real64) :: hum      = 0.0_real64  !! kPa
      real(real64) :: wind     = 0.0_real64  !! m s-1
      real(real64) :: rain     = 0.0_real64  !! mm (per-period total)
   end type meteo_detail_row_t

   type :: meteo_detail_table_t
      type(meteo_detail_row_t), allocatable :: rows(:)
      logical :: is_loaded = .false.
   contains
      procedure :: load        => meteo_detail_table_load
      procedure :: year_window => meteo_detail_table_year_window
   end type meteo_detail_table_t

   !> One rain event (sub-daily rain accumulation row).
   type :: rain_event_row_t
      real(real64) :: datetime = 0.0_real64  !! days since JD 1900 (fractional)
      real(real64) :: amount   = 0.0_real64  !! mm (raw, no unit conversion)
   end type rain_event_row_t

   type :: rain_events_table_t
      type(rain_event_row_t), allocatable :: rows(:)
      logical :: is_loaded = .false.
   contains
      procedure :: load        => rain_events_table_load
      procedure :: year_window => rain_events_table_year_window
   end type rain_events_table_t

contains

   subroutine meteo_daily_table_load(self, path, errors)
      use csv_reader_mod, only: read_csv_table
      use error_mod,      only: ERR_VALIDATION_CROSS_FIELD, ERR_VALIDATION_OUT_OF_RANGE
      class(meteo_daily_table_t), intent(inout) :: self
      character(len=*),           intent(in)    :: path
      type(error_collection_t),   intent(inout) :: errors

      real(real64), allocatable :: tbl(:,:)
      character(len=5) :: hdr(9)
      integer :: n, r
      character(len=200) :: msg

      hdr = [character(len=5) :: 'date ', 'rad  ', 'tmin ', 'tmax ', &
             'hum  ', 'wind ', 'rain ', 'etref', 'wet  ']
      call read_csv_table(trim(path), hdr, tbl, errors)
      if (errors%has_fatals()) return

      n = size(tbl, 1)
      if (allocated(self%rows)) deallocate(self%rows)
      allocate(self%rows(n))
      do r = 1, n
         self%rows(r)%date  = tbl(r, 1)
         self%rows(r)%rad   = tbl(r, 2)
         self%rows(r)%tmin  = tbl(r, 3)
         self%rows(r)%tmax  = tbl(r, 4)
         self%rows(r)%hum   = tbl(r, 5)
         self%rows(r)%wind  = tbl(r, 6)
         self%rows(r)%rain  = tbl(r, 7)
         self%rows(r)%etref = tbl(r, 8)
         self%rows(r)%wet   = tbl(r, 9)
      end do

      ! Inline validate — load owns this, never call separately.
      do r = 1, n
         if (self%rows(r)%tmin > self%rows(r)%tmax) then
            write(msg, '("meteo daily row ", I0, ": tmin > tmax")') r
            call errors%append(ERR_VALIDATION_CROSS_FIELD, trim(msg), 'meteo_csv')
            deallocate(self%rows)
            return
         end if
         if (self%rows(r)%rain < 0.0_real64) then
            write(msg, '("meteo daily row ", I0, ": rain < 0")') r
            call errors%append(ERR_VALIDATION_OUT_OF_RANGE, trim(msg), 'meteo_csv')
            deallocate(self%rows)
            return
         end if
      end do

      self%is_loaded = .true.
   end subroutine meteo_daily_table_load

   subroutine meteo_daily_table_year_window(self, year, i1, i2)
      class(meteo_daily_table_t), intent(in)  :: self
      integer,                    intent(in)  :: year
      integer,                    intent(out) :: i1, i2

      integer  :: i
      real(real64) :: t_jan1, t_dec31

      i1 = 0; i2 = 0
      if (.not. self%is_loaded) return
      if (.not. allocated(self%rows)) return

      t_jan1  = real(days_since_1900(year,  1,  1), real64)
      t_dec31 = real(days_since_1900(year, 12, 31), real64)

      ! Assumes rows are sorted ascending by date (guaranteed by csv_reader's
      ! sequential read). i1 = first in-window index, i2 = last.
      do i = 1, size(self%rows)
         if (self%rows(i)%date >= t_jan1 - 0.5_real64 .and. &
             self%rows(i)%date <= t_dec31 + 0.5_real64) then
            if (i1 == 0) i1 = i
            i2 = i
         end if
      end do
   end subroutine meteo_daily_table_year_window

   subroutine meteo_detail_table_load(self, path, errors)
      use csv_reader_mod, only: read_csv_table
      class(meteo_detail_table_t), intent(inout) :: self
      character(len=*),            intent(in)    :: path
      type(error_collection_t),    intent(inout) :: errors

      real(real64), allocatable :: tbl(:,:)
      character(len=8) :: hdr(7)
      integer :: n, r

      hdr = [character(len=8) :: 'datetime', 'record  ', 'rad     ', &
             'temp    ', 'hum     ', 'wind    ', 'rain    ']
      call read_csv_table(trim(path), hdr, tbl, errors)
      if (errors%has_fatals()) return

      n = size(tbl, 1)
      if (allocated(self%rows)) deallocate(self%rows)
      allocate(self%rows(n))
      do r = 1, n
         self%rows(r)%datetime = tbl(r, 1)
         self%rows(r)%record   = nint(tbl(r, 2))
         self%rows(r)%rad      = tbl(r, 3)
         self%rows(r)%temp     = tbl(r, 4)
         self%rows(r)%hum      = tbl(r, 5)
         self%rows(r)%wind     = tbl(r, 6)
         self%rows(r)%rain     = tbl(r, 7)
      end do

      self%is_loaded = .true.
   end subroutine meteo_detail_table_load

   subroutine meteo_detail_table_year_window(self, year, i1, i2)
      class(meteo_detail_table_t), intent(in)  :: self
      integer,                     intent(in)  :: year
      integer,                     intent(out) :: i1, i2

      integer :: i
      real(real64) :: t_jan1, t_jan1_next

      i1 = 0; i2 = 0
      if (.not. self%is_loaded) return
      if (.not. allocated(self%rows)) return

      ! Half-open interval [t_jan1, t_jan1_next) — sub-daily timestamps are
      ! fractional days, so the daily ±0.5 slack would mis-bucket Dec 31
      ! noon-to-midnight into the next year (mirrors legacy MeteoCSVDetYear).
      t_jan1      = real(days_since_1900(year,     1, 1), real64)
      t_jan1_next = real(days_since_1900(year + 1, 1, 1), real64)

      ! Assumes rows are sorted ascending by datetime (guaranteed by csv_reader's
      ! sequential read). i1 = first in-window index, i2 = last.
      do i = 1, size(self%rows)
         if (self%rows(i)%datetime >= t_jan1 .and. &
             self%rows(i)%datetime <  t_jan1_next) then
            if (i1 == 0) i1 = i
            i2 = i
         end if
      end do
   end subroutine meteo_detail_table_year_window

   subroutine rain_events_table_load(self, path, errors)
      class(rain_events_table_t), intent(inout) :: self
      character(len=*),           intent(in)    :: path
      type(error_collection_t),   intent(inout) :: errors
      ! TASK 4
   end subroutine rain_events_table_load

   subroutine rain_events_table_year_window(self, year, i1, i2)
      class(rain_events_table_t), intent(in)  :: self
      integer,                    intent(in)  :: year
      integer,                    intent(out) :: i1, i2
      i1 = 0; i2 = 0  ! TASK 4
   end subroutine rain_events_table_year_window

   ! Private helper — same epoch and convention as csv_reader's date col.
   ! Fliegel/Van Flandern formula (JD), then subtract jd1900 = 2415020.
   pure function days_since_1900(year, month, day) result(d)
      integer, intent(in) :: year, month, day
      integer :: d
      integer :: a, y, m, jd
      integer, parameter :: jd1900 = 2415020
      a = (14 - month) / 12
      y = year + 4800 - a
      m = month + 12*a - 3
      jd = day + (153*m + 2)/5 + 365*y + y/4 - y/100 + y/400 - 32045
      d = jd - jd1900
   end function days_since_1900

end module meteo_csv_mod
