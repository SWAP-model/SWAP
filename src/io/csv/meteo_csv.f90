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

   ! Stubs — real implementations land in Tasks 2-4.
   subroutine meteo_daily_table_load(self, path, errors)
      class(meteo_daily_table_t), intent(inout) :: self
      character(len=*),           intent(in)    :: path
      type(error_collection_t),   intent(inout) :: errors
      ! TASK 2
   end subroutine meteo_daily_table_load

   subroutine meteo_daily_table_year_window(self, year, i1, i2)
      class(meteo_daily_table_t), intent(in)  :: self
      integer,                    intent(in)  :: year
      integer,                    intent(out) :: i1, i2
      i1 = 0; i2 = 0  ! TASK 2
   end subroutine meteo_daily_table_year_window

   subroutine meteo_detail_table_load(self, path, errors)
      class(meteo_detail_table_t), intent(inout) :: self
      character(len=*),            intent(in)    :: path
      type(error_collection_t),    intent(inout) :: errors
      ! TASK 3
   end subroutine meteo_detail_table_load

   subroutine meteo_detail_table_year_window(self, year, i1, i2)
      class(meteo_detail_table_t), intent(in)  :: self
      integer,                     intent(in)  :: year
      integer,                     intent(out) :: i1, i2
      i1 = 0; i2 = 0  ! TASK 3
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

end module meteo_csv_mod
