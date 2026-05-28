!> Typed CSV record table for drainage open-water-level tables. (ADR 0044)
!!
!! One table type:
!!   * owl_events_table_t — 2-column per-level (date, level)
!!
!! Schema: `date` (ISO date → days-since-1900) + `level` (cm).
!! Consumer pattern: this loader is intermediate — drainage_state_init
!! converts the typed rows to the legacy `owltab(:, 2k-1:2k)` interleaved
!! (time, level) pairs that afgen() expects. See src/drainage/drainage.f90.
!!
!! Rows MUST be strictly ascending by date (afgen requires a sorted lookup
!! table). The load() method validates this and sets errors if violated.
module drainage_csv_mod
   use iso_fortran_env, only: real64
   use error_mod,       only: error_collection_t
   implicit none
   private

   public :: owl_events_table_t
   public :: owl_event_row_t

   type :: owl_event_row_t
      real(real64) :: date  = 0.0_real64  !! days since JD 1900
      real(real64) :: level = 0.0_real64  !! cm (surface water level)
   end type owl_event_row_t

   type :: owl_events_table_t
      type(owl_event_row_t), allocatable :: rows(:)
      logical :: is_loaded = .false.
   contains
      procedure :: load => owl_events_table_load
   end type owl_events_table_t

contains

   subroutine owl_events_table_load(self, path, errors)
      use csv_reader_mod, only: read_csv_table
      use error_mod,      only: ERR_VALIDATION_OUT_OF_RANGE
      class(owl_events_table_t), intent(inout) :: self
      character(len=*),          intent(in)    :: path
      type(error_collection_t),  intent(inout) :: errors

      real(real64), allocatable :: tbl(:,:)
      character(len=8) :: hdr(2)
      integer :: n, r
      character(len=200) :: msg

      ! Reset is_loaded so a failed load on a reused instance doesn't
      ! leave a stale .true. flag from a previous successful load.
      self%is_loaded = .false.

      hdr = [character(len=8) :: 'date    ', 'level   ']
      call read_csv_table(trim(path), hdr, tbl, errors)
      if (errors%has_fatals()) return

      n = size(tbl, 1)
      if (allocated(self%rows)) deallocate(self%rows)
      allocate(self%rows(n))
      do r = 1, n
         self%rows(r)%date  = tbl(r, 1)
         self%rows(r)%level = tbl(r, 2)
      end do

      ! Inline validate — rows must be sorted strictly ascending by date.
      ! afgen() requires a sorted lookup table; violating this silently
      ! produces wrong interpolation results at runtime.
      do r = 2, n
         if (self%rows(r)%date <= self%rows(r-1)%date) then
            write(msg, '("drainage owltab row ", I0, ": date not strictly ascending")') r
            call errors%append(ERR_VALIDATION_OUT_OF_RANGE, trim(msg), 'drainage_csv')
            deallocate(self%rows)
            return
         end if
      end do

      self%is_loaded = .true.
   end subroutine owl_events_table_load

end module drainage_csv_mod
