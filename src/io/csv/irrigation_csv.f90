!> Typed CSV record tables for irrigation event files. (ADR 0044)
!!
!! Two table types:
!!   * fixed_events_table_t — 4-column fixed irrigation events (date, depth_mm, conc, type)
!!   * ssdi_events_table_t  — 3-column SSDI fixed-date events (date, rate_f_mm_h, amount_f_mm)
!!
!! Schema for fixed_events:
!!   date     — ISO date → days-since-1900
!!   depth    — real64 (mm; consumer divides by 10.0 for cm)
!!   conc     — real64 (concentration M/L³)
!!   irtype   — real64 stored (0.0=sprinkler, 1.0=surface); consumer uses nint()
!!
!! Schema for ssdi_events:
!!   date     — ISO date → days-since-1900
!!   rate_f   — real64 (mm/h; consumer multiplies ×0.1×24 for cm/d)
!!   amount_f — real64 (mm; consumer multiplies ×0.1 and divides by ncomp)
!!
!! Consumer pattern: these loaders are intermediate parsing objects —
!! crop_irrigation_state_init's seed_fixed_irrigation and apply_ssdi_mode0
!! copy typed rows into fixed-size state arrays with unit conversions
!! applied at the copy site. This matches Family 1 (drainage owltab).
!!
!! Rows in ssdi_events MUST be strictly ascending by date (the SSDI loop in
!! irrigation.f90 walks forward through dates monotonically). The load()
!! method validates this and sets errors if violated.
module irrigation_csv_mod
   use iso_fortran_env, only: real64
   use error_mod,       only: error_collection_t
   implicit none
   private

   public :: fixed_events_table_t
   public :: fixed_event_row_t
   public :: ssdi_events_table_t
   public :: ssdi_event_row_t

   !> One fixed irrigation event record.
   !! irtype is stored as real64 (0.0 or 1.0); consumer applies nint().
   type :: fixed_event_row_t
      real(real64) :: date   = 0.0_real64  !! days since JD 1900
      real(real64) :: depth  = 0.0_real64  !! mm (consumer converts to cm via /10.0)
      real(real64) :: conc   = 0.0_real64  !! concentration M/L³
      real(real64) :: irtype = 0.0_real64  !! 0.0=sprinkler, 1.0=surface (consumer uses nint)
   end type fixed_event_row_t

   type :: fixed_events_table_t
      type(fixed_event_row_t), allocatable :: rows(:)
      logical :: is_loaded = .false.
   contains
      procedure :: load => fixed_events_table_load
   end type fixed_events_table_t

   !> One SSDI fixed-date event record.
   type :: ssdi_event_row_t
      real(real64) :: date     = 0.0_real64  !! days since JD 1900
      real(real64) :: rate_f   = 0.0_real64  !! mm/h (consumer converts to cm/d via ×0.1×24)
      real(real64) :: amount_f = 0.0_real64  !! mm (consumer converts to cm/compartment via ×0.1/ncomp)
   end type ssdi_event_row_t

   type :: ssdi_events_table_t
      type(ssdi_event_row_t), allocatable :: rows(:)
      logical :: is_loaded = .false.
   contains
      procedure :: load => ssdi_events_table_load
   end type ssdi_events_table_t

contains

   subroutine fixed_events_table_load(self, path, errors)
      use csv_reader_mod, only: read_csv_table
      use error_mod,      only: ERR_VALIDATION_ENUM
      class(fixed_events_table_t), intent(inout) :: self
      character(len=*),            intent(in)    :: path
      type(error_collection_t),    intent(inout) :: errors

      real(real64), allocatable :: tbl(:,:)
      character(len=8) :: hdr(4)
      integer :: n, r
      character(len=200) :: msg

      ! Reset is_loaded so a failed load on a reused instance doesn't
      ! leave a stale .true. flag from a previous successful load.
      self%is_loaded = .false.

      hdr = [character(len=8) :: 'date    ', 'depth   ', 'conc    ', 'type    ']
      call read_csv_table(trim(path), hdr, tbl, errors)
      if (errors%has_fatals()) return

      n = size(tbl, 1)
      if (allocated(self%rows)) deallocate(self%rows)
      allocate(self%rows(n))
      do r = 1, n
         self%rows(r)%date   = tbl(r, 1)
         self%rows(r)%depth  = tbl(r, 2)
         self%rows(r)%conc   = tbl(r, 3)
         self%rows(r)%irtype = tbl(r, 4)
      end do

      ! Inline validate — irtype must be 0 or 1.
      do r = 1, n
         if (nint(self%rows(r)%irtype) /= 0 .and. nint(self%rows(r)%irtype) /= 1) then
            write(msg, '("irrigation fixed_events row ", I0, ": irtype must be 0 or 1, got ", F6.1)') &
               r, self%rows(r)%irtype
            call errors%append(ERR_VALIDATION_ENUM, trim(msg), 'irrigation_csv')
            deallocate(self%rows)
            return
         end if
      end do

      self%is_loaded = .true.
   end subroutine fixed_events_table_load


   subroutine ssdi_events_table_load(self, path, errors)
      use csv_reader_mod, only: read_csv_table
      use error_mod,      only: ERR_VALIDATION_OUT_OF_RANGE
      class(ssdi_events_table_t), intent(inout) :: self
      character(len=*),           intent(in)    :: path
      type(error_collection_t),   intent(inout) :: errors

      real(real64), allocatable :: tbl(:,:)
      character(len=8) :: hdr(3)
      integer :: n, r
      character(len=200) :: msg

      ! Reset is_loaded so a failed load on a reused instance doesn't
      ! leave a stale .true. flag from a previous successful load.
      self%is_loaded = .false.

      hdr = [character(len=8) :: 'date    ', 'rate_f  ', 'amount_f']
      call read_csv_table(trim(path), hdr, tbl, errors)
      if (errors%has_fatals()) return

      n = size(tbl, 1)
      if (allocated(self%rows)) deallocate(self%rows)
      allocate(self%rows(n))
      do r = 1, n
         self%rows(r)%date     = tbl(r, 1)
         self%rows(r)%rate_f   = tbl(r, 2)
         self%rows(r)%amount_f = tbl(r, 3)
      end do

      ! Inline validate — rows must be sorted strictly ascending by date.
      ! The SSDI irrigation loop in irrigation.f90 walks forward through
      ! dates monotonically; unsorted input produces silently-wrong results.
      do r = 2, n
         if (self%rows(r)%date <= self%rows(r-1)%date) then
            write(msg, '("irrigation ssdi_events row ", I0, ": date not strictly ascending")') r
            call errors%append(ERR_VALIDATION_OUT_OF_RANGE, trim(msg), 'irrigation_csv')
            deallocate(self%rows)
            return
         end if
      end do

      self%is_loaded = .true.
   end subroutine ssdi_events_table_load

end module irrigation_csv_mod
