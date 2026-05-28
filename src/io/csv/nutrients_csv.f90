!> Typed CSV record table for nutrients amendment event files. (ADR 0044)
!!
!! One table type:
!!   * amendment_events_table_t — 4-column per-event (date, material, amount_kgha, volat_fraction)
!!
!! Schema:
!!   date           — ISO date → days-since-1900
!!   material       — integer (1..20; stored as real64, consumer applies nint())
!!   amount_kgha    — real64 (kg/ha; range [0, 500000]; consumer converts to kg/m² via ×1e-4)
!!   volat_fraction — real64 (range [0, 1])
!!
!! Consumer pattern: this loader is intermediate — load_nutrients_events in
!! nutrients_state.f90 copies typed rows into legacy wofost_soil_declarations
!! globals (MatNum, Amend, VolaFrac, TimeAmend, iamend, NuAmend, namend, isme)
!! with unit conversions applied at the copy site. The post-load sorting and
!! date-grouping logic remains in the consumer. Globals retirement is a
!! separate arc.
!!
!! No year_window constraint — consumers match events via TimeAmend/iamend
!! directly against simulation dates.
module nutrients_csv_mod
   use iso_fortran_env, only: real64
   use error_mod,       only: error_collection_t
   implicit none
   private

   public :: amendment_events_table_t
   public :: amendment_event_row_t

   !> One amendment event record.
   !! material is stored as integer (range validated to 1..20).
   !! amount_kgha is in kg/ha; consumer converts to kg/m² via ×1e-4.
   type :: amendment_event_row_t
      real(real64) :: date           = 0.0_real64  !! days since JD 1900
      integer      :: material       = 1            !! amendment material index (1..20)
      real(real64) :: amount_kgha    = 0.0_real64  !! kg/ha (consumer converts to kg/m² via ×1e-4)
      real(real64) :: volat_fraction = 0.0_real64  !! volatilisation fraction [0, 1]
   end type amendment_event_row_t

   type :: amendment_events_table_t
      type(amendment_event_row_t), allocatable :: rows(:)
      logical :: is_loaded = .false.
   contains
      procedure :: load => amendment_events_table_load
   end type amendment_events_table_t

contains

   subroutine amendment_events_table_load(self, path, errors)
      use csv_reader_mod, only: read_csv_table
      use error_mod,      only: ERR_VALIDATION_OUT_OF_RANGE
      class(amendment_events_table_t), intent(inout) :: self
      character(len=*),                intent(in)    :: path
      type(error_collection_t),        intent(inout) :: errors

      real(real64), allocatable :: tbl(:,:)
      character(len=14) :: hdr(4)
      integer :: n, r
      character(len=200) :: msg

      ! Reset is_loaded so a failed load on a reused instance doesn't
      ! leave a stale .true. flag from a previous successful load.
      self%is_loaded = .false.

      hdr = [character(len=14) :: 'date          ', 'material      ', &
                                   'amount_kgha   ', 'volat_fraction']
      call read_csv_table(trim(path), hdr, tbl, errors)
      if (errors%has_fatals()) return

      n = size(tbl, 1)
      if (allocated(self%rows)) deallocate(self%rows)
      allocate(self%rows(n))
      do r = 1, n
         self%rows(r)%date           = tbl(r, 1)
         self%rows(r)%material       = nint(tbl(r, 2))
         self%rows(r)%amount_kgha    = tbl(r, 3)
         self%rows(r)%volat_fraction = tbl(r, 4)
      end do

      ! Inline validate — material must be in [1, 20].
      do r = 1, n
         if (self%rows(r)%material < 1 .or. self%rows(r)%material > 20) then
            write(msg, '("nutrients amendment_events row ", I0, &
               &": material out of range [1, 20], got ", I0)') &
               r, self%rows(r)%material
            call errors%append(ERR_VALIDATION_OUT_OF_RANGE, trim(msg), 'nutrients_csv')
            deallocate(self%rows)
            return
         end if
      end do

      ! Inline validate — amount_kgha must be in [0, 500000].
      do r = 1, n
         if (self%rows(r)%amount_kgha < 0.0_real64 .or. &
             self%rows(r)%amount_kgha > 500000.0_real64) then
            write(msg, '("nutrients amendment_events row ", I0, &
               &": amount_kgha out of range [0, 500000], got ", ES12.4)') &
               r, self%rows(r)%amount_kgha
            call errors%append(ERR_VALIDATION_OUT_OF_RANGE, trim(msg), 'nutrients_csv')
            deallocate(self%rows)
            return
         end if
      end do

      ! Inline validate — volat_fraction must be in [0, 1].
      do r = 1, n
         if (self%rows(r)%volat_fraction < 0.0_real64 .or. &
             self%rows(r)%volat_fraction > 1.0_real64) then
            write(msg, '("nutrients amendment_events row ", I0, &
               &": volat_fraction out of range [0, 1], got ", F8.4)') &
               r, self%rows(r)%volat_fraction
            call errors%append(ERR_VALIDATION_OUT_OF_RANGE, trim(msg), 'nutrients_csv')
            deallocate(self%rows)
            return
         end if
      end do

      self%is_loaded = .true.
   end subroutine amendment_events_table_load

end module nutrients_csv_mod
