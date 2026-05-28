!> Typed CSV record tables for bottom-boundary condition files. (ADR 0044)
!!
!! Five table types, one per swbotb option:
!!   * gwl_table_t    — 2-column (date, gwl),    swbotb=1
!!   * qbot_table_t   — 2-column (date, qbot),   swbotb=2 (sw2=2) and swbotb=3 (sw4=1)
!!   * haquif_table_t — 2-column (date, haquif), swbotb=3 (sw3=2)
!!   * qhbot_table_t  — 2-column (htab, qtab),   swbotb=4 (swqhbot=2)
!!   * hbot_table_t   — 2-column (date, hbot),   swbotb=5
!!
!! Consumer pattern: these loaders are intermediate parsing objects —
!! seed_bottom_boundary in soilwater_state.f90 copies typed rows into
!! the interleaved (odd=key, even=value) state arrays gwltab/qbotab/
!! haqtab/hbotab that afgen() expects.
!!
!! All date-keyed tables require strictly ascending dates; qhbot_table_t
!! requires strictly ascending abs(htab) values (the afgen lookup key is
!! abs(gwl), matching abs(htab) stored in seed_bottom_boundary). The
!! load() method validates this in each case.
module boundary_csv_mod
   use iso_fortran_env, only: real64
   use error_mod,       only: error_collection_t
   implicit none
   private

   public :: gwl_table_t
   public :: gwl_row_t
   public :: qbot_table_t
   public :: qbot_row_t
   public :: haquif_table_t
   public :: haquif_row_t
   public :: qhbot_table_t
   public :: qhbot_row_t
   public :: hbot_table_t
   public :: hbot_row_t

   !> One groundwater-level record (swbotb=1).
   type :: gwl_row_t
      real(real64) :: date = 0.0_real64  !! days since JD 1900
      real(real64) :: gwl  = 0.0_real64  !! groundwater level (cm)
   end type gwl_row_t

   type :: gwl_table_t
      type(gwl_row_t), allocatable :: rows(:)
      logical :: is_loaded = .false.
   contains
      procedure :: load => gwl_table_load
   end type gwl_table_t

   !> One bottom-flux record (swbotb=2 sw2=2 and swbotb=3 sw4=1).
   type :: qbot_row_t
      real(real64) :: date = 0.0_real64  !! days since JD 1900
      real(real64) :: qbot = 0.0_real64  !! bottom flux (cm/d)
   end type qbot_row_t

   type :: qbot_table_t
      type(qbot_row_t), allocatable :: rows(:)
      logical :: is_loaded = .false.
   contains
      procedure :: load => qbot_table_load
   end type qbot_table_t

   !> One deep-aquifer head record (swbotb=3 sw3=2).
   type :: haquif_row_t
      real(real64) :: date   = 0.0_real64  !! days since JD 1900
      real(real64) :: haquif = 0.0_real64  !! deep-aquifer head (cm)
   end type haquif_row_t

   type :: haquif_table_t
      type(haquif_row_t), allocatable :: rows(:)
      logical :: is_loaded = .false.
   contains
      procedure :: load => haquif_table_load
   end type haquif_table_t

   !> One q(h) curve record (swbotb=4 swqhbot=2).
   !! Column 1 is htab (pressure head, cm) — NOT a date.
   !! The copy site in seed_bottom_boundary applies abs(htab) per legacy
   !! readswap.f90:1418-1419.
   type :: qhbot_row_t
      real(real64) :: htab = 0.0_real64  !! pressure head value (cm)
      real(real64) :: qtab = 0.0_real64  !! flux at that head (cm/d)
   end type qhbot_row_t

   type :: qhbot_table_t
      type(qhbot_row_t), allocatable :: rows(:)
      logical :: is_loaded = .false.
   contains
      procedure :: load => qhbot_table_load
   end type qhbot_table_t

   !> One bottom pressure-head record (swbotb=5).
   type :: hbot_row_t
      real(real64) :: date = 0.0_real64  !! days since JD 1900
      real(real64) :: hbot = 0.0_real64  !! bottom pressure head (cm)
   end type hbot_row_t

   type :: hbot_table_t
      type(hbot_row_t), allocatable :: rows(:)
      logical :: is_loaded = .false.
   contains
      procedure :: load => hbot_table_load
   end type hbot_table_t

contains

   subroutine gwl_table_load(self, path, errors)
      use csv_reader_mod, only: read_csv_table
      use error_mod,      only: ERR_VALIDATION_OUT_OF_RANGE
      class(gwl_table_t),       intent(inout) :: self
      character(len=*),         intent(in)    :: path
      type(error_collection_t), intent(inout) :: errors

      real(real64), allocatable :: tbl(:,:)
      character(len=14) :: hdr(2)
      integer :: n, r
      character(len=200) :: msg

      ! Reset is_loaded so a failed load on a reused instance doesn't
      ! leave a stale .true. flag from a previous successful load.
      self%is_loaded = .false.

      hdr = [character(len=14) :: 'date          ', 'gwl           ']
      call read_csv_table(trim(path), hdr, tbl, errors)
      if (errors%has_fatals()) return

      n = size(tbl, 1)
      if (allocated(self%rows)) deallocate(self%rows)
      allocate(self%rows(n))
      do r = 1, n
         self%rows(r)%date = tbl(r, 1)
         self%rows(r)%gwl  = tbl(r, 2)
      end do

      ! Inline validate — rows must be sorted strictly ascending by date.
      ! afgen() requires a sorted lookup table; unsorted input produces
      ! silently-wrong interpolation results at runtime.
      do r = 2, n
         if (self%rows(r)%date <= self%rows(r-1)%date) then
            write(msg, '("boundary gwl row ", I0, ": date not strictly ascending")') r
            call errors%append(ERR_VALIDATION_OUT_OF_RANGE, trim(msg), 'boundary_csv')
            deallocate(self%rows)
            return
         end if
      end do

      self%is_loaded = .true.
   end subroutine gwl_table_load


   subroutine qbot_table_load(self, path, errors)
      use csv_reader_mod, only: read_csv_table
      use error_mod,      only: ERR_VALIDATION_OUT_OF_RANGE
      class(qbot_table_t),      intent(inout) :: self
      character(len=*),         intent(in)    :: path
      type(error_collection_t), intent(inout) :: errors

      real(real64), allocatable :: tbl(:,:)
      character(len=14) :: hdr(2)
      integer :: n, r
      character(len=200) :: msg

      self%is_loaded = .false.

      hdr = [character(len=14) :: 'date          ', 'qbot          ']
      call read_csv_table(trim(path), hdr, tbl, errors)
      if (errors%has_fatals()) return

      n = size(tbl, 1)
      if (allocated(self%rows)) deallocate(self%rows)
      allocate(self%rows(n))
      do r = 1, n
         self%rows(r)%date = tbl(r, 1)
         self%rows(r)%qbot = tbl(r, 2)
      end do

      ! Inline validate — rows must be sorted strictly ascending by date.
      do r = 2, n
         if (self%rows(r)%date <= self%rows(r-1)%date) then
            write(msg, '("boundary qbot row ", I0, ": date not strictly ascending")') r
            call errors%append(ERR_VALIDATION_OUT_OF_RANGE, trim(msg), 'boundary_csv')
            deallocate(self%rows)
            return
         end if
      end do

      self%is_loaded = .true.
   end subroutine qbot_table_load


   subroutine haquif_table_load(self, path, errors)
      use csv_reader_mod, only: read_csv_table
      use error_mod,      only: ERR_VALIDATION_OUT_OF_RANGE
      class(haquif_table_t),    intent(inout) :: self
      character(len=*),         intent(in)    :: path
      type(error_collection_t), intent(inout) :: errors

      real(real64), allocatable :: tbl(:,:)
      character(len=14) :: hdr(2)
      integer :: n, r
      character(len=200) :: msg

      self%is_loaded = .false.

      hdr = [character(len=14) :: 'date          ', 'haquif        ']
      call read_csv_table(trim(path), hdr, tbl, errors)
      if (errors%has_fatals()) return

      n = size(tbl, 1)
      if (allocated(self%rows)) deallocate(self%rows)
      allocate(self%rows(n))
      do r = 1, n
         self%rows(r)%date   = tbl(r, 1)
         self%rows(r)%haquif = tbl(r, 2)
      end do

      ! Inline validate — rows must be sorted strictly ascending by date.
      do r = 2, n
         if (self%rows(r)%date <= self%rows(r-1)%date) then
            write(msg, '("boundary haquif row ", I0, ": date not strictly ascending")') r
            call errors%append(ERR_VALIDATION_OUT_OF_RANGE, trim(msg), 'boundary_csv')
            deallocate(self%rows)
            return
         end if
      end do

      self%is_loaded = .true.
   end subroutine haquif_table_load


   subroutine qhbot_table_load(self, path, errors)
      use csv_reader_mod, only: read_csv_table
      use error_mod,      only: ERR_VALIDATION_OUT_OF_RANGE
      class(qhbot_table_t),     intent(inout) :: self
      character(len=*),         intent(in)    :: path
      type(error_collection_t), intent(inout) :: errors

      real(real64), allocatable :: tbl(:,:)
      character(len=14) :: hdr(2)
      integer :: n, r
      character(len=200) :: msg

      self%is_loaded = .false.

      hdr = [character(len=14) :: 'htab          ', 'qtab          ']
      call read_csv_table(trim(path), hdr, tbl, errors)
      if (errors%has_fatals()) return

      n = size(tbl, 1)
      if (allocated(self%rows)) deallocate(self%rows)
      allocate(self%rows(n))
      do r = 1, n
         self%rows(r)%htab = tbl(r, 1)
         self%rows(r)%qtab = tbl(r, 2)
      end do

      ! qhbot stored keys are abs(htab) for afgen lookup against abs(gwl).
      ! Physical q(h) tables typically have negative htab values; the abs
      ! transformation flips their natural ordering. Validate the ordering
      ! AFTER the abs transformation, not before.
      do r = 2, n
         if (abs(self%rows(r)%htab) <= abs(self%rows(r-1)%htab)) then
            write(msg, '("boundary qhbot row ", I0, ": abs(htab) not strictly ascending")') r
            call errors%append(ERR_VALIDATION_OUT_OF_RANGE, trim(msg), 'boundary_csv')
            deallocate(self%rows)
            return
         end if
      end do

      self%is_loaded = .true.
   end subroutine qhbot_table_load


   subroutine hbot_table_load(self, path, errors)
      use csv_reader_mod, only: read_csv_table
      use error_mod,      only: ERR_VALIDATION_OUT_OF_RANGE
      class(hbot_table_t),      intent(inout) :: self
      character(len=*),         intent(in)    :: path
      type(error_collection_t), intent(inout) :: errors

      real(real64), allocatable :: tbl(:,:)
      character(len=14) :: hdr(2)
      integer :: n, r
      character(len=200) :: msg

      self%is_loaded = .false.

      hdr = [character(len=14) :: 'date          ', 'hbot          ']
      call read_csv_table(trim(path), hdr, tbl, errors)
      if (errors%has_fatals()) return

      n = size(tbl, 1)
      if (allocated(self%rows)) deallocate(self%rows)
      allocate(self%rows(n))
      do r = 1, n
         self%rows(r)%date = tbl(r, 1)
         self%rows(r)%hbot = tbl(r, 2)
      end do

      ! Inline validate — rows must be sorted strictly ascending by date.
      do r = 2, n
         if (self%rows(r)%date <= self%rows(r-1)%date) then
            write(msg, '("boundary hbot row ", I0, ": date not strictly ascending")') r
            call errors%append(ERR_VALIDATION_OUT_OF_RANGE, trim(msg), 'boundary_csv')
            deallocate(self%rows)
            return
         end if
      end do

      self%is_loaded = .true.
   end subroutine hbot_table_load

end module boundary_csv_mod
