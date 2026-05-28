!> Typed CSV record tables for soil initial-condition profiles. (ADR 0044)
!!
!! Two table types:
!!   * h_profile_table_t   — 2-column initial pressure-head profile (z, h)
!!   * cml_profile_table_t — 2-column initial solute concentration profile (z, cml)
!!
!! Both are only meaningful when swinco=3 (warm-restart from CSV).
!! cml_profile_table_t is additionally gated on swsolu=1.
!!
!! h_profile_table_t is a FULL state migration: the loaded table is stored
!! in state%soilwater%h_init and consumed by soilhydraulics.f90 directly.
!! It replaces the legacy mutation of config%soil%initial%z_init (W3 fix)
!! and eliminates the duplicate CSV read previously in swap_mod.f90 (W4 fix).
!!
!! cml_profile_table_t is a LOADER-ONLY intermediate: seed_state_from_config.f90
!! uses it to populate the existing state%solute%cml_init and zc_init arrays.
!!
!! No ascending-order validation is applied: the legacy code did not validate
!! depth ordering for these profiles, and soilhydraulics.f90's afgen consumer
!! applies abs() to depth values before interpolation (so ordering is delegated
!! to the caller's mesh). Consistent with the no-validation legacy contract.
module soil_init_csv_mod
   use iso_fortran_env, only: real64
   use error_mod,       only: error_collection_t
   implicit none
   private

   public :: h_profile_table_t
   public :: h_profile_row_t
   public :: cml_profile_table_t
   public :: cml_profile_row_t

   !> One initial pressure-head record.
   type :: h_profile_row_t
      real(real64) :: z = 0.0_real64  !! depth (cm); typically negative
      real(real64) :: h = 0.0_real64  !! pressure head (cm)
   end type h_profile_row_t

   !> Typed table for initial pressure-head profile (swinco=3, h_file CSV).
   !! [W3 fix 2026-05-28] State-side migration: stored in state%soilwater%h_init;
   !! soilwater_state_init no longer mutates config%soil%initial%z_init.
   type :: h_profile_table_t
      type(h_profile_row_t), allocatable :: rows(:)
      logical :: is_loaded = .false.
   contains
      procedure :: load => h_profile_table_load
   end type h_profile_table_t

   !> One initial solute concentration record.
   type :: cml_profile_row_t
      real(real64) :: z   = 0.0_real64  !! depth (cm)
      real(real64) :: cml = 0.0_real64  !! mobile concentration (M/L^3)
   end type cml_profile_row_t

   !> Typed table for initial solute concentration profile (swinco=3 + swsolu=1).
   !! Loader-only: seed_state_from_config copies rows into state%solute%cml_init/zc_init.
   type :: cml_profile_table_t
      type(cml_profile_row_t), allocatable :: rows(:)
      logical :: is_loaded = .false.
   contains
      procedure :: load => cml_profile_table_load
   end type cml_profile_table_t

contains

   subroutine h_profile_table_load(self, path, errors)
      use csv_reader_mod, only: read_csv_table
      class(h_profile_table_t), intent(inout) :: self
      character(len=*),         intent(in)    :: path
      type(error_collection_t), intent(inout) :: errors

      real(real64), allocatable :: tbl(:,:)
      character(len=14) :: hdr(2)
      integer :: n, r

      ! Reset so a failed load on a reused instance doesn't leave stale .true.
      self%is_loaded = .false.

      hdr = [character(len=14) :: 'z             ', 'h             ']
      call read_csv_table(trim(path), hdr, tbl, errors)
      if (errors%has_fatals()) return

      n = size(tbl, 1)
      if (allocated(self%rows)) deallocate(self%rows)
      allocate(self%rows(n))
      do r = 1, n
         self%rows(r)%z = tbl(r, 1)
         self%rows(r)%h = tbl(r, 2)
      end do

      self%is_loaded = .true.
   end subroutine h_profile_table_load

   subroutine cml_profile_table_load(self, path, errors)
      use csv_reader_mod, only: read_csv_table
      class(cml_profile_table_t), intent(inout) :: self
      character(len=*),           intent(in)    :: path
      type(error_collection_t),   intent(inout) :: errors

      real(real64), allocatable :: tbl(:,:)
      character(len=14) :: hdr(2)
      integer :: n, r

      ! Reset so a failed load on a reused instance doesn't leave stale .true.
      self%is_loaded = .false.

      hdr = [character(len=14) :: 'z             ', 'cml           ']
      call read_csv_table(trim(path), hdr, tbl, errors)
      if (errors%has_fatals()) return

      n = size(tbl, 1)
      if (allocated(self%rows)) deallocate(self%rows)
      allocate(self%rows(n))
      do r = 1, n
         self%rows(r)%z   = tbl(r, 1)
         self%rows(r)%cml = tbl(r, 2)
      end do

      self%is_loaded = .true.
   end subroutine cml_profile_table_load

end module soil_init_csv_mod
