!> Phase 4c-b Task 7: test-side wrapper for the legacy free subroutines
!! `readwofost`, `readcropfixed`, `readgrass` in `src/io/readswap.f90`.
!!
!! These readers live outside any module and are called from
!! `cropgrowth.f90` per rotation entry during simulation. `readswap()`
!! itself never invokes them, so parity tests need this bridge to
!! populate `variables%` globals (dtsmtb, slatb, amaxtb, etc.) directly
!! from a .crp file without running the full crop growth task.
!!
!! Lives under tests/; never shipped with the production binary.
module legacy_crop_helper_mod
   use iso_fortran_env, only: real64
   implicit none
   private

   public :: read_legacy_wofost
   public :: read_legacy_cropfixed
   public :: read_legacy_grass
   public :: flatten_table

   ! Explicit interface for the free subroutine `readwofost` exposed at
   ! file scope by `src/io/readswap.f90`. gfortran does not require this
   ! for an implicit-none caller as long as the linker resolves the
   ! symbol, but declaring it here gives us argument-checking inside
   ! `read_legacy_wofost`.
   interface
      subroutine readwofost(icrop, crpfil, swhydrlift, swsoybean, mg, dvsi, &
                            dvrmax1, dvrmax2, flrfphotoveg, tmaxdvr, tmindvr, &
                            toptdvr, popt, pcrt, flphenodayl, &
                            FraDeceasedLvToSoil)
         integer,            intent(in)    :: icrop
         character(len=*),   intent(in)    :: crpfil
         integer,            intent(in)    :: swhydrlift
         integer,            intent(inout) :: swsoybean
         real(8),            intent(inout) :: mg, dvsi, dvrmax1, dvrmax2
         logical,            intent(inout) :: flrfphotoveg
         real(8),            intent(inout) :: tmaxdvr, tmindvr, toptdvr
         real(8),            intent(inout) :: popt, pcrt
         logical,            intent(inout) :: flphenodayl
         real(8),            intent(inout) :: FraDeceasedLvToSoil
      end subroutine readwofost
   end interface

contains

   !> Phase 4c-b Task 7: call legacy `readwofost` for rotation entry
   !! `icrop` with crop file `crpfil` (basename, no .crp extension).
   !! Side-effect populates the long list of `variables%` globals
   !! enumerated in the `use variables, only:` clause of `readwofost`
   !! itself. Local variables here exist only to absorb the read's
   !! output arguments — they are discarded on return.
   subroutine read_legacy_wofost(icrop, crpfil)
      integer,          intent(in) :: icrop
      character(len=*), intent(in) :: crpfil

      ! Inputs (conservative defaults matching hupselbrook potatod usage).
      integer :: swhydrlift
      integer :: swsoybean

      ! Outputs absorbed and discarded.
      real(8) :: mg, dvsi, dvrmax1, dvrmax2
      real(8) :: tmaxdvr, tmindvr, toptdvr, popt, pcrt
      real(8) :: FraDeceasedLvToSoil
      logical :: flrfphotoveg, flphenodayl

      swhydrlift = 0
      swsoybean  = 0
      ! Initialise the reals/logicals so any code path inside readwofost
      ! that does not assign them does not read uninitialised memory.
      mg = 0.0d0; dvsi = 0.0d0; dvrmax1 = 0.0d0; dvrmax2 = 0.0d0
      tmaxdvr = 0.0d0; tmindvr = 0.0d0; toptdvr = 0.0d0
      popt = 0.0d0; pcrt = 0.0d0
      FraDeceasedLvToSoil = 0.0d0
      flrfphotoveg = .false.
      flphenodayl  = .false.

      call readwofost(icrop, crpfil, swhydrlift, swsoybean, mg, dvsi, &
                      dvrmax1, dvrmax2, flrfphotoveg, tmaxdvr, tmindvr, &
                      toptdvr, popt, pcrt, flphenodayl, FraDeceasedLvToSoil)
   end subroutine read_legacy_wofost

   !> Phase 4c-b Task 7 stub. Body lands in Task 13 once the cropfixed
   !! parity test is wired up.
   subroutine read_legacy_cropfixed(icrop, crpfil)
      integer,          intent(in) :: icrop
      character(len=*), intent(in) :: crpfil
      ! TODO: Phase 4c-b Task 13 — call readcropfixed(icrop, crpfil, lcc, swhydrlift).
      ! Suppress unused-argument warnings until Task 13 fills this in.
      if (.false.) print *, icrop, crpfil
   end subroutine read_legacy_cropfixed

   !> Phase 4c-b Task 7 stub. Body lands in Task 13.
   subroutine read_legacy_grass(icrop, crpfil)
      integer,          intent(in) :: icrop
      character(len=*), intent(in) :: crpfil
      ! TODO: Phase 4c-b Task 13 — call readgrass(icrop, crpfil, ...).
      if (.false.) print *, icrop, crpfil
   end subroutine read_legacy_grass

   !> Flatten a 2-D (nrows, 2) table into a length-`max_size` flat array
   !! interleaved as [x1, y1, x2, y2, ...] padded with zeros to match
   !! legacy AFGEN storage (e.g. `real(8) :: dtsmtb(30)`).
   !!
   !! If the input has zero rows, the result is all zeros.
   pure function flatten_table(table_2d, max_size) result(flat)
      real(real64), intent(in) :: table_2d(:,:)
      integer,      intent(in) :: max_size
      real(real64)             :: flat(max_size)
      integer :: nrows, i, j

      flat = 0.0_real64

      nrows = size(table_2d, 1)
      if (nrows <= 0) return
      if (size(table_2d, 2) < 2) return

      do i = 1, nrows
         j = 2 * i - 1
         if (j + 1 > max_size) exit
         flat(j)     = table_2d(i, 1)
         flat(j + 1) = table_2d(i, 2)
      end do
   end function flatten_table

end module legacy_crop_helper_mod
