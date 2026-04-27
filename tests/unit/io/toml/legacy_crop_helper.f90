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
   public :: read_legacy_grass_with_mowing
   public :: flatten_table
   public :: reset_legacy_crop_globals

   ! Explicit interfaces for the free subroutines `readwofost`,
   ! `readcropfixed`, and `readgrass` exposed at file scope by
   ! `src/io/readswap.f90`. gfortran does not require these for an
   ! implicit-none caller as long as the linker resolves the symbol,
   ! but declaring them here gives us argument-checking inside the
   ! wrappers below.
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

      subroutine readcropfixed(icrop, crpfil, lcc, swhydrlift)
         integer,          intent(in)    :: icrop
         character(len=*), intent(in)    :: crpfil
         integer,          intent(inout) :: lcc
         integer,          intent(in)    :: swhydrlift
      end subroutine readcropfixed

      subroutine readgrass(icrop, crpfil, swharvest, dmharvest, daylastharvest, &
                           dmlastharvest, swdmmow, maxdaymow, swlossmow,        &
                           swlossgrz, swdmgrz, maxdaygrz, dmgrazing, lsdb,      &
                           tagprest, swhydrlift)
         integer,          intent(in)    :: icrop
         character(len=*), intent(in)    :: crpfil
         integer,          intent(inout) :: swharvest
         real(8),          intent(inout) :: dmharvest
         integer,          intent(inout) :: daylastharvest
         real(8),          intent(inout) :: dmlastharvest
         integer,          intent(inout) :: swdmmow
         integer,          intent(inout) :: maxdaymow
         integer,          intent(inout) :: swlossmow
         integer,          intent(inout) :: swlossgrz
         integer,          intent(inout) :: swdmgrz
         integer,          intent(inout) :: maxdaygrz
         real(8),          intent(inout) :: dmgrazing
         real(8),          intent(inout) :: lsdb(100)
         real(8),          intent(inout) :: tagprest
         integer,          intent(in)    :: swhydrlift
      end subroutine readgrass
   end interface

contains

   !> Phase 4c-b Task 14: zero out the legacy-variables-module crop tables
   !! and the scalar root growth fields so that a previous test's read does
   !! not pollute the next read's parity comparison.
   !!
   !! The legacy AFGEN tables are flat arrays sized for the maximum number
   !! of rows. A reader only writes the first `2*ifnd` slots, leaving the
   !! tail untouched. When two crop reads in the same process have different
   !! `ifnd` counts, the tail of the second read is whatever the first read
   !! left there. flatten_table on the (smaller) TOML side pads with zeros,
   !! so the parity assertion fails on the polluted tail.
   !!
   !! Similarly, scalar root-growth fields (rdi/rri/rdc) are gated on
   !! `swrd` in the legacy reader: if the current crop's swrd path does
   !! not assign them, they keep the previous read's value.
   !!
   !! This routine resets every table and switch-gated scalar that the
   !! parity tests actually compare. Test-only — never used in production.
   subroutine reset_legacy_crop_globals()
      use variables, only: &
         dtsmtb, slatb, amaxtb, tmpftb, tmnftb, frtb, fltb, fstb, fotb, &
         rdrrtb, rdrstb, rfsetb, chtb, cftb, gctb, cfeictb, rdtb, rdctb, &
         rlwtb, rdi, rri, rdc, wrtmax, cofab

      dtsmtb = 0.0d0
      slatb  = 0.0d0
      amaxtb = 0.0d0
      tmpftb = 0.0d0
      tmnftb = 0.0d0
      frtb   = 0.0d0
      fltb   = 0.0d0
      fstb   = 0.0d0
      fotb   = 0.0d0
      rdrrtb = 0.0d0
      rdrstb = 0.0d0
      rfsetb = 0.0d0
      chtb   = 0.0d0
      cftb   = 0.0d0
      gctb   = 0.0d0
      cfeictb = 0.0d0
      rdtb   = 0.0d0
      rdctb  = 0.0d0
      rlwtb  = 0.0d0
      rdi    = 0.0d0
      rri    = 0.0d0
      rdc    = 0.0d0
      wrtmax = 0.0d0
      cofab  = 0.0d0
   end subroutine reset_legacy_crop_globals

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

      call reset_legacy_crop_globals()

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

   !> Phase 4c-b Task 13: call legacy `readcropfixed` for rotation entry
   !! `icrop` with crop file `crpfil` (basename, no .crp extension).
   !! Side-effect populates the `variables%` globals enumerated in the
   !! `use variables, only:` clause of `readcropfixed` itself (kdif, kdir,
   !! hlim1..hlim4, gctb, cftb, rdctb, etc.). Local `lcc` absorbs the
   !! length-of-cropping-cycle output and is discarded. `swhydrlift = 0`
   !! is the conservative default used by hupselbrook/surfacewater cases.
   subroutine read_legacy_cropfixed(icrop, crpfil)
      integer,          intent(in) :: icrop
      character(len=*), intent(in) :: crpfil

      integer :: lcc
      integer :: swhydrlift

      call reset_legacy_crop_globals()

      lcc        = 0
      swhydrlift = 0

      call readcropfixed(icrop, crpfil, lcc, swhydrlift)
   end subroutine read_legacy_cropfixed

   !> Phase 4c-b Task 13: call legacy `readgrass` for rotation entry
   !! `icrop` with crop file `crpfil` (basename, no .crp extension).
   !! Side-effect populates the long list of `variables%` globals
   !! enumerated in the `use variables, only:` clause of `readgrass`
   !! (slatb, amaxtb, kdif, hlim1..hlim4, frtb, fltb, etc.). Local
   !! variables here only exist to absorb the read's output arguments —
   !! they are discarded on return. `swhydrlift = 0` is the conservative
   !! default used by the grassgrowth case.
   subroutine read_legacy_grass(icrop, crpfil)
      integer,          intent(in) :: icrop
      character(len=*), intent(in) :: crpfil

      ! Inputs.
      integer :: swhydrlift

      ! Outputs absorbed and discarded. Types match readgrass's
      ! declarations at src/io/readswap.f90:3470-3473.
      integer :: swharvest, swdmmow, swlossmow, swlossgrz, swdmgrz
      integer :: daylastharvest, maxdaymow, maxdaygrz
      real(8) :: dmharvest, dmlastharvest, dmgrazing, tagprest
      real(8) :: lsdb(100)

      call reset_legacy_crop_globals()

      swhydrlift     = 0
      swharvest      = 0
      swdmmow        = 0
      swlossmow      = 0
      swlossgrz      = 0
      swdmgrz        = 0
      daylastharvest = 0
      maxdaymow      = 0
      maxdaygrz      = 0
      dmharvest      = 0.0d0
      dmlastharvest  = 0.0d0
      dmgrazing      = 0.0d0
      tagprest       = 0.0d0
      lsdb           = 0.0d0

      call readgrass(icrop, crpfil, swharvest, dmharvest, daylastharvest,    &
                     dmlastharvest, swdmmow, maxdaymow, swlossmow, swlossgrz, &
                     swdmgrz, maxdaygrz, dmgrazing, lsdb, tagprest, swhydrlift)
   end subroutine read_legacy_grass

   !> Phase 4d Task 20: parity-test wrapper that exposes the readgrass OUT
   !! args (swharvest, dmharvest, daylastharvest, ..., dmgrazing, tagprest)
   !! to the caller.  These are not in the `variables` module so they
   !! cannot be reached via a `use variables` clause.  The simpler
   !! `read_legacy_grass` discards them; this wrapper returns them via
   !! `intent(out)` arguments so the caller can compare against the
   !! cropgrass_config_t mowing/grazing fields.
   subroutine read_legacy_grass_with_mowing(icrop, crpfil,                    &
                                            swharvest, dmharvest,             &
                                            daylastharvest, dmlastharvest,    &
                                            swdmmow, maxdaymow,               &
                                            swdmgrz, maxdaygrz,               &
                                            dmgrazing, tagprest)
      integer,          intent(in)  :: icrop
      character(len=*), intent(in)  :: crpfil
      integer,          intent(out) :: swharvest, swdmmow, swdmgrz
      integer,          intent(out) :: daylastharvest, maxdaymow, maxdaygrz
      real(8),          intent(out) :: dmharvest, dmlastharvest
      real(8),          intent(out) :: dmgrazing, tagprest

      integer :: swhydrlift, swlossmow, swlossgrz
      real(8) :: lsdb(100)

      call reset_legacy_crop_globals()

      swhydrlift     = 0
      swharvest      = 0
      swdmmow        = 0
      swlossmow      = 0
      swlossgrz      = 0
      swdmgrz        = 0
      daylastharvest = 0
      maxdaymow      = 0
      maxdaygrz      = 0
      dmharvest      = 0.0d0
      dmlastharvest  = 0.0d0
      dmgrazing      = 0.0d0
      tagprest       = 0.0d0
      lsdb           = 0.0d0

      call readgrass(icrop, crpfil, swharvest, dmharvest, daylastharvest,    &
                     dmlastharvest, swdmmow, maxdaymow, swlossmow, swlossgrz, &
                     swdmgrz, maxdaygrz, dmgrazing, lsdb, tagprest, swhydrlift)
   end subroutine read_legacy_grass_with_mowing

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
