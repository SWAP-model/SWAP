!> Solute-transport config — populated from `[solute]` TOML section.
!! Phase 4d Task 13: skeleton + validators. Mirrors the heat /
!! bottom_boundary sister modules: a section-not-present sentinel
!! (`swsolu=0`) short-circuits validation, and a private
!! `check_table_2d` helper guards the `pertabsolu` allocatable.
!!
!! Field naming follows legacy `.swp` keys verbatim (see `readswap.f90`
!! around the `swsolu` block: `cdrain`, `cseep`, `tscf`, `ldis`,
!! `rtheta`, `bexp`, `swdc`, `swbotbc`). The only field whose name is
!! not pulled from a single legacy global is `pertabsolu` — a per-layer
!! decomposition table; legacy stores it spread across `decpot` /
!! `fdepth` arrays, but we keep the spec's grouped name.
!!
!! `swsoltyp` is a per-spec D4 placeholder switch (no legacy global
!! today) used to discriminate solute model variants.
module solute_config_mod
   use iso_fortran_env, only: real64
   use error_mod, only: error_collection_t, ERR_VALIDATION_OUT_OF_RANGE
   use validation_mod, only: check_int_enum, check_real_range, check_nonnegative_real
   implicit none
   private

   public :: solute_config_t

   type :: solute_config_t
      integer :: swsolu  = 0          !! 0/1 — solute transport switch
      integer :: swbotbc = 0          !! Solute boundary at bottom (0..2)
      real(real64) :: cdrain = 0.0_real64
      real(real64) :: cseep  = 0.0_real64
      real(real64) :: tscf   = 0.0_real64    !! Transpiration stream conc factor
      real(real64) :: ldis   = 0.0_real64    !! Dispersion length

      ! Root-uptake (decomposition switch=1 inputs):
      real(real64) :: rtheta = 0.0_real64    !! Water-content threshold
      real(real64) :: bexp   = 0.0_real64    !! Exponent

      ! Salinity (independent of crop-side ecmax/ecslop):
      real(real64) :: ecmax  = 0.0_real64
      real(real64) :: ecslop = 0.0_real64

      ! Decomposition:
      integer :: swsoltyp = 0
      integer :: swdc     = 0

      ! Per-layer decomposition table — (depth, factor) pairs, two columns.
      real(real64), allocatable :: pertabsolu(:,:)
   contains
      procedure :: validate => solute_config_validate
      procedure :: finalize => solute_config_finalize
   end type solute_config_t

contains

   !> Local helper: verify allocatable 2D table has expected ncols and >=1 row.
   !! Skips silently when unallocated. Mirrors the same-named helper in
   !! `bottom_boundary_config_mod` and `heat_config_mod`.
   subroutine check_table_2d(table, expected_cols, label, errors)
      real(real64), allocatable, intent(in)    :: table(:,:)
      integer,                   intent(in)    :: expected_cols
      character(len=*),          intent(in)    :: label
      type(error_collection_t),  intent(inout) :: errors
      character(len=128) :: msg
      integer :: nrows, ncols
      if (.not. allocated(table)) return
      nrows = size(table, 1)
      ncols = size(table, 2)
      if (ncols /= expected_cols) then
         write(msg, '("ncols=",I0," expected ",I0)') ncols, expected_cols
         call errors%append(ERR_VALIDATION_OUT_OF_RANGE, trim(msg), label)
      end if
      if (nrows < 1) then
         call errors%append(ERR_VALIDATION_OUT_OF_RANGE, "nrows<1", label)
      end if
   end subroutine check_table_2d

   subroutine solute_config_validate(self, errors)
      class(solute_config_t),   intent(in)    :: self
      type(error_collection_t), intent(inout) :: errors

      ! Section-not-present sentinel: `swsolu=0` with all defaults means the
      ! [solute] section was absent in the TOML (reader leaves config at
      ! defaults). Skip validation so existing case TOMLs without the
      ! section still load clean. Phase 4d Tasks 17-19 will add [solute] to
      ! every per-case TOML; cases 2/4/6 keep swsolu=0 explicitly.
      if (self%swsolu == 0) return

      call check_int_enum(self%swsolu,  [0, 1],       'solute.swsolu',  errors)
      ! Legacy rdsinr enforces swbotbc in 0..2 when swsolu=1.
      call check_int_enum(self%swbotbc, [0, 1, 2],    'solute.swbotbc', errors)

      ! Concentrations and uptake/dispersion params (legacy rdsdor ranges).
      call check_real_range(self%cdrain, 0.0_real64, 100.0_real64, &
                            'solute.cdrain', errors)
      call check_real_range(self%cseep,  0.0_real64, 100.0_real64, &
                            'solute.cseep',  errors)
      call check_real_range(self%tscf,   0.0_real64,  10.0_real64, &
                            'solute.tscf',   errors)
      call check_real_range(self%ldis,   0.0_real64, 100.0_real64, &
                            'solute.ldis',   errors)

      ! Root-uptake (legacy rdsdor ranges under swdc=1):
      call check_real_range(self%rtheta, 0.0_real64,   0.4_real64, &
                            'solute.rtheta', errors)
      call check_real_range(self%bexp,   0.0_real64,   2.0_real64, &
                            'solute.bexp',   errors)

      ! Salinity scalars — non-negative; cross-check against the per-crop
      ! salinity stress block is intentionally skipped (see spec D4).
      call check_nonnegative_real(self%ecmax,  'solute.ecmax',  errors)
      call check_nonnegative_real(self%ecslop, 'solute.ecslop', errors)

      ! Decomposition switches.
      call check_int_enum(self%swsoltyp, [0, 1], 'solute.swsoltyp', errors)
      call check_int_enum(self%swdc,     [0, 1], 'solute.swdc',     errors)

      ! Per-layer decomposition table: (depth, factor) pairs.
      call check_table_2d(self%pertabsolu, 2, 'solute.pertabsolu', errors)
   end subroutine solute_config_validate

   subroutine solute_config_finalize(self, errors)
      class(solute_config_t),   intent(inout) :: self
      type(error_collection_t), intent(inout) :: errors
      ! No-op for now.
   end subroutine solute_config_finalize

end module solute_config_mod
