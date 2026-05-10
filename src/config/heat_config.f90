!> Heat-flow config — populated from `[heat]` TOML section.
!! Phase 4d Task 6: skeleton + validators only. Mirrors the bottom_boundary
!! sister module: section-not-present sentinel (swhea=0) short-circuits, and
!! a private `check_table_2d` helper guards allocatable tables.
!!
!! Field naming follows legacy `.swp` keys verbatim where possible. The one
!! exception is `tsoil_init`, which holds the initial soil-temperature
!! depth/temp table read by legacy `rdador('zh')` + `rdfdor('tsoil', ...)`
!! pair. The runtime soil-temperature array is `variables%tsoil(numnod)`,
!! so we disambiguate the input table name on the config side.
module heat_config_mod
   use iso_fortran_env, only: real64
   use error_mod, only: error_collection_t, ERR_VALIDATION_OUT_OF_RANGE
   use validation_mod, only: check_int_enum, check_real_range
   implicit none
   private

   public :: heat_config_t

   type :: heat_config_t
      integer :: swhea     = 0  !! 0/1 — heat transport switch
      integer :: swcalt    = 0  !! 1=analytical, 2=numerical (per legacy rdsinr 1..2)
      integer :: swtopbhea = 0  !! 1=air temp, 2=measured surface (per legacy rdsinr 1..2)
      integer :: swbotbhea = 0  !! 1=zero flux, 2=prescribed temp (per legacy rdsinr 1..2)

      ! Soil-texture per layer (used by SWCALT=2 numerical method).
      real(real64), allocatable :: psand(:)       !! Sand fraction per layer (0..1)
      real(real64), allocatable :: psilt(:)       !! Silt fraction per layer (0..1)
      real(real64), allocatable :: pclay(:)       !! Clay fraction per layer (0..1)
      real(real64), allocatable :: porg(:)        !! Organic-matter fraction per layer (0..1)

      ! Initial soil-temperature table — pairs of (depth_cm, temp_C).
      ! Named `tsoil_init` to avoid collision with `variables%tsoil(numnod)`
      ! (the runtime soil-temperature array). Two columns by N rows.
      real(real64), allocatable :: tsoil_init(:,:)

      ! Frost params (legacy rdsdor range -10..5 °C).
      real(real64) :: tfroststa = 0.0_real64
      real(real64) :: tfrostend = 0.0_real64

      ! Phase 0 (SS-HEAT) — promoted from legacy globals to patch the
      ! silent-zero-defaults correctness gap on swcalt=1 (analytical
      ! soil-temperature method). Previously these were never read from
      ! TOML, so TOML runs silently used uninitialised (zero) legacy globals.
      real(real64) :: ddamp  = 0.0_real64    !! Damping depth of temperature wave (L)
      real(real64) :: tmean  = 0.0_real64    !! Prescribed mean annual surface temperature (degC)
      real(real64) :: tampli = 0.0_real64    !! Amplitude of annual surface temperature wave (degC)
      real(real64) :: timref = 0.0_real64    !! Time in year with top of prescribed sine wave (T)

      ! Boundary-condition tables — 2D (n_rows, 2): col 1 = time, col 2 = temperature.
      ! The adapter flattens to interleaved 1D afgen layout: (2*k-1)=time, (2*k)=value.
      real(real64), allocatable :: temtoptab(:,:)   !! Soil-surface temperature vs time
      real(real64), allocatable :: tembtab(:,:)     !! Soil-bottom temperature vs time
   contains
      procedure :: validate => heat_config_validate
      procedure :: finalize => heat_config_finalize
   end type heat_config_t

contains

   !> Local helper: verify allocatable 2D table has expected ncols and >=2 rows.
   !! Skips silently when unallocated (caller decides whether absence is an error).
   !! Mirrors the same-named helper in `bottom_boundary_config_mod`.
   subroutine check_table_2d(table, expected_cols, min_rows, label, errors)
      real(real64), allocatable, intent(in)    :: table(:,:)
      integer,                   intent(in)    :: expected_cols
      integer,                   intent(in)    :: min_rows
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
      if (nrows < min_rows) then
         write(msg, '("nrows=",I0," < required ",I0)') nrows, min_rows
         call errors%append(ERR_VALIDATION_OUT_OF_RANGE, trim(msg), label)
      end if
   end subroutine check_table_2d

   subroutine heat_config_validate(self, errors)
      class(heat_config_t),     intent(in)    :: self
      type(error_collection_t), intent(inout) :: errors
      integer :: nsand, nclay, norg

      ! Section-not-present sentinel: `swhea=0` with all defaults means the
      ! [heat] section was absent in the TOML (reader leaves config at
      ! defaults). Skip validation so existing case TOMLs without the
      ! section still load clean. Phase 4d Tasks 17-19 will add [heat] to
      ! every per-case TOML; case 6 will keep swhea=0 explicitly.
      if (self%swhea == 0) return

      call check_int_enum(self%swhea, [0, 1], 'heat.swhea', errors)
      ! Legacy rdsinr enforces swcalt in 1..2 when swhea=1.
      call check_int_enum(self%swcalt, [1, 2], 'heat.swcalt', errors)
      call check_int_enum(self%swtopbhea, [1, 2], 'heat.swtopbhea', errors)
      call check_int_enum(self%swbotbhea, [1, 2], 'heat.swbotbhea', errors)

      call check_real_range(self%tfroststa, -10.0_real64, 5.0_real64, &
                            'heat.tfroststa', errors)
      call check_real_range(self%tfrostend, -10.0_real64, 5.0_real64, &
                            'heat.tfrostend', errors)

      ! Per-layer texture tables: if any one is allocated, all four must
      ! match in length. Ranges per legacy rdfdor: 0..1 each.
      if (allocated(self%psand) .or. allocated(self%psilt) .or. &
          allocated(self%pclay) .or. allocated(self%porg)) then
         if (.not. allocated(self%psand) .or. &
             .not. allocated(self%psilt) .or. &
             .not. allocated(self%pclay) .or. &
             .not. allocated(self%porg)) then
            call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
               "psand/psilt/pclay/porg must all be set together", 'heat')
         else
            nsand = size(self%psand)
            nclay = size(self%pclay)
            norg  = size(self%porg)
            if (nsand /= size(self%psilt) .or. nsand /= nclay .or. nsand /= norg) then
               call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
                  "psand/psilt/pclay/porg length mismatch", 'heat')
            end if
            call check_fraction_array(self%psand, 'heat.psand', errors)
            call check_fraction_array(self%psilt, 'heat.psilt', errors)
            call check_fraction_array(self%pclay, 'heat.pclay', errors)
            call check_fraction_array(self%porg,  'heat.porg',  errors)
         end if
      end if

      ! Initial soil-temperature table: (depth, temp) pairs. Need >=2 rows
      ! to interpolate; the legacy reader expects matching `zh`/`tsoil`
      ! arrays of length >= 1 but a useful column needs at least two rows.
      call check_table_2d(self%tsoil_init, 2, 2, 'heat.tsoil_init', errors)

      ! Phase 0 (SS-HEAT) — validate promoted swcalt=1 analytics fields.
      call check_real_range(self%ddamp,   0.0_real64, 1.0e6_real64,   'heat.ddamp',  errors)
      call check_real_range(self%tmean, -100.0_real64, 100.0_real64,  'heat.tmean',  errors)
      call check_real_range(self%tampli,  0.0_real64, 100.0_real64,   'heat.tampli', errors)
      call check_real_range(self%timref,  0.0_real64, 366.0_real64,   'heat.timref', errors)

      ! Boundary-condition tables: when allocated, must have exactly 2
      ! columns and at least 1 row (afgen needs at least one pair).
      ! Cross-field requirement checks (e.g. must be present when
      ! swbotbhea=2 / swtopbhea=2) are deferred to a later task.
      call check_table_2d(self%temtoptab, 2, 1, 'heat.temtoptab', errors)
      call check_table_2d(self%tembtab,   2, 1, 'heat.tembtab',   errors)
   end subroutine heat_config_validate

   !> Per-element fraction range check (0..1) for soil-texture arrays.
   subroutine check_fraction_array(arr, label, errors)
      real(real64),             intent(in)    :: arr(:)
      character(len=*),         intent(in)    :: label
      type(error_collection_t), intent(inout) :: errors
      integer :: i
      do i = 1, size(arr)
         call check_real_range(arr(i), 0.0_real64, 1.0_real64, label, errors)
      end do
   end subroutine check_fraction_array

   subroutine heat_config_finalize(self, errors)
      class(heat_config_t),     intent(inout) :: self
      type(error_collection_t), intent(inout) :: errors
      ! No-op for now.
   end subroutine heat_config_finalize

end module heat_config_mod
