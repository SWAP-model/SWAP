!> Irrigation config — Phase 4d Task 9 skeleton + Task 10 validators.
!! Two public types in one module: `irrigation_config_t` holds top-level
!! `.swp` fixed-irrigation parameters; `irrigation_schedule_t` is a nested
!! sub-type embedded in each per-crop config (cropfixed, cropwofost,
!! cropgrass) for `.crp` scheduling. Readers land in Task 11; swap_config
!! wiring in Task 12. Field names mirror the legacy `.swp`/`.crp` keys
!! verbatim.
!!
!! Validator design mirrors `bottom_boundary_config_mod` and
!! `heat_config_mod`: a section-not-present sentinel short-circuits
!! (swirfix=0 top-level; schedule=0 nested), then enum/range/branch
!! checks fan out per the legacy `rdsinr`/`rdsdor` ranges in
!! `src/io/readswap.f90` and `src/crop/irrigation.f90`.
module irrigation_config_mod
   use iso_fortran_env, only: real64
   use error_mod, only: error_collection_t, ERR_VALIDATION_OUT_OF_RANGE, &
                        ERR_VALIDATION_CROSS_FIELD
   use validation_mod, only: check_int_enum, check_int_range, check_real_range, &
                             check_ordered_pair
   implicit none
   private

   public :: irrigation_config_t
   public :: irrigation_schedule_t

   ! Top-level fixed-irrigation (.swp side).
   type :: irrigation_config_t
      integer                       :: swirfix = 0
      integer                       :: swssdi  = 0       !! SSDI sub-surface drip irrigation. 0=off (only supported value); 1 stub-errored (SS-10.5).
      character(len=:), allocatable :: irgfil
      ! Phase 4f cleanup: CSV companion file path (relative to swap.toml)
      ! holding a long-form (date, depth, conc, type) fixed-events table.
      ! Replaces the legacy `.irg` external-file HACK. Mutually exclusive
      ! with `fixed_events` (inline TOML). When set, the adapter reads
      ! the CSV via `csv_reader_mod%read_csv_table` and unpacks it
      ! into the same legacy `irdate/irdepth/irconc/irtype` arrays.
      character(len=:), allocatable :: fixed_events_file
      real(real64),     allocatable :: fixed_events(:,:)
   contains
      procedure :: validate => irrigation_config_validate
      procedure :: finalize => irrigation_config_finalize
   end type irrigation_config_t

   ! Nested per-crop scheduling (.crp side).
   type :: irrigation_schedule_t
      integer      :: schedule = 0
      integer      :: startirr_day = 0
      integer      :: startirr_month = 0
      integer      :: endirr_day = 0
      integer      :: endirr_month = 0
      real(real64) :: cirrs = 0.0_real64
      integer      :: isuas = 0
      integer      :: tcs = 0
      integer      :: dcs = 0
      real(real64), allocatable :: trel_table(:,:)
      real(real64), allocatable :: raw_table(:,:)
      real(real64), allocatable :: taw_table(:,:)
      real(real64), allocatable :: dwa_table(:,:)
      real(real64) :: irgthreshold = 0.0_real64
      real(real64), allocatable :: hcri_table(:,:)
      real(real64), allocatable :: tcri_table(:,:)
      real(real64) :: dcrit = 0.0_real64
      integer      :: swcirrthres = 0
      real(real64) :: cirrthres = 0.0_real64
      real(real64) :: perirrsurp = 0.0_real64
      integer      :: tcsfix = 0
      integer      :: irgdayfix = 0
      real(real64) :: phfieldcapacity = 0.0_real64
      real(real64), allocatable :: di_table(:,:)
      real(real64) :: raithreshold = 0.0_real64
      real(real64), allocatable :: fid_table(:,:)
      integer      :: dcslim = 0
      real(real64) :: irgdepmin = 0.0_real64
      real(real64) :: irgdepmax = 0.0_real64
   contains
      procedure :: validate => irrigation_schedule_validate
      procedure :: finalize => irrigation_schedule_finalize
   end type irrigation_schedule_t

contains

   !> Local helper: verify allocatable 2D table has expected ncols and >=1 row.
   !! Skips silently when unallocated (caller decides whether absence is an error).
   !! Mirrors the same-named helper in `bottom_boundary_config_mod`.
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

   subroutine irrigation_config_validate(self, errors)
      class(irrigation_config_t), intent(in)    :: self
      type(error_collection_t),   intent(inout) :: errors
      logical :: have_irgfil, have_table, have_csv

      ! Phase 4f-extend SS-10.5 stub-error: SSDI sub-surface drip
      ! irrigation. SSDI_irrigation(1) (irrigation.f90:580) was a
      ! legacy reader still reachable from the production runtime.
      ! SS-10.5 retired its RDinit(swpfile) call and now reads swssdi
      ! from this schema slot. No regression case authors swssdi=1; if
      ! a future case needs it, the read_ssdi_input() block must also
      ! be ported to schema (ssdi_file companion table).
      if (self%swssdi == 1) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'irrigation.swssdi=1 (sub-surface drip irrigation) not yet ' // &
            'supported in the TOML pipeline; no regression case ' // &
            'exercises it. Port read_ssdi_input() to schema if a case ' // &
            'needs this.', 'irrigation')
      end if
      call check_int_enum(self%swssdi, [0, 1], 'irrigation.swssdi', errors)

      ! Sentinel: swirfix=0 ⇒ no fixed irrigation at the .swp level. Skip
      ! everything so existing case TOMLs without an [irrigation] section
      ! still load clean.
      if (self%swirfix == 0) return

      call check_int_enum(self%swirfix, [0, 1], 'irrigation.swirfix', errors)

      if (self%swirfix == 1) then
         have_irgfil = allocated(self%irgfil)
         if (have_irgfil) have_irgfil = len_trim(self%irgfil) > 0
         have_csv = allocated(self%fixed_events_file)
         if (have_csv) have_csv = len_trim(self%fixed_events_file) > 0
         have_table = allocated(self%fixed_events)
         if (.not. have_irgfil .and. .not. have_table .and. .not. have_csv) then
            call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
               "swirfix=1 requires irgfil, fixed_events, or fixed_events_file", &
               'irrigation')
         end if
         ! Phase 4f cleanup: at most one source for the events.
         if (count([have_irgfil, have_table, have_csv]) > 1) then
            call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
               "swirfix=1 cannot set more than one of " // &
               "{irgfil, fixed_events, fixed_events_file}", 'irrigation')
         end if
         call check_table_2d(self%fixed_events, 4, 'irrigation.fixed_events', errors)
      end if
   end subroutine irrigation_config_validate

   subroutine irrigation_config_finalize(self, errors)
      class(irrigation_config_t), intent(inout) :: self
      type(error_collection_t),   intent(inout) :: errors
      ! No-op.
   end subroutine irrigation_config_finalize

   subroutine irrigation_schedule_validate(self, errors)
      class(irrigation_schedule_t), intent(in)    :: self
      type(error_collection_t),     intent(inout) :: errors

      ! Sentinel: schedule=0 ⇒ no scheduled irrigation for this crop. Mirrors
      ! the legacy `if (schedule .eq. 1)` gating in src/crop/irrigation.f90.
      if (self%schedule == 0) return

      call check_int_enum(self%schedule, [0, 1], 'irrigation_schedule.schedule', errors)
      call check_int_enum(self%isuas, [0, 1], 'irrigation_schedule.isuas', errors)
      ! Legacy rdsinr enforces tcs in 1..8; tcs=5 is obsolete and rejected at
      ! runtime by src/crop/irrigation.f90. We accept the legacy range here
      ! and let the runtime reject 5 — matches existing parity behaviour.
      call check_int_range(self%tcs, 1, 8, 'irrigation_schedule.tcs', errors)
      call check_int_enum(self%dcs, [1, 2], 'irrigation_schedule.dcs', errors)
      call check_int_enum(self%swcirrthres, [0, 1], 'irrigation_schedule.swcirrthres', errors)
      call check_int_enum(self%tcsfix, [0, 1], 'irrigation_schedule.tcsfix', errors)
      call check_int_enum(self%dcslim, [0, 1], 'irrigation_schedule.dcslim', errors)

      ! Day/month windows.
      call check_int_range(self%startirr_day,   1, 31, &
         'irrigation_schedule.startirr_day', errors)
      call check_int_range(self%startirr_month, 1, 12, &
         'irrigation_schedule.startirr_month', errors)
      call check_int_range(self%endirr_day,     1, 31, &
         'irrigation_schedule.endirr_day', errors)
      call check_int_range(self%endirr_month,   1, 12, &
         'irrigation_schedule.endirr_month', errors)

      ! Real ranges (per legacy rdsdor in src/crop/irrigation.f90).
      call check_real_range(self%cirrs, 0.0_real64, 100.0_real64, &
         'irrigation_schedule.cirrs', errors)
      call check_real_range(self%irgthreshold, 0.0_real64, 20.0_real64, &
         'irrigation_schedule.irgthreshold', errors)
      call check_real_range(self%dcrit, -100.0_real64, 0.0_real64, &
         'irrigation_schedule.dcrit', errors)
      call check_real_range(self%cirrthres, 0.0_real64, 100.0_real64, &
         'irrigation_schedule.cirrthres', errors)
      call check_real_range(self%perirrsurp, 0.0_real64, 100.0_real64, &
         'irrigation_schedule.perirrsurp', errors)
      call check_real_range(self%phfieldcapacity, -1000.0_real64, 0.0_real64, &
         'irrigation_schedule.phfieldcapacity', errors)
      call check_real_range(self%raithreshold, 0.0_real64, 1000.0_real64, &
         'irrigation_schedule.raithreshold', errors)
      call check_real_range(self%irgdepmin, 0.0_real64, 1.0e7_real64, &
         'irrigation_schedule.irgdepmin', errors)
      call check_real_range(self%irgdepmax, 0.0_real64, 1.0e7_real64, &
         'irrigation_schedule.irgdepmax', errors)

      ! Branch validation per timing-criterion (tcs).
      select case (self%tcs)
      case (1)
         if (.not. allocated(self%trel_table)) then
            call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
               "tcs=1 requires trel_table", 'irrigation_schedule')
         end if
         call check_table_2d(self%trel_table, 2, 'irrigation_schedule.trel_table', errors)
      case (2)
         if (.not. allocated(self%raw_table)) then
            call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
               "tcs=2 requires raw_table", 'irrigation_schedule')
         end if
         call check_table_2d(self%raw_table, 2, 'irrigation_schedule.raw_table', errors)
      case (3)
         if (.not. allocated(self%taw_table)) then
            call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
               "tcs=3 requires taw_table", 'irrigation_schedule')
         end if
         call check_table_2d(self%taw_table, 2, 'irrigation_schedule.taw_table', errors)
      case (4)
         if (.not. allocated(self%dwa_table)) then
            call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
               "tcs=4 requires dwa_table", 'irrigation_schedule')
         end if
         call check_table_2d(self%dwa_table, 2, 'irrigation_schedule.dwa_table', errors)
      case (7)
         if (.not. allocated(self%hcri_table)) then
            call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
               "tcs=7 requires hcri_table", 'irrigation_schedule')
         end if
         call check_table_2d(self%hcri_table, 2, 'irrigation_schedule.hcri_table', errors)
      case (8)
         if (.not. allocated(self%tcri_table)) then
            call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
               "tcs=8 requires tcri_table", 'irrigation_schedule')
         end if
         call check_table_2d(self%tcri_table, 2, 'irrigation_schedule.tcri_table', errors)
      end select

      ! Branch validation per depth-criterion (dcs).
      select case (self%dcs)
      case (1)
         if (.not. allocated(self%di_table)) then
            call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
               "dcs=1 requires di_table", 'irrigation_schedule')
         end if
         call check_table_2d(self%di_table, 2, 'irrigation_schedule.di_table', errors)
      case (2)
         if (.not. allocated(self%fid_table)) then
            call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
               "dcs=2 requires fid_table", 'irrigation_schedule')
         end if
         call check_table_2d(self%fid_table, 2, 'irrigation_schedule.fid_table', errors)
      end select

      ! tcsfix=1 + tcs=6 is rejected at runtime — mirror that here.
      if (self%tcsfix == 1) then
         call check_int_range(self%irgdayfix, 1, 366, &
            'irrigation_schedule.irgdayfix', errors)
         if (self%tcs == 6) then
            call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
               "tcsfix=1 with tcs=6 not allowed", 'irrigation_schedule')
         end if
      end if

      ! Range relationship guarded by dcslim=1.
      if (self%dcslim == 1) then
         call check_ordered_pair(self%irgdepmin, self%irgdepmax, &
            'irgdepmin', 'irgdepmax', 'irrigation_schedule', errors)
      end if
   end subroutine irrigation_schedule_validate

   subroutine irrigation_schedule_finalize(self, errors)
      class(irrigation_schedule_t), intent(inout) :: self
      type(error_collection_t),     intent(inout) :: errors
      ! No-op.
   end subroutine irrigation_schedule_finalize

end module irrigation_config_mod
