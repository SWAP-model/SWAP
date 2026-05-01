!> Surface-water management config — populated from `[surface_water]`
!! TOML section. Phase 4f-prep Task B1.
module surface_water_config_mod
   use iso_fortran_env, only: real64
   use error_mod, only: error_collection_t, ERR_VALIDATION_OUT_OF_RANGE, &
                         ERR_PARSE_MISSING_REQUIRED, ERR_VALIDATION_CROSS_FIELD
   use validation_mod, only: check_int_enum, check_int_range, check_real_range
   implicit none
   private

   public :: surface_water_config_t

   type :: surface_water_config_t
      integer      :: swsrf  = 1
      integer      :: swsec  = 1
      real(real64) :: wlact  = 0.0_real64
      real(real64) :: osswlm = 0.0_real64
      integer      :: nmper  = 0

      real(real64), allocatable :: impend(:)
      integer,      allocatable :: swman(:)
      real(real64), allocatable :: wscap(:)
      real(real64), allocatable :: wldip(:)
      integer,      allocatable :: intwl(:)

      integer      :: swqhr  = 1
      real(real64) :: sofcu  = 0.0_real64

      real(real64), allocatable :: hbweir(:)
      real(real64), allocatable :: alphaw(:)
      real(real64), allocatable :: betaw(:)
   contains
      procedure :: validate => surface_water_config_validate
      procedure :: finalize => surface_water_config_finalize
   end type surface_water_config_t

contains

   !> Cross-check: if a 1D real array is allocated, its size must match
   !! `nmper`; if not allocated, it is a missing-required error. Mirrors
   !! the per-period sizing convention in macropore_config_mod.
   subroutine check_table_1d_size_matches_real(arr, nmper, label, errors)
      real(real64), allocatable, intent(in)    :: arr(:)
      integer,                   intent(in)    :: nmper
      character(len=*),          intent(in)    :: label
      type(error_collection_t),  intent(inout) :: errors
      character(len=128) :: msg
      if (.not. allocated(arr)) then
         call errors%append(ERR_PARSE_MISSING_REQUIRED, &
                            "required when swsec=2", label)
         return
      end if
      if (size(arr) /= nmper) then
         write(msg, '("size=",I0," /= nmper=",I0)') size(arr), nmper
         call errors%append(ERR_VALIDATION_OUT_OF_RANGE, trim(msg), label)
      end if
   end subroutine check_table_1d_size_matches_real

   !> Integer-typed counterpart to check_table_1d_size_matches_real.
   subroutine check_table_1d_size_matches_int(arr, nmper, label, errors)
      integer, allocatable,     intent(in)    :: arr(:)
      integer,                  intent(in)    :: nmper
      character(len=*),         intent(in)    :: label
      type(error_collection_t), intent(inout) :: errors
      character(len=128) :: msg
      if (.not. allocated(arr)) then
         call errors%append(ERR_PARSE_MISSING_REQUIRED, &
                            "required when swsec=2", label)
         return
      end if
      if (size(arr) /= nmper) then
         write(msg, '("size=",I0," /= nmper=",I0)') size(arr), nmper
         call errors%append(ERR_VALIDATION_OUT_OF_RANGE, trim(msg), label)
      end if
   end subroutine check_table_1d_size_matches_int

   !> Per-period integer enum check.
   subroutine check_int_enum_per_period(values, allowed, label, errors)
      integer,                  intent(in)    :: values(:)
      integer,                  intent(in)    :: allowed(:)
      character(len=*),         intent(in)    :: label
      type(error_collection_t), intent(inout) :: errors
      integer :: i
      character(len=128) :: ctx
      do i = 1, size(values)
         write(ctx, '(A,"[",I0,"]")') label, i
         call check_int_enum(values(i), allowed, trim(ctx), errors)
      end do
   end subroutine check_int_enum_per_period

   subroutine surface_water_config_validate(self, errors)
      class(surface_water_config_t), intent(in)    :: self
      type(error_collection_t),      intent(inout) :: errors

      ! Stub-error for swsrf=3 (primary system) must run BEFORE the swsrf=2
      ! early-return below; otherwise swsrf=3 would short-circuit without error.
      ! The remaining stub-errors (swsec=1, swqhr=2, swman=2) run AFTER the
      ! early-return: swsec and swqhr have default values of 1, which would
      ! cause false positives on every swsrf=1 config if checked here.
      if (self%swsrf == 3) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'surface_water.swsrf=3 (primary system) not yet supported in ' // &
            'the TOML pipeline; use the legacy executable.', 'surface_water')
      end if

      call check_int_enum(self%swsrf, [1, 2], 'surface_water.swsrf', errors)
      if (self%swsrf /= 2) return

      ! Stub-errors for swsec/swqhr/swman branches not yet supported in the
      ! TOML pipeline; only reached when swsrf=2 (active surface-water mgmt).
      if (self%swsec == 1) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'surface_water.swsec=1 (input water level) not yet supported ' // &
            'in the TOML pipeline; use the legacy executable.', 'surface_water')
      end if
      if (self%swqhr == 2) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'surface_water.swqhr=2 (q-h table discharge) not yet supported ' // &
            'in the TOML pipeline; use the legacy executable.', 'surface_water')
      end if
      if (allocated(self%swman)) then
         if (any(self%swman == 2)) then
            call errors%append(ERR_VALIDATION_CROSS_FIELD, &
               'surface_water.swman=2 (automatic weir) not yet supported ' // &
               'in the TOML pipeline; use the legacy executable.', 'surface_water')
         end if
      end if

      call check_int_enum(self%swsec, [1, 2], 'surface_water.swsec', errors)
      if (self%swsec /= 2) return

      call check_real_range(self%wlact,  -1.0e4_real64, 1.0e4_real64, &
                            'surface_water.wlact',  errors)
      call check_real_range(self%osswlm, 0.0_real64,    10.0_real64, &
                            'surface_water.osswlm', errors)
      call check_int_range (self%nmper,  1, 3660, 'surface_water.nmper', errors)

      if (self%nmper >= 1) then
         call check_table_1d_size_matches_real(self%impend, self%nmper, &
                                               'surface_water.impend', errors)
         call check_table_1d_size_matches_int (self%swman,  self%nmper, &
                                               'surface_water.swman',  errors)
         call check_table_1d_size_matches_real(self%wscap,  self%nmper, &
                                               'surface_water.wscap',  errors)
         call check_table_1d_size_matches_real(self%wldip,  self%nmper, &
                                               'surface_water.wldip',  errors)
         call check_table_1d_size_matches_int (self%intwl,  self%nmper, &
                                               'surface_water.intwl',  errors)

         if (allocated(self%swman) .and. size(self%swman) == self%nmper) then
            call check_int_enum_per_period(self%swman, [1, 2], &
                                           'surface_water.swman', errors)
         end if
      end if

      call check_int_enum(self%swqhr, [1, 2], 'surface_water.swqhr', errors)
      call check_real_range(self%sofcu, 0.1_real64, 1.0e5_real64, &
                            'surface_water.sofcu', errors)

      if (self%swqhr == 1 .and. self%nmper >= 1) then
         call check_table_1d_size_matches_real(self%hbweir, self%nmper, &
                                               'surface_water.hbweir', errors)
         call check_table_1d_size_matches_real(self%alphaw, self%nmper, &
                                               'surface_water.alphaw', errors)
         call check_table_1d_size_matches_real(self%betaw,  self%nmper, &
                                               'surface_water.betaw',  errors)
      end if
   end subroutine surface_water_config_validate

   subroutine surface_water_config_finalize(self, errors)
      ! Phase 4f-prep Task D2: legacy `rddre` (src/io/readswap.f90:4598-4599)
      ! normalizes the per-period exponential-weir alpha coefficient when
      ! SWQHR=1:
      !   alphaw(i) = alphaw(i) * 8.64 * 100^(1-betaw(i)) / sofcu
      ! TOML stores the raw input value; finalize applies the same
      ! normalization so config%surface_water%alphaw matches the legacy
      ! global post-rddre. Skip when swqhr /= 1 or arrays unallocated.
      class(surface_water_config_t), intent(inout) :: self
      type(error_collection_t),      intent(inout) :: errors
      integer :: i
      if (self%swqhr /= 1) return
      if (.not. allocated(self%alphaw) .or. .not. allocated(self%betaw)) return
      if (self%sofcu <= 0.0_real64) return
      do i = 1, self%nmper
         self%alphaw(i) = self%alphaw(i) * 8.64_real64 &
            * 100.0_real64 ** (1.0_real64 - self%betaw(i)) / self%sofcu
      end do
   end subroutine surface_water_config_finalize

end module surface_water_config_mod
