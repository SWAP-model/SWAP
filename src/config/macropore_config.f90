!> Macropore config — populated from `[macropore]` TOML section.
!! Phase 4f-prep Task A1: 22-key macropore block (13 scalars + 9
!! per-layer table fields). Validators switch-gated on `swmacro`.
!! See docs/phase-4e-macroporeflow-audit.md for the field catalogue.
module macropore_config_mod
   use iso_fortran_env, only: real64
   use error_mod, only: error_collection_t, ERR_VALIDATION_OUT_OF_RANGE
   use validation_mod, only: check_int_enum, check_int_range, check_real_range
   implicit none
   private

   public :: macropore_config_t

   type :: macropore_config_t
      ! Scalars (13)
      integer      :: swmacro   = 0
      real(real64) :: z_ah      = 0.0_real64
      real(real64) :: z_ic      = 0.0_real64
      real(real64) :: z_st      = 0.0_real64
      real(real64) :: vlmpstss  = 0.0_real64
      real(real64) :: ppicss    = 0.0_real64
      integer      :: numsbdm   = 0
      real(real64) :: powm      = 0.0_real64
      real(real64) :: rzah      = 0.0_real64
      real(real64) :: spoint    = 0.0_real64
      integer      :: swpowm    = 0
      real(real64) :: dipomi    = 0.0_real64
      real(real64) :: dipoma    = 0.0_real64

      ! Per-layer table fields (9)
      integer,      allocatable :: swsoilshr(:)
      integer,      allocatable :: swshrinp(:)
      real(real64), allocatable :: thetcrmp(:)
      real(real64), allocatable :: geomfac(:)
      real(real64), allocatable :: shrpar(:,:)
      integer,      allocatable :: swsorp(:)
      real(real64), allocatable :: sorpfacparl(:)
      real(real64), allocatable :: sorpmax(:)
      real(real64), allocatable :: sorpalfa(:)
   contains
      procedure :: validate => macropore_config_validate
      procedure :: finalize => macropore_config_finalize
   end type macropore_config_t

contains

   !> Local helper: confirm a 1D allocatable has at least one element.
   !! Skips silently when unallocated (caller decides whether absence is
   !! itself an error). Mirrors the bottom_boundary_config 2D helper.
   subroutine check_array_1d_size(n, label, errors)
      integer,                  intent(in)    :: n
      character(len=*),         intent(in)    :: label
      type(error_collection_t), intent(inout) :: errors
      character(len=64) :: msg
      if (n < 1) then
         write(msg, '("size=",I0," <1")') n
         call errors%append(ERR_VALIDATION_OUT_OF_RANGE, trim(msg), label)
      end if
   end subroutine check_array_1d_size

   !> Local helper: per-layer integer enum check across a 1D allocatable.
   subroutine check_int_enum_per_layer(values, allowed, label, errors)
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
   end subroutine check_int_enum_per_layer

   subroutine macropore_config_validate(self, errors)
      class(macropore_config_t), intent(in)    :: self
      type(error_collection_t),  intent(inout) :: errors

      ! swmacro is always validated; values 0/1 cover legacy 0..2 range
      ! collapsed to "off / on" per audit doc (case 3 .swp value).
      call check_int_enum(self%swmacro, [0, 1], 'macropore.swmacro', errors)

      if (self%swmacro /= 1) return

      ! Scalars (12 active under swmacro=1; powm/rzah/spoint/swpowm are
      ! optional in legacy but ranges still apply when the config carries
      ! them).
      call check_real_range(self%z_ah,    -1000.0_real64,  0.0_real64,    'macropore.z_ah',     errors)
      call check_real_range(self%z_ic,    -1000.0_real64,  0.0_real64,    'macropore.z_ic',     errors)
      call check_real_range(self%z_st,    -1000.0_real64,  0.0_real64,    'macropore.z_st',     errors)
      call check_real_range(self%vlmpstss, 0.0_real64,     0.5_real64,    'macropore.vlmpstss', errors)
      call check_real_range(self%ppicss,   0.0_real64,     0.99_real64,   'macropore.ppicss',   errors)
      call check_int_range (self%numsbdm,  0,              50,            'macropore.numsbdm',  errors)
      call check_real_range(self%powm,     0.0_real64,     100.0_real64,  'macropore.powm',     errors)
      call check_real_range(self%rzah,     0.0_real64,     1.0_real64,    'macropore.rzah',     errors)
      call check_real_range(self%spoint,   0.0_real64,     1.0_real64,    'macropore.spoint',   errors)
      call check_int_enum  (self%swpowm,   [0, 1],                        'macropore.swpowm',   errors)
      call check_real_range(self%dipomi,   0.1_real64,     1000.0_real64, 'macropore.dipomi',   errors)
      call check_real_range(self%dipoma,   0.1_real64,     1000.0_real64, 'macropore.dipoma',   errors)

      ! Per-layer 1D arrays — if allocated, size must be >= 1.
      ! nlayers cross-check deferred to finalize stage.
      if (allocated(self%swsoilshr)) then
         call check_array_1d_size(size(self%swsoilshr), 'macropore.swsoilshr', errors)
         call check_int_enum_per_layer(self%swsoilshr, [0, 1, 2], 'macropore.swsoilshr', errors)
      end if
      if (allocated(self%swshrinp)) then
         call check_array_1d_size(size(self%swshrinp), 'macropore.swshrinp', errors)
         call check_int_enum_per_layer(self%swshrinp, [1, 2], 'macropore.swshrinp', errors)
      end if
      if (allocated(self%thetcrmp)) then
         call check_array_1d_size(size(self%thetcrmp), 'macropore.thetcrmp', errors)
      end if
      if (allocated(self%geomfac)) then
         call check_array_1d_size(size(self%geomfac), 'macropore.geomfac', errors)
      end if
      if (allocated(self%swsorp)) then
         call check_array_1d_size(size(self%swsorp), 'macropore.swsorp', errors)
         call check_int_enum_per_layer(self%swsorp, [1, 2], 'macropore.swsorp', errors)
      end if
      if (allocated(self%sorpfacparl)) then
         call check_array_1d_size(size(self%sorpfacparl), 'macropore.sorpfacparl', errors)
      end if
      if (allocated(self%sorpmax)) then
         call check_array_1d_size(size(self%sorpmax), 'macropore.sorpmax', errors)
      end if
      if (allocated(self%sorpalfa)) then
         call check_array_1d_size(size(self%sorpalfa), 'macropore.sorpalfa', errors)
      end if

      ! shrpar 2D table — width is 5 columns (SHRPARA..SHRPARE per the
      ! legacy .swp documentation; column population varies by
      ! SWSOILSHR/SWSHRINP combination).
      if (allocated(self%shrpar)) then
         call check_array_1d_size(size(self%shrpar, 1), 'macropore.shrpar', errors)
         if (size(self%shrpar, 2) /= 5) then
            call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
               "shrpar ncols /= 5", 'macropore.shrpar')
         end if
      end if
   end subroutine macropore_config_validate

   subroutine macropore_config_finalize(self, errors)
      class(macropore_config_t), intent(inout) :: self
      type(error_collection_t),  intent(inout) :: errors
      ! Placeholder for downstream cross-section consistency checks
      ! (e.g., per-layer table sizes vs nlayers from soil_config_t).
      return
   end subroutine macropore_config_finalize

end module macropore_config_mod
