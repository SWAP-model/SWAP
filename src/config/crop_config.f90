!> [crop] section config: rotation list + per-crop file references.
module crop_config_mod
   use iso_fortran_env, only: real64
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD, ERR_VALIDATION_OUT_OF_RANGE
   use validation_mod, only: check_int_enum
   implicit none
   private

   public :: crop_config_t

   type :: crop_config_t
      integer :: swcrop = 0

      !> Per-entry arrays sized to the rotation length. All same length
      !! after parse; enforced in validate.
      real(real64),     allocatable :: rotation_start(:)
      real(real64),     allocatable :: rotation_end(:)
      character(len=256), allocatable :: rotation_file(:)
      integer,          allocatable :: rotation_type(:)
   contains
      procedure :: validate => crop_config_validate
      procedure :: finalize => crop_config_finalize
   end type crop_config_t

contains

   subroutine crop_config_validate(self, errors)
      class(crop_config_t),     intent(in)    :: self
      type(error_collection_t), intent(inout) :: errors

      integer :: i, n

      call check_int_enum(self%swcrop, [0, 1], "crop.swcrop", errors)

      ! If crops are disabled, skip rotation checks.
      if (self%swcrop == 0) return

      ! Rotation must have at least one entry.
      if (.not. allocated(self%rotation_start) .or. size(self%rotation_start) == 0) then
         call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
                            "crop rotation must have at least one entry when swcrop=1", &
                            "crop.rotation")
         return
      end if

      ! All rotation arrays must have the same length.
      n = size(self%rotation_start)
      if (size(self%rotation_end)  /= n .or. &
          size(self%rotation_file) /= n .or. &
          size(self%rotation_type) /= n) then
         call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
                            "crop rotation parallel arrays differ in length", &
                            "crop.rotation")
         return
      end if

      ! Per-entry checks.
      do i = 1, n
         call check_int_enum(self%rotation_type(i), [1, 2, 3], &
                             "crop.rotation.type", errors)
         if (self%rotation_start(i) >= self%rotation_end(i)) then
            call errors%append(ERR_VALIDATION_CROSS_FIELD, &
                               "crop rotation entry: start must be before end", &
                               "crop.rotation")
         end if
      end do
   end subroutine crop_config_validate

   subroutine crop_config_finalize(self, errors)
      class(crop_config_t),     intent(inout) :: self
      type(error_collection_t), intent(inout) :: errors
   end subroutine crop_config_finalize

end module crop_config_mod
