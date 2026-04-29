!> [crop] section config: rotation list + per-crop file references.
module crop_config_mod
   use iso_fortran_env, only: real64
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD, ERR_VALIDATION_OUT_OF_RANGE
   use validation_mod, only: check_int_enum
   use cropfixed_config_mod, only: cropfixed_config_t
   use cropgrass_config_mod, only: cropgrass_config_t
   use cropwofost_config_mod, only: cropwofost_config_t
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
      integer,          allocatable :: rotation_swhydrlift(:)   ! 0 or 1, per entry

      !> Per-entry sub-configs, populated when file= reference is followed.
      !! rotation_fixed(i) is meaningful when rotation_type(i) == 1.
      !! rotation_grass(i) is meaningful when rotation_type(i) == 3.
      type(cropfixed_config_t), allocatable :: rotation_fixed(:)
      type(cropgrass_config_t), allocatable :: rotation_grass(:)
      type(cropwofost_config_t), allocatable :: rotation_wofost(:)

      !> .true. for entry i when a per-crop TOML file was successfully
      !! loaded for that rotation slot. Prevents validate from checking
      !! default-constructed sub-configs for entries whose files are absent.
      logical, allocatable :: rotation_loaded(:)
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
         if (allocated(self%rotation_swhydrlift)) then
            call check_int_enum(self%rotation_swhydrlift(i), [0, 1], &
                                "crop.rotation.swhydrlift", errors)
         end if
      end do

      ! Delegate validation to per-entry sub-configs when populated.
      ! Only validate slots that were actually loaded from a file.
      if (allocated(self%rotation_loaded)) then
         do i = 1, n
            if (.not. self%rotation_loaded(i)) cycle
            if (self%rotation_type(i) == 1) then
               if (allocated(self%rotation_fixed)) then
                  call self%rotation_fixed(i)%validate(errors)
               end if
            else if (self%rotation_type(i) == 2) then
               if (allocated(self%rotation_wofost)) then
                  call self%rotation_wofost(i)%validate(errors)
               end if
            else if (self%rotation_type(i) == 3) then
               if (allocated(self%rotation_grass)) then
                  call self%rotation_grass(i)%validate(errors)
               end if
            end if
         end do
      end if
   end subroutine crop_config_validate

   subroutine crop_config_finalize(self, errors)
      class(crop_config_t),     intent(inout) :: self
      type(error_collection_t), intent(inout) :: errors
   end subroutine crop_config_finalize

end module crop_config_mod
