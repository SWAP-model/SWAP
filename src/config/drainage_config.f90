!> [drainage] section config.
module drainage_config_mod
   use iso_fortran_env, only: real64
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   use validation_mod, only: check_int_enum, check_int_range, check_real_range, check_nonnegative_real
   implicit none
   private

   public :: drainage_config_t

   type :: drainage_config_t
      integer :: swdra    = 0
      integer :: dramet   = 0
      integer :: swdivd   = 0
      integer :: swdislay = 0
      integer :: nrlevs   = 0

      real(real64) :: altcu  = 0.0_real64
      real(real64) :: basegw = 0.0_real64
      real(real64) :: entres = 0.0_real64
      real(real64) :: shape  = 0.0_real64

      integer,      allocatable :: swdtyp(:)
      real(real64), allocatable :: zbotdr(:)
      real(real64), allocatable :: drares(:)
      real(real64), allocatable :: infres(:)
      real(real64), allocatable :: L(:)
      real(real64), allocatable :: gwlinf(:)
      real(real64), allocatable :: rdrain(:)
      real(real64), allocatable :: rinfi(:)
      real(real64), allocatable :: rentry(:)
      real(real64), allocatable :: rexit(:)
      real(real64), allocatable :: widthr(:)
      real(real64), allocatable :: taludr(:)
      integer,      allocatable :: swallo(:)
   contains
      procedure :: validate => drainage_config_validate
      procedure :: finalize => drainage_config_finalize
   end type drainage_config_t

contains

   subroutine drainage_config_validate(self, errors)
      class(drainage_config_t), intent(in)    :: self
      type(error_collection_t), intent(inout) :: errors

      call check_int_enum(self%swdra,    [0, 1, 2],      "drainage.swdra",    errors)
      call check_int_enum(self%dramet,   [0, 1, 2, 3],   "drainage.dramet",   errors)
      call check_int_enum(self%swdivd,   [0, 1],         "drainage.swdivd",   errors)
      call check_int_enum(self%swdislay, [0, 1],         "drainage.swdislay", errors)
      call check_int_range(self%nrlevs,  0, 5,           "drainage.nrlevs",   errors)

      if (self%dramet == 2 .and. self%swdivd /= 1) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
                            "swdivd must be 1 when dramet=2", &
                            "drainage")
      end if
   end subroutine drainage_config_validate

   subroutine drainage_config_finalize(self, errors)
      class(drainage_config_t), intent(inout) :: self
      type(error_collection_t), intent(inout) :: errors

      ! Mirror legacy convention from readswap.f90: single-level drainage
      ! methods (dramet 1 or 2) clobber nrlevs to 1 regardless of input.
      ! Matching this is required for parity with the legacy reader. The
      ! validator already constrains nrlevs in [0, 5]; this finalize step
      ! lands AFTER validate.
      if (self%dramet /= 3) self%nrlevs = 1
   end subroutine drainage_config_finalize

end module drainage_config_mod
