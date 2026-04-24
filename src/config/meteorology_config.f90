!> [meteorology] section config.
module meteorology_config_mod
   use iso_fortran_env, only: real64
   use error_mod, only: error_collection_t
   use validation_mod, only: check_int_enum, check_real_range, check_not_empty
   implicit none
   private

   public :: meteorology_config_t

   type :: meteorology_config_t
      character(len=:), allocatable :: metfile
      real(real64) :: lat  = 0.0_real64
      real(real64) :: alt  = 0.0_real64
      real(real64) :: altw = 2.0_real64
      integer      :: swetr    = 0
      integer      :: swdivide = 0
      integer      :: swmetdetail = 0
      integer      :: nmetdetail  = 0
      integer      :: swrain   = 0
      integer      :: swetsine = 0
      integer      :: swinter  = 0
      integer      :: swmetfilall = 0
      real(real64) :: angstroma = 0.25_real64
      real(real64) :: angstromb = 0.50_real64
   contains
      procedure :: validate => meteorology_config_validate
      procedure :: finalize => meteorology_config_finalize
   end type meteorology_config_t

contains

   subroutine meteorology_config_validate(self, errors)
      class(meteorology_config_t), intent(in)    :: self
      type(error_collection_t),    intent(inout) :: errors

      call check_real_range(self%lat, -90.0_real64, 90.0_real64, "meteorology.lat", errors)
      call check_real_range(self%alt, -500.0_real64, 9000.0_real64, "meteorology.alt", errors)
      call check_int_enum(self%swetr,       [0, 1],    "meteorology.swetr",       errors)
      call check_int_enum(self%swdivide,    [0, 1],    "meteorology.swdivide",    errors)
      call check_int_enum(self%swmetdetail, [0, 1],    "meteorology.swmetdetail", errors)
      call check_int_enum(self%swrain,      [0, 1, 2], "meteorology.swrain",      errors)
      call check_int_enum(self%swinter,     [0, 1, 2], "meteorology.swinter",     errors)
   end subroutine meteorology_config_validate

   subroutine meteorology_config_finalize(self, errors)
      class(meteorology_config_t), intent(inout) :: self
      type(error_collection_t),    intent(inout) :: errors
   end subroutine meteorology_config_finalize

end module meteorology_config_mod
