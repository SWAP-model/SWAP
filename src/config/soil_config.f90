!> [soil] section config.
module soil_config_mod
   use iso_fortran_env, only: real64
   use error_mod, only: error_collection_t
   use validation_mod, only: check_int_enum, check_real_range, check_nonnegative_real
   implicit none
   private

   public :: soil_config_t

   type :: soil_config_t
      integer :: swsophy = 0
      integer :: swhyst  = 0
      integer :: swinco  = 1
      integer :: swmacro = 0
      integer :: swscal  = 0

      real(real64) :: gwli    = 0.0_real64
      real(real64) :: pondini = 0.0_real64
      real(real64) :: pondmx  = 0.0_real64
      real(real64) :: ksatexm = 0.0_real64
      real(real64) :: rsoil   = 0.0_real64

      integer :: reva_top = 0

      integer,      allocatable :: sublay(:)
      real(real64), allocatable :: hcomp(:)
      integer,      allocatable :: ncomp(:)
      integer,      allocatable :: isoillay(:)
      real(real64), allocatable :: orgmat(:)
      real(real64), allocatable :: bdens(:)
      real(real64), allocatable :: wcontent(:)
      real(real64), allocatable :: cofgen(:,:)
   contains
      procedure :: validate => soil_config_validate
      procedure :: finalize => soil_config_finalize
   end type soil_config_t

contains

   subroutine soil_config_validate(self, errors)
      class(soil_config_t),     intent(in)    :: self
      type(error_collection_t), intent(inout) :: errors

      call check_int_enum(self%swsophy, [0, 1],       "soil.swsophy", errors)
      call check_int_enum(self%swhyst,  [0, 1, 2],    "soil.swhyst",  errors)
      call check_int_enum(self%swinco,  [1, 2, 3],    "soil.swinco",  errors)
      call check_int_enum(self%swmacro, [0, 1],       "soil.swmacro", errors)
      call check_int_enum(self%swscal,  [0, 1],       "soil.swscal",  errors)

      call check_nonnegative_real(self%pondmx,  "soil.pondmx",  errors)
      call check_nonnegative_real(self%ksatexm, "soil.ksatexm", errors)
   end subroutine soil_config_validate

   subroutine soil_config_finalize(self, errors)
      class(soil_config_t),     intent(inout) :: self
      type(error_collection_t), intent(inout) :: errors
   end subroutine soil_config_finalize

end module soil_config_mod
