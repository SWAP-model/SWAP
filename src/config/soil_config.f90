!> [soil] section config.
module soil_config_mod
   use iso_fortran_env, only: real64
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   use validation_mod, only: check_int_enum, check_int_range, &
                             check_real_range, check_nonnegative_real
   implicit none
   private

   public :: soil_config_t
   public :: soil_discretization_t
   public :: soil_frost_t

   !> Optional re-discretization of the vertical grid for output reporting.
   !! When swdiscrvert == 1, dznew(:) (sized to numnodnew) carries the
   !! per-new-node thickness in cm.
   type :: soil_discretization_t
      integer :: swdiscrvert = 0       !! 0=use existing, 1=re-discretize
      integer :: numnodnew   = 0       !! count when swdiscrvert=1
      real(real64), allocatable :: dznew(:)  !! per-new-node thickness (cm)
   contains
      procedure :: validate => soil_discretization_validate
   end type soil_discretization_t

   !> Frost-induced flow reduction parameters.
   type :: soil_frost_t
      integer      :: swfrost   = 0
      real(real64) :: tfroststa = 0.0_real64
      real(real64) :: tfrostend = 0.0_real64
      integer      :: swsublim  = 0
   contains
      procedure :: validate => soil_frost_validate
   end type soil_frost_t

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
      integer :: nrstaring = 0  !! 0=user-supplied, 1..6=Staring series

      integer,      allocatable :: sublay(:)
      real(real64), allocatable :: hcomp(:)
      integer,      allocatable :: ncomp(:)
      integer,      allocatable :: isoillay(:)
      real(real64), allocatable :: orgmat(:)
      real(real64), allocatable :: bdens(:)
      real(real64), allocatable :: wcontent(:)
      real(real64), allocatable :: cofgen(:,:)
      real(real64), allocatable :: cofani(:)  !! per-soil-layer anisotropy ratio

      type(soil_discretization_t) :: discretization
      type(soil_frost_t)          :: frost
   contains
      procedure :: validate => soil_config_validate
      procedure :: finalize => soil_config_finalize
   end type soil_config_t

contains

   subroutine soil_config_validate(self, errors)
      class(soil_config_t),     intent(in)    :: self
      type(error_collection_t), intent(inout) :: errors

      integer :: i

      call check_int_enum(self%swsophy, [0, 1],       "soil.swsophy", errors)
      call check_int_enum(self%swhyst,  [0, 1, 2],    "soil.swhyst",  errors)
      call check_int_enum(self%swinco,  [1, 2, 3],    "soil.swinco",  errors)
      call check_int_enum(self%swmacro, [0, 1],       "soil.swmacro", errors)
      call check_int_enum(self%swscal,  [0, 1],       "soil.swscal",  errors)

      call check_nonnegative_real(self%pondmx,  "soil.pondmx",  errors)
      call check_nonnegative_real(self%ksatexm, "soil.ksatexm", errors)

      call check_int_range(self%nrstaring, 0, 6, "soil.nrstaring", errors)

      if (allocated(self%cofani)) then
         do i = 1, size(self%cofani)
            call check_real_range(self%cofani(i), 0.01_real64, 100.0_real64, &
                                  "soil.cofani", errors)
         end do
      end if

      call self%discretization%validate(errors)
      call self%frost%validate(errors)
   end subroutine soil_config_validate

   subroutine soil_discretization_validate(self, errors)
      class(soil_discretization_t), intent(in)    :: self
      type(error_collection_t),     intent(inout) :: errors

      integer :: i

      call check_int_enum(self%swdiscrvert, [0, 1], &
                          "soil.discretization.swdiscrvert", errors)

      if (self%swdiscrvert == 1) then
         call check_int_range(self%numnodnew, 1, 1000, &
                              "soil.discretization.numnodnew", errors)

         if (.not. allocated(self%dznew)) then
            call errors%append(ERR_VALIDATION_CROSS_FIELD, &
                               "dznew required when swdiscrvert=1", &
                               "soil.discretization")
         else
            if (size(self%dznew) /= self%numnodnew) then
               call errors%append(ERR_VALIDATION_CROSS_FIELD, &
                                  "dznew size != numnodnew", &
                                  "soil.discretization")
            end if
            do i = 1, size(self%dznew)
               call check_real_range(self%dznew(i), 1.0e-6_real64, 1000.0_real64, &
                                     "soil.discretization.dznew", errors)
            end do
         end if
      end if
   end subroutine soil_discretization_validate

   subroutine soil_frost_validate(self, errors)
      class(soil_frost_t),      intent(in)    :: self
      type(error_collection_t), intent(inout) :: errors

      call check_int_enum(self%swfrost,  [0, 1], "soil.frost.swfrost",  errors)
      call check_int_enum(self%swsublim, [0, 1], "soil.frost.swsublim", errors)

      if (self%swfrost == 1) then
         call check_real_range(self%tfroststa, -10.0_real64, 0.0_real64, &
                               "soil.frost.tfroststa", errors)
         call check_real_range(self%tfrostend, -10.0_real64, 0.0_real64, &
                               "soil.frost.tfrostend", errors)
         if (self%tfrostend >= self%tfroststa) then
            call errors%append(ERR_VALIDATION_CROSS_FIELD, &
                               "tfrostend >= tfroststa", "soil.frost")
         end if
      end if
   end subroutine soil_frost_validate

   subroutine soil_config_finalize(self, errors)
      class(soil_config_t),     intent(inout) :: self
      type(error_collection_t), intent(inout) :: errors
   end subroutine soil_config_finalize

end module soil_config_mod
