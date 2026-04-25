!> Type 3 (grass / WOFOST grass) crop config — populated from a .crp.toml file.
!! Field set is similar to cropfixed_config_t plus grass-specific management
!! (mowing, grazing, fertilizer). Per design Q2A, no shared base type.
module cropgrass_config_mod
   use iso_fortran_env, only: real64
   use error_mod, only: error_collection_t, ERR_VALIDATION_OUT_OF_RANGE
   use validation_mod, only: check_int_enum, check_int_range, check_real_range, &
                             check_nonnegative_real, check_ordered_pair
   implicit none
   private

   public :: cropgrass_config_t

   type :: cropgrass_config_t
      ! Phenology / development
      integer      :: idev = 2                  !! For grass typically 2 (temp-sum-based)
      integer      :: lcc  = 0
      real(real64) :: tbase = 0.0_real64
      real(real64) :: tsum1 = 0.0_real64
      real(real64) :: tsum2 = 0.0_real64

      ! Light & growth
      real(real64) :: kdif = 0.0_real64
      real(real64) :: kdir = 0.0_real64
      real(real64) :: eff  = 0.0_real64
      real(real64) :: amax = 0.0_real64

      ! Tables
      real(real64), allocatable :: cftb(:)
      real(real64), allocatable :: chtb(:)
      real(real64), allocatable :: rdctb(:)

      ! Root growth
      real(real64) :: rdi = 0.0_real64
      real(real64) :: rri = 0.0_real64
      real(real64) :: rdc = 0.0_real64

      ! Water stress (Feddes) — same as fixed
      real(real64) :: hlim1 = 0.0_real64
      real(real64) :: hlim2u = 0.0_real64
      real(real64) :: hlim2l = 0.0_real64
      real(real64) :: hlim3h = 0.0_real64
      real(real64) :: hlim3l = 0.0_real64
      real(real64) :: hlim4 = 0.0_real64
      real(real64) :: adcrh = 0.0_real64
      real(real64) :: adcrl = 0.0_real64
      real(real64) :: rsc   = 0.0_real64

      ! Salinity & interception
      real(real64) :: ecmax  = 0.0_real64
      real(real64) :: ecslop = 0.0_real64
      real(real64) :: cofab  = 0.0_real64

      ! Mowing schedule (grass-specific)
      integer :: swharv = 0                     !! 0=no scheduled mow, 1=scheduled
      integer :: nmow   = 0                     !! Number of mowing events
      real(real64), allocatable :: dates_mowing(:)
      real(real64), allocatable :: lai_after_mow(:)

      ! Grazing (grass-specific)
      integer :: swgraz = 0                     !! 0=no grazing, 1=scheduled
      real(real64) :: nstart_graz = 0.0_real64  !! Start day-of-year
      real(real64) :: nstop_graz  = 0.0_real64  !! Stop day-of-year
   contains
      procedure :: validate => cropgrass_config_validate
      procedure :: finalize => cropgrass_config_finalize
   end type cropgrass_config_t

contains

   subroutine cropgrass_config_validate(self, errors)
      class(cropgrass_config_t), intent(in)    :: self
      type(error_collection_t),  intent(inout) :: errors

      call check_int_enum(self%idev,   [1, 2],    'cropgrass.idev',   errors)
      call check_int_enum(self%swharv, [0, 1],    'cropgrass.swharv', errors)
      call check_int_enum(self%swgraz, [0, 1],    'cropgrass.swgraz', errors)

      call check_real_range(self%rdi, 0.0_real64, 1000.0_real64, 'cropgrass.rdi', errors)
      call check_real_range(self%rdc, 0.0_real64, 1000.0_real64, 'cropgrass.rdc', errors)
      call check_nonnegative_real(self%rri, 'cropgrass.rri', errors)

      call check_ordered_pair(self%hlim2l, self%hlim2u, 'hlim2l', 'hlim2u', 'cropgrass', errors)
      call check_ordered_pair(self%hlim3l, self%hlim3h, 'hlim3l', 'hlim3h', 'cropgrass', errors)

      call check_nonnegative_real(self%ecmax,  'cropgrass.ecmax',  errors)
      call check_nonnegative_real(self%ecslop, 'cropgrass.ecslop', errors)

      if (self%swharv == 1) then
         call check_int_range(self%nmow, 1, 50, 'cropgrass.nmow', errors)
      end if
   end subroutine cropgrass_config_validate

   subroutine cropgrass_config_finalize(self, errors)
      class(cropgrass_config_t), intent(inout) :: self
      type(error_collection_t),  intent(inout) :: errors
   end subroutine cropgrass_config_finalize

end module cropgrass_config_mod
