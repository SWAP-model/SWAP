!> Type 1 (fixed/simple) crop config — populated from a .crp.toml file.
!! Field set is the COMMON SUBSET shared with type 3 (grass), plus type 1
!! specific fields. Phase 4c-a starts with a core schema; additional fields
!! are added during parity-test iteration as needed.
module cropfixed_config_mod
   use iso_fortran_env, only: real64
   use error_mod, only: error_collection_t, ERR_VALIDATION_OUT_OF_RANGE
   use validation_mod, only: check_int_enum, check_int_range, check_real_range, &
                             check_nonnegative_real, check_ordered_pair
   implicit none
   private

   public :: cropfixed_config_t

   type :: cropfixed_config_t
      ! Phenology
      integer      :: idev = 1                  !! 1=fixed period, 2=temperature-sum-based
      integer      :: lcc  = 0                  !! Length of crop cycle (days), used when idev=1

      ! Light & growth (when idev=2 path is taken; safe defaults for idev=1)
      real(real64) :: kdif = 0.0_real64         !! Diffuse light extinction
      real(real64) :: kdir = 0.0_real64         !! Direct light extinction

      ! Crop factor & height tables (stored as flat real arrays;
      ! pairs of (development_stage, value)).
      real(real64), allocatable :: cftb(:)      !! Crop factor table
      real(real64), allocatable :: chtb(:)      !! Crop height table

      ! Root growth
      real(real64) :: rdi = 0.0_real64          !! Initial rooting depth (cm)
      real(real64) :: rri = 0.0_real64          !! Daily root extension rate (cm/d)
      real(real64) :: rdc = 0.0_real64          !! Maximum rooting depth (cm)
      real(real64), allocatable :: rdctb(:)     !! Root density distribution table

      ! Water stress (Feddes)
      real(real64) :: hlim1 = 0.0_real64
      real(real64) :: hlim2u = 0.0_real64
      real(real64) :: hlim2l = 0.0_real64
      real(real64) :: hlim3h = 0.0_real64
      real(real64) :: hlim3l = 0.0_real64
      real(real64) :: hlim4 = 0.0_real64
      real(real64) :: adcrh = 0.0_real64
      real(real64) :: adcrl = 0.0_real64
      real(real64) :: rsc   = 0.0_real64        !! Crop resistance for ET method

      ! Salinity stress
      real(real64) :: ecmax  = 0.0_real64
      real(real64) :: ecslop = 0.0_real64

      ! Interception
      real(real64) :: cofab = 0.0_real64
   contains
      procedure :: validate => cropfixed_config_validate
      procedure :: finalize => cropfixed_config_finalize
   end type cropfixed_config_t

contains

   subroutine cropfixed_config_validate(self, errors)
      class(cropfixed_config_t), intent(in)    :: self
      type(error_collection_t),  intent(inout) :: errors

      call check_int_enum(self%idev, [1, 2], 'cropfixed.idev', errors)
      if (self%idev == 1) then
         call check_int_range(self%lcc, 1, 366, 'cropfixed.lcc', errors)
      end if

      call check_real_range(self%rdi, 0.0_real64, 1000.0_real64, 'cropfixed.rdi', errors)
      call check_real_range(self%rdc, 0.0_real64, 1000.0_real64, 'cropfixed.rdc', errors)
      call check_nonnegative_real(self%rri, 'cropfixed.rri', errors)

      ! Feddes: hlim1 (saturation, near 0) > hlim2u > hlim2l > hlim3h > hlim3l > hlim4 (wilting)
      ! Permit 0 (the validator below catches "all zero" defaults via cross-field rule).
      call check_ordered_pair(self%hlim2l, self%hlim2u, 'hlim2l', 'hlim2u', 'cropfixed', errors)
      call check_ordered_pair(self%hlim3l, self%hlim3h, 'hlim3l', 'hlim3h', 'cropfixed', errors)

      call check_nonnegative_real(self%ecmax,  'cropfixed.ecmax',  errors)
      call check_nonnegative_real(self%ecslop, 'cropfixed.ecslop', errors)
   end subroutine cropfixed_config_validate

   subroutine cropfixed_config_finalize(self, errors)
      class(cropfixed_config_t), intent(inout) :: self
      type(error_collection_t),  intent(inout) :: errors
      ! No derived fields in 4c-a.
   end subroutine cropfixed_config_finalize

end module cropfixed_config_mod
