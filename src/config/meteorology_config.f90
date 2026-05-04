!> [meteorology] section config.
module meteorology_config_mod
   use iso_fortran_env, only: real64
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   use validation_mod, only: check_int_enum, check_real_range, check_not_empty
   implicit none
   private

   public :: meteorology_config_t
   public :: meteorology_evaporation_t
   public :: meteorology_snow_t

   !> Bare-soil evaporation reduction parameters.
   type :: meteorology_evaporation_t
      integer      :: swcfbs     = 0
      real(real64) :: cfbs       = 1.0_real64
      real(real64) :: cofredbl   = 0.35_real64
      real(real64) :: cofredbo   = 0.35_real64
      real(real64) :: rsigni     = 0.5_real64   !! minimum daily rainfall (cm/d) resetting Black dry counter
      real(real64) :: cfevappond = 1.25_real64  !! ponding-layer evaporation coefficient

      !> Soil-evaporation reduction method (1=Black, 2=Boesten-Stroosnijder).
      !! Legacy SWREDU; default 1 matches all cases except case 5.
      integer :: swredu = 1
   contains
      procedure :: validate => meteorology_evaporation_validate
   end type meteorology_evaporation_t

   !> Snow accumulation / melt parameters.
   type :: meteorology_snow_t
      integer      :: swsnow   = 0
      real(real64) :: snowcoef = 0.0_real64
      real(real64) :: teprrain = 0.0_real64
      real(real64) :: teprsnow = 0.0_real64
   contains
      procedure :: validate => meteorology_snow_validate
   end type meteorology_snow_t

   type :: meteorology_config_t
      character(len=:), allocatable :: metfile
      character(len=:), allocatable :: rainfile
      character(len=:), allocatable :: rain_events_file
      character(len=:), allocatable :: detail_file
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
      type(meteorology_evaporation_t) :: evaporation
      type(meteorology_snow_t)        :: snow
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
      call check_int_enum(self%swrain,      [0, 1, 2, 3], "meteorology.swrain",      errors)
      call check_int_enum(self%swinter,     [0, 1, 2], "meteorology.swinter",     errors)
      call self%evaporation%validate(errors)
      call self%snow%validate(errors)
   end subroutine meteorology_config_validate

   subroutine meteorology_evaporation_validate(self, errors)
      class(meteorology_evaporation_t), intent(in)    :: self
      type(error_collection_t),         intent(inout) :: errors

      call check_int_enum(self%swcfbs, [0, 1], "meteorology.evaporation.swcfbs", errors)
      if (self%swcfbs == 1) then
         call check_real_range(self%cfbs, 0.5_real64, 1.5_real64, &
                               "meteorology.evaporation.cfbs", errors)
      end if
      call check_real_range(self%cofredbl, 0.0_real64, 1.0_real64, &
                            "meteorology.evaporation.cofredbl", errors)
      call check_real_range(self%cofredbo, 0.0_real64, 1.0_real64, &
                            "meteorology.evaporation.cofredbo", errors)
      call check_real_range(self%rsigni, 0.0_real64, 10.0_real64, &
                            "meteorology.evaporation.rsigni", errors)
      call check_real_range(self%cfevappond, 0.0_real64, 10.0_real64, &
                            "meteorology.evaporation.cfevappond", errors)
      call check_int_enum(self%swredu, [1, 2], "meteorology.evaporation.swredu", errors)

      ! Cross-field: swredu=2 (Boesten-Stroosnijder) requires cofredbo to
      ! be set to a non-default value (anything other than the default 0.35).
      if (self%swredu == 2 .and. abs(self%cofredbo - 0.35_real64) < 1.0e-12_real64) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'meteorology.evaporation.swredu=2 requires cofredbo to be ' // &
            'authored (non-default value expected)', &
            'meteorology.evaporation')
      end if
   end subroutine meteorology_evaporation_validate

   subroutine meteorology_snow_validate(self, errors)
      class(meteorology_snow_t), intent(in)    :: self
      type(error_collection_t),  intent(inout) :: errors

      call check_int_enum(self%swsnow, [0, 1], "meteorology.snow.swsnow", errors)
      if (self%swsnow == 1) then
         call check_real_range(self%snowcoef, 0.0_real64, 10.0_real64, &
                               "meteorology.snow.snowcoef", errors)
         call check_real_range(self%teprrain, -10.0_real64, 30.0_real64, &
                               "meteorology.snow.teprrain", errors)
         call check_real_range(self%teprsnow, -10.0_real64, 30.0_real64, &
                               "meteorology.snow.teprsnow", errors)
         if (self%teprsnow > self%teprrain) then
            call errors%append(ERR_VALIDATION_CROSS_FIELD, &
                               "teprsnow > teprrain", "meteorology.snow")
         end if
      end if
   end subroutine meteorology_snow_validate

   subroutine meteorology_config_finalize(self, errors)
      class(meteorology_config_t), intent(inout) :: self
      type(error_collection_t),    intent(inout) :: errors
   end subroutine meteorology_config_finalize

end module meteorology_config_mod
