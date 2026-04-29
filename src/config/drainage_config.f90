!> [drainage] section config.
module drainage_config_mod
   use iso_fortran_env, only: real64
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   use validation_mod, only: check_int_enum, check_int_range, check_real_range, check_nonnegative_real
   implicit none
   private

   public :: drainage_config_t
   public :: drainage_surface_runoff_t

   !> [drainage.surface_runoff] sub-section: surface-runoff and rapid-drainage
   !! controls. Phase 4f-prep extension covering the legacy .dra fields
   !! SWNRSRF, SWTOPNRSRF, SWDIVDINF, FacDpthInf and adjacent rapid-drainage
   !! parameters. All fields are switch-gated; defaults of 0 / 0.0 mean
   !! "section not exercised" and validators skip their bodies.
   type, public :: drainage_surface_runoff_t
      ! Top-level switches
      integer :: swnrsrf      = 0   !! 0=no rapid drainage, 1=yes, 2=extended
      integer :: swtopnrsrf   = 0   !! adjust discharge layer bottom (0,1)
      integer :: swdivdinf    = 0   !! division for infiltration (0,1)
      integer :: swtopdislay  = 0   !! top discharge layer adjustment (0,1)

      ! Real-valued parameters
      real(real64) :: facdpthinf   = 0.0_real64   !! depth-infiltration factor [0..1]
      real(real64) :: cofintfl     = 0.0_real64   !! interflow coefficient
      real(real64) :: expintfl     = 0.0_real64   !! interflow exponent
      real(real64) :: geofac       = 0.0_real64   !! geometry factor (3=isotropic)
      real(real64) :: gwlconv      = 0.0_real64   !! groundwater level convergence
      real(real64) :: ftopdislay   = 0.0_real64   !! fraction top discharge layer [0..1]
      real(real64) :: rsurfdeep    = 0.0_real64   !! surface resistance, deep
      real(real64) :: rsurfshallow = 0.0_real64   !! surface resistance, shallow

      ! Rapid drainage cluster (gated by swnrsrf >= 1, fully required at swnrsrf=2)
      real(real64) :: rapdrareaexp = 0.0_real64
      real(real64) :: rapdraresref = 0.0_real64
      integer      :: numlevrapdra = 0
   contains
      procedure :: validate => drainage_surface_runoff_validate
   end type drainage_surface_runoff_t

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

      type(drainage_surface_runoff_t) :: surface_runoff
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

      call self%surface_runoff%validate(errors)
   end subroutine drainage_config_validate

   subroutine drainage_surface_runoff_validate(self, errors)
      class(drainage_surface_runoff_t), intent(in)    :: self
      type(error_collection_t),         intent(inout) :: errors

      call check_int_enum(self%swnrsrf,     [0, 1, 2], "drainage.surface_runoff.swnrsrf",     errors)
      call check_int_enum(self%swtopnrsrf,  [0, 1],    "drainage.surface_runoff.swtopnrsrf",  errors)
      call check_int_enum(self%swdivdinf,   [0, 1],    "drainage.surface_runoff.swdivdinf",   errors)
      call check_int_enum(self%swtopdislay, [0, 1],    "drainage.surface_runoff.swtopdislay", errors)

      ! Rapid-drainage resistances (required when swnrsrf >= 1)
      if (self%swnrsrf >= 1) then
         call check_real_range(self%rsurfdeep,    1.0e-2_real64, 1.0e6_real64, &
                               "drainage.surface_runoff.rsurfdeep",    errors)
         call check_real_range(self%rsurfshallow, 1.0e-2_real64, 1.0e6_real64, &
                               "drainage.surface_runoff.rsurfshallow", errors)
      end if

      ! Geometry factor for discharge-layer adjustment (3 = isotropic)
      if (self%swtopnrsrf == 1) then
         call check_real_range(self%geofac, 0.0_real64, 100.0_real64, &
                               "drainage.surface_runoff.geofac", errors)
      end if

      ! Top-discharge-layer fraction in [0,1]
      if (self%swtopdislay == 1) then
         call check_real_range(self%ftopdislay, 0.0_real64, 1.0_real64, &
                               "drainage.surface_runoff.ftopdislay", errors)
      end if

      ! Depth-infiltration factor in [0,1]
      if (self%swdivdinf == 1) then
         call check_real_range(self%facdpthinf, 0.0_real64, 1.0_real64, &
                               "drainage.surface_runoff.facdpthinf", errors)
      end if

      ! Always-when-set bounds (skip the default-zero case so that an
      ! unconfigured section validates clean).
      if (self%cofintfl /= 0.0_real64) then
         call check_real_range(self%cofintfl, 0.0_real64, 1.0_real64, &
                               "drainage.surface_runoff.cofintfl", errors)
      end if
      if (self%expintfl /= 0.0_real64) then
         call check_real_range(self%expintfl, 0.0_real64, 10.0_real64, &
                               "drainage.surface_runoff.expintfl", errors)
      end if
      if (self%gwlconv /= 0.0_real64) then
         call check_real_range(self%gwlconv, 0.0_real64, 100.0_real64, &
                               "drainage.surface_runoff.gwlconv", errors)
      end if

      ! Rapid-drainage cluster (fully required at swnrsrf == 2)
      if (self%swnrsrf == 2) then
         call check_real_range(self%rapdrareaexp, 0.0_real64, 100.0_real64, &
                               "drainage.surface_runoff.rapdrareaexp", errors)
         call check_real_range(self%rapdraresref, 0.0_real64, 1.0e6_real64, &
                               "drainage.surface_runoff.rapdraresref", errors)
         call check_int_range(self%numlevrapdra, 1, 5, &
                              "drainage.surface_runoff.numlevrapdra", errors)
      end if
   end subroutine drainage_surface_runoff_validate

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
