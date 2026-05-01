!> Aggregate SWAP configuration: composes the six section types.
module swap_config_mod
   use error_mod, only: error_collection_t, ERR_VALIDATION_REQUIRED
   use general_config_mod,      only: general_config_t
   use simulation_config_mod,   only: simulation_config_t
   use meteorology_config_mod,  only: meteorology_config_t
   use drainage_config_mod,     only: drainage_config_t
   use soil_config_mod,         only: soil_config_t
   use bottom_boundary_config_mod, only: bottom_boundary_config_t
   use heat_config_mod,         only: heat_config_t
   use irrigation_config_mod,   only: irrigation_config_t
   use solute_config_mod,       only: solute_config_t
   use surface_water_config_mod, only: surface_water_config_t
   use crop_config_mod,         only: crop_config_t
   implicit none
   private

   public :: swap_config_t

   type :: swap_config_t
      type(general_config_t)     :: general
      type(simulation_config_t)  :: simulation
      type(meteorology_config_t) :: meteo
      type(drainage_config_t)    :: drain
      type(soil_config_t)        :: soil
      type(bottom_boundary_config_t) :: bottom_boundary
      type(heat_config_t)        :: heat
      type(irrigation_config_t)  :: irrigation
      type(solute_config_t)      :: solute
      type(surface_water_config_t) :: surface_water
      type(crop_config_t)        :: crop
   contains
      procedure :: validate => swap_config_validate
      procedure :: finalize => swap_config_finalize
   end type swap_config_t

contains

   subroutine swap_config_validate(self, errors)
      class(swap_config_t),     intent(in)    :: self
      type(error_collection_t), intent(inout) :: errors
      call self%general%validate(errors)
      call self%simulation%validate(errors)
      call self%meteo%validate(errors)
      call self%drain%validate(errors)
      call self%soil%validate(errors)
      call self%bottom_boundary%validate(errors)
      call self%heat%validate(errors)
      call self%irrigation%validate(errors)
      call self%solute%validate(errors)
      call self%surface_water%validate(errors)
      call self%crop%validate(errors)
      ! Cross-section rules are added here as the parity test surfaces them.
      !
      ! Cross-section gating for [soil.initial] required CSV slots.
      ! Range validation has already run inside soil%validate.
      ! Skip when soil.swinco /= 3 OR when legacy inifil path is in use
      ! (matches the transitional contract in soil_config_validate).
      if (self%soil%swinco == 3 .and. &
          (.not. allocated(self%soil%inifil) .or. &
           len_trim_safe(self%soil%inifil) == 0)) then
         if (.not. allocated(self%soil%initial%h_file) .or. &
             len_trim_safe(self%soil%initial%h_file) == 0) then
            call errors%append(ERR_VALIDATION_REQUIRED, &
               'soil.initial.h_file required when soil.swinco=3', &
               'swap_config')
         end if
         if (self%heat%swhea == 1 .and. self%heat%swcalt == 2) then
            if (.not. allocated(self%soil%initial%tsoil_file) .or. &
                len_trim_safe(self%soil%initial%tsoil_file) == 0) then
               call errors%append(ERR_VALIDATION_REQUIRED, &
                  'soil.initial.tsoil_file required when soil.swinco=3 and ' // &
                  'heat.swhea=1 and heat.swcalt=2', 'swap_config')
            end if
         end if
         if (self%solute%swsolu == 1) then
            if (.not. allocated(self%soil%initial%cml_file) .or. &
                len_trim_safe(self%soil%initial%cml_file) == 0) then
               call errors%append(ERR_VALIDATION_REQUIRED, &
                  'soil.initial.cml_file required when soil.swinco=3 and ' // &
                  'solute.swsolu=1', 'swap_config')
            end if
         end if
      end if
   end subroutine swap_config_validate

   !> Allocatable-safe wrapper around len_trim. Returns 0 when not allocated.
   pure function len_trim_safe(s) result(n)
      character(len=:), allocatable, intent(in) :: s
      integer :: n
      if (allocated(s)) then
         n = len_trim(s)
      else
         n = 0
      end if
   end function len_trim_safe

   subroutine swap_config_finalize(self, errors)
      class(swap_config_t),     intent(inout) :: self
      type(error_collection_t), intent(inout) :: errors
      if (errors%has_fatals()) return
      call self%general%finalize(errors)
      call self%simulation%finalize(errors)
      call self%meteo%finalize(errors)
      call self%drain%finalize(errors)
      call self%soil%finalize(errors)
      call self%bottom_boundary%finalize(errors)
      call self%heat%finalize(errors)
      call self%irrigation%finalize(errors)
      call self%solute%finalize(errors)
      call self%surface_water%finalize(errors)
      call self%crop%finalize(errors)
   end subroutine swap_config_finalize

end module swap_config_mod
