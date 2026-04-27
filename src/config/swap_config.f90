!> Aggregate SWAP configuration: composes the six section types.
module swap_config_mod
   use error_mod, only: error_collection_t
   use general_config_mod,      only: general_config_t
   use simulation_config_mod,   only: simulation_config_t
   use meteorology_config_mod,  only: meteorology_config_t
   use drainage_config_mod,     only: drainage_config_t
   use soil_config_mod,         only: soil_config_t
   use bottom_boundary_config_mod, only: bottom_boundary_config_t
   use heat_config_mod,         only: heat_config_t
   use irrigation_config_mod,   only: irrigation_config_t
   use solute_config_mod,       only: solute_config_t
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
      call self%crop%validate(errors)
      ! Cross-section rules are added here as the parity test surfaces them.
   end subroutine swap_config_validate

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
      call self%crop%finalize(errors)
   end subroutine swap_config_finalize

end module swap_config_mod
