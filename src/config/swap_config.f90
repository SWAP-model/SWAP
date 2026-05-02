!> Aggregate SWAP configuration: composes the six section types.
module swap_config_mod
   use iso_fortran_env, only: real64
   use error_mod, only: error_collection_t, ERR_VALIDATION_REQUIRED, &
                        ERR_VALIDATION_CROSS_SECTION
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
      ! Cross-section validation for surface-water management.
      ! drainage_config_validate enforces altcu=0 in the TOML pipeline, so all
      ! coordinates here are in the same (altcu-relative) frame.
      ! zbotdr(1+nrpri) becomes zbotdr(1) since nrpri=0 for swsrf=2 (the only branch
      ! reaching this code, courtesy of upstream stub-errors).
      ! Use minval(zbotdr) to find the deepest channel bottom defensively;
      ! NRSRF tables conventionally place deepest at index 1 but we don't
      ! rely on that ordering here.
      if (self%drain%swdra == 2 .and. self%surface_water%swsrf == 2 .and. &
          self%surface_water%swsec == 2 .and. self%surface_water%swqhr == 1 .and. &
          self%drain%nrlevs >= 1 .and. allocated(self%drain%zbotdr) .and. &
          allocated(self%surface_water%hbweir) .and. &
          allocated(self%surface_water%wldip) .and. &
          allocated(self%surface_water%swman) .and. &
          allocated(self%surface_water%wscap)) then
         block
            integer :: imper
            real(real64) :: zbottom
            zbottom = minval(self%drain%zbotdr(1:self%drain%nrlevs))
            do imper = 1, self%surface_water%nmper
               if (self%surface_water%hbweir(imper) < zbottom) then
                  call errors%append(ERR_VALIDATION_CROSS_SECTION, &
                     'hbweir below deepest channel bottom of secondary system', &
                     'swap_config')
               end if
               if (self%surface_water%swman(imper) == 1 .and. &
                   self%surface_water%wscap(imper) > 1.0e-7_real64 .and. &
                   (self%surface_water%hbweir(imper) - self%surface_water%wldip(imper)) < &
                   (zbottom + 1.0e-4_real64)) then
                  call errors%append(ERR_VALIDATION_CROSS_SECTION, &
                     'target level (hbweir - wldip) below channel bottom; ' // &
                     'supply not possible', 'swap_config')
               end if
            end do
         end block
      end if

      ! Cross-section gating for [soil.initial] required CSV slots.
      ! Range validation has already run inside soil%validate.
      if (self%soil%swinco == 3) then
         block
            logical :: have_h, have_tsoil, have_cml
            have_h = .false.
            if (allocated(self%soil%initial%h_file)) &
               have_h = len_trim(self%soil%initial%h_file) > 0
            if (.not. have_h) then
               call errors%append(ERR_VALIDATION_REQUIRED, &
                  'soil.initial.h_file required when soil.swinco=3', &
                  'swap_config')
            end if
            if (self%heat%swhea == 1 .and. self%heat%swcalt == 2) then
               have_tsoil = .false.
               if (allocated(self%soil%initial%tsoil_file)) &
                  have_tsoil = len_trim(self%soil%initial%tsoil_file) > 0
               if (.not. have_tsoil) then
                  call errors%append(ERR_VALIDATION_REQUIRED, &
                     'soil.initial.tsoil_file required when soil.swinco=3 and ' // &
                     'heat.swhea=1 and heat.swcalt=2', 'swap_config')
               end if
            end if
            if (self%solute%swsolu == 1) then
               have_cml = .false.
               if (allocated(self%soil%initial%cml_file)) &
                  have_cml = len_trim(self%soil%initial%cml_file) > 0
               if (.not. have_cml) then
                  call errors%append(ERR_VALIDATION_REQUIRED, &
                     'soil.initial.cml_file required when soil.swinco=3 and ' // &
                     'solute.swsolu=1', 'swap_config')
               end if
            end if
         end block
      end if
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
      call self%surface_water%finalize(errors)
      call self%crop%finalize(errors)
   end subroutine swap_config_finalize

end module swap_config_mod
