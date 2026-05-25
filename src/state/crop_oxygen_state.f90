!> @file crop_oxygen_state.f90
!! [GR-CROP 2026-05-25] Typed crop runtime state — Bartholomeus oxygen-stress
!! workspace + persistent SAVE state. Hosts the symbols that previously lived
!! as bare globals in `variables.f90` (the `o2_*` cluster) and as `module
!! O2_pars` aliases inside `src/crop/oxygenstress.f90`.
!!
!! Two clusters:
!!
!!  Cluster A — persistent SAVE state (set on first call to OxygenStress in
!!  `calc_ini_pars`, read on every subsequent call within the same run):
!!     ini_stress (logical), d_soil_term1(:), d_soil_term2(:), gfp100(:),
!!     capac_term(:), nmin1(:), mplus1(:)
!!  Allocated by `crop_oxygen_state_init(self, numnod)`; sized to `numnod`.
!!
!!  Cluster B — per-call workspace (previously `module O2_pars` aliases over
!!  `o2_*` bare globals):
!!     w_root, w_root_z0, soil_temp, sat_water_cont, gas_filled_porosity,
!!     d_o2inwater, d_root, d_soil, perc_org_mat, soil_density, depth,
!!     shape_factor_microbialr, root_radius, r_microbial_z0,
!!     waterfilm_thickness, bunsencoeff, c_min_micro, c_macro, ctopnode
!!  Recomputed each call; persistence not required, but kept as type fields
!!  to match the legacy O2_pars sharing pattern between OxygenStress / MICRO /
!!  MACRO / SOLVE / myfunc (the helper functions that ZBREND drives need a
!!  shared scratch slot).
!!
!!  Cluster C — config-derived workspace that OxygenStress and
!!  GET_MAX_RESP_FACTOR share within a call (also recomputed each call from
!!  state%crop%common / cropfixed runtime writes):
!!     max_resp_factor, c_mroot, f_senes, q10_root, q10_microbial,
!!     specific_resp_humus, shape_factor_rootr
!!  Seeded at entry to OxygenStress from the legacy globals (cross-file
!!  writers: cropfixed_init.f90, cropfixed_runtime.f90, cropgrass_init.f90)
!!  until those sub-arcs migrate.
!!
!! The crop_oxygen_state_t is unconditionally allocated on every run; only
!! the `swoxygen=2` Bartholomeus path actually populates it.

module crop_oxygen_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   use swap_array_dimensions, only: macp
   implicit none
   private
   public :: crop_oxygen_state_t

   type :: crop_oxygen_state_t

      ! ---------------------------------------------------------------
      ! Cluster A — persistent SAVE state (per-node, sized to macp to
      ! mirror the legacy variables.f90 fixed-size arrays exactly).
      ! ---------------------------------------------------------------
      logical      :: ini_stress = .true.   !! true until calc_ini_pars has populated the per-node tables
      real(real64) :: d_soil_term1(macp) = 0.0_real64   !! pre-calc soil diffusion term1 per node
      real(real64) :: d_soil_term2(macp) = 0.0_real64   !! pre-calc soil diffusion term2 per node
      real(real64) :: gfp100(macp)       = 0.0_real64   !! gas-filled porosity at h=-100 cm per node
      real(real64) :: capac_term(macp)   = 0.0_real64   !! water-capacity FUNC coefficient per node
      real(real64) :: nmin1(macp)        = 0.0_real64   !! VG n-1 per node
      real(real64) :: mplus1(macp)       = 0.0_real64   !! VG m+1 per node

      ! ---------------------------------------------------------------
      ! Cluster B — per-call workspace (was module O2_pars / o2_*)
      ! ---------------------------------------------------------------
      real(real64) :: w_root                = 0.0_real64 !! dry weight per root length (kg/m)
      real(real64) :: w_root_z0             = 0.0_real64 !! root weight at compartment top (kg/m3)
      real(real64) :: soil_temp             = 0.0_real64 !! soil temperature (K)
      real(real64) :: sat_water_cont        = 0.0_real64 !! saturated water content (-)
      real(real64) :: gas_filled_porosity   = 0.0_real64 !! gas-filled porosity (-)
      real(real64) :: d_o2inwater           = 0.0_real64 !! O2 diffusion in water (m2/d)
      real(real64) :: d_root                = 0.0_real64 !! diffusion in root (m2/d)
      real(real64) :: d_soil                = 0.0_real64 !! soil diffusion (m2/d)
      real(real64) :: perc_org_mat          = 0.0_real64 !! organic matter percentage (%)
      real(real64) :: soil_density          = 0.0_real64 !! soil density (kg/m3)
      real(real64) :: depth                 = 0.0_real64 !! compartment thickness (m)
      real(real64) :: shape_factor_microbialr = 0.0_real64 !! shape factor microbial resp (-)
      real(real64) :: root_radius           = 0.0_real64 !! root radius (m)
      real(real64) :: r_microbial_z0        = 0.0_real64 !! microbial respiration rate at z=0 (kg/m3/d)
      real(real64) :: waterfilm_thickness   = 0.0_real64 !! water film thickness (m)
      real(real64) :: bunsencoeff           = 0.0_real64 !! Bunsen solubility coefficient
      real(real64) :: c_min_micro           = 0.0_real64 !! min O2 for microbial resp (kg/m3)
      real(real64) :: c_macro               = 0.0_real64 !! macropore O2 concentration (kg/m3)
      real(real64) :: ctopnode              = 0.0_real64 !! top-of-compartment O2 concentration (kg/m3)

      ! ---------------------------------------------------------------
      ! Cluster C — config-derived workspace (per-call) — currently
      ! seeded from legacy globals at OxygenStress entry; cross-file
      ! writers in cropfixed_init.f90, cropfixed_runtime.f90,
      ! cropgrass_init.f90 still target the legacy globals until their
      ! sub-arcs migrate.
      ! ---------------------------------------------------------------
      real(real64) :: max_resp_factor     = 1.0_real64
      real(real64) :: c_mroot             = 0.0_real64
      real(real64) :: f_senes             = 0.0_real64
      real(real64) :: q10_root            = 0.0_real64
      real(real64) :: q10_microbial       = 0.0_real64
      real(real64) :: specific_resp_humus = 0.0_real64
      real(real64) :: shape_factor_rootr  = 0.0_real64
      ! [GR-CROP 2026-05-25] w_root_ss: dry weight of roots at soil surface
      ! [kg/m3]. Written by cropfixed_runtime (afgen of wrtb table); read by
      ! OxygenStress in the static-crop (croptype==1) branch.
      real(real64) :: w_root_ss           = 0.0_real64

   contains
      procedure :: init => crop_oxygen_state_init
   end type crop_oxygen_state_t

contains

   subroutine crop_oxygen_state_init(self)
      !! Reset the initialization flag and zero per-node tables. Defaults
      !! on the type definition already cover the zeros; this is a no-op
      !! today but is retained for the lifecycle convention.
      class(crop_oxygen_state_t), intent(inout) :: self
      self%ini_stress   = .true.
      self%d_soil_term1 = 0.0_real64
      self%d_soil_term2 = 0.0_real64
      self%gfp100       = 0.0_real64
      self%capac_term   = 0.0_real64
      self%nmin1        = 0.0_real64
      self%mplus1       = 0.0_real64
   end subroutine crop_oxygen_state_init

end module crop_oxygen_state_mod
