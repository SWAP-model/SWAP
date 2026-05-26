!> @file heat_state.f90
!! Typed state record for the heat (soil-temperature) subsystem.
!!
!! All 12 included fields are instantaneous (no flzero* gating) — flat
!! layout with no cohort sub-records. Heat is the first state migration
!! without cohort sub-records (per ADR 0034).
!!
!! Excluded:
!!   - flTemperature — config flag owned by timecontrol/core, not heat
!!     home tree. Analogous to flSolute; gates all heat compute but is
!!     set at run-init from swhea, not by temperature.f90.
!!   - swhea, swcalt, swbotbhea, swtopbhea — config flags, not state.
!!   - tfroststa, tfrostend — frost-reduction parameters; config-driven.
!!   - nheat, zh — initial-condition config; not runtime state.
!!   - ddamp, tampli, tmean, timref — Phase 0 config fields (analytical
!!     method); not state.
!!   - temtoptab, tembtab — Phase 0 config tables.
!!
!! Note: rfcp is seeded to 1.0 (legacy FrozenCond convention) by the
!! caller at init time (heat_init, Task 3). The type carries no default
!! initializer for it because it is allocatable.
!!
!! See ADR 0034 and docs/superpowers/specs/2026-05-10-state-migration-heat-design.md.

module heat_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   use heat_config_mod, only: heat_config_t
   implicit none
   private
   public :: heat_state_t

   type :: heat_state_t

      ! Per-node arrays (macp-sized) — allocated by caller (heat_init / Temperature task=1)

      real(real64), allocatable :: tsoil(:)    !! soil temperature (°C) per compartment
      real(real64), allocatable :: heacap(:)   !! heat capacity (J/cm³/K) per compartment (numerical only)
      real(real64), allocatable :: heacon(:)   !! heat conductivity (J/cm/K/d) per compartment
      real(real64), allocatable :: rfcp(:)     !! reduction factor for frozen conditions per node (0–1); seeded 1.0
      real(real64), allocatable :: fquartz(:)  !! gravimetric sand+silt fraction per node (g/g mineral)
      real(real64), allocatable :: fclay(:)    !! gravimetric clay fraction per node (g/g mineral)
      real(real64), allocatable :: forg(:)     !! gravimetric organic matter fraction per node (g/g mineral)

      ! Scalar state variables — computed each timestep (numerical mode)

      real(real64) :: tetop      = 0.0_real64  !! temperature (°C) at top of soil profile (under snow cover)
      real(real64) :: tebot      = 0.0_real64  !! temperature (°C) at bottom of soil profile
      real(real64) :: zfrostbot  = 0.0_real64  !! depth of bottom of frozen layer (L)
      real(real64) :: zfrosttop  = 0.0_real64  !! depth of top of frozen layer (L)
      integer      :: nodfrostbot = 0           !! node number of deepest frozen node

   contains
      procedure :: init => heat_state_init
   end type heat_state_t

contains

   subroutine heat_state_init(self, heat_cfg, numnod_in)
      use, intrinsic :: iso_fortran_env, only: real64
      use heat_config_mod, only: heat_config_t
      class(heat_state_t),    intent(inout) :: self
      type(heat_config_t),    intent(in)    :: heat_cfg
      integer,                intent(in)    :: numnod_in

      ! heat_cfg currently unused; signature reserved for future seed migration.
      if (.not. allocated(self%tsoil))   allocate(self%tsoil(numnod_in))
      if (.not. allocated(self%heacap))  allocate(self%heacap(numnod_in))
      if (.not. allocated(self%heacon))  allocate(self%heacon(numnod_in))
      if (.not. allocated(self%rfcp))    allocate(self%rfcp(numnod_in))
      if (.not. allocated(self%fquartz)) allocate(self%fquartz(numnod_in))
      if (.not. allocated(self%fclay))   allocate(self%fclay(numnod_in))
      if (.not. allocated(self%forg))    allocate(self%forg(numnod_in))

      self%tsoil   = 0.0_real64
      self%heacap  = 0.0_real64
      self%heacon  = 0.0_real64
      self%rfcp    = 1.0_real64    ! NB: 1.0 not 0.0 — matches legacy initial value
      self%fquartz = 0.0_real64
      self%fclay   = 0.0_real64
      self%forg    = 0.0_real64
   end subroutine heat_state_init

end module heat_state_mod
