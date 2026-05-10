!> @file soilwater_state.f90
!! Typed state record for the soil-water boundary subsystem.
!!
!! All 12 included fields are instantaneous (overwritten each step) — flat
!! layout with no cohort sub-records. This matches the heat ADR 0034 precedent.
!!
!! Top-boundary fields (from boundtop):
!!   qtop, reva, hsurf, runots, QMpLatSs, ftoph, FlRunoff
!!
!! Bottom-boundary fields (from BoundBottom):
!!   qbot, qbot_nonfrozen, hbot, gwlinp, deepgw
!!
!! Excluded (config, not runtime state):
!!   - swbotb, swqhbot, swtopb — boundary-condition switches; config flags.
!!   - pond, rsro, rsroexp — surface runoff parameters; config-driven.
!!
!! No type-bound reset() procedure — all fields are instantaneous; the
!! implementer of B-1.3 simply overwrites them each step.
!!
!! See docs/superpowers/specs/2026-05-10-state-migration-boundary-design.md.

module soilwater_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   public :: soilwater_state_t

   type :: soilwater_state_t

      ! Top-boundary fluxes / surface variables (from boundtop)

      real(real64) :: qtop     = 0.0_real64  !! top-surface flux (cm/d)
      real(real64) :: reva     = 0.0_real64  !! actual soil evaporation (cm/d)
      real(real64) :: hsurf    = 0.0_real64  !! pressure head at surface (cm)
      real(real64) :: runots   = 0.0_real64  !! runoff this step (cm)
      real(real64) :: QMpLatSs = 0.0_real64  !! lateral macropore inflow (cm/d)
      logical      :: ftoph    = .false.     !! flag: pressure-head top boundary
      logical      :: FlRunoff = .false.     !! flag: runoff potential

      ! Bottom-boundary fluxes / variables (from BoundBottom)

      real(real64) :: qbot           = 0.0_real64  !! bottom flux (cm/d)
      real(real64) :: qbot_nonfrozen = 0.0_real64  !! bottom flux pre-frost snapshot
      real(real64) :: hbot           = 0.0_real64  !! prescribed head at bottom (cm)
      real(real64) :: gwlinp         = 0.0_real64  !! prescribed gwl, swbotb=1 (cm)
      real(real64) :: deepgw         = 0.0_real64  !! deep-aquifer head, swbotb=3 (cm)

   end type soilwater_state_t

end module soilwater_state_mod
