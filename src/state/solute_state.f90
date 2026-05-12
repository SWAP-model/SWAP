!> @file solute_state.f90
!! Typed state record for the solute subsystem.
!!
!! Excluded:
!!   - AgeTracer-specific globals (12 fields) — kept in variables.f90.
!!   - `ArMpSs` — shared working buffer; stays as a global until macropore
!!     migration sorts ownership.
!!
!! `sqrap` and `samcra` are solute balance output fields zeroed in
!! initialize.f90; included here as part of the solute balance.
!!
!! Reset cadence is expressed by named procedures on the parent type:
!!   - reset_intermediate() — flzerointr gate (6 fields)
!!   - reset_cumulative()   — flzerocumu gate (9 fields)
!! The `samini = sampro` rebase is physics, not cohort policy; it stays
!! inline at the call site in solute.f90.
!!
!! Originally introduced as nested cohort sub-records in ADR 0033 (Phase B).
!! Flattened in the 2026-05-12 reset-cohort-flattening arc.

module solute_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   public :: solute_state_t

   type :: solute_state_t

      ! === per-node arrays (macp-sized; allocated by caller from config) ===
      real(real64), allocatable :: cml(:)    !! soil solute concentration (M/L3 water) in mobile region
      real(real64), allocatable :: cmsy(:)   !! dissolved + adsorbed solute concentration (M/L3 soil volume)

      ! === scalar state updated during solute time-stepping (no flag-gated reset) ===
      real(real64) :: cpond   = 0.0_real64
      real(real64) :: cdrain  = 0.0_real64
      real(real64) :: cseep   = 0.0_real64
      real(real64) :: dtsolu  = 0.0_real64

      ! === instantaneous fluxes (per-step; zeroed unconditionally inside solute(2)) ===
      real(real64) :: isqbot  = 0.0_real64
      real(real64) :: isqtop  = 0.0_real64

      ! === running totals (not reset by flzerointr/flzerocumu) ===
      real(real64) :: sampro  = 0.0_real64
      real(real64) :: samcra  = 0.0_real64
      real(real64) :: solbal  = 0.0_real64
      real(real64) :: sqrap   = 0.0_real64

      ! === intermediate (reset_intermediate / gate: flzerointr) ===
      real(real64) :: imsqprec  = 0.0_real64
      real(real64) :: imsqirrig = 0.0_real64
      real(real64) :: imsqbot   = 0.0_real64
      real(real64) :: imsqdra   = 0.0_real64
      real(real64) :: imdectot  = 0.0_real64
      real(real64) :: imrottot  = 0.0_real64

      ! === cumulative (reset_cumulative / gate: flzerocumu) ===
      !! samini is in the cumulative cohort and zeroed by reset_cumulative();
      !! the samini = sampro mass-balance rebase is physics not cohort
      !! policy and lives inline at the call site (see solute.f90).
      real(real64) :: sqprec  = 0.0_real64
      real(real64) :: sqirrig = 0.0_real64
      real(real64) :: sqbot   = 0.0_real64
      real(real64) :: sqdra   = 0.0_real64
      real(real64) :: sqsur   = 0.0_real64
      real(real64) :: dectot  = 0.0_real64
      real(real64) :: rottot  = 0.0_real64
      real(real64) :: csurf   = 0.0_real64
      real(real64) :: samini  = 0.0_real64

   contains
      procedure :: reset_intermediate => solute_reset_intermediate
      procedure :: reset_cumulative   => solute_reset_cumulative
   end type solute_state_t

contains

   !> Zero the 6 intermediate fields. Called under flzerointr.
   subroutine solute_reset_intermediate(self)
      class(solute_state_t), intent(inout) :: self
      self%imsqprec  = 0.0_real64
      self%imsqirrig = 0.0_real64
      self%imsqbot   = 0.0_real64
      self%imsqdra   = 0.0_real64
      self%imdectot  = 0.0_real64
      self%imrottot  = 0.0_real64
   end subroutine solute_reset_intermediate

   !> Zero the 9 cumulative fields. Called under flzerocumu.
   subroutine solute_reset_cumulative(self)
      class(solute_state_t), intent(inout) :: self
      self%sqprec  = 0.0_real64
      self%sqirrig = 0.0_real64
      self%sqbot   = 0.0_real64
      self%sqdra   = 0.0_real64
      self%sqsur   = 0.0_real64
      self%dectot  = 0.0_real64
      self%rottot  = 0.0_real64
      self%csurf   = 0.0_real64
      self%samini  = 0.0_real64
   end subroutine solute_reset_cumulative

end module solute_state_mod
