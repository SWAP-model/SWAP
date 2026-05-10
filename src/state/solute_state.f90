!> @file solute_state.f90
!! Typed state record for the solute subsystem. Holds the
!! solute-owned variables that legacy SWAP kept in the
!! `variables.f90` globals module.
!!
!! Excluded:
!!   - AgeTracer-specific globals (12 fields: Ageirr, Agedrain, Agepre,
!!     Agepond, Agepondm1, icAgetopupw, icAgetopdwn, icAgeBot, icAgeDra,
!!     icAgeRot, icAgeSur, AgeGwl1m) — kept in variables.f90 with
!!     provenance comments. AgeTracer is currently inert (flAgeTracer is
!!     always .false.); to reactivate, define agetracer_state_t and
!!     migrate. See discovery doc Hazard #4.
!!   - `ArMpSs` — shared working buffer co-written by solute,
!!     soilhydraulics, and boundtop (all three reset and re-derive it
!!     each timestep). It is not true subsystem state. Stays as a global
!!     until macropore migration sorts ownership. See Hazard #3.
!!
!! `sqrap` and `samcra` are solute balance output fields zeroed in
!! initialize.f90. They are not written by solute.f90 in the current
!! codebase (rapid drainage and crack entrainment paths are inactive), but
!! are included here as they are definitionally part of the solute balance.
!!
!! See ADR 0032-to-be (state-migration solute subsystem) and the
!! 2026-05-10 design spec.

module solute_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   public :: solute_state_t
   public :: solute_intermediate_t
   public :: solute_cumulative_t

   !> Solute intermediate accumulators — reset when flzerointr fires.
   type :: solute_intermediate_t
      real(real64) :: imsqprec  = 0.0_real64
      real(real64) :: imsqirrig = 0.0_real64
      real(real64) :: imsqbot   = 0.0_real64
      real(real64) :: imsqdra   = 0.0_real64
      real(real64) :: imdectot  = 0.0_real64
      real(real64) :: imrottot  = 0.0_real64
   contains
      procedure :: reset => solute_intermediate_reset
   end type solute_intermediate_t

   !> Solute cumulative balance fields — reset when flzerocumu fires.
   !! samini is in the cumulative cohort and zeroed by reset();
   !! the samini = sampro mass-balance rebase is physics not cohort
   !! policy and lives inline at the call site (see solute.f90).
   type :: solute_cumulative_t
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
      procedure :: reset => solute_cumulative_reset
   end type solute_cumulative_t

   type :: solute_state_t

      ! Per-node arrays (macp-sized) — allocated by caller from config dimensions
      real(real64), allocatable :: cml(:)    !! soil solute concentration (M/L3 water) in mobile region
      real(real64), allocatable :: cmsy(:)   !! dissolved + adsorbed solute concentration (M/L3 soil volume)

      ! Scalar state variables updated during solute time-stepping
      real(real64) :: csurf   = 0.0_real64  !! total solutes (M/L2) in ponding layer on soil surface
      real(real64) :: cpond   = 0.0_real64  !! mean solute concentration (M/L3) in ponding layer
      real(real64) :: cdrain  = 0.0_real64  !! mean solute conc in aquifer or drainage system (M/L3 water)
      real(real64) :: cseep   = 0.0_real64  !! mean solute conc in upward seepage water at bottom (M/L3 water)
      real(real64) :: dtsolu  = 0.0_real64  !! solute sub-timestep (T)

      ! Instantaneous fluxes (per-step, reset each outer timestep)
      real(real64) :: isqbot  = 0.0_real64  !! instantaneous solute flux at profile bottom (M/L2/T)
      real(real64) :: isqtop  = 0.0_real64  !! instantaneous solute flux through soil surface (M/L2/T)

      ! Cumulative balance scalars (reset at balance period start)
      real(real64) :: samini  = 0.0_real64  !! total solutes (M/L2) in profile at start of balance period
      real(real64) :: sampro  = 0.0_real64  !! total solutes (M/L2) in soil column (running sum)
      real(real64) :: samcra  = 0.0_real64  !! total solutes (M/L2) entrapped in cracks
      real(real64) :: solbal  = 0.0_real64  !! cumulative solute balance error (M/L2)
      real(real64) :: dectot  = 0.0_real64  !! cumulative solute decomposition (M/L2)
      real(real64) :: imdectot = 0.0_real64 !! intermediate (within output interval) decomposition (M/L2)
      real(real64) :: rottot  = 0.0_real64  !! cumulative solutes extracted by plant roots (M/L2)
      real(real64) :: imrottot = 0.0_real64 !! intermediate root extraction (M/L2)

      ! Cumulative source/sink fluxes
      real(real64) :: sqprec   = 0.0_real64 !! cumulative solutes in precipitation (M/L2)
      real(real64) :: imsqprec = 0.0_real64 !! intermediate solutes in precipitation (M/L2)
      real(real64) :: sqirrig  = 0.0_real64 !! cumulative solutes in irrigation water (M/L2)
      real(real64) :: imsqirrig = 0.0_real64 !! intermediate solutes in irrigation (M/L2)
      real(real64) :: sqbot    = 0.0_real64 !! cumulative solutes through profile bottom (M/L2)
      real(real64) :: imsqbot  = 0.0_real64 !! intermediate solutes through bottom (M/L2)
      real(real64) :: sqdra    = 0.0_real64 !! total solutes transported to drainage canals (M/L2)
      real(real64) :: imsqdra  = 0.0_real64 !! intermediate solutes to drainage (M/L2)
      real(real64) :: sqsur    = 0.0_real64 !! cumulative solutes transported to surface water (M/L2)
      real(real64) :: sqrap    = 0.0_real64 !! cumulative solutes in rapid drainage (M/L2)

   end type solute_state_t

contains

   subroutine solute_intermediate_reset(self)
      class(solute_intermediate_t), intent(inout) :: self
      self%imsqprec  = 0.0_real64
      self%imsqirrig = 0.0_real64
      self%imsqbot   = 0.0_real64
      self%imsqdra   = 0.0_real64
      self%imdectot  = 0.0_real64
      self%imrottot  = 0.0_real64
   end subroutine solute_intermediate_reset

   subroutine solute_cumulative_reset(self)
      class(solute_cumulative_t), intent(inout) :: self
      self%sqprec  = 0.0_real64
      self%sqirrig = 0.0_real64
      self%sqbot   = 0.0_real64
      self%sqdra   = 0.0_real64
      self%sqsur   = 0.0_real64
      self%dectot  = 0.0_real64
      self%rottot  = 0.0_real64
      self%csurf   = 0.0_real64
      self%samini  = 0.0_real64
   end subroutine solute_cumulative_reset

end module solute_state_mod
