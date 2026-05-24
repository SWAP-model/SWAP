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
   use iso_c_binding, only: c_double
   use swap_array_dimensions, only: MABBC, MAHO
   implicit none
   private
   public :: solute_state_t

   type :: solute_state_t

      ! === per-node arrays (macp-sized; allocated by caller from config) ===
      real(real64), allocatable :: cml(:)    !! soil solute concentration (M/L3 water) in mobile region
      real(real64), allocatable :: cmsy(:)   !! dissolved + adsorbed solute concentration (M/L3 soil volume)

      ! === config-snapshot scalars (snapshotted from config%solute at config_to_variables) ===
      ! Bottom BC / aquifer:
      integer      :: swbotbc = 0           !! bottom-BC type for solute concentration
      integer      :: swbr    = 0           !! mixed-reservoir breakthrough switch
      real(real64) :: daquif  = 0.0_real64  !! aquifer thickness (cm)
      real(real64) :: decsat  = 0.0_real64  !! saturated-zone decay rate (1/d)
      real(real64) :: poros   = 0.0_real64  !! aquifer porosity (-)
      ! Soil chemistry:
      real(real64) :: cref    = 0.0_real64  !! reference concentration for Freundlich adsorption (M/L3)
      real(real64) :: ddif    = 0.0_real64  !! diffusion coefficient (cm2/d)
      real(real64) :: frexp   = 0.0_real64  !! Freundlich exponent (-)
      real(real64) :: kfsat   = 0.0_real64  !! saturated-zone Freundlich coefficient (cm3/g)
      ! Temperature/moisture corrections:
      real(real64) :: gampar  = 0.0_real64  !! temperature decomposition coefficient (/C)
      real(real64) :: bexp    = 0.0_real64  !! moisture-decomposition exponent (-)
      real(real64) :: rtheta  = 0.0_real64  !! reference moisture content (-)
      ! Plant uptake:
      real(real64) :: tscf    = 0.0_real64  !! relative solute uptake by roots (-)
      ! Boundary concentrations:
      real(real64) :: cirr    = 0.0_real64  !! irrigation solute concentration (M/L3)
      real(real64) :: cpre    = 0.0_real64  !! precipitation solute concentration (M/L3)
      ! Initial-condition config (legacy multi-depth init):
      integer      :: nconc   = 0           !! number of initial-concentration depth points

      ! === per-layer config arrays (MAHO-sized) and the seepage time table ===
      real(real64), allocatable :: kf(:)        !! Freundlich coefficient per layer (cm3/g)
      real(real64), allocatable :: decpot(:)    !! potential decomposition rate per layer (1/d)
      real(real64), allocatable :: fdepth(:)    !! depth-decomposition factor per layer (-)
      real(real64), allocatable :: ldis(:)      !! dispersion length per layer (cm)
      real(real64), allocatable :: cseeptab(:)  !! seepage solute concentration table (2*MABBC, time/value pairs)

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

      !> [SS-BMI2] Solute output row buffer (SoluteOutput stream).
      !! Currently placeholder only — SoluteOutput body was deleted by ADR 0009 Phase 5+.
      real(c_double),    allocatable :: output_row(:)
      character(len=32), allocatable :: output_columns(:)
      integer                        :: output_n_cols = 0

      !> [SS-BMI2] AgeTracer output row buffer (AgeTracerOutput stream).
      !! Currently inert — flAgeTracer is always false (ADR 0032).
      real(c_double),    allocatable :: agetracer_row(:)
      character(len=32), allocatable :: agetracer_columns(:)
      integer                        :: agetracer_n_cols = 0

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
