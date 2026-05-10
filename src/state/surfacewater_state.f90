!> @file surfacewater_state.f90
!! Typed state record for the surface-water subsystem. Holds the
!! surface-water-owned variables that legacy SWAP kept in the
!! `variables.f90` globals module. Excluded:
!!   - `l(Madr)` — drainage-config (m→cm conversion moves to drainage
!!     config-load in Phase 2; not surface-water state).
!!   - `fldecdt` — replaced by `intent(out) :: request_smaller_dt`
!!     argument on `SurfaceWater(task=2)`; not a state field.
!!   - `qdra(:,:)` — moved to drainage_state_t in the drainage
!!     migration arc (ADR 0031 / SS-DRST Task 3).
!! See ADR 0030 (state-migration pilot) and the 2026-05-09 design spec
!! for rationale and field provenance.

module surfacewater_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   public :: surfacewater_state_t
   public :: surfacewater_intermediate_t
   public :: surfacewater_drainage_cumulative_t
   public :: surfacewater_reservoir_cumulative_t

   !> Intermediate accumulators — reset when flzerointr fires.
   !! Owner: drainage subsystem (gated by fldrain).
   type :: surfacewater_intermediate_t
      real(real64) :: iqdra = 0.0_real64    ! intermediate lateral drainage total (cm)
      real(real64), allocatable :: inqdra(:,:)       ! (Madr, macp)
      real(real64), allocatable :: inqdra_in(:,:)    ! (Madr, macp)
      real(real64), allocatable :: inqdra_out(:,:)   ! (Madr, macp)
   contains
      procedure :: reset => surfacewater_intermediate_reset
   end type surfacewater_intermediate_t

   !> Drainage-cumulative cohort — accumulates only when fldrain is true
   !! (i.e., swdra=1 OR swdra=2). Owner: drainage subsystem; reset() is
   !! called from Drainage(2) and from SurfaceWater(2) (which only runs
   !! under swdra=2 — but when it does, it's the canonical reset site
   !! because it runs after Drainage in the timestep). Both call sites
   !! reset the full cohort: subset-zeroing is unnecessary because all
   !! 4 fields accumulate under the same fldrain gate.
   type :: surfacewater_drainage_cumulative_t
      real(real64) :: cqdra  = 0.0_real64   ! cumulative lateral drainage, all levels (cm)
      real(real64), allocatable :: cqdrain(:)        ! (Madr) cumulative drainage per level
      real(real64), allocatable :: cqdrainin(:)      ! (Madr) cumulative infiltration per level
      real(real64), allocatable :: cqdrainout(:)     ! (Madr) cumulative drainage out per level
   contains
      procedure :: reset => surfacewater_drainage_cumulative_reset
   end type surfacewater_drainage_cumulative_t

   !> Reservoir-cumulative cohort — accumulates only when flSurfaceWater
   !! is true (i.e., swdra=2 with reservoir simulation). Owner:
   !! surface-water subsystem; reset() is called only from SurfaceWater(2).
   !! Under swdra=1 these fields never accumulate, so they need no reset
   !! site in Drainage(2).
   type :: surfacewater_reservoir_cumulative_t
      real(real64) :: cqdrd  = 0.0_real64   ! cumulative drain into reservoir (cm)
      real(real64) :: cwsupp = 0.0_real64   ! cumulative external supply (cm)
      real(real64) :: cwout  = 0.0_real64   ! cumulative outflow (cm)
   contains
      procedure :: reset => surfacewater_reservoir_cumulative_reset
   end type surfacewater_reservoir_cumulative_t

   type :: surfacewater_state_t
      ! per-step / per-day scalars
      real(real64) :: wls           = 0.0_real64    ! surface water level (cm)
      real(real64) :: wlstar        = 0.0_real64    ! target surface water level (cm)
      real(real64) :: swst          = 0.0_real64    ! storage per unit area (cm)
      real(real64) :: swstini       = 0.0_real64    ! initial storage (cm)
      real(real64) :: hwlman        = 0.0_real64    ! pressure head for target level (cm)
      real(real64) :: vtair         = 0.0_real64    ! total air volume in soil column (cm)
      real(real64) :: wlsold        = 0.0_real64    ! previous-step surface level (cm)
      real(real64) :: ZDraBas       = 0.0_real64    ! drainage basis level for rapid macropore drainage (cm)
      real(real64) :: qdrtot        = 0.0_real64    ! total lateral drainage flux (cm/d)

      logical      :: overfl        = .false.       ! automatic weir overflow flag
      logical      :: flInitDraBas  = .true.        ! init drainage basis (macropore)

      integer      :: imper         = 1             ! current management period index
      integer      :: numadj        = 0             ! count of target-level adjustments

      real(real64) :: wlsbak(4)     = 0.0_real64    ! 4-step circular buffer for oscillation detection
      real(real64) :: sttab(22, 2)  = 0.0_real64    ! pre-computed level-storage table

      ! cohort sub-records — intermediate and partitioned cumulative accumulators
      type(surfacewater_intermediate_t)          :: intermediate
      type(surfacewater_drainage_cumulative_t)   :: drainage_cumulative   ! gate: fldrain
      type(surfacewater_reservoir_cumulative_t)  :: reservoir_cumulative  ! gate: flSurfaceWater
   end type surfacewater_state_t

contains

   !> Zero every field in the intermediate cohort. Allocatable arrays are
   !! zeroed only if allocated; allocation lifecycle stays with the caller.
   subroutine surfacewater_intermediate_reset(self)
      class(surfacewater_intermediate_t), intent(inout) :: self
      self%iqdra = 0.0_real64
      if (allocated(self%inqdra))     self%inqdra     = 0.0_real64
      if (allocated(self%inqdra_in))  self%inqdra_in  = 0.0_real64
      if (allocated(self%inqdra_out)) self%inqdra_out = 0.0_real64
   end subroutine surfacewater_intermediate_reset

   !> Zero the drainage-cumulative cohort. Allocatable arrays are zeroed
   !! only if allocated; allocation lifecycle stays with the caller.
   subroutine surfacewater_drainage_cumulative_reset(self)
      class(surfacewater_drainage_cumulative_t), intent(inout) :: self
      self%cqdra = 0.0_real64
      if (allocated(self%cqdrain))    self%cqdrain    = 0.0_real64
      if (allocated(self%cqdrainin))  self%cqdrainin  = 0.0_real64
      if (allocated(self%cqdrainout)) self%cqdrainout = 0.0_real64
   end subroutine surfacewater_drainage_cumulative_reset

   !> Zero the reservoir-cumulative cohort. Three scalars; no allocatables.
   subroutine surfacewater_reservoir_cumulative_reset(self)
      class(surfacewater_reservoir_cumulative_t), intent(inout) :: self
      self%cqdrd  = 0.0_real64
      self%cwsupp = 0.0_real64
      self%cwout  = 0.0_real64
   end subroutine surfacewater_reservoir_cumulative_reset

end module surfacewater_state_mod
