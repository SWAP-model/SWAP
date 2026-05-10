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
   public :: surfacewater_cumulative_t

   !> Intermediate accumulators — reset when flzerointr fires.
   type :: surfacewater_intermediate_t
      real(real64) :: iqdra = 0.0_real64    ! intermediate lateral drainage total (cm)
      real(real64), allocatable :: inqdra(:,:)       ! (Madr, macp)
      real(real64), allocatable :: inqdra_in(:,:)    ! (Madr, macp)
      real(real64), allocatable :: inqdra_out(:,:)   ! (Madr, macp)
   contains
      procedure :: reset => surfacewater_intermediate_reset
   end type surfacewater_intermediate_t

   !> Cumulative balance fields — reset when flzerocumu fires.
   type :: surfacewater_cumulative_t
      real(real64) :: cqdrd  = 0.0_real64   ! cumulative drain into reservoir (cm)
      real(real64) :: cwsupp = 0.0_real64   ! cumulative external supply (cm)
      real(real64) :: cwout  = 0.0_real64   ! cumulative outflow (cm)
      real(real64) :: cqdra  = 0.0_real64   ! cumulative lateral drainage, all levels (cm)
      real(real64), allocatable :: cqdrain(:)
      real(real64), allocatable :: cqdrainin(:)
      real(real64), allocatable :: cqdrainout(:)
   contains
      procedure :: reset => surfacewater_cumulative_reset
   end type surfacewater_cumulative_t

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

      ! cohort sub-records — intermediate and cumulative accumulators
      type(surfacewater_intermediate_t) :: intermediate
      type(surfacewater_cumulative_t)   :: cumulative
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

   !> Zero every field in the cumulative cohort. See note above.
   subroutine surfacewater_cumulative_reset(self)
      class(surfacewater_cumulative_t), intent(inout) :: self
      self%cqdrd  = 0.0_real64
      self%cwsupp = 0.0_real64
      self%cwout  = 0.0_real64
      self%cqdra  = 0.0_real64
      if (allocated(self%cqdrain))    self%cqdrain    = 0.0_real64
      if (allocated(self%cqdrainin))  self%cqdrainin  = 0.0_real64
      if (allocated(self%cqdrainout)) self%cqdrainout = 0.0_real64
   end subroutine surfacewater_cumulative_reset

end module surfacewater_state_mod
