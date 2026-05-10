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

   type :: surfacewater_state_t
      ! per-step / per-day scalars
      real(real64) :: wls           = 0.0_real64    ! surface water level (cm)
      real(real64) :: wlstar        = 0.0_real64    ! target surface water level (cm)
      real(real64) :: swst          = 0.0_real64    ! storage per unit area (cm)
      real(real64) :: swstini       = 0.0_real64    ! initial storage (cm)
      real(real64) :: hwlman        = 0.0_real64    ! pressure head for target level (cm)
      real(real64) :: vtair         = 0.0_real64    ! total air volume in soil column (cm)
      real(real64) :: wlsold        = 0.0_real64    ! previous-step surface level (cm)
      real(real64) :: cqdrd         = 0.0_real64    ! cumulative drain into reservoir (cm)
      real(real64) :: cwsupp        = 0.0_real64    ! cumulative external supply (cm)
      real(real64) :: cwout         = 0.0_real64    ! cumulative outflow (cm)
      real(real64) :: cqdra         = 0.0_real64    ! cumulative lateral drainage, all levels (cm)
      real(real64) :: ZDraBas       = 0.0_real64    ! drainage basis level for rapid macropore drainage (cm)
      real(real64) :: iqdra         = 0.0_real64    ! intermediate lateral drainage total (cm)
      real(real64) :: qdrtot        = 0.0_real64    ! total lateral drainage flux (cm/d)

      logical      :: overfl        = .false.       ! automatic weir overflow flag
      logical      :: flInitDraBas  = .true.        ! init drainage basis (macropore)

      integer      :: imper         = 1             ! current management period index
      integer      :: numadj        = 0             ! count of target-level adjustments

      real(real64) :: wlsbak(4)     = 0.0_real64    ! 4-step circular buffer for oscillation detection
      real(real64) :: sttab(22, 2)  = 0.0_real64    ! pre-computed level-storage table

      ! per-level (Madr-sized) arrays — allocated by surfacewater_init
      real(real64), allocatable :: cqdrain(:)
      real(real64), allocatable :: cqdrainin(:)
      real(real64), allocatable :: cqdrainout(:)
      real(real64), allocatable :: inqdra(:,:)       ! (Madr, macp)
      real(real64), allocatable :: inqdra_in(:,:)    ! (Madr, macp)
      real(real64), allocatable :: inqdra_out(:,:)   ! (Madr, macp)
   end type surfacewater_state_t

end module surfacewater_state_mod
