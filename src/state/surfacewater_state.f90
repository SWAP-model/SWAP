!> @file surfacewater_state.f90
!! Typed state record for the surface-water subsystem.
!! Excluded fields: `l(Madr)` (drainage config), `fldecdt`
!! (request_smaller_dt argument), `qdra(:,:)` (drainage_state_t).
!! See ADR 0030, ADR 0033, ADR 0042-flatten-reset-cohorts.

module surfacewater_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   public :: surfacewater_state_t

   type :: surfacewater_state_t

      ! === per-step / per-day scalars (no flag-gated reset) ===
      real(real64) :: wls           = 0.0_real64    ! surface water level (cm)
      real(real64) :: wlstar        = 0.0_real64    ! target surface water level (cm)
      real(real64) :: swst          = 0.0_real64    ! storage per unit area (cm)
      real(real64) :: swstini       = 0.0_real64    ! initial storage (cm)
      real(real64) :: hwlman        = 0.0_real64    ! pressure head for target level (cm)
      real(real64) :: vtair         = 0.0_real64    ! total air volume in soil column (cm)
      real(real64) :: wlsold        = 0.0_real64    ! previous-step surface level (cm)
      real(real64) :: ZDraBas       = 0.0_real64    ! drainage basis level (cm)
      real(real64) :: qdrtot        = 0.0_real64    ! total lateral drainage flux (cm/d)

      logical      :: overfl        = .false.       ! automatic weir overflow flag
      logical      :: flInitDraBas  = .true.        ! init drainage basis (macropore)

      integer      :: imper         = 1             ! current management period index
      integer      :: numadj        = 0             ! count of target-level adjustments

      real(real64) :: wlsbak(4)     = 0.0_real64    ! 4-step circular buffer
      real(real64) :: sttab(22, 2)  = 0.0_real64    ! pre-computed level-storage table

      ! === intermediate (reset_intermediate / gate: flzerointr) ===
      real(real64) :: iqdra = 0.0_real64    ! intermediate lateral drainage total (cm)
      real(real64), allocatable :: inqdra(:,:)       ! (Madr, macp)
      real(real64), allocatable :: inqdra_in(:,:)    ! (Madr, macp)
      real(real64), allocatable :: inqdra_out(:,:)   ! (Madr, macp)

      ! === cumulative — drainage subsystem
      !     (reset_cumulative_drainage / gate: flzerocumu + fldrain)
      !     Owner: drainage subsystem. Fields accumulate under fldrain
      !     (swdra=1 OR swdra=2). ===
      real(real64) :: cqdra  = 0.0_real64   ! cumulative lateral drainage (cm)
      real(real64), allocatable :: cqdrain(:)        ! (Madr) cumulative drainage per level
      real(real64), allocatable :: cqdrainin(:)      ! (Madr) cumulative infiltration per level
      real(real64), allocatable :: cqdrainout(:)     ! (Madr) cumulative drainage out per level

      ! === cumulative — reservoir subsystem
      !     (reset_cumulative_reservoir / gate: flzerocumu + flSurfaceWater)
      !     Owner: surface-water subsystem. Fields accumulate only when
      !     flSurfaceWater is true (swdra=2 only). ===
      real(real64) :: cqdrd  = 0.0_real64   ! cumulative drain into reservoir (cm)
      real(real64) :: cwsupp = 0.0_real64   ! cumulative external supply (cm)
      real(real64) :: cwout  = 0.0_real64   ! cumulative outflow (cm)

   contains
      procedure :: reset_intermediate         => surfacewater_reset_intermediate
      procedure :: reset_cumulative_drainage  => surfacewater_reset_cumulative_drainage
      procedure :: reset_cumulative_reservoir => surfacewater_reset_cumulative_reservoir
   end type surfacewater_state_t

contains

   !> Zero the intermediate cohort — flzerointr gate.
   !! Allocatable arrays are zeroed only if allocated.
   subroutine surfacewater_reset_intermediate(self)
      class(surfacewater_state_t), intent(inout) :: self
      self%iqdra = 0.0_real64
      if (allocated(self%inqdra))     self%inqdra     = 0.0_real64
      if (allocated(self%inqdra_in))  self%inqdra_in  = 0.0_real64
      if (allocated(self%inqdra_out)) self%inqdra_out = 0.0_real64
   end subroutine surfacewater_reset_intermediate

   !> Zero the drainage-cumulative cohort — flzerocumu gate, fldrain partition.
   !! Allocatable arrays are zeroed only if allocated.
   subroutine surfacewater_reset_cumulative_drainage(self)
      class(surfacewater_state_t), intent(inout) :: self
      self%cqdra = 0.0_real64
      if (allocated(self%cqdrain))    self%cqdrain    = 0.0_real64
      if (allocated(self%cqdrainin))  self%cqdrainin  = 0.0_real64
      if (allocated(self%cqdrainout)) self%cqdrainout = 0.0_real64
   end subroutine surfacewater_reset_cumulative_drainage

   !> Zero the reservoir-cumulative cohort — flzerocumu gate, flSurfaceWater partition.
   !! Three scalars; no allocatables.
   subroutine surfacewater_reset_cumulative_reservoir(self)
      class(surfacewater_state_t), intent(inout) :: self
      self%cqdrd  = 0.0_real64
      self%cwsupp = 0.0_real64
      self%cwout  = 0.0_real64
   end subroutine surfacewater_reset_cumulative_reservoir

end module surfacewater_state_mod
