!> @file surfacewater_state.f90
!! Typed state record for the surface-water subsystem.
!! Excluded fields: `l(Madr)` (drainage config), `fldecdt`
!! (request_smaller_dt argument), `qdra(:,:)` (drainage_state_t).
!! See ADR 0030, ADR 0033, ADR 0042-flatten-reset-cohorts.

module surfacewater_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   use surface_water_config_mod, only: surface_water_config_t
   use drainage_config_mod,      only: drainage_config_t
   use error_mod,                only: fatalerr_collected
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
      procedure :: init                       => surfacewater_state_init
      procedure :: reset_intermediate         => surfacewater_reset_intermediate
      procedure :: reset_cumulative_drainage  => surfacewater_reset_cumulative_drainage
      procedure :: reset_cumulative_reservoir => surfacewater_reset_cumulative_reservoir
   end type surfacewater_state_t

contains

   !> One-time runtime initialization for surfacewater state.
   !! Replaces the surviving math from legacy `rddre` (now retired from
   !! readswap) plus the post-init block inside `SurfaceWater(task=1)`.
   !!
   !! Scope: swsrf=2, swsec=2, swqhr=1, swman=1, drainage.altcu=0 only.
   !! Other branches are guarded with fatalerr_collected (defense in
   !! depth — surface_water_config_validate rejects them upstream too).
   subroutine surfacewater_state_init(self, config_sw, config_drain, numnod)
      class(surfacewater_state_t),  intent(inout) :: self
      type(surface_water_config_t), intent(in)    :: config_sw
      type(drainage_config_t),      intent(in)    :: config_drain
      integer,                      intent(in)    :: numnod

      ! Defensive guards mirroring surface_water_config_validate.
      ! swman is allocatable; slice 1:nmper covers the active management periods.
      if (config_sw%swsrf == 3 .or. config_sw%swsec == 1 .or. config_sw%swqhr == 2) then
         call fatalerr_collected('surfacewater_state_init', &
            'swsrf=3, swsec=1, or swqhr=2 not supported on the TOML path')
         return
      end if
      if (allocated(config_sw%swman)) then
         if (any(config_sw%swman(1:config_sw%nmper) == 2)) then
            call fatalerr_collected('surfacewater_state_init', &
               'swman=2 (automatic weir) not supported on the TOML path')
            return
         end if
      end if

      ! Legacy global write retained for now: bocodre reads wlp via `use variables`
      ! for the primary surface water level (not surfacewater-state owned). A separate
      ! arc will migrate readers; this write stays until then.
      block
         use variables, only: wlp
         wlp = 0.0_real64  ! swsrf=2 has no primary system
      end block

      ! ---- L1: zero defaults ----
      self%numadj = 0
      self%wlsbak = 0.0_real64

      ! ---- L2: config-derived seeds ----
      ! wls1 = wlact - altcu (legacy rddre line; altcu=0 enforced by drainage_config_validate
      ! so this equals wlact). Inlined here — retires the wls1_init transient buffer
      ! indirection at the call site once the hoist lands in Task 5.
      self%wls    = config_sw%wlact - config_drain%altcu
      self%wlstar = self%wls

      ! Allocate per-level arrays.
      allocate(self%cqdrain    (config_drain%nrlevs));            self%cqdrain    = 0.0_real64
      allocate(self%cqdrainin  (config_drain%nrlevs));            self%cqdrainin  = 0.0_real64
      allocate(self%cqdrainout (config_drain%nrlevs));            self%cqdrainout = 0.0_real64
      allocate(self%inqdra     (config_drain%nrlevs, numnod));    self%inqdra     = 0.0_real64
      allocate(self%inqdra_in  (config_drain%nrlevs, numnod));    self%inqdra_in  = 0.0_real64
      allocate(self%inqdra_out (config_drain%nrlevs, numnod));    self%inqdra_out = 0.0_real64

   end subroutine surfacewater_state_init

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
