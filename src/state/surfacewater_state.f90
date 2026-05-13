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
   public :: swstlev_from_table

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

      ! ---- L1: zero defaults ----
      self%numadj = 0
      self%wlsbak = 0.0_real64

      ! ---- L2: config-derived seeds ----
      ! wls1 = wlact - altcu (legacy rddre line; altcu=0 enforced by drainage_config_validate
      ! so this equals wlact). Inlined directly from typed config — retires the legacy
      ! `wls1_init` transient buffer once the call-site hoist lands.
      self%wls    = config_sw%wlact - config_drain%altcu
      self%wlstar = self%wls

      ! Allocate per-level arrays.
      allocate(self%cqdrain    (config_drain%nrlevs));            self%cqdrain    = 0.0_real64
      allocate(self%cqdrainin  (config_drain%nrlevs));            self%cqdrainin  = 0.0_real64
      allocate(self%cqdrainout (config_drain%nrlevs));            self%cqdrainout = 0.0_real64
      allocate(self%inqdra     (config_drain%nrlevs, numnod));    self%inqdra     = 0.0_real64
      allocate(self%inqdra_in  (config_drain%nrlevs, numnod));    self%inqdra_in  = 0.0_real64
      allocate(self%inqdra_out (config_drain%nrlevs, numnod));    self%inqdra_out = 0.0_real64

      ! ---- L3: readswap-style shape math (sttab + swst init) ----
      !
      ! sttab(:,1): water-level rows.
      !   Row 1 = +100cm above soil surface (top).
      !   Row 2 = 0cm (soil surface).
      !   Rows 3..22: divide [0, zbotdr(1+nrpri)] into 20 compartments.
      ! For swsrf=2 (no primary system) nrpri=0, so zbotdr index is 1.
      block
         integer :: i, ilev, nrpri
         real(real64) :: wdepth, wvolum, wbreadth

         nrpri = 0

         self%sttab(1, 1) = 100.0_real64
         self%sttab(2, 1) =   0.0_real64
         do i = 3, 22
            self%sttab(i, 1) = config_drain%zbotdr(1 + nrpri) * real(i - 2, real64) / 20.0_real64
         end do

         ! sttab(:,2): storage volume per unit area (cm), summed across
         ! open-channel levels (swdtyp=0). Verbatim port from legacy rddre
         ! (readswap.f90:4878-4897). l(:) is in centimetres (D6 conversion
         ! done at TOML read time).
         do i = 1, 22
            self%sttab(i, 2) = 0.0_real64
            do ilev = 1 + nrpri, config_drain%nrlevs
               if (config_drain%swdtyp(ilev) == 0 .and. self%sttab(i, 1) > config_drain%zbotdr(ilev)) then
                  if (self%sttab(i, 1) <= 0.0_real64) then
                     ! Trapezium below soil surface
                     wdepth = self%sttab(i, 1) - config_drain%zbotdr(ilev)
                     wvolum = wdepth * (config_drain%widthr(ilev) + wdepth / config_drain%taludr(ilev))
                  else
                     ! Trapezium up to surface, plus rectangle above
                     wdepth   = -config_drain%zbotdr(ilev)
                     wvolum   = wdepth * (config_drain%widthr(ilev) + wdepth / config_drain%taludr(ilev))
                     wbreadth = config_drain%widthr(ilev) + 2.0_real64 * wdepth / config_drain%taludr(ilev)
                     wdepth   = self%sttab(i, 1)
                     wvolum   = wvolum + wbreadth * wdepth
                  end if
                  self%sttab(i, 2) = self%sttab(i, 2) + wvolum / config_drain%l(ilev)
               end if
            end do
         end do
      end block

      ! Initial storage state derived from sttab + wls.
      self%swstini = swstlev_from_table(self%sttab, self%wls)
      self%swst    = self%swstini

      ! ---- Post-init shape math (was in SurfaceWater(task=1) post-call block) ----
      ! hwlman/vtair: only output reads them; default to zero here.
      self%hwlman = 0.0_real64
      self%vtair  = 0.0_real64

      ! ZDraBas (macropore drainage basis): on the TOML pipeline only swsec=2 is
      ! reachable (validator rejects swsec=1; macropore-retire ADR 0040 makes
      ! NumLevRapDra=0 which retires the swdtyp(NumLevRapDra)=1 drain-tube
      ! branch). Defensive guard rejects anything else; in-scope assignment is
      ! ZDraBas := wlstar (the value seeded in L2).
      if (config_sw%swsec /= 2) then
         call fatalerr_collected('surfacewater_state_init', &
            'swsec /= 2 not supported on the TOML path (validator rejects)')
         return
      end if
      self%ZDraBas      = self%wlstar
      self%flInitDraBas = .false.

      ! ---- Legacy global write retained for now ----
      ! bocodre reads wlp via `use variables` for the primary surface water level
      ! (not surfacewater-state owned). A separate arc will migrate readers;
      ! this write stays until then. Block-scoped `use` confines the legacy
      ! import to these two lines.
      block
         use variables, only: wlp
         wlp = 0.0_real64  ! swsrf=2 has no primary system
      end block

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

   !> Pure-table-lookup variant of `swstlev` — operates directly on the storage
   !! table without needing a full `swap_state_t`. Used by `surfacewater_state_t%init`
   !! where only the partially-populated state is available, and by the
   !! `swstlev(state, wlev)` wrapper in `surfacewater_utils`.
   !!
   !! Lives here (rather than `surfacewater_utils`) to avoid a circular
   !! module dependency: `surfacewater_utils` already depends on
   !! `swap_state_mod` (which depends on this module).
   function swstlev_from_table(sttab, wlev) result(swstlev_r)
      implicit none
      real(real64), intent(in) :: sttab(22, 2)
      real(real64), intent(in) :: wlev
      real(real64) :: swstlev_r

      integer :: i
      real(real64) :: dwl
      character(len=200) :: messag

      if (wlev < sttab(22,1)) then
         messag = 'Surface water storage below bottom of table'
         call fatalerr_collected('swstlev_from_table', messag)
      end if
      if (wlev > sttab(1,1)) then
         messag = 'Surface water storage above top of table'
         call fatalerr_collected('swstlev_from_table', messag)
      end if

      i = 0
      do
         i = i + 1
         if (wlev >= sttab(i+1,1) .and. wlev <= sttab(i,1)) exit
      end do

      dwl = (wlev - sttab(i+1,1)) / (sttab(i,1) - sttab(i+1,1))
      swstlev_r = sttab(i+1,2) + dwl * (sttab(i,2) - sttab(i+1,2))
   end function swstlev_from_table

   !> Zero the reservoir-cumulative cohort — flzerocumu gate, flSurfaceWater partition.
   !! Three scalars; no allocatables.
   subroutine surfacewater_reset_cumulative_reservoir(self)
      class(surfacewater_state_t), intent(inout) :: self
      self%cqdrd  = 0.0_real64
      self%cwsupp = 0.0_real64
      self%cwout  = 0.0_real64
   end subroutine surfacewater_reset_cumulative_reservoir

end module surfacewater_state_mod
