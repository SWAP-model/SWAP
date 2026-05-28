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
   use swap_array_dimensions, only: MABBC, MAHO, macp
   implicit none
   private
   public :: solute_state_t

   type :: solute_state_t

      ! === per-node arrays (macp-sized; allocated by caller from config) ===
      real(real64), allocatable :: cml(:)    !! soil solute concentration (M/L3 water) in mobile region
      real(real64), allocatable :: cmsy(:)   !! dissolved + adsorbed solute concentration (M/L3 soil volume)

      ! Derived time-invariant coefficients (filled once in solute_seed).
      real(real64), allocatable :: bdenskf(:)          !! bdens*kf per node (-)
      real(real64), allocatable :: bdenskfcref(:)      !! bdenskf*cref per node (M/L3)
      real(real64), allocatable :: bdenskfsatporos(:)  !! bdens*kfsat+poros per node (-)
      real(real64), allocatable :: ddiffwcs(:)         !! ddif/thetsl**2 per node (cm2/d)
      real(real64), allocatable :: decpotfdepth(:)     !! decpot*fdepth per node (1/d)

      ! === config-snapshot scalars (snapshotted from config%solute at config_to_variables) ===
      ! Top-level solute enable switch (snapshotted here so compute routines read state, not state%cfg).
      integer      :: swsolu  = 0           !! solute simulation enable (0=off, 1=on); snapshotted from config%solute%swsolu
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
      ! Initial-condition config (legacy multi-depth init).
      integer      :: nconc   = 0           !! number of initial-concentration depth points
      real(real64), allocatable :: cml_init(:)  !! initial mobile concentration table (M/L3) or per-node profile (swinco=3)
      real(real64), allocatable :: zc_init(:)   !! depth column for cml_init table (cm)

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

   contains
      procedure :: init               => solute_state_init
      procedure :: reset_intermediate => solute_reset_intermediate
      procedure :: reset_cumulative   => solute_reset_cumulative
   end type solute_state_t

contains

   !> Seed solute state from typed config + runtime dimension args.
   !!
   !! [GR-SEED 2026-05-25 Task 3] Absorbs:
   !!   - Scalar config-snapshot assignments previously in config_to_variables.f90
   !!     (Solute block, ~75 lines).
   !!   - Per-layer array broadcasts (ldis scalar→layer(1), kf/decpot/fdepth element copy).
   !!   - 2D cseeptab → interleaved afgen layout flatten.
   !!   - Runtime per-node cml/cmsy allocation + zero-fill formerly in the free
   !!     solute_init(state) in src/solute/solute.f90.
   !! [OD Step 8 2026-05-28] Absorbs:
   !!   - swinco=3 + swsolu=1 warm-restart cml profile load from
   !!     seed_state_from_config orchestrator (final cross-subsystem write).
   subroutine solute_state_init(self, config_solute, config_soil, numnod)
      use solute_config_mod, only: solute_config_t
      use soil_config_mod,   only: soil_config_t
      class(solute_state_t), intent(inout) :: self
      type(solute_config_t), intent(in)    :: config_solute
      type(soil_config_t),   intent(in)    :: config_soil
      integer,               intent(in)    :: numnod

      integer :: i, n

      ! ------------------------------------------------------------------
      ! Runtime zero-fill (formerly free solute_init in solute.f90):
      ! allocate per-node arrays sized to numnod and seed initial profiles.
      ! ------------------------------------------------------------------
      n = numnod
      if (.not. allocated(self%cml))  allocate(self%cml(n))
      if (.not. allocated(self%cmsy)) allocate(self%cmsy(n))
      if (.not. allocated(self%bdenskf))         allocate(self%bdenskf(n))
      if (.not. allocated(self%bdenskfcref))     allocate(self%bdenskfcref(n))
      if (.not. allocated(self%bdenskfsatporos)) allocate(self%bdenskfsatporos(n))
      if (.not. allocated(self%ddiffwcs))        allocate(self%ddiffwcs(n))
      if (.not. allocated(self%decpotfdepth))    allocate(self%decpotfdepth(n))
      self%bdenskf(:)         = 0.0_real64
      self%bdenskfcref(:)     = 0.0_real64
      self%bdenskfsatporos(:) = 0.0_real64
      self%ddiffwcs(:)        = 0.0_real64
      self%decpotfdepth(:)    = 0.0_real64

      ! swinco=3 (warm restart): load the per-node initial concentration profile
      ! from the cml_file CSV (OD Step 8: folded from seed_state_from_config).
      ! Other swinco values: solute task=1 interpolates from the (zc_init, cml_init)
      ! table; the initial cml here is just the default zero.
      if (config_soil%swinco == 3 .and. config_solute%swsolu == 1) then
         if (allocated(config_soil%initial%cml_file) .and. &
             len_trim(config_soil%initial%cml_file) > 0) then
            cml_load: block
               use soil_init_csv_mod, only: cml_profile_table_t
               use error_mod,         only: error_collection_t
               type(cml_profile_table_t) :: cml_tbl
               type(error_collection_t)  :: errs
               integer :: k, nrows
               call cml_tbl%load(trim(config_soil%initial%cml_file), errs)
               call errs%abort_if_fatal()
               ! abort_if_fatal() returns only on success, so cml_tbl%is_loaded
               ! must be .true. here. Guard anyway as defense-in-depth in case
               ! the error model softens later (e.g. non-fatal recoverable errors).
               if (.not. cml_tbl%is_loaded) exit cml_load
               nrows = size(cml_tbl%rows)
               self%nconc = nrows
               if (.not. allocated(self%cml_init)) then
                  allocate(self%cml_init(macp)); self%cml_init = 0.0_real64
               end if
               if (.not. allocated(self%zc_init)) then
                  allocate(self%zc_init(macp));  self%zc_init  = 0.0_real64
               end if
               do k = 1, nrows
                  self%zc_init(k)  = cml_tbl%rows(k)%z
                  self%cml_init(k) = cml_tbl%rows(k)%cml
               end do
            end block cml_load
         end if
      end if

      if (allocated(self%cml_init)) then
         self%cml(:)  = self%cml_init(1:n)
      else
         self%cml(:)  = 0.0_real64
      end if
      self%cmsy(:) = 0.0_real64

      ! ------------------------------------------------------------------
      ! Config-snapshot scalars (formerly adapter Solute block):
      ! ------------------------------------------------------------------
      self%swsolu  = config_solute%swsolu
      self%swbotbc = config_solute%swbotbc
      self%cdrain  = config_solute%cdrain
      self%tscf    = config_solute%tscf
      self%rtheta  = config_solute%rtheta
      self%bexp    = config_solute%bexp
      self%cref    = config_solute%cref
      self%cpre    = config_solute%cpre
      self%ddif    = config_solute%ddif
      self%frexp   = config_solute%frexp
      self%gampar  = config_solute%gampar
      self%daquif  = config_solute%daquif
      self%kfsat   = config_solute%kfsat
      self%decsat  = config_solute%decsat
      self%poros   = config_solute%poros
      self%swbr    = config_solute%swbr

      ! ------------------------------------------------------------------
      ! Per-layer arrays (MAHO-sized):
      ! ------------------------------------------------------------------
      ! Dispersion length: element-wise copy when array; otherwise broadcast
      ! scalar to layer(1) (legacy rdsdor 'ldis' single-value behaviour).
      if (.not. allocated(self%ldis)) then
         allocate(self%ldis(MAHO)); self%ldis = 0.0_real64
      end if
      if (allocated(config_solute%ldis_array)) then
         do i = 1, size(config_solute%ldis_array)
            self%ldis(i) = config_solute%ldis_array(i)
         end do
      else if (config_solute%ldis > 0.0_real64) then
         self%ldis(1) = config_solute%ldis
      end if

      if (.not. allocated(self%kf)) then
         allocate(self%kf(MAHO));     self%kf     = 0.0_real64
      end if
      if (.not. allocated(self%decpot)) then
         allocate(self%decpot(MAHO)); self%decpot = 0.0_real64
      end if
      if (.not. allocated(self%fdepth)) then
         allocate(self%fdepth(MAHO)); self%fdepth = 0.0_real64
      end if
      if (allocated(config_solute%kf)) then
         do i = 1, min(size(config_solute%kf), size(self%kf))
            self%kf(i) = config_solute%kf(i)
         end do
      end if
      if (allocated(config_solute%decpot)) then
         do i = 1, min(size(config_solute%decpot), size(self%decpot))
            self%decpot(i) = config_solute%decpot(i)
         end do
      end if
      if (allocated(config_solute%fdepth)) then
         do i = 1, min(size(config_solute%fdepth), size(self%fdepth))
            self%fdepth(i) = config_solute%fdepth(i)
         end do
      end if

      ! ------------------------------------------------------------------
      ! cseeptab: flatten 2D typed config to the interleaved afgen layout.
      ! afgen(cseeptab, mabbc*2, time) reads pairs as (2*k-1)=time, (2*k)=value.
      ! ------------------------------------------------------------------
      if (.not. allocated(self%cseeptab)) then
         allocate(self%cseeptab(2*MABBC)); self%cseeptab = 0.0_real64
      end if
      if (allocated(config_solute%cseeptab)) then
         do i = 1, min(size(config_solute%cseeptab, 1), size(self%cseeptab)/2)
            self%cseeptab(2*i - 1) = config_solute%cseeptab(i, 1)   ! time
            self%cseeptab(2*i)     = config_solute%cseeptab(i, 2)   ! concentration
         end do
      end if

   end subroutine solute_state_init

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
