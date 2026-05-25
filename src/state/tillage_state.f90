!> @file tillage_state.f90
!! Typed state record for the tillage subsystem (ADR 0039, Task T-1).
!!
!! Holds runtime-state fields owned by tillage.f90:
!!
!!   Group C -- 7 per-layer allocatables (sized numlay):
!!     Rho_tillage, Rho_cons, Rho_last, K_R_cons,
!!     Rho_match, N_match, Slope_match
!!
!!   Group D -- 3 per-step scalars:
!!     sumDWC, sumAvail1, sumAvail2
!!
!!   Group E -- 3 init-once geometry/cursor integers:
!!     MaxNumSoilHo, MaxNumSoilCP, iTill
!!
!!   Group AB -- init-once derived from config%soil%tillage by
!!     state%tillage%init; read-only at runtime.
!!     [GR-CROP 2026-05-25] migrated from variables.f90 till_* SAVE-state.
!!     [GR-SEED 2026-05-25 Task 9] config seeding absorbed from apply_soil_tillage.
!!     Ntill, Ntypes, i_n_model, iRedist, Max_Z_tillage, (swtill via state%cfg%soil%swtill)
!!     Date_tillage(:), Z_tillage(:), I_tillage(:), Type_Tillage(:),
!!     iType_Tillage(:), iTT1(:), iTT2(:),
!!     TAB_Rho_tillage(:), TAB_Rho_cons(:), TAB_K_R_cons(:),
!!     TAB_Rho_match(:), TAB_N_match(:)
!!
!! state%tillage%init(config_tillage, tend, numlay) allocates the Group C
!! arrays, zeroes all fields, then (if events are present) populates
!! Group AB from the typed config. Called from swap_mod after CalcGrid.
!! The swtill==1 gate lives inside the caller (swap_mod); when swtill/=1
!! the caller passes an unallocated config_tillage -- the init handles this
!! gracefully by checking allocated(config_tillage%events).
!!
!! Note: field names drop the legacy till_ prefix; the type name
!! tillage_state_t provides the namespace. Legacy name mappings:
!!   till_Rho_tillage -> Rho_tillage   till_sumDWC        -> sumDWC
!!   till_Rho_cons    -> Rho_cons       till_sumAvail1     -> sumAvail1
!!   till_Rho_last    -> Rho_last       till_sumAvail2     -> sumAvail2
!!   till_K_R_cons    -> K_R_cons       till_MaxNumSoilHo  -> MaxNumSoilHo
!!   till_Rho_match   -> Rho_match      till_MaxNumSoilCP  -> MaxNumSoilCP
!!   till_N_match     -> N_match        till_iTill         -> iTill
!!   till_Slope_match -> Slope_match
!!
!! See ADR 0039, docs/superpowers/specs/2026-05-12-state-migration-tillage-design.md
!!     docs/superpowers/plans/2026-05-12-tillage-state-migration.md

module tillage_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   use iso_c_binding, only: c_double
   implicit none
   private
   public :: tillage_state_t

   type :: tillage_state_t

      ! Group C -- per-layer allocatables (sized numlay)
      ! Allocated by init; unallocated until then.
      real(real64), allocatable :: Rho_tillage(:)  !! bulk density after tillage (g/cm3) per layer
      real(real64), allocatable :: Rho_cons(:)     !! consolidated bulk density (g/cm3) per layer
      real(real64), allocatable :: Rho_last(:)     !! bulk density at previous event (g/cm3) per layer
      real(real64), allocatable :: K_R_cons(:)     !! consolidation rate coefficient (-) per layer
      real(real64), allocatable :: Rho_match(:)    !! reference bulk density for matching (g/cm3) per layer
      real(real64), allocatable :: N_match(:)      !! exponent for bulk-density recovery (-) per layer
      real(real64), allocatable :: Slope_match(:)  !! slope for bulk-density recovery function per layer

      ! Group D -- per-step scalars (reset each timestep by tillage.f90)
      real(real64) :: sumDWC    = 0.0_real64  !! cumulative drainage-weighted consolidation (-)
      real(real64) :: sumAvail1 = 0.0_real64  !! cumulative available water fraction, path 1 (-)
      real(real64) :: sumAvail2 = 0.0_real64  !! cumulative available water fraction, path 2 (-)

      ! Group E -- init-once geometry/cursor integers
      integer :: MaxNumSoilHo = 0  !! max number of soil horizons in tillage table
      integer :: MaxNumSoilCP = 0  !! max number of soil-consolidation-parameter rows
      integer :: iTill        = 0  !! event-table cursor (index into tillage event array)

      ! Group AB -- init-once derived from state%cfg%soil%tillage by
      ! state%tillage%init. Read-only at runtime.
      ! [GR-CROP 2026-05-25] migrated from variables.f90 till_* legacy globals.
      ! [GR-SEED 2026-05-25 Task 9] formerly populated by apply_soil_tillage.
      ! Note: the swtill switch lives on state%cfg%soil%swtill (canonical),
      ! not duplicated here.
      integer :: i_n_model     = 2       !! n-parameter treatment switch (1..3)
      integer :: iRedist       = 2       !! Redistribution type after MvG change
      integer :: Ntill         = 0       !! Number of tabulated tillage events
      integer :: Ntypes        = 0       !! Number of tillage types
      real(real64) :: Max_Z_tillage = 0.0_real64  !! Max possible depth of tillage (cm)

      ! Per-event arrays (size Ntill, or Ntill+1 for Date_tillage with sentinel):
      real(real64), allocatable :: Date_tillage(:)  !! Tillage dates (days-since-1900); Date_tillage(Ntill+1) is tend+1 sentinel
      real(real64), allocatable :: Z_tillage(:)     !! Tillage depths (cm)
      real(real64), allocatable :: I_tillage(:)     !! Tillage intensity (0-1)
      integer,      allocatable :: Type_Tillage(:)  !! Tillage type index (refers to types(:))
      integer,      allocatable :: iTT1(:)          !! First position per type in iType_Tillage
      integer,      allocatable :: iTT2(:)          !! Last position per type in iType_Tillage

      ! Per-type arrays (size Ntypes):
      integer,      allocatable :: iType_Tillage(:)    !! Tillage type identifier (sequential index)
      real(real64), allocatable :: TAB_Rho_tillage(:)  !! Bulk density after tillage per type
      real(real64), allocatable :: TAB_Rho_cons(:)     !! Consolidated bulk density per type
      real(real64), allocatable :: TAB_K_R_cons(:)     !! Consolidation rate constant per type
      real(real64), allocatable :: TAB_Rho_match(:)    !! Matching-point density per type (i_n_model=3 only)
      real(real64), allocatable :: TAB_N_match(:)      !! Matching-point n per type (i_n_model=3 only)

      !> [SS-BMI2] Tillage output row buffer (DoTillage task=3 stream).
      !! N = 5: t1900, nraida, sumDWC, sumAvail1, sumAvail2.
      !! Debug writes to units 222/226 are gated by headless.
      real(c_double),    allocatable :: output_row(:)
      character(len=32), allocatable :: output_columns(:)
      integer                        :: output_n_cols = 0

   contains
      procedure :: init => tillage_state_init

   end type tillage_state_t

contains

   !> Type-bound initialiser for tillage_state_t.
   !!
   !! [GR-SEED 2026-05-25 Task 9] Absorbs the legacy free tillage_init body
   !! (Group C zero-fill / nlay alloc) plus the body of apply_soil_tillage
   !! (Group AB config seeding, ISO date conversion, sentinel, Max_Z_tillage,
   !! iTT1/iTT2 derivation). The swtill==1 early-return gate is handled by the
   !! caller (swap_mod): the caller passes config_tillage unconditionally; if
   !! events are not allocated (swtill/=1) the Group AB seeding is skipped.
   !!
   !! parse_iso_date_to_days1900 relocated here as a private helper (was in
   !! config_to_variables_mod); it had only one caller (apply_soil_tillage).
   subroutine tillage_state_init(self, config_tillage, tend, numlay)
      use soil_config_mod, only: soil_tillage_t
      class(tillage_state_t), intent(inout) :: self
      type(soil_tillage_t),   intent(in)    :: config_tillage
      real(real64),           intent(in)    :: tend
      integer,                intent(in)    :: numlay

      integer :: i, j

      ! ------------------------------------------------------------------
      ! Piece A: legacy tillage_init body -- Group C alloc + zero-fill,
      ! Group D + E zero, Group AB scalars already set by type defaults.
      ! ------------------------------------------------------------------

      ! Group C -- allocate and zero per-layer arrays
      if (allocated(self%Rho_tillage))  deallocate(self%Rho_tillage)
      if (allocated(self%Rho_cons))     deallocate(self%Rho_cons)
      if (allocated(self%Rho_last))     deallocate(self%Rho_last)
      if (allocated(self%K_R_cons))     deallocate(self%K_R_cons)
      if (allocated(self%Rho_match))    deallocate(self%Rho_match)
      if (allocated(self%N_match))      deallocate(self%N_match)
      if (allocated(self%Slope_match))  deallocate(self%Slope_match)

      allocate(self%Rho_tillage(numlay));  self%Rho_tillage  = 0.0_real64
      allocate(self%Rho_cons(numlay));     self%Rho_cons     = 0.0_real64
      allocate(self%Rho_last(numlay));     self%Rho_last     = 0.0_real64
      allocate(self%K_R_cons(numlay));     self%K_R_cons     = 0.0_real64
      allocate(self%Rho_match(numlay));    self%Rho_match    = 0.0_real64
      allocate(self%N_match(numlay));      self%N_match      = 0.0_real64
      allocate(self%Slope_match(numlay));  self%Slope_match  = 0.0_real64

      ! Group D -- zero per-step scalars
      self%sumDWC    = 0.0_real64
      self%sumAvail1 = 0.0_real64
      self%sumAvail2 = 0.0_real64

      ! Group E -- zero init-once geometry/cursor integers
      self%MaxNumSoilHo = 0
      self%MaxNumSoilCP = 0
      self%iTill        = 0

      ! ------------------------------------------------------------------
      ! Piece B: Group AB config seeding (formerly apply_soil_tillage).
      ! Skip if events are not allocated (swtill /= 1).
      ! ------------------------------------------------------------------
      if (.not. allocated(config_tillage%events)) return
      if (size(config_tillage%events) == 0)       return

      self%i_n_model = config_tillage%i_n_model
      self%iRedist   = config_tillage%iRedist

      self%Ntill  = size(config_tillage%events)
      self%Ntypes = size(config_tillage%types)

      ! Per-event arrays (sentinel: Date_tillage(Ntill+1) = tend + 1).
      if (allocated(self%Date_tillage)) deallocate(self%Date_tillage); allocate(self%Date_tillage(self%Ntill+1))
      if (allocated(self%Z_tillage))    deallocate(self%Z_tillage);    allocate(self%Z_tillage(self%Ntill))
      if (allocated(self%I_tillage))    deallocate(self%I_tillage);    allocate(self%I_tillage(self%Ntill))
      if (allocated(self%Type_Tillage)) deallocate(self%Type_Tillage); allocate(self%Type_Tillage(self%Ntill))

      do i = 1, self%Ntill
         self%Z_tillage(i)    = config_tillage%events(i)%z
         self%I_tillage(i)    = config_tillage%events(i)%intensity
         self%Type_Tillage(i) = config_tillage%events(i)%type_id
         self%Date_tillage(i) = parse_iso_date_to_days1900(config_tillage%events(i)%date)
      end do
      self%Date_tillage(self%Ntill + 1) = tend + 1.0_real64

      ! Per-type arrays.
      if (allocated(self%iType_Tillage))   deallocate(self%iType_Tillage);   allocate(self%iType_Tillage(self%Ntypes))
      if (allocated(self%TAB_Rho_cons))    deallocate(self%TAB_Rho_cons);    allocate(self%TAB_Rho_cons(self%Ntypes))
      if (allocated(self%TAB_Rho_tillage)) deallocate(self%TAB_Rho_tillage); allocate(self%TAB_Rho_tillage(self%Ntypes))
      if (allocated(self%TAB_K_R_cons))    deallocate(self%TAB_K_R_cons);    allocate(self%TAB_K_R_cons(self%Ntypes))

      do i = 1, self%Ntypes
         self%iType_Tillage(i)   = config_tillage%types(i)%id
         self%TAB_Rho_cons(i)    = config_tillage%types(i)%rho_cons
         self%TAB_Rho_tillage(i) = config_tillage%types(i)%rho_tillage
         self%TAB_K_R_cons(i)    = config_tillage%types(i)%k_R
      end do

      if (self%i_n_model == 3) then
         if (allocated(self%TAB_Rho_match)) deallocate(self%TAB_Rho_match); allocate(self%TAB_Rho_match(self%Ntypes))
         if (allocated(self%TAB_N_match))   deallocate(self%TAB_N_match);   allocate(self%TAB_N_match(self%Ntypes))
         do i = 1, self%Ntypes
            self%TAB_Rho_match(i) = config_tillage%types(i)%rho_match
            self%TAB_N_match(i)   = config_tillage%types(i)%N_match
         end do
      end if

      self%Max_Z_tillage = maxval(self%Z_tillage(1:self%Ntill))

      ! iTT1 / iTT2: first/last position per tillage type in iType_Tillage.
      ! Replicates the loop from the legacy Read_Tillage subroutine verbatim.
      if (allocated(self%iTT1)) deallocate(self%iTT1); allocate(self%iTT1(self%Ntill)); self%iTT1 = 0
      if (allocated(self%iTT2)) deallocate(self%iTT2); allocate(self%iTT2(self%Ntill)); self%iTT2 = 0
      do j = 1, self%Ntill
         do i = 1, self%Ntypes
            if (self%iTT1(j) == 0 .and. self%iType_Tillage(i) == j) self%iTT1(j) = i
            if (self%iTT1(j) >  0 .and. self%iType_Tillage(i) == j) self%iTT2(j) = i
         end do
      end do

   end subroutine tillage_state_init

   !> Private helper: ISO 'YYYY-MM-DD' string -> days since 1900 (real64).
   !!
   !! [GR-SEED 2026-05-25 Task 9] Relocated from config_to_variables_mod where
   !! it had a single caller (apply_soil_tillage). Now a private helper here.
   !! Constructs a toml_datetime from the parsed date components and delegates
   !! to the existing parse_date_to_days1900 helper in toml_field_helpers_mod.
   function parse_iso_date_to_days1900(s) result(t)
      use tomlf, only: toml_datetime
      use toml_field_helpers_mod, only: parse_date_to_days1900
      character(len=*), intent(in) :: s
      real(real64) :: t
      type(toml_datetime) :: dtv
      integer :: y, m, d
      read(s, '(i4,1x,i2,1x,i2)') y, m, d
      dtv%date%year  = y
      dtv%date%month = m
      dtv%date%day   = d
      ! Leave dtv%time fields at default (-1) so the conversion treats it as a date-only.
      t = parse_date_to_days1900(dtv)
   end function parse_iso_date_to_days1900

end module tillage_state_mod
