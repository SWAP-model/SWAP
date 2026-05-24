!> @file tillage_state.f90
!! Typed state record for the tillage subsystem (ADR 0039, Task T-1).
!!
!! Holds runtime-state fields owned by tillage.f90:
!!
!!   Group C — 7 per-layer allocatables (sized numlay):
!!     Rho_tillage, Rho_cons, Rho_last, K_R_cons,
!!     Rho_match, N_match, Slope_match
!!
!!   Group D — 3 per-step scalars:
!!     sumDWC, sumAvail1, sumAvail2
!!
!!   Group E — 3 init-once geometry/cursor integers:
!!     MaxNumSoilHo, MaxNumSoilCP, iTill
!!
!!   Group AB — init-once derived from config%soil%tillage by
!!     apply_soil_tillage (config_to_variables); read-only at runtime.
!!     [GR-CROP 2026-05-25] migrated from variables.f90 till_* SAVE-state.
!!     Ntill, Ntypes, i_n_model, iRedist, Max_Z_tillage, (swtill via state%cfg%soil%swtill)
!!     Date_tillage(:), Z_tillage(:), I_tillage(:), Type_Tillage(:),
!!     iType_Tillage(:), iTT1(:), iTT2(:),
!!     TAB_Rho_tillage(:), TAB_Rho_cons(:), TAB_K_R_cons(:),
!!     TAB_Rho_match(:), TAB_N_match(:)
!!
!! tillage_init(tl, numlay) allocates the Group C arrays and zeroes all
!! fields. Called between atmosphere_init and DoTillage(1) at swap startup
!! (T-2 handles the call-site wiring). Group AB arrays are allocated by
!! apply_soil_tillage (config_to_variables) before tillage_init runs.
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
   public :: tillage_state_t, tillage_init

   type :: tillage_state_t

      ! Group C — per-layer allocatables (sized numlay)
      ! Allocated by tillage_init; unallocated until then.
      real(real64), allocatable :: Rho_tillage(:)  !! bulk density after tillage (g/cm3) per layer
      real(real64), allocatable :: Rho_cons(:)     !! consolidated bulk density (g/cm3) per layer
      real(real64), allocatable :: Rho_last(:)     !! bulk density at previous event (g/cm3) per layer
      real(real64), allocatable :: K_R_cons(:)     !! consolidation rate coefficient (-) per layer
      real(real64), allocatable :: Rho_match(:)    !! reference bulk density for matching (g/cm3) per layer
      real(real64), allocatable :: N_match(:)      !! exponent for bulk-density recovery (-) per layer
      real(real64), allocatable :: Slope_match(:)  !! slope for bulk-density recovery function per layer

      ! Group D — per-step scalars (reset each timestep by tillage.f90)
      real(real64) :: sumDWC    = 0.0_real64  !! cumulative drainage-weighted consolidation (-)
      real(real64) :: sumAvail1 = 0.0_real64  !! cumulative available water fraction, path 1 (-)
      real(real64) :: sumAvail2 = 0.0_real64  !! cumulative available water fraction, path 2 (-)

      ! Group E — init-once geometry/cursor integers
      integer :: MaxNumSoilHo = 0  !! max number of soil horizons in tillage table
      integer :: MaxNumSoilCP = 0  !! max number of soil-consolidation-parameter rows
      integer :: iTill        = 0  !! event-table cursor (index into tillage event array)

      ! Group AB — init-once derived from state%cfg%soil%tillage by
      ! apply_soil_tillage. Read-only at runtime.
      ! [GR-CROP 2026-05-25] migrated from variables.f90 till_* legacy globals.
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

   end type tillage_state_t

contains

   subroutine tillage_init(tl, numlay)
      !! Allocates the per-layer arrays and zeroes all scalar fields.
      !! Called between atmosphere_init and DoTillage(1) at swap startup.
      !! T-2 handles the call-site wiring in swap.f90.
      type(tillage_state_t), intent(inout) :: tl
      integer,               intent(in)    :: numlay

      ! Group C — allocate and zero per-layer arrays
      if (allocated(tl%Rho_tillage))  deallocate(tl%Rho_tillage)
      if (allocated(tl%Rho_cons))     deallocate(tl%Rho_cons)
      if (allocated(tl%Rho_last))     deallocate(tl%Rho_last)
      if (allocated(tl%K_R_cons))     deallocate(tl%K_R_cons)
      if (allocated(tl%Rho_match))    deallocate(tl%Rho_match)
      if (allocated(tl%N_match))      deallocate(tl%N_match)
      if (allocated(tl%Slope_match))  deallocate(tl%Slope_match)

      allocate(tl%Rho_tillage(numlay));  tl%Rho_tillage  = 0.0_real64
      allocate(tl%Rho_cons(numlay));     tl%Rho_cons     = 0.0_real64
      allocate(tl%Rho_last(numlay));     tl%Rho_last     = 0.0_real64
      allocate(tl%K_R_cons(numlay));     tl%K_R_cons     = 0.0_real64
      allocate(tl%Rho_match(numlay));    tl%Rho_match    = 0.0_real64
      allocate(tl%N_match(numlay));      tl%N_match      = 0.0_real64
      allocate(tl%Slope_match(numlay));  tl%Slope_match  = 0.0_real64

      ! Group D — zero per-step scalars
      tl%sumDWC       = 0.0_real64
      tl%sumAvail1    = 0.0_real64
      tl%sumAvail2    = 0.0_real64

      ! Group E — zero init-once geometry/cursor integers
      tl%MaxNumSoilHo = 0
      tl%MaxNumSoilCP = 0
      tl%iTill        = 0

      ! Group AB scalars — defaults are set on type declaration; Group AB
      ! arrays are allocated/populated separately by apply_soil_tillage
      ! (config_to_variables) before tillage_init runs in swap_mod, so we
      ! don't touch them here.

   end subroutine tillage_init

end module tillage_state_mod
