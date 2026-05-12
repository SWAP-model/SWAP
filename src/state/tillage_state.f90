!> @file tillage_state.f90
!! Typed state record for the tillage subsystem (ADR 0039, Task T-1).
!!
!! Holds 13 runtime-state fields owned by tillage.f90:
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
!! tillage_init(tl, numlay) allocates the Group C arrays and zeroes all
!! fields. Called between atmosphere_init and DoTillage(1) at swap startup
!! (T-2 handles the call-site wiring).
!!
!! Excluded (config, not runtime state):
!!   - Groups A+B (18 config-constant till_* fields already in typed
!!     soil_tillage_t) — deferred to future config-consolidation arc (D5).
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

   end subroutine tillage_init

end module tillage_state_mod
