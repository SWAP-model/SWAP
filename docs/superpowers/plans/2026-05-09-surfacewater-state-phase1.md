# Surface-Water State Migration — Phase 1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Move all 29 surface-water-owned globals into a typed `surfacewater_state_t`, threaded as part of `swap_state_t` from `swap_main` down through every surface-water entry point and through every output reader of those globals. Compute writes to typed state; output reads typed state; legacy globals preserved as a transitional dual-write so check-full byte-identical regression stays green throughout, then dual-write is dropped at the end of the phase.

**Architecture:** Two new modules (`src/state/surfacewater_state.f90`, `src/state/swap_state.f90`). `state` argument threaded through `SurfaceWater(task)`, `surfacewater_init`, `runoff()`, and their callers (`swap.f90`, `boundtop.f90`, `swapoutput.f90`, `swap_csv_output.f90`). ASSOCIATE blocks in compute bodies preserve variable names. `fldecdt` becomes an `intent(out)` flag on `SurfaceWater(task=2)` instead of a global write.

**Tech Stack:** Fortran 2008, gfortran, meson + ninja + pixi, pFUnit, check-full byte-identical regression.

**Spec:** `docs/superpowers/specs/2026-05-09-state-migration-surfacewater-design.md`

**Branch:** all commits go on `refactor/surfacewater-state`. **Do NOT merge to development until Phase 2 is complete and verification at end of Phase 2 is green.**

---

## File Structure

**Created:**
- `src/state/surfacewater_state.f90` — `surfacewater_state_t` definition (29 fields per spec D2)
- `src/state/swap_state.f90` — `swap_state_t` aggregator with one field (`surfacewater`)
- `tests/unit/state/test_surfacewater_state.pf` — pFUnit construction & defaults tests
- `tests/unit/state/test_swap_state.pf` — pFUnit aggregator tests

**Modified:**
- `src/drainage/surfacewater.f90` — `SurfaceWater(task)` signature gains `state`, `request_smaller_dt`; ASSOCIATE in `WLEVBAL`, `WBALLEV`; dual-write state and globals during transitional tasks; drop global writes at end of Phase 1
- `src/drainage/surfacewater_init.f90` — `surfacewater_init` signature gains `state`; populates state alongside globals
- `src/utils/surfacewaterutils.f90` — `runoff()` signature gains `state`
- `src/boundary/boundtop.f90` — caller of `runoff()` gets `state` passed in from its caller
- `src/core/swap.f90` — declares `type(swap_state_t) :: state`; passes through to `SurfaceWater`, `boundtop`, output routines; reads `request_smaller_dt`
- `src/core/initialize.f90` — `state` initialization (zero defaults, then `surfacewater_init` populates)
- `src/io/swap_csv_output.f90` — `set_values` reads surface-water fields from `state%surfacewater` instead of globals; `csv_out` signature gains `state`
- `src/io/swapoutput.f90` — `outdrf`, `outswb`, `outend`, `swapoutput` dispatcher signatures gain `state`; reads from `state%surfacewater`; `swstini` year-reset moves to a dedicated `surfacewater_year_reset(state%surfacewater)` callback invoked from `swap_main`
- `tests/unit/drainage/test_surfacewater_init.pf` — migrate fixtures to construct `swap_state_t` and assert on `state%surfacewater%*` fields
- `tests/unit/io/toml/test_surfacewater_parity.pf` — same migration
- meson source list — add the two new state modules

**No changes to** (in this phase):
- `src/core/variables.f90` — globals stay declared until Phase 2
- `src/drainage/drainage.f90` — `qdrain` rule relocation is Phase 2
- `src/io/toml/config_to_variables.f90` — adapter unchanged

---

### Task 1: Create `surfacewater_state_t`

Define the typed state module. Pure data type with default-initialized scalars and unallocated arrays. The arrays are allocated by `surfacewater_init` in Task 4 once `Madr` and `macp` (drainage-config dimensions) are known.

**Files:**
- Create: `src/state/surfacewater_state.f90`
- Create: `tests/unit/state/test_surfacewater_state.pf`
- Modify: `tests/unit/testSuites.inc` (add `ADD_TEST_SUITE(test_surfacewater_state_suite)`)
- Modify: meson source list (add `src/state/surfacewater_state.f90`)

- [ ] **Step 1: Locate the meson source list**

Run: `grep -rln "src/utils/array_utils.f90\|surfacewaterutils.f90" --include="meson.build" . 2>/dev/null`
Expected: paths to the meson.build files that list source. Note them for Step 5.

- [ ] **Step 2: Write the failing pFUnit tests**

Create `tests/unit/state/test_surfacewater_state.pf`:

```fortran
! Tests for surfacewater_state_t — typed state record for the
! surface-water subsystem. Phase 1 of the state-migration arc.
!
! See docs/superpowers/specs/2026-05-09-state-migration-surfacewater-design.md.

@test
subroutine test_surfacewater_state_default_scalars()
   use funit
   use iso_fortran_env, only: real64
   use surfacewater_state_mod, only: surfacewater_state_t
   type(surfacewater_state_t) :: sw

   @assertEqual(0.0_real64, sw%wls,           1.0e-12_real64)
   @assertEqual(0.0_real64, sw%wlstar,        1.0e-12_real64)
   @assertEqual(0.0_real64, sw%swst,          1.0e-12_real64)
   @assertEqual(0.0_real64, sw%swstini,       1.0e-12_real64)
   @assertEqual(0.0_real64, sw%hwlman,        1.0e-12_real64)
   @assertEqual(0.0_real64, sw%vtair,         1.0e-12_real64)
   @assertEqual(0.0_real64, sw%wlsold,        1.0e-12_real64)
   @assertEqual(0.0_real64, sw%cqdrd,         1.0e-12_real64)
   @assertEqual(0.0_real64, sw%cwsupp,        1.0e-12_real64)
   @assertEqual(0.0_real64, sw%cwout,         1.0e-12_real64)
   @assertEqual(0.0_real64, sw%cqdra,         1.0e-12_real64)
   @assertEqual(0.0_real64, sw%ZDraBas,       1.0e-12_real64)
   @assertEqual(0.0_real64, sw%iqdra,         1.0e-12_real64)
   @assertEqual(0.0_real64, sw%qdrtot,        1.0e-12_real64)
   @assertFalse(sw%overfl)
   @assertTrue(sw%flInitDraBas)
   @assertEqual(1, sw%imper)
   @assertEqual(0, sw%numadj)
end subroutine test_surfacewater_state_default_scalars

@test
subroutine test_surfacewater_state_wlsbak_default()
   use funit
   use iso_fortran_env, only: real64
   use surfacewater_state_mod, only: surfacewater_state_t
   type(surfacewater_state_t) :: sw
   integer :: i

   do i = 1, 4
      @assertEqual(0.0_real64, sw%wlsbak(i), 1.0e-12_real64)
   end do
end subroutine test_surfacewater_state_wlsbak_default

@test
subroutine test_surfacewater_state_sttab_default()
   use funit
   use iso_fortran_env, only: real64
   use surfacewater_state_mod, only: surfacewater_state_t
   type(surfacewater_state_t) :: sw
   integer :: i, j

   do j = 1, 2
      do i = 1, 22
         @assertEqual(0.0_real64, sw%sttab(i,j), 1.0e-12_real64)
      end do
   end do
end subroutine test_surfacewater_state_sttab_default

@test
subroutine test_surfacewater_state_arrays_unallocated()
   ! Per-level arrays are allocatable; default state is unallocated.
   ! surfacewater_init allocates them based on drainage config (Madr, macp).
   use funit
   use surfacewater_state_mod, only: surfacewater_state_t
   type(surfacewater_state_t) :: sw

   @assertFalse(allocated(sw%cqdrain))
   @assertFalse(allocated(sw%cqdrainin))
   @assertFalse(allocated(sw%cqdrainout))
   @assertFalse(allocated(sw%qdra))
   @assertFalse(allocated(sw%inqdra))
   @assertFalse(allocated(sw%inqdra_in))
   @assertFalse(allocated(sw%inqdra_out))
end subroutine test_surfacewater_state_arrays_unallocated

@test
subroutine test_surfacewater_state_array_allocation()
   ! Caller (surfacewater_init) is responsible for allocating per-level
   ! arrays. Verify the allocation pattern works.
   use funit
   use iso_fortran_env, only: real64
   use surfacewater_state_mod, only: surfacewater_state_t
   type(surfacewater_state_t) :: sw
   integer, parameter :: Madr = 5
   integer, parameter :: macp = 100

   allocate(sw%cqdrain(Madr));         sw%cqdrain    = 0.0_real64
   allocate(sw%cqdrainin(Madr));       sw%cqdrainin  = 0.0_real64
   allocate(sw%cqdrainout(Madr));      sw%cqdrainout = 0.0_real64
   allocate(sw%qdra(Madr, macp));      sw%qdra       = 0.0_real64
   allocate(sw%inqdra(Madr, macp));    sw%inqdra     = 0.0_real64
   allocate(sw%inqdra_in(Madr, macp)); sw%inqdra_in  = 0.0_real64
   allocate(sw%inqdra_out(Madr, macp));sw%inqdra_out = 0.0_real64

   @assertTrue(allocated(sw%cqdrain))
   @assertEqual(Madr, size(sw%cqdrain))
   @assertEqual(Madr, size(sw%qdra, 1))
   @assertEqual(macp, size(sw%qdra, 2))
end subroutine test_surfacewater_state_array_allocation
```

Add `ADD_TEST_SUITE(test_surfacewater_state_suite)` to `tests/unit/testSuites.inc` near the other state/utils suites.

- [ ] **Step 3: Run tests — confirm FAIL**

Run: `pixi run -e test test-pfunit 2>&1 | tail -10`
Expected: 5 failures referencing missing module `surfacewater_state_mod`.

- [ ] **Step 4: Implement the module**

Create `src/state/surfacewater_state.f90`:

```fortran
!> @file surfacewater_state.f90
!! Typed state record for the surface-water subsystem. Owns the 29
!! variables that legacy SWAP held in the `variables.f90` globals
!! module. See ADR 0030 (state-migration pilot) and the
!! 2026-05-09 design spec for rationale and field provenance.

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
      real(real64), allocatable :: qdra(:,:)         ! (Madr, macp)
      real(real64), allocatable :: inqdra(:,:)       ! (Madr, macp)
      real(real64), allocatable :: inqdra_in(:,:)    ! (Madr, macp)
      real(real64), allocatable :: inqdra_out(:,:)   ! (Madr, macp)
   end type surfacewater_state_t

end module surfacewater_state_mod
```

- [ ] **Step 5: Add the source to meson.build**

In whichever `meson.build` lists `src/utils/array_utils.f90` (Step 1's grep result), add `src/state/surfacewater_state.f90` next to it. If `src/state/` is a new directory, also add the directory to whichever `subdir(...)` list aggregates source folders.

- [ ] **Step 6: Run tests — confirm PASS**

Run: `pixi run -e test test-pfunit 2>&1 | tail -10`
Expected: zero failures. Test count = previous baseline + 5.

- [ ] **Step 7: Commit**

```bash
git add src/state/surfacewater_state.f90 tests/unit/state/test_surfacewater_state.pf tests/unit/testSuites.inc
# Plus the meson.build files modified in Step 5
git status
git add <the meson files you edited>
git commit -m "$(cat <<'EOF'
feat(state): SS-SWST Task 1 — surfacewater_state_t typed state record

Phase 1 foundation. Defines the surface-water subsystem's typed
state record with the 29 fields from the design spec. Per-level
arrays are allocatable; surfacewater_init will allocate them in a
later task once drainage-config dimensions (Madr, macp) are
visible.

Spec: docs/superpowers/specs/2026-05-09-state-migration-surfacewater-design.md
EOF
)"
```

---

### Task 2: Create `swap_state_t` aggregator

The top-level state record. Phase 1 of the migration adds one field (`surfacewater`); subsequent subsystem-migration arcs will add more.

**Files:**
- Create: `src/state/swap_state.f90`
- Create: `tests/unit/state/test_swap_state.pf`
- Modify: `tests/unit/testSuites.inc`
- Modify: meson source list

- [ ] **Step 1: Write the failing pFUnit tests**

Create `tests/unit/state/test_swap_state.pf`:

```fortran
! Tests for swap_state_t — top-level aggregator of all per-subsystem
! typed state records. Phase 1 starts with only the surface-water
! subsystem migrated. Each subsequent migration arc adds a field.

@test
subroutine test_swap_state_has_surfacewater()
   use funit
   use iso_fortran_env, only: real64
   use swap_state_mod, only: swap_state_t
   type(swap_state_t) :: state

   @assertEqual(0.0_real64, state%surfacewater%wls,    1.0e-12_real64)
   @assertEqual(0.0_real64, state%surfacewater%swst,   1.0e-12_real64)
   @assertEqual(1,           state%surfacewater%imper)
end subroutine test_swap_state_has_surfacewater

@test
subroutine test_swap_state_independent_instances()
   ! Confirm that two state instances are independent (no shared module
   ! storage). This is the architectural intent — each "simulation"
   ! holds its own state.
   use funit
   use iso_fortran_env, only: real64
   use swap_state_mod, only: swap_state_t
   type(swap_state_t) :: state_a, state_b

   state_a%surfacewater%wls = 12.5_real64
   state_b%surfacewater%wls = -3.1_real64

   @assertEqual(12.5_real64, state_a%surfacewater%wls, 1.0e-12_real64)
   @assertEqual(-3.1_real64, state_b%surfacewater%wls, 1.0e-12_real64)
end subroutine test_swap_state_independent_instances
```

Add `ADD_TEST_SUITE(test_swap_state_suite)` to `tests/unit/testSuites.inc`.

- [ ] **Step 2: Run tests — confirm FAIL**

Run: `pixi run -e test test-pfunit 2>&1 | tail -10`
Expected: 2 failures referencing missing module `swap_state_mod`.

- [ ] **Step 3: Implement**

Create `src/state/swap_state.f90`:

```fortran
!> @file swap_state.f90
!! Top-level aggregator of all per-subsystem typed state records.
!! Phase 1 of the state-migration arc starts with only the
!! surface-water subsystem; subsequent migration arcs (water-balance,
!! crop, drainage, atmosphere, …) add their own fields here.
!!
!! Threaded through subroutine signatures from `swap_main` down,
!! replacing the implicit shared state held in `variables.f90`.

module swap_state_mod
   use surfacewater_state_mod, only: surfacewater_state_t
   implicit none
   private
   public :: swap_state_t

   type :: swap_state_t
      type(surfacewater_state_t) :: surfacewater
      ! Subsequent migration arcs add fields here.
   end type swap_state_t

end module swap_state_mod
```

- [ ] **Step 4: Add to meson**

Add `src/state/swap_state.f90` to the same meson list as `surfacewater_state.f90`.

- [ ] **Step 5: Run tests — confirm PASS**

Run: `pixi run -e test test-pfunit 2>&1 | tail -10`
Expected: zero failures. Test count = previous + 2.

- [ ] **Step 6: Commit**

```bash
git add src/state/swap_state.f90 tests/unit/state/test_swap_state.pf tests/unit/testSuites.inc <meson files>
git commit -m "$(cat <<'EOF'
feat(state): SS-SWST Task 2 — swap_state_t aggregator

Phase 1 foundation. Aggregates per-subsystem typed state records;
starts with one field (surfacewater). Subsequent subsystem
migrations grow this type one field at a time.

Spec: docs/superpowers/specs/2026-05-09-state-migration-surfacewater-design.md
EOF
)"
```

---

### Task 3: Refactor `surfacewater_init` to dual-write state and globals

Change the entry point to take `state` and populate both the typed state and the legacy globals (transitional dual-write). check-full and existing tests stay green because globals are still authoritative for any code that hasn't been migrated yet.

**Files:**
- Modify: `src/drainage/surfacewater_init.f90`
- Modify: `src/core/swap.f90` — caller of `surfacewater_init` now constructs and passes `state`
- Modify: `tests/unit/drainage/test_surfacewater_init.pf` — migrate fixtures to construct `state`, assert on `state%surfacewater%*`

- [ ] **Step 1: Read the current `surfacewater_init` and its caller**

Run: `cat src/drainage/surfacewater_init.f90`
Run: `grep -n "surfacewater_init" src/core/swap.f90`
Note: the `state` arg should pass through `swap_main → SurfaceWater(1) → surfacewater_init` (or whatever the current call chain is). Your job in this task is just `surfacewater_init` itself; the calling chain change is in Task 4.

- [ ] **Step 2: Refactor the signature**

Change the signature of `surfacewater_init` from:

```fortran
subroutine surfacewater_init()
```

(or whatever it currently is — read the file) to:

```fortran
subroutine surfacewater_init(state)
   use swap_state_mod, only: swap_state_t
   type(swap_state_t), intent(inout) :: state
```

Then inside the body, for every place that today writes a legacy global owned by surfacewater (the 29 variables — see spec D2), add a corresponding `state%surfacewater%<field> = …` write. **Both writes happen.** The legacy global write is preserved.

Allocate the per-level arrays in `state%surfacewater` against the drainage-config dimensions:

```fortran
! after determining nrlevs (the drain-level count, == Madr in legacy parlance)
! and after numnod is known (== macp):
allocate(state%surfacewater%cqdrain(nrlevs));         state%surfacewater%cqdrain    = 0.0_real64
allocate(state%surfacewater%cqdrainin(nrlevs));       state%surfacewater%cqdrainin  = 0.0_real64
allocate(state%surfacewater%cqdrainout(nrlevs));      state%surfacewater%cqdrainout = 0.0_real64
allocate(state%surfacewater%qdra(nrlevs, numnod));    state%surfacewater%qdra       = 0.0_real64
allocate(state%surfacewater%inqdra(nrlevs, numnod));  state%surfacewater%inqdra     = 0.0_real64
allocate(state%surfacewater%inqdra_in(nrlevs, numnod));  state%surfacewater%inqdra_in  = 0.0_real64
allocate(state%surfacewater%inqdra_out(nrlevs, numnod)); state%surfacewater%inqdra_out = 0.0_real64
```

Use ASSOCIATE for clarity in the body if it helps:

```fortran
associate(sw => state%surfacewater)
   sw%wls    = ...      ! same RHS as the legacy global write
   wls       = sw%wls   ! preserve legacy global
   ...
end associate
```

- [ ] **Step 3: Update the caller in `src/core/swap.f90`**

Find the existing `call surfacewater_init()` (or `call surfacewater_init(...)`) and change it to `call surfacewater_init(state)`. The `state` variable should be declared at the `swap_main` level — that's added in Task 4. For *this* task you can assume `state` exists in scope; if it doesn't compile, defer the caller change to Task 4 and stub the call temporarily — or, easier, declare a local `type(swap_state_t) :: state` at the top of `swap_main` for now and pass it to `surfacewater_init`. Task 4 builds out the rest.

- [ ] **Step 4: Migrate `tests/unit/drainage/test_surfacewater_init.pf`**

For each existing test that constructs and asserts on globals, refactor to construct a `state` instance, call `surfacewater_init(state)`, and assert on `state%surfacewater%*` instead of the global. Keep at least one assertion against the legacy global to verify the dual-write — this is a transitional test that protects against drift while we still have both targets.

Example:

```fortran
@test
subroutine test_surfacewater_init_zeros_wlsbak()
   use funit
   use iso_fortran_env, only: real64
   use swap_state_mod, only: swap_state_t
   use variables, only: wlsbak  ! legacy
   type(swap_state_t) :: state
   integer :: i

   call setup_minimal_surface_water_globals()  ! existing helper
   call surfacewater_init(state)

   ! Assert via the typed state.
   do i = 1, 4
      @assertEqual(0.0_real64, state%surfacewater%wlsbak(i), 1.0e-12_real64)
   end do

   ! Assert via the legacy global to confirm dual-write.
   do i = 1, 4
      @assertEqual(0.0_real64, wlsbak(i), 1.0e-12_real64)
   end do
end subroutine
```

- [ ] **Step 5: Build and run**

Run: `pixi run check-full`
Expected: 5 passed, 0 failed.

Run: `pixi run -e test test-pfunit 2>&1 | tail -10`
Expected: zero failures.

- [ ] **Step 6: Commit**

```bash
git add src/drainage/surfacewater_init.f90 src/core/swap.f90 tests/unit/drainage/test_surfacewater_init.pf
git commit -m "$(cat <<'EOF'
refactor(state): SS-SWST Task 3 — surfacewater_init dual-writes state and globals

surfacewater_init signature gains state arg. Body populates both
state%surfacewater%* and the legacy globals (transitional
dual-write). Per-level arrays allocated from drainage-config
dimensions. Tests assert against both targets to guard against
drift during the transitional phase.

Spec: docs/superpowers/specs/2026-05-09-state-migration-surfacewater-design.md
EOF
)"
```

---

### Task 4: Refactor `SurfaceWater(task)` entry point

Migrate the main per-step entry point. Three signature changes: add `state`, add `config` (the typed surface-water config — already exists), add `intent(out) :: request_smaller_dt` to replace the `fldecdt` global. ASSOCIATE inside `WLEVBAL` and `WBALLEV`. Dual-write state and globals.

**Files:**
- Modify: `src/drainage/surfacewater.f90` (681 LoC; main work)
- Modify: `src/core/swap.f90` — caller declares and threads `state`; reads `request_smaller_dt`
- Modify: `src/core/initialize.f90` — `state` initialization (zero-init at top of run)

- [ ] **Step 1: Update the signature of `SurfaceWater`**

Current (read the file to confirm):

```fortran
subroutine SurfaceWater(task)
   integer, intent(in) :: task
```

New:

```fortran
subroutine SurfaceWater(task, state, request_smaller_dt)
   use swap_state_mod, only: swap_state_t
   integer,            intent(in)    :: task
   type(swap_state_t), intent(inout) :: state
   logical,            intent(out)   :: request_smaller_dt
```

`request_smaller_dt` is set to `.false.` at the top of the routine; set to `.true.` at the location(s) that today write `fldecdt = .true.`.

- [ ] **Step 2: Wrap the body in ASSOCIATE blocks**

Inside `SurfaceWater(task)` and inside the helpers it calls (`WLEVBAL`, `WBALLEV`), add ASSOCIATE blocks at appropriate scopes. Example for `WLEVBAL`:

```fortran
subroutine WLEVBAL(state, request_smaller_dt)
   use swap_state_mod, only: swap_state_t
   type(swap_state_t), intent(inout) :: state
   logical,            intent(inout) :: request_smaller_dt
   ! ... (other locals)

   associate(wls    => state%surfacewater%wls,    &
             wlstar => state%surfacewater%wlstar, &
             swst   => state%surfacewater%swst,   &
             cqdrd  => state%surfacewater%cqdrd,  &
             cwsupp => state%surfacewater%cwsupp, &
             cwout  => state%surfacewater%cwout,  &
             ! …all 29 surface-water-owned scalars + arrays referenced in this routine
             vtair  => state%surfacewater%vtair)

      ! body unchanged — variable names match
      wls = wls + ...
      ! ...

      ! when fldecdt would be set:
      if (oscillation_detected) then
         request_smaller_dt = .true.
         fldecdt = .true.   ! transitional dual-write — still set the global
      end if
   end associate

   ! after the associate, also dual-write the legacy globals from state for
   ! any borrowed-by-output variables (ensures output that still reads
   ! globals sees consistent values until Task 5 migrates output)
   wls    = state%surfacewater%wls
   wlstar = state%surfacewater%wlstar
   ! ... 29 dual-writes total

end subroutine WLEVBAL
```

The dual-write pattern: the ASSOCIATE writes to state directly (since the binding aliases the field). Then after the ASSOCIATE, copy state-fields back to legacy globals so output (still reading globals) sees the same values. This dual-write is removed in Task 6.

- [ ] **Step 3: Update the caller in `src/core/swap.f90`**

Find the existing calls to `SurfaceWater(...)` (there are typically three: task=1, task=2, task=3). Change each to:

```fortran
call SurfaceWater(1, state, request_smaller_dt)
call SurfaceWater(2, state, request_smaller_dt)
call SurfaceWater(3, state, request_smaller_dt)
```

Declare locals at the top of `swap_main`:

```fortran
type(swap_state_t) :: state
logical            :: request_smaller_dt
```

For Phase 1, the `request_smaller_dt` flag is informational. **Do NOT migrate `timecontrol.f90`'s readers of `fldecdt`** — that's Phase 2's work. The dual-write inside `WLEVBAL` keeps `fldecdt = .true.` set when `request_smaller_dt = .true.`, so timecontrol's existing readers (`src/core/timecontrol.f90:62, 584`) continue to work unchanged. swap_main may receive `request_smaller_dt` and ignore it for now; the variable's existence at the call site is the architectural payload.

In Phase 2, `fldecdt` is removed from variables.f90 and timecontrol.f90 takes `request_smaller_dt` as a parameter from swap_main.

- [ ] **Step 4: Update `initialize.f90` to zero-init `state` at start of run**

Add `state` argument to `Initialize` (or the equivalent run-init routine) and zero-init the typed state. Since Fortran derived types with default-initializers handle this automatically, the only required action may be to declare the variable at top scope and let defaults apply. Confirm by inspection.

- [ ] **Step 5: Build and run**

Run: `pixi run check-full`
Expected: 5 passed, 0 failed.

Run: `pixi run -e test test-pfunit 2>&1 | tail -10`
Expected: zero failures.

If check-full fails on any case, the dual-write is incomplete — find the missing global write inside the ASSOCIATE block and add it to the post-block dual-write list.

- [ ] **Step 6: Commit**

```bash
git add src/drainage/surfacewater.f90 src/core/swap.f90 src/core/initialize.f90
git commit -m "$(cat <<'EOF'
refactor(state): SS-SWST Task 4 — SurfaceWater(task) takes typed state

SurfaceWater signature gains state and intent(out) request_smaller_dt
(replacing the legacy fldecdt global write). ASSOCIATE blocks bind
state%surfacewater%* fields to local names matching the legacy
globals. Body unchanged. Dual-write of legacy globals preserved at
the end of each ASSOCIATE — output still reads globals until Task 5.

Spec: docs/superpowers/specs/2026-05-09-state-migration-surfacewater-design.md
EOF
)"
```

---

### Task 5: Refactor `runoff()` and migrate output reads

Two related changes that ship together because output reads call `runoff()` via the surfacewater utils.

**Files:**
- Modify: `src/utils/surfacewaterutils.f90` — `runoff()` signature gains `state`
- Modify: `src/boundary/boundtop.f90` — caller of `runoff()` passes `state` (which propagates from `swap_main` via the existing call chain)
- Modify: `src/io/swap_csv_output.f90` — `set_values` reads surface-water fields from `state%surfacewater`; `csv_out` signature gains `state`
- Modify: `src/io/swapoutput.f90` — `outdrf`, `outswb`, `outend`, dispatcher signatures gain `state`; `swstini` year-reset moves to a callback
- Modify: `src/core/swap.f90` — pass `state` to all output routines
- Modify: `tests/unit/io/toml/test_surfacewater_parity.pf` — migrate to assert on state

- [ ] **Step 1: Refactor `runoff()`**

Read `src/utils/surfacewaterutils.f90`. Locate `runoff()`. Add `state` argument:

```fortran
subroutine runoff(state, ...other-existing-args, runoff_value)
   use swap_state_mod, only: swap_state_t
   type(swap_state_t), intent(in) :: state
   ! ...
   associate(wls => state%surfacewater%wls, swst => state%surfacewater%swst)
      ! body unchanged
   end associate
end subroutine
```

Update `boundtop.f90`'s call site to pass `state`. `boundtop.f90` itself takes `state` — propagate from its caller (likely a per-step routine in `swap.f90`).

- [ ] **Step 2: Migrate `set_values` in `swap_csv_output.f90`**

`set_values` currently reads surface-water-owned globals to populate the internal `vars%value` table. Switch each to read from `state%surfacewater`:

Before:
```fortran
if (vars%name(i) == 'WLS')   vars%value(1,i) = wls
```

After:
```fortran
if (vars%name(i) == 'WLS')   vars%value(1,i) = state%surfacewater%wls
```

Add `state` argument to `csv_out` and `set_values` (and any private helpers they call that need to see surfacewater state). Update the `use variables, only: …` clause of `swap_csv_output.f90` to drop the surface-water-owned symbols (`wls`, `swst`, `qdra`, etc. — see spec D2 for the full list).

`swap_csv_output.f90` retains `use variables, only: …` for borrowed reads (water-balance, crop, etc., still on globals); only the surface-water-owned reads change.

- [ ] **Step 3: Migrate `outdrf`, `outswb`, `outend` in `swapoutput.f90`**

For each routine that reads surface-water-owned globals, add `state` as an argument and switch reads to `state%surfacewater%*`. Same `use variables, only: …` cleanup pattern as Step 2.

For `outswb`'s year-reset of `swstini` (the hazard-4 mutation): move it out of `outswb` into a dedicated routine in `surfacewater.f90`:

```fortran
subroutine surfacewater_year_reset(sw)
   type(surfacewater_state_t), intent(inout) :: sw
   sw%swstini = sw%swst
end subroutine
```

Make it `public` from `surfacewater_mod`. Call from `swap_main` at the year boundary (find the existing year-boundary in `swap.f90` near where `outswb` runs).

- [ ] **Step 4: Update the dispatcher in `swapoutput.f90`**

`subroutine swapoutput(task)` becomes `subroutine swapoutput(task, state)`. Update its caller in `swap.f90`. Internal calls to migrated subroutines pass `state`.

- [ ] **Step 5: Migrate `tests/unit/io/toml/test_surfacewater_parity.pf`**

Assertions that compare against legacy globals are a parity check — they verify the dual-write keeps state and globals in sync. This task adds new state-side assertions; legacy-global assertions stay as parity guards.

- [ ] **Step 6: Build and run**

Run: `pixi run check-full`
Expected: 5 passed, 0 failed.

Run: `pixi run -e test test-pfunit 2>&1 | tail -10`
Expected: zero failures.

- [ ] **Step 7: Commit**

```bash
git add src/utils/surfacewaterutils.f90 src/boundary/boundtop.f90 \
        src/io/swap_csv_output.f90 src/io/swapoutput.f90 \
        src/drainage/surfacewater.f90 src/core/swap.f90 \
        tests/unit/io/toml/test_surfacewater_parity.pf
git commit -m "$(cat <<'EOF'
refactor(state): SS-SWST Task 5 — runoff() and output read from typed state

runoff() in surfacewaterutils takes state arg. csv_out, swapoutput,
outdrf, outswb, outend all gain state args and read surface-water
fields from state%surfacewater instead of legacy globals. outswb's
year-reset of swstini moves to surfacewater_year_reset callback,
invoked from swap_main at the year boundary. Output layer is now
read-only against the typed state.

Compute still dual-writes globals via post-ASSOCIATE block (Task 4).
Task 6 drops the dual-write.

Spec: docs/superpowers/specs/2026-05-09-state-migration-surfacewater-design.md
EOF
)"
```

---

### Task 6: Drop the dual-write

After Task 5, all readers of surface-water-owned globals are migrated (output uses state). Compute can stop writing globals. This is the test that proves the migration is complete: drop the global writes, run check-full, assert byte-identical.

**Files:**
- Modify: `src/drainage/surfacewater.f90` — remove the post-ASSOCIATE dual-write block; remove now-unused entries from `use variables, only:`
- Modify: `src/drainage/surfacewater_init.f90` — remove dual-write of globals; remove now-unused entries from `use variables, only:`
- Modify: `src/utils/surfacewaterutils.f90` — drop legacy-global reads from `runoff()` (it now reads from state only)
- Modify: `tests/unit/drainage/test_surfacewater_init.pf`, `tests/unit/io/toml/test_surfacewater_parity.pf` — drop the parity assertions against legacy globals (they're now stale anyway)

- [ ] **Step 1: Remove the post-ASSOCIATE dual-write block from `WLEVBAL`/`WBALLEV`**

Find the block of `wls = state%surfacewater%wls` style copies after each ASSOCIATE in `surfacewater.f90`. Delete.

- [ ] **Step 2: Remove dual-write from `surfacewater_init.f90`**

Find the legacy-global writes that paralleled `state%surfacewater%*` writes. Delete.

- [ ] **Step 3: Trim the `use variables, only:` clauses**

In `surfacewater.f90`, `surfacewater_init.f90`, `surfacewaterutils.f90`: review each `use variables, only: …` clause. Remove every symbol that's now in `state%surfacewater`. Keep symbols for *borrowed* globals (the 54 the discovery doc listed) — those subsystems haven't migrated yet.

- [ ] **Step 4: Drop legacy-global parity assertions from migrated tests**

In `test_surfacewater_init.pf` and `test_surfacewater_parity.pf`, remove the assertions that read `wls`, `swst`, etc. from the global module. They're stale because compute no longer writes them. Keep assertions against `state%surfacewater%*`.

The `_parity` test was specifically a parity check between the legacy reader and the typed config; that parity is no longer the right test target. Rename the test if helpful, or repurpose it to assert state values directly.

- [ ] **Step 5: Build and run**

Run: `pixi run check-full`
Expected: 5 passed, 0 failed. **This is the proof the migration works** — output now flows entirely through typed state, and the regression cases byte-match their baselines.

If a case fails, a global write was dropped that an unmigrated reader still expects. Find the unmigrated reader (`grep -rn "<global-name>" src/`) and either migrate it or restore the dual-write for that one variable.

Run: `pixi run -e test test-pfunit 2>&1 | tail -10`
Expected: zero failures.

- [ ] **Step 6: Commit**

```bash
git add src/drainage/surfacewater.f90 src/drainage/surfacewater_init.f90 \
        src/utils/surfacewaterutils.f90 \
        tests/unit/drainage/test_surfacewater_init.pf \
        tests/unit/io/toml/test_surfacewater_parity.pf
git commit -m "$(cat <<'EOF'
refactor(state): SS-SWST Task 6 — drop dual-write; state is authoritative

Compute no longer writes legacy globals for surface-water-owned
variables. After this commit, the 29 owned globals in variables.f90
are unread and unwritten — Phase 2 deletes them. check-full
byte-identical confirms the migration is functionally complete:
output flows entirely through typed state.

Trims `use variables, only:` clauses to drop owned symbols. Borrowed
symbols stay until their owning subsystems migrate.

Spec: docs/superpowers/specs/2026-05-09-state-migration-surfacewater-design.md
EOF
)"
```

---

### Task 7: Phase 1 verification + ADR draft

Final verification at the phase boundary. Author the ADR but don't merge to development — Phase 2 follows on the same branch.

**Files:**
- Create: `docs/adr/0030-state-migration-surfacewater-pilot.md` (draft only — finalize at end of Phase 2)

- [ ] **Step 1: Run the full verification suite**

```bash
pixi run -e test test-pfunit 2>&1 | tail -10
pixi run check-full
```

Expected:
- pFUnit: zero failures, test count = baseline + 7 (5 from Task 1 + 2 from Task 2).
- check-full: `Results: 5 passed, 0 failed`.

- [ ] **Step 2: Confirm globals are unread/unwritten**

```bash
# Each owned global should appear ZERO times in src/ outside variables.f90 declarations:
for v in wls wlstar swst swstini hwlman vtair wlsold cqdrd cwsupp cwout cqdra ZDraBas iqdra qdrtot overfl flInitDraBas imper numadj wlsbak sttab cqdrain cqdrainin cqdrainout qdra inqdra inqdra_in inqdra_out fldecdt; do
  count=$(grep -rE "\b$v\b" src/ --include="*.f90" 2>/dev/null | grep -v "src/core/variables.f90" | grep -v "src/state/" | wc -l)
  echo "$count  $v"
done | sort -rn | head -30
```

Expected: every count is 0 (or very small — only references inside `! comment` lines that survived from before). Investigate any `> 0` count that's a real read or write. ZDraBas, flInitDraBas, qdra, inqdra* may legitimately remain referenced as drainage-shared array names — flag and investigate.

- [ ] **Step 3: Author ADR 0030 (draft)**

Create `docs/adr/0030-state-migration-surfacewater-pilot.md`. Don't commit yet — wait for Phase 2 to land so the consequences section can include the global-deletion stats.

Use the structure of ADR 0028 as the model. Sections: Context, Decision, Consequences, References. Reference the spec, the discovery doc, and ADRs 0024 (dtutil de-shim future direction context) and 0029 (when it lands).

- [ ] **Step 4: Commit Phase 1 closure**

No commit needed if no files changed in this task. The verification result is the artifact. Move on to Phase 2 (separate plan).

---

## Self-Review Notes

- **Spec coverage:** D1 (scope) → all tasks. D2 (state shape) → Task 1. D3 (aggregator) → Task 2. D4 (threading + ASSOCIATE) → Tasks 3+4+5. D5 (intent(out) for fldecdt) → Task 4. D9 (output coupling) → Task 5. Phase 1 phasing → all tasks.
- **D6, D7, D8 are Phase 2.** Not in this plan.
- **Dual-write is the integration safety net.** check-full byte-identical at every task is the gate. If a task's check-full fails, fix the dual-write before committing.
- **Borrowed globals stay.** `use variables, only:` clauses retain the 54 borrowed symbols throughout this phase. Other subsystems own those.
- **No physics changes.** Bodies of `WLEVBAL`, `WBALLEV`, `runoff` look identical except for the ASSOCIATE wrapper. If the diff shows substantive arithmetic changes, something has gone wrong.
