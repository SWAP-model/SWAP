# Drainage State Migration — Phase 1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Define `drainage_state_t`, aggregate it under `swap_state_t`, move `qdra`/`qdrain` out of `surfacewater_state_t` into `drainage_state_t` (the cleanup of ADR 0030's acknowledged technical debt), modernize `divdra` to assumed-shape arrays, and migrate drainage compute + output to read/write the typed state with a transitional dual-write to legacy globals.

**Architecture:** New `src/state/drainage_state.f90` module. `drainage_state_t` field added to `swap_state_t`. `qdra`/`qdrain` re-classified from surfacewater to drainage state. `divdra` signature changed from explicit-shape `(Madr, macp)` to assumed-shape `(:,:)` so the global `qdra` can eventually be deleted. Drainage's per-step compute (`Drainage`, `bocodre`, `divdra`) and per-day output dual-writes state and globals during transition; Phase 2 drops the dual-write and removes the globals.

**Tech Stack:** Fortran 2008, gfortran, meson + ninja + pixi, pFUnit, check-full byte-identical regression.

**Spec:** `docs/superpowers/specs/2026-05-10-state-migration-drainage-design.md`

**Branch:** all commits go on `refactor/surfacewater-state` (the umbrella migration branch). Continues from Phase 2 of the surfacewater pilot.

---

## Lessons-learned applied from surfacewater pilot

- **Cross-subsystem reader inventory upfront** — drainage discovery doc already populated Section 3.5 (6 external reader files) so this plan's scope is locked.
- **Dual-write transitional pattern** — drainage compute writes both state and globals during Phase 1; Phase 2 drops globals.
- **`use Variables` shadow** — drainage subroutines use `use variables, only:` clauses already (per discovery); ASSOCIATE without `dr_*` prefix should work, but verify per subroutine.
- **OutputModflow consumer** — drainage is called from OutputModflow's perturbation loop (discovery hazard #5); state_om scope inherits the new `state%drainage` field automatically (deep-copy via Fortran assignment).
- **Verify subagent commits** — always grep for actual commit SHA + working-directory check before declaring done. The surfacewater pilot lost two agent reports to wrong-project commits.

---

## File Structure

**Created:**
- `src/state/drainage_state.f90` — `drainage_state_t` definition (6 fields per spec D2)
- `tests/unit/state/test_drainage_state.pf` — pFUnit construction & defaults tests

**Modified:**
- `src/state/surfacewater_state.f90` — remove `qdra` and `qdrain` declarations
- `src/state/swap_state.f90` — add `type(drainage_state_t) :: drainage` field
- `src/drainage/drainage.f90` — `Drainage` writes `state%drainage%*` alongside legacy globals; allocate per-level arrays in init path
- `src/drainage/divdra.f90` — modernize signature to assumed-shape `qdrain(:)`, `qdra(:,:)`
- `src/drainage/surfacewater.f90`, `src/drainage/surfacewater_init.f90` — every reference to `state%surfacewater%qdra` / `state%surfacewater%qdrain` → `state%drainage%qdra` / `state%drainage%qdrain`. Stop allocating these in surfacewater_init; allocate in drainage_init instead.
- `src/heat/frozencond.f90` — `state%surfacewater%qdra` / `qdrain` references → `state%drainage%`
- `src/soil/waterbalance.f90` — same rename across `integral` and `fluxes`
- `src/soil/soilhydraulics.f90` — same rename for `qdra` reads in `headcalc`
- `src/io/swap_csv_output.f90` — `state%surfacewater%iqdra` references stay (iqdra is surfacewater-owned), but any `qdra` references rename
- `src/io/swapoutput.f90` — same rename pattern; output routines that read drainage-owned fields (`outdrf`, `outbal`, `outblc`, `outafo`, `outaun`, `outend`) take `state` and read from `state%drainage`
- `src/core/swap.f90` — Drainage call already passes state; verify
- `tests/unit/testSuites.inc` — add `ADD_TEST_SUITE(test_drainage_state_suite)`
- meson source list — add `src/state/drainage_state.f90`

**No changes to** (in this phase):
- `src/core/variables.f90` — globals stay declared until Phase 2
- `src/io/toml/config_to_variables.f90` — `geofac` fix is a Phase 2 task
- The `swdislay==2` `ztopdislay` rule in surfacewater (Phase 2 owner-relocation)

---

### Task 1: Create `drainage_state_t`

Define the typed state module. Pure data type with default-initialized scalar and unallocated arrays. Per-level arrays are allocated by drainage_init (or wherever the init lifecycle goes — Task 5 decides).

**Files:**
- Create: `src/state/drainage_state.f90`
- Create: `tests/unit/state/test_drainage_state.pf`
- Modify: `tests/unit/testSuites.inc`
- Modify: meson source list

- [ ] **Step 1: Locate the meson source list**

```bash
grep -rln "src/state/surfacewater_state.f90" --include="meson.build" . 2>/dev/null
```

Note the paths — the new module goes alongside surfacewater_state.f90 in the same lists.

- [ ] **Step 2: Write the failing pFUnit tests**

Create `tests/unit/state/test_drainage_state.pf`:

```fortran
! Tests for drainage_state_t — typed state record for the drainage
! subsystem. Phase 1 of the drainage state-migration arc.
!
! See docs/superpowers/specs/2026-05-10-state-migration-drainage-design.md.

@test
subroutine test_drainage_state_default_scalars()
   use funit
   use iso_fortran_env, only: real64
   use drainage_state_mod, only: drainage_state_t
   type(drainage_state_t) :: dr

   @assertEqual(0.0_real64, dr%qdrd, 1.0e-12_real64)
end subroutine test_drainage_state_default_scalars

@test
subroutine test_drainage_state_arrays_unallocated()
   use funit
   use drainage_state_mod, only: drainage_state_t
   type(drainage_state_t) :: dr

   @assertFalse(allocated(dr%qdrain))
   @assertFalse(allocated(dr%drainl))
   @assertFalse(allocated(dr%wetper))
   @assertFalse(allocated(dr%ztopdislay))
   @assertFalse(allocated(dr%qdra))
end subroutine test_drainage_state_arrays_unallocated

@test
subroutine test_drainage_state_array_allocation()
   ! Caller (drainage_init) allocates per-level arrays. Verify the
   ! allocation pattern works.
   use funit
   use iso_fortran_env, only: real64
   use drainage_state_mod, only: drainage_state_t
   type(drainage_state_t) :: dr
   integer, parameter :: Madr = 5
   integer, parameter :: macp = 100

   allocate(dr%qdrain(Madr));         dr%qdrain     = 0.0_real64
   allocate(dr%drainl(Madr));         dr%drainl     = 0.0_real64
   allocate(dr%wetper(Madr));         dr%wetper     = 0.0_real64
   allocate(dr%ztopdislay(Madr));     dr%ztopdislay = 0.0_real64
   allocate(dr%qdra(Madr, macp));     dr%qdra       = 0.0_real64

   @assertTrue(allocated(dr%qdrain))
   @assertEqual(Madr, size(dr%qdrain))
   @assertEqual(Madr, size(dr%qdra, 1))
   @assertEqual(macp, size(dr%qdra, 2))
end subroutine test_drainage_state_array_allocation

@test
subroutine test_drainage_state_independent_instances()
   ! Two instances must be independent (no shared module storage).
   use funit
   use iso_fortran_env, only: real64
   use drainage_state_mod, only: drainage_state_t
   type(drainage_state_t) :: dr_a, dr_b

   dr_a%qdrd = 1.5_real64
   dr_b%qdrd = -2.5_real64

   @assertEqual(1.5_real64,  dr_a%qdrd, 1.0e-12_real64)
   @assertEqual(-2.5_real64, dr_b%qdrd, 1.0e-12_real64)
end subroutine test_drainage_state_independent_instances
```

Add `ADD_TEST_SUITE(test_drainage_state_suite)` to `tests/unit/testSuites.inc` near `test_surfacewater_state_suite`.

- [ ] **Step 3: Run tests — confirm FAIL**

Run: `pixi run -e test test-pfunit 2>&1 | tail -10`
Expected: 4 failures referencing missing module `drainage_state_mod`.

- [ ] **Step 4: Implement the module**

Create `src/state/drainage_state.f90`:

```fortran
!> @file drainage_state.f90
!! Typed state record for the drainage subsystem. Owns the 6
!! drainage-flux variables that legacy SWAP held in the
!! `variables.f90` globals module.
!!
!! `qdra` and `qdrain` were temporarily classified into
!! surfacewater_state_t during the surface-water migration pilot
!! (ADR 0030) because of write-site overlap. The drainage migration
!! arc moves them here, where they architecturally belong.
!!
!! See ADR 0031 (state-migration drainage subsystem) and the
!! 2026-05-10 design spec.

module drainage_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   public :: drainage_state_t

   type :: drainage_state_t
      ! Scalar fluxes
      real(real64) :: qdrd     = 0.0_real64    ! drain-direction sub-flux (cm/d)

      ! Per-level (Madr-sized) arrays — allocated by drainage_init
      real(real64), allocatable :: qdrain(:)       ! lateral drainage flux per level (cm/d)
      real(real64), allocatable :: drainl(:)       ! drain length per level (cm)
      real(real64), allocatable :: wetper(:)       ! wetted perimeter per level (cm)
      real(real64), allocatable :: ztopdislay(:)   ! top of discharge layer per level (cm)

      ! Per-level / per-compartment array (Madr × macp)
      real(real64), allocatable :: qdra(:,:)       ! lateral drainage flux per level/compartment (cm/d)
   end type drainage_state_t

end module drainage_state_mod
```

Add `src/state/drainage_state.f90` to the meson source list.

- [ ] **Step 5: Run tests — confirm PASS**

Run: `pixi run -e test test-pfunit 2>&1 | tail -10`
Expected: zero failures. Test count = previous baseline + 4.

- [ ] **Step 6: Commit**

```bash
git add src/state/drainage_state.f90 tests/unit/state/test_drainage_state.pf tests/unit/testSuites.inc
# Plus the meson.build files modified
git status
git add <the meson files>
git commit -m "$(cat <<'EOF'
feat(state): SS-DRST Task 1 — drainage_state_t typed state record

Phase 1 foundation. Defines the drainage subsystem's typed state
record with 6 fields (qdrd, qdrain, drainl, wetper, ztopdislay,
qdra). Per-level arrays are allocatable; drainage_init will
allocate them once drainage-config dimensions (Madr, macp) are
visible.

qdra and qdrain previously lived in surfacewater_state_t (per ADR
0030's acknowledged technical debt). Task 3 moves all consumer
references and removes them from surfacewater_state_t.

Spec: docs/superpowers/specs/2026-05-10-state-migration-drainage-design.md
EOF
)"
```

---

### Task 2: Add `drainage` to `swap_state_t`

**Files:**
- Modify: `src/state/swap_state.f90`
- Modify: `tests/unit/state/test_swap_state.pf` — add an assertion that drainage field exists

- [ ] **Step 1: Add the new test**

Append to `tests/unit/state/test_swap_state.pf`:

```fortran
@test
subroutine test_swap_state_has_drainage()
   use funit
   use iso_fortran_env, only: real64
   use swap_state_mod, only: swap_state_t
   type(swap_state_t) :: state

   @assertEqual(0.0_real64, state%drainage%qdrd, 1.0e-12_real64)
   @assertFalse(allocated(state%drainage%qdrain))
   @assertFalse(allocated(state%drainage%qdra))
end subroutine test_swap_state_has_drainage
```

- [ ] **Step 2: Run tests — confirm FAIL**

Run: `pixi run -e test test-pfunit 2>&1 | tail -10`
Expected: 1 new failure (state%drainage doesn't exist).

- [ ] **Step 3: Update swap_state_t**

In `src/state/swap_state.f90`:

```fortran
module swap_state_mod
   use surfacewater_state_mod, only: surfacewater_state_t
   use drainage_state_mod,     only: drainage_state_t
   implicit none
   private
   public :: swap_state_t

   type :: swap_state_t
      type(surfacewater_state_t) :: surfacewater
      type(drainage_state_t)     :: drainage
   end type swap_state_t

end module swap_state_mod
```

- [ ] **Step 4: Run tests — confirm PASS**

Run: `pixi run -e test test-pfunit 2>&1 | tail -10`
Expected: zero failures.

- [ ] **Step 5: Commit**

```bash
git add src/state/swap_state.f90 tests/unit/state/test_swap_state.pf
git commit -m "$(cat <<'EOF'
feat(state): SS-DRST Task 2 — swap_state_t gains drainage field

swap_state_t now aggregates surfacewater_state_t + drainage_state_t.
Subsequent migrations grow the type one field at a time.

Spec: docs/superpowers/specs/2026-05-10-state-migration-drainage-design.md
EOF
)"
```

---

### Task 3: Move `qdra` and `qdrain` from `surfacewater_state_t` to `drainage_state_t`

The cleanup of ADR 0030's acknowledged technical debt. Mechanical refactor: every `state%surfacewater%qdra` becomes `state%drainage%qdra`; same for `qdrain`. Remove both from `surfacewater_state_t`. Move array allocation from `surfacewater_init` to `drainage_init` (or similar).

**Files (touched by the refactor):**
- Modify: `src/state/surfacewater_state.f90` — remove `qdra` and `qdrain` declarations
- Modify: `src/drainage/surfacewater_init.f90` — drop allocation of these arrays; drop dual-writes
- Modify: `src/drainage/surfacewater.f90` — every `state%surfacewater%qdra` → `state%drainage%qdra`; same for `qdrain`
- Modify: `src/drainage/drainage.f90` — same rename
- Modify: `src/drainage/divdra.f90` — same rename (if it references state directly; otherwise via callers)
- Modify: `src/heat/frozencond.f90` — same rename (FrozenBounds writes/reads)
- Modify: `src/soil/waterbalance.f90` — same rename in `integral` and `fluxes`
- Modify: `src/soil/soilhydraulics.f90` — same rename in `headcalc`
- Modify: `src/io/swapoutput.f90` — same rename across output routines
- Modify: `src/io/swap_csv_output.f90` — same rename in `set_values`
- Modify: tests that reference `state%surfacewater%qdra` or `qdrain`

We need to **also** initialize allocation of these arrays via a new `drainage_init(state)` entry point — see step 4 below.

- [ ] **Step 1: Find every reference**

```bash
grep -rEn "state%surfacewater%(qdra|qdrain)\b" src/ tests/ --include="*.f90" --include="*.pf"
```

Tally the file/line list. Expect ~30-50 hits across compute and output.

- [ ] **Step 2: Decide where allocation moves**

Read `src/drainage/surfacewater_init.f90` and find the allocate blocks for `qdra` and `qdrain` in `state%surfacewater%*`. Two choices:

- **(a)** Add allocation to a new `drainage_init(state)` subroutine in `src/drainage/drainage.f90` (or a new `src/drainage/drainage_init.f90`). Call from `swap_main` once at startup, before any drainage compute.
- **(b)** Allocate inline in `Drainage(task=1)` if such a task branch exists, or at the top of the first `Drainage()` invocation with an `if (.not. allocated(...))` guard.

Recommend **(a)** for symmetry with `surfacewater_init`. Create `drainage_init(state)`. Allocate the per-level arrays from drainage-config dimensions (`nrlevs`, `numnod`).

- [ ] **Step 3: Add drainage_init**

In `src/drainage/drainage.f90` (or a new file `src/drainage/drainage_init.f90`):

```fortran
subroutine drainage_init(state)
   use swap_state_mod, only: swap_state_t
   use variables, only: nrlevs, numnod
   type(swap_state_t), intent(inout) :: state

   if (.not. allocated(state%drainage%qdrain))     allocate(state%drainage%qdrain(nrlevs))
   if (.not. allocated(state%drainage%drainl))     allocate(state%drainage%drainl(nrlevs))
   if (.not. allocated(state%drainage%wetper))     allocate(state%drainage%wetper(nrlevs))
   if (.not. allocated(state%drainage%ztopdislay)) allocate(state%drainage%ztopdislay(nrlevs))
   if (.not. allocated(state%drainage%qdra))       allocate(state%drainage%qdra(nrlevs, numnod))

   state%drainage%qdrain     = 0.0_real64
   state%drainage%drainl     = 0.0_real64
   state%drainage%wetper     = 0.0_real64
   state%drainage%ztopdislay = 0.0_real64
   state%drainage%qdra       = 0.0_real64
end subroutine drainage_init
```

If you put it in a new file, add it to `src/drainage/drainage_init.f90` and update meson.

Make `drainage_init` `public` from the relevant module.

- [ ] **Step 4: Wire drainage_init from swap_main**

In `src/core/swap.f90`, find where `surfacewater_init` is called (Task 1 entry of `SurfaceWater`). Call `drainage_init(state)` BEFORE `surfacewater_init` — ensures drainage state arrays are ready before any drainage compute.

```bash
grep -n "surfacewater_init\|SurfaceWater(1" src/core/swap.f90 | head
```

Add the call at the appropriate point (likely right before or right after the `flSurfaceWater`-gated `SurfaceWater(1)` call). drainage_init can run unconditionally — its arrays are sized off drainage config, which is always present.

- [ ] **Step 5: Remove qdra and qdrain allocation from surfacewater_init**

In `src/drainage/surfacewater_init.f90`, find the allocate blocks for `state%surfacewater%qdra` and `state%surfacewater%qdrain`. Delete them (drainage_init now handles).

Also delete dual-writes of these to legacy globals (now drainage's responsibility — Task 5 wires this).

- [ ] **Step 6: Remove qdra and qdrain from surfacewater_state_t**

In `src/state/surfacewater_state.f90`, remove the field declarations:

```fortran
! REMOVE:
real(real64), allocatable :: qdra(:,:)
! and remove qdrain if it was there
```

- [ ] **Step 7: Refactor every consumer reference**

Use a sed-style mechanical rename across all the files in the inventory from Step 1:

`state%surfacewater%qdra` → `state%drainage%qdra`
`state%surfacewater%qdrain` → `state%drainage%qdrain`

Be careful with `qdrain` — the legacy global `qdrain` (no `state%` prefix) stays untouched. Only the typed-state references rename.

- [ ] **Step 8: Update tests**

`tests/unit/state/test_surfacewater_state.pf` — drop assertions on `qdra` (no longer in `surfacewater_state_t`); the corresponding assertions live in `test_drainage_state.pf` (Task 1).

- [ ] **Step 9: Build and verify**

Run: `pixi run check-full`
Expected: `Results: 5 passed, 0 failed`. The rename is mechanical — same fields, different state slice.

Run: `pixi run -e test test-pfunit 2>&1 | tail -10`
Expected: zero failures.

- [ ] **Step 10: Commit**

```bash
git add -A   # be explicit about which files; this is a wide rename
git status
# verify only the expected files are staged
git commit -m "$(cat <<'EOF'
refactor(state): SS-DRST Task 3 — move qdra/qdrain to drainage_state_t

Cleanup of ADR 0030's acknowledged technical debt. qdra and qdrain
were temporarily in surfacewater_state_t because of write-site
overlap during surface-water Phase 1; this task moves them to
drainage_state_t where they architecturally belong.

Mechanical rename: state%surfacewater%qdra → state%drainage%qdra,
same for qdrain. Allocation moves from surfacewater_init to a new
drainage_init(state) entry point called from swap_main. Both fields
are removed from surfacewater_state_t.

check-full byte-identical confirms no semantic change.

Spec: docs/superpowers/specs/2026-05-10-state-migration-drainage-design.md
EOF
)"
```

---

### Task 4: Modernize `divdra` to assumed-shape

`divdra` currently takes `qdrain(Madr)` and `qdra(Madr,macp)` as explicit-shape arguments. This forces the legacy globals to remain because callers must pass arrays with statically-known dimensions. Modernize to assumed-shape so callers can pass `state%drainage%qdrain` and `state%drainage%qdra` directly. Phase 2 then deletes the legacy globals.

**Files:**
- Modify: `src/drainage/divdra.f90`
- Modify: every caller of `divdra` (verify via grep)

- [ ] **Step 1: Inspect divdra signature**

```bash
grep -n "subroutine divdra\|^.*divdra" src/drainage/divdra.f90 | head
sed -n '1,30p' src/drainage/divdra.f90
```

Note the current arg list and the explicit-shape declarations.

- [ ] **Step 2: Find all callers**

```bash
grep -rEn "\bcall divdra\b|\bcall DIVDRA\b" src/ --include="*.f90"
```

Catalog each call site.

- [ ] **Step 3: Modernize the signature**

Change in `divdra.f90`:

```fortran
! BEFORE:
real(real64), intent(in)  :: qdrain(Madr)
real(real64), intent(out) :: qdra(Madr, macp)

! AFTER:
real(real64), intent(in)  :: qdrain(:)
real(real64), intent(out) :: qdra(:,:)
```

Inside the body, replace any `Madr` and `macp` literal references with `size(qdrain)` and `size(qdra, 2)`. Verify by reading the body that the dimensions are used only for bounds (loops, array operations) — not as part of arithmetic with assumed module-scope values.

If `Madr`/`macp` appear in the body for purposes other than the qdrain/qdra arrays' dimensions (e.g., other arrays sized by the same constants), keep them via `use variables, only: Madr, macp` and only switch the qdrain/qdra-related uses.

- [ ] **Step 4: Update callers**

Every `call divdra(..., qdrain, qdra, ...)` site now passes assumed-shape arrays. The call site syntax is unchanged for `qdrain` and `qdra` arguments — Fortran's assumed-shape passing is implicit.

The callers may currently pass legacy globals `qdrain` and `qdra`. Phase 1 keeps them passing globals (because dual-write keeps globals current); Phase 2 switches to `state%drainage%qdrain` / `state%drainage%qdra`. **For Task 4, just modernize the signature; don't change call sites yet.**

- [ ] **Step 5: Build and verify**

Run: `pixi run check-full`
Expected: 5/5 byte-identical. Assumed-shape is a binary-compatible signature change for explicit-shape callers passing the same array shape.

Run: `pixi run -e test test-pfunit 2>&1 | tail -10`
Expected: zero failures.

- [ ] **Step 6: Commit**

```bash
git add src/drainage/divdra.f90
git commit -m "$(cat <<'EOF'
refactor(drainage): SS-DRST Task 4 — divdra modernized to assumed-shape

divdra now takes qdrain(:) and qdra(:,:) as assumed-shape arrays
instead of explicit-shape (Madr, macp). This unblocks deletion of
the legacy qdra and qdrain globals from variables.f90 in Phase 2:
callers can pass state%drainage%qdra / state%drainage%qdrain
directly without needing the static dimensions.

Signature change only — callers continue passing legacy globals
during Phase 1. Phase 2 switches call sites to typed state and
deletes the globals.

Spec: docs/superpowers/specs/2026-05-10-state-migration-drainage-design.md
EOF
)"
```

---

### Task 5: Drainage compute dual-writes the remaining 4 fields

Drainage compute already writes `qdra` and `qdrain` to `state%drainage%*` (after Tasks 3-4). This task adds dual-write for the remaining 4 owned fields: `qdrd`, `drainl`, `wetper`, `ztopdislay`. Compute also writes them to legacy globals for backwards compat — Phase 2 drops the dual-write.

**Files:**
- Modify: `src/drainage/drainage.f90` — add `state%drainage%X = ...` writes paralleling the existing `X = ...` global writes for the 4 fields

- [ ] **Step 1: Find the global writes**

```bash
grep -nE "^\s+(qdrd|drainl|wetper|ztopdislay)\s*=" src/drainage/drainage.f90 | head -20
```

Note each line.

- [ ] **Step 2: Add dual-writes**

For each `<global> = <RHS>` line in drainage.f90, add a paired `state%drainage%<global> = <same RHS>` line immediately after.

If the writes happen inside an ASSOCIATE block, the cleaner approach is to add the field to the ASSOCIATE bind list and use the bare alias. Check the existing structure — Task 3 of surfacewater Phase 2 added ASSOCIATE blocks in some drainage subroutines.

- [ ] **Step 3: Build and verify**

Run: `pixi run check-full`
Expected: 5/5 byte-identical. Dual-write doesn't change observable behavior.

Run: `pixi run -e test test-pfunit 2>&1 | tail -10`
Expected: zero failures.

- [ ] **Step 4: Commit**

```bash
git add src/drainage/drainage.f90
git commit -m "$(cat <<'EOF'
refactor(state): SS-DRST Task 5 — drainage compute dual-writes 4 owned fields

drainage.f90 now writes qdrd, drainl, wetper, ztopdislay to both
state%drainage%* and the legacy globals. Phase 1 transitional —
output reads still come from globals until Task 6, then Phase 2
drops the dual-write and removes the globals.

Spec: docs/superpowers/specs/2026-05-10-state-migration-drainage-design.md
EOF
)"
```

---

### Task 6: Output reads migrate to `state%drainage`

For every output routine that reads drainage-owned globals, switch reads to `state%drainage%*`. Output routines mostly already take `state` (added during surfacewater Phase 2 Task 5).

**Files:**
- Modify: `src/io/swapoutput.f90` — `outdrf` (reads `qdrain`), `outbal`, `outblc` (read `qdrain`), `outafo`, `outaun` (read drainage state via `divdra` in some paths), `outend` (reads drainage state on year boundary)
- Modify: `src/io/swap_csv_output.f90` — `set_values` reads `qdrain` family

- [ ] **Step 1: Find all output reads of the 4 fields**

```bash
grep -rnE "\b(qdrd|drainl|wetper|ztopdislay)\b" src/io/swapoutput.f90 src/io/swap_csv_output.f90 | head -30
```

For each match, classify: read or write? Read site to migrate.

(Note: most output is per-step writes to file — they READ the globals to format output. We switch those reads.)

- [ ] **Step 2: Migrate each output read**

Replace `<global>` read with `state%drainage%<global>`. The output routines already have `state` in scope (Phase 2 of surfacewater plumbed it). Drop the symbols from `use variables, only: ...` clauses.

- [ ] **Step 3: Build and verify**

Run: `pixi run check-full`
Expected: 5/5 byte-identical.

Run: `pixi run -e test test-pfunit 2>&1 | tail -10`
Expected: zero failures.

- [ ] **Step 4: Commit**

```bash
git add src/io/swapoutput.f90 src/io/swap_csv_output.f90
git commit -m "$(cat <<'EOF'
refactor(state): SS-DRST Task 6 — output reads drainage from typed state

Output routines (outdrf, outbal, outblc, outafo, outaun, outend,
set_values) now read qdrd/drainl/wetper/ztopdislay from
state%drainage instead of legacy globals. The output layer is now
read-only against the typed state for drainage-owned variables.

Compute still dual-writes globals via Task 5; Task 7 drops the
dual-write.

Spec: docs/superpowers/specs/2026-05-10-state-migration-drainage-design.md
EOF
)"
```

---

### Task 7: Drop dual-write — state becomes authoritative

After Task 6, drainage globals are unread by output. Compute can stop writing them. This is the proof that Phase 1 of the drainage migration is functionally complete.

**Files:**
- Modify: `src/drainage/drainage.f90` — remove the dual-write of `qdrd`, `drainl`, `wetper`, `ztopdislay` to legacy globals
- Modify: tests — drop any legacy-global parity assertions for drainage-owned fields

- [ ] **Step 1: Find dual-write sites**

```bash
grep -nE "^\s+(qdrd|drainl|wetper|ztopdislay)\s*=\s*state%drainage%|^\s+state%drainage%(qdrd|drainl|wetper|ztopdislay)\s*=" src/drainage/drainage.f90 | head -20
```

For each pair (legacy `X = state%drainage%X` after the state write, OR pre-existing legacy write paired with state write), delete the legacy global write. Keep the `state%drainage%X = <RHS>` writes.

- [ ] **Step 2: Trim use variables**

In `drainage.f90`, drop `qdrd`, `drainl`, `wetper`, `ztopdislay` from any `use variables, only: ...` clause that was present (Phase 1 Task 5 may have added them; not all drainage subroutines may have them).

- [ ] **Step 3: Build and verify — THE INTEGRATION GATE**

Run: `pixi run check-full`
Expected: `Results: 5 passed, 0 failed`. **This is the proof Phase 1 works** — state is authoritative for the 4 fields.

Run: `pixi run -e test test-pfunit 2>&1 | tail -10`
Expected: zero failures.

If a case fails: a reader was missed in Task 6. Find the unmigrated reader (`grep -rn "use variables.*\b(qdrd|drainl|wetper|ztopdislay)\b" src/`) and migrate or restore that variable's dual-write.

Note: `qdra` and `qdrain` dual-writes (legacy globals + state) STAY through Phase 2 because of cross-subsystem readers and `divdra`'s callers passing globals. Don't touch them in this task.

- [ ] **Step 4: Commit**

```bash
git add src/drainage/drainage.f90
# plus any modified tests
git commit -m "$(cat <<'EOF'
refactor(state): SS-DRST Task 7 — drop dual-write; drainage state authoritative

After Task 6 migrated output readers, drainage compute no longer
writes qdrd/drainl/wetper/ztopdislay to legacy globals. State is
the sole authoritative location for those 4 fields; the globals
in variables.f90 are unwritten and unread by Phase 1's perimeter
of the codebase. Phase 2 deletes them.

qdra and qdrain dual-writes remain (legacy globals still consumed
by divdra callers in compute and by frozencond). Phase 2 finishes
their migration.

check-full byte-identical confirms Phase 1's drainage migration is
functionally complete.

Spec: docs/superpowers/specs/2026-05-10-state-migration-drainage-design.md
EOF
)"
```

---

### Task 8: Phase 1 verification + ADR 0031 draft

**Files:**
- Create: `docs/adr/0031-state-migration-drainage.md` (draft — finalized at end of Phase 2)

- [ ] **Step 1: Final verification**

```bash
pixi run -e test test-pfunit 2>&1 | tail -10
pixi run check-full
```

Expected: zero pFUnit failures, 5/5 check-full.

- [ ] **Step 2: Confirm 4 owned fields are unread/unwritten globally**

```bash
for v in qdrd drainl wetper ztopdislay; do
  count=$(grep -rEn "\b$v\b" src/ --include="*.f90" 2>/dev/null \
    | grep -v "src/core/variables.f90" \
    | grep -v "src/core/initialize.f90" \
    | grep -v "src/state/" \
    | wc -l)
  echo "$count  $v"
done
```

Expected: every count is 0 (or only inside ASSOCIATE bindings using the bare names — those are not direct global accesses). Investigate any nonzero count.

`qdrain` and `qdra` will have nonzero counts — they're still globally written by drainage and read by cross-subsystem consumers. That's expected for Phase 1.

- [ ] **Step 3: Draft ADR 0031**

Create `docs/adr/0031-state-migration-drainage.md`. Don't commit yet if Phase 2 will change the consequences — draft now, finalize after Phase 2 lands. Use ADR 0030 as the template.

Required sections:
- **Status:** accepted (Phase 1 complete; Phase 2 pending)
- **Context:** drainage as migration #2; qdra/qdrain reclassification rationale
- **Decision:** state-type definition, divdra modernization, dual-write transitional
- **Consequences:** Phase 1 outcomes (4 of 6 owned globals migrated cleanly; qdra/qdrain re-classified). Phase 2 consequences TBD.
- **References:** discovery doc, design spec, predecessor ADR 0030, plan files.

- [ ] **Step 4: Commit ADR draft**

```bash
git add docs/adr/0031-state-migration-drainage.md
git commit -m "docs(adr): ADR 0031 draft — drainage state migration (Phase 1 complete; Phase 2 pending)"
```

---

## Self-Review Notes

- **Spec coverage:** D1 (scope) → all tasks. D2 (state shape) → Task 1. D3 (aggregator) → Task 2. D4 (threading) → leveraging surfacewater Phase 2 work. D5 (qdra/qdrain reclassification) → Task 3. D6 (divdra modernization) → Task 4. D9 (geofac fix) is Phase 2.
- **Dual-write transitional pattern** keeps check-full byte-identical at every commit. Task 7 is the integration gate.
- **`qdra` and `qdrain` are already in the `state%surfacewater` namespace as of Phase 1 of surfacewater migration.** Task 3 is a rename to `state%drainage`. After Task 3, those fields don't exist in `surfacewater_state_t`. Care: don't accidentally break the surfacewater code that reads them (now via `state%drainage`).
- **`divdra` modernization is independent** — Task 4 stands alone (signature change), Phase 2 then makes use of it by passing typed state.
- **Phase 2 punch list** (carry forward from Phase 1):
  1. Move `swdislay==2` `ztopdislay` rule from `SurfaceWater(2)` into `Drainage()` (spec D7).
  2. Migrate cross-subsystem readers (frozencond's `qdrain` writeback, etc.) to read from `state%drainage`.
  3. Drop `qdra` and `qdrain` dual-writes (compute writes only state).
  4. Delete the 6 owned globals (qdrd, drainl, wetper, ztopdislay, qdrain, qdra) from variables.f90.
  5. Fix the `geofac` adapter ordering inline (spec D9).
  6. ADR 0031 finalize.
