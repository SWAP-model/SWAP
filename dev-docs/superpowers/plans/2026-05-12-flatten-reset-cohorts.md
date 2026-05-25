# Flatten Reset Cohorts Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the nested cohort sub-records (`<subsystem>_intermediate_t`, `<subsystem>_cumulative_t`, and surfacewater's split cumulative cohorts) in `surfacewater_state_t`, `solute_state_t`, and `soilwater_state_t` with flat fields on the parent types. Move `reset()` from cohort-type-bound to parent-type-bound procedures, named to preserve the gate/cadence semantics.

**Architecture:** Three per-subsystem commits in order solute → soilwater → surfacewater on the development branch, each fully verified (build clean + pFUnit green + check-full byte-identical + audit grep zero) before the next begins. Each commit moves the cohort fields up to the parent type and rebinds reset procedures with names spelling out the gate (`reset_intermediate`, `reset_cumulative`, `reset_cumulative_drainage`, `reset_cumulative_reservoir`, `reset_intermediate_per_day`). Then a closing ADR documents the arc and references / supersedes the cohort-pattern portion of ADR 0033.

**Tech Stack:** Fortran 2008+, meson + ninja build (`pixi run build-linux`), pFUnit unit tests (`pixi run test-pfunit`), Python regression harness (`pixi run check-full` runs both).

**Reference design:** `docs/superpowers/specs/2026-05-12-flatten-reset-cohorts-design.md`. Read before starting.

**Verification commands (used at every task):**
- Build: `pixi run build-linux`
- pFUnit: `pixi run test-pfunit`
- Regression-only: `pixi run regression`
- Full gate (build + pFUnit + regression): `pixi run check-full`
- Per-subsystem audit greps: see each task's "Audit grep zero-out" step.

**Convention reminder (from project memory):** check-full must pass before each subsystem commit — pFUnit alone misses global-default regressions.

---

## File Structure

### Files modified per subsystem

**Task 1 — solute:**
- Rewrite: `src/state/solute_state.f90`
- Rewrite paths in: `src/solute/solute.f90` (18 refs), `src/io/swap_csv_output.f90` (6 refs)
- Rewrite tests: `tests/unit/state/test_solute_state.pf`, `tests/unit/state/test_solute_cumulative.pf`, `tests/unit/state/test_swap_state.pf` (1 ref)

**Task 2 — soilwater:**
- Rewrite: `src/state/soilwater_state.f90`
- Rewrite paths in: `src/soil/waterbalance.f90` (67), `src/io/swap_csv_output.f90` (38), `src/core/variables.f90` (16), `src/io/swapoutput.f90` (13), `src/crop/cropgrowth.f90` (12), `src/core/initialize.f90` (10), `src/crop/irrigation.f90` (9), `src/soil/soilhydraulics.f90` (5), `src/soil/soilgrid.f90` (3), `src/crop/management_soil.f90` (3), `src/core/swap.f90` (3), `src/atmosphere/meteoday.f90` (1)
- Rewrite tests: `tests/unit/state/test_soilwater_state.pf` (11)

**Task 3 — surfacewater:**
- Rewrite: `src/state/surfacewater_state.f90`
- Rewrite paths in: `src/drainage/drainage.f90` (20), `src/soil/waterbalance.f90` (12), `src/drainage/surfacewater.f90` (9), `src/io/swapoutput.f90` (6), `src/io/swap_csv_output.f90` (6), `src/soil/soilgrid.f90` (2), `src/crop/management_soil.f90` (1), `src/drainage/surfacewater_init.f90` (allocator block — paths only, no allocation changes)
- Rewrite tests: `tests/unit/state/test_surfacewater_state.pf`, `tests/unit/state/test_surfacewater_cumulative.pf`

**Task 4 — ADR:**
- Create: `docs/adr/0042-flatten-reset-cohorts.md`

### Files NOT touched
- `docs/adr/0033-cumulative-reset-cohorts.md` — historical record, left intact (the new ADR references and supersedes the cohort-pattern portion).
- Allocation lifecycle blocks (`src/drainage/surfacewater_init.f90:119-141`, `src/drainage/drainage.f90:445-467`, `src/state/soilwater_state.f90:300-419`) — only the **paths within them** rewrite from `%intr%` / `%intermediate%` etc. to flat. Allocation consolidation is out of scope (separate follow-up arc).
- AgeTracer fields, `ArMpSs`, `sqrap`, `samcra`, instantaneous fluxes (`isqbot`, `isqtop`, `isqdra`), `samini = sampro` rebase, gates and activity flags — all unchanged.

---

## Task 1: Solute — flatten cohorts and rebind resets

**Files:**
- Modify: `src/state/solute_state.f90` (full rewrite — see step 1)
- Modify: `src/solute/solute.f90`
- Modify: `src/io/swap_csv_output.f90`
- Modify: `tests/unit/state/test_solute_state.pf`
- Modify: `tests/unit/state/test_solute_cumulative.pf`
- Modify: `tests/unit/state/test_swap_state.pf`

### Step 1.1: Rewrite `src/state/solute_state.f90` to flat structure

- [ ] **Replace the entire file contents with:**

```fortran
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
   implicit none
   private
   public :: solute_state_t

   type :: solute_state_t

      ! === per-node arrays (macp-sized; allocated by caller from config) ===
      real(real64), allocatable :: cml(:)    !! soil solute concentration (M/L3 water) in mobile region
      real(real64), allocatable :: cmsy(:)   !! dissolved + adsorbed solute concentration (M/L3 soil volume)

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
      procedure :: reset_intermediate => solute_reset_intermediate
      procedure :: reset_cumulative   => solute_reset_cumulative
   end type solute_state_t

contains

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
```

This removes the public exports `solute_intermediate_t` and `solute_cumulative_t`. Any test file or source file that imported these names by name will need its `use` statement updated in subsequent steps.

### Step 1.2: Rewrite test_solute_state.pf (remove cohort path references)

- [ ] **Replace `tests/unit/state/test_solute_state.pf` with:**

```fortran
! Tests for solute_state_t — typed state record for the solute subsystem.
! Cohort sub-records flattened 2026-05-12; reset procedures now bound on
! solute_state_t directly (reset_intermediate / reset_cumulative).

@test
subroutine test_solute_state_default_scalars()
   use funit
   use iso_fortran_env, only: real64
   use solute_state_mod, only: solute_state_t
   type(solute_state_t) :: sl

   ! cumulative fields
   @assertEqual(0.0_real64, sl%samini,   1.0e-12_real64)
   @assertEqual(0.0_real64, sl%dectot,   1.0e-12_real64)
   @assertEqual(0.0_real64, sl%rottot,   1.0e-12_real64)
   @assertEqual(0.0_real64, sl%sqprec,   1.0e-12_real64)
   @assertEqual(0.0_real64, sl%sqirrig,  1.0e-12_real64)
   @assertEqual(0.0_real64, sl%sqdra,    1.0e-12_real64)
   @assertEqual(0.0_real64, sl%sqbot,    1.0e-12_real64)
   @assertEqual(0.0_real64, sl%sqsur,    1.0e-12_real64)
   @assertEqual(0.0_real64, sl%csurf,    1.0e-12_real64)
   ! intermediate fields
   @assertEqual(0.0_real64, sl%imsqprec,  1.0e-12_real64)
   @assertEqual(0.0_real64, sl%imsqirrig, 1.0e-12_real64)
   @assertEqual(0.0_real64, sl%imsqbot,   1.0e-12_real64)
   @assertEqual(0.0_real64, sl%imsqdra,   1.0e-12_real64)
   @assertEqual(0.0_real64, sl%imdectot,  1.0e-12_real64)
   @assertEqual(0.0_real64, sl%imrottot,  1.0e-12_real64)
   ! other scalar fields
   @assertEqual(0.0_real64, sl%sampro,   1.0e-12_real64)
   @assertEqual(0.0_real64, sl%samcra,   1.0e-12_real64)
   @assertEqual(0.0_real64, sl%solbal,   1.0e-12_real64)
   @assertEqual(0.0_real64, sl%sqrap,    1.0e-12_real64)
   @assertEqual(0.0_real64, sl%cpond,    1.0e-12_real64)
   @assertEqual(0.0_real64, sl%cdrain,   1.0e-12_real64)
   @assertEqual(0.0_real64, sl%cseep,    1.0e-12_real64)
   @assertEqual(0.0_real64, sl%dtsolu,   1.0e-12_real64)
end subroutine test_solute_state_default_scalars

@test
subroutine test_solute_state_arrays_unallocated()
   use funit
   use solute_state_mod, only: solute_state_t
   type(solute_state_t) :: sl

   @assertFalse(allocated(sl%cml))
   @assertFalse(allocated(sl%cmsy))
end subroutine test_solute_state_arrays_unallocated

@test
subroutine test_solute_state_array_allocation()
   use funit
   use iso_fortran_env, only: real64
   use solute_state_mod, only: solute_state_t
   type(solute_state_t) :: sl
   integer, parameter :: numnod = 10

   allocate(sl%cml(numnod));   sl%cml  = 0.0_real64
   allocate(sl%cmsy(numnod));  sl%cmsy = 0.0_real64

   @assertTrue(allocated(sl%cml))
   @assertTrue(allocated(sl%cmsy))
   @assertEqual(numnod, size(sl%cml))
   @assertEqual(numnod, size(sl%cmsy))
end subroutine test_solute_state_array_allocation

@test
subroutine test_solute_state_independent_instances()
   use funit
   use iso_fortran_env, only: real64
   use solute_state_mod, only: solute_state_t
   type(solute_state_t) :: sl_a, sl_b

   sl_a%samini = 1.5_real64
   sl_b%samini = -2.5_real64

   @assertEqual(1.5_real64,  sl_a%samini, 1.0e-12_real64)
   @assertEqual(-2.5_real64, sl_b%samini, 1.0e-12_real64)
end subroutine test_solute_state_independent_instances
```

### Step 1.3: Rewrite test_solute_cumulative.pf (reset procedures now bound on solute_state_t)

- [ ] **Replace `tests/unit/state/test_solute_cumulative.pf` with:**

```fortran
! Tests for solute_state_t reset_intermediate() and reset_cumulative() —
! the named reset procedures introduced by the 2026-05-12 cohort-flattening
! arc. Replaces the cohort-based tests from CRR Phase B Task B2.
! No allocatable arrays in the reset-eligible scalars, so 4 tests (not 5).

@test
subroutine test_solute_reset_intermediate_zeroes_all_scalars()
   use funit
   use iso_fortran_env, only: real64
   use solute_state_mod, only: solute_state_t
   type(solute_state_t) :: sl

   sl%imsqprec  = 1.0_real64
   sl%imsqirrig = 2.0_real64
   sl%imsqbot   = 3.0_real64
   sl%imsqdra   = 4.0_real64
   sl%imdectot  = 5.0_real64
   sl%imrottot  = 6.0_real64

   call sl%reset_intermediate()

   @assertEqual(0.0_real64, sl%imsqprec,  1.0e-12_real64)
   @assertEqual(0.0_real64, sl%imsqirrig, 1.0e-12_real64)
   @assertEqual(0.0_real64, sl%imsqbot,   1.0e-12_real64)
   @assertEqual(0.0_real64, sl%imsqdra,   1.0e-12_real64)
   @assertEqual(0.0_real64, sl%imdectot,  1.0e-12_real64)
   @assertEqual(0.0_real64, sl%imrottot,  1.0e-12_real64)
end subroutine test_solute_reset_intermediate_zeroes_all_scalars

@test
subroutine test_solute_reset_intermediate_from_default_state()
   use funit
   use iso_fortran_env, only: real64
   use solute_state_mod, only: solute_state_t
   type(solute_state_t) :: sl

   ! Default state is already zero; reset must not crash.
   call sl%reset_intermediate()

   @assertEqual(0.0_real64, sl%imsqprec, 1.0e-12_real64)
end subroutine test_solute_reset_intermediate_from_default_state

@test
subroutine test_solute_reset_cumulative_zeroes_all_scalars()
   use funit
   use iso_fortran_env, only: real64
   use solute_state_mod, only: solute_state_t
   type(solute_state_t) :: sl

   sl%sqprec  = 1.0_real64
   sl%sqirrig = 2.0_real64
   sl%sqbot   = 3.0_real64
   sl%sqdra   = 4.0_real64
   sl%sqsur   = 5.0_real64
   sl%dectot  = 6.0_real64
   sl%rottot  = 7.0_real64
   sl%csurf   = 8.0_real64
   sl%samini  = 9.0_real64

   call sl%reset_cumulative()

   @assertEqual(0.0_real64, sl%sqprec,  1.0e-12_real64)
   @assertEqual(0.0_real64, sl%sqirrig, 1.0e-12_real64)
   @assertEqual(0.0_real64, sl%sqbot,   1.0e-12_real64)
   @assertEqual(0.0_real64, sl%sqdra,   1.0e-12_real64)
   @assertEqual(0.0_real64, sl%sqsur,   1.0e-12_real64)
   @assertEqual(0.0_real64, sl%dectot,  1.0e-12_real64)
   @assertEqual(0.0_real64, sl%rottot,  1.0e-12_real64)
   @assertEqual(0.0_real64, sl%csurf,   1.0e-12_real64)
   @assertEqual(0.0_real64, sl%samini,  1.0e-12_real64)
end subroutine test_solute_reset_cumulative_zeroes_all_scalars

@test
subroutine test_solute_reset_cumulative_from_default_state()
   use funit
   use iso_fortran_env, only: real64
   use solute_state_mod, only: solute_state_t
   type(solute_state_t) :: sl

   ! Default state is already zero; reset must not crash.
   call sl%reset_cumulative()

   @assertEqual(0.0_real64, sl%samini, 1.0e-12_real64)
end subroutine test_solute_reset_cumulative_from_default_state
```

### Step 1.4: Migrate call sites in `src/solute/solute.f90`

- [ ] **Find every occurrence of `state%solute%intermediate%` and `state%solute%cumulative%` and remove the `%intermediate` / `%cumulative` segment.**

Use sed (verify with grep before and after):

```bash
cd /home/zawadzkim/Code/swap

# preview the affected lines
grep -n "state%solute%intermediate%\|state%solute%cumulative%" src/solute/solute.f90

# apply the rewrite (intermediate/cumulative collapse)
sed -i 's/state%solute%intermediate%/state%solute%/g; s/state%solute%cumulative%/state%solute%/g' src/solute/solute.f90

# confirm zero remaining
grep -n "state%solute%intermediate%\|state%solute%cumulative%" src/solute/solute.f90 || echo "clean"
```

- [ ] **Find every reset call and update to the new procedure names.**

Original (locate by grep):

```fortran
if (flzerointr) call state%solute%intermediate%reset()
```

Becomes:

```fortran
if (flzerointr) call state%solute%reset_intermediate()
```

And:

```fortran
if (flzerocumu) call state%solute%cumulative%reset()
```

Becomes:

```fortran
if (flzerocumu) call state%solute%reset_cumulative()
```

Apply via sed (verify first):

```bash
grep -n "state%solute%intermediate%reset\|state%solute%cumulative%reset" src/solute/solute.f90
sed -i 's/state%solute%intermediate%reset()/state%solute%reset_intermediate()/g; s/state%solute%cumulative%reset()/state%solute%reset_cumulative()/g' src/solute/solute.f90
grep -n "state%solute%intermediate%reset\|state%solute%cumulative%reset" src/solute/solute.f90 || echo "clean"
```

- [ ] **Update ASSOCIATE block targets in `src/solute/solute.f90`.**

Per-field aliases like:

```fortran
samini => state%solute%cumulative%samini, &
csurf  => state%solute%cumulative%csurf,  &
```

Have their RHS targets rewritten by the same path-collapse sed above. After the sed, the ASSOCIATE blocks now read:

```fortran
samini => state%solute%samini, &
csurf  => state%solute%csurf,  &
```

Confirm:

```bash
grep -n "=> state%solute%intermediate%\|=> state%solute%cumulative%" src/solute/solute.f90 || echo "clean"
```

### Step 1.5: Migrate call sites in `src/io/swap_csv_output.f90`

- [ ] **Apply the same path-collapse rewrite to swap_csv_output.f90:**

```bash
cd /home/zawadzkim/Code/swap
grep -n "state%solute%intermediate%\|state%solute%cumulative%" src/io/swap_csv_output.f90
sed -i 's/state%solute%intermediate%/state%solute%/g; s/state%solute%cumulative%/state%solute%/g' src/io/swap_csv_output.f90
grep -n "state%solute%intermediate%\|state%solute%cumulative%" src/io/swap_csv_output.f90 || echo "clean"
```

### Step 1.6: Migrate `tests/unit/state/test_swap_state.pf`

- [ ] **Apply the same path-collapse rewrite to test_swap_state.pf:**

```bash
cd /home/zawadzkim/Code/swap
grep -n "state%solute%intermediate%\|state%solute%cumulative%" tests/unit/state/test_swap_state.pf
sed -i 's/state%solute%intermediate%/state%solute%/g; s/state%solute%cumulative%/state%solute%/g' tests/unit/state/test_swap_state.pf
grep -n "state%solute%intermediate%\|state%solute%cumulative%" tests/unit/state/test_swap_state.pf || echo "clean"
```

### Step 1.7: Build & test

- [ ] **Build the project:**

```bash
cd /home/zawadzkim/Code/swap
pixi run build-linux
```

Expected: clean build, no errors, no new warnings. If errors mention `intermediate` or `cumulative` as a missing component, return to Steps 1.4-1.6 and search for missed paths.

- [ ] **Run the pFUnit suite:**

```bash
pixi run test-pfunit
```

Expected: all tests pass (including the 4 new tests in `test_solute_cumulative.pf` and the 4 tests in `test_solute_state.pf`).

### Step 1.8: Audit grep zero-out (solute)

- [ ] **Confirm zero remaining cohort paths anywhere in the tree:**

```bash
cd /home/zawadzkim/Code/swap
grep -rn 'state%solute%intermediate%\|state%solute%cumulative%' src/ tests/
```

Expected: no output (zero hits).

- [ ] **Also confirm no leftover `solute_intermediate_t` or `solute_cumulative_t` references:**

```bash
grep -rn 'solute_intermediate_t\|solute_cumulative_t' src/ tests/
```

Expected: no output. The two cohort types are gone.

### Step 1.9: Run the full check-full gate

- [ ] **Run the full check (build + pFUnit + regression):**

```bash
cd /home/zawadzkim/Code/swap
pixi run check-full
```

Expected: all targets pass, regression suite produces bit-identical output to the pre-refactor baseline. Per project workflow convention, this must pass **before commit** — pFUnit alone misses global-default regressions.

If regression diff appears: stop, do not commit. Investigate by re-reading the call sites you touched and confirming you didn't accidentally rename a non-cohort field (paranoid double-check via `git diff src/`).

### Step 1.10: Commit

- [ ] **Stage and commit:**

```bash
cd /home/zawadzkim/Code/swap
git add src/state/solute_state.f90 \
        src/solute/solute.f90 \
        src/io/swap_csv_output.f90 \
        tests/unit/state/test_solute_state.pf \
        tests/unit/state/test_solute_cumulative.pf \
        tests/unit/state/test_swap_state.pf

git commit -m "$(cat <<'EOF'
refactor(state): flatten solute cohort sub-records into solute_state_t

Replaces the nested solute_intermediate_t and solute_cumulative_t
cohorts (ADR 0033 Phase B) with flat fields on solute_state_t and
named reset procedures: reset_intermediate (flzerointr gate) and
reset_cumulative (flzerocumu gate). All 25 call-site references in
solute.f90, swap_csv_output.f90, and the unit tests rewritten to the
shallower paths. samini = sampro physics rebase preserved inline at
the call site. check-full byte-identical.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Soilwater — flatten cohorts and rebind resets

**Files:**
- Modify: `src/state/soilwater_state.f90` (cohort type deletions + reset rebinding; allocation block in `soilwater_init` has path rewrites only — no allocation changes)
- Modify: `src/soil/waterbalance.f90`
- Modify: `src/io/swap_csv_output.f90`
- Modify: `src/core/variables.f90`
- Modify: `src/io/swapoutput.f90`
- Modify: `src/crop/cropgrowth.f90`
- Modify: `src/core/initialize.f90`
- Modify: `src/crop/irrigation.f90`
- Modify: `src/soil/soilhydraulics.f90`
- Modify: `src/soil/soilgrid.f90`
- Modify: `src/crop/management_soil.f90`
- Modify: `src/core/swap.f90`
- Modify: `src/atmosphere/meteoday.f90`
- Modify: `tests/unit/state/test_soilwater_state.pf`

### Step 2.1: Rewrite soilwater state module

- [ ] **Edit `src/state/soilwater_state.f90`:**

1. Delete the type definitions of `soilwater_intermediate_t` (lines 75-121 in current file) and `soilwater_cumulative_t` (lines 134-151).
2. Delete the lines that declare the cohort components inside `soilwater_state_t`:
   ```fortran
   type(soilwater_intermediate_t) :: intr
   type(soilwater_cumulative_t)   :: cumu
   ```
3. In their place, paste the field declarations from those two cohort types directly inside `soilwater_state_t`, **with section-header comments**. Insert after the existing core-flat section but before `contains` (which today is at the bottom of the type definition):

```fortran
      ! ===========================================================================
      ! INTERMEDIATE accumulators (reset_intermediate / gate: flzerointr)
      !   subsumes the per-day subset, which has its own reset under flDayStart.
      ! ===========================================================================

      ! Non-per-day per-node arrays (allocated by soilwater_init)
      real(real64), allocatable :: inq(:)        !< intra-period inter-comp flux (cm)
      real(real64), allocatable :: inqrot(:)     !< intra-period root uptake (cm)
      real(real64), allocatable :: inqssdi(:)    !< intra-period SSDI flux (cm)
      real(real64), allocatable :: iqdo(:)       !< intra-period downward flux (cm)
      real(real64), allocatable :: iqup(:)       !< intra-period upward flux (cm)
      real(real64), allocatable :: IThetaBeg(:)  !< theta at start of intr period (-)

      ! Non-per-day scalars
      real(real64) :: iqrot     = 0.0_real64
      real(real64) :: iqssdi    = 0.0_real64
      real(real64) :: iqredwet  = 0.0_real64
      real(real64) :: iqreddry  = 0.0_real64
      real(real64) :: iqredsol  = 0.0_real64
      real(real64) :: iqredfrs  = 0.0_real64
      real(real64) :: ies0      = 0.0_real64
      real(real64) :: iet0      = 0.0_real64
      real(real64) :: iew0      = 0.0_real64
      real(real64) :: iintc     = 0.0_real64
      real(real64) :: iruno     = 0.0_real64
      real(real64) :: irunoCN   = 0.0_real64
      real(real64) :: irunon    = 0.0_real64
      real(real64) :: iqbot     = 0.0_real64
      real(real64) :: iqtdo     = 0.0_real64
      real(real64) :: iqtup     = 0.0_real64
      real(real64) :: IPondBeg  = 0.0_real64
      real(real64) :: iprec     = 0.0_real64
      real(real64) :: igird     = 0.0_real64
      real(real64) :: inird     = 0.0_real64

      ! ===========================================================================
      ! PER-DAY subset (reset_intermediate_per_day / gate: flDayStart)
      !   Also zeroed as part of reset_intermediate (flzerointr subsumes flDayStart).
      ! ===========================================================================
      real(real64) :: tra           = 0.0_real64
      real(real64) :: iqredwet_day  = 0.0_real64
      real(real64) :: iqreddry_day  = 0.0_real64
      real(real64) :: iqredsol_day  = 0.0_real64
      real(real64) :: iqredfrs_day  = 0.0_real64
      real(real64) :: iptra_day     = 0.0_real64
      real(real64), allocatable :: qpotrot_day(:)
      real(real64), allocatable :: qredtot_day(:)

      ! ===========================================================================
      ! CUMULATIVE accumulators (reset_cumulative / gate: flzerocumu)
      ! ===========================================================================
      real(real64) :: cqssdi    = 0.0_real64
      real(real64) :: cqrot     = 0.0_real64
      real(real64) :: cqbot     = 0.0_real64
      real(real64) :: cqbotdo   = 0.0_real64
      real(real64) :: cqbotup   = 0.0_real64
      real(real64) :: cinund    = 0.0_real64
      real(real64) :: crunon    = 0.0_real64
      real(real64) :: crunoff   = 0.0_real64
      real(real64) :: crunoffCN = 0.0_real64
      real(real64) :: cqtdo     = 0.0_real64
      real(real64) :: cqtup     = 0.0_real64
      real(real64) :: cqprai    = 0.0_real64
      real(real64) :: cgird     = 0.0_real64
      real(real64) :: cnird     = 0.0_real64
```

4. Replace the existing `contains` block at the end of the `soilwater_state_t` definition (which currently has no procedures bound, since cohorts owned them) with:

```fortran
   contains
      procedure :: reset_intermediate         => soilwater_reset_intermediate
      procedure :: reset_intermediate_per_day => soilwater_reset_intermediate_per_day
      procedure :: reset_cumulative           => soilwater_reset_cumulative
   end type soilwater_state_t
```

5. Update the `public ::` list near the top of the module to remove `soilwater_intermediate_t` and `soilwater_cumulative_t` exports.

6. Rewrite the three reset procedures at the bottom of the module:

```fortran
   !> Zero ALL fields in the intermediate cohort — flzerointr gate.
   !! Includes the per-day subset (flzerointr subsumes flDayStart).
   subroutine soilwater_reset_intermediate(self)
      class(soilwater_state_t), intent(inout) :: self

      ! Non-per-day allocatable arrays
      if (allocated(self%inq))       self%inq       = 0.0_real64
      if (allocated(self%inqrot))    self%inqrot    = 0.0_real64
      if (allocated(self%inqssdi))   self%inqssdi   = 0.0_real64
      if (allocated(self%iqdo))      self%iqdo      = 0.0_real64
      if (allocated(self%iqup))      self%iqup      = 0.0_real64
      if (allocated(self%IThetaBeg)) self%IThetaBeg = 0.0_real64

      ! Non-per-day scalars
      self%iqrot    = 0.0_real64
      self%iqssdi   = 0.0_real64
      self%iqredwet = 0.0_real64
      self%iqreddry = 0.0_real64
      self%iqredsol = 0.0_real64
      self%iqredfrs = 0.0_real64
      self%ies0     = 0.0_real64
      self%iet0     = 0.0_real64
      self%iew0     = 0.0_real64
      self%iintc    = 0.0_real64
      self%iruno    = 0.0_real64
      self%irunoCN  = 0.0_real64
      self%irunon   = 0.0_real64
      self%iqbot    = 0.0_real64
      self%iqtdo    = 0.0_real64
      self%iqtup    = 0.0_real64
      self%IPondBeg = 0.0_real64
      self%iprec    = 0.0_real64
      self%igird    = 0.0_real64
      self%inird    = 0.0_real64

      ! Per-day fields — also zeroed by full reset (flzerointr subsumes flDayStart)
      call soilwater_reset_intermediate_per_day(self)

   end subroutine soilwater_reset_intermediate

   !> Zero only the per-day subset — flDayStart gate.
   !! Non-per-day intermediate fields are UNTOUCHED.
   subroutine soilwater_reset_intermediate_per_day(self)
      class(soilwater_state_t), intent(inout) :: self
      self%tra          = 0.0_real64
      self%iqredwet_day = 0.0_real64
      self%iqreddry_day = 0.0_real64
      self%iqredsol_day = 0.0_real64
      self%iqredfrs_day = 0.0_real64
      self%iptra_day    = 0.0_real64
      if (allocated(self%qpotrot_day)) self%qpotrot_day = 0.0_real64
      if (allocated(self%qredtot_day)) self%qredtot_day = 0.0_real64
   end subroutine soilwater_reset_intermediate_per_day

   !> Zero the 14 cumulative fields — flzerocumu gate.
   subroutine soilwater_reset_cumulative(self)
      class(soilwater_state_t), intent(inout) :: self
      self%cqssdi    = 0.0_real64
      self%cqrot     = 0.0_real64
      self%cqbot     = 0.0_real64
      self%cqbotdo   = 0.0_real64
      self%cqbotup   = 0.0_real64
      self%cinund    = 0.0_real64
      self%crunon    = 0.0_real64
      self%crunoff   = 0.0_real64
      self%crunoffCN = 0.0_real64
      self%cqtdo     = 0.0_real64
      self%cqtup     = 0.0_real64
      self%cqprai    = 0.0_real64
      self%cgird     = 0.0_real64
      self%cnird     = 0.0_real64
   end subroutine soilwater_reset_cumulative
```

7. Inside `subroutine soilwater_init(sw, numnod, nlay)` (allocation block around lines 405-419 in the current file), rewrite the per-field paths from `sw%intr%X` to `sw%X`:

Before:

```fortran
allocate(sw%intr%inqrot(numnod));    sw%intr%inqrot    = 0.0_real64
allocate(sw%intr%inqssdi(numnod));   sw%intr%inqssdi   = 0.0_real64
allocate(sw%intr%IThetaBeg(numnod)); sw%intr%IThetaBeg = 0.0_real64
allocate(sw%intr%inq(numnod+1));     sw%intr%inq       = 0.0_real64
allocate(sw%intr%iqdo(numnod+1));    sw%intr%iqdo      = 0.0_real64
allocate(sw%intr%iqup(numnod+1));    sw%intr%iqup      = 0.0_real64
allocate(sw%intr%qpotrot_day(numnod)); sw%intr%qpotrot_day = 0.0_real64
allocate(sw%intr%qredtot_day(numnod)); sw%intr%qredtot_day = 0.0_real64
```

After:

```fortran
allocate(sw%inqrot(numnod));    sw%inqrot    = 0.0_real64
allocate(sw%inqssdi(numnod));   sw%inqssdi   = 0.0_real64
allocate(sw%IThetaBeg(numnod)); sw%IThetaBeg = 0.0_real64
allocate(sw%inq(numnod+1));     sw%inq       = 0.0_real64
allocate(sw%iqdo(numnod+1));    sw%iqdo      = 0.0_real64
allocate(sw%iqup(numnod+1));    sw%iqup      = 0.0_real64
allocate(sw%qpotrot_day(numnod)); sw%qpotrot_day = 0.0_real64
allocate(sw%qredtot_day(numnod)); sw%qredtot_day = 0.0_real64
```

The duplication between this allocator and the inline allocators elsewhere is **intentionally left for a follow-up arc** — do not consolidate here.

### Step 2.2: Rewrite tests/unit/state/test_soilwater_state.pf

- [ ] **Find every `state%soilwater%intr%` and `state%soilwater%cumu%` (and uses of the cohort type names) and rewrite paths.**

```bash
cd /home/zawadzkim/Code/swap

grep -n "state%soilwater%intr%\|state%soilwater%cumu%\|soilwater_intermediate_t\|soilwater_cumulative_t" tests/unit/state/test_soilwater_state.pf

sed -i \
   -e 's/state%soilwater%intr%/state%soilwater%/g' \
   -e 's/state%soilwater%cumu%/state%soilwater%/g' \
   tests/unit/state/test_soilwater_state.pf

grep -n "state%soilwater%intr%\|state%soilwater%cumu%" tests/unit/state/test_soilwater_state.pf || echo "paths clean"
```

- [ ] **If the test file declares typed locals of `soilwater_intermediate_t` / `soilwater_cumulative_t`, refactor them to declare `soilwater_state_t` instead and operate on the flat fields directly. Read the test file after the sed to confirm what's left.** Update procedure calls: `cohort%reset()` → `sw%reset_intermediate()` / `sw%reset_intermediate_per_day()` / `sw%reset_cumulative()` depending on which cohort the test was exercising.

- [ ] **If the test file imported the cohort type names via `use`, remove those imports:**

```bash
grep -n "use soilwater_state_mod" tests/unit/state/test_soilwater_state.pf
# Manually edit any line that lists soilwater_intermediate_t or soilwater_cumulative_t
# in an "only:" clause; remove just those names from the list.
```

Confirm:

```bash
grep -n "soilwater_intermediate_t\|soilwater_cumulative_t" tests/unit/state/test_soilwater_state.pf || echo "type refs clean"
```

### Step 2.3: Migrate call sites in the source tree

- [ ] **Apply the mechanical path rewrite to all 12 source files:**

```bash
cd /home/zawadzkim/Code/swap

for f in src/soil/waterbalance.f90 \
         src/io/swap_csv_output.f90 \
         src/core/variables.f90 \
         src/io/swapoutput.f90 \
         src/crop/cropgrowth.f90 \
         src/core/initialize.f90 \
         src/crop/irrigation.f90 \
         src/soil/soilhydraulics.f90 \
         src/soil/soilgrid.f90 \
         src/crop/management_soil.f90 \
         src/core/swap.f90 \
         src/atmosphere/meteoday.f90; do
    echo "=== $f ==="
    grep -c "state%soilwater%intr%\|state%soilwater%cumu%" "$f"
    sed -i \
       -e 's/state%soilwater%intr%/state%soilwater%/g' \
       -e 's/state%soilwater%cumu%/state%soilwater%/g' \
       "$f"
    grep -c "state%soilwater%intr%\|state%soilwater%cumu%" "$f"
done
```

The second `grep -c` should print `0` for every file.

- [ ] **Update reset calls. There may be `state%soilwater%intr%reset()`, `state%soilwater%intr%reset_per_day()`, and `state%soilwater%cumu%reset()` calls.** After the path-collapse sed above, those now read `state%soilwater%reset()`, `state%soilwater%reset_per_day()`, and `state%soilwater%reset()` respectively — which is broken (the intermediate and cumulative resets now collide on the same name `reset()`). Fix the procedure names:

```bash
cd /home/zawadzkim/Code/swap

# Find the reset call sites that need disambiguation
grep -rn "state%soilwater%reset()" src/ tests/
# These were all originally either intr%reset() or cumu%reset(). You need to
# look at each call site to know which it was — recover from git:
git diff HEAD -- src/ tests/ | grep -B1 'state%soilwater%reset()'
```

The cleaner approach: do the reset-call rewrite **before** the path-collapse sed. Use the following order instead (re-run from a clean working tree if you already ran the path-collapse sed):

```bash
cd /home/zawadzkim/Code/swap
git checkout -- src/ tests/  # only if you need to redo Step 2.3 from scratch

# Step A: reset-call rewrites (cohort-specific procedure names)
for f in $(grep -rl "state%soilwater%intr%reset\|state%soilwater%cumu%reset" src/ tests/); do
    sed -i \
       -e 's/state%soilwater%intr%reset_per_day()/state%soilwater%reset_intermediate_per_day()/g' \
       -e 's/state%soilwater%intr%reset()/state%soilwater%reset_intermediate()/g' \
       -e 's/state%soilwater%cumu%reset()/state%soilwater%reset_cumulative()/g' \
       "$f"
done

# Step B: now the bare path-collapse for the rest
for f in src/soil/waterbalance.f90 \
         src/io/swap_csv_output.f90 \
         src/core/variables.f90 \
         src/io/swapoutput.f90 \
         src/crop/cropgrowth.f90 \
         src/core/initialize.f90 \
         src/crop/irrigation.f90 \
         src/soil/soilhydraulics.f90 \
         src/soil/soilgrid.f90 \
         src/crop/management_soil.f90 \
         src/core/swap.f90 \
         src/atmosphere/meteoday.f90 \
         tests/unit/state/test_soilwater_state.pf; do
    sed -i \
       -e 's/state%soilwater%intr%/state%soilwater%/g' \
       -e 's/state%soilwater%cumu%/state%soilwater%/g' \
       "$f"
done
```

(If you already ran Step 2.3's original sed and reverted: use the order above.)

- [ ] **Also fix ASSOCIATE alias RHS targets.** Per-field aliases like `sw_inq => state%soilwater%intr%inq` are already collapsed by the path-collapse sed (the RHS is rewritten). Confirm:

```bash
grep -rn "=> state%soilwater%intr%\|=> state%soilwater%cumu%" src/ tests/ || echo "associates clean"
```

### Step 2.4: Build & test

- [ ] **Build:**

```bash
cd /home/zawadzkim/Code/swap
pixi run build-linux
```

Expected: clean build. Watch for:
- `intr` / `cumu` component-not-found errors (missed path)
- Procedure-name collisions if you skipped Step 2.3 Step A above
- Unused-variable warnings for any removed local aliases (clean up)

- [ ] **Run pFUnit suite:**

```bash
pixi run test-pfunit
```

Expected: all tests pass.

### Step 2.5: Audit grep zero-out (soilwater)

- [ ] **Confirm zero remaining cohort paths and type references:**

```bash
cd /home/zawadzkim/Code/swap
grep -rn 'state%soilwater%intr%\|state%soilwater%cumu%' src/ tests/
grep -rn 'soilwater_intermediate_t\|soilwater_cumulative_t' src/ tests/
```

Expected: both produce no output.

### Step 2.6: Run check-full

- [ ] **Run the full check (build + pFUnit + regression):**

```bash
cd /home/zawadzkim/Code/swap
pixi run check-full
```

Expected: all targets pass, regression byte-identical.

If regression diff appears, the most likely cause is a missed reset-procedure rename (e.g. `intr%reset()` collapsed to `reset()` with wrong target). Check `git diff` against the call sites that previously had reset calls.

### Step 2.7: Commit

- [ ] **Stage and commit:**

```bash
cd /home/zawadzkim/Code/swap

git add src/state/soilwater_state.f90 \
        src/soil/waterbalance.f90 \
        src/io/swap_csv_output.f90 \
        src/core/variables.f90 \
        src/io/swapoutput.f90 \
        src/crop/cropgrowth.f90 \
        src/core/initialize.f90 \
        src/crop/irrigation.f90 \
        src/soil/soilhydraulics.f90 \
        src/soil/soilgrid.f90 \
        src/crop/management_soil.f90 \
        src/core/swap.f90 \
        src/atmosphere/meteoday.f90 \
        tests/unit/state/test_soilwater_state.pf

git commit -m "$(cat <<'EOF'
refactor(state): flatten soilwater cohort sub-records into soilwater_state_t

Replaces the nested soilwater_intermediate_t (with subset reset_per_day)
and soilwater_cumulative_t cohorts (ADR 0038 / 0033 pattern) with flat
fields on soilwater_state_t and three named reset procedures:
  - reset_intermediate         (flzerointr gate; subsumes per-day subset)
  - reset_intermediate_per_day (flDayStart gate; per-day subset only)
  - reset_cumulative           (flzerocumu gate)
All ~190 call-site references across waterbalance.f90, swap_csv_output.f90,
swapoutput.f90, cropgrowth.f90, initialize.f90, irrigation.f90,
soilhydraulics.f90, soilgrid.f90, management_soil.f90, swap.f90,
meteoday.f90, variables.f90, and the unit tests rewritten to the
shallower paths. Allocation block in soilwater_init has path rewrites
only; allocation consolidation deferred to follow-up arc. check-full
byte-identical.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Surfacewater — flatten cohorts and rebind resets

**Files:**
- Modify: `src/state/surfacewater_state.f90`
- Modify: `src/drainage/drainage.f90`
- Modify: `src/soil/waterbalance.f90`
- Modify: `src/drainage/surfacewater.f90`
- Modify: `src/io/swapoutput.f90`
- Modify: `src/io/swap_csv_output.f90`
- Modify: `src/soil/soilgrid.f90`
- Modify: `src/crop/management_soil.f90`
- Modify: `src/drainage/surfacewater_init.f90` (path rewrites only — no allocation changes)
- Modify: `tests/unit/state/test_surfacewater_state.pf`
- Modify: `tests/unit/state/test_surfacewater_cumulative.pf`

### Step 3.1: Rewrite `src/state/surfacewater_state.f90`

- [ ] **Replace the entire file contents with:**

```fortran
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
```

This removes the public exports for the three cohort types.

### Step 3.2: Rewrite test_surfacewater_state.pf

- [ ] **Apply the path-collapse rewrite. Then audit and rewrite cohort-typed locals to use the parent type.**

```bash
cd /home/zawadzkim/Code/swap

sed -i \
   -e 's/state%surfacewater%intermediate%/state%surfacewater%/g' \
   -e 's/state%surfacewater%drainage_cumulative%/state%surfacewater%/g' \
   -e 's/state%surfacewater%reservoir_cumulative%/state%surfacewater%/g' \
   tests/unit/state/test_surfacewater_state.pf

grep -n "state%surfacewater%intermediate%\|state%surfacewater%drainage_cumulative%\|state%surfacewater%reservoir_cumulative%" tests/unit/state/test_surfacewater_state.pf || echo "paths clean"
```

- [ ] **Read the file and update any `use surfacewater_state_mod, only: surfacewater_intermediate_t, ...` imports** — remove the cohort type names from the `only:` list.

```bash
grep -n "surfacewater_intermediate_t\|surfacewater_drainage_cumulative_t\|surfacewater_reservoir_cumulative_t" tests/unit/state/test_surfacewater_state.pf
```

For each import that lists a cohort type name, delete just that name from the list. Local variable declarations of cohort types — re-declare as `surfacewater_state_t` and update field access accordingly.

### Step 3.3: Rewrite test_surfacewater_cumulative.pf

The current file tests cohort `reset()` procedures. Rewrite to test the new parent-bound procedures.

- [ ] **Replace the entire file contents with:**

```fortran
! Tests for surfacewater_state_t reset procedures:
!   - reset_intermediate          (gate: flzerointr)
!   - reset_cumulative_drainage   (gate: flzerocumu + fldrain owner)
!   - reset_cumulative_reservoir  (gate: flzerocumu + flSurfaceWater owner)
!
! The original cohort sub-record types were flattened by the 2026-05-12
! reset-cohort-flattening arc (ADR 0042); reset procedures now bind on
! the parent state type directly.

@test
subroutine test_intermediate_reset_zeroes_scalar_iqdra()
   use funit
   use iso_fortran_env, only: real64
   use surfacewater_state_mod, only: surfacewater_state_t
   type(surfacewater_state_t) :: sw

   sw%iqdra = 42.0_real64
   call sw%reset_intermediate()
   @assertEqual(0.0_real64, sw%iqdra, 1.0e-12_real64)
end subroutine test_intermediate_reset_zeroes_scalar_iqdra

@test
subroutine test_intermediate_reset_zeroes_allocated_arrays()
   use funit
   use iso_fortran_env, only: real64
   use surfacewater_state_mod, only: surfacewater_state_t
   type(surfacewater_state_t) :: sw

   allocate(sw%inqdra(2, 3));     sw%inqdra     = 1.0_real64
   allocate(sw%inqdra_in(2, 3));  sw%inqdra_in  = 2.0_real64
   allocate(sw%inqdra_out(2, 3)); sw%inqdra_out = 3.0_real64

   call sw%reset_intermediate()

   @assertEqual(0.0_real64, sum(sw%inqdra),     1.0e-12_real64)
   @assertEqual(0.0_real64, sum(sw%inqdra_in),  1.0e-12_real64)
   @assertEqual(0.0_real64, sum(sw%inqdra_out), 1.0e-12_real64)
end subroutine test_intermediate_reset_zeroes_allocated_arrays

@test
subroutine test_intermediate_reset_safe_when_arrays_not_allocated()
   use funit
   use surfacewater_state_mod, only: surfacewater_state_t
   type(surfacewater_state_t) :: sw

   call sw%reset_intermediate()

   @assertFalse(allocated(sw%inqdra))
   @assertFalse(allocated(sw%inqdra_in))
   @assertFalse(allocated(sw%inqdra_out))
end subroutine test_intermediate_reset_safe_when_arrays_not_allocated

@test
subroutine test_dr_cum_reset_zeroes_cqdra()
   use funit
   use iso_fortran_env, only: real64
   use surfacewater_state_mod, only: surfacewater_state_t
   type(surfacewater_state_t) :: sw

   sw%cqdra = 42.0_real64
   call sw%reset_cumulative_drainage()
   @assertEqual(0.0_real64, sw%cqdra, 1.0e-12_real64)
end subroutine test_dr_cum_reset_zeroes_cqdra

@test
subroutine test_dr_cum_reset_zeroes_per_level_arrays()
   use funit
   use iso_fortran_env, only: real64
   use surfacewater_state_mod, only: surfacewater_state_t
   type(surfacewater_state_t) :: sw

   allocate(sw%cqdrain(5));    sw%cqdrain    = 5.0_real64
   allocate(sw%cqdrainin(5));  sw%cqdrainin  = 5.0_real64
   allocate(sw%cqdrainout(5)); sw%cqdrainout = 5.0_real64

   call sw%reset_cumulative_drainage()

   @assertEqual(0.0_real64, sum(sw%cqdrain),    1.0e-12_real64)
   @assertEqual(0.0_real64, sum(sw%cqdrainin),  1.0e-12_real64)
   @assertEqual(0.0_real64, sum(sw%cqdrainout), 1.0e-12_real64)
end subroutine test_dr_cum_reset_zeroes_per_level_arrays

@test
subroutine test_dr_cum_reset_safe_unallocated()
   use funit
   use surfacewater_state_mod, only: surfacewater_state_t
   type(surfacewater_state_t) :: sw

   call sw%reset_cumulative_drainage()

   @assertFalse(allocated(sw%cqdrain))
   @assertFalse(allocated(sw%cqdrainin))
   @assertFalse(allocated(sw%cqdrainout))
end subroutine test_dr_cum_reset_safe_unallocated

@test
subroutine test_res_cum_reset_zeroes_all_scalars()
   use funit
   use iso_fortran_env, only: real64
   use surfacewater_state_mod, only: surfacewater_state_t
   type(surfacewater_state_t) :: sw

   sw%cqdrd  = 1.0_real64
   sw%cwsupp = 2.0_real64
   sw%cwout  = 3.0_real64

   call sw%reset_cumulative_reservoir()

   @assertEqual(0.0_real64, sw%cqdrd,  1.0e-12_real64)
   @assertEqual(0.0_real64, sw%cwsupp, 1.0e-12_real64)
   @assertEqual(0.0_real64, sw%cwout,  1.0e-12_real64)
end subroutine test_res_cum_reset_zeroes_all_scalars
```

### Step 3.4: Migrate call sites in the source tree

Apply reset-call rewrites **before** path-collapse to avoid collisions (same pattern as Task 2). The cohort `%reset()` calls map to distinct procedure names on the flat type:

- `state%surfacewater%intermediate%reset()` → `state%surfacewater%reset_intermediate()`
- `state%surfacewater%drainage_cumulative%reset()` → `state%surfacewater%reset_cumulative_drainage()`
- `state%surfacewater%reservoir_cumulative%reset()` → `state%surfacewater%reset_cumulative_reservoir()`

There are also alias forms via `sw => state%surfacewater` in `surfacewater_init.f90`. Cover those too.

- [ ] **Reset-call rewrites (Step A):**

```bash
cd /home/zawadzkim/Code/swap

for f in $(grep -rl "state%surfacewater%intermediate%reset\|state%surfacewater%drainage_cumulative%reset\|state%surfacewater%reservoir_cumulative%reset\|sw%intermediate%reset\|sw%drainage_cumulative%reset\|sw%reservoir_cumulative%reset" src/ tests/); do
    sed -i \
       -e 's/state%surfacewater%intermediate%reset()/state%surfacewater%reset_intermediate()/g' \
       -e 's/state%surfacewater%drainage_cumulative%reset()/state%surfacewater%reset_cumulative_drainage()/g' \
       -e 's/state%surfacewater%reservoir_cumulative%reset()/state%surfacewater%reset_cumulative_reservoir()/g' \
       -e 's/sw%intermediate%reset()/sw%reset_intermediate()/g' \
       -e 's/sw%drainage_cumulative%reset()/sw%reset_cumulative_drainage()/g' \
       -e 's/sw%reservoir_cumulative%reset()/sw%reset_cumulative_reservoir()/g' \
       "$f"
done
```

- [ ] **Path collapse (Step B):**

```bash
cd /home/zawadzkim/Code/swap

for f in src/drainage/drainage.f90 \
         src/soil/waterbalance.f90 \
         src/drainage/surfacewater.f90 \
         src/io/swapoutput.f90 \
         src/io/swap_csv_output.f90 \
         src/soil/soilgrid.f90 \
         src/crop/management_soil.f90 \
         src/drainage/surfacewater_init.f90 \
         tests/unit/state/test_surfacewater_state.pf; do
    echo "=== $f ==="
    grep -c "state%surfacewater%intermediate%\|state%surfacewater%drainage_cumulative%\|state%surfacewater%reservoir_cumulative%\|sw%intermediate%\|sw%drainage_cumulative%\|sw%reservoir_cumulative%" "$f"
    sed -i \
       -e 's/state%surfacewater%intermediate%/state%surfacewater%/g' \
       -e 's/state%surfacewater%drainage_cumulative%/state%surfacewater%/g' \
       -e 's/state%surfacewater%reservoir_cumulative%/state%surfacewater%/g' \
       -e 's/sw%intermediate%/sw%/g' \
       -e 's/sw%drainage_cumulative%/sw%/g' \
       -e 's/sw%reservoir_cumulative%/sw%/g' \
       "$f"
    grep -c "state%surfacewater%intermediate%\|state%surfacewater%drainage_cumulative%\|state%surfacewater%reservoir_cumulative%\|sw%intermediate%\|sw%drainage_cumulative%\|sw%reservoir_cumulative%" "$f"
done
```

The second `grep -c` should print `0` for every file.

### Step 3.5: Build & test

- [ ] **Build:**

```bash
cd /home/zawadzkim/Code/swap
pixi run build-linux
```

Expected: clean build.

- [ ] **Run pFUnit:**

```bash
pixi run test-pfunit
```

Expected: all tests pass.

### Step 3.6: Audit grep zero-out (surfacewater)

- [ ] **Confirm zero remaining cohort paths and type references:**

```bash
cd /home/zawadzkim/Code/swap

grep -rn 'state%surfacewater%intermediate%\|state%surfacewater%drainage_cumulative%\|state%surfacewater%reservoir_cumulative%' src/ tests/
grep -rn 'sw%intermediate%\|sw%drainage_cumulative%\|sw%reservoir_cumulative%' src/ tests/
grep -rn 'surfacewater_intermediate_t\|surfacewater_drainage_cumulative_t\|surfacewater_reservoir_cumulative_t' src/ tests/
```

Expected: all three produce no output.

### Step 3.7: Run check-full

- [ ] **Run the full check:**

```bash
cd /home/zawadzkim/Code/swap
pixi run check-full
```

Expected: all targets pass, regression byte-identical.

### Step 3.8: Commit

- [ ] **Stage and commit:**

```bash
cd /home/zawadzkim/Code/swap

git add src/state/surfacewater_state.f90 \
        src/drainage/drainage.f90 \
        src/soil/waterbalance.f90 \
        src/drainage/surfacewater.f90 \
        src/io/swapoutput.f90 \
        src/io/swap_csv_output.f90 \
        src/soil/soilgrid.f90 \
        src/crop/management_soil.f90 \
        src/drainage/surfacewater_init.f90 \
        tests/unit/state/test_surfacewater_state.pf \
        tests/unit/state/test_surfacewater_cumulative.pf

git commit -m "$(cat <<'EOF'
refactor(state): flatten surfacewater cohort sub-records into surfacewater_state_t

Replaces three nested cohort sub-records (surfacewater_intermediate_t,
surfacewater_drainage_cumulative_t, surfacewater_reservoir_cumulative_t
from ADR 0033 Phase A + correction) with flat fields on
surfacewater_state_t and three named reset procedures:
  - reset_intermediate          (flzerointr gate)
  - reset_cumulative_drainage   (flzerocumu gate; fldrain partition)
  - reset_cumulative_reservoir  (flzerocumu gate; flSurfaceWater partition)
All ~57 call-site references across drainage.f90, waterbalance.f90,
surfacewater.f90, swapoutput.f90, swap_csv_output.f90, soilgrid.f90,
management_soil.f90, surfacewater_init.f90, and the unit tests
rewritten to the shallower paths. Allocation block paths-only rewritten;
allocation consolidation deferred to follow-up arc. check-full
byte-identical.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: ADR for the arc

**Files:**
- Create: `docs/adr/0042-flatten-reset-cohorts.md`

### Step 4.1: Write ADR

- [ ] **Create `docs/adr/0042-flatten-reset-cohorts.md` with this content:**

```markdown
---
title: "ADR 0042 — Flatten Reset Cohorts Across State Subsystems"
date: 2026-05-12
status: accepted
supersedes-portion-of: ADR 0033
---

# ADR 0042: Flatten Reset Cohorts Across State Subsystems

**Status:** accepted
**Date:** 2026-05-12
**Branch:** `development`

## Context

ADR 0033 introduced nested cohort sub-records (`<subsystem>_intermediate_t`,
`<subsystem>_cumulative_t`, and surfacewater's partitioned cumulative
cohorts) on three subsystem state records — surfacewater, solute, soilwater
— to make reset cadence first-class in the type system. The cohort pattern
produced two real benefits: (a) co-location of fields that zero together,
(b) type-bound `reset()` syntactically compact.

In review, the wrapper bought less than expected:

- Activity gates (`flzerointr`, `flzerocumu`, `flDayStart`, `fldrain`,
  `flSurfaceWater`) remain enforced at call sites; the cohort type does
  not enforce them.
- Field-name prefixes (`i*` / `c*`) already encode cohort role; the
  wrapper is redundant with the naming convention.
- The cohort-vs-flat boundary on the parent is inconsistent: scalar
  non-reset fields (`wls`, `imper`) live flat while reset fields hide
  one level deeper.
- Access paths in hot compute loops (waterbalance, soilgrid) traverse
  three levels.
- Allocation traverses cohort paths separately from flat-field paths,
  fragmenting an otherwise unified init step.

The asymmetric reset semantics that ADR 0033 wanted to capture
(surfacewater's partitioned cumulative; soilwater's per-day subset
reset) can be expressed equally well by **named reset procedures on the
flat parent type**. Each asymmetry becomes its own procedure; the
procedure name carries the gate. The cohort form was, in effect, "two
reset procedures wearing a type-shaped costume" (surfacewater) or "a
second procedure on one cohort" (soilwater).

## Decision

Flatten the cohort sub-records into the parent state record. Move
`reset()` from cohort-type-bound to parent-type-bound procedures, named
to preserve the gate/cadence semantics.

Per subsystem:

- `surfacewater_state_t`:
  - `reset_intermediate()` (flzerointr)
  - `reset_cumulative_drainage()` (flzerocumu + fldrain partition; owner: drainage)
  - `reset_cumulative_reservoir()` (flzerocumu + flSurfaceWater partition; owner: surface-water)

- `solute_state_t`:
  - `reset_intermediate()` (flzerointr)
  - `reset_cumulative()` (flzerocumu)

- `soilwater_state_t`:
  - `reset_intermediate()` (flzerointr; subsumes per-day subset)
  - `reset_intermediate_per_day()` (flDayStart; per-day subset only)
  - `reset_cumulative()` (flzerocumu)

Section-header comments inside the flat type group fields by cadence
(`! === intermediate (reset_intermediate / gate: flzerointr) ===`,
etc.), restoring the type-signature self-documentation that ADR 0033
§82-83 valued.

## Implementation

Three per-subsystem commits in order solute → soilwater → surfacewater.
Each commit fully verified (build + pFUnit + check-full + audit grep)
before the next.

- Solute commit: `<hash>` (Task 1)
- Soilwater commit: `<hash>` (Task 2)
- Surfacewater commit: `<hash>` (Task 3)

The cohort sub-record types are removed from the public surface of each
state module. All call-site paths collapse one level shallower
(`state%X%cohort%field` → `state%X%field`). Reset call sites rewrite to
the new procedure names. ASSOCIATE block RHS targets rewrite mechanically
under the same path-collapse rule.

## Out of scope

- **Allocation consolidation.** Duplicate inline allocation blocks
  (`surfacewater_init.f90:119-141` and `drainage.f90:445-467`, and the
  fragmented allocator paths in `soilwater_state.f90:soilwater_init`)
  are pre-existing. Deferred to a separate follow-up arc.
- ADR 0033 itself is left intact as the historical record; this ADR
  references and supersedes only the cohort-pattern decision.

## Consequences

(+) Call-site paths one level shallower; less friction in hot compute
    loops.
(+) Subset / asymmetric resets are explicit named procedures rather
    than conventions buried inside cohort types.
(+) Removes the cohort-vs-flat inconsistency on the parent (some
    fields nested, some not).
(+) Field-name prefix convention (`i*` / `c*`) is no longer redundant
    with cohort wrappers.

(–) Adding a new reset-eligible field requires editing the field
    declaration and the appropriate reset procedure body — same
    cognitive load as before, but the procedure body is no longer
    co-located with the field inside a sub-type. Mitigation:
    section-header comments inside the type group fields by cadence.
(–) Type-bound `cohort%reset()` syntactic compactness is lost. The
    parent now carries multiple named reset procedures; names are
    slightly longer but spell out the gate.
(–) The ADR 0033 §82-83 self-documentation property — `state%X%cumulative%Y`
    instantly tells the reader "this is a cumulative balance variable"
    — weakens. A grep against `state%X%fieldname` alone no longer
    reveals cadence. Section-header comments partially restore this.

## References

- ADR 0030 — state-migration surfacewater pilot
- ADR 0032 — state-migration solute
- ADR 0033 — cumulative reset cohorts (the pattern this arc walks back)
- ADR 0038 — state-migration soilwater core
- Design spec: `docs/superpowers/specs/2026-05-12-flatten-reset-cohorts-design.md`
- Plan: `docs/superpowers/plans/2026-05-12-flatten-reset-cohorts.md`
```

Fill in the three commit hashes from `git log --oneline -n 5` after Tasks 1-3 ship.

### Step 4.2: Fill in commit hashes

- [ ] **Read commit hashes:**

```bash
cd /home/zawadzkim/Code/swap
git log --oneline -n 6
```

Replace the three `<hash>` placeholders in `docs/adr/0042-flatten-reset-cohorts.md` with the actual short hashes from Tasks 1, 2, 3.

### Step 4.3: Commit the ADR

- [ ] **Stage and commit:**

```bash
cd /home/zawadzkim/Code/swap
git add docs/adr/0042-flatten-reset-cohorts.md

git commit -m "$(cat <<'EOF'
docs(adr): ADR 0042 — flatten reset cohorts across state subsystems

Records the 2026-05-12 arc that replaced nested cohort sub-records
(ADR 0033) in surfacewater_state, solute_state, and soilwater_state
with flat fields plus named reset procedures on the parent types.
References Tasks 1-3 commit hashes.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Plan Self-Review

(Run after writing the complete plan; fix issues inline before handoff.)

**Spec coverage check** — each major spec section maps to a task:
- Spec §"Architecture / What goes / What stays" → Tasks 1.1, 2.1, 3.1 (state-module rewrites)
- Spec §"Per-subsystem reset interface" → procedure names and bodies in 1.1, 2.1, 3.1
- Spec §"Field grouping" → section-header comments in 1.1, 2.1, 3.1
- Spec §"Call-site migration / Path rewrites / Reset-call rewrites / ASSOCIATE blocks / Audit greps" → Tasks 1.4-1.8, 2.3-2.5, 3.4-3.6
- Spec §"Tests" → Tasks 1.2-1.3, 2.2, 3.2-3.3
- Spec §"Verification gates" → check-full step in each task
- Spec §"Sequencing" → Tasks 1 → 2 → 3 ordered
- Spec §"Scope boundaries (out of scope: allocation, ADR 0033 edits)" → explicitly noted in Step 2.1 (item 7), Step 3.4, and Task 4 ADR text
- Spec §"Write a new ADR" → Task 4

No spec section uncovered.

**Placeholder check:** ADR placeholder `<hash>` is filled in by Step 4.2 (concrete procedure). No "TBD" / "implement later" / "similar to" placeholders elsewhere.

**Type consistency check:**
- `solute_state_t` reset procedures: `reset_intermediate`, `reset_cumulative` — same names in module (Step 1.1), tests (Step 1.3), and call sites (Step 1.4). ✓
- `soilwater_state_t` reset procedures: `reset_intermediate`, `reset_intermediate_per_day`, `reset_cumulative` — same names in module (Step 2.1), tests (Step 2.2), call sites (Step 2.3). ✓
- `surfacewater_state_t` reset procedures: `reset_intermediate`, `reset_cumulative_drainage`, `reset_cumulative_reservoir` — same names in module (Step 3.1), tests (Step 3.3), call sites (Step 3.4). ✓
- `samini`, `sampro` rebase preserved inline; no name change. ✓
