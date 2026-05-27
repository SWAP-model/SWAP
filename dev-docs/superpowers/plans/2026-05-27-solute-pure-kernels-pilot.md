# Solute Pure-Kernels Pilot Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extract `solute_step`'s physics into `pure`, state-free kernels and promote the five derived coefficient arrays to typed state — making the math unit-testable, the orchestration narrative, and fixing a latent uninitialised-read bug, all byte-identical on the current regression suite.

**Architecture:** A new `solute_kernels_mod` (`src/solute/solute_kernels.f90`) holds `pure`/`elemental` functions with explicit args and no `swap_state_t` dependency (mirroring `interception.f90`/`et.f90`). `solute_mod` `use`s it. Five time-invariant coefficient arrays move from discarded locals into `solute_state_t`, computed once in `solute_seed`, read in `solute_step`.

**Tech Stack:** Fortran (gfortran), Meson/ninja build, pFUnit unit tests, Python regression harness, `pixi` task runner.

**Design spec:** `dev-docs/superpowers/specs/2026-05-27-solute-pure-kernels-pilot-design.md`

---

## Background the implementer must know

- **Build is an explicit file list, not a glob.** New `.f90` files must be added to BOTH `meson.build` (production binary) and `tests/unit/meson.build` (pFUnit binary).
- **pFUnit registration is two-place.** A new `.pf` must be added to the `pf_files` list in `tests/unit/meson.build` AND get an `ADD_TEST_SUITE(<name>_suite)` line in `tests/unit/testSuites.inc`. An unregistered `.pf` compiles but never runs — always confirm the `OK (N tests)` count rises.
- **State-schema changes need a clean rebuild.** After editing `src/state/solute_state.f90`, run `rm -rf builddir` before building; incremental Meson does not propagate `.mod` deps across the swap_modern↔swap_legacy boundary. (Only Task 2 touches the schema.)
- **Byte-identical gate, every task.** End each task with `pixi run -e test check-fast` (build + pFUnit + hupselbrook/surfacewater/salinitystress/grassgrowth). `salinitystress` is the only solute-active fast case and is the byte-identical witness for this work.
- **The five coefficients are time-invariant** (`bdens*kf`, `bdens*kf*cref`, `bdens*kfsat+poros`, `ddif/thetsl²`, `decpot*fdepth`) — soil/solute properties, correctly computed once in seed.
- **Reference commands:**
  - Build: `pixi run -e test build-linux`
  - pFUnit only: `pixi run -e test test-pfunit`
  - One regression case: `pixi run -e test regression salinitystress`
  - Fast gate: `pixi run -e test check-fast`

---

## File Structure

- **Create** `src/solute/solute_kernels.f90` — `module solute_kernels_mod`; `pure`/`elemental` physics functions only. One responsibility: solute math with explicit args.
- **Create** `tests/unit/solute/test_solute_kernels.pf` — unit tests for every kernel function.
- **Modify** `src/state/solute_state.f90` — add 5 per-node coefficient fields + allocate them in `solute_state_init` (Task 2 only).
- **Modify** `src/solute/solute.f90` — `solute_seed` fills coefficients via kernels into state; `solute_step` reads from state and calls the Freundlich/decomposition kernels; dead locals deleted.
- **Modify** `meson.build` — register `src/solute/solute_kernels.f90`.
- **Modify** `tests/unit/meson.build` — register the kernels source and the `.pf`.
- **Modify** `tests/unit/testSuites.inc` — register `test_solute_kernels_suite`.

---

## Task 1: Create `solute_kernels_mod` with coefficient functions (TDD, no behavior change yet)

Creates the four coefficient `elemental` functions and their tests. `solute.f90` is NOT touched in this task — the module lands green and standalone first.

**Files:**
- Create: `src/solute/solute_kernels.f90`
- Create: `tests/unit/solute/test_solute_kernels.pf`
- Modify: `meson.build` (Solute section, ~line 170)
- Modify: `tests/unit/meson.build` (pfunit source list ~line 180; `pf_files` list ~line 207)
- Modify: `tests/unit/testSuites.inc`

- [ ] **Step 1: Write the failing test**

Create `tests/unit/solute/test_solute_kernels.pf`:

```fortran
! Unit tests for solute_kernels_mod — pure/elemental solute physics.
@test
subroutine test_bdenskf_coeff()
   use funit
   use iso_fortran_env, only: real64
   use solute_kernels_mod, only: bdenskf_coeff
   @assertEqual(6.0_real64, bdenskf_coeff(2.0_real64, 3.0_real64), 1.0e-12_real64)
   ! kf=0 (conservative salt): correct value is 0 — the masked case today.
   @assertEqual(0.0_real64, bdenskf_coeff(1.5_real64, 0.0_real64), 1.0e-12_real64)
end subroutine test_bdenskf_coeff

@test
subroutine test_bdenskfsatporos_coeff()
   use funit
   use iso_fortran_env, only: real64
   use solute_kernels_mod, only: bdenskfsatporos_coeff
   @assertEqual(3.4_real64, bdenskfsatporos_coeff(1.5_real64, 2.0_real64, 0.4_real64), 1.0e-12_real64)
   ! kfsat=0: reduces to porosity — a non-zero value the current code never exercises.
   @assertEqual(0.35_real64, bdenskfsatporos_coeff(1.3_real64, 0.0_real64, 0.35_real64), 1.0e-12_real64)
end subroutine test_bdenskfsatporos_coeff

@test
subroutine test_ddiffwcs_coeff()
   use funit
   use iso_fortran_env, only: real64
   use solute_kernels_mod, only: ddiffwcs_coeff
   ! 0.5 / 0.4**2 = 0.5 / 0.16 = 3.125
   @assertEqual(3.125_real64, ddiffwcs_coeff(0.5_real64, 0.4_real64), 1.0e-12_real64)
   ! ddif=0 (no molecular diffusion): correct value is 0 — masked case today.
   @assertEqual(0.0_real64, ddiffwcs_coeff(0.0_real64, 0.4_real64), 1.0e-12_real64)
end subroutine test_ddiffwcs_coeff

@test
subroutine test_decpotfdepth_coeff()
   use funit
   use iso_fortran_env, only: real64
   use solute_kernels_mod, only: decpotfdepth_coeff
   @assertEqual(0.08_real64, decpotfdepth_coeff(0.1_real64, 0.8_real64), 1.0e-12_real64)
   ! decpot=0 (no decay): correct value is 0 — masked case today.
   @assertEqual(0.0_real64, decpotfdepth_coeff(0.0_real64, 0.8_real64), 1.0e-12_real64)
end subroutine test_decpotfdepth_coeff
```

- [ ] **Step 2: Register the test, then build to verify it fails**

In `tests/unit/meson.build`, add to the `pf_files` list (after the existing entries, e.g. near the io/config block — placement within the list does not matter):

```
        'solute/test_solute_kernels.pf',
```

In `tests/unit/meson.build`, add to the pfunit source list (immediately before the `'../../src/solute/solute.f90',` line ~180):

```
        '../../src/solute/solute_kernels.f90',
```

In `tests/unit/testSuites.inc`, add at the end:

```
ADD_TEST_SUITE(test_solute_kernels_suite)
```

Run: `pixi run -e test test-pfunit`
Expected: FAIL — compile error, `Cannot open module file 'solute_kernels_mod.mod'` or unresolved `bdenskf_coeff` (the module/source does not exist yet).

- [ ] **Step 3: Create the kernels module**

Create `src/solute/solute_kernels.f90`:

```fortran
!> Pure solute physics kernels — explicit args, no swap_state_t dependency.
!> Pilot extraction from solute_step (see spec 2026-05-27). Mirrors the
!> interception.f90 / et.f90 "pure function — explicit args" pattern.
module solute_kernels_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   public :: bdenskf_coeff, bdenskfsatporos_coeff, ddiffwcs_coeff, decpotfdepth_coeff

contains

   !> Freundlich bulk sorption coefficient: bdens * kf  [-].
   elemental function bdenskf_coeff(bdens, kf) result(v)
      real(real64), intent(in) :: bdens   ! dry soil bulk density (g/cm3)
      real(real64), intent(in) :: kf      ! Freundlich coefficient (cm3/g)
      real(real64)             :: v
      v = bdens*kf
   end function bdenskf_coeff

   !> Saturated-zone sorption + porosity: bdens * kfsat + poros  [-].
   elemental function bdenskfsatporos_coeff(bdens, kfsat, poros) result(v)
      real(real64), intent(in) :: bdens   ! dry soil bulk density (g/cm3)
      real(real64), intent(in) :: kfsat   ! saturated-zone Freundlich coeff (cm3/g)
      real(real64), intent(in) :: poros   ! aquifer porosity (-)
      real(real64)             :: v
      v = bdens*kfsat + poros
   end function bdenskfsatporos_coeff

   !> Tortuosity-scaled diffusion base: ddif / thetsl**2  (cm2/d).
   elemental function ddiffwcs_coeff(ddif, thetsl) result(v)
      real(real64), intent(in) :: ddif    ! molecular diffusion coefficient (cm2/d)
      real(real64), intent(in) :: thetsl  ! saturated water content (-)
      real(real64)             :: v
      v = ddif / (thetsl**2)
   end function ddiffwcs_coeff

   !> Depth-weighted potential decomposition: decpot * fdepth  (1/d).
   elemental function decpotfdepth_coeff(decpot, fdepth) result(v)
      real(real64), intent(in) :: decpot  ! potential decomposition rate (1/d)
      real(real64), intent(in) :: fdepth  ! depth-decomposition factor (-)
      real(real64)             :: v
      v = decpot*fdepth
   end function decpotfdepth_coeff

end module solute_kernels_mod
```

In `meson.build`, in the `# Solute` section, add immediately before `'src/solute/solute.f90',`:

```
    'src/solute/solute_kernels.f90',
```

- [ ] **Step 4: Build and run tests to verify they pass**

Run: `pixi run -e test test-pfunit`
Expected: PASS — the suite count rises by 4 (the four `@test` subroutines); look for `OK (` with a higher total than before this task.

- [ ] **Step 5: Commit**

```bash
git add src/solute/solute_kernels.f90 tests/unit/solute/test_solute_kernels.pf \
        meson.build tests/unit/meson.build tests/unit/testSuites.inc
git commit -m "feat(solute): add solute_kernels_mod coefficient functions + tests

Pure/elemental coefficient kernels (bdenskf/bdenskfsatporos/ddiffwcs/
decpotfdepth) with golden-value + masked-zero unit tests. Not yet wired
into solute.f90.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 2: Promote coefficients to state, wire seed/step, delete dead locals (the bug fix)

Moves the five coefficient arrays into `solute_state_t`, computes them once in `solute_seed` via the Task-1 kernels, and makes `solute_step` read them from state instead of from uninitialised locals. **State-schema change → clean rebuild required.**

**Files:**
- Modify: `src/state/solute_state.f90` (field block ~line 31; `solute_state_init` ~line 135)
- Modify: `src/solute/solute.f90` (`solute_seed` decls ~27-28 and loop ~56-66; `solute_step` decls ~91-92 and reads at 134/197/234/239/252-257)
- Modify: `tests/unit/state/test_solute_state.pf` (add an init-allocation test)

- [ ] **Step 1: Write the failing test (allocation of new coefficient fields)**

Append to `tests/unit/state/test_solute_state.pf` (before the final `end module`/last line if present; otherwise at end of file):

```fortran
@test
subroutine test_solute_state_init_allocates_coeffs()
   use funit
   use iso_fortran_env, only: real64
   use solute_state_mod,  only: solute_state_t
   use solute_config_mod, only: solute_config_t
   type(solute_state_t)  :: sl
   type(solute_config_t) :: cfg
   integer, parameter :: numnod = 5

   call sl%init(cfg, numnod)

   @assertTrue(allocated(sl%bdenskf))
   @assertTrue(allocated(sl%bdenskfcref))
   @assertTrue(allocated(sl%bdenskfsatporos))
   @assertTrue(allocated(sl%ddiffwcs))
   @assertTrue(allocated(sl%decpotfdepth))
   @assertEqual(numnod, size(sl%bdenskf))
   @assertEqual(numnod, size(sl%decpotfdepth))
end subroutine test_solute_state_init_allocates_coeffs
```

- [ ] **Step 2: Build to verify it fails**

Run: `pixi run -e test test-pfunit`
Expected: FAIL — compile error: `sl%bdenskf` etc. are not members of `solute_state_t`.

- [ ] **Step 3: Add the five fields to `solute_state_t`**

In `src/state/solute_state.f90`, in the `! === per-node arrays` block (immediately after the `cmsy` line ~31-32), add:

```fortran
      ! Derived time-invariant coefficients (filled once in solute_seed).
      real(real64), allocatable :: bdenskf(:)          !! bdens*kf per node (-)
      real(real64), allocatable :: bdenskfcref(:)      !! bdenskf*cref per node (M/L3)
      real(real64), allocatable :: bdenskfsatporos(:)  !! bdens*kfsat+poros per node (-)
      real(real64), allocatable :: ddiffwcs(:)         !! ddif/thetsl**2 per node (cm2/d)
      real(real64), allocatable :: decpotfdepth(:)     !! decpot*fdepth per node (1/d)
```

- [ ] **Step 4: Allocate the fields in `solute_state_init`**

In `src/state/solute_state.f90`, in `solute_state_init`, immediately after the existing `allocate(self%cmsy(n))` line (~136), add:

```fortran
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
```

- [ ] **Step 5: Clean rebuild and run the state test to verify it passes**

Run: `rm -rf builddir && pixi run -e test test-pfunit`
Expected: PASS — `test_solute_state_init_allocates_coeffs` green; total `OK (N tests)` count up by 1 vs Task 1.

- [ ] **Step 6: Rewrite `solute_seed` to fill coefficients via kernels**

In `src/solute/solute.f90`, add the kernels import to `solute_seed` (after the existing `use` lines, ~14):

```fortran
      use solute_kernels_mod,    only: bdenskf_coeff, bdenskfsatporos_coeff, &
                                       ddiffwcs_coeff, decpotfdepth_coeff
```

Replace the `! Derived solute concentrations.` loop body (current lines ~56-66) with:

```fortran
         ! Derived solute concentrations + time-invariant coefficients.
         sol%samini = 0.0d0
         do i = 1, mesh%numnod
            sol%bdenskf(i)         = bdenskf_coeff(soil%bdens(mesh%layer(i)), sol%kf(mesh%layer(i)))
            sol%bdenskfcref(i)     = sol%bdenskf(i)*sol%cref
            sol%bdenskfsatporos(i) = bdenskfsatporos_coeff(soil%bdens(mesh%layer(i)), sol%kfsat, sol%poros)
            sol%cmsy(i)            = soil%theta(i)*sol%cml(i) +                          &
                                     sol%bdenskfcref(i)*(sol%cml(i)/sol%cref)**sol%frexp
            sol%samini             = sol%samini + sol%cmsy(i) * mesh%dz(i)
            sol%ddiffwcs(i)        = ddiffwcs_coeff(sol%ddif, soil%thetsl(mesh%layer(i)))
            sol%decpotfdepth(i)    = decpotfdepth_coeff(sol%decpot(mesh%layer(i)), sol%fdepth(mesh%layer(i)))
         end do
         sol%sampro = sol%samini
```

Then delete the now-unused local array declaration line from `solute_seed` (current ~27-28):

```fortran
      real(8), dimension(macp) :: thetav, diffus, dispr1, vpore2, ddiffwcs, &
                                  bdenskf, bdenskfcref, bdenskfsatporos, decpotfdepth
```

(None of `thetav/diffus/dispr1/vpore2/ddiffwcs/bdenskf/bdenskfcref/bdenskfsatporos/decpotfdepth` are used elsewhere in `solute_seed`, so remove the whole line. If the compiler reports any of these still referenced in `solute_seed`, keep only those — but it should not.)

- [ ] **Step 7: Update `solute_step` to read coefficients from state and drop dead locals**

In `src/solute/solute.f90`, in `solute_step`, change the local array declaration line (current ~91-92) to keep only the genuinely-local working arrays:

```fortran
      real(8), dimension(macp) :: thetav, diffus, dispr1, vpore2
```

Then update the five read sites to use `sol%`:
- Line ~134: `diffus(i) = ddiffwcs(i) * thetav(i)**2.33d0` → `diffus(i) = sol%ddiffwcs(i) * thetav(i)**2.33d0`
- Line ~197: `decact*bdenskfcref(i)*(...)` → `decact*sol%bdenskfcref(i)*(...)`
- Line ~195: `decact = decpotfdepth(i) * ftemp * ftheta` → `decact = sol%decpotfdepth(i) * ftemp * ftheta`
- Line ~234: `sol%cml(i) = sol%cmsy(i) / (soil%theta(i) + bdenskf(i))` → `... + sol%bdenskf(i))`
- Line ~239: `dummy = bdenskf(i)*(...)` → `dummy = sol%bdenskf(i)*(...)`
- Lines ~252,254,256,257: every `bdenskfsatporos(i)` → `sol%bdenskfsatporos(i)`

- [ ] **Step 8: Clean rebuild + byte-identical regression gate**

Run: `rm -rf builddir && pixi run -e test check-fast`
Expected: PASS — all four cases, including `salinitystress`, report "regression ok (annual stats match fixture)". (Correct coefficient values equal the previously-masked zeros for this conservative-salt case, so output is byte-identical.) pFUnit `OK (N tests)` unchanged from Step 5.

- [ ] **Step 9: Commit**

```bash
git add src/state/solute_state.f90 src/solute/solute.f90 tests/unit/state/test_solute_state.pf
git commit -m "fix(solute): promote derived coefficients to state%solute

solute_step read five derived coefficient arrays (ddiffwcs/bdenskf/
bdenskfcref/bdenskfsatporos/decpotfdepth) that were only computed (as
discarded locals) in solute_seed — an uninitialised read introduced by
the seed/step split (52ac7ff), masked by the conservative-salt
salinitystress case. Promote them to typed per-node state, filled once
in solute_seed via solute_kernels_mod. Byte-identical on the suite.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 3: Extract the Freundlich isotherm inversion

Replaces the inline cmsy→cml recovery solver (current lines ~227-244) with a tested `pure function`.

**Files:**
- Modify: `src/solute/solute_kernels.f90` (add `solute_cml_from_cmsy`)
- Modify: `tests/unit/solute/test_solute_kernels.pf` (add tests)
- Modify: `src/solute/solute.f90` (`solute_step`)

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/solute/test_solute_kernels.pf`:

```fortran
@test
subroutine test_cml_from_cmsy_linear_branch()
   use funit
   use iso_fortran_env, only: real64
   use solute_kernels_mod, only: solute_cml_from_cmsy
   real(real64) :: cml
   ! |frexp-1|<0.001 → linear: cml = cmsy/(theta+bdenskf) = 4.0/(0.3+0.5) = 5.0
   cml = solute_cml_from_cmsy(4.0_real64, 0.3_real64, 0.5_real64, 1.0_real64, 1.0_real64, 0.0_real64)
   @assertEqual(5.0_real64, cml, 1.0e-12_real64)
   ! frexp=1.0005 still inside the linear threshold
   cml = solute_cml_from_cmsy(4.0_real64, 0.3_real64, 0.5_real64, 1.0005_real64, 1.0_real64, 0.0_real64)
   @assertEqual(5.0_real64, cml, 1.0e-12_real64)
end subroutine test_cml_from_cmsy_linear_branch

@test
subroutine test_cml_from_cmsy_nonlinear_satisfies_fixedpoint()
   use funit
   use iso_fortran_env, only: real64
   use solute_kernels_mod, only: solute_cml_from_cmsy
   real(real64) :: cml, theta, bdenskf, frexp, cref, cmsy, residual
   theta = 0.25_real64; bdenskf = 1.2_real64; frexp = 0.7_real64
   cref = 1.0_real64;   cmsy = 3.0_real64
   cml = solute_cml_from_cmsy(cmsy, theta, bdenskf, frexp, cref, 1.0_real64)
   ! Converged cml must satisfy cml = cmsy / (theta + bdenskf*(cml/cref)**(frexp-1))
   residual = cml - cmsy/(theta + bdenskf*(cml/cref)**(frexp-1.0_real64))
   @assertEqual(0.0_real64, residual, 1.0e-3_real64)
   @assertTrue(cml > 0.0_real64)
end subroutine test_cml_from_cmsy_nonlinear_satisfies_fixedpoint
```

- [ ] **Step 2: Build to verify it fails**

Run: `pixi run -e test test-pfunit`
Expected: FAIL — unresolved `solute_cml_from_cmsy`.

- [ ] **Step 3: Add the function to `solute_kernels_mod`**

In `src/solute/solute_kernels.f90`, add `solute_cml_from_cmsy` to the `public` list and implement it in `contains`:

```fortran
   !> Recover mobile concentration cml from total cmsy via the Freundlich
   !> isotherm. Linear shortcut when frexp ~ 1; otherwise fixed-point iterate
   !> seeded from cml_guess. Caller handles the cmsy < vsmall zeroing.
   pure function solute_cml_from_cmsy(cmsy, theta, bdenskf, frexp, cref, cml_guess) result(cml)
      real(real64), intent(in) :: cmsy       ! total (dissolved+adsorbed) conc (M/L3 soil)
      real(real64), intent(in) :: theta      ! volumetric water content (-)
      real(real64), intent(in) :: bdenskf    ! bdens*kf (-)
      real(real64), intent(in) :: frexp      ! Freundlich exponent (-)
      real(real64), intent(in) :: cref       ! reference concentration (M/L3)
      real(real64), intent(in) :: cml_guess  ! previous cml, iteration seed (M/L3)
      real(real64)             :: cml

      real(real64), parameter :: rer    = 1.0d-3
      real(real64), parameter :: vsmall = 1.0d-15
      real(real64) :: old, dummy
      logical      :: differ

      if (abs(frexp - 1.0d0) .lt. 0.001d0) then
         cml = cmsy / (theta + bdenskf)
      else
         cml = cml_guess
         if (cml .lt. vsmall) cml = vsmall
         differ = .true.
         do while (differ)
            old   = cml
            dummy = bdenskf*(cml/cref)**(frexp - 1.0d0)
            cml   = cmsy/(theta + dummy)
            if (abs(cml - old) .lt. rer*cml) differ = .false.
         end do
      end if
   end function solute_cml_from_cmsy
```

- [ ] **Step 4: Build and run tests to verify they pass**

Run: `pixi run -e test test-pfunit`
Expected: PASS — `OK (N tests)` up by 2 vs Task 2.

- [ ] **Step 5: Wire it into `solute_step`**

In `src/solute/solute.f90`, add to the `solute_step` kernels import:

```fortran
      use solute_kernels_mod,    only: solute_cml_from_cmsy
```

Replace the recovery block (current lines ~227-244) with:

```fortran
               ! Iterate to recover cml from cmsy with the Freundlich isotherm.
               if (sol%cmsy(i) .lt. vsmall) then
                  sol%cmsy(i) = 0.0d0
                  sol%cml(i)  = 0.0d0
               else
                  sol%cml(i) = solute_cml_from_cmsy(sol%cmsy(i), soil%theta(i), &
                                  sol%bdenskf(i), sol%frexp, sol%cref, sol%cml(i))
               end if
```

The `differ`, `old`, `dummy` locals and the `rer` parameter are no longer used in `solute_step` — remove `differ` from its declarations and delete the `rer` parameter line. Keep the `vsmall` parameter (still used by the guard above) and keep `dummy`/`old` only if still referenced elsewhere in `solute_step` (they are: `dummy` is used in the dtsolu loop at ~143; `old` is not — remove `old`). Verify by compiler.

- [ ] **Step 6: Byte-identical regression gate**

Run: `pixi run -e test check-fast`
Expected: PASS — all four cases byte-identical; `salinitystress` "regression ok". (No clean rebuild needed — no schema change this task.)

- [ ] **Step 7: Commit**

```bash
git add src/solute/solute_kernels.f90 tests/unit/solute/test_solute_kernels.pf src/solute/solute.f90
git commit -m "refactor(solute): extract Freundlich inversion to solute_cml_from_cmsy

Pure fixed-point isotherm solver, unit-tested for the linear branch and
fixed-point convergence. solute_step now calls it per node; byte-identical.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 4: Extract the decomposition factors

Replaces the inline temperature/moisture decomposition formulas (current lines ~185-197) with tested `pure` functions.

**Files:**
- Modify: `src/solute/solute_kernels.f90` (add `solute_ftemp`, `solute_ftheta`, `solute_decomp_ctrans`)
- Modify: `tests/unit/solute/test_solute_kernels.pf` (add tests)
- Modify: `src/solute/solute.f90` (`solute_step`)

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/solute/test_solute_kernels.pf`:

```fortran
@test
subroutine test_solute_ftemp()
   use funit
   use iso_fortran_env, only: real64
   use solute_kernels_mod, only: solute_ftemp
   ! tsoil=20 → exp(gampar*0)=1
   @assertEqual(1.0_real64, solute_ftemp(20.0_real64, 0.1_real64, .true.), 1.0e-12_real64)
   ! tsoil>=35 capped at exp(gampar*15)
   @assertEqual(exp(0.1_real64*15.0_real64), solute_ftemp(40.0_real64, 0.1_real64, .true.), 1.0e-12_real64)
   ! temperature switch off → 0
   @assertEqual(0.0_real64, solute_ftemp(10.0_real64, 0.0693_real64, .false.), 1.0e-12_real64)
end subroutine test_solute_ftemp

@test
subroutine test_solute_ftheta()
   use funit
   use iso_fortran_env, only: real64
   use solute_kernels_mod, only: solute_ftheta
   ! theta=rtheta → ratio 1 → 1
   @assertEqual(1.0_real64, solute_ftheta(0.3_real64, 0.3_real64, 0.7_real64), 1.0e-12_real64)
   ! half of rtheta, bexp=1 → 0.5
   @assertEqual(0.5_real64, solute_ftheta(0.15_real64, 0.3_real64, 1.0_real64), 1.0e-12_real64)
   ! above rtheta → clamped at 1
   @assertEqual(1.0_real64, solute_ftheta(0.6_real64, 0.3_real64, 1.0_real64), 1.0e-12_real64)
end subroutine test_solute_ftheta

@test
subroutine test_solute_decomp_ctrans()
   use funit
   use iso_fortran_env, only: real64
   use solute_kernels_mod, only: solute_decomp_ctrans
   ! bdenskfcref=0: ctrans = decact*theta*cml = 0.5*0.3*2.0 = 0.3
   @assertEqual(0.3_real64, &
      solute_decomp_ctrans(0.5_real64, 0.3_real64, 2.0_real64, 0.0_real64, 1.0_real64, 1.0_real64), &
      1.0e-12_real64)
   ! bdenskfcref=0.4, frexp=1, cref=1: + 0.5*0.4*(2.0)**1 = 0.4 → total 0.7
   @assertEqual(0.7_real64, &
      solute_decomp_ctrans(0.5_real64, 0.3_real64, 2.0_real64, 0.4_real64, 1.0_real64, 1.0_real64), &
      1.0e-12_real64)
end subroutine test_solute_decomp_ctrans
```

- [ ] **Step 2: Build to verify it fails**

Run: `pixi run -e test test-pfunit`
Expected: FAIL — unresolved `solute_ftemp` / `solute_ftheta` / `solute_decomp_ctrans`.

- [ ] **Step 3: Add the functions to `solute_kernels_mod`**

In `src/solute/solute_kernels.f90`, add the three names to `public` and implement:

```fortran
   !> Temperature reduction factor for decomposition. Capped above 35 degC;
   !> zero when the temperature switch is off.
   pure function solute_ftemp(tsoil, gampar, fl_temperature) result(ftemp)
      real(real64), intent(in) :: tsoil           ! soil temperature (degC)
      real(real64), intent(in) :: gampar          ! temperature coefficient (/C)
      logical,      intent(in) :: fl_temperature  ! temperature simulation on/off
      real(real64)             :: ftemp
      if (fl_temperature) then
         if (tsoil .lt. 35.0d0) then
            ftemp = exp(gampar*(tsoil - 20.0d0))
         else
            ftemp = exp(gampar*15.0d0)
         end if
      else
         ftemp = 0.0d0
      end if
   end function solute_ftemp

   !> Moisture reduction factor for decomposition, clamped to 1.
   pure function solute_ftheta(theta, rtheta, bexp) result(ftheta)
      real(real64), intent(in) :: theta   ! volumetric water content (-)
      real(real64), intent(in) :: rtheta  ! reference moisture content (-)
      real(real64), intent(in) :: bexp    ! moisture-decomposition exponent (-)
      real(real64)             :: ftheta
      ftheta = min(1.0d0, (theta/rtheta)**bexp)
   end function solute_ftheta

   !> Solute transformation (decomposition) rate per node.
   pure function solute_decomp_ctrans(decact, theta, cml, bdenskfcref, cref, frexp) result(ctrans)
      real(real64), intent(in) :: decact       ! actual decomposition rate (1/d)
      real(real64), intent(in) :: theta        ! volumetric water content (-)
      real(real64), intent(in) :: cml          ! mobile concentration (M/L3)
      real(real64), intent(in) :: bdenskfcref  ! bdens*kf*cref (M/L3)
      real(real64), intent(in) :: cref         ! reference concentration (M/L3)
      real(real64), intent(in) :: frexp        ! Freundlich exponent (-)
      real(real64)             :: ctrans
      ctrans = decact*theta*cml + decact*bdenskfcref*((cml/cref)**frexp)
   end function solute_decomp_ctrans
```

- [ ] **Step 4: Build and run tests to verify they pass**

Run: `pixi run -e test test-pfunit`
Expected: PASS — `OK (N tests)` up by 3 vs Task 3.

- [ ] **Step 5: Wire into `solute_step`**

In `src/solute/solute.f90`, add to the `solute_step` kernels import:

```fortran
      use solute_kernels_mod,    only: solute_ftemp, solute_ftheta, solute_decomp_ctrans
```

Replace the decomposition block (current lines ~184-197) with:

```fortran
               ! Solute decomposition.
               ftemp  = solute_ftemp(heat%tsoil(i), sol%gampar, time%flTemperature)
               ftheta = solute_ftheta(soil%theta(i), sol%rtheta, sol%bexp)
               decact = sol%decpotfdepth(i) * ftemp * ftheta
               ctrans = solute_decomp_ctrans(decact, soil%theta(i), sol%cml(i), &
                                             sol%bdenskfcref(i), sol%cref, sol%frexp)
```

(`ftemp`, `ftheta`, `decact`, `ctrans` remain `solute_step` locals — no declaration changes.)

- [ ] **Step 6: Byte-identical regression gate (full)**

Run: `pixi run -e test check-full`
Expected: PASS — all six cases; solute-active `salinitystress` byte-identical; the two `known_divergence` cases (soilhysteresis, winter) report xfail as before. This is the phase-end gate for the pilot.

- [ ] **Step 7: Commit**

```bash
git add src/solute/solute_kernels.f90 tests/unit/solute/test_solute_kernels.pf src/solute/solute.f90
git commit -m "refactor(solute): extract decomposition factors to kernels

solute_ftemp / solute_ftheta / solute_decomp_ctrans pure functions, unit
tested (temp cap, switch-off, moisture clamp, transformation rate).
solute_step composes them; byte-identical. Completes the solute pure-
kernels pilot.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Done-When

- `solute_kernels_mod` holds 7 tested pure/elemental functions; `tests/unit/solute/test_solute_kernels.pf` is registered and runs (verified by the rising `OK (N tests)` count).
- `solute_step` reads all five coefficients from `state%solute` (no uninitialised locals); the Freundlich and decomposition math is delegated to named kernels.
- `pixi run -e test check-full` passes with `salinitystress` byte-identical and only the two pre-existing `known_divergence` xfails.
- Four commits on `development`, one per task.

## Pilot evaluation (record after Task 4, do not skip)

Per the spec's evaluation criteria, jot a short note (PR description or a memory entry) on: readability of the new `solute_step` loop, diff cost vs. clarity gain, the latent bug the pattern surfaced, and whether to transpose to `soilwater`/`crop`/`heat`.
