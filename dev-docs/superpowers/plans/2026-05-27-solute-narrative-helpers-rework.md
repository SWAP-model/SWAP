# Solute Narrative-Helpers Rework Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.

**Supersedes** the kernels approach in `2026-05-27-solute-pure-kernels-pilot.md` (commits `8c08c14`..`eda838d`), per review feedback: control flags must not be passed into compute functions, and one-/two-line micro-functions add indirection without clarity. Goal: `solute_seed`/`solute_step` read as a narrative of named phases.

**Goal:** Replace the granular pure/elemental kernels with comprehensive, state-operating helper subroutines that each represent one phase of the algorithm; lift the `flTemperature` control branch to where the computation happens; keep the genuine multi-line Freundlich solver as a `pure` function in `solute_mod`.

**Architecture:** `solute_mod` holds `solute_seed`/`solute_step` (narrative orchestrators) + private phase-helper subroutines (each `subroutine …(state, …)`, `intent(inout) :: state`) + one `public pure function solute_cml_from_cmsy`. `src/solute/solute_kernels.f90` is deleted. Per-substep working arrays (`thetav`, `dispr1`, `vpore2`) are passed as explicit args between helpers — not recomputed, not hidden.

**Tech Stack:** Fortran (gfortran), Meson/ninja, pFUnit, Python regression harness, `pixi`.

**Invariant:** BYTE-IDENTICAL on the regression suite at every commit. Preserve exact arithmetic and evaluation order — **move existing blocks verbatim into helpers**; do not retype formulas. `check-fast` after each task, `check-full` at the end. No state-schema change (the `state%solute` coefficient fields from the pilot stay), so no clean rebuild needed.

**What stays from the pilot (do NOT revert):** the five `state%solute` coefficient fields + their allocation in `solute_state_init`; the `swbr==1` quarantine guard; `solute_cml_from_cmsy`'s algorithm.

---

## Target end-state of `solute_mod`

`public :: solute_seed, solute_step, solute_cml_from_cmsy`

Private helpers (all take `type(swap_state_t), intent(inout) :: state` plus listed args):
- `determine_initial_solute_profile(state)` — AFGEN initial-profile interpolation (seed).
- `derive_solute_concentrations(state)` — per-node coefficients (`bdenskf`/`bdenskfcref`/`bdenskfsatporos`/`ddiffwcs`/`decpotfdepth`) + `cmsy` + `samini`/`sampro` (seed).
- `prepare_solute_dispersion(state, thetav, dispr1, vpore2)` — pre-loop dispersion working arrays + `sol%dtsolu` (step).
- `solute_surface_flux(state, cfluxt)` — surface solute flux per substep; returns top flux `cfluxt` (step).
- `update_solute_compartments(state, thetav, dispr1, vpore2, cfluxt, isqdra)` — the per-node loop: convective+dispersive flux, decomposition (with the `time%flTemperature` branch INLINE, read from state — not passed to a leaf fn), root uptake, lateral drainage, conservation, isotherm recovery via `solute_cml_from_cmsy` (step).
- `solute_aquifer_breakthrough(state)` — the two quarantined `swbr==1` blocks (breakthrough `cdrain`/`cseep` + `sqsur`). Internally guarded by `swbr==1`; never reached at runtime because `solute_step` fatals on `swbr==1` first. Carries the FIXME (step).
- `solute_bottom_flux(state)` — bottom-of-profile `sqbot` accumulation per substep (step).
- `solute_balance(state)` — final `sampro`/`sqprec`/`sqirrig` + `solbal` block (step).

`public pure function solute_cml_from_cmsy(cmsy, theta, bdenskf, frexp, cref, cml_guess) result(cml)` — unchanged algorithm, moved from the deleted kernels module.

---

## Task 1: Remove the kernels module; inline micro-kernel math; fold the solver into solute_mod

Reverts the granular-kernel presentation while preserving exact arithmetic. After this task `solute_kernels.f90` is gone and `solute.f90` uses no kernels module.

**Files:**
- Modify `src/solute/solute.f90`
- Delete `src/solute/solute_kernels.f90`
- Modify `meson.build`, `tests/unit/meson.build`, `tests/unit/testSuites.inc`
- Delete `tests/unit/solute/test_solute_kernels.pf`; Create `tests/unit/solute/test_solute_cml.pf`

- [ ] **Step 1: Inline the coefficient arithmetic in `solute_seed`.** In the derived-concentrations loop, replace each `*_coeff(...)` call with the direct arithmetic it wraps (verbatim from how the loop read before the pilot):
  - `sol%bdenskf(i) = soil%bdens(mesh%layer(i))*sol%kf(mesh%layer(i))`
  - `sol%bdenskfsatporos(i) = soil%bdens(mesh%layer(i))*sol%kfsat + sol%poros`
  - `sol%ddiffwcs(i) = sol%ddif / (soil%thetsl(mesh%layer(i))**2)`
  - `sol%decpotfdepth(i) = sol%decpot(mesh%layer(i))*sol%fdepth(mesh%layer(i))`
  - (`sol%bdenskfcref(i) = sol%bdenskf(i)*sol%cref` is already inline — leave it.)

- [ ] **Step 2: Inline the decomposition in `solute_step`** (lift the flag to the orchestrator). Replace the three `solute_ftemp/ftheta/decomp_ctrans` calls with:
```fortran
               ! Solute decomposition.
               if (time%flTemperature) then
                  if (heat%tsoil(i) .lt. 35.0d0) then
                     ftemp = exp(sol%gampar*(heat%tsoil(i) - 20.0d0))
                  else
                     ftemp = exp(sol%gampar*15.0d0)
                  end if
               else
                  ftemp = 0.0d0
               end if
               ftheta = min(1.0d0, (soil%theta(i)/sol%rtheta)**sol%bexp)
               decact = sol%decpotfdepth(i) * ftemp * ftheta
               ctrans = decact*soil%theta(i)*sol%cml(i) +                            &
                        decact*sol%bdenskfcref(i)*((sol%cml(i)/sol%cref)**sol%frexp)
```

- [ ] **Step 3: Move `solute_cml_from_cmsy` into `solute_mod`.** Add it to the module `contains` as a `pure function` (verbatim body from `solute_kernels.f90`), add `solute_cml_from_cmsy` to `solute_mod`'s `public` list. Remove the `use solute_kernels_mod, ...` line from both `solute_seed` and `solute_step`. The call site in `update`/`solute_step` is unchanged.

- [ ] **Step 4: Delete the kernels module + registrations.** `rm src/solute/solute_kernels.f90`; remove its line from `meson.build` (Solute section) and from `tests/unit/meson.build` (pfunit source list); remove `'solute/test_solute_kernels.pf',` from `pf_files` and `ADD_TEST_SUITE(test_solute_kernels_suite)` from `testSuites.inc`. `rm tests/unit/solute/test_solute_kernels.pf`.

- [ ] **Step 5: Add the Freundlich solver test.** Create `tests/unit/solute/test_solute_cml.pf` with the two solver tests from the pilot (linear branch + threshold both-branches, and nonlinear fixed-point residual), but `use solute_mod, only: solute_cml_from_cmsy`. Register: add `'solute/test_solute_cml.pf',` to `pf_files` and `ADD_TEST_SUITE(test_solute_cml_suite)` to `testSuites.inc`. (Keep every `@assert` single-line — the pFUnit preprocessor rejects multi-line asserts.)

- [ ] **Step 6: Verify.** `pixi run -e test check-fast` → all 4 cases "regression ok"; pFUnit `OK` count consistent (the 4 coeff tests + 3 decomp tests removed, the linear+nonlinear solver tests retained). Confirm no `solute_kernels` reference remains: `grep -rn solute_kernels src tests` → empty.

- [ ] **Step 7: Commit.**
```bash
git add -A
git commit -m "refactor(solute): remove kernels module; inline micro-kernel math

Per review: control flags don't belong in compute functions, and one-line
kernels add indirection. Inline the coefficient + decomposition arithmetic
(flTemperature branch now visible in solute_step), fold the Freundlich
solver into solute_mod as a public pure function, delete solute_kernels.f90.
Byte-identical.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 2: Decompose `solute_seed` into named phase helpers

**Files:** Modify `src/solute/solute.f90`; Modify `tests/unit/solute/` (add a state-grain test).

- [ ] **Step 1: Extract `determine_initial_solute_profile(state)`.** Move the AFGEN initial-profile block (the `if (soil%swinco .ne. 3 .and. ...)` interpolation that fills `sol%cml`) verbatim into a private subroutine `determine_initial_solute_profile(state)` with its own `associate`/`use`. Replace the block in `solute_seed` with `call determine_initial_solute_profile(state)`.

- [ ] **Step 2: Extract `derive_solute_concentrations(state)`.** Move the derived-coefficients + `cmsy` + `samini`/`sampro` block verbatim into `derive_solute_concentrations(state)`. Replace in `solute_seed` with `call derive_solute_concentrations(state)`. After this, `solute_seed`'s body is the two calls (plus its `associate` may now be unnecessary — remove unused locals/associate if the compiler flags them).

- [ ] **Step 3: Verify byte-identical.** `pixi run -e test check-fast` → 4/4 ok.

- [ ] **Step 4: Add a state-grain test** for `derive_solute_concentrations`: in a new/existing solute test, build a minimal `swap_state` (allocate `soilwater`/`mesh`/`solute` arrays for a few nodes, set `bdens`/`kf`/`theta`/`dz`/`layer`/`cref`/etc.), call `derive_solute_concentrations(state)`, and assert `sol%bdenskf`, `sol%bdenskfcref`, `sol%cmsy`, and `sol%samini` match hand-computed values. (Requires `derive_solute_concentrations` be accessible — add it to a test-visible interface: simplest is to make it `public` in `solute_mod`, or test via a thin public wrapper. Prefer `public` if that's the lightest path; otherwise note the chosen approach.) Register the suite. `pixi run -e test test-pfunit` → OK count rises.

- [ ] **Step 5: Commit.**
```bash
git add -A
git commit -m "refactor(solute): solute_seed as narrative (profile + derive helpers)

determine_initial_solute_profile + derive_solute_concentrations private
helpers; solute_seed reads as two named phases. State-grain test for the
derive helper. Byte-identical.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Task 3: Decompose `solute_step` into named phase helpers

**Files:** Modify `src/solute/solute.f90`; Modify `tests/unit/solute/` (state-grain test).

Extract these private subroutines, moving the corresponding blocks **verbatim**; the orchestrator keeps the `swbr` guard, cohort resets, the `do while` loop scaffolding, and inline boundary-conc/`isqbot` lines.

- [ ] **Step 1: `prepare_solute_dispersion(state, thetav, dispr1, vpore2)`** — move the pre-loop "Maximum solute time step" block (sets `thetav`/`dispr1`/`vpore2` working arrays and `sol%dtsolu`; `diffus`/`dispr`/`vpore`/`dummy` are locals inside the helper). Declare `thetav(macp)`, `dispr1(macp)`, `vpore2(macp)` in `solute_step` and pass them. Verify compiles.

- [ ] **Step 2: `solute_surface_flux(state, cfluxt)`** — move the "Solute flux at the soil surface" block (updates `sol%csurf`/`sol%cpond`/`sol%isqtop`, sets `cfluxt` out). `intent(out) :: cfluxt`.

- [ ] **Step 3: `update_solute_compartments(state, thetav, dispr1, vpore2, cfluxt, isqdra)`** — move the entire per-node `do i = 1, mesh%numnod ... end do` loop body (convective/dispersive flux, the inline decomposition from Task 1, root uptake, lateral drainage, `isqdra`/`sqdra`/`imsqdra`, conservation, isotherm recovery, `cfluxt = cfluxb`). `cfluxt` is `intent(inout)` (seeded by the surface flux). `isqdra` is `intent(inout)`.

- [ ] **Step 4: `solute_aquifer_breakthrough(state)`** — move BOTH quarantined `if (sol%swbr .eq. 1)` blocks (the `cdrain`/`cseep` breakthrough and the `sqsur` flux). Keep the FIXME explaining it's broken/OOB and unreachable behind the top-level guard. Orchestrator calls `call solute_aquifer_breakthrough(state)` in the loop.

- [ ] **Step 5: `solute_bottom_flux(state)`** — move the "Flux through bottom of soil profile" block (`sqbot`/`imsqbot`).

- [ ] **Step 6: `solute_balance(state)`** — move the post-loop "Solute balance components" block (`sampro` recompute + `sqprec`/`imsqprec`/`sqirrig`/`imsqirrig` + `solbal`). The "Current solute flux at bottom" `isqbot` lines stay inline in `solute_step` (or include in this helper — keep verbatim either way).

  After Steps 1–6, `solute_step` is: guard → resets → boundary conc → `prepare_solute_dispersion` → `do while { advance dt; surface_flux; update_compartments; aquifer_breakthrough; bottom_flux }` → isqbot → `solute_balance`.

- [ ] **Step 7: Verify byte-identical (full).** `pixi run -e test check-full` → 5 pass + 2 known xfail (soilhysteresis/winter); pFUnit OK.

- [ ] **Step 8: Add a state-grain test** for one representative step helper (e.g. `prepare_solute_dispersion` or `update_solute_compartments`): minimal state, call, assert a derived field / `sol%dtsolu`. Make the chosen helper test-visible (public or wrapper). Register; OK count rises.

- [ ] **Step 9: Commit.**
```bash
git add -A
git commit -m "refactor(solute): solute_step as narrative phase helpers

prepare_solute_dispersion / solute_surface_flux / update_solute_compartments
/ solute_aquifer_breakthrough / solute_bottom_flux / solute_balance. The
time-step loop now reads as named stages; per-substep working arrays passed
explicitly. Byte-identical (check-full).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

## Done-When
- `solute_kernels.f90` deleted; no `solute_kernels` reference anywhere.
- `solute_seed` and `solute_step` are narrative orchestrators calling named phase helpers; no control flag is passed into a leaf compute function; no one-line wrapper functions remain (only `solute_cml_from_cmsy`, a genuine algorithm, is a standalone pure function).
- `check-full` green (5 pass + 2 pre-existing xfail); pFUnit green with state-grain helper tests + the Freundlich solver test.
- Three commits on `development`.
