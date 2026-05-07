# [nutrients] N1 — Crop-side nutrient adapter Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Wire the existing `wofost_nutrient_t` typed config to the legacy nutrient-related state via a new per-rotation adapter, and lift the `wofost_nutrient_validate` stub-error so configs with `flcropnut=true` can be loaded and validated. Runtime gate at `tillage.f90:73` stays — N3's job.

**Architecture:** Pre-flight inventory revealed that the nutrient scalars split into two groups: 7 are module-level globals in `variables.f90` (`LRNR`, `LSNR`, `NLUE`, `RNFLV`, `RNFST`, `FRNX`, `NMXLV(30)`); 13 are locals inside `cropgrowth.f90`'s `wofost` subroutine (`NLAI`, `NMAXSO`, `NPART`, `NFIXF`, `NSLA`, `RNFRT`, `TCNT`, `DVSNLT`, `DVSNT`, `RDRNS`, `FNTRT`, `ILNMXL`, `FraHarLosOrm_lv/st/so`). Approach **B** (promote-to-module): move all 13 to `module variables` so a single adapter procedure can populate everything by name. Drop the now-shadowed local declarations. Adapter lives in `cropwofost_init.f90` and fires per-rotation when `cfg%nutrient%flcropnut=true`. Validator switches from "stub-error on flcropnut=true" to "validate the inputs."

**Tech Stack:** Fortran 2008, meson + ninja, gfortran, pFUnit. Reference patterns: `apply_soil_tillage` (ADR 0021) and `apply_irrigation_ssdi` (ADR 0022) for adapter shape; `cropwofost_init_from_config` (`src/crop/cropwofost_init.f90:39`) for the call site.

**Spec:** `docs/superpowers/specs/2026-05-07-nutrients-N1-crop-side-adapter-design.md`.

---

## File Structure

**Modified:**

| Path | Change |
|---|---|
| `src/core/variables.f90` | Add 13 names to the module's nutrient block: `NLAI`, `NMAXSO`, `NPART`, `NFIXF`, `NSLA`, `RNFRT`, `TCNT`, `DVSNLT`, `DVSNT`, `RDRNS`, `FNTRT`, `ILNMXL`, `FraHarLosOrm_lv`, `FraHarLosOrm_st`, `FraHarLosOrm_so`. (15 names total — `FraHarLosOrm_*` is 3.) |
| `src/crop/cropgrowth.f90` | Drop the 15 local declarations from `subroutine wofost(task)` (around lines 1050-1056, 1072) since the names now come via `use variables`. |
| `src/config/cropwofost_config.f90` | Replace stub-error in `wofost_nutrient_validate` with input validation rules. |
| `src/crop/cropwofost_init.f90` | Add `apply_cropwofost_nutrient(nutrient)` public procedure; wire `if (cfg%nutrient%flcropnut) call apply_cropwofost_nutrient(cfg%nutrient)` into `cropwofost_init_from_config` body. |
| `tests/unit/meson.build` | Register new `.pf` file in `pf_files`. |
| `tests/unit/testSuites.inc` | Add `ADD_TEST_SUITE(test_apply_cropwofost_nutrient_suite)`. |
| `docs/adr/index.md` | Append ADR 0025 row. |

**Created:**

| Path | Responsibility |
|---|---|
| `tests/unit/io/toml/test_apply_cropwofost_nutrient.pf` | pFUnit suite: 1 default-init/silent-validate test, 5 validator tests, 1 adapter-populates-globals test. |
| `docs/adr/0025-nutrients-N1-crop-side-adapter.md` | New ADR. |

**Notes from pre-flight inventory:**

- `module variables` already contains `lrnr`, `lsnr`, `nlue`, `nni`, `rnflv`, `rnfst`, `frnx`, `nmxlv(30)`, `anlv`, `anst`, `nmaxlv`, `nmaxst`, `nmaxrt` (lines 567-568). The new names append to that nutrient cluster.
- `cropgrowth.f90:1050-1056` declares `NLAI`, `NMAXSO`, `NPART`, `NFIXF`, `NSLA`, `RNFRT`, `TCNT`, `DVSNLT`, `DVSNT`, `RDRNS`, `FNTRT`, `ILNMXL` as local `real(8)` / `integer`. Comments at `! ,LRNR, LSNR` and `! ,NLUE` etc. confirm those are NOT local — they come from `use variables`. Promoting the 12 names plus `FraHarLosOrm_*` (3 names at line 1072) to module-level eliminates the shadow.
- `wofost` subroutine has `use variables` (wildcard, line 1014). After the local declarations are removed, the names resolve via the wildcard import.
- `cropwofost_init_from_config(cfg, icrop, FraDeceasedLvToSoil)` is the per-rotation init at `cropwofost_init.f90:39`. The adapter call goes near the end of its body (line ~462, just before `end subroutine`).

---

## Task 1: Promote nutrient locals to module variables

**Files:**
- Modify: `src/core/variables.f90` (add 15 names to the nutrient block)
- Modify: `src/crop/cropgrowth.f90` (drop 15 local declarations from `wofost` subroutine)

This is a pure refactor — no behavioural change when `flCropNut=false` (the locals were uninitialized, so was the module init's zero — same result for any path that doesn't touch them). When `flCropNut=true`, the runtime path was stub-erred at `tillage.f90:73` so unreachable in any test. After this task, the names live module-globally and a future caller (the adapter) can populate them.

- [ ] **Step 1: Read the current declarations**

```bash
sed -n '565,570p' src/core/variables.f90
sed -n '1048,1075p' src/crop/cropgrowth.f90
```

You'll see in `variables.f90` the existing nutrient block:
```fortran
      real(8)   nmxlv(30)
      real(8)   nlue,anlv,anst,nmaxlv,nmaxst,nmaxrt
      real(8)   lrnr,lsnr,nni,rnflv,rnfst,frnx
```

And in `cropgrowth.f90:wofost` the locals to be promoted:
```fortran
! --- n-p-k use
      real(8) NLAI   !,LRNR, LSNR
      real(8) NMAXSO, NPART, NFIXF !,NLUE
      real(8) NSLA, RNFRT, TCNT !,RNFLV,RNFST
      real(8) DVSNLT, DVSNT, RDRNS, FNTRT !, FRNX
      !real(8) NMAXLV,NMAXST,NMAXRT,NNI,FSTR
      !real(8) NMXLV(30)
      integer ILNMXL
      ...
      real(8) FraHarLosOrm_lv,FraHarLosOrm_st,FraHarLosOrm_so
```

- [ ] **Step 2: Add the 15 names to `module variables`**

In `src/core/variables.f90`, find the existing nutrient block (around line 565-568) and append:

```fortran
      real(8)   nmxlv(30)
      real(8)   nlue,anlv,anst,nmaxlv,nmaxst,nmaxrt
      real(8)   lrnr,lsnr,nni,rnflv,rnfst,frnx
      ! N-P-K nutrient parameters from cropwofost.nutrient (N1 of [nutrients] umbrella).
      ! Promoted from local-to-wofost-subroutine after spec brainstorming
      ! revealed they could not be reached from a config-load-time adapter
      ! while declared as locals. See ADR 0025.
      real(8)   nlai, nmaxso, npart, nfixf
      real(8)   nsla, rnfrt, tcnt
      real(8)   dvsnlt, dvsnt, rdrns, fntrt
      integer   ilnmxl
      real(8)   fraharlosorm_lv, fraharlosorm_st, fraharlosorm_so
```

(Use the actual case-style of the existing nutrient block — Fortran is case-insensitive but the file convention may be lowercase or mixed.)

- [ ] **Step 3: Drop the local declarations from `cropgrowth.f90:wofost`**

Find the `! --- n-p-k use` block (around line 1049) and remove these lines:

```fortran
! --- n-p-k use
      real(8) NLAI   !,LRNR, LSNR
      real(8) NMAXSO, NPART, NFIXF !,NLUE
      real(8) NSLA, RNFRT, TCNT !,RNFLV,RNFST
      real(8) DVSNLT, DVSNT, RDRNS, FNTRT !, FRNX
      !real(8) NMAXLV,NMAXST,NMAXRT,NNI,FSTR
      !real(8) NMXLV(30)
      integer ILNMXL
```

Keep `real(8) Fstress` and `integer nut` — those are unrelated locals.

Find the `FraHarLosOrm_lv,_st,_so` declaration (around line 1072) and remove just that line:
```fortran
      real(8) FraHarLosOrm_lv,FraHarLosOrm_st,FraHarLosOrm_so
```

Keep the surrounding declarations (`FraDeceasedLvToSoil`, `HarLosOrm_*`, `HarLosNit_*`, etc.).

- [ ] **Step 4: Build and verify**

```
pixi run -e test build-linux
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:|tests,"
pixi run -e test check-full 2>&1 | grep "Results:"
```
Expected: build clean (no shadowing warnings, no unresolved symbols); pFUnit `Ok: 1, Fail: 0` (583 tests, 1 disabled — count unchanged for this task); check-full `5 passed, 0 failed`.

If gfortran complains about a variable being used before declaration, you missed an entry. Re-grep:
```
grep -nE "\\b(NLAI|NMAXSO|NPART|NFIXF|NSLA|RNFRT|TCNT|DVSNLT|DVSNT|RDRNS|FNTRT|ILNMXL|FraHarLosOrm_lv|FraHarLosOrm_st|FraHarLosOrm_so)\\b" src/crop/cropgrowth.f90 | head
```

- [ ] **Step 5: Commit**

```bash
git add src/core/variables.f90 src/crop/cropgrowth.f90
git commit -m "$(cat <<'EOF'
refactor(variables): promote N-P-K locals to module variables

Pre-flight inventory for the [nutrients] N1 sub-arc revealed the
nutrient scalars split into two groups: 7 are already module-level
in variables.f90 (LRNR, LSNR, NLUE, RNFLV, RNFST, FRNX, NMXLV);
the other 15 are locals inside cropgrowth.f90's wofost subroutine
(NLAI, NMAXSO, NPART, NFIXF, NSLA, RNFRT, TCNT, DVSNLT, DVSNT,
RDRNS, FNTRT, ILNMXL, FraHarLosOrm_lv/st/so).

Promote the 15 locals to module variables so a single adapter
procedure can populate everything via `use variables`. The wofost
subroutine continues to read them through its existing wildcard
`use variables` clause; behaviourally identical for the
flCropNut=false path (locals were zero-init in either case;
flCropNut=true has been runtime-stub-erred at tillage.f90:73).

Pure refactor; no test count change. Adapter wiring follows in
the next commit.

Part of [nutrients] N1 — crop-side nutrient adapter (spec
docs/superpowers/specs/2026-05-07-nutrients-N1-crop-side-adapter-design.md).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Validator — drop stub-error + add input checks

**Files:**
- Modify: `src/config/cropwofost_config.f90` (rewrite `wofost_nutrient_validate`)
- Create: `tests/unit/io/toml/test_apply_cropwofost_nutrient.pf` (initial 6-test suite)
- Modify: `tests/unit/meson.build` (register `.pf` file)
- Modify: `tests/unit/testSuites.inc` (register suite)

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/io/toml/test_apply_cropwofost_nutrient.pf` with 6 test subroutines (only the validator-related ones are exercisable until Task 3 lands the adapter; placeholder for Task 3 added below):

```fortran
@test
subroutine test_validate_silent_flcropnut_false()
   use cropwofost_config_mod, only: wofost_nutrient_t
   use error_mod, only: error_collection_t
   use funit
   implicit none
   type(wofost_nutrient_t) :: cfg
   type(error_collection_t) :: errs

   ! Defaults: flcropnut=.false., all numeric defaults zero, nmxlv unallocated
   call cfg%validate(errs)
   call assertEqual(0, errs%count(), 'expected zero errors when flcropnut=false')
end subroutine


@test
subroutine test_validate_passes_clean_flcropnut_true()
   use, intrinsic :: iso_fortran_env, only: real64
   use cropwofost_config_mod, only: wofost_nutrient_t
   use error_mod, only: error_collection_t
   use funit
   implicit none
   type(wofost_nutrient_t) :: cfg
   type(error_collection_t) :: errs

   cfg%flcropnut = .true.
   ! Numeric scalars: defaults are 0.0; that's fine — N1 doesn't range-check them.
   ! What N1 requires:
   !   nmxlv non-empty, size <= 30
   !   harvest fractions in [0, 1]
   allocate(cfg%nmxlv(2))
   cfg%nmxlv = [0.06_real64, 0.04_real64]
   cfg%frahar_los_orm_lv = 0.5_real64
   cfg%frahar_los_orm_st = 0.5_real64
   cfg%frahar_los_orm_so = 0.0_real64

   call cfg%validate(errs)
   call assertEqual(0, errs%count(), 'expected zero errors for clean flcropnut=true config')
end subroutine


@test
subroutine test_validate_rejects_flcropnut_true_empty_nmxlv()
   use cropwofost_config_mod, only: wofost_nutrient_t
   use error_mod, only: error_collection_t
   use funit
   implicit none
   type(wofost_nutrient_t) :: cfg
   type(error_collection_t) :: errs

   cfg%flcropnut = .true.
   ! Leave nmxlv unallocated.
   call cfg%validate(errs)
   @assertTrue(errs%count() > 0, 'expected error for empty nmxlv')
end subroutine


@test
subroutine test_validate_rejects_flcropnut_true_oversized_nmxlv()
   use, intrinsic :: iso_fortran_env, only: real64
   use cropwofost_config_mod, only: wofost_nutrient_t
   use error_mod, only: error_collection_t
   use funit
   implicit none
   type(wofost_nutrient_t) :: cfg
   type(error_collection_t) :: errs

   cfg%flcropnut = .true.
   allocate(cfg%nmxlv(31))
   cfg%nmxlv = 0.05_real64
   call cfg%validate(errs)
   @assertTrue(errs%count() > 0, 'expected error for nmxlv size > 30')
end subroutine


@test
subroutine test_validate_rejects_flcropnut_true_out_of_range_harvest_lv()
   use, intrinsic :: iso_fortran_env, only: real64
   use cropwofost_config_mod, only: wofost_nutrient_t
   use error_mod, only: error_collection_t
   use funit
   implicit none
   type(wofost_nutrient_t) :: cfg
   type(error_collection_t) :: errs

   cfg%flcropnut = .true.
   allocate(cfg%nmxlv(1))
   cfg%nmxlv(1) = 0.05_real64
   cfg%frahar_los_orm_lv = 1.5_real64   ! out of [0, 1]
   call cfg%validate(errs)
   @assertTrue(errs%count() > 0, 'expected error for frahar_los_orm_lv out of range')
end subroutine


! Placeholder; populated in Task 3 once apply_cropwofost_nutrient exists.
! Currently asserts a trivial pass so the suite compiles cleanly.
@test
subroutine test_apply_populates_legacy_globals()
   use funit
   implicit none
   ! Body added in Task 3.
   call assertEqual(1, 1)
end subroutine
```

- [ ] **Step 2: Wire the test file into meson and run to verify failures**

Add to `tests/unit/meson.build`'s `pf_files`:
```meson
        'io/toml/test_apply_cropwofost_nutrient.pf',
```

Add to `tests/unit/testSuites.inc`:
```
ADD_TEST_SUITE(test_apply_cropwofost_nutrient_suite)
```

Run:
```
pixi run -e test build-linux
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:|FAIL|tests,"
```
Expected: build clean. The 4 validator tests with `flcropnut=true` will fail because the current `wofost_nutrient_validate` stub-errors on `flcropnut=true` — the assertions expect zero errors for the clean case and at least one error for the malformed cases, but the stub errors on every `flcropnut=true` config indiscriminately. The `silent_flcropnut_false` test will pass.

- [ ] **Step 3: Replace the stub-error with input validation rules**

In `src/config/cropwofost_config.f90`, find `subroutine wofost_nutrient_validate` (around line 350+ in the file's contains section) and replace its body:

```fortran
   subroutine wofost_nutrient_validate(self, errors)
      use, intrinsic :: iso_fortran_env, only: real64
      use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD, &
                           ERR_VALIDATION_OUT_OF_RANGE
      use validation_mod, only: check_real_range
      class(wofost_nutrient_t), intent(in)    :: self
      type(error_collection_t), intent(inout) :: errors

      if (.not. self%flcropnut) return

      ! nmxlv must be allocated and non-empty when flcropnut=true.
      if (.not. allocated(self%nmxlv) .or. size(self%nmxlv) == 0) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropwofost.nutrient.nmxlv: must be a non-empty array when flcropnut=true', &
            'cropwofost.nutrient.nmxlv')
      else if (size(self%nmxlv) > 30) then
         call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
            'cropwofost.nutrient.nmxlv: at most 30 entries (legacy NMXLV(30) cap)', &
            'cropwofost.nutrient.nmxlv')
      end if

      ! Harvest-loss fractions must be in [0, 1].
      call check_real_range(self%frahar_los_orm_lv, 0.0_real64, 1.0_real64, &
                            'cropwofost.nutrient.frahar_los_orm_lv', errors)
      call check_real_range(self%frahar_los_orm_st, 0.0_real64, 1.0_real64, &
                            'cropwofost.nutrient.frahar_los_orm_st', errors)
      call check_real_range(self%frahar_los_orm_so, 0.0_real64, 1.0_real64, &
                            'cropwofost.nutrient.frahar_los_orm_so', errors)

      ! N1 leaves the 17 numeric scalars (lrnr, lsnr, nlai, ...) unranged.
      ! The legacy reader (rdsdou) didn't enforce ranges either; future
      ! tightening can land when nutrient regression fixtures exist.
   end subroutine wofost_nutrient_validate
```

If `validation_mod` doesn't already export `check_real_range` in scope (it should — it's used by `soil_config_validate`), confirm by:
```
grep -nE "public.*check_real_range" src/validation/validation.f90
```

- [ ] **Step 4: Build and verify the validator tests now pass**

```
pixi run -e test build-linux
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:|tests,"
```
Expected: build clean; pFUnit `Ok: 1, Fail: 0` (test count grows by 6); the 5 validator tests + 1 placeholder all pass.

- [ ] **Step 5: Commit**

```bash
git add src/config/cropwofost_config.f90 tests/unit/io/toml/test_apply_cropwofost_nutrient.pf tests/unit/meson.build tests/unit/testSuites.inc
git commit -m "$(cat <<'EOF'
feat(config): drop wofost_nutrient stub-error; add input validation

Replaces the blanket "flcropnut=true not yet supported" stub-error
in wofost_nutrient_validate with input validation rules:

- nmxlv must be allocated and non-empty when flcropnut=true.
- size(nmxlv) <= 30 (legacy NMXLV(30) cap).
- frahar_los_orm_lv/st/so each in [0, 1].

The 17 numeric scalars (lrnr, lsnr, nlai, nlue, nmaxso, npart, nfixf,
nsla, rnflv, rnfrt, rnfst, tcnt, dvsnlt, dvsnt, rdrns, fntrt, frnx)
stay unranged in N1 — the legacy reader (rdsdou) didn't enforce
ranges either. Future tightening can land when nutrient regression
fixtures exist (post-N3).

Validator stays silent when flcropnut=false (the default). All five
existing regression cases are unaffected — none enable flcropnut.

pFUnit suite test_apply_cropwofost_nutrient seeded with 5 validator
tests + 1 adapter placeholder (filled in Task 3).

This sub-arc still leaves the runtime gate at tillage.f90:73 in
place; lifting that is N3's job after N2 wires the soil-side
initial state.

Part of [nutrients] N1 — crop-side nutrient adapter (spec
docs/superpowers/specs/2026-05-07-nutrients-N1-crop-side-adapter-design.md).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Add `apply_cropwofost_nutrient` adapter + unit test

**Files:**
- Modify: `src/crop/cropwofost_init.f90` (add public `apply_cropwofost_nutrient`)
- Modify: `tests/unit/io/toml/test_apply_cropwofost_nutrient.pf` (replace placeholder with real test)

- [ ] **Step 1: Replace the placeholder test with a real one**

In `tests/unit/io/toml/test_apply_cropwofost_nutrient.pf`, replace `test_apply_populates_legacy_globals` with:

```fortran
@test
subroutine test_apply_populates_legacy_globals()
   use, intrinsic :: iso_fortran_env, only: real64
   use cropwofost_config_mod, only: wofost_nutrient_t
   use cropwofost_init_mod,   only: apply_cropwofost_nutrient
   use variables, only: lrnr, lsnr, nlue, rnflv, rnfst, frnx, nmxlv,           &
                        nlai, nmaxso, npart, nfixf, nsla, rnfrt, tcnt,          &
                        dvsnlt, dvsnt, rdrns, fntrt, ilnmxl,                    &
                        fraharlosorm_lv, fraharlosorm_st, fraharlosorm_so
   use funit
   implicit none
   type(wofost_nutrient_t) :: cfg

   cfg%flcropnut = .true.
   cfg%lrnr   = 0.50_real64
   cfg%lsnr   = 0.50_real64
   cfg%nlai   = 1.00_real64
   cfg%nlue   = 1.10_real64
   cfg%nmaxso = 0.0176_real64
   cfg%npart  = 1.00_real64
   cfg%nfixf  = 0.00_real64
   cfg%nsla   = 0.50_real64
   cfg%rnflv  = 0.0040_real64
   cfg%rnfrt  = 0.0048_real64
   cfg%rnfst  = 0.0015_real64
   cfg%tcnt   = 10.00_real64
   cfg%dvsnlt = 1.00_real64
   cfg%dvsnt  = 0.80_real64
   cfg%rdrns  = 0.05_real64
   cfg%fntrt  = 0.15_real64
   cfg%frnx   = 0.50_real64
   allocate(cfg%nmxlv(3))
   cfg%nmxlv = [0.06_real64, 0.04_real64, 0.02_real64]
   cfg%frahar_los_orm_lv = 0.50_real64
   cfg%frahar_los_orm_st = 0.50_real64
   cfg%frahar_los_orm_so = 0.00_real64

   call apply_cropwofost_nutrient(cfg)

   ! Module-level globals populated from cfg
   @assertEqual(0.50_real64,   lrnr,   tolerance=1.0e-12_real64)
   @assertEqual(0.50_real64,   lsnr,   tolerance=1.0e-12_real64)
   @assertEqual(1.10_real64,   nlue,   tolerance=1.0e-12_real64)
   @assertEqual(0.0040_real64, rnflv,  tolerance=1.0e-12_real64)
   @assertEqual(0.0015_real64, rnfst,  tolerance=1.0e-12_real64)
   @assertEqual(0.50_real64,   frnx,   tolerance=1.0e-12_real64)

   @assertEqual(1.00_real64,   nlai,   tolerance=1.0e-12_real64)
   @assertEqual(0.0176_real64, nmaxso, tolerance=1.0e-12_real64)
   @assertEqual(1.00_real64,   npart,  tolerance=1.0e-12_real64)
   @assertEqual(0.00_real64,   nfixf,  tolerance=1.0e-12_real64)
   @assertEqual(0.50_real64,   nsla,   tolerance=1.0e-12_real64)
   @assertEqual(0.0048_real64, rnfrt,  tolerance=1.0e-12_real64)
   @assertEqual(10.00_real64,  tcnt,   tolerance=1.0e-12_real64)
   @assertEqual(1.00_real64,   dvsnlt, tolerance=1.0e-12_real64)
   @assertEqual(0.80_real64,   dvsnt,  tolerance=1.0e-12_real64)
   @assertEqual(0.05_real64,   rdrns,  tolerance=1.0e-12_real64)
   @assertEqual(0.15_real64,   fntrt,  tolerance=1.0e-12_real64)

   ! NMXLV(1:ilnmxl) populated; ILNMXL = size(cfg%nmxlv)
   call assertEqual(3, ilnmxl, 'ilnmxl matches size(cfg%nmxlv)')
   @assertEqual(0.06_real64, nmxlv(1), tolerance=1.0e-12_real64)
   @assertEqual(0.04_real64, nmxlv(2), tolerance=1.0e-12_real64)
   @assertEqual(0.02_real64, nmxlv(3), tolerance=1.0e-12_real64)

   ! Harvest fractions
   @assertEqual(0.50_real64, fraharlosorm_lv, tolerance=1.0e-12_real64)
   @assertEqual(0.50_real64, fraharlosorm_st, tolerance=1.0e-12_real64)
   @assertEqual(0.00_real64, fraharlosorm_so, tolerance=1.0e-12_real64)
end subroutine
```

- [ ] **Step 2: Run to verify failure**

```
pixi run -e test build-linux
```
Expected: build fails — `apply_cropwofost_nutrient` is not exported by `cropwofost_init_mod`.

- [ ] **Step 3: Add the adapter to `cropwofost_init.f90`**

In `src/crop/cropwofost_init.f90`, add to the module's `public ::` section (line 35):
```fortran
   public :: cropwofost_init_from_config
   public :: apply_cropwofost_nutrient
```

Add the new procedure inside the `contains` block, after `cropwofost_init_from_config`:

```fortran
   !> Apply the per-rotation [wofost.nutrient] config to the legacy
   !! `variables` globals. Called from cropwofost_init_from_config (or
   !! directly from an alternative entry point) when
   !! cfg%nutrient%flcropnut = .true..
   !!
   !! Replaces the deleted rdinit/rdsdou block in cropgrowth.f90's wofost
   !! subroutine (legacy readers physical deletion arc, SS-C step 2).
   !!
   !! See ADR 0025 ([nutrients] N1).
   subroutine apply_cropwofost_nutrient(cfg)
      use cropwofost_config_mod, only: wofost_nutrient_t
      use variables, only: lrnr, lsnr, nlue, rnflv, rnfst, frnx, nmxlv,            &
                           nlai, nmaxso, npart, nfixf, nsla, rnfrt, tcnt,          &
                           dvsnlt, dvsnt, rdrns, fntrt, ilnmxl,                     &
                           fraharlosorm_lv, fraharlosorm_st, fraharlosorm_so
      type(wofost_nutrient_t), intent(in) :: cfg

      integer :: n

      ! Module-level scalars (already exist in module variables)
      lrnr   = cfg%lrnr
      lsnr   = cfg%lsnr
      nlue   = cfg%nlue
      rnflv  = cfg%rnflv
      rnfst  = cfg%rnfst
      frnx   = cfg%frnx

      ! Newly-promoted module variables (Task 1)
      nlai   = cfg%nlai
      nmaxso = cfg%nmaxso
      npart  = cfg%npart
      nfixf  = cfg%nfixf
      nsla   = cfg%nsla
      rnfrt  = cfg%rnfrt
      tcnt   = cfg%tcnt
      dvsnlt = cfg%dvsnlt
      dvsnt  = cfg%dvsnt
      rdrns  = cfg%rdrns
      fntrt  = cfg%fntrt

      ! NMXLV array — copy entries; ILNMXL records the active length
      n = 0
      if (allocated(cfg%nmxlv)) n = size(cfg%nmxlv)
      ilnmxl = n
      nmxlv  = 0.0_real64
      if (n > 0) nmxlv(1:n) = cfg%nmxlv(1:n)

      ! Harvest fractions
      fraharlosorm_lv = cfg%frahar_los_orm_lv
      fraharlosorm_st = cfg%frahar_los_orm_st
      fraharlosorm_so = cfg%frahar_los_orm_so
   end subroutine apply_cropwofost_nutrient
```

- [ ] **Step 4: Build and verify the test passes**

```
pixi run -e test build-linux
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:|tests,"
pixi run -e test check-full 2>&1 | grep "Results:"
```
Expected: build clean; pFUnit `Ok: 1, Fail: 0` (count unchanged from Task 2 — the placeholder counted as one test, the real test counts as one); check-full `5 passed, 0 failed`.

- [ ] **Step 5: Commit**

```bash
git add src/crop/cropwofost_init.f90 tests/unit/io/toml/test_apply_cropwofost_nutrient.pf
git commit -m "$(cat <<'EOF'
feat(crop): apply_cropwofost_nutrient adapter — typed config -> legacy globals

Adds the apply_cropwofost_nutrient public procedure to
cropwofost_init_mod. Copies all 17 nutrient scalars + NMXLV array
+ ILNMXL active-length + 3 harvest-loss fractions from
wofost_nutrient_t to the legacy module-variables globals.

The adapter is decoupled from cropwofost_init_from_config so it can
be called either as part of per-rotation init (next commit) or
directly from a future test/diagnostic entry point.

Replaces the deleted rdinit/rdsdou block in cropgrowth.f90's wofost
subroutine (legacy readers physical deletion arc, SS-C step 2,
commit dddbada).

pFUnit suite test_apply_cropwofost_nutrient now exercises the
adapter end-to-end with a fully-populated cfg and asserts every
global lands at the expected value.

Verified: build clean; pFUnit Ok: 1, Fail: 0; check-full 5/5
(adapter not yet wired into the runtime path; that's the next
commit).

Part of [nutrients] N1 — crop-side nutrient adapter (spec
docs/superpowers/specs/2026-05-07-nutrients-N1-crop-side-adapter-design.md).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Wire `apply_cropwofost_nutrient` into `cropwofost_init_from_config`

**Files:**
- Modify: `src/crop/cropwofost_init.f90` (add the call site at the end of `cropwofost_init_from_config`)

- [ ] **Step 1: Find the end of `cropwofost_init_from_config`**

```bash
sed -n '450,465p' src/crop/cropwofost_init.f90
```

You'll see the tail (around line 462):
```fortran
      dvs     = 0.0d0
      tsum    = 0.0d0
      daycrop = 0
      nofd    = 0

   end subroutine cropwofost_init_from_config
```

- [ ] **Step 2: Insert the adapter call before `end subroutine`**

Replace the snippet above with:
```fortran
      dvs     = 0.0d0
      tsum    = 0.0d0
      daycrop = 0
      nofd    = 0

      ! [nutrients] N1: apply nutrient block to legacy globals when
      ! flcropnut=true on this rotation. Validator and the runtime
      ! gate at tillage.f90:73 still control whether nutrient physics
      ! actually runs; this just makes the typed-config inputs
      ! available if/when the gate is lifted (N3).
      if (cfg%nutrient%flcropnut) call apply_cropwofost_nutrient(cfg%nutrient)

   end subroutine cropwofost_init_from_config
```

- [ ] **Step 3: Build and verify**

```
pixi run -e test build-linux
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:|tests,"
pixi run -e test check-full 2>&1 | grep "Results:"
```
Expected: build clean; pFUnit `Ok: 1, Fail: 0`; check-full `5 passed, 0 failed`. No regression case enables `flcropnut`, so the new call is unreachable in production runtime; the test count and the regression CSV outputs are unchanged.

- [ ] **Step 4: Commit**

```bash
git add src/crop/cropwofost_init.f90
git commit -m "$(cat <<'EOF'
feat(crop): wire apply_cropwofost_nutrient into cropwofost_init_from_config

The per-rotation init now calls the adapter when the rotation's
cropwofost.nutrient.flcropnut flag is true. Replaces the legacy
rdinit/rdsdou block that lived inside cropgrowth.f90's wofost
subroutine and was deleted in SS-C step 2 (commit dddbada).

The runtime gate at tillage.f90:73 is unchanged. With no current
regression case enabling flcropnut, the new call is unreachable
in production runtime today; check-full output is unchanged.
N3 of the [nutrients] umbrella will lift that gate and add the
first regression case.

Verified: build clean; pFUnit Ok: 1, Fail: 0; check-full 5/5.

Part of [nutrients] N1 — crop-side nutrient adapter (spec
docs/superpowers/specs/2026-05-07-nutrients-N1-crop-side-adapter-design.md).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: ADR 0025 + index update

**Files:**
- Create: `docs/adr/0025-nutrients-N1-crop-side-adapter.md`
- Modify: `docs/adr/index.md` (append ADR 0025 row)

- [ ] **Step 1: Write `docs/adr/0025-nutrients-N1-crop-side-adapter.md`**

```markdown
---
title: "ADR 0025 — [nutrients] N1: crop-side nutrient adapter"
date: 2026-05-07
status: accepted
---

# ADR 0025: [nutrients] N1 — crop-side nutrient adapter

## Context

The legacy SWAP nutrient subsystem reads ~25 crop-side
parameters (`LRNR`, `LSNR`, `NLAI`, `NMXLV(30)`,
`FraHarLosOrm_*`, …) from the per-rotation `<cropfil>.crp` via
TTutil's `rdinit`/`rdsdou` block inside `cropgrowth.f90`'s
`wofost` subroutine. That block was deleted in the legacy-readers
physical-deletion arc (SS-C step 2, commit dddbada) because
`flCropNut=1` is stub-errored upstream at `tillage.f90:73` and
the block was unreachable in any modern flow.

A typed-config schema for these parameters already exists:
`wofost_nutrient_t` in `src/config/cropwofost_config.f90`
(populated by the cropwofost TOML parser, commit 230751c). What
was missing: the validator stub-errored on `flcropnut=true`, and
no adapter wired the typed config to the legacy `variables`
module globals that the WOFOST physics body reads.

## Decision

This ADR captures the first sub-arc (N1) of the `[nutrients]`
umbrella. Three sub-arcs total: N1 wires the crop-side
parameters; N2 ports the soil-side initial state; N3 lifts the
runtime gate at `tillage.f90:73` and adds a regression case.

N1 specifically:

1. **Replace `wofost_nutrient_validate`'s blanket stub-error**
   with input validation rules: `nmxlv` non-empty,
   `size(nmxlv) <= 30`, `frahar_los_orm_*` in `[0, 1]`.
   The 17 numeric scalars stay unranged in N1 (matches legacy
   `rdsdou` which also didn't enforce ranges; future tightening
   can land when nutrient regression fixtures exist).
2. **Promote 15 nutrient names from local-to-`wofost`** to
   module-level in `variables.f90` (`NLAI`, `NMAXSO`, `NPART`,
   `NFIXF`, `NSLA`, `RNFRT`, `TCNT`, `DVSNLT`, `DVSNT`,
   `RDRNS`, `FNTRT`, `ILNMXL`, `FraHarLosOrm_lv/st/so`). The
   other 7 nutrient scalars (`LRNR`, `LSNR`, `NLUE`, `RNFLV`,
   `RNFST`, `FRNX`, `NMXLV(30)`) were already module-level.
3. **Add `apply_cropwofost_nutrient(cfg)`** in
   `cropwofost_init.f90` that copies all 22 module-level
   scalars + the `NMXLV` array + `ILNMXL` from the typed config.
4. **Call the adapter** from `cropwofost_init_from_config` per
   rotation when `cfg%nutrient%flcropnut = .true.`.

## What N1 does NOT do

- Lift the runtime stub-error at `tillage.f90:73`. With no
  soil-side initial state wired (N2's job), allowing
  `flCropNut=1` to reach `DoTillage` would let
  `SoilManagement(2..7)` read uninitialised soil-pool state.
  N3 lifts the gate after N2 lands.
- Add a regression case with `flcropnut=true`. Cannot run
  end-to-end until the runtime gate is lifted (N3).
- Range-check the 17 numeric nutrient scalars. The legacy
  reader didn't either; future tightening when fixtures exist.

## Consequences

- A TOML config with `crop.rotation.<n>.cropwofost.nutrient.flcropnut = true`
  loads and validates without error (no stub).
- The `wofost` subroutine reads nutrient parameters from module
  variables that the per-rotation init has populated from the
  typed config. Behaviourally identical to the legacy
  rdinit/rdsdou flow for any populated config.
- All five existing regression cases keep `flcropnut=false`;
  check-full output is byte-identical to pre-arc baseline.
- Promoting the locals to module variables creates an
  architectural consistency: all nutrient parameters live in the
  same place, accessible to anyone with `use variables`. Future
  ADR (per the architectural direction in ADR 0024) may pass
  these as explicit arguments to `wofost`, retiring the module-
  global pattern entirely.

## Acceptance

- `grep -n "not yet supported in the TOML pipeline.*nutrient" src/config/cropwofost_config.f90` → no match.
- `grep -nE "^\s*real\(8\)\s+(NLAI|NMAXSO|NPART|NFIXF|NSLA|RNFRT|TCNT|DVSNLT|DVSNT|RDRNS|FNTRT|FraHarLosOrm)" src/crop/cropgrowth.f90` → no matches (locals were dropped).
- `grep -inE "\bnlai\b|\bnmaxso\b" src/core/variables.f90` → finds the new module-level declarations.
- `grep -n "flCropNut.*not.*allowed" src/crop/tillage.f90` → one match (runtime gate untouched).
- `pixi run -e test test-pfunit` → `Ok: 1, Fail: 0`; suite count grows by 6.
- `pixi run -e test check-full` → `5 passed, 0 failed` (byte-identical CSV outputs).

## Related

- ADR 0008 — error collection over fatalerr.
- ADR 0016 — per-rotation crop config cache (the cache from
  which `apply_cropwofost_nutrient` reads).
- ADR 0021 — tillage TOML port (mirror shape).
- ADR 0022 — SSDI TOML port (mirror shape).
- ADR 0024 — `dtutil.f90` compatibility shim — same architectural
  thread: physics receives parsed inputs.
- Future: ADR 0026 (N2: soil-side `[nutrients]` TOML block).
- Future: ADR 0027 (N3: lift runtime gate + regression case).
```

- [ ] **Step 2: Append the index row**

In `docs/adr/index.md`, after the ADR 0024 row, add:

```markdown
- [ADR 0025 — [nutrients] N1: crop-side nutrient adapter](0025-nutrients-N1-crop-side-adapter.html) — Wires wofost_nutrient_t typed config to legacy variables globals via apply_cropwofost_nutrient adapter; validator stub-error replaced with input checks. Runtime gate at tillage.f90:73 stays in place pending N2/N3.
```

- [ ] **Step 3: Final verification**

```
pixi run -e test build-linux
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:|tests,"
pixi run -e test check-full 2>&1 | grep "Results:"
echo "=== final greps ==="
grep -n "not yet supported in the TOML pipeline.*nutrient" src/config/cropwofost_config.f90 || echo "OK: stub-error removed"
grep -nE "^\s*real\(8\)\s+(NLAI|NMAXSO|NPART|NFIXF|NSLA|RNFRT|TCNT|DVSNLT|DVSNT|RDRNS|FNTRT|FraHarLosOrm)" src/crop/cropgrowth.f90 || echo "OK: locals dropped"
grep -inE "\\bnlai\\b" src/core/variables.f90 | head -3
grep -n "flCropNut.*not.*allowed" src/crop/tillage.f90
```
Expected: build clean; pFUnit Ok: 1, Fail: 0; check-full 5/5; the first two greps return "OK: ..."; `nlai` is found in variables.f90 (~ line 568+); the `flCropNut.*not.*allowed` grep still returns one match (runtime gate intact).

- [ ] **Step 4: Commit**

```bash
git add docs/adr/0025-nutrients-N1-crop-side-adapter.md docs/adr/index.md
git commit -m "$(cat <<'EOF'
docs: ADR 0025 — [nutrients] N1 crop-side nutrient adapter

Captures the first sub-arc of the [nutrients] umbrella. Wires
wofost_nutrient_t typed config to legacy variables globals via
apply_cropwofost_nutrient; replaces the wofost_nutrient_validate
blanket stub-error with input validation rules (nmxlv non-empty
and <= 30 entries; harvest fractions in [0, 1]).

Runtime gate at tillage.f90:73 is intentionally NOT lifted in
N1 — the soil-side initial state isn't wired yet (N2's job).
N3 will lift the gate and add the first regression case.

- docs/adr/0025-nutrients-N1-crop-side-adapter.md: new ADR
  with context, decision, what-N1-doesn't-do, consequences,
  acceptance criteria, and forward references to N2 (ADR 0026
  candidate) and N3 (ADR 0027 candidate).
- docs/adr/index.md: ADR 0025 row.

Closes the [nutrients] N1 spec
(docs/superpowers/specs/2026-05-07-nutrients-N1-crop-side-adapter-design.md).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Self-Review

**Spec coverage:**

| Spec section | Task |
|---|---|
| §Validator changes (drop stub; add rules) | Task 2 |
| §Adapter (`apply_cropwofost_nutrient`) | Task 3 |
| §Wire-up into `cropwofost_init_from_config` | Task 4 |
| §`tillage.f90:73` UNCHANGED | All tasks (verified by acceptance grep in Task 5) |
| §Tests (6 pFUnit suites) | Tasks 2, 3 (5 validator tests in T2 + 1 adapter test in T3) |
| §ADR 0025 + index update | Task 5 |

The spec implies the adapter populates legacy globals, which it does. What the spec did NOT explicitly call out — and what the implementation plan adds — is the **promote-to-module-variables refactor (Task 1)**. This is necessary because pre-flight inventory revealed 15 of the 22 names were locals to the `wofost` subroutine, not module-level. Without Task 1, the adapter can't reach them. The adjustment is documented in ADR 0025's "Decision" section and is consistent with the spec's underlying intent (typed-config flows to legacy globals).

**Acceptance criteria** (from spec):

- `pixi run -e test build-linux` clean — verified each task.
- `pixi run -e test test-pfunit` exits 0; suite count grows by 6 — verified Tasks 2-5.
- `pixi run -e test check-full` exits 0 with `5 passed, 0 failed` — verified Tasks 1, 3, 4, 5.
- Stub-error string gone — verified Task 5 acceptance grep.
- Hand-authored fixture exercises every validator rule — Task 2's 5 validator tests.
- ADR 0025 committed — Task 5.
- `tillage.f90:73` UNCHANGED — verified Task 5 acceptance grep (one match remaining).

**Type / signature consistency:**

- `apply_cropwofost_nutrient(cfg : wofost_nutrient_t)` — used identically in Task 3 (definition) and Task 4 (call site).
- `wofost_nutrient_t` field names (`flcropnut`, `lrnr`, `lsnr`, `nlai`, …, `nmxlv`, `frahar_los_orm_lv/st/so`) match what's in `cropwofost_config.f90`.
- Module variable names in `variables.f90` (after Task 1) — `nlai`, `nmaxso`, …, `ilnmxl`, `fraharlosorm_lv/st/so` — match what the adapter uses in Task 3.
- Note: the `wofost_nutrient_t` type uses `frahar_los_orm_lv` (snake_case with underscore) while the module variable uses `fraharlosorm_lv` (no underscore separator). This mirrors the legacy global naming convention. The adapter in Task 3 maps cleanly between the two.

**Open questions punted to implementation:**

- Exact placement of new module-variable declarations within the existing nutrient block in `variables.f90`. Task 1 Step 2 shows a sample placement; the implementer can adjust to match the file's actual nutrient-block convention if it differs.
- Whether `nmxlv = 0.0_real64` zero-fill in Task 3 is necessary (the `wofost_nutrient_t%nmxlv` allocatable starts at zero size; if cfg has fewer than 30 entries, the wofost body's `afgen(NMXLV, 30, DVS)` walks all 30 — zeros are the correct sentinel for inactive slots). This is the safer choice; can be revisited if profiling shows the zero-fill is hot.

These are pickable at execution time without re-planning.

---

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-05-07-nutrients-N1-crop-side-adapter.md`.
