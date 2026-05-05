# `.crp` port — Phase 2 (cropwofost via case 5 salinitystress) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Port case 5 (salinitystress)'s `potatod.crp` (type 2, WOFOST detailed-crop) to TOML so the new executable's runtime path opens no `.crp` ASCII file for type-2 rotations. Case 5 regression remains 5/5 green throughout.

**Architecture:** Mirrors the Phase 1 (cropfixed) port shape. The per-rotation crop config cache (`crop_config_t.rotation_wofost(:)`), the dispatch in `read_crop_toml.f90` (type=2 branch), and the `crop_config_global` module-level pointer are all already in place — Phase 2 introduces only the runtime init module (`cropwofost_init.f90`) and the dispatch wiring in `wofost(task=1)`. Schema and parser are already at Phase 4c-b level; the audit (Task 1) pins down the residual gap (soybean, bulb, swrdc, nutrient fields) and the `populated` sentinel is the only schema change that gates the runtime path.

**Tech Stack:** Fortran 2008 (gfortran 13), pFUnit unit tests, Meson + Pixi, tomlf TOML parser.

**Spec:** `docs/superpowers/specs/2026-05-02-crp-port-phase2-cropwofost-design.md`

**ADR refs:** `docs/adr/0015-strangler-narrow-scope-stub-errors.md`, `docs/adr/0016-per-rotation-crop-config-cache.md`

**Phase dependency:** Phase 1 (cropfixed) must be implemented before Phase 2. Phase 2 reuses `crop_config_global_mod` (introduced in Phase 1) without modification. If Phase 1 is not yet implemented, do NOT proceed — the `associated(crop_config_global)` guard in Phase 2's dispatch wiring depends on Phase 1's `config_to_variables` setting the pointer.

**Naming reconciliation (spec → plan):**

| Spec name | Actual code name |
|---|---|
| `rotation_cropwofost(:)` (spec) | `rotation_wofost(:)` (actual, on `crop_config_t`) |
| `populated :: logical` per-config | Added in this phase to `cropwofost_config_t` |
| `rotation_loaded(:)` | Already on `crop_config_t`; set by loader when file is read |
| `crop_config_global` | Created by Phase 1 in `src/crop/crop_config_global.f90` |

---

## File map

**Source — modify:**
- `src/config/cropwofost_config.f90` — add missing fields (soybean, bulb, swrdc, nutrient sub-type); add `populated` sentinel; add stub-error validators
- `src/io/toml/read_cropwofost_toml.f90` — parse new sections; set `config%populated = .true.` at end
- `src/crop/cropgrowth.f90` — at line ~984, dispatch on `populated` sentinel between `cropwofost_init_from_config` and legacy `readwofost`
- `meson.build` — register `src/crop/cropwofost_init.f90`
- `tests/unit/meson.build` — register new test files under `tests/unit/crop/`
- `tests/unit/testSuites.inc` — register new test suite `test_cropwofost_init`

**Source — create:**
- `src/crop/cropwofost_init.f90` — runtime init from typed config (replaces `readwofost`'s runtime side-effects on the TOML path)

**Tests — create / modify:**
- `tests/unit/config/test_cropwofost_config.pf` — extend with stub-error tests + `populated` sentinel test + case-5 supported-values pass test
- `tests/unit/io/toml/test_read_cropwofost_toml.pf` — extend: round-trip full `potatod.crp.toml` fixture; assert `populated=.true.`
- `tests/unit/io/toml/test_load_swap_config.pf` — extend: assert case 5 loads all 4 rotation_wofost slots with `populated=.true.` and `rotation_loaded(1..4)=.true.`
- `tests/unit/io/toml/test_salinitystress_parity.pf` — extend: add `test_salinitystress_wofost_init_parity` asserting `cropwofost_init_from_config` produces same globals as `readwofost`
- `tests/unit/crop/test_cropwofost_init.pf` — NEW: verifies init writes correct globals + `cumdens` hand-computed value

**Test fixture — modify (in submodule):**
- `tests/swap-cases/toml/5.salinitystress/potatod.crp.toml` — add `[irrigation_schedule]` section (case 5 `schedule=0`)
- `tests/swap-cases/toml/5.salinitystress/potatod.crp` — DELETE in Task 8 (submodule pair commit)

---

## Submodule discipline

`tests/swap-cases/` is a git submodule. Inner-commit + outer-bump pair is non-negotiable. Never one without the other. Affected tasks: Task 4 (extend potatod.crp.toml) and Task 8 (delete potatod.crp).

When committing in the submodule, use file-scoped git commands (`git commit toml/5.salinitystress/potatod.crp.toml -m "..."`) so pre-existing dirty state from other case directories is not bundled.

---

## Pre-flight commands

```bash
cd /home/zawadzkim/Code/swap
git log --oneline -3
pixi run test-pfunit 2>&1 | tail -5
pixi run regression 2>&1 | tail -10
```

Expected: `Ok: 1, Fail: 0` for unit tests, `5 passed, 0 failed` for regression. If either is red, the baseline is broken — do not proceed.

Also confirm Phase 1 is complete:

```bash
test -f src/crop/crop_config_global.f90 && echo "Phase1 global OK" || echo "MISSING: run Phase 1 first"
test -f src/crop/cropfixed_init.f90     && echo "Phase1 init OK"   || echo "MISSING: run Phase 1 first"
grep -q "cropwofost_init" src/crop/cropgrowth.f90 && echo "already wired" || echo "not yet wired (expected)"
```

---

## Task 1: Audit `readwofost`

Read `src/io/readswap.f90:2520-3435` end-to-end. Classify every non-blank non-comment line into READ / VALIDATE / NORMALIZE / RUNTIME / GUARDED. Output a classification table to `docs/phase-4f-readwofost-audit.md`. This is the foundation for Tasks 2 and 5.

**Files:**
- Create: `docs/phase-4f-readwofost-audit.md`

- [ ] **Step 1: Read `readwofost` lines 2520-3000**

```bash
sed -n '2520,3000p' src/io/readswap.f90
```

- [ ] **Step 2: Read `readwofost` lines 3000-3435**

```bash
sed -n '3000,3435p' src/io/readswap.f90
```

- [ ] **Step 3: Write the audit table**

Create `docs/phase-4f-readwofost-audit.md` with a Markdown table:

```markdown
# readwofost audit — Phase 2 (.crp port cropwofost)

| Lines       | Bucket    | Notes                                                         |
| ----------- | --------- | ------------------------------------------------------------- |
| 2520-2586   | READ      | subroutine header, USE, local decls, rdinit                   |
| 2588-2637   | READ      | swcf, cftb/chtb/cfeictb tables                                |
| 2617-2637   | GUARDED   | swcf=3 (LAI-dependent dual-coeff) block                       |
| 2639-2665   | READ      | swinter + cofab + Gash tables                                 |
| 2643-2665   | GUARDED   | swinter=2 (Gash) + swinter=3 (storage-cap)                    |
| 2667-2678   | READ      | albedo/rsc/rsw when swcf=1/3 standard values vs actual        |
| 2680-2721   | READ      | soybean variant (swsoybean=1 path)                            |
| 2680-2706   | GUARDED   | swsoybean=1 entire block                                      |
| 2712-2729   | READ      | idsl/dlo/dlc/tsumea/tsumam/dtsmtb (non-soybean)               |
| 2723-2742   | READ      | vernalisation (idsl=2 path): verndvs/vernsat/vernbase/vernrtb  |
| 2723-2742   | GUARDED   | idsl=2 vernalization — stub-error (case 5 has idsl=0)         |
| 2731-2742   | GUARDED   | swbulb=1 (bulb crops) entire block                            |
| 2744-2749   | READ      | dvsend, swharv                                                |
| 2751-2761   | READ      | initial: tdwi/laiem/rgrlai; green_area: slatb/spa/ssa/span/tbase |
| 2763-2769   | READ      | assimilation: kdif/kdir/eff/amaxtb/tmpftb/tmnftb              |
| 2771-2795   | READ      | conversion: cvl/cvo/cvr/cvs; respiration: q10/rml/rmo/rmr/rms/rfsetb |
| 2785-2795   | READ      | partitioning: frtb/fltb/fstb/fotb; death: perdl/rdrrtb/rdrstb |
| 2796-2856   | READ+GUARDED | swoxygen, Feddes block (swoxygen=1), Bartholomeus (swoxygen=2 GUARDED) |
| 2858-2867   | READ      | swWrtNonox, aeratecrit                                        |
| 2869-2897   | READ+GUARDED | swdrought, Feddes block (swdrought=1), De Jong (swdrought=2 GUARDED) |
| 2899-2923   | READ+GUARDED | salinity: swsalinity=0/1/2 (all read; =2 requires swdrought=2 GUARDED) |
| 2926-2976   | READ+GUARDED | swcompensate=0/1/2 (=1/2 GUARDED), swstressor, alphacrit, dcritrtz |
| 2978-2996   | READ      | relmf, swpotrelmf                                             |
| 2990-3031   | READ      | swrdc; rdctb; swrd=1/2/3 with rdtb/rdi+rri+rdc+swdmi2rd/rlwtb+wrtmax |
| 3033-3041   | READ      | schedule; schedule=1+swdrought=2 → hlim3h/l/4                |
| 3033-3041   | GUARDED   | schedule=1 entire scheduling path                             |
| 3043-3047   | READ      | FraDeceasedLvToSoil                                           |
| 3049-3083   | READ+GUARDED | swco2=0/1; if swco2=1 read CO2AMAXTB/EFFTB/TRATB + atmofil file (GUARDED) |
| 3085-3089   | RUNTIME   | close file; irrigation(1) call when schedule=1                |
| 3091-3121   | RUNTIME   | cumdens computation from rdctb (when swdrought=1)             |
| 3123-3227   | RUNTIME   | .END-file restart block (when swinco=3 and t1900-tstart<1e-3 and |t1900-cropstart|>=1e-3) |
| 994-1021    | READ+GUARDED | N-P-K nutrient block (second rdinit; gated on flCropNut; GUARDED) |
```

Include a summary section:

```markdown
## Summary

| Bucket   | Line count (approx) |
|----------|---------------------|
| READ     | ~280                |
| VALIDATE | ~15                 |
| NORMALIZE| ~5                  |
| RUNTIME  | ~140                |
| GUARDED  | ~180                |

## Scope for Phase 2

Active (case 5 exercises): swcf=2, swinter=1, swsoybean=0, idsl=0, swbulb=0,
swoxygen=1, swwrtnonox=1, swdrought=1, swsalinity=1, swcompensate=0, swrd=2,
schedule=0, swco2=0, flCropNut=.false., swinco=3 (but condition evaluates false
for case 5 rotations — see spec risk register).

Guarded (stub-errored): swcf=3, swinter=2/3, swsoybean=1, idsl=2 vernalization,
swbulb=1, swoxygen=2, swdrought=2, swsalinity=2, swcompensate=1/2, schedule=1,
swco2=1, flCropNut=.true..
```

- [ ] **Step 4: Commit**

```bash
git add docs/phase-4f-readwofost-audit.md
git commit -m "docs(audit): readwofost line-by-line classification for cropwofost port"
```

---

## Task 2: Extend `cropwofost_config_t` schema + stub-error validators

Add the missing fields (soybean variant, bulb crops, `swrdc`, nutrient sub-type) and the `populated` sentinel. Add stub-error validators for branches not exercised by case 5.

**Files:**
- Modify: `src/config/cropwofost_config.f90`
- Test: `tests/unit/config/test_cropwofost_config.pf`

- [ ] **Step 1: Write the failing stub-error tests**

Append to `tests/unit/config/test_cropwofost_config.pf`:

```fortran
! Phase 2 cropwofost port — stub-error validators for unsupported branches.
! ADR 0015: schema accepts values 1:1; runtime plumbing not yet ported.

@test
subroutine test_cropwofost_swsoybean_one_rejected()
   use funit
   use cropwofost_config_mod, only: cropwofost_config_t
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   type(cropwofost_config_t) :: c
   type(error_collection_t) :: errors
   integer :: i
   logical :: found
   c%soybean%swsoybean = 1
   call c%validate(errors)
   found = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD .and. &
          index(errors%items(i)%message, 'swsoybean=1') > 0) found = .true.
   end do
   @assertTrue(found)
end subroutine

@test
subroutine test_cropwofost_swbulb_one_rejected()
   use funit
   use cropwofost_config_mod, only: cropwofost_config_t
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   type(cropwofost_config_t) :: c
   type(error_collection_t) :: errors
   integer :: i
   logical :: found
   c%bulb%swbulb = 1
   call c%validate(errors)
   found = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD .and. &
          index(errors%items(i)%message, 'swbulb=1') > 0) found = .true.
   end do
   @assertTrue(found)
end subroutine

@test
subroutine test_cropwofost_flcropnut_true_rejected()
   use funit
   use cropwofost_config_mod, only: cropwofost_config_t
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   type(cropwofost_config_t) :: c
   type(error_collection_t) :: errors
   integer :: i
   logical :: found
   c%nutrient%flcropnut = .true.
   call c%validate(errors)
   found = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD .and. &
          index(errors%items(i)%message, 'flcropnut') > 0) found = .true.
   end do
   @assertTrue(found)
end subroutine

@test
subroutine test_cropwofost_swco2_one_rejected()
   use funit
   use cropwofost_config_mod, only: cropwofost_config_t
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   type(cropwofost_config_t) :: c
   type(error_collection_t) :: errors
   integer :: i
   logical :: found
   c%co2%swco2 = 1
   call c%validate(errors)
   found = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD .and. &
          index(errors%items(i)%message, 'swco2=1') > 0) found = .true.
   end do
   @assertTrue(found)
end subroutine

@test
subroutine test_cropwofost_schedule_one_rejected()
   use funit
   use cropwofost_config_mod, only: cropwofost_config_t
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   type(cropwofost_config_t) :: c
   type(error_collection_t) :: errors
   integer :: i
   logical :: found
   c%schedule%schedule = 1
   call c%validate(errors)
   found = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD .and. &
          index(errors%items(i)%message, 'schedule=1') > 0) found = .true.
   end do
   @assertTrue(found)
end subroutine

@test
subroutine test_cropwofost_swdrought_two_rejected()
   use funit
   use cropwofost_config_mod, only: cropwofost_config_t
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   type(cropwofost_config_t) :: c
   type(error_collection_t) :: errors
   integer :: i
   logical :: found
   c%drought_stress%swdrought = 2
   call c%validate(errors)
   found = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD .and. &
          index(errors%items(i)%message, 'swdrought=2') > 0) found = .true.
   end do
   @assertTrue(found)
end subroutine

@test
subroutine test_cropwofost_swoxygen_two_rejected()
   use funit
   use cropwofost_config_mod, only: cropwofost_config_t
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   type(cropwofost_config_t) :: c
   type(error_collection_t) :: errors
   integer :: i
   logical :: found
   c%oxygen_stress%swoxygen = 2
   call c%validate(errors)
   found = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD .and. &
          index(errors%items(i)%message, 'swoxygen=2') > 0) found = .true.
   end do
   @assertTrue(found)
end subroutine

@test
subroutine test_cropwofost_swinter_two_rejected()
   use funit
   use cropwofost_config_mod, only: cropwofost_config_t
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   type(cropwofost_config_t) :: c
   type(error_collection_t) :: errors
   integer :: i
   logical :: found
   c%interception%swinter = 2
   call c%validate(errors)
   found = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD .and. &
          index(errors%items(i)%message, 'swinter=2') > 0) found = .true.
   end do
   @assertTrue(found)
end subroutine

@test
subroutine test_cropwofost_swcompensate_one_rejected()
   use funit
   use cropwofost_config_mod, only: cropwofost_config_t
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   type(cropwofost_config_t) :: c
   type(error_collection_t) :: errors
   integer :: i
   logical :: found
   c%compensate%swcompensate = 1
   call c%validate(errors)
   found = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD .and. &
          index(errors%items(i)%message, 'swcompensate') > 0) found = .true.
   end do
   @assertTrue(found)
end subroutine

@test
subroutine test_cropwofost_swharv_one_rejected()
   use funit
   use cropwofost_config_mod, only: cropwofost_config_t
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   type(cropwofost_config_t) :: c
   type(error_collection_t) :: errors
   integer :: i
   logical :: found
   c%harvest%swharv = 1
   call c%validate(errors)
   found = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD .and. &
          index(errors%items(i)%message, 'swharv=1') > 0) found = .true.
   end do
   @assertTrue(found)
end subroutine

@test
subroutine test_cropwofost_swsalinity_two_rejected()
   use funit
   use cropwofost_config_mod, only: cropwofost_config_t
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   type(cropwofost_config_t) :: c
   type(error_collection_t) :: errors
   integer :: i
   logical :: found
   ! swsalinity=2 (osmotic head) requires swdrought=2 which is itself stub-errored.
   c%salinity%swsalinity = 2
   call c%validate(errors)
   found = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD .and. &
          index(errors%items(i)%message, 'swsalinity=2') > 0) found = .true.
   end do
   @assertTrue(found)
end subroutine

@test
subroutine test_cropwofost_populated_false_by_default()
   use funit
   use cropwofost_config_mod, only: cropwofost_config_t
   type(cropwofost_config_t) :: c
   @assertFalse(c%populated)
end subroutine

@test
subroutine test_cropwofost_case5_supported_values_no_stub_errors()
   use funit
   use iso_fortran_env, only: real64
   use cropwofost_config_mod, only: cropwofost_config_t
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   type(cropwofost_config_t) :: c
   type(error_collection_t) :: errors
   integer :: i
   logical :: any_stub

   ! Set switches to case-5 active values — none should trigger stub-errors.
   c%soybean%swsoybean        = 0
   c%bulb%swbulb              = 0
   c%nutrient%flcropnut       = .false.
   c%harvest%swharv           = 0
   c%harvest%dvsend           = 3.0_real64
   c%crop_factor%swcf         = 2
   c%crop_factor%albedo       = 0.19_real64
   c%crop_factor%rsc          = 207.0_real64
   c%crop_factor%rsw          = 0.0_real64
   c%phenology%idsl           = 0
   c%phenology%tsumea         = 150.0_real64
   c%phenology%tsumam         = 1550.0_real64
   c%oxygen_stress%swoxygen   = 1
   c%oxygen_stress%swwrtnonox = 1
   c%oxygen_stress%aeratecrit = 0.5_real64
   c%oxygen_stress%hlim1      = -10.0_real64
   c%oxygen_stress%hlim2u     = -25.0_real64
   c%oxygen_stress%hlim2l     = -25.0_real64
   c%drought_stress%swdrought = 1
   c%drought_stress%hlim3h    = -300.0_real64
   c%drought_stress%hlim3l    = -500.0_real64
   c%drought_stress%hlim4     = -10000.0_real64
   c%drought_stress%adcrh     = 0.5_real64
   c%drought_stress%adcrl     = 0.1_real64
   c%salinity%swsalinity      = 1
   c%salinity%saltmax         = 0.732_real64
   c%salinity%saltslope       = 0.0868_real64
   c%compensate%swcompensate  = 0
   c%interception%swinter     = 1
   c%interception%cofab       = 0.25_real64
   c%co2%swco2                = 0
   c%schedule%schedule        = 0
   c%root%swrd                = 2
   c%root%swrdc               = 0
   c%root%rdi                 = 10.0_real64
   c%root%rri                 = 1.2_real64
   c%root%rdc                 = 50.0_real64
   c%root%swdmi2rd            = 1
   c%management%swpotrelmf    = 2
   c%management%relmf         = 0.8_real64

   call c%validate(errors)
   any_stub = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD) any_stub = .true.
   end do
   @assertFalse(any_stub)
end subroutine
```

- [ ] **Step 2: Run tests — verify failure**

```bash
pixi run test-pfunit 2>&1 | grep -E "FAIL|ERROR|stub" | head -20
```

Expected: at least 12 new stub-error tests fail; `populated` default test may fail because the field doesn't exist yet.

- [ ] **Step 3: Add new sub-types and `populated` to `cropwofost_config_t`**

In `src/config/cropwofost_config.f90`:

After the `use error_mod` line, add `ERR_VALIDATION_CROSS_FIELD` to the error imports:

```fortran
   use error_mod, only: error_collection_t, ERR_VALIDATION_OUT_OF_RANGE, &
                        ERR_VALIDATION_CROSS_FIELD
```

Add three new public type names after the existing list:

```fortran
   public :: wofost_soybean_t
   public :: wofost_bulb_t
   public :: wofost_nutrient_t
```

Insert the three new sub-type definitions after the `wofost_management_t` block and before the `cropwofost_config_t` top-level block:

```fortran
   ! ------------------------------------------------------------------
   ! Soybean variant (stub-errored when swsoybean=1)
   ! ------------------------------------------------------------------
   type :: wofost_soybean_t
      integer      :: swsoybean    = 0
      real(real64) :: mg           = 0.0_real64   ! maturity group
      real(real64) :: dvsi         = 0.0_real64
      real(real64) :: dvrmax1      = 0.0_real64
      real(real64) :: dvrmax2      = 0.0_real64
      real(real64) :: tmaxdvr      = 0.0_real64
      real(real64) :: tmindvr      = 0.0_real64
      real(real64) :: toptdvr      = 0.0_real64
      logical      :: flrfphotoveg = .false.
      logical      :: flphenodayl  = .false.
      real(real64) :: popt         = 0.0_real64
      real(real64) :: pcrt         = 0.0_real64
   contains
      procedure :: validate => wofost_soybean_validate
      procedure :: finalize => wofost_soybean_finalize
   end type wofost_soybean_t

   ! ------------------------------------------------------------------
   ! Bulb crops (stub-errored when swbulb=1)
   ! ------------------------------------------------------------------
   type :: wofost_bulb_t
      integer      :: swbulb = 0
      real(real64) :: pld    = 0.0_real64   ! planting density
      real(real64) :: plwti  = 0.0_real64   ! initial weight of planting material
      real(real64) :: remoc  = 0.0_real64   ! remobilisation coefficient
      real(real64), allocatable :: fbltb(:,:)  ! fraction to bulb vs DVS
   contains
      procedure :: validate => wofost_bulb_validate
      procedure :: finalize => wofost_bulb_finalize
   end type wofost_bulb_t

   ! ------------------------------------------------------------------
   ! N-P-K nutrient (stub-errored when flcropnut=.true.)
   ! ------------------------------------------------------------------
   type :: wofost_nutrient_t
      logical      :: flcropnut = .false.
      real(real64) :: lrnr      = 0.0_real64
      real(real64) :: lsnr      = 0.0_real64
      real(real64) :: nlai      = 0.0_real64
      real(real64) :: nlue      = 0.0_real64
      real(real64) :: nmaxso    = 0.0_real64
      real(real64) :: npart     = 0.0_real64
      real(real64) :: nfixf     = 0.0_real64
      real(real64) :: nsla      = 0.0_real64
      real(real64) :: rnflv     = 0.0_real64
      real(real64) :: rnfrt     = 0.0_real64
      real(real64) :: rnfst     = 0.0_real64
      real(real64) :: tcnt      = 0.0_real64
      real(real64) :: dvsnlt    = 0.0_real64
      real(real64) :: dvsnt     = 0.0_real64
      real(real64) :: rdrns     = 0.0_real64
      real(real64) :: fntrt     = 0.0_real64
      real(real64) :: frnx      = 0.0_real64
      real(real64), allocatable :: nmxlv(:,:)  ! max N concentration in leaves vs DVS
   contains
      procedure :: validate => wofost_nutrient_validate
      procedure :: finalize => wofost_nutrient_finalize
   end type wofost_nutrient_t
```

In `wofost_root_t`, add the `swrdc` field after the `swdmi2rd` field:

```fortran
      integer      :: swrdc     = 0
```

In the `cropwofost_config_t` type, add the three new sub-types and the `populated` sentinel:

```fortran
   type :: cropwofost_config_t
      ! ... existing fields unchanged ...
      type(wofost_soybean_t)         :: soybean
      type(wofost_bulb_t)            :: bulb
      type(wofost_nutrient_t)        :: nutrient
      logical                        :: populated = .false.
   contains
      procedure :: validate => cropwofost_config_validate
      procedure :: finalize => cropwofost_config_finalize
   end type cropwofost_config_t
```

- [ ] **Step 4: Add sub-type procedure implementations**

After the existing `wofost_management_finalize` subroutine, add validate + finalize stubs for the three new types, and add the `swrdc` stub-error to `wofost_root_validate`:

```fortran
   subroutine wofost_soybean_validate(self, errors)
      class(wofost_soybean_t),  intent(in)    :: self
      type(error_collection_t), intent(inout) :: errors
      call check_int_enum(self%swsoybean, [0, 1], 'wofost.soybean.swsoybean', errors)
      if (self%swsoybean == 1) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropwofost.soybean.swsoybean=1 (soybean variant) not yet ' // &
            'supported in the TOML pipeline; use the legacy executable.', &
            'cropwofost.soybean')
      end if
   end subroutine wofost_soybean_validate

   subroutine wofost_soybean_finalize(self, errors)
      class(wofost_soybean_t),  intent(inout) :: self
      type(error_collection_t), intent(inout) :: errors
      return
   end subroutine wofost_soybean_finalize

   subroutine wofost_bulb_validate(self, errors)
      class(wofost_bulb_t),     intent(in)    :: self
      type(error_collection_t), intent(inout) :: errors
      call check_int_enum(self%swbulb, [0, 1], 'wofost.bulb.swbulb', errors)
      if (self%swbulb == 1) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropwofost.bulb.swbulb=1 (bulb crops) not yet supported ' // &
            'in the TOML pipeline; use the legacy executable.', &
            'cropwofost.bulb')
      end if
   end subroutine wofost_bulb_validate

   subroutine wofost_bulb_finalize(self, errors)
      class(wofost_bulb_t),     intent(inout) :: self
      type(error_collection_t), intent(inout) :: errors
      return
   end subroutine wofost_bulb_finalize

   subroutine wofost_nutrient_validate(self, errors)
      class(wofost_nutrient_t), intent(in)    :: self
      type(error_collection_t), intent(inout) :: errors
      if (self%flcropnut) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropwofost.nutrient.flcropnut=.true. (N-P-K nutrient model) ' // &
            'not yet supported in the TOML pipeline; use the legacy executable.', &
            'cropwofost.nutrient')
      end if
   end subroutine wofost_nutrient_validate

   subroutine wofost_nutrient_finalize(self, errors)
      class(wofost_nutrient_t), intent(inout) :: self
      type(error_collection_t), intent(inout) :: errors
      return
   end subroutine wofost_nutrient_finalize
```

In `wofost_root_validate`, add after the `swdmi2rd` check:

```fortran
      if (self%swrdc == 1) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropwofost.root.swrdc=1 not yet supported in the TOML ' // &
            'pipeline; use the legacy executable.', 'cropwofost.root')
      end if
```

In `cropwofost_config_validate`, add the top-level stub-error block at the very start (before any delegation):

```fortran
   subroutine cropwofost_config_validate(self, errors)
      class(cropwofost_config_t), intent(in)    :: self
      type(error_collection_t),   intent(inout) :: errors

      ! Phase 2 stub-errors (ADR 0015)
      if (self%drought_stress%swdrought == 2) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropwofost.drought_stress.swdrought=2 (De Jong van Lier) not yet ' // &
            'supported in the TOML pipeline; use the legacy executable.', &
            'cropwofost.drought_stress')
      end if
      if (self%oxygen_stress%swoxygen == 2) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropwofost.oxygen_stress.swoxygen=2 (Bartholomeus) not yet ' // &
            'supported in the TOML pipeline; use the legacy executable.', &
            'cropwofost.oxygen_stress')
      end if
      if (self%interception%swinter == 2) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropwofost.interception.swinter=2 (Gash forest interception) not ' // &
            'yet supported in the TOML pipeline; use the legacy executable.', &
            'cropwofost.interception')
      end if
      if (self%compensate%swcompensate /= 0) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropwofost.compensate.swcompensate /= 0 (Jarvis/Walsum) not yet ' // &
            'supported in the TOML pipeline; use the legacy executable.', &
            'cropwofost.compensate')
      end if
      if (self%harvest%swharv == 1) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropwofost.harvest.swharv=1 (DVS-based harvest timing) not yet ' // &
            'supported in the TOML pipeline; use the legacy executable.', &
            'cropwofost.harvest')
      end if
      if (self%co2%swco2 == 1) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropwofost.co2.swco2=1 (CO2 assimilation correction) not yet ' // &
            'supported in the TOML pipeline; use the legacy executable.', &
            'cropwofost.co2')
      end if
      if (self%schedule%schedule == 1) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropwofost.irrigation_schedule.schedule=1 not yet supported in ' // &
            'the TOML pipeline; use the legacy executable.', &
            'cropwofost.schedule')
      end if
      if (self%salinity%swsalinity == 2) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropwofost.salinity.swsalinity=2 (osmotic head) not yet supported ' // &
            'in the TOML pipeline; use the legacy executable.', &
            'cropwofost.salinity')
      end if

      ! Delegate to sub-types
      call self%preparation%validate(errors)
      call self%sowing%validate(errors)
      call self%germination%validate(errors)
      call self%harvest%validate(errors)
      call self%crop_factor%validate(errors)
      call self%phenology%validate(errors)
      call self%initial%validate(errors)
      call self%green_area%validate(errors)
      call self%assimilation%validate(errors)
      call self%conversion%validate(errors)
      call self%respiration%validate(errors)
      call self%partitioning%validate(errors)
      call self%death%validate(errors)
      call self%root%validate(errors)
      call self%oxygen_stress%validate(errors)
      call self%drought_stress%validate(errors)
      call self%salinity%validate(errors)
      call self%compensate%validate(errors)
      call self%interception%validate(errors)
      call self%co2%validate(errors)
      call self%management%validate(errors)
      call self%soybean%validate(errors)
      call self%bulb%validate(errors)
      call self%nutrient%validate(errors)
      call self%schedule%validate(errors)
   end subroutine cropwofost_config_validate
```

Also update `cropwofost_config_finalize` to call the new sub-type finalizers:

```fortran
      call self%soybean%finalize(errors)
      call self%bulb%finalize(errors)
      call self%nutrient%finalize(errors)
```

- [ ] **Step 5: Run tests — verify pass**

```bash
pixi run test-pfunit 2>&1 | tail -10
```

Expected: `Ok: 1, Fail: 0`. If a message substring doesn't match, adjust it to match the literal text emitted by the validator.

- [ ] **Step 6: Commit**

```bash
git add src/config/cropwofost_config.f90 tests/unit/config/test_cropwofost_config.pf
git commit -m "feat(config): cropwofost schema gap-fill + stub-error validators + populated sentinel

Phase 2 of the .crp port (cropwofost via case 5). Adds soybean, bulb, and
nutrient sub-types for schema 1:1 with readwofost. Adds swrdc field to
wofost_root_t. Adds populated sentinel. Stub-errors swdrought=2,
swoxygen=2, swinter=2, swcompensate/=0, swharv=1, swco2=1, schedule=1,
swsalinity=2, swsoybean=1, swbulb=1, flcropnut=.true. per ADR 0015."
```

---

## Task 3: Extend `read_cropwofost_toml.f90` + set `populated`

Parse the new sub-sections (soybean, bulb, nutrient) and set `config%populated = .true.` at the end of `read_cropwofost_toml`.

**Files:**
- Modify: `src/io/toml/read_cropwofost_toml.f90`
- Test: `tests/unit/io/toml/test_read_cropwofost_toml.pf`

- [ ] **Step 1: Write the failing populated test**

Append to `tests/unit/io/toml/test_read_cropwofost_toml.pf`:

```fortran
@test
subroutine test_read_cropwofost_sets_populated()
   use funit
   use tomlf, only: toml_table, toml_load
   use cropwofost_config_mod, only: cropwofost_config_t
   use read_cropwofost_toml_mod, only: read_cropwofost_toml
   use error_mod, only: error_collection_t
   type(toml_table), allocatable, target :: doc
   type(toml_table), pointer             :: doc_ptr
   type(cropwofost_config_t)             :: c
   type(error_collection_t)              :: errors

   call toml_load(doc, 'tests/unit/io/toml/fixtures/cropwofost_minimal.toml')
   doc_ptr => doc
   @assertFalse(c%populated)
   call read_cropwofost_toml(doc_ptr, c, errors)
   @assertFalse(errors%has_errors())
   @assertTrue(c%populated)
end subroutine

@test
subroutine test_read_cropwofost_potatod_salinity()
   ! Round-trip case 5's potatod.crp.toml — verify the salinity and drought
   ! fields that are the headline feature of case 5.
   use funit
   use iso_fortran_env, only: real64
   use tomlf, only: toml_table, toml_load, toml_error
   use cropwofost_config_mod, only: cropwofost_config_t
   use read_cropwofost_toml_mod, only: read_cropwofost_toml
   use error_mod, only: error_collection_t
   type(toml_table), allocatable, target :: doc
   type(toml_table), pointer             :: doc_ptr
   type(cropwofost_config_t)             :: c
   type(error_collection_t)              :: errors
   type(toml_error), allocatable         :: terr

   call toml_load(doc, &
      'tests/swap-cases/toml/5.salinitystress/potatod.crp.toml', error=terr)
   @assertFalse(allocated(terr))
   doc_ptr => doc
   call read_cropwofost_toml(doc_ptr, c, errors)
   @assertFalse(errors%has_errors())
   @assertTrue(c%populated)

   ! Headline features: salinity=1 (Maas-Hoffman)
   @assertEqual(1,          c%salinity%swsalinity)
   @assertEqual(0.732_real64,  c%salinity%saltmax,   1.0e-12_real64)
   @assertEqual(0.0868_real64, c%salinity%saltslope, 1.0e-12_real64)

   ! Drought=1 (Feddes)
   @assertEqual(1,             c%drought_stress%swdrought)
   @assertEqual(-300.0_real64, c%drought_stress%hlim3h, 1.0e-12_real64)
   @assertEqual(-500.0_real64, c%drought_stress%hlim3l, 1.0e-12_real64)
   @assertEqual(-10000.0_real64, c%drought_stress%hlim4, 1.0e-12_real64)

   ! Oxygen=1 (Feddes)
   @assertEqual(1,            c%oxygen_stress%swoxygen)
   @assertEqual(-10.0_real64, c%oxygen_stress%hlim1,  1.0e-12_real64)
   @assertEqual(-25.0_real64, c%oxygen_stress%hlim2u, 1.0e-12_real64)
   @assertEqual(-25.0_real64, c%oxygen_stress%hlim2l, 1.0e-12_real64)

   ! Root depth=2 (daily increase)
   @assertEqual(2,           c%root%swrd)
   @assertEqual(10.0_real64, c%root%rdi, 1.0e-12_real64)
   @assertEqual(1.2_real64,  c%root%rri, 1.0e-12_real64)
   @assertEqual(50.0_real64, c%root%rdc, 1.0e-12_real64)
   @assertEqual(1,           c%root%swdmi2rd)

   ! Management
   @assertEqual(2,          c%management%swpotrelmf)
   @assertEqual(0.8_real64, c%management%relmf, 1.0e-12_real64)
   @assertEqual(0.3_real64, c%management%fradeceasedlvtosoil, 1.0e-12_real64)
end subroutine
```

- [ ] **Step 2: Run tests — verify failure**

```bash
pixi run test-pfunit 2>&1 | grep -E "populated|potatod_salinity|FAIL" | head -10
```

Expected: both new tests fail (`populated` not set, parser hasn't changed yet).

- [ ] **Step 3: Add parser sections for new sub-types + set `populated`**

In `src/io/toml/read_cropwofost_toml.f90`, add the three new sub-sections inside `read_cropwofost_toml` before the `irrigation_schedule` call, and set `populated` at the very end:

```fortran
      ! ----------------------------------------------------------------
      ! Soybean variant (optional section; absent in most cases)
      ! ----------------------------------------------------------------
      call get_table(doc_root, 'soybean', soy, 'soybean', errors)
      if (associated(soy)) then
         call get_optional_int_with_default(soy, 'swsoybean',    config%soybean%swsoybean,    0,          'soybean.swsoybean',    errors)
         call get_optional_real_with_default(soy, 'mg',          config%soybean%mg,           0.0_real64, 'soybean.mg',           errors)
         call get_optional_real_with_default(soy, 'dvsi',        config%soybean%dvsi,         0.0_real64, 'soybean.dvsi',         errors)
         call get_optional_real_with_default(soy, 'dvrmax1',     config%soybean%dvrmax1,      0.0_real64, 'soybean.dvrmax1',      errors)
         call get_optional_real_with_default(soy, 'dvrmax2',     config%soybean%dvrmax2,      0.0_real64, 'soybean.dvrmax2',      errors)
         call get_optional_real_with_default(soy, 'tmaxdvr',     config%soybean%tmaxdvr,      0.0_real64, 'soybean.tmaxdvr',      errors)
         call get_optional_real_with_default(soy, 'tmindvr',     config%soybean%tmindvr,      0.0_real64, 'soybean.tmindvr',      errors)
         call get_optional_real_with_default(soy, 'toptdvr',     config%soybean%toptdvr,      0.0_real64, 'soybean.toptdvr',      errors)
         call get_optional_real_with_default(soy, 'popt',        config%soybean%popt,         0.0_real64, 'soybean.popt',         errors)
         call get_optional_real_with_default(soy, 'pcrt',        config%soybean%pcrt,         0.0_real64, 'soybean.pcrt',         errors)
      end if

      ! ----------------------------------------------------------------
      ! Bulb crops (optional section; absent in most cases)
      ! ----------------------------------------------------------------
      call get_table(doc_root, 'bulb', bul, 'bulb', errors)
      if (associated(bul)) then
         call get_optional_int_with_default(bul,  'swbulb', config%bulb%swbulb, 0,          'bulb.swbulb', errors)
         call get_optional_real_with_default(bul, 'pld',    config%bulb%pld,    0.0_real64, 'bulb.pld',    errors)
         call get_optional_real_with_default(bul, 'plwti',  config%bulb%plwti,  0.0_real64, 'bulb.plwti',  errors)
         call get_optional_real_with_default(bul, 'remoc',  config%bulb%remoc,  0.0_real64, 'bulb.remoc',  errors)
         call read_table_2d(bul, 'fbltb', config%bulb%fbltb, 2, 'bulb.fbltb', errors)
      end if

      ! ----------------------------------------------------------------
      ! Nutrient model (optional section; absent in most cases)
      ! ----------------------------------------------------------------
      call get_table(doc_root, 'nutrient', nut, 'nutrient', errors)
      if (associated(nut)) then
         call get_optional_real_with_default(nut, 'lrnr',   config%nutrient%lrnr,   0.0_real64, 'nutrient.lrnr',   errors)
         call get_optional_real_with_default(nut, 'lsnr',   config%nutrient%lsnr,   0.0_real64, 'nutrient.lsnr',   errors)
         call get_optional_real_with_default(nut, 'nlai',   config%nutrient%nlai,   0.0_real64, 'nutrient.nlai',   errors)
         call get_optional_real_with_default(nut, 'nlue',   config%nutrient%nlue,   0.0_real64, 'nutrient.nlue',   errors)
         call get_optional_real_with_default(nut, 'nmaxso', config%nutrient%nmaxso, 0.0_real64, 'nutrient.nmaxso', errors)
         call get_optional_real_with_default(nut, 'npart',  config%nutrient%npart,  0.0_real64, 'nutrient.npart',  errors)
         call get_optional_real_with_default(nut, 'nfixf',  config%nutrient%nfixf,  0.0_real64, 'nutrient.nfixf',  errors)
         call get_optional_real_with_default(nut, 'nsla',   config%nutrient%nsla,   0.0_real64, 'nutrient.nsla',   errors)
         call get_optional_real_with_default(nut, 'rnflv',  config%nutrient%rnflv,  0.0_real64, 'nutrient.rnflv',  errors)
         call get_optional_real_with_default(nut, 'rnfrt',  config%nutrient%rnfrt,  0.0_real64, 'nutrient.rnfrt',  errors)
         call get_optional_real_with_default(nut, 'rnfst',  config%nutrient%rnfst,  0.0_real64, 'nutrient.rnfst',  errors)
         call get_optional_real_with_default(nut, 'tcnt',   config%nutrient%tcnt,   0.0_real64, 'nutrient.tcnt',   errors)
         call get_optional_real_with_default(nut, 'dvsnlt', config%nutrient%dvsnlt, 0.0_real64, 'nutrient.dvsnlt', errors)
         call get_optional_real_with_default(nut, 'dvsnt',  config%nutrient%dvsnt,  0.0_real64, 'nutrient.dvsnt',  errors)
         call get_optional_real_with_default(nut, 'rdrns',  config%nutrient%rdrns,  0.0_real64, 'nutrient.rdrns',  errors)
         call get_optional_real_with_default(nut, 'fntrt',  config%nutrient%fntrt,  0.0_real64, 'nutrient.fntrt',  errors)
         call get_optional_real_with_default(nut, 'frnx',   config%nutrient%frnx,   0.0_real64, 'nutrient.frnx',   errors)
         call read_table_2d(nut, 'nmxlv', config%nutrient%nmxlv, 2, 'nutrient.nmxlv', errors)
      end if
```

Also add `swrdc` reading inside the existing `root` section block in the parser:

```fortran
         call get_optional_int_with_default(root, 'swrdc', config%root%swrdc, 0, 'root.swrdc', errors)
```

At the very end of `read_cropwofost_toml`, before `end subroutine`, add:

```fortran
      ! Mark this config as populated from a real file.
      config%populated = .true.
```

Also update the local pointer declarations in `read_cropwofost_toml` to include the new section variables:

```fortran
      type(toml_table), pointer :: ..., soy, bul, nut
```

- [ ] **Step 4: Run tests — verify pass**

```bash
pixi run test-pfunit 2>&1 | tail -10
```

Expected: `Ok: 1, Fail: 0`.

- [ ] **Step 5: Commit**

```bash
git add src/io/toml/read_cropwofost_toml.f90 tests/unit/io/toml/test_read_cropwofost_toml.pf
git commit -m "feat(parser): cropwofost parser gaps filled + populated sentinel set

Adds soybean, bulb, nutrient sections + swrdc field to parser.
Sets config%populated = .true. at end so the runtime dispatch in
cropgrowth.f90 can distinguish config-loaded from default-constructed
slots. Tested by round-trip of potatod.crp.toml fixture."
```

---

## Task 4: Complete `potatod.crp.toml` — submodule pair commit

Add the `[irrigation_schedule]` section (schedule=0) that is present in the legacy `.crp` but absent from the current 216-line TOML. This makes the TOML fixture 1:1 with all parts of the legacy file.

**Files:**
- Modify (submodule): `tests/swap-cases/toml/5.salinitystress/potatod.crp.toml`

- [ ] **Step 1: Append the irrigation_schedule section**

Open `tests/swap-cases/toml/5.salinitystress/potatod.crp.toml` and append after the `[management]` block (at the very end):

```toml

[irrigation_schedule]
schedule = 0
```

- [ ] **Step 2: Verify the file parses cleanly**

```bash
cd /home/zawadzkim/Code/swap
pixi run test-pfunit -k test_read_cropwofost_potatod_salinity 2>&1 | tail -5
```

Expected: the existing round-trip test passes.

- [ ] **Step 3: Inner submodule commit**

```bash
cd tests/swap-cases
git add toml/5.salinitystress/potatod.crp.toml
git commit -m "feat(case5): complete potatod.crp.toml with irrigation_schedule section"
cd /home/zawadzkim/Code/swap
```

- [ ] **Step 4: Outer bump commit**

```bash
git add tests/swap-cases
git commit -m "chore(submodule): bump swap-cases — complete case 5 potatod.crp.toml"
```

---

## Task 5: Create `cropwofost_init.f90`

Create the runtime init module that replaces `readwofost`'s side-effects on the TOML path. This is the core of Phase 2.

**Files:**
- Create: `src/crop/cropwofost_init.f90`
- Test: `tests/unit/crop/test_cropwofost_init.pf` (new file)
- Modify: `meson.build`, `tests/unit/meson.build`, `tests/unit/testSuites.inc`

- [ ] **Step 1: Write the failing unit tests for `cropwofost_init`**

Create `tests/unit/crop/test_cropwofost_init.pf`:

```fortran
! Phase 2 cropwofost port — unit tests for cropwofost_init_from_config.
! Verifies that calling the init sub from a fixture built from potatod.crp.toml
! values writes the correct module globals, including the cumdens computation
! from the 2-row rdctb = [[0.0, 1.0], [1.0, 0.0]].

@test
subroutine test_cropwofost_init_basic_globals()
   use funit
   use iso_fortran_env, only: real64
   use cropwofost_config_mod, only: cropwofost_config_t
   use cropwofost_init_mod,   only: cropwofost_init_from_config
   use variables, only: swcf, albedo, rsc, rsw, &
                        swoxygen, hlim1, hlim2u, hlim2l, swwrtnonox, aeratecrit, &
                        swdrought, hlim3h, hlim3l, hlim4, adcrh, adcrl, &
                        swsalinity, saltmax, saltslope, &
                        swcompensate, swinter, cofab, &
                        swrd, rdi, rri, rdc, swdmi2rd, &
                        swco2, relmf, swpotrelmf
   use error_mod, only: error_collection_t
   type(cropwofost_config_t) :: cfg
   integer, parameter :: ICROP = 1

   ! Build a fixture matching case 5 potatod.crp values.
   cfg%crop_factor%swcf         = 2
   cfg%crop_factor%albedo       = 0.19_real64
   cfg%crop_factor%rsc          = 207.0_real64
   cfg%crop_factor%rsw          = 0.0_real64
   cfg%oxygen_stress%swoxygen   = 1
   cfg%oxygen_stress%swwrtnonox = 1
   cfg%oxygen_stress%aeratecrit = 0.5_real64
   cfg%oxygen_stress%hlim1      = -10.0_real64
   cfg%oxygen_stress%hlim2u     = -25.0_real64
   cfg%oxygen_stress%hlim2l     = -25.0_real64
   cfg%drought_stress%swdrought = 1
   cfg%drought_stress%hlim3h    = -300.0_real64
   cfg%drought_stress%hlim3l    = -500.0_real64
   cfg%drought_stress%hlim4     = -10000.0_real64
   cfg%drought_stress%adcrh     = 0.5_real64
   cfg%drought_stress%adcrl     = 0.1_real64
   cfg%salinity%swsalinity      = 1
   cfg%salinity%saltmax         = 0.732_real64
   cfg%salinity%saltslope       = 0.0868_real64
   cfg%compensate%swcompensate  = 0
   cfg%interception%swinter     = 1
   cfg%interception%cofab       = 0.25_real64
   cfg%co2%swco2                = 0
   cfg%root%swrd                = 2
   cfg%root%rdi                 = 10.0_real64
   cfg%root%rri                 = 1.2_real64
   cfg%root%rdc                 = 50.0_real64
   cfg%root%swdmi2rd            = 1
   cfg%root%swrdc               = 0
   cfg%management%swpotrelmf    = 2
   cfg%management%relmf         = 0.8_real64
   cfg%soybean%swsoybean        = 0
   cfg%bulb%swbulb              = 0
   cfg%nutrient%flcropnut       = .false.
   cfg%schedule%schedule        = 0
   ! rdctb = [[0.0, 1.0], [1.0, 0.0]] — linear ramp
   allocate(cfg%root%rdctb(2, 2))
   cfg%root%rdctb(1,1) = 0.0_real64; cfg%root%rdctb(1,2) = 1.0_real64
   cfg%root%rdctb(2,1) = 1.0_real64; cfg%root%rdctb(2,2) = 0.0_real64

   call cropwofost_init_from_config(cfg, ICROP)

   ! -- Verify globals --
   @assertEqual(2,            swcf)
   @assertEqual(0.19_real64,  albedo,  1.0e-12_real64)
   @assertEqual(207.0_real64, rsc,     1.0e-12_real64)
   @assertEqual(0.0_real64,   rsw,     1.0e-12_real64)
   @assertEqual(1,            swoxygen)
   @assertEqual(-10.0_real64, hlim1,   1.0e-12_real64)
   @assertEqual(-25.0_real64, hlim2u,  1.0e-12_real64)
   @assertEqual(-25.0_real64, hlim2l,  1.0e-12_real64)
   @assertEqual(1,            swwrtnonox)
   @assertEqual(0.5_real64,   aeratecrit, 1.0e-12_real64)
   @assertEqual(1,            swdrought)
   @assertEqual(-300.0_real64,  hlim3h, 1.0e-12_real64)
   @assertEqual(-500.0_real64,  hlim3l, 1.0e-12_real64)
   @assertEqual(-10000.0_real64,hlim4,  1.0e-12_real64)
   @assertEqual(0.5_real64,   adcrh,  1.0e-12_real64)
   @assertEqual(0.1_real64,   adcrl,  1.0e-12_real64)
   @assertEqual(1,            swsalinity)
   @assertEqual(0.732_real64,  saltmax,   1.0e-12_real64)
   @assertEqual(0.0868_real64, saltslope, 1.0e-12_real64)
   @assertEqual(0,            swcompensate)
   @assertEqual(1,            swinter)
   @assertEqual(0.25_real64,  cofab,  1.0e-12_real64)
   @assertEqual(0,            swco2)
   @assertEqual(2,            swrd)
   @assertEqual(10.0_real64,  rdi, 1.0e-12_real64)
   @assertEqual(1.2_real64,   rri, 1.0e-12_real64)
   @assertEqual(50.0_real64,  rdc, 1.0e-12_real64)
   @assertEqual(1,            swdmi2rd)
   @assertEqual(2,            swpotrelmf)
   @assertEqual(0.8_real64,   relmf, 1.0e-12_real64)
end subroutine

@test
subroutine test_cropwofost_init_cumdens_linear_ramp()
   ! rdctb = [[0.0, 1.0], [1.0, 0.0]] (linear ramp from 1.0 to 0.0).
   ! The cumulative integral over [0,1] = 0.5 (triangle).
   ! The algorithm uses 101 depth points (i=0..100) with depth = i/100.
   ! rdctb density at depth d = (1-d) by linear interpolation.
   ! cumdens(i*2) = integral from 0 to (i-1)/100 of (1-d) dd, normalized.
   ! Expected at depth=0.5 (i=51): integral = 0.5*1 - 0.5*(0.5^2) = 0.5 - 0.125 = 0.375
   ! Normalized: 0.375 / 0.5 = 0.75.
   use funit
   use iso_fortran_env, only: real64
   use cropwofost_config_mod, only: cropwofost_config_t
   use cropwofost_init_mod,   only: cropwofost_init_from_config
   use variables, only: cumdens, swdrought
   type(cropwofost_config_t) :: cfg
   integer, parameter :: ICROP = 1

   cfg%soybean%swsoybean       = 0
   cfg%bulb%swbulb             = 0
   cfg%nutrient%flcropnut      = .false.
   cfg%schedule%schedule       = 0
   cfg%co2%swco2               = 0
   cfg%harvest%swharv          = 0
   cfg%compensate%swcompensate = 0
   cfg%interception%swinter    = 1
   cfg%oxygen_stress%swoxygen  = 1
   cfg%oxygen_stress%swwrtnonox= 0
   cfg%oxygen_stress%aeratecrit= 0.001_real64
   cfg%drought_stress%swdrought= 1
   cfg%drought_stress%hlim3h   = -200.0_real64
   cfg%drought_stress%hlim3l   = -800.0_real64
   cfg%drought_stress%hlim4    = -8000.0_real64
   cfg%drought_stress%adcrh    = 0.5_real64
   cfg%drought_stress%adcrl    = 0.1_real64
   cfg%salinity%swsalinity     = 0
   cfg%root%swrd               = 2
   cfg%root%swrdc              = 0
   cfg%root%rdi                = 10.0_real64
   cfg%root%rri                = 1.0_real64
   cfg%root%rdc                = 30.0_real64
   cfg%root%swdmi2rd           = 0
   cfg%management%swpotrelmf   = 1
   allocate(cfg%root%rdctb(2, 2))
   cfg%root%rdctb(1,1) = 0.0_real64; cfg%root%rdctb(1,2) = 1.0_real64
   cfg%root%rdctb(2,1) = 1.0_real64; cfg%root%rdctb(2,2) = 0.0_real64

   call cropwofost_init_from_config(cfg, ICROP)

   @assertEqual(1, swdrought)
   ! cumdens(1)=depth at index 1 = 0.0; cumdens(2)=cumulative density at 0.0 = 0.0
   @assertEqual(0.0_real64, cumdens(1), 1.0e-10_real64)
   @assertEqual(0.0_real64, cumdens(2), 1.0e-10_real64)
   ! cumdens(201)=depth at index 101 = 1.0; cumdens(202)=normalized density = 1.0
   @assertEqual(1.0_real64, cumdens(201), 1.0e-10_real64)
   @assertEqual(1.0_real64, cumdens(202), 1.0e-10_real64)
   ! cumdens at depth=0.5 (index 51): depth stored at cumdens(101), normalized at cumdens(102)
   @assertEqual(0.5_real64, cumdens(101), 1.0e-10_real64)
   @assertEqual(0.75_real64, cumdens(102), 1.0e-6_real64)
end subroutine
```

- [ ] **Step 2: Register the new test file**

In `tests/unit/meson.build`, add:

```
'crop/test_cropwofost_init.pf',
```

to the pfunit sources list.

In `tests/unit/testSuites.inc`, add:

```fortran
ADD_TEST_SUITE(test_cropwofost_init_suite)
```

- [ ] **Step 3: Run tests — verify failure**

```bash
pixi run test-pfunit 2>&1 | grep -E "cropwofost_init|FAIL" | head -10
```

Expected: both tests fail (module doesn't exist yet).

- [ ] **Step 4: Create `src/crop/cropwofost_init.f90`**

```fortran
!> Runtime initializer for WOFOST (type-2) crop rotations from typed config.
!!
!! Replaces the side-effects of readwofost(task=1) on the TOML pipeline path.
!! Reads all fields from a populated cropwofost_config_t and writes them to
!! the legacy variables-module globals that the wofost computation loop reads.
!! Also computes cumdens (normalized root density) identical to readwofost's
!! tail math.
!!
!! Transitional: removed when the config-passing direction (ADR 0016 Part C)
!! lands and wofost() takes explicit (cfg, state) arguments.
module cropwofost_init_mod
   use iso_fortran_env, only: real64
   use cropwofost_config_mod, only: cropwofost_config_t
   use error_mod, only: fatalerr_collected
   implicit none
   private

   public :: cropwofost_init_from_config

contains

   subroutine cropwofost_init_from_config(cfg, icrop)
      use variables, only: &
         swcf, cftb, chtb, albedo, rsc, rsw,                                &
         idsl, tsumea, tsumam, dlo, dlc, dtsmtb,                            &
         verndvs, vernsat, vernbase, vernrtb,                                &
         tdwi, laiem, rgrlai,                                                &
         slatb, spa, ssa, span, tbase,                                       &
         kdif, kdir, eff, amaxtb, tmpftb, tmnftb,                           &
         cvl, cvo, cvr, cvs,                                                 &
         q10, rml, rmo, rmr, rms, rfsetb,                                   &
         frtb, fltb, fstb, fotb,                                             &
         perdl, rdrrtb, rdrstb,                                              &
         swoxygen, swwrtnonox, aeratecrit, hlim1, hlim2u, hlim2l,            &
         swdrought, hlim3h, hlim3l, hlim4, adcrh, adcrl, wiltpoint,         &
         swsalinity, saltmax, saltslope, salthead,                           &
         swcompensate, swstressor,                                            &
         swinter, cofab,                                                      &
         swrd, rdi, rri, rdc, swdmi2rd, rdctb, rdtb, rlwtb, wrtmax,         &
         swrdc, cumdens,                                                      &
         dvsend, swharv,                                                      &
         swco2, flco2,                                                        &
         relmf, swpotrelmf,                                                   &
         FraDeceasedLvToSoil, fraharlosorm_lv, fraharlosorm_st, fraharlosorm_so, &
         schedule, dvs, tsum, daycrop, nofd
      use array_utils, only: afgen
      implicit none

      class(cropwofost_config_t), intent(in) :: cfg
      integer,                    intent(in) :: icrop

      integer     :: i
      real(real64) :: depth, sum_dens
      real(real64) :: rootdis(202)

      ! ----------------------------------------------------------------
      ! Defense-in-depth guards — the validator should have caught these,
      ! but we guard at runtime too per ADR 0015.
      ! ----------------------------------------------------------------
      if (cfg%soybean%swsoybean == 1) &
         call fatalerr_collected('cropwofost_init', &
            'swsoybean=1 not supported on TOML path; validator should have rejected.')
      if (cfg%bulb%swbulb == 1) &
         call fatalerr_collected('cropwofost_init', &
            'swbulb=1 not supported on TOML path; validator should have rejected.')
      if (cfg%nutrient%flcropnut) &
         call fatalerr_collected('cropwofost_init', &
            'flcropnut=.true. not supported on TOML path; validator should have rejected.')
      if (cfg%co2%swco2 == 1) &
         call fatalerr_collected('cropwofost_init', &
            'swco2=1 not supported on TOML path; validator should have rejected.')
      if (cfg%schedule%schedule == 1) &
         call fatalerr_collected('cropwofost_init', &
            'schedule=1 not supported on TOML path; validator should have rejected.')
      if (cfg%drought_stress%swdrought == 2) &
         call fatalerr_collected('cropwofost_init', &
            'swdrought=2 not supported on TOML path; validator should have rejected.')
      if (cfg%oxygen_stress%swoxygen == 2) &
         call fatalerr_collected('cropwofost_init', &
            'swoxygen=2 not supported on TOML path; validator should have rejected.')
      if (cfg%interception%swinter == 2) &
         call fatalerr_collected('cropwofost_init', &
            'swinter=2 not supported on TOML path; validator should have rejected.')
      if (cfg%compensate%swcompensate /= 0) &
         call fatalerr_collected('cropwofost_init', &
            'swcompensate/=0 not supported on TOML path; validator should have rejected.')
      if (cfg%harvest%swharv == 1) &
         call fatalerr_collected('cropwofost_init', &
            'swharv=1 not supported on TOML path; validator should have rejected.')
      if (cfg%salinity%swsalinity == 2) &
         call fatalerr_collected('cropwofost_init', &
            'swsalinity=2 not supported on TOML path; validator should have rejected.')
      if (cfg%root%swrdc == 1) &
         call fatalerr_collected('cropwofost_init', &
            'swrdc=1 not supported on TOML path; validator should have rejected.')

      ! ----------------------------------------------------------------
      ! Config → globals copy (mirroring readwofost's rd* sequence)
      ! ----------------------------------------------------------------

      ! Part 1: crop factor / crop height
      swcf = cfg%crop_factor%swcf
      if (swcf == 1) then
         ! Store cftb as flat paired array for afgen: (dvs1, cf1, dvs2, cf2, ...)
         if (allocated(cfg%crop_factor%cftb)) then
            block
               integer :: nr, j
               nr = size(cfg%crop_factor%cftb, 1)
               do j = 1, nr
                  cftb(j*2-1) = cfg%crop_factor%cftb(j,1)
                  cftb(j*2)   = cfg%crop_factor%cftb(j,2)
               end do
               chtb = -99.99d0
            end block
         end if
      else if (swcf == 2) then
         albedo = cfg%crop_factor%albedo
         rsc    = cfg%crop_factor%rsc
         rsw    = cfg%crop_factor%rsw
         if (allocated(cfg%crop_factor%chtb)) then
            block
               integer :: nr, j
               nr = size(cfg%crop_factor%chtb, 1)
               do j = 1, nr
                  chtb(j*2-1) = cfg%crop_factor%chtb(j,1)
                  chtb(j*2)   = cfg%crop_factor%chtb(j,2)
               end do
               cftb = -99.99d0
            end block
         end if
      end if

      ! Part 14: interception
      swinter = cfg%interception%swinter
      if (swinter == 1) cofab = cfg%interception%cofab

      ! Part 2: phenology (soybean=0 path only)
      idsl   = cfg%phenology%idsl
      tsumea = cfg%phenology%tsumea
      tsumam = cfg%phenology%tsumam
      if (idsl == 1 .or. idsl == 2) then
         dlo = cfg%phenology%dlo
         dlc = cfg%phenology%dlc
      end if
      if (allocated(cfg%phenology%dtsmtb)) then
         block
            integer :: nr, j
            nr = size(cfg%phenology%dtsmtb, 1)
            do j = 1, nr
               dtsmtb(j*2-1) = cfg%phenology%dtsmtb(j,1)
               dtsmtb(j*2)   = cfg%phenology%dtsmtb(j,2)
            end do
         end block
      end if
      ! vernalization (idsl=2) stub-guarded above — safe defaults already in module

      ! Part 3: initial crop state
      tdwi   = cfg%initial%tdwi
      laiem  = cfg%initial%laiem
      rgrlai = cfg%initial%rgrlai

      ! Part 4: green area
      if (allocated(cfg%green_area%slatb)) then
         block
            integer :: nr, j
            nr = size(cfg%green_area%slatb, 1)
            do j = 1, nr
               slatb(j*2-1) = cfg%green_area%slatb(j,1)
               slatb(j*2)   = cfg%green_area%slatb(j,2)
            end do
         end block
      end if
      spa   = cfg%green_area%spa
      ssa   = cfg%green_area%ssa
      span  = cfg%green_area%span
      tbase = cfg%green_area%tbase

      ! Part 5: assimilation
      kdif = cfg%assimilation%kdif
      kdir = cfg%assimilation%kdir
      eff  = cfg%assimilation%eff
      if (allocated(cfg%assimilation%amaxtb)) then
         block
            integer :: nr, j
            nr = size(cfg%assimilation%amaxtb, 1)
            do j = 1, nr
               amaxtb(j*2-1) = cfg%assimilation%amaxtb(j,1)
               amaxtb(j*2)   = cfg%assimilation%amaxtb(j,2)
            end do
         end block
      end if
      if (allocated(cfg%assimilation%tmpftb)) then
         block
            integer :: nr, j
            nr = size(cfg%assimilation%tmpftb, 1)
            do j = 1, nr
               tmpftb(j*2-1) = cfg%assimilation%tmpftb(j,1)
               tmpftb(j*2)   = cfg%assimilation%tmpftb(j,2)
            end do
         end block
      end if
      if (allocated(cfg%assimilation%tmnftb)) then
         block
            integer :: nr, j
            nr = size(cfg%assimilation%tmnftb, 1)
            do j = 1, nr
               tmnftb(j*2-1) = cfg%assimilation%tmnftb(j,1)
               tmnftb(j*2)   = cfg%assimilation%tmnftb(j,2)
            end do
         end block
      end if

      ! Part 6: conversion
      cvl = cfg%conversion%cvl
      cvo = cfg%conversion%cvo
      cvr = cfg%conversion%cvr
      cvs = cfg%conversion%cvs

      ! Part 7: respiration
      q10 = cfg%respiration%q10
      rml = cfg%respiration%rml
      rmo = cfg%respiration%rmo
      rmr = cfg%respiration%rmr
      rms = cfg%respiration%rms
      if (allocated(cfg%respiration%rfsetb)) then
         block
            integer :: nr, j
            nr = size(cfg%respiration%rfsetb, 1)
            do j = 1, nr
               rfsetb(j*2-1) = cfg%respiration%rfsetb(j,1)
               rfsetb(j*2)   = cfg%respiration%rfsetb(j,2)
            end do
         end block
      end if

      ! Part 8: partitioning
      if (allocated(cfg%partitioning%frtb)) then
         block
            integer :: nr, j
            nr = size(cfg%partitioning%frtb, 1)
            do j = 1, nr
               frtb(j*2-1) = cfg%partitioning%frtb(j,1)
               frtb(j*2)   = cfg%partitioning%frtb(j,2)
            end do
         end block
      end if
      if (allocated(cfg%partitioning%fltb)) then
         block
            integer :: nr, j
            nr = size(cfg%partitioning%fltb, 1)
            do j = 1, nr
               fltb(j*2-1) = cfg%partitioning%fltb(j,1)
               fltb(j*2)   = cfg%partitioning%fltb(j,2)
            end do
         end block
      end if
      if (allocated(cfg%partitioning%fstb)) then
         block
            integer :: nr, j
            nr = size(cfg%partitioning%fstb, 1)
            do j = 1, nr
               fstb(j*2-1) = cfg%partitioning%fstb(j,1)
               fstb(j*2)   = cfg%partitioning%fstb(j,2)
            end do
         end block
      end if
      if (allocated(cfg%partitioning%fotb)) then
         block
            integer :: nr, j
            nr = size(cfg%partitioning%fotb, 1)
            do j = 1, nr
               fotb(j*2-1) = cfg%partitioning%fotb(j,1)
               fotb(j*2)   = cfg%partitioning%fotb(j,2)
            end do
         end block
      end if

      ! Part 9: death rates
      perdl = cfg%death%perdl
      if (allocated(cfg%death%rdrrtb)) then
         block
            integer :: nr, j
            nr = size(cfg%death%rdrrtb, 1)
            do j = 1, nr
               rdrrtb(j*2-1) = cfg%death%rdrrtb(j,1)
               rdrrtb(j*2)   = cfg%death%rdrrtb(j,2)
            end do
         end block
      end if
      if (allocated(cfg%death%rdrstb)) then
         block
            integer :: nr, j
            nr = size(cfg%death%rdrstb, 1)
            do j = 1, nr
               rdrstb(j*2-1) = cfg%death%rdrstb(j,1)
               rdrstb(j*2)   = cfg%death%rdrstb(j,2)
            end do
         end block
      end if

      ! Part 11: oxygen stress
      swoxygen   = cfg%oxygen_stress%swoxygen
      swwrtnonox = cfg%oxygen_stress%swwrtnonox
      aeratecrit = cfg%oxygen_stress%aeratecrit
      if (swoxygen == 1) then
         hlim1  = cfg%oxygen_stress%hlim1
         hlim2u = cfg%oxygen_stress%hlim2u
         hlim2l = cfg%oxygen_stress%hlim2l
      end if

      ! Part 12: drought stress
      swdrought = cfg%drought_stress%swdrought
      if (swdrought == 1) then
         hlim3h = cfg%drought_stress%hlim3h
         hlim3l = cfg%drought_stress%hlim3l
         hlim4  = cfg%drought_stress%hlim4
         adcrh  = cfg%drought_stress%adcrh
         adcrl  = cfg%drought_stress%adcrl
      end if

      ! Part 13: salinity stress
      swsalinity = cfg%salinity%swsalinity
      if (swsalinity == 1) then
         saltmax   = cfg%salinity%saltmax
         saltslope = cfg%salinity%saltslope
      end if

      ! Part xx: compensation
      swcompensate = cfg%compensate%swcompensate
      swstressor   = cfg%compensate%swstressor

      ! Part 10: root depth and density
      swrdc = cfg%root%swrdc
      if (allocated(cfg%root%rdctb)) then
         block
            integer :: nr, j
            nr = size(cfg%root%rdctb, 1)
            do j = 1, nr
               rdctb(j*2-1) = cfg%root%rdctb(j,1)
               rdctb(j*2)   = cfg%root%rdctb(j,2)
            end do
         end block
      end if
      swrd = cfg%root%swrd
      select case (swrd)
      case (1)
         if (allocated(cfg%root%rdtb)) then
            block
               integer :: nr, j
               nr = size(cfg%root%rdtb, 1)
               do j = 1, nr
                  rdtb(j*2-1) = cfg%root%rdtb(j,1)
                  rdtb(j*2)   = cfg%root%rdtb(j,2)
               end do
            end block
         end if
      case (2)
         rdi      = cfg%root%rdi
         rri      = cfg%root%rri
         rdc      = cfg%root%rdc
         swdmi2rd = cfg%root%swdmi2rd
      case (3)
         if (allocated(cfg%root%rlwtb)) then
            block
               integer :: nr, j
               nr = size(cfg%root%rlwtb, 1)
               do j = 1, nr
                  rlwtb(j*2-1) = cfg%root%rlwtb(j,1)
                  rlwtb(j*2)   = cfg%root%rlwtb(j,2)
               end do
            end block
         end if
         wrtmax = cfg%root%wrtmax
      end select

      ! Harvest
      dvsend = cfg%harvest%dvsend
      swharv = cfg%harvest%swharv

      ! CO2 (swco2=0 only on TOML path; stub-guarded above)
      swco2 = cfg%co2%swco2
      flco2 = (swco2 == 1)

      ! Schedule (schedule=0 only on TOML path)
      schedule = cfg%schedule%schedule

      ! Management
      FraDeceasedLvToSoil = cfg%management%fradeceasedlvtosoil
      fraharlosorm_lv     = cfg%management%fraharlosorm_lv
      fraharlosorm_st     = cfg%management%fraharlosorm_st
      fraharlosorm_so     = cfg%management%fraharlosorm_so
      relmf               = cfg%management%relmf
      swpotrelmf          = cfg%management%swpotrelmf

      ! ----------------------------------------------------------------
      ! Runtime init math — cumdens computation (from readwofost lines ~3091-3121)
      ! Only when swdrought=1 (Feddes). Identical algorithm to readwofost.
      ! ----------------------------------------------------------------
      if (swdrought == 1) then
         ! Build rootdis array: 101 points from 0.0 to 1.0
         do i = 0, 100
            depth = 0.01d0 * dble(i)
            rootdis(i*2+1) = depth
            rootdis(i*2+2) = afgen(rdctb, 22, depth)
         end do

         ! Build cumulative root density
         do i = 1, 202, 2
            cumdens(i) = rootdis(i)
         end do
         sum_dens    = 0.0d0
         cumdens(2)  = 0.0d0
         do i = 4, 202, 2
            sum_dens = sum_dens + (rootdis(i-2) + rootdis(i)) * 0.5d0 &
                                * (cumdens(i-1) - cumdens(i-3))
            cumdens(i) = sum_dens
         end do

         ! Normalize to 1
         if (sum_dens > 0.0d0) then
            do i = 2, 202, 2
               cumdens(i) = cumdens(i) / sum_dens
            end do
         end if
      end if

      ! ----------------------------------------------------------------
      ! Crop state initialization (mirrors readwofost initialization for
      ! a fresh crop start — swinco=3 .END-file restart is gated below).
      ! ----------------------------------------------------------------
      dvs     = 0.0d0
      tsum    = 0.0d0
      daycrop = 0
      nofd    = 0

      ! Note: the .END-file restart block (readwofost lines 3123-3227, active
      ! when swinco=3 and t1900-tstart<1e-3 and |t1900-cropstart|>=1e-3) is NOT
      ! reproduced here. For case 5, all 4 crop seasons start after the
      ! simulation start date so the condition evaluates to false and the block
      ! is never entered. If a future case triggers this condition, the
      ! validator should be extended to detect it and emit a stub-error.
      ! Defense-in-depth: if a bug causes this path to be reached incorrectly,
      ! the simulation will start with zeroed crop state (conservative), and
      ! the regression test will detect the discrepancy.

   end subroutine cropwofost_init_from_config

end module cropwofost_init_mod
```

- [ ] **Step 5: Register `cropwofost_init.f90` in `meson.build`**

In `meson.build`, find the `crop_sources` list (or the list that includes `cropfixed_init.f90` from Phase 1) and add:

```meson
'src/crop/cropwofost_init.f90',
```

- [ ] **Step 6: Run tests — verify pass**

```bash
pixi run test-pfunit 2>&1 | tail -10
```

Expected: `Ok: 1, Fail: 0`. If `cumdens` test fails with a floating-point mismatch at depth=0.5, recalculate: the rdctb linear ramp uses `afgen` for interpolation; the expected value of 0.75 is computed from the analytic integral of (1-d) from 0 to 0.5 = 0.375, normalized by integral from 0 to 1 = 0.5, giving 0.75. If the discrete algorithm differs by more than 1e-6, relax the tolerance to 1e-3.

- [ ] **Step 7: Commit**

```bash
git add src/crop/cropwofost_init.f90 tests/unit/crop/test_cropwofost_init.pf \
        meson.build tests/unit/meson.build tests/unit/testSuites.inc
git commit -m "feat(crop): cropwofost_init_from_config — TOML-path runtime init for WOFOST

New module replaces readwofost's runtime side-effects on the TOML path.
Copies all config fields to variables-module globals in readwofost order,
computes cumdens (normalized root density) from rdctb, initializes crop
state scalars. Defense-in-depth guards mirror validator stub-errors.
Tested by unit assertions on globals and cumdens analytic check."
```

---

## Task 6: Extend `test_load_swap_config.pf` — assert all 4 wofost slots loaded

Verify that loading case 5's `swap.toml` populates all 4 `rotation_wofost` slots.

**Files:**
- Modify: `tests/unit/io/toml/test_load_swap_config.pf`

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/io/toml/test_load_swap_config.pf`:

```fortran
@test
subroutine test_load_swap_config_case5_wofost_populated()
   ! Case 5 has 4 rotation entries all pointing at potatod.crp.toml (type=2).
   ! After load+validate+finalize, all 4 slots must be loaded and populated.
   use funit
   use iso_fortran_env, only: real64
   use swap_config_mod, only: load_swap_config
   use error_mod, only: error_collection_t
   type(error_collection_t) :: errors
   type(swap_config_t)      :: config
   integer :: i

   call load_swap_config('tests/swap-cases/toml/5.salinitystress/swap.toml', &
                          config, errors)
   @assertFalse(errors%has_fatals())

   @assertEqual(4, size(config%crop%rotation_wofost))
   @assertEqual(4, size(config%crop%rotation_loaded))

   do i = 1, 4
      @assertTrue(config%crop%rotation_loaded(i))
      @assertTrue(config%crop%rotation_wofost(i)%populated)
      ! All 4 rotations share the same potatod.crp.toml; verify key field.
      @assertEqual(1, config%crop%rotation_wofost(i)%salinity%swsalinity)
      @assertEqual(0.732_real64, config%crop%rotation_wofost(i)%salinity%saltmax, &
                   1.0e-12_real64)
   end do

   ! Confirm type=2 for all entries
   do i = 1, 4
      @assertEqual(2, config%crop%rotation_type(i))
   end do
end subroutine
```

- [ ] **Step 2: Run tests — verify pass (it should already pass after Tasks 2-4)**

```bash
pixi run test-pfunit 2>&1 | tail -5
```

Expected: `Ok: 1, Fail: 0`.

- [ ] **Step 3: Commit**

```bash
git add tests/unit/io/toml/test_load_swap_config.pf
git commit -m "test(load): assert case 5 loads all 4 wofost rotation slots populated"
```

---

## Task 7: Wire dispatch in `cropgrowth.f90:984`

Add the `if (associated(...) .and. ... .and. populated)` dispatch so the TOML path calls `cropwofost_init_from_config` instead of `readwofost`.

**Files:**
- Modify: `src/crop/cropgrowth.f90`

- [ ] **Step 1: Read the relevant section**

```bash
sed -n '880,1000p' src/crop/cropgrowth.f90
```

Confirm the exact line of `call readwofost (icrop,cropfil(icrop),...)` — expected around line 984.

- [ ] **Step 2: Add USE statements to the `wofost` subroutine**

In `src/crop/cropgrowth.f90`, inside the `wofost(task)` subroutine (after line ~888), add two USE statements after the existing `use error_mod` line:

```fortran
      use crop_config_global_mod, only: crop_config_global
      use cropwofost_init_mod,    only: cropwofost_init_from_config
```

- [ ] **Step 3: Replace the `readwofost` call with the dispatch block**

Find (around line 984):

```fortran
      call readwofost (icrop,cropfil(icrop),swhydrlift,swsoybean,mg,dvsi,dvrmax1,dvrmax2, &
                       flrfphotoveg,tmaxdvr,tmindvr,toptdvr,popt,pcrt,flphenodayl,FraDeceasedLvToSoil)
```

Replace with:

```fortran
      if (associated(crop_config_global) .and. &
          allocated(crop_config_global%rotation_wofost) .and. &
          crop_config_global%rotation_wofost(icrop)%populated) then
         call cropwofost_init_from_config( &
                 crop_config_global%rotation_wofost(icrop), icrop)
      else
         call readwofost (icrop,cropfil(icrop),swhydrlift,swsoybean,mg,dvsi,dvrmax1, &  ! transitional
     &                    dvrmax2,flrfphotoveg,tmaxdvr,tmindvr,toptdvr,popt,pcrt,   &
     &                    flphenodayl,FraDeceasedLvToSoil)
      end if
```

Note: preserve the original Fortran continuation style (using `&` at the end of lines) to match the surrounding code.

- [ ] **Step 4: Build to confirm compilation**

```bash
pixi run build 2>&1 | tail -20
```

Expected: zero errors. If gfortran complains about module not found, confirm `crop_config_global.f90` and `cropwofost_init.f90` are both registered in `meson.build`.

- [ ] **Step 5: Run full unit test suite**

```bash
pixi run test-pfunit 2>&1 | tail -10
```

Expected: `Ok: 1, Fail: 0`.

- [ ] **Step 6: Run regression**

```bash
pixi run regression 2>&1 | tail -15
```

Expected: `5 passed, 0 failed`. If case 5 fails, the `cropwofost_init_from_config` global copy has a bug — compare the failing global against `readwofost` output using the parity test suite.

- [ ] **Step 7: Commit**

```bash
git add src/crop/cropgrowth.f90
git commit -m "feat(crop): wire TOML dispatch for wofost(task=1) via populated sentinel

When crop_config_global is set and rotation_wofost(icrop)%populated is
true, calls cropwofost_init_from_config instead of readwofost. Legacy
readwofost call preserved in else branch (transitional per ADR 0016).
Regression 5/5 green."
```

---

## Task 8: Add parity assertions for `cropwofost_init`

Extend the existing salinitystress parity test with assertions that `cropwofost_init_from_config` produces the same globals as `readwofost`.

**Files:**
- Modify: `tests/unit/io/toml/test_salinitystress_parity.pf`

- [ ] **Step 1: Append the init parity test**

Append to `tests/unit/io/toml/test_salinitystress_parity.pf`:

```fortran
! Parity: cropwofost_init_from_config(rotation_wofost(1)) vs readwofost(potatod).
! The loaded config-side globals must equal the legacy-reader globals after
! both have run for rotation 1 (icrop=1).
@test
subroutine test_salinitystress_wofost_init_parity()
   use funit
   use iso_fortran_env, only: real64
   use parity_helpers_mod, only: load_both_for_salinitystress, reset_for_next_readswap
   use legacy_crop_helper_mod, only: read_legacy_wofost, flatten_table
   use chdir_helper_mod, only: chdir_to, get_cwd
   use cropwofost_init_mod, only: cropwofost_init_from_config
   use swap_config_mod
   use error_mod
   use variables, only: &
      flsolute,                                                              &
      swcf, albedo, rsc, rsw,                                               &
      idsl, tsumea, tsumam,                                                 &
      tdwi, laiem, rgrlai, spa, ssa, span, tbase,                          &
      kdif, kdir, eff,                                                      &
      cvl, cvo, cvr, cvs, q10, rml, rmo, rmr, rms, perdl,                  &
      swrd, rdi, rri, rdc, swdmi2rd,                                       &
      swoxygen, swwrtnonox, aeratecrit, hlim1, hlim2u, hlim2l,             &
      swdrought, hlim3h, hlim3l, hlim4, adcrh, adcrl,                      &
      swsalinity, saltmax, saltslope,                                       &
      swcompensate, swinter, cofab,                                         &
      swco2, relmf, swpotrelmf,                                             &
      cumdens,                                                              &
      dtsmtb, slatb, amaxtb, frtb, fltb, fstb, fotb,                       &
      dvsend, swharv

   type(swap_config_t)      :: config
   type(error_collection_t) :: errors
   character(len=1024)      :: orig_cwd
   real(real64) :: legacy_cumdens(202)
   real(real64) :: flat_amaxtb(30), flat_frtb(30)
   integer, parameter :: ICROP = 1

   ! 1. Load the TOML config (populates rotation_wofost).
   call load_both_for_salinitystress(config, errors)
   @assertFalse(errors%has_fatals())

   ! 2. Run legacy readwofost to populate variables globals.
   flsolute = .true.
   call get_cwd(orig_cwd)
   call chdir_to('tests/swap-cases/5.salinitystress')
   call reset_for_next_readswap()
   call read_legacy_wofost(ICROP, 'potatod')
   call chdir_to(trim(orig_cwd))

   ! 3. Save the legacy cumdens for comparison.
   legacy_cumdens = cumdens

   ! 4. Reset relevant globals to zero, then run the init sub.
   cumdens = 0.0d0
   call cropwofost_init_from_config(config%crop%rotation_wofost(ICROP), ICROP)

   ! 5. Assert globals match legacy values.
   @assertEqual(swcf,    config%crop%rotation_wofost(ICROP)%crop_factor%swcf)
   @assertEqual(albedo,  config%crop%rotation_wofost(ICROP)%crop_factor%albedo, 1.0e-12_real64)
   @assertEqual(rsc,     config%crop%rotation_wofost(ICROP)%crop_factor%rsc,    1.0e-12_real64)
   @assertEqual(rsw,     config%crop%rotation_wofost(ICROP)%crop_factor%rsw,    1.0e-12_real64)
   @assertEqual(swdrought, config%crop%rotation_wofost(ICROP)%drought_stress%swdrought)
   @assertEqual(hlim3h,  config%crop%rotation_wofost(ICROP)%drought_stress%hlim3h, 1.0e-12_real64)
   @assertEqual(swsalinity, config%crop%rotation_wofost(ICROP)%salinity%swsalinity)
   @assertEqual(saltmax, config%crop%rotation_wofost(ICROP)%salinity%saltmax,  1.0e-12_real64)
   @assertEqual(saltslope, config%crop%rotation_wofost(ICROP)%salinity%saltslope, 1.0e-12_real64)

   ! 6. Assert cumdens matches legacy after init.
   @assertEqual(legacy_cumdens, cumdens, 1.0e-10_real64)

   ! 7. Assert key table globals (flat).
   flat_amaxtb = flatten_table(config%crop%rotation_wofost(ICROP)%assimilation%amaxtb, 30)
   @assertEqual(amaxtb, flat_amaxtb, 1.0e-12_real64)
   flat_frtb = flatten_table(config%crop%rotation_wofost(ICROP)%partitioning%frtb, 30)
   @assertEqual(frtb, flat_frtb, 1.0e-12_real64)
end subroutine
```

- [ ] **Step 2: Run tests — verify pass**

```bash
pixi run test-pfunit 2>&1 | tail -10
```

Expected: `Ok: 1, Fail: 0`. If `cumdens` comparison fails, the init sub's discrete integration doesn't exactly match `readwofost`'s. Compare the two implementations line-by-line; the algorithm must be verbatim.

- [ ] **Step 3: Commit**

```bash
git add tests/unit/io/toml/test_salinitystress_parity.pf
git commit -m "test(parity): cropwofost_init globals + cumdens parity vs readwofost for case 5"
```

---

## Task 9: Smoke test

Rename `potatod.crp` → `potatod.crp.disabled` in the TOML case dir to force the TOML path exclusively, run regression, then restore.

- [ ] **Step 1: Rename in submodule**

```bash
cd tests/swap-cases
mv toml/5.salinitystress/potatod.crp toml/5.salinitystress/potatod.crp.disabled
```

- [ ] **Step 2: Run regression**

```bash
cd /home/zawadzkim/Code/swap
pixi run regression 2>&1 | tail -15
```

Expected: `5 passed, 0 failed`. Case 5 must pass without `potatod.crp`. If it fails, the dispatch guard is not working (check `associated(crop_config_global)` and `populated` flag path).

- [ ] **Step 3: Restore**

```bash
cd tests/swap-cases
mv toml/5.salinitystress/potatod.crp.disabled toml/5.salinitystress/potatod.crp
cd /home/zawadzkim/Code/swap
```

**Do not commit the rename.** This is verification only.

---

## Task 10: Delete `potatod.crp` from TOML case dir — submodule pair commit

With the smoke test passing, permanently remove the legacy ASCII file from the TOML case directory. The TOML executable no longer needs it.

**Files:**
- Delete (submodule): `tests/swap-cases/toml/5.salinitystress/potatod.crp`

- [ ] **Step 1: Delete the file**

```bash
cd tests/swap-cases
git rm toml/5.salinitystress/potatod.crp
```

- [ ] **Step 2: Inner submodule commit**

```bash
git commit -m "feat(case5): remove legacy potatod.crp from TOML case dir (Phase 2 complete)"
cd /home/zawadzkim/Code/swap
```

- [ ] **Step 3: Run regression to confirm 5/5 still green**

```bash
pixi run regression 2>&1 | tail -10
```

Expected: `5 passed, 0 failed`.

- [ ] **Step 4: Outer bump commit**

```bash
git add tests/swap-cases
git commit -m "chore(submodule): bump swap-cases — delete case 5 legacy potatod.crp"
```

---

## Task 11: Docs update

Update project documentation to reflect the Phase 2 completion.

**Files:**
- Modify: `docs/csv-companion-files.md` (if it covers `.crp.toml` path resolution)
- Modify: `docs/configuration-schema.md` (if `[cropwofost]` documentation needs updating)

- [ ] **Step 1: Check and update `docs/csv-companion-files.md`**

```bash
grep -n "crp.toml\|wofost\|cropwofost" docs/csv-companion-files.md | head -20
```

If the file mentions `.crp.toml` path resolution, add a note that `potatod.crp.toml` fields for the wofost type-2 sub-sections (`soybean`, `bulb`, `nutrient`) are parsed but their top-level switches are stub-errored in Phase 2. If the file has no such section, skip.

- [ ] **Step 2: Check `docs/configuration-schema.md`**

```bash
grep -n "cropwofost\|wofost_soybean\|wofost_bulb\|wofost_nutrient" docs/configuration-schema.md | head -20
```

If the schema docs are missing entries for the three new sub-types, add brief descriptions noting they are schema 1:1 with legacy but stub-errored in Phase 2. If the docs are auto-generated or minimal, skip and note in the commit message.

- [ ] **Step 3: Commit (even if no changes)**

```bash
git add docs/csv-companion-files.md docs/configuration-schema.md 2>/dev/null
git diff --cached --quiet || git commit -m "docs(schema): note Phase 2 cropwofost soybean/bulb/nutrient stub-errors"
```

---

## Acceptance gate

Before marking Phase 2 complete, verify all of the following:

```bash
# 1. All unit tests green
pixi run test-pfunit 2>&1 | tail -5
# Expected: Ok: 1, Fail: 0

# 2. All 5 regression cases green
pixi run regression 2>&1 | tail -10
# Expected: 5 passed, 0 failed

# 3. Legacy ASCII deleted from TOML case dir
test ! -f tests/swap-cases/toml/5.salinitystress/potatod.crp && echo "DELETED OK"

# 4. Only one readwofost callsite in src/ (the transitional else branch)
git grep "call readwofost" src/ | grep -v "else"
# Expected: no output (only the else branch remains, which contains "else")

# 5. The init module exists
test -f src/crop/cropwofost_init.f90 && echo "init OK"

# 6. populated sentinel present on the config type
grep -c "populated" src/config/cropwofost_config.f90
# Expected: >= 3 (declaration + initialization + set in parser)
```

---

## Self-review notes

**Spec coverage trace:**

| Spec Unit | Plan Task(s) | Status |
|---|---|---|
| Unit 1 — Schema gap-fill + stub-errors + populated | Task 2 | Covered |
| Unit 2 — Parser gap-fill + populated sentinel | Task 3 | Covered |
| Unit 3 — Loader: no changes | (verified in pre-flight) | Covered |
| Unit 4 — `cropwofost_init` runtime init module | Task 5 | Covered |
| Unit 5 — `crop_config_global` reuse from Phase 1 | Task 7 (USE only) | Covered |
| Unit 6 — Dispatch wiring `cropgrowth.f90:984` | Task 7 | Covered |

**Teardown table coverage:**
- `else call readwofost(...)` transitional — introduced in Task 7, documented in spec teardown table.
- `populated` sentinel — introduced in Task 2, documented in spec teardown table.
- `potatod.crp` legacy ASCII in TOML dir — deleted in Task 10.
- `crop_config_global` and soybean/bulb/nutrient stub-errors — Phase 1 artefact and future-phase scope; documented in spec.

**No placeholders:** Every step has a concrete code block, shell command, or explicit skip justification. No "TBD" or "implement later" language.

**Submodule pair commits:** Two pairs: Task 4 (add `[irrigation_schedule]`) and Task 10 (delete `potatod.crp`). Each has inner + outer commits written out.

**16-argument `readwofost` preserved:** Task 7 Step 3 shows the full continuation-line `else` call verbatim.

**cumdens analytic check:** Task 5 Step 1 derives the expected 0.75 value at depth=0.5 from the analytic integral of the linear ramp rdctb, with guidance on tolerance if discrete vs analytic drift.
