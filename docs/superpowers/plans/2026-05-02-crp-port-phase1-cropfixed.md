# `.crp` port — Phase 1 (cropfixed via case 6) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Port case 6 (surfacewater)'s `grass.crp` (type 1, fixed-crop) to TOML so the new executable's runtime path opens no `.crp` ASCII file for type-1 rotations. Case 6 regression remains 5/5 green throughout.

**Architecture:** Mirrors the swap.dra port shape (extend schema 1:1 with legacy fields, push validation/init into typed config + a new init module, wire the runtime callsite to dispatch on the `rotation_loaded` sentinel, fall back to the legacy reader for unported types). Per ADR 0016, the per-rotation crop config cache is the long-term shape; per ADR 0015, runtime branches not yet supported are stub-errored at validation time. The cache infrastructure (`crop_config_t.rotation_fixed(:)`, `rotation_loaded(:)`, dispatch in `read_crop_toml.f90`) is already in place and is reused unchanged.

**Tech Stack:** Fortran 2008 (gfortran 13), pFUnit unit tests, Meson + Pixi, tomlf TOML parser.

**Spec:** [docs/superpowers/specs/2026-05-02-crp-port-phase1-cropfixed-design.md](../specs/2026-05-02-crp-port-phase1-cropfixed-design.md)

**ADR:** [docs/adr/0016-per-rotation-crop-config-cache.md](../../adr/0016-per-rotation-crop-config-cache.md)

**Naming reconciliation (spec → plan):** the spec used proposed names that turned out to be already implemented under different names. The plan uses the actual names:

| Spec name (proposal) | Actual code name |
|---|---|
| `rotation_cropfixed(:)` | `rotation_fixed(:)` (already on `crop_config_t`) |
| `populated :: logical` per-config | `rotation_loaded(:)` LOGICAL array on `crop_config_t` (already implemented) |
| Loader extension in `read_crop_toml.f90` | Already implemented; type=1 branch calls `read_cropfixed_toml` |

The spec's other naming (`crop_config_global`, `cropfixed_init`, `cropfixed_init_from_config`) is not yet implemented; the plan introduces them as written.

---

## File map

**Source — modify:**
- `src/config/cropfixed_config.f90` — extend schema 1:1 with legacy `readcropfixed`; add stub-error validators
- `src/io/toml/read_cropfixed_toml.f90` — extend parser for new fields and tables
- `src/io/toml/config_to_variables.f90` — set the new `crop_config_global` pointer at the end of the crop block
- `src/crop/cropgrowth.f90` — at line 397, dispatch on `rotation_loaded(icrop)` between `cropfixed_init_from_config` and legacy `readcropfixed`
- `meson.build` — register `src/crop/crop_config_global.f90` and `src/crop/cropfixed_init.f90`
- `tests/unit/meson.build` — register new test files
- `tests/unit/testSuites.inc` — register new test suites

**Source — create:**
- `src/crop/crop_config_global.f90` — module-level pointer to the parsed crop config (transitional; teardown when config-passing direction lands)
- `src/crop/cropfixed_init.f90` — runtime init from typed config (replaces `readcropfixed`'s runtime side-effects on the TOML path)

**Tests — create / modify:**
- `tests/unit/config/test_cropfixed_config.pf` — extend with stub-error tests + table-shape tests + new scalar range tests
- `tests/unit/io/toml/test_read_cropfixed_toml.pf` — extend with round-trip of all new scalars and tables
- `tests/unit/io/toml/test_load_swap_config.pf` — extend with case-6 assertion that all three rotation slots are loaded
- `tests/unit/crop/test_cropfixed_init.pf` — NEW (creates `tests/unit/crop/` directory) — verifies init writes correct globals + cumdens hand-computed value

**Docs — modify:**
- `docs/csv-companion-files.md` — note the `.crp.toml` resolution rule alongside `.csv`
- `docs/configuration-schema.md` — extend `[cropfixed]`/`[[crop.rotation]]` documentation if not already covering the new fields

**Test fixture — modify (in submodule):**
- `tests/swap-cases/toml/6.surfacewater/grass.crp.toml` — replace skeletal 33-line file with 1:1 reproduction of legacy 202-line `grass.crp`
- `tests/swap-cases/toml/6.surfacewater/grass.crp` — DELETE in Task 9 (submodule pair commit)

---

## Submodule discipline

`tests/swap-cases/` is a git submodule. Inner-commit + outer-bump pair is non-negotiable. Never one without the other. Affected tasks: Task 4 (author grass.crp.toml) and Task 9 (delete grass.crp).

When committing in the submodule, use file-scoped git commands (e.g. `git commit toml/6.surfacewater/grass.crp.toml -m "..."`) so pre-existing dirty state from prior work in other case directories is not bundled in.

---

## Pre-flight commands

Confirm the starting state:

```bash
cd /home/zawadzkim/Code/swap
git log --oneline -3
pixi run test-pfunit 2>&1 | tail -5
pixi run regression 2>&1 | tail -10
```

Expected: `Ok: 1, Fail: 0` for unit tests, `5 passed, 0 failed` for regression. If either is red, do not proceed — the baseline is broken.

---

## Task 1: Audit `readcropfixed`

Read `src/io/readswap.f90:2037-2517` end-to-end and produce a classification table that informs Tasks 2-7. Docs only.

**Files:**
- Create: `docs/phase-4f-readcropfixed-audit.md`

- [ ] **Step 1: Read the entire `readcropfixed` body**

Use the Read tool to load `src/io/readswap.f90` lines 2037-2517.

- [ ] **Step 2: Classify each non-blank, non-comment line into one of five buckets**

Bucket definitions (same as ADR 0015 / the swap.dra audit):

1. **READ** — `rd*` calls (`rdsdor`, `rdsinr`, `rdfdor`, `rdador`, `rdinqr`). Goes away with port.
2. **VALIDATE** — `if (...) call fatalerr(...)`. Push to validators in Task 2.
3. **NORMALIZE** — coordinate / unit conversion. Push to finalizers / adapter.
4. **RUNTIME** — non-read initialization (cumdens build, root density distribution, swrdc=0 default, etc.). Goes into `cropfixed_init` in Task 6.
5. **GUARDED** — branches our scope doesn't support (swdrought=2, swoxygen=2, swcompensate∈{1,2}, swcf=3, swharv=1, schedule=1, swinter∈{2,3}, swrd=2/3, swsalinity∈{1,2}, swgc=2 with overflow, swinco=3 .END-file path). For Phase 1, all should be stub-errored at validate time.

- [ ] **Step 3: Write the audit table**

Create `docs/phase-4f-readcropfixed-audit.md` with a brief intro paragraph + a Markdown table of the form:

```markdown
| Lines       | Bucket    | Notes                                              |
| ----------- | --------- | -------------------------------------------------- |
| 2075-2077   | READ      | Open .crp file, init parser                        |
| 2080-2087   | READ      | idev (+ tsumea/tsumam/tbase if idev=2)             |
| 2095-2105   | READ      | gctb table; SCF>1 fatal under swgc=2               |
| 2099-2104   | VALIDATE  | gctb cell > 1 under swgc=2                         |
| 2117-2147   | READ      | cftb / chtb / cfeictb tables based on swcf         |
| 2126-2135   | GUARDED   | swcf=3 wet-crop branch                             |
| 2168        | READ      | rdctb table                                        |
| 2171-2174   | READ      | swrd                                               |
| 2177-2179   | READ      | rdtb when swrd=1                                   |
| 2182-2191   | GUARDED   | swrd=2 + swdmi2rd                                  |
| 2195-2199   | VALIDATE+GUARDED | swrd=3 fatal in cropfixed                  |
| 2233-2238   | READ      | hlim1/hlim2u/hlim2l when swoxygen=1                |
| 2240-2278   | GUARDED   | swoxygen=2 entire Bartholomeus block               |
| 2298-2305   | READ      | hlim3h/hlim3l/hlim4/adcrh/adcrl when swdrought=1   |
| 2307-2320   | GUARDED   | swdrought=2 De Jong van Lier block                 |
| 2324-2347   | GUARDED   | swsalinity=1/=2 (entire flsolute block)            |
| 2386-2399   | GUARDED   | swcompensate ∈ {1, 2}                              |
| 2402-2427   | READ+GUARDED | swinter (=1 supported, =2/=3 stub-errored)      |
| 2430        | READ      | schedule                                           |
| 2432-2437   | GUARDED   | schedule=1 + swdrought=2 path                      |
| 2440-2446   | RUNTIME   | close file; irrigation(1) call when schedule=1     |
| 2449-2478   | RUNTIME   | cumdens computation (root density distribution)    |
| 2480-2514   | GUARDED   | swinco=3 .END-file path                            |
```

Cover all ~480 lines. Group consecutive same-bucket lines. Include a "Summary" section at the bottom with bucket counts.

- [ ] **Step 4: Commit**

```bash
git add docs/phase-4f-readcropfixed-audit.md
git commit -m "docs(audit): readcropfixed line-by-line classification for cropfixed port"
```

---

## Task 2: Extend `cropfixed_config_t` schema 1:1

Add the ~30 missing scalar fields and 6 table fields to `cropfixed_config_t`, plus stub-error validators for unsupported branch values. Per the spec's "schema 1:1" decision: every legacy field gets a TOML home, even when its parent switch is at a value we don't support.

**Files:**
- Modify: `src/config/cropfixed_config.f90`
- Test: `tests/unit/config/test_cropfixed_config.pf`

- [ ] **Step 1: Read the existing schema**

```bash
sed -n '1,100p' src/config/cropfixed_config.f90
```

- [ ] **Step 2: Write the failing tests**

Append the following block to `tests/unit/config/test_cropfixed_config.pf`. Each `@test` is a separate subroutine; check existing file style before pasting.

```fortran
! Phase 1 cropfixed port — stub-error validators for unsupported branches.

@test
subroutine test_cropfixed_swdrought_two_rejected()
   use funit
   use cropfixed_config_mod, only: cropfixed_config_t
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   type(cropfixed_config_t) :: c
   type(error_collection_t) :: errors
   integer :: i
   logical :: found
   c%swdrought = 2
   call c%validate(errors)
   found = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD .and. &
          index(errors%items(i)%message, 'swdrought=2') > 0) found = .true.
   end do
   @assertTrue(found)
end subroutine

@test
subroutine test_cropfixed_swoxygen_two_rejected()
   use funit
   use cropfixed_config_mod, only: cropfixed_config_t
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   type(cropfixed_config_t) :: c
   type(error_collection_t) :: errors
   integer :: i
   logical :: found
   c%swoxygen = 2
   call c%validate(errors)
   found = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD .and. &
          index(errors%items(i)%message, 'swoxygen=2') > 0) found = .true.
   end do
   @assertTrue(found)
end subroutine

@test
subroutine test_cropfixed_swcompensate_one_rejected()
   use funit
   use cropfixed_config_mod, only: cropfixed_config_t
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   type(cropfixed_config_t) :: c
   type(error_collection_t) :: errors
   integer :: i
   logical :: found
   c%swcompensate = 1
   call c%validate(errors)
   found = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD .and. &
          index(errors%items(i)%message, 'swcompensate') > 0) found = .true.
   end do
   @assertTrue(found)
end subroutine

@test
subroutine test_cropfixed_swcf_three_rejected()
   use funit
   use cropfixed_config_mod, only: cropfixed_config_t
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   type(cropfixed_config_t) :: c
   type(error_collection_t) :: errors
   integer :: i
   logical :: found
   c%swcf = 3
   call c%validate(errors)
   found = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD .and. &
          index(errors%items(i)%message, 'swcf=3') > 0) found = .true.
   end do
   @assertTrue(found)
end subroutine

@test
subroutine test_cropfixed_swharv_one_rejected()
   use funit
   use cropfixed_config_mod, only: cropfixed_config_t
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   type(cropfixed_config_t) :: c
   type(error_collection_t) :: errors
   integer :: i
   logical :: found
   c%swharv = 1
   call c%validate(errors)
   found = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD .and. &
          index(errors%items(i)%message, 'swharv=1') > 0) found = .true.
   end do
   @assertTrue(found)
end subroutine

@test
subroutine test_cropfixed_swinter_two_rejected()
   use funit
   use cropfixed_config_mod, only: cropfixed_config_t
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   type(cropfixed_config_t) :: c
   type(error_collection_t) :: errors
   integer :: i
   logical :: found
   c%swinter = 2
   call c%validate(errors)
   found = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD .and. &
          index(errors%items(i)%message, 'swinter=2') > 0) found = .true.
   end do
   @assertTrue(found)
end subroutine

@test
subroutine test_cropfixed_swrd_two_rejected()
   use funit
   use cropfixed_config_mod, only: cropfixed_config_t
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   type(cropfixed_config_t) :: c
   type(error_collection_t) :: errors
   integer :: i
   logical :: found
   c%swrd = 2
   call c%validate(errors)
   found = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD .and. &
          index(errors%items(i)%message, 'swrd') > 0) found = .true.
   end do
   @assertTrue(found)
end subroutine

@test
subroutine test_cropfixed_swsalinity_one_rejected()
   use funit
   use cropfixed_config_mod, only: cropfixed_config_t
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   type(cropfixed_config_t) :: c
   type(error_collection_t) :: errors
   integer :: i
   logical :: found
   c%swsalinity = 1
   call c%validate(errors)
   found = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD .and. &
          index(errors%items(i)%message, 'swsalinity') > 0) found = .true.
   end do
   @assertTrue(found)
end subroutine

@test
subroutine test_cropfixed_case6_supported_values_pass()
   use funit
   use cropfixed_config_mod, only: cropfixed_config_t
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   type(cropfixed_config_t) :: c
   type(error_collection_t) :: errors
   integer :: i
   logical :: any_stub
   ! Case 6 supported values for the new switches.
   c%swprep = 0
   c%swsow = 0
   c%swgerm = 0
   c%dvsend = 2.0d0
   c%swharv = 0
   c%idev = 1
   c%lcc  = 366
   c%kdif = 0.75d0
   c%kdir = 0.75d0
   c%swgc = 1
   c%swcf = 1
   c%swrd = 1
   c%swoxygen = 0
   c%swwrtnonox = 0
   c%swdrought = 1
   c%swsalinity = 0
   c%swcompensate = 0
   c%swinter = 1
   c%cofab = 0.25d0
   c%hlim3h = -200.0d0
   c%hlim3l = -800.0d0
   c%hlim4  = -8000.0d0
   c%adcrh = 0.5d0
   c%adcrl = 0.1d0
   call c%validate(errors)
   any_stub = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD) any_stub = .true.
   end do
   @assertFalse(any_stub)
end subroutine
```

- [ ] **Step 3: Run tests, verify failure**

```bash
pixi run test-pfunit 2>&1 | tail -15
```

Expected: at least 8 stub-error tests fail (the new validators don't exist yet). The "case 6 supported values" test may fail to build until the new fields are declared on the type.

- [ ] **Step 4: Extend the schema with new fields**

In `src/config/cropfixed_config.f90`, extend the `type :: cropfixed_config_t` block. Add the new scalar fields below the existing ones, grouped by section, before the `contains` line. Replace the existing type block (the `type :: cropfixed_config_t ... end type` section, lines 16-59) with this expanded version (preserving existing fields):

```fortran
   type :: cropfixed_config_t
      ! ====================================================================
      ! Existing fields (Phase 4c-a) — KEPT
      ! ====================================================================
      integer      :: idev = 1
      integer      :: lcc  = 0
      real(real64) :: kdif = 0.0_real64
      real(real64) :: kdir = 0.0_real64
      real(real64), allocatable :: cftb(:)
      real(real64), allocatable :: chtb(:)
      real(real64) :: rdi = 0.0_real64
      real(real64) :: rri = 0.0_real64
      real(real64) :: rdc = 0.0_real64
      real(real64), allocatable :: rdctb(:)
      real(real64) :: hlim1 = 0.0_real64
      real(real64) :: hlim2u = 0.0_real64
      real(real64) :: hlim2l = 0.0_real64
      real(real64) :: hlim3h = 0.0_real64
      real(real64) :: hlim3l = 0.0_real64
      real(real64) :: hlim4 = 0.0_real64
      real(real64) :: adcrh = 0.0_real64
      real(real64) :: adcrl = 0.0_real64
      real(real64) :: rsc   = 0.0_real64
      real(real64) :: ecmax  = 0.0_real64
      real(real64) :: ecslop = 0.0_real64
      real(real64) :: cofab = 0.0_real64
      type(irrigation_schedule_t) :: schedule

      ! ====================================================================
      ! Phase 1 (Phase 4f .crp port) additions
      ! ====================================================================

      ! Part 0a/b/c — preparation, sowing, germination
      integer :: swprep = 0    !! 0=no preparation, 1=preparation
      integer :: swsow  = 0    !! 0=no sowing, 1=sowing before crop growth
      integer :: swgerm = 0    !! 0=no germination, 1/2=temperature-based

      ! Part 0d — harvest
      real(real64) :: dvsend = 2.0_real64   !! Development stage at harvest
      integer      :: swharv = 0            !! 0=CROPEND-based, 1=DVS-based (stub-errored)

      ! Part 1 (idev=2 path; safe defaults for idev=1)
      real(real64) :: tsumea = 0.0_real64
      real(real64) :: tsumam = 0.0_real64
      real(real64) :: tbase  = 0.0_real64

      ! Part 3 — LAI vs SCF
      integer                   :: swgc = 1   !! 1=LAI, 2=SCF
      real(real64), allocatable :: gctb(:)    !! (dvs, lai|scf) flat pairs

      ! Part 4 — crop factor / height switch
      integer                   :: swcf = 1   !! 1=crop factor (cftb), 2=crop height (chtb), 3=wet-crop (stub-errored)
      real(real64), allocatable :: cfeictb(:) !! Wet-crop factor table (only if swcf=3 - unused in Phase 1)
      real(real64) :: albedo = 0.23_real64
      real(real64) :: rsw    = 0.0_real64

      ! Part 10 — root depth & density
      integer :: swrd     = 1   !! 1=DVS table (rdtb), 2=daily increase (stub-errored), 3=biomass (fatal-errored in legacy)
      integer :: swdmi2rd = 0   !! Only used when swrd=2 (stub-errored)
      integer :: swrdc    = 0   !! Switch development root density (legacy hard-codes to 0)
      real(real64), allocatable :: rdtb(:)    !! (dvs, rd) flat pairs (when swrd=1)

      ! Part 11 — oxygen stress
      integer :: swoxygen   = 0   !! 0=none, 1=Feddes, 2=Bartholomeus (stub-errored)
      integer :: swwrtnonox = 0   !! 0=no aerobic check, 1=check
      real(real64) :: aeratecrit = 1.0e-4_real64    !! Required when swwrtnonox=1

      ! Part 12 — drought stress
      integer :: swdrought = 1   !! 1=Feddes, 2=De Jong van Lier (stub-errored)

      ! Part 13 — salinity stress
      integer :: swsalinity = 0   !! 0=none, 1=Maas-Hoffman (stub-errored), 2=osmotic head (stub-errored)
      real(real64) :: saltmax   = 0.0_real64
      real(real64) :: saltslope = 0.0_real64
      real(real64) :: salthead  = 0.0_real64

      ! Part xx — root water uptake compensation
      integer :: swcompensate = 0   !! 0=none, 1=Jarvis (stub-errored), 2=Walsum (stub-errored)
      integer :: swstressor   = 1
      real(real64) :: alphacrit = 1.0_real64
      real(real64) :: dcritrtz  = 0.0_real64

      ! Part 14 — interception
      integer :: swinter = 1   !! 0=none, 1=Von Hoyningen-Hune (supported), 2=Gash (stub-errored), 3=storage-cap (stub-errored)

      ! Part 15 — irrigation scheduling top-level switch
      integer :: schedule_switch = 0   !! 0=no scheduling, 1=apply (stub-errored)

      ! Subordinate fields under stub-errored switches (schema 1:1; runtime never reads them)
      ! swdrought=2 fields (stub-errored at validate time):
      real(real64) :: wiltpoint  = 0.0_real64
      real(real64) :: kstem      = 0.0_real64
      real(real64) :: rxylem     = 0.0_real64
      real(real64) :: rootradius = 0.0_real64
      real(real64) :: kroot      = 0.0_real64
      real(real64) :: rootcoefa  = 0.0_real64
      real(real64) :: rooteff    = 0.0_real64
      real(real64) :: stephr     = 0.0_real64
      real(real64) :: criterhr   = 0.0_real64
      real(real64) :: taccur     = 0.0_real64
      ! swoxygen=2 fields (stub-errored at validate time):
      integer      :: swoxygentype          = 1
      integer      :: swrootradius          = 1
      integer      :: swtopsub              = 1
      integer      :: nrstaring             = 1
      real(real64) :: q10_root              = 0.0_real64
      real(real64) :: q10_microbial         = 0.0_real64
      real(real64) :: specific_resp_humus   = 0.0_real64
      real(real64) :: c_mroot               = 0.0_real64
      real(real64) :: srl                   = 0.0_real64
      real(real64) :: f_senes               = 0.0_real64
      real(real64) :: dry_mat_cont_roots    = 0.0_real64
      real(real64) :: air_filled_root_por   = 0.0_real64
      real(real64) :: spec_weight_root_tissue = 0.0_real64
      real(real64) :: var_a                 = 0.0_real64
      real(real64) :: root_radiusO2         = 1.0e-3_real64
   contains
      procedure :: validate => cropfixed_config_validate
      procedure :: finalize => cropfixed_config_finalize
   end type cropfixed_config_t
```

- [ ] **Step 5: Replace `cropfixed_config_validate` to add stub-error guards**

In the same file, replace the existing `cropfixed_config_validate` subroutine with this expanded version. Keep all existing range checks; add the new stub-error block at the top (so unsupported branches are rejected before any other validation runs).

```fortran
   subroutine cropfixed_config_validate(self, errors)
      class(cropfixed_config_t), intent(in)    :: self
      type(error_collection_t),  intent(inout) :: errors

      ! ----- Phase 1 stub-errors for unsupported runtime branches -----
      ! ADR 0015: schema accepts the value 1:1 with legacy; runtime
      ! plumbing for these branches has not been ported. Cases that
      ! need them must run via the legacy executable.
      if (self%swdrought == 2) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropfixed.swdrought=2 (De Jong van Lier) not yet supported in ' // &
            'the TOML pipeline; use the legacy executable.', 'cropfixed')
      end if
      if (self%swoxygen == 2) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropfixed.swoxygen=2 (Bartholomeus) not yet supported in ' // &
            'the TOML pipeline; use the legacy executable.', 'cropfixed')
      end if
      if (self%swcompensate /= 0) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropfixed.swcompensate /= 0 (Jarvis/Walsum compensation) not ' // &
            'yet supported in the TOML pipeline; use the legacy executable.', &
            'cropfixed')
      end if
      if (self%swcf == 3) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropfixed.swcf=3 (wet-crop factor) not yet supported in the ' // &
            'TOML pipeline; use the legacy executable.', 'cropfixed')
      end if
      if (self%swharv == 1) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropfixed.swharv=1 (DVS-based harvest timing) not yet supported ' // &
            'in the TOML pipeline; use the legacy executable.', 'cropfixed')
      end if
      if (self%swinter == 2 .or. self%swinter == 3) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropfixed.swinter=2 or 3 (Gash / storage-cap interception) ' // &
            'not yet supported in the TOML pipeline; use the legacy ' // &
            'executable.', 'cropfixed')
      end if
      if (self%swrd == 2 .or. self%swrd == 3) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropfixed.swrd != 1 (alternate root extension methods) not ' // &
            'yet supported in the TOML pipeline; use the legacy executable.', &
            'cropfixed')
      end if
      if (self%swsalinity == 1 .or. self%swsalinity == 2) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropfixed.swsalinity != 0 (Maas-Hoffman / osmotic head) not ' // &
            'yet supported in the TOML pipeline; use the legacy executable.', &
            'cropfixed')
      end if
      if (self%schedule_switch == 1) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropfixed.schedule=1 (per-crop irrigation scheduling) not yet ' // &
            'supported in the TOML pipeline; use the legacy executable.', &
            'cropfixed')
      end if

      ! ----- Top-level enum / range checks -----
      call check_int_enum(self%idev, [1, 2], 'cropfixed.idev', errors)
      if (self%idev == 1) then
         call check_int_range(self%lcc, 1, 366, 'cropfixed.lcc', errors)
      end if

      call check_int_enum(self%swgc,    [1, 2], 'cropfixed.swgc',    errors)
      call check_int_enum(self%swcf,    [1, 2, 3], 'cropfixed.swcf', errors)
      call check_int_enum(self%swrd,    [1, 2, 3], 'cropfixed.swrd', errors)
      call check_int_enum(self%swoxygen,    [0, 1, 2], 'cropfixed.swoxygen',    errors)
      call check_int_enum(self%swwrtnonox,  [0, 1],    'cropfixed.swwrtnonox',  errors)
      call check_int_enum(self%swdrought,   [1, 2],    'cropfixed.swdrought',   errors)
      call check_int_enum(self%swsalinity,  [0, 1, 2], 'cropfixed.swsalinity',  errors)
      call check_int_enum(self%swcompensate,[0, 1, 2], 'cropfixed.swcompensate',errors)
      call check_int_enum(self%swinter,     [0, 1, 2, 3], 'cropfixed.swinter',  errors)
      call check_int_enum(self%swharv,      [0, 1],    'cropfixed.swharv',      errors)
      call check_int_enum(self%swprep,      [0, 1],    'cropfixed.swprep',      errors)
      call check_int_enum(self%swsow,       [0, 1],    'cropfixed.swsow',       errors)
      call check_int_enum(self%swgerm,      [0, 1, 2], 'cropfixed.swgerm',      errors)
      call check_int_enum(self%schedule_switch, [0, 1], 'cropfixed.schedule', errors)

      call check_real_range(self%dvsend, 0.0_real64, 3.0_real64,    'cropfixed.dvsend', errors)
      call check_real_range(self%kdif, 0.0_real64,  2.0_real64,     'cropfixed.kdif',   errors)
      call check_real_range(self%kdir, 0.0_real64,  2.0_real64,     'cropfixed.kdir',   errors)
      call check_real_range(self%rdi,  0.0_real64, 1000.0_real64,   'cropfixed.rdi',    errors)
      call check_real_range(self%rdc,  0.0_real64, 1000.0_real64,   'cropfixed.rdc',    errors)
      call check_nonnegative_real(self%rri, 'cropfixed.rri', errors)

      ! Feddes ordering — applied only when the corresponding branch is active.
      if (self%swdrought == 1) then
         call check_ordered_pair(self%hlim3l, self%hlim3h, 'hlim3l', 'hlim3h', 'cropfixed', errors)
      end if

      call check_nonnegative_real(self%ecmax,  'cropfixed.ecmax',  errors)
      call check_nonnegative_real(self%ecslop, 'cropfixed.ecslop', errors)

      call self%schedule%validate(errors)
   end subroutine cropfixed_config_validate
```

- [ ] **Step 6: Run tests, verify pass**

```bash
pixi run test-pfunit 2>&1 | tail -10
```

Expected: `Ok: 1, Fail: 0`. If a test fails because the message substring doesn't exactly match, adjust the message or the substring (the test searches for things like `'swdrought=2'` — make sure the validator emits that literal substring).

- [ ] **Step 7: Commit**

```bash
git add src/config/cropfixed_config.f90 tests/unit/config/test_cropfixed_config.pf
git commit -m "feat(config): cropfixed schema 1:1 with legacy + stub-error validators

Phase 1 of the .crp port (cropfixed via case 6). Schema gains all the
fields legacy readcropfixed reads, including subordinate fields under
switches that are stub-errored. Per ADR 0015, runtime branches not yet
ported (swdrought=2, swoxygen=2, swcompensate /= 0, swcf=3, swharv=1,
swinter ∈ {2,3}, swrd /= 1, swsalinity /= 0, schedule=1) reject at
validate time with ERR_VALIDATION_CROSS_FIELD."
```

---

## Task 3: Extend `read_cropfixed_toml.f90` parser

Read all the new scalars and table fields from the .crp.toml. Verifies via round-trip of a fixture.

**Files:**
- Modify: `src/io/toml/read_cropfixed_toml.f90`
- Test: `tests/unit/io/toml/test_read_cropfixed_toml.pf`
- Fixture: `tests/unit/io/toml/fixtures/cropfixed_full.toml` (create)

- [ ] **Step 1: Write the failing round-trip test**

Append to `tests/unit/io/toml/test_read_cropfixed_toml.pf`:

```fortran
@test
subroutine test_cropfixed_toml_full_roundtrip()
   use funit
   use error_mod
   use cropfixed_config_mod, only: cropfixed_config_t
   use read_cropfixed_toml_mod, only: read_cropfixed_toml
   use tomlf, only: toml_table, toml_load, toml_error
   type(cropfixed_config_t)       :: c
   type(error_collection_t)       :: errors
   type(toml_table), allocatable, target :: doc
   type(toml_table), pointer      :: doc_ptr
   type(toml_error), allocatable  :: terr
   call toml_load(doc, &
        'tests/unit/io/toml/fixtures/cropfixed_full.toml', error=terr)
   @assertFalse(allocated(terr))
   doc_ptr => doc
   call read_cropfixed_toml(doc_ptr, c, errors)
   @assertFalse(errors%has_errors())
   ! Top-level switches
   @assertEqual(0, c%swprep)
   @assertEqual(0, c%swsow)
   @assertEqual(0, c%swgerm)
   @assertEqual(2.0d0, c%dvsend, 1.0d-12)
   @assertEqual(0, c%swharv)
   @assertEqual(1, c%idev)
   @assertEqual(366, c%lcc)
   @assertEqual(0.75d0, c%kdif, 1.0d-12)
   @assertEqual(0.75d0, c%kdir, 1.0d-12)
   @assertEqual(1, c%swgc)
   @assertEqual(1, c%swcf)
   @assertEqual(1, c%swrd)
   @assertEqual(0, c%swoxygen)
   @assertEqual(0, c%swwrtnonox)
   @assertEqual(1, c%swdrought)
   @assertEqual(0, c%swsalinity)
   @assertEqual(0, c%swcompensate)
   @assertEqual(1, c%swinter)
   @assertEqual(0.25d0, c%cofab, 1.0d-12)
   @assertEqual(0, c%schedule_switch)
   ! Tables (each is a flat array of (dvs, value) pairs).
   @assertTrue(allocated(c%gctb))
   @assertEqual(4, size(c%gctb))
   @assertEqual(0.0d0, c%gctb(1), 1.0d-12)
   @assertEqual(3.0d0, c%gctb(2), 1.0d-12)
   @assertEqual(2.0d0, c%gctb(3), 1.0d-12)
   @assertEqual(3.0d0, c%gctb(4), 1.0d-12)
   @assertTrue(allocated(c%cftb))
   @assertEqual(4, size(c%cftb))
   @assertEqual(1.0d0, c%cftb(2), 1.0d-12)
   @assertTrue(allocated(c%rdtb))
   @assertEqual(4, size(c%rdtb))
   @assertEqual(30.0d0, c%rdtb(2), 1.0d-12)
   @assertTrue(allocated(c%rdctb))
   @assertEqual(4, size(c%rdctb))
   @assertEqual(1.0d0, c%rdctb(2), 1.0d-12)
   ! Feddes — drought (swdrought=1)
   @assertEqual(-200.0d0,  c%hlim3h, 1.0d-12)
   @assertEqual(-800.0d0,  c%hlim3l, 1.0d-12)
   @assertEqual(-8000.0d0, c%hlim4,  1.0d-12)
   @assertEqual(0.5d0,     c%adcrh,  1.0d-12)
   @assertEqual(0.1d0,     c%adcrl,  1.0d-12)
end subroutine
```

- [ ] **Step 2: Create the fixture**

Create `tests/unit/io/toml/fixtures/cropfixed_full.toml`:

```toml
# Reproduces case 6 grass.crp at the schema level. All switches at supported values.

[preparation]
swprep = 0
swsow  = 0
swgerm = 0

[harvest]
dvsend = 2.0
swharv = 0

[phenology]
idev = 1
lcc  = 366

[light]
kdif = 0.75
kdir = 0.75

[lai]
swgc = 1
gctb = [0.0, 3.0, 2.0, 3.0]

[crop_factor]
swcf = 1
cftb = [0.0, 1.0, 2.0, 1.0]

[root]
swrd = 1
swdmi2rd = 0
swrdc = 0
rdtb  = [0.0, 30.0, 2.0, 30.0]
rdctb = [0.0, 1.0, 1.0, 1.0]
rdi = 30.0
rri = 0.0
rdc = 30.0

[oxygen_stress]
swoxygen   = 0
swwrtnonox = 0

[drought_stress]
swdrought = 1
hlim3h    = -200.0
hlim3l    = -800.0
hlim4     = -8000.0
adcrh     = 0.5
adcrl     = 0.1

[salinity_stress]
swsalinity = 0

[compensation]
swcompensate = 0

[interception]
swinter = 1
cofab   = 0.25

[scheduling]
schedule = 0
```

- [ ] **Step 3: Run, verify failure**

```bash
pixi run test-pfunit 2>&1 | tail -10
```

Expected: round-trip test fails because the parser doesn't yet read these sections.

- [ ] **Step 4: Extend the parser**

In `src/io/toml/read_cropfixed_toml.f90`, replace the current body of `read_cropfixed_toml` with this expanded version. The TOML section names follow the new fixture (`preparation`, `harvest`, `phenology`, `light`, `lai`, `crop_factor`, `root`, `oxygen_stress`, `drought_stress`, `salinity_stress`, `compensation`, `interception`, `scheduling`).

A new helper for parsing a flat `(dvs, value)` array of pairs is needed. Add a private subroutine `read_pair_array` to the same module; it takes a TOML table, a key name, an output allocatable, a label, and the errors collection. It uses `tomlf` array iteration to fill the output array (must be even-sized).

Replacement subroutine + new helper:

```fortran
   subroutine read_pair_array(tbl, key, arr, label, errors)
      use tomlf, only: toml_array, toml_value, get_value, len
      use error_mod, only: error_collection_t, ERR_PARSE_TYPE_MISMATCH, ERR_VALIDATION_OUT_OF_RANGE
      type(toml_table), pointer, intent(in)    :: tbl
      character(len=*),          intent(in)    :: key
      real(real64), allocatable, intent(out)   :: arr(:)
      character(len=*),          intent(in)    :: label
      type(error_collection_t),  intent(inout) :: errors
      type(toml_array), pointer :: a
      integer :: stat, i, n
      real(real64) :: v
      if (.not. associated(tbl)) return
      call get_value(tbl, key, a, requested=.false., stat=stat)
      if (stat /= 0 .or. .not. associated(a)) return
      n = len(a)
      if (mod(n, 2) /= 0) then
         call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
            'expected even-length (dvs,value) pair array', label)
         return
      end if
      allocate(arr(n))
      do i = 1, n
         call get_value(a, i, v, stat=stat)
         if (stat /= 0) then
            call errors%append(ERR_PARSE_TYPE_MISMATCH, &
               'non-real cell in pair array', label)
            deallocate(arr)
            return
         end if
         arr(i) = v
      end do
   end subroutine read_pair_array

   subroutine read_cropfixed_toml(doc, config, errors)
      type(toml_table), pointer,  intent(in)    :: doc
      type(cropfixed_config_t),   intent(inout) :: config
      type(error_collection_t),   intent(inout) :: errors

      type(toml_table), pointer :: prep, harv, ph, light, lai, cf, root, &
                                    ox, dr, sa, cp, inter, sched, irr_sched

      ! Preparation, sowing, germination
      call get_table(doc, 'preparation', prep, 'preparation', errors)
      if (associated(prep)) then
         call get_optional_int_with_default(prep, 'swprep', config%swprep, 0, 'preparation.swprep', errors)
         call get_optional_int_with_default(prep, 'swsow',  config%swsow,  0, 'preparation.swsow',  errors)
         call get_optional_int_with_default(prep, 'swgerm', config%swgerm, 0, 'preparation.swgerm', errors)
      end if

      ! Harvest
      call get_table(doc, 'harvest', harv, 'harvest', errors)
      if (associated(harv)) then
         call get_optional_real_with_default(harv, 'dvsend', config%dvsend, 2.0_real64, 'harvest.dvsend', errors)
         call get_optional_int_with_default (harv, 'swharv', config%swharv, 0,            'harvest.swharv', errors)
      end if

      ! Phenology
      call get_table(doc, 'phenology', ph, 'phenology', errors)
      if (associated(ph)) then
         call get_optional_int_with_default (ph, 'idev',   config%idev,   1, 'phenology.idev',   errors)
         call get_optional_int_with_default (ph, 'lcc',    config%lcc,    0, 'phenology.lcc',    errors)
         call get_optional_real_with_default(ph, 'tsumea', config%tsumea, 0.0_real64, 'phenology.tsumea', errors)
         call get_optional_real_with_default(ph, 'tsumam', config%tsumam, 0.0_real64, 'phenology.tsumam', errors)
         call get_optional_real_with_default(ph, 'tbase',  config%tbase,  0.0_real64, 'phenology.tbase',  errors)
      end if

      ! Light
      call get_table(doc, 'light', light, 'light', errors)
      if (associated(light)) then
         call get_optional_real_with_default(light, 'kdif', config%kdif, 0.0_real64, 'light.kdif', errors)
         call get_optional_real_with_default(light, 'kdir', config%kdir, 0.0_real64, 'light.kdir', errors)
      end if

      ! LAI / SCF table
      call get_table(doc, 'lai', lai, 'lai', errors)
      if (associated(lai)) then
         call get_optional_int_with_default(lai, 'swgc', config%swgc, 1, 'lai.swgc', errors)
         call read_pair_array(lai, 'gctb', config%gctb, 'lai.gctb', errors)
      end if

      ! Crop factor / height
      call get_table(doc, 'crop_factor', cf, 'crop_factor', errors)
      if (associated(cf)) then
         call get_optional_int_with_default(cf, 'swcf', config%swcf, 1, 'crop_factor.swcf', errors)
         call read_pair_array(cf, 'cftb', config%cftb, 'crop_factor.cftb', errors)
         call read_pair_array(cf, 'chtb', config%chtb, 'crop_factor.chtb', errors)
         call get_optional_real_with_default(cf, 'albedo', config%albedo, 0.23_real64, 'crop_factor.albedo', errors)
         call get_optional_real_with_default(cf, 'rsc',    config%rsc,    0.0_real64,  'crop_factor.rsc',    errors)
         call get_optional_real_with_default(cf, 'rsw',    config%rsw,    0.0_real64,  'crop_factor.rsw',    errors)
      end if

      ! Root
      call get_table(doc, 'root', root, 'root', errors)
      if (associated(root)) then
         call get_optional_int_with_default (root, 'swrd',     config%swrd,     1,           'root.swrd',     errors)
         call get_optional_int_with_default (root, 'swdmi2rd', config%swdmi2rd, 0,           'root.swdmi2rd', errors)
         call get_optional_int_with_default (root, 'swrdc',    config%swrdc,    0,           'root.swrdc',    errors)
         call read_pair_array(root, 'rdtb',  config%rdtb,  'root.rdtb',  errors)
         call read_pair_array(root, 'rdctb', config%rdctb, 'root.rdctb', errors)
         call get_optional_real_with_default(root, 'rdi', config%rdi, 0.0_real64, 'root.rdi', errors)
         call get_optional_real_with_default(root, 'rri', config%rri, 0.0_real64, 'root.rri', errors)
         call get_optional_real_with_default(root, 'rdc', config%rdc, 0.0_real64, 'root.rdc', errors)
      end if

      ! Oxygen stress
      call get_table(doc, 'oxygen_stress', ox, 'oxygen_stress', errors)
      if (associated(ox)) then
         call get_optional_int_with_default (ox, 'swoxygen',   config%swoxygen,   0, 'oxygen_stress.swoxygen',   errors)
         call get_optional_int_with_default (ox, 'swwrtnonox', config%swwrtnonox, 0, 'oxygen_stress.swwrtnonox', errors)
         call get_optional_real_with_default(ox, 'aeratecrit', config%aeratecrit, 1.0e-4_real64, 'oxygen_stress.aeratecrit', errors)
         call get_optional_real_with_default(ox, 'hlim1',      config%hlim1,  0.0_real64, 'oxygen_stress.hlim1',  errors)
         call get_optional_real_with_default(ox, 'hlim2u',     config%hlim2u, 0.0_real64, 'oxygen_stress.hlim2u', errors)
         call get_optional_real_with_default(ox, 'hlim2l',     config%hlim2l, 0.0_real64, 'oxygen_stress.hlim2l', errors)
      end if

      ! Drought stress
      call get_table(doc, 'drought_stress', dr, 'drought_stress', errors)
      if (associated(dr)) then
         call get_optional_int_with_default (dr, 'swdrought', config%swdrought, 1,           'drought_stress.swdrought', errors)
         call get_optional_real_with_default(dr, 'hlim3h',    config%hlim3h, 0.0_real64, 'drought_stress.hlim3h', errors)
         call get_optional_real_with_default(dr, 'hlim3l',    config%hlim3l, 0.0_real64, 'drought_stress.hlim3l', errors)
         call get_optional_real_with_default(dr, 'hlim4',     config%hlim4,  0.0_real64, 'drought_stress.hlim4',  errors)
         call get_optional_real_with_default(dr, 'adcrh',     config%adcrh,  0.0_real64, 'drought_stress.adcrh',  errors)
         call get_optional_real_with_default(dr, 'adcrl',     config%adcrl,  0.0_real64, 'drought_stress.adcrl',  errors)
      end if

      ! Salinity stress (subordinate fields read but stub-errored at validate)
      call get_table(doc, 'salinity_stress', sa, 'salinity_stress', errors)
      if (associated(sa)) then
         call get_optional_int_with_default (sa, 'swsalinity', config%swsalinity, 0,           'salinity_stress.swsalinity', errors)
         call get_optional_real_with_default(sa, 'saltmax',    config%saltmax,    0.0_real64,  'salinity_stress.saltmax',    errors)
         call get_optional_real_with_default(sa, 'saltslope',  config%saltslope,  0.0_real64,  'salinity_stress.saltslope',  errors)
         call get_optional_real_with_default(sa, 'salthead',   config%salthead,   0.0_real64,  'salinity_stress.salthead',   errors)
         call get_optional_real_with_default(sa, 'ecmax',      config%ecmax,      0.0_real64,  'salinity_stress.ecmax',      errors)
         call get_optional_real_with_default(sa, 'ecslop',     config%ecslop,     0.0_real64,  'salinity_stress.ecslop',     errors)
      end if

      ! Compensation (stub-errored at validate)
      call get_table(doc, 'compensation', cp, 'compensation', errors)
      if (associated(cp)) then
         call get_optional_int_with_default (cp, 'swcompensate', config%swcompensate, 0, 'compensation.swcompensate', errors)
         call get_optional_int_with_default (cp, 'swstressor',   config%swstressor,   1, 'compensation.swstressor',   errors)
         call get_optional_real_with_default(cp, 'alphacrit',    config%alphacrit,    1.0_real64,  'compensation.alphacrit', errors)
         call get_optional_real_with_default(cp, 'dcritrtz',     config%dcritrtz,     0.0_real64,  'compensation.dcritrtz',  errors)
      end if

      ! Interception
      call get_table(doc, 'interception', inter, 'interception', errors)
      if (associated(inter)) then
         call get_optional_int_with_default (inter, 'swinter', config%swinter, 1,           'interception.swinter', errors)
         call get_optional_real_with_default(inter, 'cofab',   config%cofab,   0.0_real64,  'interception.cofab',   errors)
      end if

      ! Scheduling top-level switch
      call get_table(doc, 'scheduling', sched, 'scheduling', errors)
      if (associated(sched)) then
         call get_optional_int_with_default(sched, 'schedule', config%schedule_switch, 0, 'scheduling.schedule', errors)
      end if

      ! Per-crop irrigation schedule (existing field - kept)
      call get_table(doc, 'irrigation_schedule', irr_sched, 'irrigation_schedule', errors)
      call read_irrigation_schedule_from_section(irr_sched, config%schedule, errors)
   end subroutine read_cropfixed_toml
```

- [ ] **Step 5: Run, verify pass**

```bash
pixi run test-pfunit 2>&1 | tail -10
```

Expected: `Ok: 1, Fail: 0` including the new round-trip test.

- [ ] **Step 6: Commit**

```bash
git add src/io/toml/read_cropfixed_toml.f90 tests/unit/io/toml/test_read_cropfixed_toml.pf tests/unit/io/toml/fixtures/cropfixed_full.toml
git commit -m "feat(io/toml): cropfixed parser reads schema 1:1 with legacy

Extends read_cropfixed_toml to populate every field that legacy
readcropfixed reads (~30 new scalars + 6 flat-pair tables: gctb, cftb,
chtb, rdtb, rdctb, cfeictb). Fixture covers all fields at case-6's
supported values."
```

---

## Task 4: Author case 6's `grass.crp.toml` 1:1

Translate legacy `tests/swap-cases/6.surfacewater/grass.crp` (202 lines) into TOML form. This is a SUBMODULE PAIR COMMIT.

**Files:**
- Modify (in submodule): `tests/swap-cases/toml/6.surfacewater/grass.crp.toml`
- Test: `tests/unit/io/toml/test_load_swap_config.pf`

- [ ] **Step 1: Write the failing assertion**

Append to `tests/unit/io/toml/test_load_swap_config.pf`:

```fortran
@test
subroutine test_load_case6_grass_crp_toml_populated()
   use funit
   use error_mod, only: error_collection_t
   use swap_config_mod, only: swap_config_t
   use load_swap_config_mod, only: load_swap_config
   type(swap_config_t)      :: c
   type(error_collection_t) :: errors
   integer :: i, n_loaded
   call load_swap_config('tests/swap-cases/toml/6.surfacewater/swap.toml', c, errors)
   @assertFalse(errors%has_errors())
   ! Case 6 has 3 type-1 rotations all referencing grass.crp.toml.
   @assertEqual(3, size(c%crop%rotation_loaded))
   n_loaded = 0
   do i = 1, 3
      if (c%crop%rotation_loaded(i)) n_loaded = n_loaded + 1
   end do
   @assertEqual(3, n_loaded)
   ! Verify all three slots got the same parsed content.
   @assertEqual(1, c%crop%rotation_fixed(1)%idev)
   @assertEqual(366, c%crop%rotation_fixed(1)%lcc)
   @assertEqual(0.75d0, c%crop%rotation_fixed(2)%kdif, 1.0d-12)
   @assertEqual(0.25d0, c%crop%rotation_fixed(3)%cofab, 1.0d-12)
   @assertTrue(allocated(c%crop%rotation_fixed(1)%gctb))
   @assertEqual(4, size(c%crop%rotation_fixed(1)%gctb))
   @assertEqual(3.0d0, c%crop%rotation_fixed(1)%gctb(2), 1.0d-12)
   @assertTrue(allocated(c%crop%rotation_fixed(1)%rdctb))
   @assertEqual(4, size(c%crop%rotation_fixed(1)%rdctb))
end subroutine
```

- [ ] **Step 2: Run, verify failure**

```bash
pixi run test-pfunit 2>&1 | tail -10
```

Expected: assertions fail because the existing skeletal `grass.crp.toml` doesn't have the full content.

- [ ] **Step 3: Replace case 6's grass.crp.toml**

Overwrite `tests/swap-cases/toml/6.surfacewater/grass.crp.toml` with this 1:1 reproduction of legacy `grass.crp`:

```toml
# Case 6 surfacewater grass crop (type 1, fixed). 1:1 reproduction of
# tests/swap-cases/6.surfacewater/grass.crp. All switches at supported values.

[preparation]
swprep = 0
swsow  = 0
swgerm = 0

[harvest]
dvsend = 2.0
swharv = 0

[phenology]
idev = 1
lcc  = 366

[light]
kdif = 0.75
kdir = 0.75

[lai]
swgc = 1
gctb = [0.0, 3.0, 2.0, 3.0]

[crop_factor]
swcf = 1
cftb = [0.0, 1.0, 2.0, 1.0]

[root]
swrd = 1
swdmi2rd = 0
swrdc = 0
rdtb  = [0.0, 30.0, 2.0, 30.0]
rdctb = [0.0, 1.0, 1.0, 1.0]
rdi = 30.0
rri = 0.0
rdc = 30.0

[oxygen_stress]
swoxygen   = 0
swwrtnonox = 0

[drought_stress]
swdrought = 1
hlim3h    = -200.0
hlim3l    = -800.0
hlim4     = -8000.0
adcrh     = 0.5
adcrl     = 0.1

[salinity_stress]
swsalinity = 0

[compensation]
swcompensate = 0

[interception]
swinter = 1
cofab   = 0.25

[scheduling]
schedule = 0
```

- [ ] **Step 4: Run, verify pass**

```bash
pixi run test-pfunit 2>&1 | tail -10
```

Expected: `Ok: 1, Fail: 0` including the new case-6 assertion.

- [ ] **Step 5: Submodule inner-commit**

Use file-scoped `git commit` to avoid bundling pre-existing dirty state from other case dirs.

```bash
cd tests/swap-cases
git commit toml/6.surfacewater/grass.crp.toml -m "feat(toml/6.surfacewater): grass.crp.toml 1:1 with legacy grass.crp

Replaces the skeletal 33-line .crp.toml with full 1:1 schema coverage
of legacy grass.crp. All switches at case-6's supported values:
swprep=0, swsow=0, swgerm=0, swharv=0, idev=1, swgc=1, swcf=1, swrd=1,
swoxygen=0, swwrtnonox=0, swdrought=1, swsalinity=0, swcompensate=0,
swinter=1, schedule=0."
cd ../..
```

- [ ] **Step 6: Outer-repo bump-and-test commit**

```bash
git add tests/swap-cases tests/unit/io/toml/test_load_swap_config.pf
git commit -m "test(toml): assert case-6 grass.crp.toml fully populated; bump submodule"
```

---

## Task 5: Module-level `crop_config_global` pointer

Add a module-level pointer that the runtime path can use to read `config%crop` without changing the dispatch signature of every legacy crop sub. **Transitional** — see ADR 0016 for teardown.

**Files:**
- Create: `src/crop/crop_config_global.f90`
- Modify: `src/io/toml/config_to_variables.f90` (set the pointer at the end of the crop block)
- Modify: `meson.build` (register the new module)

- [ ] **Step 1: Create the module**

Create `src/crop/crop_config_global.f90`:

```fortran
!> Module-level pointer to the parsed crop config.
!!
!! Transitional element of the .crp port (ADR 0016). The runtime path
!! reads through this pointer to access the per-rotation cache stored
!! in `crop_config_t.rotation_fixed(:)` etc., without changing the
!! signatures of every legacy crop subroutine. After Phase 4 of the
!! .crp port and the config-passing follow-on spec, computation subs
!! take typed config + state as explicit arguments and this module
!! goes away entirely.
!!
!! Lifecycle: set by config_to_variables at the end of its crop block;
!! valid for the duration of the simulation (the underlying config
!! lives in a local of core/swap.f90's iTask=1 block, which encloses
!! all simulation initialization). Readers must check `associated(...)`
!! defensively.
module crop_config_global_mod
   use crop_config_mod, only: crop_config_t
   implicit none
   private

   public :: crop_config_global

   type(crop_config_t), pointer :: crop_config_global => null()
end module crop_config_global_mod
```

- [ ] **Step 2: Register in `meson.build`**

In the top-level `meson.build`, find the `sources` list and add this line in the alphabetical-by-path region (around the existing `src/config/crop_config.f90` entry, but in `src/crop/` ordering — likely before `src/crop/cropgrowth.f90`):

```
    'src/crop/crop_config_global.f90',
```

Order matters: `crop_config_global.f90` must compile after `src/config/crop_config.f90` (which defines `crop_config_t`) but before `src/crop/cropgrowth.f90` and `src/io/toml/config_to_variables.f90` (which use it).

- [ ] **Step 3: Set the pointer in `config_to_variables`**

In `src/io/toml/config_to_variables.f90`, find the crop block (around line 1141, where `swCrop = config%crop%swcrop` lives). At the END of the crop block (after all the rotation-array writes; look for where `config%crop` references stop appearing), add:

```fortran
      ! Phase 1 (.crp port): expose the parsed crop config to runtime
      ! subs that need per-rotation cache access. Transitional — see
      ! ADR 0016. The pointer targets the caller's local config; valid
      ! for the duration of the simulation init.
      block
         use crop_config_global_mod, only: crop_config_global
         crop_config_global => config%crop
      end block
```

Note: the pointer needs the `target` attribute on the source. Since `config_to_variables` takes `config` as `intent(in) :: config` (or similar), Fortran allows assigning `config%crop` to a pointer only if `config` is declared with `target` in the caller. Check the caller (`core/swap.f90:140`): `type(swap_config_t) :: config`. Add `target` there:

```fortran
      type(swap_config_t), target  :: config
```

- [ ] **Step 4: Build, verify clean**

```bash
pixi run swap 2>&1 | tail -20
```

Expected: clean build (no module-not-found, no syntax errors). If the `target` attribute change in `core/swap.f90` produces a complaint, adjust to where the `target` attribute is required.

- [ ] **Step 5: Run unit + regression to confirm no behavioral change yet**

```bash
pixi run test-pfunit 2>&1 | tail -5
pixi run regression 2>&1 | tail -5
```

Expected: tests still green. Adding the pointer alone changes no runtime behavior — it just exposes data that nothing yet reads.

- [ ] **Step 6: Commit**

```bash
git add src/crop/crop_config_global.f90 src/io/toml/config_to_variables.f90 src/core/swap.f90 meson.build
git commit -m "feat(crop): crop_config_global module pointer (transitional)

ADR 0016: bridge from typed config (held in core/swap.f90 local scope)
to legacy runtime subs that don't yet take config as an argument.
Teardown trigger is the config-passing follow-on spec (post-Phase 4
of the .crp port). Readers check associated(...) defensively."
```

---

## Task 6: New `cropfixed_init` runtime module

Create the runtime sub that the cropgrowth dispatcher calls in lieu of legacy `readcropfixed`. Two halves: copy config-to-globals, then runtime init math.

**Files:**
- Create: `src/crop/cropfixed_init.f90`
- Modify: `meson.build` (register)
- Test: `tests/unit/crop/test_cropfixed_init.pf` (new directory)
- Modify: `tests/unit/meson.build` (register pf file + suite)
- Modify: `tests/unit/testSuites.inc` (register suite)

- [ ] **Step 1: Create the module**

Create `src/crop/cropfixed_init.f90`:

```fortran
!> Runtime initialization for the type-1 (cropfixed) rotation, replacing
!! the per-rotation runtime side-effects of legacy `readcropfixed`.
!!
!! Operates on `cropfixed_config_t` (already validated and populated at
!! config-load time) and writes to the same `variables`-module globals
!! that legacy `readcropfixed` writes to. After the copy, runs the
!! tail-of-readcropfixed init math (cumdens build, swrdc default).
!!
!! Scope (Phase 1, ADR 0015): only the supported subset of switches.
!! Defense-in-depth runtime guards mirror the validator stub-errors.
module cropfixed_init_mod
   use iso_fortran_env, only: real64
   use cropfixed_config_mod, only: cropfixed_config_t
   implicit none
   private

   public :: cropfixed_init_from_config

contains

   subroutine cropfixed_init_from_config(cfg, icrop)
      use variables, only: idev, lcc, tsumea, tsumam, tbase,            &
                            kdif, kdir, gctb, swgc,                       &
                            cftb, chtb, cfeictb, swcf, albedo, rsc, rsw,  &
                            rdtb, rdctb, swrd, swdmi2rd, swrdc, rdi, rri, rdc, &
                            swoxygen, swwrtnonox, aeratecrit,             &
                            hlim1, hlim2u, hlim2l,                        &
                            swdrought, hlim3h, hlim3l, hlim4, adcrh, adcrl, &
                            swsalinity, saltmax, saltslope, salthead,     &
                            swcompensate, swstressor, alphacrit, dcritrtz, &
                            swinter, cofab,                               &
                            schedule, dvsend, swharv,                     &
                            cumdens
      use array_utils,  only: afgen
      use error_mod,    only: fatalerr_collected
      type(cropfixed_config_t), intent(in) :: cfg
      integer,                  intent(in) :: icrop

      integer      :: i
      real(real64) :: depth, rootdis(202), sum

      ! ---- Defense-in-depth: stub-error gates mirror cropfixed_config_validate.
      if (cfg%swdrought == 2 .or. cfg%swoxygen == 2 .or. cfg%swcf == 3 .or. &
          cfg%swharv == 1   .or. cfg%swcompensate /= 0 .or.                 &
          cfg%swinter == 2  .or. cfg%swinter == 3 .or.                      &
          cfg%swrd /= 1     .or. cfg%swsalinity /= 0 .or. cfg%schedule_switch == 1) then
         call fatalerr_collected('cropfixed_init_from_config', &
            'Unsupported runtime branch reached on the TOML path. The validator ' // &
            'should have caught this earlier.')
         return
      end if

      ! ---- Copy scalars / switches to globals ------------------------
      idev   = cfg%idev
      lcc    = cfg%lcc
      tsumea = cfg%tsumea
      tsumam = cfg%tsumam
      tbase  = cfg%tbase

      kdif = cfg%kdif
      kdir = cfg%kdir

      swgc = cfg%swgc
      swcf = cfg%swcf
      swrd = cfg%swrd
      swdmi2rd = cfg%swdmi2rd
      swrdc    = cfg%swrdc

      swoxygen   = cfg%swoxygen
      swwrtnonox = cfg%swwrtnonox
      aeratecrit = cfg%aeratecrit
      hlim1      = cfg%hlim1
      hlim2u     = cfg%hlim2u
      hlim2l     = cfg%hlim2l

      swdrought = cfg%swdrought
      hlim3h    = cfg%hlim3h
      hlim3l    = cfg%hlim3l
      hlim4     = cfg%hlim4
      adcrh     = cfg%adcrh
      adcrl     = cfg%adcrl

      swsalinity = cfg%swsalinity
      saltmax    = cfg%saltmax
      saltslope  = cfg%saltslope
      salthead   = cfg%salthead

      swcompensate = cfg%swcompensate
      swstressor   = cfg%swstressor
      alphacrit    = cfg%alphacrit
      dcritrtz     = cfg%dcritrtz

      swinter = cfg%swinter
      cofab   = cfg%cofab

      schedule = cfg%schedule_switch

      dvsend = cfg%dvsend
      swharv = cfg%swharv

      ! ---- Reflection coefficients / crop resistance (legacy:2150-2160)
      ! Only ETref-friendly branch (swcf=1 or 3) is reachable; legacy
      ! sets albedo=0.23, rsc=70, rsw=0 in that branch. We read the
      ! values from cfg in case the user authored them; otherwise fall
      ! back to the same defaults.
      if (cfg%swcf == 1) then
         albedo = 0.23_real64
         rsc    = 70.0_real64
         rsw    = 0.0_real64
      end if

      ! ---- Copy tables ------------------------------------------------
      ! gctb (size up to 2*magrs in legacy; we copy what was authored).
      if (allocated(cfg%gctb))  call copy_pair_table(cfg%gctb,  gctb)
      if (allocated(cfg%cftb))  call copy_pair_table(cfg%cftb,  cftb)
      if (allocated(cfg%chtb))  call copy_pair_table(cfg%chtb,  chtb)
      if (allocated(cfg%cfeictb)) call copy_pair_table(cfg%cfeictb, cfeictb)
      if (allocated(cfg%rdtb))  call copy_pair_table(cfg%rdtb,  rdtb)

      ! rdctb is sized 22 in legacy; we always copy.
      if (allocated(cfg%rdctb)) then
         do i = 1, min(size(cfg%rdctb), size(rdctb))
            rdctb(i) = cfg%rdctb(i)
         end do
      end if

      ! ---- Tail-of-readcropfixed RUNTIME init: cumdens (legacy:2449-2478)
      if (cfg%swdrought == 1) then
         do i = 0, 100
            depth = 0.01_real64 * real(i, real64)
            rootdis(i*2 + 1) = depth
            rootdis(i*2 + 2) = afgen(rdctb, 22, depth)
         end do
         do i = 1, 202, 2
            cumdens(i) = rootdis(i)
         end do
         sum = 0.0_real64
         cumdens(2) = 0.0_real64
         do i = 4, 202, 2
            sum = sum + (rootdis(i-2) + rootdis(i)) * 0.5_real64 &
                      * (cumdens(i-1) - cumdens(i-3))
            cumdens(i) = sum
         end do
         do i = 2, 202, 2
            cumdens(i) = cumdens(i) / sum
         end do
      end if
   end subroutine cropfixed_init_from_config

   ! Copy a flat (dvs, value) pair source into a destination global of the
   ! same flat-array shape. Destination is a fixed-size global; copy up to
   ! min(size).
   subroutine copy_pair_table(src, dst)
      real(real64), intent(in)    :: src(:)
      real(real64), intent(inout) :: dst(:)
      integer :: i
      do i = 1, min(size(src), size(dst))
         dst(i) = src(i)
      end do
   end subroutine copy_pair_table

end module cropfixed_init_mod
```

- [ ] **Step 2: Register in `meson.build`**

In the top-level `meson.build`, add (just below the `crop_config_global.f90` line from Task 5, alphabetical-ish, before `cropgrowth.f90`):

```
    'src/crop/cropfixed_init.f90',
```

- [ ] **Step 3: Build, verify clean**

```bash
pixi run swap 2>&1 | tail -10
```

Expected: clean build. If a missing import is reported (e.g. `cumdens` not in `variables`), check the `use variables, only:` list and adjust.

- [ ] **Step 4: Write the unit test**

Create `tests/unit/crop/test_cropfixed_init.pf` (creates `tests/unit/crop/` directory):

```fortran
@test
subroutine test_cropfixed_init_writes_scalars()
   use funit
   use cropfixed_init_mod, only: cropfixed_init_from_config
   use cropfixed_config_mod, only: cropfixed_config_t
   use variables, only: idev, lcc, kdif, kdir, hlim3h, hlim3l, hlim4, &
                         adcrh, adcrl, cofab, swgc, swcf, swrd, swoxygen, &
                         swdrought, swsalinity, swcompensate, swinter
   type(cropfixed_config_t) :: cfg
   ! Case 6 grass.crp values
   cfg%idev = 1
   cfg%lcc  = 366
   cfg%kdif = 0.75d0
   cfg%kdir = 0.75d0
   cfg%swgc = 1
   cfg%swcf = 1
   cfg%swrd = 1
   cfg%swoxygen = 0
   cfg%swdrought = 1
   cfg%swsalinity = 0
   cfg%swcompensate = 0
   cfg%swinter = 1
   cfg%cofab = 0.25d0
   cfg%hlim3h = -200.0d0
   cfg%hlim3l = -800.0d0
   cfg%hlim4  = -8000.0d0
   cfg%adcrh  = 0.5d0
   cfg%adcrl  = 0.1d0
   ! rdctb required by tail-init even though we don't assert on cumdens here.
   allocate(cfg%rdctb(4))
   cfg%rdctb = [0.0d0, 1.0d0, 1.0d0, 1.0d0]

   call cropfixed_init_from_config(cfg, 1)

   @assertEqual(1,    idev)
   @assertEqual(366,  lcc)
   @assertEqual(0.75d0, kdif, 1.0d-12)
   @assertEqual(0.75d0, kdir, 1.0d-12)
   @assertEqual(1,    swgc)
   @assertEqual(1,    swcf)
   @assertEqual(1,    swrd)
   @assertEqual(0,    swoxygen)
   @assertEqual(1,    swdrought)
   @assertEqual(0,    swsalinity)
   @assertEqual(0,    swcompensate)
   @assertEqual(1,    swinter)
   @assertEqual(-200.0d0,  hlim3h, 1.0d-12)
   @assertEqual(-800.0d0,  hlim3l, 1.0d-12)
   @assertEqual(-8000.0d0, hlim4,  1.0d-12)
   @assertEqual(0.5d0,     adcrh,  1.0d-12)
   @assertEqual(0.1d0,     adcrl,  1.0d-12)
   @assertEqual(0.25d0,    cofab,  1.0d-12)
end subroutine

@test
subroutine test_cropfixed_init_cumdens_constant_density()
   use funit
   use cropfixed_init_mod, only: cropfixed_init_from_config
   use cropfixed_config_mod, only: cropfixed_config_t
   use variables, only: cumdens
   type(cropfixed_config_t) :: cfg
   ! Same case-6 fixture: rdctb = [0.0, 1.0, 1.0, 1.0] is a constant
   ! density of 1.0 everywhere. After integration and normalization,
   ! cumdens should be linear from 0 at the surface to 1 at the
   ! deepest depth. Check a few representative rows.
   cfg%idev = 1
   cfg%lcc = 366
   cfg%swrd = 1
   cfg%swdrought = 1
   allocate(cfg%rdctb(4))
   cfg%rdctb = [0.0d0, 1.0d0, 1.0d0, 1.0d0]
   call cropfixed_init_from_config(cfg, 1)
   ! cumdens(1) = 0.0 (depth at surface)
   @assertEqual(0.0d0,  cumdens(1),   1.0d-12)
   ! cumdens(2) = 0.0 (storage at surface, set explicitly)
   @assertEqual(0.0d0,  cumdens(2),   1.0d-12)
   ! cumdens(202) = 1.0 (normalized total)
   @assertEqual(1.0d0,  cumdens(202), 1.0d-12)
   ! Mid-point (i=102, depth=0.5) should be ~0.5 by symmetry of constant density.
   @assertEqual(0.5d0,  cumdens(102), 1.0d-2)
end subroutine
```

- [ ] **Step 5: Wire test in `tests/unit/meson.build`**

In `tests/unit/meson.build`, in the `pf_files` list, add (alphabetical-ish, after the `config/` entries, before `io/`):

```
        'crop/test_cropfixed_init.pf',
```

- [ ] **Step 6: Wire suite in `tests/unit/testSuites.inc`**

Append:

```
ADD_TEST_SUITE(test_cropfixed_init_suite)
```

- [ ] **Step 7: Run unit tests, verify pass**

```bash
pixi run test-pfunit 2>&1 | tail -10
```

Expected: `Ok: 1, Fail: 0` including the 2 new tests.

- [ ] **Step 8: Commit**

```bash
git add src/crop/cropfixed_init.f90 meson.build \
        tests/unit/crop/test_cropfixed_init.pf \
        tests/unit/meson.build tests/unit/testSuites.inc
git commit -m "feat(crop): cropfixed_init module replaces readcropfixed runtime math

New module copies cropfixed_config_t fields to legacy variables
module globals (idev, kdif, gctb, cftb, hlim*, rdctb, rdtb, ...) and
runs the tail-of-readcropfixed RUNTIME math (cumdens normalized
cumulative root density). Defense-in-depth runtime guards mirror the
validator stub-errors per ADR 0015.

Two pFUnit tests verify scalar copy + cumdens for case-6 constant
density rdctb (analytic check: cumdens(102) ≈ 0.5 for constant
distribution)."
```

---

## Task 7: Wire `cropgrowth.f90:397` dispatch

Replace the unconditional `call readcropfixed(...)` with a `populated`-sentinel dispatch.

**Files:**
- Modify: `src/crop/cropgrowth.f90` (line 397)

- [ ] **Step 1: Read the current callsite**

```bash
sed -n '390,405p' src/crop/cropgrowth.f90
```

- [ ] **Step 2: Replace the call**

In `src/crop/cropgrowth.f90`, find the `cropfixed` subroutine's `case (1)` block. The current line ~397 reads:

```fortran
! --- read crop data
      call readcropfixed (icrop,cropfil(icrop),lcc,swhydrlift)
```

Replace with:

```fortran
! --- read crop data: dispatch on per-rotation typed-config cache
!     (ADR 0016). Falls back to legacy reader for rotations whose
!     .crp.toml is not yet authored or for rotation types not yet
!     ported (Phase 1: only type=1 cropfixed; Phases 2/3 add types
!     2 and 3). Teardown: end of Phase 4 removes the else-branch.
      block
         use crop_config_global_mod, only: crop_config_global
         use cropfixed_init_mod, only: cropfixed_init_from_config
         logical :: use_cache
         use_cache = .false.
         if (associated(crop_config_global)) then
            if (allocated(crop_config_global%rotation_loaded)) then
               if (icrop >= 1 .and. icrop <= size(crop_config_global%rotation_loaded)) then
                  if (crop_config_global%rotation_loaded(icrop)) use_cache = .true.
               end if
            end if
         end if
         if (use_cache) then
            call cropfixed_init_from_config(crop_config_global%rotation_fixed(icrop), icrop)
            ! lcc + swhydrlift are out-arguments of legacy readcropfixed;
            ! mirror them from the typed config.
            lcc        = crop_config_global%rotation_fixed(icrop)%lcc
            swhydrlift = 0   ! Phase 1 stub-errors swdrought=2; legacy reads swhydrlift only in that branch
         else
            call readcropfixed (icrop,cropfil(icrop),lcc,swhydrlift)   ! transitional fallback
         end if
      end block
```

- [ ] **Step 3: Build, verify clean**

```bash
pixi run swap 2>&1 | tail -10
```

Expected: clean build.

- [ ] **Step 4: Run unit tests, verify green**

```bash
pixi run test-pfunit 2>&1 | tail -5
```

Expected: `Ok: 1, Fail: 0`.

- [ ] **Step 5: Run regression, verify case 6 still green**

```bash
pixi run regression 2>&1 | tail -20
```

Expected: 5/5 cases green. Case 6 now exercises the new path via `cropfixed_init_from_config`. If case 6 produces different output than reference, investigate before committing — the parity audit (Task 1) is the first place to look.

- [ ] **Step 6: Commit**

```bash
git add src/crop/cropgrowth.f90
git commit -m "refactor(crop): cropfixed task=1 dispatches on per-rotation cache

Replaces unconditional call to legacy readcropfixed with a dispatch
on rotation_loaded(icrop). When the slot is loaded (case 6 type-1
rotations after Phase 1), calls cropfixed_init_from_config; else
falls back to legacy readcropfixed (transitional, removed end of
Phase 4). lcc and swhydrlift mirrored from typed config when the
cache path is taken."
```

---

## Task 8: Smoke test — verify case 6 runs without `grass.crp`

Pre-flight before deletion: rename `grass.crp` and re-run regression. Verification-only, no commit.

- [ ] **Step 1: Rename `grass.crp` in submodule**

```bash
cd /home/zawadzkim/Code/swap/tests/swap-cases/toml/6.surfacewater && mv grass.crp grass.crp.disabled && cd /home/zawadzkim/Code/swap
```

- [ ] **Step 2: Run regression**

```bash
pixi run regression 2>&1 | tail -20
```

Expected: 5/5 green. If case 6 fails with "Cannot open grass.crp" or similar, the runtime path is not yet decoupled — investigate before proceeding to Task 9.

- [ ] **Step 3: Restore the file**

```bash
cd /home/zawadzkim/Code/swap/tests/swap-cases/toml/6.surfacewater && mv grass.crp.disabled grass.crp && cd /home/zawadzkim/Code/swap
```

- [ ] **Step 4: Verify submodule clean (re Task 8 — rename should round-trip)**

```bash
cd /home/zawadzkim/Code/swap/tests/swap-cases && git status --short && cd ..
```

The submodule may have pre-existing dirty files in OTHER case dirs (unrelated to this port). Confirm `toml/6.surfacewater/grass.crp` is NOT in the dirty list (rename round-tripped cleanly). If `grass.crp` is in the dirty list, the rename did not round-trip — investigate.

- [ ] **Step 5: No commit** — verification-only.

---

## Task 9: Delete `grass.crp` from case 6 TOML dir

Final removal. Submodule pair commit.

**Files:**
- Delete (in submodule): `tests/swap-cases/toml/6.surfacewater/grass.crp`

- [ ] **Step 1: Remove the file in the submodule**

```bash
cd /home/zawadzkim/Code/swap/tests/swap-cases && git rm toml/6.surfacewater/grass.crp && cd /home/zawadzkim/Code/swap
```

- [ ] **Step 2: Confirm parity tests still pass**

The legacy `read_legacy_cropfixed` parity tests (in `tests/unit/io/toml/legacy_crop_helper.f90` and downstream `.pf` files) chdir to `tests/swap-cases/6.surfacewater/` (the LEGACY case dir, not the toml dir). That dir's `grass.crp` is unaffected.

```bash
pixi run test-pfunit 2>&1 | tail -5
```

Expected: `Ok: 1, Fail: 0`.

- [ ] **Step 3: Run regression**

```bash
pixi run regression 2>&1 | tail -10
```

Expected: 5/5 green.

- [ ] **Step 4: Submodule inner-commit**

Use file-scoped commit to avoid bundling pre-existing dirty state.

```bash
cd /home/zawadzkim/Code/swap/tests/swap-cases
git commit toml/6.surfacewater/grass.crp -m "chore(toml/6.surfacewater): remove legacy grass.crp

Cropfixed configuration is now fully expressed in grass.crp.toml.
The legacy ASCII grass.crp remains in tests/swap-cases/6.surfacewater/
for the legacy executable and the cropfixed parity test fixture."
cd /home/zawadzkim/Code/swap
```

- [ ] **Step 5: Outer-repo bump**

```bash
git add tests/swap-cases
git commit -m "chore(submodule): bump tests/swap-cases — remove grass.crp from case 6

Phase 1 of the .crp port is complete. Case 6 surfacewater runs
end-to-end via the typed pipeline + grass.crp.toml + cropfixed_init;
no .crp ASCII file is opened on the runtime path for type-1 rotations."
```

---

## Task 10: Documentation

Update the two doc files referenced by the spec.

**Files:**
- Modify: `docs/csv-companion-files.md`
- Modify: `docs/configuration-schema.md` (only if its `[cropfixed]` schema is missing the new fields)

- [ ] **Step 1: Read current state**

```bash
grep -n "\.crp\|cropfixed\|grass\.crp" docs/csv-companion-files.md | head -10
grep -n "cropfixed\|swdrought\|swoxygen" docs/configuration-schema.md | head -10
```

- [ ] **Step 2: Edit `docs/csv-companion-files.md`**

Find the "Path resolution and staging" section (around lines 92-106). Update the bullet listing what lives in the case working directory to reflect that the legacy `*.crp` files are no longer staged for type-1 rotations after Phase 1. The exact existing wording may differ; preserve everything else and adjust the `*.crp` reference. Example post-edit form:

```markdown
- The case working directory is `tests/swap-cases/toml/<N>.<case>/`. It is
  self-contained — every file SWAP reads at runtime lives there: `swap.toml`,
  `swap.dra.toml`, `*.crp.toml`, all `*.csv` companions, and
  `swap_linux.swp.template` (staged to `swap.swp` per run). The legacy `*.crp`
  ASCII files for crop types 2 and 3 (read by sub-readers in `cropgrowth.f90`
  until Phases 2 and 3 of the `.crp` port land) also live there. Type-1
  cropfixed rotations are read entirely from `*.crp.toml` post-Phase-1.
```

If the file is structured differently, adapt to its actual structure.

- [ ] **Step 3: Edit `docs/configuration-schema.md` if needed**

Read the existing `[cropfixed]` section (or the `[crop.rotation]` per-file schema, depending on how the doc is structured). If the new fields (`swprep`, `swsow`, `swgerm`, `dvsend`, `swharv`, `tsumea`, `tsumam`, `tbase`, `swgc`, `swcf`, `swrd`, `swdmi2rd`, `swrdc`, `swoxygen`, `swwrtnonox`, `swdrought`, `swsalinity`, `swcompensate`, `swstressor`, `swinter`, the `aeratecrit`/`saltmax`/`saltslope`/`salthead`/`alphacrit`/`dcritrtz` reals, plus the table fields `gctb`/`cftb`/`chtb`/`rdtb`/`rdctb`) are missing, add them to the schema table.

If the doc is up-to-date already, skip this edit and note "configuration-schema.md unchanged" in the commit message.

- [ ] **Step 4: Commit**

```bash
git add docs/csv-companion-files.md docs/configuration-schema.md
git commit -m "docs: cropfixed port (Phase 1) docs updates

csv-companion-files.md: type-1 rotations no longer stage legacy .crp.
configuration-schema.md: cropfixed schema gains schema-1:1 fields."
```

---

## Acceptance gate

After Task 10:

- [ ] `pixi run test-pfunit 2>&1 | grep "Ok:|Fail:"` → `Ok: 1, Fail: 0`.
- [ ] `pixi run regression 2>&1 | tail -5` → `5 passed, 0 failed`.
- [ ] `git grep "call readcropfixed" src/` → returns the legacy fallback in `cropgrowth.f90` only (the `else` branch). Other reachability gone.
- [ ] `git grep "call readcropfixed" tests/` → still returns parity-test references (intentional).
- [ ] `ls tests/swap-cases/toml/6.surfacewater/grass.crp` → file does not exist.
- [ ] `tests/swap-cases/toml/6.surfacewater/grass.crp.toml` → schema-1:1 with legacy.
- [ ] `git log --oneline | head -15` → ~10-12 commits since the start of Phase 1.

---

## Self-review notes

**Spec coverage:**
- Spec Unit 1 (Schema 1:1) → Task 2.
- Spec Unit 2 (Parser extension) → Task 3.
- Spec Unit 3 (Loader extension) → already implemented; no task needed (noted in plan header).
- Spec Unit 4 (Sentinel) → already implemented as `rotation_loaded(:)` on parent; no task needed.
- Spec Unit 5 (Runtime init module) → Task 6.
- Spec Unit 6 (Module-level config reference) → Task 5.
- Spec Unit 7 (Wiring change) → Task 7.
- Spec test plan → Tasks 2/3/4/6 unit tests + smoke test (Task 8) + acceptance gate.
- Spec teardown table → covered in ADR 0016 + plan task notes (Task 5 calls `crop_config_global` transitional with explicit teardown trigger; Task 7 marks the `else` fallback transitional with the same trigger).

**Type / signature consistency:**
- `cropfixed_init_from_config(cfg, icrop)` — same signature in Tasks 6, 7, and the test in Task 6. Verified.
- `crop_config_global` is a `pointer` of type `crop_config_t` in Tasks 5 and 7. Verified.
- `rotation_loaded(:)` and `rotation_fixed(:)` — actual code names, used consistently throughout.
- `schedule_switch` (the new top-level scheduling integer) vs the existing `schedule` (irrigation_schedule_t). Renamed to avoid shadowing the existing field. Used consistently in Task 2 (validator), Task 3 (parser, key `'schedule'` under section `'scheduling'`), Task 4 (TOML fixture), Task 6 (`cropfixed_init` writes `schedule = cfg%schedule_switch` to the legacy global).

**Submodule discipline:** Tasks 4 and 9 each use file-scoped `git commit <path>` to avoid bundling pre-existing dirty state. Task 8 verifies rename round-trip cleanly without leaving submodule dirty.

**Naming reconciliation noted at top of plan** — the spec's `populated`/`rotation_cropfixed` proposal was superseded by the existing `rotation_loaded(:)`/`rotation_fixed(:)` infrastructure. Plan uses actual names; spec stays as a record of the design conversation.
