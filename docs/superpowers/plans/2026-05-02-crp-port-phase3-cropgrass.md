# `.crp` port — Phase 3 (cropgrass via cases 4 + 2) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Port case 4 (oxygenstress)'s `grassd.crp` (type 3, detailed grass) to
TOML so the new executable's runtime path opens no `.crp` ASCII file for type-3
rotations. Case 2 (grassgrowth) serves as secondary verification — both cases must
regress green throughout. The legacy `readgrass` reader stays alive as a parity
fixture.

**Architecture:** Mirrors Phase 1's shape (extend schema 1:1, push validation/init
into typed config + a new init module, wire the runtime callsite to dispatch on the
`rotation_loaded` sentinel, fall back to the legacy reader for unported types). The
per-rotation crop config cache is already in place (Phase 1). Phase 3 adds the
type=3 slot population and the `cropgrass_init_from_config` subroutine. Per ADR
0016, the dispatch wiring and `crop_config_global` pointer are reused without
modification.

**Tech Stack:** Fortran 2008 (gfortran 13), pFUnit unit tests, Meson + Pixi, tomlf
TOML parser.

**Spec:** [docs/superpowers/specs/2026-05-02-crp-port-phase3-cropgrass-design.md](../specs/2026-05-02-crp-port-phase3-cropgrass-design.md)

**ADR:** [docs/adr/0016-per-rotation-crop-config-cache.md](../../adr/0016-per-rotation-crop-config-cache.md)

**Sequencing dependency:** Phase 1 (cropfixed via case 6) must be in the repository
before Phase 3 is started. Phase 3 reuses `crop_config_global` and the
`rotation_loaded(:)` sentinel introduced by Phase 1.

**Naming reconciliation (spec → plan):**

The spec uses proposed names. Confirm actual names in `src/config/crop_config.f90`
before each task. Based on Phase 1's precedent, the following divergences are
possible:

| Spec name (proposal) | Expected actual code name |
|---|---|
| `rotation_cropgrass(:)` | `rotation_grass(:)` on `crop_config_t` |
| `populated :: logical` per-config | `rotation_loaded(:)` LOGICAL array on `crop_config_t` |
| `read_cropgrass_file_toml` | May be named differently if Phase 1 established a convention |

**Always verify actual names before writing code in Tasks 3, 7, 8, and 9.**

---

## File map

**Source — modify:**
- `src/config/cropgrass_config.f90` — extend schema 1:1; add stub-error
  validators; add `populated` sentinel
- `src/io/toml/read_cropgrass_toml.f90` — extend parser for all new fields
  and tables; add `read_cropgrass_file_toml` wrapper; set `populated=.true.`
- `src/io/toml/read_crop_toml.f90` — ensure type=3 dispatch allocates
  `rotation_grass(:)` and calls `read_cropgrass_file_toml` per slot
- `src/crop/cropgrowth.f90` — at line 2068, add `if/else` dispatch on the
  `populated` sentinel around the `readgrass` call
- `meson.build` — register `src/crop/cropgrass_init.f90`
- `tests/unit/meson.build` — register new test files
- `tests/unit/testSuites.inc` — register new test suites

**Source — create:**
- `src/crop/cropgrass_init.f90` — runtime init from typed config (replaces
  `readgrass`'s runtime side-effects on the TOML path)

**Tests — create / modify:**
- `tests/unit/config/test_cropgrass_config.pf` — extend with stub-error tests
  + table-shape tests + new scalar range tests
- `tests/unit/io/toml/test_read_cropgrass_toml.pf` — extend with round-trip of
  all new scalars and tables
- `tests/unit/io/toml/test_load_swap_config.pf` — extend with case-4 and
  case-2 assertions that all rotation slots are loaded
- `tests/unit/crop/test_cropgrass_init.pf` — NEW: verifies init writes correct
  globals + cumdens hand-computed value + swoxygen=2 globals + swcompensate=1
  globals

**Docs — modify:**
- `docs/csv-companion-files.md` — note the `.crp.toml` path-resolution rule
  (if not already done by Phase 1)

**Test fixtures — modify (in submodule, two pairs):**
- `tests/swap-cases/toml/4.oxygenstress/grassd.crp.toml` — replace 53-line
  skeleton with full 1:1 reproduction (Task 5)
- `tests/swap-cases/toml/2.grassgrowth/grassd.crp.toml` — replace 59-line
  skeleton with full 1:1 reproduction (Task 6)
- `tests/swap-cases/toml/4.oxygenstress/grassd.crp` — DELETE in Task 12a
- `tests/swap-cases/toml/2.grassgrowth/grassd.crp` — DELETE in Task 12b

---

## Submodule discipline

`tests/swap-cases/` is a git submodule. Inner-commit + outer-repo bump pair is
non-negotiable. Never one without the other. Affected tasks: Task 5 (author
case-4 TOML), Task 6 (author case-2 TOML), Task 12a (delete case-4 `.crp`),
Task 12b (delete case-2 `.crp`). That is four submodule pair-commits total.

When committing in the submodule, use file-scoped git commands (e.g.
`git commit toml/4.oxygenstress/grassd.crp.toml -m "..."`) so pre-existing
dirty state from other case directories is not bundled in.

---

## Pre-flight commands

Confirm the starting state:

```bash
cd /home/zawadzkim/Code/swap
git log --oneline -3
pixi run test-pfunit 2>&1 | tail -5
pixi run regression 2>&1 | tail -10
```

Expected: unit tests all green, 5 regression cases passed. If either is red,
do not proceed — the baseline is broken.

Also verify Phase 1 is in the repository:

```bash
git log --oneline | grep -i "phase.1\|cropfixed"
grep -l "crop_config_global" src/crop/*.f90
```

If `crop_config_global` does not exist, Phase 1 has not landed. Stop and wait
for Phase 1.

Also confirm the actual name of the grass-rotation array and sentinel:

```bash
grep -n "rotation_grass\|rotation_cropgrass\|rotation_loaded\|populated" \
    src/config/crop_config.f90 | head -30
```

Record the actual name; use it verbatim in every subsequent task.

---

## Task 1: Audit `readgrass`

Read `src/io/readswap.f90` lines 3437-4270 end-to-end and produce a
classification table. Docs only — no source changes.

**Files:**
- Create: `docs/phase-4f-readgrass-audit.md`

- [ ] **Step 1: Read the entire `readgrass` body**

```bash
sed -n '3437,4270p' src/io/readswap.f90 | cat
```

Or use the Read tool with `offset=3437, limit=833`.

- [ ] **Step 2: Classify each non-blank, non-comment line into one of five buckets**

Bucket definitions (same as Phase 1 / ADR 0015):

1. **READ** — `rd*` calls. Goes away with port.
2. **VALIDATE** — `if (...) call fatalerr(...)`. Push to validators.
3. **NORMALIZE** — unit or coordinate conversion. Push to finalizers/adapter.
4. **RUNTIME** — non-read initialization (cumdens build, dateharvest bounds
   fill, `irrigation(1)` call, etc.). Goes into `cropgrass_init`.
5. **GUARDED** — branches Phase 3 does not support: `swoxygen=2
   swoxygentype=2`, `swcompensate=2`, `swdrought=2`, `swsalinity≠0`,
   `swco2=1`, `swlossgrz=1`, `swlossmow=1`, `seqgrazmow(i)≠2`.

Note that `swoxygen=2 swoxygentype=1` (physical Bartholomeus) IS supported in
Phase 3 — classify those lines as READ/RUNTIME, not GUARDED.

Note the special structure of the mowing/grazing block: `readgrass` reads
`swharvest` twice (once in the grazing block when `swGrz`, once in the mowing
block when `swMow`). In Phase 3 only `swMow=.true.` (SEQGRAZMOW=all-2) is
supported; the grazing block's `swharvest` read is therefore in the GUARDED
bucket.

- [ ] **Step 3: Write the audit table**

Create `docs/phase-4f-readgrass-audit.md` with a brief intro + a Markdown
table:

```markdown
| Lines       | Bucket    | Notes                                          |
| ----------- | --------- | ---------------------------------------------- |
| 3437-3469   | READ      | Open .crp, argument decl, locals               |
| 3508-3557   | READ      | swcf + cftb/chtb/cfeictb tables                |
| 3537-3557   | GUARDED   | swcf=3 wet-crop branch                         |
| 3560-3585   | READ      | swinter + cofab / Gash / storage-cap           |
| 3563-3585   | GUARDED   | swinter=2 (Gash) and swinter=3 (storage-cap)   |
| 3588-3598   | READ      | albedo/rsc/rsw for swcf=2 path                 |
| 3604-3613   | READ      | tdwi, laiem, rgrlai, swtsum (+ tsumtemp etc)   |
| 3617-3620   | READ      | slatb, ssa, span, tbase                        |
| 3622-3628   | READ      | kdif, kdir, eff, amaxtb, tmpftb, tmnftb        |
| 3631-3633   | READ      | cvl, cvr, cvs                                  |
| 3636-3640   | READ      | q10, rml, rmr, rms, rfsetb                     |
| 3643-3645   | READ      | frtb, fltb, fstb                               |
| 3648-3650   | READ      | perdl, rdrrtb, rdrstb                          |
| 3652-3723   | READ+GUARDED | swoxygen block                              |
| 3679-3712   | READ      | swoxygen=2 swoxygentype=1 (Bartholomeus)        |
| 3706-3711   | GUARDED   | swoxygentype=2 (reproduction functions)         |
| 3714-3723   | READ      | swwrtnonox, aeratecrit                         |
| 3725-3754   | READ+GUARDED | swdrought block                             |
| 3740-3754   | GUARDED   | swdrought=2 De Jong van Lier                   |
| 3756-3778   | GUARDED   | swsalinity block (flsolute gate)               |
| 3780-3830   | READ+GUARDED | swcompensate block (swcompensate=2 GUARDED) |
| 3832-3873   | READ      | swrdc, rdctb, swrd (+ rdtb/rdi/rri/rdc/rlwtb) |
| 3875-3885   | READ      | relmf, swpotrelmf                              |
| 3887-3982   | READ+GUARDED | seqgrazmow + grazing block (GUARDED when swGrz)|
| 3984-4037   | READ+GUARDED | mowing block                                |
| 4039-4047   | READ+GUARDED | schedule + hlim/adcr for schedule+swdrought=2 |
| 4049-4083   | GUARDED   | swco2=1 block + CO2 tables                     |
| 4085-4089   | RUNTIME   | close file; irrigation(1) when schedule=1      |
| 4091-4121   | RUNTIME   | cumdens build (trapezium sum over rdctb)        |
| 4123-4269   | GUARDED   | swinco=3 .END-file resume path                 |
```

Include a Summary section at the bottom with bucket totals and a note
distinguishing the two-path swoxygen=2 block.

- [ ] **Step 4: Commit the audit doc**

```bash
git add docs/phase-4f-readgrass-audit.md
git commit -m "docs(audit): readgrass line-by-line classification for cropgrass port"
```

---

## Task 2: Extend `cropgrass_config_t` schema 1:1

Add the ~55 missing scalar fields and ~15 table fields to
`cropgrass_config_t`, plus stub-error validators for unsupported branch values.
Schema 1:1 philosophy: every field `readgrass` reads has a TOML home, even when
the parent switch is at a stubbed value.

**Files:**
- Modify: `src/config/cropgrass_config.f90`
- Test: `tests/unit/config/test_cropgrass_config.pf`

- [ ] **Step 1: Read the existing schema and validator**

```bash
cat -n src/config/cropgrass_config.f90
```

Note the current line count and the end of the `type :: cropgrass_config_t`
block. Record the line numbers for the `contains` line and the end of the
existing `cropgrass_config_validate` subroutine.

- [ ] **Step 2: Write failing stub-error tests**

Append to `tests/unit/config/test_cropgrass_config.pf`. Each `@test` is a
separate subroutine; check existing file style before pasting.

```fortran
! Phase 3 cropgrass port — stub-error validators for unsupported branches.

@test
subroutine test_cropgrass_swdrought2_rejected()
   use funit
   use cropgrass_config_mod, only: cropgrass_config_t
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   type(cropgrass_config_t) :: c
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
subroutine test_cropgrass_swcompensate2_rejected()
   use funit
   use cropgrass_config_mod, only: cropgrass_config_t
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   type(cropgrass_config_t) :: c
   type(error_collection_t) :: errors
   integer :: i
   logical :: found
   c%swcompensate = 2
   call c%validate(errors)
   found = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD .and. &
          index(errors%items(i)%message, 'swcompensate=2') > 0) found = .true.
   end do
   @assertTrue(found)
end subroutine

@test
subroutine test_cropgrass_swsalinity1_rejected()
   use funit
   use cropgrass_config_mod, only: cropgrass_config_t
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   type(cropgrass_config_t) :: c
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
subroutine test_cropgrass_swinter2_rejected()
   use funit
   use cropgrass_config_mod, only: cropgrass_config_t
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   type(cropgrass_config_t) :: c
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
subroutine test_cropgrass_swco2_1_rejected()
   use funit
   use cropgrass_config_mod, only: cropgrass_config_t
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   type(cropgrass_config_t) :: c
   type(error_collection_t) :: errors
   integer :: i
   logical :: found
   c%swco2 = 1
   call c%validate(errors)
   found = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD .and. &
          index(errors%items(i)%message, 'swco2=1') > 0) found = .true.
   end do
   @assertTrue(found)
end subroutine

@test
subroutine test_cropgrass_swlossmow1_rejected()
   use funit
   use cropgrass_config_mod, only: cropgrass_config_t
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   type(cropgrass_config_t) :: c
   type(error_collection_t) :: errors
   integer :: i
   logical :: found
   c%swlossmow = 1
   call c%validate(errors)
   found = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD .and. &
          index(errors%items(i)%message, 'swlossmow=1') > 0) found = .true.
   end do
   @assertTrue(found)
end subroutine

@test
subroutine test_cropgrass_swlossgrz1_rejected()
   use funit
   use cropgrass_config_mod, only: cropgrass_config_t
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   type(cropgrass_config_t) :: c
   type(error_collection_t) :: errors
   integer :: i
   logical :: found
   c%swlossgrz = 1
   call c%validate(errors)
   found = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD .and. &
          index(errors%items(i)%message, 'swlossgrz=1') > 0) found = .true.
   end do
   @assertTrue(found)
end subroutine

@test
subroutine test_cropgrass_swoxygentype2_rejected()
   use funit
   use cropgrass_config_mod, only: cropgrass_config_t
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   type(cropgrass_config_t) :: c
   type(error_collection_t) :: errors
   integer :: i
   logical :: found
   c%swoxygen     = 2
   c%swoxygentype = 2
   call c%validate(errors)
   found = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD .and. &
          index(errors%items(i)%message, 'swoxygentype=2') > 0) found = .true.
   end do
   @assertTrue(found)
end subroutine

@test
subroutine test_cropgrass_seqgrazmow_grazing_rejected()
   use funit
   use cropgrass_config_mod, only: cropgrass_config_t
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   type(cropgrass_config_t) :: c
   type(error_collection_t) :: errors
   integer :: i
   logical :: found
   ! Encode a period 1 = grazing in the sequence
   c%nseqgrazmow = 1
   allocate(c%seqgrazmow(1))
   c%seqgrazmow(1) = 1
   call c%validate(errors)
   found = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD .and. &
          index(errors%items(i)%message, 'seqgrazmow') > 0) found = .true.
   end do
   @assertTrue(found)
end subroutine

@test
subroutine test_cropgrass_case4_supported_values_pass()
   use funit
   use cropgrass_config_mod, only: cropgrass_config_t
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   type(cropgrass_config_t) :: c
   type(error_collection_t) :: errors
   integer :: i
   logical :: any_stub
   ! Case 4 supported values for the new switches.
   c%swcf         = 2
   c%swinter      = 1
   c%swtsum       = 1
   c%swrd         = 2
   c%swdmi2rd     = 1
   c%swrdc        = 0
   c%swoxygen     = 2
   c%swoxygentype = 1    ! physical — supported
   c%swwrtnonox   = 1
   c%swdrought    = 1
   c%swsalinity   = 0
   c%swcompensate = 1    ! Jarvis — supported
   c%swstressor   = 1
   c%swdmmow      = 2
   c%swharv       = 1
   c%swdmgrz      = 2
   c%swlossgrz    = 0
   c%swlossmow    = 0
   c%nseqgrazmow  = 3
   allocate(c%seqgrazmow(3))
   c%seqgrazmow   = 2    ! all mowing — supported
   c%swco2        = 0
   call c%validate(errors)
   any_stub = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD) any_stub = .true.
   end do
   @assertFalse(any_stub)
end subroutine

@test
subroutine test_cropgrass_case2_supported_values_pass()
   use funit
   use cropgrass_config_mod, only: cropgrass_config_t
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   type(cropgrass_config_t) :: c
   type(error_collection_t) :: errors
   integer :: i
   logical :: any_stub
   ! Case 2 supported values for the new switches.
   c%swcf         = 2
   c%swinter      = 1
   c%swtsum       = 1
   c%swrd         = 2
   c%swoxygen     = 1    ! Feddes — supported
   c%swwrtnonox   = 1
   c%swdrought    = 1
   c%swsalinity   = 0
   c%swcompensate = 0
   c%swharv       = 2
   c%nmow         = 30
   allocate(c%mowing_dates(30))
   c%mowing_dates = 1.0_real64  ! placeholder values (shape check only)
   c%nseqgrazmow  = 3
   allocate(c%seqgrazmow(3))
   c%seqgrazmow   = 2    ! all mowing — supported
   c%swco2        = 0
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
pixi run test-pfunit 2>&1 | tail -20
```

Expected: at least 9 stub-error tests fail (new validators + new fields not
yet declared). The "supported values" tests may fail to compile until the new
fields exist.

- [ ] **Step 4: Add new fields to `cropgrass_config_t`**

In `src/config/cropgrass_config.f90`, extend the `type :: cropgrass_config_t`
block. Preserve all existing fields. Add the following grouped after the
existing `cofab` line and before `swharv`, inserting into appropriate logical
positions.

The full extended field list (Phase 3 additions):

```fortran
      ! ====================================================================
      ! Phase 3 (.crp port) additions
      ! ====================================================================

      ! Crop state initialisation (Part 2a)
      real(real64) :: tdwi   = 1000.0_real64  !! Initial total crop dry weight [kg/ha]
      real(real64) :: laiem  = 0.63_real64    !! Leaf area index at emergence
      real(real64) :: rgrlai = 0.007_real64   !! Max relative increase in LAI

      ! Start-of-growth trigger (Part 2b)
      integer      :: swtsum    = 1            !! 0=none, 1=air-temp sum, 2=soil temp
      real(real64) :: tsumtemp  = 0.0_real64   !! Soil temp threshold (when swtsum=2)
      real(real64) :: tsumdepth = 0.0_real64   !! Soil depth for temp (when swtsum=2)
      integer      :: tsumtime  = 0            !! Consecutive days above threshold

      ! Green area (Part 3)
      real(real64) :: ssa  = 0.0_real64       !! Specific stem area [ha/kg]
      real(real64) :: span = 30.0_real64      !! Leaf life span [d]
      ! (tbase already declared above)
      real(real64), allocatable :: slatb(:)   !! Specific leaf area vs DNR (flat pairs)

      ! Assimilation tables (Part 4)
      real(real64), allocatable :: amaxtb(:)  !! Max CO2 assim rate vs DNR
      real(real64), allocatable :: tmpftb(:)  !! AMAX reduction vs avg temp
      real(real64), allocatable :: tmnftb(:)  !! AMAX reduction vs min temp

      ! Biomass conversion (Part 5)
      real(real64) :: cvl = 0.685_real64      !! Leaves
      real(real64) :: cvr = 0.694_real64      !! Roots
      real(real64) :: cvs = 0.662_real64      !! Stems

      ! Maintenance respiration (Part 6)
      real(real64) :: q10 = 2.0_real64        !! Q10 factor
      real(real64) :: rml = 0.03_real64       !! Leaves
      real(real64) :: rmr = 0.015_real64      !! Roots
      real(real64) :: rms = 0.015_real64      !! Stems
      real(real64), allocatable :: rfsetb(:)  !! Senescence reduction vs DNR

      ! Partitioning (Part 7)
      real(real64), allocatable :: frtb(:)    !! Roots fraction vs DNR
      real(real64), allocatable :: fltb(:)    !! Leaves fraction vs DNR
      real(real64), allocatable :: fstb(:)    !! Stems fraction vs DNR

      ! Death rates (Part 8)
      real(real64) :: perdl = 0.05_real64     !! Max death rate leaves under water stress
      real(real64), allocatable :: rdrrtb(:)  !! Root death rate vs DNR
      real(real64), allocatable :: rdrstb(:)  !! Stem death rate vs DNR

      ! Crop factor / height switch (Part 1)
      integer      :: swcf   = 2              !! 1=crop factor, 2=crop height, 3=wet (stub)
      real(real64) :: albedo = 0.23_real64    !! Reflection coeff (when swcf=2)
      real(real64) :: rsw    = 0.0_real64     !! Canopy resistance intercepted water

      ! Root depth and density switches (Part 9)
      integer :: swrd     = 2                 !! 1=DVS table, 2=daily increase, 3=biomass
      integer :: swdmi2rd = 0                 !! 0=assimilates, 1=DM increase
      integer :: swrdc    = 0                 !! 0=unmodified, 1=modified by water extract
      real(real64) :: wrtmax = 3000.0_real64  !! Max root weight (when swrd=3)
      real(real64), allocatable :: rdtb(:)    !! Rooting depth vs DNR (when swrd=1)
      real(real64), allocatable :: rlwtb(:)   !! Rooting depth vs root weight (swrd=3)

      ! Oxygen stress additions (Part 10 — extending existing fields)
      integer      :: swoxygen   = 1          !! 0=none, 1=Feddes, 2=Bartholomeus
      integer      :: swwrtnonox = 0          !! Check aerobic conditions for root growth
      real(real64) :: aeratecrit = 1.0e-4_real64 !! Aerobic threshold for root extension

      ! Bartholomeus physical sub-model (swoxygentype=1 supported; =2 stub)
      integer      :: swoxygentype          = 1
      real(real64) :: q10_microbial         = 2.8_real64
      real(real64) :: specific_resp_humus   = 1.6e-3_real64
      real(real64) :: srl                   = 383571.0_real64
      integer      :: swrootradius          = 2  !! 1=calculated, 2=given
      real(real64) :: dry_mat_cont_roots    = 0.075_real64
      real(real64) :: air_filled_root_por   = 0.05_real64
      real(real64) :: spec_weight_root_tissue = 1.0e3_real64
      real(real64) :: var_a                 = 4.175e-10_real64
      real(real64) :: root_radiusO2         = 0.000075_real64

      ! Drought stress switch (Part 11)
      integer :: swdrought = 1               !! 1=Feddes (supported), 2=De Jong (stub)

      ! Salinity stress switch (Part 12)
      integer :: swsalinity = 0             !! 0=none (others stub)

      ! Root water uptake compensation (Part xx)
      integer      :: swcompensate = 0      !! 0=none, 1=Jarvis (supported), 2=Walsum (stub)
      integer      :: swstressor   = 1      !! Stressors to compensate
      real(real64) :: alphacrit    = 1.0_real64  !! Jarvis critical index
      real(real64) :: dcritrtz     = 0.0_real64  !! Walsum threshold (stub)

      ! Interception switch (Part 13)
      integer :: swinter = 1                !! 0=none, 1=Von Hoyningen (supported), 2/3=stub

      ! Management general (MANAGEMENT SECTION Part 1)
      integer :: nseqgrazmow = 20           !! Number of periods in SEQGRAZMOW
      integer, allocatable :: seqgrazmow(:) !! 1=graze, 2=mow, 3=dewool (only 2 supported)
      real(real64) :: mowrest  = 700.0_real64  !! Remaining DM after mowing [kg/ha]
      real(real64) :: dewrest  = 850.0_real64  !! Remaining DM after dewooling (stub)
      integer      :: swpotrelmf = 1        !! 1=theoretical, 2=attainable yield
      real(real64) :: relmf = 1.0_real64    !! Relative management factor

      ! Mowing DM-threshold flexible table (when swdmmow=2)
      real(real64), allocatable :: dmmowtb(:)   !! DM threshold vs DNR (flat pairs)

      ! Mowing regrowth delay table (Part 3 mowing settings)
      real(real64), allocatable :: dmmowdelay(:) !! DM harvest vs delay days (flat pairs)

      ! Grazing flexible DM threshold table (when swdmgrz=2)
      real(real64), allocatable :: dmgrztb(:)   !! DM threshold vs DNR (flat pairs)

      ! Grazing livestock density tables (schema 1:1; runtime skipped when seqgrazmow=all-2)
      real(real64), allocatable :: lsda(:)       !! Actual livestock density per period
      real(real64), allocatable :: daysgrazing(:)  !! Max days grazing per LSDb entry
      real(real64), allocatable :: uptgrazing(:)   !! DM uptake per LSDb entry
      real(real64), allocatable :: lossgrazing(:)  !! DM loss per LSDb entry

      ! Treading loss tables (schema 1:1; stub-errored at validator)
      integer      :: swlossmow = 0          !! 0=no losses (others stub)
      integer      :: swlossgrz = 0          !! 0=no losses (others stub)

      ! CO2 switch (Part 14)
      integer :: swco2 = 0                  !! 0=none (swco2=1 stub)

      ! Populated sentinel (ADR 0016)
      logical :: populated = .false.         !! true once read_cropgrass_toml completes
```

- [ ] **Step 5: Extend `cropgrass_config_validate` with stub-error guards**

Add the following block at the start of `cropgrass_config_validate`, before any
existing checks:

```fortran
      ! ----- Phase 3 stub-errors for unsupported runtime branches -----
      ! ADR 0015: schema accepts these values 1:1 with legacy; runtime
      ! plumbing for these branches has not been ported. Cases that
      ! need them must run via the legacy executable.
      if (self%swdrought == 2) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropgrass.swdrought=2 (De Jong van Lier) not yet supported in ' // &
            'the TOML pipeline; use the legacy executable.', 'cropgrass')
      end if
      if (self%swcompensate == 2) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropgrass.swcompensate=2 (Walsum) not yet supported in ' // &
            'the TOML pipeline; use the legacy executable.', 'cropgrass')
      end if
      if (self%swsalinity /= 0) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropgrass.swsalinity /= 0 not yet supported in the TOML ' // &
            'pipeline; use the legacy executable.', 'cropgrass')
      end if
      if (self%swinter == 2 .or. self%swinter == 3) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropgrass.swinter=2 or 3 (Gash/storage-cap) not yet supported ' // &
            'in the TOML pipeline; use the legacy executable.', 'cropgrass')
      end if
      if (self%swco2 == 1) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropgrass.swco2=1 (CO2 correction) not yet supported in the ' // &
            'TOML pipeline; use the legacy executable.', 'cropgrass')
      end if
      if (self%swlossmow == 1) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropgrass.swlossmow=1 (treading losses during mowing) not yet ' // &
            'supported in the TOML pipeline; use the legacy executable.', &
            'cropgrass')
      end if
      if (self%swlossgrz == 1) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropgrass.swlossgrz=1 (treading losses during grazing) not yet ' // &
            'supported in the TOML pipeline; use the legacy executable.', &
            'cropgrass')
      end if
      if (self%swoxygen == 2 .and. self%swoxygentype == 2) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'cropgrass.swoxygentype=2 (reproduction functions) not yet ' // &
            'supported in the TOML pipeline; use the legacy executable.', &
            'cropgrass')
      end if
      if (allocated(self%seqgrazmow)) then
         block
            integer :: _i
            do _i = 1, size(self%seqgrazmow)
               if (self%seqgrazmow(_i) /= 2) then
                  call errors%append(ERR_VALIDATION_CROSS_FIELD, &
                     'cropgrass.seqgrazmow: grazing (=1) and dewooling (=3) ' // &
                     'periods are not yet supported in the TOML pipeline; ' // &
                     'use the legacy executable.', 'cropgrass')
                  exit
               end if
            end do
         end block
      end if
```

Note: the `block` construct uses an identifier starting with `_` to avoid
shadowing any outer variable named `i`. Use a distinct name such as `jseq`
if the compiler does not accept leading underscores.

Also add `ERR_VALIDATION_CROSS_FIELD` to the module's `use error_mod` if not
already present.

- [ ] **Step 6: Run tests, verify they pass**

```bash
pixi run test-pfunit 2>&1 | grep -E "PASS|FAIL|Error" | tail -25
```

Expected: all 11 new tests pass (stub-error tests produce the expected error
codes; "supported values" tests produce zero stub-errors). Any pre-existing
failures are a baseline problem.

- [ ] **Step 7: Commit**

```bash
git add src/config/cropgrass_config.f90 tests/unit/config/test_cropgrass_config.pf
git commit -m "feat(cropgrass-config): schema 1:1 with readgrass + stub-error validators (Phase 3)"
```

---

## Task 3: Extend `read_cropgrass_toml.f90` parser

Add parsing for all new scalar and table fields. Add the `read_cropgrass_file_toml`
wrapper. Set `populated = .true.` at the end of `read_cropgrass_toml`.

**Files:**
- Modify: `src/io/toml/read_cropgrass_toml.f90`
- Test: `tests/unit/io/toml/test_read_cropgrass_toml.pf`

- [ ] **Step 1: Read the existing parser**

```bash
cat -n src/io/toml/read_cropgrass_toml.f90
```

Note the current structure: sections handled and the `read_array_1d` helper.

- [ ] **Step 2: Write a failing round-trip test for the new fields**

In `tests/unit/io/toml/test_read_cropgrass_toml.pf`, add a test that constructs
a TOML document (inline string) with all new fields and asserts that parsing
populates the correct config fields. Key assertions:

```fortran
@test
subroutine test_cropgrass_round_trip_new_scalars()
   use funit
   use tomlf, only: toml_parse, toml_table
   use read_cropgrass_toml_mod, only: read_cropgrass_toml
   use cropgrass_config_mod, only: cropgrass_config_t
   use error_mod, only: error_collection_t

   type(toml_table), pointer  :: doc
   type(cropgrass_config_t)   :: c
   type(error_collection_t)   :: errors
   character(len=:), allocatable :: toml_text
   integer :: stat

   toml_text = '[phenology]' // new_line('a') // &
               'swtsum = 1' // new_line('a') // &
               '[green_area]' // new_line('a') // &
               'ssa = 0.0004' // new_line('a') // &
               'span = 30.0' // new_line('a') // &
               '[assimilation]' // new_line('a') // &
               'cvl = 0.685' // new_line('a') // &
               'cvr = 0.694' // new_line('a') // &
               'cvs = 0.662' // new_line('a') // &
               '[oxygen_stress]' // new_line('a') // &
               'swoxygen = 2' // new_line('a') // &
               'swwrtnonox = 1' // new_line('a') // &
               'aeratecrit = 0.5' // new_line('a') // &
               '[oxygen_stress.bartholomeus]' // new_line('a') // &
               'swoxygentype = 1' // new_line('a') // &
               'srl = 383571.0' // new_line('a') // &
               'swrootradius = 2' // new_line('a') // &
               'root_radiuso2 = 0.000075' // new_line('a') // &
               '[drought_stress]' // new_line('a') // &
               'swdrought = 1' // new_line('a') // &
               '[compensation]' // new_line('a') // &
               'swcompensate = 1' // new_line('a') // &
               'alphacrit = 0.7' // new_line('a')

   call toml_parse(doc, toml_text, stat=stat)
   @assertEqual(0, stat)
   call read_cropgrass_toml(doc, c, errors)
   @assertEqual(0, errors%count())
   @assertEqual(1, c%swtsum)
   @assertRelativelyEqual(0.0004_real64, c%ssa, tolerance=1.0e-10_real64)
   @assertRelativelyEqual(30.0_real64, c%span, tolerance=1.0e-10_real64)
   @assertRelativelyEqual(0.685_real64, c%cvl, tolerance=1.0e-10_real64)
   @assertEqual(2, c%swoxygen)
   @assertEqual(1, c%swwrtnonox)
   @assertRelativelyEqual(0.5_real64, c%aeratecrit, tolerance=1.0e-10_real64)
   @assertRelativelyEqual(383571.0_real64, c%srl, tolerance=1.0e0_real64)
   @assertEqual(2, c%swrootradius)
   @assertRelativelyEqual(0.000075_real64, c%root_radiusO2, tolerance=1.0e-12_real64)
   @assertEqual(1, c%swdrought)
   @assertEqual(1, c%swcompensate)
   @assertRelativelyEqual(0.7_real64, c%alphacrit, tolerance=1.0e-10_real64)
   @assertTrue(c%populated)
end subroutine
```

- [ ] **Step 3: Run the failing test**

```bash
pixi run test-pfunit 2>&1 | grep -E "round_trip|FAIL" | tail -10
```

- [ ] **Step 4: Extend `read_cropgrass_toml.f90`**

Add new section-parsing blocks inside `read_cropgrass_toml` for:

**[phenology] extension:**
```fortran
call get_optional_int_with_default(ph, 'swtsum',    config%swtsum,    1, 'phenology.swtsum',    errors)
if (config%swtsum == 2) then
   call get_optional_real_with_default(ph, 'tsumtemp',  config%tsumtemp,  0.0_real64, 'phenology.tsumtemp',  errors)
   call get_optional_real_with_default(ph, 'tsumdepth', config%tsumdepth, 0.0_real64, 'phenology.tsumdepth', errors)
   call get_optional_int_with_default(ph, 'tsumtime',  config%tsumtime,  0,          'phenology.tsumtime',  errors)
end if
```

**[crop_state] section (new):**
```fortran
call get_table(doc, 'crop_state', cs, 'crop_state', errors)
if (associated(cs)) then
   call get_optional_real_with_default(cs, 'tdwi',   config%tdwi,   1000.0_real64, 'crop_state.tdwi',   errors)
   call get_optional_real_with_default(cs, 'laiem',  config%laiem,  0.63_real64,   'crop_state.laiem',  errors)
   call get_optional_real_with_default(cs, 'rgrlai', config%rgrlai, 0.007_real64,  'crop_state.rgrlai', errors)
end if
```

**[green_area] section (new):**
```fortran
call get_table(doc, 'green_area', ga, 'green_area', errors)
if (associated(ga)) then
   call get_optional_real_with_default(ga, 'ssa',  config%ssa,  0.0_real64,  'green_area.ssa',  errors)
   call get_optional_real_with_default(ga, 'span', config%span, 30.0_real64, 'green_area.span', errors)
   call read_array_1d(ga, 'slatb', config%slatb, 'green_area.slatb', errors)
end if
```

**[assimilation] section (new, for biomass conversion and respiration scalars):**
```fortran
call get_table(doc, 'assimilation', assim, 'assimilation', errors)
if (associated(assim)) then
   call get_optional_real_with_default(assim, 'cvl', config%cvl, 0.685_real64, 'assimilation.cvl', errors)
   call get_optional_real_with_default(assim, 'cvr', config%cvr, 0.694_real64, 'assimilation.cvr', errors)
   call get_optional_real_with_default(assim, 'cvs', config%cvs, 0.662_real64, 'assimilation.cvs', errors)
   call get_optional_real_with_default(assim, 'q10', config%q10, 2.0_real64,   'assimilation.q10', errors)
   call get_optional_real_with_default(assim, 'rml', config%rml, 0.03_real64,  'assimilation.rml', errors)
   call get_optional_real_with_default(assim, 'rmr', config%rmr, 0.015_real64, 'assimilation.rmr', errors)
   call get_optional_real_with_default(assim, 'rms', config%rms, 0.015_real64, 'assimilation.rms', errors)
   call get_optional_real_with_default(assim, 'perdl', config%perdl, 0.05_real64, 'assimilation.perdl', errors)
   call read_array_1d(assim, 'amaxtb',  config%amaxtb,  'assimilation.amaxtb',  errors)
   call read_array_1d(assim, 'tmpftb',  config%tmpftb,  'assimilation.tmpftb',  errors)
   call read_array_1d(assim, 'tmnftb',  config%tmnftb,  'assimilation.tmnftb',  errors)
   call read_array_1d(assim, 'rfsetb',  config%rfsetb,  'assimilation.rfsetb',  errors)
   call read_array_1d(assim, 'frtb',    config%frtb,    'assimilation.frtb',    errors)
   call read_array_1d(assim, 'fltb',    config%fltb,    'assimilation.fltb',    errors)
   call read_array_1d(assim, 'fstb',    config%fstb,    'assimilation.fstb',    errors)
   call read_array_1d(assim, 'rdrrtb',  config%rdrrtb,  'assimilation.rdrrtb',  errors)
   call read_array_1d(assim, 'rdrstb',  config%rdrstb,  'assimilation.rdrstb',  errors)
end if
```

**[crop_factor] section (new):**
```fortran
call get_table(doc, 'crop_factor', cf_sec, 'crop_factor', errors)
if (associated(cf_sec)) then
   call get_optional_int_with_default(cf_sec, 'swcf',   config%swcf,   2,           'crop_factor.swcf',   errors)
   call get_optional_real_with_default(cf_sec, 'albedo', config%albedo, 0.23_real64, 'crop_factor.albedo', errors)
   call get_optional_real_with_default(cf_sec, 'rsw',    config%rsw,    0.0_real64,  'crop_factor.rsw',    errors)
   ! cftb and chtb already allocated by existing parser (check not to double-read)
end if
```

**[root] extension** (add swrd, swdmi2rd, swrdc, wrtmax, rdtb, rlwtb):
```fortran
call get_optional_int_with_default(root,  'swrd',     config%swrd,     2,          'root.swrd',     errors)
call get_optional_int_with_default(root,  'swdmi2rd', config%swdmi2rd, 0,          'root.swdmi2rd', errors)
call get_optional_int_with_default(root,  'swrdc',    config%swrdc,    0,          'root.swrdc',    errors)
call get_optional_real_with_default(root, 'wrtmax',   config%wrtmax,   3000.0_real64, 'root.wrtmax', errors)
call read_array_1d(root, 'rdtb',  config%rdtb,  'root.rdtb',  errors)
call read_array_1d(root, 'rlwtb', config%rlwtb, 'root.rlwtb', errors)
```

**[oxygen_stress] section (new):**
```fortran
call get_table(doc, 'oxygen_stress', oxy, 'oxygen_stress', errors)
if (associated(oxy)) then
   call get_optional_int_with_default(oxy, 'swoxygen',    config%swoxygen,    1, 'oxygen_stress.swoxygen',    errors)
   call get_optional_int_with_default(oxy, 'swwrtnonox',  config%swwrtnonox,  0, 'oxygen_stress.swwrtnonox',  errors)
   call get_optional_real_with_default(oxy, 'aeratecrit', config%aeratecrit, 1.0e-4_real64, 'oxygen_stress.aeratecrit', errors)

   call get_table(oxy, 'bartholomeus', barto, 'oxygen_stress.bartholomeus', errors)
   if (associated(barto)) then
      call get_optional_int_with_default(barto,  'swoxygentype',           config%swoxygentype,          1,               'oxygen_stress.bartholomeus.swoxygentype',           errors)
      call get_optional_real_with_default(barto, 'q10_microbial',          config%q10_microbial,         2.8_real64,      'oxygen_stress.bartholomeus.q10_microbial',          errors)
      call get_optional_real_with_default(barto, 'specific_resp_humus',    config%specific_resp_humus,   1.6e-3_real64,   'oxygen_stress.bartholomeus.specific_resp_humus',    errors)
      call get_optional_real_with_default(barto, 'srl',                    config%srl,                   383571.0_real64, 'oxygen_stress.bartholomeus.srl',                    errors)
      call get_optional_int_with_default(barto,  'swrootradius',           config%swrootradius,          2,               'oxygen_stress.bartholomeus.swrootradius',           errors)
      call get_optional_real_with_default(barto, 'dry_mat_cont_roots',     config%dry_mat_cont_roots,    0.075_real64,    'oxygen_stress.bartholomeus.dry_mat_cont_roots',     errors)
      call get_optional_real_with_default(barto, 'air_filled_root_por',    config%air_filled_root_por,   0.05_real64,     'oxygen_stress.bartholomeus.air_filled_root_por',    errors)
      call get_optional_real_with_default(barto, 'spec_weight_root_tissue',config%spec_weight_root_tissue,1.0e3_real64,   'oxygen_stress.bartholomeus.spec_weight_root_tissue',errors)
      call get_optional_real_with_default(barto, 'var_a',                  config%var_a,                 4.175e-10_real64,'oxygen_stress.bartholomeus.var_a',                  errors)
      call get_optional_real_with_default(barto, 'root_radiuso2',          config%root_radiusO2,         0.000075_real64, 'oxygen_stress.bartholomeus.root_radiuso2',          errors)
   end if
end if
```

**[drought_stress] section (new):**
```fortran
call get_table(doc, 'drought_stress', dro, 'drought_stress', errors)
if (associated(dro)) then
   call get_optional_int_with_default(dro, 'swdrought', config%swdrought, 1, 'drought_stress.swdrought', errors)
end if
```

**[compensation] section (new):**
```fortran
call get_table(doc, 'compensation', comp, 'compensation', errors)
if (associated(comp)) then
   call get_optional_int_with_default(comp,  'swcompensate', config%swcompensate, 0,          'compensation.swcompensate', errors)
   call get_optional_int_with_default(comp,  'swstressor',   config%swstressor,   1,          'compensation.swstressor',   errors)
   call get_optional_real_with_default(comp, 'alphacrit',    config%alphacrit,    1.0_real64, 'compensation.alphacrit',    errors)
   call get_optional_real_with_default(comp, 'dcritrtz',     config%dcritrtz,     0.0_real64, 'compensation.dcritrtz',     errors)
end if
```

**[management] section (new):**
```fortran
call get_table(doc, 'management', mgmt, 'management', errors)
if (associated(mgmt)) then
   call get_optional_real_with_default(mgmt, 'mowrest',   config%mowrest,   700.0_real64, 'management.mowrest',   errors)
   call get_optional_real_with_default(mgmt, 'dewrest',   config%dewrest,   850.0_real64, 'management.dewrest',   errors)
   call get_optional_int_with_default(mgmt,  'swpotrelmf',config%swpotrelmf,1,            'management.swpotrelmf',errors)
   call get_optional_real_with_default(mgmt, 'relmf',     config%relmf,     1.0_real64,   'management.relmf',     errors)
   call read_array_1d_int(mgmt, 'seqgrazmow', config%seqgrazmow, config%nseqgrazmow, &
                          'management.seqgrazmow', errors)
   call read_array_1d(mgmt, 'dmmowtb',     config%dmmowtb,     'management.dmmowtb',     errors)
   call read_array_1d(mgmt, 'dmmowdelay',  config%dmmowdelay,  'management.dmmowdelay',  errors)
   call read_array_1d(mgmt, 'dmgrztb',     config%dmgrztb,     'management.dmgrztb',     errors)
   call read_array_1d(mgmt, 'lsda',        config%lsda,        'management.lsda',        errors)
   call read_array_1d(mgmt, 'daysgrazing', config%daysgrazing, 'management.daysgrazing', errors)
   call read_array_1d(mgmt, 'uptgrazing',  config%uptgrazing,  'management.uptgrazing',  errors)
   call read_array_1d(mgmt, 'lossgrazing', config%lossgrazing, 'management.lossgrazing', errors)
   ! swlossmow and swlossgrz (default 0):
   call get_optional_int_with_default(mgmt, 'swlossmow', config%swlossmow, 0, 'management.swlossmow', errors)
   call get_optional_int_with_default(mgmt, 'swlossgrz', config%swlossgrz, 0, 'management.swlossgrz', errors)
end if
```

**[co2] section (new):**
```fortran
call get_table(doc, 'co2', co2_sec, 'co2', errors)
if (associated(co2_sec)) then
   call get_optional_int_with_default(co2_sec, 'swco2', config%swco2, 0, 'co2.swco2', errors)
end if
```

**At end of subroutine, set populated:**
```fortran
config%populated = .true.
```

**Add `read_array_1d_int` private helper** for integer arrays (needed by
`seqgrazmow`). Pattern mirrors `read_array_1d` but uses `integer` type.

**Add `read_cropgrass_file_toml` wrapper:**
```fortran
subroutine read_cropgrass_file_toml(path, config, errors, base_path)
   character(len=*),         intent(in)    :: path
   type(cropgrass_config_t), intent(inout) :: config
   type(error_collection_t), intent(inout) :: errors
   character(len=*),         intent(in)    :: base_path

   type(toml_table), pointer :: doc
   character(len=:), allocatable :: full_path
   integer :: stat

   full_path = trim(base_path) // '/' // trim(path)
   call toml_load(doc, full_path, stat=stat)
   if (stat /= 0 .or. .not. associated(doc)) then
      call errors%append(ERR_PARSE_IO, &
                         'failed to open ' // trim(full_path), &
                         'crop.rotation.file')
      return
   end if
   call read_cropgrass_toml(doc, config, errors)
end subroutine
```

Export `read_cropgrass_file_toml` in the module's `public` statement.

- [ ] **Step 5: Run tests, verify they pass**

```bash
pixi run test-pfunit 2>&1 | grep -E "PASS|FAIL|Error" | tail -30
```

Expected: all new round-trip tests pass and `populated=.true.` assertion passes.

- [ ] **Step 6: Commit**

```bash
git add src/io/toml/read_cropgrass_toml.f90 tests/unit/io/toml/test_read_cropgrass_toml.pf
git commit -m "feat(read-cropgrass-toml): extend parser 1:1 with readgrass + file wrapper (Phase 3)"
```

---

## Task 4: Verify loader dispatch for type=3

Confirm that `read_crop_toml.f90` already allocates `rotation_grass(:)` and
calls `read_cropgrass_file_toml` for type=3 entries. If the dispatch exists but
still calls an old signature, update it. If it does not exist, add it.

**Files:**
- Modify if needed: `src/io/toml/read_crop_toml.f90`
- Test: `tests/unit/io/toml/test_load_swap_config.pf`

- [ ] **Step 1: Read the existing loader**

```bash
cat -n src/io/toml/read_crop_toml.f90
```

Look for the rotation loop, the `select case (rotation_type(i))` or equivalent,
and the type=3 branch. Record exact line numbers.

- [ ] **Step 2: Check allocation of `rotation_grass`**

Confirm that `crop_config_t%rotation_grass(:)` (or whatever the actual name is)
is allocated to `nrot` before the loop. If not, add the allocation immediately
after the `nrot` determination:

```fortran
allocate(config%rotation_grass(nrot))  ! name confirmed in pre-flight
```

- [ ] **Step 3: Check type=3 dispatch**

Confirm the type=3 case calls `read_cropgrass_file_toml`. If it calls
the old in-place reader directly (without resolving the file), update it:

```fortran
case (3)
   call read_cropgrass_file_toml(rotation_file(i), config%rotation_grass(i), &
                                 errors, base_path=case_dir)
```

If the dispatch uses a different convention for `base_path`, match it.

- [ ] **Step 4: Write failing tests**

In `tests/unit/io/toml/test_load_swap_config.pf`, add:

```fortran
@test
subroutine test_swap_config_case4_rotation_grass_loaded()
   ! Load case 4's swap.toml, assert 10 rotation_grass slots populated.
   use funit
   use swap_config_mod, only: swap_config_t, read_swap_config
   use error_mod, only: error_collection_t
   type(swap_config_t)       :: cfg
   type(error_collection_t)  :: errors
   integer :: i

   call read_swap_config('tests/swap-cases/toml/4.oxygenstress/swap.toml', cfg, errors)
   @assertEqual(0, errors%count(), message='Loading case 4 swap.toml produced errors')
   @assertTrue(allocated(cfg%crop%rotation_grass))
   @assertEqual(10, size(cfg%crop%rotation_grass))
   do i = 1, 10
      @assertTrue(cfg%crop%rotation_grass(i)%populated, &
                  message='rotation_grass slot not populated')
   end do
end subroutine

@test
subroutine test_swap_config_case2_rotation_grass_loaded()
   ! Load case 2's swap.toml, assert all rotation_grass slots populated.
   use funit
   use swap_config_mod, only: swap_config_t, read_swap_config
   use error_mod, only: error_collection_t
   type(swap_config_t)       :: cfg
   type(error_collection_t)  :: errors
   integer :: i
   integer :: nrot

   call read_swap_config('tests/swap-cases/toml/2.grassgrowth/swap.toml', cfg, errors)
   @assertEqual(0, errors%count(), message='Loading case 2 swap.toml produced errors')
   @assertTrue(allocated(cfg%crop%rotation_grass))
   nrot = size(cfg%crop%rotation_grass)
   @assertTrue(nrot > 0)
   do i = 1, nrot
      @assertTrue(cfg%crop%rotation_grass(i)%populated, &
                  message='rotation_grass slot not populated')
   end do
end subroutine
```

Note: these tests depend on Tasks 5 and 6 (the TOML files being authored).
They will initially fail because the TOML files are still skeletal. Document
this dependency — the tests become green after Task 6.

- [ ] **Step 5: Run tests (expected: fail until Task 6)**

```bash
pixi run test-pfunit 2>&1 | grep -E "case4_rotation|case2_rotation|FAIL" | head -10
```

- [ ] **Step 6: Commit loader changes (if any)**

```bash
git add src/io/toml/read_crop_toml.f90
git commit -m "feat(read-crop-toml): allocate rotation_grass and dispatch type=3 file read (Phase 3)"
```

If no changes were needed (loader already correct), skip this commit.

---

## Task 5: Author `4.oxygenstress/grassd.crp.toml` (submodule pair commit)

Replace the 53-line skeletal TOML with a full 1:1 reproduction of the legacy
`grassd.crp` for case 4.

**Files:**
- Submodule modify: `tests/swap-cases/toml/4.oxygenstress/grassd.crp.toml`

- [ ] **Step 1: Read the legacy .crp file end-to-end**

```bash
cat tests/swap-cases/4.oxygenstress/grassd.crp
```

Map every field to its TOML section as designed in Unit 6 of the spec.

- [ ] **Step 2: Read the existing TOML stub**

```bash
cat tests/swap-cases/toml/4.oxygenstress/grassd.crp.toml
```

Note what is already present (light, root, water_stress, interception, mowing,
grazing sections).

- [ ] **Step 3: Write the full TOML file**

Replace the stub with a complete file. Structure:

```toml
# oxygenstress grass crop (type 3) — Phase 3 full port
# 1:1 reproduction of tests/swap-cases/4.oxygenstress/grassd.crp

[crop_state]
tdwi   = 1000.0
laiem  = 0.63
rgrlai = 0.007

[phenology]
swtsum = 1

[crop_factor]
swcf   = 2
albedo = 0.23
rsw    = 0.0

# chtb: DNR, CH pairs (swcf=2 crop height table)
chtb = [0.0, 12.0, 180.0, 12.0, 366.0, 12.0]

[green_area]
ssa  = 0.0004
span = 30.0
slatb = [1.0, 0.0015, 80.0, 0.0015, 300.0, 0.0020, 366.0, 0.0020]

[assimilation]
kdif  = 0.60         # (also in [light] for backward compat)
kdir  = 0.75
eff   = 0.50
amaxtb  = [1.0, 40.0, 95.0, 40.0, 200.0, 35.0, 275.0, 25.0, 366.0, 25.0]
tmpftb  = [0.0, 0.0, 5.0, 0.70, 15.0, 1.0, 25.0, 1.0, 40.0, 0.0]
tmnftb  = [0.0, 0.0, 4.0, 1.0]
cvl   = 0.685
cvr   = 0.694
cvs   = 0.662
q10   = 2.0
rml   = 0.030
rmr   = 0.015
rms   = 0.015
rfsetb  = [1.0, 1.0, 366.0, 1.0]
frtb    = [1.0, 0.3, 366.0, 0.3]
fltb    = [1.0, 0.6, 366.0, 0.6]
fstb    = [1.0, 0.4, 366.0, 0.4]
perdl   = 0.05
rdrrtb  = [1.0, 0.0, 180.0, 0.02, 366.0, 0.02]
rdrstb  = [1.0, 0.0, 180.0, 0.02, 366.0, 0.02]

[root]
rdi      = 10.0
rri      =  1.0
rdc      = 40.0
swrd     = 2
swdmi2rd = 1
swrdc    = 0
rdtb     = [1.0, 10.0, 180.0, 40.0, 366.0, 40.0]
rlwtb    = [300.0, 10.0, 2500.0, 40.0]
wrtmax   = 3000.0
rdctb    = [0.0, 1.0, 1.0, 0.0]

[oxygen_stress]
swoxygen   = 2
swwrtnonox = 1
aeratecrit = 0.5

[oxygen_stress.bartholomeus]
swoxygentype          = 1
q10_microbial         = 2.8
specific_resp_humus   = 1.6e-3
srl                   = 383571.0
swrootradius          = 2
root_radiuso2         = 0.000075

[drought_stress]
swdrought = 1

[water_stress]
hlim1  =     0.0
hlim2u =     1.0
hlim2l =    -1.0
hlim3h =  -200.0
hlim3l =  -800.0
hlim4  = -8000.0
adcrh  =     0.5
adcrl  =     0.1

[compensation]
swcompensate = 1
swstressor   = 1
alphacrit    = 0.7

[interception]
swinter = 1
cofab   = 0.25

[management]
seqgrazmow  = [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]
mowrest     = 700.0
swpotrelmf  = 1
relmf       = 0.90
# swdmmow=2: flexible DM threshold table for mowing
dmmowtb     = [120.0, 4700.0, 152.0, 3700.0, 182.0, 3200.0, 213.0, 2700.0, 366.0, 2700.0]
# dmmowdelay: DM harvest value, days of delay pairs
dmmowdelay  = [0.0, 2.0, 2000.0, 3.0, 4000.0, 4.0]
# swdmgrz=2: flexible DM threshold table for grazing
dmgrztb     = [152.0, 2400.0, 244.0, 1800.0, 366.0, 1800.0]
tagprest    = 700.0
lsda        = [21.25, 21.25, 21.25, 21.25, 21.25, 21.25, 21.25, 21.25, 21.25,
               21.25, 21.25, 21.25, 21.25, 21.25, 21.25, 21.25, 21.25, 21.25,
               21.25, 21.25]
lsdb        = [21.25]
daysgrazing = [4.0]
uptgrazing  = [16.0]
lossgrazing = [4.0]

[mowing]
swharv         = 1
nmow           = 1
swdmmow        = 2
dmharvest      = 4200.0
daylastharvest = 289.0
dmlastharvest  = 2700.0
maxdaymow      = 42

[grazing]
swgraz      = 1
nstart_graz = 152
nstop_graz  = 366
maxdaygrz   = 28
dmgrazing   = 2400.0
swdmgrz     = 2

[co2]
swco2 = 0

[irrigation_schedule]
schedule = 0
```

Note: `kdif`, `kdir`, `eff`, `amax` were on the old `[light]` section. Move
them into `[assimilation]` for Phase 3 (the parser now reads them from there).
Keep `[light]` section for backward compat if needed, or remove it and update
the parser accordingly. Confirm parser handles the consolidation correctly.

- [ ] **Step 4: Run regression for case 4 and case 2**

```bash
pixi run regression 2>&1 | tail -10
```

Expected: 5/5 green (the legacy reader still runs because the TOML path is
not yet wired — the new TOML content has no effect yet).

- [ ] **Step 5: Submodule inner commit + outer bump**

```bash
cd tests/swap-cases
git add toml/4.oxygenstress/grassd.crp.toml
git commit toml/4.oxygenstress/grassd.crp.toml \
    -m "feat(case4): grassd.crp.toml 1:1 reproduction for Phase 3 cropgrass port"
cd /home/zawadzkim/Code/swap
git add tests/swap-cases
git commit -m "chore(submodule): bump swap-cases for case 4 grassd.crp.toml authoring (Phase 3)"
```

---

## Task 6: Author `2.grassgrowth/grassd.crp.toml` (submodule pair commit)

Replace the 59-line skeletal TOML with a full 1:1 reproduction of case 2's
legacy `grassd.crp`.

**Files:**
- Submodule modify: `tests/swap-cases/toml/2.grassgrowth/grassd.crp.toml`

- [ ] **Step 1: Read the legacy .crp file**

```bash
cat tests/swap-cases/2.grassgrowth/grassd.crp
```

Key differences from case 4: `swoxygen=1`, `swcompensate=0`, `swharvest=2`
with a `dateharvest` table (30 entries, calendar dates 1980-1984 → already
converted to DOY floats in the existing TOML stub).

- [ ] **Step 2: Write the full TOML file**

Same structure as case 4 except:

```toml
[oxygen_stress]
swoxygen   = 1          # Feddes (no bartholomeus sub-section)
swwrtnonox = 1
aeratecrit = 0.5

[water_stress]
hlim1  =     0.0        # Feddes hlim1/hlim2u/hlim2l used when swoxygen=1
hlim2u =     1.0
hlim2l =    -1.0
hlim3h =  -200.0
...

[compensation]
swcompensate = 0        # No compensation

[mowing]
swharv         = 2      # Fixed-date mowing
nmow           = 30
swdmmow        = 2
...
mowing_dates = [
  127.0, 149.0, 176.0, 206.0, 232.0, 261.0, 297.0,
  104.0, 139.0, 167.0, 195.0, 217.0, 251.0, 301.0,
  131.0, 152.0, 187.0, 222.0, 286.0,
  139.0, 166.0, 194.0, 229.0, 301.0,
  136.0, 159.0, 187.0, 215.0, 255.0, 312.0,
]
```

No `[oxygen_stress.bartholomeus]` sub-section (swoxygen=1 does not need it).
No `swcompensate=1` fields.

- [ ] **Step 3: Run regression for both cases**

```bash
pixi run regression 2>&1 | tail -10
```

Expected: 5/5 green (legacy reader still running).

- [ ] **Step 4: Submodule inner commit + outer bump**

```bash
cd tests/swap-cases
git add toml/2.grassgrowth/grassd.crp.toml
git commit toml/2.grassgrowth/grassd.crp.toml \
    -m "feat(case2): grassd.crp.toml 1:1 reproduction for Phase 3 cropgrass port"
cd /home/zawadzkim/Code/swap
git add tests/swap-cases
git commit -m "chore(submodule): bump swap-cases for case 2 grassd.crp.toml authoring (Phase 3)"
```

After this task, re-run the Task 4 failing tests:

```bash
pixi run test-pfunit 2>&1 | grep -E "case4_rotation|case2_rotation" | head -10
```

Expected: both rotation-loaded tests now pass (the TOML files are authored and
the loader can parse them).

---

## Task 7: Create `src/crop/cropgrass_init.f90`

Implement `cropgrass_init_from_config`: config → module globals copy + runtime
init math + defense-in-depth guards.

**Files:**
- Create: `src/crop/cropgrass_init.f90`
- Modify: `meson.build` (register new file)
- Test: `tests/unit/crop/test_cropgrass_init.pf` (new directory + file)

- [ ] **Step 1: Read the `grass(task=1)` runtime code**

```bash
sed -n '2063,2160p' src/crop/cropgrowth.f90 | cat
```

Understand what globals are consumed immediately after the `readgrass` call
(the initialization block uses `rdm`, `dvs`, `rid`, `frtb`, `fltb`, `fstb`,
`slatb`, `wrt`, `wst`, `wlv`, `lv`, etc.). The `cropgrass_init_from_config`
subroutine must write those globals so the subsequent code is unaffected.

- [ ] **Step 2: Read the `readgrass` runtime tail**

```bash
sed -n '4085,4125p' src/io/readswap.f90 | cat
```

The runtime tail is:
1. `close(crp)` — skip in init sub.
2. `irrigation(1)` call when `schedule==1` — stub-errored.
3. `cumdens` trapezium build (lines 4091-4121) — port verbatim.
4. `swinco==3` .END-file read (lines 4123-4269) — GUARDED bucket.

- [ ] **Step 3: Write failing unit tests**

Create `tests/unit/crop/test_cropgrass_init.pf`:

```fortran
module test_cropgrass_init_mod
   use funit
   use cropgrass_config_mod, only: cropgrass_config_t
   use cropgrass_init_mod,   only: cropgrass_init_from_config
   use variables,            only: kdif, kdir, swoxygen, srl, swcompensate, &
                                   alphacrit, swharvest, dateharvest, cumdens
   implicit none

contains

@test
subroutine test_cropgrass_init_writes_kdif()
   type(cropgrass_config_t) :: cfg
   cfg%kdif = 0.60_real64
   cfg%kdir = 0.75_real64
   ! ... minimal required fields to avoid fatal errors ...
   cfg%swoxygen    = 0
   cfg%swdrought   = 1
   cfg%swcompensate = 0
   cfg%swco2       = 0
   cfg%swlossgrz   = 0
   cfg%swlossmow   = 0
   cfg%nseqgrazmow = 1
   allocate(cfg%seqgrazmow(1))
   cfg%seqgrazmow  = 2
   call cropgrass_init_from_config(cfg, icrop=1)
   @assertRelativelyEqual(0.60_real64, kdif, tolerance=1.0e-10_real64)
   @assertRelativelyEqual(0.75_real64, kdir, tolerance=1.0e-10_real64)
end subroutine

@test
subroutine test_cropgrass_init_cumdens_case4()
   ! rdctb = [(0.0, 1.0), (1.0, 0.0)] — linear decline.
   ! Trapezium sum over 100 intervals of width 0.01:
   ! integral = 0.5 (area of triangle). cumdens(2) = 0/0.5 = 0.
   ! cumdens(202) = 0.5/0.5 = 1.0 after normalization.
   type(cropgrass_config_t) :: cfg
   cfg%swdrought   = 1
   cfg%swoxygen    = 0
   cfg%swcompensate = 0
   cfg%swco2       = 0
   cfg%swlossgrz   = 0
   cfg%swlossmow   = 0
   cfg%nseqgrazmow = 1
   allocate(cfg%seqgrazmow(1))
   cfg%seqgrazmow  = 2
   allocate(cfg%rdctb(4))
   cfg%rdctb = [0.0_real64, 1.0_real64, 1.0_real64, 0.0_real64]
   call cropgrass_init_from_config(cfg, icrop=1)
   @assertRelativelyEqual(0.0_real64,  cumdens(2),   tolerance=1.0e-6_real64)
   @assertRelativelyEqual(1.0_real64,  cumdens(202), tolerance=1.0e-6_real64)
end subroutine

@test
subroutine test_cropgrass_init_swoxygen2_writes_srl()
   type(cropgrass_config_t) :: cfg
   cfg%swoxygen      = 2
   cfg%swoxygentype  = 1
   cfg%srl           = 383571.0_real64
   cfg%swrootradius  = 2
   cfg%root_radiusO2 = 0.000075_real64
   cfg%swwrtnonox    = 0
   cfg%swdrought     = 1
   cfg%swcompensate  = 0
   cfg%swco2         = 0
   cfg%swlossgrz     = 0
   cfg%swlossmow     = 0
   cfg%nseqgrazmow   = 1
   allocate(cfg%seqgrazmow(1))
   cfg%seqgrazmow    = 2
   call cropgrass_init_from_config(cfg, icrop=1)
   @assertRelativelyEqual(383571.0_real64, srl, tolerance=1.0_real64)
end subroutine

@test
subroutine test_cropgrass_init_swcompensate1_writes_alphacrit()
   type(cropgrass_config_t) :: cfg
   cfg%swcompensate = 1
   cfg%alphacrit    = 0.7_real64
   cfg%swstressor   = 1
   cfg%swoxygen     = 0
   cfg%swdrought    = 1
   cfg%swco2        = 0
   cfg%swlossgrz    = 0
   cfg%swlossmow    = 0
   cfg%nseqgrazmow  = 1
   allocate(cfg%seqgrazmow(1))
   cfg%seqgrazmow   = 2
   call cropgrass_init_from_config(cfg, icrop=1)
   @assertRelativelyEqual(0.7_real64, alphacrit, tolerance=1.0e-10_real64)
end subroutine

@test
subroutine test_cropgrass_init_swharvest2_populates_dateharvest()
   type(cropgrass_config_t) :: cfg
   cfg%swharv       = 2
   cfg%nmow         = 3
   allocate(cfg%mowing_dates(3))
   cfg%mowing_dates = [127.0_real64, 149.0_real64, 176.0_real64]
   cfg%swoxygen     = 0
   cfg%swdrought    = 1
   cfg%swcompensate = 0
   cfg%swco2        = 0
   cfg%swlossgrz    = 0
   cfg%swlossmow    = 0
   cfg%nseqgrazmow  = 1
   allocate(cfg%seqgrazmow(1))
   cfg%seqgrazmow   = 2
   call cropgrass_init_from_config(cfg, icrop=1)
   ! dateharvest is stored as real(8) t1900-based timestamps in legacy;
   ! the init sub converts DOY floats accordingly. Assert the count is correct.
   ! (The exact conversion formula must be inspected from readgrass lines 3941-3943.)
   ! At minimum, assert dateharvest(1) is populated (non-zero).
   @assertTrue(dateharvest(1) > 0.0_real64)
end subroutine

end module test_cropgrass_init_mod
```

Register `tests/unit/crop/` in `tests/unit/meson.build` and add the suite to
`tests/unit/testSuites.inc` (follow Phase 1's `test_cropfixed_init.pf`
registration as the template).

- [ ] **Step 4: Run tests, verify failure**

```bash
pixi run test-pfunit 2>&1 | grep -E "cropgrass_init|FAIL" | head -10
```

Expected: all 5 new tests fail (module does not exist yet).

- [ ] **Step 5: Implement `cropgrass_init.f90`**

Create `src/crop/cropgrass_init.f90`. Module structure:

```fortran
module cropgrass_init_mod
   use iso_fortran_env, only: real64
   use cropgrass_config_mod, only: cropgrass_config_t
   use variables  ! (large list — use only: clause matching readgrass's use statement)
   use error_mod, only: fatalerr_collected
   use array_utils, only: afgen
   implicit none
   private
   public :: cropgrass_init_from_config

contains

   subroutine cropgrass_init_from_config(cfg, icrop)
      class(cropgrass_config_t), intent(in) :: cfg
      integer,                   intent(in) :: icrop

      ! --- defense-in-depth runtime guards ---
      if (cfg%swdrought == 2 .or. cfg%swsalinity /= 0 .or. &
          cfg%swcompensate == 2 .or. cfg%swco2 == 1 .or.    &
          cfg%swlossgrz == 1 .or. cfg%swlossmow == 1) then
         call fatalerr_collected('cropgrass_init', &
            'Unsupported runtime branch reached on the TOML path. ' // &
            'The validator should have caught this earlier.')
      end if
      if (cfg%swoxygen == 2 .and. cfg%swoxygentype == 2) then
         call fatalerr_collected('cropgrass_init', &
            'swoxygen=2 swoxygentype=2 (reproduction functions) not yet ' // &
            'supported on the TOML path.')
      end if

      ! --- Config → module globals copy ---
      ! (mirror readgrass's rd* calls 1:1, grouped by section)

      ! ET-related
      swcf   = cfg%swcf
      albedo = cfg%albedo
      rsw    = cfg%rsw
      if (cfg%swcf == 2) then
         ! chtb: flat (DNR, CH) pairs in cfg → chtb module global
         ! (copy logic mirrors readgrass lines 3529-3536)
         call copy_table(cfg%chtb, chtb)
         cftb = -99.99d0
      else if (cfg%swcf == 1) then
         call copy_table(cfg%cftb, cftb)
         chtb = -99.99d0
      end if

      ! Interception
      swinter = cfg%swinter
      if (cfg%swinter == 1) then
         cofab = cfg%cofab
      end if

      ! Crop state
      tdwi   = cfg%tdwi
      laiem  = cfg%laiem
      rgrlai = cfg%rgrlai

      ! Start of growth
      swtsum = cfg%swtsum
      if (cfg%swtsum == 2) then
         tsumtemp  = cfg%tsumtemp
         tsumtime  = cfg%tsumtime
         tsumdepth = cfg%tsumdepth
      end if

      ! Green area
      call copy_table(cfg%slatb, slatb)
      ssa    = cfg%ssa
      span   = cfg%span
      tbase  = cfg%tbase

      ! Assimilation
      kdif = cfg%kdif
      kdir = cfg%kdir
      eff  = cfg%eff
      call copy_table(cfg%amaxtb,  amaxtb)
      call copy_table(cfg%tmpftb,  tmpftb)
      call copy_table(cfg%tmnftb,  tmnftb)

      ! Conversion
      cvl = cfg%cvl
      cvr = cfg%cvr
      cvs = cfg%cvs

      ! Respiration
      q10 = cfg%q10
      rml = cfg%rml
      rmr = cfg%rmr
      rms = cfg%rms
      call copy_table(cfg%rfsetb, rfsetb)

      ! Partitioning
      call copy_table(cfg%frtb, frtb)
      call copy_table(cfg%fltb, fltb)
      call copy_table(cfg%fstb, fstb)

      ! Death rates
      perdl = cfg%perdl
      call copy_table(cfg%rdrrtb, rdrrtb)
      call copy_table(cfg%rdrstb, rdrstb)

      ! Root (swrd=2 path for both cases)
      swrdc = cfg%swrdc
      call copy_table(cfg%rdctb, rdctb)
      swrd = cfg%swrd
      if (cfg%swrd == 1) then
         call copy_table(cfg%rdtb, rdtb)
      else if (cfg%swrd == 2) then
         rdi      = cfg%rdi
         rri      = cfg%rri
         rdc      = cfg%rdc
         swdmi2rd = cfg%swdmi2rd
      else if (cfg%swrd == 3) then
         call copy_table(cfg%rlwtb, rlwtb)
         wrtmax = cfg%wrtmax
      end if

      ! Oxygen stress
      swoxygen = cfg%swoxygen
      if (cfg%swoxygen == 1) then
         hlim1  = cfg%hlim1
         hlim2u = cfg%hlim2u
         hlim2l = cfg%hlim2l
      else if (cfg%swoxygen == 2) then
         ! swoxygentype=1 physical path
         q10_microbial          = cfg%q10_microbial
         specific_resp_humus    = cfg%specific_resp_humus
         srl                    = cfg%srl
         swrootradius           = cfg%swrootradius
         if (cfg%swrootradius == 1) then
            dry_mat_cont_roots    = cfg%dry_mat_cont_roots
            air_filled_root_por   = cfg%air_filled_root_por
            spec_weight_root_tissue = cfg%spec_weight_root_tissue
            var_a                 = cfg%var_a
         else
            root_radiusO2 = cfg%root_radiusO2
         end if
      end if
      swwrtnonox = cfg%swwrtnonox
      aeratecrit = cfg%aeratecrit

      ! Drought stress
      swdrought = cfg%swdrought
      if (cfg%swdrought == 1) then
         hlim3h = cfg%hlim3h
         hlim3l = cfg%hlim3l
         hlim4  = cfg%hlim4
         adcrh  = cfg%adcrh
         adcrl  = cfg%adcrl
      end if

      ! Compensation
      swcompensate = cfg%swcompensate
      if (cfg%swcompensate == 1) then
         alphacrit  = cfg%alphacrit
         swstressor = cfg%swstressor
      end if

      ! Management
      mowrest    = cfg%mowrest
      swpotrelmf = cfg%swpotrelmf
      relmf      = cfg%relmf
      if (allocated(cfg%seqgrazmow)) then
         seqgrazmow(1:cfg%nseqgrazmow) = cfg%seqgrazmow
      end if

      ! Mowing (swharvest is set by the mowing block readgrass read)
      swharvest = cfg%swharv   ! note: cropgrass_config_t uses 'swharv'
      if (cfg%swharv == 1) then
         swdmmow = cfg%swdmmow
         if (cfg%swdmmow == 1) then
            dmharvest      = cfg%dmharvest
            daylastharvest = int(cfg%daylastharvest)
            dmlastharvest  = cfg%dmlastharvest
         else if (cfg%swdmmow == 2) then
            call copy_table(cfg%dmmowtb, dmmowtb)
            maxdaymow = cfg%maxdaymow
         end if
      else if (cfg%swharv == 2) then
         ! Populate dateharvest from mowing_dates (DOY floats → t1900 timestamps)
         ! This mirrors readgrass lines 3941-3943: rdatim converts calendar dates;
         ! here cfg%mowing_dates are already DOY, so convert per simulation year.
         ! The exact conversion mirrors readgrass's dateharvest fill:
         ! dateharvest(i) = tend of simulation start year + mowing_dates(i) - 1
         ! (Implementation must inspect grass(task=1) code at lines 2073+ for
         ! how dateharvest is consumed; likely it is t1900-day-of-year float.)
         call populate_dateharvest(cfg)
         dateharvest(cfg%nmow + 1) = tend + 1.0d0
      end if

      if (allocated(cfg%dmmowdelay)) then
         call copy_table(cfg%dmmowdelay, DelayRegrowthTab)
      end if

      ! --- Runtime init math ---
      ! Build cumdens from rdctb (verbatim port of readgrass lines 4091-4121)
      if (swdrought == 1) then
         block
            integer :: i
            real(real64) :: depth, sum_val
            real(real64), dimension(202) :: rootdis
            do i = 0, 100
               depth = 0.01d0 * dble(i)
               rootdis(i*2+1) = depth
               rootdis(i*2+2) = afgen(rdctb, 22, depth)
            end do
            do i = 1, 202, 2
               cumdens(i) = rootdis(i)
            end do
            sum_val = 0.0d0
            cumdens(2) = 0.0d0
            do i = 4, 202, 2
               sum_val = sum_val + (rootdis(i-2) + rootdis(i)) * 0.5d0 &
                                   * (cumdens(i-1) - cumdens(i-3))
               cumdens(i) = sum_val
            end do
            do i = 2, 202, 2
               cumdens(i) = cumdens(i) / sum_val
            end do
         end block
      end if

   end subroutine cropgrass_init_from_config

   ! Private helper: copy a flat allocatable array to a fixed-size module global.
   subroutine copy_table(src, dst)
      real(real64), allocatable, intent(in)  :: src(:)
      real(real64),              intent(out) :: dst(:)
      integer :: n
      if (.not. allocated(src)) return
      n = min(size(src), size(dst))
      dst(1:n) = src(1:n)
   end subroutine

   ! Private helper: populate dateharvest from cfg%mowing_dates (DOY floats).
   subroutine populate_dateharvest(cfg)
      type(cropgrass_config_t), intent(in) :: cfg
      ! Implementation: the dateharvest module global stores t1900-relative
      ! timestamps. readgrass uses rdatim which reads YYYY-MM-DD strings.
      ! cropgrass_init stores DOY floats and must produce the same t1900 values.
      ! (Full implementation requires inspecting variables module for
      ! dateharvest type and the rdatim → t1900 conversion formula.)
      ! Placeholder: copy DOY values directly if dateharvest is DOY-based.
      integer :: i
      do i = 1, cfg%nmow
         dateharvest(i) = cfg%mowing_dates(i)
      end do
   end subroutine

end module cropgrass_init_mod
```

**Important implementation note on `dateharvest`:** inspect the exact type and
interpretation of `dateharvest` in `src/variables.f90` before implementing
`populate_dateharvest`. The legacy `readgrass` uses `rdatim` to read calendar
date strings (e.g. `1980-05-06`) and converts them to `t1900`-based real(8)
values. The TOML stub stores them as DOY floats. The correct conversion may
require adding the year offset for each simulation year, or storing them
differently. This is the highest-risk point of the implementation.

- [ ] **Step 6: Register in meson.build**

In `meson.build`, add `src/crop/cropgrass_init.f90` to the crop sources list.
In `tests/unit/meson.build`, add the new `tests/unit/crop/` subdirectory.
In `tests/unit/testSuites.inc`, register `TestCropgrass_init`.

- [ ] **Step 7: Run tests, verify they pass**

```bash
pixi run test-pfunit 2>&1 | grep -E "cropgrass_init|PASS|FAIL" | tail -20
```

Expected: all 5 new tests pass.

- [ ] **Step 8: Commit**

```bash
git add src/crop/cropgrass_init.f90 meson.build \
        tests/unit/meson.build tests/unit/testSuites.inc \
        tests/unit/crop/test_cropgrass_init.pf
git commit -m "feat(cropgrass-init): runtime init from typed config, mirrors readgrass tail (Phase 3)"
```

---

## Task 8: Wire `cropgrowth.f90:2068` dispatch

Add the `if/else` guard around the `readgrass` call at line 2068 of
`cropgrowth.f90`.

**Files:**
- Modify: `src/crop/cropgrowth.f90`

- [ ] **Step 1: Read the wiring site**

```bash
sed -n '2060,2080p' src/crop/cropgrowth.f90 | cat
```

Confirm the exact lines of the `readgrass` call and surrounding context.

- [ ] **Step 2: Confirm actual array and sentinel names**

```bash
grep -n "rotation_grass\|rotation_cropgrass\|rotation_loaded\|populated" \
    src/config/crop_config.f90 src/crop/crop_config_global.f90 | head -20
```

Record the actual name used. Let `ARRAY_NAME` denote the confirmed name (e.g.
`rotation_grass`) and `SENTINEL` denote how it is checked (`rotation_loaded(icrop)`
or `rotation_grass(icrop)%populated`).

- [ ] **Step 3: Add `use cropgrass_init_mod` to `grass(task)` subroutine**

At the top of the `grass(task)` subroutine's `use` list (around line 2004),
add:

```fortran
use cropgrass_init_mod, only: cropgrass_init_from_config
use crop_config_global_mod, only: crop_config_global
```

(If `crop_config_global_mod` is already used elsewhere in `cropgrowth.f90`
because Phase 1 added it, ensure no duplicate use statement.)

- [ ] **Step 4: Replace the bare `readgrass` call with the dispatch**

Replace (at line ~2068):
```fortran
call readgrass (icrop,cropfil(icrop),swharvest,dmharvest,daylastharvest,dmlastharvest,swdmmow, &
                maxdaymow,swlossmow,swlossgrz,swdmgrz, &
                maxdaygrz,dmgrazing,LSDb,tagprest,swhydrlift)
```

With:
```fortran
if (associated(crop_config_global) .and. &
    allocated(crop_config_global%ARRAY_NAME) .and. &
    SENTINEL) then
   call cropgrass_init_from_config( &
           crop_config_global%ARRAY_NAME(icrop), icrop)
else
   call readgrass (icrop, cropfil(icrop), swharvest, dmharvest,          &
                   daylastharvest, dmlastharvest, swdmmow,                 &
                   maxdaymow, swlossmow, swlossgrz, swdmgrz,               &
                   maxdaygrz, dmgrazing, LSDb, tagprest, swhydrlift)  ! transitional
end if
```

Substitute the actual array name and sentinel form. The `else` branch is
**transitional** — Phase 4 removes it.

- [ ] **Step 5: Build and run unit tests**

```bash
pixi run test-pfunit 2>&1 | tail -10
```

Expected: all existing and new tests still pass. The wiring change has no
effect on the unit test suite (unit tests don't exercise the full runtime
dispatch).

- [ ] **Step 6: Run regression**

```bash
pixi run regression 2>&1 | tail -10
```

Expected: 5/5 green. Because `grassd.crp.toml` files have been authored
(Tasks 5 and 6) and the loader now populates the cache, the `if` branch
should now be taken for cases 4 and 2. Verify by temporarily adding a
print statement or by checking that removing `grassd.crp` from a TOML dir
still passes (the smoke test in Task 11 confirms this more formally).

- [ ] **Step 7: Commit**

```bash
git add src/crop/cropgrowth.f90
git commit -m "feat(cropgrowth): dispatch grass task=1 on rotation_grass populated sentinel (Phase 3)"
```

---

## Task 9: Parity test

Write or extend `tests/unit/io/toml/test_grasscrop_parity.pf` to verify that
`cropgrass_init_from_config` produces the same module-global state as the
legacy `readgrass` on a shared fixture.

**Files:**
- Create/extend: `tests/unit/io/toml/test_grasscrop_parity.pf`

- [ ] **Step 1: Identify an appropriate fixture**

The parity test needs a minimal `grassd.crp`-like file that the legacy
`readgrass` can read. Use `tests/swap-cases/4.oxygenstress/grassd.crp` as the
fixture (it lives in the legacy case dir, not the TOML dir, so it will not be
deleted in Task 12).

- [ ] **Step 2: Write the parity test**

```fortran
@test
subroutine test_grasscrop_parity_case4_kdif()
   ! After calling readgrass on case-4's grassd.crp and calling
   ! cropgrass_init_from_config on the parsed TOML, both should write
   ! the same value to the kdif module global.
   use funit
   use variables, only: kdif
   use cropgrass_config_mod, only: cropgrass_config_t
   use cropgrass_init_mod, only: cropgrass_init_from_config

   type(cropgrass_config_t) :: cfg
   real(real64) :: kdif_from_toml, kdif_from_legacy

   ! ... (set up cfg from case-4 values, call init, capture kdif)
   cfg%kdif = 0.60_real64
   ! (minimal setup)
   call cropgrass_init_from_config(cfg, icrop=1)
   kdif_from_toml = kdif

   ! (call readgrass on the legacy fixture, capture kdif)
   ! ... (use legacy_grass_helper or direct readgrass call)
   kdif_from_legacy = kdif   ! after readgrass

   @assertRelativelyEqual(kdif_from_legacy, kdif_from_toml, tolerance=1.0e-10_real64)
end subroutine
```

The exact structure depends on how Phase 1 wired parity tests
(check `tests/unit/io/toml/test_surfacewater_crop_parity.pf` for the pattern).

- [ ] **Step 3: Run and verify**

```bash
pixi run test-pfunit 2>&1 | grep -E "grasscrop_parity|FAIL" | head -10
```

- [ ] **Step 4: Commit**

```bash
git add tests/unit/io/toml/test_grasscrop_parity.pf
git commit -m "test(cropgrass-parity): parity check cropgrass_init vs readgrass for cases 4+2 (Phase 3)"
```

---

## Task 10: Full test pass

Run all tests; confirm baseline green before the smoke test.

- [ ] **Step 1: Unit tests**

```bash
pixi run test-pfunit 2>&1 | tail -15
```

Expected: all green. Fix any regressions before proceeding.

- [ ] **Step 2: Regression tests**

```bash
pixi run regression 2>&1 | tail -15
```

Expected: 5/5 green.

- [ ] **Step 3: Spot-check case 4 output column**

The case 4 regression reference includes `treddry` and `tredwet` columns
(oxygen stress reduction factors written when `swoxygen=2`). Confirm the
column values match the reference after the wiring change:

```bash
# Inspect the case 4 expected output file for tredwet/treddry presence
grep -l "tredwet\|treddry" tests/swap-cases/regression/4.oxygenstress/ 2>/dev/null | head
```

If the reference has these columns, they must match after Phase 3.

---

## Task 11: Smoke test (both cases)

Verify the new executable path is fully decoupled from the `.crp` files in the
TOML case directories.

- [ ] **Step 1: Disable `grassd.crp` in both TOML case directories**

```bash
mv tests/swap-cases/toml/4.oxygenstress/grassd.crp \
   tests/swap-cases/toml/4.oxygenstress/grassd.crp.disabled
mv tests/swap-cases/toml/2.grassgrowth/grassd.crp \
   tests/swap-cases/toml/2.grassgrowth/grassd.crp.disabled
```

- [ ] **Step 2: Run regression**

```bash
pixi run regression 2>&1 | tail -15
```

Expected: 5/5 green. If either case 4 or case 2 fails, the `grassd.crp` file
is still being opened on the runtime path. Debug before proceeding.

- [ ] **Step 3: Restore**

```bash
mv tests/swap-cases/toml/4.oxygenstress/grassd.crp.disabled \
   tests/swap-cases/toml/4.oxygenstress/grassd.crp
mv tests/swap-cases/toml/2.grassgrowth/grassd.crp.disabled \
   tests/swap-cases/toml/2.grassgrowth/grassd.crp
```

Do not commit after this step — restoration only.

---

## Task 12a: Delete `4.oxygenstress/grassd.crp` (submodule pair commit)

**Files:**
- Submodule delete: `tests/swap-cases/toml/4.oxygenstress/grassd.crp`

- [ ] **Step 1: Delete the file**

```bash
rm tests/swap-cases/toml/4.oxygenstress/grassd.crp
```

- [ ] **Step 2: Verify regression still green**

```bash
pixi run regression 2>&1 | tail -10
```

Expected: 5/5 green. If not, the smoke test (Task 11) would have caught this.

- [ ] **Step 3: Submodule inner commit + outer bump**

```bash
cd tests/swap-cases
git add -u toml/4.oxygenstress/grassd.crp
git commit toml/4.oxygenstress/grassd.crp \
    -m "feat(case4): delete legacy grassd.crp from TOML dir (Phase 3 port complete)"
cd /home/zawadzkim/Code/swap
git add tests/swap-cases
git commit -m "chore(submodule): bump swap-cases for case 4 grassd.crp deletion (Phase 3)"
```

---

## Task 12b: Delete `2.grassgrowth/grassd.crp` (submodule pair commit)

**Files:**
- Submodule delete: `tests/swap-cases/toml/2.grassgrowth/grassd.crp`

- [ ] **Step 1: Delete the file**

```bash
rm tests/swap-cases/toml/2.grassgrowth/grassd.crp
```

- [ ] **Step 2: Verify regression still green**

```bash
pixi run regression 2>&1 | tail -10
```

Expected: 5/5 green.

- [ ] **Step 3: Submodule inner commit + outer bump**

```bash
cd tests/swap-cases
git add -u toml/2.grassgrowth/grassd.crp
git commit toml/2.grassgrowth/grassd.crp \
    -m "feat(case2): delete legacy grassd.crp from TOML dir (Phase 3 port complete)"
cd /home/zawadzkim/Code/swap
git add tests/swap-cases
git commit -m "chore(submodule): bump swap-cases for case 2 grassd.crp deletion (Phase 3)"
```

---

## Task 13: Documentation + final verification

**Files:**
- Modify: `docs/csv-companion-files.md` (if not already updated by Phase 1)

- [ ] **Step 1: Update `docs/csv-companion-files.md`**

Add or verify a note in the "Path resolution and staging" section:

> `.crp.toml` files are resolved relative to the case working directory
> (the directory containing `swap.toml`), matching the convention for CSV
> companion files. Each `[[crop.rotation]]` entry with `type = 3` loads its
> crop data from the file named in the `file = "..."` field.

If Phase 1 already added this note, verify it covers type=3 (grassgrass) as
well as type=1.

- [ ] **Step 2: Final acceptance check**

```bash
cd /home/zawadzkim/Code/swap
pixi run test-pfunit 2>&1 | tail -5
pixi run regression 2>&1 | tail -10
git grep "call readgrass" src/ | grep -v "! transitional"
ls tests/swap-cases/toml/4.oxygenstress/grassd.crp 2>/dev/null && echo "EXISTS (ERROR)" || echo "Deleted (OK)"
ls tests/swap-cases/toml/2.grassgrowth/grassd.crp  2>/dev/null && echo "EXISTS (ERROR)" || echo "Deleted (OK)"
```

Expected:
- Unit tests: all green.
- Regression: 5/5 green.
- `git grep` finds only the `else` transitional fallback in `cropgrowth.f90`
  (not a bare call).
- Both `grassd.crp` files deleted from TOML case dirs.

- [ ] **Step 3: Commit docs**

```bash
git add docs/csv-companion-files.md
git commit -m "docs(crp-port-phase3): note .crp.toml path resolution for type=3 rotations"
```

---

## Acceptance gate summary

All of the following must be true before Phase 3 is considered complete:

1. `pixi run test-pfunit` → all green (including new stub-error, round-trip,
   init-globals, cumdens, parity tests).
2. `pixi run regression` → 5/5 green (including cases 4 and 2).
3. `tests/swap-cases/toml/4.oxygenstress/grassd.crp` does not exist.
4. `tests/swap-cases/toml/2.grassgrowth/grassd.crp` does not exist.
5. `git grep "call readgrass" src/` → returns exactly the transitional `else`
   branch in `cropgrowth.f90:2068` and no other sites.
6. `git grep "call readgrass" tests/` → returns only parity-test references.
7. `crop_config_global` exists and `rotation_grass(:)` is populated for both
   case 4 (10 slots) and case 2 (all slots) with `populated=.true.`.
8. Both case 4 and case 2 produce regression output identical to their
   pre-Phase-3 reference files (no numeric change).
