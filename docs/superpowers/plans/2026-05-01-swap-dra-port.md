# swap.dra → TOML port — implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Eliminate the legacy `swap.dra` ASCII file from the TOML pipeline. After this work the new executable reads only `swap.toml`, `swap.dra.toml`, and CSV companions. Case 6 (surfacewater) runs end-to-end without `tests/swap-cases/toml/6.surfacewater/swap.dra` on disk.

**Architecture:** Replace the runtime call to `rddre` (642-line legacy reader) with a new `surfacewater_init` module that does only the runtime numerical initialization (sttab, wls1, swstini, wlsbak). Push validation into config validators (cross-section in `swap_config_validate`, within-section stubs in `surface_water_config_validate`) and coordinate normalization into finalizers (`hbweir -= altcu` in `swap_config_finalize`, `wldip = abs(wldip)` in `surface_water_config_finalize`). Author the missing `[[drainage.levels]]` table in case 6's `swap.dra.toml`.

**Tech Stack:** Fortran 2008 (gfortran 13), pFUnit unit tests, Meson + Pixi, tomlf TOML parser.

**Spec:** `docs/superpowers/specs/2026-05-01-swap-dra-port-design.md`

---

## File map

**Source (modify):**
- `src/io/toml/read_drainage_toml.f90` — extend `[[drainage.levels]]` parser to read 7 missing per-level fields
- `src/config/surface_water_config.f90` — add stub-error checks for unimplemented branches; finalize wldip
- `src/config/swap_config.f90` — add cross-section validators; add cross-section finalize for `hbweir -= altcu`
- `src/drainage/surfacewater.f90` — change line 56 from `rddre` to `surfacewater_init`
- `meson.build` — register new `surfacewater_init.f90`
- `tests/unit/meson.build` — register new test files
- `tests/unit/testSuites.inc` — register new test suites

**Source (create):**
- `src/drainage/surfacewater_init.f90` — new module with runtime initialization (`sttab`, `swstini`, `wlsbak`)

**Tests (create / modify):**
- `tests/unit/io/toml/test_read_drainage_toml.pf` — extend with per-level field round-trip
- `tests/unit/config/test_surface_water_config.pf` — extend with stub-error tests + wldip finalize
- `tests/unit/config/test_swap_config.pf` — extend with cross-section validators / finalize
- `tests/unit/drainage/test_surfacewater_init.pf` — new file (creates `tests/unit/drainage/` dir)

**Test fixture (modify):**
- `tests/swap-cases/toml/6.surfacewater/swap.dra.toml` (in submodule) — add `cofani` + `[[drainage.levels]]` table; bump `nrlevs` 0 → 2; add `[drainage.surface_runoff]` if needed
- `tests/swap-cases/toml/6.surfacewater/swap.dra` (in submodule) — DELETE in final task

**Docs (modify):**
- `docs/csv-companion-files.md` — remove the `swap.dra` (`swdra=2`) parenthetical from "Path resolution and staging"
- `docs/configuration-schema.md` — extend `[[drainage.levels]]` row to cover the 7 new fields if not already covered

---

## Submodule discipline

`tests/swap-cases/` is a git submodule. Any change inside it is committed in the inner repo first, then the outer repo bumps the SHA. The inner-commit + outer-bump pair is non-negotiable; never commit one without the other. Fixture changes live in tasks 3 and 11.

---

## Task 1: rddre line-by-line audit

Read `readswap.f90:4273-4912` end-to-end and produce a classification table that informs Tasks 4–8. No code changes in this task. The audit prevents missed invariants when extracting math.

**Files:**
- Create: `docs/phase-4f-rddre-audit.md`

- [ ] **Step 1: Read the entire `rddre` body**

Use the Read tool to load `src/io/readswap.f90` lines 4273-4912 in one or two passes.

- [ ] **Step 2: Classify each non-blank, non-comment line into a bucket**

Buckets:
1. **READ** — `rd*` calls (file parsing). Goes away with port. Discard.
2. **VALIDATE** — `if (...) call fatalerr(...)`. Push to validators.
3. **NORMALIZE** — coordinate / unit conversion (e.g. `hbweir -= altcu`). Push to finalizers.
4. **RUNTIME** — sttab build, swstini, wlsbak, wls1/wlp1 interpolation. Goes into `surfacewater_init`.
5. **GUARDED** — branches our scope doesn't support (swsrf=3, swsec=1, swqhr=2, nrman2>0). Stub at the validator level; don't port.

- [ ] **Step 3: Write the audit table**

Create `docs/phase-4f-rddre-audit.md` with a Markdown table of the form:

```markdown
| Lines       | Bucket    | Notes                                      |
| ----------- | --------- | ------------------------------------------ |
| 4329-4331   | READ      | Open swap.dra, init parser                 |
| 4334-4339   | READ      | swdivd                                     |
| 4334-4339   | VALIDATE  | swdivd=0 warning                           |
| 4348-4352   | READ      | cofani                                     |
| 4596-4599   | NORMALIZE | hbweir -= altcu, alphaw normalize          |
| 4600-4603   | VALIDATE  | hbweir < zbotdr fatal                      |
| 4862-4897   | RUNTIME   | sttab build (open-channel storage table)   |
| 4899-4906   | RUNTIME   | swstini, wlsbak                            |
| 4909        | READ      | CLOSE(DRA)                                 |
| 4477-4500   | GUARDED   | swsrf=3 primary system (swdra=2 in case 6) |
| ...         | ...       | ...                                        |
```

Cover all 642 lines. Group consecutive same-bucket lines.

- [ ] **Step 4: Commit**

```bash
git add docs/phase-4f-rddre-audit.md
git commit -m "docs(audit): rddre line-by-line classification for swap.dra port"
```

---

## Task 2: Extend `[[drainage.levels]]` parser

Read the 7 surface-water-extended per-level fields that the schema declares but the parser doesn't yet read. Validates by round-tripping a fixture through the reader.

**Files:**
- Modify: `src/io/toml/read_drainage_toml.f90` (around lines 195-206 — the levels loop)
- Test: `tests/unit/io/toml/test_read_drainage_toml.pf`
- Fixture: `tests/unit/io/toml/fixtures/drainage_levels_extended.toml` (create)

- [ ] **Step 1: Write the failing test**

Add to `tests/unit/io/toml/test_read_drainage_toml.pf`:

```fortran
@test
subroutine test_drainage_levels_parses_extended_fields()
   use funit
   use error_mod
   use drainage_config_mod, only: drainage_config_t
   use read_drainage_toml_mod, only: read_drainage_toml
   use tomlf, only: toml_table, toml_parse
   type(drainage_config_t)        :: c
   type(error_collection_t)       :: errors
   type(toml_table), allocatable  :: doc
   integer :: io
   open(newunit=io, &
        file='tests/unit/io/toml/fixtures/drainage_levels_extended.toml', &
        status='old', action='read')
   call toml_parse(doc, io)
   close(io)
   call read_drainage_toml(doc, c, errors)
   @assertFalse(errors%has_errors())
   @assertEqual(2, size(c%gwlinf))
   @assertEqual(-1000.0d0, c%gwlinf(1), 1.0d-12)
   @assertEqual(-1000.0d0, c%gwlinf(2), 1.0d-12)
   @assertEqual(461.0d0,   c%rdrain(1), 1.0d-12)
   @assertEqual(296.0d0,   c%rdrain(2), 1.0d-12)
   @assertEqual(1000.0d0,  c%rinfi(1),  1.0d-12)
   @assertEqual(500.0d0,   c%rinfi(2),  1.0d-12)
   @assertEqual(1.0d0,     c%rentry(1), 1.0d-12)
   @assertEqual(5.0d0,     c%rexit(1),  1.0d-12)
   @assertEqual(100.0d0,   c%widthr(1), 1.0d-12)
   @assertEqual(50.0d0,    c%widthr(2), 1.0d-12)
   @assertEqual(0.66d0,    c%taludr(1), 1.0d-12)
   @assertEqual(1.00d0,    c%taludr(2), 1.0d-12)
end subroutine
```

- [ ] **Step 2: Create the fixture**

Create `tests/unit/io/toml/fixtures/drainage_levels_extended.toml`:

```toml
[drainage]
swdra    = 2
dramet   = 0
swdivd   = 1
swdislay = 0
nrlevs   = 2
altcu    = 0.0

[[drainage.levels]]
swdtyp = 0
L      = 390.0
zbotdr = -115.0
gwlinf = -1000.0
rdrain = 461.0
rinfi  = 1000.0
rentry = 1.0
rexit  = 5.0
widthr = 100.0
taludr = 0.66

[[drainage.levels]]
swdtyp = 0
L      = 170.0
zbotdr = -60.0
gwlinf = -1000.0
rdrain = 296.0
rinfi  = 500.0
rentry = 1.0
rexit  = 5.0
widthr = 50.0
taludr = 1.00
```

- [ ] **Step 3: Run the test, verify failure**

```bash
pixi run unit-tests 2>&1 | grep -E "test_drainage_levels_parses_extended_fields|Failures: |Ok: " | head -10
```

Expected: failure (the assertions on `c%gwlinf(1)` etc. will trip with default 0.0).

- [ ] **Step 4: Extend the parser**

In `src/io/toml/read_drainage_toml.f90`, locate the `do i = 1, n` loop in `read_drainage_inner` (around line 195). After the existing `call get_optional_real_with_default(item, 'L', ...)` call, add (preserving the existing fields):

```fortran
               call get_optional_real_with_default(item, 'gwlinf',  config%gwlinf(i),  0.0_real64,  'drainage.levels.gwlinf',  errors)
               call get_optional_real_with_default(item, 'rdrain',  config%rdrain(i),  0.0_real64,  'drainage.levels.rdrain',  errors)
               call get_optional_real_with_default(item, 'rinfi',   config%rinfi(i),   0.0_real64,  'drainage.levels.rinfi',   errors)
               call get_optional_real_with_default(item, 'rentry',  config%rentry(i),  0.0_real64,  'drainage.levels.rentry',  errors)
               call get_optional_real_with_default(item, 'rexit',   config%rexit(i),   0.0_real64,  'drainage.levels.rexit',   errors)
               call get_optional_real_with_default(item, 'widthr',  config%widthr(i),  0.0_real64,  'drainage.levels.widthr',  errors)
               call get_optional_real_with_default(item, 'taludr',  config%taludr(i),  0.0_real64,  'drainage.levels.taludr',  errors)
```

- [ ] **Step 5: Run test, verify pass**

```bash
pixi run unit-tests 2>&1 | tail -5
```

Expected: `Ok: 1` (all green).

- [ ] **Step 6: Commit**

```bash
git add src/io/toml/read_drainage_toml.f90 tests/unit/io/toml/test_read_drainage_toml.pf tests/unit/io/toml/fixtures/drainage_levels_extended.toml
git commit -m "feat(io/toml): read surface-water-extended fields in [[drainage.levels]]"
```

---

## Task 3: Author `[[drainage.levels]]` + `cofani` in case 6

Translate the legacy `swap.dra` Part 0 (`COFANI = 10.0 10.0 10.0`) and Part 1 (NRSRF=2 table) into TOML. Verifies via the existing case-6 load test.

**Files:**
- Modify (in submodule): `tests/swap-cases/toml/6.surfacewater/swap.dra.toml`
- Test: `tests/unit/io/toml/test_load_swap_config.pf` (add a case-6 specific assertion if not present)

- [ ] **Step 1: Inspect current state of case 6's swap.dra.toml**

```bash
cat tests/swap-cases/toml/6.surfacewater/swap.dra.toml
```

- [ ] **Step 2: Write the failing assertion**

Add to `tests/unit/io/toml/test_load_swap_config.pf`:

```fortran
@test
subroutine test_load_case6_dra_levels_populated()
   use funit
   use error_mod, only: error_collection_t
   use swap_config_mod, only: swap_config_t
   use load_swap_config_mod, only: load_swap_config
   type(swap_config_t)      :: c
   type(error_collection_t) :: errors
   call load_swap_config('tests/swap-cases/toml/6.surfacewater/swap.toml', c, errors)
   @assertFalse(errors%has_errors())
   @assertEqual(2, c%drain%nrlevs)
   @assertEqual(2, size(c%drain%rdrain))
   @assertEqual(461.0d0, c%drain%rdrain(1), 1.0d-12)
   @assertEqual(296.0d0, c%drain%rdrain(2), 1.0d-12)
   @assertEqual(3, size(c%drain%cofani))
   @assertEqual(10.0d0, c%drain%cofani(1), 1.0d-12)
end subroutine
```

- [ ] **Step 3: Run, verify failure**

```bash
pixi run unit-tests 2>&1 | grep -E "test_load_case6_dra_levels|Failures: |Ok: " | tail -5
```

Expected: failure (nrlevs is 0 in current swap.dra.toml, cofani is missing).

- [ ] **Step 4: Edit case 6's swap.dra.toml**

Replace the file `tests/swap-cases/toml/6.surfacewater/swap.dra.toml` with:

```toml
# surfacewater drainage configuration
# Values from tests/swap-cases/6.surfacewater/swap.dra (legacy ASCII)
# SWDRA=2 (extended drainage with surface water management).
# Surface water management lives in swap.toml's [surface_water] block;
# this file holds the drainage-side fields (basic + per-level table).

[drainage]
swdra    = 2
dramet   = 0
swdivd   = 1
swdislay = 0
nrlevs   = 2
altcu    = 0.0
cofani   = [10.0, 10.0, 10.0]

# 2-row Part 1 table from legacy swap.dra. Columns:
# SWDTYP, L (m), ZBOTDRE (cm), GWLINF (cm), RDRAIN, RINFI, RENTRY, REXIT,
# WIDTHR (cm), TALUDR.
[[drainage.levels]]
swdtyp = 0
L      = 390.0
zbotdr = -115.0
gwlinf = -1000.0
rdrain = 461.0
rinfi  = 1000.0
rentry = 1.0
rexit  = 5.0
widthr = 100.0
taludr = 0.66

[[drainage.levels]]
swdtyp = 0
L      = 170.0
zbotdr = -60.0
gwlinf = -1000.0
rdrain = 296.0
rinfi  = 500.0
rentry = 1.0
rexit  = 5.0
widthr = 50.0
taludr = 1.00
```

- [ ] **Step 5: Run test, verify pass**

```bash
pixi run unit-tests 2>&1 | tail -5
```

Expected: `Ok: N+1` (all green).

- [ ] **Step 6: Submodule inner-commit**

```bash
cd tests/swap-cases
git add toml/6.surfacewater/swap.dra.toml
git commit -m "feat(toml/6.surfacewater): author [[drainage.levels]] table + cofani

Translates Part 0 (COFANI=10.0 10.0 10.0) and Part 1 (NRSRF=2 with
SWDTYP/L/ZBOTDRE/GWLINF/RDRAIN/RINFI/RENTRY/REXIT/WIDTHR/TALUDR) from
legacy swap.dra into typed TOML form. swap.dra remains on disk for
the parity test until the runtime path is decoupled."
cd ../..
```

- [ ] **Step 7: Outer-repo bump-and-test commit**

```bash
git add tests/swap-cases tests/unit/io/toml/test_load_swap_config.pf
git commit -m "test(toml): assert case-6 drainage.levels populated; bump submodule"
```

---

## Task 4: Stub-error validators for unimplemented branches

Reject `swsrf=3`, `swsec=1`, `swqhr=2`, and `swman[i]=2` at config validation time with a clear "not yet supported in TOML pipeline" message. Defense in depth: future-task `surfacewater_init` will guard the same conditions.

**Files:**
- Modify: `src/config/surface_water_config.f90`
- Test: `tests/unit/config/test_surface_water_config.pf`

- [ ] **Step 1: Write the failing tests**

Add to `tests/unit/config/test_surface_water_config.pf`:

```fortran
@test
subroutine test_surface_water_swsrf3_rejected()
   use funit
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   use surface_water_config_mod, only: surface_water_config_t
   type(surface_water_config_t) :: s
   type(error_collection_t)     :: errors
   integer :: i
   logical :: found
   s%swsrf = 3
   call s%validate(errors)
   found = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD .and. &
          index(errors%items(i)%message, 'swsrf=3') > 0) found = .true.
   end do
   @assertTrue(found)
end subroutine

@test
subroutine test_surface_water_swsec1_rejected()
   use funit
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   use surface_water_config_mod, only: surface_water_config_t
   type(surface_water_config_t) :: s
   type(error_collection_t)     :: errors
   integer :: i
   logical :: found
   s%swsrf = 2
   s%swsec = 1
   call s%validate(errors)
   found = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD .and. &
          index(errors%items(i)%message, 'swsec=1') > 0) found = .true.
   end do
   @assertTrue(found)
end subroutine

@test
subroutine test_surface_water_swqhr2_rejected()
   use funit
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   use surface_water_config_mod, only: surface_water_config_t
   type(surface_water_config_t) :: s
   type(error_collection_t)     :: errors
   integer :: i
   logical :: found
   s%swsrf = 2
   s%swsec = 2
   s%swqhr = 2
   call s%validate(errors)
   found = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD .and. &
          index(errors%items(i)%message, 'swqhr=2') > 0) found = .true.
   end do
   @assertTrue(found)
end subroutine

@test
subroutine test_surface_water_swman2_rejected()
   use funit
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   use surface_water_config_mod, only: surface_water_config_t
   type(surface_water_config_t) :: s
   type(error_collection_t)     :: errors
   integer :: i
   logical :: found
   ! Minimum-viable swsrf=2/swsec=2 fixture with a single swman=2 period.
   s%swsrf = 2
   s%swsec = 2
   s%swqhr = 1
   s%nmper = 1
   s%wlact  = -50.0d0
   s%osswlm =   1.0d0
   s%sofcu  = 100.0d0
   allocate(s%impend(1)); s%impend(1) = 36526.0d0
   allocate(s%swman(1));  s%swman(1)  = 2
   allocate(s%wscap(1));  s%wscap(1)  = 0.0d0
   allocate(s%wldip(1));  s%wldip(1)  = 0.0d0
   allocate(s%intwl(1));  s%intwl(1)  = 1
   allocate(s%hbweir(1)); s%hbweir(1) = -50.0d0
   allocate(s%alphaw(1)); s%alphaw(1) =   1.7d0
   allocate(s%betaw(1));  s%betaw(1)  =   1.5d0
   call s%validate(errors)
   found = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD .and. &
          index(errors%items(i)%message, 'swman=2') > 0) found = .true.
   end do
   @assertTrue(found)
end subroutine
```

- [ ] **Step 2: Run, verify failure**

```bash
pixi run unit-tests 2>&1 | grep -E "test_surface_water_sw(srf3|sec1|qhr2|man2)|Failures: |Ok: " | tail -10
```

Expected: 4 failures.

- [ ] **Step 3: Add the stub-error guards**

In `src/config/surface_water_config.f90`, edit `surface_water_config_validate`. Insert the stub-error block at the very top of the subroutine body, BEFORE the existing `call check_int_enum(self%swsrf, [1, 2], ...)` line. Position is critical: existing code has `if (self%swsrf /= 2) return` which would skip the guards if they came after.

Add `ERR_VALIDATION_CROSS_FIELD` to the existing `use error_mod, only: ...` line at the top of the module if not already imported (it's already imported in the existing module — verify with `grep ERR_VALIDATION_CROSS_FIELD src/config/surface_water_config.f90`).

Code to insert at the top of `surface_water_config_validate` (immediately after the declaration block):

```fortran
      ! Stub-errors for branches not yet supported in the TOML pipeline.
      ! These must run BEFORE the swsrf=2 early-return below; otherwise
      ! swsrf=3 would short-circuit out without raising the error.
      if (self%swsrf == 3) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'surface_water.swsrf=3 (primary system) not yet supported in ' // &
            'the TOML pipeline; use the legacy executable.', 'surface_water')
      end if
      if (self%swsec == 1) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'surface_water.swsec=1 (input water level) not yet supported ' // &
            'in the TOML pipeline; use the legacy executable.', 'surface_water')
      end if
      if (self%swqhr == 2) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'surface_water.swqhr=2 (q-h table discharge) not yet supported ' // &
            'in the TOML pipeline; use the legacy executable.', 'surface_water')
      end if
      if (allocated(self%swman)) then
         if (any(self%swman == 2)) then
            call errors%append(ERR_VALIDATION_CROSS_FIELD, &
               'surface_water.swman=2 (automatic weir) not yet supported ' // &
               'in the TOML pipeline; use the legacy executable.', 'surface_water')
         end if
      end if
```

- [ ] **Step 4: Run, verify pass**

```bash
pixi run unit-tests 2>&1 | tail -5
```

Expected: all green, including the 4 new tests.

- [ ] **Step 5: Commit**

```bash
git add src/config/surface_water_config.f90 tests/unit/config/test_surface_water_config.pf
git commit -m "feat(config): reject swsrf=3/swsec=1/swqhr=2/swman=2 at validate

Branches not yet ported to the TOML pipeline are surfaced with a clear
'use the legacy executable' message rather than running the simulator
to a wrong result. Defense-in-depth gate at validation time; the
surfacewater_init runtime layer guards the same conditions."
```

---

## Task 5: Finalize `wldip = abs(wldip)`

Match the legacy `rddre` behavior of treating `wldip` as a magnitude (legacy: `wldip(imper) = abs(wldip(imper))`). Tiny finalize step within `surface_water_config`.

**Files:**
- Modify: `src/config/surface_water_config.f90` (the `surface_water_config_finalize` sub)
- Test: `tests/unit/config/test_surface_water_config.pf`

- [ ] **Step 1: Write the failing test**

```fortran
@test
subroutine test_surface_water_finalize_wldip_abs()
   use funit
   use error_mod, only: error_collection_t
   use surface_water_config_mod, only: surface_water_config_t
   type(surface_water_config_t) :: s
   type(error_collection_t)     :: errors
   s%swsrf = 2; s%swsec = 2; s%swqhr = 1
   s%nmper = 2; s%sofcu = 100.0d0
   allocate(s%wldip(2));  s%wldip  = [-3.5d0, 2.0d0]
   allocate(s%alphaw(2)); s%alphaw = [1.7d0, 1.7d0]
   allocate(s%betaw(2));  s%betaw  = [1.5d0, 1.5d0]
   call s%finalize(errors)
   @assertEqual(3.5d0, s%wldip(1), 1.0d-12)
   @assertEqual(2.0d0, s%wldip(2), 1.0d-12)
end subroutine
```

- [ ] **Step 2: Run, verify failure**

```bash
pixi run unit-tests 2>&1 | grep -E "wldip_abs|Failures: |Ok: " | tail -5
```

Expected: failure on the `s%wldip(1)` assertion.

- [ ] **Step 3: Extend the finalizer**

In `src/config/surface_water_config.f90`, locate `surface_water_config_finalize`. Before the existing `if (self%swqhr /= 1) return` check, add:

```fortran
      ! Match legacy rddre: wldip is treated as a magnitude (negatives
      ! are silently flipped). Done before swqhr-gated normalization so
      ! it runs regardless of swqhr value.
      if (allocated(self%wldip)) then
         self%wldip = abs(self%wldip)
      end if
```

- [ ] **Step 4: Run, verify pass**

```bash
pixi run unit-tests 2>&1 | tail -5
```

Expected: all green.

- [ ] **Step 5: Commit**

```bash
git add src/config/surface_water_config.f90 tests/unit/config/test_surface_water_config.pf
git commit -m "feat(config): finalize wldip = abs(wldip) (parity with rddre)"
```

---

## Task 6: Cross-section validation in `swap_config`

Two cross-section rules from `rddre`:
- `hbweir[i] (post-altcu-subtraction) > zbotdr(1)` (weir crest above deepest channel bottom of secondary system; `1+nrpri = 1` for swsrf=2).
- `hbweir[i] - wldip[i] > zbotdr(1) + 1e-4` when `swman[i]=1 .and. wscap[i] > 1e-7`.

Both validators run on the post-finalize state (Task 7 lands the altcu subtraction in finalize). The validator therefore expects `hbweir` to already be in altcu-relative coordinates.

**Files:**
- Modify: `src/config/swap_config.f90`
- Test: `tests/unit/config/test_swap_config.pf`

- [ ] **Step 1: Write the failing tests**

```fortran
@test
subroutine test_swap_config_hbweir_below_channel_bottom_fails()
   use funit
   use swap_config_mod, only: swap_config_t
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_SECTION
   type(swap_config_t)      :: c
   type(error_collection_t) :: errors
   integer :: i
   logical :: found
   c%drain%swdra    = 2
   c%drain%nrlevs   = 1
   allocate(c%drain%zbotdr(1)); c%drain%zbotdr(1) = -60.0d0
   c%surface_water%swsrf = 2
   c%surface_water%swsec = 2
   c%surface_water%swqhr = 1
   c%surface_water%nmper = 1
   c%surface_water%sofcu = 100.0d0
   c%surface_water%osswlm = 1.0d0
   c%surface_water%wlact  = -50.0d0
   allocate(c%surface_water%impend(1)); c%surface_water%impend(1) = 36526.0d0
   allocate(c%surface_water%swman(1));  c%surface_water%swman(1)  = 1
   allocate(c%surface_water%wscap(1));  c%surface_water%wscap(1)  = 0.0d0
   allocate(c%surface_water%wldip(1));  c%surface_water%wldip(1)  = 0.0d0
   allocate(c%surface_water%intwl(1));  c%surface_water%intwl(1)  = 1
   allocate(c%surface_water%hbweir(1)); c%surface_water%hbweir(1) = -100.0d0   ! BELOW zbotdr(1)
   allocate(c%surface_water%alphaw(1)); c%surface_water%alphaw(1) =   1.7d0
   allocate(c%surface_water%betaw(1));  c%surface_water%betaw(1)  =   1.5d0
   call c%validate(errors)
   found = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_SECTION .and. &
          index(errors%items(i)%message, 'hbweir') > 0) found = .true.
   end do
   @assertTrue(found)
end subroutine

@test
subroutine test_swap_config_target_below_channel_with_supply_fails()
   use funit
   use swap_config_mod, only: swap_config_t
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_SECTION
   type(swap_config_t)      :: c
   type(error_collection_t) :: errors
   integer :: i
   logical :: found
   c%drain%swdra  = 2
   c%drain%nrlevs = 1
   allocate(c%drain%zbotdr(1)); c%drain%zbotdr(1) = -60.0d0
   c%surface_water%swsrf = 2
   c%surface_water%swsec = 2
   c%surface_water%swqhr = 1
   c%surface_water%nmper = 1
   c%surface_water%sofcu = 100.0d0
   c%surface_water%osswlm = 1.0d0
   c%surface_water%wlact  = -50.0d0
   allocate(c%surface_water%impend(1)); c%surface_water%impend(1) = 36526.0d0
   allocate(c%surface_water%swman(1));  c%surface_water%swman(1)  = 1
   allocate(c%surface_water%wscap(1));  c%surface_water%wscap(1)  = 0.5d0   ! supply attempted
   allocate(c%surface_water%wldip(1));  c%surface_water%wldip(1)  = 5.0d0
   allocate(c%surface_water%intwl(1));  c%surface_water%intwl(1)  = 1
   allocate(c%surface_water%hbweir(1)); c%surface_water%hbweir(1) = -58.0d0  ! crest-wldip = -63, below zbotdr(1)=-60
   allocate(c%surface_water%alphaw(1)); c%surface_water%alphaw(1) = 1.7d0
   allocate(c%surface_water%betaw(1));  c%surface_water%betaw(1)  = 1.5d0
   call c%validate(errors)
   found = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_SECTION .and. &
          index(errors%items(i)%message, 'target level') > 0) found = .true.
   end do
   @assertTrue(found)
end subroutine
```

- [ ] **Step 2: Run, verify failures**

```bash
pixi run unit-tests 2>&1 | grep -E "hbweir_below|target_below|Failures: |Ok: " | tail -5
```

Expected: 2 failures.

- [ ] **Step 3: Add cross-section validators**

In `src/config/swap_config.f90`, extend `swap_config_validate`. After the existing `call self%crop%validate(errors)` and BEFORE the existing soil.swinco=3 block, add:

```fortran
      ! Cross-section validation for surface-water management.
      ! Task 7 enforces drain.altcu = 0 in the TOML pipeline, so all
      ! coordinates here are already in the same (altcu-relative)
      ! frame — no inline altcu subtraction needed. Comparison uses
      ! zbotdr(1+nrpri) where nrpri = 0 for swsrf=2 (the only branch
      ! reaching this code, courtesy of the upstream stub-errors).
      if (self%drain%swdra == 2 .and. self%surface_water%swsrf == 2 .and. &
          self%surface_water%swsec == 2 .and. self%surface_water%swqhr == 1 .and. &
          self%drain%nrlevs >= 1 .and. allocated(self%drain%zbotdr) .and. &
          allocated(self%surface_water%hbweir) .and. &
          allocated(self%surface_water%wldip) .and. &
          allocated(self%surface_water%swman) .and. &
          allocated(self%surface_water%wscap)) then
         block
            integer :: imper
            real(real64) :: zbottom
            zbottom = self%drain%zbotdr(1)
            do imper = 1, self%surface_water%nmper
               if (self%surface_water%hbweir(imper) < zbottom) then
                  call errors%append(ERR_VALIDATION_CROSS_SECTION, &
                     'hbweir below deepest channel bottom of secondary system', &
                     'swap_config')
               end if
               if (self%surface_water%swman(imper) == 1 .and. &
                   self%surface_water%wscap(imper) > 1.0e-7_real64 .and. &
                   (self%surface_water%hbweir(imper) - self%surface_water%wldip(imper)) < &
                   (zbottom + 1.0e-4_real64)) then
                  call errors%append(ERR_VALIDATION_CROSS_SECTION, &
                     'target level (hbweir - wldip) below channel bottom; ' // &
                     'supply not possible', 'swap_config')
               end if
            end do
         end block
      end if
```

Make sure `ERR_VALIDATION_CROSS_SECTION` is on the existing `use error_mod, only: ...` line at the top of the module:

```fortran
   use error_mod, only: error_collection_t, ERR_VALIDATION_REQUIRED, &
                         ERR_VALIDATION_CROSS_SECTION
```

- [ ] **Step 4: Run, verify pass**

```bash
pixi run unit-tests 2>&1 | tail -5
```

Expected: all green.

- [ ] **Step 5: Commit**

```bash
git add src/config/swap_config.f90 tests/unit/config/test_swap_config.pf
git commit -m "feat(config): cross-section validators for hbweir vs zbotdr

Mirrors the legacy rddre fatalerr checks at readswap.f90:4600-4613:
weir crest above deepest channel bottom, and target level (hbweir-wldip)
above bottom when supply is attempted (swman=1, wscap>0)."
```

---

## Task 7: Stub-error for non-zero `altcu`

Legacy `rddre` subtracts `altcu` from `zbotdr`, `hbweir`, and the `wls1` initial level — `altcu` and `wlact` are local-only inside `rddre`, never propagated as module globals. The TOML pipeline doesn't yet have plumbing for these subtractions in the runtime path, and no current TOML case authors `altcu /= 0` (case 6 has `altcu = 0.0`). Rather than build the altcu plumbing speculatively, this port enforces `altcu = 0` at validate time and defers the rest to a future port. This makes Task 6's validator math unambiguous (no inline `altcu` subtraction needed).

**Files:**
- Modify: `src/config/drainage_config.f90` (add the stub-error in `drainage_config_validate`)
- Test: `tests/unit/config/test_drainage_config.pf`

- [ ] **Step 1: Write the failing test**

```fortran
@test
subroutine test_drainage_altcu_nonzero_rejected()
   use funit
   use drainage_config_mod, only: drainage_config_t
   use error_mod, only: error_collection_t, ERR_VALIDATION_CROSS_FIELD
   type(drainage_config_t)  :: d
   type(error_collection_t) :: errors
   integer :: i
   logical :: found
   d%swdra = 2
   d%altcu = 5.0d0   ! non-zero: should be rejected
   call d%validate(errors)
   found = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_CROSS_FIELD .and. &
          index(errors%items(i)%message, 'altcu') > 0) found = .true.
   end do
   @assertTrue(found)
end subroutine
```

- [ ] **Step 2: Run, verify failure**

```bash
pixi run unit-tests 2>&1 | grep -E "altcu_nonzero|Failures: |Ok: " | tail -5
```

Expected: failure (validator currently accepts any altcu).

- [ ] **Step 3: Extend `drainage_config_validate`**

In `src/config/drainage_config.f90`, after the existing `check_int_range(self%nrlevs, ...)` line, add:

```fortran
      ! Stub-error: non-zero altcu requires altcu-subtraction plumbing
      ! in the runtime adapter that the TOML port hasn't built yet.
      ! All current TOML cases author altcu = 0.0; future cases needing
      ! a non-zero altcu must extend the adapter to subtract altcu from
      ! zbotdr / hbweir / wls1 globals before this guard is removed.
      if (abs(self%altcu) > 1.0e-12_real64) then
         call errors%append(ERR_VALIDATION_CROSS_FIELD, &
            'drainage.altcu /= 0 is not yet supported in the TOML pipeline. ' // &
            'Use the legacy executable for cases authoring altcu /= 0.', &
            'drainage')
      end if
```

- [ ] **Step 4: Verify Task 6 tests still pass**

Task 6's test fixtures author `altcu = 0`, so they're unaffected by this guard. Re-run:

```bash
pixi run unit-tests 2>&1 | tail -5
```

Expected: all green.

- [ ] **Step 5: Commit**

```bash
git add src/config/drainage_config.f90 tests/unit/config/test_drainage_config.pf
git commit -m "feat(config): reject drainage.altcu /= 0 in the TOML pipeline

Locks the case-6 subset to altcu = 0 and defers full altcu plumbing
(zbotdr/hbweir/wls1 subtraction) to a future port. Until then, any
case authoring altcu /= 0 must run via the legacy executable."
```

---

## Task 8: New `surfacewater_init` module + `wls1_init` global

Surviving runtime math from `rddre`: `numadj=0`, `wlsbak(1:4)=0`, set `wls1 = wls1_init` (the adapter's pre-computed initial level), `wlstar = wls1`, build `sttab` from `widthr`/`taludr`/`zbotdr`/`nrlevs`, `swstini = swstlev(wls1)`, `swst = swstini`. Operates on already-populated module globals.

Because `wlact` is not a module global today and the adapter currently doesn't propagate it, this task adds a single new global `wls1_init` (a real scalar) to the variables module. The adapter writes `wls1_init = config%surface_water%wlact` (Task 7 enforces `altcu = 0`, so `wls1 = wlact - altcu` simplifies to `wlact`). `surfacewater_init` reads `wls1_init` to seed its OUT parameter `wls1`.

**Files:**
- Modify: `src/core/variables.f90` (add `wls1_init` global)
- Modify: `src/io/toml/config_to_variables.f90` (write `wls1_init`)
- Create: `src/drainage/surfacewater_init.f90`
- Modify: `meson.build`
- Test: `tests/unit/drainage/test_surfacewater_init.pf` (new)
- Modify: `tests/unit/meson.build`
- Modify: `tests/unit/testSuites.inc`

- [ ] **Step 1a: Add `wls1_init` to `src/core/variables.f90`**

Find the surface-water block in `variables.f90` (around line 1259, where `osswlm`, `wlstar`, `wlp` are declared as `real(8)`). Add `wls1_init` to that line:

```fortran
      real(8) osswlm,wlstar,wlp,alphaw(mamp),betaw(mamp)
      real(8) wls1_init    ! TOML pipeline: initial wls1 = wlact - altcu (Task 7 forces altcu=0)
```

- [ ] **Step 1b: Adapter writes `wls1_init`**

In `src/io/toml/config_to_variables.f90`, find the surface-water block (around line 1080-1090, the section that writes `swsrf`, `swsec`, `osswlm`). Replace the comment "wlact + osswlm are loaded by readswap into local scratch; only osswlm is a module global." with active code:

```fortran
      ! TOML pipeline: pre-compute the initial water level wls1.
      ! Legacy rddre computes wls1 = wlact - altcu inside the routine;
      ! we do the same here so surfacewater_init can read it from a
      ! module global. drainage.altcu /= 0 is rejected upstream (Task 7),
      ! so this simplifies to wlact.
      wls1_init = config%surface_water%wlact - config%drain%altcu
      osswlm    = config%surface_water%osswlm
```

(Add `wls1_init` to the `use variables` import block at the top of the adapter if it's not already pulled in via a general-purpose import.)

- [ ] **Step 1c: Create `src/drainage/surfacewater_init.f90`**

```fortran
!> Runtime initialization for the surface-water management system.
!! Replaces the surviving math from the legacy `rddre` reader: build
!! the open-channel storage table sttab, set wls1/wlstar from the
!! adapter-computed wls1_init, and zero the running buffers
!! numadj/wlsbak. All inputs come from module globals already
!! populated by the TOML adapter (config_to_variables).
!!
!! Scope: swsrf=2, swsec=2, swqhr=1, swman=1, drainage.altcu=0 only.
!! Other branches are guarded with fatalerr_collected (defense in
!! depth — the config validator rejects them upstream too).
module surfacewater_init_mod
   use iso_fortran_env, only: real64
   implicit none
   private

   public :: surfacewater_init

contains

   subroutine surfacewater_init(wls1, wlp1)
      use variables, only: nrlevs, swdtyp, zbotdr, widthr, taludr, l, &
                            wls1_init, wlstar, &
                            sttab, swstini, swst, wlsbak, numadj, &
                            swsrf, swsec, swqhr, swman, nmper
      use surfacewater_utils, only: swstlev
      use error_mod, only: fatalerr_collected
      real(real64), intent(out) :: wls1, wlp1
      ! NOTE: wls1/wlp1 here are local OUT parameters; they don't
      ! collide with any module global (verified — variables.f90
      ! has wls1_init, not wls1).

      integer      :: i, ilev
      real(real64) :: wdepth, wvolum, wbreadth
      integer      :: nrpri

      ! Defensive guards mirroring surface_water_config_validate.
      ! swman is a fixed-size array in the variables module
      ! (dimensioned mamp); we slice 1:nmper to compare authored periods.
      if (swsrf == 3 .or. swsec == 1 .or. swqhr == 2) then
         call fatalerr_collected('surfacewater_init', &
            'swsrf=3, swsec=1, or swqhr=2 not supported on the TOML path')
         return
      end if
      if (any(swman(1:nmper) == 2)) then
         call fatalerr_collected('surfacewater_init', &
            'swman=2 (automatic weir) not supported on the TOML path')
         return
      end if

      ! For swsrf=2 (no primary system) nrpri = 0.
      nrpri = 0

      numadj = 0
      do i = 1, 4
         wlsbak(i) = 0.0_real64
      end do

      ! Initial water level pre-computed by the adapter
      ! (= wlact - altcu; altcu enforced = 0 by Task 7).
      wls1   = wls1_init
      wlp1   = 0.0_real64    ! swsrf=2 has no primary system
      wlstar = wls1

      ! Build sttab — storage table indexed 1..22.
      ! Row 1 = +100cm above soil surface; row 2 = 0cm; rows 3..22
      ! divide the column [0, zbotdr(1+nrpri)] into 20 compartments.
      sttab(1, 1) = 100.0_real64
      sttab(2, 1) =   0.0_real64
      do i = 3, 22
         sttab(i, 1) = zbotdr(1 + nrpri) * (i - 2) / 20.0_real64
      end do

      ! sttab(:, 2) — storage in cm as a function of water level, summed
      ! across open-channel levels (swdtyp=0). Closed-channel levels
      ! contribute zero. Verbatim port from readswap.f90:4878-4897.
      do i = 1, 22
         sttab(i, 2) = 0.0_real64
         do ilev = 1 + nrpri, nrlevs
            if (swdtyp(ilev) == 0 .and. sttab(i, 1) > zbotdr(ilev)) then
               if (sttab(i, 1) <= 0.0_real64) then
                  ! Trapezium below soil surface
                  wdepth = sttab(i, 1) - zbotdr(ilev)
                  wvolum = wdepth * (widthr(ilev) + wdepth / taludr(ilev))
               else
                  ! Trapezium up to surface, plus rectangle above
                  wdepth   = -zbotdr(ilev)
                  wvolum   = wdepth * (widthr(ilev) + wdepth / taludr(ilev))
                  wbreadth = widthr(ilev) + 2.0_real64 * wdepth / taludr(ilev)
                  wdepth   = sttab(i, 1)
                  wvolum   = wvolum + wbreadth * wdepth
               end if
               sttab(i, 2) = sttab(i, 2) + wvolum / l(ilev)
            end if
         end do
      end do

      ! Initial storage state.
      swstini = swstlev(wls1)
      swst    = swstini
   end subroutine surfacewater_init

end module surfacewater_init_mod
```

- [ ] **Step 2: Register the new module in `meson.build`**

In the top-level `meson.build`, find the line `'src/drainage/surfacewater.f90',` (around line 131) and add right BEFORE it:

```
    'src/drainage/surfacewater_init.f90',
```

The order matters: `surfacewater.f90` will `use surfacewater_init_mod` so the init module compiles first.

Also confirm `src/core/variables.f90` is already in `sources` (it is — line ~74). The `wls1_init` declaration added in Step 1a needs no separate registration.

- [ ] **Step 3: Build to verify the module compiles**

```bash
pixi run build 2>&1 | tail -20
```

Expected: clean build (no `Module not found` or syntax errors).

- [ ] **Step 4: Write the unit test**

Create `tests/unit/drainage/test_surfacewater_init.pf`:

```fortran
@test
subroutine test_surfacewater_init_zeroes_state()
   use funit
   use surfacewater_init_mod, only: surfacewater_init
   use variables, only: numadj, wlsbak, nrlevs, swdtyp, zbotdr, widthr, &
                         taludr, l, wls1_init, swsrf, swsec, swqhr, &
                         nmper, swman
   real(8) :: wls1, wlp1
   integer :: i

   ! Minimal happy-path setup mimicking case 6's first level.
   ! swman is a fixed-size global (mamp); set element 1 directly.
   swsrf  = 2
   swsec  = 2
   swqhr  = 1
   nrlevs = 1
   nmper  = 1
   wls1_init = -77.0d0    ! adapter would write this from wlact - altcu
   swman(1)  = 1
   swdtyp(1) = 0
   zbotdr(1) = -115.0d0
   widthr(1) = 100.0d0
   taludr(1) = 0.66d0
   l(1)      = 39000.0d0   ! 390m → 39000cm

   numadj = 999
   do i = 1, 4
      wlsbak(i) = 999.0d0
   end do

   call surfacewater_init(wls1, wlp1)

   @assertEqual(0,  numadj)
   @assertEqual(0.0d0, wlsbak(1), 1.0d-12)
   @assertEqual(0.0d0, wlsbak(4), 1.0d-12)
   @assertEqual(-77.0d0, wls1, 1.0d-12)
   @assertEqual(0.0d0,   wlp1, 1.0d-12)
end subroutine

@test
subroutine test_surfacewater_init_sttab_top_rows()
   use funit
   use surfacewater_init_mod, only: surfacewater_init
   use variables, only: sttab, nrlevs, swdtyp, zbotdr, widthr, taludr, l, &
                         wls1_init, swsrf, swsec, swqhr, nmper, swman
   real(8) :: wls1, wlp1

   swsrf  = 2
   swsec  = 2
   swqhr  = 1
   nrlevs = 1
   nmper  = 1
   wls1_init = -77.0d0
   swman(1)  = 1
   swdtyp(1) = 0
   zbotdr(1) = -115.0d0
   widthr(1) = 100.0d0
   taludr(1) = 0.66d0
   l(1)      = 39000.0d0

   call surfacewater_init(wls1, wlp1)

   ! Row 1 = +100cm; row 2 = 0cm.
   @assertEqual(100.0d0, sttab(1, 1), 1.0d-12)
   @assertEqual(  0.0d0, sttab(2, 1), 1.0d-12)
   ! Rows 3..22 divide [0, zbotdr(1)] = [0, -115] into 20 compartments.
   ! Row 3 = -5.75; row 22 = -115.
   @assertEqual(-5.75d0,  sttab( 3, 1), 1.0d-12)
   @assertEqual(-115.0d0, sttab(22, 1), 1.0d-12)
end subroutine
```

- [ ] **Step 5: Wire the new test in `tests/unit/meson.build`**

In `tests/unit/meson.build`, in the `pf_files` list, add:

```
        'drainage/test_surfacewater_init.pf',
```

Insert it after the config tests, before the io/ tests (alphabetical-ish placement is fine).

- [ ] **Step 6: Wire the test suite in `tests/unit/testSuites.inc`**

Add the line (matching the convention of existing entries):

```
ADD_TEST_SUITE(test_surfacewater_init_suite)
```

- [ ] **Step 7: Run, verify pass**

```bash
pixi run unit-tests 2>&1 | tail -10
```

Expected: all green, including the 2 new tests.

- [ ] **Step 8: Commit**

```bash
git add src/drainage/surfacewater_init.f90 meson.build \
        tests/unit/drainage/test_surfacewater_init.pf \
        tests/unit/meson.build tests/unit/testSuites.inc
git commit -m "feat(drainage): surfacewater_init module replaces rddre runtime math

New module does the post-read math from the legacy rddre: zeroes
numadj/wlsbak, computes wls1=wlact-altcu, wlstar=wls1, builds the
open-channel sttab storage table, and primes swstini=swstlev(wls1).
Operates on adapter-populated module globals; no file I/O.

Guards swsrf=3, swsec=1, swqhr=2, and swman=2 with fatalerr_collected
matching the upstream config validator stub-errors."
```

---

## Task 9: Wire surfacewater.f90 to call surfacewater_init

Single-line swap. The callee changes; the caller signature does not.

**Files:**
- Modify: `src/drainage/surfacewater.f90` (lines 50-58, the `case (1)` block)

- [ ] **Step 1: Read the current callsite**

```bash
sed -n '48,60p' src/drainage/surfacewater.f90
```

- [ ] **Step 2: Replace `call rddre` with `call surfacewater_init`**

Edit `src/drainage/surfacewater.f90`:

Before (line 56):
```fortran
      call rddre (wls,wlp)
```

After:
```fortran
      call surfacewater_init (wls, wlp)
```

Add the `use` import. Find the existing `use` block at the top of the routine (likely around line 30). The module's existing imports include `use variables, ...` and `use surfacewater_utils, ...`. Add a new line just below those:

```fortran
      use surfacewater_init_mod, only: surfacewater_init
```

- [ ] **Step 3: Build, verify clean**

```bash
pixi run build 2>&1 | tail -10
```

Expected: clean build.

- [ ] **Step 4: Run unit tests, verify pass**

```bash
pixi run unit-tests 2>&1 | tail -5
```

Expected: all green. (Parity tests still pass because they call `readswap()` directly, which still calls `rddre`.)

- [ ] **Step 5: Run regression case 6**

```bash
pixi run regression-tests 2>&1 | tail -20
```

Expected: 5/5 cases green. The diff should show case 6 still passes (`swap.dra` is still on disk and still being opened by rddre via the parity test, but the runtime path no longer touches it).

- [ ] **Step 6: Commit**

```bash
git add src/drainage/surfacewater.f90
git commit -m "refactor(drainage): route surfacewater task=1 to surfacewater_init

Replaces the runtime call to legacy rddre with the new
surfacewater_init module. rddre stays alive in readswap.f90 as a
parity-test reference; nothing on the runtime path opens swap.dra
anymore."
```

---

## Task 10: Verify case 6 runs without `swap.dra` (smoke test)

Pre-flight before deletion: temporarily rename `swap.dra` and re-run regression to confirm the runtime path is fully decoupled.

**Files:** none modified (transient rename only).

- [ ] **Step 1: Rename swap.dra in submodule**

```bash
cd tests/swap-cases/toml/6.surfacewater && mv swap.dra swap.dra.disabled && cd ../../../..
```

- [ ] **Step 2: Run regression tests**

```bash
pixi run regression-tests 2>&1 | tail -30
```

Expected: 5/5 green. If case 6 fails with "Cannot open swap.dra" or similar, the runtime path is not yet decoupled — investigate before proceeding.

- [ ] **Step 3: Restore the file**

```bash
cd tests/swap-cases/toml/6.surfacewater && mv swap.dra.disabled swap.dra && cd ../../../..
```

- [ ] **Step 4: Verify submodule clean**

```bash
cd tests/swap-cases && git status && cd ..
```

Expected: working tree clean (the rename was reverted).

- [ ] **Step 5: No commit** — this is a verification-only task.

---

## Task 11: Delete `swap.dra` from case 6

Final removal of the legacy ASCII file from the TOML directory. Submodule pair commit.

**Files:**
- Delete (in submodule): `tests/swap-cases/toml/6.surfacewater/swap.dra`

- [ ] **Step 1: Remove the file in the submodule**

```bash
cd tests/swap-cases && git rm toml/6.surfacewater/swap.dra && cd ..
```

- [ ] **Step 2: Confirm parity test still passes**

The parity test chdirs to `tests/swap-cases/6.surfacewater/` (the legacy dir, not the toml dir). That dir's `swap.dra` is unaffected.

```bash
pixi run unit-tests 2>&1 | grep -E "surfacewater_parity|Failures: |Ok: " | tail -5
```

Expected: all green.

- [ ] **Step 3: Run regression**

```bash
pixi run regression-tests 2>&1 | tail -10
```

Expected: 5/5 green.

- [ ] **Step 4: Submodule inner-commit**

```bash
cd tests/swap-cases
git commit -m "chore(toml/6.surfacewater): remove legacy swap.dra

Drainage and surface-water configuration is now fully expressed in
swap.toml + swap.dra.toml. swap.dra remains in tests/swap-cases/6.surfacewater/
for the legacy executable and the surfacewater parity test fixture."
cd ..
```

- [ ] **Step 5: Outer-repo bump**

```bash
git add tests/swap-cases
git commit -m "chore(submodule): bump tests/swap-cases — remove swap.dra from case 6"
```

---

## Task 12: Documentation

Update the two doc files referenced by the spec.

**Files:**
- Modify: `docs/csv-companion-files.md`
- Modify: `docs/configuration-schema.md` (only if `[[drainage.levels]]` row coverage is incomplete)

- [ ] **Step 1: Read the current state**

```bash
grep -n "swap\.dra\|swdra=2" docs/csv-companion-files.md
grep -n "drainage.levels\|gwlinf\|rdrain" docs/configuration-schema.md | head -20
```

- [ ] **Step 2: Edit `docs/csv-companion-files.md`**

Find the "Path resolution and staging" section (around line 92-106) — specifically the bullet that lists what lives in the case working directory. Remove the parenthetical `"and where the scenario requires it the legacy ASCII companion swap.dra (swdra=2, surface-water extended drainage; pending its own port)"`.

After edit:

```markdown
- The case working directory is `tests/swap-cases/toml/<N>.<case>/`. It is
  self-contained — every file SWAP reads at runtime lives there: `swap.toml`,
  `swap.dra.toml`, `*.crp.toml`, all `*.csv` companions, the legacy `*.crp`
  crop files (read by sub-readers in `cropgrowth.f90` until Phase 4f-extend
  ports them), and `swap_linux.swp.template` (staged to `swap.swp` per run).
```

- [ ] **Step 3: Edit `docs/configuration-schema.md` if needed**

If the existing `[[drainage.levels]]` schema table doesn't list `gwlinf`, `rdrain`, `rinfi`, `rentry`, `rexit`, `widthr`, `taludr`, add them with their ranges (taken from legacy `rddre`):

| Field   | Type   | Range          | Notes                            |
| ------- | ------ | -------------- | -------------------------------- |
| gwlinf  | real   | -10000..0      | Groundwater level for max infiltration (cm) |
| rdrain  | real   | 1..1e5         | Drainage resistance (d)          |
| rinfi   | real   | 1..1e5         | Infiltration resistance (d)      |
| rentry  | real   | 0..10          | Entry resistance (d)             |
| rexit   | real   | 0..10          | Exit resistance (d)              |
| widthr  | real   | 0..10000       | Bottom width of channel (cm)     |
| taludr  | real   | 0.01..5        | Side-slope (dh/dw)               |

- [ ] **Step 4: Commit**

```bash
git add docs/csv-companion-files.md docs/configuration-schema.md
git commit -m "docs: remove swap.dra parenthetical; document [[drainage.levels]] surface-water-extended fields"
```

---

## Acceptance gate

After Task 12, verify:

- [ ] `pixi run unit-tests` → all green, including new tests in Tasks 2/4/5/6/7/8.
- [ ] `pixi run regression-tests` → 5/5 cases green.
- [ ] `git grep "call rddre" src/` → returns nothing (only test files reference `rddre`).
- [ ] `git grep "call rddre" tests/` → still returns the parity test fixtures (intentional).
- [ ] `ls tests/swap-cases/toml/6.surfacewater/swap.dra` → file does not exist.
- [ ] `tests/swap-cases/toml/6.surfacewater/swap.dra.toml` has `nrlevs=2` and `[[drainage.levels]]` blocks.

---

## Self-review notes

**Spec coverage:** All 7 spec units covered. Task 4 = stub validators (swsrf=3/swsec=1/swqhr=2/swman=2). Task 5 = wldip finalize. Task 6 = cross-section validators (hbweir vs zbotdr). Task 7 = altcu /= 0 stub-error (defers full altcu plumbing). Task 8 = `surfacewater_init` + `wls1_init` global. Task 9 = wire surfacewater.f90 to call `surfacewater_init`. Tasks 10-11 = remove swap.dra. Task 12 = docs. Plus Task 1 (audit) and Tasks 2-3 (TOML authoring) up front.

**Type / signature consistency:**
- `surfacewater_init(wls1, wlp1)` — same signature as `rddre(wls1, wlp1)`. Caller is a single line in `surfacewater.f90:56`.
- `surface_water_config_t` field names (`hbweir`, `wldip`, `swman`, `wscap`, `nmper`, `wlact`) — used identically in Tasks 4, 5, 6, 7, 8. Verified.
- `drainage_config_t%altcu`, `%nrlevs`, `%zbotdr` — used in Tasks 6, 7, 8.
- New global `wls1_init` (real scalar) added to `variables` module in Task 8 step 1a; written by adapter in step 1b; read by `surfacewater_init` in step 1c.

**Coordinate system:** Task 7 enforces `drain.altcu = 0` (validator). All TOML cases today set `altcu = 0`. The legacy rddre subtraction of `altcu` from `zbotdr`/`hbweir`/`wls1` is therefore a no-op in this port — comparisons and assignments use authored values directly. Future cases needing `altcu /= 0` will require extending the adapter to subtract `altcu` from each of the three globals; this is explicitly out of scope here.

**Submodule discipline:** Tasks 3 and 11 each have explicit inner-commit + outer-bump steps. Task 10 verifies the rename round-trips cleanly without leaving the submodule dirty.
