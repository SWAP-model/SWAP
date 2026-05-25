# Port `swap.ini` to TOML + CSV Companions — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the legacy ASCII `swap.ini` reader for the salinitystress regression case with a structured `[soil.initial]` TOML sub-section (scalars) plus 1–3 separate CSV companion files (z-indexed profile tables) read via the unified `read_csv_table`.

**Architecture:** Schema-first additive migration with three commits in transitional order: (1) add the new `soil_initial_t` type, validator, and TOML reader without removing anything; (2) author the new salinitystress case data while keeping `swap.ini` and `inifil` in place so the legacy adapter path stays green; (3) rewrite the adapter to prefer the new path; (4) atomic cleanup of `inifil` from schema, reader, case data, adapter, and the `swap.ini` file. At every commit boundary the build is clean and `pixi run -e test check-full` is 5/5 green.

**Tech Stack:** Fortran 2008 (gfortran 13), pFUnit, tomlf, pixi. All verification via `pixi run -e test build-linux`, `pixi run -e test test-pfunit`, `pixi run -e test check-full`.

**Spec:** `docs/superpowers/specs/2026-05-01-swap-ini-port-design.md`

---

## File Structure

| File | Role |
|---|---|
| `src/config/soil_config.f90` | Add `soil_initial_t` type with 7 scalar fields, 1 fixed-7 real array, 3 file-path strings; add a range validator. Final task removes the legacy `inifil` slot. |
| `src/config/swap_config.f90` | Add cross-section gating call in `swap_config_t%validate` that fires `ERR_VALIDATION_REQUIRED` per the 3 file slots when their respective gating conditions are met. |
| `src/io/toml/read_soil_toml.f90` | Read the `[soil.initial]` sub-table: 7 scalar fields, `atmin7` inline TOML array (length 7), 3 string slots. Final task removes the legacy `inifil` read. |
| `src/io/toml/config_to_variables.f90` | Replace lines 619–650 with: scalar copy + 3 conditional `read_csv_table` calls when the new schema slot is populated. Uses the existing `block` scoping pattern from the daily/rain/detail meteo CSV pre-loads. |
| `tests/unit/config/test_soil_config.pf` | Tests for the new fields, 7 range tests, and 3 cross-section gating tests (positive + negative). |
| `tests/unit/io/toml/test_read_soil_toml.pf` | Test that loads a `[soil.initial]` fixture and asserts all 11 fields populate correctly. |
| `tests/unit/io/toml/fixtures/soil_initial_full.toml` | TOML fixture covering every field of the new sub-section. |
| `tests/swap-cases/toml/5.salinitystress/swap.toml` | Replace `[soil].inifil = "swap.ini"` with the new `[soil.initial]` block (in cleanup task). |
| `tests/swap-cases/toml/5.salinitystress/salinitystress.ini.h.csv` | New file. Header `z,h`. 178 rows extracted from `swap.ini`. |
| `tests/swap-cases/toml/5.salinitystress/salinitystress.ini.tsoil.csv` | New file. Header `z,tsoil`. 178 rows. |
| `tests/swap-cases/toml/5.salinitystress/salinitystress.ini.cml.csv` | New file. Header `z,cml`. 178 rows. |
| `tests/swap-cases/toml/5.salinitystress/swap.ini` | **Deleted from submodule** in cleanup task. |
| `docs/csv-companion-files.md` | Document 3 new slots in the schema table; mention `swap.ini` removal in migration history. |
| `docs/configuration-schema.md` | Document the `[soil.initial]` schema if a `[soil]` section already exists in this doc. |

---

## Conventions (read before starting)

1. **Build after every source change:** `pixi run -e test build-linux 2>&1 | grep -E "error:" | head -5`
2. **Run targeted pFUnit filter:** `pixi run -e test test-pfunit -- --filter <test_name>` (note: this codebase's `--filter` is forwarded to meson, not pFUnit — use the full `pixi run -e test test-pfunit` and grep the testlog if needed)
3. **Submodule discipline:** Tasks 4 and 6 edit the `tests/swap-cases` submodule. Each submodule commit pairs with a single outer-repo bump.
4. **No physics edits:** changes are limited to `src/config/`, `src/io/toml/`, and the test/case data tree.
5. **Pre-existing dirty submodule files:** `tests/swap-cases` has unrelated modified `swap.toml` / untracked `*.csv` files in the working tree. Use narrow `git add <path>` only — never `git add -A` or `git add .` inside the submodule.
6. **Cross-section validation pattern:** the new validator gating that depends on `swhea`/`swcalt`/`swsolu` lives in `swap_config_t%validate` (which has visibility into all sections), not in `soil_initial_t%validate` (which only does range checks).

---

## Task 1: Add `soil_initial_t` type and range validator

**Files:**
- Modify: `src/config/soil_config.f90`
- Modify: `tests/unit/config/test_soil_config.pf`

**Goal:** Add a new nested type for `[soil.initial]` with all 11 fields and a self-contained range validator. The legacy `inifil` slot stays in place. After this task, the schema has both, but nothing else uses the new slot yet.

- [ ] **Step 1: Write failing test for field round-trip**

Append to `tests/unit/config/test_soil_config.pf`:

```fortran
@test
subroutine test_soil_initial_fields_roundtrip()
   use funit
   use soil_config_mod
   type(soil_config_t) :: c
   c%initial%swirrigate = 1
   c%initial%ssnow      = 0.5d0
   c%initial%slw        = 0.1d0
   c%initial%pond       = 0.2d0
   c%initial%ldwet      = 7.0d0
   c%initial%dt         = 1.0d-7
   c%initial%atmin7     = [1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0, 7.0d0]
   c%initial%h_file     = "case.ini.h.csv"
   c%initial%tsoil_file = "case.ini.tsoil.csv"
   c%initial%cml_file   = "case.ini.cml.csv"
   @assertEqual(1,                 c%initial%swirrigate)
   @assertEqual(0.5d0,             c%initial%ssnow,  1.0d-12)
   @assertEqual(7.0d0,             c%initial%atmin7(7), 1.0d-12)
   @assertEqual("case.ini.h.csv",  c%initial%h_file)
end subroutine
```

- [ ] **Step 2: Run — expect FAIL (compile error: `initial` not a member)**

```bash
pixi run -e test build-linux 2>&1 | grep -E "error:" | head -3
```

Expected: `Error: 'initial' at (1) is not a member of the 'soil_config_t' structure`.

- [ ] **Step 3: Add the new type and field**

Edit `src/config/soil_config.f90`. **Before** the `type :: soil_config_t` declaration (currently around line 58), insert:

```fortran
   !> Initial-state inputs consumed when soil.swinco == 3.
   !! Replaces the legacy ASCII swap.ini file. Scalars copy directly to
   !! globals (ssnow, slw, pond, pondini, ldwet, dt, atmin7); the three
   !! z-indexed profiles are read via read_csv_table as separate companion
   !! CSVs. swirrigate is metadata only — no TOML-side global consumer.
   type :: soil_initial_t
      integer       :: swirrigate = 0
      real(real64)  :: ssnow  = 0.0_real64
      real(real64)  :: slw    = 0.0_real64
      real(real64)  :: pond   = 0.0_real64
      real(real64)  :: ldwet  = 0.0_real64
      real(real64)  :: dt     = 0.0_real64
      real(real64)  :: atmin7(7) = 0.0_real64
      character(len=:), allocatable :: h_file      !! header z,h
      character(len=:), allocatable :: tsoil_file  !! header z,tsoil
      character(len=:), allocatable :: cml_file    !! header z,cml
   contains
      procedure :: validate => soil_initial_validate
   end type soil_initial_t

```

Also add `public :: soil_initial_t` near line 13 (after the existing `public ::` lines).

In the `soil_config_t` body (currently ending at line 112 with `end type soil_config_t`), add a new member just before `contains` (which is on line 109):

```fortran
      type(soil_initial_t) :: initial
```

- [ ] **Step 4: Add the range validator implementation**

Append to the `contains` section (after `soil_frost_validate` ends, around line 264):

```fortran
   subroutine soil_initial_validate(self, errors)
      class(soil_initial_t),    intent(in)    :: self
      type(error_collection_t), intent(inout) :: errors
      integer :: i

      call check_int_enum(self%swirrigate, [0, 1], 'soil.initial.swirrigate', errors)
      call check_real_range(self%ssnow, 0.0_real64, 1000.0_real64, &
                            'soil.initial.ssnow', errors)
      call check_real_range(self%slw,   0.0_real64, 1000.0_real64, &
                            'soil.initial.slw',   errors)
      call check_real_range(self%pond,  0.0_real64,  100.0_real64, &
                            'soil.initial.pond',  errors)
      call check_real_range(self%ldwet, 0.0_real64,  366.0_real64, &
                            'soil.initial.ldwet', errors)
      call check_real_range(self%dt,    1.0e-12_real64, 1.0_real64, &
                            'soil.initial.dt',    errors)
      do i = 1, 7
         call check_real_range(self%atmin7(i), -50.0_real64, 50.0_real64, &
                               'soil.initial.atmin7', errors)
      end do
   end subroutine soil_initial_validate
```

Also call it from `soil_config_validate`. Find the existing block (around line 143):

```fortran
      call self%discretization%validate(errors)
      call self%frost%validate(errors)
      call self%hydraulics%validate(errors)
   end subroutine soil_config_validate
```

Add one line:

```fortran
      call self%discretization%validate(errors)
      call self%frost%validate(errors)
      call self%hydraulics%validate(errors)
      call self%initial%validate(errors)
   end subroutine soil_config_validate
```

- [ ] **Step 5: Run — expect PASS**

```bash
pixi run -e test build-linux 2>&1 | grep -E "error:" | head -3
pixi run -e test test-pfunit 2>&1 | tail -5
```

Expected: clean build, `Ok: 1, Fail: 0`.

- [ ] **Step 6: Add range-validation tests**

Append to `tests/unit/config/test_soil_config.pf`:

```fortran
@test
subroutine test_soil_initial_dt_too_small_fails()
   use funit
   use soil_config_mod
   use error_mod
   type(soil_config_t)        :: c
   type(error_collection_t)   :: errors
   c%initial%dt = 1.0d-15  ! below 1e-12 lower bound
   call c%initial%validate(errors)
   @assertTrue(errors%has_errors())
end subroutine

@test
subroutine test_soil_initial_atmin7_out_of_range_fails()
   use funit
   use soil_config_mod
   use error_mod
   type(soil_config_t)        :: c
   type(error_collection_t)   :: errors
   c%initial%atmin7(3) = 100.0d0   ! above 50.0 upper bound
   call c%initial%validate(errors)
   @assertTrue(errors%has_errors())
end subroutine

@test
subroutine test_soil_initial_zero_defaults_valid()
   use funit
   use soil_config_mod
   use error_mod
   type(soil_config_t)        :: c
   type(error_collection_t)   :: errors
   ! All defaults; nothing set explicitly. Defaults must be valid.
   call c%initial%validate(errors)
   @assertFalse(errors%has_errors())
end subroutine
```

- [ ] **Step 7: Run — expect PASS**

```bash
pixi run -e test test-pfunit 2>&1 | tail -5
```

Expected: `Ok: 1, Fail: 0`.

- [ ] **Step 8: Commit**

```bash
git add src/config/soil_config.f90 tests/unit/config/test_soil_config.pf
git commit -m "feat(soil-config): add soil_initial_t for [soil.initial] schema

Adds 7 scalars (swirrigate, ssnow, slw, pond, ldwet, dt, atmin7) plus
3 file-path slots (h_file, tsoil_file, cml_file) for the new sub-section.
Self-contained range validator wired into soil_config_validate. The
legacy [soil].inifil slot is preserved for now; cleanup at the end of
the swap.ini port plan."
```

---

## Task 2: Cross-section validator gating

**Files:**
- Modify: `src/config/swap_config.f90`
- Modify: `tests/unit/config/test_soil_config.pf` (add cross-section tests via `swap_config_t`)

**Goal:** Add gating rules at `swap_config_t%validate` that fire `ERR_VALIDATION_REQUIRED` when `swinco=3` and the matching `*_file` slot is empty (3 rules). Ranges are already covered by Task 1.

- [ ] **Step 1: Locate the validator and its imports**

```bash
grep -n "soil%validate\|use error_mod\|ERR_VALIDATION" /home/zawadzkim/Code/swap/src/config/swap_config.f90 | head -10
```

Note the line of `call self%soil%validate(errors)` — this is the insertion point.

- [ ] **Step 2: Write failing test — `swinco=3` requires `h_file`**

Append to `tests/unit/config/test_soil_config.pf`:

```fortran
@test
subroutine test_swap_config_swinco3_missing_h_file_fails()
   use funit
   use swap_config_mod
   use error_mod, only: error_collection_t, ERR_VALIDATION_REQUIRED
   type(swap_config_t)        :: c
   type(error_collection_t)   :: errors
   integer :: i
   logical :: found
   c%soil%swinco = 3
   ! h_file is unallocated by default → required-file error must fire.
   call c%validate(errors)
   @assertTrue(errors%has_errors())
   found = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_REQUIRED) found = .true.
   end do
   @assertTrue(found)
end subroutine

@test
subroutine test_swap_config_swinco3_with_h_file_passes_init_check()
   use funit
   use swap_config_mod
   use error_mod, only: error_collection_t, ERR_VALIDATION_REQUIRED
   type(swap_config_t)        :: c
   type(error_collection_t)   :: errors
   integer :: i
   logical :: found_required
   c%soil%swinco = 3
   c%soil%initial%h_file = 'case.ini.h.csv'
   ! Other validators may emit unrelated errors (e.g. missing crops);
   ! we only assert that NO ERR_VALIDATION_REQUIRED for soil.initial.h_file fires.
   call c%validate(errors)
   found_required = .false.
   do i = 1, errors%count()
      if (errors%items(i)%code == ERR_VALIDATION_REQUIRED .and. &
          index(errors%items(i)%message, 'soil.initial.h_file') > 0) then
         found_required = .true.
      end if
   end do
   @assertFalse(found_required)
end subroutine
```

- [ ] **Step 3: Run — expect FAIL (no required-error currently fires)**

```bash
pixi run -e test test-pfunit 2>&1 | tail -10
```

Expected: `Fail: > 0` mentioning `test_swap_config_swinco3_missing_h_file_fails`.

- [ ] **Step 4: Add the gating rules**

Edit `src/config/swap_config.f90`. Find `call self%soil%validate(errors)` (around line 46) and immediately after it, add:

```fortran
      ! Cross-section gating for [soil.initial] required CSV slots.
      ! Range validation has already run inside soil%validate.
      if (self%soil%swinco == 3) then
         if (.not. allocated(self%soil%initial%h_file) .or. &
             len_trim(self%soil%initial%h_file) == 0) then
            call errors%append(ERR_VALIDATION_REQUIRED, &
               'soil.initial.h_file required when soil.swinco=3', &
               'swap_config')
         end if
         if (self%heat%swhea == 1 .and. self%heat%swcalt == 2) then
            if (.not. allocated(self%soil%initial%tsoil_file) .or. &
                len_trim(self%soil%initial%tsoil_file) == 0) then
               call errors%append(ERR_VALIDATION_REQUIRED, &
                  'soil.initial.tsoil_file required when soil.swinco=3 and ' // &
                  'heat.swhea=1 and heat.swcalt=2', 'swap_config')
            end if
         end if
         if (self%solute%swsolu == 1) then
            if (.not. allocated(self%soil%initial%cml_file) .or. &
                len_trim(self%soil%initial%cml_file) == 0) then
               call errors%append(ERR_VALIDATION_REQUIRED, &
                  'soil.initial.cml_file required when soil.swinco=3 and ' // &
                  'solute.swsolu=1', 'swap_config')
            end if
         end if
      end if
```

If `ERR_VALIDATION_REQUIRED` is not in scope, add it to the existing `use error_mod, only:` line near the top of the module.

- [ ] **Step 5: Run — expect PASS**

```bash
pixi run -e test build-linux 2>&1 | grep -E "error:" | head -3
pixi run -e test test-pfunit 2>&1 | tail -10
```

Expected: clean build, `Ok: 1, Fail: 0`.

- [ ] **Step 6: Commit**

```bash
git add src/config/swap_config.f90 tests/unit/config/test_soil_config.pf
git commit -m "feat(swap-config): cross-section gating for soil.initial CSV slots

When soil.swinco=3, h_file is always required; tsoil_file when
heat.swhea=1+swcalt=2; cml_file when solute.swsolu=1. Range checks
remain inside soil_initial_validate (unchanged)."
```

---

## Task 3: TOML reader for `[soil.initial]`

**Files:**
- Create: `tests/unit/io/toml/fixtures/soil_initial_full.toml`
- Modify: `tests/unit/io/toml/test_read_soil_toml.pf`
- Modify: `src/io/toml/read_soil_toml.f90`

**Goal:** Read the `[soil.initial]` sub-table from `swap.toml`. The legacy `inifil` read at `read_soil_toml.f90:51-53` stays untouched.

- [ ] **Step 1: Create the fixture**

Write `tests/unit/io/toml/fixtures/soil_initial_full.toml` with content:

```toml
[soil]
swinco = 3

[soil.initial]
swirrigate = 1
ssnow      = 0.0
slw        = 0.0
pond       = 0.0
ldwet      = 10.0
dt         = 1.0e-7
atmin7     = [14.1, 13.2, 13.8, 16.0, 18.4, 16.3, 16.6]
h_file     = "salinitystress.ini.h.csv"
tsoil_file = "salinitystress.ini.tsoil.csv"
cml_file   = "salinitystress.ini.cml.csv"
```

- [ ] **Step 2: Write failing test**

Append to `tests/unit/io/toml/test_read_soil_toml.pf`:

```fortran
@test
subroutine test_read_soil_initial_full()
   use funit
   use tomlf, only: toml_table, toml_load
   use soil_config_mod
   use read_soil_toml_mod
   use error_mod
   type(toml_table), allocatable, target :: doc
   type(toml_table), pointer             :: doc_ptr
   type(soil_config_t)                   :: c
   type(error_collection_t)              :: errors

   call toml_load(doc, 'tests/unit/io/toml/fixtures/soil_initial_full.toml')
   doc_ptr => doc
   call read_soil_toml(doc_ptr, c, errors)

   @assertFalse(errors%has_errors())
   @assertEqual(1,        c%initial%swirrigate)
   @assertEqual(0.0d0,    c%initial%ssnow,  1.0d-12)
   @assertEqual(10.0d0,   c%initial%ldwet,  1.0d-12)
   @assertEqual(1.0d-7,   c%initial%dt,     1.0d-15)
   @assertEqual(14.1d0,   c%initial%atmin7(1), 1.0d-6)
   @assertEqual(16.6d0,   c%initial%atmin7(7), 1.0d-6)
   @assertTrue(allocated(c%initial%h_file))
   @assertEqual('salinitystress.ini.h.csv',     c%initial%h_file)
   @assertEqual('salinitystress.ini.tsoil.csv', c%initial%tsoil_file)
   @assertEqual('salinitystress.ini.cml.csv',   c%initial%cml_file)
end subroutine
```

- [ ] **Step 3: Run — expect FAIL (reader doesn't populate the new fields yet)**

```bash
pixi run -e test test-pfunit 2>&1 | tail -10
```

Expected: failures on `test_read_soil_initial_full` (fields unallocated / zero).

- [ ] **Step 4: Add the reader logic**

Edit `src/io/toml/read_soil_toml.f90`. Identify what helpers are already imported:

```bash
grep -n "use toml_field_helpers_mod" /home/zawadzkim/Code/swap/src/io/toml/read_soil_toml.f90
```

We need `get_table`, `get_optional_string_with_default`, `get_optional_int_with_default`, `get_optional_real_with_default`. Most of these are already imported (the existing reader already uses them). Verify by grepping the existing reader.

We also need to read a fixed-length real array. Check if `toml_field_helpers_mod` exposes such a helper:

```bash
grep -n "get_optional_real_array\|read_real_array\|toml_array" /home/zawadzkim/Code/swap/src/io/toml/toml_field_helpers.f90 | head -10
```

If a `get_optional_real_array_with_default(table, key, dest_array, expected_len, default, key_path, errors)` helper exists, use it. Otherwise, read the array inline using `tomlf` directly:

```fortran
block
   use tomlf, only: toml_array, get_value, len
   type(toml_array), pointer :: arr_ptr
   integer :: ilen, k
   real(real64) :: tmp
   call get_value(initial_tbl, 'atmin7', arr_ptr, requested=.false.)
   if (associated(arr_ptr)) then
      ilen = min(len(arr_ptr), 7)
      do k = 1, ilen
         call get_value(arr_ptr, k, tmp)
         config%initial%atmin7(k) = tmp
      end do
   end if
end block
```

The exact API for `tomlf`'s array access matches the pattern used elsewhere in `read_*_toml.f90` for inline arrays — search for one that already does this:

```bash
grep -rn "toml_array\|get_value.*array" /home/zawadzkim/Code/swap/src/io/toml/ | head -10
```

In `read_soil_toml.f90`, after the existing `inifil` read (around line 53), add the `[soil.initial]` reader block. **Find** the existing flow (the reader probably uses `call get_table(doc, 'soil', sec, ...)` to obtain the `[soil]` sub-table). Add:

```fortran
      ! [soil.initial] sub-table — typed home for what swap.ini used to carry.
      block
         type(toml_table), pointer :: ini_tbl
         call get_table(sec, 'initial', ini_tbl, 'soil.initial', errors)
         if (associated(ini_tbl)) then
            call get_optional_int_with_default(ini_tbl,    'swirrigate', config%initial%swirrigate, 0, &
                                               'soil.initial.swirrigate', errors)
            call get_optional_real_with_default(ini_tbl,   'ssnow',  config%initial%ssnow,  0.0d0, &
                                                'soil.initial.ssnow',  errors)
            call get_optional_real_with_default(ini_tbl,   'slw',    config%initial%slw,    0.0d0, &
                                                'soil.initial.slw',    errors)
            call get_optional_real_with_default(ini_tbl,   'pond',   config%initial%pond,   0.0d0, &
                                                'soil.initial.pond',   errors)
            call get_optional_real_with_default(ini_tbl,   'ldwet',  config%initial%ldwet,  0.0d0, &
                                                'soil.initial.ldwet',  errors)
            call get_optional_real_with_default(ini_tbl,   'dt',     config%initial%dt,     0.0d0, &
                                                'soil.initial.dt',     errors)
            call get_optional_string_with_default(ini_tbl, 'h_file',     config%initial%h_file,     '', &
                                                  'soil.initial.h_file',     errors)
            call get_optional_string_with_default(ini_tbl, 'tsoil_file', config%initial%tsoil_file, '', &
                                                  'soil.initial.tsoil_file', errors)
            call get_optional_string_with_default(ini_tbl, 'cml_file',   config%initial%cml_file,   '', &
                                                  'soil.initial.cml_file',   errors)

            ! atmin7: fixed-length 7-element inline array. Read directly via tomlf
            ! since there's no fixed-array helper. Fewer than 7 entries → leave
            ! defaults for the missing slots.
            block
               use tomlf, only: toml_array, get_value, len
               type(toml_array), pointer :: arr_ptr
               integer :: ilen, k
               real(real64) :: tmp
               call get_value(ini_tbl, 'atmin7', arr_ptr, requested=.false.)
               if (associated(arr_ptr)) then
                  ilen = min(len(arr_ptr), 7)
                  do k = 1, ilen
                     call get_value(arr_ptr, k, tmp)
                     config%initial%atmin7(k) = tmp
                  end do
               end if
            end block
         end if
      end block
```

Confirm the pattern matches other `read_*_toml.f90` files in this repo — if any TOML-array reads already use a different idiom, follow that instead.

- [ ] **Step 5: Run — expect PASS**

```bash
pixi run -e test build-linux 2>&1 | grep -E "error:" | head -3
pixi run -e test test-pfunit 2>&1 | tail -10
```

Expected: clean build, `Ok: 1, Fail: 0`.

- [ ] **Step 6: Commit**

```bash
git add src/io/toml/read_soil_toml.f90 \
        tests/unit/io/toml/fixtures/soil_initial_full.toml \
        tests/unit/io/toml/test_read_soil_toml.pf
git commit -m "feat(soil-reader): read [soil.initial] sub-table

Reads 7 scalars + atmin7(7) inline TOML array + 3 file slots into
soil_config_t.initial. Existing [soil].inifil read is unchanged
(removed in cleanup task)."
```

---

## Task 4: Migrate salinitystress case data (additive — keep `inifil` and `swap.ini` for now)

**Files (submodule):**
- Create: `tests/swap-cases/toml/5.salinitystress/salinitystress.ini.h.csv`
- Create: `tests/swap-cases/toml/5.salinitystress/salinitystress.ini.tsoil.csv`
- Create: `tests/swap-cases/toml/5.salinitystress/salinitystress.ini.cml.csv`
- Modify: `tests/swap-cases/toml/5.salinitystress/swap.toml` — **add** `[soil.initial]` block; **keep** `[soil].inifil = "swap.ini"` for now.

**Files (outer):**
- Bump submodule pointer.

**Goal:** Author the new case data without disturbing the legacy path. After this task, the case has both `inifil = "swap.ini"` and the new `[soil.initial]` block. Adapter (still legacy in this commit) reads via `inifil`/`swap.ini`. Regression remains green.

- [ ] **Step 1: Extract h profile to CSV**

The `swap.ini` file at `tests/swap-cases/toml/5.salinitystress/swap.ini` has a `*  Soil water pressure heads` section followed by lines `z_h,h` (header) and 178 data rows. Extract them to a clean CSV.

```bash
cd /home/zawadzkim/Code/swap/tests/swap-cases/toml/5.salinitystress
awk '
  /^\*[[:space:]]*Soil water pressure heads/  { in_block = 1; next }
  /^\*/                                       { in_block = 0 }
  in_block && /^[[:space:]]*z_h,h[[:space:]]*$/ { print "z,h"; next }
  in_block && NF                               { print $0 }
' swap.ini | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' > salinitystress.ini.h.csv
wc -l salinitystress.ini.h.csv
```

Expected: 179 lines (1 header + 178 data rows).

Verify a few rows match:

```bash
head -5 salinitystress.ini.h.csv
tail -3 salinitystress.ini.h.csv
```

Expected first 5: `z,h`, `-0.5,-6.01193E+01`, `-1.5,-5.90627E+01`, `-2.5,-5.80122E+01`, `-3.5,-5.69674E+01`.
Expected last 3: `-975.0, 9.13815E+02`, `-985.0, 9.23846E+02`, `-995.0, 9.33877E+02` (or similar — confirm against `swap.ini`).

- [ ] **Step 2: Extract Tsoil profile**

```bash
awk '
  /^\*[[:space:]]*Soil temperatures/          { in_block = 1; next }
  /^\*/                                       { in_block = 0 }
  in_block && /^[[:space:]]*z_Tsoil,Tsoil[[:space:]]*$/ { print "z,tsoil"; next }
  in_block && NF                               { print $0 }
' swap.ini | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' > salinitystress.ini.tsoil.csv
wc -l salinitystress.ini.tsoil.csv
```

Expected: 179 lines.

- [ ] **Step 3: Extract Cml profile**

```bash
awk '
  /^\*[[:space:]]*Solute concentrations/      { in_block = 1; next }
  /^\*/                                       { in_block = 0 }
  in_block && /^[[:space:]]*z_Cml,Cml[[:space:]]*$/ { print "z,cml"; next }
  in_block && NF                               { print $0 }
' swap.ini | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' > salinitystress.ini.cml.csv
wc -l salinitystress.ini.cml.csv
```

Expected: 179 lines.

Note: the awk patterns above have to match the exact `*  Soil water pressure heads`, `*  Soil temperatures`, `*  Solute concentrations` headings that introduce each block. If the awk extraction fails (wrong line counts), inspect `swap.ini` and adjust the regex. The expected line count is **178 data rows + 1 header = 179 lines** for each CSV; the existing `swap.ini` has identical-length z grids for all three blocks.

- [ ] **Step 4: Add `[soil.initial]` block to swap.toml**

In `tests/swap-cases/toml/5.salinitystress/swap.toml`, locate the `[soil]` section. **Keep** the `inifil = "swap.ini"` line for now. Add a new `[soil.initial]` block near the end of the file (or anywhere — TOML order doesn't matter for parsing, but keep it near `[soil]` for human readers):

```toml
[soil.initial]
swirrigate = 1
ssnow      = 0.0
slw        = 0.0
pond       = 0.0
ldwet      = 10.0
dt         = 1.0e-7
atmin7     = [14.1, 13.2, 13.8, 16.0, 18.4, 16.3, 16.6]
h_file     = "salinitystress.ini.h.csv"
tsoil_file = "salinitystress.ini.tsoil.csv"
cml_file   = "salinitystress.ini.cml.csv"
```

The scalar values come from `swap.ini`:
- `swirrigate = 1` (line `SWIRRIGATE = 1`)
- `ssnow = 0.0` (line `Ssnow =  0.00000E+00`)
- `slw = 0.0`
- `pond = 0.0`
- `ldwet = 10.0` (line `ldwet =  0.10000E+02`)
- `dt = 1.0e-7` (line `dt =  0.10000E-06`)
- `atmin7 = [14.1, 13.2, 13.8, 16.0, 18.4, 16.3, 16.6]` (the 7-row `atmin7` block)

- [ ] **Step 5: Submodule commit**

```bash
cd /home/zawadzkim/Code/swap/tests/swap-cases
git add toml/5.salinitystress/salinitystress.ini.h.csv \
        toml/5.salinitystress/salinitystress.ini.tsoil.csv \
        toml/5.salinitystress/salinitystress.ini.cml.csv \
        toml/5.salinitystress/swap.toml
git status --short
git commit -m "test(salinitystress): add [soil.initial] block + 3 init CSV companions

Extracted from swap.ini (currently still required by the legacy adapter
path). Once the adapter switches to the new schema, swap.ini and the
old [soil].inifil slot will be removed."
```

- [ ] **Step 6: Outer-repo bump**

```bash
cd /home/zawadzkim/Code/swap
git add tests/swap-cases
git commit -m "chore(swap-cases): bump for [soil.initial] migration (additive)"
```

- [ ] **Step 7: Verify regression still green (legacy adapter path still active)**

```bash
pixi run -e test check-full 2>&1 | tail -10
```

Expected: 5 passed, 0 failed. Salinitystress still uses the legacy path (`inifil = "swap.ini"`) since the adapter rewrite hasn't happened yet.

---

## Task 5: Adapter rewrite — prefer the new path

**Files:**
- Modify: `src/io/toml/config_to_variables.f90`

**Goal:** Replace the existing `swinco==3 && allocated(inifil)` block (lines 619–650) with a two-arm conditional: if the new `h_file` slot is populated, use the CSV path; otherwise fall through to the legacy ASCII path. After this task, salinitystress runs via the new path; the legacy path is dormant but intact.

- [ ] **Step 1: Locate the existing block**

```bash
grep -n "swinco == 3\|inifil\|rdsdor.*ssnow" /home/zawadzkim/Code/swap/src/io/toml/config_to_variables.f90 | head -10
```

The block currently spans roughly lines 619–650 (`if (config%soil%swinco == 3 .and. allocated(config%soil%inifil)) then` through the matching `end if`).

- [ ] **Step 2: Replace the block**

Find the exact block (starting with the comment `! readswap.f90:1593-1626`...) and replace it with:

```fortran
      ! [soil.initial] / inifil dual-path block. The new path uses the
      ! typed soil.initial schema + per-profile CSV companions; the
      ! legacy path opens swap.ini via TTutil. Both consume swinco=3.
      if (config%soil%swinco == 3) then
         if (allocated(config%soil%initial%h_file) .and. &
             len_trim(config%soil%initial%h_file) > 0) then
            ! ----- NEW CSV-based path -----------------------------------
            ssnow   = config%soil%initial%ssnow
            slw     = config%soil%initial%slw
            pond    = config%soil%initial%pond
            pondini = pond
            ldwet   = config%soil%initial%ldwet
            dt      = config%soil%initial%dt
            atmin7(:) = config%soil%initial%atmin7(:)
            ! Legacy zeroes ssnow when swsnow != 1 (readswap.f90:1606-1613).
            if (config%meteo%snow%swsnow /= 1) ssnow = 0.0d0

            ! Mandatory: initial pressure-head profile (z, h).
            block
               use csv_reader_mod,  only: read_csv_table
               use error_mod,       only: error_collection_t
               real(8), allocatable     :: tbl(:,:)
               type(error_collection_t) :: errs
               character(len=2)         :: hdr(2)
               integer :: nrows, k
               hdr(1) = 'z '
               hdr(2) = 'h '
               call read_csv_table(trim(config%soil%initial%h_file), hdr, tbl, errs)
               call errs%abort_if_fatal()
               nrows = size(tbl, 1)
               nhead = nrows
               do k = 1, nrows
                  zi(k) = tbl(k, 1)
                  h(k)  = tbl(k, 2)
               end do
            end block

            ! Optional: initial soil temperature profile.
            if (config%heat%swhea == 1 .and. config%heat%swcalt == 2) then
               block
                  use csv_reader_mod,  only: read_csv_table
                  use error_mod,       only: error_collection_t
                  real(8), allocatable     :: tbl(:,:)
                  type(error_collection_t) :: errs
                  character(len=5)         :: hdr(2)
                  integer :: nrows, k
                  hdr(1) = 'z    '
                  hdr(2) = 'tsoil'
                  call read_csv_table(trim(config%soil%initial%tsoil_file), hdr, tbl, errs)
                  call errs%abort_if_fatal()
                  nrows = size(tbl, 1)
                  do k = 1, nrows
                     zh(k)    = tbl(k, 1)
                     tsoil(k) = tbl(k, 2)
                  end do
               end block
            end if

            ! Optional: initial concentration profile (Cml).
            if (config%solute%swsolu == 1) then
               block
                  use csv_reader_mod,  only: read_csv_table
                  use error_mod,       only: error_collection_t
                  real(8), allocatable     :: tbl(:,:)
                  type(error_collection_t) :: errs
                  character(len=3)         :: hdr(2)
                  integer :: nrows, k
                  hdr(1) = 'z  '
                  hdr(2) = 'cml'
                  call read_csv_table(trim(config%soil%initial%cml_file), hdr, tbl, errs)
                  call errs%abort_if_fatal()
                  nrows = size(tbl, 1)
                  nconc = nrows
                  do k = 1, nrows
                     zc(k)  = tbl(k, 1)
                     cml(k) = tbl(k, 2)
                  end do
               end block
            end if
            ! ----- end NEW path -----------------------------------------

         else if (allocated(config%soil%inifil) .and. &
                  len_trim(config%soil%inifil) > 0) then
            ! ----- LEGACY ASCII swap.ini path (preserved for now) -------
            block
               use swap_array_dimensions, only: macp
               integer :: ini_unit, ifnd_ini
               character(len=200) :: ini_filnam
               integer, external :: getun2
               ini_filnam = trim(config%soil%inifil)
               ini_unit = getun2(10, 90, 2)
               call rdinit(ini_unit, logf, ini_filnam)
               call rdsdor('ssnow', 0.0d0, 1000.0d0, ssnow)
               if (config%meteo%snow%swsnow /= 1) ssnow = 0.0d0
               call rdsdor('slw',   0.0d0, 1000.0d0, slw)
               call rdsdor('pond',  0.0d0,  100.0d0, pond)
               pondini = pond
               call rdador('z_h',  -1.0d5,  0.0d0, zi, macp, ifnd_ini)
               call rdfdor('h',    -1.0d10, 1.0d4, h,  macp, ifnd_ini)
               nhead = ifnd_ini
               if (config%heat%swhea == 1 .and. config%heat%swcalt == 2) then
                  call rdador('z_Tsoil', -1.0d5,  0.0d0, zh,    macp, ifnd_ini)
                  call rdfdor('Tsoil',  -50.0d0, 50.0d0, tsoil, macp, ifnd_ini)
               end if
               if (config%solute%swsolu == 1) then
                  call rdador('z_Cml', -1.0d5,    0.0d0, zc,  macp, ifnd_ini)
                  call rdfdor('Cml',    0.0d0, 1.0d6,    cml, macp, ifnd_ini)
                  nconc = ifnd_ini
               end if
               close(ini_unit)
            end block
            ! ----- end LEGACY path --------------------------------------
         end if
      end if
```

The two arms are complete: the new path covers everything the legacy reader does (scalars, h profile, optional tsoil profile, optional cml profile) plus the new fields (`ldwet`, `dt`, `atmin7`).

- [ ] **Step 3: Build**

```bash
pixi run -e test build-linux 2>&1 | grep -E "error:" | head -5
```

Expected: clean.

- [ ] **Step 4: Run regression — salinitystress now uses the new path**

```bash
pixi run -e test check-full 2>&1 | tail -10
```

Expected: 5 passed, 0 failed. Salinitystress: regression ok at 1e-2 cm tolerance. The change is purely about *where* the initial-state values come from; the values themselves are identical (CSVs are direct extractions from `swap.ini`).

If salinitystress diverges, the most likely cause is a fixture extraction error — re-check the CSV row count and a few values against `swap.ini`.

- [ ] **Step 5: Commit**

```bash
git add src/io/toml/config_to_variables.f90
git commit -m "refactor(adapter): prefer [soil.initial] CSV path for swinco=3

Two-arm conditional: when soil.initial.h_file is set, read scalars from
the typed schema and profiles from CSV companions via read_csv_table;
otherwise fall through to the legacy swap.ini reader (preserved for
case-by-case migration). Salinitystress now exercises the new path."
```

---

## Task 6: Cleanup — remove `inifil`, `swap.ini`, and the legacy adapter arm

**Files (outer repo):**
- Modify: `src/config/soil_config.f90` — remove `inifil` field
- Modify: `src/io/toml/read_soil_toml.f90` — remove `inifil` read
- Modify: `src/io/toml/config_to_variables.f90` — remove the legacy arm of the dual-path block

**Files (submodule):**
- Modify: `tests/swap-cases/toml/5.salinitystress/swap.toml` — remove `inifil = "swap.ini"` line
- Delete: `tests/swap-cases/toml/5.salinitystress/swap.ini`

**Goal:** Atomic cleanup commit that removes everything related to the legacy `inifil` slot. All five edits ship together; the build and regression must remain green.

- [ ] **Step 1: Remove `inifil` from `soil_config_t`**

In `src/config/soil_config.f90`, find the line `character(len=:), allocatable :: inifil` (around line 93) and delete it.

- [ ] **Step 2: Remove `inifil` read from the TOML reader**

In `src/io/toml/read_soil_toml.f90`, find and delete the block:

```fortran
      ! Phase 4f Task B5: SWINCO=3 inifil (path to previous-run state file).
      call get_optional_string_with_default(sec, 'inifil', config%inifil, '', &
                                            'soil.inifil', errors)
```

(3 lines, around line 51–53.)

- [ ] **Step 3: Remove the legacy arm from the adapter**

In `src/io/toml/config_to_variables.f90`, find the block from Task 5 and delete the `else if (allocated(config%soil%inifil) ...)` arm and the `LEGACY` block. The remaining structure is:

```fortran
      if (config%soil%swinco == 3) then
         if (allocated(config%soil%initial%h_file) .and. &
             len_trim(config%soil%initial%h_file) > 0) then
            ! ----- NEW CSV-based path (unchanged from Task 5) -----------
            ! ... full block from Task 5 ...
         end if
      end if
```

The outer comment "[soil.initial] / inifil dual-path block" should be updated to say "[soil.initial] block" (no longer dual-path).

- [ ] **Step 4: Outer-repo build to confirm compilation**

```bash
pixi run -e test build-linux 2>&1 | grep -E "error:" | head -5
```

Expected: clean. If anything still references `config%soil%inifil`, the compiler will catch it.

- [ ] **Step 5: Update salinitystress.toml in the submodule**

```bash
cd /home/zawadzkim/Code/swap/tests/swap-cases/toml/5.salinitystress
```

Edit `swap.toml`: locate the `[soil]` section and remove the line:

```
inifil = "swap.ini"
```

Leave the rest of `[soil]` and the `[soil.initial]` block intact.

- [ ] **Step 6: Delete swap.ini from the submodule**

```bash
cd /home/zawadzkim/Code/swap/tests/swap-cases
git rm toml/5.salinitystress/swap.ini
```

- [ ] **Step 7: Submodule commit**

```bash
cd /home/zawadzkim/Code/swap/tests/swap-cases
git add toml/5.salinitystress/swap.toml
git status --short
git commit -m "test(salinitystress): drop inifil + swap.ini, [soil.initial] is canonical

The legacy swap.ini reader path was removed from the TOML adapter; the
typed [soil.initial] schema is now the only home for initial state."
```

- [ ] **Step 8: Outer-repo bump + cleanup commit**

```bash
cd /home/zawadzkim/Code/swap
git add tests/swap-cases src/config/soil_config.f90 \
        src/io/toml/read_soil_toml.f90 src/io/toml/config_to_variables.f90
git commit -m "chore(soil): drop legacy inifil slot — [soil.initial] is canonical

Atomic cleanup removing the legacy ASCII swap.ini path:
- src/config/soil_config.f90: remove inifil field
- src/io/toml/read_soil_toml.f90: remove inifil read
- src/io/toml/config_to_variables.f90: remove legacy adapter arm
- tests/swap-cases submodule: drop swap.ini + inifil = ... line in toml"
```

- [ ] **Step 9: Verify regression**

```bash
pixi run -e test check-full 2>&1 | tail -10
```

Expected: 5 passed, 0 failed.

- [ ] **Step 10: Verify pFUnit**

```bash
pixi run -e test test-pfunit 2>&1 | tail -10
```

Expected: `Ok: 1, Fail: 0`.

- [ ] **Step 11: Confirm no stragglers**

```bash
git grep -n "inifil\|swap\.ini" -- src/ tests/regression/ 2>/dev/null
```

Expected: zero matches in `src/`. Submodule check:

```bash
git -C tests/swap-cases grep -n "inifil\|swap\.ini" -- toml/ 2>/dev/null
```

Expected: zero matches.

(Matches in legacy `readswap.f90` are expected and out of scope — the `.swp` pathway still uses `inifil`.)

---

## Task 7: Documentation

**Files:**
- Modify: `docs/csv-companion-files.md`
- Modify: `docs/configuration-schema.md` (only if the file has a `[soil]` section)

- [ ] **Step 1: Update `docs/csv-companion-files.md`**

In `docs/csv-companion-files.md`, find the schema table that lists CSV slots (it documents `gwl_file`, `haquif_file`, `events_file`, etc.). Add three rows for the new slots:

| `h_file` | `[soil.initial]` | `z,h` | `swinco = 3` |
| `tsoil_file` | `[soil.initial]` | `z,tsoil` | `swinco = 3` and `heat.swhea = 1` and `heat.swcalt = 2` |
| `cml_file` | `[soil.initial]` | `z,cml` | `swinco = 3` and `solute.swsolu = 1` |

In the "Migration history" section (around line 99), add:

```markdown
- `salinitystress.ini.{h,tsoil,cml}.csv` — replaced legacy ASCII `swap.ini` profile blocks (Phase 4f cleanup, post-CSV-meteo).
```

- [ ] **Step 2: Update `docs/configuration-schema.md` if applicable**

```bash
grep -n "\[soil\]\|\[soil\." /home/zawadzkim/Code/swap/docs/configuration-schema.md 2>/dev/null | head -5
```

If the doc has a `[soil]` section, add a `[soil.initial]` sub-section documenting the 11 fields (use the field table from the spec verbatim). Otherwise, skip this step.

- [ ] **Step 3: Commit**

```bash
git add docs/csv-companion-files.md docs/configuration-schema.md 2>/dev/null
git status --short
git commit -m "docs: document [soil.initial] schema + 3 init CSV companions"
```

If `docs/configuration-schema.md` was untouched (no `[soil]` section there), the commit is just for `csv-companion-files.md`.

---

## Task 8: Final verification

- [ ] **Step 1: Full pFUnit suite**

```bash
pixi run -e test test-pfunit 2>&1 | tail -10
```

Expected: `Ok: 1, Fail: 0`.

- [ ] **Step 2: Full regression**

```bash
pixi run -e test check-full 2>&1 | tail -15
```

Expected:
```
✓ hupselbrook:    regression ok
✓ grassgrowth:    regression ok
✓ oxygenstress:   regression ok
✓ salinitystress: regression ok
✓ surfacewater:   regression ok
Results: 5 passed, 0 failed
```

- [ ] **Step 3: Manual smoke test of salinitystress**

```bash
cd /home/zawadzkim/Code/swap/tests/swap-cases
./run_case.sh -c salinitystress -e ../../builddir/swap
ls toml/5.salinitystress/
```

Expected: case runs to completion ("Swap normal completion!"); after cleanup, the listing shows the three init CSVs (`salinitystress.ini.h.csv`, `salinitystress.ini.tsoil.csv`, `salinitystress.ini.cml.csv`) and **no** `swap.ini`.

- [ ] **Step 4: Final grep — no stragglers**

```bash
cd /home/zawadzkim/Code/swap
git grep -n "inifil" -- src/io/toml/ src/config/ tests/unit/
```

Expected: zero matches. (Legacy `readswap.f90` still references `inifil` for the `.swp` path — that's expected and out of scope.)

- [ ] **Step 5: Tag the milestone**

```bash
git tag rescue/phase-swap-ini-port
git tag --list 'rescue/*' | tail -5
```

---

## Risks and mitigations

| Risk | Mitigation |
|---|---|
| Numeric drift in salinitystress regression after the migration | The CSVs are direct extractions from `swap.ini` data rows — the awk pipeline preserves the original `*` format-spec values byte-for-byte. If drift appears, check the awk extraction for trimmed leading zeros or scientific-notation parsing differences; re-run with a more conservative extractor. |
| Atomic cleanup commit breaks the build at any intermediate point | Task 6 stages all 4 source changes + the submodule edit in one git commit pair. The transitional Task 5 (dual-path) keeps the legacy arm alive until Task 6, so each commit boundary is green. |
| `tomlf` array-read API differs from the pseudocode in Task 3 Step 4 | The `block` content uses `toml_array, get_value, len` from `tomlf` directly; this is the documented API. If the project has a wrapper helper, follow whatever pattern other `read_*_toml.f90` files use for inline arrays. The fallback (untyped `get_value(arr_ptr, k, tmp)`) is portable. |
| `swirrigate` is read but never consumed in the TOML pipeline | Documented in the spec: the field exists for round-trip completeness; the adapter has a comment "config%soil%initial%swirrigate is read but unused — no TOML-side consumer." If a future consumer needs it, the slot is already in the schema. |
| Pre-existing dirty submodule files get accidentally committed | All `git add` commands inside `tests/swap-cases` use explicit paths (`toml/5.salinitystress/...`), never `-A` or `.`. Verify with `git status --short` between steps. |

---

## What this plan does NOT cover (deferred)

- **`swap.dra` port for surfacewater** — separate spec (touches ~30 fields across 4 sections of the legacy file)
- **Crop `.crp` port** — separate spec, larger scope (rotation-aware, multi-crop)
- **Other 4 regression cases** — they use `swinco=0` or `swinco=2`, unaffected by this change
- **Removing TTutil `rdsdor` / `rdador` / `rdfdor` from `readswap.f90`** — the `.swp` pathway still uses these readers in unrelated places
- **Compatibility shim for legacy `inifil`** — none. The slot is removed outright in Task 6
