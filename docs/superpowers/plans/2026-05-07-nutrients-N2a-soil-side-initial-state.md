# [nutrients] N2a — Soil-side initial state + SorpCoef Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a top-level `[nutrients]` TOML block carrying soil-side initial pool concentrations (`FOM_t(1..8)`, `Bio_t`, `Hum_t`, `cNH4_t`, `cNO3_t`) plus the `SorpCoef` sorption coefficient. Wire values to legacy `variables`-module globals via a new `apply_nutrients` adapter. Retire the broken `SoilManagement(7)` `_nut.end` diagnostic dump.

**Architecture:** New module `nutrients_config.f90` defines `nutrients_initial_t` (12 reals) and `nutrients_config_t` (a presence flag + sorp_coef + initial sub-block). New TOML reader `read_nutrients_toml` parses the optional top-level `[nutrients]` block; absent block leaves defaults (`present=.false.`, all pools 0.0, `sorp_coef=0.0`). New adapter `apply_nutrients` in `config_to_variables.f90` runs unconditionally — sets `SorpCoef` deterministically (fixes a genuine uninitialised-variable bug discovered during brainstorm) and copies pool values to the legacy globals declared in `wofost_soil_declarations.f90`. `SoilManagement(7)` body collapses to `return` (the legacy diagnostic dump used `<project>.snp` as a template, which can't exist in the modern flow).

**Tech Stack:** Fortran 2008, meson + ninja, gfortran, pFUnit. Reference patterns: `apply_soil_tillage` (commit `c99ede5`) and `apply_irrigation_ssdi` (commit `328db06`) for adapter shape; `read_drainage_toml` (commit `0c3829a` for tillage's reader) or `read_irrigation_ssdi_toml` (commit `a13cd5b`) for top-level-block reader shape; `wofost_soil_declarations.f90` for the legacy globals' declarations.

**Spec:** `docs/superpowers/specs/2026-05-07-nutrients-N2a-soil-side-initial-state-design.md`.

---

## File Structure

**Created:**

| Path | Responsibility |
|---|---|
| `src/config/nutrients_config.f90` | `nutrients_config_mod` — `nutrients_initial_t` and `nutrients_config_t` types + `validate` + `finalize`. |
| `src/io/toml/read_nutrients_toml.f90` | Parses the optional top-level `[nutrients]` block + `[nutrients.initial]` sub-table. |
| `tests/unit/io/toml/test_nutrients_config.pf` | pFUnit suite — 11 tests covering parser, validator, adapter. |
| `docs/adr/0026-nutrients-N2a-soil-side-initial-state.md` | New ADR. |

**Modified:**

| Path | Change |
|---|---|
| `src/config/swap_config.f90` | Add `type(nutrients_config_t) :: nutrients` field on `swap_config_t`. |
| `src/io/toml/load_swap_config.f90` | Append `call read_nutrients_toml(doc_ptr, config%nutrients, errors)` after the existing `read_*_toml` chain. |
| `src/io/toml/config_to_variables.f90` | Add `apply_nutrients(cfg)` helper; call it unconditionally from the top-level `config_to_variables` adapter. |
| `src/crop/management_soil.f90` | Collapse `case (7)` body to `return` + comment; drop unused locals (`snp`, `oup`, `filnamce`, `cropext`-related) if no longer referenced after the case-7 retirement. |
| `meson.build` | Add `'src/config/nutrients_config.f90'` and `'src/io/toml/read_nutrients_toml.f90'` to production sources. |
| `tests/unit/meson.build` | Add the two new sources to `pfunit_extra_sources`; register `.pf` in `pf_files`. |
| `tests/unit/testSuites.inc` | Add `ADD_TEST_SUITE(test_nutrients_config_suite)`. |
| `docs/adr/index.md` | Append ADR 0026 row. |
| `docs/configuration-schema.md` | Document `[nutrients]` + `[nutrients.initial]` blocks. |

**Pre-flight inventory results (verified at spec-authoring):**

- `load_swap_config.f90` has the `read_*_toml` chain at lines 45-56. The new call appends after `read_output_csv_toml` (which is the last one).
- `swap_config_t` field block in `swap_config.f90` is at lines 24-33; new field appends after `surface_water`.
- `config_to_variables.f90` already has `apply_soil_tillage` (line 488) and `apply_irrigation_ssdi` (line 489) calls; `apply_nutrients` joins the same area.
- `SoilManagement(7)` body spans lines 587-657. It contains two retire-able branches: the `flCropExt` block (lines 589-608, reads/rewrites the `cropext` file) and the `.snp` → `_nut.end` template-write (the bulk of the body). Both go.
- `flCropExt` is declared in `wofost_soil_declarations.f90:119`; `cstring(10000)` at :118; `cropext` (unit-number) and `filnamce` (path) declared in `management_soil.f90:46`. After case-7 collapse, verify which become unreferenced.

---

## Task 1: `nutrients_config_t` types + validator + smoke test

**Files:**
- Create: `src/config/nutrients_config.f90`
- Modify: `src/config/swap_config.f90` (add field)
- Modify: `meson.build` (add new source to production sources)
- Modify: `tests/unit/meson.build` (add to `pfunit_extra_sources`)
- Create: `tests/unit/io/toml/test_nutrients_config.pf` (smoke + validator tests)
- Modify: `tests/unit/meson.build` (register `.pf`)
- Modify: `tests/unit/testSuites.inc` (register suite)

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/io/toml/test_nutrients_config.pf`:

```fortran
@test
subroutine test_nutrients_default_init()
   use nutrients_config_mod, only: nutrients_config_t
   use funit
   implicit none
   type(nutrients_config_t) :: cfg
   integer :: i

   ! Defaults: present=false, sorp_coef=0, all initial pools at 0
   @assertFalse(cfg%present, 'present default false')
   @assertEqual(0.0d0, cfg%sorp_coef, tolerance=1.0d-12)
   do i = 1, 8
      @assertEqual(0.0d0, cfg%initial%fom(i), tolerance=1.0d-12)
   end do
   @assertEqual(0.0d0, cfg%initial%bio,  tolerance=1.0d-12)
   @assertEqual(0.0d0, cfg%initial%hum,  tolerance=1.0d-12)
   @assertEqual(0.0d0, cfg%initial%cnh4, tolerance=1.0d-12)
   @assertEqual(0.0d0, cfg%initial%cno3, tolerance=1.0d-12)
end subroutine

@test
subroutine test_nutrients_validate_silent_when_absent()
   use nutrients_config_mod, only: nutrients_config_t
   use error_mod, only: error_collection_t
   use funit
   implicit none
   type(nutrients_config_t)  :: cfg
   type(error_collection_t)  :: errs

   ! present=.false. by default — validator must stay silent.
   call cfg%validate(errs)
   call assertEqual(0, errs%count(), 'expected zero errors when present=.false.')
end subroutine

@test
subroutine test_nutrients_validate_passes_clean()
   use, intrinsic :: iso_fortran_env, only: real64
   use nutrients_config_mod, only: nutrients_config_t
   use error_mod, only: error_collection_t
   use funit
   implicit none
   type(nutrients_config_t)  :: cfg
   type(error_collection_t)  :: errs

   cfg%present   = .true.
   cfg%sorp_coef = 0.005_real64
   cfg%initial%fom = [0.5_real64, 0.3_real64, 0.2_real64, 0.1_real64, &
                      0.5_real64, 0.3_real64, 0.2_real64, 0.1_real64]
   cfg%initial%bio  = 0.4_real64
   cfg%initial%hum  = 8.0_real64
   cfg%initial%cnh4 = 0.001_real64
   cfg%initial%cno3 = 0.005_real64

   call cfg%validate(errs)
   call assertEqual(0, errs%count(), 'expected zero errors for clean populated config')
end subroutine

@test
subroutine test_nutrients_validate_rejects_negative_sorp_coef()
   use, intrinsic :: iso_fortran_env, only: real64
   use nutrients_config_mod, only: nutrients_config_t
   use error_mod, only: error_collection_t
   use funit
   implicit none
   type(nutrients_config_t)  :: cfg
   type(error_collection_t)  :: errs

   cfg%present   = .true.
   cfg%sorp_coef = -0.1_real64
   call cfg%validate(errs)
   @assertTrue(errs%count() > 0, 'expected error for negative sorp_coef')
end subroutine

@test
subroutine test_nutrients_validate_rejects_fom_above_1000()
   use, intrinsic :: iso_fortran_env, only: real64
   use nutrients_config_mod, only: nutrients_config_t
   use error_mod, only: error_collection_t
   use funit
   implicit none
   type(nutrients_config_t)  :: cfg
   type(error_collection_t)  :: errs

   cfg%present = .true.
   cfg%initial%fom(3) = 1500.0_real64
   call cfg%validate(errs)
   @assertTrue(errs%count() > 0, 'expected error for fom(3) > 1000')
end subroutine

@test
subroutine test_nutrients_validate_rejects_negative_bio()
   use, intrinsic :: iso_fortran_env, only: real64
   use nutrients_config_mod, only: nutrients_config_t
   use error_mod, only: error_collection_t
   use funit
   implicit none
   type(nutrients_config_t)  :: cfg
   type(error_collection_t)  :: errs

   cfg%present     = .true.
   cfg%initial%bio = -0.5_real64
   call cfg%validate(errs)
   @assertTrue(errs%count() > 0, 'expected error for negative bio')
end subroutine
```

- [ ] **Step 2: Wire test + run to verify failures**

Add to `tests/unit/meson.build`'s `pf_files`:
```meson
        'io/toml/test_nutrients_config.pf',
```

Add to `tests/unit/testSuites.inc`:
```
ADD_TEST_SUITE(test_nutrients_config_suite)
```

Run:
```
pixi run -e test build-linux
```
Expected: build fails — `nutrients_config_mod` does not exist.

- [ ] **Step 3: Create `src/config/nutrients_config.f90`**

```fortran
!> @file nutrients_config.f90
!! Top-level [nutrients] config block: soil-side initial pool
!! concentrations + SorpCoef sorption coefficient.
!!
!! N2a of the [nutrients] umbrella (ADR 0026). Block is optional;
!! absent → present=.false., all defaults zero.
module nutrients_config_mod
   use, intrinsic :: iso_fortran_env, only: real64
   use error_mod, only: error_collection_t, ERR_VALIDATION_OUT_OF_RANGE
   use validation_mod, only: check_real_range
   implicit none
   private

   public :: nutrients_initial_t
   public :: nutrients_config_t

   !> Initial concentrations of soil organic-matter and N pools.
   !! Units: kg/m^2 for FOM/Bio/Hum (depth-integrated); kg/m^3 for
   !! cNH4/cNO3. Legacy ranges (rdsdor [0, 1000]) preserved.
   type :: nutrients_initial_t
      real(real64) :: fom(8) = 0.0_real64
      real(real64) :: bio    = 0.0_real64
      real(real64) :: hum    = 0.0_real64
      real(real64) :: cnh4   = 0.0_real64
      real(real64) :: cno3   = 0.0_real64
   end type nutrients_initial_t

   !> [nutrients] top-level block.
   type :: nutrients_config_t
      logical      :: present   = .false.
      real(real64) :: sorp_coef = 0.0_real64
      type(nutrients_initial_t) :: initial
   contains
      procedure :: validate => nutrients_config_validate
      procedure :: finalize => nutrients_config_finalize
   end type nutrients_config_t

contains

   subroutine nutrients_config_validate(self, errors)
      class(nutrients_config_t), intent(in)    :: self
      type(error_collection_t),  intent(inout) :: errors
      integer :: i

      if (.not. self%present) return

      ! sorp_coef must be non-negative; no upper bound (literature
      ! values vary widely; high values may be intentional).
      if (self%sorp_coef < 0.0_real64) then
         call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
            'nutrients.sorp_coef: must be >= 0', &
            'nutrients.sorp_coef')
      end if

      ! Initial pool concentrations: legacy rdsdor [0, 1000] range.
      do i = 1, 8
         call check_real_range(self%initial%fom(i), 0.0_real64, 1000.0_real64, &
                               'nutrients.initial.fom', errors)
      end do
      call check_real_range(self%initial%bio,  0.0_real64, 1000.0_real64, &
                            'nutrients.initial.bio',  errors)
      call check_real_range(self%initial%hum,  0.0_real64, 1000.0_real64, &
                            'nutrients.initial.hum',  errors)
      call check_real_range(self%initial%cnh4, 0.0_real64, 1000.0_real64, &
                            'nutrients.initial.cnh4', errors)
      call check_real_range(self%initial%cno3, 0.0_real64, 1000.0_real64, &
                            'nutrients.initial.cno3', errors)
   end subroutine nutrients_config_validate

   subroutine nutrients_config_finalize(self, errors)
      class(nutrients_config_t), intent(inout) :: self
      type(error_collection_t),  intent(inout) :: errors
      ! No finalization needed; included for interface uniformity.
      return
   end subroutine nutrients_config_finalize

end module nutrients_config_mod
```

- [ ] **Step 4: Add the field to `swap_config_t`**

In `src/config/swap_config.f90`, add to the imports near the top:
```fortran
   use nutrients_config_mod, only: nutrients_config_t
```

In the `type :: swap_config_t` block (around lines 24-33), append:
```fortran
      type(nutrients_config_t)   :: nutrients
```

- [ ] **Step 5: Register sources in meson**

In `meson.build`, add `'src/config/nutrients_config.f90'` to the `sources` list, near the other `src/config/` files (e.g., after `'src/config/surface_water_config.f90'`).

In `tests/unit/meson.build`, add to `pfunit_extra_sources`:
```meson
        '../../src/config/nutrients_config.f90',
```

- [ ] **Step 6: Build and run tests**

```
pixi run -e test build-linux
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:|tests,"
pixi run -e test check-full 2>&1 | grep "Results:"
```
Expected: build clean; pFUnit `Ok: 1, Fail: 0` (test count grows by 6); check-full `5 passed, 0 failed`.

- [ ] **Step 7: Commit**

```bash
git add src/config/nutrients_config.f90 src/config/swap_config.f90 meson.build tests/unit/meson.build tests/unit/testSuites.inc tests/unit/io/toml/test_nutrients_config.pf
git commit -m "$(cat <<'EOF'
feat(config): add nutrients_config_t types + validator

Introduces src/config/nutrients_config.f90 with two types:

- nutrients_initial_t: FOM_t(1..8), Bio_t, Hum_t, cNH4_t, cNO3_t
  (12 reals; legacy rdsdor [0, 1000] range preserved)
- nutrients_config_t: presence flag, sorp_coef, initial sub-block

Block is optional (presence flag). Validator stays silent when
present=.false.. When present=.true., enforces sorp_coef >= 0
and 12-pool [0, 1000] ranges.

swap_config_t gains the `nutrients` field. No reader, adapter,
or runtime wire-up yet — those land in subsequent commits.

pFUnit smoke + validator suite seeded with 6 tests.

Part of [nutrients] N2a (spec
docs/superpowers/specs/2026-05-07-nutrients-N2a-soil-side-initial-state-design.md).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: TOML reader for `[nutrients]`

**Files:**
- Create: `src/io/toml/read_nutrients_toml.f90`
- Modify: `src/io/toml/load_swap_config.f90` (call new reader)
- Modify: `meson.build` (production source)
- Modify: `tests/unit/meson.build` (`pfunit_extra_sources`)
- Modify: `tests/unit/io/toml/test_nutrients_config.pf` (append parser tests)

- [ ] **Step 1: Append reader tests**

Append to `tests/unit/io/toml/test_nutrients_config.pf` (use whatever in-memory or temp-file load pattern is established by other reader tests — see `tests/unit/io/toml/test_read_irrigation_ssdi_toml.pf` for the closest precedent):

```fortran
@test
subroutine test_read_nutrients_full_block()
   use, intrinsic :: iso_fortran_env, only: real64
   use nutrients_config_mod, only: nutrients_config_t
   use read_nutrients_toml_mod, only: read_nutrients_toml
   use error_mod, only: error_collection_t
   use tomlf, only: toml_table, toml_load
   use funit
   implicit none
   type(toml_table), allocatable, target :: doc
   type(toml_table), pointer :: doc_ptr
   type(nutrients_config_t)  :: cfg
   type(error_collection_t)  :: errs
   character(len=:), allocatable :: toml_text

   toml_text =                                                       &
      '[nutrients]'                                  // new_line('a') // &
      'sorp_coef = 0.005'                            // new_line('a') // &
      ''                                             // new_line('a') // &
      '[nutrients.initial]'                          // new_line('a') // &
      'fom  = [0.5, 0.3, 0.2, 0.1, 0.5, 0.3, 0.2, 0.1]' // new_line('a') // &
      'bio  = 0.4'                                   // new_line('a') // &
      'hum  = 8.0'                                   // new_line('a') // &
      'cnh4 = 0.001'                                 // new_line('a') // &
      'cno3 = 0.005'

   call toml_load(doc, toml_text)
   doc_ptr => doc

   call read_nutrients_toml(doc_ptr, cfg, errs)
   call assertEqual(0, errs%count(), 'expected zero parse errors')
   @assertTrue(cfg%present, 'expected present=true')
   @assertEqual(0.005_real64, cfg%sorp_coef,        tolerance=1.0e-12_real64)
   @assertEqual(0.5_real64,   cfg%initial%fom(1),   tolerance=1.0e-12_real64)
   @assertEqual(0.1_real64,   cfg%initial%fom(8),   tolerance=1.0e-12_real64)
   @assertEqual(0.4_real64,   cfg%initial%bio,      tolerance=1.0e-12_real64)
   @assertEqual(8.0_real64,   cfg%initial%hum,      tolerance=1.0e-12_real64)
   @assertEqual(0.001_real64, cfg%initial%cnh4,     tolerance=1.0e-12_real64)
   @assertEqual(0.005_real64, cfg%initial%cno3,     tolerance=1.0e-12_real64)
end subroutine

@test
subroutine test_read_nutrients_missing_block()
   use nutrients_config_mod, only: nutrients_config_t
   use read_nutrients_toml_mod, only: read_nutrients_toml
   use error_mod, only: error_collection_t
   use tomlf, only: toml_table, toml_load
   use funit
   implicit none
   type(toml_table), allocatable, target :: doc
   type(toml_table), pointer :: doc_ptr
   type(nutrients_config_t)  :: cfg
   type(error_collection_t)  :: errs

   call toml_load(doc, '[general]' // new_line('a') // 'project = "x"')
   doc_ptr => doc

   call read_nutrients_toml(doc_ptr, cfg, errs)
   call assertEqual(0, errs%count(), 'no errors when [nutrients] absent')
   @assertFalse(cfg%present, 'expected present=false')
end subroutine

@test
subroutine test_read_nutrients_partial_fom()
   use, intrinsic :: iso_fortran_env, only: real64
   use nutrients_config_mod, only: nutrients_config_t
   use read_nutrients_toml_mod, only: read_nutrients_toml
   use error_mod, only: error_collection_t
   use tomlf, only: toml_table, toml_load
   use funit
   implicit none
   type(toml_table), allocatable, target :: doc
   type(toml_table), pointer :: doc_ptr
   type(nutrients_config_t)  :: cfg
   type(error_collection_t)  :: errs

   call toml_load(doc,                                                &
      '[nutrients]'                       // new_line('a') //         &
      '[nutrients.initial]'               // new_line('a') //         &
      'fom = [0.5, 0.3, 0.2]')
   doc_ptr => doc

   call read_nutrients_toml(doc_ptr, cfg, errs)
   call assertEqual(0, errs%count())
   @assertEqual(0.5_real64, cfg%initial%fom(1), tolerance=1.0e-12_real64)
   @assertEqual(0.3_real64, cfg%initial%fom(2), tolerance=1.0e-12_real64)
   @assertEqual(0.2_real64, cfg%initial%fom(3), tolerance=1.0e-12_real64)
   @assertEqual(0.0_real64, cfg%initial%fom(4), tolerance=1.0e-12_real64)
   @assertEqual(0.0_real64, cfg%initial%fom(8), tolerance=1.0e-12_real64)
end subroutine

@test
subroutine test_read_nutrients_oversized_fom()
   use nutrients_config_mod, only: nutrients_config_t
   use read_nutrients_toml_mod, only: read_nutrients_toml
   use error_mod, only: error_collection_t
   use tomlf, only: toml_table, toml_load
   use funit
   implicit none
   type(toml_table), allocatable, target :: doc
   type(toml_table), pointer :: doc_ptr
   type(nutrients_config_t)  :: cfg
   type(error_collection_t)  :: errs

   call toml_load(doc,                                                &
      '[nutrients]'                                              // new_line('a') // &
      '[nutrients.initial]'                                      // new_line('a') // &
      'fom = [1, 2, 3, 4, 5, 6, 7, 8, 9]')
   doc_ptr => doc

   call read_nutrients_toml(doc_ptr, cfg, errs)
   @assertTrue(errs%count() > 0, 'expected parse error for fom size > 8')
end subroutine
```

(If `tests/unit/io/toml/test_read_irrigation_ssdi_toml.pf` uses a fixture-file pattern instead of in-memory `toml_load`, switch to that pattern. Inspect that file first.)

- [ ] **Step 2: Wire and run to verify failure**

Run:
```
pixi run -e test build-linux
```
Expected: build fails — `read_nutrients_toml_mod` doesn't exist.

- [ ] **Step 3: Create `src/io/toml/read_nutrients_toml.f90`**

```fortran
!> @file read_nutrients_toml.f90
!! Parses the optional top-level [nutrients] block into
!! nutrients_config_t. Block is optional; absent → leaves
!! defaults intact and sets present = .false..
!!
!! See ADR 0026 ([nutrients] N2a).
module read_nutrients_toml_mod

   use, intrinsic :: iso_fortran_env, only: real64
   use tomlf, only: toml_table, toml_array, get_value, len
   use toml_field_helpers_mod, only: get_optional_real_with_default
   use nutrients_config_mod, only: nutrients_config_t
   use error_mod, only: error_collection_t, ERR_PARSE_TYPE_MISMATCH
   implicit none
   private

   public :: read_nutrients_toml

contains

   subroutine read_nutrients_toml(doc, cfg, errors)
      type(toml_table), pointer, intent(in)    :: doc
      type(nutrients_config_t),  intent(inout) :: cfg
      type(error_collection_t),  intent(inout) :: errors

      type(toml_table), pointer :: nut_tbl, init_tbl
      type(toml_array), pointer :: fom_arr
      integer :: stat, n, i
      real(real64) :: v

      if (.not. associated(doc)) return

      nut_tbl => null()
      call get_value(doc, 'nutrients', nut_tbl, requested=.false., stat=stat)
      if (stat /= 0 .or. .not. associated(nut_tbl)) return

      cfg%present = .true.
      call get_optional_real_with_default(nut_tbl, 'sorp_coef', cfg%sorp_coef, &
                                          0.0_real64, 'nutrients.sorp_coef', errors)

      ! [nutrients.initial] sub-table is optional.
      init_tbl => null()
      call get_value(nut_tbl, 'initial', init_tbl, requested=.false., stat=stat)
      if (stat == 0 .and. associated(init_tbl)) then
         ! fom is a 1- to 8-element array; missing entries stay at 0.0.
         fom_arr => null()
         call get_value(init_tbl, 'fom', fom_arr, requested=.false., stat=stat)
         if (stat == 0 .and. associated(fom_arr)) then
            n = len(fom_arr)
            if (n > 8) then
               call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                  'nutrients.initial.fom: at most 8 entries (legacy FOM_t(1..8) cap)', &
                  'nutrients.initial.fom')
               n = 8
            end if
            do i = 1, n
               call get_value(fom_arr, i, v, stat=stat)
               if (stat == 0) cfg%initial%fom(i) = v
            end do
         end if

         call get_optional_real_with_default(init_tbl, 'bio',  cfg%initial%bio,  &
                                             0.0_real64, 'nutrients.initial.bio',  errors)
         call get_optional_real_with_default(init_tbl, 'hum',  cfg%initial%hum,  &
                                             0.0_real64, 'nutrients.initial.hum',  errors)
         call get_optional_real_with_default(init_tbl, 'cnh4', cfg%initial%cnh4, &
                                             0.0_real64, 'nutrients.initial.cnh4', errors)
         call get_optional_real_with_default(init_tbl, 'cno3', cfg%initial%cno3, &
                                             0.0_real64, 'nutrients.initial.cno3', errors)
      end if
   end subroutine read_nutrients_toml

end module read_nutrients_toml_mod
```

- [ ] **Step 4: Wire into `load_swap_config.f90`**

Add a `use` line near the top of the module:
```fortran
   use read_nutrients_toml_mod, only: read_nutrients_toml
```

Find the chain of `read_*_toml` calls (around lines 45-56). After the `read_output_csv_toml` call (the last one), append:
```fortran
      call read_nutrients_toml(doc_ptr, config%nutrients, errors)
```

- [ ] **Step 5: Register sources in meson**

In `meson.build`, add `'src/io/toml/read_nutrients_toml.f90'` near the other `src/io/toml/` files in the `sources` list.

In `tests/unit/meson.build`, add to `pfunit_extra_sources`:
```meson
        '../../src/io/toml/read_nutrients_toml.f90',
```

- [ ] **Step 6: Build and run tests**

```
pixi run -e test build-linux
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:|tests,"
pixi run -e test check-full 2>&1 | grep "Results:"
```
Expected: build clean; pFUnit `Ok: 1, Fail: 0` (test count grows by 4 — the 4 parser tests added in Step 1); check-full `5 passed, 0 failed`.

- [ ] **Step 7: Commit**

```bash
git add src/io/toml/read_nutrients_toml.f90 src/io/toml/load_swap_config.f90 meson.build tests/unit/meson.build tests/unit/io/toml/test_nutrients_config.pf
git commit -m "$(cat <<'EOF'
feat(io): parse [nutrients] into nutrients_config_t

Adds read_nutrients_toml_mod, called from load_swap_config after
the existing read_*_toml chain. Block is optional — when
[nutrients] is absent, reader leaves defaults intact and sets
present = .false..

[nutrients.initial] sub-table is also optional. fom array
accepts 1 to 8 elements; oversized arrays append a parse error.

pFUnit suite test_nutrients_config gains 4 parser tests
(full-block, missing-block, partial-fom, oversized-fom).

Part of [nutrients] N2a (spec
docs/superpowers/specs/2026-05-07-nutrients-N2a-soil-side-initial-state-design.md).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: `apply_nutrients` adapter + adapter test

**Files:**
- Modify: `src/io/toml/config_to_variables.f90`
- Modify: `tests/unit/io/toml/test_nutrients_config.pf` (append adapter tests)

- [ ] **Step 1: Append adapter tests**

Append to `tests/unit/io/toml/test_nutrients_config.pf`:

```fortran
@test
subroutine test_apply_nutrients_populates_globals()
   use, intrinsic :: iso_fortran_env, only: real64
   use nutrients_config_mod, only: nutrients_config_t
   use config_to_variables_mod, only: apply_nutrients
   use wofost_soil_declarations, only: FOM_t, Bio_t, Hum_t, &
                                        cNH4_t, cNO3_t, SorpCoef
   use funit
   implicit none
   type(nutrients_config_t) :: cfg
   integer :: i

   cfg%present   = .true.
   cfg%sorp_coef = 0.005_real64
   cfg%initial%fom = [0.5_real64, 0.3_real64, 0.2_real64, 0.1_real64, &
                      0.5_real64, 0.3_real64, 0.2_real64, 0.1_real64]
   cfg%initial%bio  = 0.4_real64
   cfg%initial%hum  = 8.0_real64
   cfg%initial%cnh4 = 0.001_real64
   cfg%initial%cno3 = 0.005_real64

   call apply_nutrients(cfg)

   @assertEqual(0.005_real64, SorpCoef, tolerance=1.0e-12_real64)
   do i = 1, 8
      @assertEqual(cfg%initial%fom(i), FOM_t(i), tolerance=1.0e-12_real64)
   end do
   @assertEqual(0.4_real64,   Bio_t,  tolerance=1.0e-12_real64)
   @assertEqual(8.0_real64,   Hum_t,  tolerance=1.0e-12_real64)
   @assertEqual(0.001_real64, cNH4_t, tolerance=1.0e-12_real64)
   @assertEqual(0.005_real64, cNO3_t, tolerance=1.0e-12_real64)
end subroutine

@test
subroutine test_apply_nutrients_defaults_zero()
   use, intrinsic :: iso_fortran_env, only: real64
   use nutrients_config_mod, only: nutrients_config_t
   use config_to_variables_mod, only: apply_nutrients
   use wofost_soil_declarations, only: FOM_t, Bio_t, Hum_t, &
                                        cNH4_t, cNO3_t, SorpCoef
   use funit
   implicit none
   type(nutrients_config_t) :: cfg
   integer :: i

   ! Non-zero pre-existing values to detect overwrite-with-zero.
   FOM_t(1) = 99.0_real64
   Bio_t    = 99.0_real64
   Hum_t    = 99.0_real64
   cNH4_t   = 99.0_real64
   cNO3_t   = 99.0_real64
   SorpCoef = 99.0_real64

   ! cfg%present = .false. (default); apply_nutrients still runs and resets.
   call apply_nutrients(cfg)

   @assertEqual(0.0_real64, SorpCoef, tolerance=1.0e-12_real64)
   do i = 1, 8
      @assertEqual(0.0_real64, FOM_t(i), tolerance=1.0e-12_real64)
   end do
   @assertEqual(0.0_real64, Bio_t,  tolerance=1.0e-12_real64)
   @assertEqual(0.0_real64, Hum_t,  tolerance=1.0e-12_real64)
   @assertEqual(0.0_real64, cNH4_t, tolerance=1.0e-12_real64)
   @assertEqual(0.0_real64, cNO3_t, tolerance=1.0e-12_real64)
end subroutine
```

- [ ] **Step 2: Run to verify failure**

```
pixi run -e test build-linux
```
Expected: build fails — `apply_nutrients` is not exported by `config_to_variables_mod`.

- [ ] **Step 3: Add the adapter to `config_to_variables.f90`**

Add to the module's `public ::` section:
```fortran
   public :: apply_nutrients
```

Add the new procedure inside the `contains` block, after the existing `apply_irrigation_ssdi` adapter:

```fortran
   !> Apply [nutrients] config to legacy `variables`/wofost_soil_declarations
   !! globals. Always called from config_to_variables (no flCropNut gate);
   !! the cfg%present flag is informational only — defaults are zero
   !! whether or not the user supplied a [nutrients] block.
   !!
   !! Sets SorpCoef unconditionally — fixes the genuine uninitialised-
   !! variable bug discovered during the [nutrients] N2 brainstorm.
   !!
   !! See ADR 0026 ([nutrients] N2a).
   subroutine apply_nutrients(cfg)
      use nutrients_config_mod, only: nutrients_config_t
      use wofost_soil_declarations, only: FOM_t, Bio_t, Hum_t, &
                                           cNH4_t, cNO3_t, SorpCoef
      type(nutrients_config_t), intent(in) :: cfg
      integer :: i

      SorpCoef = cfg%sorp_coef
      do i = 1, 8
         FOM_t(i) = cfg%initial%fom(i)
      end do
      Bio_t  = cfg%initial%bio
      Hum_t  = cfg%initial%hum
      cNH4_t = cfg%initial%cnh4
      cNO3_t = cfg%initial%cno3
   end subroutine apply_nutrients
```

Wire the call. Find the existing line `if (flSSDI) call apply_irrigation_ssdi(config%irrigation%ssdi)` (around line 489). Right after it, add:
```fortran
      call apply_nutrients(config%nutrients)
```

(No `if`-guard — the adapter is always called, defaults handle the absent-block case.)

- [ ] **Step 4: Build and run tests**

```
pixi run -e test build-linux
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:|tests,"
pixi run -e test check-full 2>&1 | grep "Results:"
```
Expected: build clean; pFUnit `Ok: 1, Fail: 0` (test count grows by 2); check-full `5 passed, 0 failed` (regression cases all leave `[nutrients]` absent → adapter sets pools/SorpCoef to zero, which is the same as their pre-arc uninitialised behaviour for any case that doesn't read them).

- [ ] **Step 5: Commit**

```bash
git add src/io/toml/config_to_variables.f90 tests/unit/io/toml/test_nutrients_config.pf
git commit -m "$(cat <<'EOF'
feat(io): apply_nutrients populates legacy globals from typed config

Adds the apply_nutrients helper to config_to_variables_mod,
called unconditionally (not gated on any flag) from the top-level
adapter. Copies cfg%sorp_coef + initial pool concentrations
(FOM_t(1..8), Bio_t, Hum_t, cNH4_t, cNO3_t) to the legacy
variables-module globals declared in wofost_soil_declarations.

The unconditional SorpCoef assignment fixes a genuine
uninitialised-variable bug discovered during the N2 brainstorm:
SorpCoef was used by wofost_soil_watern, wofost_soil_balancecheck,
wofost_soil_amendments, and wofost_soil_cropresidues but never
assigned anywhere in the modern src/. Default 0.0 = no sorption.

cfg%present = .false. (the default when no [nutrients] block is
supplied) sets all globals to 0.0 — same as pre-arc behaviour
for any case that doesn't reach a code path reading these
globals. check-full byte-identical.

pFUnit suite test_nutrients_config gains 2 adapter tests
(populates-globals, defaults-to-zero).

Part of [nutrients] N2a (spec
docs/superpowers/specs/2026-05-07-nutrients-N2a-soil-side-initial-state-design.md).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Retire `SoilManagement(7)` body

**Files:**
- Modify: `src/crop/management_soil.f90` (collapse case-7 body; drop unreferenced locals)

- [ ] **Step 1: Read the current case-7 boundaries**

```bash
sed -n '585,660p' src/crop/management_soil.f90
```

You'll see the `case (7)` body starting around line 587 and ending just before the `case default` (around line ~657). The body has two sections:
- A `flCropExt` block (lines 589-608) that reads/rewrites the `cropext` external file.
- A `<project>.snp` template-read + `<project>_nut.end` template-write loop (the bulk).

Both go.

- [ ] **Step 2: Replace the case-7 body**

Replace lines 587 through the line before `case default` with:

```fortran
      case (7)
      ! [nutrients] N2a (ADR 0026): the legacy diagnostic dump
      ! used <project>.snp as a template (read each line; if a
      ! known key match, write the current pool value to
      ! <project>_nut.end; else echo). Without .snp files in the
      ! modern flow, the read would fail at runtime if flCropNut
      ! were ever true.
      !
      ! The flCropExt write-back at the top of the legacy case (7)
      ! is also retired — that branch was driven by a separate
      ! mechanism (cropext) tied to the N-P-K case-1 init that's
      ! now also gone.
      !
      ! A future arc can add a clean CSV-style nutrient-pool dump
      ! if anyone asks.
      return
```

- [ ] **Step 3: Drop unreferenced locals**

After the case-7 body is gone, several locals declared at the top of `SoilManagement` may become unreferenced:

```bash
grep -nE "\\b(snp|oup|cropext|filnamce|cstring|stat|line)\\b" src/crop/management_soil.f90
```

For each name, check whether it appears in any other case body. Most likely the dead-decl list is:
- `snp` (case-7 only)
- `oup` (case-7 only)
- `cropext` (cases 5-7; **may still be used in case 5/6** — verify)
- `filnamce` (case-7 only?)
- `cstring(:)` (legacy global in `wofost_soil_declarations.f90`; declared at module level, NOT a local — leave alone unless you confirm it's also unreferenced everywhere)
- `line` (case-7 only?)

Drop only the locals that have **zero remaining references** in the file. Keep `flCropExt` (it's a module global, used by case 5/6 — verify by grep before assuming it's dead). When in doubt, leave the declaration.

If a local declaration is on a comma-separated list (`integer task, sme, smm, snp, ...`), drop only the names that became unused.

- [ ] **Step 4: Build and run**

```
pixi run -e test build-linux
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:|tests,"
pixi run -e test check-full 2>&1 | grep "Results:"
echo "=== acceptance ==="
grep -nE "snp\\b" src/crop/management_soil.f90 || echo "OK: no snp refs"
grep -nE "_nut\\.end" src/crop/management_soil.f90 || echo "OK: no _nut.end refs"
```
Expected: build clean (no warnings about unused vars); pFUnit unchanged from Task 3 (no new tests); check-full `5 passed, 0 failed`; both greps return "OK: ...".

- [ ] **Step 5: Commit**

```bash
git add src/crop/management_soil.f90
git commit -m "$(cat <<'EOF'
refactor(crop): retire SoilManagement(7) _nut.end diagnostic dump

The legacy case-7 body read <project>.snp line-by-line as a
template and wrote a partially-substituted copy to
<project>_nut.end. Without .snp files in the modern flow, the
read would fail at runtime if flCropNut were ever true.

The flCropExt write-back at the top of the legacy case (7) is
also retired — that branch was driven by a separate mechanism
(cropext) tied to the N-P-K case-1 init that's now gone.

Case body collapses to `return` + comment. Future arc can add a
clean CSV-style nutrient-pool dump if anyone asks. Unused locals
that referenced only the case-7 body (snp, oup, filnamce, etc.)
were dropped.

Verified: build clean (no unused-variable warnings); pFUnit Ok: 1,
Fail: 0; check-full 5/5 (cases 2-6 of SoilManagement still gated
on flCropNut=false, so the runtime path is unaffected).

Part of [nutrients] N2a (spec
docs/superpowers/specs/2026-05-07-nutrients-N2a-soil-side-initial-state-design.md).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: ADR 0026 + index update + schema doc

**Files:**
- Create: `docs/adr/0026-nutrients-N2a-soil-side-initial-state.md`
- Modify: `docs/adr/index.md` (append ADR 0026 row)
- Modify: `docs/configuration-schema.md` (document `[nutrients]` block)

- [ ] **Step 1: Write `docs/adr/0026-nutrients-N2a-soil-side-initial-state.md`**

```markdown
---
title: "ADR 0026 — [nutrients] N2a: soil-side initial state + SorpCoef"
date: 2026-05-07
status: accepted
---

# ADR 0026: [nutrients] N2a — soil-side initial state + SorpCoef

## Context

The `[nutrients]` umbrella reactivates the SWAP nutrient
subsystem on the modern TOML pipeline. ADR 0025 (N1) wired the
crop-side parameter set per-rotation. This ADR (N2a) ports the
soil-side initial state.

The legacy code read soil-nutrient initial concentrations from
`<project>.snp` (12 reals: `FOM_t(1..8)`, `Bio_t`, `Hum_t`,
`cNH4_t`, `cNO3_t`) inside `SoilManagement(1)`, which was
collapsed to `return` in SS-C step 3 of the legacy-readers
deletion arc. The `<project>_nut.end` diagnostic dump in
`SoilManagement(7)` then read `<project>.snp` line-by-line as a
template, writing the current pool values to `_nut.end`.

Pre-flight inventory for N2 found three legacy files:
- `<project>.snp` — initial pool concentrations (this ADR's
  scope).
- `<project>.smm` — material-property overrides (deferred
  indefinitely; hardcoded defaults in `Wofost_SoilParameters`).
- `<project>.sme` — timed soil management events / fertilizer
  applications (N2b's scope).

Plus three smaller findings:
- `SorpCoef` is genuinely uninitialised in modern `src/` — used
  by `wofost_soil_watern`, `wofost_soil_balancecheck`,
  `wofost_soil_amendments`, `wofost_soil_cropresidues` but never
  assigned anywhere. Real bug.
- `DryBD` is auto-derived from `BDENS(1)` (the soil's bulk
  density) in `Wofost_SoilParameters`. No TOML port needed.
- `<project>_nut.end` writer would fail at runtime without a
  `.snp` template file present.

## Decision

This ADR (N2a):

1. **Add a top-level `[nutrients]` TOML block.** Optional;
   `present` flag distinguishes "user supplied" from "absent —
   use defaults". Carries `sorp_coef` (top-level scalar) and an
   `[nutrients.initial]` sub-table with the 12 pool values.
2. **Default `SorpCoef = 0.0`.** This is a behavioural choice
   — legacy code left `SorpCoef` uninitialised, so the
   deterministic 0.0 default may differ numerically from any
   legacy non-deterministic baseline. 0.0 = "no sorption", which
   is the physically conservative choice and matches what most
   compilers' module-init zero would produce.
3. **`apply_nutrients` is called unconditionally** from
   `config_to_variables`, regardless of `flCropNut`. The adapter
   sets `SorpCoef` and pool values whether or not nutrients are
   enabled at runtime. When `cfg%present = .false.`, the
   defaults (zero) flow through — same as the pre-arc
   uninitialised behaviour for any case that didn't read these
   globals.
4. **Retire `SoilManagement(7)`'s `_nut.end` dump.** Body
   collapses to `return`. The `flCropExt` write-back at the top
   of the legacy case is also dropped (it depended on the same
   `cropext` mechanism). A future arc can add a clean CSV-style
   nutrient-pool dump if anyone asks.

## What N2a does NOT do

- **N2b — timed soil management events.** Without N2b, simulations
  run with **zero amendments**: natural mineralization only.
  That's a useful baseline for testing.
- **Material-property overrides** (legacy `<project>.smm`). The
  17 hardcoded materials in `Wofost_SoilParameters` (Cattle
  manure → Spruce needles, with default AppAge / OrgMatFrac /
  OrgNFrac / NH4NFrac / NO3NFrac) stay as-is. Deferred
  indefinitely.
- **Lift `tillage.f90:73`.** N3's job, after both N2a + N2b
  land. With no amendments wired (N2b pending), a
  flCropNut=true simulation would run from zero pools with no
  fertilizer inputs — useful but not the full picture.
- **Add a regression case** with `flCropNut=true`. Cannot run
  end-to-end until the runtime gate is lifted (N3).
- **Replace `_nut.end` with a CSV writer.** Future enhancement.

## Schema

```toml
[nutrients]
sorp_coef = 0.005

[nutrients.initial]
fom  = [0.5, 0.3, 0.2, 0.1, 0.5, 0.3, 0.2, 0.1]   # FOM_t(1..8)
bio  = 0.4
hum  = 8.0
cnh4 = 0.001
cno3 = 0.005
```

Both `[nutrients]` and `[nutrients.initial]` are optional;
absent → defaults (zero pools, zero sorption).

## Consequences

- A TOML config with `[nutrients]` populated loads and
  validates without error; pool values flow through to the
  legacy globals.
- `SorpCoef` now has a deterministic value (0.0 by default,
  user-overridable). Eliminates the uninitialised-variable bug.
- `<project>.snp` and `<project>_nut.end` files no longer have
  any role in the modern pipeline. The diagnostic dump path is
  retired.
- All five existing regression cases produce byte-identical
  CSV outputs (no case has flCropNut=true; no case reads the
  pool globals or SorpCoef on its active code path).
- N3 will lift `tillage.f90:73` and add a regression case with
  populated `[nutrients]`; verifying byte-identical against the
  legacy binary in N3 is the end-to-end correctness gate.

## Acceptance

- `pixi run -e test test-pfunit` → `Ok: 1, Fail: 0`; suite
  count grows by 12 (the eleven N2a tests + one smoke; final
  count depends on Task-1 vs Task-2 vs Task-3 split).
- `pixi run -e test check-full` → `5 passed, 0 failed`
  (byte-identical CSVs).
- `grep -n "_nut\.end" src/crop/management_soil.f90` → no
  matches.
- `grep -n "snp\b" src/crop/management_soil.f90` → no matches.
- `grep -n "flCropNut.*not.*allowed" src/crop/tillage.f90` →
  one match (runtime gate intact).
- ADR 0026 committed.

## Related

- ADR 0008 — error collection over fatalerr.
- ADR 0024 — `dtutil.f90` compatibility shim — same
  architectural thread.
- ADR 0025 — [nutrients] N1 (crop-side adapter).
- Future: ADR 0027 (N2b: amendments).
- Future: ADR 0028 (N3: lift runtime gate + regression case).
```

- [ ] **Step 2: Append the index row**

In `docs/adr/index.md`, after the ADR 0025 row, add:

```markdown
- [ADR 0026 — [nutrients] N2a: soil-side initial state + SorpCoef](0026-nutrients-N2a-soil-side-initial-state.html) — Top-level `[nutrients]` TOML block carries initial pool concentrations (12 reals) and SorpCoef. Adapter populates legacy globals unconditionally; SorpCoef gets a deterministic default fixing a genuine uninitialised-variable bug. SoilManagement(7) `_nut.end` dump retired (broken without .snp template files).
```

- [ ] **Step 3: Document `[nutrients]` in `configuration-schema.md`**

Open `docs/configuration-schema.md`. Find a natural location for a new top-level section (the file likely has `[soil]`, `[crop]`, etc. — append after one of those, or in a new "Nutrients" section near the bottom). Mirror the formatting of the surrounding sections.

Cover:
- **Top-level `[nutrients]`** activation (optional; absent → defaults zero).
- `sorp_coef` (real, ≥ 0, default 0.0; m³/kg).
- **`[nutrients.initial]` sub-table** (optional within `[nutrients]`):
  - `fom` (real array, 1..8 elements, each in [0, 1000], default zeros) — FOM_t(1..8) initial pool concentrations (kg/m²).
  - `bio` (real, [0, 1000], default 0.0) — Bio_t initial.
  - `hum` (real, [0, 1000], default 0.0) — Hum_t initial.
  - `cnh4` (real, [0, 1000], default 0.0) — cNH4_t initial concentration (kg/m³).
  - `cno3` (real, [0, 1000], default 0.0) — cNO3_t initial concentration.

Match the existing sub-section formatting (heading depth, table style, prose density).

- [ ] **Step 4: Final verification**

```
pixi run -e test build-linux
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:|tests,"
pixi run -e test check-full 2>&1 | grep "Results:"
echo "=== acceptance ==="
grep -nE "_nut\\.end" src/crop/management_soil.f90 || echo "OK: no _nut.end"
grep -nE "snp\\b" src/crop/management_soil.f90 || echo "OK: no snp refs"
grep -n "flCropNut.*not.*allowed" src/crop/tillage.f90
ls src/config/nutrients_config.f90 src/io/toml/read_nutrients_toml.f90
```
Expected: build clean; pFUnit Ok: 1, Fail: 0 (~12 new tests); check-full 5/5; both greps return "OK: ..."; the `flCropNut.*not.*allowed` grep returns one match (runtime gate intact); both new files exist.

- [ ] **Step 5: Commit**

```bash
git add docs/adr/0026-nutrients-N2a-soil-side-initial-state.md docs/adr/index.md docs/configuration-schema.md
git commit -m "$(cat <<'EOF'
docs: ADR 0026 — [nutrients] N2a soil-side initial state + SorpCoef

Captures the second sub-arc of the [nutrients] umbrella. Adds a
top-level [nutrients] TOML block carrying initial pool
concentrations (FOM_t(1..8), Bio_t, Hum_t, cNH4_t, cNO3_t) and
SorpCoef. Adapter populates legacy globals unconditionally —
fixes the genuine uninitialised-variable bug for SorpCoef.
SoilManagement(7) _nut.end dump retired.

Runtime gate at tillage.f90:73 stays — N3's job once both N2a
and N2b are in place.

- docs/adr/0026-nutrients-N2a-soil-side-initial-state.md: new
  ADR with context, decision, what-N2a-doesn't-do (N2b/N3/.smm
  out of scope), consequences, acceptance criteria, forward
  references to ADR 0027 (N2b) and ADR 0028 (N3).
- docs/adr/index.md: ADR 0026 row.
- docs/configuration-schema.md: documents [nutrients] +
  [nutrients.initial] fields, ranges, optionality.

Closes the [nutrients] N2a spec
(docs/superpowers/specs/2026-05-07-nutrients-N2a-soil-side-initial-state-design.md).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Self-Review

**Spec coverage:**

| Spec section | Task |
|---|---|
| §Schema (`nutrients_initial_t`, `nutrients_config_t`) | Task 1 |
| §Validation rules | Task 1 |
| §TOML reader (`read_nutrients_toml`) | Task 2 |
| §Adapter (`apply_nutrients`, unconditional call) | Task 3 |
| §Runtime change (`SoilManagement(7)` retired) | Task 4 |
| §Tests (11 pFUnit suites) | Tasks 1, 2, 3 (6 + 4 + 2 = 12; one extra is the default-init smoke test) |
| §ADR 0026 + schema doc | Task 5 |

All spec sections covered.

**Acceptance criteria** (from spec):

- `pixi run -e test build-linux` clean — verified each task.
- `pixi run -e test test-pfunit` — verified each task; final +12.
- `pixi run -e test check-full` 5/5 byte-identical — verified Tasks 1, 3, 4, 5.
- No `snp\b` / `_nut.end` in management_soil.f90 — verified Task 4.
- `tillage.f90:73` UNCHANGED — verified Tasks 4, 5 acceptance grep.
- `apply_nutrients` called unconditionally — Task 3.
- ADR 0026 committed — Task 5.

**Type / signature consistency:**

- `nutrients_initial_t` fields (`fom(8)`, `bio`, `hum`, `cnh4`, `cno3`) — used identically in Tasks 1, 2, 3.
- `nutrients_config_t` fields (`present`, `sorp_coef`, `initial`) — consistent.
- `apply_nutrients(cfg)` signature — consistent.
- `read_nutrients_toml(doc, cfg, errors)` signature — consistent.
- Module variable names in `wofost_soil_declarations` (`FOM_t`, `Bio_t`, `Hum_t`, `cNH4_t`, `cNO3_t`, `SorpCoef`) — consistent with what the adapter writes and the test asserts.

**Open questions punted to implementation:**

- Whether `tests/unit/io/toml/test_read_irrigation_ssdi_toml.pf` uses an in-memory `toml_load` pattern or a fixture-file pattern. Task 2 step 1 directs the implementer to inspect that file before writing the parser tests; either pattern is fine if applied consistently.
- The exact list of unreferenced locals after `case (7)` collapse (Task 4 Step 3). The implementer greps each candidate before dropping; conservative when in doubt.
- Whether `flCropExt` is reachable in the modern flow. Spec mentions it's "likely unreachable"; the case-5 / case-6 references should be checked. Task 4 Step 3 covers this — drop only if confirmed dead.

These are pickable at execution time without re-planning.

---

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-05-07-nutrients-N2a-soil-side-initial-state.md`.
