# Phase 4c-a — Cross-File TOML + Fixed and Grass Crops Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Bring the new TOML pipeline to full field-by-field parity for cases 2.grassgrowth, 4.oxygenstress, 6.surfacewater (full type-1 + type-3 coverage), partial parity for 1.hupselbrook (its type-1 + type-3 entries; type-2 deferred), and eliminate the three Phase 4b workarounds (date convention, swmacro shadow, nrlevs hardcode).

**Architecture:** Same composite thematic TOML reader pattern from Phase 4a/4b, extended with single-entry recursive cross-file loading (`[drainage].file = "..."` and `[[crop.rotation]].file = "..."`). Two new crop config types (`cropfixed_config_t`, `cropgrass_config_t`) and matching section readers. Three small targeted source-tree fixes bundled at the start.

**Tech Stack:** gfortran 2008, pFUnit 4.15, meson + pixi, toml-f, `iso_c_binding` for `chdir` (already added Phase 4b).

Spec: `docs/superpowers/specs/2026-04-25-phase-4c-a-crop-fixed-grass-design.md`

---

## Preamble: context every task needs

**Baseline:** Phase 4b complete at commit `9fb83b6`, tag `rescue/phase-4b-parity`. 109 pFUnit tests passing, 6/6 regression green. `swap.f90` runs the pre-drift legacy path (no state_t threading). `crop_config_t` holds rotation metadata only. Cross-file TOML loading is NOT yet supported — all sections inline in `swap.toml`.

**Branch discipline:**
- Work on `development`. No per-task feature branches in Phase 4c-a.
- One commit per task. Subject `<type>(<scope>): <what>`.
- No pushes to origin during Phase 4c-a (matches the rescue policy).
- Phase exit: fast-forward `main` to `development` locally; tag `rescue/phase-4c-a-crop-fixed-grass`.

**Working directory:** `/home/zawadzkim/Code/swap` for all tasks unless noted.

**Conventions** (unchanged from 4a/4b):
- Module name = filename stem + `_mod`.
- `implicit none` after module statement; default `private`; explicit `public ::`.
- `use iso_fortran_env, only: real64` for new code; legacy idioms can stay.
- pFUnit suite filename = `tests/unit/<domain>/test_<module>.pf`; `ADD_TEST_SUITE(test_<module>_suite)` in `testSuites.inc`; `'<domain>/test_<module>.pf'` in `pf_files`.
- **No inline `! comments` on `@assert*` lines** — `funitproc` chokes on the `!`.
- Module sources go in `pfunit_extra_sources` if outside `test_base_sources`.

**Verification after each task:**
```
pixi run -e test test-pfunit         # pFUnit must stay green; count notes
pixi run -e test check-fast          # 4/6 regression cases; <90s
```
After source-tree fixes (Stage A), also run:
```
pixi run -e test check-full          # 6/6 regression cases; ~10 min
```

**One concept worth understanding before authoring tasks:**

The Phase 4b parity test calls legacy `readswap()` from inside a pFUnit test. `readswap` has the `Get_Command_Argument` fallback to `swap.swp`, so we chdir into the case dir and stage `swap_linux.swp.template` → `swap.swp` first. The crop readers (`readcropfixed`, `readcropgrass`) live INSIDE `readswap` — they execute in sequence per rotation entry, and each one OVERWRITES the same set of `variables` globals. This means after `readswap` returns, only the LAST rotation entry's crop globals are still in `variables`. Phase 4c-a parity asserts only on the LAST entry's loaded crop config vs `variables`.

---

## File structure

### New source files (all under existing directories)

| File | Responsibility |
|---|---|
| `src/io/toml/path_helpers.f90` | `directory_of(path)`, `resolve_relative_path(base, rel)`. Pure string helpers. |
| `src/config/cropfixed_config.f90` | `cropfixed_config_t` (~50-80 fields) + `validate` + `finalize`. |
| `src/config/cropgrass_config.f90` | `cropgrass_config_t` (~80-110 fields) + `validate` + `finalize`. |
| `src/io/toml/read_cropfixed_toml.f90` | `read_cropfixed_toml(doc_root, config, errors)` — populates a `cropfixed_config_t` from a `.crp.toml` root. |
| `src/io/toml/read_cropgrass_toml.f90` | Same shape, for grass. |

### Files modified

| File | Change |
|---|---|
| `src/io/toml/toml_field_helpers.f90` | `parse_date_to_days1900` JDN constant: `2415021` → `2415020`. |
| `src/io/readswap.f90` | Delete the local `integer SwMacro` declaration (line ~20). |
| `src/config/drainage_config.f90` | Extend `drainage_config_finalize` with the `nrlevs` mirror. |
| `src/config/crop_config.f90` | Add `rotation_fixed(:)` and `rotation_grass(:)` allocatable arrays; extend `validate` to delegate per entry. |
| `src/io/toml/read_drainage_toml.f90` | Add optional `base_path` arg; if `[drainage].file` present, load that file and read `[drainage]` from it. Extract `read_drainage_inner` for shared body. |
| `src/io/toml/read_crop_toml.f90` | Add optional `base_path` arg; per rotation entry, if `file` present, load that file and dispatch to `read_cropfixed_toml` / `read_cropgrass_toml` per `type`. |
| `src/io/toml/load_swap_config.f90` | Compute `base_dir = directory_of(path)`; pass `base_path=base_dir` to drainage and crop readers. |
| `tests/unit/io/toml/test_hupselbrook_parity.pf` | Re-enable `swmacro` and `nrlevs` assertions; tighten date tolerance; add type-1 (maizes) and type-3 (grassd) sub-section assertions. |
| `tests/unit/io/toml/toml_field_helpers test`, `test_hupselbrook_roundtrip` | Update expected date constants (37255 → 37257). |
| `meson.build`, `tests/unit/meson.build`, `tests/unit/testSuites.inc` | Wire new sources and suites. |
| `docs/configuration-schema.md` | Reconcile to Phase 4b/4c-a reality. |

### New test files

| File | Tests |
|---|---|
| `tests/unit/io/toml/test_path_helpers.pf` | Direct unit tests for `directory_of` and `resolve_relative_path`. |
| `tests/unit/config/test_cropfixed_config.pf` | Validator tests for `cropfixed_config_t` (~5-8 tests). |
| `tests/unit/config/test_cropgrass_config.pf` | Validator tests for `cropgrass_config_t` (~5-8 tests). |
| `tests/unit/io/toml/test_read_cropfixed_toml.pf` | Happy + missing-section + malformed (3 tests). |
| `tests/unit/io/toml/test_read_cropgrass_toml.pf` | Same shape (3 tests). |
| `tests/unit/io/toml/test_grassgrowth_parity.pf` | Per-case parity (~6-8 tests). |
| `tests/unit/io/toml/test_oxygenstress_parity.pf` | Per-case parity. |
| `tests/unit/io/toml/test_surfacewater_parity.pf` | Per-case parity. |

### New TOML files (in the `tests/swap-cases` submodule)

Per case, under `tests/swap-cases/toml/<case>/`:
- `1.hupselbrook/`: extend `swap.toml`; new `swap.dra.toml`, `maizes.crp.toml`, `grassd.crp.toml`; placeholder `potatod.crp.toml` (Phase 4c-b populates).
- `2.grassgrowth/`: extend `swap.toml`; new `swap.dra.toml`, `grassd.crp.toml`.
- `4.oxygenstress/`: extend `swap.toml`; new `swap.dra.toml`, `grassd.crp.toml`.
- `6.surfacewater/`: extend `swap.toml`; new `swap.dra.toml`, `grass.crp.toml`.

---

## Stage A — Targeted source fixes (Tasks 1-3)

These three changes unblock the parity tests by removing Phase 4b workarounds. Each is small, verified by `check-full`.

### Task 1: Align `parse_date_to_days1900` to legacy 1-based JDN

**Files:**
- Modify: `src/io/toml/toml_field_helpers.f90`
- Modify: `tests/unit/io/toml/test_toml_field_helpers.pf`
- Modify: `tests/unit/io/toml/test_hupselbrook_roundtrip.pf` (if it has hardcoded date constants)

- [ ] **Step 1: Update the JDN constant**

In `src/io/toml/toml_field_helpers.f90`, find `parse_date_to_days1900`. The line:
```fortran
jd1900 = 2415021
```
becomes:
```fortran
jd1900 = 2415020
```

After this change, `2002-01-01` → `37257` (was `37255`); `2004-12-31` → `38352` (was `38350`). These match what the legacy ttutil reader produces.

- [ ] **Step 2: Update the field-helpers unit test**

In `tests/unit/io/toml/test_toml_field_helpers.pf`, find `test_parse_date_to_days1900`. The expected value:
```fortran
@assertEqual(37255.0_real64, days, 1.0_real64)
```
becomes:
```fortran
@assertEqual(37257.0_real64, days, 1.0_real64)
```

- [ ] **Step 3: Check the round-trip test for date constants**

```
grep -n "37255\|37256\|38350\|38351" tests/unit/io/toml/test_hupselbrook_roundtrip.pf
```

For each occurrence, increment by 2 (37255→37257, 38350→38352). The round-trip itself is symmetric (load→emit→load), so the comparison is between two `c1` and `c2` values that move together; the date constant changes only matter where a HARDCODED expected value appears. If grep returns no hits, no edit needed.

- [ ] **Step 4: Build + test**

```
pixi run -e test test-pfunit
```

Expected: 109 tests still pass. The `test_parse_date_to_days1900` test now asserts 37257 against actual 37257.

```
pixi run -e test check-full
```

Expected: 6/6 regression cases green.

- [ ] **Step 5: Commit**

```
git add src/io/toml/toml_field_helpers.f90 \
        tests/unit/io/toml/test_toml_field_helpers.pf
# also if changes to roundtrip:
# git add tests/unit/io/toml/test_hupselbrook_roundtrip.pf
git commit -m "fix(toml): align parse_date_to_days1900 to legacy 1-based JDN"
```

---

### Task 2: Remove `swmacro` shadow from `readswap.f90`

**Files:**
- Modify: `src/io/readswap.f90`

- [ ] **Step 1: Locate the line**

```
grep -n "integer SwMacro\|integer swmacro" src/io/readswap.f90
```

Expected: one hit, around line 20, in the local declaration block of subroutine `readswap`.

- [ ] **Step 2: Delete the line**

The line is exactly:
```fortran
      integer SwMacro
```

Delete it. The subroutine's `use variables` import (line 9) now resolves `swmacro` to the module global, and `rdsinr('swmacro', 0, 1, swmacro)` writes the global directly.

- [ ] **Step 3: Verify no other locals shadow other globals in this subroutine** (defensive)

```
grep -nE "^\s+(integer|real|character|logical)\s+\w" src/io/readswap.f90 | head -50
```

Look for any other local variable name that's also declared as a global in `variables`. If found, leave alone (out of scope for this task — only `swmacro` is the parity blocker).

- [ ] **Step 4: Build + full regression**

```
pixi run -e test check-full
```

Expected: **6/6 regression cases green**. The macropore case (3.macroporeflow) sets `swmacro=1` in its `.swp`; if anything was implicitly relying on `swmacro` being uninitialized (=0) in the legacy path, that case is the canary. If macropore regression goes red, the implicit-zero assumption was real and we need to investigate. **STOP and report if check-full goes red.**

- [ ] **Step 5: Commit**

```
git add src/io/readswap.f90
git commit -m "fix(io): remove integer SwMacro shadow in readswap"
```

---

### Task 3: `nrlevs` clobber in `drainage_config_t%finalize`

**Files:**
- Modify: `src/config/drainage_config.f90`
- Modify: `tests/unit/config/test_drainage_config.pf`

- [ ] **Step 1: Update finalize**

Open `src/config/drainage_config.f90`. Find the `drainage_config_finalize` subroutine. Currently it's a no-op. Replace the body with:

```fortran
   subroutine drainage_config_finalize(self, errors)
      class(drainage_config_t), intent(inout) :: self
      type(error_collection_t), intent(inout) :: errors

      ! Mirror legacy convention from readswap.f90: single-level drainage
      ! methods (dramet 1 or 2) clobber nrlevs to 1 regardless of input.
      ! Matching this is required for parity with the legacy reader. The
      ! validator already constrains nrlevs in [0, 5]; this finalize step
      ! lands AFTER validate.
      if (self%dramet /= 3) self%nrlevs = 1
   end subroutine drainage_config_finalize
```

- [ ] **Step 2: Add tests**

Append to `tests/unit/config/test_drainage_config.pf`:

```fortran
@test
subroutine test_drainage_finalize_clobbers_nrlevs_when_dramet_not_3()
   use funit
   use error_mod
   use drainage_config_mod
   type(drainage_config_t)  :: d
   type(error_collection_t) :: errors
   d%dramet = 2
   d%nrlevs = 0
   call d%finalize(errors)
   @assertEqual(1, d%nrlevs)
end subroutine

@test
subroutine test_drainage_finalize_clobbers_nrlevs_even_if_user_set()
   use funit
   use error_mod
   use drainage_config_mod
   type(drainage_config_t)  :: d
   type(error_collection_t) :: errors
   d%dramet = 1
   d%nrlevs = 4
   call d%finalize(errors)
   @assertEqual(1, d%nrlevs)
end subroutine

@test
subroutine test_drainage_finalize_preserves_nrlevs_when_dramet_is_3()
   use funit
   use error_mod
   use drainage_config_mod
   type(drainage_config_t)  :: d
   type(error_collection_t) :: errors
   d%dramet = 3
   d%nrlevs = 5
   call d%finalize(errors)
   @assertEqual(5, d%nrlevs)
end subroutine
```

No inline `!` on assertion lines. Three new tests.

- [ ] **Step 3: Run**

```
pixi run -e test test-pfunit
```

Expected: 109 + 3 = 112 tests passing.

- [ ] **Step 4: Commit**

```
git add src/config/drainage_config.f90 tests/unit/config/test_drainage_config.pf
git commit -m "feat(config): nrlevs clobber to 1 when dramet/=3 in drainage finalize"
```

---

## Stage B — Path helpers (Task 4)

### Task 4: `path_helpers_mod`

**Files:**
- Create: `src/io/toml/path_helpers.f90`
- Create: `tests/unit/io/toml/test_path_helpers.pf`
- Modify: `meson.build`, `tests/unit/meson.build`, `tests/unit/testSuites.inc`

- [ ] **Step 1: Write the module**

Create `src/io/toml/path_helpers.f90`:

```fortran
!> Path-string helpers for cross-file TOML references.
!!
!! Pure string manipulation — no filesystem calls. The loader uses
!! `directory_of` to extract a base directory from a TOML file path,
!! then `resolve_relative_path` to resolve cross-file references against
!! that base.
module path_helpers_mod
   implicit none
   private

   public :: directory_of
   public :: resolve_relative_path

contains

   !> Return the directory portion of `path`, including the trailing slash.
   !! If `path` has no slash, returns "./".
   function directory_of(path) result(dir)
      character(len=*),              intent(in)  :: path
      character(len=:), allocatable              :: dir
      integer :: i
      i = index(path, '/', back=.true.)
      if (i == 0) then
         dir = './'
      else
         dir = path(1:i)
      end if
   end function directory_of

   !> Resolve `rel` against `base`. If `rel` starts with `/` it is treated
   !! as absolute and returned verbatim. Otherwise returns `base // rel`.
   !! `base` is expected to end with a slash (use `directory_of`).
   function resolve_relative_path(base, rel) result(abs)
      character(len=*),              intent(in)  :: base, rel
      character(len=:), allocatable              :: abs
      if (len_trim(rel) > 0 .and. rel(1:1) == '/') then
         abs = trim(rel)
      else
         abs = trim(base) // trim(rel)
      end if
   end function resolve_relative_path

end module path_helpers_mod
```

- [ ] **Step 2: Tests**

Create `tests/unit/io/toml/test_path_helpers.pf`:

```fortran
@test
subroutine test_directory_of_extracts_dir()
   use funit
   use path_helpers_mod
   character(len=:), allocatable :: d
   d = directory_of('tests/swap-cases/toml/1.hupselbrook/swap.toml')
   @assertEqual('tests/swap-cases/toml/1.hupselbrook/', d)
end subroutine

@test
subroutine test_directory_of_no_slash_returns_dot()
   use funit
   use path_helpers_mod
   character(len=:), allocatable :: d
   d = directory_of('swap.toml')
   @assertEqual('./', d)
end subroutine

@test
subroutine test_resolve_relative_concats()
   use funit
   use path_helpers_mod
   character(len=:), allocatable :: p
   p = resolve_relative_path('tests/swap-cases/toml/1.hupselbrook/', 'maizes.crp.toml')
   @assertEqual('tests/swap-cases/toml/1.hupselbrook/maizes.crp.toml', p)
end subroutine

@test
subroutine test_resolve_relative_absolute_passes_through()
   use funit
   use path_helpers_mod
   character(len=:), allocatable :: p
   p = resolve_relative_path('whatever/', '/abs/path/foo.toml')
   @assertEqual('/abs/path/foo.toml', p)
end subroutine
```

- [ ] **Step 3: Wire**

In `meson.build` `sources` list, append `'src/io/toml/path_helpers.f90'` (after the existing `src/io/toml/*` entries).

In `tests/unit/meson.build`:
- Append `'../../src/io/toml/path_helpers.f90'` to `pfunit_extra_sources`.
- Append `'io/toml/test_path_helpers.pf'` to `pf_files`.

In `tests/unit/testSuites.inc`: append `ADD_TEST_SUITE(test_path_helpers_suite)`.

- [ ] **Step 4: Run**

```
pixi run -e test test-pfunit
```

Expected: 112 + 4 = 116 tests passing.

- [ ] **Step 5: Commit**

```
git add src/io/toml/path_helpers.f90 \
        tests/unit/io/toml/test_path_helpers.pf \
        meson.build tests/unit/meson.build tests/unit/testSuites.inc
git commit -m "feat(io/toml): path_helpers (directory_of, resolve_relative_path)"
```

---

## Stage C — New crop config types (Tasks 5-6)

These types start with a CORE field set covering the common categories (~30 fields each). Additional fields are added during parity test iteration in Stage G when assertions reveal a gap. This is the same iterative pattern Phase 4b used.

### Task 5: `cropfixed_config_t`

**Files:**
- Create: `src/config/cropfixed_config.f90`
- Create: `tests/unit/config/test_cropfixed_config.pf`
- Modify: `meson.build`, `tests/unit/meson.build`, `tests/unit/testSuites.inc`

- [ ] **Step 1: Write the module**

Create `src/config/cropfixed_config.f90`:

```fortran
!> Type 1 (fixed/simple) crop config — populated from a .crp.toml file.
!! Field set is the COMMON SUBSET shared with type 3 (grass), plus type 1
!! specific fields. Phase 4c-a starts with a core schema; additional fields
!! are added during parity-test iteration as needed.
module cropfixed_config_mod
   use iso_fortran_env, only: real64
   use error_mod, only: error_collection_t, ERR_VALIDATION_OUT_OF_RANGE
   use validation_mod, only: check_int_enum, check_int_range, check_real_range, &
                             check_nonnegative_real, check_ordered_pair
   implicit none
   private

   public :: cropfixed_config_t

   type :: cropfixed_config_t
      ! Phenology
      integer      :: idev = 1                  !! 1=fixed period, 2=temperature-sum-based
      integer      :: lcc  = 0                  !! Length of crop cycle (days), used when idev=1

      ! Light & growth (when idev=2 path is taken; safe defaults for idev=1)
      real(real64) :: kdif = 0.0_real64         !! Diffuse light extinction
      real(real64) :: kdir = 0.0_real64         !! Direct light extinction

      ! Crop factor & height tables (stored as flat real arrays;
      ! pairs of (development_stage, value)).
      real(real64), allocatable :: cftb(:)      !! Crop factor table
      real(real64), allocatable :: chtb(:)      !! Crop height table

      ! Root growth
      real(real64) :: rdi = 0.0_real64          !! Initial rooting depth (cm)
      real(real64) :: rri = 0.0_real64          !! Daily root extension rate (cm/d)
      real(real64) :: rdc = 0.0_real64          !! Maximum rooting depth (cm)
      real(real64), allocatable :: rdctb(:)     !! Root density distribution table

      ! Water stress (Feddes)
      real(real64) :: hlim1 = 0.0_real64
      real(real64) :: hlim2u = 0.0_real64
      real(real64) :: hlim2l = 0.0_real64
      real(real64) :: hlim3h = 0.0_real64
      real(real64) :: hlim3l = 0.0_real64
      real(real64) :: hlim4 = 0.0_real64
      real(real64) :: adcrh = 0.0_real64
      real(real64) :: adcrl = 0.0_real64
      real(real64) :: rsc   = 0.0_real64        !! Crop resistance for ET method

      ! Salinity stress
      real(real64) :: ecmax  = 0.0_real64
      real(real64) :: ecslop = 0.0_real64

      ! Interception
      real(real64) :: cofab = 0.0_real64
   contains
      procedure :: validate => cropfixed_config_validate
      procedure :: finalize => cropfixed_config_finalize
   end type cropfixed_config_t

contains

   subroutine cropfixed_config_validate(self, errors)
      class(cropfixed_config_t), intent(in)    :: self
      type(error_collection_t),  intent(inout) :: errors

      call check_int_enum(self%idev, [1, 2], 'cropfixed.idev', errors)
      if (self%idev == 1) then
         call check_int_range(self%lcc, 1, 366, 'cropfixed.lcc', errors)
      end if

      call check_real_range(self%rdi, 0.0_real64, 1000.0_real64, 'cropfixed.rdi', errors)
      call check_real_range(self%rdc, 0.0_real64, 1000.0_real64, 'cropfixed.rdc', errors)
      call check_nonnegative_real(self%rri, 'cropfixed.rri', errors)

      ! Feddes: hlim1 (saturation, near 0) > hlim2u > hlim2l > hlim3h > hlim3l > hlim4 (wilting)
      ! Permit 0 (the validator below catches "all zero" defaults via cross-field rule).
      call check_ordered_pair(self%hlim2l, self%hlim2u, 'hlim2l', 'hlim2u', 'cropfixed', errors)
      call check_ordered_pair(self%hlim3l, self%hlim3h, 'hlim3l', 'hlim3h', 'cropfixed', errors)

      call check_nonnegative_real(self%ecmax,  'cropfixed.ecmax',  errors)
      call check_nonnegative_real(self%ecslop, 'cropfixed.ecslop', errors)
   end subroutine cropfixed_config_validate

   subroutine cropfixed_config_finalize(self, errors)
      class(cropfixed_config_t), intent(inout) :: self
      type(error_collection_t),  intent(inout) :: errors
      ! No derived fields in 4c-a.
   end subroutine cropfixed_config_finalize

end module cropfixed_config_mod
```

- [ ] **Step 2: Tests**

Create `tests/unit/config/test_cropfixed_config.pf`:

```fortran
@test
subroutine test_cropfixed_default_invalid_idev2_no_lcc_check()
   use funit
   use error_mod
   use cropfixed_config_mod
   type(cropfixed_config_t) :: c
   type(error_collection_t) :: errors
   c%idev = 2
   call c%validate(errors)
   @assertFalse(errors%has_errors())
end subroutine

@test
subroutine test_cropfixed_idev1_requires_lcc_in_range()
   use funit
   use error_mod
   use cropfixed_config_mod
   type(cropfixed_config_t) :: c
   type(error_collection_t) :: errors
   c%idev = 1
   c%lcc  = 0
   call c%validate(errors)
   @assertTrue(errors%has_errors())
end subroutine

@test
subroutine test_cropfixed_idev1_lcc_180_passes()
   use funit
   use error_mod
   use cropfixed_config_mod
   type(cropfixed_config_t) :: c
   type(error_collection_t) :: errors
   c%idev = 1
   c%lcc  = 180
   call c%validate(errors)
   @assertFalse(errors%has_errors())
end subroutine

@test
subroutine test_cropfixed_invalid_idev_fails()
   use funit
   use error_mod
   use cropfixed_config_mod
   type(cropfixed_config_t) :: c
   type(error_collection_t) :: errors
   c%idev = 9
   call c%validate(errors)
   @assertTrue(errors%has_errors())
end subroutine

@test
subroutine test_cropfixed_negative_rdi_fails()
   use funit
   use error_mod
   use cropfixed_config_mod
   type(cropfixed_config_t) :: c
   type(error_collection_t) :: errors
   c%idev = 2
   c%rdi  = -1.0d0
   call c%validate(errors)
   @assertTrue(errors%has_errors())
end subroutine

@test
subroutine test_cropfixed_finalize_idempotent()
   use funit
   use error_mod
   use cropfixed_config_mod
   type(cropfixed_config_t) :: c
   type(error_collection_t) :: errors
   c%idev = 2
   call c%finalize(errors)
   call c%finalize(errors)
   @assertFalse(errors%has_errors())
end subroutine
```

- [ ] **Step 3: Wire**

`meson.build` `sources` list, append `'src/config/cropfixed_config.f90'` (after the existing `src/config/*` entries, before `swap_config.f90`).

`tests/unit/meson.build`:
- Append `'../../src/config/cropfixed_config.f90'` to `pfunit_extra_sources`.
- Append `'config/test_cropfixed_config.pf'` to `pf_files`.

`tests/unit/testSuites.inc`: append `ADD_TEST_SUITE(test_cropfixed_config_suite)`.

- [ ] **Step 4: Run**

```
pixi run -e test test-pfunit
```

Expected: 116 + 6 = 122 tests passing.

- [ ] **Step 5: Commit**

```
git add src/config/cropfixed_config.f90 \
        tests/unit/config/test_cropfixed_config.pf \
        meson.build tests/unit/meson.build tests/unit/testSuites.inc
git commit -m "feat(config): cropfixed_config_t (type 1 crop schema, core fields)"
```

---

### Task 6: `cropgrass_config_t`

**Files:**
- Create: `src/config/cropgrass_config.f90`
- Create: `tests/unit/config/test_cropgrass_config.pf`
- Modify: `meson.build`, `tests/unit/meson.build`, `tests/unit/testSuites.inc`

- [ ] **Step 1: Write the module**

`cropgrass_config_t` mirrors `cropfixed_config_t`'s common subset (Feddes, salinity, root growth, interception) plus grass-specific extensions. Per Q2 decision A, no shared base — each type declares its own complete field set.

Create `src/config/cropgrass_config.f90`:

```fortran
!> Type 3 (grass / WOFOST grass) crop config — populated from a .crp.toml file.
!! Field set is similar to cropfixed_config_t plus grass-specific management
!! (mowing, grazing, fertilizer). Per design Q2A, no shared base type.
module cropgrass_config_mod
   use iso_fortran_env, only: real64
   use error_mod, only: error_collection_t, ERR_VALIDATION_OUT_OF_RANGE
   use validation_mod, only: check_int_enum, check_int_range, check_real_range, &
                             check_nonnegative_real, check_ordered_pair
   implicit none
   private

   public :: cropgrass_config_t

   type :: cropgrass_config_t
      ! Phenology / development
      integer      :: idev = 2                  !! For grass typically 2 (temp-sum-based)
      integer      :: lcc  = 0
      real(real64) :: tbase = 0.0_real64
      real(real64) :: tsum1 = 0.0_real64
      real(real64) :: tsum2 = 0.0_real64

      ! Light & growth
      real(real64) :: kdif = 0.0_real64
      real(real64) :: kdir = 0.0_real64
      real(real64) :: eff  = 0.0_real64
      real(real64) :: amax = 0.0_real64

      ! Tables
      real(real64), allocatable :: cftb(:)
      real(real64), allocatable :: chtb(:)
      real(real64), allocatable :: rdctb(:)

      ! Root growth
      real(real64) :: rdi = 0.0_real64
      real(real64) :: rri = 0.0_real64
      real(real64) :: rdc = 0.0_real64

      ! Water stress (Feddes) — same as fixed
      real(real64) :: hlim1 = 0.0_real64
      real(real64) :: hlim2u = 0.0_real64
      real(real64) :: hlim2l = 0.0_real64
      real(real64) :: hlim3h = 0.0_real64
      real(real64) :: hlim3l = 0.0_real64
      real(real64) :: hlim4 = 0.0_real64
      real(real64) :: adcrh = 0.0_real64
      real(real64) :: adcrl = 0.0_real64
      real(real64) :: rsc   = 0.0_real64

      ! Salinity & interception
      real(real64) :: ecmax  = 0.0_real64
      real(real64) :: ecslop = 0.0_real64
      real(real64) :: cofab  = 0.0_real64

      ! Mowing schedule (grass-specific)
      integer :: swharv = 0                     !! 0=no scheduled mow, 1=scheduled
      integer :: nmow   = 0                     !! Number of mowing events
      real(real64), allocatable :: dates_mowing(:)
      real(real64), allocatable :: lai_after_mow(:)

      ! Grazing (grass-specific)
      integer :: swgraz = 0                     !! 0=no grazing, 1=scheduled
      real(real64) :: nstart_graz = 0.0_real64  !! Start day-of-year
      real(real64) :: nstop_graz  = 0.0_real64  !! Stop day-of-year
   contains
      procedure :: validate => cropgrass_config_validate
      procedure :: finalize => cropgrass_config_finalize
   end type cropgrass_config_t

contains

   subroutine cropgrass_config_validate(self, errors)
      class(cropgrass_config_t), intent(in)    :: self
      type(error_collection_t),  intent(inout) :: errors

      call check_int_enum(self%idev,   [1, 2],    'cropgrass.idev',   errors)
      call check_int_enum(self%swharv, [0, 1],    'cropgrass.swharv', errors)
      call check_int_enum(self%swgraz, [0, 1],    'cropgrass.swgraz', errors)

      call check_real_range(self%rdi, 0.0_real64, 1000.0_real64, 'cropgrass.rdi', errors)
      call check_real_range(self%rdc, 0.0_real64, 1000.0_real64, 'cropgrass.rdc', errors)
      call check_nonnegative_real(self%rri, 'cropgrass.rri', errors)

      call check_ordered_pair(self%hlim2l, self%hlim2u, 'hlim2l', 'hlim2u', 'cropgrass', errors)
      call check_ordered_pair(self%hlim3l, self%hlim3h, 'hlim3l', 'hlim3h', 'cropgrass', errors)

      call check_nonnegative_real(self%ecmax,  'cropgrass.ecmax',  errors)
      call check_nonnegative_real(self%ecslop, 'cropgrass.ecslop', errors)

      if (self%swharv == 1) then
         call check_int_range(self%nmow, 1, 50, 'cropgrass.nmow', errors)
      end if
   end subroutine cropgrass_config_validate

   subroutine cropgrass_config_finalize(self, errors)
      class(cropgrass_config_t), intent(inout) :: self
      type(error_collection_t),  intent(inout) :: errors
   end subroutine cropgrass_config_finalize

end module cropgrass_config_mod
```

- [ ] **Step 2: Tests**

Create `tests/unit/config/test_cropgrass_config.pf` with 6 tests in the same shape as `test_cropfixed_config.pf`. Cover: defaults pass with idev=2; swharv=1 needs nmow in range; swharv=0 ignores nmow; swgraz enum check; negative rdi fails; finalize idempotent.

```fortran
@test
subroutine test_cropgrass_defaults_pass()
   use funit
   use error_mod
   use cropgrass_config_mod
   type(cropgrass_config_t) :: c
   type(error_collection_t) :: errors
   call c%validate(errors)
   @assertFalse(errors%has_errors())
end subroutine

@test
subroutine test_cropgrass_swharv1_requires_nmow()
   use funit
   use error_mod
   use cropgrass_config_mod
   type(cropgrass_config_t) :: c
   type(error_collection_t) :: errors
   c%swharv = 1
   c%nmow   = 0
   call c%validate(errors)
   @assertTrue(errors%has_errors())
end subroutine

@test
subroutine test_cropgrass_swharv1_with_nmow_passes()
   use funit
   use error_mod
   use cropgrass_config_mod
   type(cropgrass_config_t) :: c
   type(error_collection_t) :: errors
   c%swharv = 1
   c%nmow   = 4
   call c%validate(errors)
   @assertFalse(errors%has_errors())
end subroutine

@test
subroutine test_cropgrass_swgraz_invalid_fails()
   use funit
   use error_mod
   use cropgrass_config_mod
   type(cropgrass_config_t) :: c
   type(error_collection_t) :: errors
   c%swgraz = 9
   call c%validate(errors)
   @assertTrue(errors%has_errors())
end subroutine

@test
subroutine test_cropgrass_negative_rdi_fails()
   use funit
   use error_mod
   use cropgrass_config_mod
   type(cropgrass_config_t) :: c
   type(error_collection_t) :: errors
   c%rdi = -1.0d0
   call c%validate(errors)
   @assertTrue(errors%has_errors())
end subroutine

@test
subroutine test_cropgrass_finalize_idempotent()
   use funit
   use error_mod
   use cropgrass_config_mod
   type(cropgrass_config_t) :: c
   type(error_collection_t) :: errors
   call c%finalize(errors)
   call c%finalize(errors)
   @assertFalse(errors%has_errors())
end subroutine
```

- [ ] **Step 3: Wire** (same pattern as Task 5).

- [ ] **Step 4: Run**

```
pixi run -e test test-pfunit
```

Expected: 122 + 6 = 128 tests passing.

- [ ] **Step 5: Commit**

```
git add src/config/cropgrass_config.f90 \
        tests/unit/config/test_cropgrass_config.pf \
        meson.build tests/unit/meson.build tests/unit/testSuites.inc
git commit -m "feat(config): cropgrass_config_t (type 3 crop schema, core fields)"
```

---

## Stage D — New crop section readers (Tasks 7-8)

### Task 7: `read_cropfixed_toml`

**Files:**
- Create: `src/io/toml/read_cropfixed_toml.f90`
- Create: `tests/unit/io/toml/test_read_cropfixed_toml.pf`
- Create: `tests/unit/io/toml/fixtures/cropfixed_minimal.toml`
- Modify: `meson.build`, `tests/unit/meson.build`, `tests/unit/testSuites.inc`

- [ ] **Step 1: Write the reader**

Create `src/io/toml/read_cropfixed_toml.f90`:

```fortran
!> Reader for a type-1 (fixed crop) .crp.toml file.
!!
!! The reader is given the ROOT table of a loaded .crp.toml document
!! (loaded by read_crop_toml when it follows a [[crop.rotation]].file
!! reference). It populates a cropfixed_config_t.
!!
!! TOML schema follows the legacy section grouping where it makes sense:
!!     [phenology]    idev, lcc
!!     [light]        kdif, kdir
!!     [tables]       cftb (array of pairs), chtb, rdctb
!!     [root]         rdi, rri, rdc
!!     [water_stress] hlim1, hlim2u, hlim2l, hlim3h, hlim3l, hlim4,
!!                    adcrh, adcrl, rsc
!!     [salinity]     ecmax, ecslop
!!     [interception] cofab
module read_cropfixed_toml_mod
   use iso_fortran_env, only: real64
   use tomlf, only: toml_table
   use cropfixed_config_mod, only: cropfixed_config_t
   use toml_field_helpers_mod, only: get_table,                    &
                                     get_optional_int_with_default, &
                                     get_optional_real_with_default
   use error_mod, only: error_collection_t
   implicit none
   private

   public :: read_cropfixed_toml

contains

   subroutine read_cropfixed_toml(doc, config, errors)
      type(toml_table), pointer,  intent(in)    :: doc
      type(cropfixed_config_t),   intent(inout) :: config
      type(error_collection_t),   intent(inout) :: errors

      type(toml_table), pointer :: ph, light, root, ws, salt, inter

      call get_table(doc, 'phenology', ph, 'phenology', errors)
      if (associated(ph)) then
         call get_optional_int_with_default(ph,  'idev', config%idev, 1, 'phenology.idev', errors)
         call get_optional_int_with_default(ph,  'lcc',  config%lcc,  0, 'phenology.lcc',  errors)
      end if

      call get_table(doc, 'light', light, 'light', errors)
      if (associated(light)) then
         call get_optional_real_with_default(light, 'kdif', config%kdif, 0.0_real64, 'light.kdif', errors)
         call get_optional_real_with_default(light, 'kdir', config%kdir, 0.0_real64, 'light.kdir', errors)
      end if

      call get_table(doc, 'root', root, 'root', errors)
      if (associated(root)) then
         call get_optional_real_with_default(root, 'rdi', config%rdi, 0.0_real64, 'root.rdi', errors)
         call get_optional_real_with_default(root, 'rri', config%rri, 0.0_real64, 'root.rri', errors)
         call get_optional_real_with_default(root, 'rdc', config%rdc, 0.0_real64, 'root.rdc', errors)
      end if

      call get_table(doc, 'water_stress', ws, 'water_stress', errors)
      if (associated(ws)) then
         call get_optional_real_with_default(ws, 'hlim1',  config%hlim1,  0.0_real64, 'ws.hlim1',  errors)
         call get_optional_real_with_default(ws, 'hlim2u', config%hlim2u, 0.0_real64, 'ws.hlim2u', errors)
         call get_optional_real_with_default(ws, 'hlim2l', config%hlim2l, 0.0_real64, 'ws.hlim2l', errors)
         call get_optional_real_with_default(ws, 'hlim3h', config%hlim3h, 0.0_real64, 'ws.hlim3h', errors)
         call get_optional_real_with_default(ws, 'hlim3l', config%hlim3l, 0.0_real64, 'ws.hlim3l', errors)
         call get_optional_real_with_default(ws, 'hlim4',  config%hlim4,  0.0_real64, 'ws.hlim4',  errors)
         call get_optional_real_with_default(ws, 'adcrh',  config%adcrh,  0.0_real64, 'ws.adcrh',  errors)
         call get_optional_real_with_default(ws, 'adcrl',  config%adcrl,  0.0_real64, 'ws.adcrl',  errors)
         call get_optional_real_with_default(ws, 'rsc',    config%rsc,    0.0_real64, 'ws.rsc',    errors)
      end if

      call get_table(doc, 'salinity', salt, 'salinity', errors)
      if (associated(salt)) then
         call get_optional_real_with_default(salt, 'ecmax',  config%ecmax,  0.0_real64, 'salt.ecmax',  errors)
         call get_optional_real_with_default(salt, 'ecslop', config%ecslop, 0.0_real64, 'salt.ecslop', errors)
      end if

      call get_table(doc, 'interception', inter, 'interception', errors)
      if (associated(inter)) then
         call get_optional_real_with_default(inter, 'cofab', config%cofab, 0.0_real64, 'inter.cofab', errors)
      end if
   end subroutine read_cropfixed_toml

end module read_cropfixed_toml_mod
```

Note: tables (cftb, chtb, rdctb) are NOT read in 4c-a — they're added during parity iteration when needed.

- [ ] **Step 2: Fixture**

Create `tests/unit/io/toml/fixtures/cropfixed_minimal.toml`:

```toml
[phenology]
idev = 1
lcc  = 168

[root]
rdi = 5.0
rri = 1.2
rdc = 100.0

[water_stress]
hlim1  = -10.0
hlim2u = -25.0
hlim2l = -200.0
hlim3h = -400.0
hlim3l = -600.0
hlim4  = -8000.0
adcrh  = 0.5
adcrl  = 0.1
rsc    = 70.0

[salinity]
ecmax  = 1.7
ecslop = 12.0

[interception]
cofab = 0.25
```

- [ ] **Step 3: Tests**

Create `tests/unit/io/toml/test_read_cropfixed_toml.pf`:

```fortran
@test
subroutine test_read_cropfixed_happy()
   use funit
   use tomlf, only: toml_table, toml_load
   use cropfixed_config_mod
   use read_cropfixed_toml_mod
   use error_mod
   type(toml_table), allocatable, target :: doc
   type(toml_table), pointer             :: doc_ptr
   type(cropfixed_config_t)              :: c
   type(error_collection_t)              :: errors

   call toml_load(doc, 'tests/unit/io/toml/fixtures/cropfixed_minimal.toml')
   doc_ptr => doc
   call read_cropfixed_toml(doc_ptr, c, errors)

   @assertFalse(errors%has_errors())
   @assertEqual(1,       c%idev)
   @assertEqual(168,     c%lcc)
   @assertEqual(5.0d0,   c%rdi,    1.0d-12)
   @assertEqual(-200.0d0, c%hlim2l, 1.0d-12)
   @assertEqual(70.0d0,  c%rsc,    1.0d-12)
   @assertEqual(0.25d0,  c%cofab,  1.0d-12)
end subroutine

@test
subroutine test_read_cropfixed_missing_section_silent()
   use funit
   use tomlf, only: toml_table, toml_load
   use cropfixed_config_mod
   use read_cropfixed_toml_mod
   use error_mod
   type(toml_table), allocatable, target :: doc
   type(toml_table), pointer             :: doc_ptr
   type(cropfixed_config_t)              :: c
   type(error_collection_t)              :: errors

   call toml_load(doc, 'tests/unit/io/toml/fixtures/general_minimal.toml')
   doc_ptr => doc
   call read_cropfixed_toml(doc_ptr, c, errors)

   @assertFalse(errors%has_errors())
end subroutine
```

- [ ] **Step 4: Wire**

`meson.build`: append `'src/io/toml/read_cropfixed_toml.f90'` to `sources`.

`tests/unit/meson.build`:
- Append `'../../src/io/toml/read_cropfixed_toml.f90'` to `pfunit_extra_sources`.
- Append `'io/toml/test_read_cropfixed_toml.pf'` to `pf_files`.

`tests/unit/testSuites.inc`: append `ADD_TEST_SUITE(test_read_cropfixed_toml_suite)`.

- [ ] **Step 5: Run**

```
pixi run -e test test-pfunit
```

Expected: 128 + 2 = 130 tests passing.

- [ ] **Step 6: Commit**

```
git add src/io/toml/read_cropfixed_toml.f90 \
        tests/unit/io/toml/test_read_cropfixed_toml.pf \
        tests/unit/io/toml/fixtures/cropfixed_minimal.toml \
        meson.build tests/unit/meson.build tests/unit/testSuites.inc
git commit -m "feat(io/toml): read_cropfixed_toml section reader"
```

---

### Task 8: `read_cropgrass_toml`

Same pattern as Task 7. Schema includes the `[mowing]` and `[grazing]` sections in addition to the shared ones.

**Files:**
- Create: `src/io/toml/read_cropgrass_toml.f90`, `tests/unit/io/toml/test_read_cropgrass_toml.pf`, `tests/unit/io/toml/fixtures/cropgrass_minimal.toml`
- Modify: `meson.build`, `tests/unit/meson.build`, `tests/unit/testSuites.inc`

- [ ] **Step 1: Reader**

```fortran
module read_cropgrass_toml_mod
   use iso_fortran_env, only: real64
   use tomlf, only: toml_table
   use cropgrass_config_mod, only: cropgrass_config_t
   use toml_field_helpers_mod, only: get_table,                    &
                                     get_optional_int_with_default, &
                                     get_optional_real_with_default
   use error_mod, only: error_collection_t
   implicit none
   private

   public :: read_cropgrass_toml

contains

   subroutine read_cropgrass_toml(doc, config, errors)
      type(toml_table), pointer,  intent(in)    :: doc
      type(cropgrass_config_t),   intent(inout) :: config
      type(error_collection_t),   intent(inout) :: errors

      type(toml_table), pointer :: ph, light, root, ws, salt, inter, mow, graz

      call get_table(doc, 'phenology', ph, 'phenology', errors)
      if (associated(ph)) then
         call get_optional_int_with_default(ph,  'idev',  config%idev,  2, 'phenology.idev',  errors)
         call get_optional_int_with_default(ph,  'lcc',   config%lcc,   0, 'phenology.lcc',   errors)
         call get_optional_real_with_default(ph, 'tbase', config%tbase, 0.0_real64, 'phenology.tbase', errors)
         call get_optional_real_with_default(ph, 'tsum1', config%tsum1, 0.0_real64, 'phenology.tsum1', errors)
         call get_optional_real_with_default(ph, 'tsum2', config%tsum2, 0.0_real64, 'phenology.tsum2', errors)
      end if

      call get_table(doc, 'light', light, 'light', errors)
      if (associated(light)) then
         call get_optional_real_with_default(light, 'kdif', config%kdif, 0.0_real64, 'light.kdif', errors)
         call get_optional_real_with_default(light, 'kdir', config%kdir, 0.0_real64, 'light.kdir', errors)
         call get_optional_real_with_default(light, 'eff',  config%eff,  0.0_real64, 'light.eff',  errors)
         call get_optional_real_with_default(light, 'amax', config%amax, 0.0_real64, 'light.amax', errors)
      end if

      call get_table(doc, 'root', root, 'root', errors)
      if (associated(root)) then
         call get_optional_real_with_default(root, 'rdi', config%rdi, 0.0_real64, 'root.rdi', errors)
         call get_optional_real_with_default(root, 'rri', config%rri, 0.0_real64, 'root.rri', errors)
         call get_optional_real_with_default(root, 'rdc', config%rdc, 0.0_real64, 'root.rdc', errors)
      end if

      call get_table(doc, 'water_stress', ws, 'water_stress', errors)
      if (associated(ws)) then
         call get_optional_real_with_default(ws, 'hlim1',  config%hlim1,  0.0_real64, 'ws.hlim1',  errors)
         call get_optional_real_with_default(ws, 'hlim2u', config%hlim2u, 0.0_real64, 'ws.hlim2u', errors)
         call get_optional_real_with_default(ws, 'hlim2l', config%hlim2l, 0.0_real64, 'ws.hlim2l', errors)
         call get_optional_real_with_default(ws, 'hlim3h', config%hlim3h, 0.0_real64, 'ws.hlim3h', errors)
         call get_optional_real_with_default(ws, 'hlim3l', config%hlim3l, 0.0_real64, 'ws.hlim3l', errors)
         call get_optional_real_with_default(ws, 'hlim4',  config%hlim4,  0.0_real64, 'ws.hlim4',  errors)
         call get_optional_real_with_default(ws, 'adcrh',  config%adcrh,  0.0_real64, 'ws.adcrh',  errors)
         call get_optional_real_with_default(ws, 'adcrl',  config%adcrl,  0.0_real64, 'ws.adcrl',  errors)
         call get_optional_real_with_default(ws, 'rsc',    config%rsc,    0.0_real64, 'ws.rsc',    errors)
      end if

      call get_table(doc, 'salinity', salt, 'salinity', errors)
      if (associated(salt)) then
         call get_optional_real_with_default(salt, 'ecmax',  config%ecmax,  0.0_real64, 'salt.ecmax',  errors)
         call get_optional_real_with_default(salt, 'ecslop', config%ecslop, 0.0_real64, 'salt.ecslop', errors)
      end if

      call get_table(doc, 'interception', inter, 'interception', errors)
      if (associated(inter)) then
         call get_optional_real_with_default(inter, 'cofab', config%cofab, 0.0_real64, 'inter.cofab', errors)
      end if

      call get_table(doc, 'mowing', mow, 'mowing', errors)
      if (associated(mow)) then
         call get_optional_int_with_default(mow, 'swharv', config%swharv, 0, 'mowing.swharv', errors)
         call get_optional_int_with_default(mow, 'nmow',   config%nmow,   0, 'mowing.nmow',   errors)
      end if

      call get_table(doc, 'grazing', graz, 'grazing', errors)
      if (associated(graz)) then
         call get_optional_int_with_default(graz, 'swgraz',      config%swgraz,      0,         'grazing.swgraz',      errors)
         call get_optional_real_with_default(graz, 'nstart_graz', config%nstart_graz, 0.0_real64, 'grazing.nstart_graz', errors)
         call get_optional_real_with_default(graz, 'nstop_graz',  config%nstop_graz,  0.0_real64, 'grazing.nstop_graz',  errors)
      end if
   end subroutine read_cropgrass_toml

end module read_cropgrass_toml_mod
```

- [ ] **Step 2: Fixture**

`tests/unit/io/toml/fixtures/cropgrass_minimal.toml`:

```toml
[phenology]
idev = 2
tsum1 = 800.0
tsum2 = 800.0

[root]
rdi = 5.0
rri = 1.0
rdc = 50.0

[water_stress]
hlim1  = -10.0
hlim2u = -25.0
hlim2l = -200.0
hlim3h = -300.0
hlim3l = -500.0
hlim4  = -8000.0

[mowing]
swharv = 1
nmow   = 5

[grazing]
swgraz = 0
```

- [ ] **Step 3: Tests** (2 tests, same shape as cropfixed)

```fortran
@test
subroutine test_read_cropgrass_happy()
   use funit
   use tomlf, only: toml_table, toml_load
   use cropgrass_config_mod
   use read_cropgrass_toml_mod
   use error_mod
   type(toml_table), allocatable, target :: doc
   type(toml_table), pointer             :: doc_ptr
   type(cropgrass_config_t)              :: c
   type(error_collection_t)              :: errors

   call toml_load(doc, 'tests/unit/io/toml/fixtures/cropgrass_minimal.toml')
   doc_ptr => doc
   call read_cropgrass_toml(doc_ptr, c, errors)

   @assertFalse(errors%has_errors())
   @assertEqual(2,    c%idev)
   @assertEqual(1,    c%swharv)
   @assertEqual(5,    c%nmow)
   @assertEqual(800.0d0, c%tsum1, 1.0d-12)
   @assertEqual(50.0d0,  c%rdc,   1.0d-12)
end subroutine

@test
subroutine test_read_cropgrass_missing_section_silent()
   use funit
   use tomlf, only: toml_table, toml_load
   use cropgrass_config_mod
   use read_cropgrass_toml_mod
   use error_mod
   type(toml_table), allocatable, target :: doc
   type(toml_table), pointer             :: doc_ptr
   type(cropgrass_config_t)              :: c
   type(error_collection_t)              :: errors

   call toml_load(doc, 'tests/unit/io/toml/fixtures/general_minimal.toml')
   doc_ptr => doc
   call read_cropgrass_toml(doc_ptr, c, errors)

   @assertFalse(errors%has_errors())
end subroutine
```

- [ ] **Step 4: Wire** + run.

Expected: 130 + 2 = 132 tests passing.

- [ ] **Step 5: Commit**

```
git add src/io/toml/read_cropgrass_toml.f90 \
        tests/unit/io/toml/test_read_cropgrass_toml.pf \
        tests/unit/io/toml/fixtures/cropgrass_minimal.toml \
        meson.build tests/unit/meson.build tests/unit/testSuites.inc
git commit -m "feat(io/toml): read_cropgrass_toml section reader"
```

---

## Stage E — Cross-file loader (Tasks 9-10)

### Task 9: Extend `read_drainage_toml` with `file=` follow

**Files:**
- Modify: `src/io/toml/read_drainage_toml.f90`
- Modify: `tests/unit/io/toml/test_read_drainage_toml.pf`
- Create: `tests/unit/io/toml/fixtures/drainage_via_file.toml`, `tests/unit/io/toml/fixtures/drainage_external.dra.toml`

- [ ] **Step 1: Refactor reader to extract inner helper**

Open `src/io/toml/read_drainage_toml.f90`. The current `read_drainage_toml(doc, config, errors)` walks `[drainage]` directly. Refactor:

1. Add an optional `base_path` argument: `read_drainage_toml(doc, config, errors, base_path)`.
2. Extract the `[drainage]`-walking body into a `read_drainage_inner(drain_table, config, errors)` private subroutine.
3. The new top-level reader: get `[drainage]` table, check for `file = "..."` key. If present and `base_path` is provided, load the referenced TOML and call `read_drainage_inner` on its `[drainage]` table. Otherwise call `read_drainage_inner` on the original table.

```fortran
module read_drainage_toml_mod
   use iso_fortran_env, only: real64
   use tomlf, only: toml_table, toml_array, toml_error, toml_load, get_value, len
   use drainage_config_mod, only: drainage_config_t
   use toml_field_helpers_mod, only: get_table, get_array_of_tables,   &
                                     get_optional_int_with_default,    &
                                     get_optional_real_with_default,   &
                                     get_optional_string_with_default
   use path_helpers_mod, only: resolve_relative_path
   use error_mod, only: error_collection_t, ERR_PARSE_MALFORMED_TOML
   implicit none
   private

   public :: read_drainage_toml

contains

   subroutine read_drainage_toml(doc, config, errors, base_path)
      type(toml_table), pointer,  intent(in)    :: doc
      type(drainage_config_t),    intent(inout) :: config
      type(error_collection_t),   intent(inout) :: errors
      character(len=*), optional, intent(in)    :: base_path

      type(toml_table), pointer             :: drain_tab, ext_tab
      type(toml_table), allocatable, target :: ext_doc
      type(toml_error), allocatable         :: terr
      character(len=:), allocatable         :: file_rel, file_abs

      call get_table(doc, 'drainage', drain_tab, 'drainage', errors)
      if (.not. associated(drain_tab)) return

      call get_optional_string_with_default(drain_tab, 'file', file_rel, '', 'drainage.file', errors)
      if (len_trim(file_rel) > 0 .and. present(base_path)) then
         file_abs = resolve_relative_path(base_path, file_rel)
         call toml_load(ext_doc, trim(file_abs), error=terr)
         if (allocated(terr)) then
            call errors%append(ERR_PARSE_MALFORMED_TOML, trim(terr%message), trim(file_abs))
            return
         end if
         call get_table(ext_doc, 'drainage', ext_tab, 'drainage', errors)
         if (.not. associated(ext_tab)) return
         call read_drainage_inner(ext_tab, config, errors)
      else
         call read_drainage_inner(drain_tab, config, errors)
      end if
   end subroutine read_drainage_toml

   subroutine read_drainage_inner(sec, config, errors)
      type(toml_table), pointer,  intent(in)    :: sec
      type(drainage_config_t),    intent(inout) :: config
      type(error_collection_t),   intent(inout) :: errors

      type(toml_table), pointer :: basic, item
      type(toml_array), pointer :: levels
      integer :: i, n, stat

      call get_optional_int_with_default(sec, 'swdra',    config%swdra,    0, 'drainage.swdra',    errors)
      call get_optional_int_with_default(sec, 'dramet',   config%dramet,   0, 'drainage.dramet',   errors)
      call get_optional_int_with_default(sec, 'swdivd',   config%swdivd,   0, 'drainage.swdivd',   errors)
      call get_optional_int_with_default(sec, 'swdislay', config%swdislay, 0, 'drainage.swdislay', errors)
      call get_optional_int_with_default(sec, 'nrlevs',   config%nrlevs,   0, 'drainage.nrlevs',   errors)
      call get_optional_real_with_default(sec, 'altcu',   config%altcu,    0.0_real64, 'drainage.altcu', errors)

      call get_table(sec, 'basic', basic, 'drainage.basic', errors)
      if (associated(basic)) then
         call get_optional_real_with_default(basic, 'basegw', config%basegw, 0.0_real64, 'drainage.basic.basegw', errors)
         call get_optional_real_with_default(basic, 'entres', config%entres, 0.0_real64, 'drainage.basic.entres', errors)
         call get_optional_real_with_default(basic, 'shape',  config%shape,  0.0_real64, 'drainage.basic.shape',  errors)
      end if

      call get_array_of_tables(sec, 'levels', levels, 'drainage.levels', errors)
      if (associated(levels)) then
         n = len(levels)
         if (n > 0) then
            ! existing per-level allocation + walk — preserve from prior implementation
            ! See git history of src/io/toml/read_drainage_toml.f90 for the verbatim block.
            allocate(config%swdtyp(n), config%zbotdr(n), config%drares(n), &
                     config%infres(n), config%L(n),      config%gwlinf(n), &
                     config%rdrain(n), config%rinfi(n),  config%rentry(n), &
                     config%rexit(n),  config%widthr(n), config%taludr(n), &
                     config%swallo(n))
            config%swdtyp  = 0
            config%zbotdr  = 0.0_real64
            config%drares  = 0.0_real64
            config%infres  = 0.0_real64
            config%L       = 0.0_real64
            config%gwlinf  = 0.0_real64
            config%rdrain  = 0.0_real64
            config%rinfi   = 0.0_real64
            config%rentry  = 0.0_real64
            config%rexit   = 0.0_real64
            config%widthr  = 0.0_real64
            config%taludr  = 0.0_real64
            config%swallo  = 0
            do i = 1, n
               call get_value(levels, i, item, stat=stat)
               if (stat /= 0 .or. .not. associated(item)) cycle
               call get_optional_int_with_default(item,  'swdtyp', config%swdtyp(i), 0,           'drainage.levels.swdtyp', errors)
               call get_optional_real_with_default(item, 'zbotdr', config%zbotdr(i), 0.0_real64,  'drainage.levels.zbotdr', errors)
               call get_optional_real_with_default(item, 'drares', config%drares(i), 0.0_real64,  'drainage.levels.drares', errors)
               call get_optional_real_with_default(item, 'infres', config%infres(i), 0.0_real64,  'drainage.levels.infres', errors)
               call get_optional_real_with_default(item, 'L',      config%L(i),      0.0_real64,  'drainage.levels.L',      errors)
            end do
         end if
      end if
   end subroutine read_drainage_inner

end module read_drainage_toml_mod
```

- [ ] **Step 2: Fixtures**

`tests/unit/io/toml/fixtures/drainage_via_file.toml`:
```toml
[drainage]
file = "drainage_external.dra.toml"
```

`tests/unit/io/toml/fixtures/drainage_external.dra.toml`:
```toml
[drainage]
swdra    = 1
dramet   = 1
swdivd   = 0
nrlevs   = 1
altcu    = 0.0

[drainage.basic]
basegw = -200.0
entres = 20.0
shape  = 0.8
```

- [ ] **Step 3: Test for cross-file path**

Append to `tests/unit/io/toml/test_read_drainage_toml.pf`:

```fortran
@test
subroutine test_read_drainage_follows_file_reference()
   use funit
   use tomlf, only: toml_table, toml_load
   use drainage_config_mod
   use read_drainage_toml_mod
   use error_mod
   type(toml_table), allocatable, target :: doc
   type(toml_table), pointer             :: doc_ptr
   type(drainage_config_t)               :: d
   type(error_collection_t)              :: errors

   call toml_load(doc, 'tests/unit/io/toml/fixtures/drainage_via_file.toml')
   doc_ptr => doc
   call read_drainage_toml(doc_ptr, d, errors, base_path='tests/unit/io/toml/fixtures/')

   @assertFalse(errors%has_errors())
   @assertEqual(1, d%swdra)
   @assertEqual(1, d%dramet)
   @assertEqual(-200.0d0, d%basegw, 1.0d-12)
end subroutine
```

- [ ] **Step 4: Run**

```
pixi run -e test test-pfunit
```

Expected: 132 + 1 = 133 tests passing. Existing inline drainage tests still pass.

- [ ] **Step 5: Commit**

```
git add src/io/toml/read_drainage_toml.f90 \
        tests/unit/io/toml/test_read_drainage_toml.pf \
        tests/unit/io/toml/fixtures/drainage_via_file.toml \
        tests/unit/io/toml/fixtures/drainage_external.dra.toml
git commit -m "feat(io/toml): read_drainage_toml follows [drainage].file references"
```

---

### Task 10: Extend `read_crop_toml` + `crop_config_t` + `load_swap_config`

The biggest single task in Phase 4c-a. Touches three files in coordination.

**Files:**
- Modify: `src/config/crop_config.f90` (add per-entry array slots + extended validate)
- Modify: `src/io/toml/read_crop_toml.f90` (per-entry file follow + dispatch)
- Modify: `src/io/toml/load_swap_config.f90` (pass base_path)
- Modify: `tests/unit/io/toml/test_read_crop_toml.pf`

- [ ] **Step 1: Extend `crop_config_t`**

In `src/config/crop_config.f90`:

```fortran
! Add to use clause:
use cropfixed_config_mod, only: cropfixed_config_t
use cropgrass_config_mod, only: cropgrass_config_t

! In type declaration, add after rotation_type:
type(cropfixed_config_t), allocatable :: rotation_fixed(:)
type(cropgrass_config_t), allocatable :: rotation_grass(:)

! Extend crop_config_validate body — after the rotation enum/order checks,
! add a per-entry delegation:
do i = 1, n
   if (self%rotation_type(i) == 1) then
      if (allocated(self%rotation_fixed)) then
         call self%rotation_fixed(i)%validate(errors)
      end if
   else if (self%rotation_type(i) == 3) then
      if (allocated(self%rotation_grass)) then
         call self%rotation_grass(i)%validate(errors)
      end if
   end if
end do
```

(Type 2 entries are skipped — Phase 4c-b adds.)

- [ ] **Step 2: Extend `read_crop_toml`**

```fortran
module read_crop_toml_mod
   use iso_fortran_env, only: real64
   use tomlf, only: toml_table, toml_array, toml_datetime, toml_error, toml_load, get_value, len
   use crop_config_mod, only: crop_config_t
   use cropfixed_config_mod, only: cropfixed_config_t
   use cropgrass_config_mod, only: cropgrass_config_t
   use read_cropfixed_toml_mod, only: read_cropfixed_toml
   use read_cropgrass_toml_mod, only: read_cropgrass_toml
   use toml_field_helpers_mod, only: get_table, get_array_of_tables,   &
                                     get_optional_int_with_default,    &
                                     get_optional_string_with_default, &
                                     parse_date_to_days1900
   use path_helpers_mod, only: resolve_relative_path
   use error_mod, only: error_collection_t, ERR_PARSE_TYPE_MISMATCH, ERR_PARSE_MALFORMED_TOML
   implicit none
   private

   public :: read_crop_toml

contains

   subroutine read_crop_toml(doc, config, errors, base_path)
      type(toml_table), pointer,  intent(in)    :: doc
      type(crop_config_t),        intent(inout) :: config
      type(error_collection_t),   intent(inout) :: errors
      character(len=*), optional, intent(in)    :: base_path

      type(toml_table), pointer             :: sec, item
      type(toml_array), pointer             :: rotation
      type(toml_table), allocatable, target :: crp_doc
      type(toml_table), pointer             :: crp_doc_ptr
      type(toml_error), allocatable         :: terr
      type(toml_datetime)                   :: dtv
      integer :: i, n, stat
      character(len=:), allocatable :: fname, file_abs

      call get_table(doc, 'crop', sec, 'crop', errors)
      if (.not. associated(sec)) return

      call get_optional_int_with_default(sec, 'swcrop', config%swcrop, 0, 'crop.swcrop', errors)

      call get_array_of_tables(sec, 'rotation', rotation, 'crop.rotation', errors)
      if (.not. associated(rotation)) return

      n = len(rotation)
      if (n == 0) return

      allocate(config%rotation_start(n), config%rotation_end(n), &
               config%rotation_file(n),  config%rotation_type(n), &
               config%rotation_fixed(n), config%rotation_grass(n))
      config%rotation_start = 0.0_real64
      config%rotation_end   = 0.0_real64
      config%rotation_file  = ""
      config%rotation_type  = 0

      do i = 1, n
         call get_value(rotation, i, item, stat=stat)
         if (stat /= 0 .or. .not. associated(item)) cycle

         call get_value(item, 'start', dtv, stat=stat)
         if (stat == 0) then
            config%rotation_start(i) = parse_date_to_days1900(dtv)
         else
            call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                               "expected date at crop.rotation.start", &
                               "crop.rotation.start")
         end if

         call get_value(item, 'end', dtv, stat=stat)
         if (stat == 0) then
            config%rotation_end(i) = parse_date_to_days1900(dtv)
         else
            call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                               "expected date at crop.rotation.end", &
                               "crop.rotation.end")
         end if

         call get_optional_string_with_default(item, 'file', fname, '', 'crop.rotation.file', errors)
         config%rotation_file(i) = fname

         call get_optional_int_with_default(item, 'type', config%rotation_type(i), 0, 'crop.rotation.type', errors)

         ! If file= is specified and base_path is available, follow the
         ! reference and dispatch to the matching crop reader.
         if (len_trim(fname) > 0 .and. present(base_path)) then
            file_abs = resolve_relative_path(base_path, trim(fname))
            call toml_load(crp_doc, trim(file_abs), error=terr)
            if (allocated(terr)) then
               call errors%append(ERR_PARSE_MALFORMED_TOML, trim(terr%message), trim(file_abs))
               cycle
            end if
            crp_doc_ptr => crp_doc
            select case (config%rotation_type(i))
            case (1)
               call read_cropfixed_toml(crp_doc_ptr, config%rotation_fixed(i), errors)
            case (3)
               call read_cropgrass_toml(crp_doc_ptr, config%rotation_grass(i), errors)
            ! type 2 (WOFOST) is Phase 4c-b
            end select
         end if
      end do
   end subroutine read_crop_toml

end module read_crop_toml_mod
```

- [ ] **Step 3: Update `load_swap_config`**

```fortran
subroutine load_swap_config(path, config, errors)
   use path_helpers_mod, only: directory_of
   character(len=*),          intent(in)    :: path
   type(swap_config_t),       intent(inout) :: config
   type(error_collection_t),  intent(inout) :: errors

   type(toml_table), allocatable, target :: doc
   type(toml_table), pointer             :: doc_ptr
   type(toml_error), allocatable         :: terr
   character(len=:), allocatable         :: base_dir

   call toml_load(doc, trim(path), error=terr)
   if (allocated(terr)) then
      call errors%append(ERR_PARSE_MALFORMED_TOML, trim(terr%message), trim(path))
      return
   end if

   doc_ptr => doc
   base_dir = directory_of(trim(path))

   call read_general_toml    (doc_ptr, config%general,    errors)
   call read_simulation_toml (doc_ptr, config%simulation, errors)
   call read_meteorology_toml(doc_ptr, config%meteo,      errors)
   call read_drainage_toml   (doc_ptr, config%drain,      errors, base_path=base_dir)
   call read_soil_toml       (doc_ptr, config%soil,       errors)
   call read_crop_toml       (doc_ptr, config%crop,       errors, base_path=base_dir)
end subroutine load_swap_config
```

Add `use path_helpers_mod, only: directory_of` at the module top.

- [ ] **Step 4: Run regression**

```
pixi run -e test test-pfunit
```

Expected: 133 tests passing. The existing crop and load tests should still pass because base_path is optional and inline path is preserved.

```
pixi run -e test check-fast
```

Expected: 4/4 regression cases green.

- [ ] **Step 5: Commit**

```
git add src/config/crop_config.f90 \
        src/io/toml/read_crop_toml.f90 \
        src/io/toml/load_swap_config.f90
git commit -m "feat(io/toml): cross-file crop loading + rotation_fixed/grass slots"
```

---

## Stage F — Per-case TOML authoring (Tasks 11-14)

For each case, ONE task that authors `swap.dra.toml` + the case's `*.crp.toml` files inside the `tests/swap-cases` submodule, extends `swap.toml` to use the new file references, and bumps the outer-repo submodule pointer. After each task, run a smoke test confirming `load_swap_config` + validate + finalize produce no fatal errors.

### Task 11: hupselbrook (case 1) TOML splits

- [ ] **Step 1: Inside the submodule**, author `tests/swap-cases/toml/1.hupselbrook/`:
  - `swap.dra.toml` extracted from `tests/swap-cases/1.hupselbrook/swap.dra` — basic-section keys (dramet, swdivd, basegw, entres, shape, etc.). For dramet=2, NO `[[drainage.levels]]` array (single-level convention).
  - `maizes.crp.toml` extracted from `tests/swap-cases/1.hupselbrook/maizes.crp` — the type-1 fields covered by `cropfixed_config_t` so far (Feddes, salinity, root, interception, idev, lcc).
  - `grassd.crp.toml` extracted from `tests/swap-cases/1.hupselbrook/grassd.crp` — the type-3 fields covered by `cropgrass_config_t` so far.
  - `potatod.crp.toml` PLACEHOLDER — single-line file: `# Placeholder for Phase 4c-b WOFOST crop type 2.`

- [ ] **Step 2: Update `swap.toml`** in the submodule. Change `[drainage]` to:
```toml
[drainage]
file = "swap.dra.toml"
```
And update each `[[crop.rotation]]` entry to add `file = "<crop>.crp.toml"`:
```toml
[[crop.rotation]]
start = 2002-05-01
end   = 2002-09-30
file  = "maizes.crp.toml"
type  = 1
```
(Existing rotation entries — preserve start/end dates, just add file= and adjust type per actual `.swp.template`.)

- [ ] **Step 3: Commit submodule + bump outer pointer**

Inside `tests/swap-cases`:
```
git checkout main
git add toml/1.hupselbrook/
git commit -m "toml(hupselbrook): split into swap.dra.toml + .crp.toml files"
```
Then in outer repo:
```
git add tests/swap-cases
git commit -m "test(io/toml): hupselbrook cross-file TOML structure"
```

- [ ] **Step 4: Smoke test**

```
pixi run -e test test-pfunit
```

The existing `test_hupselbrook_loads.pf` smoke tests should still pass (load + validate + finalize produces no fatal errors). If any fail, debug — typically a missing field or a type mismatch in the new TOMLs. Iterate inside the submodule with focused commits.

---

### Tasks 12, 13, 14: grassgrowth, oxygenstress, surfacewater TOML splits

Same pattern as Task 11.

- [ ] **Task 12: grassgrowth (case 2)** — author `swap.dra.toml` + `grassd.crp.toml`. Update `swap.toml`. Smoke test. Commit submodule + outer.
- [ ] **Task 13: oxygenstress (case 4)** — author `swap.dra.toml` + `grassd.crp.toml` (per-case copy; values may differ from grassgrowth). Update `swap.toml`. Smoke test. Commit.
- [ ] **Task 14: surfacewater (case 6)** — author `swap.dra.toml` + `grass.crp.toml` (type 1). Update `swap.toml`. Smoke test. Commit.

For each: the existing per-case smoke test in `tests/unit/io/toml/test_all_cases_smoke.pf` should keep passing.

---

## Stage G — Parity tests (Tasks 15-18)

Each task authors a per-case parity test (or extends hupselbrook's) following the Phase 4b pattern. Each test:

1. chdir into the case's legacy directory.
2. Stage `swap_linux.swp.template` → `swap.swp`.
3. Call `readswap()` (no args).
4. chdir back.
5. Call `load_swap_config('tests/swap-cases/toml/<case>/swap.toml', config, errors)`.
6. Call `config%validate(errors)` and `config%finalize(errors)`.
7. Assert `variables%foo == config%<section>%<field>` for ~30-50 fields, including the LAST rotation entry's crop data via `config%crop%rotation_fixed(N)` or `rotation_grass(N)`.

### Task 15: Update `test_hupselbrook_parity.pf`

- [ ] **Step 1: Tighten existing assertions** (post-Task 1, 2, 3 effects)

In `tests/unit/io/toml/test_hupselbrook_parity.pf`:
- Date assertions (`tstart`, `tend`): tolerance `1.5d0` → `1.0d-9`.
- Re-enable `swmacro` assertion (was dropped in Phase 4b).
- Re-enable `nrlevs` assertion (was dropped in Phase 4b).

- [ ] **Step 2: Add type-1 (maizes) sub-section assertions**

After `readswap` returns, the LAST rotation entry's crop globals are in `variables`. Hupselbrook's last entry is grass (type 3) — so the type-1 (maizes) parity is harder to test directly here.

For Phase 4c-a, assert ONLY on the last entry (which is type 3). The maizes (type 1) entry can be smoke-tested via the load_swap_config path (its `rotation_fixed(I)` slot is populated; assert non-default field values).

```fortran
@test
subroutine test_hupselbrook_last_rotation_grass_parity()
   use funit
   use variables
   use swap_config_mod
   use load_swap_config_mod
   use chdir_helper_mod
   ! ... (chdir + readswap + load_swap_config boilerplate as in existing tests)
   ! After readswap, variables%hlim1, variables%hlim2u, variables%hlim2l, etc.
   ! hold the LAST crop's globals (grassd, type 3 in hupselbrook).

   ! Find the last entry's index that's type 3
   integer :: i_last, n
   n = size(config%crop%rotation_type)
   i_last = 0
   do i = n, 1, -1
      if (config%crop%rotation_type(i) == 3) then
         i_last = i
         exit
      end if
   end do
   @assertTrue(i_last > 0)

   @assertEqual(hlim1,  config%crop%rotation_grass(i_last)%hlim1,  1.0d-9)
   @assertEqual(hlim2u, config%crop%rotation_grass(i_last)%hlim2u, 1.0d-9)
   @assertEqual(hlim2l, config%crop%rotation_grass(i_last)%hlim2l, 1.0d-9)
   ! ... add more fields as available
end subroutine
```

NOTE: the actual variable names in `variables` may differ from `cropgrass_config_t` field names. Verify each via grep before adding the assertion. Drop assertions that hit fields not in `variables` (those are parity gaps for later phases).

- [ ] **Step 3: Run, iterate, commit**

```
pixi run -e test test-pfunit
```

Expected: existing 6 hupselbrook parity tests still pass + 1 new test for grass parity = 7 hupselbrook tests. Plus the previously dropped `swmacro` and `nrlevs` assertions resume passing.

```
git commit -m "test(io/toml): tighten hupselbrook parity, add grass entry assertions"
```

---

### Task 16: `test_grassgrowth_parity.pf`

Full parity for case 2.grassgrowth (only crop type 3).

- [ ] **Step 1: Author the test**

```fortran
@test
subroutine test_grassgrowth_general_parity()
   use funit
   use variables
   use swap_config_mod
   use load_swap_config_mod
   use chdir_helper_mod
   use error_mod
   type(swap_config_t)      :: config
   type(error_collection_t) :: errors
   character(len=1024)      :: orig_cwd

   call get_cwd(orig_cwd)
   call chdir_to('tests/swap-cases/2.grassgrowth')
   call stage_swp_template('swap_linux.swp.template', 'swap')
   call readswap()
   call chdir_to(trim(orig_cwd))

   call load_swap_config('tests/swap-cases/toml/2.grassgrowth/swap.toml', config, errors)
   call config%validate(errors)
   call config%finalize(errors)
   @assertFalse(errors%has_fatals())

   @assertEqual(trim(project), trim(config%general%project))
   @assertEqual(swscre,        config%general%swscre)
end subroutine

@test
subroutine test_grassgrowth_simulation_parity()
   ! ... same setup ...
   @assertEqual(tstart,    config%simulation%tstart, 1.0d-9)
   @assertEqual(tend,      config%simulation%tend,   1.0d-9)
   @assertEqual(nprintday, config%simulation%nprintday)
end subroutine

@test
subroutine test_grassgrowth_meteorology_parity()
   ! ... 6-9 meteo assertions like hupselbrook's
end subroutine

@test
subroutine test_grassgrowth_drainage_parity()
   ! ... dramet, swdivd, swdislay, nrlevs (now passes with finalize clobber)
end subroutine

@test
subroutine test_grassgrowth_soil_parity()
   ! ... swsophy, swhyst, swinco, gwli
end subroutine

@test
subroutine test_grassgrowth_crop_rotation_parity()
   ! ... swcrop, rotation entry counts and types
end subroutine

@test
subroutine test_grassgrowth_grass_data_parity()
   ! ... LAST rotation entry's grass-specific fields
   @assertEqual(hlim1,  config%crop%rotation_grass(N)%hlim1,  1.0d-9)
   ! ... etc
end subroutine
```

- [ ] **Step 2: Wire**

`tests/unit/testSuites.inc`: append `ADD_TEST_SUITE(test_grassgrowth_parity_suite)`.
`tests/unit/meson.build` `pf_files`: append `'io/toml/test_grassgrowth_parity.pf'`.

- [ ] **Step 3: Run, iterate**

```
pixi run -e test test-pfunit
```

Iterate as in Phase 4b: each failing assertion → check if it's a TOML mismatch (fix the `.toml`), config_t gap (skip + log Phase 4d work), or genuine parity bug (fix appropriately).

- [ ] **Step 4: Commit**

```
git add tests/unit/io/toml/test_grassgrowth_parity.pf \
        tests/unit/testSuites.inc tests/unit/meson.build
git commit -m "test(io/toml): grassgrowth (case 2) full parity"
```

---

### Tasks 17, 18: oxygenstress and surfacewater parity

- [ ] **Task 17: `test_oxygenstress_parity.pf`** — same shape as grassgrowth (also type 3). Wire + run + iterate + commit.
- [ ] **Task 18: `test_surfacewater_parity.pf`** — same shape but type 1 instead of type 3. Use `rotation_fixed(N)` instead of `rotation_grass(N)`. Wire + run + iterate + commit.

---

## Stage H — Closeout (Tasks 19-21)

### Task 19: Refresh `docs/configuration-schema.md`

Reconcile the doc to Phase 4b/4c-a reality.

- [ ] **Step 1: Open the file** and walk through it. Three known divergences to fix:
  1. `[crop.rotation]` is documented as parallel string arrays. Reality: `[[crop.rotation]]` array-of-tables. Update.
  2. `[output]` placed at top level. Reality: nested under `[simulation.output]`. Update.
  3. `[meteorology].alt` and `[meteorology].altw` documented under `[meteorology.evapotranspiration]`. Reality: at section root. Update.

- [ ] **Step 2: Add a new "cross-file references" subsection** documenting `[drainage].file` and `[[crop.rotation]].file` patterns and the path-resolution rule (relative to the referencing file's directory).

- [ ] **Step 3: Add subsection sketches for `cropfixed` and `cropgrass`**:
```markdown
### `*.crp.toml` — type 1 (fixed crop)

Sections: [phenology], [light], [root], [water_stress], [salinity], [interception].
See src/config/cropfixed_config.f90 for the authoritative field set.

### `*.crp.toml` — type 3 (grass)

Sections: as type 1 plus [mowing], [grazing].
See src/config/cropgrass_config.f90.
```

- [ ] **Step 4: Cross-link** to `docs/toml-format-guide.md` for conventions.

- [ ] **Step 5: Build FORD**

```
pixi run -e docs docs-build
```

Expected: clean build.

- [ ] **Step 6: Commit**

```
git add docs/configuration-schema.md
git commit -m "docs(schema): reconcile to Phase 4b/4c-a reality + cross-file"
```

---

### Task 20: Coverage rebaseline

```
pixi run -e coverage coverage-report
```

Update `docs/coverage-baseline.md` with Phase 4c-a numbers. Append a new section after Phase 4b's; note the new modules (`cropfixed`, `cropgrass`, `path_helpers`) and their coverage levels.

```
git add docs/coverage-baseline.md
git commit -m "docs(baselines): record Phase 4c-a coverage numbers"
```

Append a check-full record to `tests/regression/baselines/`:

```
git add tests/regression/baselines/
git commit -m "docs(baselines): record Phase 4c-a check-full output"
```

---

### Task 21: Phase 4c-a tag and fast-forward

- [ ] **Step 1: Final check-full**

```
pixi run -e test check-full
```

Expected: 6/6 cases green.

- [ ] **Step 2: Verify "no physics changed" constraint**

```
git diff rescue/phase-4b-parity..HEAD -- \
   src/atmosphere/ src/soil/ src/crop/ src/boundary/ \
   src/drainage/ src/macropore/ src/solute/ src/heat/ src/utils/
```

Expected: small diff. Only `src/io/readswap.f90` (the swmacro shadow removal) and `src/config/drainage_config.f90` (finalize update) changes intersect with these paths — `src/io/` is in scope, `src/config/` is in scope. No physics algorithm changes.

- [ ] **Step 3: Fast-forward main**

```
git checkout main
git merge --ff-only development
git checkout development
```

- [ ] **Step 4: Tag**

```
git tag rescue/phase-4c-a-crop-fixed-grass
```

No push.

- [ ] **Step 5: Verify**

```
git log --oneline -n 1 main
git tag --list 'rescue/phase-*'
```

Expected: 7 rescue tags, all on consistent commits.

---

## Self-Review

### Spec coverage

Walking the spec sections and pointing each to a task:

- **In scope: Cross-file TOML loading** → Tasks 9, 10
- **In scope: cropfixed_config_t** → Task 5
- **In scope: cropgrass_config_t** → Task 6
- **In scope: read_cropfixed_toml + read_cropgrass_toml** → Tasks 7, 8
- **In scope: path_helpers** → Task 4
- **In scope: crop_config_t extension** → Task 10 step 1
- **In scope: read_drainage_toml + read_crop_toml + load_swap_config extensions** → Tasks 9, 10
- **In scope: parse_date_to_days1900 alignment** → Task 1
- **In scope: swmacro shadow removal** → Task 2
- **In scope: nrlevs finalize mirror** → Task 3
- **In scope: per-case TOML files for cases 1, 2, 4, 6** → Tasks 11, 12, 13, 14
- **In scope: per-case parity tests** → Tasks 15 (extends), 16, 17, 18
- **In scope: schema doc refresh** → Task 19
- **In scope: coverage rebaseline + tag** → Tasks 20, 21

All spec scope items are covered. No gaps found.

### Placeholder scan

- "// add more fields as available" at task 15 step 2 — this is contextual guidance for an iterative test development pattern, not a placeholder. Acceptable.
- "// existing per-level allocation + walk — preserve from prior implementation" at task 9 step 1 — the verbatim block IS provided right below; the comment is just an explanation. Acceptable.

No actual placeholders.

### Type consistency

- `cropfixed_config_t` and `cropgrass_config_t` field names verified consistent across Tasks 5, 6, 7, 8, 10.
- `rotation_fixed(:)` and `rotation_grass(:)` array names consistent in Tasks 10, 15, 16, 17, 18.
- `read_cropfixed_toml` / `read_cropgrass_toml` signatures consistent in Tasks 7, 8, 10.
- `directory_of` and `resolve_relative_path` used consistently in Tasks 4, 9, 10.
- `base_path` argument name consistent across reader extensions.
