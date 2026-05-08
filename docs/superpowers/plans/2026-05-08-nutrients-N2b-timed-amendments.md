# [nutrients] N2b — Timed amendments via CSV companion Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an optional `events_file` field on `nutrients_config_t` pointing to a CSV companion of timed soil management events (date, material, amount_kgha, volat_fraction). Extend `apply_nutrients` with CSV staging + per-row validation + sort-by-date + same-day grouping logic. Runtime path (`SoilManagement(3)` → `Wofost_SoilAmendents`) unchanged.

**Architecture:** Schema gains one allocatable string field. TOML reader parses the optional `events_file` path. The adapter's existing `apply_nutrients` (added in N2a, commit `5fb0e0b`) gets a new private helper `apply_nutrients_events` that — when `events_file` is populated — loads the CSV via `read_csv_table` (auto-converts column-1 ISO dates to days-since-1900), validates each row's material/amount/volat ranges, sorts in-place by date, populates `MatNum(:)`, `Amend(:)` (kg/ha → kg/m²), `VolaFrac(:)`, then groups same-day dosages into `TimeAmend(:)`/`NuAmend(:)`/`iAmend(:,:)`/`namend`/`isme`. Sort + group algorithm is a verbatim port of the deleted pre-collapse `SoilManagement(1)` legacy block.

**Tech Stack:** Fortran 2008, meson + ninja, gfortran, pFUnit. Reference patterns: `apply_ssdi_mode0` (commit `328db06`) for the CSV-staging shape; the existing `apply_nutrients` (commit `5fb0e0b`) for the procedure to extend; `wofost_soil_declarations.f90` for the legacy globals (`maxamn=1000`, `MatNum(maxamn)`, `Amend(maxamn)`, `VolaFrac(maxamn)`, `TimeAmend(maxamn)`, `NuAmend(maxamn)`, `iamend(maxamn,maxamn)`, `namend`, `isme`).

**Spec:** `docs/superpowers/specs/2026-05-08-nutrients-N2b-timed-amendments-design.md`.

---

## File Structure

**Created:**

| Path | Responsibility |
|---|---|
| `tests/unit/io/toml/fixtures/nutrients_events_small.csv` | 3-row CSV fixture for adapter tests; deliberately includes one out-of-order row so the sort test is meaningful + one same-day pair so the grouping test is meaningful. |
| `docs/adr/0027-nutrients-N2b-timed-amendments.md` | New ADR. |

**Modified:**

| Path | Change |
|---|---|
| `src/config/nutrients_config.f90` | Add `character(len=:), allocatable :: events_file` field on `nutrients_config_t`. No new validator rules in this sub-arc (per-row CSV validation lives in the adapter). |
| `src/io/toml/read_nutrients_toml.f90` | Parse optional `events_file` string from `[nutrients]`. |
| `src/io/toml/config_to_variables.f90` | Extend `apply_nutrients` body with `call apply_nutrients_events(cfg)` at the tail. Add the new private `apply_nutrients_events(cfg)` procedure. |
| `tests/unit/io/toml/test_nutrients_config.pf` | +6 tests (parser × 2, adapter × 4 — see Task 4). |
| `docs/adr/index.md` | Append ADR 0027 row. |
| `docs/configuration-schema.md` | Document `events_file` and CSV column format under `[nutrients]`. |

**Pre-flight inventory results (verified at spec-authoring):**

- `nutrients_config_t` is at `src/config/nutrients_config.f90` line 29; `events_file` field appends after `sorp_coef` (line 32).
- `wofost_soil_declarations.f90` exports `maxamn` (parameter, value 1000) at line 11; the legacy globals (`MatNum`, `Amend`, etc.) at lines 61-73.
- `pathwork` is set early in `config_to_variables` from `[general.paths].work` (existing mechanism).
- Existing `apply_nutrients` is in `src/io/toml/config_to_variables.f90`; check the tail for the right insertion point. Find by `grep -n "subroutine apply_nutrients" src/io/toml/config_to_variables.f90`.
- The CSV fixture mirrors the `tests/unit/io/toml/fixtures/ssdi_events_small.csv` shape (commit Task-4 of SSDI port). Same dir, same naming convention.

---

## Task 1: Schema field + reader

**Files:**
- Modify: `src/config/nutrients_config.f90` (add `events_file` field)
- Modify: `src/io/toml/read_nutrients_toml.f90` (parse `events_file`)
- Modify: `tests/unit/io/toml/test_nutrients_config.pf` (append 2 parser tests)

- [ ] **Step 1: Append the parser tests**

Append to `tests/unit/io/toml/test_nutrients_config.pf`:

```fortran
@test
subroutine test_read_nutrients_with_events_file()
   use nutrients_config_mod, only: nutrients_config_t
   use read_nutrients_toml_mod, only: read_nutrients_toml
   use error_mod, only: error_collection_t
   use tomlf, only: toml_table, toml_loads
   use funit
   implicit none
   type(toml_table), allocatable, target :: doc
   type(toml_table), pointer :: doc_ptr
   type(nutrients_config_t)  :: cfg
   type(error_collection_t)  :: errs
   character(len=*), parameter :: toml_text =                       &
      '[nutrients]'                              // new_line('a') // &
      'events_file = "amendments.csv"'

   call toml_loads(doc, toml_text)
   doc_ptr => doc

   call read_nutrients_toml(doc_ptr, cfg, errs)
   call assertEqual(0, errs%count(), 'expected zero parse errors')
   @assertTrue(cfg%present, 'expected present=true')
   @assertTrue(allocated(cfg%events_file), 'events_file allocated')
   @assertEqual('amendments.csv', trim(cfg%events_file))
end subroutine

@test
subroutine test_read_nutrients_no_events_file()
   use nutrients_config_mod, only: nutrients_config_t
   use read_nutrients_toml_mod, only: read_nutrients_toml
   use error_mod, only: error_collection_t
   use tomlf, only: toml_table, toml_loads
   use funit
   implicit none
   type(toml_table), allocatable, target :: doc
   type(toml_table), pointer :: doc_ptr
   type(nutrients_config_t)  :: cfg
   type(error_collection_t)  :: errs
   character(len=*), parameter :: toml_text =                       &
      '[nutrients]'                              // new_line('a') // &
      'sorp_coef = 0.005'

   call toml_loads(doc, toml_text)
   doc_ptr => doc

   call read_nutrients_toml(doc_ptr, cfg, errs)
   call assertEqual(0, errs%count())
   @assertFalse(allocated(cfg%events_file), 'events_file unallocated when absent')
end subroutine
```

- [ ] **Step 2: Run to verify failures**

```
pixi run -e test build-linux
```
Expected: build fails — `cfg%events_file` is not a member of `nutrients_config_t`.

- [ ] **Step 3: Add the field to `nutrients_config_t`**

In `src/config/nutrients_config.f90`, find the `nutrients_config_t` type (around line 29) and add the new field:

```fortran
   type :: nutrients_config_t
      logical      :: present   = .false.
      real(real64) :: sorp_coef = 0.0_real64
      character(len=:), allocatable :: events_file       ! relative to pathwork; CSV companion (N2b, ADR 0027)
      type(nutrients_initial_t) :: initial
   contains
      procedure :: validate => nutrients_config_validate
      procedure :: finalize => nutrients_config_finalize
   end type nutrients_config_t
```

The validator stays as-is. Per-row CSV validation lives in the adapter (Task 3).

- [ ] **Step 4: Wire the field into the reader**

In `src/io/toml/read_nutrients_toml.f90`, add to the imports:
```fortran
   use toml_field_helpers_mod, only: get_optional_real_with_default, &
                                     get_optional_string_with_default
```
(The first import already exists; add the second.)

Inside `read_nutrients_toml`, just after the `sorp_coef` parse (and before the `init_tbl` lookup), add:

```fortran
      call get_optional_string_with_default(nut_tbl, 'events_file', &
                                            cfg%events_file, '', &
                                            'nutrients.events_file', errors)
```

If `get_optional_string_with_default` doesn't allocate the string when the key is absent (i.e., leaves `cfg%events_file` unallocated), that's the desired behaviour. Verify by inspecting the helper's signature in `src/io/toml/toml_field_helpers.f90`. If the helper unconditionally allocates to the default value (empty string), the test `test_read_nutrients_no_events_file` will fail with `events_file allocated to empty string` — adjust the assertion to `len_trim(cfg%events_file) == 0` instead of `.not. allocated`.

- [ ] **Step 5: Build + run tests**

```
pixi run -e test build-linux
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:|tests,"
pixi run -e test check-full 2>&1 | grep "Results:"
```
Expected: build clean; pFUnit `Ok: 1, Fail: 0` (test count grows by 2); check-full `5 passed, 0 failed`.

- [ ] **Step 6: Commit**

```bash
git add src/config/nutrients_config.f90 src/io/toml/read_nutrients_toml.f90 tests/unit/io/toml/test_nutrients_config.pf
git commit -m "$(cat <<'EOF'
feat(config): add events_file field to nutrients_config_t

Adds an optional `events_file` (allocatable string) field on
nutrients_config_t. The TOML reader parses it from
[nutrients].events_file; absent or empty → unallocated.

The field is the CSV companion path for timed soil management
events (fertilizer applications). The adapter staging + sort +
group logic lands in the next commit (Task 3).

Per-row CSV validation lives in the adapter, not in
nutrients_config_validate. Mirrors the SSDI mode-0 pattern
(commit 328db06).

Verified: build clean; pFUnit Ok: 1, Fail: 0 (test count +2);
check-full 5/5.

Part of [nutrients] N2b (spec
docs/superpowers/specs/2026-05-08-nutrients-N2b-timed-amendments-design.md).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: CSV fixture

**Files:**
- Create: `tests/unit/io/toml/fixtures/nutrients_events_small.csv`

- [ ] **Step 1: Create the fixture file**

Create `tests/unit/io/toml/fixtures/nutrients_events_small.csv`:

```
date,material,amount_kgha,volat_fraction
2003-06-20,1,25000.0,0.10
2003-04-15,10,100.0,0.05
2003-06-20,3,500.0,0.02
```

The first two rows are deliberately out-of-order (June before April) so the sort test in Task 3 is meaningful. The two rows dated `2003-06-20` (materials 1 and 3) become the same-day group so the grouping test is meaningful.

After sort: row 1 = `2003-04-15` (mat 10), rows 2-3 = `2003-06-20` (mat 1, then mat 3 — same-day group).

- [ ] **Step 2: Commit**

```bash
git add tests/unit/io/toml/fixtures/nutrients_events_small.csv
git commit -m "$(cat <<'EOF'
test(fixtures): add nutrients_events_small.csv for N2b adapter tests

Three-row fixture with deliberate ordering to exercise both the
sort-by-date and same-day grouping logic in the upcoming
apply_nutrients_events adapter:

- Rows 1-2 are out-of-order (June before April) — sort test.
- Rows 1, 3 are same-date (2003-06-20) — grouping test.

Used by Task 3 / 4 of [nutrients] N2b (spec
docs/superpowers/specs/2026-05-08-nutrients-N2b-timed-amendments-design.md).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: `apply_nutrients_events` adapter

**Files:**
- Modify: `src/io/toml/config_to_variables.f90` (extend `apply_nutrients` + add private `apply_nutrients_events`)

- [ ] **Step 1: Read the current `apply_nutrients` boundaries**

```bash
grep -n "subroutine apply_nutrients" src/io/toml/config_to_variables.f90
```

You'll see two matches: `subroutine apply_nutrients` (the public entry, added in N2a — commit `5fb0e0b`) and the matching `end subroutine apply_nutrients`.

- [ ] **Step 2: Extend `apply_nutrients` to call the new helper**

Inside `apply_nutrients(cfg)` (in `src/io/toml/config_to_variables.f90`), at the very end of the procedure body (before `end subroutine apply_nutrients`), add:

```fortran
      ! N2b (ADR 0027): stage timed amendments from the CSV companion.
      call apply_nutrients_events(cfg)
```

- [ ] **Step 3: Add the `apply_nutrients_events` procedure**

Add the new private procedure inside the `contains` block of `config_to_variables_mod`, immediately after `apply_nutrients`:

```fortran
   !> Stage timed soil management events from the CSV companion
   !! at cfg%events_file. Sets the legacy globals consumed by
   !! SoilManagement(3) and Wofost_SoilAmendents.
   !!
   !! Default (no events_file or empty CSV): namend = 0, isme = 1.
   !!
   !! See ADR 0027 ([nutrients] N2b).
   subroutine apply_nutrients_events(cfg)
      use, intrinsic :: iso_fortran_env, only: real64
      use csv_reader_mod, only: read_csv_table
      use error_mod, only: error_collection_t, fatalerr_collected
      use nutrients_config_mod, only: nutrients_config_t
      use variables, only: pathwork
      use wofost_soil_declarations, only: MatNum, Amend, VolaFrac, &
                                            TimeAmend, NuAmend, iamend, &
                                            namend, isme, maxamn
      type(nutrients_config_t), intent(in) :: cfg

      real(real64), allocatable :: tbl(:,:)
      type(error_collection_t)  :: errs
      character(len=300) :: csvpath
      character(len=14)  :: hdr(4)
      integer :: i, j, n
      real(real64) :: tmp_date, tmp_amount, tmp_volat
      real(real64) :: tmp_mat_real

      ! Default: no amendments. Reset legacy globals to a known state.
      namend = 0
      isme   = 1

      if (.not. allocated(cfg%events_file)) return
      if (len_trim(cfg%events_file) == 0)   return

      hdr(1) = 'date          '
      hdr(2) = 'material      '
      hdr(3) = 'amount_kgha   '
      hdr(4) = 'volat_fraction'
      csvpath = trim(pathwork) // trim(cfg%events_file)
      call read_csv_table(trim(csvpath), hdr, tbl, errs)
      call errs%abort_if_fatal()

      n = 0
      if (allocated(tbl)) n = size(tbl, 1)
      if (n < 1) return     ! Empty CSV: no amendments. Not an error.
      if (n > maxamn) then
         call fatalerr_collected('apply_nutrients_events', &
            'CSV row count exceeds maxamn (1000)')
         return
      end if

      ! Per-row validation
      do i = 1, n
         if (nint(tbl(i, 2)) < 1 .or. nint(tbl(i, 2)) > 20) then
            call fatalerr_collected('apply_nutrients_events', &
               'material out of range [1, 20]')
            return
         end if
         if (tbl(i, 3) < 0.0_real64 .or. tbl(i, 3) > 500000.0_real64) then
            call fatalerr_collected('apply_nutrients_events', &
               'amount_kgha out of range [0, 500000]')
            return
         end if
         if (tbl(i, 4) < 0.0_real64 .or. tbl(i, 4) > 1.0_real64) then
            call fatalerr_collected('apply_nutrients_events', &
               'volat_fraction out of range [0, 1]')
            return
         end if
      end do

      ! Sort by date (in-place bubble sort, mirrors deleted SoilManagement(1)).
      ! Acceptable O(n^2) given n <= 1000 and this runs once at config-load.
      do i = 1, n - 1
         do j = i + 1, n
            if (tbl(i, 1) > tbl(j, 1)) then
               tmp_date     = tbl(i, 1); tbl(i, 1) = tbl(j, 1); tbl(j, 1) = tmp_date
               tmp_mat_real = tbl(i, 2); tbl(i, 2) = tbl(j, 2); tbl(j, 2) = tmp_mat_real
               tmp_amount   = tbl(i, 3); tbl(i, 3) = tbl(j, 3); tbl(j, 3) = tmp_amount
               tmp_volat    = tbl(i, 4); tbl(i, 4) = tbl(j, 4); tbl(j, 4) = tmp_volat
            end if
         end do
      end do

      ! Populate per-event legacy globals
      do i = 1, n
         MatNum(i)   = nint(tbl(i, 2))
         Amend(i)    = 1.0e-4_real64 * tbl(i, 3)   ! kg/ha -> kg/m^2
         VolaFrac(i) = tbl(i, 4)
      end do

      ! Group dosages per date (mirrors deleted SoilManagement(1)).
      j = 1
      NuAmend(j)   = 1
      TimeAmend(j) = tbl(1, 1)
      iamend(1, 1) = 1
      do i = 2, n
         if (abs(tbl(i, 1) - tbl(i - 1, 1)) < 1.0e-3_real64) then
            NuAmend(j) = NuAmend(j) + 1
         else
            j = j + 1
            NuAmend(j)   = 1
            TimeAmend(j) = tbl(i, 1)
         end if
         iamend(j, NuAmend(j)) = i
      end do

      namend = j
      isme   = 1
   end subroutine apply_nutrients_events
```

Notes on global names: in `wofost_soil_declarations.f90` the array is named lowercase `iamend(maxamn,maxamn)` (per line 61). Use the actual case-style of that declaration; if `iAmend` (mixed-case) is what's exported, swap accordingly. Verify by `grep -nE "^\s*integer.*iamend" src/crop/wofost_soil_declarations.f90` before pasting.

If `apply_nutrients_events` declared as `private` causes pFUnit imports in Task 4 to fail, add it to the module's `public ::` list (the test imports it via `use config_to_variables_mod, only: apply_nutrients_events`).

- [ ] **Step 4: Build and verify (no new tests yet — that's Task 4)**

```
pixi run -e test build-linux
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:|tests,"
pixi run -e test check-full 2>&1 | grep "Results:"
```
Expected: build clean; pFUnit unchanged from Task 1 (still +2 from baseline; the new procedure has no callers in tests yet); check-full `5 passed, 0 failed` (regression cases all leave `events_file` unallocated → adapter takes the no-op branch and resets `namend=0, isme=1`, which is the same as their pre-arc state).

- [ ] **Step 5: Commit**

```bash
git add src/io/toml/config_to_variables.f90
git commit -m "$(cat <<'EOF'
feat(io): apply_nutrients_events — CSV staging + sort + group

Adds the apply_nutrients_events private procedure to
config_to_variables_mod, called from the existing apply_nutrients
adapter (added in N2a, commit 5fb0e0b).

When cfg%events_file is allocated and non-empty:
- read_csv_table loads (date, material, amount_kgha,
  volat_fraction); column 1 auto-converts ISO dates to
  days-since-1900 per ADR 0012.
- Per-row validation (material in [1, 20], amount in
  [0, 500000], volat in [0, 1]; tightening from legacy's
  typo'd [0, 500000] for volat).
- In-place bubble sort by date (mirrors deleted
  pre-collapse SoilManagement(1) algorithm verbatim).
- Populates per-event MatNum, Amend (kg/ha -> kg/m^2),
  VolaFrac legacy globals.
- Groups same-day dosages into TimeAmend, NuAmend, iamend;
  sets namend, isme = 1.

When events_file is unallocated/empty: namend = 0, isme = 1
(no amendments — natural-mineralization-only baseline).

Adapter unit tests follow in next commit. check-full
byte-identical (regression cases all leave events_file
unallocated; their reset to namend=0/isme=1 matches their
pre-arc no-amendment-init state).

Part of [nutrients] N2b (spec
docs/superpowers/specs/2026-05-08-nutrients-N2b-timed-amendments-design.md).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Adapter pFUnit tests

**Files:**
- Modify: `tests/unit/io/toml/test_nutrients_config.pf` (append 4 adapter tests)

- [ ] **Step 1: Append the 4 adapter tests**

Append to `tests/unit/io/toml/test_nutrients_config.pf`:

```fortran
@test
subroutine test_apply_events_csv_populates_legacy_globals()
   use, intrinsic :: iso_fortran_env, only: real64
   use nutrients_config_mod, only: nutrients_config_t
   use config_to_variables_mod, only: apply_nutrients_events
   use variables, only: pathwork
   use wofost_soil_declarations, only: MatNum, Amend, VolaFrac, &
                                         TimeAmend, NuAmend, iamend, &
                                         namend, isme
   use funit
   implicit none
   type(nutrients_config_t) :: cfg
   integer :: dy_apr15, dy_jun20

   pathwork = 'tests/unit/io/toml/fixtures/'
   cfg%events_file = 'nutrients_events_small.csv'

   call apply_nutrients_events(cfg)

   ! 3 rows, two share the date 2003-06-20 → 2 amendment groups.
   call assertEqual(2, namend, 'two distinct amendment dates')
   call assertEqual(1, isme,   'isme reset to 1')

   ! After sort: 2003-04-15 first (mat 10), 2003-06-20 second group (mat 1, mat 3).
   ! NuAmend(1) = 1 (just the 2003-04-15 event)
   ! NuAmend(2) = 2 (the two 2003-06-20 events)
   call assertEqual(1, NuAmend(1))
   call assertEqual(2, NuAmend(2))

   ! Sorted order: row 1 = mat 10, row 2 = mat 1, row 3 = mat 3
   call assertEqual(10, MatNum(1))
   call assertEqual(1,  MatNum(2))
   call assertEqual(3,  MatNum(3))

   ! Amend = kg/ha * 1.0e-4 = kg/m^2
   @assertEqual(0.01_real64,    Amend(1), tolerance=1.0e-12_real64)  ! 100   * 1e-4
   @assertEqual(2.5_real64,     Amend(2), tolerance=1.0e-12_real64)  ! 25000 * 1e-4
   @assertEqual(0.05_real64,    Amend(3), tolerance=1.0e-12_real64)  ! 500   * 1e-4

   ! VolaFrac round-trips exactly
   @assertEqual(0.05_real64, VolaFrac(1), tolerance=1.0e-12_real64)
   @assertEqual(0.10_real64, VolaFrac(2), tolerance=1.0e-12_real64)
   @assertEqual(0.02_real64, VolaFrac(3), tolerance=1.0e-12_real64)

   ! Days-since-1900 for 2003-04-15 and 2003-06-20:
   !   2003-04-15 = 37726 (per the SSDI test fixture's same date)
   !   2003-06-20 = 37792 (66 days later)
   dy_apr15 = 37726
   dy_jun20 = 37792
   @assertEqual(real(dy_apr15, real64), TimeAmend(1), tolerance=1.0e-6_real64)
   @assertEqual(real(dy_jun20, real64), TimeAmend(2), tolerance=1.0e-6_real64)

   ! iamend pointers
   call assertEqual(1, iamend(1, 1), 'group 1: row 1')
   call assertEqual(2, iamend(2, 1), 'group 2: row 2')
   call assertEqual(3, iamend(2, 2), 'group 2: row 3')
end subroutine

@test
subroutine test_apply_events_no_file_zero_amendments()
   use nutrients_config_mod, only: nutrients_config_t
   use config_to_variables_mod, only: apply_nutrients_events
   use wofost_soil_declarations, only: namend, isme
   use funit
   implicit none
   type(nutrients_config_t) :: cfg

   ! Pre-existing non-zero values to detect overwrite.
   namend = 99
   isme   = 99

   ! cfg%events_file is unallocated by default.
   call apply_nutrients_events(cfg)

   call assertEqual(0, namend, 'namend reset to 0')
   call assertEqual(1, isme,   'isme reset to 1')
end subroutine

@test
subroutine test_apply_events_empty_file()
   use nutrients_config_mod, only: nutrients_config_t
   use config_to_variables_mod, only: apply_nutrients_events
   use variables, only: pathwork
   use wofost_soil_declarations, only: namend, isme
   use funit
   implicit none
   type(nutrients_config_t) :: cfg

   ! Empty-string events_file is treated like absent.
   pathwork = 'tests/unit/io/toml/fixtures/'
   cfg%events_file = ''
   namend = 99
   isme   = 99

   call apply_nutrients_events(cfg)

   call assertEqual(0, namend)
   call assertEqual(1, isme)
end subroutine

@test
subroutine test_apply_events_explicit_state_after_run()
   use, intrinsic :: iso_fortran_env, only: real64
   use nutrients_config_mod, only: nutrients_config_t
   use config_to_variables_mod, only: apply_nutrients_events
   use variables, only: pathwork
   use wofost_soil_declarations, only: namend, isme, NuAmend
   use funit
   implicit none
   type(nutrients_config_t) :: cfg

   ! Run on the fixture, then re-run with no file — verify reset.
   pathwork = 'tests/unit/io/toml/fixtures/'
   cfg%events_file = 'nutrients_events_small.csv'
   call apply_nutrients_events(cfg)
   call assertEqual(2, namend, 'after fixture run: 2 groups')

   ! Re-init: deallocate events_file, expect reset.
   if (allocated(cfg%events_file)) deallocate(cfg%events_file)
   call apply_nutrients_events(cfg)
   call assertEqual(0, namend, 'after no-file run: namend reset to 0')
   call assertEqual(1, isme,   'after no-file run: isme reset to 1')
end subroutine
```

The day-number values (`37726`, `37792`) match the days-since-1900 convention used elsewhere in the codebase (verified in `apply_irrigation_ssdi`'s tests, commit `328db06`). If `read_csv_table`'s date parsing produces a slightly different epoch by 1 day, adjust the literals; the equality should hold within tolerance `1.0e-6`.

Validation-rejection tests (material out-of-range, etc.) are deliberately omitted — the adapter calls `fatalerr_collected` on bad input, which means `error stop` exit-1 from a pFUnit context. Catching that requires harness gymnastics that aren't worth the test value. Manual verification at N3 against fixtures with deliberately-bad data is sufficient.

- [ ] **Step 2: Run to verify the tests pass**

```
pixi run -e test build-linux
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:|tests,"
pixi run -e test check-full 2>&1 | grep "Results:"
```
Expected: build clean; pFUnit `Ok: 1, Fail: 0` (test count grows by 4); check-full `5 passed, 0 failed`.

If any test fails on a date-numeric assertion, double-check the `read_csv_table` ISO-date convention by looking at the SSDI adapter test (`tests/unit/io/toml/test_apply_irrigation_ssdi.pf`) and copy the exact day-number used there for `2003-04-15`.

- [ ] **Step 3: Commit**

```bash
git add tests/unit/io/toml/test_nutrients_config.pf
git commit -m "$(cat <<'EOF'
test(io): pFUnit adapter tests for apply_nutrients_events

Four tests exercising the CSV-staging path:

- test_apply_events_csv_populates_legacy_globals: 3-row
  fixture (one same-day pair); asserts namend=2, NuAmend(1)=1,
  NuAmend(2)=2, MatNum/Amend/VolaFrac correctly populated and
  sorted, TimeAmend matches days-since-1900, iamend pointers
  correct.
- test_apply_events_no_file_zero_amendments: events_file
  unallocated; namend=0, isme=1.
- test_apply_events_empty_file: events_file = ''; same as
  unallocated.
- test_apply_events_explicit_state_after_run: re-running with
  deallocated events_file resets namend/isme.

Validation-rejection tests (out-of-range material, etc.) are
deliberately omitted — fatalerr_collected aborts via error stop,
which is hard to catch in pFUnit. Manual verification at N3.

Verified: build clean; pFUnit Ok: 1, Fail: 0 (count +4);
check-full 5/5.

Part of [nutrients] N2b (spec
docs/superpowers/specs/2026-05-08-nutrients-N2b-timed-amendments-design.md).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: ADR 0027 + index update + schema doc

**Files:**
- Create: `docs/adr/0027-nutrients-N2b-timed-amendments.md`
- Modify: `docs/adr/index.md` (append ADR 0027 row)
- Modify: `docs/configuration-schema.md` (document `events_file` + CSV format)

- [ ] **Step 1: Write `docs/adr/0027-nutrients-N2b-timed-amendments.md`**

```markdown
---
title: "ADR 0027 — [nutrients] N2b: timed amendments via CSV companion"
date: 2026-05-08
status: accepted
---

# ADR 0027: [nutrients] N2b — timed amendments via CSV companion

## Context

ADR 0025 (N1) wired the crop-side nutrient parameters; ADR 0026
(N2a) wired the soil-side initial state and `SorpCoef`. The
remaining gap before the runtime gate at `tillage.f90:73` can
be lifted is the timed soil management events — fertilizer
applications and manure spreading driven by the legacy
`<project>.sme` file.

Pre-flight inventory found:
- `<project>.sme` legacy schema: `smedate` (date), `MatNum`
  (1..20), `Dosagekgha` (0..500000), `VolatFraction` (legacy
  range 0..500000 — typo; physically a fraction in [0, 1]).
  Up to `maxamn = 1000` rows.
- The pre-collapse `SoilManagement(1)` block (deleted in SS-C
  step 3, commit `67f2545`) sorted events by date, converted
  kg/ha → kg/m², and grouped same-day dosages into
  `TimeAmend(:)`, `NuAmend(:)`, `iamend(:,:)`.
- The runtime path (`SoilManagement(3)` → `Wofost_SoilAmendents`)
  consumes those grouped globals; no code change needed there.

User feedback during brainstorming: "the toml would explode" —
inline `[[nutrients.events]]` arrays were rejected in favour
of a CSV companion. Same shape as SSDI events (ADR 0022).

## Decision

`nutrients_config_t` gains an optional
`character(len=:), allocatable :: events_file` field. The TOML
reader parses `[nutrients].events_file = "amendments.csv"`.
When non-empty, the adapter resolves the path against
`pathwork` and stages a 4-column CSV
(`date,material,amount_kgha,volat_fraction`) via
`read_csv_table` (column 1's ISO date auto-converts to
days-since-1900 per ADR 0012).

`apply_nutrients_events` (new private helper inside
`config_to_variables_mod`, called from the existing
`apply_nutrients`) does:

1. Default to `namend = 0`, `isme = 1` (no amendments).
2. If `events_file` allocated and non-empty: load CSV,
   validate per-row ranges, sort in-place by date (bubble
   sort), populate `MatNum(:)`, `Amend(:)` (kg/ha → kg/m²),
   `VolaFrac(:)`, then group same-day dosages into
   `TimeAmend(:)`, `NuAmend(:)`, `iamend(:,:)`, set `namend`,
   `isme = 1`.
3. Sort + group algorithm is a verbatim port of the deleted
   pre-collapse `SoilManagement(1)` legacy block.

The `volat_fraction` validator range tightens from legacy's
typo'd `[0, 500000]` to physically meaningful `[0, 1]`.

## Schema

```toml
[nutrients]
sorp_coef   = 0.005
events_file = "amendments.csv"

[nutrients.initial]
fom  = [0.5, 0.3, 0.2, 0.1, 0.5, 0.3, 0.2, 0.1]
bio  = 0.4
hum  = 8.0
cnh4 = 0.001
cno3 = 0.005
```

```
date,material,amount_kgha,volat_fraction
2003-04-15,10,100.0,0.05
2003-06-20,1,25000.0,0.10
2003-06-20,3,500.0,0.02
```

## Consequences

- With N1 + N2a + N2b in place, every soil-side and crop-side
  nutrient input flows through TOML + CSV. No legacy file
  reads remain on the nutrient-init path.
- `apply_nutrients` is called unconditionally per N2a; the
  `events_file` branch is also unconditional but takes the
  no-op path when the field is unallocated. check-full
  byte-identical.
- N3 (next sub-arc) lifts `tillage.f90:73`'s runtime
  stub-error and adds a regression case with populated
  `[nutrients]`. Verifying byte-identical against the legacy
  binary is N3's correctness gate.
- `volat_fraction` range tightening is a behavioural change.
  Legacy values outside `[0, 1]` were physically meaningless;
  ADR 0027 documents the intentional tightening.

## What N2b does NOT do

- **Lift `tillage.f90:73`.** N3.
- **Add a regression case.** N3.
- **Material-property overrides** (legacy `<project>.smm`).
  Out of scope across the whole `[nutrients]` umbrella.
- **`error_collection_t` consolidation.** Adapter validation
  routes through `fatalerr_collected` — matches SSDI mode-0
  precedent. Consolidation is a separate architectural arc.
- **Range-validate `events_file` shape at validate-time.**
  Per-row CSV validation lives in the adapter (file isn't
  parsed at config-validate time).

## Acceptance

- `pixi run -e test test-pfunit` → `Ok: 1, Fail: 0`; suite
  count grows by 6.
- `pixi run -e test check-full` → `5 passed, 0 failed`
  (byte-identical CSVs).
- `grep -n "subroutine apply_nutrients_events" src/io/toml/config_to_variables.f90`
  → one match.
- CSV fixture committed at
  `tests/unit/io/toml/fixtures/nutrients_events_small.csv`.
- `tillage.f90:73` UNCHANGED — `grep -n "flCropNut.*not.*allowed" src/crop/tillage.f90` returns one match.
- ADR 0027 committed.

## Related

- ADR 0008 — error collection over fatalerr.
- ADR 0012 — CSV companion input files.
- ADR 0022 — SSDI events CSV companion (mirror shape).
- ADR 0024 — `dtutil.f90` shim — same architectural thread.
- ADR 0025 — [nutrients] N1 (crop-side adapter).
- ADR 0026 — [nutrients] N2a (soil-side initial state).
- Future: ADR 0028 (N3: lift runtime gate + regression case).
```

- [ ] **Step 2: Append the index row**

In `docs/adr/index.md`, after the ADR 0026 row, add:

```markdown
- [ADR 0027 — [nutrients] N2b: timed amendments via CSV companion](0027-nutrients-N2b-timed-amendments.html) — Optional `events_file` on `[nutrients]` references a 4-column CSV (`date,material,amount_kgha,volat_fraction`); adapter stages, sorts, groups same-day dosages, populates legacy globals. `volat_fraction` range tightened to [0, 1] (legacy typo was [0, 500000]). Runtime path unchanged.
```

(Read the file first to see the exact existing format and adapt if needed.)

- [ ] **Step 3: Document `events_file` in `configuration-schema.md`**

Open `docs/configuration-schema.md`. Find the existing `[nutrients]` section (added in N2a, commit `0052588`). Append:

- **`events_file`** (string, optional): relative path to a CSV companion of timed soil management events. Resolved against `[general.paths].work` (`pathwork` global). Absent or empty → no amendments.
- **CSV columns:**
  - `date` — ISO `YYYY-MM-DD`. Strictly ascending after sort (same-day rows are grouped).
  - `material` — integer in `[1, 20]`. Indexes the hardcoded materials in `Wofost_SoilParameters` (1=Cattle manure, 10=Mineral N fertilizer, etc.).
  - `amount_kgha` — real in `[0, 500000]` (kg/ha).
  - `volat_fraction` — real in `[0, 1]` (volatilization fraction).
- **Limits:** up to 1000 rows (`maxamn` parameter in `wofost_soil_declarations`).
- **Sorting + grouping:** the adapter sorts by date and groups same-day rows into single amendment events.

Match the existing sub-section formatting.

- [ ] **Step 4: Final verification**

```
pixi run -e test build-linux
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:|tests,"
pixi run -e test check-full 2>&1 | grep "Results:"
echo "=== acceptance ==="
grep -n "subroutine apply_nutrients_events" src/io/toml/config_to_variables.f90 || echo "FAIL"
ls tests/unit/io/toml/fixtures/nutrients_events_small.csv
ls docs/adr/0027-nutrients-N2b-timed-amendments.md
grep -n "flCropNut.*not.*allowed" src/crop/tillage.f90
```
Expected: build clean; pFUnit `Ok: 1, Fail: 0` (count +6 over the N2a baseline); check-full 5/5; subroutine grep returns one match; both files exist; the `flCropNut.*not.*allowed` grep returns one match (runtime gate intact).

- [ ] **Step 5: Commit**

```bash
git add docs/adr/0027-nutrients-N2b-timed-amendments.md docs/adr/index.md docs/configuration-schema.md
git commit -m "$(cat <<'EOF'
docs: ADR 0027 — [nutrients] N2b timed amendments via CSV companion

Captures the third sub-arc of the [nutrients] umbrella. Adds an
optional events_file field on nutrients_config_t pointing at a
4-column CSV companion (date, material, amount_kgha,
volat_fraction). Adapter stages, validates per-row, sorts by
date, groups same-day dosages, populates legacy globals.

Sort + group algorithm is a verbatim port of the deleted
pre-collapse SoilManagement(1) legacy block. Runtime path
(SoilManagement(3) -> Wofost_SoilAmendents) unchanged.

volat_fraction range tightened to [0, 1] (legacy typo was
[0, 500000]).

Runtime gate at tillage.f90:73 stays — N3's job.

- docs/adr/0027-nutrients-N2b-timed-amendments.md: new ADR.
- docs/adr/index.md: ADR 0027 row.
- docs/configuration-schema.md: documents events_file + CSV
  format under [nutrients].

With N1 + N2a + N2b complete, every soil-side and crop-side
nutrient input flows through TOML + CSV. Closes the
[nutrients] N2b spec
(docs/superpowers/specs/2026-05-08-nutrients-N2b-timed-amendments-design.md).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Self-Review

**Spec coverage:**

| Spec section | Task |
|---|---|
| §Schema (`events_file` field) | Task 1 |
| §TOML shape (parser) | Task 1 |
| §CSV companion shape | Task 2 (fixture) + Task 3 (adapter reads it) |
| §Adapter (`apply_nutrients_events`) | Task 3 |
| §Tests (6 pFUnit) | Tasks 1, 4 (2 + 4 = 6) |
| §ADR 0027 + index + schema doc | Task 5 |

All spec sections covered. Validation-rejection tests omitted by design (catching `fatalerr_collected` from pFUnit is harness-fragile; manual verification at N3).

**Acceptance criteria** (from spec):

- `pixi run -e test build-linux` clean — verified each task.
- `pixi run -e test test-pfunit` +6 tests — verified Tasks 1, 4, 5.
- `pixi run -e test check-full` 5/5 byte-identical — verified Tasks 1, 3, 4.
- `subroutine apply_nutrients_events` exists — Task 3.
- CSV fixture committed — Task 2.
- `tillage.f90:73` UNCHANGED — Task 5 acceptance grep.
- ADR 0027 committed — Task 5.

**Type / signature consistency:**

- `events_file` field name (`character(len=:), allocatable`) consistent across Tasks 1, 3, 4.
- `apply_nutrients_events(cfg)` signature consistent.
- Module-variable names (`MatNum`, `Amend`, `VolaFrac`, `TimeAmend`, `NuAmend`, `iamend`, `namend`, `isme`, `maxamn`) match `wofost_soil_declarations.f90` (verified at spec-authoring; the case-style is lowercase per the file).
- CSV header strings (`date`, `material`, `amount_kgha`, `volat_fraction`) consistent in Tasks 2, 3, 4.

**Open questions punted to implementation:**

- Whether `get_optional_string_with_default` leaves the LHS unallocated when the key is absent (Task 1 Step 4 documents the contingency: adjust `test_read_nutrients_no_events_file`'s assertion if the helper unconditionally allocates).
- Exact case of `iamend` global (line 61 of `wofost_soil_declarations.f90` — lowercase per pre-flight grep, but verify at impl time).
- Whether `apply_nutrients_events` needs to be public or stays private (Task 3 Step 3 notes it depends on Task 4's pFUnit imports — the test suite imports it via `only:` so it must be public).
- Day-number literals (`37726`, `37792`) for the fixture dates — Task 4 Step 2 says cross-check against the SSDI test if assertions fail.

These are pickable at execution time without re-planning.

---

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-05-08-nutrients-N2b-timed-amendments.md`.
