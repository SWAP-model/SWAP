# SSDI TOML port — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Port the SSDI subsystem from TTutil-based `swap.swp` + `ssdi_file` reads to a typed `[irrigation.ssdi]` TOML block + CSV companion for fixed-mode events; retire `read_ssdi_input` and `src/io/checkdate.f90`.

**Architecture:** Add three nested types to `irrigation_config_t` (`irrigation_ssdi_fixed_t`, `irrigation_ssdi_scheduled_t`, `irrigation_ssdi_t`) holding the parameters previously read by `read_ssdi_input`. Mode is discriminated by an explicit `schedule = 0|1` field with mode-specific sub-tables. Fixed-mode event tables (up to 366 rows) live in a CSV companion. Add an `apply_irrigation_ssdi` helper to `config_to_variables.f90` that reads the CSV (mode 0) or copies the typed sub-block (mode 1), populates the legacy globals (`ssdi_date`, `ssdi_rate_f`, `ssdi_amount_f`, `ifnd`, `nod_ssdi`, `ssdi_sched_type`, `ssdi_threshold`, `ssdi_threshold_z`, `ssdi_amount`, `ssdi_appl_rate`, `sw_interval`, `days_interval`, `nirri`, `dt_SSDI_event`, `qssdi`), and runs the deferred date-window check that used to live in the deleted `checkdate`. Delete `read_ssdi_input`, collapse `SSDI_irrigation(1)`, and delete `src/io/checkdate.f90`.

**Tech Stack:** Fortran 2008, meson + ninja build, toml-f for parsing, pFUnit for unit tests, `csv_reader_mod` for CSV parsing (auto-converts column-1 ISO dates to days-since-1900). Reference patterns: `src/io/toml/read_irrigation_toml.f90` (existing irrigation reader, supports both inline + CSV companion), `src/io/toml/config_to_variables.f90` (look at the `fixed_events_file` block — `swirfix == 1 .and. allocated(...fixed_events_file)`, ~lines 184-220 — for the canonical CSV-staging adapter pattern).

**Spec:** `docs/superpowers/specs/2026-05-07-ssdi-toml-port-design.md`.

---

## File Structure

**Created:**

| Path | Responsibility |
|---|---|
| `src/io/toml/read_irrigation_ssdi_toml.f90` | New module `read_irrigation_ssdi_toml_mod` parsing `[irrigation.ssdi]` (scalars + sub-tables `[irrigation.ssdi.fixed]` / `[irrigation.ssdi.scheduled]`). Sibling to `read_irrigation_toml.f90`. |
| `tests/unit/io/toml/test_read_irrigation_ssdi_toml.pf` | Schema parsing tests (mode 0, mode 1, missing block, optional fields). |
| `tests/unit/io/toml/test_irrigation_ssdi_validate.pf` | Validator rule tests. |
| `tests/unit/io/toml/test_apply_irrigation_ssdi.pf` | Adapter tests for both modes (mode 0 with CSV fixture; mode 1 inline). |
| `tests/unit/io/toml/fixtures/ssdi_events_small.csv` | 3-row CSV fixture for the adapter test. |
| `docs/adr/0022-ssdi-toml-port.md` | New ADR. |

**Modified:**

| Path | Change |
|---|---|
| `src/config/irrigation_config.f90` | Add `irrigation_ssdi_fixed_t`, `irrigation_ssdi_scheduled_t`, `irrigation_ssdi_t`; add `ssdi` field on `irrigation_config_t`; extend `irrigation_config_validate`. |
| `src/io/toml/read_irrigation_toml.f90` | After parsing the `[irrigation]` table, call `read_irrigation_ssdi_toml`. |
| `src/io/toml/config_to_variables.f90` | Add `apply_irrigation_ssdi(ssdi)` helper; call it from `config_to_variables` when `flSSDI` is true. Drop the existing `swpfile`/`logf` HACK comments at the swssdi line if they're now dead. |
| `src/crop/irrigation.f90` | Delete `subroutine read_ssdi_input`; collapse `SSDI_irrigation(1)` body to a single `return`; drop unused locals (`swp`, `swpfile`, `logf`, `ssdi_file`, related TTutil reads). |
| `meson.build` | Add `'src/io/toml/read_irrigation_ssdi_toml.f90'`; drop `'src/io/checkdate.f90'`. |
| `tests/unit/meson.build` | Add the new reader to `pfunit_extra_sources`; drop `'../../src/io/checkdate.f90'`; add the three new `.pf` files to `pf_files`. |
| `tests/unit/testSuites.inc` | Register the three new test suites. |
| `docs/configuration-schema.md` | Document `[irrigation.ssdi]` (both modes) + the events-CSV format. |
| `docs/adr/index.md` | Add ADR 0022 row. |

**Deleted:**

| Path | Reason |
|---|---|
| `src/io/checkdate.f90` | Only caller (`read_ssdi_input`) is gone; date-window check moves into `apply_irrigation_ssdi` validator pass, expressed against `error_collection_t`. |

---

## Task 1: Schema types

**Files:**
- Modify: `src/config/irrigation_config.f90` (add type definitions before `irrigation_config_t`; add `ssdi` field on `irrigation_config_t`; export new types)
- Test: `tests/unit/io/toml/test_apply_irrigation_ssdi.pf` (default-init test added here; suite grows in Task 4)
- Modify: `tests/unit/meson.build` (register new `.pf` file)
- Modify: `tests/unit/testSuites.inc` (register the new suite — Task 1 of the tillage arc set this precedent)

- [ ] **Step 1: Write the failing test**

Create `tests/unit/io/toml/test_apply_irrigation_ssdi.pf`:

```fortran
@test
subroutine test_irrigation_ssdi_default_init()
   use irrigation_config_mod, only: irrigation_config_t
   use funit
   implicit none
   type(irrigation_config_t) :: cfg

   ! Default-init values for the new ssdi sub-type.
   call assertEqual(0, cfg%ssdi%schedule, 'schedule default')
   @assertEqual(0.0d0, cfg%ssdi%ssdi_z(1), tolerance=1.0d-12)
   @assertEqual(0.0d0, cfg%ssdi%ssdi_z(2), tolerance=1.0d-12)
   @assertFalse(allocated(cfg%ssdi%fixed%events_file), 'events_file not allocated by default')
   call assertEqual(0, cfg%ssdi%scheduled%sched_type, 'sched_type default')
   call assertEqual(1, cfg%ssdi%scheduled%days_interval, 'days_interval default')
end subroutine
```

- [ ] **Step 2: Wire the test file and run to verify it fails**

Add to `tests/unit/meson.build`'s `pf_files` list:
```meson
        'io/toml/test_apply_irrigation_ssdi.pf',
```

Add to `tests/unit/testSuites.inc`:
```
ADD_TEST_SUITE(test_apply_irrigation_ssdi_suite)
```

Run:
```
pixi run -e test build-linux
```
Expected: build fails — `cfg%ssdi` is not a member of `irrigation_config_t`.

- [ ] **Step 3: Add the new types in `src/config/irrigation_config.f90`**

Insert the three types **before** the `type :: irrigation_config_t` declaration at line 26 (immediately after the existing types). Add the new types to the module's `public` exports too (look for the existing `public ::` block).

```fortran
   public :: irrigation_ssdi_fixed_t
   public :: irrigation_ssdi_scheduled_t
   public :: irrigation_ssdi_t
```

```fortran
   !> Fixed-mode SSDI configuration. Date table lives in a CSV file
   !! pointed to by events_file (parsed by the adapter, not the TOML
   !! reader).
   type :: irrigation_ssdi_fixed_t
      character(len=:), allocatable :: events_file   !! relative path to ssdi_events.csv
   end type irrigation_ssdi_fixed_t

   !> Scheduled-trigger-mode SSDI configuration.
   type :: irrigation_ssdi_scheduled_t
      integer      :: sched_type      = 0          !! 1=Tred, 2=presh, 3=watc
      real(real64) :: threshold       = 0.0_real64 !! semantics determined by sched_type
      real(real64) :: threshold_depth = 0.0_real64 !! cm; required when sched_type > 1
      real(real64) :: ssdi_amount     = 0.0_real64 !! mm/event
      real(real64) :: ssdi_appl_rate  = 0.0_real64 !! mm/h
      integer      :: sw_interval     = 0          !! 0=daily, 1=interval-gated
      integer      :: days_interval   = 1          !! 1..366; required when sw_interval=1
   end type irrigation_ssdi_scheduled_t

   !> [irrigation.ssdi] block container. Mode discriminator: schedule.
   type :: irrigation_ssdi_t
      integer      :: schedule  = 0                  !! 0=fixed-date, 1=scheduled-trigger
      real(real64) :: ssdi_z(2) = 0.0_real64         !! cm; both elements equal for single-depth
      type(irrigation_ssdi_fixed_t)     :: fixed
      type(irrigation_ssdi_scheduled_t) :: scheduled
   end type irrigation_ssdi_t
```

Add the field on `irrigation_config_t` (alongside other nested members):
```fortran
      type(irrigation_ssdi_t) :: ssdi
```

- [ ] **Step 4: Build and run the test to verify it passes**

```
pixi run -e test build-linux
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:|tests,"
```
Expected: build clean; pFUnit `Ok: 1, Fail: 0` (tests count grows by 1).

- [ ] **Step 5: Commit**

```bash
git add src/config/irrigation_config.f90 tests/unit/io/toml/test_apply_irrigation_ssdi.pf tests/unit/meson.build tests/unit/testSuites.inc
git commit -m "$(cat <<'EOF'
feat(config): add irrigation_ssdi_t types to irrigation_config_t

Adds irrigation_ssdi_fixed_t, irrigation_ssdi_scheduled_t, and
irrigation_ssdi_t containers to irrigation_config_mod, plus the
`ssdi` field on irrigation_config_t. Defaults are off-state (no
events_file allocated; schedule=0; sched_type=0). pFUnit smoke
test asserts default-init.

No reader, validator, or adapter wired yet — those land in
subsequent tasks.

Part of SSDI TOML port (spec
docs/superpowers/specs/2026-05-07-ssdi-toml-port-design.md).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Validator extension

**Files:**
- Modify: `src/config/irrigation_config.f90` (extend `irrigation_config_validate`)
- Test: `tests/unit/io/toml/test_irrigation_ssdi_validate.pf` (new)
- Modify: `tests/unit/meson.build`, `tests/unit/testSuites.inc` (register the new suite)

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/io/toml/test_irrigation_ssdi_validate.pf`. Write **all of the following test subroutines concretely** in the file. Each follows the same shape as the tillage validate tests; reference `tests/unit/io/toml/test_soil_tillage_validate.pf` for the exact pattern (build a `soil_config_t` analog → call `validate` → assert error count or absence).

Tests to include:

1. `test_validate_swssdi1_schedule_out_of_range` — `cfg%swssdi=1`, `cfg%ssdi%schedule = 5` → expect error.
2. `test_validate_mode0_requires_events_file` — `schedule=0`, `events_file` empty → expect error.
3. `test_validate_mode0_rejects_populated_scheduled_block` — `schedule=0`, `events_file = "x.csv"`, `scheduled%sched_type = 1` → expect error (fields must stay at defaults).
4. `test_validate_mode1_rejects_populated_events_file` — `schedule=1`, `events_file = "x.csv"` → expect error.
5. `test_validate_mode1_sched_type_in_range` — `schedule=1`, `sched_type=4` → expect error.
6. `test_validate_mode1_threshold_range_per_sched_type` — three sub-tests:
   - `sched_type=1, threshold=2.0` → expect error (range [0,1])
   - `sched_type=2, threshold=1.0` → expect error (range [-1e7,0])
   - `sched_type=3, threshold=-0.1` → expect error (range [0,1])
7. `test_validate_mode1_threshold_depth_required_when_sched_type_gt_1` — `sched_type=2`, `threshold_depth=10.0` → expect error (must be in [-100, 0]).
8. `test_validate_mode1_amount_appl_rate_in_range` — `ssdi_amount=200.0` (mm) → expect error.
9. `test_validate_mode1_days_interval_required_when_sw_interval_1` — `sw_interval=1`, `days_interval=0` → expect error.
10. `test_validate_ssdi_z_ordering` — `ssdi_z=[-50.0d0, -30.0d0]` (top below bottom) → expect error.
11. `test_validate_swssdi0_skips_block` — `swssdi=0`, no ssdi block → expect zero new errors.
12. `test_validate_clean_mode0_no_errors` — `swssdi=1`, `schedule=0`, `events_file = "ssdi.csv"`, `ssdi_z=[-30.0d0, -50.0d0]`, scheduled struct at defaults → expect zero errors.
13. `test_validate_clean_mode1_no_errors` — `swssdi=1`, `schedule=1`, `events_file` empty, scheduled fully populated in-range → expect zero errors.

Skeleton for the first two:

```fortran
@test
subroutine test_validate_swssdi1_schedule_out_of_range()
   use irrigation_config_mod, only: irrigation_config_t, irrigation_config_validate
   use error_mod, only: error_collection_t
   use funit
   implicit none
   type(irrigation_config_t) :: cfg
   type(error_collection_t)  :: errs

   cfg%swssdi = 1
   cfg%ssdi%schedule = 5    ! out of [0, 1]
   cfg%ssdi%fixed%events_file = ''     ! intentional — separate failing rule but check schedule first
   call irrigation_config_validate(cfg, errs)
   @assertTrue(errs%count() > 0, 'expected error for schedule=5')
end subroutine

@test
subroutine test_validate_mode0_requires_events_file()
   use irrigation_config_mod, only: irrigation_config_t, irrigation_config_validate
   use error_mod, only: error_collection_t
   use funit
   implicit none
   type(irrigation_config_t) :: cfg
   type(error_collection_t)  :: errs

   cfg%swssdi = 1
   cfg%ssdi%schedule = 0
   ! events_file deliberately left at default (unallocated) → expect error
   cfg%ssdi%ssdi_z = [-30.0d0, -50.0d0]
   call irrigation_config_validate(cfg, errs)
   @assertTrue(errs%count() > 0, 'expected error for schedule=0 with no events_file')
end subroutine
```

Write the rest in the same shape.

- [ ] **Step 2: Wire and run to verify failures**

Add to `tests/unit/meson.build`'s `pf_files`:
```meson
        'io/toml/test_irrigation_ssdi_validate.pf',
```

Add to `tests/unit/testSuites.inc`:
```
ADD_TEST_SUITE(test_irrigation_ssdi_validate_suite)
```

Run:
```
pixi run -e test build-linux
pixi run -e test test-pfunit 2>&1 | grep -E "FAIL|Fail:"
```
Expected: build clean; tests fail (validator hasn't been extended).

- [ ] **Step 3: Extend `irrigation_config_validate` in `src/config/irrigation_config.f90`**

Find `subroutine irrigation_config_validate` and `end subroutine irrigation_config_validate`. Insert a new `if (self%swssdi == 1)` block just before the `end subroutine` line. Use the same error-code constants as the tillage validator: `ERR_VALIDATION_OUT_OF_RANGE` for ranges, `ERR_VALIDATION_CROSS_FIELD` for cross-field invariants. Reference the tillage block in `src/config/soil_config.f90` for the exact pattern.

```fortran
      ! [irrigation.ssdi] validation — only fires when swssdi = 1.
      if (self%swssdi == 1) then
         ! Mode discriminator
         call check_int_range(self%ssdi%schedule, 0, 1, &
                              'irrigation.ssdi.schedule', errors)

         ! Shared: ssdi_z range + ordering (top above or equal bottom; cm depths are negative)
         call check_real_range(self%ssdi%ssdi_z(1), -100.0_real64, 0.0_real64, &
                               'irrigation.ssdi.ssdi_z', errors)
         call check_real_range(self%ssdi%ssdi_z(2), -100.0_real64, 0.0_real64, &
                               'irrigation.ssdi.ssdi_z', errors)
         if (self%ssdi%ssdi_z(1) < self%ssdi%ssdi_z(2)) then
            call errors%append(ERR_VALIDATION_CROSS_FIELD, &
                               'ssdi_z(1) must be >= ssdi_z(2) (top above or equal bottom; cm depths are negative)', &
                               'irrigation.ssdi.ssdi_z')
         end if

         ! Mode-specific
         select case (self%ssdi%schedule)
         case (0)
            ! Fixed mode: events_file required, scheduled-block at defaults
            if (.not. allocated(self%ssdi%fixed%events_file) .or. &
                len_trim(self%ssdi%fixed%events_file) == 0) then
               call errors%append(ERR_VALIDATION_CROSS_FIELD, &
                                  'must be non-empty when schedule = 0', &
                                  'irrigation.ssdi.fixed.events_file')
            end if
            if (self%ssdi%scheduled%sched_type /= 0 .or. &
                self%ssdi%scheduled%threshold /= 0.0_real64 .or. &
                self%ssdi%scheduled%ssdi_amount /= 0.0_real64 .or. &
                self%ssdi%scheduled%ssdi_appl_rate /= 0.0_real64) then
               call errors%append(ERR_VALIDATION_CROSS_FIELD, &
                                  '[irrigation.ssdi.scheduled] must be at defaults when schedule = 0 (silently-ignored fields are a footgun)', &
                                  'irrigation.ssdi.scheduled')
            end if

         case (1)
            ! Scheduled mode: events_file empty, scheduled-block populated
            if (allocated(self%ssdi%fixed%events_file) .and. &
                len_trim(self%ssdi%fixed%events_file) > 0) then
               call errors%append(ERR_VALIDATION_CROSS_FIELD, &
                                  'must be empty when schedule = 1', &
                                  'irrigation.ssdi.fixed.events_file')
            end if

            call check_int_range(self%ssdi%scheduled%sched_type, 1, 3, &
                                 'irrigation.ssdi.scheduled.sched_type', errors)
            select case (self%ssdi%scheduled%sched_type)
            case (1)
               call check_real_range(self%ssdi%scheduled%threshold, 0.0_real64, 1.0_real64, &
                                     'irrigation.ssdi.scheduled.threshold', errors)
            case (2)
               call check_real_range(self%ssdi%scheduled%threshold, -1.0e7_real64, 0.0_real64, &
                                     'irrigation.ssdi.scheduled.threshold', errors)
            case (3)
               call check_real_range(self%ssdi%scheduled%threshold, 0.0_real64, 1.0_real64, &
                                     'irrigation.ssdi.scheduled.threshold', errors)
            end select
            if (self%ssdi%scheduled%sched_type > 1) then
               call check_real_range(self%ssdi%scheduled%threshold_depth, -100.0_real64, 0.0_real64, &
                                     'irrigation.ssdi.scheduled.threshold_depth', errors)
            end if
            call check_real_range(self%ssdi%scheduled%ssdi_amount, 0.0_real64, 100.0_real64, &
                                  'irrigation.ssdi.scheduled.ssdi_amount', errors)
            call check_real_range(self%ssdi%scheduled%ssdi_appl_rate, 0.0_real64, 100.0_real64, &
                                  'irrigation.ssdi.scheduled.ssdi_appl_rate', errors)
            call check_int_range(self%ssdi%scheduled%sw_interval, 0, 1, &
                                 'irrigation.ssdi.scheduled.sw_interval', errors)
            if (self%ssdi%scheduled%sw_interval == 1) then
               call check_int_range(self%ssdi%scheduled%days_interval, 1, 366, &
                                    'irrigation.ssdi.scheduled.days_interval', errors)
            end if
         end select
      end if
```

You may need to add `use validation_mod, only: check_int_range, check_real_range` if it's not already in scope (it is in `soil_config.f90`; verify the same import pattern exists in `irrigation_config.f90`).

- [ ] **Step 4: Build and run tests**

```
pixi run -e test build-linux
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:|tests,"
```
Expected: build clean; `Ok: 1, Fail: 0`.

- [ ] **Step 5: Commit**

```bash
git add src/config/irrigation_config.f90 tests/unit/io/toml/test_irrigation_ssdi_validate.pf tests/unit/meson.build tests/unit/testSuites.inc
git commit -m "$(cat <<'EOF'
feat(config): validate [irrigation.ssdi] when swssdi=1

Extends irrigation_config_validate with rules for the new SSDI
block: schedule discriminator [0..1], ssdi_z range and ordering,
mode-mismatch (events_file required for mode 0; scheduled-block
at defaults; events_file empty for mode 1), sched_type [1..3],
threshold range per sched_type (mode-specific), threshold_depth
required when sched_type > 1, ssdi_amount/appl_rate ranges,
sw_interval enum, days_interval required when sw_interval=1.

Date-window check (CSV dates against [tstart, tend]) is deferred
to the adapter (Task 4) since validator can't see the staged CSV
rows.

Validator stays silent when swssdi=0.

pFUnit suite test_irrigation_ssdi_validate exercises every rule.

Part of SSDI TOML port (spec
docs/superpowers/specs/2026-05-07-ssdi-toml-port-design.md).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: TOML reader

**Files:**
- Create: `src/io/toml/read_irrigation_ssdi_toml.f90`
- Modify: `src/io/toml/read_irrigation_toml.f90` (call new sibling)
- Modify: `meson.build`, `tests/unit/meson.build` (`pfunit_extra_sources`)
- Test: `tests/unit/io/toml/test_read_irrigation_ssdi_toml.pf` (new)
- Modify: `tests/unit/meson.build`, `tests/unit/testSuites.inc` (register the new suite)

- [ ] **Step 1: Write the failing test**

Create `tests/unit/io/toml/test_read_irrigation_ssdi_toml.pf`. Reference `tests/unit/io/toml/test_read_soil_tillage_toml.pf` for the in-memory TOML loading idiom (the tillage spec settled on the project's standard pattern; mirror it).

Three tests:

1. `test_read_mode0_full_block`: TOML containing
   ```toml
   [irrigation]
   swssdi = 1
   [irrigation.ssdi]
   schedule = 0
   ssdi_z   = [-30.0, -50.0]
   [irrigation.ssdi.fixed]
   events_file = "ssdi_events.csv"
   ```
   Assertions: `ssdi%schedule == 0`, `ssdi%ssdi_z(1) == -30.0d0`, `ssdi%ssdi_z(2) == -50.0d0`, `allocated(ssdi%fixed%events_file)`, `trim(ssdi%fixed%events_file) == "ssdi_events.csv"`, `ssdi%scheduled%sched_type == 0` (default).

2. `test_read_mode1_full_block`: TOML containing
   ```toml
   [irrigation]
   swssdi = 1
   [irrigation.ssdi]
   schedule = 1
   ssdi_z   = [-30.0, -50.0]
   [irrigation.ssdi.scheduled]
   sched_type = 1
   threshold = 0.7
   ssdi_amount = 10.0
   ssdi_appl_rate = 5.0
   sw_interval = 0
   ```
   Assertions: every populated scheduled-block field round-trips correctly; `events_file` not allocated; `days_interval == 1` (default).

3. `test_read_missing_block_keeps_defaults`: TOML with `[irrigation]` and `swssdi = 0`, no `[irrigation.ssdi]` block. Assertions: `ssdi%schedule == 0`, `.not. allocated(ssdi%fixed%events_file)`, `ssdi%scheduled%sched_type == 0`.

- [ ] **Step 2: Wire and run to verify failure**

Add to `tests/unit/meson.build`'s `pf_files`:
```meson
        'io/toml/test_read_irrigation_ssdi_toml.pf',
```

Add to `tests/unit/testSuites.inc`:
```
ADD_TEST_SUITE(test_read_irrigation_ssdi_toml_suite)
```

Run `pixi run -e test build-linux`. Expected: build fails — `read_irrigation_ssdi_toml_mod` doesn't exist.

- [ ] **Step 3: Create `src/io/toml/read_irrigation_ssdi_toml.f90`**

```fortran
!> @file read_irrigation_ssdi_toml.f90
!! Parses the optional [irrigation.ssdi] block into
!! irrigation_config_t%ssdi. Block omitted when irrigation.swssdi = 0;
!! reader leaves defaults intact in that case. Sibling to
!! read_irrigation_toml.f90 to keep that file focused.
module read_irrigation_ssdi_toml_mod

   use, intrinsic :: iso_fortran_env, only: real64
   use tomlf, only: toml_table, toml_array, get_value, len
   use toml_field_helpers_mod, only: get_optional_int_with_default,       &
                                     get_optional_real_with_default,      &
                                     get_optional_string_with_default
   use irrigation_config_mod, only: irrigation_ssdi_t,                    &
                                    irrigation_ssdi_fixed_t,              &
                                    irrigation_ssdi_scheduled_t
   use error_mod, only: error_collection_t, ERR_PARSE_TYPE_MISMATCH
   implicit none
   private

   public :: read_irrigation_ssdi_toml

contains

   !> Populate `ssdi` from the optional [irrigation.ssdi] table on `irr_sec`.
   !! Missing block is benign — leaves defaults.
   subroutine read_irrigation_ssdi_toml(irr_sec, ssdi, errors)
      type(toml_table), pointer, intent(in)    :: irr_sec
      type(irrigation_ssdi_t),   intent(inout) :: ssdi
      type(error_collection_t),  intent(inout) :: errors

      type(toml_table), pointer :: ssdi_tbl, fixed_tbl, sched_tbl
      type(toml_array), pointer :: zarr
      integer :: stat, n, k
      real(real64) :: zval

      if (.not. associated(irr_sec)) return

      ssdi_tbl => null()
      call get_value(irr_sec, 'ssdi', ssdi_tbl, requested=.false., stat=stat)
      if (stat /= 0 .or. .not. associated(ssdi_tbl)) return

      call get_optional_int_with_default(ssdi_tbl, 'schedule', ssdi%schedule, 0, &
                                         'irrigation.ssdi.schedule', errors)

      ! ssdi_z is a 2-element array; both elements equal for single-depth.
      zarr => null()
      call get_value(ssdi_tbl, 'ssdi_z', zarr, requested=.false., stat=stat)
      if (stat == 0 .and. associated(zarr)) then
         n = len(zarr)
         if (n == 1) then
            call get_value(zarr, 1, zval, stat=stat)
            if (stat == 0) then
               ssdi%ssdi_z(1) = zval
               ssdi%ssdi_z(2) = zval
            end if
         else if (n == 2) then
            do k = 1, 2
               call get_value(zarr, k, zval, stat=stat)
               if (stat == 0) ssdi%ssdi_z(k) = zval
            end do
         else
            call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                               'ssdi_z must be a 1- or 2-element array', &
                               'irrigation.ssdi.ssdi_z')
         end if
      end if

      ! [irrigation.ssdi.fixed]
      fixed_tbl => null()
      call get_value(ssdi_tbl, 'fixed', fixed_tbl, requested=.false., stat=stat)
      if (stat == 0 .and. associated(fixed_tbl)) then
         call get_optional_string_with_default(fixed_tbl, 'events_file', &
                                               ssdi%fixed%events_file, '', &
                                               'irrigation.ssdi.fixed.events_file', errors)
      end if

      ! [irrigation.ssdi.scheduled]
      sched_tbl => null()
      call get_value(ssdi_tbl, 'scheduled', sched_tbl, requested=.false., stat=stat)
      if (stat == 0 .and. associated(sched_tbl)) then
         call read_scheduled_fields(sched_tbl, ssdi%scheduled, errors)
      end if
   end subroutine read_irrigation_ssdi_toml


   subroutine read_scheduled_fields(sec, sched, errors)
      type(toml_table), pointer,  intent(in)    :: sec
      type(irrigation_ssdi_scheduled_t), intent(inout) :: sched
      type(error_collection_t),   intent(inout) :: errors

      call get_optional_int_with_default(sec, 'sched_type', sched%sched_type, 0, &
                                         'irrigation.ssdi.scheduled.sched_type', errors)
      call get_optional_real_with_default(sec, 'threshold', sched%threshold, 0.0_real64, &
                                          'irrigation.ssdi.scheduled.threshold', errors)
      call get_optional_real_with_default(sec, 'threshold_depth', sched%threshold_depth, 0.0_real64, &
                                          'irrigation.ssdi.scheduled.threshold_depth', errors)
      call get_optional_real_with_default(sec, 'ssdi_amount', sched%ssdi_amount, 0.0_real64, &
                                          'irrigation.ssdi.scheduled.ssdi_amount', errors)
      call get_optional_real_with_default(sec, 'ssdi_appl_rate', sched%ssdi_appl_rate, 0.0_real64, &
                                          'irrigation.ssdi.scheduled.ssdi_appl_rate', errors)
      call get_optional_int_with_default(sec, 'sw_interval', sched%sw_interval, 0, &
                                         'irrigation.ssdi.scheduled.sw_interval', errors)
      call get_optional_int_with_default(sec, 'days_interval', sched%days_interval, 1, &
                                         'irrigation.ssdi.scheduled.days_interval', errors)
   end subroutine read_scheduled_fields

end module read_irrigation_ssdi_toml_mod
```

- [ ] **Step 4: Wire into `read_irrigation_toml.f90`**

Add a `use` line near the top:
```fortran
   use read_irrigation_ssdi_toml_mod, only: read_irrigation_ssdi_toml
```

After the existing `[irrigation]` parsing in `subroutine read_irrigation_toml` (find by name; the existing parser already uses a `sec` pointer to the irrigation table — confirm before pasting), add:
```fortran
      call read_irrigation_ssdi_toml(sec, config%ssdi, errors)
```

- [ ] **Step 5: Add the file to meson sources**

In `meson.build` production source list (where `read_irrigation_toml.f90` is registered), add:
```meson
    'src/io/toml/read_irrigation_ssdi_toml.f90',
```

In `tests/unit/meson.build` `pfunit_extra_sources`, add:
```meson
        '../../src/io/toml/read_irrigation_ssdi_toml.f90',
```

- [ ] **Step 6: Build and run tests to verify pass**

```
pixi run -e test build-linux
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:|tests,"
```
Expected: build clean; `Ok: 1, Fail: 0`.

- [ ] **Step 7: Commit**

```bash
git add src/io/toml/read_irrigation_ssdi_toml.f90 src/io/toml/read_irrigation_toml.f90 meson.build tests/unit/meson.build tests/unit/testSuites.inc tests/unit/io/toml/test_read_irrigation_ssdi_toml.pf
git commit -m "$(cat <<'EOF'
feat(io): parse [irrigation.ssdi] into irrigation_config_t%ssdi

Adds read_irrigation_ssdi_toml_mod, called from read_irrigation_toml
after the existing [irrigation] parsing. Block is optional — when
[irrigation.ssdi] is absent, reader leaves defaults intact. Both
mode sub-tables ([irrigation.ssdi.fixed], [irrigation.ssdi.scheduled])
are read unconditionally; the validator (Task 2) enforces
mode-mismatch.

ssdi_z accepts a 1- or 2-element array; single-depth case populates
both array elements with the same value.

pFUnit suite test_read_irrigation_ssdi_toml covers mode 0, mode 1,
and missing block.

Part of SSDI TOML port (spec
docs/superpowers/specs/2026-05-07-ssdi-toml-port-design.md).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: `apply_irrigation_ssdi` adapter

**Files:**
- Modify: `src/io/toml/config_to_variables.f90` (add `apply_irrigation_ssdi`; call from `config_to_variables` when `flSSDI`)
- Test: `tests/unit/io/toml/test_apply_irrigation_ssdi.pf` (extend with adapter tests)
- Create: `tests/unit/io/toml/fixtures/ssdi_events_small.csv` (3-row CSV fixture)

This task is the most substantive. It does:
- Resolves `ssdi_z` array → `nod_ssdi(1:2)` layer indices via `zbotcp` walk.
- Mode 0: reads CSV via `read_csv_table` (column 1 auto-converts ISO dates → days-since-1900); validates row count, ascending dates, date-window (replaces `checkdate`); copies into `ssdi_date`/`ssdi_rate_f`/`ssdi_amount_f`; computes initial `nirri` from `tstart`; converts mm/h → cm/d, mm → cm; spreads `amount_f` over compartments.
- Mode 1: copies typed config to legacy globals; sets `nirri = 1`, `dt_SSDI_event = 1.0d0`, `days_counter = 366`; converts units; spreads `ssdi_amount` over compartments.
- Both: `qssdi = 0.0d0`.

Reference: `src/io/toml/config_to_variables.f90:184-220` (the `swirfix == 1 .and. allocated(...fixed_events_file)` block) — that's the canonical CSV-staging adapter pattern in this codebase. Path resolution: `trim(pathwork) // trim(events_file)`.

Pre-step: read `subroutine SSDI_irrigation` case(1) body in `src/crop/irrigation.f90` end-to-end (currently around lines 350-460) — that's the legacy code being replaced. The adapter must produce **byte-identical legacy globals** for swssdi=1 cases.

- [ ] **Step 1: Create the CSV fixture**

Create `tests/unit/io/toml/fixtures/ssdi_events_small.csv`:

```
date,rate_f,amount_f
2003-04-15,5.0,10.0
2003-05-02,5.0,8.0
2003-06-10,4.0,6.0
```

- [ ] **Step 2: Write the failing tests**

Append to `tests/unit/io/toml/test_apply_irrigation_ssdi.pf` (created in Task 1; already has the default-init test):

```fortran
@test
subroutine test_apply_mode0_populates_globals_from_csv()
   use irrigation_config_mod, only: irrigation_config_t
   use config_to_variables_mod, only: apply_irrigation_ssdi
   use variables, only: ssdi_date_irr, ssdi_rate_f_irr, ssdi_amount_f_irr,    &
                        ifnd, nod_ssdi_irr, qssdi, NumNod, zbotcp,            &
                        tend, tstart, t1900, pathwork, mairg
   use, intrinsic :: iso_fortran_env, only: real64
   use funit
   implicit none
   type(irrigation_config_t) :: cfg
   integer :: dy_apr15

   ! Minimal grid for ssdi_z resolution.
   NumNod = 4
   if (allocated(zbotcp)) deallocate(zbotcp); allocate(zbotcp(4))
   zbotcp = [-10.0d0, -25.0d0, -40.0d0, -60.0d0]

   ! Simulation window covers the fixture dates (2003-04-15 .. 2003-06-10).
   tstart = 37621.0_real64   ! 2003-01-01 (days since 1900)
   tend   = 37985.0_real64   ! 2003-12-31
   t1900  = tstart
   pathwork = 'tests/unit/io/toml/fixtures/'   ! relative to project root

   cfg%swssdi = 1
   cfg%ssdi%schedule = 0
   cfg%ssdi%ssdi_z = [-25.0d0, -40.0d0]
   cfg%ssdi%fixed%events_file = 'ssdi_events_small.csv'

   call apply_irrigation_ssdi(cfg%ssdi)

   @assertEqual(3, ifnd, 'three rows in fixture')
   ! Date columns parsed to days-since-1900; first row is 2003-04-15.
   dy_apr15 = 37726   ! days from 1900-01-01 to 2003-04-15
   @assertEqual(real(dy_apr15, real64), ssdi_date_irr(1), tolerance=1.0d-6)
   ! Unit conversion: rate_f mm/h -> cm/d  (5 * 0.1 * 24 = 12.0)
   @assertEqual(12.0_real64, ssdi_rate_f_irr(1), tolerance=1.0d-9)
   ! Unit conversion + spread: amount_f mm -> cm, then divided by (nod_ssdi(2)-nod_ssdi(1)+1)
   ! ssdi_z = [-25, -40] -> nod_ssdi(1)=2, nod_ssdi(2)=3 (2 nodes); 10mm -> 1.0cm / 2 = 0.5
   @assertEqual(0.5_real64, ssdi_amount_f_irr(1), tolerance=1.0d-9)
   @assertEqual(2, nod_ssdi_irr(1))
   @assertEqual(3, nod_ssdi_irr(2))
end subroutine

@test
subroutine test_apply_mode1_copies_scheduled_block()
   use irrigation_config_mod, only: irrigation_config_t
   use config_to_variables_mod, only: apply_irrigation_ssdi
   use variables, only: ssdi_sched_type_irr, ssdi_threshold_irr,        &
                        ssdi_threshold_z_irr, ssdi_amount_irr,          &
                        ssdi_appl_rate_irr, sw_interval_irr,            &
                        days_interval_irr, dt_SSDI_event, nirri,        &
                        qssdi, NumNod, zbotcp, nod_ssdi_irr
   use, intrinsic :: iso_fortran_env, only: real64
   use funit
   implicit none
   type(irrigation_config_t) :: cfg

   NumNod = 4
   if (allocated(zbotcp)) deallocate(zbotcp); allocate(zbotcp(4))
   zbotcp = [-10.0d0, -25.0d0, -40.0d0, -60.0d0]

   cfg%swssdi = 1
   cfg%ssdi%schedule = 1
   cfg%ssdi%ssdi_z = [-25.0d0, -40.0d0]
   cfg%ssdi%scheduled%sched_type = 1
   cfg%ssdi%scheduled%threshold  = 0.7d0
   cfg%ssdi%scheduled%ssdi_amount    = 10.0d0
   cfg%ssdi%scheduled%ssdi_appl_rate = 5.0d0
   cfg%ssdi%scheduled%sw_interval    = 0
   cfg%ssdi%scheduled%days_interval  = 1

   call apply_irrigation_ssdi(cfg%ssdi)

   @assertEqual(1,     ssdi_sched_type_irr)
   @assertEqual(0.7d0, ssdi_threshold_irr,  tolerance=1.0d-12)
   ! mm -> cm, then spread over (nod_ssdi(2)-nod_ssdi(1)+1) = 2 compartments
   @assertEqual(0.5d0, ssdi_amount_irr,     tolerance=1.0d-9)
   ! mm/h -> cm/d
   @assertEqual(12.0d0, ssdi_appl_rate_irr, tolerance=1.0d-9)
   @assertEqual(0,     sw_interval_irr)
   ! Regression-fix invariant from d8a88d6: dt_SSDI_event default is 1.0
   @assertEqual(1.0d0, dt_SSDI_event,       tolerance=1.0d-12)
   @assertEqual(1,     nirri)
   @assertEqual(2, nod_ssdi_irr(1))
   @assertEqual(3, nod_ssdi_irr(2))
end subroutine
```

Note on global names: this codebase uses `_irr`-suffixed legacy globals for SSDI parameters (see `src/crop/irrigation.f90:343-349` and the `variables` module). The actual exact list is:
- `ssdi_date_irr`, `ssdi_rate_f_irr`, `ssdi_amount_f_irr` — fixed-mode arrays
- `ssdi_sched_type_irr`, `ssdi_threshold_irr`, `ssdi_threshold_z_irr` — scheduled-mode scalars
- `ssdi_amount_irr`, `ssdi_appl_rate_irr`, `sw_interval_irr`, `days_interval_irr` — scheduled-mode scalars
- `nod_ssdi_irr(2)`, `nirri_ssdi_irr` — shared
- `dt_SSDI_event`, `qssdi`, `nirri`, `qssdi(:)` — non-suffixed globals

Verify exact names by `grep -n "_irr\b" src/core/variables.f90 | head -30` before writing the tests. If any name above is wrong, fix in-place.

- [ ] **Step 3: Run to verify failures**

```
pixi run -e test build-linux
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:|tests,"
```
Expected: build fails or tests fail (`apply_irrigation_ssdi` doesn't exist).

- [ ] **Step 4: Add `apply_irrigation_ssdi` to `src/io/toml/config_to_variables.f90`**

Add to the module's `public` exports:
```fortran
   public :: apply_irrigation_ssdi
```

Add to module-level `use`:
```fortran
   use irrigation_config_mod, only: irrigation_ssdi_t
```

Insert near the end of `contains` (after the existing tillage adapter from ADR 0021):

```fortran
   !> Apply [irrigation.ssdi] config to legacy `variables` globals.
   !! Called from config_to_variables when flSSDI is true. Reads
   !! the events CSV (mode 0) or copies the scheduled sub-block
   !! (mode 1), allocates and populates the legacy globals,
   !! resolves ssdi_z to layer indices, and runs the deferred
   !! date-window check (replaces the deleted checkdate).
   !!
   !! Replaces SSDI_irrigation(1) and read_ssdi_input (deletion in
   !! Task 5).
   subroutine apply_irrigation_ssdi(ssdi)
      use, intrinsic :: iso_fortran_env, only: real64
      use irrigation_config_mod, only: irrigation_ssdi_t
      use error_mod, only: fatalerr_collected
      use variables, only: NumNod, zbotcp, tstart, tend, mairg,           &
                           pathwork,                                      &
                           ! shared / mode-0 globals
                           ifnd, nirri, qssdi, dt_SSDI_event,             &
                           nod_ssdi_irr, ssdi_date_irr, ssdi_rate_f_irr,  &
                           ssdi_amount_f_irr,                             &
                           ! mode-1 globals
                           ssdi_sched_type_irr, ssdi_threshold_irr,       &
                           ssdi_threshold_z_irr, ssdi_amount_irr,         &
                           ssdi_appl_rate_irr, sw_interval_irr,           &
                           days_interval_irr, days_counter_irr
      type(irrigation_ssdi_t), intent(in) :: ssdi

      integer  :: i, j, nod_top, nod_bot, ncomp
      real(real64) :: spread

      ! Resolve ssdi_z(1:2) -> layer indices via zbotcp walk.
      ! Mirrors src/crop/irrigation.f90 SSDI_irrigation(1) lines ~414-422.
      do j = 1, 2
         i = 1
         do while (zbotcp(i) > (ssdi%ssdi_z(j) + 1.0e-5_real64))
            i = i + 1
            if (i > NumNod) exit
         end do
         nod_ssdi_irr(j) = i
      end do
      nod_top = nod_ssdi_irr(1)
      nod_bot = nod_ssdi_irr(2)
      ncomp   = nod_bot - nod_top + 1
      if (ncomp < 1) then
         call fatalerr_collected('apply_irrigation_ssdi', &
                                 'ssdi_z resolves to zero compartments — check ssdi_z vs grid')
      end if

      qssdi = 0.0_real64

      select case (ssdi%schedule)
      case (0)
         call apply_ssdi_mode0(ssdi, ncomp)
      case (1)
         call apply_ssdi_mode1(ssdi, ncomp)
      end select
   end subroutine apply_irrigation_ssdi


   !> Mode-0 (fixed-date): stage CSV; populate ssdi_*_f_irr; deferred
   !! date-window validation; initial nirri entry-point from tstart.
   subroutine apply_ssdi_mode0(ssdi, ncomp)
      use, intrinsic :: iso_fortran_env, only: real64
      use csv_reader_mod, only: read_csv_table
      use error_mod, only: error_collection_t, fatalerr_collected
      use irrigation_config_mod, only: irrigation_ssdi_t
      use variables, only: pathwork, ifnd, nirri, tstart, tend, mairg,        &
                           ssdi_date_irr, ssdi_rate_f_irr, ssdi_amount_f_irr
      type(irrigation_ssdi_t), intent(in) :: ssdi
      integer,                 intent(in) :: ncomp

      real(real64), allocatable :: tbl(:,:)
      type(error_collection_t)  :: errs
      character(len=300) :: csvpath
      character(len=8)   :: hdr(3)
      integer :: i, n, nirri_init
      logical :: any_in_window, window_in_dates

      hdr(1) = 'date    '
      hdr(2) = 'rate_f  '
      hdr(3) = 'amount_f'
      csvpath = trim(pathwork) // trim(ssdi%fixed%events_file)
      call read_csv_table(trim(csvpath), hdr, tbl, errs)
      call errs%abort_if_fatal()

      n = 0
      if (allocated(tbl)) n = size(tbl, 1)
      if (n < 1) then
         call fatalerr_collected('apply_irrigation_ssdi', &
                                 'mode 0: events CSV has no rows')
         return
      end if
      if (n > mairg) then
         call fatalerr_collected('apply_irrigation_ssdi', &
                                 'mode 0: events CSV exceeds mairg rows')
         return
      end if

      ! Allocate legacy globals (sized to mairg, the legacy convention).
      if (allocated(ssdi_date_irr))     deallocate(ssdi_date_irr)
      if (allocated(ssdi_rate_f_irr))   deallocate(ssdi_rate_f_irr)
      if (allocated(ssdi_amount_f_irr)) deallocate(ssdi_amount_f_irr)
      allocate(ssdi_date_irr(mairg))
      allocate(ssdi_rate_f_irr(mairg))
      allocate(ssdi_amount_f_irr(mairg))
      ssdi_date_irr     = 0.0_real64
      ssdi_rate_f_irr   = 0.0_real64
      ssdi_amount_f_irr = 0.0_real64

      ! Copy + ascending check + unit conversions + compartment spread.
      do i = 1, n
         ssdi_date_irr(i) = tbl(i, 1)
         if (i > 1 .and. ssdi_date_irr(i) <= ssdi_date_irr(i-1)) then
            call fatalerr_collected('apply_irrigation_ssdi', &
                                    'mode 0: ssdi_date not strictly ascending')
            return
         end if
         ! mm/h -> cm/d
         ssdi_rate_f_irr(i)   = tbl(i, 2) * 0.1_real64 * 24.0_real64
         ! mm -> cm, then spread over `ncomp` compartments
         ssdi_amount_f_irr(i) = (tbl(i, 3) * 0.1_real64) / real(ncomp, real64)
      end do
      ifnd = n

      ! Date-window check (replaces the deleted checkdate). Pass iff
      ! at least one date in [tstart, tend], OR [tstart, tend] is
      ! contained in [date(1), date(n)].
      any_in_window  = .false.
      do i = 1, n
         if (ssdi_date_irr(i) >= tstart - 1.0e-6_real64 .and. &
             ssdi_date_irr(i) <= tend   + 1.0e-6_real64) then
            any_in_window = .true.
            exit
         end if
      end do
      window_in_dates = (ssdi_date_irr(1) <= tstart + 1.0e-6_real64 .and. &
                        ssdi_date_irr(n) >= tend   - 1.0e-6_real64)
      if (.not. any_in_window .and. .not. window_in_dates) then
         call fatalerr_collected('apply_irrigation_ssdi', &
                                 'mode 0: no ssdi_date within simulation period')
      end if

      ! Initial entry-point nirri based on tstart.
      nirri_init = 1
      do i = 1, n - 1
         if (tstart >= ssdi_date_irr(i)) nirri_init = i
      end do
      if (tstart >= ssdi_date_irr(n)) nirri_init = n
      nirri = nirri_init
   end subroutine apply_ssdi_mode0


   !> Mode-1 (scheduled-trigger): copy scheduled sub-block to legacy
   !! globals; set nirri/dt_SSDI_event/days_counter defaults
   !! (preserves the d8a88d6 regression-fix invariant for
   !! dt_SSDI_event = 1.0).
   subroutine apply_ssdi_mode1(ssdi, ncomp)
      use, intrinsic :: iso_fortran_env, only: real64
      use irrigation_config_mod, only: irrigation_ssdi_t
      use variables, only: nirri, dt_SSDI_event,                          &
                           ssdi_sched_type_irr, ssdi_threshold_irr,       &
                           ssdi_threshold_z_irr, ssdi_amount_irr,         &
                           ssdi_appl_rate_irr, sw_interval_irr,           &
                           days_interval_irr, days_counter_irr
      type(irrigation_ssdi_t), intent(in) :: ssdi
      integer,                 intent(in) :: ncomp

      ssdi_sched_type_irr  = ssdi%scheduled%sched_type
      ssdi_threshold_irr   = ssdi%scheduled%threshold
      ssdi_threshold_z_irr = ssdi%scheduled%threshold_depth

      ! Unit conversions: mm -> cm, mm/h -> cm/d, then spread.
      ssdi_amount_irr      = (ssdi%scheduled%ssdi_amount * 0.1_real64) / real(ncomp, real64)
      ssdi_appl_rate_irr   = ssdi%scheduled%ssdi_appl_rate * 0.1_real64 * 24.0_real64

      sw_interval_irr      = ssdi%scheduled%sw_interval
      if (ssdi%scheduled%sw_interval == 0) then
         days_interval_irr = 1
      else
         days_interval_irr = ssdi%scheduled%days_interval
      end if
      days_counter_irr     = 366

      nirri          = 1
      dt_SSDI_event  = 1.0_real64
   end subroutine apply_ssdi_mode1
```

Wire the call. Find the existing `flSSDI = (config%irrigation%swssdi == 1)` line in `config_to_variables` (around line 485 — line numbers may shift). Right after that line, add:
```fortran
      if (flSSDI) call apply_irrigation_ssdi(config%irrigation%ssdi)
```

Note on global names: the actual `variables` module may use `nod_ssdi(2)` (no `_irr` suffix) in some places — check `src/crop/irrigation.f90:343-349` for the exact mapping. In that file the locals are reset *from* the `_irr` versions (`nod_ssdi = nod_ssdi_irr`). The adapter writes the `_irr` names; the runtime later reads them into the unsuffixed locals. **Do not rename — verify by grep what the actual storage names are and write to those.**

- [ ] **Step 5: Build and run tests + check-full**

```
pixi run -e test build-linux
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:|tests,"
pixi run -e test check-full 2>&1 | grep -E "Results:"
```
Expected: build clean; pFUnit `Ok: 1, Fail: 0`; check-full `5 passed, 0 failed` (regression cases all have `swssdi=0`, the new adapter call short-circuits — they're unchanged).

- [ ] **Step 6: Commit**

```bash
git add src/io/toml/config_to_variables.f90 tests/unit/io/toml/test_apply_irrigation_ssdi.pf tests/unit/io/toml/fixtures/ssdi_events_small.csv
git commit -m "$(cat <<'EOF'
feat(io): apply_irrigation_ssdi populates legacy globals from typed config

Adds the apply_irrigation_ssdi helper to config_to_variables_mod,
called when flSSDI = .true.. Resolves ssdi_z(1:2) to nod_ssdi(1:2)
layer indices via the zbotcp walk; mode 0 reads the events CSV
(via read_csv_table per ADR 0012; column 1 auto-converts ISO dates
to days-since-1900); mode 1 copies the scheduled sub-block; both
modes apply legacy unit conversions (mm/h -> cm/d, mm -> cm, then
spread over the in-zone compartments) and initialize qssdi = 0.

The deferred date-window check (at least one date in [tstart,tend]
OR window contained in [date(1), date(ifnd)]) replaces the
to-be-deleted checkdate.f90 — same logic, expressed against
fatalerr_collected.

Mode 1 sets dt_SSDI_event = 1.0 (preserves the d8a88d6 regression-fix
invariant) and the entry-point nirri = 1.

Replaces SSDI_irrigation(1) and read_ssdi_input (deletion in next
commit). check-full regression suite unchanged (all five cases have
swssdi=0; adapter short-circuits).

Part of SSDI TOML port (spec
docs/superpowers/specs/2026-05-07-ssdi-toml-port-design.md).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: Retire `read_ssdi_input` + collapse `SSDI_irrigation(1)`

**Files:**
- Modify: `src/crop/irrigation.f90` (delete `subroutine read_ssdi_input`, collapse `SSDI_irrigation(1)` body, drop unused locals/imports)

- [ ] **Step 1: Verify the adapter is doing the work**

Confirm the global state populated by `apply_irrigation_ssdi` covers everything `SSDI_irrigation(1)` did:
- ssdi_z resolution → `nod_ssdi_irr(1:2)` ✓
- CSV staging → `ssdi_date_irr`, `ssdi_rate_f_irr`, `ssdi_amount_f_irr`, `ifnd` (mode 0) ✓
- Scheduled-block copy → `ssdi_sched_type_irr`, `ssdi_threshold_irr`, etc. (mode 1) ✓
- `nirri` initial entry-point ✓
- `dt_SSDI_event = 1.0` ✓
- `qssdi = 0.0` ✓
- Date-window check ✓

Anything that **stays** in `SSDI_irrigation(1)`:
- Reading the `_irr` snapshots into the per-day working locals at the top of `SSDI_irrigation` (lines ~343-349) — those are needed every day, not just at init.

So `case (1)` body collapses; the subroutine itself stays (still has cases 2 and 9).

- [ ] **Step 2: Delete `subroutine read_ssdi_input`**

In `src/crop/irrigation.f90`, find by `grep -n "subroutine read_ssdi_input\|end subroutine read_ssdi_input" src/crop/irrigation.f90`. Delete the entire body inclusive of the header comment line above it.

- [ ] **Step 3: Collapse `SSDI_irrigation` `case (1)` to a single comment + return**

Find by `grep -n "case (1)" src/crop/irrigation.f90` (in the SSDI_irrigation subroutine — there's another `case (1)` in `irrigation()`). Replace the entire case-1 block (currently the swp open + rdinit + rdsinr + rdscha + read_ssdi_input call + nirri loop + nod_ssdi resolution + amount spread) with:

```fortran
   case (1)
      ! [irrigation.ssdi] init was performed at config-load time by
      ! apply_irrigation_ssdi (config_to_variables.f90). Per ADR 0022,
      ! this case is now a no-op; the runtime reads the staged
      ! _irr snapshots into per-day locals at the top of
      ! SSDI_irrigation (above the select case).
      return
```

- [ ] **Step 4: Drop unused locals + imports**

After deletion, run:
```bash
grep -nE "swpfile|logf\b|ssdi_file|RDinit|getun2|rdscha|rdsinr|rdinqr" src/crop/irrigation.f90
```
Expected: zero matches inside the SSDI_irrigation subroutine. If any local declarations (`integer :: swp`, `character(len=...) :: ssdi_file`) became unused, drop them. Trim the `use variables, only: ...` line at the top of the SSDI_irrigation subroutine to drop `swpfile`, `logf`, `ssdi_file` if they're no longer referenced.

Also remove the `if (swssdi == 0) return` self-check at `SSDI_irrigation` case(2) (line ~634; per ADR 0020 commit's "harmless defense-in-depth" note — `flSSDI` gating is now the single source of truth).

- [ ] **Step 5: Build and run**

```
pixi run -e test build-linux
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:|tests,"
pixi run -e test check-full 2>&1 | grep -E "Results:"
echo "=== acceptance greps ==="
grep -n "subroutine read_ssdi_input" src/ || echo "OK: no read_ssdi_input"
grep -n "swpfile" src/crop/irrigation.f90 || echo "OK: no swpfile in irrigation.f90"
```
Expected: build clean; pFUnit unchanged; check-full 5/5; both greps return "OK: ...".

Note: `swpfile` may still appear in `src/io/toml/config_to_variables.f90` (the `logf`-opening HACK) — that cleanup is out of scope for this spec; deferred to a follow-on commit per the spec's "Out of scope" section. Verify via `grep -rn "swpfile" src/` returning only the config_to_variables.f90 hits.

- [ ] **Step 6: Commit**

```bash
git add src/crop/irrigation.f90
git commit -m "$(cat <<'EOF'
refactor(crop): retire read_ssdi_input; collapse SSDI_irrigation(1)

Removes the last TTutil-based reads in irrigation.f90:
- Deletes subroutine read_ssdi_input (was the swap.swp +
  ssdi_file reader).
- Collapses SSDI_irrigation(1) body to a return — the case-1
  init responsibilities (CSV staging, ssdi_z -> nod_ssdi
  resolution, nirri entry-point, dt_SSDI_event default,
  qssdi=0) all moved to apply_irrigation_ssdi (previous commit,
  config_to_variables.f90).
- Drops swpfile, logf, ssdi_file from `use variables` and the
  associated locals.
- Removes the harmless `if (swssdi == 0) return` self-check at
  case(2) — flSSDI gating is now the single source of truth.

apply_irrigation_ssdi is the sole SSDI-init path; runs at
config-load time. SSDI_irrigation cases 2 and 9 (per-day apply
and end-of-event reset) stay unchanged.

Acceptance:
- `grep -n "subroutine read_ssdi_input" src/` -> no matches.
- `grep -n "swpfile" src/crop/irrigation.f90` -> no matches.
- check-full 5/5 (all five regression cases have swssdi=0; the
  adapter call short-circuits).

The remaining swpfile references in src/ are confined to
src/io/toml/config_to_variables.f90 (the logf-opening HACK,
out of scope here; ADR 0022 follow-up).

Part of SSDI TOML port (spec
docs/superpowers/specs/2026-05-07-ssdi-toml-port-design.md).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: Retire `src/io/checkdate.f90`

**Files:**
- Delete: `src/io/checkdate.f90`
- Modify: `meson.build` (drop production source)
- Modify: `tests/unit/meson.build` (drop test source)

- [ ] **Step 1: Verify zero callers remain**

```bash
grep -rn "\bcheckdate\b" src/ tests/ --include="*.f90" --include="*.F90" --include="*.pf" 2>/dev/null
```
Expected: only matches inside `src/io/checkdate.f90` itself (the subroutine definition and its file-header comments). Confirm there are no other callers — Task 5 deleted the only call site (in `read_ssdi_input`).

If any caller appears outside `src/io/checkdate.f90`, **stop** — that's a regression that needs investigation. Check `read_ssdi_input` was actually deleted in Task 5 and the date-window check was relocated to `apply_irrigation_ssdi` in Task 4.

- [ ] **Step 2: Drop from meson sources**

In `meson.build` (production sources), remove this line:
```meson
    'src/io/checkdate.f90',
```

In `tests/unit/meson.build` `pfunit_extra_sources`, remove this line:
```meson
        '../../src/io/checkdate.f90',
```

- [ ] **Step 3: Delete the file**

```bash
git rm src/io/checkdate.f90
```

- [ ] **Step 4: Build and run**

```
pixi run -e test build-linux
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:|tests,"
pixi run -e test check-full 2>&1 | grep -E "Results:"
echo "=== final greps ==="
grep -rn "checkdate" src/ tests/ 2>/dev/null && echo "FAIL" || echo "OK: no checkdate"
ls src/io/checkdate.f90 2>/dev/null && echo "FAIL" || echo "OK: file gone"
```
Expected: build clean; pFUnit unchanged; check-full 5/5; both greps return "OK: ...".

- [ ] **Step 5: Commit**

```bash
git add meson.build tests/unit/meson.build
git commit -m "$(cat <<'EOF'
refactor(io): retire src/io/checkdate.f90 — no remaining callers

Deletes the date-range validator that was extracted from the
deleted readswap.f90 in SS-C step 5 (commit fa1aaec). It existed
solely to support read_ssdi_input's date-window check; with that
caller deleted in the previous commit, checkdate has zero callers.

The equivalent date-window logic now lives in
apply_irrigation_ssdi (config_to_variables.f90), expressed against
fatalerr_collected — same behaviour, same ±1.d-6 tolerance, same
"any date in window OR window contained in dates" semantics.

Acceptance:
- `grep -rn "checkdate" src/ tests/` -> no matches.
- src/io/checkdate.f90 does not exist.
- check-full 5/5.

Closes ADR 0019's "two stub-readers" footnote: with ADR 0021
(tillage) and ADR 0022 (SSDI) landed, the only remaining TTutil
runtime calls in src/ are utility (rdsets/rdfrom in swap_main.f90
+ rddtmp in swapoutput.f90) — not data readers.

Part of SSDI TOML port (spec
docs/superpowers/specs/2026-05-07-ssdi-toml-port-design.md).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: ADR 0022 + schema doc

**Files:**
- Create: `docs/adr/0022-ssdi-toml-port.md`
- Modify: `docs/adr/index.md` (add ADR 0022 row)
- Modify: `docs/configuration-schema.md` (document `[irrigation.ssdi]` + CSV format)

- [ ] **Step 1: Write `docs/adr/0022-ssdi-toml-port.md`**

```markdown
---
title: "ADR 0022 — SSDI parameters ported to [irrigation.ssdi] TOML block"
date: 2026-05-07
status: accepted
---

# ADR 0022: SSDI parameters ported to `[irrigation.ssdi]` TOML block

## Context

ADR 0019's "Update 2026-05-06" closed the umbrella physical
deletion of the legacy fixed-format reader, but two TTutil-based
reads survived as foot-noted exceptions: `Read_Tillage` (closed
by ADR 0021) and `read_ssdi_input` (closed by this ADR).

This ADR ports the SSDI subsystem and retires the standalone
`src/io/checkdate.f90` along with it.

## Decision

SSDI parameters move to a new `[irrigation.ssdi]` sub-table with
an explicit `schedule = 0|1` discriminator. Mode-specific
parameters live in disjoint sub-tables:

- `[irrigation.ssdi.fixed]` (`schedule = 0`) — the date-rate-amount
  table is held in a CSV companion (`events_file = "..."`); CSV-only,
  no inline path. Tables can hit 366 rows and inline TOML becomes
  unreadable at that scale.
- `[irrigation.ssdi.scheduled]` (`schedule = 1`) — trigger-based
  parameters (`sched_type`, thresholds, amount, application rate,
  interval gating) inline.

`swssdi` stays at `[irrigation].swssdi` (least disruption to
existing TOML cases; parallel to the soil.tillage / soil.swtill
pattern from ADR 0021).

`read_ssdi_input` is deleted. A new helper
`apply_irrigation_ssdi` in `config_to_variables.f90` populates
the legacy `_irr`-suffixed globals (`ssdi_date_irr`,
`ssdi_rate_f_irr`, `ssdi_amount_f_irr`, `ssdi_sched_type_irr`,
`ssdi_threshold_irr`, `ssdi_threshold_z_irr`, `ssdi_amount_irr`,
`ssdi_appl_rate_irr`, `sw_interval_irr`, `days_interval_irr`,
`nod_ssdi_irr`, plus `ifnd`, `nirri`, `dt_SSDI_event`, `qssdi`)
at config-load time.

`src/io/checkdate.f90` is deleted; the date-window check lives
in `apply_irrigation_ssdi`'s deferred-validation tail, expressed
against `fatalerr_collected` instead of `fatalerr` + `STOP`.

## Schema

```toml
[irrigation]
swssdi = 1

[irrigation.ssdi]
schedule = 0          # 0 fixed-date, 1 scheduled-trigger
ssdi_z   = [-30.0, -50.0]   # cm; both equal for single-depth

# When schedule = 0:
[irrigation.ssdi.fixed]
events_file = "ssdi_events.csv"

# When schedule = 1:
[irrigation.ssdi.scheduled]
sched_type      = 1   # 1=Tred, 2=presh, 3=watc
threshold       = 0.7
threshold_depth = -50.0   # cm; required when sched_type > 1
ssdi_amount     = 10.0    # mm
ssdi_appl_rate  = 5.0     # mm/h
sw_interval     = 0       # 0|1
# days_interval required only when sw_interval = 1
```

CSV companion (`ssdi_events.csv`):
```
date,rate_f,amount_f
2003-04-15,5.0,10.0
2003-05-02,5.0,8.0
```

When `swssdi = 0`, the entire `[irrigation.ssdi]` block can be
omitted.

## Consequences

- `swpfile` global pointer no longer read by `irrigation.f90`. The
  remaining `swpfile`/`logf` plumbing in
  `src/io/toml/config_to_variables.f90` becomes orphan; cleanup
  is the next step (separate, smaller follow-on commit).
- `read_ssdi_input` deleted; `apply_irrigation_ssdi` is the sole
  SSDI-init path.
- `src/io/checkdate.f90` deleted; `fatalerr_collected` is the
  sole validation channel for the date-window check.
- With ADR 0021 + ADR 0022 landed, the only TTutil runtime calls
  remaining in `src/` are utility (`rdsets`/`rdfrom` in
  `swap_main.f90` + `rddtmp` in `swapoutput.f90`) — not data
  readers.

## Tests

pFUnit suites under `tests/unit/io/toml/`:
- `test_read_irrigation_ssdi_toml`
- `test_irrigation_ssdi_validate`
- `test_apply_irrigation_ssdi`

CSV fixture: `tests/unit/io/toml/fixtures/ssdi_events_small.csv`.

No `swssdi = 1` regression case is added; deferred until a
known-good legacy comparator is available.

## Acceptance

- `grep -n "subroutine read_ssdi_input" src/` -> no matches.
- `grep -n "swpfile" src/crop/irrigation.f90` -> no matches.
- `grep -rn "checkdate" src/` -> no matches.
- `src/io/checkdate.f90` does not exist.
- `pixi run -e test test-pfunit` -> `Ok: 1, Fail: 0`.
- `pixi run -e test check-full` -> `5 passed, 0 failed`.
```

- [ ] **Step 2: Add `[irrigation.ssdi]` to `docs/configuration-schema.md`**

Open `docs/configuration-schema.md`, find the existing `[irrigation]` section. Append a new `### [irrigation.ssdi]` sub-section using the same shape as the surrounding sub-sections. Cover:

- **Activation:** present only when `[irrigation].swssdi = 1`.
- **Top-level fields:**
  - `schedule` (integer, 0..1, default 0) — mode discriminator.
  - `ssdi_z` (real array of 1 or 2 elements, each in [-100.0, 0.0] cm) — application depth or interval. For single-depth, supply one element.
- **`[irrigation.ssdi.fixed]`** (required when `schedule = 0`):
  - `events_file` (string) — relative path to events CSV.
- **`[irrigation.ssdi.scheduled]`** (required when `schedule = 1`):
  - `sched_type` (integer, 1..3) — 1=Tred, 2=presh, 3=watc.
  - `threshold` (real) — range depends on `sched_type` (1 → [0,1]; 2 → [-1e7, 0]; 3 → [0,1]).
  - `threshold_depth` (real, [-100, 0] cm) — required when `sched_type > 1`.
  - `ssdi_amount` (real, [0, 100] mm).
  - `ssdi_appl_rate` (real, [0, 100] mm/h).
  - `sw_interval` (integer, 0..1).
  - `days_interval` (integer, 1..366) — required when `sw_interval = 1`.
- **CSV format** (`events_file`):
  - Columns: `date` (ISO YYYY-MM-DD), `rate_f` (mm/h), `amount_f` (mm).
  - Up to 366 rows; dates strictly ascending; at least one row in `[start_date, end_date]` OR window contained in date range.

Match the existing sub-section formatting (heading depth, table style, prose density).

- [ ] **Step 3: Update `docs/adr/index.md`**

Append a row for ADR 0022 — match the existing format. Likely:
```markdown
| [0022](0022-ssdi-toml-port.md) | SSDI parameters ported to `[irrigation.ssdi]` TOML block; checkdate.f90 retired | 2026-05-07 | accepted |
```

- [ ] **Step 4: Final verification**

```
pixi run -e test build-linux
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:|tests,"
pixi run -e test check-full 2>&1 | grep -E "Results:"
echo "=== acceptance greps ==="
grep -n "subroutine read_ssdi_input" src/ || echo "OK: no read_ssdi_input"
grep -n "swpfile" src/crop/irrigation.f90 || echo "OK: no swpfile in irrigation.f90"
grep -rn "checkdate" src/ 2>/dev/null && echo "FAIL: checkdate refs remain" || echo "OK: no checkdate"
ls src/io/checkdate.f90 2>/dev/null && echo "FAIL: file still exists" || echo "OK: file deleted"
```
Expected: all green; greps return "OK: ...".

- [ ] **Step 5: Commit**

```bash
git add docs/adr/0022-ssdi-toml-port.md docs/configuration-schema.md docs/adr/index.md
git commit -m "$(cat <<'EOF'
docs: ADR 0022 — SSDI parameters ported to [irrigation.ssdi]

Captures the decision to port the SSDI subsystem from TTutil
swap.swp + ssdi_file reads to a typed [irrigation.ssdi] sub-table
with an explicit schedule = 0|1 discriminator; CSV companion for
fixed-mode events; src/io/checkdate.f90 retired.

Closes the second half of ADR 0019's "two stub-readers also
survive" footnote (ADR 0021 covered tillage). With ADR 0021 +
ADR 0022 landed, the only remaining TTutil runtime calls in src/
are utility (rdsets/rdfrom + rddtmp), not data readers.

- docs/adr/0022-ssdi-toml-port.md: new ADR.
- docs/configuration-schema.md: documents [irrigation.ssdi] both
  modes + the events-CSV format.
- docs/adr/index.md: ADR 0022 row.

Closes the SSDI TOML port spec
(docs/superpowers/specs/2026-05-07-ssdi-toml-port-design.md).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Self-Review

**Spec coverage:**

| Spec section | Task |
|---|---|
| §Schema (`irrigation_ssdi_*` types) | Task 1 |
| §Validation rules (per-block) | Task 2 |
| §TOML reader (`read_irrigation_ssdi_toml`) | Task 3 |
| §CSV reader (used inside adapter) | Task 4 |
| §Adapter (`apply_irrigation_ssdi`, both modes, deferred date-window check) | Task 4 |
| §Runtime change (delete `read_ssdi_input`, collapse `SSDI_irrigation(1)`) | Task 5 |
| §Retire `checkdate.f90` | Task 6 |
| §Tests (3 pFUnit suites + CSV fixture) | Tasks 1, 2, 3, 4 (one suite + the fixture per task) |
| §ADR 0022 + schema doc | Task 7 |

All spec sections covered.

**Acceptance criteria** (from spec):
- `pixi run -e test build-linux` clean — verified each task.
- `pixi run -e test test-pfunit` — verified each task.
- `pixi run -e test check-full` — verified Tasks 4-7.
- `grep -n "subroutine read_ssdi_input" src/` — verified Tasks 5, 7.
- `grep -rn "swpfile" src/` — verified Task 5 (note: the
  config_to_variables.f90 plumbing is intentionally out of scope per
  spec §Non-goals).
- `grep -rn "checkdate" src/` — verified Task 6.
- `src/io/checkdate.f90` does not exist — verified Task 6.
- ADR 0022 committed — Task 7.
- `docs/configuration-schema.md` updated — Task 7.

**Type / signature consistency:**
- `irrigation_ssdi_fixed_t` field `events_file` — used identically in Tasks 1, 2, 3, 4.
- `irrigation_ssdi_scheduled_t` fields (`sched_type`, `threshold`, `threshold_depth`, `ssdi_amount`, `ssdi_appl_rate`, `sw_interval`, `days_interval`) — consistent.
- `irrigation_ssdi_t` fields (`schedule`, `ssdi_z`, `fixed`, `scheduled`) — consistent.
- `apply_irrigation_ssdi(ssdi)` signature in Task 4 matches the call site.
- `apply_ssdi_mode0(ssdi, ncomp)` and `apply_ssdi_mode1(ssdi, ncomp)` private helpers introduced in Task 4 only.

**Open questions punted to implementation:**
- Exact `_irr` suffix list — Task 4 Step 2 instructs the implementer to verify by `grep -n "_irr\b" src/core/variables.f90 | head -30` before writing the tests.
- Whether `path_helpers_mod` (or the bare `pathwork` pattern from existing CSV adapters) is the right path-resolution helper — the existing `fixed_events_file` block in `config_to_variables.f90:184-220` uses `trim(pathwork) // trim(filename)`, so Task 4 follows that.
- The CSV reader's exact `expected_header` argument format — `read_csv_table(path, expected_header(:), table, errors)` where `expected_header` is `character(len=*)` and column 1 being `'date'` triggers ISO-date auto-conversion. Padding strings to a uniform length is a Fortran idiom; Task 4 uses 8-char headers. Implementer can adjust if the CSV reader is stricter.

These are pickable at execution time without re-planning.

---

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-05-07-ssdi-toml-port.md`.
