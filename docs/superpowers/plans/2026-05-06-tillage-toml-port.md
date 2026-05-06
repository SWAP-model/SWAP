# Tillage TOML port — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Port the tillage subsystem from TTutil-based `swap.swp` reads to a typed `[soil.tillage]` TOML block; retire `Read_Tillage`; remove `tillage.f90`'s dependence on `swpfile`.

**Architecture:** Add three nested types to `soil_config_t` (`soil_tillage_event_t`, `soil_tillage_type_t`, `soil_tillage_t`) holding the parameters previously read by `Read_Tillage`. Extend the TOML reader and validator. Add an `apply_soil_tillage` helper to `config_to_variables.f90` that allocates and populates the legacy globals (`Date_tillage`, `Z_tillage`, `Type_tillage`, `iType_Tillage`, `TAB_Rho_cons`, `TAB_Rho_tillage`, `TAB_K_R_cons`, `TAB_Rho_match`, `TAB_N_match`, `Ntill`, `Ntypes`, `Max_Z_tillage`, `iTT1`, `iTT2`) at config-load time. Delete `subroutine Read_Tillage` and the `call Read_Tillage` line in `DoTillage(1)`.

**Tech Stack:** Fortran 2008, meson + ninja build, toml-f for parsing, pFUnit for unit tests, gfortran. Reference patterns: `src/io/toml/read_drainage_toml.f90` (array-of-tables), `src/io/toml/read_irrigation_toml.f90` (date parsing via `parse_date_to_days1900(toml_datetime)`), `src/io/toml/config_to_variables.f90:482-485` (existing `flTillage` setup).

**Spec:** `docs/superpowers/specs/2026-05-06-tillage-toml-port-design.md`.

---

## File Structure

**Created:**

| Path | Responsibility |
|---|---|
| `src/io/toml/read_soil_tillage_toml.f90` | New module `read_soil_tillage_toml_mod` that parses `[soil.tillage]` (scalars + `[[events]]` + `[[types]]` arrays-of-tables) into `soil_tillage_t`. Sibling to `read_soil_toml.f90` to keep that file from growing. |
| `tests/unit/io/toml/test_read_soil_tillage_toml.pf` | pFUnit suite asserting schema parses correctly: scalars, events array, types array, optional `rho_match`/`N_match` when `i_n_model=3`, missing-block default. |
| `tests/unit/io/toml/test_soil_tillage_validate.pf` | pFUnit suite asserting every validator rule fires correctly. |
| `tests/unit/io/toml/test_apply_soil_tillage.pf` | pFUnit suite asserting `apply_soil_tillage` populates the legacy globals correctly. |
| `docs/adr/0021-tillage-toml-port.md` | ADR capturing the decision. |

**Modified:**

| Path | Change |
|---|---|
| `src/config/soil_config.f90` | Add `soil_tillage_event_t`, `soil_tillage_type_t`, `soil_tillage_t` type definitions; add `tillage` field on `soil_config_t`; extend `soil_config_validate` with the new rules. |
| `src/io/toml/read_soil_toml.f90` | After parsing the `[soil]` table, call `read_soil_tillage_toml` to populate `config%soil%tillage`. |
| `src/io/toml/config_to_variables.f90` | Add `apply_soil_tillage(tillage)` helper; call it from `config_to_variables` when `flTillage` is true. |
| `src/crop/tillage.f90` | Delete `subroutine Read_Tillage` (lines 423-508) and the `call Read_Tillage` line in `DoTillage(1)` (line 64). Drop `swpfile` and `logf` from imports if no longer used. Optionally drop the `if (swtill /= 1) return` self-checks at lines 53/72 (defense-in-depth, redundant under `flTillage` gating). |
| `meson.build` | Add `'src/io/toml/read_soil_tillage_toml.f90'` to production sources. |
| `tests/unit/meson.build` | Add `'../../src/io/toml/read_soil_tillage_toml.f90'` to `pfunit_extra_sources`; add the three new `.pf` files to `pf_files`. |
| `docs/configuration-schema.md` | Document the `[soil.tillage]` block. |

---

## Task 1: Schema types

**Files:**
- Modify: `src/config/soil_config.f90` (add type definitions before `soil_config_t`; add `tillage` field on `soil_config_t`)
- Test: `tests/unit/io/toml/test_apply_soil_tillage.pf` (default-init test added here; the suite will grow with adapter tests in Task 4)
- Modify: `tests/unit/meson.build` (register the new `.pf` file)

- [ ] **Step 1: Write the failing test**

Create `tests/unit/io/toml/test_apply_soil_tillage.pf`:

```fortran
@test
subroutine test_soil_tillage_default_init()
   use soil_config_mod, only: soil_config_t
   use funit, only: assertEqual, assertFalse, assertEqual
   implicit none
   type(soil_config_t) :: cfg

   ! Default-init values for the new tillage sub-type.
   call assertEqual(2, cfg%tillage%i_n_model, 'i_n_model default')
   call assertEqual(2, cfg%tillage%iRedist,   'iRedist default')
   @assertFalse(allocated(cfg%tillage%events), 'events not allocated by default')
   @assertFalse(allocated(cfg%tillage%types),  'types not allocated by default')
end subroutine
```

- [ ] **Step 2: Wire the test file into meson and run to verify it fails**

Add to `tests/unit/meson.build`'s `pf_files` list:

```meson
        'io/toml/test_apply_soil_tillage.pf',
```

Run:
```
pixi run -e test build-linux
```
Expected: build fails — `cfg%tillage` is not a member of `soil_config_t`.

- [ ] **Step 3: Add the new types in `src/config/soil_config.f90`**

Insert the three types **before** the `type :: soil_config_t` declaration at line 79 (immediately after `end type soil_initial_t`):

```fortran
   !> Single tillage event row read from `[[soil.tillage.events]]`.
   type :: soil_tillage_event_t
      character(len=10) :: date      = ''       !! ISO YYYY-MM-DD (parsed to days-since-1900 by adapter)
      real(real64)      :: z         = 0.0_real64
      real(real64)      :: intensity = 0.0_real64
      integer           :: type_id   = 0
   end type soil_tillage_event_t

   !> Single tillage type row read from `[[soil.tillage.types]]`.
   type :: soil_tillage_type_t
      integer      :: id          = 0
      real(real64) :: rho_cons    = 0.0_real64
      real(real64) :: rho_tillage = 0.0_real64
      real(real64) :: k_R         = 0.0_real64
      real(real64) :: rho_match   = -99.0_real64   !! used only when i_n_model = 3
      real(real64) :: N_match     = -99.0_real64   !! used only when i_n_model = 3
   end type soil_tillage_type_t

   !> `[soil.tillage]` block container.
   type :: soil_tillage_t
      integer :: i_n_model = 2
      integer :: iRedist   = 2
      type(soil_tillage_event_t), allocatable :: events(:)
      type(soil_tillage_type_t),  allocatable :: types(:)
   end type soil_tillage_t
```

Add the field on `soil_config_t` (alongside `discretization`, `frost`, `hydraulics`, `initial`):

```fortran
      type(soil_tillage_t)        :: tillage
```

- [ ] **Step 4: Build and run the test to verify it passes**

```
pixi run -e test build-linux
pixi run -e test test-pfunit
```
Expected: build clean; pFUnit `Ok: 1, Fail: 0`.

- [ ] **Step 5: Commit**

```bash
git add src/config/soil_config.f90 tests/unit/io/toml/test_apply_soil_tillage.pf tests/unit/meson.build
git commit -m "$(cat <<'EOF'
feat(config): add soil_tillage_t types to soil_config_t

Adds soil_tillage_event_t, soil_tillage_type_t, and soil_tillage_t
containers to soil_config_mod, plus the `tillage` field on
soil_config_t. Defaults are off-state (no events/types allocated;
i_n_model=2; iRedist=2). pFUnit smoke test asserts default-init.

No reader, validator, or adapter wired yet — those land in
subsequent tasks.

Part of tillage TOML port (spec
docs/superpowers/specs/2026-05-06-tillage-toml-port-design.md).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Validator extension

**Files:**
- Modify: `src/config/soil_config.f90:134-192` (extend `soil_config_validate`)
- Test: `tests/unit/io/toml/test_soil_tillage_validate.pf` (new)
- Modify: `tests/unit/meson.build` (register the new `.pf` file)

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/io/toml/test_soil_tillage_validate.pf` with one test per rule. Sample skeleton (write all rules in one file):

```fortran
@test
subroutine test_validate_swtill1_requires_events()
   use soil_config_mod, only: soil_config_t, soil_config_validate
   use error_mod, only: error_collection_t
   use funit
   implicit none
   type(soil_config_t)     :: cfg
   type(error_collection_t) :: errs

   cfg%swtill = 1
   cfg%tillage%i_n_model = 2
   cfg%tillage%iRedist   = 2
   ! events left unallocated → expect error
   call soil_config_validate(cfg, errs)
   @assertTrue(errs%count() > 0, 'expected error for swtill=1 with no events')
end subroutine

@test
subroutine test_validate_event_type_id_must_match_a_type()
   use soil_config_mod, only: soil_config_t, soil_tillage_event_t, soil_tillage_type_t, soil_config_validate
   use error_mod, only: error_collection_t
   use funit
   implicit none
   type(soil_config_t)     :: cfg
   type(error_collection_t) :: errs

   cfg%swtill = 1
   allocate(cfg%tillage%events(1))
   cfg%tillage%events(1) = soil_tillage_event_t(date='2003-04-15', z=30.0d0, &
                                                intensity=1.0d0, type_id=2)
   allocate(cfg%tillage%types(1))
   cfg%tillage%types(1) = soil_tillage_type_t(id=1, rho_cons=1500.0d0, &
                                              rho_tillage=1100.0d0, k_R=0.05d0)
   call soil_config_validate(cfg, errs)
   @assertTrue(errs%count() > 0, 'expected error for events(1).type_id not in types')
end subroutine

@test
subroutine test_validate_dates_strictly_ascending()
   ! ... same shape, two events with date(2) <= date(1), expect error
end subroutine

@test
subroutine test_validate_intensity_in_range()
   ! ... events(1).intensity = 1.5, expect error
end subroutine

@test
subroutine test_validate_i_n_model_3_requires_rho_match_N_match()
   ! ... cfg%tillage%i_n_model = 3, types(1).rho_match left default (-99), expect error
end subroutine

@test
subroutine test_validate_swtill0_skips_tillage_block()
   ! ... cfg%swtill = 0, no events/types, expect zero new errors
end subroutine

@test
subroutine test_validate_clean_config_no_errors()
   ! ... cfg%swtill = 1, one event referencing one type, all in range, expect zero errors
end subroutine
```

(Write all seven concrete test subroutines; for brevity here only the first two are spelled out. Mirror them following the same shape.)

- [ ] **Step 2: Wire the test file and run to verify failures**

Add to `tests/unit/meson.build`'s `pf_files` list:
```meson
        'io/toml/test_soil_tillage_validate.pf',
```

Run:
```
pixi run -e test build-linux
pixi run -e test test-pfunit 2>&1 | grep -E "FAIL|Fail:"
```
Expected: build clean (the test fixtures compile against the existing types from Task 1); tests fail because the validator hasn't been extended yet — every "expected error" assertion fails because `errs%count() == 0`.

- [ ] **Step 3: Extend `soil_config_validate` in `src/config/soil_config.f90`**

Insert, just before `end subroutine soil_config_validate` at line 192, a new block guarded by `if (self%swtill == 1)`:

```fortran
      ! [soil.tillage] validation — only fires when swtill = 1.
      ! Depth check (events.z vs |zbotcp(NumNod)|) is deferred to the
      ! adapter (apply_soil_tillage) since NumNod isn't known here.
      if (self%swtill == 1) then
         call check_int_range(self%tillage%i_n_model, 1, 3, &
                              "soil.tillage.i_n_model", errors)
         call check_int_range(self%tillage%iRedist,   0, 2, &
                              "soil.tillage.iRedist",   errors)

         if (.not. allocated(self%tillage%events) .or. &
             size(self%tillage%events) == 0) then
            call errors%append(ERR_VALIDATION_EMPTY, &
                               'must be non-empty when soil.swtill = 1', &
                               'soil.tillage.events')
         end if
         if (.not. allocated(self%tillage%types) .or. &
             size(self%tillage%types) == 0) then
            call errors%append(ERR_VALIDATION_EMPTY, &
                               'must be non-empty when soil.swtill = 1', &
                               'soil.tillage.types')
         end if

         if (allocated(self%tillage%events)) then
            do i = 1, size(self%tillage%events)
               associate (ev => self%tillage%events(i))
                  call check_real_range(ev%intensity, 0.0_real64, 1.0_real64, &
                                        'soil.tillage.events.intensity', errors)
                  if (ev%type_id < 1) then
                     call errors%append(ERR_VALIDATION_RANGE, &
                                        'must be >= 1', &
                                        'soil.tillage.events.type_id')
                  end if
                  ! membership: type_id must appear in tillage.types(:).id
                  if (allocated(self%tillage%types)) then
                     if (.not. any(self%tillage%types%id == ev%type_id)) then
                        call errors%append(ERR_VALIDATION_RANGE, &
                                           'type_id not present in soil.tillage.types', &
                                           'soil.tillage.events.type_id')
                     end if
                  end if
                  ! strictly ascending dates
                  if (i > 1) then
                     if (ev%date <= self%tillage%events(i-1)%date) then
                        call errors%append(ERR_VALIDATION_RANGE, &
                                           'dates must be strictly ascending', &
                                           'soil.tillage.events.date')
                     end if
                  end if
               end associate
            end do
         end if

         if (allocated(self%tillage%types)) then
            do i = 1, size(self%tillage%types)
               associate (ty => self%tillage%types(i))
                  call check_real_range(ty%rho_cons,    100.0_real64, 3000.0_real64, &
                                        'soil.tillage.types.rho_cons',    errors)
                  call check_real_range(ty%rho_tillage, 100.0_real64, 3000.0_real64, &
                                        'soil.tillage.types.rho_tillage', errors)
                  call check_real_range(ty%k_R,         1.0e-4_real64, 10.0_real64, &
                                        'soil.tillage.types.k_R',         errors)
                  if (self%tillage%i_n_model == 3) then
                     call check_real_range(ty%rho_match, 100.0_real64, 3000.0_real64, &
                                           'soil.tillage.types.rho_match', errors)
                     call check_real_range(ty%N_match,   1.001_real64, 10.0_real64, &
                                           'soil.tillage.types.N_match',   errors)
                  end if
               end associate
            end do
         end if
      end if
```

Note: the `i` loop variable already exists in `soil_config_validate` (used by the cofani check) — reuse it.

If `ERR_VALIDATION_EMPTY` / `ERR_VALIDATION_RANGE` aren't already exported by `error_mod`, use the existing parse-error or generic-validation code that `check_real_range` itself uses; cross-reference `src/error/error.f90` for the canonical names. Validator helper signatures (`check_int_range`, `check_real_range`, `errors%append`) are already used elsewhere in this file — copy verbatim.

- [ ] **Step 4: Build and run tests to verify they pass**

```
pixi run -e test build-linux
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:"
```
Expected: build clean; `Ok: 1, Fail: 0`.

Also confirm the date-window check (events within `simulation.start_date`/`end_date`) is **not** in this task — it requires cross-section access to `swap_config_t` and lives in the adapter or a higher-level cross-validator. Track as a follow-up if needed; the per-block validator only enforces local rules.

- [ ] **Step 5: Commit**

```bash
git add src/config/soil_config.f90 tests/unit/io/toml/test_soil_tillage_validate.pf tests/unit/meson.build
git commit -m "$(cat <<'EOF'
feat(config): validate [soil.tillage] when swtill=1

Extends soil_config_validate with rules for the new tillage block:
non-empty events/types, type_id set-membership, ascending dates,
i_n_model/iRedist enums, intensity 0..1, rho_cons/rho_tillage 100..3000,
k_R 1e-4..10, optional rho_match/N_match required and ranged when
i_n_model=3. Depth-against-grid check deferred to apply_soil_tillage
(NumNod isn't available here). Cross-section date-window check
deferred to a higher-level validator.

Validator stays silent when swtill=0.

pFUnit suite test_soil_tillage_validate exercises every rule.

Part of tillage TOML port (spec
docs/superpowers/specs/2026-05-06-tillage-toml-port-design.md).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: TOML reader

**Files:**
- Create: `src/io/toml/read_soil_tillage_toml.f90`
- Modify: `src/io/toml/read_soil_toml.f90` (call the new sibling reader after parsing `[soil]`)
- Modify: `meson.build` (add the new file to production sources)
- Modify: `tests/unit/meson.build` (add to `pfunit_extra_sources`)
- Test: `tests/unit/io/toml/test_read_soil_tillage_toml.pf` (new)
- Modify: `tests/unit/meson.build` (register the new `.pf` file)

- [ ] **Step 1: Write the failing test**

Create `tests/unit/io/toml/test_read_soil_tillage_toml.pf`:

```fortran
@test
subroutine test_read_full_tillage_block()
   use soil_config_mod, only: soil_config_t
   use load_swap_config_mod, only: load_swap_config
   use error_mod, only: error_collection_t
   use funit
   implicit none
   type(soil_config_t)      :: soil
   type(error_collection_t) :: errs
   character(len=*), parameter :: toml_text = &
'[general]'                                  // new_line('a') // &
'project = "tillage_test"'                   // new_line('a') // &
'[soil]'                                     // new_line('a') // &
'swtill = 1'                                 // new_line('a') // &
'[soil.tillage]'                             // new_line('a') // &
'i_n_model = 2'                              // new_line('a') // &
'iRedist   = 2'                              // new_line('a') // &
'[[soil.tillage.events]]'                    // new_line('a') // &
'date      = 2003-04-15'                     // new_line('a') // &
'z         = 30.0'                           // new_line('a') // &
'intensity = 1.0'                            // new_line('a') // &
'type_id   = 1'                              // new_line('a') // &
'[[soil.tillage.events]]'                    // new_line('a') // &
'date      = 2004-04-20'                     // new_line('a') // &
'z         = 30.0'                           // new_line('a') // &
'intensity = 1.0'                            // new_line('a') // &
'type_id   = 1'                              // new_line('a') // &
'[[soil.tillage.types]]'                     // new_line('a') // &
'id          = 1'                            // new_line('a') // &
'rho_cons    = 1500.0'                       // new_line('a') // &
'rho_tillage = 1100.0'                       // new_line('a') // &
'k_R         = 0.05'

   ! Use whatever in-memory loader pattern the existing parity tests use;
   ! e.g. write to a temp file, call load_swap_config, fetch config%soil.
   ! Reference: tests/unit/io/toml/test_read_irrigation_toml.pf for a
   ! comparable pattern.
   !
   ! ... load via temp-file or an in-memory hook ...
   ! soil = config%soil

   @assertEqual(2, soil%tillage%i_n_model)
   @assertEqual(2, soil%tillage%iRedist)
   @assertTrue(allocated(soil%tillage%events))
   @assertEqual(2, size(soil%tillage%events))
   @assertEqual('2003-04-15', soil%tillage%events(1)%date)
   @assertEqual(30.0d0,       soil%tillage%events(1)%z, tolerance=1.0d-9)
   @assertEqual(1,            soil%tillage%events(1)%type_id)
   @assertTrue(allocated(soil%tillage%types))
   @assertEqual(1, size(soil%tillage%types))
   @assertEqual(1500.0d0, soil%tillage%types(1)%rho_cons,    tolerance=1.0d-9)
   @assertEqual(1100.0d0, soil%tillage%types(1)%rho_tillage, tolerance=1.0d-9)
   @assertEqual(0.05d0,   soil%tillage%types(1)%k_R,         tolerance=1.0d-12)
end subroutine

@test
subroutine test_read_missing_block_keeps_defaults()
   ! Same shape but no [soil.tillage] block; assert defaults intact:
   ! soil%tillage%i_n_model == 2 (default)
   ! .not. allocated(soil%tillage%events)
   ! .not. allocated(soil%tillage%types)
end subroutine

@test
subroutine test_read_optional_rho_match_N_match_when_i_n_model_3()
   ! TOML with i_n_model=3 and rho_match/N_match present in types(1);
   ! assert they read into the typed config.
end subroutine
```

Read the existing parity-test loader pattern in `tests/unit/io/toml/test_hupselbrook_loads.pf` (or a sibling) for the exact temp-file vs in-memory loading idiom this codebase uses.

- [ ] **Step 2: Wire the test file and run to verify it fails**

Add to `tests/unit/meson.build`'s `pf_files` list:
```meson
        'io/toml/test_read_soil_tillage_toml.pf',
```

Run:
```
pixi run -e test build-linux
```
Expected: build fails — `read_soil_tillage_toml` module / subroutine doesn't exist yet, or the test compiles but the new fields are never populated and the assertions fail.

- [ ] **Step 3: Create the reader module**

Create `src/io/toml/read_soil_tillage_toml.f90`. Pattern to copy: `src/io/toml/read_drainage_toml.f90`'s `read_levels` (array-of-tables handling) plus `src/io/toml/read_irrigation_toml.f90`'s date parsing (`parse_date_to_days1900` / `toml_datetime`).

```fortran
!> @file read_soil_tillage_toml.f90
!! Parses the optional [soil.tillage] block into soil_config_t%tillage.
!! Block is omitted when soil.swtill = 0; reader leaves defaults intact
!! in that case. Sibling to read_soil_toml.f90 to keep that file focused.
module read_soil_tillage_toml_mod

   use, intrinsic :: iso_fortran_env, only: real64
   use tomlf, only: toml_table, toml_array, toml_datetime, get_value, len
   use toml_field_helpers_mod, only: get_table, get_array_of_tables, &
                                     get_optional_int_with_default,  &
                                     get_optional_real_with_default
   use soil_config_mod, only: soil_tillage_t, soil_tillage_event_t, &
                              soil_tillage_type_t
   use error_mod, only: error_collection_t, ERR_PARSE_TYPE_MISMATCH
   implicit none
   private

   public :: read_soil_tillage_toml

contains

   !> Populate `tillage` from the optional [soil.tillage] table on `soil_sec`.
   !! Missing block is benign — leaves defaults.
   subroutine read_soil_tillage_toml(soil_sec, tillage, errors)
      type(toml_table), pointer, intent(in)    :: soil_sec
      type(soil_tillage_t),      intent(inout) :: tillage
      type(error_collection_t),  intent(inout) :: errors

      type(toml_table), pointer :: till_tbl, item
      type(toml_array), pointer :: events_arr, types_arr
      integer :: i, n, stat
      character(len=:), allocatable :: date_str

      if (.not. associated(soil_sec)) return

      ! [soil.tillage] is optional.
      till_tbl => null()
      call get_value(soil_sec, 'tillage', till_tbl, requested=.false., stat=stat)
      if (stat /= 0 .or. .not. associated(till_tbl)) return

      call get_optional_int_with_default(till_tbl, 'i_n_model', tillage%i_n_model, 2, &
                                         'soil.tillage.i_n_model', errors)
      call get_optional_int_with_default(till_tbl, 'iRedist',   tillage%iRedist,   2, &
                                         'soil.tillage.iRedist',   errors)

      ! [[soil.tillage.events]]
      call get_array_of_tables(till_tbl, 'events', events_arr, &
                               'soil.tillage.events', errors)
      if (associated(events_arr)) then
         n = len(events_arr)
         if (n > 0) then
            allocate(tillage%events(n))
            do i = 1, n
               item => null()
               call get_value(events_arr, i, item, stat=stat)
               if (stat /= 0 .or. .not. associated(item)) cycle
               call read_event_row(item, tillage%events(i), errors)
            end do
         else
            allocate(tillage%events(0))
         end if
      end if

      ! [[soil.tillage.types]]
      call get_array_of_tables(till_tbl, 'types', types_arr, &
                               'soil.tillage.types', errors)
      if (associated(types_arr)) then
         n = len(types_arr)
         if (n > 0) then
            allocate(tillage%types(n))
            do i = 1, n
               item => null()
               call get_value(types_arr, i, item, stat=stat)
               if (stat /= 0 .or. .not. associated(item)) cycle
               call read_type_row(item, tillage%types(i), errors)
            end do
         else
            allocate(tillage%types(0))
         end if
      end if
   end subroutine read_soil_tillage_toml


   subroutine read_event_row(row, ev, errors)
      type(toml_table), pointer, intent(in)    :: row
      type(soil_tillage_event_t), intent(out)  :: ev
      type(error_collection_t),  intent(inout) :: errors
      type(toml_datetime) :: dtv
      integer             :: stat

      ! Date — TOML local-date literal converted to days-since-1900 by the
      ! adapter; here we keep the ISO string for round-trip fidelity.
      call get_value(row, 'date', dtv, stat=stat)
      if (stat /= 0) then
         call errors%append(ERR_PARSE_TYPE_MISMATCH, &
                            "missing or non-date 'date' key", &
                            'soil.tillage.events')
      else
         write(ev%date, '(i4.4,"-",i2.2,"-",i2.2)') &
            dtv%date%year, dtv%date%month, dtv%date%day
      end if

      call get_optional_real_with_default(row, 'z',         ev%z,         0.0_real64, &
                                          'soil.tillage.events.z',         errors)
      call get_optional_real_with_default(row, 'intensity', ev%intensity, 0.0_real64, &
                                          'soil.tillage.events.intensity', errors)
      call get_optional_int_with_default(row,  'type_id',   ev%type_id,   0, &
                                         'soil.tillage.events.type_id',   errors)
   end subroutine read_event_row


   subroutine read_type_row(row, ty, errors)
      type(toml_table), pointer, intent(in)    :: row
      type(soil_tillage_type_t), intent(out)   :: ty
      type(error_collection_t),  intent(inout) :: errors
      call get_optional_int_with_default(row,  'id',          ty%id,          0, &
                                         'soil.tillage.types.id',          errors)
      call get_optional_real_with_default(row, 'rho_cons',    ty%rho_cons,    0.0_real64, &
                                          'soil.tillage.types.rho_cons',    errors)
      call get_optional_real_with_default(row, 'rho_tillage', ty%rho_tillage, 0.0_real64, &
                                          'soil.tillage.types.rho_tillage', errors)
      call get_optional_real_with_default(row, 'k_R',         ty%k_R,         0.0_real64, &
                                          'soil.tillage.types.k_R',         errors)
      call get_optional_real_with_default(row, 'rho_match',   ty%rho_match,   -99.0_real64, &
                                          'soil.tillage.types.rho_match',   errors)
      call get_optional_real_with_default(row, 'N_match',     ty%N_match,     -99.0_real64, &
                                          'soil.tillage.types.N_match',     errors)
   end subroutine read_type_row

end module read_soil_tillage_toml_mod
```

Note: this stores the date as an ISO string in `event_t%date`. The adapter (Task 4) parses that back into days-since-1900. This is a deliberate split — keeps the reader concerned with on-disk shape, the adapter concerned with global-state semantics. If you'd rather store a `real(real64)` directly in the event type, change `event_t%date` to `real(real64) :: t1900` and call `parse_date_to_days1900(dtv)` here. Document the choice in the module header.

- [ ] **Step 4: Wire the new module into `read_soil_toml.f90`**

In `src/io/toml/read_soil_toml.f90`, add a `use` line:
```fortran
   use read_soil_tillage_toml_mod, only: read_soil_tillage_toml
```

And after the existing `[soil]` parsing, add (after the `initial`/`hydraulics`/`frost`/`discretization` sub-table calls):
```fortran
   call read_soil_tillage_toml(sec, config%tillage, errors)
```

- [ ] **Step 5: Add the file to meson sources**

In `meson.build`, in the production source list (where `read_soil_toml.f90` is registered), add:
```meson
    'src/io/toml/read_soil_tillage_toml.f90',
```

In `tests/unit/meson.build`, in `pfunit_extra_sources`, add:
```meson
        '../../src/io/toml/read_soil_tillage_toml.f90',
```

- [ ] **Step 6: Build and run tests to verify pass**

```
pixi run -e test build-linux
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:"
```
Expected: build clean; `Ok: 1, Fail: 0`.

- [ ] **Step 7: Commit**

```bash
git add src/io/toml/read_soil_tillage_toml.f90 src/io/toml/read_soil_toml.f90 meson.build tests/unit/meson.build tests/unit/io/toml/test_read_soil_tillage_toml.pf
git commit -m "$(cat <<'EOF'
feat(io): parse [soil.tillage] into soil_config_t%tillage

Adds read_soil_tillage_toml_mod, called from read_soil_toml after
the existing [soil] sub-table walks (initial / hydraulics / frost /
discretization). Block is optional — when [soil.tillage] is absent,
reader leaves defaults intact. Date in `event.date` is kept as ISO
string here; adapter parses to days-since-1900.

Pattern: mirrors read_drainage_toml's array-of-tables walk plus
read_irrigation_toml's date handling.

pFUnit suite test_read_soil_tillage_toml covers full block, missing
block, and i_n_model=3 with optional rho_match/N_match.

Part of tillage TOML port (spec
docs/superpowers/specs/2026-05-06-tillage-toml-port-design.md).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: `apply_soil_tillage` adapter

**Files:**
- Modify: `src/io/toml/config_to_variables.f90` (add `apply_soil_tillage` helper; call from `config_to_variables`)
- Test: `tests/unit/io/toml/test_apply_soil_tillage.pf` (extend with adapter tests)

- [ ] **Step 1: Write the failing tests**

Append to the existing `tests/unit/io/toml/test_apply_soil_tillage.pf` (created in Task 1):

```fortran
@test
subroutine test_apply_populates_per_event_globals()
   use soil_config_mod, only: soil_config_t, soil_tillage_event_t, soil_tillage_type_t
   use config_to_variables_mod, only: apply_soil_tillage
   use variables, only: Date_tillage, Z_tillage, I_tillage, Type_tillage, &
                        Ntill, Max_Z_tillage, NumNod, zbotcp, tend
   use funit
   implicit none
   type(soil_config_t) :: cfg

   ! Minimal grid for the deferred z-range check.
   NumNod = 2
   if (allocated(zbotcp)) deallocate(zbotcp); allocate(zbotcp(2))
   zbotcp = [-50.0d0, -100.0d0]
   tend = 36898.0d0   ! 2004-12-31

   cfg%swtill = 1
   allocate(cfg%tillage%events(2))
   cfg%tillage%events(1) = soil_tillage_event_t(date='2003-04-15', z=30.0d0, &
                                                intensity=1.0d0, type_id=1)
   cfg%tillage%events(2) = soil_tillage_event_t(date='2004-04-20', z=30.0d0, &
                                                intensity=0.8d0, type_id=1)
   allocate(cfg%tillage%types(1))
   cfg%tillage%types(1) = soil_tillage_type_t(id=1, rho_cons=1500.0d0, &
                                              rho_tillage=1100.0d0, k_R=0.05d0)

   call apply_soil_tillage(cfg%tillage)

   @assertEqual(2, Ntill)
   @assertEqual(30.0d0, Z_tillage(1), tolerance=1.0d-9)
   @assertEqual(1.0d0,  I_tillage(1), tolerance=1.0d-9)
   @assertEqual(1,      Type_tillage(1))
   ! Sentinel slot: Date_tillage(Ntill+1) = tend + 1
   @assertEqual(tend + 1.0d0, Date_tillage(Ntill+1), tolerance=1.0d-9)
   @assertEqual(30.0d0, Max_Z_tillage, tolerance=1.0d-9)
end subroutine

@test
subroutine test_apply_populates_per_type_globals()
   ! Same setup; assert iType_Tillage(1) == 1, TAB_Rho_cons(1) == 1500, etc.
end subroutine

@test
subroutine test_apply_iTT1_iTT2_indices_when_two_types()
   ! events of types 1 and 2; types(1)..(3) covering both ids;
   ! assert iTT1/iTT2 first/last positions per the legacy Read_Tillage logic.
end subroutine

@test
subroutine test_apply_z_outside_grid_raises_error()
   ! events(1).z = 200.0 with zbotcp(NumNod) = -100 → outside, expect error.
   ! This exercises the deferred z-range validation.
end subroutine

@test
subroutine test_apply_i_n_model_3_allocates_match_arrays()
   ! cfg%tillage%i_n_model = 3; types(1).rho_match = 1200, N_match = 2.0;
   ! assert TAB_Rho_match and TAB_N_match are allocated and populated.
end subroutine
```

Reference for the legacy semantics: read the deleted `subroutine Read_Tillage` body from git
(`git show HEAD~N:src/crop/tillage.f90`) — concretely, the loop that
fills `iTT1`/`iTT2`:
```fortran
do j = 1, Ntill
   do i = 1, Ntypes
      if (iTT1(j) == 0 .and. iType_Tillage(i) == j) iTT1(j) = i
      if (iTT1(j) >  0 .and. iType_Tillage(i) == j) iTT2(j) = i
   end do
end do
```

- [ ] **Step 2: Run to verify failures**

```
pixi run -e test build-linux
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:"
```
Expected: tests fail because `apply_soil_tillage` doesn't exist or globals are unallocated.

- [ ] **Step 3: Add `apply_soil_tillage` to `src/io/toml/config_to_variables.f90`**

Add a `use` for the new types:
```fortran
   use soil_config_mod, only: soil_config_t, soil_tillage_t
```

Add `apply_soil_tillage` to the module's `public` exports.

Insert this subroutine in the module's `contains` section (place near the end, after the other private adapters):

```fortran
   !> Apply [soil.tillage] config to legacy `variables` globals.
   !! Called from config_to_variables when flTillage is true. Allocates
   !! per-event and per-type arrays, populates them from the typed config,
   !! parses event dates to days-since-1900, sets the Ntill+1 sentinel,
   !! computes Max_Z_tillage and the iTT1/iTT2 first/last-position
   !! indices. Replaces the deleted Read_Tillage subroutine.
   subroutine apply_soil_tillage(tillage)
      use, intrinsic :: iso_fortran_env, only: real64
      use soil_config_mod, only: soil_tillage_t
      use variables, only: Date_tillage, Z_tillage, I_tillage, Type_tillage,        &
                           iType_Tillage, TAB_Rho_cons, TAB_Rho_tillage,            &
                           TAB_K_R_cons, TAB_Rho_match, TAB_N_match,                &
                           iTT1, iTT2, Ntill, Ntypes, Max_Z_tillage,                &
                           i_n_model, iRedist, NumNod, zbotcp, tend
      use error_mod, only: error_collection_t, ERR_VALIDATION_RANGE
      ! NOTE: the global error collector is the one already used by the
      ! adapter; reference how config_to_variables passes its errors arg
      ! around and follow that pattern.
      type(soil_tillage_t), intent(in) :: tillage

      integer :: i, j

      i_n_model = tillage%i_n_model
      iRedist   = tillage%iRedist

      Ntill  = size(tillage%events)
      Ntypes = size(tillage%types)

      ! Per-event arrays (sentinel: Date_tillage(Ntill+1) = tend + 1)
      if (allocated(Date_tillage)) deallocate(Date_tillage); allocate(Date_tillage(Ntill+1))
      if (allocated(Z_tillage))    deallocate(Z_tillage);    allocate(Z_tillage(Ntill))
      if (allocated(I_tillage))    deallocate(I_tillage);    allocate(I_tillage(Ntill))
      if (allocated(Type_tillage)) deallocate(Type_tillage); allocate(Type_tillage(Ntill))

      do i = 1, Ntill
         Z_tillage(i)    = tillage%events(i)%z
         I_tillage(i)    = tillage%events(i)%intensity
         Type_tillage(i) = tillage%events(i)%type_id
         Date_tillage(i) = parse_iso_date_to_days1900(tillage%events(i)%date)
      end do
      Date_tillage(Ntill + 1) = tend + 1.0_real64

      ! Deferred z-range check (validator can't see NumNod / zbotcp).
      do i = 1, Ntill
         if (Z_tillage(i) < 0.0_real64 .or. &
             Z_tillage(i) > abs(zbotcp(NumNod))) then
            ! Append to the in-flight error collection. Use whatever the
            ! enclosing config_to_variables call passes; a singleton or
            ! pass-through, depending on how other adapter helpers do it.
            ! Match the pattern of, e.g., apply_soil_initial (search the
            ! file for an analogous deferred check).
         end if
      end do

      ! Per-type arrays
      if (allocated(iType_Tillage))   deallocate(iType_Tillage);   allocate(iType_Tillage(Ntypes))
      if (allocated(TAB_Rho_cons))    deallocate(TAB_Rho_cons);    allocate(TAB_Rho_cons(Ntypes))
      if (allocated(TAB_Rho_tillage)) deallocate(TAB_Rho_tillage); allocate(TAB_Rho_tillage(Ntypes))
      if (allocated(TAB_K_R_cons))    deallocate(TAB_K_R_cons);    allocate(TAB_K_R_cons(Ntypes))

      do i = 1, Ntypes
         iType_Tillage(i)   = tillage%types(i)%id
         TAB_Rho_cons(i)    = tillage%types(i)%rho_cons
         TAB_Rho_tillage(i) = tillage%types(i)%rho_tillage
         TAB_K_R_cons(i)    = tillage%types(i)%k_R
      end do

      if (i_n_model == 3) then
         if (allocated(TAB_Rho_match)) deallocate(TAB_Rho_match); allocate(TAB_Rho_match(Ntypes))
         if (allocated(TAB_N_match))   deallocate(TAB_N_match);   allocate(TAB_N_match(Ntypes))
         do i = 1, Ntypes
            TAB_Rho_match(i) = tillage%types(i)%rho_match
            TAB_N_match(i)   = tillage%types(i)%N_match
         end do
      end if

      Max_Z_tillage = maxval(Z_tillage(1:Ntill))

      ! iTT1 / iTT2 first/last-position-per-type indices, mirroring the
      ! deleted Read_Tillage:
      if (allocated(iTT1)) deallocate(iTT1); allocate(iTT1(Ntill)); iTT1 = 0
      if (allocated(iTT2)) deallocate(iTT2); allocate(iTT2(Ntill)); iTT2 = 0
      do j = 1, Ntill
         do i = 1, Ntypes
            if (iTT1(j) == 0 .and. iType_Tillage(i) == j) iTT1(j) = i
            if (iTT1(j) >  0 .and. iType_Tillage(i) == j) iTT2(j) = i
         end do
      end do
   end subroutine apply_soil_tillage


   !> ISO 'YYYY-MM-DD' -> days since 1900 (real(real64)).
   !! Mirrors the legacy SWAP date convention used elsewhere; reuse the
   !! helper from toml_field_helpers if a `parse_date_string_to_days1900`
   !! already exists; otherwise this implements it via toml_datetime.
   function parse_iso_date_to_days1900(s) result(t)
      use tomlf, only: toml_datetime
      use toml_field_helpers_mod, only: parse_date_to_days1900
      character(len=*), intent(in) :: s
      real(real64) :: t
      type(toml_datetime) :: dtv
      integer :: y, m, d
      read(s, '(i4,1x,i2,1x,i2)') y, m, d
      dtv%date%year  = y
      dtv%date%month = m
      dtv%date%day   = d
      t = parse_date_to_days1900(dtv)
   end function parse_iso_date_to_days1900
```

If the existing adapter's deferred-check pattern uses a module-scope or pass-through `error_collection_t`, copy that exact pattern — search `config_to_variables.f90` for an analogous "I check something that needs the grid" comment. The validator's per-block rules already caught everything that doesn't depend on the grid, so the only deferred check here is z-range.

In `config_to_variables` (the existing top-level subroutine in the same file), wire the call. Right after the existing `flTillage = (config%soil%swtill == 1)` line at lines 482-485, add:
```fortran
      if (flTillage) call apply_soil_tillage(config%soil%tillage)
```

- [ ] **Step 4: Build and run tests to verify pass**

```
pixi run -e test build-linux
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:"
pixi run -e test check-full 2>&1 | grep -E "Results:"
```
Expected: build clean; pFUnit `Ok: 1, Fail: 0`; check-full `5 passed, 0 failed` (all five regression cases still have `swtill=0`, so the new adapter call short-circuits).

- [ ] **Step 5: Commit**

```bash
git add src/io/toml/config_to_variables.f90 tests/unit/io/toml/test_apply_soil_tillage.pf
git commit -m "$(cat <<'EOF'
feat(io): apply_soil_tillage populates legacy globals from typed config

Adds the apply_soil_tillage helper to config_to_variables_mod, called
when flTillage = .true.. Allocates per-event (Date_tillage(Ntill+1),
Z_tillage, I_tillage, Type_tillage) and per-type (iType_Tillage,
TAB_Rho_cons, TAB_Rho_tillage, TAB_K_R_cons, optional TAB_Rho_match
and TAB_N_match) arrays, copies values from soil_config_t%tillage,
parses events.date (ISO YYYY-MM-DD) to days-since-1900 via the existing
parse_date_to_days1900 helper, sets the Ntill+1 sentinel
(tend + 1), computes Max_Z_tillage, and populates iTT1/iTT2.

Replaces the still-extant Read_Tillage swpfile-based reader (deletion
in next commit). check-full regression suite unchanged (all five cases
have swtill=0, adapter short-circuits).

Part of tillage TOML port (spec
docs/superpowers/specs/2026-05-06-tillage-toml-port-design.md).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: Retire `Read_Tillage`

**Files:**
- Modify: `src/crop/tillage.f90` (delete `subroutine Read_Tillage` at lines 423-508; remove the `call Read_Tillage` line in `DoTillage(1)` near line 64; clean up `swpfile`/`logf` imports)

- [ ] **Step 1: Verify the adapter is doing the work**

Confirm Task 4 already runs `apply_soil_tillage` before any `DoTillage(1)` call site. Trace:

```
src/io/toml/config_to_variables.f90:flTillage = ... ; if (flTillage) call apply_soil_tillage(...)   # config-load
src/core/swap.f90:165:                  if (flTillage) call DoTillage(1)                            # iTask=1 init
src/crop/tillage.f90:DoTillage(1):      ... allocations ... call Read_Tillage ... stub-checks ...
```

So at `DoTillage(1)` entry, the per-event/per-type globals are already populated (by `apply_soil_tillage`). `Read_Tillage` becomes a no-op.

- [ ] **Step 2: Delete the `call Read_Tillage` line in `DoTillage(1)`**

In `src/crop/tillage.f90`, find the case(1) block (around line 56-90) and remove:
```fortran
      ! read input data
      call Read_Tillage
      
      ! no tillage required; leave DoTillage immediately
      if (swtill == 0) return
```

The second `if (swtill == 0) return` was a belt-and-suspenders against `Read_Tillage` having returned without doing anything (legacy disable path). It's now dead code: `DoTillage(1)` is reached only when `flTillage = (swtill == 1)`. Remove both lines.

The per-NumLay allocations (`Rho_tillage(NumLay)`, `Rho_cons(NumLay)`, etc.) and the stub-error checks below stay — they're independent of `Read_Tillage`.

- [ ] **Step 3: Delete `subroutine Read_Tillage` entirely**

In `src/crop/tillage.f90`, delete the entire subroutine body from `subroutine Read_Tillage` (line 423) through `end subroutine Read_Tillage` (line 508), inclusive. Also delete the `! ******** Read_Tillage ********` comment header above it.

- [ ] **Step 4: Drop `swpfile` and `logf` from imports if unused**

Check the top-of-file `use variables, only:` lines:
```bash
grep -n "swpfile\|logf" src/crop/tillage.f90
```
After deletion, both should be unreferenced. Edit the `use variables, only: ...` line to drop them.

If the `! Phase 4f-extend SS-B (ADR 0020)` comments at lines 53/72 referenced `swpfile`, remove or update them to point at ADR 0021.

- [ ] **Step 5: Build and run tests**

```
pixi run -e test build-linux
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:"
pixi run -e test check-full 2>&1 | grep -E "Results:"
```
Expected: build clean (no unresolved references to `Read_Tillage` anywhere); pFUnit `Ok: 1, Fail: 0`; check-full `5 passed, 0 failed`.

Final acceptance grep:
```
grep -n "subroutine Read_Tillage\b" src/
grep -n "swpfile" src/crop/tillage.f90
```
Both should return zero matches.

- [ ] **Step 6: Optionally drop the case-2/3 self-checks**

In `src/crop/tillage.f90:53` and `tillage.f90:72`, the comments per ADR 0020 commit message note "defense-in-depth; not load-bearing post-ADR 0020 but harmless". Now that `apply_soil_tillage` is the single source of truth and `flTillage` gates entry, these `if (swtill /= 1) return` self-checks can be removed for consistency. This is optional clean-up — leave for a follow-up if you want this commit smaller.

- [ ] **Step 7: Commit**

```bash
git add src/crop/tillage.f90
git commit -m "$(cat <<'EOF'
refactor(crop): retire Read_Tillage; tillage.f90 no longer reads swpfile

Removes the last TTutil-based read in tillage.f90:
- Deletes subroutine Read_Tillage (was lines 423-508).
- Removes the `call Read_Tillage` and the now-dead
  `if (swtill == 0) return` from DoTillage(1).
- Drops swpfile and logf from `use variables` since neither is
  referenced any more.

apply_soil_tillage (in config_to_variables.f90, previous commit) is
now the sole tillage-init path; it runs at config-load time, before
DoTillage(1) is called. The remaining DoTillage(1) body
(per-NumLay allocations + stub-error checks) is unchanged.

Acceptance:
- `grep -n "subroutine Read_Tillage\\b" src/` → no matches.
- `grep -n "swpfile" src/crop/tillage.f90` → no matches.
- check-full 5/5 (all five regression cases have swtill=0; the
  adapter call short-circuits).

Part of tillage TOML port (spec
docs/superpowers/specs/2026-05-06-tillage-toml-port-design.md).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: ADR 0021 + schema doc

**Files:**
- Create: `docs/adr/0021-tillage-toml-port.md`
- Modify: `docs/adr/index.md` (add the new ADR to the index, if there's an index)
- Modify: `docs/configuration-schema.md` (document the `[soil.tillage]` block)

- [ ] **Step 1: Write `docs/adr/0021-tillage-toml-port.md`**

```markdown
---
title: "ADR 0021 — Tillage parameters ported to [soil.tillage] TOML block"
date: 2026-05-06
status: accepted
---

# ADR 0021: Tillage parameters ported to `[soil.tillage]` TOML block

## Context

ADR 0019's "Update 2026-05-06" closed the umbrella physical deletion
of the legacy fixed-format reader, but two TTutil-based reads
survived as foot-noted exceptions: `Read_Tillage` and
`read_ssdi_input`. Both still opened the staged `swap.swp` via
`RDinit` to load parameters that had no TOML schema yet.

This ADR ports the tillage subsystem; SSDI follows in ADR 0022.

## Decision

Tillage parameters move to a new `[soil.tillage]` sub-table. The
`swtill` switch stays at `[soil].swtill` (least disruption to
existing TOML cases). Per-event and per-type tables are inline TOML
arrays-of-tables (`[[soil.tillage.events]]`, `[[soil.tillage.types]]`)
— small, tightly coupled, and benefit from same-file readability.

`Read_Tillage` is deleted. A new helper `apply_soil_tillage` in
`config_to_variables.f90` populates the legacy globals
(`Date_tillage`, `Z_tillage`, `Type_tillage`, `iType_Tillage`,
`TAB_Rho_cons`, `TAB_Rho_tillage`, `TAB_K_R_cons`, optional
`TAB_Rho_match` / `TAB_N_match`, plus `Ntill`, `Ntypes`,
`Max_Z_tillage`, `iTT1`, `iTT2`) at config-load time.

## Schema

```toml
[soil]
swtill = 1   # unchanged location

[soil.tillage]
i_n_model = 2   # 1..3
iRedist   = 2   # 0..2

[[soil.tillage.events]]
date      = 2003-04-15
z         = 30.0
intensity = 1.0
type_id   = 1

[[soil.tillage.types]]
id          = 1
rho_cons    = 1500.0
rho_tillage = 1100.0
k_R         = 0.05
# rho_match / N_match required only when i_n_model = 3
```

When `swtill = 0`, the `[soil.tillage]` block can be omitted.

## Consequences

- `swpfile` global pointer no longer read by `tillage.f90`.
- `Read_Tillage` deleted; `apply_soil_tillage` is the sole
  tillage-init path.
- Pattern reusable for SSDI (next ADR).
- `flCropNut` un-stub-erring + nutrient reactivation remains a
  separate arc.

## Tests

pFUnit suites under `tests/unit/io/toml/`:
- `test_read_soil_tillage_toml`
- `test_soil_tillage_validate`
- `test_apply_soil_tillage`

No `swtill = 1` regression case is added in this ADR; deferred
until a known-good legacy comparator is available.

## Acceptance

- `grep -n "subroutine Read_Tillage\\b" src/` → no matches.
- `grep -n "swpfile" src/crop/tillage.f90` → no matches.
- `pixi run -e test test-pfunit` → `Ok: 1, Fail: 0`.
- `pixi run -e test check-full` → `5 passed, 0 failed`.
```

- [ ] **Step 2: Add `[soil.tillage]` to `docs/configuration-schema.md`**

Find the existing `[soil]` section and append a `[soil.tillage]` sub-section. Document each field: name, type, range, default, and "required only when …" notes for `rho_match` / `N_match`. Use the same shape as the existing `[soil.discretization]` / `[soil.hydraulics]` sub-sections.

- [ ] **Step 3: Update `docs/adr/index.md` (if present)**

Append:
```markdown
| [0021](0021-tillage-toml-port.md) | Tillage parameters ported to `[soil.tillage]` TOML block | 2026-05-06 | accepted |
```

(Format: match the existing index columns.)

- [ ] **Step 4: Final verification**

```
pixi run -e test build-linux
pixi run -e test test-pfunit 2>&1 | grep -E "^Ok:|^Fail:"
pixi run -e test check-full 2>&1 | grep -E "Results:"
echo "=== acceptance greps ==="
grep -n "subroutine Read_Tillage\b" src/ || echo "OK: no Read_Tillage"
grep -n "swpfile" src/crop/tillage.f90 || echo "OK: no swpfile in tillage.f90"
```
Expected: all green; both greps return "OK: ...".

- [ ] **Step 5: Commit**

```bash
git add docs/adr/0021-tillage-toml-port.md docs/configuration-schema.md docs/adr/index.md
git commit -m "$(cat <<'EOF'
docs: ADR 0021 — tillage parameters ported to [soil.tillage]

Captures the decision to port the tillage subsystem from TTutil
swap.swp reads to a typed [soil.tillage] sub-table. Closes the
half of ADR 0019's "two stub-readers also survive" footnote that
covered tillage; ADR 0022 will close SSDI.

- docs/adr/0021-tillage-toml-port.md: new ADR with rationale,
  schema, consequences, and acceptance checks.
- docs/configuration-schema.md: documents [soil.tillage] fields,
  ranges, optionality.
- docs/adr/index.md: ADR 0021 row.

Closes the tillage TOML port spec
(docs/superpowers/specs/2026-05-06-tillage-toml-port-design.md).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Self-Review

**Spec coverage:**

| Spec section | Task |
|---|---|
| §Schema (`soil_tillage_event_t`, `soil_tillage_type_t`, `soil_tillage_t`) | Task 1 |
| §Validation rules (per-block) | Task 2 |
| §Reader (`read_soil_tillage_toml`) | Task 3 |
| §Adapter (`apply_soil_tillage`, deferred z-range check) | Task 4 |
| §Runtime change (delete `Read_Tillage`, drop `swpfile` from `tillage.f90`) | Task 5 |
| §Tests (3 pFUnit suites) | Tasks 1-4 (one suite added per task as the surface area lands) |
| §ADR 0021 | Task 6 |
| §`configuration-schema.md` update | Task 6 |

All spec sections covered.

**Acceptance criteria** (from spec):

- `pixi run -e test build-linux` clean — verified each task.
- `pixi run -e test test-pfunit` — verified each task.
- `pixi run -e test check-full` — verified Tasks 4, 5, 6.
- `grep -n "subroutine Read_Tillage" src/` — verified Task 5.
- `grep -n "swpfile" src/crop/tillage.f90` — verified Task 5.
- Hand-authored fixture exercises every validator rule — Task 2's `test_soil_tillage_validate.pf`.
- ADR 0021 committed — Task 6.
- `docs/configuration-schema.md` updated — Task 6.

**Type / signature consistency:**

- `soil_tillage_event_t` has fields `date` (string), `z`, `intensity`, `type_id` — used identically in Tasks 1, 2, 3, 4.
- `soil_tillage_type_t` fields `id`, `rho_cons`, `rho_tillage`, `k_R`, `rho_match`, `N_match` — consistent.
- `apply_soil_tillage(tillage)` signature in Task 4 matches the call site wired in Task 4.
- `parse_iso_date_to_days1900` introduced in Task 4 only; no earlier task references it.

**Open questions punted to implementation:**

- Cross-section validator for date-window (events within
  `simulation.start_date`/`end_date`). Per-block validator can't see
  this; spec mentions it but defers to a higher-level validator. If
  no higher-level validator exists, add a dedicated cross-validator
  pass; otherwise skip in Task 2 and address as a follow-up.
- Whether the deferred z-range check appends to the in-flight
  `error_collection_t` or aborts via `fatalerr_collected`. Match
  whatever existing adapter helpers in the same file do
  (Task 4 Step 3 notes this).
- Whether to drop the `if (swtill /= 1) return` self-checks at
  `tillage.f90:53,72` in Task 5 or defer. Optional in the plan.

These are pickable at execution time without re-planning.

---

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-05-06-tillage-toml-port.md`.
