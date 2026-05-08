---
title: "[nutrients] N2b — Timed amendments via CSV companion (ADR 0027 candidate)"
author: Mateusz Zawadzki
date: 2026-05-08
status: draft
---

# [nutrients] N2b — Timed soil management events via CSV companion

Third sub-arc of the `[nutrients]` umbrella. Adds a CSV
companion `events_file` to the `[nutrients]` block carrying
timed soil management events (fertilizer applications,
manure spreading, …). Extends `apply_nutrients` with CSV
staging + sort + group logic mirroring the deleted legacy
`SoilManagement(1)` block. Runtime path
(`SoilManagement(3)` calling `Wofost_SoilAmendents`) is
unchanged.

**Does NOT lift the runtime gate at `tillage.f90:73`.** N3's
job after both N2a + N2b land.

## Goal

Reach a state where:

1. A `[nutrients]` block can carry an optional
   `events_file = "amendments.csv"` pointing at a CSV
   companion with timed amendments (date, material,
   amount_kgha, volat_fraction).
2. `apply_nutrients` reads the CSV at config-load time,
   validates row contents, sorts by date, groups same-day
   dosages, populates the legacy globals (`MatNum(:)`,
   `Amend(:)` in kg/m², `VolaFrac(:)`, `TimeAmend(:)`,
   `NuAmend(:)`, `iAmend(:,:)`, `namend`, `isme`).
3. The runtime path (`SoilManagement(3)` calling
   `Wofost_SoilAmendents`) consumes the populated globals
   without code change. Once `flCropNut=true` (N3), this
   just works.
4. Absent `events_file` → `namend = 0`, `isme = 1`, no
   amendments. Natural-mineralization-only baseline (same
   as N2a alone).
5. `pixi run -e test test-pfunit` passes; new pFUnit tests
   cover parser, validator, adapter (incl. sort + group).
   `pixi run -e test check-full` is unchanged (no case
   enables flCropNut yet).
6. ADR 0027 captures the decision.

## Non-goals

- **Inline events array.** The user explicitly preferred CSV
  ("the toml would explode"). Multi-year operational fertilizer
  schedules can run to dozens of events; CSV scales without
  visual noise. Mirrors the SSDI events decision (ADR 0022).
- **Material-property overrides** (legacy `<project>.smm`).
  Out-of-scope across the whole `[nutrients]` umbrella; the 17
  hardcoded materials in `Wofost_SoilParameters` stay.
- **Lift `tillage.f90:73`.** N3's job. With N2a + N2b in place,
  the runtime is finally fully wired; N3 is the lift-and-test.
- **Add a regression case.** N3.
- **Tighten ranges on the 17 nutrient scalars** (still in
  cropwofost.nutrient — unranged per N1's decision).

## Schema

### Update to `src/config/nutrients_config.f90`

`nutrients_config_t` gains one new field:

```fortran
type :: nutrients_config_t
   logical :: present = .false.
   real(real64) :: sorp_coef = 0.0_real64
   character(len=:), allocatable :: events_file       ! NEW; relative to pathwork
   type(nutrients_initial_t) :: initial
contains
   procedure :: validate => nutrients_config_validate
   procedure :: finalize => nutrients_config_finalize
end type
```

No new typed-record for the events themselves — CSV row
staging happens at adapter time, not at config-load time.
Mirrors `apply_ssdi_mode0`'s pattern (commit `328db06`).

### TOML shape

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

`events_file` is **optional**. Absent or empty string → no
amendments at runtime; `namend = 0`. The adapter still runs
(unconditionally per N2a) but takes the no-events branch.

Resolution: `trim(pathwork) // trim(events_file)` — same
convention as `apply_ssdi_mode0` (mode-0 CSV staging) and the
existing `fixed_events_file` adapter for `swirfix=1` irrigation.

### CSV companion shape

```
date,material,amount_kgha,volat_fraction
2003-04-15,10,100.0,0.05
2003-06-20,1,25000.0,0.10
2003-06-20,3,500.0,0.02
```

Columns:
- `date` — ISO `YYYY-MM-DD`. `read_csv_table` (per ADR 0012)
  auto-converts to days-since-1900 real(real64) when
  column 1's header is `'date'`.
- `material` — integer in [1, 20]. Maps to `MatNum`; the
  17 hardcoded materials in `Wofost_SoilParameters` are at
  indices 1-17 (Cattle manure → Spruce needles); indices
  18-20 are unused slots (`maxmat=20`).
- `amount_kgha` — real in [0, 500000] kg/ha (legacy
  `rdfdor` range preserved).
- `volat_fraction` — real in [0, 1] (tightened from legacy's
  typo'd `[0, 500000]` — physically a fraction).

Same-day rows are valid; the adapter groups them.

Up to `maxamn = 1000` rows (legacy parameter cap in
`wofost_soil_declarations.f90`).

## Validation

In `nutrients_config_validate`, no new rules added in this
sub-arc beyond what N2a already enforces. The CSV row
validation lives in the adapter — `events_file` is just a
file path at config-validate time.

**Adapter-time validation (deferred from validator):**

| Rule | Where it fires |
|---|---|
| `events_file` resolves to an existing file | `apply_nutrients`, after `read_csv_table` (which logs a clear "cannot open CSV" error) |
| Row count `<= maxamn` (1000) | `apply_nutrients` |
| `material(i)` in [1, 20] | `apply_nutrients`, per-row |
| `amount_kgha(i)` in [0, 500000] | `apply_nutrients`, per-row |
| `volat_fraction(i)` in [0, 1] | `apply_nutrients`, per-row |
| Optional cross-section: dates within `[tstart, tend]` | not enforced in N2b — Wofost_SoilAmendents simply won't fire for out-of-window dates; future enhancement |

All adapter-side errors route through `fatalerr_collected`
(the established pattern from `apply_ssdi_mode0`).

## Adapter — extend `apply_nutrients`

`src/io/toml/config_to_variables.f90`'s existing
`apply_nutrients(cfg)` (added in N2a, commit `5fb0e0b`) gets
a new private helper `apply_nutrients_events(cfg)` and a
call to it. The shape:

```fortran
subroutine apply_nutrients(cfg)
   ...
   ! N2a logic: SorpCoef + initial pools (unchanged)
   SorpCoef = cfg%sorp_coef
   ! ... pool population ...

   ! N2b: timed amendments
   call apply_nutrients_events(cfg)
end subroutine

subroutine apply_nutrients_events(cfg)
   use, intrinsic :: iso_fortran_env, only: real64
   use csv_reader_mod, only: read_csv_table
   use error_mod, only: error_collection_t, fatalerr_collected
   use nutrients_config_mod, only: nutrients_config_t
   use variables, only: pathwork
   use wofost_soil_declarations, only: MatNum, Amend, VolaFrac,           &
                                         TimeAmend, NuAmend, iAmend,        &
                                         namend, isme, maxamn
   type(nutrients_config_t), intent(in) :: cfg

   real(real64), allocatable :: tbl(:,:)
   type(error_collection_t)  :: errs
   character(len=300) :: csvpath
   character(len=14)  :: hdr(4)
   integer :: i, j, n
   real(real64) :: tmp_date, tmp_amount, tmp_volat
   integer :: tmp_mat

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
            tmp_date    = tbl(i, 1); tbl(i, 1) = tbl(j, 1); tbl(j, 1) = tmp_date
            tmp_mat     = nint(tbl(i, 2))
            tbl(i, 2)   = tbl(j, 2); tbl(j, 2) = real(tmp_mat, real64)
            tmp_amount  = tbl(i, 3); tbl(i, 3) = tbl(j, 3); tbl(j, 3) = tmp_amount
            tmp_volat   = tbl(i, 4); tbl(i, 4) = tbl(j, 4); tbl(j, 4) = tmp_volat
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
   iAmend(1, 1) = 1
   do i = 2, n
      if (abs(tbl(i, 1) - tbl(i - 1, 1)) < 1.0e-3_real64) then
         NuAmend(j) = NuAmend(j) + 1
      else
         j = j + 1
         NuAmend(j)   = 1
         TimeAmend(j) = tbl(i, 1)
      end if
      iAmend(j, NuAmend(j)) = i
   end do

   namend = j
   isme   = 1
end subroutine apply_nutrients_events
```

The bubble-sort is acceptable for this use case (n ≤ 1000,
runs once at config-load). The grouping algorithm is a
verbatim port of the deleted pre-collapse `SoilManagement(1)`
logic — preserves bit-identical TimeAmend / NuAmend / iAmend
construction.

`maxamn` is exported by `wofost_soil_declarations` (parameter
= 1000). Confirm at implementation-time that
`wofost_soil_declarations` exports `maxamn` (the file declares
it as `Integer, Parameter :: maxamn = 1000` at line 11; verify
it's visible from the adapter via `use
wofost_soil_declarations, only: maxamn`).

## Runtime — no change

`SoilManagement(3)` already does:

```fortran
case (3)
   if (abs(TimeAmend(isme) + 1.0d0 - t1900) .lt. 1.d-3) then
      call Wofost_SoilAmendents
      isme = isme + 1
   endif
```

Once `flCropNut=true` (N3 lifts the gate), this fires once per
amendment-group date. `Wofost_SoilAmendents` reads `MatNum`,
`Amend`, `VolaFrac`, `iAmend`, `NuAmend`, `isme` (all populated
by the adapter) and applies the amendments to the soil pools.
Zero runtime code change in this sub-arc.

A subtle behaviour: the legacy check is `TimeAmend(isme) + 1
== t1900`, i.e., "today is the day **after** the amendment
date". This is intentional — the legacy applies amendments on
the day after the listed date (a one-day delay reflecting
typical agronomic timing). Preserve verbatim.

`isme` overflow protection: when `isme > namend`,
`TimeAmend(isme)` reads off the end. The legacy code didn't
guard against this; in practice `isme` only increments past
`namend` when no further amendments are scheduled, and the
condition `TimeAmend(isme) + 1 == t1900` is unlikely to fire
on a stale value. Leave as-is unless N3 testing surfaces a
concrete bug.

## Tests

Extend `tests/unit/io/toml/test_nutrients_config.pf`:

| Test | Asserts |
|---|---|
| `test_read_nutrients_with_events_file` | TOML with `events_file = "x.csv"`; `cfg%events_file` allocated to `"x.csv"`, `cfg%present == .true.` |
| `test_read_nutrients_no_events_file` | TOML with `[nutrients]` but no `events_file`; `cfg%events_file` unallocated |
| `test_apply_events_csv_populates_legacy_globals` | CSV fixture (3 rows, one same-day pair); call adapter; assert `namend == 2`, `NuAmend(1) == 2`, `NuAmend(2) == 1`, `Amend(:)` in kg/m² (1.0e-4 of CSV value), `MatNum(:)` correct, `TimeAmend(:)` sorted, `isme == 1` |
| `test_apply_events_no_file_zero_amendments` | `cfg%events_file` unallocated; `namend == 0`, `isme == 1` |
| `test_apply_events_rejects_material_out_of_range` | CSV with `material = 25` (column 2); expect adapter abort via `fatalerr_collected` |
| `test_apply_events_rejects_volat_fraction_above_1` | CSV with `volat_fraction = 1.5`; expect abort |

CSV fixture: `tests/unit/io/toml/fixtures/nutrients_events_small.csv`:
```
date,material,amount_kgha,volat_fraction
2003-06-20,1,25000.0,0.10
2003-04-15,10,100.0,0.05
2003-06-20,3,500.0,0.02
```

Note the first two rows are deliberately out-of-order so the
sort test is meaningful. After sort: 2003-04-15 (mat 10),
2003-06-20 (mat 1, then mat 3 — both grouped under
TimeAmend(2)).

`check-full`: 5/5 unchanged. No regression case sets
`flCropNut=true`.

## ADR 0027

`docs/adr/0027-nutrients-N2b-timed-amendments.md`:

- **Decision.** `events_file` (optional) on `[nutrients]`
  references a CSV companion with `date,material,amount_kgha,volat_fraction`
  columns. Adapter stages, validates per-row, sorts by date,
  groups same-day dosages, populates legacy globals.
- **Rationale.** CSV chosen over inline TOML because typical
  multi-year fertilizer schedules can grow to 50+ events;
  inline becomes unreadable at that scale. Same shape as SSDI
  events (ADR 0022).
- **Range tightening.** `volat_fraction` validator range
  changes from legacy's typo'd `[0, 500000]` to physically
  meaningful `[0, 1]`. Documented as a behavioural change.
- **Bubble-sort.** Acceptable at n ≤ 1000; runs once at
  config-load; mirrors the deleted legacy SoilManagement(1)
  algorithm.
- **Consequences.** With N2a + N2b in place, all soil-side
  nutrient inputs flow through TOML + CSV. N3 can lift the
  runtime gate.

`docs/adr/index.md` updated. `docs/configuration-schema.md`
gains an `events_file` row + CSV-format documentation under
`[nutrients]`.

## Acceptance criteria

- [ ] `pixi run -e test build-linux` clean.
- [ ] `pixi run -e test test-pfunit` exits 0; suite count
      grows by 6.
- [ ] `pixi run -e test check-full` → `5 passed, 0 failed`
      byte-identical CSV outputs.
- [ ] `grep -n "subroutine apply_nutrients_events" src/io/toml/config_to_variables.f90`
      → one match.
- [ ] CSV fixture committed at
      `tests/unit/io/toml/fixtures/nutrients_events_small.csv`.
- [ ] `tillage.f90:73` UNCHANGED — runtime stub-error stays.
- [ ] ADR 0027 committed.

## Commit cadence

~5 commits, each independently buildable + green:

1. `feat(config): add events_file field to nutrients_config_t`
2. `feat(io): parse events_file in read_nutrients_toml`
3. `feat(io): apply_nutrients_events — CSV staging + sort + group`
4. `test: pFUnit suite extensions for events`
5. `docs: ADR 0027 + index + schema doc`

(Commits 3 + 4 may merge at plan time. Plan-time decision.)

## Risk

- **CSV path resolution.** `pathwork` is set early in
  `config_to_variables` from `[general.paths].work` and is
  used by other CSV adapters (`apply_ssdi_mode0`, the irrigation
  fixed-events block). Reuse the same convention. Not a new
  failure mode.
- **`isme` overflow** when amendments run past simulation end
  (or start before). Legacy code didn't guard against it; the
  equality check `TimeAmend(isme) + 1 == t1900` is unlikely to
  match stale values. Documented as known-good for now; N3
  testing may surface a concrete bug requiring a guard.
- **`maxamn = 1000` cap.** A user with > 1000 events must
  either (a) split into multiple simulations or (b) we raise
  the cap (a one-line change in `wofost_soil_declarations`).
  Not gated on N2b.
- **`volat_fraction` range tightening** (legacy typo'd
  `[0, 500000]` → N2b's `[0, 1]`) is a behavioural change.
  Any input file in legacy use with `volat_fraction` outside
  `[0, 1]` would now fail validation. ADR 0027 documents the
  change; legacy values outside `[0, 1]` were physically
  meaningless anyway.
- **Adapter validation aborts via `fatalerr_collected`** rather
  than appending to `error_collection_t` for a full-pass
  validation. This matches the SSDI-mode-0 precedent.
  `error_collection_t` consolidation across adapters is a
  separate architectural arc (out of scope here).
