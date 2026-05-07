---
title: "SSDI TOML port (ADR 0022 candidate)"
author: Mateusz Zawadzki
date: 2026-05-07
status: draft
---

# SSDI TOML port — design spec

Sequel to `2026-05-06-tillage-toml-port-design.md`. ADR 0019's
"Update 2026-05-06" closed the umbrella legacy-reader physical
deletion, but two TTutil-based reads survived as foot-noted
exceptions: `Read_Tillage` and `read_ssdi_input`. ADR 0021 (already
landed) ported tillage. This spec ports SSDI — the second and final
half — and retires `src/io/checkdate.f90` along with it.

After this lands, `swpfile` has no readers anywhere in `src/`.

## Goal

Reach a state where:

1. SSDI parameters (`ssdi_schedule`, `ssdi_sched_type`, `ssdi_z`,
   `ssdi_amount`, `ssdi_appl_rate`, `ssdi_threshold`,
   `ssdi_threshold_z`, `sw_interval`, `days_interval`, plus the
   per-event `ssdi_date` / `ssdi_rate_f` / `ssdi_amount_f` table
   for fixed mode) are loaded from `[irrigation.ssdi]` in
   `swap.toml` via the typed pipeline.
2. Fixed-mode event tables (up to 366 rows per simulation) live in
   a CSV companion file, referenced from the TOML.
3. `read_ssdi_input` (in `src/crop/irrigation.f90`) is deleted; an
   `apply_irrigation_ssdi` helper in
   `src/io/toml/config_to_variables.f90` is the sole entry point
   for SSDI-init.
4. `irrigation.f90` no longer references `swpfile`, `logf`,
   `ssdi_file`, `RDinit`, `RDsinr`, `RDsdor`, `RDfdor`, `RDatim`,
   or `RDinar`. All TTutil dependencies in this file go away.
5. `src/io/checkdate.f90` is deleted (no remaining callers — the
   date-window check moves into the validator/adapter, expressed
   against `error_collection_t`).
6. Validators reject malformed `[irrigation.ssdi]` blocks (mode
   mismatch, out-of-range threshold per sched_type, missing
   `days_interval` when `sw_interval=1`, …). Errors collected via
   `error_collection_t`, not `STOP`.
7. pFUnit suites cover schema parsing, validation, and adapter
   behaviour. check-full stays 5/5 (no `swssdi=1` regression case
   added; deferred, same rationale as tillage).
8. ADR 0022 captures the decision; `configuration-schema.md`
   documents the new block + CSV format.

End-state: `grep -rn "swpfile" src/` returns no matches.

## Non-goals

- **Adding a `swssdi=1` regression case.** All five existing cases
  default to `swssdi=0` and are unaffected by this work. Authoring
  a sixth synthetic-SSDI case is a separate, follow-on task —
  preferably after we have a known-good `swssdi=1` path to compare
  the legacy binary against.
- **Removing `swpfile`/`logf` plumbing in
  `config_to_variables.f90`.** With no remaining readers, the
  pointers are dead weight. Cleanest as a follow-on commit; called
  out in ADR 0022's "next step" but out of scope for this spec.
- **Refactoring `SSDI_irrigation(2)`/`SSDI_irrigation(9)` runtime
  logic.** This spec is parameter-load-only. Case(2) (per-day event
  apply) and case(9) (end-of-event reset) stay unchanged.
- **Nutrient reactivation.** Separate umbrella; covered by its own
  spec arc.
- **Dual inline+CSV support for events.** Tables can hit 366 rows;
  inline becomes unreadable. CSV-only.

## Schema

### New types in `src/config/irrigation_config.f90`

```fortran
!> Fixed-mode SSDI configuration. Date table lives in a CSV file
!! pointed to by events_file (parsed by the adapter, not the TOML
!! reader, so the typed-config layer stays file-agnostic).
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

`irrigation_config_t` gains one new field:

```fortran
type(irrigation_ssdi_t) :: ssdi
```

`swssdi` stays at `[irrigation].swssdi` (its current location). The
`[irrigation.ssdi]` block is **optional**: when `swssdi = 0` the
block can be omitted entirely.

### TOML shape — fixed-date mode (schedule=0)

```toml
[irrigation]
swssdi = 1

[irrigation.ssdi]
schedule = 0
ssdi_z   = [-30.0, -50.0]   # depth interval in cm; for single-depth, both equal

[irrigation.ssdi.fixed]
events_file = "ssdi_events.csv"
```

CSV companion (`ssdi_events.csv`):
```
date,rate_f,amount_f
2003-04-15,5.0,10.0
2003-05-02,5.0,8.0
...
```

Columns:
- `date` — ISO `YYYY-MM-DD`, strictly ascending across rows.
- `rate_f` — mm/h, in [0.0, 100.0].
- `amount_f` — mm, in [0.0, 100.0].

Up to `mairg = 366` rows.

### TOML shape — scheduled-trigger mode (schedule=1)

```toml
[irrigation.ssdi]
schedule = 1
ssdi_z   = [-30.0, -50.0]

[irrigation.ssdi.scheduled]
sched_type      = 1            # 1=Tred, 2=presh, 3=watc
threshold       = 0.7          # semantics per sched_type
threshold_depth = -50.0        # cm; required when sched_type > 1
ssdi_amount     = 10.0         # mm
ssdi_appl_rate  = 5.0          # mm/h
sw_interval     = 0            # 0|1
# days_interval required only when sw_interval = 1
```

When `swssdi = 0`, the entire `[irrigation.ssdi]` block can be
omitted.

## Validation rules

In `irrigation_config_validate`, fired only when `self%swssdi == 1`:

| Rule | Error context |
|---|---|
| `schedule` in [0, 1] | `irrigation.ssdi.schedule` |
| `ssdi_z(i)` in [-100.0, 0.0] for `i = 1, 2` | `irrigation.ssdi.ssdi_z` |
| `ssdi_z(1) >= ssdi_z(2)` (top above or equal to bottom; cm increases downward as a magnitude) | `irrigation.ssdi.ssdi_z` |
| `schedule == 0` ⇒ `len_trim(fixed.events_file) > 0` | `irrigation.ssdi.fixed.events_file` |
| `schedule == 0` ⇒ scheduled-block fields all at defaults (hard error if user populated both — silently-ignored fields are a footgun) | `irrigation.ssdi.scheduled` |
| `schedule == 1` ⇒ `fixed.events_file` empty | `irrigation.ssdi.fixed.events_file` |
| `schedule == 1` ⇒ `scheduled.sched_type` in [1, 3] | `irrigation.ssdi.scheduled.sched_type` |
| `schedule == 1` ⇒ `scheduled.threshold` in mode-specific range: `sched_type=1` → [0, 1]; `=2` → [-1.0e7, 0]; `=3` → [0, 1] | `irrigation.ssdi.scheduled.threshold` |
| `schedule == 1 .and. sched_type > 1` ⇒ `scheduled.threshold_depth` in [-100, 0] | `irrigation.ssdi.scheduled.threshold_depth` |
| `schedule == 1` ⇒ `scheduled.ssdi_amount` in [0, 100] (mm) | `irrigation.ssdi.scheduled.ssdi_amount` |
| `schedule == 1` ⇒ `scheduled.ssdi_appl_rate` in [0, 100] (mm/h) | `irrigation.ssdi.scheduled.ssdi_appl_rate` |
| `schedule == 1` ⇒ `scheduled.sw_interval` in [0, 1] | `irrigation.ssdi.scheduled.sw_interval` |
| `schedule == 1 .and. sw_interval == 1` ⇒ `scheduled.days_interval` in [1, 366] | `irrigation.ssdi.scheduled.days_interval` |

**Deferred to adapter** (validator can't see the simulation window
or the staged CSV rows):

| Rule | Where it fires |
|---|---|
| CSV row count `>= 1` (mode 0) | `apply_irrigation_ssdi`, after CSV staging |
| CSV dates strictly ascending (mode 0) | `apply_irrigation_ssdi` |
| At least one CSV date in `[tstart, tend]`, OR `[tstart, tend] ⊆ [dates(1), dates(ifnd)]` (mode 0) | `apply_irrigation_ssdi`; this replaces the deleted `checkdate` call |
| `ssdi_z` resolved against `zbotcp(NumNod)` | `apply_irrigation_ssdi` |

All errors are appended to `error_collection_t`; the
finalize/abort checkpoint surfaces them at once (ADR 0008).

## CSV reader

The fixed-mode events CSV is parsed by the existing
`src/io/csv_reader.f90:read_csv_table` (per ADR 0012). Mirror the
adapter pattern from `src/drainage/surfacewater_init.f90` for
analogous date-table handling.

CSV reading happens at **adapter time** (parallel to how meteo CSV
is read), not at TOML-parse time — keeps the typed-config layer
file-agnostic.

The reader fills new staging fields on `irrigation_ssdi_fixed_t`:
either extend that type with `dates(:)`, `rate_f(:)`, `amount_f(:)`,
`ifnd` (preferred — keeps staging colocated with config), OR thread
locals through the adapter (alternative, but loses the
config-as-snapshot property).

## TOML reader

New module `src/io/toml/read_irrigation_ssdi_toml.f90`, sibling to
`read_irrigation_toml.f90`. Pattern mirrors `read_drainage_toml`'s
sub-table walk plus the existing soil-tillage reader (ADR 0021).
Called from `read_irrigation_toml` after the existing `[irrigation]`
parsing.

The block is optional: if absent, parser leaves
`config%irrigation%ssdi` at its default-initialized state (schedule=0,
files empty, scheduled struct at defaults).

## Adapter

`src/io/toml/config_to_variables.f90` gets a new helper:

```fortran
subroutine apply_irrigation_ssdi(ssdi)
   type(irrigation_ssdi_t), intent(in) :: ssdi
   ! ...
end subroutine
```

Called from `config_to_variables` when `flSSDI` is true (right
after the existing `flSSDI = (config%irrigation%swssdi == 1)`
line).

`apply_irrigation_ssdi` does:

1. Resolves `ssdi_z` array → `nod_ssdi(1:2)` layer indices via the
   `zbotcp` walk currently in `SSDI_irrigation(1)` lines 414-422.
2. Mode 0 (fixed):
   - Reads the CSV via `read_csv_table` (ADR 0012); stages
     `dates(:)`, `rate_f(:)`, `amount_f(:)`, `ifnd`.
   - Validates row count, ascending dates, date-window check
     (replaces the deleted `checkdate` call).
   - Allocates legacy globals `ssdi_date(mairg)`, `ssdi_rate_f(mairg)`,
     `ssdi_amount_f(mairg)`; copies values; sets `ifnd` global.
   - Computes initial entry-point `nirri` based on `t1900` (current
     case(1) lines 363-371).
   - Converts units to legacy: `rate_f` mm/h → cm/d
     (`* 0.1d0 * 24.0d0`); `amount_f` mm → cm (`* 0.1d0`).
   - Spreads `amount_f` over the `nod_ssdi(2) - nod_ssdi(1) + 1`
     compartments (current case(1) lines 425-426).
3. Mode 1 (scheduled):
   - Copies `sched_type`, `threshold`, `threshold_z`,
     `ssdi_amount`, `ssdi_appl_rate`, `sw_interval`,
     `days_interval` to legacy globals.
   - Sets `nirri = 1`, `dt_SSDI_event = 1.0d0`, `days_counter = 366`
     (current case(1) lines 372-380; preserves the `d8a88d6`
     regression-fix invariant).
   - Converts units: `ssdi_amount` mm → cm; `ssdi_appl_rate` mm/h
     → cm/d.
   - Spreads `ssdi_amount` over compartments
     (`amount = amount / dble(nod_ssdi(2) - nod_ssdi(1) + 1)`).
4. Both modes: `qssdi = 0.0d0`.

Names match the legacy globals' `ssdi_*_irr` variants where
relevant — see `src/crop/irrigation.f90:343-349` for the existing
mapping (`days_counter_irr`, `nirri_ssdi_irr`, `ssdi_date_irr`,
`ssdi_rate_f_irr`, `ssdi_amount_f_irr`).

## Runtime change in `src/crop/irrigation.f90`

- **Delete `subroutine read_ssdi_input`** (currently lines
  494-554). All TTutil reads + the `checkdate` call go with it.
- **Collapse `SSDI_irrigation(1)`** to a no-op (or delete the case
  body). Everything it used to do —`getun2`/`rdinit`/`rdsinr` of
  swap.swp, the `ssdi_file` open, the call to `read_ssdi_input`,
  the `ssdi_z` → `nod_ssdi` walk, the `nirri` initial entry
  computation, the unit conversions — is now done in
  `apply_irrigation_ssdi` at config-load time.
- The case-2 self-check at `irrigation.f90:634` (per ADR 0020
  commit's "harmless defense-in-depth" note) stays or gets removed
  — small follow-up either way; remove for consistency with the
  tillage retirement.
- The file's `use variables, only: ...` line drops `swpfile`,
  `logf`, `ssdi_file`, `swp` (and any other formerly-required
  TTutil-side names) if no longer referenced.

After this lands, `grep -n "swpfile\|rdinit\|RDinit" src/crop/irrigation.f90`
returns zero matches.

## Retire `src/io/checkdate.f90`

`checkdate` was extracted from the deleted `readswap.f90` in SS-C
step 5 specifically because `read_ssdi_input` called it. With
`read_ssdi_input` gone:

- The date-window check moves into `apply_irrigation_ssdi`'s
  deferred-validation tail. Same logic: at least one CSV date in
  `[tstart, tend]` OR `[tstart, tend] ⊆ [dates(1), dates(ifnd)]`.
  Now expressed against `error_collection_t` (`%append` with
  `ERR_VALIDATION_CROSS_FIELD`), not `fatalerr` + `STOP`.
- **Delete `src/io/checkdate.f90`.**
- Drop `'src/io/checkdate.f90'` from `meson.build` production
  sources.
- Drop `'../../src/io/checkdate.f90'` from
  `tests/unit/meson.build`'s `pfunit_extra_sources`.
- After deletion: `grep -rn "checkdate" src/` returns no matches.

## Tests

Three new pFUnit suites under `tests/unit/io/toml/`:

| File | Asserts |
|---|---|
| `test_read_irrigation_ssdi_toml.pf` | TOML parses correctly: schedule=0 with `events_file`, schedule=1 with `[scheduled]` block, missing block default, optional `days_interval` only meaningful when `sw_interval=1` |
| `test_irrigation_ssdi_validate.pf` | Validator catches every rule: schedule out of range, mode-mismatch (events_file populated when schedule=1, etc.), threshold out-of-range per sched_type, ssdi_z ordering, missing days_interval when sw_interval=1 |
| `test_apply_irrigation_ssdi.pf` | Adapter populates legacy globals: build minimal `swap_config_t` for both modes (mode 0 with a small CSV fixture), call `apply_irrigation_ssdi`, assert `ssdi_date`, `ssdi_rate_f`, `ssdi_amount_f`, `ifnd`, `nod_ssdi(1:2)`, `qssdi`, `dt_SSDI_event`, `nirri` match expectations; for mode 1 assert `sched_type`, `threshold`, etc. |

CSV fixture: `tests/unit/io/toml/fixtures/ssdi_events_small.csv`
(~3 rows).

`tests/unit/meson.build` lists the new `.pf` files +
`testSuites.inc` registers the suites. No changes to existing
parity suites.

`pixi run -e test check-full` is unchanged: 5/5 in ~42s; all
existing cases keep `swssdi=0` and the new adapter call
short-circuits.

## ADR

`docs/adr/0022-ssdi-toml-port.md` (new):

- **Decision.** SSDI parameters move to `[irrigation.ssdi]`
  sub-table with explicit `schedule = 0|1` discriminator and a
  CSV companion for fixed-mode events. Legacy `read_ssdi_input`
  retired. `checkdate.f90` retired (date-window check expressed
  against `error_collection_t`).
- **Rationale.** Closes the second half of ADR 0019's
  "two stub-readers also survive" footnote (tillage was the
  first; ADR 0021). SSDI was the last `swpfile` reader. CSV
  chosen for events because tables can hit 366 rows; explicit
  `schedule` discriminator chosen because the two modes have
  largely disjoint parameter sets.
- **Consequences.**
  - `swpfile` global pointer no longer read anywhere in `src/`.
    Cleanup of `swpfile`/`logf` plumbing in
    `config_to_variables.f90` is the next step.
  - `read_ssdi_input` deleted; `apply_irrigation_ssdi` is the
    sole SSDI-init path.
  - `checkdate.f90` deleted; `error_collection_t` is the sole
    validation channel.
  - With ADR 0021 + ADR 0022 landed, the only remaining TTutil
    runtime calls in `src/` are utility (`rdsets`/`rdfrom` in
    `swap_main.f90` + `rddtmp` in `swapoutput.f90`) — see ADR
    0019 §"Production runtime" for context.

`docs/configuration-schema.md` gains an `[irrigation.ssdi]`
section documenting fields, ranges, the events-CSV format, and
both mode shapes.

## Acceptance criteria

- [ ] `pixi run -e test build-linux` clean.
- [ ] `pixi run -e test test-pfunit` exits 0; suite count grows
      by 3 new suites; existing 555 tests unchanged.
- [ ] `pixi run -e test check-full` exits 0 with `5 passed,
      0 failed`.
- [ ] `grep -n "subroutine read_ssdi_input" src/` returns zero
      matches.
- [ ] `grep -rn "swpfile" src/` returns zero matches.
- [ ] `grep -rn "checkdate" src/` returns zero matches.
- [ ] `src/io/checkdate.f90` does not exist.
- [ ] A hand-authored fixture exercises every validator rule.
- [ ] ADR 0022 committed.
- [ ] `docs/configuration-schema.md` updated.

## Commit cadence

~6-7 commits, each independently buildable + green:

1. `feat(config): add irrigation_ssdi_t types + validator`
2. `feat(io): parse [irrigation.ssdi] in TOML reader`
3. `feat(io): apply_irrigation_ssdi adapter helper (mode-0 CSV staging + mode-1 copy)`
4. `refactor(crop): retire read_ssdi_input; collapse SSDI_irrigation(1)`
5. `refactor(io): retire src/io/checkdate.f90 — no remaining callers`
6. `test: pFUnit suites for irrigation.ssdi schema + adapter`
7. `docs: ADR 0022 + configuration-schema.md update`

The plan (writing-plans) will refine ordering — commit (5)
depends on (4) being merged first; (6) lands tests alongside
each implementation step rather than as a single tail commit.

## Risk

- **CSV fixture path resolution.** Adapter must resolve
  `events_file` relative to the TOML directory, not the binary's
  cwd. Same pattern as meteo / drainage CSV companions; reuse
  `path_helpers_mod`.
- **Date-window check semantics.** Legacy `checkdate` accepted
  EITHER "at least one date in window" OR "window ⊆ date range".
  The validator must preserve both branches (see source at
  `src/io/checkdate.f90:21-31` before deletion).
- **Mode-discriminator validation false negatives.** A user could
  populate `[irrigation.ssdi.scheduled]` while setting
  `schedule = 0`. Validator should at minimum warn (or error)
  rather than silently ignore the populated fields. Decided as
  hard error in §Validation rules above.
- **`ssdi_z` ordering convention.** Legacy stores depths as
  negative cm (downward-negative). `ssdi_z(1) >= ssdi_z(2)` means
  top is closer to the surface (less negative). Document this in
  the configuration-schema.md TOML example.
