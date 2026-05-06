---
title: "Tillage TOML port (ADR 0021 candidate)"
author: Mateusz Zawadzki
date: 2026-05-06
status: draft
---

# Tillage TOML port — design spec

Sequel to `2026-05-05-legacy-readers-physical-deletion-design.md`.
That umbrella deleted the legacy fixed-format reader from `src/`, but
left two small TTutil-based reads alive: `Read_Tillage` (this spec)
and `read_ssdi_input` (next spec). Both still open the staged
`swap.swp` via `RDinit` to load parameters that have no TOML schema
yet.

This spec ports the tillage subsystem to a typed `[soil.tillage]`
TOML block, retires `Read_Tillage`, and removes `tillage.f90`'s
dependence on `swpfile`. SSDI follows in the next spec; both close
ADR 0019's "two stub-readers also survive" footnote.

## Goal

Reach a state where:

1. Tillage parameters (`i_n_model`, `iRedist`, `Date_tillage`,
   `Z_tillage`, `I_tillage`, `Type_tillage`, `iType_Tillage`,
   `Rho_cons`, `Rho_tillage`, `k_R`, optional `Rho_match`/`N_match`)
   are loaded from `[soil.tillage]` in `swap.toml` via the typed
   pipeline.
2. `Read_Tillage` (in `src/crop/tillage.f90`) is deleted; an
   `apply_soil_tillage` helper in
   `src/io/toml/config_to_variables.f90` is the sole entry point
   for tillage-init.
3. `tillage.f90` no longer references `swpfile`, `logf`, `RDinit`,
   `RDsinr`, `RDfdor`, `RDftim`, `RDfinr`, or `RDinne`. The file's
   only TTutil dependency goes away.
4. Validators reject malformed `[soil.tillage]` blocks (empty
   events when `swtill=1`, dates outside simulation window,
   `type_id` not present in `types`, …). Errors are collected via
   `error_collection_t`, not `STOP`.
5. pFUnit suites cover schema parsing, validation, and adapter
   behaviour. check-full stays 5/5 (no swtill=1 regression case is
   added; that's deferred).
6. ADR 0021 captures the decision; `configuration-schema.md`
   documents the new block.

End-state: `swpfile` has only one remaining reader (SSDI), retired
in the sibling spec.

## Non-goals

- **SSDI TOML port.** `read_ssdi_input` remains as-is until the
  next spec. The `swpfile`/`logf` plumbing in
  `config_to_variables.f90` stays put.
- **Adding a tillage regression case.** All five existing cases
  default to `swtill=0` and are unaffected by this work. Authoring
  a sixth synthetic-tillage case is a separate, follow-on task —
  preferably after we have a known-good `swtill=1` path to compare
  the legacy binary against.
- **Refactoring `DoTillage(2)`/`DoTillage(3)` runtime logic.** This
  spec is parameter-load-only. The physics/state-update entry
  points stay unchanged.
- **Removing the `flCropNut`-stub-error guard at
  `tillage.f90:79`.** That belongs to the nutrient-reactivation
  arc.

## Schema

### New types in `src/config/soil_config.f90`

```fortran
type :: soil_tillage_event_t
   character(len=10) :: date              ! ISO YYYY-MM-DD
   real(8)           :: z          = 0.0d0  ! depth, cm (0 .. |zbotcp(NumNod)|)
   real(8)           :: intensity  = 0.0d0  ! 0..1
   integer           :: type_id    = 0      ! >= 1; references soil_tillage_type_t.id
end type

type :: soil_tillage_type_t
   integer :: id          = 0               ! 1..N, matches event.type_id
   real(8) :: rho_cons    = 0.0d0           ! 100..3000 (kg/m^3)
   real(8) :: rho_tillage = 0.0d0           ! 100..3000
   real(8) :: k_R         = 0.0d0           ! 1.0e-4..10.0
   real(8) :: rho_match   = -99.0d0         ! optional, used iff i_n_model=3
   real(8) :: N_match     = -99.0d0         ! optional, used iff i_n_model=3
end type

type :: soil_tillage_t
   integer :: i_n_model = 2                 ! 1..3
   integer :: iRedist   = 2                 ! 0..2
   type(soil_tillage_event_t), allocatable :: events(:)
   type(soil_tillage_type_t),  allocatable :: types(:)
end type
```

`soil_config_t` gains one new field:

```fortran
type(soil_tillage_t) :: tillage
```

`swtill` stays at `[soil].swtill` (its current location). The
`[soil.tillage]` block is **optional**: when `swtill = 0` the block
can be omitted entirely.

### TOML shape

```toml
[soil]
swtill = 1               # unchanged location

[soil.tillage]
i_n_model = 2            # 1..3
iRedist   = 2            # 0..2

[[soil.tillage.events]]
date      = 2003-04-15
z         = 30.0
intensity = 1.0
type_id   = 1

[[soil.tillage.events]]
date      = 2004-04-20
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

## Validation rules

In `soil_config_validate`, fired only when `self%swtill == 1`:

| Rule | Error context |
|---|---|
| `tillage.events` allocated and non-empty | `soil.tillage.events` |
| `tillage.types` allocated and non-empty | `soil.tillage.types` |
| `i_n_model` in [1, 3] | `soil.tillage.i_n_model` |
| `iRedist` in [0, 2] | `soil.tillage.iRedist` |
| `events(i).intensity` in [0, 1] | `soil.tillage.events[i].intensity` |
| `events(i).type_id >= 1` | `soil.tillage.events[i].type_id` |
| `events(i).type_id` ∈ {`types(j).id`} (set membership) | `soil.tillage.events[i].type_id` |
| `events` strictly date-ascending | `soil.tillage.events` |
| `events(i).date` ∈ [`simulation.start_date`, `simulation.end_date`] | `soil.tillage.events[i].date` |
| `types(i).rho_cons` in [100, 3000] | `soil.tillage.types[i].rho_cons` |
| `types(i).rho_tillage` in [100, 3000] | `soil.tillage.types[i].rho_tillage` |
| `types(i).k_R` in [1.0e-4, 10.0] | `soil.tillage.types[i].k_R` |
| `i_n_model == 3 ⇒ rho_match` in [100, 3000] for every type | `soil.tillage.types[i].rho_match` |
| `i_n_model == 3 ⇒ N_match` in [1.001, 10.0] for every type | `soil.tillage.types[i].N_match` |

**Deferred to adapter** (validator can't see `NumNod` / `zbotcp`):

| Rule | Where it fires |
|---|---|
| `events(i).z` in [0, `\|zbotcp(NumNod)\|`] | `apply_soil_tillage`, after grid is built |

All errors are appended to `error_collection_t`; the
`finalize`/abort checkpoint surfaces them at once (ADR 0008).

## Reader

Either extend `src/io/toml/read_soil_toml.f90` with a new
`read_soil_tillage_toml` helper, or create a sibling
`src/io/toml/read_soil_tillage_toml.f90`. Pattern follows
`read_meteorology_toml`'s array-of-tables walks for
`[[meteorology.rainfall_events]]`. The file structure question
(extend vs sibling) is decided in the implementation plan;
preference is sibling (keeps `read_soil_toml.f90` from growing).

The block is optional: if absent, the parser leaves
`config%soil%tillage` at its default-initialized state (empty
events/types, scalars at defaults).

## Adapter

`src/io/toml/config_to_variables.f90` gets a new helper:

```fortran
subroutine apply_soil_tillage(tillage)
   type(soil_tillage_t), intent(in) :: tillage
   ! ...
end subroutine
```

Called from the existing TOML→globals init pass when `flTillage`
is true (i.e., right after `flTillage = (config%soil%swtill == 1)`).
The call site replaces — and removes the need for — the
`call DoTillage(1)` ⇒ `Read_Tillage` chain at swap.f90:165.

`apply_soil_tillage` does:

1. Validates `events(i).z` against `|zbotcp(NumNod)|` (deferred
   rule). Errors → `error_collection_t`.
2. Allocates the legacy globals at the right size:
   - `Date_tillage(Ntill + 1)` (sentinel slot)
   - `Z_tillage(Ntill)`, `I_tillage(Ntill)`, `Type_tillage(Ntill)`
   - `iType_Tillage(Ntypes)`, `TAB_Rho_cons(Ntypes)`,
     `TAB_Rho_tillage(Ntypes)`, `TAB_K_R_cons(Ntypes)`
   - When `i_n_model == 3`: `TAB_Rho_match(Ntypes)`,
     `TAB_N_match(Ntypes)`
3. Copies values from the typed config; parses `events(i).date`
   into `t1900` representation via the same date helper used by
   `simulation.start_date`.
4. Sets the sentinel: `Date_tillage(Ntill + 1) = tend + 1.0d0`.
5. Computes `Max_Z_tillage = maxval(Z_tillage(1:Ntill))`.
6. Computes `iTT1`/`iTT2` first/last-position indices per current
   `Read_Tillage` body.
7. Sets module-level scalars: `i_n_model`, `iRedist`, `Ntill`,
   `Ntypes`.

## Runtime change in `src/crop/tillage.f90`

- **Delete `subroutine Read_Tillage`** in its entirety.
- **Remove the `Max_Z_tillage = 0.0d0` reset** at the top of
  whatever was the case-1 init body that called Read_Tillage —
  `apply_soil_tillage` writes `Max_Z_tillage` directly.
- **Drop the `if (swtill /= 1) return` self-checks** at lines 53
  and 72 (defense-in-depth per ADR 0020 commit, now redundant —
  `flTillage` is the single source of truth).
- The file's `use variables, only: ...` line drops `swpfile` and
  `logf` from its imports if they're no longer needed.
- The `flCropNut` stub-error guard at `tillage.f90:79` stays
  (separate concern, owned by the nutrient-reactivation arc).

After this lands, `grep -n "swpfile" src/crop/tillage.f90` returns
zero matches.

## Tests

Three new pFUnit suites under `tests/unit/io/toml/`:

| File | Asserts |
|---|---|
| `test_read_soil_tillage_toml.pf` | TOML parses correctly: scalars, events array, types array, optional rho_match/N_match when `i_n_model=3`; missing-block case (swtill=0) leaves defaults |
| `test_soil_tillage_validate.pf` | Validator catches every rule above: empty events/types when swtill=1, type_id not in types, dates non-ascending, dates outside simulation window, out-of-range scalars, missing rho_match/N_match when i_n_model=3 |
| `test_apply_soil_tillage.pf` | Adapter populates globals correctly: build minimal `swap_config_t` with two events + one type, set `NumNod`/`zbotcp` to satisfy depth check, call `apply_soil_tillage`, assert `Ntill==2`, `Ntypes==1`, `Date_tillage(3) == tend+1`, `Max_Z_tillage == maxval(Z_tillage)`, `iTT1`/`iTT2` correct |

`tests/unit/meson.build` lists the new `.pf` files. No changes to
existing parity suites (they don't exercise tillage).

`pixi run -e test check-full` is unchanged: 5/5 in ~42s; all
existing cases keep `swtill=0`.

## ADR

`docs/adr/0021-tillage-toml-port.md` (new):

- **Decision.** Tillage parameters move to `[soil.tillage]`
  sub-table. Legacy `swpfile`-based `Read_Tillage` retired.
- **Rationale.** ADR 0019/0020 closed the legacy-reader runtime
  path. Tillage was one of two remaining swpfile readers (the
  other is SSDI, deferred to its own ADR). Inline TOML over CSV:
  tables are small (1-10 events, 1-5 types), tightly coupled
  (events reference types by id), benefit from same-file
  readability.
- **Consequences.**
  - `swpfile` global pointer no longer read by `tillage.f90`.
  - `Read_Tillage` deleted; `apply_soil_tillage` is the sole
    tillage-init path.
  - Pattern reusable for SSDI (next ADR).
  - When `flCropNut` is reactivated (separate arc), the existing
    stub-error guard in `DoTillage` at `tillage.f90:79` will be
    revisited.

`docs/configuration-schema.md` gains a `[soil.tillage]` section
documenting fields + ranges.

## Acceptance criteria

- [ ] `pixi run -e test build-linux` clean.
- [ ] `pixi run -e test test-pfunit` exits 0; suite count grows
      by 3 new suites; existing 540 tests unchanged.
- [ ] `pixi run -e test check-full` exits 0 with `5 passed,
      0 failed`.
- [ ] `grep -n "subroutine Read_Tillage" src/` returns zero
      matches.
- [ ] `grep -n "swpfile" src/crop/tillage.f90` returns zero
      matches.
- [ ] A hand-authored `swap.toml` snippet under
      `tests/swap-cases/toml/` (or in the new
      `test_apply_soil_tillage.pf` fixture) exercises every
      validator rule above.
- [ ] ADR 0021 committed.
- [ ] `docs/configuration-schema.md` updated.

## Commit cadence

~5-6 commits, each independently buildable + green:

1. `feat(config): add soil_tillage_t types + validator`
2. `feat(io): parse [soil.tillage] in TOML reader`
3. `feat(io): apply_soil_tillage adapter helper`
4. `refactor(crop): retire Read_Tillage; tillage.f90 no longer reads swpfile`
5. `test: pFUnit suites for soil.tillage schema + adapter`
6. `docs: ADR 0021 + configuration-schema.md update`

The plan (writing-plans) will refine ordering — the validator
might need to land alongside (1) or (2), and (4) depends on (3)
working end-to-end.

## Risk

- **Date parsing.** `events(i).date` arrives as ISO string;
  `apply_soil_tillage` must convert to the `t1900` real(8) form
  used by `Date_tillage(:)`. The drainage / surface-water adapters
  already do this; reuse the same helper.
- **Empty `tillage` block when `swtill=0`.** TOML library must
  tolerate the missing `[soil.tillage]` table; verified by the
  default-init test.
- **`zbotcp` not set at validate-time.** Mitigation: split the
  z-range check into the adapter (deferred validation pattern,
  matches drainage/levels already in use).
