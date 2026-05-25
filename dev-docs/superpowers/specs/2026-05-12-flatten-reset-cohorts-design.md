---
title: "Flatten Reset Cohorts — Design Spec"
date: 2026-05-12
status: draft
supersedes-portion-of: ADR 0033
---

# Flatten Reset Cohorts Across State Subsystems

## Context

Three subsystem state records — `surfacewater_state_t`, `solute_state_t`, `soilwater_state_t` — currently group their flag-gated reset fields into nested cohort sub-records (`<subsystem>_intermediate_t`, `<subsystem>_cumulative_t`, plus surfacewater's two-way cumulative split). The cohort pattern was introduced by ADR 0033 to make reset cadence first-class in the type system and to express asymmetric reset semantics (different gates, different subsets, different owners) as types rather than comments.

In review, the cohort wrapper buys less than expected:

- Gates remain enforced at call sites (`if (flzerointr) call ...reset()`); the type does not enforce them.
- Field-name prefixes (`i*` = intermediate, `c*` = cumulative) already encode cohort role; the wrapper is redundant with the naming convention.
- The cohort boundary is inconsistent: scalar non-reset fields (`wls`, `swst`, `imper`) live flat on the parent while reset fields hide one level deeper.
- Access paths in hot compute paths (waterbalance, soilgrid) traverse three levels (`state%X%cohort%field`).
- Allocation traverses cohort paths (e.g. `state%surfacewater%intermediate%inqdra` and `state%surfacewater%drainage_cumulative%cqdrain` are allocated as if they were two unrelated arrays).

The asymmetric reset semantics that ADR 0033 wanted to capture (surfacewater's partitioned cumulative; soilwater's per-day subset reset) can be expressed equally well by **named reset procedures on the flat parent type**. Each asymmetry becomes its own procedure; the procedure name carries the gate. That is more honest than the cohort form: surfacewater's two cumulative cohorts are already "two reset procedures wearing a type-shaped costume", and soilwater's `reset_per_day` is already a second procedure on one cohort.

## Decision

Flatten the cohort sub-records into the parent state record. Move `reset()` from cohort-type-bound to parent-type-bound procedures, named to preserve the gate/cadence semantics that the cohort type names used to carry. All three subsystems migrate in one arc.

## Architecture

### What goes

- `surfacewater_intermediate_t`, `surfacewater_drainage_cumulative_t`, `surfacewater_reservoir_cumulative_t` — deleted.
- `solute_intermediate_t`, `solute_cumulative_t` — deleted.
- `soilwater_intermediate_t`, `soilwater_cumulative_t` — deleted.
- Nested components (`%intermediate`, `%cumulative`, `%drainage_cumulative`, `%reservoir_cumulative`, `%intr`, `%cumu`) — removed from the parent types.

### What stays

- All field names, types, shapes, semantics — unchanged. Only the path changes.
- Call-site gating (`if (flzerointr)`, `if (flzerocumu)`, `if (flDayStart)`) — unchanged, still in caller's hands.
- Activity gates (`fldrain`, `flSurfaceWater`) — unchanged, still enforced by which call site invokes which procedure.
- Allocation lifecycle (`surfacewater_init.f90` and inline allocators in `drainage.f90`) — unchanged in this arc. Cleanup deferred to a follow-up.
- Instantaneous-flux inline zeroing (`isqbot`, `isqtop`, `isqdra` in `solute(task=2)`) — unchanged.
- `samini = sampro` rebase in `solute(2)` — unchanged.
- ASSOCIATE blocks at call sites — preserved; only the targets in `ASSOCIATE (alias => target, ...)` headers are rewritten one level shallower.

## Per-subsystem reset interface

All reset procedures are type-bound on the parent state record. Procedure bodies are identical to today's cohort `reset()` bodies, just relocated. Allocatable arrays remain `if (allocated(...))`-guarded.

### `surfacewater_state_t`

| Procedure | Fields zeroed | Call sites | Gate |
|---|---|---|---|
| `reset_intermediate()` | `iqdra`, `inqdra`, `inqdra_in`, `inqdra_out` | `Drainage(2)`, `SurfaceWater(2)` | `flzerointr` |
| `reset_cumulative_drainage()` | `cqdra`, `cqdrain`, `cqdrainin`, `cqdrainout` | `Drainage(2)`, `SurfaceWater(2)` | `flzerocumu` (owner: drainage; fields accumulate under `fldrain`) |
| `reset_cumulative_reservoir()` | `cqdrd`, `cwsupp`, `cwout` | `SurfaceWater(2)` only | `flzerocumu` (owner: surface-water; fields accumulate under `flSurfaceWater`) |

### `solute_state_t`

| Procedure | Fields zeroed | Call site | Gate |
|---|---|---|---|
| `reset_intermediate()` | `imsqprec`, `imsqirrig`, `imsqbot`, `imsqdra`, `imdectot`, `imrottot` | `solute(2)` | `flzerointr` |
| `reset_cumulative()` | `sqprec`, `sqirrig`, `sqbot`, `sqdra`, `sqsur`, `dectot`, `rottot`, `csurf`, `samini` | `solute(2)` | `flzerocumu` |

### `soilwater_state_t`

| Procedure | Fields zeroed | Gate |
|---|---|---|
| `reset_intermediate()` | All intermediate fields (main scalars + per-node arrays + per-day subset) | `flzerointr` (subsumes `flDayStart` effect) |
| `reset_intermediate_per_day()` | Per-day subset only (8 fields: `tra`, `i*_day` scalars, `qpotrot_day`, `qredtot_day`) | `flDayStart` (when `flzerointr` is false) |
| `reset_cumulative()` | 14 cumulative fields | `flzerocumu` |

The full `reset_intermediate()` explicitly zeros the per-day subset as part of the full zero. It does not delegate to `reset_intermediate_per_day()` — matches current behavior and avoids procedure-call indirection. Exact intermediate field cardinality to be confirmed during plan-writing against the current `soilwater_intermediate_reset` body.

### Naming rule

`reset_<cadence>[_<partition>]` where `<cadence>` is `intermediate` or `cumulative` and `<partition>` is an activity-gate qualifier when one exists (`drainage`, `reservoir`, `per_day`).

## Field grouping inside the flat type

To preserve the type-signature self-documentation that ADR 0033 valued, each flat type groups its fields with section-header comments by reset cadence/gate. The grouping is a reading aid only — no runtime structure.

Illustrative shape of `surfacewater_state_t`:

```fortran
type :: surfacewater_state_t
   ! === per-step / per-day scalars (no flag-gated reset) ===
   real(real64) :: wls           = 0.0_real64
   real(real64) :: wlstar        = 0.0_real64
   real(real64) :: swst          = 0.0_real64
   real(real64) :: swstini       = 0.0_real64
   real(real64) :: hwlman        = 0.0_real64
   real(real64) :: vtair         = 0.0_real64
   real(real64) :: wlsold        = 0.0_real64
   real(real64) :: ZDraBas       = 0.0_real64
   real(real64) :: qdrtot        = 0.0_real64
   logical      :: overfl        = .false.
   logical      :: flInitDraBas  = .true.
   integer      :: imper         = 1
   integer      :: numadj        = 0
   real(real64) :: wlsbak(4)     = 0.0_real64
   real(real64) :: sttab(22, 2)  = 0.0_real64

   ! === intermediate (reset_intermediate / gate: flzerointr) ===
   real(real64) :: iqdra = 0.0_real64
   real(real64), allocatable :: inqdra(:,:)
   real(real64), allocatable :: inqdra_in(:,:)
   real(real64), allocatable :: inqdra_out(:,:)

   ! === cumulative — drainage subsystem
   !     (reset_cumulative_drainage / gate: flzerocumu + fldrain) ===
   real(real64) :: cqdra = 0.0_real64
   real(real64), allocatable :: cqdrain(:)
   real(real64), allocatable :: cqdrainin(:)
   real(real64), allocatable :: cqdrainout(:)

   ! === cumulative — reservoir subsystem
   !     (reset_cumulative_reservoir / gate: flzerocumu + flSurfaceWater) ===
   real(real64) :: cqdrd  = 0.0_real64
   real(real64) :: cwsupp = 0.0_real64
   real(real64) :: cwout  = 0.0_real64
contains
   procedure :: reset_intermediate          => surfacewater_reset_intermediate
   procedure :: reset_cumulative_drainage   => surfacewater_reset_cumulative_drainage
   procedure :: reset_cumulative_reservoir  => surfacewater_reset_cumulative_reservoir
end type surfacewater_state_t
```

Same convention applies to `solute_state_t` and `soilwater_state_t`. Each section header names: (a) which reset procedure zeros it, (b) which flag(s) gate it.

## Call-site migration

### Path rewrites (mechanical)

| Before | After |
|---|---|
| `state%surfacewater%intermediate%iqdra` | `state%surfacewater%iqdra` |
| `state%surfacewater%intermediate%inqdra` | `state%surfacewater%inqdra` |
| `state%surfacewater%drainage_cumulative%cqdra` | `state%surfacewater%cqdra` |
| `state%surfacewater%drainage_cumulative%cqdrain` | `state%surfacewater%cqdrain` |
| `state%surfacewater%reservoir_cumulative%cqdrd` | `state%surfacewater%cqdrd` |
| `state%solute%intermediate%imsqprec` | `state%solute%imsqprec` |
| `state%solute%cumulative%sqprec` | `state%solute%sqprec` |
| `state%soilwater%intr%X` | `state%soilwater%X` |
| `state%soilwater%cumu%X` | `state%soilwater%X` |

### Reset-call rewrites

```fortran
! before (surfacewater example)
if (flzerointr) call state%surfacewater%intermediate%reset()
if (flzerocumu) then
   call state%surfacewater%drainage_cumulative%reset()
   call state%surfacewater%reservoir_cumulative%reset()
end if

! after
if (flzerointr) call state%surfacewater%reset_intermediate()
if (flzerocumu) then
   call state%surfacewater%reset_cumulative_drainage()
   call state%surfacewater%reset_cumulative_reservoir()
end if
```

### ASSOCIATE blocks

Preserved. Aliases (`sw_cqdrd`, `sw_iqdra`, `sl => state%solute`, `sw => state%surfacewater`) all still resolve; only the target paths in the `ASSOCIATE (alias => target, ...)` headers move one level shallower. No call-site body changes.

### Audit greps

Each subsystem migration is gated by zero hits on:

- surfacewater:
  - `grep -rn 'state%surfacewater%intermediate%\|state%surfacewater%drainage_cumulative%\|state%surfacewater%reservoir_cumulative%' src/ tests/`
  - `grep -rn 'sw%intermediate%\|sw%drainage_cumulative%\|sw%reservoir_cumulative%' src/ tests/` (alias `sw => state%surfacewater` exists in `surfacewater_init.f90`)
- solute:
  - `grep -rn 'state%solute%intermediate%\|state%solute%cumulative%' src/ tests/`
  - Per-field ASSOCIATE aliases (e.g. `samini => state%solute%cumulative%samini`) are caught by the `state%solute%cumulative%` grep on their RHS — no separate alias-form grep required today. If a block alias `sl => state%solute` is introduced later, add `sl%intermediate%` / `sl%cumulative%` to the audit.
- soilwater:
  - `grep -rn 'state%soilwater%intr%\|state%soilwater%cumu%' src/ tests/`
  - Per-field ASSOCIATE aliases (e.g. `sw_inq => state%soilwater%intr%inq`) are caught by the `state%soilwater%intr%` grep on their RHS — no separate alias-form grep required today.

Per the ADR 0033 §113 lesson, all ASSOCIATE alias RHS targets must be rewritten alongside direct accesses. The greps above operate on RHS path patterns, which catches both direct uses and per-field ASSOCIATE bindings.

### Estimated call-site counts

- surfacewater: ~70-90 (ADR 0030/0033 reported 113 before bridge-sync cleanups already landed)
- solute: ~31 (ADR 0033 §91)
- soilwater: to be inventoried during plan-writing

## Tests

### Existing pFUnit coverage

- `tests/unit/state/test_surfacewater_state.pf` — `intermediate%reset()` (scalar zero, allocatable zero, safe-when-unallocated)
- `tests/unit/state/test_surfacewater_cumulative.pf` — `drainage_cumulative%reset()`, `reservoir_cumulative%reset()`
- Solute and soilwater equivalents — 4 + N tests per ADR 0033

### Migration plan

- Update test bodies to invoke `call state%X%reset_<cadence>()` against the parent type instead of `call state%X%<cohort>%reset()`.
- Same field assertions, same shape, same test count.
- Drop any `use` imports of the deleted cohort type names.

### `test_surfacewater_parity.pf`

The TOML adapter parity test follows the same path-rewrite rules as production code. No semantic test changes.

### No new tests required

Reset semantics are unchanged; the procedure binding moves from sub-type to parent-type. The existing tests prove the same field-zeroing contract.

## Verification gates

Each subsystem migration commit is gated by:

1. **Build clean** — `meson compile -C builddir`, no new warnings.
2. **pFUnit green** — full unit-test suite passes (~650+ tests currently).
3. **check-full byte-identical** — regression suite produces bit-identical CSV/binary output to the pre-refactor baseline.
4. **Audit grep zero** — the audit greps for the migrated subsystem return 0 hits.

Per the project workflow convention, check-full must be verified **before** committing each subsystem migration. pFUnit alone misses global-default regressions.

## Sequencing

All three subsystems land on one branch / PR, with per-subsystem commits for reviewability. Commit order on that branch is:

1. **solute** — simplest (no asymmetry: one intermediate cohort + one cumulative cohort, no allocatables in cohorts).
2. **soilwater** — one asymmetry (per-day subset reset within intermediate).
3. **surfacewater** — two asymmetries (split cumulative by activity gate).

Each step's audit-grep and verification gate must pass before moving to the next. The per-subsystem commits keep the diff reviewable.

## Scope boundaries

### In scope

- Flatten cohort sub-records in surfacewater_state, solute_state, soilwater_state.
- Move `reset()` procedures to parent type with `reset_<cadence>[_<partition>]` naming.
- Rewrite call-site paths and ASSOCIATE block targets.
- Update existing pFUnit tests to new procedure paths.
- Write a new ADR (next-numbered) recording this arc and referencing/superseding the cohort-pattern portion of ADR 0033.

### Out of scope

- **Allocation consolidation.** The duplication between `surfacewater_init.f90:119-141` and `drainage.f90:445-467` (and equivalent in other subsystems) is pre-existing and visible. Deferred to a separate small follow-up arc immediately after this one.
- Any other state subsystem (drainage, heat, atmosphere, tillage, timecontrol, crop-uptake, boundary). Their state types didn't adopt the cohort pattern.
- Field renaming, type changes, or reshuffling fields between subsystems.
- Editing ADR 0033 itself. The historical decision document remains intact; the new ADR references it.

## Consequences

(+) Reset cadence still legible from the type — section-header comments and procedure names spell out the gate.
(+) Call-site paths one level shallower; less friction in hot compute loops (waterbalance, soilgrid).
(+) Subset/asymmetric resets are explicit named procedures instead of conventions buried inside a cohort.
(+) Removes the cohort-vs-flat inconsistency on the parent (some fields nested, some not).
(+) Field naming convention (`i*` / `c*` prefixes) is no longer redundant with cohort wrappers.

(–) Adding a new reset-eligible field requires a new field declaration *and* a line in the appropriate reset procedure — same cognitive load as before, but the procedure body is no longer co-located with the field in the type definition. Mitigation: the section-header comments make the grouping visible at the type definition.
(–) Type-bound `reset()` syntactic compactness (`call cohort%reset()`) is lost; the parent now carries multiple named reset procedures. Names are slightly longer but spell out the gate.
(–) ADR 0033 §82-83 self-documentation property weakens: a reader of a single field-access line no longer sees the cadence in the path. Section-header comments in the type definition partially restore this, but a grep against `state%X%fieldname` alone no longer reveals cadence.

## References

- ADR 0030 — state-migration surfacewater pilot
- ADR 0032 — state-migration solute
- ADR 0033 — cumulative reset cohorts (the pattern this arc walks back)
- ADR 0038 — state-migration soilwater core
- `src/state/surfacewater_state.f90`, `src/state/solute_state.f90`, `src/state/soilwater_state.f90` — current cohort definitions
