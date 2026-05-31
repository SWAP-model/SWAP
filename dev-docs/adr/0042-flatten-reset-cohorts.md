---
title: "ADR 0042 — Flatten Reset Cohorts Across State Subsystems"
date: 2026-05-12
status: accepted
supersedes-portion-of: ADR 0033
---

# ADR 0042: Flatten Reset Cohorts Across State Subsystems

**Status:** accepted
**Date:** 2026-05-12
**Branch:** `development`

## Context

ADR 0033 introduced nested cohort sub-records (`<subsystem>_intermediate_t`,
`<subsystem>_cumulative_t`, and surfacewater's partitioned cumulative
cohorts) on three subsystem state records — surfacewater, solute, soilwater
— to make reset cadence first-class in the type system. The cohort pattern
produced two real benefits: (a) co-location of fields that zero together,
(b) type-bound `reset()` syntactically compact.

In review, the wrapper bought less than expected:

- Activity gates (`flzerointr`, `flzerocumu`, `flDayStart`, `fldrain`,
  `flSurfaceWater`) remain enforced at call sites; the cohort type does
  not enforce them.
- Field-name prefixes (`i*` / `c*`) already encode cohort role; the
  wrapper is redundant with the naming convention.
- The cohort-vs-flat boundary on the parent is inconsistent: scalar
  non-reset fields (`wls`, `imper`) live flat while reset fields hide
  one level deeper.
- Access paths in hot compute loops (waterbalance, soilgrid) traverse
  three levels.
- Allocation traverses cohort paths separately from flat-field paths,
  fragmenting an otherwise unified init step.

The asymmetric reset semantics that ADR 0033 wanted to capture
(surfacewater's partitioned cumulative; soilwater's per-day subset
reset) can be expressed equally well by **named reset procedures on the
flat parent type**. Each asymmetry becomes its own procedure; the
procedure name carries the gate. The cohort form was, in effect, "two
reset procedures wearing a type-shaped costume" (surfacewater) or "a
second procedure on one cohort" (soilwater).

## Decision

Flatten the cohort sub-records into the parent state record. Move
`reset()` from cohort-type-bound to parent-type-bound procedures, named
to preserve the gate/cadence semantics.

Per subsystem:

- `surfacewater_state_t`:
  - `reset_intermediate()` (flzerointr)
  - `reset_cumulative_drainage()` (flzerocumu + fldrain partition; owner: drainage)
  - `reset_cumulative_reservoir()` (flzerocumu + flSurfaceWater partition; owner: surface-water)

- `solute_state_t`:
  - `reset_intermediate()` (flzerointr)
  - `reset_cumulative()` (flzerocumu)

- `soilwater_state_t`:
  - `reset_intermediate()` (flzerointr; subsumes per-day subset)
  - `reset_intermediate_per_day()` (flDayStart; per-day subset only)
  - `reset_cumulative()` (flzerocumu)

Section-header comments inside the flat type group fields by cadence
(`! === intermediate (reset_intermediate / gate: flzerointr) ===`,
etc.), restoring the type-signature self-documentation that ADR 0033
§82-83 valued.

## Implementation

Three per-subsystem commits in order solute → soilwater → surfacewater.
Each commit fully verified (build + pFUnit + check-full + audit grep)
before the next.

- Solute commit: `5c65fae` (Task 1)
- Soilwater commit: `dcd5b90` (Task 2)
- Surfacewater commit: `46ba744` (Task 3)

The cohort sub-record types are removed from the public surface of each
state module. All call-site paths collapse one level shallower
(`state%X%cohort%field` → `state%X%field`). Reset call sites rewrite to
the new procedure names. ASSOCIATE block RHS targets rewrite mechanically
under the same path-collapse rule.

## Out of scope

- **Allocation consolidation.** Duplicate inline allocation blocks
  (`surfacewater_init.f90:119-141` and `drainage.f90:445-467`, and the
  fragmented allocator paths in `soilwater_state.f90:soilwater_init`)
  are pre-existing. Deferred to a separate follow-up arc.
- ADR 0033 itself is left intact as the historical record; this ADR
  references and supersedes only the cohort-pattern decision.

## Consequences

(+) Call-site paths one level shallower; less friction in hot compute
    loops.
(+) Subset / asymmetric resets are explicit named procedures rather
    than conventions buried inside cohort types.
(+) Removes the cohort-vs-flat inconsistency on the parent (some
    fields nested, some not).
(+) Field-name prefix convention (`i*` / `c*`) is no longer redundant
    with cohort wrappers.

(–) Adding a new reset-eligible field requires editing the field
    declaration and the appropriate reset procedure body — same
    cognitive load as before, but the procedure body is no longer
    co-located with the field inside a sub-type. Mitigation:
    section-header comments inside the type group fields by cadence.
(–) Type-bound `cohort%reset()` syntactic compactness is lost. The
    parent now carries multiple named reset procedures; names are
    slightly longer but spell out the gate.
(–) The ADR 0033 §82-83 self-documentation property — `state%X%cumulative%Y`
    instantly tells the reader "this is a cumulative balance variable"
    — weakens. A grep against `state%X%fieldname` alone no longer
    reveals cadence. Section-header comments partially restore this.

## References

- ADR 0030 — state-migration surfacewater pilot
- ADR 0032 — state-migration solute
- ADR 0033 — cumulative reset cohorts (the pattern this arc walks back)
- ADR 0038 — state-migration soilwater core
