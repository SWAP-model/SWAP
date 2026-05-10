---
title: "ADR 0033 — Cumulative Reset Cohorts"
date: 2026-05-10
status: accepted (Phase A complete; Phase B pending)
---

# ADR 0033: Cumulative Reset Cohorts

**Status:** accepted (Phases A + B + C complete)
**Date:** 2026-05-10
**Branch:** `refactor/cumulative-reset-cohorts`

## Context

After three subsystem state-migration arcs (surfacewater ADR 0030, drainage ADR 0031, solute ADR 0032), each subsystem state record carries a flat list of fields including its cumulative balance and intermediate accumulator variables. Reset logic for those variables lives in inline `if (flzerointr)` / `if (flzerocumu)` blocks scattered across multiple files (`surfacewater.f90`, `drainage.f90`, `solute.f90`). The migration arcs produced several incidents around these reset blocks:

- Forgotten resets when adding new cumulative fields.
- Duplicated reset blocks zeroing non-identical subsets across files (surfacewater + drainage; the asymmetry was real and caught only at the integration gate of Task A5).
- Silent drift between dual-write blocks during multi-task migrations.

The goal is to make reset cadence first-class in the type system: each subsystem state type carries cohort sub-records — one per reset cadence — and the only method on each cohort is `reset()`. Direct field writes for accumulation everywhere else; no setters, no auto-resets.

## Decision

Each subsystem state record gains two nested sub-records:

- `<subsystem>_intermediate_t` — fields that reset on `flzerointr`.
- `<subsystem>_cumulative_t` — fields that reset on `flzerocumu`.

Each cohort declares one type-bound procedure: `reset()`. The procedure zeroes every field it owns, with `allocated()` guards on per-level/per-node allocatable arrays.

Field paths grow one level deeper:
- Before: `state%surfacewater%cqdra`
- After:  `state%surfacewater%cumulative%cqdra`

Accumulation sites mutate fields directly: `state%X%cumulative%Y = state%X%cumulative%Y + delta`. No `add_*` methods. Special non-zero rebases (e.g. solute's `samini = sampro`, in Phase B) stay inline at the flag-gated call site, immediately after the corresponding `reset()` call, with an explanatory comment.

## Phase A outcomes (surfacewater, completed 2026-05-10)

- `src/state/surfacewater_state.f90` adds `surfacewater_intermediate_t` (4 fields: `iqdra`, `inqdra(:,:)`, `inqdra_in(:,:)`, `inqdra_out(:,:)`) and `surfacewater_cumulative_t` (7 fields: `cqdrd`, `cwsupp`, `cwout`, `cqdra`, `cqdrain(:)`, `cqdrainin(:)`, `cqdrainout(:)`).
- `surfacewater_state_t` carries them as nested components: `intermediate` and `cumulative`. The 11 fields are no longer flat in `surfacewater_state_t`.
- 113 reader/writer sites across `surfacewater.f90`, `drainage.f90`, `surfacewater_init.f90`, `waterbalance.f90`, `swapoutput.f90`, `swap_csv_output.f90`, `soilgrid.f90`, `management_soil.f90`, and `tests/unit/state/test_surfacewater_state.pf` migrated to the nested paths. ASSOCIATE block aliases (`sw_cqdrd`, etc.) keep the same local names but point at the new nested paths.
- 5 pFUnit tests for the two `reset()` procedures cover scalar zeroing, allocatable-array zeroing, and safe-when-not-allocated semantics.
- The 22-line inline reset block in `surfacewater.f90:SurfaceWater(task=2)` collapses into two `call state%surfacewater%X%reset()` invocations.
- check-full byte-identical at every commit; pFUnit 646 tests passing.

### Cumulative-cohort ownership asymmetry (Task A5 finding)

The `flzerocumu` reset block at `drainage.f90:Drainage()` zeros a STRICT SUBSET of the cumulative cohort: only `cqdra`, `cqdrain(:)`, `cqdrainin(:)`, `cqdrainout(:)` (4 of 7 cohort fields). The other 3 (`cqdrd`, `cwsupp`, `cwout`) are accumulated by `SurfaceWater(task=2)`, which runs AFTER `Drainage()` in `swap_main`. If drainage called the full cohort `reset()`, it would silently zero `cqdrd`/`cwsupp`/`cwout` mid-accumulation, corrupting the surface-water balance calculation downstream.

**Decision:** hybrid pattern.

- `surfacewater.f90:SurfaceWater(task=2)` calls the full `cumulative%reset()` — the canonical owner of the full cohort.
- `drainage.f90:Drainage()` keeps element-by-element zeroing of its 4-field subset, with an `allocated()` guard, with a comment pointing to this ADR section. It does NOT call the type-bound `reset()`.
- `intermediate%reset()` is symmetric — both sites zero the identical 4-field set, so both call the full cohort reset.

**Lesson:** when a cohort field is mid-accumulation across multiple subsystem entry points, the cohort `reset()` must be called from the LATEST point in the call chain that the cohort is touched, not earlier sites. The pattern doesn't require unanimous full-reset; subset zeroing at intermediate sites is legitimate and the type system can't enforce it. ADR readers should treat "every reset block becomes a cohort reset() call" as a heuristic, not a rule.

## Consequences

(+) Reset logic is co-located with the data; pFUnit can drive a populated cohort through `reset()` without SWAP scaffolding (see `tests/unit/state/test_surfacewater_cumulative.pf`).
(+) Adding a new cumulative field is one line in the cohort type; no need to remember to zero it in N reset blocks.
(+) Reset cadence is now part of the type signature — readers see `intermediate` vs `cumulative` and immediately know reset semantics.
(+) Field paths self-document: `state%X%cumulative%Y` says "this is a cumulative balance variable in subsystem X."
(–) Field paths grow one level deeper. ASSOCIATE blocks mitigate this for compute-heavy bodies.
(–) Asymmetric subset-resets (the Task A5 finding) cannot be enforced by the cohort type; they remain a per-call-site policy.

## Phase B outcomes (solute, completed 2026-05-10)

- `src/state/solute_state.f90` adds `solute_intermediate_t` (6 fields: `imsqprec`, `imsqirrig`, `imsqbot`, `imsqdra`, `imdectot`, `imrottot`) and `solute_cumulative_t` (9 fields: `sqprec`, `sqirrig`, `sqbot`, `sqdra`, `sqsur`, `dectot`, `rottot`, `csurf`, `samini`).
- `solute_state_t` carries them as nested `intermediate` and `cumulative` components. The 15 fields are no longer flat in `solute_state_t`.
- 31 reader/writer sites across `solute.f90`, `swapoutput.f90`, `swap_csv_output.f90`, and the unit tests migrated to nested paths. Two ASSOCIATE blocks in solute.f90 (case 1 and case 2) keep their local alias names but point at the new nested paths. Two additional ASSOCIATE blocks in swapoutput.f90 (`sl => state%solute`) discovered during compile (initial grep missed them because they alias `sl%`, not `state%solute%`) — both retargeted.
- 4 pFUnit tests added for the two `reset()` procedures (no allocatables in solute cohorts, so 4 instead of 5 — the allocation-guard test is N/A).
- The 15-line inline reset block in `solute(task=2)` collapses into two `call state%solute%X%reset()` invocations + an inline `samini = sampro` rebase preserved at the call site. The rebase resolves through the case(2) ASSOCIATE alias to `state%solute%cumulative%samini`. Comment at the call site explains that the rebase is physics (mass-balance anchoring), not cohort policy.
- `isqbot`, `isqtop`, `isqdra` zeroing remains inline inside `solute(task=2)` — these are instantaneous per-step fluxes (reset every step unconditionally, not flag-gated). They are owned by `solute_state_t` directly, not by a cohort.
- check-full byte-identical at every commit; pFUnit 650 tests passing.

### Phase B finding: solute has NO ownership-asymmetry hazard

Unlike Phase A's surfacewater/drainage finding, solute's reset block is the sole writer of all 15 cohort fields — there is no second site that zeros a strict subset. The full cohort `reset()` is a clean drop-in.

## Migration playbook update (Phase C — see ADR 0033 + plan)

The cohort pattern is now part of the migration playbook. Future subsystem migration discoveries classify owned fields by reset cadence (instantaneous / intermediate / cumulative) upfront, and design specs pre-define cohort sub-records as part of the state-type design. See discovery-template Section 2 update and design-template addendum (committed in CRR Phase C).

## Consequences (final)

(+) Reset logic is co-located with the data; pFUnit can drive a populated cohort through `reset()` without SWAP scaffolding.
(+) Adding a new cumulative field is one line in the cohort type; no need to remember to zero it in N reset blocks.
(+) Reset cadence is now part of the type signature — readers see `intermediate` vs `cumulative` and immediately know reset semantics.
(+) Field paths self-document: `state%X%cumulative%Y` says "this is a cumulative balance variable in subsystem X."
(–) Field paths grow one level deeper. ASSOCIATE blocks mitigate this for compute-heavy bodies.
(–) Asymmetric subset-resets (the Task A5 finding for surfacewater/drainage) cannot be enforced by the cohort type; they remain a per-call-site policy. Solute has no asymmetric reset, so the Phase A pattern works clean there.
(–) Two additional alias forms (`sl => state%solute`, `sw => state%surfacewater`) require greps tailored beyond `state%X%` — Phase B Task B3 caught this when initial inventory missed two ASSOCIATE blocks. Future migrations should grep for ALL alias forms during the audit step.

## References

- Plan: `docs/superpowers/plans/2026-05-10-cumulative-reset-refactor.md`
- Phase A commit chain: `50ad0a2` → `55fa08f` → `45d4514` → `e43fac3` → `e5b8f7d` → `fba1288`
- Phase B commit chain: `8208382` → `d4c0cb5` → `ed30706` → `fea3cf4`
- Predecessor state migrations: ADR 0030 (surfacewater), ADR 0031 (drainage), ADR 0032 (solute).
