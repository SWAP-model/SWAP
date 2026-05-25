---
title: "ADR 0028 — [nutrients] N3 — Runtime Activation"
date: 2026-05-08
status: accepted
---

# ADR 0028: [nutrients] N3 — Runtime Activation

**Status:** accepted
**Date:** 2026-05-08
**Sub-arc of:** [nutrients] umbrella (N1 → N2a → N2b → **N3**)

## Context

After ADRs 0025/0026/0027, all nutrient typed-config flows to the
legacy globals via three adapters: apply_cropwofost_nutrient (per
rotation, N1), apply_nutrients (top-level [nutrients], N2a), and
apply_nutrients_events (CSV companion, N2b). The legacy global
flCropNut, which gates nutrient physics inside the daily loop and
inside DoTillage, was never assigned anywhere on the modern path —
it defaulted to .false. and stayed there. Two stub-errors held the
gate closed: a validator at cropwofost_init.f90:85 and a runtime
guard at tillage.f90:73.

## Decision

Drive flCropNut per-rotation from cfg%nutrient%flcropnut, inside
cropwofost_init_from_config. Remove both stub-errors. Defer the
regression fixture authorship — no legacy nutrient-enabled .crp
exists to mirror, and the parity baseline question is a separate
workstream.

Per-rotation assignment is correct because cropwofost_init_from_config
runs at every rotation start, before that rotation's first daily step.
A sequence of rotations with mixed flcropnut values therefore toggles
the global correctly. This mirrors the per-rotation crop config cache
pattern (ADR 0016).

An OR-reduction across all rotations was rejected — it would enable
nutrient gates for rotations whose typed config explicitly disabled
them.

## Consequences

- The [nutrients] umbrella is complete on the TOML path. A rotation
  with cropwofost.nutrient.flcropnut = true now runs the WOFOST
  nutrient subsystem end-to-end.
- No new dependencies introduced — cropwofost_init.f90 is
  TTutil-clean (modernized in N1) and the change is one assignment
  plus the existing N1 adapter call.
- Correctness of nutrient outputs is not validated by this arc. No
  legacy parity fixture exists; first regression case is a separate
  workstream.
- Stale "flCropNut=1 is stub-errored upstream" comments across
  swap.f90, management_soil.f90, cropwofost_init.f90 and
  cropgrowth.f90 have been refreshed.

## References

- ADR 0025 — N1 crop-side adapter
- ADR 0026 — N2a soil-side initial state
- ADR 0027 — N2b timed amendments
- ADR 0023 / 0024 — TTutil retirement context
- ADR 0016 — per-rotation crop config cache (the precedent for the
  per-rotation assignment pattern)
- Spec: docs/superpowers/specs/2026-05-08-nutrients-N3-runtime-activation-design.md
