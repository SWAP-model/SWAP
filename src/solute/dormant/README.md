# src/solute/dormant — dormant solute features

Solute-side features whose **legacy bodies are preserved** but which
have **no live dispatch site** in the TOML-pipeline build. Files in
this directory are deliberately **excluded from the meson build**:
they exist as source-controlled reference for future reactivation.

See `src/soil/dormant/README.md` for the rationale (kept short here
to avoid duplication): in brief, dormant modules own their own
imports so retiring a bare global from `variables.f90` no longer
needs to consider dead-routine references, and the reactivation
checklist lives next to the body.

## Current dormant modules

| File | Feature | Reactivation prerequisite |
|------|---------|----------------------------|
| `agetracer.f90` | `AgeTracer` — groundwater-age tracking after Goode (1996), "Direct simulation of groundwater age" (WRR vol. 32, p. 289-296). Ageing equation solved alongside the main solute transport step. | Define `agetracer_state_t` (host the 12 Age* globals — Ageirr, Agedrain, Agepre, Agepond, Agepondm1, Agegwl1m, icAgetopupw, icAgetopdwn, icAgeBot, icAgeDra, icAgeRot, icAgeSur). Resolve the `cml` dual-use hazard (AgeTracer overwrites `state%solute%cml` at the end of task=2). Add a `flAgeTracer` config gate (currently never assigned anywhere). Restore the dispatch site (legacy: `call AgeTracer(1, state)` from solute init, `call AgeTracer(2, state)` from solute step). Update `outage`'s flAgeTracer guard in `swapoutput.f90`. |

## Build exclusion

These files are **not listed** in the top-level `meson.build` sources
list. They will not be compiled. To reactivate:

1. Work through the reactivation prerequisite for that file (above).
2. Add the file path to `meson.build`.
3. Restore the `use <module>_mod, only: <routines>` line wherever the
   dispatch site lives (for AgeTracer: `src/core/swap_mod.f90` had
   two `use agetracer_mod, only: AgeTracer` imports at lines 112 and
   487 — both were dead since ADR 0032 retired `flAgeTracer`).
4. Run `pixi run check-fast` — broken imports will surface
   immediately because the dormant body's `use variables, only:`
   lines reference globals that may have been retired in the
   meantime.

## Discovery hazards documented at extraction time

- **`cml` dual-use**: `state%solute%cml(:)` is the solute compute
  array. AgeTracer's task=2 finalization overwrites it with the
  age-tracer ageing concentrations so that the standard output
  channel picks them up. This creates a hidden coupling between
  the two computations. On reactivation, either add a separate
  `state%agetracer%cml_age(:)` field, or formalize the sequencing
  rule (and document it loudly).

- **Stale globals**: 22 bare globals (the Age* family + cml dual-use
  partners) were retired by commit `f1115bf` (ADR 0032) — they are
  no longer in `variables.f90`. The dormant body below still
  references them, so it WILL NOT COMPILE in the current tree.
  This is by design: reactivation must re-home them on `state` (see
  ADR 0032 + the 2026-05-10 solute migration discovery spec).
