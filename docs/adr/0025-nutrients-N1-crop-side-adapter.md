---
title: "ADR 0025 — [nutrients] N1: crop-side nutrient adapter"
date: 2026-05-07
status: accepted
---

# ADR 0025: [nutrients] N1 — crop-side nutrient adapter

## Context

The legacy SWAP nutrient subsystem reads ~25 crop-side
parameters (`LRNR`, `LSNR`, `NLAI`, `NMXLV(30)`,
`FraHarLosOrm_*`, …) from the per-rotation `<cropfil>.crp` via
TTutil's `rdinit`/`rdsdou` block inside `cropgrowth.f90`'s
`wofost` subroutine. That block was deleted in the legacy-readers
physical-deletion arc (SS-C step 2, commit dddbada) because
`flCropNut=1` is stub-errored upstream at `tillage.f90:73` and
the block was unreachable in any modern flow.

A typed-config schema for these parameters already exists:
`wofost_nutrient_t` in `src/config/cropwofost_config.f90`
(populated by the cropwofost TOML parser, commit 230751c). What
was missing: the validator stub-errored on `flcropnut=true`, and
no adapter wired the typed config to the legacy `variables`
module globals that the WOFOST physics body reads.

## Decision

This ADR captures the first sub-arc (N1) of the `[nutrients]`
umbrella. Three sub-arcs total: N1 wires the crop-side
parameters; N2 ports the soil-side initial state; N3 lifts the
runtime gate at `tillage.f90:73` and adds a regression case.

N1 specifically:

1. **Replace `wofost_nutrient_validate`'s blanket stub-error**
   with input validation rules: `nmxlv` non-empty,
   `size(nmxlv) <= 30`, `frahar_los_orm_*` in `[0, 1]`.
   The 17 numeric scalars stay unranged in N1 (matches legacy
   `rdsdou` which also didn't enforce ranges; future tightening
   can land when nutrient regression fixtures exist).
2. **Promote 15 nutrient names from local-to-`wofost`** to
   module-level in `variables.f90` (`NLAI`, `NMAXSO`, `NPART`,
   `NFIXF`, `NSLA`, `RNFRT`, `TCNT`, `DVSNLT`, `DVSNT`,
   `RDRNS`, `FNTRT`, `ILNMXL`, `FraHarLosOrm_lv/st/so`). The
   other 7 nutrient scalars (`LRNR`, `LSNR`, `NLUE`, `RNFLV`,
   `RNFST`, `FRNX`, `NMXLV(30)`) were already module-level.
3. **Add `apply_cropwofost_nutrient(cfg)`** in
   `cropwofost_init.f90` that copies all 22 module-level
   scalars + the `NMXLV` array + `ILNMXL` from the typed config.
4. **Call the adapter** from `cropwofost_init_from_config` per
   rotation when `cfg%nutrient%flcropnut = .true.`.

Implementation also flattened `wofost_nutrient_t%nmxlv` from
`(:,:)` to `(:)` to match the legacy `NMXLV(30)` layout, and
added `frahar_los_orm_lv/st/so` fields to the typed config that
were missing pre-N1. The TOML reader was updated to match.

## What N1 does NOT do

- Lift the runtime stub-error at `tillage.f90:73`. With no
  soil-side initial state wired (N2's job), allowing
  `flCropNut=1` to reach `DoTillage` would let
  `SoilManagement(2..7)` read uninitialised soil-pool state.
  N3 lifts the gate after N2 lands.
- Add a regression case with `flcropnut=true`. Cannot run
  end-to-end until the runtime gate is lifted (N3).
- Range-check the 17 numeric nutrient scalars. The legacy
  reader didn't either; future tightening when fixtures exist.

## Consequences

- A TOML config with `crop.rotation.<n>.cropwofost.nutrient.flcropnut = true`
  loads and validates without error (no stub).
- The `wofost` subroutine reads nutrient parameters from module
  variables that the per-rotation init has populated from the
  typed config. Behaviourally identical to the legacy
  rdinit/rdsdou flow for any populated config.
- All five existing regression cases keep `flcropnut=false`;
  check-full output is byte-identical to pre-arc baseline.
- Promoting the locals to module variables creates an
  architectural consistency: all nutrient parameters live in the
  same place, accessible to anyone with `use variables`. Future
  ADR (per the architectural direction in ADR 0024) may pass
  these as explicit arguments to `wofost`, retiring the module-
  global pattern entirely.

## Acceptance

- `grep -n "not yet supported in the TOML pipeline.*nutrient" src/config/cropwofost_config.f90` → no match.
- `grep -nE "^\s*real\(8\)\s+(NLAI|NMAXSO|NPART|NFIXF|NSLA|RNFRT|TCNT|DVSNLT|DVSNT|RDRNS|FNTRT|FraHarLosOrm)" src/crop/cropgrowth.f90` → no matches (locals were dropped).
- `grep -inE "\bnlai\b|\bnmaxso\b" src/core/variables.f90` → finds the new module-level declarations.
- `grep -n "flCropNut.*not.*allowed" src/crop/tillage.f90` → one match (runtime gate untouched).
- `pixi run -e test test-pfunit` → `Ok: 1, Fail: 0`; suite count grows by 6.
- `pixi run -e test check-full` → `5 passed, 0 failed` (byte-identical CSV outputs).

## Related

- ADR 0008 — error collection over fatalerr.
- ADR 0016 — per-rotation crop config cache (the cache from
  which `apply_cropwofost_nutrient` reads).
- ADR 0021 — tillage TOML port (mirror shape).
- ADR 0022 — SSDI TOML port (mirror shape).
- ADR 0024 — `dtutil.f90` compatibility shim — same architectural
  thread: physics receives parsed inputs.
- Future: ADR 0026 (N2: soil-side `[nutrients]` TOML block).
- Future: ADR 0027 (N3: lift runtime gate + regression case).
