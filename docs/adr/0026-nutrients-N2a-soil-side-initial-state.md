---
title: "ADR 0026 — [nutrients] N2a: soil-side initial state + SorpCoef"
date: 2026-05-07
status: accepted
---

# ADR 0026: [nutrients] N2a — soil-side initial state + SorpCoef

## Context

The `[nutrients]` umbrella reactivates the SWAP nutrient
subsystem on the modern TOML pipeline. ADR 0025 (N1) wired the
crop-side parameter set per-rotation. This ADR (N2a) ports the
soil-side initial state.

The legacy code read soil-nutrient initial concentrations from
`<project>.snp` (12 reals: `FOM_t(1..8)`, `Bio_t`, `Hum_t`,
`cNH4_t`, `cNO3_t`) inside `SoilManagement(1)`, which was
collapsed to `return` in SS-C step 3 of the legacy-readers
deletion arc. The `<project>_nut.end` diagnostic dump in
`SoilManagement(7)` then read `<project>.snp` line-by-line as a
template, writing the current pool values to `_nut.end`.

Pre-flight inventory for N2 found three legacy files:
- `<project>.snp` — initial pool concentrations (this ADR's
  scope).
- `<project>.smm` — material-property overrides (deferred
  indefinitely; hardcoded defaults in `Wofost_SoilParameters`).
- `<project>.sme` — timed soil management events / fertilizer
  applications (N2b's scope).

Plus three smaller findings:
- `SorpCoef` is genuinely uninitialised in modern `src/` — used
  by `wofost_soil_watern`, `wofost_soil_balancecheck`,
  `wofost_soil_amendments`, `wofost_soil_cropresidues` but never
  assigned anywhere. Real bug.
- `DryBD` is auto-derived from `BDENS(1)` (the soil's bulk
  density) in `Wofost_SoilParameters`. No TOML port needed.
- `<project>_nut.end` writer would fail at runtime without a
  `.snp` template file present.

## Decision

This ADR (N2a):

1. **Add a top-level `[nutrients]` TOML block.** Optional;
   `present` flag distinguishes "user supplied" from "absent —
   use defaults". Carries `sorp_coef` (top-level scalar) and an
   `[nutrients.initial]` sub-table with the 12 pool values.
2. **Default `SorpCoef = 0.0`.** This is a behavioural choice
   — legacy code left `SorpCoef` uninitialised, so the
   deterministic 0.0 default may differ numerically from any
   legacy non-deterministic baseline. 0.0 = "no sorption", which
   is the physically conservative choice and matches what most
   compilers' module-init zero would produce.
3. **`apply_nutrients` is called unconditionally** from
   `config_to_variables`, regardless of `flCropNut`. The adapter
   sets `SorpCoef` and pool values whether or not nutrients are
   enabled at runtime. When `cfg%present = .false.`, the
   defaults (zero) flow through — same as the pre-arc
   uninitialised behaviour for any case that didn't read these
   globals.
4. **Retire `SoilManagement(7)`'s `_nut.end` dump.** Body
   collapses to `return`. The `flCropExt` write-back at the top
   of the legacy case is also dropped (it depended on the same
   `cropext` mechanism). A future arc can add a clean CSV-style
   nutrient-pool dump if anyone asks. Note: the `cropext`
   module global itself was kept — `case (6)` still references
   it.

## What N2a does NOT do

- **N2b — timed soil management events.** Without N2b, simulations
  run with **zero amendments**: natural mineralization only.
  That's a useful baseline for testing.
- **Material-property overrides** (legacy `<project>.smm`). The
  17 hardcoded materials in `Wofost_SoilParameters` (Cattle
  manure → Spruce needles, with default AppAge / OrgMatFrac /
  OrgNFrac / NH4NFrac / NO3NFrac) stay as-is. Deferred
  indefinitely.
- **Lift `tillage.f90:73`.** N3's job, after both N2a + N2b
  land. With no amendments wired (N2b pending), a
  flCropNut=true simulation would run from zero pools with no
  fertilizer inputs — useful but not the full picture.
- **Add a regression case** with `flCropNut=true`. Cannot run
  end-to-end until the runtime gate is lifted (N3).
- **Replace `_nut.end` with a CSV writer.** Future enhancement.

## Schema

```toml
[nutrients]
sorp_coef = 0.005

[nutrients.initial]
fom  = [0.5, 0.3, 0.2, 0.1, 0.5, 0.3, 0.2, 0.1]   # FOM_t(1..8)
bio  = 0.4
hum  = 8.0
cnh4 = 0.001
cno3 = 0.005
```

Both `[nutrients]` and `[nutrients.initial]` are optional;
absent → defaults (zero pools, zero sorption).

## Consequences

- A TOML config with `[nutrients]` populated loads and
  validates without error; pool values flow through to the
  legacy globals.
- `SorpCoef` now has a deterministic value (0.0 by default,
  user-overridable). Eliminates the uninitialised-variable bug.
- `<project>.snp` and `<project>_nut.end` files no longer have
  any role in the modern pipeline. The diagnostic dump path is
  retired.
- All five existing regression cases produce byte-identical
  CSV outputs (no case has flCropNut=true; no case reads the
  pool globals or SorpCoef on its active code path).
- N3 will lift `tillage.f90:73` and add a regression case with
  populated `[nutrients]`; verifying byte-identical against the
  legacy binary in N3 is the end-to-end correctness gate.

## Acceptance

- `pixi run -e test test-pfunit` → `Ok: 1, Fail: 0` (601 tests,
  1 disabled — was 583 pre-N2a; +18 from N1 + N2a).
- `pixi run -e test check-full` → `5 passed, 0 failed`
  (byte-identical CSVs).
- `grep -nE "_nut\.end" src/crop/management_soil.f90` → matches
  only in the retirement-comment block (no executable code).
- `grep -n "snp\b" src/crop/management_soil.f90` → matches only
  in the retirement-comment block.
- `grep -n "flCropNut.*not.*allowed" src/crop/tillage.f90` →
  one match (runtime gate intact).
- ADR 0026 committed.

## Related

- ADR 0008 — error collection over fatalerr.
- ADR 0024 — `dtutil.f90` compatibility shim — same
  architectural thread.
- ADR 0025 — [nutrients] N1 (crop-side adapter).
- Future: ADR 0027 (N2b: amendments).
- Future: ADR 0028 (N3: lift runtime gate + regression case).
