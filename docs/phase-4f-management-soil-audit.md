# Phase 4f-extend SS-7 — management_soil.f90 audit

**Date:** 2026-05-05
**Sub-spec:** SS-7 (legacy reader retirement umbrella)
**Outcome:** **CLOSED — unreachable in TOML runtime; gate (`flCropNut`) defaults to `.false.` and is not set by any TOML adapter or working case. No code work required.**

## Empirical evidence

Verified at HEAD `9ed4a00` (2026-05-05):

1. **All 8 `call SoilManagement(*)` sites in `src/core/swap.f90`** (lines
   191, 317, 324, 327, 334, 371, 408) are guarded by `if (flCropNut)`.

2. **`flCropNut` writes:** the only writes anywhere in `src/` are:
   - `src/io/readswap.f90:474` sets it to `.false.` (legacy reader,
     not in production call graph)
   - `src/crop/tillage.f90:79` — `if (flCropNut) call fatalerr_collected
     ('DoTillage', 'flCropNut = 1 not (yet) allowed')` — an explicit
     stub-error: even the legacy code path rejects `flCropNut = .true.`
     mid-tillage.
   No TOML adapter writes `flCropNut`.

3. **Default value:** Fortran logical defaults to `.false.`, so under
   the TOML pipeline `flCropNut` is `.false.` for the entire
   simulation.

4. **No working case authors a nutrient block** in either TOML or the
   legacy `.swp.template` files.

Conclusion: `management_soil.f90`'s three TTutil `rdinit` calls (lines
129, 142, 197 — `.smm`, `.sme`, `.snp` companion files) are unreachable
from the modern binary's call graph.

## What survives where

| Symbol | Where | Status |
|---|---|---|
| `SoilManagement(task=1)` init block (lines ~31-220 with three rdinit calls) | `src/crop/management_soil.f90` | Dead in production runtime; reachable only if `flCropNut` is somehow set true (which is also explicitly rejected at `tillage.f90:79`). |
| `SoilManagement(task=2..7)` per-day/per-event branches | `src/crop/management_soil.f90` | Same gate — dead in TOML runtime. |
| `flCropNut` global | `src/core/variables.f90:564` | Live but always `.false.` under the TOML pipeline. |
| `tillage.f90:79` `flCropNut=1 not allowed` rejection | `src/crop/tillage.f90` | Active belt-and-suspenders rejection at the tillage code path. |

## Retirement gate

The dead `case(1)` block in `management_soil.f90` retires alongside
`readswap.f90` itself in SS-11 (umbrella closeout). Until then it is
no-op in TOML runtime.

## Optional follow-ups (not required to close SS-7)

- **Stub-error a TOML attempt to enable nutrients.** No TOML schema
  field even exists for `flCropNut` today; any future port should
  start with the SS-4 / SS-5 pattern (validator that rejects until
  per-rotation `[crop.<rotation>.nutrients]` plus `.smm` / `.sme` /
  `.snp` schema slots are wired).
- **Promote the `tillage.f90` ad-hoc rejection** into a typed
  validator on the eventual nutrient config slot.
