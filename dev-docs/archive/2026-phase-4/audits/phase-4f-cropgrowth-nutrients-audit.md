# Phase 4f-extend SS-8 — cropgrowth.f90 nutrient block audit

**Date:** 2026-05-05
**Sub-spec:** SS-8 (legacy reader retirement umbrella)
**Outcome:** **CLOSED — unreachable in TOML runtime; gate (`flCropNut`) defaults to `.false.` and is not set by any TOML adapter or working case. No code work required.**

## Empirical evidence

Verified at HEAD `9ed4a00` (2026-05-05):

1. **The `rdinit` at `src/crop/cropgrowth.f90:1155`** is inside an
   `if (flCropNut) then` block (line 1151). Two other `flCropNut` gates
   exist at lines 351 and 1281 in the same file.

2. **`flCropNut` writes:** identical situation to SS-7 — only
   `readswap.f90:474` (legacy, sets `.false.`) and `tillage.f90:79`
   (rejects `.true.`). No TOML adapter writes the flag.

3. **Default value:** `.false.`. The `cropgrowth.f90:1151` gate is
   never true under the TOML pipeline.

4. **No working case authors nutrient parameters** in TOML or legacy.

Conclusion: the `LRNR/LSNR/NLAI/NLUE/...` nutrient block read at
`cropgrowth.f90:1155+` from the per-rotation `.crp` file is
unreachable from the modern binary.

## What survives where

| Symbol | Where | Status |
|---|---|---|
| `cropgrowth.f90:1151-1280` nutrient `rdinit` block + N-P-K state machine | `src/crop/cropgrowth.f90` | Dead in production runtime; gated on `flCropNut`. |
| `cropgrowth.f90:351` and `:1281` `flCropNut` checks | `src/crop/cropgrowth.f90` | Same gate; dead in TOML mode. |

## Shared gate with SS-7

SS-7 (`management_soil.f90`) and SS-8 (`cropgrowth.f90` nutrients)
share the `flCropNut` runtime gate. If a future case requires nutrient
simulation, **both must be ported simultaneously** (and `tillage.f90`'s
`fatalerr_collected` rejection lifted) — they are one logical feature
spread across three files.

## Retirement gate

The dead nutrient block in `cropgrowth.f90` retires alongside
`management_soil.f90` and `readswap.f90` in SS-11 (umbrella closeout).
Per the shared-gate note above, that retirement is one logical change.

## Optional follow-ups (not required to close SS-8)

See SS-7 audit (`phase-4f-management-soil-audit.md`) — the
recommendations are identical and apply to the combined nutrient
feature, not to either file in isolation.
