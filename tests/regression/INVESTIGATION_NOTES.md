# Regression fixtures — investigation notes

**Status:** open, bumped to a follow-on spec
**Opened:** 2026-04-22 (Rescue Phase 1, Task 1)

## What changed

The rescue committed to **gfortran-only** during Phases 1–4 (see `docs/adr/0001-gfortran-first.md` — written in Task 8 of the Phase 1 plan). Phase 0's baseline regression numbers were recorded under **ifx** (Intel), because `pixi.toml` silently hardcoded `FC=ifx` in the production configure tasks and `meson.build` had an Intel-specific flag path including `-init=zero`.

When the compiler was swapped to gfortran (plus `-finit-local-zero` to match Intel's zero-initialization semantics), five of six cases matched the existing ifx-based fixtures cleanly. Two cases diverged:

- **`macropore`**: `DRAINAGE` drifts by ~21 at year 1998 (out of ~100), ~3 at year 1999, with corresponding small `GWL` deltas (0.3–2.0). This is large enough to be a real behavioral difference, not floating-point noise.
- **`oxygenstress`**: the pre-existing `MOWDM` deviation (max 85 in year 1995) was already visible under ifx; same magnitude and same year under gfortran — so this deviation is NOT compiler-driven. It is a pre-existing physics or fixture issue that has been tolerated from well before the rescue.

## What we did

Rather than picking a side, we kept both reference sets:

- `*_expected.json` — historical ifx-produced values. **Unchanged.** These document what the upstream-Intel-compiled reference produced.
- `*_expected_gfortran.json` — regenerated from the gfortran+finit-local-zero build on 2026-04-22. These are what the regression harness actually compares against under the current compiler policy.

The harness `CASES` dict points at the `_gfortran` files only. A future compiler-policy change would add a selector — not in scope for the rescue.

## What needs investigation (not part of the rescue)

The open question is whether the **macropore DRAINAGE divergence** is a physics problem or a compiler-flag artefact. Specifically:

1. Does gfortran `-finit-local-zero` cover every initialization path that ifx's `-init=zero` covered? (`SAVE`d module variables, allocated arrays, derived-type components, and COMMON blocks are all separately-controllable dimensions.)
2. Do additional FP-model flags (`-ffp-contract=off`, `-fno-unsafe-math-optimizations`, `-fno-fast-math`) move gfortran's output closer to ifx's? If so, the divergence is numerical compiler choice. If not, something in the macropore code is genuinely compiler-sensitive (stale pointer, undefined-order-of-evaluation arithmetic, etc.).
3. For oxygenstress MOWDM: is the 1995-spike a crop-growth/mowing-schedule physics bug that the fixture has always masked, or a numerical artefact of a scheduled event? Since gfortran and ifx both hit the same number, this is almost certainly physics, not compiler.

None of these block the rescue. They are tracked here as open items to revisit during Phase 4 module cleanup or as a dedicated physics-audit follow-on spec.

## Reproducing the ifx reference (if needed for an investigation)

The legacy SWAP 4.2.0 Intel-compiled Linux binary is preserved at `tests/reference/swap420`. Running it via `pixi run swap-ref` generates output comparable to the historical `*_expected.json` fixtures (modulo the MOWDM deviation the fixtures themselves encode).
