# Phase 4f-extend SS-6 — irrigation.f90 audit

**Date:** 2026-05-05
**Sub-spec:** SS-6 (legacy reader retirement umbrella)
**Outcome:** **CLOSED — already retired by Phase 4f strangler swap. No code work required.**

## Why this audit exists

The umbrella spec (`2026-05-04-legacy-reader-retirement-design.md`)
described SS-6 as: *"Port the embedded irrigation block (lines 79–326)
currently read via `rdinit` from each `.crp` file."* That description
predated the Phase 4f strangler swap and assumed the legacy reader was
still in the production call graph. Empirical investigation (this audit)
shows the runtime gap is already closed; the only remaining work is
record-keeping.

## Empirical evidence

Verified at HEAD `88bd433` (2026-05-05):

1. **`call readswap()` is absent from production runtime.**
   `grep -rn "call readswap\b" src/ --include='*.f90'` returns one match —
   a comment in `readswap.f90` itself. Per the umbrella status table:
   *"Phase 4f strangler-swap | DONE | `readswap()` no longer called from
   `swap.f90`"*.

2. **`call irrigation(1)` exists only inside `readswap.f90`.**
   Three call sites (lines 2444, 3087, 4087), all guarded by
   `if (schedule .eq. 1)` where `schedule` is read from the legacy `.crp`
   file. With `readswap()` out of the production call graph, none of
   these sites is reachable from the modern binary.

3. **TOML case workdirs contain no `.crp` files.**
   `ls tests/swap-cases/toml/*/*.crp` → no matches. The modern binary,
   when run from a TOML case workdir (per `tests/swap-cases/run_case.sh`),
   has no `.crp` to open even if `irrigation.f90:79` were reached.

4. **Regression is green TOML-only.**
   `pixi run -e test check-full` → `5 passed, 0 failed`. The five
   non-macropore cases (1, 2, 4, 5, 6) all run successfully against
   their TOML inputs without ever invoking `irrigation.f90`'s `case(1)`
   reader.

## What survives where

| Symbol | Where | Status |
|---|---|---|
| `irrigation.f90` `case(1)` (lines ~70–303) — TTutil block | `src/crop/irrigation.f90` | Dead in production runtime; reachable only via parity-test harness which intentionally exercises legacy reading. |
| `irrigation.f90` `case(2)` (per-day calc, lines ~305+) | `src/crop/irrigation.f90` | Live: invoked from `src/core/swap.f90:246` as `if (flIrrigate) call irrigation(2)`. |
| `irrigation_config_t` (top-level `[irrigation]` schema) | `src/config/irrigation_config.f90` | Live, full field set already present. |
| `read_irrigation_toml.f90` | `src/io/toml/read_irrigation_toml.f90` | Live, populates the schema. |
| Adapter `swirfix` + `fixed_events`/`fixed_events_file` copy | `src/io/toml/config_to_variables.f90:933-981` | Live; case 5 uses fixed_events_file end-to-end. |
| Per-rotation `cirrs / cirrthres / dcrit / isuas / perirrsurp / raithreshold / swcirrthres` adapter copy | `src/io/toml/config_to_variables.f90:982-986` | NOT wired (comment notes "Phase 4g territory; the strangler adapter does not touch them"). However, `flIrrigate` stays at default `.false.` for cases that don't author scheduling, so `irrigation(2)` doesn't fire and these globals don't matter for current cases. |

## What "scheduling" means in TOML mode today

No working TOML case authors `[scheduling]` (`schedule = 1`). The legacy
`.crp` files DO author `SCHEDULE = 1` for every crop, but those files
are not read by the modern binary. So in TOML mode:

- `flIrrigate` defaults to `.false.`
- `irrigation(2)` is gated by `if (flIrrigate)` in `swap.f90:246` — never fires
- Scheduling-related globals (`tcs`, `dcs`, `irgthreshold`, …) are uninitialized but unused
- `swirfix=1` cases use the fixed-events adapter path independently

If a future TOML case adds `schedule = 1`, the per-rotation scheduler
parameters need to be wired to legacy globals. Until then this is
unused machinery.

## Retirement gate

The dead `case(1)` block in `irrigation.f90` retires alongside
`readswap.f90` itself in SS-11 (umbrella closeout). At that point:
- Both `case(1)` block and the three `call irrigation(1)` sites in
  `readswap.f90` go.
- `irrigation.f90` keeps only `case(2)` (per-day calc, still live).
- The `flIrrigate` initialization moves from `irrigation(1)` into the
  TOML adapter (it's currently set to `.true.` only by the TTutil
  reader at `irrigation.f90:86` when SCHEDULE=1).

## Optional follow-ups (not required to close SS-6)

- **Stub-error a TOML `schedule = 1` until per-rotation wiring lands.**
  Mirror the SS-4 / SS-5 pattern: `meteorology_config_validate` style
  validator that rejects `schedule = 1` with `ERR_VALIDATION_CROSS_FIELD`
  pointing at this audit. Prevents silent-broken-runs if anyone authors
  scheduling in TOML before the per-rotation wiring exists.
- **Per-rotation schema slot.** The umbrella spec's original SS-6 design
  (per-rotation `crop.<rotation>.irrigation` sub-section) is still the
  right shape if/when scheduling needs TOML support. Defer until a
  case actually needs it.
