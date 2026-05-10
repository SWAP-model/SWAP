# Known Issue: SS-BND Phase 0 — swbotb=2/4/8 Regression Coverage Gap

**Date:** 2026-05-10
**Status:** known-issue
**Related arc:** SS-BND (boundary-conditions state migration)
**Filed at:** B-0.5 (read-only audit; no code changes)

---

## Context

Phase 0 (Tasks B-0.1–B-0.3) closed 8 config gaps in
`bottom_boundary_config_t`, grouped by bottom-boundary gate:

| Gate fields | swbotb path |
|---|---|
| sinmax / sinamp / sinave | swbotb=2 (sine-wave gwl) |
| cofqha / cofqhb / cofqhc / swcofqhc | swbotb=4 (q(h) exp curve) |
| hplate | swbotb=8 (lysimeter) |

B-0.4 verified the adapter (`config_to_variables.f90`) correctly writes
the 8 legacy globals at runtime. However, none of these three paths are
exercised by byte-identical regression.

---

## check-full Coverage (Audited B-0.5)

`swbotb` values from `tests/swap-cases/toml/<case>/swap.toml`:

| Case | swbotb | Path |
|---|---|---|
| 2.grassgrowth | 1 | free drainage (Darcy) |
| 3.macroporeflow | 3 | prescribed head (gwl table) |
| 4.oxygenstress | 3 | prescribed head (gwl table) |
| 5.salinitystress | 3 | prescribed head (gwl table) |
| 6.surfacewater | 3 | prescribed head (gwl table) |
| 1.hupselbrook | 6 | seepage face |

**swbotb=2, swbotb=4, swbotb=8: zero check-full coverage.**

---

## Risk Assessment

**Current regression risk: low.** The 8 fields are covered by pFUnit
unit tests (B-0.4). They populate legacy globals only; `BoundBottom`
logic is unchanged from pre-arc behavior.

**User-facing risk: moderate.** A future reader or adapter regression
on swbotb=2/4/8 paths will not be caught by the byte-identical suite.
Any user running a TOML simulation with these bottom-boundary modes
relies on parsing correctness that is not regression-guarded.

---

## Recommended Fixture Authoring

Out of scope for this arc. Deferred per plan B-0.5 (analogous to the
2026-05-08 nutrient-regression-case deferral). Suggested future cases:

- `tests/swap-cases/toml/7.sinegwl/` — swbotb=2
- `tests/swap-cases/toml/8.qhcurve/` — swbotb=4
- `tests/swap-cases/toml/9.lysimeter/` — swbotb=8

Each requires reference output from a known-good binary before inclusion
in the check-full harness.

---

## Cross-References

| Artifact | Path |
|---|---|
| Discovery | `docs/superpowers/specs/2026-05-10-state-migration-boundary-discovery.md` |
| Design | `docs/superpowers/specs/2026-05-10-state-migration-boundary-design.md` |
| Plan | `docs/superpowers/plans/2026-05-10-boundary-state-migration.md` |
| Successor ADR | `docs/adr/0035-state-migration-boundary.md` (drafted at B-2.7) |
| pFUnit config tests | `tests/unit/config/test_boundary_phase0_fields.pf` |
