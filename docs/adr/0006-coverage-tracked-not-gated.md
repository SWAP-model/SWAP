---
title: "ADR 0006 — Coverage is tracked, not gated"
---

# ADR 0006 — Coverage is tracked, not gated

Status: accepted (2026-04-24, during rescue Phase 3)

## Context

During Phase 3 of the rescue spec we produce the first coverage baseline for
the SWAP modernization tree. Coverage tools (`gcov`, `gcovr`) are available
and cheap to run via `pixi run -e coverage coverage-report`. The question is
whether any coverage number becomes a gate on Phase 4 work — for example,
"PR cannot merge if line coverage drops below 50%".

## Decision

Coverage is **tracked**, not **gated**. The baseline in
`docs/coverage-baseline.md` is a reference point for Phase 4 refactors. No
CI check blocks a change on coverage; no pixi task fails on a coverage target.

## Consequences

Positive:

- Phase 4 can move fast on risky refactors (crop consolidation, I/O
  consolidation, `fatalerr` replacement) without chasing coverage percentages
  during the physics-preserving window.
- Authors can add characterization tests where they matter (pure routines,
  state lifecycles, TOML readers) without artificially padding counts on
  trivial getters.
- The baseline document stays short and actionable.

Negative:

- A future drop in coverage is only caught in review, not mechanically.
  Mitigation: the baseline table is per-domain, so a material regression in
  one domain is visible at a glance in the next coverage re-run.

## Revisit trigger

When rescue exits at `rescue/complete` (end of Phase 4), re-evaluate whether
coverage should become a gate in the compartment-state follow-on spec. At
that point the code is stabilised enough that a "never goes down" ratchet
may be cheap to add.
