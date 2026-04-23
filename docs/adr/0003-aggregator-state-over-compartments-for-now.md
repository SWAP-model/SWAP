# ADR 0003 — Aggregator state now, compartment state later

Status: accepted (2026-04-23, during rescue Phase 2)

## Context

The legacy SWAP codebase used a global `variables` module that every physics routine read and wrote. This blocks parallelization, breaks test isolation, and makes the dependency graph opaque.

The rescue spec's original forward plan (captured in the now-deleted `docs/next-refactoring-task.md`) was to skip directly to compartment-based state (ponding, canopy, snow, soil, saturated zone, etc.) where process functions explicitly pass water/heat/solute between compartments. This is the long-term target.

At the rescue baseline `e256bc0` the codebase had progressed partway through an aggregator pattern (one `swap_state_t` composed of domain-specific `*_state_t` types) without completing it. Multiple modules still `use variables` and rely on legacy globals.

## Decision

The rescue stays on the aggregator pattern through Phase 4 exit. Compartment-based state is a follow-on spec.

Rationale:

1. Aggregator is partially done; finishing it is bounded work.
2. Compartment state is a larger physics-aware redesign; attempting it mid-rescue risks physics regressions we cannot easily detect.
3. The aggregator already enables test isolation and state passing; it is sufficient for Phase 4's fix-and-clean work.
4. Compartment state will require TDD-first rewrite per compartment with regression checks; that belongs to its own dedicated spec.

## Consequences

- **Positive**: the rescue has a bounded, achievable goal. Phase 4 closes the pre-compartment gates (every module has clean `config_t` / `initial_t` / `state_t` separation) rather than attempting a physics-level redesign.
- **Positive**: Phase 3 test coverage builds on the aggregator; state-type lifecycle tests ported during Phase 4 remain useful when compartment state lands later.
- **Negative**: aggregator is not the final state. Phase 4 produces code that gets rewritten again when compartment state lands. Some work is "throwaway".
- **Negative**: `src/core/variables.f90` and `src/core/swap_state_sync.f90` stay in the tree until the aggregator is complete. They are flagged as legacy in `docs/code-style.md`.

## When to revisit

After `rescue/complete`. A dedicated compartment-state spec and plan supersede this ADR.
