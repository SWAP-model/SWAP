# ADR 0053 — Compartment exchange records

**Status:** Accepted (2026-07-04) — architecture; per-surface implementation
arcs follow (see design doc
[`2026-07-03-arc-t2a-exchange-records-design.md`](../2026-07-03-arc-t2a-exchange-records-design.md)).
**Implements:** ADR 0051 (core-vs-component boundary) — the mechanism that makes a
compartment's coupling *explicit* so its implementation can be swapped (near-term
target: external WOFOST). **Relates to:** ADR 0043 (state threaded, no hidden
state), ADR 0050 (BMI/XMI substrate; the bottom-boundary XMI arrays are already a
de-facto exchange record).

## Context

The state migration removed the bare globals but left the coupling: any compute
routine can still read/write any sibling compartment's `state%X` fields
(§1 of the 2026-07-03 review — e.g. crop reads `soil%tra`, soil-hydraulics reads
`heat%tsoil`). This implicit all-touch-all coupling blocks (a) swapping a
compartment's implementation (you'd have to chase every reach-across site) and
(b) per-component testing/threading. A code-verified inventory found the runtime
crossings concentrate on ~6 surfaces, soil-water being the hub.

## Decision

Introduce **typed, directional exchange records — one per coupling surface** — as
the single channel between two compartments. A compartment reads/writes only its
side of the relevant record and never touches a sibling's state directly.

- **Ownership:** a `swap_exchange_t` aggregate held as `state%exchange` on
  `swap_state_t`, threaded like every other state record (incremental; may evolve
  to first-class MF6-style objects later).
- **Copy semantics:** a record holds the crossing *values* (copied at the
  boundary), not pointers into sibling state — so the surface is a real contract
  and an *external* component (no shared memory) is a drop-in producer/consumer.
- **Directional + cadence-aware:** `from_<producer>` / `to_<consumer>` groups with
  an explicit update order (e.g. crop fields write-once-per-day; actual
  transpiration finalizes at the day boundary).
- **One combined crop water-relations record** (crop ↔ atmosphere+soil-water) —
  the full crop↔SWAP interface (~18 dynamic fields: ET/canopy inputs + root sink +
  stress feedback), because relative transpiration inherently spans an atmosphere
  quantity (`ptra`) and a soil quantity (`tra`), and an external crop must feed
  SWAP's ET, not just the sink. Crop config-*constants* stay config-passed.
- **Scope now:** the four *clean* surfaces — `crop_water`, `heat_soil`,
  `drain_soil`, `atmos_soil`. **Deferred:** `solute` (a many-to-one consumer →
  a driver-bundle, not a bilateral pair) and `surfacewater` (bidirectional
  pond/runoff loop). **Excluded (not coupling):** the `atmo%cumu`/`intr` water-
  balance accumulation, which is output accounting.

## Consequences

- Each coupling surface becomes a readable, minimal, typed contract; a
  compartment's entire dependency on the rest of SWAP is its record(s).
- **Enables the external-WOFOST plug:** once `crop_water` is the only crop↔rest
  channel, `crop_step()` is a swappable producer — native crop, in-house WOFOST,
  or external WOFOST over BMI/XMI (structurally identical to the running
  SWAP↔MODFLOW 6 coupling; daily cadence; the record is the contract).
- Improves testability (a compartment tests against a small record, not a whole
  `swap_state_t`) and is a prerequisite for safe threaded execution.
- **Byte-identical rollout:** each surface is cut over by populating the record
  with a pure copy of the value its reach-across line currently reads, then
  repointing that line — output unchanged; `check-fast` + `check-bindings` green
  per surface.

## Alternatives considered

- **One `swap_exchange_t` with all crossings in a single record.** Rejected —
  recreates all-depend-on-all. One record per surface keeps each contract minimal.
- **First-class exchange objects off `swap_state_t` (full MF6 split) now.**
  Deferred — more signature churn than the incremental `state%exchange` form,
  which can evolve into it.
- **Keep reaching into sibling state (status quo).** Rejected — it is the
  residual coupling that blocks component swap and threading.

## Sequencing

`crop_water` (pattern-setter + WOFOST-enabler) → `heat_soil` → `drain_soil` →
`atmos_soil`; then a separate pass for the deferred solute/surfacewater surfaces.
Details, field tables, and the update protocol: the T2-A design doc.
</content>
