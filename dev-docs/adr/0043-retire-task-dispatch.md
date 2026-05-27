---
title: "ADR 0043 — Retire Integer Task-Dispatch Lifecycle Multiplexing"
date: 2026-05-27
status: accepted
branch: development
---

# ADR 0043: Retire Integer Task-Dispatch Lifecycle Multiplexing

**Status:** accepted
**Date:** 2026-05-27
**Branch:** `development`

## Context

Many SWAP 4.2.0 compute subroutines carry the legacy idiom of an integer
lifecycle selector: `subroutine X(task, …)` (also spelled `iTask`, `swtask`)
with an internal `select case (task)` whose branches are *lifecycle phases*,
not algorithm options:

- `case (1)` — initialization
- `case (2)` — per-step / rate computation
- `case (3)` — terminal / output

The orchestrator then drives them with magic-int calls
(`call SoilWater(1, state)` … `call SoilWater(2, state)` …), which the reader
must decode against the callee's case labels. This is the opposite of the
canonical pattern already in place — `swap_init` / `swap_run_step` /
`swap_close` (free module procedures, time loop in `swap_main`), with each
subsystem exposing a type-bound `state%X%init(...)`, and per-step work as free
procedures called with no integer (`call Drainage(state)`).

Reading the dispatcher bodies revealed that "initialization" is genuinely
**two phases**, already latent in the code but conflated by the `case (1)`
label:

1. **Construct** — allocate arrays, zero-fill, snapshot config scalars/tables
   into state. Order-free, idempotent, takes no sibling subsystem. This is
   exactly what the type-bound `state%X%init(cfg, …)` procedures already do
   (e.g. `heat_state_init`, `solute_state_init`).
2. **Seed** — compute derived initial state from *sibling subsystems that must
   already be seeded*. E.g. `Solute(1)` computes `cmsy`/`samini` from
   `soil%theta`, `soil%bdens`, `mesh%z`; `Temperature(1)` sets the tsoil
   *profile* and `fquartz/fclay/forg` from `soil%thetas`/`soil%orgmat` and is
   deliberately ordered *after* `SoilWater(1)`. This is order-dependent
   (soilwater → heat → solute) and lives in the old `case (1)` body.

These cannot be merged: folding the seed into the type-bound `init` would
either drop the derived computation (`heat_state_init` takes no `soil`) or
break the ordering. They are distinct phases.

## Decision

Retire the integer task-dispatch idiom. For each such subroutine `X(task, …)`:

- **Phase 1 (construct) stays.** The type-bound `state%X%init(...)` is the
  canonical, order-free constructor; its call site in `swap_init_body` is left
  exactly where it is. We commit to it as *the* init mechanism.
- **`case (1)` → free `x_seed(state)`** (Phase 2): the old body, integer
  dropped, called at the *same site* the dispatcher's `(1)` was called — so
  cross-subsystem ordering is preserved byte-for-byte. It takes the whole
  `state` because it reads sibling subsystems.
- **Step / sub-step cases → named free procedures** taking `state`
  (`x_step`, `x_update`, `x_finalize`, `x_lateral`/`x_balance`, …). The integer
  is dropped; the phase becomes the name.
- **Operation pairs that merely *resembled* lifecycle** (save/restore, lazy
  table build/lookup) get honest names, not lifecycle verbs:
  `soilwater_save_state` / `soilwater_restore_state`;
  `matricflux_build_table` / `matric_flux`.

### Multi-instance invariant

In-process parallel multi-instance execution (multiple `swap_state_t` running
concurrently in one process) is on the near-term roadmap. To support it, this
ADR establishes a binding invariant for all extracted procedures:

> **Every procedure takes `state` explicitly and touches nothing outside it —
> no module-level variables, no `SAVE` locals.**

This makes the procedures re-entrant: N instances each mutate only their own
`state`. Replacing `X(int, state)` dispatch with `x_phase(state)` procedures
is itself a step *toward* this invariant. The Tier 1 + Tier 2 subsystems
covered by the first retirement arc already satisfy it (verified: no module
vars, no `SAVE` in `soilhydraulics`, `solute`, `temperature`, `tillage`,
`irrigation`, `surfacewater`).

### Explicitly out of scope (NOT the anti-pattern)

Genuine option-switches, where the integer encodes an algorithm/option choice
rather than a lifecycle phase, are correct and stay:

- `select case (model)` in `soil/WC_K_models_04_11.f90` (water-content /
  conductivity model choice).
- `reduceva_apply(task_flag)` in `atmosphere/et.f90` — `1=daily / 2=sub-daily`
  computation mode (mild redundancy with its `daily` arg; not lifecycle).
- `SharedSimulation(task)` — file open/read/write/close; dormant.

## Consequences

- `swap_init_body` and `swap_run_step` become self-documenting: a sequence of
  named `proc(state)` calls instead of magic ints.
- Output is byte-identical against the `swap420gf` reference at every step —
  the retirement is a pure rename + signature change with no relocation of
  initialization, so behavior is preserved by construction.
- The arc is scoped to **Tier 1 + Tier 2** (10 independent dispatchers). The
  **crop cluster** (`CropGrowth`, `CropFixed`/`Wofost`/`Grass`,
  `ArableLandGerm`, `SoilManagement`, …) is deferred: its dispatchers call
  each other with magic ints and carry phase-sequence ordering needing its own
  design.

### Forward note — the real parallelism blocker

The literal obstacle to parallel in-process runs is **not** the dispatch
idiom but module-global state in the deferred crop cluster:

- `crop_config_global` — a module-global `pointer`
  (`src/crop/crop_config_global.f90:23`), reassigned on every
  `crop_state_init` (`src/state/crop_state.f90:85`, flagged transitional under
  ADR 0016). Initialising instance A then B leaves the pointer aimed at B;
  A's crop routines then read B's config.
- `SAVE` locals in `src/crop/cropgrowth.f90:601,678` and
  `src/crop/cropgrass_runtime.f90:94` — persistent state shared across all
  instances calling those routines.

Eliminating these (route reads through `state%cfg%crop` / move persistence
into `state`) is the prerequisite for multi-instance and is tracked as a
distinct effort, overlapping the crop-cluster Tier 3 retirement. A possible
end state is a single top-level `state%init(config)` constructor wrapping
Phase 1 (all subsystems) → Phase 2 (seeds, in order) → first-day compute, so
spinning up N instances is a loop that cannot mis-sequence.
