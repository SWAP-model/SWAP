---
title: "ADR 0016 — Per-rotation crop config cache and config-passing direction"
date: 2026-05-02
status: accepted
---

# ADR 0016: Per-rotation crop config cache and config-passing direction

## Context

The `.crp` port (Phases 1-4 of Phase 4f) replaces three legacy crop-data
readers (`readcropfixed`, `readwofost`, `readgrass`) with the typed-config
pipeline (ADR 0007). These readers differ from the swap.dra port in two
important ways:

1. **They are called per-rotation, not once at simulation init.** A SWAP case
   has N crop rotations (N = 1 for simple cases, N = 5 for grassgrowth, N = 3
   for surfacewater, etc.). Each rotation entry references a `.crp` file via
   `cropfil(icrop)`; the legacy reader is invoked at the start of each
   rotation by `cropgrowth.f90`'s task=1 dispatch. The file may be a different
   physical file per entry, or the same file referenced multiple times.

2. **The three crop modes (cropfixed / cropwofost / cropgrass) are
   independent.** A case may mix types in one rotation (case 1 hupselbrook
   uses all three). Each mode has its own legacy reader, its own typed config
   type (`cropfixed_config_t`, `cropwofost_config_t`, `cropgrass_config_t`),
   and its own runtime initialization.

The yearly meteorology CSV pattern (ADR 0013) faced a structurally similar
problem: many years of meteo data, accessed at runtime by the legacy
`MeteoCSVDetYear` reader. The chosen solution was to load all years into an
in-memory cache at config-load time and dispatch into the cache from the
runtime per timestep.

The same pattern fits the per-rotation crop case: load all rotation `.crp.toml`
files at config-load time, store the parsed content in parallel arrays on
`crop_config_t`, and dispatch into the cache from `cropgrowth.f90`'s
per-rotation init.

A separate but related question is how the runtime accesses the cache. The
existing pattern across the codebase is that legacy computation subroutines
read inputs via `use variables, only: <globals>`. The typed-config adapter
(`config_to_variables`) writes the same globals so the runtime is unaware of
the typed pipeline's existence. This decoupling has served the strangler
migration well, but it constrains future work: testability is poor (mutating
shared state per test), and parallelism is impossible (globals are not
thread-safe). The natural endpoint is for runtime subs to take config + state
as explicit arguments.

## Decision

### Part A — Per-rotation crop config cache

Add a parallel-array cache to `crop_config_t`, one allocatable array per
crop mode:

```fortran
type :: crop_config_t
   ! ... existing rotation metadata ...
   type(cropfixed_config_t),  allocatable :: rotation_cropfixed(:)
   type(cropwofost_config_t), allocatable :: rotation_cropwofost(:)
   type(cropgrass_config_t),  allocatable :: rotation_cropgrass(:)
end type
```

Each array is sized to the number of rotations. Per rotation entry `i`, the
loader populates exactly one slot based on `rotation_type(i)`:

- `rotation_type(i) = 1` → `rotation_cropfixed(i)` populated; the other two
  arrays' slot `i` remain at their default-initialized (unpopulated) state.
- `rotation_type(i) = 2` → `rotation_cropwofost(i)` populated.
- `rotation_type(i) = 3` → `rotation_cropgrass(i)` populated.

Each crop config type carries a `populated :: logical = .false.` sentinel
set true at the end of its in-place TOML reader. The runtime dispatch checks
the sentinel.

### Part B — Module-level pointer for runtime access (transitional)

A new module `src/crop/crop_config_global_mod` provides a public pointer:

```fortran
type(crop_config_t), pointer :: crop_config_global => null()
```

`config_to_variables`, at the end of its crop block, sets
`crop_config_global => config%crop`. `cropgrowth.f90` (and any other runtime
sub that needs per-rotation crop content) reads through this pointer.

The pointer is **transitional** — see Part C.

### Part C — Future direction: pass config + state explicitly

Once the per-rotation cache is in place across all three crop modes
(end of Phase 4), the next architectural step is to **eliminate the
module-global mutation pattern entirely**. Computation subroutines —
starting with the crop-related ones, eventually generalizing — should
take typed config and typed state as explicit arguments:

```fortran
! Today (post-Phase 4):
subroutine cropfixed(task)
   use variables, only: idev, kdif, gctb, ...
   ! ... reads/mutates module globals ...
end subroutine

! Future direction:
subroutine cropfixed(task, cfg, state)
   class(cropfixed_config_t), intent(in)    :: cfg
   class(crop_state_t),       intent(inout) :: state
   ! ... pure on (cfg, state) ...
end subroutine
```

This eliminates `crop_config_global`, the `populated` sentinels, and the
need for a module-level singleton. It also unlocks:
- Per-rotation isolation (no shared state between rotations).
- Parallel evaluation of independent cases.
- Simple unit testing (compose any config + state, call the sub, check the
  output state).
- Python bindings that pass config dictionaries through to the runtime
  without any global setup.

This direction is **out of scope for the `.crp` port phases**; it is the
natural cleanup once all crop readers are ported and the `crop_config_global`
pointer is the only thing the runtime needs. A follow-on spec sequences this
work after Phase 4.

## Consequences

**Positive:**

- The per-rotation cache makes the typed pipeline complete for crop input
  data: every field a rotation needs is loaded at config-load time, validated
  before the simulation starts, and indexed by rotation at runtime.
- The pattern is identical for all three crop modes, so Phases 2 and 3 are
  pure repetition of Phase 1's plumbing with different schemas.
- The cache positions the codebase for the eventual config-passing future
  without committing to it now. Phase 4 ends with the legacy fallback
  removed and the cache as the sole source of crop runtime data.
- The future direction (Part C) is recorded explicitly so that, after
  Phase 4, the next architectural step is documented rather than rediscovered.

**Negative:**

- Memory overhead: the cache holds parsed crop content for all rotations
  even when the simulation is currently in rotation 1 of N. For typical
  cases this is negligible (a few KB per rotation), but very-large-N cases
  (many decades of grass rotations) might see noticeable memory growth.
  Mitigation: the parsed content is on the order of single-digit KB per
  rotation, so even 1000 rotations is ~10 MB. Acceptable.
- The `crop_config_global` pointer introduces a singleton-like coupling
  during the transition. Marked transitional in Part C; teardown trigger is
  the config-passing follow-on spec.
- The `populated` sentinel adds a field to each crop config type that
  exists only to drive runtime dispatch. Marked transitional; same teardown
  trigger as the global pointer.

**Neutral:**

- The cache layout (parallel arrays per type, one populated slot per
  rotation) is a Fortran-idiomatic choice. A more polymorphic alternative
  (one allocatable array of `class(crop_config_base_t)`) was considered and
  rejected: Fortran's `select type` ergonomics are clunky compared to
  parallel arrays, and mixing types in one array makes per-type allocation
  fiddly. Parallel arrays are cheap when most rotations are one type
  (case 6: all type-1; case 2: all type-3). For mixed-type cases (case 1),
  exactly one slot is populated per rotation across the three arrays.

## Teardown plan

Strangler-fig discipline (per project convention): every transitional
element introduced by this ADR has a documented removal trigger.

| Element | Trigger | Removal | Replacement |
|---|---|---|---|
| Per-rotation legacy reader fallback (`if/else` in `cropgrowth.f90`) | End of Phase 4 (`.crp` port) — all three crop types ported | The `if (associated(crop_config_global) .and. ... .and. populated)` guards around `readcropfixed`/`readwofost`/`readgrass` calls | Unconditional `*_init_from_config` calls; the dispatch is on `rotation_type(icrop)` only |
| Legacy `readcropfixed`/`readwofost`/`readgrass` subs in `readswap.f90` | Parity tests rewritten against fixture values | The three subroutines (~2200 lines combined) | Fixture-based parity tests storing expected `variables`-module state |
| `crop_config_global` module pointer | Config-passing follow-on spec lands (Part C) | The `crop_config_global_mod` module entirely | `cropfixed_init_from_config(cfg, ...)` and downstream subs take `cfg`/`state` as explicit arguments, threaded from `core/swap.f90` |
| `populated :: logical` sentinel on crop config types | Same as `crop_config_global` | The `populated` field on each crop config type | Dispatch happens at the call site based on `rotation_type(icrop)`; no runtime sentinel needed because there is no fallback |
| `parallel arrays` (`rotation_cropfixed`, `_cropwofost`, `_cropgrass`) | Not transitional. The cache shape itself is the long-term design. The arrays remain even after the config-passing direction lands; only the access pattern changes. | — | — |

## Pattern reference

This ADR generalizes the per-year cache pattern documented in ADR 0013 (CSV
meteorology) to per-rotation crop data. The implementations:

- `src/config/crop_config.f90` — adds the three parallel arrays.
- `src/io/toml/read_crop_toml.f90` — populates the right slot per rotation.
- `src/crop/crop_config_global_mod.f90` (new) — module-level pointer.
- `src/crop/cropfixed_init.f90` (Phase 1, new) — runtime sub reading from
  the cache and writing to legacy globals.
- `src/crop/cropwofost_init.f90` (Phase 2, future).
- `src/crop/cropgrass_init.f90` (Phase 3, future).
- `src/crop/cropgrowth.f90` — dispatch site at the per-rotation init point.
