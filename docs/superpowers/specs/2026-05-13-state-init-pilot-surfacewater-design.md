---
title: "State Init Pilot — Surfacewater — Design Spec"
date: 2026-05-13
status: draft
---

# State Init Pilot — Surfacewater

## Context

After the 2026-05-12 cohort-flattening arc, the state records on `swap_state_t` have a clean shape (flat fields + named reset procedures) but the **initialization story** remains tangled. A single subsystem's init logic is currently spread across four places:

1. **Config validate + finalize** (clean): `config/<section>_config.f90` modules with `validate`/`finalize` procedures producing config-self math (e.g. `alphaw` normalization, `abs(wldip)`, `nrlevs` clobber).
2. **Per-subsystem `<X>_init` procedures** — sometimes in `src/state/`, sometimes in subsystem modules (`src/drainage/surfacewater_init.f90`, `solute_init` in `src/solute/`, `heat_init` in `src/heat/`). Inconsistent signatures and scope.
3. **`config_to_variables.f90`** (1789 LoC adapter) — populates legacy variables.f90 globals AND transient `*_init_buf` module-level buffers as a bridge between config-load time and state-allocation time.
4. **Inline seeds in `swap_init` body** — `state%timecontrol%iyear = tc_iyear_init_buf`, `state%soilwater%pondini = pondini_init_buf`, the swinco=3 atmosphere block, etc.

The same `state%X` gets touched in 4-5 hops: config schema → `config_to_variables` → transient buffer → `<X>_init` allocation → manual seed in `swap_init`. The transient buffer indirection exists because state isn't allocated yet when `config_to_variables` runs.

A separate but related smell: init logic is buried inside `<Subsystem>(task=1, state)` Phase-1 dispatchers. `swap_init` calls `SurfaceWater(1, state, ...)`, which internally calls `surfacewater_init(state)`. Init order, dependencies, and conditional activation become invisible at the driver level.

## Decision: target architecture

Adopt the modern scientific-modeling-framework pattern: **top-level orchestration with state-bound init implementation**.

```
1.  load_swap_config(file, config)            ← S1 read
2.  config%validate(errors)                   ← S2 validate
3.  config%finalize(errors)                   ← S3 config-self math
                                                  (alphaw normalization, etc.)

4.  For each subsystem in dependency order:
       call state%X%init(config%X, [dims])    ← S4 state init
                                                  alloc + zero + seed from config
                                                  + readswap-style "make-shape" math
                                                  (sttab, swst init, hatm default, etc.)
       — driver-side flag guards for optional subsystems
       — core components (timecontrol, soilwater, atmosphere) init unconditionally

5.  Physics modules read config (read-only) and mutate state during their Phase-1
    and runtime — unchanged. cofgen, mfluxtable, etc. all populated here.
    The `<Subsystem>(task=1, state)` dispatcher idiom remains for physics setup work
    that genuinely belongs to a Phase-1 step.
```

After S4 every state is **well-formed but physically inert** — arrays allocated, defaults set, readswap-style shape math computed. S5 (physics) writes the heavy stuff.

### What "state init" includes (S4)

- **L0 — allocation:** allocate arrays sized from typed config + grid dimensions.
- **L1 — zero defaults:** explicit zero/false initialization of scalars, logicals, integers (Fortran's per-component init takes care of declarations; init body covers anything the declaration didn't).
- **L2 — config seed:** copy directly from typed config slice (`state%surfacewater%wls = config%surface_water%wlact - config%drain%altcu`). No transient buffers.
- **L3 — readswap-style shape math:** deterministic transformations from config to state shape that put data in the form physics expects. Examples: `sttab` (open-channel storage table) built from drainage geometry; `swst` initial storage computed from `sttab` + `wls`; `hatm = -2.75e5` (legacy soilhydraulics default); `ZDraBas` macropore drainage basis derived from drainage config.

State `init` does **not** include:

- **Physics math** — VG-Mualem fitting (cofgen), JvL matric-flux integration (mfluxtable), Hooghoudt/Ernst drainage equations. These live in physics modules (`SoilHydraulics(1)`, `MatricFlux(1)`, `bocodrb`, etc.) and read config + mutate state during their normal Phase-1 / runtime execution.
- **Validation or config-self math** — those are S2 / S3.

### Implementation pattern

Each state record gets one type-bound `init` procedure on the state module:

```fortran
type :: surfacewater_state_t
   ! ... fields (per the 2026-05-12 flatten-cohorts arc) ...
contains
   procedure :: init                        => surfacewater_state_init
   procedure :: reset_intermediate          => surfacewater_reset_intermediate
   procedure :: reset_cumulative_drainage   => surfacewater_reset_cumulative_drainage
   procedure :: reset_cumulative_reservoir  => surfacewater_reset_cumulative_reservoir
end type surfacewater_state_t
```

The implementation lives in the state module (`src/state/<X>_state.f90`). It imports the typed config slice(s) it needs and any utility modules its math requires (e.g., `surfacewater_utils` for `swstlev`). State modules grow modestly but stay coherent — they own their full lifecycle.

### Driver-side activation gating

The `init` procedure body is **flag-agnostic** — when called, it inits unconditionally. Activity gating lives at the driver:

```fortran
! In swap_init S4 pass:

! Core components — always init, no flag guard:
call state%timecontrol%init(config%simulation)
call state%soilwater%init(config%soil, numnod, numlay)
call state%atmosphere%init(config%soil, config%meteo)

! Optional subsystems — driver-side flag guards (matches existing convention):
if (flTillage)       call state%tillage%init(config%soil, numlay)
if (flTemperature)   call state%heat%init(config%heat, dims)
if (flDrainage)      call state%drainage%init(config%drain, dims)
if (flSurfaceWater)  call state%surfacewater%init(config%surface_water, config%drain, numnod)
if (flSolute)        call state%solute%init(config%solute, dims)
```

Rationale: matches SWAP's existing convention (`if (flSolute) call solute_init(state)` is already there); keeps init bodies simple and unit-testable in isolation; makes activation order explicit at the orchestration site.

## Pilot scope: surfacewater only

This spec covers the **pilot subsystem** only. Surfacewater is chosen because it exhibits all four "spread" loci from the Context section (cohort allocation, transient buffer, inline seed, readswap-style derived math) — making it the most instructive case to validate the pattern before rolling out to other subsystems.

### In scope

- **Define `surfacewater_state_t%init` on the state module.** Type-bound procedure in `src/state/surfacewater_state.f90`. Signature:
  ```fortran
  subroutine surfacewater_state_init(self, config_sw, config_drain, numnod)
     class(surfacewater_state_t),   intent(inout) :: self
     type(surface_water_config_t),  intent(in)    :: config_sw
     type(drainage_config_t),       intent(in)    :: config_drain
     integer,                       intent(in)    :: numnod
  end subroutine surfacewater_state_init
  ```
- **Body content** consolidates everything currently in:
  - `src/drainage/surfacewater_init.f90:27-145` (allocation + zero + config seed + sttab build + swst init)
  - The post-init block in `src/drainage/surfacewater.f90:66-93` inside `SurfaceWater(task=1)`: `hwlman = 0`, `vtair = 0`, and the `ZDraBas` initialization from drainage config and `state%surfacewater%wlstar`.
- **Body reads from typed config slices, not legacy globals.** Previously: `use variables, only: nrlevs, numnod, swdtyp, zbotdr, widthr, taludr, l, wls1_init, wlp, swsrf, swsec, swqhr, swman, nmper`. After: read these from `config_sw` and `config_drain` directly.
- **Inline the `wls1_init` computation.** Previously `config_to_variables.f90:1133` computes `wls1_init = config%surface_water%wlact - config%drain%altcu` into a module-level transient buffer; `surfacewater_init.f90:72` reads it. New `state%surfacewater%init` body computes inline:
  ```fortran
  self%wls    = config_sw%wlact - config_drain%altcu
  self%wlstar = self%wls
  ```
- **Hoist the call site.** Remove `call surfacewater_init(state)` from `src/drainage/surfacewater.f90:66` (inside `SurfaceWater(task=1)`). Add `if (flSurfaceWater) call state%surfacewater%init(config%surface_water, config%drain, numnod)` directly to the S4 pass in `swap_init` (`src/core/swap_mod.f90`). Also remove the post-init block at `surfacewater.f90:70-93` since `init` now subsumes it. `SurfaceWater(task=1)` case becomes a stub (`return`) with a comment.
- **Retire the `wlp = 0` global write.** This was a write to a legacy global made by `surfacewater_init` (since `wlp` is the primary surface-water level and surfacewater's branch is `swsrf=2` with no primary). After hoisting, this write must remain *for now* (other physics still reads `wlp` via `use variables`) — keep the write inside the new state init body until a follow-up arc retires the global. Mark with a comment.
- **Delete `src/drainage/surfacewater_init.f90`.** The module is now empty. Update `meson.build` to drop the file.
- **`wls1_init` transient buffer is left orphaned by this pilot.** The new `init` body inlines `wlact - altcu` directly from typed config and no longer reads `wls1_init`. After `surfacewater_init.f90` is deleted, the write at `config_to_variables.f90:1133` and the declaration in `variables.f90` become orphans (the value is computed and stored but never read). Removing them requires touching `config_to_variables.f90`, which is explicitly out of scope for this pilot (a separate arc will retire the transient buffers + adapter rows together). Audit: `grep -rn 'wls1_init' src/` after the pilot should show only the orphan write and declaration; surfacewater code no longer references it.

### Out of scope

- **All other subsystems** — solute, soilwater, atmosphere, tillage, heat, drainage, timecontrol. Their existing `_init` procedures stay. They migrate in a separate arc once the pattern is validated.
- **`config_to_variables.f90` retirement** — the surface_water section of the adapter still populates legacy globals (`swsrf`, `nmper`, `swdtyp`, etc.) that other physics modules read (`bocodre`, `SurfaceWater(2)`, output). Stopping those writes requires migrating the readers — a separate arc that the user plans to dispatch with a different agent.
- **Allocation consolidation** between `surfacewater_init` and the drainage.f90:445-467 inline allocator block — pre-existing and explicitly deferred from the cohort-flattening arc; deferred further.
- **Module-level pointer binding cleanup** (`bind_cofgen_target`, `bind_state_targets`, `bind_tc_target`) — a separate architectural concern.
- **`SurfaceWater(task=1)` deletion** — the case stays as a stub. Whether to remove the task=1 case entirely is left for a small follow-up cleanup.
- **Cross-subsystem state derivation patterns** (e.g., a hypothetical `state%heat%init` that depends on `state%soilwater%cofgen`) — surfacewater is self-contained, doesn't surface this question.
- **`Initialize()` retirement** — the legacy globals init routine. Out of scope.

## Architecture details

### Module dependencies after the pilot

`src/state/surfacewater_state.f90` gains these `use` statements:

```fortran
use surface_water_config_mod, only: surface_water_config_t
use drainage_config_mod,      only: drainage_config_t
use surfacewater_utils,       only: swstlev
use error_mod,                only: fatalerr_collected
```

`swstlev` is a small (~30 LoC) helper in `surfacewater_utils` that linearly interpolates the `sttab` storage table; it already exists. The new `init` calls it for the `swst = swstlev(self, wls1)` step.

`surface_water_config_mod` and `drainage_config_mod` are pure type-definition modules — no compute dependencies, no module-import cycles introduced.

### Validation defenses

The current `surfacewater_init` body has two `fatalerr_collected` guards:

```fortran
if (swsrf == 3 .or. swsec == 1 .or. swqhr == 2) then
   call fatalerr_collected('surfacewater_init', ...)
   return
end if
if (any(swman(1:nmper) == 2)) then
   call fatalerr_collected('surfacewater_init', ...)
   return
end if
```

These are belt-and-braces — `surface_water_config_validate` already rejects these branches. The new `init` retains both guards (defense in depth — a future change to the validator shouldn't silently break this).

### Call ordering in `swap_init`

The S4 pass for surfacewater requires:
- `CalcGrid()` has run (to populate `numnod` — currently sourced from `variables`; the new signature takes it as an explicit argument).
- `config` is validated + finalized.

The pilot inserts the call at the same orchestration point where `if (flSurfaceWater) call SurfaceWater(1, ...)` lives today (`swap_mod.f90:161`), but as an explicit S4 call before the physics pass. Approximate ordering:

```fortran
! ... existing code through CalcGrid() and current `_init` calls ...

! NEW — S4 surface-water init (hoisted out of SurfaceWater(1)):
if (flSurfaceWater) call state%surfacewater%init(config%surface_water, config%drain, numnod)

! ... rest of existing Phase-1 chain, with SurfaceWater(1)'s case(1) now a stub:
if (flSurfaceWater) call SurfaceWater(1, state, request_smaller_dt)   ! no-op stub
```

(Whether to keep the `SurfaceWater(1, ...)` call at all is a small follow-up; the pilot keeps it as a stub to minimize blast radius.)

### Tests

The existing pFUnit test `tests/unit/state/test_surfacewater_state.pf` covers default-zero scalars and array-unallocated state. Add new tests covering the `init` procedure:

- `test_surfacewater_init_allocates_cohort_arrays` — call `init` with a minimal `config_sw` + `config_drain` + `numnod`; assert `cqdrain`, `cqdrainin`, `cqdrainout`, `inqdra`, `inqdra_in`, `inqdra_out` are allocated to the right sizes.
- `test_surfacewater_init_sets_wls_from_config` — assert `state%wls == config_sw%wlact - config_drain%altcu` after init.
- `test_surfacewater_init_builds_sttab` — assert the storage table has the expected shape (rows 1 = +100, row 2 = 0, rows 3-22 dividing `[0, zbotdr(1)]` into 20 compartments).
- `test_surfacewater_init_computes_swst_from_sttab` — assert `state%swst == state%swstini == swstlev(state, state%wls)`.
- `test_surfacewater_init_sets_zdrabas_for_swsec_2` — assert ZDraBas equals `state%wlstar` when `swsec == 2`.

Plus the existing `tests/unit/io/toml/test_surfacewater_parity.pf` continues to exercise the full TOML→state pipeline.

## Verification gates

- **Build clean** — `pixi run build-linux`, no new warnings.
- **pFUnit green** — full suite passes, including the new init tests.
- **check-full byte-identical** — regression suite (hupselbrook, surfacewater, grassgrowth, salinitystress, oxygenstress) produces bit-identical output to the pre-pilot baseline. **This is the critical gate** — surfacewater is the namesake case and any drift will surface here immediately.
- **Audit greps** — `grep -rn 'surfacewater_init' src/ tests/` returns only references in this design's references / commit messages, not in source code. `grep -rn 'wls1_init' src/` returns zero hits.

Per project workflow convention, run `check-full` **before each commit**, not just at end-of-arc. pFUnit alone misses global-default regressions.

## Out of scope (recap)

- Other subsystems' state inits.
- `config_to_variables.f90` retirement.
- Allocation consolidation (drainage.f90 inline allocator).
- Module-level pointer-binding cleanup.
- `SurfaceWater(task=1)` dispatcher case removal.
- `Initialize()` legacy globals retirement.

## Consequences

(+) Single explicit `state%surfacewater%init(...)` call in `swap_init`. The init story for surfacewater becomes readable at the orchestration level instead of buried inside `SurfaceWater(task=1)`.
(+) Retires the `wls1_init` transient buffer. The strangler-fig bridge for surfacewater is gone.
(+) Establishes the pilot pattern (type-bound `init` on state module; driver-side flag guard; encapsulated body; reads typed config slices) for the remaining 6 subsystems' future migration.
(+) `surfacewater_init.f90` (a separate 148-LoC module that exists solely to host this work) is retired — code moves to where it logically belongs.
(+) Init body becomes unit-testable in isolation: pFUnit can construct typed config inputs and exercise the init without the full SWAP global-flag environment.

(–) `surfacewater_state.f90` grows from 121 to ~250 LoC (adding init body + new imports + post-init blocks). The state module is no longer "pure data"; it owns its lifecycle.
(–) Adds two cross-section config dependencies (`drainage_config_t` alongside `surface_water_config_t`) to the state module's import list. Explicit dependency surfacing — but a slight increase in the module-import graph.
(–) The dispatch idiom `SurfaceWater(task=1)` becomes inconsistent with peer subsystems for one commit's worth of code (case(1) is now a stub; case(2) and case(3) are real physics). This is intentional pilot behavior — the inconsistency motivates eventual cleanup of the task=1 dispatcher across all subsystems.
(–) Leaves an orphaned `wls1_init` write + declaration as discussed above. The clean-up requires touching `config_to_variables.f90`, which is reserved for a separate follow-up arc.

## References

- 2026-05-12 cohort-flattening spec: `docs/superpowers/specs/2026-05-12-flatten-reset-cohorts-design.md`
- ADR 0030 — state-migration surfacewater pilot (original migration)
- ADR 0033 — cumulative reset cohorts (superseded by ADR 0042 for the cohort pattern)
- ADR 0042 — flatten reset cohorts (immediately preceding arc)
- Current `surfacewater_init`: `src/drainage/surfacewater_init.f90`
- Current call site: `src/drainage/surfacewater.f90:66` (inside `SurfaceWater(task=1)`)
- Hoisting target: `src/core/swap_mod.f90` S4 pass
