# Drainage State-Type Migration — Design

**Date:** 2026-05-10
**Status:** accepted (pending implementation)
**Migration #:** 2 of N (subsystem-by-subsystem state migration umbrella)
**Branch:** `refactor/surfacewater-state` (continuing the umbrella migration branch)
**Discovery:** `docs/superpowers/specs/2026-05-09-state-migration-drainage-discovery.md`
**Predecessor playbook:** `docs/superpowers/specs/2026-05-09-state-migration-surfacewater-design.md` + ADR 0030
**Successor ADR:** 0031

## Goal

Migrate the drainage subsystem off `variables.f90` module globals into a typed `drainage_state_t` aggregated under `swap_state_t`. After this arc, drainage-owned variables (including `qdra` and `qdrain` — currently mis-classified into `surfacewater_state_t` because of historical write-site overlap) live in `drainage_state_t`. `divdra` is modernized to take assumed-shape arrays so the legacy `qdra(Madr,macp)` global can be deleted. The pre-existing `geofac` adapter ordering bug at `config_to_variables.f90:431` is fixed inline.

## Out of scope

- The remaining ~6–8 subsystems. They follow this same playbook in subsequent arcs.
- Promoting any deferred config fields (no Phase 0 in this arc — discovery did not surface any unmodeled configs).
- Refactoring `divdra`'s internals beyond the signature change.
- Resolving the deeper `geofac` schema duplication (two TOML fields mapping to one legacy global). The inline fix preserves the existing semantics; full schema reconciliation is a separate doc-level decision.

## Context

The surface-water migration (Phase 1 + Phase 2, commits `acae2af..2e64a87`) established the playbook: discovery → design → plan → execute, with state-type definition, argument-threaded `swap_state_t`, ASSOCIATE in compute, dual-write transitional pattern, and check-full byte-identical at every commit. Drainage is migration #2.

The drainage discovery doc cataloged 6 owned globals, 52 borrowed, 6 external reader files, 9 cross-subsystem hazards. The owned set is much smaller than surfacewater's (29 owned), but several hazards involve the now-shipped surfacewater migration:

- `qdra` and `qdrain` are currently in `surfacewater_state_t`, written by drainage compute. This is a mis-classification we now correct.
- `divdra` keeps `qdra` alive as a legacy global because it takes the array as an explicit-shape `(Madr, macp)` argument — this is the reason a single global survived surfacewater Task 11.
- `ztopdislay` has near-duplicate redistribution code in both `drainage()` and `SurfaceWater(2)`.
- `frozencond.f90:FrozenBounds` writes `qdrain` AND calls `DIVDRA` with global arrays (3rd writer of `qdra`).

## Decisions

### D1. Pilot scope: drainage subsystem home tree

The four home files: `drainage.f90`, `divdra.f90`, `drainage_config.f90`, `read_drainage_toml.f90`. Plus call-site adjustments in callers (`swap.f90`, `surfacewater.f90`, `frozencond.f90`, `swapoutput.f90`, `swap_csv_output.f90`, etc.). No incidental refactoring of other subsystems beyond what the migration requires.

### D2. State-type shape

Single derived type `drainage_state_t` in `src/state/drainage_state.f90`:

```fortran
module drainage_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   public :: drainage_state_t

   type :: drainage_state_t
      real(real64) :: qdrd     = 0.0_real64    ! drain-direction sub-flux (cm/d)
      ! Per-level (Madr-sized) arrays, allocated by drainage_init from
      ! drainage config dimensions
      real(real64), allocatable :: qdrain(:)       ! (Madr) lateral drainage flux per level (cm/d)
      real(real64), allocatable :: drainl(:)       ! (Madr) drain length per level (cm)
      real(real64), allocatable :: wetper(:)       ! (Madr) wetted perimeter per level (cm)
      real(real64), allocatable :: ztopdislay(:)   ! (Madr) top of discharge layer per level (cm)
      ! Per-level / per-compartment array
      real(real64), allocatable :: qdra(:,:)       ! (Madr, macp) lateral drainage flux per level/compartment (cm/d)
   end type drainage_state_t

end module drainage_state_mod
```

Field names match legacy globals so ASSOCIATE blocks in compute bodies preserve variable names.

`qdra` and `qdrain` are **moved** from `surfacewater_state_t` to `drainage_state_t`. Per Q1 decision: drainage is the natural owner; surface-water consumes via `state%drainage%qdra` etc. This is a renaming refactor in surfacewater code (no semantic change — the same fields just live in a different state slice).

### D3. State-type aggregation: `swap_state_t` grows by one field

`src/state/swap_state.f90` gains:

```fortran
type :: swap_state_t
   type(surfacewater_state_t) :: surfacewater
   type(drainage_state_t)     :: drainage   ! NEW
end type swap_state_t
```

### D4. Argument-threading: extend the existing pattern

`Drainage(state)` already takes state (added in surfacewater Phase 2 Task 3 to thread surfacewater reads). The signature stays — the routine now writes its own `state%drainage%*` instead of writing `state%surfacewater%qdra`/`qdrain` and the legacy globals.

`divdra` (Q2 decision) gets a new signature: take `qdra` as an assumed-shape `qdra(:,:)` argument instead of explicit-shape `qdra(Madr,macp)`. This is the change that lets us delete the `qdra` global entirely. Verify that `divdra`'s body is compatible — assumed-shape arrays carry their own bounds via Fortran intrinsic `size()`, but indexed accesses inside the body are unchanged.

ASSOCIATE in compute bodies: bare-name aliases where `use variables, only:` doesn't shadow; `dr_*` prefix where it does. Same pattern as surfacewater used (`sw_*`).

### D5. Move `qdra` and `qdrain` from `surfacewater_state_t` to `drainage_state_t`

This is the cleanup of the technical debt acknowledged in ADR 0030. The mechanical move:

- Delete `qdrain` and `qdra` declarations from `surfacewater_state_t`.
- Add them (with the same shapes and semantics) to `drainage_state_t`.
- Replace every `state%surfacewater%qdra` / `state%surfacewater%qdrain` reference with `state%drainage%qdra` / `state%drainage%qdrain` across the codebase. Mostly mechanical — grep + edit.
- The allocation logic in `surfacewater_init` for these arrays also moves to a new `drainage_init` (or wherever drainage initialization happens — TBD by inspection of drainage's lifecycle entry points).

### D6. Modernize `divdra` to assumed-shape

`subroutine divdra(numnod, nrlevs, dz, ksatfit, …, qdrain, qdra, …)` currently has `qdrain(Madr)` and `qdra(Madr, macp)` as explicit-shape arguments. Change to assumed-shape:

```fortran
real(real64), intent(in)  :: qdrain(:)
real(real64), intent(out) :: qdra(:,:)
```

Inside the body, `Madr` and `macp` constants are replaced with `size(qdrain)` and `size(qdra, 2)`. This removes the dependency on the global `Madr`/`macp` parameter values.

After the modernization, `divdra` callers pass `state%drainage%qdrain` and `state%drainage%qdra` directly — no more legacy global staging.

### D7. Resolve `ztopdislay` dual-write

Drainage and SurfaceWater both write `ztopdislay`. Per Q3 decision: drainage is the single owner. Two paths:

- The `swdislay == 1` path runs in `Drainage()` (drainage.f90:494-537).
- The `swdislay == 2` path runs in `SurfaceWater(2)` (surfacewater.f90:146-204).

Move the `swdislay == 2` path into drainage. SurfaceWater stops writing `ztopdislay`. The two paths converge into a single drainage routine that handles both modes. Verify by reading the two blocks that they're truly near-duplicates (the discovery doc said so) and consolidating doesn't change behavior.

If consolidation is non-trivial (the two blocks differ in subtle ways), keep them separate but make drainage the writer in both cases — surface-water just calls into drainage at the right point.

### D8. Resolve `frozencond.f90:FrozenBounds` `qdrain` writeback + `DIVDRA` call

`FrozenBounds` modifies `qdrain` after drainage has computed it (frozen-soil correction), then calls `DIVDRA` again with the corrected `qdrain` to redistribute. After the migration:

- `FrozenBounds` takes `state` as `intent(inout)` (it already does — added during surfacewater Phase 2 Task 5).
- `qdrain` modifications go into `state%drainage%qdrain` (replacing global writes).
- The `DIVDRA` call passes `state%drainage%qdrain` and `state%drainage%qdra` (assumed-shape post-D6).
- Drop legacy global writes.

### D9. Fix `geofac` adapter ordering bug

At `src/io/toml/config_to_variables.f90:431`, the unconditional `geofac = config%drain%surface_runoff%geofac` overwrites the `ipos=5` write at line 306. The TOML schema has two distinct `geofac` fields (drainage geometry vs surface-runoff), but they map to a single legacy global.

Fix: gate line 431 on `if (config%drain%ipos /= 5)`:

```fortran
if (config%drain%ipos /= 5) then
   geofac = config%drain%surface_runoff%geofac
end if
```

This preserves the existing behavior (line 306's `ipos=5` write wins when applicable; the surface-runoff write wins otherwise). The deeper schema duplication is documented in the spec but not resolved here.

A test fixture exercising both branches should be present or added.

### D10. Lifecycle entry point: drainage_init

If drainage has no current init routine (it might be initialized inline in `Drainage(task=1)` if it has a task switch, or initialized by `surfacewater_init` for the per-level arrays), introduce `drainage_init(state)` as a public entry point analogous to `surfacewater_init`. Allocates `state%drainage` per-level arrays from drainage config dimensions.

Verify by inspection: if `surfacewater_init` currently allocates `qdrain` etc. (because surfacewater_state_t hosted them), that allocation moves to `drainage_init`.

## Phasing

The arc decomposes into two phases with check-full as the gate:

- **Phase 1 — State type, threading, dual-write, output reads.** Define `drainage_state_t`, add to `swap_state_t`, thread state through drainage entry points, dual-write `state%drainage%*` and legacy globals during compute, redirect output reads to typed state, drop dual-write from drainage home files. Move `qdra`/`qdrain` from surfacewater_state_t to drainage_state_t (mechanical refactor with check-full at every step). Modernize `divdra` signature.
- **Phase 2 — Owner-rule relocation, cross-subsystem reader migration, global removal.** Move `swdislay==2` `ztopdislay` rule into drainage. Migrate cross-subsystem readers (frozencond, surfacewater, output routines) to read from `state%drainage`. Drop legacy `qdrain`, `qdra`, `qdrd`, `drainl`, `wetper`, `ztopdislay` globals from `variables.f90`. Fix `geofac` adapter ordering inline. ADR 0031 finalize.

Each phase ships independently with check-full byte-identical and pFUnit green as the integration gates.

## Testing

- **pFUnit:** existing drainage tests stay green (test the discovery's Section 9 to enumerate). Migrate them to assert against `state%drainage%*` where applicable. Add a new pFUnit test for `drainage_state_t` (default values + array allocation lifecycle), matching `test_surfacewater_state.pf`. Add a test for the `geofac` adapter fix covering both `ipos==5` and `ipos!=5` cases.
- **check-full:** byte-identical against `swap.ok` for all 5 cases. Integration gate.

## Non-goals (explicit)

- Rewriting drainage compute logic. ASSOCIATE preserves variable names; the bodies of `bocodre`, `divdra`, etc. should not change semantically.
- Adding new physics features.
- Resolving the schema-level `geofac` duplication.

## ADR 0031

`docs/adr/0031-state-migration-drainage.md` records:
- Decision: drainage is migration #2; subsystems migrate one at a time on the umbrella branch.
- The `qdra`/`qdrain` re-classification (moved from surfacewater_state_t to drainage_state_t).
- The `divdra` assumed-shape modernization.
- The `ztopdislay` single-owner resolution.
- The `geofac` adapter ordering inline fix.
- Cross-references to the discovery doc, ADR 0030, and the playbook artifacts.
