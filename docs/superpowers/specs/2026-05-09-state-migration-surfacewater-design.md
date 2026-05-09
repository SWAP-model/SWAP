# Surface-Water State-Type Migration — Design

**Date:** 2026-05-09
**Status:** accepted (pending implementation)
**Pilot for:** state-type migration playbook (first subsystem migration arc; pattern carries over to subsequent subsystems)
**Discovery input:** `docs/superpowers/specs/2026-05-08-state-migration-surfacewater-discovery.md`
**Successor ADR:** 0030

## Goal

Migrate the surface-water subsystem off `variables.f90` module globals into a typed `surfacewater_state_t`, threaded through subroutine signatures from `swap_main` down. After this arc, the subsystem owns no globals — all owned variables live in the typed state. This is the pilot for a subsystem-by-subsystem migration that ultimately empties `variables.f90` and exposes a `swap_state_t` as the binary's primary stateful boundary (the surface a future Python embedding will see).

## Out of scope

- The other ~6–10 subsystems. They follow this playbook in subsequent arcs.
- Deleting `swapoutput.f90`. That's a sister arc — it removes hazards (3) and (4) below for free, but doesn't gate the state migration. Coordinated, not bundled.
- C interoperability layer / Python bindings. Future arc.
- Re-entrancy beyond what argument-threading naturally provides (single-instance is acceptable for MVP).

## Context

After 9 weeks of TOML-pipeline work, the input side is fully typed: TOML → typed config → adapter → globals. The output side is mixed (CSV path uses an internal typed struct; the older ASCII path goes globals → file directly). Compute reads/writes globals throughout. `variables.f90` is 1314 LoC of flat module-level globals (~1200 of them) — a single shared mutable state surface that prevents multi-instance use, blocks Python embedding via a clean state object, and couples subsystems implicitly.

The plan is to migrate subsystem-by-subsystem onto typed state types, threaded as arguments. Each subsystem's compute, output, and config converge on a single state slice. `variables.f90` shrinks per arc until empty.

The discovery doc catalogs surface-water as a peripheral, well-bounded subsystem (1433 LoC, 29 owned globals, 54 borrowed globals, 6 entry points). Surfacing the cross-subsystem hazards revealed 8 issues; this design resolves each.

## Decisions

### D1. Pilot scope: surface-water subsystem only

The four home files (`surfacewater.f90`, `surfacewater_init.f90`, `surfacewaterutils.f90`, the typed config + TOML reader) plus call-site adjustments wherever those files are invoked. No incidental refactoring of other subsystems.

### D2. State-type shape

Single derived type `surfacewater_state_t` aggregating all 29 owned globals, declared in a new module `src/state/surfacewater_state.f90`. Field names match the legacy global names (case-preserving) so `ASSOCIATE` blocks in compute bodies don't churn variable names.

```fortran
module surfacewater_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   public :: surfacewater_state_t

   type :: surfacewater_state_t
      ! per-step / per-day state
      real(real64) :: wls           = 0.0_real64    ! surface water level (cm)
      real(real64) :: wlstar        = 0.0_real64    ! target surface water level (cm)
      real(real64) :: swst          = 0.0_real64    ! storage per unit area (cm)
      real(real64) :: swstini       = 0.0_real64    ! initial storage (cm)
      real(real64) :: hwlman        = 0.0_real64    ! pressure head for target level (cm)
      real(real64) :: vtair         = 0.0_real64    ! total air volume (cm) — discovery §8 flagged this as misfiled in variables.f90 (declared in soil-water region but written exclusively by surface-water); it IS surface-water-owned
      real(real64) :: wlsold        = 0.0_real64    ! previous-step surface level (cm)
      real(real64) :: cqdrd         = 0.0_real64    ! cumulative drain into reservoir (cm)
      real(real64) :: cwsupp        = 0.0_real64    ! cumulative external supply (cm)
      real(real64) :: cwout         = 0.0_real64    ! cumulative outflow (cm)
      real(real64) :: cqdra         = 0.0_real64    ! cumulative lateral drainage (cm)
      real(real64) :: ZDraBas       = 0.0_real64    ! drainage basis level (cm)
      real(real64) :: iqdra         = 0.0_real64    ! intermediate lateral drainage total (cm)
      real(real64) :: qdrtot        = 0.0_real64    ! total lateral drainage flux (cm/d)
      logical      :: overfl        = .false.       ! automatic weir overflow flag
      logical      :: flInitDraBas  = .true.        ! init drainage basis (macropore)
      integer      :: imper         = 1             ! current management period index
      integer      :: numadj        = 0             ! count of target-level adjustments
      real(real64) :: wlsbak(4)     = 0.0_real64    ! 4-step circular buffer for oscillation detection

      ! per-level (Madr-sized) arrays — allocatable, sized at init from drainage config
      real(real64), allocatable :: cqdrain(:)        ! cumulative drainage per level
      real(real64), allocatable :: cqdrainin(:)      ! cumulative infiltration per level
      real(real64), allocatable :: cqdrainout(:)     ! cumulative outflow per level
      real(real64), allocatable :: qdra(:,:)         ! lateral drainage flux per level/compartment

      ! per-level intermediates (Madr × macp arrays)
      real(real64), allocatable :: inqdra(:,:)
      real(real64), allocatable :: inqdra_in(:,:)
      real(real64), allocatable :: inqdra_out(:,:)

      ! pre-computed level-storage table
      real(real64) :: sttab(22, 2)  = 0.0_real64
   end type surfacewater_state_t
end module surfacewater_state_mod
```

`l(Madr)` (drain spacing) is **not** in this struct — it's a drainage-config field, not surface-water state. See D6.

`fldecdt` and `qdrain(Madr)` are **not** in this struct — they belong to other subsystems. See D5 and D7.

### D3. State-type aggregation: build `swap_state_t` from day 1

Create `src/state/swap_state.f90` declaring:

```fortran
module swap_state_mod
   use surfacewater_state_mod, only: surfacewater_state_t
   implicit none
   private
   public :: swap_state_t

   type :: swap_state_t
      type(surfacewater_state_t) :: surfacewater
      ! Subsequent subsystem migrations add fields here.
      ! Non-migrated subsystems continue to use variables.f90 globals.
   end type swap_state_t
end module swap_state_mod
```

Each subsequent subsystem migration adds one field. No placeholders for un-migrated subsystems.

### D4. Argument-threading from `swap_main` down

`swap_main` declares one `type(swap_state_t) :: state` and passes it to subsystem dispatchers. Surface-water's entry points become:

```fortran
subroutine SurfaceWater(task, state, config, request_smaller_dt)
   integer,                   intent(in)    :: task
   type(swap_state_t),        intent(inout) :: state
   type(swap_config_t),       intent(in)    :: config
   logical,                   intent(out)   :: request_smaller_dt   ! see D5
end subroutine
```

Inside the routine: `ASSOCIATE(sw => state%surfacewater)` makes the body read/write `sw%wls` etc. with minimal churn. Borrowed reads use `state%water%theta`, `state%timecontrol%t1900`, etc. — but during this arc those subsystems aren't migrated yet, so borrowed reads still come from `variables.f90` globals via narrow `use variables, only: …` clauses. The `use variables` clause shrinks as more subsystems migrate; for surfacewater specifically, the *owned* globals leave (~29 symbols) but borrowed (~54) stay until their owners migrate.

`runoff()` (in `surfacewaterutils`) becomes `runoff(sw, …other-args, runoff_value)` — explicit state argument. `boundtop.f90` (its caller) gets `state` passed in by *its* caller, threading down from `swap_main`.

### D5. `fldecdt` becomes an `intent(out)` flag

Hazard (5) resolution. The "decrease timestep" signal leaves the state struct entirely. `SurfaceWater(task=2, …)` adds `intent(out) :: request_smaller_dt`. `swap_main` reads it after the call and updates `dt` accordingly. Self-documenting at the call site, no shared-state coupling.

The legacy global `fldecdt` is removed from variables.f90 in Phase 2 (alongside the rest of the owned globals); during Phase 1 it stays as a transitional shim that the adapter writes into.

### D6. `l(Madr)` m→cm conversion moves to config-load

Hazard (6) resolution. Drainage TOML config currently expresses `l` in metres; the legacy reader and `surfacewater_init` together convert to cm in-place. After this arc, the drainage typed config either (a) declares `l` in cm directly with the TOML reader doing the conversion at load, or (b) carries an explicit `l_unit` field. Choose (a) — simpler, no per-call ambiguity. The conversion happens once, in `read_drainage_toml.f90` (or wherever the drainage config's `apply_*` adapter sits). Subsystems read pre-converted values.

The legacy global `l(Madr)` keeps existing for now (drainage subsystem hasn't migrated yet); only its *units* contract changes (from "value depends on whether init has run" to "always cm post-load"). Surfacewater stops mutating it.

### D7. `qdrain(Madr)` rule moves to drainage owner

Hazard (1) resolution. The `if (gwl > 998) qdrain(:) = 0.0` rule moves from `surfacewater.f90:117-120` to `drainage.f90/bocodre`. Surface-water no longer references `qdrain` at all. Drainage stays globally-driven for now (it hasn't migrated), so the rule is implemented against globals — but the *ownership* is corrected. When drainage migrates later, the rule moves cleanly into `drainage_state_t`'s update routine.

This is a small per-step change in `bocodre`; we verify check-full byte-identical after the relocation.

### D8. `wls ↔ runots` cycle: stays explicit, no architectural change

Hazard (7) resolution. SWAP's existing timestep-iteration loop already converges this cycle implicitly — making the dependencies explicit doesn't break the convergence. After argument-threading: `runoff(state%surfacewater, state%boundary%runots, …)` is called from boundtop, modifying `runots`; surfacewater's per-step uses the updated `runots`. The cycle is identical to today; it just stops hiding behind module globals.

### D9. Output coupling: read-only state-aware paths

Hazards (3) and (4) (output reaching back into state, compute called from output) live inside `swapoutput.f90`. We do **not** delete swapoutput.f90 in this arc — that's a sister arc. But we DO update its surface-water reads to go through the typed state when the relevant fields have moved. Concretely:
- `outdrf`, `outswb`, `outend` open `use variables, only: …` for surfacewater-owned globals → those `use` clauses break once we remove the globals. We replace each with reads from the typed state, threaded in: `subroutine outswb(task, state)` and so on.
- The `outswb` year-end mutation of `swstini` (hazard 4) moves to a dedicated `surfacewater_year_reset(state%surfacewater)` callback invoked from the daily/yearly loop in `swap_main`. The output routine becomes read-only.
- The `SurfaceWater(2)` call from inside `swapoutput.f90` (hazard 3, sensitivity loop) is removed — the sensitivity loop is dead in the modern flow (all outputs are file-based; this is leftover from interactive-debug code). Verify by running with `swend=2` after deletion.

CSV output (`set_values` in `swap_csv_output.f90`) gets the same treatment for the surface-water-related variables it currently reads from globals.

### D10. Config: 12 missing fields — DEFERRED (separate arc)

The discovery doc Section 8 item (8) lists 12 fields used by stub-errored branches (`swman=2`, `swsec=1`, `swqhr=2`) that the legacy reader populated but `surface_water_config_t` doesn't declare: `wlsman`, `gwlcrit`, `nphase`, `nodhd`, `dropr`, `vcrit`, `hcrit`, `hqhtab`, `qqhtab` plus three siblings the discovery doc enumerates in detail.

**Status: deferred.** On a fresh pass, none of the 29 owned state-type fields overlap with these 12 config inputs, and the regression cases never exercise `swman=2` / `swqhr=2` / `swsec=1` (all stub-errored). Promoting these fields therefore does not unblock Phase 1 or Phase 2. It also drags in TOML-schema design (2D `mamp × mamte` tables — nested arrays vs sub-tables vs CSV companion) that the spec doesn't pin down.

Track as its own future arc, kicked off whenever someone needs to enable `swman=2` / `swqhr=2`. Until then, the stub-errors keep these branches inactive and `surface_water_config_t` stays at its current shape.

## Phasing

The arc decomposes into two phases with check-full as the gate between them:

- **Phase 1 — State type + per-step migration.** Define `surfacewater_state_t` and `swap_state_t`. Thread `state` through all surfacewater entry points and call sites. Move `fldecdt` to `intent(out)`. Apply ASSOCIATE inside compute bodies. Update output-side reads to use threaded state. Verify check-full at end of phase.
- **Phase 2 — Owner-rule relocation.** Move `qdrain` zeroing rule to drainage's `bocodre`. Move `l(Madr)` m→cm conversion to drainage config load. Remove dead `SurfaceWater(2)` call from `swapoutput.f90`. Remove `fldecdt`, the 29 owned globals, and any incidentally-orphaned declarations from `variables.f90`.

(Phase 0 — promote 12 missing config fields — is deferred. See D10.)

Each phase ships independently with check-full byte-identical and pFUnit green as the integration gates.

## Testing

- **pFUnit:** existing surfacewater tests (`test_surfacewater_init.pf`, `test_surface_water_config.pf`, `test_read_surface_water_toml.pf`, `test_surfacewater_parity.pf`) must stay green. Migrate them to use the typed state where they reference globals. Add new tests for the lifecycle entry points (`SurfaceWater(task=1)` initializes state correctly; `task=2` updates state and signals `request_smaller_dt` correctly under oscillation; `task=3` finalizes cleanly).
- **check-full:** byte-identical against `swap.ok` for all 5 cases. This is the integration gate at every phase boundary.
- **Hand-verification:** none required — surfacewater is exercised by the `surfacewater` regression case under `tests/swap-cases/toml/6.surfacewater/`.

## Non-goals (explicitly)

- Rewriting compute logic. ASSOCIATE preserves variable names exactly; the body of `WLEVBAL` etc. should not change semantically.
- Adding new physics features.
- Generalizing the state-type machinery before we have a second subsystem migration to learn from. YAGNI — `swap_state_t` starts with one field, grows as we go.

## ADR 0030

`docs/adr/0030-state-migration-surfacewater-pilot.md` records:
- The decision to thread typed state from `swap_main` down (not module-level intermediate).
- The single-owner principle for cross-subsystem variables.
- The hazard-by-hazard resolution table.
- The migration playbook structure (discovery → design → plan → execute) for subsequent subsystems.

## Reusable playbook artifacts produced by this arc

For subsequent subsystems (water-balance, crop, drainage, atmosphere, …):

- The discovery template (already populated for surfacewater).
- The state-type module pattern (`src/state/<subsystem>_state.f90`).
- The argument-threading pattern (entry points take `state`, `config`, optional `intent(out)` signals).
- The hazard taxonomy (single-owner, intent-out signal, config-load conversion, ownership relocation).
- The phasing pattern (config-completion → state-migration → owner-rule-relocation).

Each subsequent subsystem follows the same four-phase shape: discovery → design → plan → execute. Each gets its own ADR. Each adds one field to `swap_state_t`.
