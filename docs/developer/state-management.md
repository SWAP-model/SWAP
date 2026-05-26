---
title: State management
author: SWAP modernization team
---

# State management

## Context

The legacy entry point for shared model data is `src/core/variables.f90`, a
single Fortran module that declares well over a thousand global mutable
symbols — scalars, arrays, and flags that every physics routine `use`s and
writes into as a side effect. This design has three concrete problems for the
rescue. First, the dependency graph between routines is invisible: a caller
cannot tell from a subroutine's argument list what state it reads or
produces, because both happen through `use variables`. Second, the module
makes the code test-hostile; a pFUnit test cannot set up the subset of state
a routine actually needs without dragging in the full global set, nor can it
isolate two independent tests that happen to touch overlapping globals.
Third, the globals preclude any form of multi-instance or thread-parallel
execution — two SWAP runs in the same process would race on the same
storage, so BMI coupling, ensemble drivers, and future GPU offload are all
blocked until the globals are gone.

The modernization replaces this with explicit state passing. Each physics
domain owns a derived type, all domain types are composed into a single
aggregate `swap_state_t`, and routines receive the slice of state they need
by reference with explicit `intent`. Phase 2 of the rescue establishes the
aggregator as the canonical pattern; Phase 4 removes the last readers of the
legacy `variables` module. The architectural decision to prefer a single
aggregator over per-compartment state trees is expected to be recorded in
ADR 0003 (added by a later Phase 2 task).

## Three type categories

Conceptually, every datum carried by SWAP falls into one of three
lifetime/mutability categories. The convention below is the target for Phase
3–4; the current rescue-era code has all three merged into each domain's
`*_state_t`, and the split into three concrete types is a refinement that
will be introduced incrementally.

**`config_t` — read-only configuration.** Fields populated once from the
configuration files (the `.swp` TOML or legacy input, the `.dra` drainage
input, the `.crp` crop parameterizations). Method selectors (`swetr`,
`swinter`, `dramet`), Angstrom coefficients, file paths, and similar knobs
live here. Downstream physics routines take `config` with `intent(in)` and
are forbidden from writing to it. Pinning immutability at the argument list
makes configuration a real type-level contract rather than a convention.

**`initial_t` — derived-at-init, frozen thereafter.** Quantities computed
once at start-of-simulation from `config_t` but not part of the dynamic
time-step loop. The numerical grid (compartment thickness, node depths,
layer-to-node mapping) built by `CalcGrid_state` is the canonical example:
it depends on configuration but is constant for the whole run. `initial_t`
is built during Task 1 and read with `intent(in)` during Task 2.

**`state_t` — mutated each step.** Water content, pressure head, leaf area
index, accumulated fluxes, solute concentrations — anything the time loop
updates. This is the only category that legitimately moves.

A given domain may have all three (soil will, once fully refactored: config
from `.swp`, initial grid, dynamic water/solute state), two
(atmosphere essentially has config and state but no meaningful initial),
or just one (the boundary module presently exposes only a single state).
Today the aggregator's fields are all `*_state_t` for uniformity; the
config/initial split is applied field-by-field inside them and is
mechanically extractable into separate types later.

## Canonical module layout

For a domain named `foo`, the source file is `src/foo/foo_state.f90`. It
contains a single module `foo_state_mod` with the following public surface,
matching the existing implementations in `src/atmosphere/atmosphere_state.f90`
and `src/drainage/drainage_state.f90`:

```fortran
module foo_state_mod
    implicit none
    private

    public :: foo_state_t
    public :: foo_state_init
    public :: foo_state_finalize
    public :: foo_state_reset_cumulative
    public :: foo_state_reset_intermediate

    type :: foo_state_t
        ! scalar fields with default initializers
        ! allocatable arrays sized in init
    end type foo_state_t

contains
    subroutine foo_state_init(state, ...)
        type(foo_state_t), intent(inout) :: state
        ! allocate arrays, call foo_state_reset_all
    end subroutine foo_state_init

    subroutine foo_state_finalize(state)
        type(foo_state_t), intent(inout) :: state
        ! deallocate arrays; optionally call foo_state_reset_all
    end subroutine foo_state_finalize
end module foo_state_mod
```

The reset trio that the established domains (atmosphere, drainage, soil)
all expose carries a specific meaning:

- `foo_state_reset_all(state)` — private helper called by both `_init` and
  `_finalize`; zeroes every field. Used so the defaults from the type
  declaration are not relied on implicitly.
- `foo_state_reset_cumulative(state)` — zeroes fields that accumulate over
  a reporting period (`cgrai`, `cevap`, `cqdrain(:)` …). Invoked when
  `flzerocumu` fires in the time loop.
- `foo_state_reset_intermediate(state)` — zeroes intra-day intermediate
  accumulators (`inrai`, `ievap`, `inqdra(:,:)` …). Invoked when
  `flzerointr` fires.

Bridges to the legacy globals — `foo_state_from_variables(state, ...)` and
`foo_state_to_variables(state, ...)` — live in `src/core/swap_state_sync.f90`
rather than in each domain module (see Legacy bridge below), because they
need to `use variables` and we want the domain modules to stay free of that
coupling.

## Aggregation

`src/core/swap_state_mod.f90` declares `type :: swap_state_t`, which
composes one instance of each domain's state type as a named field:

```fortran
type :: swap_state_t
    type(time_state_t)         :: time
    type(soil_state_t)         :: soil
    type(atmosphere_state_t)   :: atm
    type(crop_state_t)         :: crop
    type(drainage_state_t)     :: drain
    type(boundary_state_t)     :: boundary
    type(macropore_state_t)    :: macro
    type(solute_state_t)       :: solute
    type(heat_state_t)         :: heat
    type(snow_state_t)         :: snow
    type(surfacewater_state_t) :: surfwater
    type(wofost_soil_state_t)  :: wofost_soil
    type(oxygenstress_state_t) :: oxystress
    type(irrigation_state_t)   :: irrig
    type(tillage_state_t)      :: tillage
    integer :: numnod = 0
    integer :: numlay = 0
    integer :: nrlevs = 0
    integer :: ncrop = 0
    logical :: initialized = .false.
end type swap_state_t
```

The driver (`src/core/swap_main.f90`) holds the `swap_state_t` instance and
passes it by reference to the three-phase entry point `swap(iCaller, iTask,
state, toswap, fromswap)`. Physics routines receive the slice they actually
need — `SoilWater_state(state, task)` takes the full aggregator for now,
but the target after Phase 3 is that leaf-level routines take only
`state%soil` or `state%drain`, keeping dependency edges explicit and
minimal rather than passing "all of the state to all of the routines."

## Scoped aliasing with ASSOCIATE

Long `%`-chains through the aggregator read poorly. Ported physics
routines open an `associate` block at the top of the body and bind each
sub-record they touch to a short canonical alias (`soil`, `drai`, `time`,
`atmo`, …). The same name is used for the same sub-record across every
file, so `soil%gwl` means `state%soilwater%gwl` no matter which routine
you are reading.

See [state-aliasing.md](state-aliasing.html) for the canonical alias
table, a worked drainage example, and the rules on when to use the
pattern versus extracting a helper subroutine.

## Lifecycle contract

A `swap_state_t` instance goes through exactly one well-defined lifecycle
per simulation:

1. **Init — exactly once.** The main program calls `swap(iTask=1, state, …)`,
   which reads configuration, sizes the grid, and then invokes
   `swap_state_init(state, numnod, numlay, nrlevs, ncrop)`. That routine
   calls each domain's `*_state_init` in turn (see lines 724–775 of
   `swap_state_mod.f90`). After init, `state%initialized` is `.true.` and
   every allocatable is allocated.

2. **Step — arbitrarily many times.** `swap(iTask=2, state, …)` advances
   one day. Within the day, the time-control logic may raise `flzerocumu`
   or `flzerointr`, causing `*_state_reset_cumulative` /
   `*_state_reset_intermediate` to fire on the relevant domains.

3. **Finalize — exactly once.** `swap(iTask=3, state, …)` calls
   `swap_state_finalize(state)`, which calls each domain's `*_state_finalize`
   to deallocate arrays. `state%initialized` is set back to `.false.`.

The invariant: at any observable program point, `state` is either
un-initialized, fully-initialized, or finalized — never partial. Calling a
physics routine before init, or between finalize and a subsequent init, is
a programmer error and produces undefined behaviour. During the rescue the
contract is by convention; Phase 4 adds checked invariants (assertion on
`state%initialized`, shape checks on allocatables).

## How to add a new domain state

To introduce a new domain `foo`:

1. **Create the module.** Add `src/foo/foo_state.f90` with
   `module foo_state_mod`, `type :: foo_state_t`, and the public
   `foo_state_init(state, …)` and `foo_state_finalize(state)` routines.
   If the domain accumulates within-day or across-day fluxes, add
   `foo_state_reset_intermediate` and `foo_state_reset_cumulative` mirroring
   the atmosphere/drainage pattern. The `_init` routine should allocate
   arrays and then delegate zeroing to a private `foo_state_reset_all`.

2. **Compose it into the aggregator.** In `src/core/swap_state_mod.f90`,
   add `use foo_state_mod, only: foo_state_t, foo_state_init,
   foo_state_finalize`, declare `type(foo_state_t) :: foo` as a field of
   `swap_state_t`, add `call foo_state_init(state%foo, …)` in
   `swap_state_init`, and add `call foo_state_finalize(state%foo)` in
   `swap_state_finalize`.

3. **Wire the three-phase driver.** In `src/core/swap.f90`, the
   aggregator's init inside the `iTask == 1` block already covers
   allocation — but if the domain also has a legacy `*_state(state, task)`
   entry point of its own (as `SoilWater_state`, `SurfaceWater_state`,
   `MacroPore_state`, etc. do), add a `call Foo_state(state, 1)` in the
   Task 1 block, a `call Foo_state(state, 2)` in the Task 2 loop, and
   `call Foo_state(state, 3)` in the Task 3 block if it needs a per-run
   wrap-up distinct from `_finalize`.

4. **Wire the build.** Add `src/foo/foo_state.f90` to the main
   `meson.build` sources list, and to the `test_base_sources` array in
   `tests/unit/meson.build` so pFUnit-driven unit tests that touch the
   new state compile against it.

5. **Optional — legacy bridge.** If the domain still shares data with the
   legacy `variables` module, add `foo_state_from_variables(state, …)` and
   `foo_state_to_variables(state, …)` procedures to
   `src/core/swap_state_sync.f90`, following the many existing examples
   there. These disappear at Phase 4 exit.

6. **Test.** When Phase 4 re-introduces `.pf` files, add a pFUnit
   lifecycle test under `tests/unit/` that exercises
   init → reset → finalize in sequence and asserts no leaks and no
   leftover non-zero fields.

## Legacy bridge

`src/core/swap_state_sync.f90` is rescue-era scaffolding. It is a single
~3300-line module that declares pair-wise copy procedures
(`state_from_variables(state)` / `state_to_variables(state)` at the top
level, plus per-domain `*_state_from_variables` / `*_state_to_variables`
and the narrower `*_outputs_from_variables` output-only variants). Each
procedure copies fields between the aggregator and the legacy
`variables.f90` globals.

The bridge exists because the port is gradual: some physics routines have
been updated to read from and write to `state`, while others still read
and write the globals. The convention during the dynamic phase is that
the driver pushes `state` into `variables` before calling a
not-yet-ported routine, and pulls changes back after. Routines that have
been ported skip the bridge and take `state` directly.

`swap_state_sync.f90` is temporary. Its line count is an explicit measure
of remaining legacy coupling; it shrinks monotonically as physics
routines are ported. It is deleted in its entirety at Phase 4 exit, along
with `src/core/variables.f90` itself.

## Reference

See [architecture.md](architecture.md) for the higher-level picture of
how `swap_state_t` fits into the three-phase control flow and the driver.
When ADR 0003 is added (a later Phase 2 task), it will record the decision
to use this single-aggregator pattern in preference to the originally
sketched per-compartment state tree.
