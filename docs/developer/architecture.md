---
title: Architecture
author: SWAP modernization team
---

# Architecture

## Overview

SWAP (Soil–Water–Atmosphere–Plant) is a process-based model that simulates
transport of water, heat, and solutes in a one-dimensional vertical soil
column coupled to the atmosphere at the top and to a groundwater / drainage
system at the bottom. A crop canopy mediates the exchange: it intercepts
precipitation, partitions available energy into potential transpiration and
soil evaporation, and withdraws water from the root zone. The spatial
discretization is a stack of numerical compartments along the vertical axis;
no lateral extent is resolved explicitly. Lateral drainage to ditches and
tile drains is parameterized through aggregated flux expressions rather than
through a 2-D/3-D grid.

The physics covered includes Richards-equation unsaturated soil-water flow,
soil heat conduction with optional freezing effects, advection–dispersion
solute transport with simple reactive terms, optional preferential flow
through a macropore continuum, surface-water management for polder-style
systems, and crop growth via three alternative crop models. The outer
integration loop runs on a daily time step — meteo input and crop updates
advance one calendar day at a time — while the water and solute solvers
internally subdivide the day into adaptive sub-steps whenever convergence
requires it. This gives a clean separation between the agronomic cadence
(daily) and the numerical cadence (adaptive, sub-daily), with the daily
loop serving as the integration envelope for output and crop processes.

## Three-phase control flow

The model is factored into a thin driver and a three-phase entry point so
that SWAP can also be compiled as a shared library and invoked from a host
program. The driver is `src/core/swap_main.f90`, which owns the rerun loop:
it reads `reruns.dat`, iterates over the rerun sets, and for each set calls
the three-phase entry point `swap(iCaller, iTask, state, toswap, fromswap)`
defined in `src/core/swap.f90` with `iTask` equal to 1, 2, and 3 in turn.
The `swap_state_t` container (see below) is passed in and out on every
call, so the entry point is re-entrant with respect to its state.

Task 1 is initialization. It reads the input configuration via the
typed TOML pipeline (`ReadSwapToml_state` and the
`config_to_variables` adapter) — the legacy fixed-format reader was
deleted in 2026-05-06, so a `.toml` argument is now mandatory. It
initializes the
legacy `variables` module, builds the numerical grid via `CalcGrid_state`,
allocates the aggregated state via `swap_state_init`, and initializes each
physics domain (soil water, surface water, macropores, soil temperature,
snow, solute, age tracer, soil management) by calling the domain's
`*_state` routine with task=1. Output files are opened and headers are
written at the end of task 1.

Task 2 is the dynamic step — the per-day time loop. Within each day the
physics routines are called in a strict order so that each consumer reads
state already updated by its upstream producers: read meteo for the day
(`ReadMeteoDay_state`), evaluate crop growth phenology, compute irrigation,
process meteo (including daytime/nighttime splitting), apply tillage,
evaluate snow, compute frozen-soil conductivity reductions, compute
potential and actual root water extraction (`RootExtraction_state`),
determine the bottom boundary condition (`BoundBottom_state`), compute
drainage fluxes (`Drainage_state`), advance surface-water balance, solve
the Richards equation for soil water (`SoilWater_state`, which uses the
Thomas tridiagonal algorithm `tridag` in `src/utils/numericalsolvers.f90`),
step soil temperature, step solute and optional age-tracer transport,
handle macropore preferential flow as part of the soil-water sub-step,
and finally write output. Sub-daily time-step reductions are driven by the
`fldtreduce` flag inside the soil-water loop. Task 3 is closure: final
output, iteration-time statistics, and deallocation.

## State aggregation

The legacy code kept essentially all simulation state in `src/core/variables.f90`,
a single module of hundreds of global variables that every physics routine
read from and wrote to. That design made multi-instance execution, thread
safety, and any future GPU offloading impossible, and it made data flow
between physics routines invisible at the call site. The modernization
replaces this with an explicit, composed state container defined in
`src/core/swap_state_mod.f90`: the top-level `swap_state_t` derived type
is a record of domain-specific sub-states, each of which lives in its own
module and carries only the variables that a single domain owns.

The composed sub-types are `soil_state_t` (from
`src/soil/soil_state.f90`), `atmosphere_state_t`
(`src/atmosphere/atmosphere_state.f90`), `boundary_state_t`
(`src/boundary/boundary_state.f90`), `drainage_state_t` and
`surfacewater_state_t` (`src/drainage/`), `heat_state_t`
(`src/heat/heat_state.f90`), `solute_state_t` (`src/solute/solute_state.f90`),
`macropore_state_t` (`src/macropore/macropore_state.f90`),
`cropgrowth_state_t` (`src/crop/cropgrowth_state.f90`), plus a
`time_state_t` and an `io_handles_t` that keep time/control data and file
handles separate from physical state. Physics routines now receive
`state` as an argument and access only the sub-record they need.

Because not every physics routine has been converted yet, rescue-era
scaffolding in `src/core/swap_state_sync.f90` provides a bidirectional
bridge between `swap_state_t` and the legacy `variables` module. The
dynamic step begins with `state_to_variables(state)` to restore legacy
globals from the aggregated state and ends with `state_from_variables(state)`
to copy changes back. This bridge is intentionally temporary and is slated
for removal once all routines take `state` directly. For the full pattern
— the `config_t` / `initial_t` / `state_t` separation, the use of
`ASSOCIATE` for aliasing inside physics routines, and the lifecycle
contract between `*_state_init` and `*_state_finalize` — see
[state-management.md](state-management.md).

## I/O layer

**As of Phase 4f-extend (2026-05-05):** TOML is the only runtime input
path. Configuration enters via `src/io/toml/load_swap_config.f90`
(`load_swap_config`), which parses `swap.toml` plus its companion files
(`swap.dra.toml`, `*.crp.toml`, CSV companions) through the vendored
toml-f library, populates a typed `swap_config_t` (`src/config/`), and
hands off to `src/io/toml/config_to_variables.f90`
(`config_to_variables`) which copies the validated config into the
legacy `variables` module globals that the physics solvers read.

Per-rotation crop configuration loads through three init modules —
`src/crop/cropfixed_init.f90`, `src/crop/cropwofost_init.f90`,
`src/crop/cropgrass_init.f90` — that read from the cached
`crop_config_global` (populated by `read_crop_toml.f90`) at the moment
each rotation activates.

Meteo data is CSV: a multi-year `.csv` file declared via
`[meteorology.temporal].file`, optionally accompanied by a sub-daily
`detail_file` and a rain-events `events_file`. The reader is
`src/io/readmeteo.f90` (`ReadMeteoYear`, `ReadRainEvents`,
`MeteoCSVYear`, `MeteoCSVDetYear`); all TTutil-based per-year `.YYY`
and all-years `.met` reader paths were retired in SS-5 (ADR 0014).

The legacy fixed-format reader `src/io/readswap.f90` and its per-crop
helpers (`readcropfixed`, `readwofost`, `readgrass`, `rddre`,
`readarablelandgerm`) were **physically deleted** on 2026-05-06 as the
follow-on to ADR 0019. The case-1 init blocks they fed
(`irrigation(1)`, `SoilManagement(1)`, the cropgrowth nutrient block,
`ArableLandGerm(1)`) were collapsed to no-op returns or removed
outright. Parity test suites under `tests/unit/io/toml/` switched to
literal-value assertions and no longer drive any legacy reader. See
ADR 0019's "Update 2026-05-06" section for the deletion record.

Two TTutil utility calls survive in the production runtime — they are
not data readers:
- `swap_main.f90` — `rdsets` / `rdfrom` for the optional `reruns.dat`
  parameter sweep mechanism (auto-skipped if file missing).
- `swapoutput.f90` — `rddtmp` for TTutil scratch-file cleanup at exit.

Two stub-readers also survive — `Read_Tillage` (`tillage.f90:433`) and
`SSDI_irrigation(1)` (`irrigation.f90:580`) both call
`RDinit(unit, 0, swpfile)` to look up `swtill` / `swssdi`. Both default
to 0 in every regression case; the swpfile-pointing hack at
`config_to_variables.f90:1171-1186` keeps the file open succeeding.
Tracked as the SS-10.5 follow-up.

Output is CSV-first. `src/io/swap_csv_output.f90` provides the shared CSV
writer primitives used by the domain-specific output modules:
`src/io/swapoutput.f90` for the soil-water, irrigation, temperature,
solute, snow, crop, and surface-water writers, and
`src/io/macroporeoutput.f90` for macropore-specific output. Per-domain
writers are called from task 2 on days where `flOutput`, `flOutputShort`,
`flHarvestDay`, or `flIrrigationOutput` indicates the user requested
output. For the full TOML schema — every section and key that
`ReadSwapToml_state` and `ReadDrainageToml_state` recognize — see
[configuration-schema.md](configuration-schema.md).

## Subdirectory map

- `src/atmosphere/` — Atmospheric forcing and canopy boundary processes:
  precipitation partitioning, canopy interception, reference and actual
  evapotranspiration, snow accumulation and melt, and day- and sub-day
  meteo processing. Files: `atmosphere_constants.f90`,
  `atmosphere_state.f90`, `et.f90`, `interception.f90`, `meteoday.f90`,
  `meteodt.f90`, `precipitation.f90`, `snow.f90`. Exports
  `atmosphere_state_t`.

- `src/boundary/` — Top and bottom boundary conditions for the soil
  column: `boundtop.f90` (surface fluxes consistent with ponding and
  runoff logic) and `boundbottom.f90` (groundwater/pressure-head or
  free-drainage options). Exports `boundary_state_t` alongside
  constants in `boundary_constants.f90`.

- `src/core/` — Entry points and cross-cutting infrastructure. The
  driver `swap_main.f90`, the three-phase entry point `swap.f90`, the
  time controller `timecontrol.f90`, one-shot initializer
  `initialize.f90`, the aggregated state module `swap_state_mod.f90`,
  the rescue-era legacy bridge pair `variables.f90` +
  `swap_state_sync.f90`, structured logging in `swap_log.f90`, global
  `constants.f90`, and array-dimension parameters in `arrays.f90` plus
  the `arrays.fi` / `params.fi` / `description.fi` include files. This
  subdirectory is the only one that `use`s every domain's state module.

- `src/crop/` — Crop growth and crop–soil interaction. Three crop
  models coexist: fixed (LAI/crop-factor tables), grass, and WOFOST.
  Files include `cropgrowth.f90`, `cropgrowth_state.f90`,
  `irrigation.f90`, `management_soil.f90`, `oxygenstress.f90`,
  `rootextraction.f90`, `tillage.f90`, the WOFOST nutrient and
  soil-interaction files `wofostnut.f90` and `wofost_soil_*.f90`. This
  is the largest domain by file count and accounts for most of the
  model's optional complexity. Exports `cropgrowth_state_t`.

- `src/drainage/` — Drainage fluxes to ditches and tile drains
  (`drainage.f90`), the `divdra.f90` routine that distributes total
  drainage over drainage levels, and the polder-style surface-water
  reservoir in `surfacewater.f90`. Exports `drainage_state_t` and
  `surfacewater_state_t`.

- `src/error/` — Empty at the rescue baseline. Placeholder for the
  Phase 4 error-handling module that will replace the current
  `fatalerr` / ad-hoc `call exit` pattern.

- `src/heat/` — Soil temperature transport (`temperature.f90`) and the
  frozen-soil thermal-conductivity reduction (`frozencond.f90`).
  Exports `heat_state_t`.

- `src/io/` — All file I/O. Canonical TOML readers `readswaptoml.f90`
  and `readdrainagetoml.f90`; the meteo reader `readmeteo.f90`; the
  CSV output primitives in `swap_csv_output.f90`; the domain output
  writers in `swapoutput.f90` and `macroporeoutput.f90`; and a small
  `checkdate.f90` (date-range validator extracted from the deleted
  legacy reader, still used by `read_ssdi_input` for the swssdi=1
  path pending ADR 0021).

- `src/macropore/` — Preferential flow through the macropore continuum:
  rate calculations (`macrorate.f90`), the integrator (`macropore.f90`),
  and the state carrier (`macropore_state.f90`). Exports
  `macropore_state_t`.

- `src/soil/` — Soil hydraulics and the soil-water solver. Van
  Genuchten / Mualem / PDI constitutive models in
  `WC_K_models_04_11.f90`, tabulated-property support in
  `sptabulated.f90`, grid construction and discretization in
  `soilgrid.f90`, the main Richards solver wrapping in
  `soilhydraulics.f90`, and the water-balance integrator
  `waterbalance.f90`. Exports `soil_state_t`.

- `src/solute/` — Solute transport: advection, hydrodynamic dispersion,
  and simple first-order / linear-sorption reactions. Files:
  `solute.f90`, `solute_state.f90`. Exports `solute_state_t`.

- `src/utils/` — Cross-cutting low-level helpers that multiple domains
  depend on: array utilities (`arrayutils.f90`), I/O utilities
  (`ioutils.f90`), numerical solvers (`numericalsolvers.f90`, including
  the Thomas tridiagonal algorithm `tridag` and the band-matrix routines
  `bandec` / `banbks`), soil-hydraulics helpers
  (`soilhydraulicsutils.f90`), surface-water helpers
  (`surfacewaterutils.f90`), and the shared simulation / exchange
  helpers (`sharedsimulation.f90`, `sharedexchange.f90`).

## Dependency direction

The module graph is intended to flow inward toward the core: physics
subdirectories depend on `src/core/` (for the aggregated state types,
time control, and constants) and on `src/utils/` (for numerical solvers
and low-level helpers). Physics subdirectories do not depend on each
other directly. When a module in one domain needs a value computed by
a module in another domain — soil water needing root-zone extraction,
say, or drainage needing the current surface-water level — that value
lives in the corresponding sub-state of `swap_state_t`, and the
consumer reads it from there rather than `use`-importing the producer.

The consequence is that `src/core/swap_state_mod.f90` is the only module
in the codebase that imports every domain's `*_state_t` type. Domain
modules import only their own state type plus whatever they inherit from
`src/core/` and `src/utils/`. This keeps the dependency graph flat,
makes inter-domain data flow visible at the `state%` accessor level
instead of buried in a global module, and leaves the door open for
multi-instance execution and parallel runs where each instance owns its
own `swap_state_t`. Residual coupling through `variables` exists only
for the duration of the rescue and is isolated to `swap_state_sync.f90`.

## Rescue caveat

This repository is currently under a rescue-and-stabilize workflow, and
parts of the architecture described above are still in transition. In
particular the legacy `variables` module coexists with `swap_state_t`,
and follow-on rescue work continues. The legacy fixed-format reader
code that previously lived in `src/io/readswap.f90` and the case-1
init blocks it fed were physically deleted on 2026-05-06 (see ADR
0019's closing update); the production runtime has been TOML-only
since 2026-05-05. For the full modernization history — capstone
summaries and ADRs — see `dev-docs/` in the repository
(`phase-4-modernization-summary.md`,
`post-phase-4-modernization-summary.md`, `adr/`).
