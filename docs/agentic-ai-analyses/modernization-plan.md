# SWAP codebase survey and modernization options

## What is here today
- **Scope**: SWAP v4.2.0 Fortran 90/77 mix; executable-oriented; heavy global state in `variables` module and fixed-size arrays from `arrays.fi`; ttutil provides text I/O, logging, unit management, and RD/RW helpers.
- **Top-level driver**: `swap_main.f90` calls `swap()` three times (init/dynamic/close). rerun handling via `reruns.dat`; logging to `reruns.log`; exits process after completion.
- **Core time loop**: `swap.f90` orchestrates initialization, day-level logic, sub-daily time steps, and finalization. Exchange structs `swap_input`/`swap_output` exist for DLL/BMI-style calls. Dynamic run loops until `flrunend` is set in `TimeControl`.

### Core computation (keep as intact as possible)
These routines implement physics/biology; they expect globals from `variables` and work mostly in memory:
- **Soil water and flow**: `soilwater.f90`, `boundbottom.f90`, `boundtop.f90`, `tridag.f90`, `calcgrid.f90`, `calcgwl.f90`, `frozencond.f90`, `macropore.f90`, `macrorate.f90`, `macroporeoutput.f90` (state + reporting), `frozencond.f90`.
- **Atmosphere & surface**: `surfacewater.f90`, `fluxes.f90`, `penmon.f90`, `meteodt.f90`, `meteoday.f90`, `snow.f90`, `temperature.f90`.
- **Crop & management**: `cropgrowth.f90`, `rootextraction.f90`, `tillage.f90`, `management_soil.f90`, `irrigation.f90`, `wofost_*` soil/crop nutrient extensions.
- **Solute & tracers**: `solute.f90`, `age tracer` hooks via `sharedexchange` placeholders, `wofostnut.f90` nutrient coupling.
- **Coupled helpers**: `timecontrol.f90`, `integral.f90`, `functions.f90`, `calcgrid.f90`, `checkmassbal.f90`, `hysteresis.f90`, `oxygenstress.f90`, `divdra.f90`, `drainage.f90`.

### Input, configuration, and orchestration
- **Configuration parsing**: `readswap.f90` (main .swp parsing; sets almost all flags, sizes, file names; builds output schedules), `arrays.fi` (dimensions), `params.fi`/`description.fi` for constants and text labels.
- **Meteorological forcing**: `readmeteo.f90` reads annual/detailed meteo files (and optional rain event files), populates arrays for day/time-step use.
- **Shared simulation hooks**: `sharedsimulation.f90` handles command-line shared runs through temporary files and blocking sleep; `sharedexchange.f90` placeholders for exchanging data.
- **Exchange wrappers**: `swap_exchange` (types), optional DLL/BMI handling inside `swap.f90`; prototypes `swap_bmi.f90`, `swap_bmi_original.f90.bak`; `pyswap_interface.f90` exists for Python linking.

### Output and file interaction
- **Standard outputs**: `swapoutput.f90`, `soilwateroutput.f90`, `swap_csv_output.f90`, `macro*output`, `temperatureoutput`, `soluteoutput`, `age tracer` outputs, `irrigationoutput`, `surfacewateroutput`, `snowoutput`.
- **Run bookkeeping**: `swap_main` writes `Swap.ok`, rerun logs, opens/closes files; `checkmassbal` writes balance diagnostics.
- **ttutil dependence**: ubiquitous `rds*`, `rdf*`, `rd*` routines for parsing, `fopens`, `getun/getun2` for units, `mess*` for logging, `sleep`, `rdinqr` presence checks, `dtdpar/dtardp` date conversions. `rdmodulettutil.f90` exposes variable metadata. ttutil is currently baked in via include files and Fortran 77 sources in `ttutil/src`.

### Structural observations
- **Global state**: `variables.f90` declares hundreds of scalars/arrays; mutation is implicit. This prevents reentrancy and thread safety and ties logic to file-driven initialization.
- **Static allocation**: array bounds fixed by `arrays.fi`; limits parallelism and complicates dynamic sizing.
- **Side-effect heavy**: many routines open/read/write files directly; initialization assumes presence of CLI arguments and working directory layout; sleep-based shared-simulation is brittle.
- **Separation of concerns**: computation routines mostly clean, but tightly coupled to globals and flags set by file readers; exchange layer is minimal and optional.

## Modernization goals (as stated)
- Provide a library callable from Python, avoiding executable-only usage.
- Support memory-only mode (no text I/O): Python supplies configuration and forcing; outputs returned in memory.
- Enable orchestration for parallel/multi-run scenarios and future GPU work.
- Keep core computation logic intact where possible.

## Modernization options
### Option A – Wrapper-first, minimal intrusion
- Treat existing core as a single-instance engine.
- Implement a stable C/Fortran API (or BMI) that:
  - Initializes from provided in-memory config/forcing structures (bypassing `readswap`, `readmeteo` by supplying pre-filled `variables` fields).
  - Accepts per-step meteorology and management via `swap_input` and returns state via `swap_output` arrays (extend types as needed).
  - Redirects output routines to callback writers or in-memory buffers; gate file writers behind flags.
- Add a thin "state reset" routine that zeroes globals to defaults between runs.
- Use f2py/ctypes/cffi to expose the API to Python; manage a Python-side lock to avoid concurrent use.
- Pros: fast to deliver; keeps physics untouched. Cons: still single-instance, globals remain, file parsing code largely unused but present.

### Option B – Modular state refactor (recommended phased path)
- Introduce a `type(swap_state)` encapsulating what is currently in `variables` (begin with time/meta/meteo; later crop/soil). Pass it explicitly into routines instead of using module globals. Start at boundaries: `TimeControl`, `ReadMeteoDay`, `SoilWater`, `CropGrowth`.
- Create an `io_bridge` layer that maps between external config structures (from Python JSON/dicts) and the Fortran state, replacing ttutil readers. Maintain a compatibility path that still uses `readswap` for legacy runs.
- Replace output procedures with strategy-style callbacks: interface `SwapOutputSink` that can be a file writer or memory collector. Keep CSV writer behind a flag.
- Make array sizing dynamic: parameterize sizes from input and allocate at runtime; retire `arrays.fi` where possible.
- Separate orchestration: move rerun loop and CLI handling out of library; keep `swap_main` as a thin legacy wrapper calling the new API.
- Prepare for reentrancy: no `save` variables in computation unless stored in the state object; avoid `stop/exit`; return error codes.
- Pros: positions for parallel runs and GPU experimentation; clearer API boundaries. Cons: more effort; requires careful staged migration.

### Option C – Full interface rewrite with new driver
- Write a new driver module that mirrors BMI/EML-style APIs (initialize, update_until, finalize) using the refactored state, leaving legacy drivers untouched. Defer deeper physics changes; focus on clean interfaces and testable boundaries.
- Use this to validate memory-only workflows and incremental replacement of ttutil.

## Suggested phased plan
1) **Inventory & tests**: catalog current I/O touch points and create regression tests for a small case using the existing executable and a library-call harness (Python calling current `swap` via `swap_exchange`).
2) **API skeleton**: finalize Fortran module exposing `swap_init(state, config)`, `swap_step(state, forcing, outputs)`, `swap_finalize(state)`. Add an in-memory logging/error collector.
3) **Decouple inputs**: implement a new config loader (Fortran or Python-side) that fills `swap_state` without ttutil. Keep `readswap` for legacy CLI runs behind feature flags.
4) **Redirect outputs**: introduce callback hooks; wrap legacy writers with guards so library runs can disable file I/O.
5) **State extraction**: move time/meteo/crop/soil arrays from `variables` into the state type progressively; update core routines to accept the state (initially via `use`, then via explicit arguments).
6) **Memory-only path**: make meteo ingestion callable each step from Python; bypass `ReadMeteoYear/Day` file reads when external forcing is provided.
7) **Thread-safety & parallelism**: ensure no global `save` state remains; add a test matrix for concurrent runs (multiprocessing first, then threading if safe).
8) **GPU/acceleration hooks**: identify hotspots (`SoilWater`, `Penman-Monteith`, `macropore`) and prepare data-layout changes (contiguous arrays, SoA) for potential offloading; keep this deferred until API stabilizes.

## Module impact summary
- **Mostly stable (compute)**: soil water/heat/solute/crop modules listed above; changes limited to argument lists and state access.
- **Refactor required**: `swap.f90` (driver sequencing), `swap_main.f90` (CLI-only path), `variables.f90` (state encapsulation), `readswap.f90`/`readmeteo.f90` (new input sources), `swapoutput.f90` and related writers (callback plumbing), `sharedsimulation.f90`/`sharedexchange.f90` (replace with structured interfaces), `swap_bmi.f90` (align with new API), `pyswap_interface.f90` (Python entry points).
- **External dependency**: ttutil can remain for legacy text parsing but should be isolated behind a thin adapter; long term replace with native parsing or Python-side provisioning.

## Next steps (practical)
- Decide between Option A (fast wrapper) and Option B (staged refactor). A can be a stepping stone while building tests for B.
- Define the minimal `swap_state` and `swap_forcing` derived types and map them to existing globals.
- Prototype a Python-driven run that sets forcing via `swap_input` and reads `swap_output` to validate the orchestration loop without file I/O.
- Introduce feature flags to disable file writing and to skip ttutil initialization when running in library mode.
