---
title: Next refactoring step
author: Mateusz Zawadzki
---

# Current status
I have divided the states now into:
- `config` and `init`, which are variables read from the configuration files
- `state` which are variables that are being processed and passed between subroutines.

Currently the `state`s are under the config and init modules that correspond to the sections of the configuration files.

# The next refactoring step
## New states
The next step in the refactoring should be to create compartment state modules organized by physical storage locations rather than configuration file sections


- ponding_state_t - surface water with runoff tracking
- canopy_state_t - interception
- snow_state_t - frozen storage
- soil_state_t - unsaturated zone (+ heat + solute)
- saturated_zone_state_t - groundwater (+ heat + solute)
- macropore_state_t - preferential flow (to consider)
- crop_state_t - living biomass
- drainage_state_t - artificial drainage
- irrigation_state_t - scheduled irrigation events
- surface_water_state_t - managed surface water (optional)
- boundary_reader_state_t - bottom boundary file I/O (if swbbcfile=1)

Non-storage:

- meteo_forcing_t - current weather (forcing)
- meteo_reader_state_t - file I/O (infrastructure)

Then we will start vectorizing the subroutines and functions to pass the "water" between these states. This will dramatically improve the redability of the code and make the parallelization and gpu acceleration even easier. 

## Refactoring processes

Refactor process functions to explicitly pass water/heat/solute between compartments:

```fortran
    call compute_infiltration(ponding, soil, flux)  ! ponding → soil
    call compute_interception(canopy, forcing, throughfall)  ! gross → canopy → net
    call compute_drainage(soil, drainage, flux)  ! soil → drains
```
Separate computed fluxes from stored state - **fluxes become function outputs, not state variables**

Benefits:

- Clear water balance: each compartment tracks its storage + cumulative in/out
- Easier parallelization: independent compartments can be solved in parallel
- GPU acceleration: array operations on compartment profiles
- Testing: mock individual compartments and test processes in isolation
- Maintainability: each compartment enforces its own physical constraints

Before compartment refactoring, ensure:

- [ ] All modules have clean config_t, initial_t, state_t separation
- [ ] TOML parser correctly populates both config and initial
- [ ] Main program uses init(state, config, initial) pattern
- [ ] All process functions use config as intent(in)
- [ ] Model runs correctly with current sectional organization

Then proceed with compartment refactoring incrementally (one compartment at a time).

# Some further ideas for improvement
- [ ] (MINOR) include a disable module guard in the initialization. Currently everything is always initialized and finalized, but this could be reduced.