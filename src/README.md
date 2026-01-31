# SWAP Source Code Structure

This directory contains the SWAP model source code organized into thematic/functional modules.

## Module Overview

```
src/
├── core/          # Core framework, state management, main program
├── io/            # Input file parsing and output generation
├── soil/          # Soil physics, water flow, numerical solution
├── atmosphere/    # Meteorology, evapotranspiration, snow
├── crop/          # Crop growth, nutrients, irrigation, management
├── drainage/      # Lateral drainage and surface water
├── boundary/      # Top and bottom boundary conditions
├── macropore/     # Preferential flow through macropores
├── solute/        # Solute transport
├── heat/          # Soil temperature and frost
└── utils/         # Shared simulation utilities
```

## Module Descriptions

| Module | Description | Key Files |
|--------|-------------|-----------|
| **core** | Central framework providing state variables, initialization, time control, and the main SWAP subroutine | `variables.f90`, `swap.f90`, `swap_main.f90` |
| **io** | All input/output operations including reading the main input file and writing results | `readswap.f90`, `swapoutput.f90`, `swap_csv_output.f90` |
| **soil** | Soil water flow simulation using Richards equation with Mualem-van Genuchten hydraulics | `headcalc.f90`, `functions.f90`, `soilwater.f90` |
| **atmosphere** | Atmospheric boundary: Penman-Monteith ET, interception, snow dynamics | `penmon.f90`, `meteoday.f90`, `snow.f90` |
| **crop** | Crop growth models (simple to WOFOST), root water uptake, irrigation, tillage | `cropgrowth.f90`, `rootextraction.f90`, `irrigation.f90` |
| **drainage** | Lateral drainage (Hooghoudt/Ernst equations) and surface water management | `drainage.f90`, `surfacewater.f90` |
| **boundary** | Top (infiltration/evaporation) and bottom (groundwater) boundary conditions | `boundtop.f90`, `boundbottom.f90` |
| **macropore** | Dual-domain preferential flow through shrinkage cracks and biopores | `macropore.f90`, `macrorate.f90` |
| **solute** | Advection-dispersion equation for reactive solute transport | `solute.f90` |
| **heat** | Soil heat flow and frozen soil conditions | `temperature.f90`, `frozencond.f90` |
| **utils** | Utilities for coupling with external models | `sharedexchange.f90`, `sharedsimulation.f90` |

## Dependencies

The modules have the following dependency structure:

```
                      ┌─────────┐
                      │  core   │
                      └────┬────┘
           ┌───────────────┼───────────────┐
           ▼               ▼               ▼
      ┌────────┐      ┌────────┐      ┌────────┐
      │   io   │      │  soil  │      │  crop  │
      └────────┘      └───┬────┘      └───┬────┘
                          │               │
        ┌─────────────────┼───────────────┤
        ▼                 ▼               ▼
   ┌─────────┐      ┌──────────┐    ┌───────────┐
   │ boundary│      │ drainage │    │ macropore │
   └─────────┘      └──────────┘    └───────────┘
        │                 │               │
        └─────────────────┼───────────────┘
                          ▼
              ┌───────────────────────┐
              │ atmosphere │  heat    │
              │   solute   │  utils   │
              └───────────────────────┘
```

All modules depend on `core/variables.f90` for shared state.

## Building

From the repository root:

```bash
# Using pixi (recommended)
pixi run build-linux

# Manual build
meson setup builddir
meson compile -C builddir
```

## Testing

Unit tests are organized in `tests/unit/` with a subfolder for each module.

```bash
# Run all tests
pixi run test-linux

# Run specific test case
pixi run test-linux-hupselbrook
```

## Adding New Files

1. Create the file in the appropriate module folder
2. Add the file path to `meson.build` in the `sources` list
3. Consider the compilation order (dependencies must come first)
4. Update the module's README.md

## Include Files

Common include files are located in `core/`:
- `arrays.fi` - Array dimension parameters
- `params.fi` - Numerical constants
- `description.fi` - Model version information
