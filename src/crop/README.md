# Crop Module

This folder contains routines for crop growth simulation, root water uptake, nutrient dynamics, and agricultural management.

## Files

### Crop Growth

| File | Description |
|------|-------------|
| `cropgrowth.f90` | Main crop growth simulation (~4500 lines). Implements multiple crop models: simple crop, detailed crop with fixed development (DVS), WOFOST-based dynamic crop growth. Handles crop calendar, emergence, growth stages, LAI development, root growth, and harvest. |
| `rootextraction.f90` | Root water uptake calculations. Implements multiple reduction functions: Feddes et al. (1978) for drought/wet stress, De Jong van Lier (2008) microscopic approach, compensation mechanisms (Jarvis 1989). |
| `oxygenstress.f90` | Oxygen stress model for waterlogged conditions. Implements Bartholomeus et al. (2008) approach calculating oxygen diffusion, respiration, and resulting stress factors (~1800 lines). |

### WOFOST Soil Nutrient Model

| File | Description |
|------|-------------|
| `wofostnut.f90` | Crop nutrient uptake and stress. Calculates nitrogen nutrition index (NNI) and nutrient-limited growth. Based on LINTUL4 and WOFOST nutrient modules. |
| `wofost_soil_declarations.f90` | Variable declarations for the WOFOST soil module including organic matter pools, mineral nitrogen, and soil amendment parameters. |
| `wofost_soil_interface.f90` | Interface routines connecting SWAP soil water with WOFOST soil nutrient dynamics. |
| `wofost_soil_parameters.f90` | Parameter definitions for organic matter decomposition, nitrification, denitrification rate constants. |
| `wofost_soil_rateconstants.f90` | Temperature and moisture response functions for soil biological processes. |
| `wofost_soil_orgmatn.f90` | Organic matter and nitrogen dynamics. Models fresh organic matter (FOM), biomass (BIO), humus (HUM) pools and N mineralization. |
| `wofost_soil_watern.f90` | Soil water and nitrogen interaction. Calculates N uptake based on water uptake and concentration. |
| `wofost_soil_amendments.f90` | Handles organic amendments (manure, crop residues) and their decomposition. |
| `wofost_soil_cropresidues.f90` | Crop residue incorporation and decomposition after harvest. |
| `wofost_soil_balancecheck.f90` | Mass balance verification for organic matter and nitrogen pools. |

### Management

| File | Description |
|------|-------------|
| `irrigation.f90` | Irrigation scheduling and application. Supports fixed scheduling, threshold-based triggering (soil moisture, pressure head), and deficit-based irrigation with various criteria. |
| `tillage.f90` | Soil tillage simulation. Models bulk density changes, hydraulic property modifications, and soil structure evolution after tillage events. |
| `management_soil.f90` | Soil management operations including organic matter applications, fertilizer additions, and coordination with WOFOST soil module. |

## Root Water Uptake

### Feddes Reduction Function
Uptake reduction based on pressure head thresholds:
- h1 (anaerobiosis point)
- h2 (field capacity, optimal range start)
- h3h/h3l (wilting point, high/low transpiration)
- h4 (permanent wilting point)

### Salinity and Frost Stress
Additional reduction factors for:
- Osmotic stress from high EC
- Low temperature effects on root function

## Dependencies

- `core/variables.f90` - Crop state variables
- `soil/` - Soil water state for stress calculations
- `atmosphere/` - Meteorological data for growth calculations
- Crop parameter files (`.crp` format)
