# Atmosphere Module

This folder contains routines for atmospheric processes including evapotranspiration calculation, meteorological data processing, and snow dynamics.

## Files

| File | Description |
|------|-------------|
| `penmon.f90` | Penman-Monteith evapotranspiration calculations. Computes potential evaporation from bare soil (ES0), potential transpiration from dry canopy (ET0), and potential evaporation from wet canopy (EW0). Includes detailed radiation balance calculations and aerodynamic resistance computations. |
| `meteoday.f90` | Daily meteorological data processing. Contains the `MeteoVars` module and routines for processing daily weather, calculating interception (Von Hoyningen-Hüne/Braden, Gash, Rutter methods), and reducing soil evaporation. |
| `meteodt.f90` | Sub-daily (time step) meteorological processing. Handles detailed rainfall events, distributes ET as sine wave during the day, and processes weather data at iteration time scale. |
| `snow.f90` | Snow accumulation and melt simulation. Models snow pack formation, temperature-driven melting, rain-on-snow events, and liquid water storage in snow. |

## Key Processes

### Penman-Monteith Equation
Reference evapotranspiration is calculated using the FAO Penman-Monteith approach:

$$ET_0 = \frac{0.408 \Delta (R_n - G) + \gamma \frac{900}{T+273} u_2 (e_s - e_a)}{\Delta + \gamma (1 + 0.34 u_2)}$$

### Interception Models

Three interception methods are available:
1. **Von Hoyningen-Hüne/Braden** - Simple empirical model for agricultural crops
2. **Gash** - Analytical model based on rainfall characteristics
3. **Rutter** - Detailed canopy water balance model for forests

### Snow Melt
Temperature-index approach with corrections for:
- Degree-day melting
- Rain-on-snow contributions
- Sublimation losses

## Dependencies

- `core/variables.f90` - Meteorological state variables
- Weather input files (`.met` format)

## Refactoring
### Modularization
- [ ] 