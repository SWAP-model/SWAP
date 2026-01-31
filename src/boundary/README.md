# Boundary Module

This folder contains routines for handling the top and bottom boundary conditions of the soil column.

## Files

| File | Description |
|------|-------------|
| `boundtop.f90` | Top boundary condition calculations. Determines infiltration, evaporation, and ponding at the soil surface. Handles soil evaporation flux based on atmospheric demand vs. soil supply, surface runoff, and interaction with macropores at the soil surface. |
| `boundbottom.f90` | Bottom boundary condition calculations. Implements multiple options for the lower boundary of the soil column including prescribed groundwater level, prescribed flux, free drainage, lysimeter, and regional aquifer coupling. |

## Top Boundary Conditions

The soil surface can experience:

### Infiltration
- Rainfall minus interception
- Irrigation applications
- Runon from surrounding areas
- Snowmelt contributions

### Evaporation
Actual soil evaporation is limited by:
1. Atmospheric demand (potential evaporation)
2. Soil water supply (hydraulic conductivity at surface)
3. Reduction functions for dry conditions

### Ponding and Runoff
- Surface storage up to maximum ponding depth
- Runoff when ponding exceeds threshold
- Runoff-curve number method (optional)

## Bottom Boundary Options

| SWBOTB | Description |
|--------|-------------|
| -1 | Prescribed groundwater level |
| 1 | Groundwater level from input table |
| 2 | Prescribed regional bottom flux |
| 3 | Calculate flux from regional aquifer head |
| 4 | Calculate flux from average groundwater level |
| 5 | Pressure head at bottom from soil-air interface |
| 6 | Free drainage (unit hydraulic gradient) |
| 7 | Zero flux (impermeable layer) |
| 8 | Free outflow at soil-air interface |

### Regional Aquifer Coupling
For options 3-4, the bottom flux is calculated as:
$$q_{bot} = \frac{h_{aquifer} - h_{bottom}}{c_{aquifer}}$$

Where $c_{aquifer}$ is the vertical resistance of the aquitard.

## Frozen Soil

Both boundaries interact with frozen soil conditions:
- Reduced infiltration when surface is frozen
- Modified bottom flux for frozen deep layers
- Coordinated with `heat/frozencond.f90`

## Dependencies

- `core/variables.f90` - Boundary state variables
- `soil/functions.f90` - Hydraulic conductivity calculations
- `atmosphere/` - Precipitation and evaporation fluxes
- `macropore/` - Surface water entry into macropores
