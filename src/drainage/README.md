# Drainage Module

This folder contains routines for lateral drainage, surface water management, and the interaction between groundwater and drainage systems.

## Files

| File | Description |
|------|-------------|
| `drainage.f90` | Main lateral drainage calculations. Implements Hooghoudt and Ernst drainage equations for pipe drains and open ditches. Supports up to 5 drainage levels with different characteristics. Calculates drainage/infiltration fluxes based on groundwater level and drainage resistance. |
| `divdra.f90` | Distribution of drainage fluxes over soil compartments. Calculates model discharge layers and distributes lateral fluxes vertically through the saturated zone based on transmissivity. Handles both drainage outflow and infiltration from surface water. |
| `surfacewater.f90` | Surface water management and extended drainage. Simulates water levels in primary and secondary surface water systems, weir management, pump operation, and surface water-groundwater interaction. Handles runoff and surface ponding. |

## Drainage Equations

### Hooghoudt Equation
For tile drains with equivalent depth correction:
$$q_d = \frac{8 K_b d_e \Delta h + 4 K_a (\Delta h)^2}{L^2}$$

Where:
- $q_d$ = drainage flux (m/d)
- $K_a$, $K_b$ = hydraulic conductivity above/below drain
- $d_e$ = equivalent depth (Hooghoudt correction)
- $\Delta h$ = midpoint water table height above drain
- $L$ = drain spacing (m)

### Ernst Equation
For layered soils with radial and horizontal flow components:
$$q_d = \frac{\Delta h}{W_h + W_r + W_e}$$

Where $W_h$, $W_r$, $W_e$ are horizontal, radial, and entry resistances.

## Drainage Levels

SWAP supports up to 5 drainage levels:
1. Primary system (main ditches, rivers)
2. Secondary system (field ditches)
3. Tertiary system (tile drains)
4. Quaternary system (additional drains)
5. Rapid macropore drainage

Each level can have:
- Different depths and spacings
- Individual drainage/infiltration resistances
- Time-variable surface water levels

## Surface Water Options

| Mode | Description |
|------|-------------|
| Basic | Prescribed water levels in drainage system |
| Extended | Dynamic simulation of surface water levels |
| Multi-level | Different water levels per drainage system |

## Dependencies

- `core/variables.f90` - Drainage state variables
- `soil/calcgwl.f90` - Groundwater level calculation
- `macropore/` - Rapid drainage through macropores
