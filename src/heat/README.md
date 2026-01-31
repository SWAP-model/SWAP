# Heat Module

This folder contains routines for soil temperature simulation and frost dynamics.

## Files

| File | Description |
|------|-------------|
| `temperature.f90` | Main soil temperature simulation. Solves the heat conduction equation with temperature-dependent thermal properties. Supports both analytical (sine wave) and numerical solutions. Accounts for soil composition (minerals, organic matter, water, air) in calculating thermal conductivity and heat capacity. |
| `frozencond.f90` | Frozen soil conditions. Determines reduction factors for hydraulic conductivity in frozen soil, calculates frost depth, and handles freeze-thaw dynamics. Links temperature state to water flow restrictions. |

## Heat Transport Equation

The soil heat flow equation:
$$C_h \frac{\partial T}{\partial t} = \frac{\partial}{\partial z}\left[\lambda \frac{\partial T}{\partial z}\right]$$

Where:
- $T$ = soil temperature (°C)
- $C_h$ = volumetric heat capacity (J/cm³/K)
- $\lambda$ = thermal conductivity (W/m/K)

## Thermal Properties

### Heat Capacity
$$C_h = \sum_i \phi_i C_i$$

Where $\phi_i$ are volume fractions of:
- Mineral particles
- Organic matter
- Water
- Air

### Thermal Conductivity
Uses the De Vries method accounting for particle shape factors and the continuous medium (water or air dominated).

## Solution Methods

### Analytical (SWCALT = 1)
Sine wave approximation assuming harmonic temperature variation:
$$T(z,t) = T_{mean} + A_0 \exp(-z/d) \sin(\omega t - z/d + \phi)$$

Where $d$ is the damping depth.

### Numerical (SWCALT = 2)
Implicit finite difference solution with:
- Upper BC: air temperature or energy balance
- Lower BC: specified temperature or zero gradient

## Frost Effects

### Hydraulic Conductivity Reduction
$$K_{frozen} = K_{unfrozen} \cdot f_{frost}(T)$$

Linear reduction between:
- $T_{frost,start}$ (typically 0°C)
- $T_{frost,end}$ (typically -3°C)

### Frost Depth
Calculated by interpolating between nodes where temperature crosses 0°C.

### Impact on Transpiration
Root water uptake is reduced or eliminated in frozen soil layers.

## Dependencies

- `core/variables.f90` - Temperature state variables
- `soil/` - Soil composition for thermal properties
- `atmosphere/` - Air temperature boundary condition
- `crop/` - Frost effects on root uptake
