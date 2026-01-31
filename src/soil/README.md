# Soil Module

This folder contains routines for soil physics, water flow calculations, and the numerical solution of the Richards equation.

## Files

### Hydraulic Properties

| File | Description |
|------|-------------|
| `functions.f90` | Core hydraulic functions: water content (`watcon`), hydraulic conductivity (`hconduc`), moisture capacity (`moiscap`), mean conductivity calculations (arithmetic, geometric, harmonic means). |
| `WC_K_models_04_11.f90` | Advanced water retention and conductivity models including Mualem-van Genuchten (unimodal and bimodal), PDI model, and film flow corrections. |
| `sptabulated.f90` | Spline-based tabulated soil physical relations. Contains TSPACK routines for smooth interpolation of measured retention curves. Also includes the `doln` module for log-normal transformations. |
| `hysteresis.f90` | Hysteresis in soil water retention. Tracks wetting/drying history and adjusts van Genuchten parameters for scanning curves. |

### Numerical Solution

| File | Description |
|------|-------------|
| `headcalc.f90` | Main iterative solver for pressure head calculations. Implements the implicit finite difference solution of Richards equation with Newton-Raphson iteration. |
| `tridag.f90` | Tridiagonal matrix solver (Thomas algorithm). Solves the linearized system of equations arising from the implicit discretization. |
| `fluxes.f90` | Calculates water fluxes between soil compartments from the solved pressure head distribution. |

### Grid and State

| File | Description |
|------|-------------|
| `soilwater.f90` | Soil water state variable calculations. Initializes soil hydraulic properties and manages soil water storage state. |
| `calcgrid.f90` | Calculates computational grid parameters (node positions, distances, layer assignments) from input discretization. |
| `calcgwl.f90` | Groundwater level calculation. Searches for water table and perched water tables in the soil profile. |
| `watstor.f90` | Updates water storage in the soil profile. Simple integration of volumetric water content over compartments. |
| `convertdiscrvert.f90` | Converts vertical discretization between different grid configurations (for ANIMO/PEARL coupling). |
| `integral.f90` | Calculates intermediate and cumulative water fluxes for mass balance tracking. |
| `checkmassbal.f90` | Mass balance checking routines. Validates water conservation per period for model coupling. |

## Key Equations

### Richards Equation
The unsaturated flow is described by:
$$\frac{\partial \theta}{\partial t} = \frac{\partial}{\partial z}\left[K(h)\left(\frac{\partial h}{\partial z} + 1\right)\right] - S(h)$$

### Mualem-van Genuchten
Water retention:
$$\theta(h) = \theta_r + \frac{\theta_s - \theta_r}{[1 + |\alpha h|^n]^m}$$

Hydraulic conductivity:
$$K(h) = K_s S_e^l [1-(1-S_e^{1/m})^m]^2$$

## Dependencies

- `core/variables.f90` - State variables and parameters
- `core/arrays.fi` - Array dimension parameters
