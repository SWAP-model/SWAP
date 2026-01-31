# Solute Module

This folder contains routines for simulating solute transport in the soil, with focus on simple conservative and reactive tracers.

## Files

| File | Description |
|------|-------------|
| `solute.f90` | Main solute transport simulation. Implements convection-dispersion equation for a single reactive solute. Includes linear adsorption (Freundlich isotherm), first-order decay, root uptake, and volatilization. Handles boundary conditions at soil surface and bottom. |

## Transport Equation

The convection-dispersion equation (CDE) solved is:

$$\frac{\partial (\theta c + \rho_b s)}{\partial t} = \frac{\partial}{\partial z}\left[D \frac{\partial c}{\partial z} - q c\right] - \mu (\theta c + \rho_b s) - S_r c$$

Where:
- $c$ = solute concentration in soil water (mg/cm³)
- $s$ = adsorbed concentration (mg/g)
- $\theta$ = volumetric water content (-)
- $\rho_b$ = bulk density (g/cm³)
- $D$ = dispersion coefficient (cm²/d)
- $q$ = water flux (cm/d)
- $\mu$ = decay rate constant (1/d)
- $S_r$ = root water uptake rate (1/d)

## Adsorption

Freundlich isotherm:
$$s = K_f \cdot c^{n_f}$$

For $n_f = 1$ (linear isotherm):
$$R = 1 + \frac{\rho_b K_d}{\theta}$$

Where $R$ is the retardation factor.

## Decay/Transformation

First-order decay with temperature and moisture corrections:
$$\mu_{act} = \mu_{ref} \cdot f_T(T) \cdot f_\theta(\theta)$$

## Boundary Conditions

### Top Boundary
- Incoming solute with rain/irrigation
- Volatilization losses
- Surface application (instantaneous or constant flux)

### Bottom Boundary
- Convective outflow
- Zero gradient condition

## Applications

The solute module can simulate:
- Conservative tracers (Cl, Br)
- Age tracers (for residence time)
- Simple pesticides (with decay)
- Salts (for salinity effects on crops)

For detailed pesticide fate, SWAP is typically coupled with PEARL model.

## Dependencies

- `core/variables.f90` - Solute state variables
- `soil/` - Water content and fluxes
- `crop/` - Root uptake of solutes
