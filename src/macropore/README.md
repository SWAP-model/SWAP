# Macropore Module

This folder contains routines for simulating preferential flow through macropores, shrinkage cracks, and biopores.

## Files

| File | Description |
|------|-------------|
| `macropore.f90` | Main macropore simulation routine (~1800 lines). Implements the MACRO-concept for preferential flow including macropore geometry, water entry at surface and from matrix, lateral exchange between macropores and matrix, and rapid drainage through the macropore network. |
| `macrorate.f90` | Rate calculations for macropore water flow. Computes water fluxes into and through macropores, absorption into matrix, and drainage to tile drains. Handles the water balance of the macropore domain. |
| `macroporeoutput.f90` | Output routines for macropore-specific variables. Writes detailed macropore state and flux information for analysis and debugging. |

## Macropore Concept

SWAP uses a dual-domain approach:

### Two Domains
1. **Main Bypass Domain** - Large continuous cracks/pores that can transport water rapidly to depth
2. **Internal Catchment Domain** - Smaller pores that primarily infiltrate laterally into the matrix

### Macropore Types
- **Shrinkage cracks** - Develop in clay soils during drying
- **Biopores** - Created by roots and soil fauna
- **Structural pores** - Aggregation-related voids

## Key Processes

### Surface Entry
Water enters macropores at the surface when:
- Rainfall intensity > infiltration capacity
- Ponding occurs
- Fraction depends on macropore area at surface

### Vertical Transport
- Rapid downward flow through air-filled macropores
- Flow velocity depends on macropore geometry
- Can reach tile drains directly (rapid drainage)

### Lateral Exchange
Matrix-macropore exchange by:
- **Absorption** - Capillary uptake from macropore walls
- **Infiltration** - Pressure-driven flow when macropore is water-filled
- Based on sorptivity and hydraulic conductivity

### Rapid Drainage
Direct contribution to tile drains when:
- Macropore water level is above drain
- Connected macropore network exists

## Macropore Geometry

$$V_{mac}(z) = V_{mac,surface} \cdot f(z)$$

Where $f(z)$ describes the depth distribution of macropore volume, typically exponentially decreasing with depth.

## Dependencies

- `core/variables.f90` - Macropore state variables
- `soil/` - Matrix water content for exchange
- `drainage/` - Rapid drainage to tile drains
- `boundary/` - Surface water entry
