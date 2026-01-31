# Utils Module

This folder contains utility routines for shared simulations and model coupling.

## Files

| File | Description |
|------|-------------|
| `sharedexchange.f90` | Data exchange routines for coupled simulations. Contains `FromSwap` and `ToSwap` subroutines for writing SWAP results to and reading data from shared files. Used when SWAP is coupled with other models (WOFOST, ANIMO, PEARL). Currently contains stub implementations. |
| `sharedsimulation.f90` | Shared simulation coordination. Handles process synchronization when SWAP runs alongside other models. Manages shared data files, process IDs, and timing coordination for multi-model simulations. |

## Coupling Concept

SWAP can be coupled with other models in several ways:

### Loose Coupling
- Models run sequentially
- Data exchange via files
- Used for: SWAP → ANIMO, SWAP → PEARL

### Tight Coupling
- Models run simultaneously
- Data exchange via shared memory/files
- Real-time synchronization required
- Used for: SWAP ↔ WOFOST (internal), SWAP ↔ external crop models

## Data Exchange Variables

### From SWAP (output to other models)
- Soil water content profiles
- Root zone water status
- Temperature profiles
- Water fluxes
- Groundwater levels

### To SWAP (input from other models)
- Crop LAI and height
- Root depth
- Potential transpiration
- Crop-specific parameters

## Implementation Notes

These modules are designed as extension points for coupling SWAP with external models. The current implementation contains the interface structure but requires customization for specific coupling applications.

## Dependencies

- `core/variables.f90` - Access to model state
- External model interfaces (application-specific)
