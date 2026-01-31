# Crop Unit Tests

This folder contains unit tests for crop growth and management routines.

## Test Scope

Tests in this folder should cover:

### Crop Growth
- **cropgrowth.f90**
  - Development stage calculations
  - LAI development
  - Root growth patterns
  
- **rootextraction.f90**
  - Feddes reduction function
  - Compensation mechanisms
  - Stress factor calculations
  
- **oxygenstress.f90**
  - Oxygen diffusion
  - Respiration calculations
  - Stress reduction factors

### WOFOST Modules
- **wofostnut.f90**
  - Nutrient uptake
  - NNI calculation
  
- **wofost_soil_*.f90**
  - Organic matter decomposition
  - Nitrogen mineralization
  - Mass balance verification

### Management
- **irrigation.f90**
  - Threshold triggering
  - Application calculations
  
- **tillage.f90**
  - Bulk density changes
  - Hydraulic property updates

## Test Files

| File | Description |
|------|-------------|
| (to be created) | Tests for root uptake reduction |
| (to be created) | Tests for crop development |
| (to be created) | Tests for irrigation scheduling |

## Running Tests

```bash
# From the swap directory
pixi run test-unit
```

## Test Cases

- Known stress scenarios (drought, waterlogging)
- Published crop growth data
- Tillage experiment data
