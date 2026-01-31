# I/O Module

This folder contains all input/output routines for reading input files and writing output files.

## Files

| File | Description |
|------|-------------|
| `readswap.f90` | Reads the main SWAP input file (`.swp`). Parses all model configuration including soil properties, drainage parameters, boundary conditions, crop settings, and output options. This is the primary input parser (~5000 lines). |
| `readmeteo.f90` | Reads meteorological input data from weather files (`.met`). Handles both daily and detailed (sub-daily) weather data formats. Includes routines for reading rain event data and parameter validation. |
| `swapoutput.f90` | Main output controller for traditional SWAP output files. Manages opening, writing, and closing of various output files (`.bal`, `.blc`, `.vap`, etc.). Contains extensive output formatting routines (~4700 lines). |
| `swap_csv_output.f90` | Modern CSV output module providing flexible, user-configurable output. Supports selection of output variables via input file, handles nodal and profile outputs, and produces machine-readable CSV format. |

## Input File Formats

SWAP uses the FSE (Fortran Simulation Environment) input format with:
- Keyword-based parameter assignment
- Tables for time-dependent data
- Comments marked with `*` or `!`

## Output File Types

| Extension | Content |
|-----------|---------|
| `.csv` | User-defined variables in CSV format |
| `.bal` | Water balance output |
| `.blc` | Detailed water balance components |
| `.vap` | Vertical profile output |
| `.log` | Simulation log and warnings |

## Dependencies

- `core/variables.f90` - Access to all model state variables
- Relies on TTUTIL library for FSE-format parsing (`rdparm`, `rdtable`, etc.)
