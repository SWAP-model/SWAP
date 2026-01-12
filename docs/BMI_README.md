# SWAP BMI (Basic Model Interface)

This module provides a standardized interface for running SWAP as a library from Python or other languages.

## Overview

The BMI (Basic Model Interface) is a standard developed by the Community Surface Dynamics Modeling System (CSDMS) for model interoperability. This implementation enables:

- **Zero file I/O overhead**: Set parameters directly in memory
- **Batch processing**: Run hundreds of instances efficiently
- **Parameter estimation**: Easy integration with optimization libraries
- **Model coupling**: Standard interface for linking with other models
- **Future GPU acceleration**: Foundation for CUDA implementation

## Architecture

### Fortran Module (`swap_bmi.f90`)

The BMI module exposes SWAP functionality through C-compatible functions using `iso_c_binding`. Key features:

- **Two initialization modes**:
  - `bmi_initialize(config_file)`: Traditional file-based initialization
  - `bmi_initialize_memory()`: Skip file I/O for batch processing
  
- **Zero-copy data access**: Direct access to Fortran arrays via pointers
- **Day-by-day control**: `update_day()` for stepping through time
- **Parameter setters**: Direct modification of hydraulic properties and meteorology

### Python Wrapper (`pyswap/bmi.py`)

Pythonic interface using ctypes:

```python
from pyswap.bmi import SwapBMI
import numpy as np

with SwapBMI() as swap:
    # Initialize without files
    swap.initialize_memory()
    
    # Set parameters
    swap.set_ksat(np.full(30, 10.0))  # cm/day
    swap.set_alpha(np.full(30, 0.02))  # 1/cm
    swap.set_n(np.full(30, 1.5))
    
    # Run simulation
    for day in range(365):
        swap.set_meteo_day(
            rainfall=5.0,
            etref=4.0,
            tmin=15.0,
            tmax=25.0
        )
        swap.update_day()
        
        # Get results (zero-copy numpy arrays)
        theta = swap.get_theta()
        gwl = swap.get_gwl()
```

## Compilation

### Linux with Intel Fortran

```bash
cd /path/to/swap
meson setup build --buildtype=release
meson compile -C build
```

### Linux with GCC

```bash
cd /path/to/swap
FC=gfortran meson setup build --buildtype=release
meson compile -C build
```

### Windows Cross-Compilation (MinGW)

```bash
cd /path/to/swap
meson setup builddir_windows --cross-file=cross_mingw.txt
meson compile -C builddir_windows
```

After compilation, the shared library will be in:
- Linux: `build/libswap.so`
- Windows: `builddir_windows/swap.dll`

## Installation

### Install Python Package

```bash
cd /path/to/pySWAP
pip install -e .
```

### Install SWAP Library (System-wide)

```bash
cd /path/to/swap
meson install -C build
```

Or manually copy the shared library to a directory in your library path:

```bash
# Linux
sudo cp build/libswap.so /usr/local/lib/
sudo ldconfig

# Or set LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/path/to/swap/build:$LD_LIBRARY_PATH
```

## BMI Functions

### Lifecycle

- `initialize(config_file)`: Initialize from .swp file
- `initialize_memory()`: Initialize without files (for batch mode)
- `update_day()`: Run one day
- `finalize()`: Clean up and close

### Time Management

- `get_current_time()`: Days since 1900
- `get_time_step()`: Time step in days
- `get_time_units()`: Returns "days since 1900"

### Grid Information

- `get_grid_size()`: Number of soil nodes
- `get_var_grid(var_name)`: Grid ID for variable

### State Variables (Getters)

Zero-copy access to arrays:
- `get_theta()`: Soil moisture (m³/m³)
- `get_h()`: Pressure head (cm)
- `get_dz()`: Layer thickness (cm)
- `get_z()`: Depth of nodes (cm)
- `get_rwu()`: Root water uptake per node (cm/day)

Scalar outputs:
- `get_gwl()`: Groundwater level (cm below surface)
- `get_tact()`: Actual transpiration (cm/day)
- `get_tpot()`: Potential transpiration (cm/day)

### Parameters (Setters)

Hydraulic properties (per node):
- `set_ksat(array)`: Saturated hydraulic conductivity (cm/day)
- `set_alpha(array)`: Van Genuchten alpha (1/cm)
- `set_n(array)`: Van Genuchten n (-)
- `set_theta_res(array)`: Residual water content (m³/m³)
- `set_theta_sat(array)`: Saturated water content (m³/m³)

Meteorology (daily):
- `set_rainfall(value)`: Daily rainfall (mm)
- `set_etref(value)`: Reference ET (cm/day)
- `set_tmin(value)`: Minimum temperature (°C)
- `set_tmax(value)`: Maximum temperature (°C)
- `set_meteo_day(rain, et, tmin, tmax)`: Set all at once

State initialization:
- `set_theta(array)`: Initial soil moisture profile

## Usage Examples

See [examples/bmi_examples.py](../../pySWAP/examples/bmi_examples.py) for comprehensive examples:

1. **Basic usage**: File-based initialization
2. **Memory mode**: No file I/O, all parameters in memory
3. **Parameter study**: Vary Ksat and observe response
4. **Batch processing**: Monte Carlo simulation (50+ instances)
5. **Calibration**: Simple parameter estimation

### Quick Start

```python
from pyswap.bmi import SwapBMI
import numpy as np

# Example: Parameter sensitivity
ksat_values = [5, 10, 20, 50]  # cm/day
results = {}

for ksat in ksat_values:
    with SwapBMI() as swap:
        swap.initialize_memory()
        
        # Set soil properties
        n_nodes = swap.get_grid_size()
        swap.set_ksat(np.full(n_nodes, ksat))
        swap.set_alpha(np.full(n_nodes, 0.02))
        swap.set_n(np.full(n_nodes, 1.5))
        swap.set_theta(np.full(n_nodes, 0.3))
        
        # Run 30 days
        gwl_series = []
        for day in range(30):
            swap.set_meteo_day(10.0, 4.0, 15.0, 25.0)
            swap.update_day()
            gwl_series.append(swap.get_gwl())
        
        results[ksat] = gwl_series

# Analyze results...
```

## Performance Tips

### Memory Mode vs File Mode

For batch processing, use `initialize_memory()`:

```python
# SLOW (file I/O for each instance)
for params in parameter_set:
    with SwapBMI() as swap:
        swap.initialize("config.swp")  # Reads files
        # ...

# FAST (no file I/O)
for params in parameter_set:
    with SwapBMI() as swap:
        swap.initialize_memory()  # Skip files
        swap.set_ksat(params['ksat'])
        # ...
```

### Zero-Copy Arrays

Arrays are returned as numpy views of Fortran memory (no copying):

```python
theta = swap.get_theta()  # Zero-copy view
mean_theta = theta.mean()  # Operations on view

# If you need to store for later, explicitly copy:
theta_stored = swap.get_theta().copy()
```

### Parallel Processing

Each `SwapBMI` instance is independent and can be run in parallel:

```python
from multiprocessing import Pool

def run_simulation(params):
    with SwapBMI() as swap:
        swap.initialize_memory()
        swap.set_ksat(params['ksat'])
        # ... run simulation
        return results

with Pool(8) as pool:
    results = pool.map(run_simulation, parameter_sets)
```

## Integration with Other Tools

### xmipy (Advanced)

For more advanced BMI features, use `xmipy`:

```bash
pip install xmipy
```

```python
from xmipy import XmiWrapper

swap = XmiWrapper(lib_path="libswap.so")
swap.initialize("config.swp")
# ... use BMI functions
```

### Model Coupling

BMI enables coupling SWAP with other models:

```python
# Pseudo-code for SWAP-MODFLOW coupling
swap = SwapBMI()
modflow = ModflowBMI()

while time < end_time:
    # Exchange data
    recharge = swap.get_gwl()  # or appropriate variable
    modflow.set_recharge(recharge)
    
    # Update both models
    swap.update_day()
    modflow.update_stress_period()
    
    # Get feedback
    head = modflow.get_head()
    swap.set_bottom_boundary(head)
```

## Troubleshooting

### Library Not Found

```python
# Error: libswap.so not found
```

Solutions:
1. Set `LD_LIBRARY_PATH`:
   ```bash
   export LD_LIBRARY_PATH=/path/to/swap/build:$LD_LIBRARY_PATH
   ```

2. Install system-wide:
   ```bash
   sudo meson install -C build
   ```

3. Specify path in Python:
   ```python
   swap = SwapBMI(lib_path="/absolute/path/to/libswap.so")
   ```

### Symbol Not Found

```python
# Error: undefined symbol: bmi_initialize
```

Check that `swap_bmi.f90` was compiled:
```bash
nm -D build/libswap.so | grep bmi_
```

Should show all BMI functions.

### Segmentation Fault

Common causes:
- Array size mismatch (use `get_grid_size()` first)
- Calling functions before `initialize()`
- Incorrect parameter types

Enable debugging:
```bash
meson setup build --buildtype=debug -Db_sanitize=address
```

## Extending the BMI

To add new variables:

1. **In `swap_bmi.f90`**:
   ```fortran
   function bmi_get_value_ptr_myvar(ptr) result(status) bind(C)
       type(c_ptr), intent(out) :: ptr
       integer(c_int) :: status
       
       ptr = c_loc(myvar)  ! myvar from variables module
       status = BMI_SUCCESS
   end function
   ```

2. **In `pyswap/bmi.py`**:
   ```python
   def get_myvar(self) -> np.ndarray:
       ptr = ctypes.c_void_p()
       self.lib.bmi_get_value_ptr_myvar(ctypes.byref(ptr))
       n = self.get_grid_size()
       return np.ctypeslib.as_array(
           (ctypes.c_double * n).from_address(ptr.value)
       )
   ```

## References

- [CSDMS BMI Documentation](https://csdms.colorado.edu/wiki/BMI)
- [BMI Specification](https://bmi.readthedocs.io/)
- [xmipy Documentation](https://github.com/Deltares/xmipy)
- [SWAP Documentation](https://swap.nl/)

## License

Same as SWAP: GPL-3.0

## Citation

If you use SWAP BMI in your research, please cite:

```
Zawadzki, M. (2024). SWAP Basic Model Interface (BMI) for Python integration.
SWAP 4.2.0. https://github.com/SWAP-model/swap
```
