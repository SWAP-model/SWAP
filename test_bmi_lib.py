"""
Test the SWAP BMI library compilation and basic functionality.
This is a minimal test that doesn't require a full SWAP configuration.
"""

import ctypes
import os
import sys
from pathlib import Path

# Find the library
lib_paths = [
    Path("/home/zawadzkim/Code/swap/builddir/libswap.so"),
    Path("../swap/builddir/libswap.so"),
    Path("./builddir/libswap.so"),
]

lib_path = None
for path in lib_paths:
    if path.exists():
        lib_path = str(path.absolute())
        break

if not lib_path:
    print("ERROR: Could not find libswap.so")
    print("Searched in:")
    for path in lib_paths:
        print(f"  - {path.absolute()}")
    sys.exit(1)

print(f"Found library: {lib_path}")

# Load the library
try:
    lib = ctypes.CDLL(lib_path)
    print("✓ Successfully loaded library")
except Exception as e:
    print(f"✗ Failed to load library: {e}")
    sys.exit(1)

# Check for key BMI functions
functions_to_check = [
    "initialize",
    "initialize_memory",
    "finalize",
    "update",
    "update_day",
    "get_current_time",
    "get_grid_size",
    "get_value_ptr_theta",
    "get_value_gwl",
]

print("\nChecking for BMI functions:")
missing = []
for func_name in functions_to_check:
    try:
        func = getattr(lib, func_name)
        print(f"  ✓ {func_name}")
    except AttributeError:
        print(f"  ✗ {func_name} - NOT FOUND")
        missing.append(func_name)

if missing:
    print(f"\n❌ {len(missing)} functions missing")
    sys.exit(1)
else:
    print(f"\n✅ All {len(functions_to_check)} required BMI functions found")

# Try to set up function signatures
print("\nSetting up function signatures:")
try:
    # Initialize
    lib.initialize.argtypes = [ctypes.c_char_p]
    lib.initialize.restype = ctypes.c_int
    
    lib.initialize_memory.argtypes = []
    lib.initialize_memory.restype = ctypes.c_int
    
    # Lifecycle
    lib.finalize.argtypes = []
    lib.finalize.restype = ctypes.c_int
    
    lib.update_day.argtypes = []
    lib.update_day.restype = ctypes.c_int
    
    # Time
    lib.get_current_time.argtypes = [ctypes.POINTER(ctypes.c_double)]
    lib.get_current_time.restype = ctypes.c_int
    
    # Grid
    lib.get_grid_size.argtypes = [ctypes.POINTER(ctypes.c_int)]
    lib.get_grid_size.restype = ctypes.c_int
    
    # Getters
    lib.get_value_ptr_theta.argtypes = [ctypes.POINTER(ctypes.c_void_p)]
    lib.get_value_ptr_theta.restype = ctypes.c_int
    
    lib.get_value_gwl.argtypes = [ctypes.POINTER(ctypes.c_double)]
    lib.get_value_gwl.restype = ctypes.c_int
    
    print("  ✓ Function signatures configured")
except Exception as e:
    print(f"  ✗ Failed to set function signatures: {e}")
    sys.exit(1)

# Test basic memory initialization (without actual SWAP config files)
print("\nTesting basic functionality:")
try:
    # Note: initialize_memory() may fail without proper SWAP initialization
    # This is expected - we're just testing if the function is callable
    status = lib.initialize_memory()
    if status == 0:
        print("  ✓ initialize_memory() called successfully (status=0)")
        
        # Try to get grid size
        grid_size = ctypes.c_int()
        status = lib.get_grid_size(ctypes.byref(grid_size))
        if status == 0:
            print(f"  ✓ get_grid_size() = {grid_size.value}")
        else:
            print(f"  ⚠ get_grid_size() returned status {status}")
        
        # Finalize
        status = lib.finalize()
        if status == 0:
            print("  ✓ finalize() called successfully")
    else:
        print(f"  ⚠ initialize_memory() returned status {status}")
        print("    (This is expected without proper SWAP configuration)")
        
except Exception as e:
    print(f"  ⚠ Runtime test failed: {e}")
    print("    (This is expected without proper SWAP configuration)")

print("\n" + "="*60)
print("SUMMARY: Library compilation and symbol export successful!")
print("="*60)
print("\nNext steps:")
print("1. Test with pyswap.bmi module: python -c 'from pyswap.bmi import SwapBMI'")
print("2. Run examples: python examples/bmi_examples.py")
print("3. Test with actual SWAP configuration files")
