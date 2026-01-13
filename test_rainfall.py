"""
Test SWAP BMI with actual Hupselbrook simulation and compute annual rainfall.
"""

import ctypes
import sys
from pathlib import Path
from collections import defaultdict

# Find the library
lib_path = Path("/home/zawadzkim/Code/swap/builddir/libswap.so")
if not lib_path.exists():
    print(f"ERROR: Could not find {lib_path}")
    sys.exit(1)

print(f"Loading library: {lib_path}")
lib = ctypes.CDLL(str(lib_path))

# Set up function signatures
lib.initialize.argtypes = [ctypes.c_char_p]
lib.initialize.restype = ctypes.c_int

lib.update_day.argtypes = []
lib.update_day.restype = ctypes.c_int

lib.finalize.argtypes = []
lib.finalize.restype = ctypes.c_int

lib.get_current_time.argtypes = [ctypes.POINTER(ctypes.c_double)]
lib.get_current_time.restype = ctypes.c_int

lib.get_value_cgrai.argtypes = [ctypes.POINTER(ctypes.c_double)]
lib.get_value_cgrai.restype = ctypes.c_int

lib.get_value_cnrai.argtypes = [ctypes.POINTER(ctypes.c_double)]
lib.get_value_cnrai.restype = ctypes.c_int

lib.get_value_year.argtypes = [ctypes.POINTER(ctypes.c_int)]
lib.get_value_year.restype = ctypes.c_int

lib.get_value_daynr.argtypes = [ctypes.POINTER(ctypes.c_int)]
lib.get_value_daynr.restype = ctypes.c_int

lib.get_value_gwl.argtypes = [ctypes.POINTER(ctypes.c_double)]
lib.get_value_gwl.restype = ctypes.c_int

# Find Hupselbrook swap.swp
swap_file = None
search_paths = [
    Path("/home/zawadzkim/Code/swap/swap_org/cases/1.hupselbrook/swap.swp"),
    Path("./swap_org/cases/1.hupselbrook/swap.swp"),
    Path("./tests/cases/1.hupselbrook/swap.swp"),
]

for path in search_paths:
    if path.exists():
        swap_file = str(path.absolute())
        break

if not swap_file:
    print("ERROR: Could not find Hupselbrook.swp")
    print("Searched:")
    for p in search_paths:
        print(f"  {p}")
    sys.exit(1)

print(f"Found SWAP file: {swap_file}")
print()

# Initialize SWAP
print("Initializing SWAP...")
status = lib.initialize(swap_file.encode('utf-8'))
if status != 0:
    print(f"ERROR: initialize() failed with status {status}")
    sys.exit(1)
print("✓ SWAP initialized")

# Track rainfall by year
rainfall_by_year = defaultdict(lambda: {'start': 0.0, 'end': 0.0})
prev_year = None
prev_cgrai = 0.0
day_count = 0

print("\nRunning simulation and tracking rainfall...")
print("=" * 60)

try:
    while True:
        # Get current state
        year = ctypes.c_int()
        daynr = ctypes.c_int()
        cgrai = ctypes.c_double()
        cnrai = ctypes.c_double()
        gwl = ctypes.c_double()
        time = ctypes.c_double()
        
        lib.get_value_year(ctypes.byref(year))
        lib.get_value_daynr(ctypes.byref(daynr))
        lib.get_value_cgrai(ctypes.byref(cgrai))
        lib.get_value_cnrai(ctypes.byref(cnrai))
        lib.get_value_gwl(ctypes.byref(gwl))
        lib.get_current_time(ctypes.byref(time))
        
        year_val = year.value
        
        # Track start of year
        if prev_year is None or year_val != prev_year:
            if prev_year is not None:
                # Save end value for previous year
                rainfall_by_year[prev_year]['end'] = prev_cgrai
            # Save start value for new year
            rainfall_by_year[year_val]['start'] = cgrai.value
            rainfall_by_year[year_val]['year'] = year_val
            
            if prev_year != year_val and prev_year is not None:
                annual_rain = rainfall_by_year[prev_year]['end'] - rainfall_by_year[prev_year]['start']
                print(f"Year {prev_year}: {annual_rain:8.2f} cm  (cumulative: {rainfall_by_year[prev_year]['end']:8.2f} cm)")
        
        prev_year = year_val
        prev_cgrai = cgrai.value
        day_count += 1
        
        # Update one day
        status = lib.update_day()
        if status != 0:
            break
            
except KeyboardInterrupt:
    print("\nSimulation interrupted by user")

# Handle last year
if prev_year is not None:
    rainfall_by_year[prev_year]['end'] = prev_cgrai
    annual_rain = rainfall_by_year[prev_year]['end'] - rainfall_by_year[prev_year]['start']
    print(f"Year {prev_year}: {annual_rain:8.2f} cm  (cumulative: {rainfall_by_year[prev_year]['end']:8.2f} cm)")

print("=" * 60)
print(f"\nSimulation completed: {day_count} days")

# Print summary
print("\n" + "=" * 60)
print("ANNUAL RAINFALL SUMMARY")
print("=" * 60)
print(f"{'Year':<8} {'Rainfall (cm)':<15} {'Rainfall (mm)':<15}")
print("-" * 60)

total_rainfall = 0.0
years = sorted([y for y in rainfall_by_year.keys()])

for year in years:
    data = rainfall_by_year[year]
    annual_rain = data['end'] - data['start']
    total_rainfall += annual_rain
    print(f"{year:<8} {annual_rain:>12.2f}     {annual_rain*10:>12.1f}")

if years:
    avg_rainfall = total_rainfall / len(years)
    print("-" * 60)
    print(f"{'Mean':<8} {avg_rainfall:>12.2f}     {avg_rainfall*10:>12.1f}")
    print(f"{'Total':<8} {total_rainfall:>12.2f}     {total_rainfall*10:>12.1f}")
    print(f"\nNumber of years: {len(years)}")
print("=" * 60)

# Finalize
print("\nFinalizing SWAP...")
lib.finalize()
print("✓ SWAP finalized")
print("\n✅ Analysis complete!")
