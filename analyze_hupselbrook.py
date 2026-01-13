"""
Run SWAP Hupselbrook simulation and analyze rainfall and groundwater levels per year.
"""

import subprocess
import pandas as pd
from pathlib import Path
import sys

# Paths
swap_dir = Path("/home/zawadzkim/Code/swap/swap_org/cases/1.hupselbrook")
swap_exe = Path("/home/zawadzkim/Code/swap/builddir/swap")

if not swap_exe.exists():
    print(f"ERROR: SWAP executable not found at {swap_exe}")
    sys.exit(1)

if not swap_dir.exists():
    print(f"ERROR: Hupselbrook directory not found at {swap_dir}")
    sys.exit(1)

print(f"Running SWAP simulation in: {swap_dir}")
print(f"Using executable: {swap_exe}")
print()

# Run SWAP
try:
    result = subprocess.run(
        [str(swap_exe)],
        cwd=swap_dir,
        capture_output=True,
        text=True,
        timeout=30
    )
    
    if result.returncode != 0:
        print("SWAP execution failed!")
        print("STDOUT:", result.stdout)
        print("STDERR:", result.stderr)
        sys.exit(1)
    
    print("✓ SWAP simulation completed successfully")
    print()
    
except subprocess.TimeoutExpired:
    print("ERROR: SWAP simulation timed out (>30s)")
    sys.exit(1)
except Exception as e:
    print(f"ERROR running SWAP: {e}")
    sys.exit(1)

# Find CSV output file
csv_files = list(swap_dir.glob("*.csv"))
if not csv_files:
    print("ERROR: No CSV file found after SWAP run")
    print("Available files:")
    for f in swap_dir.iterdir():
        if f.is_file():
            print(f"  {f.name}")
    sys.exit(1)

csv_file = csv_files[0]
print(f"Reading CSV output: {csv_file.name}")
print()

# Read CSV
try:
    df = pd.read_csv(csv_file)
    print(f"✓ Loaded {len(df)} rows from CSV")
    print(f"Columns: {', '.join(df.columns.tolist())}")
    print()
except Exception as e:
    print(f"ERROR reading CSV: {e}")
    sys.exit(1)

# Check for required columns
if 'YEAR' not in df.columns:
    # Try to extract year from date column
    date_cols = [c for c in df.columns if 'date' in c.lower() or 'time' in c.lower()]
    if date_cols:
        print(f"Note: Using date column '{date_cols[0]}' to extract year")
        df['YEAR'] = pd.to_datetime(df[date_cols[0]]).dt.year
    else:
        print("ERROR: Cannot find YEAR or date column in CSV")
        print(f"Available columns: {df.columns.tolist()}")
        sys.exit(1)

# Find rain column
rain_cols = [c for c in df.columns if 'rain' in c.lower() or 'prec' in c.lower()]
if not rain_cols:
    print("ERROR: No rainfall column found in CSV")
    print(f"Available columns: {df.columns.tolist()}")
    sys.exit(1)

rain_col = rain_cols[0]
print(f"Using rainfall column: {rain_col}")

# Find GWL column
gwl_cols = [c for c in df.columns if 'gwl' in c.lower() or 'groundwater' in c.lower()]
if gwl_cols:
    gwl_col = gwl_cols[0]
    print(f"Using GWL column: {gwl_col}")
else:
    gwl_col = None
    print("Note: No GWL column found")

print()
print("=" * 70)
print("ANNUAL SUMMARY")
print("=" * 70)

# Group by year
yearly = df.groupby('YEAR').agg({
    rain_col: 'sum',
    gwl_col: 'mean' if gwl_col else lambda x: None
})

print(f"\n{'Year':<8} {'Rainfall (cm)':<18} {'Rainfall (mm)':<18} {'Mean GWL (cm)':<15}")
print("-" * 70)

for year, row in yearly.iterrows():
    rain_cm = row[rain_col]
    rain_mm = rain_cm * 10
    
    if gwl_col:
        gwl_val = row[gwl_col]
        print(f"{year:<8} {rain_cm:>15.2f}    {rain_mm:>15.1f}    {gwl_val:>12.2f}")
    else:
        print(f"{year:<8} {rain_cm:>15.2f}    {rain_mm:>15.1f}")

print("-" * 70)

# Overall statistics
total_rain = yearly[rain_col].sum()
mean_rain = yearly[rain_col].mean()

print(f"{'Mean':<8} {mean_rain:>15.2f}    {mean_rain*10:>15.1f}", end="")
if gwl_col:
    mean_gwl = yearly[gwl_col].mean()
    print(f"    {mean_gwl:>12.2f}")
else:
    print()

print(f"{'Total':<8} {total_rain:>15.2f}    {total_rain*10:>15.1f}")
print(f"\nNumber of years: {len(yearly)}")
print("=" * 70)

print("\n✅ Analysis complete!")
