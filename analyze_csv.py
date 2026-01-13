#!/usr/bin/env python3
"""Analyze SWAP CSV output for annual rainfall and GWL statistics."""

import pandas as pd
import sys

csv_file = '/home/zawadzkim/Code/swap/tests/cases/1.hupselbrook/result_output.csv'

try:
    # Skip the first 6 header lines (metadata)
    df = pd.read_csv(csv_file, skiprows=6)
    print(f"✓ Loaded {len(df)} rows from {csv_file}")
    print(f"Columns: {list(df.columns)}\n")
    
    # Check if we have required columns
    if 'DATETIME' not in df.columns and 'Date' not in df.columns:
        print("Available columns:", df.columns.tolist())
        sys.exit(1)
    
    # Extract year from datetime
    date_col = 'DATETIME' if 'DATETIME' in df.columns else 'Date'
    df['Year'] = pd.to_datetime(df[date_col], errors='coerce').dt.year
    
    # Variables to analyze
    vars_to_analyze = ['RAIN', 'IRRIG', 'INTERC', 'RUNOFF', 'EPOT', 'EACT', 
                       'DRAINAGE', 'QBOTTOM', 'GWL', 'TPOT', 'TACT', 'DSTOR']
    
    # Calculate annual statistics - sum for fluxes, mean for state variables (GWL)
    agg_dict = {}
    for var in vars_to_analyze:
        if var in df.columns:
            # GWL is a state variable, use mean; others are fluxes, use sum
            agg_dict[var] = 'mean' if var == 'GWL' else 'sum'
    
    annual = df.groupby('Year').agg(agg_dict).round(2)
    
    print("=" * 160)
    print("ANNUAL STATISTICS (cm)")
    print("=" * 160)
    
    # Header
    header = f"{'Year':<6}"
    for var in vars_to_analyze:
        if var in annual.columns:
            header += f" {var:>10}"
    print(header)
    print("-" * 160)
    
    # Data rows
    for year, row in annual.iterrows():
        line = f"{int(year):<6}"
        for var in vars_to_analyze:
            if var in annual.columns:
                line += f" {row[var]:>10.2f}"
        print(line)
    
    print("-" * 160)
    
    # Mean row
    line = f"{'Mean':<6}"
    for var in vars_to_analyze:
        if var in annual.columns:
            line += f" {annual[var].mean():>10.2f}"
    print(line)
    
    # Total row (for flux variables only)
    line = f"{'Total':<6}"
    for var in vars_to_analyze:
        if var in annual.columns:
            if var == 'GWL':
                line += f" {'-':>10}"
            else:
                line += f" {annual[var].sum():>10.2f}"
    print(line)
    print("=" * 160)
    
except Exception as e:
    print(f"ERROR: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
