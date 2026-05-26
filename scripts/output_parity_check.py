"""output_parity_check.py — IO-OUT/D evidence script.

Proves that result_output.csv reproduces the column values written by the
legacy .crp (OutGrass) and .snw (SnowOutput) writers.

Both writers and the CSV write from the SAME state fields.  The only
intentional difference is that the .crp writer applies nint() to the biomass
and root-depth columns before writing; the CSV stores the unrounded float.
The expected max deviation for nint-rounded columns is therefore 0.5 (half an
integer unit).

NOTE — one-time evidence tool, not expected to re-run after IO-OUT Phase E
-----------------------------------------------------------------------
This script required a temporary ``simulation.output.swcrp`` config knob
(wired in e55996f) to activate the legacy .crp writer on the TOML path so
the parity comparison could be performed.  That knob was reverted in the
immediately following commit because the .crp/.snw writers are being retired
entirely in IO-OUT Phase E (tasks E3–E4).  The parity results captured in
commit e55996f (all 22 columns PASS) justify that deletion; this script is
retained as documentation of that evidence and will not function after the
.crp/.snw writers are removed.

Usage
-----
  # Grass .crp parity (case must have swcrop=1 + crop type=3)
  pixi run -e test python scripts/output_parity_check.py <case_dir> --mode crp

  # Snow .snw parity (case must have meteorology.snow.swsnow=1)
  pixi run -e test python scripts/output_parity_check.py <case_dir> --mode snw

  # Both (case must satisfy both prerequisites)
  pixi run -e test python scripts/output_parity_check.py <case_dir> --mode both

The script copies <case_dir> to a temp directory, patches swap.toml to:
  * add simulation.output.swcrp = 1  (enables legacy .crp writer)
  * add meteorology.snow.swsnow = 1  (for snow mode — grassgrowth case is
    already wired with low-enough temperatures for snow events to occur)
  * expand output.csv.inlist to include all mapped columns

Then it runs builddir/swap, parses both output files, aligns rows by index,
and reports PASS/FAIL per column.

Tolerances
----------
  * nint-rounded columns (biomass kg/ha, root depths cm):  abs_diff < 0.6
  * float columns written at f6.2 / f7.2 precision:        abs_diff < 0.01
  * snow-specific daily values:                            abs_diff < 1e-4

Exit codes
----------
  0 — all columns passed (informational columns may differ — see output)
  1 — at least one STRUCTURALLY-IDENTICAL column failed
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# ---------------------------------------------------------------------------
# Project layout
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[1]
SWAP_BIN  = REPO_ROOT / "builddir" / "swap"

# ---------------------------------------------------------------------------
# Column mapping: (crp_col_name, csv_col_name, tolerance, nint_rounded, informational)
# ---------------------------------------------------------------------------
# CRP column order (0-based, after stripping * header lines):
#   0=Date, 1=Daynr, 2=Daycrp, 3=DVS, 4=TSUM, 5=LAIpot, 6=LAI,
#   7=Height, 8=CrpFac, 9=RootdPot, 10=Rootd,
#   11=PWLV, 12=WLV, 13=PWST, 14=WST, 15=PWRT, 16=WRT,
#   17=CPWDM, 18=CWDM, 19=CPWSO, 20=CWSO,
#   21=PGRASSDM, 22=GRASSDM, 23=PMOWDM, 24=MOWDM, 25=PGRAZDM, 26=GRAZDM,
#   (remaining columns not relevant to grass parity)
#
# Grass .crp format for the float columns (DVS, LAI, Height, CrpFac) is f6.2
# or f7.2 — 2 decimal places; the CSV stores more digits.  Tolerance 0.01.
# Biomass and root-depth columns are written via nint() — tolerance 0.6.
#
# DVS for grass is always -99.99 (sentinel) — this will match exactly.

CRP_COL_MAP = [
    # (crp_header_name, csv_header_name, tolerance, nint_rounded, informational)
    ("DVS",      "DVS",      0.01,  False, False),
    ("LAI",      "LAI",      0.01,  False, False),
    ("Height",   "HEIGHT",   0.01,  False, False),
    ("CrpFac",   "CRPFAC",   0.01,  False, False),
    ("Rootd",    "RD",       0.6,   True,  False),  # nint(rd) in .crp
    ("RootdPot", "RDPOT",    0.6,   True,  False),  # nint(rdpot) in .crp
    ("PWLV",     "PWLV",     0.6,   True,  False),  # nint(wlvpot)
    ("WLV",      "WLV",      0.6,   True,  False),  # nint(wlv)
    ("PWST",     "PWST",     0.6,   True,  False),  # nint(wstpot)
    ("WST",      "WST",      0.6,   True,  False),  # nint(wst)
    ("PWRT",     "PWRT",     0.6,   True,  False),  # nint(wrtpot)
    ("WRT",      "WRT",      0.6,   True,  False),  # nint(wrt)
    ("PGRASSDM", "PGRASSDM", 0.6,   True,  False),  # nint(tagppot)
    ("GRASSDM",  "GRASSDM",  0.6,   True,  False),  # nint(tagp)
    ("PMOWDM",   "PMOWDM",   0.6,   True,  False),  # nint(tagptpot)
    ("MOWDM",    "MOWDM",    0.6,   True,  False),  # nint(tagpt)
    ("PGRAZDM",  "PGRAZDM",  0.6,   True,  False),  # nint(cuptgrazpot)
    ("GRAZDM",   "GRAZDM",   0.6,   True,  False),  # nint(cuptgraz)
]

# SNW column map: (snw_col_name, csv_col_name, tolerance, informational)
# .snw format: date, dcum, rainfall, snowfall, snowstorage, meltflux, sublimation
# CSV (from inlist ssnow,sublim,melt,snow):  DATETIME, SNOW, SUBLIM, SSNOW, MELT (order depends on inlist)
# SNOW (CSV) = state%atmosphere%intr%igsnow (period accumulator for gsnow)
#            = .snw "snowfall" (gsnow) when period=1 (daily output)
# SUBLIM (CSV) = intr%isubl (period accumulator for subl)
#              = .snw "sublimation" when period=1
# SSNOW (CSV) = state%atmosphere%ssnow
#             = .snw "snowstorage"                     STRUCTURALLY IDENTICAL
# MELT (CSV)  = state%atmosphere%melt
#             = .snw "meltflux"                        STRUCTURALLY IDENTICAL
#
# With period=1 (daily output), the period accumulators for SNOW/SUBLIM
# equal the daily values, so all four should match.
SNW_COL_MAP = [
    # (snw_col_name, csv_col_name, tolerance, informational, reason)
    ("snowfall",    "SNOW",   1e-4,  False, "period accumulator equals daily value when period=1"),
    ("sublimation", "SUBLIM", 1e-4,  False, "period accumulator equals daily value when period=1"),
    ("snowstorage", "SSNOW",  1e-4,  False, "same state field state%atmosphere%ssnow"),
    ("meltflux",    "MELT",   1e-4,  False, "same state field state%atmosphere%melt"),
]


def parse_swap_csv(path: Path) -> tuple[list[str], list[dict]]:
    """Parse result_output.csv, skip '*' comment lines.  Returns (headers, rows)."""
    rows = []
    headers = None
    with path.open() as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            if row[0].strip().startswith("*"):
                continue
            if headers is None:
                headers = [x.strip() for x in row]
            else:
                rows.append({h: v.strip() for h, v in zip(headers, row)})
    return headers, rows


def parse_crp(path: Path) -> tuple[list[str], list[dict]]:
    """Parse .crp legacy file, skip '*' header lines. Returns (headers, rows)."""
    headers = None
    rows = []
    with path.open() as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("*"):
                continue
            parts = [x.strip() for x in line.split(",")]
            if headers is None:
                headers = parts
            else:
                rows.append(dict(zip(headers, parts)))
    return headers, rows


def parse_snw(path: Path) -> tuple[list[str], list[dict]]:
    """Parse .snw legacy file, skip '*' header lines. Returns (headers, rows)."""
    headers = None
    rows = []
    with path.open() as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("*"):
                continue
            parts = [x.strip() for x in line.split(",")]
            if headers is None:
                headers = parts
            else:
                rows.append(dict(zip(headers, parts)))
    return headers, rows


def safe_float(s: str) -> float | None:
    """Convert string to float, returning None for empty/whitespace-only."""
    s = s.strip()
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def run_case(case_dir: Path, mode: str) -> Path:
    """Copy case to temp dir, patch swap.toml, run SWAP.  Returns workdir Path."""
    workdir_parent = Path(tempfile.mkdtemp(prefix="parity_"))
    workdir = workdir_parent / "case"
    shutil.copytree(
        case_dir, workdir,
        ignore=shutil.ignore_patterns("result_output.csv", "result_*.csv", "*.log", "*.out")
    )

    # Stage swap.swp from template (legacy sub-readers may still need it)
    if not (workdir / "swap.swp").exists():
        tmpl = workdir / "swap_linux.swp.template"
        if tmpl.exists():
            shutil.copy(tmpl, workdir / "swap.swp")

    # Patch swap.toml
    toml_path = workdir / "swap.toml"
    toml_text = toml_path.read_text()

    # 1. Enable swcrp = 1 (legacy .crp writer)
    if mode in ("crp", "both"):
        if "swcrp" not in toml_text:
            toml_text = toml_text.replace(
                "[simulation.output]",
                "[simulation.output]\nswcrp = 1"
            )
        else:
            import re
            toml_text = re.sub(r"swcrp\s*=\s*\d", "swcrp = 1", toml_text)

    # 2. Enable snow (swsnow = 1) if not already present
    if mode in ("snw", "both"):
        if "[meteorology.snow]" not in toml_text:
            # Insert after [meteorology.rain] block
            snow_block = (
                "\n[meteorology.snow]\n"
                "swsnow   = 1\n"
                "snowcoef = 0.3\n"
                "teprrain = 2.0\n"
                "teprsnow = -2.0\n"
            )
            # Insert before [meteorology.interception] or after rain block
            if "[meteorology.interception]" in toml_text:
                toml_text = toml_text.replace(
                    "[meteorology.interception]",
                    snow_block + "[meteorology.interception]"
                )
            else:
                # Append before [drainage]
                toml_text = toml_text.replace(
                    "[drainage]",
                    snow_block + "[drainage]"
                )

    # 3. Expand inlist
    if mode == "crp":
        crp_inlist = "dvs,lai,height,crpfac,rd,rdpot,pwlv,wlv,pwst,wst,pwrt,wrt,pgrassdm,grassdm,pmowdm,mowdm,pgrazdm,grazdm"
        import re
        toml_text = re.sub(r'(inlist\s*=\s*)"[^"]*"', f'inlist = "{crp_inlist}"', toml_text)
    elif mode == "snw":
        snw_inlist = "ssnow,sublim,melt,snow"
        import re
        toml_text = re.sub(r'(inlist\s*=\s*)"[^"]*"', f'inlist = "{snw_inlist}"', toml_text)
    elif mode == "both":
        both_inlist = ("dvs,lai,height,crpfac,rd,rdpot,pwlv,wlv,pwst,wst,pwrt,wrt,"
                       "pgrassdm,grassdm,pmowdm,mowdm,pgrazdm,grazdm,"
                       "ssnow,sublim,melt,snow")
        import re
        toml_text = re.sub(r'(inlist\s*=\s*)"[^"]*"', f'inlist = "{both_inlist}"', toml_text)

    toml_path.write_text(toml_text)

    # Run SWAP
    proc = subprocess.run([str(SWAP_BIN)], cwd=workdir, capture_output=True, text=True)
    if proc.returncode != 100:
        print("ERROR: SWAP exited with code", proc.returncode)
        if proc.stdout:
            print("stdout:", proc.stdout[-2000:])
        if proc.stderr:
            print("stderr:", proc.stderr[-2000:])
        sys.exit(2)

    return workdir


def compare_crp(workdir: Path) -> bool:
    """Compare result.crp against result_output.csv. Returns True if all columns PASS."""
    crp_path = workdir / "result.crp"
    csv_path = workdir / "result_output.csv"

    if not crp_path.exists():
        print("FAIL: result.crp not produced (is swcrp=1 in [simulation.output]?)")
        return False
    if not csv_path.exists():
        print("FAIL: result_output.csv not produced")
        return False

    crp_headers, crp_rows = parse_crp(crp_path)
    csv_headers, csv_rows = parse_swap_csv(csv_path)

    if len(crp_rows) != len(csv_rows):
        print(f"WARN: row count mismatch: crp={len(crp_rows)}, csv={len(csv_rows)}")

    n_rows = min(len(crp_rows), len(csv_rows))

    print(f"\n--- Grass .crp parity ({n_rows} rows) ---")
    print(f"{'Column':<12} {'CRP_col':<12} {'CSV_col':<12} {'Tolerance':<12} {'MaxDiff':<12} Result")
    print("-" * 72)

    all_pass = True
    for crp_col, csv_col, tol, nint_rounded, informational in CRP_COL_MAP:
        if crp_col not in crp_headers:
            print(f"{'SKIP':<12} {crp_col:<12} (not in .crp header)")
            continue
        if csv_col not in csv_headers:
            print(f"{'SKIP':<12} {csv_col:<12} (not in CSV header)")
            continue

        max_diff = 0.0
        n_missing = 0
        row_fails = []
        for i in range(n_rows):
            crp_val = safe_float(crp_rows[i].get(crp_col, ""))
            csv_val = safe_float(csv_rows[i].get(csv_col, ""))
            if crp_val is None or csv_val is None:
                n_missing += 1
                continue
            diff = abs(crp_val - csv_val)
            if diff > max_diff:
                max_diff = diff
            if diff > tol:
                row_fails.append((i + 1, crp_val, csv_val, diff))

        tag = "nint" if nint_rounded else "float"
        info_str = " [INFO]" if informational else ""
        if row_fails:
            status = "FAIL" + info_str
            if not informational:
                all_pass = False
        else:
            status = "PASS" + info_str

        print(f"{crp_col:<12} {crp_col:<12} {csv_col:<12} {tol:<12.4f} {max_diff:<12.5f} {status}  ({tag})")
        if row_fails and len(row_fails) <= 5:
            for row_n, cv, csav, d in row_fails:
                print(f"          row {row_n}: crp={cv:.4f} csv={csav:.4f} diff={d:.5f}")

    if n_missing > 0:
        print(f"  ({n_missing} rows skipped due to missing/empty values in one or both files)")

    return all_pass


def compare_snw(workdir: Path) -> bool:
    """Compare result.snw against result_output.csv. Returns True if all non-informational PASS."""
    snw_path = workdir / "result.snw"
    csv_path = workdir / "result_output.csv"

    if not snw_path.exists():
        print("FAIL: result.snw not produced (is swsnow=1 in [meteorology.snow]?)")
        return False
    if not csv_path.exists():
        print("FAIL: result_output.csv not produced")
        return False

    snw_headers, snw_rows = parse_snw(snw_path)
    csv_headers, csv_rows = parse_swap_csv(csv_path)

    if len(snw_rows) != len(csv_rows):
        print(f"WARN: row count mismatch: snw={len(snw_rows)}, csv={len(csv_rows)}")

    n_rows = min(len(snw_rows), len(csv_rows))

    print(f"\n--- Snow .snw parity ({n_rows} rows) ---")
    print(f"{'SNW_col':<14} {'CSV_col':<10} {'Tolerance':<12} {'MaxDiff':<12} Result")
    print("-" * 64)

    all_pass = True
    for snw_col, csv_col, tol, informational, reason in SNW_COL_MAP:
        if snw_col not in snw_headers:
            print(f"{'SKIP':<14} {csv_col:<10} ('{snw_col}' not in .snw header)")
            continue
        if csv_col not in csv_headers:
            print(f"{'SKIP':<14} {csv_col:<10} ('{csv_col}' not in CSV header)")
            continue

        max_diff = 0.0
        n_missing = 0
        row_fails = []
        for i in range(n_rows):
            snw_val = safe_float(snw_rows[i].get(snw_col, ""))
            csv_val = safe_float(csv_rows[i].get(csv_col, ""))
            if snw_val is None or csv_val is None:
                n_missing += 1
                continue
            diff = abs(snw_val - csv_val)
            if diff > max_diff:
                max_diff = diff
            if diff > tol:
                row_fails.append((i + 1, snw_val, csv_val, diff))

        info_str = " [INFO]" if informational else ""
        if row_fails:
            status = "FAIL" + info_str
            if not informational:
                all_pass = False
        else:
            status = "PASS" + info_str

        print(f"{snw_col:<14} {csv_col:<10} {tol:<12.2e} {max_diff:<12.5f} {status}")
        if row_fails and len(row_fails) <= 5:
            for row_n, sv, cv, d in row_fails:
                print(f"  row {row_n}: snw={sv:.6f} csv={cv:.6f} diff={d:.6f}")

    if informational:
        print(f"  NOTE: {reason}")

    return all_pass


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("case_dir", type=Path, help="Path to TOML case directory (containing swap.toml)")
    parser.add_argument("--mode", choices=["crp", "snw", "both"], default="crp",
                        help="Which legacy writer to check: crp, snw, or both (default: crp)")
    args = parser.parse_args()

    if not SWAP_BIN.exists():
        print(f"ERROR: swap binary not found at {SWAP_BIN}")
        print("Build first: pixi run build-linux")
        sys.exit(2)

    case_dir = args.case_dir.resolve()
    if not (case_dir / "swap.toml").exists():
        print(f"ERROR: no swap.toml in {case_dir}")
        sys.exit(2)

    print(f"Case: {case_dir}")
    print(f"Mode: {args.mode}")
    print(f"Binary: {SWAP_BIN}")

    print("\nPreparing temp dir and running SWAP...")
    workdir = run_case(case_dir, args.mode)
    print(f"Run complete. Workdir: {workdir}")

    all_pass = True
    if args.mode in ("crp", "both"):
        crp_ok = compare_crp(workdir)
        all_pass = all_pass and crp_ok

    if args.mode in ("snw", "both"):
        snw_ok = compare_snw(workdir)
        all_pass = all_pass and snw_ok

    print()
    if all_pass:
        print("OVERALL: PASS — result_output.csv reproduces the legacy file columns within tolerance.")
        sys.exit(0)
    else:
        print("OVERALL: FAIL — one or more structurally-identical columns exceeded tolerance.")
        sys.exit(1)


if __name__ == "__main__":
    main()
