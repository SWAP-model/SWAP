"""Regression checks for SWAP output CSV files.


Runs test cases in isolated temp directories, aggregates the
`result_output.csv` files, and compares annual stats against stored fixtures.
Fails with a non-zero exit if values differ beyond tolerance.
"""


import csv
import json
import math
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import NamedTuple


TESTS_DIR = Path(__file__).resolve().parent.parent
SWAP_BIN = Path(__file__).resolve().parents[2] / "builddir" / "swap"
TOL = 1e-2  # cm tolerance on aggregated values



class CaseConfig(NamedTuple):
    """Configuration for a test case."""
    name: str
    case_dir: str  # relative path from tests/cases
    fixture: str   # fixture filename in tests/regression
    flux_vars: list[str]  # summed annually
    state_vars: list[str]  # averaged annually
    cumul_vars: list[str] = []  # last value per year (for cumulative outputs)



# Registered test cases
CASES = {
    "hupselbrook": CaseConfig(
        name="hupselbrook",
        case_dir="1.hupselbrook",
        fixture="hupselbrook_expected.json",
        flux_vars=["RAIN", "IRRIG", "INTERC", "RUNOFF", "EPOT", "EACT",
                   "DRAINAGE", "QBOTTOM", "TPOT", "TACT", "DSTOR"],
        state_vars=["GWL"],
    ),
    "macropore": CaseConfig(
        name="macropore",
        case_dir="3.macroporeflow",
        fixture="macropore_expected.json",
        flux_vars=["DRAINAGE"],
        state_vars=["GWL"],
    ),
    "grassgrowth": CaseConfig(
        name="grassgrowth",
        case_dir="2.grassgrowth",
        fixture="grassgrowth_expected.json",
        flux_vars=[],
        state_vars=[],
        cumul_vars=["PGRASSDM", "GRASSDM", "PMOWDM", "MOWDM"],
    ),
    "oxygenstress": CaseConfig(
        name="oxygenstress",
        case_dir="4.oxygenstress",
        fixture="oxygenstress_expected.json",
        flux_vars=[],
        state_vars=["TREDDRY", "TREDWET"],
        cumul_vars=["PGRASSDM", "GRASSDM", "PMOWDM", "MOWDM"],
    ),
    "salinitystress": CaseConfig(
        name="salinitystress",
        case_dir="5.salinitystress",
        fixture="salinitystress_expected.json",
        flux_vars=[],
        state_vars=["TREDDRY", "TREDWET", "TREDSOL", "CPWSO", "CWSO",
                    "CONC[-5.0]", "CONC[-25.0]", "CONC[-55.0]"],
    ),
    "surfacewater": CaseConfig(
        name="surfacewater",
        case_dir="6.surfacewater",
        fixture="surfacewater_expected.json",
        flux_vars=[],
        state_vars=["GWL", "POND"],
    ),
}



def load_fixture(path: Path):
    with path.open() as f:
        return json.load(f)



def aggregate(csv_path: Path, flux_vars: list[str], state_vars: list[str], cumul_vars: list[str] = None):
    """Aggregate daily CSV into annual stats.

    - flux_vars: summed annually
    - state_vars: averaged annually  
    - cumul_vars: last value per year (for cumulative outputs)
    """
    if cumul_vars is None:
        cumul_vars = []
    all_vars = flux_vars + state_vars + cumul_vars
    years = {}
    with csv_path.open() as f:
        reader = csv.reader(f)
        # skip header lines starting with '*'
        for row in reader:
            if row and row[0].startswith("*"):
                continue
            headers = row
            break
        data = list(csv.DictReader(f, fieldnames=headers))


    for rec in data:
        if not rec.get("DATETIME"):
            continue
        year = rec["DATETIME"].split("-")[0]
        yr = years.setdefault(year, {k: [] for k in all_vars})
        for k in all_vars:
            if k in rec and rec[k]:
                yr[k].append(float(rec[k]))


    annual = {}
    for year, vals in years.items():
        annual[year] = {}
        for k in flux_vars:
            if vals[k]:
                annual[year][k] = round(sum(vals[k]), 2)
        for k in state_vars:
            if vals[k]:
                annual[year][k] = round(sum(vals[k]) / len(vals[k]), 2)
        for k in cumul_vars:
            if vals[k]:
                annual[year][k] = round(vals[k][-1], 2)  # last value


    # totals/means across years
    totals = {}
    means = {}
    n_years = len(annual)
    for k in flux_vars:
        totals[k] = round(sum(annual[y].get(k, 0.0) for y in annual), 2)
        means[k] = round(totals[k] / n_years, 2)
    for k in state_vars:
        means[k] = round(sum(annual[y].get(k, 0.0) for y in annual) / n_years, 2)
    for k in cumul_vars:
        means[k] = round(sum(annual[y].get(k, 0.0) for y in annual) / n_years, 2)


    return annual, totals, means



def compare(expected, actual_years, actual_totals, actual_means):
    """Compare expected vs actual values and collect all mismatches."""
    mismatches = []

    def check_block(block_name, exp_block, act_block):
        for year_or_key, exp_vals in exp_block.items():
            if isinstance(exp_vals, dict):
                act_vals = act_block.get(year_or_key, {})
                for var, exp_val in exp_vals.items():
                    act_val = act_vals.get(var)
                    matches = act_val is not None and math.isclose(act_val, exp_val, abs_tol=TOL)
                    if not matches:
                        mismatches.append({
                            "block": block_name,
                            "year": year_or_key,
                            "var": var,
                            "expected": exp_val,
                            "actual": act_val,
                            "diff": abs(act_val - exp_val) if act_val is not None else None
                        })
            else:
                act_val = act_block.get(year_or_key)
                matches = act_val is not None and math.isclose(act_val, exp_vals, abs_tol=TOL)
                if not matches:
                    mismatches.append({
                        "block": block_name,
                        "year": year_or_key,
                        "var": "-",
                        "expected": exp_vals,
                        "actual": act_val,
                        "diff": abs(act_val - exp_vals) if act_val is not None else None
                    })


    check_block("years", expected["years"], actual_years)
    if "total" in expected:
        check_block("total", expected["total"], actual_totals)
    if "mean" in expected:
        check_block("mean", expected["mean"], actual_means)

    if mismatches:
        # Build comparison table
        lines = ["\n  Mismatches found (tolerance={:.0e}):".format(TOL)]
        lines.append("  {:>6} {:>12} {:>14} {:>14} {:>12}".format(
            "Year", "Variable", "Expected", "Actual", "Diff"))
        lines.append("  " + "-" * 60)
        for m in mismatches:
            diff_str = f"{m['diff']:.4f}" if m['diff'] is not None else "N/A"
            act_str = f"{m['actual']:.4f}" if m['actual'] is not None else "None"
            lines.append("  {:>6} {:>12} {:>14.4f} {:>14} {:>12}".format(
                m['year'], m['var'], m['expected'], act_str, diff_str))

        raise AssertionError("\n".join(lines))



def run_case(case: CaseConfig) -> bool:
    """Run a single test case. Returns True on success, False on failure."""
    case_dir = TESTS_DIR / "cases" / case.case_dir
    fixture_path = TESTS_DIR / "regression" / case.fixture


    if not case_dir.exists():
        print(f"✗ {case.name}: case directory not found at {case_dir}")
        return False


    if not fixture_path.exists():
        print(f"✗ {case.name}: fixture not found at {fixture_path}")
        return False


    expected = load_fixture(fixture_path)


    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)

        # Copy case files, excluding any pre-existing output files
        shutil.copytree(
            case_dir, 
            tmp / "case",
            ignore=shutil.ignore_patterns(
                'result_output.csv',
                'result_*.csv',
                '*.log',
                'output.*',
                '*.out'
            )
        )
        workdir = tmp / "case"


        # ensure template is named swap.swp
        swap_file = workdir / "swap.swp"
        if not swap_file.exists():
            template = workdir / "swap_linux.swp.template"
            if template.exists():
                shutil.copy(template, swap_file)
            else:
                # try other common names
                for alt in workdir.glob("*.swp"):
                    shutil.copy(alt, swap_file)
                    break


        # Record time before running to verify output is fresh
        before_run = time.time()


        # run swap
        proc = subprocess.run([str(SWAP_BIN)], cwd=workdir, capture_output=True, text=True)

        # Check exit code (swap_main exits with 100 on success)
        if proc.returncode != 100:
            print(f"✗ {case.name}: swap failed with exit code {proc.returncode} (expected 100)")
            if proc.stdout:
                print(f"stdout:\n{proc.stdout}")
            if proc.stderr:
                print(f"stderr:\n{proc.stderr}")
            return False


        csv_path = workdir / "result_output.csv"
        if not csv_path.exists():
            print(f"✗ {case.name}: result_output.csv not produced")
            return False


        # Verify the CSV was created by this run (not a pre-existing file)
        if csv_path.stat().st_mtime < before_run:
            print(f"✗ {case.name}: result_output.csv exists but was not created by this run")
            print(f"   File timestamp: {csv_path.stat().st_mtime}, run started at: {before_run}")
            return False


        try:
            cumul_vars = case.cumul_vars if hasattr(case, 'cumul_vars') else []
            annual, totals, means = aggregate(csv_path, case.flux_vars, case.state_vars, cumul_vars)
            compare(expected, annual, totals, means)
        except AssertionError as e:
            print(f"✗ {case.name}: {e}")
            return False


    print(f"✓ {case.name}: regression ok (annual stats match fixture)")
    return True



def main():
    if not SWAP_BIN.exists():
        raise SystemExit(f"swap binary not found at {SWAP_BIN}; build first (pixi run build-linux)")


    # Parse command line to select cases
    args = sys.argv[1:]
    if args:
        selected = []
        for arg in args:
            if arg in CASES:
                selected.append(CASES[arg])
            else:
                print(f"Unknown case: {arg}. Available: {', '.join(CASES.keys())}")
                sys.exit(1)
    else:
        selected = list(CASES.values())


    passed = 0
    failed = 0
    for case in selected:
        if run_case(case):
            passed += 1
        else:
            failed += 1


    print(f"\n{passed} passed, {failed} failed")
    if failed:
        sys.exit(1)



if __name__ == "__main__":
    main()