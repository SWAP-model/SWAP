"""Regression checks for SWAP output CSV files.

Runs test cases in isolated temp directories, aggregates the
`result_output.csv` files, and compares annual stats against stored fixtures.
Fails with a non-zero exit if values differ beyond tolerance.
"""

import csv
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor
import concurrent.futures
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

    input_files: dict[str, str] = {}


# Registered test cases
CASES = {
    "hupselbrook": CaseConfig(
        name="hupselbrook",
        case_dir="1.hupselbrook",
        fixture="hupselbrook_expected_gfortran.json",
        flux_vars=["RAIN", "IRRIG", "INTERC", "RUNOFF", "EPOT", "EACT",
                   "DRAINAGE", "QBOTTOM", "TPOT", "TACT", "DSTOR"],
        state_vars=["GWL"],
        input_files={
                "swp": "swap_linux.swp.template",
                "metfile": "283.csv",
                "cropfiles": ["grassd.crp", "maizes.crp", "potatod.crp"],
                "drafile": "swap.dra"
        }
    ),
    # Case 3 (macroporeflow) excluded from check-full per ADR 0011.
    # The case directory still exists at tests/swap-cases/3.macroporeflow/
    # for archival; macropore_config_t is orphan infrastructure per ADR 0010.
    # When the future macropore phase reactivates the module, this entry
    # comes back.
    "grassgrowth": CaseConfig(
        name="grassgrowth",
        case_dir="2.grassgrowth",
        fixture="grassgrowth_expected_gfortran.json",
        flux_vars=[],
        state_vars=[],
        cumul_vars=["PGRASSDM", "GRASSDM", "PMOWDM", "MOWDM"],
    ),
    "oxygenstress": CaseConfig(
        name="oxygenstress",
        case_dir="4.oxygenstress",
        fixture="oxygenstress_expected_gfortran.json",
        flux_vars=[],
        state_vars=["TREDDRY", "TREDWET"],
        cumul_vars=["PGRASSDM", "GRASSDM", "PMOWDM", "MOWDM"],
    ),
    "salinitystress": CaseConfig(
        name="salinitystress",
        case_dir="5.salinitystress",
        fixture="salinitystress_expected_gfortran.json",
        flux_vars=[],
        state_vars=["TREDDRY", "TREDWET", "TREDSOL", "CPWSO", "CWSO",
                    "CONC[-5.0]", "CONC[-25.0]", "CONC[-55.0]"],
    ),
    "surfacewater": CaseConfig(
        name="surfacewater",
        case_dir="6.surfacewater",
        fixture="surfacewater_expected_gfortran.json",
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

def load_case(case_name: str):
    """load case with pyswap
    Currently in development. Problem with pyswap is now that it cannot auto-detect all config files and they have to by specified manually. Also, pyswap does not support
    the detailed rain files and meteo files with .YYY extension.
    """
    import pyswap as psp
    case_dir = TESTS_DIR / "swap-cases" / case_name
    if not case_dir.exists():
        raise FileNotFoundError(f"Case directory not found: {case_dir}")
    
    meta = psp.components.Metadata(
        project="SWAP Regression Tests",
        author="Test Author",
        email="test@email.com",
        institution="Test Institution",
        description=f"Test case for {case_name}",
        swap_ver="4.2.0"
    )

    files = {
        
    }

    met = psp.load_met(case_dir / "met.csv")
    grassd = psp.load_crp(case_dir / "grassd.crp")
    maizes = psp.load_crp(case_dir / "maizes.crp")
    potatod = psp.load_crp(case_dir / "potatod.crp")
    drainage = psp.load_dra(case_dir / "drainage.dra")
    ml: psp.Model = psp.load_swp(case_dir / "swap_linux.swp.template", meta)

    ml.crop.cropfiles = {
        "grassd": grassd,
        "maizes": maizes,
        "potatod": potatod
    }
    ml.lateraldrainage.drafile = drainage
    ml.meteorology.metfile = met

    return ml

def _run_and_aggregate(case: CaseConfig):
    """Run a single SWAP case in a temp dir and return aggregated stats.

    Returns a tuple ``(annual, totals, means)``. Raises RuntimeError on any
    runtime/output failure so callers can surface the message cleanly.
    """
    case_dir = TESTS_DIR / "swap-cases" / case.case_dir
    if not case_dir.exists():
        raise RuntimeError(f"case directory not found at {case_dir}")

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

        # Phase 4f: SWAP now requires swap.toml as the canonical entry
        # point. Stage it from tests/swap-cases/toml/<case>/ alongside
        # cross-file siblings (swap.dra.toml, *.crp.toml). Cases without
        # a populated TOML directory still fail loudly — Phase 4f-extend
        # ports them one by one. Mirrors run_case.sh --toml.
        toml_src = TESTS_DIR / "swap-cases" / "toml" / case.case_dir
        if (toml_src / "swap.toml").exists():
            shutil.copy(toml_src / "swap.toml", workdir / "swap.toml")
            for extra in ("swap.dra.toml",):
                if (toml_src / extra).exists():
                    shutil.copy(toml_src / extra, workdir / extra)
            for crp in toml_src.glob("*.crp.toml"):
                shutil.copy(crp, workdir / crp.name)

        # Record time before running to verify output is fresh
        before_run = time.time()

        # run swap
        proc = subprocess.run([str(SWAP_BIN)], cwd=workdir, capture_output=True, text=True)

        # Check exit code (swap_main exits with 100 on success)
        if proc.returncode != 100:
            detail = f"exit code {proc.returncode} (expected 100)"
            if proc.stdout:
                detail += f"\nstdout:\n{proc.stdout}"
            if proc.stderr:
                detail += f"\nstderr:\n{proc.stderr}"
            raise RuntimeError(f"swap failed: {detail}")

        csv_path = workdir / "result_output.csv"
        if not csv_path.exists():
            raise RuntimeError("result_output.csv not produced")

        # Verify the CSV was created by this run (not a pre-existing file)
        if csv_path.stat().st_mtime < before_run:
            raise RuntimeError(
                f"result_output.csv exists but was not created by this run "
                f"(file mtime {csv_path.stat().st_mtime} < run start {before_run})"
            )

        cumul_vars = case.cumul_vars if hasattr(case, 'cumul_vars') else []
        return aggregate(csv_path, case.flux_vars, case.state_vars, cumul_vars)


def run_case(case: CaseConfig) -> tuple[bool, float]:
    """Run a single test case. Returns (success, execution_time) tuple."""
    start_time = time.perf_counter()

    fixture_path = TESTS_DIR / "regression" / case.fixture

    if not fixture_path.exists():
        print(f"✗ {case.name}: fixture not found at {fixture_path}")
        elapsed = time.perf_counter() - start_time
        return False, elapsed

    expected = load_fixture(fixture_path)

    try:
        annual, totals, means = _run_and_aggregate(case)
    except RuntimeError as exc:
        print(f"✗ {case.name}: {exc}")
        elapsed = time.perf_counter() - start_time
        return False, elapsed

    try:
        compare(expected, annual, totals, means)
    except AssertionError as e:
        print(f"✗ {case.name}: {e}")
        elapsed = time.perf_counter() - start_time
        return False, elapsed

    elapsed = time.perf_counter() - start_time
    print(f"✓ {case.name}: regression ok (annual stats match fixture) [{elapsed:.2f}s]")
    return True, elapsed


def _regen_one_case(case: CaseConfig) -> Path:
    """Regenerate a gfortran-baseline fixture for a single case.

    Writes ``{case.name}_expected_gfortran.json`` next to the other fixtures
    and returns the output path.
    """
    annual, totals, means = _run_and_aggregate(case)
    payload = {
        "years": annual,
        "total": totals,
        "mean": means,
    }
    out_path = TESTS_DIR / "regression" / f"{case.name}_expected_gfortran.json"
    with out_path.open("w") as f:
        json.dump(payload, f, indent=2, sort_keys=True)
        f.write("\n")
    print(f"✓ {case.name}: wrote {out_path.name}")
    return out_path


def main():
    overall_start = time.perf_counter()

    # Parse command line
    args = sys.argv[1:]
    regenerate = False
    if args and args[0] == "--regenerate-fixtures":
        regenerate = True
        args = args[1:]

    if not SWAP_BIN.exists():
        raise SystemExit(f"swap binary not found at {SWAP_BIN}; build first (pixi run build-linux)")

    # Parse command line to select cases
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

    if regenerate:
        print(f"Regenerating fixtures for {len(selected)} case(s) as *_expected_gfortran.json ...\n")
        regen_count = 0
        for case in selected:
            try:
                _regen_one_case(case)
                regen_count += 1
            except Exception as exc:
                print(f"✗ {case.name}: regeneration failed: {exc}")
        print(f"\nRegenerated {regen_count}/{len(selected)} fixture(s).")
        if regen_count != len(selected):
            sys.exit(1)
        return

    # Determine number of workers (defaults to CPU count)
    max_workers = min(len(selected), os.cpu_count() or 1)
    
    print(f"Running {len(selected)} test case(s) with {max_workers} worker(s)...\n")
    
    # Run cases in parallel
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        future_to_case = {executor.submit(run_case, case): case for case in selected}
        
        passed = 0
        failed = 0
        timings = []
        
        # Process results as they complete
        for future in concurrent.futures.as_completed(future_to_case):
            case = future_to_case[future]
            try:
                success, elapsed = future.result()
                timings.append((case.name, elapsed))
                if success:
                    passed += 1
                else:
                    failed += 1
            except Exception as exc:
                print(f'✗ {case.name} generated an exception: {exc}')
                failed += 1

    overall_elapsed = time.perf_counter() - overall_start
    
    # Print summary with timing information
    print(f"\n{'='*60}")
    print(f"Results: {passed} passed, {failed} failed")
    print(f"Total execution time: {overall_elapsed:.2f}s")
    
    if timings:
        print(f"\nIndividual test timings:")
        for name, elapsed in sorted(timings, key=lambda x: x[1], reverse=True):
            print(f"  {name:20s} {elapsed:6.2f}s")
    
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
