#!/usr/bin/env python3
"""Regression check for SWAP output CSV (Hupselbrook case).

Runs the bundled test case in an isolated temp directory, aggregates the
`result_output.csv` file, and compares annual stats against a stored fixture.
Fails with a non-zero exit if values differ beyond tolerance.
"""

import csv
import json
import math
import shutil
import subprocess
import tempfile
from pathlib import Path

CASE_DIR = Path(__file__).resolve().parent.parent / "cases" / "1.hupselbrook"
FIXTURE = Path(__file__).resolve().parent / "hupselbrook_expected.json"
SWAP_BIN = Path(__file__).resolve().parents[2] / "builddir" / "swap"
TOL = 1e-2  # cm tolerance on aggregated values

# Columns treated as sums; GWL is averaged per year
FLUX_VARS = [
    "RAIN",
    "IRRIG",
    "INTERC",
    "RUNOFF",
    "EPOT",
    "EACT",
    "DRAINAGE",
    "QBOTTOM",
    "TPOT",
    "TACT",
    "DSTOR",
]
STATE_VARS = ["GWL"]


def load_fixture(path: Path):
    with path.open() as f:
        return json.load(f)


def aggregate(csv_path: Path):
    """Aggregate monthly CSV into annual stats (sum for flux, mean for GWL)."""
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
        yr = years.setdefault(year, {k: [] for k in FLUX_VARS + STATE_VARS})
        for k in FLUX_VARS:
            if k in rec and rec[k]:
                yr[k].append(float(rec[k]))
        for k in STATE_VARS:
            if k in rec and rec[k]:
                yr[k].append(float(rec[k]))

    annual = {}
    for year, vals in years.items():
        annual[year] = {}
        for k in FLUX_VARS:
            if vals[k]:
                annual[year][k] = round(sum(vals[k]), 2)
        for k in STATE_VARS:
            if vals[k]:
                annual[year][k] = round(sum(vals[k]) / len(vals[k]), 2)

    # totals/means across years
    totals = {}
    means = {}
    n_years = len(annual)
    for k in FLUX_VARS:
        totals[k] = round(sum(annual[y].get(k, 0.0) for y in annual), 2)
        means[k] = round(totals[k] / n_years, 2)
    for k in STATE_VARS:
        means[k] = round(sum(annual[y].get(k, 0.0) for y in annual) / n_years, 2)

    return annual, totals, means


def compare(expected, actual_years, actual_totals, actual_means):
    def check_block(block_name, exp_block, act_block):
        for year_or_key, exp_vals in exp_block.items():
            if isinstance(exp_vals, dict):
                act_vals = act_block.get(year_or_key, {})
                for var, exp_val in exp_vals.items():
                    act_val = act_vals.get(var)
                    if act_val is None or not math.isclose(act_val, exp_val, abs_tol=TOL):
                        raise AssertionError(
                            f"Mismatch {block_name}:{year_or_key}:{var}: expected {exp_val}, got {act_val}"
                        )
            else:
                act_val = act_block.get(year_or_key)
                if act_val is None or not math.isclose(act_val, exp_vals, abs_tol=TOL):
                    raise AssertionError(
                        f"Mismatch {block_name}:{year_or_key}: expected {exp_vals}, got {act_val}"
                    )

    check_block("years", expected["years"], actual_years)
    check_block("total", expected["total"], actual_totals)
    check_block("mean", expected["mean"], actual_means)


def main():
    if not SWAP_BIN.exists():
        raise SystemExit(f"swap binary not found at {SWAP_BIN}; build first (pixi run build-linux)")

    expected = load_fixture(FIXTURE)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        # copy case files to temp sandbox
        shutil.copytree(CASE_DIR, tmp / "case", dirs_exist_ok=True)
        workdir = tmp / "case"

        # ensure template is named swap.swp
        swap_file = workdir / "swap.swp"
        if not swap_file.exists():
            shutil.copy(workdir / "swap_linux.swp.template", swap_file)

        # run swap
        proc = subprocess.run([str(SWAP_BIN)], cwd=workdir, capture_output=True, text=True)
        if proc.returncode not in (0, 100):  # swap_main exits with 100 on success
            raise SystemExit(
                f"swap failed with code {proc.returncode}\nstdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
            )

        csv_path = workdir / "result_output.csv"
        if not csv_path.exists():
            raise SystemExit("result_output.csv not produced")

        annual, totals, means = aggregate(csv_path)
        compare(expected, annual, totals, means)

    print("✓ regression ok (hupselbrook annual stats match fixture)")


if __name__ == "__main__":
    main()
