"""Generate regression reference fixtures from the gfortran-compiled SWAP 4.2.0.

The reference binary `tests/reference/swap420gf` is the *unmodified* 4.2.0 source
compiled with the modern build's gfortran flags (see tests/reference/README.md). It
reads the **legacy ASCII** inputs in `tests/swap-cases/<N>.<case>/`. This script runs
it per case, aggregates `result_output.csv` with the same `aggregate()` and the same
flux/state/cumul variable sets the regression harness uses, and writes
`<case>_reference_gf.json` next to the other fixtures.

The regression harness (`test_output_regression.py`) compares the modern build's
output against these reference fixtures — so a divergence means the modern build's
*physics* differs from 4.2.0 (the compiler is held constant).

Usage:
    python3 tests/regression/regen_reference.py            # all registered cases
    python3 tests/regression/regen_reference.py hupselbrook # one case
"""

import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from test_output_regression import CASES, aggregate, TESTS_DIR

REF_BIN = TESTS_DIR / "reference" / "swap420gf"
# Committed legacy inputs per case; everything else in the dir is generated output.
LEGACY_INPUTS = ("*.swp.template", "*.crp", "*.dra", "*.met", "*.csv", "*.ini",
                 "*.bbc", "*.dat", "*.irg")


def _run_reference_and_aggregate(case):
    """Run swap420gf on a temp copy of the case's legacy ASCII dir, aggregate."""
    from test_output_regression import case_legacy_dir
    legacy_dir = case_legacy_dir(case)
    if not legacy_dir.exists():
        raise RuntimeError(f"legacy dir not found: {legacy_dir}")

    with tempfile.TemporaryDirectory() as tmpdir:
        work = Path(tmpdir) / "case"
        work.mkdir()
        # Copy only committed input files (avoid dragging in stale outputs).
        for pat in LEGACY_INPUTS:
            for f in legacy_dir.glob(pat):
                shutil.copy(f, work / f.name)
        # Stage swap.swp from the linux template (the binary reads swap.swp).
        template = work / "swap_linux.swp.template"
        if template.exists():
            shutil.copy(template, work / "swap.swp")

        proc = subprocess.run([str(REF_BIN)], cwd=work,
                              capture_output=True, text=True)
        # swap420 returns 100 on normal completion.
        if proc.returncode not in (0, 100):
            raise RuntimeError(
                f"swap420gf failed (rc={proc.returncode})\n"
                f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}")

        csv_path = work / "result_output.csv"
        if not csv_path.exists():
            raise RuntimeError("swap420gf produced no result_output.csv")

        return aggregate(csv_path, case.flux_vars, case.state_vars, case.cumul_vars)


def regen_one(case):
    annual, totals, means = _run_reference_and_aggregate(case)
    out = TESTS_DIR / "regression" / f"{case.name}_reference_gf.json"
    with out.open("w") as f:
        json.dump({"years": annual, "total": totals, "mean": means},
                  f, indent=2, sort_keys=True)
        f.write("\n")
    print(f"✓ {case.name}: wrote {out.name}")


def main():
    if not REF_BIN.exists():
        raise SystemExit(f"reference binary not found: {REF_BIN}")
    args = sys.argv[1:]
    selected = [CASES[a] for a in args] if args else list(CASES.values())
    for a in args:
        if a not in CASES:
            raise SystemExit(f"unknown case: {a}. known: {', '.join(CASES)}")
    failed = 0
    for case in selected:
        try:
            regen_one(case)
        except Exception as exc:
            print(f"✗ {case.name}: {exc}")
            failed += 1
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
