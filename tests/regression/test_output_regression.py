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

    # Non-empty => the modern build is KNOWN to diverge from the 4.2.0
    # reference for this case (a documented, not-yet-fixed physics regression).
    # The harness reports it as an expected-divergence (xfail): it does not
    # fail the suite, but prints the mismatch so the gap stays visible. The
    # string is a short reason / pointer to INVESTIGATION_NOTES.md.
    known_divergence: str = ""

    input_files: dict[str, str] = {}

    # local=True => case inputs live in the MAIN repo under
    # tests/regression/cases/<case_dir>/{legacy,toml}/ (not the swap-cases
    # submodule). Used for the generated crop-switch cases.
    local: bool = False

    # Non-empty => the option is gated / its compute not yet restored, so the
    # MODERN build is expected to fatal-error on this case. The harness treats
    # that runtime error as xfail (a pending restoration target); the swap420gf
    # fixture is still generated. A successful modern run flips it to xpass.
    pending_restore: str = ""


# Registered test cases
CASES = {
    "hupselbrook": CaseConfig(
        name="hupselbrook",
        case_dir="1.hupselbrook",
        fixture="hupselbrook_reference_gf.json",
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
    # [MACRO-RETIRE 2026-05-12] Case 3 (macroporeflow) retired per ADR 0040.
    # Macropore physics deleted from rescue branch; legacy SWAP 4.2.0
    # implementation preserved on branch legacy/swap-4.2.0. The case dir
    # still exists in the tests/swap-cases submodule for archival. A future
    # macropore feature arc will re-introduce this entry.
    "grassgrowth": CaseConfig(
        name="grassgrowth",
        case_dir="2.grassgrowth",
        fixture="grassgrowth_reference_gf.json",
        flux_vars=[],
        state_vars=[],
        cumul_vars=["PGRASSDM", "GRASSDM", "PMOWDM", "MOWDM"],
    ),
    "oxygenstress": CaseConfig(
        name="oxygenstress",
        case_dir="4.oxygenstress",
        fixture="oxygenstress_reference_gf.json",
        flux_vars=[],
        state_vars=["TREDDRY", "TREDWET"],
        cumul_vars=["PGRASSDM", "GRASSDM", "PMOWDM", "MOWDM"],
    ),
    "salinitystress": CaseConfig(
        name="salinitystress",
        case_dir="5.salinitystress",
        fixture="salinitystress_reference_gf.json",
        flux_vars=[],
        state_vars=["TREDDRY", "TREDWET", "TREDSOL", "CPWSO", "CWSO",
                    "CONC[-5.0]", "CONC[-25.0]", "CONC[-55.0]"],
    ),
    "surfacewater": CaseConfig(
        name="surfacewater",
        case_dir="6.surfacewater",
        fixture="surfacewater_reference_gf.json",
        flux_vars=[],
        state_vars=["GWL", "POND"],
    ),
    # Clone of hupselbrook with hysteresis active (SWHYST=1). Exercises the
    # soil-water-retention hysteresis path, dormant in all 6 base cases.
    # Hysteresis shifts GWL/DRAINAGE/storage vs the base, so the standard
    # water-balance columns assert the feature (see option-triage 2026-05-27).
    "soilhysteresis": CaseConfig(
        name="soilhysteresis",
        case_dir="7.soilhysteresis",
        fixture="soilhysteresis_reference_gf.json",
        flux_vars=["RAIN", "INTERC", "RUNOFF", "EPOT", "EACT",
                   "DRAINAGE", "QBOTTOM", "TPOT", "TACT", "DSTOR"],
        state_vars=["GWL"],
        # [FIX-ADAPTIVEDT 2026-06-11] Root cause (drainage first-step no-op) is
        # fixed; swinter3 + swcf3_maize (same desync) now match 4.2.0. This case
        # shares the identical root cause and is expected to xpass, but could not
        # be re-run locally (private tests/swap-cases submodule unavailable).
        # Re-run check-full with swap-cases present and REMOVE this flag once it
        # reports xpass. See INVESTIGATION_NOTES.md 2026-06-11.
        known_divergence="adaptive-dt threshold desync vs 4.2.0 (~1.9cm GWL, "
                         "amplified by hysteresis) — root cause FIXED 2026-06-11; "
                         "pending xpass re-confirmation with swap-cases; "
                         "see INVESTIGATION_NOTES.md",
    ),
    # Winter cluster: snow accumulation/melt (SWSNOW=1) + frost-reduced soil
    # water flow (SWFROST=1) + snow sublimation (SWSUBLIM=1). Clone of
    # hupselbrook; SNOWINCO=0 so snow accumulates from sub-zero precip days.
    # Asserts SNOW/SSNOW (snow storage) alongside the water balance.
    "winter": CaseConfig(
        name="winter",
        case_dir="8.winter",
        fixture="winter_reference_gf.json",
        flux_vars=["RAIN", "INTERC", "RUNOFF", "EPOT", "EACT",
                   "DRAINAGE", "QBOTTOM", "TPOT", "TACT", "DSTOR"],
        state_vars=["GWL", "SNOW"],  # SSNOW is identically 0 here — not asserted
        # Snow path reproduces 4.2.0 exactly (SNOW max 0.574 both); residual
        # ~0.04cm DRAINAGE/RUNOFF drift is the shared adaptive-dt desync (below).
        # [FIX-ADAPTIVEDT 2026-06-11] Same root cause as soilhysteresis/swinter3
        # (drainage first-step no-op) — now FIXED. Expected to xpass; could not be
        # re-run locally (private tests/swap-cases submodule unavailable). Re-run
        # check-full with swap-cases present and REMOVE this flag once it reports
        # xpass. See INVESTIGATION_NOTES.md 2026-06-11.
        known_divergence="adaptive-dt threshold desync vs 4.2.0 (~0.04cm "
                         "DRAINAGE/RUNOFF); snow path matches exactly — root cause "
                         "FIXED 2026-06-11; pending xpass re-confirmation with "
                         "swap-cases; see INVESTIGATION_NOTES.md",
    ),
}

# --- Crop-switch coverage cases (hupselbrook + one setting), generated by
#     gen_switch_cases.py, inputs under tests/regression/cases/<name>/.
_SW_FLUX = ["RAIN", "IRRIG", "INTERC", "RUNOFF", "EPOT", "EACT",
            "DRAINAGE", "QBOTTOM", "TPOT", "TACT", "DSTOR"]
_SW_STATE = ["GWL"]


def _switch_case(name, **kw):
    return CaseConfig(name=name, case_dir=name, fixture=f"{name}_reference_gf.json",
                      flux_vars=_SW_FLUX, state_vars=_SW_STATE, local=True, **kw)


CASES.update({
    # shipped, byte-identical -> passing
    "swrd2":         _switch_case("swrd2"),
    "swharv1":       _switch_case("swharv1"),
    "swcompensate1": _switch_case("swcompensate1"),
    "swcompensate2": _switch_case("swcompensate2"),
    "swinter2":      _switch_case("swinter2"),
    "swcf3":         _switch_case("swcf3"),
    # [FIX-ADAPTIVEDT 2026-06-11] Previously xfail ("0.01cm GWL FP artifact").
    # Root cause was the drainage first-step no-op (see INVESTIGATION_NOTES.md
    # 2026-06-11), not FP rounding: the wasted first step doubled dt prematurely
    # and the 0.01cm GWL was the residual desync. Now byte-identical to 4.2.0.
    "swcf3_maize":   _switch_case("swcf3_maize"),
    # gated / divergent -> pending restoration target (modern errors -> xfail)
    "swsalinity1":   _switch_case("swsalinity1",
                                  pending_restore="cropfixed swsalinity=1 gated "
                                  "(salinity->cml feedback divergence); see INVESTIGATION_NOTES.md"),
    "swoxygen2":     _switch_case("swoxygen2",
                                  pending_restore="wofost swoxygen=2 gated (type-2 "
                                  "oxygen-stress 0.01cm TACT divergence); see INVESTIGATION_NOTES.md"),
    # deleted/dormant compute, restoration targets (modern fatal-errors -> xfail)
    # [FIX-ADAPTIVEDT 2026-06-11] Previously xfail (adaptive-dt desync). The
    # drainage first-step no-op fix removed the premature dt-doubling, so the
    # adapted-Rutter run now matches 4.2.0 byte-for-byte under adaptive dt too.
    "swinter3":      _switch_case("swinter3"),
    "swdrought2":    _switch_case("swdrought2",
                                  pending_restore="swdrought=2 (De Jong van Lier) compute deleted; "
                                  "recover jongvanlier.f90 from 5c82f0a^"),
    # NOTE: swsalinity=2 (osmotic head) is intentionally NOT a case — swap420gf
    # itself SIGSEGVs on it (matricflux) in the hupselbrook config, so no oracle
    # fixture can be produced. See INVESTIGATION_NOTES.md.
})


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
                # total/mean blocks are flat {var: value} (no per-year nesting),
                # so year_or_key IS the variable name; there is no specific year.
                act_val = act_block.get(year_or_key)
                matches = act_val is not None and math.isclose(act_val, exp_vals, abs_tol=TOL)
                if not matches:
                    mismatches.append({
                        "block": block_name,
                        "year": "-",
                        "var": year_or_key,
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
        # Build comparison table. Columns: Block (years/total/mean), Year
        # (calendar year for the years block, "-" for run-level total/mean),
        # Variable, Expected, Actual, Diff.
        row_fmt = "  {:<6} {:>6} {:>13} {:>14} {:>14} {:>12}"
        header = row_fmt.format("Block", "Year", "Variable",
                                "Expected", "Actual", "Diff")
        lines = ["\n  Mismatches found (tolerance={:.0e}):".format(TOL),
                 header,
                 "  " + "-" * (len(header) - 2)]
        for m in mismatches:
            exp_str = f"{m['expected']:.4f}"
            act_str = f"{m['actual']:.4f}" if m['actual'] is not None else "None"
            diff_str = f"{m['diff']:.4f}" if m['diff'] is not None else "N/A"
            lines.append(row_fmt.format(
                m['block'], str(m['year']), m['var'], exp_str, act_str, diff_str))

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

def case_legacy_dir(case: CaseConfig) -> Path:
    """Directory of legacy ASCII inputs for a case (swap420gf reads these)."""
    if getattr(case, "local", False):
        return TESTS_DIR / "regression" / "cases" / case.case_dir / "legacy"
    return TESTS_DIR / "swap-cases" / case.case_dir


def case_toml_dir(case: CaseConfig) -> Path:
    """Directory of TOML inputs for a case (the modern build reads these)."""
    if getattr(case, "local", False):
        return TESTS_DIR / "regression" / "cases" / case.case_dir / "toml"
    return TESTS_DIR / "swap-cases" / "toml" / case.case_dir


def _run_and_aggregate(case: CaseConfig):
    """Run a single SWAP case in a temp dir and return aggregated stats.

    Returns a tuple ``(annual, totals, means)``. Raises RuntimeError on any
    runtime/output failure so callers can surface the message cleanly.
    """
    toml_dir = case_toml_dir(case)
    if not toml_dir.exists() or not (toml_dir / "swap.toml").exists():
        raise RuntimeError(
            f"TOML case directory not found at {toml_dir} "
            f"(missing dir or swap.toml). Phase 0 of CSV meteo finalization "
            f"made the TOML dir the sole source of truth for regression."
        )

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)

        # Phase 0: TOML dir is the self-contained source. Per case, it
        # contains swap.toml, swap.dra.toml, *.crp.toml, *.csv companions,
        # legacy *.crp crop files, swap_linux.swp.template, and (where the
        # scenario requires it) legacy ASCII companions like swap.dra
        # (swdra=2 → rddre()) and swap.ini (swinco=3 → rdinit()).
        # Everything SWAP needs at runtime lives here; nothing is read
        # from the legacy <N>.<case>/ dirs anymore.
        shutil.copytree(
            toml_dir,
            tmp / "case",
            ignore=shutil.ignore_patterns(
                'result_output.csv',
                'result_*.csv',
                '*.log',
                'output.*',
                '*.out',
            ),
        )
        workdir = tmp / "case"

        # Stage swap.swp from the in-dir template (legacy crop sub-readers
        # still call RDinit(swpfile)).
        swap_file = workdir / "swap.swp"
        if not swap_file.exists():
            template = workdir / "swap_linux.swp.template"
            if template.exists():
                shutil.copy(template, swap_file)

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


def run_case(case: CaseConfig) -> tuple[str, float]:
    """Run a single test case.

    Returns ``(status, execution_time)`` where ``status`` is one of:
      - ``"pass"``  — modern output matches the 4.2.0 reference.
      - ``"fail"``  — unexpected mismatch / error (fails the suite).
      - ``"xfail"`` — expected divergence (case.known_divergence set and it
                      diverged): reported, but does NOT fail the suite.
      - ``"xpass"`` — case.known_divergence set but the case now MATCHES; the
                      flag should be removed. Treated as a non-fatal warning.
    """
    start_time = time.perf_counter()

    fixture_path = TESTS_DIR / "regression" / case.fixture

    if not fixture_path.exists():
        print(f"✗ {case.name}: fixture not found at {fixture_path}")
        return "fail", time.perf_counter() - start_time

    expected = load_fixture(fixture_path)

    try:
        annual, totals, means = _run_and_aggregate(case)
    except RuntimeError as exc:
        elapsed = time.perf_counter() - start_time
        if getattr(case, "pending_restore", ""):
            print(f"⚠ {case.name}: PENDING RESTORE — {case.pending_restore} "
                  f"(modern build cannot run this yet) [xfail]")
            return "xfail", elapsed
        print(f"✗ {case.name}: {exc}")
        return "fail", elapsed

    try:
        compare(expected, annual, totals, means)
    except AssertionError as e:
        elapsed = time.perf_counter() - start_time
        if case.known_divergence:
            print(f"⚠ {case.name}: KNOWN DIVERGENCE — {case.known_divergence} "
                  f"[xfail, not a suite failure]{e}")
            return "xfail", elapsed
        print(f"✗ {case.name}: {e}")
        return "fail", elapsed

    elapsed = time.perf_counter() - start_time
    if case.known_divergence:
        print(f"⚠ {case.name}: now MATCHES 4.2.0 but known_divergence is still "
              f"set — remove the flag. [xpass] [{elapsed:.2f}s]")
        return "xpass", elapsed
    if getattr(case, "pending_restore", ""):
        print(f"⚠ {case.name}: modern build now RUNS this and matches 4.2.0 — "
              f"the option is restored; remove pending_restore. [xpass] [{elapsed:.2f}s]")
        return "xpass", elapsed
    print(f"✓ {case.name}: regression ok (annual stats match fixture) [{elapsed:.2f}s]")
    return "pass", elapsed


def _regen_one_case(case: CaseConfig) -> Path:
    """Snapshot the *modern* build's output as a golden-master baseline.

    DIAGNOSTIC ONLY. Writes ``{case.name}_expected_gfortran.json`` (the modern
    build's own output). This is NOT the regression reference — the harness
    compares against ``{case.name}_reference_gf.json``, generated from the
    gfortran-compiled SWAP 4.2.0 by ``regen_reference.py``. Regenerating the
    reference from the modern build would make the fidelity check trivially
    pass, so the two are kept distinct. Use this only to record where the
    modern build currently stands relative to 4.2.0.
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
        xfailed = 0
        xpassed = 0
        timings = []

        # Process results as they complete
        for future in concurrent.futures.as_completed(future_to_case):
            case = future_to_case[future]
            try:
                status, elapsed = future.result()
                timings.append((case.name, elapsed))
                if status == "pass":
                    passed += 1
                elif status == "xfail":
                    xfailed += 1
                elif status == "xpass":
                    xpassed += 1
                else:
                    failed += 1
            except Exception as exc:
                print(f'✗ {case.name} generated an exception: {exc}')
                failed += 1

    overall_elapsed = time.perf_counter() - overall_start

    # Print summary with timing information
    print(f"\n{'='*60}")
    summary = f"Results: {passed} passed, {failed} failed"
    if xfailed:
        summary += f", {xfailed} known-divergence (xfail)"
    if xpassed:
        summary += f", {xpassed} unexpectedly-matching (xpass — remove flag)"
    print(summary)
    print(f"Total execution time: {overall_elapsed:.2f}s")
    
    if timings:
        print(f"\nIndividual test timings:")
        for name, elapsed in sorted(timings, key=lambda x: x[1], reverse=True):
            print(f"  {name:20s} {elapsed:6.2f}s")
    
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
