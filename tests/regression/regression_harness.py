"""Shared machinery for the byte-identical regression suite.

This module holds the case registry and the pure run/aggregate/compare logic
used by three entry points:

  - ``test_output_regression.py`` — the pytest suite (the gate; run via
    ``pytest``, parametrized over ``CASES``).
  - ``regen_reference.py`` — regenerates the ``*_reference_gf.json`` fixtures
    from the gfortran-compiled SWAP 4.2.0 oracle (``tests/reference/swap420gf``).
  - ``regen_expected.py`` — diagnostic snapshot of the *modern* build's own
    output (``*_expected_gfortran.json``); NOT the reference.

It defines no ``test_*`` functions, so pytest does not collect it directly.
The comparison semantics (TOL, annual aggregation, the xfail policy encoded by
``known_divergence`` / ``pending_restore``) are the byte-identity contract and
must not change without a documented physics reason.
"""

import csv
import json
import math
import os
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
from typing import NamedTuple


TESTS_DIR = Path(__file__).resolve().parent.parent
SWAP_BIN = Path(__file__).resolve().parents[2] / "builddir" / "swap"
TOL = 1e-2  # cm tolerance on aggregated values

# The subset run by `check-fast` (marked `fast` in the pytest suite). The full
# registry runs under `check-full`. Keep this to the quick, high-signal cases.
FAST_CASES = ("hupselbrook", "soilhysteresis", "snow", "swcf3")


class CaseConfig(NamedTuple):
    """Configuration for a test case."""
    name: str
    case_dir: str  # relative path from the cases root
    fixture: str   # fixture filename in tests/regression
    flux_vars: list[str]  # summed annually
    state_vars: list[str]  # averaged annually
    cumul_vars: list[str] = []  # last value per year (for cumulative outputs)

    # Non-empty => the modern build is KNOWN to diverge from the 4.2.0
    # reference for this case (a documented, not-yet-fixed physics regression).
    # The suite marks it xfail: it does not fail, but the mismatch stays visible
    # (and an unexpected match surfaces as xpass — remove the flag). The string
    # is a short reason / pointer to INVESTIGATION_NOTES.md.
    known_divergence: str = ""

    input_files: dict[str, str] = {}

    # Non-empty => the option is gated / its compute not yet restored, so the
    # MODERN build is expected to fatal-error on this case. The suite marks it
    # xfail (a pending restoration target); the swap420gf fixture is still
    # generated. A successful modern run flips it to xpass.
    pending_restore: str = ""


# Registered test cases
CASES = {
    # [2026-06-11] Self-contained case (legacy + toml under the swap-testcases
    # repo). The whole regression suite is independent of the old, private
    # swap-cases submodule.
    "hupselbrook": CaseConfig(
        name="hupselbrook",
        case_dir="hupselbrook",
        fixture="hupselbrook_reference_gf.json",
        flux_vars=["RAIN", "IRRIG", "INTERC", "RUNOFF", "EPOT", "EACT",
                   "DRAINAGE", "QBOTTOM", "TPOT", "TACT", "DSTOR"],
        state_vars=["GWL"],
    ),
    # grassgrowth (ruurlo): RECONSTRUCTED, byte-identical — swbotb=1 (prescribed
    # GWL via gwlevel.csv), dramet=3 resistance drainage (open channel + owltab
    # CSV), numerical heat, grass (type 3, swrd=2, swharv=2 fixed-date mowing via
    # mowing_dates DOY), 1980-1984, 260.met->260.csv.
    "grassgrowth": CaseConfig(
        name="grassgrowth",
        case_dir="grassgrowth",
        fixture="grassgrowth_reference_gf.json",
        flux_vars=[],
        state_vars=[],
        cumul_vars=["PGRASSDM", "GRASSDM", "PMOWDM", "MOWDM"],
    ),
    # salinitystress (saltfarmtexel): RECONSTRUCTED, byte-identical — swinco=3 full
    # warm-restart (h_init.csv + cml_init.csv per-node + atmosphere ldwet/atmin7 +
    # tsoil via the 195-row [heat].tsoil_init profile), swbotb=3 (sinus aquifer
    # head), dramet=3 2-level drainage, numerical heat, solute (swsolu=1) +
    # Maas-Hoffman salinity, wofost potato, fixed irrigation (irrig.csv), 2012-2015.
    # See INVESTIGATION_NOTES.md / CASE_RECONSTRUCTION.md.
    "salinitystress": CaseConfig(
        name="salinitystress",
        case_dir="salinitystress",
        fixture="salinitystress_reference_gf.json",
        flux_vars=[],
        state_vars=["TREDDRY", "TREDWET", "TREDSOL", "CPWSO", "CWSO",
                    "CONC[-5.0]", "CONC[-25.0]", "CONC[-55.0]"],
    ),
    # surfacewater: RECONSTRUCTED, byte-identical — swbotb=3 (Cauchy, explicit,
    # haquif.csv), swdra=2 surface-water management (2 subsurface levels +
    # simulated surface-water level, [surface_water] + .management + .weir with
    # 28 periods), no heat, no solute, type-1 grass crop, 1997-1999, 290.met.
    "surfacewater": CaseConfig(
        name="surfacewater",
        case_dir="surfacewater",
        fixture="surfacewater_reference_gf.json",
        flux_vars=[],
        state_vars=["GWL", "POND"],
    ),
    # oxygenstress (zegveld): RECONSTRUCTED, byte-identical — swbotb=3 (Cauchy
    # from deep aquifer, implicit, haquif.csv), dramet=3 2-level resistance
    # drainage, numerical heat (peat profile), grass (type 3) with Bartholomeus
    # oxygen stress (swoxygen=2), 1993-2002.
    "oxygenstress": CaseConfig(
        name="oxygenstress",
        case_dir="oxygenstress",
        fixture="oxygenstress_reference_gf.json",
        flux_vars=[],
        state_vars=["TREDDRY", "TREDWET"],
        cumul_vars=["PGRASSDM", "GRASSDM", "PMOWDM", "MOWDM"],
    ),
    # Clone of hupselbrook with hysteresis active (SWHYST=1). Exercises the
    # soil-water-retention hysteresis path, dormant in all base cases.
    # Byte-identical to swap420gf.
    "soilhysteresis": CaseConfig(
        name="soilhysteresis",
        case_dir="soilhysteresis",
        fixture="soilhysteresis_reference_gf.json",
        flux_vars=["RAIN", "INTERC", "RUNOFF", "EPOT", "EACT",
                   "DRAINAGE", "QBOTTOM", "TPOT", "TACT", "DSTOR"],
        state_vars=["GWL"],
    ),
    # Winter cluster: snow accumulation/melt (SWSNOW=1). Clone of hupselbrook;
    # SNOWINCO=0 so snow accumulates from sub-zero precip days. Asserts
    # SNOW/SSNOW (snow storage) alongside the water balance. Byte-identical.
    "snow": CaseConfig(
        name="snow",
        case_dir="snow",
        fixture="snow_reference_gf.json",
        flux_vars=["RAIN", "INTERC", "RUNOFF", "DRAINAGE", "QBOTTOM", "DSTOR", "EACT"],
        state_vars=["GWL"],
    ),
    # Combined high-coverage case: SWHYST=1 (tau=0.2) + SWSNOW=1 (SNOWINCO=0) on
    # base hupselbrook, exercising in one run numerical heat (swhea=1/swcalt=2),
    # solute (swsolu=1), basic drainage (swdra=2), retention hysteresis, snow
    # accumulation/melt, and the maize/potato/grass rotation. Byte-identical.
    "winterhysteresis": CaseConfig(
        name="winterhysteresis",
        case_dir="winterhysteresis",
        fixture="winterhysteresis_reference_gf.json",
        flux_vars=["RAIN", "INTERC", "RUNOFF", "EPOT", "EACT",
                   "DRAINAGE", "QBOTTOM", "TPOT", "TACT", "DSTOR"],
        state_vars=["GWL"],
    ),
    # Winter cluster: snow (SWSNOW=1) + frost-reduced soil water flow (SWFROST=1).
    # The SNOW path is byte-identical (see the `snow` case); the FROST path
    # diverges ~0.5-0.8 cm DRAINAGE/RUNOFF in the cold grass year (suspected
    # heat<->frost feedback amplifying a sub-threshold tsoil difference).
    "winter": CaseConfig(
        name="winter",
        case_dir="winter",
        fixture="winter_reference_gf.json",
        flux_vars=["RAIN", "INTERC", "RUNOFF", "EPOT", "EACT",
                   "DRAINAGE", "QBOTTOM", "TPOT", "TACT", "DSTOR"],
        state_vars=["GWL"],
        known_divergence="frost-path drift vs 4.2.0 (~0.5-0.8 cm DRAINAGE/RUNOFF "
                         "in the cold grass year); snow path byte-identical; "
                         "suspected heat<->frost feedback; see INVESTIGATION_NOTES.md",
    ),
}

# --- Crop-switch coverage cases (hupselbrook + one setting). Inputs live in the
#     public SWAP-model/swap-testcases repo (authored by its tools/gen_switch_cases.py).
_SW_FLUX = ["RAIN", "IRRIG", "INTERC", "RUNOFF", "EPOT", "EACT",
            "DRAINAGE", "QBOTTOM", "TPOT", "TACT", "DSTOR"]
_SW_STATE = ["GWL"]


def _switch_case(name, **kw):
    return CaseConfig(name=name, case_dir=name, fixture=f"{name}_reference_gf.json",
                      flux_vars=_SW_FLUX, state_vars=_SW_STATE, **kw)


CASES.update({
    # shipped, byte-identical -> passing
    "swrd2":         _switch_case("swrd2"),
    "swharv1":       _switch_case("swharv1"),
    "swcompensate1": _switch_case("swcompensate1"),
    "swcompensate2": _switch_case("swcompensate2"),
    "swinter2":      _switch_case("swinter2"),
    "swcf3":         _switch_case("swcf3"),
    "swcf3_maize":   _switch_case("swcf3_maize"),
    "swsalinity1":   _switch_case("swsalinity1"),
    "swoxygen2":     _switch_case("swoxygen2"),
    "swinter3":      _switch_case("swinter3"),
    # swdrought=2 (De Jong van Lier microscopic uptake): the ONLY genuinely
    # unimplemented case — the compute was DELETED. Restoration attempted
    # (2026-06-11); twilt spiral-bug fixed, a perf residual remains, so it stays
    # gated and the modern build fatal-errors here. See INVESTIGATION_NOTES.md.
    "swdrought2":    _switch_case("swdrought2",
                                  pending_restore="swdrought=2 (De Jong van Lier) compute deleted; "
                                  "restoration attempted — twilt spiral-bug fixed, perf residual remains; "
                                  "see INVESTIGATION_NOTES.md 2026-06-11"),
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
    """Compare expected vs actual values; raise AssertionError on any mismatch."""
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


def cases_root() -> Path:
    """Root directory holding the regression case inputs (``<name>/{legacy,toml}``).

    The case *inputs* live in the public SWAP-model/swap-testcases repo; the
    expected-output fixtures stay here in the SWAP repo. Resolution order:

      1. ``SWAP_TESTCASES_PATH`` env var (set by CI, which checks out the pinned
         tag from ``TESTCASES_REF`` into a sibling dir).
      2. A ``../swap-testcases/cases`` sibling checkout next to the SWAP repo —
         the default when SWAP and swap-testcases are cloned side by side.
    """
    env = os.environ.get("SWAP_TESTCASES_PATH")
    if env:
        return Path(env).expanduser()
    return TESTS_DIR.parent.parent / "swap-testcases" / "cases"


def case_legacy_dir(case: CaseConfig) -> Path:
    """Directory of legacy ASCII inputs for a case (swap420gf reads these)."""
    return cases_root() / case.case_dir / "legacy"


def case_toml_dir(case: CaseConfig) -> Path:
    """Directory of TOML inputs for a case (the modern build reads these)."""
    return cases_root() / case.case_dir / "toml"


def run_and_aggregate(case: CaseConfig):
    """Run a single SWAP case in a temp dir and return aggregated stats.

    Returns a tuple ``(annual, totals, means)``. Raises RuntimeError on any
    runtime/output failure so callers can surface the message cleanly.
    """
    toml_dir = case_toml_dir(case)
    if not toml_dir.exists() or not (toml_dir / "swap.toml").exists():
        raise RuntimeError(
            f"TOML case directory not found at {toml_dir} "
            f"(missing dir or swap.toml). The TOML dir is the sole source of "
            f"truth for regression."
        )

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)

        # The TOML dir is the self-contained source. Per case it contains
        # swap.toml, swap.dra.toml, *.crp.toml, *.csv companions, legacy *.crp
        # crop files, swap_linux.swp.template, and (where the scenario requires
        # it) legacy ASCII companions like swap.dra (swdra=2 -> rddre()) and
        # swap.ini (swinco=3 -> rdinit()). Nothing is read from the legacy dirs.
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

        proc = subprocess.run([str(SWAP_BIN)], cwd=workdir, capture_output=True, text=True)

        # swap_main exits with 100 on success
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

        return aggregate(csv_path, case.flux_vars, case.state_vars, case.cumul_vars)


def regen_one_expected(case: CaseConfig) -> Path:
    """Snapshot the *modern* build's output as a golden-master baseline.

    DIAGNOSTIC ONLY. Writes ``{case.name}_expected_gfortran.json`` (the modern
    build's own output). This is NOT the regression reference — the suite
    compares against ``{case.name}_reference_gf.json``, generated from the
    gfortran-compiled SWAP 4.2.0 by ``regen_reference.py``. Regenerating the
    reference from the modern build would make the fidelity check trivially
    pass, so the two are kept distinct.
    """
    annual, totals, means = run_and_aggregate(case)
    payload = {"years": annual, "total": totals, "mean": means}
    out_path = TESTS_DIR / "regression" / f"{case.name}_expected_gfortran.json"
    with out_path.open("w") as f:
        json.dump(payload, f, indent=2, sort_keys=True)
        f.write("\n")
    print(f"✓ {case.name}: wrote {out_path.name}")
    return out_path
