"""Dev validator: prove a re-enabled cropfixed switch is oracle-correct.

For a given scenario it builds a maizes (type-1) crop variant that exercises one
previously-stub-errored switch, runs it through BOTH the legacy oracle
(tests/reference/swap420gf, classic ASCII input) and the modern build
(builddir/swap, TOML input), aggregates the hupselbrook flux/state vars the
regression harness uses, and reports any divergence beyond TOL.

This is a development aid (not a committed regression case). Usage:
    python tests/regression/_switch_validate.py swrd2
"""
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from test_output_regression import aggregate

TESTS = Path(__file__).resolve().parent.parent
CLASSIC = TESTS / "swap-cases" / "1.hupselbrook"
TOMLDIR = TESTS / "swap-cases" / "toml" / "1.hupselbrook"
REF_BIN = TESTS / "reference" / "swap420gf"
MODERN_BIN = Path(__file__).resolve().parents[2] / "builddir" / "swap"

FLUX = ["RAIN", "IRRIG", "INTERC", "RUNOFF", "EPOT", "EACT",
        "DRAINAGE", "QBOTTOM", "TPOT", "TACT", "DSTOR"]
STATE = ["GWL"]
TOL = 1e-2


def _patch(path: Path, subs):
    """Apply a list of (regex, replacement) substitutions to a text file."""
    txt = path.read_text()
    for pat, rep in subs:
        new = re.sub(pat, rep, txt, flags=re.MULTILINE)
        if new == txt:
            raise SystemExit(f"patch no-op in {path.name}: /{pat}/")
        txt = new
    path.write_text(txt)


def run_legacy(crp_subs, crp="maizes"):
    with tempfile.TemporaryDirectory() as d:
        w = Path(d)
        for pat in ("*.swp.template", "*.crp", "*.dra", "*.met", "*.csv", "*.ini"):
            for f in CLASSIC.glob(pat):
                shutil.copy(f, w / f.name)
        shutil.copy(w / "swap_linux.swp.template", w / "swap.swp")
        _patch(w / f"{crp}.crp", crp_subs)
        p = subprocess.run([str(REF_BIN)], cwd=w, capture_output=True, text=True)
        out = w / "result_output.csv"
        if not out.exists():
            raise SystemExit(f"LEGACY failed rc={p.returncode}\n{p.stdout}\n{p.stderr}")
        return aggregate(out, FLUX, STATE)


def run_modern(toml_subs, crp="maizes"):
    with tempfile.TemporaryDirectory() as d:
        w = Path(d) / "case"
        shutil.copytree(TOMLDIR, w, ignore=shutil.ignore_patterns('result_*.csv'))
        if not (w / "swap.swp").exists():
            shutil.copy(w / "swap_linux.swp.template", w / "swap.swp")
        _patch(w / f"{crp}.crp.toml", toml_subs)
        p = subprocess.run([str(MODERN_BIN)], cwd=w, capture_output=True, text=True)
        out = w / "result_output.csv"
        if p.returncode != 100 or not out.exists():
            raise SystemExit(f"MODERN failed rc={p.returncode}\n{p.stdout}\n{p.stderr}")
        return aggregate(out, FLUX, STATE)


# --- scenarios: (classic-crp patches, toml patches) ----------------------
# toml sections for stub-errored switches are absent in maizes.crp.toml, so we
# append them after the last line (cofab = 0.25).
_COFAB = r"^cofab\s*=\s*0\.25"
SCENARIOS = {
    # swrd=2: classic already has RDI/RRI/RDC/SWDMI2RD=1; toml needs swrd flip
    # + swdmi2rd=1 added (toml default is 0, classic maizes is 1).
    "swrd2": (
        [(r"^  SWRD = 1\b", "  SWRD = 2")],
        [(r"^swrd = 1\b", "swrd = 2\nswdmi2rd = 1")],
    ),
    # swharv=1: DVS-based harvest timing. dvsend=2.0 in both (default).
    "swharv1": (
        [(r"^  SWHARV = 0\b", "  SWHARV = 1")],
        [(r"^swharv = 0\b", "swharv = 1")],
    ),
    # swsalinity=1: Maas-Hoffman. hupselbrook has SWSOLU=1 so the branch is live.
    "swsal1": (
        [(r"^  SWSALINITY = 0\b", "  SWSALINITY = 1")],
        [(_COFAB, "cofab   = 0.25\n\n[salinity_stress]\nswsalinity = 1\n"
                  "saltmax = 3.0\nsaltslope = 0.1\n")],
    ),
    # swinter=2: Gash forest interception. Classic maizes already lists the
    # T/PFREE/PSTEM/SCANOPY/AVPREC/AVEVAP table; toml gains the 5 flat tables.
    "swinter2": (
        [(r"^  SWINTER = 1\b", "  SWINTER = 2")],
        [(r"^swinter = 1\b", "swinter = 2"),
         (_COFAB, "cofab   = 0.25\n"
                  "pfreetb   = [0.0, 0.9, 365.0, 0.9]\n"
                  "pstemtb   = [0.0, 0.05, 365.0, 0.05]\n"
                  "scanopytb = [0.0, 0.4, 365.0, 0.4]\n"
                  "avprectb  = [0.0, 6.0, 365.0, 6.0]\n"
                  "avevaptb  = [0.0, 1.5, 365.0, 1.5]\n")],
    ),
    # wofost swoxygen=2 (Bartholomeus physical): potatod already carries the full
    # param block, so just flip the switch.
    "wof_swoxygen2": (
        [(r"^  SWOXYGEN = 1\b", "  SWOXYGEN = 2")],
        [(r"(^swoxygen\s+= )1\b", r"\g<1>2")],
    ),
    # ===== WOFOST (potatod) restorations =====
    "wof_swharv1": (
        [(r"^  SWHARV = 0\b", "  SWHARV = 1")],
        [(r"^swharv = 0\b", "swharv = 1")],
    ),
    "wof_swcomp1": (
        [(r"^  SWCOMPENSATE = 0\b", "  SWCOMPENSATE = 1"),
         (r"^  SWSTRESSOR = 3\b", "  SWSTRESSOR = 2"),
         (r"^  ALPHACRIT = 1\.0\b", "  ALPHACRIT = 0.7")],
        [(r"^swcompensate = 0\b", "swcompensate = 1"),
         (r"^swstressor   = 3\b", "swstressor   = 2"),
         (r"^alphacrit    = 1\.0\b", "alphacrit    = 0.7")],
    ),
    "wof_swcomp2": (
        [(r"^  SWCOMPENSATE = 0\b", "  SWCOMPENSATE = 2"),
         (r"^  SWSTRESSOR = 3\b", "  SWSTRESSOR = 2")],
        [(r"^swcompensate = 0\b", "swcompensate = 2"),
         (r"^swstressor   = 3\b", "swstressor   = 2")],
    ),
    "wof_swinter2": (
        [(r"^  SWINTER = 1\b", "  SWINTER = 2")],
        [(r"^swinter = 1\b", "swinter = 2")],
    ),
    # ===== GRASS (grassd) restorations =====
    # grassd uses deprecated SWJARVIS=4; swap for explicit SWCOMPENSATE=2 Walsum.
    "grs_swcomp2": (
        [(r"^  SWJARVIS = 4\b",
          "  SWCOMPENSATE = 2\n  SWSTRESSOR = 2\n  DCRITRTZ = 16.0")],
        [(r"^swcompensate = 1\b", "swcompensate = 2"),
         (r"^swstressor   = 1\b", "swstressor   = 2")],
    ),
    "grs_swinter2": (
        [(r"^  SWINTER =  1\b", "  SWINTER =  2"),
         (r"(^  COFAB =  0\.25.*$)",
          r"\1" + "\n\n     T  PFREE  PSTEM  SCANOPY  AVPREC  AVEVAP\n"
          "   0.0    0.9   0.05      0.4     6.0     1.5\n"
          " 365.0    0.9   0.05      0.4     6.0     1.5\n* End of table")],
        [(r"^swinter = 1\b", "swinter = 2"),
         (r"(^cofab = 0\.25.*$)",
          r"\1" + "\ngashtb = [\n  [0.0, 0.9, 0.05, 0.4, 6.0, 1.5],\n"
          "  [365.0, 0.9, 0.05, 0.4, 6.0, 1.5],\n]")],
    ),
    # wofost swcf=3: LAI-indexed cf/cfeic/ch (replacing the swcf=2 DVS/CF/CH table).
    "wof_swcf3": (
        [(r"^  SWCF = 2\b", "  SWCF = 3"),
         (r" DVS   CF    CH\n 0\.0  1\.0   1\.0\n 1\.0  1\.1  40\.0\n"
          r" 2\.0  1\.1  50\.0\n\* End of table",
          " LAI   CF   CFEIC    CH\n 0.0  1.0    0.9    1.0\n"
          " 10.0 1.0    0.9   50.0\n* End of table")],
        [(r"^swcf   = 2\b", "swcf   = 3"),
         (r"cftb = \[\n  \[0\.0, 1\.0\],\n  \[1\.0, 1\.1\],\n  \[2\.0, 1\.1\],\n\]\n"
          r"chtb = \[\n  \[0\.0,  1\.0\],\n  \[1\.0, 40\.0\],\n  \[2\.0, 50\.0\],\n\]",
          "cftb = [\n  [0.0, 1.0],\n  [10.0, 1.0],\n]\n"
          "cfeictb = [\n  [0.0, 0.9],\n  [10.0, 0.9],\n]\n"
          "chtb = [\n  [0.0, 1.0],\n  [10.0, 50.0],\n]")],
    ),
    # grass swcf=3: LAI-indexed cf/cfeic/ch (replacing the swcf=2 DNR/CH table).
    "grs_swcf3": (
        [(r"^  SWCF = 2\b", "  SWCF = 3"),
         (r"    DNR       CH     CF\n    0\.0     12\.0    1\.0\n"
          r"  180\.0     12\.0    1\.0\n  366\.0     12\.0    1\.0\n\* End of table",
          "    LAI    CF   CFEIC     CH\n    0.0   1.0    0.9    12.0\n"
          "   10.0   1.0    0.9    12.0\n* End of table")],
        [(r"^swcf   = 2\b", "swcf   = 3"),
         (r"^chtb = \[0\.0, 12\.0, 180\.0, 12\.0, 366\.0, 12\.0\]",
          "cftb = [0.0, 1.0, 10.0, 1.0]\n"
          "cfeictb = [0.0, 0.9, 10.0, 0.9]\n"
          "chtb = [0.0, 12.0, 10.0, 12.0]")],
    ),
    # diagnostic: swsalinity=1 but saltslope=0 -> branch runs, alpsol always 1.
    "swsal1_noslope": (
        [(r"^  SWSALINITY = 0\b", "  SWSALINITY = 1"),
         (r"^  SALTSLOPE = 0\.1\b", "  SALTSLOPE = 0.0")],
        [(_COFAB, "cofab   = 0.25\n\n[salinity_stress]\nswsalinity = 1\n"
                  "saltmax = 3.0\nsaltslope = 0.0\n")],
    ),
    # swcompensate=1: Jarvis. alphacrit<1 + swstressor=2 (drought) engages the
    # compensation math (maize has summer drought stress).
    "swcomp1": (
        [(r"^  SWCOMPENSATE = 0\b", "  SWCOMPENSATE = 1"),
         (r"^  SWSTRESSOR = 3\b", "  SWSTRESSOR = 2"),
         (r"^  ALPHACRIT = 1\.0\b", "  ALPHACRIT = 0.7")],
        [(_COFAB, "cofab   = 0.25\n\n[compensation]\nswcompensate = 1\n"
                  "swstressor = 2\nalphacrit = 0.7\n")],
    ),
    # swcf=1: plain crop factor (already allowed) — isolation probe for the
    # crop-factor ET path vs legacy on maize.
    "swcf1": (
        [(r"^  SWCF = 2\b", "  SWCF = 1")],
        [(r"^swcf   = 2\b",
          "swcf   = 1\n"
          "cftb   = [0.0, 0.8, 0.3, 0.8, 0.5, 0.9, 0.7, 1.0, 1.0, 1.1, 1.4, 1.2, 2.0, 1.2]")],
    ),
    # swcf=3: wet-crop factor. Classic gains a CFW column; toml gains cftb (the
    # CF column maize already lists) + cfeictb (constant 0.9 wet factor).
    "swcf3": (
        [(r"^  SWCF = 2\b", "  SWCF = 3"),
         (r" DVS   CF     CH\n 0\.0  0\.8    1\.0\n 0\.3  0\.8   15\.0\n"
          r" 0\.5  0\.9   40\.0\n 0\.7  1\.0  140\.0\n 1\.0  1\.1  170\.0\n"
          r" 1\.4  1\.2  180\.0\n 2\.0  1\.2  175\.0",
          " DVS   CF     CH   CFW\n 0.0  0.8    1.0   0.9\n 0.3  0.8   15.0   0.9\n"
          " 0.5  0.9   40.0   0.9\n 0.7  1.0  140.0   0.9\n 1.0  1.1  170.0   0.9\n"
          " 1.4  1.2  180.0   0.9\n 2.0  1.2  175.0   0.9")],
        [(r"^swcf   = 2\b",
          "swcf   = 3\n"
          "cftb   = [0.0, 0.8, 0.3, 0.8, 0.5, 0.9, 0.7, 1.0, 1.0, 1.1, 1.4, 1.2, 2.0, 1.2]\n"
          "cfeictb = [0.0, 0.9, 2.0, 0.9]")],
    ),
    # swcompensate=2: Walsum. dcritrtz=5.0 (classic), swstressor=2 (drought).
    "swcomp2": (
        [(r"^  SWCOMPENSATE = 0\b", "  SWCOMPENSATE = 2"),
         (r"^  SWSTRESSOR = 3\b", "  SWSTRESSOR = 2")],
        [(_COFAB, "cofab   = 0.25\n\n[compensation]\nswcompensate = 2\n"
                  "swstressor = 2\ndcritrtz = 5.0\n")],
    ),
}


# scenario -> crop basename to patch (default maizes / type-1).
SCENARIO_CROP = {
    "wof_swharv1": "potatod", "wof_swcomp1": "potatod", "wof_swcomp2": "potatod",
    "wof_swinter2": "potatod",
    "grs_swcomp2": "grassd", "grs_swinter2": "grassd", "grs_swcf3": "grassd",
    "wof_swcf3": "potatod", "wof_swoxygen2": "potatod",
}


def main():
    if len(sys.argv) != 2 or sys.argv[1] not in SCENARIOS:
        raise SystemExit(f"usage: {sys.argv[0]} <{'|'.join(SCENARIOS)}>")
    crp = SCENARIO_CROP.get(sys.argv[1], "maizes")
    crp_subs, toml_subs = SCENARIOS[sys.argv[1]]
    leg_a, leg_t, leg_m = run_legacy(crp_subs, crp)
    mod_a, mod_t, mod_m = run_modern(toml_subs, crp)

    print(f"\n=== {sys.argv[1]}: legacy(swap420gf) vs modern ===")
    worst = 0.0
    for k in sorted(set(leg_t) | set(mod_t)):
        lv, mv = leg_t.get(k, float('nan')), mod_t.get(k, float('nan'))
        d = abs(lv - mv)
        worst = max(worst, d)
        flag = "  <-- DIFF" if d > TOL else ""
        print(f"  {k:10s} legacy={lv:14.6f} modern={mv:14.6f} |d|={d:.2e}{flag}")
    for k in sorted(set(leg_m) | set(mod_m)):
        lv, mv = leg_m.get(k, float('nan')), mod_m.get(k, float('nan'))
        d = abs(lv - mv)
        worst = max(worst, d)
        flag = "  <-- DIFF" if d > TOL else ""
        print(f"  {k:10s} legacy={lv:14.6f} modern={mv:14.6f} |d|={d:.2e}{flag}")
    print(f"worst |d| = {worst:.2e}  ->  {'PASS' if worst <= TOL else 'FAIL'}")
    sys.exit(0 if worst <= TOL else 1)


if __name__ == "__main__":
    main()
