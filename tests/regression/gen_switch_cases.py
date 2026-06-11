"""Generate permanent regression cases for the crop-switch options.

Each case is hupselbrook with ONE setting (+ its required parameters) changed, in
BOTH legacy ASCII and TOML form, written under tests/regression/cases/<name>/
{legacy,toml}/. Fixtures are then produced from swap420gf by regen_reference.py
and the cases registered in test_output_regression.py.

The patch tuples mirror the validated scenarios in _switch_validate.py.

Usage:  python tests/regression/gen_switch_cases.py
"""
import re
import shutil
import sys
from pathlib import Path

TESTS = Path(__file__).resolve().parent.parent
# [2026-06-11] Source from the self-contained LOCAL hupselbrook case (no longer
# the swap-cases submodule, whose toml/ tree was lost from the remote).
HUP_LEGACY = TESTS / "regression" / "cases" / "hupselbrook" / "legacy"
HUP_TOML = TESTS / "regression" / "cases" / "hupselbrook" / "toml"
OUT = TESTS / "regression" / "cases"

_COFAB = r"^cofab\s+=\s+0\.25"


def _patch(path: Path, subs):
    txt = path.read_text()
    for pat, rep in subs:
        new = re.sub(pat, rep, txt, flags=re.MULTILINE)
        if new == txt:
            raise SystemExit(f"patch no-op in {path.name}: /{pat}/")
        txt = new
    path.write_text(txt)


# Each case: { 'crp': basename, 'legacy': [(pat,rep)...], 'toml': [(pat,rep)...] }
# or a list of such dicts to patch multiple crops in one case.
GASH_CRP = ("\n\n     T  PFREE  PSTEM  SCANOPY  AVPREC  AVEVAP\n"
            "   0.0    0.9   0.05      0.4     6.0     1.5\n"
            " 365.0    0.9   0.05      0.4     6.0     1.5\n* End of table")
GASH_TOML = ("cofab   = 0.25\n"
             "pfreetb   = [0.0, 0.9, 365.0, 0.9]\n"
             "pstemtb   = [0.0, 0.05, 365.0, 0.05]\n"
             "scanopytb = [0.0, 0.4, 365.0, 0.4]\n"
             "avprectb  = [0.0, 6.0, 365.0, 6.0]\n"
             "avevaptb  = [0.0, 1.5, 365.0, 1.5]\n")

CASES = {
    # ---- shipped, byte-identical -> passing regressions ----
    "swrd2": [dict(crp="maizes",
                   legacy=[(r"^  SWRD = 1\b", "  SWRD = 2")],
                   toml=[(r"^swrd = 1\b", "swrd = 2\nswdmi2rd = 1")])],
    "swharv1": [dict(crp="maizes",
                     legacy=[(r"^  SWHARV = 0\b", "  SWHARV = 1")],
                     toml=[(r"^swharv = 0\b", "swharv = 1")]),
                dict(crp="potatod",
                     legacy=[(r"^  SWHARV = 0\b", "  SWHARV = 1")],
                     toml=[(r"(^swharv\s+= )0\b", r"\g<1>1")])],
    "swcompensate1": [dict(crp="maizes",
                           legacy=[(r"^  SWCOMPENSATE = 0\b", "  SWCOMPENSATE = 1"),
                                   (r"^  SWSTRESSOR = 3\b", "  SWSTRESSOR = 2"),
                                   (r"^  ALPHACRIT = 1\.0\b", "  ALPHACRIT = 0.7")],
                           toml=[(_COFAB, "cofab   = 0.25\n\n[compensation]\nswcompensate = 1\n"
                                          "swstressor = 2\nalphacrit = 0.7\n")])],
    "swcompensate2": [dict(crp="maizes",
                           legacy=[(r"^  SWCOMPENSATE = 0\b", "  SWCOMPENSATE = 2"),
                                   (r"^  SWSTRESSOR = 3\b", "  SWSTRESSOR = 2")],
                           toml=[(_COFAB, "cofab   = 0.25\n\n[compensation]\nswcompensate = 2\n"
                                          "swstressor = 2\ndcritrtz = 5.0\n")]),
                      dict(crp="grassd",
                           legacy=[(r"^  SWJARVIS = 4\b",
                                    "  SWCOMPENSATE = 2\n  SWSTRESSOR = 2\n  DCRITRTZ = 16.0")],
                           toml=[(r"^swcompensate = 1\b", "swcompensate = 2"),
                                 (r"^swstressor   = 1\b", "swstressor   = 2")])],
    "swinter2": [dict(crp="maizes",
                      legacy=[(r"^  SWINTER = 1\b", "  SWINTER = 2")],
                      toml=[(r"^swinter = 1\b", "swinter = 2"), (_COFAB, GASH_TOML)]),
                 dict(crp="potatod",
                      legacy=[(r"^  SWINTER = 1\b", "  SWINTER = 2")],
                      toml=[(r"(^swinter\s+= )1\b", r"\g<1>2")])],
    # ---- gated, pending restoration -> registered xfail (pending_restore) ----
    "swsalinity1": [dict(crp="maizes",
                         legacy=[(r"^  SWSALINITY = 0\b", "  SWSALINITY = 1")],
                         toml=[(_COFAB, "cofab   = 0.25\n\n[salinity_stress]\nswsalinity = 1\n"
                                        "saltmax = 3.0\nsaltslope = 0.1\n")])],
    "swoxygen2": [dict(crp="potatod",
                       legacy=[(r"^  SWOXYGEN = 1\b", "  SWOXYGEN = 2")],
                       toml=[(r"(^swoxygen\s+= )1\b", r"\g<1>2")])],
    # ---- shipped swcf=3 (LAI-indexed wet crop factor) — wofost+grass, byte-identical ----
    # maize swcf=3 is split out (swcf3_maize) because it carries the documented
    # pre-existing crop-factor + Penman-Monteith 0.01cm GWL FP artifact.
    "swcf3": [
        dict(crp="potatod",
             legacy=[(r"^  SWCF = 2\b", "  SWCF = 3"),
                     (r" DVS   CF    CH\n 0\.0  1\.0   1\.0\n 1\.0  1\.1  40\.0\n"
                      r" 2\.0  1\.1  50\.0\n\* End of table",
                      " LAI   CF   CFEIC    CH\n 0.0  1.0    0.9    1.0\n"
                      " 10.0 1.0    0.9   50.0\n* End of table")],
             toml=[(r"^swcf   = 2\b", "swcf   = 3"),
                   (r"cftb = \[\n  \[0\.0, 1\.0\],\n  \[1\.0, 1\.1\],\n  \[2\.0, 1\.1\],\n\]\n"
                    r"chtb = \[\n  \[0\.0,  1\.0\],\n  \[1\.0, 40\.0\],\n  \[2\.0, 50\.0\],\n\]",
                    "cftb = [\n  [0.0, 1.0],\n  [10.0, 1.0],\n]\n"
                    "cfeictb = [\n  [0.0, 0.9],\n  [10.0, 0.9],\n]\n"
                    "chtb = [\n  [0.0, 1.0],\n  [10.0, 50.0],\n]")]),
        dict(crp="grassd",
             legacy=[(r"^  SWCF = 2\b", "  SWCF = 3"),
                     (r"    DNR       CH     CF\n    0\.0     12\.0    1\.0\n"
                      r"  180\.0     12\.0    1\.0\n  366\.0     12\.0    1\.0\n\* End of table",
                      "    LAI    CF   CFEIC     CH\n    0.0   1.0    0.9    12.0\n"
                      "   10.0   1.0    0.9    12.0\n* End of table")],
             toml=[(r"^swcf   = 2\b", "swcf   = 3"),
                   (r"(^chtb = \[0\.0, 12\.0, 180\.0, 12\.0, 366\.0, 12\.0\])",
                    "cftb = [0.0, 1.0, 10.0, 1.0]\n"
                    "cfeictb = [0.0, 0.9, 10.0, 0.9]\n"
                    "chtb = [0.0, 12.0, 10.0, 12.0]")]),
    ],
    "swcf3_maize": [
        dict(crp="maizes",
             legacy=[(r"^  SWCF = 2\b", "  SWCF = 3"),
                     (r" DVS   CF     CH\n 0\.0  0\.8    1\.0\n 0\.3  0\.8   15\.0\n"
                      r" 0\.5  0\.9   40\.0\n 0\.7  1\.0  140\.0\n 1\.0  1\.1  170\.0\n"
                      r" 1\.4  1\.2  180\.0\n 2\.0  1\.2  175\.0",
                      " DVS   CF     CH   CFW\n 0.0  0.8    1.0   0.9\n 0.3  0.8   15.0   0.9\n"
                      " 0.5  0.9   40.0   0.9\n 0.7  1.0  140.0   0.9\n 1.0  1.1  170.0   0.9\n"
                      " 1.4  1.2  180.0   0.9\n 2.0  1.2  175.0   0.9")],
             toml=[(r"^swcf   = 2\b",
                    "swcf   = 3\n"
                    "cftb   = [0.0, 0.8, 0.3, 0.8, 0.5, 0.9, 0.7, 1.0, 1.0, 1.1, 1.4, 1.2, 2.0, 1.2]\n"
                    "cfeictb = [0.0, 0.9, 2.0, 0.9]")]),
    ],
}

# JvL drought-stress block (swdrought=2), values from legacy grass.crp.
_JVL = ("\n  WILTPOINT = -20000.0\n  KSTEM = 1.03d-4\n  RXYLEM = 0.02\n"
        "  ROOTRADIUS = 0.05\n  KROOT = 3.5d-5\n  ROOTCOEFA = 0.53\n"
        "  SWHYDRLIFT = 0\n  ROOTEFF = 1.0\n  STEPHR = 1.0\n"
        "  CRITERHR = 0.001\n  TACCUR = 0.001")

CASES.update({
    # ---- not-yet-restored compute (deleted/dormant) -> pending restoration ----
    "swinter3": [dict(crp="maizes",
                      legacy=[(r"^  SWINTER = 1\b", "  SWINTER = 3"),
                              (r"(^  COFAB = 0\.25\b.*$)",
                               r"\1" + "\n  FIMIN = 0.04\n  SICCAPLAI = 0.0042")],
                      toml=[(r"^swinter = 1\b", "swinter = 3"),
                            (_COFAB, "cofab   = 0.25\nfimin = 0.04\nsiccaplai = 0.0042\n")])],
    "swdrought2": [dict(crp="maizes",
                        legacy=[(r"^  SWDROUGHT = 1\b", "  SWDROUGHT = 2"),
                                (r"(^  ADCRL = 0\.1\b.*$)", r"\1" + _JVL)],
                        toml=[(r"^swdrought = 1\b", "swdrought = 2")])],
    # swsalinity=2 (osmotic head) is intentionally omitted: swap420gf itself
    # SIGSEGVs on it (matricflux) in this config, so no oracle fixture exists.
})

# Cases whose modern run currently fatal-errors (option gated / compute not yet
# restored). Registered with pending_restore so the suite treats the modern error
# as xfail; the swap420gf fixture is still generated as a restoration target.
PENDING = {"swsalinity1", "swoxygen2", "swinter3", "swdrought2"}


def gen_one(name, specs):
    legacy_dir = OUT / name / "legacy"
    toml_dir = OUT / name / "toml"
    for d in (legacy_dir, toml_dir):
        if d.exists():
            shutil.rmtree(d)
        d.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(HUP_LEGACY, legacy_dir)
    shutil.copytree(HUP_TOML, toml_dir,
                    ignore=shutil.ignore_patterns('result_output.csv', 'Swap.ok'))
    for spec in specs:
        _patch(legacy_dir / f"{spec['crp']}.crp", spec["legacy"])
        _patch(toml_dir / f"{spec['crp']}.crp.toml", spec["toml"])
    print(f"  generated {name}: {[s['crp'] for s in specs]}")


def main():
    sel = sys.argv[1:] or list(CASES)
    OUT.mkdir(exist_ok=True)
    for name in sel:
        gen_one(name, CASES[name])
    print(f"generated {len(sel)} case(s) under {OUT}")


if __name__ == "__main__":
    main()
