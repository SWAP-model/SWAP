"""Generate system-level regression cases (hydrology/meteo switches).

Companion to gen_switch_cases.py, which patches CROP files. This one patches the
top-level swap.swp / swap.toml (and could patch swap.dra) to flip system switches
like hysteresis (SWHYST), snow (SWSNOW), and frost (SWFROST). Each case is base
hupselbrook with ONE feature group enabled, in BOTH legacy ASCII and TOML form,
written under tests/regression/cases/<name>/{legacy,toml}/. Fixtures come from
swap420gf via regen_reference.py; cases are registered in test_output_regression.py.

Why these exist: after the swap-cases toml/ tree was lost from the remote, the
former cases 7.soilhysteresis / 8.winter could no longer be regenerated from the
submodule. These reconstructions are self-contained and reproducible from the
local hupselbrook base.

KEY GOTCHA (byte-identity): when SWHYST=1 the legacy .swp carries `TAU = 0.2`,
but the modern TOML reader DEFAULTS `tau` to 0.0 when absent — a degenerate value
that flips the wetting/drying scanning curve every step. The TOML MUST set
`tau = 0.2` explicitly. (Worth a modern validator check: require tau>0 when
swhyst!=0.)

Usage:  python tests/regression/gen_system_cases.py [name ...]
"""
import re
import shutil
import sys
from pathlib import Path

TESTS = Path(__file__).resolve().parent.parent
HUP = TESTS / "regression" / "cases" / "hupselbrook"   # local self-contained base
OUT = TESTS / "regression" / "cases"

# A patch = (filename_in_dir, [(regex, replacement), ...]). The legacy swp lives
# in swap_linux.swp.template; the modern config in swap.toml.
SWP = "swap_linux.swp.template"
TOML = "swap.toml"

# Snow block appended to the modern TOML (inserted before [soil]).
_SNOW_TOML = ("\n[meteorology.snow]\n"
              "swsnow   = 1\nsnowcoef = 0.3\nteprrain = 2.0\nteprsnow = -2.0\n")
# Frost block appended to the modern TOML (inserted before [crop]).
_FROST_TOML = ("\n[soil.frost]\n"
               "swfrost   = 1\ntfroststa = 0.0\ntfrostend = -1.0\nswsublim  = 0\n")


def _ins_before(text, anchor, block):
    if anchor not in text:
        raise SystemExit(f"anchor not found: {anchor!r}")
    return text.replace(anchor, block + anchor, 1)


CASES = {
    # SWHYST=1 retention hysteresis. tau=0.2 mandatory in TOML (see module note).
    "soilhysteresis": {
        "legacy": {SWP: [(r"^  SWHYST = 0\b", "  SWHYST = 1")]},
        "toml":   {TOML: [(r"(?m)^swhyst  = 0$", "swhyst  = 1\ntau     = 0.2")]},
    },
    # SWSNOW=1 snow accumulation/melt, SNOWINCO=0 (build snow from sub-zero precip).
    "snow": {
        "legacy": {SWP: [(r"^  SWSNOW = 0\b", "  SWSNOW = 1"),
                         (r"^  SNOWINCO = 22\.0\b", "  SNOWINCO = 0.0")]},
        "toml":   {TOML: [("INSERT_SNOW", None)]},   # handled specially below
    },
    # SWSNOW=1 + SWFROST=1 winter cluster. Frost path is a KNOWN_DIVERGENCE
    # (heat<->frost feedback); kept for coverage + investigation.
    "winter": {
        "legacy": {SWP: [(r"^  SWSNOW = 0\b", "  SWSNOW = 1"),
                         (r"^  SWFROST = 0\b", "  SWFROST = 1"),
                         (r"^  SNOWINCO = 22\.0\b", "  SNOWINCO = 0.0")]},
        "toml":   {TOML: [("INSERT_SNOW", None), ("INSERT_FROST", None)]},
    },
    # Combined high-coverage: hysteresis + snow (both byte-identical individually).
    "winterhysteresis": {
        "legacy": {SWP: [(r"^  SWHYST = 0\b", "  SWHYST = 1"),
                         (r"^  SWSNOW = 0\b", "  SWSNOW = 1"),
                         (r"^  SNOWINCO = 22\.0\b", "  SNOWINCO = 0.0")]},
        "toml":   {TOML: [(r"(?m)^swhyst  = 0$", "swhyst  = 1\ntau     = 0.2"),
                          ("INSERT_SNOW", None)]},
    },
}


def _apply(path: Path, subs):
    txt = path.read_text()
    for pat, rep in subs:
        if pat == "INSERT_SNOW":
            txt = _ins_before(txt, "\n[soil]\n", _SNOW_TOML)
        elif pat == "INSERT_FROST":
            txt = _ins_before(txt, "\n[crop]\n", _FROST_TOML)
        else:
            new = re.sub(pat, rep, txt, flags=re.MULTILINE)
            if new == txt:
                raise SystemExit(f"patch no-op in {path.name}: /{pat}/")
            txt = new
    path.write_text(txt)


def gen_one(name, spec):
    dst = OUT / name
    if dst.exists():
        shutil.rmtree(dst)
    for side in ("legacy", "toml"):
        shutil.copytree(HUP / side, dst / side,
                        ignore=shutil.ignore_patterns("result_output.csv", "Swap.ok", "observed"))
        for fname, subs in spec.get(side, {}).items():
            _apply(dst / side / fname, subs)
    print(f"  generated {name}")


def main():
    sel = sys.argv[1:] or list(CASES)
    for name in sel:
        if name not in CASES:
            raise SystemExit(f"unknown system case: {name}. known: {', '.join(CASES)}")
        gen_one(name, CASES[name])
    print(f"generated {len(sel)} system case(s); now run regen_reference.py for fixtures")


if __name__ == "__main__":
    main()
