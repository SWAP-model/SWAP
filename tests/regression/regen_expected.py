"""Diagnostic snapshot of the *modern* build's own regression output.

Writes ``<case>_expected_gfortran.json`` — the modern build's aggregated output
for each case. This is NOT the regression reference: the pytest suite compares
against ``<case>_reference_gf.json`` (produced from the gfortran-compiled SWAP
4.2.0 oracle by ``regen_reference.py``). Snapshotting the modern build as its
own reference would make the fidelity check trivially pass, so the two are kept
distinct. Use this only to record where the modern build currently stands
relative to 4.2.0.

Usage:
    python3 tests/regression/regen_expected.py            # all registered cases
    python3 tests/regression/regen_expected.py hupselbrook # one case
"""

import sys

from regression_harness import CASES, SWAP_BIN, regen_one_expected


def main():
    if not SWAP_BIN.exists():
        raise SystemExit(f"swap binary not found at {SWAP_BIN}; build first "
                         f"(pixi run build-linux)")
    args = sys.argv[1:]
    for a in args:
        if a not in CASES:
            raise SystemExit(f"unknown case: {a}. known: {', '.join(CASES)}")
    selected = [CASES[a] for a in args] if args else list(CASES.values())

    print(f"Snapshotting {len(selected)} case(s) as *_expected_gfortran.json ...\n")
    failed = 0
    for case in selected:
        try:
            regen_one_expected(case)
        except Exception as exc:
            print(f"✗ {case.name}: {exc}")
            failed += 1
    print(f"\nDone: {len(selected) - failed}/{len(selected)} written.")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
