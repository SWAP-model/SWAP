"""Byte-identical regression suite (pytest).

Each registered case runs the modern SWAP build in an isolated temp dir,
aggregates its ``result_output.csv`` into annual stats, and compares them
against the stored ``*_reference_gf.json`` fixture (produced from the
gfortran-compiled SWAP 4.2.0 oracle by ``regen_reference.py``). A case fails
if any aggregated value differs beyond ``TOL``.

Run it with pytest — there is no standalone CLI:

    pixi run -e test check-fast     # the `fast` subset, xdist-parallel
    pixi run -e test check-full     # every case
    pytest tests/regression -k hupselbrook        # one case
    pytest tests/regression -m fast               # the fast subset

Two policies are expressed as xfail markers, driven by the case registry in
``regression_harness.py``:

  - ``known_divergence`` — the modern build is knowingly off from 4.2.0; the
    case xfails (visible, non-fatal). If it starts matching, it xpasses — a
    signal to remove the flag.
  - ``pending_restore`` — the option's compute is gated/deleted, so the modern
    build fatal-errors (RuntimeError); the case xfails until restored.

The fixtures are the assertion of record; regenerate them only for a
documented physics change (``regen_reference.py``), never from the modern build.
"""

import pytest

from regression_harness import (
    CASES,
    FAST_CASES,
    SWAP_BIN,
    TESTS_DIR,
    compare,
    load_fixture,
    run_and_aggregate,
)


def _case_param(case):
    """Build a pytest param for one case, attaching fast/xfail markers."""
    marks = []
    if case.name in FAST_CASES:
        marks.append(pytest.mark.fast)
    if case.known_divergence:
        # Divergence surfaces as an assertion mismatch; xfail non-strict so an
        # unexpected match reports as xpass rather than failing the suite.
        marks.append(pytest.mark.xfail(
            reason=f"known divergence: {case.known_divergence}",
            strict=False, raises=AssertionError))
    if case.pending_restore:
        # The modern build is expected to fatal-error (RuntimeError) here.
        marks.append(pytest.mark.xfail(
            reason=f"pending restore: {case.pending_restore}",
            strict=False, raises=RuntimeError))
    return pytest.param(case, id=case.name, marks=marks)


@pytest.fixture(scope="session", autouse=True)
def _require_swap_binary():
    """Fail loudly (once) if the modern build is missing, rather than per case."""
    if not SWAP_BIN.exists():
        pytest.fail(
            f"swap binary not found at {SWAP_BIN}; build first "
            f"(pixi run build-linux)", pytrace=False)


@pytest.mark.parametrize("case", [_case_param(c) for c in CASES.values()])
def test_regression(case):
    """Modern output must match the 4.2.0 reference fixture within TOL."""
    fixture_path = TESTS_DIR / "regression" / case.fixture
    if not fixture_path.exists():
        pytest.fail(f"fixture not found at {fixture_path}", pytrace=False)

    expected = load_fixture(fixture_path)
    annual, totals, means = run_and_aggregate(case)  # RuntimeError => pending_restore xfail
    compare(expected, annual, totals, means)         # AssertionError => mismatch / known_divergence xfail
