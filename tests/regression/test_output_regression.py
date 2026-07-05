"""Byte-identical regression suite (pytest, live double-run).

Each registered case runs **two** engines in isolated temp dirs and compares
their aggregated annual stats within ``TOL``:

  - the modern SWAP build on the TOML inputs, and
  - the selected *reference* engine on its input variant (default ``swap420gf``,
    the gfortran-compiled SWAP 4.2.0 oracle, reading the legacy ASCII inputs).

There are no stored expected fixtures — the reference side is regenerated every
run. Pick the reference with ``SWAP_REGRESSION_REF`` (see the ``REFERENCES``
registry in ``regression_harness.py``); the machinery is ready for a future
released SWAP to serve as the reference too.

Run it with pytest:

    pixi run -e test check-fast     # the `fast` subset, xdist-parallel
    pixi run -e test check-full     # every case
    pytest tests/regression -k hupselbrook        # one case
    pytest tests/regression -m fast               # the fast subset

Two policies are expressed as xfail markers from the case registry:

  - ``known_divergence`` — the modern build is knowingly off from the reference;
    the case xfails (visible, non-fatal). If it starts matching, it xpasses — a
    signal to remove the flag.
  - ``pending_restore`` — the option's compute is gated/deleted, so the modern
    build fatal-errors (RuntimeError); the case xfails until restored.
"""

import pytest

from regression_harness import (
    CASES,
    FAST_CASES,
    SWAP_BIN,
    active_reference,
    compare,
    run_and_aggregate,
    run_reference_and_aggregate,
)

# The reference is selected once (SWAP_REGRESSION_REF), shared across cases.
REFERENCE = active_reference()


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
def _require_binaries():
    """Fail loudly (once) if either engine is missing, rather than per case."""
    if not SWAP_BIN.exists():
        pytest.fail(
            f"modern swap binary not found at {SWAP_BIN}; build first "
            f"(pixi run build-linux)", pytrace=False)
    if not REFERENCE.binary.exists():
        pytest.fail(
            f"reference '{REFERENCE.name}' binary not found at "
            f"{REFERENCE.binary}", pytrace=False)


@pytest.mark.parametrize("case", [_case_param(c) for c in CASES.values()])
def test_regression(case):
    """Modern output must match the reference engine's output within TOL."""
    # Modern first: a pending_restore case fatal-errors here (RuntimeError) and
    # xfails before a reference run is spent.
    m_annual, m_totals, m_means = run_and_aggregate(case)
    # Live "expected" side.
    r_annual, r_totals, r_means = run_reference_and_aggregate(case, REFERENCE)
    expected = {"years": r_annual, "total": r_totals, "mean": r_means}
    compare(expected, m_annual, m_totals, m_means)  # AssertionError => known_divergence xfail
