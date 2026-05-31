"""Coupled SWAP <-> MODFLOW 6 smoke test (qualitative).

A short coupled run (ncol=10, nper=30) must:
  * complete without error,
  * yield finite final heads,
  * show a lateral gradient (heads[0] != heads[-1]) between the two channels,
  * write the final-water-table plot.

The definition of done is qualitative (sensible, nonzero exchange + a
recharge-shaped water table), not a strict numeric assertion.

Runs with zero args (pixi run -e test test-coupling); positional args override
the defaults (used by the meson registration / manual invocation):

    python test_coupled_smoke.py [LIBSWAP_XMI [LIBMF6 [CASE_DIR]]]

Or as a meson test (libmf6 path discovered at configure time): `coupled-smoke`.
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _defaults import default_case_dir, default_libmf6, default_libswap  # noqa: E402
from run_coupled import main  # noqa: E402


def run_smoke(libswap: Path, libmf6: Path, case_dir: Path) -> None:
    run = Path(tempfile.mkdtemp(prefix="swapmf_smoke_")) / "run"
    main(libswap, libmf6, case_dir, run, ncol=10, nper=30)

    hds = np.load(run / "heads.npy")
    final = hds[-1]
    assert np.all(np.isfinite(final)), f"non-finite final heads: {final}"
    assert final[0] != final[-1], f"no lateral gradient: {final}"
    assert (run / "coupled_gwl.png").exists(), "plot not written"
    # Exchange sanity: the interior water table must differ from the linear
    # channel-to-channel interpolation (i.e. recharge actually moved the head).
    linear = np.linspace(final[0], final[-1], final.size)
    assert not np.allclose(final[1:-1], linear[1:-1], atol=1e-6), \
        "interior heads match the no-recharge linear profile (no coupling effect)"
    print("coupled smoke OK; final water table:", np.array2string(final, precision=3))


if __name__ == "__main__":
    libswap = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else default_libswap()
    libmf6 = Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else default_libmf6()
    case_dir = Path(sys.argv[3]).resolve() if len(sys.argv) > 3 else default_case_dir()
    run_smoke(libswap, libmf6, case_dir)
