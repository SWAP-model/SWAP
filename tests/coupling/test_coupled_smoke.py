"""Coupled SWAP <-> MODFLOW 6 smoke (qualitative).

A short coupled run (ncol=10, nper=30) must complete without error, yield
finite final heads, show a lateral gradient between the two channels, write the
water-table plot, and show recharge actually moving the interior head (the
interior differs from the no-recharge linear channel-to-channel profile).

The definition of done is qualitative, not a strict numeric assertion. Requires
libmf6.so (the libmf6 fixture skips when it is absent).
"""
from __future__ import annotations

import numpy as np

from run_coupled import main


def test_coupled_smoke(swap_lib, libmf6, coupled_case_dir, tmp_path):
    run = tmp_path / "run"
    main(swap_lib, libmf6, coupled_case_dir, run, ncol=10, nper=30)

    hds = np.load(run / "heads.npy")
    final = hds[-1]
    assert np.all(np.isfinite(final)), f"non-finite final heads: {final}"
    assert final[0] != final[-1], f"no lateral gradient: {final}"
    assert (run / "coupled_gwl.png").exists(), "plot not written"
    linear = np.linspace(final[0], final[-1], final.size)
    assert not np.allclose(final[1:-1], linear[1:-1], atol=1e-6), \
        "interior heads match the no-recharge linear profile (no coupling effect)"
