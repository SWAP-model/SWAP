"""MODFLOW 6 standalone sanity check (no SWAP coupling).

Builds a small 10-column model via flopy and asserts a lateral gradient between
the two channels. Proves MODFLOW 6 itself runs in the environment before the
coupled smoke attempts two-way exchange. Skips if the mf6 executable is absent.
"""
import shutil

import pytest

flopy = pytest.importorskip("flopy")

from build_modflow import build


def test_modflow_standalone(tmp_path):
    if shutil.which("mf6") is None:
        pytest.skip("mf6 executable not on PATH (ships in the pixi test env)")

    ws = tmp_path / "mf"
    sim = build(ws, ncol=10, nper=3)
    ret, _ = sim.run_simulation(silent=True)
    assert ret, "MODFLOW standalone run failed"
    hds = flopy.utils.HeadFile(ws / "swapmf.hds").get_data().squeeze()
    assert hds[0] != hds[-1], "expected a lateral gradient between the two channels"
