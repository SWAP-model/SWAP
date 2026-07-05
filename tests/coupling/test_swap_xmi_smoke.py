"""SWAP-only XMI contract smoke (no MODFLOW).

Loads libswap via xmipy, runs a few coupled-style steps, exercises
get_value_ptr on the three exchange arrays, and asserts shapes + finite values.
"""
import shutil

import numpy as np
from xmipy import XmiWrapper


def test_swap_xmi_contract(swap_lib, coupled_case_dir, tmp_path):
    work = tmp_path / "run"
    work.mkdir()
    # The case is self-contained: stage every file (swap.toml + meteo + crop
    # files + gwl csv + ensemble.txt, which ships =3).
    for src in coupled_case_dir.iterdir():
        if src.is_file():
            shutil.copy(src, work / src.name)

    swap = XmiWrapper(lib_path=str(swap_lib), working_directory=str(work))
    swap.initialize(str(work / "swap.toml"))
    try:
        gwl = swap.get_value_ptr("gwl")
        qbot_volume = swap.get_value_ptr("qbot_volume")
        storage_coef = swap.get_value_ptr("storage_coef")
        assert gwl.shape == (3,), gwl.shape
        assert storage_coef.shape == (3,)
        assert np.allclose(storage_coef, 0.15)

        # Three coupled-style daily steps with injected heads per column.
        for _ in range(3):
            gwl[:] = np.array([-0.5, -0.75, -1.0])  # metres
            swap.prepare_time_step(1.0)
            swap.prepare_solve(0)
            swap.solve(0)
            swap.finalize_solve(0)
            swap.finalize_time_step()
            assert np.all(np.isfinite(qbot_volume)), qbot_volume

        # Snapshot recharge BEFORE finalize(): get_value_ptr arrays are a
        # zero-copy view onto the Fortran allocatables, which finalize()
        # deallocates (reading the view afterwards is a use-after-free).
        recharge = qbot_volume.copy()
        assert np.all(np.isfinite(recharge))
    finally:
        swap.finalize()
