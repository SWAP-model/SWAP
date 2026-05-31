"""SWAP-only XMI contract smoke test (no MODFLOW). Loads libswap_xmi via xmipy,
runs a few coupled-style steps, exercises get_value_ptr on the 3 exchange
arrays, asserts shapes and finite values."""
import sys, shutil
from pathlib import Path
import numpy as np
from xmipy import XmiWrapper

lib_path = Path(sys.argv[1])           # libswap_xmi.so
case_dir = Path(sys.argv[2])           # .../tests/coupling/hupselbrook_coupled
work     = Path(sys.argv[3])           # writable work dir (meson build dir)

work.mkdir(parents=True, exist_ok=True)
# The case is self-contained: stage every file in it (swap.toml + meteo + crop
# files + gwl csv + ensemble.txt). ensemble.txt ships =3 in the case dir.
for src in case_dir.iterdir():
    if src.is_file():
        shutil.copy(src, work / src.name)

swap = XmiWrapper(lib_path=str(lib_path), working_directory=str(work))
swap.initialize(str(work / "swap.toml"))

gwl          = swap.get_value_ptr("gwl")
qbot_volume  = swap.get_value_ptr("qbot_volume")
storage_coef = swap.get_value_ptr("storage_coef")
assert gwl.shape == (3,), gwl.shape
assert storage_coef.shape == (3,)
assert np.allclose(storage_coef, 0.15)

# Three coupled-style daily steps with distinct injected heads per column.
for day in range(3):
    gwl[:] = np.array([-0.5, -0.75, -1.0])      # metres
    swap.prepare_time_step(1.0)
    swap.prepare_solve(0)
    converged = swap.solve(0)
    swap.finalize_solve(0)
    swap.finalize_time_step()
    assert np.all(np.isfinite(qbot_volume)), qbot_volume

# Snapshot the recharge BEFORE finalize(): the get_value_ptr arrays are a
# zero-copy view onto the Fortran ensemble's allocatables, which finalize()
# deallocates -> reading the view afterwards is a use-after-free.
recharge = qbot_volume.copy()
swap.finalize()
print("SWAP XMI smoke OK:", recharge)
