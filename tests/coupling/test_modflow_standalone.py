import sys, tempfile
from pathlib import Path
import flopy
from build_modflow import build

ws = Path(tempfile.mkdtemp()) / "mf"
sim = build(ws, ncol=10, nper=3)
ret, _ = sim.run_simulation(silent=True)
assert ret, "MODFLOW standalone run failed"
hds = flopy.utils.HeadFile(ws / "swapmf.hds").get_data().squeeze()
assert hds[0] != hds[-1], "expected a lateral gradient between the two channels"
print("MODFLOW standalone OK, heads:", hds)
