"""Build a 1-layer x 1-row x N-col unconfined transient MODFLOW 6 model:
two CHD 'channels' at the ends, RCHA on the interior cells (coupling target),
K grounded to SWAP's ~12.5 cm/d soil (=0.125 m/d). Lengths in metres, time in days."""
from pathlib import Path
import flopy

def build(ws: Path, ncol: int = 10, nper: int = 1096, h0: float = -0.5, h1: float = -1.0):
    ws.mkdir(parents=True, exist_ok=True)
    sim = flopy.mf6.MFSimulation(sim_name="swapmf", sim_ws=str(ws), exe_name="mf6")
    tdis = flopy.mf6.ModflowTdis(sim, time_units="days",
                                 perioddata=[(1.0, 1, 1.0)] * nper, nper=nper)
    ims = flopy.mf6.ModflowIms(sim, complexity="SIMPLE", outer_maximum=50,
                               inner_maximum=100, linear_acceleration="BICGSTAB")
    gwf = flopy.mf6.ModflowGwf(sim, modelname="swapmf", newtonoptions="NEWTON",
                               save_flows=True)
    delr = 10.0
    dis = flopy.mf6.ModflowGwfdis(gwf, nlay=1, nrow=1, ncol=ncol,
                                  delr=delr, delc=10.0, top=0.0, botm=-10.0,
                                  length_units="meters")
    flopy.mf6.ModflowGwfic(gwf, strt=-0.75)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=0.125)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=1e-5, sy=0.15, transient={0: True})
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[
        [(0, 0, 0), h0], [(0, 0, ncol - 1), h1]])
    rch_cells = {0: [[(0, 0, j), 0.0] for j in range(1, ncol - 1)]}
    rch = flopy.mf6.ModflowGwfrch(gwf, stress_period_data=rch_cells, pname="RCHA",
                                  maxbound=ncol - 2)
    oc = flopy.mf6.ModflowGwfoc(gwf, head_filerecord="swapmf.hds",
                                budget_filerecord="swapmf.cbc",
                                saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")])
    sim.write_simulation()
    return sim

if __name__ == "__main__":
    import sys
    build(Path(sys.argv[1]), ncol=int(sys.argv[2]) if len(sys.argv) > 2 else 10,
          nper=int(sys.argv[3]) if len(sys.argv) > 3 else 1096)
    print("wrote MODFLOW 6 model")
