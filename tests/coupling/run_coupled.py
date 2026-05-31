"""Orchestrate a coupled SWAP <-> MODFLOW 6 run.

Build a two-channel MODFLOW 6 cross-section (flopy), stage the SWAP coupled
work dir, run the loop via the vendored ``SwapMod`` driver (MODFLOW leads the
clock, SWAP solves one column-ensemble step per day), then collect MODFLOW
heads and plot the final water table.

Usage:
    python run_coupled.py LIBSWAP_XMI LIBMF6 CASE_DIR RUN_DIR [NCOL] [NPER]

  LIBSWAP_XMI : path to libswap_xmi.so  (the XMI kernel; NOT libswap_bmi.so)
  LIBMF6      : path to libmf6.so
  CASE_DIR    : .../tests/coupling/hupselbrook_coupled (self-contained source
                of the SWAP coupled config + meteo + crop files + gwl csv)
  RUN_DIR     : writable run directory (created)
  NCOL        : MODFLOW columns (default 10); interior recharge cells = NCOL-2
  NPER        : daily stress periods (default 1096 = full 2002-2004)

Units / sign conventions (see README.md):
  * gwl(i)          MODFLOW head, metres, datum = surface = MODFLOW top = 0.
  * qbot_volume(i)  SWAP recharge DEPTH over the step, metres, positive down.
                    The driver divides by delt -> m/d rate for the RCHA package
                    (MODFLOW RCHA `recharge` is a flux rate L/T; MF6 multiplies
                    by cell area internally, so no area factor here).
  * storage_coef(i) specific yield (-), fixed 0.15 placeholder for the smoke.
"""
from __future__ import annotations

import os
import shutil
import sys
import tomllib
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent.resolve()
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "imod_coupler_fork"))

from build_modflow import build  # noqa: E402
from config import BaseConfig, SwapModConfig  # noqa: E402  (vendored fork)
from swapmod import SwapMod  # noqa: E402  (vendored fork)

SWAP_STAGE_FILES = [
    "swap.toml",
    "swap.gwl.csv",
    "283.csv",
    "grassd.crp.toml",
    "maizes.crp.toml",
    "potatod.crp.toml",
]


def stage_run_dir(libswap: Path, libmf6: Path, case_dir: Path, run_dir: Path,
                  ncol: int, nper: int) -> tuple[Path, int]:
    """Build the MODFLOW model + SWAP work dir + coupler config. Returns
    (config_dir, n_couple)."""
    run_dir.mkdir(parents=True, exist_ok=True)
    n_couple = ncol - 2  # interior recharge cells == SWAP column count

    # 1. MODFLOW model in run_dir/mf
    build(run_dir / "mf", ncol=ncol, nper=nper)

    # 2. SWAP work dir in run_dir/swap
    swap_ws = run_dir / "swap"
    swap_ws.mkdir(exist_ok=True)
    for f in SWAP_STAGE_FILES:
        src = case_dir / f
        if src.exists():
            shutil.copy(src, swap_ws / f)
    # ensemble.txt: SWAP column count == MODFLOW interior recharge cells.
    (swap_ws / "ensemble.txt").write_text(f"{n_couple}\n")

    # 3. Coupler config with discovered absolute dll paths.
    cfg_text = (HERE / "imod_coupler.toml").read_text() \
        .replace("<ABS_PATH_TO>/libmf6.so", str(libmf6)) \
        .replace("<ABS_PATH_TO>/libswap_xmi.so", str(libswap))
    (run_dir / "imod_coupler.toml").write_text(cfg_text)
    return run_dir, n_couple


def run_coupled(run_dir: Path, sanity_steps: int = 5) -> None:
    """Run the coupled loop from run_dir (chdir'd by SwapModConfig)."""
    cfg_dir = run_dir.resolve()
    data = tomllib.loads((cfg_dir / "imod_coupler.toml").read_text())
    base = BaseConfig(log_level="INFO", timing=False)
    smcfg = SwapModConfig(config_dir=cfg_dir, **data)  # chdirs to cfg_dir

    driver = SwapMod(base, smcfg)
    driver.initialize()
    print(f"[run_coupled] initialized; max_iter={int(driver.max_iter)} "
          f"delt(start) n_couple={driver.swap_volume.shape[0]} "
          f"mf6_head.shape={driver.mf6_head.shape} "
          f"mf6_recharge.shape={driver.mf6_recharge.shape}")

    step = 0
    while driver.get_current_time() < driver.get_end_time():
        driver.update()
        step += 1
        if step <= sanity_steps:
            heads_int = driver.mf6_head[driver._rch_idx()]
            print(
                f"[step {step:4d}] "
                f"qbot_volume(m)={_fmt(driver.swap_volume)}  "
                f"recharge(m/d)={_fmt(driver.mf6_recharge)}  "
                f"interior_head(m)={_fmt(heads_int)}  "
                f"chd=({driver.mf6_head[0]:.3f},{driver.mf6_head[-1]:.3f})"
            )
    print(f"[run_coupled] loop done after {step} steps")
    driver.finalize()


def _fmt(a: np.ndarray) -> str:
    a = np.asarray(a)
    return "[" + " ".join(f"{x:+.4e}" for x in a[: min(4, a.size)]) + \
        ("...]" if a.size > 4 else "]")


def collect_and_plot(run_dir: Path) -> np.ndarray:
    import flopy
    hds = flopy.utils.HeadFile(run_dir / "mf" / "swapmf.hds").get_alldata()
    hds = np.asarray(hds).squeeze()  # (nper, ncol)
    np.save(run_dir / "heads.npy", hds)
    final = hds[-1]
    print(f"[collect] heads shape={hds.shape}  final water table={np.array2string(final, precision=3)}")
    _plot(hds, run_dir / "coupled_gwl.png")
    return hds


def _plot(hds: np.ndarray, out: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    final = hds[-1]
    ncol = final.size
    x = np.arange(ncol)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(x, final, marker="o", color="#1f77b4", label="final head")
    ax.scatter([0, ncol - 1], [final[0], final[-1]], color="crimson", zorder=5,
               label="CHD channels")
    ax.axhline(0.0, color="grey", lw=0.7, ls=":")
    ax.set_xlabel("MODFLOW column")
    ax.set_ylabel("groundwater head [m] (datum = surface)")
    ax.set_title("Coupled SWAP<->MODFLOW6: final water table\n"
                 "(two channels + SWAP recharge on interior cells)")
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"[plot] wrote {out}")


def main(libswap: Path, libmf6: Path, case_dir: Path, run_dir: Path,
         ncol: int = 10, nper: int = 1096) -> None:
    cwd0 = Path.cwd()
    try:
        stage_run_dir(libswap, libmf6, case_dir, run_dir, ncol, nper)
        run_coupled(run_dir)
        collect_and_plot(run_dir)
    finally:
        os.chdir(cwd0)
    print("[run_coupled] coupled run complete")


if __name__ == "__main__":
    main(
        Path(sys.argv[1]).resolve(),
        Path(sys.argv[2]).resolve(),
        Path(sys.argv[3]).resolve(),
        Path(sys.argv[4]).resolve(),
        ncol=int(sys.argv[5]) if len(sys.argv) > 5 else 10,
        nper=int(sys.argv[6]) if len(sys.argv) > 6 else 1096,
    )
