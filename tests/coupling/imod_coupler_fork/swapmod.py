"""SwapMod: the coupling between SWAP and MODFLOW 6.

Vendored from ``imod_coupler/drivers/swapmod/swapmod.py`` (``swap`` branch) and
adapted to be self-contained:

* imports the local :class:`~.swap_wrapper.SwapWrapper` fork (SWAP-native
  exchange variable names);
* imports the local :class:`~.driver.Driver` ABC and local ``BaseConfig`` /
  ``Coupling`` / ``SwapModConfig`` instead of the full imod_coupler package;
* uses the stdlib ``logging`` module instead of ``loguru``;
* replaces ``imod_coupler.logging.exchange_collector.ExchangeCollector`` with a
  trivial local no-op stub; and
* vendors a minimal ``Mf6Wrapper`` (the MODFLOW 6 side of the coupling) locally,
  since only the SWAP wrapper is forked in this task.

The driver loop: MODFLOW 6 leads the clock, SWAP solves once per step, and the
exchanges move groundwater head (m, MF6 -> SWAP) and recharge volume + storage
coefficient (SWAP -> MF6).
"""

from __future__ import annotations

import logging
from typing import Any

from numpy.typing import NDArray
from xmipy import XmiWrapper

try:  # imported as a package: ``import imod_coupler_fork``
    from .config import BaseConfig, Coupling, SwapModConfig
    from .driver import Driver
    from .swap_wrapper import SwapWrapper
except ImportError:  # imported flat: dir on sys.path, ``import swapmod``
    from config import BaseConfig, Coupling, SwapModConfig
    from driver import Driver
    from swap_wrapper import SwapWrapper

logger = logging.getLogger(__name__)


class ExchangeCollector:
    """No-op stand-in for ``imod_coupler.logging.exchange_collector``.

    The real collector records exchanged arrays to file; this fork drops that
    behaviour to avoid pulling in the full imod_coupler logging stack.
    """

    @classmethod
    def from_file(cls, output_config_file: Any) -> "ExchangeCollector":
        return cls()

    def finalize(self) -> None:
        pass


class Mf6Wrapper(XmiWrapper):
    """Minimal MODFLOW 6 XMI wrapper for the SWAP coupling.

    Vendored helper methods mirror the upstream
    ``imod_coupler.kernelwrappers.mf6_wrapper.Mf6Wrapper`` surface that
    :meth:`SwapMod.couple` relies on.

    MODFLOW 6 stores every memory-manager variable address in UPPERCASE; the
    coupling config supplies the model/package names in their case-as-written
    (e.g. ``swapmf``), so we uppercase the full address before handing it to
    the BMI ``get_value_ptr`` (validated against libmf6 7.x).
    """

    def get_value_ptr(self, name: str) -> NDArray[Any]:  # type: ignore[override]
        return super().get_value_ptr(name.upper())

    def get_head(self, mf6_flowmodel_key: str) -> NDArray[Any]:
        return self.get_value_ptr(f"{mf6_flowmodel_key}/X")

    def get_recharge(
        self, mf6_flowmodel_key: str, mf6_pkg_key: str
    ) -> NDArray[Any]:
        return self.get_value_ptr(f"{mf6_flowmodel_key}/{mf6_pkg_key}/RECHARGE")[:]

    def get_recharge_nodes(
        self, mf6_flowmodel_key: str, mf6_pkg_key: str
    ) -> NDArray[Any]:
        """Zero-copy pointer to the recharge package's 1-based node list.

        MODFLOW stores the boundary cell list as 1-based node numbers in
        ``<MODEL>/<PKG>/NODELIST`` and only fills it during prepare_time_step,
        so callers must convert to 0-based (``- 1``) at use time, after the
        first prepare.
        """
        return self.get_value_ptr(f"{mf6_flowmodel_key}/{mf6_pkg_key}/NODELIST")

    def get_storage(self, mf6_flowmodel_key: str) -> NDArray[Any]:
        return self.get_value_ptr(f"{mf6_flowmodel_key}/STO/SS")

    def has_sc1(self, mf6_flowmodel_key: str) -> bool:
        return bool(self.get_value_ptr(f"{mf6_flowmodel_key}/STO/ISTOR_COEF")[0])

    def get_area(self, mf6_flowmodel_key: str) -> NDArray[Any]:
        return self.get_value_ptr(f"{mf6_flowmodel_key}/DIS/AREA")

    def get_top(self, mf6_flowmodel_key: str) -> NDArray[Any]:
        return self.get_value_ptr(f"{mf6_flowmodel_key}/DIS/TOP")

    def get_bot(self, mf6_flowmodel_key: str) -> NDArray[Any]:
        return self.get_value_ptr(f"{mf6_flowmodel_key}/DIS/BOT")

    def max_iter(self) -> Any:
        return self.get_value_ptr("SLN_1/MXITER")[0]


class SwapMod(Driver):
    """The driver coupling SWAP and MODFLOW 6."""

    base_config: BaseConfig  # the parsed information from the configuration file
    swapmod_config: SwapModConfig  # the parsed config specific to SwapMod
    coupling: Coupling  # the coupling information

    timing: bool  # true, when timing is enabled
    mf6: Mf6Wrapper  # the MODFLOW 6 XMI kernel
    swap: SwapWrapper  # the SWAP XMI kernel

    max_iter: NDArray[Any]  # max. nr outer iterations in MODFLOW kernel
    delt: float  # time step from MODFLOW 6 (leading)

    mf6_head: NDArray[Any]  # the hydraulic head array in the coupled model
    mf6_recharge: NDArray[Any]  # the coupled recharge array from the RCH package
    mf6_storage: NDArray[Any]  # the specific storage array (ss)
    mf6_has_sc1: bool  # when true, mf6 storage is a storage coefficient (sc1)
    mf6_area: NDArray[Any]  # cell area (size:nodes)
    mf6_top: NDArray[Any]  # top of cell (size:nodes)
    mf6_bot: NDArray[Any]  # bottom of cell (size:nodes)

    swap_head: NDArray[Any]  # internal SWAP groundwater head
    swap_volume: NDArray[Any]  # unsaturated zone flux (as a volume!)
    swap_storage: NDArray[Any]  # SWAP storage coefficients (MODFLOW's sc1)

    def __init__(self, base_config: BaseConfig, swapmod_config: SwapModConfig):
        """Constructs the `SwapMod` object."""
        self.base_config = base_config
        self.swapmod_config = swapmod_config
        self.coupling = swapmod_config.coupling[
            0
        ]  # Adapt as soon as we have multimodel support

    def initialize(self) -> None:
        self.mf6 = Mf6Wrapper(
            lib_path=self.swapmod_config.kernels.modflow6.dll,
            lib_dependency=self.swapmod_config.kernels.modflow6.dll_dep_dir,
            working_directory=self.swapmod_config.kernels.modflow6.work_dir,
            timing=self.base_config.timing,
        )
        self.swap = SwapWrapper(
            lib_path=self.swapmod_config.kernels.swap.dll,
            lib_dependency=self.swapmod_config.kernels.swap.dll_dep_dir,
            working_directory=self.swapmod_config.kernels.swap.work_dir,
            timing=self.base_config.timing,
        )
        # Print output to stdout
        self.mf6.set_int("ISTDOUTTOFILE", 0)
        self.mf6.initialize()
        self.swap.initialize()
        self.log_version()
        if self.coupling.output_config_file is not None:
            self.exchange_logger = ExchangeCollector.from_file(
                self.coupling.output_config_file
            )
        else:
            self.exchange_logger = ExchangeCollector()
        self.couple()

    def log_version(self) -> None:
        logger.info(f"MODFLOW version: {self.mf6.get_version()}")
        logger.info(f"SWAP version: {self.swap.get_version()}")

    def couple(self) -> None:
        """Couple Modflow and SWAP.

        The MODFLOW head/storage arrays are full-grid (size = nodes); the
        recharge package and the SWAP column ensemble cover only the interior
        recharge cells (size = nbound = ncol-2 for the two-channel demo). We
        derive the interior-cell index map from the recharge package's CELLID
        / NODELIST so the SWAP<->MODFLOW head and storage exchanges address
        exactly the recharge cells (whereas RECHARGE is already nbound-sized).
        """
        self.mf6_head = self.mf6.get_head(self.coupling.mf6_model)
        self.mf6_recharge = self.mf6.get_recharge(
            self.coupling.mf6_model, self.coupling.mf6_swap_recharge_pkg
        )
        self.mf6_storage = self.mf6.get_storage(self.coupling.mf6_model)
        self.mf6_has_sc1 = self.mf6.has_sc1(self.coupling.mf6_model)
        self.mf6_area = self.mf6.get_area(self.coupling.mf6_model)
        self.mf6_top = self.mf6.get_top(self.coupling.mf6_model)
        self.mf6_bot = self.mf6.get_bot(self.coupling.mf6_model)
        self.max_iter = self.mf6.max_iter()
        # Zero-copy pointer to the recharge package's 1-based node list. MF6
        # only populates it during prepare_time_step, so we keep the pointer
        # and dereference (0-based) lazily inside the exchange routines.
        self.rch_nodelist = self.mf6.get_recharge_nodes(
            self.coupling.mf6_model, self.coupling.mf6_swap_recharge_pkg
        )

        self.swap_head = self.swap.get_head_ptr()
        self.swap_volume = self.swap.get_volume_ptr()
        self.swap_storage = self.swap.get_storage_ptr()

    def update(self) -> None:
        # Prepare the MODFLOW time step FIRST: this populates the recharge
        # package's NODELIST (the interior-cell index map) that the head and
        # storage exchanges below depend on. (we cannot set the timestep yet
        # in Modflow -> pass the dummy value 0.0.)
        self.mf6.prepare_time_step(0.0)

        # heads to SWAP
        self.exchange_mod2swap()

        self.delt = self.mf6.get_time_step()
        self.swap.prepare_time_step(self.delt)
        self.swap.prepare_solve(0)
        self.swap.solve(0)
        self.swap.finalize_solve(0)
        self.exchange_swap2mod()

        # convergence loop
        self.mf6.prepare_solve(1)
        for kiter in range(1, self.max_iter + 1):
            has_converged = self.do_iter(1)
            if has_converged:
                logger.debug(f"MF6-SWAP converged in {kiter} iterations")
                break
        self.mf6.finalize_solve(1)

        self.mf6.finalize_time_step()
        self.swap.finalize_time_step()

    def finalize(self) -> None:
        self.mf6.finalize()
        self.swap.finalize()
        self.exchange_logger.finalize()

    def get_current_time(self) -> float:
        return self.mf6.get_current_time()

    def get_end_time(self) -> float:
        return self.mf6.get_end_time()

    def exchange_swap2mod(self) -> None:
        """Exchange SWAP to Modflow.

        Storage maps SWAP's per-column specific yield onto the recharge cells'
        MODFLOW storage entries; recharge is SWAP's per-day percolation depth
        (m) converted to a rate (m/d) for the RCHA package (already
        nbound-sized). No area factor: MODFLOW RCHA `recharge` is a flux rate
        (L/T) that MF6 multiplies by cell area internally.
        """
        self.mf6_storage[self._rch_idx()] = self.swap_storage[:]

        # Divide recharge volume (m over the step) by delta time -> rate (m/d).
        tled = 1 / self.delt
        self.mf6_recharge[:] = tled * self.swap_volume[:]

    def exchange_mod2swap(self) -> None:
        """Exchange Modflow to SWAP (head of the recharge cells, m)."""
        self.swap_head[:] = self.mf6_head[self._rch_idx()]

    def _rch_idx(self) -> NDArray[Any]:
        """0-based grid-node indices of the recharge cells (resolved lazily).

        ``rch_nodelist`` is MF6's 1-based NODELIST, only populated during
        prepare_time_step; convert to 0-based at use time.
        """
        return self.rch_nodelist[:] - 1

    def do_iter(self, sol_id: int) -> bool:
        """Execute a single iteration."""
        has_converged = self.mf6.solve(sol_id)
        return has_converged

    def report_timing_totals(self) -> None:
        total_mf6 = self.mf6.report_timing_totals()
        total_swap = self.swap.report_timing_totals()
        total = total_mf6 + total_swap
        logger.info(f"Total elapsed time in numerical kernels: {total:0.4f} seconds")
