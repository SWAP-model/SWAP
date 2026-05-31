"""Self-contained fork of the imod_coupler ``swap`` branch swapmod driver.

This package vendors a minimal subset of Deltares' imod_coupler (the ``swap``
branch) so that the SWAP <-> MODFLOW 6 coupling driver can be imported and run
*without* installing the full imod_coupler package. Two adaptations are made
relative to upstream:

* the kernel wrapper (:class:`~.swap_wrapper.SwapWrapper`) exposes SWAP-native
  ``get_value_ptr`` names (``gwl`` / ``qbot_volume`` / ``storage_coef``); and
* the ``Driver`` ABC, ``BaseConfig`` and ``ExchangeCollector`` are vendored as
  minimal local stand-ins, and ``loguru`` is replaced by the stdlib ``logging``
  module, so the only third-party imports are ``xmipy``, ``numpy`` and
  ``pydantic``.
"""

from __future__ import annotations

try:  # imported as a package
    from .driver import Driver
    from .swapmod import SwapMod
except ImportError:  # imported flat (dir on sys.path)
    from driver import Driver
    from swapmod import SwapMod

__all__ = ["Driver", "SwapMod"]
