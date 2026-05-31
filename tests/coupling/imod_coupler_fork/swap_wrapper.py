"""Fork of imod_coupler swap-branch SwapWrapper using SWAP-native var names."""
import os
from pathlib import Path
from typing import Any, Union

import numpy as np
from numpy.typing import NDArray
from xmipy import XmiWrapper

#: Default coupled SWAP config filename staged in the SWAP working dir. The
#: SWAP XMI kernel reads its ensemble column count from ``ensemble.txt`` next
#: to this config; both live in the kernel's working_directory.
DEFAULT_SWAP_CONFIG = "swap_coupled.toml"


class SwapWrapper(XmiWrapper):
    def __init__(self, lib_path, lib_dependency=None, working_directory=None, timing=False):
        super().__init__(lib_path, lib_dependency, working_directory, timing)

    def initialize(self, config_file: Union[str, os.PathLike[Any]] = "") -> None:
        """Initialize SWAP, defaulting to the staged coupled config.

        The imod_coupler driver calls ``initialize()`` with no argument; the
        SWAP XMI kernel needs the config path (it reads ``ensemble.txt`` beside
        it for the column count). xmipy chdirs to the working_directory first,
        so a bare filename resolves there.
        """
        if not config_file:
            config_file = DEFAULT_SWAP_CONFIG
        super().initialize(config_file)

    def get_head_ptr(self) -> NDArray[np.float64]:
        return self.get_value_ptr("gwl")            # was 'dhgwmod'

    def get_volume_ptr(self) -> NDArray[np.float64]:
        return self.get_value_ptr("qbot_volume")    # was 'dvsim'

    def get_storage_ptr(self) -> NDArray[np.float64]:
        return self.get_value_ptr("storage_coef")   # was 'dsc1sim'
