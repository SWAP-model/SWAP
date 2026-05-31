"""Fork of imod_coupler swap-branch SwapWrapper using SWAP-native var names."""
from pathlib import Path
from typing import Union

import numpy as np
from numpy.typing import NDArray
from xmipy import XmiWrapper


class SwapWrapper(XmiWrapper):
    def __init__(self, lib_path, lib_dependency=None, working_directory=None, timing=False):
        super().__init__(lib_path, lib_dependency, working_directory, timing)

    def get_head_ptr(self) -> NDArray[np.float64]:
        return self.get_value_ptr("gwl")            # was 'dhgwmod'

    def get_volume_ptr(self) -> NDArray[np.float64]:
        return self.get_value_ptr("qbot_volume")    # was 'dvsim'

    def get_storage_ptr(self) -> NDArray[np.float64]:
        return self.get_value_ptr("storage_coef")   # was 'dsc1sim'
