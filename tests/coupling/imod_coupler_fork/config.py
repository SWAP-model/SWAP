"""Config models for the vendored swapmod driver.

Vendored from ``imod_coupler/drivers/swapmod/config.py`` (``swap`` branch) and
adapted from pydantic v1 (``@validator``) to pydantic v2 (``@field_validator``).
The minimal local :class:`BaseConfig` replaces ``imod_coupler.config.BaseConfig``
so the fork imports without the full imod_coupler package.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any, List, Optional

from pydantic import BaseModel, DirectoryPath, FilePath, field_validator


class BaseConfig(BaseModel):
    """Minimal stand-in for ``imod_coupler.config.BaseConfig``."""

    log_level: str = "INFO"
    timing: bool = False


class Kernel(BaseModel):
    dll: FilePath
    dll_dep_dir: Optional[DirectoryPath] = None
    work_dir: DirectoryPath

    @field_validator("dll")
    @classmethod
    def resolve_dll(cls, dll: FilePath) -> FilePath:
        return dll.resolve()

    @field_validator("dll_dep_dir")
    @classmethod
    def resolve_dll_dep_dir(
        cls, dll_dep_dir: Optional[DirectoryPath]
    ) -> Optional[DirectoryPath]:
        if dll_dep_dir is not None:
            dll_dep_dir = dll_dep_dir.resolve()
        return dll_dep_dir


class Kernels(BaseModel):
    modflow6: Kernel
    swap: Kernel


class Coupling(BaseModel):
    mf6_model: str  # the MODFLOW 6 model that will be coupled
    mf6_swap_recharge_pkg: str  # the recharge package that will be used for coupling
    output_config_file: Optional[FilePath] = None

    @field_validator("output_config_file")
    @classmethod
    def resolve_file_path(cls, file_path: Optional[FilePath]) -> Optional[FilePath]:
        if file_path is not None:
            file_path = file_path.resolve()
        return file_path


class SwapModConfig(BaseModel):
    kernels: Kernels
    coupling: List[Coupling]

    def __init__(self, config_dir: Path, **data: Any) -> None:
        """Model for the SwapMod config validated by pydantic.

        The validation expects the current working directory at config-file
        level, so it is changed during initialization.

        Args:
            config_dir (Path): Directory where the config file resides.
        """
        os.chdir(config_dir)
        super().__init__(**data)

    @field_validator("coupling")
    @classmethod
    def restrict_coupling_count(cls, coupling: List[Coupling]) -> List[Coupling]:
        if len(coupling) == 0:
            raise ValueError("At least one coupling has to be defined.")
        if len(coupling) > 1:
            raise ValueError("Multi-model coupling is not yet supported.")
        return coupling
