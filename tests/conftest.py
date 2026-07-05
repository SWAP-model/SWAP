"""Shared pytest fixtures for the Python test harnesses.

Covers the C-ABI binding smokes (tests/bmi, tests/cffi-demo) and the coupling
smokes (tests/coupling). The byte-identical regression suite (tests/regression)
resolves its own paths via regression_harness and does not use these.
"""
from __future__ import annotations

import ctypes.util
import os
import shutil
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
RUNNABLE_MODEL = REPO_ROOT / "tests" / "unit" / "fixtures" / "runnable_model"
COUPLED_CASE = REPO_ROOT / "tests" / "coupling" / "hupselbrook_coupled"


@pytest.fixture(scope="session")
def swap_lib() -> Path:
    """Path to the built libswap.so (SWAP_LIB overrides). Skip if not built."""
    env = os.environ.get("SWAP_LIB")
    lib = Path(env) if env else REPO_ROOT / "builddir" / "libswap.so"
    if not lib.exists():
        pytest.skip(f"libswap.so not built at {lib} (run: pixi run build-linux)")
    return lib


@pytest.fixture(scope="session")
def libmf6() -> Path:
    """Discover libmf6.so (sys.prefix/lib, $CONDA_PREFIX/lib, loader path).

    Skips when absent — MODFLOW 6 ships in the pixi test env, so coupled tests
    are expected to skip outside it rather than fail.
    """
    candidates = [Path(sys.prefix) / "lib" / "libmf6.so"]
    conda = os.environ.get("CONDA_PREFIX")
    if conda:
        candidates.append(Path(conda) / "lib" / "libmf6.so")
    for c in candidates:
        if c.exists():
            return c
    found = ctypes.util.find_library("mf6")
    if found:
        return Path(found)
    pytest.skip("libmf6.so not found (ships in the pixi test env: pixi run -e test ...)")


@pytest.fixture
def bindings_case(tmp_path) -> Path:
    """Stage a self-contained runnable case into a temp dir; return its path.

    Defaults to the runnable_model unit fixture (in-repo, one-year model);
    override with SWAP_BINDINGS_CASE to smoke a different case.
    """
    env = os.environ.get("SWAP_BINDINGS_CASE")
    src = Path(env) if env else RUNNABLE_MODEL
    dst = tmp_path / "case"
    shutil.copytree(src, dst)
    return dst


@pytest.fixture(scope="session")
def coupled_case_dir() -> Path:
    """The self-contained SWAP<->MODFLOW6 coupling case directory."""
    return COUPLED_CASE
