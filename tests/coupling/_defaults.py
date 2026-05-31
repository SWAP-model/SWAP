"""Shared default-path discovery for the coupling smoke scripts.

Lets ``test_swap_xmi_smoke.py`` / ``test_coupled_smoke.py`` run with zero CLI
args (the convenient ``pixi run -e test test-xmi`` / ``test-coupling`` path)
while still honouring explicit positional args (the meson-registered + manual
invocations). Paths are derived from the script location so the smokes work
from any CWD.
"""
from __future__ import annotations

import ctypes.util
import os
import sys
from pathlib import Path

# tests/coupling/_defaults.py -> parents[2] == repo root
REPO_ROOT = Path(__file__).resolve().parents[2]
CASE_DIR = REPO_ROOT / "tests" / "coupling" / "hupselbrook_coupled"
LIBSWAP_XMI = REPO_ROOT / "builddir" / "libswap_xmi.so"


def default_libswap() -> Path:
    """Path to the built SWAP XMI kernel; clear error if not built yet."""
    if not LIBSWAP_XMI.exists():
        sys.exit(
            f"error: libswap_xmi.so not found at {LIBSWAP_XMI}\n"
            "       build it first:  pixi run build-linux"
        )
    return LIBSWAP_XMI


def default_libmf6() -> Path:
    """Discover libmf6.so: sys.prefix/lib, then $CONDA_PREFIX/lib, then the
    system loader search path."""
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
    searched = "\n".join(f"         {c}" for c in candidates)
    sys.exit(
        "error: libmf6.so not found. Searched:\n"
        f"{searched}\n"
        "         ctypes.util.find_library('mf6') -> None\n"
        "       it ships in the pixi test env; run via 'pixi run -e test ...'"
    )


def default_case_dir() -> Path:
    return CASE_DIR
