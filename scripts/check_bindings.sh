#!/usr/bin/env bash
# check_bindings.sh — the C-ABI binding gate for binding/driver arcs.
#
# Drives the one shared library (builddir/libswap.so, ADR 0050) through every
# consumer contract in one command, so a binding change is verified end-to-end:
#   1. BMI single-column        (tests/bmi/test_bmi_smoke.py)
#   2. CAPI in-memory           (tests/cffi-demo/test_capi_smoke.py)
#   3. XMI ensemble smoke        (tests/coupling/test_swap_xmi_smoke.py)
#   4. SWAP<->MODFLOW6 coupled   (tests/coupling/test_coupled_smoke.py; auto-
#                                 skipped by the libmf6 fixture if absent)
#
# These are pytest tests (shared fixtures in tests/conftest.py resolve the
# library and stage a self-contained case). This script just points pytest at
# them with the built library and runs them together. The exe byte-identical
# regression + pFUnit stay in `pixi run -e test check-fast` (run that too).
#
# Env overrides: SWAP_LIB (library path), SWAP_BINDINGS_CASE (case staged for
# the BMI/CAPI smokes; defaults to the runnable_model unit fixture).
set -euo pipefail

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
lib="${SWAP_LIB:-$repo/builddir/libswap.so}"

[ -f "$lib" ] || { echo "ERROR: $lib not built (run: pixi run build-linux)"; exit 1; }

export SWAP_LIB="$lib"
exec pytest "$repo/tests/bmi" "$repo/tests/cffi-demo" "$repo/tests/coupling" -v
