#!/usr/bin/env bash
# check_bindings.sh — the C-ABI six-path gate for binding/driver arcs (T0.5).
#
# Drives the one shared library (builddir/libswap.so, ADR 0050) through every
# consumer contract, so a binding change is verified end-to-end in one command:
#   1. BMI single-column       (tests/bmi/hello_swap.py — get_value_double)
#   2. CAPI in-memory           (tests/cffi-demo/run_ensemble.py — swap_view_array)
#   3. XMI ensemble smoke        (tests/coupling/test_swap_xmi_smoke.py)
#   4. SWAP<->MODFLOW6 coupled   (tests/coupling/test_coupled_smoke.py; skipped
#                                 if libmf6.so is absent)
# The exe byte-identical regression + pFUnit stay in `pixi run -e test check-fast`
# (run that too before committing). The meson bmi/cffi-demo suites reference the
# un-checked-out tests/swap-cases submodule; this script stages the self-contained
# regression case instead, into a temp dir (never the committed tree).
set -euo pipefail

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
lib="$repo/builddir/libswap.so"
case_src="$repo/tests/regression/cases/hupselbrook/toml"
py="python"

[ -f "$lib" ] || { echo "ERROR: $lib not built (run: pixi run build-linux)"; exit 1; }

work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
cp -r "$case_src"/. "$work/case/"

fail=0
run() { echo "=== $1 ==="; shift; if "$@"; then echo "  PASS"; else echo "  FAIL"; fail=1; fi; }

run "1. BMI hello_swap"       bash -c "cd '$work/case' && $py '$repo/tests/bmi/hello_swap.py' '$lib'"
run "2. CAPI run_ensemble"    $py "$repo/tests/cffi-demo/run_ensemble.py" "$work/case" "$lib"
run "3. XMI smoke"            $py "$repo/tests/coupling/test_swap_xmi_smoke.py" "$lib"
if [ -f "$repo/.pixi/envs/test/lib/libmf6.so" ]; then
   run "4. SWAP<->MODFLOW6 coupled smoke" $py "$repo/tests/coupling/test_coupled_smoke.py" "$lib"
else
   echo "=== 4. coupled smoke: SKIPPED (libmf6.so absent) ==="
fi

echo
if [ "$fail" -eq 0 ]; then echo "check-bindings: all paths OK"; else echo "check-bindings: FAILURES above"; fi
exit "$fail"
