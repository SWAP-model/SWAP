#!/usr/bin/env bash
# Wrapper for the unit-swap-tests pFUnit binary.
#
# Why this exists: the binary aggregates ~50 test suites into one process,
# but several suites (hupselbrook_parity, grassgrowth_parity, …) invoke
# legacy `readswap()`, which mutates `variables` module globals. Cleanup
# helpers like reset_for_next_readswap() handle most of that, but some
# state still leaks between suites and a later suite fatal-errors via
# RDDATA on what should be valid input. Each suite passes when run in
# isolation; only the cross-suite sequence is broken.
#
# Until the legacy reader is fully retired (umbrella spec SS-11), this
# wrapper invokes the binary once per suite (`-f <suite>` + `--tap
# <file>`), aggregates TAP results, and exits non-zero on any failure.
#
# Args:
#   $1 — path to the unit-swap-tests binary (passed by meson)

set -u

BINARY="${1:-./builddir/tests/unit/unit-swap-tests}"
SUITES_INC="tests/unit/testSuites.inc"

if [ ! -x "$BINARY" ]; then
   echo "run_pfunit.sh: binary not found or not executable: $BINARY" >&2
   exit 2
fi
if [ ! -f "$SUITES_INC" ]; then
   echo "run_pfunit.sh: $SUITES_INC not found (run from project root)" >&2
   exit 2
fi

TAP_DIR=$(mktemp -d)
trap 'rm -rf "$TAP_DIR"' EXIT

total_ok=0
total_fail=0
fail_lines=()

# Extract suite names. Lines look like: ADD_TEST_SUITE(test_foo_suite)
# Skip lines that start with `!` (Fortran comment — disabled suites).
suites=$(sed -n 's/^ADD_TEST_SUITE(\([^)]*\))$/\1/p' "$SUITES_INC")

for suite in $suites; do
   tap="$TAP_DIR/$suite.tap"
   "$BINARY" --tap "$tap" -f "$suite" < /dev/null > /dev/null 2>&1 || true
   if [ ! -s "$tap" ]; then
      fail_lines+=("$suite (no TAP output — binary crashed)")
      total_fail=$((total_fail + 1))
      continue
   fi
   ok=$(grep -c "^ ok " "$tap" 2>/dev/null || true)
   fail=$(grep -c "^ not ok " "$tap" 2>/dev/null || true)
   total_ok=$((total_ok + ${ok:-0}))
   total_fail=$((total_fail + ${fail:-0}))
   if [ "${fail:-0}" -gt 0 ]; then
      while IFS= read -r line; do
         fail_lines+=("${line# not ok - }")
      done < <(grep "^ not ok " "$tap")
   fi
done

echo
echo "================================================================"
echo "pFUnit summary: $total_ok passed, $total_fail failed"
echo "================================================================"

if [ ${#fail_lines[@]} -gt 0 ]; then
   echo
   echo "Failures:"
   for f in "${fail_lines[@]}"; do
      echo "  - $f"
   done
   exit 1
fi
exit 0
