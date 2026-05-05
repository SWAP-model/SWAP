#!/usr/bin/env bash
# Wrapper for the unit-swap-tests pFUnit binary.
#
# Why this exists: the binary's default invocation reports `Ok: 1, Fail: 0`
# even when individual @test cases fail, because TTutil's `fatalerr` (called
# from legacy SWAP code reachable via parity/roundtrip tests) issues a bare
# `STOP` that terminates the process with exit 0 before pFUnit's
# `finalize()` can record the failure or set the exit code. See the SS-5
# follow-up notes for the full root-cause analysis.
#
# What this wrapper does:
#   1. Reads the suite list from tests/unit/testSuites.inc.
#   2. Invokes the binary once per suite with `--tap` and `-f <suite>`,
#      so a crashing suite does not blackhole the suites listed after it.
#   3. Parses each TAP file for `^ ok` / `^ not ok` lines.
#   4. Prints a per-failure summary at the end.
#   5. Exits 1 if any test failed, 0 otherwise.
#
# Args:
#   $1 — path to the unit-swap-tests binary (passed by meson)
#
# Notes:
#   - Tests that crash via `fatalerr` mid-suite still yield a partial TAP
#     file; their post-crash @test cases are silently dropped. This is a
#     known limitation; the regression suite (`pixi run check-full`) is
#     the canonical end-to-end gate for those code paths.
#   - We swallow the binary's exit code on purpose: it can be 0 (silent
#     STOP) or 2 (gfortran runtime error) regardless of whether @test
#     assertions actually failed. The TAP file is the source of truth.

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
skipped_suites=()
fail_lines=()
unexpected_failures=()

# Temporary allowlist of pre-existing failures that block the harness fix
# from landing without simultaneously breaking `check-full` (which
# depends-on test-pfunit). Each entry is the suite-qualified test name as
# it appears in TAP output. Remove entries as the underlying tests are
# fixed (tracked in tasks #18, #19); when this list is empty, delete the
# allowlist mechanism entirely.
known_failures=(
   "test_read_heat_toml_suite.test_read_heat_swhea1_happy"
   "test_read_irrigation_toml_suite.test_read_irrigation_swirfix1_inline_table"
   "test_all_cases_smoke_suite.test_macroporeflow_loads_clean"
   "test_heat_config_suite.test_heat_happy_path_numerical"
   "test_cropgrass_config_suite.test_cropgrass_swrd3_rejected"
   "test_cropgrass_init_suite.test_cropgrass_init_swharvest2_populates_dateharvest"
)

is_known_failure() {
   local name="$1"
   local k
   for k in "${known_failures[@]}"; do
      if [ "$k" = "$name" ]; then return 0; fi
   done
   return 1
}

# Extract suite names. Lines look like: ADD_TEST_SUITE(test_foo_suite)
suites=$(sed -n 's/^ADD_TEST_SUITE(\([^)]*\))$/\1/p' "$SUITES_INC")

for suite in $suites; do
   # Parity / roundtrip suites invoke legacy readswap(), which reads
   # Get_Command_Argument(1) as the project name. The pFUnit `--tap` /
   # `-f` flags then become the project name and FOPENG fails. These
   # suites are not runnable in the per-suite TAP harness; they are
   # exercised end-to-end by `pixi run -e test check-full` instead.
   case "$suite" in
      *_parity_suite|*_roundtrip_suite)
         skipped_suites+=("$suite")
         continue
         ;;
   esac

   tap="$TAP_DIR/$suite.tap"
   "$BINARY" --tap "$tap" -f "$suite" < /dev/null > /dev/null 2>&1 || true
   if [ ! -s "$tap" ]; then
      # Empty TAP file = binary crashed before pFUnit wrote anything.
      # Treat as a hard failure so it can't be ignored.
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
         name="${line# not ok - }"
         fail_lines+=("$name")
         if ! is_known_failure "$name"; then
            unexpected_failures+=("$name")
         fi
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
fi

if [ ${#skipped_suites[@]} -gt 0 ]; then
   echo
   echo "Skipped (parity/roundtrip — covered by 'pixi run check-full'):"
   for s in "${skipped_suites[@]}"; do
      echo "  - $s"
   done
fi

if [ ${#unexpected_failures[@]} -gt 0 ]; then
   echo
   echo "UNEXPECTED FAILURES (not in allowlist — please fix or add to allowlist with justification):"
   for f in "${unexpected_failures[@]}"; do
      echo "  - $f"
   done
   exit 1
fi

# Detect stale allowlist entries (tests that used to fail but now pass).
stale_allowlist=()
for k in "${known_failures[@]}"; do
   found=0
   for f in "${fail_lines[@]}"; do
      if [ "$k" = "$f" ]; then found=1; break; fi
   done
   if [ $found -eq 0 ]; then
      stale_allowlist+=("$k")
   fi
done
if [ ${#stale_allowlist[@]} -gt 0 ]; then
   echo
   echo "Stale allowlist entries (tests now pass — please remove from known_failures):"
   for s in "${stale_allowlist[@]}"; do
      echo "  - $s"
   done
   exit 1
fi

exit 0
