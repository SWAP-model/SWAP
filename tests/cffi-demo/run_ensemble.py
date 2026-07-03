"""SS-BMI2: in-memory workflow demo.

Loads SWAP TOML config in Python; passes it to SWAP via cffi
(no filesystem touches for config); runs SWAP in headless mode;
inspects state via the capi accessor surface; finalizes.

Usage:
    python run_ensemble.py <case_dir> <path_to_libswap.so>

Example:
    python tests/cffi-demo/run_ensemble.py \\
           tests/swap-cases/toml/1.hupselbrook \\
           builddir/libswap.so

Pure-cffi proof of concept — pyswap rewrite (separate arc) replaces
this with a proper class-based API.
"""

import pathlib
import sys

import struct

import cffi

HERE = pathlib.Path(__file__).resolve().parent


def _read_doubles(ffi, ptr, n):
    """Read n doubles from a cffi double* into a Python list."""
    raw = bytes(ffi.buffer(ptr, n * 8))
    return list(struct.unpack(f"{n}d", raw))


def main(case_dir: pathlib.Path, lib_path: pathlib.Path) -> int:
    ffi = cffi.FFI()
    ffi.cdef((HERE / "swap_bmi.h").read_text())
    ffi.cdef((HERE / "swap_capi.h").read_text())
    lib = ffi.dlopen(str(lib_path))

    # Load TOML config in Python (no filesystem read on the Fortran side)
    toml_path = case_dir / "swap.toml"
    toml_text = toml_path.read_text()
    toml_buf = toml_text.encode("utf-8")
    print(f"Loaded TOML from {toml_path} ({len(toml_buf)} bytes)")

    # Wire up headless + initialize
    assert lib.swap_set_headless(1) == 0, "swap_set_headless failed"

    # Stay in the case directory so any relative paths in the TOML
    # (e.g. drainage .csv references) resolve correctly.
    import os
    os.chdir(case_dir)

    rc = lib.swap_initialize_from_toml_string(toml_buf, len(toml_buf))
    assert rc == 0, f"swap_initialize_from_toml_string returned {rc}"

    # Read time bounds via BMI
    t_start = ffi.new("double *")
    t_end   = ffi.new("double *")
    lib.get_start_time(t_start)
    lib.get_end_time(t_end)
    print(f"Run window: t1900 = {t_start[0]:.4f} -> {t_end[0]:.4f}")

    # Read initial soil moisture profile via capi zero-copy accessor
    ptr = ffi.new("double **")
    n = ffi.new("int *")
    rc = lib.swap_view_array(b"theta\0", ptr, n)
    assert rc == 0, f"swap_view_array(theta) returned {rc}"
    theta = _read_doubles(ffi, ptr[0], n[0])
    print(f"Initial theta profile: {n[0]} nodes, theta[0]={theta[0]:.4f}, "
          f"theta[-1]={theta[-1]:.4f}")

    # Run loop until simulation ends
    flag = ffi.new("double *")
    n_steps = 0
    max_steps = 200000
    while n_steps < max_steps:
        rc = lib.update()
        assert rc == 0, f"update returned {rc}"
        n_steps += 1
        lib.swap_get_scalar(b"flRunEnd\0", flag)
        if flag[0] > 0.5:
            break

    print(f"Completed {n_steps} update steps")

    # Read final t1900 + theta
    t_cur = ffi.new("double *")
    lib.get_current_time(t_cur)
    print(f"Final t1900 = {t_cur[0]:.4f}")
    theta_final = _read_doubles(ffi, ptr[0], n[0])
    print(f"Final theta[0] = {theta_final[0]:.4f}")
    assert 0.0 < theta_final[0] < 1.0, f"theta[0] out of (0, 1): {theta_final[0]}"

    # Read water balance summary struct
    summary = ffi.new("swap_water_balance_t *")
    lib.swap_get_water_balance(summary)
    print(f"Water balance summary:")
    print(f"  rain        = {summary.rain:.4f}")
    print(f"  transp_act  = {summary.transp_act:.4f}")
    print(f"  runoff      = {summary.runoff:.4f}")
    print(f"  percolation = {summary.percolation:.4f}")

    assert lib.finalize() == 0
    print("BMI Phase 2 demo: OK")
    return 0


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("usage: run_ensemble.py <case_dir> <path_to_libswap.so>",
              file=sys.stderr)
        sys.exit(2)
    case = pathlib.Path(sys.argv[1]).resolve()
    lib  = pathlib.Path(sys.argv[2]).resolve()
    sys.exit(main(case, lib))
