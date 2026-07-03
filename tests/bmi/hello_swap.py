"""SS-DRV Phase 1: BMI hello-world.

Loads libswap.so via cffi, drives the hupselbrook simulation
through one initialize / update / finalize cycle, and asserts the
two sentinel get_value_double variables return physically plausible
data. End-to-end proof that the C binding path works.
"""

import pathlib
import sys

import cffi

HERE = pathlib.Path(__file__).resolve().parent


def main() -> int:
    if len(sys.argv) < 2:
        print("usage: hello_swap.py <path-to-libswap.so>", file=sys.stderr)
        return 2

    so_path = sys.argv[1]
    header_text = (HERE / "swap_bmi.h").read_text()

    ffi = cffi.FFI()
    ffi.cdef(header_text)
    lib = ffi.dlopen(so_path)

    rc = lib.initialize(b"swap.toml\0", 0)
    assert rc == 0, f"initialize returned {rc}"

    rc = lib.update()
    assert rc == 0, f"update returned {rc}"

    t = ffi.new("double *")
    assert lib.get_current_time(t) == 0
    print(f"current_time after 1 step = {t[0]}")

    dt = ffi.new("double *")
    assert lib.get_time_step(dt) == 0
    print(f"time_step                = {dt[0]}")

    buf = ffi.new("double[500]")
    assert lib.get_value_double(b"soil_water_content\0", 500, buf) == 0
    theta0 = buf[0]
    print(f"theta[0]                 = {theta0}")
    assert 0.0 < theta0 < 1.0, f"theta[0]={theta0} not in (0, 1)"

    # Unknown variable name must return rc != 0
    rc = lib.get_value_double(b"nonexistent_variable\0", 1, buf)
    assert rc != 0, "expected nonzero rc for unknown variable name"

    rc = lib.finalize()
    assert rc == 0, f"finalize returned {rc}"

    print("BMI hello-world: OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
