"""CAPI in-memory smoke: load SWAP TOML config in Python and pass it to SWAP
via cffi (no filesystem read on the Fortran side), run in headless mode to
completion, and inspect state through the capi zero-copy accessor surface and
the water-balance struct.

Pure-cffi proof of concept — pyswap (a separate arc) provides the class-based
API for real use.
"""

import pathlib
import struct

import cffi

HERE = pathlib.Path(__file__).resolve().parent


def _read_doubles(ffi, ptr, n):
    """Read n doubles from a cffi double* into a Python list."""
    raw = bytes(ffi.buffer(ptr, n * 8))
    return list(struct.unpack(f"{n}d", raw))


def test_capi_in_memory_run(swap_lib, bindings_case, monkeypatch):
    ffi = cffi.FFI()
    ffi.cdef((HERE / "swap_bmi.h").read_text())
    ffi.cdef((HERE / "swap_capi.h").read_text())
    lib = ffi.dlopen(str(swap_lib))

    # Load TOML config in Python (no filesystem read on the Fortran side).
    toml_buf = (bindings_case / "swap.toml").read_text().encode("utf-8")

    assert lib.swap_set_headless(1) == 0, "swap_set_headless failed"

    # Stay in the case dir so relative paths in the TOML resolve.
    monkeypatch.chdir(bindings_case)

    rc = lib.swap_initialize_from_toml_string(toml_buf, len(toml_buf))
    assert rc == 0, f"swap_initialize_from_toml_string returned {rc}"

    t_start = ffi.new("double *")
    t_end = ffi.new("double *")
    lib.get_start_time(t_start)
    lib.get_end_time(t_end)
    assert t_end[0] > t_start[0]

    # Initial soil moisture profile via capi zero-copy accessor.
    ptr = ffi.new("double **")
    n = ffi.new("int *")
    rc = lib.swap_view_array(b"theta\0", ptr, n)
    assert rc == 0, f"swap_view_array(theta) returned {rc}"
    assert n[0] > 0

    # Run loop until the simulation ends.
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
    assert n_steps < max_steps, "simulation did not reach flRunEnd"

    theta_final = _read_doubles(ffi, ptr[0], n[0])
    assert 0.0 < theta_final[0] < 1.0, f"theta[0] out of (0, 1): {theta_final[0]}"

    # Water-balance summary struct is readable.
    summary = ffi.new("swap_water_balance_t *")
    lib.swap_get_water_balance(summary)

    assert lib.finalize() == 0
