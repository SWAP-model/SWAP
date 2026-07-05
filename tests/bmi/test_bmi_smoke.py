"""BMI hello-world smoke: load libswap via cffi and drive one
initialize / update / finalize cycle over the runnable_model case, asserting
the sentinel get_value_double variables return physically plausible data.
End-to-end proof that the C binding path works.
"""

import pathlib

import cffi

HERE = pathlib.Path(__file__).resolve().parent


def test_bmi_hello_world(swap_lib, bindings_case, monkeypatch):
    ffi = cffi.FFI()
    ffi.cdef((HERE / "swap_bmi.h").read_text())
    lib = ffi.dlopen(str(swap_lib))

    # initialize() reads swap.toml from CWD; run inside the staged case dir.
    monkeypatch.chdir(bindings_case)

    rc = lib.initialize(b"swap.toml\0", 0)
    assert rc == 0, f"initialize returned {rc}"

    rc = lib.update()
    assert rc == 0, f"update returned {rc}"

    t = ffi.new("double *")
    assert lib.get_current_time(t) == 0

    dt = ffi.new("double *")
    assert lib.get_time_step(dt) == 0

    buf = ffi.new("double[500]")
    assert lib.get_value_double(b"soil_water_content\0", 500, buf) == 0
    theta0 = buf[0]
    assert 0.0 < theta0 < 1.0, f"theta[0]={theta0} not in (0, 1)"

    # Unknown variable name must return rc != 0
    rc = lib.get_value_double(b"nonexistent_variable\0", 1, buf)
    assert rc != 0, "expected nonzero rc for unknown variable name"

    rc = lib.finalize()
    assert rc == 0, f"finalize returned {rc}"
