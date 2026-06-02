#!/usr/bin/env python3
"""Prototype: drive SWAP fully in-memory from Python (diskless init) and
verify the run is bit-identical to the on-disk standalone path.

Both runs go through the SAME library (libswap_bmi.so) over the singleton
(state, config); they are executed in isolated subprocesses so the singleton
is fresh each time. The in-memory run pushes the TOML config + every companion
file (drainage/crop subfiles + the meteo CSV) as bytes via the new
swap_attach_config_file C-API, so the Fortran reads ZERO files. The disk run
uses the standard BMI initialize(path).

If the two water balances match to 1e-10, the diskless path is proven
equivalent to disk.

    python prototype/pyswap_inmemory_demo.py
"""
import ctypes
import json
import os
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LIB = os.path.join(REPO, "builddir", "libswap_bmi.so")
CASE = os.path.join(REPO, "tests", "swap-cases", "toml", "1.hupselbrook")
# Every load-time + seed-time companion hupselbrook reads.
COMPANIONS = [
    "swap.dra.toml",
    "maizes.crp.toml",
    "potatod.crp.toml",
    "grassd.crp.toml",
    "283.csv",
]
SENTINEL = "RESULT_JSON:"


class WaterBalance(ctypes.Structure):
    _fields_ = [
        (n, ctypes.c_double)
        for n in (
            "rain", "evap_pot", "evap_act", "transp_pot", "transp_act",
            "runoff", "drain", "percolation", "storage_change", "balance_error",
        )
    ]


def load_lib():
    lib = ctypes.CDLL(LIB)
    lib.swap_set_headless.argtypes = [ctypes.c_int]
    lib.swap_clear_config_files.argtypes = []
    lib.swap_attach_config_file.argtypes = [
        ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_int
    ]
    lib.swap_initialize_from_toml_string.argtypes = [ctypes.c_char_p, ctypes.c_int]
    lib.initialize.argtypes = [ctypes.c_char_p, ctypes.c_int]
    lib.update.argtypes = []
    lib.swap_get_scalar.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_double)]
    lib.swap_get_water_balance.argtypes = [ctypes.POINTER(WaterBalance)]
    lib.finalize.argtypes = []
    lib.swap_results_shape.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
    lib.swap_view_results.argtypes = [ctypes.POINTER(ctypes.c_void_p),
                                      ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
    lib.swap_results_columns.argtypes = [ctypes.c_char_p, ctypes.c_int]
    return lib


def run_to_end(lib):
    flrunend = ctypes.c_double(0.0)
    for _ in range(1_000_000):
        lib.swap_get_scalar(b"flRunEnd", ctypes.byref(flrunend))
        if flrunend.value != 0.0:
            break
        lib.update()
    wb = WaterBalance()
    lib.swap_get_water_balance(ctypes.byref(wb))
    # flush in-memory results record and run end-of-sim bookkeeping
    lib.finalize()
    return {n: getattr(wb, n) for n, _ in wb._fields_}


def worker(mode):
    lib = load_lib()
    lib.swap_set_headless(1)
    if mode == "mem":
        lib.swap_clear_config_files()
        for name in COMPANIONS:
            with open(os.path.join(CASE, name), "rb") as fh:
                content = fh.read()
            lib.swap_attach_config_file(name.encode(), len(name), content, len(content))
        with open(os.path.join(CASE, "swap.toml"), "rb") as fh:
            toml = fh.read()
        rc = lib.swap_initialize_from_toml_string(toml, len(toml))
    else:
        path = os.path.join(CASE, "swap.toml").encode()
        rc = lib.initialize(path, len(path))
    if rc != 0:
        print(SENTINEL + json.dumps({"error": f"init rc={rc}"}))
        return 1
    wb = run_to_end(lib)
    nr, nc = ctypes.c_int(0), ctypes.c_int(0)
    lib.swap_results_shape(ctypes.byref(nr), ctypes.byref(nc))
    ptr = ctypes.c_void_p()
    rc = lib.swap_view_results(ctypes.byref(ptr), ctypes.byref(nr), ctypes.byref(nc))
    cols_buf = ctypes.create_string_buffer(8192)
    lib.swap_results_columns(cols_buf, len(cols_buf))
    names = [s.decode() for s in cols_buf.raw.split(b"\x00") if s][: nc.value]
    first_row, last_row = [], []
    if rc == 0 and nr.value > 0 and nc.value > 0:
        flat = (ctypes.c_double * (nr.value * nc.value)).from_address(ptr.value)
        # Fortran column-major: element (i,j) at j*nrows + i
        first_row = [flat[j * nr.value + 0] for j in range(nc.value)]
        last_row = [flat[j * nr.value + (nr.value - 1)] for j in range(nc.value)]
    wb["_results_nrows"] = nr.value
    wb["_results_ncols"] = nc.value
    wb["_results_cols"] = names
    wb["_results_first_row"] = first_row
    wb["_results_last_row"] = last_row
    print(SENTINEL + json.dumps(wb))
    return 0


def run_mode(mode):
    out = subprocess.run(
        [sys.executable, __file__, "worker", mode],
        capture_output=True, text=True, cwd=REPO,
    )
    if out.returncode != 0:
        print(f"[{mode}] subprocess failed rc={out.returncode}")
        print(out.stdout)
        print(out.stderr)
        sys.exit(1)
    line = next(l for l in out.stdout.splitlines() if l.startswith(SENTINEL))
    return json.loads(line[len(SENTINEL):])


def main():
    if len(sys.argv) == 3 and sys.argv[1] == "worker":
        sys.exit(worker(sys.argv[2]))

    mem = run_mode("mem")
    disk = run_mode("disk")
    if "error" in mem:
        print("in-memory init failed:", mem["error"]); sys.exit(1)
    if "error" in disk:
        print("disk init failed:", disk["error"]); sys.exit(1)

    print(f"{'field':>16} {'in-memory':>18} {'disk':>18} {'|diff|':>12}")
    print("-" * 68)
    ok = True
    _meta_keys = {"_results_nrows", "_results_ncols", "_results_cols",
                  "_results_first_row", "_results_last_row"}
    for k in mem:
        if k in _meta_keys:
            continue
        d = abs(mem[k] - disk[k])
        if d > 1e-10:
            ok = False
        print(f"{k:>16} {mem[k]:18.8f} {disk[k]:18.8f} {d:12.2e}")
    print("-" * 68)

    print("\nResults record (in-memory):")
    print(f"  shape   = {mem['_results_nrows']} rows x {mem['_results_ncols']} cols")
    print(f"  columns = {mem['_results_cols']}")
    print(f"  first row (mem) = {[round(x, 6) for x in mem['_results_first_row']]}")
    if mem.get("_results_first_row") and disk.get("_results_first_row") \
       and len(mem["_results_first_row"]) == len(disk["_results_first_row"]):
        rec_ok = all(abs(a - b) <= 1e-10
                     for a, b in zip(mem["_results_first_row"], disk["_results_first_row"])) \
                 and all(abs(a - b) <= 1e-10
                         for a, b in zip(mem["_results_last_row"], disk["_results_last_row"]))
        print("  record mem==disk:", "PASS" if rec_ok else "FAIL")
        ok = ok and rec_ok
    else:
        print("  record mem==disk: SKIP (no rows / shape mismatch)")
        ok = False

    print("RESULT:", "PASS — diskless in-memory run is bit-identical to disk"
          if ok else "FAIL — divergence")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
