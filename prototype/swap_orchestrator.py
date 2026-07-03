#!/usr/bin/env python3
"""SWAP column orchestrator — drive many SWAP columns from Python via netCDF.

This is the first slice of the gridded/ensemble orchestration layer described
in dev-docs/2026-06-11-python-orchestration-plan.md. It runs N independent SWAP
columns from a single Python process by:

  1. Reading a *grid spec* (an xarray.Dataset / netCDF file, or a plain list of
     ColumnSpec) that describes, per column, a base TOML config plus optional
     per-column parameter overrides.
  2. Running each column **in-memory** through the SWAP BMI C-ABI
     (libswap.so) — no per-column scratch files, no subprocess shell-out.
  3. Running the columns **in parallel**, one OS process per column, via a
     multiprocessing pool. Process isolation sidesteps the remaining Tier-3
     module-global state (crop_config_global, a few SAVE locals) that still
     blocks true in-process threading — see ADR 0043. Each worker loads a fresh
     copy of the library, so columns never share Fortran state.
  4. Collecting every column's in-memory results record (zero-copy numpy views
     of the Fortran results matrix) into one tidy xarray.Dataset indexed by
     (column, time), and writing it to netCDF.

The *same* machinery serves two use cases the maintainers called out:

  * **Gridded simulation** — `column` is a spatial cell; overrides carry per-cell
    soil / boundary / meteo parameters. Couple to MODFLOW by feeding the
    `groundwater_level` override and reading `bottom_flux` back each step (that
    tighter loop is the XMI ensemble path, swap_ensemble_mod.f90; this module is
    the loose / offline driver).
  * **PEST / sensitivity / calibration** — `column` is a parameter *sample*;
    overrides carry the sampled values. Run a Morris/Sobol design or a PEST
    run-manager batch as one orchestrated sweep and read the objective columns
    straight out of the returned Dataset.

Design intent: the orchestration logic here is pure Python over a thin, stable
C-ABI. It does not import pyswap. pyswap's role (see the plan doc) is to *build*
the base TOML + companions from its Pydantic components; this module *runs* them.

Smoke test:  python prototype/swap_orchestrator.py --selftest
"""
from __future__ import annotations

import argparse
import ctypes
import io
import os
import sys
from dataclasses import dataclass, field
from multiprocessing import get_context
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

REPO = Path(__file__).resolve().parent.parent
DEFAULT_LIB = REPO / "builddir" / "libswap.so"

# Companion-file extensions SWAP loads at init/seed time alongside swap.toml.
# Anything matching these in the case dir is attached to the in-memory store.
_COMPANION_GLOBS = ("*.crp.toml", "*.dra.toml", "*.csv", "*.bbc.toml", "*.ini")


# --------------------------------------------------------------------------- #
# Column specification
# --------------------------------------------------------------------------- #
@dataclass
class ColumnSpec:
    """One SWAP column to run.

    name        : label used as the `column` coordinate value.
    base_dir    : directory holding swap.toml + companion files.
    overrides   : nested {table: {key: value}} edits applied to the base TOML
                  before the run (e.g. {"simulation.numerical": {"dtmax": 0.02}}
                  or {"soil.initial": {"gwli": -120.0}}). Dotted table paths are
                  supported as keys ("a.b.c"); values are written verbatim.
    companions  : explicit companion filenames; auto-detected from base_dir when
                  None.
    """

    name: str
    base_dir: Path
    overrides: Mapping[str, Mapping[str, Any]] = field(default_factory=dict)
    companions: Sequence[str] | None = None

    def resolved_companions(self) -> list[str]:
        if self.companions is not None:
            return list(self.companions)
        found: list[str] = []
        for pat in _COMPANION_GLOBS:
            found += [p.name for p in sorted(Path(self.base_dir).glob(pat))]
        return found


@dataclass
class ColumnResult:
    """Outcome of one column run."""

    name: str
    ok: bool
    error: str = ""
    times: np.ndarray = field(default_factory=lambda: np.empty(0))
    columns: list[str] = field(default_factory=list)
    data: np.ndarray = field(default_factory=lambda: np.empty((0, 0)))
    water_balance: dict[str, float] = field(default_factory=dict)


# --------------------------------------------------------------------------- #
# TOML override
# --------------------------------------------------------------------------- #
def apply_overrides(base_toml: str, overrides: Mapping[str, Mapping[str, Any]]) -> str:
    """Return base_toml with the given table/key edits applied.

    Uses tomlkit so the rest of the document (comments, ordering, arrays-of-
    tables) is preserved byte-for-byte — only the targeted scalars change.
    """
    if not overrides:
        return base_toml
    import tomlkit

    doc = tomlkit.parse(base_toml)
    for table_path, kv in overrides.items():
        node: Any = doc
        for part in table_path.split("."):
            if part not in node:
                node[part] = tomlkit.table()
            node = node[part]
        for key, value in kv.items():
            node[key] = value
    return tomlkit.dumps(doc)


# --------------------------------------------------------------------------- #
# Water-balance struct mirrored from swap_capi_mod.f90 swap_water_balance_t
# --------------------------------------------------------------------------- #
class _WaterBalance(ctypes.Structure):
    _fields_ = [
        (n, ctypes.c_double)
        for n in (
            "rain", "evap_pot", "evap_act", "transp_pot", "transp_act",
            "runoff", "drain", "percolation", "storage_change", "balance_error",
        )
    ]


def _load_lib(lib_path: Path) -> ctypes.CDLL:
    lib = ctypes.CDLL(str(lib_path))
    lib.swap_set_headless.argtypes = [ctypes.c_int]
    lib.swap_clear_config_files.argtypes = []
    lib.swap_attach_config_file.argtypes = [
        ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_int
    ]
    lib.swap_initialize_from_toml_string.argtypes = [ctypes.c_char_p, ctypes.c_int]
    lib.update.argtypes = []
    lib.finalize.argtypes = []
    lib.swap_get_scalar.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_double)]
    lib.swap_get_water_balance.argtypes = [ctypes.POINTER(_WaterBalance)]
    lib.swap_results_shape.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
    lib.swap_view_results.argtypes = [
        ctypes.POINTER(ctypes.c_void_p), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)
    ]
    lib.swap_view_results_times.argtypes = [ctypes.POINTER(ctypes.c_void_p), ctypes.POINTER(ctypes.c_int)]
    lib.swap_results_columns.argtypes = [ctypes.c_char_p, ctypes.c_int]
    return lib


def _run_one(spec: ColumnSpec, lib_path: Path) -> ColumnResult:
    """Run a single column to completion in-memory. Executed inside a worker
    process; loads its own copy of the library."""
    try:
        lib = _load_lib(lib_path)
        lib.swap_set_headless(1)
        lib.swap_clear_config_files()

        base_dir = Path(spec.base_dir)
        for name in spec.resolved_companions():
            content = (base_dir / name).read_bytes()
            lib.swap_attach_config_file(name.encode(), len(name), content, len(content))

        toml_text = (base_dir / "swap.toml").read_text()
        toml_text = apply_overrides(toml_text, spec.overrides)
        # Opt in to the in-memory results record (zero-copy readout below).
        if "results_in_memory" not in toml_text:
            toml_text += "\n[output.csv]\nresults_in_memory = 1\n"
        buf = toml_text.encode()

        rc = lib.swap_initialize_from_toml_string(buf, len(buf))
        if rc != 0:
            return ColumnResult(spec.name, ok=False, error=f"init rc={rc}")

        flrunend = ctypes.c_double(0.0)
        for _ in range(5_000_000):
            lib.swap_get_scalar(b"flRunEnd", ctypes.byref(flrunend))
            if flrunend.value != 0.0:
                break
            lib.update()

        wb = _WaterBalance()
        lib.swap_get_water_balance(ctypes.byref(wb))
        water_balance = {n: getattr(wb, n) for n, _ in wb._fields_}

        lib.finalize()  # flush the in-memory results record

        nr, nc = ctypes.c_int(0), ctypes.c_int(0)
        lib.swap_results_shape(ctypes.byref(nr), ctypes.byref(nc))
        nrows, ncols = nr.value, nc.value

        cols_buf = ctypes.create_string_buffer(1 << 16)
        lib.swap_results_columns(cols_buf, len(cols_buf))
        col_names = [s.decode() for s in cols_buf.raw.split(b"\x00") if s][:ncols]

        data = np.empty((nrows, ncols))
        times = np.empty(nrows)
        if nrows > 0 and ncols > 0:
            ptr = ctypes.c_void_p()
            lib.swap_view_results(ctypes.byref(ptr), ctypes.byref(nr), ctypes.byref(nc))
            # Fortran column-major (nrows x ncols): element (i,j) at j*nrows+i.
            flat = np.ctypeslib.as_array(
                (ctypes.c_double * (nrows * ncols)).from_address(ptr.value)
            )
            data = flat.reshape((ncols, nrows)).T.copy()

            tptr = ctypes.c_void_p()
            ntimes = ctypes.c_int(0)
            lib.swap_view_results_times(ctypes.byref(tptr), ctypes.byref(ntimes))
            if tptr.value:
                times = np.ctypeslib.as_array(
                    (ctypes.c_double * ntimes.value).from_address(tptr.value)
                ).copy()

        return ColumnResult(
            spec.name, ok=True, times=times, columns=col_names,
            data=data, water_balance=water_balance,
        )
    except Exception as exc:  # noqa: BLE001 — report, don't crash the pool
        return ColumnResult(spec.name, ok=False, error=f"{type(exc).__name__}: {exc}")


# Module-level worker shim (must be picklable for spawn).
def _worker(args):
    spec, lib_path = args
    return _run_one(spec, lib_path)


def run_columns(
    specs: Sequence[ColumnSpec],
    lib_path: Path = DEFAULT_LIB,
    processes: int | None = None,
) -> list[ColumnResult]:
    """Run every column, in parallel (one process each). Order is preserved."""
    lib_path = Path(lib_path)
    if not lib_path.exists():
        raise FileNotFoundError(f"SWAP BMI library not found: {lib_path}")
    payload = [(s, lib_path) for s in specs]
    if processes == 1 or len(specs) == 1:
        return [_worker(p) for p in payload]
    ctx = get_context("spawn")  # fresh interpreter per worker → no shared state
    with ctx.Pool(processes=processes) as pool:
        return pool.map(_worker, payload)


# --------------------------------------------------------------------------- #
# netCDF / xarray adapters
# --------------------------------------------------------------------------- #
def results_to_dataset(results: Sequence[ColumnResult]):
    """Pack column results into one xarray.Dataset indexed by (column, time).

    Output vars are the union of SWAP result columns; the time axis is taken
    from the first OK column (all columns share the simulation calendar).
    Water-balance scalars are added as (column,) data vars prefixed `wb_`.
    """
    import xarray as xr

    ok = [r for r in results if r.ok and r.data.size]
    if not ok:
        raise RuntimeError(
            "no column produced results; errors: "
            + "; ".join(f"{r.name}: {r.error}" for r in results if not r.ok)
        )
    ref = ok[0]
    ntime = ref.data.shape[0]
    var_names = ref.columns
    col_coord = [r.name for r in results]

    data_vars: dict[str, Any] = {}
    for j, var in enumerate(var_names):
        arr = np.full((len(results), ntime), np.nan)
        for i, r in enumerate(results):
            if r.ok and r.data.shape == (ntime, len(var_names)):
                arr[i, :] = r.data[:, j]
        data_vars[var] = (("column", "time"), arr)

    wb_keys = next((r.water_balance.keys() for r in ok), [])
    for key in wb_keys:
        data_vars[f"wb_{key}"] = (
            ("column",),
            np.array([r.water_balance.get(key, np.nan) if r.ok else np.nan for r in results]),
        )
    data_vars["ok"] = (("column",), np.array([r.ok for r in results]))

    # SWAP result times are days since 1900-01-01 (the t1900 axis).
    time_coord = ref.times[:ntime] if ref.times.size >= ntime else np.arange(ntime)
    coords = {
        "column": col_coord,
        "time": ("time", time_coord, {"units": "days since 1900-01-01", "axis": "T"}),
    }
    return xr.Dataset(data_vars, coords=coords)


def specs_from_dataset(grid, base_dir: Path, override_vars: Sequence[str] | None = None):
    """Build ColumnSpec list from an xarray grid spec.

    The grid must have a `column` dimension. Each data variable listed in
    `override_vars` (default: all 1-D vars over `column`) becomes a per-column
    override, keyed by its attribute `swap_path` (e.g. "simulation.numerical")
    and the var name as the TOML key. Lets a calibration design or a gridded
    parameter field map straight onto SWAP inputs.
    """
    columns = [str(c) for c in grid["column"].values]
    if override_vars is None:
        override_vars = [v for v in grid.data_vars if grid[v].dims == ("column",)]
    specs = []
    for i, col in enumerate(columns):
        ov: dict[str, dict[str, Any]] = {}
        for v in override_vars:
            table = grid[v].attrs.get("swap_path")
            if not table:
                continue
            key = grid[v].attrs.get("swap_key", v)
            ov.setdefault(table, {})[key] = grid[v].values[i].item()
        specs.append(ColumnSpec(name=col, base_dir=base_dir, overrides=ov))
    return specs


def write_netcdf(path: Path, results: Sequence[ColumnResult]) -> None:
    results_to_dataset(results).to_netcdf(path)


# --------------------------------------------------------------------------- #
# Self-test
# --------------------------------------------------------------------------- #
def _selftest() -> int:
    """Run a tiny gridded sweep against a local case and assert it works."""
    case = REPO / "tests" / "regression" / "cases" / "swcf3" / "toml"
    if not (case / "swap.toml").exists():
        print(f"SELFTEST SKIP: case not found at {case}")
        return 0
    if not DEFAULT_LIB.exists():
        print(f"SELFTEST SKIP: build libswap.so first ({DEFAULT_LIB})")
        return 0

    # 3 columns: baseline + two with a perturbed initial groundwater level.
    # This is exactly the shape of a 1-parameter sensitivity sweep.
    specs = [
        ColumnSpec("gwli_-75", case),
        ColumnSpec("gwli_-100", case, overrides={"soil.initial": {"gwli": -100.0}}),
        ColumnSpec("gwli_-120", case, overrides={"soil.initial": {"gwli": -120.0}}),
    ]
    print(f"Running {len(specs)} columns in parallel via {DEFAULT_LIB.name} ...")
    results = run_columns(specs)

    for r in results:
        status = "ok" if r.ok else f"FAIL ({r.error})"
        wb = r.water_balance
        extra = ""
        if r.ok:
            extra = (f"  rows={r.data.shape[0]} cols={len(r.columns)} "
                     f"drain={wb.get('drain', float('nan')):.3f} "
                     f"balerr={wb.get('balance_error', float('nan')):.2e}")
        print(f"  {r.name:12s} {status}{extra}")

    assert all(r.ok for r in results), "a column failed"
    assert all(r.data.shape[0] > 0 for r in results), "a column produced no rows"
    assert all(abs(r.water_balance["balance_error"]) < 1e-3 for r in results), \
        "water balance did not close"
    # The override must take effect: a different initial groundwater level
    # changes the early-transient GWL (the cumulative balance washes the
    # initial condition out over a multi-year run, so assert on the transient).
    gwl_idx = results[0].columns.index("GWL")
    first_gwl = [r.data[0, gwl_idx] for r in results]
    assert any(abs(g - first_gwl[0]) > 1.0 for g in first_gwl[1:]), \
        "gwli override had no effect — overrides not applied?"
    print(f"  override check: first-row GWL per column = "
          f"{[round(g, 2) for g in first_gwl]}")

    ds = results_to_dataset(results)
    out = REPO / "prototype" / "_selftest_orchestrator.nc"
    ds.to_netcdf(out)
    reloaded = _reload_check(out)
    print(f"netCDF written: {out.name}  dims={dict(ds.sizes)}  vars={len(ds.data_vars)}")
    print(f"netCDF reloaded OK: columns={list(reloaded['column'].values)}")
    out.unlink(missing_ok=True)
    print("SELFTEST PASS")
    return 0


def _reload_check(path: Path):
    import xarray as xr

    with xr.open_dataset(path) as ds:
        ds.load()
        return ds


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selftest", action="store_true", help="run the bundled smoke test")
    ap.add_argument("--lib", type=Path, default=DEFAULT_LIB, help="path to libswap.so")
    args = ap.parse_args(argv)
    if args.selftest:
        return _selftest()
    ap.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
