# SWAP &harr; MODFLOW 6 coupling demo

A runnable example of **SWAP coupled to MODFLOW 6**: a two-channel
lateral-flow cross-section where the saturated groundwater (MODFLOW) and the
unsaturated soil columns (SWAP) exchange head and recharge once per day. SWAP
runs as an *ensemble* of N homogeneous columns (one per MODFLOW interior cell)
behind an XMI shared-library kernel; a thin fork of imod_coupler's `swapmod`
driver runs the loop.

This is a **qualitative smoke** demo: the goal is a completed coupled run with
sensible, nonzero exchange values and a recharge-driven water-table shape, not
a byte-exact numeric target.

---

## What runs

```
 MODFLOW 6 (libmf6.so)                 SWAP XMI (libswap_xmi.so)
 1 layer x 1 row x ncol, unconfined    ensemble of (ncol-2) columns,
 two CHD channels (cells 0, ncol-1)    one shared swap_config_t,
 RCHA on the (ncol-2) interior cells   swbotb=1 prescribed-GWL bottom
        ^   recharge (m/d)                    |   head (m)
        |                                     v
        +-----------  swapmod driver  --------+
        (MODFLOW leads the clock; SWAP solves one day per step)
```

MODFLOW leads the clock. Each day the driver:
1. `mf6.prepare_time_step` (populates the RCHA package's `NODELIST`),
2. injects MODFLOW interior-cell heads into SWAP (`gwl`, m),
3. SWAP solves the day for every column,
4. reads SWAP recharge (`qbot_volume`, m/step) and storage (`storage_coef`),
   writes them onto the MODFLOW RCHA `recharge` (as a rate) and `STO`,
5. iterates MODFLOW to convergence.

---

## Prerequisites

Everything is in the **pixi `test` environment**:

- `xmipy` (drives the XMI C-ABI),
- `flopy` + a MODFLOW 6 binary / `libmf6.so` (conda-forge `modflow6`),
- the built SWAP XMI kernel `builddir/libswap_xmi.so`
  (`pixi run build-linux`).

`libmf6.so` ships in the env at `.pixi/envs/test/lib/libmf6.so`.

> **Load `libswap_xmi.so`, NOT `libswap_bmi.so`.** Only the XMI kernel exports
> `get_value_ptr` / `solve` and the ensemble exchange variables.

---

## How to run

Build the kernel, then run the orchestrator from the repo root:

```bash
pixi run build-linux

# short wiring run (ncol=10 -> 8 interior columns, 30 days):
pixi run -e test python tests/coupling/run_coupled.py \
    builddir/libswap_xmi.so \
    .pixi/envs/test/lib/libmf6.so \
    tests/swap-cases/toml/1.hupselbrook \
    /tmp/swapmf_run 10 30

# full run (2002-2004, 1096 days):
pixi run -e test python tests/coupling/run_coupled.py \
    builddir/libswap_xmi.so .pixi/envs/test/lib/libmf6.so \
    tests/swap-cases/toml/1.hupselbrook /tmp/swapmf_full 10 1096
```

`run_coupled.py` builds the MODFLOW model in `RUN_DIR/mf`, stages the SWAP work
dir in `RUN_DIR/swap`, writes `RUN_DIR/imod_coupler.toml` with the discovered
absolute DLL paths, runs the loop, saves `heads.npy`, and plots the final water
table to `coupled_gwl.png`.

### Smoke test

```bash
pixi run -e test python tests/coupling/test_coupled_smoke.py \
    builddir/libswap_xmi.so .pixi/envs/test/lib/libmf6.so \
    tests/swap-cases/toml/1.hupselbrook
```

It is also registered as the meson test `coupled-smoke` (suite `coupling`),
which discovers `libmf6.so` at `.pixi/envs/test/lib/libmf6.so` at configure
time; if that path is absent the test is not registered (run the script
manually). The SWAP-only XMI contract test is `swap-xmi-smoke` (no MODFLOW).

```bash
pixi run meson test -C builddir --suite coupling -v
```

> **Test-env binding caveat.** Meson resolves `python3` at *configure* time via
> `find_installation('python3')`, which binds the **default** pixi env's
> interpreter. That interpreter lacks `xmipy`/`flopy` (they live only in the
> `test` env), so `meson test` runs both coupling tests under a python without
> the deps and they error. This affects `swap-xmi-smoke` and `coupled-smoke`
> equally. **Run the smoke scripts directly via `pixi run -e test python ...`
> (above)** — that is the supported path and what the qualitative DoD is
> verified against.

---

## Exchange variables and unit / sign conventions

The SWAP XMI kernel exposes three rank-1 exchange arrays (length = column
count = `ncol-2`), addressed by their **SWAP-native** names:

| variable        | dir.      | meaning                                              |
|-----------------|-----------|------------------------------------------------------|
| `gwl`           | MF6 -> SWAP | groundwater head, **metres**, datum = surface = MODFLOW `top` = 0 (below surface is negative). Injected as SWAP's prescribed-GWL bottom boundary. |
| `qbot_volume`   | SWAP -> MF6 | recharge **depth over the day, metres**, **positive = into groundwater** (SWAP `qbot < 0` downward percolation, sign-flipped). |
| `storage_coef`  | SWAP -> MF6 | specific yield (-), fixed `0.15` placeholder for the smoke. |

**Recharge unit/sign (the #1 correctness point).** The driver converts the
per-day depth to a rate: `mf6_recharge[:] = qbot_volume[:] / delt` (m/d). The
MODFLOW RCHA `recharge` variable is a **flux rate (L/T)** that MF6 multiplies
by cell area internally, so there is **no area factor** in the exchange. Sign:
positive `qbot_volume` (downward percolation) becomes positive RCHA recharge,
which **mounds** the water table. Verified in the full run: the interior
mean head rose from -0.75 m to -0.45 m, peaking near -0.31 m -- a recharge dome
between the two channels.

---

## SWAP coupled config deltas

`tests/swap-cases/toml/1.hupselbrook/swap_coupled.toml` is the standalone
hupselbrook config with two deltas for coupling:

- `[bottom_boundary] swbotb = 1` -- prescribed groundwater level; in coupled
  mode `flcoupled_gwl` makes BoundBottom use the injected `gwl` instead of the
  `gwl_file` table (the table still seeds the initial level / passes
  validation).
- `[drainage] swdra = 0` -- lateral drainage off (MODFLOW now carries lateral
  groundwater flow between the channels).

## `ensemble.txt` rule

The SWAP XMI kernel reads its **column count** from `ensemble.txt` (a single
integer) staged next to the config in the SWAP work dir. It **must equal the
MODFLOW interior recharge-cell count = `ncol - 2`**, because the driver's
`[:] = [:]` exchanges require SWAP array length == RCHA `NBOUND`. `run_coupled.py`
writes `ensemble.txt` automatically.

---

## Vendored driver fork

`imod_coupler_fork/` is a **vendored, self-contained fork** of imod_coupler's
`swap`-branch `swapmod` driver:

- `swap_wrapper.py` -- `SwapWrapper(XmiWrapper)` with SWAP-native exchange
  names (`gwl` / `qbot_volume` / `storage_coef`), and an `initialize()` that
  defaults to `swap_coupled.toml` (the driver calls `initialize()` with no
  argument).
- `swapmod.py` -- `SwapMod(Driver)` plus a minimal `Mf6Wrapper`. Fork-specific
  adaptations validated against this `libmf6` (7.x) build:
  - **`Mf6Wrapper.get_value_ptr` uppercases** the address -- MODFLOW 6 stores
    every memory-manager variable in UPPERCASE, while the config supplies the
    model name as written (`swapmf`).
  - the MODFLOW head (`SWAPMF/X`, size = `ncol`) and storage (`SWAPMF/STO/SS`,
    size = `ncol`) arrays are **full-grid**; the recharge package and the SWAP
    ensemble cover only the interior cells. The driver maps between them via
    the RCHA package's `NODELIST` (0-based), so head/storage exchanges address
    exactly the recharge cells (`RECHARGE` is already interior-sized).
  - `prepare_time_step` is called **before** the head exchange each step, since
    MF6 only populates `NODELIST` during prepare.
- `config.py` -- pydantic-v2 `SwapModConfig` (`kernels.{modflow6,swap}` +
  `coupling[]`); `driver.py` -- the `Driver` ABC.

The validated MODFLOW 6 BMI addresses (model `swapmf`):

| purpose       | address                  |
|---------------|--------------------------|
| head          | `SWAPMF/X`               |
| recharge      | `SWAPMF/RCHA/RECHARGE`   |
| recharge cells| `SWAPMF/RCHA/NODELIST`   |
| storage (SS)  | `SWAPMF/STO/SS`          |
| has sc1       | `SWAPMF/STO/ISTOR_COEF`  |
| cell area     | `SWAPMF/DIS/AREA`        |
| cell top/bot  | `SWAPMF/DIS/TOP` / `/BOT`|
| max outer iter| `SLN_1/MXITER`           |

---

## Files

- `build_modflow.py` -- flopy two-channel CHD + RCHA MODFLOW 6 builder.
- `run_coupled.py` -- orchestrator (build + stage + run + plot).
- `test_swap_xmi_smoke.py` -- SWAP-only XMI contract test (`swap-xmi-smoke`).
- `test_coupled_smoke.py` -- coupled qualitative smoke (`coupled-smoke`).
- `imod_coupler.toml` -- coupler config (DLL paths rewritten at runtime).
- `imod_coupler_fork/` -- vendored driver + wrappers + config models.
- `ensemble.txt` -- example column-count sidecar.
