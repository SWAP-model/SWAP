# src/

Source code for the SWAP modernization.

Architecture: [../docs/architecture.md](../docs/architecture.md).
State pattern: [../docs/state-management.md](../docs/state-management.md).
Code style: [../docs/code-style.md](../docs/code-style.md).
Contributing: [../docs/contributing.md](../docs/contributing.md).

## Layout (ADR 0048)

Feature-first for behavior, layer-first for the shared data model.

**Foundation & framing**

| Directory | What it owns |
|---|---|
| [core/](core/README.md) | Foundation leaves only: constants, array dims, date utils, logging, array helpers, numerical solvers. |
| [driver/](driver/README.md) | Lifecycle orchestration (`swap_mod`) and the `swap` executable (`swap_main`, `swap_ensemble_mod`). |
| [bindings/](bindings/README.md) | C-ABI facades: `swap_capi`, `swap_bmi`, `swap_xmi`, `bmi_constants`. |
| [error/](error/README.md) | Error collection + `fatalerr`. |
| [validation/](validation/README.md) | Config validation helpers. |

**Shared data model & edges**

| Directory | What it owns |
|---|---|
| [config/](config/README.md) | Read-only typed config records (`*_config_t`), the parsed TOML result. |
| [state/](state/README.md) | Mutable typed state records (`*_state_t`), aggregated under `swap_state_t`. |
| [io/](io/README.md) | TOML readers and CSV readers/writers — all file I/O. |

**Compute feature folders**

| Directory | What it owns |
|---|---|
| [atmosphere/](atmosphere/README.md) | Precipitation, interception, ET, snow, meteo. |
| [soilwater/](soilwater/README.md) | Soil hydraulics, grid, water balance, top/bottom boundary conditions. |
| [heat/](heat/README.md) | Soil temperature, frozen-soil conductivity. |
| [solute/](solute/README.md) | Solute transport. |
| [drainage/](drainage/README.md) | Drainage flux, surface-water dynamics. |
| [timecontrol/](timecontrol/README.md) | Time loop, calendar, schedule, dt control. |
| [crop/](crop/README.md) | Crop growth dispatch + cross-mode (rootextraction, oxygenstress) + irrigation + tillage, with `fixed/`, `grass/`, `wofost/` sub-packages. All WOFOST nutrient dynamics live in `crop/wofost/`. |

## Dependency direction

Compute feature folders depend on the foundation (`core/`) and the shared data
model (`state/`, `config/`); they do not depend on each other directly —
inter-domain communication is through `swap_state_t`. `driver/` sits at the apex
(depends on everything); `bindings/` wraps `driver/` with the C-ABI.

## Conventions

See [../docs/code-style.md](../docs/code-style.md). The linter config is `.fprettify.rc` at the repo root; run `pixi run lint-check` against files you change before committing.
