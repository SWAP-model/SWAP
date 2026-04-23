# src/

Source code for the SWAP modernization.

Architecture: [../docs/architecture.md](../docs/architecture.md).
State pattern: [../docs/state-management.md](../docs/state-management.md).
Code style: [../docs/code-style.md](../docs/code-style.md).
Contributing: [../docs/contributing.md](../docs/contributing.md).

## Subdirectories

| Directory | What it owns |
|---|---|
| [atmosphere/](atmosphere/README.md) | Precipitation, interception, ET, snow, meteo. |
| [boundary/](boundary/README.md) | Top / bottom boundary conditions. |
| [core/](core/README.md) | Entry points, aggregated state, time control, legacy bridge. |
| [crop/](crop/README.md) | Crop growth (fixed / grass / WOFOST), irrigation, tillage, rootextraction. |
| [drainage/](drainage/README.md) | Drainage flux, surface water state. |
| [error/](error/README.md) | (Empty — Phase 4 placeholder.) |
| [heat/](heat/README.md) | Soil temperature, frozen soil conductivity. |
| [io/](io/README.md) | TOML and legacy readers, CSV output writers. |
| [macropore/](macropore/README.md) | Macropore flow, rates, output. |
| [soil/](soil/README.md) | Soil hydraulics, grid, water balance. |
| [solute/](solute/README.md) | Solute transport. |
| [utils/](utils/README.md) | Arrays, solvers, shared helpers. |

## Dependency direction

Physics subdirs depend on `core/` (state aggregator, time control) and `utils/` (low-level helpers). They do not depend on each other directly; inter-domain communication is through `swap_state_t`.

## Conventions

See [../docs/code-style.md](../docs/code-style.md). The linter config is `.fprettify.rc` at the repo root; run `pixi run lint-check` against files you change before committing.
