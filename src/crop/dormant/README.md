# src/crop/dormant — dormant crop-module features

Crop-module features whose **legacy bodies are preserved** but which
have **no live dispatch site** in the TOML-pipeline build. The files
in this directory are deliberately **excluded from the meson build**:
they exist as source-controlled reference for future reactivation.

Why not delete? Each feature is a documented capability of the
legacy SWAP code (alternative drought/oxygen formulations, microscopic
uptake, etc.) that the TOML migration did not re-wire. Keeping the
bodies here — instead of leaving them as unreferenced public exports
in active modules — has two benefits:

1. The main modules' `use variables, only:` blocks stop being held
   hostage by dead routines. Each dormant module owns its own
   imports; retiring a bare global from `variables.f90` no longer
   needs to consider dormant-routine references.
2. The reactivation checklist for each feature is co-located with
   the body, not buried in an ADR.

## Current dormant modules

| File | Feature | Reactivation prerequisite |
|------|---------|----------------------------|
| `jongvanlier.f90` | De Jong van Lier (2013) microscopic root-water uptake — `swdrought = 2` path; contains `JongvanLier` (outer Newton-Raphson on `hleaf`/`Tactual`) and `JongvanLierLoop` (per-node soil-root pressure-head balance + matric-flux extraction). Reads `criterhr`, `stephr`, `kroot`, `kstem`, `rxylem`, `rootradius`, `rootcoefa`, `rooteff`, `flhydrlift`, `twilt`, `wiltpoint`. | Restore the two `fatalerr_collected` dispatch stubs in `rootextraction.f90` (search for `'swdrought=2'`): (1) replace the stub at the top of `RootExtraction` with `call JongvanLier(state)`, (2) restore `alpdry = soil%alpJvLier` in the per-node combination loop. Wire `criterhr`/`stephr`/`kroot`/`kstem`/`rxylem`/`rootradius`/`rootcoefa`/`rooteff` through a new `crop_jvl_t` sub-record on `state%cfg%crop` (`cropfixed_config_t` has slots but no TOML reader). Drop the dormant module's `use variables` blocks. Add this file to both `meson.build` and `tests/unit/meson.build` under the Crop section. Restore the retired legacy declarations (`rootcoefa`, `rooteff`, `rootradius`, `kstem`) and orphan stubs (`CriterHr`, `Kroot`, `Rxylem`, `StepHr`) — or move them directly to `state%cfg%crop`. |
| `oxygenrepro.f90` | Bartholomeus reproduction-function oxygen stress — `swoxygen=2` + `swoxygentype=2` path; contains `oxygen_dat` (loads the 4 hard-coded 18×6 coefficient tables for top/sub × slope/intercept) and `OxygenReproFunction` (per-node `rwu_factor` as a polynomial in soil temperature, depth, and mean gas-filled porosity above the node). Reads `state%mesh%zbotcp`. | Restore the dispatch site in `RootExtraction` (the `fatalerr_collected('RootExtraction', 'swoxygen=2/swoxygentype=2 (OxygenReproFunction) is dormant — …')` stub near rootextraction.f90:155) with the original `oxygen_dat` + `OxygenReproFunction` calls. Add the two 6-element arrays (`OxygenSlope`, `OxygenIntercept`) to `state%crop%oxygen` (or a new `crop_oxygen_repro_state_t`) — the legacy `variables.f90` declarations for them were already retired with the Bartholomeus migration. Wire `swtopsub` and `nrstaring` from the per-rotation config (cropfixed/cropwofost already have schema slots under their `swoxygen=2` sub-record, but the validator stub-errors that branch today). Add this file to `meson.build` and `tests/unit/meson.build` under the Crop section, and add `use oxygenrepro_dormant_mod, only: oxygen_dat, OxygenReproFunction` to `rootextraction.f90`. Add a TOML regression fixture before relying on the path — the existing 6 regression cases do not cover it. |

## Build exclusion

These files are **not listed** in the top-level `meson.build` sources
list. They will not be compiled. To reactivate a feature:

1. Work through the reactivation prerequisite for that file (above).
2. Add the file path to `meson.build`.
3. Add a `use <module>_mod, only: <routines>` line wherever you wire
   the dispatch site.
4. Run `pixi run check-fast` — broken imports will surface immediately
   because the dormant body's `use variables, only:` lines reference
   globals that may have been retired in the meantime.

## Style note

The dormant bodies are preserved verbatim from the active state at the
moment of extraction (Pattern 1 sub-record associate applied, stale
migration breadcrumbs stripped, bare-global imports retained). They
are **not** kept in lockstep with subsequent refactors of the active
modules — the reactivation step is the right time to harmonize them
again.
