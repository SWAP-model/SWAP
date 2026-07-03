# src/soil/dormant — dormant water-balance features

Soil-water features whose **legacy bodies are preserved** but which
have **no live dispatch site** in the TOML-pipeline build. The files in
this directory are deliberately **excluded from the meson build**:
they exist as source-controlled reference for future reactivation.

Why not delete? Each feature is a documented capability of the
legacy SWAP code (varying gate switches, coupling targets, perched-
table heuristics) that the TOML migration did not re-wire. Keeping
the bodies here — instead of leaving them as unreferenced public
exports in active modules — has two benefits:

1. The main modules' `use variables, only:` blocks stop being held
   hostage by dead routines. Each dormant module owns its own
   imports; retiring a bare global from `variables.f90` no longer
   needs to consider dormant-routine references.
2. The reactivation checklist for each feature is co-located with
   the body, not buried in an ADR.

> **Staleness note (2026-07-03, arc T1-B).** Reactivation steps below may name
> scaffolding that has since been deleted: `variables.f90` / bare globals (gone —
> runtime symbols live on typed `state%…` records) and the `state%cfg` pointer
> (retired, ADR 0046 — use a `config%…` argument). Read those as "wire through the
> typed `state`/`config` records"; harmonizing the frozen body is part of
> reactivation.

## Current dormant modules

| File | Feature | Reactivation prerequisite |
|------|---------|----------------------------|
| `watertable.f90` | `watertable()` — perched-water-table search with `CritUndSatVol` volume-integral heuristic (April 2008). The helper `level()` is **also** in this file (verbatim copy) so the dormant module compiles self-contained when reactivated; the **live** copy of `level()` remains in `src/soil/waterbalance.f90` because `calcgwl` calls it. On reactivation, drop the dormant `level` and import from `soilwaterbalance_mod` instead. | Add a config gate (`CritUndSatVol > 0`) and a dispatch site that calls `watertable` from `calcgwl` or the output coupling. Migrate `CritUndSatVol` from `variables.f90` to `state%cfg%soil` or `state%soilwater`. |
| `checkmassbal.f90` | Per-period mass-balance audit for ANIMO/PEARL output coupling (June 2003). Writes `<outfil>.dwb.csv` listing per-subsystem deviations that exceed `CritDevMasBal`. The macropore branches (Section 4, IQExcMtxDm1/2 accumulation, MPDOM1/MPDOM2 writes) have been stripped — ADR 0040 retired macropores permanently. | Add a config gate (e.g., `output_csv.dwb=true`), a dispatch site that calls `checkmassbal` at every output period, and migrate the 6 globals it reads (`NumNodNew`, `outfil`, `pathwork`, `DZNew`, `CritDevMasBal`, `dev_cmb`) from `variables.f90` to `state%mesh` (regrid temps) and `state%cfg%output`/`state%output` (paths, criterion, file handle). |
| `regrid.f90` | `ConvertDiscrVert()` — mid-simulation vertical re-discretization (project state from a fine compute grid onto a coarser output grid) gated by `SwDiscrvert==1` | Restore the `SwDiscrvert==1` dispatch site (likely in the output coupling section of `swap_mod.f90`). Migrate `SwDiscrvert/numnodNew/dzNew` from `variables.f90` to `state%mesh` (regrid-temporary fields). The macropore-retired-zero arrays the `swop==2` branch consumes can stay deleted — `swop==2` is the macropore path which is permanently retired (ADR 0040). |
| `sptabulated.f90` | `module doTSPACK` + `module TSPACK` (Renka tension-spline library, ~3,370 lines) + `EvalTabulatedFunction` + `PreProcTabulatedFunction` — the entire `swsophy=1` tabulated soil-hydraulics path. All 6 regression cases run `swsophy=0` (analytical MvG); the tabulated branch had no live dispatch site since the TOML migration. | Restore the 5 `call EvalTabulatedFunction` invocations in `src/utils/soilhydraulicsutils.f90` (`watcon`, `moiscap`, `dhconduc`, `hconduc`, `prhead`) — signatures preserved as comments at each `fatalerr_collected` stub. Wire `numtablay`/`ientrytablay`/`sptablay` through a new `soil_hydraulics_tab_t` sub-record on `state%cfg%soil` (no config home today — they were read by the retired `.swp` reader). Expand layer→node in `config_to_variables` or `soilwater_init`, writing directly to `state%soilwater%{numtab, ientrytab, sptab}`. Restore the `PreProcTabulatedFunction` dispatch (no live caller today either). Decide whether `do_ln_trans` (currently a compile-time parameter in `src/soil/do_ln_trans.f90`) should become a runtime switch. Add this file to both `meson.build` and `tests/unit/meson.build`. |

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
