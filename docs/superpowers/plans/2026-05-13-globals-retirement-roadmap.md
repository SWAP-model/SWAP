---
title: "Globals Retirement — Subsystem-by-Subsystem Roadmap"
date: 2026-05-13
status: draft
context: post-SS-BMI2 — retire 33 use variables consumers, then delete config_to_variables + variables.f90 + initialize.f90
---

# Globals Retirement Roadmap

## Why now

Three foundational arcs shipped this week — SS-DRV, SS-TCM, SS-BMI2 — and the working prototype is byte-for-byte regression-clean. The remaining tech-debt visible in `swap_mod.f90`'s `swap_init_from_loaded_config` (transient `*_init_buf` reads, `bind_*_target` pointer-binding calls, inline swinco=3 seeds) is *symptomatic*. The underlying disease:

- **33 files** still `use variables`.
- **~300+ distinct symbols** flow through that legacy bare-globals interface.
- **4 038 lines** of legacy adapter code keep the bridge alive (`config_to_variables.f90` 1802 + `variables.f90` 1391 + `initialize.f90` 845).
- The adapter populates bare globals → other code reads from bare globals → `config_to_variables` stays load-bearing → transient buffers exist to bridge it → the seeds in `swap_mod` persist.

Migrating each reader cluster to read directly from `state` / `config` slices unwinds the whole chain. After all readers move:
- `config_to_variables` becomes dead code (delete).
- `variables.f90` declarations have no consumers (delete).
- `initialize.f90` zero-fills nothing (delete).
- Transient buffers in `config_to_variables` retire as a natural side effect.
- The strangler-fig leftovers in `swap_mod` retire as a natural side effect.

Net deletion target: ~4 000 lines of legacy plumbing.

## Strategy: subsystem-by-subsystem reader migration, then final retirement arc

8 reader-migration arcs (one per source cluster) + 1 final retirement arc. Each reader arc is sized like SS-DRV / SS-TCM / SS-BMI2 — 1 to 3 days of focused work. Total estimated work: 2–3 weeks.

## Pattern: what each reader-migration arc does

For one source cluster (e.g. `src/heat/`):

1. **Inventory.** For each file in the cluster, run:
   ```bash
   grep -A 10 "^\s*use variables" <file>
   ```
   List every symbol imported (whether via `only:` or bare `use variables`).
2. **Classify each symbol.** Each falls into one of:
   - **Belongs on `state%<X>`** (runtime mutable data) → already migrated or migrates as part of this arc (e.g. `tsoil` → `state%heat%tsoil`).
   - **Belongs on `config%<X>`** (read-only after init) → reader takes `config` as arg (or pulls from `state%<X>` if the state-init-pilot pattern has already cached the config value into state).
   - **Genuine cross-subsystem global** (`flMacroPore`, `dt_SSDI_event`) → migrate to a more appropriate state subrecord, or leave behind for the final retirement arc to address with a per-symbol fix.
3. **Migrate each reader.** Drop the symbol from `use variables, only:`, replace the bare-name reads with `state%<X>%<symbol>` or `config%<X>%<symbol>`. Where state is not in scope, thread it through the signature.
4. **Confirm `use variables` removed entirely** from each file in the cluster (the goal is `grep -n "use variables" <file>` returns empty).
5. **Verification gates per memory `feedback_per_task_regression_gate.md`:**
   - `pixi run build-linux` (clean rebuild when state schema touched, per `feedback_state_schema_clean_rebuild.md`)
   - `pixi run test-pfunit`
   - `pixi run -e test python tests/regression/test_output_regression.py hupselbrook surfacewater salinitystress grassgrowth`
   - Byte-for-byte 4/4 regression preserved.

When a needed field doesn't exist on `state%<X>` yet (legacy global that hasn't been migrated), the arc has two options:
- **Migrate the field into state as part of this arc** (small SS-* style schema addition + dual-write transition + reader cutover + global retirement, all in one). Best when the field is owned by this cluster.
- **Defer the field** to a follow-on arc; the reader continues importing it via `use variables` for now. Acceptable for cross-subsystem fields that don't naturally belong to this cluster.

## Arcs in suggested order

### Arc 1 — `utils/` (small, focused pilot)
**Files:** `src/utils/soilhydraulicsutils.f90`, `src/utils/surfacewaterutils.f90`. ~10 symbols total via clean `only:` clauses. **Estimated effort: 1–2 days.**

**Why first:** smallest scope, cleanest existing imports, validates the pattern after surfacewater_state pilot. **Note:** this arc does NOT retire the `bind_*_target` pattern — that's a deeper architectural cleanup that lands when the state pointers can be threaded as explicit arguments (deferred to Arc 9 retirement or a later atomic fix).

### Arc 2 — `boundary/` (Richards-equation hot path)
**Files:** `src/boundary/boundbottom.f90`, `src/boundary/boundtop.f90`. Bare `use variables` (no `only:`). Read `dtmin`, `dtmax` (already migrated to state%timecontrol), `h`, `pondmx`, etc. **Estimated effort: 1–2 days.**

**Risk:** these are called many times per timestep — regression is sensitive to any solver-input mistake. Mitigated by clean rebuild + check-fast after every commit.

### Arc 3 — `heat/` (mostly self-contained)
**Files:** `src/heat/frozencond.f90`, `src/heat/temperature.f90`. Bare `use variables`. Gated by `swfrost==1` so risk on the regression suite is concentrated in oxygenstress / salinitystress. **Estimated effort: 1–2 days.**

`heat_init` lives in `src/heat/temperature.f90` rather than `src/state/heat_state.f90` — this arc relocates it as a type-bound `state%heat%init(config%heat, dims)` per the surfacewater_state pilot pattern.

### Arc 4 — `atmosphere/` (medium-size, meteo-heavy)
**Files:** `src/atmosphere/et.f90`, `src/atmosphere/interception.f90`, `src/atmosphere/meteoday.f90`, `src/atmosphere/meteodt.f90`. **Estimated effort: 2–3 days.**

Includes the swinco=3 atmosphere warm-restart seeding currently inline in `swap_mod`. After this arc, `swap_mod` calls `state%atmosphere%init(config%soil, config%meteo)` and the inline block is deleted.

### Arc 5 — `soil/` (depends on utils being done)
**Files:** `src/soil/soilgrid.f90`, `src/soil/soilhydraulics.f90`, `src/soil/waterbalance.f90`, `src/soil/WC_K_models_04_11.f90`. **Estimated effort: 2–3 days.**

Soil hydraulics calls soilhydraulicsutils — utils must be migrated first so the soil-side migration can read state/config consistently. `WC_K_models_04_11.f90` holds the `bind_cofgen_target` pointer pattern — this arc may retire it by threading state directly.

### Arc 6 — `drainage/` (medium)
**Files:** `src/drainage/divdra.f90`, `src/drainage/drainage.f90`, `src/drainage/surfacewater.f90`. **Estimated effort: 2 days.**

Drainage already has `drainage_init(state, config)` and surfacewater has `state%surfacewater%init` from the pilot. This arc cleans up the runtime readers (`SurfaceWater(task=2|3)`, `DIVDRA`).

### Arc 7 — `io/` (output formatters)
**Files:** `src/io/readmeteo.f90`, `src/io/swap_csv_output.f90`, `src/io/swapoutput.f90`. **Estimated effort: 2–3 days.**

Output formatters read many config flags (output cadence, format switches). Most of `state%timecontrol` already covers what they need (the SS-BMI2 sweep). This arc finishes by migrating any leftover legacy globals these files still read.

### Arc 8 — `crop/` (biggest, last)
**Files:** `src/crop/cropfixed_init.f90`, `cropgrass_init.f90`, `cropwofost_init.f90`, `cropgrowth.f90`, `irrigation.f90`, `oxygenstress.f90`, `rootextraction.f90`, `tillage.f90`, plus `wofost_soil_parameters.f90`. **Estimated effort: 3–5 days.**

The crop subsystem is the biggest and most coupled. May need internal decomposition: separate sub-arcs for crop-init vs. crop-runtime, or per crop type (fixed / grass / wofost). Decide during the brainstorm for this arc.

**Note:** there is no `crop_state_t` today — crop runtime fields live partly in `state%soilwater`, partly as bare globals in `variables.f90`. A `crop_state_t` may need to be designed here (or an "amenities" arc precedes it). This is the only arc that introduces a brand-new state subrecord.

### Arc 9 — Final retirement (the payoff)
**Files:** `src/io/toml/config_to_variables.f90` (delete, ~1802 lines), `src/core/variables.f90` (delete, ~1391 lines), `src/core/initialize.f90` (delete, ~845 lines), `src/core/swap_mod.f90` (drop transient-buffer reads + bind_*_target calls + swinco=3 inline block), `meson.build` (drop deleted files). **Estimated effort: 1–2 days.**

Preconditions verified by grep:
- `grep -rn "use variables" src/ --include="*.f90"` returns empty (every reader migrated).
- `grep -rn "tc_iyear_init_buf\|h_init_buf\|pondini_init_buf\|pond_init_buf" src/` returns no live readers.
- `grep -rn "bind_cofgen_target\|bind_state_targets\|bind_tc_target" src/` returns no live callers.

After this arc: `swap_init_from_loaded_config` collapses to ~30 lines of clean S4 orchestration (state%X%init calls in dependency order, plus the Phase-1 physics calls). The strangler-fig is gone.

## Suggested dependency ordering

```
Arc 1 (utils) ──┐
                ├─→ Arc 5 (soil)  ──┐
                │                   ├─→ Arc 8 (crop) ──→ Arc 9 (retirement)
Arc 2 (boundary)│                   │
Arc 3 (heat)    ├──→ Arc 4 (atmosphere) ──→ Arc 6 (drainage)
Arc 7 (io)──────┘                                       
```

Arcs 1, 2, 3, 7 can run in any order (mostly independent). Arc 4 (atmosphere) is best after at least one of those is done so the state-init-pilot pattern has another data point. Arc 5 (soil) follows utils. Arc 6 (drainage) follows atmosphere. Arc 8 (crop) is last because it's the biggest and most coupled, possibly introducing `crop_state_t`. Arc 9 closes the loop.

## Verification discipline (applies to every arc)

Per memory `feedback_per_task_regression_gate.md` + `feedback_state_schema_clean_rebuild.md`:

1. **Every implementer subagent** ends with `rm -rf builddir && pixi run build-linux && pixi run test-pfunit && pixi run -e test python tests/regression/test_output_regression.py hupselbrook surfacewater salinitystress grassgrowth` (4/4 byte-for-byte required).
2. **Every arc** ends with `pixi run check-full` (5/5).
3. **Each commit migrates one file or one tight cohort.** Clean bisection lane preserved.

## Memory & ADR consequences

When arc 9 ships:
- Memory `feedback_state_schema_clean_rebuild.md` becomes obsolete IF the Meson cross-static-library `.mod` deps are also fixed in arc 9 (see BMI roadmap A1). Retire the memory then.
- Memory `feedback_per_task_regression_gate.md` stays.
- ADR candidate: a new ADR documenting the post-cleanup state architecture — "all state through `swap_state_t`, all config through `swap_config_t`, no bare globals." Worth writing as the cap on the arc.

## What this arc DOES NOT change

- BMI / cffi surfaces (`swap_bmi_mod`, `swap_capi_mod`) — they keep working throughout. They depend only on `state` and `config`, which is exactly what this arc cements.
- The Python demo `tests/cffi-demo/run_ensemble.py` — unchanged.
- Physics — every value computed by every formula stays bit-identical. Only data-access paths change.

## Files & artifacts inventory

**This roadmap:** `docs/superpowers/plans/2026-05-13-globals-retirement-roadmap.md` (this file).

**Per-arc deliverables** (when picked up):
- Spec → `docs/superpowers/specs/YYYY-MM-DD-globals-<cluster>-design.md`
- Plan → `docs/superpowers/plans/YYYY-MM-DD-globals-<cluster>.md`

**State init pilot precedent:**
- Spec: `docs/superpowers/specs/2026-05-13-state-init-pilot-surfacewater-design.md`
- Implementation: `src/state/surfacewater_state.f90` (the `init` type-bound procedure)
- Driver call: `src/core/swap_mod.f90:172` (`call state%surfacewater%init(config%surface_water, config%drain, numnod)`)

**Adapter reference (for arc 9):**
- `src/io/toml/config_to_variables.f90` — 1802 lines, target of retirement.
- `src/core/variables.f90` — 1391 lines, target of retirement.
- `src/core/initialize.f90` — 845 lines, target of retirement.
