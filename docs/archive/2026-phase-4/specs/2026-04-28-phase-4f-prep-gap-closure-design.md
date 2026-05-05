---
title: "Phase 4f-prep — Schema Gap Closure (Top-N)"
author: Mateusz Zawadzki
date: 2026-04-28
status: draft
---

# Phase 4f-prep: Schema Gap Closure (Top-N)

A focused phase that closes the highest-priority real-domain (G) gaps from the Phase 4e variables audit, sized to unblock Phase 4f's strangler-fig of `readswap()`. Phase 4f-prep delivers two new typed configs (`macropore_config_t`, `surface_water_config_t`) and extensions to five existing configs covering the top-10 blockers from the audit. End-state: real-domain G count drops from 118 to ≤30; Phase 4f opens with high confidence the new schema covers what physics actually reads.

## Background

Phase 4e produced [`docs/phase-4f-config-to-variables-audit.md`](../../phase-4f-config-to-variables-audit.md). The audit classified every `variables%` field referenced by execution paths as **C** (covered by `swap_config_t`), **R** (runtime state), **G** (gap — read from input but no config slot), or **RETIRED** (legacy output switches deprecated per ADR 0009). The post-Phase-4e count is 203 C / 806 R / 189 G / 18 RETIRED. Of the 189 G, 71 are in a heuristic "Other / uncategorised" overflow bucket; the remaining 118 are real-domain gaps in drainage, soil, meteo, crop, time/control, macropore, bottom boundary, heat, irrigation, and solute.

The Phase 4f strangler-fig replaces `call ReadSwap()` in `src/core/swap.f90:124` with a TOML load + `config_to_variables(config)` adapter. That adapter copies every (C) field; (R) fields stay alone; **(G) fields are blockers** — the adapter has nothing to copy from. Phase 4f-prep closes enough G gaps that the strangler-fig can land without runtime "uninitialized field" surprises.

Per the brainstorming session: scope is the **top-N priority subset** (Q1 option A) — close the audit's top-10 plus their direct neighbors, not the full long tail. Macropore coverage is **full** (Q2 option B) — all 22 keys, one full new `macropore_config_t`. Per-case TOMLs and parity tests are extended **only on cases that exercise the new section** (Q3 option A) — symmetry across cases is not a goal.

## Scope

### In

- **Two new typed configs**: `macropore_config_t` (full 22-key coverage) and `surface_water_config_t` (covers `swsrf`/`swsec` cluster + surface-runoff sub-block). Both Phase 4d-deferred.
- **Extensions to five existing configs**:
  - `simulation_config_t` — `[simulation.numerical]` sub-section: `dt`, `dtmin`, `dtmax`, `MaxIt`, `MaxBackTr`, `taccur`.
  - `meteorology_config_t` — `[meteorology.snow]` (`swsnow` + sub-fields) + `[meteorology.rainfall]` extension (`cnref`, `cndry`, `cnwet`, `cofred`, `cfeic`, `cfevappond`).
  - `soil_config_t` — `[soil.discretization]` (`swdiscrvert`, `dznew(:)`, `numnodnew`); `[soil.anisotropy]` (`cofani(:)` per layer); `[soil.frost]` (`swfrost`, `tfroststa`, `tfrostend`, `swsublim`); `nrstaring` top-level.
  - `drainage_config_t` — `[drainage.surface_runoff]` sub-section (`swnrsrf`, `SwTopnrsrf`, `swdivdinf`, `FacDpthInf`).
  - `crop_config_t` — `[[crop.rotation]].swhydrlift`. Verify `croptype` audit-flag is a false positive for the existing `rotation_type(:)` field; if real, address.
- **Per-case TOML extensions** only where the section is exercised by the case's legacy `.swp`. Macropore section → only case 3. Surface-water section → only case 6. Snow/frost → wherever those switches are 1 in legacy. Numerical sub-section → wherever non-default values are set.
- **Parity test extensions** only on the same cases as TOML extensions.
- **Tests-first per task** — each new field/section gets validator unit tests and reader unit tests before per-case parity assertions are added.
- **Closeout**: schema doc bump, coverage rebaseline, regression baseline log, tag `rescue/phase-4f-prep-gap-closure`.

### Out

- **Resolving the 71 "Other / uncategorised" G entries.** Those need re-categorization more than new fields. Stay flagged for opportunistic resolution during Phase 4f if any cause runtime issues.
- **Wiring the new pipeline into `swap.f90`.** That's Phase 4f.
- **Deletion of legacy output writer code paths** per ADR 0009. Phase 5+.
- **Macropore-specific physics changes.** Phase 4f-prep is purely additive to the schema — physics paths stay on the legacy reader until Phase 4f.

## Critical decisions

### D1 — `macropore_config_t` shape (full, 22 keys)

Per the macroporeflow audit, the keys split into:
- **13 scalars** (switch + thresholds): `swmacro`, `z_ah`, `z_ic`, `z_st`, `vlmpstss`, `ppicss`, `numsbdm`, `powm`, `rzah`, `spoint`, `swpowm`, `dipomi`, `dipoma`.
- **9 per-layer table fields**: `swsoilshr(:)`, `swshrinp(:)`, `thetcrmp(:)`, `geomfac(:)`, `shrpar(:,:)` (4 cols), `swsorp(:)`, `sorpfacparl(:)`, `sorpmax(:)`, `sorpalfa(:)`.

The per-layer fields are sized by `nlayers` from `soil_config_t`. The `shrpar` matrix has 4 columns (a, b, c, d coefficients) per layer. Validator checks: array sizes consistent with `nlayers`; switch-gated bounds when `swmacro=1`.

Top-level TOML section: `[macropore]`. No cross-file reference (the macropore block in legacy `.swp` is inline; we mirror).

### D2 — `surface_water_config_t` shape

The audit's drainage section flagged 20 G entries; ~10 are surface-water-specific. Initial fields: `swsrf`, `swsec`, `nmper`, `dropr`, `wlsman`, `gwlcrit`, `nphase`, `wscap`, `wldip`, plus `RapDra*` table.

This is genuinely surface-water management (drainage by surface flow + weir control), not the same as `[drainage]`'s subsurface drainage. Separate type, separate top-level section: `[surface_water]`.

### D3 — Switch-gated validators

Most new fields are only required when their gating switch is 1. The validator pattern from Phase 4d (e.g., `bottom_boundary_config_t` per-`swbotb`-branch) is the template:

```fortran
if (self%swmacro == 1) then
   call check_real_range(self%z_ah, -100.0d0, 0.0d0, 'macropore.z_ah', errors)
   ! ... etc ...
end if
```

When the switch is 0, validators skip. This matches Q3 option A's behavior — cases that don't exercise the section have no `[macropore]` block in their TOML; the reader leaves `swmacro = 0` (default); the validator skips its body entirely.

### D4 — `croptype` audit flag verification

The audit flagged `croptype` as a top-10 G entry but `crop_config_t` already has `rotation_type(:)`. This is likely a false positive (auditor's exact-name match missed the alias). **Task A0 verifies this** before any implementation: grep for `croptype` usage; if it's solely an alias for `rotation_type`, document and skip; if it's a genuinely distinct field, add it.

### D5 — `swhydrlift` per-rotation-entry

The legacy crop sub-readers (`readwofost`, `readcropfixed`, `readgrass`) take a `swhydrlift` argument. The Phase 4c-b legacy_crop_helper hardcodes it to 0. The `[[crop.rotation]]` array entries should expose this as an optional per-entry field defaulting to 0.

### D6 — Test ordering (per user note)

The user explicitly asked for tests that "make work easier" — i.e., tests that catch bugs early, not just coverage padding. Three layers:

1. **Unit tests** for each new validator stanza and reader function. **Authored before** the implementation (TDD-style for the small additions; not full TDD for the large schema additions, but at least the validator unit tests land in the same commit as the validator code).
2. **Per-case parity tests** when the new section is exercised by a regression case. These catch "schema disagrees with legacy `.swp` semantics" issues early.
3. **`pixi run -e test check-fast`** after each task's commit, **`check-full` after each Part**. The 6/6 regression suite is the integration safety net — physics behavior must not change.

The standing parity-helper pattern (`reset_for_next_readswap()` + `load_both_for_<case>`) carries forward unchanged.

## Critical files

### Create

| File | Responsibility |
|---|---|
| `src/config/macropore_config.f90` | New `macropore_config_t` type + validators + finalize. |
| `src/config/surface_water_config.f90` | New `surface_water_config_t` type. |
| `src/io/toml/read_macropore_toml.f90` | Section reader. |
| `src/io/toml/read_surface_water_toml.f90` | Section reader. |
| `tests/unit/config/test_macropore_config.pf` | Validator unit tests. |
| `tests/unit/config/test_surface_water_config.pf` | Validator unit tests. |
| `tests/unit/io/toml/test_read_macropore_toml.pf` | Reader unit tests. |
| `tests/unit/io/toml/test_read_surface_water_toml.pf` | Reader unit tests. |
| `tests/regression/baselines/phase-4f-prep-gap-closure.log` | Closeout regression log. |

### Modify

| File | Change |
|---|---|
| `src/config/swap_config.f90` | Add `macropore` and `surface_water` sub-fields. |
| `src/config/simulation_config.f90` | Add numerical sub-section + validators + tests. |
| `src/config/meteorology_config.f90` | Add snow + rainfall extensions + tests. |
| `src/config/soil_config.f90` | Add discretization + anisotropy + frost + nrstaring + tests. |
| `src/config/drainage_config.f90` | Add surface_runoff sub-section + tests. |
| `src/config/crop_config.f90` | Add `swhydrlift` per-rotation-entry + tests. |
| `src/io/toml/load_swap_config.f90` | Wire the two new section readers + extension reads. |
| `src/io/toml/read_simulation_toml.f90` | Read new numerical sub-section. |
| `src/io/toml/read_meteorology_toml.f90` | Read new snow + extended rainfall. |
| `src/io/toml/read_soil_toml.f90` | Read new sub-sections. |
| `src/io/toml/read_drainage_toml.f90` | Read new surface_runoff sub-section. |
| `src/io/toml/read_crop_toml.f90` | Read per-rotation `swhydrlift`. |
| `tests/unit/config/test_*_config.pf` | Extend with new field tests for each modified config. |
| `tests/unit/io/toml/test_read_*_toml.pf` | Extend with new field reader tests. |
| `tests/unit/io/toml/test_*_parity.pf` | Per-case parity assertions where the section is exercised. |
| `tests/swap-cases/toml/<case>/swap.toml` | Per-case TOML extensions where exercised (submodule). |
| `tests/unit/meson.build`, `tests/unit/testSuites.inc` | Register new test files. |
| `meson.build` | Register new source files. |
| `docs/configuration-schema.md` | Document new sections + extensions. |
| `docs/phase-4f-config-to-variables-audit.md` | Reclassify closed gaps from G to C. |
| `docs/coverage-baseline.md` | Phase 4f-prep section. |

## Phase 4f sketch (informational; detailed plan after 4f-prep closes)

Phase 4f opens with the variables audit at ≤30 G entries (real-domain). The strangler-fig is then:

1. Author `src/io/toml/config_to_variables.f90` with a single subroutine that copies every (C) field plus the (RETIRED) zero-forcing per ADR 0009.
2. Replace `call ReadSwap()` in `src/core/swap.f90:124` with TOML load + `config_to_variables(config)` call.
3. Verify check-full 6/6 green via the new path.
4. `git rm src/io/readswap.f90`; update meson.
5. Tag `rescue/phase-4f-strangler-swp`.

## Verification gates

After each Part: `pixi run -e test test-pfunit` clean; `pixi run -e test check-fast` 4/4 green.
After each Part of meaningful size: `pixi run -e test check-full` 6/6 green.
Final: `git diff rescue/phase-4e-error-prep..HEAD -- src/` shows ONLY new/modified config files + readers + load_swap_config wiring + tests. No physics-algorithm changes.

End-state success: real-domain G count in `docs/phase-4f-config-to-variables-audit.md` is ≤30. Audit doc reclassification commit lands as part of closeout.

## File count summary

- Create: ~10 new source/test files (2 config types + 2 readers + 4 unit-test files + 2 misc).
- Modify: ~20 existing source/test files (config extensions, reader extensions, parity test extensions, audit + schema docs).
- Plus per-case TOML edits in submodule (1-2 cases per new section, ~6-10 file edits across the submodule).
- LoC delta: estimated +2,500 / -100 (mostly additive; small deletions in audit doc as G entries become C).
