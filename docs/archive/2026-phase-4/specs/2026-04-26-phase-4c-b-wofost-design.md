---
title: "Phase 4c-b — WOFOST Crop (Type 2) + Salinitystress Case"
author: Mateusz Zawadzki
date: 2026-04-26
status: draft
---

# Phase 4c-b: WOFOST Crop (Type 2) + Salinitystress Case

The second half of Phase 4c. Closes the crop-rotation matrix opened by 4c-a by adding the WOFOST type 2 reader, completing hupselbrook's potatod entry, and bringing 5.salinitystress to full parity end-to-end. After 4c-b, all five non-macropore cases have complete `.swp` + `.dra` + `.crp` parity against the legacy `variables` globals.

## Background

4c-a delivered cross-file TOML loading and types 1 (fixed) and 3 (grass). Hupselbrook's `potatod.crp` (type 2, WOFOST) was stubbed with a placeholder `potatod.crp.toml`; salinitystress was deferred entirely. The dispatch in `read_crop_toml.f90` already has a `select case (rotation_type(i))` block; 4c-b adds the `case (2)` arm and the matching reader.

4c-a's crop parity used **hardcoded values** extracted from the `.crp` files (Phase 4b "Option B"). The reason: the legacy crop sub-readers — `readwofost`, `readgrass`, `readcropfixed` in `readswap.f90` — are called from `cropgrowth.f90` per-rotation-entry during simulation, *not* from the top-level `readswap()` invoked in pFUnit. So `variables%dtsmtb`, `variables%hlim1`, etc. stay at zero after `readswap()` returns. Hardcoded extraction proved "TOML matches the .crp on paper" but did not exercise the legacy reader.

For 4c-b we upgrade. `readswap.f90` is a flat file of free subroutines (no module wrapper), and `readwofost(icrop, crpfil, swhydrlift, swsoybean, mg, dvsi, dvrmax1, dvrmax2, flrfphotoveg, …)` at line 2520 is callable directly from any code unit, including pFUnit. The 4c-b parity test therefore calls `readwofost()` from the test body to populate `variables%` globals, then loads the same crop file via TOML, and compares actually-read values on both sides. Every WOFOST scalar and table goes through real legacy-vs-new parity. No hardcoded shortcut for type 2.

The same upgrade applies retroactively to types 1 and 3 — `readcropfixed` (line 2037) and `readgrass` (line 3437) are equally callable. 4c-b will tighten 4c-a's hardcoded crop assertions to legacy-reader-driven assertions for hupselbrook's maizes/grassd entries and the grassgrowth/oxygenstress/surfacewater cases as a closeout step.

## Scope

### In

- **`cropwofost_config_t`** in `src/config/cropwofost_config.f90`. Single top-level type composed of small grouped sub-types (Option B, user-confirmed), all field names mirroring the legacy `.crp` variable names verbatim:
  - `wofost_preparation_t` — `swprep`, `zprep`, `hprep`, `maxprepdelay`
  - `wofost_sowing_t` — `swsow`, `zsow`, `hsow`, `ztempsow`, `tempsow`, `maxsowdelay`
  - `wofost_germination_t` — `swgerm`, `tsumemeopt`, `tbasem`, `teffmx`, `hdrygerm`, `hwetgerm`, `zgerm`, `agerm`
  - `wofost_harvest_t` — `dvsend`, `swharv`
  - `wofost_cropfactor_t` — `swcf`, `albedo`, `rsc`, `rsw`, plus `cftb(:,:)` / `chtb(:,:)` (CF or CH vs DVS)
  - `wofost_phenology_t` — `idsl`, `tsumea`, `tsumam`, `dlo`, `dlc`, `vernsat`, `vernbase`, `verndvs`, plus tables `dtsmtb(:,:)`, `verntb(:,:)`
  - `wofost_initial_t` — `tdwi`, `laiem`, `rgrlai`
  - `wofost_greenarea_t` — `spa`, `ssa`, `span`, `tbase`, plus `slatb(:,:)`
  - `wofost_assimilation_t` — `kdif`, `kdir`, `eff`, plus tables `amaxtb(:,:)`, `tmpftb(:,:)`, `tmnftb(:,:)`
  - `wofost_conversion_t` — `cvl`, `cvo`, `cvr`, `cvs`
  - `wofost_respiration_t` — `q10`, `rml`, `rmo`, `rmr`, `rms`, plus `rfsetb(:,:)`
  - `wofost_partitioning_t` — `frtb(:,:)`, `fltb(:,:)`, `fstb(:,:)`, `fotb(:,:)`
  - `wofost_death_t` — `perdl`, plus `rdrrtb(:,:)`, `rdrstb(:,:)`
  - `wofost_root_t` — `swrd`, `rdi`, `rri`, `rdc`, `swdmi2rd`, `wrtmax`, plus tables `rdtb(:,:)`, `rlwtb(:,:)`, `rdctb(:,:)`
  - `wofost_oxygen_stress_t` — `swoxygen`, `swwrtnonox`, `aeratecrit`, `hlim1`, `hlim2u`, `hlim2l`, `q10_microbial`, `specific_resp_humus`, `srl`, `swrootradius`, `dry_mat_cont_roots`, `air_filled_root_por`, `spec_weight_root_tissue`, `var_a`, `root_radiusO2`
  - `wofost_drought_stress_t` — `swdrought`, `hlim3h`, `hlim3l`, `hlim4`, `adcrh`, `adcrl`
  - `wofost_salinity_t` — `swsalinity`, `saltmax`, `saltslope`, `salthead`
  - `wofost_compensate_t` — `swcompensate`, `swstressor`, `alphacrit`, `dcritrtz`
  - `wofost_interception_t` — `swinter`, `cofab`, plus `gashtb(:,:)` for SWINTER=2 (PFREE/PSTEM/SCANOPY/AVPREC/AVEVAP vs T)
  - `wofost_co2_t` — `swco2`, `atmofil`, plus tables `co2amaxtb(:,:)`, `co2efftb(:,:)`, `co2tratb(:,:)`
  - `wofost_management_t` — `fraharlosorm_lv`, `fraharlosorm_st`, `fraharlosorm_so`, `fradeceasedlvtosoil`, `swpotrelmf`, `relmf`
  - **Out of 4c-b**: nitrogen-use block (LINTUL4 fields, all currently commented in the .crp template; no test case sets them), and the irrigation-scheduling section (Part 1–3 of the irrigation block; this maps to a future `irrigation_config_t` in Phase 4d).
- **`read_cropwofost_toml`** in `src/io/toml/read_cropwofost_toml.f90`. Reads grouped TOML sections matching the sub-type layout (`[preparation]`, `[sowing]`, `[germination]`, `[harvest]`, `[crop_factor]`, `[phenology]`, `[initial]`, `[green_area]`, `[assimilation]`, `[conversion]`, `[respiration]`, `[partitioning]`, `[death]`, `[root]`, `[oxygen_stress]`, `[drought_stress]`, `[salinity]`, `[compensate]`, `[interception]`, `[co2]`, `[management]`). Tables encoded as TOML arrays of arrays (e.g., `dtsmtb = [[0.0, 0.0], [2.0, 0.0], [13.0, 11.0], [30.0, 28.0]]`). Each table is a 2-column or 6-column 2D `real(8)` allocatable. Validators check column count, row count ≤ 15 (legacy AFGEN limit), and monotone first-column where the legacy reader requires it.
- **Crop config dispatch** in `read_crop_toml.f90`: extend the `select case (rotation_type(i))` block with `case (2) ; call read_cropwofost_toml(crp_doc_ptr, config%rotation_wofost(i), errors)`.
- **`crop_config_t` extension**: add `type(cropwofost_config_t), allocatable :: rotation_wofost(:)`. Validate dispatches by type when `rotation_loaded(i)` is true.
- **TOML files in `tests/swap-cases/toml/`**:
  - `1.hupselbrook/potatod.crp.toml` — replaces the 4c-a placeholder; full WOFOST scalars + tables.
  - `5.salinitystress/swap.toml` — new top-level (case has no Phase 4a TOML yet).
  - `5.salinitystress/swap.dra.toml` — new drainage file.
  - `5.salinitystress/potatod.crp.toml` — WOFOST with `swsalinity = 1` to exercise the Maas-Hoffman branch (legacy file uses `swsalinity = 1`).
- **Parity tests**:
  - `test_salinitystress_parity.pf` — new file, ~7 subtests mirroring the 4c-a layout (timing, meteorology, drainage, soil, crop scalars, crop tables). Calls `readswap()` for `.swp`/`.dra` parity, then `readwofost()` for `.crp` parity.
  - `test_hupselbrook_parity.pf` — extend with WOFOST sub-test that calls `readwofost()` for the potatod rotation entry and asserts every scalar + table.
  - **Tighten 4c-a parity**: in `test_hupselbrook_parity.pf` (maizes), `test_grassgrowth_parity.pf`, `test_oxygenstress_parity.pf`, `test_surfacewater_parity.pf`, replace hardcoded crop-scalar assertions with `readcropfixed()` / `readgrass()` calls, then assert `variables%foo` against `config%crop%rotation_fixed(i)%foo` / `rotation_grass(i)%foo`. Hardcoded values become test fixtures only when the legacy reader doesn't populate the field (rare, document each).
- **Schema doc**: extend `docs/configuration-schema.md` with the `[[crop.rotation]]` type-2 file shape, table encoding convention, and the full sub-section list.
- **Coverage rebaseline** + tag `rescue/phase-4c-b-wofost`.

### Out

- **Macropore case** (3.macroporeflow) — out for all crop parity work in 4c-x.
- **Bottom-boundary, heat, solute, irrigation configs** — Phase 4d.
- **WOFOST nitrogen-use (LINTUL4)** — none of the regression cases activate it; deferred to a later phase if/when a test case requires it.
- **WOFOST irrigation-scheduling block** — folds into Phase 4d's `irrigation_config_t`.
- **Wiring new TOML path into runtime** — still legacy-readswap-only. Phase 4 proper.

## Critical decisions

### D1 — Sub-config layering (user-confirmed Option B)

`cropwofost_config_t` is a thin envelope holding ~21 grouped sub-types. Each sub-type owns a focused validator. Field names inside each sub-type match legacy `.crp` variable names verbatim — `dtsmtb`, `slatb`, `amaxtb`, `hlim1`, `swoxygen`, etc. — so a reader of the legacy `.crp` can map 1:1 to TOML. The grouping is purely organizational; no field is renamed for "clarity".

Why this matters: the 4c-a precedent had `cropfixed_config_t` and `cropgrass_config_t` as flat structs because each was small enough (~80 and ~110 fields). WOFOST at ~150 scalars + ~16 tables would be too unwieldy flat — validators would be 200+ lines, and the TOML file structure would have no internal grouping signals. The grouped sub-type approach mirrors the legacy `.crp` Part 0–15 sectioning and keeps each validator under ~30 lines.

### D2 — Table encoding in TOML

Legacy AFGEN tables are 2-column (`X Y` rows) with up to 15 records. The interception table is 6-column (`T PFREE PSTEM SCANOPY AVPREC AVEVAP`). All tables encoded as TOML arrays of arrays:

```toml
[phenology]
idsl   = 0
tsumea = 150.0
tsumam = 1550.0
dtsmtb = [
  [ 0.0,  0.0],
  [ 2.0,  0.0],
  [13.0, 11.0],
  [30.0, 28.0],
]
```

Reader allocates `real(8) :: dtsmtb(:,:)` with `(nrows, ncols)`. Validator checks `ncols == 2` (or 6 for `gashtb`), `nrows >= 2`, `nrows <= 15`, and (for tables where the legacy AFGEN demands it) strict monotonicity of column 1.

### D3 — Parity strategy: call legacy `readwofost()` directly

The 4c-a hardcoded-values workaround was a constraint of the call graph (`readswap()` doesn't reach crop sub-readers). In 4c-b we lift that by calling `readwofost(icrop=1, crpfil='potatod', ...)` directly from the parity test body after `readswap()`. The legacy reader populates `variables%` globals; the new TOML reader populates `config%crop%rotation_wofost(i)`. Comparison is real legacy vs real new on every field.

This requires figuring out the right argument set for `readwofost`. The signature has 8 named args after `(icrop, crpfil)` for switches and outputs — the test passes the same values the live caller in `cropgrowth.f90` passes. We will document the full call inside `chdir_helper.f90` or a new `legacy_crop_helper.f90` so each parity test can do `call read_legacy_wofost(icrop, crpfil)` cleanly.

The same helper will expose `read_legacy_cropfixed(icrop, crpfil)` and `read_legacy_grass(icrop, crpfil)` for the 4c-a tightening pass.

### D4 — `rotation_wofost(:)` parallel-array shape (matches 4c-a)

The 4c-a precedent stores per-entry sub-configs in parallel arrays keyed by the rotation index:
```fortran
type(cropfixed_config_t), allocatable :: rotation_fixed(:)
type(cropgrass_config_t), allocatable :: rotation_grass(:)
logical,                  allocatable :: rotation_loaded(:)
```

4c-b adds:
```fortran
type(cropwofost_config_t), allocatable :: rotation_wofost(:)
```

The dispatch in `validate` and `read_crop_toml` extends with a `case (2)` arm. Allocations are sized to `nrot` even though only one slot per type is populated per rotation entry — keeps indexing simple at the cost of a few unused allocatables.

### D5 — Salinitystress case TOML authoring

5.salinitystress has no Phase 4a TOML — 4c-b authors the full set from scratch (`swap.toml`, `swap.dra.toml`, `potatod.crp.toml`). The `.swp` and `.dra` parity follow the 4c-a pattern (extract from legacy files, assert against `variables` globals). The crop file activates `swsalinity = 1` so the Maas-Hoffman validator branch is exercised.

### D6 — 4c-a tightening pass (closeout, not blocker)

After 4c-b's WOFOST work is done and salinitystress is green, we sweep through the four 4c-a tests and replace hardcoded crop-scalar assertions with legacy-reader-driven ones. This is closeout work: it does not gate the 4c-b tag. If the tightening surfaces a parity bug, fix it; if it surfaces a *new* readswap-shape limitation, document and defer to Phase 4d.

## Critical files

### Create

- `src/config/cropwofost_config.f90` — top-level + 21 sub-types, ~1100 LoC.
- `src/io/toml/read_cropwofost_toml.f90` — section reader, ~700 LoC.
- `tests/unit/config/test_cropwofost_config.pf` — type construction, validator coverage.
- `tests/unit/io/toml/test_read_cropwofost_toml.pf` — table parsing, error paths.
- `tests/unit/io/toml/legacy_crop_helper.f90` — wrappers around `readwofost`, `readcropfixed`, `readgrass`.
- `tests/unit/io/toml/test_salinitystress_parity.pf` — full case-5 parity.
- `tests/swap-cases/toml/5.salinitystress/swap.toml` — new (submodule).
- `tests/swap-cases/toml/5.salinitystress/swap.dra.toml` — new (submodule).
- `tests/swap-cases/toml/5.salinitystress/potatod.crp.toml` — new (submodule).
- `tests/swap-cases/toml/1.hupselbrook/potatod.crp.toml` — replaces 4c-a placeholder (submodule).

### Modify

- `src/config/crop_config.f90` — add `rotation_wofost(:)`; extend validate dispatch.
- `src/io/toml/read_crop_toml.f90` — add `case (2)` arm.
- `tests/unit/io/toml/test_hupselbrook_parity.pf` — add WOFOST sub-test; tighten maizes assertions.
- `tests/unit/io/toml/test_grassgrowth_parity.pf` — tighten grass assertions to use legacy reader.
- `tests/unit/io/toml/test_oxygenstress_parity.pf` — tighten grass assertions.
- `tests/unit/io/toml/test_surfacewater_parity.pf` — tighten fixed-crop assertions.
- `tests/unit/meson.build` — register new sources + tests.
- `tests/unit/testSuites.inc` — register new pFUnit suites.
- `meson.build` — register `cropwofost_config.f90` + `read_cropwofost_toml.f90`.
- `docs/configuration-schema.md` — type-2 schema section.

## Tasks (high-level — full breakdown in plan file)

### Part A — WOFOST config + reader

1. Author `cropwofost_config.f90` skeleton with all 21 sub-types and field declarations. Defer validators to Task 2.
2. Author validators per sub-type. Per-section pFUnit tests for the validator surface.
3. Author `read_cropwofost_toml.f90`. Per-section pFUnit tests for happy-path + at least one error path each.
4. Wire into `crop_config.f90` (`rotation_wofost(:)`) and `read_crop_toml.f90` (`case (2)` arm). Validate dispatch.
5. Build clean. `pixi run -e test test-pfunit` green.

### Part B — Hupselbrook potatod entry

6. Author `tests/swap-cases/toml/1.hupselbrook/potatod.crp.toml` (replace placeholder).
7. Author `legacy_crop_helper.f90` with `read_legacy_wofost` wrapper. Verify it populates `variables%` globals as expected via a smoke pFUnit test.
8. Extend `test_hupselbrook_parity.pf` with WOFOST sub-test. Iterate until green.

### Part C — Salinitystress case end-to-end

9. Audit 5.salinitystress legacy files (`.swp`, `.dra`, `.crp`). Document any case-specific quirks.
10. Author `5.salinitystress/swap.toml` + `swap.dra.toml`.
11. Author `5.salinitystress/potatod.crp.toml` with `swsalinity = 1`.
12. Author `test_salinitystress_parity.pf`: timing, meteo, drainage, soil, crop scalars, crop tables. Iterate until green.

### Part D — 4c-a tightening + closeout

13. Extend `legacy_crop_helper.f90` with `read_legacy_cropfixed`, `read_legacy_grass`.
14. Tighten `test_grassgrowth_parity.pf`, `test_oxygenstress_parity.pf`, `test_surfacewater_parity.pf`, and the maizes/grassd entries in `test_hupselbrook_parity.pf` to call legacy readers. Drop hardcoded values.
15. Schema doc update + coverage rebaseline + tag `rescue/phase-4c-b-wofost`.

## Risks & mitigations

| Risk | Mitigation |
|---|---|
| `readwofost` signature requires runtime state already initialized by `readswap` (e.g., `mg`, `dvsi` come from `.swp`-section reads) | Verify by tracing a green run with `gdb` or print statements; if so, the parity test calls `readswap()` first (already does for `.swp` parity) before `readwofost()` |
| AFGEN tables in legacy code use 1D flat storage `dtsmtb(30)` with sentinel-terminated records, not 2D `dtsmtb(:,:)` | Reader unflattens on read; comparison-helper in `legacy_crop_helper` re-flattens for assertions, OR compares row-by-row up to the populated count |
| Salinitystress case has unusual `.swp` settings that don't map to existing config sections | Audit step (Task 9) catches early; either extend the section or document the deferred field for 4d |
| Tightening pass (Task 14) surfaces silent 4c-a parity bugs | Each surfaced bug becomes a focused commit with the legacy-reader assertion as the test; no different from any normal parity tightening |
| WOFOST scalar parity fails on values that legacy normalizes (e.g., unit conversions in `readwofost`) | Each failure is a teaching moment about the legacy reader; if non-trivial, document with a `READING NOTE` comment in the TOML file or the test |

## Verification

After Part A:
```
pixi run -e test test-pfunit            # all green; new tests pass
pixi run -e test check-full             # 6/6 green (no physics changes)
```

After Part B:
```
pixi run -e test test-pfunit            # hupselbrook WOFOST sub-test green
```

After Part C:
```
pixi run -e test test-pfunit            # 5.salinitystress parity green
```

After Part D:
```
pixi run -e test test-pfunit            # all 5 non-macropore cases full parity
pixi run -e test check-full             # 6/6 green
pixi run -e coverage coverage-report    # rebaselined
```

Final tag: `rescue/phase-4c-b-wofost`.

## File count summary

- Create: ~10 files (`cropwofost_config.f90`, `read_cropwofost_toml.f90`, helper, 2 test_*.pf, 4 TOML files in submodule).
- Modify: ~10 files (crop_config.f90, read_crop_toml.f90, 4 parity tests, meson + testSuites).
- LoC estimate: ~2500 lines added (mostly mechanical: field declarations, validator stanzas, TOML readers, parity assertions).

## Phase 4d preview (informational; separate spec to follow)

Phase 4d will add the remaining configs needed before retiring `readswap.f90`:

- `bottom_boundary_config_t` — switches and fields under `[bottom_boundary]`.
- `heat_config_t` — `swhea`, frost params, soil-temperature initialization.
- `irrigation_config_t` — fixed-irrigation table from `.swp`, plus the WOFOST-side scheduling block (TCS/DCS, thresholds, depth tables) factored out of 4c-b's deferred section.
- `solute_config_t` — `swsolu`, decomposition, root uptake, dispersion.
- Per-event `mowing` and `grazing` tables for grass type 3 (deferred from 4c-a's `cropgrass_config_t` `swharv=0` workaround).

Each gets its own config_t + reader + per-case TOML extensions + parity-test extensions. Ordering: alphabetical by config name (bottom_boundary, heat, irrigation, mowing/grazing, solute), one focused commit per config_t + reader + tests, then per-case TOML wiring + parity extensions at the end. Final tag: `rescue/phase-4d-remaining-configs`.
