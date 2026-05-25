# Phase 4c-b — WOFOST Crop (Type 2) + Salinitystress Case Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `cropwofost_config_t` (Option B sub-config layering, ~21 grouped sub-types) + `read_cropwofost_toml`. Bring 1.hupselbrook's potatod entry and 5.salinitystress to full parity using `readwofost()` called directly from pFUnit. Tighten 4c-a's hardcoded crop assertions across all four prior parity tests by calling `readcropfixed()` / `readgrass()` directly.

**Architecture:** Same composite reader pattern. New `cropwofost_config_t` is a top-level type holding 21 grouped sub-types — every field name verbatim from legacy `.crp` (no renames). Tables encoded as TOML arrays of arrays; reader stores `(:,:)` 2D; parity helper converts to/from legacy flat `(30)` storage. New `legacy_crop_helper.f90` exposes thin wrappers around `readwofost`, `readcropfixed`, `readgrass` for use by parity tests.

**Tech Stack:** Unchanged from 4c-a — gfortran 2008, pFUnit 4.15, meson + pixi, toml-f, `iso_c_binding` for chdir.

Spec: `docs/superpowers/specs/2026-04-26-phase-4c-b-wofost-design.md`

---

## Preamble: context every task needs

**Baseline:** Phase 4c-a complete at commit `0ac5c4a`, tag `rescue/phase-4c-a-crop-fixed-grass`. 155 pFUnit tests passing, 6/6 regression green. `swap.toml` cross-file loading active (`[drainage].file`, `[[crop.rotation]].file`). `crop_config_t` holds `rotation_fixed(:)`, `rotation_grass(:)`, `rotation_loaded(:)` parallel arrays with type-2 dispatch as a no-op. The `case (2)` arm in `read_crop_toml.f90` is the wedge 4c-b expands.

**Branch discipline:**
- Work on `development`. No per-task feature branches.
- One commit per task. Subject `<type>(<scope>): <what>`.
- No pushes to origin during 4c-b (matches rescue policy).
- Phase exit: fast-forward `main` to `development` locally; tag `rescue/phase-4c-b-wofost`.

**Working directory:** `/home/zawadzkim/Code/swap` for all tasks unless noted (submodule work is inside `tests/swap-cases/`).

**Conventions** (unchanged from 4a/4b/4c-a):
- Module name = filename stem + `_mod`.
- `implicit none` after module statement; default `private`; explicit `public ::`.
- `use iso_fortran_env, only: real64` for new code; legacy idioms can stay.
- pFUnit suite filename = `tests/unit/<domain>/test_<module>.pf`; `ADD_TEST_SUITE(test_<module>_suite)` in `testSuites.inc`; `'<domain>/test_<module>.pf'` in `pf_files`.
- **No inline `! comments` on `@assert*` lines** — `funitproc` chokes on the `!`.
- Module sources go in `pfunit_extra_sources` if outside `test_base_sources`.

**Verification after each task:**
```
pixi run -e test test-pfunit         # pFUnit must stay green; count notes
pixi run -e test check-fast          # 4/6 regression cases; <90s
```
After Part A and again at closeout, also run:
```
pixi run -e test check-full          # 6/6 regression cases; ~6 min
```

**Two concepts worth understanding before authoring tasks:**

1. **Legacy AFGEN table storage.** Each table (e.g., `dtsmtb`, `slatb`, `amaxtb`) is declared in `src/core/variables.f90` as `real(8) :: dtsmtb(30)` — a flat array. The legacy `rdador` reader fills it with X/Y pairs interleaved (`dtsmtb(1)=x1, dtsmtb(2)=y1, dtsmtb(3)=x2, dtsmtb(4)=y2, …`) and writes the count to a separate `ifnd` integer. The reader stops at the count; subsequent slots are zero. The new reader stores tables as `real(8), allocatable :: dtsmtb(:,:)` with shape `(nrows, 2)`. Parity comparison must convert one shape to the other — done in `legacy_crop_helper.f90` via a `flatten_table(src_2d) result(flat)` helper.

2. **`readwofost` argument signature.** Line 2520 of `readswap.f90`:
   ```
   subroutine readwofost (icrop,crpfil,swhydrlift,swsoybean,mg,dvsi, &
                          dvrmax1,dvrmax2,flrfphotoveg, …)
   ```
   The first two are inputs (`icrop`, `crpfil`); the rest are switches and outputs that the live caller in `cropgrowth.f90` passes from the rotation entry. Task 7 (the helper module) traces the live call site, captures the exact arg list, and exposes a thin wrapper `read_legacy_wofost(icrop, crpfil)` that supplies the rest from the just-loaded `swap_config_t` (or hardcoded defaults that match `cropgrowth.f90` for rotation-entry context). This is the most fragile bit of 4c-b — verify by stepping through a single rotation entry of the live hupselbrook run with print statements before authoring the wrapper.

---

## File structure

### New source files

| File | Responsibility |
|---|---|
| `src/config/cropwofost_config.f90` | `cropwofost_config_t` + 21 grouped sub-types + per-sub-type `validate` + top-level `finalize`. |
| `src/io/toml/read_cropwofost_toml.f90` | `read_cropwofost_toml(doc_root, config, errors)` — populates `cropwofost_config_t` from `.crp.toml` root. |
| `tests/unit/io/toml/legacy_crop_helper.f90` | `read_legacy_wofost(icrop, crpfil)`, `read_legacy_cropfixed(icrop, crpfil)`, `read_legacy_grass(icrop, crpfil)`. Plus `flatten_table(table_2d)` helper. |
| `tests/unit/config/test_cropwofost_config.pf` | Type construction + per-sub-type validator coverage. |
| `tests/unit/io/toml/test_read_cropwofost_toml.pf` | Section reader happy-path + error paths (missing keys, malformed tables, validator failures). |
| `tests/unit/io/toml/test_salinitystress_parity.pf` | Full case-5 parity (timing, meteo, drainage, soil, crop scalars, crop tables). |
| `tests/swap-cases/toml/1.hupselbrook/potatod.crp.toml` | Replaces 4c-a placeholder. (Submodule.) |
| `tests/swap-cases/toml/5.salinitystress/swap.toml` | New top-level. (Submodule.) |
| `tests/swap-cases/toml/5.salinitystress/swap.dra.toml` | New drainage. (Submodule.) |
| `tests/swap-cases/toml/5.salinitystress/potatod.crp.toml` | New WOFOST with `swsalinity = 1`. (Submodule.) |

### Modified source files

| File | Change |
|---|---|
| `src/config/crop_config.f90` | Add `type(cropwofost_config_t), allocatable :: rotation_wofost(:)`. Extend validate dispatch with `case (2)`. |
| `src/io/toml/read_crop_toml.f90` | Replace `case (2)` no-op with `call read_cropwofost_toml(crp_doc_ptr, config%rotation_wofost(i), errors)`. |
| `tests/unit/io/toml/test_hupselbrook_parity.pf` | Add WOFOST sub-test for potatod entry; tighten maizes/grassd assertions to use legacy_crop_helper. |
| `tests/unit/io/toml/test_grassgrowth_parity.pf` | Tighten grass assertions to `read_legacy_grass`. |
| `tests/unit/io/toml/test_oxygenstress_parity.pf` | Tighten grass assertions. |
| `tests/unit/io/toml/test_surfacewater_parity.pf` | Tighten fixed-crop assertions. |
| `tests/unit/meson.build` | Register `cropwofost_config.f90`, `read_cropwofost_toml.f90`, `legacy_crop_helper.f90` in `pfunit_extra_sources`; new `.pf` files in `pf_files`. |
| `tests/unit/testSuites.inc` | Register new pFUnit suites. |
| `meson.build` | Register `cropwofost_config.f90` + `read_cropwofost_toml.f90` in `sources`. |
| `docs/configuration-schema.md` | Add type-2 section with full sub-section list and table-encoding example. |

---

## Tasks

### Part A — WOFOST config + reader

- [ ] **Task 1 — `cropwofost_config_t` skeleton (no validators).**
  Author `src/config/cropwofost_config.f90` with the top-level type and all 21 sub-types declared as listed in the spec D1 section. Field names verbatim from legacy `.crp` Part 0–15. Tables declared as `real(8), allocatable :: <name>(:,:)`. No validators yet — only the type structure, default initializers (zero/false), and a top-level `finalize` that calls each sub-type's `finalize` (each sub-type's `finalize` is a no-op for now). Register in `meson.build` `sources` and `tests/unit/meson.build` `pfunit_extra_sources`.
  **Verify:** `pixi run -e test test-pfunit` builds clean; existing 155 tests pass.
  **Commit:** `feat(config): add cropwofost_config_t skeleton (Phase 4c-b Task 1)`.

- [ ] **Task 2 — `cropwofost_config_t` validators.**
  Add `validate` to each sub-type. Each validator checks `int_enum`, `range`, and `ordered_pair` rules per the legacy `.crp` value-range comments. Top-level `validate` dispatches to each sub-type. Plus `tests/unit/config/test_cropwofost_config.pf` with one test per sub-type validator covering happy-path + at least one failure mode.
  **Verify:** new tests pass; existing tests still pass.
  **Commit:** `feat(config): add cropwofost_config_t validators (Phase 4c-b Task 2)`.

- [ ] **Task 3 — `read_cropwofost_toml`.**
  Author `src/io/toml/read_cropwofost_toml.f90`. Reads grouped sections (`[preparation]`, `[sowing]`, `[germination]`, `[harvest]`, `[crop_factor]`, `[phenology]`, `[initial]`, `[green_area]`, `[assimilation]`, `[conversion]`, `[respiration]`, `[partitioning]`, `[death]`, `[root]`, `[oxygen_stress]`, `[drought_stress]`, `[salinity]`, `[compensate]`, `[interception]`, `[co2]`, `[management]`). Use existing `toml_field_helpers` for scalars. Tables read via `get_value(section, "tablename", arr_of_arr)` then unflattened into `(:,:)`. Plus `tests/unit/io/toml/test_read_cropwofost_toml.pf` with happy-path read of a synthetic full WOFOST TOML + error tests for missing required keys and malformed tables.
  **Verify:** new tests pass.
  **Commit:** `feat(io/toml): add read_cropwofost_toml (Phase 4c-b Task 3)`.

- [ ] **Task 4 — Wire dispatch in `crop_config_t` and `read_crop_toml`.**
  Add `type(cropwofost_config_t), allocatable :: rotation_wofost(:)` to `crop_config_t`. Allocate to `nrot` size in the same place `rotation_fixed` and `rotation_grass` are allocated. Extend `crop_config_t%validate` to dispatch by `rotation_type(i) == 2` to `rotation_wofost(i)%validate` when `rotation_loaded(i)`. In `read_crop_toml.f90`, replace the `case (2)` no-op with `call read_cropwofost_toml(crp_doc_ptr, config%rotation_wofost(i), errors)` and set `rotation_loaded(i) = .true.`.
  **Verify:** `test-pfunit` green; `check-fast` 4/4 green.
  **Commit:** `feat(config,io/toml): wire WOFOST dispatch into crop config (Phase 4c-b Task 4)`.

- [ ] **Task 5 — Part A check-full.**
  No code change. Run `pixi run -e test check-full`; verify 6/6 regression green (no physics changes — purely config additions). Commit only the baselines log update.
  **Verify:** `tests/regression/baselines/phase-4c-b-task-5-check-full.log` recorded.
  **Commit:** `docs(baselines): record phase 4c-b Part A check-full output`.

### Part B — Hupselbrook potatod entry

- [ ] **Task 6 — Author `1.hupselbrook/potatod.crp.toml`.**
  Inside `tests/swap-cases/` submodule. Replace the 4c-a placeholder with a full WOFOST TOML matching `tests/swap-cases/1.hupselbrook/potatod.crp` (704-line file). All 21 sub-sections populated. Tables verbatim from the .crp. Bump outer-repo submodule pointer.
  **Verify:** `test_hupselbrook_loads.pf` (4c-a) still passes — proves the new TOML is load-able and the WOFOST reader doesn't choke on the `[[crop.rotation]]` type-2 entry.
  **Commit (submodule):** `feat(toml/hupselbrook): full potatod WOFOST TOML (Phase 4c-b Task 6)`.
  **Commit (outer):** `chore(swap-cases): bump submodule for hupselbrook potatod WOFOST TOML`.

- [ ] **Task 7 — `legacy_crop_helper.f90` with `read_legacy_wofost`.**
  Author `tests/unit/io/toml/legacy_crop_helper.f90`. Module exposes:
  - `subroutine read_legacy_wofost(icrop, crpfil)` — calls `readwofost(icrop, crpfil, …)` with the right args. Trace the live call site in `cropgrowth.f90` to determine the exact arg values for a rotation-entry call.
  - `subroutine read_legacy_cropfixed(icrop, crpfil)` and `subroutine read_legacy_grass(icrop, crpfil)` — stubs for now (Part D fills them in).
  - `pure function flatten_table(table_2d) result(flat)` — converts `(:,:)` to flat `(:)` interleaved X/Y so parity tests can compare against legacy `dtsmtb(30)` shape.
  Register in `pfunit_extra_sources`. Add a smoke pFUnit test `tests/unit/io/toml/test_legacy_crop_helper.pf` that calls `read_legacy_wofost(1, 'potatod')` after staging `swap_linux.swp.template` and `readswap()`, and asserts that `variables%dtsmtb(1)` is non-zero (proves legacy reader fired).
  **Verify:** smoke test passes.
  **Commit:** `feat(test/io/toml): add legacy_crop_helper for crop-reader parity (Phase 4c-b Task 7)`.

- [ ] **Task 8 — Extend `test_hupselbrook_parity.pf` with WOFOST sub-test.**
  Add a new `@test` subroutine `test_hupselbrook_parity_potatod_wofost`. Pattern:
  1. Stage `.swp.template` → `.swp` and `chdir` into `tests/swap-cases/1.hupselbrook`.
  2. `call readswap()` to populate `.swp`/`.dra`-side globals.
  3. `call read_legacy_wofost(icrop=3, crpfil='potatod')` (rotation entry index TBD — verify which entry potatod is in hupselbrook's rotation table).
  4. `chdir` back; `call load_swap_config('tests/swap-cases/toml/1.hupselbrook/swap.toml', config, errors)`.
  5. `@assertEqual` on every WOFOST scalar (~80 fields) and every table (compare via `flatten_table` helper).
  Iterate until green: TOML edits, reader bugs, validator gaps, missing fields all surface here.
  **Verify:** new sub-test green; existing 4c-a hupselbrook sub-tests still green.
  **Commit:** `test(parity): hupselbrook potatod WOFOST full parity (Phase 4c-b Task 8)`.

### Part C — Salinitystress case end-to-end

- [ ] **Task 9 — Audit `5.salinitystress` legacy files.**
  Read `tests/swap-cases/5.salinitystress/swap_linux.swp.template`, `swap.dra`, `potatod.crp`. Document case-specific quirks in `docs/phase-4c-b-salinitystress-audit.md`: which switches are set unusually (e.g., `swsalinity = 1` in the .crp), which sections present that other cases skip, anything that requires a new field or validator. No code change; this is a research task.
  **Verify:** audit doc exists and is referenced from later tasks.
  **Commit:** `docs(phase-4c-b): salinitystress legacy-file audit (Phase 4c-b Task 9)`.

- [ ] **Task 10 — Author `5.salinitystress/swap.toml` + `swap.dra.toml`.**
  Inside submodule. Mirror 4c-a authoring pattern (general/simulation/meteorology sections from `.swp`; drainage from `.dra` with cross-file reference). Bump outer-repo submodule pointer.
  **Verify:** `test_all_cases_smoke.pf` (4a) extended to include case 5 — passes.
  **Commit (submodule):** `feat(toml/salinitystress): swap.toml + swap.dra.toml (Phase 4c-b Task 10)`.
  **Commit (outer):** `chore(swap-cases): bump submodule for salinitystress swap+dra TOML`.

- [ ] **Task 11 — Author `5.salinitystress/potatod.crp.toml`.**
  Inside submodule. WOFOST TOML matching `tests/swap-cases/5.salinitystress/potatod.crp` (418-line file). `swsalinity = 1` activates Maas-Hoffman path. Differences from hupselbrook potatod tracked in the Task 9 audit doc. Bump submodule pointer.
  **Verify:** TOML loads via `load_swap_config`; validator passes.
  **Commit (submodule):** `feat(toml/salinitystress): potatod WOFOST TOML (Phase 4c-b Task 11)`.
  **Commit (outer):** `chore(swap-cases): bump submodule for salinitystress potatod TOML`.

- [ ] **Task 12 — Author `test_salinitystress_parity.pf`.**
  Pattern follows 4c-a `test_oxygenstress_parity.pf`: 7 subtests covering timing, meteorology, drainage, soil, crop scalars, crop tables, and a smoke "config validates clean". Crop side calls `read_legacy_wofost`. Iterate until green.
  **Verify:** all 7 sub-tests green; pFUnit count = previous + 7.
  **Commit:** `test(parity): salinitystress (case 5) full parity (Phase 4c-b Task 12)`.

### Part D — 4c-a tightening + closeout

- [ ] **Task 13 — Fill in `read_legacy_cropfixed` and `read_legacy_grass`.**
  Replace the Task 7 stubs with real implementations. Trace `cropgrowth.f90` for the exact arg lists at each rotation-entry call. Smoke-test each via the existing `test_legacy_crop_helper.pf`.
  **Verify:** smoke tests for each helper pass.
  **Commit:** `feat(test/io/toml): finish legacy_crop_helper fixed/grass wrappers (Phase 4c-b Task 13)`.

- [ ] **Task 14 — Tighten 4c-a parity tests.**
  Per file:
  - `test_hupselbrook_parity.pf` (maizes + grassd entries): replace hardcoded crop-scalar assertions with `read_legacy_cropfixed` / `read_legacy_grass` calls; assert `variables%foo == config%crop%rotation_fixed(N)%foo`.
  - `test_grassgrowth_parity.pf`: same, for grassd.
  - `test_oxygenstress_parity.pf`: same, for grassd.
  - `test_surfacewater_parity.pf`: same, for grass (type-1 fixed).
  Hardcoded values become test fixtures only when the legacy reader doesn't populate the field — document each with a `! READING NOTE: legacy <foo> not populated by readcropfixed; hardcoded from .crp` comment in the test.
  Each file gets its own commit.
  **Verify:** all 4 tests green; total pFUnit count unchanged or +N for new sub-tests.
  **Commit (split per file):** four commits, e.g.:
  - `test(parity): tighten hupselbrook crop parity to use legacy readers (Phase 4c-b Task 14a)`
  - `test(parity): tighten grassgrowth parity to use legacy readgrass (Phase 4c-b Task 14b)`
  - `test(parity): tighten oxygenstress parity to use legacy readgrass (Phase 4c-b Task 14c)`
  - `test(parity): tighten surfacewater parity to use legacy readcropfixed (Phase 4c-b Task 14d)`

- [ ] **Task 15 — Closeout: schema doc + coverage rebaseline + tag.**
  - Extend `docs/configuration-schema.md` with the type-2 section: full sub-section list, table encoding convention, validator summary.
  - Run `pixi run -e coverage coverage-report`. Update `docs/coverage-baseline.md` with 4c-b numbers.
  - Run `pixi run -e test check-full` — confirm 6/6 green. Save log under `tests/regression/baselines/phase-4c-b-wofost.log`.
  - Fast-forward `main` to `development` locally. Tag `rescue/phase-4c-b-wofost`. No push.
  **Verify:** tag exists; baselines updated; check-full + test-pfunit both green.
  **Commit:** `docs(baselines,schema): phase 4c-b closeout (Phase 4c-b Task 15)`.

---

## Risk register (operational)

| Risk | Mitigation |
|---|---|
| `readwofost` arg signature requires runtime values not yet read by `readswap()` | Task 7 traces the live call site before authoring the wrapper; if values are unavailable, hard-code defaults that match the live caller for a rotation entry |
| TOML table arrays-of-arrays don't parse cleanly in toml-f | If toml-f's array-of-array support is shaky, fall back to two parallel arrays per table (`dtsmtb_x = [...]; dtsmtb_y = [...]`) with a reader that zips them |
| WOFOST scalar parity fails because legacy reader applies unit conversions | Document each unit conversion in a `! READING NOTE` test comment; if conversions are non-trivial, add a `wofost_normalize` step in `read_cropwofost_toml` to match legacy behavior |
| Salinitystress case has fields that don't fit existing config sections | Task 9 audit catches early; either extend the section or document and defer to 4d |
| Tightening surfaces silent 4c-a parity bugs | Each surfaced bug becomes a focused commit with the legacy-reader assertion as the regression test |
| `check-full` red after Part A | No physics path changes in Part A — if red, a meson wiring or `use` chain issue; bisect by toggling new sources off |

---

## Verification gates

After Task 5 (Part A done): `check-full` 6/6 green, `test-pfunit` green with new tests added.
After Task 8 (Part B done): hupselbrook full parity green including potatod WOFOST.
After Task 12 (Part C done): all 5 non-macropore cases parity-green.
After Task 14 (Part D tightening done): no hardcoded crop assertions remain in any 4c-a parity test (search: `! READING NOTE.*hardcoded`).
After Task 15 (closeout): `rescue/phase-4c-b-wofost` tagged; coverage rebaselined.
