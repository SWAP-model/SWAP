# Phase 4d — Bottom Boundary, Heat, Irrigation, Solute Configs + Grass Mowing/Grazing Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Land four new config types (`bottom_boundary`, `heat`, `irrigation`, `solute`) plus extend `cropgrass_config_t` with per-event mowing/grazing tables. Roll out `[bottom_boundary]`, `[heat]`, `[irrigation]`, `[solute]` sections to all five non-macropore cases. Extend each parity test with assertions for the new fields.

**Architecture:** Same composite reader pattern from 4a/b/c. New types live in `src/config/`; new readers in `src/io/toml/`. `irrigation_config_t` (top-level, fixed-irrigation) + `irrigation_schedule_t` (nested per-crop sub-type) — Option A per user decision. Each per-crop config (`cropfixed`, `cropgrass`, `cropwofost`) gains a `schedule :: irrigation_schedule_t` field.

**Tech Stack:** Unchanged from 4c — gfortran 2008, pFUnit 4.15, meson + pixi, toml-f, `iso_c_binding` chdir.

Spec: `docs/superpowers/specs/2026-04-27-phase-4d-remaining-configs-design.md`

---

## Preamble: context every task needs

**Baseline:** Phase 4c-b complete at commit `2127e26`, tag `rescue/phase-4c-b-wofost`. 221 pFUnit tests passing, 6/6 regression green. WOFOST type-2 fully wired with parity for hupselbrook potatod and salinitystress. The 4c-a parity tests use `legacy_crop_helper` wrappers (`read_legacy_cropfixed`, `read_legacy_grass`, `read_legacy_wofost`).

**Branch discipline:**
- Work on `development`. No per-task feature branches.
- One commit per task. Subject `<type>(<scope>): <what>`.
- No pushes to origin during 4d (matches rescue policy).
- Phase exit: fast-forward `main` to `development` locally; tag `rescue/phase-4d-remaining-configs`.

**Working directory:** `/home/zawadzkim/Code/swap` for outer-repo work, `tests/swap-cases/` for submodule TOML authoring.

**Conventions** (unchanged):
- Module name = filename stem + `_mod`.
- `implicit none` after module statement; default `private`; explicit `public ::`.
- `use iso_fortran_env, only: real64` for new code.
- pFUnit suite filename = `tests/unit/<domain>/test_<module>.pf`; register in `pf_files` and `testSuites.inc`.
- **No inline `! comments` on `@assert*` lines** — `funitproc` chokes on the `!`.
- Submodule edits use the file-write-only pattern from 4c-b — subagent writes files, user batch-commits at end.

**Verification after each task:**
```
pixi run -e test test-pfunit         # F count must stay 0
pixi run -e test check-fast          # 4/6 cases; <90s
```
After Part F (grass extension) and again at closeout, also:
```
pixi run -e test check-full          # 6/6 cases; ~6 min
```

**Three concepts worth understanding before authoring tasks:**

1. **`SWBOTB` legacy storage.** The 8 SWBOTB branches use mostly-disjoint sets of legacy globals. SWBOTB=3 (Cauchy with regional aquifer head) populates `shape, hdrain, rimlay, aqave, aqamp, aqomeg`. Other branches populate fewer. The new `bottom_boundary_config_t` declares all fields as nullable (allocatable for tables, scalar with sensible defaults for the rest); the validator gates field requirements on `swbotb` value.

2. **Solute vs salinity disjointness.** Crop-side salinity stress (`config%crop%rotation_wofost(N)%salinity`) and solute-side parameters (`config%solute%cdrain` etc.) are independent. Validators do NOT cross-check. `salinity.swsalinity=1` in the crop config asks the crop water-uptake routine to apply Maas-Hoffman. `solute.swsolu=1` in the solute config asks the solute transport routine to actually transport salt. They can be on/off independently.

3. **Irrigation Option A nesting.** `swap_config_t%irrigation` holds top-level fixed-irrigation. Each per-crop type (`cropfixed_config_t%schedule`, `cropwofost_config_t%schedule`, `cropgrass_config_t%schedule`) holds the per-crop scheduling block. Both are populated independently from the `.swp` and `.crp` files respectively. The validator on `swap_config_t` calls `config%irrigation%validate()`; the validator on each crop type calls `self%schedule%validate()`. No cross-type coupling.

---

## File structure

### New source files

| File | Responsibility |
|---|---|
| `src/config/bottom_boundary_config.f90` | `bottom_boundary_config_t` — all SWBOTB branches. |
| `src/config/heat_config.f90` | `heat_config_t` — SWHEA + frost + initial-temperature table. |
| `src/config/irrigation_config.f90` | `irrigation_config_t` (top-level) + `irrigation_schedule_t` (nested). Both in same file (related). |
| `src/config/solute_config.f90` | `solute_config_t` — solute transport, dispersion, salinity. |
| `src/io/toml/read_bottom_boundary_toml.f90` | Section reader for `[bottom_boundary]`. |
| `src/io/toml/read_heat_toml.f90` | Section reader for `[heat]`. |
| `src/io/toml/read_irrigation_toml.f90` | Section reader for top-level `[irrigation]`. |
| `src/io/toml/read_solute_toml.f90` | Section reader for `[solute]`. |
| `tests/unit/config/test_bottom_boundary_config.pf` | Validator coverage. |
| `tests/unit/config/test_heat_config.pf` | Validator coverage. |
| `tests/unit/config/test_irrigation_config.pf` | Validator coverage (top-level + nested schedule). |
| `tests/unit/config/test_solute_config.pf` | Validator coverage. |
| `tests/unit/io/toml/test_read_bottom_boundary_toml.pf` | Section-reader tests. |
| `tests/unit/io/toml/test_read_heat_toml.pf` | Section-reader tests. |
| `tests/unit/io/toml/test_read_irrigation_toml.pf` | Section-reader tests. |
| `tests/unit/io/toml/test_read_solute_toml.pf` | Section-reader tests. |
| `docs/phase-4d-case-audit.md` | Case audit (Task 1 output). |

### Modified source files

| File | Change |
|---|---|
| `src/config/swap_config.f90` | Add `bottom_boundary`, `heat`, `irrigation`, `solute` fields. Extend validate/finalize dispatches. |
| `src/config/cropgrass_config.f90` | Extend with mowing/grazing per-event tables; add `schedule :: irrigation_schedule_t`. |
| `src/config/cropwofost_config.f90` | Add `schedule :: irrigation_schedule_t`. |
| `src/config/cropfixed_config.f90` | Add `schedule :: irrigation_schedule_t`. |
| `src/io/toml/read_cropgrass_toml.f90` | Extend for new tables + `[irrigation_schedule]` block. |
| `src/io/toml/read_cropwofost_toml.f90` | Extend for `[irrigation_schedule]` block. |
| `src/io/toml/read_cropfixed_toml.f90` | Extend for `[irrigation_schedule]` block. |
| `src/io/toml/load_swap_config.f90` | Call new top-level readers. |
| `tests/unit/io/toml/test_hupselbrook_parity.pf` | Add bottom_boundary/heat/irrigation/solute assertions. |
| `tests/unit/io/toml/test_grassgrowth_parity.pf` | Same + grass mowing/grazing assertions. |
| `tests/unit/io/toml/test_oxygenstress_parity.pf` | Same + grass mowing/grazing. |
| `tests/unit/io/toml/test_surfacewater_parity.pf` | Add new sections. |
| `tests/unit/io/toml/test_salinitystress_parity.pf` | Add new sections (case 5 has SWBOTB=3, SWSOLU=1). |
| `tests/unit/meson.build`, `tests/unit/testSuites.inc`, `meson.build` | Register new sources/tests. |
| `docs/configuration-schema.md` | Add four new sections + grass mowing extension. |

### Submodule edits (`tests/swap-cases/toml/`)

For each of cases 1, 2, 4, 5, 6:
- `swap.toml` — extend with `[bottom_boundary]`, `[heat]`, `[irrigation]`, `[solute]` sections.
- (Cases 2, 4): `grassd.crp.toml` — extend with mowing/grazing tables.
- (Cases that use scheduled irrigation per Task 1 audit): per-crop TOMLs — add `[irrigation_schedule]` blocks.

---

## Tasks

### Part A — Audit + foundation

- [ ] **Task 1 — Case audit.**
  Read each of the 5 non-macropore cases' `.swp`, `.dra`, `.crp` files. Record per-case values for `swbotb`, `swhea`, `swsolu`, `swirfix`, `schedule` (per crop), and any other 4d-relevant switch. Document in `docs/phase-4d-case-audit.md` with one table per switch.
  **Verify:** doc exists; covers all 5 cases.
  **Commit:** `docs(phase-4d): case audit (Phase 4d Task 1)`.

### Part B — Bottom boundary

- [ ] **Task 2 — `bottom_boundary_config_t` skeleton.**
  Author `src/config/bottom_boundary_config.f90`. Top-level type with all SWBOTB-1..8 fields declared (most as scalars with sensible defaults, tables as `real(real64), allocatable :: foo(:,:)`). Stub `validate` and `finalize`. Register in `meson.build` `sources` and `tests/unit/meson.build` `pfunit_extra_sources`.
  **Verify:** `test-pfunit` green.
  **Commit:** `feat(config): add bottom_boundary_config_t skeleton (Phase 4d Task 2)`.

- [ ] **Task 3 — `bottom_boundary_config_t` validators.**
  Add `validate` with branches per `swbotb` value. Pass-through `validate` for `swbotb=4` (free drainage — no params). Add `tests/unit/config/test_bottom_boundary_config.pf` covering each branch's happy path + at least one failure mode per branch (~16 tests).
  **Verify:** new tests pass.
  **Commit:** `feat(config): add bottom_boundary_config_t validators (Phase 4d Task 3)`.

- [ ] **Task 4 — `read_bottom_boundary_toml`.**
  Author `src/io/toml/read_bottom_boundary_toml.f90`. Reads `[bottom_boundary]` section + optional sub-sections per branch. Plus `tests/unit/io/toml/test_read_bottom_boundary_toml.pf` (~6 tests: happy paths for SWBOTB=2, =3, =4 + 3 error paths).
  **Verify:** new tests pass.
  **Commit:** `feat(io/toml): add read_bottom_boundary_toml (Phase 4d Task 4)`.

- [ ] **Task 5 — Wire into `swap_config_t` + `load_swap_config`.**
  Add `type(bottom_boundary_config_t) :: bottom_boundary` to `swap_config_t`. Extend validate/finalize dispatch. In `load_swap_config.f90`, call `read_bottom_boundary_toml` after the soil reader. Add a smoke test in `test_load_swap_config.pf` exercising the dispatch.
  **Verify:** `test-pfunit` green; `check-fast` 4/4 green.
  **Commit:** `feat(config,io/toml): wire bottom_boundary into swap_config (Phase 4d Task 5)`.

### Part C — Heat

- [ ] **Task 6 — `heat_config_t` skeleton + validators.**
  Single commit: skeleton + validators + tests. ~12 fields, ~6 validator tests. Register in meson + testSuites.
  **Commit:** `feat(config): add heat_config_t (Phase 4d Task 6)`.

- [ ] **Task 7 — `read_heat_toml` + tests.**
  ~5 reader tests (SWHEA=0 silent, SWHEA=1 happy, frost params, soil-temperature table).
  **Commit:** `feat(io/toml): add read_heat_toml (Phase 4d Task 7)`.

- [ ] **Task 8 — Wire into `swap_config_t`.**
  Same pattern as Task 5.
  **Commit:** `feat(config,io/toml): wire heat into swap_config (Phase 4d Task 8)`.

### Part D — Irrigation

- [ ] **Task 9 — `irrigation_config_t` + `irrigation_schedule_t`.**
  Both types in `src/config/irrigation_config.f90`. Top-level type ~10 fields; schedule sub-type ~25 fields with table allocatables for TCS=1..8 branches. Stub validators.
  **Commit:** `feat(config): add irrigation_config_t and irrigation_schedule_t skeletons (Phase 4d Task 9)`.

- [ ] **Task 10 — Validators + tests.**
  Both validators with branch coverage. ~12 validator tests.
  **Commit:** `feat(config): add irrigation validators (Phase 4d Task 10)`.

- [ ] **Task 11 — Top-level `read_irrigation_toml` + per-crop reader extensions.**
  Author `src/io/toml/read_irrigation_toml.f90`. Extend `read_cropwofost_toml.f90`, `read_cropfixed_toml.f90`, `read_cropgrass_toml.f90` to read `[irrigation_schedule]` block when present. Tests for both (top-level and per-crop).
  **Commit:** `feat(io/toml): add read_irrigation_toml and per-crop schedule readers (Phase 4d Task 11)`.

- [ ] **Task 12 — Wire into `swap_config_t` + per-crop types.**
  Add `irrigation` field at top level; add `schedule :: irrigation_schedule_t` to each per-crop type. Update validate dispatches. Smoke test.
  **Commit:** `feat(config,io/toml): wire irrigation top-level + per-crop schedule (Phase 4d Task 12)`.

### Part E — Solute

- [ ] **Task 13 — `solute_config_t` + validators + reader + tests.**
  Single commit: skeleton + validators + reader. ~25 fields. ~10 tests across config + reader.
  **Commit:** `feat(config,io/toml): add solute_config_t and reader (Phase 4d Task 13)`.

- [ ] **Task 14 — Wire into `swap_config_t`.**
  **Commit:** `feat(config,io/toml): wire solute into swap_config (Phase 4d Task 14)`.

### Part F — Grass mowing/grazing extension

- [ ] **Task 15 — Extend `cropgrass_config_t` with per-event tables.**
  Add `mowing_dates(:)`, `mowing_heights(:)`, `nmow`, `swdmmow`, `dmharvest`, `daylastharvest`, `dmlastharvest`, `maxdaymow`, `nstart_graz`, `nstop_graz`, `maxdaygrz`, `dmgrazing`, `swdmgrz`, `lsdb_default`. Validator branches for `swharv=0/1/2`. Update existing `test_cropgrass_config.pf` with new branch coverage.
  **Commit:** `feat(config): extend cropgrass_config_t with per-event mowing/grazing tables (Phase 4d Task 15)`.

- [ ] **Task 16 — Extend `read_cropgrass_toml.f90` for new tables.**
  Read mowing/grazing tables conditionally. Update `test_read_cropgrass_toml.pf`.
  **Commit:** `feat(io/toml): extend read_cropgrass_toml for mowing/grazing tables (Phase 4d Task 16)`.

### Part G — Per-case TOML rollout

- [ ] **Task 17 — Extend per-case `swap.toml` with new sections (5 cases).**
  For each of cases 1, 2, 4, 5, 6, edit the case's `swap.toml` to add `[bottom_boundary]`, `[heat]`, `[irrigation]`, `[solute]` sections. Values from Task 1 audit. Submodule edits — DO NOT commit submodule yourself; user batch-commits at end. Update outer-repo `test_all_cases_smoke.pf` if needed (should still pass).
  **Verify:** smoke tests green for all 5 cases; F count 0.
  **Commit (outer-repo only — no submodule changes yet, just code that exercises them later):** none in this task; defer to user-driven submodule commit + outer-repo bump after Task 19.

- [ ] **Task 18 — Extend cases 2 and 4 `grassd.crp.toml` with mowing/grazing tables.**
  Submodule edit. Same submodule-commit pattern as Task 17.

- [ ] **Task 19 — Extend per-crop TOMLs with `[irrigation_schedule]` blocks (case-specific).**
  Per Task 1 audit, identify which crops have scheduled irrigation. Add `[irrigation_schedule]` blocks to those `.crp.toml` files. Submodule edit.

- [ ] **User commits Tasks 17-19 submodule changes + outer-repo bump.**
  Single user-driven step. Subagent reports the file list; user runs:
  ```
  cd tests/swap-cases
  git add toml/<paths>
  git commit -m "feat(toml): per-case extensions for Phase 4d sections (Phase 4d Tasks 17-19)"
  cd ..
  git add tests/swap-cases
  git commit -m "chore(swap-cases): bump submodule for Phase 4d per-case extensions"
  ```

### Part H — Parity-test extensions

- [ ] **Task 20 — Extend each parity test with new-section assertions.**
  Walk all 5 parity tests (`test_hupselbrook_parity.pf`, `test_grassgrowth_parity.pf`, `test_oxygenstress_parity.pf`, `test_surfacewater_parity.pf`, `test_salinitystress_parity.pf`). For each, add sub-tests for `bottom_boundary`, `heat`, `irrigation`, `solute`. Use existing `load_both_for_<case>` helpers (extended if needed). For grass cases, add mowing/grazing assertions; for crop cases, add per-crop irrigation_schedule assertions if the case uses scheduled irrigation. ~30-50 new assertions per case.
  **Verify:** all parity tests green; F count 0; check-full 6/6 green.
  **Commit (one per file, 5 total):**
  - `test(parity): add 4d-section assertions to hupselbrook (Phase 4d Task 20a)`
  - `test(parity): add 4d-section assertions + grass extensions to grassgrowth (Phase 4d Task 20b)`
  - `test(parity): add 4d-section assertions + grass extensions to oxygenstress (Phase 4d Task 20c)`
  - `test(parity): add 4d-section assertions to surfacewater (Phase 4d Task 20d)`
  - `test(parity): add 4d-section assertions to salinitystress (Phase 4d Task 20e)`

### Part I — Closeout

- [ ] **Task 21 — Schema doc extension.**
  Extend `docs/configuration-schema.md` with the four new sections + grass mowing/grazing extension.
  **Commit:** `docs(schema): document Phase 4d sections (Phase 4d Task 21)`.

- [ ] **Task 22 — Coverage rebaseline.**
  Run `pixi run -e coverage coverage-report`. Update `docs/coverage-baseline.md` with 4d numbers.
  **Commit:** `docs(baselines): rebaseline coverage at Phase 4d (Phase 4d Task 22)`.

- [ ] **Task 23 — check-full final log + tag + main fast-forward.**
  Run `pixi run -e test check-full`. Save log to `tests/regression/baselines/phase-4d-remaining-configs.log`. Fast-forward `main` to `development` locally. Tag `rescue/phase-4d-remaining-configs`. No push.
  **Commit:** `docs(baselines): record Phase 4d check-full output (Phase 4d Task 23)`.

---

## Risk register (operational)

| Risk | Mitigation |
|---|---|
| SWBOTB=3 case 5 has unit-conversion subtleties (e.g., aquifer head in cm vs m) | Task 1 audit records exact values; Task 4 reader passes them through; Task 20e parity catches any unit drift |
| Irrigation `irgfil` is a path reference — case may have its own `swap.irg` file | Task 11 reader reads the path; not Phase 4d's job to follow it (deferred to a future phase if needed). Validator just checks the file path is non-empty when `swirfix=1` |
| Grass mowing tables in case 2 currently have `swharv=0` workaround that needs removal | Task 15 explicitly removes the workaround when extending; commit message says so |
| Parity test expansion (Task 20) surfaces silent legacy-reader bugs | Each surfaced bug becomes a focused commit with a regression test |
| Submodule edits across 3 tasks (17-19) bunch up; user has to commit them all at once | Subagent reports a clean file list + commit script for user; mirrors 4c-b's pattern |
| Coverage stays flat at 53.1% (gcovr issue) | Pre-existing, document only |
| `cropgrass_config_t` extension may break 4c-a/b validator tests | Task 15 explicitly updates `test_cropgrass_config.pf` validators |

---

## Verification gates

After Task 5: bottom_boundary wired, smoke test green.
After Task 8: heat wired.
After Task 12: irrigation wired (both halves).
After Task 14: solute wired.
After Task 16: grass per-event tables wired.
After Tasks 17-19 + user commit: all 5 cases load `swap_config_t` cleanly with the new sections.
After Task 20: all 5 parity tests pass with new-section assertions.
After Task 23: tag `rescue/phase-4d-remaining-configs` exists; main fast-forwarded.
