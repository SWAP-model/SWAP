# Phase 4f-prep — Schema Gap Closure (Top-N) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close the top-N priority real-domain gaps from the Phase 4e audit so Phase 4f's strangler-fig of `readswap()` opens with ≤30 remaining G entries. Add `macropore_config_t` (full 22-key) and `surface_water_config_t` (new); extend five existing configs; per-case TOMLs and parity tests only where the new section is exercised.

**Architecture:** Same patterns as Phase 4d. Each new config: skeleton → validators (with switch-gating per Phase 4d's `bottom_boundary` precedent) → reader → wire into `swap_config_t` + `load_swap_config`. Each existing-config extension: extend type → extend validator → extend reader → extend tests. Per-case TOML: only cases whose legacy `.swp` exercises the new section.

**Tech Stack:** Unchanged — gfortran 2008, pFUnit 4.15, meson + pixi, toml-f.

Spec: `docs/superpowers/specs/2026-04-28-phase-4f-prep-gap-closure-design.md`

---

## Preamble: context every task needs

**Baseline:** Phase 4e complete at commit `c4fe7eb`, tag `rescue/phase-4e-error-prep`. main and development aligned. 317 pFUnit tests, 6/6 regression green. Submodule at `952341e` (clean).

**Branch discipline:**
- Work on `development`. No feature branches.
- One commit per task. Subject `<type>(<scope>): <what> (Phase 4f-prep Task N)`.
- No pushes during 4f-prep.
- Phase exit: fast-forward main; tag `rescue/phase-4f-prep-gap-closure`.

**Working directory:** `/home/zawadzkim/Code/swap`. Submodule edits in `tests/swap-cases/`.

**Conventions** (unchanged from 4d/4e):
- Module name = filename stem + `_mod`.
- `implicit none` after module statement; default `private`; explicit `public ::`.
- `use iso_fortran_env, only: real64` for new code.
- pFUnit: `'<domain>/test_<module>.pf'` in `pf_files`; `ADD_TEST_SUITE(test_<module>_suite)` in `testSuites.inc`; module sources in `pfunit_extra_sources` if outside `test_base_sources`.

**Verification after each task:**
```
pixi run -e test test-pfunit
./builddir/tests/unit/unit-swap-tests > /tmp/v.log 2>&1
echo "exit=$?"
grep -E "^[\.F ]+$|^\.+$" /tmp/v.log | tr -dc 'F' | wc -c
pixi run -e test check-fast
```
After each Part of meaningful size:
```
pixi run -e test check-full
```

**Two concepts worth understanding before authoring tasks:**

1. **Switch-gated validators.** Most new fields are required only when their gating switch is 1. Pattern from Phase 4d's `bottom_boundary_config_t`: validator body wraps each branch in `if (self%switch == N)`. Cases that don't exercise the section have `switch = 0` (default) and their TOML has no section block — the reader leaves the type at defaults; the validator skips its body.

2. **Per-case parity-test extension only when exercised.** Adding a `[macropore]` schema does NOT require updating all 6 parity tests — only `test_macroporeflow_parity.pf` (case 3 is the only case with `swmacro=1`). Other cases stay green because their TOML has no `[macropore]` block, the reader doesn't read one, and the validator skips. The parity assertion catches mismatches only where the section is actually used.

---

## File structure

(See spec section "Critical files" for the full file list. Plan's task ordering below references the same files.)

---

## Tasks

### Part A — Audit verification + macropore_config_t

- [ ] **Task A0 — Verify `croptype` audit flag.**
  Research-only. Run `grep -rn "\bcroptype\b" src/ docs/` and inspect each match. The audit listed `croptype` as a top-10 G entry, but `crop_config_t` already has `rotation_type(:)`. Determine: is `croptype` a distinct legacy field, or just an alias the auditor missed? Document in 1-2 paragraphs at the top of `docs/phase-4f-config-to-variables-audit.md` (extending the existing audit doc, not a new file). If false-positive, mark croptype as C with note. If real, update Task scope and proceed.
  **Commit:** `docs(phase-4f-prep): verify croptype audit-flag (Phase 4f-prep Task A0)`.

- [ ] **Task A1 — `macropore_config_t` skeleton + validators + tests.**
  Create `src/config/macropore_config.f90` with the 13 scalars + 9 per-layer table fields per spec D1. Write validators with switch-gating on `swmacro`. Per-layer table sizes consistent with `nlayers`.
  Author `tests/unit/config/test_macropore_config.pf` covering: (a) `swmacro=0` validates clean with everything default; (b) `swmacro=1` requires scalar bounds; (c) array-size mismatch surfaces; (d) full happy-path with all 22 fields populated.
  Register in `meson.build` `sources` and `tests/unit/meson.build`.
  **Verify:** test-pfunit green; new tests pass; existing 317 tests still pass.
  **Commit:** `feat(config): add macropore_config_t (Phase 4f-prep Task A1)`.

- [x] **Task A2 — SKIPPED per ADR 0010.** No `read_macropore_toml`; the new TOML pipeline does not parse `[macropore]`. `macropore_config_t` stays as orphan infrastructure for future macropore work.

- [x] **Task A3 — SKIPPED per ADR 0010.** No wiring of `macropore_config_t` into `swap_config_t`. Case 3 (3.macroporeflow) continues to run on the legacy `readswap()` path; Phase 4f's strangler-fig keeps `readswap.f90` as a fallback for macropore cases.

### Part B — surface_water_config_t

- [ ] **Task B1 — `surface_water_config_t` skeleton + validators + tests.**
  Same shape as A1, smaller. Create `src/config/surface_water_config.f90` with `swsrf`, `swsec`, `nmper`, `dropr`, `wlsman`, `gwlcrit`, `nphase`, `wscap`, `wldip`, plus the `RapDra*` table.
  Switch-gating on `swsrf`/`swsec`.
  Author `tests/unit/config/test_surface_water_config.pf` (~6 tests).
  **Commit:** `feat(config): add surface_water_config_t (Phase 4f-prep Task B1)`.

- [ ] **Task B2 — `read_surface_water_toml` + tests.**
  `src/io/toml/read_surface_water_toml.f90`. Reads `[surface_water]` section.
  `tests/unit/io/toml/test_read_surface_water_toml.pf`.
  **Commit:** `feat(io/toml): add read_surface_water_toml (Phase 4f-prep Task B2)`.

- [ ] **Task B3 — Wire `surface_water` into `swap_config_t` + `load_swap_config`.**
  Add to `swap_config_t`. Wire reader.
  **Commit:** `feat(config,io/toml): wire surface_water into swap_config (Phase 4f-prep Task B3)`.

- [ ] **Task B4 — Part A+B check-full baseline.**
  Run `check-full`; record at `tests/regression/baselines/phase-4f-prep-task-b4-check-full.log`. Expect 6/6 green (purely additive; no physics changes).
  **Commit:** `docs(baselines): record phase 4f-prep Part A+B check-full (Phase 4f-prep Task B4)`.

### Part C — Extensions to existing configs

- [ ] **Task C1 — `simulation_config_t.numerical` sub-section.**
  Add `[simulation.numerical]` fields: `dt`, `dtmin`, `dtmax`, `MaxIt`, `MaxBackTr`, `taccur`. Extend validator with bound checks. Extend `read_simulation_toml`. Extend `tests/unit/config/test_simulation_config.pf` and `tests/unit/io/toml/test_read_simulation_toml.pf`.
  **Commit:** `feat(config,io/toml): add simulation.numerical sub-section (Phase 4f-prep Task C1)`.

- [ ] **Task C2 — `meteorology_config_t` snow + rainfall extension.**
  Add `[meteorology.snow]` (`swsnow` + sub-fields) and extend `[meteorology.rainfall]` (`cnref`, `cndry`, `cnwet`, `cofred`, `cfeic`, `cfevappond`).
  Extend validator + reader + per-section tests.
  **Commit:** `feat(config,io/toml): add meteorology snow + extend rainfall (Phase 4f-prep Task C2)`.

- [ ] **Task C3 — `soil_config_t` extensions.**
  Add `[soil.discretization]` (`swdiscrvert`, `dznew(:)`, `numnodnew`); `[soil.anisotropy]` (`cofani(:)`); `[soil.frost]` (`swfrost`, `tfroststa`, `tfrostend`, `swsublim`); `nrstaring` top-level.
  Extend validator + reader + tests.
  **Commit:** `feat(config,io/toml): extend soil with discretization+anisotropy+frost (Phase 4f-prep Task C3)`.

- [ ] **Task C4 — `drainage_config_t.surface_runoff` sub-section.**
  Add `[drainage.surface_runoff]` (`swnrsrf`, `SwTopnrsrf`, `swdivdinf`, `FacDpthInf`).
  Extend validator + reader + tests.
  **Commit:** `feat(config,io/toml): add drainage.surface_runoff sub-section (Phase 4f-prep Task C4)`.

- [ ] **Task C5 — `crop_config_t.rotation` extensions.**
  Add `swhydrlift` per-rotation-entry. If A0 confirmed `croptype` is a real distinct field, add it here too.
  Extend validator + reader + tests.
  **Commit:** `feat(config,io/toml): extend crop.rotation with swhydrlift (Phase 4f-prep Task C5)`.

- [ ] **Task C6 — Part C check-full.**
  Record at `tests/regression/baselines/phase-4f-prep-task-c6-check-full.log`.
  **Commit:** `docs(baselines): record phase 4f-prep Part C check-full (Phase 4f-prep Task C6)`.

### Part D — Per-case TOML extensions + parity coverage

Per Q3 option A: only cases that exercise the new section get TOML extensions and parity-test assertions.

- [x] **Task D1 — SKIPPED per ADR 0010.** No `[macropore]` block in case 3's TOML; no parity assertions for macropore fields. Case 3's existing parity test (Phase 4e Task C3, commit `37e79c6`) already scopes to the schema-covered subset and stays green.

- [ ] **Task D2 — Surface-water TOML in case 6 + parity assertions.**
  Same pattern, case 6.

- [ ] **Task D3 — Snow / frost TOML extensions where exercised.**
  Audit which cases have `swsnow=1`, `swfrost=1`, `swsublim=1`. Extend their TOMLs with the new sub-sections + parity assertions in their `test_*_parity.pf`.

- [ ] **Task D4 — Numerical sub-section TOML extensions where exercised.**
  Audit which cases have non-default `dt`/`dtmax`/etc. Extend.

- [ ] **Task D5 — Soil discretization / anisotropy where exercised.**
  Same pattern.

- [ ] **Task D6 — Drainage surface_runoff where exercised.**
  Same pattern.

- [ ] **User commits Tasks D1-D6 submodule changes + outer-repo bumps.**

- [ ] **Task D7 — Part D check-full.**
  Record at `tests/regression/baselines/phase-4f-prep-task-d7-check-full.log`.

### Part E — Closeout

- [ ] **Task E1 — Update audit doc.**
  In `docs/phase-4f-config-to-variables-audit.md`: reclassify all the closed gaps from G to C with their new schema paths. Update summary counts. Final real-domain G count target: ≤30.
  **Commit:** `docs(phase-4f): reclassify closed gaps after Phase 4f-prep (Phase 4f-prep Task E1)`.

- [ ] **Task E2 — Schema doc extension.**
  Document the new `[macropore]`, `[surface_water]`, `[simulation.numerical]`, `[meteorology.snow]`, `[soil.discretization|anisotropy|frost]`, `[drainage.surface_runoff]` sections in `docs/configuration-schema.md`.
  **Commit:** `docs(schema): document Phase 4f-prep extensions (Phase 4f-prep Task E2)`.

- [ ] **Task E3 — Coverage rebaseline.**
  Run `pixi run -e coverage coverage-report`. Update `docs/coverage-baseline.md` with Phase 4f-prep section.
  **Commit:** `docs(coverage): rebaseline at Phase 4f-prep (Phase 4f-prep Task E3)`.

- [ ] **Task E4 — Closeout: check-full + tag + main fast-forward.**
  Final check-full. Save log. Tag `rescue/phase-4f-prep-gap-closure`. Fast-forward main.

---

## Risks & mitigations

| Risk | Mitigation |
|---|---|
| `croptype` turns out to be real distinct field, growing C5 scope | A0 verifies first; if real, scope grows by 1-2 fields and 1 commit |
| Per-layer table sizing in `macropore_config_t` requires `nlayers` from `soil_config_t` not yet finalized at validate time | Schema-cross-check happens at finalize; the macropore validator can defer per-layer size checks to its own finalize stage if needed |
| `surface_water_config_t` overlaps with existing `drainage_config_t` legacy semantics | The legacy `.swp` block names distinguish them; if conflict surfaces, scope down B-section to only the cleanly-distinct fields and document |
| Per-case TOML extension iterates and surfaces hidden TOML/reader bugs | That's the point — parity tests catch them. Each fix is a separate focused commit |
| 30-gap target is overshot (more gaps than expected stay open) | Acceptable; document at closeout. Phase 4f opens with whatever G count we reach; gaps surfaced at runtime become focused fix commits inside Phase 4f |

---

## Verification gates

After each Part: `test-pfunit` green; new test count grows by the expected number; F count 0; check-fast 4/4 green.
After Parts B/C/D: `check-full` 6/6 green.
Final: real-domain G count in audit doc ≤30; tag exists; main fast-forwarded.

## File count summary

- Create: ~10 new files (2 config types + 2 readers + 4 unit-test files + 2 misc docs/baselines).
- Modify: ~20 files (existing configs + readers + parity tests + audit doc + schema doc).
- LoC delta: ~+2500 / -100.
