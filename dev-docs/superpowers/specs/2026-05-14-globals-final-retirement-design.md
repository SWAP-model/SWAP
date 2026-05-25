---
title: "GR-FINAL — Globals Retirement Arc 9 (Final retirement; delete the 3 adapter files)"
date: 2026-05-14
status: approved
context: globals-retirement Arc 9 — capstone; delete variables.f90 + config_to_variables.f90 + initialize.f90
---

# GR-FINAL Design

## Goal

Migrate all 32 remaining `use variables` consumers, retire the swap_mod strangler-fig leftovers (transient buffers + swinco=3 inline block), and delete the 3 adapter files:

- `src/core/variables.f90` (~1414 lines)
- `src/io/toml/config_to_variables.f90` (~1775 lines)
- `src/core/initialize.f90` (~851 lines)

Total deletion target: **~4040 lines of legacy plumbing**.

End state: all simulation data flows through typed `state%X` (runtime mutable) and `config%X` (read-only after init). No bare globals. The strangler-fig fully demolished.

## Architecture

Five phases:

- **Phase A** — Preparation: schema completeness audit, retirement-readiness inventory, baseline.
- **Phase B** — Major cluster migrations: readmeteo, swapoutput (split into 2), cropgrowth + crop runtime, timecontrol_mod, soil/drainage, small consumers.
- **Phase C** — Strangler-fig elimination: swap_mod transient buffers + swinco=3 inline block + config_to_variables.f90 incremental shrinkage.
- **Phase D** — Global retirement: iterative variables.f90 declaration deletion + initialize.f90 zero-fill cleanup. Compile errors surface any missed reader → iterate.
- **Phase E** — Final adapter deletion: delete 3 adapter files; meson cleanup; final ADR.

## Verification Discipline

Per `feedback_per_task_regression_gate.md` + `feedback_state_schema_clean_rebuild.md`:

- **Every implementer subagent ends with VG:**
  ```bash
  rm -rf builddir && pixi run build-linux \
    && pixi run test-pfunit \
    && pixi run -e test python tests/regression/test_output_regression.py \
         hupselbrook surfacewater salinitystress grassgrowth
  ```
  Expected: clean build, pFUnit 741, regression **4/4 byte-for-byte**.

- **Every phase ends with:** `pixi run check-full` 5/5 byte-for-byte.

- **Controller quality bar (per user directive):** subagents do NOT defer due to missing infrastructure. Controller intervenes inline with state/config arg threading, schema extensions, or sibling migrations.

## Phase A — Preparation

### A1: Pre-flight baseline

Lock baseline. Run check-full 5/5; empty marker commit.

### A2: Audit retirement readiness

Categorize all ~150 remaining declarations in `src/core/variables.f90` into:

1. **Write-only-by-adapter** (no remaining readers) — directly deletable after adapter writes drop.
2. **Read-only-by-readers** (consumers still read; need migration).
3. **Both** (need migration + adapter-rewrite).
4. **Config-loaded mirror** (state mirror is sourced via adapter; convert to direct config-sourcing then delete legacy global).

Produce inventory note: `docs/superpowers/notes/2026-05-14-gr-final-retirement-inventory.md` with per-symbol classification.

### A3: Add `max_resp_factor` to config

Verify if `max_resp_factor` is in any config_t. If not, add to `config%crop%fixed` (or wherever consistent with legacy). Update adapter populator.

### A4: Phase A close marker

check-full 5/5; commit marker.

## Phase B — Major cluster migrations

### B1: readmeteo bare blanket (line 17)

`src/io/readmeteo.f90` has `use variables` at line 17 with implicit access to many globals. Convert to narrow imports / state refs.

The bare site in `readmeteo` (top of file) governs MANY subroutines. The replacement must thread `state` through each subroutine. Many already take state from prior arc work — verify.

Expected: empty `use variables` line after migration (or one narrow deferral with comment).

### B2: readmeteo CSV writers (lines 267, 366, 461, 554)

Four narrow sites. Each writes per-year meteo data:
- Line 267: `raincsv_dat, nraincsv` — rain CSV input data
- Lines 366, 461: `arad, atmn, atmx, ahum, awin, arai, aetr, wet, ad, am, ...` — daily meteo arrays (already on state%atmosphere)
- Line 554: `metcsv_det, nmetcsv_det` — meteo detail (sub-daily) CSV

Migrate each writer to populate `state%atmosphere%X` directly. After this, the legacy `arad/atmn/atmx/...` globals become unused. Retirement happens in Phase D1.

### B3: swapoutput.f90 — first half (sites 1-7)

`src/io/swapoutput.f90` (1962L, 15 use-variables sites). Split across 2 commits to bound risk.

Sites covered in this commit (audit during impl):
- OutSwap (main driver)
- OutInc (water balance)
- OutBalSwbal (cumulative balance)
- OutDrf (drainage)
- OutAge (tracer)
- OutInco (initial conditions)
- OutMassBal (mass balance)

Each writer routes its reads to state references. Replacement table (key symbols):

| Was | Now |
|---|---|
| Mesh: numnod/dz/z/disnod/ztopcp/zbotcp/layer | `state%mesh%X` |
| Atmosphere arrays | `state%atmosphere%X` |
| Crop fields | `state%crop%X` or `state%crop%common/wofost/grass/fixed%X` |
| Nutrients | `state%nutrients%X` |
| Cumulative balance globals (cgrai/cnrai/etc.) | already on `state%atmosphere%X` |

Update callers of each writer; verify state arg threading.

### B4: swapoutput.f90 — second half (sites 8-15)

Remaining writers:
- OutTem (temperature)
- OutSnow (snow)
- OutCropFixed
- OutWofost
- OutCropGrass
- OutBma (biomass)
- OutVap (vapor)
- OutCsvFmt (CSV format helpers)

Same pattern as B3. After B3+B4, swapoutput.f90 should have NO bare `use variables` (modulo documented narrow deferrals for things like `logf` or symbols not yet on state).

### B5: cropgrowth.f90 narrow imports (10 sites)

`src/crop/cropgrowth.f90` was narrowed in GR-CROP Phase B but kept 10 narrow `use variables, only: ...` sites for config-loaded crop parameters (idev/tsumea/tbase/lrnr/etc.).

Substantial config threading required: each subroutine (CropFixed, CropWofost, CropGrass, sub-subroutines like ArableLandGerm/FacCO2/etc.) needs `config` arg (most already have it from prior arcs — verify per-subroutine).

Per-symbol routing:
- `idev/tsumea/tsumam/tbase` → `config%crop%fixed%X`
- `lrnr/lsnr/nlue/rnflv/rnfst/frnx/nmxlv` → `config%crop%wofost%nutrient%X` or `config%crop%wofost%X`
- `c_mroot/f_senes/q10_root/q10_microbial/shape_factor_rootr/specific_resp_humus` → `config%crop%fixed%X` (per the A3 finding)
- `max_resp_factor` → `config%crop%fixed%max_resp_factor` (added in A3 if missing)
- Other narrow imports — classify per-symbol.

After this, cropgrowth.f90 has NO `use variables`.

### B6: Other crop runtime files

Per-file narrow import migration:
- `src/crop/oxygenstress.f90` (4 sites) — config threading; close any residual GR-ATM deferrals
- `src/crop/rootextraction.f90` (4 sites) — config threading
- `src/crop/irrigation.f90` (1 site at line 320: mairg, irrigevent, qssdi, qssdisum, dt_SSDI_event — these may stay if irrigation-event-specific arrays aren't in state)
- `src/crop/tillage.f90` (1 site)
- `src/crop/management_soil.f90` (1 narrow site)

After this, crop runtime files have zero `use variables` clauses.

### B7: Crop init files

- `src/crop/cropfixed_init.f90` (1 narrow site at line 28)
- `src/crop/cropgrass_init.f90` (2 narrow sites at lines 34, 433)
- `src/crop/cropwofost_init.f90` (2 narrow sites at lines 42, 503)
- `src/crop/wofost_soil_parameters.f90` (1 narrow site at line 14)

Mostly config threading for `idev/tsumea/tsumam/tbase/lrnr/etc.` per A3 + B5.

### B8: timecontrol_mod.f90 (5 sites)

`src/core/timecontrol_mod.f90` has 5 `use variables` sites — likely residual narrow imports for time-related config (tcum, daycum, daynr, etc.). Most are already on `state%timecontrol`. Migrate remainders + close any deferrals.

### B9: Soil cluster

- `src/soil/soilhydraulics.f90` (4 sites)
- `src/soil/waterbalance.f90` (3 sites)
- `src/soil/soilgrid.f90` (residuals)

Config threading for soil params; some may close boundtop's deferred `swkmean/swredu/flrunon/runonarr` if threaded properly.

### B10: Drainage cluster

- `src/drainage/drainage.f90` (4 sites)
- `src/drainage/surfacewater.f90`, `divdra.f90` residuals

Config + state threading.

### B11: Small consumers

Per-file audit of remaining `use variables` sites:
- `src/boundary/boundbottom.f90` narrow deferrals (logf + SwBotb3ResVert + gwltab/qbotab/haqtab/hbotab)
- `src/boundary/boundtop.f90` narrow deferrals (logf, nird, swkmean/swredu/flrunon/runonarr, swpondmx/pondmxtab)
- `src/heat/temperature.f90`, `frozencond.f90` residuals
- `src/solute/solute.f90`, `agetracer.f90` residuals
- `src/atmosphere/meteoday.f90`, `meteodt.f90`, `et.f90`, `interception.f90` residual narrow imports
- `config_to_variables.f90` self-reads (these will retire naturally during Phase C+D)

For each, determine if the symbol is state, config, or write-only. Migrate or document deferral. By end of B11, the only `use variables` should be `config_to_variables.f90` self-reads (Phase C+D will handle them).

### B12: Phase B close marker

check-full 5/5; commit marker.

## Phase C — Strangler-fig elimination

### C1: Retire swap_mod transient buffer reads

`src/core/swap_mod.f90` lines 120-142 + 297-310: `h_init_buf`, `pondini_init_buf`, `pond_init_buf`, `tc_iyear_init_buf`, `tc_imonth_init_buf`, `tc_dt_init_buf` are transient buffers in `config_to_variables_mod`. swap_mod reads them to seed state.

Replace each read with direct `config%X` sourcing:
- `tc_iyear_init_buf` → `config%simulation%iyear` (verify path)
- `tc_imonth_init_buf` → `config%simulation%imonth`
- `tc_dt_init_buf` → `config%timecontrol%dt_default` or similar
- `pondini_init_buf` → `config%soil%water%pondini`
- `pond_init_buf` → similar
- `h_init_buf` → `config%soil%initial%h_init_array`

Verify each config path exists. If missing, add to config schema in this task.

Then remove the buffer declarations from `config_to_variables_mod` (and the `use config_to_variables_mod` line in swap_mod if it becomes empty).

### C2: Retire swap_mod swinco=3 inline block

Block at ~line 297-312 conditionally seeds atmosphere state when `swinco==3`. Per the prior arcs, the atmosphere state mirrors track legacy through simulation now. This block may be either:
- Redundant (Phase A dual-writes already cover seeding) → fully retire
- Still needed for the warm-restart case where atmosphere data must be applied differently → simplify to read from `config%X` directly

Inspect carefully + decide per-line. If retire, drop the entire `if (config%soil%swinco == 3) then ... end if` block.

### C3: Drop adapter dual-write blocks for write-only globals

In swap_mod, the `[SS-GR-ATM A6/A7/A9/A10-A12]` and `[SS-GR-CROP A14-A17]` blocks mirror legacy globals into state. As Phase D deletes those legacy globals, these dual-writes become dead code.

This is an incremental task — drops one block at a time as Phase D retires the corresponding globals. Likely several sub-commits.

### C4: Shrink config_to_variables.f90 — drop adapter writes for retired globals

For each global being retired in Phase D, the corresponding adapter line in `config_to_variables.f90` becomes dead code. Drop it.

Goal: by end of Phase D, `config_to_variables.f90` has minimal content (or is empty), enabling Phase E2's deletion.

### C5: Phase C close marker

check-full 5/5; commit marker.

## Phase D — Global retirement (incremental)

Iterative; each task targets a category of globals. After each, run VG. Compile errors surface missed readers — patch in same task.

### D1: Retire atmosphere multi-consumer globals

Symbols: `arad/atmn/atmx/ahum/awin/arai/aetr/wet/atav/tav/ad/am/daynrfirst/daynrlast/atmin7/nofd/teprrain/teprsnow/Tav/tavd/rh` plus rain CSV / meteo detail arrays.

Preconditions: B1+B2+B3+B4+B11 have migrated atmosphere consumers.

Tombstone declarations in `variables.f90`; drop zero-fills in `initialize.f90`; drop legacy populator in `config_to_variables.f90`.

### D2: Retire crop runtime globals

All crop fields that have state mirrors and no remaining bare-global readers after B5+B6+B7.

### D3: Retire nutrient globals (if any)

Verify if `variables.f90` has any nutrient declarations (most are in `Wofost_Soil_Declarations` module). If yes, retire.

### D4: Retire timecontrol globals

After B8, any remaining bare timecontrol fields retire.

### D5: Retire soil/drainage runtime globals

After B9-B10.

### D6: Retire residual scattered globals

After B11. Catch-all for everything else.

### D7: Phase D close marker — variables.f90 sanity check

`grep -nE "^[[:space:]]*\b(real\(8\)|integer|logical|character).*::" src/core/variables.f90` shows empty (or only deprecated-comments). If anything remains, identify why and either retire it or extend Phase D.

## Phase E — Final adapter deletion

### E1: Verify variables.f90 has no remaining declarations

Sanity grep. If any remain, return to Phase D.

### E2: Delete `src/core/variables.f90`

- `git rm src/core/variables.f90`
- Remove from `meson.build` legacy sources
- Remove from `tests/unit/meson.build` if listed
- Iterative build cycle — patch any remaining `use variables` lines (should be only in the about-to-delete `config_to_variables.f90` and `initialize.f90`)

### E3: Delete `src/io/toml/config_to_variables.f90`

- `git rm src/io/toml/config_to_variables.f90`
- Remove from `meson.build`
- Remove `use config_to_variables_mod` line in `swap_mod.f90` (verify the buffer reads were retired in C1)
- Iterative build cycle

### E4: Delete `src/core/initialize.f90`

- `git rm src/core/initialize.f90`
- Remove from `meson.build`
- Remove `call initialize(...)` from `swap_mod.f90` (likely nothing left to initialize after Phase D)
- Iterative build cycle

### E5: meson.build + tests/unit/meson.build cleanup

Drop entries for the 3 deleted files. Verify no broken references.

### E6: Final ADR — post-retirement state architecture

`docs/superpowers/decisions/0044-state-architecture-final.md` documenting:
- The completed strangler-fig demolition
- The two-paradigm data flow: `state%X` (runtime mutable, typed) + `config%X` (read-only after init, typed)
- The state subrecord taxonomy (mesh, soilwater, atmosphere, heat, drainage, surfacewater, solute, tillage, timecontrol, crop, nutrients)
- Lessons learned: Phase A.5 dual-write coverage; zero-resets count; controller intervention pattern
- Migration history: GR-UTILS → GR-BH → GR-ATM → GR-CROP → GR-FINAL

### E7: Final verification + arc-complete marker

- check-full 5/5 byte-for-byte
- BMI + cffi-demo passing
- pFUnit all-pass
- `wc -l src/core/variables.f90` → file not found (deleted)
- `wc -l src/io/toml/config_to_variables.f90` → file not found
- `wc -l src/core/initialize.f90` → file not found
- Arc-complete marker commit.

## Risk & Mitigation

| Risk | Mitigation |
|---|---|
| swapoutput's 15 sites cascade failures | Split into 2 commits (B3/B4); VG after each |
| Crop config threading creates ripple through 30+ subroutines | Controller intervenes inline; trust the audit |
| swap_mod swinco=3 block has subtle semantics | Inspect carefully in C2; don't retire prematurely |
| variables.f90 deletion (E2) exposes compile errors | Iterative patch cycle; expected ~5-15 patches |
| config_to_variables.f90 has self-reads (8 sites) | Will retire naturally during Phase C+D as legacy populator lines drop |
| Phase D order matters — must retire globals AFTER their readers migrate | Each Phase D task explicitly states its preconditions |
| Adapter deletion in E might surface previously-hidden init dependencies | E4 deletion is last; verify `initialize.f90` is empty (only zero-fills, no actual logic) before deleting |

## What this Arc DOES NOT Change

- BMI / cffi surfaces (`swap_bmi_mod`, `swap_capi_mod`) — continue working throughout. Their dependency is on `state` and `config`, NOT on legacy globals.
- The Python demo `tests/cffi-demo/run_ensemble.py` — unchanged.
- Physics — every value computed by every formula stays bit-identical. Only data-access paths change.
- The state subrecord schemas — they're already established; this arc doesn't add new subrecords (except maybe `max_resp_factor` in A3).

## Effort

- **Total tasks:** 35 meta-tasks (A1-A4 + B1-B12 + C1-C5 + D1-D7 + E1-E7)
- **Estimated commits:** 70-140 (each task typically 2-4 commits)
- **Estimated effort:** 15-25 days subagent-driven
- **Phases gated by check-full 5/5 byte-for-byte**

## Memory & ADR Consequences

- ADR 0044 (NEW): post-retirement state architecture — final cap on the rescue arc series.
- After arc closes, update memory `project_state_rescue_complete_2026-05-12.md` with full retirement summary.
- Retirement of `feedback_state_schema_clean_rebuild.md` — once variables.f90 is deleted and `.mod` cross-deps from `swap_modern`/`swap_legacy` boundary are addressed, this memory becomes obsolete. Evaluate at arc close.

## Files & Artifacts

- **This spec:** `docs/superpowers/specs/2026-05-14-globals-final-retirement-design.md`
- **Plan (next):** `docs/superpowers/plans/2026-05-14-globals-final-retirement.md`
- **Roadmap context:** `docs/superpowers/plans/2026-05-13-globals-retirement-roadmap.md` (Arc 9 finale)
- **Retirement inventory note (produced in A2):** `docs/superpowers/notes/2026-05-14-gr-final-retirement-inventory.md`
- **Prior arcs:** GR-UTILS (`26123ad`), GR-BH (`a65bdf3`), GR-ATM (`57b00b3`), GR-CROP (`476efba`)
