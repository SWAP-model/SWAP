---
title: "Phase 4e — Unified Error Handling + Phase 4f Prep"
author: Mateusz Zawadzki
date: 2026-04-27
status: draft
---

# Phase 4e: Unified Error Handling + Phase 4f Prep

A small, mostly mechanical phase that completes the error-handling story (`error_collection_t` everywhere) and prepares the runway for Phase 4f's strangler-fig replacement of `readswap()`. Phase 4e ships three things: error rewrite across surviving files, an adapter audit doc that tells Phase 4f what `config_to_variables(config)` must assign, and a macropore-case TOML so all 6 cases run via the new path before we rewire the main driver.

## Background

Phases 4a–4d delivered the typed config + TOML loader pipeline and proved field-by-field parity against `variables` globals for 5 of 6 regression cases. The next step is to actually wire the new pipeline into `swap.f90` execution (Phase 4f) — but two prep items would make 4f's work safer:

1. **Error idiom inconsistency.** New code paths (`src/config/`, `src/io/toml/`, `src/validation/`) use `error_collection_t%append + abort_if_fatal`. Legacy paths (`src/io/readmeteo.f90`, `src/io/swapoutput.f90`, every physics module) use `call fatalerr(routine, message)` which prompts stdin → fatal Fortran runtime error in non-TTY contexts (this is what the parity-helper FOPENG cleanup mitigates). After 4f, the new TOML reader is the user-facing failure surface, but the surviving legacy code still calls `fatalerr` from physics paths. Unifying the two means all errors flow through one channel — testable, redirectable, recoverable.

2. **No comprehensive `variables%` ↔ `swap_config_t` map.** Phase 4d added bottom_boundary, heat, irrigation, solute configs but didn't audit completeness. There may be `variables%foo` fields that no `*_config_t` covers — Phase 4f's `config_to_variables(config)` adapter would either need to leave those at default or surface them as a gap. An audit doc is cheap insurance.

A third item is independent: case 3 (3.macroporeflow) still has no TOML. It's been deferred since Phase 4c-a, with a short-term rationale (macropore is a special physics path and its `.crp` shape differs). Phase 4f's verification gate is "all 6 cases green via the new path" — that's only achievable if macropore has TOML. Doing it as part of 4e (alongside the other prep items) keeps 4f's scope tight.

## Scope

### In

- **Error replacement (Part A)**: replace `call fatalerr(routine, msg)` with `error_collection_t%append(...)` + `abort_if_fatal()` in every file that will SURVIVE Phases 4f and 4g. Specifically:
  - SKIP `src/io/readswap.f90` (4f deletes it).
  - SKIP the crop sub-reader sections (`readcropfixed`, `readwofost`, `readgrass`) inside `src/io/readswap.f90` (4g deletes them — these are inside readswap.f90 so the file-level skip covers it).
  - REPLACE everywhere else that has `call fatalerr`: `src/io/swapoutput.f90` (44), `src/crop/cropgrowth.f90` (19, but careful — only the parts that aren't crop sub-readers), `src/crop/tillage.f90` (16), `src/io/swap_csv_output.f90` (11), `src/io/readmeteo.f90` (9), and the physics modules in the long tail (~30 files, mostly 1–7 calls each).
  - The legacy `subprojects/ttutil/src/fatalerr.f90` itself stays untouched (vendored).
- **Adapter audit (Part B)**: produce `docs/phase-4f-config-to-variables-audit.md`. For every `variables%` field referenced from `swap.f90` execution paths (Initialize, dynamic loop, closure), classify it as one of:
  - **(C) covered**: populated by some `*_config_t` field today; document the path.
  - **(R) runtime**: state computed during simulation, never read from input — stays in `variables`, the adapter doesn't touch it.
  - **(G) gap**: legacy reader populates it, no `*_config_t` covers it, but execution paths read it → must be added before Phase 4f. Likely to be empty or near-empty after 4d, but verify.
- **Macropore case TOML (Part C)**: author `tests/swap-cases/toml/3.macroporeflow/swap.toml` + `swap.dra.toml` + any `.crp.toml` files referenced. Plus `test_macroporeflow_parity.pf` mirroring the 4d-extended pattern. NO new `*_config_t` is needed — macropore physics is unchanged; the file is just a TOML transcription of the existing `.swp.template`.
- **Closeout (Part D)**: schema doc bump (one note about the audit doc), coverage rebaseline if numbers move, baseline log, tag `rescue/phase-4e-error-prep`.

### Out

- **Replacing `fatalerr` inside `readswap.f90`** — those calls disappear in 4f.
- **Replacing `fatalerr` inside the crop sub-readers in cropgrowth.f90** — disappear in 4g.
- **Phase 4f itself** (the strangler-fig). Spec sketch is in this doc; detailed plan comes after 4e closes.
- **Phase 4g itself** (the .dra and .crp strangler-fig). Spec sketch only.
- **`fatalerr` deletion**: even after 4e, `subprojects/ttutil/src/fatalerr.f90` still exists because the file is vendored and the implementation is referenced from `readswap.f90` (until 4f deletes that). Phase 4f deletion is what eventually orphans it.
- **Logging redirection**: `error_collection_t` summarize prints to stdout; `swap_log` writes to the .log file. Bridging them is Phase 4f territory.

## Critical decisions

### D1 — Error replacement signature pattern

Today's code:
```fortran
if (idsl < 0 .or. idsl > 2) then
   call fatalerr('readwofost', 'illegal value for idsl')
end if
```

After 4e:
```fortran
if (idsl < 0 .or. idsl > 2) then
   call errors%append(ERR_VALIDATION_OUT_OF_RANGE, 'illegal value for idsl', 'readwofost')
   call errors%abort_if_fatal()
end if
```

Two side effects worth flagging:
- The replacement requires the surrounding subroutine to have access to an `errors` argument. Most physics subroutines today take no error argument. **The mechanical fix:** add `type(error_collection_t), intent(inout) :: errors` as a new argument and propagate it up the call chain. This is a bigger ripple than just rewriting the call sites — every caller of every modified subroutine needs to pass an `errors` it owns.
- For routines deep in the call stack where threading `errors` would touch dozens of files, an alternative is a module-level `errors` instance in `error_mod` (a singleton) that any subroutine can append to. Cheaper signature-wise, slightly less testable. **Decision: use the singleton pattern for physics paths, threaded `errors` for I/O and configuration paths.** The 4e plan groups commits accordingly.

### D2 — `abort_if_fatal()` placement

Where does the abort fire? Two options:
- **Inline at every call site**: `errors%append; errors%abort_if_fatal()`. Mirrors `fatalerr`'s behavior 1:1.
- **Deferred to a synchronization point**: `errors%append` only; `abort_if_fatal()` fires at the top of `swap.f90`'s main loop iterations, between `Initialize` and the dynamic loop, etc.

**Decision: inline at every call site for Phase 4e.** Behavior parity with `fatalerr` is the priority; refactoring the abort discipline is post-Phase-4 work.

### D3 — Audit doc format (Part B)

Single markdown file `docs/phase-4f-config-to-variables-audit.md`. Structure:

```
# `variables%` field audit for Phase 4f config_to_variables adapter

Each row: legacy global, status (C/R/G), TOML path (if C), notes.

| variable | status | source / target | notes |
|---|---|---|---|
| project | C | config%general%project | string copy |
| swscre | C | config%general%swscre | int copy |
| tstart | C | config%simulation%tstart | already days-since-1900 from parse_date_to_days1900 |
| t1900 | R | runtime — TimeControl(1) initializes | not in adapter |
| swfrost | C | config%bottom_boundary?? | NEEDS CHECK — see notes |
| ... | | | |
```

Method: grep `^use variables, only:` across `src/` and produce the union of all field names referenced by physics/I/O subroutines. Cross-check each against the `*_config_t` field list. Anything that isn't covered is either (R) genuinely runtime or (G) a gap.

### D4 — Macropore TOML scope

Macropore's `.swp.template` references `.crp` and possibly `.dra` files. The `.crp` is the same shape as other cases (likely type 1 fixed). Macropore-specific input (the macropore physics parameters) lives in the `.swp` itself under a section that today's `swap_config_t` may not cover.

**Risk:** authoring macropore TOML may surface that `[macropore]` section coverage is missing — i.e., a Phase 4d-pluss gap.

**Decision:** if a macropore-specific section is missing from the schema, scope down Part C to "TOML for everything except macropore-physics-specific keys; load+validate clean; parity test asserts only the schema-covered fields." Defer the missing `macropore_config_t` to a separate task (Phase 4f-prep or later). The check-full regression will still pass because macropore physics runs through legacy `readswap.f90` until 4f anyway.

### D5 — `fatalerr` shim during transition

Some physics subroutines may have so many call sites that threading `errors` everywhere is painful. For those, a **shim**: keep the routine signature unchanged, but add a thin wrapper `subroutine fatalerr_collected(routine, msg)` that internally calls `error_mod`'s singleton append. Replace `call fatalerr` with `call fatalerr_collected` mechanically. This preserves the legacy abort behavior while routing through the new collection. **Use this for `src/utils/numericalsolvers.f90`, `src/soil/sptabulated.f90`, and any deep-call-chain physics file.** Don't use it for I/O or configuration code where threading `errors` is clean.

## Critical files

### Modify (Phase 4e)

- `src/error/error.f90` — possibly add `fatalerr_collected` shim per D5.
- ~30 files under `src/` containing `call fatalerr` (excluding `readswap.f90` and crop sub-readers).
- `tests/unit/error/test_error.pf` — add a test for `fatalerr_collected` if added.
- `tests/swap-cases/toml/3.macroporeflow/` — new TOML files (in submodule).
- `tests/unit/io/toml/test_macroporeflow_parity.pf` — new parity test.
- `tests/unit/io/toml/parity_helpers.f90` — add `load_both_for_macroporeflow` helper.
- `tests/unit/meson.build`, `tests/unit/testSuites.inc` — register new test.
- `docs/phase-4f-config-to-variables-audit.md` — new audit doc (Part B output).
- `docs/coverage-baseline.md` — Phase 4e section.
- `docs/configuration-schema.md` — small note pointing to audit doc.
- `tests/regression/baselines/phase-4e-error-prep.log` — new.

### Create

(All in the modify list above.)

### Delete

None in 4e. Phase 4f deletes `readswap.f90` and orphans `subprojects/ttutil/src/fatalerr.f90` (vendored; we don't touch).

## Phase 4f sketch (informational; detailed plan after 4e closes)

**Goal:** Replace `call ReadSwap()` in `swap.f90:124` with the new TOML pipeline. Verify check-full 6/6 green via the new path. Tag `rescue/phase-4f-strangler-swp`.

**Steps:**

1. Author `src/io/toml/config_to_variables.f90` with `subroutine config_to_variables(config)` that copies every (C)-classified field from the audit. Each field gets a one-line assignment.
2. In `swap.f90`, replace:
   ```
   call ReadSwap()
   ```
   with:
   ```
   call load_swap_config('swap.toml', config, errors)
   call config%validate(errors)
   call config%finalize(errors)
   call errors%abort_if_fatal()
   call config_to_variables(config)
   ```
3. Per-case: rename or remove `swap_linux.swp.template` from each case dir, OR keep it and run an integration test that demonstrates `swap.toml` is what the binary reads.
4. Build, run check-full. Iterate.
5. Once green: `git rm src/io/readswap.f90`, update meson.
6. Add `bbcfil` deprecation warning in `read_bottom_boundary_toml.f90` per user note.
7. Tag.

**Risks:** Many. The audit's job (Part B) is to surface them ahead of time.

## Phase 4g sketch (informational)

**Goal:** Replace the per-rotation-entry `readwofost` / `readcropfixed` / `readgrass` calls in `cropgrowth.f90` with calls into the existing `*_config_t` data. Verify 6/6 green. Tag `rescue/phase-4g-strangler-crop-dra`.

**Steps:** mirror 4f pattern but for crop and drainage sub-readers. The cross-file TOML loading already works; what's missing is `config_to_variables_per_rotation_entry(rotation_idx)` that copies the right `cropwofost_config_t` / `cropfixed_config_t` / `cropgrass_config_t` slice into the legacy globals.

**Note on bottom-up vs top-down:** could also be done before 4f. Doing 4f first proves the strangler-fig pattern works on the contained `.swp` reader; 4g extends it. Reverse order is harder because the .swp reader sets up state that the crop sub-readers depend on.

## Verification gates

Phase 4e per-Part:
- After Part A: `pixi run -e test test-pfunit` green; `pixi run -e test check-full` 6/6.
- After Part B: audit doc exists, every `variables%` field classified.
- After Part C: macropore parity test green; smoke load test for case 3 green.
- After Part D: tagged.

Final: `git diff rescue/phase-4d-remaining-configs..HEAD -- src/` shows error-replacement diffs and zero physics-algorithm changes. (The `git diff` is the audit-trail proof that 4e didn't accidentally change physics.)

## File count summary

- Modify: ~30 source files (error replacement) + 8 test/doc files.
- Create: ~8 files (3 macropore TOMLs in submodule, parity test, audit doc, baseline log, schema doc snippet, helper extension).
- LoC delta: ~+1,000 / ~-500 (replace lines roughly 2:1; net +500 because error_collection_t calls are slightly verbose).
