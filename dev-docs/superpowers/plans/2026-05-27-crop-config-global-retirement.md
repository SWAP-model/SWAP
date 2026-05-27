# crop_config_global Retirement Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Delete the `crop_config_global` module pointer by routing its runtime reads through `state%cfg%crop`, removing the first named multi-instance blocker in the crop subsystem.

**Architecture:** `crop_config_global` and `state%cfg%crop` are the *same object* (verified: `state%crop%init(config%crop,…)` + `state%cfg => config`, both wired in `swap_init_body` before any crop runtime call). The pointer is read-only at runtime (the `rotation_loaded` cache is written only at config-load via the real `config%`). So each `crop_config_global%X` read → `crop_cfg%X` (reusing the existing per-file associate alias) or `state%cfg%crop%X` (for associate-statement targets) is byte-identical by construction.

**Tech Stack:** Fortran 2008 (gfortran), Meson via `pixi`, byte-identical regression vs `swap420gf`.

---

## Refactoring discipline (read first)

Behavior-preserving access-path migration. The invariant is the test:

```
pixi run -e test check-fast
```
(build-linux + pFUnit + byte-compare 4 cases: hupselbrook, surfacewater, salinitystress, grassgrowth). **Every task ends with `check-fast` 4/4 clean before commit.**

> **Known pre-existing xfails — not task regressions:** `soilhysteresis` and `winter` fail against `swap420gf` (adaptive-dt desync, commit `cfc447e`), and are NOT in `check-fast`'s 4 cases. So per-task `check-fast` must be fully clean; they appear only in the final `check-full` and are expected.

**The rewrite rule (per read):**
- Files with an in-scope `associate (crop_cfg => state%cfg%crop)`: `crop_config_global%X` → `crop_cfg%X`. Verify the `crop_cfg` alias is actually in scope at the read (the associate block encloses it); if not, extend the block or use `state%cfg%crop%X`.
- `associate (alias => crop_config_global%…)` *statements* (only in `cropgrowth_helpers.f90`): rewrite the target to `state%cfg%crop%…` (an associate target cannot reference a sibling alias defined in the same `associate`).
- Remove the file's local `use crop_config_global_mod, only: crop_config_global` import once that file has zero `crop_config_global` references.

**Do NOT** change `read_crop_toml.f90` (config-load writes to `rotation_loaded` stay). **Do NOT** touch dispatch. Branch `development`; development-only commits.

**Ordering:** Tasks 1–5 rewrite readers (module still exists, still set in `crop_state_init`, just unread — compiles fine). Task 6 removes the set-site + the module file once readers are zero.

---

## Task 1: cropfixed_runtime.f90

**Files:** Modify `src/crop/cropfixed_runtime.f90` (reads ~73–80; existing alias `crop_cfg => state%cfg%crop` at `:55`; local import at `:68`).

- [ ] **Step 1: Rewrite the reads.** Grep first: `grep -n "crop_config_global" src/crop/cropfixed_runtime.f90`. For each *code* reference (not comments), replace `crop_config_global%` with `crop_cfg%`. Confirm the `associate(... crop_cfg => state%cfg%crop ...)` opened at `:55` encloses every read (it wraps the `cropfixed` body). Examples: `crop_config_global%rotation_loaded` → `crop_cfg%rotation_loaded`; `crop_config_global%rotation_fixed(crop%common%icrop)` → `crop_cfg%rotation_fixed(crop%common%icrop)`.
- [ ] **Step 2: Remove the dead import.** Delete line `:68` `use crop_config_global_mod, only: crop_config_global`.
- [ ] **Step 3: Confirm clean.** `grep -n "crop_config_global" src/crop/cropfixed_runtime.f90` → only comments (or none).
- [ ] **Step 4: Build.** `pixi run -e test build-linux` → clean.
- [ ] **Step 5: Regression.** `pixi run -e test check-fast` → 4/4 byte-identical.
- [ ] **Step 6: Commit.**
```bash
git add src/crop/cropfixed_runtime.f90
git commit -m "refactor(crop): cropfixed_runtime reads via crop_cfg, drop crop_config_global"
```

---

## Task 2: cropgrass_runtime.f90

**Files:** Modify `src/crop/cropgrass_runtime.f90` (reads ~118–131; existing alias `crop_cfg => state%cfg%crop` at `:101`; local import at `:113`).

- [ ] **Step 1: Rewrite the reads.** `grep -n "crop_config_global" src/crop/cropgrass_runtime.f90`. Replace each code `crop_config_global%` with `crop_cfg%`. NOTE line ~131 is `associate (cfg => crop_config_global%rotation_grass(crop%common%icrop))` — this is an associate *target*; rewrite it to `associate (cfg => crop_cfg%rotation_grass(crop%common%icrop))` **only if** `crop_cfg` (from the `:101` associate) is in scope and lexically *outer* to this nested associate (a nested associate target CAN reference an alias from an enclosing associate). If the `:101` associate does not enclose line 131, use `state%cfg%crop%rotation_grass(crop%common%icrop)` instead. Verify by reading the block structure.
- [ ] **Step 2: Remove the dead import.** Delete `:113` `use crop_config_global_mod, only: crop_config_global`.
- [ ] **Step 3: Confirm clean.** `grep -n "crop_config_global" src/crop/cropgrass_runtime.f90` → comments only.
- [ ] **Step 4: Build.** `pixi run -e test build-linux` → clean.
- [ ] **Step 5: Regression.** `pixi run -e test check-fast` → 4/4 (the `grassgrowth` case exercises this directly — strong test).
- [ ] **Step 6: Commit.**
```bash
git add src/crop/cropgrass_runtime.f90
git commit -m "refactor(crop): cropgrass_runtime reads via crop_cfg, drop crop_config_global"
```

---

## Task 3: cropwofost_runtime.f90

**Files:** Modify `src/crop/cropwofost_runtime.f90` (reads ~189–196; existing alias `crop_cfg => state%cfg%crop` at `:171`; local import at `:184`).

- [ ] **Step 1: Rewrite the reads.** `grep -n "crop_config_global" src/crop/cropwofost_runtime.f90`. Replace each code `crop_config_global%` with `crop_cfg%` (e.g. `crop_config_global%rotation_wofost(crop%common%icrop)` → `crop_cfg%rotation_wofost(crop%common%icrop)`). Confirm the `:171` associate encloses the reads.
- [ ] **Step 2: Remove the dead import.** Delete `:184` `use crop_config_global_mod, only: crop_config_global`.
- [ ] **Step 3: Confirm clean.** `grep -n "crop_config_global" src/crop/cropwofost_runtime.f90` → comments only.
- [ ] **Step 4: Build.** `pixi run -e test build-linux` → clean.
- [ ] **Step 5: Regression.** `pixi run -e test check-fast` → 4/4.
- [ ] **Step 6: Commit.**
```bash
git add src/crop/cropwofost_runtime.f90
git commit -m "refactor(crop): cropwofost_runtime reads via crop_cfg, drop crop_config_global"
```

---

## Task 4: cropgrowth_helpers.f90 (associate-target rewrites)

**Files:** Modify `src/crop/cropgrowth_helpers.f90` (3 associate-target reads at ~123, ~153, ~188 inside `ArableLandGerm`; local import at `:86`). These reads are `associate(prep|sow|germ => crop_config_global%rotation_wofost(crop%icrop)%…)` — there is NO `crop_cfg` alias here, and an associate target can't reference a sibling alias, so spell out `state%cfg%crop`.

- [ ] **Step 1: Rewrite the 3 associate targets.** `grep -n "crop_config_global" src/crop/cropgrowth_helpers.f90`.
  - `:123` `associate(prep => crop_config_global%rotation_wofost(crop%icrop)%preparation)` → `associate(prep => state%cfg%crop%rotation_wofost(crop%icrop)%preparation)`
  - `:153` `associate(sow => crop_config_global%rotation_wofost(crop%icrop)%sowing)` → `associate(sow => state%cfg%crop%rotation_wofost(crop%icrop)%sowing)`
  - `:188` `associate(germ => crop_config_global%rotation_wofost(crop%icrop)%germination)` → `associate(germ => state%cfg%crop%rotation_wofost(crop%icrop)%germination)`
  `state` is the `ArableLandGerm(task, tsoil, state)` argument and is in scope; `crop%icrop` is whatever alias was already used in the original target (keep it identical). Handle any additional code references the grep finds the same way (→ `state%cfg%crop%…`).
- [ ] **Step 2: Remove the dead import.** Delete `:86` `use crop_config_global_mod, only: crop_config_global`.
- [ ] **Step 3: Confirm clean.** `grep -n "crop_config_global" src/crop/cropgrowth_helpers.f90` → comments only.
- [ ] **Step 4: Build.** `pixi run -e test build-linux` → clean.
- [ ] **Step 5: Regression.** `pixi run -e test check-fast` → 4/4. (These branches are dead in the TOML pipeline per the in-code comments, so this is compile + no-change confirmation.)
- [ ] **Step 6: Commit.**
```bash
git add src/crop/cropgrowth_helpers.f90
git commit -m "refactor(crop): cropgrowth_helpers reads via state%cfg%crop, drop crop_config_global"
```

---

## Task 5: cropgrowth.f90 (largest; two nested imports)

**Files:** Modify `src/crop/cropgrowth.f90` (~16 reads at 187–210 and 356–357; existing alias `crop_cfg => state%cfg%crop` at `:89`; TWO local imports at `:174` and `:354`).

- [ ] **Step 1: Rewrite the reads.** `grep -n "crop_config_global" src/crop/cropgrowth.f90`. Replace each code `crop_config_global%` with `crop_cfg%`. Examples: `crop_config_global%rotation_loaded(state%crop%common%icrop)` → `crop_cfg%rotation_loaded(state%crop%common%icrop)`; `crop_config_global%rotation_wofost(state%crop%common%icrop)%germination` (associate target at ~:210) → `crop_cfg%rotation_wofost(state%crop%common%icrop)%germination` (the `:89` associate `crop_cfg` is outer to this nested associate, so referencing it is legal); `crop_config_global%rotation_wofost(...)%bulb%pld` (~:356) → `crop_cfg%...%bulb%pld`. Confirm the `:89` associate encloses all read sites (it wraps the `CropGrowth` body).
- [ ] **Step 2: Remove BOTH dead imports.** Delete `:174` and `:354` `use crop_config_global_mod, only: crop_config_global` (these were scoped inside case blocks; once reads use `crop_cfg`, both are dead).
- [ ] **Step 3: Confirm clean.** `grep -n "crop_config_global" src/crop/cropgrowth.f90` → comments only.
- [ ] **Step 4: Build.** `pixi run -e test build-linux` → clean.
- [ ] **Step 5: Regression.** `pixi run -e test check-fast` → 4/4.
- [ ] **Step 6: Commit.**
```bash
git add src/crop/cropgrowth.f90
git commit -m "refactor(crop): cropgrowth reads via crop_cfg, drop crop_config_global"
```

---

## Task 6: Delete the module + set-site

**Files:** Modify `src/state/crop_state.f90` (set-site `:85`, import `:52`); Delete `src/crop/crop_config_global.f90`.

Precondition: Tasks 1–5 done → the only remaining references are the set-site and import in `crop_state.f90` and the module file itself. Verify: `grep -rn "crop_config_global" src/ --include=*.f90` should show ONLY `src/state/crop_state.f90` (2 lines) and `src/crop/crop_config_global.f90`.

- [ ] **Step 1: Remove the set-site + import.** In `src/state/crop_state.f90`: delete line `:85` `crop_config_global => crop_cfg` (and its `! 4. Legacy crop_config_global pointer …` comment at `:84`), and delete the import `:52` `use crop_config_global_mod, only: crop_config_global`. NOTE: `crop_cfg` is declared `type(crop_config_t), target, intent(in)` — the `target` attribute was there for this pointer association. It is now unused; leave it (harmless) or remove `, target` if no other pointer targets it (grep `=> crop_cfg` in the file first; if none remain, dropping `target` is safe but optional — prefer leaving it to minimize the diff).
- [ ] **Step 2: Delete the module file.** `git rm src/crop/crop_config_global.f90`. Then check the Meson build list: `grep -rn "crop_config_global" meson.build src/*/meson.build 2>/dev/null` — if the file is listed in any `meson.build` sources array, remove that entry.
- [ ] **Step 3: Confirm fully gone.** `grep -rn "crop_config_global" src/ --include=*.f90` → zero matches (comments may remain in unrelated files; if any reference the deleted module by name in a way that misleads, update it).
- [ ] **Step 4: Clean rebuild.** This touches `src/state/crop_state.f90`, so per the state-schema rule do a clean build: `rm -rf builddir && pixi run -e test build-linux` → clean (incremental Meson does not propagate `.mod` deps across the swap_modern↔swap_legacy boundary).
- [ ] **Step 5: Full regression.** `pixi run -e test check-full` → 5 passed / 0 failed / 2 known xfails (`soilhysteresis`, `winter`). Confirm those are the ONLY failures.
- [ ] **Step 6: Commit.**
```bash
git add src/state/crop_state.f90 src/crop/crop_config_global.f90 meson.build
git commit -m "refactor(crop): delete crop_config_global module + set-site (multi-instance blocker retired)"
```

---

## Task 7: Spec close-out

**Files:** Modify `dev-docs/superpowers/specs/2026-05-27-crop-config-global-retirement-design.md`.

- [ ] **Step 1:** Set front-matter `status: complete`; add a one-line result note (commits, check-full result).
- [ ] **Step 2: Commit.**
```bash
git add dev-docs/superpowers/specs/2026-05-27-crop-config-global-retirement-design.md
git commit -m "docs(spec): mark crop_config_global retirement complete"
```

---

## Self-review notes (for the executor)

- **Spec coverage:** all 5 reader files (spec table) → Tasks 1–5; module + set-site deletion (spec "Then delete the global") → Task 6; success criteria (`grep` zero, file deleted, `check-full`) → Task 6 Steps 3/5. Access-form rule (crop_cfg alias vs spelled-out for associate targets) is in the discipline section and applied per task.
- **Ordering safety:** module stays defined and set through Tasks 1–5 (readers removed first); only Task 6 removes the producer + module. Every intermediate state compiles.
- **State-schema rule:** Task 6 touches `crop_state.f90` → mandatory `rm -rf builddir` before its build (Step 4).
- **Line numbers** are 2026-05-27 anchors; grep each file first (Step 1 of each task) since edits don't shift across files but the counts are approximate.
- **No new behavior** — the test is byte-identical regression; there are no unit tests to write (TDD's "failing test" maps to "regression must stay green").
