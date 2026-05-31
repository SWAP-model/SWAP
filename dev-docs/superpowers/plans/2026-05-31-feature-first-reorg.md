# Feature-First Source-Tree Reorganization — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reorganize `src/` into a feature-first-for-behavior / layer-first-for-data tree per ADR 0048, byte-identical at every commit.

**Architecture:** Pure `git mv` + build-path edits. No `.f90` *content* changes — Fortran module names are path-independent, and ninja computes module dependencies automatically, so source-list order within a target is irrelevant. Each commit moves one cohesive group and updates the two build files in lockstep: root `meson.build` (source lists + `src_inc`) and `tests/unit/meson.build` (explicit `../../src/...` paths).

**Tech Stack:** Fortran (gfortran), Meson/ninja, pixi, pFUnit.

---

## Invariants (apply to EVERY task)

- **The "test" is byte-identical regression**, not new unit tests. There is no new test code in this plan.
- After any source move you MUST `rm -rf builddir` before rebuilding — incremental Meson does not propagate `.mod` deps across the `swap_modern`↔`swap_legacy` boundary (CLAUDE.md).
- Each task ends green on: `pixi run -e test check-fast` (4/4 byte-identical) **and** a clean full build of every target (exe + `libswap_bmi.so` + `libswap_xmi.so`).
- The pFUnit test count must not change. Note the baseline in Task 0 and re-confirm each task.
- `git mv` (not `mv`) so history follows the files.
- Branch: `reorg-feature-first` off `development`. Commit per task. Do **not** touch `origin/main`.
- When editing `meson.build` source lists, substitute only **quoted full paths** (`'src/old/x.f90'` → `'src/new/x.f90'`). Comments mentioning old paths are historical record — leave them.

### The two build files

- **`meson.build`** (root): the `sources` list (→ `swap_legacy`), the three
  `modern_*_sources` lists (→ `swap_modern` / `swap_modern_xmi`), the
  `executable('swap', sources: …)` line, and the `src_inc =
  include_directories(...)` list. `src_inc` MUST NOT name a directory that does
  not exist (Meson errors at configure), so renamed/removed dirs must be dropped
  from it.
- **`tests/unit/meson.build`**: 134 explicit `../../src/<dir>/<file>.f90` paths.

---

## Task 0: Branch + baseline

**Files:** none moved.

- [ ] **Step 1: Create the branch**

```bash
cd /home/zawadzkim/Code/swap
git checkout development
git checkout -b reorg-feature-first
```

- [ ] **Step 2: Establish the green baseline**

```bash
pixi run clean && pixi run build-linux
pixi run -e test check-fast
```

Expected: build succeeds; `check-fast` reports 4/4 byte-identical PASS. Record the pFUnit `OK (N tests)` count printed by the run — call it `N_BASE`. Every subsequent task must reproduce `N_BASE`.

- [ ] **Step 3: Confirm the shared-library targets build**

```bash
ninja -C builddir swap_bmi swap_xmi
```

Expected: `libswap_bmi.so` and `libswap_xmi.so` link with no errors.

No commit (no changes yet).

---

## Task 1: Split `core/` → `driver/` + `bindings/`

Move orchestration and the C-ABI facades out of `core/`, leaving `core/` as the foundation floor. These files live in the `modern_*_sources` lists and the executable source — **not** the `sources` list.

**Files:**
- Move: `src/core/swap_mod.f90` → `src/driver/swap_mod.f90`
- Move: `src/core/swap_main.f90` → `src/driver/swap_main.f90`
- Move: `src/core/swap_ensemble_mod.f90` → `src/driver/swap_ensemble_mod.f90`
- Move: `src/core/swap_capi_mod.f90` → `src/bindings/swap_capi_mod.f90`
- Move: `src/core/swap_bmi_mod.f90` → `src/bindings/swap_bmi_mod.f90`
- Move: `src/core/swap_xmi_mod.f90` → `src/bindings/swap_xmi_mod.f90`
- Move: `src/core/bmi_constants_mod.f90` → `src/bindings/bmi_constants_mod.f90`
- Modify: `meson.build` (`modern_core_sources`, `modern_sources`, `modern_xmi_sources`, `executable(...)`, `src_inc`)

- [ ] **Step 1: Move the files**

```bash
cd /home/zawadzkim/Code/swap
mkdir -p src/driver src/bindings
git mv src/core/swap_mod.f90           src/driver/swap_mod.f90
git mv src/core/swap_main.f90          src/driver/swap_main.f90
git mv src/core/swap_ensemble_mod.f90  src/driver/swap_ensemble_mod.f90
git mv src/core/swap_capi_mod.f90      src/bindings/swap_capi_mod.f90
git mv src/core/swap_bmi_mod.f90       src/bindings/swap_bmi_mod.f90
git mv src/core/swap_xmi_mod.f90       src/bindings/swap_xmi_mod.f90
git mv src/core/bmi_constants_mod.f90  src/bindings/bmi_constants_mod.f90
```

- [ ] **Step 2: Update `meson.build` source lists**

```bash
sed -i \
  -e "s|'src/core/swap_mod.f90'|'src/driver/swap_mod.f90'|" \
  -e "s|'src/core/swap_ensemble_mod.f90'|'src/driver/swap_ensemble_mod.f90'|" \
  -e "s|'src/core/bmi_constants_mod.f90'|'src/bindings/bmi_constants_mod.f90'|" \
  -e "s|'src/core/swap_capi_mod.f90'|'src/bindings/swap_capi_mod.f90'|" \
  -e "s|'src/core/swap_bmi_mod.f90'|'src/bindings/swap_bmi_mod.f90'|" \
  -e "s|'src/core/swap_xmi_mod.f90'|'src/bindings/swap_xmi_mod.f90'|" \
  -e "s|sources: 'src/core/swap_main.f90'|sources: 'src/driver/swap_main.f90'|" \
  meson.build
```

- [ ] **Step 3: Add the new dirs to `src_inc`**

In `meson.build`, edit the `src_inc = include_directories(...)` call to append `'src/driver', 'src/bindings'` to the list. After editing, the call should include those two entries alongside the existing ones.

- [ ] **Step 4: Verify no stale `core/` paths remain for moved files**

```bash
grep -nE "src/core/(swap_mod|swap_main|swap_ensemble_mod|swap_capi_mod|swap_bmi_mod|swap_xmi_mod|bmi_constants_mod)\.f90" meson.build tests/unit/meson.build
```

Expected: **no output** (all references updated; `tests/unit/meson.build` does not reference these app/facade files).

- [ ] **Step 5: Clean rebuild + check-fast + lib targets**

```bash
pixi run clean && pixi run build-linux
ninja -C builddir swap_bmi swap_xmi
pixi run -e test check-fast
```

Expected: full build + both `.so` link; `check-fast` 4/4 byte-identical; pFUnit count == `N_BASE`.

- [ ] **Step 6: Commit**

```bash
git add -A
git commit -m "refactor(tree): split core/ into core/ (leaves) + driver/ + bindings/

Move orchestration (swap_mod, swap_main, swap_ensemble_mod) to src/driver/
and the C-ABI facades (swap_capi/bmi/xmi, bmi_constants) to src/bindings/.
core/ now holds only foundation leaves. Pure git mv + meson path edits; no
use-statement changes. Byte-identical (ADR 0048).

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 2: Rename `soil/` → `soilwater/`, absorb `boundary/`

**Files:**
- Move: `src/soil/` → `src/soilwater/` (whole dir, incl. `dormant/`)
- Move: `src/boundary/boundtop.f90` → `src/soilwater/boundtop.f90`
- Move: `src/boundary/boundbottom.f90` → `src/soilwater/boundbottom.f90`
- Modify: `meson.build` (`sources`, `src_inc`), `tests/unit/meson.build`

- [ ] **Step 1: Move the files**

```bash
cd /home/zawadzkim/Code/swap
git mv src/soil src/soilwater
git mv src/boundary/boundtop.f90    src/soilwater/boundtop.f90
git mv src/boundary/boundbottom.f90 src/soilwater/boundbottom.f90
rmdir src/boundary
```

- [ ] **Step 2: Update `meson.build` `sources`**

```bash
sed -i \
  -e "s|'src/soil/waterbalance.f90'|'src/soilwater/waterbalance.f90'|" \
  -e "s|'src/soil/soilhydraulics.f90'|'src/soilwater/soilhydraulics.f90'|" \
  -e "s|'src/soil/WC_K_models_04_11.f90'|'src/soilwater/WC_K_models_04_11.f90'|" \
  -e "s|'src/boundary/boundtop.f90'|'src/soilwater/boundtop.f90'|" \
  -e "s|'src/boundary/boundbottom.f90'|'src/soilwater/boundbottom.f90'|" \
  meson.build
```

- [ ] **Step 3: Update `src_inc`**

In `meson.build`'s `src_inc` list: replace `'src/soil'` with `'src/soilwater'` and **remove** `'src/boundary'`.

- [ ] **Step 4: Update `tests/unit/meson.build`**

```bash
sed -i \
  -e "s|\.\./\.\./src/soil/|../../src/soilwater/|g" \
  -e "s|\.\./\.\./src/boundary/|../../src/soilwater/|g" \
  tests/unit/meson.build
```

- [ ] **Step 5: Verify no stale paths remain**

```bash
grep -nE "src/soil/|src/boundary/" meson.build tests/unit/meson.build | grep -vE "^\s*#|src/soilwater"
```

Expected: **no output** (remaining `src/soil*` hits, if any, are `src/soilwater/...` or comments).

- [ ] **Step 6: Clean rebuild + verify**

```bash
pixi run clean && pixi run build-linux
ninja -C builddir swap_bmi swap_xmi
pixi run -e test check-fast
```

Expected: builds; `check-fast` 4/4 byte-identical; pFUnit == `N_BASE`.

- [ ] **Step 7: Commit**

```bash
git add -A
git commit -m "refactor(tree): rename soil/ -> soilwater/, fold in boundary/

boundtop/boundbottom are soil-water boundary conditions, not a separate
domain. Pure git mv + meson path edits. Byte-identical (ADR 0048).

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 3: Dissolve `utils/`

`utils/` was never one thing: two true leaf-math files → `core/`; two state-touching files → their feature folders (killing the `utils → state` inversion).

**Files:**
- Move: `src/utils/arrayutils.f90` → `src/core/arrayutils.f90`
- Move: `src/utils/numericalsolvers.f90` → `src/core/numericalsolvers.f90`
- Move: `src/utils/soilhydraulicsutils.f90` → `src/soilwater/soilhydraulicsutils.f90`
- Move: `src/utils/surfacewaterutils.f90` → `src/drainage/surfacewaterutils.f90`
- Modify: `meson.build` (`sources`, `src_inc`), `tests/unit/meson.build`

> Depends on Task 2 (needs `src/soilwater/` to exist).

- [ ] **Step 1: Move the files**

```bash
cd /home/zawadzkim/Code/swap
git mv src/utils/arrayutils.f90         src/core/arrayutils.f90
git mv src/utils/numericalsolvers.f90   src/core/numericalsolvers.f90
git mv src/utils/soilhydraulicsutils.f90 src/soilwater/soilhydraulicsutils.f90
git mv src/utils/surfacewaterutils.f90  src/drainage/surfacewaterutils.f90
rmdir src/utils
```

- [ ] **Step 2: Update `meson.build` `sources`**

```bash
sed -i \
  -e "s|'src/utils/arrayutils.f90'|'src/core/arrayutils.f90'|" \
  -e "s|'src/utils/numericalsolvers.f90'|'src/core/numericalsolvers.f90'|" \
  -e "s|'src/utils/soilhydraulicsutils.f90'|'src/soilwater/soilhydraulicsutils.f90'|" \
  -e "s|'src/utils/surfacewaterutils.f90'|'src/drainage/surfacewaterutils.f90'|" \
  meson.build
```

- [ ] **Step 3: Update `src_inc`**

In `meson.build`'s `src_inc` list: **remove** `'src/utils'`.

- [ ] **Step 4: Update `tests/unit/meson.build`**

```bash
sed -i \
  -e "s|\.\./\.\./src/utils/arrayutils.f90|../../src/core/arrayutils.f90|g" \
  -e "s|\.\./\.\./src/utils/numericalsolvers.f90|../../src/core/numericalsolvers.f90|g" \
  -e "s|\.\./\.\./src/utils/soilhydraulicsutils.f90|../../src/soilwater/soilhydraulicsutils.f90|g" \
  -e "s|\.\./\.\./src/utils/surfacewaterutils.f90|../../src/drainage/surfacewaterutils.f90|g" \
  tests/unit/meson.build
```

- [ ] **Step 5: Verify no stale paths remain**

```bash
grep -nE "src/utils/" meson.build tests/unit/meson.build | grep -vE "^\s*#"
```

Expected: **no output**.

- [ ] **Step 6: Clean rebuild + verify**

```bash
pixi run clean && pixi run build-linux
ninja -C builddir swap_bmi swap_xmi
pixi run -e test check-fast
```

Expected: builds; `check-fast` 4/4 byte-identical; pFUnit == `N_BASE`.

- [ ] **Step 7: Commit**

```bash
git add -A
git commit -m "refactor(tree): dissolve utils/

Leaf math (arrayutils, numericalsolvers) -> core/; state-operating helpers
(soilhydraulicsutils -> soilwater/, surfacewaterutils -> drainage/). Removes
the utils->state layering inversion. Byte-identical (ADR 0048).

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 4: Extract `timecontrol/` from `core/`

`timecontrol_mod.f90` is a compute subsystem (pairs with `state/timecontrol_state`). It lives in the `modern_core_sources` list.

**Files:**
- Move: `src/core/timecontrol_mod.f90` → `src/timecontrol/timecontrol_mod.f90`
- Modify: `meson.build` (`modern_core_sources`, `src_inc`), `tests/unit/meson.build` (if referenced)

- [ ] **Step 1: Move the file**

```bash
cd /home/zawadzkim/Code/swap
mkdir -p src/timecontrol
git mv src/core/timecontrol_mod.f90 src/timecontrol/timecontrol_mod.f90
```

- [ ] **Step 2: Update `meson.build`**

```bash
sed -i "s|'src/core/timecontrol_mod.f90'|'src/timecontrol/timecontrol_mod.f90'|" meson.build
```

In `src_inc`: append `'src/timecontrol'`.

- [ ] **Step 3: Update `tests/unit/meson.build` if it references the file**

```bash
sed -i "s|\.\./\.\./src/core/timecontrol_mod.f90|../../src/timecontrol/timecontrol_mod.f90|g" tests/unit/meson.build
```

(The `sed` is a no-op if the path is absent — safe either way.)

- [ ] **Step 4: Verify no stale paths remain**

```bash
grep -nE "src/core/timecontrol_mod.f90" meson.build tests/unit/meson.build | grep -vE "^\s*#"
```

Expected: **no output**.

- [ ] **Step 5: Clean rebuild + verify**

```bash
pixi run clean && pixi run build-linux
ninja -C builddir swap_bmi swap_xmi
pixi run -e test check-fast
```

Expected: builds; `check-fast` 4/4 byte-identical; pFUnit == `N_BASE`.

- [ ] **Step 6: Commit**

```bash
git add -A
git commit -m "refactor(tree): extract timecontrol/ from core/

timecontrol_mod is a compute subsystem (pairs with state/timecontrol_state),
not a foundation leaf. Byte-identical (ADR 0048).

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 5: Restructure `crop/` into `fixed/` `grass/` `wofost/`

Split the three growth models into subfolders; move all WOFOST nutrient dynamics under `wofost/`. The dispatcher and cross-mode files (`cropgrowth`, `cropgrowth_helpers`, `rootextraction`, `oxygenstress`) and — for now — `irrigation`, `tillage` stay at `crop/` root. `dormant/` stays.

**Files:**
- Move → `src/crop/fixed/`: `cropfixed_init.f90`, `cropfixed_runtime.f90`
- Move → `src/crop/grass/`: `cropgrass_init.f90`, `cropgrass_runtime.f90`
- Move → `src/crop/wofost/`: `cropwofost_init.f90`, `cropwofost_runtime.f90`, `wofostnut.f90`, `wofost_soil_declarations.f90`, `wofost_soil_interface.f90`, `wofost_soil_parameters.f90`, `wofost_soil_rateconstants.f90`, `wofost_soil_orgmatn.f90`, `wofost_soil_watern.f90`, `wofost_soil_amendments.f90`, `wofost_soil_cropresidues.f90`, `wofost_soil_balancecheck.f90`, `management_soil.f90`
- Modify: `meson.build` (`sources`, `src_inc`), `tests/unit/meson.build`

- [ ] **Step 1: Move the files**

```bash
cd /home/zawadzkim/Code/swap
mkdir -p src/crop/fixed src/crop/grass src/crop/wofost

git mv src/crop/cropfixed_init.f90    src/crop/fixed/cropfixed_init.f90
git mv src/crop/cropfixed_runtime.f90 src/crop/fixed/cropfixed_runtime.f90

git mv src/crop/cropgrass_init.f90    src/crop/grass/cropgrass_init.f90
git mv src/crop/cropgrass_runtime.f90 src/crop/grass/cropgrass_runtime.f90

git mv src/crop/cropwofost_init.f90    src/crop/wofost/cropwofost_init.f90
git mv src/crop/cropwofost_runtime.f90 src/crop/wofost/cropwofost_runtime.f90
git mv src/crop/wofostnut.f90          src/crop/wofost/wofostnut.f90
git mv src/crop/management_soil.f90    src/crop/wofost/management_soil.f90
for f in declarations interface parameters rateconstants orgmatn watern amendments cropresidues balancecheck; do
  git mv "src/crop/wofost_soil_${f}.f90" "src/crop/wofost/wofost_soil_${f}.f90"
done
```

- [ ] **Step 2: Update `meson.build` `sources`**

```bash
sed -i \
  -e "s|'src/crop/cropfixed_init.f90'|'src/crop/fixed/cropfixed_init.f90'|" \
  -e "s|'src/crop/cropfixed_runtime.f90'|'src/crop/fixed/cropfixed_runtime.f90'|" \
  -e "s|'src/crop/cropgrass_init.f90'|'src/crop/grass/cropgrass_init.f90'|" \
  -e "s|'src/crop/cropgrass_runtime.f90'|'src/crop/grass/cropgrass_runtime.f90'|" \
  -e "s|'src/crop/cropwofost_init.f90'|'src/crop/wofost/cropwofost_init.f90'|" \
  -e "s|'src/crop/cropwofost_runtime.f90'|'src/crop/wofost/cropwofost_runtime.f90'|" \
  -e "s|'src/crop/wofostnut.f90'|'src/crop/wofost/wofostnut.f90'|" \
  -e "s|'src/crop/management_soil.f90'|'src/crop/wofost/management_soil.f90'|" \
  -e "s|'src/crop/wofost_soil_declarations.f90'|'src/crop/wofost/wofost_soil_declarations.f90'|" \
  -e "s|'src/crop/wofost_soil_interface.f90'|'src/crop/wofost/wofost_soil_interface.f90'|" \
  -e "s|'src/crop/wofost_soil_parameters.f90'|'src/crop/wofost/wofost_soil_parameters.f90'|" \
  -e "s|'src/crop/wofost_soil_rateconstants.f90'|'src/crop/wofost/wofost_soil_rateconstants.f90'|" \
  -e "s|'src/crop/wofost_soil_orgmatn.f90'|'src/crop/wofost/wofost_soil_orgmatn.f90'|" \
  -e "s|'src/crop/wofost_soil_watern.f90'|'src/crop/wofost/wofost_soil_watern.f90'|" \
  -e "s|'src/crop/wofost_soil_amendments.f90'|'src/crop/wofost/wofost_soil_amendments.f90'|" \
  -e "s|'src/crop/wofost_soil_cropresidues.f90'|'src/crop/wofost/wofost_soil_cropresidues.f90'|" \
  -e "s|'src/crop/wofost_soil_balancecheck.f90'|'src/crop/wofost/wofost_soil_balancecheck.f90'|" \
  meson.build
```

- [ ] **Step 3: Update `src_inc`**

In `meson.build`'s `src_inc` list: append `'src/crop/fixed', 'src/crop/grass', 'src/crop/wofost'`.

- [ ] **Step 4: Update `tests/unit/meson.build`**

```bash
sed -i \
  -e "s|\.\./\.\./src/crop/cropfixed_init.f90|../../src/crop/fixed/cropfixed_init.f90|g" \
  -e "s|\.\./\.\./src/crop/cropfixed_runtime.f90|../../src/crop/fixed/cropfixed_runtime.f90|g" \
  -e "s|\.\./\.\./src/crop/cropgrass_init.f90|../../src/crop/grass/cropgrass_init.f90|g" \
  -e "s|\.\./\.\./src/crop/cropgrass_runtime.f90|../../src/crop/grass/cropgrass_runtime.f90|g" \
  -e "s|\.\./\.\./src/crop/cropwofost_init.f90|../../src/crop/wofost/cropwofost_init.f90|g" \
  -e "s|\.\./\.\./src/crop/cropwofost_runtime.f90|../../src/crop/wofost/cropwofost_runtime.f90|g" \
  -e "s|\.\./\.\./src/crop/wofostnut.f90|../../src/crop/wofost/wofostnut.f90|g" \
  -e "s|\.\./\.\./src/crop/management_soil.f90|../../src/crop/wofost/management_soil.f90|g" \
  -e "s|\.\./\.\./src/crop/wofost_soil_\([a-z]*\)\.f90|../../src/crop/wofost/wofost_soil_\1.f90|g" \
  tests/unit/meson.build
```

- [ ] **Step 5: Verify only the intended `crop/` files remain at root**

```bash
# Source-list references that should NOT have moved (still at crop/ root):
grep -nE "src/crop/(cropgrowth|cropgrowth_helpers|rootextraction|oxygenstress|irrigation|tillage)\.f90" meson.build
# Any leftover flat references to moved files (should be empty):
grep -nE "src/crop/(cropfixed|cropgrass|cropwofost|wofostnut|management_soil|wofost_soil)[A-Za-z_]*\.f90" meson.build tests/unit/meson.build | grep -vE "src/crop/(fixed|grass|wofost)/" | grep -vE "^\s*#"
```

Expected: first grep lists the 6 stay-put files; second grep produces **no output**.

- [ ] **Step 6: Clean rebuild + verify**

```bash
pixi run clean && pixi run build-linux
ninja -C builddir swap_bmi swap_xmi
pixi run -e test check-fast
```

Expected: builds; `check-fast` 4/4 byte-identical; pFUnit == `N_BASE`.

- [ ] **Step 7: Commit**

```bash
git add -A
git commit -m "refactor(tree): split crop/ into fixed/ grass/ wofost/

Three growth models get their own subfolders; all WOFOST nutrient dynamics
(wofost_soil_*, wofostnut, management_soil) move under crop/wofost/. Dispatcher
and cross-mode files (cropgrowth, rootextraction, oxygenstress) plus irrigation
and tillage stay at crop/ root for now. Byte-identical (ADR 0048).

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 6: Documentation + full-suite gate

**Files:**
- Modify: `dev-docs/adr/0048-feature-first-source-reorg.md` (status + Verification)
- Modify: `CLAUDE.md` (Architecture section — directory names)
- Modify: `src/README.md` (if it describes the tree)

- [ ] **Step 1: Run the full regression suite**

```bash
pixi run -e test check-full
```

Expected: 5/5 non-xfail byte-identical (pre-existing `winter` + `soilhysteresis` xfails unchanged). Record the result for the ADR.

- [ ] **Step 2: Flip ADR 0048 to Accepted and fill Verification**

In `dev-docs/adr/0048-feature-first-source-reorg.md`:
- Change `**Status:** Proposed (2026-05-31) — flips to Accepted …` to `**Status:** Accepted (2026-05-31)`.
- Replace the "Verification" body with the actual results: `check-fast` 4/4, `check-full` 5/5 byte-identical, pFUnit `N_BASE` unchanged, exe + `libswap_bmi.so` + `libswap_xmi.so` all link.

- [ ] **Step 3: Update `CLAUDE.md` Architecture section**

In `CLAUDE.md`, update the "Architecture (current)" and "Where to look" wording so the directory references match the new tree (`src/core/` = leaves only; orchestration in `src/driver/`; C-ABI in `src/bindings/`; `src/soilwater/`; `src/timecontrol/`; `crop/{fixed,grass,wofost}/`; `utils/` and `boundary/` no longer exist). Do not change the non-negotiables.

- [ ] **Step 4: Update `src/README.md` if it enumerates folders**

```bash
grep -nE "boundary|utils|src/core|soil/" src/README.md
```

If it describes the old layout, edit it to match the new tree. If it does not mention directories, skip.

- [ ] **Step 5: Final clean build + commit**

```bash
pixi run clean && pixi run build-linux
ninja -C builddir swap_bmi swap_xmi
pixi run -e test check-fast
git add -A
git commit -m "docs(tree): finalize ADR 0048 + update CLAUDE.md/README for new layout

check-full 5/5 byte-identical. ADR 0048 accepted.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Self-Review notes

- **Spec coverage:** every ADR 0048 move maps to a task — `core` split (T1),
  `soil`→`soilwater`+`boundary` (T2), `utils` dissolution (T3), `timecontrol`
  (T4), `crop` split incl. WOFOST nutrients (T5), docs + full gate (T6).
  `config/`/`state/`/`io/`/`error/`/`validation/` are intentionally untouched.
- **Out-of-scope items** (`state→io` inversion, `wofost_soil_interface` hidden
  `SAVE`, promoting irrigation/tillage) are explicitly NOT in this plan — they
  are separate follow-ups recorded in the ADR.
- **Ordering dependency:** T3 must follow T2 (`soilhydraulicsutils` needs
  `src/soilwater/`). All other tasks are order-independent but are sequenced
  core→features→crop→docs for readable history.
- **Risk control:** every task is `git mv` + quoted-path edits + clean rebuild +
  byte-identical `check-fast`; any task that breaks the build or diverges output
  is reverted in isolation without affecting the others.
