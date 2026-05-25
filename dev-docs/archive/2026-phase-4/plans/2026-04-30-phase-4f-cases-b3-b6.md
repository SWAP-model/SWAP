# Phase 4f Tasks B3–B6: Strangler-Swap parity for grassgrowth, oxygenstress, salinitystress, surfacewater

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Get each of the four remaining non-macropore regression cases to pass `pixi run -e test regression <case>` with the 1e-2 cm tolerance, by authoring the case's TOML files and closing any new schema/adapter gaps surfaced. After all four pass, full `pixi run -e test check-full` is green and the phase is tagged.

**Architecture:** The strangler-swap pattern landed for hupselbrook (commits `6207e83 → abf76bd`). Each remaining case follows the same iteration pattern as hupselbrook B2: (a) audit legacy `.swp` + `.dra` keys against current TOML, (b) close gaps in schema / reader / adapter / case TOML / case `.crp.toml` until regression passes. Crop legacy readers stay in place — `.crp` files are still staged via `run_case.sh --toml`. The per-rotation crop adapter is deferred to Phase 4f-extend.

**Tech Stack:** Fortran (gfortran 13), tomlf, ttutil, meson, pFUnit, pixi, Python (regression). Existing tooling: `pixi run -e test build-linux`, `pixi run -e test regression <case>`, `pixi run -e test check-full`.

---

## File Structure

For every case (B3–B6) the same set of files is touched:

- `tests/swap-cases/toml/<N>.<case>/swap.toml` — the case's top-level TOML, **inside the `tests/swap-cases` submodule**. Modify to author every key the case's `.swp` lists.
- `tests/swap-cases/toml/<N>.<case>/swap.dra.toml` — the case's drainage TOML, **inside the submodule**. Modify to author every key the case's `.dra` lists.
- `src/config/*_config.f90` — extend typed configs only when a case authors a key the schema doesn't model.
- `src/io/toml/read_*_toml.f90` — extend section readers only when adding a schema field.
- `src/io/toml/config_to_variables.f90` — adapter copies. Add HACK markers for narrow scalars; convert HACKs to schema slots only when a Phase 4f-extend pass tackles them.

Reference materials:
- `docs/phase-4f-hupselbrook-parity-audit.md` — the audit document produced for hupselbrook lists every category of gap and the resolution pattern.
- `src/io/readswap.f90` — the legacy reader. Every `rd*` call is a potential schema entry. Each conditional block (e.g. `if (swdra .ne. 0)`, `if (swsolu .eq. 1)`) marks a switch-gated section that may or may not apply to a given case.
- The `adapter` (`src/io/toml/config_to_variables.f90`) header docstring documents the HACK marker convention.
- Case `.swp` / `.dra` / `.crp` source-of-truth files are at `tests/swap-cases/<N>.<case>/swap_linux.swp.template`, `swap.dra`, and `<crop>.crp`.

Each case task ends with a single submodule commit (case TOML changes) + outer-repo commits (schema / reader / adapter changes if any).

---

## Conventions & invariants

These hold for every task below:

1. **Submodule discipline.** Edits to `tests/swap-cases/toml/<case>/*.toml` are committed inside the submodule first, then the outer-repo bump is committed separately. Do not use `git add -A`; name files explicitly.
2. **HACK marker pattern** in `config_to_variables.f90`: a one-line comment `! HACK Phase 4f-extend: <one-line description>` immediately precedes the assignment. Every grep-discoverable HACK must be removable by Phase 4f-extend.
3. **Verification command** for regression is `pixi run -e test regression <case>` — runs the case-only target, ~1 second after a successful build.
4. **Build cycle:** `pixi run -e test build-linux` from repo root. If the build fails, fix the compilation error before re-running regression.
5. **Cropfile staging:** `run_case.sh --toml` (auto-enabled by the auto-detect default) stages `swap.toml + swap.dra.toml + *.crp.toml + swap.swp`. The legacy crop sub-readers in `cropgrowth.f90` consume the staged `.crp` files. Do NOT remove `swap.swp` staging.
6. **No physics edits.** Fixes go in `src/config/*_config.f90`, `src/io/toml/read_*_toml.f90`, `src/io/toml/config_to_variables.f90`, or the case TOML. Do not modify `src/soil/`, `src/atmosphere/`, `src/crop/`, `src/drainage/` (non-config), `src/heat/`, `src/utils/`.
7. **Audit document per case:** for each B3–B6 task, produce `docs/phase-4f-<case>-parity-audit.md` mirroring the hupselbrook audit format. The doc is the gap inventory; the commits close them.
8. **One commit per gap closure.** Don't batch unrelated fixes; mirror the hupselbrook pattern (`shape` fix, `flCropReadFile/Open + rdmax`, `rsigni + cfevappond` were each their own commit). If multiple gaps land in one logical schema extension, that single extension is one commit.
9. **Don't push** anywhere; tagging is local-only at closeout.
10. **Macropore (case 3) is excluded** per ADR 0011. Do not run or include it.

---

## Task 1: Phase 4f Task B3 — grassgrowth parity

**Goal:** `pixi run -e test regression grassgrowth` passes 1e-2 cm tolerance.

**Files:**
- Audit: produce `docs/phase-4f-grassgrowth-parity-audit.md`
- Modify (submodule): `tests/swap-cases/toml/2.grassgrowth/swap.toml`, `tests/swap-cases/toml/2.grassgrowth/swap.dra.toml`
- Modify (outer repo, only if new gaps): `src/config/*_config.f90`, `src/io/toml/read_*_toml.f90`, `src/io/toml/config_to_variables.f90`
- Verification: `tests/swap-cases/2.grassgrowth/swap_linux.swp.template`, `tests/swap-cases/2.grassgrowth/swap.dra`

- [ ] **Step 1: Sanity-check baseline state**

```bash
cd /home/zawadzkim/Code/swap
pixi run -e test regression grassgrowth 2>&1 | tail -25
```

Expected: case fails, possibly with `FATAL ... Could not open file 'swap.toml'` or with a numerical diff table. Record the failure mode — that is the starting state.

- [ ] **Step 2: Author the audit document**

Walk every line of `tests/swap-cases/2.grassgrowth/swap_linux.swp.template` and `tests/swap-cases/2.grassgrowth/swap.dra`. For each `KEY = value` line, classify against current TOML coverage. Save as `docs/phase-4f-grassgrowth-parity-audit.md` with the same section structure as `docs/phase-4f-hupselbrook-parity-audit.md`:

- General + Output
- Meteorology
- Crop
- Irrigation
- Soil water
- Drainage / .dra
- Bottom boundary
- Heat
- Solute

Status legend: OK / TOML / READER / ADAPTER / FINALIZE / N/A.

End with a "Highest-priority gaps" section ordered by expected water-balance impact. Commit:

```bash
git add docs/phase-4f-grassgrowth-parity-audit.md
git commit -m "docs(phase-4f): grassgrowth parity audit (Phase 4f Task B3)"
```

- [ ] **Step 3: Close gaps iteratively**

For each gap in priority order:

- If TOML-missing: author the value in `tests/swap-cases/toml/2.grassgrowth/swap.toml` or `swap.dra.toml`. Submodule commit:
  ```bash
  cd tests/swap-cases
  git add toml/2.grassgrowth/swap.toml  # or swap.dra.toml
  git commit -m "feat(toml/grassgrowth): <one-line description> (Phase 4f Task B3)"
  cd ../..
  git add tests/swap-cases
  git commit -m "chore(swap-cases): bump for <description> (Phase 4f Task B3)"
  ```

- If READER-missing or ADAPTER-missing: extend schema in `src/config/*_config.f90` and reader in `src/io/toml/read_*_toml.f90` if structurally new; otherwise wire in `src/io/toml/config_to_variables.f90` with a HACK marker. Outer-repo commit:
  ```bash
  git add src/config/*_config.f90 src/io/toml/read_*_toml.f90 src/io/toml/config_to_variables.f90
  git commit -m "feat(toml,<section>): add <field> to <schema|adapter> (Phase 4f Task B3)"
  ```

- If FINALIZE-mismatch: port the legacy `readswap.f90` finalize logic to the adapter (mirror `populate_outdatint_monthly` pattern at `src/io/toml/config_to_variables.f90:798`). Outer-repo commit:
  ```bash
  git add src/io/toml/config_to_variables.f90
  git commit -m "fix(adapter): port <finalize> for grassgrowth (Phase 4f Task B3)"
  ```

After each commit, build and run regression:
```bash
pixi run -e test build-linux 2>&1 | grep -E "Error|error:" | head -3
pixi run -e test regression grassgrowth 2>&1 | tail -25
```

If the build fails, fix the compilation error before continuing. If regression still has diffs, return to Step 3 with the next-priority gap.

- [ ] **Step 4: Verify regression passes**

```bash
pixi run -e test regression grassgrowth
```

Expected output:
```
✓ grassgrowth: regression ok (annual stats match fixture)
Results: 1 passed, 0 failed
```

If still failing after the audit gaps are closed, augment the audit doc with the residual diff and dispatch a focused diagnostic agent (mirroring the hupselbrook approach used for `shape`, `rsigni`, `swmonth` etc.).

---

## Task 2: Phase 4f Task B4 — oxygenstress parity

**Goal:** `pixi run -e test regression oxygenstress` passes.

**Files:** identical structure to Task 1, with paths under `4.oxygenstress`.

Oxygenstress is the case where `SwOxygen` (oxygen-stress switch) flips into a non-default mode. Bulk density (`bdens`) becomes load-bearing because the legacy code at `src/crop/cropgrowth.f90:3668` aborts if BDENS values aren't realistic when `SwOxygen=2`. Confirm the .swp's value flows into the global early in the audit.

- [ ] **Step 1: Sanity-check baseline state**

```bash
pixi run -e test regression oxygenstress 2>&1 | tail -25
```

- [ ] **Step 2: Author the audit document**

Walk `tests/swap-cases/4.oxygenstress/swap_linux.swp.template` and `swap.dra`. Save as `docs/phase-4f-oxygenstress-parity-audit.md`. **Pay extra attention to `SwOxygen`, `BDENS`, and any `SOX_*` (oxygen-related) keys.**

Commit:
```bash
git add docs/phase-4f-oxygenstress-parity-audit.md
git commit -m "docs(phase-4f): oxygenstress parity audit (Phase 4f Task B4)"
```

- [ ] **Step 3: Close gaps iteratively**

Same iteration pattern as Task 1 Step 3, with `oxygenstress` substituted for `grassgrowth` and case path `4.oxygenstress`. Commit message format:
```
feat(toml/oxygenstress): <description> (Phase 4f Task B4)
chore(swap-cases): bump for <description> (Phase 4f Task B4)
feat(toml,<section>): <field> (Phase 4f Task B4)
fix(adapter): <fix> for oxygenstress (Phase 4f Task B4)
```

- [ ] **Step 4: Verify regression passes**

```bash
pixi run -e test regression oxygenstress
```

Expected: `✓ oxygenstress: regression ok`.

---

## Task 3: Phase 4f Task B5 — salinitystress parity

**Goal:** `pixi run -e test regression salinitystress` passes.

**Files:** identical structure, paths under `5.salinitystress`.

Salinitystress is the case where the `[solute]` block does real work. `swsolu`, `cdrain`, `cseep`, `tscf`, `ldis`, the per-layer `cml` initial concentrations, and possibly `swsp`/`frexp`/`cref`/`kf` may all be load-bearing. Drainage may use a different DRAMET than hupselbrook.

- [ ] **Step 1: Sanity-check baseline state**

```bash
pixi run -e test regression salinitystress 2>&1 | tail -25
```

- [ ] **Step 2: Author the audit document**

Walk `tests/swap-cases/5.salinitystress/swap_linux.swp.template` and `swap.dra`. Save as `docs/phase-4f-salinitystress-parity-audit.md`. **Pay extra attention to the [solute] section: every legacy `rd*` call inside `if (swsolu .eq. 1)` blocks in `src/io/readswap.f90`.**

Commit:
```bash
git add docs/phase-4f-salinitystress-parity-audit.md
git commit -m "docs(phase-4f): salinitystress parity audit (Phase 4f Task B5)"
```

- [ ] **Step 3: Close gaps iteratively**

Same iteration pattern. Commit format:
```
feat(toml/salinitystress): <description> (Phase 4f Task B5)
chore(swap-cases): bump for <description> (Phase 4f Task B5)
feat(toml,solute): <field> (Phase 4f Task B5)
fix(adapter): <fix> for salinitystress (Phase 4f Task B5)
```

- [ ] **Step 4: Verify regression passes**

```bash
pixi run -e test regression salinitystress
```

Expected: `✓ salinitystress: regression ok`.

---

## Task 4: Phase 4f Task B6 — surfacewater parity

**Goal:** `pixi run -e test regression surfacewater` passes.

**Files:** identical structure, paths under `6.surfacewater`.

Surfacewater is the case where the `swsrf` / `swsec` / `nmper` / `wscap` / `wldip` / `intwl` block is fully populated — this is the largest bock of cross-table inputs. `swqhr` may switch the discharge formula; `OWLTAB` time-series tables surface here. `nmper`, the per-period management arrays (`swman`, `wscap`, `wldip`, `intwl`, `impend`) and the `qhrtab`/`qtab` tables all matter.

- [ ] **Step 1: Sanity-check baseline state**

```bash
pixi run -e test regression surfacewater 2>&1 | tail -25
```

- [ ] **Step 2: Author the audit document**

Walk `tests/swap-cases/6.surfacewater/swap_linux.swp.template` and `swap.dra`. Save as `docs/phase-4f-surfacewater-parity-audit.md`. **Pay extra attention to the surface-water management blocks (Part 4A–4D in the legacy `.dra`) and the `OWLTAB` per-level water-level time series.**

Commit:
```bash
git add docs/phase-4f-surfacewater-parity-audit.md
git commit -m "docs(phase-4f): surfacewater parity audit (Phase 4f Task B6)"
```

- [ ] **Step 3: Close gaps iteratively**

Same iteration pattern. Commit format:
```
feat(toml/surfacewater): <description> (Phase 4f Task B6)
chore(swap-cases): bump for <description> (Phase 4f Task B6)
feat(toml,surface_water): <field> (Phase 4f Task B6)
fix(adapter): <fix> for surfacewater (Phase 4f Task B6)
```

- [ ] **Step 4: Verify regression passes**

```bash
pixi run -e test regression surfacewater
```

Expected: `✓ surfacewater: regression ok`.

---

## Task 5: Closeout — full check-full + tag

**Goal:** All 5 non-macropore cases green, phase tagged, main fast-forwarded.

**Files:**
- Modify: `pixi.toml` (no changes expected, but verify check-fast / check-full task targets)
- Modify (submodule, only if changes are present): pointer bump
- Verification: full check-full

- [ ] **Step 1: Run full check-full**

```bash
cd /home/zawadzkim/Code/swap
pixi run -e test check-full 2>&1 | tail -30
```

Expected output:
```
✓ hupselbrook: regression ok
✓ grassgrowth: regression ok
✓ oxygenstress: regression ok
✓ salinitystress: regression ok
✓ surfacewater: regression ok
Results: 5 passed, 0 failed
```

(`macroporeflow` is excluded per ADR 0011 — `tests/regression/test_output_regression.py` already drops it.)

If any case fails: return to the relevant Task N (1–4) Step 4 to close the residual gap. Each case is independent; a residual diff in case 5 doesn't require revisiting case 2.

- [ ] **Step 2: Run pFUnit suite to confirm unit tests still green**

```bash
pixi run -e test test-pfunit 2>&1 | tail -10
```

Expected:
```
1/1 swap:unit-pfunit / unit-swap-tests OK
Ok: 1
```

- [ ] **Step 3: Tag the phase**

```bash
git tag rescue/phase-4f-strangler-swp
git tag --list rescue/*
```

Expected: the new tag appears in the list. Tag is local — do not push.

- [ ] **Step 4: Fast-forward main to development**

Confirm with the user before this step — main is the long-lived branch. If the user approves:

```bash
git checkout main
git merge --ff-only development
git checkout development
```

Expected: clean fast-forward, no merge commit. If non-fast-forward, abort and investigate (the user may have committed to main directly).

- [ ] **Step 5: Update phase status doc**

Locate the phase status doc (e.g. `docs/PHASE_4F_STATUS.md` or similar — search with `ls docs/PHASE_*` first). If one exists, update it; if not, append a one-liner to `docs/superpowers/specs/2026-04-29-phase-4f-strangler-swp-design.md`:

```markdown
## Closeout

Phase 4f-strangler-swp tagged at <commit-sha> on <date>.
All 5 non-macropore regression cases pass via the new TOML pipeline.

Outstanding work for Phase 4f-extend:
- Replace adapter HACK markers with proper schema slots
  (grep `HACK Phase 4f-extend` in src/io/toml/config_to_variables.f90)
- Build the per-rotation crop adapter (cropfixed → wofost → grass)
- Drop `swap.swp` staging from `run_case.sh`
- Macropore case re-port (case 3, ADR 0011)
```

Commit:
```bash
git add docs/superpowers/specs/2026-04-29-phase-4f-strangler-swp-design.md
git commit -m "docs(phase-4f): closeout — all non-macropore cases green"
```

---

## Risks & mitigations

| Risk | Mitigation |
|---|---|
| A case's `.swp` exposes a fundamentally new section (e.g. SWSCAL, FILENAMESOPHY) | Schema-extend, document the gap in the case audit, HACK marker if narrow. Don't physics-port — the legacy reader is the parity target. |
| Schema extension for case N breaks earlier case's regression | Re-run `pixi run -e test regression <earlier-case>` after each schema extension; rollback the extension if numbers move. The validators are switch-gated, so additive fields shouldn't impact unauthored cases — but verify. |
| Auditor misses a finalize transformation in `readswap.f90` | The audit document is required precisely to surface these. If a numerical diff persists after every authored key has been wired, search `readswap.f90` for non-`rd*` assignments inside the relevant section block (e.g. `paramvg(7,lay) = 1 - 1/npar`, `l(1) = 100*lm2`). |
| Submodule commit / outer-repo bump get out of sync | Always `cd tests/swap-cases && git commit ... && cd ../.. && git add tests/swap-cases && git commit ...` as a pair. Verify with `git submodule status` after each pair. |
| HACK markers proliferate without follow-up | The phase docstring at the top of `config_to_variables.f90` documents the convention; closeout Step 5 lists them as Phase 4f-extend work. Acceptable for now. |
| User wants a different commit cadence (e.g. squash per case) | Confirm at task start. Default = one commit per gap closure to mirror hupselbrook B2. |

---

## Verification

After each Task N (1–4):
```bash
pixi run -e test regression <case>
```
Must show `✓ <case>: regression ok`.

Final closeout (Task 5):
```bash
pixi run -e test check-full      # 5/5 green
pixi run -e test test-pfunit     # OK
git tag --list rescue/*          # rescue/phase-4f-strangler-swp present
```
