# Rescue Phase 0 — Baseline Reset Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Establish an unambiguous git topology anchored on the last known-green commit (`e256bc0`), preserve all potentially-salvageable drifted work as archive branches, capture the original SWAP 4.2.0 reference implementation as an orphan branch, and record the verified baseline timings and tolerances.

**Architecture:** Local-only git operations (no pushes during this phase). Reset `main` to `e256bc0`, create `development` from it, archive the drifted branches by renaming under `archive/` prefix, create `legacy/swap-4.2.0` as a curated orphan branch from the `swap_org/` tree. Verify the full regression suite green at documented tolerances and tag `rescue/phase-0-baseline`.

**Tech Stack:** git, pixi, meson, gfortran, pFUnit, pytest (via `pixi run regression`).

**Spec:** [docs/superpowers/specs/2026-04-22-rescue-and-stabilize-design.md](../specs/2026-04-22-rescue-and-stabilize-design.md) — Phase 0

**Starting state at the time of execution:**
- `HEAD` detached at `e256bc0` (the green baseline).
- `main` at `7587ca3` locally and on `origin/main` (drift-era tip).
- Branches present: `main`, `swaplib`, `swaplib-simple`, plus the detached HEAD.
- The session-start drifted commit `cfb4aee` is dangling (no branch points to it).
- Working tree has drift-era modifications to tracked files (from the pre-reset state) and untracked directories (`fpm_install/`, `builddir/`, spec file written during brainstorming).
- `swap_org/` present on disk but untracked; 36 MB.

**End state after this plan:**
- Local `main` = `e256bc0`. `origin/main` **unchanged** (still `7587ca3`).
- `development` branched from `main` at `e256bc0`, plus one commit adding the rescue spec.
- Branches `archive/wip-drifted`, `archive/main-pre-rescue`, `archive/swaplib`, `archive/swaplib-simple` exist locally.
- Branches `swaplib` and `swaplib-simple` deleted locally (archives preserve their content).
- Orphan branch `legacy/swap-4.2.0` contains a curated snapshot of `swap_org/`.
- Tag `rescue/phase-0-baseline` points at the spec commit on `development`.
- `docs/superpowers/specs/2026-04-22-baseline-record.md` records exact timings and accepted tolerances observed on this machine.

---

## Task 1: Pre-flight reconnaissance

**Files:** none created or modified.

- [ ] **Step 1: Confirm we are detached at `e256bc0`**

Run:

```bash
git rev-parse HEAD
git status --short | head -30
git branch -a
```

Expected: `HEAD` resolves to `e256bc0a43e349ad3692f6c97fa98d06e07943d8`. `git branch -a` shows `* (HEAD detached at e256bc0)`, `main`, `swaplib`, `swaplib-simple`, `remotes/origin/main`.

If `HEAD` is not at `e256bc0`, stop and investigate before proceeding. All subsequent tasks assume this baseline.

- [ ] **Step 2: Confirm the dangling drifted commit is still reachable**

Run:

```bash
git log --oneline cfb4aee -1
```

Expected: outputs `cfb4aee refactor: begin to work on the output writing.` (or similar title). If this fails with "bad object," the commit has already been garbage-collected and we cannot archive it — proceed anyway (nothing we can do), and note in the baseline record that `wip-drifted` could not be recovered.

- [ ] **Step 3: Record the drifted main tip**

Run:

```bash
git log --oneline 7587ca3 -1
git log --oneline main -1
```

Expected: `main` points at `7587ca3 include the new updated workflow excluding the eindows build`. Record this SHA — it is what `archive/main-pre-rescue` will point to.

- [ ] **Step 4: Confirm GPL v2 LICENSE at repo root**

Run:

```bash
head -3 LICENSE
wc -l LICENSE
```

Expected: first few lines name the GPL or reference it, length ~339 lines (standard GPL v2 text). If the file is missing or is not GPL v2, note this — Phase 2 will add/correct it, but the legacy branch creation in Task 10 depends on GPL v2 being the licensing basis.

---

## Task 2: Create safety archive branches

**Files:** none modified. Pure git ref manipulation.

- [ ] **Step 1: Archive the dangling drifted HEAD**

Run:

```bash
git branch archive/wip-drifted cfb4aee
```

Expected: silent success. If this fails because `cfb4aee` is unreachable (garbage collected), skip with a note; continue to the next step.

- [ ] **Step 2: Archive the drift-era main tip**

Run:

```bash
git branch archive/main-pre-rescue 7587ca3
```

Expected: silent success.

- [ ] **Step 3: Archive the BMI / Python-binding experiment branch**

Run:

```bash
git branch archive/swaplib swaplib
```

Expected: silent success.

- [ ] **Step 4: Archive the state-sync experiment branch**

Run:

```bash
git branch archive/swaplib-simple swaplib-simple
```

Expected: silent success.

- [ ] **Step 5: Verify all four archives exist and point where expected**

Run:

```bash
for b in archive/wip-drifted archive/main-pre-rescue archive/swaplib archive/swaplib-simple; do
  echo "=== $b ==="
  git log --oneline -1 "$b"
done
```

Expected: four lines of output, one per archive branch, each showing the commit it was branched from.

- [ ] **Step 6: Commit — nothing to commit in this task** (git branch operations don't require commits). Skip.

---

## Task 3: Reset local `main` to baseline and create `development`

**Files:** working tree returns to exact `e256bc0` state for tracked files. Untracked files (including the rescue spec) preserved.

- [ ] **Step 1: Reset `main` pointer to `e256bc0` (local only)**

Run:

```bash
git branch -f main e256bc0
```

Expected: silent success. This moves the local `main` pointer from `7587ca3` to `e256bc0` without touching the working tree or `origin/main`.

- [ ] **Step 2: Verify `main` now points at the baseline**

Run:

```bash
git log --oneline -1 main
git log --oneline -1 origin/main
```

Expected: `main` shows `e256bc0 refactor: cleaning up tests and docs`. `origin/main` still shows `7587ca3 include the new updated workflow...` (unchanged). If `origin/main` also changed, stop — something pushed accidentally.

- [ ] **Step 3: Create `development` branch at the baseline**

Run:

```bash
git branch development e256bc0
```

Expected: silent success.

- [ ] **Step 4: Discard modifications to tracked files so a branch checkout can succeed**

Run:

```bash
git checkout -- .
git status --short | grep -v "^??" | head
```

Expected: second command prints nothing (no modified tracked files remain). The untracked (`??`) entries — including the rescue spec at `docs/superpowers/specs/2026-04-22-rescue-and-stabilize-design.md` — are preserved.

- [ ] **Step 5: Check out `development`**

Run:

```bash
git checkout development
git rev-parse HEAD
git branch --show-current
```

Expected: `HEAD` resolves to `e256bc0...`. Current branch is `development`. Untracked files still present.

- [ ] **Step 6: Commit — nothing to commit** (branch pointer manipulation only). Skip.

---

## Task 4: Remove drift-era untracked clutter

**Files:** delete `fpm_install/`, `build/`, `builddir/` at repo root (all rebuildable; none in git at this baseline).

- [ ] **Step 1: Confirm each target is NOT tracked (safety check)**

Run:

```bash
for d in fpm_install build builddir; do
  if git ls-files --error-unmatch "$d" 2>/dev/null | head -1; then
    echo "DANGER: $d is tracked — do not delete"
  else
    echo "OK: $d not tracked, safe to delete"
  fi
done
```

Expected: three "OK" lines. If any says DANGER, stop and investigate.

- [ ] **Step 2: Delete the three directories**

Run:

```bash
rm -rf fpm_install build builddir
ls -d fpm_install build builddir 2>&1 | head -3
```

Expected: second command says each "cannot access ... No such file or directory".

- [ ] **Step 3: Commit — nothing to commit yet** (these were untracked; deletion doesn't produce a commit). Skip.

---

## Task 5: Commit the rescue spec to `development`

**Files:**
- Add: `docs/superpowers/specs/2026-04-22-rescue-and-stabilize-design.md` (already on disk).
- Directory `docs/superpowers/plans/` also already exists; this plan file will be included.

- [ ] **Step 1: Confirm the spec file is present**

Run:

```bash
ls -la docs/superpowers/specs/2026-04-22-rescue-and-stabilize-design.md
wc -l docs/superpowers/specs/2026-04-22-rescue-and-stabilize-design.md
```

Expected: file present, ~200–260 lines.

- [ ] **Step 2: Confirm this plan file is present**

Run:

```bash
ls -la docs/superpowers/plans/2026-04-22-rescue-phase-0-baseline.md
```

Expected: file present.

- [ ] **Step 3: Stage spec and plan**

Run:

```bash
git add docs/superpowers/specs/2026-04-22-rescue-and-stabilize-design.md
git add docs/superpowers/plans/2026-04-22-rescue-phase-0-baseline.md
git status --short
```

Expected: two `A` entries for the two files, nothing else.

- [ ] **Step 4: Commit**

Run:

```bash
git commit -m "docs: add rescue & stabilize spec and phase 0 plan

Spec captures the Rescue & Stabilize workflow that returns the
modernization branch to a clean, tested, documented baseline from
which compartment-state refactoring, Python bindings, and GPU work
proceed as follow-on specs. Plan implements Phase 0: baseline reset
and archival of drifted branches."
```

Expected: commit succeeds. Note the new commit SHA.

- [ ] **Step 5: Verify the commit**

Run:

```bash
git log --oneline -2
```

Expected: newest commit is the one just made; its parent is `e256bc0`.

---

## Task 6: Baseline verification — compile

**Files:** none modified. `pixi run build-linux` will create a `builddir/`.

- [ ] **Step 1: Build**

Run:

```bash
pixi run build-linux 2>&1 | tail -30
```

Expected: ends with a successful meson/ninja completion. An executable `swap` is produced in `builddir/`.

- [ ] **Step 2: Verify the binary exists**

Run:

```bash
ls -la builddir/swap
file builddir/swap
```

Expected: file exists, identified as an ELF executable.

- [ ] **Step 3: Commit — nothing to commit** (build artifacts are gitignored). Skip.

---

## Task 7: Baseline verification — unit tests

**Files:** none modified.

- [ ] **Step 1: Identify the correct pFUnit task name**

Run:

```bash
grep -nE "pfunit|test-pfunit|_compile-pfunit" pixi.toml | head -15
```

Expected: lines referencing `_configure-pfunit`, `_compile-pfunit`, `test-pfunit`. Use whichever task name runs the pFUnit suites directly. At baseline this is `test-pfunit`.

- [ ] **Step 2: Run pFUnit suites**

Run:

```bash
pixi run test-pfunit 2>&1 | tail -40
```

Expected: all suites green. If any suite fails at this baseline, the failure was already present at `e256bc0` and must be documented as a pre-existing tolerance (same category as the MOWDM oxygenstress deviation) — do not skip investigation.

- [ ] **Step 3: Record the unit-test result**

Note for Task 9: "pFUnit at e256bc0: N suites green" (or document specific pre-existing failures).

---

## Task 8: Baseline verification — regression (full)

**Files:** none modified. Results accumulate under `tests/swap-cases/*/result*` (gitignored).

- [ ] **Step 1: Initialize submodule if needed**

Run:

```bash
git submodule status tests/swap-cases
```

If the status line starts with `-` (not initialized) run:

```bash
git submodule update --init tests/swap-cases
```

Expected: submodule ready.

- [ ] **Step 2: Run the full regression suite**

Run:

```bash
time pixi run regression 2>&1 | tee /tmp/phase0-regression.log | tail -30
```

Expected, based on the previously observed baseline run:

- `surfacewater`: ok, ~3s
- `hupselbrook`: ok, ~10s
- `salinitystress`: ok, ~27s
- `grassgrowth`: ok, ~27s
- `oxygenstress`: **FAIL on MOWDM** with deviations of 0.1 to 85.0 (pre-existing, accepted tolerance)
- `macropore`: ok, ~344s
- Total wall time ~590s; "Results: 5 passed, 1 failed" where the failure is the MOWDM deviation only.

If results diverge beyond this pattern, stop and investigate — we are not at the expected baseline.

- [ ] **Step 3: Preserve the log**

Run:

```bash
mkdir -p tests/regression/baselines
cp /tmp/phase0-regression.log tests/regression/baselines/phase-0-baseline.log
```

Expected: file copied.

- [ ] **Step 4: Commit — defer** to Task 9, which commits the full baseline record together with this log.

---

## Task 9: Record baseline timings and tolerances

**Files:**
- Create: `docs/superpowers/specs/2026-04-22-baseline-record.md`
- Keep: `tests/regression/baselines/phase-0-baseline.log` (from Task 8 Step 3)

- [ ] **Step 1: Write the baseline record**

Create `docs/superpowers/specs/2026-04-22-baseline-record.md` with this exact content (replace placeholder timings with actual observed values from Task 8 if they differ):

````markdown
---
title: Baseline record at rescue/phase-0-baseline
date: 2026-04-22
commit: e256bc0a43e349ad3692f6c97fa98d06e07943d8
status: accepted
---

# Phase 0 Baseline Record

This file locks in the observed state of the `e256bc0` baseline at the moment the Rescue & Stabilize workflow started. Future phases verify against these numbers and flag regressions.

## Git topology

- `main` = `development` = `e256bc0`
- `origin/main` untouched (still `7587ca3`)
- Archives: `archive/wip-drifted` (`cfb4aee`), `archive/main-pre-rescue` (`7587ca3`), `archive/swaplib`, `archive/swaplib-simple`
- Orphan: `legacy/swap-4.2.0` created in Task 10

## Compile

- `pixi run build-linux` succeeds with gfortran.
- Produces `builddir/swap`.

## Unit tests (pFUnit)

- `pixi run test-pfunit`: all suites green at `e256bc0`.
- Any failures observed at this commit are recorded below; none expected.

## Regression (full)

Full suite wall time: ~590 s. Run log preserved at `tests/regression/baselines/phase-0-baseline.log`.

| Case           | Status | Time (s) | Notes                                            |
|----------------|--------|---------:|--------------------------------------------------|
| surfacewater   | pass   | ~3       |                                                  |
| hupselbrook    | pass   | ~10      |                                                  |
| salinitystress | pass   | ~27      |                                                  |
| grassgrowth    | pass   | ~27      |                                                  |
| oxygenstress   | fail   | ~178     | MOWDM deviations (pre-existing, accepted)        |
| macropore      | pass   | ~344     | Accepted; perf-regression follow-on spec tracks it |

## Accepted tolerances

1. **`oxygenstress` MOWDM deviations** — fixture expected vs actual differ by up to 85.0 in year 1995 (specific values in the run log). These deviations existed before the rescue and are accepted; do not treat as a regression during Phase 0–4.
2. **`macropore` runtime ~344 s** — approximately double the pre-drift ~170 s baseline. Accepted for this spec; addressed in a separate performance follow-on spec between Phase 4 and the compartment-state refactor.

## Fast test set

Phases 1–4 use `pixi run check-fast` (added in Phase 1) for everyday iteration:

- surfacewater, hupselbrook, salinitystress, grassgrowth (~67 s total).

Full set is re-run before every phase tag via `pixi run check-full`.
````

- [ ] **Step 2: Stage and commit**

Run:

```bash
git add docs/superpowers/specs/2026-04-22-baseline-record.md tests/regression/baselines/phase-0-baseline.log
git commit -m "docs: record phase 0 baseline timings and accepted tolerances

Locks in the observed e256bc0 state: build, unit tests, and full
regression suite timings. Documents the two known deviations
(oxygenstress MOWDM and macropore runtime) as accepted tolerances
for the duration of the rescue spec."
```

Expected: commit succeeds.

---

## Task 10: Create `legacy/swap-4.2.0` orphan branch with curated `swap_org/` content

**Files:** operates primarily on `swap_org/` via a staging directory. On `development` the working tree is unchanged at the end.

The orphan branch is created through a temporary git worktree so we never mutate `development`'s working tree or risk checking `swap_org/` files into it.

- [ ] **Step 1: Stage curated content outside the repo**

Run:

```bash
rm -rf /tmp/swap-legacy-staging
mkdir /tmp/swap-legacy-staging
cd swap_org && tar cf - \
  --exclude='bin' \
  --exclude='linux/swap420' \
  --exclude='linux/*.o' \
  --exclude='linux/*.mod' \
  --exclude='*.zip' \
  --exclude='rsoftware' \
  . | tar xf - -C /tmp/swap-legacy-staging
cd ..
ls /tmp/swap-legacy-staging/
du -sh /tmp/swap-legacy-staging
```

Expected: staging directory contains `source_swp_4.2.0/`, `source_ttutil_4.27/`, `cases/`, `compiler_settings/`, `doc/`, `license/`, `linux/` (whatever is left after exclusions — likely just build scripts), `readme_4.2.0.txt`, `xdata/`. No `bin/`, no `*.zip`, no compiled binary, no `rsoftware/`. Size substantially smaller than the 36 MB original.

- [ ] **Step 2: Review what was kept and what was excluded**

Run:

```bash
find /tmp/swap-legacy-staging -name '*.zip' -o -name 'swap420' -o -name '*.o' -o -name '*.mod' | head
du -sh /tmp/swap-legacy-staging/*
```

Expected: first command prints nothing. Second shows per-directory sizes; no single directory is disproportionately large.

- [ ] **Step 3: Create orphan branch via a temporary worktree**

Run:

```bash
git worktree add --detach /tmp/swap-legacy-wt
cd /tmp/swap-legacy-wt
git checkout --orphan legacy/swap-4.2.0
git rm -rf --cached . >/dev/null 2>&1 || true
find . -mindepth 1 -maxdepth 1 ! -name '.git' -exec rm -rf {} +
cp -a /tmp/swap-legacy-staging/. .
ls -la | head
```

Expected: `ls` shows the curated content at the worktree root; no leftover files from `development` besides `.git`.

- [ ] **Step 4: Commit the orphan root**

Run:

```bash
cd /tmp/swap-legacy-wt
git add .
git commit -m "legacy: import SWAP 4.2.0 reference implementation (GPL v2)

Curated snapshot of swap_org/ captured as an orphan branch so the
reference implementation lives in the repository rather than only on
the author's machine. Excluded: compiled binaries, .zip archives,
rsoftware/, build artifacts. Preserved: source, ttutil source, test
case inputs, docs, compiler settings, xdata, and the LICENSE files
required by GPL v2 §1."
git log --oneline -1
cd /home/zawadzkim/Code/swap
```

Expected: one commit on `legacy/swap-4.2.0` with no parent. `git log --oneline -1` shows exactly one line.

- [ ] **Step 5: Remove the temporary worktree**

Run:

```bash
git worktree remove /tmp/swap-legacy-wt
rm -rf /tmp/swap-legacy-staging
git branch -a | grep legacy
```

Expected: worktree removed cleanly; `legacy/swap-4.2.0` still present in branch listing.

- [ ] **Step 6: Verify the orphan branch has exactly one commit with no parent**

Run:

```bash
git log --oneline legacy/swap-4.2.0
git rev-list --parents legacy/swap-4.2.0
```

Expected: single commit. `git rev-list --parents` shows the commit SHA followed by nothing (no parent SHA).

- [ ] **Step 7: Confirm we're back on `development`**

Run:

```bash
git branch --show-current
git status --short | head -10
```

Expected: on `development`, clean working tree (except `swap_org/` which remains untracked on this branch and should be added to `.gitignore` in Phase 1).

- [ ] **Step 8: Commit — no commit needed on `development`** (the orphan branch is self-contained). Skip.

---

## Task 11: Delete superseded local branches

**Files:** none.

- [ ] **Step 1: Confirm archives cover the branches about to be deleted**

Run:

```bash
git log --oneline -1 archive/swaplib
git log --oneline -1 archive/swaplib-simple
diff <(git log --oneline swaplib) <(git log --oneline archive/swaplib) && echo "swaplib archive matches"
diff <(git log --oneline swaplib-simple) <(git log --oneline archive/swaplib-simple) && echo "swaplib-simple archive matches"
```

Expected: both "matches" messages print. If either diff is non-empty, stop — the archive is not identical to the branch and deletion would lose history.

- [ ] **Step 2: Delete `swaplib`**

Run:

```bash
git branch -D swaplib
```

Expected: "Deleted branch swaplib (was <sha>)."

- [ ] **Step 3: Delete `swaplib-simple`**

Run:

```bash
git branch -D swaplib-simple
```

Expected: "Deleted branch swaplib-simple (was <sha>)."

- [ ] **Step 4: Verify final branch set**

Run:

```bash
git branch
git branch --list 'archive/*'
git branch --list 'legacy/*'
```

Expected:
- Local branches: `main`, `development` (current).
- `archive/*`: four entries (`archive/main-pre-rescue`, `archive/swaplib`, `archive/swaplib-simple`, `archive/wip-drifted` — the last only if Task 2 Step 1 succeeded).
- `legacy/*`: one entry (`legacy/swap-4.2.0`).

- [ ] **Step 5: Commit — nothing to commit** (branch deletion only). Skip.

---

## Task 12: Tag `rescue/phase-0-baseline`

**Files:** none.

- [ ] **Step 1: Confirm `development` HEAD is the spec + baseline-record commits**

Run:

```bash
git log --oneline development | head -3
```

Expected: top entry is the "docs: record phase 0 baseline timings" commit, parent is the "docs: add rescue & stabilize spec" commit, grandparent is `e256bc0`.

- [ ] **Step 2: Create the tag on `development` HEAD**

Run:

```bash
git tag -a rescue/phase-0-baseline -m "Phase 0 complete: baseline reset at e256bc0.

main and development anchored at e256bc0; archive branches for
wip-drifted, pre-rescue main, swaplib, and swaplib-simple preserved;
legacy/swap-4.2.0 orphan branch captures SWAP 4.2.0 reference tree
under GPL v2; baseline timings recorded in
docs/superpowers/specs/2026-04-22-baseline-record.md."
```

Expected: silent success.

- [ ] **Step 3: Verify the tag**

Run:

```bash
git tag -l 'rescue/*'
git show --stat rescue/phase-0-baseline | head -20
```

Expected: tag exists and points at the expected commit.

- [ ] **Step 4: Commit — tags are not commits. Skip.**

---

## Task 13 (OPTIONAL — requires explicit user confirmation): Push archive, legacy, and tag to `origin`

Pushes are visible to others (`origin` is the public SWAP-model/SWAP repository). This task must not run without the user's explicit go-ahead at execution time.

**Files:** none. Network operations only.

- [ ] **Step 1: Ask the user**

Prompt the user:

> "Phase 0 is complete locally. Push `archive/*`, `legacy/swap-4.2.0`, and the `rescue/phase-0-baseline` tag to `origin` (github.com/SWAP-model/SWAP)? Local `main` and `development` will NOT be pushed. Y/N?"

If the answer is no or unclear, skip this task. Phase 0 is still complete.

- [ ] **Step 2: Push archives and legacy**

Only if the user answered yes:

```bash
git push origin \
  archive/wip-drifted \
  archive/main-pre-rescue \
  archive/swaplib \
  archive/swaplib-simple \
  legacy/swap-4.2.0 \
  rescue/phase-0-baseline
```

Expected: all refs pushed successfully. If the push fails because the remote rejects, stop and report.

- [ ] **Step 3: Verify the remote has the new refs**

Run:

```bash
git ls-remote origin | grep -E 'archive/|legacy/|rescue/' | head
```

Expected: the newly-pushed refs visible on `origin`.

- [ ] **Step 4: Commit — nothing to commit.** Skip.

---

## Plan self-review

**Spec coverage.** Phase 0 of the spec lists 10 numbered steps (inventory, wip-drifted archive, main-pre-rescue archive, main reset, development creation, legacy/swap-4.2.0 orphan, delete swaplib/swaplib-simple, verify, commit spec, tag). Mapping:

- Spec step 1 (inventory) → Task 1 + Task 2 Step 1 (also archives anything unique).
- Spec step 2 (`archive/wip-drifted`) → Task 2 Step 1.
- Spec step 3 (`archive/main-pre-rescue`) → Task 2 Step 2.
- Spec step 4 (reset local `main`, no push) → Task 3 Steps 1–2.
- Spec step 5 (create `development`, commits go there) → Task 3 Steps 3–5.
- Spec step 6 (`legacy/swap-4.2.0` with curation and license confirmation) → Task 10 (curation executed; license re-verified in Task 1 Step 4).
- Spec step 7 (delete `swaplib` / `swaplib-simple`) → Task 11.
- Spec step 8 (verify `check-fast` + `check-full`) → Tasks 6–8. `check-fast` and `check-full` pixi tasks don't exist yet; they are introduced in Phase 1. Phase 0 uses `build-linux` + `test-pfunit` + `regression` (full) instead, which covers the same ground. Deviation noted here so the Phase 1 plan does not assume they're already present.
- Spec step 9 (commit spec to `development`) → Task 5.
- Spec step 10 (tag `rescue/phase-0-baseline`) → Task 12.

Additional items this plan adds that the spec implies but doesn't enumerate: pre-flight verification (Task 1), cleaning drift-era untracked dirs (Task 4), recording baseline timings as a committed artifact (Task 9), and the optional push step (Task 13).

**Placeholder scan.** No "TBD", "TODO", or "fill in later" instructions. Every command is concrete.

**Type / name consistency.** Branch names match the spec throughout: `main`, `development`, `archive/wip-drifted`, `archive/main-pre-rescue`, `archive/swaplib`, `archive/swaplib-simple`, `legacy/swap-4.2.0`. Tag name matches: `rescue/phase-0-baseline`. Spec file path matches: `docs/superpowers/specs/2026-04-22-rescue-and-stabilize-design.md`.

**Known risks at execution time:**

1. `cfb4aee` may be garbage-collected before Task 2 runs. Task 1 Step 2 checks for this; if unreachable, `archive/wip-drifted` cannot be created and the step is skipped with a note. The user authorised treating the drifted work as expendable.
2. The regression suite takes ~10 minutes; running it is the longest single step in this plan.
3. `origin/main` must not change during execution. Task 3 Step 2 verifies this explicitly.
4. Curated `swap_org/` size depends on what's inside `cases/` and `doc/`; if unexpectedly large (>50 MB), revisit the curation list before committing the orphan branch.
5. If the `docs/superpowers/specs/2026-04-22-baseline-record.md` placeholder timings differ from observed values in Task 8, the executor edits them before committing in Task 9 Step 1. The file is a record, not a contract.
