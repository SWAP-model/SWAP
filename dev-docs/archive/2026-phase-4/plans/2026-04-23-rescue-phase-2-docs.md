# Rescue Phase 2 — Canonical Documentation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the drift-era document scatter under `docs/` with a small coherent set of canonical documents, wire up `fprettify` so the linter actually enforces what the style guide says, add per-subdir READMEs, fix the root LICENSE per the Phase 0 finding, and refresh the FORD project file so `pixi run -e docs docs-build` produces a clean site.

**Architecture:** Purge drift-era pages wholesale (git history on `archive/wip-drifted` preserves them). Write eight top-level canonical docs and three new ADRs. Configure `.fprettify.rc` to match `code-style.md` so documentation and linter stay in lockstep. Per-subdir READMEs for each `src/<domain>/`. Replace the root LICENSE (currently LGPL v2.1) with GPL v2 per ADR 0005. Update `docs.md` (FORD project file) so the new layout renders cleanly. End with `check-full` green and tag `rescue/phase-2-docs`.

**Tech Stack:** Markdown, FORD 7.0, fprettify 0.3.x, pixi, git.

**Spec:** [docs/superpowers/specs/2026-04-22-rescue-and-stabilize-design.md](../specs/2026-04-22-rescue-and-stabilize-design.md) — Phase 2

**Starting state at execution:**
- `development` HEAD = `17a49f2` (Phase 1 closeout). `main` = `17a49f2`. `origin/main` = `7587ca3` (untouched).
- Tags: `rescue/phase-0-baseline` (`ad4411d`), `rescue/phase-1-infra` (`17a49f2`).
- `docs/adr/` has `0001-gfortran-first.md`, `0002-single-builddir.md` (Phase 1).
- `docs/superpowers/` has specs and plans (preserved).
- `docs/pfunit-vendoring.md` present (absorbed into `build-and-test.md` in Task 6, then deleted).
- `docs/code-style-guide.md` present (superseded by `code-style.md` in Task 8).
- `docs/` also has drift-era files to purge (Task 1 lists them).
- Root `LICENSE` is LGPL v2.1 (Phase 0 finding; corrected in Task 13).
- `fprettify` available in pixi; `pixi run lint` exists but has no config file.
- FORD project file at repo root as `docs.md`.
- `check-fast` passes 4/4, `check-full` passes 6/6.

**End state:**
- `docs/` = `index.md`, `architecture.md`, `state-management.md`, `configuration-schema.md`, `build-and-test.md`, `dependency-management.md`, `code-style.md`, `contributing.md`, `adr/0001…0005`, plus preserved `superpowers/`, `public/`, `api/`.
- `.fprettify.rc` at repo root encodes the style rules `code-style.md` declares.
- 12 `src/<domain>/README.md` files + refreshed `src/README.md`.
- Root `LICENSE` = GPL v2; `NOTICE.md` clarifies the dual-license situation.
- `docs.md` excludes `docs/superpowers/` from FORD output; `pixi run -e docs docs-build` succeeds.
- Tag `rescue/phase-2-docs` on `main` = `development`.

**Deferred out of Phase 2:**
- `docs/api/` tracking policy (currently tracked + gitignored — inconsistent; fix in a later chore).
- Per-procedure `!>` docstrings inside Fortran source (Phase 4, as modules are touched).
- pFUnit proper vendoring.
- Markdown linter / link checker for docs (new deps — YAGNI during rescue).

**Conventions:**
- One concern per commit.
- Every commit keeps `pixi run -e test check-fast` green.
- No pushes to any remote during Phase 2.
- Before tagging: clean build + `check-full` green + `pixi run -e docs docs-build` succeeds.

---

## Task 1: Purge drift-era documentation

**Files:**
- Delete: `docs/refactoring.md`, `docs/refactoring-diary.md`, `docs/next-refactoring-task.md`, `docs/newschema.md`, `docs/state_management_architecture.md`, `docs/state_management_pattern.md`, `docs/style_guide.md`, `docs/architecture/functions.md`, `docs/architecture/index.md`, `docs/architecture/modules.md`, `docs/architecture/subroutines.md`.
- Preserve (for later tasks): `docs/code-style-guide.md` (Task 8), `docs/index.md` (Task 2 rewrites in place), `docs/pfunit-vendoring.md` (Task 6).

- [ ] **Step 1: Confirm starting state**

```bash
cd /home/zawadzkim/Code/swap
git rev-parse HEAD               # 17a49f2
git log --oneline -1 origin/main  # 7587ca3
git branch --show-current          # development
```

- [ ] **Step 2: Delete the files**

```bash
cd /home/zawadzkim/Code/swap
git rm docs/refactoring.md \
       docs/refactoring-diary.md \
       docs/next-refactoring-task.md \
       docs/newschema.md \
       docs/state_management_architecture.md \
       docs/state_management_pattern.md \
       docs/style_guide.md \
       docs/architecture/functions.md \
       docs/architecture/index.md \
       docs/architecture/modules.md \
       docs/architecture/subroutines.md
rmdir docs/architecture
```

- [ ] **Step 3: Commit**

```bash
git commit -m "docs: purge drift-era documentation pages

Remove seven top-level drift-era markdown files and the auto-
generated-looking docs/architecture/ subdirectory. Content worth
keeping is rewritten in later Phase 2 tasks (index.md stays and
gets rewritten in place in Task 2; architecture.md is new in
Task 3; state-management.md is new in Task 4; code-style.md is new
in Task 8 and absorbs code-style-guide.md).

Original content preserved on archive/wip-drifted if ever needed."
```

- [ ] **Step 4: Verify check-fast**

```bash
pixi run -e test check-fast 2>&1 | tail -5
```

Expected: 4/4 pass.

---

## Task 2: Rewrite `docs/index.md`

**Files:**
- Modify: `docs/index.md`

- [ ] **Step 1: Overwrite `docs/index.md` with this content**

```markdown
---
title: SWAP documentation
author: SWAP Team
---

# SWAP — Soil Water Atmosphere Plant

SWAP is a one-dimensional vertical simulation model for transport processes (water, heat, solutes) in the soil–plant–atmosphere continuum at field scale. This repository holds the modernization of SWAP 4.2.0 under active rescue-and-stabilize work.

## Start here

1. [Architecture overview](architecture.html) — three-phase control flow, state aggregation, module boundaries.
2. [Build and test](build-and-test.html) — build, fast/full test split, regression harness, fixture policy.
3. [Contributing](contributing.html) — branch model, commit conventions, verification gates.

## Reference

- [State management](state-management.html) — `config_t` / `initial_t` / `state_t`, ASSOCIATE, lifecycle.
- [Configuration schema](configuration-schema.html) — TOML input reference.
- [Dependency management](dependency-management.html) — subprojects, pixi deps, pFUnit, how to bump.
- [Code style](code-style.html) — Fortran 2008 conventions, names, intent; matched by `.fprettify.rc`.
- [Architecture decision records](adr/) — the non-obvious choices and why.

## Current status

The modernization is in a rescue-and-stabilize workflow; see `docs/superpowers/specs/` and `docs/superpowers/plans/` for the internal planning. Status of the rescue is tracked in git tags named `rescue/phase-N-*`.

## API reference

FORD-generated API reference at [api/index.html](api/index.html). Run `pixi run -e docs docs-build` to regenerate.
```

- [ ] **Step 2: Commit**

```bash
cd /home/zawadzkim/Code/swap
git add docs/index.md
git commit -m "docs: rewrite index.md as a short FORD landing page

Links to the eight top-level canonical docs added in Phase 2
(architecture, build-and-test, contributing, state-management,
configuration-schema, dependency-management, code-style, ADRs).
Removes drift-era stub content."
```

---

## Task 3: Write `docs/architecture.md`

**Files:**
- Create: `docs/architecture.md`

- [ ] **Step 1: Reconnaissance**

```bash
cd /home/zawadzkim/Code/swap
grep -nE "^\s*subroutine|^\s*module" src/core/swap.f90 | head -20
grep -nE "^\s*subroutine|^\s*module" src/core/swap_main.f90 | head -10
head -80 src/core/swap_state_mod.f90
```

Record the three-phase task=1/2/3 structure and `swap_state_t`'s member fields.

- [ ] **Step 2: Write `docs/architecture.md`**

Sections (prose, ~200–400 words each):

1. **Overview** — what SWAP simulates; 1-D vertical column; daily time stepping.
2. **Three-phase control flow** — `swap(iTask)` with 1 (init) / 2 (dynamic step) / 3 (closure). Physics ordering inside task=2: meteo → crop → irrigation → root-extraction → bottom BC → drainage → surface-water → soilwater (Richards via TDMA) → heat → solute → macropore → output. `swap_main.f90` drives reruns; `swap.f90` holds `swap(iTask)`.
3. **State aggregation** — `swap_state_t` in `src/core/swap_state_mod.f90` composes domain states (atmosphere, boundary, drainage, surfacewater, heat, solute, soil, macropore, crop). Explicit state passing replaces the legacy global `variables` module. Link to `state-management.md`.
4. **I/O layer** — TOML via toml-f is canonical (`src/io/readswaptoml.f90`, `readdrainagetoml.f90`). Legacy fixed-format (`readswap.f90`) coexists. Output via `swap_csv_output.f90` + domain writers. Link to `configuration-schema.md`.
5. **Subdirectory map** — one paragraph per `src/<domain>/` (atmosphere, boundary, core, crop, drainage, error, heat, io, macropore, soil, solute, utils). Note `src/error/` is an empty Phase 4 placeholder.
6. **Dependency direction** — core → domain states → physics routines → utils. Domain subdirs don't import each other directly; they communicate via `swap_state_t`. `src/core/swap_state_mod.f90` is the only module aware of all domain types.
7. **Rescue caveat** — see `docs/superpowers/specs/` for the rescue spec and `docs/adr/` for architectural decisions.

- [ ] **Step 3: Commit**

```bash
cd /home/zawadzkim/Code/swap
git add docs/architecture.md
git commit -m "docs: add architecture.md — three-phase flow, state aggregator, I/O"
```

---

## Task 4: Write `docs/state-management.md`

**Files:**
- Create: `docs/state-management.md`

- [ ] **Step 1: Reconnaissance**

```bash
cd /home/zawadzkim/Code/swap
head -100 src/atmosphere/atmosphere_state.f90
grep -nE "ASSOCIATE|associate|_state_init|_state_finalize" src/atmosphere/atmosphere_state.f90 src/drainage/drainage_state.f90 | head -20
head -50 src/core/swap_state_mod.f90
head -30 src/core/swap_state_sync.f90
```

- [ ] **Step 2: Write the document**

Sections:

1. **Context** — legacy `variables` globals problem; aggregator pattern solution.
2. **Three type categories** — `config_t` (intent(in) everywhere), `initial_t` (frozen after init), `state_t` (mutated in dynamic loop). A domain may have one, two, or all three.
3. **Canonical module layout** — for domain `foo`: `src/foo/foo_state.f90` with `module foo_state_mod`, `type foo_state_t`, `foo_state_init`, `foo_state_finalize`, optional `foo_state_reset_*`, optional `foo_state_from_variables` / `foo_state_to_variables` bridges.
4. **Aggregation** — `swap_state_t` in `src/core/swap_state_mod.f90` composes everyone.
5. **ASSOCIATE** — example:
   ```fortran
   associate(h => state%soil%h, dz => state%soil%dz)
       ! ... block body ...
   end associate
   ```
   Scoped only; does not persist across calls.
6. **Lifecycle contract** — init exactly once per simulation; finalize at shutdown; resets between timesteps/days. Never partial state.
7. **How to add a new domain state** — create the module, add to `swap_state_t`, wire init/finalize in `src/core/swap.f90` task=1 / task=3, add to `tests/unit/meson.build`'s `test_base_sources`.
8. **Legacy bridge** — `src/core/swap_state_sync.f90` is a rescue-era shim between aggregator state and the legacy `variables` module. Temporary; removed in Phase 4.
9. **Reference** — see `docs/adr/0003-aggregator-state-over-compartments-for-now.md` for why this pattern instead of compartment-based state.

- [ ] **Step 3: Commit**

```bash
cd /home/zawadzkim/Code/swap
git add docs/state-management.md
git commit -m "docs: add state-management.md — config/initial/state, ASSOCIATE, lifecycle"
```

---

## Task 5: Write `docs/configuration-schema.md`

**Files:**
- Create: `docs/configuration-schema.md`

- [ ] **Step 1: Reconnaissance — enumerate the real schema from the readers**

```bash
cd /home/zawadzkim/Code/swap
grep -nE "get_value|get_array|section" src/io/readswaptoml.f90 | head -60
grep -nE "get_value|get_array|section" src/io/readdrainagetoml.f90 | head -40
find tests/swap-cases -name '*.toml' -not -path '*/result*' | head
```

Readers are the authoritative schema. Note: TOML crop reader does NOT exist at baseline; crop configs are still fixed-format `.crp`.

- [ ] **Step 2: Write the document**

Sections:

1. **Input file layout** — `.swp` (main, TOML), `.dra` (drainage, TOML), `.crp` (crop, fixed-format at baseline), `.met` (meteo, fixed-format), optional `.bbc` (bottom boundary).
2. **TOML conventions** — sections, subsections, arrays-of-tables, scalar types. toml-f is stricter than YAML.
3. **`.swp` reference** — enumerate every section and key from `readswaptoml.f90`'s `get_value`/`get_array` calls. For each: type, required-or-optional, description, example value. Structure as a table per section.
4. **`.dra` reference** — same treatment from `readdrainagetoml.f90`.
5. **`.crp` reference (baseline)** — placeholder: "TOML crop reader is a Phase 4 deliverable. At baseline crop configs use fixed-format `.crp`; see `readcropfixed` / `readgrass` / `readwofost` subroutines for the legacy format. This section will be filled in when Phase 4 adds `readcrop_toml.f90`."
6. **Example** — a minimal valid `.swp` that runs `hupselbrook`. If a TOML variant exists under `tests/swap-cases/*/`, reference it; otherwise inline the smallest example compatible with the reader.
7. **Validation** — missing required keys → read-time error. Unknown keys → silently ignored today (Phase 4 adds strict validation).

- [ ] **Step 3: Commit**

```bash
cd /home/zawadzkim/Code/swap
git add docs/configuration-schema.md
git commit -m "docs: add configuration-schema.md — TOML input reference

Section-by-section reference for .swp and .dra TOML inputs,
extracted from readswaptoml.f90 / readdrainagetoml.f90 (authoritative
schema sources). Crop TOML is a Phase 4 deliverable — baseline
crop configs still use fixed-format .crp."
```

---

## Task 6: Write `docs/build-and-test.md`; absorb `pfunit-vendoring.md`

**Files:**
- Create: `docs/build-and-test.md`
- Delete: `docs/pfunit-vendoring.md`
- Modify: `.github/agents/swap-fortran.agent.md` (update the pFUnit-vendoring pointer)

- [ ] **Step 1: Read pfunit-vendoring.md for content to absorb**

```bash
cat /home/zawadzkim/Code/swap/docs/pfunit-vendoring.md
```

- [ ] **Step 2: Write `docs/build-and-test.md`**

Sections:

1. **Quick start** — `pixi run -e test check-fast`. Build + pFUnit (empty at baseline) + four fast cases.
2. **Prerequisites** — pixi handles everything (python, gfortran, meson, ninja, fprettify, pyswap, FORD).
3. **Build** — `pixi run build-linux` → `builddir/swap`. Single `builddir/` per ADR 0002. `pixi run clean` wipes.
4. **Fast vs full test protocol** — the canonical gates:
   - `check-fast`: build + pFUnit + 4 cases (hupselbrook, surfacewater, salinitystress, grassgrowth). Target <90s. Every iteration.
   - `check-full`: + oxygenstress + macropore. Target ~10 min. Before every phase tag.
   - `regression [case-name ...]`: any subset by name.
   - `regression --regenerate-fixtures`: rewrite `*_expected_gfortran.json`; use only after a deliberate physics change; see `contributing.md` for the policy.
5. **Regression harness internals** — `tests/regression/test_output_regression.py`. Each case in isolated temp dir. Aggregates annual flux sums / state means / cumulative values. Tolerance 1e-2. Cases defined in top-of-file `CASES` dict.
6. **Fixture policy** — `*_expected.json` = historical ifx reference (kept, not compared). `*_expected_gfortran.json` = what harness uses. `tests/regression/INVESTIGATION_NOTES.md` documents the macropore DRAINAGE divergence and oxygenstress MOWDM deviation.
7. **`tests/swap-cases/` submodule** — pinned to `main`. Dirties after each regression run (harness renames input files); outer-repo pointer unchanged. `git submodule update tests/swap-cases` resets it.
8. **pFUnit** — (absorb pfunit-vendoring.md content here):
   - `tests/pFUnit/` is a gitlink (mode 160000, SHA pinned) NOT declared in `.gitmodules`. So `git clone` of the outer repo does NOT auto-populate it.
   - Install dir consumed via `PFUNIT_ROOT` (set in pixi `[activation]`).
   - To rebuild if lost: `cd tests/pFUnit && mkdir build && cd build && cmake -DCMAKE_INSTALL_PREFIX=install_gfortran .. && make -j && make install`.
   - Full vendoring as a meson subproject is deferred.
9. **Cross-compilation** — `pixi run build-windows-cross` best-effort via mingw.
10. **Cleaning** — `pixi run clean` removes build dirs.
11. **Common problems** — "Unsupported Fortran compiler" (Intel sourced; see ADR 0001); "pFUnit not found" (rebuild); "No suitable tests defined." (expected at baseline).

- [ ] **Step 3: Delete pfunit-vendoring.md**

```bash
cd /home/zawadzkim/Code/swap
git rm docs/pfunit-vendoring.md
```

- [ ] **Step 4: Update agent pointer**

Edit `.github/agents/swap-fortran.agent.md`. Replace:
```
6. `docs/pfunit-vendoring.md` — the pFUnit gitlink peculiarity and how to rebuild the install.
```
with:
```
6. `docs/build-and-test.md` — build/test workflow, fast-vs-full gates, pFUnit gitlink peculiarity with rebuild steps.
```

- [ ] **Step 5: Commit**

```bash
git add docs/build-and-test.md .github/agents/swap-fortran.agent.md
git status --short
git commit -m "docs: add build-and-test.md; absorb pfunit-vendoring.md

Single canonical build/test doc covering quick start, fast vs full
protocol, regression harness internals, fixture policy, pFUnit
(gitlink peculiarity + rebuild steps absorbed from deleted
pfunit-vendoring.md), submodule dirtying, cleaning, common
problems. Agent definition updated."
```

---

## Task 7: Write `docs/dependency-management.md`

**Files:**
- Create: `docs/dependency-management.md`

- [ ] **Step 1: Reconnaissance**

```bash
cd /home/zawadzkim/Code/swap
ls subprojects/
for f in subprojects/*.wrap; do echo "=== $f ==="; cat "$f"; done
grep -nE "^\[dependencies\]|^\[feature" pixi.toml | head
head -50 pixi.toml
```

- [ ] **Step 2: Write the document**

Sections:

1. **Two layers** — pixi (toolchain); meson subprojects (Fortran libs).
2. **pixi deps** — enumerate from `[dependencies]`, `[feature.test.dependencies]`, `[feature.docs.dependencies]`:
   - `meson`, `gfortran`, `ninja`, `python 3.11` — build chain + test harness runtime.
   - `fprettify` — Fortran formatter (see `code-style.md` and `.fprettify.rc`).
   - `pytest`, `pytest-xdist`, `pytest-timeout`, `pandas`, `pyswap` — test feature.
   - `ford` — docs feature.
3. **Meson subprojects** — list each from `subprojects/*.wrap`:
   - `ttutil` — legacy utilities; vendored from SWAP-model/ttutil.
   - `toml-f` — pure-Fortran TOML parser.
   - `test-drive` — test harness helpers (indirect).
4. **pFUnit** — cross-reference `build-and-test.md` for the gitlink state.
5. **Bumping versions** — per layer:
   - pixi: edit `pixi.toml`, `pixi install`, commit both `pixi.toml` and `pixi.lock`.
   - Subprojects: edit `subprojects/<name>.wrap`, clear subproject, reconfigure, run `check-full`, commit.
   - pFUnit: rebuild install, bump gitlink SHA if needed.
6. **No new deps during rescue** — one variable at a time; proposals wait for Phase 4 exit.

- [ ] **Step 3: Commit**

```bash
cd /home/zawadzkim/Code/swap
git add docs/dependency-management.md
git commit -m "docs: add dependency-management.md"
```

---

## Task 8: Write `docs/code-style.md` (replaces `code-style-guide.md`)

**Files:**
- Create: `docs/code-style.md`
- Delete: `docs/code-style-guide.md`

- [ ] **Step 1: Read existing code-style-guide.md for anything worth preserving**

```bash
cat /home/zawadzkim/Code/swap/docs/code-style-guide.md
```

- [ ] **Step 2: Write `docs/code-style.md`**

Sections:

1. **Language target** — Fortran 2008; `-std=legacy` tolerates old idioms in legacy modules; new code is 2008 free-form.
2. **Module conventions** —
   - One module per file. `module foo_mod` / `end module foo_mod`.
   - `implicit none` immediately after `module`.
   - `private` default; explicit `public ::` for exports.
   - No COMMON, no `.fi` includes in new code.
3. **Types and intent** —
   - Derived types for state/config/initial: `type :: foo_state_t`.
   - Every argument has explicit `intent(in|out|inout)`.
   - `real(kind=8)` acceptable during rescue; `real(kind=real64)` from `iso_fortran_env` preferred in new code.
4. **Naming** —
   - Types: `<domain>_state_t`, `<domain>_config_t`, `<domain>_initial_t`.
   - Modules: `<domain>_mod` (→ `foo_state_mod`).
   - Procedures: `<domain>_state_init`, `<domain>_state_finalize`, `<domain>_process_<verb>`.
   - Variables: snake_case. Prefer full words unless it's a physics convention (`theta`, `psi`, `K`).
5. **Control flow** — `select case` over `if/elseif` chains. Early returns fine. No `goto` in new code.
6. **ASSOCIATE** — scoped aliasing for repeated `state%field` references; not a substitute for extracting a helper.
7. **Comments / docstrings** — FORD consumes `!>` above procedures/types. Every new public procedure gets a one-line summary + arg descriptions. Explain WHY, not WHAT.
8. **Formatting** — **`pixi run lint` (fprettify) is authoritative.** The rules are encoded in `.fprettify.rc`; see Task 9 for the config. If `lint` disagrees with something written here, `.fprettify.rc` wins — open an issue to reconcile rather than ignore the formatter.
9. **Tests** — new public procedures come with a pFUnit test; see `build-and-test.md` and `contributing.md`.
10. **Legacy rule** — `src/core/variables.f90` and `src/core/swap_state_sync.f90` are rescue-era. Do not extend; new features wrap through `swap_state_t`.

- [ ] **Step 3: Delete old guide**

```bash
cd /home/zawadzkim/Code/swap
git rm docs/code-style-guide.md
```

- [ ] **Step 4: Commit**

```bash
git add docs/code-style.md
git status --short
git commit -m "docs: add code-style.md; delete drift-era code-style-guide.md

Canonical Fortran-style reference matching the .fprettify.rc
config added in the next task (Task 9). When pixi run lint
disagrees with this doc, the config wins — open an issue to
reconcile rather than ignore the formatter."
```

---

## Task 9: Configure `.fprettify.rc` to match `code-style.md`

**Goal:** The linter in pixi (`pixi run lint`) must reflect what `code-style.md` declares. Today `pixi run lint` runs `fprettify` with no config file, so it uses defaults that may not match our style. Create `.fprettify.rc` at repo root and extend `lint` to target Fortran sources with sensible settings.

**Files:**
- Create: `.fprettify.rc`
- Modify: `pixi.toml` (extend `lint` task to target src recursively, add `lint-check` dry-run variant)
- Modify: `docs/code-style.md` (reference the exact `.fprettify.rc` location and explain dry-run)

- [ ] **Step 1: Inspect current lint task and fprettify defaults**

```bash
cd /home/zawadzkim/Code/swap
grep -n "fprettify\|lint " pixi.toml
pixi run fprettify --help 2>&1 | head -40
```

Record: default indent (3?), default whitespace/case style, recognized config file name (`.fprettify.rc`).

- [ ] **Step 2: Write `.fprettify.rc`**

Write to `/home/zawadzkim/Code/swap/.fprettify.rc`:

```ini
# fprettify configuration for the SWAP modernization.
# Matches docs/code-style.md. If the formatter disagrees with the
# doc, this file wins — open an issue to reconcile rather than
# ignore the formatter.
[fprettify]
indent = 4
line-length = 132
whitespace = 2
whitespace-comma = true
whitespace-assignment = true
whitespace-decl = true
whitespace-relational = true
whitespace-logical = true
whitespace-plusminus = true
whitespace-multdiv = false
whitespace-print = true
whitespace-type = true
whitespace-intrinsics = true
case = [1, 1, 1, 1]
strict-indent = false
enable-decl = true
enable-replacements = false
disable-indent-mod = true
disable-whitespace-mod = false
```

Notes on the chosen values:

- `indent = 4` — matches modern Fortran style and the visible spacing in most `src/` files.
- `line-length = 132` — Fortran 2008 free-form maximum; matches the project's convention of allowing wide lines over continuation gymnastics.
- `case = [1, 1, 1, 1]` — lower-case keywords, intrinsics, procedures, user-defined. Consistent with snake_case naming.
- `enable-replacements = false` — do NOT auto-rewrite operators (e.g., `.eq.` → `==`). Legacy SWAP code has plenty of old-style comparisons; mechanical rewriting invites spurious physics-affecting diffs. When Phase 4 touches a module, the engineer can enable replacements manually file-by-file.
- `disable-indent-mod = true` — don't re-indent module scopes; some legacy files have hand-crafted alignment worth preserving in rescue-era commits.

If `pixi run fprettify --help` reveals options with different names than expected, adjust. The goal is: the committed `.fprettify.rc` produces zero diff when applied to the current `src/` tree (i.e., `fprettify --diff` produces no output for files not actively being rewritten), and any style choice a contributor reads in `code-style.md` is enforceable with `pixi run lint`.

- [ ] **Step 3: Extend the `lint` task in `pixi.toml`**

Find the current `lint` task line (under `[tasks]`):

```toml
lint = { cmd = "fprettify" }
```

Replace with two tasks:

```toml
# Format Fortran sources in place per .fprettify.rc.
lint = { cmd = "fprettify --config-file .fprettify.rc -r src tests/unit" }

# Dry-run: check what lint would change without writing. Use in CI
# or locally before committing to see diffs.
lint-check = { cmd = "fprettify --config-file .fprettify.rc -d -r src tests/unit" }
```

If `fprettify` does not accept `-r` for recursive or the flag differs, adjust to the correct syntax (visible via `fprettify --help`).

- [ ] **Step 4: Dry-run to see what would change**

```bash
cd /home/zawadzkim/Code/swap
pixi run lint-check 2>&1 | tee /tmp/phase2-task9-lintcheck.log | tail -30
```

Possible outcomes:

- **Zero diff** — perfect; the rescue baseline already conforms to the chosen style. Proceed to Step 7.
- **Small diff, mechanical (whitespace only)** — proceed to Step 5 to apply.
- **Large diff, touching hundreds of files** — STOP. The chosen config is too aggressive for the rescue. Tighten the config to produce minimal diff and try again (common culprits: `enable-replacements`, `strict-indent`, `case`). This is NOT a Phase 2 physics-reformat exercise.
- **Diff that reformats physics-sensitive patterns** — STOP. Adjust config so it does not touch that.

- [ ] **Step 5: Apply lint if the diff is small and safe**

If Step 4's diff is small/mechanical:

```bash
cd /home/zawadzkim/Code/swap
pixi run lint 2>&1 | tee /tmp/phase2-task9-lintapply.log | tail
git diff --stat
```

Expected: modifications to `.f90` files in `src/` and possibly `tests/unit/`. Review the diff: it should be whitespace, line continuation, keyword case. No content changes (no removed/added logic, no renamed identifiers, no operator swaps).

If the diff contains anything beyond whitespace/case, STOP — reconfigure `.fprettify.rc` more conservatively.

- [ ] **Step 6: Verify regression still green after lint**

```bash
cd /home/zawadzkim/Code/swap
rm -rf builddir
pixi run -e test check-full 2>&1 | tail -30
```

Expected: 6/6 pass. If not, the linter re-wrote something that changes behaviour — revert and tighten the config.

If `check-full` is green, proceed.

- [ ] **Step 7: Update `docs/code-style.md` Section 8 to point at `.fprettify.rc` and explain `lint-check`**

Replace Section 8 of `docs/code-style.md` with:

```markdown
8. **Formatting** — `pixi run lint` (fprettify) is authoritative. The rules are encoded in `.fprettify.rc` at the repository root. Run `pixi run lint-check` first to see diffs without applying; `pixi run lint` rewrites files in place.

If the formatter disagrees with something written here, `.fprettify.rc` wins — open an issue to reconcile rather than ignore the formatter. The linter is deliberately conservative during the rescue: `enable-replacements` is off so legacy operator forms (`.eq.`, `.le.`) are not auto-rewritten, and `disable-indent-mod` is on so hand-crafted module-scope alignment is preserved.
```

- [ ] **Step 8: Commit**

```bash
cd /home/zawadzkim/Code/swap
git add .fprettify.rc pixi.toml docs/code-style.md
# Also stage any whitespace-only diffs produced by Step 5 if applicable:
git add -u src/ tests/unit/
git status --short
git commit -m "$(cat <<'EOF'
style(lint): configure fprettify to match code-style.md

Add .fprettify.rc at repo root encoding the style rules that
docs/code-style.md documents. Extend the pixi lint task to target
src/ and tests/unit/ recursively; add lint-check as a dry-run
variant.

Config is deliberately conservative during the rescue:
- enable-replacements = false (do NOT rewrite .eq. to == etc.;
  legacy operator forms are left alone until Phase 4 touches them).
- disable-indent-mod = true (hand-crafted module-scope alignment
  preserved).
- case = [1,1,1,1] (lower-case keywords/intrinsics/procedures/
  user-defined, matching snake_case naming).

If whitespace-only diffs were produced by applying lint to the
current tree, they are included here — verified check-full still
6/6 green, so no behavioural change.

docs/code-style.md Section 8 now points at .fprettify.rc and
documents the lint-check dry-run workflow.
EOF
)"
```

---

## Task 10: Write `docs/contributing.md`

**Files:**
- Create: `docs/contributing.md`

- [ ] **Step 1: Write the document**

Sections:

1. **You are here** — repo in rescue-and-stabilize. Read the spec at `docs/superpowers/specs/...` before substantive changes.
2. **Branch model** —
   - `main` advances only on phase completion.
   - `development` — in-phase commits.
   - `archive/*` — pre-rescue snapshots.
   - `legacy/swap-4.2.0` — orphan with upstream SWAP 4.2.0.
   - Phase 4+: per-change feature branches off `development`.
   - `origin/main` is NOT advanced during the rescue.
3. **Commit conventions** — conventional-commits prefixes (`feat:`, `fix:`, `docs:`, `refactor:`, `chore:`, `test:`, `style:`), optional scopes, imperative mood, body wrapped at 72 cols.
4. **Verification gates** —
   - `pixi run -e test check-fast` before every merge to `development`.
   - `pixi run lint-check` before every commit (or `pixi run lint` to auto-fix).
   - `pixi run -e test check-full` before every phase tag.
   - Tag names: `rescue/phase-N-<shortname>`.
5. **Fixture policy** —
   - Harness compares against `*_expected_gfortran.json`.
   - `*_expected.json` = historical ifx reference (never delete).
   - Deliberate physics changes: `python tests/regression/test_output_regression.py --regenerate-fixtures`. Document before/after numeric diff for at least one affected variable in the commit message.
   - Never regenerate silently to make a failing test pass without understanding why.
6. **Accepted tolerances** — macropore DRAINAGE divergence; oxygenstress MOWDM deviation. See `tests/regression/INVESTIGATION_NOTES.md`.
7. **Style** — run `pixi run lint` before committing. `.fprettify.rc` is authoritative.
8. **Adding new code** — cross-references `state-management.md` (Section 7), `configuration-schema.md`, `build-and-test.md`.
9. **Agents** — `.github/agents/swap-fortran.agent.md` specifies required context for automated agents.
10. **When in doubt** — match surrounding code; read adjacent ADRs; open a discussion before guessing.

- [ ] **Step 2: Commit**

```bash
cd /home/zawadzkim/Code/swap
git add docs/contributing.md
git commit -m "docs: add contributing.md — branch model, commits, fixture policy, gates"
```

---

## Task 11: ADR 0003 — aggregator-state-over-compartments-for-now

**Files:**
- Create: `docs/adr/0003-aggregator-state-over-compartments-for-now.md`

- [ ] **Step 1: Write the ADR**

Standard ADR structure (Context / Decision / Consequences / Re-visit). Key points:

- Legacy `variables` module is the problem: globally mutable state.
- Original plan (in deleted `next-refactoring-task.md`): compartment-based state.
- At baseline: aggregator pattern is partially done.
- Decision: finish the aggregator through Phase 4; compartment is a follow-on spec.
- Rationale: bounded work; compartment needs physics-aware TDD; aggregator already enables testing.
- Consequences: some throwaway work; `variables.f90` + `swap_state_sync.f90` stay until Phase 4 exit (flagged legacy in `code-style.md`).

Full text (write verbatim):

```markdown
# ADR 0003 — Aggregator state now, compartment state later

Status: accepted (2026-04-23, during rescue Phase 2)

## Context

The legacy SWAP codebase used a global `variables` module that every physics routine read and wrote. This blocks parallelization, breaks test isolation, and makes the dependency graph opaque.

The rescue spec's original forward plan (captured in the now-deleted `docs/next-refactoring-task.md`) was to skip directly to compartment-based state (ponding, canopy, snow, soil, saturated zone, etc.) where process functions explicitly pass water/heat/solute between compartments. This is the long-term target.

At the rescue baseline `e256bc0` the codebase had progressed partway through an aggregator pattern (one `swap_state_t` composed of domain-specific `*_state_t` types) without completing it. Multiple modules still `use variables` and rely on legacy globals.

## Decision

The rescue stays on the aggregator pattern through Phase 4 exit. Compartment-based state is a follow-on spec.

Rationale:

1. Aggregator is partially done; finishing it is bounded work.
2. Compartment state is a larger physics-aware redesign; attempting it mid-rescue risks physics regressions we cannot easily detect.
3. The aggregator already enables test isolation and state passing; it is sufficient for Phase 4's fix-and-clean work.
4. Compartment state will require TDD-first rewrite per compartment with regression checks; that belongs to its own dedicated spec.

## Consequences

- **Positive**: the rescue has a bounded, achievable goal. Phase 4 closes the pre-compartment gates (every module has clean `config_t` / `initial_t` / `state_t` separation) rather than attempting a physics-level redesign.
- **Positive**: Phase 3 test coverage builds on the aggregator; state-type lifecycle tests ported during Phase 4 remain useful when compartment state lands later.
- **Negative**: aggregator is not the final state. Phase 4 produces code that gets rewritten again when compartment state lands. Some work is "throwaway".
- **Negative**: `src/core/variables.f90` and `src/core/swap_state_sync.f90` stay in the tree until the aggregator is complete. They are flagged as legacy in `docs/code-style.md`.

## When to revisit

After `rescue/complete`. A dedicated compartment-state spec and plan supersede this ADR.
```

- [ ] **Step 2: Commit**

```bash
cd /home/zawadzkim/Code/swap
git add docs/adr/0003-aggregator-state-over-compartments-for-now.md
git commit -m "docs(adr): ADR 0003 — aggregator state now, compartment state later"
```

---

## Task 12: ADR 0004 — pFUnit for Fortran unit tests

**Files:**
- Create: `docs/adr/0004-pfunit-for-unit-tests.md`

- [ ] **Step 1: Write the ADR**

```markdown
# ADR 0004 — pFUnit for Fortran unit tests

Status: accepted (2026-04-23, during rescue Phase 2; wiring finalized in Phase 1)

## Context

Options considered at rescue start:

- **pFUnit 4** — mature, CMake-based, preprocessor-generated Fortran test drivers. De-facto standard in Fortran scientific computing (NASA, GFDL, ECMWF use it).
- **test-drive** (Fortran-Lang) — lightweight, pure Fortran, no preprocessor. Easier to vendor but less feature-complete (no parameterized tests, no fixtures in the pFUnit sense).
- **Handroll** — a few dozen lines of `call assert_equal(...)` helpers. Maximum control, minimum ergonomics.

At the rescue baseline, `tests/unit/meson.build` was already wired for pFUnit 4.15, though all suites were orphaned against drift-era modules (Phase 1 removed them). The scaffolding (generator, driver, `testSuites.inc`) remained functional.

## Decision

Use pFUnit 4 for Fortran unit tests throughout the rescue and beyond.

## Consequences

- **Positive**: ecosystem standard; conventions transfer to and from other scientific Fortran projects.
- **Positive**: rich fixtures and parameterized tests — useful for state lifecycles and physics routines across parameter combinations.
- **Positive**: wiring already in place; adding a new test is `.pf` file + `ADD_TEST_SUITE()` line.
- **Negative**: pFUnit is CMake-based, not meson-native. The rescue lives with a gitlink checkout at `tests/pFUnit/` (see `docs/build-and-test.md`). Proper vendoring is deferred.
- **Negative**: preprocessor adds Python compile-time dependency; pixi manages this, so net friction is zero.

## Alternative if circumstances change

test-drive as fallback. Conversion cost: rewriting `.pf` suites in test-drive's Fortran-native API.
```

- [ ] **Step 2: Commit**

```bash
cd /home/zawadzkim/Code/swap
git add docs/adr/0004-pfunit-for-unit-tests.md
git commit -m "docs(adr): ADR 0004 — pFUnit for Fortran unit tests"
```

---

## Task 13: ADR 0005 — licensing audit, and remediate root LICENSE

**Files:**
- Create: `docs/adr/0005-licensing-audit.md`
- Modify: `LICENSE` (replace with GPL v2 text)
- Create: `NOTICE.md`

This is split into ADR-write + license-remediate in one commit because the ADR would be hollow without the fix immediately following.

- [ ] **Step 1: Inspect current root LICENSE and upstream GPL v2**

```bash
cd /home/zawadzkim/Code/swap
head -10 LICENSE
wc -l LICENSE
git show legacy/swap-4.2.0:source_swp_4.2.0/LICENSE | head -20
```

Expected: root LICENSE is LGPL v2.1 (179 lines, first line references lgpl-2.1.en.html). Upstream SWAP 4.2.0 LICENSE declares GPL v2 (a short notice file, not the full license text).

- [ ] **Step 2: Obtain the full GPL v2 text**

The short upstream notice is not sufficient as a root LICENSE (it's a declaration, not the license text). The full GPL v2 is 339 lines of the canonical https://www.gnu.org/licenses/old-licenses/gpl-2.0.txt content.

Pragmatic options:

1. If network is available: `curl -sSL https://www.gnu.org/licenses/old-licenses/gpl-2.0.txt -o LICENSE`.
2. If network is not available: STOP and report BLOCKED asking the controller to paste the GPL v2 text. The subagent must NOT invent or paraphrase license text.

Verify the downloaded content begins with `GNU GENERAL PUBLIC LICENSE / Version 2, June 1991` and is approximately 339 lines.

```bash
head -3 LICENSE
wc -l LICENSE
```

- [ ] **Step 3: Write `NOTICE.md`**

Write to `/home/zawadzkim/Code/swap/NOTICE.md`:

```markdown
# License notices

## SWAP modernization (this repository)

Licensed under the **GNU General Public License Version 2**. See `LICENSE` for the full text.

This is a derivative work of SWAP 4.2.0 by Wageningen University & Research (WUR), originally distributed under GPL v2. Per GPL v2 §2(b), the derivative work is licensed on the same terms.

The canonical upstream SWAP 4.2.0 tree is preserved on the `legacy/swap-4.2.0` orphan branch in this repository for reference.

## TTUTIL 4.27 (bundled)

Portions of this repository incorporate TTUTIL 4.27 by WUR, distributed under the **GNU Lesser General Public License Version 2.1**. TTUTIL source files carry their own LGPL v2.1 notice at build time via `subprojects/ttutil/`, and are preserved at `legacy/swap-4.2.0:source_ttutil_4.27/`.

## History

Prior to rescue Phase 2 (ADR 0005), the root `LICENSE` file in this repository contained the LGPL v2.1 text — an inherited misnomer. ADR 0005 documents the correction. Users who pulled the repository during that period hold a copy under the notice that was present at the time.
```

- [ ] **Step 4: Write ADR 0005**

Write to `/home/zawadzkim/Code/swap/docs/adr/0005-licensing-audit.md`:

```markdown
# ADR 0005 — Licensing audit: root LICENSE vs upstream GPL v2

Status: accepted (2026-04-23, during rescue Phase 2); remediation committed in the same task as this ADR

## Context

Phase 0 baseline verification flagged an inconsistency:

- Repository-root `LICENSE` was the **GNU Lesser General Public License v2.1**. 179 lines; first line referenced `lgpl-2.1.en.html`.
- Upstream SWAP 4.2.0 (preserved on `legacy/swap-4.2.0` at `source_swp_4.2.0/LICENSE`) declares the **GNU General Public License v2**. TTUTIL 4.27 is separately LGPL v2.1.

The modernization tree in `src/` is a derivative work of the upstream GPL v2 SWAP code. Under GPL v2 §2(b) a derivative work must be licensed "as a whole at no charge to all third parties under the terms of this License" — i.e., GPL v2 (or v2+). It cannot be redistributed under LGPL v2.1 without permission from every copyright holder of the original SWAP work. For SWAP that means Wageningen University & Research (WUR), not the current modernization maintainer.

The mismatch is almost certainly accidental inheritance — TTUTIL's LGPL text was copied into the root `LICENSE` slot at some point, or the LGPL v2.1 text was dropped in as a placeholder. There is no evidence of an explicit relicensing agreement with WUR.

## Decision

The repository's license is corrected to GPL v2 to match upstream SWAP 4.2.0's license. The TTUTIL subset, when included verbatim or in derivative form, retains its LGPL v2.1 notice in its own file (preserved at build time via `subprojects/ttutil/` and at rest in `legacy/swap-4.2.0`).

Remediation committed in the same task as this ADR:

1. Root `LICENSE` replaced with the full GPL v2 text from https://www.gnu.org/licenses/old-licenses/gpl-2.0.txt.
2. `NOTICE.md` added clarifying that SWAP portions are GPL v2 and bundled TTUTIL portions are LGPL v2.1.

## Consequences

- **Positive**: the modernization repo complies with GPL v2 §2(b) going forward. No legal ambiguity for downstream users.
- **Positive**: TTUTIL's LGPL v2.1 is preserved at the per-file level, where it was always intended to live.
- **Negative**: anyone who previously relied on the root-level LGPL v2.1 notice to justify linking SWAP into proprietary code was mistaken. This ADR + fix makes the actual licensing explicit and forecloses that interpretation.
- **Negative**: the correction is a historical fix; downstream consumers who pulled the repo during the LGPL-v2.1-at-root period received it under that notice. We cannot retroactively change their copy's terms. Going forward, every new clone sees GPL v2.

## When to revisit

Closed by the remediation commit. If upstream SWAP ever relicenses (later GPL version, dual-license), that triggers ADR 0006.
```

- [ ] **Step 5: Commit all three together**

```bash
cd /home/zawadzkim/Code/swap
git add LICENSE NOTICE.md docs/adr/0005-licensing-audit.md
git status --short
git commit -m "$(cat <<'EOF'
chore(license): replace root LICENSE with GPL v2; add NOTICE.md; ADR 0005

Phase 0 discovered the repo-root LICENSE was LGPL v2.1 while
upstream SWAP 4.2.0 (legacy/swap-4.2.0:source_swp_4.2.0/LICENSE)
declares GPL v2. A derivative work of GPL v2 must be licensed
GPL v2; the mismatch was almost certainly accidental inheritance
from TTUTIL.

ADR 0005 records the finding, decision, and rationale.

Remediation:
- LICENSE now contains the full GPL v2 text from
  https://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
- NOTICE.md documents the two-license situation (SWAP = GPL v2;
  bundled TTUTIL = LGPL v2.1 with per-file headers inside
  subprojects/ttutil/ and preserved on legacy/swap-4.2.0).

This is a corrective change, not a policy change. Downstream users
who pulled the repo during the LGPL-in-root-LICENSE period hold
their copy under that notice; going forward every new clone sees
GPL v2.
EOF
)"
```

- [ ] **Step 6: Verify check-fast**

```bash
pixi run -e test check-fast 2>&1 | tail -5
```

Expected: 4/4 pass.

---

## Task 14: Per-subdir READMEs and refresh `src/README.md`

**Files:**
- Modify: `src/README.md`
- Create: `src/atmosphere/README.md`, `src/boundary/README.md`, `src/core/README.md`, `src/crop/README.md`, `src/drainage/README.md`, `src/error/README.md`, `src/heat/README.md`, `src/io/README.md`, `src/macropore/README.md`, `src/soil/README.md`, `src/solute/README.md`, `src/utils/README.md`.

Single commit, 13 files.

- [ ] **Step 1: Inventory per-subdir content**

For each subdir:

```bash
cd /home/zawadzkim/Code/swap
for d in atmosphere boundary core crop drainage error heat io macropore soil solute utils; do
  echo "=============== src/$d ==============="
  ls src/$d/
  echo "--- exports ---"
  grep -nE "^\s*public\s*::|^\s*public$" src/$d/*.f90 2>/dev/null | head -20
  echo "--- external uses ---"
  grep -hnE "^\s*use\s+" src/$d/*.f90 2>/dev/null | awk '{print $2}' | sort -u | head -20
done
```

Record per subdir: file list, public exports, external `use` modules (filter to other `src/` domains or `swap_*`).

- [ ] **Step 2: Write per-subdir READMEs**

Template — each README has these three sections exactly:

```markdown
# src/<domain>

## Responsibility

<one-paragraph description>

## Public interface

- `<domain>_state_t` (type) — <brief>
- `<domain>_state_init(state, config)` — <brief>
- ... <other exports>

## Dependencies

Imports from:
- `src/core/` — <modules>
- `src/utils/` — <modules>
- ...

Does NOT depend on: <siblings if non-obvious>
```

Per-subdir content (guidance — tune to the actual inventory):

- **atmosphere** — precipitation, interception, ET, snow, meteo day/timestep. Exports `atmosphere_state_t` + ET/interception procedures. Depends on `src/core/`, `src/utils/`.
- **boundary** — top/bottom boundary conditions. Exports `boundary_state_t` + `boundbottom` / `boundtop`.
- **core** — entry points (`swap_main.f90`, `swap.f90`, `initialize.f90`, `timecontrol.f90`), aggregated state (`swap_state_mod.f90`), legacy `variables` module + rescue bridge (`swap_state_sync.f90`), logging, constants. Most densely dependent.
- **crop** — crop growth, WOFOST soil interaction, irrigation, tillage, oxygenstress, rootextraction, management. Exports `cropgrowth_state_t`. Most complex domain.
- **drainage** — drainage calculations, surface water state.
- **error** — Empty at rescue baseline. Placeholder for Phase 4 error-handling module.
- **heat** — soil temperature, frozen soil conductivity.
- **io** — TOML readers, legacy fixed-format reader, meteo reader, CSV output.
- **macropore** — macropore flow, rates, output.
- **soil** — soil hydraulics, soil grid, water balance, tabulated soil properties.
- **solute** — solute transport.
- **utils** — array utilities, numerical solvers, hydraulics helpers, surface-water helpers, shared exchange, shared simulation utilities.

Each README ~30–60 lines.

- [ ] **Step 3: Refresh `src/README.md`**

Overwrite with:

```markdown
# src/

Source code for the SWAP modernization.

Architecture: [../docs/architecture.md](../docs/architecture.md).
State pattern: [../docs/state-management.md](../docs/state-management.md).
Code style: [../docs/code-style.md](../docs/code-style.md).
Contributing: [../docs/contributing.md](../docs/contributing.md).

## Subdirectories

| Directory | What it owns |
|---|---|
| [atmosphere/](atmosphere/README.md) | Precipitation, interception, ET, snow, meteo. |
| [boundary/](boundary/README.md) | Top / bottom boundary conditions. |
| [core/](core/README.md) | Entry points, aggregated state, time control, legacy bridge. |
| [crop/](crop/README.md) | Crop growth (fixed / grass / WOFOST), irrigation, tillage, rootextraction. |
| [drainage/](drainage/README.md) | Drainage flux, surface water state. |
| [error/](error/README.md) | (Empty — Phase 4 placeholder.) |
| [heat/](heat/README.md) | Soil temperature, frozen soil conductivity. |
| [io/](io/README.md) | TOML and legacy readers, CSV output writers. |
| [macropore/](macropore/README.md) | Macropore flow, rates, output. |
| [soil/](soil/README.md) | Soil hydraulics, grid, water balance. |
| [solute/](solute/README.md) | Solute transport. |
| [utils/](utils/README.md) | Arrays, solvers, shared helpers. |

## Dependency direction

Physics subdirs depend on `core/` (state aggregator, time control) and `utils/` (low-level helpers). They do not depend on each other directly; inter-domain communication is through `swap_state_t`.

## Conventions

See [../docs/code-style.md](../docs/code-style.md). The linter config is `.fprettify.rc` at the repo root; run `pixi run lint-check` before committing.
```

- [ ] **Step 4: Commit**

```bash
cd /home/zawadzkim/Code/swap
git add src/README.md src/*/README.md
git status --short
git commit -m "docs: per-subdir READMEs for every src/<domain>

13 files: src/README.md refreshed + one README per subdir
(atmosphere, boundary, core, crop, drainage, error, heat, io,
macropore, soil, solute, utils). Each has three sections:
Responsibility, Public interface, Dependencies.
src/error/README.md notes the directory is empty at rescue
baseline and is a Phase 4 placeholder."
```

- [ ] **Step 5: Verify check-fast**

```bash
pixi run -e test check-fast 2>&1 | tail -5
```

Expected: 4/4 pass.

---

## Task 15: Update FORD config and verify docs build

**Files:**
- Modify: `docs.md`

- [ ] **Step 1: Inspect current docs.md**

```bash
cat /home/zawadzkim/Code/swap/docs.md
```

- [ ] **Step 2: Update `exclude_dir` to also skip `docs/superpowers`**

Find the `exclude_dir:` line(s). Add `./docs/superpowers` so FORD does not publish internal planning docs as pages.

Result should look like:

```
exclude_dir: ./subprojects
             ./docs/superpowers
```

- [ ] **Step 3: Build docs**

```bash
cd /home/zawadzkim/Code/swap
pixi run -e docs docs-build 2>&1 | tee /tmp/phase2-task15-ford.log | tail -30
```

Expected: FORD parses `docs.md`, picks up the new top-level docs, produces `docs/api/index.html`. Warnings about missing source docstrings are acceptable (Phase 4 concern); hard errors are not.

- [ ] **Step 4: Spot check**

```bash
ls docs/api/page/ 2>/dev/null
```

Expected: page files for each top-level doc (`architecture.html`, `state-management.html`, `build-and-test.html`, etc.) and for ADRs. If `superpowers/` appears, the exclude failed — fix `docs.md` syntax.

- [ ] **Step 5: Commit**

```bash
cd /home/zawadzkim/Code/swap
git add docs.md
git status --short
```

Only `docs.md` should be staged. `docs/api/` FORD output is not committed here.

```bash
git commit -m "docs(ford): exclude docs/superpowers from FORD output

exclude_dir skips docs/superpowers/ — internal rescue specs and
plans are not user-facing FORD pages. All other FORD settings
unchanged. pixi run -e docs docs-build succeeds and renders the
new top-level docs plus five ADRs."
```

- [ ] **Step 6: Verify check-fast**

```bash
pixi run -e test check-fast 2>&1 | tail -5
```

Expected: 4/4 pass.

---

## Task 16: Phase 2 closeout — check-full, FF main, tag

**Files:**
- Create: `tests/regression/baselines/phase-2-docs.log`

- [ ] **Step 1: Clean build + full verification**

```bash
cd /home/zawadzkim/Code/swap
rm -rf builddir
time pixi run -e test check-full 2>&1 | tee /tmp/phase2-check-full.log | tail -30
```

Expected: 6/6 pass, ~6 minutes. If anything fails, STOP.

- [ ] **Step 2: Preserve log**

```bash
cp /tmp/phase2-check-full.log tests/regression/baselines/phase-2-docs.log
git add tests/regression/baselines/phase-2-docs.log
git commit -m "docs(baselines): record phase 2 check-full output"
```

- [ ] **Step 3: Fast-forward main**

```bash
git branch -f main development
git log --oneline -1 main
git log --oneline -1 development
git log --oneline -1 origin/main   # still 7587ca3
```

- [ ] **Step 4: Tag**

```bash
git tag -a rescue/phase-2-docs main -m "$(cat <<'EOF'
Phase 2 complete: canonical documentation + linter config.

Drift-era docs/ pages purged. New canonical docs:
  - docs/index.md, architecture.md, state-management.md,
    configuration-schema.md, build-and-test.md,
    dependency-management.md, code-style.md, contributing.md
  - docs/adr/0003-aggregator-state-over-compartments-for-now.md
  - docs/adr/0004-pfunit-for-unit-tests.md
  - docs/adr/0005-licensing-audit.md

.fprettify.rc encodes the style that code-style.md documents.
pixi lint + lint-check tasks now target src/ and tests/unit/.

Licensing remediation: root LICENSE replaced with GPL v2 per
ADR 0005; NOTICE.md documents the dual-license situation.

Per-subdir READMEs in every src/<domain>/.

FORD config excludes docs/superpowers/. pixi run -e docs docs-build
produces a clean site.

Deferred: docs/api/ tracking policy; per-procedure docstrings
inside Fortran source (Phase 4); pFUnit full vendoring; markdown
linting and link checking for docs.

check-full green at documented tolerances. main = development;
origin/main untouched at 7587ca3.
EOF
)"
```

- [ ] **Step 5: Final summary**

```bash
echo "=== branches ==="
git branch
git branch --list 'archive/*'
git branch --list 'legacy/*'
echo "=== tags ==="
git tag -l 'rescue/*'
echo "=== HEAD / main / origin ==="
git log --oneline -1 main
git log --oneline -1 development
git log --oneline -1 origin/main
echo "=== Phase 2 commits ==="
git log --oneline rescue/phase-1-infra..rescue/phase-2-docs
```

Expected: `main == development`; `origin/main` = `7587ca3`; three `rescue/*` tags; ~16 Phase 2 commits.

---

## Plan self-review

**Spec coverage.** Phase 2 spec rows mapped to tasks:

| Spec row | Plan task |
|---|---|
| Purge drift-era pages | Task 1 |
| Write index.md | Task 2 |
| architecture.md | Task 3 |
| state-management.md | Task 4 |
| configuration-schema.md | Task 5 |
| build-and-test.md | Task 6 |
| dependency-management.md | Task 7 |
| code-style.md | Task 8 |
| contributing.md | Task 10 |
| ADRs (0001–0004 + 0005 licensing) | 0001/0002 in Phase 1; 0003 Task 11; 0004 Task 12; 0005 + remediation Task 13 |
| Per-subdir READMEs | Task 14 |
| Update docs.md (FORD) | Task 15 |
| Phase tag | Task 16 |

Added beyond spec: Task 9 (`.fprettify.rc` and `lint-check` task) — user-requested addition so what `code-style.md` documents is what `pixi run lint` enforces.

**Placeholder scan.** No "TBD"/"TODO" in step instructions. Content-generation tasks include section outlines and reconnaissance commands grounding the subagent in actual repo state. One deliberate stop-condition: Task 13 Step 2 — if network is unavailable for fetching the GPL v2 text, the subagent reports BLOCKED rather than inventing license text.

**Type / name consistency.** File names match across tasks (`docs/architecture.md`, `docs/adr/000{3,4,5}-*.md`, `.fprettify.rc`, `LICENSE`, `NOTICE.md`, `src/*/README.md`, `docs.md`). Task-cross-references are all internal to this plan (no broken pointers).

**Known risks:**

1. Task 9 (linter): the `.fprettify.rc` defaults chosen may produce a larger diff than expected. Step 4 explicitly stops the task if the diff is large or physics-affecting. Tighter conservative config produced if so.
2. Task 13 (license): needs network access for the GPL v2 text. BLOCKED if unavailable — controller pastes text in a follow-up.
3. Task 15 (FORD): FORD may complain about something in the new layout. Fix is `docs.md` tweaks.
4. Per-subdir READMEs (Task 14): risk of inconsistency across 13 files; the subagent drafts all 13 in one pass using the template.
